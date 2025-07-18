# =============================================
# Projet : Bayesian reanalysis of NIVAS
# Script : 02_models.R
# Auteur : A. Naudet
# Date création : 2025-05-13
# =============================================

# Loading -----------------------------------------------------------------
library(here)
source(here("R", "setup.R"))
source(here("R", "data_preparation.R"))
source(here("R", "priors_specification.R"))
# set_cmdstan_path("~/R/cmdstan/cmdstan-2.36.0/")

# Fit models -------------------------------------------------------

## Primary outcome - Day7 reintubation ----

model_list_reint <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group")

## Primary outcome - Day7 reintubation - multivariate model ----
explic_reint_multi1 <- c("group + center + age_inf_60 + copd + surg_site + epidural + sepsis + cardio_isch + chronic_cardio_faill + bmi_sup_30")
model_list_reint_d7_multivariate1 <- fit_brm_opti(
  data = df,
  outcome = "reintubation_d7",
  explic = explic_reint_multi1
)

# # Intercept models for reintubation
model_list_reint_inter1 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter1)
model_list_reint_inter2 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter2)
model_list_reint_inter3 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter3)



# Death on day 30

model_list_mortality_d30 <- fit_brm_opti(df, "mortality_d30", "group")
explic_reint_mortality30 <- c("group + center + surg_site + epidural + copd + sofa + age + igs")

model_list_mortalityd30_multivariate1 <- fit_brm_opti(
  data = df,
  outcome = "mortality_d30",
  explic = explic_reint_mortality30
)
# Health care associated infection
model_list_hcai_d7 <- fit_brm_opti(df, "hcai_d7", "group")
# pneumonia to day-7
model_list_pneumonia_d7 <- fit_brm_opti(df, "pneumonia_d7", "group")


# Results -----------------------------------------------------------------

# Day-7 reintubation

analyse_outcome("reintubation_d7", df, "group", model_list_reint)



analyse_outcome("reintubation_d7", df, explic_reint_multi1, model_list_reint_d7_multivariate1)
analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter1)
analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter2)
analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter3)

# Death -------------------------------------------------------------------
analyse_outcome("mortality_d30", df, "group", model_list_mortality_d30)
analyse_outcome("mortality_d30", df, "group", model_list_mortalityd30_multivariate1)
analyse_outcome(
  outcome = "mortality_d30",
  df = df,
  predicteur = explic_reint_mortality30,
  m_list = model_list_mortalityd30_multivariate1,
  ARR = T,
  verbose = T
)

# hcai --------------------------------------------------------------------
analyse_outcome("hcai_d7", df, "group", model_list_hcai_d7)
# pneumonia ---------------------------------------------------------------
analyse_outcome("pneumonia_d7", df, "group", model_list_pneumonia_d7)


# Convergence -------------------------------------------------------------

check_convergence(model_list_reint)
check_convergence(model_list_reint_d7_multivariate1)
check_convergence(model_list_reint_inter1)
check_convergence(model_list_reint_inter2)
check_convergence(model_list_mortality_d30)
check_convergence(model_list_hcai_d7)
check_convergence(model_list_pneumonia_d7)


# Fig. 2 ------------------------------------------------------------------
#
#

# Extract model data
draws_reint_d7_weakprior <- get_draws(model_list_reint[[1]], "groupNIV")
draws_reint_d7_weakprior <- as.data.frame(draws_reint_d7_weakprior)

# Convert ARR to percentage for easier interpretation
draws_reint_d7_weakprior$arr_pct <- draws_reint_d7_weakprior$arr * 100

# Calculate key probabilities based on ARR thresholds
p_no_benefit <- mean(draws_reint_d7_weakprior$arr_pct < 0) # P(ARR < 0%)
p_weak_benefit <- mean(draws_reint_d7_weakprior$arr_pct >= 2) # P(ARR ≥ 2%)
p_moderate_benefit <- mean(draws_reint_d7_weakprior$arr_pct >= 5) # P(ARR ≥ 5%)
p_strong_benefit <- mean(draws_reint_d7_weakprior$arr_pct >= 10) # P(ARR ≥ 10%)

# Calculate density on ARR percentage
dens <- density(draws_reint_d7_weakprior$arr_pct)
dens_df <- data.frame(x = dens$x, y = dens$y)



# Define clinically relevant ARR thresholds (in %)
threshold_neutral <- 0 # No benefit threshold
threshold_weak <- 2 # Weak benefit threshold
threshold_moderate <- 5 # Moderate benefit threshold
threshold_strong <- 10 # Strong benefit threshold

# Create mutually exclusive zones based on ARR magnitude
dens_df <- dens_df %>%
  mutate(
    zone = case_when(
      x < threshold_neutral ~ "no_meaningful_benefit", # ARR < 0%
      x >= threshold_neutral & x < threshold_weak ~ "weak_benefit", # 0% ≤ ARR < 2%
      x >= threshold_weak & x < threshold_moderate ~ "moderate_benefit", # 2% ≤ ARR < 5%
      x >= threshold_moderate & x < threshold_strong ~ "good_benefit", # 5% ≤ ARR < 10%
      x >= threshold_strong ~ "strong_benefit" # ARR ≥ 10%
    )
  )



# Corrected color scheme as requested
colors_publication <- c(
  "no_meaningful_benefit" = "#179299", # No meaningful benefit < 0%
  "weak_benefit" = "#40a02b", # Weak 0-2%
  "moderate_benefit" = "#df8e1d", # Moderate 2-5%
  "good_benefit" = "#fe640b", # Good 5-10%
  "strong_benefit" = "#d20f39" # Strong >10%
)


fig2_arr <- ggplot(dens_df, aes(x = x)) +

  # Filled areas
  geom_ribbon(
    aes(ymin = 0, ymax = y, fill = zone),
    alpha = 0.8
  ) +

  # Density line
  geom_line(
    aes(y = y),
    color = "black",
    linewidth = 1.0
  ) +

  # Reference lines for ARR thresholds
  geom_vline(
    xintercept = threshold_neutral,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = threshold_weak,
    linetype = "dotted",
    color = "black",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = threshold_moderate,
    linetype = "dotted",
    color = "black",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept = threshold_strong,
    linetype = "dotted",
    color = "black",
    linewidth = 0.6
  ) +

  # Manual colors (no legend)
  scale_fill_manual(
    values = colors_publication,
    guide = "none" # Remove legend
  ) +


  # Publication theme
  theme_minimal() +
  theme(
    # Text - ALL BLACK
    text = element_text(color = "black", family = "Arial"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Axes
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5),

    # Grid
    panel.grid.major = element_line(color = "#e6e9ef", linewidth = 0.3),
    panel.grid.minor = element_blank(),

    # Margins
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10, unit = "pt")
  ) +

  # Labels
  labs(
    x = "Absolute Risk Reduction (%)",
    y = "Posterior Density"
  ) +

  # Single axis with ARR percentage scale
  scale_x_continuous(
    breaks = sort(c(seq(-10, 30, by = 5), 2)),
    labels = paste0(sort(c(seq(-10, 30, by = 5), 2)), "%")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))



save_arr_plot <- function(plot_obj, filename, width = 7, height = 5, dpi = 600) {
  ggsave(
    plot = plot_obj,
    filename = filename,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = "png",
    bg = "white"
  )
}


# Display plot
print(fig2_arr)

# Save plot
save_arr_plot(fig2_arr, fig_path("fig2_arr_final.png"))


#
#
#
