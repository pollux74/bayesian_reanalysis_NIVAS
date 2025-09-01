# =============================================
# Projet : Bayesian reanalysis of NIVAS
# Script : 02_models.R
# Author : A. Naudet--Lasserre
# Date last modification : 2025-09-01
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
explic_reint_multi1 <- c("group + age_inf_60 + copd + surg_site + epidural + sepsis + cardio_isch + chronic_cardio_faill + bmi_sup_30 + (1 | center) ")

model_list_reint_d7_multivariate1 <- fit_brm_opti(
  data = df,
  outcome = "reintubation_d7",
  explic = explic_reint_multi1
)

explic_reint_multi2 <- c("group + age_inf_60 + copd + surg_site + epidural + sepsis + cardio_isch + chronic_cardio_faill + bmi_sup_30")
model_list_reint_d7_multivariate2 <- fit_brm_opti(
  data = df,
  outcome = "reintubation_d7",
  explic = explic_reint_multi2
)

# Intercept models for reintubation
model_list_reint_inter1 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter1)
model_list_reint_inter2 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter2)
model_list_reint_inter3 <- fit_brm_opti(data = df, outcome = "reintubation_d7", explic = "group", intercept = P_reintubation_inter3)

## Day-30 reintubation
model_list_reint30 <- fit_brm_opti(data = df, outcome = "reintubation_d30", explic = "group")

# DAyt-30 mortality
model_list_mortality_d30 <- fit_brm_opti(df, outcome = "mortality_d30", explic = "group")
explic_reint_mortality30 <- c("group + surg_site + epidural + copd + sofa + age + igs")

model_list_mortalityd30_multivariate1 <- fit_brm_opti(
  data = df,
  outcome = "mortality_d30",
  explic = explic_reint_mortality30
)
# Day-90 mortality

model_list_mortality_d90 <- fit_brm_opti(df, "mortality_d90", explic = "group")
# Health care associated infection
model_list_hcai_d7 <- fit_brm_opti(df, "hcai_d7", explic = "group")
model_list_hcai_d30 <- fit_brm_opti(df, "hcai_d30", explic = "group")
# pneumonia
model_list_pneumonia_d7 <- fit_brm_opti(df, "pneumonia_d7", explic = "group")
model_list_pneumonia_d30 <- fit_brm_opti(df, "pneumonia_d30", explic = "group")


# LOS

model_list_hosLOS <- fit_brm_opti(df, outcome = "hosLOS", "group")
model_list_icuLOS <- fit_brm_opti(df, outcome = "icuLOS", "group")

# Results -----------------------------------------------------------------

# Day-7 reintubation

analyse_outcome("reintubation_d7", df, "group", model_list_reint)

analyse_outcome("reintubation_d7", df, explic_reint_multi1, model_list_reint_d7_multivariate1)
analyse_outcome("reintubation_d7", df, explic_reint_multi2, m_list = model_list_reint_d7_multivariate2)

# Comparison of models adjustement
## WAIC comparison
waic_comparison <- model_weights(
  model_list_reint_d7_multivariate1[[1]], # random effect (1 | center)
  model_list_reint_d7_multivariate2[[1]], # no adjustement for study center
  weights = "waic"
)
names(waic_comparison) <- c("Random_effect", "Without_center_adjustement")
print(waic_comparison)

analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter1)
analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter2)
analyse_outcome("reintubation_d7", df, "group", model_list_reint_inter3)

# Day-30 reintubation
analyse_outcome(outcome = "reintubation_d30", df, "group", model_list_reint30)

# Day-30 Mortality -------------------------------------------------------------------
analyse_outcome(outcome = "mortality_d30", df, "group", model_list_mortality_d30)
analyse_outcome(outcome = "mortality_d30", df, "group", model_list_mortalityd30_multivariate1)
analyse_outcome(
  outcome = "mortality_d30",
  df = df,
  predicteur = explic_reint_mortality30,
  m_list = model_list_mortalityd30_multivariate1,
  ARR = T,
  verbose = T
)

analyse_outcome("mortality_d90", df, "group", model_list_mortality_d90)

# hcai --------------------------------------------------------------------
analyse_outcome("hcai_d7", df, "group", model_list_hcai_d7)
analyse_outcome("hcai_d30", df, "group", model_list_hcai_d30)
# pneumonia ---------------------------------------------------------------
analyse_outcome("pneumonia_d7", df, "group", model_list_pneumonia_d7)
analyse_outcome("pneumonia_d30", df, "group", model_list_pneumonia_d30)

# LOS
analyse_outcome("hosLOS", df, "group", model_list_hosLOS)
analyse_outcome("icuLOS", df, "group", model_list_icuLOS)


# Convergence -------------------------------------------------------------

check_convergence(model_list_reint)
check_convergence(model_list_reint_d7_multivariate1)
check_convergence(model_list_reint_inter1)
check_convergence(model_list_reint_inter2)
check_convergence(model_list_reint_inter3)
check_convergence(model_list_mortality_d30)
check_convergence(model_list_mortality_d90)
check_convergence(model_list_hcai_d7)
check_convergence(model_list_pneumonia_d7)
check_convergence(model_list_pneumonia_d30)
check_convergence(model_list_hcai_d30)
check_convergence(model_list_reint30)
check_convergence(model_list_icuLOS)
check_convergence(model_list_hosLOS)


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

# Subgroup analysis -------------------------------------------------------

# Age subgroups
df_young <- df[df$age_inf_60 == "Oui" & !is.na(df$age_inf_60), ]
df_old <- df[df$age_inf_60 == "Non" & !is.na(df$age_inf_60), ]

model_young <- fit_brm_opti(df_young, "reintubation_d7", "group")
model_old <- fit_brm_opti(df_old, "reintubation_d7", "group")

# Surgery site subgroups
df_below <- df[df$surg_site == "Sous" & !is.na(df$surg_site), ]
df_above <- df[df$surg_site == "Sus" & !is.na(df$surg_site), ]

model_below <- fit_brm_opti(df_below, "reintubation_d7", "group")
model_above <- fit_brm_opti(df_above, "reintubation_d7", "group")

# Epidural subgroups
df_epi_yes <- df[df$epidural == "Oui" & !is.na(df$epidural), ]
df_epi_no <- df[df$epidural == "Non" & !is.na(df$epidural), ]

model_epi_yes <- fit_brm_opti(df_epi_yes, "reintubation_d7", "group")
model_epi_no <- fit_brm_opti(df_epi_no, "reintubation_d7", "group")




# Extract results from subgroup models
extract_results <- function() {
  results <- list()

  # Helper function
  get_model_results <- function(model_list, name, subgroup, level, n) {
    if (!is.null(model_list[[1]])) {
      model_results <- compute_posterior_metrics(
        model_list = model_list[1],
        prior_names = name,
        group_name = "groupNIV"
      )

      return(data.frame(
        subgroup = subgroup,
        level = level,
        median_or = model_results$summary$median_or,
        ci_lower = model_results$summary$ci_lower,
        ci_upper = model_results$summary$ci_upper,
        n_patients = n
      ))
    }
    return(NULL)
  }

  # Overall analysis
  if (exists("model_list_reint") && !is.null(model_list_reint[[1]])) {
    overall <- get_model_results(model_list_reint, "Overall", "Overall", "All patients", nrow(df))
    if (!is.null(overall)) results[["overall"]] <- overall
  }

  # Age subgroups
  young <- get_model_results(model_young, "Young", "Age", "< 60 years", 112)
  if (!is.null(young)) results[["young"]] <- young

  old <- get_model_results(model_old, "Old", "Age", "≥ 60 years", 181)
  if (!is.null(old)) results[["old"]] <- old

  # Surgery site
  below <- get_model_results(model_below, "Below", "Surgery site", "Below diaphragm", 109)
  if (!is.null(below)) results[["below"]] <- below

  above <- get_model_results(model_above, "Above", "Surgery site", "Above diaphragm", 184)
  if (!is.null(above)) results[["above"]] <- above

  # Epidural
  epi_yes <- get_model_results(model_epi_yes, "EpiYes", "Epidural", "Yes", 44)
  if (!is.null(epi_yes)) results[["epi_yes"]] <- epi_yes

  epi_no <- get_model_results(model_epi_no, "EpiNo", "Epidural", "No", 249)
  if (!is.null(epi_no)) results[["epi_no"]] <- epi_no

  # Combine all results
  final_data <- do.call(rbind, results)
  rownames(final_data) <- NULL

  return(final_data)
}

# Create classic forest plot
create_forest_plot <- function(data) {
  # Prepare data
  plot_data <- data %>%
    mutate(
      label = paste0(level, " (n=", n_patients, ")"),
      or_text = sprintf("%.2f (%.2f-%.2f)", median_or, ci_lower, ci_upper),

      # Set display order
      order = case_when(
        subgroup == "Overall" ~ 1,
        subgroup == "Age" & level == "< 60 years" ~ 2,
        subgroup == "Age" & level == "≥ 60 years" ~ 3,
        subgroup == "Surgery site" & level == "Below diaphragm" ~ 4,
        subgroup == "Surgery site" & level == "Above diaphragm" ~ 5,
        subgroup == "Epidural" & level == "Yes" ~ 6,
        subgroup == "Epidural" & level == "No" ~ 7
      ),

      # Group labels for facets
      group_name = case_when(
        subgroup == "Overall" ~ "Overall",
        subgroup == "Age" ~ "Age groups",
        subgroup == "Surgery site" ~ "Surgical site",
        subgroup == "Epidural" ~ "Epidural analgesia"
      )
    ) %>%
    arrange(order)

  # Create plot
  p <- ggplot(plot_data, aes(y = reorder(label, desc(order)))) +

    # Reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +

    # Point estimates and CIs
    geom_pointrange(
      aes(x = median_or, xmin = ci_lower, xmax = ci_upper),
      size = 0.8,
      linewidth = 1.2,
      shape = 21,
      fill = "white",
      stroke = 1,
      color = "#1F78B4" # Blue color
    ) +

    # OR values as text
    geom_text(
      aes(x = max(plot_data$ci_upper) * 1.1, label = or_text),
      size = 3.2,
      hjust = 0,
      family = "Arial"
    ) +

    # Log scale
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
      labels = c("0.1", "0.25", "0.5", "1.0", "2.0", "4.0", "8.0"),
      limits = c(0.1, max(plot_data$ci_upper) * 1.3)
    ) +

    # Facets by subgroup
    facet_grid(
      group_name ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +

    # Theme similar to previous version
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      # Remove y-axis title
      axis.title.y = element_blank(),

      # Customize x-axis
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 12),

      # Customize y-axis labels
      axis.text.y = element_text(size = 12, hjust = 0),

      # Customize facet labels
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 12),
      strip.placement = "outside",

      # Remove minor grid, keep major x grid
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),

      # Add panel borders
      panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5),

      # Adjust margins
      plot.margin = margin(10, 100, 10, 10)
    ) +

    # Labels
    labs(x = "OR (95% CrI)")

  return(p)
}

# Run the analysis
forest_data <- extract_results()
forest_plot <- create_forest_plot(forest_data)

# Display plot
print(forest_plot)

ggsave(
  plot = forest_plot,
  filename = fig_path("forest_plot_classic.tiff"),
  width = 7, # Standard journal width in inches
  height = 5, # Height in inches (ratio ~1.4)
  dpi = 600, # 600 dpi for print quality
  compression = "lzw", # Lossless compression
  device = "tiff" # Explicit format specification
)

# Show results table
print(forest_data)
