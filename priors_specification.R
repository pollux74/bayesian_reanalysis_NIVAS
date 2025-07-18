# =============================================
# Project : Bayesian reanalysis of NIVAS
# Script : prior_specification.R
# Author : A. Naudet--Lasserre
# Creation date : 2025-05-13
# =============================================

# I. Original NIVAS trial power assumptions ----------------------------------

# Sample size and event rate assumptions from original protocol
nivas_power_params <- list(
  n_niv = 150L, # Sample size NIV group (integer for clarity)
  n_o2 = 150L, # Sample size O2 group
  prop_niv = 0.40, # Expected mortality rate NIV group
  prop_o2 = 0.65, # Expected mortality rate O2 group (reference)
  alpha = 0.05, # Type I error (for transparency)
  power = 0.80 # Statistical power (for transparency)
)

# Display assumptions
cat("=== NIVAS Trial Power Assumptions ===\n")
cat(
  "NIV group: n =", nivas_power_params$n_niv,
  ", expected mortality =", nivas_power_params$prop_niv, "\n"
)
cat(
  "O2 group:  n =", nivas_power_params$n_o2,
  ", expected mortality =", nivas_power_params$prop_o2, "\n\n"
)

# Contingency table construction

# Calculate expected cell counts
nivas_power_table <- with(nivas_power_params, {
  matrix(
    data = c(
      prop_niv * n_niv,     (1 - prop_niv) * n_niv, # NIV: deaths, alive
      prop_o2 * n_o2,       (1 - prop_o2) * n_o2 # O2:  deaths, alive
    ),
    nrow = 2,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(
      Treatment = c("NIV", "O2_standard"),
      Outcome = c("Death", "Alive")
    )
  )
})

# Display contingency table
cat("=== Expected Contingency Table (Power Calculation) ===\n")
print(nivas_power_table)



# Calculate OR with confidence interval
or_power_calc <- epitools::oddsratio(
  x = nivas_power_table,
  method = "wald",
  conf.level = 0.95
)

# Extract key statistics
or_estimate <- or_power_calc$measure[2, "estimate"]
or_lower <- or_power_calc$measure[2, "lower"]
or_upper <- or_power_calc$measure[2, "upper"]
log_or_estimate <- log(or_estimate)

# Display results
cat("=== Odds Ratio from Power Assumptions ===\n")
cat("OR estimate:", round(or_estimate, 3), "\n")
cat("95% CI: [", round(or_lower, 3), ", ", round(or_upper, 3), "]\n")
cat("Log OR:", round(log_or_estimate, 3), "\n")


# II. Priors definitions -------------------------------------------------------

## A. Minimally informative and skeptical priors for all outcomes ----------------------------

# Define a minimally informative prior as per Goligher et al. (2018 & 2023)
prior_minimal <- data.frame(
  prior_mean = log(1), # Centered on the absence of treatment effect
  prior_sd = 10, # Large sd to reflect minimal information
  interpretation = "Representing no priors beliefs.
    Equivalent to frequentist analysis"
)

prior_skeptical <- data.frame(
  prior_mean = 0, # log(OR = 1) = no treatment effect
  prior_sd = 0.355, #
  interpretation = "95% probabily for the OR to fall between 0.5 and 2"
)

## B. Day-7 Reintubation ----
priors_reintubation_d7 <- list()
### 1. Minimally informative prior ----

priors_reintubation_d7[["Minimally informative"]] <- prior_minimal

### 2.  Enthusiastic priorr ----

priors_reintubation <- data.frame(
  prior_name = c("Weakly informative", "Skeptical", "Optimistic", "Pessimistic"),
  prior_mean = c(0, 0, log(0.3589744), log(1 / 0.3589744)),
  prior_sd = c(5, 0.355, 0.9884903, 0.9884903)
)
# Standard deviation to allow 15% chance of harm (OR > 1) (Zampieri et al.2021)
enthusiastique_sd <- calibrate_prior_sd(
  mean = log(or_estimate),
  target_p = 0.85,
  threshold = 0 # 5% probability that log(OR) > 0 (NIV harmful vs O2)
)

# Belief in NIV protective effect with 5% chance of harm (OR > 1)
priors_reintubation_d7[["Enthusiastic"]] <- data.frame(
  prior_mean = log(or_estimate), # Logarithmic OR for reintubation
  prior_sd = enthusiastique_sd,
  interpretation = "Belief in NIV protection with 15% chance of harm"
)
### 3. Pessimistic ----
# Build as the opposite of the Enthusiatic prior
priors_reintubation_d7[["Pessimistic"]] <- data.frame(
  prior_mean = -log(or_estimate), # Logarithmic OR for reintubation
  prior_sd = enthusiastique_sd,
  interpretation = "Assumes NIV si worse than oxygen therapy but allows 15% chance of benefit."
)

### 4. Skeptical prior ----

# Centered on null with small probability of large protective effect
priors_reintubation_d7[["Skeptical"]] <- prior_skeptical


priors_reintubation_d7 <- bind_rows(priors_reintubation_d7, .id = "prior_type")
print(priors_reintubation_d7)

# Round numerical columns to 1 decimal place and rename columns
priors_reintubation_d7_clean <- priors_reintubation_d7 %>%
  mutate(
    prior_mean = round(prior_mean, 2),
    prior_sd = round(prior_sd, 2)
  ) %>%
  rename(
    "Prior Belief" = prior_type,
    "Assumed Mean of Logarithm of OR" = prior_mean,
    "Assumed SD of Logarithm of OR" = prior_sd,
    "Clinical Interpretation" = interpretation
  )

priors_reintubation_d7_clean


## C. Day-30 Mortality ----
priors_mortality_d30 <- list()

### 1. Minimally informative prior ----
priors_mortality_d30[["Minimally informative"]] <- prior_minimal

### 2. Skeptical prior ----
priors_mortality_d30[["Skeptical"]] <- prior_skeptical

priors_mortality_d30 <- bind_rows(priors_mortality_d30, .id = "prior_type")
priors_mortality_d30_clean <- priors_mortality_d30 %>%
  mutate(
    prior_mean = round(prior_mean, 2),
    prior_sd = round(prior_sd, 2)
  ) %>%
  rename(
    "Prior Belief" = prior_type,
    "Assumed Mean of Logarithm of OR" = prior_mean,
    "Assumed SD of Logarithm of OR" = prior_sd,
    "Clinical Interpretation" = interpretation
  )

print(priors_mortality_d30_clean)

# print(priors_death)
priors_hcai_d7 <- priors_mortality_d30
priors_pneumonia_d7 <- priors_mortality_d30

# Visualisation des priors -------------------------------------------------------------------------
# Usage with your data (with legend)
p1 <- plot_priors_distributions(priors_reintubation_d7)
p1
# Version without legend if needed

# Pour des figures de travail (plus légères)


ggsave(
  plot = p1,
  filename = paste0(here("outputs", "figures"), "/prior_distributions.tiff"), # TIFF pour publication
  width = 7, # Largeur en pouces (standard journal)
  height = 5, # Hauteur en pouces (ratio ~1.4)
  dpi = 600, # 600 dpi pour impression
  device = "tiff" # Expliciter le format
)


# Intercept prior (eSupllemental method) ---------------------------------------------------------

# Helper functions ----------------------------------------------------------

logit_center <- function(x) {
  log(x / (1 - x))
}

logit_CI <- function(x) {
  round(log(x / (1 - x)), 2)
}

intercept_sd <- function(x) {
  abs(x) / 1.96
}

# Prior specifications ------------------------------------------------------

# Clinical baseline probabilities and bounds
baseline_center1 <- 0.65
baseline_lower1 <- 0.35
baseline_center2 <- 0.50
baseline_lower2 <- 0.20
baseline_center3 <- 0.50
baseline_lower3 <- 0.30

# Calculate logit bounds and standard deviations
intercept_ci1 <- logit_CI(baseline_lower1)
intercept_ci2 <- logit_CI(baseline_lower2)
intercept_ci3 <- logit_CI(baseline_lower3)

sd_intercept1 <- intercept_sd(logit_center(baseline_center1) - intercept_ci1)
sd_intercept2 <- intercept_sd(logit_center(baseline_center2) - intercept_ci2)
sd_intercept3 <- intercept_sd(logit_center(baseline_center3) - intercept_ci3)

# Prior dataframes for brms
P_reintubation_inter1 <- data.frame(
  prior_type = "Prior Reintubation d7  Intercept 1 ",
  prior_mean = logit_center(baseline_center1),
  prior_sd = sd_intercept1
)

P_reintubation_inter2 <- data.frame(
  prior_type = "Prior Reintubation d7  Intercept 2 ",
  prior_mean = logit_center(baseline_center2),
  prior_sd = sd_intercept2
)

P_reintubation_inter3 <- data.frame(
  prior_type = "Prior Reintubation d7  Intercept 3 ",
  prior_mean = logit_center(baseline_center3),
  prior_sd = sd_intercept3
)

# Display results
print(P_reintubation_inter1)
print(P_reintubation_inter2)
print(P_reintubation_inter3)
