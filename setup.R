# =============================================
# Project : Bayesian reanalysis of NIVAS
# Script : setup.R
# Author : A. Naudet--Lasserre
# Creation date : 2025-05-13
# =============================================

# Reproducibility settings
set.seed(42) # For any future random operations
options(digits = 4)

# Project validation ------------------------------------------------------

# Verify script is run from project root (parent of R/)
if (!dir.exists("R")) {
  stop(
    "Run this script from project root directory (should contain R/ folder)\n",
    "Current directory: ", getwd(),
    call. = FALSE
  )
}
cat("✓ Project setup initiated from:", basename(getwd()), "\n")

# Package loading ---------------------------------------------------------

packages <- c(
  "ggplot2", "bayesplot", "gridExtra", "patchwork", "ggdist", # Visualization
  "brms",
  # "cmdstanr", # Bayesian
  "tidyverse", "here" # Data & utils
)

# Install missing packages
missing <- packages[!packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))
cat("✓ Loaded", length(packages), "packages\n")

# Configuration -----------------------------------------------------------

# Stan setup for bayesian modeling
options(brms.backend = "cmdstanr")

# Create output directories if needed
# output_dirs <- c("outputs/figures", "outputs/tables", "outputs/models")
# for (dir_path in output_dirs) {
#   if (!dir.exists(here::here(dir_path))) {
#     dir.create(here::here(dir_path), recursive = TRUE)
#   }
# }

# Load custom functions
if (file.exists(here::here("R", "functions.R"))) {
  source(here::here("R", "functions.R"))
}

cat("✓ Setup complete\n")
