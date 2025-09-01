# =============================================
# Project : Bayesian reanalysis of NIVAS
# Script : data_preparation.R
# Author : A. Naudet--Lasserre
# Date last modification : 2025-09-01
# =============================================

# Load data base(not available for patient data privacy concern) ----------------------------------------------------------
df <- read.csv(here("Data", "nivas_data_raw.csv"), stringsAsFactors = FALSE, sep = ";")
str(df)

# Data_management ---------------------------------------------------------

## Primary outcome: reintubation at day-7 ----
df$reintubation_d7 <- factor(df$intub_J7, levels = c("0", "1"), labels = c("No", "Yes"))
## fay-30 reintubation ---
df$reintubation_d30 <- factor(df$intub_J30, levels = c("0", "1"), labels = c("No", "Yes"))
## Treatment group ----
df$group <- factor(df$GROUPE)
df$group <- recode(
  df$group,
  "BRAS A" = "0", # Oxygen group
  "BRAS B" = "1" # NIV group
)
df$group <- factor(df$group, levels = c("0", "1"), labels = c("O2_standard", "NIV"))
table(df$reintubation_d7, df$group)
## Age ----
df$age <- as.numeric(df$AGE)
summary(df$age)
## Age >60 ans
df$age_inf_60 <- ifelse(df$age < 60, 1, 0)
df$age_inf_60 <- factor(df$age_inf_60, levels = c("0", "1"), labels = c("Non", "Oui"))
## death to day 30 ----
df$mortality_d30 <- factor(df$suividcJ30, levels = c("0", "1"), labels = c("non", "oui"))
table(df$mortality_d30, df$suividcJ30)
## death to day 90 ----
df$mortality_d90 <- factor(df$dcbilanj90, levels = c("0", "1"), labels = c("non", "oui"))
table(df$mortality_d90, df$dcbilanj90)
## day-7 HCAI
df$hcai_d7 <- factor(df$INFNOSO_ON3, levels = c("0", "1"), labels = c("non", "oui"))
table(df$group, df$hcai_d7, deparse.level = 2)
## day-30 HCAI
df$hcai_d30 <- factor(df$INFNOSOglobal_ON5, levels = c("0", "1"), labels = c("non", "oui"))
table(df$group, df$hcai_d30, deparse.level = 2)
## Day-7 Pneumonia
df$pneumonia_d7 <- factor(df$PNEUMONIE_ON3bis, levels = c("0", "1"), labels = c("Non", "Oui"))
table(df$pneumonia_d7, df$PNEUMONIE_ON3bis, deparse.level = 2)
## Day-30 pneumonia0 ----
df$pneumonia_d30 <- factor(df$PNEUMONIE_ON5global, levels = c("0", "1"), labels = c("Non", "Oui"))
table(df$pneumonia_d30, df$PNEUMONIE_ON5global, deparse.level = 2)

## Hospital LOS
df$hosLOS <- as.numeric(df$Dureetotalehospit_jours)
summary(df$hosLOS)
## ICU LOS
df$icuLOS <- as.numeric(df$nbjourreaj30)
summary(df$icuLOS)


## BMI ----
df$bmi <- as.numeric(gsub(",", ".", df$BMI))
df$bmi_sup_30 <- ifelse(df$bmi > 30, 1, 0)
summary(df$bmi)
## COPD ----
df$copd <- factor(df$BRONCHOPNEUMO_ON, levels = c("0", "1"), labels = c("Non", "Oui"))
## Ischemic heart disease ----
df$cardio_isch <- factor(df$CORO_ON, levels = c("0", "1"), labels = c("Non", "Oui"))
# Chronique heart faillure
df$chronic_cardio_faill <- factor(df$INSUFCARDIAQUE_ON, levels = c("0", "1"), labels = c("Non", "Oui"))
## Sepsis ----
df$sepsis <- factor(df$SEPSIS_ON, levels = c("0", "1"), labels = c("Non", "Oui"))
## Study site
df$center <- factor(df$CENTRE)
levels(df$center)
## Surgical site
df$surg_site <- factor(df$OPERATION_T, levels = c("1", "2"), labels = c("Sus", "Sous"))
levels(df$surg_site)
table(df$surg_site, df$OPERATION_T)
## Epidural anesthessia
df$epidural <- factor(df$ANALGESIEAPD, levels = c("0", "1"), labels = c("Non", "Oui"))
## SOFA score
df$sofa <- as.numeric(df$SOFAOTAL)
table(df$sofa)
## IGS score
df$igs <- as.numeric(df$SCOREIGS)
table(df$igs)

# # Synthetic dataset for bivariate analysis reproduction ----------------------
#
# This code creates a synthetic dataset that reproduces the key statistics
# from our study. The original data cannot be shared due to privacy requirements.
#
# Dataset: 293 patients, 2 treatment groups
# - Group 1 (NIV): n=148 | Group 0 (O2): n=145
# - Variables: treatment group, reintubation, age, death, infections, pneumonia

# Uncomment to generate the dataset:

# df <- data.frame(
#   group = c(rep("1", 148), rep("0", 145)), # 1=NIV, 0=O2 therapy
#   reintubation = c(rep(1, 49), rep(0, 99), rep(1, 66), rep(0, 79)),
#   age = c(rnorm(293, 62.1, 13.8)), # Mean age: 62.1 years
#   death = c(rep(1, 15), rep(0, 133), rep(1, 22), rep(0, 122)),
#   hcai = c(rep(1, 27), rep(0, 121), rep(1, 44), rep(0, 101)),
#   pneumonia = c(rep(1, 15), rep(0, 133), rep(1, 32), rep(0, 113))
# )

# df$group <- factor(df$group, levels = c("0", "1"))
