##################################################
## Project: MSc HDS Project
## Date: 26 January 2026
## Author: Autumn O'Donnell
##################################################

## Section: 0.0 Load Libraries
##################################################
library(readxl)
library(dplyr)
library(table1)
library(survival)
library(survminer)

## Section: 0.1 Load Dataset
##################################################
SECOMBIT <- read_excel("data/SECOMBIT 4yr OS_raw_data_v2.xlsx")

## Section: 1.0 Clean Data
##################################################
# Clean the data to improve readability
clean_df <- SECOMBIT %>%
  mutate(
    
    # ID
    id = ordine, 
    
    # Overall Survival
    os_time   = OS,
    os_status = status,
    
    # Progression-Free Survival
    pfs_time   = `PFS TOTAL`,
    pfs_status = `Progr total`,
    
    # Treatment arm (A as baseline)
    Arm = factor(ARM, levels = c("A", "B", "C")),
    
    # Number of sites (1–2 as baseline)
    Sites = factor(sites, levels = c("1-2", ">=3"),
                   labels = c("1 - 2", ">= 3")),
    
    # LDH (Normal as baseline)
    LDH = factor(ULN_LDH, levels = c("normal", "elevated"),
                 labels = c("Normal", "Elevated")),
    
    # TMB (<10 as baseline)
    TMB = factor(TMB, levels = c("<10", ">=10"),
                 labels = c("< 10", ">= 10")),
    
    # JAK (Wild Type as baseline)
    JAK = factor(JAK, levels = c("wt", "mut"),
                 labels = c("Wild Type (Normal)", "Mutated"))
  ) %>%
  select(id, Arm, Sites, LDH, TMB, JAK, pfs_time, pfs_status, os_time, os_status)


## Section: 2.0  Summary table with table1
##################################################
table1(~ Sites + LDH + TMB + JAK | Arm,
       data = clean_df,
       overall = FALSE)


## Section: 3.0   Kaplan-Meier plots
##################################################

# OS KM plot
fit_os <- survfit(Surv(os_time, os_status) ~ Arm, data = clean_df)

ggsurvplot(fit_os, data = clean_df,
           xlab = "Time (in months)",
           ylab = "Survival Probability",
           title = "Overall Survival (OS)",
           #risk.table = TRUE
           )

# PFS KM plot
fit_pfs <- survfit(Surv(pfs_time, pfs_status) ~ Arm, data = clean_df)

ggsurvplot(fit_pfs, data = clean_df,
           xlab = "Time (in months)",
           ylab = "Survival Probability",
           title = "Progression-Free Survival (PFS)",
           #risk.table = TRUE
           )

## Section: 3.1   Log-rank tests by Arm
##################################################

logrank_os <- survdiff(Surv(os_time, os_status) ~ Arm, data = clean_df)
logrank_pfs <- survdiff(Surv(pfs_time, pfs_status) ~ Arm, data = clean_df)

# Print results
logrank_os
logrank_pfs


## Section: 4.0    Cox model
##################################################

# OS Cox model
cox_os <- coxph(Surv(os_time, os_status) ~ Arm + Sites + LDH + TMB + JAK,
                data = clean_df, ties = "efron")
summary(cox_os)

# Test proportional hazards for OS
ph_test_os <- cox.zph(cox_os)
ph_test_os

# PFS Cox model
cox_pfs <- coxph(Surv(pfs_time, pfs_status) ~ Arm + Sites + LDH + TMB + JAK,
                 data = clean_df, ties = "efron")
summary(cox_pfs)

# Test proportional hazards for PFS
ph_test_pfs <- cox.zph(cox_pfs)
ph_test_pfs
