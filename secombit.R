##################################################
## Project: MSc HDS Project
## Date: 26 January 2026
## Author: Autumn O'Donnell
##################################################

## Section: 0.0 Load Libraries
##################################################
library(dplyr)
library(table1)
library(survival)
library(survminer)

## Section: 0.1 Load Dataset
##################################################
SECOMBIT <- read.csv("Secombit data.csv", check.names = FALSE)
str(SECOMBIT)

## Helper: Convert common status encodings to 0/1
to_event01 <- function(x) {
  # convert factors safely
  if (is.factor(x)) x <- as.character(x)
  x_chr <- trimws(tolower(as.character(x)))
  
  # if already numeric-like, coerce
  suppressWarnings({
    x_num <- as.numeric(x_chr)
  })
  
  # If numeric conversion worked for most values, use it
  if (mean(is.na(x_num)) < 0.5) {
    # standardize: anything >0 becomes 1
    return(ifelse(is.na(x_num), NA, ifelse(x_num > 0, 1, 0)))
  }
  
  # Otherwise map common labels
  event_labels    <- c("1", "event", "dead", "death", "progressed", "progression", "yes", "true")
  censored_labels <- c("0", "censored", "alive", "no", "false")
  
  out <- ifelse(x_chr %in% event_labels, 1,
                ifelse(x_chr %in% censored_labels, 0, NA))
  return(out)
}

## Section: 1.0 Clean Data
##################################################
clean_df <- SECOMBIT %>%
  transmute(
    # ID
    id = ordine,
    
    # Treatment arm (A baseline)
    Arm = factor(ARM, levels = c("A", "B", "C")),
    
    # Sites (1–2 baseline)
    Sites = factor(sites, levels = c("1-2", ">=3"),
                   labels = c("1 - 2", ">= 3")),
    
    # LDH (Normal baseline)
    LDH = factor(ULN_LDH, levels = c("normal", "elevated"),
                 labels = c("Normal", "Elevated")),
    
    # TMB (<10 baseline)
    TMB = factor(TMB, levels = c("<10", ">=10"),
                 labels = c("< 10", ">= 10")),
    
    # JAK (wt baseline)
    JAK = factor(JAK, levels = c("wt", "mut"),
                 labels = c("Wild Type (Normal)", "Mutated")),
    
    # Times (force numeric)
    os_time  = suppressWarnings(as.numeric(OS)),
    pfs_time = suppressWarnings(as.numeric(`PFS TOTAL`)),
    
    # Status (force 0/1)
    os_status  = to_event01(status),
    pfs_status = to_event01(`Progr total`)
  ) %>%
  # Optional: drop rows missing key survival inputs
  filter(!is.na(os_time), !is.na(os_status),
         !is.na(pfs_time), !is.na(pfs_status))

## Quick sanity checks (highly recommended)
table(clean_df$os_status, useNA = "ifany")
table(clean_df$pfs_status, useNA = "ifany")
summary(clean_df$os_time)
summary(clean_df$pfs_time)

## Section: 2.0 Summary table with table1
##################################################
# nicer labels in table1
label(clean_df$Sites) <- "Sites"
label(clean_df$LDH)   <- "LDH"
label(clean_df$TMB)   <- "TMB"
label(clean_df$JAK)   <- "JAK"

table1(~ Sites + LDH + TMB + JAK | Arm,
       data = clean_df,
       overall = FALSE)

## Section: 3.0 Kaplan-Meier plots
##################################################
fit_os <- survfit(Surv(os_time, os_status) ~ Arm, data = clean_df)

ggsurvplot(
  fit_os, data = clean_df,
  xlab = "Time (months)",
  ylab = "Survival probability",
  title = "Overall Survival (OS)",
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE
)

fit_pfs <- survfit(Surv(pfs_time, pfs_status) ~ Arm, data = clean_df)

ggsurvplot(
  fit_pfs, data = clean_df,
  xlab = "Time (months)",
  ylab = "Survival probability",
  title = "Progression-Free Survival (PFS)",
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE
)

## Section: 3.1 Log-rank tests
##################################################
logrank_os  <- survdiff(Surv(os_time, os_status) ~ Arm, data = clean_df)
logrank_pfs <- survdiff(Surv(pfs_time, pfs_status) ~ Arm, data = clean_df)

logrank_os
logrank_pfs

## Section: 4.0 Cox model
##################################################
cox_os <- coxph(Surv(os_time, os_status) ~ Arm + Sites + LDH + TMB + JAK,
                data = clean_df, ties = "efron")
summary(cox_os)

ph_test_os <- cox.zph(cox_os)
ph_test_os
# plot(ph_test_os)  # uncomment to visualize

cox_pfs <- coxph(Surv(pfs_time, pfs_status) ~ Arm + Sites + LDH + TMB + JAK,
                 data = clean_df, ties = "efron")
summary(cox_pfs)

ph_test_pfs <- cox.zph(cox_pfs)
ph_test_pfs
# plot(ph_test_pfs)

