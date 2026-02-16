##################################################
## Project: MSc HDS Project
## Date: 13 February 2026
## Author: Autumn Johnson
##################################################

## Install Package (Only once if required)
##################################################
# install.packages("WinRatio")
# install.packages("BuyseTest")
# install.packages("msm")

## Section: 0.0 Load Libraries
##################################################
library(readxl)
library(dplyr)
library(table1)
library(survival)
library(survminer)
library(WinRatio) # For Win ratio
library(BuyseTest) # For Net Benefit
library(msm) # For Multi-state model

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
                 labels = c("Wild Type (Normal)", "Mutated")),
    
  ) %>%
  select(id, Arm, Sites, LDH, TMB, JAK, pfs_time, pfs_status, os_time, os_status)


## Section: 5.0    Win Ratio
##################################################

## Section: 5.1    B vs. A
##################################################
wr_ba <- winratio(id = "id",
                  trt = "Arm",
                  outcomes = list(
                    outc1 = c("os_status", "s", "os_time"),
                    outc2 = c("pfs_status", "s", "pfs_time")
                  ),
                  fu = "os_time",
                  data = clean_df %>% filter(Arm != "C") %>%
                    filter(!is.na(pfs_time)| !is.na(os_time)) %>%
                    mutate(pfs_time = round(pfs_time, 1),
                           os_time = round(os_time, 1))
)

summary(wr_ba)

## Section: 5.2    C vs. A
##################################################
wr_ca <- winratio(id = "id",
                  trt = "Arm",
                  outcomes = list(
                    outc1 = c("os_status", "s", "os_time"),
                    outc2 = c("pfs_status", "s", "pfs_time")
                  ),
                  fu = "os_time",
                  data = clean_df %>% filter(Arm != "B") %>%
                    filter(!is.na(pfs_time)| !is.na(os_time)) %>%
                    mutate(pfs_time = round(pfs_time, 1),
                           os_time = round(os_time, 1))
)

summary(wr_ca)

## Section: 5.3    C vs. B
##################################################
wr_cb <- winratio(id = "id",
                  trt = "Arm",
                  outcomes = list(
                    outc1 = c("os_status", "s", "os_time"),
                    outc2 = c("pfs_status", "s", "pfs_time")
                  ),
                  fu = "os_time",
                  data = clean_df %>% filter(Arm != "A") %>%
                    filter(!is.na(pfs_time)| !is.na(os_time)) %>%
                    mutate(pfs_time = round(pfs_time, 1),
                           os_time = round(os_time, 1))
)

summary(wr_cb)

# Adjusted p-values for multiple testing
pvals <- c(
  "B_vs_A" = wr_ab$p.value,
  "C_vs_B" = wr_bc$p.value,
  "C_vs_A" = wr_ac$p.value
)

p.adjust(pvals, method = "bonferroni")

## Section: 6.0    Net Benefit
##################################################

## Section: 6.1    B vs. A
##################################################
nb_df_BA <- clean_df %>%
  filter(Arm %in% c("A","B")) %>%
  filter(!is.na(os_time) | !is.na(pfs_time)) %>%   # keep if at least one endpoint present
  mutate(
    trt = factor(Arm, levels = c("A","B"))
  )

nb_BA <- BuyseTest(
  trt ~ 
    tte(os_time, os_status, threshold = 0) +   # priority 1: OS
    tte(pfs_time, pfs_status, threshold = 0),  # priority 2: PFS
  data = nb_df_BA,
  hierarchical = TRUE
)

summary(nb_BA)

## Section: 6.2    C vs. A
##################################################
nb_df_CA <- clean_df %>%
  filter(Arm %in% c("A","C")) %>%
  filter(!is.na(os_time) | !is.na(pfs_time)) %>%
  mutate(trt = factor(Arm, levels = c("A","C")))   # A control, C treatment

nb_CA <- BuyseTest(
  trt ~ 
    tte(os_time, os_status, threshold = 0) +
    tte(pfs_time, pfs_status, threshold = 0),
  data = nb_df_CA,
  hierarchical = TRUE
)

summary(nb_CA)

## Section: 6.3    C vs. B
##################################################
nb_df_CB <- clean_df %>%
  filter(Arm %in% c("B","C")) %>%
  filter(!is.na(os_time) | !is.na(pfs_time)) %>%
  mutate(trt = factor(Arm, levels = c("B","C")))   # B control, C treatment

nb_CB <- BuyseTest(
  trt ~ 
    tte(os_time, os_status, threshold = 0) +
    tte(pfs_time, pfs_status, threshold = 0),
  data = nb_df_CB,
  hierarchical = TRUE
)

summary(nb_CB)


## Section: 7.0    Multi-state model
##################################################

## Section: 7.1    Define states
##################################################
# State 1: Stable
# State 2: Progression
# State 3: Death

states <- c("Stable", "Progression", "Death")

## Section: 7.2    Define transition intensity matrix (Q)
##################################################
Q <- matrix(0, nrow = 3, ncol = 3, dimnames = list(states, states))

Q["Stable","Progression"] <- 0.1
Q["Stable","Death"]       <- 0.05
Q["Progression","Death"]  <- 0.1
# Other transitions remain 0 (impossible)


## Section: 7.3    Prepare data in msm long format
##################################################
msm_df <- clean_df %>%
  select(id, Arm, os_time, os_status, pfs_time, pfs_status) %>%
  rowwise() %>%
  do({
    id <- .$id
    arm <- .$Arm
    
    # Start in Stable at time 0
    df_list <- list(
      data.frame(id = id, time = 0, state = 1, Arm = arm)
    )
    
    # Progression, only if before death
    if(!is.na(.$pfs_time) && .$pfs_status == 1 && 
       (is.na(.$os_time) || .$pfs_time < .$os_time)){
      df_list <- c(df_list, list(data.frame(id = id, time = .$pfs_time, state = 2, Arm = arm)))
    }
    
    # Death, if occurred
    if(!is.na(.$os_time) && .$os_status == 1){
      df_list <- c(df_list, list(data.frame(id = id, time = .$os_time, state = 3, Arm = arm)))
    }
    
    bind_rows(df_list)
  }) %>%
  ungroup() %>%
  arrange(id, time)


## Section: 7.4    Multi-state  model
##################################################
# Arm as covariate affecting all transitions
msm_fit <- msm(
  state ~ time,
  subject = id,
  data = msm_df,
  qmatrix = Q,
  covariates = ~ Arm,
  exacttimes = TRUE
)

# Summary
summary(msm_fit)

## Section: 7.5    Transition probabilities
##################################################
# Probabilities from Stable at time t
pt <- pmatrix.msm(msm_fit, t = 12)   # e.g., 12 months
pt

# Example: plot transition probabilities from Stable over time
times <- seq(0, 24, by = 1)  # 0–24 months
probs <- sapply(times, function(t) pmatrix.msm(msm_fit, t = t)[1, ])
matplot(times, t(probs), type = "l", lty = 1, col = 1:3,
        xlab = "Time (months)", ylab = "Probability",
        main = "State occupation probabilities from Stable")
legend("topright", legend = states, col = 1:3, lty = 1)
