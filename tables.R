library(table1)
library(here)

df <- read.csv("Secombit data.csv")

summary_stat <- table1(~ sites + ULN_LDH + TMB + JAK | ARM, data = df)


