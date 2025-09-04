# 02_models_time_invariant.R
# Cause-specific Cox model for first interruption (with competing risks represented in event codes)

library(readr)
library(dplyr)
library(survival)

time_invariant_df <- read_csv("data/time_invariant_df.csv", show_col_types = FALSE)

# Event coding: 1=interruption, 2=competing (death/LTFU/transfer)
# Fit cause-specific Cox for interruption (event==1)
time_invariant_df <- time_invariant_df %>%
  mutate(event_interruption = ifelse(event == 1, 1, 0))

cox1 <- coxph(Surv(time, event_interruption) ~ AgeARTStartCategory + Gender + DSD + MMD +
                BaselineCD4Category + YearStartART + Adherence + WHOStage + NCD,
              data = time_invariant_df)

summary(cox1)
saveRDS(cox1, "data/cox_first_interrupt.rds")