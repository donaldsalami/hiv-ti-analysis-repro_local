# 03_models_time_varying.R
# Discrete-time pooled logistic model for interval risk (time-varying)

library(readr)
library(dplyr)
library(stats)

duration_df <- read_csv("data/duration_df.csv", show_col_types = FALSE)

# outcome is 'status' indicating interruption in the interval (0/1)
# logit with lagged MMD and log interval length as covariates
duration_df <- duration_df %>% mutate(log_interval = log(pmax(TimeReturnToCare, 1)))

glm1 <- glm(status ~ lag_MMD + AgeARTStartCategory + Gender + DSD + Adherence + WHOStage +
              ViralLoadCategory + log_interval,
            data = duration_df,
            family = binomial())

summary(glm1)
saveRDS(glm1, "data/pooled_logistic_timevarying.rds")