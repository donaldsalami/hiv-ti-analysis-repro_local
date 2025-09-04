# 05_figures_tables.R
# Example plots/tables (base R and ggplot2 optional)

library(readr)
library(dplyr)
library(ggplot2)

time_invariant_df <- read_csv("data/time_invariant_df.csv", show_col_types = FALSE)

# Kaplan-Meier for any-event endpoint (for illustration)
library(survival)
km <- survfit(Surv(time, event==1) ~ DSD, data = time_invariant_df)

# Quick ggplot of median times by DSD
summary_df <- time_invariant_df %>%
  group_by(DSD) %>%
  summarise(mean_time = mean(time), .groups="drop")

p <- ggplot(summary_df, aes(x=DSD, y=mean_time)) + geom_col() + labs(y="Mean time to event (days)")
ggsave("data/fig_mean_time_by_DSD.png", width = 6, height = 4, dpi = 120)