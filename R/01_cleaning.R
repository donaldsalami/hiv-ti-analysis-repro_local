# 01_cleaning.R
# Example cleaning for synthetic data produced by data-raw/generate_dummy_data.R

library(dplyr)
library(readr)

sub_df <- read_csv("data/sub_df.csv", show_col_types = FALSE)

# Basic recoding examples
sub_df <- sub_df %>%
  mutate(
    Gender = factor(Gender, levels = c("Male","Female")),
    DSD = factor(DSD, levels = c("StandardCare","FastTrack","CommunityART")),
    MMD = factor(MMD, ordered = TRUE, levels = c("0-2","3-5","6+"))
  )

write_csv(sub_df, "data/sub_df_clean.csv")
message("Saved: data/sub_df_clean.csv")