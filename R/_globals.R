# R/_globals.R
# Helper paths and options
options(stringsAsFactors = FALSE)
DATA_DIR    <- "data"
FIG_DIR     <- "figures"
OUT_DIR     <- "outputs"

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)