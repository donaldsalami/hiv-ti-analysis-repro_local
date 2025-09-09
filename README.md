# HIV Treatment Interruption Analysis — Reproducible Repository (Synthetic Data)

This repository organizes R scripts used to analyze first and recurrent ART treatment interruptions,
with **fully synthetic (dummy) data** so results are reproducible without sharing identifiable data.

## What’s inside

- `R/`
  - `01_first_interruption_time_invariant.R`: Aalen–Johansen, CIF plots, facility-frailty Cox
  - `02_multistate_time_varying.R`: multistate Cox with `strata(trans)` + `cluster(patientUID)`, `msfit()` → `probtrans()`
  - `03_mediation_time_varying.R`: discrete-time pooled model, DSD×MMD interaction, CMAverse natural effects + CDEs (visit-based)
  - `04_mediation_time_invariant.R`: baseline Cox mediation (CMAverse), ordinal MMD mediator
  - `99_session_info.R`: session info for reproducibility
- `data-raw/`
  - `generate_dummy_data.R`: now aligned with all factor levels used by mediation scripts
- `data/`
  - `Time_invariant_df_Final_4Aug25`
  - `Time_dependent_iit_df_Final_4Aug25.csv`
  - `README.md`
- `figures/`, `outputs/`: created by scripts

## Usage

1. Run `source("data-raw/generate_dummy_data.R")`.
2. Run scripts in order: `R/01_first_interruption_time_invariant.R`, `R/02_multistate_time_varying.R`.
3. For mediation:
   - Visit-based (time-varying): `R/03_mediation_time_varying.R`
   - Baseline (time-invariant): `R/04_mediation_time_invariant.R`

> Note: CMAverse is required for mediation (`install.packages("CMAverse")`). For parallel runs, also install `future` and `furrr`.
