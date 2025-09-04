# R/02_multistate_time_varying.R
# Cleaned from user's script — relative paths, consistent outputs

source("R/_globals.R")
suppressPackageStartupMessages({
  library(survival); library(mstate); library(readr); library(dplyr)
  library(tidyr); library(purrr); library(stringr); library(ggplot2)
})

infile <- file.path(DATA_DIR, "Time_dependent_iit_df_Final_4Aug25.csv")
if (!file.exists(infile)) stop("Missing input: ", infile, "\nRun data-raw/generate_dummy_data.R first.")
time_varying_df <- readr::read_csv(infile, show_col_types = FALSE)

# Map states to numeric
state_map <- c("State1"=1, "State2"=2, "State3"=3, "State4"=4, "State5"=5)
time_varying_df <- time_varying_df %>%
  mutate(State_from = as.integer(state_map[State_from]),
         State_to   = as.integer(state_map[State_to]))

# Transition keys (plot set)
transition_map <- list("1_2"=1,"2_3"=2,"2_4"=3,"2_5"=4,"3_4"=5,"3_5"=6,"4_3"=7,"4_5"=8)

mutistate_df <- time_varying_df %>%
  mutate(trans_key = paste(State_from, State_to, sep="_"),
         trans = sapply(trans_key, function(k) { if (!is.null(transition_map[[k]])) transition_map[[k]] else NA_integer_ }),
         status = ifelse(!is.na(trans) & status==0, 1, status))

# Factors + references
facs <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation","NCD",
          "Adherence","WHOStage","ViralLoadCategory","StabilityAssessment","BaselineRegimen",
          "RegimenLine","MMD","DifferentiatedCare","PwP_countCategory","YearStartART",
          "CovidLockDown","MissedAppointmentCumCategory","DefaultAppointmentCumCategory","TimeReturnToCareCategory")

mutistate_df <- mutistate_df %>% mutate(across(all_of(facs), as.factor))

mutistate_df <- within(mutistate_df, {
  AgeARTStartCategory  <- relevel(AgeARTStartCategory, ref="0-14")
  Gender               <- relevel(Gender, ref="Female")
  EducationLevel       <- relevel(EducationLevel, ref="Tertiary")
  MaritalStatus        <- relevel(MaritalStatus, ref="Other")
  Occupation           <- relevel(Occupation, ref="Employed")
  NCD                  <- relevel(NCD, ref="0")
  Adherence            <- relevel(Adherence, ref="Good")
  WHOStage             <- relevel(WHOStage, ref="1")
  ViralLoadCategory    <- relevel(ViralLoadCategory, ref="<200")
  StabilityAssessment  <- relevel(StabilityAssessment, ref="Stable")
  BaselineRegimen      <- relevel(BaselineRegimen, ref="Other")
  RegimenLine          <- relevel(RegimenLine, ref="Other")
  MMD                  <- relevel(MMD, ref="<3 Months")
  DifferentiatedCare   <- relevel(DifferentiatedCare, ref="Standard Care")
  PwP_countCategory    <- relevel(PwP_countCategory, ref="0")
  YearStartART         <- relevel(YearStartART, ref="2017")
  CovidLockDown        <- relevel(CovidLockDown, ref="0")
  MissedAppointmentCumCategory  <- relevel(MissedAppointmentCumCategory, ref="0")
  DefaultAppointmentCumCategory <- relevel(DefaultAppointmentCumCategory, ref="0")
  TimeReturnToCareCategory      <- relevel(TimeReturnToCareCategory, ref="0 days")
})

mutistate_df$trans <- as.factor(mutistate_df$trans)
mutistate_df$patientUID <- as.factor(mutistate_df$patientUID)
mutistate_df$SiteCode <- as.factor(mutistate_df$SiteCode)

ms_covariates <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation","NCD",
                   "Adherence","WHOStage","ViralLoadCategory","StabilityAssessment","BaselineRegimen",
                   "RegimenLine","MMD","DifferentiatedCare","PwP_countCategory","YearStartART",
                   "CovidLockDown","MissedAppointmentCumCategory","DefaultAppointmentCumCategory","TimeReturnToCareCategory")

# Global multistate Cox with cluster + strata
cox_frailty_ms <- coxph(
  as.formula(paste("Surv(tstart, tstop, status) ~",
                   paste(c(ms_covariates,"cluster(patientUID)","strata(trans)"), collapse=" + "))),
  data = mutistate_df, method="efron", model=TRUE
)

# Transitions structure (no absorbing from A to D directly in this plotting set)
tmat <- transMat(
  x = list(
    c(2),      # A -> B
    c(3,4,5),  # B -> R,C,D
    c(4,5),    # R -> C,D
    c(3,5),    # C -> R,D
    c()        # D absorbing
  ),
  names = c("A","B","R","C","D")
)

# Build representative newdata (sample per trans)
build_newdata_for_msfit <- function(data, model, n_per_trans = 5) {
  trans_levels <- levels(model.frame(model)$`strata(trans)`)
  data <- data %>% filter(!is.na(trans)) %>% mutate(trans = factor(trans, levels = trans_levels))
  newdata <- data %>% group_by(trans) %>% slice_sample(n = n_per_trans, replace = TRUE) %>% ungroup()
  newdata$trans <- factor(newdata$trans, levels = trans_levels)
  newdata$strata <- factor(as.character(newdata$trans), levels = trans_levels)
  newdata
}

template_rows <- build_newdata_for_msfit(mutistate_df, cox_frailty_ms, n_per_trans = 5)

# msfit + probtrans
msf <- msfit(cox_frailty_ms, trans=tmat, newdata=template_rows, variance=TRUE)
pt  <- probtrans(msf, predt = 0)

# Plot state occupation from A
png(file.path(FIG_DIR, "MultiState_State_Occupation_Probabilities.png"), width=900, height=650, res=120)
plot(pt, from=1, type="filled",
     xlab="Time (months)", ylab="Probability",
     col=c("green","red","blue","purple","gray"), lwd=2, xaxt="n",
     main="State Occupation Probabilities")
time_days <- pt[[1]]$time
month_ticks <- pretty(time_days/30.44)
axis(1, at = month_ticks*30.44, labels = month_ticks)
legend("topright", legend=c("Treatment Initiated (A)","First Interruption (B)",
                            "Re-engaged (R)","Subsequent Interruption (C)","Discontinued (D)"),
       fill=c("green","red","blue","purple","gray"), border="black", bty="n", cex=0.8)
dev.off()

# Extract state occupation at 1,3,6,12,24 months
months <- c(1,3,6,12,24)
time_points <- months * 30.44
get_probs <- function(pt, t) {
  idx <- which.min(abs(pt[[1]]$time - t))
  probs <- pt[[1]][idx, grepl("^pstate", names(pt[[1]]))]
  as.numeric(probs)
}
prob_matrix <- t(sapply(time_points, function(t) get_probs(pt, t)))
rownames(prob_matrix) <- paste0(months, " months")
colnames(prob_matrix) <- c("A","B","R","C","D")

readr::write_csv(as.data.frame(prob_matrix), file.path(OUT_DIR, "StateOccupationProbabilities_1_3_6_12_24.csv"))

# By MMD profiles
mmd_levels <- levels(mutistate_df$MMD)
mmd_profiles <- do.call("rbind", lapply(mmd_levels, function(m) {x <- template_rows; x$MMD <- m; x$Profile <- m; x}))
mmd_profiles$strata <- factor(as.character(mmd_profiles$trans), levels = levels(model.frame(cox_frailty_ms)$`strata(trans)`))

pt_list <- list()
for (m in mmd_levels) {
  this_data <- dplyr::filter(mmd_profiles, Profile == m)
  msf_m <- msfit(cox_frailty_ms, trans=tmat, newdata=this_data, variance=TRUE)
  pt_list[[m]] <- probtrans(msf_m, predt=0)
}

# Summaries at 6,12,24 months with CIs if available
extract_summary <- function(pt, group_name, group_type, time_points) {
  df <- as.data.frame(pt[[1]]); df$time_months <- df$time / 30.44
  prob_cols <- grep("^pstate", names(df), value = TRUE)
  se_cols <- paste0("se", gsub("pstate", "", prob_cols))
  if (!all(se_cols %in% names(df))) return(NULL)
  out <- list()
  for (t in time_points) {
    row <- df[which.min(abs(df$time_months - t)), , drop = FALSE]
    probs <- as.numeric(row[prob_cols]); ses <- as.numeric(row[se_cols])
    lowers <- pmax(0, probs - 1.96*ses); uppers <- pmin(1, probs + 1.96*ses)
    out[[as.character(t)]] <- data.frame(GroupType = group_type, Profile = group_name,
                                         TimeMonths = round(row$time/30.44, 1),
                                         State = paste0("State ", seq_along(probs)),
                                         Mean = round(probs, 3),
                                         Lower95CI = round(lowers, 3),
                                         Upper95CI = round(uppers, 3))
  }
  dplyr::bind_rows(out)
}

target_months <- c(6,12,24)
mmd_summary <- dplyr::bind_rows(lapply(mmd_levels, function(m) extract_summary(pt_list[[m]], m, "MMD", target_months)))
readr::write_csv(mmd_summary, file.path(OUT_DIR, "Transition_Probabilities_by_MMD_6_12_24.csv"))

# By Differentiated Care profiles
dcare_levels <- levels(mutistate_df$DifferentiatedCare)
dcare_profiles <- do.call("rbind", lapply(dcare_levels, function(d) {x <- template_rows; x$DifferentiatedCare <- d; x$Profile <- d; x}))
dcare_profiles$strata <- factor(as.character(dcare_profiles$trans), levels = levels(model.frame(cox_frailty_ms)$`strata(trans)`))

pt_list_dcare <- list()
for (d in dcare_levels) {
  this_data <- dplyr::filter(dcare_profiles, Profile == d)
  msf_d <- msfit(cox_frailty_ms, trans=tmat, newdata=this_data, variance=TRUE)
  pt_list_dcare[[d]] <- probtrans(msf_d, predt=0)
}
dcare_summary <- dplyr::bind_rows(lapply(dcare_levels, function(d) extract_summary(pt_list_dcare[[d]], d, "DifferentiatedCare", target_months)))
readr::write_csv(dcare_summary, file.path(OUT_DIR, "Transition_Probabilities_by_DSD_6_12_24.csv"))

# Transition-specific models (per trans) with minimal row threshold
fit_transition_model <- function(k) {
  df_k <- dplyr::filter(mutistate_df, trans == k)
  if (nrow(df_k) < 50) return(NULL)
  formula_k <- as.formula(paste("Surv(tstart, tstop, status) ~", paste(ms_covariates, collapse = " + ")))
  survival::coxph(formula_k, data = df_k, method = "efron")
}
models_by_trans <- purrr::map(1:8, fit_transition_model)
names(models_by_trans) <- paste0("trans", 1:8)

extract_hr_table <- function(model, transition_id) {
  if (is.null(model)) return(NULL)
  est <- summary(model)$coefficients; conf <- summary(model)$conf.int
  tibble::tibble(Transition = transition_id, Covariate = rownames(est),
                 HR = round(conf[, "exp(coef)"], 2),
                 Lower95CI = round(conf[, "lower .95"], 2),
                 Upper95CI = round(conf[, "upper .95"], 2),
                 p_value = round(est[, "Pr(>|z|)"], 3),
                 `HR (95% CI)` = paste0(round(conf[, "exp(coef)"], 2), " (",
                                        round(conf[, "lower .95"], 2), ", ",
                                        round(conf[, "upper .95"], 2), ")"))
}
results_hr <- purrr::map2_dfr(models_by_trans, names(models_by_trans), extract_hr_table)
readr::write_csv(results_hr, file.path(OUT_DIR, "Multistate_Transition_Models.csv"))