# R/01_first_interruption_time_invariant.R
# Cleaned from user's script — uses relative paths and exports to figures/ and outputs/

# ---- Setup ----
source("R/_globals.R")

suppressPackageStartupMessages({
  library(survival); library(mstate); library(readr); library(dplyr); library(lubridate)
  library(ggplot2); library(tidyr); library(stringr); library(purrr); library(broom)
  library(gtsummary)
})

# ---- Data ----
infile <- file.path(DATA_DIR, "Time_invariant_df_Final_4Aug25.csv")
if (!file.exists(infile)) {
  stop("Missing input: ", infile, "\nRun data-raw/generate_dummy_data.R first.")
}
Time_invariant_df <- readr::read_csv(infile, show_col_types = FALSE)

time_invariant_df <- Time_invariant_df %>%
  select(patientUID, SiteCode, status, tstop,
         AgeARTStartCategory, Gender, EducationLevel, MaritalStatus, Occupation, NCD,
         Adherence, WHOStage, ViralLoadCategory, StabilityAssessment, BaselineRegimen,
         RegimenLine, MMD, DifferentiatedCare, PwP_countCategory, YearStartART,
         CovidLockDown, MissedAppointment, DefaultAppointment)

# Factors and reference levels
facs <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation","NCD",
          "Adherence","WHOStage","ViralLoadCategory","StabilityAssessment","BaselineRegimen",
          "RegimenLine","MMD","DifferentiatedCare","PwP_countCategory","YearStartART",
          "CovidLockDown","MissedAppointment","DefaultAppointment")

time_invariant_df <- time_invariant_df %>% mutate(across(all_of(facs), as.factor))

time_invariant_df <- within(time_invariant_df, {
  AgeARTStartCategory  <- relevel(AgeARTStartCategory, ref = "25-34")
  Gender               <- relevel(Gender, ref = "Female")
  EducationLevel       <- relevel(EducationLevel, ref = "None")
  MaritalStatus        <- relevel(MaritalStatus, ref = "Married Monogamous")
  Occupation           <- relevel(Occupation, ref = "Unknown")
  NCD                  <- relevel(NCD, ref = "0")
  Adherence            <- relevel(Adherence, ref = "Good")
  WHOStage             <- relevel(WHOStage, ref = "1")
  ViralLoadCategory    <- relevel(ViralLoadCategory, ref = "<200")
  StabilityAssessment  <- relevel(StabilityAssessment, ref = "Stable")
  BaselineRegimen      <- relevel(BaselineRegimen, ref = "2NRTIs + INSTI")
  RegimenLine          <- relevel(RegimenLine, ref = "FirstLine")
  MMD                  <- relevel(MMD, ref = "<3 Months")
  DifferentiatedCare   <- relevel(DifferentiatedCare, ref = "Standard Care")
  PwP_countCategory    <- relevel(PwP_countCategory, ref = "0")
  YearStartART         <- relevel(YearStartART, ref = "2017")
  CovidLockDown        <- relevel(CovidLockDown, ref = "0")
  MissedAppointment    <- relevel(MissedAppointment, ref = "0")
  DefaultAppointment   <- relevel(DefaultAppointment, ref = "0")
})

# time in months + event factor
time_invariant_df <- time_invariant_df %>% mutate(time = tstop/30)
time_invariant_df$event <- factor(time_invariant_df$status, levels = 0:2,
                                  labels = c("censor","interruption","competing"))

# Aalen–Johansen via survfit(Surv(time,event))
sf <- survfit(Surv(time, event) ~ 1, data = time_invariant_df, id = time_invariant_df$patientUID)

# Probability-in-state stacked plot
sfsum <- summary(sf)
states <- sfsum$states
ord <- match(c("(s0)", "interruption", "competing"), states)
time <- c(0, sfsum$time)
P <- rbind(c(1,0,0), as.matrix(sfsum$pstate)[, ord, drop = FALSE])
colnames(P) <- c("(s0)", "interruption", "competing")
p_s0  <- P[, "(s0)"]; p_int <- P[, "interruption"]; p_cmp <- P[, "competing"]

png(file.path(FIG_DIR, "TimeInvariant_Stacked_probability_plot.png"), width=900, height=650, res=120)
plot(NA, xlim = range(time), ylim = c(0, 1), xaxs = "i", yaxs = "i",
     xlab = "Months since ART start", ylab = "Probability in state",
     main = "Probability in State Over Time (Aalen–Johansen)")
for (i in 1:(length(time)-1)) {
  t0 <- time[i]; t1 <- time[i+1]
  y0 <- 0; y1 <- p_s0[i+1]; rect(t0,y0,t1,y1,col="#009E73", border=NA)
  y0 <- y1; y1 <- y1 + p_int[i+1]; rect(t0,y0,t1,y1,col="#e34a33", border=NA)
  y0 <- y1; y1 <- 1; rect(t0,y0,t1,y1,col="#1A85FF", border=NA)
}
box(); legend("right", fill=c("#009E73","#e34a33","#1A85FF"),
              legend=c("(s0) at risk","interruption","competing"), bty="n")
dev.off()

# CIF-style curves
png(file.path(FIG_DIR, "TimeInvariant_CIF_curves.png"), width=900, height=650, res=120)
plot(NA, xlim = range(time), ylim = c(0, 1), xaxs = "i", yaxs = "i",
     xlab = "Months since ART start", ylab = "Cumulative incidence",
     main = "Cumulative Incidence by Event Type (Aalen–Johansen)")
lines(time, p_int, type="s", lwd=2, col="#e34a33")
lines(time, p_cmp, type="s", lwd=2, col="#1A85FF")
legend("right", col=c("#e34a33","#1A85FF"), lwd=2, bty="n",
       legend=c("interruption","competing"))
dev.off()

# Transition probabilities at 1,3,6,12 months
t_req <- c(1,3,6,12)
ss <- summary(sf, times = t_req, extend = TRUE)
state_name <- "interruption"
stopifnot(state_name %in% colnames(ss$pstate))
res <- data.frame(month = ss$time, p_interruption = ss$pstate[, state_name])
if (!is.null(ss$lower) && !is.null(ss$upper)) {
  res$lower <- ss$lower[, state_name]; res$upper <- ss$upper[, state_name]
}
res$`p_interruption(%)` <- round(100 * res$p_interruption, 1)
if ("lower" %in% names(res)) res$`lower(%)` <- round(100 * res$lower, 1)
if ("upper" %in% names(res)) res$`upper(%)` <- round(100 * res$upper, 1)

readr::write_csv(res, file.path(OUT_DIR, "TimeInvariant_prob_interruption_1_3_6_12m.csv"))

# Cox model with facility-level frailty
covariates <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation","NCD",
                "Adherence","WHOStage","ViralLoadCategory","StabilityAssessment","BaselineRegimen",
                "RegimenLine","MMD","DifferentiatedCare","PwP_countCategory","YearStartART",
                "CovidLockDown","MissedAppointment","DefaultAppointment")

cox_formula <- as.formula(paste("Surv(time, status == 1) ~",
                                paste(c(covariates, "frailty(SiteCode)"), collapse=" + ")))
cox_frailty <- coxph(cox_formula, data = time_invariant_df, method="efron")
zph <- cox.zph(cox_frailty)

# Export table
tbl <- tbl_regression(cox_frailty, exponentiate=TRUE)
df_tbl <- as.data.frame(tbl) %>% mutate(`HR(95% CI)` = paste0(`**HR**`, " (", `**95% CI**`, ")"))
readr::write_csv(df_tbl, file.path(OUT_DIR, "Time_invariant_Cox_Model.csv"))