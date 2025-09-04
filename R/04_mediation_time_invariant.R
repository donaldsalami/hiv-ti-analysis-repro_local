# R/04_mediation_time_invariant.R
# Baseline mediation with Cox outcome and ordinal MMD mediator (CMAverse)

source("R/_globals.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(purrr); library(forcats); library(CMAverse); library(future); library(furrr)
})

infile <- file.path(DATA_DIR, "Baseline_iit_df_Final_4Aug25.csv")
stopifnot(file.exists(infile))
Baseline_iit_df <- readr::read_csv(infile, show_col_types = FALSE)

df <- Baseline_iit_df %>%
  mutate(iit_event = if_else(as.numeric(status)==1, 1L, 0L))

df <- df %>%
  select(patientUID, iit_event, ART_durationMths,
         AgeARTStartCategory, Gender, EducationLevel, MaritalStatus, Occupation, NCD,
         Adherence, WHOStage, ViralLoadCategory, StabilityAssessment, BaselineRegimen,
         RegimenLine, MMD, DifferentiatedCare, PwP_countCategory, YearStartART,
         CovidLockDown, MissedAppointment, DefaultAppointment)

fac_vars <- c("MissedAppointment","DefaultAppointment","MaritalStatus","RegimenLine",
              "DifferentiatedCare","StabilityAssessment","Adherence","BaselineRegimen",
              "AgeARTStartCategory","Gender","EducationLevel","YearStartART","WHOStage",
              "ViralLoadCategory","MMD","CovidLockDown","PwP_countCategory","NCD")
df <- df %>% mutate(across(all_of(fac_vars), \(x) if (!is.factor(x)) factor(x) else x))

df <- within(df, {
  MissedAppointment   <- relevel(MissedAppointment,   ref = "0")
  DefaultAppointment  <- relevel(DefaultAppointment,  ref = "0")
  MaritalStatus       <- relevel(MaritalStatus,       ref = "Other")
  RegimenLine         <- relevel(RegimenLine,         ref = "Other")
  BaselineRegimen     <- relevel(BaselineRegimen,     ref = "Other")
  DifferentiatedCare  <- relevel(DifferentiatedCare,  ref = "Standard Care")
  StabilityAssessment <- relevel(StabilityAssessment, ref = "Stable")
  Adherence           <- relevel(Adherence,           ref = "Good")
  EducationLevel      <- relevel(EducationLevel,      ref = "Tertiary")
  AgeARTStartCategory <- relevel(AgeARTStartCategory, ref = "0-14")
  YearStartART        <- relevel(YearStartART,        ref = "2017")
  PwP_countCategory   <- relevel(PwP_countCategory,   ref = "0")
  WHOStage            <- relevel(WHOStage,            ref = "1")
  NCD                 <- relevel(NCD,                 ref = "0")
  Gender              <- relevel(Gender,              ref = "Female")
  ViralLoadCategory   <- relevel(ViralLoadCategory,   ref = "<200")
  MMD                 <- relevel(MMD,                 ref = "<3 Months")
  CovidLockDown       <- relevel(CovidLockDown,       ref = "0")
})

# Ordered factors and menus
df$AgeARTStartCategory <- factor(df$AgeARTStartCategory, levels=c("0-14","15-24","25-34","35-44","45-59","60+"), ordered=TRUE)
df$WHOStage            <- factor(df$WHOStage, levels=c("1","2","3","4"), ordered=TRUE)
df$YearStartART        <- factor(df$YearStartART, levels=as.character(2017:2023), ordered=TRUE)
df$PwP_countCategory   <- factor(df$PwP_countCategory, levels=c("0","1","2","3 or more"), ordered=TRUE)
df$ViralLoadCategory   <- factor(df$ViralLoadCategory, levels=c("<200","200-999",">1000","Unknown"), ordered=TRUE)
df$StabilityAssessment <- factor(df$StabilityAssessment, levels=c("Stable","Unstable"), ordered=TRUE)
df$NCD                 <- factor(df$NCD, levels=c("0","1"), ordered=TRUE)

# Settings
FAST_MODE   <- FALSE
PARALLEL    <- TRUE
SEED        <- 2025
WORKERS     <- max(1, parallel::detectCores() - 2)
set.seed(SEED)

if (FAST_MODE) { NBOOT_MAIN <- 60; NBOOT_CDE <- 30; BASEC_MODE <- "core" } else { NBOOT_MAIN <- 150; NBOOT_CDE <- 80; BASEC_MODE <- "full" }

MMD_LEVELS <- c("<3 Months","3-5 Months","6+ Months")
DSD_LEVELS <- c("Community ART Distribution HCW Led","Community ART Distribution Peer Led",
                "Facility ART Distribution Group","Fast Track")

basec_full <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation","NCD",
                "Adherence","WHOStage","ViralLoadCategory","StabilityAssessment","BaselineRegimen",
                "RegimenLine","PwP_countCategory","YearStartART","CovidLockDown")
basec_core <- c("AgeARTStartCategory","Gender","WHOStage","ViralLoadCategory","YearStartART")
basec_vars <- if (BASEC_MODE=="core") basec_core else basec_full

df <- df %>% mutate(ART_durationMths = as.numeric(ART_durationMths),
                    iit_event = as.integer(iit_event),
                    DifferentiatedCare = factor(DifferentiatedCare, levels=c("Standard Care", DSD_LEVELS)),
                    MMD = factor(MMD, levels=MMD_LEVELS, ordered=TRUE))

make_subset <- function(one_level) {
  df %>% mutate(A = if_else(DifferentiatedCare==one_level,"DSD",
                            if_else(DifferentiatedCare=="Standard Care","SC",NA_character_)),
                A = factor(A, levels=c("SC","DSD")),
                M = factor(MMD, levels=MMD_LEVELS, ordered=TRUE)) %>%
    filter(!is.na(A)) %>%
    tidyr::drop_na(ART_durationMths, iit_event, A, M, dplyr::all_of(basec_vars))
}

prune_basec <- function(dat, covars, event="iit_event", min_per_cell=3L) {
  keep <- covars[covars %in% names(dat)]
  if (!length(keep)) return(list(dat=dat, basec=character(0)))
  ok <- vapply(keep, function(v) {
    x <- dat[[v]]
    if (!is.factor(x)) {
      vv <- try(stats::var(as.numeric(x), na.rm=TRUE), silent=TRUE)
      return(is.numeric(vv) && is.finite(vv) && !is.na(vv) && vv>0)
    }
    x <- droplevels(as.factor(x)); if (nlevels(x)<2) return(FALSE)
    ev <- dat[[event]]; if (all(is.na(ev))) return(FALSE)
    tab <- try(table(x, ev, useNA="no"), silent=TRUE)
    if (inherits(tab,"try-error")) return(FALSE)
    isTRUE(all(rowSums(tab)>=min_per_cell) && all(colSums(tab)>=min_per_cell))
  }, logical(1))
  keep <- keep[isTRUE(ok)]
  if (!length(keep)) return(list(dat=dat, basec=character(0)))
  dat[keep] <- lapply(dat[keep], function(z) if (is.factor(z)) droplevels(z) else z)
  list(dat=dat, basec=keep)
}

cmest_once <- function(d, basec, mval_scalar, EM=FALSE, nboot=60) {
  CMAverse::cmest(data=d, model="rb", estimation="imputation", inference="bootstrap", nboot=nboot,
                  outcome="ART_durationMths", event="iit_event", yreg="coxph",
                  mediator="M", mreg=list("ordinal"), exposure="A", basec=if(length(basec)) basec else NULL,
                  EMint=EM, astar="SC", a="DSD", mval=list(mval_scalar))
}

safe_fit <- function(expr, dsd_label, extra_cols=list()) {
  out <- tryCatch(expr, error=function(e) e, warning=function(w) w)
  if (inherits(out,"error") || inherits(out,"warning")) {
    return(tibble(dsd_contrast=as.character(dsd_label), effect=NA_character_, est=NA_real_,
                  se=NA_real_, lo=NA_real_, hi=NA_real_, p=NA_real_, .fail=TRUE,
                  error_msg=as.character(conditionMessage(out)), !!!extra_cols))
  }
  tibble(dsd_contrast=as.character(dsd_label), effect=as.character(names(out$effect.pe)),
         est=as.numeric(out$effect.pe), se=as.numeric(out$effect.se),
         lo=as.numeric(out$effect.ci.low), hi=as.numeric(out$effect.ci.high),
         p=as.numeric(out$effect.pval), .fail=FALSE, error_msg=NA_character_, !!!extra_cols)
}
normalize_effect <- function(eff) { sub("^r","", tolower(eff)) }

run_natural <- function(nboot=150, mval_default="<3 Months") {
  purrr::map_dfr(DSD_LEVELS, function(lev) {
    dsub <- make_subset(lev)
    if (nlevels(droplevels(dsub$A)) < 2 || length(unique(dsub$iit_event)) < 2)
      return(tibble(dsd_contrast=lev, kind="natural", effect=NA, est=NA, se=NA, lo=NA, hi=NA, p=NA,
                    .fail=TRUE, error_msg="No exposure or event variation"))
    tab <- with(dsub, table(A,M))
    if (any(tab==0)) return(tibble(dsd_contrast=lev, kind="natural", effect=NA, est=NA, se=NA, lo=NA, hi=NA, p=NA,
                                   .fail=TRUE, error_msg="Empty A×M cell"))
    pr <- prune_basec(dsub, basec_vars, event="iit_event", min_per_cell=3L)
    mlev <- levels(dsub$M); mval_use <- if (mval_default %in% mlev) mval_default else mlev[1]
    safe_fit(cmest_once(pr$dat, pr$basec, mval_scalar=mval_use, EM=FALSE, nboot=nboot),
             dsd_label=lev, extra_cols=list(kind="natural"))
  }) %>% mutate(effect=normalize_effect(effect)) %>%
    filter(.fail | effect %in% c("te","pnde","tnde","pnie","tnie","pm"))
}

run_cde <- function(nboot=80, m_levels=MMD_LEVELS) {
  grid <- tidyr::expand_grid(dsd=DSD_LEVELS, m=m_levels)
  purrr::pmap_dfr(list(grid$dsd, grid$m), function(lev, mlev) {
    dsub <- make_subset(lev)
    tab <- with(dsub, table(A,M))
    if (any(tab==0) || !(mlev %in% levels(dsub$M)))
      return(tibble(dsd_contrast=lev, mmd_cde_at=mlev, kind="cde",
                    effect=NA, est=NA, se=NA, lo=NA, hi=NA, p=NA, .fail=TRUE,
                    error_msg="Empty A×M cell or mediator level missing"))
    pr <- prune_basec(dsub, basec_vars, event="iit_event", min_per_cell=3L)
    safe_fit(cmest_once(pr$dat, pr$basec, mval_scalar=mlev, EM=FALSE, nboot=nboot),
             dsd_label=lev, extra_cols=list(kind="cde", mmd_cde_at=mlev))
  }) %>% mutate(effect=normalize_effect(effect)) %>% filter(.fail | effect=="cde")
}

natural_tbl <- run_natural(nboot=NBOOT_MAIN) %>% arrange(dsd_contrast)
cde_tbl     <- run_cde(nboot=NBOOT_CDE, m_levels=MMD_LEVELS) %>% arrange(dsd_contrast, mmd_cde_at)

readr::write_csv(natural_tbl, file.path(OUT_DIR, "natural_tbl_timeinvariant.csv"))
readr::write_csv(cde_tbl,     file.path(OUT_DIR, "cde_tbl_timeinvariant.csv"))