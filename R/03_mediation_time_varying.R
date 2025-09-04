# R/03_mediation_time_varying.R
# Time-varying mediation (visit-based pooled model) + DSD×MMD interaction + CMAverse natural effects/CDEs

source("R/_globals.R")
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(splines); library(sandwich)
  library(car); library(emmeans); library(purrr); library(tibble)
  library(CMAverse)
})

# ---- Load data ----
infile <- file.path(DATA_DIR, "Time_dependent_iit_df_Final_4Aug25.csv")
stopifnot(file.exists(infile))
df <- readr::read_csv(infile, show_col_types = FALSE)

# ---- Factor preparation ----
fac_vars <- c("MissedAppointmentCumCategory","DefaultAppointmentCumCategory",
              "MaritalStatus","RegimenLine","DifferentiatedCare","StabilityAssessment",
              "Adherence","BaselineRegimen","AgeARTStartCategory","Gender","EducationLevel",
              "YearStartART","WHOStage","ViralLoadCategory","MMD","CovidLockDown",
              "PwP_countCategory","NCD","TimeReturnToCareCategory")
df <- df %>% mutate(across(all_of(fac_vars), \(x) if (!is.factor(x)) factor(x) else x))

df <- within(df, {
  MaritalStatus       <- relevel(MaritalStatus,       ref = "Other")
  RegimenLine         <- relevel(RegimenLine,         ref = "Other")
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
  MissedAppointmentCumCategory  <- relevel(MissedAppointmentCumCategory,  ref = "0")
  DefaultAppointmentCumCategory <- relevel(DefaultAppointmentCumCategory, ref = "0")
  TimeReturnToCareCategory <- relevel(TimeReturnToCareCategory, ref = "0 days")
})

# Explicit ordering where needed
df$AgeARTStartCategory <- factor(df$AgeARTStartCategory, levels=c("0-14","15-24","25-34","35-44","45-59","60+"), ordered=TRUE)
df$WHOStage            <- factor(df$WHOStage, levels=c("1","2","3","4"), ordered=TRUE)
df$YearStartART        <- factor(df$YearStartART, levels=as.character(2017:2023), ordered=TRUE)
df$PwP_countCategory   <- factor(df$PwP_countCategory, levels=c("0","1","2","3 or more"), ordered=TRUE)
df$ViralLoadCategory   <- factor(df$ViralLoadCategory, levels=c("<200","200-999",">1000","Unknown"), ordered=TRUE)
df$MissedAppointmentCumCategory   <- factor(df$MissedAppointmentCumCategory, levels=c("0","1-3","4-6",">7"), ordered=TRUE)
df$DefaultAppointmentCumCategory  <- factor(df$DefaultAppointmentCumCategory, levels=c("0","1-3","4-6",">7"), ordered=TRUE)
df$StabilityAssessment <- factor(df$StabilityAssessment, levels=c("Stable","Unstable"), ordered=TRUE)
df$NCD                 <- factor(df$NCD, levels=c("0","1"), ordered=TRUE)
df$TimeReturnToCareCategory <- factor(df$TimeReturnToCareCategory, levels=c("0 days","1–90 days","91–180 days","181+ days"), ordered=TRUE)

# ---- Build interval features ----
df <- df %>%
  mutate(tstart=as.numeric(tstart), tstop=as.numeric(tstop), status=as.integer(status),
         int_len = pmax(1e-8, tstop - tstart), log_int_len = log(int_len)) 

bh <- as.data.frame(ns(df$tstart, df = 4)); colnames(bh) <- paste0("tbase", 1:4)
df <- bind_cols(df, bh)

# Lag mediator
df <- df %>% arrange(patientUID, tstart, tstop) %>% group_by(patientUID) %>%
  mutate(M_lag = dplyr::lag(MMD, default = dplyr::first(MMD))) %>% ungroup()

basec_tv <- c("tbase1","tbase2","tbase3","tbase4","log_int_len",
              "AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation",
              "NCD","Adherence","WHOStage","ViralLoadCategory","StabilityAssessment",
              "BaselineRegimen","RegimenLine","PwP_countCategory","YearStartART","CovidLockDown",
              "MissedAppointmentCumCategory","DefaultAppointmentCumCategory","TimeReturnToCareCategory")
fac_only <- setdiff(basec_tv, c("tbase1","tbase2","tbase3","tbase4","log_int_len"))
df[fac_only] <- lapply(df[fac_only], function(x) if (!is.factor(x)) factor(x) else x)

DSD_LEVELS <- c("Community ART Distribution HCW Led","Community ART Distribution Peer Led",
                "Facility ART Distribution Group","Fast Track")
MMD_LEVELS <- c("<3 Months","3-5 Months","6+ Months")

# ---- A×M interaction tables ----
interaction_table_one <- function(dsd_level, data = df, covars = basec_tv,
                                  mmd_lvls = MMD_LEVELS, link = c("cloglog","logit"),
                                  rg_limit = 5000) {
  link <- match.arg(link)
  dsub <- data %>% filter(DifferentiatedCare %in% c("Standard Care", dsd_level)) %>%
    mutate(Y = as.integer(status),
           A = factor(if_else(DifferentiatedCare==dsd_level, "DSD","SC"), levels=c("SC","DSD")),
           M = droplevels(factor(M_lag, levels = mmd_lvls, ordered = FALSE))) %>%
    tidyr::drop_na(Y, A, M) %>% droplevels()

  # Drop constant covariates
  has_2lvls <- function(x) {
    if (is.factor(x) || is.character(x)) length(unique(x[!is.na(x)])) >= 2
    else stats::sd(as.numeric(x), na.rm=TRUE) > 0
  }
  keep_cov <- covars[covars %in% names(dsub)]
  keep_cov <- keep_cov[sapply(dsub[keep_cov], has_2lvls)]

  rhs <- if (length(keep_cov)) paste(c("A * M", keep_cov), collapse=" + ") else "A * M"
  fml <- stats::as.formula(paste("Y ~", rhs))
  fam <- binomial(link = "cloglog")
  fit <- glm(fml, data = dsub, family = fam)
  V <- sandwich::vcovCL(fit, cluster = dsub$patientUID)

  nm_ok <- colnames(V)
  V_ok <- V; cf_ok <- coef(fit)[nm_ok]
  inter_terms_ok <- grep("^A.*:M", nm_ok, value = TRUE)

  if (length(inter_terms_ok) == 0) {
    omnibus <- tibble::tibble(dsd_contrast = dsd_level, test = "A × M (omnibus)",
                              chisq = NA_real_, df = NA_integer_, p_interaction = NA_real_)
  } else {
    L <- diag(length(nm_ok))[match(inter_terms_ok, nm_ok), , drop = FALSE]
    wtest <- car::linearHypothesis(model = fit, hypothesis.matrix = L, rhs = rep(0, nrow(L)),
                                   coef. = cf_ok, vcov. = V_ok, test = "Chisq")
    omnibus <- tibble::tibble(dsd_contrast = dsd_level, test = "A × M (omnibus)",
                              chisq = unname(wtest$Chisq[2]), df = unname(wtest$Df[2]),
                              p_interaction = signif(unname(wtest$`Pr(>Chisq)`[2]), 3))
  }

  fac_covs <- keep_cov[sapply(dsub[keep_cov], is.factor)]
  nuisance_vars <- setdiff(fac_covs, c("A","M"))
  rg <- try(emmeans::ref_grid(fit, vcov.=V_ok, data=dsub, type="link",
                              cov.reduce=mean, nuisance=nuisance_vars, rg.limit=rg_limit), silent=TRUE)
  scale_label <- "HR"
  if (inherits(rg, "try-error")) {
    per_level <- tibble::tibble(dsd_contrast = dsd_level, M_level = levels(dsub$M),
                                log_ratio = NA_real_, SE=NA_real_, ratio=NA_real_,
                                lo=NA_real_, hi=NA_real_, p=NA_real_, scale=scale_label)
    return(list(omnibus=omnibus, per_level=per_level))
  }

  emm <- emmeans::emmeans(rg, ~ A | M)
  comp <- emmeans::contrast(emm, method="revpairwise", by="M")
  tab <- as.data.frame(comp)
  per_level <- tab %>% transmute(dsd_contrast = dsd_level, M_level = as.character(M),
                                 log_ratio = estimate, SE = SE, ratio = exp(estimate),
                                 lo = exp(estimate - 1.96*SE), hi = exp(estimate + 1.96*SE),
                                 p = p.value, scale = scale_label)
  list(omnibus=omnibus, per_level=per_level)
}

res_list <- lapply(DSD_LEVELS, interaction_table_one)
omnibus_tbl <- bind_rows(lapply(res_list, `[[`, "omnibus"))
perM_tbl    <- bind_rows(lapply(res_list, `[[`, "per_level"))

readr::write_csv(omnibus_tbl, file.path(OUT_DIR, "omnibus_interaction_by_DSD.csv"))
readr::write_csv(perM_tbl,  file.path(OUT_DIR, "conditional_effects_by_MMD_and_DSD.csv"))

# ---- CMAverse time-varying natural effects & CDEs ----
# Build tv-ready data with lagged M and time basis (already above)
tv_df <- df

MMD_LEVELS <- c("<3 Months","3-5 Months","6+ Months")
DSD_LEVELS <- c("Community ART Distribution HCW Led","Community ART Distribution Peer Led",
                "Facility ART Distribution Group","Fast Track")
basec_tv <- c("AgeARTStartCategory","Gender","EducationLevel","MaritalStatus","Occupation",
              "NCD","Adherence","WHOStage","ViralLoadCategory","StabilityAssessment",
              "BaselineRegimen","RegimenLine","PwP_countCategory","YearStartART","CovidLockDown",
              "MissedAppointmentCumCategory","DefaultAppointmentCumCategory","TimeReturnToCareCategory",
              "tbase1","tbase2","tbase3","tbase4","log_int_len")

# prune helper
prune_covars <- function(dsub, covars, event="status", min_per_cell=3L) {
  keep <- covars[covars %in% names(dsub)]
  if (!length(keep)) return(list(dat=dsub, covs=NULL))
  ok <- vapply(keep, function(v) {
    x <- dsub[[v]]
    if (!is.factor(x)) return(stats::sd(as.numeric(x), na.rm=TRUE) > 0)
    x <- droplevels(x); if (nlevels(x) < 2) return(FALSE)
    tab <- table(x, dsub[[event]]); all(rowSums(tab)>=min_per_cell) && all(colSums(tab)>=min_per_cell)
  }, logical(1))
  keep <- keep[isTRUE(ok)]
  list(dat=dsub, covs=if (length(keep)) keep else NULL)
}

cmest_tv <- function(d, basec, nboot=60, EM=FALSE, mlevel=NULL) {
  CMAverse::cmest(data=d, model="rb", estimation="imputation", inference="bootstrap", nboot=nboot,
                  outcome="status", yreg="logistic", mediator="M", mreg=list("ordinal"),
                  exposure="A", basec=basec, EMint=EM, astar="SC", a="DSD",
                  mval=list(if (is.null(mlevel)) "<3 Months" else mlevel))
}

safe_fit <- function(expr, extra_cols=list()) {
  out <- tryCatch(expr, error=function(e) e)
  if (inherits(out,"error"))
    return(tibble(effect=NA_character_, est=NA_real_, se=NA_real_, lo=NA_real_, hi=NA_real_,
                  p=NA_real_, .fail=TRUE, error_msg=conditionMessage(out), !!!extra_cols))
  tibble(effect=names(out$effect.pe), est=as.numeric(out$effect.pe), se=as.numeric(out$effect.se),
         lo=as.numeric(out$effect.ci.low), hi=as.numeric(out$effect.ci.high),
         p=as.numeric(out$effect.pval), .fail=FALSE, error_msg=NA_character_, !!!extra_cols)
}

normalize_effect <- function(eff) { sub("^r","", tolower(eff)) }

tv_nat_by_dsd <- function(tv_df, nboot=60, EM=FALSE, mref="<3 Months") {
  purrr::map_dfr(DSD_LEVELS, function(lev) {
    dsub <- tv_df %>%
      filter(DifferentiatedCare %in% c("Standard Care", lev)) %>%
      transmute(patientUID,
                A = factor(if_else(DifferentiatedCare==lev, "DSD","SC"), levels=c("SC","DSD")),
                M = factor(M_lag, levels=MMD_LEVELS, ordered=TRUE),
                status = as.integer(status),
                across(all_of(basec_tv))) %>%
      drop_na(A,M,status)
    if (any(table(dsub$A,dsub$M)==0))
      return(tibble(dsd_contrast=lev, effect=NA, est=NA, se=NA, lo=NA, hi=NA, p=NA,
                    .fail=TRUE, error_msg="Empty A×M cell", kind="natural"))
    pr <- prune_covars(dsub, basec_tv)
    safe_fit(cmest_tv(pr$dat, pr$covs, nboot=nboot, EM=EM, mlevel=mref),
             extra_cols=list(dsd_contrast=lev, kind="natural"))
  }) %>% mutate(effect=normalize_effect(effect)) %>%
    filter(.fail | effect %in% c("te","pnde","tnde","pnie","tnie","pm"))
}

tv_cde_by_dsd <- function(tv_df, nboot=40, EM=FALSE, m_levels=MMD_LEVELS) {
  grid <- tidyr::expand_grid(dsd=DSD_LEVELS, m=m_levels)
  purrr::pmap_dfr(list(grid$dsd, grid$m), function(lev, mlev) {
    dsub <- tv_df %>%
      filter(DifferentiatedCare %in% c("Standard Care", lev)) %>%
      transmute(patientUID,
                A = factor(if_else(DifferentiatedCare==lev, "DSD","SC"), levels=c("SC","DSD")),
                M = factor(M_lag, levels=MMD_LEVELS, ordered=TRUE),
                status = as.integer(status),
                across(all_of(basec_tv))) %>%
      drop_na(A,M,status)
    if (any(table(dsub$A,dsub$M)==0) || !(mlev %in% levels(dsub$M)))
      return(tibble(dsd_contrast=lev, mmd_cde_at=mlev, effect=NA, est=NA, se=NA, lo=NA, hi=NA, p=NA,
                    .fail=TRUE, error_msg="Empty A×M or M level missing", kind="cde"))
    pr <- prune_covars(dsub, basec_tv)
    safe_fit(cmest_tv(pr$dat, pr$covs, nboot=nboot, EM=EM, mlevel=mlev),
             extra_cols=list(dsd_contrast=lev, mmd_cde_at=mlev, kind="cde"))
  }) %>% mutate(effect=normalize_effect(effect)) %>% filter(.fail | effect=="cde")
}

nat_tbl <- tv_nat_by_dsd(df, nboot = 60, EM = FALSE) %>% arrange(dsd_contrast)
cde_tbl <- tv_cde_by_dsd(df, nboot = 40, EM = FALSE) %>% arrange(dsd_contrast, mmd_cde_at)

readr::write_csv(nat_tbl, file.path(OUT_DIR, "natural_tbl_timevarying.csv"))
readr::write_csv(cde_tbl, file.path(OUT_DIR, "cde_tbl_timevarying.csv"))