# data-raw/generate_dummy_data.R
# Synthetic data generator aligned to time-invariant/varying & mediation scripts

set.seed(123)

N <- 6000L
sites <- sprintf("S%03d", sample(1:150, N, replace=TRUE))

# Factor levels (union across scripts)
AgeARTStartCategory <- sample(c("0-14","15-24","25-34","35-44","45-59","60+"),
                              N, replace=TRUE, prob=c(0.05,0.15,0.45,0.25,0.07,0.03))
Gender              <- sample(c("Female","Male"), N, replace=TRUE, prob=c(0.55,0.45))
EducationLevel      <- sample(c("Tertiary","Secondary","Primary","None","Other"),
                              N, replace=TRUE, prob=c(0.12,0.45,0.32,0.05,0.06))
MaritalStatus       <- sample(c("Other","Single","Married","Divorced/Widowed"),
                              N, replace=TRUE, prob=c(0.05,0.45,0.45,0.05))
Occupation          <- sample(c("Employed","Unemployed","Student","Other"),
                              N, replace=TRUE, prob=c(0.35,0.45,0.10,0.10))
NCD                 <- sample(c("0","1"), N, replace=TRUE, prob=c(0.8,0.2))
Adherence           <- sample(c("Good","Fair","Poor"), N, replace=TRUE, prob=c(0.65,0.25,0.10))
WHOStage            <- sample(c("1","2","3","4"), N, replace=TRUE, prob=c(0.35,0.35,0.20,0.10))
ViralLoadCategory   <- sample(c("<200","200-999",">1000","Unknown"),
                              N, replace=TRUE, prob=c(0.70,0.15,0.10,0.05))
StabilityAssessment <- sample(c("Stable","Unstable"), N, replace=TRUE, prob=c(0.7,0.3))
BaselineRegimen     <- sample(c("Other","TDF/3TC/DTG","ABC/3TC/DTG"),
                              N, replace=TRUE, prob=c(0.1,0.7,0.2))
RegimenLine         <- sample(c("Other","1L","2L","3L+"),
                              N, replace=TRUE, prob=c(0.05,0.85,0.09,0.01))
MMD                 <- sample(c("<3 Months","3-5 Months","6+ Months"),
                              N, replace=TRUE, prob=c(0.25,0.35,0.40))
DifferentiatedCare  <- sample(c("Standard Care",
                                "Community ART Distribution HCW Led",
                                "Community ART Distribution Peer Led",
                                "Facility ART Distribution Group",
                                "Fast Track"),
                              N, replace=TRUE, prob=c(0.32,0.20,0.15,0.08,0.25))
PwP_countCategory   <- sample(c("0","1","2","3 or more"),
                              N, replace=TRUE, prob=c(0.70,0.20,0.07,0.03))
YearStartART        <- sample(as.character(2017:2023), N, replace=TRUE)
CovidLockDown       <- sample(c("0","1"), N, replace=TRUE, prob=c(0.8,0.2))
MissedAppointment   <- sample(c("0","1"), N, replace=TRUE, prob=c(0.80,0.20))
DefaultAppointment  <- sample(c("0","1"), N, replace=TRUE, prob=c(0.88,0.12))

# Time-to-event in months (cap at ~89 months), then days
rate_base <- 1/18
rate_adj  <- ifelse(MMD=="6+ Months", 0.7, ifelse(MMD=="3-5 Months", 0.85, 1.1)) *
             ifelse(DifferentiatedCare=="Standard Care", 1.1, 0.9)
time_m    <- rexp(N, rate = rate_base * rate_adj)
time_m    <- pmin(time_m, 89)
tstop     <- round(time_m * 30.44, 1)

# status: 0=censor, 1=interruption, 2=competing
status <- sample(c(0,1,2), N, replace=TRUE, prob=c(0.25, 0.60, 0.15))

# Time-invariant data + ART_durationMths for mediation
Time_invariant_df <- data.frame(
  patientUID = seq_len(N),
  SiteCode = sites,
  status = status,
  tstop = tstop,
  ART_durationMths = round(time_m, 2),
  AgeARTStartCategory, Gender, EducationLevel, MaritalStatus, Occupation,
  NCD, Adherence, WHOStage, ViralLoadCategory, StabilityAssessment,
  BaselineRegimen, RegimenLine, MMD, DifferentiatedCare, PwP_countCategory,
  YearStartART, CovidLockDown, MissedAppointment, DefaultAppointment
)

# ---------- Time-varying person-period ----------
K <- sample(3:6, N, replace=TRUE)
tv_list <- vector("list", N)
for (i in seq_len(N)) {
  k <- K[i]
  tlen <- sample(c(14, 30, 60, 90), k, replace=TRUE, prob=c(0.2,0.4,0.3,0.1))
  tstart <- c(0, head(cumsum(tlen), -1))
  tstop  <- cumsum(tlen)

  state_from <- integer(k); state_to <- integer(k); cur <- 1L
  for (j in seq_len(k)) {
    state_from[j] <- cur
    if (cur == 1) {
      cur <- 2
    } else if (cur == 2) {
      cur <- sample(c(3,4,5), 1, prob=c(0.55,0.30,0.15))
    } else if (cur == 3) {
      cur <- sample(c(4,5), 1, prob=c(0.65,0.35))
    } else if (cur == 4) {
      cur <- sample(c(3,5), 1, prob=c(0.70,0.30))
    } else if (cur == 5) {
      cur <- 5
    }
    state_to[j] <- cur
  }

  status_iv <- rbinom(k, 1, prob=0.8)

  tv_list[[i]] <- data.frame(
    patientUID = i,
    SiteCode = sample(sites, 1),
    status = status_iv,
    tstart = tstart,
    tstop = tstop,
    State = paste0("State", state_to),
    State_from = paste0("State", state_from),
    State_to   = paste0("State", state_to),
    AgeARTStartCategory = sample(AgeARTStartCategory, k, replace=TRUE),
    Gender = sample(Gender, k, replace=TRUE),
    EducationLevel = sample(EducationLevel, k, replace=TRUE),
    MaritalStatus = sample(MaritalStatus, k, replace=TRUE),
    Occupation = sample(Occupation, k, replace=TRUE),
    NCD = sample(NCD, k, replace=TRUE),
    Adherence = sample(Adherence, k, replace=TRUE),
    WHOStage = sample(WHOStage, k, replace=TRUE),
    ViralLoadCategory = sample(ViralLoadCategory, k, replace=TRUE),
    StabilityAssessment = sample(StabilityAssessment, k, replace=TRUE),
    BaselineRegimen = sample(BaselineRegimen, k, replace=TRUE),
    RegimenLine = sample(RegimenLine, k, replace=TRUE),
    MMD = sample(MMD, k, replace=TRUE),
    DifferentiatedCare = sample(DifferentiatedCare, k, replace=TRUE),
    PwP_countCategory = sample(PwP_countCategory, k, replace=TRUE),
    YearStartART = sample(YearStartART, k, replace=TRUE),
    CovidLockDown = sample(CovidLockDown, k, replace=TRUE),
    MissedAppointmentCumCategory = sample(c("0","1-3","4-6",">7"), k, replace=TRUE, prob=c(0.70,0.20,0.07,0.03)),
    DefaultAppointmentCumCategory = sample(c("0","1-3","4-6",">7"), k, replace=TRUE, prob=c(0.85,0.10,0.04,0.01)),
    TimeReturnToCareCategory = sample(c("0 days","1–90 days","91–180 days","181+ days"),  # note en-dash
                                      k, replace=TRUE, prob=c(0.55,0.25,0.15,0.05))
  )
}

Time_dependent_iit_df <- do.call(rbind, tv_list)

if (!dir.exists("data")) dir.create("data", recursive = TRUE)
readr::write_csv(Time_invariant_df, "data/Time_invariant_df_Final_4Aug25.csv")
readr::write_csv(Time_dependent_iit_df, "data/Time_dependent_iit_df_Final_4Aug25.csv")

message("Wrote synthetic: data/Time_invariant_df_Final_4Aug25.csv and data/Time_dependent_iit_df_Final_4Aug25.csv")