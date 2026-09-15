#<<<<<<<<<<<########################################################>>>>>>>>>>>#
#<<<<<<<<<<<########## METEOR Target trial emulation ###############>>>>>>>>>>>#
#<<<<<<<<<<<########################################################>>>>>>>>>>>#

# R version 4.5.1


#==============================================================================#
# This script runs the analyses of the article: 
#
# Effect of biological DMARDs on physical function in patients with 
# axial spondyloarthritis: Results from a strategy trial emulation
#
# Table 1, Table 2 & Supplementary Table S3 were done in Stata V15.1
# Figure 1 & Figure S1 were done with DAGitty and drawn with Excalidraw
#==============================================================================#


#########################################################################
####################### Libraries  ######################################
#########################################################################

# --- Data handling -----------------------------------------------------
library(data.table)
library(dplyr)
# --- Estimators --------------------------------------------------------
library(ltmle)    # Longitudinal targeted minimum loss-based estimation
library(ipw)      # Fits numerator and denominator models, 
                  # returns the cumulative weight within patient.
library(survey)   # Marginal structural model with robust (sandwich) SE
library(gfoRmula) # Parametric g-formula for time-varying treatments.
library(geepack)  # GEE regression
# --- Output ------------------------------------------------------------
library(flextable)
library(ggplot2)
library(scales)
library(writexl)

#########################################################################
####################### Time-varying ltmle ##############################
#######################      Table 3        #############################
#########################################################################

#########################################################################
########>>>>>>>>>>>>>>>>>>> ltmle package (estimand: always vs never)
#########################################################################

packageVersion("ltmle") # 1.3.0

#########################################################################
# Estimator: Longitudinal targeted minimum loss-based estimation (LTMLE)
#
# Estimand: effect on BASFI at 12 months of "always treat with bDMARDs" versus
#           "never treat", sustained over the 6- and 12-month decisions.
#           This static contrast is the same estimand as the g-formula and MSM.
#
# Data:  Wide format, one row per patient (meteor12ltmlemocondswitch.csv),
#        n = 352, no missing values.
#
#        Three time points: baseline, 6 months, 12 months.
#
#-------------------------------------------------------------------------------
# Variables; The last column gives the equivalent in the long DB used in gformula and MSM
#
#
#   variable       source    role                     time point   long-file equivalent
#   -------------  --------  -----------------------  -----------  --------------------
#   id             database  patient index            -            id
#   mny            database  baseline confounder      baseline     mny
#   asasmri        database  baseline confounder      baseline     asasmri
#   hla            database  baseline confounder      baseline     hla
#   sex            database  baseline confounder      baseline     sex
#   age            database  baseline confounder      baseline     age
#   comorbbin      database  baseline confounder      baseline     comorbbin
#   ibdbl          database  baseline confounder      baseline     ibdbl
#   pertvt1        database  baseline confounder      baseline     pertv at t0 = 0
#   emmtvt1        database  baseline confounder      baseline     emmtv at t0 = 0
#   comedtvt1      database  baseline confounder      baseline     comedtv at t0 = 0
#   L0d            database  baseline confounder      baseline     asdastotal at t0 = 0
#   L0e            database  baseline confounder      baseline     basfitotal at t0 = 0
#   A0             database  treatment node           baseline     bionewnext at t0 = 0
#   L1a            database  time-varying confounder  6 months     pertv at t0 = 1
#   L1b            database  time-varying confounder  6 months     emmtv at t0 = 1
#   L1c            database  time-varying confounder  6 months     comedtv at t0 = 1
#   L1d            database  time-varying confounder  6 months     asdastotal at t0 = 1
#   L1e            database  time-varying confounder  6 months     basfitotal at t0 = 1
#   L1s            derived   switching, see below     6 months     switchnext at t0 = 1
#   A1             database  treatment node           6 months     bionewnext at t0 = 1
#   L2a            database  time-varying confounder  12 months    pertv at t0 = 2
#   L2b            database  time-varying confounder  12 months    emmtv at t0 = 2
#   L2c            database  time-varying confounder  12 months    comedtv at t0 = 2
#   L2d            database  time-varying confounder  12 months    asdastotal at t0 = 2
#   Y2             database  outcome                  12 months    basfitotalt3
#   crpelevatedt1  database  subgroup variable        baseline     crpelevatedt1
#
#   L1s is switchnextt2, renamed in Step 0. 
#   id and crpelevatedt1 are not nodes, so ltmle treats them as baseline covariates;
#   neither appears in any Qform or gform, so neither enters a regression.
#   crpelevatedt1 is used only to define the CRP+ and CRP- subgroups.
#
#-------------------------------------------------------------------------------
#   Labels
#
#   mny           modified New York criterion fulfilled
#   asasmri       ASAS-positive MRI of the sacroiliac joints
#   hla           HLA-B27 positive
#   sex           1 = male, 0 = female
#   age           age at baseline, years
#   comorbbin     any comorbidity
#   ibdbl         inflammatory bowel disease at baseline
#   crpelevatedt1 elevated CRP at baseline
#   pertvt1, L1a, L2a   peripheral manifestation: arthritis or enthesitis
#   emmtvt1, L1b, L2b   extra-musculoskeletal manifestation: uveitis or
#                       psoriasis. IBD is excluded
#   comedtvt1, L1c, L2c comedication: NSAID, csDMARD or glucocorticoids
#   L0d, L1d, L2d       ASDAS
#   L0e, L1e            BASFI, 0-10
#   A0, A1              bDMARD decision, "next" convention below
#   L1s                 bDMARD switch decision, "next" convention below
#   Y2                  BASFI at 12 months
#
#-------------------------------------------------------------------------------
# CONVENTIONS
#
#    Column order is the causal order. ltmle defines the parents of each node
#    as every column to its left. 
#
#    In the database switchnextt2 sits after A1. Step 0 renames it to
#    L1s and relocates it before A1, so the 6-month switching decision precedes
#    the 6-month bDMARD decision. 
#
#    The "next" convention. A0 is the decision
#    taken at baseline, covering baseline to 6 months, i.e. bDMARD at 6 months.
#    A1 is the decision taken at 6 months, covering 6 to 12 months, i.e. bDMARD
#    at 12 months. L1s is the switching decision taken at 6 months. There is no
#    third decision node: the decision taken at 12 months would follow the
#    outcome. This is why every rule returns two elements, c(A0, A1), or three
#    when switching is intervened on, c(A0, L1s, A1).
#
#    Note on suffixes. The long file uses a tN suffix meaning visit number,
#    1-indexed, so pertvt1 is baseline. This file uses node names whose digit
#    is the time point, 0-indexed, so L1d is 6 months. The three columns that
#    keep their long-file names here, pertvt1, emmtvt1, comedtvt1 and
#    crpelevatedt1, are baseline under both readings. 
#
#    L1s is treated two ways. In the stop strategies it stays an ordinary
#    covariate, an Lnode. In the switch strategies it becomes an intervention
#    node, an Anode with its own treatment model gSW. The comparator arm is
#    therefore not numerically identical between the two families.
#
#    Fitted models. ltmle groups consecutive L and Y nodes
#    that are not separated by a treatment node into a block, and fits one Q
#    regression at the first node of each block. With L1s as an Lnode there are
#    two blocks, {L1a ... L1s} and {L2a ... Y2}, so exactly two Q regressions
#    are fitted:
#
#      Q at L2a : E[ BASFI at 12M | baseline, A0, L1a-L1e, L1s, A1 ]
#                 the outcome regression
#      Q at L1a : E[ Q(L2a) | baseline, A0 ]
#
#    together with the treatment models, g at A0 and g at A1 (and g at L1s in
#    the switch strategies). Formulas supplied for any other node are dropped
#    with a message at run time; where they are given they only state the
#    assumed structure.
#
#    This is the defining contrast with the parametric g-formula. The g-formula
#    must model every post-treatment confounder because it simulates them
#    forward one visit at a time. LTMLE regresses them out sequentially and
#    never models them, so the 12-month confounders L2a to L2d enter no
#    regression at all.
#===============================================================================

#####>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")

meteor12ltmlemocondswitch <- read.csv(file="meteor12ltmlemocondswitch.csv", header=TRUE, sep=",")
names(meteor12ltmlemocondswitch)[names(meteor12ltmlemocondswitch) == "switchnextt2"] <- "L1s"
meteor12ltmlemocondswitch <- meteor12ltmlemocondswitch %>%
  relocate(L1s, .before = A1)

#####>>>>>> Strategies
never <- function(row) c(0,0)
always <- function(row) c(1,1) 


#####>>>>>> Specify the 2 outcome (L1a and L2a) and the 2 treatment (A0 and A1) models 

Lnodes <- c("L1a", "L1b", "L1c", "L1d", "L1e", "L1s", "L2a", "L2b", "L2c", "L2d")
Anodes <- c("A0","A1")
Ynodes <- c("Y2")
Qform <- c("L1a" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0", ## Only this OUTCOME model is fitted to not adjust away all the L1 post-treatment confounders (this is not fitting L1a), the outcome here is the targeted outcome coming from the L2a model
           "L1b" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0", ## ignored, only here to inform DAG relations
           "L1c" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0", ## ignored, only here to inform DAG relations
           "L1d" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0 + L1a + L1b + L1c", ## ignored, only here to inform DAG relations
           "L1e" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0 + L1a + L1b + L1c + L1d", ## ignored, only here to inform DAG relations
           "L1s" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e      + L1a + L1b + L1c + L1d + L1e", ## ignored, only here to inform DAG relations
           
           "L2a" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1", ## Only this OUTCOME model is fitted to not adjust away all the L2 post-treatment confounders (this is not fitting L2a) and A0 is included here to allow A0 influence Y as in gformula. The outcome here is the observed Y
           "L2b" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e      + L1a + L1b + L1c + L1d + L1e + L1s + A1", ## ignored, only here to inform DAG relations
           "L2c" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e      + L1a + L1b + L1c + L1d + L1e + L1s + A1", ## ignored, only here to inform DAG relations
           "L2d" = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e      + L1a + L1b + L1c + L1d + L1e + L1s + A1 + L2a + L2b + L2c", ## ignored, only here to inform DAG relations
           "Y2"  = "Q.kplus1 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e + A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1 + L2a + L2b + L2c + L2d") ## ignored, only here to inform DAG relations
gform <- c("A0" = "A0 ~ mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e",
           "A1" = "A1 ~ asasmri + hla + age + sex + comorbbin + ibdbl + A0*mny*L1d + L1a + L1b + L1c + L1d + L1e")
EY.11 <- ltmle(meteor12ltmlemocondswitch, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=always, estimate.time = FALSE)
EY.00 <- ltmle(meteor12ltmlemocondswitch, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes, 
               Qform = Qform, gform = gform, 
               rule=never, estimate.time = FALSE)

#####>>>>>> ATE 
ATE <- ltmle(meteor12ltmlemocondswitch, Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes,
             Qform = Qform, gform = gform,
             rule = list(always, never), estimate.time = FALSE)
print(summary(ATE))


#####>>>>>> ATE, PO and e-value (Table 3)

#####>>>>>> Extract estimates
sATE <- summary(ATE)                      
get_est <- function(x) c(est = x$estimate, se = x$std.dev, lo = x$CI[1], hi = x$CI[2])

ate  <- get_est(sATE$effect.measures$ATE)
po11 <- get_est(sATE$effect.measures$treatment)   # always treated
po00 <- get_est(sATE$effect.measures$control)     # never treated


#####>>>>>> e-value for a continuous outcome 
sd_y <- sd(meteor12ltmlemocondswitch$Y2, na.rm = TRUE)
evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * est / sd_y)
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
ev <- evalue_md(ate["est"], sd_y)


#####>>>>>> LTMLE data for Table 3
t3 <- data.frame(
  Method            = "LTMLE",
  Treatment_effect  = sprintf("%.2f (%.2f; %.2f)", ate["est"], ate["lo"], ate["hi"]),
  SE                = sprintf("%.2f", ate["se"]),
  Always_treated_PO = sprintf("%.2f (%.2f)", po11["est"], po11["se"]),
  Never_treated_PO  = sprintf("%.2f (%.2f)", po00["est"], po00["se"]),
  e_value           = sprintf("%.1f", ev),
  stringsAsFactors  = FALSE
)

t3_num <- data.frame(Method = "LTMLE", t(c(ate, po11 = po11, po00 = po00, ev)))
print(t3)

setwd("C:/Users/alexa/OneDrive/work/Projects/Causal_axSpA/METEOR/Data/Main manuscript/0_Final code/Tables/")
writexl::write_xlsx(t3, "Table3ltmle.xlsx")


#####>>>>>> Diagnostics (Supplementary Table S4)

s  <- summary(ATE)$effect.measures$ATE
d  <- meteor12ltmlemocondswitch

# Followers: patients whose observed treatment matches "always treated"
ab   <- t(apply(d, 1, always))
foll <- sum(d$A0 == ab[,1] & d$A1 == ab[,2])

# Truncation: % of patients with cumulative g bounded in either arm
trunc <- 100 * mean(apply(ATE$cum.g != ATE$cum.g.unbounded, 1, any))

# Influence-curve SD
IC <- ATE$IC; if (is.matrix(IC) && ncol(IC) == 2) IC <- IC[,1] - IC[,2]

S4 <- data.frame(
  `Treatment effect (95% CI)` = sprintf("%.2f (%.2f; %.2f)", s$estimate, s$CI[1], s$CI[2]),
  `SE`                        = sprintf("%.2f", s$std.dev),
  `Followers N (%)`           = sprintf("%d (%.0f)", foll, 100 * foll / nrow(d)),
  `Truncation (%)`            = sprintf("%.1f", trunc),
  `Influence curve (SD)`      = sprintf("%.1f", sd(IC)),
  check.names = FALSE)

cat("\nSupplementary Table S4 (N =", nrow(d), ")\n\n")
print(S4, row.names = FALSE, right = FALSE)

setwd("C:/mypath/Tables/")
writexl::write_xlsx(S4, "TableS4ltmle.xlsx")


#########################################################################
#########################   Time-varying MSM  ###########################
#######################         Table 3       ###########################
#########################################################################


#########################################################################
#####>>>>>>>>>>>>>>>> ipw and survey packages (estimand: always vs never)
#########################################################################


packageVersion("ipw") # 1.2.1.1
packageVersion("survey") # 4.5


#===============================================================================
# Estimator: Marginal structural model, inverse probability of treatment weights
#
# Estimand: effect on BASFI at 12 months of "always treat with bDMARDs" versus
#           "never treat", sustained over the 6- and 12-month decisions. Here
#           the contrast is cumA = 2 versus cumA = 0 in a model saturated in
#           cumulative treatment, which is the same estimand the g-formula and
#           LTMLE analyses target.
#
# Data:  meteor12longnext.csv, long format, 3 rows per patient, n = 352.
#        t0 = 0  baseline    t0 = 1  6 months    t0 = 2  12 months
#
#    The database is restricted to 12-month completers and no variable used
#    below has a missing value, apart from the two decisions at t0 = 2.
#
#    Rows are used in two ways. The weight models are fitted on
#    dec = d[d$t0 < 2, ], the two decision visits. The outcome is taken from
#    the t0 == 2 rows only.
#
#-------------------------------------------------------------------------------
# VARIABLES. Every other column in the file is ignored, including
# switchnext, which this estimator does not condition on.
#
#   variable         source    role                     t0=0                 t0=1                 t0=2
#   ---------------  --------  -----------------------  -------------------  -------------------  -------------------
#   id               database  patient index            patient index        patient index        patient index
#   t0               database  time index               time index           time index           time index
#   mny              database  baseline confounder      measured             repeated             not used
#   asasmri          database  baseline confounder      measured             repeated             not used
#   hla              database  baseline confounder      measured             repeated             not used
#   sex              database  baseline confounder      measured             repeated             not used
#   age              database  baseline confounder      measured             repeated             not used
#   comorbbin        database  baseline confounder      measured             repeated             not used
#   ibdbl            database  baseline confounder      measured             repeated             not used
#   pertv            database  time-varying confounder  baseline value       6-month value        not used
#   emmtv            database  time-varying confounder  baseline value       6-month value        not used
#   comedtv          database  time-varying confounder  baseline value       6-month value        not used
#   asdastotal       database  time-varying confounder  baseline value       6-month value        not used
#   basfitotal       database  time-varying confounder  baseline value       6-month value        not used
#   bionewnext       database  treatment                decision BL-6M       decision 6M-12M      NA, row dropped
#   basfitotalt3     database  outcome                  NA                   NA                   BASFI at 12 months
#   bionewnext_lag1  derived   treatment history        0                    baseline decision    not used
#
#   Patient-level, one row per patient in `out` (and in `chk` for the base-R verification):
#
#   cumA    derived  exposure in the MSM: bionewnext summed over t0 = 0 and 1,
#                    so 0, 1 or 2, relabelled Never / Once / Always
#   uw, sw  derived  cumulative unstabilised and stabilised weights, carried
#                    from the t0 == 1 row where ipwtm's running product is
#                    complete
#   uw_t, sw_t  derived  the same weights truncated at the 1st and 99th
#                    percentiles
#
#-------------------------------------------------------------------------------
#   Labels
#
#   mny           modified New York criterion fulfilled
#   asasmri       ASAS-positive MRI of the sacroiliac joints
#   hla           HLA-B27 positive
#   sex           1 = male, 0 = female
#   age           age at baseline, years
#   comorbbin     any comorbidity
#   ibdbl         inflammatory bowel disease at baseline
#   pertv         peripheral manifestation: arthritis or enthesitis
#   emmtv         extra-musculoskeletal manifestation: uveitis or psoriasis.
#                 IBD is excluded
#   comedtv       comedication: NSAID, csDMARD or glucocorticoids
#   asdastotal    ASDAS
#   basfitotal    BASFI, 0-10
#   bionewnext    bDMARD decision, "next" convention below
#   basfitotalt3  BASFI at 12 months, i.e. basfitotal at t0 == 2
#
#-------------------------------------------------------------------------------
# CONVENTIONS
#
#    bionewnext in row t is the decision taken at visit t, covering the
#    interval (t, t+1]. So bionewnext at t0 = 0 is bDMARD at 6 months, at
#    t0 = 1 is bDMARD at 12 months, and at t0 = 2 is NA because it follows the
#    outcome. There are therefore exactly two decisions, and dec drops the
#    t0 == 2 rows.
#
#    Because the decision in row t is taken at visit t, the covariates that
#    precede it are the ones in the same row. The denominator model conditions
#    on same-row covariates for this reason; using lagged ones would leave each
#    decision's own visit unadjusted.
#
#    THE outcome is defined only at t0 = 2 and is merged in from those rows, so
#    sd(out$basfitotalt3) is taken over the same 352 values as
#    sd(basfitotalt3, na.rm = TRUE) in the g-formula script and
#    "sum $Y if t0 == 2" in Stata.
#
#    The denominator model is the same specification as the
#    bionewnext covariate model in the g-formula script, so that any difference
#    between the two estimates is attributable to the estimator rather
#    than to the adjustment set. The one addition here is factor(t0), which
#    lets the two decisions have different intercepts.
#===============================================================================

#####>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")
d <- read.csv("meteor12longnext.csv", header = TRUE, sep = ",") # Open dataset
d <- d[order(d$id, d$t0), ] # Reorder ID and Time (t0)

#####>>>>>> Step 0 - lagged treatment

## ---------------------------------------------------------------------------
## Lags reaching before baseline are set to 0, matching gfoRmula's
## baselags = FALSE, so all three methods use the same history coding.
## ---------------------------------------------------------------------------

d$bionewnext_lag1 <- ave(d$bionewnext, d$id, FUN = function(x) c(0, head(x, -1)))
d$bionewnext_lag1[d$t0 < 1] <- 0
dec <- d[d$t0 < 2, ]           # the two treatment decisions (excludes last visit)

#####>>>>>> Step 1 and 2 - weights

## ---------------------------------------------------------------------------
## Denominator: pooled logistic regression over both decisions, conditioning on
##   the baseline confounders and on the post-treatment confounders measured at
##   the visit at which each decision is taken, plus the three-way interaction
##   between previous bDMARD use, ASDAS and mNY. factor(t0) lets the intercept
##   differ between the two decisions.
## Numerator: treatment history and visit only. Baseline covariates are NOT
##   included, because the marginal structural model has no baseline terms;
##   including them would change the estimand.
##
## ipwtm returns the cumulative product of numerator/denominator within patient,
## which is exactly the weight the marginal structural model needs.
## ---------------------------------------------------------------------------

w <- ipwtm(
  exposure    = bionewnext,
  family      = "binomial",
  link        = "logit",
  numerator   = ~ bionewnext_lag1 + factor(t0),
  denominator = ~ bionewnext_lag1 * asdastotal * mny + basfitotal +
    pertv + emmtv + comedtv +
    asasmri + hla + sex + age + comorbbin + ibdbl + factor(t0),
  id          = id,
  timevar     = t0,
  type        = "all",
  data        = dec)

dec$sw <- w$ipw.weights

## Unstabilised weights, for comparison: same denominator, numerator = 1
wu <- ipwtm(exposure = bionewnext, family = "binomial", link = "logit",
            numerator = ~ 1,
            denominator = ~ bionewnext_lag1 * asdastotal * mny + basfitotal +
              pertv + emmtv + comedtv +
              asasmri + hla + sex + age + comorbbin + ibdbl +
              factor(t0),
            id = id, timevar = t0, type = "all", data = dec)
dec$uw <- wu$ipw.weights

#####>>>>>> Step 3 Prepare the MSM

## ---------------------------------------------------------------------------
## carry the final cumulative weight and cumulative treatment to the
## 12-month row, where the outcome is defined
## ---------------------------------------------------------------------------

last <- dec[dec$t0 == 1, c("id", "sw", "uw")]
cum  <- aggregate(bionewnext ~ id, data = dec, FUN = sum)
names(cum)[2] <- "cumA"

out <- merge(d[d$t0 == 2, c("id", "basfitotalt3")], last, by = "id")
out <- merge(out, cum, by = "id")
out$cumA <- factor(out$cumA, levels = c(0, 1, 2),
                   labels = c("Never", "Once", "Always"))

## Truncation at the 1st and 99th percentiles.

trim <- function(x, ty = 2) {
  q <- quantile(x, c(0.01, 0.99), type = ty)
  pmin(pmax(x, q[1]), q[2])
}
out$sw_t <- trim(out$sw)
out$uw_t <- trim(out$uw)

cat(sprintf("cumA: never %d, once %d, always %d\n",
            sum(out$cumA == "Never"), sum(out$cumA == "Once"),
            sum(out$cumA == "Always")))
cat(sprintf("stabilised cumulative weights: mean %.3f  min %.3f  max %.3f\n",
            mean(out$sw), min(out$sw), max(out$sw)))

#####>>>>>> Step 4  MSM

## ---------------------------------------------------------------------------
## Saturated in cumulative treatment, so the coefficient on "Always" is the
## contrast of always versus never treating. svyglm gives the robust standard
## error, treating the weights as known.
## ---------------------------------------------------------------------------

fit_msm <- function(wt) {
  des <- svydesign(ids = ~ id, weights = as.formula(paste("~", wt)), data = out)
  m   <- svyglm(basfitotalt3 ~ cumA, design = des)
  b   <- coef(m)["cumAAlways"]
  se  <- sqrt(diag(vcov(m)))["cumAAlways"]
  po0 <- coef(m)["(Intercept)"]
  c(effect = b, se = se, lo = b - 1.96*se, hi = b + 1.96*se,
    po_never = po0, po_always = po0 + b)
}

res <- t(sapply(c("uw", "sw", "uw_t", "sw_t"), fit_msm))
print(round(res, 4))


#####>>>>>> ATE, PO and e-value (Table 3)

wt_used <- "sw_t" # stabilised, trimmed at 1st/99th pct
des <- svydesign(ids = ~ id, weights = as.formula(paste("~", wt_used)), data = out)
m   <- svyglm(basfitotalt3 ~ cumA, design = des)
b  <- coef(m); V <- vcov(m)
i0 <- "(Intercept)"; iA <- "cumAAlways"
ci <- confint(m, iA) 
ate_msm  <- c(est = unname(b[iA]), se = unname(sqrt(V[iA, iA])), lo = ci[1], hi = ci[2])
po00_msm <- c(est = unname(b[i0]), se = unname(sqrt(V[i0, i0])))
po11_msm <- c(est = unname(b[i0] + b[iA]),
              se  = unname(sqrt(V[i0, i0] + V[iA, iA] + 2 * V[i0, iA])))

## E-value: same sd_y as the LTMLE and g-formula scripts (the 352 t0 == 2 values)
evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * est / sd_y)
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
sd_y   <- sd(out$basfitotalt3)
ev_msm <- evalue_md(ate_msm["est"], sd_y)

t3_msm <- data.frame(
  Method            = "MSM",
  Treatment_effect  = sprintf("%.2f (%.2f; %.2f)", ate_msm["est"], ate_msm["lo"], ate_msm["hi"]),
  SE                = sprintf("%.2f", ate_msm["se"]),
  Always_treated_PO = sprintf("%.2f (%.2f)", po11_msm["est"], po11_msm["se"]),
  Never_treated_PO  = sprintf("%.2f (%.2f)", po00_msm["est"], po00_msm["se"]),
  e_value           = sprintf("%.1f", ev_msm),
  stringsAsFactors  = FALSE
)
print(t3_msm)

setwd("C:/mypath/Tables/")
writexl::write_xlsx(t3_msm, "Table3msm.xlsx")


#####>>>>>> Diagnostics (Supplementary Table S5)
n_pt      <- nrow(out)
followers <- sum(out$cumA == "Always")    

## Drop (%): patients whose covariate profile perfectly predicts treatment.
## glm keeps such rows rather than deleting them (unlike Stata's logit), so
## detect them as fitted probabilities numerically at 0 or 1, and separately
## count any decision row the model did drop.

stopifnot(length(fitted(w$den.mod)) == nrow(dec))   # no list wise deletion
dec$ps      <- fitted(w$den.mod)
drop_ids    <- unique(dec$id[dec$ps < 1e-8 | dec$ps > 1 - 1e-8])
drop_pct    <- 100 * length(drop_ids) / n_pt

## IPTW: cumulative stabilised weights after trimming, and how many were trimmed
wv     <- out[[wt_used]]
n_trim <- sum(out$sw != out$sw_t)

s5 <- data.frame(
  Treatment_effect = t2_msm$Treatment_effect,
  SE               = t2_msm$SE,
  Followers        = sprintf("%d (%.0f)", followers, 100 * followers / n_pt),
  Drop             = sprintf("%.1f", drop_pct),
  IPTW             = sprintf("%.2f (%.2f; %.2f)", mean(wv), min(wv), max(wv)),
  N_trimmed        = n_trim,
  stringsAsFactors = FALSE
)
print(s5)
writexl::write_xlsx(s5, "TableS5msm.xlsx")


###############################################################################
##  Balancing diagnostics (supplementary Table S7 and S8)
###############################################################################

# Reference: Jackson JW, Diagnostics for confounding of
#   time-varying and other joint exposures, Epidemiology 2016;27:859-869

#===============================================================================
#   Covariate balance between treated and untreated at each bDMARD decision,
#   before and after weighting, computed within strata of treatment history.
#
#   Each table reports Jackson's diagnostic 1 (unweighted columns: time-varying
#   confounding in the study population) and diagnostic 3 (weighted columns:
#   residual time-varying confounding in the weighted population) side by side.
#
#   S7  first decision, taken at baseline (bDMARD at 6 months). Everyone shares
#       the same history at this point, nobody has been treated yet, so there is
#       a single stratum and the table is the familiar treated-versus-untreated
#       comparison.
#
#   S8  second decision, taken at 6 months (bDMARD at 12 months). Patients now
#       arrive with one of two histories, treated or untreated over the first
#       interval, so the table has one panel per history.
#
#-------------------------------------------------------------------------------
# Why check balance of the second decision stratified on history
#
#   The post-weighting balance check, Diagnostic 3, compares
#   prior covariate means across exposure groups among persons who have followed
#   a particular exposure trajectory up to that point in time. The stratification
#   on history is part of the definition of the diagnostic.
#
#   Two reasons to stratify on treatment history:
#
#   1. Stabilized weights only balance conditionally. The numerator of the
#      stabilised weight conditions on treatment history, so the weighted
#      pseudo-population has the 12-month decision independent of the 6-month
#      covariates given the 6-month decision, not marginally. A pooled
#      treated-versus-untreated comparison at the second decision therefore has
#      no interpretation as residual confounding, whatever value it takes.
#
#   2. The covariates measured at 6 months are affected by the first decision.
#      Pooling the two histories makes the second-decision groups differ in
#      their composition of first-decision treatment, so the pooled comparison
#      partly reflects the effect of the drug on ASDAS and BASFI between
#      baseline and 6 months. Within a history stratum both groups received the
#      same first-interval treatment, so that path is closed and a difference
#      is interpretable as confounding.
#
#   The covariates included at each decision are the full covariate history
#   preceding that decision, not only the baseline ones. Restricting the table
#   to baseline covariates would drop the time-varying confounders that motivate
#   the marginal structural model in the first place.
#
#-------------------------------------------------------------------------------

#####>>>>>> Treatment history, built from bionewnext

## bionewnext in row t is the decision taken at visit t and covering (t, t+1],
## so bionewnext at t0 == 0 is the bDMARD decision for the first interval and at
## t0 == 1 for the second. History at a given decision is therefore the sequence
## of bionewnext values in the STRICTLY EARLIER rows:
##
##   decision at t0 = 0 : history is empty, one stratum, all 352 patients
##   decision at t0 = 1 : history is bionewnext at t0 = 0, two strata
##
## hist_lab is that sequence written out, so the panels are self-labelling and
## the same code would extend unchanged to a third decision.

hist0 <- dec[dec$t0 == 0, c("id", "bionewnext")]
names(hist0) <- c("id", "A_int1")

hist_lab <- function(a) ifelse(a == 1, "bDMARD in interval 1", "no bDMARD in interval 1")

#####>>>>>> Helpers

is_binary <- function(x) all(x %in% c(0, 1))

wmean <- function(x, w) sum(w * x) / sum(w)

## Weighted sample variance; with w = 1 it reduces to sum((x - m)^2) / (n - 1),
## so the unweighted columns are the ordinary sample variances.
wvar <- function(x, w) {
  m  <- wmean(x, w)
  s1 <- sum(w)
  s2 <- sum(w^2)
  (s1 / (s1^2 - s2)) * sum(w * (x - m)^2)
}

## Unweighted pooled SD within the panel: the denominator Jackson specifies,
## computed once per covariate and reused for the weighted columns.
pooled_sd <- function(x, a) {
  if (is_binary(x)) {
    p1 <- mean(x[a == 1]); p0 <- mean(x[a == 0])
    sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  } else {
    sqrt((var(x[a == 1]) + var(x[a == 0])) / 2)
  }
}

smd <- function(x, a, w = NULL, denom = NULL) {
  if (is.null(w)) w <- rep(1, length(x))
  if (is.null(denom)) denom <- pooled_sd(x, a)
  m1 <- wmean(x[a == 1], w[a == 1])
  m0 <- wmean(x[a == 0], w[a == 0])
  c(treated = m1, untreated = m0, smd = (m1 - m0) / denom)
}

balance_table <- function(dat, vars, a = "A", w = "wt", w_untrimmed = "wt_raw") {
  A <- dat[[a]]
  W <- dat[[w]]
  has_raw <- !is.null(dat[[w_untrimmed]])
  rows <- lapply(names(vars), function(v) {
    x  <- dat[[v]]
    sd <- pooled_sd(x, A)           # unweighted, fixed across columns
    un <- smd(x, A, denom = sd)
    wt <- smd(x, A, W, denom = sd)
    row <- data.frame(
      Variable         = unname(vars[v]),
      Treated_orig     = sprintf("%.2f", un["treated"]),
      Untreated_orig   = sprintf("%.2f", un["untreated"]),
      SMD_unweighted   = sprintf("%.2f", un["smd"]),
      Treated_pseudo   = sprintf("%.2f", wt["treated"]),
      Untreated_pseudo = sprintf("%.2f", wt["untreated"]),
      SMD_weighted     = sprintf("%.2f", wt["smd"]),
      stringsAsFactors = FALSE)
    if (has_raw)
      row$SMD_weighted_untrimmed <- sprintf("%.2f", smd(x, A, dat[[w_untrimmed]], denom = sd)["smd"])
    row
  })
  do.call(rbind, rows)
}

## Approximate sampling SE of an SMD, as in the header. Evaluated at d = 0,
## which is the relevant case when asking whether an observed value is
## distinguishable from no imbalance.
smd_se <- function(n1, n0, d = 0) sqrt((n1 + n0)/(n1 * n0) + d^2/(2 * (n1 + n0)))

report_panel <- function(dat, vars, label, min_arm = 10) {
  A <- dat[["A"]]
  n1 <- sum(A == 1); n0 <- sum(A == 0)
  cat(sprintf("\n%s\n", strrep("-", 78)))
  cat(sprintf("%s\n", label))
  cat(sprintf("  original population: treated %d, untreated %d, total %d\n", n1, n0, n1 + n0))
  cat(sprintf("  pseudo-population:   treated %.0f, untreated %.0f, total %.0f (mean weight %.3f)\n",
              sum(dat$wt[A == 1]), sum(dat$wt[A == 0]), sum(dat$wt), mean(dat$wt)))
  
  if (min(n1, n0) < 2) {
    cat("  NOT ESTIMABLE: fewer than two patients in one arm.\n")
    return(invisible(NULL))
  }
  
  tab <- balance_table(dat, vars)
  tab$Imbalanced <- ifelse(abs(as.numeric(tab$SMD_weighted)) >= 0.10, "*", "")
  
  se <- smd_se(n1, n0)
  cat(sprintf("  approximate SE of an SMD here: %.2f (smaller arm n = %d)\n", se, min(n1, n0)))
  if (min(n1, n0) < min_arm) {
    cat("  UNINFORMATIVE: the smaller arm is too small for the SMD to carry evidence.\n")
    cat("  Report this panel descriptively and do not apply the 0.10 threshold to it.\n")
  } else if (se >= 0.10) {
    cat("  CAUTION: the SE reaches the 0.10 convention, so individual rows near the\n")
    cat("  threshold cannot be distinguished from sampling variation.\n")
  }
  
  print(tab, row.names = FALSE)
  invisible(tab)
}

#####>>>>>> Supplementary Table S7: first decision, single history stratum

b6 <- dec[dec$t0 == 0, ]
b6$A      <- b6$bionewnext
b6$wt     <- trim(b6$sw)
b6$wt_raw <- b6$sw
stopifnot(nrow(b6) == length(unique(dec$id)), !anyNA(b6$wt))

vars_s7 <- c(asdastotal = "ASDAS",
             basfitotal = "BASFI",
             pertv      = "Per features",
             emmtv      = "EMM",
             comedtv    = "Comedication",
             age        = "Age (years)",
             sex        = "Male sex",
             comorbbin  = "Comorbidities",
             mny        = "mNY",
             asasmri    = "BME MRI-SIJ",
             hla        = "HLA-B27",
             ibdbl      = "IBD")

cat("\n================ Supplementary Table S7: bDMARD decision at baseline ================\n")
cat("History before this decision is empty, so there is a single stratum.\n")
s7 <- report_panel(b6, vars_s7, "All patients (no prior bDMARD decision)")

#####>>>>>> Supplementary Table S8: second decision, one panel per history

## Baseline values and 6-month values side by side, one row per patient. The
## time-fixed covariates, the treatment and the weight come from the t0 == 1
## row, which is the row the second decision is taken in and the row the
## denominator model conditions on.

bl <- dec[dec$t0 == 0, c("id", "asdastotal", "basfitotal", "pertv", "emmtv", "comedtv")]
names(bl) <- c("id", "asdas_bl", "basfi_bl", "per_bl", "emm_bl", "comed_bl")

m6 <- dec[dec$t0 == 1, c("id", "asdastotal", "basfitotal", "pertv", "emmtv", "comedtv",
                         "age", "sex", "comorbbin", "mny", "asasmri", "hla", "ibdbl",
                         "bionewnext")]
names(m6) <- c("id", "asdas_6m", "basfi_6m", "per_6m", "emm_6m", "comed_6m",
               "age", "sex", "comorbbin", "mny", "asasmri", "hla", "ibdbl", "A")

b12 <- merge(bl, m6, by = "id")
b12 <- merge(b12, out[, c("id", "sw", "sw_t")], by = "id")  # cumulative weights
b12 <- merge(b12, hist0, by = "id")                        # first-interval decision
names(b12)[names(b12) == "sw_t"] <- "wt"
names(b12)[names(b12) == "sw"]   <- "wt_raw"
stopifnot(nrow(b12) == nrow(bl), !anyNA(b12$wt), !anyNA(b12$A_int1))

vars_s8 <- c(asdas_bl  = "ASDAS BL",
             basfi_bl  = "BASFI BL",
             per_bl    = "Per features BL",
             emm_bl    = "EMM BL",
             comed_bl  = "Comedication BL",
             asdas_6m  = "ASDAS 6M",
             basfi_6m  = "BASFI 6M",
             per_6m    = "Per features 6M",
             emm_6m    = "EMM 6M",
             comed_6m  = "Comedication 6M",
             age       = "Age (years)",
             sex       = "Male sex",
             comorbbin = "Comorbidities",
             mny       = "mNY",
             asasmri   = "BME MRI-SIJ",
             hla       = "HLA-B27",
             ibdbl     = "IBD")

cat("\n\n=============== Supplementary Table S8: bDMARD decision at 6 months ===============\n")
cat("One panel per treatment history over the first interval.\n")

s8 <- list()
for (h in c(0, 1)) {
  panel <- b12[b12$A_int1 == h, ]
  s8[[as.character(h)]] <- report_panel(panel, vars_s8, hist_lab(h))
}

#####>>>>>> Export

setwd("C:/mypath/Tables/")
sheets <- list(S7 = s7)
if (!is.null(s8[["0"]])) sheets$S8_no_bDMARD_int1  <- s8[["0"]]
if (!is.null(s8[["1"]])) sheets$S8_bDMARD_int1     <- s8[["1"]]
writexl::write_xlsx(sheets, "TableS7S8byhistory.xlsx")

###############################################################################
##  Manual code for double checking
###############################################################################

den <- glm(bionewnext ~ bionewnext_lag1 * asdastotal * mny + basfitotal +
             pertv + emmtv + comedtv + asasmri + hla + sex + age +
             comorbbin + ibdbl + factor(t0), family = binomial, data = dec)
num <- glm(bionewnext ~ bionewnext_lag1 + factor(t0),
           family = binomial, data = dec)
pd <- predict(den, type = "response")
pn <- predict(num, type = "response")
gd <- pd*dec$bionewnext + (1 - pd)*(1 - dec$bionewnext)
gn <- pn*dec$bionewnext + (1 - pn)*(1 - dec$bionewnext)

chk <- aggregate(cbind(uw = 1/gd, sw = gn/gd) ~ id, data = dec, FUN = prod)
chk <- merge(chk, cum, by = "id")
chk <- merge(d[d$t0 == 2, c("id", "basfitotalt3")], chk, by = "id")
chk$uw_t <- trim(chk$uw); chk$sw_t <- trim(chk$sw)

robust <- function(wt) {
  m <- lm(basfitotalt3 ~ factor(cumA), data = chk, weights = chk[[wt]])
  X <- model.matrix(m); u <- residuals(m); ww <- chk[[wt]]
  br <- solve(t(X*ww) %*% X); mt <- t(X*ww*u) %*% (X*ww*u)
  n <- nrow(X); k <- ncol(X)
  V <- br %*% mt %*% br * n/(n - k) # HC1, as Stata's robust
  c(effect = coef(m)[3], se = sqrt(V[3, 3]),
    po_never = coef(m)[1], po_always = coef(m)[1] + coef(m)[3])
}
cat("\nBase-R verification:\n")
print(round(t(sapply(c("uw", "sw", "uw_t", "sw_t"), robust)), 4))


#########################################################################
####################### Time-varying gformula ###########################
#######################      Table 3          ###########################
#########################################################################


#########################################################################
############>>>>>>>>>>>>>>>> gfoRmula package (estimand: always vs never)
#########################################################################

packageVersion("gfoRmula") # 1.1.1

#===============================================================================
# Estimator: Time-varying parametric G-formula
#
# Estimand: effect on BASFI at 12 months of "always treat with bDMARDs" versus
#           "never treat", sustained over the 6- and 12-month decisions.
#
# Data:  meteor12longnext.csv, long format, 3 rows per patient, n = 352.
#        t0 = 0  baseline    t0 = 1  6 months    t0 = 2  12 months
#
#    The file is restricted to 12-month completers and no variable used
#    below has a missing value, apart from the two decisions at t0 = 2.
#
#-------------------------------------------------------------------------------
#
#   variable      source    role                     t0=0                 t0=1                 t0=2
#   ------------  --------  -----------------------  -------------------  -------------------  -------------------
#   id            database  patient index            patient index        patient index        patient index
#   t0            database  time index               time index           time index           time index
#   mny           database  baseline confounder      measured             repeated             repeated
#   asasmri       database  baseline confounder      measured             repeated             repeated
#   hla           database  baseline confounder      measured             repeated             repeated
#   sex           database  baseline confounder      measured             repeated             repeated
#   age           database  baseline confounder      measured             repeated             repeated
#   comorbbin     database  baseline confounder      measured             repeated             repeated
#   ibdbl         database  baseline confounder      measured             repeated             repeated
#   pertv         database  time-varying confounder  baseline value       6-month value        12-month value
#   emmtv         database  time-varying confounder  baseline value       6-month value        12-month value
#   comedtv       database  time-varying confounder  baseline value       6-month value        12-month value
#   asdastotal    database  time-varying confounder  baseline value       6-month value        12-month value
#   basfitotal    database  time-varying confounder  baseline value       6-month value        12-month value
#   bionewnext    database  treatment                decision BL-6M       decision 6M-12M      NA
#   switchnext    database  treatment-related        decision BL-6M       decision 6M-12M      NA
#   basfitotalt3  database  outcome                  NA                   NA                   BASFI at 12 months
#   lag1_VAR      derived   history                  0                    baseline value       6-month value
#   lag2_VAR      derived   history                  0                    0                    baseline value
#
#   lag1_ and lag2_ exist for each of the seven time-varying confounders in
#   covnames. gformula() builds them internally from histories = c(lagged) and
#   histvars; baselags = FALSE is what zero-fills them before baseline.
#
#-------------------------------------------------------------------------------
#   Labels
#
#   mny           modified New York criterion fulfilled
#   asasmri       ASAS-positive MRI of the sacroiliac joints
#   hla           HLA-B27 positive
#   sex           1 = male, 0 = female
#   age           age at baseline, years
#   comorbbin     any comorbidity
#   ibdbl         inflammatory bowel disease at baseline
#   pertv         peripheral manifestation: arthritis or enthesitis
#   emmtv         extra-musculoskeletal manifestation: uveitis or psoriasis.
#                 IBD is excluded
#   comedtv       comedication: NSAID, csDMARD or glucocorticoids
#   asdastotal    ASDAS
#   basfitotal    BASFI, 0-10
#   bionewnext    bDMARD decision, "next" convention below
#   switchnext    bDMARD switch decision, "next" convention below
#   basfitotalt3  BASFI at 12 months, i.e. basfitotal at t0 == 2
#
#-------------------------------------------------------------------------------
# CONVENTIONS
#
#    bionewnext and switchnext in row t are the decisions taken at visit t,
#    covering the interval (t, t+1]. So bionewnext at t0 = 0
#    is bDMARD at 6 months, at t0 = 1 is bDMARD at 12 months, and at t0 = 2 is
#    NA because it follows the outcome. switchnext is identically 0 at
#    t0 = 0, since a first bDMARD is not a switch, and events occur at t0 = 1.
#    This is why interventions are specified as c(0, 0, NA) and c(1, 1, NA):
#    nothing is imposed at the third time point.
#
#    The outcome column is NA at t0 = 0 and t0 = 1, unlike the covariates,
#    which is why sd(basfitotalt3, na.rm = TRUE) here and "sum $Y if t0 == 2"
#    in Stata are computed over the same 352 values. It is also what
#    outcome_type = 'continuous_eof' expects.
#===============================================================================


###################>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")
meteor12 <- read.csv(file="meteor12longnext.csv", header=TRUE, sep=",")
setDT(meteor12)


###################>>>>>> Run gformula

set.seed(as.numeric(Sys.time()))
random_seed <- sample(1:1000000, 1)  # Pick a random integer
ncores <- parallel::detectCores()
ncores <- ncores-1

outcome_type <- 'continuous_eof'
id <- 'id'
time_name <- 't0'
covnames <- c('pertv', 'emmtv', 'comedtv',  'asdastotal', 'basfitotal','bionewnext', 'switchnext')
outcome_name <- 'basfitotalt3'
covtypes <- c('binary', 'binary', 'binary', 'normal', 'normal','binary', 'binary')
histories <- c(lagged)
histvars <- list(c('bionewnext','switchnext', 'pertv', 'emmtv', 'comedtv', 'asdastotal', 'basfitotal'))
covparams <- list(covmodels = c(pertv        ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext + lag1_switchnext + t0,
                                emmtv        ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext + lag1_switchnext + t0,
                                comedtv      ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext + lag1_switchnext + t0,
                                asdastotal   ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext + lag1_switchnext + t0 + pertv + emmtv + comedtv,
                                basfitotal   ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext + lag1_switchnext + t0 + pertv + emmtv + comedtv + asdastotal,
                                bionewnext   ~       asasmri + hla + sex + age + comorbbin              + ibdbl                                                                                                                                                                + lag1_bionewnext*asdastotal*mny            + pertv + emmtv + comedtv              + basfitotal, 
                                switchnext   ~ mny + asasmri + hla + sex + age + comorbbin              + ibdbl                                                                                   + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal                                             + pertv + emmtv + comedtv + asdastotal + basfitotal))
ymodel <-                       basfitotalt3 ~ mny + asasmri + hla + sex + age + comorbbin + lag2_pertv + ibdbl + lag2_emmtv + lag2_comedtv + lag2_asdastotal + lag2_basfitotal + lag2_bionewnext + lag1_pertv + lag1_emmtv + lag1_comedtv + lag1_asdastotal + lag1_basfitotal + lag1_bionewnext                           + pertv + emmtv + comedtv + asdastotal
intvars <- list('bionewnext', 'bionewnext')
interventions <- list(list(c(static, c(0, 0, NA))),
                      list(c(static, c(1, 1, NA))))
int_descript <- c('Never treat', 'Always treat')
nsimul   <- 1000
nsamples <-1000
gform_cont_eof <- gformula(obs_data = meteor12,
                           outcome_type = outcome_type, id = id,
                           time_name = time_name, covnames = covnames,
                           outcome_name = outcome_name, covtypes = covtypes,
                           covparams = covparams, ymodel = ymodel,
                           intvars = intvars, interventions = interventions,
                           int_descript = int_descript, histories = histories,
                           histvars = histvars, basecovs = c("hla", "mny", "asasmri", "sex", "age", "comorbbin", "arthtvt1", "enthtvt1", "pertvt1" , "psotvt1", "ibdbl", "aautvt1", "emmtvt1" , "nsaidnewt1", "csdmardnewt1", "gcnewt1", "comedtvt1"),
                           nsimul = nsimul, ref_int =1, seed =random_seed, nsamples =nsamples,
                           sim_data_b=FALSE, baselags = FALSE, parallel = TRUE, ncores = ncores
                           ,boot_diag = TRUE)
gform_cont_eof



###################>>>>>> Summary table and diagnostics

#####>>>>>> Manual ATE (bootstrap percentile 95% CI)
all_contrast_boot <- gform_cont_eof$bootests[["Always treat"]] -gform_cont_eof$bootests[["Never treat"]]
all_ATE     <- gform_cont_eof$result[Interv.==2, `g-form mean`] - gform_cont_eof$result[Interv.==1, `g-form mean`]
all_se_ATE  <- sd(all_contrast_boot, na.rm = TRUE)
all_ci      <- quantile(all_contrast_boot, c(.025,.975), na.rm = TRUE)
all_ci_low  <- all_ci[1]
all_ci_high <- all_ci[2]

#####>>>>>> Manual Regime-specific estimates (bootstrap percentile 95% CI)
Obs_Mean        <- mean(meteor12$basfitotalt3, na.rm = TRUE)
all_usual       <- gform_cont_eof$result[Interv.==0, `g-form mean`]
all_common      <- gform_cont_eof$result[Interv.==1, `g-form mean`]
all_always      <- gform_cont_eof$result[Interv.==2, `g-form mean`]
all_usual_boot  <- gform_cont_eof$bootests[["Natural course"]]
all_common_boot <- gform_cont_eof$bootests[["Never treat"]]
all_always_boot <- gform_cont_eof$bootests[["Always treat"]]
all_u_ci <- quantile(all_usual_boot,  c(.025,.975), na.rm=TRUE)
all_c_ci <- quantile(all_common_boot, c(.025,.975), na.rm=TRUE)
all_a_ci <- quantile(all_always_boot, c(.025,.975), na.rm=TRUE)

#####>>>>>> Followers (Always treated)
wide_treat_all <- meteor12[t0 %in% c(0,1), .(
  bionewnext_t0 = bionewnext[t0 == 0],
  bionewnext_t1 = bionewnext[t0 == 1]
), by = id]

all_followers_N <- wide_treat_all[bionewnext_t0 == 1 & bionewnext_t1 == 1, .N]
all_followers_Perc <- (all_followers_N / (nrow(meteor12)/3)) * 100

#####>>>>>>  E-value for the point estimate, on the standardised scale
sd_Y <- sd(meteor12$basfitotalt3, na.rm = TRUE)
RR   <- exp(0.91 * all_ATE / sd_Y)
RR_p <- if (RR < 1) 1 / RR else RR
evalue_point <- RR_p + sqrt(RR_p * (RR_p - 1))


#####>>>>>> Summary table
all_summary_row <- data.table(
  Patients     = "All",
  N            = nrow(meteor12)/3,
  Followers_N = all_followers_N,
  Followers_Perc = sprintf("%.1f", all_followers_Perc),
  nsimul       = nsimul,
  nsamples     = nsamples,
  ATE          = sprintf("%.2f", all_ATE),
  SE_ATE       = sprintf("%.2f", all_se_ATE),
  CI_ATE       = sprintf("(%.2f; %.2f)", all_ci_low, all_ci_high),
  Always       = sprintf("%.2f", all_always),
  SE_Always    = sprintf("%.2f", sd(all_always_boot, na.rm = TRUE)),
  CI_Always    = sprintf("(%.2f; %.2f)", all_a_ci[1], all_a_ci[2]),
  Never        = sprintf("%.2f", all_common),
  SE_Never     = sprintf("%.2f", sd(all_common_boot, na.rm = TRUE)),
  CI_Never     = sprintf("(%.2f; %.2f)", all_c_ci[1], all_c_ci[2]),
  Obs_Mean     = sprintf("%.2f", Obs_Mean),
  Usual        = sprintf("%.2f", all_usual),
  SE_Usual     = sprintf("%.2f", sd(all_usual_boot, na.rm = TRUE)),
  CI_Usual     = sprintf("(%.2f; %.2f)", all_u_ci[1], all_u_ci[2]),
  evalue_point = sprintf("%.2f", evalue_point))

#####>>>>>> Print table
all_ft <- flextable(all_summary_row)
all_ft <- align(all_ft, align = "right", part = "all")
all_ft <- autofit(all_ft)
all_ft



#####>>>>>> ATE, PO and e-value (Table 2) 

# 1. Create the formatted data frame
table3_data <- data.table(
  ATE          = sprintf("%.2f", all_ATE),
  CI_ATE       = sprintf("(%.2f; %.2f)", all_ci_low, all_ci_high),
  SE_ATE       = sprintf("%.2f", all_se_ATE),
  Always       = sprintf("%.2f", all_always),
  SE_Always    = sprintf("%.2f", sd(all_always_boot, na.rm = TRUE)),
  Never        = sprintf("%.2f", all_common),
  SE_Never     = sprintf("%.2f", sd(all_common_boot, na.rm = TRUE)),
  evalue_point = sprintf("%.2f", evalue_point)
)

# 2. Export the data frame directly to Excel
setwd("C:/mypath/Tables/")
writexl::write_xlsx(table3_data, "Table3gformula.xlsx")

# 3. Create the flextable for viewer display/Word output
table3gformula <- flextable(table3_data)
table3gformula


#####>>>>>> Diagnostics (Supplementary Table S6)

n_pt <- uniqueN(meteor12$id) 

s6 <- data.table(
  Treatment_effect = sprintf("%.2f (%.2f; %.2f)", all_ATE, all_ci_low, all_ci_high),
  SE               = sprintf("%.2f", all_se_ATE),
  Followers        = sprintf("%d (%.0f)", all_followers_N,
                             100 * all_followers_N / n_pt),
  Observed_mean    = sprintf("%.2f", Obs_Mean),
  Estimated_mean   = sprintf("%.2f (%.2f)", all_usual,
                             sd(all_usual_boot, na.rm = TRUE)))

s6_ft <- flextable(s6); s6_ft
writexl::write_xlsx(s6, "TableS6gformula.xlsx")


################################################################################
############################ Time-varying GEE ##################################
################################ Table 3 #######################################
################################################################################

## Estimator: Time-varying GEE
##
## Estimand: The longitudinal effect of a bDMARD during the preceding interval
## on BASFI measured at the end of that interval, adjusted for the autoregressor
## and the confounder history in the previous visit. 
## 
## This is not the sustained-regime contrast of the LTMLE, MSM and g-formula. 
## It is the traditional approach included in Table 2 for comparison.
##
## Data: meteor12longnext.csv, long format, 3 rows per patient, n = 352.
##       t0 = 0 baseline, t0 = 1 six months, t0 = 2 twelve months.
##       Only t0 = 1 and t0 = 2 contribute, because the treatment is missing at
##       baseline by construction: 704 observations, two rows per patient.
##
## CONVENTIONS
##
##   The treatment has to be recoded. The long file stores bDMARD under the
##   "next" convention required by gfoRmula: bionewnext in row t0 is the decision
##   taken at visit t0, covering (t0, t0+1]. The GEE needs the opposite, the
##   bDMARD covering the interval that ends at the visit where the outcome is
##   measured. The two are one visit apart:
##
##       bionew(t0) = bionewnext(t0 - 1)
##
##       t0 = 0 (BL)   no previous row  -> NA, row drops out
##       t0 = 1 (6M)   bionewnext at BL -> bDMARD started between BL and 6M
##       t0 = 2 (12M)  bionewnext at 6M -> bDMARD started between 6M and 12M
##
##   Lagged confounders: ASDAS, BASFI, pertv and emmtv are measured at the visit,
##   so they enter at t0 - 1 to precede the treatment interval. comedtv does NOT
##   get a lag: it is built from nsaidnew, csdmardnew and gcnew, which carry the
##   same lagged-by-design convention as bionew, so the same-row value already
##   refers to the treatment interval.
##
##   The outcome is basfitotal at each visit, not basfitotalt3. basfitotalt3 is
##   defined only at t0 = 2 and would leave a single row per patient.
################################################################################

###################>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")
gee_d <- read.csv("meteor12longnext.csv", header = TRUE, sep = ",")
gee_d <- gee_d[order(gee_d$id, gee_d$t0), ]


###################>>>>>> Step 0: treatment and lagged confounders

## One visit back within patient. The first row of each patient becomes NA,
## which is what makes the baseline row drop out of every model below.
lag1 <- function(x, id) ave(x, id, FUN = function(z) c(NA, head(z, -1)))

gee_d$bionew          <- lag1(gee_d$bionewnext, gee_d$id)   # bDMARD covering the interval ending here
gee_d$lag1_asdastotal <- lag1(gee_d$asdastotal, gee_d$id)
gee_d$lag1_basfitotal <- lag1(gee_d$basfitotal, gee_d$id)
gee_d$lag1_pertv      <- lag1(gee_d$pertv,      gee_d$id)
gee_d$lag1_emmtv      <- lag1(gee_d$emmtv,      gee_d$id)

gee_g <- gee_d[gee_d$t0 > 0, ]
stopifnot(nrow(gee_g) == 2 * length(unique(gee_d$id)), !anyNA(gee_g$bionew))
cat(sprintf("GEE sample: %d observations, %d patients\n",
            nrow(gee_g), length(unique(gee_g$id))))


###################>>>>>> Step 1: model terms

A <- "bionew"                                        # treatment
Y <- "basfitotal"                                    # outcome at the visit
P <- "lag1_asdastotal"                               # mediator/confounder
D <- "lag1_basfitotal"                               # autoregressor
L <- c("lag1_pertv", "lag1_emmtv", "comedtv")        # comedtv unlagged
W <- c("age", "sex", "comorbbin", "mny", "asasmri", "hla", "ibdbl")   # BL time-fixed

f_crude <- as.formula(paste(Y, "~", A, "+", D))
f_part  <- as.formula(paste(Y, "~", A, "+", P, "+", D, "+", paste(W, collapse = " + ")))
f_full  <- as.formula(paste(Y, "~", A, "+", P, "+", D, "+",
                            paste(L, collapse = " + "), "+", paste(W, collapse = " + ")))


#####>>>>>> Step 2: autoregressive time-lagged GEE

fit_gee <- function(f)
  geeglm(f, id = id, waves = t0, data = gee_g,
         family = gaussian, corstr = "exchangeable")

m_crude <- fit_gee(f_crude)
m_part  <- fit_gee(f_part)
m_full  <- fit_gee(f_full)

gee_row <- function(m, label) {
  s  <- summary(m)$coefficients
  b  <- s[A, "Estimate"]
  se <- s[A, "Std.err"]
  data.frame(Model = label,
             effect = b, se = se, lo = b - 1.96 * se, hi = b + 1.96 * se,
             stringsAsFactors = FALSE)
}
print(do.call(rbind, list(gee_row(m_crude, "Crude"),
                          gee_row(m_part,  "Without other time-varying"),
                          gee_row(m_full,  "Fully adjusted"))), digits = 3)

cat(sprintf("\nWorking correlation (exchangeable alpha) = %.4f\n",
            summary(m_full)$corr[1, 1]))


#####>>>>>> ATE, marginal outcomes and e-value (Table 2)

s_full  <- summary(m_full)$coefficients
b_gee   <- unname(s_full[A, "Estimate"])
se_gee  <- unname(s_full[A, "Std.err"])
ate_gee <- c(est = b_gee, se = se_gee,
             lo = b_gee - 1.96 * se_gee, hi = b_gee + 1.96 * se_gee)

mf   <- model.frame(m_full)
nd0  <- nd1 <- gee_g[rownames(mf), ]
nd0[[A]] <- 0
nd1[[A]] <- 1
po0_gee <- mean(predict(m_full, newdata = nd0))
po1_gee <- mean(predict(m_full, newdata = nd1))

## E-value on the same scale as the other methods, using the SD of the outcome
## in the estimation sample.
evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * est / sd_y)
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
sd_y_gee <- sd(mf[[Y]])
ev_gee   <- evalue_md(ate_gee["est"], sd_y_gee)

t3_gee <- data.frame(
  Method            = "GEE",
  Treatment_effect  = sprintf("%.2f (%.2f; %.2f)", ate_gee["est"], ate_gee["lo"], ate_gee["hi"]),
  SE                = sprintf("%.2f", ate_gee["se"]),
  Always_treated_PO = "NA",
  Never_treated_PO  = "NA",
  e_value           = sprintf("%.1f", ev_gee),
  stringsAsFactors  = FALSE
)
print(t3_gee)
cat(sprintf("\nPredictive margins: untreated %.2f, treated %.2f (difference %.2f)\n",
            po0_gee, po1_gee, po1_gee - po0_gee))

setwd("C:/mypath/Tables/")
writexl::write_xlsx(t3_gee, "Table3gee.xlsx")



#########################################################################
####################### Time-varying ltmle ##############################
####################  Supplementary Table S9 #######################3###
#######################      Figure 3      ##############################
#########################################################################


###############################################################################
## Time-varying LTMLE — All strategies, stratified on CRP
###############################################################################

###################>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")
d <- read.csv(file = "meteor12ltmlemocondswitch.csv", header = TRUE, sep = ",")
names(d)[names(d) == "switchnextt2"] <- "L1s"
d <- d %>% relocate(L1s, .before = A1)
setDT(d)

###################>>>>>> Models

## ---------------------------------------------------------------------------
## Models (identical across strategies and identical to LTMLE above) 
## ---------------------------------------------------------------------------
B <- "mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e"

Qform <- c("L1a" = paste("Q.kplus1 ~", B, "+ A0"), ## Model Targeted outcome coming from the model below (ignores L1 confounders)
           "L2a" = paste("Q.kplus1 ~", B, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1")) ## A0 added to allow A0 to influence Y ## Model Observed outcome (ignores L2 confounders)

## Only these two are fitted. ltmle groups consecutive L/Y nodes
## formulas for the other nodes would be dropped anyway

gA0 <- paste("A0 ~", B)
gA1 <- "A1 ~ asasmri + hla + age + sex + comorbbin + ibdbl + A0*mny*L1d + L1a + L1b + L1c + L1d + L1e"
gSW <- paste("L1s ~", B, "+ L1a + L1b + L1c + L1d + L1e")

GB <- c(0.01, 0.99) # gbounds

## Only treated patients can switch. Declaring this deterministic stops the
## switching model from assigning positive probability where switching is
## structurally impossible.
det.g <- function(data, current.node, nodes) {
  if (!"L1s" %in% names(data)) return(NULL)
  if (current.node != which(names(data) == "L1s")) return(NULL)
  list(is.deterministic = data$A0 == 0, prob1 = 0)
}

###################>>>>>> Strategies

## ---------------------------------------------------------------------------
## Thresholds: LDA 2.1, inactive disease 1.3
## clinically important improvement 1.1 on the change score.
## ---------------------------------------------------------------------------

STRAT <- list(
  list(lab = "Static (always treated)",    sw = FALSE,
       f = function(r) c(1, 1)),
  list(lab = "Dynamic 1 (stop if no LDA)", sw = FALSE,
       f = function(r) c(1, as.numeric(r[["L1d"]] <  2.1))),
  list(lab = "Dynamic 2 (stop if no ID)",  sw = FALSE,
       f = function(r) c(1, as.numeric(r[["L1d"]] <  1.3))),
  list(lab = "Dynamic 3 (switch if no LDA)", sw = TRUE,
       f = function(r) c(1, as.numeric(r[["L1d"]] >= 2.1), 1)),
  list(lab = "Dynamic 4 (switch if no ID)",  sw = TRUE,
       f = function(r) c(1, as.numeric(r[["L1d"]] >= 1.3), 1)),
  list(lab = "Dynamic 5 (stop if no CII)", sw = FALSE,
       f = function(r) c(1, as.numeric(r[["L1d"]] - r[["L0d"]] <  -1.1))),
  list(lab = "Dynamic 6 (switch if no CII)", sw = TRUE,
       f = function(r) c(1, as.numeric(r[["L1d"]] - r[["L0d"]] >= -1.1), 1)))

SUB <- list("All"  = d,
            "CRP+" = d[crpelevatedt1 == 1],
            "CRP-" = d[crpelevatedt1 == 0])

###################>>>>>> Function to run the strategy

## ---------------------------------------------------------------------------
## One strategy in one subgroup
## ---------------------------------------------------------------------------
run_one <- function(dat, s) {
  
  Anodes <- if (s$sw) c("A0", "L1s", "A1") else c("A0", "A1")
  Lnodes <- if (s$sw) c("L1a","L1b","L1c","L1d","L1e","L2a","L2b","L2c","L2d")
  else      c("L1a","L1b","L1c","L1d","L1e","L1s","L2a","L2b","L2c","L2d")
  gform  <- if (s$sw) c("A0" = gA0, "L1s" = gSW, "A1" = gA1)
  else      c("A0" = gA0, "A1" = gA1)
  ref    <- if (s$sw) function(r) c(0, 0, 0) else function(r) c(0, 0)
  
  args <- list(data = dat, Anodes = Anodes, Lnodes = Lnodes, Ynodes = "Y2",
               Qform = Qform, gform = gform, rule = list(s$f, ref),
               gbounds = GB, estimate.time = FALSE)
  if (s$sw) args$deterministic.g.function <- det.g
  fit <- do.call(ltmle, args)
  
  e <- summary(fit)$effect.measures
  
  ## Followers: the rule applied row by row, compared with what was observed.
  ab   <- t(apply(dat, 1, s$f))
  foll <- if (s$sw) mean(dat$A0 == ab[,1] & dat$L1s == ab[,2] & dat$A1 == ab[,3])
  else      mean(dat$A0 == ab[,1] & dat$A1 == ab[,2])
  
  ## Truncation, reported two ways. In Figure 3 the switch rows correspond to
  ## the strategy arm alone and the static row to the union of both arms
  tm <- fit$cum.g != fit$cum.g.unbounded
  tr_strategy <- 100 * mean(apply(tm[,,1, drop = FALSE], 1, any))
  tr_either   <- 100 * mean(apply(tm, 1, any))
  
  IC <- fit$IC
  if (is.matrix(IC) && ncol(IC) == 2) IC <- IC[,1] - IC[,2]
  
  ## E-value on the standardised scale (VanderWeele & Ding)
  sdY <- sd(dat$Y2, na.rm = TRUE)
  ev  <- function(b) { rr <- exp(0.91 * abs(b / sdY)); rr + sqrt(rr * (rr - 1)) }
  ci_near <- if (abs(e$ATE$CI[1]) < abs(e$ATE$CI[2])) e$ATE$CI[1] else e$ATE$CI[2]
  
  data.table(n = nrow(dat),
             ate = e$ATE$estimate, se = e$ATE$std.dev,
             lo = e$ATE$CI[1], hi = e$ATE$CI[2],
             po_strategy = e$treatment$estimate, se_strategy = e$treatment$std.dev,
             po_never    = e$control$estimate,   se_never    = e$control$std.dev,
             followers_pct = 100 * foll,
             trunc_strategy_pct = tr_strategy, trunc_either_pct = tr_either,
             sd_IC = sd(IC),
             evalue = ev(e$ATE$estimate), evalue_ci = ev(ci_near))
}

###################>>>>>> Run

results <- rbindlist(lapply(names(SUB), function(g)
  rbindlist(lapply(STRAT, function(s)
    cbind(subgroup = g, strategy = s$lab,
          suppressMessages(run_one(SUB[[g]], s)))))))

results[, ate_ci := sprintf("%.2f (%.2f; %.2f)", ate, lo, hi)]
results[, po_strategy_lab := sprintf("%.2f (%.2f)", po_strategy, se_strategy)]
results[, po_never_lab    := sprintf("%.2f (%.2f)", po_never,    se_never)]

###################>>>>>> Full table

for (g in names(SUB)) {
  r <- results[subgroup == g]
  cat("\n", strrep("-", 108), "\n", sep = "")
  cat(sprintf("LTMLE. Subgroup: %s (n = %d). Comparator: never treated\n", g, r$n[1]))
  cat(strrep("-", 108), "\n", sep = "")
  cat(sprintf("%-29s %-22s %5s %14s %14s %6s %7s %5s %5s\n",
              "Strategy", "Effect (95% CI)", "SE", "Strategy PO",
              "Never PO", "Foll%", "Trunc%", "sdIC", "E"))
  for (i in seq_len(nrow(r)))
    cat(sprintf("%-29s %-22s %5.2f %14s %14s %6.1f %7.1f %5.1f %5.1f\n",
                r$strategy[i], r$ate_ci[i], r$se[i], r$po_strategy_lab[i],
                r$po_never_lab[i], r$followers_pct[i],
                r$trunc_strategy_pct[i], r$sd_IC[i], r$evalue[i]))
}
cat("\nTrunc%: strategy arm. results$trunc_either_pct holds the union of both arms.\n")

## Full numeric table for export
print(results)
setwd("C:/mypath/Tables/")
writexl::write_xlsx(results, "TableS9 & Figure 3.xlsx")




###################>>>>>> Figure 3


fsize <- 4.0
CLIP  <- c(-2.5, 1.0)    # CI whiskers clipped to the plotted range

# ================== 1. Build the data from `results` ==================

strat_order <- c(                       # order in `results` (STRAT labels)
  "Static (always treated)",
  "Dynamic 1 (stop if no LDA)",
  "Dynamic 2 (stop if no ID)",
  "Dynamic 3 (switch if no LDA)",
  "Dynamic 4 (switch if no ID)",
  "Dynamic 5 (stop if no CII)",
  "Dynamic 6 (switch if no CII)")

strat_label <- c(                       # capitalised labels for the figure
  "Static (Always treated)",
  "Dynamic 1 (Stop if no LDA)",
  "Dynamic 2 (Stop if no ID)",
  "Dynamic 3 (Switch if no LDA)",
  "Dynamic 4 (Switch if no ID)",
  "Dynamic 5 (Stop if no CII)",
  "Dynamic 6 (Switch if no CII)")
names(strat_label) <- strat_order

df <- as.data.frame(results) %>%
  filter(subgroup %in% c("CRP+", "CRP-")) %>%
  mutate(
    strat_base            = strategy,
    ate_est               = ate,
    ci_low_trunc          = pmax(lo, CLIP[1]),
    ci_high_trunc         = pmin(hi, CLIP[2]),
    ate_ci_label          = ate_ci,
    ate_se_label          = sprintf("%.2f", se),
    ey11_label            = po_strategy_lab,
    ey00_label            = po_never_lab,
    percent_followers_all = sprintf("%.1f", followers_pct),
    pct_gtrunc            = sprintf("%.1f", trunc_strategy_pct),
    sd_IC                 = sprintf("%.2f", sd_IC),
    evalue_point          = sprintf("%.1f", evalue))

rows <- list()
for (s in strat_order) {
  rp  <- df %>% filter(strat_base == s, subgroup == "CRP+")
  rm  <- df %>% filter(strat_base == s, subgroup == "CRP-")
  hdr <- rp[1, ] %>% mutate(row_type = "header", row_label = strat_label[[s]])
  rp  <- rp %>% mutate(row_type = "crp", row_label = " CRP+")
  rm  <- rm %>% mutate(row_type = "crp", row_label = " CRP-")
  rows <- c(rows, list(hdr, rp, rm))
}

df2 <- bind_rows(rows)

n_rows         <- nrow(df2)
df2$plot_order <- seq_len(n_rows)
df2$row_id     <- n_rows + 1L - df2$plot_order
df2$row_fac    <- factor(df2$row_id, levels = seq_len(n_rows))

df_bg <- df2 %>%
  group_by(group = (row_id - 1) %/% 3) %>%
  summarise(
    xmin = min(as.numeric(row_fac)) - 0.5,
    xmax = max(as.numeric(row_fac)) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(fill = rep(c("grey97", "grey90"), length.out = n()))

df_data <- df2 %>% filter(row_type == "crp")
df_hdr  <- df2 %>% filter(row_type == "header")

y_cols      <- c(se = 2.50, ey11 = 3.60, ey00 = 4.80, follow = 6.20,
                 pct_gtrunc = 7.40, sd_IC = 8.60, evalue_point = 9.60)
y_label_pos <- -6.7

header_shift <- +0.55

# ================== 2. Forest plot  ==================

forest_plot <- ggplot(df2, aes(x = row_fac, y = ate_est)) +
  
  geom_rect(
    data = df_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 1
  ) +
  scale_fill_identity() +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  
  geom_linerange(
    data = df_data,
    aes(ymin = ci_low_trunc, ymax = ci_high_trunc),
    linewidth = 0.6,
    color = "black"
  ) +
  
  geom_point(
    data = df_data,
    shape = 18,
    size = 3.5,
    color = "black"
  ) +
  
  geom_text(
    data = df_data,
    aes(label = ate_ci_label),
    hjust = 0.5,
    vjust = -0.68,
    size = fsize
  ) +
  
  # Right columns
  geom_text(data = df_data, aes(y = y_cols["se"],           label = ate_se_label),          size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["ey11"],         label = ey11_label),            size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["ey00"],         label = ey00_label),            size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["follow"],       label = percent_followers_all), size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["pct_gtrunc"],   label = pct_gtrunc),            size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["sd_IC"],        label = sd_IC),                 size = fsize, vjust = 0.5) +
  geom_text(data = df_data, aes(y = y_cols["evalue_point"], label = evalue_point),          size = fsize, vjust = 0.5) +
  
  # Right column headers
  annotate(
    "text",
    x = n_rows + 0.8 + header_shift,
    y = y_cols,
    label = c("SE", "Treated\nPO (SE)", "Never treated\nPO (SE)", "Followers\n%",
              "Truncation\n %", "Infl curve\nSD", "e-value"),
    size = fsize,
    fontface = "bold",
    lineheight = 0.82,
    vjust = 0.5
  ) +
  
  # Left: strategy headers and CRP labels
  geom_text(
    data = df_hdr,
    aes(y = y_label_pos, label = row_label),
    hjust = 0,
    size = fsize,
    fontface = "bold"
  ) +
  geom_text(
    data = df_data,
    aes(y = y_label_pos, label = row_label),
    hjust = 0,
    size = fsize
  ) +
  
  annotate(
    "text",
    x = n_rows + 0.8 + header_shift,
    y = y_label_pos,
    label = "Strategy",
    size = fsize,
    fontface = "bold",
    hjust = 0
  ) +
  
  scale_y_continuous(
    limits = c(-7.8, 11.0),
    breaks = c(-2.5, -2, -1, 0, 1),
    oob    = scales::oob_keep
  ) +
  scale_x_discrete(drop = FALSE) +
  
  coord_flip(clip = "off") +
  
  labs(
    title   = NULL,
    y       = NULL,
    x       = NULL,
    caption = "\u2190  Better function               Worse function  \u2192"
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    panel.background = element_rect(fill = "gray95", color = NA),
    plot.background  = element_rect(fill = "gray95", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.ticks       = element_blank(),
    axis.line        = element_blank(),
    plot.caption     = element_text(
      size   = 12,
      face   = "italic",
      colour = "grey25",
      hjust  = 0.40,
      margin = margin(t = 10)
    ),
    plot.margin      = margin(40, 150, 10, 10)
  )

forest_plot

ggsave(filename = "Figure 3.png",
       path = "C:/mypath/Figures/",
       width = 17, height = 10, dpi = 600)

ggsave(filename = "Figure 3.emf",
       path = "C:/mypath/Figures/",
       width = 17, height = 10, dpi = 600)




################################################################################
####################### Sensitivity analyses ###################################
################################################################################

################################################################################
#####>>>>>>>>>>>>> LTMLE with censoring nodes (IPCW) (estimand: always vs never)
#####>>>>>>>>>>>>> Population: eligible patients with complete BL data
################################################################################

###################>>>>>> Not used to export data

###################>>>>>> Prepare the data

setwd("C:/mypath/Datasets/")
d <- read.csv("meteor12ltmlemocondswitchMI.csv")
names(d)[names(d) == "switchnextt2"] <- "L1s"


###################>>>>>> Numbers for the flow diagram (Supplementary Figure S2)

basecov  <- c("mny","asasmri","hla","sex","age","comorbbin",
              "pertvt1","ibdbl","emmtvt1","comedtvt1","L0d","L0e")
scores   <- c("L0d","L0e","L1d","L1e","L2d","Y2")   # ASDAS and BASFI, 3 visits
othercov <- setdiff(basecov, scores)

cat("\n--- SUPPLEMENTARY FIGURE S2 ---\n")
cat("Eligible patients (3 visits, bDMARD-naive, high disease activity):", nrow(d), "\n")

##### Main analysis: complete case
cc      <- complete.cases(d[, scores])            # complete ASDAS + BASFI, 3 visits
cc_full <- cc & complete.cases(d[, othercov])     # + complete baseline covariates
cat("\n[Main analysis, complete case]\n")
cat("  excluded, missing ASDAS or BASFI in >=1 visit :", sum(!cc), "\n")
cat("  complete ASDAS and BASFI at all three visits  :", sum(cc), "\n")
cat("  excluded, missing a baseline covariate        :", sum(cc & !cc_full), "\n")
cat("  ANALYSED (g-formula, MSM, LTMLE, GEE)         :", sum(cc_full), "\n")
cat("    elevated CRP / normal CRP                   :",
    sum(d$crpelevatedt1[cc_full] == 1, na.rm = TRUE), "/",
    sum(d$crpelevatedt1[cc_full] == 0, na.rm = TRUE), "\n")

##### Sensitivity analysis: censoring-weighted
ok   <- complete.cases(d[, basecov])
noB  <- is.na(d$L0e)
noC  <- !complete.cases(d[, setdiff(basecov, "L0e")])
s    <- d[ok, ]
u6   <- !is.na(s$L1d) & !is.na(s$L1e)
u12  <- u6 & !is.na(s$L2d) & !is.na(s$Y2)
late <- !u6 & !is.na(s$L2d) & !is.na(s$Y2)        # 12M data despite 6M censoring

cat("\n[Sensitivity analysis, censoring-weighted]\n")
cat("  excluded, missing baseline BASFI              :", sum(noB), "\n")
cat("  excluded, missing a baseline covariate        :", sum(noC & !noB), "\n")
cat("  excluded, total                               :", sum(noB | noC), "\n")
cat("  complete baseline data                        :", nrow(s), "\n")
cat("    both 6M and 12M                             :", sum(u12), "\n")
cat("    6M only                                     :", sum(u6 & !u12), "\n")
cat("    baseline only                               :", sum(!u6), "\n")
cat("      of which 12M recorded but 6M missing      :", sum(late), "\n")
cat("  ANALYSED (LTMLE with censoring weights)       :", nrow(s), "\n")
cat("    elevated CRP / normal CRP                   :",
    sum(s$crpelevatedt1 == 1, na.rm = TRUE), "/",
    sum(s$crpelevatedt1 == 0, na.rm = TRUE), "\n\n")

## internal consistency
stopifnot(sum(!cc) + sum(cc & !cc_full) + sum(cc_full) == nrow(d),
          sum(noB | noC) + nrow(s) == nrow(d),
          sum(u12) + sum(u6 & !u12) + sum(!u6) == nrow(s))

d <- s   # keep last: everything above needs the full 761


###################>>>>>> Censoring nodes

## uncensored = that visit's required nodes observed; monotonicity enforced
unc6  <- !is.na(d$L1d) & !is.na(d$L1e)
unc12 <- !is.na(d$L2d) & !is.na(d$Y2)
unc12[!unc6] <- FALSE

d$C1 <- BinaryToCensoring(is.censored = !unc6)
d$C2 <- BinaryToCensoring(is.censored = !unc12)

## all data after a censoring event must be NA
d[!unc6,  c("L1a","L1b","L1c","L1d","L1e","L1s","A1",
            "L2a","L2b","L2c","L2d","Y2")] <- NA
d[!unc12, c("L2a","L2b","L2c","L2d","Y2")] <- NA
d$C2[!unc6] <- NA


###################>>>>>> Strategies and models

never  <- function(row) c(0, 0)
always <- function(row) c(1, 1)

## column order = causal order; C1 precedes the 6M nodes it censors.
## id is NOT included: any column before the first A node is treated by ltmle
## as a baseline covariate.
nc <- c(basecov, "A0", "C1", "L1a","L1b","L1c","L1d","L1e","L1s",
        "A1", "C2", "L2a","L2b","L2c","L2d", "Y2")
Lnodes <- c("L1a","L1b","L1c","L1d","L1e","L1s","L2a","L2b","L2c","L2d")
Anodes <- c("A0","A1")
Cnodes <- c("C1","C2")
Ynodes <- c("Y2")

base <- "mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e"

Qform <- c("L1a" = paste("Q.kplus1 ~", base, "+ A0"),
           "L2a" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"))
gform <- c("A0" = paste("A0 ~", base),
           "C1" = paste("C1 ~", base, "+ A0"),
           "A1" = "A1 ~ asasmri + hla + age + sex + comorbbin + ibdbl + A0*mny*L1d + L1a + L1b + L1c + L1d + L1e",
           "C2" = paste("C2 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"))


###################>>>>>> Estimation

## Pooled SD for the e-value. Computed after the censoring block, so it is taken
## over the 352 patients whose outcome survives monotone censoring - the same
## 352 as the main analysis, which puts the S10 e-values on the Table 2 scale.
## Computing it before the block would add the 26 late-recorded outcomes.
sd_y <- sd(d$Y2, na.rm = TRUE)
stopifnot(sum(!is.na(d$Y2)) == 352)

evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * est / sd_y)
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}

###################>>>>>> Diagnostics

trunc_pct <- function(f) {
  100 * mean(apply(f$cum.g != f$cum.g.unbounded, 1,
                   function(x) any(x, na.rm = TRUE)))
}

sd_ic <- function(f) {
  IC <- f$IC
  if (is.matrix(IC) && ncol(IC) == 2) IC <- IC[, 1] - IC[, 2]
  sd(IC)
}

run <- function(z, lab, sd_y) {
  f  <- ltmle(z[, nc], Anodes = Anodes, Cnodes = Cnodes, Lnodes = Lnodes,
              Ynodes = Ynodes, Qform = Qform, gform = gform,
              rule = list(always, never), estimate.time = FALSE)
  em <- summary(f)$effect.measures
  a  <- em$ATE; p11 <- em$treatment; p00 <- em$control
  n  <- nrow(z)
  
  cat(sprintf("%-12s n=%3d  ATE=%+.3f (%.3f; %.3f)  SE=%.3f\n",
              lab, n, a$estimate, a$CI[1], a$CI[2], a$std.dev))
  
  list(
    label        = lab,
    n            = n,
    est          = a$estimate, se = a$std.dev, lo = a$CI[1], hi = a$CI[2],
    po11         = p11$estimate, se11 = p11$std.dev,
    po00         = p00$estimate, se00 = p00$std.dev,
    evalue       = evalue_md(a$estimate, sd_y),
    followers    = sum(z$A0 == 1 & z$A1 == 1, na.rm = TRUE),
    n_unc6       = sum(z$C1 == "uncensored"),
    cens6        = sum(z$C1 == "censored"),
    cens12       = sum(z$C2 == "censored", na.rm = TRUE),
    trunc        = trunc_pct(f),
    w_mean       = mean(1 / f$cum.g, na.rm = TRUE),
    w_max        = max(1 / f$cum.g, na.rm = TRUE),
    sd_ic        = sd_ic(f),
    fit          = f)
}

cat("--- IPCW LTMLE: always vs never treated, BASFI at 12 months ---\n")
r <- list(
  run(d,                                "All",          sd_y),
  run(d[which(d$crpelevatedt1 == 1), ], "CRP elevated", sd_y),
  run(d[which(d$crpelevatedt1 == 0), ], "CRP normal",   sd_y))
names(r) <- sapply(r, `[[`, "label")


###################>>>>>> Supplementary Table S10 (effect estimates)

s10 <- do.call(rbind, lapply(r, function(x) data.frame(
  Population        = x$label,
  N                 = x$n,
  Treatment_effect  = sprintf("%.2f (%.2f; %.2f)", x$est, x$lo, x$hi),
  SE                = sprintf("%.2f", x$se),
  Always_treated_PO = sprintf("%.2f (%.2f)", x$po11, x$se11),
  Never_treated_PO  = sprintf("%.2f (%.2f)", x$po00, x$se00),
  e_value           = sprintf("%.1f", x$evalue),
  stringsAsFactors  = FALSE)))
rownames(s10) <- NULL
print(s10)


###################>>>>>> Diagnostics table

## Followers is reported two ways: out of the whole IPCW population, and out of
## the patients still uncensored at 6 months. Only the second is comparable with
## the 102 (29%) in S4, because a patient censored at 6M has A1 missing and so
## can never be counted as a follower of "always treated".
s10diag <- do.call(rbind, lapply(r, function(x) data.frame(
  Population         = x$label,
  N                  = x$n,
  Followers          = sprintf("%d (%.0f)", x$followers, 100 * x$followers / x$n),
  Followers_of_unc6  = sprintf("%d (%.0f)", x$followers, 100 * x$followers / x$n_unc6),
  Censored_6M        = sprintf("%d (%.0f)", x$cens6,  100 * x$cens6  / x$n),
  Censored_12M       = sprintf("%d (%.0f)", x$cens12, 100 * x$cens12 / x$n),
  Truncation         = sprintf("%.1f", x$trunc),
  Weights            = sprintf("%.2f (max %.1f)", x$w_mean, x$w_max),
  Influence_curve    = sprintf("%.1f", x$sd_ic),
  stringsAsFactors   = FALSE)))
rownames(s10diag) <- NULL
print(s10diag)


################################################################################
#####>>>>>>>>>>>>> LTMLE, alternative definitions of "treated with a bDMARD"
#####>>>>>>>>>>>>> Population: completers (N = 352), same as the main analysis
################################################################################

## Nodes, Q models, g models and the two strategies are identical to the main
## analysis. Only the exposure columns A0 and A1 are redefined, so any change in
## the estimate is attributable to the definition of treatment. 

###################>>>>>> Not used to export data


###################>>>>>> Prepare the data
setwd("C:/mypath/Datasets/")
WIDE <- "meteor12ltmlemocondswitch.csv"    # one row per patient, n = 352
LONG <- "originalfulllong761.csv"          # long file, only for the extra exposure columns


w <- read.csv(WIDE, header = TRUE, sep = ",")
names(w)[names(w) == "switchnextt2"] <- "L1s"
w <- w %>% relocate(L1s, .before = A1)

## Column order is the causal order. Anything merged in below must be kept out
## of the ltmle call; node.cols is what does that. id is excluded because any
## column preceding the first A node is treated by ltmle as a baseline covariate.
node.cols <- c("mny","asasmri","crpelevatedt1","hla","sex","age","comorbbin",
               "pertvt1","ibdbl","emmtvt1","comedtvt1","L0d","L0e","A0",
               "L1a","L1b","L1c","L1d","L1e","L1s","A1","L2a","L2b","L2c","L2d","Y2")
stopifnot(all(node.cols %in% names(w)))


###################>>>>>> Alternative exposure definitions, from the long file

## bionew1mnext : bDMARD with >= 1 month of exposure in the interval
## bionew3mnext : bDMARD with >= 3 months of exposure in the interval
## expdays1/2   : exposurebio1interval at the 6M and 12M visits, i.e. days of
##                bDMARD exposure within interval 1 (BL-6M) and interval 2 (6M-12M)
## bionew6mnext is excluded: only 11 of the 352 are always treated
## under it, which is too few to estimate the always-treated potential outcome.

long <- read.csv(LONG)

alt  <- subset(long, t == 0,
               select = c("id","bionewnextt1","bionewnextt2",
                          "bionew1mnextt1","bionew1mnextt2",
                          "bionew3mnextt1","bionew3mnextt2"))
exp1 <- subset(long, t == 6,  select = c("id","exposurebio1interval"))
exp2 <- subset(long, t == 12, select = c("id","exposurebio1interval"))
names(exp1)[2] <- "expdays1"
names(exp2)[2] <- "expdays2"

extra <- merge(merge(alt, exp1, by = "id"), exp2, by = "id")
n0    <- nrow(w)
w     <- merge(w, extra, by = "id", all.x = TRUE, sort = FALSE)

## The merge must not change the sample, and the long file's primary exposure
## must reproduce A0/A1 exactly. If it does not, the "next" convention differs
## between the two files and every alternative definition is misaligned.
stopifnot(nrow(w) == n0,
          !anyNA(w[, setdiff(names(extra), "id")]),
          all(w$A0 == w$bionewnextt1),
          all(w$A1 == w$bionewnextt2))
cat("Merged from the long file:", setdiff(names(extra), "id"), "\n")


###################>>>>>> Strategies and models

never  <- function(row) c(0, 0)
always <- function(row) c(1, 1)

Lnodes <- c("L1a","L1b","L1c","L1d","L1e","L1s","L2a","L2b","L2c","L2d")
Anodes <- c("A0","A1")
Ynodes <- c("Y2")

## ltmle groups consecutive L/Y nodes not separated by a treatment node into one
## block and fits a single Q regression per block, so only the L1a and L2a
## formulas are used. The rest are kept because they state the assumed structure
## and mirror Supplementary Table S1; ltmle drops them with a message.
base <- "mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e"

Qform <- c("L1a" = paste("Q.kplus1 ~", base, "+ A0"),
           "L1b" = paste("Q.kplus1 ~", base, "+ A0"),
           "L1c" = paste("Q.kplus1 ~", base, "+ A0"),
           "L1d" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c"),
           "L1e" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d"),
           "L1s" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e"),
           "L2a" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"),
           "L2b" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"),
           "L2c" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"),
           "L2d" = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1 + L2a + L2b + L2c"),
           "Y2"  = paste("Q.kplus1 ~", base, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1 + L2a + L2b + L2c + L2d"))

gform <- c("A0" = paste("A0 ~", base),
           "A1" = "A1 ~ asasmri + hla + age + sex + comorbbin + ibdbl + A0*mny*L1d + L1a + L1b + L1c + L1d + L1e")


###################>>>>>> Estimation

sd_y <- sd(w$Y2, na.rm = TRUE)

evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * est / sd_y)
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}

## Diagnostics defined exactly as in the Supplementary Table S4 script.
trunc_pct <- function(f) {
  100 * mean(apply(f$cum.g != f$cum.g.unbounded, 1,
                   function(x) any(x, na.rm = TRUE)))
}

sd_ic <- function(f) {
  IC <- f$IC
  if (is.matrix(IC) && ncol(IC) == 2) IC <- IC[, 1] - IC[, 2]
  sd(IC)
}

run_def <- function(a_cols, lab, sd_y) {
  z     <- w
  z$A0  <- z[[a_cols[1]]]
  z$A1  <- z[[a_cols[2]]]
  
  f  <- ltmle(z[, node.cols], Anodes = Anodes, Lnodes = Lnodes, Ynodes = Ynodes,
              Qform = Qform, gform = gform,
              rule = list(always, never), estimate.time = FALSE)
  em <- summary(f)$effect.measures
  a  <- em$ATE; p11 <- em$treatment; p00 <- em$control
  n  <- nrow(z)
  
  alw <- z$A0 == 1 & z$A1 == 1
  nev <- z$A0 == 0 & z$A1 == 0
  
  cat(sprintf("%-32s n=%3d  ATE=%+.3f (%.3f; %.3f)  SE=%.3f\n",
              lab, n, a$estimate, a$CI[1], a$CI[2], a$std.dev))
  
  list(
    label     = lab,
    n         = n,
    est       = a$estimate, se = a$std.dev, lo = a$CI[1], hi = a$CI[2],
    po11      = p11$estimate, se11 = p11$std.dev,
    po00      = p00$estimate, se00 = p00$std.dev,
    evalue    = evalue_md(a$estimate, sd_y),
    followers = sum(alw),
    n_never   = sum(nev),
    ## Median days of bDMARD within each interval, always vs never treated.
    ## Interval 1 is a partial-exposure period by construction (the drug is
    ## started part-way through it), so the two intervals are not summed.
    dose      = sprintf("%.0f / %.0f vs %.0f / %.0f",
                        median(z$expdays1[alw]), median(z$expdays2[alw]),
                        median(z$expdays1[nev]), median(z$expdays2[nev])),
    trunc     = trunc_pct(f),
    sd_ic     = sd_ic(f),
    fit       = f)
}

defs <- list(
  "Any bDMARD recorded (primary)" = c("A0", "A1"),
  "At least 1 month of exposure"  = c("bionew1mnextt1", "bionew1mnextt2"),
  "At least 3 months of exposure" = c("bionew3mnextt1", "bionew3mnextt2"))

cat("--- LTMLE: always vs never treated, BASFI at 12 months, by exposure definition ---\n")
rd <- Map(run_def, defs, names(defs), MoreArgs = list(sd_y = sd_y))


###################>>>>>> Rows for Supplementary Table S10

s10_def <- do.call(rbind, lapply(rd, function(x) data.frame(
  Analysis          = x$label,
  N                 = x$n,
  Treatment_effect  = sprintf("%.2f (%.2f; %.2f)", x$est, x$lo, x$hi),
  SE                = sprintf("%.2f", x$se),
  Always_treated_PO = sprintf("%.2f (%.2f)", x$po11, x$se11),
  Never_treated_PO  = sprintf("%.2f (%.2f)", x$po00, x$se00),
  e_value           = sprintf("%.1f", x$evalue),
  stringsAsFactors  = FALSE)))
rownames(s10_def) <- NULL
print(s10_def)


###################>>>>>> Rows for the S10 diagnostics table

s10_def_diag <- do.call(rbind, lapply(rd, function(x) data.frame(
  Analysis         = x$label,
  N                = x$n,
  Always_treated   = sprintf("%d (%.0f)", x$followers, 100 * x$followers / x$n),
  Never_treated    = sprintf("%d (%.0f)", x$n_never,   100 * x$n_never   / x$n),
  Exposure_days    = x$dose,
  Truncation       = sprintf("%.1f", x$trunc),
  Influence_curve  = sprintf("%.1f", x$sd_ic),
  stringsAsFactors = FALSE)))
rownames(s10_def_diag) <- NULL
print(s10_def_diag)



################################################################################
####################### Sensitivity analyses ###################################
####################  Supplementary  Table S10 #################################
################################################################################

################################################################################
#####>>>>>>>>>>>>> All strategies (static + dynamic), stratified on CRP, under
#####>>>>>>>>>>>>> the two sensitivity analyses:
#####>>>>>>>>>>>>>   A. censoring weights (IPCW), N = 452
#####>>>>>>>>>>>>>   B. alternative definitions of "treated", N = 352
################################################################################

## Strategies, Q models and g models are identical to Figure 3. Analysis A adds
## censoring nodes and widens the population; analysis B redefines A0 and A1
## only. 84 fits (A: 7 strategies x 3 subgroups; B: 7 x 3 x 3 definitions).
##
## The primary-definition rows of analysis B reproduce Supplementary Table S9
## and Figure 3 exactly and are kept as the reference row of each block.

setwd("C:/mypath/Datasets/")

WIDE <- "meteor12ltmlemocondswitch.csv"     # completers, n = 352
MI   <- "meteor12ltmlemocondswitchMI.csv"   # all eligible, n = 761
LONG <- "originalfulllong761.csv"           # long file, for the extra exposure columns

basecov <- c("mny","asasmri","hla","sex","age","comorbbin",
             "pertvt1","ibdbl","emmtvt1","comedtvt1","L0d","L0e")


###################>>>>>> Strategies

## Rules are exactly as in Figure 3. 
STRAT <- list(
  list(lab = "Static (always treated)",      sw = FALSE, f = function(r) c(1, 1)),
  list(lab = "Dynamic 1 (stop if no LDA)",   sw = FALSE, f = function(r) c(1, as.numeric(r[["L1d"]] <  2.1))),
  list(lab = "Dynamic 2 (stop if no ID)",    sw = FALSE, f = function(r) c(1, as.numeric(r[["L1d"]] <  1.3))),
  list(lab = "Dynamic 3 (switch if no LDA)", sw = TRUE,  f = function(r) c(1, as.numeric(r[["L1d"]] >= 2.1), 1)),
  list(lab = "Dynamic 4 (switch if no ID)",  sw = TRUE,  f = function(r) c(1, as.numeric(r[["L1d"]] >= 1.3), 1)),
  list(lab = "Dynamic 5 (stop if no CII)",   sw = FALSE, f = function(r) c(1, as.numeric(r[["L1d"]] - r[["L0d"]] <  -1.1))),
  list(lab = "Dynamic 6 (switch if no CII)", sw = TRUE,  f = function(r) c(1, as.numeric(r[["L1d"]] - r[["L0d"]] >= -1.1), 1)))


###################>>>>>> Models

B <- "mny + asasmri + hla + age + sex + comorbbin + pertvt1 + ibdbl + emmtvt1 + comedtvt1 + L0d + L0e"

## ltmle groups consecutive L/Y nodes into one block and fits one Q regression
## per block, so only L1a and L2a are fitted. The block boundaries are the same
## whether L1s is an A node or an L node, and with or without C nodes.
Qform <- c("L1a" = paste("Q.kplus1 ~", B, "+ A0"),
           "L2a" = paste("Q.kplus1 ~", B, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1"))

gA0 <- paste("A0 ~", B)
gA1 <- "A1 ~ asasmri + hla + age + sex + comorbbin + ibdbl + A0*mny*L1d + L1a + L1b + L1c + L1d + L1e"
gSW <- paste("L1s ~", B, "+ L1a + L1b + L1c + L1d + L1e")
gC1 <- paste("C1 ~", B, "+ A0")
gC2 <- paste("C2 ~", B, "+ A0 + L1a + L1b + L1c + L1d + L1e + L1s + A1")

GB <- c(0.01, 0.99)   # gbounds

## Only a patient who actually received a bDMARD can switch. A0 is checked
## because ltmle sets it to the counterfactual value, so the constraint lifts
## under a strategy that assigns treatment; anyb0 is the primary exposure
## indicator, which stays fixed when A0 is redefined in analysis B.

det.g <- function(data, current.node, nodes) {
  if (!"L1s" %in% names(data)) return(NULL)
  if (current.node != which(names(data) == "L1s")) return(NULL)
  list(is.deterministic = data$A0 == 0 & data$anyb0 == 0, prob1 = 0)
}

###################>>>>>> Helpers

## E-value using the SD of the outcome in the analysed sample. This is the
## convention behind Table 2 and Figure 3: the SDs differ by subgroup 
evalue_md <- function(est, sd_y) {
  rr <- exp(0.91 * abs(est / sd_y))
  rr + sqrt(rr * (rr - 1))
}

## apply() over a frame containing the C1/C2 censoring factors coerces every
## column to character, which breaks the arithmetic in the Dynamic 5 and 6
## rules. Loop over rows instead, as ltmle does internally.
apply_rule <- function(dat, f, k)
  t(vapply(seq_len(nrow(dat)), function(i) as.numeric(f(dat[i, ])), numeric(k)))


###################>>>>>> One strategy in one sample

run_one <- function(dat, s, cols, ipcw = FALSE) {
  
  Anodes <- if (s$sw) c("A0","L1s","A1") else c("A0","A1")
  Lnodes <- if (s$sw) c("L1a","L1b","L1c","L1d","L1e","L2a","L2b","L2c","L2d")
  else      c("L1a","L1b","L1c","L1d","L1e","L1s","L2a","L2b","L2c","L2d")
  Cnodes <- if (ipcw) c("C1","C2") else NULL
  
  gform <- c("A0" = gA0)
  if (ipcw) gform <- c(gform, "C1" = gC1)
  if (s$sw) gform <- c(gform, "L1s" = gSW)
  gform <- c(gform, "A1" = gA1)
  if (ipcw) gform <- c(gform, "C2" = gC2)
  
  ref <- if (s$sw) function(r) c(0, 0, 0) else function(r) c(0, 0)
  
  args <- list(data = as.data.frame(dat)[, cols], Anodes = Anodes, Cnodes = Cnodes,
               Lnodes = Lnodes, Ynodes = "Y2", Qform = Qform, gform = gform,
               rule = list(s$f, ref), gbounds = GB, estimate.time = FALSE)
  if (s$sw) args$deterministic.g.function <- det.g
  f <- suppressWarnings(do.call(ltmle, args))
  
  e <- summary(f)$effect.measures
  n <- nrow(dat)
  
  ## Followers: rule applied row by row against what was observed. Censored
  ## patients have NA in A1 / L1s and count as non-followers, so the percentage
  ## is also given against those still uncensored at 6 months.
  ab  <- apply_rule(as.data.frame(dat)[, cols], s$f, if (s$sw) 3 else 2)
  hit <- if (s$sw) (dat$A0 == ab[, 1]) & (dat$L1s == ab[, 2]) & (dat$A1 == ab[, 3])
  else      (dat$A0 == ab[, 1]) & (dat$A1 == ab[, 2])
  foll  <- sum(hit %in% TRUE)
  n_obs <- if (ipcw) sum(dat$C1 == "uncensored") else n
  
  ## Truncation, two ways: strategy arm alone, and the union of both arms.
  ## na.rm is required here - cum.g carries NA after a censoring event.
  tm   <- f$cum.g != f$cum.g.unbounded
  anyT <- function(x) any(x, na.rm = TRUE)
  
  data.frame(
    strategy = s$lab, n = n, n_obs = n_obs,
    ate = e$ATE$estimate, se = e$ATE$std.dev, lo = e$ATE$CI[1], hi = e$ATE$CI[2],
    po_strategy = e$treatment$estimate, se_strategy = e$treatment$std.dev,
    po_never    = e$control$estimate,   se_never    = e$control$std.dev,
    followers   = foll,
    trunc_strategy_pct = 100 * mean(apply(tm[, , 1, drop = FALSE], 1, anyT)),
    trunc_either_pct   = 100 * mean(apply(tm, 1, anyT)),
    sd_IC_asS4 = { IC <- f$IC; sd(IC[, 1] - IC[, 2]) },
    sd_IC_ate  = e$ATE$std.dev * sqrt(n),
    evalue     = evalue_md(e$ATE$estimate, sd(dat$Y2, na.rm = TRUE)),
    stringsAsFactors = FALSE)
}


################################################################################
#####>>>>>> ANALYSIS A. Censoring weights (IPCW)
################################################################################

a <- read.csv(MI)
names(a)[names(a) == "switchnextt2"] <- "L1s"
a <- a[complete.cases(a[, basecov]), ]
stopifnot(nrow(a) == 452)
a$anyb0 <- a$A0

unc6  <- !is.na(a$L1d) & !is.na(a$L1e)
unc12 <- !is.na(a$L2d) & !is.na(a$Y2)
unc12[!unc6] <- FALSE

a$C1 <- BinaryToCensoring(is.censored = !unc6)
a$C2 <- BinaryToCensoring(is.censored = !unc12)
a[!unc6,  c("L1a","L1b","L1c","L1d","L1e","L1s","A1",
            "L2a","L2b","L2c","L2d","Y2")] <- NA
a[!unc12, c("L2a","L2b","L2c","L2d","Y2")] <- NA
a$C2[!unc6] <- NA

## Causal order; C1 precedes the 6M nodes it censors. id is excluded because any
## column before the first A node is treated by ltmle as a baseline covariate.
cols_A <- c(basecov, "crpelevatedt1", "anyb0", "A0", "C1",
            "L1a","L1b","L1c","L1d","L1e","L1s", "A1", "C2",
            "L2a","L2b","L2c","L2d", "Y2")

SUB_A <- list("All"  = a,
              "CRP+" = a[which(a$crpelevatedt1 == 1), ],
              "CRP-" = a[which(a$crpelevatedt1 == 0), ])

cat("--- A. IPCW, all strategies ---\n")
resA <- do.call(rbind, lapply(names(SUB_A), function(g)
  do.call(rbind, lapply(STRAT, function(s) {
    cat(sprintf("  %-5s %-30s\n", g, s$lab))
    cbind(analysis = "IPCW (N=452)", subgroup = g,
          run_one(SUB_A[[g]], s, cols_A, ipcw = TRUE))
  }))))


################################################################################
#####>>>>>> ANALYSIS B. Alternative definitions of "treated"
################################################################################

b <- read.csv(WIDE)
names(b)[names(b) == "switchnextt2"] <- "L1s"

long <- read.csv(LONG)
alt  <- subset(long, t == 0,
               select = c("id","bionewnextt1","bionewnextt2",
                          "bionew1mnextt1","bionew1mnextt2",
                          "bionew3mnextt1","bionew3mnextt2"))
n0 <- nrow(b)
b  <- merge(b, alt, by = "id", all.x = TRUE, sort = FALSE)

## The long file's primary exposure must reproduce A0/A1 exactly, or the "next"
## convention differs between the files and every alternative definition is
## shifted by one interval.
stopifnot(nrow(b) == n0, all(b$A0 == b$bionewnextt1), all(b$A1 == b$bionewnextt2))

cols_B <- c(basecov, "crpelevatedt1", "anyb0", "A0",
            "L1a","L1b","L1c","L1d","L1e","L1s","A1",
            "L2a","L2b","L2c","L2d","Y2")

DEFS <- list(
  "Any bDMARD recorded (primary)" = c("A0", "A1"),
  "At least 1 month of exposure"  = c("bionew1mnextt1", "bionew1mnextt2"),
  "At least 3 months of exposure" = c("bionew3mnextt1", "bionew3mnextt2"))

cat("--- B. Alternative treatment definitions, all strategies ---\n")
resB <- do.call(rbind, lapply(names(DEFS), function(dn) {
  z    <- b
  z$anyb0 <- z$bionewnextt1          # primary exposure, before A0 is redefined
  z$A0 <- z[[DEFS[[dn]][1]]]
  z$A1 <- z[[DEFS[[dn]][2]]]
  SUB  <- list("All"  = z,
               "CRP+" = z[which(z$crpelevatedt1 == 1), ],
               "CRP-" = z[which(z$crpelevatedt1 == 0), ])
  do.call(rbind, lapply(names(SUB), function(g)
    do.call(rbind, lapply(STRAT, function(s) {
      cat(sprintf("  %-30s %-5s %-30s\n", dn, g, s$lab))
      cbind(analysis = dn, subgroup = g, run_one(SUB[[g]], s, cols_B, ipcw = FALSE))
    }))))
}))

res <- rbind(resA, resB)


###################>>>>>> Rows for Supplementary Table S10

s10_strat <- data.frame(
  Analysis         = res$analysis,
  Subgroup         = res$subgroup,
  Strategy         = res$strategy,
  N                = res$n,
  Treatment_effect = sprintf("%.2f (%.2f; %.2f)", res$ate, res$lo, res$hi),
  SE               = sprintf("%.2f", res$se),
  Strategy_PO      = sprintf("%.2f (%.2f)", res$po_strategy, res$se_strategy),
  Never_treated_PO = sprintf("%.2f (%.2f)", res$po_never,    res$se_never),
  e_value          = sprintf("%.1f", res$evalue),
  stringsAsFactors = FALSE)

###################>>>>>> Rows for the S10 diagnostics table

s10_strat_diag <- data.frame(
  Analysis            = res$analysis,
  Subgroup            = res$subgroup,
  Strategy            = res$strategy,
  N                   = res$n,
  Followers           = sprintf("%d (%.0f)", res$followers, 100 * res$followers / res$n),
  Followers_of_obs    = sprintf("%d (%.0f)", res$followers, 100 * res$followers / res$n_obs),
  Truncation          = sprintf("%.0f", res$trunc_strategy_pct),
  Truncation_either   = sprintf("%.0f", res$trunc_either_pct),
  Influence_curve     = sprintf("%.1f", res$sd_IC_asS4),
  Influence_curve_ATE = sprintf("%.1f", res$sd_IC_ate),
  stringsAsFactors    = FALSE)

print(s10_strat)
print(s10_strat_diag)



###################>>>>>> Export

setwd("C:/mypath/Tables/")
write_xlsx(list(S10_strategies             = s10_strat,
                S10_strategies_diagnostics = s10_strat_diag,
                numeric                    = res),
           "TableS10sensanalysis.xlsx")
