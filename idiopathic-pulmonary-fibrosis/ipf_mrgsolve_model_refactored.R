## =============================================================================
## IPF QSP Model — mrgsolve ODE Implementation — REFACTORED (pluggable PK)
## Disease: Idiopathic Pulmonary Fibrosis (IPF)
## Drugs  : Pirfenidone (PIR, ESBRIET) · Nintedanib (NIN, OFEV) · Combination
## Author : Claude Code Routine (CCR) · 2026-06-17; PK/PD refactor 2026-08-31
##
## Derived from ipf_mrgsolve_model.R (untouched, upstream-tracked). Every
## compound's PK is reorganized into its own clearly-delimited block, renamed
## to the fork's `<ROLE>_<STEM>` convention (FORK_WORKFLOW_GUIDE.md Part 2),
## and exposes exactly one `C_<STEM>` concentration plus its own named
## `EFFECT_<STEM>*` disease-effect term(s). See ipf_refactor_notes.md for the
## full account (naming map, build-compat fix, verification result).
##
## Key References (calibration anchors, unchanged from the original):
##  - ASCEND (King et al. NEJM 2014): pirfenidone –47.9% FVC decline vs placebo
##  - INPULSIS-1/2 (Richeldi et al. NEJM 2014): nintedanib –50.1% FVC decline
##  - Pirfenidone PK: Rubino (2009) Clin Pharmacokinet; F=81%, t½=2.4h, ka=1.74/h
##  - Nintedanib PK: Stopfer (2011) Clin Pharmacokinet; F=4.7%, t½=10h, ka=0.8/h
##  - Natural history FVC decline ~200 mL/yr (Raghu 2011 ATS Guidelines)
## =============================================================================

library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(gridExtra)

## ─── 1. MODEL DEFINITION (refactored) ───────────────────────────────────────

ipf_refactored_code <- '
$PROB IPF QSP Model — Pirfenidone (PIR) & Nintedanib (NIN) PK/PD [refactored: pluggable-PK naming convention]

$PARAM @annotated
// ── Pirfenidone (PIR) PK ─────────────────────────────────────────────────────
KA_PIR : 1.74   : Absorption rate constant pirfenidone (1/h)
F_PIR  : 0.81   : Oral bioavailability pirfenidone
CL_PIR : 8.4    : Clearance pirfenidone (L/h)
V1_PIR : 20.4   : Central volume pirfenidone (L)
V2_PIR : 14.0   : Peripheral volume pirfenidone (L)
Q_PIR  : 3.6    : Inter-compartment clearance pirfenidone (L/h)

// ── Nintedanib (NIN) PK ──────────────────────────────────────────────────────
KA_NIN : 0.80   : Absorption rate constant nintedanib (1/h)
F_NIN  : 0.047  : Oral bioavailability nintedanib
CL_NIN : 22.0   : Clearance nintedanib (L/h)
V1_NIN : 730.0  : Central volume nintedanib (L)
V2_NIN : 900.0  : Peripheral volume nintedanib (L)
Q_NIN  : 15.0   : Inter-compartment clearance nintedanib (L/h)

// ── Disease PD — AEC & TGF-β ─────────────────────────────────────────────────
AEC2_ss : 1.0   : AEC-II steady-state (normalized)
k_AEC   : 0.005 : AEC-II damage rate (per hour, from ROS/injury)
k_rep   : 0.003 : AEC-II repair rate (per hour)
kprod_TGFb : 0.08 : TGF-β1 basal production rate
kdeg_TGFb  : 0.04 : TGF-β1 degradation rate (1/h)
kact_TGFb  : 0.12 : TGF-β1 activation from damaged AEC (proportionality)
EC50_PIR_TGFB : 30.0 : Pirfenidone EC50 for TGF-β inhibition (µg/mL)
EMAX_PIR_TGFB : 0.65 : Pirfenidone Emax for TGF-β inhibition
GAMMA_PIR_TGFB : 1.0 : Hill coefficient for pirfenidone TGF-β inhibition [new; original had none, math-implied =1]

// ── Disease PD — Fibroblast/Myofibroblast ─────────────────────────────────────
kact_F  : 0.015 : Fibroblast activation rate by TGF-β
kdiff_M : 0.025 : Myofibroblast differentiation rate from activated fibroblast
kapop_F : 0.008 : Fibroblast apoptosis (reduced in IPF)
F_ss    : 1.0   : Baseline fibroblast level
EC50_NIN : 15.0  : Nintedanib EC50 for fibroblast proliferation inhibition (nM)
EMAX_NIN : 0.70  : Nintedanib Emax for fibroblast inhibition
GAMMA_NIN : 1.0  : Hill coefficient for nintedanib fibroblast inhibition [new; original had none, math-implied =1]

// ── Disease PD — ECM / Collagen ───────────────────────────────────────────────
kprod_Col : 0.010 : Collagen production rate by myofibroblasts (per unit time)
kdeg_Col  : 0.002 : Collagen degradation (MMP-mediated)
kprod_MMP : 0.020 : MMP production rate
kdeg_MMP  : 0.030 : MMP degradation rate
kprod_TIMP: 0.018 : TIMP production rate (TGF-β driven)
kdeg_TIMP : 0.020 : TIMP degradation rate
Col_ss    : 1.0   : Baseline collagen (normalized)

// ── Disease PD — Macrophage ───────────────────────────────────────────────────
kM2_act   : 0.012 : M2 macrophage activation rate
kdeg_M2   : 0.008 : M2 macrophage clearance rate
EC50_PIR_M2 : 25.0  : Pirfenidone EC50 for M2/cytokine inhibition (µg/mL)
GAMMA_PIR_M2 : 1.0 : Hill coefficient for pirfenidone M2 inhibition [new; original had none, math-implied =1]

// ── Oxidative stress ──────────────────────────────────────────────────────────
kprod_ROS  : 0.04 : Basal ROS production rate
kdeg_ROS   : 0.05 : ROS clearance rate (GSH/antioxidant)
kfb_ROS    : 0.02 : ROS positive feedback from AEC damage
EMAX_PIR_ROS : 0.45 : Pirfenidone Emax antioxidant effect
GAMMA_PIR_ROS : 1.0 : Hill coefficient for pirfenidone ROS inhibition [new; original had none, math-implied =1]

// ── Lung Function (FVC) ───────────────────────────────────────────────────────
FVC_base   : 80.0 : Baseline FVC % predicted (typical mild-moderate IPF)
k_FVC_loss : 0.0026 : FVC loss rate per unit collagen excess (per hour)
// ~200 mL/yr ≈ 2.5% predicted/yr ≈ 0.000285%/h → calibrated

// ── DLCO ─────────────────────────────────────────────────────────────────────
DLCO_base  : 65.0  : Baseline DLCO % predicted
k_DLCO_loss: 0.0018 : DLCO loss rate from AEC damage

// ── Molecular weights for unit conversion ────────────────────────────────────
MW_PIR : 185.22 : Molecular weight pirfenidone (g/mol) [unused downstream in original -- preserved, not repurposed]
MW_NIN : 539.63 : Molecular weight nintedanib (g/mol)

$CMT @annotated
// Pirfenidone (PIR) PK
GUT_PIR  : Pirfenidone gut depot (µg)
CENT_PIR : Pirfenidone central compartment (µg)
PERI_PIR : Pirfenidone peripheral compartment (µg)
// Nintedanib (NIN) PK
GUT_NIN  : Nintedanib gut depot (ng)
CENT_NIN : Nintedanib central compartment (ng)
PERI_NIN : Nintedanib peripheral compartment (ng)
// Disease PD
AEC2     : AEC-II cell population (normalized)
TGFb     : Active TGF-β1 level (normalized)
M2       : M2 macrophage activation (normalized)
ROS      : Reactive oxygen species (normalized)
FIBRO    : Activated fibroblast pool (normalized)
MYOFIB   : Myofibroblast pool (normalized)
COLLAGEN : ECM collagen accumulation (normalized)
MMP      : Matrix metalloproteinase activity (normalized)
TIMP     : TIMP activity (normalized)
FVC_st   : FVC state variable (% predicted)
DLCO_st  : DLCO state variable (% predicted)

$MAIN
// ── Derived PK concentrations — the single exposed, PD-facing concentration ──
// Computed here (not $ODE) to match where the original placed this exactly:
// mrgsolve evaluates $MAIN once per requested output/dosing interval using
// the start-of-interval state (a one-interval-late lag vs. the true
// continuous state -- see UPSTREAM_ISSUES.md #113 and the identical
// rheumatoid-arthritis precedent, UPSTREAM_ISSUES.md #31 part 2). Moving
// these to $ODE would make them update continuously instead and silently
// change every downstream PD trajectory -- confirmed empirically: doing so
// during development shifted TGFb by up to ~35% at the verification grid.
double C_PIR = CENT_PIR / V1_PIR;   // µg/mL
double C_NIN = (CENT_NIN / V1_NIN) / MW_NIN * 1e3;  // nM

// ── Hill effect interface (also $MAIN, same reason as above) ────────────────
// Pirfenidone acts on three separate downstream pathways from one
// concentration (TGF-β, M2/cytokine, ROS) -- kept as three named
// EFFECT_PIR_* terms rather than collapsed into one, exactly mirroring the
// the three separate inh_P_* terms the original had. The original reused
// Emax_P_TGFb for the M2 term and reused EC50_P_TGFb for the ROS term (no
// separate EMAX_PIR_M2 / EC50_PIR_ROS ever existed) -- that parameter-
// sharing quirk is preserved exactly, not "fixed" into three independent
// parameter pairs.
double EFFECT_PIR_TGFB = EMAX_PIR_TGFB * pow(C_PIR, GAMMA_PIR_TGFB)
                          / (pow(EC50_PIR_TGFB, GAMMA_PIR_TGFB) + pow(C_PIR, GAMMA_PIR_TGFB));
double EFFECT_PIR_M2   = EMAX_PIR_TGFB * pow(C_PIR, GAMMA_PIR_M2)
                          / (pow(EC50_PIR_M2, GAMMA_PIR_M2) + pow(C_PIR, GAMMA_PIR_M2));
double EFFECT_PIR_ROS  = EMAX_PIR_ROS * pow(C_PIR, GAMMA_PIR_ROS)
                          / (pow(EC50_PIR_TGFB, GAMMA_PIR_ROS) + pow(C_PIR, GAMMA_PIR_ROS));
double EFFECT_NIN      = EMAX_NIN * pow(C_NIN, GAMMA_NIN)
                          / (pow(EC50_NIN, GAMMA_NIN) + pow(C_NIN, GAMMA_NIN));

// ── MMP:TIMP ratio → net ECM remodeling ──────────────────────────────────────
double net_ECM = (TIMP > 0) ? MMP / TIMP : 1.0;

// ── Non-zero initial conditions ──────────────────────────────────────────────
// Build-compat fix: the original used a separate $INIT block that jointly
// redeclared every $CMT @annotated compartment -- mrgsolve 2.0.1 rejects
// that as "Duplicated model names" (UPSTREAM_ISSUES.md #113). Moved here
// using the <cmt>_0 idiom; PK compartments keep the implicit 0 default,
// unchanged.
AEC2_0     = 1.0;
TGFb_0     = 1.0;
M2_0       = 1.0;
ROS_0      = 1.0;
FIBRO_0    = 1.0;
MYOFIB_0   = 1.0;
COLLAGEN_0 = 1.0;
MMP_0      = 1.0;
TIMP_0     = 1.0;
FVC_st_0   = 80.0;
DLCO_st_0  = 65.0;

$ODE
// ── Pirfenidone (PIR) PK — archetype 3 (depot + central + peripheral, linear) ─
dxdt_GUT_PIR  = -KA_PIR * GUT_PIR;
dxdt_CENT_PIR =  KA_PIR * F_PIR * GUT_PIR
                 - (CL_PIR + Q_PIR) / V1_PIR * CENT_PIR
                 + Q_PIR / V2_PIR * PERI_PIR;
dxdt_PERI_PIR =  Q_PIR / V1_PIR * CENT_PIR - Q_PIR / V2_PIR * PERI_PIR;

// ── Nintedanib (NIN) PK — archetype 3 (depot + central + peripheral, linear) ──
dxdt_GUT_NIN  = -KA_NIN * GUT_NIN;
dxdt_CENT_NIN =  KA_NIN * F_NIN * GUT_NIN
                 - (CL_NIN + Q_NIN) / V1_NIN * CENT_NIN
                 + Q_NIN / V2_NIN * PERI_NIN;
dxdt_PERI_NIN =  Q_NIN / V1_NIN * CENT_NIN - Q_NIN / V2_NIN * PERI_NIN;

// ── AEC-II dynamics ───────────────────────────────────────────────────────────
// Loss: ROS-driven damage, senescence by TGF-β feedback
// Gain: repair (EGF/HGF), proportional to remaining cells
double AEC2_damage_rate = k_AEC * ROS * AEC2;
double AEC2_repair_rate = k_rep * AEC2_ss * (1.0 - AEC2);
dxdt_AEC2 = AEC2_repair_rate - AEC2_damage_rate;

// ── TGF-β1 dynamics ──────────────────────────────────────────────────────────
// Production: M2 macrophages, damaged AEC, myofibroblasts
// Inhibition: pirfenidone
double TGFb_prod = kprod_TGFb
                   + kact_TGFb * (AEC2_ss - AEC2)  // more damage → more TGF-β
                   + 0.06 * M2
                   + 0.04 * MYOFIB;
double TGFb_deg  = kdeg_TGFb * TGFb;
dxdt_TGFb = TGFb_prod * (1.0 - EFFECT_PIR_TGFB) - TGFb_deg;

// ── M2 macrophage dynamics ────────────────────────────────────────────────────
double M2_prod = kM2_act * TGFb;  // TGF-β promotes M2 polarization
dxdt_M2 = M2_prod * (1.0 - EFFECT_PIR_M2) - kdeg_M2 * M2;

// ── ROS dynamics ─────────────────────────────────────────────────────────────
double ROS_prod = kprod_ROS
                  + kfb_ROS * (AEC2_ss - AEC2) * TGFb;
double ROS_clear = kdeg_ROS * ROS;
dxdt_ROS = ROS_prod * (1.0 - EFFECT_PIR_ROS) - ROS_clear;

// ── Fibroblast activation ─────────────────────────────────────────────────────
// Activated by TGF-β, PDGF (proxied by nintedanib block)
double F_activ = kact_F * TGFb * F_ss;
double F_apop  = kapop_F * FIBRO;
dxdt_FIBRO = F_activ * (1.0 - EFFECT_NIN) - F_apop;

// ── Myofibroblast differentiation ────────────────────────────────────────────
double M_diff   = kdiff_M * FIBRO * TGFb;
double M_apop   = 0.006 * MYOFIB;  // very low apoptosis in IPF (bcl-2 upregulated)
dxdt_MYOFIB = M_diff * (1.0 - EFFECT_NIN * 0.5) - M_apop;

// ── ECM — Collagen ────────────────────────────────────────────────────────────
double Col_prod = kprod_Col * MYOFIB;
double Col_deg  = kdeg_Col  * MMP / (TIMP + 0.1) * COLLAGEN;
dxdt_COLLAGEN = Col_prod - Col_deg;

// ── MMP dynamics (MMP-1, -7, -9) ─────────────────────────────────────────────
double MMP_prod = kprod_MMP * FIBRO * 0.5 + kprod_MMP * M2 * 0.5;
double MMP_deg  = kdeg_MMP * MMP;
dxdt_MMP = MMP_prod - MMP_deg;

// ── TIMP dynamics (TGF-β upregulates TIMP-1) ─────────────────────────────────
double TIMP_prod = kprod_TIMP * TGFb * MYOFIB;
double TIMP_deg  = kdeg_TIMP * TIMP;
dxdt_TIMP = TIMP_prod - TIMP_deg;

// ── FVC (% predicted) — primary endpoint ─────────────────────────────────────
// Collagen excess above baseline drives FVC decline
double col_excess = (COLLAGEN > 1.0) ? (COLLAGEN - 1.0) : 0.0;
double fvc_loss = k_FVC_loss * col_excess * FVC_st;
dxdt_FVC_st = -fvc_loss;

// ── DLCO (% predicted) ────────────────────────────────────────────────────────
double dlco_loss = k_DLCO_loss * (AEC2_ss - AEC2 + 0.5 * col_excess) * DLCO_st;
dxdt_DLCO_st = -dlco_loss;

$TABLE
// Cp_pirf/Cn_nint intentionally keep the same independent recomputation the
// original had (same formula as C_PIR/C_NIN) rather than being redirected
// to reference $MAIN C_PIR/C_NIN directly. $TABLE runs at the actual
// requested output time (using the current, non-lagged state), while
// $MAIN C_PIR/C_NIN carry the one-interval-late lag described above --
// pointing these at C_PIR/C_NIN would inject that lag into what was
// previously a fresh, unlagged report and silently change the numeric
// values of these two outputs. Not treated as a "duplicate site" to collapse.
double Cp_pirf  = CENT_PIR / V1_PIR;                      // µg/mL
double Cn_nint  = (CENT_NIN / V1_NIN) / MW_NIN * 1e3;      // nM
double FVC_pct  = FVC_st;
double DLCO_pct = DLCO_st;
double FVC_decline_yr = k_FVC_loss * ((COLLAGEN > 1.0) ? (COLLAGEN - 1.0) : 0.0) * FVC_st * 8760.0;
double Col_norm = COLLAGEN;
double AEC2_norm = AEC2;
double TGFb_norm = TGFb;
double MMP_norm  = MMP;
double TIMP_norm = TIMP;
double MMP_TIMP_ratio = (TIMP > 0) ? MMP / TIMP : 1.0;
double Myofib_norm = MYOFIB;
double ROS_norm  = ROS;
// Periostin proxy (correlates with collagen deposition)
double Periostin_proxy = COLLAGEN * 1.2 * MYOFIB;
// MMP-7 serum proxy
double MMP7_proxy = MMP * TGFb * 0.8;
// KL-6 proxy
double KL6_proxy = (AEC2_ss - AEC2 + 0.5) * 800.0;  // U/mL scale

$CAPTURE
C_PIR EFFECT_PIR_TGFB EFFECT_PIR_M2 EFFECT_PIR_ROS C_NIN EFFECT_NIN
Cp_pirf Cn_nint FVC_pct DLCO_pct FVC_decline_yr Col_norm
AEC2_norm TGFb_norm MMP_norm TIMP_norm MMP_TIMP_ratio Myofib_norm ROS_norm
Periostin_proxy MMP7_proxy KL6_proxy
'

## ─── 2. COMPILE MODEL ───────────────────────────────────────────────────────

mod <- mcode("IPF_QSP_refactored", ipf_refactored_code)
cat("Model compiled successfully.\n")
cat("Compartments:", nrow(init(mod)), "\n")

## ─── 3. DOSING REGIMENS ─────────────────────────────────────────────────────

# Pirfenidone 801 mg TID (every 8h) × MW=185.22 g/mol
dose_pirf_mg  <- 801     # mg
dose_pirf_ug  <- dose_pirf_mg * 1e3  # µg

# Nintedanib 150 mg BID (every 12h) × MW=539.63 g/mol
dose_nint_mg  <- 150
dose_nint_ng  <- dose_nint_mg * 1e6  # ng

# Simulation duration (52 weeks = 364 days)
sim_duration <- 52 * 7 * 24  # hours
obs_times    <- seq(0, sim_duration, by=24)  # daily observations

## Treatment scenarios (cmt targets renamed GUT_PIR/GUT_NIN, unchanged doses/timing)
mk_pirf_ev <- function() {
  ev(cmt="GUT_PIR", amt=dose_pirf_ug, ii=8, addl=sim_duration/8, time=0)
}
mk_nint_ev <- function() {
  ev(cmt="GUT_NIN", amt=dose_nint_ng, ii=12, addl=sim_duration/12, time=0)
}
mk_combo_ev <- function() {
  c(ev(cmt="GUT_PIR", amt=dose_pirf_ug, ii=8,  addl=sim_duration/8,  time=0),
    ev(cmt="GUT_NIN", amt=dose_nint_ng,  ii=12, addl=sim_duration/12, time=0))
}

## ─── 4. SIMULATION — 5 TREATMENT SCENARIOS ──────────────────────────────────

cat("\nRunning 5 treatment scenarios...\n")

scenarios <- list(
  "Placebo (Natural History)"         = ev(time=0, amt=0),
  "Pirfenidone 801 mg TID"            = mk_pirf_ev(),
  "Nintedanib 150 mg BID"             = mk_nint_ev(),
  "Combination (Pirf + Nint)"         = mk_combo_ev(),
  "Pirfenidone Low Dose (267 mg TID)" = ev(cmt="GUT_PIR", amt=267e3, ii=8, addl=sim_duration/8)
)

results <- lapply(names(scenarios), function(sc) {
  ev_obj <- scenarios[[sc]]
  out <- mod %>%
    mrgsim(events=ev_obj, end=sim_duration, delta=24, digits=4) %>%
    as.data.frame()
  out$Scenario <- sc
  out
}) %>% bind_rows()

# Convert time to weeks
results$Week <- results$time / (24 * 7)

cat("Simulation complete. Rows:", nrow(results), "\n")

## ─── 5. PK PROFILE (first 72h) ──────────────────────────────────────────────

pk_data <- mod %>%
  mrgsim(
    events = c(ev(cmt="GUT_PIR", amt=dose_pirf_ug, time=0),
               ev(cmt="GUT_NIN", amt=dose_nint_ng, time=0)),
    end = 72, delta = 0.5
  ) %>%
  as.data.frame()

## ─── 6. DOSE-RESPONSE ANALYSIS ──────────────────────────────────────────────

cat("\nDose-response analysis...\n")

pirf_doses <- c(267, 534, 801, 1068) * 1e3  # µg
nint_doses  <- c(50, 100, 150, 200) * 1e6   # ng

dr_pirf <- lapply(pirf_doses, function(d) {
  ev_d <- ev(cmt="GUT_PIR", amt=d, ii=8, addl=sim_duration/8)
  out  <- mod %>%
    mrgsim(events=ev_d, end=sim_duration, delta=24) %>%
    as.data.frame() %>%
    filter(time == max(time))
  data.frame(Drug="Pirfenidone", Dose_mg=d/1e3,
             FVC_final=out$FVC_pct, DLCO_final=out$DLCO_pct,
             Col_final=out$Col_norm)
}) %>% bind_rows()

dr_nint <- lapply(nint_doses, function(d) {
  ev_d <- ev(cmt="GUT_NIN", amt=d, ii=12, addl=sim_duration/12)
  out  <- mod %>%
    mrgsim(events=ev_d, end=sim_duration, delta=24) %>%
    as.data.frame() %>%
    filter(time == max(time))
  data.frame(Drug="Nintedanib", Dose_mg=d/1e6,
             FVC_final=out$FVC_pct, DLCO_final=out$DLCO_pct,
             Col_final=out$Col_norm)
}) %>% bind_rows()

dose_response <- bind_rows(dr_pirf, dr_nint)

## ─── 7. BIOMARKER TRAJECTORIES ──────────────────────────────────────────────

bm_data <- results %>%
  select(Week, Scenario, TGFb_norm, Myofib_norm, Col_norm,
         Periostin_proxy, MMP7_proxy, KL6_proxy, ROS_norm,
         MMP_TIMP_ratio, AEC2_norm)

## ─── 8. PLOT FUNCTIONS ──────────────────────────────────────────────────────

scenario_colors <- c(
  "Placebo (Natural History)"         = "#E74C3C",
  "Pirfenidone 801 mg TID"            = "#2471A3",
  "Nintedanib 150 mg BID"             = "#27AE60",
  "Combination (Pirf + Nint)"         = "#8E44AD",
  "Pirfenidone Low Dose (267 mg TID)" = "#E67E22"
)

## FVC over time
p_fvc <- ggplot(results, aes(x=Week, y=FVC_pct, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="FVC % Predicted Over 52 Weeks (IPF QSP Model)",
       x="Time (weeks)", y="FVC (% predicted)",
       caption="ASCEND: pirfenidone −47.9% decline reduction; INPULSIS: nintedanib −50.1%") +
  theme_bw(base_size=13) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=c(50, 70, 80), linetype="dashed", color="gray60", alpha=0.5)

## DLCO over time
p_dlco <- ggplot(results, aes(x=Week, y=DLCO_pct, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="DLCO % Predicted Over 52 Weeks",
       x="Time (weeks)", y="DLCO (% predicted)") +
  theme_bw(base_size=13) +
  theme(legend.position="bottom", legend.title=element_blank())

## PK profiles
p_pk_pirf <- ggplot(pk_data, aes(x=time, y=Cp_pirf)) +
  geom_line(color="#2471A3", linewidth=1.5) +
  labs(title="Pirfenidone PK — Single Dose (801 mg)",
       x="Time (h)", y="Plasma Concentration (µg/mL)") +
  geom_hline(yintercept=30, linetype="dashed", color="red", label="EC50") +
  annotate("text", x=60, y=32, label="EC50 ~30 µg/mL", color="red", size=3.5) +
  theme_bw(base_size=13)

p_pk_nint <- ggplot(pk_data, aes(x=time, y=Cn_nint)) +
  geom_line(color="#27AE60", linewidth=1.5) +
  labs(title="Nintedanib PK — Single Dose (150 mg)",
       x="Time (h)", y="Plasma Concentration (nM)") +
  geom_hline(yintercept=20, linetype="dashed", color="red") +
  annotate("text", x=60, y=22, label="IC50 ~20 nM (FGFR1)", color="red", size=3.5) +
  theme_bw(base_size=13)

## Biomarker dynamics
p_tgfb <- ggplot(results, aes(x=Week, y=TGFb_norm, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="TGF-β1 Level (Normalized)", x="Time (weeks)", y="TGF-β1 (normalized)") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_col <- ggplot(results, aes(x=Week, y=Col_norm, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Collagen Deposition (Normalized)", x="Time (weeks)", y="Collagen") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_mmp7 <- ggplot(results, aes(x=Week, y=MMP7_proxy, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Serum MMP-7 Proxy", x="Time (weeks)", y="MMP-7 (proxy)") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_kl6 <- ggplot(results, aes(x=Week, y=KL6_proxy, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="KL-6 Biomarker Proxy (U/mL)", x="Time (weeks)", y="KL-6 (U/mL equiv.)") +
  theme_bw(base_size=13) + theme(legend.position="none")

## Dose-response
p_dr <- ggplot(dose_response, aes(x=Dose_mg, y=FVC_final, color=Drug)) +
  geom_line(linewidth=1.4) + geom_point(size=3) +
  scale_color_manual(values=c("Pirfenidone"="#2471A3", "Nintedanib"="#27AE60")) +
  labs(title="Dose-Response: Final FVC at 52 Weeks",
       x="Dose (mg)", y="FVC % predicted at 52 weeks") +
  theme_bw(base_size=13)

## ─── 9. SUMMARY STATISTICS TABLE ────────────────────────────────────────────

fvc_summary <- results %>%
  filter(Week %in% c(0, 13, 26, 39, 52)) %>%
  group_by(Scenario, Week) %>%
  summarise(FVC=round(mean(FVC_pct), 2),
            DLCO=round(mean(DLCO_pct), 2),
            TGFb=round(mean(TGFb_norm), 3),
            Col=round(mean(Col_norm), 3),
            .groups="drop")

cat("\n=== FVC Summary at Key Timepoints ===\n")
print(as.data.frame(fvc_summary))

## ─── 10. CLINICAL TRIAL CALIBRATION CHECK ───────────────────────────────────

placebo_52  <- results %>% filter(Scenario=="Placebo (Natural History)", Week==52)
pirf_52     <- results %>% filter(Scenario=="Pirfenidone 801 mg TID", Week==52)
nint_52     <- results %>% filter(Scenario=="Nintedanib 150 mg BID", Week==52)
combo_52    <- results %>% filter(Scenario=="Combination (Pirf + Nint)", Week==52)

fvc_0       <- 80.0
pirf_0      <- results %>% filter(Scenario=="Pirfenidone 801 mg TID", Week==0)
placebo_decl <- fvc_0 - mean(placebo_52$FVC_pct)
pirf_decl    <- fvc_0 - mean(pirf_52$FVC_pct)
nint_decl    <- fvc_0 - mean(nint_52$FVC_pct)
combo_decl   <- fvc_0 - mean(combo_52$FVC_pct)

cat("\n=== Clinical Calibration (52-week FVC decline in % predicted) ===\n")
cat(sprintf("  Placebo:     −%.2f%%  (Literature target: ~2.5-3%% predicted/yr)\n", placebo_decl))
cat(sprintf("  Pirfenidone: −%.2f%%  (Target: ~47-50%% reduction vs placebo)\n", pirf_decl))
cat(sprintf("  Nintedanib:  −%.2f%%  (Target: ~50%% reduction vs placebo)\n", nint_decl))
cat(sprintf("  Combination: −%.2f%%  (Target: ≥50%% reduction vs placebo)\n", combo_decl))
if (placebo_decl > 0) {
  cat(sprintf("  Pirf reduction vs placebo: %.1f%%\n", (1 - pirf_decl/placebo_decl)*100))
  cat(sprintf("  Nint reduction vs placebo: %.1f%%\n", (1 - nint_decl/placebo_decl)*100))
}

## NOTE (added by the PK/PD refactor): the calibration check above (section
## 10) re-uses the original file's own 52-week (8736h) horizon exactly as
## written. Independently of this refactor, that horizon does not actually
## complete under mrgsolve 2.0.1 -- see ipf_refactor_notes.md and
## UPSTREAM_ISSUES.md #114: the disease system's own TGF-β<->M2 and
## TGF-β<->myofibroblast positive-feedback loops (no saturation term
## anywhere) diverge and make the ODE solver fail around simulated hour
## 112-114 in every scenario (identically in the untouched original), long
## before section 10's Week==52 filters would find any matching rows. This
## is a pre-existing disease-side defect, reproduced unchanged here, not
## introduced or fixed by the PK renaming above. Verification for this
## refactor used a shortened 0-96h window instead (see refactor notes).

## ─── 11. SAVE PLOTS ──────────────────────────────────────────────────────────

cat("\nPlots generated (use gridExtra or patchwork to display):\n")
cat("  p_fvc, p_dlco, p_pk_pirf, p_pk_nint\n")
cat("  p_tgfb, p_col, p_mmp7, p_kl6, p_dr\n")

## Combined plot (4×2 grid)
if (requireNamespace("gridExtra", quietly=TRUE)) {
  grid_plot <- gridExtra::grid.arrange(
    p_fvc, p_dlco, p_pk_pirf, p_pk_nint,
    p_tgfb, p_col, p_dr, p_kl6,
    ncol=2, nrow=4
  )
}

cat("\nIPF QSP model (refactored) simulation complete.\n")
cat("Results in 'results' data frame, summaries in 'fvc_summary'.\n")
