## =============================================================================
## COPD QSP Model — mrgsolve ODE Implementation (PK/PD REFACTORED)
## =============================================================================
## Disease: Chronic Obstructive Pulmonary Disease (COPD)
##
## This is a PK/PD refactor of copd_mrgsolve_model.R — pluggable-PK naming
## convention applied to all four modeled compounds (LAMA/LABA/ICS/PDE4I).
## See copd_refactor_notes.md for full details: archetype per compound,
## the pre-existing build defect found and fixed here (disclosed), and the
## verification result. The original file is untouched.
##
## Compartments (19 ODE states):
##   PK  — LAMA (Tiotropium), LABA (Salmeterol), ICS (Budesonide), PDE4I (Roflumilast)
##   PD  — Inflammation mediators (IL-8, NE, CRP, Eosinophils)
##         FEV1 dynamics (natural decline + drug effect + exacerbation)
##         Emphysema / ECM index
##         Exacerbation risk (Poisson rate)
##         RV function / PVR
##
## Treatment scenarios:
##   1. Placebo / Natural history
##   2. LAMA monotherapy (Tiotropium 18 µg/day)
##   3. LABA/LAMA dual therapy (Salmeterol 50 µg bid + Tiotropium 18 µg/day)
##   4. ICS/LABA (Budesonide 400 µg bid + Salmeterol 50 µg bid)
##   5. Triple therapy LAMA/LABA/ICS
##   6. Triple + PDE4I (Roflumilast 500 µg/day; severe COPD + chronic bronchitis)
##
## Key clinical trials calibrated:
##   UPLIFT (Tiotropium), TORCH (ICS/LABA), POET-COPD (LAMA vs LABA),
##   SUMMIT (FP/VI), IMPACT (FF/UMEC/VI triple), ETHOS (BUD/GLYCO/FORM),
##   FLAME (Indacaterol/Glyco), EINSTEIN (Roflumilast)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# =============================================================================
# mrgsolve model specification (refactored: pluggable PK naming convention)
# =============================================================================
code_copd_refactored <- '
$PROB
COPD QSP ODE Model — LAMA/LABA/ICS/PDE4I PK + Disease PD (refactored)

$PARAM @annotated
// ── LAMA PK (Tiotropium) ──────────────────────────────────────────────────
KA_LAMA      : 0.5      : LAMA lung absorption rate constant (1/h) — inhaled
F_LUNG_LAMA  : 0.20     : Lung fraction deposited
F_LAMA       : 0.03     : Systemic bioavailability fraction (lung depot -> central)
CL_LAMA      : 3.92     : Tiotropium systemic CL (L/h) — UPLIFT PK
V1_LAMA      : 36.0     : Central volume (L) tiotropium
Q_LAMA       : 0.5      : Inter-compartmental CL (L/h)
V2_LAMA      : 8.0      : Peripheral volume (L)
EC50_LAMA    : 0.10     : EC50 for M3 bronchodilation (ng/mL lung)
EMAX_LAMA    : 0.18     : Max FEV1 gain from LAMA (fraction of pred)
GAMMA_LAMA   : 1        : Hill coefficient (no explicit Hill term in original; implicit 1)

// ── LABA PK (Salmeterol) ─────────────────────────────────────────────────
KA_LABA      : 1.2      : LABA lung absorption rate (1/h)
F_LUNG_LABA  : 0.20     : Lung fraction deposited
F_LABA       : 0.05     : Systemic fraction absorbed (lung depot -> central)
CL_LABA      : 23.0     : Salmeterol systemic CL (L/h) — Pope 1999
V1_LABA      : 45.0     : Central volume (L)
Q_LABA       : 2.0      : Intercompartmental CL (L/h)
V2_LABA      : 12.0     : Peripheral volume (L)
EC50_LABA    : 0.08     : EC50 for β2 bronchodilation (ng/mL lung)
EMAX_LABA    : 0.15     : Max FEV1 gain from LABA (fraction of pred)
GAMMA_LABA   : 1        : Hill coefficient (no explicit Hill term in original; implicit 1)

// ── ICS PK (Budesonide) ──────────────────────────────────────────────────
KA_ICS       : 0.60     : ICS lung absorption (1/h)
F_LUNG_ICS   : 0.22     : Lung fraction deposited
F_ICS        : 0.10     : Systemic fraction absorbed (gut first-pass 90%)
CL_ICS       : 78.0     : Budesonide systemic CL (L/h) — Thorsson 1994
V1_ICS       : 180.0    : Central volume (L)
Q_ICS        : 8.0      : Intercompartmental CL (L/h)
V2_ICS       : 60.0     : Peripheral volume (L)
EC50_ICS     : 0.15     : EC50 for anti-inflammatory GR effect (ng/mL lung)
EMAX_ICS     : 0.30     : Max fractional reduction in IL8 / exacerbation rate
GAMMA_ICS    : 1.5      : Hill coefficient for ICS PD
Eos_thresh   : 300.0    : Blood eosinophil threshold for ICS response (cells/µL)

// ── PDE4I PK (Roflumilast) ───────────────────────────────────────────────
KA_PDE4I     : 0.70     : Oral absorption rate (1/h) — declared, unused (see refactor notes)
F_PDE4I      : 0.80     : Oral bioavailability — declared, unused (see refactor notes)
CL_PDE4I     : 9.1      : Apparent oral CL (L/h) — Bethke 2007
V1_PDE4I     : 210.0    : Central volume (L)
EMAX_PDE4I   : 0.20     : Max reduction in exacerbation rate
EC50_PDE4I   : 2.5      : EC50 ng/mL for cAMP-mediated anti-inflammation
GAMMA_PDE4I  : 1        : Hill coefficient (no explicit Hill term in original; implicit 1)

// ── Inflammation PD ───────────────────────────────────────────────────────
ksyn_IL8    : 15.0     : IL-8 baseline production (pg/mL/h)
kout_IL8    : 0.12     : IL-8 elimination rate (1/h)   → t½ ~5.8h
IL80        : 125.0    : Baseline sputum IL-8 (pg/mL) — Donnelly 1993 (renamed from IL8_0, see $MAIN note)
ksyn_NE     : 4.0      : NE baseline production (µg/mL/h)
kout_NE     : 0.08     : NE elimination (1/h)           → t½ ~8.7h
NE_0        : 50.0     : Baseline sputum NE (µg/mL)
ksyn_CRP    : 0.10     : CRP synthesis rate (mg/L/h)
kout_CRP    : 0.020    : CRP elimination (1/h)           → t½ ~35h
CRP0        : 5.0      : Baseline serum CRP (mg/L) (renamed from CRP_0, see $MAIN note; unused elsewhere in original)
ksyn_Eos    : 2.0      : Eosinophil turnover rate (cells/µL/h)
kout_Eos    : 0.003    : Eosinophil removal rate (1/h)
Eos0        : 200.0    : Baseline blood eosinophils (cells/µL) (renamed from Eos_0, see $MAIN note; unused elsewhere in original DSL, overridden from the R-side dose_response scenario)

// ── FEV1 Disease Dynamics ─────────────────────────────────────────────────
FEV10       : 55.0     : Baseline FEV1 % predicted (GOLD 2-3 phenotype) (renamed from FEV1_0, see $MAIN note)
k_FEV1_dec  : 0.000055 : Natural FEV1 decline rate (fraction/h) → ~48 mL/yr
k_AE_FEV1   : 0.015    : FEV1 decline per moderate/severe exacerbation
k_inflam_FEV1 : 0.0002 : Inflammation-driven FEV1 decline amplifier

// ── Emphysema / Structural ────────────────────────────────────────────────
Emph0       : 20.0     : Baseline emphysema index (% LAA on CT) (renamed from Emph_0, see $MAIN note; unused elsewhere in original)
k_emph_prog : 0.000020 : Emphysema progression rate (fraction/h) [irreversible]

// ── Exacerbation Rate ─────────────────────────────────────────────────────
lambda_AE_0 : 1.8      : Baseline annual exacerbation rate (events/yr)
k_AE        : 0.000205 : Rate constant for exacerbation accumulation (events/h)
k_post_AE   : 0.50     : Post-exacerbation increase in rate (multiplicative)

// ── Pulmonary Vascular ────────────────────────────────────────────────────
PVR0        : 280.0    : Baseline PVR (dyn·s·cm⁻⁵) — elevated in GOLD 3-4 (renamed from PVR_0, see $MAIN note)
mPAP_0      : 26.0     : Baseline mPAP (mmHg)
k_PVR_prog  : 0.000015 : PVR progression rate with FEV1↓
CO_0        : 4.5      : Cardiac output at rest (L/min)
PAWP_0      : 8.0      : Pulmonary artery wedge pressure (mmHg)

// ── Dose Flags ───────────────────────────────────────────────────────────
DOSE_LAMA   : 0        : LAMA dosing on (1) or off (0)
DOSE_LABA   : 0        : LABA dosing on (1) or off (0)
DOSE_ICS    : 0        : ICS dosing on (1) or off (0)
DOSE_PDE4I  : 0        : PDE4I dosing on (1) or off (0)

$CMT @annotated
// LAMA PK
GUT_LAMA    : LAMA lung depot (ng)
CENT_LAMA   : LAMA central plasma (ng/mL * V1)
PERI_LAMA   : LAMA peripheral

// LABA PK
GUT_LABA    : LABA lung depot (ng)
CENT_LABA   : LABA central plasma
PERI_LABA   : LABA peripheral

// ICS PK
GUT_ICS     : ICS lung depot (ng)
CENT_ICS    : ICS central plasma
PERI_ICS    : ICS peripheral

// PDE4I PK
CENT_PDE4I  : PDE4I central (ng/mL * V1)

// PD — Inflammation
IL8         : Sputum IL-8 (pg/mL)
NE_sput     : Sputum NE (µg/mL)
CRP         : Serum CRP (mg/L)
Eos         : Blood eosinophils (cells/µL)

// Disease structural states
FEV1        : FEV1 % predicted
Emph        : Emphysema index (% LAA)
PVR         : Pulmonary vascular resistance (dyn·s·cm⁻⁵)

// Cumulative exacerbation counter
AE_cum      : Cumulative exacerbations
AE_rate_ann : Annualized exacerbation rate

$MAIN
// Initial values (all PK compartments default to mrgsolve’s implicit 0 and
// are not restated here). Uses the modern `<CMT>_0 = value;` idiom in $MAIN
// rather than a separate $INIT block -- see the pre-existing build defects
// disclosed in the refactor notes: (1) the original’’’s own $INIT block (bare
// `NAME = value` assignments) redeclares every $CMT compartment a second
// time under mrgsolve 2.0.1, which treats $INIT as an alternative
// compartment-declaring block, not a companion to $CMT (same defect class as
// UPSTREAM_ISSUES.md #38/#42/#51); (2) six of the original’’’s own baseline
// $PARAM names (IL8_0/CRP_0/Eos_0/FEV1_0/Emph_0/PVR_0) collide with
// mrgsolve’’’s own auto-reserved per-compartment `<CMT>_0` initial-value
// symbol for the identically-named compartments IL8/CRP/Eos/FEV1/Emph/PVR
// (same defect class as UPSTREAM_ISSUES.md #60) -- those six params are
// renamed IL80/CRP0/Eos0/FEV10/Emph0/PVR0 below (dropping the underscore,
// same convention as #60’’’s fix) to free up the `_0` symbol for its
// mrgsolve-reserved initial-value use. Both fixes are syntax-only and
// non-numeric: every initial value below is identical to the original’’’s
// own $INIT values, and every renamed param keeps its original value.
IL8_0         = 125;
NE_sput_0     = 50;
CRP_0         = 5;
Eos_0         = 200;
FEV1_0        = 55;
Emph_0        = 20;
PVR_0         = 280;
AE_rate_ann_0 = 1.8;

$GLOBAL
// C_<STEM>/EFFECT_<STEM> live here, not in $PARAM: mrgsolve 2.0.1 compiles
// $PARAM members as read-only references inside $ODE, so a value that must
// be recomputed every timestep from state cannot also be a $PARAM (confirmed
// against this exact file with a minimal reproduction: a $PARAM `C_X = 0.0`
// reassigned in $ODE fails with "assignment of read-only reference C_X" —
// same constraint already documented in several sibling refactors, e.g.
// breast-cancer/bc_refactor_notes.md). Declared here instead and exposed via
// $CAPTURE, so every covariate is still discoverable in simulation output.
double C_LAMA;
double C_LAMA_PLASMA;
double EFFECT_LAMA;
double C_LABA;
double C_LABA_PLASMA;
double EFFECT_LABA;
double C_ICS;
double C_ICS_PLASMA;
double EFFECT_ICS;
double C_PDE4I;
double EFFECT_PDE4I;

$ODE
// ── LAMA PK (Archetype 3: depot + central + peripheral, linear) ───────────
// C_LAMA is the lung (epithelial-lining-fluid) concentration -- the ONLY
// concentration PD reads for this compound (see refactor notes: genuinely
// distinct tissue site from systemic plasma, same allowance the guide gives
// TMDD/tissue-compartment archetypes). It is derived directly from the lung
// depot state, not from CENT_LAMA/V1_LAMA -- no new compartment was needed
// since the depot already represents the lung deposit.
C_LAMA = GUT_LAMA / (F_LUNG_LAMA * 200.0); // ng/mL lung ELF
dxdt_GUT_LAMA  = -KA_LAMA * GUT_LAMA;
dxdt_CENT_LAMA =  KA_LAMA * GUT_LAMA * F_LAMA
                  - (CL_LAMA + Q_LAMA) * (CENT_LAMA / V1_LAMA)
                  + Q_LAMA * (PERI_LAMA / V2_LAMA);
dxdt_PERI_LAMA =  Q_LAMA * (CENT_LAMA / V1_LAMA)
                  - Q_LAMA * (PERI_LAMA / V2_LAMA);
C_LAMA_PLASMA = CENT_LAMA / V1_LAMA; // systemic plasma -- informational only, never read by PD

// ── LABA PK (Archetype 3) ──────────────────────────────────────────────────
C_LABA = GUT_LABA / (F_LUNG_LABA * 200.0);
dxdt_GUT_LABA  = -KA_LABA * GUT_LABA;
dxdt_CENT_LABA =  KA_LABA * GUT_LABA * F_LABA
                  - (CL_LABA + Q_LABA) * (CENT_LABA / V1_LABA)
                  + Q_LABA * (PERI_LABA / V2_LABA);
dxdt_PERI_LABA =  Q_LABA * (CENT_LABA / V1_LABA)
                  - Q_LABA * (PERI_LABA / V2_LABA);
C_LABA_PLASMA = CENT_LABA / V1_LABA;

// ── ICS PK (Archetype 3) ───────────────────────────────────────────────────
C_ICS = GUT_ICS / (F_LUNG_ICS * 200.0);
dxdt_GUT_ICS  = -KA_ICS * GUT_ICS;
dxdt_CENT_ICS =  KA_ICS * GUT_ICS * F_ICS
                 - (CL_ICS + Q_ICS) * (CENT_ICS / V1_ICS)
                 + Q_ICS * (PERI_ICS / V2_ICS);
dxdt_PERI_ICS =  Q_ICS * (CENT_ICS / V1_ICS)
                 - Q_ICS * (PERI_ICS / V2_ICS);
C_ICS_PLASMA = CENT_ICS / V1_ICS;

// ── PDE4I PK (Archetype 1: no depot, single compartment, linear) ──────────
// KA_PDE4I/F_PDE4I (declared above) are not used here: dosing goes directly
// into CENT_PDE4I as a bolus (see the original’s own dose_PDE4i event and
// refactor notes) -- these are dead parameters carried over unchanged.
// PDE4I is oral, so its exposed concentration is the ordinary central/
// plasma concentration (no separate lung tissue site, unlike the three
// inhaled compounds above).
C_PDE4I = CENT_PDE4I / V1_PDE4I;
dxdt_CENT_PDE4I = -CL_PDE4I * C_PDE4I;

// ── Drug Effect Functions (named Hill interface) ──────────────────────────
// LAMA bronchodilation (Emax model on lung conc.)
EFFECT_LAMA = DOSE_LAMA * EMAX_LAMA * pow(C_LAMA, GAMMA_LAMA) /
              (pow(EC50_LAMA, GAMMA_LAMA) + pow(C_LAMA, GAMMA_LAMA));

// LABA bronchodilation
EFFECT_LABA = DOSE_LABA * EMAX_LABA * pow(C_LABA, GAMMA_LABA) /
              (pow(EC50_LABA, GAMMA_LABA) + pow(C_LABA, GAMMA_LABA));

// ICS anti-inflammatory — Eos-guided response (Hill equation)
double Eos_factor = (Eos > Eos_thresh) ? 1.0 : (Eos / Eos_thresh);
EFFECT_ICS = DOSE_ICS * EMAX_ICS *
             pow(C_ICS, GAMMA_ICS) /
             (pow(EC50_ICS, GAMMA_ICS) + pow(C_ICS, GAMMA_ICS)) * Eos_factor;

// PDE4I anti-inflammatory (cAMP-mediated NF-kB suppression)
EFFECT_PDE4I = DOSE_PDE4I * EMAX_PDE4I * pow(C_PDE4I, GAMMA_PDE4I) /
               (pow(EC50_PDE4I, GAMMA_PDE4I) + pow(C_PDE4I, GAMMA_PDE4I));

// ── Inflammation ODEs (disease-side; EFFECT_* combined only here) ────────
// IL-8: CS/TLR → NF-κB driven; suppressed by ICS + PDE4I
double IL8_stim_factor = (FEV1 < 50.0) ? 1.5 : 1.0; // severe → higher IL-8
dxdt_IL8    = ksyn_IL8 * IL8_stim_factor * (1.0 - EFFECT_ICS * 0.5 - EFFECT_PDE4I * 0.4)
              - kout_IL8 * IL8;

// NE: neutrophil-driven; suppressed by ICS
dxdt_NE_sput = ksyn_NE * (IL8 / IL80) * (1.0 - EFFECT_ICS * 0.3)
               - kout_NE * NE_sput;

// CRP: systemic, driven by IL-6 proxy (correlated with IL-8)
dxdt_CRP    = ksyn_CRP * (IL8 / IL80) * (1.0 - EFFECT_ICS * 0.4)
              - kout_CRP * CRP;

// Eosinophils: Th2/ILC2 driven; ICS has a SECOND, non-Hill, flag-only
// effect here (not part of EFFECT_ICS’s concentration-driven interface --
// depends only on the DOSE_ICS switch, not on C_ICS -- kept exactly as the
// original; see refactor notes).
double Eos_stim = (DOSE_ICS == 1.0) ? 0.3 : 1.0; // ICS depletes Eos
dxdt_Eos    = ksyn_Eos * Eos_stim - kout_Eos * Eos;

// ── FEV1 Dynamics ─────────────────────────────────────────────────────────
// Natural decline: 30-50 mL/yr in COPD (converted to fraction/h)
// Drug treatment provides bronchodilatory gain (not on decline rate directly)
// Exacerbations cause step-decrements
// AE rate (events/h → converted)
double AE_rate_h = (lambda_AE_0 * (1.0 - EFFECT_ICS*0.35 - EFFECT_PDE4I*0.15))
                   / (24.0 * 365.0);
// FEV1 bronchodilatory gain (plateau effect when both LAMA+LABA used)
// Combined bronchodilatory gain capped at sum (slightly sub-additive)
double E_BD = EFFECT_LAMA + EFFECT_LABA - 0.3 * EFFECT_LAMA * EFFECT_LABA;
// FEV1 target with treatment
double FEV1_target = FEV10 * (1.0 + E_BD);

// dFEV1/dt: slow trend toward target (bronchodilation) minus disease decline
double inflam_penalty = k_inflam_FEV1 * (IL8 / IL80 - 1.0) * FEV1;
dxdt_FEV1   = 0.05 * (FEV1_target - FEV1)          // approach BD plateau (slow)
              - k_FEV1_dec * FEV1                   // natural decline
              - inflam_penalty                      // inflammation driven decline
              - k_AE_FEV1 * AE_rate_h * FEV1;      // exacerbation step-loss

// ── Emphysema (irreversible progression) ──────────────────────────────────
dxdt_Emph  = k_emph_prog * (NE_sput / NE_0) *
             (1.0 - EFFECT_ICS * 0.05) *       // ICS minimal effect on emphysema
             (100.0 - Emph);              // asymptote at 100%

// ── Pulmonary Vascular Resistance ────────────────────────────────────────
// PVR increases as FEV1 declines (hypoxia-driven)
double FEV1_frac = (FEV1 > 0.1 ? FEV1 : 0.1) / FEV10;
dxdt_PVR   = k_PVR_prog * (1.0 / FEV1_frac - 1.0) * PVR0
             - k_PVR_prog * 0.1 * (PVR - PVR0);

// ── Exacerbation Accumulation ─────────────────────────────────────────────
dxdt_AE_cum      = AE_rate_h;
dxdt_AE_rate_ann = 0.001 * (lambda_AE_0 * (1.0 - EFFECT_ICS*0.35 - EFFECT_PDE4I*0.15)
                             * (1.0 + k_post_AE * AE_cum/100.0)
                             - AE_rate_ann);

$TABLE
double mPAP     = CO_0 * PVR / 80.0 + PAWP_0;
double DAS_FEV1 = FEV1;                        // % predicted
double GOLD_stage = (FEV1 >= 80) ? 1 :
                    (FEV1 >= 50) ? 2 :
                    (FEV1 >= 30) ? 3 : 4;
double CAT_approx = 10.0 + (55.0 - FEV1) * 0.35 + (CRP / 5.0) * 1.2;
if(CAT_approx > 40) CAT_approx = 40.0;
if(CAT_approx < 0)  CAT_approx = 0.0;

// Oxygen saturation proxy (Severinghaus-inspired rough model)
double SpO2 = 98.0 - 0.12 * (PVR / 280.0 - 1.0) * 5.0
              - 0.15 * (100.0 - FEV1) / 50.0;
if(SpO2 > 99.0) SpO2 = 99.0;
if(SpO2 < 85.0) SpO2 = 85.0;

$CAPTURE
C_LAMA C_LAMA_PLASMA C_LABA C_LABA_PLASMA C_ICS C_ICS_PLASMA C_PDE4I
EFFECT_LAMA EFFECT_LABA EFFECT_ICS EFFECT_PDE4I E_BD
mPAP GOLD_stage CAT_approx SpO2
'

mod_copd_refactored <- mcode("COPD_QSP_refactored", code_copd_refactored)

# =============================================================================
# DOSING REGIMENS (identical amounts/timing to the original; compartment
# names updated to the refactored convention)
# =============================================================================
# LAMA: Tiotropium 18 µg once daily
# — lung dose = 18 µg × F_lung(0.20) = 3.6 µg = 3600 ng into GUT_LAMA
dose_LAMA  <- ev(cmt = "GUT_LAMA", amt = 3600, ii = 24, addl = 363, time = 0)

# LABA: Salmeterol 50 µg twice daily
# — lung dose = 50 µg × F_lung(0.20) = 10 µg = 10000 ng bid
dose_LABA  <- ev(cmt = "GUT_LABA", amt = 10000, ii = 12, addl = 727, time = 0)

# ICS: Budesonide 400 µg twice daily
# — lung dose = 400 µg × F_lung(0.22) = 88 µg = 88000 ng bid
dose_ICS   <- ev(cmt = "GUT_ICS", amt = 88000, ii = 12, addl = 727, time = 0)

# PDE4I: Roflumilast 500 µg once daily oral
# — systemic dose = 500 µg × F(0.80) = 400 µg = 400000 ng → CENT_PDE4I
dose_PDE4i <- ev(cmt = "CENT_PDE4I", amt = 400000, ii = 24, addl = 363, time = 0)

# Simulation time: 1 year = 8760 hours
sim_time <- seq(0, 8760, by = 12)

# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================
scenarios <- list(
  "1. Placebo (Natural History)" = list(
    params = list(DOSE_LAMA=0, DOSE_LABA=0, DOSE_ICS=0, DOSE_PDE4I=0),
    ev = NULL
  ),
  "2. LAMA Monotherapy\n(Tiotropium 18µg qd — UPLIFT)" = list(
    params = list(DOSE_LAMA=1, DOSE_LABA=0, DOSE_ICS=0, DOSE_PDE4I=0),
    ev = dose_LAMA
  ),
  "3. LABA/LAMA Dual\n(Salmeterol+Tiotropium — FLAME/POET-COPD)" = list(
    params = list(DOSE_LAMA=1, DOSE_LABA=1, DOSE_ICS=0, DOSE_PDE4I=0),
    ev = c(dose_LAMA, dose_LABA)
  ),
  "4. ICS/LABA Dual\n(Budesonide+Salmeterol — TORCH)" = list(
    params = list(DOSE_LAMA=0, DOSE_LABA=1, DOSE_ICS=1, DOSE_PDE4I=0),
    ev = c(dose_LABA, dose_ICS)
  ),
  "5. Triple Therapy\n(LAMA/LABA/ICS — IMPACT/ETHOS)" = list(
    params = list(DOSE_LAMA=1, DOSE_LABA=1, DOSE_ICS=1, DOSE_PDE4I=0),
    ev = c(dose_LAMA, dose_LABA, dose_ICS)
  ),
  "6. Triple + Roflumilast\n(Severe COPD + Chronic Bronchitis)" = list(
    params = list(DOSE_LAMA=1, DOSE_LABA=1, DOSE_ICS=1, DOSE_PDE4I=1),
    ev = c(dose_LAMA, dose_LABA, dose_ICS, dose_PDE4i)
  )
)

# =============================================================================
# RUN SIMULATIONS
# =============================================================================
run_scenario <- function(name, params_list, evs) {
  mod_tmp <- param(mod_copd_refactored, params_list)
  if (is.null(evs)) {
    out <- mrgsim(mod_tmp, tgrid = sim_time, delta = 12)
  } else {
    out <- mrgsim(mod_tmp, events = evs, tgrid = sim_time, delta = 12)
  }
  as.data.frame(out) %>%
    mutate(scenario = name, time_day = time / 24, time_wk = time / 168)
}

results_list <- mapply(
  FUN      = run_scenario,
  name     = names(scenarios),
  params_list = lapply(scenarios, `[[`, "params"),
  evs      = lapply(scenarios, `[[`, "ev"),
  SIMPLIFY = FALSE
)
results <- bind_rows(results_list)

# =============================================================================
# PLOTS
# =============================================================================
scen_colors <- c(
  "1. Placebo (Natural History)"                           = "#616161",
  "2. LAMA Monotherapy\n(Tiotropium 18µg qd — UPLIFT)"    = "#1565C0",
  "3. LABA/LAMA Dual\n(Salmeterol+Tiotropium — FLAME/POET-COPD)" = "#2E7D32",
  "4. ICS/LABA Dual\n(Budesonide+Salmeterol — TORCH)"     = "#AD1457",
  "5. Triple Therapy\n(LAMA/LABA/ICS — IMPACT/ETHOS)"     = "#E65100",
  "6. Triple + Roflumilast\n(Severe COPD + Chronic Bronchitis)" = "#6A1B9A"
)

theme_copd <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#37474F", color = NA),
    strip.text       = element_text(color = "white", face = "bold"),
    legend.position  = "bottom",
    legend.key.width = unit(1.2, "cm"),
    panel.grid.minor = element_blank()
  )

# A) FEV1 % predicted over time
p_fev1 <- ggplot(results, aes(x = time_day, y = FEV1, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(30, 50, 80), linetype = "dashed",
             color = c("#C62828","#F57F17","#1B5E20"), alpha = 0.6) +
  annotate("text", x = 370, y = c(28, 48, 82),
           label = c("GOLD 4 threshold","GOLD 3 threshold","GOLD 1 threshold"),
           size = 2.8, hjust = 1, color = c("#C62828","#F57F17","#1B5E20")) +
  scale_color_manual(values = scen_colors) +
  labs(title = "FEV1 Trajectory (1 year)",
       x = "Day", y = "FEV1 (% predicted)",
       color = "Treatment Scenario") +
  theme_copd

# B) Sputum IL-8
p_IL8 <- ggplot(results, aes(x = time_day, y = IL8, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 125, linetype = "dashed", color = "#666666", alpha = 0.7) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Sputum IL-8 (pg/mL)",
       x = "Day", y = "IL-8 (pg/mL)", color = NULL) +
  theme_copd + theme(legend.position = "none")

# C) Serum CRP
p_CRP <- ggplot(results, aes(x = time_day, y = CRP, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Serum CRP (mg/L)",
       x = "Day", y = "CRP (mg/L)", color = NULL) +
  theme_copd + theme(legend.position = "none")

# D) Blood eosinophils
p_Eos <- ggplot(results, aes(x = time_day, y = Eos, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(100, 300), linetype = "dashed",
             color = c("#B71C1C","#1B5E20"), alpha = 0.6) +
  annotate("text", x = 370, y = c(95, 310),
           label = c("Eos 100 (ICS benefit threshold)","Eos 300 (ICS response)"),
           size = 2.6, hjust = 1, color = c("#B71C1C","#1B5E20")) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Blood Eosinophils (cells/µL)",
       x = "Day", y = "Eosinophils (cells/µL)", color = NULL) +
  theme_copd + theme(legend.position = "none")

# E) Cumulative exacerbations
p_AE <- ggplot(results, aes(x = time_day, y = AE_cum, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Cumulative Exacerbation Count (1 yr)",
       x = "Day", y = "Cumulative exacerbations", color = NULL) +
  theme_copd + theme(legend.position = "none")

# F) Emphysema progression
p_emph <- ggplot(results, aes(x = time_day, y = Emph, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Emphysema Index (% LAA on HRCT)",
       x = "Day", y = "Emphysema index (%)", color = NULL) +
  theme_copd + theme(legend.position = "none")

# G) PVR
p_PVR <- ggplot(results, aes(x = time_day, y = PVR, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 240, linetype = "dashed", color = "#E53935", alpha = 0.6) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Pulmonary Vascular Resistance (dyn·s·cm⁻⁵)",
       x = "Day", y = "PVR", color = NULL) +
  theme_copd + theme(legend.position = "none")

# H) CAT score proxy
p_CAT <- ggplot(results, aes(x = time_day, y = CAT_approx, color = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(10, 20, 30), linetype = "dashed",
             color = c("#43A047","#F9A825","#E53935"), alpha = 0.7) +
  scale_color_manual(values = scen_colors) +
  labs(title = "CAT Score Approximation",
       x = "Day", y = "CAT score", color = NULL) +
  theme_copd + theme(legend.position = "none")

# Combine panels
fig1 <- p_fev1 / (p_IL8 | p_CRP | p_Eos) /
        (p_AE  | p_emph | p_PVR) /
        p_CAT
fig1 <- fig1 + plot_annotation(
  title    = "COPD QSP Model — 6 Treatment Scenarios over 1 Year (refactored)",
  subtitle = "Scenarios: Placebo · LAMA · LABA+LAMA · ICS+LABA · Triple · Triple+PDE4I",
  caption  = "Model calibrated to UPLIFT / TORCH / FLAME / IMPACT / ETHOS trial data.\nICS benefit gated by blood eosinophil count (Eos ≥300 → full benefit).",
  theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                   plot.subtitle = element_text(size = 10),
                   plot.caption  = element_text(size = 8, hjust = 0))
)

print(fig1)

# =============================================================================
# SUMMARY TABLE — END-OF-YEAR VALUES
# =============================================================================
summary_tbl <- results %>%
  filter(time_day >= 364.5) %>%
  group_by(scenario) %>%
  summarise(
    FEV1_yr1         = round(mean(FEV1), 1),
    GOLD_stage_yr1   = round(mean(GOLD_stage), 0),
    IL8_yr1          = round(mean(IL8), 1),
    CRP_yr1          = round(mean(CRP), 2),
    Eos_yr1          = round(mean(Eos), 0),
    AE_cum_yr1       = round(max(AE_cum), 2),
    AE_rate_yr1      = round(mean(AE_rate_ann), 2),
    Emph_yr1         = round(mean(Emph), 1),
    PVR_yr1          = round(mean(PVR), 0),
    CAT_yr1          = round(mean(CAT_approx), 1),
    SpO2_yr1         = round(mean(SpO2), 1),
    .groups = "drop"
  ) %>%
  rename(
    "Scenario"           = scenario,
    "FEV1 % pred"        = FEV1_yr1,
    "GOLD Stage"         = GOLD_stage_yr1,
    "IL-8 (pg/mL)"       = IL8_yr1,
    "CRP (mg/L)"         = CRP_yr1,
    "Eos (cells/µL)"     = Eos_yr1,
    "Cum AEs"            = AE_cum_yr1,
    "AE rate/yr"         = AE_rate_yr1,
    "Emphysema %"        = Emph_yr1,
    "PVR (dyn)"          = PVR_yr1,
    "CAT score"          = CAT_yr1,
    "SpO2 (%)"           = SpO2_yr1
  )

print(summary_tbl, n = Inf)

# =============================================================================
# DOSE-RESPONSE: ICS dose effect on sputum IL-8 at 24 weeks
# =============================================================================
ICS_doses_ug <- c(0, 50, 100, 200, 400, 800, 1600)
dose_response <- lapply(ICS_doses_ug, function(d) {
  lung_dose <- d * 0.22 * 1000   # ng
  dose_ev   <- ev(cmt = "GUT_ICS", amt = lung_dose, ii = 12,
                  addl = (12 * 7 * 24 / 12 - 1), time = 0)  # 12 weeks bid
  out <- mrgsim(
    param(mod_copd_refactored, list(DOSE_LAMA=0, DOSE_LABA=0, DOSE_ICS=1, DOSE_PDE4I=0,
                         Eos0=350)),  # eosinophilic phenotype (renamed from Eos_0, see model $MAIN note)
    events = dose_ev,
    tgrid  = seq(0, 2016, by = 24),   # 12 weeks
    delta  = 24
  ) %>% as.data.frame()
  data.frame(
    ICS_dose_ug = d,
    IL8_wk12    = mean(tail(out$IL8, 5)),
    Eos_wk12    = mean(tail(out$Eos, 5))
  )
}) %>% bind_rows()

p_DR <- ggplot(dose_response, aes(x = ICS_dose_ug, y = IL8_wk12)) +
  geom_line(color = "#880E4F", linewidth = 1.2) +
  geom_point(color = "#880E4F", size = 3) +
  scale_x_log10(
    breaks = c(50, 100, 200, 400, 800, 1600),
    labels = c("50","100","200","400","800","1600")
  ) +
  labs(
    title    = "ICS Dose-Response: Sputum IL-8 at 12 Weeks (Eosinophilic Phenotype)",
    subtitle = "Eosinophilic COPD (Eos ≥350 cells/µL); ICS active via GR–NF-κB transrepression",
    x        = "ICS Dose (µg/day, log scale)",
    y        = "Sputum IL-8 at Week 12 (pg/mL)",
    caption  = "EC50 for ICS PD = 0.15 ng/mL lung; Hill = 1.5. Calibrated to Bleecker 2017 (IMPACT)."
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_DR)

cat("\nCOPD mrgsolve QSP Model (refactored) — Simulation Complete\n")
cat("6 scenarios × 1 year × 19 ODE states\n")
cat("Clinical trial calibrations: UPLIFT / TORCH / FLAME / IMPACT / ETHOS\n")
