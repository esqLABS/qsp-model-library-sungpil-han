## ============================================================
## Familial Hypercholesterolemia (FH) QSP Model
## REFACTORED sibling of fh_mrgsolve_model.R -- pluggable-PK naming
## convention applied to the file's three census-tracked compounds:
## Rosuvastatin (STATIN stem), Evolocumab (PCSK9I stem), Ezetimibe (EZE
## stem). Lomitapide and Bempedoic acid are flag-driven (no PK
## compartments of their own) and are byte-for-byte identical to the
## original everywhere they appear. All numeric parameter values are
## copied verbatim from the original; only PK/PD naming and structure
## are reorganized, plus:
##   - five pre-existing mrgsolve-2.0.1 build defects (unrelated to any
##     compound's own PK/PD shape) found and fixed with syntax-only,
##     non-numeric changes so this file actually compiles through the
##     qspserver mrgsolve_api -- logged in translations/UPSTREAM_ISSUES.md
##     and disclosed in full in fh_refactor_notes.md.
##   - three new, math-implied Hill constants per compound
##     (GAMMA_STATIN=1, GAMMA_EZE=1, EC50_EZE promoted out of an embedded
##     literal) -- renames/promotions, not refits; see notes.
## The original (fh_mrgsolve_model.R) is untouched.
## Refactor added 2026-09-02.
## ============================================================
##
## Disease:   Familial Hypercholesterolemia (LDLR/APOB/PCSK9)
## Drugs:     Statins (Rosuvastatin), PCSK9 inhibitors (Evolocumab),
##            Ezetimibe, Inclisiran (unmodeled flag only),
##            Bempedoic acid, Lomitapide
## Scenarios: HetFH baseline, statin mono, statin+EZE,
##            PCSK9i mono, statin+PCSK9i, lomitapide (HomFH)
##
## References: Stecula et al. 2019, Harrington et al. 2023,
##             FOURIER trial 2017, ORION trials 2020-2022,
##             Watts et al. 2020 (FH global guidelines)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ----------------------------------------------------------
## Model code
## ----------------------------------------------------------
fh_code <- '
$PARAM
// === Statin PK (Rosuvastatin-like, oral, mg -> ug/mL) ===
// Bespoke structure (not a plain archetype): GUT_STATIN feeds BOTH
// CENT_STATIN (systemic) and LIV_STATIN (hepatic site of action) in
// PARALLEL as two independent first-pass splits of the same absorbed
// dose -- there is no Q-mediated exchange between CENT_STATIN and
// LIV_STATIN, unlike a true central/peripheral pair. This is the
// original structure exactly, renamed only. See refactor notes.
KA_STATIN     = 0.45,    // GI absorption rate constant (h^-1) (was ka_S)
F_STATIN      = 0.20,    // oral bioavailability (hepatic first-pass ~80%) (was F_S)
V1_STATIN     = 2.0,     // apparent volume, PER KG body weight -- combined with BW at each use site, matching the original convention (was Vd_S)
CL_STATIN     = 0.60,    // systemic (non-hepatic) clearance (L/h) (was CL_S)
CL_HEP_STATIN = 0.45,    // hepatic clearance (L/h) -- predominant route (was CL_H)
DOSE_STATIN   = 0.0,     // statin daily dose (mg), 0 = no statin; NOT read by $ODE -- dosing is via GUT_STATIN events, this is R-side bookkeeping only, same as the original (was DOSE_S)
EC50_STATIN   = 0.08,    // statin hepatic conc for 50% HMGCR inhibition (ug/mL) (was EC50_SI)
EMAX_STATIN   = 0.92,    // maximum fraction HMGCR inhibition by statin (was Emax_SI)
GAMMA_STATIN  = 1,       // Hill coefficient [new; no explicit Hill exponent in the original, so gamma=1 -- a rename, not a refit]

// === PCSK9 inhibitor PK (Evolocumab-like, subcutaneous mAb) -- archetype 4, TMDD ===
// The TMDD "free target" is PCSK9_pl, the disease-side circulating PCSK9
// pool (declared under PD compartments below, not renamed to this
// compound’s stem -- see refactor notes for why).
KA_PCSK9I    = 0.0245,  // SC absorption rate (h^-1) ~ bioavail/t_abs (was ka_P)
F_PCSK9I     = 0.72,    // SC bioavailability evolocumab (was F_P)
V1_PCSK9I    = 3.1,     // central volume (L) (was Vc_P)
V2_PCSK9I    = 2.6,     // peripheral volume (L) (was Vp_P)
CL_PCSK9I    = 0.0060,  // linear clearance evolocumab (L/h) (was CL_P)
Q_PCSK9I     = 0.0043,  // intercompartment clearance (L/h) (was Q_P)
KTMDD_PCSK9I = 0.00015, // target-mediated (PCSK9-bound) clearance (L/h/nmol) (was kTMDD)
DOSE_PCSK9I  = 0.0,     // PCSK9i dose (mg), 0 = none [420 mg q4w]; NOT read by $ODE -- dosing is via GUT_PCSK9I events, R-side bookkeeping only (was DOSE_P)
MW_PCSK9I    = 144000,  // molecular weight evolocumab (Da) (was MW_P)
KON_PCSK9I   = 0.012,   // PCSK9i binding to free PCSK9 (1/nM/h) (was kon_PK9)
KOFF_PCSK9I  = 1e-5,    // PCSK9i-PCSK9 dissociation (h^-1) (was koff_PK9)
// Hill-interface parameters for EFFECT_PCSK9I. NOT a fit: EFFECT_PCSK9I
// below is the exact ODE-state fractional PCSK9 target engagement
// (COMPLEX_PCSK9I vs. total PCSK9 currently in the system) -- the same
// TMDD-occupancy pattern used for tocilizumab elsewhere in this corpus
// (rheumatoid-arthritis). EC50_PCSK9I is the analytically derived
// quasi-steady-state Kd, (KOFF_PCSK9I+kout_PK9)/KON_PCSK9I -- exposed for
// discoverability/future closed-form use, NOT consumed by EFFECT_PCSK9I
// itself (same treatment as EC50_TCZ in the RA precedent). See notes.
EMAX_PCSK9I  = 1,          // occupancy saturates at 1 by definition [new]
EC50_PCSK9I  = 2.000833,   // = (KOFF_PCSK9I+kout_PK9)/KON_PCSK9I, nM [new, derived, not fit]
GAMMA_PCSK9I = 1,          // [new]

// === Ezetimibe PK (oral, mg) -- archetype 3 without a peripheral compartment (depot + central only) ===
KA_EZE   = 0.35,   // absorption rate (h^-1) (was ka_EZE)
F_EZE    = 0.35,   // oral bioavailability (active + inactive EH recycling)
V1_EZE   = 2.5,    // volume of distribution, PER KG body weight (was Vd_EZE)
CL_EZE   = 0.15,   // systemic clearance (L/h)
DOSE_EZE = 0.0,    // ezetimibe dose (mg, typically 10 mg/day); NOT read by $ODE -- dosing is via GUT_EZE events, R-side bookkeeping only
EC50_EZE = 0.05,   // ezetimibe conc for 50% intestinal cholesterol-absorption block (mg/L) [new; was an embedded literal 0.05 inside the original ratio -- a promotion, not a refit]
EMAX_EZE = 0.55,   // max fraction reduction in intestinal cholesterol absorption (was fEZE_abs)
GAMMA_EZE = 1,     // Hill coefficient [new; no explicit Hill exponent in the original, so gamma=1]

// === HMGCR/LDLR turnover === (disease-side, unchanged)
kout_HMG = 0.28,   // HMGCR turnover rate (h^-1) ~2.5 h half-life
kin_LR   = 0.025,  // LDLR synthesis rate (h^-1, relative to baseline)
kout_LR  = 0.020,  // LDLR degradation baseline (h^-1)
EC50_LD  = 0.30,   // HMGCR inhib fraction for 50% LDLR upregulation
Emax_LD  = 2.5,    // max fold-increase in LDLR synthesis due to SREBP2

// === PCSK9 dynamics === (disease-side, unchanged)
kin_PK9  = 10.0,   // PCSK9 synthesis rate (ng/mL/h) -> steady-state ~400 ng/mL
kout_PK9 = 0.024,  // PCSK9 clearance (h^-1), t1/2 ~29 h -- also the COMPLEX_PCSK9I clearance rate (same role as KDEG in a classic TMDD receptor pool)
EC50_PK9_LR = 200, // PCSK9 conc for 50% LDLR degradation (ng/mL)
kP9_LDLR = 0.010,  // PCSK9-mediated extra LDLR degradation (per unit PCSK9)

// === Lipoprotein dynamics (mg/dL) === (disease-side, unchanged)
LDL0_het = 280,    // HetFH baseline LDL-C (mg/dL)
LDL0_hom = 680,    // HomFH baseline LDL-C (mg/dL)
VLDL0    = 32,     // baseline VLDL-C (mg/dL)
IDL0     = 15,     // baseline IDL-C (mg/dL)
HDL0     = 46,     // baseline HDL-C (mg/dL)
TG0      = 155,    // baseline TG (mg/dL)
PCSK9_0  = 400,    // baseline plasma PCSK9 (ng/mL)
LDLR_0   = 1.0,    // baseline LDLR (normalized, =1)

// VLDL kinetics
kVLDL_s  = 0.040,  // VLDL -> IDL conversion (h^-1)
kIDL_s   = 0.085,  // IDL -> LDL conversion (h^-1)
kLDL_cl  = 0.0065, // LDL clearance baseline (h^-1)
kLDL_cl_LR = 0.028,// Additional LDL clearance per LDLR unit (h^-1)
kHDL_eq  = 0.012,  // HDL equilibration rate
fLomit   = 0.80,   // max MTP inhibition by lomitapide (-> down VLDL)
DOSE_LOMT= 0.0,    // Lomitapide dose flag (1 = on)
DOSE_BEMP= 0.0,    // Bempedoic acid flag (1 = on, adds ~18% LDL reduction)
fBemp    = 0.18,   // max additional LDL-C reduction by bempedoic acid

// === Genetic parameters ===
LDLR_fxn = 0.50,   // LDLR residual function (0.5=HetFH, 0.05=HomFH, 1=normal)

// === Body parameters ===
BW       = 75.0,   // body weight (kg)
Vmix     = 4.0     // plasma mixing volume for LP (L)

$CMT
// Statin PK (bespoke parallel-split structure -- see $PARAM header comment)
GUT_STATIN         // Statin GI (mg) (was GUT_S)
CENT_STATIN        // Statin central plasma (ug) (was CENT_S)
LIV_STATIN         // Statin liver (ug) -- primary site of action (was LIV_S)

// PCSK9 inhibitor PK (archetype 4, TMDD)
GUT_PCSK9I         // PCSK9i subcutaneous depot (mg) (was SC_PCSK9I)
CENT_PCSK9I        // PCSK9i central plasma (mg)
PERI_PCSK9I        // PCSK9i peripheral compartment (mg)
COMPLEX_PCSK9I     // PCSK9i-PCSK9 complex (nmol/L equiv) (was COMP_PK9)

// Ezetimibe PK (archetype 3 minus peripheral: depot + central only)
GUT_EZE            // Ezetimibe GI (mg)
CENT_EZE           // Ezetimibe plasma (mg)

// PD compartments (disease-side, unchanged)
HMGCR_rel          // HMGCR relative activity (1 = normal)
LDLR_rel           // LDLR surface level (1 = baseline for genotype)
PCSK9_pl           // Plasma PCSK9 (ng/mL) -- disease-side free/total circulating target; deliberately NOT renamed to a PCSK9I_ stem, see refactor notes

// Lipoprotein compartments (mg/dL)
VLDL_C             // VLDL-C
IDL_C              // IDL-C
LDL_C              // LDL-C
HDL_C              // HDL-C
TG_C               // Triglycerides

$MAIN
// Initial conditions, set here via the <cmt>_0 idiom rather than a
// separate $INIT block. This is a mandatory, syntax-only build fix, not
// a modeling change -- see refactor notes and the new UPSTREAM_ISSUES.md
// entry for why: this build (a) evaluates $INIT as one list() call with
// no access to $PARAM values (fails on any $INIT line that references a
// sibling $PARAM, e.g. "PCSK9_pl = PCSK9_0"), and (b) independently
// rejects any model where $CMT and a separate $INIT block jointly
// declare the same compartments ("Duplicated model names"). Every value
// below is copied verbatim from the original’s own $INIT block.
GUT_STATIN_0     = 0;
CENT_STATIN_0    = 0;
LIV_STATIN_0     = 0;
GUT_PCSK9I_0     = 0;
CENT_PCSK9I_0    = 0;
PERI_PCSK9I_0    = 0;
COMPLEX_PCSK9I_0 = 0;
GUT_EZE_0        = 0;
CENT_EZE_0       = 0;
HMGCR_rel_0      = 1.0;
LDLR_rel_0       = 1.0;
PCSK9_pl_0       = PCSK9_0;
VLDL_C_0         = VLDL0;
IDL_C_0          = IDL0;
LDL_C_0          = LDL0_het;
HDL_C_0          = HDL0;
TG_C_0           = TG0;

$ODE

// -------------------------------------------------------
// [1-3] STATIN PK (bespoke parallel-split structure)
// -------------------------------------------------------
double C_STATIN     = LIV_STATIN / (V1_STATIN * BW * 0.26);  // hepatic conc (ug/mL) -- SITE OF ACTION, the single canonical PD-facing concentration
double C_STATIN_SYS = CENT_STATIN / (V1_STATIN * BW);        // systemic conc (ug/mL) -- reporting only; no PD equation reads this (was C_sys_S / C_statin_sys, now a single site)

dxdt_GUT_STATIN  = -KA_STATIN * GUT_STATIN;
dxdt_CENT_STATIN =  KA_STATIN * F_STATIN * GUT_STATIN - (CL_STATIN / (V1_STATIN * BW)) * CENT_STATIN
                    - (CL_HEP_STATIN / (V1_STATIN * BW * 0.26)) * CENT_STATIN;
dxdt_LIV_STATIN  =  KA_STATIN * F_STATIN * GUT_STATIN * 0.80   // first-pass hepatic extraction ~80%
                    - (CL_HEP_STATIN / (V1_STATIN * BW * 0.26)) * LIV_STATIN;

// -------------------------------------------------------
// [4-7] PCSK9 INHIBITOR PK (archetype 4, TMDD)
// -------------------------------------------------------
double C_PCSK9I  = CENT_PCSK9I / V1_PCSK9I;          // mg/L -- canonical exposed concentration, single site (was C_pcsk9i_mgL, computed only in $TABLE; now the one site both PK and reporting read)
double PCSK9I_nM = C_PCSK9I / (MW_PCSK9I * 1e-6);    // nM; algebraically identical to CENT_PCSK9I/(MW_PCSK9I*1e-6*V1_PCSK9I) -- derived from the canonical C_PCSK9I instead of being recomputed independently from CENT_PCSK9I a second time (was PCSK9i_nM)

dxdt_GUT_PCSK9I     = -KA_PCSK9I * GUT_PCSK9I;
dxdt_CENT_PCSK9I    =  KA_PCSK9I * F_PCSK9I * GUT_PCSK9I
                       - (CL_PCSK9I / V1_PCSK9I) * CENT_PCSK9I
                       - (Q_PCSK9I  / V1_PCSK9I) * CENT_PCSK9I
                       + (Q_PCSK9I  / V2_PCSK9I) * PERI_PCSK9I
                       - KTMDD_PCSK9I * PCSK9_pl * CENT_PCSK9I;
dxdt_PERI_PCSK9I    = (Q_PCSK9I / V1_PCSK9I) * CENT_PCSK9I
                      - (Q_PCSK9I / V2_PCSK9I) * PERI_PCSK9I;
dxdt_COMPLEX_PCSK9I =  KON_PCSK9I * PCSK9I_nM * (PCSK9_pl / 1000.0)
                       - KOFF_PCSK9I * COMPLEX_PCSK9I
                       - kout_PK9 * COMPLEX_PCSK9I;  // complex is cleared

// Hill interface (TMDD): fractional PCSK9 target engagement, an EXACT
// ODE-state ratio -- not a fit. Total currently-circulating PCSK9
// (free+bound) is PCSK9_pl + COMPLEX_PCSK9I*1000, using the same *1000
// unit-equivalence factor the original already applies wherever it
// combines these two state variables (see PCSK9_free below). This
// diagnostic does NOT feed back into any disease equation -- the
// disease side still consumes the native mass-action reduction in
// PCSK9_free exactly as the original did. See refactor notes.
double OCC_PCSK9I    = (COMPLEX_PCSK9I * 1000.0) / (PCSK9_pl + COMPLEX_PCSK9I * 1000.0);
double EFFECT_PCSK9I = OCC_PCSK9I;

// -------------------------------------------------------
// [8-9] EZETIMIBE PK (archetype 3 without a peripheral compartment)
// -------------------------------------------------------
double C_EZE = CENT_EZE / (V1_EZE * BW);  // mg/L -- canonical exposed concentration, single site

dxdt_GUT_EZE  = -KA_EZE * GUT_EZE;
dxdt_CENT_EZE =  KA_EZE * F_EZE * GUT_EZE
                 - (CL_EZE / (V1_EZE * BW)) * CENT_EZE;

double EFFECT_EZE = EMAX_EZE * pow(C_EZE, GAMMA_EZE)
                     / (pow(EC50_EZE, GAMMA_EZE) + pow(C_EZE, GAMMA_EZE));

// -------------------------------------------------------
// [10] HMGCR ACTIVITY (inhibited by hepatic statin)
// -------------------------------------------------------
double EFFECT_STATIN = EMAX_STATIN * pow(C_STATIN, GAMMA_STATIN)
                       / (pow(EC50_STATIN, GAMMA_STATIN) + pow(C_STATIN, GAMMA_STATIN));
double kin_HMG = kout_HMG;   // at steady state HMGCR_rel = 1
double Bemp_eff = DOSE_BEMP * fBemp * 0.55;  // bempedoic acid ACLY->HMGCR effect

dxdt_HMGCR_rel = kin_HMG * (1.0 - EFFECT_STATIN - Bemp_eff) - kout_HMG * HMGCR_rel;

// -------------------------------------------------------
// [11] LDLR SURFACE EXPRESSION
// -------------------------------------------------------
// SREBP2 activated when HMGCR is inhibited (down intracellular cholesterol)
double HMGCR_inh_frac = 1.0 - HMGCR_rel;   // fraction inhibited 0->1
double LDLR_up = Emax_LD * HMGCR_inh_frac / (EC50_LD + HMGCR_inh_frac); // fold-increase

// PCSK9 degrades LDLR (concentration-dependent); single canonical site --
// was recomputed a second time in $TABLE under the same name, which this
// build rejects outright (duplicate "double PCSK9_free" declaration);
// now declared once here and captured directly (see refactor notes).
double PCSK9_free = PCSK9_pl - COMPLEX_PCSK9I * 1000.0; // free PCSK9 (ng/mL)
if(PCSK9_free < 0) PCSK9_free = 0;
double PCSK9_effect_LDLR = kP9_LDLR * PCSK9_free / (EC50_PK9_LR + PCSK9_free); // extra degradation

// Genotype scaling
double LDLR_max = LDLR_fxn; // peak possible = genotype-scaled
double kin_LR_eff = kin_LR * LDLR_max * (1.0 + LDLR_up);
double kout_LR_eff = kout_LR + PCSK9_effect_LDLR;

dxdt_LDLR_rel = kin_LR_eff - kout_LR_eff * LDLR_rel;

// -------------------------------------------------------
// [12] PLASMA PCSK9
// -------------------------------------------------------
// Statin increases PCSK9 transcription (SREBP2 feedback)
double PCSK9_statin_feedback = 1.0 + 0.40 * HMGCR_inh_frac; // statins increase PCSK9 ~40%
double PCSK9_inclisiran_red  = 1.0;  // flag for inclisiran (simplified); unused, kept as original

// PCSK9i binds PCSK9
double PCSK9_pl_use = PCSK9_pl;
if(PCSK9_pl_use < 0) PCSK9_pl_use = 0;

dxdt_PCSK9_pl = kin_PK9 * PCSK9_statin_feedback
                - kout_PK9 * PCSK9_pl_use
                - KON_PCSK9I * PCSK9I_nM * PCSK9_pl_use
                + KOFF_PCSK9I * COMPLEX_PCSK9I * 1000.0;

// -------------------------------------------------------
// [13-17] LIPOPROTEIN DYNAMICS
// -------------------------------------------------------
// Ezetimibe reduces intestinal cholesterol input -> down VLDL substrate
double EZE_effect_VLDL = 1.0 - EFFECT_EZE;
// Lomitapide reduces VLDL assembly
double LOMT_eff = 1.0 - DOSE_LOMT * fLomit;
// VLDL production rate
double k_VLDL_prod = kVLDL_s * VLDL0 * EZE_effect_VLDL * LOMT_eff;

dxdt_VLDL_C = k_VLDL_prod
              - kVLDL_s * VLDL_C   // VLDL -> IDL lipolysis
              - (HMGCR_inh_frac * 0.20) * VLDL_C; // statins modestly reduce VLDL

// IDL
dxdt_IDL_C  = kVLDL_s * VLDL_C
              - kIDL_s  * IDL_C;   // IDL -> LDL

// LDL clearance: baseline + LDLR-mediated + PCSK9i enhanced
double kLDL_total = kLDL_cl + kLDL_cl_LR * LDLR_rel;
// Bempedoic acid additional LDL-C reduction through HMGCR
double bemp_LDL_add = DOSE_BEMP * fBemp;

dxdt_LDL_C  = kIDL_s * IDL_C
              - kLDL_total * LDL_C
              - bemp_LDL_add * kLDL_cl * LDL_C;

// HDL: statins increase HDL ~6%, CETP inhibition not modelled separately
double HDL_target = HDL0 * (1.0 + 0.06 * HMGCR_inh_frac);
dxdt_HDL_C  = kHDL_eq * (HDL_target - HDL_C);

// TG: statins reduce TG, fibrates not modelled
double TG_target = TG0 * (1.0 - 0.15 * HMGCR_inh_frac);
dxdt_TG_C   = kHDL_eq * (TG_target - TG_C);

$TABLE
// Derived PD (disease-side, unchanged). Statin/PCSK9i/Ezetimibe
// concentrations and effects are declared once in $ODE (see above) and
// captured directly below -- no longer recomputed a second time here
// under different names (was C_statin_sys/C_statin_liv/C_pcsk9i_mgL/
// C_eze_mgL; see refactor notes for the "normalize duplicate
// concentration sites" fix this satisfies).
double LDL_reduction_pct  = (LDL0_het > 0) ? (LDL0_het - LDL_C) / LDL0_het * 100.0 : 0;
double NonHDL_C = LDL_C + VLDL_C + IDL_C;
double TC       = LDL_C + VLDL_C + HDL_C + IDL_C;
double LDLR_pct = LDLR_rel * 100.0;

// CVD risk surrogates (simplified Framingham-style)
double CVD_risk_10yr = 0.01 * exp(0.0035 * LDL_C - 0.01 * HDL_C);

// LDL goals (ESC 2019 very-high risk: <55 mg/dL; high risk: <70 mg/dL)
double LDL_goal_55  = (LDL_C <= 55) ? 1 : 0;
double LDL_goal_70  = (LDL_C <= 70) ? 1 : 0;

$CAPTURE
C_STATIN, C_STATIN_SYS, C_PCSK9I, C_EZE,
EFFECT_STATIN, EFFECT_PCSK9I, EFFECT_EZE, OCC_PCSK9I,
LDLR_pct, PCSK9_free,
NonHDL_C, TC, LDL_reduction_pct, CVD_risk_10yr,
LDL_goal_55, LDL_goal_70
'

## ----------------------------------------------------------
## Compile model
## ----------------------------------------------------------
mod <- mcode("FH_QSP_REFACTORED", fh_code)
cat("Model compiled successfully. Compartments:", length(Init(mod)), "\n")

## ----------------------------------------------------------
## Helper: build event data for daily oral + q4w SC dosing
## ----------------------------------------------------------
build_ev <- function(dose_statin_mg  = 0,
                     dose_pcsk9i_mg  = 0,    # 0 or 420 mg q4w
                     dose_eze_mg     = 0,
                     n_days          = 365) {
  evlist <- list()

  # Statin: once daily oral (add to GUT_STATIN)
  if (dose_statin_mg > 0) {
    ev_s <- ev(amt = dose_statin_mg, cmt = "GUT_STATIN", ii = 24, addl = n_days - 1)
    evlist[["statin"]] <- ev_s
  }

  # PCSK9i: q4w subcutaneous (every 28 days)
  if (dose_pcsk9i_mg > 0) {
    n_inj <- floor(n_days / 28)
    ev_p <- ev(amt = dose_pcsk9i_mg, cmt = "GUT_PCSK9I", ii = 28 * 24, addl = n_inj - 1)
    evlist[["pcsk9i"]] <- ev_p
  }

  # Ezetimibe: once daily oral
  if (dose_eze_mg > 0) {
    ev_e <- ev(amt = dose_eze_mg, cmt = "GUT_EZE", ii = 24, addl = n_days - 1)
    evlist[["eze"]] <- ev_e
  }

  if (length(evlist) == 0) return(NULL)
  Reduce(c, evlist)
}

## ----------------------------------------------------------
## Treatment Scenarios
## ----------------------------------------------------------
scenarios <- list(
  list(
    label    = "1. HetFH -- No Treatment",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = NULL,
    DOSE_PCSK9I = 0, DOSE_STATIN = 0, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#E74C3C"
  ),
  list(
    label    = "2. HetFH -- Rosuvastatin 40 mg/d",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40),
    DOSE_PCSK9I = 0, DOSE_STATIN = 40, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#E67E22"
  ),
  list(
    label    = "3. HetFH -- Rosuvastatin 40 mg + Ezetimibe 10 mg",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40, dose_eze_mg = 10),
    DOSE_PCSK9I = 0, DOSE_STATIN = 40, DOSE_EZE = 10,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#F39C12"
  ),
  list(
    label    = "4. HetFH -- Evolocumab 420 mg q4w (PCSK9i mono)",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_pcsk9i_mg = 420),
    DOSE_PCSK9I = 420, DOSE_STATIN = 0, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#27AE60"
  ),
  list(
    label    = "5. HetFH -- Rosuvastatin 40 mg + Evolocumab 420 mg q4w",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40, dose_pcsk9i_mg = 420),
    DOSE_PCSK9I = 420, DOSE_STATIN = 40, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#2E86C1"
  ),
  list(
    label    = "6. HomFH -- Lomitapide + Statin + Evolocumab",
    LDLR_fxn = 0.05,
    LDL0     = 680,
    ev_      = build_ev(dose_statin_mg = 40, dose_pcsk9i_mg = 420),
    DOSE_PCSK9I = 420, DOSE_STATIN = 40, DOSE_EZE = 0,
    DOSE_LOMT= 1, DOSE_BEMP = 0,
    color    = "#8E44AD"
  )
)

## ----------------------------------------------------------
## Run simulations
## ----------------------------------------------------------
sim_results <- lapply(scenarios, function(sc) {
  # Update initial LDL for genotype
  init_vals <- Init(mod)
  init_vals["LDL_C"]    <- sc$LDL0
  init_vals["LDLR_rel"] <- sc$LDLR_fxn
  init_vals["PCSK9_pl"] <- 400

  params_update <- list(
    LDLR_fxn    = sc$LDLR_fxn,
    LDL0_het    = sc$LDL0,
    DOSE_STATIN = sc$DOSE_STATIN,
    DOSE_PCSK9I = sc$DOSE_PCSK9I,
    DOSE_EZE    = sc$DOSE_EZE,
    DOSE_LOMT   = sc$DOSE_LOMT,
    DOSE_BEMP   = sc$DOSE_BEMP
  )

  m2 <- mod %>%
    param(params_update) %>%
    init(init_vals)

  out <- if (!is.null(sc$ev_)) {
    mrgsim(m2, sc$ev_, end = 365 * 24, delta = 12, carry_out = "evid")
  } else {
    mrgsim(m2, end = 365 * 24, delta = 12)
  }

  as.data.frame(out) %>%
    mutate(Scenario = sc$label, Color = sc$color,
           time_days = time / 24)
})

results_df <- bind_rows(sim_results)

## ----------------------------------------------------------
## Plot 1: LDL-C over 52 weeks
## ----------------------------------------------------------
p1 <- ggplot(results_df, aes(time_days, LDL_C, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "navy", alpha = 0.7) +
  geom_hline(yintercept = 55, linetype = "dotted", color = "darkred", alpha = 0.7) +
  annotate("text", x = 350, y = 75,  label = "70 mg/dL (ESC high risk)",
           color = "navy", size = 3, hjust = 1) +
  annotate("text", x = 350, y = 60,  label = "55 mg/dL (ESC very high risk)",
           color = "darkred", size = 3, hjust = 1) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Familial Hypercholesterolemia -- LDL-C Time Course",
       x = "Time (days)", y = "LDL-C (mg/dL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 2))

## ----------------------------------------------------------
## Plot 2: LDLR Surface Expression
## ----------------------------------------------------------
p2 <- ggplot(results_df, aes(time_days, LDLR_pct, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Hepatic LDLR Surface Expression",
       x = "Time (days)", y = "LDLR (%  of normal)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 3: Plasma PCSK9 levels
## ----------------------------------------------------------
p3 <- ggplot(results_df, aes(time_days, PCSK9_pl, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Plasma PCSK9 Concentration",
       x = "Time (days)", y = "PCSK9 (ng/mL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 4: Lipid Panel Summary at Week 52
## ----------------------------------------------------------
summary_52w <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    LDL_C    = mean(LDL_C),
    HDL_C    = mean(HDL_C),
    TG       = mean(TG_C),
    NonHDL_C = mean(NonHDL_C),
    .groups  = "drop"
  ) %>%
  pivot_longer(cols = c(LDL_C, HDL_C, TG, NonHDL_C),
               names_to = "Lipid", values_to = "mgdL")

p4 <- ggplot(summary_52w, aes(Scenario, mgdL, fill = Lipid)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Lipid Panel at Week 52",
       x = "", y = "Concentration (mg/dL)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

## ----------------------------------------------------------
## Plot 5: % LDL Reduction waterfall
## ----------------------------------------------------------
ldl_reduc <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario, Color) %>%
  summarise(ldl_pct = mean(LDL_reduction_pct), .groups = "drop") %>%
  arrange(desc(ldl_pct))

p5 <- ggplot(ldl_reduc, aes(reorder(Scenario, ldl_pct), ldl_pct, fill = Scenario)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", ldl_pct)), hjust = -0.1, size = 4) +
  scale_fill_manual(values = setNames(ldl_reduc$Color, ldl_reduc$Scenario)) +
  coord_flip(ylim = c(-50, 90)) +
  labs(title = "LDL-C Reduction at Week 52\n(vs. HetFH untreated baseline 280 mg/dL)",
       x = "", y = "% LDL-C Reduction") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 6: HMGCR Activity over time
## ----------------------------------------------------------
p6 <- ggplot(results_df, aes(time_days, HMGCR_rel * 100, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "HMGCR Relative Activity",
       x = "Time (days)", y = "HMGCR Activity (% baseline)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Combine & save
## ----------------------------------------------------------
combined_plot <- (p1 + p5) / (p2 + p3) / (p4 + p6) +
  plot_annotation(
    title    = "Familial Hypercholesterolemia (FH) -- QSP Simulation Dashboard (Refactored)",
    subtitle = "mrgsolve ODE model: PK/PD of statin, PCSK9 inhibitor, ezetimibe, lomitapide",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 12))
  )

ggsave("fh_simulation_dashboard_refactored.png", combined_plot,
       width = 18, height = 16, dpi = 150)
cat("Dashboard saved: fh_simulation_dashboard_refactored.png\n")

## ----------------------------------------------------------
## Steady-state sensitivity: LDL reduction vs statin dose
## ----------------------------------------------------------
statin_doses <- c(5, 10, 20, 40, 80)

ss_results <- lapply(statin_doses, function(d) {
  m2 <- mod %>%
    param(LDLR_fxn = 0.50, LDL0_het = 280, DOSE_STATIN = d) %>%
    init(LDL_C = 280, LDLR_rel = 0.50)
  ev_ <- build_ev(dose_statin_mg = d)
  out <- mrgsim(m2, ev_, end = 180 * 24, delta = 24)
  tail(as.data.frame(out), 1) %>%
    mutate(Dose = d, LDL_reduction = LDL_reduction_pct)
})
ss_df <- bind_rows(ss_results)

cat("\nSteady-state LDL-C reduction by statin dose (HetFH):\n")
cat(sprintf("  Dose %3d mg/d -> LDL-C = %5.1f mg/dL (down %.0f%%)\n",
            ss_df$Dose, ss_df$LDL_C, ss_df$LDL_reduction))

## ----------------------------------------------------------
## Print summary table
## ----------------------------------------------------------
cat("\n=== LDL-C Summary at Week 52 ===\n")
summ <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    LDL_C_mgdL = round(mean(LDL_C), 1),
    HDL_C_mgdL = round(mean(HDL_C), 1),
    TG_mgdL    = round(mean(TG_C), 1),
    PCSK9_ngmL = round(mean(PCSK9_pl), 0),
    LDLR_pct   = round(mean(LDLR_pct), 0),
    LDL_red_pct= round(mean(LDL_reduction_pct), 1),
    LDL_goal55 = round(mean(LDL_goal_55) * 100),
    LDL_goal70 = round(mean(LDL_goal_70) * 100),
    .groups    = "drop"
  )
print(as.data.frame(summ))
