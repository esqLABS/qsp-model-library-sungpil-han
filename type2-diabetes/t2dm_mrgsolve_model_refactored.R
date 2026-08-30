## =============================================================================
##  Type 2 Diabetes Mellitus — mrgsolve QSP Model (REFACTORED)
##  Disease: T2DM  |  Version: 1.0-refactored  |  Date: 2026-08-30
##
##  This is a PK/PD naming-convention refactor of `t2dm_mrgsolve_model.R`
##  (fork workflow, FORK_WORKFLOW_GUIDE.md Part 2). Every one of the seven
##  compounds' PK/PD blocks was renamed to the fork's pluggable-PK convention
##  (GUT_<STEM>/CENT_<STEM>/PERI_<STEM>, C_<STEM>, EFFECT_<STEM>, etc.) with
##  the disease-side glucose/insulin/glucagon/GLP-1/beta-cell/IR/FFA/HbA1c/
##  eGFR/UACR machinery left otherwise byte-for-behavior identical. See
##  `t2dm_refactor_notes.md` for the full archetype-by-archetype rationale,
##  the three build-compatibility fixes applied (also logged in
##  translations/UPSTREAM_ISSUES.md #99-#100), and the verification results.
##
##  Model structure (unchanged from the original):
##    - 7 drug PK compartments (Metformin, Empagliflozin, Semaglutide,
##      Sitagliptin, Glimepiride, Insulin Degludec, Pioglitazone)
##    - Glucose-insulin minimal model with tissue compartment (Bergman ext.)
##    - Glucagon dynamics
##    - β-cell mass/function trajectories (adapted from Topp et al., 2000)
##    - Incretin (GLP-1) compartment with DPP-4 degradation
##    - Hepatic & peripheral insulin resistance indices
##    - FFA / adiposity / body weight
##    - HbA1c, eGFR, UACR endpoints
##
##  Treatment scenarios (7, unchanged):
##    1. No treatment (diet/exercise baseline)
##    2. Metformin 1000 mg BID
##    3. Metformin + Empagliflozin 10 mg QD
##    4. Metformin + Semaglutide 1 mg SC weekly
##    5. Triple: Metformin + Empagliflozin + Semaglutide
##    6. Metformin + Insulin Degludec 20 U QD
##    7. Metformin + Sitagliptin 100 mg QD (DPP-4i)
##
##  Parameter calibration notes (unchanged from the original):
##    - Glucose-insulin: Bergman minimal model (Bergman 1989, Diabetes)
##    - β-cell mass: Topp et al. (2000), J Theor Biol 244:501
##    - SGLT2i PD: Ferrannini et al. (2012), Diabetes Care
##    - GLP-1RA PK: Lau et al. (2015), Pharm Res (semaglutide)
##    - Metformin PK: Gong et al. (2012), Clin Pharmacokinet
##    - DECLARE-TIMI 58 (Wiviott 2019): empagliflozin → HbA1c ↓0.7%
##    - SUSTAIN-6 (Marso 2016): semaglutide → HbA1c ↓1.4%, wt ↓4.5 kg
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# =============================================================================
# MODEL DEFINITION
# =============================================================================
t2dm_model_code <- '
$PROB T2DM QSP Model — Multi-Drug PK/PD (refactored: pluggable PK naming convention)

$PARAM
// ---- Patient baseline ----
BW      = 90,    // body weight (kg)
HbA1c0  = 9.0,   // baseline HbA1c (%)
Gp0     = 180,   // baseline plasma glucose (mg/dL)
Ip0     = 20,    // baseline plasma insulin (mU/L)
Gc0     = 150,   // baseline glucagon (pg/mL)
eGFR0   = 80,    // baseline eGFR (mL/min/1.73m2)
UACR0   = 50,    // baseline UACR (ug/mg)

// ---- Drug dosing switches (1=on, 0=off) ----
USE_MET  = 0,    // Metformin
USE_EMPA = 0,    // Empagliflozin
USE_SEMA = 0,    // Semaglutide
USE_DPP4 = 0,    // Sitagliptin
USE_SU   = 0,    // Glimepiride
USE_INS  = 0,    // Insulin Degludec
USE_PIOG = 0,    // Pioglitazone

// ---- Metformin PK (Gong 2012, Clin Pharmacokinet) ----
// Archetype 3 (depot + central + peripheral, linear). Original parameterized
// distribution as CL/V elimination + k12/k21 micro-rate-constants; rewritten
// to CL/Q/V1/V2 per the naming guide (Q_MET = k12_met*V1_MET = 0.15*300 = 45;
// V2_MET = Q_MET/k21_met = 45/0.08 = 562.5) -- algebraically identical to the
// original for every input, not a refit (see refactor notes).
KA_MET  = 1.2,   // absorption rate (h-1)
F_MET   = 0.55,  // bioavailability (applied at the dosing-amount level, not
                 // inside the ODE -- matches original: build_doses() pre-
                 // multiplies the dose by F_MET rather than the ODE using it)
V1_MET  = 300,   // central vol (L/70kg)
CL_MET  = 35,    // clearance (L/h)
Q_MET   = 45,    // inter-compartmental CL (L/h), derived = k12_met * V1_MET
V2_MET  = 562.5, // peripheral vol (L), derived = Q_MET / k21_met
MW_MET  = 165.6, // molecular weight (declared, unused in $ODE -- same in original)

// ---- Empagliflozin PK (Macha 2013, Clin Pharmacokinet) ----
// Archetype 1 (no depot, single compartment, linear). KA_EMPA/F_EMPA are
// declared (as in the original) but not referenced in $ODE -- empagliflozin
// is dosed directly into CENT_EMPA with F already folded into the dose
// amount by the R-side build_doses(), same as the original.
KA_EMPA = 0.60,  // absorption (h-1) -- unused in $ODE, kept for fidelity
F_EMPA  = 0.86,  // -- unused in $ODE, kept for fidelity
V1_EMPA = 73,    // (L/70kg)
CL_EMPA = 10.6,

// ---- Semaglutide PK (Lau 2015, Pharm Res) ----
// Archetype 3 variant: depot + central, no peripheral compartment (linear).
KA_SEMA = 0.0085,// SC absorption (h-1), t1/2abs~5d
F_SEMA  = 0.89,  // -- unused in $ODE (also unused in the original own dosing calc), kept for fidelity
V1_SEMA = 12.5,  // (L)
CL_SEMA = 0.053, // very low CL (t1/2~168h)

// ---- Sitagliptin PK (He 2009, Br J Clin Pharmacol) ----
// Archetype 1 (no depot, single compartment, linear).
KA_DPP4 = 0.80,  // -- unused in $ODE, kept for fidelity
F_DPP4  = 0.87,  // -- unused in $ODE, kept for fidelity
V1_DPP4 = 198,
CL_DPP4 = 12.4,

// ---- Glimepiride PK (Massi-Benedetti 1996) ----
// Archetype 1 (no depot, single compartment, linear).
KA_SU   = 0.50,  // -- unused in $ODE, kept for fidelity
F_SU    = 1.00,  // -- unused in $ODE, kept for fidelity
V1_SU   = 12.6,
CL_SU   = 3.1,

// ---- Insulin Degludec PK (Kurtzhals 2011) ----
// Archetype 3 variant: depot + central, no peripheral (linear). NOTE: the
// original never couples this compound concentration/PK to any disease
// PD term (see refactor notes + UPSTREAM_ISSUES.md) -- no EFFECT_INS exists
// because the original has none to rename.
KA_INS  = 0.048, // from SC hexamer depot (h-1) -- a single first-order
                 // absorption step; the original does not model a
                 // multi-hexamer/two-step depot despite the "hexamer" name
V1_INS  = 8.0,
CL_INS  = 0.96,
F_INS   = 0.91,  // -- unused in $ODE, kept for fidelity

// ---- Pioglitazone PK (Eckland 2000) ----
// Archetype 1 (no depot, single compartment, linear).
KA_PIOG = 0.70,  // -- unused in $ODE, kept for fidelity
F_PIOG  = 0.83,  // -- unused in $ODE, kept for fidelity
V1_PIOG = 89,
CL_PIOG = 5.6,

// ---- Metformin PD -- Hill interface (two effect axes, EGP suppression and
// mild Rd stimulation) -- rename only: original hardcoded these as bare
// numeric literals (0.30, 0.12, 800) directly inside the effect formula
// rather than as named $PARAM entries; pulled out here with the same
// values, no numeric change. Shared EC50/gamma, distinct Emax per axis. ----
EMAX_MET_EGP = 0.30, // max fractional EGP suppression (was inline literal)
EMAX_MET_RD  = 0.12, // max fractional Rd increase (was inline literal)
EC50_MET     = 800,  // ng/mL metformin (was inline literal)
GAMMA_MET    = 1,    // original had no explicit Hill coefficient

// ---- Glucose-Insulin Dynamics (Bergman + extensions) [disease, unchanged] ----
SI       = 8e-4, // insulin sensitivity (dL/mU/h, normal ~1e-3)
Sg       = 0.01, // glucose effectiveness (h-1)
EGP0     = 2.4,  // basal EGP (mg/kg/min)
Rd0      = 2.4,  // basal Rd (mg/kg/min)
p2       = 0.05, // remote compartment transfer (h-1)
Vg       = 1.5,  // glucose distribution vol (dL/kg -> scaled by BW)
Vi       = 0.05, // insulin vol (L/kg)

// ---- Insulin secretion [disease, unchanged] ----
beta_sens = 0.6, // GSIS sensitivity (beta-cell function)
beta_M0  = 1.0,  // normalized baseline beta-cell mass
k_prolif = 1e-4, // beta-cell proliferation constant (d-1)
k_apop   = 1e-3, // beta-cell apoptosis rate (d-1) -- elevated in T2DM
Gth      = 90,   // glucose threshold for GSIS (mg/dL)
phi_max  = 900,  // max insulin secretion rate (mU/L/h)

// ---- Glucagon dynamics [disease, unchanged] ----
kout_Gc  = 0.3,  // glucagon turnover (h-1)
Gc_Gp50  = 100,  // Gp at 50% glucagon suppression (mg/dL)
EGP_Gc   = 0.015,// glucagon effect on EGP (mg/kg/min per pg/mL above basal)

// ---- GLP-1 dynamics [disease, unchanged] ----
kin_GLP1 = 15,   // basal GLP-1 production (pmol/L/h)
kout_GLP1= 4.0,  // GLP-1 elimination (h-1)
kDPP4    = 3.5,  // DPP-4 degradation rate (h-1)
GLP1_Gp  = 0.05, // GLP-1 response to meal glucose (pmol/L per mg/dL)

// ---- DPP-4 inhibition (Sitagliptin) -- Hill interface ----
EMAX_DPP4  = 0.80,  // max DPP-4 inhibition (was Imax_dpp4)
EC50_DPP4  = 100,   // ng/mL (was IC50_dpp4*1000 = 0.10*1000, folded into one
                    // named constant -- same number the original computed
                    // at the point of use, not a new value)
GAMMA_DPP4 = 1,     // original had no explicit Hill coefficient

// ---- SGLT2i PD (Empagliflozin) -- Hill interface ----
TmG_base  = 340, // max tubular glucose reabsorption (mg/min/1.73m2) [disease]
RGT0      = 180, // renal glucose threshold (mg/dL) [disease]
EMAX_EMPA  = 0.55,  // max reduction in TmG (was Imax_empa)
EC50_EMPA  = 30000, // ng/mL (was IC50_empa*1000 = 30*1000, folded)
GAMMA_EMPA = 1,     // original had no explicit Hill coefficient
GFR_val   = 100, // individual GFR (mL/min) for UGE calc [disease]

// ---- SU (Glimepiride) PD -- Hill interface ----
EMAX_SU   = 2.5, // fold-increase in insulin secretion (SU max) (was Emax_su)
EC50_SU   = 50,  // ng/mL Glimepiride (was EC50_su)
GAMMA_SU  = 1.5, // Hill coefficient (was n_su)

// ---- GLP-1RA / Semaglutide PD -- Hill interface (two effect axes) ----
EMAX_SEMA_INS  = 0.6, // max increase in GSIS (fold) (was Emax_sema_ins)
EC50_SEMA      = 5,   // nmol/L semaglutide (was EC50_sema, shared by both axes)
GAMMA_SEMA_INS = 1.5, // Hill coefficient on the GSIS axis (was the hardcoded exponent 1.5)
EMAX_SEMA_GC   = 0.40,// max glucagon suppression (was a hardcoded literal
                      // 0.40 inside the E_sema_Gc formula -- pulled out as an
                      // explicit named parameter, same value, not a new claim)
GAMMA_SEMA_GC  = 1,   // original glucagon-suppression term had no explicit
                      // Hill coefficient (implicit gamma=1 ratio)
EMAX_SEMA_WT   = 4.5, // max weight loss (kg) from GLP-1RA at 1yr (was Emax_sema_wt)
                      // -- time-based bespoke term, NOT part of the
                      // concentration-driven Hill interface (see refactor notes)
K_SEMA_WT      = 0.004,// weight loss rate constant (h-1) (was ksema_wt)

// ---- Insulin Resistance dynamics [disease, unchanged] ----
IR_H0    = 1.0,  // hepatic IR index (1=normal, >1=resistant)
IR_P0    = 1.0,  // peripheral IR index
k_IR_FFA = 0.005,// FFA-driven IR increase rate
k_IR_rec = 0.001,// natural IR recovery rate

// ---- FFA / adiposity [disease, unchanged] ----
FFA0     = 0.6,  // baseline plasma FFA (mmol/L)
kFFA_rel = 0.15, // lipolysis rate
kFFA_up  = 0.20, // peripheral FFA uptake/clearance
FFA_Ip50 = 15,   // Ip for 50% lipolysis inhibition (mU/L)
BW_loss_rate = 0.0,// body weight change from lifestyle (kg/day)

// ---- HbA1c kinetics [disease, unchanged] ----
kHbA1c   = 0.0084,// HbA1c equilibration rate (h-1, t1/2~60 days)
HbA1c_ss = 0.165, // HbA1c-to-mean glucose slope (% per mg/dL above 90)

// ---- Renal / complication dynamics [disease, unchanged] ----
k_eGFR_decline = 2e-5,  // eGFR decline rate per unit HbA1c excess (mL/min/day)
k_UACR_rise    = 0.004, // UACR rise rate per unit HbA1c excess
k_eGFR_empa    = 1.5e-5,// nephroprotective effect of SGLT2i on eGFR decline
k_UACR_empa    = 0.003, // SGLT2i -> UACR reduction

// ---- PPARg / Pioglitazone PD -- Hill interface ----
EMAX_PIOG    = 0.35,  // max IR reduction (35%) (was Emax_piog_IR)
EC50_PIOG    = 200,   // ng/mL pioglitazone (was EC50_piog)
GAMMA_PIOG   = 1,     // original had no explicit Hill coefficient
EMAX_PIOG_WT = 3.0,   // max weight GAIN (kg, fluid retention) (was Emax_piog_wt)
                      // -- time-based bespoke term, NOT part of the
                      // concentration-driven Hill interface (see refactor notes)

// ---- Body weight [disease, unchanged] ----
BW_rate  = 0     // net weight change rate (kg/day, baseline)

$CMT
// Drug PK
GUT_MET CENT_MET PERI_MET     // Metformin (gut, central, peripheral)
CENT_EMPA                      // Empagliflozin plasma
GUT_SEMA CENT_SEMA             // Semaglutide SC depot, plasma
CENT_DPP4                      // Sitagliptin plasma
CENT_SU                        // Glimepiride plasma
GUT_INS CENT_INS               // Insulin degludec SC depot, plasma
CENT_PIOG                      // Pioglitazone plasma

// Glucose-Insulin
Gp Gt Ip X_action             // Plasma glucose, tissue glucose, insulin, remote action

// Endocrine
Gc                             // Glucagon (pg/mL)
GLP1                           // Active GLP-1 (pmol/L)
beta_mass                      // beta-cell mass (normalized)

// Intermediate
IR_H IR_P                      // Hepatic and peripheral insulin resistance
FFA                            // Plasma FFA (mmol/L)
BW_t                           // Body weight (kg)

// Endpoints
HbA1c_cmpt                     // HbA1c (%)
eGFR_cmpt                      // eGFR (mL/min/1.73m2)
UACR_cmpt                      // UACR (ug/mg creatinine)

$MAIN
// ---- Compartment initial values ----
// Build-compat fix: the original declared these same values via a separate
// $INIT block, but mrgsolve 2.0.1 treats $CMT + $INIT jointly naming the
// same compartments as a duplicate declaration ("Duplicated model names",
// a hard compile error) -- confirmed on the untouched original. Moved here
// as the standard mrgsolve <cmt>_0 initial-condition idiom instead; same
// values, no numeric change. See UPSTREAM_ISSUES.md.
GUT_MET_0 = 0;
CENT_MET_0 = 0;
PERI_MET_0 = 0;
CENT_EMPA_0 = 0;
GUT_SEMA_0 = 0;
CENT_SEMA_0 = 0;
CENT_DPP4_0 = 0;
CENT_SU_0 = 0;
GUT_INS_0 = 0;
CENT_INS_0 = 0;
CENT_PIOG_0 = 0;
Gp_0 = 180;
Gt_0 = 180;
Ip_0 = 20;
X_action_0 = 0.016;
Gc_0 = 150;
GLP1_0 = 5;
beta_mass_0 = 0.7;         // T2DM: ~50-70% of normal beta-cell mass
IR_H_0 = 2.5;              // elevated hepatic IR
IR_P_0 = 2.0;              // elevated peripheral IR
FFA_0 = 0.8;               // elevated FFA in T2DM
BW_t_0 = 90;               // baseline body weight
HbA1c_cmpt_0 = 9.0;
eGFR_cmpt_0 = 80;
UACR_cmpt_0 = 50;

// ---- Plasma drug concentrations (the exposed C_<STEM> covariate site) ----
double C_MET   = CENT_MET / V1_MET * 1000;   // ng/mL
double C_EMPA  = CENT_EMPA / V1_EMPA * 1000; // ng/mL
double C_SEMA  = CENT_SEMA / V1_SEMA;        // nmol/L (MW~4114 -> approx nmol/L; no *1000, matches original)
double C_DPP4  = CENT_DPP4 / V1_DPP4 * 1000; // ng/mL
double C_SU    = CENT_SU / V1_SU * 1000;     // ng/mL
double C_INS   = CENT_INS / V1_INS * 1000;   // ng/mL (exposed per naming convention; see
                                              // refactor notes -- the original never reads
                                              // this concentration in any disease-PD term)
double C_PIOG  = CENT_PIOG / V1_PIOG * 1000; // ng/mL

// ---- Drug effect functions (Hill interface, EFFECT_<STEM>) ----
// Metformin: AMPK -> reduces EGP (max 30%), mild Rd increase. Rename only --
// original was already a plain Emax*C/(C+EC50) ratio (gamma=1 implicit).
double EFFECT_MET_EGP = USE_MET * EMAX_MET_EGP * C_MET / (C_MET + EC50_MET);
double EFFECT_MET_RD  = USE_MET * EMAX_MET_RD  * C_MET / (C_MET + EC50_MET);

// Empagliflozin: SGLT2 inhibition -> UGE increase. Rename only.
double EFFECT_EMPA = USE_EMPA * EMAX_EMPA * C_EMPA / (C_EMPA + EC50_EMPA);
double TmG    = TmG_base * (1 - EFFECT_EMPA);   // reduced reabsorption capacity
double UGE    = fmax(0.0, (GFR_val * Gp/100 - TmG) / (24*60)); // mg/min -> mg/h per dL
// Urinary glucose excretion (mg/dL plasma equivalent loss per h)
double UGE_loss = fmax(0.0, GFR_val * Gp / 100 - TmG) * 60 / BW_t; // mg/kg/h

// Semaglutide: GLP-1R agonism -- amplifies insulin, suppresses glucagon,
// reduces weight. Two concentration-driven Hill axes (rename only) plus a
// third, time-only weight-loss term the original models as a function of
// elapsed time rather than concentration -- kept exactly as-is, bespoke,
// outside the Hill interface (see refactor notes).
double EFFECT_SEMA_INS = USE_SEMA * EMAX_SEMA_INS * pow(C_SEMA, GAMMA_SEMA_INS)
                          / (pow(EC50_SEMA, GAMMA_SEMA_INS) + pow(C_SEMA, GAMMA_SEMA_INS));
double EFFECT_SEMA_GC  = USE_SEMA * EMAX_SEMA_GC * pow(C_SEMA, GAMMA_SEMA_GC)
                          / (pow(EC50_SEMA, GAMMA_SEMA_GC) + pow(C_SEMA, GAMMA_SEMA_GC));
double BW_SEMA_EFFECT = USE_SEMA * EMAX_SEMA_WT * (1 - exp(-K_SEMA_WT * TIME));

// DPP-4 inhibition: prevents GLP-1 degradation (extends GLP-1 half-life ~2x). Rename only.
double EFFECT_DPP4 = USE_DPP4 * EMAX_DPP4 * C_DPP4 / (C_DPP4 + EC50_DPP4);
double keff_DPP4 = kDPP4 * (1 - EFFECT_DPP4); // NB: computed but not read below --
                                               // the ODE recomputes kDPP4*(1-EFFECT_DPP4)
                                               // directly, exactly as the original did
                                               // with kDPP4*(1-DPP4_inh); preserved as-is
                                               // (harmless pre-existing dead intermediate)

// Sulfonylurea (Glimepiride): increases insulin secretion via SUR1 closure. Rename only.
double EFFECT_SU = USE_SU * EMAX_SU * pow(C_SU, GAMMA_SU) / (pow(C_SU, GAMMA_SU) + pow(EC50_SU, GAMMA_SU));

// Pioglitazone: PPARg -> reduces IR (concentration-driven Hill term, rename
// only) + increases adiponectin/fluid retention (time-only bespoke term,
// same treatment as the semaglutide weight axis above).
double EFFECT_PIOG_IR = USE_PIOG * EMAX_PIOG * C_PIOG / (C_PIOG + EC50_PIOG);
double BW_PIOG_GAIN = USE_PIOG * EMAX_PIOG_WT * (1 - exp(-0.003 * TIME));

// ---- GLP-1 effective level (endogenous + semaglutide contribution) ----
double GLP1_total = GLP1 + USE_SEMA * C_SEMA * 50; // sema in pmol/L equiv

// ---- Insulin secretion (GSIS) ----
double GSIS_base = phi_max * beta_mass * beta_sens
                   * pow(fmax(0.0, Gp - Gth), 1.5)
                   / (pow(fmax(0.0, Gp - Gth), 1.5) + pow(90.0, 1.5));
double GSIS_incretin = GSIS_base * (1 + 0.35 * GLP1_total/(GLP1_total + 10.0)
                                    + EFFECT_SEMA_INS);
double GSIS_su       = GSIS_incretin * (1 + EFFECT_SU);
double dIp_sec       = GSIS_su;  // total secretion into portal

// ---- EGP: suppressed by insulin (via X_action), glucagon drives it up ----
double EGP_ins_supp  = 1.0 / (1.0 + 3.0 * X_action);
double EGP_Gc_drive  = 1.0 + EGP_Gc * fmax(0.0, Gc - Gc0);
double EGP_IR_drive  = IR_H / 2.5; // elevated with hepatic IR
double EGP_val       = EGP0 * EGP_ins_supp * EGP_Gc_drive * EGP_IR_drive
                        * (1 - EFFECT_MET_EGP);
// cap to physiological range
EGP_val = fmax(0.5, fmin(EGP_val, 8.0));

// ---- Rd: insulin-stimulated glucose disposal ----
double SI_eff  = SI / IR_P;      // peripheral IR reduces SI
double Rd_ins  = SI_eff * X_action * Gt;
double Rd_val  = Rd0 + Rd_ins + EFFECT_MET_RD * Rd0;

// ---- Glucagon equation parameters ----
double Gc_ss   = Gc0 * (Gc_Gp50 / fmax(Gp, 50.0)) * (1 - EFFECT_SEMA_GC)
                * (1.0 / (1.0 + 0.02 * fmax(0.0, Ip - Ip0)));

// ---- beta-cell mass dynamics ----
double beta_prolif_rate = k_prolif * beta_mass * fmax(0.0, Gp - 90) / 90.0;
double beta_apop_rate   = k_apop * beta_mass
                          * (1.0 + 0.5 * fmax(0.0, FFA - 0.6))
                          * (1.0 + 0.3 * fmax(0.0, Gp - 180.0)/180.0);
// GLP-1 (incretin) protects beta-cell
double GLP1_beta_prot   = 1.0 - 0.3 * GLP1_total / (GLP1_total + 10.0);
beta_apop_rate *= GLP1_beta_prot;

// ---- FFA: anti-lipolytic effect of insulin ----
double Ip_antilipol  = 1.0 / (1.0 + Ip / FFA_Ip50);
double FFA_release_r = kFFA_rel * Ip_antilipol * BW_t / 90.0;
double FFA_uptake_r  = kFFA_up * FFA;

// ---- IR dynamics ----
// Hepatic IR driven by FFA, glucotoxicity, fat accumulation
double IR_H_drive = k_IR_FFA * fmax(0.0, FFA - 0.5) + k_IR_FFA * fmax(0.0, Gp - 150)/300.0;
double IR_H_recov = k_IR_rec * (IR_H - 1.0);
// Metformin reduces hepatic IR
double IR_H_met   = USE_MET * 0.002 * EFFECT_MET_EGP;
// Pioglitazone reduces both IR_H and IR_P
double dIRH_piog  = EFFECT_PIOG_IR * 0.003 * IR_H;
double dIRP_piog  = EFFECT_PIOG_IR * 0.003 * IR_P;

// ---- Body weight dynamics ----
double BW_UGE_loss   = UGE_loss / 40.0 * 0.001;   // caloric loss from glucosuria (kg/d)
double BW_sema_rate  = USE_SEMA * K_SEMA_WT * (EMAX_SEMA_WT - BW_SEMA_EFFECT) / 24.0;
double BW_piog_rate  = USE_PIOG * 0.003 * (EMAX_PIOG_WT - BW_PIOG_GAIN) / 24.0;
double dBW           = BW_rate - BW_UGE_loss - BW_sema_rate + BW_piog_rate;

// ---- HbA1c equilibration to current mean glucose ----
double Gp_eq         = fmax(Gp, 70.0);
double HbA1c_target  = 5.0 + HbA1c_ss * fmax(0.0, Gp_eq - 90.0);
double dHbA1c        = kHbA1c * (HbA1c_target - HbA1c_cmpt);

// ---- eGFR decline ----
double HbA1c_excess  = fmax(0.0, HbA1c_cmpt - 7.0);
double eGFR_decline  = k_eGFR_decline * HbA1c_excess * 24.0; // per day
double eGFR_protect  = USE_EMPA * k_eGFR_empa * 24.0 * EFFECT_EMPA;
double deGFR         = -(eGFR_decline - eGFR_protect);

// ---- UACR dynamics ----
double UACR_drive    = k_UACR_rise * HbA1c_excess * 24.0;
double UACR_protect  = USE_EMPA * k_UACR_empa * 24.0 * EFFECT_EMPA;
double dUACR         = UACR_drive - UACR_protect;

$ODE
// ===== Drug PK ODEs =====
// Metformin
dxdt_GUT_MET  = -KA_MET * GUT_MET;
dxdt_CENT_MET =  KA_MET*GUT_MET - (CL_MET + Q_MET)/V1_MET * CENT_MET + Q_MET/V2_MET * PERI_MET;
dxdt_PERI_MET =  Q_MET/V1_MET * CENT_MET - Q_MET/V2_MET * PERI_MET;

// Empagliflozin
dxdt_CENT_EMPA  = -CL_EMPA/V1_EMPA * CENT_EMPA;

// Semaglutide
dxdt_GUT_SEMA = -KA_SEMA * GUT_SEMA;
dxdt_CENT_SEMA  = KA_SEMA * GUT_SEMA - CL_SEMA/V1_SEMA * CENT_SEMA;

// Sitagliptin
dxdt_CENT_DPP4 = -CL_DPP4/V1_DPP4 * CENT_DPP4;

// Glimepiride
dxdt_CENT_SU    = -CL_SU/V1_SU * CENT_SU;

// Insulin Degludec
dxdt_GUT_INS  = -KA_INS * GUT_INS;
dxdt_CENT_INS   = KA_INS * GUT_INS - CL_INS/V1_INS * CENT_INS;

// Pioglitazone
dxdt_CENT_PIOG  = -CL_PIOG/V1_PIOG * CENT_PIOG;

// ===== Glucose-Insulin ODEs (Bergman extended) =====
// Plasma glucose (mg/dL)
dxdt_Gp = (EGP_val - Rd_val) / (Vg * BW_t) * 10 - Sg * (Gp - 90) - UGE_loss / (Vg * BW_t);

// Tissue glucose
dxdt_Gt = Sg * (Gp - Gt) - Rd_ins / (Vg * BW_t) * 5;

// Plasma insulin (endogenous, mU/L) -- NB: this term reuses CL_INS/V1_INS
// (renamed from CL_ins/Vc_ins) as the clearance constant for *endogenous*
// plasma insulin, regardless of USE_INS or dosing. This is the original
// own construction (confirmed: Cins/INS_C is never otherwise read by any
// disease equation), preserved exactly under the renamed identifiers --
// see refactor notes and UPSTREAM_ISSUES.md.
dxdt_Ip = (dIp_sec - CL_INS * Ip / V1_INS) / Vi / BW_t * 10;

// Remote insulin action compartment
dxdt_X_action = -p2 * X_action + p2 * SI * Ip;

// ===== Glucagon (pg/mL) =====
dxdt_Gc = kout_Gc * (Gc_ss - Gc);

// ===== Active GLP-1 (pmol/L) =====
dxdt_GLP1 = kin_GLP1 + GLP1_Gp * fmax(0.0, Gp - 90) - (kDPP4 * (1 - EFFECT_DPP4) + kout_GLP1) * GLP1;

// ===== beta-cell mass (normalized, 0-1) =====
dxdt_beta_mass = (beta_prolif_rate - beta_apop_rate) / (24.0 * 365.25);  // per year

// ===== Insulin Resistance indices =====
dxdt_IR_H = IR_H_drive - IR_H_recov - IR_H_met - dIRH_piog;
dxdt_IR_P = k_IR_FFA * fmax(0.0, FFA - 0.5)
             - k_IR_rec * (IR_P - 1.0) - dIRP_piog
             + 0.001 * fmax(0.0, Gp - 150)/300.0;

// ===== FFA (mmol/L) =====
dxdt_FFA = FFA_release_r - FFA_uptake_r;

// ===== Body weight (kg) =====
dxdt_BW_t = dBW;

// ===== Endpoints =====
dxdt_HbA1c_cmpt = dHbA1c;
dxdt_eGFR_cmpt  = deGFR;
dxdt_UACR_cmpt  = dUACR;

$CAPTURE
C_MET C_EMPA C_SEMA C_DPP4 C_SU C_INS C_PIOG
EGP_val Rd_val UGE_loss
EFFECT_MET_EGP EFFECT_MET_RD EFFECT_EMPA EFFECT_SEMA_INS EFFECT_SEMA_GC EFFECT_DPP4 EFFECT_SU EFFECT_PIOG_IR
beta_prolif_rate beta_apop_rate
'

# ============================================================================
# COMPILE MODEL
# ============================================================================
mod <- mcode("T2DM_QSP_refactored", t2dm_model_code)

# ============================================================================
# DOSE EVENTS BUILDER (compartment names renamed to match the naming
# convention; dose amounts, intervals, and durations unchanged)
# ============================================================================
build_doses <- function(scenario, start_h = 0, days = 365) {
  events <- list()
  total_h <- days * 24

  if (scenario %in% c("metformin","combo_empa","combo_sema","triple","insulin","dpp4i")) {
    # Metformin 1000 mg BID (F=55% -> 550 mg reaches gut cmpt)
    # Represent as mg bioavailable entering GUT_MET
    doses_met <- ev(amt = 1000 * 0.55, cmt = "GUT_MET", ii = 12, addl = 2*days - 1, time = start_h)
    events <- c(events, list(doses_met))
  }

  if (scenario %in% c("combo_empa","triple")) {
    # Empagliflozin 10 mg QD (F=86%)
    doses_empa <- ev(amt = 10 * 0.86 * 1e6 / 450.9, cmt = "CENT_EMPA",
                     ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_empa))
  }

  if (scenario %in% c("combo_sema","triple")) {
    # Semaglutide 1 mg SC weekly
    doses_sema <- ev(amt = 1000 / 4113.6 * 1000, cmt = "GUT_SEMA",
                     ii = 168, addl = ceiling(days/7) - 1, time = start_h)
    events <- c(events, list(doses_sema))
  }

  if (scenario == "dpp4i") {
    # Sitagliptin 100 mg QD
    doses_dpp4 <- ev(amt = 100 * 0.87 * 1e6 / 407.5, cmt = "CENT_DPP4",
                     ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_dpp4))
  }

  if (scenario == "insulin") {
    # Insulin degludec 20 U QD (1 U = ~0.0347 mg)
    doses_ins <- ev(amt = 20 * 0.0347 * 0.91 * 1000 / 6103, cmt = "GUT_INS",
                    ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_ins))
  }

  if (length(events) == 0) return(ev(amt = 0, cmt = "GUT_MET", time = 0))
  Reduce(c, events)
}

# ============================================================================
# SIMULATION FUNCTION
# ============================================================================
run_scenario <- function(scenario_name, use_flags, days = 365) {
  params_update <- list(
    USE_MET  = use_flags[["met"]],
    USE_EMPA = use_flags[["empa"]],
    USE_SEMA = use_flags[["sema"]],
    USE_DPP4 = use_flags[["dpp4"]],
    USE_SU   = use_flags[["su"]],
    USE_INS  = use_flags[["ins"]],
    USE_PIOG = use_flags[["piog"]]
  )

  dose_ev <- build_doses(scenario_name, days = days)
  sim_times <- seq(0, days * 24, by = 6)  # every 6 hours

  out <- mod %>%
    param(params_update) %>%
    ev(dose_ev) %>%
    mrgsim(end = days * 24, delta = 6) %>%
    as.data.frame()

  out$scenario <- scenario_name
  out
}

# ============================================================================
# SCENARIO DEFINITIONS (7 scenarios, unchanged)
# ============================================================================
scenarios <- list(
  list(
    name  = "No treatment",
    key   = "no_tx",
    flags = list(met=0, empa=0, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin 1000 BID",
    key   = "metformin",
    flags = list(met=1, empa=0, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Empagliflozin",
    key   = "combo_empa",
    flags = list(met=1, empa=1, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Semaglutide",
    key   = "combo_sema",
    flags = list(met=1, empa=0, sema=1, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Triple (Met+Empa+Sema)",
    key   = "triple",
    flags = list(met=1, empa=1, sema=1, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Insulin Degludec",
    key   = "insulin",
    flags = list(met=1, empa=0, sema=0, dpp4=0, su=0, ins=1, piog=0)
  ),
  list(
    name  = "Metformin + Sitagliptin",
    key   = "dpp4i",
    flags = list(met=1, empa=0, sema=0, dpp4=1, su=0, ins=0, piog=0)
  )
)

# ============================================================================
# RUN ALL SCENARIOS
# ============================================================================
cat("Running T2DM QSP simulations (refactored)...\n")
sim_results <- map_dfr(scenarios, function(s) {
  cat(sprintf("  [%s] %s\n", s$key, s$name))
  run_scenario(s$key, s$flags, days = 365)
})

# Convert time to days
sim_results <- sim_results %>%
  mutate(time_days = time / 24,
         scenario = factor(scenario, levels = map_chr(scenarios, "key"),
                          labels = map_chr(scenarios, "name")))

cat("Simulation complete. N rows:", nrow(sim_results), "\n")

# ============================================================================
# SUMMARY TABLE AT 52 WEEKS
# ============================================================================
summary_52w <- sim_results %>%
  filter(abs(time_days - 365) < 1) %>%
  group_by(scenario) %>%
  slice_tail(n = 1) %>%
  summarise(
    HbA1c_52w     = round(mean(HbA1c_cmpt), 2),
    dHbA1c        = round(mean(HbA1c_cmpt) - 9.0, 2),
    Gp_52w        = round(mean(Gp), 1),
    BW_52w        = round(mean(BW_t), 1),
    dBW           = round(mean(BW_t) - 90.0, 1),
    eGFR_52w      = round(mean(eGFR_cmpt), 1),
    UACR_52w      = round(mean(UACR_cmpt), 0),
    beta_mass_52w = round(mean(beta_mass), 3),
    .groups = "drop"
  )

cat("\n=== 52-Week Summary ===\n")
print(summary_52w, n = 10)

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

# Color palette for 7 scenarios
scen_colors <- c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3",
  "#ff7f00","#a65628","#f781bf"
)

# Plot 1: HbA1c trajectories
p_hba1c <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                  aes(x = time_days, y = HbA1c_cmpt, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 7.0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scen_colors) +
  labs(title = "HbA1c Trajectory Over 1 Year",
       x = "Time (days)", y = "HbA1c (%)", color = "Treatment") +
  theme_bw() + theme(legend.position = "right")

# Plot 2: Plasma glucose
p_glucose <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                    aes(x = time_days, y = Gp, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 126, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Plasma Glucose Over Time",
       x = "Time (days)", y = "Plasma Glucose (mg/dL)", color = "Treatment") +
  theme_bw()

# Plot 3: Body weight
p_bw <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
               aes(x = time_days, y = BW_t, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Body Weight Over Time",
       x = "Time (days)", y = "Body Weight (kg)", color = "Treatment") +
  theme_bw()

# Plot 4: β-cell mass
p_beta <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                 aes(x = time_days, y = beta_mass, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "β-cell Mass (Normalized)",
       x = "Time (days)", y = "β-cell Mass (relative)", color = "Treatment") +
  theme_bw()

# Plot 5: eGFR
p_egfr <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                 aes(x = time_days, y = eGFR_cmpt, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
  scale_color_manual(values = scen_colors) +
  labs(title = "eGFR Over Time (Renal Protection)",
       x = "Time (days)", y = "eGFR (mL/min/1.73m²)", color = "Treatment") +
  theme_bw()

# Plot 6: Urinary Glucose Excretion (SGLT2i effect)
p_uge <- ggplot(sim_results %>%
                  filter(scenario %in% c("Metformin + Empagliflozin","Triple (Met+Empa+Sema)"),
                         time_days %% 1 < 0.3),
                aes(x = time_days, y = UGE_loss * 24, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors[c(3,5)]) +
  labs(title = "Urinary Glucose Excretion (SGLT2i)",
       x = "Time (days)", y = "UGE (mg/kg/day)", color = "Treatment") +
  theme_bw()

# Plot 7: GLP-1 with and without DPP-4 inhibition
p_glp1 <- ggplot(sim_results %>%
                   filter(scenario %in% c("No treatment","Metformin + Sitagliptin",
                                          "Metformin + Semaglutide"),
                          time_days < 7),
                 aes(x = time_days, y = GLP1, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors[c(1,7,4)]) +
  labs(title = "Active GLP-1 Levels (First Week)",
       x = "Time (days)", y = "Active GLP-1 (pmol/L)", color = "Treatment") +
  theme_bw()

# Summary bar chart
p_summary <- summary_52w %>%
  ggplot(aes(x = reorder(scenario, dHbA1c), y = dHbA1c, fill = scenario)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", dHbA1c)), hjust = 1.1, color = "white") +
  scale_fill_manual(values = scen_colors) +
  coord_flip() +
  labs(title = "ΔHbA1c at 52 Weeks",
       x = NULL, y = "ΔHbA1c (%)") +
  theme_bw()

# Display plots
print(p_hba1c)
print(p_glucose)
print(p_bw)
print(p_beta)
print(p_egfr)
print(summary_52w)

# ============================================================================
# CLINICAL CALIBRATION CHECKPOINTS
# ============================================================================
cat("\n=== Clinical Calibration Validation ===\n")
calib_checks <- tibble::tribble(
  ~Drug,             ~Trial,           ~Endpoint,    ~Observed,  ~Simulated,
  "Metformin",       "UKPDS-34",       "ΔHbA1c(%)",  "-1.4",
    as.character(round(filter(summary_52w, grepl("Metformin 1000", scenario))$dHbA1c[1], 1)),
  "Empagliflozin",   "EMPA-REG",       "ΔHbA1c(%)",  "-0.7",
    as.character(round(filter(summary_52w, grepl("Empa", scenario) & !grepl("Sema", scenario))$dHbA1c[1], 1)),
  "Semaglutide",     "SUSTAIN-6",      "ΔHbA1c(%)",  "-1.4",
    as.character(round(filter(summary_52w, grepl("Sema", scenario) & !grepl("Empa", scenario))$dHbA1c[1], 1)),
  "Semaglutide",     "SUSTAIN-6",      "ΔBW(kg)",    "-4.5",
    as.character(round(filter(summary_52w, grepl("Sema", scenario) & !grepl("Empa", scenario))$dBW[1], 1)),
  "Sitagliptin",     "TECOS",          "ΔHbA1c(%)",  "-0.7",
    as.character(round(filter(summary_52w, grepl("Sitagliptin", scenario))$dHbA1c[1], 1))
)
print(calib_checks)

cat("\nModel ready (refactored). Use Shiny app for interactive exploration.\n")
