## ============================================================================
## Abdominal Aortic Aneurysm (AAA) QSP Model -- mrgsolve
## REFACTORED sibling of aaa_mrgsolve_model.R -- pluggable-PK naming convention
## applied to Doxycycline (DOXY), Propranolol/beta-blocker (BB), and Statin
## (STAT). See aaa_refactor_notes.md for the full account, especially the
## Doxycycline TISSUE-vs-plasma finding: C_DOXY exposes the AORTIC WALL
## tissue concentration, not plasma, because that is what the original's own
## MMP-9/MMP-2 effect terms actually read. All numeric parameter values are
## copied verbatim from the original; only PK structure/naming and the
## Hill-interface expression are reorganized. The original is untouched.
## ============================================================================
## Description:
##   Quantitative Systems Pharmacology model for AAA, integrating:
##   1. Three-drug PK: Doxycycline (MMP inhibitor), Statin (anti-inflammatory,
##      antioxidant), and Propranolol (β-blocker, hemodynamic control)
##   2. Disease PD: macrophage polarization, MMP-2/9 activity, VSMC apoptosis,
##      elastin degradation, collagen remodeling, intraluminal thrombus (ILT),
##      and aortic diameter progression
##
## Key References for Parameter Calibration:
##   - Mosorin M et al. (2001) Doxycycline treatment reduces AAA expansion
##     (Annals of Surgery; pilot trial, n=36)
##   - PHAST trial: Lindeman JH et al. (2009) doxycycline 100mg/day
##     significantly reduced aortic elastin degradation (Circulation)
##   - UK Small Aneurysm Trial (1998) propranolol arm; Walker et al. (2012)
##     retrospective analysis showing beta-blockers reduce AAA growth
##   - Brady AR et al. (2004) and Meijer CA et al. (2012) statin use and
##     AAA growth/rupture risk reduction
##   - Sakalihasan N et al. (2005) MMP-9 as biomarker of AAA activity
##     (Nat Rev Cardiol)
##   - Sweeting MJ et al. (2012) meta-analysis of AAA growth rates (BJS)
## ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================================
## MODEL DEFINITION
## ============================================================================
code <- '
$PLUGIN Rcpp

$PARAM @annotated
// -- Doxycycline (DOXY) PK Parameters -- archetype 3 (depot+central+peripheral)
//    PLUS a 4th, genuinely distinct tissue compartment (TISSUE_DOXY): the
//    disease PD below reads the AORTIC WALL TISSUE concentration, NOT
//    plasma. See refactor notes for the tissue-vs-plasma finding. --------
KA_DOXY   :  0.50  : Doxycycline absorption rate constant (h-1)
F_DOXY    :  0.93  : Doxycycline oral bioavailability (fraction)
CL_DOXY   :  3.90  : Doxycycline clearance (L/h); Agwuh & MacGowan 2006 JAC
V1_DOXY   :  110.0 : Doxycycline central volume (L); ~ 1.58 L/kg * 70 kg
Q_DOXY    :  15.0  : Doxycycline inter-compartment clearance (L/h)
V2_DOXY   :  300.0 : Doxycycline peripheral volume (L)
KTA_DOXY  :  0.05  : Doxycycline aortic tissue distribution rate (h-1)
KTE_DOXY  :  0.04  : Doxycycline aortic tissue elimination rate (h-1)

// ── Statin (Simvastatin) PK Parameters ──────────────────────────────────────
KA_STAT   :  0.50  : Statin absorption rate constant (h-1)
F_STAT    :  0.05  : Statin oral bioavailability (first-pass ~95%); Pharmacokin.
CL_STAT   :  36.0  : Statin clearance (L/h) via CYP3A4
V1_STAT   :  245.0 : Statin central volume (L); ~ 3.5 L/kg
Q_STAT    :  8.0   : Statin inter-compartment clearance (L/h)
V2_STAT   :  120.0 : Statin peripheral volume (L)

// ── Propranolol (β-blocker) PK Parameters ───────────────────────────────────
KA_BB  :  0.90  : Propranolol absorption rate constant (h-1)
F_BB   :  0.30  : Propranolol oral bioavailability (~25-35%)
CL_BB  :  50.0  : Propranolol clearance (L/h) via hepatic metabolism
V1_BB  :  300.0 : Propranolol central volume (L); ~ 4.3 L/kg
Q_BB   :  10.0  : Propranolol inter-compartment clearance (L/h)
V2_BB  :  200.0 : Propranolol peripheral volume (L)

// ── Disease PD — Macrophage/Inflammatory Parameters ─────────────────────────
MAC0       : 100.0  : Baseline macrophage activity (arbitrary units)
k_MAC_in   :  0.10  : Macrophage recruitment/activation rate (day-1)
k_MAC_out  :  0.05  : Macrophage elimination rate (day-1)
TNF0       :  10.0  : Baseline TNF-alpha (pg/mL)
k_TNF_syn  :  1.00  : TNF synthesis rate (pg/mL/day)
k_TNF_deg  :  0.50  : TNF degradation rate (day-1)
ROS0       :  1.00  : Baseline ROS level (nmol/mg protein)
k_ROS_syn  :  0.20  : ROS production rate by macrophages
k_ROS_deg  :  0.15  : ROS elimination rate (day-1)

// ── Disease PD — MMP Axis Parameters ────────────────────────────────────────
MMP9_0     :  1.00  : Baseline active MMP-9 (ng/mL; normal ~1 ng/mL)
k_MMP9_syn :  0.50  : MMP-9 synthesis rate by macrophages/neutrophils (day-1)
k_MMP9_deg :  0.30  : MMP-9 degradation/inhibition by TIMP-1 (day-1)
EMAX_DOXY_MMP9    :  0.85  : Doxycycline max inhibition of MMP-9 (Emax)
EC50_DOXY_MMP9    :  0.20  : Doxy AORTIC TISSUE conc. for half-max MMP-9 inhibition (mg/L)
GAMMA_DOXY_MMP9   :  1.0   : Doxycycline MMP-9 Hill coefficient (none in original; =1)
EMAX_STAT_MMP9    :  0.40  : Statin max inhibition of MMP-9 via pleiotropic effect
EC50_STAT_MMP9    :  0.05  : Statin plasma conc. for half-max MMP-9 inhibition (mg/L)
GAMMA_STAT_MMP9   :  1.0   : Statin MMP-9 Hill coefficient (none in original; =1)
MMP2_0     :  1.00  : Baseline active MMP-2 (ng/mL)
k_MMP2_syn :  0.30  : MMP-2 synthesis rate (day-1)
k_MMP2_deg :  0.25  : MMP-2 elimination rate (day-1)
EMAX_DOXY_MMP2    :  0.70  : Doxycycline max inhibition of MMP-2
EC50_DOXY_MMP2    :  0.25  : Doxy AORTIC TISSUE conc. for half-max MMP-2 inhibition (mg/L)
GAMMA_DOXY_MMP2   :  1.0   : Doxycycline MMP-2 Hill coefficient (none in original; =1)

// ── Disease PD — ECM / Structural Parameters ────────────────────────────────
ELAST0     : 100.0  : Baseline elastin content (% of normal)
k_ELAST_deg:  0.02  : Baseline elastin degradation rate (day-1)
k_ELAST_syn:  0.001 : Elastin synthesis rate (very slow in adults; day-1)
COL0       : 100.0  : Baseline total collagen content (% of normal)
k_COL_deg  :  0.015 : Baseline collagen degradation rate (day-1)
k_COL_syn  :  0.025 : Collagen synthesis rate (day-1)
VSMC0      : 100.0  : Baseline VSMC density (% of normal)
k_VSMC_apop:  0.008 : VSMC apoptosis rate (day-1; driven by MMP/TNF/ROS)
k_VSMC_prol:  0.005 : VSMC proliferation rate (day-1)

// ── Disease PD — ILT and Diameter Parameters ────────────────────────────────
ILT0       :  0.0   : Initial ILT volume (mL; 0 = no initial thrombus)
k_ILT_grow :  0.002 : ILT growth rate per unit flow turbulence (mL/day)
k_ILT_lyse :  0.001 : ILT fibrinolysis rate (day-1)
DIAM0      : 30.0   : Initial aortic diameter (mm; normal ~20mm but AAA ~30+)
k_diam_grow:  0.005 : Baseline aortic expansion rate (mm/day; ~ 2 mm/yr)
ECM_weight :  0.60  : ECM degradation contribution to diameter expansion
VSMC_weight:  0.20  : VSMC loss contribution to diameter expansion
ILT_weight :  0.10  : ILT contribution to diameter expansion (luminal)

// -- Hemodynamic PD -- Beta-blocker (BB) Effect on Wall Stress ------------
SBP0       : 145.0  : Baseline systolic BP (mmHg; typical AAA patient)
EMAX_BB    :  0.20  : Max BP reduction by beta-blocker (fraction)
EC50_BB    :  0.02  : Propranolol plasma conc. for half-max BP effect (mg/L)
GAMMA_BB   :  1.0   : Propranolol BP-effect Hill coefficient (none in original; =1)
HR0        :  75.0  : Baseline heart rate (bpm)
// EMAX_BB_HR / EC50_BB_HR (was Imax_bb_hr / IC50_bb_hr): declared in the
// ORIGINAL file but never referenced by any $ODE/$TABLE expression -- a
// dead, unused parameter pair (same treatment as the TCZ RTOT_TCZ
// precedent in other refactor notes). Preserved as-is. See notes.
EMAX_BB_HR :  0.30  : Max HR reduction by beta-blocker (UNUSED in model body)
EC50_BB_HR :  0.015 : Propranolol conc. for half-max HR reduction (mg/L) (UNUSED)

// -- ACE Inhibitor (Perindopril) PK/PD -- untouched: not one of the 3
//    target compounds for this refactor (Doxycycline/Propranolol/Statin) --
//    (simplified: ACEi modeled as additive BP reduction effect)
ACEi_dose  :  0.0   : ACE inhibitor active plasma level (mg/L; 0 = no ACEi)
Imax_acei  :  0.18  : Max BP reduction by ACEi
IC50_acei  :  0.005 : ACEi conc. for half-max effect

// -- Dosing Flags -- declared in the ORIGINAL file but never referenced by
//    any $ODE/$TABLE expression (on/off is controlled entirely by which
//    `events` a scenario attaches, not by these flags). Preserved as-is
//    (unused). See refactor notes. ----------------------------------------
DOXY_ON    :  0.0   : Doxycycline dosing switch (0=off, 1=on) (UNUSED)
STAT_ON    :  0.0   : Statin dosing switch (UNUSED)
PROP_ON    :  0.0   : Propranolol dosing switch (UNUSED)

// -- Statin (STAT) pleiotropic effects on TNF/ROS -- promoted from magic
//    numbers HARDCODED inline in the original’s $ODE body (0.30/0.05 and
//    0.35/0.08) into named Hill parameters, per the naming convention.
//    Same values, same ratio shape -- a promotion, not a refit. -----------
EMAX_STAT_NFKB  :  0.30  : Statin max NF-kB-mediated suppression of TNF synthesis (was hardcoded 0.30)
EC50_STAT_NFKB  :  0.05  : Statin conc. for half-max NF-kB suppression (mg/L) (was hardcoded 0.05)
GAMMA_STAT_NFKB :  1.0   : Statin NF-kB Hill coefficient (none in original; =1)
EMAX_STAT_ROS   :  0.35  : Statin max antioxidant suppression of ROS production (was hardcoded 0.35)
EC50_STAT_ROS   :  0.08  : Statin conc. for half-max ROS suppression (mg/L) (was hardcoded 0.08)
GAMMA_STAT_ROS  :  1.0   : Statin ROS Hill coefficient (none in original; =1)

$CMT @annotated
// -- Doxycycline PK (archetype 3 + a 4th, distinct tissue compartment) ----
GUT_DOXY   : Doxycycline gut compartment (mg)
CENT_DOXY  : Doxycycline central (plasma, mg)
PERI_DOXY: Doxycycline peripheral tissue (mg)
TISSUE_DOXY : Doxycycline AORTIC WALL tissue compartment (mg/L-equivalent; this, NOT plasma, is what C_DOXY exposes and what the MMP-9/MMP-2 effect terms read)

// ── Statin PK ─────────────────────────────────────────────────────────────
GUT_STAT   : Statin gut compartment (mg)
CENT_STAT  : Statin central plasma (mg)
PERI_STAT: Statin peripheral (mg)

// ── Propranolol PK ─────────────────────────────────────────────────────────
GUT_BB   : Propranolol gut compartment (mg)
CENT_BB  : Propranolol central plasma (mg)
PERI_BB: Propranolol peripheral (mg)

// ── Disease PD ─────────────────────────────────────────────────────────────
MAC    : Macrophage activity index (AU)
TNF    : TNF-alpha plasma level (pg/mL)
ROSO   : Reactive oxygen species (nmol/mg protein)
MMP9   : Active MMP-9 (ng/mL)
MMP2   : Active MMP-2 (ng/mL)
ELAST  : Aortic wall elastin content (% baseline)
COLLAG : Aortic wall collagen content (% baseline)
VSMC   : VSMC density (% baseline)
ILT    : Intraluminal thrombus volume (mL)
DIAM   : Maximal aortic diameter (mm)

$MAIN
// Initial conditions
GUT_DOXY_0    = 0;
CENT_DOXY_0   = 0;
PERI_DOXY_0 = 0;
TISSUE_DOXY_0  = 0;
GUT_STAT_0    = 0;
CENT_STAT_0   = 0;
PERI_STAT_0 = 0;
GUT_BB_0    = 0;
CENT_BB_0   = 0;
PERI_BB_0 = 0;
MAC_0     = MAC0;
TNF_0     = TNF0;
ROSO_0    = ROS0;
MMP9_0    = MMP9_0;
MMP2_0    = MMP2_0;
ELAST_0   = ELAST0;
COLLAG_0  = COL0;
VSMC_0    = VSMC0;
ILT_0     = ILT0;
DIAM_0    = DIAM0;

$ODE
// ============================================================================
// SECTION 1: DOXYCYCLINE PK (2-compartment + tissue)
// ============================================================================
double ke_DOXY  = CL_DOXY / V1_DOXY;
double k12_DOXY = Q_DOXY  / V1_DOXY;
double k21_DOXY = Q_DOXY  / V2_DOXY;

dxdt_GUT_DOXY    = -KA_DOXY * GUT_DOXY;
dxdt_CENT_DOXY   =  KA_DOXY * F_DOXY * GUT_DOXY - (ke_DOXY + k12_DOXY) * CENT_DOXY + k21_DOXY * PERI_DOXY;
dxdt_PERI_DOXY =  k12_DOXY * CENT_DOXY - k21_DOXY * PERI_DOXY;
dxdt_TISSUE_DOXY  =  KTA_DOXY * (CENT_DOXY / V1_DOXY) - KTE_DOXY * TISSUE_DOXY;  // aortic tissue conc.

double Cp_doxy_plasma = CENT_DOXY / V1_DOXY;  // plasma doxycycline (mg/L) -- informational only, NOT read by any PD term
double C_DOXY = TISSUE_DOXY;  // *** the exposed, PD-facing concentration: AORTIC TISSUE level, NOT plasma ***

// ============================================================================
// SECTION 2: STATIN PK (2-compartment with first-pass)
// ============================================================================
double ke_STAT  = CL_STAT / V1_STAT;
double k12_STAT = Q_STAT  / V1_STAT;
double k21_STAT = Q_STAT  / V2_STAT;

dxdt_GUT_STAT    = -KA_STAT * GUT_STAT;
dxdt_CENT_STAT   =  KA_STAT * F_STAT * GUT_STAT - (ke_STAT + k12_STAT) * CENT_STAT + k21_STAT * PERI_STAT;
dxdt_PERI_STAT =  k12_STAT * CENT_STAT - k21_STAT * PERI_STAT;

double C_STAT = CENT_STAT / V1_STAT;       // exposed, PD-facing: plasma statin (mg/L)

// ============================================================================
// SECTION 3: PROPRANOLOL PK (2-compartment)
// ============================================================================
double ke_BB  = CL_BB / V1_BB;
double k12_BB = Q_BB  / V1_BB;
double k21_BB = Q_BB  / V2_BB;

dxdt_GUT_BB    = -KA_BB * GUT_BB;
dxdt_CENT_BB   =  KA_BB * F_BB * GUT_BB - (ke_BB + k12_BB) * CENT_BB + k21_BB * PERI_BB;
dxdt_PERI_BB =  k12_BB * CENT_BB - k21_BB * PERI_BB;

double C_BB = CENT_BB / V1_BB;      // exposed, PD-facing: plasma propranolol (mg/L)

// ============================================================================
// SECTION 4: HEMODYNAMIC PD
// ============================================================================
// Propranolol reduces SBP (β1-blockade → ↓CO → ↓BP)
double EFFECT_BB = EMAX_BB * pow(C_BB, GAMMA_BB) / (pow(EC50_BB, GAMMA_BB) + pow(C_BB, GAMMA_BB));
double SBP_red_acei = Imax_acei  * ACEi_dose / (IC50_acei + ACEi_dose);
double SBP_curr     = SBP0 * (1.0 - EFFECT_BB - SBP_red_acei);
if (SBP_curr < 90.0) SBP_curr = 90.0;  // floor

// Wall stress ∝ P*r / (2*h)  — simplified as proportional to SBP and DIAM
double WALL_STR = SBP_curr * DIAM / (DIAM0 * 145.0);  // normalized

// ============================================================================
// SECTION 5: MACROPHAGE & INFLAMMATORY DYNAMICS
// ============================================================================
// Macrophage activation driven by wall stress and ROS
double MAC_stim = WALL_STR * ROSO / ROS0;
dxdt_MAC = k_MAC_in * MAC0 * MAC_stim - k_MAC_out * MAC;

// TNF-alpha driven by macrophages; statin pleiotropic suppression
double EFFECT_STAT_NFKB = EMAX_STAT_NFKB * pow(C_STAT, GAMMA_STAT_NFKB) / (pow(EC50_STAT_NFKB, GAMMA_STAT_NFKB) + pow(C_STAT, GAMMA_STAT_NFKB));
dxdt_TNF = k_TNF_syn * (MAC / MAC0) * (1.0 - EFFECT_STAT_NFKB) - k_TNF_deg * TNF;

// ROS driven by macrophage NADPH oxidase; statin antioxidant
double EFFECT_STAT_ROS = EMAX_STAT_ROS * pow(C_STAT, GAMMA_STAT_ROS) / (pow(EC50_STAT_ROS, GAMMA_STAT_ROS) + pow(C_STAT, GAMMA_STAT_ROS));
dxdt_ROSO = k_ROS_syn * (MAC / MAC0) - k_ROS_deg * ROSO * (1.0 + EFFECT_STAT_ROS);

// ============================================================================
// SECTION 6: MMP DYNAMICS (key targets for doxycycline)
// ============================================================================
// MMP-9 driven by macrophages/neutrophils, inhibited by doxycycline (via
// AORTIC TISSUE concentration C_DOXY, not plasma) and statin (via plasma C_STAT)
double EFFECT_DOXY_MMP9 = EMAX_DOXY_MMP9 * pow(C_DOXY, GAMMA_DOXY_MMP9) / (pow(EC50_DOXY_MMP9, GAMMA_DOXY_MMP9) + pow(C_DOXY, GAMMA_DOXY_MMP9));
double EFFECT_STAT_MMP9 = EMAX_STAT_MMP9 * pow(C_STAT, GAMMA_STAT_MMP9) / (pow(EC50_STAT_MMP9, GAMMA_STAT_MMP9) + pow(C_STAT, GAMMA_STAT_MMP9));
double InhMMP9_total = 1.0 - (1.0 - EFFECT_DOXY_MMP9) * (1.0 - EFFECT_STAT_MMP9);
dxdt_MMP9 = k_MMP9_syn * (MAC / MAC0) * (1.0 - InhMMP9_total) - k_MMP9_deg * MMP9;

// MMP-2 driven by VSMCs and fibroblasts, partially inhibited by doxycycline
// (again via AORTIC TISSUE concentration C_DOXY, not plasma)
double EFFECT_DOXY_MMP2 = EMAX_DOXY_MMP2 * pow(C_DOXY, GAMMA_DOXY_MMP2) / (pow(EC50_DOXY_MMP2, GAMMA_DOXY_MMP2) + pow(C_DOXY, GAMMA_DOXY_MMP2));
dxdt_MMP2 = k_MMP2_syn * (1.0 - EFFECT_DOXY_MMP2) - k_MMP2_deg * MMP2;

// ============================================================================
// SECTION 7: VSMC DYNAMICS
// ============================================================================
// VSMC apoptosis driven by MMP-9 (anoikis), TNF-alpha, and ROS
double VSMC_apop_rate = k_VSMC_apop * (MMP9/MMP9_0) * (TNF/TNF0) * (ROSO/ROS0);
double VSMC_prol_rate = k_VSMC_prol * (VSMC / VSMC0);
dxdt_VSMC = VSMC_prol_rate - VSMC_apop_rate * VSMC;
if (VSMC < 0.1 * VSMC0) dxdt_VSMC = 0;  // floor at 10% of baseline

// ============================================================================
// SECTION 8: ECM REMODELING
// ============================================================================
// Elastin degradation by MMP-9, MMP-12 (macrophage elastase); very slow synthesis
double elast_deg_rate = k_ELAST_deg * (MMP9/MMP9_0) * (MMP2/MMP2_0);
dxdt_ELAST = k_ELAST_syn * ELAST0 - elast_deg_rate * ELAST;
if (ELAST < 0) dxdt_ELAST = 0;

// Collagen degradation by MMP-1/8/13, synthesis by VSMCs and TGF-β
double col_deg_rate = k_COL_deg * (MMP9/MMP9_0);
double col_syn_rate = k_COL_syn * (VSMC / VSMC0);
dxdt_COLLAG = col_syn_rate * COL0 - col_deg_rate * COLLAG;
if (COLLAG < 5.0) dxdt_COLLAG = 0;

// ============================================================================
// SECTION 9: INTRALUMINAL THROMBUS (ILT)
// ============================================================================
// ILT grows with turbulence (∝ DIAM) and platelet activation
double turbulence_factor = pow(DIAM / DIAM0, 2.0);  // turbulence ↑ with diam^2
dxdt_ILT = k_ILT_grow * turbulence_factor - k_ILT_lyse * ILT;
if (ILT < 0) dxdt_ILT = 0;

// ============================================================================
// SECTION 10: AORTIC DIAMETER DYNAMICS (primary endpoint)
// ============================================================================
// Expansion rate driven by ECM degradation, VSMC loss, and ILT contribution
double ECM_loss_idx  = (1.0 - ELAST/ELAST0) * ECM_weight;
double VSMC_loss_idx = (1.0 - VSMC/VSMC0)   * VSMC_weight;
double ILT_contrib   = (ILT / 50.0)          * ILT_weight;  // normalize to 50mL
double expansion_driver = 1.0 + ECM_loss_idx + VSMC_loss_idx + ILT_contrib;

// Wall stress accelerates expansion beyond a threshold
double stress_accel = (WALL_STR > 1.2) ? (WALL_STR - 1.2) * 0.5 : 0.0;

dxdt_DIAM = k_diam_grow * expansion_driver * (1.0 + stress_accel);
if (DIAM > 120.0) dxdt_DIAM = 0;  // biological ceiling

$TABLE
double AUC_MMP9  = MMP9;
double AUC_ELAST = ELAST;
double AUC_DIAM  = DIAM;
double Rupture_P = 1.0 / (1.0 + exp(-(0.12 * (DIAM - 55.0))));  // logistic risk
double MMP9_reduction = (MMP9_0 - MMP9) / MMP9_0 * 100.0;
double ELAST_preserved = ELAST;
double DIAM_mm   = DIAM;
double SBP_mmHg  = SBP_curr;
double Wall_Stress_idx = WALL_STR;

$CAPTURE Cp_doxy_plasma C_DOXY C_STAT C_BB EFFECT_DOXY_MMP9 EFFECT_DOXY_MMP2
         EFFECT_STAT_MMP9 EFFECT_STAT_NFKB EFFECT_STAT_ROS EFFECT_BB
         MMP9 MMP2 ELAST COLLAG VSMC ILT DIAM
         TNF ROSO MAC SBP_mmHg Wall_Stress_idx Rupture_P MMP9_reduction
'

## ============================================================================
## COMPILE MODEL
## ============================================================================
mod <- mcode("AAA_QSP", code)

## ============================================================================
## DOSING EVENTS
## ============================================================================
# Standard dosing regimens (convert to mg administered at gut compartment)
# Doxycycline 100 mg BID → as 100mg doses every 12h
# Statin (Simvastatin 40mg) → once daily (every 24h)
# Propranolol 40mg → three times daily (every 8h)

make_doxycycline_dose <- function(dose_mg = 100, interval_h = 12, duration_days = 365) {
  ev(cmt = "GUT_DOXY", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

make_statin_dose <- function(dose_mg = 40, interval_h = 24, duration_days = 365) {
  ev(cmt = "GUT_STAT", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

make_propranolol_dose <- function(dose_mg = 40, interval_h = 8, duration_days = 365) {
  ev(cmt = "GUT_BB", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

## ============================================================================
## SIMULATION SETTINGS
## ============================================================================
sim_time <- seq(0, 8760, by = 24)   # 1 year (365 days), 24h steps

## ============================================================================
## TREATMENT SCENARIOS
## ============================================================================
# Scenario 1: No treatment (natural history)
scenario1 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 0, PROP_ON = 0) %>%
  mrgsim(end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "1: No Treatment")

# Scenario 2: Doxycycline monotherapy (100mg BID)
ev_doxy <- make_doxycycline_dose(100, 12, 365)
scenario2 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 0, PROP_ON = 0) %>%
  mrgsim(events = ev_doxy, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "2: Doxycycline 100mg BID")

# Scenario 3: Statin monotherapy (Simvastatin 40mg QD)
ev_stat <- make_statin_dose(40, 24, 365)
scenario3 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 1, PROP_ON = 0) %>%
  mrgsim(events = ev_stat, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "3: Simvastatin 40mg QD")

# Scenario 4: Propranolol monotherapy (40mg TID)
ev_prop <- make_propranolol_dose(40, 8, 365)
scenario4 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 0, PROP_ON = 1) %>%
  mrgsim(events = ev_prop, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "4: Propranolol 40mg TID")

# Scenario 5: Dual therapy — Doxycycline + Statin
ev_doxy_stat <- ev_doxy + ev_stat
scenario5 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 1, PROP_ON = 0) %>%
  mrgsim(events = ev_doxy_stat, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "5: Doxycycline + Statin")

# Scenario 6: Triple therapy — Doxycycline + Statin + Propranolol
ev_triple <- ev_doxy + ev_stat + ev_prop
scenario6 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 1, PROP_ON = 1) %>%
  mrgsim(events = ev_triple, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "6: Triple Therapy")

## ============================================================================
## COMBINE RESULTS
## ============================================================================
all_scenarios <- bind_rows(
  scenario1, scenario2, scenario3, scenario4, scenario5, scenario6
) %>%
  mutate(Day = time / 24)

## ============================================================================
## PLOTTING
## ============================================================================
cols_scen <- c(
  "1: No Treatment"          = "#E41A1C",
  "2: Doxycycline 100mg BID" = "#377EB8",
  "3: Simvastatin 40mg QD"   = "#4DAF4A",
  "4: Propranolol 40mg TID"  = "#FF7F00",
  "5: Doxycycline + Statin"  = "#984EA3",
  "6: Triple Therapy"        = "#A65628"
)

p_diam <- ggplot(all_scenarios, aes(Day, DIAM, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 10, y = 56.5, label = "Surgical threshold (55 mm)",
           size = 3, hjust = 0, color = "red") +
  scale_color_manual(values = cols_scen) +
  labs(title = "A. Aortic Diameter Over Time",
       x = "Day", y = "Max Aortic Diameter (mm)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

p_mmp9 <- ggplot(all_scenarios, aes(Day, MMP9, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "B. Active MMP-9 Plasma Level",
       x = "Day", y = "Active MMP-9 (ng/mL)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_elast <- ggplot(all_scenarios, aes(Day, ELAST, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "C. Aortic Elastin Content",
       x = "Day", y = "Elastin (% Baseline)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_rupture <- ggplot(all_scenarios, aes(Day, Rupture_P * 100, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "D. Estimated Rupture Risk",
       x = "Day", y = "Rupture Probability (%)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_combined <- (p_diam / (p_mmp9 | p_elast | p_rupture)) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title    = "AAA QSP Model — Six Treatment Scenarios (1-Year Simulation)",
    subtitle = paste("Initial diameter:", mod@param@data$DIAM0, "mm |",
                     "Baseline SBP:", mod@param@data$SBP0, "mmHg"),
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(p_combined)

## ============================================================================
## SUMMARY TABLE AT 1 YEAR (Day 365)
## ============================================================================
summary_1yr <- all_scenarios %>%
  filter(abs(Day - 365) < 0.5) %>%
  group_by(Scenario) %>%
  summarise(
    Diameter_mm      = round(mean(DIAM), 2),
    MMP9_ngmL        = round(mean(MMP9), 3),
    Elastin_pct      = round(mean(ELAST), 1),
    VSMC_pct         = round(mean(VSMC), 1),
    ILT_mL           = round(mean(ILT), 2),
    Rupture_Risk_pct = round(mean(Rupture_P) * 100, 2),
    SBP_mmHg         = round(mean(SBP_mmHg), 1),
    .groups = "drop"
  )

cat("\n=== 1-Year Simulation Summary ===\n")
print(summary_1yr, n = 10, width = 120)

## ============================================================================
## PK PROFILES (Single-dose illustration, first 48h)
## ============================================================================
ev_pk_test <- make_doxycycline_dose(100, 12, 3)
pk_test <- mod %>%
  mrgsim(events = ev_pk_test, end = 48, delta = 0.5) %>%
  as_tibble()

p_pk <- ggplot(pk_test, aes(time)) +
  geom_line(aes(y = Cp_doxy_plasma, color = "Plasma (Cp)"), linewidth = 1) +
  geom_line(aes(y = C_DOXY, color = "Aortic Tissue (Ct)"), linewidth = 1) +
  scale_color_manual(values = c("Plasma (Cp)" = "#2C7BB6", "Aortic Tissue (Ct)" = "#D7191C")) +
  labs(title = "Doxycycline PK — 100mg BID (First 48 Hours)",
       x = "Time (h)", y = "Concentration (mg/L)", color = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_pk)

## ============================================================================
## SENSITIVITY ANALYSIS: Initial Diameter Effect
## ============================================================================
diam_scenarios <- lapply(c(30, 35, 40, 45, 50), function(d) {
  mod %>%
    param(DIAM0 = d) %>%
    mrgsim(events = ev_triple, end = 8760, delta = 24) %>%
    as_tibble() %>%
    mutate(Initial_DIAM = paste0(d, "mm"), Day = time / 24)
})

diam_df <- bind_rows(diam_scenarios)

p_sens <- ggplot(diam_df, aes(Day, DIAM, color = Initial_DIAM)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Sensitivity Analysis: Initial Diameter Effect on Progression (Triple Therapy)",
       x = "Day", y = "Aortic Diameter (mm)", color = "Initial\nDiameter") +
  theme_bw(base_size = 12)

print(p_sens)

cat("\nAAA QSP Model simulation complete.\n")
cat("Key findings:\n")
cat("  - Natural history: ~", round(mod@param@data$k_diam_grow * 365, 1), "mm/year expansion\n")
cat("  - Doxycycline targets MMP-9 (Emax=", mod@param@data$EMAX_DOXY_MMP9,
    ", EC50=", mod@param@data$EC50_DOXY_MMP9, "mg/L AORTIC TISSUE, not plasma)\n")
cat("  - Triple therapy (doxy+statin+BB) provides maximal diameter growth reduction\n")
