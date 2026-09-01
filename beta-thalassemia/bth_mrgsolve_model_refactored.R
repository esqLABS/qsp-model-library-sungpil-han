## ============================================================
## Beta-Thalassemia QSP Model — mrgsolve Implementation
## ============================================================
## Disease: Beta-Thalassemia (Transfusion-Dependent [TDT] &
##          Non-Transfusion-Dependent [NTDT])
##
## Model scope:
##   - Erythropoiesis cascade (HSC → mature RBC): 9 compartments
##   - EPO feedback loop
##   - Iron metabolism & hepcidin-ERFE axis: 5 compartments
##   - PK: Luspatercept (2-cpt SC), Deferasirox (1-cpt PO),
##          Hydroxyurea (1-cpt PO)
##   - Transfusion dosing events
##   - 22 ODEs total (>15 requirement satisfied)
##
## Clinical calibration references:
##   - BELIEVE trial (luspatercept TDT, N Engl J Med 2020)
##   - BEYOND trial (luspatercept NTDT, N Engl J Med 2022)
##   - ESCALATOR trial (deferasirox, Blood 2008)
##   - Musallam et al. Lancet 2012 (NTDT natural history)
##   - Cazzola et al. Blood 2016 (ERFE/hepcidin biology)
##
## Author: QSP Library (CCR Auto-generated, 2026-06-25)
## ============================================================
##
## PK/PD REFACTOR (see FORK_WORKFLOW_GUIDE.md Part 2 and
## bth_refactor_notes.md for full detail):
##   All three compounds -- Luspatercept (L -> stem LUSP, spelled out
##   because the census's bare "L" abbreviation is not distinctive enough
##   to be a usable stem), Deferasirox (DFX), Hydroxyurea (HU) -- refactored
##   per their rows in driver-patches/data/compound_perturbation_census.md
##   ("Normalize duplicate concentration sites, then redirect").
##   Archetype 3 (depot+central+peripheral, linear) for luspatercept;
##   archetype 3's no-peripheral variant (depot+central, linear) for
##   deferasirox and hydroxyurea. Naming: GUT_<STEM>/CENT_<STEM>/
##   PERI_<STEM>/KA_<STEM>/F_<STEM>/CL_<STEM>/V1_<STEM>/V2_<STEM>/Q_<STEM>/
##   C_<STEM>/EC50_<STEM>/EMAX_<STEM>/GAMMA_<STEM>/EFFECT_<STEM>.
##
##   THE DUPLICATE-CONCENTRATION-SITE PROBLEM: each compound's plasma
##   concentration was computed twice -- once cached in $ODE (`C_LUSPAT`,
##   `C_DFX`, `C_HU`) for use by the disease-effect terms, and once more,
##   completely independently, in $TABLE (`C_Luspa`, `C_Deferasirox`,
##   `C_HU_calc`) purely for reporting. Luspatercept's disease effect was
##   ALSO independently re-derived a second time in $TABLE (`IE_frac_eff`),
##   separately from the $ODE-computed effect that the erythropoiesis ODEs
##   actually use.
##
##   FIX: `C_LUSP`/`C_DFX`/`C_HU` and `EFFECT_LUSP`/`EFFECT_DFX`/`EFFECT_HU`
##   are now declared as $GLOBAL preprocessor macros (`#define`), not
##   $ODE-local doubles -- so the identical definition is textually
##   substituted at every use site in BOTH $ODE and $TABLE. This is true
##   single-source-of-truth normalization (one definition, many use sites)
##   rather than two independently-typed re-derivations that merely happen
##   to agree today. It also sidesteps the dose-instant stale-value
##   reporting artifact documented in FORK_WORKFLOW_GUIDE.md ("The
##   dose-instant reporting artifact") the same way
##   clostridioides-difficile-infection's refactor found: a macro has no
##   per-block-evaluation cache to go stale, because it is not a variable.
##   $TABLE's old `C_Luspa`/`C_Deferasirox`/`C_HU_calc` local doubles are
##   gone (superseded by directly `$CAPTURE`-ing the `C_<STEM>` macros);
##   `IE_frac_eff` is kept (it is a genuinely distinct disease-level
##   reporting quantity, not just a duplicate concentration) but its RHS
##   now reads the canonical `EFFECT_LUSP` macro instead of re-deriving the
##   ratio inline -- see bth_refactor_notes.md for the important caveat that
##   `IE_frac_eff` was ALREADY, in the original, a luspatercept-only partial
##   recompute (it never included hydroxyurea's compounding reduction that
##   `ie_eff2` in $ODE applies) -- that pre-existing behavior is preserved
##   exactly, not "fixed".
##
##   Luspatercept's and hydroxyurea's effect terms were already bare
##   Emax/EC50 ratios (`Emax_L*C/(EC50_L+C)`, `Emax_HU*C/(EC50_HU+C)`) --
##   promoting them to the guide's explicit Hill form
##   (EMAX*C^GAMMA/(EC50^GAMMA+C^GAMMA), GAMMA=1) is a RENAME, not a refit
##   (GAMMA_LUSP=1/GAMMA_HU=1 are new named params; pow(x,1)=x so no
##   numeric change). Deferasirox's chelation effect
##   (`k_DFX_Fe*C_DFX`, renamed `EFFECT_DFX = KFE_DFX*C_DFX`) is
##   deliberately KEPT LINEAR/unbounded, not forced into the saturating
##   Hill template -- the original's own mechanism has no saturation term,
##   so imposing one would silently cap a mechanism the original leaves
##   uncapped (same reasoning as alkaptonuria's nitisinone term, see
##   alkaptonuria/aku_refactor_notes.md).
##
##   Verified against the original (rerun via the qspserver mrgsolve_api
##   service, http://localhost:8007) on the original file's own dosing
##   scenarios (2/3/4/5, plus the orphaned `hu_dose` event object -- see
##   bth_refactor_notes.md for why hydroxyurea needed its own check): every
##   $CAPTURE'd/compartment output matched the original EXACTLY (max abs
##   diff 0.0) -- see bth_refactor_notes.md for full detail.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL CODE
## ============================================================

code <- '
$PLUGIN nm-vars
// [PK/PD refactor -- build-compatibility fix, not a PK/PD change] the
// original also declares `autodiff`, which this mrgsolve build does not
// provide ("plugin autodiff could not be found"); dropped here per the
// syntax-only-fix policy in FORK_WORKFLOW_GUIDE.md. Logged as
// translations/UPSTREAM_ISSUES.md #131 (also covers two more pre-existing,
// unrelated build defects fixed below: $CAPTURE duplicating $CMT
// compartment names, and pmax() used inside $ODE).

$PARAM
// -------------------------------------------------------
// Erythropoiesis — baseline (non-thalassemia healthy adult)
// -------------------------------------------------------
k_BFU_prod  = 1.5      // BFU-E production from HSC (cells/uL/day)
k_BFU_diff  = 0.28     // BFU-E → CFU-E differentiation rate (/day)
k_CFU_diff  = 0.42     // CFU-E → ProEB differentiation (/day)
k_pro_diff  = 0.50     // ProEB → BasoEB (/day)
k_baso_diff = 0.45     // BasoEB → PolyEB (/day)
k_poly_diff = 0.40     // PolyEB → OrthoEB (/day)
k_ortho_diff= 0.55     // OrthoEB → Reticulocyte (/day)
k_retic_mat = 2.0      // Reticulocyte → Mature RBC maturation (/day)
k_rbc_elim  = 0.0083   // Mature RBC elimination: t½ ~ 120d (/day)
                        // [BTH: effectively shorter due to hemolysis]

// -------------------------------------------------------
// Ineffective erythropoiesis (IE) — BTH penalty
// ie_frac: fraction of erythroblasts that undergo apoptosis
// (0 = normal; 0.7 = severe BTH ≈ 70% IE)
// -------------------------------------------------------
ie_frac     = 0.70     // IE fraction (β0/β0 severe: 0.70–0.90)

// -------------------------------------------------------
// Hemoglobin
// -------------------------------------------------------
HGB_per_RBC = 28e-6    // g Hb per RBC (28 pg = 28e-6 µg; convert later)
BLOOD_VOL   = 4.5      // blood volume (L) for Hb calc
RBC_ref     = 4.5e6    // normal RBC count (cells/µL) → Hb ~14 g/dL

// -------------------------------------------------------
// EPO feedback
// -------------------------------------------------------
EPO_base    = 15.0     // basal EPO (IU/L) — healthy
EPO_max     = 8000.0   // max EPO (IU/L) — severe BTH
HGB_target  = 14.0     // set-point Hb (g/dL) for EPO feedback
kEPO_elim   = 0.35     // EPO elimination (/day, t½ ~2h → ~8h effective)
EC50_EPO    = 8.0      // EC50 of EPO on BFU-E proliferation (IU/L)
Emax_EPO    = 3.5      // Emax EPO effect on BFU-E (fold-increase)

// -------------------------------------------------------
// ERFE — Erythroferrone
// -------------------------------------------------------
k_ERFE_prod = 0.15     // ERFE production proportional to erythroblast pool
k_ERFE_elim = 0.35     // ERFE elimination (/day, t½ ~2d)

// -------------------------------------------------------
// Hepcidin
// -------------------------------------------------------
HEPCIDIN_base   = 25.0 // basal hepcidin (nM) healthy
k_HEPC_prod     = 0.08 // hepcidin synthesis (/day)
k_HEPC_elim     = 0.35 // hepcidin elimination (/day)
IC50_ERFE_HEPC  = 2.0  // ERFE IC50 for hepcidin suppression (relative units)
Emax_ERFE_HEPC  = 0.90 // max suppression of hepcidin by ERFE (90%)

// -------------------------------------------------------
// Iron metabolism
// -------------------------------------------------------
k_Fe_absorb  = 1.2    // dietary Fe absorption coefficient (µmol/L/day)
                       // modulated by hepcidin/ferroportin
k_Fe_RBC     = 0.25   // Fe incorporation into RBCs (/day)
k_Fe_liver   = 0.04   // plasma Fe → liver storage (/day)
k_Fe_release = 0.01   // liver Fe → plasma (basal mobilization, /day)
Fe_plasma_0  = 18.0   // normal plasma Fe (µmol/L, Tf-bound)
FERRITIN_0   = 100.0  // normal serum ferritin (µg/L)
k_FERR_form  = 0.005  // liver Fe → ferritin conversion
k_FERR_elim  = 0.04   // ferritin elimination (/day)
LIC_0        = 1.5    // normal LIC (mg Fe/g dw)
LIC_max      = 40.0   // max LIC capacity

// Cardiac iron
k_FeCARD_in  = 0.002  // NTBI → cardiac iron uptake (fraction/day)
k_FeCARD_out = 0.008  // cardiac Fe clearance (/day)

// -------------------------------------------------------
// Luspatercept PK (2-cpt SC)
// Ref: Platzbecker et al. NEJM 2017 (MDS), Viprakasit NEJM 2020
// F=60%, t½_terminal=~11 days, CL=0.35 L/day, Vc=~8L, Vp=~20L
// Q=0.9 L/day
// [PK/PD refactor] renamed to the GUT_/CENT_/PERI_/KA_/F_/CL_/V1_/V2_/Q_
// <STEM> convention: Ka_L->KA_LUSP, F_L->F_LUSP, CL_L->CL_LUSP,
// Vc_L->V1_LUSP, Vp_L->V2_LUSP, Q_L->Q_LUSP. No value changed.
// -------------------------------------------------------
KA_LUSP = 0.27        // SC absorption (/day)
F_LUSP  = 0.60        // bioavailability
CL_LUSP = 0.35        // clearance (L/day / 70 kg normalised)
V1_LUSP = 8.0         // central volume (L)
V2_LUSP = 20.0        // peripheral volume (L)
Q_LUSP  = 0.90        // inter-compartmental clearance (L/day)
WT      = 70.0        // body weight (kg) for dose calculation

// Luspatercept PD
// [PK/PD refactor] EC50_L->EC50_LUSP, Emax_L->EMAX_LUSP (rename, same
// values). GAMMA_LUSP=1.0 is a NEW named param -- the original had no
// explicit Hill exponent (bare Emax*C/(EC50+C) ratio); pow(x,1)=x so this
// is not a refit.
EC50_LUSP   = 0.5        // EC50 for late erythropoiesis improvement (µg/mL)
EMAX_LUSP   = 0.65       // max reduction in IE fraction by luspatercept
GAMMA_LUSP  = 1.0        // Hill coefficient (new; original had none explicit)
// Effect: IE_effective = ie_frac * (1 - EFFECT_LUSP), see $GLOBAL

// -------------------------------------------------------
// Deferasirox PK (1-cpt oral)
// F=70%, Tmax~3.5h, t½~8-16h, CL~0.20 L/h/kg → 14 L/day/70kg
// [PK/PD refactor] Ka_DFX->KA_DFX, V_DFX->V1_DFX (rename, same values).
// k_DFX_Fe->KFE_DFX (rename; kept LINEAR, not forced into an Emax/EC50
// Hill form -- see $GLOBAL comment and bth_refactor_notes.md).
// -------------------------------------------------------
KA_DFX   = 1.50       // absorption (/day → converts from /h=0.18)
F_DFX    = 0.70
CL_DFX   = 14.0       // L/day
V1_DFX   = 100.0      // L (apparent Vd, protein-bound 99%)
KFE_DFX  = 0.12       // DFX → Fe chelation efficiency on LIC (/day per µg/mL)

// -------------------------------------------------------
// Hydroxyurea PK (1-cpt oral)
// F=80%, t½~4h, Vd=0.6 L/kg → 42L, CL=7.3 L/h = 175 L/day
// [PK/PD refactor] Ka_HU->KA_HU, V_HU->V1_HU, Emax_HU->EMAX_HU (rename,
// same values). GAMMA_HU=1.0 is a NEW named param, same reasoning as
// GAMMA_LUSP above.
// -------------------------------------------------------
KA_HU    = 6.0         // absorption (/day, fast: /h = 0.25)
F_HU     = 0.80
CL_HU    = 175.0
V1_HU    = 42.0
EC50_HU  = 10.0        // HU concentration for HbF induction (µg/mL)
EMAX_HU  = 0.40        // max HbF increase (fraction of total Hb)
GAMMA_HU = 1.0         // Hill coefficient (new; original had none explicit)

// -------------------------------------------------------
// Transfusion parameters
// (handled via addl/rate dosing events in R code)
// -------------------------------------------------------
TX_FE_per_unit = 225.0 // Fe per unit pRBC (µmol → ~200 mg Fe)
                        // 1 unit pRBC = 200 mL, 1 mg Fe/mL → ~200 mg
                        // Convert: 200 mg / 56 g/mol = 3571 µmol → per L blood vol

$GLOBAL
// [PK/PD refactor] Canonical, single-source-of-truth concentration/effect
// definitions for the three PK compounds, per
// driver-patches/data/compound_perturbation_census.md’s
// "Normalize duplicate concentration sites, then redirect" classification
// for beta-thalassemia|DFX, beta-thalassemia|HU, beta-thalassemia|L.
// Declared as $GLOBAL preprocessor macros (not $ODE-local doubles) so the
// SAME definition textually substitutes at every use site in both $ODE
// and $TABLE -- true single-site normalization, not two independently
// hand-written re-derivations kept in sync by convention. This also
// sidesteps the dose-instant stale-value reporting artifact documented in
// FORK_WORKFLOW_GUIDE.md ("The dose-instant reporting artifact"), the same
// fix clostridioides-difficile-infection’s refactor found: a macro has no
// per-block-evaluation cache to go stale.
#define C_LUSP  (CENT_LUSP / V1_LUSP)   // luspatercept plasma conc (ug/mL); was LUSPAT_C1/Vc_L, independently re-derived as C_LUSPAT ($ODE) and C_Luspa ($TABLE)
#define C_DFX   (CENT_DFX  / V1_DFX)    // deferasirox plasma conc (ug/mL); was DFX_CENT/V_DFX, independently re-derived as C_DFX ($ODE) and C_Deferasirox ($TABLE)
#define C_HU    (CENT_HU   / V1_HU)     // hydroxyurea plasma conc (ug/mL); was HU_CENT/V_HU, independently re-derived as C_HU ($ODE) and C_HU_calc ($TABLE)

// EFFECT_LUSP: exact rename of the original bare ratio
// IE_effect_L = Emax_L*C_LUSPAT/(EC50_L+C_LUSPAT); GAMMA_LUSP=1 is new
// (pow(x,1)=x, no numeric change).
#define EFFECT_LUSP (EMAX_LUSP * pow(C_LUSP, GAMMA_LUSP) / (pow(EC50_LUSP, GAMMA_LUSP) + pow(C_LUSP, GAMMA_LUSP)))

// EFFECT_DFX: exact rename of the original k_DFX_Fe*C_DFX. Deliberately
// kept LINEAR/unbounded -- the original’s own chelation-removal mechanism
// has no saturation term, so forcing it into the guide’s saturating
// Emax/EC50 Hill template would silently cap a mechanism the original
// leaves uncapped (same reasoning as alkaptonuria/aku_refactor_notes.md’s
// nitisinone term).
#define EFFECT_DFX  (KFE_DFX * C_DFX)

// EFFECT_HU: exact rename of the original bare ratio
// HU_HbF = Emax_HU*C_HU/(EC50_HU+C_HU); GAMMA_HU=1 is new (pow(x,1)=x, no
// numeric change).
#define EFFECT_HU   (EMAX_HU * pow(C_HU, GAMMA_HU) / (pow(EC50_HU, GAMMA_HU) + pow(C_HU, GAMMA_HU)))

$INIT
// Erythropoiesis compartments (healthy steady-state approximations)
BFU_E    = 5.0        // BFU-E progenitors (cells/µL)
CFU_E    = 3.5
PRO_E    = 7.0
BASO_E   = 14.0
POLY_E   = 14.0
ORTHO_E  = 14.0
RETIC    = 50.0       // reticulocytes (cells/µL)
RBC_MAT  = 4500.0     // mature RBCs (cells/µL, × 10^3 scale)
EPO_CMT  = 15.0       // EPO (IU/L)

// ERFE / Hepcidin
ERFE_CMT    = 1.0     // relative units
HEPC_CMT    = 25.0    // hepcidin (nM)

// Iron compartments
FE_PL    = 18.0       // plasma iron (µmol/L)
FE_LIV   = 1.5        // liver iron content (mg Fe/g dw)
FERR_CMT = 100.0      // serum ferritin (µg/L)
FE_CARD  = 0.1        // cardiac iron (relative, T2* inverse proxy)

// Luspatercept PK -- [PK/PD refactor] LUSPAT_SC->GUT_LUSP, LUSPAT_C1->CENT_LUSP, LUSPAT_C2->PERI_LUSP
GUT_LUSP  = 0.0
CENT_LUSP = 0.0
PERI_LUSP = 0.0

// Deferasirox PK -- [PK/PD refactor] DFX_GUT->GUT_DFX, DFX_CENT->CENT_DFX
GUT_DFX  = 0.0
CENT_DFX = 0.0

// Hydroxyurea PK -- [PK/PD refactor] HU_GUT->GUT_HU, HU_CENT->CENT_HU
GUT_HU   = 0.0
CENT_HU  = 0.0

$ODE
// -------------------------------------------------------
// [PK/PD refactor] Concentration and disease-effect quantities for the
// three PK compounds now come from the $GLOBAL macros above (C_LUSP,
// C_DFX, C_HU, EFFECT_LUSP, EFFECT_DFX, EFFECT_HU) -- this is the single
// site those quantities are defined; every read below (and every read in
// $TABLE) uses the identical macro-expanded expression. No local doubles
// re-derive them here any more.
// -------------------------------------------------------

// -------------------------------------------------------
// Luspatercept effect on IE fraction
// (reduces late-stage apoptosis, improves ortho-EB output)
// -------------------------------------------------------
double ie_eff       = ie_frac * (1.0 - EFFECT_LUSP);  // effective IE fraction

// -------------------------------------------------------
// EPO effect on BFU-E proliferation
// -------------------------------------------------------
double EPO_stim = 1.0 + Emax_EPO * EPO_CMT / (EC50_EPO + EPO_CMT);

// -------------------------------------------------------
// ERFE → Hepcidin suppression
// -------------------------------------------------------
double ERFE_sup_HEPC = 1.0 - Emax_ERFE_HEPC * ERFE_CMT / (IC50_ERFE_HEPC + ERFE_CMT);

// -------------------------------------------------------
// Ferroportin activity (hepcidin-dependent)
// norm: HEPC_CMT = 25 nM → FPN activity = 1.0
// -------------------------------------------------------
double FPN_act = HEPCIDIN_base / (HEPCIDIN_base + HEPC_CMT);  // 0→1 (high hepc = low FPN)
double Fe_abs  = k_Fe_absorb * (1.0 - FPN_act);              // reduced absorption when hepc high

// -------------------------------------------------------
// Hydroxyurea effect: HbF induction (reduces effective ie_frac)
// (HU partially compensates α/β imbalance via γ-globin)
// -------------------------------------------------------
double ie_eff2 = ie_eff * (1.0 - 0.3 * EFFECT_HU);  // partial IE reduction by HbF

// -------------------------------------------------------
// Hemoglobin (g/dL) — derived
// Hb = RBC_MAT (cells/uL) * Hb_per_RBC / 100
// Using simplified: Hb = RBC_MAT / RBC_ref * 14.0
// -------------------------------------------------------
double HGB = (RBC_MAT / RBC_ref) * 14.0;

// -------------------------------------------------------
// Total erythroblast pool for ERFE production
// -------------------------------------------------------
double ERYTHRO_POOL = PRO_E + BASO_E + POLY_E + ORTHO_E;

// -------------------------------------------------------
// ERFE derived from erythroblast pool
// -------------------------------------------------------
double ERFE_prod = k_ERFE_prod * ERYTHRO_POOL;

// -------------------------------------------------------
// EPO secretion — feedback based on Hb
// Hill-type: EPO rises as Hb falls
// -------------------------------------------------------
double EPO_prod = EPO_base + (EPO_max - EPO_base) * pow(HGB_target, 4) /
                              (pow(HGB_target, 4) + pow(HGB + 0.01, 4));

// -------------------------------------------------------
// ERYTHROBLAST ODEs (linear cascade with IE losses)
// -------------------------------------------------------
dxdt_BFU_E   = k_BFU_prod * EPO_stim                    // production (EPO-driven)
              - k_BFU_diff * BFU_E;                       // differentiation

dxdt_CFU_E   = k_BFU_diff * BFU_E
              - k_CFU_diff * CFU_E;

dxdt_PRO_E   = k_CFU_diff * CFU_E
              - k_pro_diff * PRO_E
              - ie_eff2 * k_pro_diff * PRO_E;             // IE apoptosis at early stage

dxdt_BASO_E  = k_pro_diff  * PRO_E * (1.0 - ie_eff2 * 0.3)  // some survive
              - k_baso_diff * BASO_E
              - ie_eff2 * k_baso_diff * BASO_E;

dxdt_POLY_E  = k_baso_diff * BASO_E * (1.0 - ie_eff2 * 0.3)
              - k_poly_diff * POLY_E
              - ie_eff2 * 0.5 * k_poly_diff * POLY_E;    // IE less at late stage

dxdt_ORTHO_E = k_poly_diff * POLY_E  * (1.0 - ie_eff2 * 0.15)
              - k_ortho_diff * ORTHO_E;

dxdt_RETIC   = k_ortho_diff * ORTHO_E
              - k_retic_mat * RETIC;

dxdt_RBC_MAT = k_retic_mat * RETIC
              - k_rbc_elim  * RBC_MAT;

// -------------------------------------------------------
// EPO ODE
// -------------------------------------------------------
dxdt_EPO_CMT = EPO_prod - kEPO_elim * EPO_CMT;

// -------------------------------------------------------
// ERFE ODE
// -------------------------------------------------------
dxdt_ERFE_CMT = ERFE_prod - k_ERFE_elim * ERFE_CMT;

// -------------------------------------------------------
// HEPCIDIN ODE
// -------------------------------------------------------
dxdt_HEPC_CMT = k_HEPC_prod * HEPCIDIN_base * ERFE_sup_HEPC
              - k_HEPC_elim * HEPC_CMT;

// -------------------------------------------------------
// IRON compartments
// -------------------------------------------------------
// Plasma iron — absorption, RBC use, liver storage
dxdt_FE_PL   = Fe_abs                                    // dietary absorption
              + k_Fe_release * FE_LIV                     // release from liver
              - k_Fe_RBC    * FE_PL                       // RBC synthesis uptake
              - k_Fe_liver  * FE_PL;                      // liver storage

// Liver iron content — storage + transfusion iron input (via dosing event FE_LIV)
dxdt_FE_LIV  = k_Fe_liver * FE_PL                        // from plasma
              - k_Fe_release * FE_LIV                     // remobilization
              - EFFECT_DFX  * FE_LIV;                     // chelation removal [was DFX_Fe_remov * FE_LIV]

// Serum ferritin — derived from liver iron
dxdt_FERR_CMT= k_FERR_form * FE_LIV                      // proportional to LIC
              - k_FERR_elim * FERR_CMT;

// Cardiac iron — driven by NTBI when LIC > threshold
double NTBI_gen = fmax(0.0, (FE_LIV - 7.0) * 0.05);     // NTBI only when LIC > 7 -- [PK/PD refactor -- build-compat fix] pmax (R syntax) -> fmax (C++); see UPSTREAM_ISSUES.md #131
dxdt_FE_CARD  = k_FeCARD_in  * NTBI_gen
               - k_FeCARD_out * FE_CARD;

// -------------------------------------------------------
// LUSPATERCEPT PK (2-cpt SC) -- [PK/PD refactor] renamed compartments/params only
// -------------------------------------------------------
dxdt_GUT_LUSP  = -KA_LUSP * GUT_LUSP;

dxdt_CENT_LUSP = KA_LUSP * F_LUSP * GUT_LUSP
                - (CL_LUSP + Q_LUSP) / V1_LUSP * CENT_LUSP
                + Q_LUSP / V2_LUSP * PERI_LUSP;

dxdt_PERI_LUSP = Q_LUSP / V1_LUSP * CENT_LUSP
                - Q_LUSP / V2_LUSP * PERI_LUSP;

// -------------------------------------------------------
// DEFERASIROX PK (1-cpt oral) -- [PK/PD refactor] renamed compartments/params only
// -------------------------------------------------------
dxdt_GUT_DFX  = -KA_DFX * GUT_DFX;
dxdt_CENT_DFX =  KA_DFX * F_DFX * GUT_DFX
               - CL_DFX / V1_DFX * CENT_DFX;

// -------------------------------------------------------
// HYDROXYUREA PK (1-cpt oral) -- [PK/PD refactor] renamed compartments/params only
// -------------------------------------------------------
dxdt_GUT_HU   = -KA_HU * GUT_HU;
dxdt_CENT_HU  =  KA_HU * F_HU * GUT_HU
               - CL_HU / V1_HU * CENT_HU;

$TABLE
// Report calculated quantities
double Hb_gdL     = (RBC_MAT / RBC_ref) * 14.0;
double EPO_IUL    = EPO_CMT;
double ERFE_rel   = ERFE_CMT;
double HEPC_nM    = HEPC_CMT;
double LIC_mgFe   = FE_LIV;
double FERR_ugL   = FERR_CMT;
double CARD_T2star = 50.0 / (FE_CARD + 0.01);  // proxy: T2* inversely proportional
double Retic_pct  = 100.0 * RETIC / (RETIC + RBC_MAT / 100.0);
// [PK/PD refactor] C_Luspa/C_Deferasirox/C_HU_calc removed -- superseded by
// directly $CAPTURE-ing the canonical C_LUSP/C_DFX/C_HU macros (see
// $GLOBAL and $CAPTURE below), which are now the single site those
// concentrations are defined, reused verbatim here.
// IE_frac_eff kept as its own distinct reporting quantity (NOT the same
// thing as EFFECT_LUSP: it applies EFFECT_LUSP to ie_frac). Its RHS now
// reads the canonical EFFECT_LUSP macro instead of re-deriving the ratio
// inline; this is exactly the original quantity, unchanged numerically.
// NOTE (preserved, not "fixed"): this was already, in the original, a
// luspatercept-only partial recompute -- it never included hydroxyurea’s
// additional reduction that $ODE-side ie_eff2 applies. See
// bth_refactor_notes.md.
double IE_frac_eff = ie_frac * (1.0 - EFFECT_LUSP);

$CAPTURE Hb_gdL EPO_IUL ERFE_rel HEPC_nM LIC_mgFe FERR_ugL CARD_T2star
         C_LUSP C_DFX C_HU EFFECT_LUSP EFFECT_DFX EFFECT_HU
         Retic_pct IE_frac_eff
// [PK/PD refactor -- build-compat fix] the original also listed
// BFU_E CFU_E PRO_E BASO_E POLY_E ORTHO_E RETIC RBC_MAT FE_PL FE_CARD here;
// those are $CMT-implicit compartments (declared via $INIT), and mrgsolve
// 2.0.1 rejects a compartment name inside $CAPTURE ("compartment should
// not be in $CAPTURE"). Removed -- all ten are still reported, unchanged,
// as compartments (mrgsolve concatenates compartments + captures into
// outputPaths, so nothing becomes unobservable). See UPSTREAM_ISSUES.md #131.
'

## ============================================================
## COMPILE MODEL
## ============================================================
mod_bth <- mcode("beta_thalassemia_qsp_refactored", code)

## ============================================================
## DOSING EVENTS
## ============================================================

## --- Luspatercept: 1.0 mg/kg SC q21d (BELIEVE trial starting dose)
## [PK/PD refactor] cmt "LUSPAT_SC" -> "GUT_LUSP"
luspat_dose_mgkg <- 1.0
WT_kg <- 70
luspat_dose_total <- luspat_dose_mgkg * WT_kg  # mg → µg*1000 in model units
luspatercept_ev <- ev(
  ID    = 1,
  amt   = luspat_dose_total * 1000,  # µg
  cmt   = "GUT_LUSP",
  addl  = 11,     # 12 doses = ~36 weeks
  ii    = 21      # every 21 days
)

## --- Deferasirox: 30 mg/kg/day PO daily
## [PK/PD refactor] cmt "DFX_GUT" -> "GUT_DFX"
dfx_dose <- ev(
  ID    = 1,
  amt   = 30 * WT_kg,   # mg
  cmt   = "GUT_DFX",
  addl  = 364,
  ii    = 1
)

## --- Hydroxyurea: 20 mg/kg/day PO daily
## [PK/PD refactor] cmt "HU_GUT" -> "GUT_HU"
hu_dose <- ev(
  ID    = 1,
  amt   = 20 * WT_kg,
  cmt   = "GUT_HU",
  addl  = 364,
  ii    = 1
)

## --- Transfusions: every 21 days, 2 units pRBC
## Transfusion modelled as a bolus to FE_LIV (iron load)
## and an immediate Hb boost (modelled via RBC_MAT reset in init or as bolus)
## Iron per unit: ~200 mg Fe = 3571 µmol; normalised to LIC units
## Approx: 2 units → LIC increase of ~0.3 mg/g dw / event
tx_iron_ev <- ev(
  ID    = 1,
  amt   = 0.30,    # LIC increase per transfusion event (mg Fe/g dw)
  cmt   = "FE_LIV",
  addl  = 17,      # 18 transfusion events over year
  ii    = 21
)

## ============================================================
## SCENARIOS
## ============================================================

## Helper: run simulation and tag with scenario name
run_scenario <- function(params = list(), events = NULL, days = 365, label = "Scenario") {
  m <- mod_bth
  if (length(params) > 0) m <- param(m, params)
  if (is.null(events)) {
    out <- mrgsim(m, end = days, delta = 1)
  } else {
    out <- mrgsim(m, events = events, end = days, delta = 1)
  }
  as_tibble(out) %>% mutate(Scenario = label)
}

## Scenario 1: Natural history — severe TDT (no treatment)
sc1_params <- list(ie_frac = 0.80, EPO_base = 50)
sc1 <- run_scenario(
  params = sc1_params,
  days   = 365,
  label  = "1. Natural History (TDT, no Tx)"
)

## Scenario 2: Regular transfusions only (q21d, no chelation)
sc2_params <- list(ie_frac = 0.80, EPO_base = 50)
sc2 <- run_scenario(
  params = sc2_params,
  events = tx_iron_ev,
  days   = 365,
  label  = "2. Transfusions Only (no chelation)"
)

## Scenario 3: Transfusions + Deferasirox 30 mg/kg/day
sc3_events <- ev_seq(tx_iron_ev, dfx_dose)
sc3_params  <- list(ie_frac = 0.80, EPO_base = 50)
sc3 <- run_scenario(
  params = sc3_params,
  events = sc3_events,
  days   = 365,
  label  = "3. Transfusions + Deferasirox"
)

## Scenario 4: Luspatercept monotherapy (NTDT, milder)
sc4_params <- list(ie_frac = 0.55, EPO_base = 30)
sc4 <- run_scenario(
  params = sc4_params,
  events = luspatercept_ev,
  days   = 252,   # 12 cycles × 21 days
  label  = "4. Luspatercept (NTDT, 1.0 mg/kg q21d)"
)

## Scenario 5: Luspatercept + Transfusions + Deferasirox (TDT)
sc5_params  <- list(ie_frac = 0.80, EPO_base = 50)
sc5_events  <- ev_seq(luspatercept_ev, tx_iron_ev, dfx_dose)
sc5 <- run_scenario(
  params = sc5_params,
  events = sc5_events,
  days   = 252,
  label  = "5. Luspatercept + Tx + Deferasirox (TDT)"
)

## Scenario 6: Gene Therapy (simulated by restoring normal parameters)
## After engraftment (~3-6 months post-infusion), BTH phenotype largely corrected
sc6_params <- list(ie_frac = 0.05, EPO_base = 15)   # near-normal
sc6 <- run_scenario(
  params = sc6_params,
  days   = 365,
  label  = "6. Gene Therapy (beti-cel, engrafted)"
)

## Combine all scenarios
all_sc <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6)

## ============================================================
## PLOTS
## ============================================================

theme_set(theme_bw(base_size = 12))

## Hemoglobin over time
p_hb <- ggplot(all_sc, aes(x = time, y = Hb_gdL, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 9.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 200, y = 10.0, label = "Pre-Tx target ≥9.5 g/dL", size = 3) +
  labs(title = "Beta-Thalassemia QSP: Hemoglobin Response",
       x = "Time (days)", y = "Hemoglobin (g/dL)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Liver Iron Content
p_lic <- ggplot(all_sc, aes(x = time, y = LIC_mgFe, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 15, linetype = "dotted", color = "red") +
  annotate("text", x = 200, y = 8.0,  label = "LIC 7 mg/g (chelation threshold)", size = 3) +
  annotate("text", x = 200, y = 16.0, label = "LIC 15 mg/g (high risk)", size = 3, color = "red") +
  labs(title = "Liver Iron Content (LIC)", x = "Time (days)", y = "LIC (mg Fe/g dw)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## EPO and ERFE
p_epo_erfe <- all_sc %>%
  select(time, Scenario, EPO_IUL, ERFE_rel) %>%
  pivot_longer(c(EPO_IUL, ERFE_rel), names_to = "Biomarker", values_to = "Value") %>%
  ggplot(aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Biomarker, scales = "free_y",
             labeller = labeller(Biomarker = c(EPO_IUL = "EPO (IU/L)",
                                               ERFE_rel = "ERFE (rel. units)"))) +
  labs(title = "EPO & Erythroferrone Dynamics", x = "Time (days)", y = "") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Hepcidin over time
p_hepc <- ggplot(all_sc, aes(x = time, y = HEPC_nM, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "gray40") +
  labs(title = "Hepcidin Dynamics", x = "Time (days)", y = "Hepcidin (nM)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Luspatercept PK (Scenario 4 & 5)
## [PK/PD refactor] y = C_Luspa -> y = C_LUSP (canonical $CAPTURE name)
p_luspat_pk <- all_sc %>%
  filter(grepl("Luspatercept|Gene", Scenario)) %>%
  ggplot(aes(x = time, y = C_LUSP, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Luspatercept Plasma Concentration",
       x = "Time (days)", y = "Luspatercept (µg/mL)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

## Cardiac iron proxy (T2*)
p_t2star <- ggplot(all_sc, aes(x = time, y = CARD_T2star, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
  annotate("text", x = 300, y = 19, label = "T2* < 20 ms: cardiac iron overload risk", size = 3, color = "red") +
  labs(title = "Cardiac T2* (Iron Proxy)", x = "Time (days)", y = "Cardiac T2* (ms, higher=safer)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom")

## Clinical calibration summary
cat("\n==============================\n")
cat("CLINICAL CALIBRATION SUMMARY\n")
cat("==============================\n")
cat("\nScenario steady-state endpoints (day 250+):\n")
sum_tab <- all_sc %>%
  filter(time > 240) %>%
  group_by(Scenario) %>%
  summarise(
    Hb_mean     = round(mean(Hb_gdL,   na.rm=TRUE), 2),
    LIC_mean    = round(mean(LIC_mgFe, na.rm=TRUE), 2),
    EPO_mean    = round(mean(EPO_IUL,  na.rm=TRUE), 1),
    ERFE_mean   = round(mean(ERFE_rel, na.rm=TRUE), 2),
    HEPC_mean   = round(mean(HEPC_nM,  na.rm=TRUE), 2),
    T2star_mean = round(mean(CARD_T2star,na.rm=TRUE),1),
    .groups = "drop"
  )
print(sum_tab)

cat("\n--- BELIEVE Trial Calibration (Sc 5 vs reported) ---\n")
cat("Observed Hb increase (BELIEVE): +1.0 to +2.0 g/dL vs model\n")
sc5_dhb <- mean(sc5$Hb_gdL[sc5$time > 200]) - sc5$Hb_gdL[1]
cat(sprintf("  Model ΔHb (Sc 5 vs baseline): +%.2f g/dL\n", sc5_dhb))

cat("\n--- BEYOND Trial Calibration (Sc 4 NTDT) ---\n")
cat("Observed TI rate (BEYOND): 77.7% vs placebo 0%\n")
sc4_dhb <- mean(sc4$Hb_gdL[sc4$time > 200]) - sc4$Hb_gdL[1]
cat(sprintf("  Model ΔHb (Sc 4 NTDT luspatercept): +%.2f g/dL\n", sc4_dhb))

cat("\n--- Deferasirox LIC reduction (Sc 3 vs 2) ---\n")
sc3_lic <- mean(sc3$LIC_mgFe[sc3$time > 300])
sc2_lic <- mean(sc2$LIC_mgFe[sc2$time > 300])
cat(sprintf("  LIC with chelation: %.2f vs without: %.2f mg/g\n", sc3_lic, sc2_lic))

## Print plots
print(p_hb)
print(p_lic)
print(p_epo_erfe)
print(p_hepc)
print(p_luspat_pk)
print(p_t2star)

## ============================================================
## PARAMETER SENSITIVITY ANALYSIS
## ============================================================
cat("\n--- Sensitivity: IE fraction effect on Hb ---\n")
ie_vals <- seq(0.1, 0.90, by = 0.10)
sens_hb <- sapply(ie_vals, function(ie) {
  m2 <- param(mod_bth, ie_frac = ie, EPO_base = 15 + 50*ie)
  out2 <- mrgsim(m2, end = 365, delta = 1)
  mean(as_tibble(out2)$Hb_gdL[300:365])
})
sens_df <- data.frame(IE_frac = ie_vals, Hb_ss = sens_hb)
cat("IE Fraction vs Steady-State Hb:\n")
print(sens_df)

p_sens <- ggplot(sens_df, aes(x = IE_frac, y = Hb_ss)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(title = "Sensitivity: IE Fraction vs Steady-State Hemoglobin",
       x = "Ineffective Erythropoiesis Fraction",
       y = "Steady-State Hb (g/dL)") +
  theme_bw()
print(p_sens)
