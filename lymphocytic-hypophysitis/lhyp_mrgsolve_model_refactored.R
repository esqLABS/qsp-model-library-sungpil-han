## =============================================================================
## Lymphocytic Hypophysitis (LyH) — QSP mrgsolve Model
## File: lhyp_mrgsolve_model_refactored.R
##
## Disease: Autoimmune lymphocytic infiltration of the pituitary gland
## causing multi-axis hypopituitarism (ACTH, TSH, FSH/LH, GH, PRL, ADH)
##
## Model scope:
##   - Immune dynamics: naive T-cells, effector T (Th1/Th17), Tregs, B-cells,
##     plasma cells, anti-pituitary antibodies (APA)
##   - Pituitary: inflammatory mass volume, functional cell mass (0–1)
##   - HPA axis: ACTH → cortisol (negative feedback, diurnal rhythm)
##   - HPT axis: TSH → free T4 (negative feedback)
##   - HPG axis: FSH, LH, estradiol (simplified)
##   - GH/IGF-1 axis: GH → IGF-1 (negative feedback)
##   - Prolactin: elevated by stalk effect
##   - ADH: posterior pituitary involvement → diabetes insipidus
##   - Drug PK: prednisolone (2-compartment oral), azathioprine (1-compartment),
##     rituximab (1-compartment IV)
##
## Parameters calibrated to:
##   - Honegger J et al. (2015) Endocrine: LyH hormone deficiency rates
##   - Caputo M et al. (2019) Clin Exp Immunol: cytokine profiles
##   - Fleseriu M et al. (2021) Pituitary: treatment outcomes
##   - Faje A et al. (2014) J Clin Endocrinol Metab: checkpoint inhibitor LyH
##
## 22 ODE compartments, 5 treatment scenarios, R6+ compatible
## =============================================================================
##
## [refactor] Per FORK_WORKFLOW_GUIDE.md Part 2 (pluggable PK, named Hill
## interface). Per `driver-patches/data/compound_perturbation_census.md`,
## this file has exactly two rows, BOTH classified "Normalize duplicate
## concentration sites, then redirect":
##   - Azathioprine (stem AZA) -- target: purine synthesis
##   - Rituximab    (stem RTX) -- target: CD20 (checked for TMDD -- the
##     original models rituximab's PK as plain 1-compartment linear
##     clearance with an algebraic Hill-ratio effect term, NOT
##     target-mediated drug disposition; no receptor-binding ODE state
##     exists anywhere in the file for CD20)
## Prednisolone is a third compound in this file but has NO row of its own
## in the census and is completely out of scope: every `Pred_*`/`*_pred`/
## `E_pred`/`Cpred` name, value, and formula is left untouched.
##
## Both AZA and RTX had the identical defect: their plasma concentration
## was computed twice, under two different local names, in two different
## mrgsolve blocks that cannot see each other's locals ($MAIN and $TABLE):
##   AZA: `double Caza = AZA_plasma / Vc_aza;` ($MAIN) vs
##        `double Caza_obs = AZA_plasma / Vc_aza;` ($TABLE)
##   RTX: `double Crtx = RTX_plasma / Vc_rtx;` ($MAIN) vs
##        `double Crtx_obs = RTX_plasma / Vc_rtx;` ($TABLE)
## Normalized to one $GLOBAL-predeclared name each, C_AZA/C_RTX, assigned
## from the identical formula at both of the original's two call sites
## (matching the established pattern in hereditary-spherocytosis,
## cervical-cancer, breast-cancer, AMD, membranous-nephropathy, and
## stable-angina refactors -- $TABLE cannot read a $MAIN local or a value
## $MAIN wrote earlier in the same interval without risking a stale read,
## so the formula is kept at both of the original's own evaluation sites,
## just under one shared name instead of two independently-maintained
## copies). See lhyp_refactor_notes.md for full detail.
##
## Archetypes:
##   AZA: Archetype 3 minus peripheral (depot + central, linear) --
##        GUT_AZA -> CENT_AZA, no peripheral compartment in the original.
##   RTX: Archetype 1 (no depot, single compartment, linear elimination) --
##        CENT_RTX only.
## Both Hill effect terms were already the canonical Emax*C/(EC50+C) shape
## in the original (implicit gamma=1) -- renames, not refits.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- MODEL CODE (mrgsolve format) -------------------------------------------
code <- '
$PROB Lymphocytic Hypophysitis QSP Model
  Autoimmune pituitary disease with multi-axis endocrine disruption.
  22 ODE compartments, prednisolone/azathioprine/rituximab PK/PD.

$PARAM
  // ---- Prednisolone PK (out of scope -- no census row, untouched) ----
  ka     = 1.20,   // absorption rate constant (h^-1)
  Vc     = 45.0,   // central volume (L)
  Vp     = 120.0,  // peripheral volume (L)
  CL_P   = 8.50,   // total clearance (L/h)
  Q_P    = 25.0,   // inter-compartmental clearance (L/h)
  F_pred = 0.82,   // oral bioavailability

  // ---- Azathioprine PK (refactored: stem AZA) ----
  KA_AZA = 0.80,   // absorption rate constant (h^-1) [was ka_aza]
  V1_AZA = 30.0,   // central volume of distribution (L) [was Vc_aza]
  CL_AZA = 12.0,   // clearance (L/h) - includes TPMT metabolism [was CL_aza]
  F_AZA  = 0.50,   // oral bioavailability [was F_aza] -- NOTE: unused inside
                   // $ODE in the original (bioavailability was applied
                   // externally, by scaling the dosed amount in the R-side
                   // aza_events() helper, not inside the DSL). Preserved
                   // exactly as unused here too -- see refactor notes.

  // ---- Rituximab PK (refactored: stem RTX) ----
  V1_RTX = 3.1,    // central volume (L) - typical for IgG antibody [was Vc_rtx]
  CL_RTX = 0.015,  // clearance (L/h) [was CL_rtx]
  MW_RTX = 145000, // molecular weight (g/mol, not used directly) [was MW_rtx]

  // ---- Immune Cell Dynamics ----
  // Naive T cells
  kprolif_Tn  = 0.050, // naive T proliferation rate (day^-1)
  kdeath_Tn   = 0.020, // naive T death rate (day^-1)
  Tn_baseline = 1000,  // steady-state naive T (AU)

  // Effector T cells (Th1 + Th17 aggregate)
  kact_T      = 0.12,  // naive→effector activation (AU^-1 day^-1)
  kdeath_Te   = 0.08,  // effector T death rate (day^-1)
  kinact_T    = 0.05,  // effector→inactive T rate (day^-1)

  // Regulatory T cells
  kprolif_Tr  = 0.030, // Treg proliferation (day^-1)
  kdeath_Tr   = 0.025, // Treg death rate (day^-1)
  kinduct_Tr  = 0.015, // induction from naive (day^-1, TGF-β driven)

  // B cells and plasma cells
  kprolif_Bn  = 0.030, // naive B proliferation (day^-1)
  kdeath_Bn   = 0.015, // naive B death (day^-1)
  kact_B      = 0.04,  // B cell activation (day^-1)
  kdeath_Bp   = 0.020, // plasma cell death (day^-1)
  Bn_baseline = 500,   // steady-state naive B (AU)

  // Anti-pituitary antibodies (APA)
  kprod_APA   = 0.025, // APA production per plasma cell (day^-1)
  kdeg_APA    = 0.008, // APA degradation (day^-1)

  // ---- Pituitary Dynamics ----
  kpit_inflam = 0.015, // pituitary inflammation rate (per immune unit)
  kpit_repair = 0.004, // pituitary repair rate (day^-1)
  kpit_fibs   = 0.002, // fibrosis accumulation (day^-1)

  // ---- HPA Axis ----
  kprod_ACTH  = 0.50,  // ACTH production baseline (pg/mL/h)
  kdeg_ACTH   = 0.35,  // ACTH degradation (h^-1), t1/2 ~2h
  kprod_Cort  = 2.00,  // cortisol production per ACTH unit (nmol/L/h per pg/mL)
  kdeg_Cort   = 0.12,  // cortisol degradation (h^-1), t1/2 ~1.5h
  IC50_Cort   = 400,   // cortisol IC50 for ACTH feedback (nmol/L)
  n_Hill_Cort = 2.0,   // Hill coefficient for cortisol feedback
  Emax_Cort   = 0.90,  // maximal cortisol feedback on ACTH

  // ---- HPT Axis ----
  kprod_TSH   = 0.80,  // TSH production baseline (mIU/L/h)
  kdeg_TSH    = 0.15,  // TSH degradation (h^-1), t1/2 ~1h in disease
  kprod_fT4   = 0.50,  // fT4 production per TSH (pmol/L/h per mIU/L)
  kdeg_fT4    = 0.06,  // fT4 degradation (h^-1), t1/2 ~7 days
  IC50_fT4    = 15.0,  // fT4 IC50 for TSH feedback (pmol/L)
  Emax_fT4    = 0.85,  // maximal fT4 feedback on TSH

  // ---- GH/IGF-1 Axis ----
  kprod_GH    = 1.50,  // GH production baseline (ng/mL/h, pulsatile)
  kdeg_GH     = 0.40,  // GH degradation (h^-1), t1/2 ~20 min
  kprod_IGF1  = 0.30,  // IGF-1 production per GH (ng/mL/h per ng/mL)
  kdeg_IGF1   = 0.045, // IGF-1 degradation (h^-1), t1/2 ~15h
  IC50_IGF1   = 200,   // IGF-1 IC50 for GH feedback (ng/mL)

  // ---- HPG Axis (simplified) ----
  kprod_FSH   = 0.30,  // FSH production baseline (IU/L/h)
  kdeg_FSH    = 0.10,  // FSH degradation (h^-1)
  kprod_LH    = 0.40,  // LH production baseline (IU/L/h)
  kdeg_LH     = 0.20,  // LH degradation (h^-1)
  kprod_E2    = 5.0,   // Estradiol production per FSH (pmol/L/h per IU/L)
  kdeg_E2     = 0.08,  // Estradiol degradation (h^-1)
  IC50_E2     = 200,   // E2 IC50 for FSH/LH feedback (pmol/L)

  // ---- Prolactin ----
  kprod_PRL   = 0.50,  // PRL production baseline (ng/mL/h)
  kdeg_PRL    = 0.12,  // PRL degradation (h^-1)
  k_stalk_PRL = 0.30,  // stalk effect coefficient (elevates PRL)

  // ---- ADH (posterior pituitary) ----
  kprod_ADH   = 0.40,  // ADH production baseline (pg/mL/h)
  kdeg_ADH    = 0.20,  // ADH degradation (h^-1)
  k_post_loss = 0.80,  // posterior involvement coefficient (0=none, 1=severe)

  // ---- Drug Effect Parameters ----
  EC50_pred   = 0.50,  // prednisolone EC50 for immune suppression (mg/L) [out of scope]
  Emax_pred   = 0.85,  // max prednisolone effect (0-1) [out of scope]
  EC50_AZA    = 0.20,  // azathioprine EC50 for T/B suppression (mg/L equiv.) [was EC50_aza]
  EMAX_AZA    = 0.70,  // max azathioprine immune effect [was Emax_aza]
  GAMMA_AZA   = 1.0,   // Hill coefficient (new -- original had no explicit
                       // exponent; pow(x,1)=x, rename not a refit)
  EC50_RTX    = 0.05,  // rituximab EC50 for B-cell depletion (mg/L) [was EC50_rtx]
  EMAX_RTX    = 0.90,  // max rituximab B-cell depletion effect [was Emax_rtx]
  GAMMA_RTX   = 1.0,   // Hill coefficient (new -- same rationale as GAMMA_AZA)

  // ---- Baseline Steady-State Values ----
  ACTH_ss     = 22.0,  // pg/mL (normal 6-50 pg/mL)
  Cort_ss     = 500,   // nmol/L (normal 200-700 nmol/L morning)
  TSH_ss      = 2.5,   // mIU/L (normal 0.4-4.0 mIU/L)
  fT4_ss      = 15.0,  // pmol/L (normal 12-22 pmol/L)
  GH_ss       = 2.0,   // ng/mL (normal 0.5-5 ng/mL)
  IGF1_ss     = 180,   // ng/mL (normal 100-300 ng/mL)
  FSH_ss      = 5.0,   // IU/L
  LH_ss       = 5.0,   // IU/L
  E2_ss       = 200,   // pmol/L
  PRL_ss      = 15.0,  // ng/mL
  ADH_ss      = 3.0,   // pg/mL
  PitFunc_ss  = 1.0    // normalized pituitary function (0-1)

// [build-compat fix -- see lhyp_refactor_notes.md / UPSTREAM_ISSUES.md]
// The original declares all 25 compartments in a separate $CMT block AND
// re-declares every one of them in $INIT (NAME = value). mrgsolve 2.0.1
// treats $INIT as its own compartment-declaring block, not a companion to
// $CMT, so using both for the same names is rejected at build time with
// "Duplicated model names" -- reproduces identically from the untouched
// original, unrelated to AZA/RTX’s own PK. Settled-policy syntax-only fix
// (matching the identical defect class already fixed this way in
// lymphangioleiomyomatosis/lam_mrgsolve_model_refactored.R, issue #117):
// the $CMT block is deleted here and $INIT alone (below) declares all 25
// compartments, in the exact same order the original’s $CMT block used --
// declares no new compartment, changes no numeric value, and leaves
// 1-based dosing/output compartment indices unchanged.

$GLOBAL
// Pluggable-PK refactor (AZA, RTX only -- see file header and
// lhyp_refactor_notes.md). mrgsolve 2.0.1 hoists every block-local
// "double NAME" into one shared anonymous namespace, so the same name
// declared as a local in two different DSL blocks collides -- and $TABLE
// cannot read a $MAIN local at all. The single, normalized exposed
// concentration/effect for each compound is predeclared once here and
// assigned (no "double" keyword) at each of the original’s own
// evaluation sites below; never re-declared.
// Discoverability fix: C_AZA/C_RTX removed from this bare forward-declare.
// The $TABLE evaluation site (the one that determines what $CAPTURE
// reports) below is now a genuine `double C_<STEM> = <expr>;` initializing
// statement instead, so downstream tooling can pattern-match it as a
// single contiguous statement; the $MAIN evaluation site keeps its bare
// assignment, referring to the same file-scope member. EFFECT_AZA/
// EFFECT_RTX were not in scope for this fix and keep the original
// forward-declare mechanism.
double EFFECT_AZA, EFFECT_RTX;

$MAIN
  // ---- Drug concentrations ----
  double Cpred = Pred_central / Vc;       // mg/L [out of scope, untouched]
  // AZA/RTX: single normalized site (see $GLOBAL above). The original
  // computed this identical formula independently here in $MAIN (as
  // "Caza"/"Crtx") AND again in $TABLE (as "Caza_obs"/"Crtx_obs") --
  // redirected onto one shared name, recomputed identically at both of
  // the original’s own call sites (kept in $MAIN here, matching where
  // the original placed it, per the guide’s "keep a calculation in the
  // block the original used it in" rule).
  C_AZA = CENT_AZA / V1_AZA;    // mg/L [was: double Caza = AZA_plasma / Vc_aza;]
  C_RTX = CENT_RTX / V1_RTX;    // mg/L [was: double Crtx = RTX_plasma / Vc_rtx;]

  // ---- Drug effect functions ----
  double E_pred = Emax_pred * Cpred / (EC50_pred + Cpred);   // [out of scope, untouched]
  // EFFECT_AZA/EFFECT_RTX: canonical named Hill interface for each
  // compound -- a rename, not a refit. The original’s E_aza/E_rtx ratios
  // were already exactly this shape (implicit gamma=1, now GAMMA_AZA/
  // GAMMA_RTX = 1 explicitly).
  EFFECT_AZA = EMAX_AZA * pow(C_AZA, GAMMA_AZA) / (pow(EC50_AZA, GAMMA_AZA) + pow(C_AZA, GAMMA_AZA));
  EFFECT_RTX = EMAX_RTX * pow(C_RTX, GAMMA_RTX) / (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX));

  // ---- Combined immunosuppression ----
  double ImmunoSupp = 1.0 - (1.0 - (1.0 - E_pred)) * (1.0 - EFFECT_AZA * 0.5);

  // ---- Pituitary function (clamped 0–1) ----
  double Pfunc = (PitFunc > 1.0) ? 1.0 : (PitFunc < 0.0) ? 0.0 : PitFunc;
  double Pinfl = (PitInf  < 0.0) ? 0.0 : PitInf;

  // ---- HPA axis ----
  // Cortisol negative feedback on ACTH (Hill equation)
  double Cort_FB = 1.0 - Emax_Cort * pow(Cortisol, n_Hill_Cort) /
                         (pow(IC50_Cort, n_Hill_Cort) + pow(Cortisol, n_Hill_Cort));
  double ACTH_prod = kprod_ACTH * Cort_FB * Pfunc;

  // Prednisolone suppresses ACTH (glucocorticoid effect via GR)
  ACTH_prod = ACTH_prod * (1.0 - E_pred * 0.70);

  // ---- HPT axis ----
  double fT4_FB = 1.0 - Emax_fT4 * fT4 / (IC50_fT4 + fT4);
  double TSH_prod = kprod_TSH * fT4_FB * Pfunc;

  // ---- GH axis ----
  double IGF1_FB = 1.0 - 0.80 * IGF1 / (IC50_IGF1 + IGF1);
  double GH_prod = kprod_GH * IGF1_FB * Pfunc;

  // ---- HPG axis ----
  double E2_FB = 1.0 - 0.80 * E2 / (IC50_E2 + E2);
  double FSH_prod = kprod_FSH * E2_FB * Pfunc;
  double LH_prod  = kprod_LH  * E2_FB * Pfunc;

  // ---- Prolactin (stalk effect increases PRL) ----
  // Stalk involvement reduces dopamine → increases PRL
  double stalk_factor = 1.0 + k_stalk_PRL * Pinfl;
  double PRL_prod = kprod_PRL * stalk_factor;
  // Functional lactotrophs still needed for PRL
  double PRL_clearance = kdeg_PRL * Pfunc;

  // ---- ADH (posterior pituitary involvement) ----
  double ADH_prod = kprod_ADH * (1.0 - k_post_loss * Pinfl / (1.0 + Pinfl));

  // ---- Immune cell kinetics ----
  double Imm_act = (Te + APA * 0.01);   // aggregate immune activation driving pit damage

  double pit_damage = kpit_inflam * Imm_act * (1.0 - ImmunoSupp);
  double pit_recovery = kpit_repair * Pfunc;

$ODE
  // ---- Prednisolone PK (out of scope -- untouched) ----
  dxdt_Pred_gut     = -ka * Pred_gut;
  dxdt_Pred_central =  ka * Pred_gut
                     - (CL_P / Vc) * Pred_central
                     - (Q_P  / Vc) * Pred_central
                     + (Q_P  / Vp) * Pred_periph;
  dxdt_Pred_periph  =  (Q_P  / Vc) * Pred_central
                     - (Q_P  / Vp) * Pred_periph;

  // ---- Azathioprine PK (refactored: stem AZA) ----
  dxdt_GUT_AZA  = -KA_AZA * GUT_AZA;
  dxdt_CENT_AZA =  KA_AZA * GUT_AZA - (CL_AZA / V1_AZA) * CENT_AZA;

  // ---- Rituximab PK (IV, 1-compartment; refactored: stem RTX) ----
  dxdt_CENT_RTX = -(CL_RTX / V1_RTX) * CENT_RTX;

  // ---- Immune Cell Dynamics (time in days, convert: divide h params by 24) ----
  // Naive T cells (de novo thymic production + homeostatic proliferation)
  dxdt_Tn =   kprolif_Tn / 24.0 * Tn_baseline
            - kdeath_Tn  / 24.0 * Tn
            - kact_T     / 24.0 * Tn * (1.0 - ImmunoSupp);

  // Effector T cells
  dxdt_Te =   kact_T    / 24.0 * Tn * (1.0 - ImmunoSupp)
            - kdeath_Te / 24.0 * Te * (1.0 + E_pred * 2.5)  // GC-induced apoptosis
            - kinact_T  / 24.0 * Te;

  // Regulatory T cells (induced by TGF-β, enhanced by prednisolone)
  dxdt_Tr =   kinduct_Tr / 24.0 * Tn * (1.0 + E_pred * 1.5)
            + kprolif_Tr / 24.0 * Tr
            - kdeath_Tr  / 24.0 * Tr;

  // Naive B cells
  dxdt_Bn =   kprolif_Bn / 24.0 * Bn_baseline
            - kdeath_Bn  / 24.0 * Bn
            - kact_B     / 24.0 * Bn * (1.0 - ImmunoSupp) * (1.0 - EFFECT_RTX);

  // Plasma cells
  dxdt_Bp =   kact_B    / 24.0 * Bn * (1.0 - ImmunoSupp) * (1.0 - EFFECT_RTX)
            - kdeath_Bp / 24.0 * Bp * (1.0 + E_pred * 1.5);

  // Anti-pituitary antibodies
  dxdt_APA =   kprod_APA / 24.0 * Bp
             - kdeg_APA  / 24.0 * APA;

  // ---- Pituitary Dynamics ----
  // Inflammatory volume (MRI visible mass, relative)
  dxdt_PitInf  =   pit_damage / 24.0
               - kpit_fibs / 24.0 * Pinfl  // some turns to fibrosis
               - (kpit_repair + E_pred * 0.02) / 24.0 * Pinfl;

  // Functional pituitary cell mass (0–1)
  dxdt_PitFunc  = - pit_damage  / 24.0
                  + pit_recovery / 24.0;

  // ---- HPA Axis ----
  dxdt_ACTH    = ACTH_prod   - kdeg_ACTH * ACTH;
  dxdt_Cortisol = kprod_Cort * ACTH - kdeg_Cort * Cortisol;

  // ---- HPT Axis ----
  dxdt_TSH  = TSH_prod   - kdeg_TSH * TSH;
  dxdt_fT4  = kprod_fT4 * TSH - kdeg_fT4 * fT4;

  // ---- GH / IGF-1 Axis ----
  dxdt_GH   = GH_prod    - kdeg_GH * GH;
  dxdt_IGF1 = kprod_IGF1 * GH - kdeg_IGF1 * IGF1;

  // ---- HPG Axis ----
  dxdt_FSH  = FSH_prod   - kdeg_FSH * FSH;
  dxdt_LH   = LH_prod    - kdeg_LH  * LH;
  dxdt_E2   = kprod_E2 * FSH - kdeg_E2 * E2;

  // ---- Prolactin ----
  dxdt_PRL  = PRL_prod - kdeg_PRL * PRL;

  // ---- ADH ----
  dxdt_ADH  = ADH_prod - kdeg_ADH * ADH;

$INIT
  // Drug compartments start at zero
  Pred_gut = 0, Pred_central = 0, Pred_periph = 0,
  GUT_AZA  = 0, CENT_AZA     = 0,
  CENT_RTX = 0,

  // Immune cells at healthy steady state
  Tn  = 1000, Te = 50,   Tr  = 200,
  Bn  = 500,  Bp = 20,   APA = 10,

  // Pituitary: LyH starts with partial involvement
  PitInf  = 0.3,  // ~1 SD above normal (mild lymphocytic infiltration)
  PitFunc = 0.85, // 15% functional loss at baseline (early LyH)

  // Hormone axes (normal healthy steady state for disease onset context)
  ACTH     = 22.0,  Cortisol = 500,
  TSH      = 2.5,   fT4      = 15.0,
  GH       = 2.0,   IGF1     = 180,
  FSH      = 5.0,   LH       = 5.0,   E2  = 200,
  PRL      = 15.0,  ADH      = 3.0

$TABLE
  double Cpred_obs = Pred_central / Vc;   // prednisolone plasma (mg/L) [out of scope]
  // AZA/RTX: reassign the single normalized site here too (TABLE cannot
  // see the $MAIN assignment for this same record without risking a
  // stale read -- see $GLOBAL/refactor notes), identical formula to
  // $MAIN above, matching the original’s own Caza_obs/Crtx_obs
  // recomputation. Captured directly below instead of under a second,
  // separately-named "_obs" variable.
  // Discoverability fix: converted to genuine `double NAME = expr;`
  // initializing statements (identical formula/value) so this, the site
  // that actually feeds $CAPTURE, is literal-text-discoverable.
  double C_AZA = CENT_AZA / V1_AZA;
  double C_RTX = CENT_RTX / V1_RTX;

  // Pituitary function score (composite, 0=all deficient, 1=normal)
  double PFS = (PitFunc < 0) ? 0 : (PitFunc > 1) ? 1 : PitFunc;

  // Cortisol/ACTH ratio (indicator of adrenal reserve)
  double CA_ratio = (ACTH > 0.1) ? Cortisol / ACTH : 999;

  // Relative hormone levels vs baseline
  double ACTH_pct  = ACTH    / 22.0   * 100;
  double Cort_pct  = Cortisol/ 500.0  * 100;
  double TSH_pct   = TSH     / 2.5    * 100;
  double fT4_pct   = fT4     / 15.0   * 100;
  double GH_pct    = GH      / 2.0    * 100;
  double IGF1_pct  = IGF1    / 180.0  * 100;
  double FSH_pct   = FSH     / 5.0    * 100;
  double PRL_ratio = PRL     / 15.0;  // >3 = clinically significant elevation

  double ADA_score = (ACTH < 10) ? 1 : 0;  // secondary adrenal insufficiency flag

  // [build-compat fix -- see lhyp_refactor_notes.md / UPSTREAM_ISSUES.md]
  // The original ends $TABLE with three bare, unheadered "capture NAME1
  // NAME2 ..." lines (lowercase, no preceding $CAPTURE marker), which
  // mrgsolve 2.0.1 tries to compile as a call to an undeclared C++
  // function "capture(...)" and rejects -- the same defect class already
  // logged for bronchiectasis/bex_mrgsolve_model.R (UPSTREAM_ISSUES.md
  // #46). Settled-policy syntax-only fix: each bare "capture" is given a
  // real "$CAPTURE" header (one block per original line, same names, same
  // order) -- no symbol added, removed, or reordered.
  $CAPTURE Cpred_obs PFS CA_ratio
  $CAPTURE ACTH_pct Cort_pct TSH_pct fT4_pct GH_pct IGF1_pct FSH_pct
  $CAPTURE PRL_ratio ADA_score C_AZA EFFECT_AZA C_RTX EFFECT_RTX
'

## ---- Compile the model -------------------------------------------------------
mod <- mcode("LymphoHypophysitis", code, quiet = TRUE)

## ---- Helper: single dose event builder --------------------------------------
## Prednisolone: oral, BID dosing (out of scope, untouched)
pred_events <- function(daily_mg, duration_days, start_day = 0) {
  dose_per_admin <- daily_mg / 2
  times <- seq(start_day * 24, (start_day + duration_days) * 24 - 12, by = 12)
  ev(ID = 1, time = times, amt = dose_per_admin * 0.82,
     cmt = "Pred_gut", evid = 1)
}

## Azathioprine: oral, OD
aza_events <- function(daily_mg, duration_days, start_day = 30) {
  times <- seq(start_day * 24, (start_day + duration_days) * 24 - 24, by = 24)
  ev(ID = 1, time = times, amt = daily_mg * 0.50,
     cmt = "GUT_AZA", evid = 1)
}

## Rituximab: IV bolus every 6 months
rtx_events <- function(dose_mg, n_infusions = 2, interval_days = 180) {
  times <- seq(0, (n_infusions - 1) * interval_days * 24, by = interval_days * 24)
  ev(ID = 1, time = times, amt = dose_mg,
     cmt = "CENT_RTX", evid = 1)
}

## ---- Simulation parameters --------------------------------------------------
sim_days <- 730   # 2 years
sim_times <- seq(0, sim_days * 24, by = 4)  # hourly output every 4 h

## ============================================================================
## SCENARIO 1: No treatment (natural history)
## ============================================================================
out_s1 <- mod %>%
  mrgsim(end = sim_days * 24, delta = 4) %>%
  as.data.frame() %>%
  mutate(Scenario = "S1: No Treatment (Natural History)",
         time_d = time / 24)

## ============================================================================
## SCENARIO 2: High-dose prednisolone (60 mg/day × 4 weeks),
##             then taper (40→30→20→10 mg/day, 2 weeks each)
## ============================================================================
e_s2 <- pred_events(60,  28,  0) +
         pred_events(40,  14,  28) +
         pred_events(30,  14,  42) +
         pred_events(20,  14,  56) +
         pred_events(10,  14,  70) +
         pred_events(5,   180, 84)   # maintenance 5 mg/day for 6 months

out_s2 <- mod %>%
  mrgsim(events = e_s2, end = sim_days * 24, delta = 4) %>%
  as.data.frame() %>%
  mutate(Scenario = "S2: High-Dose Prednisolone + Taper",
         time_d = time / 24)

## ============================================================================
## SCENARIO 3: Prednisolone + Azathioprine (steroid-sparing)
## ============================================================================
e_s3 <- pred_events(60, 28, 0) +
         pred_events(40, 14, 28) +
         pred_events(20, 14, 42) +
         pred_events(10, 30, 56) +
         pred_events(5, 500, 86) +      # low-dose maintenance
         aza_events(150, 600, 30)       # AZA 150 mg/day starting month 1

out_s3 <- mod %>%
  mrgsim(events = e_s3, end = sim_days * 24, delta = 4) %>%
  as.data.frame() %>%
  mutate(Scenario = "S3: Pred + Azathioprine (Steroid-Sparing)",
         time_d = time / 24)

## ============================================================================
## SCENARIO 4: Rituximab (checkpoint inhibitor-associated LyH)
## ============================================================================
e_s4 <- rtx_events(1000, n_infusions = 2, interval_days = 180) +
         pred_events(20, 30, 0)    # short-course prednisolone during first infusion

out_s4 <- mod %>%
  mrgsim(events = e_s4, end = sim_days * 24, delta = 4) %>%
  as.data.frame() %>%
  mutate(Scenario = "S4: Rituximab (Anti-CD20) Therapy",
         time_d = time / 24)

## ============================================================================
## SCENARIO 5: Standard low-dose prednisolone only (maintenance)
##             — suboptimal treatment
## ============================================================================
e_s5 <- pred_events(10, 730, 0)  # 10 mg/day throughout

out_s5 <- mod %>%
  mrgsim(events = e_s5, end = sim_days * 24, delta = 4) %>%
  as.data.frame() %>%
  mutate(Scenario = "S5: Low-Dose Prednisolone Maintenance",
         time_d = time / 24)

## ============================================================================
## COMBINED RESULTS
## ============================================================================
all_scenarios <- bind_rows(out_s1, out_s2, out_s3, out_s4, out_s5)

## ---- Publication-ready plots ------------------------------------------------

theme_qsp <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#1a1a2e", color = NA),
    panel.background = element_rect(fill = "#16213e", color = NA),
    panel.grid.major = element_line(color = "#2d3561", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    text             = element_text(color = "#eaeaea"),
    axis.text        = element_text(color = "#aaaaaa"),
    strip.text       = element_text(color = "#80deea", face = "bold"),
    legend.background = element_rect(fill = "#1a1a2e", color = NA),
    legend.text      = element_text(color = "#eaeaea"),
    plot.title       = element_text(color = "#80deea", face = "bold", size = 13)
  )

scenario_colors <- c(
  "S1: No Treatment (Natural History)"       = "#ef5350",
  "S2: High-Dose Prednisolone + Taper"        = "#42a5f5",
  "S3: Pred + Azathioprine (Steroid-Sparing)" = "#66bb6a",
  "S4: Rituximab (Anti-CD20) Therapy"         = "#ffa726",
  "S5: Low-Dose Prednisolone Maintenance"     = "#ab47bc"
)

## Plot 1: Pituitary Function Score
p1 <- all_scenarios %>%
  ggplot(aes(x = time_d, y = PFS * 100, color = Scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Pituitary Function Score Over Time",
       x = "Time (days)", y = "Pituitary Function (%)") +
  geom_hline(yintercept = 70, linetype = "dashed", color = "yellow", alpha = 0.5) +
  annotate("text", x = 10, y = 68, label = "Hypopituitarism threshold", color = "yellow", size = 3) +
  theme_qsp

## Plot 2: HPA Axis (ACTH and Cortisol)
p2 <- all_scenarios %>%
  select(time_d, Scenario, ACTH, Cortisol) %>%
  pivot_longer(c(ACTH, Cortisol), names_to = "hormone", values_to = "level") %>%
  ggplot(aes(x = time_d, y = level, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~hormone, scales = "free_y",
             labeller = labeller(hormone = c(ACTH = "ACTH (pg/mL)", Cortisol = "Cortisol (nmol/L)"))) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "HPA Axis Dynamics", x = "Time (days)", y = "Concentration") +
  theme_qsp

## Plot 3: HPT and GH Axes
p3 <- all_scenarios %>%
  select(time_d, Scenario, TSH, fT4, GH, IGF1) %>%
  pivot_longer(c(TSH, fT4, GH, IGF1), names_to = "hormone", values_to = "level") %>%
  ggplot(aes(x = time_d, y = level, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~hormone, scales = "free_y") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "HPT and GH/IGF-1 Axes", x = "Time (days)", y = "Level") +
  theme_qsp

## Plot 4: Immune dynamics and APA
p4 <- all_scenarios %>%
  select(time_d, Scenario, Te, Tr, APA, Bp) %>%
  pivot_longer(c(Te, Tr, APA, Bp), names_to = "cell", values_to = "level") %>%
  ggplot(aes(x = time_d, y = level, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~cell, scales = "free_y",
             labeller = labeller(cell = c(Te = "Effector T cells (AU)",
                                          Tr = "Regulatory T cells (AU)",
                                          APA = "Anti-Pituitary Abs (AU)",
                                          Bp = "Plasma cells (AU)"))) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Immune Cell and Antibody Dynamics",
       x = "Time (days)", y = "Level") +
  theme_qsp

## Plot 5: Prolactin and ADH (posterior involvement)
p5 <- all_scenarios %>%
  select(time_d, Scenario, PRL, ADH, FSH, LH) %>%
  pivot_longer(c(PRL, ADH, FSH, LH), names_to = "hormone", values_to = "level") %>%
  ggplot(aes(x = time_d, y = level, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~hormone, scales = "free_y") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Prolactin, ADH, FSH, LH Dynamics",
       x = "Time (days)", y = "Level") +
  theme_qsp

## Plot 6: Drug concentrations
p6 <- all_scenarios %>%
  filter(Scenario %in% c("S2: High-Dose Prednisolone + Taper",
                          "S3: Pred + Azathioprine (Steroid-Sparing)",
                          "S4: Rituximab (Anti-CD20) Therapy")) %>%
  select(time_d, Scenario, Cpred_obs, C_AZA, C_RTX) %>%
  pivot_longer(c(Cpred_obs, C_AZA, C_RTX), names_to = "drug", values_to = "conc") %>%
  mutate(drug = recode(drug,
                        Cpred_obs = "Prednisolone (mg/L)",
                        C_AZA     = "Azathioprine (mg/L)",
                        C_RTX     = "Rituximab (mg/L)")) %>%
  ggplot(aes(x = time_d, y = conc, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~drug, scales = "free_y") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Drug Plasma Concentrations (PK)",
       x = "Time (days)", y = "Concentration") +
  theme_qsp

## ---- Summary statistics at key timepoints ------------------------------------
summary_table <- all_scenarios %>%
  filter(time_d %in% c(0, 30, 90, 180, 365, 730)) %>%
  group_by(Scenario, time_d) %>%
  summarise(
    Pit_Func_pct  = round(mean(PFS * 100), 1),
    ACTH_pgmL     = round(mean(ACTH), 1),
    Cortisol_nmolL = round(mean(Cortisol), 0),
    TSH_mIUL      = round(mean(TSH), 2),
    fT4_pmolL     = round(mean(fT4), 1),
    GH_ngmL       = round(mean(GH), 2),
    IGF1_ngmL     = round(mean(IGF1), 0),
    PRL_ngmL      = round(mean(PRL), 1),
    ADH_pgmL      = round(mean(ADH), 2),
    APA_AU        = round(mean(APA), 1),
    .groups = "drop"
  )

print(summary_table)

## ---- Print plots -------------------------------------------------------------
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

message("\n=== Lymphocytic Hypophysitis QSP Model (PK/PD refactored) ===")
message("Compiled with mrgsolve, 22 ODE compartments")
message("5 treatment scenarios simulated over 2 years")
message("Outputs: pituitary function, hormone axes, immune dynamics, drug PK")
message("See lhyp_shiny_app.R for interactive dashboard")
message("See lhyp_refactor_notes.md for the PK/PD refactor (AZA, RTX) detail")
