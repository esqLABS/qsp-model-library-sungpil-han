## =============================================================================
## Cervical Cancer QSP Model — mrgsolve ODE Implementation (PK/PD REFACTORED)
## HPV-driven Cervical Squamous Cell Carcinoma / Adenocarcinoma
## =============================================================================
## [refactor] Per FORK_WORKFLOW_GUIDE.md Part 2 (pluggable PK, named Hill
## interface). Five compounds refactored, per the existing rows in
## driver-patches/data/compound_perturbation_census.md, all classified
## "Redirect concentration (clean single site)": Bevacizumab (BEV), Cisplatin
## (CIS), Paclitaxel (PAC), Pembrolizumab (PEM), and TV-ADC / tisotumab
## vedotin (TVADC). No other compounds exist in this file. Full rationale,
## per-compound archetype, and verification results are in
## cc_refactor_notes.md.
##
## Renamed (values unchanged from the original):
##   CIS_C1/CIS_C2                 -> CENT_CIS/PERI_CIS
##   PAC_C1/PAC_C2                 -> CENT_PAC/PERI_PAC
##   BEV_C1/BEV_C2                 -> CENT_BEV/PERI_BEV
##   PEMBRO_C1/PEMBRO_C2           -> CENT_PEM/PERI_PEM (harmonized to the
##     "PEM" stem the original's own $PARAM block already used --
##     CL_PEM/V1_PEM/Q_PEM/V2_PEM -- resolving the original's own PEM- vs
##     PEMBRO- naming inconsistency, exactly the chaos this refactor removes)
##   TV_ADC_C1/TV_ADC_C2           -> CENT_TVADC/PERI_TVADC
##   CL_TV/V1_TV/Q_TV/V2_TV        -> CL_TVADC/V1_TVADC/Q_TVADC/V2_TVADC
##     (renamed off the bare "TV" stem, which collided in meaning -- not in
##     symbol name -- with the unrelated TV = tumor-volume compartment)
##   MMAE_free                     -> MMAE_TVADC
##   k_dec_MMAE                    -> k_dec_MMAE_TVADC
##   Pt_DNA                        -> ADDUCT_CIS
##   BEV_effect                    -> VEGF_bind_BEV (rename only, same formula)
##   eff_ICI                       -> EFFECT_PEM (rename only, same formula)
##   pac_eff, adc_eff              -> EFFECT_PAC, EFFECT_TVADC (Hill form made
##     explicit; see below)
## New (not in original; named Hill-interface parameters that make explicit a
## shape the original already had implicitly as a plain ratio -- rename, not
## a refit; also FR_MMAE_TVADC, naming a literal 0.4 fraction that was
## hardcoded in the original's dxdt_MMAE_free line):
##   EMAX_PAC=1, EC50_PAC=100.0, GAMMA_PAC=1     (was PAC_C1/(PAC_C1+100.0))
##   EMAX_BEV=1, EC50_BEV=10.0,  GAMMA_BEV=1     (was BEV_C1/(BEV_C1+10.0))
##   EMAX_TVADC=1, EC50_TVADC=1.0, GAMMA_TVADC=1 (was MMAE_free/(MMAE_free+1))
##   FR_MMAE_TVADC = 0.4
##
## Exposed concentration convention note: this original's own dosing already
## divides the bolus dose by V1 before adding it to the central compartment
## (e.g. dose_cisplatin: amt = dose_mg/V1), and its own $TABLE captured the
## central compartment directly as "the concentration" (CIS_Conc = CIS_C1,
## no further /V1). So C_<STEM> here is a direct alias of CENT_<STEM>, *not*
## CENT_<STEM>/V1_<STEM> as the guide's canonical template shows -- dividing
## again would silently change the original's numeric behavior. This is
## disclosed in full in cc_refactor_notes.md.
##
## No TMDD anywhere in this file: Bevacizumab and Pembrolizumab (the two
## compounds explicitly checked) are both plain linear 2-compartment PK with
## no receptor/complex compartment -- confirmed against the actual equations,
## not assumed from real-world biology. TV-ADC is bespoke (archetype 2 base
## + an extra, non-archetypal free-payload (MMAE) compartment) -- see notes.
## =============================================================================
##
## Compartments (19 ODEs, unchanged count/roles from the original -- renamed only):
##   1  CENT_CIS   — Cisplatin central compartment
##   2  PERI_CIS   — Cisplatin peripheral compartment
##   3  CENT_PAC   — Paclitaxel central compartment (recurrent/interval CT)
##   4  PERI_PAC   — Paclitaxel peripheral compartment
##   5  CENT_BEV   — Bevacizumab central compartment
##   6  PERI_BEV   — Bevacizumab peripheral compartment
##   7  CENT_PEM   — Pembrolizumab central compartment
##   8  PERI_PEM   — Pembrolizumab peripheral compartment
##   9  CENT_TVADC — Tisotumab vedotin (ADC) central compartment
##  10  PERI_TVADC — Tisotumab vedotin (ADC) peripheral compartment
##  11  MMAE_TVADC — Free intratumoral MMAE payload (relative)
##  12  VEGF       — Free VEGF-A (ng/mL)
##  13  ADDUCT_CIS — Platinum-DNA adduct burden (relative, 0-1)
##  14  RT_SF      — Cumulative radiation surviving fraction (log-scale damage)
##  15  TV         — Tumor volume (cm³, Gompertz growth + multi-modal kill)
##  16  SCCAg      — SCC-Ag serum biomarker (ng/mL)
##  17  HPVload    — HPV viral load (relative, log copies)
##  18  CD8T       — CD8+ T effector cells (relative)
##  19  PDL1_exp   — Tumor PD-L1 expression (relative, CPS-like)
##
## Key References (calibration): unchanged from the original -- see
## cc_mrgsolve_model.R's header for the full citation list (Rose 1999 NEJM
## GOG-120, Green 2001 Lancet, Reece 1987, Eifel 2004 JCO, Fowler 1989,
## Pötter 2018 EMBRACE, Tewari 2014 NEJM GOG-240 PMID 24552320, Lorusso 2024
## Lancet / Monk 2023 JCO KEYNOTE-A18, Colombo 2021 NEJM KEYNOTE-826 PMID
## 34534429, Ahamadi 2017 CPT:PSP, Coleman 2021 Lancet Oncol innovaTV 204,
## Vergote 2024 JCO/NEJM Evid innovaTV 301, Lu 2008, Gaarenstroom 2000,
## Kasibhatla 2007).
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---------------------------------------------------------------
## mrgsolve model specification
## ---------------------------------------------------------------
code_cc <- '
$PROB Cervical Cancer QSP Model (PK/PD refactored)
Cisplatin-based CCRT +/- Pembrolizumab (KEYNOTE-A18) +/- Bevacizumab (GOG-240)
+/- Tisotumab Vedotin ADC (innovaTV 301) for recurrent/metastatic disease
[refactor] CIS/PAC/BEV/PEM/TVADC PK+PD blocks renamed to the forks
plumbable-PK convention (C_<STEM> exposed concentration, EFFECT_<STEM> named
disease effect). No other compounds exist in this file. See
cc_refactor_notes.md for full rationale and verification.

$PARAM @annotated
// ── Cisplatin PK (2-compartment, weekly during CCRT) ─────────────
CL_CIS   : 30.0  : Cisplatin (total Pt) clearance (L/h; Reece 1987)
V1_CIS   : 15.0  : Central Vd (L)
Q_CIS    : 20.0  : Inter-compartmental clearance (L/h)
V2_CIS   : 30.0  : Peripheral Vd (L)

// ── Paclitaxel PK (2-compartment, recurrent interval chemo) ──────
CL_PAC   : 13.2  : Paclitaxel total CL (L/h)
V1_PAC   : 6.5   : Central Vd (L)
Q_PAC    : 7.0   : Inter-comp CL (L/h)
V2_PAC   : 113.0 : Peripheral Vd (L)

// ── Bevacizumab PK (2-compartment, IV, GOG-240 15mg/kg q3w) ──────
CL_BEV   : 0.207 : Bevacizumab clearance (L/day; Lu 2008)
V1_BEV   : 2.91  : Central Vd (L)
Q_BEV    : 0.469 : Inter-comp CL (L/day)
V2_BEV   : 1.91  : Peripheral Vd (L)

// ── Pembrolizumab PK (2-compartment, IV, 200mg q3w) ──────────────
CL_PEM   : 0.213 : Pembrolizumab clearance (L/day; Ahamadi 2017)
V1_PEM   : 3.34  : Central Vd (L)
Q_PEM    : 0.638 : Inter-comp CL (L/day)
V2_PEM   : 2.68  : Peripheral Vd (L)

// ── Tisotumab vedotin ADC PK (2-compartment, IV, 2.0 mg/kg q3w) ──
CL_TVADC : 0.51  : TV-ADC clearance (L/day; conjugate)
V1_TVADC : 3.35  : Central Vd (L)
Q_TVADC  : 0.55  : Inter-comp CL (L/day)
V2_TVADC : 2.75  : Peripheral Vd (L)
k_dec_MMAE_TVADC: 5.0  : MMAE release/decay rate (1/day)
FR_MMAE_TVADC: 0.4     : Fraction of conjugate clearance flux released as free MMAE payload

// ── VEGF dynamics ─────────────────────────────────────────────────
VEGF0    : 0.20  : Baseline free VEGF-A (ng/mL)
ksyn_VEGF: 1.8   : VEGF production from tumor (scaled to TV)
kdeg_VEGF: 8.0   : VEGF degradation rate (1/day)
kbind_BEV: 45.0  : Bevacizumab-VEGF binding rate

// ── Platinum-DNA adducts ──────────────────────────────────────────
k_adduct : 0.02  : Rate of Pt-DNA adduct formation (1/(µg/mL·h))
k_repair : 0.10  : Adduct repair rate (1/h; NER activity)

// ── Radiotherapy LQ model ─────────────────────────────────────────
alpha_rt : 0.30  : LQ alpha (1/Gy; cervix SCC, Fowler 1989)
beta_rt  : 0.030 : LQ beta (1/Gy^2); alpha/beta = 10 Gy
k_radiosens: 1.6 : Cisplatin radiosensitization multiplier on alpha
k_reoxy  : 0.05  : Reoxygenation recovery rate between fractions (1/day)
k_repop  : 0.008 : Accelerated repopulation growth offset (1/day, from day 21)

// ── Tumor growth (Gompertz model) ─────────────────────────────────
TV0      : 40.0  : Initial tumor volume (cm³; FIGO IIB-IIIB bulky)
kg       : 0.010 : Gompertz growth rate (1/day)
TV_max   : 1500.0: Carrying capacity (cm³)
k_kill_RT: 1.0   : Radiation kill scaling (per unit cumulative damage)
k_kill_Pt: 0.005 : Platinum kill rate constant (1/(relative adduct·day))
k_kill_pac: 0.003: Paclitaxel kill rate constant (1/day, Imax-scaled)
k_kill_ICI: 0.004: Pembrolizumab-enhanced CD8 kill rate (1/day)
k_kill_ADC: 0.006: Tisotumab vedotin kill rate (1/day, Imax-scaled)
k_kill_bev: 0.0015: Anti-angiogenic growth-inhibitory contribution (1/day)

// ── Paclitaxel Hill effect (named; rename of a plain ratio) ────────
EMAX_PAC : 1.0   : Max paclitaxel cytotoxic effect fraction (implicit in original)
EC50_PAC : 100.0 : EC50 for paclitaxel cytotoxic effect (ng/mL)
GAMMA_PAC: 1.0   : Hill coefficient (original had no explicit Hill term)

// ── Bevacizumab Hill effect (named; rename of a plain ratio) ───────
EMAX_BEV : 1.0   : Max anti-angiogenic growth-inhibitory effect fraction
EC50_BEV : 10.0  : EC50 for anti-angiogenic effect (mg/L)
GAMMA_BEV: 1.0   : Hill coefficient (original had no explicit Hill term)

// ── TV-ADC payload Hill effect (named; rename of a plain ratio) ────
EMAX_TVADC : 1.0 : Max ADC payload cytotoxic effect fraction
EC50_TVADC : 1.0 : EC50 for free MMAE payload effect (relative)
GAMMA_TVADC: 1.0 : Hill coefficient (original had no explicit Hill term)

// ── SCC-Ag dynamics ───────────────────────────────────────────────
// [refactor/build-fix] SCCAg_0 -> SCCAg0: renamed only to avoid colliding
// with the auto-generated "<CMT>_0" initial-value idiom used for the
// SCCAg compartment below (see $MAIN and cc_refactor_notes.md). This
// parameter is otherwise unused in $ODE/$MAIN/$TABLE, same as in the
// original -- rename only, no numeric or behavioral change.
SCCAg0   : 8.0   : Baseline SCC-Ag (ng/mL; advanced SCC, nl <1.5)
ksyn_SCC : 0.10  : SCC-Ag production per unit tumor (ng/mL/cm³/day)
kdeg_SCC : 0.25  : SCC-Ag degradation rate (1/day; t½≈2.8 days)

// ── HPV viral load dynamics ───────────────────────────────────────
// [refactor/build-fix] HPVload_0 -> HPVload0: same reason as SCCAg0 above
// (collision with the HPVload compartment’s "<CMT>_0" init idiom).
HPVload0 : 5.0   : Baseline HPV viral load (log10 copies, relative)
k_HPV_prod: 0.08 : HPV production proportional to tumor (1/day)
k_HPV_clear: 0.05: HPV clearance rate (1/day, immune + treatment-driven)
k_HPV_ICI : 2.0  : ICI-boosted clearance multiplier

// ── CD8+ T cell dynamics ──────────────────────────────────────────
// [refactor/build-fix] CD8T_0 -> CD8T0: same reason as SCCAg0 above
// (collision with the CD8T compartment’s "<CMT>_0" init idiom).
CD8T0    : 1.0   : Baseline CD8+ T (relative)
k_CD8_in : 0.12  : CD8+ influx rate (1/day)
k_CD8_out: 0.10  : CD8+ efflux rate (1/day)
k_exhaust: 0.35  : T cell exhaustion by tumor load & PD-L1
k_ICI    : 2.5   : ICI boost to CD8+ (fold increase)

// ── PD-L1 expression dynamics ─────────────────────────────────────
PDL1_0   : 1.0   : Baseline tumor PD-L1 expression (relative, CPS-like)
k_PDL1_up: 0.02  : IFN-gamma-driven adaptive PD-L1 upregulation (1/day)
PDL1_max : 3.0   : Maximum relative PD-L1 expression

// ── Scenario flags (0=off, 1=on) ──────────────────────────────────
ICI_flag  : 0    : Pembrolizumab present (0/1)
RT_flag   : 0    : Concurrent EBRT/brachytherapy active (0/1)
CPS_high  : 1    : PD-L1 CPS >= 1 status (1=eligible, 0=low)

$CMT @annotated
CENT_CIS  : Cisplatin central (µg/mL)
PERI_CIS  : Cisplatin peripheral (µg/mL)
CENT_PAC  : Paclitaxel central (ng/mL)
PERI_PAC  : Paclitaxel peripheral (ng/mL)
CENT_BEV  : Bevacizumab central (mg/L)
PERI_BEV  : Bevacizumab peripheral (mg/L)
CENT_PEM  : Pembrolizumab central (µg/mL)
PERI_PEM  : Pembrolizumab peripheral (µg/mL)
CENT_TVADC: Tisotumab vedotin central (µg/mL)
PERI_TVADC: Tisotumab vedotin peripheral (µg/mL)
MMAE_TVADC: Free intratumoral MMAE (relative)
VEGF      : Free VEGF-A (ng/mL)
ADDUCT_CIS: Platinum-DNA adducts (relative 0-1)
RT_SF     : Cumulative radiation damage (-log survFrac, relative)
TV        : Tumor volume (cm³)
SCCAg     : SCC-Ag serum (ng/mL)
HPVload   : HPV viral load (relative log10 copies)
CD8T      : CD8+ T cell relative level
PDL1_exp  : Tumor PD-L1 expression (relative)

$GLOBAL
// [refactor] exposed concentration / named-effect doubles for the five
// refactored compounds, predeclared here rather than in $PARAM: mrgsolve
// 2.0.1 compiles $PARAM members as read-only references inside $ODE, so a
// value that must be recomputed every timestep from state cannot also live
// in $PARAM (same constraint documented in the breast-cancer, AMD, and
// membranous-nephropathy refactors). These are visible in every
// simulation’s output via $CAPTURE and in /model_manifest’s outputPaths.
double C_CIS, C_PAC, C_BEV, C_PEM, C_TVADC;
double EFFECT_CIS, EFFECT_PAC, EFFECT_BEV, EFFECT_PEM, EFFECT_TVADC;

$MAIN
// [refactor] eff_ICI -> EFFECT_PEM (rename only, identical formula).
// Pembrolizumabs disease effect in this model is a binary treatment-flag
// switch (ICI_flag x CPS_high) -- it is NOT driven by its own PK
// concentration/occupancy (C_PEM is exposed below but not consumed by any
// disease equation). This is a property of the original, preserved as-is,
// not a new limitation introduced by the refactor -- see cc_refactor_notes.md.
EFFECT_PEM = ICI_flag * CPS_high;

// [build-fix] Pre-existing upstream defect: the original declares $CMT
// (compartment names) *and* a separate $INIT block assigning starting
// values to the same names -- mrgsolve 2.0.1 treats $INIT as its own
// compartment-declaring block, so using both together redeclares every
// compartment twice ("Duplicated model names"). Fixed here, syntax-only,
// by deleting $INIT and moving its assignments to the modern "<CMT>_0"
// idiom below (values unchanged). Disclosed in full, with the exact error
// text, in cc_refactor_notes.md and logged as a new UPSTREAM_ISSUES.md
// entry; the checked-in original still has $CMT+$INIT and does not build.
CENT_CIS_0=0; PERI_CIS_0=0;
CENT_PAC_0=0; PERI_PAC_0=0;
CENT_BEV_0=0; PERI_BEV_0=0;
CENT_PEM_0=0; PERI_PEM_0=0;
CENT_TVADC_0=0; PERI_TVADC_0=0; MMAE_TVADC_0=0;
VEGF_0=0.20;
ADDUCT_CIS_0=0;
RT_SF_0=0;
TV_0=40;
SCCAg_0=8.0;
HPVload_0=5.0;
CD8T_0=1.0;
PDL1_exp_0=1.0;

$ODE
// ── Cisplatin 2-compartment (archetype 2; renamed only) ────────────
// [refactor] C_CIS is a direct alias of CENT_CIS, not CENT_CIS/V1_CIS --
// the original already doses amt=dose_mg/V1 and captures the central
// compartment directly as "the concentration" with no further division;
// see header note and cc_refactor_notes.md.
C_CIS = CENT_CIS;
dxdt_CENT_CIS = -(CL_CIS/V1_CIS)*CENT_CIS - (Q_CIS/V1_CIS)*CENT_CIS
               + (Q_CIS/V2_CIS)*PERI_CIS;
dxdt_PERI_CIS = (Q_CIS/V1_CIS)*CENT_CIS - (Q_CIS/V2_CIS)*PERI_CIS;

// ── Paclitaxel 2-compartment (archetype 2; renamed only) ────────────
C_PAC = CENT_PAC;
dxdt_CENT_PAC = -(CL_PAC/V1_PAC)*CENT_PAC - (Q_PAC/V1_PAC)*CENT_PAC
               + (Q_PAC/V2_PAC)*PERI_PAC;
dxdt_PERI_PAC = (Q_PAC/V1_PAC)*CENT_PAC - (Q_PAC/V2_PAC)*PERI_PAC;

// ── Bevacizumab 2-compartment (archetype 2; renamed only) ───────────
C_BEV = CENT_BEV;
double VEGF_bind_BEV = kbind_BEV * C_BEV * VEGF;
dxdt_CENT_BEV = -(CL_BEV/V1_BEV)*CENT_BEV - (Q_BEV/V1_BEV)*CENT_BEV
               + (Q_BEV/V2_BEV)*PERI_BEV;
dxdt_PERI_BEV = (Q_BEV/V1_BEV)*CENT_BEV - (Q_BEV/V2_BEV)*PERI_BEV;

// ── Pembrolizumab 2-compartment (archetype 2; renamed only) ─────────
C_PEM = CENT_PEM;
dxdt_CENT_PEM = -(CL_PEM/V1_PEM)*CENT_PEM - (Q_PEM/V1_PEM)*CENT_PEM
                + (Q_PEM/V2_PEM)*PERI_PEM;
dxdt_PERI_PEM = (Q_PEM/V1_PEM)*CENT_PEM - (Q_PEM/V2_PEM)*PERI_PEM;

// ── Tisotumab vedotin ADC: archetype-2 conjugate PK + bespoke payload ─
// [refactor] bespoke, non-archetypal 3rd compartment (MMAE_TVADC): fed by a
// named fraction (FR_MMAE_TVADC, was a hardcoded 0.4) of the conjugates own
// clearance flux, first-order payload decay. Not TMDD (no receptor pool);
// does not match archetypes 1-4 -- kept as-is, see cc_refactor_notes.md.
C_TVADC = CENT_TVADC;
dxdt_CENT_TVADC = -(CL_TVADC/V1_TVADC)*CENT_TVADC - (Q_TVADC/V1_TVADC)*CENT_TVADC
                  + (Q_TVADC/V2_TVADC)*PERI_TVADC;
dxdt_PERI_TVADC = (Q_TVADC/V1_TVADC)*CENT_TVADC - (Q_TVADC/V2_TVADC)*PERI_TVADC;
dxdt_MMAE_TVADC = FR_MMAE_TVADC * (CL_TVADC/V1_TVADC) * CENT_TVADC
                  - k_dec_MMAE_TVADC * MMAE_TVADC;
// [build-fix] "if(MMAE_TVADC < 0) MMAE_TVADC = 0;" removed -- see the
// build-defect note above dxdt_VEGF below and cc_refactor_notes.md.

// ── VEGF dynamics (unchanged, disease-side; reads C_BEV) ─────────────
double VEGF_prod = ksyn_VEGF * (TV / TV0);
dxdt_VEGF = VEGF_prod - kdeg_VEGF * VEGF - VEGF_bind_BEV;
// [build-fix] Pre-existing upstream defect, independent of the $CMT+$INIT
// one above: mrgsolve 2.0.1 passes every $CMT state into $ODE as a
// `const double&`, so a direct "if(X < val) X = val;" clamp on a bare
// compartment name is a hard compile error ("assignment of read-only
// reference"). It was also always a behavioral no-op even in mrgsolve
// versions where it compiled: only dxdt_* feeds the integrator, so
// reassigning the state variable itself inside $ODE can never affect the
// next solver step or the reported trajectory (same reasoning as
// UPSTREAM_ISSUES.md #36, FLUID_EX/PR_FRAC in age-related-macular-
// degeneration). All nine such clamps in this file (MMAE_TVADC, VEGF,
// ADDUCT_CIS, RT_SF, CD8T, PDL1_exp, TV, SCCAg, HPVload) were deleted
// outright here -- confirmed dead, not a numeric change. Logged as a new
// UPSTREAM_ISSUES.md entry; the checked-in original still has all nine
// and does not build.

// ── Platinum-DNA adducts (was Pt_DNA; cisplatins own PD state) ───────
// [refactor] EFFECT_CIS = ADDUCT_CIS: a linear (non-Hill) PD turnover state
// already in the original -- no saturating nonlinearity to approximate, so
// this is a direct alias/rename, not a fit. See cc_refactor_notes.md.
dxdt_ADDUCT_CIS = k_adduct * C_CIS - k_repair * ADDUCT_CIS;
EFFECT_CIS = ADDUCT_CIS;

// ── Radiation cumulative damage (LQ model, continuous approximation) ─
// Effective alpha increases with concurrent cisplatin sensitization
double alpha_eff = alpha_rt * (1.0 + (k_radiosens - 1.0) * (ADDUCT_CIS / (ADDUCT_CIS + 0.3)));
double dose_rate = RT_flag * 2.0;   // 2 Gy/fraction-day equivalent, active only when RT_flag=1
double rt_damage_rate = alpha_eff * dose_rate + beta_rt * dose_rate * dose_rate;
dxdt_RT_SF = rt_damage_rate - k_reoxy * RT_SF;

// ── CD8+ T cell dynamics (unchanged, disease-side; reads EFFECT_PEM) ──
double ICI_effect = 1.0 + EFFECT_PEM * (k_ICI - 1.0);
double exhaustion  = k_exhaust * (TV / TV_max) * (PDL1_exp / PDL1_0);
dxdt_CD8T = k_CD8_in * ICI_effect - k_CD8_out * CD8T - exhaustion * CD8T;

// ── PD-L1 adaptive expression (IFN-gamma-like feedback from CD8T) ─────
dxdt_PDL1_exp = k_PDL1_up * CD8T * (1 - PDL1_exp/PDL1_max);

// ── Tumor volume (Gompertz + multi-modal kill) ────────────────────────
// [refactor] EFFECT_PAC/EFFECT_BEV/EFFECT_TVADC: each already a plain Hill
// ratio in the original (pac_eff, kill_bevs inline ratio, adc_eff) -- pulled
// out as named Emax/EC50/Hill terms with EMAX=1/GAMMA=1 where the original
// had no explicit Emax/Hill (rename, not a fit; see cc_refactor_notes.md).
double grow_term = kg * TV * log(TV_max / TV);
double repop_term = RT_flag * k_repop * TV;              // accelerated repopulation during prolonged RT
double kill_RT   = k_kill_RT * (1 - exp(-RT_SF)) * TV;
double kill_Pt   = k_kill_Pt * EFFECT_CIS * TV;
EFFECT_PAC = EMAX_PAC * pow(C_PAC, GAMMA_PAC) / (pow(EC50_PAC, GAMMA_PAC) + pow(C_PAC, GAMMA_PAC));
double kill_pac  = k_kill_pac * EFFECT_PAC * TV;
double kill_ICI  = k_kill_ICI * EFFECT_PEM * CD8T * TV;
EFFECT_TVADC = EMAX_TVADC * pow(MMAE_TVADC, GAMMA_TVADC) / (pow(EC50_TVADC, GAMMA_TVADC) + pow(MMAE_TVADC, GAMMA_TVADC));
double kill_ADC  = k_kill_ADC * EFFECT_TVADC * TV;
EFFECT_BEV = EMAX_BEV * pow(C_BEV, GAMMA_BEV) / (pow(EC50_BEV, GAMMA_BEV) + pow(C_BEV, GAMMA_BEV));
double kill_bev  = k_kill_bev * EFFECT_BEV * TV;
dxdt_TV = grow_term + repop_term - kill_RT - kill_Pt - kill_pac - kill_ICI - kill_ADC - kill_bev;

// ── SCC-Ag (turnover, proportional to tumor) ──────────────────────────
double SCC_prod = ksyn_SCC * TV;
dxdt_SCCAg = SCC_prod - kdeg_SCC * SCCAg;

// ── HPV viral load (production from tumor, clearance via immune/tx) ──
double HPV_clear_eff = k_HPV_clear * (1.0 + EFFECT_PEM * (k_HPV_ICI - 1.0));
dxdt_HPVload = k_HPV_prod * (TV/TV0) - HPV_clear_eff * HPVload;
// [build-fix] the "if(TV<0.01) TV=0.01;", "if(SCCAg<0.1) SCCAg=0.1;", and
// "if(HPVload<0) HPVload=0;" clamps (plus the six more above -- MMAE_TVADC,
// VEGF, ADDUCT_CIS, RT_SF, CD8T, PDL1_exp) were all removed as part of the
// same pre-existing, always-dead-clamp defect; see the note above
// dxdt_VEGF and cc_refactor_notes.md.

$TABLE
// [refactor/build-fix] C_<STEM>/EFFECT_<STEM> are recomputed here, directly
// from the live $CMT state, rather than trusting the copies $ODE last set.
// Verification found that at a timestamp where a dose event and a
// requested output row coincide (e.g. t=0 with a bolus dose there), the
// qspserver mrgsolve_api run reports that row using the state *after* the
// event but *before* $ODE is re-invoked -- so a $GLOBAL double only ever
// written inside $ODE (as C_CIS/EFFECT_CIS/... were) can read one step
// stale at exactly that row (confirmed: CIS_Conc read 0 instead of the
// correct 4.5333 at the post-dose t=0 row using the $ODE-only version).
// The original never has this problem because "capture CIS_Conc = CIS_C1;"
// reads the bare compartment directly, every time, in $TABLE. Recomputing
// each C_<STEM>/EFFECT_<STEM> again here (redundant with the identical
// $ODE computation used for cross-compound coupling during integration,
// e.g. VEGF_bind_BEV, kill_Pt) closes that gap and reproduces the
// original’s own behavior exactly -- confirmed by the verification run
// in cc_refactor_notes.md (max abs diff 0.0 at every dose timestamp,
// including simultaneous dose/observation rows).
C_CIS = CENT_CIS;
C_PAC = CENT_PAC;
C_BEV = CENT_BEV;
C_PEM = CENT_PEM;
C_TVADC = CENT_TVADC;
EFFECT_CIS = ADDUCT_CIS;
EFFECT_PAC = EMAX_PAC * pow(C_PAC, GAMMA_PAC) / (pow(EC50_PAC, GAMMA_PAC) + pow(C_PAC, GAMMA_PAC));
EFFECT_BEV = EMAX_BEV * pow(C_BEV, GAMMA_BEV) / (pow(EC50_BEV, GAMMA_BEV) + pow(C_BEV, GAMMA_BEV));
EFFECT_TVADC = EMAX_TVADC * pow(MMAE_TVADC, GAMMA_TVADC) / (pow(EC50_TVADC, GAMMA_TVADC) + pow(MMAE_TVADC, GAMMA_TVADC));
EFFECT_PEM = ICI_flag * CPS_high;

// [refactor] legacy output names kept identical to the original (for
// direct, same-name verification) but now computed from the renamed state.
double CIS_Conc    = C_CIS;
double PAC_Conc    = C_PAC;
double BEV_Conc    = C_BEV;
double PEMBRO_Conc = C_PEM;
double TVADC_Conc  = C_TVADC;
double MMAE_lvl    = MMAE_TVADC;
double VEGF_free   = VEGF;
double PtDNA_rel   = ADDUCT_CIS;
double RT_damage   = RT_SF;
double TumorVol    = TV;
double SCCAg_lvl   = SCCAg;
double HPV_rel     = HPVload;
double CD8T_rel    = CD8T;
double PDL1_rel    = PDL1_exp;
double TV_change   = (TV - TV0) / TV0 * 100;

$CAPTURE
CIS_Conc PAC_Conc BEV_Conc PEMBRO_Conc TVADC_Conc MMAE_lvl VEGF_free PtDNA_rel RT_damage
TumorVol SCCAg_lvl HPV_rel CD8T_rel PDL1_rel TV_change
C_CIS C_PAC C_BEV C_PEM C_TVADC
EFFECT_CIS EFFECT_PAC EFFECT_BEV EFFECT_PEM EFFECT_TVADC
'

## ---------------------------------------------------------------
## Compile the model
## ---------------------------------------------------------------
mod_cc <- mcode("cc_model_refactored", code_cc)

## ---------------------------------------------------------------
## Dosing event functions
## ---------------------------------------------------------------

## Cisplatin: 40 mg/m² weekly x5-6 during CCRT (BSA ~1.7 m² -> 68 mg)
dose_cisplatin <- function(n_doses = 6, interval_d = 7, V1 = 15) {
  dose_mg <- 68
  ev(cmt="CENT_CIS", amt=dose_mg/V1, time=seq(0, (n_doses-1)*interval_d, by=interval_d))
}

## EBRT + brachytherapy: continuous RT_flag "on" via parameter update is
## simulated as a fixed active window (day 0-49, ~5 weeks EBRT + boost)
## handled via idata/param switching per-scenario below (RT_flag=1 during window)

## Paclitaxel: 175 mg/m² q3w interval chemo (recurrent setting, ~300 mg)
dose_paclitaxel <- function(n_cycles = 6, interval_d = 21, start_d = 0, V1 = 6.5) {
  dose_mg <- 300
  ev(cmt="CENT_PAC", amt=dose_mg*1000/V1, time=seq(start_d, start_d+(n_cycles-1)*interval_d, by=interval_d))
}

## Bevacizumab: 15 mg/kg q3w IV (GOG-240) = ~1050 mg per dose
dose_bevacizumab <- function(start_d = 0, n_doses = 20, interval_d = 21, V1 = 2.91) {
  dose_mg <- 1050
  times <- seq(start_d, start_d + (n_doses-1)*interval_d, by=interval_d)
  ev(cmt="CENT_BEV", amt=dose_mg/V1, time=times)
}

## Pembrolizumab: 200 mg q3w IV (KEYNOTE-A18 / KEYNOTE-826)
dose_pembrolizumab <- function(start_d = 0, n_doses = 35, interval_d = 21, V1 = 3.34) {
  dose_mg <- 200
  times <- seq(start_d, start_d + (n_doses-1)*interval_d, by=interval_d)
  ev(cmt="CENT_PEM", amt=dose_mg/V1, time=times)
}

## Tisotumab vedotin: 2.0 mg/kg q3w IV (innovaTV 301; ~140mg for 70kg)
dose_tisotumab <- function(start_d = 0, n_doses = 20, interval_d = 21, V1 = 3.35) {
  dose_mg <- 140
  times <- seq(start_d, start_d + (n_doses-1)*interval_d, by=interval_d)
  ev(cmt="CENT_TVADC", amt=dose_mg/V1, time=times)
}

## ---------------------------------------------------------------
## Treatment Scenarios
## ---------------------------------------------------------------
sim_time <- seq(0, 730, by=1)  # 2-year simulation (days)

## S1: No treatment (natural history / untreated progression)
mod_S1 <- mod_cc %>% param(RT_flag=0, ICI_flag=0, CPS_high=1)
out_S1  <- mrgsim(mod_S1, end=730, delta=1)

## S2: Cisplatin-based concurrent chemoradiation (CCRT), RTOG-90-01 style
##     EBRT+brachy over ~7 weeks (day 0-49), weekly cisplatin x6
mod_S2 <- mod_cc %>% param(RT_flag=1, ICI_flag=0, CPS_high=1)
ev_S2 <- dose_cisplatin(6)
out_S2 <- mrgsim(mod_S2, events=ev_S2, end=730, delta=1)

## S3: CCRT + Pembrolizumab (KEYNOTE-A18 regimen: concurrent + maintenance
##     up to ~2 years); Pembro starts with CCRT, continues as maintenance
mod_S3 <- mod_cc %>% param(RT_flag=1, ICI_flag=1, CPS_high=1)
ev_S3 <- ev_seq(dose_cisplatin(6), dose_pembrolizumab(0, n_doses=35))
out_S3 <- mrgsim(mod_S3, events=ev_S3, end=730, delta=1)

## S4: Recurrent/metastatic — Chemo + Bevacizumab (GOG-240 regimen:
##     paclitaxel + cisplatin + bevacizumab)
mod_S4 <- mod_cc %>% param(RT_flag=0, ICI_flag=0, CPS_high=1)
ev_S4 <- ev_seq(dose_cisplatin(6, interval_d=21), dose_paclitaxel(6), dose_bevacizumab(0, n_doses=20))
out_S4 <- mrgsim(mod_S4, events=ev_S4, end=730, delta=1)

## S5: Recurrent/metastatic — Tisotumab vedotin monotherapy after
##     progression on platinum (innovaTV 301: post-platinum 2L+)
mod_S5 <- mod_cc %>% param(RT_flag=0, ICI_flag=0, CPS_high=1)
ev_S5 <- dose_tisotumab(0, n_doses=20)
out_S5 <- mrgsim(mod_S5, events=ev_S5, end=730, delta=1)

## S6: Recurrent/metastatic 1st line — Chemo + Bevacizumab + Pembrolizumab
##     (KEYNOTE-826 triplet, PD-L1 CPS>=1 population)
mod_S6 <- mod_cc %>% param(RT_flag=0, ICI_flag=1, CPS_high=1)
ev_S6 <- ev_seq(
  dose_cisplatin(6, interval_d=21),
  dose_paclitaxel(6),
  dose_bevacizumab(0, n_doses=20),
  dose_pembrolizumab(0, n_doses=35)
)
out_S6 <- mrgsim(mod_S6, events=ev_S6, end=730, delta=1)

## ---------------------------------------------------------------
## Summary: 24-month PFS proxy and key endpoints
## ---------------------------------------------------------------
summarize_scenario <- function(out, label) {
  df <- as.data.frame(out)
  pfs_d <- df %>% filter(TumorVol > 80) %>% pull(time) %>% min()
  pfs_d <- if(is.infinite(pfs_d)) ">730" else round(pfs_d)
  sccag_nadir <- min(df$SCCAg_lvl)
  sccag_nadir_t <- df$time[which.min(df$SCCAg_lvl)]
  tv_min <- min(df$TumorVol)
  best_resp <- round((tv_min - 40) / 40 * 100, 1)
  hpv_final <- round(tail(df$HPV_rel, 1), 3)
  data.frame(
    Scenario      = label,
    PFS_days      = pfs_d,
    SCCAg_nadir   = round(sccag_nadir, 2),
    SCCAg_nadir_t = round(sccag_nadir_t),
    BestResp_pct  = best_resp,
    HPV_final     = hpv_final
  )
}

summary_table <- rbind(
  summarize_scenario(out_S1, "S1: Untreated (natural history)"),
  summarize_scenario(out_S2, "S2: Cisplatin CCRT (RTOG-90-01 style)"),
  summarize_scenario(out_S3, "S3: CCRT + Pembrolizumab (KEYNOTE-A18)"),
  summarize_scenario(out_S4, "S4: Chemo+Bevacizumab, R/M (GOG-240)"),
  summarize_scenario(out_S5, "S5: Tisotumab vedotin, R/M (innovaTV 301)"),
  summarize_scenario(out_S6, "S6: Chemo+Bev+Pembro, R/M (KEYNOTE-826)")
)
print(summary_table)

## ---------------------------------------------------------------
## Visualization
## ---------------------------------------------------------------

combine_sims <- function(..., labels = NULL) {
  sims <- list(...)
  lapply(seq_along(sims), function(i) {
    df <- as.data.frame(sims[[i]])
    df$Scenario <- if(!is.null(labels)) labels[i] else paste0("S", i)
    df
  }) %>% bind_rows()
}

scenario_labels <- c(
  "S1: Untreated",
  "S2: Cisplatin CCRT",
  "S3: CCRT+Pembrolizumab",
  "S4: Chemo+Bev (GOG-240)",
  "S5: Tisotumab vedotin",
  "S6: Chemo+Bev+Pembro (KEYNOTE-826)"
)

all_sims <- combine_sims(out_S1, out_S2, out_S3, out_S4, out_S5, out_S6,
                         labels = scenario_labels)

## --- Plot 1: Tumor Volume over time ---
p1 <- ggplot(all_sims, aes(x=time, y=TumorVol, color=Scenario)) +
  geom_line(size=0.9) +
  geom_hline(yintercept=80, linetype="dashed", color="gray50") +
  annotate("text", x=680, y=90, label="PD threshold (2×BL)", size=2.8, color="gray50") +
  labs(title="Tumor Volume (cm³) — 6 Treatment Scenarios",
       x="Day", y="Tumor Volume (cm³)") +
  scale_color_brewer(palette="Set1") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7))

## --- Plot 2: SCC-Ag serum ---
p2 <- ggplot(all_sims, aes(x=time, y=SCCAg_lvl, color=Scenario)) +
  geom_line(size=0.9) +
  geom_hline(yintercept=1.5, linetype="dashed", color="darkgreen") +
  annotate("text", x=680, y=2, label="ULN 1.5 ng/mL", size=2.8, color="darkgreen") +
  labs(title="SCC-Ag Serum Level (ng/mL)",
       x="Day", y="SCC-Ag (ng/mL)") +
  scale_y_log10() +
  scale_color_brewer(palette="Set1") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

## --- Plot 3: Radiation damage + HPV viral load (S2/S3) ---
p3_rt <- ggplot(as.data.frame(out_S3), aes(x=time, y=RT_damage)) +
  geom_line(color="#607D8B") +
  labs(title="Cumulative Radiation Damage (S3, CCRT window)",
       x="Day", y="RT damage (relative)") +
  coord_cartesian(xlim=c(0,120)) +
  theme_bw(base_size=10)

p3_hpv <- ggplot(as.data.frame(out_S3), aes(x=time, y=HPV_rel)) +
  geom_line(color="#C0392B") +
  labs(title="HPV Viral Load Decline (S3, CCRT+Pembrolizumab)",
       x="Day", y="HPV load (relative log10)") +
  theme_bw(base_size=10)

## --- Plot 4: Drug PK (Cisplatin) ---
p4 <- ggplot(as.data.frame(out_S2), aes(x=time)) +
  geom_line(aes(y=CIS_Conc), color="#FF8F00") +
  labs(title="Cisplatin Central PK (S2, weekly CCRT)",
       x="Day", y="Cisplatin (µg/mL)") +
  coord_cartesian(xlim=c(0,60)) +
  theme_bw(base_size=10)

## --- Plot 5: VEGF suppression (bevacizumab scenarios) ---
p5 <- ggplot(all_sims %>% filter(Scenario %in% c(
    "S2: Cisplatin CCRT",
    "S4: Chemo+Bev (GOG-240)",
    "S6: Chemo+Bev+Pembro (KEYNOTE-826)"
  )), aes(x=time, y=VEGF_free, color=Scenario)) +
  geom_line(size=0.9) +
  labs(title="Free VEGF-A (Bevacizumab Scenarios)",
       x="Day", y="Free VEGF-A (ng/mL)") +
  scale_color_brewer(palette="Set2") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7))

## --- Combined figure ---
main_fig <- (p1 | p2) / (p4 | p3_rt) / (p5 | p3_hpv)
print(main_fig + plot_annotation(
  title    = "Cervical Cancer QSP Model — Simulation Results (PK/PD refactored)",
  subtitle = "HPV-driven Cervical SCC · 6 Treatment Scenarios · 2-Year Projection",
  caption  = "Calibrated to RTOG-90-01, GOG-240, KEYNOTE-A18/826, innovaTV 204/301"
))

## ---------------------------------------------------------------
## Key Parameter Calibration Notes
## ---------------------------------------------------------------
## Cisplatin (weekly CCRT):
##   - Standard 40 mg/m² weekly x5-6 concurrent with pelvic RT
##   - Radiosensitizer: Green JA 2001 Lancet meta-analysis (concurrent
##     platinum-based CT + RT improves OS, absolute benefit ~6% at 5yr)
##   - Rose PG 1999 NEJM (GOG-120): cisplatin-containing regimens superior
##     to hydroxyurea in locally advanced disease
##
## Radiotherapy (EBRT + brachytherapy):
##   - LQ model alpha/beta = 10 Gy for cervical SCC (Fowler 1989 Br J Radiol)
##   - EBRT 45-50 Gy/25 fx pelvis (± extended field to para-aortic nodes)
##   - Brachytherapy boost to total EQD2 >= 85 Gy to HR-CTV
##   - EMBRACE-I/II (Pötter R 2018/2021 Lancet Oncol/Radiother Oncol):
##     image-guided adaptive brachytherapy improves local control
##   - RTOG 90-01 (Eifel PJ 2004 JCO; Morris M 1999 NEJM): concurrent
##     cisplatin+5FU+RT superior to extended-field RT alone
##
## Bevacizumab (recurrent/metastatic):
##   - GOG-240 (Tewari KS et al. 2014 NEJM, PMID 24552320): chemo+bev
##     improved OS 16.8 vs 13.3 mo (HR 0.71) in recurrent/persistent/
##     metastatic cervical cancer
##   - t½≈20 days (IgG1); CL=0.207 L/day (Lu 2008 Cancer Chemother Pharmacol)
##   - 15 mg/kg q3w -> Cmax≈360 µg/mL
##
## Pembrolizumab:
##   - KEYNOTE-A18 (Lorusso D et al. 2024 Lancet; Monk BJ et al. 2023 JCO):
##     pembrolizumab + CCRT improved PFS (HR 0.70) vs CCRT alone in
##     high-risk locally advanced cervical cancer (FIGO 2014 IB2-IVA)
##   - KEYNOTE-826 (Colombo N et al. 2021 NEJM, PMID 34534429): pembro +
##     chemo ± bevacizumab improved OS in 1st-line recurrent/metastatic
##     disease (PD-L1 CPS>=1 population showed greatest benefit)
##   - PK: Ahamadi M et al. 2017 CPT Pharmacometrics Syst Pharmacol;
##     t½≈22 days, linear PK at approved doses, near-saturating receptor
##     occupancy at 200mg q3w
##
## Tisotumab vedotin (ADC):
##   - innovaTV 204 (Coleman RL et al. 2021 Lancet Oncol, PMID 34310922):
##     single-agent activity in previously treated recurrent/metastatic
##     disease (ORR ~24%)
##   - innovaTV 301 (Vergote I et al. 2024 J Clin Oncol/NEJM Evid):
##     confirmed OS benefit vs investigator-choice chemo in 2L+ setting
##     (median OS 11.5 vs 9.5 mo, HR 0.70)
##   - Mechanism: anti-tissue factor (TF) mAb conjugated to MMAE via
##     protease-cleavable linker; bystander killing via membrane-permeable
##     MMAE independent of TF expression on neighboring cells
##   - Key AEs: ocular toxicity (conjunctivitis/keratitis), bleeding events
##
## Tumor growth / biomarker calibration:
##   - Gompertz model, doubling time consistent with locally advanced
##     bulky cervical SCC pre-treatment growth kinetics
##   - SCC-Ag half-life ~2.8 days (Gaarenstroom KN 2000 Int J Gynecol
##     Cancer); tracks tumor burden and recurrence
##   - HPV viral load decline post-CCRT reflects both direct cytoreduction
##     and immune-mediated clearance (enhanced under ICI exposure)
## =============================================================================
