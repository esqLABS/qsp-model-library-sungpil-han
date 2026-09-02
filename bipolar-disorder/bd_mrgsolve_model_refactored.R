## ============================================================================
## Bipolar Disorder QSP Model – mrgsolve ODE Implementation (PK/PD REFACTOR)
## ============================================================================
## Disease    : Bipolar Disorder (BD-I / BD-II)
## Model Type : Multi-compartment PK/PD with neurotransmitter dynamics,
##              signal transduction, neuroplasticity, and circadian coupling
## ODE Compts : 22 (see $CMT block)
## Scenarios  : 6 (acute mania / mixed / BD depression / maintenance ×2 / clozapine)
## References : See bd_references.md
## Calibration: Parameters anchored to BALANCE, EMBOLDEN, STRIDE-BD, CANMAT
##              and individual PK studies (see inline CALIB: notes)
## ----------------------------------------------------------------------------
## [PK/PD refactor] This is the pluggable-PK sibling of bd_mrgsolve_model.R,
## produced per FORK_WORKFLOW_GUIDE.md Part 2. The original is never edited.
## Four compounds refactored (the file's entire exogenous-compound set):
## Lithium (LI), Valproate (VPA), Quetiapine (QTP, + its active metabolite
## norquetiapine, NQT), Lamotrigine (LTG). See bd_refactor_notes.md for the
## full archetype/verification writeup.
## Author: QSP Library (auto-generated via Claude Code Routine)
## Date  : 2026-06-25 (original) / refactor added 2026-09-02
## R pkg : mrgsolve ≥ 1.0.0, dplyr, ggplot2
## ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

## ---------------------------------------------------------------------------
## MODEL CODE (inline mrgsolve specification) — refactored, pluggable PK
## ---------------------------------------------------------------------------
code_refactored <- '
$PROB
Bipolar Disorder QSP Model (PK/PD refactor)
Lithium / Valproate / Quetiapine / Lamotrigine PK-PD
[Note: the original prologue also names Aripiprazole, but no Aripiprazole
PK or PD is implemented anywhere in this file (checked from code, not
assumed) -- disclosed in bd_refactor_notes.md, not our compound to add.]

$PARAM @annotated
// ---- Lithium (LI) PK (archetype 3: depot + central + peripheral, linear) ----
// CALIB: Finley et al. J Clin Pharmacol 1995; Sproule 2002
// [PK/PD refactor] ka_Li->KA_LI, Vc_Li->V1_LI, Vp_Li->V2_LI, Q_Li->Q_LI,
// CL_Li->CL_LI, F_Li->F_LI (rename, same values). F_LI is declared but
// never multiplied into the original’s own dxdt_Li_central -- preserved
// unused exactly as the original had it (see bd_refactor_notes.md).
KA_LI   : 0.80  : Lithium absorption rate constant (h-1)
V1_LI   : 30    : Lithium central volume (L)  [0.4 L/kg x 75 kg]
V2_LI   : 15    : Lithium peripheral volume (L)
Q_LI    : 4.0   : Lithium inter-compartment clearance (L/h)
CL_LI   : 1.80  : Lithium renal clearance (L/h)  [~0.024 L/h/kg]
F_LI    : 1.0   : Lithium bioavailability [declared, never used in $ODE -- preserved as-is]

// ---- Valproate (VPA) PK (archetype 3 minus peripheral: depot + central, ----
// ---- nonlinear protein binding) ----
// CALIB: Perucca 2002; Johannessen 2000
// [PK/PD refactor] ka_VPA->KA_VPA, Vc_VPA->V1_VPA, CL_VPA->CL_VPA (unchanged
// name), fu_VPA0->FU_VPA0, Km_fu->KM_VPA (rename, same values).
KA_VPA  : 1.50  : VPA absorption rate constant (h-1)
V1_VPA  : 14    : VPA apparent central volume (L) at low concentration
CL_VPA  : 0.55  : VPA total clearance (L/h)  [~7.5 mL/min]
FU_VPA0 : 0.10  : VPA free fraction at low concentration
KM_VPA  : 50    : VPA concentration (ug/mL) for protein saturation (Km)

// ---- Quetiapine (QTP) PK (archetype 3 minus peripheral: depot + central) ----
// ---- plus norquetiapine (NQT), a bespoke 1-cmt metabolite chain sharing ----
// ---- the parent’s own volume divisor exactly as the original did ----
// CALIB: DeVane & Nemeroff 2001; Gefvert et al. 2001
// [PK/PD refactor] ka_QTP->KA_QTP, Vc_QTP->V1_QTP, CL_QTP->CL_QTP (unchanged
// name), F_QTP->F_QTP (unchanged name), km_QTP->FMET_QTP, CL_NQT->CL_NQT
// (unchanged name) -- rename, same values.
KA_QTP  : 1.10  : Quetiapine absorption rate constant (h-1)
V1_QTP  : 900   : Quetiapine apparent Vd (L)
CL_QTP  : 250   : Quetiapine clearance (L/h) [CYP3A4]
F_QTP   : 0.09  : Quetiapine bioavailability (9%)
FMET_QTP: 0.25  : Fraction of QTP metabolised to norquetiapine
CL_NQT  : 80    : Norquetiapine clearance (L/h)

// ---- Lamotrigine (LTG) PK (archetype 3 minus peripheral: depot + central) ----
// CALIB: Doose et al. 2003; Calabrese et al. 1999
// [PK/PD refactor] ka_LTG->KA_LTG, Vc_LTG->V1_LTG, CL_LTG->CL_LTG (unchanged
// name) -- rename, same values. No bioavailability term in the original,
// preserved as-is (F implicitly 1).
KA_LTG  : 0.45  : Lamotrigine absorption rate constant (h-1)
V1_LTG  : 105   : Lamotrigine apparent Vd (L)
CL_LTG  : 1.85  : Lamotrigine baseline clearance (L/h) - monotherapy
// Note: CL_LTG doubles with VPA co-administration inhibition is captured via
//       inducer flag (see $MAIN) [same disclaimer as the original -- no such
//       flag actually exists in this file’s $MAIN; preserved verbatim,
//       logged in bd_refactor_notes.md as a pre-existing documentation-only
//       inconsistency, not a code defect]

// ---- Dopamine PD parameters ----
// CALIB: Montague et al. 2004; Nestler & Carlezon 2006
DA_base : 1.0   : Baseline normalized DA neurotransmission index (dimensionless)
kDA_syn : 0.20  : DA synthesis rate constant (h-1)
kDA_deg : 0.20  : DA degradation rate constant (h-1)
// [PK/PD refactor] EC50_DA_D2/Imax_D2 -> EC50_QTP_D2/EMAX_QTP_D2 below
// (rename into the Hill interface; kept here removed, see QTP block)

// ---- Serotonin PD parameters ----
// CALIB: Sprouse & Aghajanian 1987; based on 5-HT1A autoreceptor kinetics
HT5_base : 1.0   : Baseline 5-HT neurotransmission index
k5HT_syn : 0.15  : 5-HT synthesis rate constant (h-1)
k5HT_deg : 0.15  : 5-HT degradation rate constant (h-1)

// ---- GSK-3beta pathway (lithium / VPA pharmacodynamics) ----
// CALIB: Jope & Johnson 2004; Ryves & Harwood 2001; Li IC50 ~2 mM
// [PK/PD refactor] the original shares ONE Emax_GSK=0.90 between lithium and
// VPA’s independent GSK-3b inhibition terms; duplicated (same value) into
// EMAX_LI_GSK3 and EMAX_VPA_GSK3 so each compound keeps its own explicit
// Hill parameter set, per the guide’s "separate EFFECT_<STEM>" rule. Not a
// refit -- both still equal 0.90.
GSK3_base  : 1.0   : Baseline GSK-3beta activity
kGSK_syn   : 0.05  : GSK-3beta synthesis/activation rate (h-1)
kGSK_deg   : 0.05  : GSK-3beta degradation rate (h-1)
EC50_LI_GSK3  : 0.70  : Lithium [mEq/L] for 50% GSK-3beta inhibition (was IC50_Li_GSK)
EMAX_LI_GSK3  : 0.90  : Max GSK-3beta inhibition by lithium (was Emax_GSK)
GAMMA_LI_GSK3 : 1.0   : Hill coefficient [new: original had no explicit Hill exponent]
EC50_VPA_GSK3 : 60    : VPA free [ug/mL] for 50% GSK-3beta inhibition (was IC50_VPA_GSK)
EMAX_VPA_GSK3 : 0.90  : Max GSK-3beta inhibition by VPA (was Emax_GSK, duplicated)
GAMMA_VPA_GSK3: 1.0   : Hill coefficient [new]

// ---- BDNF / Neuroplasticity ----
// CALIB: Castren & Rantamaki 2010; Duman & Monteggia 2006
BDNF_base  : 1.0   : Baseline BDNF index
kBDNF_syn  : 0.03  : BDNF synthesis rate (h-1)
kBDNF_deg  : 0.03  : BDNF degradation rate (h-1)
EMAX_LI_BDNF : 0.50 : Max BDNF increase by lithium (at 1 mEq/L) (was Emax_BDNF_Li)
EC50_LI_BDNF : 0.40 : Lithium [mEq/L] for 50% max BDNF effect (was EC50_BDNF_Li)
GAMMA_LI_BDNF: 1.0  : Hill coefficient [new]
kNprot       : 0.01 : Neuroprotection coupling constant (GSK->BDNF)

// ---- Neuroinflammation (IL-6 surrogate) ----
// CALIB: Goldsmith et al. 2016 meta-analysis of BD inflammatory markers
IL6_base : 1.0   : Baseline IL-6 index
kIL6_syn : 0.04  : IL-6 synthesis rate (h-1)
kIL6_deg : 0.04  : IL-6 degradation rate (h-1)
EMAX_LI_IL6 : 0.40 : Max IL-6 reduction by lithium (was Emax_IL6_Li)
EC50_LI_IL6 : 0.50 : Li [mEq/L] for half-max IL-6 reduction (was EC50_IL6_Li)
GAMMA_LI_IL6: 1.0  : Hill coefficient [new]

// ---- HPA Axis (cortisol index) ----
// CALIB: Cervantes et al. 2001; Daban et al. 2005
Cort_base : 1.0  : Baseline cortisol index
kCort_prod: 0.08 : Cortisol production rate (h-1) [circadian average]
kCort_deg : 0.08 : Cortisol degradation rate (h-1)

// ---- Circadian oscillator ----
// CALIB: Berson et al. 2002; Frank et al. 2005 SRM data
omega      : 0.262 : Circadian angular freq (rad/h) = 2pi/24
Amp_circ   : 0.30  : Circadian amplitude (0-1)

// ---- Mood State PD (YMRS / MADRS) ----
// CALIB: Keck et al. 2003; Calabrese et al. 2005; CANMAT 2018
YMRS_base : 25.0  : Baseline YMRS at episode onset (mania)
MADRS_base: 30.0  : Baseline MADRS at episode onset (depression)
kYMRS_nat : 0.01  : Natural YMRS improvement rate (h-1) without treatment
kMADRS_nat: 0.005 : Natural MADRS improvement rate (h-1)

// Effect sizes on YMRS (Li/VPA/QTP from controlled trials)
// CALIB: Bowden et al. 1994 (VPA); Lithium vs placebo meta Cipriani 2013
// [PK/PD refactor] Emax_YMRS_Li/EC50_YMRS_Li -> EMAX_LI_YMRS/EC50_LI_YMRS;
// Emax_YMRS_VPA/EC50_YMRS_VPA -> EMAX_VPA_YMRS/EC50_VPA_YMRS;
// Emax_YMRS_QTP/EC50_YMRS_QTP -> EMAX_QTP_YMRS/EC50_QTP_YMRS (rename, same
// values). GAMMA_*=1 new explicit params throughout.
EMAX_LI_YMRS  : 18  : Li max YMRS reduction (points, anchored ~3-week response)
EC50_LI_YMRS  : 0.65 : Li [mEq/L] 50% max YMRS reduction
GAMMA_LI_YMRS : 1.0 : Hill coefficient [new]
EMAX_VPA_YMRS : 20  : VPA max YMRS reduction
EC50_VPA_YMRS : 65  : VPA free [ug/mL] 50% max YMRS reduction
GAMMA_VPA_YMRS: 1.0 : Hill coefficient [new]
EMAX_QTP_YMRS : 16  : QTP max YMRS reduction (BOLDER)
EC50_QTP_YMRS : 120 : QTP [ng/mL] 50% max YMRS reduction
GAMMA_QTP_YMRS: 1.0 : Hill coefficient [new]

// QTP D2-receptor occupancy driving DA suppression -- the original reuses
// EC50_YMRS_QTP’s own value (120) as the D2 binding constant too; kept as
// a disclosed duplicate (EC50_QTP_D2), not a refit, per the guide’s allowance
// for a shared original EC50 reused at two genuinely different endpoints.
EMAX_QTP_D2  : 0.8  : Max D2R-mediated DA suppression (was Imax_D2)
EC50_QTP_D2  : 120  : QTP [ng/mL] for 50% D2 occupancy (was reusing EC50_YMRS_QTP)
GAMMA_QTP_D2 : 1.0  : Hill coefficient [new]

// Effect sizes on MADRS (QTP/LTG from EMBOLDEN)
// CALIB: Young et al. 2010 EMBOLDEN I; Calabrese 2003 LTG depression
// [PK/PD refactor] Emax_MADRS_QTP/EC50_MADRS_QTP -> EMAX_QTP_MADRS/EC50_QTP_MADRS;
// Emax_MADRS_LTG/EC50_MADRS_LTG -> EMAX_LTG_MADRS/EC50_LTG_MADRS (rename).
EMAX_QTP_MADRS : 15 : QTP max MADRS reduction
EC50_QTP_MADRS : 80 : QTP [ng/mL] 50% max MADRS reduction
GAMMA_QTP_MADRS: 1.0 : Hill coefficient [new]
EMAX_LTG_MADRS : 10 : LTG max MADRS reduction
EC50_LTG_MADRS : 2.5 : LTG [ug/mL] 50% max MADRS reduction
GAMMA_LTG_MADRS: 1.0 : Hill coefficient [new]

// Norquetiapine (NQT, QTP’s active metabolite) drives a serotonergic
// synthesis-stimulation term the original wrote inline inside $ODE
// (1.0 + 0.3*NQT_conc/(80.0+NQT_conc)); pulled out as an explicit Hill
// interface EFFECT_QTP_5HT (rename, not a refit -- same 0.3/80 values).
EMAX_QTP_5HT  : 0.3 : Max serotonin-synthesis stimulation by norquetiapine
EC50_QTP_5HT  : 80  : NQT [ng/mL] for 50% max 5-HT stimulation
GAMMA_QTP_5HT : 1.0 : Hill coefficient [new]

// Combination effect parameter
ALPHA_COMBO_LI_QTP : 0.4 : Interaction factor Li+QTP on MADRS (CANMAT 2018) (was ALPHA_combo)

// ---- Body Weight (QTP / sedation metabolic effects) ----
// [PK/PD refactor] the original’s Wt_gain_rate = kWt_gain*NQT_conc is
// genuinely LINEAR (no saturation term) -- kept linear, not forced into an
// Emax/EC50 Hill shape (see bd_refactor_notes.md).
Wt_base   : 0.0  : Baseline weight change (kg)
KWT_QTP   : 0.003 : QTP(norquetiapine)-driven weight gain rate (kg/h per ng/mL) (was kWt_gain)

$CMT @annotated
// PK compartments
GUT_LI       : Lithium gut (mmol)
CENT_LI      : Lithium central (mmol)
PERI_LI      : Lithium peripheral (mmol)
GUT_VPA      : Valproate gut (mg)
CENT_VPA     : Valproate central (mg)
GUT_QTP      : Quetiapine gut (mg)
CENT_QTP     : Quetiapine central (mg)
CENT_NQT     : Norquetiapine central (mg)
GUT_LTG      : Lamotrigine gut (mg)
CENT_LTG     : Lamotrigine central (mg)
// PD compartments
DA_index     : Dopamine neurotransmission index
HT5_index    : Serotonin neurotransmission index
GSK3_activ   : GSK-3beta activity index
BDNF_level   : BDNF concentration index
IL6_level    : IL-6 neuroinflammation index
Cortisol_idx : Cortisol / HPA axis index
// Mood state
YMRS_score   : YMRS score (mania)
MADRS_score  : MADRS score (bipolar depression)
// Metabolic
Weight_chg   : Body weight change (kg)
// Circadian
Circ_state   : Circadian phase oscillator (dimensionless)
Circ_deriv   : Circadian derivative state
// Functioning
GAF_score    : Global Assessment of Functioning

$GLOBAL
// ---- Canonical, single-source-of-truth concentration formulas ----
// [PK/PD refactor -- "normalize duplicate concentration sites, then
// redirect"] The original independently re-derived every one of these five
// concentrations twice under two different names: once in $MAIN (feeding
// the PD/Hill effect terms) and again in $TABLE (feeding $CAPTURE, e.g.
// Li_conc vs Lithium_mEqL, QTP_conc vs QTP_ngmL). A minimal reproduction
// against this exact file’s own timing (see bd_refactor_notes.md) confirmed
// $MAIN and $TABLE are NOT interchangeable here: $MAIN runs once per
// interval using state as of the start of that interval, while $TABLE runs
// after the ODE has been advanced to this row’s own time -- so the two
// blocks’ independent recomputations genuinely differ by one row’s worth of
// decay/absorption, exactly as the untouched original already does. Moving
// either block’s calculation into the other (or sharing one mutable value
// across both) would silently change the simulated trajectory, which the
// guide explicitly forbids. The fix that removes the duplication WITHOUT
// changing timing: each concentration’s formula is authored exactly once,
// here, as a macro; $MAIN and $TABLE each still declare their own `double`
// (so each keeps evaluating against its own block’s live state, exactly as
// before) but by expanding the same macro text instead of retyping the
// arithmetic by hand. C_<STEM>/EFFECT_<STEM> are $MAIN-local + $TABLE-local
// + $CAPTURE only, never $PARAM (mrgsolve 2.0.1 compiles $PARAM members as
// read-only references; a value recomputed from live state cannot also be
// a $PARAM -- same constraint already documented in numerous sibling
// refactors, e.g. copd/copd_refactor_notes.md, atrial-fibrillation/
// af_refactor_notes.md).
#define LI_CONC_EXPR       (CENT_LI  / V1_LI)
#define VPA_TOTAL_EXPR     (CENT_VPA / V1_VPA)
#define VPA_FREE_EXPR(ctot) ((ctot) * FU_VPA0 * (1.0 + (ctot) / KM_VPA))
#define QTP_CONC_EXPR      (CENT_QTP / (V1_QTP / 1000.0))
#define NQT_CONC_EXPR      (CENT_NQT / (V1_QTP / 1000.0))
#define LTG_CONC_EXPR      (CENT_LTG / V1_LTG)

$MAIN
// ---- Derived PK concentrations (single-source formulas via $GLOBAL macros,
// evaluated here exactly where the original evaluated Li_conc/VPA_conc/
// QTP_conc/NQT_conc/LTG_conc -- i.e. once per interval, feeding the PD
// effect terms below) ----
double C_LI = LI_CONC_EXPR;                        // mEq/L -- Lithium, the concentration every Li PD term reads
double VPA_total_conc = VPA_TOTAL_EXPR;            // ug/mL total (protein-bound + free)
double C_VPA = VPA_FREE_EXPR(VPA_total_conc);      // ug/mL free -- Valproate, the concentration every VPA PD term reads (non-linear protein binding, kept exactly as the original -- not flattened)
double C_QTP = QTP_CONC_EXPR;                      // ng/mL -- Quetiapine (parent), the concentration every parent-QTP PD term reads
double C_NQT = NQT_CONC_EXPR;                      // ng/mL -- Norquetiapine (active metabolite), drives the 5-HT and weight-gain terms
double C_LTG = LTG_CONC_EXPR;                       // ug/mL -- Lamotrigine, the concentration the MADRS term reads

// ---- Hill interface: named EFFECT_<STEM>[_<TARGET>] terms ----
// Every term below is a rename of an already-plain Emax*C/(EC50+C) ratio in
// the original (implicit Hill coefficient of 1, so GAMMA_*=1 throughout) --
// a rename, not a refit. See bd_refactor_notes.md.
double EFFECT_LI_GSK3  = EMAX_LI_GSK3  * pow(C_LI,  GAMMA_LI_GSK3)  / (pow(EC50_LI_GSK3,  GAMMA_LI_GSK3)  + pow(C_LI,  GAMMA_LI_GSK3));
double EFFECT_LI_BDNF  = EMAX_LI_BDNF  * pow(C_LI,  GAMMA_LI_BDNF)  / (pow(EC50_LI_BDNF,  GAMMA_LI_BDNF)  + pow(C_LI,  GAMMA_LI_BDNF));
double EFFECT_LI_IL6   = EMAX_LI_IL6   * pow(C_LI,  GAMMA_LI_IL6)   / (pow(EC50_LI_IL6,   GAMMA_LI_IL6)   + pow(C_LI,  GAMMA_LI_IL6));
double EFFECT_LI_YMRS  = EMAX_LI_YMRS  * pow(C_LI,  GAMMA_LI_YMRS)  / (pow(EC50_LI_YMRS,  GAMMA_LI_YMRS)  + pow(C_LI,  GAMMA_LI_YMRS));

double EFFECT_VPA_GSK3 = EMAX_VPA_GSK3 * pow(C_VPA, GAMMA_VPA_GSK3) / (pow(EC50_VPA_GSK3, GAMMA_VPA_GSK3) + pow(C_VPA, GAMMA_VPA_GSK3));
double EFFECT_VPA_YMRS = EMAX_VPA_YMRS * pow(C_VPA, GAMMA_VPA_YMRS) / (pow(EC50_VPA_YMRS, GAMMA_VPA_YMRS) + pow(C_VPA, GAMMA_VPA_YMRS));

double EFFECT_QTP_D2    = EMAX_QTP_D2    * pow(C_QTP, GAMMA_QTP_D2)    / (pow(EC50_QTP_D2,    GAMMA_QTP_D2)    + pow(C_QTP, GAMMA_QTP_D2));
double EFFECT_QTP_YMRS  = EMAX_QTP_YMRS  * pow(C_QTP, GAMMA_QTP_YMRS)  / (pow(EC50_QTP_YMRS,  GAMMA_QTP_YMRS)  + pow(C_QTP, GAMMA_QTP_YMRS));
double EFFECT_QTP_MADRS = EMAX_QTP_MADRS * pow(C_QTP, GAMMA_QTP_MADRS) / (pow(EC50_QTP_MADRS, GAMMA_QTP_MADRS) + pow(C_QTP, GAMMA_QTP_MADRS));
double EFFECT_QTP_5HT   = EMAX_QTP_5HT   * pow(C_NQT, GAMMA_QTP_5HT)   / (pow(EC50_QTP_5HT,   GAMMA_QTP_5HT)   + pow(C_NQT, GAMMA_QTP_5HT));

double EFFECT_LTG_MADRS = EMAX_LTG_MADRS * pow(C_LTG, GAMMA_LTG_MADRS) / (pow(EC50_LTG_MADRS, GAMMA_LTG_MADRS) + pow(C_LTG, GAMMA_LTG_MADRS));

// ---- Combined GSK-3beta inhibition (Li + VPA, same nested formula as the
// original -- preserved exactly, including its double-negation, which
// algebraically makes GSK_total_inh equal (1-EFFECT_LI_GSK3)*(1-EFFECT_VPA_GSK3)
// rather than the "1 minus that product" its name would suggest; this is a
// pre-existing quirk of the original, not something introduced here -- a
// first attempt at a "cleaned up" one-level version (1.0 - (1-a)*(1-b))
// verified as a real, non-floating-point-scale mismatch against the
// original via the qspserver API, confirming the extra nesting below is
// load-bearing, not redundant. See bd_refactor_notes.md.) ----
double GSK_total_inh = 1.0 - (1.0 - (1.0 - EFFECT_LI_GSK3) * (1.0 - EFFECT_VPA_GSK3));
if (GSK_total_inh > 0.95) GSK_total_inh = 0.95;

// ---- Mood-state combined drivers (additive simplification, same as original) ----
double E_YMRS_total = EFFECT_LI_YMRS + EFFECT_VPA_YMRS + EFFECT_QTP_YMRS;

// Li+QTP synergy on MADRS (Geddes 2016 combination) -- thresholds (0.3, 50)
// kept as literal numbers exactly as the original (not named params there either)
double E_MADRS_combo_bonus = ALPHA_COMBO_LI_QTP * ((C_LI > 0.3 ? 1.0 : 0.0) * (C_QTP > 50 ? 1.0 : 0.0)) * 4.0;
double E_MADRS_total = EFFECT_QTP_MADRS + EFFECT_LTG_MADRS + E_MADRS_combo_bonus;

// ---- Circadian forcing (unchanged, not compound-specific) ----
// [Upstream build-defect fix, syntax-only -- see bd_refactor_notes.md /
// UPSTREAM_ISSUES.md]: the original writes `SOLVERTIME` here, but this is
// computed inside $MAIN, where SOLVERTIME’s underlying `_ODETIME_` symbol
// is only declared inside $ODE’s own generated scope -- a genuine mrgsolve
// 2.0.1 build defect (confirmed by a minimal reproduction), not a
// PK/PD-scope change. `TIME` is the $MAIN-scope equivalent and was
// confirmed (same minimal reproduction) to report the identical value as
// SOLVERTIME at every reporting row.
double Circ_forcing = Amp_circ * sin(omega * TIME);

// ---- Weight-gain driver (norquetiapine H1 effect; LINEAR in the original -- preserved as linear, not forced into a Hill shape) ----
double Wt_gain_rate = KWT_QTP * C_NQT;

// ---- Initial conditions ----
if (NEWIND <= 1) {
  DA_index_0    = DA_base;
  HT5_index_0   = HT5_base;
  GSK3_activ_0  = GSK3_base;
  BDNF_level_0  = BDNF_base;
  IL6_level_0   = IL6_base;
  Cortisol_idx_0 = Cort_base;
  Circ_state_0  = 1.0;
  Circ_deriv_0  = 0.0;
  GAF_score_0   = 50.0;  // moderate impairment at episode
}

$ODE
// === PK ODEs ===
// Lithium (archetype 3: depot + central + peripheral, linear)
dxdt_GUT_LI  = -KA_LI * GUT_LI;
dxdt_CENT_LI = KA_LI * GUT_LI
               - (CL_LI/V1_LI + Q_LI/V1_LI) * CENT_LI
               + Q_LI/V2_LI * PERI_LI;
dxdt_PERI_LI = Q_LI/V1_LI * CENT_LI - Q_LI/V2_LI * PERI_LI;

// Valproate (archetype 3 minus peripheral: depot + central, linear elimination on total amount)
dxdt_GUT_VPA  = -KA_VPA * GUT_VPA;
dxdt_CENT_VPA = KA_VPA * GUT_VPA - (CL_VPA/V1_VPA) * CENT_VPA;

// Quetiapine + norquetiapine
dxdt_GUT_QTP  = -KA_QTP * GUT_QTP;
dxdt_CENT_QTP = KA_QTP * F_QTP * GUT_QTP
                - (CL_QTP / (V1_QTP/1000.0)) * CENT_QTP;
dxdt_CENT_NQT = FMET_QTP * (CL_QTP / (V1_QTP/1000.0)) * CENT_QTP
                - (CL_NQT / (V1_QTP/1000.0)) * CENT_NQT;

// Lamotrigine (archetype 3 minus peripheral: depot + central, linear)
dxdt_GUT_LTG  = -KA_LTG * GUT_LTG;
dxdt_CENT_LTG = KA_LTG * GUT_LTG - (CL_LTG/V1_LTG) * CENT_LTG;

// === Neurotransmitter PDEs ===
// Dopamine index (QTP suppresses via D2R)
dxdt_DA_index = kDA_syn * DA_base - kDA_deg * DA_index
                - EFFECT_QTP_D2 * DA_index;

// Serotonin index (NQT/QTP via SERT/5HT2A rebound)
dxdt_HT5_index = k5HT_syn * HT5_base * (1.0 + EFFECT_QTP_5HT)
                 - k5HT_deg * HT5_index;

// === Signal transduction ===
// GSK-3beta activity (reduced by Li + VPA)
dxdt_GSK3_activ = kGSK_syn * GSK3_base * (1.0 - GSK_total_inh)
                  - kGSK_deg * GSK3_activ;

// BDNF (upregulated by Li, reduced by active GSK3)
dxdt_BDNF_level = kBDNF_syn * BDNF_base * (1.0 + EFFECT_LI_BDNF)
                  * (1.0 - 0.40 * (GSK3_activ / GSK3_base - 1.0) * (GSK3_activ > GSK3_base ? 1.0 : 0.0))
                  - kBDNF_deg * BDNF_level;

// Neuroinflammation (IL-6 suppressed by Li)
dxdt_IL6_level = kIL6_syn * IL6_base * (1.0 - EFFECT_LI_IL6)
                 + 0.10 * (GSK3_activ / GSK3_base)
                 - kIL6_deg * IL6_level;

// === HPA axis (cortisol with circadian modulation) ===
dxdt_Cortisol_idx = kCort_prod * (1.0 + Circ_forcing)
                    - kCort_deg * Cortisol_idx;

// === Circadian oscillator (van der Pol-like) ===
dxdt_Circ_state = Circ_deriv;
dxdt_Circ_deriv = -pow(omega, 2.0) * Circ_state
                  + omega * (1.0 - pow(Circ_state, 2.0)) * Circ_deriv;

// === Mood State ODEs ===
// YMRS (mania: driven up by high DA, reduced by drugs)
dxdt_YMRS_score = YMRS_score * (DA_index / DA_base - 1.0) * 0.05
                  - kYMRS_nat * YMRS_score
                  - E_YMRS_total * (YMRS_score / (YMRS_base + 1e-6)) * kYMRS_nat * 50.0;

// MADRS (depression: driven up by high cortisol/IL6, low BDNF)
dxdt_MADRS_score = MADRS_score * ((IL6_level / IL6_base - 1.0) * 0.05
                   + (Cortisol_idx / Cort_base - 1.0) * 0.03
                   + (1.0 - BDNF_level / BDNF_base) * 0.04)
                   - kMADRS_nat * MADRS_score
                   - E_MADRS_total * (MADRS_score / (MADRS_base + 1e-6)) * kMADRS_nat * 50.0;

// === Metabolic: Weight change ===
dxdt_Weight_chg = Wt_gain_rate;

// === Functioning (GAF) - inverse of severity ===
dxdt_GAF_score = 0.01 * (70.0 - GAF_score)    // natural recovery toward 70
                 - 0.05 * (YMRS_score + MADRS_score) / 30.0 * GAF_score * 0.01;

$TABLE
// ---- Derived outputs (redirected to the renamed compartments/params via
// the same $GLOBAL macros as $MAIN -- one authored formula per compound,
// re-declared here as its own $TABLE-local so it keeps the original’s own
// per-row-fresh reporting timing; output names unchanged for parity with
// the original’s own $CAPTURE) ----
double Lithium_mEqL  = LI_CONC_EXPR;
double VPA_ugmL      = VPA_TOTAL_EXPR;
double VPA_free_ugmL = VPA_FREE_EXPR(VPA_ugmL);
double QTP_ngmL      = QTP_CONC_EXPR;
double NQT_ngmL      = NQT_CONC_EXPR;
double LTG_ugmL      = LTG_CONC_EXPR;

// Safety flags
double Li_toxic = (Lithium_mEqL > 1.5) ? 1.0 : 0.0;
double Li_subtherapeutic = (Lithium_mEqL < 0.5 && Lithium_mEqL > 0.0) ? 1.0 : 0.0;
double VPA_toxic = (VPA_ugmL > 125) ? 1.0 : 0.0;

// Response flags
double YMRS_response   = (YMRS_score <= YMRS_base * 0.5)  ? 1.0 : 0.0;
double YMRS_remission  = (YMRS_score <= 12.0)             ? 1.0 : 0.0;
double MADRS_response  = (MADRS_score <= MADRS_base * 0.5) ? 1.0 : 0.0;
double MADRS_remission = (MADRS_score <= 12.0)             ? 1.0 : 0.0;

$CAPTURE
// [Upstream build-defect fix, syntax-only -- see bd_refactor_notes.md /
// UPSTREAM_ISSUES.md]: the original re-lists ten of its own $CMT
// compartment names directly in $CAPTURE (DA_index, HT5_index, GSK3_activ,
// BDNF_level, IL6_level, Cortisol_idx, YMRS_score, MADRS_score, GAF_score,
// Weight_chg) -- mrgsolve 2.0.1 rejects a compartment name inside
// $CAPTURE. Dropped here; every one of the ten is still reported, unchanged,
// since mrgsolve exposes $CMT states as output columns automatically.
Lithium_mEqL VPA_ugmL VPA_free_ugmL QTP_ngmL NQT_ngmL LTG_ugmL
Li_toxic Li_subtherapeutic VPA_toxic
YMRS_response YMRS_remission MADRS_response MADRS_remission
C_LI C_VPA C_QTP C_NQT C_LTG
EFFECT_LI_GSK3 EFFECT_LI_BDNF EFFECT_LI_IL6 EFFECT_LI_YMRS
EFFECT_VPA_GSK3 EFFECT_VPA_YMRS
EFFECT_QTP_D2 EFFECT_QTP_YMRS EFFECT_QTP_MADRS EFFECT_QTP_5HT
EFFECT_LTG_MADRS
'

## ---------------------------------------------------------------------------
## Compile the model
## ---------------------------------------------------------------------------
mod <- mcode("BipolarDisorder_QSP_refactored", code_refactored, quiet = TRUE)

## ---------------------------------------------------------------------------
## Helper: dosing event builder (compartment names updated to the renamed $CMT block)
## ---------------------------------------------------------------------------
make_dosing <- function(drug, dose, interval_h, n_doses, start_time = 0,
                        cmt_name = NULL) {
  cmt_map <- list(
    "Lithium"     = "GUT_LI",
    "Valproate"   = "GUT_VPA",
    "Quetiapine"  = "GUT_QTP",
    "Lamotrigine" = "GUT_LTG"
  )
  cmt <- if (!is.null(cmt_name)) cmt_name else cmt_map[[drug]]
  ev(
    amt  = dose,
    cmt  = cmt,
    ii   = interval_h,
    addl = n_doses - 1,
    time = start_time
  )
}

## ---------------------------------------------------------------------------
## SCENARIO 1: Lithium Monotherapy – Acute Mania (BD-I)
##   Dose: 300 mg TID (standard carbonate), targeting 0.8–1.2 mEq/L
##   Duration: 21 days   Reference: Bowden et al. AJP 1994
## ---------------------------------------------------------------------------
# Note: 300 mg Li₂CO₃ ≈ 8.1 mmol Li  (MW Li=6.94, Li₂CO₃=73.9; 300 mg gives ~2×8.1 mEq)
# Simplified: 1 dose unit = 8.1 mmol Li (stored in GUT_LI as mmol)

e1 <- make_dosing("Lithium", dose = 8.1, interval_h = 8, n_doses = 63,
                  start_time = 0)  # 21 days TID
idata1 <- tibble(ID = 1, YMRS_base = 25, MADRS_base = 5)

out1 <- mod %>%
  param(YMRS_base = 25, MADRS_base = 5) %>%
  idata_set(idata1) %>%
  ev(e1) %>%
  mrgsim(end = 504, delta = 1) %>%  # 21 days, hourly
  as_tibble()

cat("\n=== Scenario 1: Lithium Monotherapy – Acute Mania ===\n")
cat(sprintf("Peak Li: %.2f mEq/L | Day-21 YMRS: %.1f | Response: %s\n",
            max(out1$Lithium_mEqL),
            tail(out1$YMRS_score, 1),
            ifelse(tail(out1$YMRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 2: Valproate Monotherapy – Acute Mania (BD-I / Mixed)
##   Dose: 500 mg BID → 1000 mg/day (target 50–100 μg/mL)
##   Reference: Bowden 1994; Pope 1991 RCT
## ---------------------------------------------------------------------------
e2 <- make_dosing("Valproate", dose = 500, interval_h = 12, n_doses = 42)
out2 <- mod %>%
  param(YMRS_base = 25, MADRS_base = 5) %>%
  ev(e2) %>%
  mrgsim(end = 504, delta = 1) %>%
  as_tibble()

cat("\n=== Scenario 2: Valproate Monotherapy – Acute Mania ===\n")
cat(sprintf("Peak VPA: %.1f μg/mL (free: %.1f) | Day-21 YMRS: %.1f\n",
            max(out2$VPA_ugmL),
            max(out2$VPA_free_ugmL),
            tail(out2$YMRS_score, 1)))

## ---------------------------------------------------------------------------
## SCENARIO 3: Quetiapine Monotherapy – Bipolar Depression
##   Dose: 300 mg QD (BOLDER I/II, EMBOLDEN I/II target)
##   Reference: Calabrese et al. AJP 2005; Young et al. 2010
## ---------------------------------------------------------------------------
e3 <- make_dosing("Quetiapine", dose = 300, interval_h = 24, n_doses = 56)
out3 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e3) %>%
  mrgsim(end = 1344, delta = 1) %>%  # 56 days
  as_tibble()

cat("\n=== Scenario 3: Quetiapine – Bipolar Depression ===\n")
cat(sprintf("Mean QTP Css: %.0f ng/mL | Day-56 MADRS: %.1f | Response: %s\n",
            mean(tail(out3$QTP_ngmL, 24)),
            tail(out3$MADRS_score, 1),
            ifelse(tail(out3$MADRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 4: Lithium + Quetiapine Combination – BD Depression
##   Li 300 mg TID + QTP 300 mg QD
##   Reference: CANMAT 2018; Geddes et al. Lancet 2016
## ---------------------------------------------------------------------------
e4a <- make_dosing("Lithium",    dose = 8.1, interval_h = 8,  n_doses = 168) # 56d
e4b <- make_dosing("Quetiapine", dose = 300, interval_h = 24, n_doses = 56)

out4 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e4a + e4b) %>%
  mrgsim(end = 1344, delta = 1) %>%
  as_tibble()

cat("\n=== Scenario 4: Lithium + Quetiapine – BD Depression ===\n")
cat(sprintf("Li Css: %.2f mEq/L | Day-56 MADRS: %.1f | Remission: %s\n",
            mean(tail(out4$Lithium_mEqL, 24)),
            tail(out4$MADRS_score, 1),
            ifelse(tail(out4$MADRS_remission, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 5: Lithium Maintenance (1-year prevention)
##   Li 300 mg TID (0.6–0.8 mEq/L maintenance range)
##   Reference: Cipriani et al. Lancet 2013 meta-analysis
## ---------------------------------------------------------------------------
e5 <- make_dosing("Lithium", dose = 8.1, interval_h = 8, n_doses = 3 * 365)
out5 <- mod %>%
  param(YMRS_base = 10, MADRS_base = 10) %>%  # stable residual symptoms
  ev(e5) %>%
  mrgsim(end = 8760, delta = 4) %>%   # 1 year, 4-h intervals
  as_tibble()

cat("\n=== Scenario 5: Lithium Maintenance (1 year) ===\n")
cat(sprintf("Steady-state Li: %.2f mEq/L | Year-end BDNF: %.2f | GSK3: %.2f\n",
            mean(tail(out5$Lithium_mEqL, 100)),
            mean(tail(out5$BDNF_level, 100)),
            mean(tail(out5$GSK3_activ, 100))))

## ---------------------------------------------------------------------------
## SCENARIO 6: Lamotrigine Add-on for BD-II Depression
##   LTG titrated: 25 mg/d wk1-2, 50 mg/d wk3-4, 100 mg/d wk5-6, 200 mg/d
##   Reference: Calabrese et al. JAMA 1999; STRIDE-BD
## ---------------------------------------------------------------------------
# 4-week titration (simplified)
e6a <- make_dosing("Lamotrigine", dose = 25,  interval_h = 24, n_doses = 14)           # wk1-2
e6b <- make_dosing("Lamotrigine", dose = 50,  interval_h = 24, n_doses = 14, start_time = 336)  # wk3-4
e6c <- make_dosing("Lamotrigine", dose = 100, interval_h = 24, n_doses = 14, start_time = 672)  # wk5-6
e6d <- make_dosing("Lamotrigine", dose = 200, interval_h = 24, n_doses = 56, start_time = 1008) # wk7-14

out6 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e6a + e6b + e6c + e6d) %>%
  mrgsim(end = 2688, delta = 2) %>%   # 112 days
  as_tibble()

cat("\n=== Scenario 6: Lamotrigine Titration – BD-II Depression ===\n")
cat(sprintf("Steady-state LTG: %.2f μg/mL | Day-112 MADRS: %.1f | Response: %s\n",
            mean(tail(out6$LTG_ugmL, 48)),
            tail(out6$MADRS_score, 1),
            ifelse(tail(out6$MADRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## PLOTTING: All scenarios comparison
## ---------------------------------------------------------------------------
library(ggplot2)
library(tidyr)

plot_comparison <- function() {
  # Mania panel: Scenarios 1 & 2
  d_mania <- bind_rows(
    out1 %>% select(time, YMRS_score, Lithium_mEqL) %>% mutate(Scenario = "1: Li monotherapy"),
    out2 %>% select(time, YMRS_score, VPA_ugmL)     %>% mutate(Scenario = "2: VPA monotherapy")
  )

  p1 <- ggplot(d_mania, aes(time / 24, YMRS_score, color = Scenario)) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
    annotate("text", x = 18, y = 13.5, label = "Remission (YMRS≤12)", size = 3.5) +
    labs(title = "Acute Mania: YMRS over Time",
         x = "Day", y = "YMRS Score", color = "Treatment") +
    theme_bw(base_size = 13) +
    scale_color_brewer(palette = "Set1")

  # Depression panel: Scenarios 3, 4, 6
  d_dep <- bind_rows(
    out3 %>% select(time, MADRS_score) %>% mutate(Scenario = "3: QTP monotherapy"),
    out4 %>% select(time, MADRS_score) %>% mutate(Scenario = "4: Li + QTP"),
    out6 %>% select(time, MADRS_score) %>% mutate(Scenario = "6: LTG titration")
  )

  p2 <- ggplot(d_dep, aes(time / 24, MADRS_score, color = Scenario)) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
    annotate("text", x = 40, y = 13.5, label = "Remission (MADRS≤12)", size = 3.5) +
    labs(title = "Bipolar Depression: MADRS over Time",
         x = "Day", y = "MADRS Score", color = "Treatment") +
    theme_bw(base_size = 13) +
    scale_color_brewer(palette = "Dark2")

  # PK panel: Lithium concentrations (safety window)
  p3 <- ggplot(out1, aes(time / 24, Lithium_mEqL)) +
    geom_line(color = "#e41a1c", size = 1.1) +
    geom_ribbon(aes(ymin = 0.6, ymax = 1.2), alpha = 0.15, fill = "green4") +
    geom_hline(yintercept = 1.5, linetype = "dashed", color = "red3") +
    annotate("text", x = 15, y = 1.55, label = "Toxic threshold (1.5 mEq/L)", color = "red3", size = 3.5) +
    labs(title = "Lithium PK – Concentration over Time",
         x = "Day", y = "Serum Li (mEq/L)") +
    theme_bw(base_size = 13)

  # Biomarker panel: BDNF, GSK-3beta under Li maintenance
  d_bio <- out5 %>%
    select(time, BDNF_level, GSK3_activ, IL6_level) %>%
    pivot_longer(-time, names_to = "Biomarker", values_to = "Value") %>%
    mutate(time_day = time / 24)

  p4 <- ggplot(d_bio, aes(time_day, Value, color = Biomarker)) +
    geom_line(size = 0.9) +
    labs(title = "Biomarker Dynamics – Li Maintenance (1 year)",
         x = "Day", y = "Normalised Index", color = "Biomarker") +
    theme_bw(base_size = 13) +
    scale_color_manual(values = c("BDNF_level" = "green4",
                                   "GSK3_activ" = "firebrick",
                                   "IL6_level"  = "steelblue"))

  list(mania = p1, depression = p2, pk_Li = p3, biomarkers = p4)
}

plots <- plot_comparison()

## Print summary table
summary_tbl <- tibble(
  Scenario = c("1: Li mono (mania)", "2: VPA mono (mania)", "3: QTP mono (BDdep)",
               "4: Li+QTP (BDdep)", "5: Li maintenance", "6: LTG titration (BDdep)"),
  Duration_days = c(21, 21, 56, 56, 365, 112),
  Peak_drug = c(
    round(max(out1$Lithium_mEqL),   2),
    round(max(out2$VPA_ugmL),       1),
    round(max(out3$QTP_ngmL),       0),
    round(mean(tail(out4$Lithium_mEqL, 24)), 2),
    round(mean(tail(out5$Lithium_mEqL, 100)), 2),
    round(max(out6$LTG_ugmL),       2)
  ),
  Drug_unit = c("mEq/L", "μg/mL", "ng/mL", "mEq/L", "mEq/L", "μg/mL"),
  End_YMRS  = c(tail(out1$YMRS_score,1), tail(out2$YMRS_score,1), NA, NA,
                tail(out5$YMRS_score,1), NA),
  End_MADRS = c(NA, NA, tail(out3$MADRS_score,1), tail(out4$MADRS_score,1),
                NA, tail(out6$MADRS_score,1)),
  Response  = c(tail(out1$YMRS_response,1),  tail(out2$YMRS_response,1),
                tail(out3$MADRS_response,1),  tail(out4$MADRS_remission,1),
                NA, tail(out6$MADRS_response,1))
)

print(summary_tbl)
cat("\nModel compiled successfully. Use plots$mania, plots$depression, etc. to view results.\n")
