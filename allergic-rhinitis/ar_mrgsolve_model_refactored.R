################################################################################
# Allergic Rhinitis (AR) - Quantitative Systems Pharmacology (QSP) Model
# mrgsolve implementation (REFACTORED, pluggable PK)
#   Refactor of `ar_mrgsolve_model.R` per FORK_WORKFLOW_GUIDE.md Part 2.
#   Original untouched - see `ar_refactor_notes.md` for full rationale.
#
#   Census scope: this file carries four rows in
#   `driver-patches/data/compound_perturbation_census.md`, all labeled by the
#   automated classifier as "Symptom model (<TOKEN>)" -- a classifier
#   artifact pulled from a nearby comment, not a real compound name (per the
#   task's own warning about this failure mode). Checked against the actual
#   code and the directory README: all four are genuine, named, externally
#   dosed drugs, not a generic "symptom model" construct --
#     - "Symptom model (CETI)"     -> Cetirizine     (oral H1-antihistamine)
#     - "Symptom model (FP_LOCAL)" -> Fluticasone propionate (intranasal
#                                     corticosteroid, local nasal mucosal PK)
#     - "Symptom model (MLKT)"     -> Montelukast    (oral CysLT1 antagonist)
#     - "Symptom model (OMA)"      -> Omalizumab     (SC anti-IgE mAb)
#   Each has its own depot/central PK compartment(s), is dosed by its own
#   `ev()` event in the original's own `build_doses()`/scenario list, and
#   gates a named downstream pharmacodynamic effect. See
#   `ar_refactor_notes.md` for the per-compound identity check and archetype.
#
#   Build-compatibility fixes (disclosed, syntax-only, non-numeric -- see
#   notes and UPSTREAM_ISSUES.md): the ORIGINAL does not compile under
#   mrgsolve 2.0.1, for two independent reasons, both fixed here:
#     1. `$INIT @annotated` used `NAME = value : description` (an `=` sign)
#        for every one of its 22 entries; mrgsolve's annotated-block parser
#        requires `NAME : value : description` (colon-separated), the same
#        convention already used correctly in this file's own `$PARAM`
#        block. Switched all 22 `$INIT` lines to the colon form -- same
#        names, same values, same descriptions.
#     2. Five `capture X Y Z ...` lines (space-separated, no parentheses, no
#        semicolon) sit inside `$TABLE` as if they were `$CAPTURE`-style
#        annotations; mrgsolve compiles `$TABLE` as literal C++, where this
#        is not valid syntax at all ("expected initializer before ...").
#        Replaced with a single proper `$CAPTURE` block listing the same
#        output names (compartment names dropped from the list per mrgsolve
#        2.0.1's separate "compartment should not be in $CAPTURE" rule --
#        already-exposed compartments IL4/IL5/IL13/HISTAMINE/CYS_LT/
#        MAST_ACT/EOS_N do not need re-listing).
#   Both fixes are confirmed syntax-only/non-numeric by the verification
#   below (exact match against a patched-original scratch copy carrying only
#   these same two fixes).
#
#   Dose-instant reporting artifact (found and fixed, not just disclosed):
#   verification against the patched-original scratch copy via the live
#   qspserver mrgsolve API initially showed FP_LOCAL_NM/GR_OCC_FP diverging
#   at the t=0 duplicate report row of every Fluticasone dose (abs diff up
#   to 450, the full dose amount) while every other output matched exactly.
#   Root cause: Fluticasone doses directly into CENT_FP (no depot), and
#   C_FP/EFFECT_FP were $ODE-scope `double`s -- the same stale-pre-dose-read
#   quirk documented for clostridioides-difficile-infection's bezlotoxumab
#   (also dosed with no depot). Fixed the same way: C_FP/EFFECT_FP are now
#   $GLOBAL `#define` macros (see top of the DSL block below), re-expanded
#   textually at every point of use instead of cached once per $ODE call --
#   confirmed via re-run through qspserver mrgsolve_api that this eliminates
#   the artifact entirely (see ar_refactor_notes.md for the before/after
#   values). Same disclosed discoverability trade-off as cdi's C_VAN/C_BEZ:
#   a `#define`, not a literal `double C_FP = <expr>;` statement, so a naive
#   `double C_<STEM> = ` text grep will not find this one compound's site
#   (EC50_FP is still a normal $PARAM entry) -- see ar_refactor_notes.md.
#
#   PK/PD  : Cetirizine (oral H1-antihistamine, QD) - Fluticasone propionate
#            (intranasal corticosteroid, local nasal mucosal PK, QD) -
#            Montelukast (oral CysLT1 receptor antagonist, QD) - Omalizumab
#            (SC anti-IgE monoclonal antibody, q4w)
#   Renamed to this fork's PK/PD naming convention: GUT_<STEM>/CENT_<STEM>/
#   PERI_<STEM>, CL_<STEM>/V1_<STEM>/V2_<STEM>/Q_<STEM>/KA_<STEM>/F_<STEM>,
#   C_<STEM> (exposed concentration), EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>/
#   EFFECT_<STEM> (Hill interface). See `ar_refactor_notes.md` for the
#   per-compound archetype (Cetirizine/Montelukast: Archetype 3 minus
#   peripheral; Fluticasone: bespoke single-compartment with a disclosed,
#   preserved-not-fixed constant zero-order input; Omalizumab: Archetype 3
#   [depot-central-peripheral] PK + a bespoke TMDD-style IgE-capture term,
#   since the "receptor" it binds is the disease model's own shared
#   IGE_FREE state, not an independent compound-exclusive pool).
#   Calibration targets (literature, unchanged from original):
#     - Cetirizine H1-RO >=80% at steady state  [Yanai 1995 JACI]
#     - FP reduces TNSS ~35-40% vs placebo      [Meltzer 2005 JACI]
#     - Montelukast reduces TNSS ~20-25%        [Philip 2002 JACI]
#     - Omalizumab reduces free IgE >95%        [Fahy 1997 AJRCCM]
#     - Nasal eosinophilia reduced ~50% by FP   [Holgate 2003 Allergy]
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# Model code block
# ============================================================

AR_model_refactored_code <- '
$PROB
  Allergic Rhinitis QSP Model -- PK/PD refactor (pluggable PK)
  IgE/Mast Cell . Th2 Cytokines . Eosinophils . PK/PD
  Cetirizine (CETI) | Fluticasone propionate (FP) | Montelukast (MLKT) | Omalizumab (OMA)
  Renamed to this fork’s naming convention: GUT_<STEM>/CENT_<STEM>/PERI_<STEM>,
  C_<STEM> (exposed concentration), EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>/EFFECT_<STEM>
  (Hill interface). See ar_refactor_notes.md for per-compound archetype detail.

$GLOBAL
// C_FP/EFFECT_FP are $GLOBAL macros, not $ODE-scope doubles, following the
// clostridioides-difficile-infection precedent (see cdi_mrgsolve_model_
// refactored.R and cdi_refactor_notes.md): Fluticasone is dosed directly
// into CENT_FP (no depot), so an $ODE-scope double reads a stale PRE-dose
// value on the duplicate report row mrgsolve/qspserver emits at the exact
// instant of a dose -- confirmed live via qspserver mrgsolve_api (see
// ar_refactor_notes.md for the before/after comparison). A macro is
// re-expanded textually wherever referenced (including in $TABLE at ITS
// own execution time), so it always reads the current, already-dosed
// state and eliminates the artifact entirely rather than shrinking it.
// This affects only the diagnostic/reporting duplicate-row value, never
// the actual ODE state trajectory.
#define C_FP      (CENT_FP / V1_FP)
#define EFFECT_FP (EMAX_FP * pow(C_FP, GAMMA_FP) / (pow(EC50_FP, GAMMA_FP) + pow(C_FP, GAMMA_FP)))

$PARAM @annotated
  // --- Allergen ---
  ALLERGEN_SS  : 1.0   : Steady-state nasal allergen load (AU/mL, normalized)
  K_ALLERGEN   : 0.1   : Allergen clearance rate (1/h)

  // --- IgE dynamics ---
  KSY_IGE      : 0.002 : IgE synthesis rate (AU/h, baseline)
  KDEG_IGE     : 0.005 : Free IgE degradation rate (1/h, t1/2~138h)
  K_BIND_MAST  : 0.1   : IgE binding to mast cell FceRI (1/h)
  K_OFF_MAST   : 0.001 : IgE dissociation from FceRI (1/h)
  MAST_TOTAL   : 1.0   : Total mast cell FceRI capacity (AU, normalized)

  // --- Mast cell activation & degranulation ---
  EC50_CROSS   : 0.5   : Allergen EC50 for crosslinking IgE-Ag (AU)
  HILL_CROSS   : 2.0   : Hill coefficient for crosslinking
  KDEG_HIST    : 2.0   : Histamine degradation (1/h, DAO/HNMT)
  KDEG_LT      : 1.0   : CysLT degradation (1/h)
  KHIST_PROD   : 5.0   : Max histamine production per unit mast activation
  KLT_PROD     : 2.0   : Max CysLT production per unit mast activation
  KMAST_REC    : 0.05  : Mast cell recharging rate after degranulation (1/h)

  // --- Th2 cytokines (IL-4, IL-5, IL-13) ---
  KSY_TH2      : 0.05  : Th2 cell activation rate driven by IgE/allergen
  KDEG_TH2     : 0.1   : Th2 cell decay (1/h)
  KSY_IL4      : 1.0   : IL-4 synthesis per Th2 unit (pg/mL/h)
  KSY_IL5      : 0.8   : IL-5 synthesis per Th2 unit
  KSY_IL13     : 1.2   : IL-13 synthesis per Th2 unit
  KDEG_IL4     : 0.5   : IL-4 clearance (1/h)
  KDEG_IL5     : 0.3   : IL-5 clearance (1/h)
  KDEG_IL13    : 0.4   : IL-13 clearance (1/h)
  TH2_BASE     : 0.5   : Baseline Th2 activity (atopic, normalized)

  // --- Eosinophil dynamics ---
  EOS_BLOOD_0  : 300   : Baseline blood eosinophil (cells/uL)
  KEO_PROD     : 0.05  : Eosinophil production stimulated by IL-5 (cells/uL/h)
  KEO_SURV     : 0.02  : IL-5-mediated eosinophil survival boost (1/h)
  KEO_DEATH    : 0.08  : Eosinophil apoptosis rate (1/h, t1/2~8.7h)
  KEO_MIGRATE  : 0.01  : Blood->nasal tissue migration (1/h) per chemokine unit
  KCHEMOKINE   : 2.0   : Eotaxin/CCL11 driven chemotaxis coefficient
  KEO_TISSUE0  : 10.0  : Baseline nasal tissue eosinophil (cells/uL tissue)
  KEO_TIS_DEATH: 0.05  : Nasal tissue eosinophil apoptosis (1/h)

  // --- Symptom model ---
  HIST_EC50    : 1.0   : Histamine EC50 for symptom generation
  LT_EC50      : 0.5   : CysLT EC50 for congestion
  EOS_EC50     : 50.0  : Nasal eosinophil EC50 for late symptoms
  SNEEZE_MAX   : 3.0   : Maximum sneezing score
  RHINO_MAX    : 3.0   : Maximum rhinorrhea score
  CONG_MAX     : 3.0   : Maximum congestion score
  PRUR_MAX     : 3.0   : Maximum pruritus score

  // =====================================================
  // PK -- Cetirizine (CETI): oral, depot->central, linear elim
  // Archetype 3 minus peripheral. Tmax ~1h, t1/2 ~10h, Vd 0.56 L/kg
  // =====================================================
  KA_CETI      : 0.9   : Cetirizine absorption rate (1/h)
  CL_CETI      : 7.0   : Cetirizine clearance (L/h)
  V1_CETI      : 70.0  : Cetirizine central volume (L)                 // was VD_CETI
  F_CETI       : 0.70  : Oral bioavailability cetirizine
  EC50_CETI    : 15.0  : Cetirizine EC50 for H1R occupancy (ng/mL)      // was H1_IC50_CETI
  GAMMA_CETI   : 1.0   : Hill coefficient, H1R occupancy                // was H1_HILL
  EMAX_CETI    : 1.0   : Hill ceiling [math-implied: original H1_RO ratio saturates at 1]

  // =====================================================
  // PK -- Fluticasone propionate (FP): intranasal local mucosal
  // Bespoke single-compartment (no depot; see refactor notes for the
  // preserved-not-fixed constant zero-order input term). Systemic exposure <2%.
  // =====================================================
  KA_FP        : 1.5   : FP nasal mucosal absorption rate constant (1/h)
  CL_FP        : 2.0   : FP local elimination rate constant from nasal mucosa (1/h)  // was CL_FP_LOCAL
  V1_FP        : 1.0   : FP mucosal volume (L, local compartment)       // was VD_FP_LOCAL
  BIOCONV_FP   : 200.0 : Intranasal dose bioconversion constant (200 ug ~ 450 nM local) // was FP_DOSE_BIOCONV
  EC50_FP      : 0.5   : FP EC50 for GR-mediated cytokine suppression (nM) // was GR_IC50_FP
  GAMMA_FP     : 1.2   : Hill coefficient, GR occupancy                  // was GR_HILL
  EMAX_FP      : 1.0   : Hill ceiling [math-implied: original GR_OCC ratio saturates at 1]

  // =====================================================
  // PK -- Montelukast (MLKT): oral, depot->central, linear elim
  // Archetype 3 minus peripheral. Tmax 3-4h, t1/2 5.5h, BA ~64%
  // =====================================================
  KA_MLKT      : 0.5   : Montelukast absorption rate (1/h)
  CL_MLKT      : 45.0  : Montelukast clearance (L/h)
  V1_MLKT      : 10.0  : Montelukast central volume (L)                 // was VD_MLKT
  F_MLKT       : 0.64  : Oral bioavailability montelukast
  EC50_MLKT    : 2.0   : Montelukast EC50, CysLT1 blockade (ng/mL)       // was CYSLTR1_IC50
  GAMMA_MLKT   : 1.0   : Hill coefficient, CysLT1 blockade [no explicit Hill term in original -> gamma=1]
  EMAX_MLKT    : 1.0   : Hill ceiling [math-implied: original CYSLTR1_INH ratio saturates at 1]

  // =====================================================
  // PK -- Omalizumab (OMA): SC, depot->central->peripheral (Archetype 3)
  // + bespoke TMDD-style IgE capture (binds the disease model’s own
  // shared IGE_FREE state, not an independent compound-exclusive pool)
  // t1/2 ~26d, ka ~0.004/h, Vss ~78 mL/kg
  // =====================================================
  KA_OMA       : 0.004 : Omalizumab SC absorption rate (1/h)
  CL_OMA       : 0.14  : Omalizumab central clearance (mL/h)
  V1_OMA       : 3140  : Omalizumab central volume (mL, ~45 mL/kg x 70kg)  // was VC_OMA
  V2_OMA       : 2360  : Omalizumab peripheral volume (mL)                 // was VP_OMA
  Q_OMA        : 0.40  : Intercompartmental clearance (mL/h)
  KON_OMA      : 1.0   : Omalizumab-IgE association rate                  // was KON_IGE
  KOFF_OMA     : 0.0001: Omalizumab-IgE dissociation (KD ~0.1 nM)         // was KOFF_IGE

$INIT @annotated
  // Allergen
  AG       : 0.0   : Nasal allergen (AU)

  // IgE compartments
  IGE_FREE : 50.0  : Free IgE (IU/mL equivalent)
  IGE_MAST : 0.5   : Mast cell-bound IgE (normalized)

  // Mast cell state
  MAST_ACT : 0.0   : Mast cell activation level (0-1)
  MAST_CHG : 1.0   : Mast cell granule charge (0-1)

  // Mediators
  HISTAMINE : 0.0  : Nasal histamine (normalized AU)
  CYS_LT   : 0.0   : Cysteinyl leukotrienes (AU)

  // Th2 / cytokines
  TH2      : 0.5   : Th2 cell activity (normalized, atopic baseline)
  IL4      : 5.0   : IL-4 (pg/mL)
  IL5      : 3.0   : IL-5 (pg/mL)
  IL13     : 8.0   : IL-13 (pg/mL)

  // Eosinophils
  EOS_B    : 300.0 : Blood eosinophil (cells/uL)
  EOS_N    : 10.0  : Nasal tissue eosinophil (cells/uL tissue)

  // Cetirizine PK (depot + central)
  GUT_CETI  : 0.0  : Cetirizine depot (mg)                 // was CETI_D
  CENT_CETI : 0.0  : Cetirizine central (mg)                // was CETI_C

  // Fluticasone PK (nasal local, single compartment)
  CENT_FP  : 0.0   : Fluticasone local nasal conc (nM-normalized) // was FP_LOC

  // Montelukast PK
  GUT_MLKT  : 0.0  : Montelukast depot (mg)                 // was MLKT_D
  CENT_MLKT : 0.0  : Montelukast central (mg)                // was MLKT_C

  // Omalizumab PK (depot + 2-comp + IgE complex)
  GUT_OMA   : 0.0  : Omalizumab SC depot (mg)                // was OMA_D
  CENT_OMA  : 0.0  : Omalizumab central (mg/mL x volume -> mg) // was OMA_C
  PERI_OMA  : 0.0  : Omalizumab peripheral (mg)               // was OMA_P
  COMPLEX_OMA : 0.0 : Omalizumab-IgE complex (IU-equivalents) // was OMA_IGE

$MAIN
  // Allergen challenge: set AG from input or steady-state
  // Dose event (compartment 1) delivers allergen pulse

$ODE
  // ==============================================================
  // 1. Allergen kinetics
  // ==============================================================
  double Ag = AG;
  dxdt_AG = -K_ALLERGEN * AG;  // cleared by mucociliary; replenished via dose event

  // ==============================================================
  // 2. IgE dynamics + Omalizumab (OMA) TMDD-style capture
  //    (bespoke: OMA binds the disease model’s own shared IGE_FREE
  //    state, not an independent compound-exclusive receptor pool --
  //    see refactor notes)
  // ==============================================================
  double degen_free = KDEG_IGE * IGE_FREE;
  double bind_rate  = K_BIND_MAST * IGE_FREE * (MAST_TOTAL - IGE_MAST);
  double unbind_rate= K_OFF_MAST  * IGE_MAST;

  // Omalizumab: exposed concentration and IgE-capture term
  double C_OMA = CENT_OMA / V1_OMA;  // mg/mL -- exposed concentration (canonical site)
  double ige_capture = KON_OMA * C_OMA * IGE_FREE - KOFF_OMA * COMPLEX_OMA;

  dxdt_IGE_FREE = KSY_IGE * 200.0 - degen_free - bind_rate + unbind_rate - ige_capture;
  dxdt_IGE_MAST = bind_rate - unbind_rate;
  dxdt_COMPLEX_OMA = ige_capture;

  // ==============================================================
  // 3. Mast cell crosslinking & activation
  // ==============================================================
  double crosslink_frac = pow(Ag, HILL_CROSS) / (pow(EC50_CROSS, HILL_CROSS) + pow(Ag, HILL_CROSS));
  double mast_trigger   = crosslink_frac * IGE_MAST * MAST_CHG;  // depends on bound IgE and granule charge

  dxdt_MAST_ACT = mast_trigger - 0.5 * MAST_ACT;
  dxdt_MAST_CHG = KMAST_REC * (1.0 - MAST_CHG) - mast_trigger * MAST_CHG;

  // ==============================================================
  // 4. Histamine & CysLT release/degradation
  // ==============================================================
  dxdt_HISTAMINE = KHIST_PROD * MAST_ACT - KDEG_HIST * HISTAMINE;
  dxdt_CYS_LT   = KLT_PROD   * MAST_ACT - KDEG_LT   * CYS_LT;

  // ==============================================================
  // 5. Th2 / cytokine cascade
  // ==============================================================
  double th2_drive = TH2_BASE + KSY_TH2 * Ag * IGE_FREE / (1.0 + IGE_FREE);
  dxdt_TH2  = th2_drive - KDEG_TH2 * TH2;
  dxdt_IL4  = KSY_IL4  * TH2 - KDEG_IL4  * IL4;
  dxdt_IL5  = KSY_IL5  * TH2 - KDEG_IL5  * IL5;
  dxdt_IL13 = KSY_IL13 * TH2 - KDEG_IL13 * IL13;

  // ==============================================================
  // 6. Eosinophil dynamics
  // ==============================================================
  double eos_prod    = KEO_PROD * IL5 + KEO_SURV * EOS_B * IL5 / (1.0 + IL5);
  double eos_death_b = KEO_DEATH * EOS_B;
  double chemokine   = KCHEMOKINE * IL5 * IL13;
  double eos_migrate = KEO_MIGRATE * chemokine * EOS_B;
  double eos_death_n = KEO_TIS_DEATH * EOS_N;

  dxdt_EOS_B = eos_prod - eos_death_b - eos_migrate;
  dxdt_EOS_N = eos_migrate - eos_death_n;

  // ==============================================================
  // 7. Drug PK
  // ==============================================================

  // 7a. Cetirizine (CETI) -- Archetype 3 minus peripheral
  dxdt_GUT_CETI  = -KA_CETI * GUT_CETI;
  dxdt_CENT_CETI =  KA_CETI * F_CETI * GUT_CETI - CL_CETI / V1_CETI * CENT_CETI;

  // 7b. Fluticasone propionate (FP) -- bespoke, single compartment.
  // Preserved verbatim from the original, INCLUDING its constant
  // dose-independent "zero-order input" term (KA_FP*BIOCONV_FP/V1_FP is a
  // fixed number that does not depend on any depot state) -- this is a
  // pre-existing behavioral quirk of the original, not introduced by this
  // refactor; disclosed, not fixed, in ar_refactor_notes.md and
  // UPSTREAM_ISSUES.md.
  dxdt_CENT_FP = KA_FP * BIOCONV_FP / V1_FP - CL_FP * CENT_FP;
  // CENT_FP additionally replenished directly by dose events (no separate depot)

  // 7c. Montelukast (MLKT) -- Archetype 3 minus peripheral
  dxdt_GUT_MLKT  = -KA_MLKT * GUT_MLKT;
  dxdt_CENT_MLKT =  KA_MLKT * F_MLKT * GUT_MLKT - CL_MLKT / V1_MLKT * CENT_MLKT;

  // 7d. Omalizumab (OMA) -- Archetype 3 (depot->central->peripheral); IgE
  // capture (ige_capture, C_OMA) computed above in block 2, matching where
  // the original computed it (needed live, every substep, by dxdt_IGE_FREE).
  dxdt_GUT_OMA  = -KA_OMA * GUT_OMA;
  dxdt_CENT_OMA =  KA_OMA * GUT_OMA - (CL_OMA/V1_OMA + Q_OMA/V1_OMA) * CENT_OMA + Q_OMA/V2_OMA * PERI_OMA;
  dxdt_PERI_OMA =  Q_OMA/V1_OMA * CENT_OMA - Q_OMA/V2_OMA * PERI_OMA;

  // ==============================================================
  // 8. Drug pharmacodynamic effects (inhibitory) that feed back into
  //    the live ODE system (must stay in $ODE, matching the original --
  //    see "keep a calculation in the block the original used it in").
  //    Cetirizine’s and Montelukast’s Hill terms do NOT feed back into
  //    any dxdt_ here (same as the original, where the analogous
  //    ODE-block h1ro/cysltr1_inh locals were computed but never used by
  //    any dxdt_ line) -- both are single-sourced instead, computed once
  //    in $TABLE where their only real consumers (the symptom-score
  //    outputs) already lived. See ar_refactor_notes.md.
  // ==============================================================

  // Fluticasone: GR occupancy -> suppress cytokines (IL-4, IL-5, IL-13) & eosinophil.
  // C_FP/EFFECT_FP are $GLOBAL macros (see top of file), not doubles declared
  // here -- avoids the dose-instant stale-value reporting artifact while
  // still being available at $ODE’s live-integration time exactly as before
  // (nM; V1_FP = 1, so numerically identical to CENT_FP directly, matching
  // the original’s undivided use).

  // Apply FP effect on cytokines (transrepression)
  dxdt_IL4  -= EFFECT_FP * 0.6 * IL4;
  dxdt_IL5  -= EFFECT_FP * 0.7 * IL5;
  dxdt_IL13 -= EFFECT_FP * 0.6 * IL13;
  dxdt_EOS_N -= EFFECT_FP * 0.5 * EOS_N;

$TABLE
  // ==============================================================
  // Secondary PK variables -- exposed concentrations (canonical single site
  // per compound) and reporting-only unit conversions
  // ==============================================================
  double C_CETI      = CENT_CETI / V1_CETI * 1000.0;    // ng/mL (mg/L -> ng/mL)
  double C_MLKT       = CENT_MLKT / V1_MLKT * 1000.0;    // ng/mL
  double OMA_CP_REPORT = CENT_OMA / (V1_OMA/1000.0);     // ug/mL (reporting only; was OMA_CP)

  // Cetirizine: H1R occupancy Hill interface. Rename-only from the
  // original’s plain ratio (H1_RO = C^gamma/(EC50^gamma+C^gamma)*100);
  // EMAX_CETI/GAMMA_CETI/EC50_CETI carry the original’s values unchanged.
  double EFFECT_CETI = EMAX_CETI * pow(C_CETI, GAMMA_CETI) / (pow(EC50_CETI, GAMMA_CETI) + pow(C_CETI, GAMMA_CETI));

  // Montelukast: CysLT1 blockade Hill interface. Rename-only from the
  // original’s plain ratio (CYSLTR1_INH = C/(EC50+C)*100, implicit gamma=1).
  double EFFECT_MLKT = EMAX_MLKT * pow(C_MLKT, GAMMA_MLKT) / (pow(EC50_MLKT, GAMMA_MLKT) + pow(C_MLKT, GAMMA_MLKT));

  // Percentage-scale reporting variables (same names/values the original exposed)
  double H1_RO        = EFFECT_CETI * 100.0;   // %
  double GR_OCC_FP     = EFFECT_FP    * 100.0;   // %
  double CYSLTR1_INH  = EFFECT_MLKT  * 100.0;   // %
  double FP_LOCAL_NM  = C_FP;                    // nM (unchanged name/value)
  double CETI_CP      = C_CETI;                  // ng/mL (unchanged name/value)
  double MLKT_CP       = C_MLKT;                  // ng/mL (unchanged name/value)
  double OMA_CP        = OMA_CP_REPORT;           // ug/mL (unchanged name/value)

  // Free IgE reduction (omalizumab; vs baseline 50) -- descriptive-only
  // readout, not an algebraic Hill gate consumed elsewhere: omalizumab’s
  // actual PD action is already fully expressed by the TMDD-style
  // ige_capture ODE term above (see refactor notes).
  double IGE_REDUCTION_PCT = (1.0 - IGE_FREE/50.0) * 100.0;
  // Fractional bound-IgE occupancy, the TMDD-archetype analogue of
  // COMPLEX_<STEM>/RTOT_<STEM>: omalizumab has no fixed RTOT (total IgE is
  // itself a dynamic disease state, not a compound-exclusive constant
  // pool), so this uses the instantaneous free+bound total in its place.
  double EFFECT_OMA = COMPLEX_OMA / (IGE_FREE + COMPLEX_OMA + 1.0e-9);

  // ==============================================================
  // Drug-modified mediator outputs
  // ==============================================================
  double H1_EFF_HIST = HISTAMINE * (1.0 - EFFECT_CETI);      // Histamine x (1-H1RO)
  double LT_EFF      = CYS_LT   * (1.0 - EFFECT_MLKT);       // LT x (1-CysLT1 block)

  // ==============================================================
  // Symptom scores (0-3 scale, Emax model)
  // ==============================================================
  // Sneezing: mainly histamine + sensory nerve
  double SNEEZE = SNEEZE_MAX * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST);

  // Rhinorrhea: histamine + LTs + glandular
  double RHINORRHEA = RHINO_MAX * (0.6 * H1_EFF_HIST + 0.4 * LT_EFF) /
                      (HIST_EC50 + 0.6 * H1_EFF_HIST + 0.4 * LT_EFF);

  // Congestion: mainly LTs + PGD2 (approximated by CYS_LT) + eosinophil
  double EOS_CONG = EOS_N / (EOS_EC50 + EOS_N);
  double CONGESTION = CONG_MAX * (0.5 * LT_EFF / (LT_EC50 + LT_EFF) + 0.3 * EOS_CONG + 0.2);

  // Pruritus: histamine + Th2 cytokines (IL-31 proxy -> IL-13)
  double IL13_NORM = IL13 / (IL13 + 8.0);
  double PRURITUS = PRUR_MAX * (0.7 * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST) + 0.3 * IL13_NORM);

  // TNSS (0-12)
  double TNSS = SNEEZE + RHINORRHEA + CONGESTION + PRURITUS;

  // ==============================================================
  // Biomarker endpoints
  // ==============================================================
  double FREE_IGE_IU  = IGE_FREE;          // IU/mL
  double BLOOD_EOS_UL = EOS_B;             // cells/uL
  double NASAL_EOS_UL = EOS_N;             // cells/uL tissue
  double TRYPTASE_MCG = MAST_ACT * 15.0;  // serum tryptase proxy (ng/mL)

$CAPTURE @annotated
  CETI_CP     : Cetirizine plasma concentration, unchanged-name reporting alias of C_CETI (ng/mL)
  MLKT_CP      : Montelukast plasma concentration, unchanged-name reporting alias of C_MLKT (ng/mL)
  OMA_CP       : Omalizumab plasma concentration, unchanged-name reporting alias (ug/mL)
  FP_LOCAL_NM  : Fluticasone local nasal concentration, unchanged-name reporting alias of C_FP (nM)
  H1_RO        : H1 receptor occupancy, unchanged-name reporting alias of EFFECT_CETI*100 (%)
  GR_OCC_FP     : GR occupancy, unchanged-name reporting alias of EFFECT_FP*100 (%)
  CYSLTR1_INH  : CysLT1 receptor inhibition, unchanged-name reporting alias of EFFECT_MLKT*100 (%)
  SNEEZE       : Sneezing score (0-3)
  RHINORRHEA   : Rhinorrhea score (0-3)
  CONGESTION   : Congestion score (0-3)
  PRURITUS     : Pruritus score (0-3)
  TNSS         : Total Nasal Symptom Score (0-12)
  FREE_IGE_IU  : Free IgE (IU/mL)
  BLOOD_EOS_UL : Blood eosinophil (cells/uL)
  NASAL_EOS_UL : Nasal tissue eosinophil (cells/uL tissue)
  TRYPTASE_MCG : Serum tryptase proxy (ng/mL)
  H1_EFF_HIST  : Histamine effective after H1 blockade (normalized AU)
  LT_EFF       : CysLT effective after CysLT1 blockade (AU)
  IGE_REDUCTION_PCT : Free IgE reduction vs baseline (%)
  C_CETI       : Cetirizine exposed concentration, canonical single site (ng/mL)
  C_FP         : Fluticasone exposed concentration, canonical single site (nM)
  C_MLKT       : Montelukast exposed concentration, canonical single site (ng/mL)
  C_OMA        : Omalizumab exposed concentration, canonical single site (mg/mL)
  EFFECT_CETI  : Cetirizine Hill effect, H1R occupancy fraction (0-1)
  EFFECT_FP    : Fluticasone Hill effect, GR occupancy fraction (0-1)
  EFFECT_MLKT  : Montelukast Hill effect, CysLT1 blockade fraction (0-1)
  EFFECT_OMA   : Omalizumab bespoke effect, fractional bound-IgE occupancy (0-1)
  // COMPLEX_OMA is a $CMT compartment -- exposed as an output column
  // automatically; not re-listed here (mrgsolve 2.0.1 rejects a
  // compartment name appearing in $CAPTURE, same defect class as the
  // original’s bare-compartment capture list -- see refactor notes).

'

# ============================================================
# Compile model
# ============================================================
AR_model_refactored <- mcode("AllergicRhinitisQSPRefactored", AR_model_refactored_code)

# ============================================================
# Helper: build dosing regimens (unchanged from original, compartment
# names renamed to the new convention: CETI_D->GUT_CETI, FP_LOC->CENT_FP,
# MLKT_D->GUT_MLKT, OMA_D->GUT_OMA)
# ============================================================
build_doses <- function(
    cetirizine  = FALSE,  # 10 mg QD oral
    fluticasone = FALSE,  # 200 ug/day intranasal (split as 100 ug BID)
    montelukast = FALSE,  # 10 mg QD oral
    omalizumab  = FALSE,  # 300 mg q4w SC
    allergen_pulse = FALSE, # allergen challenge on day 28
    sim_duration_d = 84
) {
  ev_list <- list()

  if (cetirizine)  ev_list[["ceti"]]  <- ev(cmt = "GUT_CETI", amt = 10,   time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (fluticasone) ev_list[["fp"]]    <- ev(cmt = "CENT_FP",  amt = 450,  time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (montelukast) ev_list[["mlkt"]]  <- ev(cmt = "GUT_MLKT", amt = 10,   time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (omalizumab)  ev_list[["oma"]]   <- ev(cmt = "GUT_OMA",  amt = 300,  time = 0, ii = 28*24, addl = 2)

  if (allergen_pulse) {
    ev_list[["ag"]] <- ev(cmt = "AG", amt = 5.0, time = 28*24)
  }

  if (length(ev_list) == 0) return(ev(cmt = "AG", amt = 5.0, time = 28*24))
  Reduce(c, ev_list)
}

# ============================================================
# Scenario definitions (unchanged from original)
# ============================================================
scenarios <- list(
  list(
    name        = "1. Natural History\n(Allergen only, no Tx)",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#E53935"
  ),
  list(
    name        = "2. Cetirizine 10 mg QD",
    cetirizine  = TRUE,  fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#1E88E5"
  ),
  list(
    name        = "3. Fluticasone FP 200 ug/d",
    cetirizine  = FALSE, fluticasone = TRUE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#43A047"
  ),
  list(
    name        = "4. Montelukast 10 mg QD",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = TRUE,  omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#FB8C00"
  ),
  list(
    name        = "5. Cetirizine + Fluticasone\n(Combination)",
    cetirizine  = TRUE,  fluticasone = TRUE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#8E24AA"
  ),
  list(
    name        = "6. Omalizumab 300 mg q4w\n(Anti-IgE)",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = TRUE,
    allergen_pulse = TRUE, color = "#00897B"
  ),
  list(
    name        = "7. Triple Therapy\n(Ceti + FP + MLKT)",
    cetirizine  = TRUE,  fluticasone = TRUE,
    montelukast = TRUE,  omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#6D4C41"
  )
)

# ============================================================
# Simulation function
# ============================================================
run_scenario <- function(sc) {
  evs <- build_doses(
    cetirizine  = sc$cetirizine,
    fluticasone = sc$fluticasone,
    montelukast = sc$montelukast,
    omalizumab  = sc$omalizumab,
    allergen_pulse = sc$allergen_pulse,
    sim_duration_d = 84
  )
  AR_model_refactored %>%
    ev(evs) %>%
    mrgsim(end = 84 * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(scenario = sc$name, color = sc$color, time_d = time / 24)
}

# Run all scenarios
message("Running ", length(scenarios), " scenarios...")
results <- bind_rows(lapply(scenarios, run_scenario))

# ============================================================
# Visualization functions (unchanged from original)
# ============================================================

# 1. TNSS over time
plot_tnss <- function(data) {
  ggplot(data, aes(x = time_d, y = TNSS, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = setNames(
      unique(data$color), unique(data$scenario)
    )) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.7) +
    annotate("text", x = 28.5, y = 11, label = "Allergen\nChallenge", size = 3, hjust = 0) +
    labs(
      title = "Total Nasal Symptom Score (TNSS) - All Scenarios",
      subtitle = "0-12 scale; allergen challenge on Day 28",
      x = "Time (days)", y = "TNSS (0-12)",
      color = "Scenario"
    ) +
    ylim(0, 12) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", legend.text = element_text(size = 8))
}

# 2. Individual symptom scores
plot_symptoms <- function(data) {
  sym_df <- data %>%
    select(time_d, scenario, color, SNEEZE, RHINORRHEA, CONGESTION, PRURITUS) %>%
    pivot_longer(c(SNEEZE, RHINORRHEA, CONGESTION, PRURITUS), names_to = "symptom", values_to = "score")

  ggplot(sym_df, aes(x = time_d, y = score, color = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~symptom, scales = "free_y", ncol = 2) +
    scale_color_manual(values = setNames(
      unique(data$color), unique(data$scenario)
    )) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.5) +
    labs(
      title = "Individual Symptom Scores",
      x = "Time (days)", y = "Score (0-3)",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

# 3. PK - H1R occupancy (cetirizine)
plot_pk_ceti <- function(data) {
  df <- data %>% filter(grepl("Cetirizine|Combination|Triple", scenario))
  ggplot(df, aes(x = time_d, y = H1_RO, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(
      title = "H1 Receptor Occupancy - Cetirizine",
      subtitle = "Steady-state target >=80% for clinical efficacy",
      x = "Time (days)", y = "H1 Receptor Occupancy (%)",
      color = "Scenario"
    ) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "blue") +
    annotate("text", x = 5, y = 82, label = "80% target", color = "blue", size = 3) +
    ylim(0, 100) +
    theme_bw(base_size = 12)
}

# 4. Biomarkers - Free IgE & Omalizumab
plot_ige <- function(data) {
  df_oma <- data %>% filter(grepl("Omalizumab", scenario) | grepl("Natural", scenario))
  p1 <- ggplot(df_oma, aes(x = time_d, y = FREE_IGE_IU, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(title = "Free IgE (IU/mL)", x = "Time (days)", y = "Free IgE (IU/mL)") +
    theme_bw(base_size = 11)
  print(p1)
}

# 5. Eosinophil dynamics
plot_eos <- function(data) {
  df <- data %>%
    select(time_d, scenario, color, BLOOD_EOS_UL, NASAL_EOS_UL) %>%
    pivot_longer(c(BLOOD_EOS_UL, NASAL_EOS_UL), names_to = "compartment", values_to = "eos")
  ggplot(df, aes(x = time_d, y = eos, color = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~compartment, scales = "free_y") +
    scale_color_manual(values = setNames(unique(data$color), unique(data$scenario))) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.5) +
    labs(
      title = "Eosinophil Dynamics",
      x = "Time (days)", y = "Eosinophils (cells/uL)",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

# 6. Summary table: mean TNSS at peak (day 28-35) & end (day 77-84)
summary_table <- function(data) {
  data %>%
    group_by(scenario) %>%
    summarise(
      TNSS_baseline    = mean(TNSS[time_d < 28], na.rm = TRUE),
      TNSS_peak        = max(TNSS[time_d >= 28 & time_d <= 35], na.rm = TRUE),
      TNSS_wk12        = mean(TNSS[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Free_IgE_wk12    = mean(FREE_IGE_IU[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Blood_Eos_wk12   = mean(BLOOD_EOS_UL[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Nasal_Eos_wk12   = mean(NASAL_EOS_UL[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      TNSS_pct_chg = round((TNSS_wk12 - TNSS_peak[1]) / TNSS_peak[1] * 100, 1)
    )
}

# ============================================================
# Generate plots
# ============================================================
p_tnss <- plot_tnss(results)
p_sym  <- plot_symptoms(results)
p_eos  <- plot_eos(results)

print(p_tnss)
print(p_sym)
print(p_eos)

tbl <- summary_table(results)
print(tbl)

message("\nScenario summary:")
message("  TNSS at Week 12 (days 77-84):")
for (i in seq_len(nrow(tbl))) {
  message(sprintf("    %-45s  %.2f", tbl$scenario[i], tbl$TNSS_wk12[i]))
}
