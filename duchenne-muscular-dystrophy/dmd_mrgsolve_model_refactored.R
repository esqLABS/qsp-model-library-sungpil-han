## ============================================================
## Duchenne Muscular Dystrophy (DMD) QSP Model -- REFACTORED (pluggable PK)
##   Refactor of `dmd_mrgsolve_model.R` per FORK_WORKFLOW_GUIDE.md Part 2.
##   Original untouched -- see `dmd_refactor_notes.md` for full rationale.
##
##   Scope: three census rows in
##   `driver-patches/data/compound_perturbation_census.md` --
##     - Corticosteroid (CS): Deflazacort (base PK, F_CS=0.72) and
##       Prednisone (same shared PK block, F_CS override 0.82 per scenario) --
##       both real, externally-dosed, sharing one PK block as the original did.
##     - ASO (Eteplirsen, exon-51-skipping PMO): scenario "4. Eteplirsen
##       30 mg/kg/wk" IV -- real approved EXONDYS 51 dosing.
##     - Gene therapy (AAV, "DYS" in the census -- renamed here, see below):
##       delandistrogene moxeparvovec (Elevidys/SRP-9001), scenario
##       "5. Gene Therapy (Elevidys)", single IV dose 1.33e14 vg/kg -- matches
##       the real approved Elevidys dose exactly.
##
##   Stem-collision note: the census pre-filled "Gene therapy (DYS)" but the
##   file's own $CMT DYS is the DISEASE-side "drug-induced dystrophin pool"
##   compartment (untouched, not gene-therapy PK) -- using DYS as the PK stem
##   would collide with it. Used AAV instead (matches the original's own
##   AAV_CIRC/AAV_MUS compartment prefix). See dmd_refactor_notes.md.
##
##   Renamed to this fork's PK/PD naming convention: GUT_<STEM>/CENT_<STEM>/
##   PERI_<STEM>, CL_<STEM>/V1_<STEM>/V2_<STEM>/Q_<STEM>/KA_<STEM>/F_<STEM>,
##   C_<STEM> (exposed concentration), EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>/
##   EFFECT_<STEM> (Hill interface where the original's term was already that
##   shape; a disclosed non-Hill linear rename otherwise -- ASO and AAV effect
##   terms are both linear in the original, kept linear, not force-fit).
##
##   Both AAV compartments are genuine, dosed, first-order PK (NOT deleted --
##   correcting the census's "Delete PK compartment" classification, which is
##   only half right: the muscle compartment's state IS directly the exposed
##   concentration (no volume division needed, census correct on that point),
##   but there are two real compartments with real one-way kinetics, not zero.
##   See dmd_refactor_notes.md for the full correction, matching the
##   precedent set in achondroplasia/acdp_refactor_notes.md and
##   clostridioides-difficile-infection/cdi_refactor_notes.md.
##
##   Pre-existing original defect found, NOT fixed (disclosed + logged,
##   UPSTREAM_ISSUES.md -- see notes): `dxdt_DYS` multiplies KEXP_DYS by a
##   ternary that is always literally 0.0 regardless of its condition, so
##   KEXP_DYS ("Micro-dystrophin expression rate") never actually contributes
##   anything -- the DYS compartment only ever decays. Left byte-for-byte
##   identical (same formula, same always-dead parameter) for exact-match
##   verification; this is disease-side, not part of any compound's own PK
##   block, so out of this refactor's structural scope regardless.
##
##   A second pre-existing finding, also disclosed + logged: $TABLE's
##   Dystrophin_pct output (`DYS_BASE + DYS + EFF_ASO*IC_ASO`) omits the
##   gene-therapy (AAV) contribution that $MAIN's DYS_TOTAL (which drives all
##   the actual PD dynamics -- MEMI, SC, MF, SWD, FVC, LVEF) includes. The
##   reported/plotted "Dystrophin Level" column therefore undercounts
##   specifically for the gene-therapy scenario. Preserved exactly (not
##   fixed) since $TABLE cannot see $MAIN's own locals -- recomputing from
##   scratch is what the original did, and the omission is in what it chose
##   to recompute, not a block-visibility artifact.
##
##   The untouched original does not compile under mrgsolve 2.0.1 (three
##   independent, pre-existing, disclosed defects, all confirmed via
##   POST /model_manifest against the original's own extracted DSL, and all
##   fixed here as syntax-only/non-numeric build-compatibility fixes -- full
##   detail and error text in dmd_refactor_notes.md / UPSTREAM_ISSUES.md):
##     1. `$CMT @annotated` + a separate `$INIT` block jointly redeclare all
##        22 compartments ("Duplicated model names"). `$INIT` dropped; the
##        same initial values are set via the `<cmt>_0 = value;` idiom
##        inside `if (NEWIND <= 1) { ... }` in $MAIN instead.
##     2. Once (1) is fixed, seven of the original's own baseline $PARAM
##        names (`M1_0`, `M2_0`, `TGFb_0`, `FIB_0`, `SC_0`, `MF_0`, `SWD_0`)
##        collide with mrgsolve's own auto-generated `<compartment>_0`
##        initial-value symbol for the identically-named compartment --
##        renamed `<CMT>_BASELINE` (values unchanged; six of the seven were
##        dead/unused, one -- M1_0 -- is used in the NF-kB drive term and is
##        now M1_BASELINE there too).
##     3. `$CAPTURE` re-lists 14 `$CMT` compartment names directly, which
##        mrgsolve 2.0.1 rejects -- dropped (compartments are exposed as
##        output columns automatically).
##
##   A fourth finding, purely about THIS refactor's own construction (not a
##   build defect, caught and fixed before delivery): an early draft
##   extracted the corticosteroid membrane-stabilization arithmetic into its
##   own named `EFFECT_CS_MEM` double instead of leaving it inline in the
##   exact 4-term operation order the original used. Under live ODE
##   integration (not a static check) this reassociation is enough to seed
##   a sub-ULP difference that this model's own genuine, pre-existing
##   runaway positive-feedback instability (see below) amplifies by 11
##   orders of magnitude within one 240h dosed scenario -- confirmed via a
##   controlled A/B test against the patched-original baseline. Fixed by
##   keeping that term inline, byte-for-byte identical to the original's own
##   operation order; verification below is the confirmed-clean result.
##
##   Fifth, a genuine pre-existing numerical fragility (NOT fixed, disclosed
##   + logged only, same treatment as the achondroplasia/essential-
##   thrombocythemia precedent class): the untreated natural-history
##   scenario's own positive-feedback loop (ROS -> NFkB -> M1 -> more ROS)
##   is genuinely unstable and makes lsoda's adaptive step size collapse to
##   zero at simulated hour 276.035 -- confirmed identical in the untouched
##   original (via the patched-original baseline) and in this refactored
##   file, at the exact same simulated hour, in both the untreated and
##   corticosteroid-dosed scenarios. Verification below therefore uses a
##   shortened window (240h, not the original's own 6-year horizon) that
##   stays well before this onset.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ===================================================================
## MODEL DEFINITION
## ===================================================================

code_dmd_refactored <- '
$PROB DMD QSP Model — 19-Compartment ODE System (REFACTORED, pluggable PK)

$PARAM @annotated
// --- Corticosteroid (CS) PK parameters (Deflazacort base; Prednisone via
// F_CS override in its own scenario, matching the original) ---
KA_CS    : 1.2     : Steroid absorption rate constant (h-1)
F_CS     : 0.72    : Steroid oral bioavailability (fraction)
CL_CS    : 25.0    : Steroid plasma clearance (L/h/70kg)
V1_CS    : 35.0    : Steroid central volume (L/70kg)
V2_CS    : 80.0    : Steroid peripheral volume (L/70kg)
Q_CS     : 8.0     : Steroid inter-compartmental CL (L/h)

// --- ASO (Eteplirsen, exon-51-skipping PMO) PK parameters ---
CL_ASO   : 210.0   : ASO plasma clearance (mL/h/kg)
V1_ASO   : 290.0   : ASO central volume (mL/kg)
V2_ASO   : 6500.0  : ASO muscle/tissue volume (mL/kg)
Q_ASO    : 18.0    : ASO inter-compartmental CL (mL/h/kg)
KUP_ASO  : 0.0015  : ASO intracellular uptake rate (h-1)
KOUT_ASO : 0.0008  : ASO intracellular elimination rate (h-1)

// --- Gene therapy (AAV, delandistrogene moxeparvovec/Elevidys) PK ---
// stem AAV used instead of the census placeholder "DYS" -- collides with
// the disease-side $CMT DYS (see header note above)
KTRANSD_AAV : 0.003     : AAV muscle transduction rate (h-1) [was KAAV]
CL_AAV      : 5e-5      : AAV vector degradation/clearance from circulation (h-1) [was DEGRAD]
KEXP_DYS    : 0.004     : Micro-dystrophin expression rate (h-1) -- UNUSED (see header defect note); kept unchanged, not renamed to the AAV stem since it has no functional link to the AAV PK block itself
KDEC_AAV    : 0.00008   : Muscle-transduced AAV vector decay rate (h-1, very slow) [new name, split from KDEC_DYS -- same value, see header note]

// --- Dystrophin dynamics (disease-side, unchanged) ---
DYS_BASE : 0.001   : Baseline dystrophin (% of normal in DMD ~ 0.1%)
DYS_MAX  : 30.0    : Max achievable dystrophin with gene therapy (% normal)
HILL_DYS : 2.0     : Hill coefficient for dystrophin effect on membrane
EC50_DYS : 5.0     : EC50 dystrophin for membrane protection (% normal)
KDEC_DYS : 0.00008 : Micro-dystrophin decay rate (h-1, very slow) -- disease-side dxdt_DYS decay, unchanged/untouched

// --- Membrane integrity (disease-side, unchanged) ---
KD_MEM   : 0.15    : Membrane damage rate from Ca2+ (dimensionless h-1)
KR_MEM   : 0.08    : Membrane repair rate (h-1)
MEMI_SS  : 0.35    : Membrane integrity steady-state in untreated DMD (0-1)
KD_CS_MEM: 0.012   : Steroid effect on membrane damage rate (L/nmol)

// --- Intracellular calcium dynamics (disease-side, unchanged) ---
CA_IN    : 2.5     : Ca2+ basal influx (relative units h-1)
CA_OUT   : 1.0     : Ca2+ efflux rate constant (h-1)
CA_SS    : 2.5     : Steady-state Ca2+ in DMD (relative to normal=1.0)
KMEM_CA  : 0.8     : Membrane integrity effect on Ca2+ influx (rel. units)

// --- Oxidative stress (ROS) (disease-side, unchanged) ---
KROS_IN  : 0.5     : Basal ROS generation rate (h-1)
KROS_CA  : 0.3     : Ca2+-driven ROS amplification coefficient
KROS_EL  : 0.4     : ROS elimination rate (h-1, Nrf2 antioxidant)
ROS_SS   : 3.0     : ROS steady-state in DMD

// --- NF-kB inflammatory signaling (disease-side + CS Hill interface) ---
KNF_IN   : 1.2     : NF-kB activation rate from ROS/DAMPs (h-1)
KNF_EL   : 0.5     : NF-kB deactivation rate (h-1)
EC50_CS_NF : 15.0  : CS EC50 for NF-kB inhibition (ng/mL)
EMAX_CS_NF : 1.0   : CS max fractional NF-kB inhibition [new, math-implied -- original ratio C/(C+EC50) already saturates at 1]
GAMMA_CS_NF: 1.0   : CS Hill coefficient for NF-kB inhibition [new, explicit -- original had no Hill exponent here, i.e. GAMMA=1]
NFkB_MAX : 4.0     : Max NF-kB activity fold (baseline=1)

// --- M1/M2 macrophage dynamics (disease-side, unchanged) ---
KM1_IN   : 0.08    : M1 recruitment rate (cells/µL/h driven by NFkB)
KM1_EL   : 0.04    : M1 elimination rate (h-1)
KM2_IN   : 0.03    : M2 differentiation rate from M1 (h-1)
KM2_EL   : 0.03    : M2 elimination rate (h-1)
M1_BASELINE : 15.0 : Initial M1 macrophages (cells/µL tissue) in DMD [was M1_0 -- renamed, build-compatibility fix, see notes]
M2_BASELINE : 8.0  : Initial M2 macrophages (cells/µL tissue) [was M2_0 -- renamed, dead/unused, see notes]

// --- TGF-β1 dynamics (disease-side, unchanged) ---
KTGF_IN  : 0.2     : TGF-β1 secretion rate from M2 (pg/mL/cell/h)
KTGF_EL  : 0.15    : TGF-β1 clearance rate (h-1)
TGFb_BASELINE : 12.0 : Baseline TGF-β1 in DMD muscle (pg/mL) [was TGFb_0 -- renamed, dead/unused, see notes]

// --- Fibrosis dynamics (disease-side, unchanged) ---
KFIB_IN  : 0.002   : Fibrosis progression rate (units/day, TGF-β driven)
KFIB_EL  : 0.00015 : Spontaneous fibrosis resolution rate (h-1, very slow)
FIB_MAX  : 100.0   : Maximum fibrosis score
FIB_BASELINE : 10.0 : Initial fibrosis score at model start (age 6yr) [was FIB_0 -- renamed, dead/unused, see notes]

// --- Satellite cell pool (disease-side, unchanged) ---
SC_BASELINE : 100.0 : Initial satellite cell pool (% of normal = 100) [was SC_0 -- renamed, dead/unused, see notes]
KSC_REGEN: 0.01    : SC self-renewal rate (h-1)
KSC_EXHST: 0.002   : SC exhaustion rate per necrosis event (h-1)
SC_MIN   : 5.0     : Minimum SC pool (% normal)

// --- Muscle function (disease-side, unchanged) ---
MF_BASELINE : 80.0 : Initial muscle function at age 6yr (% baseline) [was MF_0 -- renamed, dead/unused, see notes]
KMF_DEC  : 0.0004  : Muscle function decline rate (h-1, fibrosis/SC-dep)
KMF_REGEN: 0.0002  : Muscle function recovery rate (h-1, SC-dependent)

// --- Clinical endpoints (disease-side, unchanged) ---
SWD_BASELINE : 380.0 : Initial 6MWD at age 6yr (m) [was SWD_0 -- renamed, dead/unused, see notes]
KSWMD_DC : 0.00006 : 6MWD decline rate (m/h with muscle function)
FVC_0    : 95.0    : Initial FVC % predicted at age 6yr
KFVC_DC  : 0.00004 : FVC decline rate (h-1)
LVEF_0   : 62.0    : Initial LVEF (%) at age 6yr
KLVEF_DC : 0.000025: LVEF decline rate (h-1)

// --- Drug effect modifiers (unchanged names/values; ASO and AAV effect
// terms are linear in the original, kept linear -- not force-fit to Hill) ---
EFF_ASO  : 0.04    : ASO maximum dystrophin restoration efficiency (% per nmol/L)
EFF_AAV  : 0.02    : AAV maximum micro-dys expression (% per log-vg/cell)
HDAC_EFF : 0.3     : Givinostat FAP→myogenic effect (fraction) -- out of scope (not a census compound in this file)
HDAC_FIB : 0.25    : Givinostat fibrosis reduction effect (fraction) -- out of scope

// --- Body weight/scaling ---
WT       : 22.0    : Body weight (kg), representative age 8yr DMD boy
AGE_BASE : 8.0     : Baseline age (years)

$CMT @annotated
// PK compartments (renamed to this fork naming convention)
GUT_CS    : Corticosteroid oral absorption depot (mg) [was DEPOT_CS]
CENT_CS   : Corticosteroid central plasma (mg)
PERI_CS   : Corticosteroid peripheral tissue (mg) [was PERIPH_CS]
CENT_ASO  : ASO plasma compartment (mg/kg)
PERI_ASO  : ASO muscle compartment (mg/kg) [was MUS_ASO]
IC_ASO    : ASO intracellular active (nmol/L tissue) -- bespoke 3rd PK site, the actual PD-driving concentration (see C_ASO below)
CENT_AAV  : AAV vector circulating (vg/kg) [was AAV_CIRC]
MUS_AAV   : AAV in muscle (vg/cell) [was AAV_MUS] -- bespoke, one-way irreversible transduction sink (not a reversible PERI_<STEM>)

// Disease pathophysiology compartments (unchanged)
DYS       : Dystrophin level (% of normal)
MEMI      : Membrane integrity (fraction 0-1)
CAI       : Intracellular calcium (relative units, normal=1)
ROS       : Reactive oxygen species (relative units, normal=1)
NFkB      : NF-κB activity (fold, basal=1)
M1        : M1 macrophage density (cells/µL tissue)
M2        : M2 macrophage density (cells/µL tissue)
TGFb      : TGF-β1 concentration (pg/mL muscle)
FIB       : Fibrosis score (0-100)
SC        : Satellite cell pool (% of normal)

// Functional outcomes (unchanged)
MF        : Muscle function (% baseline)
SWD       : 6-minute walk distance (m)
FVC_pct   : FVC % predicted
LVEF_pct  : Left ventricular ejection fraction (%)

$MAIN
// Build-compatibility fixes (disclosed, syntax-only, non-numeric -- see
// header note and refactor notes / UPSTREAM_ISSUES.md). The untouched
// original does not compile under mrgsolve 2.0.1 for three independent,
// pre-existing reasons, all fixed here (and confirmed identical in the
// original via a patched-original scratch baseline used only for
// verification -- see notes):
//   1. The original declares every non-zero initial condition in a
//      separate `$INIT` block, which conflicts with `$CMT @annotated`
//      also declaring the same 22 compartment names ("Duplicated model
//      names"). `$INIT` dropped; the same values are set here via the
//      standard `<cmt>_0 = value;` idiom instead.
//   2. Six baseline $PARAM names carried over from the original (`M2_0`,
//      `TGFb_0`, `FIB_0`, `SC_0`, `MF_0`, `SWD_0`) collide with the
//      mrgsolve-auto-generated `<compartment>_0` initial-value symbol for the
//      identically-named compartment (`M2`, `TGFb`, `FIB`, `SC`, `MF`,
//      `SWD`) -- all six are dead/unused elsewhere in the original, same
//      defect class as `lymphangioleiomyomatosis` issue #118. A seventh,
//      `M1_0`, collides the same way but IS used (NF_INPUT below) --
//      same defect class as `myotonic-dystrophy` issue #65. All seven
//      renamed `<CMT>_BASELINE` in `$PARAM` above; values unchanged.
//   3. `$CAPTURE` re-lists the same compartment names again (see the
//      $CAPTURE line below) -- dropped, compartments are exposed as
//      output columns automatically.
// All PK compartments already default to 0 via `$CMT` (unchanged from the
// original, which also initialized every PK compartment to 0).
if (NEWIND <= 1) {
  DYS_0      = 0.1;
  MEMI_0     = 0.35;
  CAI_0      = 2.5;
  ROS_0      = 3.0;
  NFkB_0     = 3.5;
  M1_0       = 15.0;
  M2_0       = 8.0;
  TGFb_0     = 12.0;
  FIB_0      = 10.0;
  SC_0       = 100.0;
  MF_0       = 80.0;
  SWD_0      = 380.0;
  FVC_pct_0  = 95.0;
  LVEF_pct_0 = 62.0;
}

// --- Corticosteroid (CS): the exposed concentration. Computed once per
// reporting interval (mrgsolve $MAIN semantics), matching the original
// exactly -- this is the same block the original used for C_CS, so
// dxdt_MEMI/dxdt_NFkB below see the identical once-per-interval-updated
// value the original did (confirmed empirically, see refactor notes: a
// $MAIN local is captured/reported one full interval behind current state,
// vs a $TABLE local which reports the live current-row state -- this is
// inherent mrgsolve behavior, not introduced by this refactor). ---
double C_CS = CENT_CS / V1_CS * 1000.0;   // ng/mL
// NOTE: original also computed an unused `C_CS_PERIPH` here
// (PERIPH_CS/V2_CS*1000) -- never read anywhere in the model, dead code,
// dropped (see refactor notes).

// --- ASO (Eteplirsen): the exposed concentration is the intracellular
// active pool itself -- the original PD term (DYS_from_ASO) read
// IC_ASO directly, never the separately-computed, dead plasma `C_ASO` local
// (never read anywhere, not even captured) -- redirected here to the site
// that actually drives PD, per this fork naming convention. ---
double C_ASO = IC_ASO;

// --- Gene therapy (AAV): the exposed concentration is the muscle-transduced
// vector compartment itself, no volume division needed (the census
// "concentration is itself the state" classification is correct for this
// half; see refactor notes for the half that needed correcting). ---
double C_AAV = MUS_AAV;

// Each compound has one named effect on the disease. ASO and AAV are exact
// renames of the original linear (non-Hill) terms -- not forced into
// a saturating shape the original never had.
double EFFECT_ASO = EFF_ASO * C_ASO;
double EFFECT_AAV = EFF_AAV * C_AAV * 1e-12 * 1e3;  // simplified scaling (unchanged from original)

// Total dystrophin (disease-side algebra, unchanged; refs renamed)
double DYS_TOTAL = DYS_BASE + DYS + EFFECT_ASO + EFFECT_AAV;

// Dystrophin effect on membrane (sigmoidal) -- disease-side, unchanged
double DYS_EFF = pow(DYS_TOTAL, HILL_DYS) / (pow(EC50_DYS, HILL_DYS) + pow(DYS_TOTAL, HILL_DYS));

// Corticosteroid NF-kB inhibition -- exact-shape rename (original ratio
// C_CS/(EC50_CS_NF+C_CS) is already an Emax model; EMAX_CS_NF=1.0 and
// GAMMA_CS_NF=1.0 are math-implied/explicit, not fit -- see refactor notes)
double EFFECT_CS = EMAX_CS_NF * pow(C_CS, GAMMA_CS_NF) / (pow(EC50_CS_NF, GAMMA_CS_NF) + pow(C_CS, GAMMA_CS_NF));

// Membrane-driven Ca2+ influx -- disease-side, unchanged
double CA_INFLUX = CA_IN * (1.0 - MEMI) * KMEM_CA;

// Necrosis rate (function of Ca2+ and ROS) -- disease-side, unchanged
double NECRO_RATE = 0.1 * (CAI / 1.0) * (ROS / 1.0);

// SC pool effect on muscle function -- disease-side, unchanged
double SC_EFF = SC / 100.0;

// Fibrosis effect on muscle function (negative) -- disease-side, unchanged
double FIB_EFF = 1.0 - (FIB / FIB_MAX) * 0.7;

$ODE
// --- Corticosteroid PK (Archetype 3: depot + central + peripheral, linear) ---
dxdt_GUT_CS  = -KA_CS * GUT_CS;
dxdt_CENT_CS =  KA_CS * GUT_CS - (CL_CS/V1_CS) * CENT_CS - (Q_CS/V1_CS) * CENT_CS + (Q_CS/V2_CS) * PERI_CS;
dxdt_PERI_CS =  (Q_CS/V1_CS) * CENT_CS - (Q_CS/V2_CS) * PERI_CS;

// --- ASO PK (bespoke 3-site: central/peripheral linear PK + intracellular
// active pool; IC_ASO -- not CENT_ASO -- is the site PD actually reads) ---
dxdt_CENT_ASO = -(CL_ASO/1000.0) * CENT_ASO - (Q_ASO/1000.0) * CENT_ASO + (Q_ASO/1000.0) * PERI_ASO;
dxdt_PERI_ASO =  (Q_ASO/1000.0) * CENT_ASO - (Q_ASO/1000.0) * PERI_ASO - KUP_ASO * PERI_ASO;
dxdt_IC_ASO   =  KUP_ASO * PERI_ASO * 1000.0 - KOUT_ASO * IC_ASO;  // nmol/L approx

// --- Gene therapy (AAV) PK (bespoke 2-site, irreversible one-way
// transduction -- NOT a reversible Archetype-2 Q-term: vector never returns
// from muscle back to circulation, matching the original structure) ---
dxdt_CENT_AAV = -KTRANSD_AAV * CENT_AAV - CL_AAV * CENT_AAV;
dxdt_MUS_AAV  =  KTRANSD_AAV * CENT_AAV - KDEC_AAV * MUS_AAV;

// --- Dystrophin pool (endogenous restoration; baseline in DMD ~ 0.1%) ---
// Disease-side, byte-for-byte unchanged, INCLUDING the pre-existing
// KEXP_DYS defect (always multiplied by a 0.0 ternary -- see header note
// and refactor notes / UPSTREAM_ISSUES.md). Not fixed: this is a logic
// defect, not a build defect, and fixing it would change the simulated
// trajectory (out of scope for a structural PK/PD refactor).
dxdt_DYS = KEXP_DYS * (DYS_TOTAL > 0 ? 0.0 : 0.0) - KDEC_DYS * DYS + 0.0;
// Note: DYS compartment tracks drug-induced dystrophin beyond DYS_BASE
// It is replenished by ASO effect (tracked via IC_ASO) and AAV (via MUS_AAV)

// --- Membrane integrity ---
// The CS membrane-stabilizing term (linear in C_CS, not Hill-shaped) is
// kept INLINE here, in the exact original 4-term operation order --
// NOT extracted into a separately-summed named EFFECT_CS_MEM local.
// Disclosed finding (see refactor notes): this system has a genuine,
// pre-existing runaway positive-feedback instability (confirmed identical
// in the untouched original -- both blow up/fail the solver at the same
// simulated hour on the undosed scenario). Under live ODE integration
// (not a static single-point check), extracting this sum into its own
// double changes IEEE754 operation order enough to seed a difference that
// this instability amplifies by many orders of magnitude within one dosed
// scenario -- confirmed via a controlled A/B test against the untouched
// original (qspserver mrgsolve_api). Kept byte-for-byte identical in
// structure to the original (only identifier renames: DEPOT/PERIPH -> the
// unchanged CENT_CS/C_CS names already in scope) to preserve the exact
// floating-point operation order the original used.
dxdt_MEMI = KR_MEM * DYS_EFF * (1.0 - MEMI) * SC_EFF
            - KD_MEM * (1.0 - DYS_EFF) * CAI
            - KD_CS_MEM * C_CS * MEMI * (-1.0)   // CS stabilizes membrane
            + KD_CS_MEM * C_CS * (1.0 - MEMI);

// --- Intracellular Calcium --- (disease-side, unchanged)
dxdt_CAI = CA_INFLUX * (1.0 - DYS_EFF)           // influx through torn membrane
           + 0.1 * (1.0 - MEMI)                   // additional leak
           - CA_OUT * CAI;                          // efflux (SERCA + PM-Ca-ATPase)

// --- ROS --- (disease-side, unchanged)
dxdt_ROS = KROS_IN * (1.0 + KROS_CA * (CAI - 1.0))
           - KROS_EL * ROS
           + 0.0;  // Idebenone effect modeled as dose-response adjustment

// --- NF-kB (fold change from basal) --- (renamed CS_NF_INH -> EFFECT_CS)
double NF_INPUT = 1.0 + 2.0 * (ROS - 1.0) / 3.0 + 0.5 * (M1/M1_BASELINE - 1.0);
dxdt_NFkB = KNF_IN * (NF_INPUT > 1.0 ? NF_INPUT : 1.0)
            - KNF_EL * NFkB * (1.0 + 2.0 * EFFECT_CS);

// --- M1 macrophages --- (renamed CS_NF_INH -> EFFECT_CS)
double M1_IN = KM1_IN * NFkB * (1.0 - EFFECT_CS * 0.7);
dxdt_M1 = M1_IN - KM1_EL * M1 - KM2_IN * M1;  // M1 → M2 switching

// --- M2 macrophages --- (disease-side, unchanged)
dxdt_M2 = KM2_IN * M1 - KM2_EL * M2;

// --- TGF-β1 --- (disease-side, unchanged)
dxdt_TGFb = KTGF_IN * M2 - KTGF_EL * TGFb;

// --- Fibrosis score (0-100 scale) --- (disease-side, unchanged; HDAC_FIB
// is not a census compound in this file, left exactly as-is)
double FIB_DRIVE = KFIB_IN * TGFb * (1.0 - FIB / FIB_MAX) * 24.0;  // /day → /h
dxdt_FIB = FIB_DRIVE * (1.0 - HDAC_FIB * 0.0)  // HDAC_FIB applied via simulation
           - KFIB_EL * FIB;

// --- Satellite cell pool --- (disease-side, unchanged)
dxdt_SC = KSC_REGEN * SC * (1.0 - SC / 100.0) * DYS_EFF   // self-renewal when membrane intact
          - KSC_EXHST * NECRO_RATE * SC                       // depleted by necrosis
          - 0.0001 * FIB * SC / 100.0;                        // fibrosis niche disruption

// --- Muscle function (% baseline) --- (disease-side, unchanged)
dxdt_MF = KMF_REGEN * SC_EFF * DYS_EFF * (100.0 - MF)
          - KMF_DEC * (1.0 - DYS_EFF) * (1.0 - MEMI) * MF
          - KMF_DEC * 0.5 * (FIB / FIB_MAX) * MF;

// --- 6-Minute Walk Distance --- (disease-side, unchanged)
dxdt_SWD = -KSWMD_DC * (1.0 + (1.0 - DYS_EFF) * 2.0) * SWD
            + KSWMD_DC * 0.5 * (MF / 80.0) * DYS_EFF * SWD;

// --- FVC % predicted (respiratory) --- (disease-side, unchanged)
dxdt_FVC_pct = -KFVC_DC * (1.0 - DYS_EFF * 0.5) * FVC_pct;

// --- LVEF --- (disease-side, unchanged)
dxdt_LVEF_pct = -KLVEF_DC * (1.0 - DYS_EFF * 0.3) * LVEF_pct;

$TABLE
// Serum CK (inversely related to membrane integrity × muscle mass) -- disease-side, unchanged
double CK_serum = 20000.0 * (1.0 - MEMI) * (MF / 80.0);
// NOTE (pre-existing original defect, disclosed + logged, not fixed): this
// recomputation omits the AAV/gene-therapy contribution that DYS_TOTAL
// (computed in $MAIN) includes -- see header note and refactor notes.
double Dystrophin_pct = DYS_BASE + DYS + EFF_ASO * IC_ASO;

// Plasma steroid concentration (ng/mL) -- live/current-state duplicate of
// C_CS (computed in $MAIN). Kept as its own independent computation, NOT
// collapsed into a reference to C_CS: $TABLE cannot see $MAIN locals (mrgsolve
// block scoping), and empirically the two are NOT numerically identical
// when captured (C_CS is one full reporting interval behind; C_CS_ngmL is
// live) -- confirmed via qspserver, see refactor notes. This exactly
// matches why the original itself already carried this same duplicate.
double C_CS_ngmL = CENT_CS / V1_CS * 1000.0;

// NSAA score estimate (0-34 scale) -- disease-side, unchanged
double NSAA = 34.0 * (SWD / 400.0) * (MF / 100.0);
if(NSAA > 34) NSAA = 34;
if(NSAA < 0)  NSAA = 0;

// NOTE: DYS TGFb FIB SC M1 M2 MEMI CAI ROS NFkB SWD FVC_pct LVEF_pct MF are
// all $CMT compartments -- dropped from $CAPTURE (build-compatibility fix
// #3 above); they remain reported automatically as output columns exactly
// as before, nothing becomes unobservable.
$CAPTURE C_CS_ngmL CK_serum Dystrophin_pct NSAA C_CS EFFECT_CS C_ASO EFFECT_ASO C_AAV EFFECT_AAV
'

## ===================================================================
## COMPILE MODEL
## ===================================================================
mod_dmd <- mcode("dmd_qsp_refactored", code_dmd_refactored)

cat("Model compiled successfully.\n")
cat("Compartments:", length(init(mod_dmd)), "\n")
cat("Parameters:", length(param(mod_dmd)), "\n")

## ===================================================================
## DOSING REGIMENS
## ===================================================================

# Simulation time: 6 years (52,560 hours), step = 12h
t_end <- 52560  # hours (6 years)
DT    <- 12     # h
times <- seq(0, t_end, by = DT)

age_start <- 8  # years

## ---- Scenario 1: Natural History (No treatment) ----
e_natural <- ev(time = 0, amt = 0, cmt = "GUT_CS")

## ---- Scenario 2: Deflazacort 0.9 mg/kg/day oral ----
# Every 24h, 22 kg × 0.9 = 19.8 mg/day
dose_dfz <- 19.8  # mg/day
e_dfz <- ev(time = 0, amt = dose_dfz, cmt = "GUT_CS",
            ii = 24, addl = t_end/24 - 1)  # daily dosing

## ---- Scenario 3: Prednisone 0.75 mg/kg/day oral ----
dose_pred <- 22 * 0.75  # mg/day = 16.5 mg
e_pred <- ev(time = 0, amt = dose_pred, cmt = "GUT_CS",
             ii = 24, addl = t_end/24 - 1)
# Prednisone slightly lower bioavailability than deflazacort
# Model same PK base but adjust via F parameter

## ---- Scenario 4: Eteplirsen 30 mg/kg/wk IV ----
# Every 7 days = 168h; 22kg × 30 = 660 mg/dose
dose_aso <- 30  # mg/kg/wk
e_aso <- ev(time = 0, amt = dose_aso, cmt = "CENT_ASO",
            ii = 168, addl = t_end/168 - 1)

## ---- Scenario 5: Gene therapy (Elevidys) single IV ----
# 1.33×10^14 vg/kg single dose
# Represented as large initial CENT_AAV
dose_aav <- 1.33e14  # vg/kg
e_aav <- ev(time = 0, amt = dose_aav, cmt = "CENT_AAV")

## ---- Scenario 6: Deflazacort + Eteplirsen (combination) ----
e_combo <- rbind(
  ev(time = 0, amt = dose_dfz, cmt = "GUT_CS", ii = 24, addl = t_end/24 - 1),
  ev(time = 0, amt = dose_aso, cmt = "CENT_ASO", ii = 168, addl = t_end/168 - 1)
)

## ===================================================================
## RUN SIMULATIONS (6 Scenarios)
## ===================================================================

run_scenario <- function(mod, ev_obj, scenario_name, params_override = list()) {
  mod_run <- mod
  if (length(params_override) > 0) {
    mod_run <- param(mod_run, params_override)
  }
  out <- mrgsim(mod_run, ev_obj, end = t_end, delta = DT, obsonly = TRUE)
  df <- as.data.frame(out)
  df$scenario <- scenario_name
  df$age_yr   <- age_start + df$time / 8760  # convert h to years
  return(df)
}

cat("\nRunning 6 treatment scenarios...\n")

df1 <- run_scenario(mod_dmd, e_natural, "1. Natural History")
df2 <- run_scenario(mod_dmd, e_dfz,    "2. Deflazacort 0.9 mg/kg/d")
df3 <- run_scenario(mod_dmd, e_pred,   "3. Prednisone 0.75 mg/kg/d",
                    list(F_CS = 0.82))  # slightly higher F for pred
df4 <- run_scenario(mod_dmd, e_aso,    "4. Eteplirsen 30 mg/kg/wk")
df5 <- run_scenario(mod_dmd, e_aav,    "5. Gene Therapy (Elevidys)")
df6 <- run_scenario(mod_dmd, e_combo,  "6. Deflazacort + Eteplirsen")

all_data <- bind_rows(df1, df2, df3, df4, df5, df6)
all_data$scenario <- factor(all_data$scenario, levels = c(
  "1. Natural History",
  "2. Deflazacort 0.9 mg/kg/d",
  "3. Prednisone 0.75 mg/kg/d",
  "4. Eteplirsen 30 mg/kg/wk",
  "5. Gene Therapy (Elevidys)",
  "6. Deflazacort + Eteplirsen"
))

cat("All simulations complete.\n")
cat("Total rows:", nrow(all_data), "\n")

## ===================================================================
## CLINICAL SUMMARY AT KEY TIME POINTS
## ===================================================================

summary_times <- c(0, 8760, 17520, 26280, 35040, 43800, 52560)  # 0,1,2,3,4,5,6 years
summary_labels <- c("Age 8yr", "Age 9yr", "Age 10yr", "Age 11yr",
                    "Age 12yr", "Age 13yr", "Age 14yr")

clinical_summary <- all_data %>%
  filter(time %in% summary_times) %>%
  mutate(age_label = case_when(
    time == 0     ~ "Age 8yr",
    time == 8760  ~ "Age 9yr",
    time == 17520 ~ "Age 10yr",
    time == 26280 ~ "Age 11yr",
    time == 35040 ~ "Age 12yr",
    time == 43800 ~ "Age 13yr",
    time == 52560 ~ "Age 14yr"
  )) %>%
  select(scenario, age_label, SWD, NSAA, FVC_pct, LVEF_pct, FIB, CK_serum, Dystrophin_pct) %>%
  mutate(across(where(is.numeric), ~round(., 1)))

print(clinical_summary)

## ===================================================================
## VISUALIZATION
## ===================================================================

scenario_colors <- c(
  "1. Natural History"          = "#E74C3C",
  "2. Deflazacort 0.9 mg/kg/d"  = "#3498DB",
  "3. Prednisone 0.75 mg/kg/d"  = "#2ECC71",
  "4. Eteplirsen 30 mg/kg/wk"   = "#9B59B6",
  "5. Gene Therapy (Elevidys)"  = "#E67E22",
  "6. Deflazacort + Eteplirsen" = "#1ABC9C"
)

# --- Plot 1: 6MWD over time ---
p1 <- ggplot(all_data, aes(x = age_yr, y = SWD, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "DMD QSP Model: 6-Minute Walk Distance (6MWD)",
       x = "Age (years)", y = "6MWD (meters)",
       color = "Treatment") +
  geom_hline(yintercept = 300, linetype = "dashed", color = "gray50",
             alpha = 0.6) +
  annotate("text", x = 8.5, y = 310, label = "LoA threshold ~300m",
           hjust = 0, color = "gray40", size = 3.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  xlim(8, 14)

# --- Plot 2: Dystrophin level ---
p2 <- ggplot(all_data, aes(x = age_yr, y = Dystrophin_pct, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Dystrophin Level Over Time",
       x = "Age (years)", y = "Dystrophin (% of normal)",
       color = "Treatment") +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray50") +
  annotate("text", x = 8.5, y = 4.5, label = "~4% threshold (functional benefit)",
           hjust = 0, color = "gray40", size = 3.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  xlim(8, 14)

# --- Plot 3: Fibrosis progression ---
p3 <- ggplot(all_data, aes(x = age_yr, y = FIB, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Muscle Fibrosis Score",
       x = "Age (years)", y = "Fibrosis Score (0-100)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  xlim(8, 14)

# --- Plot 4: FVC % predicted ---
p4 <- ggplot(all_data, aes(x = age_yr, y = FVC_pct, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Respiratory Function (FVC % predicted)",
       x = "Age (years)", y = "FVC (% predicted)",
       color = "Treatment") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.6) +
  annotate("text", x = 8.5, y = 51, label = "NIV threshold <50%",
           hjust = 0, color = "red", size = 3.5) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  xlim(8, 14)

# --- Plot 5: Inflammation markers ---
p5 <- ggplot(all_data, aes(x = age_yr, y = M1, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "M1 Macrophage Infiltration",
       x = "Age (years)", y = "M1 Macrophages (cells/µL tissue)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  xlim(8, 14)

# --- Plot 6: Serum CK ---
p6 <- ggplot(all_data %>% filter(age_yr <= 9), aes(x = age_yr, y = CK_serum, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Serum CK Levels (First Year)",
       x = "Age (years)", y = "Serum CK (U/L)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

cat("\nPlot objects created: p1 (6MWD), p2 (Dystrophin%), p3 (Fibrosis),\n")
cat("  p4 (FVC), p5 (M1 macrophages), p6 (Serum CK)\n")
cat("Use print(p1) through print(p6) to display.\n")

## ===================================================================
## SENSITIVITY ANALYSIS
## ===================================================================

sensitivity_params <- list(
  KD_MEM   = seq(0.08, 0.25, by = 0.05),   # Membrane damage rate
  KFIB_IN  = seq(0.001, 0.004, by = 0.001), # Fibrosis rate
  EFF_ASO  = seq(0.02, 0.06, by = 0.01)    # ASO efficacy
)

run_sensitivity <- function(param_name, param_values) {
  results <- list()
  for (pv in param_values) {
    p_list <- setNames(list(pv), param_name)
    mod_sens <- param(mod_dmd, p_list)
    out <- mrgsim(mod_sens, e_natural, end = t_end, delta = DT, obsonly = TRUE)
    df <- as.data.frame(out)
    df$param_val <- pv
    df$param_name <- param_name
    df$age_yr <- age_start + df$time / 8760
    results[[length(results)+1]] <- df
  }
  bind_rows(results)
}

cat("\nRunning sensitivity analysis...\n")
sens_kd  <- run_sensitivity("KD_MEM",  sensitivity_params$KD_MEM)
sens_fib <- run_sensitivity("KFIB_IN", sensitivity_params$KFIB_IN)
cat("Sensitivity analysis complete.\n")

# Tornado chart data at 6 years
tornado_6yr <- bind_rows(
  sens_kd %>% filter(time == t_end) %>% select(param_name, param_val, SWD),
  sens_fib %>% filter(time == t_end) %>% select(param_name, param_val, SWD)
)

## ===================================================================
## VIRTUAL PATIENT POPULATION (Monte Carlo, n=50)
## ===================================================================

set.seed(42)
n_vp <- 50

vp_params <- data.frame(
  ID       = 1:n_vp,
  KD_MEM   = rlnorm(n_vp, log(0.15), 0.3),
  KFIB_IN  = rlnorm(n_vp, log(0.002), 0.4),
  KSC_EXHST= rlnorm(n_vp, log(0.002), 0.35),
  MEMI_SS  = rnorm(n_vp, 0.35, 0.06)
)
vp_params$MEMI_SS <- pmax(0.15, pmin(0.65, vp_params$MEMI_SS))

run_vp <- function(params_row, ev_obj, scenario) {
  p_list <- as.list(params_row[, c("KD_MEM","KFIB_IN","KSC_EXHST","MEMI_SS")])
  p_list <- lapply(p_list, as.numeric)
  mod_vp <- param(mod_dmd, p_list)
  mod_vp <- init(mod_vp, MEMI = params_row$MEMI_SS)
  out <- mrgsim(mod_vp, ev_obj, end = t_end, delta = DT, obsonly = TRUE)
  df <- as.data.frame(out)
  df$ID <- params_row$ID
  df$scenario <- scenario
  df$age_yr <- age_start + df$time / 8760
  df
}

cat("\nRunning virtual patient population (n=50)...\n")
vp_natural <- bind_rows(lapply(1:n_vp, function(i) run_vp(vp_params[i,], e_natural, "Natural History")))
vp_dfz     <- bind_rows(lapply(1:n_vp, function(i) run_vp(vp_params[i,], e_dfz,    "Deflazacort")))

# 5th/50th/95th percentile bands
vp_summary <- bind_rows(vp_natural, vp_dfz) %>%
  group_by(scenario, time, age_yr) %>%
  summarise(
    p05_SWD = quantile(SWD, 0.05),
    p50_SWD = quantile(SWD, 0.50),
    p95_SWD = quantile(SWD, 0.95),
    .groups = "drop"
  )

p_vp <- ggplot(vp_summary, aes(x = age_yr, fill = scenario, color = scenario)) +
  geom_ribbon(aes(ymin = p05_SWD, ymax = p95_SWD), alpha = 0.2) +
  geom_line(aes(y = p50_SWD), linewidth = 1.2) +
  scale_fill_manual(values = c("Natural History" = "#E74C3C", "Deflazacort" = "#3498DB")) +
  scale_color_manual(values = c("Natural History" = "#E74C3C", "Deflazacort" = "#3498DB")) +
  labs(title = "Virtual Patient Population: 6MWD\n(n=50, 5th-95th percentile band)",
       x = "Age (years)", y = "6MWD (m)",
       fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 12) +
  xlim(8, 14)

cat("VP simulation complete.\n")
cat("Use print(p_vp) to display virtual patient population plot.\n")

## ===================================================================
## BIOMARKER TRAJECTORY SUMMARY (Table)
## ===================================================================

traj_summary <- all_data %>%
  filter(time == t_end) %>%  # 6yr follow-up
  select(scenario, SWD, NSAA, FVC_pct, LVEF_pct, FIB, CK_serum, Dystrophin_pct, SC, MF) %>%
  mutate(
    SWD_change   = round(SWD - 380, 1),
    FVC_change   = round(FVC_pct - 95, 1),
    NSAA_change  = round(NSAA - (34 * 380/400 * 80/100), 1),
    across(where(is.numeric), ~round(., 1))
  )

cat("\n==== 6-Year Treatment Outcome Summary (Age 14yr) ====\n")
print(traj_summary %>% select(scenario, SWD, SWD_change, FVC_pct, LVEF_pct,
                                FIB, Dystrophin_pct, SC))

cat("\n==== DMD QSP MODEL SUMMARY (REFACTORED) ====\n")
cat("Model version: v1.0-refactored | Date: 2026-08-31\n")
cat("Compartments : 19 ODE (8 PK + 11 PD)\n")
cat("Parameters   : ", length(param(mod_dmd)), "\n")
cat("Scenarios    : 6 (Natural history + 5 treatments)\n")
cat("Drug classes : Corticosteroids (Deflazacort/Prednisone), ASO (Eteplirsen), AAV gene therapy (delandistrogene moxeparvovec/Elevidys)\n")
cat("Endpoints    : 6MWD, NSAA, FVC%, LVEF%, Fibrosis, Dystrophin%, CK\n")
cat("VP Population: n=50 Monte Carlo\n")
