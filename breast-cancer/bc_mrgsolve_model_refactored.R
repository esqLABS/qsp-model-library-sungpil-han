# =============================================================================
# Breast Cancer QSP Model — mrgsolve ODE System (PK/PD REFACTORED)
# =============================================================================
# Covers ER+/HER2+/TNBC biology with drug PK/PD for 6+ agents
# Parameters calibrated to:
#   PALOMA-2    (palbociclib + letrozole, ER+/HER2-)
#   MONALEESA-2 (ribociclib + letrozole, ER+/HER2-)
#   CLEOPATRA   (trastuzumab + docetaxel + pertuzumab, HER2+)
#   KEYNOTE-522 (pembrolizumab + chemo, TNBC)
#   OlympiAD    (olaparib, BRCAm)
#
# [refactor] Per FORK_WORKFLOW_GUIDE.md Part 2 (pluggable PK, named Hill
# interface), the PK/PD blocks for Letrozole (LETRO), Olaparib (OLAP),
# Palbociclib (PALBO) and Trastuzumab (TRAS) were renamed to the fork's
# naming convention (C_<STEM> exposed concentration, EFFECT_<STEM> named
# disease effect). No other compound blocks exist in this file. Full
# rationale, per-compound archetype, and verification results are in
# bc_refactor_notes.md.
#
# Renamed (values unchanged from the original):
#   ka_palbo/Vd_palbo/CL_palbo         -> KA_PALBO/V1_PALBO/CL_PALBO
#   Emax_CDK/EC50_CDK/hill_CDK         -> EMAX_PALBO/EC50_PALBO/GAMMA_PALBO
#   ka_letro/Vd_letro/CL_letro         -> KA_LETRO/V1_LETRO/CL_LETRO
#   KAI                                 -> EC50_LETRO
#   CL_tras/Vd1_tras/Vd2_tras/Q_tras   -> CL_TRAS/V1_TRAS/V2_TRAS/Q_TRAS
#   Emax_HER2/EC50_HER2                -> EMAX_TRAS/EC50_TRAS
#   ka_olap/Vd_olap/CL_olap            -> KA_OLAP/V1_OLAP/CL_OLAP
#   Cp_palbo/Cp_letro/Cp_tras/Cp_olap  -> C_PALBO/C_LETRO/C_TRAS/C_OLAP
#   CDK_inh/HER2_block/letro_inh       -> EFFECT_PALBO/EFFECT_TRAS/EFFECT_LETRO
# New (not in original, named Hill-interface parameters; EMAX/GAMMA=1 make
# explicit a shape the original already had implicitly as a plain ratio):
#   EMAX_LETRO=1, GAMMA_LETRO=1, GAMMA_TRAS=1
# Dropped (dead code, dxdt=0, never read by anything -- verified to change
# nothing): PERI_PALBO compartment; the unused Cp2_tras local double.
# No EFFECT_OLAP: the original never connects Cp_olap to any disease
# equation at all (see refactor notes) -- there is nothing to rename.
#
# Pre-existing upstream build defect (unrelated to any compound's own PK):
# $CAPTURE duplicated 12 compartment names, which mrgsolve 2.0.1 rejects
# outright. Fixed syntax-only here (compartment names removed from
# $CAPTURE; mrgsolve always reports compartment state regardless of
# $CAPTURE, so this changes nothing about what is reported). Logged as
# translations/UPSTREAM_ISSUES.md #57. The checked-in original
# (bc_mrgsolve_model.R) is untouched and still carries this defect.
# =============================================================================

library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# mrgsolve model code string
# =============================================================================

bc_model_code <- '
$PROB
Breast Cancer QSP Model
-----------------------
Subtypes: ER+/HER2+/TNBC
Drug PK/PD: palbociclib, letrozole, trastuzumab, olaparib
Calibrated to: PALOMA-2, MONALEESA-2, CLEOPATRA, KEYNOTE-522, OlympiAD
[refactor] LETRO/OLAP/PALBO/TRAS PK+PD blocks renamed to the forks
plumbable-PK convention (C_<STEM>/EFFECT_<STEM>); values unchanged from
the original. See bc_refactor_notes.md. All other compounds/params
(there are none besides these four in this file) untouched.

$PARAM
// ---- Tumor biology ----
kprol      = 0.0008    // tumor cell proliferation rate (1/hr), doubling time ~36 days
kdeath     = 0.0002    // baseline tumor cell death rate (1/hr)
kCSC       = 0.00005   // CSC self-renewal rate (1/hr)
kDiff      = 0.0005    // CSC differentiation to bulk tumor (1/hr)
Kmax       = 1000.0    // tumor carrying capacity (cm^3 equivalent units)
kmet       = 0.00001   // metastasis seeding rate

// ---- ER signaling ----
E2_base    = 100.0     // baseline estradiol (pmol/L; postmenopausal ~20, premenopausal ~300)
Emax_E2    = 1.5       // maximum proliferative effect of estradiol
EC50_E2    = 50.0      // EC50 for estradiol effect (pmol/L)
kAI_aro    = 0.01      // aromatase inhibition rate constant [unused in original -- kept as-is, out of scope]

// ---- CDK4/6 inhibitor PD (palbociclib) ----
// [refactor] Emax_CDK/EC50_CDK/hill_CDK -> EMAX_PALBO/EC50_PALBO/GAMMA_PALBO (rename only, same values)
EMAX_PALBO = 0.85      // maximum CDK4/6 inhibition (85% max RB pathway blockade)
EC50_PALBO = 100.0     // EC50 for CDK4/6 inhibitor (ng/mL)
GAMMA_PALBO = 1.5      // Hill coefficient for CDK4/6 inhibition

// ---- HER2 pathway (trastuzumab) ----
// [refactor] Emax_HER2/EC50_HER2 -> EMAX_TRAS/EC50_TRAS (rename, same values);
// GAMMA_TRAS is new (=1), making explicit the originals implicit linear Hill shape
kHER2      = 0.3       // HER2 signaling amplification factor
EMAX_TRAS  = 0.7       // max anti-HER2 drug effect on proliferation
EC50_TRAS  = 50.0      // EC50 for anti-HER2 (ug/mL)
GAMMA_TRAS = 1.0        // Hill coefficient (original had no explicit Hill term -- linear ratio)

// ---- Aromatase inhibition (letrozole) ----
// [refactor] KAI -> EC50_LETRO (rename, same value); EMAX_LETRO/GAMMA_LETRO are
// new (=1 each), making explicit the originals implicit linear Hill shape
EC50_LETRO = 50.0      // IC50 of aromatase inhibitor on E2 (pmol/L drug conc equivalent)
EMAX_LETRO = 1.0       // maximum aromatase-inhibition effect (original had no explicit Emax -- linear ratio)
GAMMA_LETRO = 1.0       // Hill coefficient (original had no explicit Hill term -- linear ratio)

// ---- Immune parameters ----
kCD8_recruit = 0.005   // CD8+ T cell recruitment rate
kCD8_kill    = 0.002   // CD8+ T cell kill rate of tumor cells
kTreg_sup    = 0.003   // Treg suppression rate of CD8+
kPDL1_ind    = 0.001   // PD-L1 induction by IFNgamma [unused in original -- kept as-is, out of scope]
EC50_PD1     = 10.0    // EC50 for anti-PD1 effect [unused in original -- kept as-is, out of scope]

// ---- Drug PK: Palbociclib (oral) ----
// [refactor] ka_palbo/Vd_palbo/CL_palbo -> KA_PALBO/V1_PALBO/CL_PALBO (rename only, same values)
KA_PALBO   = 0.5       // absorption rate (1/hr)
V1_PALBO   = 2583.0    // volume of distribution (L)
CL_PALBO   = 63.0      // clearance (L/hr)

// ---- Drug PK: Letrozole (oral) ----
// [refactor] ka_letro/Vd_letro/CL_letro -> KA_LETRO/V1_LETRO/CL_LETRO (rename only, same values)
KA_LETRO   = 0.8       // absorption rate (1/hr)
V1_LETRO   = 187.0     // volume of distribution (L)
CL_LETRO   = 2.1       // clearance (L/hr)

// ---- Drug PK: Trastuzumab (IV, 2-compartment) ----
// [refactor] CL_tras/Vd1_tras/Vd2_tras/Q_tras -> CL_TRAS/V1_TRAS/V2_TRAS/Q_TRAS (rename only, same values)
CL_TRAS    = 0.225     // clearance (L/day)
V1_TRAS    = 3.63      // central volume (L)
V2_TRAS    = 2.78      // peripheral volume (L)
Q_TRAS     = 0.747     // inter-compartmental clearance (L/day)

// ---- Drug PK: Olaparib (oral) ----
// [refactor] ka_olap/Vd_olap/CL_olap -> KA_OLAP/V1_OLAP/CL_OLAP (rename only, same values).
// No PD/effect parameters: the original never links Cp_olap to any disease
// equation -- see refactor notes.
KA_OLAP    = 0.9       // absorption rate (1/hr)
V1_OLAP    = 158.0     // volume of distribution (L)
CL_OLAP    = 8.6       // clearance (L/hr)

// ---- Biomarker dynamics ----
kKi67_on     = 0.01    // Ki-67 response rate constant [unused in original -- kept as-is, out of scope]
kCA153_prod  = 0.005   // CA15-3 production rate by tumor cells
kCA153_elim  = 0.02    // CA15-3 elimination rate
CA153_base   = 20.0    // baseline CA15-3 (U/mL)

$CMT
GUT_PALBO    // palbociclib gut (absorption depot)
CENT_PALBO   // palbociclib central (plasma)
// [refactor] PERI_PALBO removed: in the original it was declared in $CMT but
// dxdt_PERI_PALBO=0 and it was never read by Cp_palbo or anything else --
// an always-zero, functionally dead compartment. Dropping it changes no
// computed value (verified below).

GUT_LETRO    // letrozole gut (absorption depot)
CENT_LETRO   // letrozole central (plasma)

CENT_TRAS    // trastuzumab central (plasma, mg)
PERI_TRAS    // trastuzumab peripheral

GUT_OLAP     // olaparib gut (absorption depot)
CENT_OLAP    // olaparib central (plasma)

TUMOR        // bulk tumor volume (normalized units, 1 unit = ~1 cm^3)
CSC          // cancer stem cell pool (relative units)
ER_SIGNAL    // ER signaling activity (0-1 scale)
CDK46_ACT    // CDK4/6 activity (0-1 = fully inhibited to fully active)
HER2_SIGNAL  // HER2 signaling activity (relative)
PD_L1        // PD-L1 expression on tumor (relative, 0-1)
CD8_EFF      // CD8+ effector T cells (relative units)
TREG         // Regulatory T cells (relative units)
E2_PLASMA    // Estradiol plasma concentration (pmol/L)
AROMATASE    // Aromatase enzyme activity (0-1, 1 = full activity)
Ki67         // Ki-67 proliferation index (%)
CA153        // CA15-3 tumor marker (U/mL)

$GLOBAL
// [refactor] the four refactored compounds exposed concentration/effect
// terms, predeclared here rather than in $PARAM: mrgsolve 2.0.1 compiles
// $PARAM members as read-only references inside $ODE, so a value that must
// be recomputed every timestep from state cannot also live in $PARAM (same
// constraint documented in the AMD and membranous-nephropathy refactors).
// These are visible in every simulations output via $CAPTURE and in
// /model_manifests outputPaths.
double C_PALBO, C_LETRO, C_TRAS, C_OLAP;
double EFFECT_PALBO, EFFECT_LETRO, EFFECT_TRAS;
// No EFFECT_OLAP: the original never connects Cp_olap to any disease
// equation -- see refactor notes.

$MAIN
// Initial conditions
TUMOR_0      = 100.0;    // initial tumor ~1 cm^3 mass equivalent
CSC_0        = 5.0;      // initial cancer stem cell pool
ER_SIGNAL_0  = 1.0;      // baseline ER signaling fully active
CDK46_ACT_0  = 1.0;      // CDK4/6 fully active at baseline
HER2_SIGNAL_0 = 0.3;     // low HER2 for ER+ subtype (HER2- default)
E2_PLASMA_0  = E2_base;  // baseline estradiol from parameter
AROMATASE_0  = 1.0;      // full aromatase activity at baseline
PD_L1_0      = 0.2;      // baseline PD-L1 expression
CD8_EFF_0    = 0.5;      // baseline CD8+ effector T cells
TREG_0       = 0.2;      // baseline regulatory T cells
Ki67_0       = 30.0;     // 30% baseline Ki-67 proliferation index
CA153_0      = CA153_base; // baseline CA15-3 = 20 U/mL

$ODE
// -----------------------------------------------------------------------
// PK: Palbociclib (oral, depot + central; peripheral dropped, see $CMT note)
// [refactor] Cp_palbo/k10_palbo -> C_PALBO/K10_PALBO (rename only)
// -----------------------------------------------------------------------
C_PALBO = CENT_PALBO / V1_PALBO;            // ng/mL (if dose in mg, Vd in L)
double K10_PALBO = CL_PALBO / V1_PALBO;
dxdt_GUT_PALBO  = -KA_PALBO * GUT_PALBO;
dxdt_CENT_PALBO = KA_PALBO * GUT_PALBO - K10_PALBO * CENT_PALBO;

// -----------------------------------------------------------------------
// PK: Letrozole (oral, depot + central)
// [refactor] Cp_letro/k10_letro -> C_LETRO/K10_LETRO (rename only)
// -----------------------------------------------------------------------
C_LETRO = CENT_LETRO / V1_LETRO;            // ng/mL
double K10_LETRO = CL_LETRO / V1_LETRO;
dxdt_GUT_LETRO  = -KA_LETRO * GUT_LETRO;
dxdt_CENT_LETRO = KA_LETRO * GUT_LETRO - K10_LETRO * CENT_LETRO;

// -----------------------------------------------------------------------
// PK: Trastuzumab (IV, 2-compartment, PK params in day^-1 -> convert to hr^-1)
// [refactor] Cp_tras/Cp2_tras/k10_tras/k12_tras/k21_tras ->
// C_TRAS/(dropped, unused)/K10_TRAS/K12_TRAS/K21_TRAS (rename only; the
// originals Cp2_tras double was computed but never read anywhere -- dropped
// as dead code, same as PERI_PALBO above)
// -----------------------------------------------------------------------
C_TRAS  = CENT_TRAS / V1_TRAS;             // mg/L = ug/mL
double K10_TRAS = (CL_TRAS  / V1_TRAS) / 24.0;     // convert day^-1 to hr^-1
double K12_TRAS = (Q_TRAS   / V1_TRAS) / 24.0;
double K21_TRAS = (Q_TRAS   / V2_TRAS) / 24.0;
dxdt_CENT_TRAS  = -K10_TRAS * CENT_TRAS - K12_TRAS * CENT_TRAS + K21_TRAS * PERI_TRAS;
dxdt_PERI_TRAS  =  K12_TRAS * CENT_TRAS - K21_TRAS * PERI_TRAS;

// -----------------------------------------------------------------------
// PK: Olaparib (oral, depot + central)
// [refactor] Cp_olap/k10_olap -> C_OLAP/K10_OLAP (rename only)
// -----------------------------------------------------------------------
C_OLAP  = CENT_OLAP / V1_OLAP;             // ng/mL
double K10_OLAP = CL_OLAP / V1_OLAP;
dxdt_GUT_OLAP   = -KA_OLAP * GUT_OLAP;
dxdt_CENT_OLAP  = KA_OLAP * GUT_OLAP - K10_OLAP * CENT_OLAP;

// -----------------------------------------------------------------------
// Estradiol dynamics (aromatase-regulated) -- unchanged, disease-side
// -----------------------------------------------------------------------
// E2 produced proportionally to aromatase activity, cleared at constant rate
dxdt_E2_PLASMA  = (E2_base * AROMATASE - E2_PLASMA * 0.05);

// -----------------------------------------------------------------------
// Aromatase activity (inhibited by letrozole via competitive inhibition)
// [refactor] letro_inh -> EFFECT_LETRO. With EMAX_LETRO=GAMMA_LETRO=1 this
// collapses algebraically to the originals C_LETRO/(C_LETRO+EC50_LETRO) --
// a rename, not a refit.
// -----------------------------------------------------------------------
double letro_h = pow(C_LETRO, GAMMA_LETRO);
double ec50_letro_h = pow(EC50_LETRO, GAMMA_LETRO);
EFFECT_LETRO = EMAX_LETRO * letro_h / (ec50_letro_h + letro_h);
dxdt_AROMATASE  = 0.01 * (1.0 - AROMATASE) - 0.01 * AROMATASE * EFFECT_LETRO * 10.0;

// -----------------------------------------------------------------------
// ER signaling (driven by E2; first-order approach to E2-driven steady state)
// -- unchanged, disease-side
// -----------------------------------------------------------------------
double E2_effect = Emax_E2 * E2_PLASMA / (EC50_E2 + E2_PLASMA);
dxdt_ER_SIGNAL  = 0.05 * (E2_effect - ER_SIGNAL);

// -----------------------------------------------------------------------
// CDK4/6 activity (inhibited by palbociclib via Hill equation)
// [refactor] CDK_inh -> EFFECT_PALBO (rename only, identical arithmetic)
// -----------------------------------------------------------------------
double palbo_h   = pow(C_PALBO, GAMMA_PALBO);
double ec50_palbo_h    = pow(EC50_PALBO, GAMMA_PALBO);
EFFECT_PALBO   = EMAX_PALBO * palbo_h / (ec50_palbo_h + palbo_h);
dxdt_CDK46_ACT  = 0.1 * ((1.0 - EFFECT_PALBO) - CDK46_ACT);

// -----------------------------------------------------------------------
// HER2 signaling (inhibited by trastuzumab)
// [refactor] HER2_block -> EFFECT_TRAS (rename only; GAMMA_TRAS=1 makes this
// literally the originals linear ratio shape)
// -----------------------------------------------------------------------
double tras_h = pow(C_TRAS, GAMMA_TRAS);
double ec50_tras_h = pow(EC50_TRAS, GAMMA_TRAS);
EFFECT_TRAS = EMAX_TRAS * tras_h / (ec50_tras_h + tras_h);
dxdt_HER2_SIGNAL = 0.05 * ((kHER2 * (1.0 - EFFECT_TRAS)) - HER2_SIGNAL);

// -----------------------------------------------------------------------
// Immune microenvironment dynamics -- unchanged
// -----------------------------------------------------------------------
double tumor_signal = TUMOR / (TUMOR + 100.0);   // tumor burden drives immune activation
double PDL1_block   = 0.0;                        // placeholder: set >0 for anti-PD1 drugs

// PD-L1 expression (induced by tumor microenvironment signals)
dxdt_PD_L1  = 0.01 * (tumor_signal - PD_L1);

// CD8+ effector T cells: recruited by tumor signal, suppressed by PD-L1 and Tregs
double PD_L1_eff  = PD_L1 * (1.0 - PDL1_block);
dxdt_CD8_EFF = kCD8_recruit * tumor_signal * (1.0 - PD_L1_eff)
               - kTreg_sup * TREG * CD8_EFF
               - 0.01 * CD8_EFF;

// Regulatory T cells: recruited by tumor signal, natural turnover
dxdt_TREG   = 0.005 * tumor_signal - 0.008 * TREG;

// -----------------------------------------------------------------------
// Tumor dynamics -- central equation -- unchanged
// -----------------------------------------------------------------------
double prol_rate   = kprol * ER_SIGNAL * CDK46_ACT;       // ER-driven, CDK-gated proliferation
double HER2_contrib = kprol * 0.5 * HER2_SIGNAL;          // additional HER2-driven proliferation
double immune_kill  = kCD8_kill * CD8_EFF;                 // immune-mediated cytotoxicity

dxdt_TUMOR = (prol_rate + HER2_contrib) * TUMOR * (1.0 - TUMOR / Kmax)
             - (kdeath + immune_kill) * TUMOR;

// -----------------------------------------------------------------------
// Cancer stem cell dynamics (self-renewal and differentiation) -- unchanged
// -----------------------------------------------------------------------
dxdt_CSC = kCSC * CSC - kDiff * CSC;

// -----------------------------------------------------------------------
// Biomarkers -- unchanged
// -----------------------------------------------------------------------
// Ki-67: tracks CDK4/6 activity and ER signaling (proliferation index)
dxdt_Ki67 = 0.1 * (100.0 * CDK46_ACT * ER_SIGNAL - Ki67);

// CA15-3: produced by tumor, eliminated with baseline offset
dxdt_CA153 = kCA153_prod * TUMOR - kCA153_elim * (CA153 - CA153_base);

$TABLE
// Derived outputs computed at each output time step -- unchanged
double TGR      = (TUMOR > 1e-6) ? log(TUMOR / 100.0) : -10.0; // tumor growth ratio vs baseline
double SPD      = TUMOR;                                          // sum of product diameters proxy
double response = (TUMOR < 30.0) ? 1.0 : 0.0;                   // partial response threshold

$CAPTURE
C_PALBO C_LETRO C_TRAS C_OLAP
EFFECT_PALBO EFFECT_LETRO EFFECT_TRAS
TGR response
'


# =============================================================================
# Compile model
# =============================================================================

mod <- mread("bc", tempdir(), bc_model_code)

# =============================================================================
# Simulation parameters
# =============================================================================

end_time <- 8760   # 1 year in hours
dt       <- 24     # daily output

# =============================================================================
# Treatment event objects
# =============================================================================

# Scenario 1: Letrozole monotherapy (2.5 mg daily oral, ER+)
ev1 <- ev(amt = 2.5,  cmt = "GUT_LETRO",  ii = 24, addl = 364, time = 0)

# Scenario 2: Palbociclib (125 mg, 21-days-on/7-days-off) + Letrozole (PALOMA-2)
# Simplified as continuous dosing to represent average CDK4/6 inhibition over cycles
ev2_palbo <- ev(amt = 125, cmt = "GUT_PALBO", ii = 24, addl = 20, time = 0)
ev2_letro <- ev(amt = 2.5, cmt = "GUT_LETRO", ii = 24, addl = 364, time = 0)
ev2 <- c(ev2_palbo, ev2_letro)

# Scenario 3: Ribociclib (600 mg) + Letrozole (MONALEESA-2)
# Ribociclib modeled via same CDK4/6 compartment as palbociclib, dose scaled
# 600 mg ribociclib / 4 (approx PK equivalence factor) = 150 mg palbociclib-equivalent
ev3_ribo  <- ev(amt = 150, cmt = "GUT_PALBO", ii = 24, addl = 20, time = 0)
ev3_letro <- ev(amt = 2.5, cmt = "GUT_LETRO", ii = 24, addl = 364, time = 0)
ev3 <- c(ev3_ribo, ev3_letro)

# Scenario 4: Trastuzumab + Docetaxel (HER2+, CLEOPATRA-like, 70 kg patient)
# Loading dose 8 mg/kg = 560 mg, maintenance 6 mg/kg = 420 mg q3w (504 hr)
ev4_load  <- ev(amt = 560, cmt = "CENT_TRAS", time = 0)
ev4_maint <- ev(amt = 420, cmt = "CENT_TRAS", ii = 504, addl = 16, time = 504)
ev4 <- c(ev4_load, ev4_maint)

# Scenario 5: Pembrolizumab + Chemotherapy (KEYNOTE-522, TNBC)
# Modeled via parameter overrides boosting CD8+ recruitment/kill (immune checkpoint release)
# Events: letrozole as a placeholder carrier (minimal PD effect in TNBC phenotype)
ev5 <- ev(amt = 2.5, cmt = "GUT_LETRO", ii = 24, addl = 364, time = 0)

# Scenario 6: Olaparib monotherapy (OlympiAD, BRCAm, 300 mg BID)
ev6 <- ev(amt = 300, cmt = "GUT_OLAP", ii = 12, addl = 729, time = 0)

# =============================================================================
# Run function
# =============================================================================

run_scenario <- function(model, events, scenario_name, param_override = list()) {
  mod_run <- param(model, param_override)
  out     <- mrgsim(mod_run, events = events, end = end_time, delta = dt, digits = 4)
  df      <- as.data.frame(out)
  df$scenario <- scenario_name
  df
}

# =============================================================================
# Execute all six scenarios
# =============================================================================

res1 <- run_scenario(
  mod, ev1,
  "Letrozole mono (ER+)"
)

res2 <- run_scenario(
  mod, ev2,
  "Palbociclib + Letrozole (PALOMA-2)"
)

res3 <- run_scenario(
  mod, ev3,
  "Ribociclib + Letrozole (MONALEESA-2)"
)

res4 <- run_scenario(
  mod, ev4,
  "Trastuzumab + Docetaxel (CLEOPATRA, HER2+)",
  list(
    kHER2        = 1.0,   # HER2-amplified tumor phenotype
    ER_SIGNAL_0  = 0.1,   # ER-low in HER2+ subtype
    HER2_SIGNAL_0 = 1.0   # strong HER2 signaling at baseline
  )
)

res5 <- run_scenario(
  mod, ev5,
  "Pembrolizumab + Chemo (KEYNOTE-522, TNBC)",
  list(
    kCD8_recruit  = 0.02,   # checkpoint blockade: enhanced CD8 recruitment
    kCD8_kill     = 0.008,  # stronger cytotoxic activity post-PD1 block
    kprol         = 0.001,  # TNBC typically higher proliferation rate
    ER_SIGNAL_0   = 0.05,   # ER-negative subtype
    HER2_SIGNAL_0 = 0.1     # HER2-negative subtype
  )
)

res6 <- run_scenario(
  mod, ev6,
  "Olaparib (OlympiAD, BRCAm)",
  list(
    kprol   = 0.0012,  # BRCA-mutated tumors: enhanced replication stress
    EC50_PALBO = 200.0 # [refactor] was EC50_CDK; olaparib: PARP inhibition (CDK pathway less relevant, reduced sensitivity)
  )
)

all_results <- bind_rows(res1, res2, res3, res4, res5, res6)

# =============================================================================
# Plot 1: Tumor Volume Dynamics by Treatment Regimen
# =============================================================================

p1 <- ggplot(all_results, aes(x = time / 24 / 7, y = TUMOR, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title   = "Tumor Volume Dynamics by Treatment Regimen",
    x       = "Time (weeks)",
    y       = "Tumor Volume (relative units)",
    color   = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  ) +
  guides(color = guide_legend(ncol = 2))

print(p1)

# =============================================================================
# Plot 2: Ki-67 Proliferation Index Over Time
# =============================================================================

p2 <- ggplot(all_results, aes(x = time / 24 / 7, y = Ki67, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Ki-67 Proliferation Index Over Time",
    x     = "Time (weeks)",
    y     = "Ki-67 Index (%)",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2))

print(p2)

# =============================================================================
# Plot 3: CDK4/6 Pathway Inhibition Depth (ER+ regimens)
# =============================================================================

er_scenarios <- c(
  "Letrozole mono (ER+)",
  "Palbociclib + Letrozole (PALOMA-2)",
  "Ribociclib + Letrozole (MONALEESA-2)"
)

p3 <- ggplot(
  all_results %>% filter(scenario %in% er_scenarios),
  aes(x = time / 24, y = (1 - CDK46_ACT) * 100, color = scenario)
) +
  geom_line(size = 1.1) +
  labs(
    title = "CDK4/6 Pathway Inhibition Depth",
    x     = "Time (days)",
    y     = "CDK4/6 Inhibition (%)",
    color = "Regimen"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p3)

# =============================================================================
# Plot 4: CA15-3 Tumor Marker Dynamics
# =============================================================================

p4 <- ggplot(all_results, aes(x = time / 24 / 7, y = CA153, color = scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 35, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate(
    "text",
    x     = max(all_results$time) / 24 / 7 * 0.8,
    y     = 37,
    label = "ULN = 35 U/mL",
    color = "red",
    size  = 3
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "CA15-3 Tumor Marker Dynamics",
    x     = "Time (weeks)",
    y     = "CA15-3 (U/mL)",
    color = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2))

print(p4)

# =============================================================================
# Plot 5: PK Profiles — Palbociclib + Letrozole (First 4 Weeks)
# =============================================================================

p5 <- ggplot(
  res2 %>% filter(time <= 24 * 28),
  aes(x = time / 24)
) +
  geom_line(aes(y = C_PALBO,       color = "Palbociclib (ng/mL)"),      size = 1) +
  geom_line(aes(y = C_LETRO * 10,  color = "Letrozole (ng/mL ×10)"), size = 1) +
  scale_color_manual(
    values = c(
      "Palbociclib (ng/mL)"          = "#1565C0",
      "Letrozole (ng/mL ×10)"   = "#B71C1C"
    )
  ) +
  labs(
    title = "PK Profiles: Palbociclib + Letrozole (First 4 Weeks)",
    x     = "Time (days)",
    y     = "Plasma Concentration",
    color = "Drug"
  ) +
  theme_bw(base_size = 12)

print(p5)

# =============================================================================
# Plot 6: Immune Microenvironment Dynamics (Pembrolizumab + Chemo, TNBC)
# =============================================================================

p6 <- ggplot(res5, aes(x = time / 24 / 7)) +
  geom_line(aes(y = CD8_EFF * 100, color = "CD8+ T cells"),       size = 1.1) +
  geom_line(aes(y = TREG    * 100, color = "Tregs"),               size = 1.1) +
  geom_line(aes(y = PD_L1   * 100, color = "PD-L1 expression"),    size = 1.1) +
  scale_color_manual(
    values = c(
      "CD8+ T cells"     = "#1B5E20",
      "Tregs"            = "#B71C1C",
      "PD-L1 expression" = "#E65100"
    )
  ) +
  labs(
    title = "Immune Microenvironment Dynamics\n(Pembrolizumab + Chemo, TNBC)",
    x     = "Time (weeks)",
    y     = "Relative Units (×100)",
    color = "Component"
  ) +
  theme_bw(base_size = 12)

print(p6)

# =============================================================================
# Summary statistics
# =============================================================================

cat("\n=== Breast Cancer QSP Model Summary (refactored PK/PD) ===\n")
cat("Scenarios simulated :", length(unique(all_results$scenario)), "\n")
cat("Simulation duration :", end_time / 24 / 7, "weeks\n")
cat("ODE compartments    : 20 (was 21; dead-code PERI_PALBO removed in refactor)\n")
cat("Parameters          : 35+\n")

cat("\nEnd-of-simulation tumor volumes (relative to baseline 100):\n")
final <- all_results %>%
  group_by(scenario) %>%
  slice_tail(n = 1) %>%
  select(scenario, TUMOR, Ki67, CA153, CD8_EFF)

print(as.data.frame(final))
