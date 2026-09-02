##############################################################################
# GERD QSP Model — mrgsolve implementation (PK/PD REFACTORED sibling)
# Disease: Gastroesophageal Reflux Disease (GERD)
# Model: 20-compartment ODE system
#   • PK: PPI (omeprazole default / esomeprazole scenario override),
#          H2RA (famotidine), P-CAB (vonoprazan), Prokinetic (domperidone)
#   • PD: H+/K+-ATPase turnover, gastric pH, esophageal acid exposure,
#          mucosal damage/healing, symptom score
# Calibration: Miner 2003 (omeprazole), Hunt 1984 (famotidine),
#              Ashida 2016 (vonoprazan VOYAGE), Kahrilas 2008 (esomeprazole)
# Scenarios: 6 treatment arms (unchanged from original)
# Author: Claude Code Routine (CCR) — 2026-06-18 (original), refactor 2026-09-02
#
# ── PK/PD REFACTOR — see gerd_refactor_notes.md for the full account ────────
# This sibling never edits gerd_mrgsolve_model.R; it is a fresh file that:
#   1. Isolates each of the four compounds' PK into its own GUT_<STEM>/
#      CENT_<STEM> compartment pair (naming convention from
#      FORK_WORKFLOW_GUIDE.md Part 2), renamed from PPI_GUT/PPI_CENT-style
#      to GUT_PPI/CENT_PPI-style, V_<STEM> -> V1_<STEM>.
#   2. Exposes exactly one concentration per compound, `C_<STEM>`, computed
#      once (in $MAIN, matching where the original computed its own
#      CP_<STEM> and consumed it for PD) and reused everywhere else —
#      never recomputed a second time from raw state.
#   3. Confirmed real duplicate-concentration-site defect, corpus-wide
#      classifier under-labelled it: the original independently computed
#      each of CP_PPI/CP_H2RA/CP_PCAB/CP_PROK TWICE under the SAME name
#      (once in $MAIN as a `double`, once in $TABLE as a `capture`
#      assignment) -- this is not just wasteful duplication, it is a
#      genuine mrgsolve 2.0.1 build defect (see notes: "redefinition of
#      capture ::CP_PPI" / "previously declared double ::CP_PPI"). census
#      only flagged PPI as "duplicate"; H2RA/P-CAB/Prokinetic have the
#      identical defect and are corrected here (see gerd_refactor_notes.md
#      and the census update).
#   4. `LES_pressure` in the original's own $TABLE was ALSO an independent
#      third re-derivation of the prokinetic LES effect (same formula as
#      $MAIN's PROK_EFF_LES, hardcoded EC50=0.005 baked in twice) --
#      consolidated to reference the single $MAIN-computed
#      EFFECT_PROK_LES.
#   5. Each compound's effect on disease is a named `EFFECT_<STEM>` Hill
#      term (rename only, gamma=1 where the original had no explicit Hill
#      coefficient; PPI's own explicit HILL_PPI is preserved as
#      GAMMA_PPI). Two prokinetic effects (LES pressure, gastric emptying)
#      keep separate EFFECT_PROK_LES / EFFECT_PROK_EMP terms sharing one
#      promoted EC50_PROK (was a hardcoded 0.005 literal in two places in
#      the original) -- same shape, no refit.
#   6. Two pre-existing upstream mrgsolve-2.0.1 build defects found in the
#      ORIGINAL file (unrelated to any single compound's own PK), fixed
#      syntax-only in THIS file only (never in the checked-in original),
#      logged as translations/UPSTREAM_ISSUES.md #<N> (see notes):
#        a. `$CMT @annotated` + a separate `$INIT` block jointly
#           redeclaring the same 17 compartment names ("Duplicated model
#           names" build failure). Fixed by dropping `$INIT` and setting
#           non-zero initial values via the `<cmt>_0 = value;` idiom in
#           $MAIN (matching this guide's own Archetype 4 example).
#        b. The $MAIN/$TABLE double-declaration of CP_PPI/CP_H2RA/
#           CP_PCAB/CP_PROK described in point 3 above.
#   7. `C_<STEM>`/`EFFECT_<STEM>` are $MAIN-local doubles, exposed via a
#      bare `$CAPTURE` list (not `$PARAM`) -- matching the corpus's
#      confirmed, disclosed mrgsolve 2.0.1 constraint that a `$PARAM`
#      member cannot also be reassigned each step in $MAIN/$ODE (same
#      precedent as acute-intermittent-porphyria/sepsis/amd).
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─── Model code ──────────────────────────────────────────────────────────────
gerd_refactored_code <- '
$PROB GERD QSP Model (refactored) — PPI/H2RA/P-CAB/Prokinetic PK-PD, pluggable PK

$PARAM @annotated
// ── Drug Doses (informational only -- actual dosing is via ev() events in
//    the R scenarios below; these were never read by $MAIN/$ODE/$TABLE in
//    the original either) ───────────────────────────────────────────────────
DOSE_PPI   : 20  : PPI dose (mg, omeprazole equivalent)
DOSE_H2RA  : 40  : H2RA dose (mg, famotidine equivalent)
DOSE_PCAB  : 20  : P-CAB dose (mg, vonoprazan)
DOSE_PROK  : 10  : Prokinetic dose (mg, domperidone)

// ── PPI PK Parameters (Omeprazole; esomeprazole scenario overrides F_PPI/CL_PPI) ──
KA_PPI     : 0.80  : PPI absorption rate constant (h-1)
CL_PPI     : 18.0  : PPI apparent clearance (L/h)
V1_PPI     : 20.0  : PPI central volume of distribution (L)
F_PPI      : 0.65  : PPI bioavailability (CYP2C19 EM)
CYP2C19    : 1.0   : CYP2C19 phenotype multiplier (EM=1, IM=0.6, PM=0.25)

// ── H2RA PK Parameters (Famotidine) ─────────────────────────────────────────
KA_H2RA    : 1.50  : H2RA absorption rate constant (h-1)
CL_H2RA    : 14.0  : H2RA clearance (L/h, renal dominant)
V1_H2RA    : 90.0  : H2RA central volume of distribution (L)
F_H2RA     : 0.50  : H2RA oral bioavailability

// ── P-CAB PK Parameters (Vonoprazan) ─────────────────────────────────────────
KA_PCAB    : 2.00  : P-CAB absorption rate constant (h-1)
CL_PCAB    : 8.0   : P-CAB clearance (L/h, CYP3A4)
V1_PCAB    : 300.0 : P-CAB central volume of distribution (L, lipophilic)
F_PCAB     : 0.60  : P-CAB oral bioavailability

// ── Prokinetic PK (Domperidone) ───────────────────────────────────────────────
KA_PROK    : 1.20  : Prokinetic absorption rate (h-1)
CL_PROK    : 60.0  : Prokinetic clearance (L/h)
V1_PROK    : 450.0 : Prokinetic central volume of distribution (L)
F_PROK     : 0.15  : Prokinetic bioavailability (first-pass)

// ── H+/K+-ATPase (Proton Pump) Dynamics ────────────────────────────────────
PUMP_TOTAL : 100.0 : Total H+/K+-ATPase pool (normalized units)
PUMP_ACT0  : 30.0  : Baseline active pump fraction (%)
K_ACT      : 0.30  : Pump activation rate (h-1)
K_DEACT    : 0.30  : Pump deactivation rate (h-1)
K_SYN_PUMP : 3.0   : Pump synthesis rate (units/h)
K_DEG_PUMP : 0.03  : Pump degradation rate (h-1)

// ── Hill interface: Emax/EC50/Hill-coefficient per compound (renamed from
//    the original IC50_*/HILL_*/EMAX_* -- values unchanged, this is a
//    rename, not a refit; GAMMA_* = 1 added explicit where the original had
//    no Hill coefficient of its own) ────────────────────────────────────────
EC50_PPI   : 0.15  : PPI EC50 for covalent pump binding (mg/L) -- was IC50_PPI
GAMMA_PPI  : 1.5   : PPI Hill coefficient -- was HILL_PPI
EMAX_PPI   : 1.0   : PPI maximal pump-binding fraction (new explicit; original had no scaling factor, i.e. implicit Emax=1)
EC50_PCAB  : 0.08  : P-CAB EC50 for pump inhibition (mg/L) -- was IC50_PCAB
GAMMA_PCAB : 1.0   : P-CAB Hill coefficient (new explicit; original was a plain ratio)
EMAX_PCAB  : 1.0   : P-CAB maximal pump-inhibition fraction (new explicit; original had no scaling factor)
EC50_H2RA  : 0.04  : H2RA EC50 for H2 receptor block (mg/L) -- was IC50_H2RA
GAMMA_H2RA : 1.0   : H2RA Hill coefficient (new explicit; original was a plain ratio)
EMAX_H2RA  : 0.60  : H2RA maximal inhibition (60% of stimulated acid) -- same name/value as original
EC50_PROK  : 0.005 : Prokinetic EC50, shared by both LES and emptying effects (new explicit; was a hardcoded 0.005 literal baked into two separate formulas in the original)
GAMMA_PROK : 1.0   : Prokinetic Hill coefficient (new explicit; original was a plain ratio)
EMAX_PROK_LES : 5.0  : Prokinetic max LES pressure increase (mmHg) -- was K_PROK_LES
EMAX_PROK_EMP : 0.30 : Prokinetic max gastric emptying rate increase (h-1) -- was K_PROK_EMP

// ── Gastric Acid Secretion ──────────────────────────────────────────────────
ACID_BASE  : 3.5   : Baseline gastric acid secretion rate (mmol/h)
ACID_MAX   : 12.0  : Maximal acid secretion rate (mmol/h, pentagastrin)
HIST_STIM  : 0.35  : Histamine contribution to acid secretion (weights H2RA into combined ACID_INH)
GAST_STIM  : 0.45  : Gastrin contribution to acid secretion
ACH_STIM   : 0.20  : Acetylcholine contribution
K_NEUT_INT : 0.5   : Intragastric acid neutralization rate (h-1)
PH_BUFF    : 1.5   : Buffer capacity constant
VOL_GAS    : 200.0 : Gastric volume (mL)

// ── LES & Reflux ─────────────────────────────────────────────────────────────
LES_P0     : 22.0  : Baseline LES pressure (mmHg)
TLESR_BASE : 6.0   : Baseline TLESR rate (events/h)
GASTRIC_EMP: 0.50  : Gastric emptying rate constant (h-1)
REFLUX_K   : 0.05  : Acid reflux to esophageal acid exposure conversion

// ── Esophageal Mucosal Dynamics ──────────────────────────────────────────────
MUC_HEAL0  : 0.30  : Baseline mucosal healing rate (units/day)
MUC_INJ_K  : 0.40  : Mucosal injury rate constant (per AET unit)
K_HEAL     : 0.20  : Mucosal healing rate constant (h-1)
GRADE_K    : 0.015 : Damage->LA grade conversion constant
MAX_DAMAGE : 100.0 : Maximum damage score (100 = LA Grade D)

// ── Symptom Score ────────────────────────────────────────────────────────────
SYM_K      : 0.05  : Symptom sensitivity (per AET unit)
SYM_HEAL_K : 0.30  : Symptom relief rate with acid suppression

$CMT @annotated
// PPI PK
GUT_PPI    : PPI gut compartment (mg) -- was PPI_GUT
CENT_PPI   : PPI central compartment (mg) -- was PPI_CENT

// H2RA PK
GUT_H2RA   : H2RA gut compartment (mg) -- was H2RA_GUT
CENT_H2RA  : H2RA central compartment (mg) -- was H2RA_CENT

// P-CAB PK
GUT_PCAB   : P-CAB gut compartment (mg) -- was PCAB_GUT
CENT_PCAB  : P-CAB central compartment (mg) -- was PCAB_CENT

// Prokinetic PK
GUT_PROK   : Prokinetic gut compartment (mg) -- was PROK_GUT
CENT_PROK  : Prokinetic central compartment (mg) -- was PROK_CENT

// H+/K+-ATPase Pool
PUMP_INACT : Inactive pump pool (tubulovesicular)
PUMP_ACT   : Active pump pool (secretory canaliculus)
PUMP_INH   : Irreversibly inhibited pump pool (PPIs only)

// Gastric Acid & pH
ACID_RATE  : Gastric acid secretion rate (mmol/h)
GAS_pH     : Intragastric pH (transformed)

// Esophageal Compartments
AET        : Acid exposure time (% time pH<4)
MUC_DMG    : Mucosal damage score (0-100)
MUC_HEAL   : Mucosal integrity index (0-100)

// Symptom
SYM_SCORE  : Symptom score (GERD-Q equivalent, 0-18)

// Barrett progression (long-term)
BE_RISK    : Cumulative Barrett risk index (0-1)

$MAIN
// Upstream build-defect fix (a): the original paired $CMT @annotated with a
// separate $INIT block re-declaring the same 17 names ("Duplicated model
// names" -- confirmed live against qspserver, see gerd_refactor_notes.md).
// Non-zero initial values are set here instead, via mrgsolve’s own
// <cmt>_0 idiom (this guide’s own Archetype 4 example) -- no numeric change,
// PPI/H2RA/P-CAB/Prokinetic compartments still start at 0 same as original.
PUMP_INACT_0 = 70;
PUMP_ACT_0   = 30;
PUMP_INH_0   = 0;
ACID_RATE_0  = 3.5;
GAS_pH_0     = 1.8;
AET_0        = 15.0;
MUC_DMG_0    = 25.0;
MUC_HEAL_0   = 75.0;
SYM_SCORE_0  = 8.0;
BE_RISK_0    = 0.0;

// ── Exposed concentrations (single canonical site per compound; computed
//    once here, matching the block the original used for its own PD read --
//    never recomputed again in $TABLE, see upstream build-defect fix (b)
//    in the header notes) ─────────────────────────────────────────────────
double C_PPI  = CENT_PPI  / V1_PPI;
double C_H2RA = CENT_H2RA / V1_H2RA;
double C_PCAB = CENT_PCAB / V1_PCAB;
double C_PROK = CENT_PROK / V1_PROK;

// PPI: Emax covalent inhibition of active pumps (rename of INH_PPI, same shape)
double EFFECT_PPI = EMAX_PPI * pow(C_PPI, GAMMA_PPI)
                    / (pow(EC50_PPI, GAMMA_PPI) + pow(C_PPI, GAMMA_PPI));

// P-CAB: reversible ionic block (rename of INH_PCAB, same shape, gamma=1)
double EFFECT_PCAB = EMAX_PCAB * pow(C_PCAB, GAMMA_PCAB)
                     / (pow(EC50_PCAB, GAMMA_PCAB) + pow(C_PCAB, GAMMA_PCAB));

// H2RA: competitive antagonism at H2 receptor (rename of INH_H2RA, same shape, gamma=1)
double EFFECT_H2RA = EMAX_H2RA * pow(C_H2RA, GAMMA_H2RA)
                     / (pow(EC50_H2RA, GAMMA_H2RA) + pow(C_H2RA, GAMMA_H2RA));

// Combined acid suppression (PPI+PCAB irreversible/reversible; H2RA partial) --
// combined only at this point of use, each compound’s own EFFECT_<STEM> kept
// separate above, per the guide’s multi-drug-one-pathway rule.
double ACID_INH = 1.0 - (1.0 - EFFECT_PPI) * (1.0 - EFFECT_PCAB) * (1.0 - EFFECT_H2RA * HIST_STIM);
if(ACID_INH > 0.99) ACID_INH = 0.99;

// Prokinetic: two independent named effects sharing one EC50_PROK (was a
// 0.005 literal hardcoded separately in $MAIN and again in $TABLE’s
// LES_pressure re-derivation -- promoted to one explicit shared parameter,
// same value, same shape, no refit).
double EFFECT_PROK_LES = EMAX_PROK_LES * pow(C_PROK, GAMMA_PROK)
                          / (pow(EC50_PROK, GAMMA_PROK) + pow(C_PROK, GAMMA_PROK));
double EFFECT_PROK_EMP = EMAX_PROK_EMP * pow(C_PROK, GAMMA_PROK)
                          / (pow(EC50_PROK, GAMMA_PROK) + pow(C_PROK, GAMMA_PROK));

double LES_P = LES_P0 + EFFECT_PROK_LES;

// TLESR rate — inversely proportional to LES pressure
double TLESR = TLESR_BASE * (LES_P0 / LES_P);

// Prokinetic effect on gastric emptying
double EMP_RATE = GASTRIC_EMP + EFFECT_PROK_EMP;

// Acid secretion (effective)
double ACID_EFF = ACID_BASE + (ACID_MAX - ACID_BASE) * (1.0 - ACID_INH);

// pH calculation from acid secretion rate
double pH_calc = 1.0 + PH_BUFF * exp(-ACID_EFF / 2.0);
if(pH_calc < 1.0) pH_calc = 1.0;
if(pH_calc > 7.0) pH_calc = 7.0;

// Reflux acid volume → AET
double REFLUX_RATE = TLESR * REFLUX_K * (pH_calc < 4.0 ? 1.0 : exp(-pH_calc + 4.0));
double AET_INST = REFLUX_RATE * 100.0;  // % time pH<4

// Mucosal damage dynamics
double DMG_RATE = MUC_INJ_K * (AET / 100.0) * (1.0 - MUC_HEAL / MAX_DAMAGE);
double HEAL_RATE = K_HEAL * (MUC_HEAL0 / 24.0) * (1.0 + ACID_INH);

$ODE
// ── PPI PK ──────────────────────────────────────────────────────────────────
dxdt_GUT_PPI  = -KA_PPI * GUT_PPI;
dxdt_CENT_PPI = KA_PPI * F_PPI * GUT_PPI * CYP2C19
                - (CL_PPI / V1_PPI) / CYP2C19 * CENT_PPI;

// ── H2RA PK ─────────────────────────────────────────────────────────────────
dxdt_GUT_H2RA  = -KA_H2RA * GUT_H2RA;
dxdt_CENT_H2RA = KA_H2RA * F_H2RA * GUT_H2RA - (CL_H2RA / V1_H2RA) * CENT_H2RA;

// ── P-CAB PK ─────────────────────────────────────────────────────────────────
dxdt_GUT_PCAB  = -KA_PCAB * GUT_PCAB;
dxdt_CENT_PCAB = KA_PCAB * F_PCAB * GUT_PCAB - (CL_PCAB / V1_PCAB) * CENT_PCAB;

// ── Prokinetic PK ─────────────────────────────────────────────────────────────
dxdt_GUT_PROK  = -KA_PROK * GUT_PROK;
dxdt_CENT_PROK = KA_PROK * F_PROK * GUT_PROK - (CL_PROK / V1_PROK) * CENT_PROK;

// ── Proton Pump Pool Dynamics ────────────────────────────────────────────────
// Synthesis into inactive pool; translocation to active; irreversible PPI binding
dxdt_PUMP_INACT = K_SYN_PUMP - K_DEG_PUMP * PUMP_INACT - K_ACT * PUMP_INACT;
dxdt_PUMP_ACT   = K_ACT * PUMP_INACT - K_DEACT * PUMP_ACT
                  - K_DEG_PUMP * PUMP_ACT
                  - EFFECT_PPI * 0.5 * PUMP_ACT;  // covalent PPI binding
dxdt_PUMP_INH   = EFFECT_PPI * 0.5 * PUMP_ACT - K_DEG_PUMP * PUMP_INH; // degraded at normal rate

// ── Gastric Acid & pH ────────────────────────────────────────────────────────
double PUMP_FRAC = PUMP_ACT / PUMP_ACT0;
dxdt_ACID_RATE  = 0.1 * (ACID_EFF * PUMP_FRAC - ACID_RATE);  // half-life ~7h
dxdt_GAS_pH     = 0.5 * (pH_calc - GAS_pH);                   // slow equilibration

// ── Esophageal Acid Exposure Time ────────────────────────────────────────────
// AET drifts toward instantaneous AET with clearance
dxdt_AET = 0.2 * (AET_INST - AET);

// ── Mucosal Damage / Healing ─────────────────────────────────────────────────
dxdt_MUC_DMG  = DMG_RATE - HEAL_RATE;
if(MUC_DMG < 0) dxdt_MUC_DMG = 0;
if(MUC_DMG > MAX_DAMAGE) dxdt_MUC_DMG = -0.1;

dxdt_MUC_HEAL = -DMG_RATE + HEAL_RATE;
if(MUC_HEAL < 0) dxdt_MUC_HEAL = 0;
if(MUC_HEAL > MAX_DAMAGE) dxdt_MUC_HEAL = 0;

// ── Symptom Score ─────────────────────────────────────────────────────────────
double SYM_TARGET = 15.0 * (AET / 30.0) + 3.0 * (MUC_DMG / MAX_DAMAGE);
dxdt_SYM_SCORE = 0.1 * (SYM_TARGET - SYM_SCORE);

// ── Barrett Risk (long-term cumulative) ────────────────────────────────────────
dxdt_BE_RISK = 0.0001 * AET * (1.0 - BE_RISK);

$TABLE
// Disease-side outputs unchanged from the original, except LES_pressure --
// the original independently re-derived the prokinetic LES effect a THIRD
// time here (same formula as $MAIN’s PROK_EFF_LES, with the EC50 0.005
// literal baked in yet again); redirected to reference the single
// $MAIN-computed EFFECT_PROK_LES instead (see header notes point 4).
capture ACID_INH_pct = ACID_INH * 100;
capture pH       = GAS_pH;
capture AET_pct  = AET;
capture DMG      = MUC_DMG;
capture HEAL     = MUC_HEAL;
capture SYM      = SYM_SCORE;
capture PUMP_active = PUMP_ACT;
capture PUMP_inhibited = PUMP_INH;
capture LES_pressure = LES_P0 + EFFECT_PROK_LES;
capture Barrett  = BE_RISK;

$CAPTURE @annotated
// Canonical single-site concentrations and Hill effects per compound --
// $MAIN-local doubles + bare $CAPTURE (not $PARAM): declaring these as
// $PARAM would make them read-only in $MAIN, conflicting with recomputing
// them every step (same confirmed mrgsolve 2.0.1 constraint documented for
// acute-intermittent-porphyria/sepsis/amd in this corpus).
C_PPI  : PPI plasma concentration (mg/L) -- single canonical site, was CP_PPI (duplicated $MAIN/$TABLE)
C_H2RA : H2RA plasma concentration (mg/L) -- single canonical site, was CP_H2RA (duplicated $MAIN/$TABLE)
C_PCAB : P-CAB plasma concentration (mg/L) -- single canonical site, was CP_PCAB (duplicated $MAIN/$TABLE)
C_PROK : Prokinetic plasma concentration (mg/L) -- single canonical site, was CP_PROK (duplicated $MAIN/$TABLE)
EFFECT_PPI : PPI covalent pump-inhibition fraction (0-1) -- was INH_PPI
EFFECT_PCAB : P-CAB reversible pump-inhibition fraction (0-1) -- was INH_PCAB
EFFECT_H2RA : H2RA H2-receptor-block fraction (0-1) -- was INH_H2RA
EFFECT_PROK_LES : Prokinetic LES-pressure-increase effect (mmHg) -- was PROK_EFF_LES
EFFECT_PROK_EMP : Prokinetic gastric-emptying-rate-increase effect (h-1) -- new name, was inline in EMP_RATE
'

mod <- mcode("GERD_QSP_REFACTORED", gerd_refactored_code)

# ─── Treatment Scenarios (unchanged from the original; compartment names
#     updated to the refactored GUT_<STEM> convention) ──────────────────────
# Dosing events for 8 weeks (56 days, 1344 h), QD morning dosing

make_events <- function(drug, dose, interval = 24, duration = 56*24) {
  times <- seq(0, duration - 1, by = interval)
  ev(amt = dose, cmt = drug, time = times)
}

scenarios <- list(
  "No Treatment (Control)" = list(
    mrgsolve::ev(amt = 0, cmt = "GUT_PPI", time = 0),
    param = list(DOSE_PPI = 0, DOSE_H2RA = 0, DOSE_PCAB = 0, DOSE_PROK = 0)
  ),
  "Omeprazole 20 mg QD (Standard PPI)" = list(
    ev = make_events("GUT_PPI", 20),
    param = list(DOSE_PPI = 20)
  ),
  "Esomeprazole 40 mg QD (High-dose PPI)" = list(
    ev = make_events("GUT_PPI", 40),
    param = list(DOSE_PPI = 40, F_PPI = 0.73, CL_PPI = 12)
  ),
  "Vonoprazan 20 mg QD (P-CAB)" = list(
    ev = make_events("GUT_PCAB", 20),
    param = list(DOSE_PCAB = 20)
  ),
  "Famotidine 40 mg BID (H2RA)" = list(
    ev = c(make_events("GUT_H2RA", 40, interval = 12)),
    param = list(DOSE_H2RA = 40)
  ),
  "Eso 40 mg QD + Domperidone 10 mg TID\n(PPI + Prokinetic)" = list(
    ev = c(make_events("GUT_PPI", 40),
           make_events("GUT_PROK", 10, interval = 8)),
    param = list(DOSE_PPI = 40, F_PPI = 0.73, CL_PPI = 12, DOSE_PROK = 10)
  )
)

# ─── Simulate all scenarios ──────────────────────────────────────────────────
sim_list <- lapply(seq_along(scenarios), function(i) {
  scen <- scenarios[[i]]
  nm   <- names(scenarios)[i]

  m2 <- mod
  if (!is.null(scen$param)) m2 <- param(m2, scen$param)

  ev_obj <- if (inherits(scen[[1]], "ev")) scen[[1]] else scen$ev

  out <- mrgsim(m2, ev_obj,
                start = 0, end = 56 * 24, delta = 0.5,
                carry_out = "evid")
  out <- as.data.frame(out)
  out$Scenario <- nm
  out
})

all_sims <- bind_rows(sim_list)

# ─── Summary Statistics (Week 8 endpoint) ───────────────────────────────────
summary_stats <- all_sims %>%
  filter(time == max(time)) %>%
  group_by(Scenario) %>%
  summarise(
    pH_mean         = round(mean(pH, na.rm = TRUE), 2),
    AET_pct_mean    = round(mean(AET_pct, na.rm = TRUE), 1),
    Acid_inh_pct    = round(mean(ACID_INH_pct, na.rm = TRUE), 1),
    Mucosal_damage  = round(mean(DMG, na.rm = TRUE), 1),
    Symptom_score   = round(mean(SYM, na.rm = TRUE), 1),
    Barrett_risk    = round(mean(Barrett, na.rm = TRUE), 4),
    .groups = "drop"
  )

print(summary_stats)

# ─── Plots ───────────────────────────────────────────────────────────────────
cols <- c(
  "No Treatment (Control)"                      = "#E53935",
  "Omeprazole 20 mg QD (Standard PPI)"          = "#FB8C00",
  "Esomeprazole 40 mg QD (High-dose PPI)"       = "#8E24AA",
  "Vonoprazan 20 mg QD (P-CAB)"                 = "#00897B",
  "Famotidine 40 mg BID (H2RA)"                 = "#1E88E5",
  "Eso 40 mg QD + Domperidone 10 mg TID\n(PPI + Prokinetic)" = "#43A047"
)

# Thin data for plotting
plot_data <- all_sims %>% filter(time %% 6 == 0)

p1 <- ggplot(plot_data, aes(x = time / 24, y = pH, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "A. Intragastric pH", x = "Time (days)", y = "pH") +
  geom_hline(yintercept = 4, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 4.2, label = "pH 4 threshold", size = 3) +
  theme_classic() + theme(legend.position = "none")

p2 <- ggplot(plot_data, aes(x = time / 24, y = AET_pct, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "B. Acid Exposure Time (%)", x = "Time (days)", y = "AET (%)") +
  geom_hline(yintercept = 6, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 7, label = "Lyon 2.0 cutoff (6%)", size = 3) +
  theme_classic() + theme(legend.position = "none")

p3 <- ggplot(plot_data, aes(x = time / 24, y = DMG, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "C. Mucosal Damage Score", x = "Time (days)", y = "Damage (0-100)") +
  theme_classic() + theme(legend.position = "none")

p4 <- ggplot(plot_data, aes(x = time / 24, y = SYM, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "D. Symptom Score (GERD-Q)", x = "Time (days)", y = "Score (0-18)") +
  geom_hline(yintercept = 8, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 8.5, label = "GERD-Q ≥8 = GERD", size = 3) +
  theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 7))

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title    = "GERD QSP Model — 8-week Treatment Simulation (refactored)",
    subtitle = "PPI vs H2RA vs P-CAB vs Prokinetic combination",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined)

# ─── CYP2C19 Phenotype Sensitivity Analysis ──────────────────────────────────
cyp_scenarios <- list(
  "Ultra-Rapid (UM, CYP2C19×2.0)"  = list(CYP2C19 = 2.0),
  "Extensive (EM, CYP2C19×1.0)"    = list(CYP2C19 = 1.0),
  "Intermediate (IM, CYP2C19×0.6)" = list(CYP2C19 = 0.6),
  "Poor (PM, CYP2C19×0.25)"        = list(CYP2C19 = 0.25)
)

cyp_sims <- lapply(seq_along(cyp_scenarios), function(i) {
  m2 <- param(mod, cyp_scenarios[[i]])
  out <- mrgsim(m2, make_events("GUT_PPI", 20),
                start = 0, end = 56 * 24, delta = 0.5)
  out <- as.data.frame(out)
  out$Phenotype <- names(cyp_scenarios)[i]
  out
})

cyp_data <- bind_rows(cyp_sims) %>% filter(time %% 6 == 0)

p_cyp <- ggplot(cyp_data, aes(x = time / 24, y = pH, colour = Phenotype)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "CYP2C19 Phenotype Effect on Gastric pH\n(Omeprazole 20 mg QD)",
    x        = "Time (days)",
    y        = "Intragastric pH",
    colour   = "CYP2C19 Phenotype"
  ) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_colour_manual(values = c("#C62828", "#1976D2", "#388E3C", "#F57F17")) +
  theme_classic()

print(p_cyp)

# ─── Dose-Response at Week 8 ─────────────────────────────────────────────────
doses <- c(5, 10, 20, 40, 80)
dr_ppi <- lapply(doses, function(d) {
  out <- mrgsim(mod, make_events("GUT_PPI", d),
                start = 0, end = 56 * 24, delta = 1)
  out <- as.data.frame(out)
  data.frame(Dose = d, AET = tail(out$AET_pct, 1), pH = tail(out$pH, 1),
             Drug = "PPI (Omeprazole)")
})

dr_pcab <- lapply(doses, function(d) {
  out <- mrgsim(mod, make_events("GUT_PCAB", d),
                start = 0, end = 56 * 24, delta = 1)
  out <- as.data.frame(out)
  data.frame(Dose = d, AET = tail(out$AET_pct, 1), pH = tail(out$pH, 1),
             Drug = "P-CAB (Vonoprazan)")
})

dr_data <- bind_rows(c(dr_ppi, dr_pcab))

p_dr <- ggplot(dr_data, aes(x = Dose, y = AET, colour = Drug, group = Drug)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_x_log10() +
  labs(
    title  = "Dose-Response: AET (%) at Week 8",
    x      = "Dose (mg, log scale)",
    y      = "Acid Exposure Time (%)",
    colour = "Drug Class"
  ) +
  geom_hline(yintercept = 6, linetype = "dashed") +
  scale_colour_manual(values = c("#8E24AA", "#00897B")) +
  theme_classic()

print(p_dr)

message("\n=== GERD QSP Model (refactored) simulation complete ===")
message("Key Week-8 endpoints:")
print(summary_stats)
