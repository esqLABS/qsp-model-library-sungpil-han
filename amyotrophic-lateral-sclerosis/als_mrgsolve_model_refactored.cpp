$PROB
ALS QSP Model — REFACTORED (pluggable PK)
Amyotrophic Lateral Sclerosis — Motor Neuron Degeneration, Protein Aggregation,
Excitotoxicity, Neuroinflammation, Drug PK/PD

Refactor sibling of als_mrgsolve_model.R. Four compounds isolated to the
naming convention (RIL/EDA/TOF/PB): Riluzole, Edaravone, Tofersen,
Phenylbutyrate/AMX0035. See als_refactor_notes.md for the full account,
especially the Tofersen CSF-vs-plasma finding (C_TOF exposes the CSF
concentration, not plasma) and the two pre-existing mrgsolve 2.0.1 build
defects fixed syntax-only in this sibling (never in the original).

$PARAM @annotated
// ── Motor Neuron Parameters ──────────────────────────────────────
k_MN_death    : 0.00080  : baseline MN death rate constant (1/day)
k_MN_death_U  : 0.00060  : upper MN death rate constant  (1/day)
MN_upper_REF  : 1.0      : baseline upper motor neuron reference (normalized) (was MN_upper_0)
MN_lower_REF  : 1.0      : baseline lower motor neuron reference (normalized) (was MN_lower_0)

// ── SOD1 Protein Dynamics ────────────────────────────────────────
k_SOD1_syn    : 0.50     : SOD1 synthesis rate (AU/day)
k_SOD1_deg    : 0.10     : SOD1 wild-type degradation rate (1/day)
k_SOD1_mis    : 0.08     : SOD1 misfolding rate constant (1/day)
k_SOD1_clr    : 0.015    : misfolded SOD1 clearance (UPS/autophagy, 1/day)
f_SOD1_mut    : 0.0      : fraction mutant (0=sporadic, 1=SOD1-ALS)

// ── TDP-43 Dynamics ──────────────────────────────────────────────
k_TDP_export  : 0.020    : TDP-43 nuclear→cytoplasmic export (1/day)
k_TDP_import  : 0.150    : TDP-43 cytoplasmic→nuclear import (1/day)
k_TDP_agg     : 0.010    : TDP-43 cytoplasmic aggregation (1/day)
k_TDP_clr     : 0.008    : TDP-43 aggregate clearance (1/day)
TDP43_nuc_REF : 10.0     : baseline nuclear TDP-43 reference (AU) (was TDP43_nuc_0)

// ── Glutamate Excitotoxicity ──────────────────────────────────────
Glu_base      : 1.0      : baseline synaptic glutamate (AU)
k_Glu_rel     : 2.0      : glutamate release rate (AU/day)
k_EAAT2       : 1.80     : EAAT2 baseline uptake constant (1/day)
k_EAAT2_ALS   : 0.40     : EAAT2 ALS reduction factor (dimensionless <1)
k_Ca_entry    : 0.50     : Ca2+ influx rate constant per Glu excess
k_Ca_efflux   : 1.20     : Ca2+ efflux/buffering rate (1/day)
Ca_i_REF      : 0.10     : baseline intracellular Ca2+ reference (μM) (was Ca_i_0)

// ── Oxidative Stress ─────────────────────────────────────────────
k_ROS_prod    : 0.40     : ROS production rate per Ca2+ & SOD1mis (AU/day)
k_ROS_scav    : 0.60     : ROS scavenging by GSH (1/day per AU_GSH)
k_GSH_syn     : 0.50     : GSH synthesis rate (AU/day)
k_GSH_cons    : 0.12     : GSH consumption (1/day, basal turnover)
ROS_REF       : 0.50     : baseline ROS reference (AU) (was ROS_0)
GSH_REF       : 4.0      : baseline GSH reference (AU) (was GSH_0)

// ── Mitochondrial Function ───────────────────────────────────────
k_Mito_dam    : 0.06     : mitochondrial damage by ROS (1/day per AU_ROS)
k_Mito_rep    : 0.03     : mitochondrial repair/mitophagy (1/day)
Mito_REF      : 1.0      : baseline mitochondrial function reference (normalized) (was Mito_0)

// ── Neuroinflammation ────────────────────────────────────────────
k_Mic_act     : 0.12     : microglial activation rate (1/day)
k_Mic_res     : 0.060    : microglial resolution rate (1/day)
k_TNFa_prod   : 0.60     : TNF-α production by M1 microglia (AU/day)
k_TNFa_deg    : 0.35     : TNF-α degradation rate (1/day)
k_IL1b_prod   : 0.40     : IL-1β production rate (AU/day)
k_IL1b_deg    : 0.25     : IL-1β degradation (1/day)
Mic_0         : 0.10     : baseline microglial activation (0–1)

// ── Trophic Support ──────────────────────────────────────────────
k_BDNF_prod   : 0.25     : BDNF production (AU/day)
k_BDNF_deg    : 0.18     : BDNF degradation (1/day)
BDNF_REF      : 1.0      : baseline BDNF reference (AU) (was BDNF_0)
BDNF_prot_wt  : 0.30     : max neuroprotection weight from BDNF (dimensionless)

// ── Biomarkers ───────────────────────────────────────────────────
k_NfL_rel     : 0.10     : NfL release rate per MN death (AU/day)
k_NfL_clr     : 0.025    : NfL clearance from CSF (1/day)
NfL_CSF_REF   : 5.0      : baseline CSF NfL reference (pg/mL AU) (was NfL_CSF_0)

// ── Clinical Endpoints ───────────────────────────────────────────
ALSFRS_REF    : 48.0     : initial ALSFRS-R total score reference (was ALSFRS_0)
FVC_REF       : 100.0    : initial FVC % predicted reference (was FVC_0)
k_FVC_dec     : 0.0060   : FVC decline rate per unit TNFa & Mic (1/day)

// ── Riluzole PK (archetype 3: depot+central+peripheral) ──────────
F_RIL         : 0.60     : riluzole oral bioavailability
KA_RIL        : 0.80     : absorption rate (1/h) (was ka_RIL)
CL_RIL        : 28.0     : clearance (L/h)
V1_RIL        : 245.0    : central volume (L)
Q_RIL         : 15.0     : inter-compartmental CL (L/h)
V2_RIL        : 112.0    : peripheral volume (L)
// Riluzole Hill interface (rename of original’s plain C/(C+IC50) ratio)
EC50_RIL      : 0.50     : conc. for half-max glutamate-release inhibition (μg/mL) (was IC50_RIL)
EMAX_RIL      : 0.60     : max inhibition of Glu release (fraction) (was Emax_RIL)
GAMMA_RIL     : 1.0      : Hill coefficient (none in original; =1)

// ── Edaravone PK (archetype 1: single compartment, linear) ───────
CL_EDA        : 18.0     : edaravone clearance (L/h)
V1_EDA        : 120.0    : edaravone volume (L) (was V_EDA)
// Edaravone Hill interface
EC50_EDA      : 1.20     : conc. for half-max ROS scavenging (μg/mL) (was IC50_EDA)
EMAX_EDA      : 0.70     : max ROS reduction (fraction) (was Emax_EDA)
GAMMA_EDA     : 1.0      : Hill coefficient (none in original; =1)

// ── Tofersen PK (archetype 3: depot+central(plasma)+peripheral(CSF)) ──
KA_TOF        : 0.030    : tofersen SC absorption rate (1/h) (was ka_TOF)
CL_TOF        : 0.50     : tofersen plasma clearance (L/h)
V1_TOF        : 15.0     : central (plasma) volume (L)
Q_TOF         : 0.30     : distribution to CSF (L/h)
V2_TOF        : 5.0      : CSF volume (L)
// Tofersen Hill interface — reads CSF concentration (C_TOF), NOT plasma;
// this is already a plain Emax/EC50 ratio in the original despite Tofersen
// being an antisense oligonucleotide (see refactor notes) — a rename, not
// a mechanistic-kinetics fit.
EC50_TOF      : 0.10     : EC50 for SOD1 mRNA knockdown (μg/mL CSF)
EMAX_TOF      : 0.80     : max SOD1 mRNA knockdown (fraction) (was Emax_TOF)
GAMMA_TOF     : 1.0      : Hill coefficient (none in original; =1)

// ── Phenylbutyrate PK (archetype 3 minus peripheral: depot+central) ──
F_PB          : 0.85     : phenylbutyrate bioavailability
KA_PB         : 1.20     : PB absorption rate (1/h) (was ka_PB)
CL_PB         : 12.0     : PB clearance (L/h)
V1_PB         : 50.0     : PB volume (L) (was V_PB)
// PB Hill interface — two effects share the same EC50/gamma (as the
// original does: one EC50_PB reused for both), distinct EMAX each.
EC50_PB       : 50.0     : EC50 ER-stress/mitochondrial effect (μmol/L PB)
GAMMA_PB      : 1.0      : Hill coefficient (none in original; =1)
EMAX_PB       : 0.50     : max ER stress / CHOP reduction (was Emax_PB) — computed but NOT read by any $ODE below, same as the original (see refactor notes: pre-existing dead effect term)
EMAX_PB_MITO  : 0.30     : max mitochondrial protection by PB+TUDCA (was Emax_mito_PB)

$CMT @annotated
// ── Drug PK ──────────────────────────────────────────────────────
GUT_RIL    : Riluzole oral depot (mg)
CENT_RIL   : Riluzole central compartment (mg)
PERI_RIL   : Riluzole peripheral compartment (mg)
CENT_EDA   : Edaravone central compartment (mg)
GUT_TOF    : Tofersen SC depot (mg)
CENT_TOF   : Tofersen plasma compartment (mg)
PERI_TOF   : Tofersen CSF compartment (mg) — the PD-facing tissue site
GUT_PB     : Phenylbutyrate oral depot (mg)
CENT_PB    : Phenylbutyrate plasma compartment (mg)

// ── Disease Biology ───────────────────────────────────────────────
MN_upper   : Upper motor neurons (corticospinal, normalized)
MN_lower   : Lower motor neurons (spinal/bulbar, normalized)
SOD1_wt    : Wild-type SOD1 protein (AU)
SOD1_mis   : Misfolded SOD1 aggregates (AU)
TDP43_nuc  : Nuclear TDP-43 (AU)
TDP43_cyto : Cytoplasmic TDP-43 (AU)
Glut_syn   : Synaptic glutamate (AU)
Ca_i       : Intracellular calcium (μM)
ROS        : Reactive oxygen species (AU)
GSH        : Glutathione (AU)
Mito       : Mitochondrial integrity (normalized 0–1)
Mic_act    : Microglial activation state (0–1)
TNFa       : TNF-α (AU)
BDNF       : BDNF trophic factor (AU)

// ── Biomarkers & Clinical ─────────────────────────────────────────
NfL_CSF    : CSF neurofilament light chain (pg/mL AU)
ALSFRS     : ALSFRS-R total score (0–48)
FVC        : FVC % predicted

$MAIN
// ── Initial conditions (nonzero only; zero is mrgsolve’s own default for
// the PK compartments) — set via the `<CMT>_0 = value;` idiom rather than
// a $CMT-annotation value column or a separate $INIT block. Both of those
// alternatives were tried and rejected: a $CMT value column collides with
// mrgsolve 2.0.1’s own bare-vs-annotated duplicate-declaration check when
// paired with $INIT (see refactor notes, defect 1), and — found only once
// that was fixed — the qspserver mrgsolve_api runner does not honor a
// $CMT-annotation value column at all (confirmed empirically: a compartment
// declared `NAME : 48.0 : desc` still reported 0 at t=0 through this API;
// only an explicit `NAME_0 = 48.0;` assignment in $MAIN is honored).
MN_upper_0   = 1.0;
MN_lower_0   = 1.0;
SOD1_wt_0    = 5.0;
TDP43_nuc_0  = 10.0;
TDP43_cyto_0 = 0.5;
Glut_syn_0   = 1.0;
Ca_i_0       = 0.1;
ROS_0        = 0.5;
GSH_0        = 4.0;
Mito_0       = 1.0;
Mic_act_0    = 0.1;
TNFa_0       = 0.2;
BDNF_0       = 1.0;
NfL_CSF_0    = 5.0;
ALSFRS_0     = 48.0;
FVC_0        = 100.0;

// ── Exposed, PD-facing concentrations (one per compound) ─────────
double C_RIL  = CENT_RIL / V1_RIL;         // riluzole μg/mL
double C_EDA  = CENT_EDA / V1_EDA;         // edaravone μg/mL
double C_TOF  = PERI_TOF / V2_TOF;         // tofersen CSF μg/mL — PD-facing site
double Cp_TOF_plasma = CENT_TOF / V1_TOF;  // tofersen plasma μg/mL — informational only, NOT read by any PD term
double C_PB   = (CENT_PB / V1_PB) * 1000.0; // PB μmol/L (MW≈122)

// ── Hill-interface drug effects ───────────────────────────────────
double EFFECT_RIL      = EMAX_RIL      * pow(C_RIL, GAMMA_RIL)     / (pow(EC50_RIL, GAMMA_RIL)     + pow(C_RIL, GAMMA_RIL)     + 1e-9);
double EFFECT_EDA      = EMAX_EDA      * pow(C_EDA, GAMMA_EDA)     / (pow(EC50_EDA, GAMMA_EDA)     + pow(C_EDA, GAMMA_EDA)     + 1e-9);
double EFFECT_TOF      = EMAX_TOF      * pow(C_TOF, GAMMA_TOF)     / (pow(EC50_TOF, GAMMA_TOF)     + pow(C_TOF, GAMMA_TOF)     + 1e-9);
double EFFECT_PB       = EMAX_PB       * pow(C_PB,  GAMMA_PB)      / (pow(EC50_PB,  GAMMA_PB)      + pow(C_PB,  GAMMA_PB)      + 1e-9); // dead: not read below, matches original
double EFFECT_PB_MITO  = EMAX_PB_MITO  * pow(C_PB,  GAMMA_PB)      / (pow(EC50_PB,  GAMMA_PB)      + pow(C_PB,  GAMMA_PB)      + 1e-9);

// ── SOD1 aggregation burden ──────────────────────────────────────
double total_SOD1  = SOD1_wt + SOD1_mis + 1e-9;
double frac_SOD1mis = SOD1_mis / total_SOD1;
double SOD1_burden  = f_SOD1_mut * frac_SOD1mis;  // 0 in sporadic ALS

// ── TDP-43 pathology ─────────────────────────────────────────────
double frac_TDP_cyto = TDP43_cyto / (TDP43_nuc + TDP43_cyto + 1e-9);

// ── Excitotoxicity normalized ────────────────────────────────────
double Glu_excess = fmax(0.0, Glut_syn - Glu_base);

// ── ROS & inflammation burden ────────────────────────────────────
double ROS_norm  = ROS / (ROS_REF + 1e-9);
double TNFa_norm = TNFa / 1.0;

// ── Net motor neuron death rate (multi-hit) ──────────────────────
double MN_death_rate = k_MN_death * (
    1.0
  + 3.0 * SOD1_burden         // SOD1 tox
  + 2.0 * frac_TDP_cyto       // TDP-43 pathology
  + 1.5 * ROS_norm             // oxidative stress
  + 1.2 * Glu_excess           // excitotoxicity
  + 0.8 * TNFa_norm            // neuro-inflammation
  - BDNF_prot_wt * (BDNF / BDNF_REF)  // trophic protection
);
MN_death_rate = fmax(0.0001, MN_death_rate); // floor to avoid negatives

// ── EAAT2 reduction in ALS (inflammation-driven) ─────────────────
double EAAT2_eff = k_EAAT2 * (k_EAAT2_ALS + (1.0 - k_EAAT2_ALS) / (1.0 + Mic_act));

// ── SOD1 synthesis (reduced by tofersen) ─────────────────────────
double SOD1_syn_eff = k_SOD1_syn * (1.0 - EFFECT_TOF);

// ── NfL release proportional to MN death ─────────────────────────
double NfL_release = k_NfL_rel * MN_death_rate * (MN_upper + MN_lower);

// ── Microglial activation stimulus ───────────────────────────────
double Mic_stim = 0.5 * SOD1_burden
                + 0.3 * frac_TDP_cyto
                + 0.3 * ROS_norm
                + 0.5 * (MN_upper_REF - MN_upper)
                + 0.5 * (MN_lower_REF - MN_lower);

$ODE
// ─── DRUG PK ODEs ────────────────────────────────────────────────
// Riluzole — archetype 3 (depot+central+peripheral)
dxdt_GUT_RIL  = -KA_RIL * GUT_RIL;
dxdt_CENT_RIL = F_RIL * KA_RIL * GUT_RIL
                - (CL_RIL + Q_RIL) / V1_RIL * CENT_RIL
                + Q_RIL / V2_RIL * PERI_RIL;
dxdt_PERI_RIL = Q_RIL / V1_RIL * CENT_RIL
                - Q_RIL / V2_RIL * PERI_RIL;

// Edaravone — archetype 1 (single compartment, IV bolus/infusion via event)
dxdt_CENT_EDA = -(CL_EDA / V1_EDA) * CENT_EDA;

// Tofersen — archetype 3 (depot+central(plasma)+peripheral(CSF))
dxdt_GUT_TOF  = -KA_TOF * GUT_TOF;
dxdt_CENT_TOF = KA_TOF * GUT_TOF
                - (CL_TOF + Q_TOF) / V1_TOF * CENT_TOF
                + Q_TOF / V2_TOF * PERI_TOF;
dxdt_PERI_TOF = Q_TOF / V1_TOF * CENT_TOF
                - Q_TOF / V2_TOF * PERI_TOF;

// Phenylbutyrate (AMX0035) — archetype 3 minus peripheral (depot+central)
dxdt_GUT_PB  = -KA_PB * GUT_PB;
dxdt_CENT_PB = F_PB * KA_PB * GUT_PB - (CL_PB / V1_PB) * CENT_PB;

// ─── DISEASE BIOLOGY ODEs ────────────────────────────────────────
// Motor Neurons
dxdt_MN_upper  = -k_MN_death_U * MN_death_rate * MN_upper;
dxdt_MN_lower  = -k_MN_death   * MN_death_rate * MN_lower;

// SOD1 protein dynamics
dxdt_SOD1_wt   = SOD1_syn_eff
                 - k_SOD1_deg * SOD1_wt
                 - k_SOD1_mis * SOD1_wt * f_SOD1_mut;
dxdt_SOD1_mis  = k_SOD1_mis * SOD1_wt * f_SOD1_mut
                 - k_SOD1_clr * SOD1_mis;

// TDP-43 dynamics (nuclear ↔ cytoplasmic shuttle)
dxdt_TDP43_nuc  = -k_TDP_export * TDP43_nuc
                  + k_TDP_import * TDP43_cyto;
dxdt_TDP43_cyto =  k_TDP_export * TDP43_nuc
                  - k_TDP_import * TDP43_cyto
                  - k_TDP_agg   * TDP43_cyto;

// Synaptic glutamate
dxdt_Glut_syn  = k_Glu_rel * (1.0 - EFFECT_RIL) * MN_lower
                 - EAAT2_eff * Glut_syn;

// Intracellular calcium
dxdt_Ca_i      = k_Ca_entry * Glu_excess
                 - k_Ca_efflux * (Ca_i - Ca_i_REF);

// Reactive oxygen species
double ROS_prod = k_ROS_prod * (1.0 + SOD1_burden) * Ca_i / Mito;
double ROS_scav = k_ROS_scav * GSH * ROS + EFFECT_EDA * ROS;
dxdt_ROS       = ROS_prod - ROS_scav;

// Glutathione
dxdt_GSH       = k_GSH_syn - k_GSH_cons * GSH - k_ROS_scav * GSH * ROS;

// Mitochondrial function (0–1, 1=healthy)
dxdt_Mito      = k_Mito_rep * (1.0 - Mito) * (1.0 + EFFECT_PB_MITO)
                 - k_Mito_dam * ROS * Mito;

// Microglial activation
dxdt_Mic_act   = k_Mic_act * Mic_stim * (1.0 - Mic_act)
                 - k_Mic_res * Mic_act;

// TNF-α
dxdt_TNFa      = k_TNFa_prod * Mic_act - k_TNFa_deg * TNFa;

// BDNF trophic support
dxdt_BDNF      = k_BDNF_prod * (1.0 - 0.35 * TNFa_norm)
                 - k_BDNF_deg * BDNF;

// ─── BIOMARKER & CLINICAL ODEs ───────────────────────────────────
// CSF NfL (elevated with neurodegeneration)
dxdt_NfL_CSF   = NfL_release - k_NfL_clr * NfL_CSF;

// ALSFRS-R total score (declines proportionally to MN loss)
double MN_frac  = 0.5 * (MN_upper + MN_lower) / (0.5 * (MN_upper_REF + MN_lower_REF));
dxdt_ALSFRS    = -k_MN_death * ALSFRS * (1.8 - MN_frac) * (1.0 + 0.5 * Mic_act);

// FVC % predicted (respiratory drive declines with lower MN loss + inflammation)
dxdt_FVC       = -k_FVC_dec * FVC * (1.0 + TNFa_norm + Mic_act);

$TABLE
// C_RIL/C_EDA/C_TOF/Cp_TOF_plasma/C_PB/EFFECT_* are already declared as
// `double` locals in $MAIN above (matching where the original computed
// Cp_RIL/Cp_EDA/Ccsf_TOF/Cp_PB/E_RIL/E_EDA/E_TOF/E_PB/E_mito_PB, and
// needed there so the disease ODEs read them with the same timing the
// original used). They are refreshed here with a bare re-ASSIGNMENT
// (no `double`/`capture` keyword — re-declaring under `capture NAME =
// expr;` redefines the same symbol mrgsolve already hoists from $MAIN,
// and mrgsolve 2.0.1 rejects that outright: confirmed live,
// “redefinition of ‘capture {anonymous}::C_RIL’ ... previously declared
// here: ‘double {anonymous}::C_RIL’” — the same defect class as
// UPSTREAM_ISSUES.md’s “$CAPTURE duplicates compartment names”, just
// against a $MAIN double instead of a $CMT compartment).
//
// This reassignment is not merely a compile workaround — it fixes a
// real reporting-timing bug found live via the qspserver mrgsolve_api:
// $MAIN evaluates once per interval using state from the *start* of
// that interval (this guide’s own documented semantics), so a value
// read directly from the $MAIN double lags the true per-row state by
// one reporting step whenever a dose lands on a reporting time (proven
// with a minimal reproduction: a single-compartment test model showed
// C_X permanently offset by one row against its own CENT_X/V1_X, not
// just a self-healing dose-instant blip). $TABLE runs after $ODE has
// integrated up to the row actually being reported, so reassigning here
// (same formula, same variable) picks up the correct, contemporaneous
// value — this is exactly why the original recomputes Cp_RIL/Cp_EDA/
// Ccsf_TOF/Cp_PB_uM/E_RIL_cap/E_EDA_cap/E_TOF_cap in $TABLE instead of
// just exposing its own $MAIN locals directly. See refactor notes.
C_RIL          = CENT_RIL / V1_RIL;
C_EDA          = CENT_EDA / V1_EDA;
C_TOF          = PERI_TOF / V2_TOF;
Cp_TOF_plasma  = CENT_TOF / V1_TOF;
C_PB           = (CENT_PB / V1_PB) * 1000.0;
EFFECT_RIL     = EMAX_RIL * pow(C_RIL, GAMMA_RIL) / (pow(EC50_RIL, GAMMA_RIL) + pow(C_RIL, GAMMA_RIL) + 1e-9);
EFFECT_EDA     = EMAX_EDA * pow(C_EDA, GAMMA_EDA) / (pow(EC50_EDA, GAMMA_EDA) + pow(C_EDA, GAMMA_EDA) + 1e-9);
EFFECT_TOF     = EMAX_TOF * pow(C_TOF, GAMMA_TOF) / (pow(EC50_TOF, GAMMA_TOF) + pow(C_TOF, GAMMA_TOF) + 1e-9);
EFFECT_PB      = EMAX_PB * pow(C_PB, GAMMA_PB) / (pow(EC50_PB, GAMMA_PB) + pow(C_PB, GAMMA_PB) + 1e-9);
EFFECT_PB_MITO = EMAX_PB_MITO * pow(C_PB, GAMMA_PB) / (pow(EC50_PB, GAMMA_PB) + pow(C_PB, GAMMA_PB) + 1e-9);

capture MN_total    = MN_upper + MN_lower;
capture MN_pct      = 100.0 * (MN_upper + MN_lower) / 2.0;
capture SOD1_frac   = SOD1_mis / (SOD1_wt + SOD1_mis + 1e-9);
capture TDP_cyto_frac = TDP43_cyto / (TDP43_nuc + TDP43_cyto + 1e-9);
capture ROS_norm_out = ROS / ROS_REF;
capture MN_death_out = k_MN_death * (1.0 + 3.0 * f_SOD1_mut * SOD1_mis / (SOD1_wt + SOD1_mis + 1e-9));

$CAPTURE C_RIL C_EDA C_TOF Cp_TOF_plasma C_PB
         EFFECT_RIL EFFECT_EDA EFFECT_TOF EFFECT_PB EFFECT_PB_MITO