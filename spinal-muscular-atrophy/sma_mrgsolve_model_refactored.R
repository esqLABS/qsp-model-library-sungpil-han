## =============================================================================
## Spinal Muscular Atrophy (SMA) -- mrgsolve QSP Model
## PK/PD REFACTOR of sma_mrgsolve_model.R -- pluggable-PK naming convention
## (see ../FORK_WORKFLOW_GUIDE.md Part 2 and sma_refactor_notes.md)
## =============================================================================
## Disease: Spinal Muscular Atrophy (5q-SMA)
## Model covers:
##   1. SMN2 pre-mRNA alternative splicing (exon 7 inclusion dynamics)
##   2. SMN protein homeostasis (FL-SMN pool)
##   3. Motor neuron survival/degeneration (alpha-MN pool)
##   4. Neuromuscular junction maturation
##   5. Skeletal muscle mass dynamics
##   6. Clinical endpoint proxies (CMAP, HFMSE/CHOP-INTEND, FVC)
##   7. Nusinersen (NUS) intrathecal PK -- bespoke CSF-cascade + saturable
##      CNS-uptake transport (does not fit any of the guide's archetypes 1-4)
##   8. Risdiplam (RIS) oral PK -- archetype 3 minus peripheral (depot +
##      central, linear), plus a bespoke PD-silent CNS-partition compartment
##   9. Onasemnogene abeparvovec / Zolgensma (ZOL) IV AAV9 gene-therapy PK --
##      bespoke single concentration-valued compartment feeding an
##      irreversible transduction -> transgene-transcription cascade
##
## Calibration references:
##   - Darras et al. (2019) NEJM -- ENDEAR trial (nusinersen type I)
##   - Mercuri et al. (2018) NEJM -- CHERISH trial (nusinersen type II/III)
##   - Baranello et al. (2021) NEJM -- firefish/sunfish (risdiplam)
##   - Day et al. (2021) NEJM -- STR1VE (Zolgensma)
##   - Kletzl et al. (2019) J Pharmacokinet Pharmacodyn -- nusinersen PK
##   - Poirier et al. (2021) CPT:PSP -- risdiplam PK
##   - Al-Zaidy et al. (2019) Mol Ther -- Zolgensma biodistribution
##
## Usage:
##   library(mrgsolve); library(tidyverse)
##   source("sma_mrgsolve_model_refactored.R")
##   mod <- sma_model()
##   out <- mrgsim(mod, ev_nusinersen(), delta=1, end=730)
##
## Naming convention applied (see FORK_WORKFLOW_GUIDE.md Part 2):
##   Nusinersen  stem NUS: CENT_NUS/PERI_NUS/TISSUE_NUS, V1_NUS/V2_NUS/
##               V_CNS_NUS, Q_NUS, CL_NUS/CL_CNS_NUS, VMAX_NUS/KM_NUS,
##               EMAX_NUS/EC50_NUS/GAMMA_NUS, C_NUS, EFFECT_NUS
##   Risdiplam   stem RIS: GUT_RIS/CENT_RIS/CNSTISSUE_RIS, KA_RIS, F_RIS,
##               V1_RIS, CL_RIS, KP_RIS, EMAX_RIS/EC50_RIS/GAMMA_RIS,
##               C_RIS, EFFECT_RIS
##   Onasemnogene/Zolgensma stem ZOL: CENT_ZOL/MNLOAD_ZOL/MRNA_ZOL, CL_ZOL,
##               K_TRANS_ZOL/K_TXN_ZOL/K_MRNA_DEG_ZOL/K_PROT_ZOL/EFF_TG_ZOL,
##               C_ZOL, EFFECT_ZOL (bespoke -- not a Hill function, see notes)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ─────────────────────────────────────────────────────────────
## 1. Model Definition
## ─────────────────────────────────────────────────────────────
sma_model <- function() {
  mrgsolve::mcode("sma_qsp_refactored", '
$PROB
SMA QSP model: SMN biology, motor neuron degeneration,
nusinersen / risdiplam / onasemnogene-abeparvovec (Zolgensma) PK-PD
(PK/PD refactor: pluggable-PK naming convention, see FORK_WORKFLOW_GUIDE.md)

$PARAM @annotated
// --- SMN2 Splicing Parameters (disease-side, compound-agnostic) ---
E7I_base    : 0.10 : Baseline exon-7 inclusion from SMN2 (fraction, ~10%)
E7I_max     : 0.90 : Maximum achievable exon-7 inclusion fraction

// --- SMN mRNA & Protein Dynamics ---
k_SMN2_txn  : 1.0  : SMN2 transcription rate constant (relative units/day)
k_FL_deg    : 0.693: FL-SMN mRNA degradation rate (t1/2 ~1 day; /day)
k_d7_deg    : 2.079: SMN-Delta7 mRNA degradation rate (t1/2 ~8 h; /day)
k_prot_syn  : 2.0  : SMN protein synthesis rate per FL-mRNA unit
k_prot_deg  : 0.231: SMN protein degradation rate (t1/2 ~3 days; /day)
SMN_thresh  : 0.30 : SMN protein threshold (fraction of normal) below which MN death begins

// --- Motor Neuron Pool ---
MN0         : 1.0  : Initial motor neuron pool (normalized to 1.0 = 100%)
k_MN_death  : 0.002: Daily MN death rate in SMA (fraction/day)
k_MN_rescue : 0.04 : Rate at which SMN restoration rescues MNs (/day)
k_MN_spont  : 0.0001: Spontaneous MN death rate (aging; /day)
MN_min      : 0.05 : Minimum MN pool fraction (5% irreducible)
d_MN_hill   : 2.0  : Hill coefficient for SMN-dependent MN death rate

// --- Neuromuscular Junction ---
k_NMJ_mat   : 0.01 : Rate of NMJ maturation (/day)
k_NMJ_det   : 0.005: Rate of NMJ deterioration with MN loss (/day)
NMJ_max     : 1.0  : Maximum NMJ maturity score

// --- Skeletal Muscle ---
Muscle0     : 1.0  : Initial muscle mass (normalized)
k_muscle_at : 0.004: Daily muscle atrophy rate (denervation; /day)
k_muscle_gr : 0.002: Daily muscle recovery/growth rate (/day)
muscle_min  : 0.10 : Minimum muscle mass fraction

// --- Clinical Endpoint Scaling ---
CMAP_max    : 10.0 : Maximum CMAP amplitude (mV, normal adult)
CHOP_max    : 64.0 : Maximum CHOP-INTEND score
HFMSE_max   : 66.0 : Maximum HFMSE score
FVC_max     : 100.0: Maximum FVC % predicted
RULM_max    : 37.0 : Maximum RULM score

// --- Disease Subtype Parameters ---
SMN2_copies : 2.0  : Number of SMN2 gene copies (1-4)
disease_onset_age: 0.0: Age at symptom onset (months, for natural history)

// ============================================================
// Nusinersen (NUS) -- intrathecal ASO, bespoke CSF-cascade PK
// (bulk CSF flow lumbar->cervical + saturable Michaelis-Menten
//  uptake into CNS tissue; no archetype 1-4 fits this transport
//  shape -- see refactor notes)
// ============================================================
V1_NUS      : 30.0 : CSF lumbar ("central"/dosing-site) volume (mL)
V2_NUS      : 70.0 : CSF cervical+cranial ("peripheral") volume (mL)
V_CNS_NUS   : 1500.0: CNS tissue volume (g, for ng/g) [bespoke, no table slot]
Q_NUS       : 0.43 : CSF bulk flow rate, lumbar->cervical, one-way (mL/min = 619 mL/day)
CL_NUS      : 0.05 : ASO clearance from CSF pools (mL/min = 72 mL/day)
CL_CNS_NUS  : 0.001: ASO elimination from CNS tissue (fraction/day) [bespoke]
KA_CNS_NUS  : 0.05 : Uptake rate from CSF to CNS tissue (/day) [UNUSED -- see refactor notes]
VMAX_NUS    : 50.0 : Max uptake capacity from CSF to CNS (ng/day)
KM_NUS      : 5.0  : Michaelis constant for CNS uptake (ng/mL)
EMAX_NUS    : 0.60 : Nusinersen max fractional increase in exon-7 inclusion (0-1)
EC50_NUS    : 5.0  : Nusinersen CNS tissue EC50 for exon-7 inclusion (ng/g)
GAMMA_NUS   : 1.5  : Hill coefficient, nusinersen splicing effect

// ============================================================
// Risdiplam (RIS) -- oral small molecule
// Archetype 3 minus peripheral (depot + central, linear) for the
// plasma PK that PD actually reads, plus one bespoke, PD-silent
// "CNS tissue" compartment preserved from the original for fidelity
// (see refactor notes -- it is never read by any effect equation,
// in either the original or here)
// ============================================================
KA_RIS      : 1.5  : Oral absorption rate constant (/h, converted to /day in ODE)
F_RIS       : 0.99 : Oral bioavailability
V1_RIS      : 70.0 : Volume of distribution (L)
CL_RIS      : 1.0  : Clearance (L/h, converted to /day in ODE)
KP_RIS      : 0.5  : Brain:plasma partition ratio [bespoke, drives the PD-silent CNS compartment only]
FU_RIS      : 0.07 : Unbound fraction in plasma [UNUSED -- see refactor notes]
EMAX_RIS    : 0.50 : Risdiplam max fractional increase in exon-7 inclusion
EC50_RIS    : 80.0 : Risdiplam plasma EC50 for splicing effect (ng/mL)
GAMMA_RIS   : 1.2  : Hill coefficient, risdiplam splicing effect

// ============================================================
// Onasemnogene abeparvovec / Zolgensma (ZOL) -- IV AAV9 gene therapy
// Bespoke: single-compartment linear elimination (concentration-
// valued state, no separate volume term -- same "state IS the
// concentration" deviation as e.g. neonatal-hyperbilirubinemia\'s
// SNMP/PB/UDCA) feeding an irreversible transduction -> transgene
// transcription cascade. Not a Hill/occupancy relationship -- see
// refactor notes for why this is bespoke rather than Archetype 4.
// ============================================================
CL_ZOL      : 2.31 : Vector genome clearance rate from plasma (/day, t1/2 ~7h)
K_TRANS_ZOL : 0.001: Rate of motor neuron transduction (vg/nucleus/day) [bespoke]
K_TXN_ZOL   : 0.5  : Transgene transcription rate (relative units/day)
K_MRNA_DEG_ZOL: 0.693: Transgene mRNA degradation rate (t1/2 ~1 day; /day)
K_PROT_ZOL  : 2.0  : Transgene protein synthesis rate [feeds only a dead-code toggle -- see refactor notes]
EFF_TG_ZOL  : 1.0  : Transgene expression efficiency (0-1)
AAV9_AB_BLOCK_ZOL: 0.0 : Anti-AAV9 antibody neutralization (0=no block, 1=full block)

$CMT @annotated
// Nusinersen (NUS) PK -- bespoke CSF cascade
CENT_NUS    : CSF lumbar nusinersen, dosing/"central" site (ng)
PERI_NUS    : CSF cervical+cranial nusinersen, "peripheral" site (ng)
TISSUE_NUS  : CNS tissue nusinersen (ng/g equivalent) -- the PD-reading site

// Risdiplam (RIS) PK -- archetype 3 minus peripheral, plus bespoke tissue cmt
GUT_RIS     : GI tract risdiplam (mg)
CENT_RIS    : Plasma risdiplam (mg) -- the PD-reading site
CNSTISSUE_RIS: CNS-partitioned risdiplam (mg) -- PD-silent, preserved for fidelity

// Onasemnogene / Zolgensma (ZOL) PK + transduction cascade -- bespoke
CENT_ZOL    : Plasma vector genome, concentration-valued (vg/mL * Vdist)
MNLOAD_ZOL  : Transduced motor-neuron vg load (vg/nucleus)
MRNA_ZOL    : Transgene-derived FL-SMN mRNA (relative units) -- the PD-reading quantity

// Disease biology (unchanged)
FL_SMN_mRNA : Full-length SMN mRNA (relative units)
dSMN_mRNA   : SMN-delta7 mRNA (relative units)
SMN_prot    : Full-length SMN protein pool (normalized, 0-1)
MN_pool     : Alpha motor neuron pool (normalized, 0-1)
NMJ_score   : NMJ maturation score (0-1)
Muscle_mass : Skeletal muscle mass (normalized, 0-1)

// Cumulative endpoints
AUC_SMN     : AUC of SMN protein (for efficacy assessment)
MN_lost     : Cumulative MN loss (fraction)

$GLOBAL
// File-scope globals for the exposed concentrations and Hill/EFFECT
// interfaces -- assigned (no "double") in $ODE, reused as-is in
// $TABLE/$CAPTURE. Declaring them once here (rather than "double X = ..."
// inside both $ODE and $TABLE) avoids the mrgsolve pitfall where a
// block-local double of the same name in two DSL blocks is hoisted into
// one shared anonymous C++ namespace and collides -- the same mechanism
// behind this file\'s own pre-existing $TABLE/$CAPTURE build defect fixed
// below (see refactor notes and UPSTREAM_ISSUES.md).
double C_NUS, C_RIS, C_ZOL, EFFECT_NUS, EFFECT_RIS, EFFECT_ZOL;

$MAIN
// SMN2 copy-number scaling (more SMN2 copies = more baseline FL-SMN)
double SMN2_scale = SMN2_copies / 2.0;  // normalized to 2 copies

// Exon-7 inclusion at t=0 (drug-naive)
double E7I_ss = E7I_base * SMN2_scale;

// Initial SMN mRNA at steady state (no drug)
double FL_mRNA_ss = k_SMN2_txn * E7I_ss / k_FL_deg;
double dSMN_mRNA_ss = k_SMN2_txn * (1.0 - E7I_ss) / k_d7_deg;

// Initial SMN protein (normalized: normal = 1.0 for SMN2_copies=4)
double SMN_prot_ss = k_prot_syn * FL_mRNA_ss / k_prot_deg;

// Initial conditions
FL_SMN_mRNA_0  = FL_mRNA_ss;
dSMN_mRNA_0    = dSMN_mRNA_ss;
SMN_prot_0     = SMN_prot_ss;
MN_pool_0      = MN0;
NMJ_score_0    = 0.8 * MN0;   // NMJ roughly proportional to MN at baseline
Muscle_mass_0  = Muscle0;
AUC_SMN_0      = 0.0;
MN_lost_0      = 0.0;

$ODE
// -----------------------------------------------------------
// Nusinersen (NUS) PK -- intrathecal, bespoke CSF cascade
// CENT_NUS (dosing site, CSF lumbar) --Q_NUS(one-way)--> PERI_NUS (CSF cervical)
// CENT_NUS --saturable Vmax/Km uptake--> TISSUE_NUS (CNS tissue, the PD site)
// -----------------------------------------------------------
double conc_CENT_NUS = CENT_NUS / V1_NUS;
double conc_PERI_NUS = PERI_NUS / V2_NUS;

double flow_NUS   = (Q_NUS * 1440.0) * conc_CENT_NUS;      // ng/day (Q in mL/min -> mL/day)
double CLday_NUS  = CL_NUS * 1440.0;                       // mL/day
double uptake_NUS = (VMAX_NUS * conc_CENT_NUS) / (KM_NUS + conc_CENT_NUS);  // saturable CNS uptake

dxdt_CENT_NUS   = - flow_NUS - CLday_NUS * conc_CENT_NUS - uptake_NUS;
dxdt_PERI_NUS   =   flow_NUS - CLday_NUS * conc_PERI_NUS;
dxdt_TISSUE_NUS =   uptake_NUS - CL_CNS_NUS * TISSUE_NUS;

// C_NUS is the single concentration PD reads (CNS tissue, ng/g)
C_NUS = TISSUE_NUS / V_CNS_NUS;

// -----------------------------------------------------------
// Risdiplam (RIS) PK -- oral, depot + central (archetype 3 minus
// peripheral); CNSTISSUE_RIS is a second, PD-silent compartment
// preserved unchanged from the original (bespoke, non-standard
// Kp-scaled input/output -- see refactor notes)
// -----------------------------------------------------------
double KAday_RIS = KA_RIS * 24.0;   // convert /h to /day
double CLday_RIS = CL_RIS * 24.0;   // convert /h to /day

dxdt_GUT_RIS      = - KAday_RIS * GUT_RIS;
dxdt_CENT_RIS     =   KAday_RIS * F_RIS * GUT_RIS - CLday_RIS * (CENT_RIS / V1_RIS);
dxdt_CNSTISSUE_RIS =  KAday_RIS * KP_RIS * GUT_RIS - CLday_RIS * KP_RIS * (CNSTISSUE_RIS / V1_RIS);

// C_RIS is the single concentration PD reads (plasma, ng/mL)
C_RIS = (CENT_RIS / V1_RIS) * 1000.0;   // mg/L -> ng/mL

// -----------------------------------------------------------
// Onasemnogene / Zolgensma (ZOL) PK -- IV AAV9 vector, single
// concentration-valued compartment (bespoke, no explicit volume
// term in the original -- preserved, not flattened), feeding an
// irreversible transduction -> transgene transcription cascade
// -----------------------------------------------------------
dxdt_CENT_ZOL   = - CL_ZOL * CENT_ZOL;
dxdt_MNLOAD_ZOL =   K_TRANS_ZOL * CENT_ZOL * (1.0 - AAV9_AB_BLOCK_ZOL);
dxdt_MRNA_ZOL   =   K_TXN_ZOL * EFF_TG_ZOL * MNLOAD_ZOL - K_MRNA_DEG_ZOL * MRNA_ZOL;

// C_ZOL is the single concentration PD would read (vg/mL-equivalent);
// EFFECT_ZOL is the disease-facing quantity -- not a Hill/occupancy
// function (none exists mechanistically in the original), but the
// transgene mRNA pool itself, exactly as the original fed it straight
// into SMN synthesis (see refactor notes).
C_ZOL = CENT_ZOL;
EFFECT_ZOL = MRNA_ZOL;

// Pre-existing dead-code toggle in the original, preserved verbatim
// (SMN_from_ZOL is computed but never consumed by any dxdt/output --
// see refactor notes; not fixed, per "log what you find, don\'t fix").
double SMN_from_ZOL = K_PROT_ZOL > 0 ? MRNA_ZOL : 0.0;

// -----------------------------------------------------------
// Drug effects on exon-7 inclusion -- each compound\'s own named
// Hill interface, combined only at the point disease equations use them
// -----------------------------------------------------------
// Nusinersen effect (ISS-N1 blockade increases E7I) -- rename only,
// original was already this exact Hill shape
EFFECT_NUS = EMAX_NUS * pow(C_NUS, GAMMA_NUS) /
             (pow(EC50_NUS, GAMMA_NUS) + pow(C_NUS, GAMMA_NUS));

// Risdiplam effect (SRSF1/Tra2b enhancement) -- rename only
EFFECT_RIS = EMAX_RIS * pow(C_RIS, GAMMA_RIS) /
             (pow(EC50_RIS, GAMMA_RIS) + pow(C_RIS, GAMMA_RIS));

// Combined exon-7 inclusion rate (max at E7I_max) -- independent
// combination of NUS and RIS, unchanged from the original
double E7I_current = E7I_base + (E7I_max - E7I_base) * (EFFECT_NUS + EFFECT_RIS -
                     EFFECT_NUS * EFFECT_RIS);
E7I_current = fmin(E7I_current, E7I_max);
E7I_current = fmax(E7I_current, E7I_base);

// Scale by SMN2 copy number
E7I_current *= SMN2_scale;

// -----------------------------------------------------------
// SMN mRNA Dynamics
// -----------------------------------------------------------
dxdt_FL_SMN_mRNA = k_SMN2_txn * E7I_current - k_FL_deg * FL_SMN_mRNA;
dxdt_dSMN_mRNA   = k_SMN2_txn * (1.0 - fmin(E7I_current, 1.0)) - k_d7_deg * dSMN_mRNA;

// -----------------------------------------------------------
// SMN Protein Dynamics -- Zolgensma\'s contribution now reads the
// named EFFECT_ZOL interface (numerically identical to the
// original\'s direct A_tg_mRNA/MRNA_ZOL read, since EFFECT_ZOL is
// declared as its exact identity above)
// -----------------------------------------------------------
double SMN_synthesis = k_prot_syn * (FL_SMN_mRNA + EFFECT_ZOL);
dxdt_SMN_prot = SMN_synthesis - k_prot_deg * SMN_prot;
double SMN_norm = fmax(SMN_prot, 0.0);  // bounded at 0

// -----------------------------------------------------------
// Motor Neuron Pool Dynamics
// -----------------------------------------------------------
double SMN_ratio = SMN_norm / SMN_thresh;
double death_rate_factor = 1.0 / (1.0 + pow(SMN_ratio, d_MN_hill));
double k_death_effective = k_MN_death * death_rate_factor + k_MN_spont;

double rescue_factor = fmax(0.0, (SMN_norm - SMN_thresh) / (1.0 - SMN_thresh));
double k_rescue_effective = k_MN_rescue * rescue_factor;

dxdt_MN_pool = k_rescue_effective * (MN0 - MN_pool) - k_death_effective * MN_pool;
dxdt_MN_pool = fmax(dxdt_MN_pool, -(MN_pool - MN_min));  // don\'t go below minimum

// -----------------------------------------------------------
// NMJ Maturation
// -----------------------------------------------------------
dxdt_NMJ_score = k_NMJ_mat * MN_pool * (NMJ_max - NMJ_score)
                 - k_NMJ_det * (1.0 - MN_pool) * NMJ_score;

// -----------------------------------------------------------
// Skeletal Muscle Mass
// -----------------------------------------------------------
double innervation = NMJ_score * MN_pool;
dxdt_Muscle_mass = k_muscle_gr * innervation * (Muscle0 - Muscle_mass)
                   - k_muscle_at * (1.0 - innervation) * Muscle_mass;
dxdt_Muscle_mass = fmax(dxdt_Muscle_mass, -(Muscle_mass - muscle_min));

// -----------------------------------------------------------
// Cumulative Endpoints
// -----------------------------------------------------------
dxdt_AUC_SMN = SMN_prot;
dxdt_MN_lost = fmax(0.0, -dxdt_MN_pool);

$TABLE
// Derived clinical variables. NOTE: the original\'s own $TABLE re-declared
// local doubles CMAP/HFMSE/RULM with the SAME name it then passed to
// "capture NAME = NAME;" -- that self-referential form makes mrgsolve
// re-declare the symbol a second time as a `capture`, colliding with the
// `double` already declared for it and failing to compile (confirmed via
// qspserver /model_manifest on the untouched original; logged in
// UPSTREAM_ISSUES.md). FVC/CHOP_INTEND did NOT hit this because their own
// local variable names (FVC_pct, CHOP) already differed from their
// capture names. Fixed here, syntax-only, by giving CMAP/HFMSE/RULM their
// own differently-named locals too (CMAP_val/HFMSE_val/RULM_val) -- same
// pattern the original already used successfully for FVC/CHOP_INTEND.
double CMAP_val  = CMAP_max * MN_pool * NMJ_score;
double FVC_pct   = FVC_max * Muscle_mass * 0.85 + 15.0;  // respiratory muscle proxy
double HFMSE_val = HFMSE_max * pow(Muscle_mass, 1.5) * MN_pool;
double CHOP      = CHOP_max * Muscle_mass * NMJ_score;
double RULM_val  = RULM_max * Muscle_mass * MN_pool;

// Exon-7 inclusion (approximation from current drug levels). The
// original recomputed this independently from C_CNS_nusinersen and
// C_plasma_risdiplam using the SAME Hill formulas and the SAME
// underlying state values as EFFECT_NUS/EFFECT_RIS above -- so this is
// mathematically identical to reusing those two globals directly, which
// is what is done here (no numeric change, just removing a redundant
// second Hill evaluation of the same quantity). NOTE: this reporting
// formula pre-existing in the original is NOT the same as the actual
// disease-driving E7I_current in $ODE -- it omits both the
// EFFECT_NUS*EFFECT_RIS complement term and the SMN2_scale factor. This
// pre-existing inconsistency between the reported "E7_inclusion" output
// and the quantity that actually drives disease dynamics is preserved
// as-is (not fixed) and disclosed in the refactor notes and
// UPSTREAM_ISSUES.md.
// Refresh the exposed concentration/effect globals directly from
// state, exactly at THIS reported row -- found during this refactor\'s
// own verification (see refactor notes) that reusing an $ODE-block-
// local double (conc_CENT_NUS) for a directly-dosed compartment\'s
// concentration in $TABLE can read a value from the solver\'s last
// internal derivative evaluation rather than the state exactly at this
// reported time -- invisible for a continuously-varying compartment,
// but a full dose-sized artifact exactly at the row where a bolus
// lands on the compartment being read (e.g. C_CSF_lumbar right at a
// nusinersen dose). Reassigning (no "double", they are already
// $GLOBAL) guarantees freshness regardless of that internal solver
// detail, mirroring what the original always did by recomputing every
// $TABLE quantity directly from state rather than reusing an $ODE
// local.
double CENT_NUS_conc_out = CENT_NUS / V1_NUS;
C_NUS = TISSUE_NUS / V_CNS_NUS;
C_RIS = (CENT_RIS / V1_RIS) * 1000.0;
C_ZOL = CENT_ZOL;
EFFECT_NUS = EMAX_NUS * pow(C_NUS, GAMMA_NUS) /
             (pow(EC50_NUS, GAMMA_NUS) + pow(C_NUS, GAMMA_NUS));
EFFECT_RIS = EMAX_RIS * pow(C_RIS, GAMMA_RIS) /
             (pow(EC50_RIS, GAMMA_RIS) + pow(C_RIS, GAMMA_RIS));
EFFECT_ZOL = MRNA_ZOL;

double E7I_out = E7I_base + (E7I_max - E7I_base) * (EFFECT_NUS + EFFECT_RIS);

capture CMAP      = CMAP_val;
capture FVC       = FVC_pct;
capture HFMSE     = HFMSE_val;
capture CHOP_INTEND = CHOP;
capture RULM      = RULM_val;
capture E7_inclusion = E7I_out;
capture SMN_protein  = SMN_prot;
capture MN_fraction  = MN_pool;
capture NMJ_maturity = NMJ_score;
capture C_CSF_lumbar = CENT_NUS_conc_out;
capture C_CNS_nusinersen = C_NUS;
capture C_plasma_risdiplam = C_RIS;
capture Transgene_mRNA = MRNA_ZOL;

$CAPTURE
// C_NUS/C_RIS/C_ZOL/EFFECT_NUS/EFFECT_RIS/EFFECT_ZOL are already declared
// ($GLOBAL) and assigned ($ODE) above -- listed bare here (no "= expr")
// to avoid the exact re-declaration collision fixed above for
// CMAP/HFMSE/RULM, this time between $GLOBAL and $TABLE scope (same
// underlying mrgsolve mechanism, same fix: never write
// "capture NAME = NAME;" for a name already declared as a double
// elsewhere -- use a bare $CAPTURE entry instead).
C_NUS C_RIS C_ZOL EFFECT_NUS EFFECT_RIS EFFECT_ZOL
')
}

## ─────────────────────────────────────────────────────────────
## 2. Dosing Event Builders
##    (compartment order is unchanged from the original -- cmt indices
##    1/4/7 still address the nusinersen/risdiplam/Zolgensma dosing
##    sites, now named CENT_NUS/GUT_RIS/CENT_ZOL)
## ─────────────────────────────────────────────────────────────

# Nusinersen (intrathecal): ENDEAR/CHERISH loading + maintenance
# cmt=1 (CENT_NUS, was A_CSF_L), intrathecal 12 mg = 12000 ng injected
ev_nusinersen <- function(start_day = 0) {
  loading <- c(0, 14, 28, 63)
  maintenance <- seq(63 + 120, 730, by = 120)  # every 4 months
  days <- unique(c(loading + start_day, maintenance + start_day))
  days <- days[days <= 730]
  mrgsolve::ev(
    ID    = 1,
    time  = days,
    amt   = 12000,    # 12 mg = 12000 ng in CSF lumbar (cmt=1, CENT_NUS)
    cmt   = 1,
    evid  = 1,
    addl  = 0
  )
}

# Risdiplam oral 5 mg daily (fixed dose adult)
# amt in mg, cmt=4 (GUT_RIS, was A_gut_RIS)
ev_risdiplam <- function(start_day = 0, end_day = 730, dose_mg = 5) {
  mrgsolve::ev(
    ID    = 1,
    time  = seq(start_day, end_day, by = 1),
    amt   = dose_mg,
    cmt   = 4,
    evid  = 1
  )
}

# Risdiplam weight-based (pediatric) 0.2 mg/kg
ev_risdiplam_peds <- function(weight_kg = 15, start_day = 0, end_day = 730) {
  ev_risdiplam(start_day, end_day, dose_mg = 0.2 * weight_kg)
}

# Zolgensma single IV dose (1.1 x 10^14 vg/kg, 15 kg child = 1.65 x 10^15 vg)
# Expressed as normalized units; cmt=7 (CENT_ZOL, was A_plasma_ZOL)
ev_zolgensma <- function(dose_vg = 1.65e15) {
  mrgsolve::ev(
    ID   = 1,
    time = 0,
    amt  = dose_vg / 1e14,  # normalized to 10^14 units
    cmt  = 7,
    evid = 1
  )
}

## ─────────────────────────────────────────────────────────────
## 3. Simulation Scenarios
## ─────────────────────────────────────────────────────────────
run_scenarios <- function() {
  mod <- sma_model()

  # Parameter sets for disease subtypes
  params_type1 <- list(SMN2_copies = 2, MN0 = 1.0, k_MN_death = 0.004)
  params_type2 <- list(SMN2_copies = 3, MN0 = 1.0, k_MN_death = 0.002)
  params_type3 <- list(SMN2_copies = 4, MN0 = 1.0, k_MN_death = 0.001)

  scenarios <- list(
    # ── Scenario 1: Untreated SMA Type I natural history ──────
    list(
      name   = "SMA Type I — No Treatment",
      params = params_type1,
      events = mrgsolve::ev(time = 9999, amt = 0, cmt = 1)  # no drug
    ),

    # ── Scenario 2: Nusinersen in SMA Type I (ENDEAR) ─────────
    list(
      name   = "SMA Type I — Nusinersen",
      params = params_type1,
      events = ev_nusinersen(start_day = 0)
    ),

    # ── Scenario 3: Risdiplam in SMA Type II (SUNFISH) ────────
    list(
      name   = "SMA Type II — Risdiplam",
      params = params_type2,
      events = ev_risdiplam(dose_mg = 5)
    ),

    # ── Scenario 4: Zolgensma in presymptomatic SMA type I ───
    list(
      name   = "Presymptomatic SMA — Zolgensma",
      params = modifyList(params_type1, list(MN0 = 0.95)),
      events = ev_zolgensma(dose_vg = 1.65e15)
    ),

    # ── Scenario 5: Nusinersen late start (type II, 2 yrs) ───
    list(
      name   = "SMA Type II — Nusinersen Late Start",
      params = params_type2,
      events = ev_nusinersen(start_day = 365)  # start after 1 year
    ),

    # ── Scenario 6: Risdiplam pediatric weight-based ─────────
    list(
      name   = "SMA Type II — Risdiplam Pediatric (15 kg)",
      params = params_type2,
      events = ev_risdiplam_peds(weight_kg = 15)
    )
  )

  results <- lapply(seq_along(scenarios), function(i) {
    sc  <- scenarios[[i]]
    mod2 <- mrgsolve::param(mod, sc$params)
    out  <- mrgsolve::mrgsim(mod2, sc$events, delta = 1, end = 730, obsonly = TRUE)
    df   <- as.data.frame(out)
    df$Scenario <- sc$name
    df
  })

  bind_rows(results)
}

## ─────────────────────────────────────────────────────────────
## 4. Visualization
## ─────────────────────────────────────────────────────────────
plot_results <- function(df) {
  pal <- c(
    "SMA Type I — No Treatment"                     = "#B71C1C",
    "SMA Type I — Nusinersen"                       = "#4CAF50",
    "SMA Type II — Risdiplam"                       = "#2196F3",
    "Presymptomatic SMA — Zolgensma"                = "#9C27B0",
    "SMA Type II — Nusinersen Late Start"           = "#FF9800",
    "SMA Type II — Risdiplam Pediatric (15 kg)"    = "#00BCD4"
  )
  lt <- c(
    "SMA Type I — No Treatment"                     = "dashed",
    "SMA Type I — Nusinersen"                       = "solid",
    "SMA Type II — Risdiplam"                       = "solid",
    "Presymptomatic SMA — Zolgensma"                = "solid",
    "SMA Type II — Nusinersen Late Start"           = "dotdash",
    "SMA Type II — Risdiplam Pediatric (15 kg)"    = "dotted"
  )

  p1 <- ggplot(df, aes(time, SMN_protein, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "SMN Protein Pool Over Time",
         x = "Days", y = "SMN Protein (normalized)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p2 <- ggplot(df, aes(time, MN_fraction, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "Motor Neuron Pool Over Time",
         x = "Days", y = "MN Pool (fraction of baseline)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p3 <- ggplot(df, aes(time, CMAP, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "CMAP Amplitude Over Time",
         x = "Days", y = "CMAP (mV)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p4 <- ggplot(df, aes(time, CHOP_INTEND, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "CHOP-INTEND Score Over Time",
         x = "Days", y = "CHOP-INTEND (0–64)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p5 <- ggplot(df, aes(time, FVC, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "FVC % Predicted Over Time",
         x = "Days", y = "FVC (% predicted)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p6 <- ggplot(df, aes(time, E7_inclusion, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "Exon 7 Inclusion Rate (SMN2)",
         x = "Days", y = "Exon 7 Inclusion (fraction)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
}

## ─────────────────────────────────────────────────────────────
## 5. Sensitivity Analysis
## ─────────────────────────────────────────────────────────────
sensitivity_analysis <- function() {
  mod <- sma_model()

  # Vary SMN2 copy number (1–4)
  sa_smn2 <- lapply(1:4, function(copies) {
    m <- mrgsolve::param(mod, list(SMN2_copies = copies, k_MN_death = 0.003))
    ev_null <- mrgsolve::ev(time = 9999, amt = 0, cmt = 1)
    df <- as.data.frame(mrgsolve::mrgsim(m, ev_null, delta = 1, end = 365))
    df$SMN2_copies <- copies
    df
  }) |> bind_rows()

  # Vary nusinersen EC50 (EC50_NUS -- unchanged name, was already this)
  ec50_vals <- c(2.5, 5.0, 10.0, 20.0)
  sa_ec50 <- lapply(ec50_vals, function(ec50) {
    m <- mrgsolve::param(mod, list(EC50_NUS = ec50, SMN2_copies = 2))
    df <- as.data.frame(mrgsolve::mrgsim(m, ev_nusinersen(), delta = 1, end = 730))
    df$EC50_NUS <- ec50
    df
  }) |> bind_rows()

  list(smn2_copies = sa_smn2, ec50_nusinersen = sa_ec50)
}

## ─────────────────────────────────────────────────────────────
## 6. Virtual Population (PopPK)
## ─────────────────────────────────────────────────────────────
virtual_population <- function(n = 100, seed = 42) {
  set.seed(seed)
  mod <- sma_model()

  # Sample interindividual variability
  # NOTE: `Emax_NUS` was renamed `EMAX_NUS` by this refactor (naming
  # convention EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>) -- updated here to
  # match; `EC50_NUS`/`k_prot_deg`/`k_MN_death`/`SMN2_copies` were already
  # spelled this way in the original and are unchanged.
  idata <- data.frame(
    ID           = 1:n,
    SMN2_copies  = sample(2:3, n, replace = TRUE, prob = c(0.6, 0.4)),
    k_MN_death   = rlnorm(n, log(0.002), 0.3),
    EC50_NUS     = rlnorm(n, log(5.0), 0.4),
    EMAX_NUS     = rnorm(n, 0.60, 0.05) |> pmin(0.85) |> pmax(0.35),
    k_prot_deg   = rlnorm(n, log(0.231), 0.2)
  )

  ev_nus <- ev_nusinersen()
  out <- mrgsolve::mrgsim(mod, idata = idata, events = ev_nus,
                          delta = 7, end = 730, obsonly = TRUE)
  as.data.frame(out)
}

## ─────────────────────────────────────────────────────────────
## 7. Main: Run Everything
## ─────────────────────────────────────────────────────────────
if (FALSE) {  # set TRUE to run interactively
  library(mrgsolve); library(dplyr); library(ggplot2); library(cowplot)

  cat("Building SMA QSP model (refactored)...\n")
  mod <- sma_model()
  print(mod)

  cat("Running 6 treatment scenarios...\n")
  df_all <- run_scenarios()

  cat("Plotting results...\n")
  p <- plot_results(df_all)
  print(p)

  cat("Running sensitivity analysis...\n")
  sa <- sensitivity_analysis()
  print(sa$smn2_copies |> filter(time == 365) |>
    group_by(SMN2_copies) |> summarise(SMN_prot = last(SMN_prot), MN = last(MN_pool)))

  cat("Generating virtual population (n=100)...\n")
  vpc <- virtual_population(100)
  cat("VPC summary:\n")
  print(vpc |> filter(time == 365) |>
    group_by(ID) |> slice(n()) |>
    summarise(across(c(MN_fraction, CMAP, CHOP_INTEND, FVC),
                     list(median = median, q5 = ~quantile(., 0.05), q95 = ~quantile(., 0.95)))))
}
