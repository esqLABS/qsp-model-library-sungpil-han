## ============================================================
##  Fibromyalgia (FM) — mrgsolve QSP Model — PK/PD-refactored sibling
##  (see fm_refactor_notes.md; original untouched at fm_mrgsolve_model.R)
##
##  Compartments:
##    Drug PK  : Duloxetine (DUL), Pregabalin (PRE),
##               Milnacipran (MIL), Amitriptyline (TCA)
##    Peripheral: DRG sensitization, NGF, PGE2
##    Spinal   : Dorsal-horn WDR, Substance P (CSF),
##               NMDA-receptor state, Wind-up, LTP (central sensitization)
##    Brain    : Synaptic NE, 5-HT, descending inhibition (DPMS)
##    HPA Axis : CRH, ACTH, Cortisol, negative feedback
##    ANS      : SNS tone, HRV
##    Sleep    : Sleep pressure (adenosine), SWS depth
##    Immune   : Microglia activation, IL-1beta (spinal)
##    Outcomes : Pain score, FIQ, fatigue, depression
##  Total: 31 ODE compartments (9 drug PK + 22 disease; the original file's own
##  header comment claimed "10" PK / "30" total -- both off by one against its
##  own $CMT list, a pre-existing doc-only inconsistency, not fixed here since
##  it is cosmetic and outside this refactor's scope).
##
##  PK/PD refactor summary (see fm_refactor_notes.md for full detail):
##    - Every compound's PK isolated into its own GUT_<STEM>/CENT_<STEM>/
##      PERI_<STEM> compartments (naming convention from FORK_WORKFLOW_GUIDE.md
##      Part 2). DUL keeps its peripheral compartment (Archetype 3, full);
##      PRE/MIL/TCA have no peripheral compartment in the original and keep
##      none here (Archetype 3 without PERI_<STEM>).
##    - Each compound exposes exactly one concentration site, C_<STEM>,
##      computed once in $ODE (the block the original itself used); the
##      original's duplicate "Cp_DUL / Cp_DUL_out"-style pair (one name in
##      $ODE, a second name recomputing the identical ratio in $TABLE) is
##      normalized to the single canonical name, referenced directly in
##      $TABLE/$CAPTURE rather than recomputed under a second name.
##    - Every compound's effect on the disease system is an explicit named
##      Hill term, EFFECT_<STEM_PATHWAY>, built from EMAX_/EC50_/GAMMA_
##      parameters carrying the original's values unchanged (GAMMA=1 where
##      the original had no explicit Hill coefficient) -- a rename of an
##      already-plain-ratio effect term, not a refit. DUL/MIL/TCA each act
##      on two distinct disease-facing pathways (SERT and NET reuptake
##      inhibition) plus, for TCA alone, a third (antihistaminergic/
##      anticholinergic sedation feeding sleep), so each gets its own
##      pathway-qualified EFFECT_<STEM>_<PATHWAY> rather than one collapsed
##      EFFECT_<STEM> that would hide which pathway drives which disease
##      equation. PRE has one pathway (alpha2-delta channel block) and one
##      EFFECT_PRE term.
##    - Pre-existing upstream build defect found and fixed *only* in this
##      sibling (never in the original): $CAPTURE re-listed 18 names that are
##      already $CMT compartments, which this mrgsolve build's validObject()
##      rejects outright. Logged as UPSTREAM_ISSUES.md #144; fixed here by
##      simply not re-listing compartment names in $CAPTURE (mrgsolve reports
##      compartment state automatically) -- a syntax-only change, confirmed
##      numerically inert by the verification run in fm_refactor_notes.md.
## ============================================================

library(mrgsolve)

fm_code <- '
$PROB
Fibromyalgia QSP — Duloxetine / Pregabalin / Milnacipran / Amitriptyline PK-PD
31-compartment ODE model including:
  Drug PK (4 drugs, pluggable per-compound PK blocks), Peripheral sensitization,
  Central sensitization (spinal/brain), HPA axis, ANS, Sleep dynamics,
  Neuroinflammation, Clinical outcomes

$PARAM @annotated
// ---- Duloxetine (DUL) PK -- Archetype 3, depot + central + peripheral ----
KA_DUL  : 0.80   : Duloxetine oral absorption rate (1/h)
CL_DUL  : 54.0   : Duloxetine total clearance (L/h)
V1_DUL  : 1640   : Duloxetine central volume (L)
Q_DUL   : 18.0   : Duloxetine inter-compartmental clearance (L/h)
V2_DUL  : 820    : Duloxetine peripheral volume (L)
F_DUL   : 0.50   : Duloxetine bioavailability

// ---- Pregabalin (PRE) PK -- Archetype 3 without peripheral ----
KA_PRE  : 1.30   : Pregabalin absorption rate (1/h, rapid absorption)
CL_PRE  : 6.8    : Pregabalin renal clearance (L/h, CLcr-based)
V1_PRE  : 42.0   : Pregabalin central volume (L)
F_PRE   : 0.90   : Pregabalin bioavailability

// ---- Milnacipran (MIL) PK -- Archetype 3 without peripheral ----
KA_MIL  : 0.90   : Milnacipran absorption rate (1/h)
CL_MIL  : 50.0   : Milnacipran clearance (L/h)
V1_MIL  : 300    : Milnacipran central volume (L)
F_MIL   : 0.85   : Milnacipran bioavailability

// ---- Amitriptyline (TCA) PK -- Archetype 3 without peripheral ----
KA_TCA  : 0.70   : Amitriptyline absorption rate (1/h)
CL_TCA  : 40.0   : Amitriptyline clearance (L/h)
V1_TCA  : 1500   : Amitriptyline central volume (L)
F_TCA   : 0.48   : Amitriptyline bioavailability

// ---- Hill interface: SERT-reuptake-inhibition pathway (DUL, MIL, TCA) ----
EC50_DUL_SERT   : 0.003 : Duloxetine EC50, SERT-inhibition pathway (mg/L)
EMAX_DUL_SERT   : 1.0   : Duloxetine Emax, SERT-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_DUL_SERT  : 1.0   : Duloxetine Hill coefficient, SERT-inhibition pathway (none in original; =1, rename not a fit)
EC50_MIL_SERT   : 0.015 : Milnacipran EC50, SERT-inhibition pathway (mg/L)
EMAX_MIL_SERT   : 1.0   : Milnacipran Emax, SERT-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_MIL_SERT  : 1.0   : Milnacipran Hill coefficient, SERT-inhibition pathway (none in original; =1, rename not a fit)
EC50_TCA_SERT   : 0.010 : Amitriptyline EC50, SERT-inhibition pathway (mg/L)
EMAX_TCA_SERT   : 1.0   : Amitriptyline Emax, SERT-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_TCA_SERT  : 1.0   : Amitriptyline Hill coefficient, SERT-inhibition pathway (none in original; =1, rename not a fit)

// ---- Hill interface: NET-reuptake-inhibition pathway (DUL, MIL, TCA) ----
EC50_DUL_NET    : 0.012 : Duloxetine EC50, NET-inhibition pathway (mg/L)
EMAX_DUL_NET    : 1.0   : Duloxetine Emax, NET-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_DUL_NET   : 1.0   : Duloxetine Hill coefficient, NET-inhibition pathway (none in original; =1, rename not a fit)
EC50_MIL_NET    : 0.018 : Milnacipran EC50, NET-inhibition pathway (mg/L)
EMAX_MIL_NET    : 1.0   : Milnacipran Emax, NET-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_MIL_NET   : 1.0   : Milnacipran Hill coefficient, NET-inhibition pathway (none in original; =1, rename not a fit)
EC50_TCA_NET    : 0.025 : Amitriptyline EC50, NET-inhibition pathway (mg/L)
EMAX_TCA_NET    : 1.0   : Amitriptyline Emax, NET-inhibition pathway (dimensionless; original Emax_SNRI, shared, split per stem/pathway -- rename not a refit)
GAMMA_TCA_NET   : 1.0   : Amitriptyline Hill coefficient, NET-inhibition pathway (none in original; =1, rename not a fit)

// ---- Hill interface: Amitriptyline sedation/SWS-support pathway ----
EC50_TCA_SEDATION  : 0.05 : Amitriptyline EC50, sedative SWS-support pathway (mg/L)
EMAX_TCA_SEDATION  : 1.0  : Amitriptyline Emax, sedative SWS-support pathway (dimensionless; implicit in original plain ratio -- rename not a refit)
GAMMA_TCA_SEDATION : 1.0  : Amitriptyline Hill coefficient, sedative SWS-support pathway (none in original; =1, rename not a fit)

// ---- Hill interface: Pregabalin alpha2-delta calcium-channel pathway ----
EC50_PRE  : 0.05  : Pregabalin EC50, alpha2-delta channel-block pathway (mg/L)
EMAX_PRE  : 0.70  : Pregabalin Emax, alpha2-delta channel-block pathway (max fractional block)
GAMMA_PRE : 1.0   : Pregabalin Hill coefficient, alpha2-delta channel-block pathway (none in original; =1, rename not a fit)

$PARAM
// ---- Peripheral sensitization ----
kprod_NGF   = 0.05   // baseline NGF synthesis h^-1
kdeg_NGF    = 0.10   // NGF degradation h^-1
kact_TRPV1  = 0.20   // PGE2/NGF -> TRPV1 sensitization
kdeg_PGE2   = 0.30   // h^-1
kprod_PGE2  = 0.08   // h^-1 baseline
kstim_DRG   = 0.40   // DRG firing rate constant
kdeg_DRG    = 0.50   // DRG decay h^-1

// ---- Spinal / Central sensitization ----
kprod_SP    = 0.15   // Substance P production h^-1
kdeg_SP     = 0.20   // h^-1
kWU         = 0.05   // wind-up accumulation h^-1
kWU_decay   = 0.02   // wind-up decay h^-1
kLTP        = 0.08   // LTP induction rate h^-1
kLTP_decay  = 0.005  // LTP spontaneous decay h^-1
kNMDA_act   = 0.12   // NMDA receptor activation by Glu/SP
kNMDA_decay = 0.08   // NMDA inactivation h^-1
Emax_inhib  = 0.70   // max inhibitory effect of descending DPMS

// ---- Synaptic transmitters (brain) ----
ksyn_NE     = 0.30   // NE synthesis/pool h^-1
ksyn_5HT    = 0.25   // 5-HT synthesis h^-1
kdeg_NE     = 0.40   // NE reuptake/degradation h^-1
kdeg_5HT    = 0.35   // 5-HT reuptake/degradation h^-1
kdesc_NE    = 0.20   // descending NE -> spinal inhibition
kdesc_5HT   = 0.15   // descending 5-HT -> spinal inhibition

// ---- HPA axis ----
kprod_CRH   = 0.50   // h^-1
kdeg_CRH    = 0.80   // h^-1
kprod_ACTH  = 0.40   // h^-1
kdeg_ACTH   = 0.60   // h^-1
kprod_CORT  = 0.30   // h^-1
kdeg_CORT   = 0.20   // h^-1
kfb_CORT    = 0.60   // cortisol negative feedback strength

// ---- ANS ----
kSNS_base   = 0.60   // baseline SNS tone
kSNS_stress = 0.20   // stress-to-SNS coupling
kSNS_decay  = 0.40   // h^-1
kHRV_base   = 0.80   // baseline HRV (normalized)

// ---- Sleep ----
kaden_prod  = 0.15   // adenosine production h^-1
kaden_clear = 0.12   // adenosine clearance h^-1
kSWS_drive  = 0.20   // adenosine -> SWS drive
kSWS_decay  = 0.25   // SWS decay h^-1
kSWS_pain_inh = 0.30 // pain -> SWS inhibition coefficient
kSWS_TCA    = 0.15   // TCA sedation -> SWS support

// ---- Neuroinflammation ----
kprod_MG    = 0.10   // microglia activation h^-1
kdeg_MG     = 0.12   // h^-1
kprod_IL1b  = 0.18   // IL-1beta production h^-1
kdeg_IL1b   = 0.25   // h^-1
kMG_cortisol = 0.15  // cortisol -> microglia suppression

// ---- Clinical outcomes ----
k_pain_LTP  = 0.40   // LTP -> pain score scaling
k_pain_SP   = 0.20   // CSF SP -> pain score
k_FIQ_pain  = 0.35   // pain -> FIQ
k_FIQ_sleep = 0.25   // sleep deprivation -> FIQ
k_FIQ_dep   = 0.20   // depression -> FIQ
k_fatigue   = 0.30   // sleep loss -> fatigue
k_dep_LTP   = 0.15   // central sensitization -> depression
pain_base   = 5.0    // baseline NRS pain score
FIQ_base    = 55.0   // baseline FIQR score
fatigue_base = 60.0  // baseline fatigue VAS

// ---- Dosing flags ----
use_DUL = 0     // 0=off, 1=on
use_PRE = 0
use_MIL = 0
use_TCA = 0

$CMT
// Drug PK compartments (9)
GUT_DUL CENT_DUL PERI_DUL     // duloxetine: gut, central, peripheral
GUT_PRE CENT_PRE              // pregabalin: gut, central
GUT_MIL CENT_MIL               // milnacipran: gut, central
GUT_TCA CENT_TCA                // amitriptyline: gut, central

// Peripheral sensitization (3)
NGF PGE2 DRG_act

// Spinal (6)
SP_csf NMDA_state WindUp LTP_cs NE_syn SHT_syn

// HPA axis (3)
CRH ACTH CORT

// ANS + Sleep (4)
SNS_tone SWS_depth Adenosine DPMS

// Neuroinflammation (2)
MG_act IL1b_sp

// Clinical outcomes (4) -- tracked as ODEs for smoothing
Pain_score FIQ_score Fatigue_VAS Depression_score

$MAIN
// ----- PK initial steady-state approximations -----
// Pain begins at FM steady state (untreated)
if(NEWIND <= 1) {
  // start at baseline disease state
}

$ODE
// ============================================================
// Drug PK (each compound isolated in its own block; single exposed
// concentration C_<STEM> is the only site PD equations below read)
// ============================================================
// -- Duloxetine (Archetype 3: depot + central + peripheral) --
double ka_DUL_eff = use_DUL * KA_DUL;
dxdt_GUT_DUL  = -ka_DUL_eff * GUT_DUL;
dxdt_CENT_DUL =  ka_DUL_eff * GUT_DUL * F_DUL
                 - (CL_DUL/V1_DUL) * CENT_DUL
                 - (Q_DUL/V1_DUL)  * CENT_DUL
                 + (Q_DUL/V2_DUL)  * PERI_DUL;
dxdt_PERI_DUL =  (Q_DUL/V1_DUL)   * CENT_DUL
                 - (Q_DUL/V2_DUL)  * PERI_DUL;
double C_DUL  = CENT_DUL / V1_DUL;

// -- Pregabalin (Archetype 3 without peripheral) --
double ka_PRE_eff = use_PRE * KA_PRE;
dxdt_GUT_PRE  = -ka_PRE_eff * GUT_PRE;
dxdt_CENT_PRE =  ka_PRE_eff * GUT_PRE * F_PRE - (CL_PRE/V1_PRE) * CENT_PRE;
double C_PRE  = CENT_PRE / V1_PRE;

// -- Milnacipran (Archetype 3 without peripheral) --
double ka_MIL_eff = use_MIL * KA_MIL;
dxdt_GUT_MIL  = -ka_MIL_eff * GUT_MIL;
dxdt_CENT_MIL =  ka_MIL_eff * GUT_MIL * F_MIL - (CL_MIL/V1_MIL) * CENT_MIL;
double C_MIL  = CENT_MIL / V1_MIL;

// -- Amitriptyline (Archetype 3 without peripheral) --
double ka_TCA_eff = use_TCA * KA_TCA;
dxdt_GUT_TCA  = -ka_TCA_eff * GUT_TCA;
dxdt_CENT_TCA =  ka_TCA_eff * GUT_TCA * F_TCA - (CL_TCA/V1_TCA) * CENT_TCA;
double C_TCA  = CENT_TCA / V1_TCA;

// ============================================================
// Drug PD — named Hill terms per compound per pathway, combined only
// at the point the disease equations actually use them
// ============================================================
// SERT-inhibition pathway (DUL, MIL, TCA each contribute independently)
double EFFECT_DUL_SERT = EMAX_DUL_SERT * pow(C_DUL, GAMMA_DUL_SERT)
                          / (pow(EC50_DUL_SERT, GAMMA_DUL_SERT) + pow(C_DUL, GAMMA_DUL_SERT));
double EFFECT_MIL_SERT = EMAX_MIL_SERT * pow(C_MIL, GAMMA_MIL_SERT)
                          / (pow(EC50_MIL_SERT, GAMMA_MIL_SERT) + pow(C_MIL, GAMMA_MIL_SERT));
double EFFECT_TCA_SERT = EMAX_TCA_SERT * pow(C_TCA, GAMMA_TCA_SERT)
                          / (pow(EC50_TCA_SERT, GAMMA_TCA_SERT) + pow(C_TCA, GAMMA_TCA_SERT));
double inh_SERT = EFFECT_DUL_SERT + EFFECT_MIL_SERT + EFFECT_TCA_SERT;
inh_SERT = (inh_SERT > 1.0) ? 1.0 : inh_SERT;

// NET-inhibition pathway (DUL, MIL, TCA each contribute independently)
double EFFECT_DUL_NET = EMAX_DUL_NET * pow(C_DUL, GAMMA_DUL_NET)
                         / (pow(EC50_DUL_NET, GAMMA_DUL_NET) + pow(C_DUL, GAMMA_DUL_NET));
double EFFECT_MIL_NET = EMAX_MIL_NET * pow(C_MIL, GAMMA_MIL_NET)
                         / (pow(EC50_MIL_NET, GAMMA_MIL_NET) + pow(C_MIL, GAMMA_MIL_NET));
double EFFECT_TCA_NET = EMAX_TCA_NET * pow(C_TCA, GAMMA_TCA_NET)
                         / (pow(EC50_TCA_NET, GAMMA_TCA_NET) + pow(C_TCA, GAMMA_TCA_NET));
double inh_NET  = EFFECT_DUL_NET + EFFECT_MIL_NET + EFFECT_TCA_NET;
inh_NET = (inh_NET > 1.0) ? 1.0 : inh_NET;

// Alpha2-delta channel block (pregabalin, single pathway)
double EFFECT_PRE = EMAX_PRE * pow(C_PRE, GAMMA_PRE)
                     / (pow(EC50_PRE, GAMMA_PRE) + pow(C_PRE, GAMMA_PRE));

// TCA sedation -> SWS support (third, independent pathway for TCA)
double EFFECT_TCA_SEDATION = EMAX_TCA_SEDATION * pow(C_TCA, GAMMA_TCA_SEDATION)
                              / (pow(EC50_TCA_SEDATION, GAMMA_TCA_SEDATION) + pow(C_TCA, GAMMA_TCA_SEDATION));

// ============================================================
// Peripheral Sensitization
// ============================================================
// NGF dynamics (NGF elevated in FM muscle/skin)
dxdt_NGF  = kprod_NGF * 1.5          // elevated production in FM
             - kdeg_NGF * NGF;

// PGE2 dynamics (mast cell / COX pathway)
dxdt_PGE2 = kprod_PGE2 * (1.0 + IL1b_sp * 0.5)  // IL-1beta amplifies PGE2
             - kdeg_PGE2 * PGE2;

// DRG afferent activity (A-delta/C fibers)
double TRPV1_sens = kact_TRPV1 * (PGE2 + NGF * 0.5);
dxdt_DRG_act = kstim_DRG * TRPV1_sens
               - kdeg_DRG * DRG_act;

// ============================================================
// Spinal Dorsal Horn — Central Sensitization
// ============================================================
// Substance P (CSF proxy)
dxdt_SP_csf = kprod_SP * DRG_act * (1.0 - EFFECT_PRE * 0.6)
              - kdeg_SP * SP_csf;

// NMDA receptor state (0–1 scale, 0=resting)
dxdt_NMDA_state = kNMDA_act * (SP_csf + DRG_act * 0.5) * (1.0 - NMDA_state)
                  - kNMDA_decay * NMDA_state;

// Wind-up (WU): repetitive C-fiber activity -> temporal summation
dxdt_WindUp = kWU * DRG_act * NMDA_state
              - kWU_decay * WindUp;

// Long-term potentiation / central sensitization index (0–1)
dxdt_LTP_cs = kLTP * WindUp * (1.0 + IL1b_sp * 0.4)
              - kLTP_decay * LTP_cs
              - LTP_cs * DPMS * Emax_inhib; // descending inhibition

// ============================================================
// Supraspinal — Synaptic Monoamines & Descending Inhibition
// ============================================================
// Synaptic NE (descending)
dxdt_NE_syn  = ksyn_NE  * (1.0 + inh_NET  * 2.0)  // NET block -> up NE
               - kdeg_NE * (1.0 - inh_NET) * NE_syn;

// Synaptic 5-HT (descending)
dxdt_SHT_syn = ksyn_5HT * (1.0 + inh_SERT * 2.0)  // SERT block -> up 5-HT
               - kdeg_5HT * (1.0 - inh_SERT) * SHT_syn;

// Descending Pain Modulating System (DPMS, 0–1 scale)
dxdt_DPMS    = kdesc_NE * NE_syn + kdesc_5HT * SHT_syn
               - 0.30 * DPMS      // intrinsic turnover
               - DPMS * SNS_tone * 0.10;  // SNS stress partially offsets

// ============================================================
// HPA Axis
// ============================================================
// CRH (hypothalamus) — stress + pain driven
dxdt_CRH  = kprod_CRH * (1.0 + LTP_cs * 0.5 + SNS_tone * 0.3)
             - kdeg_CRH  * CRH
             - kfb_CORT  * CORT * CRH;  // negative feedback

// ACTH (pituitary)
dxdt_ACTH = kprod_ACTH * CRH - kdeg_ACTH * ACTH;

// Cortisol (adrenal)
dxdt_CORT  = kprod_CORT * ACTH - kdeg_CORT * CORT;

// ============================================================
// ANS — Sympathetic tone
// ============================================================
dxdt_SNS_tone = kSNS_base + kSNS_stress * (LTP_cs + 0.5 * (1.0 - CORT * 0.5))
                - kSNS_decay * SNS_tone;

// ============================================================
// Sleep — SWS depth and adenosine
// ============================================================
// Adenosine (sleep pressure)
dxdt_Adenosine = kaden_prod - kaden_clear * Adenosine
                 + SNS_tone * 0.05;  // stress slows clearance

// Slow-wave sleep depth (0–1)
dxdt_SWS_depth = kSWS_drive * Adenosine
                 - kSWS_decay * SWS_depth
                 - kSWS_pain_inh * LTP_cs * SWS_depth   // pain disrupts SWS
                 + kSWS_TCA * EFFECT_TCA_SEDATION;        // TCA sedation

// ============================================================
// Neuroinflammation
// ============================================================
// Microglia activation (spinal)
dxdt_MG_act  = kprod_MG * DRG_act * (1.0 + WindUp * 0.5)
               - kdeg_MG * MG_act
               - kMG_cortisol * CORT * MG_act;

// IL-1beta (spinal)
dxdt_IL1b_sp = kprod_IL1b * MG_act - kdeg_IL1b * IL1b_sp;

// ============================================================
// Clinical Outcomes — smoothed ODE representations
// ============================================================
// Pain NRS (0–10): driven by LTP, SP, offset by DPMS + drug effect
double pain_target = pain_base
                     + k_pain_LTP * LTP_cs * 6.0        // max +6
                     + k_pain_SP  * SP_csf * 2.0
                     - DPMS * 4.0                         // descending inhibition
                     - (inh_NET + inh_SERT) * 2.0;       // SNRI analgesic
pain_target = (pain_target < 0) ? 0 : (pain_target > 10 ? 10 : pain_target);
dxdt_Pain_score    = 0.15 * (pain_target - Pain_score);

// FIQR score (0–100): composite
double FIQ_target  = FIQ_base
                     + k_FIQ_pain  * (Pain_score - 5.0) * 8.0
                     + k_FIQ_sleep * (1.0 - SWS_depth) * 20.0
                     + k_FIQ_dep   * Depression_score * 0.5;
FIQ_target = (FIQ_target < 0) ? 0 : (FIQ_target > 100 ? 100 : FIQ_target);
dxdt_FIQ_score     = 0.10 * (FIQ_target - FIQ_score);

// Fatigue VAS (0–100)
double fatigue_tgt = fatigue_base
                     + k_fatigue * (1.0 - SWS_depth) * 30.0
                     - inh_NET * 20.0;    // NE up reduces fatigue
fatigue_tgt = (fatigue_tgt < 0) ? 0 : (fatigue_tgt > 100 ? 100 : fatigue_tgt);
dxdt_Fatigue_VAS   = 0.12 * (fatigue_tgt - Fatigue_VAS);

// Depression (PHQ-9 scaled 0–27)
double dep_target  = 8.0
                     + k_dep_LTP * LTP_cs * 12.0
                     + (1.0 - SWS_depth) * 4.0
                     - (inh_SERT + inh_NET) * 5.0;
dep_target = (dep_target < 0) ? 0 : (dep_target > 27 ? 27 : dep_target);
dxdt_Depression_score = 0.08 * (dep_target - Depression_score);

$TABLE
// C_DUL/C_PRE/C_MIL/C_TCA and EFFECT_* are already in scope here ($ODE and
// $TABLE compile into the same function in this mrgsolve build; redeclaring
// them under the same name a second time is a C++ redefinition error, not a
// scoping choice -- confirmed empirically, matching the precedent already
// established for this exact situation elsewhere in this corpus, e.g.
// acne-vulgaris). Only the two output-only percent metrics below are new
// here, preserved from the original including its own asymmetry (they omit
// the TCA contribution and the Emax_SNRI/EMAX_*_SERT multiplier that the
// live $ODE inh_SERT/inh_NET do include -- an original inconsistency kept
// unchanged for byte-identical verification; see fm_refactor_notes.md).
double inh_SERT_pct = 100 * (C_DUL/(EC50_DUL_SERT+C_DUL) + C_MIL/(EC50_MIL_SERT+C_MIL));
double inh_NET_pct  = 100 * (C_DUL/(EC50_DUL_NET+C_DUL)  + C_MIL/(EC50_MIL_NET+C_MIL));
double Ca_block_pct = 100 * EMAX_PRE * C_PRE/(EC50_PRE+C_PRE);

$CAPTURE @annotated
C_DUL   : Duloxetine plasma concentration (mg/L) -- canonical single site, redirectable
C_PRE   : Pregabalin plasma concentration (mg/L) -- canonical single site, redirectable
C_MIL   : Milnacipran plasma concentration (mg/L) -- canonical single site, redirectable
C_TCA   : Amitriptyline plasma concentration (mg/L) -- canonical single site, redirectable
inh_SERT_pct : SERT inhibition, DUL+MIL only, percent (matches original output metric, TCA omitted)
inh_NET_pct  : NET inhibition, DUL+MIL only, percent (matches original output metric, TCA omitted)
Ca_block_pct : Pregabalin alpha2-delta calcium-channel block, percent
EFFECT_DUL_SERT : Duloxetine SERT-inhibition Hill effect (0-1)
EFFECT_DUL_NET  : Duloxetine NET-inhibition Hill effect (0-1)
EFFECT_MIL_SERT : Milnacipran SERT-inhibition Hill effect (0-1)
EFFECT_MIL_NET  : Milnacipran NET-inhibition Hill effect (0-1)
EFFECT_TCA_SERT : Amitriptyline SERT-inhibition Hill effect (0-1)
EFFECT_TCA_NET  : Amitriptyline NET-inhibition Hill effect (0-1)
EFFECT_TCA_SEDATION : Amitriptyline sedative SWS-support Hill effect (0-1)
EFFECT_PRE      : Pregabalin alpha2-delta Hill effect (0-1)
'

## ============================================================
## Build and compile model
## ============================================================
fm_mod <- mcode("fibromyalgia_qsp", fm_code)

## ============================================================
## Initial conditions — FM steady state (untreated)
## ============================================================
FM_init <- init(fm_mod,
  NGF      = 1.5,    # elevated in FM
  PGE2     = 1.2,
  DRG_act  = 0.8,
  SP_csf   = 1.8,    # high CSF SP
  NMDA_state = 0.4,
  WindUp   = 0.3,
  LTP_cs   = 0.55,   # significant central sensitization
  NE_syn   = 0.7,
  SHT_syn  = 0.6,
  DPMS     = 0.35,   # blunted descending inhibition
  CRH      = 0.9,
  ACTH     = 0.85,
  CORT     = 0.8,
  SNS_tone = 0.75,   # elevated SNS
  Adenosine = 1.2,
  SWS_depth = 0.35,  # poor SWS
  MG_act   = 0.65,
  IL1b_sp  = 0.55,
  Pain_score   = 6.5,
  FIQ_score    = 68,
  Fatigue_VAS  = 72,
  Depression_score = 12
)

## ============================================================
## Dosing Events
## ============================================================
# Duloxetine 60 mg QD
dose_DUL_60 <- ev(amt = 60, ii = 24, addl = 83, cmt = "GUT_DUL",
                  param = list(use_DUL = 1))

# Pregabalin 150 mg BID
dose_PRE_150bid <- ev(amt = 150, ii = 12, addl = 167, cmt = "GUT_PRE",
                      param = list(use_PRE = 1))

# Milnacipran 50 mg BID
dose_MIL_50bid <- ev(amt = 50, ii = 12, addl = 167, cmt = "GUT_MIL",
                     param = list(use_MIL = 1))

# Amitriptyline 25 mg QHS
dose_TCA_25 <- ev(amt = 25, ii = 24, addl = 83, cmt = "GUT_TCA", time = 22,
                  param = list(use_TCA = 1))

# Combination: DUL + PRE
dose_combo <- ev_seq(dose_DUL_60, dose_PRE_150bid)

## ============================================================
## Simulation Scenarios
## ============================================================
sim_time <- seq(0, 84 * 24, by = 1)  # 12 weeks in hours

# Scenario 1: Untreated FM baseline
sim_base <- fm_mod %>%
  init(FM_init) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Untreated FM", time_days = time / 24)

# Scenario 2: Duloxetine 60 mg QD
sim_DUL <- fm_mod %>%
  init(FM_init) %>%
  param(use_DUL = 1) %>%
  ev(dose_DUL_60) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Duloxetine 60 mg QD", time_days = time / 24)

# Scenario 3: Pregabalin 150 mg BID
sim_PRE <- fm_mod %>%
  init(FM_init) %>%
  param(use_PRE = 1) %>%
  ev(dose_PRE_150bid) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Pregabalin 150mg BID", time_days = time / 24)

# Scenario 4: Milnacipran 50 mg BID
sim_MIL <- fm_mod %>%
  init(FM_init) %>%
  param(use_MIL = 1) %>%
  ev(dose_MIL_50bid) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Milnacipran 50mg BID", time_days = time / 24)

# Scenario 5: Duloxetine + Pregabalin combination
sim_COMBO <- fm_mod %>%
  init(FM_init) %>%
  param(use_DUL = 1, use_PRE = 1) %>%
  ev(dose_combo) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "DUL + PRE Combo", time_days = time / 24)

# Scenario 6: Low-dose Amitriptyline (sleep-targeted)
sim_TCA <- fm_mod %>%
  init(FM_init) %>%
  param(use_TCA = 1) %>%
  ev(dose_TCA_25) %>%
  mrgsim(end = 84 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Amitriptyline 25mg QHS", time_days = time / 24)

# Combine all scenarios
sim_all <- bind_rows(sim_base, sim_DUL, sim_PRE, sim_MIL, sim_COMBO, sim_TCA)

## ============================================================
## Visualization
## ============================================================
library(ggplot2)
library(dplyr)
library(tidyr)

colors_scen <- c(
  "Untreated FM"         = "#E74C3C",
  "Duloxetine 60 mg QD"  = "#3498DB",
  "Pregabalin 150mg BID" = "#2ECC71",
  "Milnacipran 50mg BID" = "#9B59B6",
  "DUL + PRE Combo"      = "#E67E22",
  "Amitriptyline 25mg QHS" = "#1ABC9C"
)

# Plot 1: Pain score over time
p1 <- ggplot(sim_all, aes(time_days, Pain_score, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "NRS Pain Score (0-10) — 12-Week Treatment",
       x = "Time (days)", y = "Pain NRS", color = "Scenario") +
  theme_bw(14) + ylim(0, 10) +
  geom_hline(yintercept = c(3, 5), linetype = "dashed", color = "grey50")

# Plot 2: FIQ score
p2 <- ggplot(sim_all, aes(time_days, FIQ_score, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "FIQR Score (0-100) — Functional Impact",
       x = "Time (days)", y = "FIQR", color = "Scenario") +
  theme_bw(14)

# Plot 3: SWS depth and fatigue
p3 <- ggplot(sim_all, aes(time_days, SWS_depth, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Slow-Wave Sleep Depth (normalized)",
       x = "Time (days)", y = "SWS Depth (a.u.)", color = "Scenario") +
  theme_bw(14)

# Plot 4: Central sensitization (LTP)
p4 <- ggplot(sim_all, aes(time_days, LTP_cs, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Central Sensitization Index (Spinal LTP)",
       x = "Time (days)", y = "LTP_cs (0-1)", color = "Scenario") +
  theme_bw(14)

# Plot 5: PK — duloxetine plasma concentration
p5 <- ggplot(filter(sim_DUL, time_days <= 14),
             aes(time_days, C_DUL)) +
  geom_line(color = "#3498DB", linewidth = 1) +
  labs(title = "Duloxetine PK — Plasma Concentration (first 14 days)",
       x = "Time (days)", y = "Cp (mg/L)") +
  theme_bw(14)

# Plot 6: PD biomarker — CSF Substance P
p6 <- ggplot(sim_all, aes(time_days, SP_csf, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = colors_scen) +
  labs(title = "CSF Substance P (biomarker of central sensitization)",
       x = "Time (days)", y = "SP_csf (a.u.)", color = "Scenario") +
  theme_bw(14)

# Print summary at week 12
summary_wk12 <- sim_all %>%
  filter(time_days >= 83 & time_days <= 84) %>%
  group_by(scenario) %>%
  summarise(
    Pain_NRS     = round(mean(Pain_score), 2),
    FIQR         = round(mean(FIQ_score), 1),
    Fatigue_VAS  = round(mean(Fatigue_VAS), 1),
    Depression   = round(mean(Depression_score), 1),
    SWS_depth    = round(mean(SWS_depth), 3),
    LTP_cs       = round(mean(LTP_cs), 3),
    SP_csf       = round(mean(SP_csf), 3),
    DPMS         = round(mean(DPMS), 3),
    .groups = "drop"
  ) %>%
  arrange(Pain_NRS)

cat("\n=== FM QSP Model — 12-week Treatment Outcome Summary ===\n")
print(summary_wk12)

## ============================================================
## Responder analysis: >=30% pain reduction
## ============================================================
baseline_pain <- mean(filter(sim_base, time_days <= 1)$Pain_score)

resp_30 <- sim_all %>%
  filter(time_days >= 83) %>%
  group_by(scenario) %>%
  summarise(
    pain_wk12 = mean(Pain_score),
    pct_change = 100 * (pain_wk12 - baseline_pain) / baseline_pain,
    responder_30 = pct_change <= -30,
    responder_50 = pct_change <= -50,
    .groups = "drop"
  )

cat("\n=== Responder Analysis ===\n")
print(resp_30)

## ============================================================
## Output plots
## ============================================================
library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)
grid.arrange(p5, p6, nrow = 1)

message("FM QSP model simulation complete.")
