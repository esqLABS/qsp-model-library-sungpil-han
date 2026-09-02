## ============================================================
## Behcet's Disease (BD) — QSP mrgsolve Model — REFACTORED (pluggable PK)
## ============================================================
## Description:
##   Fork refactor of bd_mrgsolve_model.R (never edited — see
##   bd_refactor_notes.md). Same disease biology (neutrophil activation,
##   Th1/Th17 balance, cytokine network, vascular inflammation, multi-organ
##   clinical endpoints) and the same five drugs (colchicine, prednisolone,
##   adalimumab, apremilast, canakinumab), reorganized so each compound's
##   PK lives in its own named block, exposes exactly one concentration
##   (C_<STEM>) and one or more named Hill effects (EFFECT_<STEM>_*), per
##   the fork's naming convention. No disease-side equation, baseline
##   parameter, or numeric value is changed.
##
## Refactor summary (see bd_refactor_notes.md for full detail):
##   - Colchicine  (COL)  : archetype 3 (depot + central + peripheral)
##   - Prednisolone(PRED) : archetype 2 (central + peripheral, no depot)
##   - Adalimumab  (ADA)  : archetype 2 (central + peripheral, no depot)
##   - Apremilast  (APR)  : archetype 1 (single compartment, no depot)
##   - Canakinumab (CAN)  : archetype 2 (central + peripheral, no depot)
##   - The original's own $MAIN (PD-facing) and $TABLE (capture-facing)
##     each computed the same concentration formula under a different
##     name (Cp_col vs Cp_COL, etc.) -- textually a duplicate, but NOT a
##     numeric one: mrgsolve evaluates $MAIN with interval-start state and
##     $TABLE with interval-end state, so the two differ by one reporting
##     step. Verification caught this when an initial attempt collapsed
##     them into one value, introducing a reproducible one-row lag (see
##     bd_refactor_notes.md). Resolved by keeping each computed exactly
##     where the original had it: C_<STEM> as a single, canonical
##     `double C_<STEM> = <expr>;` in $MAIN (the value the Hill/
##     EFFECT_<STEM> interface reads, PD-facing, matching the original's
##     Cp_col/Cp_pred/etc. placement exactly), and $TABLE keeping its own
##     independent, fresh recomputation for the Cp_<STEM> reporting
##     columns -- same as the original did, just renamed.
##   - Neither Adalimumab nor Canakinumab models receptor-binding kinetics
##     in the original (both are a plain Emax*C/(EC50+C) ratio on plasma
##     concentration, no receptor/complex compartments) — so neither
##     needed archetype 4 (TMDD); confirmed from the code, not assumed.
##
## Author: QSP Auto-generator (Claude Code Routine) — fork refactor pass
## Date:   2026-09-02
## ============================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ============================================================
## Model Code
## ============================================================
bd_code_refactored <- '
$PROB
Behcet Disease QSP Model — PK/PD with Immune Dynamics (refactored, pluggable PK)
20 ODE compartments (8 PK + 12 PD)

$PARAM @annotated
// ---- Colchicine PK (archetype 3: depot + central + peripheral) ----
KA_COL   : 1.2   : Colchicine absorption rate constant (1/h)
F_COL    : 0.44  : Colchicine oral bioavailability (fraction)
CL_COL   : 17.0  : Colchicine clearance (L/h)
V1_COL   : 28.0  : Colchicine central volume (L)
V2_COL   : 5500  : Colchicine tissue volume (L) — high Vd=21L/kg
Q_COL    : 30.0  : Colchicine inter-compartmental CL (L/h)

// ---- Prednisolone PK (archetype 2: central + peripheral, no depot) ----
// KA_PRED/F_PRED are carried over unused, exactly as in the original: the
// original declares an absorption rate and oral bioavailability but never
// gives prednisolone a depot compartment, dosing straight into the
// central compartment via a rate=-2 flag instead — see refactor notes.
CL_PRED  : 8.4   : Prednisolone CL (L/h)
V1_PRED  : 22.0  : Prednisolone central volume (L)
V2_PRED  : 50.0  : Prednisolone tissue volume (L)
Q_PRED   : 10.0  : Prednisolone intercomp CL (L/h)
F_PRED   : 0.80  : Prednisolone bioavailability (declared, unused — see notes)
KA_PRED  : 1.0   : Prednisolone absorption rate (1/h) (declared, unused — see notes)

// ---- Adalimumab (anti-TNF) PK (archetype 2: central + peripheral, no depot) ----
// F_ADA is likewise carried over unused — no depot compartment models the
// SC absorption it would apply to; dosed straight into central via rate=-2.
CL_ADA   : 0.012 : Adalimumab CL (L/h) — t1/2 ~2 wk
V1_ADA   : 2.8   : Adalimumab central volume (L)
V2_ADA   : 3.4   : Adalimumab peripheral volume (L)
Q_ADA    : 0.003 : Adalimumab intercomp CL (L/h)
F_ADA    : 0.64  : Adalimumab SC bioavailability (declared, unused — see notes)

// ---- Apremilast PK (archetype 1: single compartment, no depot) ----
// KA_APR/F_APR are carried over unused — the original declares an oral
// absorption rate and bioavailability but doses AAPR (now CENT_APR)
// directly as an instantaneous bolus with no depot in $ODE at all.
KA_APR   : 0.58  : Apremilast absorption rate (1/h) (declared, unused — see notes)
CL_APR   : 10.0  : Apremilast clearance (L/h)
V1_APR   : 87.0  : Apremilast volume of distribution (L)
F_APR    : 0.73  : Apremilast bioavailability (declared, unused — see notes)

// ---- Canakinumab PK (archetype 2: central + peripheral, no depot) ----
// F_CAN is likewise carried over unused, same pattern as ADA above.
CL_CAN   : 0.007 : Canakinumab CL (L/h) — t1/2 ~26 d
V1_CAN   : 3.0   : Canakinumab central volume (L)
V2_CAN   : 3.2   : Canakinumab peripheral volume (L)
Q_CAN    : 0.0015: Canakinumab intercomp CL (L/h)
F_CAN    : 0.66  : Canakinumab SC bioavailability (declared, unused — see notes)

// ---- Disease Baseline Parameters ----
NEU0     : 1.0   : Baseline neutrophil activation (normalized)
TH1_BASE0  : 1.0   : Baseline Th1 cells (normalized) (renamed from TH1_0 -- build-compatibility fix, see notes)
TH17_BASE0 : 1.0   : Baseline Th17 cells (normalized) (renamed from TH17_0 -- build-compatibility fix, see notes)
TREG_BASE0 : 1.0   : Baseline Treg cells (normalized) (renamed from TREG_0 -- build-compatibility fix, see notes)
TNFA0    : 1.0   : Baseline TNF-alpha (normalized)
IL1B0    : 1.0   : Baseline IL-1beta (normalized)
IL6_0    : 1.0   : Baseline IL-6 (normalized)
IL17A0   : 1.0   : Baseline IL-17A (normalized)
EA0      : 0.5   : Baseline endothelial activation (normalized)

// ---- Disease Kinetic Rate Constants ----
kNEU_in  : 0.05  : Neutrophil activation synthesis rate (1/h)
kNEU_out : 0.05  : Neutrophil activation decay rate (1/h)
kTH1_in  : 0.03  : Th1 proliferation/activation rate (1/h)
kTH1_out : 0.03  : Th1 decay rate (1/h)
kTH17_in : 0.03  : Th17 proliferation rate (1/h)
kTH17_out: 0.03  : Th17 decay rate (1/h)
kTREG_in : 0.02  : Treg generation rate (1/h)
kTREG_out: 0.02  : Treg decay rate (1/h)
kTNFA_syn: 0.08  : TNF-alpha synthesis rate (1/h)
kTNFA_deg: 0.08  : TNF-alpha degradation rate (1/h)
kIL1B_syn: 0.06  : IL-1beta synthesis rate (1/h)
kIL1B_deg: 0.06  : IL-1beta degradation rate (1/h)
kIL6_syn : 0.07  : IL-6 synthesis rate (1/h)
kIL6_deg : 0.07  : IL-6 degradation rate (1/h)
kIL17_syn: 0.05  : IL-17A synthesis rate (1/h)
kIL17_deg: 0.05  : IL-17A decay rate (1/h)
kEA_on   : 0.04  : Endothelial activation rate (1/h)
kEA_off  : 0.04  : Endothelial deactivation rate (1/h)

// ---- Cytokine Cross-talk Amplification ----
a_TNFA_TH1 : 0.3  : TNF-alpha amplification of Th1
a_TNFA_NEU : 0.4  : TNF-alpha priming of neutrophils
a_IL17_NEU : 0.3  : IL-17A driving neutrophil recruitment
a_TH17_IL17: 0.5  : Th17 → IL-17A production
a_TH1_TNFA : 0.4  : Th1 → TNF-alpha production
a_NEU_TNFA : 0.3  : Neutrophil → TNF-alpha
a_NEU_IL1B : 0.3  : Neutrophil → IL-1B
a_IL1B_IL6 : 0.3  : IL-1beta driving IL-6
a_TNFA_IL6 : 0.2  : TNF-alpha driving IL-6
a_IL6_TH17 : 0.2  : IL-6 driving Th17 polarization
a_TREG_inh : 0.3  : Treg suppression of Th17/Th1

// ---- Organ Manifestation Parameters ----
kOUL_on  : 0.02  : Oral ulcer onset rate (/h)
kOUL_off : 0.01  : Oral ulcer healing rate (/h)
kOCI_on  : 0.015 : Ocular inflammation onset rate (/h)
kOCI_off : 0.008 : Ocular inflammation resolution (/h)
kEA_vul  : 0.03  : EA driving vascular/ocular pathology (/h)

// ---- Drug PD Potency Parameters (Hill interface) ----
// EMAX_/GAMMA_ pulled out explicit — the original’s effect terms were
// already a plain ratio (colchicine: no Emax at all, i.e. Emax=1) or a
// plain Emax*C/(EC50+C) (n=1); this is a rename, not a refit.
EC50_COL_NEU  : 0.8  : Colchicine EC50 on neutrophil migration (ng/mL)
EMAX_COL_NEU  : 1.0  : Colchicine max effect on neutrophil migration (original had no Emax term, i.e. Emax=1)
GAMMA_COL_NEU : 1.0  : Colchicine Hill coefficient, neutrophil effect (original had no explicit Hill term)
EC50_COL_IL1  : 1.5  : Colchicine EC50 on IL-1B via NLRP3 (ng/mL)
EMAX_COL_IL1  : 1.0  : Colchicine max effect on IL-1B (original had no Emax term, i.e. Emax=1)
GAMMA_COL_IL1 : 1.0  : Colchicine Hill coefficient, IL-1B effect (original had no explicit Hill term)
EC50_PRED_TNF : 50.0 : Prednisolone EC50 for TNF-alpha suppression (ng/mL)
EC50_PRED_IL6 : 30.0 : Prednisolone EC50 for IL-6 suppression (ng/mL)
EMAX_PRED     : 0.80 : Max effect of prednisolone (shared by TNF and IL-6 effects, as in the original)
GAMMA_PRED_TNF: 1.0  : Prednisolone Hill coefficient, TNF-alpha effect (original had no explicit Hill term)
GAMMA_PRED_IL6: 1.0  : Prednisolone Hill coefficient, IL-6 effect (original had no explicit Hill term)
EC50_ADA_TNF  : 1500 : Adalimumab EC50 for TNF neutralization (ng/mL)
EMAX_ADA      : 0.90 : Max effect of adalimumab on TNF-alpha
GAMMA_ADA_TNF : 1.0  : Adalimumab Hill coefficient, TNF-alpha effect (original had no explicit Hill term)
EC50_APR_TNF  : 200  : Apremilast EC50 for TNF-alpha (ng/mL)
EC50_APR_IL17 : 300  : Apremilast EC50 for IL-17A (ng/mL)
EMAX_APR      : 0.65 : Max Emax of apremilast (shared by TNF and IL-17A effects, as in the original)
GAMMA_APR_TNF : 1.0  : Apremilast Hill coefficient, TNF-alpha effect (original had no explicit Hill term)
GAMMA_APR_IL17: 1.0  : Apremilast Hill coefficient, IL-17A effect (original had no explicit Hill term)
EC50_CAN_IL1  : 800  : Canakinumab EC50 for IL-1beta (ng/mL)
EMAX_CAN      : 0.92 : Canakinumab max effect on IL-1beta
GAMMA_CAN_IL1 : 1.0  : Canakinumab Hill coefficient, IL-1beta effect (original had no explicit Hill term)

// ---- Disease Severity (HLA-B51 effect) ----
HLAB51_factor: 1.4  : HLA-B51 positive multiplier on disease severity

$CMT @annotated
// PK Compartments (same 10 slots, same order, as the original — cmt
// index-based dosing events from the original file work unchanged)
GUT_COL  : Colchicine gut depot (mg)
CENT_COL : Colchicine central (mg)
PERI_COL : Colchicine tissue (mg)
CENT_PRED: Prednisolone central (mg)
PERI_PRED: Prednisolone tissue (mg)
CENT_ADA : Adalimumab central (mg)
PERI_ADA : Adalimumab peripheral (mg)
CENT_APR : Apremilast central (mg)
CENT_CAN : Canakinumab central (mg)
PERI_CAN : Canakinumab peripheral (mg)

// PD Compartments — Immune Cells
NEU      : Neutrophil activation state (normalized)
TH1      : Th1 cell activity (normalized)
TH17     : Th17 cell activity (normalized)
TREG     : Regulatory T cell activity (normalized)

// PD Compartments — Cytokines
TNFA     : TNF-alpha (normalized)
IL1B     : IL-1beta (normalized)
IL6C     : IL-6 (normalized)
IL17A    : IL-17A (normalized)

// PD Compartments — Organ Manifestations
EA       : Endothelial activation (normalized)
OUL      : Oral ulcer activity index
OCI      : Ocular inflammation index
BDCAF    : BDCAF composite score

$MAIN
// --- Single exposed concentration per compound (ng/mL), PD-facing ---
// Pluggable-PK interface (fork refactor — see behcet-disease/
// bd_refactor_notes.md). Each is a single, contiguous `double C_<STEM> =
// <expr>;` initializing statement (never a bare forward declaration
// assigned later — that pattern is invisible to this corpus’s own
// downstream discoverability tooling, which pattern-matches the source
// text directly; see the guide’s "discoverable by downstream tooling"
// section). This is the one canonical site each compound’s concentration
// is computed at for the disease-effect (EFFECT_<STEM>) calculations
// below, matching exactly where and how the original computed
// Cp_col/Cp_pred/etc. for the same purpose in its own $MAIN (not moved
// into $ODE — see the guide’s warning on why that changes the simulated
// trajectory). $TABLE’s own Cp_<STEM> reporting columns are a separate,
// independent recomputation, not a reference to these — see the $TABLE
// block below for why (verification-caught reporting-lag finding,
// bd_refactor_notes.md); no $GLOBAL indirection is needed since nothing
// outside $MAIN/$ODE reads these.
double C_COL  = CENT_COL  / V1_COL  * 1000;
double C_PRED = CENT_PRED / V1_PRED * 1000;
double C_ADA  = CENT_ADA  / V1_ADA  * 1000;
double C_APR  = CENT_APR  / V1_APR  * 1000;
double C_CAN  = CENT_CAN  / V1_CAN  * 1000;

// === Drug Effect Functions (Hill equation interface) ===
// Colchicine effects (rename only — original was Cp/(IC50+Cp), i.e.
// Emax=1, n=1)
double EFFECT_COL_NEU = EMAX_COL_NEU * pow(C_COL, GAMMA_COL_NEU)
                 / (pow(EC50_COL_NEU, GAMMA_COL_NEU) + pow(C_COL, GAMMA_COL_NEU));
double EFFECT_COL_IL1 = EMAX_COL_IL1 * pow(C_COL, GAMMA_COL_IL1)
                 / (pow(EC50_COL_IL1, GAMMA_COL_IL1) + pow(C_COL, GAMMA_COL_IL1));

// Prednisolone effects (rename only — original was Emax*Cp/(EC50+Cp))
double EFFECT_PRED_TNF = EMAX_PRED * pow(C_PRED, GAMMA_PRED_TNF)
                  / (pow(EC50_PRED_TNF, GAMMA_PRED_TNF) + pow(C_PRED, GAMMA_PRED_TNF));
double EFFECT_PRED_IL6 = EMAX_PRED * pow(C_PRED, GAMMA_PRED_IL6)
                  / (pow(EC50_PRED_IL6, GAMMA_PRED_IL6) + pow(C_PRED, GAMMA_PRED_IL6));

// Adalimumab effects (anti-TNF; rename only — plain Emax*Cp/(EC50+Cp),
// no receptor-binding ODE in the original, so no TMDD here)
double EFFECT_ADA_TNF = EMAX_ADA * pow(C_ADA, GAMMA_ADA_TNF)
                 / (pow(EC50_ADA_TNF, GAMMA_ADA_TNF) + pow(C_ADA, GAMMA_ADA_TNF));

// Apremilast effects (PDE4 inhibitor → cAMP ↑; rename only)
double EFFECT_APR_TNF  = EMAX_APR * pow(C_APR, GAMMA_APR_TNF)
                  / (pow(EC50_APR_TNF, GAMMA_APR_TNF) + pow(C_APR, GAMMA_APR_TNF));
double EFFECT_APR_IL17 = EMAX_APR * pow(C_APR, GAMMA_APR_IL17)
                  / (pow(EC50_APR_IL17, GAMMA_APR_IL17) + pow(C_APR, GAMMA_APR_IL17));

// Canakinumab effects (anti-IL-1B; rename only — plain Emax*Cp/(EC50+Cp),
// no receptor-binding ODE in the original, so no TMDD here)
double EFFECT_CAN_IL1 = EMAX_CAN * pow(C_CAN, GAMMA_CAN_IL1)
                 / (pow(EC50_CAN_IL1, GAMMA_CAN_IL1) + pow(C_CAN, GAMMA_CAN_IL1));

// === Baseline disease synthesis rates with HLA-B51 === (unchanged)
double DIS = HLAB51_factor; // disease severity multiplier

// Cytokine network drivers
// TNF-alpha synthesis: from Th1, neutrophils, macrophages
double TNFA_syn = kTNFA_syn * DIS * (1.0 + a_TH1_TNFA*(TH1 - 1.0) + a_NEU_TNFA*(NEU - 1.0));

// IL-1B synthesis: NLRP3 activation (driven by NEU, EA)
double IL1B_syn = kIL1B_syn * DIS * (1.0 + a_NEU_IL1B*(NEU - 1.0));

// IL-6 synthesis: driven by IL-1B and TNF-alpha
double IL6_syn  = kIL6_syn  * DIS * (1.0 + a_IL1B_IL6*(IL1B - 1.0) + a_TNFA_IL6*(TNFA - 1.0));

// IL-17A synthesis: from Th17 cells
double IL17_syn = kIL17_syn * DIS * (1.0 + a_TH17_IL17*(TH17 - 1.0));

// Neutrophil activation: driven by IL-8 (proxy=TNFA + IL17A), suppressed by colchicine
double NEU_syn  = kNEU_in * DIS * (1.0 + a_TNFA_NEU*(TNFA - 1.0) + a_IL17_NEU*(IL17A - 1.0));

// Th1 polarization: driven by TNF-alpha, IL-12
double TH1_syn  = kTH1_in * DIS * (1.0 + a_TNFA_TH1*(TNFA - 1.0));

// Th17 polarization: driven by IL-6, IL-23; suppressed by Treg
double TH17_syn = kTH17_in * DIS * (1.0 + a_IL6_TH17*(IL6C - 1.0)) / (1.0 + a_TREG_inh*TREG);

// Treg: induced by TGF-beta (inverse of disease inflammation)
double TREG_syn = kTREG_in / (1.0 + 0.2*(TNFA - 1.0));

// Endothelial activation: driven by TNF-alpha, IL-1B, IL-17A
double EA_on_rate = kEA_on * DIS * (TNFA + IL1B + IL17A) / 3.0;

// Oral ulcer onset: driven by IL-17A, TNFA, neutrophils
double OUL_on_rate = kOUL_on * DIS * (IL17A*0.4 + TNFA*0.4 + NEU*0.2);

// Ocular inflammation: driven by endothelial activation
double OCI_on_rate = kOCI_on * DIS * EA;

$ODE
// ---- COLCHICINE PK ----
dxdt_GUT_COL  = -KA_COL * GUT_COL;
dxdt_CENT_COL =  KA_COL * GUT_COL * F_COL
                - (CL_COL + Q_COL) / V1_COL * CENT_COL
                + Q_COL / V2_COL * PERI_COL;
dxdt_PERI_COL =  Q_COL / V1_COL * CENT_COL
                - Q_COL / V2_COL * PERI_COL;

// ---- PREDNISOLONE PK (2-CMT, no depot — dosed via RATE, as original) ----
dxdt_CENT_PRED =  -(CL_PRED + Q_PRED) / V1_PRED * CENT_PRED
                  + Q_PRED / V2_PRED * PERI_PRED;
dxdt_PERI_PRED =  Q_PRED / V1_PRED * CENT_PRED
                 - Q_PRED / V2_PRED * PERI_PRED;

// ---- ADALIMUMAB PK ----
dxdt_CENT_ADA =  -(CL_ADA + Q_ADA) / V1_ADA * CENT_ADA
                 + Q_ADA / V2_ADA * PERI_ADA;
dxdt_PERI_ADA =  Q_ADA / V1_ADA * CENT_ADA
                - Q_ADA / V2_ADA * PERI_ADA;

// ---- APREMILAST PK ----
dxdt_CENT_APR = -CL_APR / V1_APR * CENT_APR;

// ---- CANAKINUMAB PK ----
dxdt_CENT_CAN =  -(CL_CAN + Q_CAN) / V1_CAN * CENT_CAN
                 + Q_CAN / V2_CAN * PERI_CAN;
dxdt_PERI_CAN =  Q_CAN / V1_CAN * CENT_CAN
                - Q_CAN / V2_CAN * PERI_CAN;

// ---- NEUTROPHIL ACTIVATION ----
// Drug effects: colchicine inhibits neutrophil migration
dxdt_NEU = NEU_syn - kNEU_out * NEU * (1.0 + EFFECT_COL_NEU);

// ---- Th1 CELLS ----
// Drug effects: prednisolone
dxdt_TH1 = TH1_syn - kTH1_out * TH1 * (1.0 + EFFECT_PRED_TNF*0.5);

// ---- Th17 CELLS ----
// Drug effects: prednisolone, apremilast (via cAMP)
dxdt_TH17 = TH17_syn - kTH17_out * TH17 * (1.0 + EFFECT_PRED_TNF*0.3 + EFFECT_APR_IL17*0.4);

// ---- TREG CELLS ----
// Treg induced over time as inflammation resolves
dxdt_TREG = TREG_syn - kTREG_out * TREG;

// ---- TNF-alpha ----
// Drug effects: adalimumab neutralizes TNF-alpha; pred reduces synthesis; apremilast PDE4 pathway
dxdt_TNFA = TNFA_syn * (1.0 - EFFECT_ADA_TNF) * (1.0 - EFFECT_PRED_TNF) * (1.0 - EFFECT_APR_TNF*0.6)
            - kTNFA_deg * TNFA;

// ---- IL-1beta ----
// Drug effects: colchicine (NLRP3 inhibition), canakinumab (neutralization), prednisolone
dxdt_IL1B = IL1B_syn * (1.0 - EFFECT_COL_IL1) * (1.0 - EFFECT_CAN_IL1) * (1.0 - EFFECT_PRED_TNF*0.4)
            - kIL1B_deg * IL1B;

// ---- IL-6 ----
// Drug effects: prednisolone reduces IL-6 synthesis
dxdt_IL6C = IL6_syn * (1.0 - EFFECT_PRED_IL6)
            - kIL6_deg * IL6C;

// ---- IL-17A ----
// Drug effects: apremilast (PDE4 inh), prednisolone
dxdt_IL17A = IL17_syn * (1.0 - EFFECT_APR_IL17) * (1.0 - EFFECT_PRED_TNF*0.3)
             - kIL17_deg * IL17A;

// ---- ENDOTHELIAL ACTIVATION ----
// Driven by TNFA, IL1B, IL17A; suppressed by treatment
dxdt_EA = EA_on_rate * (1.0 - EFFECT_ADA_TNF*0.5 - EFFECT_PRED_TNF*0.3) - kEA_off * EA;

// ---- ORAL ULCER ACTIVITY INDEX ----
// Onset driven by cytokines/neutrophils; healing suppressed by high EA
dxdt_OUL = OUL_on_rate * (1.0 - EFFECT_PRED_TNF*0.4 - EFFECT_ADA_TNF*0.4 - EFFECT_APR_TNF*0.3)
           - kOUL_off * OUL;

// ---- OCULAR INFLAMMATION INDEX ----
// Driven by endothelial activation; responds well to anti-TNF/steroid
dxdt_OCI = OCI_on_rate * (1.0 - EFFECT_ADA_TNF*0.6 - EFFECT_PRED_TNF*0.5)
           - kOCI_off * OCI;

// ---- BDCAF COMPOSITE SCORE ----
// BDCAF = weighted sum of organ manifestations (dynamic)
dxdt_BDCAF = 0.25 * OUL + 0.25 * OCI + 0.25 * EA + 0.25 * (TNFA + IL1B + IL6C + IL17A)/4.0
             - 0.1 * BDCAF;

$TABLE
// NOT redirected to $MAIN’s C_<STEM> -- verification caught a real,
// reproducible one-row reporting lag when this was first tried as a
// straight alias (see bd_refactor_notes.md for the full trace). mrgsolve
// calls $MAIN once per record using the compartment state as of the
// *start* of that record’s interval (the same "stale" value the disease
// $ODE’s Hill effects deliberately act on all interval), then advances the
// ODE, then calls $TABLE using the state as of the *end* of that interval
// (after any dose at this exact record has been applied) -- so the
// original’s own $TABLE recomputation was not a true duplicate of $MAIN’s,
// despite using the identical formula: same formula, two different points
// in time. Collapsing them (aliasing here to $MAIN’s C_<STEM>) reproduced
// exactly this lag when checked against the original via the qspserver
// API. Kept as the original had it: an independent, fresh, end-of-
// interval recomputation from compartment state, purely for reporting.
double Cp_COL  = CENT_COL  / V1_COL  * 1000;  // ng/mL
double Cp_PRED = CENT_PRED / V1_PRED * 1000;  // ng/mL
double Cp_ADA  = CENT_ADA  / V1_ADA  * 1000;  // ng/mL
double Cp_APR  = CENT_APR  / V1_APR  * 1000;  // ng/mL
double Cp_CAN  = CENT_CAN  / V1_CAN  * 1000;  // ng/mL

// Derived biomarkers
double Oral_Ulcer_Score    = OUL;
double Ocular_Inflam_Score = OCI;
double Disease_Activity    = BDCAF;
double Endoth_Activation   = EA;
double TNFA_level          = TNFA;
double IL1B_level          = IL1B;
double IL6_level           = IL6C;
double IL17A_level         = IL17A;
double Neutrophil_Act      = NEU;
double Th17_Activity       = TH17;
double Th1_Activity        = TH1;

$CAPTURE @annotated
Cp_COL  : Colchicine plasma concentration (ng/mL)
Cp_PRED : Prednisolone plasma concentration (ng/mL)
Cp_ADA  : Adalimumab plasma concentration (ng/mL)
Cp_APR  : Apremilast plasma concentration (ng/mL)
Cp_CAN  : Canakinumab plasma concentration (ng/mL)
C_COL   : Colchicine — single exposed concentration (ng/mL)
C_PRED  : Prednisolone — single exposed concentration (ng/mL)
C_ADA   : Adalimumab — single exposed concentration (ng/mL)
C_APR   : Apremilast — single exposed concentration (ng/mL)
C_CAN   : Canakinumab — single exposed concentration (ng/mL)
EFFECT_COL_NEU  : Colchicine effect on neutrophil migration (fraction)
EFFECT_COL_IL1  : Colchicine effect on IL-1B/NLRP3 (fraction)
EFFECT_PRED_TNF : Prednisolone effect on TNF-alpha (fraction)
EFFECT_PRED_IL6 : Prednisolone effect on IL-6 (fraction)
EFFECT_ADA_TNF  : Adalimumab effect on TNF-alpha (fraction)
EFFECT_APR_TNF  : Apremilast effect on TNF-alpha (fraction)
EFFECT_APR_IL17 : Apremilast effect on IL-17A (fraction)
EFFECT_CAN_IL1  : Canakinumab effect on IL-1beta (fraction)
'

## NOTE (fork build-compatibility fix -- see bd_refactor_notes.md and
## UPSTREAM_ISSUES.md): the original's own $CAPTURE @annotated block also
## listed NEU, TH1, TH17, TREG, TNFA, IL1B, IL6C, IL17A, EA, OUL, OCI,
## BDCAF -- but every one of those is already a $CMT compartment, and
## mrgsolve 2.0.1 rejects a $CAPTURE entry that duplicates a compartment
## name ("compartment should not be in $CAPTURE"). This is a pre-existing
## defect in the original, confirmed by feeding the original's own DSL
## (unmodified) through the same qspserver mrgsolve_api /model_manifest
## call: it fails with the identical error. Removing the 12 duplicated
## names from $CAPTURE here changes nothing numeric -- mrgsolve reports
## every $CMT compartment in its output by default whether or not it is
## also named in $CAPTURE, so NEU/TH1/.../BDCAF remain in every run's
## output columns exactly as before.
##
## A second, independent build-compatibility fix was needed once the
## above was fixed: the original's own $PARAM block declares TH1_0,
## TH17_0, and TREG_0 (baseline placeholders, never referenced by any
## $MAIN/$ODE/$TABLE calculation -- confirmed by grep) with names that
## collide exactly with mrgsolve's own auto-generated "<CMT>_0" initial-
## value symbol for compartments TH1, TH17, and TREG, producing a C++
## "conflicting declaration" compile error. Every sibling baseline param
## in the same block (NEU0, TNFA0, IL1B0, IL6_0, IL17A0, EA0) already
## avoids this by omitting the underscore or using a different compartment
## name -- only these three kept the colliding underscored form. Renamed
## to TH1_BASE0/TH17_BASE0/TREG_BASE0 here; since none of the three is
## read anywhere, this is a pure syntax fix with no numeric effect.

## ============================================================
## Compile Model
## ============================================================
bd_mod_refactored <- mcode("BehcetDisease_QSP_refactored", bd_code_refactored)

## ============================================================
## Initial Conditions (Active Behcet's Disease, HLA-B51+)
## Same values as the original, compartment names renamed to match.
## ============================================================
bd_init_refactored <- list(
  GUT_COL = 0,
  CENT_COL = 0, PERI_COL = 0,
  CENT_PRED = 0, PERI_PRED = 0,
  CENT_ADA = 0, PERI_ADA = 0,
  CENT_APR = 0,
  CENT_CAN = 0, PERI_CAN = 0,
  NEU   = 2.5,   # Elevated neutrophil activation (active disease)
  TH1   = 2.0,   # Elevated Th1
  TH17  = 2.8,   # Markedly elevated Th17 (Th17-dominant BD)
  TREG  = 0.5,   # Reduced Treg (immune dysregulation)
  TNFA  = 3.0,   # Elevated TNF-alpha
  IL1B  = 2.5,   # Elevated IL-1beta
  IL6C  = 2.0,   # Elevated IL-6
  IL17A = 2.8,   # Elevated IL-17A
  EA    = 2.0,   # Elevated endothelial activation
  OUL   = 3.0,   # Active oral ulcers
  OCI   = 2.0,   # Active ocular inflammation
  BDCAF = 8.0    # Moderate-severe BDCAF (scale 0-12)
)

## ============================================================
## Treatment Scenarios
## Identical to the original: same numeric cmt indices work unchanged
## because the refactored $CMT block preserves the exact same 10-slot PK
## ordering (see bd_refactor_notes.md).
## ============================================================
time_grid <- seq(0, 2160, by = 4)  # 90 days, 4-hour intervals

## Scenario 1: Untreated (disease natural course)
ev_untreated <- ev(cmt = 1, amt = 0, time = 0)

## Scenario 2: Colchicine monotherapy (0.5 mg BID = 1 mg/day)
ev_colch <- ev(
  data.frame(
    time = c(seq(0, 2156, by = 12)),
    cmt  = 1,   # GUT_COL
    amt  = 0.5, # mg per dose
    evid = 1
  )
)

## Scenario 3: Prednisolone monotherapy (40 mg/day, oral)
ev_pred <- ev(
  data.frame(
    time = c(seq(0, 2156, by = 24)),
    cmt  = 4,   # CENT_PRED
    amt  = 40,  # mg/day
    evid = 1,
    rate = -2   # immediate absorption (infusion rate flag for oral)
  )
)

## Scenario 4: Adalimumab (anti-TNF, SC 40 mg Q2W)
ev_ada <- ev(
  data.frame(
    time = seq(0, 2160, by = 336),  # Q2W = every 14 days = 336 h
    cmt  = 6,    # CENT_ADA
    amt  = 40,   # mg
    evid = 1,
    rate = -2
  )
)

## Scenario 5: Apremilast (30 mg BID, oral)
ev_apremilast <- ev(
  data.frame(
    time = c(seq(0, 2156, by = 12)),
    cmt  = 8,    # CENT_APR
    amt  = 30,   # mg per dose
    evid = 1
  )
)

## Scenario 6: Canakinumab (SC 150 mg Q8W)
ev_can <- ev(
  data.frame(
    time = seq(0, 2160, by = 1344),  # Q8W = 8*7*24 = 1344 h
    cmt  = 9,    # CENT_CAN
    amt  = 150,  # mg
    evid = 1,
    rate = -2
  )
)

## Scenario 7: Colchicine + Prednisolone combination
ev_combo1 <- rbind(
  as.data.frame(ev_colch),
  as.data.frame(ev_pred)
) %>% arrange(time)
ev_combo1_obj <- as.ev(ev_combo1)

## Scenario 8: Adalimumab + Apremilast combination (for refractory BD)
ev_combo2 <- rbind(
  as.data.frame(ev_ada),
  as.data.frame(ev_apremilast)
) %>% arrange(time)
ev_combo2_obj <- as.ev(ev_combo2)

## ============================================================
## Run Simulations
## ============================================================
run_scenario <- function(events, label, n_subj = 1) {
  bd_mod_refactored %>%
    init(bd_init_refactored) %>%
    mrgsim(events = events, end = 2160, delta = 4, carry_out = "time") %>%
    as.data.frame() %>%
    mutate(
      scenario   = label,
      time_days  = time / 24
    )
}

# Run all scenarios
cat("Running Behcet's Disease QSP simulations (refactored)...\n")
results <- list(
  run_scenario(ev_untreated,   "1_Untreated"),
  run_scenario(ev_colch,       "2_Colchicine"),
  run_scenario(ev_pred,        "3_Prednisolone"),
  run_scenario(ev_ada,         "4_Adalimumab_anti-TNF"),
  run_scenario(ev_apremilast,  "5_Apremilast_PDE4i"),
  run_scenario(ev_can,         "6_Canakinumab_anti-IL1B"),
  run_scenario(ev_combo1_obj,  "7_Colch+Pred_Combo"),
  run_scenario(ev_combo2_obj,  "8_Ada+Aprem_Combo_Refractory")
)
sim_all <- bind_rows(results)
cat("Simulations complete. Rows:", nrow(sim_all), "\n")

## ============================================================
## Visualization
## ============================================================
# Color palette
scen_colors <- c(
  "1_Untreated"                = "#e63946",
  "2_Colchicine"               = "#f4a261",
  "3_Prednisolone"             = "#2a9d8f",
  "4_Adalimumab_anti-TNF"      = "#457b9d",
  "5_Apremilast_PDE4i"         = "#52b788",
  "6_Canakinumab_anti-IL1B"    = "#9d4edd",
  "7_Colch+Pred_Combo"         = "#f48c06",
  "8_Ada+Aprem_Combo_Refractory" = "#1d3557"
)

## Plot 1: Disease Activity (BDCAF) — All Scenarios
p1 <- ggplot(sim_all, aes(x = time_days, y = BDCAF, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "BDCAF Disease Activity Score",
       x = "Time (days)", y = "BDCAF Score") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray50", alpha = 0.7) +
  annotate("text", x = 85, y = 3.3, label = "Remission threshold", size = 2.5, color = "gray50") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8))

## Plot 2: TNF-alpha dynamics
p2 <- ggplot(sim_all, aes(x = time_days, y = TNFA, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "TNF-α (Normalized)", x = "Time (days)", y = "TNF-α") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 3: IL-1beta dynamics
p3 <- ggplot(sim_all, aes(x = time_days, y = IL1B, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "IL-1β (Normalized)", x = "Time (days)", y = "IL-1β") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 4: IL-17A dynamics
p4 <- ggplot(sim_all, aes(x = time_days, y = IL17A, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "IL-17A (Normalized)", x = "Time (days)", y = "IL-17A") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 5: Oral Ulcer Activity
p5 <- ggplot(sim_all, aes(x = time_days, y = OUL, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "Oral Ulcer Activity Index", x = "Time (days)", y = "Oral Ulcer Score") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 6: Ocular Inflammation
p6 <- ggplot(sim_all, aes(x = time_days, y = OCI, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "Ocular Inflammation Index", x = "Time (days)", y = "Ocular Score") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 7: Endothelial Activation
p7 <- ggplot(sim_all, aes(x = time_days, y = EA, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "Endothelial Activation", x = "Time (days)", y = "EA (normalized)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## Plot 8: Th17 Activity
p8 <- ggplot(sim_all, aes(x = time_days, y = TH17, color = scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = scen_colors, name = "Treatment") +
  labs(title = "Th17 Cell Activity", x = "Time (days)", y = "Th17 (normalized)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

## PK Plot: Adalimumab concentration-time profile
pk_ada <- sim_all %>% filter(scenario == "4_Adalimumab_anti-TNF")
p_pk_ada <- ggplot(pk_ada, aes(x = time_days, y = Cp_ADA)) +
  geom_line(color = "#457b9d", linewidth = 1) +
  geom_hline(yintercept = 1500, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 70, y = 1700, label = "EC50 for TNF neutralization", size = 2.5, color = "red") +
  labs(title = "Adalimumab PK (SC 40 mg Q2W)", x = "Time (days)", y = "Adalimumab (ng/mL)") +
  theme_minimal(base_size = 11)

## Combined dashboard
dashboard <- (p1 | p_pk_ada) / (p2 | p3 | p4) / (p5 | p6 | p7 | p8)
print(dashboard)

## ============================================================
## Summary Table: Endpoint Reduction at Day 90
## ============================================================
summary_tbl <- sim_all %>%
  filter(time_days >= 88 & time_days <= 90) %>%
  group_by(scenario) %>%
  summarise(
    BDCAF_d90    = mean(BDCAF, na.rm = TRUE),
    OralUlcer_d90 = mean(OUL,   na.rm = TRUE),
    Ocular_d90   = mean(OCI,   na.rm = TRUE),
    TNFA_d90     = mean(TNFA,  na.rm = TRUE),
    IL1B_d90     = mean(IL1B,  na.rm = TRUE),
    IL17A_d90    = mean(IL17A, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BDCAF_reduce_pct  = (8.0 - BDCAF_d90) / 8.0 * 100,
    OUL_reduce_pct    = (3.0 - OralUlcer_d90) / 3.0 * 100,
    Remission         = ifelse(BDCAF_d90 < 3, "Yes", "No")
  )

cat("\n=== Behcet's Disease QSP: Day-90 Outcomes (refactored) ===\n")
print(summary_tbl, n = 20, width = 120)

## ============================================================
## Virtual Patient Analysis: HLA-B51 Status
## ============================================================
cat("\n--- Simulating HLA-B51+ vs HLA-B51- patients under adalimumab ---\n")

# HLA-B51 negative (HLAB51_factor = 1.0)
bd_neg <- bd_mod_refactored %>% param(HLAB51_factor = 1.0) %>%
  init(bd_init_refactored) %>%
  mrgsim(events = ev_ada, end = 2160, delta = 4) %>%
  as.data.frame() %>%
  mutate(HLAB51 = "HLA-B51 Negative", time_days = time / 24)

# HLA-B51 positive (HLAB51_factor = 1.4)
bd_pos <- bd_mod_refactored %>% param(HLAB51_factor = 1.4) %>%
  init(bd_init_refactored) %>%
  mrgsim(events = ev_ada, end = 2160, delta = 4) %>%
  as.data.frame() %>%
  mutate(HLAB51 = "HLA-B51 Positive", time_days = time / 24)

hlab51_sim <- bind_rows(bd_neg, bd_pos)

p_hlab51 <- ggplot(hlab51_sim, aes(x = time_days, y = BDCAF, color = HLAB51)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("HLA-B51 Positive" = "#e63946", "HLA-B51 Negative" = "#457b9d")) +
  labs(title = "HLA-B51 Genotype Impact on Adalimumab Response",
       subtitle = "Behcet's Disease QSP Model (refactored)",
       x = "Time (days)", y = "BDCAF Score", color = "Genotype") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(p_hlab51)

## ============================================================
## Dose-Response Analysis: Adalimumab
## ============================================================
doses_ada <- c(10, 20, 40, 80)  # mg SC Q2W

dose_response <- purrr::map_dfr(doses_ada, function(d) {
  ev_d <- ev(data.frame(time = seq(0, 2160, by = 336), cmt = 6, amt = d, evid = 1, rate = -2))
  bd_mod_refactored %>%
    init(bd_init_refactored) %>%
    mrgsim(events = ev_d, end = 2160, delta = 24) %>%
    as.data.frame() %>%
    filter(time == 2160) %>%
    mutate(dose_mg = d, time_days = time / 24)
})

cat("\n=== Adalimumab Dose-Response (Day 90) ===\n")
print(dose_response %>% select(dose_mg, BDCAF, OUL, OCI, TNFA))

cat("\nBehcet's Disease QSP Model (refactored) — Analysis Complete\n")
cat("Key insights:\n")
cat("  1. Anti-TNF (adalimumab) shows strong effect on ocular BD\n")
cat("  2. Colchicine primarily controls mucocutaneous flares\n")
cat("  3. Canakinumab (anti-IL-1B) effective for refractory oral ulcers\n")
cat("  4. Apremilast achieves good oral ulcer control without immunosuppression\n")
cat("  5. HLA-B51+ patients require more aggressive therapy\n")
cat("  6. Combination therapy (adalimumab + apremilast) superior for multi-organ BD\n")
