# ============================================================================
# VEXAS Syndrome QSP Model — mrgsolve specification (REFACTORED for pluggable PK)
#
# Disease: VEXAS (Vacuoles, E1 enzyme, X-linked, Autoinflammatory, Somatic)
# Mechanism: somatic UBA1 (codon Met41) mutation in HSC -> loss of cytosolic
#            UBA1b -> defective cytoplasmic ubiquitination -> ER stress, UPR,
#            mitochondrial dysfunction, NF-kB / NLRP3 / type-I IFN
#            hyperactivation -> autoinflammation + myelodysplastic cytopenias.
#
# This file is a sibling of the original `vexas_mrgsolve_model.R` (untouched,
# never edited) produced by the fork's PK/PD refactor pass — see
# FORK_WORKFLOW_GUIDE.md Part 2 and `vexas_refactor_notes.md` for the full
# rationale. In one sentence: every compound's PK/effect block is renamed to
# the fork's naming convention (GUT_<STEM>/CENT_<STEM>/CL_<STEM>/V1_<STEM>/
# KA_<STEM>/F_<STEM>, C_<STEM>, EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>,
# EFFECT_<STEM>) with NO change to PK archetype, NO refit, and NO change to
# any disease-biology equation. Compounds covered (per
# driver-patches/data/compound_perturbation_census.md): ANA (anakinra), AZA
# (azacitidine), CAN (canakinumab), PRED (prednisone), RUX (ruxolitinib),
# TOC (tocilizumab).
#
# Verified byte-for-byte against the original via the qspserver mrgsolve_api
# container (POST /model_manifest, POST /run_simulation) across the
# original's own five dosing scenarios (2-6) plus one bespoke canakinumab
# dosing test (the original file's own scenarios() never doses CAN — see
# notes): max relative deviation across every $CAPTURE output and every
# $CMT state, full time grid, all six scenarios = 0.0 (numerically exact;
# the one non-zero figure recorded, ~6.5e-13, is floating-point noise on a
# concentration that is analytically zero in that scenario). This is a pure
# structural reorganization / rename, not a refit — see
# vexas_refactor_notes.md for the full verification table.
#
# 23 ODE compartments cover:
#   - HSC clone dynamics + VAF (UBA1 mutant fraction)
#   - Misfolded protein burden / ER stress / mitochondrial ROS
#   - Core cytokines: IL-6, IL-1b, TNF-a, IFN-a, CXCL8, CCL2
#   - Acute-phase: CRP, ferritin
#   - Hematology: Hb (macrocytic anemia), Platelets, Neutrophils
#   - Clinical: fever index, skin/chondritis activity, VTE risk
#   - Steroid HPA suppression
#   - Drug PK: prednisone, tocilizumab, anakinra, canakinumab, ruxolitinib
#             (all archetype-3-minus-peripheral: depot + central, linear
#             elimination), azacitidine (same PK archetype, but see notes
#             on its disconnected PD)
#
# Therapy scenarios in `scenarios()` (unchanged dosing from the original,
# only compartment names updated):
#   1. Untreated natural history
#   2. Prednisone 1 mg/kg/d -> taper
#   3. Tocilizumab 162 mg SC q1w + low-dose prednisone
#   4. Anakinra 100 mg SC q24h
#   5. Ruxolitinib 10 mg PO BID + 10 mg/d prednisone
#   6. Azacitidine 75 mg/m^2 SC d1-7 q28 + supportive
#   7. Allogeneic HSCT (instantaneous clone replacement)
#
# Original author: QSP-routine 2026-06-30  |  Refactor: PK/PD fork pass
# Compatible with mrgsolve >= 1.5 (also verified against mrgsolve 2.0.1 via
# the qspserver mrgsolve_api container — see vexas_refactor_notes.md)
# ============================================================================

library(mrgsolve)

vexas_code_refactored <- '
$PROB
# VEXAS syndrome QSP model (23 ODEs) -- REFACTORED for pluggable PK
# Original author: QSP-routine 2026-06-30
# Refactor: renamed each compound PK/effect block to the fork naming
# convention (GUT_<STEM>/CENT_<STEM>/CL_<STEM>/V1_<STEM>/KA_<STEM>/F_<STEM>,
# C_<STEM>, EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>, EFFECT_<STEM>); no PK
# archetype change, no refit. See vexas_refactor_notes.md.
# Time units: hours.  Concentration units indicated per compartment.

$PARAM @annotated
// ---- Disease genetics / clone ----
VAF0       :  0.55 : Baseline UBA1 mutant variant allele fraction (monocyte)
k_clone    : 0.0008: Clonal expansion rate (1/h) [age-driven]
k_HSCT     : 0.0   : HSCT-induced clone eradication rate (1/h)
k_Aza      : 0.0   : Azacitidine clone-modulating rate (1/h) -- NOT driven by C_AZA anywhere in this file, unchanged from original; see refactor notes

// ---- Proteostasis ----
kp_misf    : 0.05  : Misfolded protein production proportional to VAF (au/h)
kd_misf    : 0.02  : Baseline misfolded clearance (1/h)
hill_misf  : 2     : Hill coefficient (UPR/inflammation)

// ---- Mitochondrial / ER stress ----
k_ROS_in   : 0.03  : ROS generation by misfolded burden (au/h)
k_ROS_out  : 0.30  : ROS scavenging (1/h)
k_ERst_in  : 0.04  : ER-stress accumulation (au/h)
k_ERst_out : 0.10  : ER-stress resolution (1/h)

// ---- Cytokine synthesis (NF-kB / inflammasome / IFN) ----
k_IL1b_syn : 0.25  : IL-1b synthesis rate (pg/mL/h)
k_IL1b_deg : 0.40  : IL-1b degradation (1/h, t1/2 ~= 1.7h)
k_IL6_syn  : 0.50  : IL-6 synthesis rate (pg/mL/h)
k_IL6_deg  : 0.40  : IL-6 degradation (1/h)
k_TNF_syn  : 0.30  : TNF-a synthesis (pg/mL/h)
k_TNF_deg  : 0.50  : TNF-a degradation (1/h, t1/2 ~= 80 min)
k_IFN_syn  : 0.15  : IFN-a synthesis (au/h)
k_IFN_deg  : 0.35  : IFN-a degradation (1/h)
k_CXCL8_syn: 0.20  : CXCL8/IL-8 synthesis (pg/mL/h)
k_CXCL8_deg: 0.45  : CXCL8 degradation (1/h)
k_CCL2_syn : 0.10  : CCL2/MCP-1 synthesis (pg/mL/h)
k_CCL2_deg : 0.30  : CCL2 degradation (1/h)

// ---- Acute-phase reactants ----
k_CRP_syn  : 0.08  : CRP synthesis driven by IL-6 (mg/L/h)
EC50_IL6_CRP: 35   : IL-6 EC50 for CRP production (pg/mL)
k_CRP_deg  : 0.038 : CRP first-order elimination (1/h, t1/2 ~= 18h)
k_FER_syn  : 0.05  : Ferritin synthesis from IL-6/IL-1b (ng/mL/h)
k_FER_deg  : 0.012 : Ferritin elimination (1/h, t1/2 ~= 60h)
EC50_IL6_FER: 50   : IL-6 EC50 for ferritin (pg/mL)

// ---- Hematology ----
Hb0        : 13.0  : Baseline normal Hb (g/dL)
k_Hb_in    : 0.018 : Hb production rate (g/dL / h)
k_Hb_loss  : 0.0012: Hb intrinsic loss/turnover (1/h ~ 25-d t1/2)
Imax_Hb    : 0.65  : Max % suppression of erythropoiesis by inflammation
EC50_Hb    : 80    : Combined inflam index for half-max Hb suppression
PLT0       : 250   : Baseline platelets (x10^9/L)
k_PLT_in   : 0.30  : Platelet production rate (x10^9/L/h)
k_PLT_loss : 0.0058: Platelet elim (1/h, t1/2 ~= 5d)
Imax_PLT   : 0.55  : Max % suppression of platelet production
EC50_PLT   : 70    : Inflam EC50 for PLT suppression
ANC0       : 4.5   : Baseline ANC (x10^9/L)
k_ANC_in   : 0.025 : ANC production
k_ANC_loss : 0.10  : ANC turnover (1/h)
Imax_ANC   : 0.40  : Max % suppression
EC50_ANC   : 70    : Inflam EC50 for ANC

// ---- Clinical activity ----
k_Fev_syn  : 0.05  : Fever index buildup rate
k_Fev_dec  : 0.10  : Fever resolution (1/h)
k_Skin_syn : 0.03  : Skin/chondritis activity buildup
k_Skin_dec : 0.05  : Skin activity resolution
k_VTE_syn  : 0.0008: VTE risk buildup (1/h) driven by IL-6+ferritin
k_VTE_dec  : 0.005 : VTE risk decay (1/h)
EC50_VTE   : 60    : Inflam EC50 for VTE risk

// ---- HPA axis (steroid suppression) ----
HPA0       : 1.0   : Baseline endogenous cortisol output (au)
k_HPA_supp : 0.04  : Suppression of HPA by exogenous GC (1/h per mg/mL drug)
k_HPA_rec  : 0.005 : HPA recovery rate (1/h)

// ---- Drug PK: PREDNISONE (oral, depot+central, linear -- archetype 3 minus peripheral) ----
KA_PRED    : 0.9   : Absorption rate constant (1/h) (was ka_PRED)
F_PRED     : 0.85  : Bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_PRED    : 16    : Clearance (L/h, ~3h t1/2)
V1_PRED    : 60    : Central volume (L) (was V_PRED)
EMAX_PRED  : 0.6   : Max fractional suppression of NF-kB activity by prednisone (was a hardcoded 0.6 coefficient in NFkB_act)
EC50_PRED  : 0.05  : Prednisone EC50 for NF-kB suppression (mg/L) (was a hardcoded 0.05 in NFkB_act)
GAMMA_PRED : 1     : Hill coefficient (original had no explicit Hill exponent => 1)

// ---- Drug PK: TOCILIZUMAB (SC, depot+central, linear -- archetype 3 minus peripheral) ----
KA_TOC     : 0.012 : SC absorption rate constant (1/h ~4d t1/2 abs) (was ka_TOC)
F_TOC      : 0.80  : SC bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_TOC     : 0.013 : Linear clearance (L/h)
V1_TOC     : 5.0   : Central volume (L) (was V_TOC)
IL6RMAX_TOC: 5.0   : Apparent IL-6R density (mg/L) (was IL6Rmax; declared in original, never referenced anywhere in $ODE -- see refactor notes; preserved unused)
EMAX_TOC   : 0.95  : Max fractional suppression of IL-6 signaling by TOC (was a hardcoded 0.95 coefficient on occ_TOC)
EC50_TOC   : 0.5   : TOC EC50 = apparent TOC-IL6R Kd (mg/L) (was Kd_TOC)
GAMMA_TOC  : 1     : Hill coefficient (occ_TOC = C/(Kd+C) already exactly Hill-shaped, gamma=1)

// ---- Drug PK: ANAKINRA (SC, depot+central, linear -- archetype 3 minus peripheral) ----
KA_ANA     : 0.60  : Absorption rate constant (1/h) (was ka_ANA)
F_ANA      : 0.95  : SC bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_ANA     : 2.8   : Clearance (L/h, t1/2 ~=3-4h)
V1_ANA     : 12    : Central volume (L) (was V_ANA)
EMAX_ANA   : 0.90  : Max fractional suppression of IL-1b by ANA (was a hardcoded 0.90 coefficient on occ_ANA)
EC50_ANA   : 1.0   : ANA EC50 = apparent ANA-IL1R Kd (mg/L) (was Kd_ANA)
GAMMA_ANA  : 1     : Hill coefficient (occ_ANA = C/(Kd+C) already exactly Hill-shaped, gamma=1)

// ---- Drug PK: CANAKINUMAB (SC mAb, depot+central, linear -- archetype 3 minus peripheral) ----
KA_CAN     : 0.005 : Absorption rate constant (1/h) (was ka_CAN)
F_CAN      : 0.66  : Bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_CAN     : 0.005 : Clearance (L/h, t1/2 ~=26d)
V1_CAN     : 6.0   : Central volume (L) (was V_CAN)
EMAX_CAN   : 0.85  : Max fractional suppression of IL-1b by CAN (was a hardcoded 0.85 coefficient on occ_CAN)
EC50_CAN   : 0.3   : CAN EC50 = apparent Kd (mg/L) (was Kd_CAN)
GAMMA_CAN  : 1     : Hill coefficient (occ_CAN = C/(Kd+C) already exactly Hill-shaped, gamma=1)

// ---- Drug PK: RUXOLITINIB (oral, depot+central, linear -- archetype 3 minus peripheral) ----
KA_RUX     : 1.6   : Absorption rate constant (1/h) (was ka_RUX)
F_RUX      : 0.95  : Bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_RUX     : 22    : Clearance (L/h, t1/2 ~=3h)
V1_RUX     : 75    : Central volume (L) (was V_RUX)
EMAX_RUX   : 0.70  : Max fractional suppression of IFN/JAK-driven signaling by RUX (was a hardcoded 0.70 coefficient on occ_RUX)
EC50_RUX   : 100   : RUX EC50 = apparent JAK Kd (ug/L) (was Kd_RUX)
GAMMA_RUX  : 1     : Hill coefficient (occ_RUX = C/(Kd+C) already exactly Hill-shaped, gamma=1)

// ---- Drug PK: AZACITIDINE (SC, depot+central, linear -- archetype 3 minus peripheral) ----
KA_AZA     : 1.8   : Absorption rate constant (1/h) (was ka_AZA)
F_AZA      : 0.89  : SC bioavailability (declared in original, never referenced in any dxdt -- see refactor notes; preserved unused)
CL_AZA     : 250   : Clearance (L/h, t1/2 ~=0.7h)
V1_AZA     : 76    : Central volume (L) (was V_AZA)
// No EMAX_AZA/EC50_AZA/GAMMA_AZA/EFFECT_AZA: unchanged from original, C_AZA
// is exposed (below) but is not read by any PD equation anywhere in this
// file -- dxdt_VAF clone-suppression term uses the free-standing k_Aza
// parameter above instead, which nothing here ever sets from C_AZA. Not
// invented/fixed here; see refactor notes and UPSTREAM_ISSUES.md.

// ---- Body / dosing helpers ----
WT         : 75    : Body weight (kg)
BSA        : 1.85  : Body surface area (m^2)

$CMT @annotated
// Disease/biology compartments (unchanged from original)
VAF        : UBA1 mutant variant allele fraction (unitless 0-1)
MISF       : Misfolded protein burden (au)
ROS        : ROS / mitochondrial stress (au)
ERST       : ER stress / UPR (au)
IL1B       : IL-1b (pg/mL)
IL6        : IL-6 (pg/mL)
TNFa       : TNF-a (pg/mL)
IFNa       : IFN-a (au)
CXCL8      : CXCL8 (pg/mL)
CCL2       : CCL2 (pg/mL)
CRP        : CRP (mg/L)
FER        : Ferritin (ng/mL)
HB         : Hb (g/dL)
PLT        : Platelet count (x10^9/L)
ANC        : Absolute neutrophil count (x10^9/L)
FEV        : Fever activity index (0-100)
SKIN       : Skin/chondritis activity index (0-100)
VTE        : VTE risk index (0-100)
HPA        : HPA axis output (fraction of baseline)
// Drug PK depots & central compartments (renamed to convention; same order as original)
GUT_PRED   : Prednisone GI depot (mg) (was AGUT_PRED)
CENT_PRED  : Prednisone central compartment (mg) (was CC_PRED)
GUT_TOC    : Tocilizumab SC depot (mg) (was ADEP_TOC)
CENT_TOC   : Tocilizumab central compartment (mg) (was CC_TOC)
GUT_ANA    : Anakinra SC depot (mg) (was ADEP_ANA)
CENT_ANA   : Anakinra central compartment (mg) (was CC_ANA)
GUT_CAN    : Canakinumab SC depot (mg) (was ADEP_CAN)
CENT_CAN   : Canakinumab central compartment (mg) (was CC_CAN)
GUT_RUX    : Ruxolitinib GI depot (mg) (was AGUT_RUX)
CENT_RUX   : Ruxolitinib central compartment (mg) (was CC_RUX)
GUT_AZA    : Azacitidine SC depot (mg) (was ADEP_AZA)
CENT_AZA   : Azacitidine central compartment (mg) (was CC_AZA)

$MAIN
VAF_0   = VAF0;
MISF_0  = 30.0;
ROS_0   = 5.0;
ERST_0  = 6.0;
IL1B_0  = 8.0;
IL6_0   = 120.0;
TNFa_0  = 25.0;
IFNa_0  = 8.0;
CXCL8_0 = 60.0;
CCL2_0  = 280.0;
CRP_0   = 110.0;
FER_0   = 1800.0;
HB_0    = 9.5;
PLT_0   = 130.0;
ANC_0   = 5.5;
FEV_0   = 40.0;
SKIN_0  = 35.0;
VTE_0   = 25.0;
HPA_0   = 1.0;

$ODE
// ---- Drug concentrations: the single redirect point per compound ----
double C_PRED  = CENT_PRED / V1_PRED;        // mg/L prednisone
double C_TOC   = CENT_TOC  / V1_TOC;         // mg/L
double C_ANA   = CENT_ANA  / V1_ANA;         // mg/L
double C_CAN   = CENT_CAN  / V1_CAN;         // mg/L
double C_RUX   = CENT_RUX  / V1_RUX * 1000;  // ug/L (unit conversion preserved from original)
double C_AZA   = CENT_AZA  / V1_AZA;         // mg/L

// ---- Drug PD: named Hill-interface effect terms (renamed only -- each is
// already exactly the original algebraic ratio/coefficient, gamma fixed
// at 1 because the original had no explicit Hill exponent on any of these) ----
double EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));
double EFFECT_TOC  = EMAX_TOC  * pow(C_TOC,  GAMMA_TOC)  / (pow(EC50_TOC,  GAMMA_TOC)  + pow(C_TOC,  GAMMA_TOC));
double EFFECT_ANA  = EMAX_ANA  * pow(C_ANA,  GAMMA_ANA)  / (pow(EC50_ANA,  GAMMA_ANA)  + pow(C_ANA,  GAMMA_ANA));
double EFFECT_CAN  = EMAX_CAN  * pow(C_CAN,  GAMMA_CAN)  / (pow(EC50_CAN,  GAMMA_CAN)  + pow(C_CAN,  GAMMA_CAN));
double EFFECT_RUX  = EMAX_RUX  * pow(C_RUX,  GAMMA_RUX)  / (pow(EC50_RUX,  GAMMA_RUX)  + pow(C_RUX,  GAMMA_RUX));

// ---- Effective free cytokine signaling (identical shape to original;
// occ_<drug>*coefficient replaced 1:1 by EFFECT_<STEM>) ----
double IL6_eff   = IL6 * (1.0 - EFFECT_TOC);
double IL1_eff   = IL1B * (1.0 - EFFECT_ANA - EFFECT_CAN);
double IFN_eff   = IFNa * (1.0 - EFFECT_RUX);
// Original: "(1 + C_PRED/(0.05+C_PRED) > 0 ? (1.0 - 0.6*C_PRED/(0.05+C_PRED)) : 1.0)".
// Since C_PRED >= 0, the ratio is in [0,1) so the guard is always true; this
// simplifies exactly (no numeric change) to 1 - EFFECT_PRED.
double NFkB_act  = 1.0 - EFFECT_PRED;

// Combined inflammation index (for hematology suppression)
double INFLAM    = IL6_eff + 2.0*IL1_eff + 0.5*TNFa + 0.3*IFN_eff;

// ---- Disease biology ODEs (unchanged) ----
dxdt_VAF  = (k_clone*VAF*(1-VAF)) - k_HSCT*VAF - k_Aza*VAF;

dxdt_MISF = kp_misf*100*VAF - kd_misf*MISF;
dxdt_ROS  = k_ROS_in*MISF - k_ROS_out*ROS;
dxdt_ERST = k_ERst_in*MISF - k_ERst_out*ERST;

// Cytokines: NF-kB driven (steroids down), inflammasome (IL-1b), IFN
dxdt_IL1B  = k_IL1b_syn*(ROS/(ROS + 5))*100*NFkB_act - k_IL1b_deg*IL1B;
dxdt_IL6   = k_IL6_syn *(ERST/(ERST + 6))*100*NFkB_act - k_IL6_deg*IL6;
dxdt_TNFa  = k_TNF_syn *(MISF/(MISF + 40))*100*NFkB_act - k_TNF_deg*TNFa;
dxdt_IFNa  = k_IFN_syn *(ROS/(ROS + 4))*100*(1.0 - EFFECT_RUX) - k_IFN_deg*IFNa;
dxdt_CXCL8 = k_CXCL8_syn*(IL1_eff/(IL1_eff + 5))*100 - k_CXCL8_deg*CXCL8;
dxdt_CCL2  = k_CCL2_syn *(IL6_eff/(IL6_eff + 40))*100 - k_CCL2_deg*CCL2;

// Acute-phase
dxdt_CRP   = k_CRP_syn*(IL6_eff/(EC50_IL6_CRP + IL6_eff))*100 - k_CRP_deg*CRP;
dxdt_FER   = k_FER_syn*((IL6_eff+IL1_eff)/(EC50_IL6_FER + IL6_eff + IL1_eff))*100 - k_FER_deg*FER;

// Hematology -- suppression by inflammation
double sup_Hb  = Imax_Hb  * INFLAM / (EC50_Hb + INFLAM);
double sup_PLT = Imax_PLT * INFLAM / (EC50_PLT + INFLAM);
double sup_ANC = Imax_ANC * INFLAM / (EC50_ANC + INFLAM);

dxdt_HB  = k_Hb_in *(1 - sup_Hb )*Hb0  - k_Hb_loss *HB;
dxdt_PLT = k_PLT_in*(1 - sup_PLT)*PLT0 - k_PLT_loss*PLT;
dxdt_ANC = k_ANC_in*(1 - sup_ANC)*ANC0 - k_ANC_loss*ANC;

// Clinical activity scores
dxdt_FEV  = k_Fev_syn *(IL1_eff/(IL1_eff + 10) + IL6_eff/(IL6_eff + 50))*50 - k_Fev_dec *FEV;
dxdt_SKIN = k_Skin_syn*(CXCL8/(CXCL8 + 30) + IL1_eff/(IL1_eff + 10))*50 - k_Skin_dec*SKIN;
dxdt_VTE  = k_VTE_syn *((IL6_eff+FER/100)/(EC50_VTE + IL6_eff + FER/100))*100 - k_VTE_dec*VTE;

// HPA axis
dxdt_HPA  = -k_HPA_supp*C_PRED*HPA + k_HPA_rec*(1 - HPA);

// ---- Drug PK ODEs (renamed only; F_<STEM> not applied, matching the
// original -- see refactor notes) ----
dxdt_GUT_PRED  = -KA_PRED*GUT_PRED;
dxdt_CENT_PRED =  KA_PRED*GUT_PRED - CL_PRED/V1_PRED*CENT_PRED;

dxdt_GUT_TOC   = -KA_TOC*GUT_TOC;
dxdt_CENT_TOC  =  KA_TOC*GUT_TOC - CL_TOC/V1_TOC*CENT_TOC;

dxdt_GUT_ANA   = -KA_ANA*GUT_ANA;
dxdt_CENT_ANA  =  KA_ANA*GUT_ANA - CL_ANA/V1_ANA*CENT_ANA;

dxdt_GUT_CAN   = -KA_CAN*GUT_CAN;
dxdt_CENT_CAN  =  KA_CAN*GUT_CAN - CL_CAN/V1_CAN*CENT_CAN;

dxdt_GUT_RUX   = -KA_RUX*GUT_RUX;
dxdt_CENT_RUX  =  KA_RUX*GUT_RUX - CL_RUX/V1_RUX*CENT_RUX;

dxdt_GUT_AZA   = -KA_AZA*GUT_AZA;
dxdt_CENT_AZA  =  KA_AZA*GUT_AZA - CL_AZA/V1_AZA*CENT_AZA;

$CAPTURE @annotated
C_PRED     : Prednisone plasma conc (mg/L)
C_TOC      : Tocilizumab plasma conc (mg/L)
C_ANA      : Anakinra plasma conc (mg/L)
C_CAN      : Canakinumab plasma conc (mg/L)
C_RUX      : Ruxolitinib plasma conc (ug/L)
C_AZA      : Azacitidine plasma conc (mg/L)
EFFECT_PRED: Prednisone fractional suppression of NF-kB activity (was inline NFkB_act coefficient, not previously captured)
EFFECT_TOC : Tocilizumab fractional suppression of IL-6 signaling (was occ_TOC*0.95)
EFFECT_ANA : Anakinra fractional suppression of IL-1b signaling (was occ_ANA*0.90)
EFFECT_CAN : Canakinumab fractional suppression of IL-1b signaling (was occ_CAN*0.85, not previously captured)
EFFECT_RUX : Ruxolitinib fractional suppression of IFN/JAK signaling (was occ_RUX*0.70)
INFLAM     : Composite inflammation index
IL6_eff    : Free/biologically-active IL-6 (pg/mL)
IL1_eff    : Free IL-1b (pg/mL)
IFN_eff    : Free IFN-a (au)
'

# ---- Compile ----------------------------------------------------------------
mod_vexas_refactored <- mcode("vexas_qsp_refactored", vexas_code_refactored)

# ---- Helpers ---------------------------------------------------------------
prednisone_taper <- function(start_mg = 60, taper_weeks = 12, end_mg = 5) {
  weeks <- 0:taper_weeks
  doses <- pmax(seq(start_mg, end_mg, length.out = length(weeks)), end_mg)
  data.frame(time = weeks*168, cmt = "GUT_PRED", amt = doses, evid = 1, ii = 24, addl = 6)
}

scenarios <- function() {
  list(
    "1_untreated"   = NULL,
    "2_pred_taper"  = prednisone_taper(60, 12, 5),
    "3_toci_lowGC"  = rbind(
      data.frame(time = seq(0, 84*24, by = 7*24), cmt = "GUT_TOC", amt = 162, evid = 1, ii = 0, addl = 0),
      data.frame(time = seq(0, 84*24, by = 24),  cmt = "GUT_PRED", amt = 10, evid = 1, ii = 0, addl = 0)
    ),
    "4_anakinra"    = data.frame(time = seq(0, 84*24, by = 24), cmt = "GUT_ANA", amt = 100, evid = 1, ii = 0, addl = 0),
    "5_ruxo_pred"   = rbind(
      data.frame(time = seq(0, 84*24, by = 12), cmt = "GUT_RUX", amt = 10, evid = 1, ii = 0, addl = 0),
      data.frame(time = seq(0, 84*24, by = 24), cmt = "GUT_PRED", amt = 10, evid = 1, ii = 0, addl = 0)
    ),
    "6_azacitidine" = data.frame(time = c(outer(seq(0, 6)*24, seq(0, 168, by = 28*24), `+`)),
                                  cmt = "GUT_AZA", amt = 75*1.85, evid = 1, ii = 0, addl = 0),
    "7_HSCT"        = data.frame(time = 24, cmt = "VAF", amt = 0, evid = 2, # set k_HSCT high transiently
                                  ii = 0, addl = 0)  # use $TABLE override in user code
  )
}

# ---- Example run ------------------------------------------------------------
if (interactive()) {
  library(dplyr); library(ggplot2); library(tidyr)
  evts <- scenarios()[["3_toci_lowGC"]]
  out  <- mod_vexas_refactored %>%
    ev(evts) %>%
    mrgsim(end = 84*24, delta = 6) %>% as.data.frame()
  out %>% select(time, IL6, IL6_eff, CRP, FER, HB, PLT) %>%
    pivot_longer(-time) %>%
    ggplot(aes(time/24, value)) + geom_line() + facet_wrap(~name, scales = "free_y") +
    labs(x = "Days", title = "VEXAS (refactored) -- Tocilizumab + low-dose prednisone")
}

# ---- Parameter sources / calibration notes (unchanged from original) -------
# - UBA1 VAF ranges: Beck et al. NEJM 2020 (median ~75% monocyte);
#   Georgin-Lavialle et al. AJH 2024.
# - IL-6 elevation: median 80-250 pg/mL in active disease (Patel 2022).
# - CRP: typically 80-200 mg/L during flare.
# - Ferritin: 1000-10,000 ng/mL; HLH-like crises >10,000.
# - Hb baseline 8-10 g/dL macrocytic (MCV >100); PLT 80-150x10^9/L.
# - Tocilizumab PK: V~5L, CL~10-15 mL/h, t1/2 11-13 d (SC, EMA/FDA labels).
# - Anakinra PK: t1/2 4-6h, F~95% SC (Kineret label).
# - Canakinumab PK: t1/2 26 d, F~66% SC (Ilaris label).
# - Ruxolitinib PK: t1/2 3h, CL 22 L/h (Jakafi label).
# - Azacitidine PK: t1/2 ~0.7h (Vidaza label, oral CC-486 differs).
# - Prednisone PK: F~0.85, CL 6-8 L/h/70kg, t1/2 3h (active prednisolone).
# - 5-yr survival ~50% untreated, ~60-70% on biologic+HSCT (Hines 2023).
# - Allogeneic HSCT modeled as instantaneous VAF->0 (Hadjadj 2024 cohort).
#
# ---- Refactor-specific notes -------------------------------------------
# See vexas_refactor_notes.md for: archetype classification per compound,
# the qspserver mrgsolve_api verification table (all 6 scenarios, max
# relative deviation), and two upstream defects logged (not fixed) in
# translations/UPSTREAM_ISSUES.md entries #85-86 (dead F_<drug>
# bioavailability parameters across all six compounds; azacitidine PK
# fully disconnected from PD).
# ============================================================================
