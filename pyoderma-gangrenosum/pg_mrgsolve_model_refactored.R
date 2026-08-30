# =============================================================================
# Pyoderma Gangrenosum (PG) — mrgsolve QSP Model — REFACTORED
# =============================================================================
# PK/PD refactor per FORK_WORKFLOW_GUIDE.md Part 2 ("pluggable PK, named Hill
# interface"). Original: pg_mrgsolve_model.R (untouched, ground truth).
#
# Scope of this refactor: all six compounds this file models are refactored
# (Adalimumab, Anakinra, Cyclosporine, Infliximab, Prednisone, Ustekinumab) --
# there are no other compounds in this file to leave untouched. Every
# disease-side compartment/equation (cytokine turnover, neutrophil/Th17/Treg/
# macrophage dynamics, MMP-9/ROS/ulcer kinetics, CRP/calprotectin, clinical
# endpoints) is copied verbatim from the original; only each compound's own
# PK block and its named effect term were renamed to the fork's convention.
#
# Per-compound archetype (determined from this file's own equations, not
# assumed from any other refactor of the same drug elsewhere in the corpus):
#   - Adalimumab (ADA): Archetype 3 (depot + central + peripheral, linear).
#   - Infliximab (IFX): Archetype 2 (no depot [IV], 2-cmt, linear).
#   - Anakinra (ANA):   Archetype 3 minus peripheral (depot + central, linear).
#   - Cyclosporine (CSA): Archetype 3 minus peripheral (depot + central, linear).
#   - Ustekinumab (UST): Archetype 3 minus peripheral (depot + central,
#     linear). Checked specifically for TMDD (per task instruction, since UST
#     is an anti-IL-12/23 antibody where receptor-binding kinetics would be
#     plausible) -- there are NO receptor/complex compartments or on/off-rate
#     parameters anywhere in this file for UST; its own $ODE is a plain
#     one-depot linear-elimination compartment (`dxdt_UST_CENT = KA_UST *
#     F1_UST * UST_SC - (CL_UST/V_UST) * UST_CENT`) and its effect on disease
#     is a plain algebraic Emax ratio on CONC_UST. Confirmed NOT TMDD in this
#     file; kept as the plain-PK archetype it actually is.
#   - Prednisone (PRED): Archetype 3 minus peripheral, minus bioavailability
#     (the original declares no `F1_PRED`-equivalent parameter and applies
#     none to the absorbed dose -- so no `F_PRED` was invented; see guide's
#     "never invent or default a PK parameter" rule).
#
# All six compounds' effect terms were already exact Hill ratios
# (Emax*C/(EC50+C), implicit gamma=1) in the original -- every EFFECT_<STEM>
# below is a RENAME, not a refit. See pg_refactor_notes.md for the full
# naming map, the shared-Emax finding (logged as a cross-compound parameter
# leak, same class as UPSTREAM_ISSUES.md #84), the build-compatibility fix
# required to compile this file at all under mrgsolve 2.0.1 (12 `$PARAM`
# annotation lines missing a description field), and the qspserver
# mrgsolve_api verification results.
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

pg_qsp_code <- '
$PROB
# Pyoderma Gangrenosum QSP model (refactored PK/PD interface)
# v1.0 — Drug PK/PD, neutrophilic inflammation, ulcer dynamics, healing
# Time unit: days

$PARAM @annotated
// -------- Adalimumab (ADA) — Archetype 3 (depot+central+peripheral) --------
KA_ADA       :  0.30   : ADA absorption rate (1/d)
F_ADA        :  0.64   : ADA SC bioavailability                          : was F1_ADA
CL_ADA       :  0.31   : ADA clearance (L/d)
V1_ADA       :  4.7    : ADA central volume (L)
Q_ADA        :  0.45   : ADA inter-compartmental Q (L/d)
V2_ADA       :  2.6    : ADA peripheral volume (L)

// -------- Infliximab (IFX) — Archetype 2 (no depot [IV], 2-cmt) --------
CL_IFX       :  0.32   : IFX clearance (L/d)
V1_IFX       :  3.5    : IFX central volume (L)
Q_IFX        :  0.40   : IFX Q (L/d)
V2_IFX       :  2.0    : IFX peripheral volume (L)

// -------- Anakinra (ANA) — Archetype 3 minus peripheral --------
KA_ANA       :  6.0    : ANA absorption rate (1/d)
F_ANA        :  0.95   : ANA SC bioavailability                          : was F1_ANA
CL_ANA       :  17.0   : ANA clearance (L/d)
V1_ANA       :  20.0   : ANA volume of distribution (L)                  : was V_ANA

// -------- Cyclosporine (CSA) — Archetype 3 minus peripheral --------
KA_CSA       :  1.4    : CSA absorption rate (1/d)
F_CSA        :  0.30   : CSA bioavailability                             : was F1_CSA
CL_CSA       :  27.0   : CSA clearance (L/d)
V1_CSA       :  85.0   : CSA volume of distribution (L)                  : was V_CSA

// -------- Ustekinumab (UST) — Archetype 3 minus peripheral (not TMDD, see header) --------
KA_UST       :  0.20   : UST absorption rate (1/d)
F_UST        :  0.57   : UST bioavailability                             : was F1_UST
CL_UST       :  0.45   : UST clearance (L/d)
V1_UST       :  4.6    : UST volume of distribution (L)                  : was V_UST

// -------- Prednisone (PRED) — Archetype 3 minus peripheral, no F term (original has none) --------
KA_PRED      :  4.0    : PRED absorption rate (1/d)
CL_PRED      :  90.0   : PRED clearance (L/d)
V1_PRED      :  35.0   : PRED volume of distribution (L)                 : was V_PRED

// -------- PD: cytokine homeostasis (turnover) --------
kin_TNF      :  3.0    : TNF synthesis (pg/mL/d)
kout_TNF     :  3.0    : TNF degradation (1/d)
kin_IL1      :  4.0    : IL-1b synthesis baseline
kout_IL1     :  4.0    : IL-1b degradation (1/d)
kin_IL17     :  2.0    : IL-17A synthesis baseline
kout_IL17    :  2.0    : IL-17A degradation (1/d)
kin_IL6      :  3.0    : IL-6 baseline
kout_IL6     :  3.0    : IL-6 degradation
kin_IL23     :  1.5    : IL-23 baseline
kout_IL23    :  1.5    : IL-23 degradation
kin_IL8      :  6.0    : IL-8 baseline
kout_IL8     :  6.0    : IL-8 degradation

// -------- Disease forcing (PG-specific) --------
DSEV         :  3.0    : Disease severity multiplier (1=remission, 3-5=active)
F_pathergy   :  0.5    : Pathergy event toggle (0/1)
F_comorb_IBD :  0.6    : IBD comorbidity coefficient

// -------- Inflammatory drivers (cross-cytokine) --------
k_TNF_drive   : 0.20   : disease drives TNF
k_IL1_drive   : 0.25   : disease drives IL-1
k_IL17_drive  : 0.20   : Th17/IL-23 drives IL-17
k_IL6_drive   : 0.15   : TNF/IL-1 drives IL-6
k_IL23_drive  : 0.10   : DC drives IL-23
k_IL8_drive   : 0.40   : TNF+IL-1+IL-17 drives IL-8

// -------- Cell dynamics --------
kin_Neu      : 0.6     : neutrophil influx baseline (1e9/L/d)
kout_Neu     : 0.5     : neutrophil clearance (1/d)
k_NETform    : 0.5     : NET formation per neutrophil unit
kout_NET     : 0.8     : NET degradation/clearance
kin_Th17     : 0.05    : Th17 baseline (au/d)
kout_Th17    : 0.05    : Th17 clearance (1/d)
kin_Treg     : 0.06    : Treg baseline (au/d)
kout_Treg    : 0.06    : Treg clearance (1/d)
kin_M1       : 0.20    : M1 macrophage baseline (au/d)
kout_M1      : 0.20    : M1 macrophage clearance (1/d)

// -------- Tissue / ulcer kinetics --------
k_ulcer_grow : 0.30    : ulcer growth from inflammation (cm2/d)
k_heal_max   : 0.18    : maximum healing rate (cm2/d)
EC50_heal_drug : 0.20  : drug effect on healing rate
UlcerArea_0  : 12.0    : baseline ulcer area cm2
k_MMP        : 0.5     : MMP-9 synthesis rate
kout_MMP     : 0.4     : MMP-9 degradation
k_ROS        : 0.6     : ROS gen
kout_ROS     : 0.5     : ROS clearance (1/d)

// -------- Endpoints --------
k_CRP_syn    : 0.6     : CRP from IL-6 (mg/L per IL-6 unit)
kout_CRP     : 0.3     : CRP clearance (1/d)
k_Calp_syn   : 0.3     : calprotectin synthesis rate
kout_Calp    : 0.3     : calprotectin clearance (1/d)

// -------- Drug effect potencies: Hill interface (EC50/EMAX/GAMMA per compound) --------
// The original used ONE shared "Emax_drug = 0.92" for all six compounds
// (see pg_refactor_notes.md, "shared EMAX finding") -- split here into six
// independent EMAX_<STEM> parameters, each carrying the identical original
// value, so each compound stays independently driveable per the fork
// convention (no numeric change: overriding none of them reproduces the
// original`s single shared value exactly).
EC50_ADA     : 0.20    : ADA EC50 on TNF (ug/mL)                    : was IC50_ADA_TNF
EMAX_ADA     : 0.92    : ADA max effect                             : was shared Emax_drug
GAMMA_ADA    : 1       : ADA Hill coefficient (original had none, added as 1)
EC50_IFX     : 0.10    : IFX EC50 on TNF (ug/mL)                    : was IC50_IFX_TNF
EMAX_IFX     : 0.92    : IFX max effect                             : was shared Emax_drug
GAMMA_IFX    : 1       : IFX Hill coefficient (original had none, added as 1)
EC50_ANA     : 0.50    : ANA EC50 on IL-1 axis (ug/mL)              : was IC50_ANA_IL1
EMAX_ANA     : 0.92    : ANA max effect                             : was shared Emax_drug
GAMMA_ANA    : 1       : ANA Hill coefficient (original had none, added as 1)
EC50_UST     : 0.50    : UST EC50 on IL-23 (ug/mL)                  : was IC50_UST_IL23
EMAX_UST     : 0.92    : UST max effect                             : was shared Emax_drug
GAMMA_UST    : 1       : UST Hill coefficient (original had none, added as 1)
EC50_CSA     : 0.10    : CSA EC50 on Th17 (ng/mL)                   : was IC50_CSA_Th17
EMAX_CSA     : 0.92    : CSA max effect                             : was shared Emax_drug
GAMMA_CSA    : 1       : CSA Hill coefficient (original had none, added as 1)
EC50_PRED    : 5.0     : PRED EC50 on NF-kB-dependent cytokines (ng/mL) : was IC50_PRED
EMAX_PRED    : 0.92    : PRED max effect                            : was shared Emax_drug
GAMMA_PRED   : 1       : PRED Hill coefficient (original had none, added as 1)

$CMT @annotated
GUT_ADA   : adalimumab SC depot (mg)                 : was ADA_SC
CENT_ADA  : adalimumab central (mg)                  : was ADA_CENT
PERI_ADA  : adalimumab peripheral (mg)                : was ADA_PERI
CENT_IFX  : infliximab central (mg)                  : was IFX_CENT
PERI_IFX  : infliximab peripheral (mg)                : was IFX_PERI
GUT_ANA   : anakinra SC depot (mg)                   : was ANA_SC
CENT_ANA  : anakinra central (mg)                    : was ANA_CENT
GUT_CSA   : cyclosporine gut (mg)                    : was CSA_GUT
CENT_CSA  : cyclosporine central (mg)                : was CSA_CENT
GUT_UST   : ustekinumab SC depot (mg)                : was UST_SC
CENT_UST  : ustekinumab central (mg)                 : was UST_CENT
GUT_PRED  : prednisone gut (mg)                      : was PRED_GUT
CENT_PRED : prednisone central (mg)                  : was PRED_CENT
TNFa      : TNF-alpha (pg/mL)
IL1b      : IL-1beta (pg/mL)
IL17A     : IL-17A (pg/mL)
IL6       : IL-6 (pg/mL)
IL23      : IL-23 (pg/mL)
IL8       : IL-8 (pg/mL)
Neutroph  : tissue neutrophil index (au)
NET       : NET burden (au)
Th17      : Th17 index (au)
Treg      : Treg index (au)
M1        : macrophage M1 index (au)
MMP9_act  : active MMP-9 (au)
ROS_les   : lesional ROS (au)
Ulcer     : ulcer area (cm2)
Healed    : cumulative healed fraction (0-1)
CRP       : C-reactive protein (mg/L)
Calprot   : serum calprotectin (ng/mL)

$MAIN
CENT_ADA_0 = 0;
TNFa_0     = 3.0;
IL1b_0     = 4.0;
IL17A_0    = 2.0;
IL6_0      = 3.0;
IL23_0     = 1.5;
IL8_0      = 6.0;
Neutroph_0 = 1.0;
NET_0      = 0.5;
Th17_0     = 1.0;
Treg_0     = 1.0;
M1_0       = 1.0;
MMP9_act_0 = 1.0;
ROS_les_0  = 1.0;
Ulcer_0    = UlcerArea_0;
Healed_0   = 0.0;
CRP_0      = 6.0;
Calprot_0  = 100.0;

$ODE
// ===== Drug PK =====
double C_ADA  = CENT_ADA / V1_ADA;     // mg/L = ug/mL
double C_IFX  = CENT_IFX / V1_IFX;
double C_ANA  = CENT_ANA / V1_ANA;
double C_CSA  = CENT_CSA / V1_CSA * 1000.0; // ng/mL (mg/L*1000)
double C_UST  = CENT_UST / V1_UST;
double C_PRED = CENT_PRED / V1_PRED * 1000.0; // ng/mL

dxdt_GUT_ADA   = -KA_ADA  * GUT_ADA;
dxdt_CENT_ADA  =  KA_ADA  * F_ADA * GUT_ADA
                 - (CL_ADA/V1_ADA) * CENT_ADA
                 - (Q_ADA/V1_ADA) * CENT_ADA + (Q_ADA/V2_ADA) * PERI_ADA;
dxdt_PERI_ADA  =  (Q_ADA/V1_ADA) * CENT_ADA - (Q_ADA/V2_ADA) * PERI_ADA;

dxdt_CENT_IFX  = -(CL_IFX/V1_IFX) * CENT_IFX
                 - (Q_IFX/V1_IFX) * CENT_IFX + (Q_IFX/V2_IFX) * PERI_IFX;
dxdt_PERI_IFX  =  (Q_IFX/V1_IFX) * CENT_IFX - (Q_IFX/V2_IFX) * PERI_IFX;

dxdt_GUT_ANA   = -KA_ANA  * GUT_ANA;
dxdt_CENT_ANA  =  KA_ANA  * F_ANA * GUT_ANA - (CL_ANA/V1_ANA) * CENT_ANA;

dxdt_GUT_CSA   = -KA_CSA  * GUT_CSA;
dxdt_CENT_CSA  =  KA_CSA  * F_CSA * GUT_CSA - (CL_CSA/V1_CSA) * CENT_CSA;

dxdt_GUT_UST   = -KA_UST  * GUT_UST;
dxdt_CENT_UST  =  KA_UST  * F_UST * GUT_UST - (CL_UST/V1_UST) * CENT_UST;

dxdt_GUT_PRED  = -KA_PRED * GUT_PRED;
dxdt_CENT_PRED =  KA_PRED * GUT_PRED - (CL_PRED/V1_PRED) * CENT_PRED;

// ===== PD: drug effects (Hill interface; renames of the original`s already-
// Hill-shaped Emax ratios -- gamma=1 in every case, added explicitly since
// the original had no Hill exponent term) =====
double EFFECT_ADA  = EMAX_ADA  * pow(C_ADA,  GAMMA_ADA)  / (pow(EC50_ADA,  GAMMA_ADA)  + pow(C_ADA,  GAMMA_ADA));
double EFFECT_IFX  = EMAX_IFX  * pow(C_IFX,  GAMMA_IFX)  / (pow(EC50_IFX,  GAMMA_IFX)  + pow(C_IFX,  GAMMA_IFX));
double EFFECT_ANA  = EMAX_ANA  * pow(C_ANA,  GAMMA_ANA)  / (pow(EC50_ANA,  GAMMA_ANA)  + pow(C_ANA,  GAMMA_ANA));
double EFFECT_UST  = EMAX_UST  * pow(C_UST,  GAMMA_UST)  / (pow(EC50_UST,  GAMMA_UST)  + pow(C_UST,  GAMMA_UST));
double EFFECT_CSA  = EMAX_CSA  * pow(C_CSA,  GAMMA_CSA)  / (pow(EC50_CSA,  GAMMA_CSA)  + pow(C_CSA,  GAMMA_CSA));
double EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));

// Combined TNF suppression -- dead/unused in the original (dxdt_TNFa below
// recomputes the same product inline instead of reading this variable);
// preserved verbatim, harmless.
double SUP_TNF  = 1.0 - (1.0 - (1-EFFECT_ADA)*(1-EFFECT_IFX)*(1-EFFECT_PRED)) ;
double SUP_IL1  = (1 - EFFECT_ANA) * (1 - EFFECT_PRED);
double SUP_IL17 = (1 - EFFECT_UST) * (1 - EFFECT_CSA);
double SUP_IL23 = (1 - EFFECT_UST);
double SUP_IL6  = (1 - EFFECT_PRED) * (1 - EFFECT_IFX * 0.6);

// Disease-driven forcing
double F_dis = DSEV * (1.0 + 0.5*F_pathergy + 0.4*F_comorb_IBD);

// ===== Cytokine ODEs (synthesis driven; drug-suppressed) =====
dxdt_TNFa  = (kin_TNF  + k_TNF_drive  * F_dis * (Neutroph + Th17 + M1))
               * (1 - EFFECT_ADA)*(1 - EFFECT_IFX)*(1 - EFFECT_PRED)
             - kout_TNF * TNFa;

dxdt_IL1b  = (kin_IL1  + k_IL1_drive  * F_dis * (Neutroph + M1)) * SUP_IL1
             - kout_IL1 * IL1b;

dxdt_IL17A = (kin_IL17 + k_IL17_drive * F_dis * Th17) * SUP_IL17
             - kout_IL17 * IL17A;

dxdt_IL6   = (kin_IL6  + k_IL6_drive  * F_dis * (TNFa + IL1b)/10.0) * SUP_IL6
             - kout_IL6 * IL6;

dxdt_IL23  = (kin_IL23 + k_IL23_drive * F_dis * M1) * SUP_IL23
             - kout_IL23 * IL23;

dxdt_IL8   = (kin_IL8  + k_IL8_drive  * F_dis * (TNFa + IL1b + IL17A)/15.0)
               * (1 - EFFECT_PRED*0.6)
             - kout_IL8 * IL8;

// ===== Cell ODEs =====
double rec_neu = IL8 / 6.0;            // recruitment scaled by IL-8
dxdt_Neutroph = kin_Neu + k_IL8_drive * 0.1 * rec_neu * F_dis - kout_Neu * Neutroph;
dxdt_NET      = k_NETform * Neutroph * (1 + ROS_les) - kout_NET * NET;
dxdt_Th17     = kin_Th17 + 0.1 * IL23/1.5 * F_dis * SUP_IL17 - kout_Th17 * Th17;
dxdt_Treg     = kin_Treg - kout_Treg * Treg - 0.02 * Th17;   // Th17 suppresses Treg
dxdt_M1       = kin_M1   + 0.05 * (TNFa + IL1b)/7.0 * F_dis - kout_M1 * M1;

// ===== Tissue: MMP-9, ROS, ulcer =====
dxdt_MMP9_act = k_MMP * (Neutroph + 0.3 * M1) - kout_MMP * MMP9_act;
dxdt_ROS_les  = k_ROS * Neutroph - kout_ROS * ROS_les;

// Ulcer growth from MMP-9 / ROS / NET ; healing accelerated by drugs (cytokine suppression)
double inflam_drive = (MMP9_act + 0.6*ROS_les + 0.5*NET);
double drug_heal_eff = EFFECT_ADA + EFFECT_IFX + EFFECT_ANA + EFFECT_UST + EFFECT_CSA + 0.5*EFFECT_PRED;
double heal_rate = k_heal_max * drug_heal_eff / (EC50_heal_drug + drug_heal_eff);

dxdt_Ulcer = k_ulcer_grow * inflam_drive * (Ulcer / (Ulcer + 5)) - heal_rate * Ulcer;
dxdt_Healed = heal_rate * Ulcer / UlcerArea_0;   // fraction healed

// ===== Endpoints =====
dxdt_CRP     = k_CRP_syn * IL6 - kout_CRP * CRP;
dxdt_Calprot = k_Calp_syn * (Neutroph + 0.5*NET) - kout_Calp * Calprot;

$TABLE
double Pain_VAS   = std::min(10.0, 1.5 + 0.4*Ulcer + 0.1*TNFa);
double PARACELSUS = std::min(60.0, 5 + 1.5*Ulcer + 0.5*Pain_VAS + 2.0*(Neutroph>2 ? 1 : 0)*5
                                  + 0.8*M1 + 0.4*Th17);
double DLQI       = std::min(30.0, 4 + 0.6*Ulcer + 0.3*Pain_VAS);
double HiSCRpseudo = (Ulcer < 0.5 * UlcerArea_0) ? 1.0 : 0.0;  // >=50% reduction
double CompleteHeal = (Ulcer < 0.05 * UlcerArea_0) ? 1.0 : 0.0;

$CAPTURE @annotated
C_ADA  : Adalimumab conc (ug/mL)                    : was CONC_ADA
C_IFX  : Infliximab conc (ug/mL)                     : was CONC_IFX
C_ANA  : Anakinra conc (ug/mL)                       : was CONC_ANA
C_CSA  : Cyclosporine conc (ng/mL)                   : was CONC_CSA
C_UST  : Ustekinumab conc (ug/mL)                    : was CONC_UST
C_PRED : Prednisone conc (ng/mL)                      : was CONC_PRED
EFFECT_ADA  : ADA Hill effect on TNF (0-1)
EFFECT_IFX  : IFX Hill effect on TNF (0-1)
EFFECT_ANA  : ANA Hill effect on IL-1 axis (0-1)
EFFECT_UST  : UST Hill effect on IL-23 (0-1)
EFFECT_CSA  : CSA Hill effect on Th17 (0-1)
EFFECT_PRED : PRED Hill effect on NF-kB-dependent cytokines (0-1)
Pain_VAS  : VAS pain 0-10
PARACELSUS : PARACELSUS score
DLQI      : DLQI score
HiSCRpseudo : >=50% ulcer reduction
CompleteHeal : ulcer cleared
'

# -----------------------------------------------------------------------------
# Compile model
# -----------------------------------------------------------------------------
pg_qsp_model <- mcode("pg_qsp_refactored", pg_qsp_code)

# -----------------------------------------------------------------------------
# Helper: simulate a scenario over 168 days (24 weeks) with a defined dose regimen
# (compartment names updated to the fork's naming convention; same doses/timing)
# -----------------------------------------------------------------------------
simulate_pg <- function(model, scenario, tmax = 168, ...) {
  ev <- switch(scenario,
    # 1) Standard of care: Prednisone 60 mg/d x14d -> taper x 12 weeks
    "Prednisone_SOC" = ev(amt = 60, cmt = "GUT_PRED", ii = 1, addl = 13) +
                       ev(amt = 40, cmt = "GUT_PRED", time = 14, ii = 1, addl = 13) +
                       ev(amt = 20, cmt = "GUT_PRED", time = 28, ii = 1, addl = 27) +
                       ev(amt = 10, cmt = "GUT_PRED", time = 56, ii = 1, addl = 55),
    # 2) Cyclosporine: 4 mg/kg/d (~ 280 mg) x6 mo
    "Cyclosporine"   = ev(amt = 140, cmt = "GUT_CSA", ii = 0.5, addl = 335),  # BID
    # 3) Infliximab induction 5 mg/kg @0,2,6 wk then q8w (assume 70 kg = 350 mg)
    "Infliximab"     = ev(amt = 350, cmt = "CENT_IFX", time = 0) +
                       ev(amt = 350, cmt = "CENT_IFX", time = 14) +
                       ev(amt = 350, cmt = "CENT_IFX", time = 42) +
                       ev(amt = 350, cmt = "CENT_IFX", time = 98) +
                       ev(amt = 350, cmt = "CENT_IFX", time = 154),
    # 4) Adalimumab: 80 mg loading -> 40 mg q1w
    "Adalimumab"     = ev(amt = 80, cmt = "GUT_ADA", time = 0) +
                       ev(amt = 40, cmt = "GUT_ADA", time = 7, ii = 7, addl = 22),
    # 5) Anakinra 100 mg/d SC
    "Anakinra"       = ev(amt = 100, cmt = "GUT_ANA", ii = 1, addl = 167),
    # 6) Ustekinumab 90 mg SC q12w (after 0 wk + 4 wk loading)
    "Ustekinumab"    = ev(amt = 90, cmt = "GUT_UST", time = 0) +
                       ev(amt = 90, cmt = "GUT_UST", time = 28) +
                       ev(amt = 90, cmt = "GUT_UST", time = 112),
    # 7) Combo: Prednisone short course + Cyclosporine maintenance
    "Combo_PRED_CSA" = ev(amt = 60, cmt = "GUT_PRED", ii = 1, addl = 13) +
                       ev(amt = 30, cmt = "GUT_PRED", time = 14, ii = 1, addl = 27) +
                       ev(amt = 140, cmt = "GUT_CSA", time = 0, ii = 0.5, addl = 335),
    # 8) Combo: Infliximab + Methotrexate-like Th17 suppression simulated via CsA low-dose
    "Combo_IFX_low_CsA" = ev(amt = 350, cmt = "CENT_IFX", time = 0) +
                          ev(amt = 350, cmt = "CENT_IFX", time = 14) +
                          ev(amt = 350, cmt = "CENT_IFX", time = 42) +
                          ev(amt = 350, cmt = "CENT_IFX", time = 98) +
                          ev(amt = 70, cmt = "GUT_CSA", ii = 0.5, addl = 335),
    # 9) No treatment (natural history)
    "NoTreatment"    = ev(amt = 0, cmt = "GUT_PRED", time = 0)
  )

  out <- model %>%
    ev(ev) %>%
    mrgsim(end = tmax, delta = 0.5)

  as.data.frame(out) %>% mutate(scenario = scenario)
}

# -----------------------------------------------------------------------------
# Run 9 scenarios (use any subset as needed)
# -----------------------------------------------------------------------------
scenarios <- c("NoTreatment", "Prednisone_SOC", "Cyclosporine",
               "Infliximab", "Adalimumab", "Anakinra",
               "Ustekinumab", "Combo_PRED_CSA", "Combo_IFX_low_CsA")

if (interactive() || identical(Sys.getenv("RUN_PG_SIM"), "1")) {
  res <- bind_rows(lapply(scenarios, function(s)
    simulate_pg(pg_qsp_model, s, tmax = 168)))

  message("Endpoint summary @ Week 24:")
  res %>% filter(abs(time - 168) < 0.6) %>%
    select(scenario, Ulcer, Healed, PARACELSUS, DLQI, Pain_VAS, CRP) %>%
    print()
}

# -----------------------------------------------------------------------------
# Virtual population helper — sample CL_ADA, IC50, baseline ulcer area, comorbidities
# -----------------------------------------------------------------------------
make_vpop <- function(n = 200, seed = 42) {
  set.seed(seed)
  data.frame(
    ID            = seq_len(n),
    CL_ADA        = rlnorm(n, log(0.31), 0.25),
    CL_IFX        = rlnorm(n, log(0.32), 0.30),
    CL_CSA        = rlnorm(n, log(27),   0.25),
    UlcerArea_0   = pmax(2, rlnorm(n, log(10), 0.5)),
    DSEV          = pmax(1, rnorm(n, 3.0, 0.7)),
    F_pathergy    = rbinom(n, 1, 0.45),
    F_comorb_IBD  = rbinom(n, 1, 0.40)
  )
}

# Example: run virtual population on adalimumab
run_vpop_scenario <- function(model, scenario = "Adalimumab",
                              n = 50, tmax = 168) {
  vp <- make_vpop(n = n)
  bind_rows(lapply(seq_len(n), function(i) {
    m <- model %>% param(as.list(vp[i, c("CL_ADA","CL_IFX","CL_CSA","UlcerArea_0",
                                         "DSEV","F_pathergy","F_comorb_IBD")]))
    simulate_pg(m, scenario, tmax = tmax) %>% mutate(ID = i)
  }))
}

cat("[PG QSP refactored] Model compiled.  Scenarios available:\n  ",
    paste(scenarios, collapse = "\n   "),
    "\nCall simulate_pg(pg_qsp_model, '<scenario>') to generate trajectories.\n")
