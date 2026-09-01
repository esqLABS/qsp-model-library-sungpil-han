# =============================================================================
# Cluster Headache (CH) — Quantitative Systems Pharmacology (QSP) ODE model
# mrgsolve format
# REFACTORED sibling of ch_mrgsolve_model.R -- pluggable-PK naming convention
# applied to all seven of the file's compounds: Sumatriptan (SUMA), Zolmitriptan
# (ZOL), Verapamil (VERA), Lithium (LI), Topiramate (TOPI), Galcanezumab
# (GALCA), Prednisolone (PRED). See ch_refactor_notes.md for the full account.
# All numeric parameter values are copied verbatim from the original; only PK
# structure/naming and the Hill-interface expression are reorganized. The
# original (ch_mrgsolve_model.R) is untouched.
#
# Build-compat note: the original does not compile under mrgsolve 2.0.1 as
# written (`EMAX_TOPI :  0.35` is missing the required third @annotated
# description field). Per this fork's settled build-compat policy, a purely
# descriptive third field was added to that one line (changing nothing about
# its name, value, or use in any equation) so this file actually compiles and
# runs. Logged as UPSTREAM_ISSUES.md #93. The checked-in original is left
# untouched and still carries the defect.
#
# Scope
#   - Hypothalamic generator (circadian gate) driving trigeminovascular attacks
#   - Attack-rate state controlled by CGRP & PACAP tone; modulated by treatments
#   - PK/PD of 7 drugs : sumatriptan SC, zolmitriptan IN, verapamil PR,
#                        lithium, topiramate, galcanezumab (CGRP mAb), prednisone
#   - Devices/non-PK : O2 (effect compartment), GON block (decaying effect)
#   - Clinical endpoints : attacks/week, mean weekly attack rate, response
#
# Citations used to anchor parameters (see ch_references.md)
#   - Goadsby et al. NEJM 2019 (galcanezumab eCH) -- mean attack reduction 8.7 vs 5.2
#   - Cohen et al. JAMA 2009 (high-flow O2) -- 78% pain-free @15 min
#   - Ekbom et al. NEJM 1991 (SC sumatriptan) -- 74% relief @15 min
#   - Leone et al. Neurology 2000 (verapamil RCT) -- 240 mg/d efficacy
#   - Steiner et al. Cephalalgia 1997 (lithium vs verapamil)
#   - Wei DY, Goadsby PJ. Curr Pain Headache Rep 2021 (pharmacology review)
#
# Compartments (renamed to the fork's pluggable-PK convention; 23 main ODE +
# auxiliary, same count and same 1-based ordering as the original):
#   1  GUT_SUMA       sumatriptan SC depot           (mg)   [was SUMA_DEPOT]
#   2  CENT_SUMA      sumatriptan central             (mg)   [was SUMA_CENT]
#   3  GUT_ZOL        zolmitriptan IN depot            (mg)   [was ZOL_DEPOT]
#   4  CENT_ZOL       zolmitriptan central             (mg)   [was ZOL_CENT]
#   5  PERI_ZOL       zolmitriptan peripheral          (mg)   [was ZOL_PER]
#   6  GUT_VERA       verapamil PR depot               (mg)   [was VERA_DEPOT]
#   7  CENT_VERA      verapamil central                (mg)   [was VERA_CENT]
#   8  PERI_VERA      verapamil peripheral             (mg)   [was VERA_PER]
#   9  GUT_LI         lithium depot                    (mg)   [was LI_DEPOT]
#   10 CENT_LI        lithium central                  (mg)   [was LI_CENT]
#   11 GUT_TOPI       topiramate depot                 (mg)   [was TOPI_DEPOT]
#   12 CENT_TOPI      topiramate central               (mg)   [was TOPI_CENT]
#   13 GUT_GALCA      galcanezumab SC depot            (mg)   [was GALCA_SC]
#   14 CENT_GALCA     galcanezumab central             (mg)   [was GALCA_CENT]
#   15 GUT_PRED       prednisone depot                 (mg)   [was PRED_DEPOT]
#   16 CENT_PRED      prednisolone central             (mg)   [was PRED_CENT]
#   17 HYPO_DRIVE     hypothalamic excitability (0=quiet 1=active) -- unchanged
#   18 CGRP           trigeminal CGRP tone (a.u., reference 1.0) -- unchanged
#   19 PACAP          PACAP / VIP tone (a.u.) -- unchanged
#   20 PIAL           pial vasodilation effect site (a.u.) -- unchanged
#   21 ATTACK_HZ      instantaneous attack hazard /h -- unchanged
#   22 CUM_ATTACKS    cumulative attacks count -- unchanged
#   23 BOUT_TIMER     within-bout time (d) for episodic CH gate -- unchanged
#   24 O2_EFFECT      bolus O2 effect compartment (transient) -- unchanged
#   25 GON_EFF        GON block effect compartment (decays over weeks) -- unchanged
#
# Exposed pluggable-PK interface (naming convention; see ch_refactor_notes.md):
#   C_SUMA/EFFECT_SUMA, C_ZOL/EFFECT_ZOL, C_VERA/EFFECT_VERA, C_LI/EFFECT_LI,
#   C_TOPI/EFFECT_TOPI, C_GALCA/EFFECT_GALCA, C_PRED/EFFECT_PRED -- each a
#   rename of the original's own plain Hill ratio, not a refit (GAMMA_*=1
#   added explicitly where the original's Hill coefficient was implicit).
#   Predeclared in $GLOBAL, not $PARAM: mrgsolve 2.0.1 compiles $PARAM members
#   as read-only references inside $ODE, so a value recomputed every timestep
#   from state cannot also live in $PARAM (same constraint documented in the
#   AAA/AMD/breast-cancer refactors). All fourteen are still visible in every
#   simulation's output via $CAPTURE and in /model_manifest's outputPaths.
#
# Five+ scenarios (see "events" examples at bottom) -- dosing, compartment
# indices (1-based) and timing all unchanged from the original:
#   S0 No treatment
#   S1 O2 + sumatriptan abortive
#   S2 Verapamil 240->480 mg/d preventive
#   S3 Verapamil + lithium (chronic CH)
#   S4 Galcanezumab 300 mg SC monthly
#   S5 Prednisone bridge + verapamil (induction)
#   S6 GON block + verapamil
# =============================================================================

library(mrgsolve)
library(tidyverse)

ch_code <- '
$PROB Cluster Headache QSP — hypothalamic-trigeminovascular axis with 7 drugs (REFACTORED: pluggable-PK naming convention)

$PARAM @annotated
// ---- sumatriptan (SUMA) SC PK -- archetype 3 minus peripheral (depot+central)
KA_SUMA   :  2.5   : SC absorption rate /h
CL_SUMA   :  18    : Sumatriptan CL (L/h)
V1_SUMA   :  9     : Central V (L)
EC50_SUMA :  20    : Conc giving 50% acute abortive effect (ng/mL)
EMAX_SUMA :  0.92  : Max fraction of attack hazard reduction
GAMMA_SUMA:  1.0   : Hill coefficient (none in original; =1)

// ---- zolmitriptan (ZOL) IN PK -- archetype 3 (depot+central+peripheral); IN
//      absorption is a single first-order depot in the original, no biphasic
//      nasal kinetics modeled -- see refactor notes
KA_ZOL    :  1.5   : IN absorption /h
CL_ZOL    :  9     : CL (L/h)
V1_ZOL    :  90    : Central V (L)
Q_ZOL     :  6     : intercompartmental (L/h)
V2_ZOL    :  60    : peripheral V (L)
EMAX_ZOL  :  0.80  : abortive Emax
EC50_ZOL  :  6     : ng/mL
GAMMA_ZOL :  1.0   : Hill coefficient (none in original; =1)

// ---- verapamil (VERA) PR PK -- archetype 3 (depot+central+peripheral)
KA_VERA   :  0.20  : absorption /h
CL_VERA   :  60    : apparent CL (L/h)
V1_VERA   :  300   : central V (L)
Q_VERA    :  35    : intercompartmental (L/h)
V2_VERA   :  500   : peripheral V (L)
EC50_VERA :  90    : ng/mL effective for prevention
EMAX_VERA :  0.55  : max fractional attack-rate reduction
GAMMA_VERA:  1.0   : Hill coefficient (none in original; =1)

// ---- lithium (LI) PK -- archetype 3 minus peripheral (depot+central);
//      CL_LI is renally scaled by CrCL in $MAIN (preserved from original)
KA_LI     :  1.5   : absorption /h
CL_LI     :  1.5   : CL (L/h, ~25 mL/min)
V1_LI     :  45    : V (L)
EC50_LI   :  0.6   : mEq/L
EMAX_LI   :  0.50  : preventive
GAMMA_LI  :  1.0   : Hill coefficient (none in original; =1)

// ---- topiramate (TOPI) PK -- archetype 3 minus peripheral (depot+central)
KA_TOPI   :  0.6   : /h
CL_TOPI   :  1.2   : L/h
V1_TOPI   :  60    : L
EC50_TOPI :  5     : ug/mL
EMAX_TOPI :  0.35  : Max fractional attack-rate reduction (preventive); description field added, see refactor notes / UPSTREAM_ISSUES #93
GAMMA_TOPI:  1.0   : Hill coefficient (none in original; =1)

// ---- galcanezumab (GALCA) PK -- SC mAb; original models this as LINEAR
//      one-compartment-with-depot PK, NOT target-mediated drug disposition
//      (no free-receptor/complex compartment anywhere in the file) --
//      archetype 3 minus peripheral (depot+central)
KA_GALCA   :  0.012 : SC absorption /h (~ F 0.75 t_abs days)
CL_GALCA   :  0.008 : L/h ~ 0.19 L/d
V1_GALCA   :  7.5   : L
EC50_GALCA :  1.5   : ug/mL drug at site, IC50 of CGRP binding
EMAX_GALCA :  0.65  : max attack-rate reduction (eCH trial scale)
GAMMA_GALCA:  1.0   : Hill coefficient (none in original; =1)

// ---- prednisolone (PRED) PK (active metabolite of prednisone prodrug) --
//      archetype 3 minus peripheral (depot+central)
KA_PRED   :  2.5   : /h (rapid)
CL_PRED   :  6     : L/h
V1_PRED   :  40    : L
EC50_PRED :  20    : ng/mL
EMAX_PRED :  0.70  : transitional bridge
GAMMA_PRED:  1.0   : Hill coefficient (none in original; =1)

// ---- disease parameters (hypothalamic generator) -- untouched, disease-side
BASE_HYPO :  0.05  : tonic hypothalamic drive baseline (quiet remission)
BOUT_AMP  :  0.95  : in-bout drive amplitude
KIN_HYPO  :  0.020 : drive build-up /h
KOUT_HYPO :  0.010 : drive decay /h
KIN_CGRP  :  0.15  : CGRP production rate driven by hypothalamus
KOUT_CGRP :  0.30  : CGRP elimination /h
KIN_PACAP :  0.08  : PACAP production /h
KOUT_PACAP:  0.20  : PACAP elimination /h
PIAL_HALF :  3     : pial response half-time (h)
ATTACK_K0 :  0.001 : baseline attack hazard /h (~ 0.024 /d ~ no attacks)
ATTACK_KSAT:  0.18 : max hazard /h ( ~ 4 attacks/day cap )
CGRP_SET  :  1.0   : reference CGRP tone (mAb at this point ~50% effect)
PACAP_SET :  1.0   : reference PACAP tone
KO2_ON    :  3.0   : O2 effect time-constant on (/h)
KO2_OFF   :  6.0   : O2 effect decay /h
GON_HL    :  336   : GON-block effect half-life (h ~ 2 weeks)

// ---- circadian / bout gate
CIRC_AMP  :  0.35  : 24-h amplitude on hazard
CIRC_PHASE:  3     : night peak (~03:00 hr-of-day)
BOUT_LEN  :  42    : days in bout
REMISSION :  240   : days in remission
BOUT_ON   :  1     : start in active bout (1) or not (0)

// ---- patient covariates
CrCL      :  100   : creatinine clearance (mL/min)
WT        :  78    : body weight (kg)
SEX       :  1     : 1=M 0=F
SMOKER    :  1     : 1 active 0 no
CHRONIC   :  0     : 1 chronic CH (no remission)

$CMT @annotated
GUT_SUMA    : sumatriptan SC depot
CENT_SUMA   : sumatriptan central
GUT_ZOL     : zolmitriptan IN depot
CENT_ZOL    : zolmitriptan central
PERI_ZOL    : zolmitriptan peripheral
GUT_VERA    : verapamil PR depot
CENT_VERA   : verapamil central
PERI_VERA   : verapamil peripheral
GUT_LI      : lithium depot
CENT_LI     : lithium central
GUT_TOPI    : topiramate depot
CENT_TOPI   : topiramate central
GUT_GALCA   : galcanezumab SC depot
CENT_GALCA  : galcanezumab central
GUT_PRED    : prednisone depot
CENT_PRED   : prednisolone central
HYPO_DRIVE  : hypothalamic excitability
CGRP        : trigeminal CGRP tone (a.u.)
PACAP       : PACAP tone (a.u.)
PIAL        : pial vasodilation effect site
ATTACK_HZ   : instantaneous attack hazard /h
CUM_ATTACKS : cumulative attacks
BOUT_TIMER  : days within current bout
O2_EFFECT   : O2 effect site
GON_EFF     : GON block effect site

$GLOBAL
// [refactor] exposed effect terms for the 7 refactored compounds,
// predeclared here rather than in $PARAM: mrgsolve 2.0.1 compiles $PARAM
// members as read-only references inside $ODE, so a value that must be
// recomputed every timestep from state cannot also live in $PARAM (same
// constraint documented in the AAA/AMD/breast-cancer refactors). These are
// visible in every simulation output via $CAPTURE and in
// /model_manifest\'s outputPaths.
// [discoverability fix] C_SUMA/C_ZOL/C_VERA/C_LI/C_TOPI/C_GALCA/C_PRED were
// previously predeclared here too (bare-assigned once in $ODE below); moved
// to single-site `double C_<STEM> = <expr>;` initializers in $TABLE instead
// so each is literal-text-discoverable as `double C_<STEM> = ...;` by
// downstream tooling that regexes the corpus for compound PK. Predeclaring
// them here AND initializing with `double` in $TABLE would collide
// ("ambiguous reference") since mrgsolve auto-declares a member from every
// `double NAME = ...;` site found anywhere in the text.
double EFFECT_SUMA, EFFECT_ZOL, EFFECT_VERA, EFFECT_LI, EFFECT_TOPI, EFFECT_GALCA, EFFECT_PRED;

$MAIN
double CL_LI_eff = CL_LI * (CrCL/100);   // renal scaling for lithium (renamed from cl_li)
double bout_drive = BOUT_ON * BOUT_AMP;
// chronic CH never enters remission gate
if(CHRONIC == 1) bout_drive = BOUT_AMP;

$ODE
// concentrations (exposed, PD-facing)
C_SUMA  = CENT_SUMA * 1000 / V1_SUMA;          // ng/mL (mg->ug & ng adj.)
C_ZOL   = CENT_ZOL  * 1000 / V1_ZOL;           // ng/mL
C_VERA  = CENT_VERA * 1000 / V1_VERA;          // ng/mL
C_LI    = CENT_LI   / V1_LI;                   // mEq/L (mg -> mEq approx)
C_TOPI  = CENT_TOPI / V1_TOPI;                 // ug/mL
C_GALCA = CENT_GALCA / V1_GALCA;               // ug/mL
C_PRED  = CENT_PRED * 1000 / V1_PRED;          // ng/mL

// ---- PK ODEs
dxdt_GUT_SUMA  = -KA_SUMA  * GUT_SUMA;
dxdt_CENT_SUMA =  KA_SUMA  * GUT_SUMA - (CL_SUMA / V1_SUMA) * CENT_SUMA;

dxdt_GUT_ZOL   = -KA_ZOL  * GUT_ZOL;
dxdt_CENT_ZOL  =  KA_ZOL  * GUT_ZOL
                  - (CL_ZOL / V1_ZOL) * CENT_ZOL
                  - (Q_ZOL  / V1_ZOL) * CENT_ZOL
                  + (Q_ZOL  / V2_ZOL) * PERI_ZOL;
dxdt_PERI_ZOL  =  (Q_ZOL / V1_ZOL) * CENT_ZOL - (Q_ZOL / V2_ZOL) * PERI_ZOL;

dxdt_GUT_VERA  = -KA_VERA * GUT_VERA;
dxdt_CENT_VERA =  KA_VERA * GUT_VERA
                  - (CL_VERA / V1_VERA) * CENT_VERA
                  - (Q_VERA  / V1_VERA) * CENT_VERA
                  + (Q_VERA  / V2_VERA) * PERI_VERA;
dxdt_PERI_VERA =  (Q_VERA / V1_VERA) * CENT_VERA - (Q_VERA / V2_VERA) * PERI_VERA;

dxdt_GUT_LI    = -KA_LI * GUT_LI;
dxdt_CENT_LI   =  KA_LI * GUT_LI - (CL_LI_eff / V1_LI) * CENT_LI;

dxdt_GUT_TOPI  = -KA_TOPI * GUT_TOPI;
dxdt_CENT_TOPI =  KA_TOPI * GUT_TOPI - (CL_TOPI / V1_TOPI) * CENT_TOPI;

dxdt_GUT_GALCA  = -KA_GALCA * GUT_GALCA;
dxdt_CENT_GALCA =  KA_GALCA * GUT_GALCA - (CL_GALCA / V1_GALCA) * CENT_GALCA;

dxdt_GUT_PRED  = -KA_PRED * GUT_PRED;
dxdt_CENT_PRED =  KA_PRED * GUT_PRED - (CL_PRED / V1_PRED) * CENT_PRED;

// ---- circadian modulation (24-h cycle, fits peak in early hours)
double hr     = fmod(SOLVERTIME, 24.0);
double circ   = 1.0 + CIRC_AMP * cos(2*M_PI*(hr - CIRC_PHASE) / 24.0);

// ---- bout timer (episodic vs chronic)
dxdt_BOUT_TIMER = 1.0/24.0;   // accrue days; reset via $TABLE if needed (here simple)

// ---- hypothalamic drive (target circ * bout_drive)
double hypo_target = BASE_HYPO + circ * bout_drive;
dxdt_HYPO_DRIVE = KIN_HYPO * (hypo_target - HYPO_DRIVE);

// ---- drug effects on hypothalamic drive (preventive); Hill interface,
//      rename not refit -- every original term was already EMAX*C/(C+EC50)
EFFECT_VERA   = EMAX_VERA  * pow(C_VERA,  GAMMA_VERA)  / (pow(EC50_VERA,  GAMMA_VERA)  + pow(C_VERA,  GAMMA_VERA));
EFFECT_LI     = EMAX_LI    * pow(C_LI,    GAMMA_LI)    / (pow(EC50_LI,    GAMMA_LI)    + pow(C_LI,    GAMMA_LI));
EFFECT_TOPI   = EMAX_TOPI  * pow(C_TOPI,  GAMMA_TOPI)  / (pow(EC50_TOPI,  GAMMA_TOPI)  + pow(C_TOPI,  GAMMA_TOPI));
EFFECT_GALCA  = EMAX_GALCA * pow(C_GALCA, GAMMA_GALCA) / (pow(EC50_GALCA, GAMMA_GALCA) + pow(C_GALCA, GAMMA_GALCA));
EFFECT_PRED   = EMAX_PRED  * pow(C_PRED,  GAMMA_PRED)  / (pow(EC50_PRED,  GAMMA_PRED)  + pow(C_PRED,  GAMMA_PRED));
double E_gon    = GON_EFF;                                // 0..1 fraction
double preventive = 1 - (1-EFFECT_VERA)*(1-EFFECT_LI)*(1-EFFECT_TOPI)*(1-EFFECT_GALCA)*(1-EFFECT_PRED)*(1-E_gon);
if(preventive > 0.95) preventive = 0.95;

// ---- trigeminovascular tone driven by hypothalamus, dampened by preventive
double drive_eff = HYPO_DRIVE * (1 - preventive);
dxdt_CGRP  = KIN_CGRP  * drive_eff - KOUT_CGRP  * CGRP;
dxdt_PACAP = KIN_PACAP * drive_eff - KOUT_PACAP * PACAP;

// ---- pial response (effect compartment of CGRP+PACAP)
double pial_drive = 0.6*CGRP + 0.4*PACAP;
dxdt_PIAL  = (log(2)/PIAL_HALF) * (pial_drive - PIAL);

// ---- O2 effect decay
dxdt_O2_EFFECT = -KO2_OFF * O2_EFFECT;

// ---- GON block decay (single exponential, k = ln2/HL)
dxdt_GON_EFF = -(log(2)/GON_HL) * GON_EFF;

// ---- attack hazard (CGRP-driven, capped, acute abortive abates)
EFFECT_SUMA = EMAX_SUMA * pow(C_SUMA, GAMMA_SUMA) / (pow(EC50_SUMA, GAMMA_SUMA) + pow(C_SUMA, GAMMA_SUMA));
EFFECT_ZOL  = EMAX_ZOL  * pow(C_ZOL,  GAMMA_ZOL)  / (pow(EC50_ZOL,  GAMMA_ZOL)  + pow(C_ZOL,  GAMMA_ZOL));
double E_O2_acute   = O2_EFFECT;                          // 0..1
double acute_abort  = 1 - (1 - EFFECT_SUMA)*(1 - EFFECT_ZOL)*(1 - E_O2_acute);
if(acute_abort > 0.99) acute_abort = 0.99;

double hz_raw = ATTACK_KSAT * (PIAL / (PIAL + 1.0));
double hz_now = (ATTACK_K0 + hz_raw) * (1 - acute_abort);
dxdt_ATTACK_HZ = (hz_now - ATTACK_HZ) * 5.0;   // smoothing

dxdt_CUM_ATTACKS = ATTACK_HZ;                  // integrate hazard

$TABLE
capture ConcSuma    = C_SUMA;
capture ConcZol     = C_ZOL;
capture ConcVera    = C_VERA;
capture ConcLi      = C_LI;
capture ConcTopi    = C_TOPI;
capture ConcGalca   = C_GALCA;
capture ConcPred    = C_PRED;
capture HypoDrive   = HYPO_DRIVE;
capture CGRPtone    = CGRP;
capture PACAPtone   = PACAP;
capture Pialtone    = PIAL;
capture HazardPerH  = ATTACK_HZ;
capture AttacksWeek = ATTACK_HZ * 168.0;
capture Preventive  = preventive;
capture AcuteAbort  = acute_abort;

// [discoverability fix] single contiguous `double C_<STEM> = <expr>;`
// initializers, identical formulas to the $ODE bare assignments above, so
// each compound is literal-text-discoverable by downstream tooling. By the
// time $TABLE runs, $ODE has already bare-assigned the current-timestep
// value into these same names (no longer $GLOBAL-predeclared -- see
// $GLOBAL comment), so this recompute is not stale and the ConcXxx
// captures above (which read C_<STEM> before this point) already saw the
// correct current value.
double C_SUMA  = CENT_SUMA * 1000 / V1_SUMA;
double C_ZOL   = CENT_ZOL  * 1000 / V1_ZOL;
double C_VERA  = CENT_VERA * 1000 / V1_VERA;
double C_LI    = CENT_LI   / V1_LI;
double C_TOPI  = CENT_TOPI / V1_TOPI;
double C_GALCA = CENT_GALCA / V1_GALCA;
double C_PRED  = CENT_PRED * 1000 / V1_PRED;

$CAPTURE
C_SUMA C_ZOL C_VERA C_LI C_TOPI C_GALCA C_PRED
EFFECT_SUMA EFFECT_ZOL EFFECT_VERA EFFECT_LI EFFECT_TOPI EFFECT_GALCA EFFECT_PRED
'

# ---- compile ----------------------------------------------------------------
ch_mod <- mcode("cluster_headache_qsp_refactored", ch_code)

# ============================================================================
# Helper builders for the 6 standard scenarios (compartment names renamed;
# dosing amounts/timing/routes unchanged from the original)
# ============================================================================

build_scenario <- function(scenario = "S0") {
  ev <- ev(amt = 0, cmt = 1, time = 0)   # placeholder

  if (scenario == "S0") {
    # no drug
    return(ev)
  }
  if (scenario == "S1") {
    # 1 attack at day 1 morning: 6 mg sumatriptan SC + O2 bolus (set effect to 0.78)
    e1 <- ev(time = 24 + 4,  amt = 6,    cmt = "GUT_SUMA")
    e2 <- ev(time = 24 + 4,  amt = 0.78, cmt = "O2_EFFECT", evid = 1)  # bolus
    return(c(e1, e2))
  }
  if (scenario == "S2") {
    # verapamil 240 mg PR BID for 14 d, then 480 mg/d (i.e. 240 BID)
    e1 <- ev(amt = 240, cmt = "GUT_VERA", ii = 12, addl = 27, time = 0)
    return(e1)
  }
  if (scenario == "S3") {
    # verapamil + lithium 600 mg/d split BID
    e1 <- ev(amt = 240, cmt = "GUT_VERA", ii = 12, addl = 55, time = 0)
    e2 <- ev(amt = 300, cmt = "GUT_LI",   ii = 12, addl = 55, time = 0)
    return(c(e1, e2))
  }
  if (scenario == "S4") {
    # galcanezumab 300 mg SC monthly x 3
    e1 <- ev(amt = 300, cmt = "GUT_GALCA", ii = 24*28, addl = 2, time = 0)
    return(e1)
  }
  if (scenario == "S5") {
    # prednisone 60 mg/d x 5 d -> taper -> verapamil 240 mg BID
    e1 <- ev(amt = 60, cmt = "GUT_PRED", ii = 24, addl = 4, time = 0)
    e2 <- ev(amt = 40, cmt = "GUT_PRED", ii = 24, addl = 2, time = 24*5)
    e3 <- ev(amt = 20, cmt = "GUT_PRED", ii = 24, addl = 2, time = 24*8)
    e4 <- ev(amt = 240, cmt = "GUT_VERA", ii = 12, addl = 50, time = 0)
    return(c(e1, e2, e3, e4))
  }
  if (scenario == "S6") {
    # GON block (one-shot effect 0.65) + verapamil
    e1 <- ev(time = 0, amt = 0.65, cmt = "GON_EFF", evid = 1)
    e2 <- ev(amt = 240, cmt = "GUT_VERA", ii = 12, addl = 27, time = 0)
    return(c(e1, e2))
  }
  ev
}

# ============================================================================
# Example simulation: run all six scenarios (12-week horizon)
# ============================================================================

simulate_all <- function(weeks = 12) {
  scenarios <- c("S0","S1","S2","S3","S4","S5","S6")
  out_list <- list()
  for (sc in scenarios) {
    ev_sc <- build_scenario(sc)
    out <- ch_mod %>%
      ev(ev_sc) %>%
      mrgsim(end = 24*7*weeks, delta = 1) %>%
      as_tibble() %>%
      mutate(scenario = sc)
    out_list[[sc]] <- out
  }
  bind_rows(out_list)
}

# ============================================================================
# Quick sanity helpers ---------------------------------------------------------
# (run with: source("ch_mrgsolve_model_refactored.R"); simulate_all() %>% ...)
# ============================================================================
