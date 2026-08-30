# =============================================================================
#  Essential Tremor (본태성 떨림) — QSP model for mrgsolve
#  ---------------------------------------------------------------------------
#  48 differential equations.  22 therapeutic scenarios.
#
#  THE ONE IDEA THIS FILE IS BUILT ON
#  ---------------------------------------------------------------------------
#  A tremor is not a LEVEL of anything.  It is a LIMIT CYCLE.  So the state
#  variable is not "tremor" — it is the AMPLITUDE ENVELOPE of an oscillator
#  sitting just past a supercritical Hopf bifurcation:
#
#        G_total = G0*(1+PROG) * [ w_C*Phi_C + w_P*Phi_P ]
#        Phi_C   = (a_O*phi_olive + a_R*phi_cblthal) * phi_thal * phi_ctx
#        Phi_P   = phi_spindle * phi_nmj
#        mu      = G_total - 1
#        dr/dt   = (1/tauA) * ( r*(mu - beta*r^2) + eps ) / (1 + |mu|)
#        r*      = sqrt(mu/beta)     for mu > 0 ;  ~0 for mu < 0
#
#  Three consequences do all the work, and none of them is asserted anywhere —
#  they are algebraic:
#
#   (1) AMPLITUDE comes from GAIN, FREQUENCY comes from DELAY.  Drugs enter the
#       phi_* terms (gain) and never touch tau_loop (delay), so every drug in
#       this file changes amplitude and leaves frequency alone — which is what
#       is observed.  Frequency falls with age and Purkinje loss instead,
#       because those lengthen the loop.
#
#   (2) Because r* = sqrt(mu), the SAME pharmacology abolishes a mild tremor
#       and barely touches a severe one.  The responder/non-responder split in
#       every beta-blocker trial ever run is one equation crossing threshold,
#       not two populations of biology.  And for the same reason two drugs are
#       SUPRA-additive, because sqrt is concave and the second drug acts where
#       the map is steepest.
#
#   (3) The therapeutic CEILING is set by TOPOLOGY, not potency.  The inferior
#       olive is one of two PARALLEL branches of the central loop (weight a_O);
#       the Vim relay is a SERIES element.  So a perfect Cav3.1 blocker cannot
#       do better than (a_O*0 + a_R) = a_R, while a Vim lesion multiplies the
#       whole central limb by ~0.08.  That is why surgery beats drugs, and it
#       is a statement about wiring, not about drug quality.
#
#  A fourth result is about measurement rather than biology: rating scales are
#  LOGARITHMIC in amplitude (Elble, Brain 2006), rating = a + 2*log10(A).  A
#  50% accelerometric reduction is therefore 0.60 points per item, forever, at
#  any baseline.  The "discordance" between accelerometry and TETRAS in the
#  literature is not a discordance.  It is the logarithm.
#
#  VERIFICATION
#  ---------------------------------------------------------------------------
#  Every equation below was independently re-implemented in dependency-free
#  Python (et_verify.py) and integrated with a hand-written RK4.  That exercise
#  found and fixed five real defects in this model, each marked [FIXED] in
#  place.  It also refuted the author's prior expectation about combination
#  therapy (see scenario 7).  et_verify.py also checks that the slow-envelope
#  reduction is legitimate by integrating the full 2-D oscillator at 5.5 Hz and
#  comparing its peak amplitude with r* (agreement 0.06-2.0% over mu 0.15-5).
#
#  CAVEAT.  This is an educational/research model.  It is not fitted to patient
#  data, not validated, and must not be used for clinical decisions.
# =============================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# Essential Tremor QSP — oscillator formulation
# 48 ODEs: 27 PK/effect-site, 9 adaptation/disease, 3 tremor envelopes,
#          6 organ-system, plus 3 slow markers.

$PARAM @annotated
// ---------------- oscillator core -----------------------------------------
G0      : 1.60  : Untreated loop gain (mild 1.15, moderate 1.60, severe 6.0)
WC      : 0.60  : Central loop weight (WC+WP=1)
WP      : 0.40  : Peripheral stretch-reflex loop weight
AO      : 0.35  : Olivary share of central drive (PARALLEL branch)
AR      : 0.65  : Cerebello-thalamic share (AO+AR=1)
BETA    : 1.00  : Hopf saturation coefficient
TAUA    : 0.35  : h : tremor envelope time constant
// [REFACTOR build-compat] renamed EPS -> EPS_NOISE: "EPS" is a reserved word
// under mrgsolve 2.0.1 (the untouched original fails to compile with
// "Reserved words in model names: EPS"); see et_refactor_notes.md and
// UPSTREAM_ISSUES.md #102. Non-numeric, name-only.
EPS_NOISE     : 1e-4  : Noise seeding of the limit cycle
KPROG   : 5.71e-6 : 1/h : progression of loop gain (G0 doubles in 20 y)
KEXC    : 0.0   : Hypothesis switch: oscillation drives cerebellar injury
KDEGPC  : 1.427e-6 : 1/h : Purkinje integrity decline

// ---------------- effector gain shares ------------------------------------
HDG     : 0.55  : Neck effector gain share (0.85 = head-tremor phenotype)
VXG     : 0.45  : Laryngeal effector gain share
WC_HD   : 0.85  : Neck loop central weight
WP_HD   : 0.15  : Neck loop peripheral weight
NECKBTX : 0.0   : Neck botulinum effect (0-1)

// ---------------- biomechanics --------------------------------------------
J0      : 0.00256 : kg.m2 : hand inertia about the wrist
KST     : 7.30  : N.m/rad : wrist stiffness (f0 = 8.50 Hz unloaded)
ZETA0   : 0.35  : damping RATIO, held constant under load (co-contraction)
BD      : 0.0955 : N.m.s/rad : wrist damping coefficient (documentation only)
MLOAD   : 0.0   : kg : added mass (weight-loading test)
LLOAD   : 0.10  : m : moment arm of the added mass
KAMP    : 1.902 : cm : amplitude calibration
KNZ     : 0.021 : cm : physiological-tremor noise gain
KCAT    : 12.0  : Adrenergic amplification of physiological tremor
FNZ     : 0.55  : beta2-dependent fraction of physiological tremor
J0_HD   : 0.0250 : kg.m2 : head inertia
KST_HD  : 8.80  : N.m/rad : neck stiffness (f0 = 2.99 Hz)
BD_HD   : 0.329 : N.m.s/rad
KAMP_HD : 3.10  : deg : head amplitude calibration
ASYM    : 0.80  : non-dominant / dominant amplitude ratio

// ---------------- loop delay / frequency ----------------------------------
TAU0    : 0.18182 : s : loop delay at age 50 (-> 5.50 Hz)
AGE     : 60.0  : y
KAGE    : 0.110 : per decade above 50 (-> -0.055 Hz/y)
KPC     : 0.42  : Purkinje loss lengthens the loop
KLOADN  : 0.012 : mass loading effect on the CENTRAL delay (small on purpose)
KVISD   : 0.055 : visuomotor delay added during intention

// ---------------- peripheral loop / beta2 ---------------------------------
FB2     : 0.60  : beta2-dependent fraction of peripheral loop gain
ADR0    : 0.20  : nM : resting adrenaline equivalent
KD_AG   : 0.35  : nM : adrenaline at spindle beta2 (resting occupancy 0.36)
STRESS  : 0.0   : multiplier on adrenergic drive
CAFF    : 0.0   : caffeine
THYRO   : 0.0   : thyrotoxicosis
SALB    : 0.0   : nM : beta2-agonist equivalent
SF      : 3.0   : NMJ safety factor
HNMJ    : 2.5   : supralinearity of release on SNAP-25

// ---------------- disease / adaptation ------------------------------------
KG_OL   : 0.70  : GABA-A potentiation reaching the olive
KREBF   : 2.20  : acute counter-adaptation -> rebound gain
KREBS   : 1.20  : chronic counter-adaptation -> baseline shift
KAF     : 0.90  : fast adaptation gain
TAUF_ON : 1.0   : h : acute tolerance onset (Mellanby, fast)
TAUF_OFF: 5.0   : h : acute tolerance offset (slow -> this asymmetry IS the rebound)
KAS     : 0.90  : slow adaptation gain (ethanol only)
TAUS    : 720.0 : h
HARM    : 0.0   : harmaline drive on the olive
KRR     : 3.70  : thalamic re-routing gain
REROUTE_MAX : 0.45 : re-routing capacity (patient-level, 1.0 = relapser)
TAU_RR  : 13140 : h

// ---------------- Vim lesion / DBS ----------------------------------------
LESION  : 0.0   : lesion switch
VLES    : 0.0   : mm3 : lesion volume
V50L    : 45.0  : mm3 : lesion efficacy V50
HLES    : 2.5   : lesion efficacy Hill
V50A    : 260.0 : mm3 : lesion ataxia V50
AXL     : 1.00  : lesion ataxia Emax
DBSON   : 0.0   : DBS switch
FSTIM   : 130.0 : Hz
VTA     : 250.0 : mm3 : volume of tissue activated
F50D    : 80.0  : Hz : frequency threshold
HDBS    : 4.0   : frequency Hill
V50D    : 90.0  : mm3
// [REFACTOR build-compat] added the missing ": description" annotation field
// below: the untouched original's "EMAX_DBS: 0.90" (two colon-separated
// fields) fails mrgsolve 2.0.1's @annotated $PARAM parser ("improper
// annotation format"), which requires exactly three fields. Non-numeric,
// description-only; see et_refactor_notes.md and UPSTREAM_ISSUES.md #102.
EMAX_DBS: 0.90 : DBS max block fraction
KENT    : 0.30  : low-frequency entrainment gain
F50E    : 60.0  : Hz : entrainment cutoff
BILAT   : 0.0   : bilateral flag

// ---------------- propranolol --------------------------------------------

// ---------------- atenolol / nadolol -------------------------------------

// ---------------- primidone / phenobarbital / PEMA -----------------------

// ---------------- topiramate / gabapentin / benzodiazepine ---------------

// ---------------- ethanol / 1-octanol ------------------------------------

// ---------------- T-type calcium blocker --------------------------------

// ---------------- botulinum toxin ---------------------------------------

// ---------------- organ systems -----------------------------------------

// ---------------- task / state ------------------------------------------


$PARAM
// Drug PK/PD parameters.  Kept in a plain $PARAM rather than @annotated
// because the annotated form requires exactly one parameter per line.
  MW_PRP = 259.3,
  FU_PRP = 0.10,
  KA_PRP = 0.35,
  V1_PRP = 250,
  V2_PRP = 300,
  CL_PRP = 60,
  Q_PRP = 60,
  F0_PRP = 0.20,
  FMX_PRP = 0.25,
  FD50_PRP = 150,
  DPRP = 0.0,
  KPUU_PRP = 1.50,
  KE0_PRP = 1.20,
  KI_PRP_B1 = 1.80,
  KI_PRP_B2 = 0.60,
  EMAX_PRP_CBL = 0.25,
  EC50_PRP_CBL = 150.0,
  KB2REG = 0.55,
  TAU_B2REG = 336.0,
  MW_ATN = 266.3,
  FU_ATN = 0.95,
  KA_ATN = 0.80,
  V1_ATN = 70,
  CL_ATN = 10,
  F_ATN = 0.50,
  KI_ATN_B1 = 30.0,
  KI_ATN_B2 = 1000.0,
  MW_NAD = 309.4,
  FU_NAD = 0.75,
  KA_NAD = 0.40,
  V1_NAD = 140,
  CL_NAD = 18,
  F_NAD = 0.32,
  KI_NAD_B1 = 3.00,
  KI_NAD_B2 = 1.20,
  MW_PRM = 218.25,
  MW_PB = 232.24,
  MW_PEM = 190.20,
  KA_PRM = 1.50,
  F_PRM = 0.92,
  V1_PRM = 40,
  CL_PRM = 2.80,
  FM_PB = 0.25,
  FM_PEM = 0.45,
  V1_PB = 45,
  CL_PB = 0.32,
  V1_PEM = 45,
  CL_PEM = 2.00,
  KP_PRM = 0.90,
  KP_PB = 0.70,
  KE0_PRM = 1.00,
  KE0_PB = 0.35,
  EMAX_PRM_GABA = 0.62,
  EC50_PRM_GABA = 28.0,
  EMAX_PB_GABA = 0.55,
  EC50_PB_GABA = 220.0,
  EMAX_PEM_GABA = 0.10,
  EC50_PEM_GABA = 400.0,
  EMAX_PRM_NACH = 0.12,
  EC50_PRM_NACH = 30.0,
  KIND = 1.30,
  TAU_IND = 336.0,
  MW_TOP = 339.4,
  KA_TOP = 1.00,
  F_TOP = 0.80,
  V1_TOP = 60,
  CL_TOP = 1.40,
  KP_TOP = 0.90,
  EMAX_TOP_CTX = 0.28,
  EC50_TOP_CTX = 14.0,
  EMAX_TOP_GABA = 0.15,
  EC50_TOP_GABA = 30.0,
  MW_GBP = 171.2,
  KA_GBP = 1.20,
  V1_GBP = 60,
  CL_GBP = 10,
  KP_GBP = 0.15,
  EMAX_GBP = 0.18,
  EC50_GBP = 3.0,
  EMAX_BZD = 0.30,
  BZDLEV = 0.0,
  V1_ETH = 47.6,
  KA_ETH = 4.00,
  VMAX_ETH = 7.14,
  KM_ETH = 0.08,
  KE0_ETH = 6.00,
  KP_ETH = 1.00,
  EMAX_ETH_GABA = 0.65,
  EC50_ETH_GABA = 0.45,
  EMAX_ETH_OL = 0.75,
  EC50_ETH_OL = 0.30,
  EMAX_ETH_CTX = 0.10,
  EC50_ETH_CTX = 0.90,
  KPCTX_ETH = 1.00,
  EC50_INTOX_ETH = 0.80,
  ETHCHR = 0.0,
  MW_OCT = 130.2,
  KA_OCT = 2.00,
  V1_OCT = 50,
  CL_OCT = 25,
  KP_OCT = 1.00,
  EMAX_OCT_GABA = 0.35,
  EC50_OCT_GABA = 120.0,
  EMAX_OCT_OL = 0.60,
  EC50_OCT_OL = 60.0,
  KPCTX_OCT = 0.15,
  EC50_INTOX_OCT = 250.0,
  MW_TTB = 400.0,
  KA_TTB = 1.20,
  F_TTB = 0.70,
  V1_TTB = 200,
  CL_TTB = 15,
  FU_TTB = 0.15,
  KPUU_TTB = 1.50,
  KE0_TTB = 0.80,
  EMAX_TTB_OL = 0.85,
  EC50_TTB_OL = 80.0,
  KDEG_BTX = 0.35,
  KCL_BTX = 0.00328,
  KR_SNAP = 3.21e-4,
  FSPILL = 0.15,
  HR0 = 72,
  EHR = 0.31,
  SBP0 = 132,
  ESBP = 14,
  FEV10 = 3.20,
  ASTHMA = 0.0,
  EFEV_A = 0.30,
  EFEV_N = 0.04,
  TAU_FEV = 2.0,
  EC50_SED_PB = 90.0,
  EC50_SED_PRM = 45.0,
  TAU_SEDTOL = 240.0,
  KSEDTOL = 0.80,
  HCO30 = 24.0,
  EHCO3 = 5.0,
  EC50_HCO3 = 20.0,
  TAU_HCO3 = 48.0,
  BW0 = 78.0,
  EBW = 0.09,
  EC50_BW = 25.0,
  TAU_BW = 1440.0,
  ECOG = 45.0,
  EC50_COG = 30.0,
  TAU_COG = 336.0,
  BMD0 = 1.00,
  KBMD = 1.9e-6,
  ALT0 = 25.0,
  EALT = 45.0,
  TAU_ALT = 720.0,
  TASK_INT = 0.0,
  KVIS = 0.22,
  FATIGUE = 0.0,
  KFAT = 0.10,
  // ---- refactor-added: explicit Hill exponents (gamma=1, none explicit
  //      in the original) and promoted magic-number EC50s (same values as
  //      the original's inline literals, now named/discoverable) ----
  GAMMA_PRP_CBL = 1.0,
  GAMMA_PRM_GABA = 1.0,
  GAMMA_PB_GABA = 1.0,
  GAMMA_PEM_GABA = 1.0,
  GAMMA_PRM_NACH = 1.0,
  GAMMA_TOP_CTX = 1.0,
  GAMMA_TOP_GABA = 1.0,
  GAMMA_GBP = 1.0,
  GAMMA_ETH_GABA = 1.0,
  GAMMA_ETH_OL = 1.0,
  GAMMA_ETH_CTX = 1.0,
  GAMMA_OCT_GABA = 1.0,
  GAMMA_OCT_OL = 1.0,
  GAMMA_TTB_OL = 1.0,
  EC50_SED_PRP = 300.0,
  EC50_SED_GBP = 4.0,
  EC50_SED_TOP = 30.0,
  EC50_SED_ETH = 0.90,
  EC50_ATAX_PB = 120.0,
  EC50_ATAX_ETH = 1.20

$GLOBAL
// mrgsolve compiles $MAIN, $ODE and $TABLE as separate function bodies, so a
// variable declared `double x = ...` inside $ODE is NOT visible to $TABLE and
// cannot be named in $CAPTURE.  Almost everything in this model is shared that
// way -- the gain terms are computed in $ODE and $TABLE is what turns them into
// clinical scales -- so every shared quantity is declared once at file scope
// and only assigned in the blocks below.
#define TCLAMP(x) ((x) < 0.0 ? 0.0 : ((x) > 4.0 ? 4.0 : (x)))
// [REFACTOR build-compat] added F_PRP to this declaration list: the untouched
// original computes and assigns F_PRP in $MAIN ("F_PRP = F0_PRP + ...") but
// never declares it anywhere, so it fails to compile ("'F_PRP' was not
// declared in this scope"). This is inside propranolol's own block (the
// refactor's own scope), so it is fixed here as part of the archetype work,
// not as a foreign patch; still logged as a pre-existing defect in the
// checked-in original -- see et_refactor_notes.md and UPSTREAM_ISSUES.md #102.
double PC0, G00, MU00, GH0, GV0, C_PRP, F_PRP;
double CF_PRP, C_ATN, CF_ATN, C_NAD, CF_NAD, C_PRM;
double CU_PRM, C_PB, CU_PB, C_PEM, CU_PEM, C_TOP;
double CU_TOP, CB_TOP, C_GBP, CU_GBP, CB_GBP, C_ETH;
double C_OCT, CU_OCT, CB_OCT, C_TTB, CF_TTB, IQ_B2;
double IQ_B1, OCCB2, OCCB1, AG, OCC_AG, OCC_AG0;
double RAG, PHI_SPIN, dnm, PHI_NMJ, GRIP, PHI_P;
double EFFECT_PRM_GABA, EFFECT_PB_GABA, EFFECT_PEM_GABA, EFFECT_ETH_GABA, EFFECT_OCT_GABA, EFFECT_TOP_GABA;
double P_BZD, P_RAW, P_EFF, P_ETHR, reb_f, reb_s;
double REB, EFFECT_TTB_OL, ETH_OL, PHI_OL, PHI_CBL, LES_EFF;
double DBSF, DBSV, DBSB, eterm, ENTR, PHI_TH;
double PHI_CTX, PHI_C, G0T, GTOT, MU, GHD;
double MU_HD, GVX, MU_VX, JREL, TAU_L, FNEUR;
double JTOT, F0M, ZETA, rN, HN_, HP_;
double NZ, A_PHYS, A_LC, A_UL, F_OBS, F0_HD;
double Z_HD, rH_, H_HD, A_HD, A_VX, KE_PRM;
double tgtF, SED_RAW, rsq, supp, iT, FEV_SS;
double HCSS, BWSS, CGSS, HR, SBP, SED;
double ATAXd, ATAX, INTOXd, INTOX, T_R, T_L;
double T_HD, T_VX, UL_T, SPI, HW, LL;
double FACE, TONG, STAND, TETRAS_PS, GRIPLOSS, ADLf;
double TETRAS_ADL, FTM, SPIRAL, BAINF, QUESTd, QUEST;
// ---- refactor-added: per-drug named Hill/quotient effect terms that
//      were previously computed inline (never their own named variable)
double EFFECT_PRP_B2, EFFECT_PRP_B1, EFFECT_ATN_B2, EFFECT_ATN_B1, EFFECT_NAD_B2, EFFECT_NAD_B1;
double EFFECT_PRP_CBL, EFFECT_PRM_NACH, EFFECT_GBP, EFFECT_ETH_CTX, EFFECT_TOP_CTX, EFFECT_ETH_OL, EFFECT_OCT_OL;

$CMT @annotated
// ---- propranolol / atenolol / nadolol
GUT_PRP : Propranolol gut (mg)
CENT_PRP : Propranolol central (mg)
PERI_PRP : Propranolol peripheral (mg)
CE_PRP : Propranolol free brain (nM)
GUT_ATN : Atenolol gut (mg)
CENT_ATN : Atenolol central (mg)
GUT_NAD : Nadolol gut (mg)
CENT_NAD : Nadolol central (mg)
// ---- primidone family
GUT_PRM : Primidone gut (mg)
CENT_PRM : Primidone central (mg)
CENT_PB  : Phenobarbital central (mg)
CENT_PEM : PEMA central (mg)
CE_PRM : Primidone brain (uM)
CE_PB  : Phenobarbital brain (uM)
// ---- other orals
GUT_TOP : Topiramate gut (mg)
CENT_TOP : Topiramate central (mg)
GUT_GBP : Gabapentin gut (mg)
CENT_GBP : Gabapentin central (mg)
// ---- ethanol / octanol
GUT_ETH : Ethanol gut (g)
CENT_ETH : Ethanol body water (g)
CE_ETH : Ethanol brain (g/L)
GUT_OCT : Octanol gut (mg)
CENT_OCT : Octanol central (mg)
// ---- T-type blocker
GUT_TTB : T-type blocker gut (mg)
CENT_TTB : T-type blocker central (mg)
CE_TTB : T-type blocker free brain (nM)
// ---- botulinum
A_BTXT : BoNT depot, tremor-dominant muscle (U)
A_BTXG : BoNT depot, grip muscle (U)
SNAPT  : Functional SNAP-25, tremor muscle (fraction)
SNAPG  : Functional SNAP-25, grip muscle (fraction)
// ---- adaptation / disease
ADAPTF  : Fast GABA-A counter-adaptation
ADAPTS  : Slow GABA-A tolerance
B2REG   : beta-adrenoceptor up-regulation
IND     : Hepatic enzyme induction
SEDTOL  : Sedation tolerance
PCINT   : Purkinje cell structural integrity
DNDIS   : Dentate nucleus disinhibition
REROUTE : Thalamic re-routing (habituation)
PROG    : Loop-gain progression
// ---- tremor envelopes (the oscillator)
R_UL : Upper-limb tremor envelope
R_HD : Head tremor envelope
R_VX : Voice tremor envelope
// ---- organ systems
FEV1 : Forced expiratory volume (L)
HCO3 : Serum bicarbonate (mmol/L)
BW   : Body weight (kg)
BMD  : Bone mineral density index
ALT  : Alanine aminotransferase (U/L)
COG  : Word-finding impairment (0-100)

$MAIN
SNAPT_0  = 1.0;   SNAPG_0 = 1.0;
PCINT_0  = 1.0;
FEV1_0   = FEV10; HCO3_0  = HCO30;
BW_0     = BW0;   BMD_0   = BMD0;  ALT_0 = ALT0;
// seat the tremor envelopes at the untreated fixed point so that t=0 is a
// steady state rather than a transient
PC0  = AO + AR;                       // = 1 by construction
G00  = G0*(WC*PC0 + WP*1.0);
MU00 = G00 - 1.0;
R_UL_0 = (MU00 > 0) ? sqrt(MU00/BETA) : 1e-3;
GH0  = G0*HDG*(WC_HD*PC0 + WP_HD*1.0) - 1.0;
R_HD_0 = (GH0  > 0) ? sqrt(GH0/BETA)  : 1e-3;
GV0  = G0*VXG*(0.90*PC0 + 0.10*1.0) - 1.0;
R_VX_0 = (GV0  > 0) ? sqrt(GV0/BETA)  : 1e-3;
// dose-dependent first-pass escape for propranolol (saturable CYP2D6/1A2)
F_PRP = F0_PRP + FMX_PRP*DPRP/(FD50_PRP + DPRP);
F_GUT_PRP = F_PRP;
F_GUT_ATN = F_ATN;
F_GUT_NAD = F_NAD;
F_GUT_PRM = F_PRM;
F_GUT_TOP = F_TOP;

$ODE
// =========================================================================
// A. CONCENTRATIONS
// =========================================================================
C_PRP  = CENT_PRP/V1_PRP;                       // mg/L total
CF_PRP = C_PRP*FU_PRP*1e6/MW_PRP;             // nM free
C_ATN  = CENT_ATN/V1_ATN;
CF_ATN = C_ATN*FU_ATN*1e6/MW_ATN;
C_NAD  = CENT_NAD/V1_NAD;
CF_NAD = C_NAD*FU_NAD*1e6/MW_NAD;
C_PRM  = CENT_PRM/V1_PRM;  CU_PRM = C_PRM*1e3/MW_PRM;   // uM
C_PB   = CENT_PB/V1_PB;    CU_PB  = C_PB*1e3/MW_PB;
C_PEM  = CENT_PEM/V1_PEM;  CU_PEM = C_PEM*1e3/MW_PEM;
C_TOP  = CENT_TOP/V1_TOP;  CU_TOP = C_TOP*1e3/MW_TOP;
CB_TOP = KP_TOP*CU_TOP;
C_GBP  = CENT_GBP/V1_GBP;  CU_GBP = C_GBP*1e3/MW_GBP;
CB_GBP = KP_GBP*CU_GBP;
C_ETH  = CENT_ETH/V1_ETH;                        // g/L
C_OCT  = CENT_OCT/V1_OCT;  CU_OCT = C_OCT*1e3/MW_OCT;
CB_OCT = KP_OCT*CU_OCT;
C_TTB  = CENT_TTB/V1_TTB;
CF_TTB = C_TTB*FU_TTB*1e6/MW_TTB;

// =========================================================================
// B. beta-ADRENOCEPTOR OCCUPANCY  (Gaddum competitive antagonism)
//    ONE number, OCCB2, will turn out to carry both the tremor benefit and
//    the asthma contraindication.  That is not a coincidence in this model,
//    it is the same term appearing in two places.
// =========================================================================
// Each drug's own named competitive-occupancy quotient (bespoke: a linear
// Ci/Ki term feeding a SHARED Gaddum competitive-antagonism equilibrium at
// the same receptor, not an independent per-drug Hill saturation -- see
// refactor notes). Combined into the shared occupancy only at this point.
EFFECT_PRP_B2 = CF_PRP/KI_PRP_B2;
EFFECT_ATN_B2 = CF_ATN/KI_ATN_B2;
EFFECT_NAD_B2 = CF_NAD/KI_NAD_B2;
EFFECT_PRP_B1 = CF_PRP/KI_PRP_B1;
EFFECT_ATN_B1 = CF_ATN/KI_ATN_B1;
EFFECT_NAD_B1 = CF_NAD/KI_NAD_B1;
IQ_B2 = EFFECT_PRP_B2 + EFFECT_ATN_B2 + EFFECT_NAD_B2;
IQ_B1 = EFFECT_PRP_B1 + EFFECT_ATN_B1 + EFFECT_NAD_B1;
OCCB2 = IQ_B2/(1.0 + IQ_B2);
OCCB1 = IQ_B1/(1.0 + IQ_B1);
// [FIXED — defect 3 found by et_verify.py] B2REG was originally multiplying the
// AGONIST concentration, which made atenolol make tremor WORSE on treatment
// (+2.5%).  Up-regulation is a receptor-NUMBER effect, so it belongs on the
// beta2-mediated gain term, where it is silent under blockade and only appears
// as withdrawal rebound once the antagonist washes out.
AG      = ADR0*(1.0 + STRESS + CAFF + THYRO) + SALB;
OCC_AG  = AG/(KD_AG*(1.0 + IQ_B2) + AG);
OCC_AG0 = ADR0/(KD_AG + ADR0);
RAG     = OCC_AG/OCC_AG0;

// =========================================================================
// C. PERIPHERAL STRETCH-REFLEX LOOP
// =========================================================================
PHI_SPIN = 1.0 - FB2 + FB2*(1.0 + B2REG)*RAG;
dnm      = 1.0 - exp(-SF);
PHI_NMJ  = (1.0 - exp(-SF*pow(SNAPT, HNMJ)))/dnm;
GRIP     = 100.0*(1.0 - exp(-SF*pow(SNAPG, HNMJ)))/dnm;
PHI_P    = PHI_SPIN*PHI_NMJ;

// =========================================================================
// D. GABA-A POTENTIATION  (primidone + phenobarbital + PEMA + ethanol +
//    octanol + topiramate + benzodiazepine), with counter-adaptation
// =========================================================================
EFFECT_PRM_GABA = EMAX_PRM_GABA*pow(CE_PRM,GAMMA_PRM_GABA)/(pow(EC50_PRM_GABA,GAMMA_PRM_GABA) + pow(CE_PRM,GAMMA_PRM_GABA));
EFFECT_PB_GABA  = EMAX_PB_GABA *pow(CE_PB,GAMMA_PB_GABA)  /(pow(EC50_PB_GABA,GAMMA_PB_GABA)   + pow(CE_PB,GAMMA_PB_GABA));
EFFECT_PEM_GABA = EMAX_PEM_GABA*pow(CU_PEM,GAMMA_PEM_GABA)/(pow(EC50_PEM_GABA,GAMMA_PEM_GABA) + pow(CU_PEM,GAMMA_PEM_GABA));
EFFECT_ETH_GABA = EMAX_ETH_GABA*pow(CE_ETH,GAMMA_ETH_GABA)/(pow(EC50_ETH_GABA,GAMMA_ETH_GABA) + pow(CE_ETH,GAMMA_ETH_GABA));
EFFECT_OCT_GABA = EMAX_OCT_GABA*pow(CB_OCT,GAMMA_OCT_GABA)/(pow(EC50_OCT_GABA,GAMMA_OCT_GABA) + pow(CB_OCT,GAMMA_OCT_GABA));
EFFECT_TOP_GABA = EMAX_TOP_GABA*pow(CB_TOP,GAMMA_TOP_GABA)/(pow(EC50_TOP_GABA,GAMMA_TOP_GABA) + pow(CB_TOP,GAMMA_TOP_GABA));
P_BZD = EMAX_BZD*BZDLEV;
P_RAW = EFFECT_PRM_GABA + EFFECT_PB_GABA + EFFECT_PEM_GABA + EFFECT_ETH_GABA + EFFECT_OCT_GABA + EFFECT_TOP_GABA + P_BZD;
if(P_RAW > 0.92) P_RAW = 0.92;
P_EFF = P_RAW/(1.0 + ADAPTF + ADAPTS);
// [FIXED — defect 4] Written as REB = 1+KREBF*ADAPTF, the adaptation blunted
// the acute suppression (PHI_CBL exceeded 1 at peak ethanol) and, with a single
// symmetric tau, never outlived the drug at all: the rebound came out +0.1%.
// The rebound is the UNOPPOSED part of the adaptation, so subtract the drive.
P_ETHR = EFFECT_ETH_GABA;  // identical formula/inputs -- reuse, don't recompute
reb_f = ADAPTF - KAF*P_RAW;  if(reb_f < 0) reb_f = 0.0;
reb_s = ADAPTS - KAS*P_ETHR; if(reb_s < 0) reb_s = 0.0;
REB   = 1.0 + KREBF*reb_f + KREBS*reb_s;      // >1 = rebound above baseline

// =========================================================================
// E. THE TWO PARALLEL BRANCHES OF THE CENTRAL LOOP
//    This is where the T-type ceiling lives.  phi_olive can be driven to zero
//    and Phi_C still cannot fall below a_R.
// =========================================================================
EFFECT_TTB_OL = EMAX_TTB_OL*pow(CE_TTB,GAMMA_TTB_OL)/(pow(EC50_TTB_OL,GAMMA_TTB_OL) + pow(CE_TTB,GAMMA_TTB_OL));
EFFECT_ETH_OL = EMAX_ETH_OL*pow(CE_ETH,GAMMA_ETH_OL)/(pow(EC50_ETH_OL,GAMMA_ETH_OL) + pow(CE_ETH,GAMMA_ETH_OL));
EFFECT_OCT_OL = EMAX_OCT_OL*pow(CB_OCT,GAMMA_OCT_OL)/(pow(EC50_OCT_OL,GAMMA_OCT_OL) + pow(CB_OCT,GAMMA_OCT_OL));
// combined only at this point of use -- each drug's own EFFECT_<STEM>_OL
// stays independently driveable; same shared ceiling cap the original
// applied to their pre-combined sum
ETH_OL = EFFECT_ETH_OL + EFFECT_OCT_OL;
if(ETH_OL > 0.95) ETH_OL = 0.95;
EFFECT_PRP_CBL = EMAX_PRP_CBL*pow(CE_PRP,GAMMA_PRP_CBL)/(pow(EC50_PRP_CBL,GAMMA_PRP_CBL) + pow(CE_PRP,GAMMA_PRP_CBL));
PHI_OL = (1.0 - EFFECT_TTB_OL)*(1.0 - ETH_OL)*(1.0 - KG_OL*P_EFF)
                *(1.0 + KG_OL*(REB - 1.0))*(1.0 + HARM);
PHI_CBL = (1.0 + DNDIS)*REB*(1.0 - P_EFF)
                 *(1.0 - EFFECT_PRP_CBL);

// =========================================================================
// F. THE SERIES ELEMENT — Vim relay.  Surgery acts HERE, which is why it is
//    not subject to the parallel-branch ceiling.
// =========================================================================
LES_EFF = 0.0;
if(VLES > 0) LES_EFF = LESION*(pow(VLES,HLES)/(pow(V50L,HLES) + pow(VLES,HLES)));
DBSF = pow(FSTIM,HDBS)/(pow(F50D,HDBS) + pow(FSTIM,HDBS));
DBSV = VTA/(V50D + VTA);
DBSB = EMAX_DBS*DBSF*DBSV*DBSON;
eterm = (F50E - FSTIM)/F50E;  if(eterm < 0) eterm = 0.0;
ENTR = 1.0 + KENT*eterm*DBSON;        // low-frequency stim can ENTRAIN
PHI_TH = (1.0 - 0.92*LES_EFF)*(1.0 - DBSB)*ENTR*(1.0 + KRR*REROUTE);
if(PHI_TH < 0.02) PHI_TH = 0.02;

// =========================================================================
// G. CORTICAL LIMB, LOOP GAIN, BIFURCATION
// =========================================================================
EFFECT_TOP_CTX  = EMAX_TOP_CTX*pow(CB_TOP,GAMMA_TOP_CTX)/(pow(EC50_TOP_CTX,GAMMA_TOP_CTX) + pow(CB_TOP,GAMMA_TOP_CTX));
EFFECT_PRM_NACH = EMAX_PRM_NACH*pow(CE_PRM,GAMMA_PRM_NACH)/(pow(EC50_PRM_NACH,GAMMA_PRM_NACH) + pow(CE_PRM,GAMMA_PRM_NACH));
EFFECT_GBP      = EMAX_GBP*pow(CB_GBP,GAMMA_GBP)/(pow(EC50_GBP,GAMMA_GBP) + pow(CB_GBP,GAMMA_GBP));
EFFECT_ETH_CTX  = EMAX_ETH_CTX*pow(CE_ETH,GAMMA_ETH_CTX)/(pow(EC50_ETH_CTX,GAMMA_ETH_CTX) + pow(CE_ETH,GAMMA_ETH_CTX));
PHI_CTX = (1.0 - EFFECT_TOP_CTX) * (1.0 - EFFECT_PRM_NACH) * (1.0 - EFFECT_GBP) * (1.0 - EFFECT_ETH_CTX)
               * (1.0 + KFAT*FATIGUE);
PHI_C = (AO*PHI_OL + AR*PHI_CBL)*PHI_TH*PHI_CTX;
G0T   = G0*(1.0 + PROG);
GTOT  = G0T*(WC*PHI_C + WP*PHI_P)*(1.0 + KVIS*TASK_INT);
MU    = GTOT - 1.0;
GHD   = G0T*HDG*(WC_HD*PHI_C + WP_HD*PHI_P)*(1.0 - 0.55*NECKBTX);
MU_HD = GHD - 1.0;
GVX   = G0T*VXG*(0.90*PHI_C + 0.10*PHI_P);
MU_VX = GVX - 1.0;

// =========================================================================
// H. FREQUENCY (delay) AND MECHANICS — note that NO drug appears here
// =========================================================================
JREL  = MLOAD*LLOAD*LLOAD/J0;
TAU_L = TAU0*(1.0 + KAGE*(AGE-50.0)/10.0 + KPC*(1.0-PCINT)
                     + KVISD*TASK_INT + KLOADN*JREL);
FNEUR = 1.0/TAU_L;
JTOT  = J0 + MLOAD*LLOAD*LLOAD;
F0M   = (1.0/(2.0*M_PI))*sqrt(KST/JTOT);
// [FIXED — defect 7] Holding the damping COEFFICIENT fixed makes zeta fall as
// 1/sqrt(J), which sharpened the resonance under load and DOUBLED ET amplitude
// (1.91 -> 3.70 cm at +500 g).  Loading provokes co-contraction, raising B with
// K, so the damping RATIO is what is held constant.
ZETA  = ZETA0;
rN    = FNEUR/F0M;
HN_   = 1.0/sqrt(pow(1.0-rN*rN,2.0) + pow(2.0*ZETA*rN,2.0));
HP_   = 1.0/(2.0*ZETA);
// [FIXED — defect 2] NZ = 1+KCAT*(RAG-1) returned -0.73 on propranolol, i.e. a
// negative physiological tremor.  Only FNZ of resting physiological tremor is
// beta2-dependent, so blockade can at most remove that fraction.
NZ = (RAG > 1.0) ? (1.0 + KCAT*(RAG-1.0)) : (1.0 - FNZ*(1.0-RAG));
if(NZ < 0.30) NZ = 0.30;
A_PHYS = KNZ*NZ*HP_;
A_LC   = KAMP*R_UL*HN_;
A_UL   = sqrt(A_LC*A_LC + A_PHYS*A_PHYS);
F_OBS  = (A_LC*FNEUR + A_PHYS*F0M)/(A_LC + A_PHYS);
F0_HD  = (1.0/(2.0*M_PI))*sqrt(KST_HD/J0_HD);
Z_HD   = ZETA0;
rH_    = FNEUR/F0_HD;
H_HD   = 1.0/sqrt(pow(1.0-rH_*rH_,2.0) + pow(2.0*Z_HD*rH_,2.0));
A_HD   = KAMP_HD*R_HD*H_HD;
A_VX   = 6.0*R_VX;

// =========================================================================
// I. PK / PD DIFFERENTIAL EQUATIONS
// =========================================================================
dxdt_GUT_PRP = -KA_PRP*GUT_PRP;
dxdt_CENT_PRP =  KA_PRP*GUT_PRP - CL_PRP/V1_PRP*CENT_PRP
               - Q_PRP*(CENT_PRP/V1_PRP - PERI_PRP/V2_PRP);
dxdt_PERI_PRP =  Q_PRP*(CENT_PRP/V1_PRP - PERI_PRP/V2_PRP);
dxdt_CE_PRP =  KE0_PRP*(KPUU_PRP*CF_PRP - CE_PRP);
dxdt_GUT_ATN = -KA_ATN*GUT_ATN;
dxdt_CENT_ATN =  KA_ATN*GUT_ATN - CL_ATN/V1_ATN*CENT_ATN;
dxdt_GUT_NAD = -KA_NAD*GUT_NAD;
dxdt_CENT_NAD =  KA_NAD*GUT_NAD - CL_NAD/V1_NAD*CENT_NAD;

KE_PRM = CL_PRM/V1_PRM;
dxdt_GUT_PRM = -KA_PRM*GUT_PRM;
dxdt_CENT_PRM =  KA_PRM*GUT_PRM - KE_PRM*CENT_PRM;
dxdt_CENT_PB  =  FM_PB *KE_PRM*CENT_PRM*(MW_PB /MW_PRM) - CL_PB /V1_PB *CENT_PB;
dxdt_CENT_PEM =  FM_PEM*KE_PRM*CENT_PRM*(MW_PEM/MW_PRM) - CL_PEM/V1_PEM*CENT_PEM;
dxdt_CE_PRM =  KE0_PRM*(KP_PRM*CU_PRM - CE_PRM);
dxdt_CE_PB  =  KE0_PB *(KP_PB *CU_PB  - CE_PB);

dxdt_GUT_TOP = -KA_TOP*GUT_TOP;
dxdt_CENT_TOP =  KA_TOP*GUT_TOP - CL_TOP*(1.0 + KIND*IND)/V1_TOP*CENT_TOP;
dxdt_GUT_GBP = -KA_GBP*GUT_GBP;
dxdt_CENT_GBP =  KA_GBP*GUT_GBP - CL_GBP/V1_GBP*CENT_GBP;

dxdt_GUT_ETH = -KA_ETH*GUT_ETH;
dxdt_CENT_ETH =  KA_ETH*GUT_ETH - VMAX_ETH*C_ETH/(KM_ETH + C_ETH);  // saturable
dxdt_CE_ETH =  KE0_ETH*(KP_ETH*C_ETH - CE_ETH);
dxdt_GUT_OCT = -KA_OCT*GUT_OCT;
dxdt_CENT_OCT =  KA_OCT*GUT_OCT - CL_OCT/V1_OCT*CENT_OCT;

dxdt_GUT_TTB = -KA_TTB*GUT_TTB;
dxdt_CENT_TTB =  KA_TTB*GUT_TTB - CL_TTB/V1_TTB*CENT_TTB;
dxdt_CE_TTB =  KE0_TTB*(KPUU_TTB*CF_TTB - CE_TTB);

dxdt_A_BTXT = -KDEG_BTX*A_BTXT;
dxdt_A_BTXG = -KDEG_BTX*A_BTXG;
dxdt_SNAPT  =  KR_SNAP*(1.0 - SNAPT) - KCL_BTX*A_BTXT*SNAPT;
dxdt_SNAPG  =  KR_SNAP*(1.0 - SNAPG) - KCL_BTX*A_BTXG*SNAPG;

// =========================================================================
// J. ADAPTATION AND DISEASE
// =========================================================================
tgtF = KAF*P_RAW;
dxdt_ADAPTF = (tgtF - ADAPTF)/((tgtF > ADAPTF) ? TAUF_ON : TAUF_OFF);
dxdt_ADAPTS = (KAS*P_ETHR - ADAPTS)/TAUS;
dxdt_B2REG  = (KB2REG*OCCB2 - B2REG)/TAU_B2REG;
dxdt_IND    = (CU_PB/(CU_PB + 60.0) - IND)/TAU_IND;

SED_RAW = 100.0*(0.55*CE_PB/(EC50_SED_PB + CE_PB)
                      + 0.45*CE_PRM/(EC50_SED_PRM + CE_PRM)
                      + 0.60*CE_ETH/(EC50_SED_ETH + CE_ETH)
                      + 0.30*CB_TOP/(EC50_SED_TOP + CB_TOP)
                      + 0.35*CB_GBP/(EC50_SED_GBP + CB_GBP)
                      + 0.25*CE_PRP/(EC50_SED_PRP + CE_PRP));
if(SED_RAW > 100.0) SED_RAW = 100.0;
dxdt_SEDTOL = (KSEDTOL*SED_RAW/100.0 - SEDTOL)/TAU_SEDTOL;

// [FIXED — defect 5] R_UL^2 overflowed once the envelope diverged; cap it.
rsq = R_UL*R_UL;  if(rsq > 1e6) rsq = 1e6;
dxdt_PCINT   = -KDEGPC*(1.0 + KEXC*rsq);
dxdt_DNDIS   = (2.4*(1.0 - PCINT) - DNDIS)/8760.0;
supp  = ((LES_EFF > 0.2) || (DBSB > 0.2)) ? 1.0 : 0.0;
dxdt_REROUTE = (REROUTE_MAX*supp - REROUTE)/TAU_RR;
dxdt_PROG    = KPROG*(1.0 + KEXC*rsq);

// =========================================================================
// K. THE OSCILLATOR — slow-envelope (Hopf normal form) reduction
//    [FIXED — defect 1] Written as dr/dt = (1/TAUA) r (mu - beta r^2), the
//    linearised relaxation time near r* is TAUA/(2 mu), so it SHRINKS with
//    severity: at G0=12 the system was stiff enough that RK4 returned NaN.
//    Dividing by (1+|mu|) leaves the fixed point r* = sqrt(mu/beta) exactly
//    unchanged and caps the relaxation time at TAUA/2 for every severity.
//    et_verify.py confirms r* reproduces the peak amplitude of the full 2-D
//    oscillator to 0.06-2.0% over mu = 0.15 to 5.
// =========================================================================
iT = 1.0/TAUA;
dxdt_R_UL = iT*(R_UL*(MU    - BETA*R_UL*R_UL) + EPS_NOISE)/(1.0 + fabs(MU));
dxdt_R_HD = iT*(R_HD*(MU_HD - BETA*R_HD*R_HD) + EPS_NOISE)/(1.0 + fabs(MU_HD));
dxdt_R_VX = iT*(R_VX*(MU_VX - BETA*R_VX*R_VX) + EPS_NOISE)/(1.0 + fabs(MU_VX));

// =========================================================================
// L. ORGAN SYSTEMS
// =========================================================================
FEV_SS = FEV10*(1.0 - EFEV_A*ASTHMA*OCCB2 - EFEV_N*OCCB2);
dxdt_FEV1 = (FEV_SS - FEV1)/TAU_FEV;
HCSS = HCO30 - EHCO3*CB_TOP/(EC50_HCO3 + CB_TOP);
dxdt_HCO3 = (HCSS - HCO3)/TAU_HCO3;
BWSS = BW0*(1.0 - EBW*CB_TOP/(EC50_BW + CB_TOP));
dxdt_BW   = (BWSS - BW)/TAU_BW;
dxdt_BMD  = -KBMD*IND*BMD;
dxdt_ALT  = (ALT0 + EALT*ETHCHR - ALT)/TAU_ALT;
CGSS = ECOG*CB_TOP/(EC50_COG + CB_TOP);
dxdt_COG  = (CGSS - COG)/TAU_COG;

$TABLE
// =========================================================================
// M. CLINICAL ENDPOINTS
//    Rating scales are LOGARITHMIC in amplitude (Elble 2006):
//        rating = 2.0 + 2.0*log10(A_cm),  so 1 point = 3.16-fold amplitude.
//    Everything downstream of this line is why accelerometry and TETRAS
//    disagree in every published trial.
// =========================================================================
HR   = HR0*(1.0 - EHR*OCCB1);
SBP  = SBP0 - ESBP*OCCB1;
SED  = SED_RAW/(1.0 + SEDTOL);
ATAXd = 0.50*CE_PB/(EC50_ATAX_PB + CE_PB) + 0.60*CE_ETH/(EC50_ATAX_ETH + CE_ETH)
             + 0.35*DBSV*DBSON*(1.0 + BILAT)/2.0;
if(VLES > 0) ATAXd += AXL*(VLES/(V50A + VLES));
ATAX = 100.0*(ATAXd > 1.0 ? 1.0 : ATAXd);
INTOXd = CE_ETH*KPCTX_ETH/(EC50_INTOX_ETH + CE_ETH*KPCTX_ETH)
              + CB_OCT*KPCTX_OCT/(EC50_INTOX_OCT + CB_OCT*KPCTX_OCT);
INTOX = 100.0*(INTOXd > 1.0 ? 1.0 : INTOXd);

T_R  = TCLAMP(2.0 + 2.0*log10(A_UL      > 1e-6 ? A_UL      : 1e-6));
T_L  = TCLAMP(2.0 + 2.0*log10(A_UL*ASYM > 1e-6 ? A_UL*ASYM : 1e-6));
T_HD = TCLAMP(2.0 + 2.0*log10(A_HD/3.0  > 1e-6 ? A_HD/3.0  : 1e-6));
T_VX = (A_VX > 1e-6) ? TCLAMP(1.2*log10(A_VX) + 1.5) : 0.0;
UL_T = 3.0*T_R + 3.0*T_L;                               // 3 manoeuvres/arm
SPI  = TCLAMP(0.85*T_R + 0.15) + TCLAMP(0.85*T_L + 0.15);
HW   = TCLAMP(0.90*T_R);
LL   = 2.0*TCLAMP(0.55*T_R);
FACE = TCLAMP(0.5*T_HD);
TONG = TCLAMP(0.5*T_HD);
STAND= TCLAMP(0.60*T_R + 0.40*T_HD);
TETRAS_PS = UL_T + SPI + HW + LL + T_HD + FACE + TONG + T_VX + STAND;
GRIPLOSS  = (100.0 - GRIP)/100.0;
ADLf = 0.24*T_R + 0.35*GRIPLOSS;
if(ADLf > 1.0) ADLf = 1.0;
if(ADLf < 0.0) ADLf = 0.0;
TETRAS_ADL = 48.0*ADLf;
FTM    = 2.25*TETRAS_PS;
SPIRAL = TCLAMP(0.85*T_R + 0.15);
BAINF  = (2.5*T_R > 10.0) ? 10.0 : 2.5*T_R;
QUESTd = 0.52*ADLf + 0.13*SED/100.0 + 0.12*ATAX/100.0
              + 0.09*COG/100.0 + 0.08*GRIPLOSS + 0.06*INTOX/100.0;
QUEST  = 100.0*(QUESTd > 1.0 ? 1.0 : QUESTd);

$CAPTURE @annotated
GTOT   : Total loop gain
MU     : Bifurcation parameter (G-1)
PHI_C  : Central limb gain factor
PHI_P  : Peripheral limb gain factor
PHI_OL : Olivary branch factor
PHI_CBL: Cerebello-thalamic branch factor
PHI_TH : Vim relay factor (series element)
PHI_SPIN : Spindle gain factor
PHI_NMJ  : Neuromuscular transmission factor
A_UL   : cm : Upper-limb tremor amplitude
A_LC   : cm : Limit-cycle component
A_PHYS : cm : Physiological tremor floor
A_HD   : deg : Head tremor amplitude
A_VX   : pct : Voice F0 modulation
FNEUR  : Hz : Central oscillation frequency
F0M    : Hz : Limb mechanical resonance
F_OBS  : Hz : Observed dominant peak frequency
T_R    : Dominant upper-limb TETRAS item
TETRAS_PS  : TETRAS performance subscale (0-64)
TETRAS_ADL : TETRAS ADL subscale (0-48)
FTM    : Fahn-Tolosa-Marin (0-144)
SPIRAL : Archimedes spiral rating (0-4)
BAINF  : Bain-Findley scale (0-10)
QUEST  : QUEST quality of life (0-100, higher worse)
OCCB2  : beta2 occupancy (efficacy AND the asthma contraindication)
OCCB1  : beta1 occupancy
P_RAW  : Raw GABA-A potentiation
P_EFF  : Adaptation-corrected GABA-A potentiation
EFFECT_PRM_GABA  : GABA-A potentiation from primidone itself
EFFECT_PB_GABA   : GABA-A potentiation from phenobarbital
EFFECT_PEM_GABA  : GABA-A potentiation from PEMA (primidone metabolite)
EFFECT_ETH_GABA  : GABA-A potentiation from ethanol
EFFECT_OCT_GABA  : GABA-A potentiation from 1-octanol
EFFECT_TOP_GABA  : GABA-A potentiation from topiramate
REB    : Rebound gain factor
EFFECT_TTB_OL : T-type/Cav3.1 block fraction (CX-8998-class blocker)
EFFECT_ETH_OL : Ethanol's own olivary suppression term (pre-cap)
EFFECT_OCT_OL : 1-octanol's own olivary suppression term (pre-cap)
EFFECT_PRP_CBL  : Propranolol's cerebello-thalamic (central) effect
EFFECT_PRM_NACH : Primidone's cortical sodium-channel effect
EFFECT_GBP      : Gabapentin's cortical effect
EFFECT_ETH_CTX  : Ethanol's cortical effect
EFFECT_TOP_CTX  : Topiramate's cortical effect
EFFECT_PRP_B2 : Propranolol's own beta2 receptor-occupancy quotient
EFFECT_PRP_B1 : Propranolol's own beta1 receptor-occupancy quotient
EFFECT_ATN_B2 : Atenolol's own beta2 receptor-occupancy quotient
EFFECT_ATN_B1 : Atenolol's own beta1 receptor-occupancy quotient
EFFECT_NAD_B2 : Nadolol's own beta2 receptor-occupancy quotient
EFFECT_NAD_B1 : Nadolol's own beta1 receptor-occupancy quotient
LES_EFF: Effective Vim lesion fraction
DBSB   : DBS block fraction
GRIP   : pct : Grip strength
HR     : bpm : Heart rate
SBP    : mmHg : Systolic blood pressure
SED    : Sedation index (0-100)
ATAX   : Ataxia index (0-100)
INTOX  : Intoxication index (0-100)
C_PRP  : mg/L : Propranolol plasma
C_ATN  : mg/L : Atenolol plasma
C_NAD  : mg/L : Nadolol plasma
C_PRM  : mg/L : Primidone plasma
C_PB   : mg/L : Phenobarbital plasma
C_PEM  : mg/L : PEMA plasma
C_ETH  : g/L : Blood ethanol
C_TOP  : mg/L : Topiramate plasma
C_GBP  : mg/L : Gabapentin plasma
C_OCT  : mg/L : 1-octanol plasma
C_TTB  : mg/L : T-type blocker plasma
CB_TOP : uM : Topiramate brain-equivalent (algebraic partition)
CB_GBP : uM : Gabapentin brain-equivalent (algebraic partition)
CB_OCT : uM : 1-octanol brain-equivalent (algebraic partition)
'

mod <- mcode("essential_tremor", code, atol = 1e-8, rtol = 1e-8)

# =============================================================================
#  HELPERS
# =============================================================================
W <- 24 * 7                      # one week in hours

# oral regimen -> mrgsolve event object
reg <- function(cmt, amt, ii = 24, dur_wk = 24, start = 0) {
  ev(amt = amt, cmt = cmt, ii = ii, addl = ceiling(dur_wk * W / ii) - 1,
     time = start)
}

# run and return the mean of every output over the final dosing interval.
# Reading a single trough understates a once-daily drug: propranolol 160 mg
# reads -31% at trough but -37% on the interval mean, and trial assessments
# are not troughs.
final_mean <- function(m, e = NULL, wk = 24, last = 24, ...) {
  out <- mrgsim(m, events = e, end = wk * W, delta = 1, ...) %>% as.data.frame()
  tail(out, last + 1) %>% summarise(across(where(is.numeric), mean))
}

# =============================================================================
#  SCENARIOS
# =============================================================================
cat("\n=== ET QSP: 22 scenarios ===================================\n")

## --- S0  patient phenotypes -------------------------------------------------
## G0 is the untreated loop gain.  It is the only thing that distinguishes a
## patient who will be a spectacular responder from one who will not.
PH <- list(mild = 1.15, moderate = 1.60, severe = 6.00, very_severe = 12.0)

s0 <- lapply(names(PH), function(nm) {
  final_mean(param(mod, G0 = PH[[nm]])) %>% mutate(pheno = nm)
}) %>% bind_rows()
print(s0 %>% select(pheno, GTOT, MU, A_UL, F_OBS, TETRAS_PS))
BASE <- s0[s0$pheno == "moderate", ]

## --- S1-S5  beta-blockade, and why beta1-selectivity fails ------------------
## Propranolol acts on beta2 receptors on the muscle SPINDLE.  Atenolol is
## beta1-selective (Ki_b2 1000 nM) AND barely enters the brain; nadolol is
## non-selective and peripherally restricted.  The ordering of the three is
## therefore a prediction, not an input.
s1 <- bind_rows(
  final_mean(param(mod, DPRP =  60), reg("GUT_PRP",  60)) %>% mutate(rx="propranolol 60"),
  final_mean(param(mod, DPRP = 120), reg("GUT_PRP", 120)) %>% mutate(rx="propranolol 120"),
  final_mean(param(mod, DPRP = 160), reg("GUT_PRP", 160)) %>% mutate(rx="propranolol 160"),
  final_mean(param(mod, DPRP = 240), reg("GUT_PRP", 240)) %>% mutate(rx="propranolol 240"),
  final_mean(param(mod, DPRP = 320), reg("GUT_PRP", 320)) %>% mutate(rx="propranolol 320"),
  final_mean(mod, reg("GUT_ATN", 100)) %>% mutate(rx="atenolol 100"),
  final_mean(mod, reg("GUT_NAD", 120)) %>% mutate(rx="nadolol 120")
) %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL,
             dTETRAS = TETRAS_PS - BASE$TETRAS_PS)
print(s1 %>% select(rx, OCCB2, A_UL, dA_pct, TETRAS_PS, dTETRAS, HR))

## --- S6  THE LOG-SCALE ARTEFACT ---------------------------------------------
## Nothing here is a new mechanism.  The point is that a -37% amplitude change
## is -0.41 points on an item, because the scale is a logarithm.  Reporting
## "% change in TETRAS" is a category error.
z <- s1[s1$rx == "propranolol 160", ]
cat(sprintf("\n[log artefact] amplitude %+.1f%%  ->  TETRAS %+.2f pts (%+.1f%%), item %+.2f pts\n",
            z$dA_pct, z$dTETRAS, 100*z$dTETRAS/BASE$TETRAS_PS, z$T_R - BASE$T_R))

## --- S7  SEVERITY DEPENDENCE — the sqrt law ---------------------------------
## Identical drug, identical receptor occupancy, four patients.  r* = sqrt(mu)
## means the mild patient is pushed BELOW the bifurcation and the severe one is
## not.  This is the responder/non-responder split, derived.
s7 <- lapply(names(PH), function(nm) {
  b <- final_mean(param(mod, G0 = PH[[nm]]))
  t <- final_mean(param(mod, G0 = PH[[nm]], DPRP = 160), reg("GUT_PRP", 160))
  data.frame(pheno = nm, A0 = b$A_UL, A_on = t$A_UL,
             dA_pct = 100*(t$A_UL - b$A_UL)/b$A_UL,
             dTETRAS = t$TETRAS_PS - b$TETRAS_PS)
}) %>% bind_rows()
print(s7)

## --- S8-S10  primidone: the parent molecule is the active moiety ------------
## Primidone t1/2 ~10 h, phenobarbital t1/2 ~100 h.  If phenobarbital were the
## active moiety the effect could not appear on day 1.  The model is told only
## the two potencies and the two half-lives; the time course is an output.
s8 <- bind_rows(
  final_mean(mod, reg("GUT_PRM", 250/3, ii = 8)) %>% mutate(rx="primidone 250"),
  final_mean(mod, reg("GUT_PRM", 500/3, ii = 8)) %>% mutate(rx="primidone 500"),
  final_mean(mod, reg("GUT_PRM", 750/3, ii = 8)) %>% mutate(rx="primidone 750")
) %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL,
             pct_parent = 100*EFFECT_PRM_GABA/(EFFECT_PRM_GABA + EFFECT_PB_GABA))
print(s8 %>% select(rx, C_PRM, C_PB, EFFECT_PRM_GABA, EFFECT_PB_GABA, pct_parent, A_UL, dA_pct, SED))

s9 <- mrgsim(mod, events = reg("GUT_PRM", 250/3, ii = 8), end = 24*W, delta = 2) %>%
  as.data.frame() %>% filter(time %in% c(2,6,12,24,72,336,4032))
cat("\n[primidone onset] the effect is present before phenobarbital exists:\n")
print(s9 %>% select(time, C_PRM, C_PB, EFFECT_PRM_GABA, EFFECT_PB_GABA, A_UL))

## the phenobarbital level that would be needed to match primidone 250 mg/d
p <- as.list(param(mod))
tgt <- s8$P_RAW[s8$rx == "primidone 250"]
need <- p$EC50_PB_GABA * tgt / max(p$EMAX_PB_GABA - tgt, 1e-9)
cat(sprintf("[phenobarbital equivalence] needs %.0f uM brain = %.1f mg/L plasma; therapeutic range is 10-40 mg/L\n",
            need, need/p$KP_PB*p$MW_PB/1e3))

## --- S11  COMBINATION IS SUPRA-ADDITIVE -------------------------------------
## The author expected combination therapy to be SUB-additive (two drugs each
## cutting gain, amplitude going as sqrt).  et_verify.py refuted that: sqrt is
## CONCAVE, so pushing mu toward zero moves onto the steepest part of the map
## and the second drug buys MORE than it would alone.  Reported as the model
## found it, not as it was expected.
s11 <- final_mean(param(mod, DPRP = 160),
                  c(reg("GUT_PRP", 160), reg("GUT_PRM", 250/3, ii = 8)))
cat(sprintf("\n[combination] propranolol %+.1f%% | primidone %+.1f%% | sum %+.1f%% | observed %+.1f%%\n",
            s1$dA_pct[s1$rx=="propranolol 160"], s8$dA_pct[s8$rx=="primidone 250"],
            s1$dA_pct[s1$rx=="propranolol 160"] + s8$dA_pct[s8$rx=="primidone 250"],
            100*(s11$A_UL - BASE$A_UL)/BASE$A_UL))

## --- S12-S13  ethanol: suppression, then rebound ---------------------------
## The rebound is not an extra assumption.  ADAPTF has tau = 8 h and ethanol
## is gone in ~4 h, so the counter-adaptation necessarily outlives the drug.
s12 <- mrgsim(mod, events = ev(amt = 28, cmt = "GUT_ETH"), end = 40, delta = 0.25) %>%
  as.data.frame()
cat(sprintf("\n[ethanol 2 drinks] nadir %.3f cm (%+.1f%%) at %.2f h; rebound peak %+.1f%% at %.1f h\n",
            min(s12$A_UL), 100*(min(s12$A_UL)-BASE$A_UL)/BASE$A_UL,
            s12$time[which.min(s12$A_UL)],
            100*(max(s12$A_UL[s12$time>6])-BASE$A_UL)/BASE$A_UL,
            s12$time[s12$time>6][which.max(s12$A_UL[s12$time>6])]))

s13 <- lapply(1:4, function(nd) {
  o <- mrgsim(mod, events = ev(amt = 14*nd, cmt = "GUT_ETH"), end = 48, delta = 0.25) %>%
    as.data.frame()
  data.frame(drinks = nd,
             nadir_pct = 100*(min(o$A_UL)-BASE$A_UL)/BASE$A_UL,
             rebound_pct = 100*(max(o$A_UL[o$time>6])-BASE$A_UL)/BASE$A_UL,
             peak_INTOX = max(o$INTOX))
}) %>% bind_rows()
print(s13)

## --- S14  daily ethanol: the self-medication trap --------------------------
## ADAPTS (tau 30 d) both blunts the drug and shifts baseline gain UP, so the
## sober-morning tremor gets worse while the evening relief stays.  That is an
## escalation gradient, generated rather than asserted.
s14 <- mrgsim(mod, events = ev(amt = 42, cmt = "GUT_ETH", ii = 24, addl = 89,
                               time = 19), end = 91*24, delta = 1) %>%
  as.data.frame() %>% filter(time %% 24 == 6)      # sober each morning
cat("\n[daily ethanol, sober-morning tremor]\n")
print(s14 %>% filter(time/24 %in% c(1,7,14,30,60,89)) %>%
        mutate(day = time %/% 24) %>% select(day, ADAPTS, GTOT, A_UL))

## --- S15  1-octanol: tremor benefit per unit intoxication -----------------
s15 <- bind_rows(
  {o <- mrgsim(mod, events=ev(amt=28,  cmt="GUT_ETH"), end=24, delta=0.25) %>% as.data.frame()
   data.frame(agent="ethanol 2 drinks", dA=100*(min(o$A_UL)-BASE$A_UL)/BASE$A_UL, INTOX=max(o$INTOX))},
  {o <- mrgsim(mod, events=ev(amt=560, cmt="GUT_OCT"), end=24, delta=0.25) %>% as.data.frame()
   data.frame(agent="1-octanol 8 mg/kg", dA=100*(min(o$A_UL)-BASE$A_UL)/BASE$A_UL, INTOX=max(o$INTOX))},
  {o <- mrgsim(mod, events=ev(amt=1120,cmt="GUT_OCT"), end=24, delta=0.25) %>% as.data.frame()
   data.frame(agent="1-octanol 16 mg/kg", dA=100*(min(o$A_UL)-BASE$A_UL)/BASE$A_UL, INTOX=max(o$INTOX))}
) %>% mutate(benefit_per_intox = abs(dA)/pmax(INTOX, 0.01))
print(s15)

## --- S16  T-TYPE BLOCKER: A CEILING SET BY TOPOLOGY -----------------------
## The headline result.  Cav3.1 sits on ONE of two parallel branches, so even a
## perfect blocker leaves a_R of the central drive standing.  In the harmaline
## rat the oscillation IS the olive (a_O = 1) and the same drug is curative.
## The species discordance that killed CX-8998's rationale is ONE parameter.
s16 <- bind_rows(
  final_mean(mod, reg("GUT_TTB", 100), wk = 12) %>% mutate(setting="human ET, 100 mg/d"),
  final_mean(param(mod, EMAX_TTB_OL = 1.0, EC50_TTB_OL = 1.0), reg("GUT_TTB", 100), wk = 12) %>%
    mutate(setting="human ET, perfect Cav3 block")
) %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL)
print(s16 %>% select(setting, PHI_OL, PHI_C, GTOT, A_UL, dA_pct))

rat  <- param(mod, AO = 1.0, AR = 0.0, WC = 0.95, WP = 0.05, G0 = 2.2)
rb   <- final_mean(rat, wk = 12)
rt   <- final_mean(rat, reg("GUT_TTB", 100), wk = 12)
cat(sprintf("[harmaline rat, a_O=1] same drug: %+.1f%% (human ET: %+.1f%%)\n",
            100*(rt$A_UL - rb$A_UL)/rb$A_UL, s16$dA_pct[1]))

cat("[ceiling sweep] the best a perfect Cav3.1 blocker can ever do vs a_O:\n")
for (ao in c(0.20, 0.35, 0.50, 0.65, 0.80, 1.00)) {
  pb <- param(mod, AO = ao, AR = 1 - ao)
  b  <- final_mean(pb, wk = 8)
  t  <- final_mean(param(pb, EMAX_TTB_OL = 1.0, EC50_TTB_OL = 1.0), reg("GUT_TTB", 100), wk = 8)
  cat(sprintf("   a_O = %.2f -> ceiling %+6.1f %%\n", ao, 100*(t$A_UL-b$A_UL)/b$A_UL))
}

## --- S17  MRgFUS: one lesion volume, two outcomes -------------------------
## Benefit saturates (V50 = 45 mm3) while ataxia does not (V50 = 260 mm3).
## The therapeutic volume window is therefore derived, not chosen.
s17 <- lapply(c(20,40,60,90,120,180,250,400), function(v) {
  final_mean(param(mod, LESION = 1, VLES = v), wk = 8) %>% mutate(V = v)
}) %>% bind_rows() %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL)
print(s17 %>% select(V, LES_EFF, GTOT, A_UL, dA_pct, ATAX))

## --- S18  DBS frequency is a switch, not a dial --------------------------
## Hill 4 on frequency with f50 = 80 Hz reproduces the clinical >100 Hz rule,
## and the entrainment term reproduces the fact that LOW-frequency Vim
## stimulation can make tremor worse.
s18 <- lapply(c(10,30,50,60,80,100,130,185,250), function(f) {
  final_mean(param(mod, DBSON = 1, FSTIM = f), wk = 8) %>% mutate(f = f)
}) %>% bind_rows() %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL)
print(s18 %>% select(f, DBSB, PHI_TH, GTOT, A_UL, dA_pct))

## --- S19  habituation: who relapses is a capacity, not a duration --------
s19 <- lapply(c(0.45, 0.75, 0.95, 1.00), function(rm) {
  lapply(c(1, 3, 5), function(yr) {
    final_mean(param(mod, LESION = 1, VLES = 120, REROUTE_MAX = rm),
               wk = yr*52) %>% mutate(REROUTE_MAX = rm, yr = yr)
  }) %>% bind_rows()
}) %>% bind_rows() %>% mutate(dA_pct = 100*(A_UL - BASE$A_UL)/BASE$A_UL)
print(s19 %>% select(REROUTE_MAX, yr, REROUTE, PHI_TH, A_UL, dA_pct))

## --- S20  BOTULINUM: the window is PRECISION, not dose -------------------
## Two SNAP-25 pools, one in the tremorogenic wrist muscles and one in the
## finger flexors that carry grip.  f_spill is the only thing that separates
## them.  Halving the dose scales BOTH down and opens nothing.
s20 <- lapply(list(c(100,0.15), c(100,0.45), c(50,0.45), c(150,0.15)), function(x) {
  d <- x[1]; fs <- x[2]
  o <- mrgsim(param(mod, FSPILL = fs),
              events = c(ev(amt = d*(1-fs), cmt = "A_BTXT"),
                         ev(amt = d*fs,     cmt = "A_BTXG")),
              end = 6*W, delta = 24) %>% as.data.frame()
  k <- which.min(o$A_UL)
  data.frame(dose = d, f_spill = fs, SNAPT = o$SNAPT[k], SNAPG = o$SNAPG[k],
             dA_pct = 100*(o$A_UL[k]-BASE$A_UL)/BASE$A_UL,
             grip = o$GRIP[k], QUEST = o$QUEST[k])
}) %>% bind_rows()
print(s20)

## --- S21  WEIGHT LOADING: the model reproduces the diagnostic test -------
## ET frequency is central (delay), EPT frequency IS the mechanical resonance.
## Adding 500 g moves f0 as 1/sqrt(J) and barely touches the central delay, so
## the standard clinical discriminator falls out of the equations.
s21 <- bind_rows(
  final_mean(mod, wk = 8) %>% mutate(cond="ET unloaded"),
  final_mean(param(mod, MLOAD = 0.5), wk = 8) %>% mutate(cond="ET +500 g"),
  final_mean(param(mod, G0 = 0.42, THYRO = 9.0), wk = 8) %>% mutate(cond="EPT unloaded"),
  final_mean(param(mod, G0 = 0.42, THYRO = 9.0, MLOAD = 0.5), wk = 8) %>% mutate(cond="EPT +500 g")
)
print(s21 %>% select(cond, MU, F0M, FNEUR, F_OBS, A_LC, A_PHYS, A_UL))

## --- S22  head tremor emerges when the NECK effector crosses G = 1 -------
s22 <- lapply(list(c(0.55,1.6), c(0.85,1.6), c(0.85,2.4), c(0.85,6.0)), function(x) {
  b <- final_mean(param(mod, HDG = x[1], G0 = x[2]), wk = 8)
  t <- final_mean(param(mod, HDG = x[1], G0 = x[2], DPRP = 160),
                  reg("GUT_PRP", 160), wk = 8)
  data.frame(HDG = x[1], G0 = x[2], MU_HD = b$MU_HD, A_HD = b$A_HD,
             arm_dA = 100*(t$A_UL-b$A_UL)/b$A_UL,
             head_dA = ifelse(b$A_HD > 1e-3, 100*(t$A_HD-b$A_HD)/b$A_HD, NA))
}) %>% bind_rows()
print(s22)

## --- S23  the contraindication IS the efficacy --------------------------
## OCCB2 appears in exactly two places: the spindle gain and the airway.  So
## the tremor benefit and the FEV1 fall in an asthmatic are the same number.
s23 <- bind_rows(
  final_mean(param(mod, ASTHMA = 1, DPRP = 160), reg("GUT_PRP", 160), wk = 8) %>% mutate(rx="propranolol 160"),
  final_mean(param(mod, ASTHMA = 1), reg("GUT_NAD", 120), wk = 8) %>% mutate(rx="nadolol 120"),
  final_mean(param(mod, ASTHMA = 1), reg("GUT_ATN", 100), wk = 8) %>% mutate(rx="atenolol 100")
) %>% mutate(dA_pct = 100*(A_UL-BASE$A_UL)/BASE$A_UL,
             dFEV1_pct = 100*(FEV1 - 3.20)/3.20)
print(s23 %>% select(rx, OCCB2, dA_pct, FEV1, dFEV1_pct))

## --- S24  progression: sqrt law -> fastest worsening EARLY --------------
## dA/dt = (dG/dt)/(2 sqrt(mu)) is largest when mu is smallest.  A constant
## mechanistic progression rate therefore produces a DECELERATING amplitude,
## which is the clinical picture patients describe.
s24 <- lapply(c(1,2,5,10,15,20), function(yr) {
  final_mean(mod, wk = yr*52) %>% mutate(yr = yr)
}) %>% bind_rows()
print(s24 %>% select(yr, PROG, GTOT, MU, A_UL, TETRAS_PS, F_OBS))

## --- S25  is early suppression disease-modifying? (explicit hypothesis) --
## KEXC couples oscillation amplitude back onto cerebellar injury.  Default 0.
## Switched on, it makes early tremor control disease-modifying — an open
## question, flagged as a hypothesis and NOT as a result.
s25 <- bind_rows(
  final_mean(param(mod, KEXC = 0.00), wk = 10*52) %>% mutate(arm="KEXC 0, untreated"),
  final_mean(param(mod, KEXC = 0.30), wk = 10*52) %>% mutate(arm="KEXC 0.3, untreated"),
  final_mean(param(mod, KEXC = 0.30, DPRP = 160), reg("GUT_PRP", 160, dur_wk = 10*52),
             wk = 10*52) %>% mutate(arm="KEXC 0.3, propranolol 10 y")
)
print(s25 %>% select(arm, PROG, GTOT, A_UL, TETRAS_PS))

cat("\n=== end ====================================================\n")
