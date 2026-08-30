# =============================================================================
# Idiopathic Multicentric Castleman Disease (iMCD) — mrgsolve QSP Model
# REFACTORED: pluggable PK / named Hill interface (FORK_WORKFLOW_GUIDE.md Part 2)
#
# Derived from imcd_mrgsolve_model.R. Original untouched (upstream-tracked);
# every compound's PK is reorganized into its own clearly-delimited block,
# renamed to the fork's `<ROLE>_<STEM>` convention, and exposes exactly one
# `C_<STEM>` concentration plus (where the original had one) one `EFFECT_<STEM>`
# disease-effect term. See imcd_refactor_notes.md for the full account.
#
# 19 ODE compartments, 7 treatment scenarios (unchanged from the original)
# Calibration: CONCERT trial (van Rhee 2014 Lancet Oncol) for siltuximab
#              ACTEMRA studies (Nishimoto 2005 Blood) for tocilizumab
#              TAFRO subtype (Iwaki 2016, Fajgenbaum 2019 JCI) for sirolimus
#
# Author: QSP CCR Library · CDCN 2017 diagnostic criteria
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

code <- '
$PROB
# iMCD QSP Model (Idiopathic Multicentric Castleman Disease) -- refactored PK
- IL-6 cytokine storm + lymphadenopathy + acute phase
- Drugs : siltuximab (SILT), tocilizumab (TOCZ), sirolimus (SIRO),
          rituximab (RTX), anakinra (ANA), ruxolitinib (RUX),
          CHOP-doxorubicin (DOXO), CHOP-cyclophosphamide (CYC), prednisone (PRED)
- Endpoints : CRP, IL-6 (total/free), Hb, IgG, LN tumor burden,
              VEGF, platelet, CDCN response, TAFRO severity, OS hazard

$PARAM @annotated
// ---------------- siltuximab (SILT, anti-IL-6 mAb, FDA 2014) -- Archetype 2 --
CL_SILT   :  0.232  : siltuximab CL (L/d) van Rhee 2014 PopPK
V1_SILT   :  4.47   : Vc (L)
V2_SILT   :  2.49   : Vp (L)
Q_SILT    :  0.527  : Q (L/d)
EC50_SILT :  1e-3   : IL-6 binding Kd (nM approx 1 pM) Kurzrock 2013 (renamed from Kd_SILT)
EMAX_SILT :  1      : Hill ceiling (math-implied -- original had no explicit Emax)
GAMMA_SILT:  1      : Hill coefficient (original had no explicit Hill exponent)
MW_SILT   :  145000 : Da
DOSE_SILT :  11     : mg/kg q3w label dose
// ---------------- tocilizumab (TOCZ, anti-IL-6R mAb) -- Archetype 2 --------
CL_TOCZ   :  0.20   : (L/d) Frey 2010
V1_TOCZ   :  4.1    :
V2_TOCZ   :  2.5    :
Q_TOCZ    :  0.50   :
EC50_TOCZ :  2.5     : nM (IL-6R) (renamed from Kd_TOCZ)
EMAX_TOCZ :  1      : Hill ceiling (math-implied)
GAMMA_TOCZ:  1      : Hill coefficient (original had no explicit exponent)
DOSE_TOCZ :  8      : mg/kg q2w
// ---------------- sirolimus (SIRO, mTORC1) -- Archetype 3 (depot+cent+peri) --
CL_SIRO   :  8.2    : L/h  -> 197 L/d MacDonald 2000
V1_SIRO   :  12     : L
V2_SIRO   :  100    : L
Q_SIRO    :  10     : L/h
F_SIRO    :  0.14   : oral bioavail
KA_SIRO   :  2.0    : 1/h (renamed from Ka_SIRO)
DOSE_SIRO :  2      : mg PO QD (target trough 6-14 ng/mL)
EC50_SIRO :  8      : ng/mL
EMAX_SIRO :  1      : Hill ceiling (math-implied)
GAMMA_SIRO:  1      : Hill coefficient (original had no explicit exponent)
// ---------------- rituximab (RTX, anti-CD20) -- Archetype 2 ---------------
CL_RTX    :  0.32   : L/d
V1_RTX    :  3.3    :
V2_RTX    :  4.5    :
Q_RTX     :  0.42   :
DOSE_RTX  :  375    : mg/m2 weekly x4 induction
KILL_RTX  :  0.030  : 1/(d * (mg/L)) CD20 lysis rate constant (renamed from RTX_kill; linear, not Hill-saturating -- bespoke effect, see notes)
// ---------------- anakinra (ANA, IL-1Ralpha) -- Archetype 3 minus peripheral --
CL_ANA    :  9.0    : L/h
V1_ANA    :  20     : L (renamed from V_ANA)
F_ANA     :  0.95   : SC
KA_ANA    :  0.40   : 1/h (renamed from Ka_ANA)
DOSE_ANA  :  100    : mg SC QD
// ---------------- ruxolitinib (RUX, JAK1/2) -- Archetype 3 minus peripheral --
CL_RUX    :  18     : L/h
V1_RUX    :  72     : L (renamed from V_RUX)
F_RUX     :  0.95   :
KA_RUX    :  3.0    : 1/h [declared but UNUSED in the original -- its own ODE reads Ka_SIRO instead; preserved unfixed]
DOSE_RUX  :  20     : mg PO BID
EC50_RUX  :  300    : nM JAK1/2 (renamed from IC50_RUX)
EMAX_RUX  :  1      : Hill ceiling (math-implied)
GAMMA_RUX :  1      : Hill coefficient (original had no explicit exponent)
// ---------------- CHOP doxorubicin (DOXO) -- Archetype 1 ------------------
CL_DOXO   :  45     : L/h
V1_DOXO   :  500    : L (renamed from V_DOXO)
DOSE_DOXO :  50     : mg/m2 q21d
EMAX_DOXO :  0.05   : max fractional cytotoxic-kill contribution (renamed from literal 0.05)
EC50_DOXO :  0.1    : mg/L (renamed from literal 0.1)
GAMMA_DOXO:  1       : Hill coefficient (original had no explicit exponent)
// ---------------- CHOP cyclophosphamide (CYC) -- Archetype 1 --------------
CL_CYC    :  6      : L/h
V1_CYC    :  35     : L (renamed from V_CYC)
DOSE_CYC  :  750    : mg/m2 q21d
EMAX_CYC  :  0.04   : max fractional cytotoxic-kill contribution (renamed from literal 0.04)
EC50_CYC  :  5      : mg/L (renamed from literal 5)
GAMMA_CYC :  1       : Hill coefficient (original had no explicit exponent)
// ---------------- Prednisone (PRED) -- Archetype 3 minus peripheral -------
CL_PRED   :  5.6    : L/h
V1_PRED   :  35     : L (renamed from V_PRED)
KA_PRED   :  2.5    : 1/h (renamed from Ka_PRED)
F_PRED    :  0.85   :
EC50_PRED :  40     : ng/mL
EMAX_PRED :  1      : Hill ceiling (math-implied)
GAMMA_PRED:  1       : Hill coefficient (original had no explicit exponent)
DOSE_PRED :  60     : mg PO QD (approx 1 mg/kg)
//
// ---------------- IL-6 axis & turnover --------------------------------
kIL6_base : 0.20    : baseline IL-6 production (pg/mL/d) -- healthy
kIL6_iMCD : 2.5     : iMCD over-production multiplier
kIL6_deg  : 8       : IL-6 elimination 1/d (half-life ~2 h, but in serum apparent)
sIL6R_b   : 50      : ng/mL baseline soluble IL-6R
//
// ---------------- Lymph node / plasmablast ----------------------------
kgrow_LN  : 0.012   : 1/d growth rate of LN size (cm composite)
LN_max    : 25      : cm composite
kshrink_LN: 0.005   : 1/d natural turnover
//
kprol_PB  : 0.20    : plasmablast proliferation 1/d (under IL-6)
kdeath_PB : 0.10    : 1/d basal apoptosis
PB_base   : 0.05    : baseline plasmablast fraction (% of LN)
//
kprol_Bmem: 0.08    : memory B turnover
kdeath_Bmem: 0.04   :
//
// ---------------- Acute phase ----------------------------------------
kCRP_in   : 12      : mg/L/d max CRP synthesis under IL-6
kCRP_out  : 0.55    : 1/d (CRP half-life ~ 19 h)
EC50_CRP  : 5       : pg/mL IL-6
//
kHb_in    : 0.16    : g/dL/d homeostasis
kHb_out   : 0.012   : 1/d (RBC lifespan ~120 d)
Hb_max    : 14.5    : g/dL
Hepcidin_EC50 : 30  : pg/mL IL-6 driving hepcidin
//
kIgG_in   : 0.6     : g/dL/d
kIgG_out  : 0.035   : 1/d (IgG t1/2 ~21 d)
IgG_max   : 5.0     : g/dL polyclonal max under disease
//
// ---------------- VEGF / TAFRO --------------------------------------
kVEGF_in  : 50      : pg/mL/d
kVEGF_out : 1.5     : 1/d
EC50_VEGF : 20      : pg/mL IL-6
//
kAnasarca_in : 0.20 :  units/d
kAnasarca_out: 0.10 :  1/d
EC50_VEGFan  : 200  : pg/mL VEGF
//
// ---------------- Platelet ------------------------------------------
kPlt_in   : 25      : x10^9/L per d
kPlt_out  : 0.10    : 1/d
Plt_max   : 400     :
Plt_TAFRO_kill : 0.15 : platelet consumption rate driven by IL-6/VEGF
//
// ---------------- mTORC1 / TAFRO disease driver ---------------------
kmTOR_in  : 0.10    : arbitrary signal/d
kmTOR_out : 0.05    : 1/d
//
// ---------------- Hazard / safety ----------------------------------
h0_OS     : 5e-5    : 1/d baseline 5-yr OS hazard ~ 0.09
beta_CRP  : 0.0008  : CRP hazard coefficient

$CMT @annotated
CENT_SILT : siltuximab central (mg)
PERI_SILT : siltuximab peripheral (mg)
CENT_TOCZ : tocilizumab central (mg)
PERI_TOCZ : tocilizumab peripheral (mg)
GUT_SIRO  : sirolimus gut (mg)
CENT_SIRO : sirolimus central (mg)
PERI_SIRO : sirolimus peripheral (mg)
CENT_RTX  : rituximab central (mg)
PERI_RTX  : rituximab peripheral (mg)
GUT_ANA   : anakinra SC depot (mg)
CENT_ANA  : anakinra central (mg)
GUT_RUX   : ruxolitinib gut (mg)
CENT_RUX  : ruxolitinib central (mg)
CENT_DOXO : doxorubicin central (mg)
CENT_CYC  : cyclophosphamide central (mg)
GUT_PRED  : prednisone gut (mg)
CENT_PRED : prednisone central (mg)
IL6_T     : IL-6 serum (pg/mL)  -- TOTAL (free + bound)
IL6_F     : IL-6 free  (pg/mL)
LN        : lymph node composite size (cm)
PB        : plasmablast burden  (fraction LN)
Bmem      : memory B-cell pool (% baseline)
CRP       : C-reactive protein (mg/L)
Hb        : Hemoglobin (g/dL)
IgG       : Polyclonal IgG (g/dL)
VEGF      : VEGF-A (pg/mL)
Anasarca  : Anasarca / 3rd-space fluid score (0-10)
Plt       : Platelet x10^9/L
mTOR      : mTORC1 activity (0-1)
HAZ       : Cumulative OS hazard (death by AE/disease)

$MAIN
// initial conditions
IL6_T_0   = 60;            // pg/mL active disease (CONCERT median)
IL6_F_0   = 30;
LN_0      = 8;             // cm composite (significant lymphadenopathy)
PB_0      = 0.12;          // 12% of LN are plasmablasts
Bmem_0    = 100;
CRP_0     = 120;           // mg/L active iMCD
Hb_0      = 9.5;           // g/dL anemia of inflammation
IgG_0     = 4.0;           // hypergammaglobulinemia
VEGF_0    = 800;           // pg/mL
Anasarca_0= 2.5;
Plt_0     = 300;
mTOR_0    = 0.55;
HAZ_0     = 0;

$ODE
// ============================================================
// ---------- Siltuximab (SILT) PK -- Archetype 2 ----------
// ============================================================
double C_SILT     = CENT_SILT / V1_SILT;              // mg/L
double k10_SILT   = CL_SILT / V1_SILT;
double k12_SILT   = Q_SILT   / V1_SILT;
double k21_SILT   = Q_SILT   / V2_SILT;
dxdt_CENT_SILT = -(k10_SILT + k12_SILT) * CENT_SILT + k21_SILT * PERI_SILT;
dxdt_PERI_SILT =  k12_SILT * CENT_SILT - k21_SILT * PERI_SILT;

double C_SILT_nM   = C_SILT / MW_SILT * 1e9;                          // nM
double EFFECT_SILT = pow(C_SILT_nM, GAMMA_SILT) /
                     (pow(EC50_SILT, GAMMA_SILT) + pow(C_SILT_nM, GAMMA_SILT)); // fraction of IL-6 neutralized
double SILT_free_frac = 1.0 - EMAX_SILT * EFFECT_SILT;                // free (unbound) IL-6 fraction

// ============================================================
// ---------- Tocilizumab (TOCZ) PK -- Archetype 2 ----------
// ============================================================
double C_TOCZ     = CENT_TOCZ / V1_TOCZ;
double k10_TOCZ   = CL_TOCZ / V1_TOCZ;
double k12_TOCZ   = Q_TOCZ  / V1_TOCZ;
double k21_TOCZ   = Q_TOCZ  / V2_TOCZ;
dxdt_CENT_TOCZ = -(k10_TOCZ + k12_TOCZ) * CENT_TOCZ + k21_TOCZ * PERI_TOCZ;
dxdt_PERI_TOCZ =  k12_TOCZ * CENT_TOCZ - k21_TOCZ * PERI_TOCZ;

// NOTE: original hardcodes the same 145000 Da (IgG1) literal for tocilizumab
// molar conversion as it does for siltuximab, rather than a named MW_TOCZ
// parameter -- preserved exactly, not promoted to a new parameter.
double C_TOCZ_nM   = C_TOCZ / 145000 * 1e9;                           // nM
double EFFECT_TOCZ = pow(C_TOCZ_nM, GAMMA_TOCZ) /
                     (pow(EC50_TOCZ, GAMMA_TOCZ) + pow(C_TOCZ_nM, GAMMA_TOCZ)); // fraction of IL-6R blocked
double TOCZ_block  = 1.0 - EMAX_TOCZ * EFFECT_TOCZ;                   // fraction of IL-6R signal getting through

// ============================================================
// ---------- Sirolimus (SIRO) PK -- Archetype 3 ----------
// ============================================================
double k10_SIRO = (CL_SIRO * 24) / V1_SIRO;     // /d
double k12_SIRO = (Q_SIRO  * 24) / V1_SIRO;
double k21_SIRO = (Q_SIRO  * 24) / V2_SIRO;
double C_SIRO   = CENT_SIRO / V1_SIRO * 1000;    // ng/mL
dxdt_GUT_SIRO  = -KA_SIRO * GUT_SIRO;
dxdt_CENT_SIRO =  KA_SIRO * GUT_SIRO * F_SIRO - (k10_SIRO + k12_SIRO) * CENT_SIRO + k21_SIRO * PERI_SIRO;
dxdt_PERI_SIRO =  k12_SIRO * CENT_SIRO - k21_SIRO * PERI_SIRO;

double EFFECT_SIRO = pow(C_SIRO, GAMMA_SIRO) / (pow(EC50_SIRO, GAMMA_SIRO) + pow(C_SIRO, GAMMA_SIRO));
double Sirol_block  = 1.0 - EMAX_SIRO * EFFECT_SIRO;

// ============================================================
// ---------- Rituximab (RTX) PK -- Archetype 2 ----------
// ============================================================
double k10_RTX = CL_RTX / V1_RTX;
double k12_RTX = Q_RTX  / V1_RTX;
double k21_RTX = Q_RTX  / V2_RTX;
double C_RTX   = CENT_RTX / V1_RTX;
dxdt_CENT_RTX = -(k10_RTX + k12_RTX) * CENT_RTX + k21_RTX * PERI_RTX;
dxdt_PERI_RTX =  k12_RTX * CENT_RTX - k21_RTX * PERI_RTX;

// Bespoke (not Hill-shaped): the original CD20-lysis kill term is a plain
// first-order rate on concentration, no saturating EC50 anywhere.
double EFFECT_RTX = KILL_RTX * C_RTX;

// ============================================================
// ---------- Anakinra (ANA) PK -- Archetype 3 minus peripheral ----------
// (no EFFECT_ANA: the original never couples ANA_C/Cana to any disease/PD
//  equation anywhere in the file -- see refactor notes / UPSTREAM_ISSUES)
// ============================================================
double k10_ANA = (CL_ANA * 24) / V1_ANA;
double C_ANA   = CENT_ANA / V1_ANA * 1000;            // ng/mL
dxdt_GUT_ANA = -KA_ANA * 24 * GUT_ANA;
dxdt_CENT_ANA  =  KA_ANA * 24 * GUT_ANA * F_ANA - k10_ANA * CENT_ANA;

// ============================================================
// ---------- Ruxolitinib (RUX) PK -- Archetype 3 minus peripheral ----------
// PRESERVED BUG (not fixed, disclosed): the original own RUX_GUT/RUX_C
// ODEs read Ka_SIRO (sirolimus own absorption rate), not Ka_RUX (which is
// declared but never referenced in $ODE at all). KA_RUX above is therefore
// carried forward exactly as the original dead-but-declared parameter;
// KA_SIRO is the value that actually governs ruxolitinib absorption here,
// same as in the original. See refactor notes / UPSTREAM_ISSUES.
// ============================================================
double k10_RUX = (CL_RUX * 24) / V1_RUX;
double C_RUX   = CENT_RUX / V1_RUX * 1000;
dxdt_GUT_RUX = -KA_SIRO * 24 * GUT_RUX;
dxdt_CENT_RUX   =  KA_SIRO * 24 * GUT_RUX * F_RUX - k10_RUX * CENT_RUX;

double EFFECT_RUX = pow(C_RUX, GAMMA_RUX) / (pow(EC50_RUX, GAMMA_RUX) + pow(C_RUX, GAMMA_RUX));
double JAK_block   = 1.0 - EMAX_RUX * EFFECT_RUX;

// ============================================================
// ---------- CHOP doxorubicin (DOXO) PK -- Archetype 1 ----------
// ============================================================
double k10_DOXO = (CL_DOXO * 24) / V1_DOXO;
double C_DOXO = CENT_DOXO / V1_DOXO;
dxdt_CENT_DOXO = -k10_DOXO * CENT_DOXO;
double EFFECT_DOXO = EMAX_DOXO * pow(C_DOXO, GAMMA_DOXO) / (pow(EC50_DOXO, GAMMA_DOXO) + pow(C_DOXO, GAMMA_DOXO));

// ============================================================
// ---------- CHOP cyclophosphamide (CYC) PK -- Archetype 1 ----------
// ============================================================
double k10_CYC = (CL_CYC * 24) / V1_CYC;
double C_CYC  = CENT_CYC / V1_CYC;
dxdt_CENT_CYC = -k10_CYC * CENT_CYC;
double EFFECT_CYC = EMAX_CYC * pow(C_CYC, GAMMA_CYC) / (pow(EC50_CYC, GAMMA_CYC) + pow(C_CYC, GAMMA_CYC));

// ============================================================
// ---------- Prednisone (PRED) PK -- Archetype 3 minus peripheral ----------
// ============================================================
double k10_PRED = (CL_PRED * 24) / V1_PRED;
double C_PRED = CENT_PRED / V1_PRED * 1000;
dxdt_GUT_PRED = -KA_PRED * 24 * GUT_PRED;
dxdt_CENT_PRED   =  KA_PRED * 24 * GUT_PRED * F_PRED - k10_PRED * CENT_PRED;

double EFFECT_PRED = pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));
double Pred_block   = 1.0 - EMAX_PRED * EFFECT_PRED;

// ============================================================
// ---------- Disease side: IL-6 axis (untouched math, renamed refs only) ----
// ============================================================
double DiseaseDrive = kIL6_iMCD * (LN / LN_0) * (1 + 0.6 * (mTOR - 0.3));
double IL6_synth    = kIL6_base + DiseaseDrive * 1.0;       // pg/mL/d

// dead in the original (computed, never used downstream) -- preserved as-is
double IL6_F_calc  = IL6_T * SILT_free_frac;

dxdt_IL6_T = IL6_synth * Pred_block - kIL6_deg * IL6_T;
dxdt_IL6_F = IL6_synth * Pred_block * SILT_free_frac - kIL6_deg * IL6_F;

double IL6_signal = IL6_F * TOCZ_block * JAK_block;          // effective signal at gp130 + post-receptor

// ---------- mTORC1 ----------
dxdt_mTOR = kmTOR_in * (1 + 0.4 * (IL6_signal / 50)) * Sirol_block - kmTOR_out * mTOR;

// ---------- Plasmablast & B-cell ----------
double IL6_prol = IL6_signal / (IL6_signal + 10);            // saturating
double Bortz_kill = 0;       // placeholder if bortezomib added
double cyto_kill = EFFECT_DOXO + EFFECT_CYC;
dxdt_PB = kprol_PB * (IL6_prol + 0.6 * mTOR) * (1 - PB / 0.8) - (kdeath_PB + cyto_kill) * PB;

dxdt_Bmem = kprol_Bmem * (1 - Bmem/120) - (kdeath_Bmem + EFFECT_RTX) * Bmem;

// ---------- Lymph node size ----------
double LN_drive = 0.7 * (PB / 0.10) + 0.4 * mTOR + 0.2 * (IL6_signal / 30);
dxdt_LN = kgrow_LN * LN_drive * (1 - LN/LN_max) - (kshrink_LN + 0.5 * cyto_kill) * LN;

// ---------- Acute phase ----------
dxdt_CRP = kCRP_in * IL6_signal / (IL6_signal + EC50_CRP) - kCRP_out * CRP;

// hepcidin / anemia
double Hepc = IL6_signal / (IL6_signal + Hepcidin_EC50);
dxdt_Hb  = kHb_in * (1 - Hepc) * (1 - Hb / Hb_max) - kHb_out * Hb;

// IgG polyclonal
dxdt_IgG = kIgG_in * (PB / 0.10) * (1 - IgG / IgG_max) - (kIgG_out + 0.1 * EFFECT_RTX) * IgG;

// VEGF / Anasarca / TAFRO axis (sirolimus reduces)
dxdt_VEGF = kVEGF_in * IL6_signal / (IL6_signal + EC50_VEGF) * (0.7 + 0.6 * mTOR) - kVEGF_out * VEGF;
dxdt_Anasarca = kAnasarca_in * VEGF / (VEGF + EC50_VEGFan) - kAnasarca_out * Anasarca;

// Platelet -- TAFRO consumption
double TAFRO_consum = Plt_TAFRO_kill * (VEGF / 800) * (mTOR / 0.5);
dxdt_Plt = kPlt_in * (1 - Plt/Plt_max) - (kPlt_out + TAFRO_consum) * Plt;

// ---------- Hazard ----------
dxdt_HAZ = h0_OS * (1 + beta_CRP * CRP);

$TABLE
double IL6_serum = IL6_T;
double IL6_free  = IL6_F;
double LN_size   = LN;
double CRP_lab   = CRP;
double Hb_lab    = Hb;
double IgG_lab   = IgG;
double VEGF_lab  = VEGF;
double Plt_lab   = Plt;
double Anasarca_lab = Anasarca;
double mTOR_act  = mTOR;
double Survival  = exp(-HAZ);

// CDCN response composite (lower is better)
double CDCN_resp = 0.25*(IL6_serum/60) + 0.25*(CRP/120) + 0.20*(LN/8) +
                   0.10*(1 - Hb/12) + 0.10*(IgG/4) + 0.10*(Anasarca/2.5);

$CAPTURE C_SILT EFFECT_SILT C_TOCZ EFFECT_TOCZ C_SIRO EFFECT_SIRO C_RTX EFFECT_RTX
         C_ANA C_RUX EFFECT_RUX C_DOXO EFFECT_DOXO C_CYC EFFECT_CYC C_PRED EFFECT_PRED
         IL6_serum IL6_free LN_size CRP_lab Hb_lab IgG_lab VEGF_lab Plt_lab Anasarca_lab mTOR_act
         Survival CDCN_resp
'

mod <- mcode("imcd_refactored", code)

# =============================================================================
# Treatment scenarios (identical dosing to the original; compartment names
# updated to the refactored convention: SILT_C->CENT_SILT, TOCZ_C->CENT_TOCZ,
# SIRO_GUT->GUT_SIRO, RTX_C->CENT_RTX, ANA_SC->GUT_ANA, DOXO_C->CENT_DOXO,
# CYC_C->CENT_CYC, PRED_GUT->GUT_PRED)
# =============================================================================
sim_one <- function(label, ev) {
  mrgsim(mod, events = ev, end = 365, delta = 1) %>%
    as_tibble() %>% mutate(scenario = label)
}

# ---- Scenario 1: Untreated natural history -----
ev_none <- ev(amt = 0, cmt = "CENT_SILT", time = 0)
s1 <- sim_one("S1: Untreated", ev_none)

# ---- Scenario 2: Siltuximab 11 mg/kg IV q3w (CONCERT trial label) -----
# 70 kg -> 770 mg per dose
ev_silt <- ev(amt = 770, cmt = "CENT_SILT", time = seq(0, 360, 21))
s2 <- sim_one("S2: Siltuximab 11 mg/kg q3w", ev_silt)

# ---- Scenario 3: Tocilizumab 8 mg/kg IV q2w (Nishimoto regimen) -----
ev_tocz <- ev(amt = 560, cmt = "CENT_TOCZ", time = seq(0, 360, 14))
s3 <- sim_one("S3: Tocilizumab 8 mg/kg q2w", ev_tocz)

# ---- Scenario 4: Sirolimus PO QD (TAFRO Fajgenbaum 2019 JCI) -----
ev_siro <- ev(amt = 2, cmt = "GUT_SIRO", ii = 1, addl = 364)
s4 <- sim_one("S4: Sirolimus 2 mg QD (TAFRO)", ev_siro)

# ---- Scenario 5: Rituximab 375 mg/m2 weekly x4 + prednisone -----
ev_rtx_pred <- bind_rows(
  ev(amt = 750, cmt = "CENT_RTX", time = c(0, 7, 14, 21)),     # 1.5 m2 ~750 mg
  ev(amt = 60,  cmt = "GUT_PRED", ii = 1, addl = 27, time = 0)
)
s5 <- sim_one("S5: Rituximab + Prednisone", ev_rtx_pred)

# ---- Scenario 6: CHOP-like + Siltuximab (induction) -----
ev_chop_silt <- bind_rows(
  ev(amt = 75 * 1.5,  cmt = "CENT_DOXO", time = seq(0, 6*21, 21)),
  ev(amt = 750 * 1.5, cmt = "CENT_CYC",  time = seq(0, 6*21, 21)),
  ev(amt = 100,       cmt = "GUT_PRED", ii = 1, addl = 4, time = seq(0, 6*21, 21)),
  ev(amt = 770,       cmt = "CENT_SILT", time = seq(0, 360, 21))
)
s6 <- sim_one("S6: CHOP + Siltuximab", ev_chop_silt)

# ---- Scenario 7: Triple therapy (Sirolimus + Siltuximab + Anakinra) for refractory TAFRO -----
ev_triple <- bind_rows(
  ev(amt = 770, cmt = "CENT_SILT", time = seq(0, 360, 21)),
  ev(amt = 2,   cmt = "GUT_SIRO", ii = 1, addl = 364),
  ev(amt = 100, cmt = "GUT_ANA", ii = 1, addl = 364)
)
s7 <- sim_one("S7: Siltuximab + Sirolimus + Anakinra", ev_triple)

# Combine and quick plot
all_sim <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

# Example plot
if (interactive()) {
  ggplot(all_sim, aes(time, CRP_lab, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(title = "iMCD QSP (refactored PK): CRP trajectory by treatment scenario",
         x = "Day", y = "CRP (mg/L)") +
    theme_bw()
}

# =============================================================================
# Calibration notes (representative literature anchors, unchanged from original)
# =============================================================================
# 1) CONCERT trial (van Rhee 2014 Lancet Oncol): siltuximab 11 mg/kg q3w
#    durable tumor + symptomatic response 34% vs 0% placebo (P=0.0012).
#    CRP normalized to < 10 mg/L by week 6; IL-6 increased due to clearance
#    block (total IL-6 rises because free IL-6 is sequestered).
#
# 2) Nishimoto 2005 Blood: tocilizumab 8 mg/kg q2w —
#    CRP < 10 mg/L by week 2; IgG, fibrinogen, ESR normalize by week 4;
#    lymph node shrinkage by 50% at 6 months in ~60% of patients.
#
# 3) Fajgenbaum 2019 JCI Insight: sirolimus rescues 3/3 refractory TAFRO patients,
#    target trough 6-14 ng/mL — rationale for mTORC1 / Tfh hyperactivation.
#
# 4) van Rhee 2018 Blood Adv consensus treatment guidelines — 1st line siltuximab
#    (or tocilizumab where siltuximab unavailable). Sirolimus for refractory.
#
# 5) Dispenzieri & Fajgenbaum 2020 Blood: 5-year OS ≈ 65% (TAFRO worse, IPL better).
# =============================================================================
