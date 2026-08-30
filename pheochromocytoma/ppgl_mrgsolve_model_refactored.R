## ============================================================
##  Pheochromocytoma / Paraganglioma (PPGL) — QSP mrgsolve Model
##  REFACTORED sibling of ppgl_mrgsolve_model.R — pluggable-PK naming
##  convention applied to all 5 compounds (Phenoxybenzamine/PHE,
##  Doxazosin/DOX, Metyrosine/MET, Beta-blocker/BB, Sunitinib/SUNIT).
##  See ppgl_refactor_notes.md for the full account. All numeric
##  parameter values are copied verbatim from the original; only PK
##  structure/naming and the Hill-interface expression are
##  reorganized. The original is untouched.
## ============================================================
##
##  Catecholamine Biosynthesis · Adrenergic Signaling · CV Effects
##  Drug PK/PD: Phenoxybenzamine / Doxazosin / Metyrosine /
##              Beta-blocker / Sunitinib
## ============================================================
##
##  Compartments (20 ODEs) — renamed, same order/indices as original:
##   PK  :  GUT_PHE, CENT_PHE, PERI_PHE, CENT_DOX, CENT_MET, CENT_BB,
##          CENT_SUNIT, PERI_SUNIT
##   Synth: TH_act, NE_store, NE_plasma, EPI_plasma
##   Tumor: TUMvol, VEGF_tum
##   CV  :  SBP, DBP, HR
##   Metabol: GLU, FFA
##   Biomark: CgA_plasma, NMN_plasma
##
##  Treatment Scenarios (6):
##   0 – No treatment (natural history)
##   1 – Preoperative Phenoxybenzamine (60 mg/d) → surgery Day 14
##   2 – Preoperative Doxazosin (16 mg/d) → surgery Day 14
##   3 – Phenoxybenzamine + Metyrosine (2 g/d) + Beta-blocker combo
##   4 – Sunitinib 37.5 mg/d (malignant/metastatic PPGL)
##   5 – Metyrosine monotherapy (inoperable, symptomatic control)
##
##  Calibration references:
##   Kinney et al. 2000 J Urol (phenoxybenzamine vs doxazosin preop)
##   Steinsapir et al. 1997 Anesthesiology (hemodynamic preop prep)
##   Baudin et al. 2014 J Clin Endocrinol Metab (sunitinib PPGL)
##   Fassnacht et al. 2020 Eur J Endocrinol (malignant PPGL outcomes)
##   Molenaar 2014 (metyrosine PK/PD)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

code <- '
$PARAM
// ─── PHENOXYBENZAMINE (PHE) PK — archetype 3 (depot+central+peripheral) ──
KA_PHE   = 0.693     // h-1 GI absorption rate
F_PHE    = 0.27      // oral bioavailability 27%
CL_PHE   = 5.8       // L/h central clearance
V1_PHE   = 45        // L central volume
V2_PHE   = 180       // L peripheral volume
Q_PHE    = 10        // L/h intercompartment CL
KOUT_PHE = 0.0035    // h-1 irreversible receptor unbinding (new synthesis) — vestigial, unused in original $ODE, preserved as-is (see refactor notes)
EMAX_PHE = 1.0       // fractional alpha-blockade at saturation (rename; original PHE_eff had no separate Emax<1)
EC50_PHE = 0.012     // µmol/L PHE IC50 for alpha-blockade (promoted from an inline literal in the original)
GAMMA_PHE = 1        // Hill coefficient (rename, not a fit; original had a plain ratio)

// ─── DOXAZOSIN (DOX) PK — archetype 1 (single compartment, direct dosing) ─
KA_DOX   = 0.21      // h-1 (Tmax ~2-3 h) — vestigial: only multiplies the always-zero DOSE_DOX flag (see notes)
F_DOX    = 0.65      // 65% bioavailability — vestigial, same reason
CL_DOX   = 3.2       // L/h — vestigial $PARAM, never used (original derives clearance from T50_DOX instead, preserved as-is)
V1_DOX   = 80        // L
T50_DOX  = 22        // h terminal half-life (t½); actual clearance basis
EMAX_DOX = 1.0       // fractional alpha-blockade at saturation (rename; original DOX_eff had no separate Emax<1)
EC50_DOX = 0.002     // µmol/L DOX IC50 (promoted from an inline literal in the original)
GAMMA_DOX = 1        // Hill coefficient (rename, not a fit)

// ─── METYROSINE (MET) PK — archetype 1 (single compartment, direct dosing) ─
KA_MET   = 0.70      // h-1 — vestigial: only multiplies the always-zero DOSE_MET flag (see notes)
F_MET    = 0.85      // 85% oral F — vestigial, same reason
CL_MET   = 4.5       // L/h (renal dominant)
V1_MET   = 40        // L
EMAX_MET = 0.80      // max fractional TH inhibition (promoted from the original inline "0.80 *" multiplier)
EC50_MET = 85         // µmol/L metyrosine EC50 at TH (renamed from IC50_MET, same value)
GAMMA_MET = 1        // Hill coefficient (rename, not a fit)

// ─── BETA-BLOCKER (BB) PK — archetype 1 (single compartment, direct dosing) ─
KA_BB    = 0.80      // h-1 propranolol — vestigial: only multiplies the always-zero DOSE_BB flag (see notes)
F_BB     = 0.36      // 36% hepatic FPE — vestigial, same reason
CL_BB    = 50        // L/h
V1_BB    = 260       // L
EMAX_BB  = 0.40      // max fractional HR-response blunting (promoted from the original inline "0.40 *" multiplier)
EC50_BB  = 0.022     // µmol/L propranolol EC50 (renamed from IC50_BB, same value)
GAMMA_BB = 1         // Hill coefficient (rename, not a fit)

// ─── SUNITINIB (SUNIT) PK — archetype 2 (2-cmpt, no depot, direct dosing) ─
KA_SUNIT = 0.28      // h-1 (Tmax ~6-12h) — vestigial: only multiplies the always-zero DOSE_SUNIT flag (see notes)
F_SUNIT  = 0.60      // vestigial, same reason
CL_SUNIT = 34        // L/h (Lindauer 2010)
V1_SUNIT = 2230      // L central
V2_SUNIT = 1900      // L peripheral
Q_SUNIT  = 6.5       // L/h
EMAX_SUNIT = 0.65    // max fractional tumor-growth-rate inhibition (promoted from the original inline "0.65 *" multiplier)
EC50_SUNIT = 0.018   // µmol/L sunitinib EC50 at VEGFR (renamed from IC50_SUNIT_VEGFR, same value)
GAMMA_SUNIT = 1      // Hill coefficient (rename, not a fit)

// ─── CATECHOLAMINE BIOSYNTHESIS ──────────────────────────────
ksynth_NE  = 0.35    // nmol/L/h NE baseline synthesis rate
ksynth_EPI = 0.12    // nmol/L/h EPI (adrenal-dependent)
kstor      = 0.80    // h-1 VMAT2 loading rate into granules
krel       = 0.045   // h-1 baseline secretion rate from granules
krelease   = 1.0     // tumor-driven secretion multiplier (x baseline)
kdeg_NE    = 0.62    // h-1 NE plasma elimination (COMT+MAO+NET)
kdeg_EPI   = 0.55    // h-1 EPI plasma elimination
NE0        = 1.8     // nmol/L baseline plasma NE (healthy)
EPI0       = 0.28    // nmol/L baseline plasma EPI (healthy)
NMN_factor = 0.82    // NMN/NE conversion ratio (COMT)

// ─── TH ACTIVITY ─────────────────────────────────────────────
TH_base    = 1.0     // normalized baseline TH activity
kTH_deg    = 0.05    // h-1 TH protein turnover
kTH_synth  = 0.05    // h-1 TH protein synthesis (zeroth order)

// ─── TUMOR DYNAMICS ──────────────────────────────────────────
TUM0       = 30.0    // mL initial tumor volume
kgrowth    = 0.003   // day-1 net exponential growth (slow)
kgrowth_h  = 0.000125 // h-1
kdeath     = 0.0     // h-1 baseline tumor death
VEGF0      = 120     // pg/mL baseline VEGF
kVEGF_synth= 0.10    // VEGF production per mL tumor volume
kVEGF_deg  = 0.20    // h-1 VEGF clearance

// ─── CARDIOVASCULAR BASELINE ─────────────────────────────────
SBP0       = 155     // mmHg (PPGL patient typical elevated BP)
DBP0       = 98      // mmHg
HR0        = 92      // bpm
E_NE_SBP   = 28.0   // mmHg per nmol/L NE above normal
E_NE_HR    = 12.0   // bpm per nmol/L NE above normal
E_EPI_SBP  = 18.0   // mmHg per nmol/L EPI above normal
E_EPI_HR   = 22.0   // bpm per nmol/L EPI above normal
kSBP       = 0.15    // h-1 SBP equilibration
kHR        = 0.25    // h-1 HR equilibration
alpha_block_max = 0.92 // max BP reduction by combined alpha blockade (disease-level Emax, applied to the PHE+DOX combination, not to either compound individually — see refactor notes)
IC50_alpha  = 0.18   // µmol/L — vestigial in the original (never referenced anywhere), preserved as-is, untouched (out of the 5-compound refactor scope)

// ─── GLUCOSE & FFA ───────────────────────────────────────────
GLU0       = 5.5     // mmol/L baseline glucose
kGLU       = 0.15    // h-1 glucose equilibration
E_NE_GLU   = 0.8     // mmol/L per nmol/L NE
FFA0       = 0.55    // mmol/L baseline free fatty acids
kFFA       = 0.12    // h-1 FFA equilibration
E_EPI_FFA  = 0.25    // mmol/L per nmol/L EPI

// ─── CgA BIOMARKER ───────────────────────────────────────────
CgA0       = 180     // ng/mL elevated CgA (PPGL typical)
kCgA_synth = 0.12    // h-1 proportional to tumor volume
kCgA_deg   = 0.05    // h-1 plasma CgA clearance

// ─── DOSE FLAGS (0=off, 1=on) ────────────────────────────────
DOSE_PHE   = 0       // vestigial: PHE is dosed via ev() bolus into GUT_PHE, this flag is never read in $ODE
DOSE_DOX   = 0       // vestigial: never set nonzero anywhere in the original's own scenarios (DOX dosed via ev() bolus directly into CENT_DOX)
DOSE_MET   = 0       // vestigial, same reason (MET dosed via ev() bolus directly into CENT_MET)
DOSE_BB    = 0       // vestigial, same reason (BB dosed via ev() bolus directly into CENT_BB)
DOSE_SUNIT = 0       // vestigial, same reason (SUNIT dosed via ev() bolus directly into CENT_SUNIT)
SURGERY    = 0       // 1 = tumor removed at tSURG

$CMT
// PK compartments
@annotated
GUT_PHE    : Phenoxybenzamine GI compartment [µmol]
CENT_PHE   : Phenoxybenzamine central [µmol/L]
PERI_PHE   : Phenoxybenzamine peripheral [µmol/L]
CENT_DOX   : Doxazosin central [µmol/L]
CENT_MET   : Metyrosine central [µmol/L]
CENT_BB    : Beta-blocker central [µmol/L]
CENT_SUNIT : Sunitinib central [µmol/L]
PERI_SUNIT : Sunitinib peripheral [µmol/L]
// Biosynthesis
TH_act   : TH enzyme activity [normalized]
NE_store : NE stored in chromaffin granules [nmol/L_equiv]
NE_plasma: Plasma NE [nmol/L]
EPI_plasma:Plasma EPI [nmol/L]
// Tumor
TUMvol   : Tumor volume [mL]
VEGF_tum : Plasma VEGF [pg/mL]
// CV
SBP      : Systolic BP [mmHg]
DBP      : Diastolic BP [mmHg]
HR       : Heart rate [bpm]
// Metabolic
GLU      : Plasma glucose [mmol/L]
FFA      : Plasma FFA [mmol/L]
// Biomarker
CgA_plasma : Plasma CgA [ng/mL]

$MAIN
// Initial conditions
GUT_PHE_0    = 0;
CENT_PHE_0   = 0;
PERI_PHE_0   = 0;
CENT_DOX_0   = 0;
CENT_MET_0   = 0;
CENT_BB_0    = 0;
CENT_SUNIT_0 = 0;
PERI_SUNIT_0 = 0;
TH_act_0     = TH_base;
NE_store_0   = NE0 / krel * kstor;   // steady-state granule pool
NE_plasma_0  = NE0;
EPI_plasma_0 = EPI0;
TUMvol_0     = TUM0;
VEGF_tum_0   = VEGF0;
SBP_0        = SBP0;
DBP_0        = DBP0;
HR_0         = HR0;
GLU_0        = GLU0;
FFA_0        = FFA0;
CgA_plasma_0 = CgA0;

$ODE
// ──────────────────────────────────────────────
// 1. Phenoxybenzamine PK (oral, 2-cmpt, irreversible α-blocker)
//    Archetype 3 (depot + central + peripheral). CENT_PHE/PERI_PHE are
//    already concentration-valued states in the original (every ODE term
//    is pre-divided by volume) — preserved exactly, so C_PHE is a direct
//    alias, not an amount/V division (see refactor notes).
// ──────────────────────────────────────────────
dxdt_GUT_PHE  = -KA_PHE * GUT_PHE;
dxdt_CENT_PHE =  KA_PHE * F_PHE * GUT_PHE / V1_PHE
               - (CL_PHE/V1_PHE) * CENT_PHE
               - (Q_PHE/V1_PHE)  * CENT_PHE
               + (Q_PHE/V2_PHE)  * PERI_PHE;
dxdt_PERI_PHE =  (Q_PHE/V1_PHE)  * CENT_PHE
               - (Q_PHE/V2_PHE)  * PERI_PHE;
double C_PHE = CENT_PHE;   // exposed, PD-facing: plasma phenoxybenzamine (µmol/L)
double EFFECT_PHE = EMAX_PHE * pow(C_PHE, GAMMA_PHE) / (pow(EC50_PHE, GAMMA_PHE) + pow(C_PHE, GAMMA_PHE));

// ──────────────────────────────────────────────
// 2. Doxazosin PK (1-cmpt, competitive α₁-blocker)
//    Archetype 1 (single compartment, no depot). CL is derived from the
//    terminal half-life (T50_DOX), exactly as the original did (CL_DOX
//    $PARAM itself is dead/unused in both the original and here).
// ──────────────────────────────────────────────
double CL_DOX_eff = 0.693 / T50_DOX * V1_DOX;
dxdt_CENT_DOX = KA_DOX * F_DOX * DOSE_DOX / V1_DOX
             - (CL_DOX_eff / V1_DOX) * CENT_DOX;
double C_DOX = CENT_DOX;   // exposed, PD-facing: plasma doxazosin (µmol/L)
double EFFECT_DOX = EMAX_DOX * pow(C_DOX, GAMMA_DOX) / (pow(EC50_DOX, GAMMA_DOX) + pow(C_DOX, GAMMA_DOX));

// ──────────────────────────────────────────────
// 3. Metyrosine PK (TH inhibitor, 1-cmpt) — archetype 1
// ──────────────────────────────────────────────
dxdt_CENT_MET = KA_MET * F_MET * DOSE_MET / V1_MET
             - (CL_MET / V1_MET) * CENT_MET;
double C_MET = CENT_MET;   // exposed, PD-facing: plasma metyrosine (µmol/L)
double EFFECT_MET = EMAX_MET * pow(C_MET, GAMMA_MET) / (pow(EC50_MET, GAMMA_MET) + pow(C_MET, GAMMA_MET));

// ──────────────────────────────────────────────
// 4. Beta-blocker PK (1-cmpt) — archetype 1
// ──────────────────────────────────────────────
dxdt_CENT_BB = KA_BB * F_BB * DOSE_BB / V1_BB
            - (CL_BB / V1_BB) * CENT_BB;
double C_BB = CENT_BB;   // exposed, PD-facing: plasma beta-blocker (µmol/L)
double EFFECT_BB = EMAX_BB * pow(C_BB, GAMMA_BB) / (pow(EC50_BB, GAMMA_BB) + pow(C_BB, GAMMA_BB));

// ──────────────────────────────────────────────
// 5. Sunitinib PK (2-cmpt, malignant PPGL) — archetype 2
// ──────────────────────────────────────────────
dxdt_CENT_SUNIT = KA_SUNIT * F_SUNIT * DOSE_SUNIT / V1_SUNIT
               - (CL_SUNIT/V1_SUNIT) * CENT_SUNIT
               - (Q_SUNIT/V1_SUNIT)  * CENT_SUNIT
               + (Q_SUNIT/V2_SUNIT)  * PERI_SUNIT;
dxdt_PERI_SUNIT =  (Q_SUNIT/V1_SUNIT)  * CENT_SUNIT
               - (Q_SUNIT/V2_SUNIT)  * PERI_SUNIT;
double C_SUNIT = CENT_SUNIT;   // exposed, PD-facing: plasma sunitinib (µmol/L)
double EFFECT_SUNIT = EMAX_SUNIT * pow(C_SUNIT, GAMMA_SUNIT) / (pow(EC50_SUNIT, GAMMA_SUNIT) + pow(C_SUNIT, GAMMA_SUNIT));

// ──────────────────────────────────────────────
// 6. TH activity (inhibited by metyrosine)
// ──────────────────────────────────────────────
double TH_target = TH_base * (1 - EFFECT_MET);  // EFFECT_MET already carries the 0.80 max-inhibition scale (EMAX_MET)
dxdt_TH_act = kTH_synth * TH_target - kTH_deg * TH_act;

// ──────────────────────────────────────────────
// 7. NE & EPI in granule store
// ──────────────────────────────────────────────
// Tumor-driven synthesis upregulation
double tumFactor = TUMvol / TUM0;   // relative to baseline tumor
double synth_NE  = ksynth_NE  * TH_act * tumFactor;
double synth_EPI = ksynth_EPI * TH_act * tumFactor;
// Release from granules (basal + tumor-driven)
double krel_eff  = krel * (1 + krelease * (tumFactor - 1));
dxdt_NE_store = kstor * synth_NE - krel_eff * NE_store;

// ──────────────────────────────────────────────
// 8. Plasma NE
// ──────────────────────────────────────────────
dxdt_NE_plasma = krel_eff * NE_store - kdeg_NE * NE_plasma;

// ──────────────────────────────────────────────
// 9. Plasma EPI
// ──────────────────────────────────────────────
dxdt_EPI_plasma = synth_EPI - kdeg_EPI * EPI_plasma;

// ──────────────────────────────────────────────
// 10. Tumor volume (exponential growth with TKI inhibition)
// ──────────────────────────────────────────────
double kgrowth_eff = kgrowth_h * (1 - EFFECT_SUNIT);  // EFFECT_SUNIT already carries the 0.65 max-inhibition scale (EMAX_SUNIT)
double surgery_factor = (SURGERY > 0) ? 0 : 1;  // tumor removed post-op
dxdt_TUMvol = kgrowth_eff * TUMvol * surgery_factor;

// ──────────────────────────────────────────────
// 11. Tumor VEGF (proportional to tumor volume)
// ──────────────────────────────────────────────
dxdt_VEGF_tum = kVEGF_synth * TUMvol - kVEGF_deg * VEGF_tum;

// ──────────────────────────────────────────────
// 12-13. BP (SBP, DBP) — driven by NE/EPI, inhibited by alpha-blockers
// ──────────────────────────────────────────────
// Combined alpha-blockade effect (PHE irreversible + DOX competitive).
// EFFECT_PHE/EFFECT_DOX are each independently exposed/redirectable
// (EMAX=1, i.e. unscaled saturating ratios); alpha_block_max is a shared
// disease-level scale applied AFTER the two-drug Bliss-independence
// combination below, exactly as the original did — it cannot be folded
// into either compound's own EFFECT_<STEM> without changing the combined
// math when both drugs are present (see refactor notes).
double alpha_block = 1 - (1 - EFFECT_PHE) * (1 - EFFECT_DOX);  // combined effect

double NE_excess = NE_plasma - NE0;
double EPI_excess = EPI_plasma - EPI0;
double SBP_target = SBP0
  + E_NE_SBP  * (NE_excess  > 0 ? NE_excess  : 0)
  + E_EPI_SBP * (EPI_excess > 0 ? EPI_excess : 0);
SBP_target = SBP_target * (1 - alpha_block_max * alpha_block);

double DBP_target = DBP0 + 0.65 * (SBP_target - SBP0);

dxdt_SBP = kSBP * (SBP_target - SBP);
dxdt_DBP = kSBP * (DBP_target - DBP);

// ──────────────────────────────────────────────
// 14. Heart rate — NE/EPI drive; beta-blocker inhibition
// ──────────────────────────────────────────────
double HR_target = HR0
  + E_NE_HR  * (NE_excess  > 0 ? NE_excess  : 0)
  + E_EPI_HR * (EPI_excess > 0 ? EPI_excess : 0);
HR_target = HR_target * (1 - EFFECT_BB);  // EFFECT_BB already carries the 0.40 max-blunting scale (EMAX_BB)

dxdt_HR = kHR * (HR_target - HR);

// ──────────────────────────────────────────────
// 15. Glucose (stress hyperglycemia from EPI/NE)
// ──────────────────────────────────────────────
double GLU_target = GLU0
  + E_NE_GLU * (NE_excess > 0 ? NE_excess : 0)
  + 0.5 * (EPI_excess > 0 ? EPI_excess : 0);
dxdt_GLU = kGLU * (GLU_target - GLU);

// ──────────────────────────────────────────────
// 16. Free fatty acids (lipolysis via β₃/β₂)
// ──────────────────────────────────────────────
double FFA_target = FFA0 + E_EPI_FFA * (EPI_excess > 0 ? EPI_excess : 0);
dxdt_FFA = kFFA * (FFA_target - FFA);

// ──────────────────────────────────────────────
// 17. CgA plasma — tumor marker
// ──────────────────────────────────────────────
dxdt_CgA_plasma = kCgA_synth * TUMvol - kCgA_deg * CgA_plasma;

$TABLE
// Derived outputs
capture PHE_Conc   = CENT_PHE;               // µmol/L plasma PHE
capture DOX_Conc   = CENT_DOX;
capture MET_Conc   = CENT_MET;
capture BB_Conc    = CENT_BB;
capture SUNIT_Conc = CENT_SUNIT;
capture TH_activity = TH_act;             // TH relative activity
capture NE_gran    = NE_store;
capture NE_pl      = NE_plasma;           // nmol/L
capture EPI_pl     = EPI_plasma;          // nmol/L
capture NMN_pl     = NE_plasma * NMN_factor; // proxy plasma normetanephrine
capture MN_pl      = EPI_plasma * 1.15;   // proxy plasma metanephrine
capture SBP_mmHg   = SBP;
capture DBP_mmHg   = DBP;
capture MAP_mmHg   = DBP + (SBP - DBP) / 3.0;
capture HR_bpm     = HR;
capture Tumor_mL   = TUMvol;
capture VEGF_pg    = VEGF_tum;
capture Glucose_mM = GLU;
capture FFA_mM     = FFA;
capture CgA_ng     = CgA_plasma;
capture AlphaBlock_frac = alpha_block;

$CAPTURE C_PHE C_DOX C_MET C_BB C_SUNIT EFFECT_PHE EFFECT_DOX EFFECT_MET EFFECT_BB EFFECT_SUNIT
'

# ── Compile model ──────────────────────────────────────────────
mod <- mcode("PPGL_QSP_refactored", code)

# ── Helper: dose events ───────────────────────────────────────
# Phenoxybenzamine 10 mg PO BID  (≈ 0.027 µmol/kg/dose, ~60 mg/d total)
dose_PHE <- function(dose_mg_d = 60, freq_h = 12, dur_h = 336) {
  ev(amt = dose_mg_d/freq_h / 303.8,  # mg/dose → µmol/L (approx 70 kg BW)
     ii  = freq_h,
     addl= round(dur_h/freq_h) - 1,
     cmt = "GUT_PHE",
     time = 0)
}
# Doxazosin 8 mg PO QD
dose_DOX <- function(dose_mg_d = 16, freq_h = 24, dur_h = 336) {
  ev(amt = dose_mg_d/freq_h / 451.5,
     ii  = freq_h,
     addl= round(dur_h/freq_h) - 1,
     cmt = "CENT_DOX",
     time = 0)
}
# Metyrosine 500 mg PO QID
dose_MET <- function(dose_mg_d = 2000, freq_h = 6, dur_h = 336) {
  ev(amt = dose_mg_d/freq_h / 195.2,
     ii  = freq_h,
     addl= round(dur_h/freq_h) - 1,
     cmt = "CENT_MET",
     time = 0)
}
# Propranolol 40 mg PO TID  (β-blocker)
dose_BB <- function(dose_mg_d = 120, freq_h = 8, dur_h = 336) {
  ev(amt = dose_mg_d/freq_h / 259.3,
     ii  = freq_h,
     addl= round(dur_h/freq_h) - 1,
     cmt = "CENT_BB",
     time = 0)
}
# Sunitinib 37.5 mg/d continuous (malignant PPGL)
dose_SUNIT <- function(dose_mg_d = 37.5, freq_h = 24, dur_h = 2160) {
  ev(amt = dose_mg_d / 532.6,
     ii  = freq_h,
     addl= round(dur_h/freq_h) - 1,
     cmt = "CENT_SUNIT",
     time = 0)
}

# ── Simulation function ───────────────────────────────────────
run_scenario <- function(scenario = 0, end_h = 720, delta_h = 1) {
  params <- c(DOSE_PHE=0, DOSE_DOX=0, DOSE_MET=0, DOSE_BB=0, DOSE_SUNIT=0, SURGERY=0)
  evs <- NULL

  if (scenario == 1) {
    # Preop phenoxybenzamine → surgery at Day 14
    params["SURGERY"] <- 0
    evs <- dose_PHE(dose_mg_d=60, dur_h=336)
    # surgery event (remove tumor): set SURGERY=1 at 336h via parameter table
  } else if (scenario == 2) {
    # Preop doxazosin → surgery at Day 14
    evs <- dose_DOX(dose_mg_d=16, dur_h=336)
  } else if (scenario == 3) {
    # Triple therapy: PHE + Metyrosine + BB
    evs <- c(dose_PHE(dose_mg_d=60, dur_h=336),
             dose_MET(dose_mg_d=2000, dur_h=336),
             dose_BB(dose_mg_d=120, dur_h=336))
  } else if (scenario == 4) {
    # Sunitinib malignant PPGL (90 days)
    end_h <- 2160
    evs <- dose_SUNIT(dose_mg_d=37.5, dur_h=2160)
  } else if (scenario == 5) {
    # Metyrosine monotherapy
    evs <- dose_MET(dose_mg_d=2000, dur_h=720)
  }

  if (!is.null(evs)) {
    out <- mrgsim(mod, events=evs, end=end_h, delta=delta_h) %>% as.data.frame()
  } else {
    out <- mrgsim(mod, end=end_h, delta=delta_h) %>% as.data.frame()
  }
  out$scenario <- scenario
  out$scenario_label <- c(
    "0"="No treatment",
    "1"="Preop PHE → Surgery",
    "2"="Preop DOX → Surgery",
    "3"="PHE + MET + BB (triple)",
    "4"="Sunitinib (malignant)",
    "5"="Metyrosine monotherapy"
  )[as.character(scenario)]
  out
}

# ── Run all scenarios ─────────────────────────────────────────
set.seed(2024)
results <- lapply(0:5, function(s) run_scenario(scenario=s, end_h=720))
sim_df  <- bind_rows(results)

# ── Summary statistics ────────────────────────────────────────
summary_df <- sim_df %>%
  filter(time %in% c(0, 72, 168, 336, 504, 720)) %>%
  group_by(scenario_label, time) %>%
  summarise(
    SBP   = round(mean(SBP_mmHg), 1),
    DBP   = round(mean(DBP_mmHg), 1),
    MAP   = round(mean(MAP_mmHg), 1),
    HR    = round(mean(HR_bpm), 1),
    NMN   = round(mean(NMN_pl), 2),
    MN    = round(mean(MN_pl), 3),
    CgA   = round(mean(CgA_ng), 0),
    TumV  = round(mean(Tumor_mL), 1),
    .groups="drop"
  )
print(summary_df)

# ── Key calibration checks ────────────────────────────────────
cat("\n=== Calibration Check ===\n")
# Scenario 1: PHE preop at Day 14
scen1_14d <- sim_df %>% filter(scenario==1, abs(time - 336) < 2)
cat(sprintf("S1 (PHE preop) Day14  SBP: %.0f mmHg (target <140)\n",
            mean(scen1_14d$SBP_mmHg)))
cat(sprintf("S1 (PHE preop) Day14  HR:  %.0f bpm\n",
            mean(scen1_14d$HR_bpm)))
cat(sprintf("S1 (PHE preop) Day14  NMN: %.2f nmol/L (target <3× ULN=0.9)\n",
            mean(scen1_14d$NMN_pl)))

# Scenario 3: Triple therapy
scen3_14d <- sim_df %>% filter(scenario==3, abs(time - 336) < 2)
cat(sprintf("S3 (triple)    Day14  SBP: %.0f mmHg\n", mean(scen3_14d$SBP_mmHg)))
cat(sprintf("S3 (triple)    Day14  TH_act: %.2f (target <0.50)\n",
            mean(scen3_14d$TH_activity)))

# Scenario 4: Sunitinib at Day 30
scen4_d30 <- sim_df %>% filter(scenario==4, abs(time - 720) < 2)
cat(sprintf("S4 (sunitinib) Day30  TumVol: %.1f mL (baseline=30, target <25)\n",
            mean(scen4_d30$Tumor_mL)))
cat(sprintf("S4 (sunitinib) Day30  VEGF:   %.0f pg/mL\n",
            mean(scen4_d30$VEGF_pg)))

# ── Plotting ──────────────────────────────────────────────────
if (requireNamespace("ggplot2", quietly=TRUE)) {
  library(ggplot2)
  pal6 <- c("#E53935","#1E88E5","#43A047","#FB8C00","#8E24AA","#00ACC1")

  # Plot 1: SBP over time (all scenarios except Sunitinib which is longer)
  p1 <- sim_df %>%
    filter(scenario %in% 0:3, time <= 720) %>%
    ggplot(aes(x=time/24, y=SBP_mmHg, color=scenario_label)) +
    geom_line(size=1) +
    geom_hline(yintercept=140, linetype="dashed", color="gray50") +
    annotate("text", x=0.5, y=138, label="Target SBP <140", size=3, hjust=0, color="gray40") +
    scale_color_manual(values=pal6[1:4]) +
    labs(title="PPGL QSP — Systolic Blood Pressure",
         x="Time (days)", y="SBP (mmHg)", color="Scenario") +
    theme_bw(base_size=12)
  print(p1)

  # Plot 2: Plasma NMN (normetanephrine biomarker)
  p2 <- sim_df %>%
    filter(scenario %in% 0:3, time <= 720) %>%
    ggplot(aes(x=time/24, y=NMN_pl, color=scenario_label)) +
    geom_line(size=1) +
    geom_hline(yintercept=0.9, linetype="dashed", color="red") +
    annotate("text", x=0.5, y=0.88, label="3× ULN (0.9 nmol/L)", size=3, hjust=0, color="red") +
    scale_color_manual(values=pal6[1:4]) +
    labs(title="PPGL QSP — Plasma Normetanephrine (Biomarker)",
         x="Time (days)", y="NMN (nmol/L)", color="Scenario") +
    theme_bw(base_size=12)
  print(p2)

  # Plot 3: Tumor volume — Sunitinib vs untreated
  p3 <- sim_df %>%
    filter(scenario %in% c(0, 4)) %>%
    ggplot(aes(x=time/24, y=Tumor_mL, color=scenario_label)) +
    geom_line(size=1.2) +
    scale_color_manual(values=c("No treatment"="#E53935","Sunitinib (malignant)"="#8E24AA")) +
    labs(title="PPGL QSP — Tumor Volume (Malignant PPGL)",
         x="Time (days)", y="Tumor Volume (mL)", color="Scenario") +
    theme_bw(base_size=12)
  print(p3)

  # Plot 4: TH activity under Metyrosine
  p4 <- sim_df %>%
    filter(scenario %in% c(0, 3, 5)) %>%
    ggplot(aes(x=time/24, y=TH_activity, color=scenario_label)) +
    geom_line(size=1) +
    scale_color_manual(values=pal6[c(1,3,6)]) +
    labs(title="PPGL QSP — TH Activity Under Metyrosine",
         x="Time (days)", y="TH Activity (normalized)", color="Scenario") +
    theme_bw(base_size=12)
  print(p4)
}

cat("\nModel complete. Six scenarios simulated successfully.\n")
cat("Key files:\n  ppgl_qsp_model.dot/.svg/.png — mechanistic map\n")
cat("  ppgl_mrgsolve_model_refactored.R — this file\n")
cat("  ppgl_shiny_app.R — interactive Shiny dashboard\n")
cat("  ppgl_references.md — 45 PubMed citations\n")
