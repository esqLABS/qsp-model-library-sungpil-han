## ============================================================
## Takayasu Arteritis (TA) — QSP mrgsolve ODE Model — REFACTORED (pluggable PK)
## Author : QSP Disease Model Library (CCR) — refactor pass
## Date   : 2026-08-30
## ============================================================
## Key References:
##   Nakaoka 2018 (Lancet): tocilizumab in TA (TAKT trial)
##   Hellmich 2020: EULAR recommendations for large-vessel vasculitis
##   Tombetti 2019 (Nat Rev Rheum): pathogenesis review
##   Hatemi 2022: ITAS disease activity scoring
##   Mekinian 2012: PET-CT monitoring of TA
## ------------------------------------------------------------
## REFACTOR SCOPE (see ta_refactor_notes.md for full detail):
##   Only the four compounds' own PK blocks and their downstream Hill
##   effect terms were touched (renamed to the fork's pluggable-PK
##   naming convention: GUT_/CENT_/PERI_/C_/EFFECT_/EMAX_/EC50_/GAMMA_
##   per stem PRED, TCZ, MTX, IFX). Every disease-biology compartment
##   (IL6, sIL6R, IL6_cmplx, TNF, TH1, TH17, TREG, VWI, ST) and every
##   biomarker (CRP, PET, VWT) plus $TABLE's NIH_SCORE/ITAS_SCORE/
##   RESPONSE_FLAG are otherwise byte-identical in formula to the
##   original -- only the drug-effect variable names feeding into them
##   were substituted 1:1 (Inh_PRED->EFFECT_PRED, Occ_TCZ->EFFECT_TCZ,
##   Inh_MTX->EFFECT_MTX, Inh_IFX->EFFECT_IFX).
##
##   BUILD-COMPATIBILITY FIX (disclosed, non-numeric; see notes and
##   UPSTREAM_ISSUES.md): the original does not compile under mrgsolve
##   2.0.1 at all (three independent, pre-existing defects, all outside
##   any compound's own PK/PD math). The same three syntax-only fixes
##   verification needed are carried into this delivered file:
##     1. $PARAM and $CMT (both non-annotated blocks) used C-style
##        block comments, which this mrgsolve build's non-annotated
##        block parser cannot handle -- stripped (comments only, no
##        value/name touched).
##     2. $INIT (non-annotated `name = value` block) redeclared the same
##        compartment names $CMT already declares ("Duplicated model
##        names") -- moved to the standard $MAIN `<CMT>_0 = value;` idiom.
##     3. The original's three inline `capture A B C` statements inside
##        $TABLE are not valid multi-name capture syntax under this
##        mrgsolve build (C++ compile error) -- consolidated into one
##        ordinary $CAPTURE block, same names, same values.
##   None of these three changes the value of any parameter, compartment,
##   or captured output -- confirmed by verification (see notes).
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ------------------------------------------------------------
## 1. Model Definition
## ------------------------------------------------------------
ta_model_code <- '
$PROB
  Takayasu Arteritis QSP Model (refactored: pluggable PK naming)
  20 ODE compartments: Drug PK (Prednisone/Prednisolone, Tocilizumab,
  Methotrexate, Infliximab) + Disease Biology (IL-6, sIL-6R, TNF-alpha,
  Th1, Th17, Treg, Vessel Wall Inflammation, Stenosis Index) +
  Biomarkers (CRP, ESR, PET-SUV, MRI-VWT)

$PARAM
  
  KA_PRED  = 1.50   
  V1_PRED  = 28.0   
  V2_PRED  = 56.0   
  Q_PRED   = 40.0   
  CL_PRED  = 18.0   
  KE0_PRED = 0.25   
  F_PRED   = 0.82   

  
  KA_TCZ   = 0.0083 
  V1_TCZ   = 3.5    
  V2_TCZ   = 2.0    
  Q_TCZ    = 0.30   
  CL_TCZ   = 0.18   
  CLMM_TCZ = 0.080  
  KM_TCZ   = 1.20   
  F_TCZ    = 0.80   

  
  KA_MTX     = 0.90   
  V1_MTX     = 18.0   
  CL_MTX     = 5.0    
  KPOLY_MTX  = 0.12   
  KDEPOLY_MTX= 0.018  
  F_MTX      = 0.70

  
  V1_IFX   = 3.0    
  V2_IFX   = 1.8    
  Q_IFX    = 0.22   
  CL_IFX   = 0.16   

  
  ksyn_IL6 = 0.30   
  kdeg_IL6 = 0.25   
  IL6_base = 1.2    

  
  ksyn_sR  = 3.0    
  kdeg_sR  = 0.015  
  sR_base  = 200.0  
  kon_IL6R = 0.001  
  koff_IL6R= 0.002  

  
  ksyn_TNF = 0.10
  kdeg_TNF = 0.30
  TNF_base = 0.5    

  
  
  ksyn_Th1 = 0.05   
  kdeg_Th1 = 0.02
  Th1_base = 25.0   

  
  ksyn_Th17= 0.03
  kdeg_Th17= 0.02
  Th17_base= 15.0

  
  ksyn_Treg= 0.02
  kdeg_Treg= 0.018
  Treg_base= 8.0

  
  ksyn_VWI = 0.008  
  kdeg_VWI = 0.004
  VWI_base = 2.0    

  
  kprog_ST = 0.0005 
  kreg_ST  = 0.0001 
  ST_base  = 0.0    

  
  ksyn_CRP = 0.60   
  kdeg_CRP = 0.025  
  CRP_base = 5.0    

  
  kESR     = 2.5    
  ESR_base = 20.0   

  
  ksyn_PET = 0.0015
  kdeg_PET = 0.003
  PET_base = 1.5    

  
  kVWT     = 0.002  
  VWT_base = 1.5    

  
  EMAX_PRED   = 0.85   
  EC50_PRED   = 0.10   
  GAMMA_PRED  = 1.5

  EMAX_TCZ    = 0.95   
  EC50_TCZ    = 0.50   
  GAMMA_TCZ   = 1.2

  EMAX_MTX    = 0.55   
  EC50_MTX    = 0.05   
  GAMMA_MTX   = 1.0

  EMAX_IFX    = 0.90   
  EC50_IFX    = 0.80   
  GAMMA_IFX   = 1.3

  
  amp_IL6_Th1  = 0.006  
  amp_IL6_Th17 = 0.004  
  amp_VWI_IL6  = 0.20   
  amp_VWI_TNF  = 0.15   
  amp_VWI_Th17 = 0.08   
  amp_ST_VWI   = 0.40   
  amp_IL6_TNF  = 0.05   
  inh_Treg_Th1 = 0.025  
  inh_Treg_Th17= 0.030  

  
  WT    = 65.0   
  DOSE_PRED = 0  
  DOSE_TCZ  = 0  
  DOSE_MTX  = 0  
  DOSE_IFX  = 0  

$CMT
  
  GUT_PRED CENT_PRED PERI_PRED EFF_PRED
  
  GUT_TCZ CENT_TCZ PERI_TCZ
  
  GUT_MTX CENT_MTX POLY_MTX
  
  CENT_IFX PERI_IFX
  
  IL6 sIL6R IL6_cmplx TNF
  TH1 TH17 TREG
  VWI ST
  
  CRP PET VWT

$MAIN
  GUT_PRED_0  = 0;   CENT_PRED_0 = 0;   PERI_PRED_0 = 0;   EFF_PRED_0 = 0;
  GUT_TCZ_0   = 0;   CENT_TCZ_0  = 0;   PERI_TCZ_0  = 0;
  GUT_MTX_0   = 0;   CENT_MTX_0  = 0;   POLY_MTX_0  = 0;
  CENT_IFX_0  = 0;   PERI_IFX_0  = 0;
  IL6_0       = 1.2; sIL6R_0     = 200; IL6_cmplx_0 = 0; TNF_0 = 0.5;
  TH1_0       = 25;  TH17_0      = 15;  TREG_0 = 8;
  VWI_0       = 2.0; ST_0        = 0.0;
  CRP_0       = 5.0; PET_0       = 1.5; VWT_0 = 1.5;

$ODE
  /* ========== Drug PK ========== */

  /* -- Prednisone gut -> prednisolone plasma (Archetype 3) -- */
  double dGUT_PRED  = -KA_PRED * GUT_PRED;
  double dCENT_PRED = KA_PRED * F_PRED * GUT_PRED / V1_PRED
                     - (CL_PRED/V1_PRED) * CENT_PRED
                     - (Q_PRED/V1_PRED) * CENT_PRED
                     + (Q_PRED/V2_PRED) * PERI_PRED;
  double dPERI_PRED = (Q_PRED/V1_PRED) * CENT_PRED
                     - (Q_PRED/V2_PRED) * PERI_PRED;
  double dEFF_PRED  = KE0_PRED * (CENT_PRED - EFF_PRED);

  /* -- Tocilizumab SC -> plasma (Archetype 3 + parallel linear/MM clearance) -- */
  double dGUT_TCZ = -KA_TCZ * GUT_TCZ;
  double CL_TCZ_tot = CL_TCZ + CLMM_TCZ * CENT_TCZ / (KM_TCZ + CENT_TCZ);
  double dCENT_TCZ  = KA_TCZ * F_TCZ * GUT_TCZ / V1_TCZ
                     - (CL_TCZ_tot/V1_TCZ) * CENT_TCZ
                     - (Q_TCZ/V1_TCZ) * CENT_TCZ
                     + (Q_TCZ/V2_TCZ) * PERI_TCZ;
  double dPERI_TCZ  = (Q_TCZ/V1_TCZ) * CENT_TCZ
                     - (Q_TCZ/V2_TCZ) * PERI_TCZ;

  /* -- MTX oral -> plasma -> polyglutamates (Archetype 3 minus peripheral, plus POLY_MTX) -- */
  double dGUT_MTX  = -KA_MTX * GUT_MTX;
  double dCENT_MTX = KA_MTX * F_MTX * GUT_MTX / V1_MTX
                    - (CL_MTX/V1_MTX) * CENT_MTX
                    - KPOLY_MTX * CENT_MTX;
  double dPOLY_MTX = KPOLY_MTX * CENT_MTX - KDEPOLY_MTX * POLY_MTX;

  /* -- Infliximab IV -> plasma (Archetype 2, no depot) -- */
  double dCENT_IFX = -(CL_IFX/V1_IFX) * CENT_IFX
                    - (Q_IFX/V1_IFX) * CENT_IFX
                    + (Q_IFX/V2_IFX) * PERI_IFX;
  double dPERI_IFX = (Q_IFX/V1_IFX) * CENT_IFX
                    - (Q_IFX/V2_IFX) * PERI_IFX;

  /* ========== Exposed concentrations (pluggable-PK interface) ========== */
  /* Each C_<STEM> is exactly the value the original’s own effect term read --
     no volume division was ever applied at the point of use in the original,
     for any of the four compounds; preserved verbatim, not "fixed" (see notes).
     fmax(..., 0.0) guards a NaN fragility found while verifying (logged in
     UPSTREAM_ISSUES.md and ta_refactor_notes.md): every one of these
     compartments can numerically undershoot to a tiny negative value between
     doses during long adaptive-step integration, and pow() on a negative base
     with the non-integer GAMMA_PRED/GAMMA_TCZ/GAMMA_IFX exponents below is NaN
     in C++, which then propagates through the whole coupled ODE system for
     the rest of the run. This is inside these four compounds’ own scope, so
     fixed here per the guide’s point 4, not merely logged -- confirmed by
     verification to change nothing except averting that blowup (see notes). */
  double C_PRED = fmax(EFF_PRED, 0.0);   /* effect-site (biophase) level -- the value Inh_PRED actually read, not raw plasma CENT_PRED */
  double C_TCZ  = fmax(CENT_TCZ, 0.0);   /* compartment read directly, undivided -- matches original Occ_TCZ */
  double C_MTX  = fmax(POLY_MTX, 0.0);   /* active polyglutamate pool -- matches original Inh_MTX; GAMMA_MTX=1 so not NaN-exposed, guarded for consistency */
  double C_IFX  = fmax(CENT_IFX, 0.0);   /* compartment read directly, undivided -- matches original Inh_IFX */

  /* ========== Drug PD (Hill interface: rename, not a fit) ========== */

  /* Prednisolone inhibits IL-6, TNF, Th1, Th17 via GR */
  double EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) /
                    (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));

  /* Tocilizumab blocks IL-6R (trans-signaling & signaling) */
  double EFFECT_TCZ  = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ) /
                    (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));

  /* MTX-PG inhibits T cell proliferation */
  double EFFECT_MTX  = EMAX_MTX * pow(C_MTX, GAMMA_MTX) /
                    (pow(EC50_MTX, GAMMA_MTX) + pow(C_MTX, GAMMA_MTX));

  /* Infliximab neutralizes TNF-alpha */
  double EFFECT_IFX  = EMAX_IFX * pow(C_IFX, GAMMA_IFX) /
                    (pow(EC50_IFX, GAMMA_IFX) + pow(C_IFX, GAMMA_IFX));

  /* ========== Disease Compartments ========== */

  /* IL-6 ODE
     Synthesis driven by Th1, Th17, TNF, amplified by disease;
     Inhibited by prednisolone (GR) and blocked by TCZ (feedback rise sIL-6R)
     TCZ leads to paradoxical serum IL-6 rise (sIL-6R)
  */
  double IL6_syn = ksyn_IL6
                   + amp_IL6_Th1  * TH1
                   + amp_IL6_Th17 * TH17
                   + amp_IL6_TNF  * TNF;
  double IL6_deg = kdeg_IL6 * (1 + EFFECT_TCZ * 0.2);  /* slight increase in IL-6 half-life blockade */
  double dIL6    = IL6_syn * (1 - EFFECT_PRED)
                   - IL6_deg * IL6
                   - kon_IL6R * IL6 * sIL6R
                   + koff_IL6R * IL6_cmplx;

  /* Soluble IL-6R: TCZ increases sIL-6R as free receptor accumulates */
  double sIL6R_factor = 1 + 2.5 * EFFECT_TCZ;  /* TCZ blocks membrane IL-6R, shedding increases */
  double dSIL6R  = ksyn_sR * sIL6R_factor
                   - kdeg_sR * sIL6R
                   - kon_IL6R * IL6 * sIL6R
                   + koff_IL6R * IL6_cmplx;

  double dIL6_cmplx = kon_IL6R * IL6 * sIL6R
                      - koff_IL6R * IL6_cmplx
                      - kdeg_IL6 * IL6_cmplx;

  /* TNF-alpha ODE */
  double dTNF = ksyn_TNF + 0.03 * TH1
                - kdeg_TNF * TNF * (1 - EFFECT_IFX) * (1 - 0.6 * EFFECT_PRED);

  /* Th1 ODE: driven by IL-6, IFN-gamma circuit; inhibited by Treg and drugs */
  double dTH1 = ksyn_Th1 * (1 + 0.05 * IL6)
                - kdeg_Th1 * TH1 * (1 - EFFECT_MTX) * (1 - EFFECT_PRED * 0.5)
                - inh_Treg_Th1 * TREG * TH1;

  /* Th17 ODE: driven by IL-6 + IL-23 axis, TCZ markedly suppresses */
  double dTH17 = ksyn_Th17 * (1 + 0.04 * IL6) * (1 - EFFECT_TCZ * 0.8)
                 - kdeg_Th17 * TH17 * (1 - EFFECT_MTX * 0.7) * (1 - EFFECT_PRED * 0.4)
                 - inh_Treg_Th17 * TREG * TH17;

  /* Treg ODE: prednisolone and TCZ partially restore Treg */
  double dTREG = ksyn_Treg * (1 + 0.3 * EFFECT_PRED + 0.2 * EFFECT_TCZ)
                 - kdeg_Treg * TREG;

  /* Vessel Wall Inflammation Index [0-10]
     Driven by IL-6, TNF, Th17; suppressed by all drugs via Emax
  */
  double VWI_drive = amp_VWI_IL6 * IL6
                     + amp_VWI_TNF * TNF
                     + amp_VWI_Th17 * TH17;
  double Drug_inh_VWI = 1 - (1 - EFFECT_PRED) * (1 - EFFECT_TCZ * 0.9)
                              * (1 - EFFECT_IFX * 0.7) * (1 - EFFECT_MTX * 0.3);
  double dVWI = ksyn_VWI * VWI_drive * (1 - Drug_inh_VWI)
                - kdeg_VWI * VWI;

  /* Stenosis Index [0-100 %] — driven by cumulative VWI */
  double dST = kprog_ST * amp_ST_VWI * VWI
               - kreg_ST * ST;

  /* CRP (mg/L) — produced proportional to IL-6; TCZ rapidly normalizes */
  double dCRP = ksyn_CRP * IL6 * (1 - EFFECT_TCZ * 0.95)
                - kdeg_CRP * CRP;

  /* PET-CT SUVmax — correlates with VWI and vessel wall metabolic activity */
  double dPET = ksyn_PET * VWI
                - kdeg_PET * PET;

  /* MRI Vessel Wall Thickness (mm) — driven by cumulative stenosis/VWI */
  double dVWT = kVWT * VWI - 0.001 * VWT;

  /* ========== Assign DES ========== */
  dxdt_GUT_PRED  = dGUT_PRED;
  dxdt_CENT_PRED = dCENT_PRED;
  dxdt_PERI_PRED = dPERI_PRED;
  dxdt_EFF_PRED  = dEFF_PRED;
  dxdt_GUT_TCZ   = dGUT_TCZ;
  dxdt_CENT_TCZ  = dCENT_TCZ;
  dxdt_PERI_TCZ  = dPERI_TCZ;
  dxdt_GUT_MTX   = dGUT_MTX;
  dxdt_CENT_MTX  = dCENT_MTX;
  dxdt_POLY_MTX  = dPOLY_MTX;
  dxdt_CENT_IFX  = dCENT_IFX;
  dxdt_PERI_IFX  = dPERI_IFX;
  dxdt_IL6       = dIL6;
  dxdt_sIL6R     = dSIL6R;
  dxdt_IL6_cmplx = dIL6_cmplx;
  dxdt_TNF       = dTNF;
  dxdt_TH1       = dTH1;
  dxdt_TH17      = dTH17;
  dxdt_TREG      = dTREG;
  dxdt_VWI       = dVWI;
  dxdt_ST        = dST;
  dxdt_CRP       = dCRP;
  dxdt_PET       = dPET;
  dxdt_VWT       = dVWT;

$TABLE
  double NIH_SCORE = 0;
  /* NIH Disease Activity Score: 0-20 scale
     based on new/worsening features:
     systemic symptoms (2pts), ESR (2pts), angiography (2pts),
     ischemic symptoms (6pts), BP difference (2pts)
  */
  double ESR_now = ESR_base + kESR * (CRP - CRP_base);
  NIH_SCORE = 2 * (CRP > 20 ? 1 : CRP/20) +  /* systemic inflammation */
              2 * (ESR_now > 40 ? 1 : ESR_now/40) +  /* ESR */
              3 * (VWI/10) +                   /* vascular inflammation */
              4 * (ST/50) +                    /* stenosis extent */
              3 * (PET/4);                     /* PET activity */
  if(NIH_SCORE > 20) NIH_SCORE = 20;

  double ITAS_SCORE = 0;
  /* ITAS 2010 simplification */
  ITAS_SCORE = 1.5 * (CRP > 10 ? 1 : 0) +
               1.5 * (VWI > 5 ? 1 : 0) +
               3.0 * (ST > 20 ? 1 : 0) +
               2.0 * (PET > 2.5 ? 1 : 0);

  double RESPONSE_FLAG = NIH_SCORE < 4 ? 1 : 0;  /* low disease activity */
  double CP_PRED  = CENT_PRED;
  double CP_TCZ   = CENT_TCZ;
  double CP_MTX   = CENT_MTX;
  double CP_IFX   = CENT_IFX;

$CAPTURE ESR_now NIH_SCORE ITAS_SCORE RESPONSE_FLAG CP_PRED CP_TCZ CP_MTX CP_IFX Drug_inh_VWI C_PRED C_TCZ C_MTX C_IFX EFFECT_PRED EFFECT_TCZ EFFECT_MTX EFFECT_IFX
'

mod <- mcode("TakayasuArteritis", ta_model_code)

## ------------------------------------------------------------
## 2. Helper Functions: Dosing Events
## ------------------------------------------------------------

## Prednisone 1 mg/kg/day orally (continuous via events)
pred_events <- function(dose_mg = 65, duration_days = 365,
                        taper_to = 10, taper_start = 60) {
  times_init  <- seq(0, taper_start * 24, by = 24)
  times_taper <- seq((taper_start + 7) * 24, duration_days * 24, by = 24)
  dose_taper  <- seq(dose_mg, taper_to,
                     length.out = length(times_taper))
  ev_init  <- ev(time = times_init,  amt = dose_mg, cmt = "GUT_PRED")
  ev_taper <- ev(time = times_taper, amt = dose_taper, cmt = "GUT_PRED")
  c(ev_init, ev_taper)
}

## Tocilizumab 162 mg SC q2w (standard TA dosing)
tcz_events <- function(start_day = 0, n_doses = 26) {
  times <- seq(start_day * 24, by = 14 * 24, length.out = n_doses)
  ev(time = times, amt = 162, cmt = "GUT_TCZ")
}

## Methotrexate 15 mg/week orally
mtx_events <- function(start_day = 0, duration_days = 365) {
  times <- seq(start_day * 24, duration_days * 24, by = 7 * 24)
  ev(time = times, amt = 15, cmt = "GUT_MTX")
}

## Infliximab 5 mg/kg IV: 0, 2, 6 weeks then q6w
ifx_events <- function(start_day = 0, wt_kg = 65, n_maint = 8) {
  dose <- 5 * wt_kg
  induction_times <- c(0, 2, 6) * 7 * 24 + start_day * 24
  maint_times     <- seq((6 + 6) * 7 * 24, by = 6 * 7 * 24,
                         length.out = n_maint) + start_day * 24
  ev(time = c(induction_times, maint_times),
     amt  = dose, cmt = "CENT_IFX", rate = -2)  # 2h infusion
}

## ------------------------------------------------------------
## 3. Simulation Parameters
## ------------------------------------------------------------
sim_end <- 365 * 24   # 1 year in hours
dt      <- 4           # 4-hour output step
times   <- seq(0, sim_end, by = dt)

base_params <- list(
  ksyn_IL6 = 0.45,   # active TA: elevated IL-6 synthesis
  ksyn_TNF = 0.20,
  ksyn_Th1 = 0.10,
  ksyn_Th17= 0.07,
  VWI      = 6.0,    # active disease starting VWI
  IL6      = 8.0,    # elevated baseline
  TNF      = 2.0,
  TH1      = 50.0,
  TH17     = 35.0,
  CRP      = 45.0,
  PET      = 3.5,
  VWT      = 4.5
)

## ------------------------------------------------------------
## 4. Treatment Scenarios
## ------------------------------------------------------------

## Scenario 1: No treatment (natural history)
ev_none <- ev(time = 0, amt = 0, cmt = "GUT_PRED")

## Scenario 2: Prednisone monotherapy (1 mg/kg/day → taper to 10 mg/day)
ev_pred <- pred_events(dose_mg = 65, duration_days = 365,
                       taper_to = 10, taper_start = 60)

## Scenario 3: Prednisone + Methotrexate (standard first-line combo)
ev_pred_mtx <- c(
  pred_events(dose_mg = 65, taper_to = 10, taper_start = 60),
  mtx_events(start_day = 0)
)

## Scenario 4: Prednisone + Tocilizumab (TAKT trial regimen)
ev_pred_tcz <- c(
  pred_events(dose_mg = 65, taper_to = 7.5, taper_start = 90),
  tcz_events(start_day = 0)
)

## Scenario 5: Prednisone + Infliximab (refractory TA)
ev_pred_ifx <- c(
  pred_events(dose_mg = 65, taper_to = 10, taper_start = 60),
  ifx_events(start_day = 0)
)

scenarios <- list(
  "1. No Treatment (Natural History)" = ev_none,
  "2. Prednisone Monotherapy"         = ev_pred,
  "3. Prednisone + Methotrexate"      = ev_pred_mtx,
  "4. Prednisone + Tocilizumab (TAKT)"= ev_pred_tcz,
  "5. Prednisone + Infliximab (Refract.)"= ev_pred_ifx
)

## ------------------------------------------------------------
## 5. Run All Scenarios
## ------------------------------------------------------------
run_scenario <- function(ev_obj, scen_name) {
  idata <- as.data.frame(base_params)
  out <- mod %>%
    param(base_params) %>%
    init(IL6 = base_params$IL6,
         TNF = base_params$TNF,
         TH1 = base_params$TH1,
         TH17 = base_params$TH17,
         VWI = base_params$VWI,
         CRP = base_params$CRP,
         PET = base_params$PET,
         VWT = base_params$VWT) %>%
    ev(ev_obj) %>%
    mrgsim(end = sim_end, delta = dt, carry_out = "evid") %>%
    as_tibble() %>%
    mutate(time_days = time / 24, Scenario = scen_name)
  out
}

results <- bind_rows(lapply(names(scenarios), function(nm) {
  run_scenario(scenarios[[nm]], nm)
}))

## ------------------------------------------------------------
## 6. Visualization
## ------------------------------------------------------------
theme_ta <- theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#E65100", color = "white"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

scen_colors <- c(
  "1. No Treatment (Natural History)"       = "#D32F2F",
  "2. Prednisone Monotherapy"               = "#F57C00",
  "3. Prednisone + Methotrexate"            = "#388E3C",
  "4. Prednisone + Tocilizumab (TAKT)"      = "#1565C0",
  "5. Prednisone + Infliximab (Refract.)"   = "#6A1B9A"
)

p1 <- ggplot(results, aes(time_days, IL6, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "A. Serum IL-6 (pg/mL)", x = "Time (days)", y = "IL-6 (pg/mL)") +
  theme_ta

p2 <- ggplot(results, aes(time_days, CRP, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 12, label = "CRP = 10 mg/L", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "B. CRP (mg/L)", x = "Time (days)", y = "CRP (mg/L)") +
  theme_ta

p3 <- ggplot(results, aes(time_days, NIH_SCORE, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 4.5, label = "Remission threshold", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "C. NIH Disease Activity Score (0-20)", x = "Time (days)", y = "NIH Score") +
  theme_ta

p4 <- ggplot(results, aes(time_days, VWI, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "D. Vessel Wall Inflammation Index (0-10)", x = "Time (days)", y = "VWI Score") +
  theme_ta

p5 <- ggplot(results, aes(time_days, ST, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "E. Arterial Stenosis Index (%)", x = "Time (days)", y = "Stenosis (%)") +
  theme_ta

p6 <- ggplot(results, aes(time_days, PET, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 330, y = 2.7, label = "PET activity threshold", size = 3) +
  scale_color_manual(values = scen_colors) +
  labs(title = "F. PET-CT FDG SUVmax (Vascular)", x = "Time (days)", y = "SUVmax") +
  theme_ta

p_pk1 <- ggplot(
  results %>% filter(grepl("TCZ", Scenario)),
  aes(time_days, CP_TCZ, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "G. Tocilizumab Plasma Concentration",
       x = "Time (days)", y = "Concentration (mg/L)") +
  theme_ta

p_th <- ggplot(
  results %>% select(time_days, Scenario, TH1, TH17, TREG) %>%
    pivot_longer(c(TH1, TH17, TREG), names_to = "Cell", values_to = "Count"),
  aes(time_days, Count, color = Cell, linetype = Scenario)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = c(TH1 = "#D32F2F", TH17 = "#1565C0", TREG = "#2E7D32"),
                     labels = c("Th1 Cells", "Th17 Cells", "Treg Cells")) +
  labs(title = "H. T Cell Populations (cells/µL)", x = "Time (days)", y = "Count (cells/µL)") +
  theme_ta

## Combined dashboard
dashboard <- (p1 + p2) / (p3 + p4) / (p5 + p6)
dashboard_full <- dashboard + plot_annotation(
  title = "Takayasu Arteritis QSP Model — Treatment Scenario Comparison",
  subtitle = "5 scenarios: natural history vs. prednisolone ± MTX / TCZ / IFX",
  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
)

## ------------------------------------------------------------
## 7. Clinical Trial Calibration Reference
## ------------------------------------------------------------
calibration_table <- data.frame(
  Trial          = c("TAKT (Nakaoka 2018 Lancet)",
                     "NIHON-BVAS Pred monotherapy",
                     "Abisror 2013 JRA",
                     "Comarmond 2012 Medicine",
                     "Hoffman 1994 Ann Intern Med"),
  Intervention   = c("Tocilizumab + Pred vs Pred",
                     "Pred 1 mg/kg taper",
                     "Infliximab (TNF)",
                     "Infliximab rescue",
                     "Methotrexate + Pred"),
  Key_Outcome    = c("Time to relapse HR 0.41 (TCZ arm)",
                     "~70% initial remission rate",
                     "Remission 93% refractory TA",
                     "Response 67% refractory",
                     "Steroid-sparing; 72% remission"),
  Model_Parameter= c("EMAX_TCZ=0.95, EC50_TCZ=0.50",
                     "EMAX_PRED=0.85, taper over 60 days",
                     "EMAX_IFX=0.90, EC50_IFX=0.80",
                     "EFFECT_IFX 0-0.90 dose-range",
                     "EMAX_MTX=0.55, EC50_MTX=0.05"),
  stringsAsFactors = FALSE
)

cat("==================================================\n")
cat("  Takayasu Arteritis QSP Model — Calibration\n")
cat("==================================================\n")
print(calibration_table, row.names = FALSE)

cat("\n--- End-of-year summary (day 365) ---\n")
summary_365 <- results %>%
  filter(abs(time_days - 365) < 0.25) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  select(Scenario, IL6, CRP, VWI, ST, PET, NIH_SCORE, ITAS_SCORE) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))
print(as.data.frame(summary_365), row.names = FALSE)

print(dashboard_full)

## ------------------------------------------------------------
## 8. Sensitivity Analysis: IL-6 synthesis rate vs. TCZ efficacy
## ------------------------------------------------------------
cat("\n--- Sensitivity: ksyn_IL6 vs. 1-year NIH score (TCZ scenario) ---\n")
ksyn_range <- seq(0.20, 0.80, by = 0.15)
sens_results <- lapply(ksyn_range, function(k) {
  p_override <- modifyList(base_params, list(ksyn_IL6 = k, IL6 = k / 0.25 * 1.2))
  out <- mod %>%
    param(p_override) %>%
    init(IL6 = p_override$IL6, CRP = 45, VWI = 6, TH1 = 50, TH17 = 35) %>%
    ev(scenarios[["4. Prednisone + Tocilizumab (TAKT)"]]) %>%
    mrgsim(end = sim_end, delta = dt) %>%
    as_tibble() %>%
    filter(abs(time / 24 - 365) < 0.5) %>%
    slice(1) %>%
    mutate(ksyn_IL6 = k)
})
sens_df <- bind_rows(sens_results)
cat(sprintf("  ksyn_IL6=%.2f → NIH_SCORE=%.2f\n",
            sens_df$ksyn_IL6, sens_df$NIH_SCORE))
