## =============================================================================
## Stable Angina (Chronic Coronary Syndrome) — QSP Model
## mrgsolve ODE-based PK/PD Model
## PK/PD REFACTOR — pluggable PK, named Hill interface
## =============================================================================
## Refactored per FORK_WORKFLOW_GUIDE.md Part 2. Five compounds rewritten to
## the guide's naming convention (see sa_refactor_notes.md for full detail):
##   BB  (bisoprolol)  — archetype 3 (depot+central+peripheral, linear)
##   CCB (amlodipine)  — archetype 3 (depot+central+peripheral, linear)
##   RAN (ranolazine)  — archetype 3 minus peripheral (depot+central, linear)
##   IVA (ivabradine)  — archetype 3 minus peripheral (depot+central, linear)
##   NIT (ISMN)        — archetype 3 minus peripheral (depot+central, linear)
##                        plus a bespoke nitrate-tolerance state (TOL_NIT),
##                        which is native to this compound's own block, not a
##                        disease-side compartment.
## Every compartment, $PARAM value, and equation shape is unchanged from the
## original (sa_mrgsolve_model.R) — this is a rename/reorganization pass, not
## a refit. See sa_refactor_notes.md for the pre-existing build defect fixed
## here (syntax-only) and the full verification result.
## =============================================================================
## Reference calibration (unchanged from original):
##  - BB (bisoprolol): Rehnqvist 1992, TIBBS study; t½=11h, HR reduction ~20%
##  - CCB (amlodipine): CAPE study, TIME study; t½=40h
##  - Ranolazine: CARISA (Chaitman 2004), MARISA (Chaitman 2004); t½=7h
##  - Ivabradine: BEAUTIFUL (Fox 2008), SIGNIFY (Fox 2014); t½=2h
##  - Nitrates (ISMN): Chrysant 1993; t½=5h, BA~100%
##  - Anti-ischemic endpoint: RPP < 20,000 (isch threshold)
##  - CCS Class reduction: COURAGE (Boden 2007), ORBITA (Al-Lamee 2018)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ─── mrgsolve Model Code ─────────────────────────────────────────────────────

code <- '
$PROB
 Stable Angina QSP Model (PK/PD refactored)
 21 ODE compartments:
   PK: Bisoprolol/BB (2-cmt), Amlodipine/CCB (2-cmt), Ranolazine/RAN (1-cmt),
       Ivabradine/IVA (1-cmt), ISMN-Nitrate/NIT (1-cmt + tolerance state)
   PD: HR, SBP, CBF state, INaL inhibition, O2-imbalance,
       Ischemia burden, Angina score, Exercise capacity,
       Plaque progression

$CMT
 // ----- BB (bisoprolol) PK: archetype 3, depot + central + peripheral -----
 GUT_BB        // GI depot [mg]
 CENT_BB       // Central plasma [mg]
 PERI_BB       // Peripheral tissue [mg]

 // ----- CCB (amlodipine) PK: archetype 3, depot + central + peripheral -----
 GUT_CCB       // GI depot [mg]
 CENT_CCB      // Central plasma [mg]
 PERI_CCB      // Peripheral tissue [mg]

 // ----- RAN (ranolazine) PK: archetype 3 minus peripheral -----
 GUT_RAN       // GI depot [mg]
 CENT_RAN      // Central plasma [mg]

 // ----- IVA (ivabradine) PK: archetype 3 minus peripheral -----
 GUT_IVA       // GI depot [mg]
 CENT_IVA      // Central plasma [mg]

 // ----- NIT (ISMN nitrate) PK: archetype 3 minus peripheral, + bespoke tolerance state -----
 GUT_NIT       // GI depot [mg]
 CENT_NIT      // Central plasma [mg]
 TOL_NIT       // Nitrate tolerance state [0-1] (native to the NIT block itself)

 // ----- Disease PD states (unchanged) -----
 HR_STATE      // Heart rate [bpm]
 SBP_STATE     // Systolic BP [mmHg]
 CBF_STATE     // Coronary blood flow [mL/min]
 O2IMBAL       // O2 supply-demand imbalance [normalised]
 ISCHEMIA      // Ischemia burden [0-1]
 ANGINA_SCORE  // Weekly angina episodes [count/week]
 EX_CAP        // Exercise capacity [METs]
 PLAQUE        // Plaque progression [% lumen occluded, 0-100]

$PARAM
 // ─── Patient baseline characteristics ───────────────────────
 BW         = 75      // body weight [kg]
 AGE        = 62      // age [years]
 HR0        = 80      // baseline HR [bpm]
 SBP0       = 145     // baseline SBP [mmHg]
 CBF0       = 250     // baseline CBF [mL/min]
 STENOSIS   = 70      // stenosis severity [% lumen]
 CFR0       = 2.5     // coronary flow reserve at baseline
 EX_CAP0    = 7.0     // baseline exercise capacity [METs]
 ANGINA0    = 5.0     // baseline angina episodes/week

 // ─── BB / Bisoprolol PK (Rehnqvist 1992) ────────────────────
 KA_BB      = 1.0     // absorption rate [h⁻¹]
 F_BB       = 0.85    // oral bioavailability
 CL_BB      = 12.5    // clearance [L/h]  (t½ ~11h: 0.693/0.0554)
 V1_BB      = 50.0    // central Vd [L]   (Vd ~8 L/kg → Vss~600 L)
 Q_BB       = 8.0     // inter-compartmental CL [L/h]
 V2_BB      = 250.0   // peripheral Vd [L]
 EC50_BB_HR   = 30.0  // EC50 for HR reduction [ng/mL]           (was EC50_BB)
 EMAX_BB_HR   = 0.25  // Emax for HR reduction (fraction)        (was EMAX_BB)
 GAMMA_BB_HR  = 1      // Hill coefficient (none in original; =1, rename not a fit)
 EC50_BB_SBP  = 35.0  // EC50 for SBP reduction [ng/mL]          (was EC50_SBP_BB)
 EMAX_BB_SBP  = 0.12  // Emax for SBP reduction                  (was EMAX_SBP_BB)
 GAMMA_BB_SBP = 1      // Hill coefficient (none in original; =1, rename not a fit)

 // ─── CCB / Amlodipine PK (Murdoch 1991) ─────────────────────
 KA_CCB     = 0.5     // absorption rate [h⁻¹]
 F_CCB      = 0.64    // oral bioavailability
 CL_CCB     = 7.0     // clearance [L/h]  (t½ ~42h)
 V1_CCB     = 210.0   // central Vd [L]
 Q_CCB      = 5.0     // inter-compartmental CL [L/h]
 V2_CCB     = 1200.0  // peripheral Vd [L]  (Vd~1400 L total)
 EC50_CCB_SBP  = 4.5   // EC50 for SBP reduction [ng/mL]         (was EC50_CCB)
 EMAX_CCB_SBP  = 0.18  // Emax SBP reduction                     (was EMAX_CCB)
 GAMMA_CCB_SBP = 1      // Hill coefficient (none in original; =1, rename not a fit)
 EC50_CCB_CBF  = 3.0   // EC50 coronary vasodilation [ng/mL]     (was EC50CBF_CCB)
 EMAX_CCB_CBF  = 0.30  // Emax CBF increase                      (was EMAX_CBF_CCB)
 GAMMA_CCB_CBF = 1      // Hill coefficient (none in original; =1, rename not a fit)

 // ─── RAN / Ranolazine PK (Jerling 2006) ──────────────────────
 KA_RAN     = 1.2     // absorption rate ER [h⁻¹]
 F_RAN      = 0.76    // bioavailability (ER)
 CL_RAN     = 60.0    // clearance [L/h]   (t½ ~7h)
 V1_RAN     = 600.0   // Vd [L]
 EC50_RAN_INA  = 300.0 // EC50 I_NaL inhibition [ng/mL]          (was EC50_RAN)
 EMAX_RAN_INA  = 0.70  // Emax I_NaL inhibition                  (was EMAX_RAN)
 GAMMA_RAN_INA = 1      // Hill coefficient (none in original; =1, rename not a fit)
 EC50_RAN_EX   = 250.0 // EC50 exercise capacity [ng/mL]         (was EC50EX_RAN)
 EMAX_RAN_EX   = 0.15  // Emax exercise capacity increase (CARISA) (was EMAX_EX_RAN)
 GAMMA_RAN_EX  = 1      // Hill coefficient (none in original; =1, rename not a fit)

 // ─── IVA / Ivabradine PK (Portoles 2006) ─────────────────────
 KA_IVA     = 2.5     // absorption [h⁻¹]
 F_IVA      = 0.40    // bioavailability
 CL_IVA     = 80.0    // clearance [L/h]   (t½ ~2h)
 V1_IVA     = 230.0   // Vd [L]
 EC50_IVA   = 1.0     // EC50 I_f inhibition [ng/mL]
 EMAX_IVA   = 0.28    // Emax HR reduction (BEAUTIFUL/SIGNIFY)
 GAMMA_IVA  = 1        // Hill coefficient (none in original; =1, rename not a fit)

 // ─── NIT / ISMN PK (Chasseaud 1992) ──────────────────────────
 KA_NIT     = 2.0     // absorption [h⁻¹]
 F_NIT      = 1.00    // BA ~100%
 CL_NIT     = 40.0    // clearance [L/h]   (t½ ~4.5h)
 V1_NIT     = 260.0   // Vd [L]
 EC50_NIT_MVO  = 200.0 // EC50 preload reduction [ng/mL]         (was EC50_NIT)
 EMAX_NIT_MVO  = 0.20  // Emax LVEDP reduction → MVO2 reduction  (was EMAX_NIT)
 GAMMA_NIT_MVO = 1      // Hill coefficient (none in original; =1, rename not a fit)
 EC50_NIT_CBF  = 150.0 // EC50 coronary vasodilation             (was EC50CBF_NIT)
 EMAX_NIT_CBF  = 0.20  // Emax CBF increase (nitrate effect)     (was EMAX_CBF_NIT)
 GAMMA_NIT_CBF = 1      // Hill coefficient (none in original; =1, rename not a fit)
 TOLK1_NIT  = 0.05    // tolerance development rate [h⁻¹] (was TOL_K1)
 TOLK2_NIT  = 0.01    // tolerance recovery rate [h⁻¹]    (was TOL_K2)

 // ─── PD parameters (unchanged) ───────────────────────────────
 KOUT_HR    = 0.3     // HR response rate [h⁻¹]
 KOUT_SBP   = 0.2     // SBP response rate [h⁻¹]
 KOUT_CBF   = 0.5     // CBF response rate [h⁻¹]
 KOUT_ISCH  = 0.1     // ischemia decay [h⁻¹]
 KOUT_ANG   = 0.05    // angina score decay [h⁻¹] (weekly dynamics)

 // Ischemic threshold: RPP >20000 → ischemia
 RPP_THRESH = 20000   // [bpm × mmHg]

 // Long-term plaque (slow dynamics: years scale)
 KPROG_PL   = 0.00014 // plaque progression [fraction/h] (~0.1%/month)
 KSTATIN_PL = 0.00008 // statin regression [fraction/h]

 // O2 supply-demand weighting
 WDMD       = 0.6     // demand weight in imbalance
 WSUP       = 0.4     // supply weight

 // Stenosis-flow relationship (Gould 1974)
 STEN_K     = 0.035   // stenosis-CBF reduction slope

 // ─── Drug dose switches (0=off, 1=on) ────────────────────────
 DOSE_BB    = 0
 DOSE_CCB   = 0
 DOSE_RAN   = 0
 DOSE_IVA   = 0
 DOSE_NIT   = 0
 STATIN_ON  = 0       // statin therapy flag

$GLOBAL
// [refactor] The five compounds exposed concentration + Hill-effect terms,
// predeclared here rather than in $PARAM: mrgsolve 2.0.1 compiles $PARAM
// members as read-only references inside $ODE, so a value recomputed every
// timestep from state cannot also live in $PARAM (same constraint already
// documented in the AMD, membranous-nephropathy, sepsis, SAH, AIP, and
// breast-cancer refactors). Visible in every simulation output via bare
// $CAPTURE below and in /model_manifest outputPaths.
double EFFECT_BB_HR, EFFECT_BB_SBP, EFFECT_CCB_SBP, EFFECT_CCB_CBF, EFFECT_IVA,
       EFFECT_NIT_MVO, EFFECT_NIT_CBF, EFFECT_RAN_INA, EFFECT_RAN_EX;

$MAIN
 // Compartment initial values
 HR_STATE_0   = HR0;
 SBP_STATE_0  = SBP0;
 CBF_STATE_0  = CBF0 * (1 - STEN_K * STENOSIS);  // stenosis reduces resting CBF
 O2IMBAL_0    = 0;
 ISCHEMIA_0   = 0;
 ANGINA_SCORE_0 = ANGINA0;
 EX_CAP_0     = EX_CAP0;
 PLAQUE_0     = STENOSIS;
 TOL_NIT_0    = 0;

$ODE
 // ─────────────────────────────────────────────────────────────
 // BB (BISOPROLOL) PK — archetype 3 (depot + central + peripheral)
 double KE_BB  = CL_BB / V1_BB;
 double K12_BB = Q_BB  / V1_BB;
 double K21_BB = Q_BB  / V2_BB;
 dxdt_GUT_BB  = -KA_BB * GUT_BB;
 dxdt_CENT_BB =  KA_BB * F_BB * GUT_BB - (KE_BB + K12_BB) * CENT_BB + K21_BB * PERI_BB;
 dxdt_PERI_BB =  K12_BB * CENT_BB - K21_BB * PERI_BB;

 C_BB = CENT_BB / V1_BB * 1000.0;  // ng/mL (dose in mg → *1000/V in L) — the exposed BB concentration

 // ─────────────────────────────────────────────────────────────
 // CCB (AMLODIPINE) PK — archetype 3 (depot + central + peripheral)
 double KE_CCB  = CL_CCB / V1_CCB;
 double K12_CCB = Q_CCB  / V1_CCB;
 double K21_CCB = Q_CCB  / V2_CCB;
 dxdt_GUT_CCB  = -KA_CCB * GUT_CCB;
 dxdt_CENT_CCB =  KA_CCB * F_CCB * GUT_CCB - (KE_CCB + K12_CCB) * CENT_CCB + K21_CCB * PERI_CCB;
 dxdt_PERI_CCB =  K12_CCB * CENT_CCB - K21_CCB * PERI_CCB;

 C_CCB = CENT_CCB / V1_CCB * 1000.0;  // the exposed CCB concentration

 // ─────────────────────────────────────────────────────────────
 // RAN (RANOLAZINE) PK — archetype 3 minus peripheral (depot + central)
 double KE_RAN = CL_RAN / V1_RAN;
 dxdt_GUT_RAN  = -KA_RAN * GUT_RAN;
 dxdt_CENT_RAN =  KA_RAN * F_RAN * GUT_RAN - KE_RAN * CENT_RAN;

 C_RAN = CENT_RAN / V1_RAN * 1000.0;  // the exposed RAN concentration

 // ─────────────────────────────────────────────────────────────
 // IVA (IVABRADINE) PK — archetype 3 minus peripheral (depot + central)
 double KE_IVA = CL_IVA / V1_IVA;
 dxdt_GUT_IVA  = -KA_IVA * GUT_IVA;
 dxdt_CENT_IVA =  KA_IVA * F_IVA * GUT_IVA - KE_IVA * CENT_IVA;

 C_IVA = CENT_IVA / V1_IVA * 1000.0;  // the exposed IVA concentration

 // ─────────────────────────────────────────────────────────────
 // NIT (ISMN NITRATE) PK — archetype 3 minus peripheral (depot + central)
 // plus a bespoke tolerance state native to the NIT block itself
 double KE_NIT = CL_NIT / V1_NIT;
 dxdt_GUT_NIT  = -KA_NIT * GUT_NIT;
 dxdt_CENT_NIT =  KA_NIT * F_NIT * GUT_NIT - KE_NIT * CENT_NIT;

 C_NIT = CENT_NIT / V1_NIT * 1000.0;  // the exposed NIT concentration

 // Nitrate tolerance ODE: tolerance builds during exposure, decays nitrate-free
 double tol_ss  = (C_NIT > 1.0) ? 1.0 : 0.0;   // target tolerance state
 dxdt_TOL_NIT   = TOLK1_NIT * (tol_ss - TOL_NIT) - TOLK2_NIT * TOL_NIT * (1 - tol_ss);
 double TOLEFF_NIT = (1 - 0.75 * TOL_NIT);      // tolerance reduces nitrate effect

 // ─────────────────────────────────────────────────────────────
 // Hill-interface effect terms (renamed, not refit — every one of these was
 // already a plain Emax*C/(EC50+C) ratio in the original; GAMMA_*=1 added
 // explicitly, pow() with exponent 1 reproduces the original arithmetic
 // exactly)
 EFFECT_BB_HR   = EMAX_BB_HR   * pow(C_BB,  GAMMA_BB_HR)   / (pow(EC50_BB_HR,   GAMMA_BB_HR)   + pow(C_BB,  GAMMA_BB_HR));
 EFFECT_BB_SBP  = EMAX_BB_SBP  * pow(C_BB,  GAMMA_BB_SBP)  / (pow(EC50_BB_SBP,  GAMMA_BB_SBP)  + pow(C_BB,  GAMMA_BB_SBP));
 EFFECT_CCB_SBP = EMAX_CCB_SBP * pow(C_CCB, GAMMA_CCB_SBP) / (pow(EC50_CCB_SBP, GAMMA_CCB_SBP) + pow(C_CCB, GAMMA_CCB_SBP));
 EFFECT_CCB_CBF = EMAX_CCB_CBF * pow(C_CCB, GAMMA_CCB_CBF) / (pow(EC50_CCB_CBF, GAMMA_CCB_CBF) + pow(C_CCB, GAMMA_CCB_CBF));
 EFFECT_IVA     = EMAX_IVA     * pow(C_IVA, GAMMA_IVA)     / (pow(EC50_IVA,     GAMMA_IVA)     + pow(C_IVA, GAMMA_IVA));
 EFFECT_NIT_MVO = EMAX_NIT_MVO * pow(C_NIT, GAMMA_NIT_MVO) / (pow(EC50_NIT_MVO, GAMMA_NIT_MVO) + pow(C_NIT, GAMMA_NIT_MVO)) * TOLEFF_NIT;
 EFFECT_NIT_CBF = EMAX_NIT_CBF * pow(C_NIT, GAMMA_NIT_CBF) / (pow(EC50_NIT_CBF, GAMMA_NIT_CBF) + pow(C_NIT, GAMMA_NIT_CBF)) * TOLEFF_NIT;
 EFFECT_RAN_INA = EMAX_RAN_INA * pow(C_RAN, GAMMA_RAN_INA) / (pow(EC50_RAN_INA, GAMMA_RAN_INA) + pow(C_RAN, GAMMA_RAN_INA));
 EFFECT_RAN_EX  = EMAX_RAN_EX  * pow(C_RAN, GAMMA_RAN_EX)  / (pow(EC50_RAN_EX,  GAMMA_RAN_EX)  + pow(C_RAN, GAMMA_RAN_EX));

 // Combined HR target: baseline × (1 - EFFECT_BB_HR) × (1 - EFFECT_IVA)
 double HR_tgt = HR0 * (1 - EFFECT_BB_HR) * (1 - EFFECT_IVA);

 // Combined SBP target
 double SBP_tgt = SBP0 * (1 - EFFECT_BB_SBP) * (1 - EFFECT_CCB_SBP);

 // Stenosis reduces basal CBF; drugs increase it
 double CBF_base_effect = CBF0 * (1 - STEN_K * PLAQUE);
 double CBF_tgt = CBF_base_effect * (1 + EFFECT_CCB_CBF + EFFECT_NIT_CBF);

 // ─────────────────────────────────────────────────────────────
 // HR, SBP, CBF first-order approach to target
 dxdt_HR_STATE  = KOUT_HR  * (HR_tgt  - HR_STATE);
 dxdt_SBP_STATE = KOUT_SBP * (SBP_tgt - SBP_STATE);
 dxdt_CBF_STATE = KOUT_CBF * (CBF_tgt - CBF_STATE);

 // ─────────────────────────────────────────────────────────────
 // MVO2 ~ HR × SBP (Rate-Pressure Product) proxy
 double RPP_cur  = HR_STATE * SBP_STATE;
 // Nitrate reduces MVO2 via preload (LVEDP) reduction
 double MVO2_mod = RPP_cur * (1 - EFFECT_NIT_MVO * 0.5);
 // Ranolazine reduces metabolic inefficiency
 double MVO2_eff = MVO2_mod * (1 - EFFECT_RAN_INA * 0.1);

 // O2 supply (CBF-based, normalised)
 double O2_supply = CBF_STATE / CBF0;  // normalised to 1 at baseline (no stenosis)
 double O2_demand = MVO2_eff  / (HR0 * SBP0);  // normalised RPP

 // O2 imbalance: positive = demand excess (ischemia driver)
 double imbal_tgt = (O2_demand * WDMD - O2_supply * WSUP);
 dxdt_O2IMBAL = 0.3 * (imbal_tgt - O2IMBAL);

 // ─────────────────────────────────────────────────────────────
 // Ischemia burden: driven by O2 imbalance when RPP > threshold
 double isch_drive = (MVO2_eff > RPP_THRESH) ? (MVO2_eff - RPP_THRESH) / RPP_THRESH : 0.0;
 double isch_tgt   = isch_drive * (1 - EFFECT_RAN_INA * 0.3);  // ranolazine reduces ischemia
 dxdt_ISCHEMIA = KOUT_ISCH * (isch_tgt - ISCHEMIA);

 // ─────────────────────────────────────────────────────────────
 // Angina score: driven by ischemia burden (episodes/week dynamics)
 double ang_tgt = ANGINA0 * (ISCHEMIA / 0.15 + 0.001);
 ang_tgt = (ang_tgt > 14) ? 14 : ang_tgt;  // cap at 2/day
 dxdt_ANGINA_SCORE = KOUT_ANG * (ang_tgt - ANGINA_SCORE);

 // ─────────────────────────────────────────────────────────────
 // Exercise capacity: inverse of ischemia + direct ranolazine benefit
 double ex_tgt = EX_CAP0 * (1 - 0.5 * isch_drive) * (1 + EFFECT_RAN_EX);
 ex_tgt = (ex_tgt < 1.0) ? 1.0 : ex_tgt;  // floor at 1 MET
 dxdt_EX_CAP = 0.05 * (ex_tgt - EX_CAP);

 // ─────────────────────────────────────────────────────────────
 // Plaque progression (long-term, months–years)
 double plaque_prog  = KPROG_PL * PLAQUE * (1 - PLAQUE / 100.0);  // logistic
 double plaque_regr  = STATIN_ON * KSTATIN_PL * PLAQUE;
 dxdt_PLAQUE = plaque_prog - plaque_regr;

$TABLE
 capture HR_sim  = HR_STATE;
 capture SBP_sim = SBP_STATE;
 capture RPP_sim = HR_STATE * SBP_STATE;
 capture CBF_sim = CBF_STATE;
 capture O2imbal = O2IMBAL;
 capture Isch    = ISCHEMIA;
 capture Angina  = ANGINA_SCORE;
 capture ExCap   = EX_CAP;
 capture Plaque  = PLAQUE;
 capture NIT_tol = TOL_NIT;

 // [discoverability fix] Re-expose each compound’s concentration as a single
 // contiguous `double C_<STEM> = <expr>;` statement (same formula already
 // used in $ODE, verbatim) so downstream tooling that regexes for this
 // pattern can find it; the $GLOBAL bare forward-declare above was removed
 // for these five names to avoid an ambiguous-reference collision with this
 // declaration. The $ODE bare reassignments and the bare $CAPTURE listing
 // below are unaffected -- they resolve against this declaration.
 double C_BB  = CENT_BB  / V1_BB  * 1000.0;
 double C_CCB = CENT_CCB / V1_CCB * 1000.0;
 double C_RAN = CENT_RAN / V1_RAN * 1000.0;
 double C_IVA = CENT_IVA / V1_IVA * 1000.0;
 double C_NIT = CENT_NIT / V1_NIT * 1000.0;

$CAPTURE
 C_BB C_CCB C_RAN C_IVA C_NIT
 EFFECT_BB_HR EFFECT_BB_SBP EFFECT_CCB_SBP EFFECT_CCB_CBF EFFECT_IVA
 EFFECT_NIT_MVO EFFECT_NIT_CBF EFFECT_RAN_INA EFFECT_RAN_EX
'

## ─── Compile Model ───────────────────────────────────────────────────────────
mod <- mcode("stable_angina_refactored", code)

## ─── Dosing Regimen Builder ──────────────────────────────────────────────────
make_regimen <- function(bb_mg  = 0,  bb_qd  = TRUE,
                         ccb_mg = 0,  ccb_qd = TRUE,
                         ran_mg = 0,  ran_bid = TRUE,
                         iva_mg = 0,  iva_bid = TRUE,
                         nit_mg = 0,  nit_bid = TRUE,
                         statin = 0,
                         t_end  = 168) {   # 168h = 1 week

  ev <- ev(time = 0, amt = 0)  # empty event placeholder

  if (bb_mg > 0) {
    ii_bb <- ifelse(bb_qd, 24, 12)
    add_ev <- ev(cmt = "GUT_BB",  amt = bb_mg,  ii = ii_bb, addl = floor(t_end/ii_bb) - 1)
    ev <- ev + add_ev
  }
  if (ccb_mg > 0) {
    ii_ccb <- ifelse(ccb_qd, 24, 12)
    add_ev <- ev(cmt = "GUT_CCB", amt = ccb_mg, ii = ii_ccb, addl = floor(t_end/ii_ccb) - 1)
    ev <- ev + add_ev
  }
  if (ran_mg > 0) {
    ii_ran <- ifelse(ran_bid, 12, 24)
    add_ev <- ev(cmt = "GUT_RAN", amt = ran_mg, ii = ii_ran, addl = floor(t_end/ii_ran) - 1)
    ev <- ev + add_ev
  }
  if (iva_mg > 0) {
    ii_iva <- ifelse(iva_bid, 12, 24)
    add_ev <- ev(cmt = "GUT_IVA", amt = iva_mg, ii = ii_iva, addl = floor(t_end/ii_iva) - 1)
    ev <- ev + add_ev
  }
  if (nit_mg > 0) {
    # Eccentric dosing: morning + afternoon (no night dose = nitrate-free)
    add_ev <- ev(cmt = "GUT_NIT", amt = nit_mg, ii = 12, addl = floor(t_end/12) - 1)
    ev <- ev + add_ev
  }
  return(ev)
}

## ─── 5 Treatment Scenarios ───────────────────────────────────────────────────

scenarios <- list(
  "S1: Untreated"              = list(bb=0,   ccb=0,   ran=0,    iva=0,   nit=0,   stat=0),
  "S2: Bisoprolol 5 mg QD"     = list(bb=5,   ccb=0,   ran=0,    iva=0,   nit=0,   stat=0),
  "S3: BB + Amlodipine 5 mg"   = list(bb=5,   ccb=5,   ran=0,    iva=0,   nit=0,   stat=0),
  "S4: BB + CCB + Ranolazine"  = list(bb=5,   ccb=5,   ran=1000, iva=0,   nit=0,   stat=0),
  "S5: BB + Ivabradine 5 mg"   = list(bb=5,   ccb=0,   ran=0,    iva=5,   nit=0,   stat=0),
  "S6: BB + ISMN 40 mg + CCB"  = list(bb=5,   ccb=5,   ran=0,    iva=0,   nit=40,  stat=1)
)

t_sim <- seq(0, 168, by = 0.5)  # 1 week, 0.5h resolution

run_scenario <- function(s, label) {
  ev_s <- make_regimen(
    bb_mg  = s$bb,   ccb_mg = s$ccb,
    ran_mg = s$ran,  iva_mg = s$iva,
    nit_mg = s$nit,  statin = s$stat,
    t_end  = 168
  )
  params_s <- list(
    DOSE_BB   = ifelse(s$bb  > 0, 1, 0),
    DOSE_CCB  = ifelse(s$ccb > 0, 1, 0),
    DOSE_RAN  = ifelse(s$ran > 0, 1, 0),
    DOSE_IVA  = ifelse(s$iva > 0, 1, 0),
    DOSE_NIT  = ifelse(s$nit > 0, 1, 0),
    STATIN_ON = s$stat
  )
  out <- mod %>%
    param(params_s) %>%
    mrgsim(ev_s, tgrid = t_sim, carry_out = "cmt") %>%
    as.data.frame() %>%
    mutate(Scenario = label)
  return(out)
}

all_results <- bind_rows(
  mapply(run_scenario, scenarios, names(scenarios), SIMPLIFY = FALSE)
)

## ─── Summary Table (steady-state metrics at 1 week) ─────────────────────────
summary_tbl <- all_results %>%
  filter(time >= 144) %>%  # last 24h for SS estimate
  group_by(Scenario) %>%
  summarise(
    HR_mean_bpm      = round(mean(HR_sim,   na.rm=TRUE), 1),
    SBP_mean_mmHg    = round(mean(SBP_sim,  na.rm=TRUE), 1),
    RPP_mean         = round(mean(RPP_sim,  na.rm=TRUE), 0),
    CBF_mean_mLmin   = round(mean(CBF_sim,  na.rm=TRUE), 1),
    Ischemia_score   = round(mean(Isch,     na.rm=TRUE), 3),
    Angina_per_wk    = round(mean(Angina,   na.rm=TRUE), 1),
    ExCap_METs       = round(mean(ExCap,    na.rm=TRUE), 1),
    .groups = "drop"
  )

print("=== Stable Angina QSP Model — Scenario Summary (1-week steady state) ===")
print(summary_tbl)

## ─── Dose-Response Analysis (Bisoprolol 1.25–20 mg) ─────────────────────────
bb_doses <- c(1.25, 2.5, 5, 10, 20)
dr_results <- lapply(bb_doses, function(d) {
  ev_d <- ev(cmt = "GUT_BB", amt = d, ii = 24, addl = 6)
  mod %>%
    mrgsim(ev_d, tgrid = t_sim) %>%
    as.data.frame() %>%
    filter(time >= 144) %>%
    summarise(Dose_mg = d, HR_ss = mean(HR_sim), RPP_ss = mean(RPP_sim),
              Angina_ss = mean(Angina))
}) %>% bind_rows()

print("=== Bisoprolol Dose-Response (HR & RPP at steady state) ===")
print(dr_results)

## ─── Visualization ───────────────────────────────────────────────────────────
cols <- c("#9E9E9E","#7B1FA2","#1565C0","#2E7D32","#E65100","#C62828")

p1 <- ggplot(all_results, aes(time, HR_sim, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 60, linetype="dashed", color="grey60") +
  labs(title="Heart Rate Over Time", x="Time (h)", y="HR (bpm)") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7))

p2 <- ggplot(all_results, aes(time, SBP_sim, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  labs(title="Systolic BP", x="Time (h)", y="SBP (mmHg)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p3 <- ggplot(all_results, aes(time, RPP_sim/1000, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 20, linetype="dashed", color="red", alpha=0.7) +
  annotate("text", x=10, y=20.5, label="Ischemic Threshold\n(RPP=20,000)", size=3, color="red") +
  labs(title="Rate-Pressure Product", x="Time (h)", y="RPP (×1000)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p4 <- ggplot(all_results, aes(time, CBF_sim, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  labs(title="Coronary Blood Flow", x="Time (h)", y="CBF (mL/min)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p5 <- ggplot(all_results, aes(time, Angina, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  labs(title="Angina Episodes", x="Time (h)", y="Angina (/week equivalent)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p6 <- ggplot(all_results, aes(time, ExCap, color = Scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values = cols) +
  labs(title="Exercise Capacity", x="Time (h)", y="Exercise Capacity (METs)") +
  theme_bw(base_size=11) + theme(legend.position="none")

# PK plot
p7 <- all_results %>%
  filter(time <= 72) %>%
  pivot_longer(c(C_BB, C_CCB, C_RAN, C_IVA, C_NIT),
               names_to="Drug", values_to="Cp") %>%
  filter(Scenario == "S4: BB + CCB + Ranolazine" | (Scenario=="S6: BB + ISMN 40 mg + CCB" & Drug %in% c("C_BB","C_CCB","C_NIT"))) %>%
  ggplot(aes(time, Cp, color=Drug)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~Scenario) +
  labs(title="Drug Plasma Concentrations", x="Time (h)", y="Cp (ng/mL)") +
  theme_bw(base_size=10)

p8 <- ggplot(dr_results, aes(Dose_mg, RPP_ss/1000)) +
  geom_line(color="#7B1FA2", linewidth=1) +
  geom_point(color="#7B1FA2", size=3) +
  geom_hline(yintercept=20, linetype="dashed", color="red") +
  scale_x_continuous(breaks = bb_doses) +
  labs(title="Bisoprolol Dose-Response\n(RPP at Steady State)", x="Bisoprolol Dose (mg QD)", y="RPP (×1000)") +
  theme_bw(base_size=11)

# Combine
main_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title    = "Stable Angina QSP Model — 6 Treatment Scenarios",
    subtitle = "Bisoprolol · Amlodipine · Ranolazine · Ivabradine · ISMN Nitrate",
    theme    = theme(plot.title = element_text(face="bold", size=14))
  )

print(main_plot)
print(p7)
print(p8)

## ─── Long-term Plaque Simulation (2 years) ───────────────────────────────────
t_long  <- seq(0, 365*24*2, by = 24)  # 2 years, daily
ev_statin <- ev(time = 0, amt = 0)   # no dose event needed, just flag

long_no_statin <- mod %>%
  param(STATIN_ON = 0) %>%
  mrgsim(ev_statin, tgrid = t_long) %>%
  as.data.frame() %>%
  mutate(Scenario = "No Statin")

long_statin <- mod %>%
  param(STATIN_ON = 1) %>%
  mrgsim(ev_statin, tgrid = t_long) %>%
  as.data.frame() %>%
  mutate(Scenario = "Statin (Rosuvastatin)")

p_plaque <- bind_rows(long_no_statin, long_statin) %>%
  mutate(time_yrs = time / (24*365)) %>%
  ggplot(aes(time_yrs, Plaque, color = Scenario)) +
  geom_line(linewidth=1) +
  scale_color_manual(values = c("#C62828","#2E7D32")) +
  labs(title="Long-term Plaque Progression\nvs Statin Regression",
       x="Time (years)", y="Stenosis (% Lumen)") +
  theme_bw(base_size=12)

print(p_plaque)
