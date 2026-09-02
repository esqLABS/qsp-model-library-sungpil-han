## =============================================================================
## Essential Hypertension — mrgsolve QSP/PK-PD Model (PK/PD REFACTOR)
## =============================================================================
## Disease:   Hypertension (본태성 고혈압)
## Model:     Multi-compartment ODE system with RAAS, SNS, endothelial, renal,
##            and cardiac sub-models + PK for 5 antihypertensive drug classes
## ODE states: 18 compartments (the original header comment claims 22; the
##            actual $CMT block — original and here — has 9 drug-PK + 2 RAAS
##            + 5 haemodynamic + 2 remodeling = 18; logged as an original
##            documentation inconsistency, not fixed upstream, see notes)
## Scenarios:  6 treatment scenarios (unchanged from the original)
##
## PK/PD REFACTOR (see FORK_WORKFLOW_GUIDE.md Part 2 and
## eh_refactor_notes.md for the full account). Per the census row for this
## file (driver-patches/data/compound_perturbation_census.md), 5 compounds
## were renamed to the forks pluggable-PK convention. Confirmed real
## identities from the original code/README (the census ACEI/ARB/BB/CCB
## labels are class abbreviations, not the actual modeled molecule):
##   - ACEI stem -> ramipril (dosed prodrug) / ramiprilat (the active
##     metabolite actually tracked by the ODE compartments)
##   - ARB  stem -> losartan (dosed prodrug) / EXP3174 (the active
##     metabolite actually tracked by the ODE compartments)
##   - CCB  stem -> amlodipine
##   - BB   stem -> bisoprolol
##   - HCTZ stem -> hydrochlorothiazide (already a specific real drug; not
##     a class label)
## All five are archetype 2 (no depot, two compartments, linear CL/Q/V) in
## the original ODE system EXCEPT HCTZ, which is archetype 1 (single
## compartment, linear elimination). Renamed compartments and parameters to
## the fork convention: <STEM>_C -> CENT_<STEM>, <STEM>_P -> PERI_<STEM>,
## the exposed concentration to C_<STEM>, and the original's plain-ratio
## effect terms (ACE_inhib, AT1R_block, BB_block, VGCC_block, NCC_inhib)
## renamed to EFFECT_<STEM>, with EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>
## pulled out as explicit named parameters (rename, not a refit — see
## eh_refactor_notes.md). All PK/PD parameter VALUES are copied verbatim
## from the original.
##
## Found in passing, disclosed rather than fixed (not in scope for a
## PK-reorganization refactor): KA_<STEM>, F_<STEM>, and DOSE_<STEM> are
## declared in $PARAM for every one of the 5 compounds but are never
## referenced anywhere in the original's $MAIN/$ODE/$TABLE — there is no
## first-order absorption depot in this model at all; bioavailability is
## applied externally, as a hardcoded numeric literal duplicating F_<STEM>'s
## value, inside the R-side dosing helper (make_dose("ACEI_C", 10 * 0.28 *
## ...)) rather than inside the DSL. Preserved unused, verbatim, for
## parameter-value fidelity to the original; does not affect any simulated
## trajectory since they are dead parameters in both files.
##
## Parameter calibration references (unchanged from original):
##   - ACEI PK:  Breslin et al., Clin Pharmacokinet 2003
##   - ARB PK:   Gottwald et al., Clin Pharmacokinet 2002
##   - CCB PK:   Faulkner et al., J Cardiovasc Pharmacol 1986
##   - BB PK:    Leopold et al., Eur J Clin Pharmacol 1986
##   - HCTZ PK:  Beermann, Eur J Clin Pharmacol 1984
##   - RAAS PD:  Mager et al., J Pharmacokinet Pharmacodyn 2003
##   - BP model: Pruijm et al., Am J Hypertens 2013
##   - LVH model:Devereux et al., J Am Coll Cardiol 1989
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## =============================================================================
## 1. MODEL DEFINITION
## =============================================================================

eh_code_refactored <- '

$PROB Essential Hypertension QSP Model v1.0 — refactored (pluggable PK)
  18 ODE states | 5 drug classes | 6 treatment scenarios

$PARAM @annotated
  // -- ACE Inhibitor (ramipril prodrug -> ramiprilat active metabolite) PK --
  KA_ACEI  : 1.2    : Absorption rate constant (1/h), ramipril (declared, unused: no depot in this model)
  F_ACEI   : 0.28   : Oral bioavailability (fraction), ramipril to ramiprilat (declared, unused: applied externally at dose time)
  V1_ACEI  : 8.0    : Central volume of distribution (L), ramiprilat
  V2_ACEI  : 32.0   : Peripheral volume of distribution (L), ramiprilat
  CL_ACEI  : 6.5    : Clearance (L/h), ramiprilat (renal)
  Q_ACEI   : 3.0    : Intercompartmental clearance (L/h)
  DOSE_ACEI: 0.0    : Dose flag for ACE inhibitor (0 = off) (declared, unused)
  EMAX_ACEI : 1.0   : Emax of ramiprilat ACE inhibition (implicit in original; =1, rename not a fit)
  EC50_ACEI : 0.005 : EC50 of ramiprilat for ACE (mg/L equiv.) (renamed from IC50_ACEI)
  GAMMA_ACEI: 1.0   : Hill coefficient, ACE inhibition (none in original; =1, rename not a fit)

  // -- ARB (losartan prodrug -> EXP3174 active metabolite) PK --
  KA_ARB   : 1.1    : Absorption rate constant (1/h), losartan (declared, unused: no depot in this model)
  F_ARB    : 0.33   : Oral bioavailability (fraction), EXP3174 (declared, unused: applied externally at dose time)
  V1_ARB   : 14.0   : Central Vc (L), EXP3174
  V2_ARB   : 45.0   : Peripheral Vp (L), EXP3174
  CL_ARB   : 5.2    : Clearance (L/h), EXP3174
  Q_ARB    : 3.5    : Intercompartmental clearance (L/h)
  DOSE_ARB : 0.0    : Dose flag for ARB (0 = off) (declared, unused)
  EMAX_ARB : 1.0    : Emax of EXP3174 AT1R blockade (implicit in original; =1, rename not a fit)
  EC50_ARB : 0.02   : EC50 of EXP3174 for AT1R (mg/L equiv.) (renamed from IC50_ARB)
  GAMMA_ARB: 1.0    : Hill coefficient, AT1R blockade (none in original; =1, rename not a fit)

  // -- CCB (amlodipine) PK --
  KA_CCB   : 0.25   : Absorption rate constant (1/h), amlodipine (slow) (declared, unused: no depot in this model)
  F_CCB    : 0.64   : Oral bioavailability (fraction), amlodipine (declared, unused: applied externally at dose time)
  V1_CCB   : 21.0   : Central Vc (L), amlodipine
  V2_CCB   : 400.0  : Peripheral Vp (L), amlodipine (large Vd ~21 L/kg)
  CL_CCB   : 3.5    : Clearance (L/h), amlodipine (hepatic CYP3A4)
  Q_CCB    : 8.0    : Intercompartmental clearance (L/h)
  DOSE_CCB : 0.0    : Dose flag for CCB (0 = off) (declared, unused)
  EMAX_CCB : 1.0    : Emax of amlodipine L-VGCC blockade (implicit in original; =1, rename not a fit)
  EC50_CCB : 0.003  : EC50 of amlodipine for L-VGCC (mg/L equiv.) (pulled from an inline literal in the original, not a named parameter there)
  GAMMA_CCB: 1.0    : Hill coefficient, L-VGCC blockade (none in original; =1, rename not a fit)

  // -- Beta-blocker (bisoprolol) PK --
  KA_BB    : 1.3    : Absorption rate constant (1/h), bisoprolol (declared, unused: no depot in this model)
  F_BB     : 0.80   : Oral bioavailability (fraction), bisoprolol (declared, unused: applied externally at dose time)
  V1_BB    : 12.0   : Central Vc (L), bisoprolol
  V2_BB    : 100.0  : Peripheral Vp (L), bisoprolol
  CL_BB    : 9.0    : Clearance (L/h), bisoprolol (renal + hepatic)
  Q_BB     : 4.5    : Intercompartmental clearance (L/h)
  DOSE_BB  : 0.0    : Dose flag for beta-blocker (0 = off) (declared, unused)
  EMAX_BB  : 1.0    : Emax of bisoprolol beta1-AR blockade (implicit in original; =1, rename not a fit)
  EC50_BB  : 0.10   : EC50 of bisoprolol for beta1-AR (mg/L equiv.) (renamed from IC50_BB)
  GAMMA_BB : 1.0    : Hill coefficient, beta1-AR blockade (none in original; =1, rename not a fit)

  // -- Thiazide Diuretic (HCTZ) PK --
  KA_HCTZ  : 1.5    : Absorption rate constant (1/h), HCTZ (declared, unused: no depot in this model)
  F_HCTZ   : 0.70   : Oral bioavailability (fraction), HCTZ (declared, unused: applied externally at dose time)
  V1_HCTZ  : 4.0    : Central Vc (L), HCTZ
  CL_HCTZ  : 18.0   : Clearance (L/h), HCTZ (renal)
  DOSE_HCTZ: 0.0    : Dose flag for HCTZ (0 = off) (declared, unused)
  EMAX_HCTZ : 1.0   : Emax of HCTZ NCC inhibition (implicit in original; =1, rename not a fit)
  EC50_HCTZ : 0.02  : EC50 of HCTZ for NCC (mg/L equiv.) (pulled from an inline literal in the original, not a named parameter there)
  GAMMA_HCTZ: 1.0   : Hill coefficient, NCC inhibition (none in original; =1, rename not a fit)

  // -- RAAS PD Parameters --
  ANGII0   : 15.0   : Baseline AngII (pg/mL), normal = 8-25 pg/mL
  KPROD_AII: 0.15   : AngII production rate constant (1/h)
  KDEG_AII : 0.15   : AngII degradation rate constant (1/h)
  ALDO0    : 180.0  : Baseline aldosterone (pmol/L), normal 110-860
  KPROD_AL : 0.12   : Aldosterone production rate (1/h)
  KDEG_AL  : 0.22   : Aldosterone degradation rate (1/h)

  // -- Sympathetic Tone --
  SNS0     : 1.0    : Baseline sympathetic tone (normalized, 1 = normal)
  KRET_SNS : 0.30   : SNS tone return-to-baseline rate (1/h)

  // -- Nitric Oxide / Endothelial --
  NO0      : 1.0    : Baseline NO index (normalized)
  KPROD_NO : 0.50   : NO production rate (1/h)
  KDEG_NO  : 0.50   : NO degradation rate (1/h)

  // -- TPR dynamics --
  TPR0     : 1.0    : Baseline TPR (normalized, =1 -> MAP_baseline)
  KTPR_RET : 0.08   : TPR return-to-baseline rate constant (1/h)
  // Contributions to TPR change
  ALPHA_AII: 0.40   : AngII contribution coefficient to TPR
  ALPHA_SNS: 0.30   : SNS tone contribution coefficient to TPR
  ALPHA_NO : 0.20   : NO-mediated vasodilation coefficient (reduces TPR)
  ALPHA_CCB: 0.25   : CCB-mediated vasodilation coefficient

  // -- Cardiac PD --
  HR0      : 70.0   : Baseline heart rate (bpm)
  SV0      : 70.0   : Baseline stroke volume (mL)
  ALPHA_HR : 0.20   : SNS effect coefficient on HR
  BETA_BB  : 0.30   : Beta-blocker reduction in HR (fraction)

  // -- Plasma Volume / Na --
  PV0      : 3.2    : Baseline plasma volume (L)
  KPVRET   : 0.06   : Plasma volume return rate (1/h)
  HCTZ_PV  : 0.10   : HCTZ-induced plasma volume reduction coefficient

  // -- Mean Arterial Pressure --
  MAP0     : 100.0  : Baseline MAP (mmHg) - pre-treatment hypertensive patient
  PP0      : 50.0   : Baseline pulse pressure (mmHg)

  // -- LV Hypertrophy (chronic remodeling) --
  LVM0     : 210.0  : Baseline LV mass (g) - hypertensive (normal < 200 g)
  KLVM_ON  : 0.002  : LV hypertrophy growth rate (g/h per mmHg above 93)
  KLVM_RET : 0.0005 : LV mass regression rate (1/h) with treatment
  MAP_THRESH: 93.0  : MAP threshold above which LVH progresses

  // -- eGFR --
  EGFR0    : 72.0   : Baseline eGFR (mL/min/1.73m2), mildly reduced
  KEGFR_DEC: 0.0001 : eGFR decline rate (/h per mmHg above MAP_THRESH)
  KEGFR_RET: 0.0003 : eGFR improvement rate with MAP control

$CMT @annotated
  // Drug PK compartments (9 compartments, renamed to the fork convention)
  CENT_ACEI : ACE inhibitor (ramiprilat) central (mg)
  PERI_ACEI : ACE inhibitor (ramiprilat) peripheral (mg)
  CENT_ARB  : ARB (EXP3174) central (mg)
  PERI_ARB  : ARB (EXP3174) peripheral (mg)
  CENT_CCB  : CCB (amlodipine) central (mg)
  PERI_CCB  : CCB (amlodipine) peripheral (mg)
  CENT_BB   : Beta-blocker (bisoprolol) central (mg)
  PERI_BB   : Beta-blocker (bisoprolol) peripheral (mg)
  CENT_HCTZ : Thiazide diuretic (HCTZ) central (mg)

  // RAAS/PD compartments (2 compartments)
  ANGII    : Angiotensin II concentration (pg/mL)
  ALDO     : Aldosterone concentration (pmol/L)

  // Cardiovascular/hemodynamic states (5 compartments)
  SNS_T    : Sympathetic tone (normalized)
  NO_IDX   : Nitric oxide index (normalized)
  TPR_N    : Total peripheral resistance (normalized)
  CO_L     : Cardiac output (L/min)
  PV_L     : Plasma volume (L)

  // Chronic remodeling compartments (2 compartments)
  LVM_G    : LV mass (g)
  EGFR_ML  : eGFR (mL/min/1.73m2)

$MAIN
  // Initial conditions - hypertensive patient (untreated)
  ANGII_0  = ANGII0;
  ALDO_0   = ALDO0;
  SNS_T_0  = SNS0;
  NO_IDX_0 = NO0;
  TPR_N_0  = TPR0;
  CO_L_0   = (MAP0 / 80.0);  // CO such that MAP = CO x TPR_norm x 80 ~ 100
  PV_L_0   = PV0;
  LVM_G_0  = LVM0;
  EGFR_ML_0= EGFR0;

  // Drug depot doses enter via events - initial = 0
  CENT_ACEI_0 = 0; PERI_ACEI_0 = 0;
  CENT_ARB_0  = 0; PERI_ARB_0  = 0;
  CENT_CCB_0  = 0; PERI_CCB_0  = 0;
  CENT_BB_0   = 0; PERI_BB_0   = 0;
  CENT_HCTZ_0 = 0;

$ODE
  // -- Drug PK (archetype 2: no depot, two compartments, linear) --

  // ACE inhibitor (ramiprilat, active metabolite, 2-comp)
  double C_ACEI = CENT_ACEI / V1_ACEI;  // mg/L (the single canonical concentration site)
  dxdt_CENT_ACEI = -CL_ACEI/V1_ACEI * CENT_ACEI - Q_ACEI/V1_ACEI * CENT_ACEI
                   + Q_ACEI/V2_ACEI * PERI_ACEI;
  dxdt_PERI_ACEI =  Q_ACEI/V1_ACEI * CENT_ACEI - Q_ACEI/V2_ACEI * PERI_ACEI;

  // ARB (EXP3174, active metabolite, 2-comp)
  double C_ARB = CENT_ARB / V1_ARB;
  dxdt_CENT_ARB = -CL_ARB/V1_ARB * CENT_ARB - Q_ARB/V1_ARB * CENT_ARB
                  + Q_ARB/V2_ARB * PERI_ARB;
  dxdt_PERI_ARB =  Q_ARB/V1_ARB * CENT_ARB - Q_ARB/V2_ARB * PERI_ARB;

  // CCB (amlodipine, 2-comp)
  double C_CCB = CENT_CCB / V1_CCB;
  dxdt_CENT_CCB = -CL_CCB/V1_CCB * CENT_CCB - Q_CCB/V1_CCB * CENT_CCB
                  + Q_CCB/V2_CCB * PERI_CCB;
  dxdt_PERI_CCB =  Q_CCB/V1_CCB * CENT_CCB - Q_CCB/V2_CCB * PERI_CCB;

  // Beta-blocker (bisoprolol, 2-comp)
  double C_BB = CENT_BB / V1_BB;
  dxdt_CENT_BB = -CL_BB/V1_BB * CENT_BB - Q_BB/V1_BB * CENT_BB
                 + Q_BB/V2_BB * PERI_BB;
  dxdt_PERI_BB =  Q_BB/V1_BB * CENT_BB - Q_BB/V2_BB * PERI_BB;

  // HCTZ (archetype 1: single compartment, renal clearance predominant)
  double C_HCTZ = CENT_HCTZ / V1_HCTZ;
  dxdt_CENT_HCTZ = -CL_HCTZ/V1_HCTZ * CENT_HCTZ;

  // -- Hill interface: one named EFFECT_<STEM> per compound --------------
  // Every one of the original’s effect terms was already a plain ratio
  // C/(C+IC50); this is a rename (EMAX/EC50/GAMMA pulled out explicitly),
  // not a refit. GAMMA_<STEM> = 1 throughout (no Hill coefficient in the
  // original).
  double EFFECT_ACEI = EMAX_ACEI * pow(C_ACEI, GAMMA_ACEI) / (pow(EC50_ACEI, GAMMA_ACEI) + pow(C_ACEI, GAMMA_ACEI));   // ACE inhibition fraction, 0->1
  double EFFECT_ARB  = EMAX_ARB  * pow(C_ARB,  GAMMA_ARB)  / (pow(EC50_ARB,  GAMMA_ARB)  + pow(C_ARB,  GAMMA_ARB));   // AT1R blockade fraction, 0->1
  double EFFECT_BB   = EMAX_BB   * pow(C_BB,   GAMMA_BB)   / (pow(EC50_BB,   GAMMA_BB)   + pow(C_BB,   GAMMA_BB));   // beta1-AR blockade fraction, 0->1
  double EFFECT_CCB  = EMAX_CCB  * pow(C_CCB,  GAMMA_CCB)  / (pow(EC50_CCB,  GAMMA_CCB)  + pow(C_CCB,  GAMMA_CCB));   // L-VGCC blockade fraction, 0->1
  double EFFECT_HCTZ = EMAX_HCTZ * pow(C_HCTZ, GAMMA_HCTZ) / (pow(EC50_HCTZ, GAMMA_HCTZ) + pow(C_HCTZ, GAMMA_HCTZ));  // NCC inhibition fraction, 0->1

  // -- RAAS ODEs -----------------------------------------------------------

  // AngII production reduced by ACE inhibition & AT1R blockade (feedback)
  // AT1R blockade removes negative feedback -> reactive hyperreninaemia
  // ACE inhibition reduces AngII conversion
  double ACE_activity = 1.0 - EFFECT_ACEI;   // fraction of ACE still active
  double AngII_feedback = 1.0 + 0.8 * EFFECT_ARB;  // reactive rise in AngI/AngII
  double AngII_prod = KPROD_AII * ANGII0 * ACE_activity * AngII_feedback;
  double AngII_deg  = KDEG_AII * ANGII;

  dxdt_ANGII = AngII_prod - AngII_deg;

  // Aldosterone driven by AngII (AT1R stimulates adrenal cortex)
  // Both ACEI and ARB reduce aldosterone (but ARB blocks AT1R more directly)
  double Aldo_stim = (ANGII / ANGII0) * (1.0 - EFFECT_ARB);
  dxdt_ALDO = KPROD_AL * ALDO0 * Aldo_stim - KDEG_AL * ALDO;

  // -- Sympathetic Tone ODE -------------------------------------------------
  // Beta-blockers reduce effective sympathetic cardiac output
  // SNS tone itself unchanged but downstream effects modulated
  double SNS_input = SNS0 * (1.0 - 0.15 * EFFECT_BB);  // reduced effect
  dxdt_SNS_T = KRET_SNS * (SNS_input - SNS_T);

  // -- Nitric Oxide Index ODE -------------------------------------------------
  // AngII activates NADPH oxidase -> ROS -> quenches NO
  // ACEI prevents bradykinin degradation -> up-eNOS activation -> up-NO
  double NO_AngII_suppress = 0.3 * (ANGII / ANGII0 - 1.0);  // AngII reduces NO
  double NO_ACEI_boost     = 0.4 * EFFECT_ACEI;             // bradykinin effect
  double NO_target = NO0 * (1.0 - NO_AngII_suppress + NO_ACEI_boost);
  if (NO_target < 0.1) NO_target = 0.1;
  if (NO_target > 2.5) NO_target = 2.5;
  dxdt_NO_IDX = KPROD_NO * NO_target - KDEG_NO * NO_IDX;

  // -- TPR ODE ---------------------------------------------------------------
  // TPR driven by: AngII vasoconstriction, SNS tone, offset by NO and CCB
  double TPR_AngII_effect = ALPHA_AII * (ANGII / ANGII0 - 1.0);
  double TPR_SNS_effect   = ALPHA_SNS * (SNS_T / SNS0 - 1.0);
  double TPR_NO_effect    = ALPHA_NO  * (NO_IDX / NO0 - 1.0);
  double TPR_CCB_effect   = ALPHA_CCB * EFFECT_CCB;
  double TPR_HCTZ_effect  = 0.08 * EFFECT_HCTZ;  // volume depletion -> mild vasoconstriction

  double TPR_target = TPR0 * (1.0
                       + TPR_AngII_effect
                       + TPR_SNS_effect
                       - TPR_NO_effect * (NO_IDX / NO0)
                       - TPR_CCB_effect
                       - EFFECT_ARB * 0.15
                       - EFFECT_ACEI  * 0.12);
  if (TPR_target < 0.3) TPR_target = 0.3;
  dxdt_TPR_N = KTPR_RET * (TPR_target - TPR_N);

  // -- Cardiac Output ODE ---------------------------------------------------
  // HR modulated by SNS and beta-blocker; SV by preload (PV) and afterload (TPR)
  double HR_current = HR0 * (1.0 + ALPHA_HR*(SNS_T/SNS0 - 1.0) - BETA_BB*EFFECT_BB);
  double SV_current = SV0 * (PV_L / PV0) * (1.0 - 0.12*(TPR_N - 1.0));
  if (HR_current < 40) HR_current = 40;
  if (SV_current < 20) SV_current = 20;
  double CO_target = (HR_current * SV_current) / 1000.0;  // L/min
  dxdt_CO_L = 0.5 * (CO_target - CO_L);

  // -- Plasma Volume ODE ----------------------------------------------------
  // HCTZ causes Na+ excretion -> plasma volume reduction
  // Aldosterone retention -> volume expansion
  double Aldo_vol_effect  = 0.08 * (ALDO / ALDO0 - 1.0);
  double HCTZ_vol_effect  = HCTZ_PV * EFFECT_HCTZ;
  double PV_target = PV0 * (1.0 + Aldo_vol_effect - HCTZ_vol_effect);
  if (PV_target < 1.5) PV_target = 1.5;
  dxdt_PV_L = KPVRET * (PV_target - PV_L);

  // -- MAP (algebraic, updated each step) ------------------------------------
  // MAP = CO x TPR x scaling_factor
  // Normalization: at baseline CO=1.25 L/min, TPR_N=1.0 -> MAP=100 mmHg
  double MAP_calc = CO_L * TPR_N * 80.0;  // 80 = mmHg.min/L normalization

  // -- LV Mass ODE (slow remodeling, weeks-months timescale) ------------------
  double LVM_stimulus = (MAP_calc > MAP_THRESH) ? (MAP_calc - MAP_THRESH) : 0.0;
  double LVM_regression = (MAP_calc < MAP_THRESH) ? KLVM_RET * (LVM_G - 180.0) : 0.0;
  dxdt_LVM_G = KLVM_ON * LVM_stimulus - LVM_regression;

  // -- eGFR ODE (slow remodeling) ---------------------------------------------
  double EGFR_decline = (MAP_calc > MAP_THRESH)
                         ? KEGFR_DEC * (MAP_calc - MAP_THRESH) * EGFR_ML
                         : 0.0;
  double EGFR_recovery = (MAP_calc <= MAP_THRESH)
                          ? KEGFR_RET * (EGFR0 - EGFR_ML)
                          : 0.0;
  dxdt_EGFR_ML = -EGFR_decline + EGFR_recovery;

$TABLE
  // -- Derived outputs -------------------------------------------------------
  double MAP_out = CO_L * TPR_N * 80.0;
  double ART_STIFF_proxy = 1.0 + 0.005 * (LVM_G - LVM0);
  // Arterial stiffness proxy increases with LVM growth
  double PP_out  = PP0 * (1.0 + 0.8*(ART_STIFF_proxy - 1.0));

  double SBP_out = MAP_out + (PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0))) * 2.0/3.0;
  double DBP_out = MAP_out - (PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0))) * 1.0/3.0;

  // mrgsolve forward-declares every "double NAME = expr;" seen anywhere in
  // $ODE/$TABLE at file scope, so C_<STEM>/EFFECT_<STEM> are already
  // declared from their $ODE initializers above -- re-assign here (no
  // "double") rather than re-declaring, which would be a redefinition
  // error. This bare-reassignment pattern is what actually lets a value be
  // shared cleanly between blocks. The original recomputed CENT_BB/V1_BB
  // independently in BOTH its HR_out expression and its BB_pct expression (a
  // genuine duplicate concentration site within $TABLE, on top of the
  // separate $ODE site) -- normalized here into the single C_BB below,
  // referenced by both.
  C_ACEI = CENT_ACEI / V1_ACEI;
  C_ARB  = CENT_ARB  / V1_ARB;
  C_CCB  = CENT_CCB  / V1_CCB;
  C_BB   = CENT_BB   / V1_BB;
  C_HCTZ = CENT_HCTZ / V1_HCTZ;

  double HR_out  = HR0 * (1.0 + ALPHA_HR*(SNS_T/SNS0 - 1.0) - BETA_BB*(C_BB/(C_BB + EC50_BB)));
  double SV_out  = (CO_L / HR_out) * 1000.0;  // mL

  double ACE_pct = 100.0 * C_ACEI / (C_ACEI + EC50_ACEI);
  double AT1_pct = 100.0 * C_ARB  / (C_ARB  + EC50_ARB );
  double BB_pct  = 100.0 * C_BB   / (C_BB   + EC50_BB  );
  double CCB_pct = 100.0 * C_CCB  / (C_CCB  + EC50_CCB );
  double HCTZ_pct= 100.0 * C_HCTZ / (C_HCTZ + EC50_HCTZ);

  EFFECT_ACEI = EMAX_ACEI * pow(C_ACEI, GAMMA_ACEI) / (pow(EC50_ACEI, GAMMA_ACEI) + pow(C_ACEI, GAMMA_ACEI));
  EFFECT_ARB  = EMAX_ARB  * pow(C_ARB,  GAMMA_ARB)  / (pow(EC50_ARB,  GAMMA_ARB)  + pow(C_ARB,  GAMMA_ARB));
  EFFECT_BB   = EMAX_BB   * pow(C_BB,   GAMMA_BB)   / (pow(EC50_BB,   GAMMA_BB)   + pow(C_BB,   GAMMA_BB));
  EFFECT_CCB  = EMAX_CCB  * pow(C_CCB,  GAMMA_CCB)  / (pow(EC50_CCB,  GAMMA_CCB)  + pow(C_CCB,  GAMMA_CCB));
  EFFECT_HCTZ = EMAX_HCTZ * pow(C_HCTZ, GAMMA_HCTZ) / (pow(EC50_HCTZ, GAMMA_HCTZ) + pow(C_HCTZ, GAMMA_HCTZ));

  capture MAP   = MAP_out;
  capture SBP   = SBP_out;
  capture DBP   = DBP_out;
  capture PP    = PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0));
  capture HR    = HR_out;
  capture CO    = CO_L;
  capture TPR   = TPR_N;
  capture PV    = PV_L;
  capture AngII = ANGII;
  capture Aldo  = ALDO;
  capture NO    = NO_IDX;
  capture LVM   = LVM_G;
  capture eGFR  = EGFR_ML;

  capture Cp_ACEI = C_ACEI;
  capture Cp_ARB  = C_ARB;
  capture Cp_CCB  = C_CCB;
  capture Cp_BB   = C_BB;
  capture Cp_HCTZ = C_HCTZ;

  capture ACE_inhib_pct = ACE_pct;
  capture AT1R_block_pct = AT1_pct;
  capture BB_block_pct   = BB_pct;
  capture VGCC_block_pct = CCB_pct;
  capture NCC_inhib_pct  = HCTZ_pct;

$CAPTURE @annotated
  MAP  : Mean arterial pressure (mmHg)
  SBP  : Systolic blood pressure (mmHg)
  DBP  : Diastolic blood pressure (mmHg)
  PP   : Pulse pressure (mmHg)
  HR   : Heart rate (bpm)
  CO   : Cardiac output (L/min)
  TPR  : Total peripheral resistance (normalized)
  PV   : Plasma volume (L)
  AngII: Angiotensin II (pg/mL)
  Aldo : Aldosterone (pmol/L)
  NO   : Nitric oxide index
  LVM  : LV mass (g)
  eGFR : eGFR (mL/min/1.73m2)
  Cp_ACEI : Ramiprilat plasma concentration (mg/L) (same site as C_ACEI)
  Cp_ARB  : EXP3174 plasma concentration (mg/L) (same site as C_ARB)
  Cp_CCB  : Amlodipine plasma concentration (mg/L) (same site as C_CCB)
  Cp_BB   : Bisoprolol plasma concentration (mg/L) (same site as C_BB)
  Cp_HCTZ : HCTZ plasma concentration (mg/L) (same site as C_HCTZ)
  ACE_inhib_pct  : ACE inhibition (%) (same quantity as EFFECT_ACEI x 100)
  AT1R_block_pct : AT1R blockade (%) (same quantity as EFFECT_ARB x 100)
  BB_block_pct   : beta1-AR blockade (%) (same quantity as EFFECT_BB x 100)
  VGCC_block_pct : L-VGCC blockade (%) (same quantity as EFFECT_CCB x 100)
  NCC_inhib_pct  : NCC inhibition (%) (same quantity as EFFECT_HCTZ x 100)
  C_ACEI : Ramiprilat plasma concentration (mg/L) - canonical driveable site
  C_ARB  : EXP3174 plasma concentration (mg/L) - canonical driveable site
  C_CCB  : Amlodipine plasma concentration (mg/L) - canonical driveable site
  C_BB   : Bisoprolol plasma concentration (mg/L) - canonical driveable site
  C_HCTZ : HCTZ plasma concentration (mg/L) - canonical driveable site
  EFFECT_ACEI : ACE inhibition fraction (0-1)
  EFFECT_ARB  : AT1R blockade fraction (0-1)
  EFFECT_BB   : beta1-AR blockade fraction (0-1)
  EFFECT_CCB  : L-VGCC blockade fraction (0-1)
  EFFECT_HCTZ : NCC inhibition fraction (0-1)
'

## =============================================================================
## 2. COMPILE MODEL
## =============================================================================

mod_eh_refactored <- mread("essential_hypertension_refactored", tempdir(), eh_code_refactored)
cat("Refactored model compiled successfully.\n")
cat("ODE compartments:", length(cmt(mod_eh_refactored)), "\n")


## =============================================================================
## 3. DOSING REGIMENS (6 TREATMENT SCENARIOS) -- identical to the original
## =============================================================================

SIM_DURATION <- 24 * 7 * 24  # 4032 hours = 24 weeks
SIM_DELTA    <- 1             # 1-hour intervals

# Helper: build dosing event for a drug
make_dose <- function(drug_cmt, dose_mg, interval_h = 24,
                      duration = SIM_DURATION) {
  ev(ID = 1, amt = dose_mg, cmt = drug_cmt,
     ii = interval_h, addl = floor(duration / interval_h) - 1)
}

## Scenario 1 - No Treatment (Untreated Hypertension)
dose_S1 <- ev(ID = 1, amt = 0, cmt = "CENT_ACEI", time = 0)

## Scenario 2 - ACE Inhibitor Monotherapy (Ramipril 10 mg QD)
dose_S2 <- make_dose("CENT_ACEI", 10 * 0.28 * 1e3 * 1e-3) # mg (ramiprilat equiv.)

## Scenario 3 - ARB Monotherapy (Losartan 100 mg QD)
dose_S3 <- make_dose("CENT_ARB", 100 * 0.33)  # mg EXP3174 equiv. (33 mg)

## Scenario 4 - CCB Monotherapy (Amlodipine 10 mg QD)
dose_S4 <- make_dose("CENT_CCB", 10 * 0.64)   # mg (6.4 mg)

## Scenario 5 - Beta-Blocker Monotherapy (Bisoprolol 10 mg QD)
dose_S5 <- make_dose("CENT_BB",  10 * 0.80)   # mg (8 mg)

## Scenario 6 - Triple Therapy (ACEI + CCB + Thiazide) - Standard 1st-line
dose_S6 <- make_dose("CENT_ACEI",  5 * 0.28) +
           make_dose("CENT_CCB",   5 * 0.64) +
           make_dose("CENT_HCTZ", 12.5 * 0.70)


## =============================================================================
## 4. RUN SIMULATIONS
## =============================================================================

run_scenario <- function(mod, dose_ev, scenario_name) {
  out <- mod %>%
    ev(dose_ev) %>%
    mrgsim(end = SIM_DURATION, delta = SIM_DELTA) %>%
    as_tibble() %>%
    mutate(scenario = scenario_name,
           time_wk  = time / (24 * 7))
  out
}

results <- bind_rows(
  run_scenario(mod_eh_refactored, dose_S1, "S1: No Treatment"),
  run_scenario(mod_eh_refactored, dose_S2, "S2: ACEI (Ramipril 10 mg)"),
  run_scenario(mod_eh_refactored, dose_S3, "S3: ARB (Losartan 100 mg)"),
  run_scenario(mod_eh_refactored, dose_S4, "S4: CCB (Amlodipine 10 mg)"),
  run_scenario(mod_eh_refactored, dose_S5, "S5: BB (Bisoprolol 10 mg)"),
  run_scenario(mod_eh_refactored, dose_S6, "S6: Triple Therapy")
)

cat("Simulation complete. Rows:", nrow(results), "\n")

cat("\n=== SIMULATION COMPLETE ===\n")
