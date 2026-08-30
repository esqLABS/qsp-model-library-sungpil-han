## ============================================================
## Bronchiectasis QSP Model — mrgsolve ODE Implementation
## REFACTORED PK sibling of bex_mrgsolve_model.R
## ============================================================
## Disease:   Bronchiectasis (non-CF and CF)
## Pathways:  Vicious Cycle (Infection → Neutrophilic Inflammation →
##            Protease Imbalance → Airway Damage → MCC Failure → Re-infection)
## Drugs:     Azithromycin (anti-inflammatory macrolide)
##            Inhaled Tobramycin (anti-pseudomonal)
##            Ciprofloxacin (systemic anti-pseudomonal)
##            Dornase Alfa / Hypertonic Saline (mucoactive)
##            Ivacaftor + Elexacaftor/Tezacaftor (CFTR modulators, CF-BEX)
## Key Refs:  Chalmers et al. EMBARCE 2012; Wong et al. BLESS 2012;
##            Serisier et al. AZISAST 2013; Barker et al. RESPIRE 2018;
##            EMBARCE: Lancet 2012;380:660–667
## Calibration: FEV1 decline ~40mL/yr (severe); exacerbation rate ~3/yr (BSI≥9)
## ============================================================
##
## THIS FILE vs. the original (bex_mrgsolve_model.R):
## Per the fork's PK/PD refactor guide (see FORK_WORKFLOW_GUIDE.md, Part 2),
## the three compounds this model gives their own PK compartments —
## Azithromycin (AZM), Ciprofloxacin (CIPRO), Tobramycin (TOBRA) — have been
## renamed to the guide's naming convention (KA_/F_/CL_/V1_/V2_/Q_<STEM>,
## C_<STEM> exposed concentration, EFFECT_<STEM> Hill interface) and
## reorganized so each compound's PK is a clearly delimited block exposing
## exactly one concentration variable read by its own effect term. No
## equation and no parameter value was changed anywhere in the disease-PD
## core, in the DNase block, or in any of the three refactored compounds'
## own PK/PD arithmetic (every rename is either a 1:1 identifier swap or a
## GAMMA=1 Hill-formula rewrite that is arithmetically identical to the
## original's plain ratio, pow(x,1)==x). See bex_refactor_notes.md for the
## full per-compound account and the verification results.
##
## This file ALSO carries one syntax-only, non-numeric, non-behavioral
## build-compatibility fix that the original needs to compile under
## mrgsolve 2.0.1 at all (the original's $TABLE block ends in three bare
## lowercase "capture ..." lines with no "$CAPTURE" block header, and the
## first of those lines also lists six $CMT compartment names that
## mrgsolve's real $CAPTURE syntax refuses to duplicate) — see
## bex_refactor_notes.md and translations/UPSTREAM_ISSUES.md for the exact
## defect and the exact fix, and confirmation this changes no simulated
## value.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL DEFINITION
## ============================================================
code <- '
$PROB Bronchiectasis QSP Model — Vicious Cycle + Multi-drug PD (refactored PK: AZM/CIPRO/TOBRA pluggable)

$PARAM
// ---- Bacterial Dynamics ----
Bmax        = 1e9,   // Max bacterial load (CFU/mL sputum)
B0          = 1e7,   // Baseline PA load (CFU/mL)
kgrow       = 0.8,   // PA net growth rate (d-1)
kbiofilm    = 0.15,  // Biofilm formation rate constant (d-1)
BF0         = 0.3,   // Baseline biofilm fraction (0-1)
kBF_clear   = 0.05,  // Biofilm clearance rate (d-1, spontaneous)

// ---- Innate Immunity / Neutrophil ----
Nbase       = 2.5e6, // Baseline sputum neutrophil count (cells/mL)
kN_recruit  = 0.6,   // IL-8-driven neutrophil recruitment rate (d-1)
kN_death    = 0.3,   // Neutrophil apoptosis rate (d-1)
kN_bact     = 1e-9,  // Neutrophil bacterial killing efficacy (mL/cell/d)

// ---- IL-8 (CXCL8) Dynamics ----
IL8_base    = 5.0,   // Baseline IL-8 (ng/mL sputum)
kIL8_prod   = 2.0,   // NF-kB-driven IL-8 production rate (ng/mL/d)
kIL8_deg    = 0.5,   // IL-8 degradation rate (d-1)

// ---- Neutrophil Elastase Dynamics ----
NE_base     = 0.8,   // Baseline NE activity (µg/mL sputum)
kNE_prod    = 1.5e-7,// NE release per neutrophil (µg/cell/d)
kNE_deg     = 0.4,   // NE degradation/inhibition rate (d-1)
A1AT_base   = 100,   // Alpha-1 antitrypsin (µg/mL, inhibitory capacity)

// ---- Mucociliary Clearance ----
MCC0        = 1.0,   // Baseline MCC (normalized, 1=normal)
kMCC_NE     = 0.003, // NE-mediated MCC impairment (mL/µg/d)
kMCC_recov  = 0.05,  // MCC recovery rate (d-1)
MCC_min     = 0.1,   // Minimum MCC

// ---- Airway Damage Score ----
AD0         = 0.3,   // Baseline airway damage (0-1, 0=normal, 1=destroyed)
kAD_NE      = 0.002, // NE-driven airway damage rate (mL/µg/d)
kAD_repair  = 0.01,  // Structural repair rate (d-1, very slow)
AD_max      = 0.95,  // Max possible damage

// ---- FEV1 ----
FEV1_0      = 1800,  // Baseline FEV1 (mL, moderate BEX)
kFEV1_AD    = 500,   // FEV1 decrement per damage unit (mL/unit)
FEV1_min    = 400,   // Minimum FEV1 (mL)

// ---- Exacerbation State ----
Ex_threshold = 0.7,  // Bacterial load threshold triggering exacerbation (normalized)
kEx_resolve  = 0.25, // Exacerbation resolution rate (d-1) with treatment

// ---- Sputum Albumin (Inflammation marker) ----
ALB_sput_base = 0.5, // Baseline sputum albumin (mg/mL)

// ---- PK: Azithromycin (AZM) — refactored, archetype 3 (depot+central+peripheral) ----
// Oral 250mg 3x/week or 500mg daily
// Original names in comments: Ka_AZM, F_AZM, Vd_AZM, CL_AZM, Vp_AZM, Q_AZM
KA_AZM      = 0.5,   // (was Ka_AZM) Absorption rate (h-1, as labeled by original — see refactor notes on the day/hour axis question)
F_AZM       = 0.37,  // Bioavailability (37%)
V1_AZM      = 31,    // (was Vd_AZM) Apparent Vd used directly in the elimination/inter-cpt terms below (L/kg per original comment; NOT multiplied by body weight there — see refactor notes)
CL_AZM      = 22,    // Clearance (L/h) → t½ ~ 68h
V2_AZM      = 800,   // (was Vp_AZM) Peripheral Vd (tissue) (L)
Q_AZM       = 5,     // Inter-compartmental CL (L/h)
// PD of AZM (anti-inflammatory) — Hill interface, rename only (original was a plain Emax*C/(EC50+C) ratio, so GAMMA=1 exactly reproduces it)
EC50_AZM_IL8 = 0.05, // AZM lung conc producing 50% IL-8 suppression (mg/L)
EMAX_AZM_IL8 = 0.50, // (was Emax_AZM) Max IL-8/LTB4 suppression fraction
GAMMA_AZM_IL8 = 1,   // no explicit Hill coefficient in the original; added as an explicit rename target
EC50_AZM_QS  = 0.02, // AZM conc for 50% QS inhibition (mg/L)
EMAX_AZM_QS  = 0.35, // (was Emax_AZM_QS) Max quorum sensing inhibition
GAMMA_AZM_QS = 1,    // no explicit Hill coefficient in the original; added as an explicit rename target

// ---- PK: Inhaled Tobramycin (TOBRA) — refactored, bespoke 2-cpt one-way structure ----
// 300mg BID via TOBI Podhaler. LUNG_TOBRA is the dosed, local, PD-exposed
// compartment (bespoke role, not "GUT_"/"CENT_" in the usual sense — see
// refactor notes); CENT_TOBRA is a one-way, informational-only systemic sink.
// Original names in comments: F_Tobra_inh, Vd_Tobra, CL_Tobra, ke_Tobra_lung
F_TOBRA        = 0.10, // (was F_Tobra_inh) fraction of local lung dose reaching the systemic sink
V1_TOBRA       = 0.26, // (was Vd_Tobra) systemic distribution volume (L/kg per original comment; used directly, see notes)
CL_TOBRA       = 5,    // (was CL_Tobra) Systemic CL (L/h, renal) — systemic sink only, not PD-exposed
KE_TOBRA_LUNG  = 0.8,  // (was ke_Tobra_lung) local lung elimination rate constant (bespoke role — no plain CL_/V slot fits a local depot with its own elimination)
// Local lung PD — already exactly Emax*C^gamma/(EC50^gamma+C^gamma) in the original: rename only
MIC_PA_TOBRA   = 4,    // (was MIC_PA_Tobra) MIC of PA to tobramycin (mg/L) — decorative only, unused by any ODE/TABLE expression in the original or here
EC50_TOBRA     = 16,   // (was EC50_Tobra) Tobramycin conc (lung, mg/L) for 50% bacterial kill
EMAX_TOBRA     = 0.85, // (was Emax_Tobra) Max additional bacterial kill fraction
GAMMA_TOBRA    = 1.5,  // (was Hill_Tobra) already an explicit Hill coefficient in the original

// ---- PK: Ciprofloxacin (CIPRO) — refactored, archetype 3 without peripheral (depot+central only) ----
// 500mg BID oral
// Original names in comments: Ka_Cipro, F_Cipro, Vd_Cipro, CL_Cipro, ELF_ratio
KA_CIPRO     = 1.2,  // (was Ka_Cipro) Absorption rate (h-1)
F_CIPRO      = 0.70, // (was F_Cipro) Bioavailability 70%
V1_CIPRO     = 2.5,  // (was Vd_Cipro) Vd (L/kg × 70kg, used directly below)
CL_CIPRO     = 25,   // (was CL_Cipro) CL (L/h) → t½ ~ 4.2h
// Lung ELF penetration ratio ~0.6
ELF_RATIO_CIPRO = 0.60, // (was ELF_ratio)
// PD - AUC/MIC driven — already exactly Emax*C/(EC50+C) in the original: rename only
MIC_PA_CIPRO = 0.5,  // (was MIC_PA_Cipro) MIC ciprofloxacin (mg/L; susceptible) — decorative only, unused by any ODE/TABLE expression in the original or here
EMAX_CIPRO   = 0.80, // (was Emax_Cipro)
EC50_CIPRO   = 1.0,  // (was EC50_Cipro) Lung ELF conc for 50% kill (mg/L)
GAMMA_CIPRO  = 1,    // no explicit Hill coefficient in the original; added as an explicit rename target

// ---- PK: Dornase Alfa (inhaled) — untouched, not one of the 3 scoped compounds ----
// 2.5mg daily → local digestion of eDNA
Kd_DNase     = 0.3,  // DNase elimination from airways (h-1)
Emax_DNase   = 0.60, // Max mucus viscosity reduction
EC50_DNase   = 0.5,  // DNase lung conc (mg/L) for 50% effect

// ---- Simulation switches ----
AZM_on       = 0,    // 1 = AZM therapy active
TIP_on       = 0,    // 1 = Inhaled Tobramycin
Cipro_on     = 0,    // 1 = Oral Ciprofloxacin
DNase_on     = 0,    // 1 = Dornase Alfa

$CMT
// PK compartments
GUT_AZM      // AZM GI compartment (mg)
CENT_AZM     // AZM central plasma (mg)
PERI_AZM     // AZM peripheral tissue (mg)

LUNG_TOBRA   // (was LUNG_Tobra) Inhaled tobramycin lung depot (mg) — dosed, local, PD-exposed
CENT_TOBRA   // (was CENT_Tobra) Tobramycin systemic sink (mg) — one-way inflow only, informational

GUT_CIPRO    // (was GUT_Cipro) Ciprofloxacin GI (mg)
CENT_CIPRO   // (was CENT_Cipro) Ciprofloxacin plasma (mg)

LUNG_DNase   // Dornase alfa lung (mg)

// Disease state compartments
BACT         // Planktonic bacterial load (CFU/mL, ×1e7)
BIOFILM      // Biofilm fraction (0-1)
NEUT         // Airway neutrophil count (×1e6 cells/mL)
IL8          // Sputum IL-8 (ng/mL)
NE           // Neutrophil elastase activity (µg/mL)
MCC          // Mucociliary clearance index (0-1)
AD           // Airway damage score (0-1)
EXAC         // Exacerbation state (0 = quiescent, 1 = active)

$MAIN
// Derived concentrations
// AZM: exposed C_AZM is the LUNG concentration (what both AZM PD terms read),
// C_AZM_PLASMA is informational only (kept for reporting, matches original
// C_AZM_plasma; note it uses V1_AZM*70 while the ODE below uses V1_AZM alone
// -- a pre-existing internal inconsistency in the original, preserved as-is,
// see refactor notes; not fixed here).
double C_AZM_PLASMA = CENT_AZM / (V1_AZM * 70);  // mg/L
double C_AZM        = C_AZM_PLASMA * 10;         // exposed conc (lung ~10x plasma, macrolide)

// TOBRA: exposed C_TOBRA is the local lung concentration (what the PD term
// reads); C_TOBRA_SYS is informational only, and — exactly as in the
// original (C_Tobra_sys) — is never read by anything else in this model.
double C_TOBRA     = LUNG_TOBRA / 60;            // mg/L (lung ~60mL ELF equiv) -- exposed
double C_TOBRA_SYS = CENT_TOBRA / (V1_TOBRA * 70); // informational only, unused (dead, same as original)

// CIPRO: exposed C_CIPRO is the ELF (lung) concentration (what the PD term
// reads); C_CIPRO_PLASMA is informational only.
double C_CIPRO_PLASMA = CENT_CIPRO / (V1_CIPRO * 70);
double C_CIPRO        = C_CIPRO_PLASMA * ELF_RATIO_CIPRO; // exposed

double C_DNase_lung  = LUNG_DNase / 30;             // mg/L

// Initial conditions
BACT_0 = B0 / 1e7;      // normalized to ×1e7
NEUT_0 = Nbase / 1e6;
IL8_0  = IL8_base;
NE_0   = NE_base;
MCC_0  = MCC0;
AD_0   = AD0;
BIOFILM_0 = BF0;
EXAC_0    = 0;

$ODE
// ============================================================
// DRUG PK ODEs
// ============================================================
// AZM — archetype 3 (depot + central + peripheral), renamed identifiers only.
// (CL_AZM/V1_AZM) + (Q_AZM/V1_AZM) below is (CL_AZM+Q_AZM)/V1_AZM, arithmetically
// identical to what the original computed as two separate terms.
dxdt_GUT_AZM   = -KA_AZM * GUT_AZM;
dxdt_CENT_AZM  = KA_AZM * GUT_AZM * F_AZM
                 - (CL_AZM/V1_AZM) * CENT_AZM
                 - (Q_AZM/V1_AZM) * CENT_AZM
                 + (Q_AZM/V2_AZM) * PERI_AZM;
dxdt_PERI_AZM  = (Q_AZM/V1_AZM) * CENT_AZM
                 - (Q_AZM/V2_AZM) * PERI_AZM;

// Inhaled Tobramycin — bespoke 2-cpt, one-way (LUNG_TOBRA -> CENT_TOBRA),
// renamed identifiers only; arithmetic identical to the original.
dxdt_LUNG_TOBRA  = -KE_TOBRA_LUNG * LUNG_TOBRA
                   - (F_TOBRA * KE_TOBRA_LUNG) * LUNG_TOBRA; // partial absorption
dxdt_CENT_TOBRA  = F_TOBRA * KE_TOBRA_LUNG * LUNG_TOBRA
                   - (CL_TOBRA / (V1_TOBRA * 70)) * CENT_TOBRA;

// Ciprofloxacin — archetype 3 without a peripheral compartment (depot +
// central only), renamed identifiers only.
dxdt_GUT_CIPRO   = -KA_CIPRO * GUT_CIPRO;
dxdt_CENT_CIPRO  = KA_CIPRO * GUT_CIPRO * F_CIPRO
                   - (CL_CIPRO / (V1_CIPRO * 70)) * CENT_CIPRO;

// DNase — untouched, not one of the 3 scoped compounds
dxdt_LUNG_DNase  = -Kd_DNase * LUNG_DNase;

// ============================================================
// DRUG PD — Effect calculations (Hill interface: EFFECT_<STEM>, from C_<STEM>)
// ============================================================
// AZM anti-inflammatory: IL-8 suppression (rename only; GAMMA=1 reproduces
// the original plain Emax*C/(EC50+C) ratio exactly, pow(x,1)==x)
double EFFECT_AZM_IL8 = EMAX_AZM_IL8 * pow(C_AZM, GAMMA_AZM_IL8)
                        / (pow(EC50_AZM_IL8, GAMMA_AZM_IL8) + pow(C_AZM, GAMMA_AZM_IL8));
// AZM quorum sensing inhibition → reduces biofilm growth (rename only)
double EFFECT_AZM_QS  = EMAX_AZM_QS * pow(C_AZM, GAMMA_AZM_QS)
                        / (pow(EC50_AZM_QS, GAMMA_AZM_QS) + pow(C_AZM, GAMMA_AZM_QS));

// Tobramycin — direct bactericidal on planktonic PA (already Hill-shaped in
// the original; rename only, GAMMA_TOBRA keeps the originally fitted value 1.5)
double EFFECT_TOBRA = EMAX_TOBRA * pow(C_TOBRA, GAMMA_TOBRA) /
                      (pow(EC50_TOBRA, GAMMA_TOBRA) + pow(C_TOBRA, GAMMA_TOBRA));

// Ciprofloxacin — concentration-dependent kill, AUC/MIC (rename only;
// GAMMA_CIPRO=1 reproduces the original plain ratio exactly)
double EFFECT_CIPRO = EMAX_CIPRO * pow(C_CIPRO, GAMMA_CIPRO)
                      / (pow(EC50_CIPRO, GAMMA_CIPRO) + pow(C_CIPRO, GAMMA_CIPRO));

// Total bactericidal effect (additive, capped at 0.95)
double E_bact_kill = 1 - (1 - EFFECT_TOBRA * TIP_on) * (1 - EFFECT_CIPRO * Cipro_on);
if(E_bact_kill > 0.95) E_bact_kill = 0.95;

// DNase effect — mucus viscosity reduction → MCC enhancement (untouched)
double E_DNase_MCC = Emax_DNase * C_DNase_lung / (EC50_DNase + C_DNase_lung);

// ============================================================
// DISEASE STATE ODEs
// ============================================================

// Bacterial load (normalized, units ×1e7 CFU/mL)
// Growth logistic, MCC-dependent clearance, antibiotic kill
double BACT_n   = BACT / (Bmax / 1e7);         // normalized 0-1
double kBact_MCC = kN_bact * NEUT * 1e6;        // innate killing (mL/cell/d × cells/mL → d-1)
double kBact_drug = E_bact_kill * 2.0;           // drug-mediated kill (d-1 max)

dxdt_BACT = kgrow * BACT * (1 - BACT_n)         // logistic growth
            - kBact_MCC * BACT                   // neutrophil killing
            - MCC * 0.4 * BACT                   // MCC-dependent clearance
            - kBact_drug * BACT;                 // antibiotic effect

// Biofilm fraction (0-1)
double Biofilm_growth = kbiofilm * BACT_n * (1 - BIOFILM) * (1 - EFFECT_AZM_QS * AZM_on);
double Biofilm_clear  = kBF_clear * BIOFILM * MCC;
dxdt_BIOFILM = Biofilm_growth - Biofilm_clear;

// IL-8 (ng/mL sputum)
// Produced by NF-kB in response to bacteria; suppressed by AZM
double IL8_prod = kIL8_prod * BACT_n * (1 - EFFECT_AZM_IL8 * AZM_on);
dxdt_IL8 = IL8_prod + 0.3 * IL8_base          // baseline production
           - kIL8_deg * (IL8 - IL8_base);

// Neutrophil count (×1e6 cells/mL sputum)
double IL8_stim = IL8 / (IL8 + 5);            // IL-8 saturation function
dxdt_NEUT = kN_recruit * IL8_stim * (Nbase/1e6)
            - kN_death * NEUT
            + 0.1 * (Nbase/1e6 - NEUT);       // homeostatic return

// Neutrophil Elastase (µg/mL sputum)
double NE_inhib_frac = A1AT_base / (A1AT_base + 50); // fractional inhibition
dxdt_NE = kNE_prod * NEUT * 1e6 * (1 - NE_inhib_frac)
          - kNE_deg * (NE - NE_base);

// Mucociliary Clearance Index (0-1; 1=normal)
// Damaged by NE; enhanced by DNase
double MCC_target = MCC0 - kMCC_NE * NE + E_DNase_MCC * DNase_on;
if(MCC_target < MCC_min) MCC_target = MCC_min;
if(MCC_target > 1.0) MCC_target = 1.0;
dxdt_MCC = kMCC_recov * (MCC_target - MCC);

// Airway Damage Score (0-1)
double AD_input = kAD_NE * NE * (1 - AD/AD_max);
double AD_repair = kAD_repair * (1 - AD);
dxdt_AD = AD_input - AD_repair;

// Exacerbation state (quasi-binary, 0-1)
// Triggered when BACT normalized > threshold
double Ex_drive = (BACT_n > Ex_threshold) ? (BACT_n - Ex_threshold) * 5 : 0;
dxdt_EXAC = Ex_drive * (1 - EXAC) - kEx_resolve * EXAC;

$TABLE
// Derived outputs for simulation output
double FEV1 = FEV1_0 - kFEV1_AD * AD;
if(FEV1 < FEV1_min) FEV1 = FEV1_min;
double FEV1_pct = FEV1 / 2400 * 100;        // % predicted (reference 2400 mL)

double BACT_cfu = BACT * 1e7;               // back to absolute CFU/mL
double log10_BACT = log10(BACT_cfu + 1);

double NEUT_abs = NEUT * 1e6;               // absolute neutrophils/mL

double Sputum_purulence = IL8 / (IL8 + 10); // 0-1 index

double AZM_lung_conc = C_AZM;       // was C_AZM_lung
double Tobra_lung_conc = C_TOBRA;   // was C_Tobra_lung
double Cipro_ELF_conc  = C_CIPRO;   // was C_Cipro_ELF

double Exac_active = (EXAC > 0.5) ? 1 : 0;

// NOTE on the two $CAPTURE fixes below (both syntax-only, disclosed in
// bex_refactor_notes.md): (1) the original writes three bare lowercase
// "capture ..." lines with no "$CAPTURE"/"$CMT"-style block header at all,
// which is not valid mrgsolve DSL and fails to compile; (2) once given a
// real block header, the first line is revealed to also list six names
// (IL8, NE, MCC, AD, BIOFILM, EXAC) that duplicate existing $CMT
// compartment names, which mrgsolve 2.0.1 also refuses to build (they are
// always in the output regardless). Both fixes only touch which line/
// symbols are captured, never a value or equation.
$CAPTURE FEV1 FEV1_pct log10_BACT
$CAPTURE Exac_active AZM_lung_conc Tobra_lung_conc Cipro_ELF_conc
$CAPTURE NEUT_abs Sputum_purulence
$CAPTURE C_AZM C_TOBRA C_CIPRO EFFECT_AZM_IL8 EFFECT_AZM_QS EFFECT_TOBRA EFFECT_CIPRO
'

## Compile model
mod <- mcode("bronchiectasis_qsp", code)

## ============================================================
## SIMULATION PARAMETERS & SCENARIOS
## ============================================================

## Observation times (days)
times_long <- seq(0, 365*2, by = 1)    # 2 years daily
times_short <- seq(0, 14, by = 0.5)    # 2-week exacerbation episode

## Initial conditions (moderate-severe BEX, BSI 9-14)
init_vals <- init(mod,
  BACT      = 0.7,   # 7×1e6 CFU/mL (moderate PA colonization)
  BIOFILM   = 0.40,
  NEUT      = 5.0,   # 5×1e6 cells/mL (elevated)
  IL8       = 20,    # ng/mL (elevated)
  NE        = 2.5,   # µg/mL (elevated)
  MCC       = 0.45,  # Impaired MCC
  AD        = 0.50,  # Moderate structural damage
  EXAC      = 0
)

## Helper — dosing event table
## NOTE: cmt = "LUNG_TOBRA" / "GUT_CIPRO" here (was "LUNG_Tobra" / "GUT_Cipro"
## in the original) — renamed to match the refactored compartment names above.
make_events <- function(dose_AZM = 0, dose_TIP = 0, dose_Cipro = 0, dose_DNase = 0,
                        days = 365, freq_AZM = 3, freq_TIP = 0.5,
                        freq_Cipro = 2, freq_DNase = 1) {
  evs <- c()
  if(dose_AZM > 0) {
    # AZM: e.g. 250mg 3x/week → every ~2.3 days
    t_azm <- seq(0, days, by = 7/freq_AZM)
    evs <- c(evs, ev(amt = dose_AZM, cmt = "GUT_AZM", time = t_azm))
  }
  if(dose_TIP > 0) {
    # Tobramycin inhaled BID = twice daily
    t_tip <- seq(0, days, by = 1/freq_TIP)
    evs <- c(evs, ev(amt = dose_TIP, cmt = "LUNG_TOBRA", time = t_tip))
  }
  if(dose_Cipro > 0) {
    # Ciprofloxacin 500mg BID
    t_cip <- seq(0, days, by = 1/freq_Cipro)
    evs <- c(evs, ev(amt = dose_Cipro, cmt = "GUT_CIPRO", time = t_cip))
  }
  if(dose_DNase > 0) {
    t_dn <- seq(0, days, by = 1)
    evs <- c(evs, ev(amt = dose_DNase, cmt = "LUNG_DNase", time = t_dn))
  }
  if(length(evs) == 0) return(ev(amt = 0, cmt = 1, time = Inf))
  do.call(c, evs)
}

## ============================================================
## SCENARIO 1: No Treatment (Natural history)
## ============================================================
cat("Running Scenario 1: No treatment...\n")
s1_params <- param(mod, AZM_on = 0, TIP_on = 0, Cipro_on = 0, DNase_on = 0)
s1 <- mrgsim(s1_params, init = init_vals,
             tgrid = times_long, carry_out = "evid") %>%
  as_tibble() %>%
  mutate(scenario = "No Treatment")

## ============================================================
## SCENARIO 2: Azithromycin 250mg 3×/week (EMBRACE protocol)
## ============================================================
cat("Running Scenario 2: Azithromycin maintenance...\n")
s2_params <- param(mod, AZM_on = 1, TIP_on = 0, Cipro_on = 0, DNase_on = 0)
s2_ev <- make_events(dose_AZM = 250, days = 365*2)
s2 <- mrgsim(s2_params, init = init_vals, events = s2_ev,
             tgrid = times_long) %>%
  as_tibble() %>%
  mutate(scenario = "Azithromycin 250mg 3×/wk")

## ============================================================
## SCENARIO 3: Inhaled Tobramycin 300mg BID (RESPIRE protocol)
## 28 days on / 28 days off cycling
## ============================================================
cat("Running Scenario 3: Inhaled Tobramycin (on/off cycle)...\n")
s3_params <- param(mod, AZM_on = 0, TIP_on = 1, Cipro_on = 0, DNase_on = 0)
# On-cycles: months 1,3,5... Off-cycles: months 2,4,6...
on_days <- unlist(lapply(seq(0, 365*2 - 56, by = 56), function(d) seq(d, d+27)))
t_tip_on <- sort(unique(unlist(lapply(on_days, function(d) c(d, d + 0.5)))))
t_tip_on <- t_tip_on[t_tip_on <= 365*2]
s3_ev <- ev(amt = 300, cmt = "LUNG_TOBRA", time = t_tip_on)
s3 <- mrgsim(s3_params, init = init_vals, events = s3_ev,
             tgrid = times_long) %>%
  as_tibble() %>%
  mutate(scenario = "Inhaled Tobramycin 300mg BID (cycled)")

## ============================================================
## SCENARIO 4: Combination — AZM + Inhaled Tobramycin
## ============================================================
cat("Running Scenario 4: AZM + Inhaled Tobramycin combination...\n")
s4_params <- param(mod, AZM_on = 1, TIP_on = 1, Cipro_on = 0, DNase_on = 0)
s4_ev <- do.call(c, list(make_events(dose_AZM = 250, days = 365*2),
                          ev(amt = 300, cmt = "LUNG_TOBRA", time = t_tip_on)))
s4 <- mrgsim(s4_params, init = init_vals, events = s4_ev,
             tgrid = times_long) %>%
  as_tibble() %>%
  mutate(scenario = "AZM + Inhaled Tobramycin")

## ============================================================
## SCENARIO 5: AZM + Inhaled Tobramycin + Dornase Alfa
## (Intensive maintenance — CF-like protocol)
## ============================================================
cat("Running Scenario 5: AZM + TIP + Dornase Alfa...\n")
s5_params <- param(mod, AZM_on = 1, TIP_on = 1, Cipro_on = 0, DNase_on = 1)
s5_ev_dn  <- make_events(dose_AZM = 250, dose_TIP = 300, dose_DNase = 2.5,
                          days = 365*2, freq_TIP = 0.5)
s5 <- mrgsim(s5_params, init = init_vals, events = s5_ev_dn,
             tgrid = times_long) %>%
  as_tibble() %>%
  mutate(scenario = "AZM + TIP + Dornase Alfa")

## ============================================================
## SCENARIO 6: Acute Exacerbation Treatment
## Oral Ciprofloxacin 500mg BID × 14 days
## ============================================================
cat("Running Scenario 6: Acute exacerbation — Ciprofloxacin 500mg BID × 14d...\n")
init_exac <- init(mod,
  BACT = 2.0, BIOFILM = 0.50, NEUT = 12.0, IL8 = 60, NE = 6.0,
  MCC = 0.25, AD = 0.60, EXAC = 0.85
)
s6_params <- param(mod, AZM_on = 0, TIP_on = 0, Cipro_on = 1, DNase_on = 0)
s6_ev <- make_events(dose_Cipro = 500, days = 14, freq_Cipro = 2)
s6 <- mrgsim(s6_params, init = init_exac, events = s6_ev,
             tgrid = times_short) %>%
  as_tibble() %>%
  mutate(scenario = "Acute Exacerbation: Cipro 500mg BID ×14d")

## ============================================================
## COMBINE SCENARIOS 1-5 (Long-term)
## ============================================================
results_long <- bind_rows(s1, s2, s3, s4, s5)

## ============================================================
## PLOTS
## ============================================================

## Theme
theme_bex <- theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        strip.background = element_rect(fill = "#2196F3"),
        strip.text = element_text(color = "white", face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

col_scen <- c(
  "No Treatment"                          = "#E53935",
  "Azithromycin 250mg 3×/wk"             = "#43A047",
  "Inhaled Tobramycin 300mg BID (cycled)" = "#1E88E5",
  "AZM + Inhaled Tobramycin"              = "#FB8C00",
  "AZM + TIP + Dornase Alfa"             = "#8E24AA"
)

## Figure 1 — FEV1 over 2 years
p1 <- ggplot(results_long, aes(x = time/365, y = FEV1_pct, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_scen) +
  labs(title = "Figure 1 — FEV1 (% predicted) Over 2-Year Treatment Period",
       x = "Time (years)", y = "FEV1 % predicted",
       color = "Treatment Scenario") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.8, y = 52, label = "GOLD Severe BEX cutoff", size = 3) +
  ylim(30, 85) + theme_bex

## Figure 2 — Bacterial Load
p2 <- ggplot(results_long, aes(x = time/365, y = log10_BACT, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_scen) +
  labs(title = "Figure 2 — Sputum Bacterial Load (log10 CFU/mL)",
       x = "Time (years)", y = "log10[CFU/mL]",
       color = "Treatment Scenario") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "gray") +
  annotate("text", x = 1.5, y = 5.15, label = "1×10⁵ CFU/mL", size = 3) +
  theme_bex

## Figure 3 — Sputum IL-8
p3 <- ggplot(results_long, aes(x = time/365, y = IL8, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_scen) +
  labs(title = "Figure 3 — Sputum IL-8 (ng/mL) — Inflammatory Burden",
       x = "Time (years)", y = "IL-8 (ng/mL)",
       color = "Treatment Scenario") +
  theme_bex

## Figure 4 — Neutrophil Elastase
p4 <- ggplot(results_long, aes(x = time/365, y = NE, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_scen) +
  labs(title = "Figure 4 — Sputum Neutrophil Elastase Activity (µg/mL)",
       x = "Time (years)", y = "NE Activity (µg/mL)",
       color = "Treatment Scenario") +
  geom_hline(yintercept = 0.8, linetype = "dotted") +
  theme_bex

## Figure 5 — Mucociliary Clearance & Airway Damage
p5 <- results_long %>%
  select(time, MCC, AD, scenario) %>%
  pivot_longer(c(MCC, AD), names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
    "MCC" = "MCC Index (0-1, 1=normal)",
    "AD"  = "Airway Damage Score (0-1)")) %>%
  ggplot(aes(x = time/365, y = value, color = scenario, linetype = variable)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_scen) +
  scale_linetype_manual(values = c("solid","dashed")) +
  labs(title = "Figure 5 — Mucociliary Clearance & Airway Damage Score",
       x = "Time (years)", y = "Index (0-1)",
       color = "Treatment", linetype = "Variable") +
  theme_bex

## Figure 6 — Acute Exacerbation: Cipro PK/PD
p6 <- ggplot(s6, aes(x = time)) +
  geom_line(aes(y = log10_BACT * 2, color = "Bacterial Load (log10 CFU/mL, ×2)"), linewidth = 1) +
  geom_line(aes(y = Cipro_ELF_conc * 5, color = "Cipro ELF Conc (mg/L, ×5)"), linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = EXAC * 10, color = "Exacerbation Score (×10)"), linewidth = 0.8, linetype = "dotted") +
  scale_color_manual(values = c("Bacterial Load (log10 CFU/mL, ×2)" = "#E53935",
                                "Cipro ELF Conc (mg/L, ×5)" = "#1E88E5",
                                "Exacerbation Score (×10)" = "#FB8C00")) +
  labs(title = "Figure 6 — Acute Exacerbation: Ciprofloxacin 500mg BID × 14 days",
       x = "Time (days)", y = "Normalized Value (see legend for scale factors)",
       color = "Variable") +
  theme_bex

## Print all figures
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

## ============================================================
## SUMMARY TABLE
## ============================================================
summary_tbl <- results_long %>%
  filter(time %in% c(0, 180, 365, 730)) %>%
  group_by(scenario, time) %>%
  summarise(
    FEV1_pct_mean    = round(mean(FEV1_pct), 1),
    log10_BACT_mean  = round(mean(log10_BACT), 2),
    IL8_mean         = round(mean(IL8), 1),
    NE_mean          = round(mean(NE), 2),
    MCC_mean         = round(mean(MCC), 2),
    AD_mean          = round(mean(AD), 2),
    ExacActive_pct   = round(mean(Exac_active) * 100, 0),
    .groups = "drop"
  ) %>%
  mutate(Time_Label = case_when(
    time == 0   ~ "Baseline",
    time == 180 ~ "6 months",
    time == 365 ~ "1 year",
    time == 730 ~ "2 years"
  ))

cat("\n===== Bronchiectasis QSP Model Summary =====\n")
print(as.data.frame(summary_tbl[, c("scenario","Time_Label","FEV1_pct_mean",
                                     "log10_BACT_mean","IL8_mean","NE_mean",
                                     "MCC_mean","AD_mean","ExacActive_pct")]))

## ============================================================
## DOSE-RESPONSE: AZM dose titration (125 / 250 / 500 mg 3×/wk)
## ============================================================
cat("\nRunning AZM dose-response analysis...\n")
dr_results <- purrr::map_dfr(c(0, 125, 250, 500), function(dose) {
  p_dr <- param(mod, AZM_on = ifelse(dose > 0, 1, 0))
  ev_dr <- if(dose > 0) make_events(dose_AZM = dose, days = 365) else
    ev(amt = 0, cmt = 1, time = Inf)
  mrgsim(p_dr, init = init_vals, events = ev_dr,
         tgrid = seq(0, 365, by = 7)) %>%
    as_tibble() %>%
    mutate(AZM_dose = paste0(dose, "mg 3×/wk"))
})

p_dr <- ggplot(dr_results %>% filter(time %in% c(90, 180, 270, 365)),
               aes(x = factor(time), y = FEV1_pct, fill = AZM_dose)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Blues") +
  labs(title = "AZM Dose-Response: FEV1 (% predicted) at 3/6/9/12 months",
       x = "Time (days)", y = "FEV1 % predicted", fill = "AZM Dose") +
  theme_bex

print(p_dr)

cat("\n[Model complete] Bronchiectasis QSP mrgsolve model execution finished.\n")
