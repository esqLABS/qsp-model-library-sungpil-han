## =============================================================================
## Ovarian Cancer (HGSOC) QSP Model — mrgsolve ODE Implementation (REFACTORED)
## High-Grade Serous Ovarian Carcinoma
## Pluggable-PK refactor of oc_mrgsolve_model.R — see oc_refactor_notes.md.
## Original untouched; this is a new sibling file (fork workflow guide, Part 2).
## =============================================================================
## Compartments (18 ODEs):
##   1  CENT_CAR  — Carboplatin central compartment      (was CAR_C1)
##   2  PERI_CAR  — Carboplatin peripheral compartment    (was CAR_C2)
##   3  CENT_PAC  — Paclitaxel central compartment        (was PAC_C1)
##   4  PERI_PAC  — Paclitaxel peripheral compartment      (was PAC_C2)
##   5  PERI2_PAC — Paclitaxel deep peripheral compartment (was PAC_C3)
##   6  GUT_OLA   — Olaparib gut absorption compartment    (was OLA_gut)
##   7  CENT_OLA  — Olaparib central compartment           (was OLA_C1)
##   8  PERI_OLA  — Olaparib peripheral compartment        (was OLA_C2)
##   9  CENT_NIRA — Niraparib central compartment          (was NIRA_C1)
##  10  PERI_NIRA — Niraparib peripheral compartment       (was NIRA_C2)
##  11  CENT_BEV  — Bevacizumab central compartment        (was BEV_C1)
##  12  PERI_BEV  — Bevacizumab peripheral compartment     (was BEV_C2)
##  13  VEGF      — Free VEGF-A (ng/mL)
##  14  TV        — Tumor volume (cm³, Gompertz growth)
##  15  CA125     — CA-125 serum (U/mL)
##  16  Pt_DNA    — Platinum-DNA adducts (relative, 0-1)
##  17  CD8T      — CD8+ T effector cells (relative)
##  18  HRD       — Effective HRD damage accumulation (0-1)
##
## Key References (calibration):
##   - Carboplatin PK: Chatelut 1995 JNCI; CL=GFR×0.134+0.00571×BW
##   - Paclitaxel PK: Gianni 1995 JCO; non-linear Michaelis-Menten
##   - Olaparib PK: Doherty 2014 Clin Pharmacokinet (300mg BID)
##   - Niraparib PK: Sandhu 2013 JCO; Benitez-Llambay 2020
##   - Bevacizumab PK: Lu 2008 Cancer Chemother Pharmacol
##   - PARP inhibitor efficacy: SOLO-1, PRIMA, PAOLA-1 calibration
##   - Tumor growth: Oza 2015 (Gompertz logistic OC model)
##   - CA-125 dynamics: Rustin 2007 (CALYPSO pooled analysis)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---------------------------------------------------------------
## mrgsolve model specification
## ---------------------------------------------------------------
code_oc <- '
$PROB Ovarian Cancer (HGSOC) QSP Model — refactored (pluggable PK)
Carboplatin+Paclitaxel +/- Bevacizumab +/- PARP Inhibitor Maintenance

$PARAM @annotated
// ── Carboplatin PK (Archetype 2: no depot, 2-cmt, linear) ─────────
CL_CAR   : 4.2   : Carboplatin clearance (L/h; Chatelut 1995)
V1_CAR   : 15.0  : Central Vd (L)
Q_CAR    : 6.5   : Inter-compartmental clearance (L/h)
V2_CAR   : 35.0  : Peripheral Vd (L)

// ── Paclitaxel PK (Archetype 2 extended: no depot, 3-cmt, linear) ──
CL_PAC   : 13.2  : Paclitaxel total CL (L/h; Gianni 1995)
V1_PAC   : 6.5   : Central Vd (L)
Q2_PAC   : 7.0   : Rapid inter-comp CL (L/h)
V2_PAC   : 113.0 : Peripheral Vd 2 (L)
Q3_PAC   : 2.2   : Slow inter-comp CL (L/h)
V3_PAC   : 1088.0: Peripheral Vd 3 (deep, L)
EMAX_PAC : 1      : Max fractional G2/M-arrest effect (rename, implicit in original ratio)
EC50_PAC : 100.0  : Paclitaxel EC50 for G2/M arrest (ng/mL)
GAMMA_PAC: 1      : Hill coefficient (original had no explicit Hill term)

// ── Olaparib PK (Archetype 3: depot+central+peripheral, oral, linear) ──
KA_OLA   : 1.74  : Olaparib oral absorption rate (1/h)
F_OLA    : 0.77  : Olaparib oral bioavailability
CL_OLA   : 8.9   : Olaparib apparent CL/F (L/h; Doherty 2014)
V1_OLA   : 67.0  : Olaparib central Vd/F (L)
Q_OLA    : 4.2   : Inter-comp CL (L/h)
V2_OLA   : 100.0 : Peripheral Vd/F (L)
EMAX_OLA : 0.8    : Max fractional PARP-trapping contribution (rename, implicit in original ratio)
EC50_OLA : 500    : Olaparib EC50 for PARP trapping (ng/mL)
GAMMA_OLA: 1      : Hill coefficient (original had no explicit Hill term)

// ── Niraparib PK (bespoke: original lacks a real depot — see notes) ──
KA_NIRA  : 0.36  : Niraparib "absorption" rate constant (1/h; original defect, see refactor notes)
F_NIRA   : 0.73  : Niraparib bioavailability
CL_NIRA  : 16.2  : Niraparib apparent CL/F (L/h)
V1_NIRA  : 537.0 : Niraparib central Vd/F (L)
Q_NIRA   : 5.0   : Inter-comp CL (L/h)
V2_NIRA  : 537.0 : Peripheral Vd/F (L)
EMAX_NIRA : 0.7   : Max fractional PARP-trapping contribution (rename, implicit in original ratio)
EC50_NIRA : 2000  : Niraparib EC50 for PARP trapping (ng/mL)
GAMMA_NIRA: 1     : Hill coefficient (original had no explicit Hill term)

// ── Bevacizumab PK (Archetype 2: no depot, 2-cmt, linear, day-scale) ──
CL_BEV   : 0.207 : Bevacizumab clearance (L/day; Lu 2008)
V1_BEV   : 2.91  : Central Vd (L)
Q_BEV    : 0.469 : Inter-comp CL (L/day)
V2_BEV   : 1.91  : Peripheral Vd (L)

// ── VEGF dynamics ─────────────────────────────────────────────────
VEGF0    : 0.15  : Baseline free VEGF-A (ng/mL; healthy ~0.15)
ksyn_VEGF: 1.5   : VEGF production from tumor (scaled to TV)
kdeg_VEGF: 10.0  : VEGF degradation rate (1/day)
kbind_BEV: 50.0  : Bevacizumab-VEGF binding rate

// ── Platinum-DNA adducts ──────────────────────────────────────────
k_adduct : 0.015 : Rate of Pt-DNA adduct formation (1/(µg/L·h))
k_repair : 0.08  : Adduct repair rate (1/h; NER activity)
HRD_sens  : 1.0  : HRD sensitivity multiplier (1=HRD+; 0.4=HRD-)

// ── HRD accumulation (PARPi effect) ──────────────────────────────
k_HRD_in : 0.1   : HRD damage accrual from PARPi trapping (1/h)
k_HRD_out: 0.02  : HRD repair rate (1/h)

// ── Tumor growth (Gompertz model) ─────────────────────────────────
TV0      : 50.0  : Initial tumor volume (cm³; FIGO III debulked)
kg       : 0.008 : Gompertz growth rate (1/day)
TV_max   : 3000.0: Carrying capacity (cm³)
k_kill_Pt: 0.004 : Platinum kill rate constant (1/(relative adduct·day))
k_kill_T : 0.001 : CD8+ T cell kill rate constant (1/(rel cell·day))
k_kill_Pi: 0.003 : PARPi-induced kill rate (1/day, when HRD high)

// ── CA-125 dynamics ───────────────────────────────────────────────
CA125_BASE : 300.0 : Baseline CA-125 (U/mL; advanced OC) (renamed from CA125_0, see build-compat note)
ksyn_CA125: 3.0  : CA-125 production per unit tumor (U/mL/cm³/day)
kdeg_CA125: 0.03 : CA-125 degradation rate (1/day; t½≈23 days)

// ── CD8+ T cell dynamics ──────────────────────────────────────────
CD8T_BASE : 1.0   : Baseline CD8+ T (relative) (renamed from CD8T_0, see build-compat note)
k_CD8_in : 0.1   : CD8+ influx rate (1/day)
k_CD8_out: 0.1   : CD8+ efflux rate (1/day)
k_exhaust: 0.3   : T cell exhaustion by tumor load
k_ICI    : 2.0   : ICI boost to CD8+ (fold increase)

// ── Scenario flags (0=off, 1=on) ──────────────────────────────────
ICI_flag  : 0    : Immune checkpoint inhibitor (0/1)
BRCAmut   : 1    : BRCA mutation status (1=mut, 0=wt)
HRD_pos   : 1    : HRD positive status (1=HRD+, 0=HRD-)

$CMT @annotated
CENT_CAR  : Carboplatin central (µg/mL)
PERI_CAR  : Carboplatin peripheral (µg/mL)
CENT_PAC  : Paclitaxel central (ng/mL)
PERI_PAC  : Paclitaxel peripheral (ng/mL)
PERI2_PAC : Paclitaxel deep peripheral (ng/mL)
GUT_OLA   : Olaparib gut compartment (mg)
CENT_OLA  : Olaparib central (ng/mL)
PERI_OLA  : Olaparib peripheral (ng/mL)
CENT_NIRA : Niraparib central (ng/mL)
PERI_NIRA : Niraparib peripheral (ng/mL)
CENT_BEV  : Bevacizumab central (mg/L)
PERI_BEV  : Bevacizumab peripheral (mg/L)
VEGF      : Free VEGF-A (ng/mL)
TV        : Tumor volume (cm³)
CA125     : CA-125 serum (U/mL)
Pt_DNA    : Platinum-DNA adducts (relative 0-1)
CD8T      : CD8+ T cell relative level
HRD       : HRD damage accumulation (0-1)

$GLOBAL
double C_CAR, C_PAC, C_OLA, C_NIRA, C_BEV;
double EFFECT_PAC, EFFECT_OLA, EFFECT_NIRA, VEGF_bind_BEV;

$MAIN
// Initial conditions (build-compat: `<CMT>_0` idiom in $MAIN, no separate
// $INIT block — see refactor notes, defect #1)
CENT_CAR_0  = 0;
PERI_CAR_0  = 0;
CENT_PAC_0  = 0;
PERI_PAC_0  = 0;
PERI2_PAC_0 = 0;
GUT_OLA_0   = 0;
CENT_OLA_0  = 0;
PERI_OLA_0  = 0;
CENT_NIRA_0 = 0;
PERI_NIRA_0 = 0;
CENT_BEV_0  = 0;
PERI_BEV_0  = 0;
VEGF_0      = VEGF0;
TV_0        = TV0;
CA125_0     = CA125_BASE;
Pt_DNA_0    = 0;
CD8T_0      = CD8T_BASE;
HRD_0       = 0;

// Effective HRD sensitivity
double eff_HRD = HRD_pos * HRD_sens + (1 - HRD_pos) * 0.3;

$ODE
// ── Carboplatin: Archetype 2 (no depot, 2-cmt, linear) ────────────
// C_CAR is a direct alias, not CENT/V1: the original's own dosing
// function already divides the bolus dose by V1 before adding it to
// the central compartment, and the original's own $TABLE captured the
// central compartment directly as "the concentration" — see notes.
C_CAR = CENT_CAR;
dxdt_CENT_CAR = -(CL_CAR + Q_CAR)/V1_CAR*CENT_CAR + Q_CAR/V2_CAR*PERI_CAR;
dxdt_PERI_CAR =  Q_CAR/V1_CAR*CENT_CAR - Q_CAR/V2_CAR*PERI_CAR;

// ── Paclitaxel: Archetype 2 extended (no depot, 3-cmt, linear) ────
// Despite the header comment claiming "Michaelis-Menten nonlinear" PK
// (Gianni 1995), the actual dxdt equations below are purely linear
// (constant CL_PAC, no Km/Vmax saturation term anywhere) — see notes.
C_PAC = CENT_PAC;
dxdt_CENT_PAC  = -(CL_PAC + Q2_PAC + Q3_PAC)/V1_PAC*CENT_PAC
                  + Q2_PAC/V2_PAC*PERI_PAC + Q3_PAC/V3_PAC*PERI2_PAC;
dxdt_PERI_PAC  =  Q2_PAC/V1_PAC*CENT_PAC - Q2_PAC/V2_PAC*PERI_PAC;
dxdt_PERI2_PAC =  Q3_PAC/V1_PAC*CENT_PAC - Q3_PAC/V3_PAC*PERI2_PAC;
// G2/M arrest — rename only, identical Hill-1 ratio as the original's
// inline `pac_eff = PAC_C1/(PAC_C1+100.0)`
EFFECT_PAC = EMAX_PAC * pow(C_PAC, GAMMA_PAC) / (pow(EC50_PAC, GAMMA_PAC) + pow(C_PAC, GAMMA_PAC));

// ── Olaparib: Archetype 3 (depot+central+peripheral, oral, linear) ─
C_OLA = CENT_OLA;
dxdt_GUT_OLA  = -KA_OLA * GUT_OLA;
dxdt_CENT_OLA =  (F_OLA*KA_OLA*GUT_OLA)/V1_OLA - (CL_OLA + Q_OLA)/V1_OLA*CENT_OLA
                  + Q_OLA/V2_OLA*PERI_OLA;
dxdt_PERI_OLA =  Q_OLA/V1_OLA*CENT_OLA - Q_OLA/V2_OLA*PERI_OLA;
// PARP-trapping contribution — rename only, identical guarded ratio as
// the original's inline `if(OLA_C1>0.1) parp_trap += 0.8*OLA_C1/(OLA_C1+500)`
EFFECT_OLA = 0.0;
if(C_OLA > 0.1) EFFECT_OLA = EMAX_OLA * pow(C_OLA, GAMMA_OLA) / (pow(EC50_OLA, GAMMA_OLA) + pow(C_OLA, GAMMA_OLA));

// ── Niraparib: bespoke — original lacks a real depot compartment ──
// (KA_NIRA/F_NIRA parameters exist and a term applies them, but the
// term reads NIRA_C1/CENT_NIRA itself, not a separate depot state —
// this is preserved exactly as an original defect; see refactor notes)
C_NIRA = CENT_NIRA;
dxdt_CENT_NIRA = (F_NIRA*KA_NIRA*CENT_NIRA)/V1_NIRA - (CL_NIRA + Q_NIRA)/V1_NIRA*CENT_NIRA
                  + Q_NIRA/V2_NIRA*PERI_NIRA;
dxdt_PERI_NIRA =  Q_NIRA/V1_NIRA*CENT_NIRA - Q_NIRA/V2_NIRA*PERI_NIRA;
// PARP-trapping contribution — rename only, identical guarded ratio as
// the original's inline `if(NIRA_C1>0.1) parp_trap += 0.7*NIRA_C1/(NIRA_C1+2000)`
EFFECT_NIRA = 0.0;
if(C_NIRA > 0.1) EFFECT_NIRA = EMAX_NIRA * pow(C_NIRA, GAMMA_NIRA) / (pow(EC50_NIRA, GAMMA_NIRA) + pow(C_NIRA, GAMMA_NIRA));

// ── Bevacizumab: Archetype 2 (no depot, 2-cmt, linear, day-scale) ──
// No TMDD anywhere in this file's own equations: plain linear PK, no
// receptor/complex compartment. VEGF neutralization is a first-order
// mass-action sink (rename of `BEV_effect`), not a Hill/Emax ratio —
// see "No EFFECT_BEV" in the refactor notes for why there is no
// disease-facing Hill term to expose for this compound in this file.
C_BEV = CENT_BEV;
VEGF_bind_BEV = kbind_BEV * C_BEV * VEGF;
dxdt_CENT_BEV = -(CL_BEV + Q_BEV)/V1_BEV*CENT_BEV + Q_BEV/V2_BEV*PERI_BEV;
dxdt_PERI_BEV =  Q_BEV/V1_BEV*CENT_BEV - Q_BEV/V2_BEV*PERI_BEV;

// ── VEGF dynamics ─────────────────────────────────────────────────
double VEGF_prod = ksyn_VEGF * (TV / TV0);
dxdt_VEGF = VEGF_prod - kdeg_VEGF * VEGF - VEGF_bind_BEV;

// ── Platinum-DNA adducts ──────────────────────────────────────────
// HRD patients repair adducts more slowly (eff_HRD reduces repair)
dxdt_Pt_DNA = k_adduct * C_CAR - k_repair * (1 - 0.6 * eff_HRD) * Pt_DNA;

// ── HRD damage (PARP inhibitor trapping + endogenous) ─────────────
// Combined only at the point the disease equation uses it — each
// compound's EFFECT_<STEM> stays separate up to here.
dxdt_HRD = k_HRD_in * (eff_HRD * (EFFECT_OLA + EFFECT_NIRA)) - k_HRD_out * HRD;

// ── CD8+ T cell dynamics ──────────────────────────────────────────
double ICI_effect = 1.0 + ICI_flag * (k_ICI - 1.0);
double exhaustion  = k_exhaust * TV / TV_max;
dxdt_CD8T = k_CD8_in * ICI_effect - k_CD8_out * CD8T - exhaustion * CD8T;

// ── Tumor volume (Gompertz + drug kill) ───────────────────────────
double grow_term = kg * TV * log(TV_max / TV);
// Platinum kill (adduct-dependent)
double kill_Pt   = k_kill_Pt * Pt_DNA * TV;
// G2/M arrest (paclitaxel-dependent)
double kill_pac  = 0.6 * k_kill_Pt * EFFECT_PAC * TV;
// PARPi synthetic lethality
double parp_kill = k_kill_Pi * HRD * TV;
// CD8+ T cell killing
double kill_CD8  = k_kill_T * CD8T * TV;
dxdt_TV = grow_term - kill_Pt - kill_pac - parp_kill - kill_CD8;

// ── CA-125 (turnover, proportional to tumor) ─────────────────────
double CA125_prod = ksyn_CA125 * TV;
dxdt_CA125 = CA125_prod - kdeg_CA125 * CA125;

$TABLE
// Recompute every C_<STEM>/EFFECT_<STEM> directly from live $CMT state
// here too (redundant with $ODE) to avoid a stale $GLOBAL read at a
// timestep where a dose event and a requested output row coincide —
// see the cervical-cancer/cc_refactor_notes.md precedent for why.
C_CAR  = CENT_CAR;
C_PAC  = CENT_PAC;
C_OLA  = CENT_OLA;
C_NIRA = CENT_NIRA;
C_BEV  = CENT_BEV;
EFFECT_PAC = EMAX_PAC * pow(C_PAC, GAMMA_PAC) / (pow(EC50_PAC, GAMMA_PAC) + pow(C_PAC, GAMMA_PAC));
EFFECT_OLA = 0.0;
if(C_OLA > 0.1) EFFECT_OLA = EMAX_OLA * pow(C_OLA, GAMMA_OLA) / (pow(EC50_OLA, GAMMA_OLA) + pow(C_OLA, GAMMA_OLA));
EFFECT_NIRA = 0.0;
if(C_NIRA > 0.1) EFFECT_NIRA = EMAX_NIRA * pow(C_NIRA, GAMMA_NIRA) / (pow(EC50_NIRA, GAMMA_NIRA) + pow(C_NIRA, GAMMA_NIRA));
VEGF_bind_BEV = kbind_BEV * C_BEV * VEGF;

capture CAR_Conc  = C_CAR;
capture PAC_Conc  = C_PAC;
capture OLA_Conc  = C_OLA;
capture NIRA_Conc = C_NIRA;
capture BEV_Conc  = C_BEV;
capture VEGF_free = VEGF;
capture TumorVol  = TV;
capture CA125_lvl = CA125;
capture PtDNA_rel = Pt_DNA;
capture HRD_dmg   = HRD;
capture CD8T_rel  = CD8T;
capture TV_change = (TV - TV0) / TV0 * 100;  // % change from baseline
capture EFFECT_PAC_out  = EFFECT_PAC;
capture EFFECT_OLA_out  = EFFECT_OLA;
capture EFFECT_NIRA_out = EFFECT_NIRA;
capture VEGF_bind_BEV_out = VEGF_bind_BEV;
'

## ---------------------------------------------------------------
## Compile the model
## ---------------------------------------------------------------
mod_oc <- mcode("oc_model", code_oc)

## ---------------------------------------------------------------
## Dosing event functions
## ---------------------------------------------------------------

## Carboplatin: AUC-based (Calvert). GFR=90 mL/min → AUC 6
## Dose (mg) = AUC × (GFR + 25) = 6 × (90+25) = 690 mg
## Approximate concentration: Dose/V1 = 690/15 = 46 µg/mL peak

dose_carboplatin <- function(n_cycles = 6, interval_d = 21, V1 = 15) {
  dose_mg <- 690
  ev(cmt="CENT_CAR", amt=dose_mg/V1, time=seq(0, (n_cycles-1)*interval_d, by=interval_d))
}

## Paclitaxel: 175 mg/m² → 300 mg for 1.72 m² BSA, over 3h
## Peak conc ≈ 300 mg / (6.5 L) ≈ 46,000 ng/mL (non-linear PK adjusts)
dose_paclitaxel <- function(n_cycles = 6, interval_d = 21, V1 = 6.5) {
  dose_mg <- 300
  ev(cmt="CENT_PAC", amt=dose_mg*1000/V1, time=seq(0, (n_cycles-1)*interval_d, by=interval_d))
}

## Olaparib: 300 mg BID (twice daily) starting at time 0
dose_olaparib <- function(start_d = 0, dur_d = 365) {
  times <- seq(start_d*24, (start_d + dur_d)*24 - 12, by=12)
  ev(cmt="GUT_OLA", amt=300, time=times)
}

## Niraparib: 300 mg QD
dose_niraparib <- function(start_d = 0, dur_d = 365) {
  times <- seq(start_d*24, (start_d + dur_d)*24 - 24, by=24)
  ev(cmt="CENT_NIRA", amt=300, time=times)
}

## Bevacizumab: 15 mg/kg q3w IV = ~1050 mg per dose
## Concentration: 1050 mg / V1(L=2.91) = 361 mg/L
dose_bevacizumab <- function(start_d = 0, n_doses = 22, interval_d = 21, V1 = 2.91) {
  dose_mg <- 1050
  times <- seq(start_d, start_d + (n_doses-1)*interval_d, by=interval_d)
  ev(cmt="CENT_BEV", amt=dose_mg/V1, time=times)
}

## ---------------------------------------------------------------
## Treatment Scenarios
## ---------------------------------------------------------------
sim_time <- seq(0, 730, by=1)  # 2-year simulation (days)

## S1: No treatment (natural progression)
mod_S1 <- mod_oc %>% param(BRCAmut=1, HRD_pos=1, ICI_flag=0)
out_S1  <- mrgsim(mod_S1, end=730, delta=1)

## S2: Carboplatin + Paclitaxel × 6 cycles (standard 1st line)
ev_S2 <- ev_seq(dose_carboplatin(6), dose_paclitaxel(6))
out_S2 <- mrgsim(mod_S1, events=ev_S2, end=730, delta=1)

## S3: Carboplatin + Paclitaxel + Bevacizumab × 6 → Bevacizumab maintenance
ev_S3 <- ev_seq(
  dose_carboplatin(6),
  dose_paclitaxel(6),
  dose_bevacizumab(0, n_doses=22)
)
out_S3 <- mrgsim(mod_S1, events=ev_S3, end=730, delta=1)

## S4: Carbo+Pacli × 6 → Olaparib maintenance (BRCA mutant, SOLO-1)
## Olaparib starts day 126 (after 6 cycles)
ev_S4 <- ev_seq(
  dose_carboplatin(6),
  dose_paclitaxel(6),
  dose_olaparib(start_d=126, dur_d=604)
)
mod_S4 <- mod_oc %>% param(BRCAmut=1, HRD_pos=1, ICI_flag=0)
out_S4 <- mrgsim(mod_S4, events=ev_S4, end=730, delta=1)

## S5: Carbo+Pacli × 6 → Niraparib maintenance (all-comers, PRIMA)
ev_S5 <- ev_seq(
  dose_carboplatin(6),
  dose_paclitaxel(6),
  dose_niraparib(start_d=126, dur_d=604)
)
mod_S5 <- mod_oc %>% param(BRCAmut=0, HRD_pos=1, ICI_flag=0)
out_S5 <- mrgsim(mod_S5, events=ev_S5, end=730, delta=1)

## S6: Carbo+Pacli+Bev × 6 → Olaparib+Bev maintenance (PAOLA-1)
ev_S6 <- ev_seq(
  dose_carboplatin(6),
  dose_paclitaxel(6),
  dose_bevacizumab(0, n_doses=6),      # 6 cycles with chemo
  dose_bevacizumab(126, n_doses=16),   # maintenance Bev alone from cycle 7
  dose_olaparib(start_d=126, dur_d=604)
)
mod_S6 <- mod_oc %>% param(BRCAmut=0, HRD_pos=1, ICI_flag=0)
out_S6 <- mrgsim(mod_S6, events=ev_S6, end=730, delta=1)

## ---------------------------------------------------------------
## Summary: 24-month PFS proxy and key endpoints
## ---------------------------------------------------------------
summarize_scenario <- function(out, label) {
  df <- as.data.frame(out)
  # PFS: time when TV > 2× baseline (progressive disease surrogate)
  pfs_d <- df %>% filter(TumorVol > 100) %>% pull(time) %>% min()
  pfs_d <- if(is.infinite(pfs_d)) ">730" else round(pfs_d)
  # CA-125 nadir
  ca125_nadir <- min(df$CA125_lvl)
  ca125_nadir_t <- df$time[which.min(df$CA125_lvl)]
  # Best tumor response
  tv_min <- min(df$TumorVol)
  best_resp <- round((tv_min - 50) / 50 * 100, 1)
  data.frame(
    Scenario      = label,
    PFS_days      = pfs_d,
    CA125_nadir   = round(ca125_nadir, 1),
    CA125_nadir_t = round(ca125_nadir_t),
    BestResp_pct  = best_resp
  )
}

summary_table <- rbind(
  summarize_scenario(out_S1, "S1: Untreated"),
  summarize_scenario(out_S2, "S2: Carbo+Pacli ×6"),
  summarize_scenario(out_S3, "S3: Carbo+Pacli+Bev → Bev maint"),
  summarize_scenario(out_S4, "S4: Carbo+Pacli → Olaparib maint (BRCA+)"),
  summarize_scenario(out_S5, "S5: Carbo+Pacli → Niraparib maint (HRD+)"),
  summarize_scenario(out_S6, "S6: Carbo+Pacli+Bev → Ola+Bev maint (PAOLA-1)")
)
print(summary_table)

## ---------------------------------------------------------------
## Visualization
## ---------------------------------------------------------------

## Helper: combine outputs
combine_sims <- function(..., labels = NULL) {
  sims <- list(...)
  lapply(seq_along(sims), function(i) {
    df <- as.data.frame(sims[[i]])
    df$Scenario <- if(!is.null(labels)) labels[i] else paste0("S", i)
    df
  }) %>% bind_rows()
}

scenario_labels <- c(
  "S1: Untreated",
  "S2: Carbo+Pacli",
  "S3: Carbo+Pacli+Bev→Bev",
  "S4: Carbo+Pacli→Olaparib (BRCA+)",
  "S5: Carbo+Pacli→Niraparib (HRD+)",
  "S6: Carbo+Pacli+Bev→Ola+Bev (PAOLA-1)"
)

all_sims <- combine_sims(out_S1, out_S2, out_S3, out_S4, out_S5, out_S6,
                         labels = scenario_labels)

## --- Plot 1: Tumor Volume over time ---
p1 <- ggplot(all_sims, aes(x=time, y=TumorVol, color=Scenario)) +
  geom_line(size=0.9) +
  geom_hline(yintercept=100, linetype="dashed", color="gray50") +
  annotate("text", x=680, y=110, label="PD threshold (2×BL)", size=2.8, color="gray50") +
  labs(title="Tumor Volume (cm³) — 6 Treatment Scenarios",
       x="Day", y="Tumor Volume (cm³)") +
  scale_y_continuous(limits=c(0, 3000)) +
  scale_color_brewer(palette="Set1") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7))

## --- Plot 2: CA-125 serum ---
p2 <- ggplot(all_sims, aes(x=time, y=CA125_lvl, color=Scenario)) +
  geom_line(size=0.9) +
  geom_hline(yintercept=35, linetype="dashed", color="darkgreen") +
  annotate("text", x=680, y=40, label="ULN 35 U/mL", size=2.8, color="darkgreen") +
  labs(title="CA-125 Serum Level (U/mL)",
       x="Day", y="CA-125 (U/mL)") +
  scale_y_log10() +
  scale_color_brewer(palette="Set1") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

## --- Plot 3: PARPi concentration + HRD ---
p3_ola <- ggplot(as.data.frame(out_S4), aes(x=time)) +
  geom_line(aes(y=OLA_Conc), color="#E91E63") +
  labs(title="Olaparib Central Conc (S4)", x="Day", y="Olaparib (ng/mL)") +
  theme_bw(base_size=10)

p3_hrd <- ggplot(as.data.frame(out_S4), aes(x=time, y=HRD_dmg)) +
  geom_line(color="#880E4F") +
  labs(title="HRD Damage Accumulation (S4 Olaparib maint.)",
       x="Day", y="HRD damage (0–1)") +
  theme_bw(base_size=10)

## --- Plot 4: Drug PK (Carboplatin) ---
p4 <- ggplot(as.data.frame(out_S2), aes(x=time)) +
  geom_line(aes(y=CAR_Conc), color="#FF8F00") +
  labs(title="Carboplatin Central PK (S2)",
       x="Day", y="Carboplatin (µg/mL)") +
  coord_cartesian(xlim=c(0,180)) +
  theme_bw(base_size=10)

## --- Plot 5: VEGF suppression (bevacizumab scenarios) ---
p5 <- ggplot(all_sims %>% filter(Scenario %in% c(
    "S2: Carbo+Pacli",
    "S3: Carbo+Pacli+Bev→Bev",
    "S6: Carbo+Pacli+Bev→Ola+Bev (PAOLA-1)"
  )), aes(x=time, y=VEGF_free, color=Scenario)) +
  geom_line(size=0.9) +
  labs(title="Free VEGF-A (Bevacizumab Scenarios)",
       x="Day", y="Free VEGF-A (ng/mL)") +
  scale_color_brewer(palette="Set2") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7))

## --- Combined figure ---
main_fig <- (p1 | p2) / (p4 | p3_ola) / (p5 | p3_hrd)
print(main_fig + plot_annotation(
  title    = "Ovarian Cancer QSP Model — Simulation Results",
  subtitle = "High-Grade Serous OC · 6 Treatment Scenarios · 2-Year Projection",
  caption  = "Calibrated to SOLO-1, PRIMA, PAOLA-1 clinical trials"
))

## ---------------------------------------------------------------
## Key Parameter Calibration Notes
## ---------------------------------------------------------------
## Carboplatin:
##   - Calvert formula: Dose=AUC×(GFR+25); AUC target 5-6 mg/mL·min
##   - CL primarily renal; t½ alpha=1.1h, beta=6h (Egorin 1984 Cancer Res)
##   - DNA adduct formation peaks 1-2h post-infusion (Bajorin 1992)
##
## Paclitaxel:
##   - Non-linear PK: Michaelis-Menten, Km≈2.17µM (Gianni 1995 JCO)
##   - t½ alpha=0.34h, beta=1.3h, gamma=27h (3-compartment model)
##   - CYP3A4/CYP2C8 metabolism; P-gp efflux
##
## Olaparib:
##   - tmax≈1.5h, t½≈11.9h (300mg BID; Doherty 2014)
##   - Geometric mean Cmax=5.5µM at steady state
##   - CYP3A4 major metabolizer; F≈77%
##
## Niraparib:
##   - tmax≈3h, t½≈36h (QD dosing; Sandhu 2013)
##   - Large Vd (1074 L); CL≈16.2 L/h
##   - Dose reduction to 200mg for BW<77kg or platelets<150k
##
## Bevacizumab:
##   - t½≈20 days (IgG1 antibody; Lu 2008 Cancer Chemother Pharmacol)
##   - CL=0.207 L/day (mainly catabolism/target-mediated)
##   - 15 mg/kg q3w → Cmax≈360 µg/mL
##
## Tumor growth calibration:
##   - Gompertz model: Oza 2015 (ICON7); doubling time ~60 days untreated
##   - CA-125 t½≈23 days (Rustin 1996 JCO)
##   - SOLO-1: olaparib maint. PFS median not reached (vs 13.8mo ctrl)
##   - PRIMA: niraparib PFS 13.8mo (HRD+) vs 8.2mo (ctrl) (Gonzalez-Martin 2019)
##   - PAOLA-1: Ola+Bev PFS 22.1mo vs 16.6mo (HRD+) (Ray-Coquard 2019)
