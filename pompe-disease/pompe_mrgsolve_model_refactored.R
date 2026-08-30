## =================================================================
## Pompe Disease (GSDII) — Quantitative Systems Pharmacology
## mrgsolve ODE model (22 compartments)
##
## Disease: GAA enzyme deficiency → lysosomal glycogen accumulation →
##          autophagic buildup → cardiomyopathy (IOPD) + skeletal myopathy
##          + diaphragmatic / respiratory failure (LOPD)
##
## Drugs implemented (PK + PD coupling):
##   • Alglucosidase alfa (Lumizyme/Myozyme) 20 mg/kg IV q2w
##   • Avalglucosidase alfa (Nexviazyme, COMET trial) 20 mg/kg IV q2w
##   • Cipaglucosidase alfa + Miglustat (Pombiliti+Opfolda, PROPEL)
##   • AAV9-hGAA gene-therapy bolus (research stage)
##   • Rituximab-based ITI (CRIM- IOPD prophylaxis)
##
## Endpoints:
##   • LV mass index (g/m^2)  – IOPD primary
##   • FVC upright (% pred)   – LOPD primary (COMET)
##   • 6-minute walk distance (m)
##   • GMFM-88, ventilator-free survival, anti-GAA ADA titre, Hex4
##
## Parameter values approximate adult (70 kg) LOPD unless noted; IOPD
## variant flagged by IOPD_FLAG = 1. All values illustrative; see
## `pompe_references.md` for source ranges.
##
## [refactor] Per FORK_WORKFLOW_GUIDE.md Part 2 (pluggable PK, named Hill
## interface), the PK/PD blocks for Alglucosidase alfa (ALGLU), Avalglucosidase
## alfa (AVAL), Cipaglucosidase alfa (CIPA) and Miglustat (MIG) were renamed to
## the fork's naming convention (C_<STEM> exposed concentration, EFFECT_<STEM>
## named disease effect). AAV9-hGAA gene therapy and Rituximab ITI are out of
## scope for this refactor and are completely untouched. Full rationale,
## per-compound archetype, and verification results are in
## pompe_refactor_notes.md.
##
## Renamed (values unchanged from the original):
##   ALGLU_C/ALGLU_P                      -> CENT_ALGLU/PERI_ALGLU
##   AVAL_C/AVAL_P                        -> CENT_AVAL/PERI_AVAL
##   CIPA_C/CIPA_P                        -> CENT_CIPA/PERI_CIPA
##   MIG_A/MIG_C, V_MIG                   -> GUT_MIG/CENT_MIG, V1_MIG
##   MIG_STAB/MIG_EC50                    -> EMAX_MIG/EC50_MIG
##   Cp_alglu/Cp_aval/Cp_cipa/Cp_mig       -> C_ALGLU/C_AVAL/C_CIPA/C_MIG
##   Cp_cipa_stab                         -> C_CIPA_STAB
##   uptake_alglu/uptake_aval/uptake_cipa -> EFFECT_ALGLU/EFFECT_AVAL/EFFECT_CIPA
##   mig_eff                              -> EFFECT_MIG
## New (named Hill-interface parameters completing the interface; GAMMA=1
## makes explicit a shape the original already had implicitly as a plain
## ratio): GAMMA_MIG=1 (new $PARAM). EMAX_ALGLU/EC50_ALGLU/GAMMA_ALGLU,
## EMAX_AVAL/EC50_AVAL/GAMMA_AVAL, EMAX_CIPA/EC50_CIPA/GAMMA_CIPA are computed
## in $ODE (not $PARAM) as exact algebraic derivations of the shared,
## pre-existing VMAX_UPT/KM_UPT/RHO_AVAL/RHO_CIPA parameters -- this
## guarantees bit-exact reproduction of the original's arithmetic and avoids
## a second, drifting copy of the same numbers; exposed via $CAPTURE instead.
##
## Pre-existing upstream build defects (unrelated to any compound's own PK,
## both fixed syntax-only here, logged as translations/UPSTREAM_ISSUES.md #71):
##  1. 11 $PARAM @annotated lines had no description field (mrgsolve 2.0.1
##     requires "name : value : description"); a description was added to
##     each (7 belong to AVAL/CIPA's own PK block, folded into this refactor;
##     4 are disease-side DIAPH_LOSS/DIAPH_GAIN/SMWT_MAX/SMWT_MIN, a pure
##     build-compat fix needed for the whole file to compile at all).
##  2. $TABLE ended in 13 bare `capture NAME;` lines (no `$CAPTURE` header,
##     one bare identifier per line). These compile without error under this
##     mrgsolve build but are silently non-functional -- none of the 13 names
##     were actually retrievable as simulation output (confirmed: requesting
##     any of them via /run_simulation fails "not a compartment or captured
##     item"). Replaced with real $CAPTURE header blocks; same names, same
##     values, now actually exposed. The checked-in original
##     (pompe_mrgsolve_model.R) is untouched and still carries both defects.
## =================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

pompe_code <- '
$PROB
# Pompe Disease QSP v1.0
# 22-compartment ODE; adult LOPD baseline, IOPD switchable.
# [refactor] ALGLU/AVAL/CIPA/MIG PK+PD blocks renamed to the fork
# plumbable-PK convention (C_<STEM>/EFFECT_<STEM>); values unchanged from
# the original. AAV/RTX blocks are out of scope and untouched. See
# pompe_refactor_notes.md. Two pre-existing build defects (missing $PARAM
# annotations, non-functional bare capture lines) fixed syntax-only --
# translations/UPSTREAM_ISSUES.md #71.

$PARAM @annotated
// ---- Patient & disease ----
WT        :  70    : Body weight (kg)
BSA       :  1.80  : Body surface area (m^2)
IOPD_FLAG :   0    : 1 = infantile-onset, 0 = LOPD adult
CRIM_NEG  :   0    : 1 = CRIM-negative (higher ADA risk)
GAA_BASE  :   0.10 : Residual GAA activity (fraction of normal)
SEX_M     :   1    : 1 = male, 0 = female (covariate placeholder)

// ---- Alglucosidase alfa PK (2-cmt) ----
KA_ALGLU  :   0.0  : Not used (IV infusion)
CL_ALGLU  :  21    : Clearance L/d (~0.27 mL/min/kg * 70 * 1440)
V1_ALGLU  :   3.5  : Central volume (L) ~50 mL/kg
Q_ALGLU   :   2.0  : Intercompartmental clearance (L/d)
V2_ALGLU  :   4.5  : Peripheral volume (L)
M6P_ALGLU :   2.5  : Mannose-6-P per mole alglucosidase (mol/mol)

// ---- Avalglucosidase alfa PK ----
CL_AVAL   :  17    : Clearance (L/d) lower than alglu
V1_AVAL   :   3.3  : Central V (L)
Q_AVAL    :   2.0  : Intercompartmental clearance (L/d)
V2_AVAL   :   4.5  : Peripheral volume (L)
M6P_AVAL  :  37.5  : ~15-fold higher M6P content

// ---- Cipaglucosidase alfa PK ----
CL_CIPA   :  22    : Clearance (L/d)
V1_CIPA   :   3.6  : Central volume (L)
Q_CIPA    :   1.8  : Intercompartmental clearance (L/d)
V2_CIPA   :   4.5  : Peripheral volume (L)
M6P_CIPA  :  25    : Mannose-6-P content relative to alglucosidase (mol/mol)

// ---- Miglustat (stabiliser) ----  [renamed: V_MIG->V1_MIG, MIG_STAB->EMAX_MIG, MIG_EC50->EC50_MIG]
KA_MIG    :   3.5  : Miglustat absorption (1/d)
CL_MIG    :  84    : Clearance (L/d)
V1_MIG    : 100    : Central volume (L)
F_MIG     :   0.83 : Bioavailability
EMAX_MIG  :   0.35 : Plasma stabilisation factor (Imax) [renamed from MIG_STAB, value unchanged]
EC50_MIG  :   1.5  : Plasma stabilisation EC50 (mg/L) [renamed from MIG_EC50, value unchanged]
GAMMA_MIG :   1.0  : Hill coefficient for stabilisation (new; original had no explicit exponent, implicit gamma=1)

// ---- M6P-receptor uptake (Michaelis–Menten) ----
VMAX_UPT  :   8.0  : Maximal tissue uptake rate (mg/d)
KM_UPT    :   0.7  : Km (mg/L) for CI-MPR
RHO_AVAL  :   3.0  : Relative uptake potency Aval vs Alglu (M6P x15 → uptake gain x3)
RHO_CIPA  :   1.8  : Relative uptake potency Cipa vs Alglu
ADA_KI    :  10    : ADA neutralising IC50 (titre units)
ADA_BLOCK :   0.95 : Maximal ADA-mediated block fraction

// ---- Tissue lysosomal pool ----
KOUT_TISS :   0.07 : Lysosomal enzyme degradation (1/d)
GAA_TURN  :   0.05 : Endogenous GAA recovery rate (1/d)
GLYC_SS_M :  10    : Steady-state lysosomal glycogen muscle (a.u.)
GLYC_SS_C :  15    : Steady-state cardiac glycogen (a.u., IOPD ref)
KIN_GLYC  :   0.20 : Glycogen accumulation rate (1/d)
KOUT_GLYC :   0.04 : Baseline glycogen clearance (1/d)
KGAA_GLYC :   1.20 : GAA-driven glycogen hydrolysis (1/d / nmol)
HX4_GAIN  :   0.30 : Hex4 generation gain
KOUT_HX4  :   0.6  : Urinary Hex4 clearance (1/d)

// ---- ADA dynamics ----
KP_AB     :   0.012 : Plasmablast generation per dose (1/d)
KP_AB_CRIM:   0.045 : CRIM- plasmablast surge (1/d)
KOUT_AB   :   0.05  : ADA decay (1/d, ~14 d t1/2)
ADA_AMP   :   1.0   : Amplification factor on antigen exposure

// ---- Disease physiology ----
MM_BASE   :  20.0   : Muscle mass index baseline (kg)
KM_LOSS   :   0.0012 : Muscle loss rate (1/d) per a.u. glycogen
KM_GAIN   :   0.0006 : Muscle regeneration rate per ERT delivery (1/d)
DIAPH_BASE:   1.0    : Diaphragm function (1 = normal)
DIAPH_LOSS:   0.0010 : Diaphragm function loss rate (1/d)
DIAPH_GAIN:   0.0006 : Diaphragm function gain rate per delivered enzyme (1/d)
FVC_MAX   :  95      : FVC upright max (% predicted)
FVC_MIN   :  20      : Floor
LVMI_BASE :  60      : LV mass index baseline (g/m^2, adult)
LVMI_IOPD : 250      : IOPD presenting LVMI
KLV_GAIN  :   0.05   : LV hypertrophy rate per cardiac glycogen
KLV_REVERSE:  0.04   : Reverse remodel rate per delivered enzyme
SMWT_MAX  : 600    : 6MWT ceiling (m)
SMWT_MIN  : 100    : 6MWT floor (m)

// ---- Endpoint thresholds ----
VENT_FVC_THR : 35    : FVC < 35% triggers vent_failure hazard
VENT_HAZARD  :  0.0008 : Daily hazard for vent_failure if FVC < threshold

// ---- AAV gene therapy ----
AAV_DOSE  :   0      : Single bolus (vg/kg, set externally)
AAV_kexp  :   0.0035 : Expression rise (1/d)
AAV_DECAY :   0.00010: Vector dilution / immune loss (1/d)
AAV_GAIN  :   0.40   : Max contribution to tissue GAA (fraction of normal)
AAV_NAB   :   1.0    : Pre-existing AAV NAb modifier (0–1, 0 blocks)

// ---- Rituximab ITI (for CRIM-) ----
RTX_KIN   :   0.0    : Set externally per dose
RTX_KOUT  :   0.04   : Rituximab decay (1/d, t1/2 ~ 20 d)
RTX_ADA_K :   0.20   : Maximal inhibition of plasmablast generation
RTX_EC50  :  10      : Rituximab EC50 (mg/L)

$CMT @annotated
CENT_ALGLU   : alglucosidase central (mg)
PERI_ALGLU   : alglucosidase peripheral (mg)
CENT_AVAL    : avalglucosidase central (mg)
PERI_AVAL    : avalglucosidase peripheral (mg)
CENT_CIPA    : cipaglucosidase central (mg)
PERI_CIPA    : cipaglucosidase peripheral (mg)
GUT_MIG     : miglustat absorption depot (mg)
CENT_MIG     : miglustat central (mg)
GAA_M     : muscle lysosomal GAA pool (a.u.)
GAA_C     : cardiac lysosomal GAA pool (a.u.)
GAA_D     : diaphragm GAA pool (a.u.)
GLYC_M    : muscle glycogen (a.u.)
GLYC_C    : cardiac glycogen (a.u.)
GLYC_D    : diaphragm glycogen (a.u.)
HEX4      : plasma Hex4 biomarker (a.u.)
ADA_T     : anti-GAA ADA titre (units)
LVMI      : LV mass index (g/m^2)
MM_IDX    : muscle mass index (kg)
DIAPH_F   : diaphragm function (0-1)
FVC_UP    : FVC upright (% predicted)
AAV_X     : AAV-driven tissue GAA expression (a.u.)
RTX_C     : Rituximab central (mg)

$MAIN
// initial conditions tuned for IOPD vs LOPD presentation
double GLYC_init_M  = GLYC_SS_M * (1 - GAA_BASE) + 1e-6;
double GLYC_init_C  = (IOPD_FLAG > 0.5 ? GLYC_SS_C : 2.0) * (1 - GAA_BASE);
double GLYC_init_D  = GLYC_SS_M * (1 - GAA_BASE);
double LVMI_init    = (IOPD_FLAG > 0.5 ? LVMI_IOPD : LVMI_BASE);

GAA_M_0   = GAA_BASE * 10.0;
GAA_C_0   = GAA_BASE * 10.0;
GAA_D_0   = GAA_BASE * 10.0;
GLYC_M_0  = GLYC_init_M;
GLYC_C_0  = GLYC_init_C;
GLYC_D_0  = GLYC_init_D;
HEX4_0    = 8.0 * (1 - GAA_BASE);
ADA_T_0   = 0.0;
LVMI_0    = LVMI_init;
MM_IDX_0  = MM_BASE * (IOPD_FLAG > 0.5 ? 0.4 : 0.85);
DIAPH_F_0 = DIAPH_BASE * (IOPD_FLAG > 0.5 ? 0.30 : 0.80);
FVC_UP_0  = FVC_MAX * (IOPD_FLAG > 0.5 ? 0.35 : 0.65);
AAV_X_0   = 0.0;
RTX_C_0   = 0.0;

$ODE
// ---- Plasma concentrations ----
double C_ALGLU = CENT_ALGLU / V1_ALGLU;
double C_AVAL  = CENT_AVAL  / V1_AVAL;
double C_CIPA  = CENT_CIPA  / V1_CIPA;
double C_MIG   = CENT_MIG   / V1_MIG;
double Cp_rtx   = RTX_C   / V1_ALGLU;   // approx

// Miglustat stabilisation (Imax on Cipa CL) -- named Hill interface, GAMMA_MIG=1
// so pow(x,1)==x exactly (IEEE-754): identical arithmetic to the original
// plain ratio, just made an explicit 3-parameter Hill term per the fork
// naming convention.
double EFFECT_MIG  = EMAX_MIG * pow(C_MIG, GAMMA_MIG) / (pow(EC50_MIG, GAMMA_MIG) + pow(C_MIG, GAMMA_MIG));
double C_CIPA_STAB = C_CIPA * (1 + EFFECT_MIG);

// ADA-mediated neutralisation
double ada_block = ADA_BLOCK * pow(ADA_T,2) / (pow(ADA_KI,2) + pow(ADA_T,2));

// CI-MPR uptake (Michaelis–Menten, sum across drugs) -- named Hill interface.
// EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM> are computed here (not declared as
// independent $PARAM entries) because they are exact algebraic derivations
// of the shared, pre-existing VMAX_UPT/KM_UPT/RHO_AVAL/RHO_CIPA parameters
// -- computing them from those parameters every step guarantees bit-exact
// reproduction of the original arithmetic and avoids a second, drifting
// copy of the same numbers. See pompe_refactor_notes.md.
double EMAX_ALGLU = VMAX_UPT;
double EC50_ALGLU = KM_UPT;
double GAMMA_ALGLU = 1.0;
double EMAX_AVAL  = VMAX_UPT * RHO_AVAL;
double EC50_AVAL  = KM_UPT / RHO_AVAL;
double GAMMA_AVAL = 1.0;
double EMAX_CIPA  = VMAX_UPT * RHO_CIPA;
double EC50_CIPA  = KM_UPT;
double GAMMA_CIPA = 1.0;

double EFFECT_ALGLU = EMAX_ALGLU * pow(C_ALGLU, GAMMA_ALGLU) / (pow(EC50_ALGLU, GAMMA_ALGLU) + pow(C_ALGLU, GAMMA_ALGLU)) * (1 - ada_block);
double EFFECT_AVAL  = EMAX_AVAL  * pow(C_AVAL,  GAMMA_AVAL)  / (pow(EC50_AVAL,  GAMMA_AVAL)  + pow(C_AVAL,  GAMMA_AVAL))  * (1 - 0.5*ada_block);
double EFFECT_CIPA  = EMAX_CIPA  * pow(C_CIPA_STAB, GAMMA_CIPA) / (pow(EC50_CIPA, GAMMA_CIPA) + pow(C_CIPA_STAB, GAMMA_CIPA)) * (1 - 0.5*ada_block);

double tissue_supply = EFFECT_ALGLU + EFFECT_AVAL + EFFECT_CIPA;
double aav_supply    = AAV_GAIN * AAV_X * AAV_NAB;

// Drug PK ODEs (IV infusions/bolus drive in vivo dosing events)
dxdt_CENT_ALGLU = - (CL_ALGLU / V1_ALGLU) * CENT_ALGLU - Q_ALGLU/V1_ALGLU * CENT_ALGLU + Q_ALGLU/V2_ALGLU * PERI_ALGLU;
dxdt_PERI_ALGLU =   Q_ALGLU/V1_ALGLU * CENT_ALGLU - Q_ALGLU/V2_ALGLU * PERI_ALGLU;
dxdt_CENT_AVAL  = - (CL_AVAL / V1_AVAL) * CENT_AVAL - Q_AVAL/V1_AVAL * CENT_AVAL + Q_AVAL/V2_AVAL * PERI_AVAL;
dxdt_PERI_AVAL  =   Q_AVAL/V1_AVAL * CENT_AVAL - Q_AVAL/V2_AVAL * PERI_AVAL;
dxdt_CENT_CIPA  = - (CL_CIPA / V1_CIPA) * CENT_CIPA - Q_CIPA/V1_CIPA * CENT_CIPA + Q_CIPA/V2_CIPA * PERI_CIPA;
dxdt_PERI_CIPA  =   Q_CIPA/V1_CIPA * CENT_CIPA - Q_CIPA/V2_CIPA * PERI_CIPA;
dxdt_GUT_MIG   = - KA_MIG * GUT_MIG;
dxdt_CENT_MIG   =   KA_MIG * GUT_MIG * F_MIG - (CL_MIG / V1_MIG) * CENT_MIG;
dxdt_RTX_C   = - RTX_KOUT * RTX_C;

// Tissue lysosomal enzyme pools (muscle, cardiac, diaphragm)
double rec_GAA   = GAA_TURN * (1 - GAA_BASE);  // small endogenous recovery only if GAA_BASE>0
dxdt_GAA_M   = 0.55 * tissue_supply + 0.6 * aav_supply + rec_GAA - KOUT_TISS * GAA_M;
dxdt_GAA_C   = 0.25 * tissue_supply + 0.3 * aav_supply + rec_GAA - KOUT_TISS * GAA_C;
dxdt_GAA_D   = 0.20 * tissue_supply + 0.1 * aav_supply + rec_GAA - KOUT_TISS * GAA_D;

// Lysosomal glycogen dynamics
dxdt_GLYC_M = KIN_GLYC - (KOUT_GLYC + KGAA_GLYC * GAA_M / 10.0) * GLYC_M;
dxdt_GLYC_C = KIN_GLYC * (IOPD_FLAG > 0.5 ? 1.5 : 0.4)
              - (KOUT_GLYC + KGAA_GLYC * GAA_C / 10.0) * GLYC_C;
dxdt_GLYC_D = KIN_GLYC - (KOUT_GLYC + KGAA_GLYC * GAA_D / 10.0) * GLYC_D;

// Hex4 biomarker (sum proxy)
double hex4_gen = HX4_GAIN * (GLYC_M + 0.5*GLYC_D);
dxdt_HEX4 = hex4_gen - KOUT_HX4 * HEX4;

// ADA dynamics (per-dose Ag triggers via mtime/event)
double ag_drive = (C_ALGLU + C_AVAL + C_CIPA) * ADA_AMP;
double rtx_eff  = RTX_ADA_K * Cp_rtx / (RTX_EC50 + Cp_rtx);
double KP_AB_eff = (CRIM_NEG > 0.5 ? KP_AB_CRIM : KP_AB) * (1 - rtx_eff);
dxdt_ADA_T   = KP_AB_eff * ag_drive - KOUT_AB * ADA_T;

// AAV expression rise then dilution
dxdt_AAV_X   = AAV_kexp * (AAV_DOSE > 0 ? 1.0 - AAV_X : -AAV_X)
               - AAV_DECAY * AAV_X;

// LVMI dynamics (IOPD primary)
dxdt_LVMI = KLV_GAIN * GLYC_C - KLV_REVERSE * (GAA_C + aav_supply);

// Skeletal muscle mass index
dxdt_MM_IDX = -KM_LOSS * GLYC_M * MM_IDX + KM_GAIN * (GAA_M + aav_supply) * (MM_BASE - MM_IDX);

// Diaphragm function (0–1)
dxdt_DIAPH_F = -DIAPH_LOSS * GLYC_D * DIAPH_F + DIAPH_GAIN * (GAA_D + aav_supply) * (1.0 - DIAPH_F);

// FVC upright tracks diaphragm function and muscle mass
double FVC_target = FVC_MIN + (FVC_MAX - FVC_MIN) * (0.7 * DIAPH_F + 0.3 * MM_IDX / MM_BASE);
dxdt_FVC_UP = 0.05 * (FVC_target - FVC_UP);

$TABLE
// 6MWT depends on muscle + diaphragm + FVC
double SMWT = SMWT_MIN + (SMWT_MAX - SMWT_MIN) * (0.5*MM_IDX/MM_BASE + 0.3*DIAPH_F + 0.2*FVC_UP/FVC_MAX);
if (SMWT < SMWT_MIN) SMWT = SMWT_MIN;
if (SMWT > SMWT_MAX) SMWT = SMWT_MAX;

// Composite endpoints
double VENT_RISK = (FVC_UP < VENT_FVC_THR ? VENT_HAZARD * (VENT_FVC_THR - FVC_UP) : 0.0);
double SF36_PCS  = 30 + 25 * (MM_IDX/MM_BASE) + 15 * (FVC_UP/FVC_MAX);
double NTproBNP  = 100 + 25 * LVMI;
double EF_LV     = 65 - 0.05 * fmax(0, LVMI - LVMI_BASE);
double CK        = 200 + 80 * GLYC_M;

$CAPTURE C_ALGLU C_AVAL C_CIPA C_MIG Cp_rtx SMWT VENT_RISK SF36_PCS NTproBNP EF_LV CK tissue_supply ada_block
$CAPTURE EFFECT_ALGLU EFFECT_AVAL EFFECT_CIPA EFFECT_MIG C_CIPA_STAB
$CAPTURE EMAX_ALGLU EC50_ALGLU GAMMA_ALGLU EMAX_AVAL EC50_AVAL GAMMA_AVAL EMAX_CIPA EC50_CIPA GAMMA_CIPA
'

pompe_mod <- mcode("pompe_qsp", pompe_code)

## ----------------------------------------------------------------
## Convenience scenario runner
## ----------------------------------------------------------------

#' Run a Pompe-disease QSP scenario
#' @param scenario  one of c("no_tx","alglu","aval","cipa_mig","aav_gt","alglu_iti")
#' @param iopd      logical; 1 = infantile, 0 = LOPD
#' @param years     simulation duration (years)
#' @return long-format data.frame with simulated outputs
pompe_run <- function(scenario = "alglu", iopd = FALSE, years = 3) {
  if (!requireNamespace("mrgsolve", quietly = TRUE))
    stop("mrgsolve package required")

  param <- list(IOPD_FLAG = as.integer(iopd),
                CRIM_NEG  = as.integer(iopd) * 0.2,  # ~20% IOPD CRIM- assumption
                AAV_DOSE  = 0)

  ev <- mrgsolve::ev()

  q2w <- function(amt_per_kg, days, drug_cmt) {
    n <- floor(days / 14)
    times <- (0:n) * 14
    mrgsolve::ev(time = times, amt = amt_per_kg * 70, rate = -1, cmt = drug_cmt, evid = 1)
  }

  if (scenario == "alglu")     ev <- q2w(20, years*365, "CENT_ALGLU")
  if (scenario == "aval")      ev <- q2w(20, years*365, "CENT_AVAL")
  if (scenario == "cipa_mig") {
    ev_cipa <- q2w(20, years*365, "CENT_CIPA")
    n   <- years*365
    ev_mig <- mrgsolve::ev(time = seq(0, n, 1), amt = 195, cmt = "GUT_MIG", evid = 1)
    ev <- c(ev_cipa, ev_mig)
  }
  if (scenario == "aav_gt") {
    param$AAV_DOSE <- 1
    ev <- mrgsolve::ev(time = 0, amt = 0, cmt = "AAV_X", evid = 0)
  }
  if (scenario == "alglu_iti") {
    ev_a  <- q2w(20, years*365, "CENT_ALGLU")
    ev_rt <- mrgsolve::ev(time = c(-14, -7, 0, 7), amt = 700, cmt = "RTX_C", evid = 1)
    ev    <- c(ev_a, ev_rt)
    param$CRIM_NEG <- 1
  }
  if (scenario == "no_tx") {
    ev <- mrgsolve::ev(time = 0, amt = 0, cmt = "CENT_ALGLU", evid = 0)
  }

  out <- pompe_mod |>
    mrgsolve::param(param) |>
    mrgsolve::ev(ev) |>
    mrgsolve::mrgsim(end = years * 365, delta = 1) |>
    as.data.frame()

  out$scenario <- scenario
  out
}

## ----------------------------------------------------------------
## Calibration notes (anchors)
## ----------------------------------------------------------------
## • Adult LOPD baseline 6MWT ~ 300 m (Kishnani 2019 NEJM COMET)
## • Alglucosidase alfa CL ~ 0.27 mL/min/kg, V_ss ~ 100 mL/kg (Hahn 2008)
## • Avalglucosidase alfa ~3-fold higher M6P-tagged binding & uptake (Pena 2019)
## • Cipaglucosidase+Miglustat (PROPEL): +14 m 6MWT vs alglucosidase (Schoser Lancet Neurol 2021)
## • IOPD untreated mortality <12 mo, ERT brings 1-yr ventilator-free survival to ~88% (Kishnani 2007)
## • CRIM- IOPD HSAT ≥51,200 → loss of clinical response (Banugaria 2011, Messinger 2012 ITI)
## • Hex4 declines on ERT (Young 2009, An 2005)

## Example: out <- pompe_run("alglu", iopd = FALSE, years = 3)
##          ggplot(out, aes(time, FVC_UP)) + geom_line()
