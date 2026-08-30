## ============================================================
## Alcoholic Liver Disease (ALD) — QSP Model (mrgsolve)
## 22 ODE Compartments | 7 Treatment Scenarios
## Mechanistic scope: Ethanol PK → Oxidative Stress → Gut-Liver Axis
##   → Kupffer Cell/NLRP3 → Neutrophil Infiltration
##   → Hepatocyte Death/Regeneration → Fibrosis → Drug PK/PD
## Parameters calibrated to STOPAH (Thursz 2015 NEJM), EASL 2018
##
## ---- REFACTORED sibling of ald_mrgsolve_model.R ----
## Scope of this refactor: G-CSF (GCSF), NAC, and Prednisolone (PRED) only —
## their PK compartments, PK parameters, and PD (Hill) effect interfaces are
## renamed to the fork's pluggable-PK naming convention
## (`FORK_WORKFLOW_GUIDE.md`, Part 2). Pentoxifylline, anakinra, and every
## disease-side equation not driven by one of the three in-scope compounds
## are byte-for-byte copies of the original, aside from a small set of
## disclosed, non-numeric mrgsolve-2.0.1 build-compatibility fixes that were
## required just to get the ORIGINAL file to compile at all (see
## `ald_refactor_notes.md` and `translations/UPSTREAM_ISSUES.md`).
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- mrgsolve model code ----
ald_code <- '
$PARAM @annotated
// --- Ethanol Metabolism ---
k_etoh_abs   : 1.2    : Ethanol absorption rate (gut→blood, /h)
k_etoh_elim  : 0.15   : Zero-order ethanol elimination (BAC units/h)
k_ADH        : 0.9    : ADH-mediated metabolism (/h)
k_CYP2E1_bas : 0.05   : Basal CYP2E1 activity
k_CYP2E1_ind : 0.15   : CYP2E1 induction rate constant (chronic)
km_etoh_CYP  : 50.0   : Km ethanol for CYP2E1 (mg/dL)
k_AA_clear   : 2.0    : Acetaldehyde clearance by ALDH2 (/h)

// --- Oxidative Stress ---
k_ROS_prod   : 0.08   : Basal ROS production (from CYP2E1 activity)
k_ROS_clear  : 0.5    : ROS clearance by GSH/SOD (/h)
GSH0         : 5.0    : Baseline GSH (mM)
k_GSH_synth  : 0.3    : GSH synthesis rate (mM/h)
k_GSH_depl   : 0.08   : GSH depletion by ROS (/h)

// --- Gut-Liver Axis ---
LPS0         : 1.0    : Baseline gut LPS (relative units)
k_LPS_prod   : 0.02   : LPS influx from dysbiosis
k_LPS_clear  : 0.15   : Hepatic LPS clearance (/h)
k_perm_etoh  : 0.005  : Ethanol → intestinal permeability increase

// --- Kupffer Cell Activation ---
KC0          : 1.0    : Baseline Kupffer cell activity (rel)
k_KC_act     : 0.3    : KC activation rate by LPS/AA (/h)
k_KC_res     : 0.05   : KC resolution rate (/h)
EC50_KC_LPS  : 2.0    : LPS EC50 for KC activation

// --- Cytokines ---
TNF0         : 1.0    : Baseline TNF-α (rel)
k_TNF_prod   : 0.5    : TNF production rate by KC
k_TNF_clear  : 0.8    : TNF clearance (/h)
IL1B0        : 1.0    : Baseline IL-1β (rel)
k_IL1B_prod  : 0.4    : IL-1β production (KC/NLRP3)
k_IL1B_clear : 0.7    : IL-1β clearance (/h)

// --- Neutrophil Infiltration ---
NEUT0        : 1.0    : Baseline liver neutrophils (rel)
k_neut_rec   : 0.2    : Neutrophil recruitment rate
k_neut_clear : 0.15   : Neutrophil clearance (/h)
EC50_neut    : 3.0    : CXCL1/8 EC50 for neutrophil recruitment

// --- Hepatocyte Dynamics ---
H0           : 1.0    : Healthy hepatocyte fraction (baseline = 1)
k_Hdeath_TNF : 0.04   : TNF-mediated hepatocyte death rate
k_Hdeath_ROS : 0.03   : ROS-mediated hepatocyte death rate
k_Hdeath_neut: 0.02   : Neutrophil-mediated hepatocyte death
k_Hregen     : 0.005  : Hepatocyte regeneration rate (/h)
k_Hregen_max : 0.02   : Max regeneration rate (GCSF-enhanced)

// --- ALT Kinetics ---
ALT0         : 30.0   : Baseline ALT (IU/L)
k_ALT_rel    : 0.8    : ALT release from injured hepatocytes
k_ALT_clear  : 0.02   : ALT clearance from serum (/h, t½ ~35h)

// --- Bilirubin & INR ---
BILI0        : 1.2    : Baseline bilirubin (mg/dL)
k_bili_prod  : 0.012  : Bilirubin production rate
k_bili_conj  : 0.04   : Hepatic bilirubin conjugation/excretion
INR0         : 1.0    : Baseline INR
k_clot_synth : 0.03   : Clotting factor synthesis rate (/h)
k_clot_clear : 0.025  : Clotting factor clearance (/h)

// --- Fibrosis ---
F0           : 0.3    : Baseline fibrosis score (early-mod ALD)
k_fib        : 0.0003 : Fibrosis progression rate
k_fib_regress: 0.0001 : Fibrosis regression rate
EC50_fib_TGF : 2.0    : TGF-β EC50 for fibrosis progression

// --- Prednisolone (PRED) PK ---
// [REFACTORED: renamed from ka_pred/F_pred/CL_pred/Vc_pred/Q_pred/Vp_pred
// to the fork naming convention; same values, archetype 3 (depot + central
// + peripheral, linear).]
KA_PRED      : 1.5    : Absorption rate (/h)
F_PRED       : 0.82   : Oral bioavailability
CL_PRED      : 8.5    : Clearance (L/h)
V1_PRED      : 35.0   : Central volume (L)
Q_PRED       : 15.0   : Inter-compartment clearance (L/h)
V2_PRED      : 50.0   : Peripheral volume (L)

// --- NAC PK ---
// [REFACTORED: Vc_NAC -> V1_NAC, k_NAC_tissue -> K_TISSUE_NAC (kept as an
// added elimination-like term, not part of the naming table but suffixed
// with the NAC stem per convention); archetype 1 (single compartment,
// linear elimination plus this extra tissue-loss term).]
CL_NAC       : 12.0   : NAC clearance (L/h)
V1_NAC       : 30.0   : NAC central volume (L)
K_TISSUE_NAC : 0.4    : NAC tissue distribution / extra elimination (/h)

// --- G-CSF (GCSF) PK ---
// [REFACTORED: ka_GCSF -> KA_GCSF (kept, still unused -- see refactor notes),
// Vc_GCSF -> V1_GCSF; archetype 1 (single compartment, linear elimination).]
KA_GCSF      : 0.5    : G-CSF SC absorption (/h) [unused in $ODE, see notes]
CL_GCSF      : 0.8    : G-CSF clearance (L/h, receptor-mediated)
V1_GCSF      : 4.5    : G-CSF central volume (L)

// --- PD Effect Parameters ---
// [REFACTORED: Emax_pred/EC50_pred -> EMAX_PRED/EC50_PRED, Emax_NAC/EC50_NAC
// -> EMAX_NAC/EC50_NAC, Emax_GCSF/EC50_GCSF -> EMAX_GCSF/EC50_GCSF; each had
// an implicit linear (Hill n=1) ratio in the original, so GAMMA_*=1 is
// added explicitly per the guide (rename, not a refit). Pentoxifylline and
// anakinra own Emax/EC50 are untouched, out of scope.]
EMAX_PRED    : 0.75   : Max NF-κB inhibition by prednisolone
EC50_PRED    : 150.0  : EC50 prednisolone (ng/mL)
GAMMA_PRED   : 1      : Hill coefficient, prednisolone NF-κB inhibition (fixed; original had none)
EMAX_NAC     : 0.65   : Max GSH restoration by NAC
EC50_NAC     : 80.0   : EC50 NAC plasma (μg/mL)
GAMMA_NAC    : 1      : Hill coefficient, NAC GSH restoration (fixed; original had none)
EMAX_GCSF    : 1.8    : Max neutrophil fold-increase by G-CSF
EC50_GCSF    : 5.0    : EC50 G-CSF (ng/mL)
GAMMA_GCSF   : 1      : Hill coefficient, G-CSF neutrophil stimulation (fixed; original had none)
Emax_pento   : 0.40   : Max TNF inhibition by pentoxifylline
EC50_pento   : 600.0  : EC50 pentoxifylline (ng/mL)
Emax_anakin  : 0.70   : Max IL-1β blockade by anakinra
EC50_anakin  : 2000.0 : EC50 anakinra (ng/mL)

$INIT
ETOH    = 0.0   // Blood ethanol (mg/dL)
AA      = 0.0   // Acetaldehyde (μM)
ROS     = 1.0   // Reactive oxygen species (rel)
GSH     = 5.0   // Glutathione (mM)
LPS     = 1.0   // Gut LPS in portal (rel)
KC      = 1.0   // Kupffer cell activation (rel)
TNF     = 1.0   // TNF-α (rel)
IL1B    = 1.0   // IL-1β (rel)
NEUT    = 1.0   // Liver neutrophils (rel)
H       = 1.0   // Healthy hepatocyte fraction
ALT     = 30.0  // Serum ALT (IU/L)
BILI    = 1.2   // Serum bilirubin (mg/dL)
INR     = 1.0   // INR
F       = 0.3   // Fibrosis score (0–4)
GUT_PRED  = 0.0 // Prednisolone gut [was PRED_gut]
CENT_PRED = 0.0 // Prednisolone central (ng/mL) [was PRED_C]
PERI_PRED = 0.0 // Prednisolone peripheral [was PRED_P]
CENT_NAC  = 0.0 // NAC plasma (μg/mL) [was NAC_C]
CENT_GCSF = 0.0 // G-CSF plasma (ng/mL) [was GCSF_C]
PTX_C   = 0.0   // Pentoxifylline plasma (ng/mL)
ANK_C   = 0.0   // Anakinra plasma (ng/mL)

$GLOBAL
#define MELD_calc (3.78*log(BILI+0.01) + 11.2*log(INR+0.01) + 9.57*log(CREA+0.01) + 6.43)
// CREA fixed at patient baseline; simplified MELD here uses ALT proxy
// Full implementation: link creatinine compartment

$MAIN
// Ethanol input via dosing events (ETOH dose = total g, approximate BAC)

// [REFACTORED: exposed concentrations, one per compound, per the fork
// naming convention. CENT_PRED/CENT_NAC/CENT_GCSF already store the
// PD-facing concentration directly in the original (see refactor notes for
// why the prednisolone dxdt_ terms are not a unit-consistent amount/V
// pair), so C_<STEM> is an identity alias, not a /V1 division -- dividing
// again here would change the numeric value the Hill terms below read.]
double C_PRED = CENT_PRED;   // ng/mL, identity (see refactor notes)
double C_NAC  = CENT_NAC;    // μg/mL, identity
double C_GCSF = CENT_GCSF;   // ng/mL, identity

// Drug PD: NF-κB inhibition by prednisolone
// [REFACTORED: eff_pred_NFkB -> EFFECT_PRED, explicit Hill form
// (GAMMA_PRED=1 makes pow(x,1)=x, i.e. this is a rename, not a refit).]
double EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));

// Drug PD: IL-1β blockade by anakinra [unchanged, out of scope]
double eff_anakin_IL1B = Emax_anakin * ANK_C / (EC50_anakin + ANK_C);

// Drug PD: TNF-α inhibition by pentoxifylline [unchanged, out of scope]
double eff_pento_TNF = Emax_pento * PTX_C / (EC50_pento + PTX_C);

// Drug PD: GSH restoration by NAC
// [REFACTORED: eff_NAC_GSH -> EFFECT_NAC, explicit Hill form (GAMMA_NAC=1,
// rename not a refit).]
double EFFECT_NAC = EMAX_NAC * pow(C_NAC, GAMMA_NAC) / (pow(EC50_NAC, GAMMA_NAC) + pow(C_NAC, GAMMA_NAC));

// Drug PD: G-CSF neutrophil stimulation
// [REFACTORED: the original eff_GCSF_neut = 1.0 + Emax_GCSF*C/(EC50+C) is
// "1 + a plain Hill ratio"; EFFECT_GCSF now holds just the Hill ratio
// (0..EMAX_GCSF) per the guide naming convention, and the "+1" baseline
// fold is re-applied at the one point of use in $ODE (dxdt_NEUT), matching
// "combine only at the point the disease equations actually use them."]
double EFFECT_GCSF = EMAX_GCSF * pow(C_GCSF, GAMMA_GCSF) / (pow(EC50_GCSF, GAMMA_GCSF) + pow(C_GCSF, GAMMA_GCSF));

// Effective ROS clearance (GSH-dependent)
double ROS_clearance = k_ROS_clear * (GSH / GSH0);

// KC activation by LPS and acetaldehyde (NLRP3 signal)
double KC_drive = k_KC_act * LPS / (EC50_KC_LPS + LPS) * (1.0 + 0.3 * AA / 5.0);

// Neutrophil drive by KC-derived CXCL1 (proportional to KC)
double neut_drive = k_neut_rec * KC / (EC50_neut / 3.0 + KC);

// Hepatocyte death (combined TNF, ROS, neutrophil)
double Hdeath = H * (k_Hdeath_TNF * TNF + k_Hdeath_ROS * ROS + k_Hdeath_neut * NEUT);
// Clamp H > 0
double H_safe = (H > 0.01) ? H : 0.01;

// Regeneration rate limited by fibrosis and boosted by G-CSF
// [REFACTORED: (eff_GCSF_neut - 1.0)/Emax_GCSF -> EFFECT_GCSF/EMAX_GCSF,
// an exact algebraic identity: eff_GCSF_neut - 1.0 == EFFECT_GCSF by
// construction above, so this term is numerically unchanged.]
double regen_rate = k_Hregen * (1.0 - F / 4.0) + k_Hregen_max * EFFECT_GCSF / EMAX_GCSF;

// Fibrosis driven by KC/TGF-β (proxy: KC activation)
double TGFb_proxy = KC;
double fib_prog = k_fib * TGFb_proxy / (EC50_fib_TGF + TGFb_proxy) * (1.0 - H_safe);
double fib_regress = k_fib_regress * H_safe;

$ODE
// Ethanol (BAC, mg/dL): absorbed from dose events, eliminated linearly
dxdt_ETOH = -k_etoh_elim * ETOH - k_ADH * ETOH;

// Acetaldehyde: produced by ADH + CYP2E1, cleared by ALDH2
double AA_prod = (k_ADH * ETOH + k_CYP2E1_bas * ETOH * ETOH / (km_etoh_CYP + ETOH));
dxdt_AA = AA_prod - k_AA_clear * AA;

// ROS: produced by CYP2E1 activity, cleared by GSH antioxidant defense
double ROS_prod = k_ROS_prod * (1.0 + k_CYP2E1_ind * ETOH / 50.0) + 0.05 * KC;
dxdt_ROS = ROS_prod - ROS_clearance * ROS;

// Glutathione: synthesized (NAC-boosted), depleted by ROS
double GSH_synth_eff = k_GSH_synth * (1.0 + EFFECT_NAC);
dxdt_GSH = GSH_synth_eff - k_GSH_depl * ROS * GSH - 0.02 * AA * GSH;

// Gut LPS: elevated by ethanol-induced dysbiosis
double LPS_prod = LPS0 * k_LPS_prod * (1.0 + k_perm_etoh * ETOH);
dxdt_LPS = LPS_prod - k_LPS_clear * LPS;

// Kupffer cell activation
double KC_inhib = EFFECT_PRED;  // prednisolone suppresses KC NF-κB
dxdt_KC = KC_drive * (1.0 - KC_inhib) - k_KC_res * KC;

// TNF-α
double TNF_inhib = EFFECT_PRED + eff_pento_TNF - EFFECT_PRED * eff_pento_TNF;
dxdt_TNF = k_TNF_prod * KC * (1.0 - TNF_inhib) - k_TNF_clear * (TNF - TNF0);

// IL-1β (NLRP3-driven; also blocked by anakinra)
double IL1B_inhib = EFFECT_PRED + eff_anakin_IL1B - EFFECT_PRED * eff_anakin_IL1B;
dxdt_IL1B = k_IL1B_prod * KC * (1.0 - IL1B_inhib) - k_IL1B_clear * (IL1B - IL1B0);

// Liver neutrophils
// [REFACTORED: eff_GCSF_neut -> (1.0 + EFFECT_GCSF), same combined value.]
dxdt_NEUT = neut_drive * (1.0 + EFFECT_GCSF) - k_neut_clear * NEUT;

// Healthy hepatocyte fraction
dxdt_H = -Hdeath + regen_rate * (1.0 - H_safe);

// Serum ALT
double ALT_release_rate = k_ALT_rel * Hdeath * ALT0 * 20.0;
dxdt_ALT = ALT_release_rate - k_ALT_clear * (ALT - ALT0);

// Bilirubin: inversely proportional to hepatocyte function
double conj_capacity = k_bili_conj * H_safe;
dxdt_BILI = k_bili_prod - conj_capacity * BILI;

// INR: clotting factor synthesis depends on hepatocyte function
double clot_prod = k_clot_synth * H_safe;
double clot_clear = k_clot_clear;
// INR rises as clotting factor falls (simplified: dINR/dt = base_prod_deficit)
dxdt_INR = (k_clot_synth - clot_prod) / k_clot_synth * 0.01 - 0.002 * (INR - INR0) * H_safe;

// Fibrosis (Laennec score, 0–4)
dxdt_F = fib_prog - fib_regress;

// ---- Drug PK ----
// Prednisolone (2-compartment oral, dose in mg)
// [REFACTORED: PRED_gut/PRED_C/PRED_P -> GUT_PRED/CENT_PRED/PERI_PRED,
// ka_pred/F_pred/CL_pred/Vc_pred/Q_pred/Vp_pred -> KA_PRED/F_PRED/CL_PRED/
// V1_PRED/Q_PRED/V2_PRED; equations byte-identical in structure, renamed
// only (archetype 3).]
dxdt_GUT_PRED  = -KA_PRED * GUT_PRED;
dxdt_CENT_PRED =  KA_PRED * F_PRED * GUT_PRED / V1_PRED * 1000.0
                  - (CL_PRED / V1_PRED) * CENT_PRED
                  - (Q_PRED  / V1_PRED) * CENT_PRED
                  + (Q_PRED  / V2_PRED) * PERI_PRED;
dxdt_PERI_PRED =  (Q_PRED / V1_PRED) * CENT_PRED - (Q_PRED / V2_PRED) * PERI_PRED;

// NAC (IV infusion, 1-compartment simplified)
// [REFACTORED: NAC_C -> CENT_NAC, Vc_NAC -> V1_NAC, k_NAC_tissue -> K_TISSUE_NAC.]
dxdt_CENT_NAC = -(CL_NAC / V1_NAC) * CENT_NAC - K_TISSUE_NAC * CENT_NAC;

// G-CSF (SC, 1-compartment)
// [REFACTORED: GCSF_C -> CENT_GCSF, Vc_GCSF -> V1_GCSF.]
dxdt_CENT_GCSF = -(CL_GCSF / V1_GCSF) * CENT_GCSF;

// Pentoxifylline (1-compartment, oral) [unchanged, out of scope]
dxdt_PTX_C = -(0.6 + 0.5) * PTX_C;   // CL/V = 0.6/h + k12=0.5

// Anakinra (SC, 1-compartment, t½~4h → ke~0.17/h) [unchanged, out of scope]
dxdt_ANK_C = -0.17 * ANK_C;

$TABLE
double MELD  = 3.78*log(BILI + 0.01) + 11.2*log(INR + 0.01) + 9.57*log(1.0 + 0.01) + 6.43;
double DF    = 4.6 * (INR - 1.0) * 14.0 + BILI;  // Maddrey DF (simplified)
double ABIC  = 40.0 * 0.1 + BILI * 0.08 + INR * 0.8 + 1.0 * 0.3; // example ABIC
double logit_d90 = -3.5 + 0.18 * MELD;            // logistic approximation
double prob_d90 = 1.0 / (1.0 + exp(-logit_d90));  // 90-day mortality probability
capture MELD_out  = MELD;
capture DF_out    = DF;
capture prob_d90_out  = prob_d90;
capture ALT_out   = ALT;
capture BILI_out  = BILI;
capture INR_out   = INR;
capture H_out     = H * 100.0;   // percent hepatocytes viable
capture GSH_out   = GSH;
capture ROS_out   = ROS;
capture KC_out    = KC;
capture NEUT_out  = NEUT;
capture F_out     = F;
capture TNF_out   = TNF;
capture IL1B_out  = IL1B;
// [REFACTORED: PRED_C_out/NAC_C_out/GCSF_C_out -> C_PRED_out/C_NAC_out/
// C_GCSF_out (same values, renamed to the exposed-concentration convention),
// plus new EFFECT_*_out captures for qspserver /model_manifest discoverability.
// Recomputed fresh here from the $CMT state (not read back from the $MAIN
// locals C_PRED/EFFECT_PRED/...) because mrgsolve evaluates $MAIN once at
// the START of each interval, before that ODE integration step --
// capturing the $MAIN member directly in $TABLE reports it one output
// interval late (confirmed empirically against this model: with delta=1h
// dosing, $MAIN-sourced C_PRED_out trailed CENT_PRED by exactly one 1h
// step). Recomputing from CENT_PRED/CENT_NAC/CENT_GCSF (the $CMT state,
// current as of the $TABLE evaluation time) avoids that lag; the actual
// disease equations (KC_inhib, TNF_inhib, IL1B_inhib, GSH_synth_eff,
// dxdt_NEUT, regen_rate) still consume the $MAIN-timed EFFECT_PRED/
// EFFECT_NAC/EFFECT_GCSF above, unchanged, exactly matching the original
// $MAIN-then-$ODE evaluation order -- only these *reporting* captures
// are recomputed to avoid the lag, see refactor notes.]
double C_PRED_tbl    = CENT_PRED;
double C_NAC_tbl     = CENT_NAC;
double C_GCSF_tbl    = CENT_GCSF;
double EFFECT_PRED_tbl = EMAX_PRED * pow(C_PRED_tbl, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED_tbl, GAMMA_PRED));
double EFFECT_NAC_tbl  = EMAX_NAC  * pow(C_NAC_tbl,  GAMMA_NAC)  / (pow(EC50_NAC,  GAMMA_NAC)  + pow(C_NAC_tbl,  GAMMA_NAC));
double EFFECT_GCSF_tbl = EMAX_GCSF * pow(C_GCSF_tbl, GAMMA_GCSF) / (pow(EC50_GCSF, GAMMA_GCSF) + pow(C_GCSF_tbl, GAMMA_GCSF));
capture C_PRED_out    = C_PRED_tbl;
capture C_NAC_out     = C_NAC_tbl;
capture C_GCSF_out    = C_GCSF_tbl;
capture EFFECT_PRED_out = EFFECT_PRED_tbl;
capture EFFECT_NAC_out  = EFFECT_NAC_tbl;
capture EFFECT_GCSF_out = EFFECT_GCSF_tbl;
'

## Compile the model
mod <- mcode("ALD_QSP_refactored", ald_code, quiet = TRUE)

## ---- Helper: build event table ----
## [REFACTORED: cmt = "PRED_gut"/"NAC_C"/"GCSF_C" -> "GUT_PRED"/"CENT_NAC"/
## "CENT_GCSF" to match the renamed compartments above. This function is
## dead code in the original file too (never called by build_events() below
## from the executed "Run all 7 scenarios" pipeline further down) -- kept
## dead here as well, renamed only for consistency with the rest of this
## sibling file.]
build_events <- function(scenario, duration_days = 90) {
  evt <- eventd()

  if (scenario %in% c("S1","S2","S3","S4","S5","S6","S7")) {
    # Add daily ethanol exposure for scenarios without abstinence
    if (scenario %in% c("S1")) {
      # S1: active drinking (120 g/day ethanol → BAC ~80 mg/dL pulse)
      for (d in seq(0, duration_days - 1)) {
        evt <- add(evt, evd(time = d * 24, amt = 80, cmt = "ETOH", rate = 4))
      }
    }
  }

  if (scenario == "S2") {
    # S2: Abstinence + supportive care (no active drug, just stops drinking)
    # No ethanol events; natural history with withdrawal
  }

  if (scenario %in% c("S3", "S5")) {
    # Prednisolone 40 mg QD x 28 days (dose in gut compartment, mg units)
    for (d in seq(0, 27)) {
      evt <- add(evt, evd(time = d * 24, amt = 40, cmt = "GUT_PRED"))
    }
  }

  if (scenario %in% c("S4", "S5")) {
    # NAC IV: 150 mg/kg (70kg) = 10500 mg day 1, then 50 mg/kg x 4 days
    # Simplified: IV bolus loading → input to CENT_NAC directly (amt in μg/mL × L = mg)
    evt <- add(evt, evd(time = 0,  amt = 350, cmt = "CENT_NAC", rate = 17.5)) # 20h infusion
    for (d in seq(1, 4)) {
      evt <- add(evt, evd(time = d * 24, amt = 117, cmt = "CENT_NAC", rate = 4.9))
    }
  }

  if (scenario == "S6") {
    # G-CSF 5 μg/kg x 5 days SC (350 μg/day, 70 kg)
    for (d in seq(0, 4)) {
      evt <- add(evt, evd(time = d * 24, amt = 350, cmt = "CENT_GCSF",
                          rate = 350 / Vc_GCSF_dose))
    }
  }

  if (scenario == "S7") {
    # Prednisolone 40 mg x 28 days + Anakinra 100 mg SC daily x 28 days
    for (d in seq(0, 27)) {
      evt <- add(evt, evd(time = d * 24, amt = 40,    cmt = "GUT_PRED"))
      evt <- add(evt, evd(time = d * 24, amt = 100000, cmt = "ANK_C",
                          rate = 100000 / 4))  # 100 mg → 100000 ng over 4h
    }
  }

  evt
}

## ---- Simplified simulation function ----
run_sim <- function(scenario, duration_days = 90) {
  tfinal <- duration_days * 24

  # Set scenario-specific parameters
  params_override <- list()

  if (scenario == "S1") {
    # Active alcohol use (no treatment)
    params_override <- list()
  } else if (scenario == "S2") {
    # Abstinence only — reduce ethanol input, natural recovery
    params_override <- list(k_etoh_elim = 0.30)
  } else if (scenario == "S3") {
    # Prednisolone 40 mg x 28d
    params_override <- list()
  } else if (scenario == "S4") {
    # NAC + supportive
    params_override <- list()
  } else if (scenario == "S5") {
    # Prednisolone + NAC (GET protocol)
    params_override <- list()
  } else if (scenario == "S6") {
    # G-CSF 5 μg/kg x 5 days
    params_override <- list()
  } else if (scenario == "S7") {
    # Prednisolone + Anakinra (investigational)
    params_override <- list()
  }

  m <- do.call(param, c(list(mod), params_override))

  # Initial conditions for severe AH (MELD ~25)
  init_severe <- init(mod,
    ETOH = if (scenario == "S1") 60 else 0,
    AA   = 2.0,
    ROS  = 3.5,
    GSH  = 2.0,  # depleted
    LPS  = 4.0,
    KC   = 3.5,
    TNF  = 4.0,
    IL1B = 5.0,
    NEUT = 3.0,
    H    = 0.55,
    ALT  = 180,
    BILI = 10.0,
    INR  = 1.8,
    F    = 1.5
  )

  m <- init(m, init_severe)

  sim <- mrgsim(m, end = tfinal, delta = 1,
                carry.out = c("time"),
                outvars = c("MELD_out","DF_out","prob_d90_out",
                            "ALT_out","BILI_out","INR_out",
                            "H_out","GSH_out","ROS_out",
                            "KC_out","NEUT_out","F_out",
                            "TNF_out","IL1B_out",
                            "C_PRED_out","C_NAC_out","C_GCSF_out"))

  df <- as.data.frame(sim)
  df$time_days <- df$time / 24
  df$scenario  <- scenario
  df
}

## ---- Run all 7 scenarios ----
scenarios <- list(
  S1 = "Active Drinking (No Rx)",
  S2 = "Abstinence Only",
  S3 = "Prednisolone 40mg x28d",
  S4 = "NAC IV (GET protocol)",
  S5 = "Prednisolone + NAC",
  S6 = "G-CSF 5μg/kg x5d",
  S7 = "Prednisolone + Anakinra"
)

## Simplified simulation (no events for now, use steady parameters)
cat("Running 7 ALD QSP scenarios (refactored PK/PD interface)...\n")

results_list <- list()
for (sc in names(scenarios)) {
  cat(" Scenario:", sc, "-", scenarios[[sc]], "\n")

  # Direct parameter simulation (event-based dosing simplified)
  tfinal <- 90 * 24

  # Parameter sets per scenario
  # [REFACTORED: Emax_pred/Emax_NAC/Emax_GCSF -> EMAX_PRED/EMAX_NAC/EMAX_GCSF;
  # Emax_pento/Emax_anakin/k_etoh_elim unchanged, out of scope.]
  plist <- list(
    S1 = list(EMAX_PRED=0, EMAX_NAC=0, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.05),   # continued drinking
    S2 = list(EMAX_PRED=0, EMAX_NAC=0, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.30),   # abstinence → faster ethanol clearance
    S3 = list(EMAX_PRED=0.75, EMAX_NAC=0, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.20),
    S4 = list(EMAX_PRED=0, EMAX_NAC=0.65, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.20),
    S5 = list(EMAX_PRED=0.75, EMAX_NAC=0.65, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.20),
    S6 = list(EMAX_PRED=0, EMAX_NAC=0, EMAX_GCSF=1.8, Emax_pento=0, Emax_anakin=0,
              k_etoh_elim=0.20),
    S7 = list(EMAX_PRED=0.75, EMAX_NAC=0, EMAX_GCSF=1, Emax_pento=0, Emax_anakin=0.70,
              k_etoh_elim=0.20)
  )

  p_sc <- plist[[sc]]

  # Build steady-state drug concentration proxies
  # Prednisolone steady-state ~200 ng/mL with 40mg QD
  PRED_C_ss  <- if (p_sc$EMAX_PRED > 0) 200 else 0
  NAC_C_ss   <- if (p_sc$EMAX_NAC > 0) 150 else 0
  GCSF_C_ss  <- 0  # pulsatile; effect handled via EMAX_GCSF change
  ANK_C_ss   <- if (p_sc$Emax_anakin > 0) 3500 else 0

  m_sc <- do.call(param, c(list(mod), p_sc))

  m_sc <- init(m_sc,
    AA   = if (sc == "S1") 3.0 else 1.0,
    ROS  = 3.5,
    GSH  = 2.0,
    LPS  = 4.0,
    KC   = 3.5,
    TNF  = 4.0,
    IL1B = 5.0,
    NEUT = 3.0,
    H    = 0.55,
    ALT  = 180,
    BILI = 10.0,
    INR  = 1.8,
    F    = 1.5,
    CENT_PRED = PRED_C_ss,
    CENT_NAC  = NAC_C_ss,
    ANK_C  = ANK_C_ss,
    ETOH   = if (sc == "S1") 60 else 0
  )

  sim <- mrgsim(m_sc, end = tfinal, delta = 4)
  df  <- as.data.frame(sim)
  df$time_days <- df$time / 24
  df$scenario  <- sc
  df$label     <- scenarios[[sc]]
  results_list[[sc]] <- df
}

results <- bind_rows(results_list)

## ---- Plotting ----
sc_colors <- c(S1="#D32F2F", S2="#1976D2", S3="#7B1FA2", S4="#388E3C",
               S5="#F57C00", S6="#00796B", S7="#5D4037")
sc_labels <- unlist(scenarios)

p1 <- ggplot(results, aes(time_days, MELD_out, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey50") +
  annotate("text", x = 85, y = 21, label = "MELD 20", size = 3) +
  labs(x = "Time (days)", y = "MELD Score",
       title = "MELD Score over 90 Days",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p2 <- ggplot(results, aes(time_days, ALT_out, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "grey60") +
  labs(x = "Time (days)", y = "ALT (IU/L)",
       title = "Serum ALT Kinetics",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p3 <- ggplot(results, aes(time_days, BILI_out, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "Bilirubin (mg/dL)",
       title = "Serum Bilirubin",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p4 <- ggplot(results, aes(time_days, H_out, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "Viable Hepatocytes (%)",
       title = "Hepatocyte Viability",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p5 <- ggplot(results, aes(time_days, KC_out, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "KC Activation (rel)",
       title = "Kupffer Cell / Inflammation",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p6 <- ggplot(results, aes(time_days, prob_d90_out * 100, color = scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "Est. 90-day Mortality (%)",
       title = "90-day Mortality Risk",
       color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

panel <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title = "Alcoholic Liver Disease (ALD) — QSP Model (refactored PK/PD)\n7-Scenario Simulation Panel",
    subtitle = "Severe AH baseline: MELD ~24, ALT 180 IU/L, Bilirubin 10 mg/dL",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("ald_simulation_panel_refactored.png", panel, width = 14, height = 14,
       dpi = 150, path = "alcoholic-liver-disease/")
cat("Panel saved: alcoholic-liver-disease/ald_simulation_panel_refactored.png\n")

## ---- GSH & Oxidative Stress sub-analysis ----
results_ox <- results %>%
  select(time_days, scenario, label, GSH_out, ROS_out, TNF_out, IL1B_out) %>%
  pivot_longer(c(GSH_out, ROS_out, TNF_out, IL1B_out),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
    GSH_out  = "Glutathione (mM)",
    ROS_out  = "ROS (rel.)",
    TNF_out  = "TNF-α (rel.)",
    IL1B_out = "IL-1β (rel.)"
  ))

p_ox <- ggplot(results_ox %>% filter(time_days <= 28),
               aes(time_days, value, color = scenario)) +
  geom_line(size = 0.9) +
  facet_wrap(~variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "Level", color = "Scenario",
       title = "Oxidative Stress & Cytokine Dynamics (Day 0–28)",
       subtitle = "First 28 days of treatment window") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E0F7FA"))

ggsave("ald_oxidative_cytokine_refactored.png", p_ox, width = 12, height = 8,
       dpi = 150, path = "alcoholic-liver-disease/")
cat("Saved: ald_oxidative_cytokine_refactored.png\n")

## ---- Fibrosis & Mortality trajectories ----
p_fib_mort <- results %>%
  select(time_days, scenario, label, F_out, prob_d90_out) %>%
  pivot_longer(c(F_out, prob_d90_out), names_to = "var", values_to = "val") %>%
  mutate(var = recode(var,
    F_out    = "Fibrosis Score (Laennec 0–4)",
    prob_d90_out = "90-day Mortality Prob."
  )) %>%
  ggplot(aes(time_days, val, color = scenario)) +
  geom_line(size = 0.9) +
  facet_wrap(~var, scales = "free_y") +
  scale_color_manual(values = sc_colors, labels = sc_labels) +
  labs(x = "Time (days)", y = "", color = "Scenario",
       title = "Fibrosis & Mortality Outcomes",
       subtitle = "90-day follow-up across 7 treatment scenarios") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#FBE9E7"))

ggsave("ald_fibrosis_mortality_refactored.png", p_fib_mort, width = 12, height = 6,
       dpi = 150, path = "alcoholic-liver-disease/")
cat("Saved: ald_fibrosis_mortality_refactored.png\n")

cat("\n=== Simulation Complete ===\n")
cat("Day-90 outcomes by scenario:\n")
results %>%
  filter(time_days >= 89) %>%
  group_by(scenario, label) %>%
  slice_tail(n = 1) %>%
  select(scenario, label, MELD_out, ALT_out, BILI_out, INR_out, F_out, prob_d90_out) %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  print(n = 20)
