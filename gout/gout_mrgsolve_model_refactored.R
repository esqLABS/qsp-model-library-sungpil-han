## ============================================================
## Gout QSP Model — mrgsolve ODE Implementation (REFACTORED)
## Disease: Gout (Hyperuricemia & Crystal-Induced Arthritis)
## Version: 1.0-refactored  Date: 2026-08-30
##
## This is a fork-owned sibling of gout_mrgsolve_model.R, refactored for
## pluggable PK per driver-patches/FORK_WORKFLOW_GUIDE.md Part 2. Only the
## nine compound blocks (ALLO, OXY, FEBU, PROB, LESI, COLCH, INDO, ANA,
## CANA) are renamed/restructured to the guide's naming convention
## (GUT_<STEM>/CENT_<STEM>/PERI_<STEM>/C_<STEM>/EFFECT_<STEM>, etc.).
## The disease-side biology (purine metabolism, urate distribution, renal
## handling, crystal/NLRP3/cytokine cascade, pain, tophus, eGFR) is
## byte-identical to the original. See gout_refactor_notes.md for the
## full archetype-per-compound breakdown, the ALLO/OXY parent-metabolite
## finding, and the qspserver mrgsolve_api verification results.
##
## Calibrated against (same as original):
##   - Becker et al. 2010 (CONFIRMS, febuxostat vs allopurinol)
##   - Sundy et al. 2011 (canakinumab CANTOS)
##   - Terkeltaub et al. 2010 (colchicine AGREE)
##   - Saag et al. 2017 (lesinurad CLEAR studies)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## Model Code (refactored)
## ============================================================
gout_model_code_refactored <- '
$PROB
Gout QSP Model v1.0 (refactored)
Purine metabolism, urate kinetics, crystal inflammation, drug PK/PD
Refactored for pluggable PK per driver-patches/FORK_WORKFLOW_GUIDE.md:
compound blocks (ALLO/OXY/FEBU/PROB/LESI/COLCH/INDO/ANA/CANA) renamed to
the guide'"'"'s naming convention; disease-side biology is byte-identical to
the original.

$PARAM @annotated
// --- Purine Production ---
kprod_UA    : 0.15 : Baseline uric acid production rate (mg/dL/h)
kprod_diet  : 0.03 : Dietary purine contribution (mg/dL/h)
kXO_max     : 1.0  : Maximum XO enzyme capacity (normalized)
km_XO       : 0.5  : Michaelis constant for XO (normalized substrate)
fruct_coeff : 0.2  : Fructose-driven purine production coefficient

// --- Urate Distribution ---
Vc_UA       : 14.0 : Central volume for urate (L)
Vp_UA       : 28.0 : Peripheral volume for urate (L)
Vsyn_UA     : 0.5  : Synovial volume for urate (L)
ktp_UA      : 0.08 : Central to peripheral urate transfer (h^-1)
kpt_UA      : 0.04 : Peripheral to central urate transfer (h^-1)
ksyn_in     : 0.06 : Central to synovial urate transfer (h^-1)
ksyn_out    : 0.05 : Synovial to central urate return (h^-1)

// --- Renal Excretion ---
GFR         : 120  : Glomerular filtration rate (mL/min)
FEurate0    : 0.08 : Baseline fractional excretion of urate
kURAT1      : 0.88 : URAT1 reabsorption fraction (0-1)
kOAT_sec    : 0.12 : OAT1/3 secretion fraction
kABCG2_r    : 0.04 : ABCG2 renal secretion fraction

// --- Gut ABCG2 ---
kgut_ABCG2  : 0.015 : Intestinal ABCG2 urate secretion (h^-1)

// --- Crystal Formation ---
sUA_sat     : 6.8  : Saturation threshold for MSU crystal formation (mg/dL)
kcryst      : 0.002 : Crystal formation rate constant (h^-1 per mg/dL excess)
kdissolve   : 0.001 : Crystal dissolution rate constant (h^-1)
kflare_cryst: 0.005 : Crystal shedding-induced flare rate

// --- Inflammasome (NLRP3) & IL-1b ---
kNLRP3      : 0.8  : NLRP3 activation rate per crystal concentration
kIL1b_prod  : 0.5  : IL-1b production rate (pg/mL/h)
kIL1b_deg   : 0.3  : IL-1b clearance rate (h^-1)
IL1b0       : 1.0  : Baseline IL-1b (pg/mL)

// --- TNF-alpha ---
kTNFa_prod  : 0.3  : TNF-alpha production rate (pg/mL/h)
kTNFa_deg   : 0.4  : TNF-alpha degradation rate (h^-1)
TNFa0       : 5.0  : Baseline TNF-alpha (pg/mL)

// --- Neutrophil Influx ---
kPMN_rec    : 0.2  : Neutrophil recruitment rate
kPMN_deg    : 0.15 : Neutrophil clearance rate (h^-1)
PMN0        : 1.0  : Baseline neutrophil influx (normalized)

// --- Pain & Inflammation Score ---
kpain_IL1b  : 0.04 : Pain contribution from IL-1b
kpain_PGE2  : 0.03 : Pain contribution from PGE2
kpain_res   : 0.1  : Pain resolution rate (h^-1)
pain_max    : 10.0 : Maximum NRS pain score

// --- PK: Allopurinol (GUT_ALLO/CENT_ALLO) -> Oxypurinol (CENT_OXY/PERI_OXY) ---
// Allopurinol has no direct PD effect in this model (matches the
// urolithiasis ALLO/OXP precedent) -- only its metabolite oxypurinol
// inhibits xanthine oxidase. C_ALLO is informational only.
KA_ALLO     : 0.9  : Allopurinol absorption rate constant (h^-1)
F_ALLO      : 0.80 : Allopurinol oral bioavailability
CL_ALLO     : 9.0  : Allopurinol clearance (L/h)
V1_ALLO     : 30.0 : Allopurinol central volume (L)
KCONV_OXY   : 0.5  : Allopurinol to oxypurinol conversion rate (h^-1)

// --- PK: Oxypurinol (allopurinol'"'"'s own active metabolite) ---
CL_OXY      : 0.5  : Oxypurinol renal clearance (L/h)
V1_OXY      : 45.0 : Oxypurinol central volume (L)
V2_OXY      : 30.0 : Oxypurinol peripheral volume (L)
KCP_OXY     : 0.5  : Oxypurinol central-to-peripheral transfer (h^-1, amount-based; not Q/V1/V2-consistent, see notes)
KPC_OXY     : 0.2  : Oxypurinol peripheral-to-central transfer (h^-1, amount-based; not Q/V1/V2-consistent, see notes)
EMAX_OXY    : 1.0  : Max fractional XO inhibition by oxypurinol (new, explicit; original was bare C/(C+Ki))
EC50_OXY    : 0.001: Oxypurinol EC50 (=original Ki_Oxy) for XO inhibition (mg/L)
GAMMA_OXY   : 1.0  : Hill coefficient, oxypurinol XO inhibition (new, explicit; original had no exponent)

// --- PK: Febuxostat (CL/Q/V-consistent 2-compartment; Q_FEBU derived so
//     Q/V1=ktp_Febu=0.3 and Q/V2=kpt_Febu=0.15 both hold exactly) ---
KA_FEBU     : 1.2  : Febuxostat absorption rate constant (h^-1)
F_FEBU      : 0.49 : Febuxostat oral bioavailability
CL_FEBU     : 4.0  : Febuxostat total clearance (L/h)
V1_FEBU     : 12.0 : Febuxostat central volume (L)
V2_FEBU     : 24.0 : Febuxostat peripheral volume (L)
Q_FEBU      : 3.6  : Febuxostat inter-compartmental clearance (L/h; derived = ktp_Febu*Vc_Febu = kpt_Febu*Vp_Febu, both 3.6 -- exact, not fitted)
EMAX_FEBU   : 1.0  : Max fractional XO inhibition by febuxostat (new, explicit)
EC50_FEBU   : 0.000001 : Febuxostat EC50 (=original Ki_Febu) for XO (mg/L, very potent)
GAMMA_FEBU  : 1.0  : Hill coefficient, febuxostat XO inhibition (new, explicit)

// --- PK: Probenecid (CL/Q/V-consistent 2-compartment; Q_PROB derived,
//     Q/V1=ktp_Prob=0.4 and Q/V2=kpt_Prob=0.2 both hold exactly, =4.0) ---
KA_PROB     : 1.5  : Probenecid absorption rate (h^-1)
F_PROB      : 1.0  : Probenecid bioavailability
CL_PROB     : 3.5  : Probenecid clearance (L/h)
V1_PROB     : 10.0 : Probenecid central volume (L)
V2_PROB     : 20.0 : Probenecid peripheral volume (L)
Q_PROB      : 4.0  : Probenecid inter-compartmental clearance (L/h; derived = ktp_Prob*Vc_Prob = kpt_Prob*Vp_Prob, both 4.0 -- exact, not fitted)
EMAX_PROB   : 1.0  : Max fractional URAT1 inhibition by probenecid (new, explicit; original ratio form C/IC50/(1+C/IC50) is algebraically C/(C+IC50))
EC50_PROB   : 5.0  : Probenecid EC50 (=original IC50_Prob) for URAT1 (mg/L)
GAMMA_PROB  : 1.0  : Hill coefficient, probenecid URAT1 inhibition (new, explicit)

// --- PK: Lesinurad (no peripheral compartment in the original) ---
KA_LESI     : 2.1  : Lesinurad absorption rate (h^-1)
F_LESI      : 0.95 : Lesinurad bioavailability
CL_LESI     : 8.0  : Lesinurad clearance (L/h)
V1_LESI     : 14.0 : Lesinurad central volume (L)
EMAX_LESI   : 1.0  : Max fractional URAT1 inhibition by lesinurad (new, explicit)
EC50_LESI   : 0.1  : Lesinurad EC50 (=original IC50_Lesi) for URAT1 (mg/L)
GAMMA_LESI  : 1.0  : Hill coefficient, lesinurad URAT1 inhibition (new, explicit)

// --- PK: Colchicine (amount-based transfer constants; NOT Q/V1/V2-
//     consistent -- ktp_Colch*Vc_Colch=210 != kpt_Colch*Vp_Colch=168, see notes) ---
KA_COLCH    : 1.8  : Colchicine absorption rate (h^-1)
F_COLCH     : 0.45 : Colchicine oral bioavailability
CL_COLCH    : 20.0 : Colchicine total clearance (L/h)
V1_COLCH    : 100.0: Colchicine central volume (L)
V2_COLCH    : 400.0: Colchicine peripheral volume (L)
KCP_COLCH   : 2.1  : Colchicine central-to-peripheral transfer (h^-1, amount-based)
KPC_COLCH   : 0.42 : Colchicine peripheral-to-central transfer (h^-1, amount-based)
EMAX_COLCH  : 0.85 : Maximum NLRP3 inhibition by colchicine (=original Emax_Colch, rename only)
EC50_COLCH  : 0.0003: Colchicine EC50 (=original IC50_Colch) for NLRP3 (mg/L)
GAMMA_COLCH : 1.0  : Hill coefficient, colchicine NLRP3 inhibition (new, explicit)

// --- PK: Indomethacin (NSAID); PD effect reads the tissue/synovial
//     compartment (C_INDO_PERI), not plasma (C_INDO, kept informational),
//     same "PD reads the non-plasma site" pattern as
//     abdominal-aortic-aneurysm/Doxycycline. Amount-based transfer, NOT
//     Q/V1/V2-consistent (ktp_Indo*Vc_Indo=9.6 != kpt_Indo*Vt_Indo=6) ---
KA_INDO     : 1.5  : Indomethacin absorption rate (h^-1)
F_INDO      : 0.90 : Indomethacin bioavailability
CL_INDO     : 8.0  : Indomethacin clearance (L/h)
V1_INDO     : 16.0 : Indomethacin central volume (L)
V2_INDO     : 20.0 : Indomethacin tissue/synovial volume (L) (=original Vt_Indo)
KCP_INDO    : 0.6  : Indomethacin central-to-tissue transfer (h^-1, amount-based)
KPC_INDO    : 0.3  : Indomethacin tissue-to-central transfer (h^-1, amount-based)
EMAX_INDO   : 1.0  : Max fractional COX inhibition by indomethacin (new, explicit)
EC50_INDO   : 0.002: Indomethacin EC50 (=original IC50_Indo) for COX (mg/L)
GAMMA_INDO  : 1.0  : Hill coefficient, indomethacin COX inhibition (new, explicit)

// --- PK: Anakinra (IL-1Ra); SC depot named GUT_ANA per convention
//     (absorption depot role, not necessarily oral). Amount-based
//     transfer, NOT Q/V1/V2-consistent (ktp_Ana*Vc_Ana=1.2 != kpt_Ana*Vp_Ana=0.4) ---
KA_ANA      : 0.6  : Anakinra SC absorption rate (h^-1)
F_ANA       : 0.95 : Anakinra SC bioavailability
CL_ANA      : 8.0  : Anakinra total clearance (L/h)
V1_ANA      : 8.0  : Anakinra central volume (L)
V2_ANA      : 5.0  : Anakinra peripheral volume (L)
KCP_ANA     : 0.15 : Anakinra central-to-peripheral transfer (h^-1, amount-based)
KPC_ANA     : 0.08 : Anakinra peripheral-to-central transfer (h^-1, amount-based)
EMAX_ANA    : 1.0  : Max fractional IL-1R blockade by anakinra (new, explicit)
EC50_ANA    : 0.5  : Anakinra EC50 (=original IC50_Ana) for IL-1R blockade (mg/L)
GAMMA_ANA   : 1.0  : Hill coefficient, anakinra IL-1R blockade (new, explicit)

// --- PK: Canakinumab (SC, mass-action neutralization of the disease'"'"'s
//     own IL-1b state, not a canakinumab-owned receptor pool -- kept as
//     genuine drug-target binding per the guide'"'"'s TMDD guidance, same
//     "shared disease state stays unrenamed" treatment as VEGF_FREE in
//     age-related-macular-degeneration/amd_refactor_notes.md).
//     No algebraic EFFECT_CANA exists in the original -- neutralization
//     is expressed entirely by the mass-action ODE below, so
//     EFFECT_CANA below is a diagnostic-only fractional-occupancy
//     readout, not fed back into any disease equation (matches the
//     "no fit possible/attempted, multi-step ODE chain kept unflattened"
//     precedent for EFFECT_VIT in the same amd notes). Amount-based
//     transfer, NOT Q/V1/V2-consistent. ---
KA_CANA     : 0.009: Canakinumab SC absorption rate (h^-1)
F_CANA      : 0.70 : Canakinumab SC bioavailability
CL_CANA     : 0.23 : Canakinumab clearance (L/h)
V1_CANA     : 6.0  : Canakinumab central volume (L)
V2_CANA     : 3.5  : Canakinumab peripheral volume (L)
KCP_CANA    : 0.005: Canakinumab central-to-peripheral transfer (h^-1, amount-based)
KPC_CANA    : 0.003: Canakinumab peripheral-to-central transfer (h^-1, amount-based)
KON_CANA    : 1.0  : Canakinumab-IL1b association (L/mg/h)
KOFF_CANA   : 0.00001: Canakinumab-IL1b dissociation (h^-1)

// --- Tophus Dynamics ---
ktoph_form  : 0.0001: Tophus crystal deposition rate (cm3/h per unit crystal)
ktoph_diss  : 0.00005: Tophus dissolution rate (cm3/h per unit sUA lowering)
Toph0       : 0.0  : Initial tophus volume (cm3, 0 for non-tophaceous)

// --- Disease Progression ---
k_joint_dmg : 0.0001: Cumulative joint damage rate (units/h per flare)
k_eGFR_loss : 0.00002: eGFR decline per unit chronic hyperuricemia (h^-1)
eGFR0       : 90.0 : Baseline eGFR (mL/min/1.73m2)

// --- Covariates (patient-specific) ---
BW          : 80.0 : Body weight (kg)
AGE         : 50.0 : Age (years)
SEX         : 1.0  : Sex (1=male, 0=female)
RACE_AFRO   : 0.0  : African ancestry (1=yes; up risk)
FOOD_score  : 0.5  : Diet purine score (0-1)
ETOH        : 0.0  : Alcohol intake (drinks/day)

$CMT @annotated
// Purine/Urate
A_UA_gut    : Urate precursor gut absorption depot
A_UA_cent   : Urate central compartment (plasma, mg)
A_UA_peri   : Urate peripheral tissue compartment (mg)
A_UA_syn    : Urate in synovial fluid (mg)
A_Crystal   : MSU crystal pool in joint (mg, normalized)
A_Tophus    : Tophus volume (cm3)

// Inflammation
A_IL1b      : Active IL-1b concentration (pg/mL)
A_TNFa      : TNF-alpha (pg/mL)
A_PMN       : Neutrophil influx to joint (normalized)
A_Pain      : Acute gout pain score (NRS 0-10)
A_JointDmg  : Cumulative joint damage (arbitrary units)
A_eGFR      : eGFR (mL/min/1.73m2)

// PK: Allopurinol / Oxypurinol
GUT_ALLO    : Allopurinol gut depot
CENT_ALLO   : Allopurinol central
CENT_OXY    : Oxypurinol central (allopurinol'"'"'s own active metabolite)
PERI_OXY    : Oxypurinol peripheral

// PK: Febuxostat
GUT_FEBU    : Febuxostat gut depot
CENT_FEBU   : Febuxostat central
PERI_FEBU   : Febuxostat peripheral

// PK: Probenecid
GUT_PROB    : Probenecid gut depot
CENT_PROB   : Probenecid central
PERI_PROB   : Probenecid peripheral

// PK: Lesinurad
GUT_LESI    : Lesinurad gut depot
CENT_LESI   : Lesinurad central

// PK: Colchicine
GUT_COLCH   : Colchicine gut depot
CENT_COLCH  : Colchicine central
PERI_COLCH  : Colchicine peripheral

// PK: Indomethacin
GUT_INDO    : Indomethacin gut depot
CENT_INDO   : Indomethacin central
PERI_INDO   : Indomethacin tissue/synovial (PD-relevant site)

// PK: Anakinra
GUT_ANA     : Anakinra SC depot
CENT_ANA    : Anakinra central
PERI_ANA    : Anakinra peripheral

// PK: Canakinumab
GUT_CANA    : Canakinumab SC depot
CENT_CANA   : Canakinumab central
PERI_CANA   : Canakinumab peripheral
COMPLEX_CANA: IL-1b:Canakinumab complex (canakinumab-owned; A_IL1b itself is the shared disease state, kept unrenamed)

$MAIN
// Disease-side concentrations (unchanged)
double C_UA     = A_UA_cent / Vc_UA;         // serum urate (mg/dL)
double C_UA_syn = A_UA_syn / Vsyn_UA;        // synovial urate (mg/dL)
double C_Crystal = A_Crystal;                // crystal concentration

// Exposed compound concentrations (single point an external covariate
// could later substitute in for each compound)
double C_ALLO   = CENT_ALLO / V1_ALLO;      // allopurinol (mg/L) -- informational only, no direct PD effect
double C_OXY    = CENT_OXY  / V1_OXY;       // oxypurinol (mg/L)
double C_FEBU   = CENT_FEBU / V1_FEBU;      // febuxostat (mg/L)
double C_PROB   = CENT_PROB / V1_PROB;      // probenecid (mg/L)
double C_LESI   = CENT_LESI / V1_LESI;      // lesinurad (mg/L)
double C_COLCH  = CENT_COLCH / V1_COLCH;    // colchicine (mg/L)
double C_INDO   = CENT_INDO / V1_INDO;      // indomethacin plasma (mg/L) -- informational only, PD reads tissue below
double C_INDO_PERI = PERI_INDO / V2_INDO;   // indomethacin tissue/synovial (mg/L) -- the site COX inhibition actually reads
double C_ANA    = CENT_ANA / V1_ANA;        // anakinra (mg/L)
double C_CANA   = CENT_CANA / V1_CANA;      // canakinumab (mg/L)

// XO inhibition (combined oxypurinol + febuxostat) -- named Hill effect
// per compound, combined only at the point the disease equation (XO_activity) uses them
double EFFECT_OXY  = EMAX_OXY  * pow(C_OXY,  GAMMA_OXY)  / (pow(EC50_OXY,  GAMMA_OXY)  + pow(C_OXY,  GAMMA_OXY));
double EFFECT_FEBU = EMAX_FEBU * pow(C_FEBU, GAMMA_FEBU) / (pow(EC50_FEBU, GAMMA_FEBU) + pow(C_FEBU, GAMMA_FEBU));
double XO_inhib      = 1.0 - (1.0 - EFFECT_OXY) * (1.0 - EFFECT_FEBU);
double XO_activity   = kXO_max * (1.0 - XO_inhib);

// URAT1 inhibition (probenecid + lesinurad)
double EFFECT_PROB = EMAX_PROB * pow(C_PROB, GAMMA_PROB) / (pow(EC50_PROB, GAMMA_PROB) + pow(C_PROB, GAMMA_PROB));
double EFFECT_LESI = EMAX_LESI * pow(C_LESI, GAMMA_LESI) / (pow(EC50_LESI, GAMMA_LESI) + pow(C_LESI, GAMMA_LESI));
double URAT1_inhib   = 1.0 - (1.0 - EFFECT_PROB) * (1.0 - EFFECT_LESI);
double kURAT1_eff    = kURAT1 * (1.0 - URAT1_inhib);

// Renal clearance of urate
double CLr_UA = GFR * 0.001 * 60.0 * (1.0 - kURAT1_eff + kOAT_sec + kABCG2_r); // L/h

// Gut excretion (ABCG2)
double ABCG2_Q141K_effect = 1.0; // set < 1.0 for Q141K polymorphism patients
double kgut_eff = kgut_ABCG2 * ABCG2_Q141K_effect;

// NLRP3 inhibition by colchicine
double EFFECT_COLCH = EMAX_COLCH * pow(C_COLCH, GAMMA_COLCH) / (pow(EC50_COLCH, GAMMA_COLCH) + pow(C_COLCH, GAMMA_COLCH));

// COX inhibition by NSAID (reads the tissue/synovial concentration)
double EFFECT_INDO = EMAX_INDO * pow(C_INDO_PERI, GAMMA_INDO) / (pow(EC50_INDO, GAMMA_INDO) + pow(C_INDO_PERI, GAMMA_INDO));

// IL-1b effective concentration (accounting for Cana neutralization)
double IL1b_free  = A_IL1b;
double IL1b_bound = COMPLEX_CANA;

// Canakinumab: diagnostic-only fractional IL-1b neutralization (NOT fed
// back into any disease equation -- the real effect is already fully
// expressed by the mass-action dxdt_A_IL1b/dxdt_COMPLEX_CANA terms below;
// no single-line Hill effect exists in the original to rename or fit)
double EFFECT_CANA = COMPLEX_CANA / (COMPLEX_CANA + A_IL1b + 1.0e-9);

// IL-1R blockade by anakinra
double EFFECT_ANA = EMAX_ANA * pow(C_ANA, GAMMA_ANA) / (pow(EC50_ANA, GAMMA_ANA) + pow(C_ANA, GAMMA_ANA));

// PGE2 proxy (proportional to COX activity * inflammatory stimulus)
double PGE2 = (1.0 - EFFECT_INDO) * (IL1b_free / (IL1b_free + 10.0));

// Crystal formation: driven by supersaturation in synovial fluid
double dCryst_form = (C_UA_syn > sUA_sat) ?
    kcryst * (C_UA_syn - sUA_sat) : 0.0;
double dCryst_diss = kdissolve * A_Crystal * (1.0 / (1.0 + C_UA_syn / sUA_sat));

// Flare trigger (crystal concentration drives NLRP3)
double NLRP3_act = kNLRP3 * A_Crystal / (1.0 + A_Crystal) * (1.0 - EFFECT_COLCH);

// IL-1b production with IL-1Ra (anakinra) effect
double IL1b_prod = kIL1b_prod * NLRP3_act * (1.0 - EFFECT_ANA);

// Uric acid production rate (diet + endogenous + alcohol + fructose)
double UA_prod = (kprod_UA + kprod_diet * (1.0 + FOOD_score + ETOH * fruct_coeff))
                 * XO_activity * (BW / 70.0);

// Pain: driven by IL-1b, PGE2, PMN activation
double pain_drive = kpain_IL1b * IL1b_free + kpain_PGE2 * PGE2 * 10.0
                    + 0.02 * A_PMN;

// Initial conditions
A_UA_cent_0  = 6.0 * Vc_UA;    // sUA ~6 mg/dL at baseline
A_UA_peri_0  = 6.0 * Vp_UA * 0.8;
A_UA_syn_0   = 5.5 * Vsyn_UA;
A_IL1b_0     = IL1b0;
A_TNFa_0     = TNFa0;
A_PMN_0      = PMN0;
A_Pain_0     = 0.0;
A_eGFR_0     = eGFR0;
A_Tophus_0   = Toph0;

$ODE
// =============================================================
// Urate Distribution
// =============================================================
dxdt_A_UA_gut  = 0.0; // Fed from dosing events (not dietary)
dxdt_A_UA_cent = UA_prod                                 // production
               - CLr_UA * C_UA                           // renal excretion
               - kgut_eff * A_UA_cent                    // gut secretion
               - ktp_UA * A_UA_cent                      // to peripheral
               + kpt_UA * A_UA_peri                      // from peripheral
               - ksyn_in * A_UA_cent                     // to synovial
               + ksyn_out * A_UA_syn;                    // from synovial

dxdt_A_UA_peri = ktp_UA * A_UA_cent - kpt_UA * A_UA_peri;

dxdt_A_UA_syn  = ksyn_in * A_UA_cent
               - ksyn_out * A_UA_syn
               - dCryst_form * Vsyn_UA                   // crystal deposition removes UA
               + dCryst_diss * Vsyn_UA;                  // dissolution releases UA

// =============================================================
// Crystal & Tophus
// =============================================================
dxdt_A_Crystal = dCryst_form - dCryst_diss
               - kflare_cryst * A_Crystal * A_PMN;       // crystal clearance by PMN

dxdt_A_Tophus  = ktoph_form * A_Crystal
               - ktoph_diss * A_Tophus * fmax(0.0, 6.0 - C_UA); // dissolve if sUA<6

// =============================================================
// NLRP3 Inflammasome & Cytokines
// =============================================================
dxdt_A_IL1b    = IL1b_prod
               - kIL1b_deg * A_IL1b
               - KON_CANA * A_IL1b * CENT_CANA + KOFF_CANA * COMPLEX_CANA;

dxdt_A_TNFa    = kTNFa_prod * NLRP3_act * (1.0 - 0.8 * EFFECT_INDO)
               - kTNFa_deg * A_TNFa;

dxdt_A_PMN     = kPMN_rec * (IL1b_free / (IL1b_free + 5.0) + 0.5 * A_TNFa / (A_TNFa + 20.0))
               - kPMN_deg * A_PMN;

// =============================================================
// Pain & Disease Outcomes
// =============================================================
dxdt_A_Pain    = pain_drive - kpain_res * A_Pain;

dxdt_A_JointDmg = k_joint_dmg * A_PMN * A_Crystal;

dxdt_A_eGFR    = -k_eGFR_loss * fmax(0.0, C_UA - 6.0) * A_eGFR;

// =============================================================
// PK: Allopurinol / Oxypurinol
// =============================================================
dxdt_GUT_ALLO  = -KA_ALLO * GUT_ALLO;
dxdt_CENT_ALLO = KA_ALLO * F_ALLO * GUT_ALLO
               - (CL_ALLO / V1_ALLO) * CENT_ALLO
               - KCONV_OXY * CENT_ALLO;
dxdt_CENT_OXY  = KCONV_OXY * CENT_ALLO
               - (CL_OXY / V1_OXY) * CENT_OXY
               - KCP_OXY * CENT_OXY
               + KPC_OXY * PERI_OXY;
dxdt_PERI_OXY  = KCP_OXY * CENT_OXY - KPC_OXY * PERI_OXY;

// =============================================================
// PK: Febuxostat
// =============================================================
dxdt_GUT_FEBU  = -KA_FEBU * GUT_FEBU;
dxdt_CENT_FEBU = KA_FEBU * F_FEBU * GUT_FEBU
               - (CL_FEBU + Q_FEBU) / V1_FEBU * CENT_FEBU
               + Q_FEBU / V2_FEBU * PERI_FEBU;
dxdt_PERI_FEBU = Q_FEBU / V1_FEBU * CENT_FEBU - Q_FEBU / V2_FEBU * PERI_FEBU;

// =============================================================
// PK: Probenecid
// =============================================================
dxdt_GUT_PROB  = -KA_PROB * GUT_PROB;
dxdt_CENT_PROB = KA_PROB * F_PROB * GUT_PROB
               - (CL_PROB + Q_PROB) / V1_PROB * CENT_PROB
               + Q_PROB / V2_PROB * PERI_PROB;
dxdt_PERI_PROB = Q_PROB / V1_PROB * CENT_PROB - Q_PROB / V2_PROB * PERI_PROB;

// =============================================================
// PK: Lesinurad
// =============================================================
dxdt_GUT_LESI  = -KA_LESI * GUT_LESI;
dxdt_CENT_LESI = KA_LESI * F_LESI * GUT_LESI
               - (CL_LESI / V1_LESI) * CENT_LESI;

// =============================================================
// PK: Colchicine
// =============================================================
dxdt_GUT_COLCH  = -KA_COLCH * GUT_COLCH;
dxdt_CENT_COLCH = KA_COLCH * F_COLCH * GUT_COLCH
                - (CL_COLCH / V1_COLCH) * CENT_COLCH
                - KCP_COLCH * CENT_COLCH
                + KPC_COLCH * PERI_COLCH;
dxdt_PERI_COLCH = KCP_COLCH * CENT_COLCH - KPC_COLCH * PERI_COLCH;

// =============================================================
// PK: Indomethacin
// =============================================================
dxdt_GUT_INDO  = -KA_INDO * GUT_INDO;
dxdt_CENT_INDO = KA_INDO * F_INDO * GUT_INDO
               - (CL_INDO / V1_INDO) * CENT_INDO
               - KCP_INDO * CENT_INDO
               + KPC_INDO * PERI_INDO;
dxdt_PERI_INDO = KCP_INDO * CENT_INDO - KPC_INDO * PERI_INDO;

// =============================================================
// PK: Anakinra (SC)
// =============================================================
dxdt_GUT_ANA   = -KA_ANA * GUT_ANA;
dxdt_CENT_ANA  = KA_ANA * F_ANA * GUT_ANA
               - (CL_ANA / V1_ANA) * CENT_ANA
               - KCP_ANA * CENT_ANA
               + KPC_ANA * PERI_ANA;
dxdt_PERI_ANA  = KCP_ANA * CENT_ANA - KPC_ANA * PERI_ANA;

// =============================================================
// PK: Canakinumab (SC, target-mediated disposition -- binds the shared
// disease IL-1b state directly on amount, not concentration, matching
// the original'"'"'s own dxdt_A_IL1b/dxdt_A_IL1b_Cana terms exactly)
// =============================================================
dxdt_GUT_CANA    = -KA_CANA * GUT_CANA;
dxdt_CENT_CANA   = KA_CANA * F_CANA * GUT_CANA
                 - (CL_CANA / V1_CANA) * CENT_CANA
                 - KCP_CANA * CENT_CANA
                 + KPC_CANA * PERI_CANA
                 - KON_CANA * CENT_CANA * A_IL1b
                 + KOFF_CANA * COMPLEX_CANA;
dxdt_PERI_CANA   = KCP_CANA * CENT_CANA - KPC_CANA * PERI_CANA;
dxdt_COMPLEX_CANA = KON_CANA * CENT_CANA * A_IL1b
                 - KOFF_CANA * COMPLEX_CANA
                 - (CL_CANA / V1_CANA) * COMPLEX_CANA; // catabolism of complex

$TABLE
capture sUA      = A_UA_cent / Vc_UA;
capture sUA_syn  = A_UA_syn / Vsyn_UA;
capture Crystal  = A_Crystal;
capture Tophus   = A_Tophus;
capture IL1b_f   = A_IL1b;
capture TNFa_f   = A_TNFa;
capture PMN      = A_PMN;
capture Pain     = fmin(A_Pain, pain_max);
capture JntDmg   = A_JointDmg;
capture eGFR_sim = A_eGFR;
capture XO_inh   = XO_inhib * 100.0;  // % XO inhibition
capture URAT1_inh = URAT1_inhib * 100.0;
capture C_ALLO_OUT = C_ALLO;
capture C_OXY_OUT  = C_OXY;
capture C_FEBU_OUT = C_FEBU;
capture C_PROB_OUT = C_PROB;
capture C_LESI_OUT = C_LESI;
capture C_COLCH_OUT = C_COLCH;
capture C_INDO_OUT  = C_INDO;
capture C_INDO_PERI_OUT = C_INDO_PERI;
capture C_ANA_OUT  = C_ANA;
capture C_CANA_OUT = C_CANA;
capture EFFECT_OXY_OUT   = EFFECT_OXY;
capture EFFECT_FEBU_OUT  = EFFECT_FEBU;
capture EFFECT_PROB_OUT  = EFFECT_PROB;
capture EFFECT_LESI_OUT  = EFFECT_LESI;
capture EFFECT_COLCH_OUT = EFFECT_COLCH;
capture EFFECT_INDO_OUT  = EFFECT_INDO;
capture EFFECT_ANA_OUT   = EFFECT_ANA;
capture EFFECT_CANA_OUT  = EFFECT_CANA;
capture FEurate    = CLr_UA * C_UA / (GFR * 0.001 * 60.0 * C_UA + 0.0001) * 100.0;
'

## Compile model (distinct model name to avoid a soloc cache collision
## with the original if both are ever compiled in the same R session)
mod <- mread_cache("gout_qsp_refactored", tempdir(), gout_model_code_refactored)

## ============================================================
## SCENARIO DEFINITIONS
## ============================================================

## Helper: build dosing event
dose_ev <- function(drug_cmt, amount, ii, addl, start=0) {
    ev(cmt=drug_cmt, amt=amount, ii=ii, addl=addl, time=start, rate=0)
}

# Simulation time: 52 weeks (8736 h)
sim_end  <- 52 * 7 * 24   # hours
obs_times <- seq(0, sim_end, by=24)

## ============================================================
## Scenario 1: Untreated Hyperuricemia (baseline)
## ============================================================
scen1 <- mod %>%
    param(FOOD_score=0.7, ETOH=2.0) %>%   # high-purine diet + alcohol
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="1_Untreated")

## ============================================================
## Scenario 2: Allopurinol 300 mg/day (standard urate-lowering)
## ============================================================
e2 <- ev(cmt="GUT_ALLO", amt=300, ii=24, addl=sim_end/24 - 1, time=0)

scen2 <- mod %>%
    ev(e2) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="2_Allopurinol300")

## ============================================================
## Scenario 3: Febuxostat 80 mg/day (non-purine XO inhibitor)
## ============================================================
e3 <- ev(cmt="GUT_FEBU", amt=80, ii=24, addl=sim_end/24 - 1, time=0)

scen3 <- mod %>%
    ev(e3) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="3_Febuxostat80")

## ============================================================
## Scenario 4: Combination — Allopurinol + Lesinurad
##   (suboptimal responders needing dual therapy)
## ============================================================
e4a <- ev(cmt="GUT_ALLO",  amt=300, ii=24, addl=sim_end/24 - 1, time=0)
e4b <- ev(cmt="GUT_LESI",  amt=200, ii=24, addl=sim_end/24 - 1, time=0)

scen4 <- mod %>%
    ev(e4a + e4b) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="4_Allo_Lesinurad")

## ============================================================
## Scenario 5: Acute Gout Flare Treatment — Colchicine
##   Low-dose regimen: 1.2mg then 0.6mg 1h later
## ============================================================
e5_flare <- ev(cmt="A_Crystal", amt=5, time=0)   # induce crystal flare
e5a <- ev(cmt="GUT_COLCH", amt=1.2, time=0)
e5b <- ev(cmt="GUT_COLCH", amt=0.6, time=1)
e5c <- ev(cmt="GUT_COLCH", amt=0.6, ii=12, addl=7, time=12) # maintenance 5 days

scen5 <- mod %>%
    ev(e5_flare + e5a + e5b + e5c) %>%
    mrgsim(end=14*24, delta=2) %>%
    as_tibble() %>%
    mutate(Scenario="5_Colchicine_acute")

## ============================================================
## Scenario 6: Acute Flare — Indomethacin 50mg TID
## ============================================================
e6_flare <- ev(cmt="A_Crystal", amt=5, time=0)
e6 <- ev(cmt="GUT_INDO", amt=50, ii=8, addl=20, time=0) # 50mg q8h × 7 days

scen6 <- mod %>%
    ev(e6_flare + e6) %>%
    mrgsim(end=14*24, delta=2) %>%
    as_tibble() %>%
    mutate(Scenario="6_Indomethacin_acute")

## ============================================================
## Scenario 7: Biologic — Canakinumab 150mg SC (refractory flares)
## ============================================================
e7_flare <- ev(cmt="A_Crystal", amt=5, time=0)
e7 <- ev(cmt="GUT_CANA", amt=150, time=0)

scen7 <- mod %>%
    ev(e7_flare + e7) %>%
    mrgsim(end=90*24, delta=12) %>%
    as_tibble() %>%
    mutate(Scenario="7_Canakinumab")

## ============================================================
## Scenario 8: Febuxostat 80mg + Flare prophylaxis (Colchicine 0.5mg/day)
## ============================================================
e8a <- ev(cmt="GUT_FEBU",  amt=80,  ii=24, addl=sim_end/24-1, time=0)
e8b <- ev(cmt="GUT_COLCH", amt=0.5, ii=24, addl=sim_end/24-1, time=0)

scen8 <- mod %>%
    ev(e8a + e8b) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="8_Febu_ColchProphylaxis")

## ============================================================
## RESULTS SUMMARY
## ============================================================
cat("\n=== Gout QSP Model (refactored) — 52-week Outcomes Summary ===\n")
summary_all <- list(scen1, scen2, scen3, scen4, scen8) %>%
    bind_rows() %>%
    filter(time == max(time)) %>%
    select(Scenario, sUA, XO_inh, URAT1_inh, Crystal, Tophus, eGFR_sim)

print(summary_all)

## ============================================================
## VISUALIZATION
## ============================================================
chronic_data <- bind_rows(scen1, scen2, scen3, scen4, scen8) %>%
    mutate(week = time / 168)

acute_data <- bind_rows(scen5, scen6, scen7) %>%
    mutate(day = time / 24)

p1 <- ggplot(chronic_data, aes(x=week, y=sUA, color=Scenario)) +
    geom_line(linewidth=1.0) +
    geom_hline(yintercept=6.0, linetype="dashed", color="red") +
    geom_hline(yintercept=5.0, linetype="dotted", color="blue") +
    labs(title="Serum Urate Over 52 Weeks",
         subtitle="Red dashed = target <6 mg/dL; Blue dotted = target <5 mg/dL (tophaceous)",
         x="Week", y="sUA (mg/dL)") +
    theme_bw() + theme(legend.position="bottom")

p2 <- ggplot(chronic_data, aes(x=week, y=Crystal, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="MSU Crystal Pool Over 52 Weeks",
         x="Week", y="Crystal Burden (normalized)") +
    theme_bw() + theme(legend.position="bottom")

p3 <- ggplot(acute_data, aes(x=day, y=Pain, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Acute Gout Pain Score (NRS 0-10)",
         x="Day", y="NRS Pain Score") +
    theme_bw() + theme(legend.position="bottom")

p4 <- ggplot(chronic_data, aes(x=week, y=IL1b_f, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Synovial IL-1b Over 52 Weeks",
         x="Week", y="IL-1b (pg/mL)") +
    theme_bw() + theme(legend.position="bottom")

p5 <- ggplot(chronic_data, aes(x=week, y=Tophus, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Tophus Volume Regression Over 52 Weeks",
         x="Week", y="Tophus Volume (cm3)") +
    theme_bw() + theme(legend.position="bottom")

p6 <- ggplot(chronic_data, aes(x=week, y=eGFR_sim, color=Scenario)) +
    geom_line(linewidth=1.0) +
    geom_hline(yintercept=60, linetype="dashed", color="orange") +
    labs(title="eGFR Trajectory Over 52 Weeks",
         subtitle="Orange = CKD stage 3 threshold",
         x="Week", y="eGFR (mL/min/1.73m2)") +
    theme_bw() + theme(legend.position="bottom")

## ============================================================
## SENSITIVITY ANALYSIS: sUA target achievement
## ============================================================
cat("\n=== sUA Target Achievement (<6 mg/dL) at Week 24 ===\n")
target_week24 <- bind_rows(scen2, scen3, scen4, scen8) %>%
    filter(abs(time - 24*7*24) < 24) %>%
    group_by(Scenario) %>%
    summarise(
        sUA_wk24    = mean(sUA),
        target_met  = mean(sUA) < 6.0,
        XO_inh_pct  = mean(XO_inh),
        .groups="drop"
    )
print(target_week24)

cat("\nModel compilation and simulation complete (refactored).\n")
