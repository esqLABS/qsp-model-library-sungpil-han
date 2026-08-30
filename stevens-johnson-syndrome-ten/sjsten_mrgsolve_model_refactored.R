## =========================================================================
## Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)
## Quantitative Systems Pharmacology (QSP) Model — mrgsolve
## ** PK/PD-refactored sibling ** — see sjsten_refactor_notes.md
##
## Disease scope:
##   Drug-induced severe cutaneous adverse reaction (SCAR).
##   Captures: culprit drug exposure → drug-specific CD8⁺ CTL expansion
##   → cytotoxic effectors (granulysin, sFasL, perforin/GzmB, TNF-α, IFN-γ)
##   → keratinocyte apoptosis/necroptosis → epidermal detachment (BSA%)
##   → SCORTEN-driven mortality hazard
##   → re-epithelialization kinetics under treatment.
##
## Treatment scenarios:
##   1) Supportive care only (drug withdrawal + ICU)
##   2) IVIG 1 g/kg/d × 4
##   3) Cyclosporine 3 mg/kg PO BID × 10 d (Valeyrie-Allanore 2010)
##   4) Etanercept 25 mg SC × 2 (Wang 2018 JCI RCT)
##   5) Methylprednisolone 1 mg/kg IV pulse
##   6) Combination: Cyclosporine + Etanercept
##   7) JAK inhibitor (tofacitinib 5 mg BID) — investigational
##
## Calibration references:
##   - SCORTEN: Bastuji-Garin 2000 J Invest Dermatol (mortality vs score)
##   - Granulysin kinetics: Chung 2008 Nat Med; Fujita 2011 J Allergy Clin Immunol
##   - Etanercept: Wang 2018 J Clin Invest (mortality 8.3% vs 47.6%, p<0.01)
##   - Cyclosporine: Valeyrie-Allanore 2010 Br J Dermatol (mortality SMR 0.43)
##   - IVIG: Viard 1998 Science (mechanism); Roujeau 2017 meta (mortality NS)
##
## Units:
##   Time = days; concentrations = mg/L (drug) or pg/mL (cytokines);
##   BSA detachment = % body surface area.
## =========================================================================

[PROB]
SJS/TEN QSP model — drug → CD8⁺ CTL → granulysin/FasL/TNF → keratinocyte
apoptosis → BSA detachment → SCORTEN → mortality + re-epithelialization.
Refactored for pluggable PK: every compound exposes a single named
C_<STEM> concentration and, where the original's own effect term is
Hill/Emax-shaped, a single named EFFECT_<STEM>. See sjsten_refactor_notes.md.

[PLUGIN] Rcpp

[PARAM]
// ===== Patient covariates =====
WT        = 65,      // kg
AGE       = 45,      // years
HLA_RISK  = 1,       // 0 = absent, 1 = present (e.g. HLA-B*15:02 with CBZ)
PREDOSE   = 1,       // continued culprit drug? 0/1

// ===== Culprit drug PK (illustrative: carbamazepine-like) [stem: DRUG] =====
CL_DRUG   = 4.0,     // L/h
V1_DRUG   = 70,      // L
KA_DRUG   = 0.5,     // 1/h

// ===== Antigen presentation =====
KON_HLA   = 0.05,    // drug x HLA -> antigen formation
KOFF_HLA  = 0.10,
HLA_FOLD  = 50,      // OR-like risk multiplier if HLA_RISK=1

// ===== T-cell activation =====
EMAX_TACT = 5.0,
EC50_TACT = 2.0,     // antigen units
KDEG_T    = 0.20,    // 1/d
TACT_BASE = 0.05,    // basal activation

// ===== Cytokine kinetics =====
KSYN_TNF  = 8.0,     // pg/mL/d per activated CTL
KDEG_TNF  = 6.0,     // 1/d
KSYN_IFN  = 12.0,
KDEG_IFN  = 4.0,
KSYN_IL15 = 4.0,
KDEG_IL15 = 1.5,
KSYN_GNLY = 50.0,    // ng/mL/d
KDEG_GNLY = 0.5,
KSYN_FASL = 4.0,     // pg/mL/d
KDEG_FASL = 1.0,
KSYN_HMGB = 20.0,
KDEG_HMGB = 2.0,

// ===== Keratinocyte apoptosis kinetics =====
EMAX_APOP = 1.5,     // max KC death rate (1/d)
EC50_GNLY = 25.0,    // ng/mL
EC50_TNF  = 50.0,    // pg/mL
EC50_FASL = 8.0,     // pg/mL
KC0       = 100,     // arbitrary epidermal mass (100 = full)
KREEPI    = 0.07,    // 1/d basal re-epithelialization
KREEPI_EGF = 0.05,
EGF_STIM  = 1.0,     // baseline EGF/KGF drive

// ===== BSA detachment dynamics =====
BSA_MAX   = 100,     // %
KBSA_LOSS = 0.10,
KBSA_REEPI = 0.05,

// ===== SCORTEN-based mortality =====
BASELINE_HAZ = 0.005, // per day
SCORTEN_BETA = 0.55,  // ln-OR per point
AGE_THRESH = 40,
HR_THRESH  = 120,
BUN_VAL    = 9.0,    // mmol/L baseline
GLC_VAL    = 8.0,    // mmol/L
HCO3_VAL   = 22,
CA_VAL     = 1.0,    // 1 = present malignancy

// ===== Treatment scenarios (toggle 0/1, dose mg) — Hill params renamed =====
SCEN_IVIG  = 0,
DOSE_IVIG  = 1.0,    // g/kg/d
DAYS_IVIG  = 4,
EMAX_IVIG  = 0.5,    // was EFF_IVIG_FASL: 50% sFasL neutralization
EC50_IVIG  = 0.5,    // was inline literal 0.5 in the original's ratio
GAMMA_IVIG = 1.0,    // original had no explicit Hill exponent (implicit 1)

SCEN_CSA   = 0,
DOSE_CSA   = 3.0,    // mg/kg/d divided BID
EMAX_CSA   = 0.6,    // was EFF_CSA_T: 60% T activation inhibition
EC50_CSA   = 0.05,   // was inline literal 0.05
GAMMA_CSA  = 1.0,

SCEN_ETAN  = 0,
DOSE_ETAN  = 25,     // mg SC
EMAX_ETAN  = 0.8,    // was EFF_ETAN_TNF
EC50_ETAN  = 0.3,    // was inline literal 0.3
GAMMA_ETAN = 1.0,

SCEN_INFL  = 0,
DOSE_INFL  = 5,      // mg/kg
EMAX_INFL  = 0.85,   // was EFF_INFL_TNF
EC50_INFL  = 0.5,    // was inline literal 0.5
GAMMA_INFL = 1.0,

SCEN_PRED  = 0,
DOSE_PRED  = 1.0,    // mg/kg
EMAX_PRED  = 0.5,    // was EFF_PRED_T
EC50_PRED  = 0.05,   // was inline literal 0.05
GAMMA_PRED = 1.0,

SCEN_TPE   = 0,      // therapeutic plasma exchange — not a modeled compound, untouched
EFF_TPE_GNLY = 0.7,

SCEN_JAKI  = 0,
DOSE_JAKI  = 5,      // mg BID
EMAX_JAKI  = 0.7,    // was EFF_JAKI_IFN
EC50_JAKI  = 0.1,    // was inline literal 0.1
GAMMA_JAKI = 1.0,

SCEN_WD    = 1,      // drug withdrawal toggle (1 = stopped on day 0)
DAY_WD     = 0,

// ===== Drug PK for biologics/immunomodulators [renamed to V1_<STEM> convention] =====
CL_IVIG    = 0.05,   // L/d
V1_IVIG    = 5,
KA_ETAN    = 0.4,    // 1/d SC
CL_ETAN    = 0.07,
V1_ETAN    = 7.6,
CL_INFL    = 0.27,
V1_INFL    = 4,
CL_CSA     = 5,      // L/h
V1_CSA     = 100,
KA_CSA     = 1.4,
CL_PRED    = 12,     // L/h
V1_PRED    = 50,
KA_JAKI    = 6.0,
CL_JAKI    = 22,
V1_JAKI    = 96

[CMT]
CENT_DRUG   // 1  culprit drug central (amount) [stem: DRUG]
GUT_DRUG    // 2  culprit drug depot
Ag_HLA      // 3  drug-HLA complex (antigen units) -- disease state, unchanged
T_act       // 4  activated CD8 CTL clones
TNF         // 5  pg/mL
IFNg        // 6  pg/mL
IL15        // 7  pg/mL
GNLY        // 8  granulysin ng/mL
sFasL       // 9  pg/mL
HMGB1       // 10 ng/mL
KC_alive    // 11 epidermal mass
BSA_loss    // 12 % BSA detachment
Cum_hazard  // 13 cumulative mortality hazard
Surv        // 14 survival probability
CENT_IVIG   // 15 mg/L [stem: IVIG] -- bespoke: state IS the concentration (see notes)
GUT_ETAN    // 16 [stem: ETAN]
CENT_ETAN   // 17 mg/L -- bespoke
CENT_INFL   // 18 mg/L [stem: INFL] -- bespoke, no depot (IV)
GUT_CSA     // 19 [stem: CSA]
CENT_CSA    // 20 mg/L -- bespoke
CENT_PRED   // 21 mg/L [stem: PRED] -- bespoke, no depot (IV)
GUT_JAKI    // 22 [stem: JAKI]
CENT_JAKI   // 23 mg/L -- bespoke

[MAIN]
// initial conditions
CENT_DRUG_0 = 0;
KC_alive_0 = KC0;
BSA_loss_0 = 0;
Surv_0 = 1;

[ODE]
// ---------- 0. Culprit drug PK (if continued) [stem: DRUG, Archetype 3 minus peripheral] ----------
// NOTE: drug_input is computed but never used anywhere else in the original
// file (PREDOSE/SCEN_WD/DAY_WD have no downstream effect) -- preserved
// verbatim as dead code, not fixed; see "Anything else flagged" in
// sjsten_refactor_notes.md and UPSTREAM_ISSUES.md.
double drug_input = (PREDOSE > 0.5 && (SCEN_WD < 0.5 || SOLVERTIME < DAY_WD)) ? 1.0 : 0.0;
dxdt_GUT_DRUG  = -KA_DRUG*GUT_DRUG;
dxdt_CENT_DRUG =  KA_DRUG*GUT_DRUG - (CL_DRUG/V1_DRUG)*CENT_DRUG;

// the exposed concentration
double C_DRUG = CENT_DRUG / V1_DRUG; // mg/L

// ---------- 1. Antigen presentation (drug + HLA) ----------
double hla_amp = (HLA_RISK > 0.5) ? HLA_FOLD : 1.0;
// bespoke: the culprit drug's own driving term is linear mass-action
// (KON_HLA * hla_amp * C_DRUG), not a saturating Hill/Emax shape, so no
// EMAX_DRUG/EC50_DRUG/GAMMA_DRUG were invented for it (see notes). Named
// per the guide's "one named function of concentration" principle anyway.
double EFFECT_DRUG = hla_amp * C_DRUG;
dxdt_Ag_HLA = KON_HLA*EFFECT_DRUG - KOFF_HLA*Ag_HLA;

// ---------- 2. Treatment factors ----------
// bespoke: each of these six compartments already holds a concentration
// directly (no CENT/V division in the original's own effect terms, even
// though CL/V is used correctly as the elimination rate constant) -- see
// "C_<STEM> is an identity of the compartment" in sjsten_refactor_notes.md.
double C_CSA  = CENT_CSA;
double C_PRED = CENT_PRED;
double C_ETAN = CENT_ETAN;
double C_INFL = CENT_INFL;
double C_IVIG = CENT_IVIG;
double C_JAKI = CENT_JAKI;

double EFFECT_CSA  = EMAX_CSA  * pow(C_CSA,  GAMMA_CSA)  / (pow(EC50_CSA,  GAMMA_CSA)  + pow(C_CSA,  GAMMA_CSA));
double EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED, GAMMA_PRED) + pow(C_PRED, GAMMA_PRED));
double EFFECT_ETAN = EMAX_ETAN * pow(C_ETAN, GAMMA_ETAN) / (pow(EC50_ETAN, GAMMA_ETAN) + pow(C_ETAN, GAMMA_ETAN));
double EFFECT_INFL = EMAX_INFL * pow(C_INFL, GAMMA_INFL) / (pow(EC50_INFL, GAMMA_INFL) + pow(C_INFL, GAMMA_INFL));
double EFFECT_IVIG = EMAX_IVIG * pow(C_IVIG, GAMMA_IVIG) / (pow(EC50_IVIG, GAMMA_IVIG) + pow(C_IVIG, GAMMA_IVIG));
double EFFECT_JAKI = EMAX_JAKI * pow(C_JAKI, GAMMA_JAKI) / (pow(EC50_JAKI, GAMMA_JAKI) + pow(C_JAKI, GAMMA_JAKI));

double f_csa  = 1.0 - (SCEN_CSA  > 0.5 ? EFFECT_CSA  : 0.0);
double f_pred = 1.0 - (SCEN_PRED > 0.5 ? EFFECT_PRED : 0.0);
double f_tnf_etan = 1.0 - (SCEN_ETAN > 0.5 ? EFFECT_ETAN : 0.0);
double f_tnf_infl = 1.0 - (SCEN_INFL > 0.5 ? EFFECT_INFL : 0.0);
double f_tnf      = f_tnf_etan * f_tnf_infl;
double f_fasl     = 1.0 - (SCEN_IVIG > 0.5 ? EFFECT_IVIG : 0.0);
double f_ifn      = 1.0 - (SCEN_JAKI > 0.5 ? EFFECT_JAKI : 0.0);
double f_gnly     = 1.0 - (SCEN_TPE  > 0.5 ? EFF_TPE_GNLY * exp(-(SOLVERTIME-1)/2.0) : 0.0); // TPE: not a compound, untouched

// ---------- 3. T-cell activation ----------
double drive = EMAX_TACT * Ag_HLA / (Ag_HLA + EC50_TACT);
dxdt_T_act   = drive*f_csa*f_pred + TACT_BASE - KDEG_T*T_act;

// ---------- 4. Cytokines (driven by T_act, IL-15 amplifies CTL) ----------
dxdt_TNF   = KSYN_TNF *T_act*f_tnf - KDEG_TNF *TNF;
dxdt_IFNg  = KSYN_IFN *T_act*f_ifn - KDEG_IFN *IFNg;
dxdt_IL15  = KSYN_IL15*T_act       - KDEG_IL15*IL15;
dxdt_GNLY  = KSYN_GNLY*T_act*f_gnly - KDEG_GNLY*GNLY;
dxdt_sFasL = KSYN_FASL*T_act       - KDEG_FASL*sFasL*(1.0 - (1.0 - f_fasl));
dxdt_HMGB1 = KSYN_HMGB*(KC0 - KC_alive)/KC0 - KDEG_HMGB*HMGB1;

// ---------- 5. Keratinocyte apoptosis ----------
double e_gnly = GNLY  / (GNLY  + EC50_GNLY);
double e_tnf  = TNF   / (TNF   + EC50_TNF);
double e_fasl = sFasL / (sFasL + EC50_FASL);
double k_apop = EMAX_APOP * (1.0 - (1-e_gnly)*(1-e_tnf)*(1-e_fasl));

dxdt_KC_alive = -k_apop*KC_alive + KREEPI*(KC0 - KC_alive)*EGF_STIM;

// ---------- 6. BSA detachment ----------
double damage_frac = (KC0 - KC_alive) / KC0;
dxdt_BSA_loss = KBSA_LOSS*damage_frac*(BSA_MAX - BSA_loss)
                - KBSA_REEPI*BSA_loss*(KC_alive/KC0);

// ---------- 7. SCORTEN-based mortality hazard ----------
double scorten = 0;
if (AGE > AGE_THRESH) scorten += 1;
if (BSA_loss > 10)    scorten += 1;
scorten += (BUN_VAL > 10)  ? 1 : 0;
scorten += (HCO3_VAL< 20)  ? 1 : 0;
scorten += (GLC_VAL > 14)  ? 1 : 0;
scorten += (CA_VAL  > 0.5) ? 1 : 0;
scorten += (HR_THRESH>120) ? 1 : 0;

double rx_haz_mult = 1.0;
if (SCEN_ETAN > 0.5) rx_haz_mult *= 0.18; // Wang 2018 reduction
if (SCEN_CSA  > 0.5) rx_haz_mult *= 0.43; // Valeyrie-Allanore
if (SCEN_IVIG > 0.5) rx_haz_mult *= 0.80;

double hazard = BASELINE_HAZ * exp(SCORTEN_BETA*scorten) * rx_haz_mult;
dxdt_Cum_hazard = hazard;
dxdt_Surv       = -hazard*Surv;

// ---------- 8. Biologic / immunomodulator PK ----------
double r_ivig = (SCEN_IVIG > 0.5 && SOLVERTIME < DAYS_IVIG) ? DOSE_IVIG*WT/V1_IVIG : 0.0; // edge case: zero-order infusion built into the ODE, not an ev()
dxdt_CENT_IVIG = r_ivig - (CL_IVIG/V1_IVIG)*CENT_IVIG;

dxdt_GUT_ETAN  = -KA_ETAN*GUT_ETAN;
dxdt_CENT_ETAN =  KA_ETAN*GUT_ETAN - (CL_ETAN/V1_ETAN)*CENT_ETAN;

dxdt_CENT_INFL = - (CL_INFL/V1_INFL)*CENT_INFL;

dxdt_GUT_CSA  = -KA_CSA*GUT_CSA;
dxdt_CENT_CSA =  KA_CSA*GUT_CSA - (CL_CSA/V1_CSA)*CENT_CSA;

dxdt_CENT_PRED = - (CL_PRED/V1_PRED)*CENT_PRED;

dxdt_GUT_JAKI  = -KA_JAKI*GUT_JAKI;
dxdt_CENT_JAKI =  KA_JAKI*GUT_JAKI - (CL_JAKI/V1_JAKI)*CENT_JAKI;

[TABLE]
double SCORTEN = (AGE>AGE_THRESH) + (BSA_loss>10)
               + (BUN_VAL>10) + (HCO3_VAL<20) + (GLC_VAL>14)
               + (CA_VAL>0.5) + (HR_THRESH>120);
double PredMort = 1 - Surv;
double Re_epi = 100*(KC_alive/KC0);

[CAPTURE]
C_DRUG SCORTEN PredMort Re_epi EFFECT_DRUG
C_CSA C_PRED C_ETAN C_INFL C_IVIG C_JAKI
EFFECT_CSA EFFECT_PRED EFFECT_ETAN EFFECT_INFL EFFECT_IVIG EFFECT_JAKI

/*
=========================================================================
Example R driver — load model, run scenarios, plot
(renamed to match the refactored compartment/parameter names above)
=========================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

mod <- mread("sjsten_mrgsolve_model_refactored.R")

scenarios <- list(
  Supportive  = list(),
  IVIG        = list(SCEN_IVIG=1),
  Cyclosporine= list(SCEN_CSA=1),
  Etanercept  = list(SCEN_ETAN=1),
  Methylpred  = list(SCEN_PRED=1),
  CSA_ETAN    = list(SCEN_CSA=1, SCEN_ETAN=1),
  JAK_inhib   = list(SCEN_JAKI=1)
)

ev_drug <- ev(amt=400, ii=12, addl=2, cmt="GUT_DRUG")    # culprit
ev_ivig <- ev(amt=0, time=0, cmt="CENT_IVIG")            # via rate term
ev_etan <- ev(amt=25, ii=3, addl=1, cmt="GUT_ETAN")
ev_csa  <- ev(amt=100, ii=0.5, addl=20, cmt="GUT_CSA")   # ~3 mg/kg BID
ev_pred <- ev(amt=60, ii=1,  addl=5, cmt="CENT_PRED")
ev_jaki <- ev(amt=5,  ii=0.5, addl=20, cmt="GUT_JAKI")

simulate <- function(par) {
  e <- ev_drug
  if (isTRUE(par$SCEN_ETAN==1)) e <- c(e, ev_etan)
  if (isTRUE(par$SCEN_CSA ==1)) e <- c(e, ev_csa)
  if (isTRUE(par$SCEN_PRED==1)) e <- c(e, ev_pred)
  if (isTRUE(par$SCEN_JAKI==1)) e <- c(e, ev_jaki)
  mod %>% param(par) %>% ev(e) %>% mrgsim(end=30, delta=0.1) %>% as_tibble()
}

out <- bind_rows(lapply(names(scenarios), function(n)
  simulate(scenarios[[n]]) %>% mutate(scenario=n)))

ggplot(out, aes(time, BSA_loss, color=scenario)) + geom_line() +
  labs(x="Day", y="BSA detachment (%)", title="SJS/TEN — treatment effect")

ggplot(out, aes(time, PredMort, color=scenario)) + geom_line() +
  labs(x="Day", y="Predicted mortality")
*/
