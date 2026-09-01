## =============================================================================
## Cholelithiasis (Gallstone Disease) — QSP Model
## mrgsolve ODE-based PK/PD Simulation
## REFACTORED for pluggable PK: see chol_refactor_notes.md. Original untouched at
## chol_mrgsolve_model.R; this file is a derived sibling, not an edit of it.
## =============================================================================
## Reference parameters calibrated from:
##   - Bachrach WH & Hofmann AF (1982) Ursodeoxycholic acid in the treatment
##     of cholesterol cholelithiasis. Dig Dis Sci 27:737-761
##   - Paumgartner G & Beuers U (2002) Ursodeoxycholic acid in cholestatic
##     liver disease. Hepatology 36:525-531
##   - Jazrawi RP et al (1992) Gut 33:381-386
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- mrgsolve model code block ----------------------------------------------
code <- '
$PROB Cholelithiasis QSP Model - UDCA/Statin PK/PD + Gallstone Dynamics (refactored, pluggable PK)

$PARAM
// ---- UDCA PK Parameters ----
DOSE_UDCA   = 750,   // total daily UDCA dose (mg/day), split TID
BWT         = 70,    // body weight (kg)
KA_UDCA     = 0.80,  // absorption rate constant UDCA (h^-1)
F_UDCA      = 0.50,  // oral bioavailability UDCA (fraction)
V1_UDCA     = 12.0,  // volume of distribution UDCA (L/kg) -> 840 L for 70 kg
CL_UDCA     = 40.0,  // total clearance UDCA (L/h)
KEHC_UDCA   = 0.35,  // enterohepatic cycling rate constant (h^-1)
FRAC_HEP_UDCA  = 0.70,  // hepatic first-pass extraction fraction
FRAC_BILE_UDCA = 0.80,  // fraction of hepatic UDCA secreted into bile

// ---- Statin PK Parameters (Simvastatin) ----
DOSE_STAT   = 0,     // statin daily dose (mg), 0 = off; 40 = typical
KA_STAT     = 0.60,  // absorption rate constant statin (h^-1)
F_STAT      = 0.05,  // oral bioavailability statin (first-pass ~95% extraction)
V1_STAT     = 3.2,   // volume of distribution statin (L/kg)
CL_STAT     = 70.0,  // total clearance statin (L/h)
EC50_STAT   = 0.008, // EC50 for HMGCR inhibition (mg/L)
EMAX_STAT   = 0.85,  // maximum HMGCR inhibition by statin (fraction)
GAMMA_STAT  = 1,     // Hill coefficient (math-implied by original’s plain ratio; new, not a fit)

// ---- Ezetimibe PK Parameters ----
DOSE_EZET   = 0,     // ezetimibe daily dose (mg), 0 = off; 10 = typical
KA_EZET     = 0.50,
F_EZET      = 0.35,
V1_EZET     = 4.0,
CL_EZET     = 8.0,
EC50_EZET   = 0.006, // EC50 for NPC1L1 inhibition (mg/L)
EMAX_EZET   = 0.80,  // maximum intestinal cholesterol absorption inhibition (never applied to any state in the original -- see refactor notes)
GAMMA_EZET  = 1,     // Hill coefficient (math-implied; new, not a fit)

// ---- Bile Acid Pool Dynamics ----
BA_synth0   = 0.52,  // baseline bile acid synthesis rate (g/day ~ 0.022 g/h)
BA_pool0    = 3.5,   // baseline total BA pool (g)
kBA_EHC     = 0.42,  // BA enterohepatic cycling rate (h^-1; ~6-10 cycles/day)
kBA_fecal   = 0.014, // BA fecal loss rate constant (h^-1)
K_BA_UDCA   = 0.25,  // UDCA effect coefficient on BA pool synthesis (was E_UDCA_BA)
KD_FXR      = 8.0,   // FXR activation concentration for BA feedback (umol/L)

// ---- Hepatic Cholesterol ----
CHOL_h0     = 15.0,  // baseline hepatic free cholesterol pool (mmol)
k_CHOL_syn  = 1.50,  // hepatic cholesterol synthesis rate (mmol/h)
k_CHOL_deg  = 0.10,  // hepatic cholesterol utilization/degradation (h^-1)
E_STAT_CHOL = 0.60,  // fractional reduction in cholesterol synthesis by statin -- unused in the original (dead parameter), preserved verbatim, not renamed
E_FXR_CHOL  = 0.20,  // FXR-mediated reduction in biliary cholesterol secretion (baseline disease term, not UDCA-dose-dependent)

// ---- Biliary Composition & Saturation ----
k_CHOL_bil  = 0.080, // rate of biliary cholesterol secretion (mmol/h)
k_PL_bil    = 0.18,  // rate of biliary phospholipid secretion (mmol/h)
PL_bil0     = 12.0,  // baseline biliary PL (mmol/L in GB bile)
BA_bil0     = 35.0,  // baseline biliary BA concentration (mmol/L in GB bile)
CHOL_bil0   = 4.2,   // baseline biliary cholesterol (mmol/L)
// CSI = CHOL_bil / (0.1875*BA_bil + 0.1429*PL_bil)

// ---- Gallbladder Volume & Motility ----
GB_vol0     = 30.0,  // fasting gallbladder volume (mL)
GB_vol_min  = 5.0,   // minimum GB volume (emptied, mL)
k_GB_fill   = 0.025, // GB filling rate constant (h^-1)
CCK_peak    = 1.0,   // normalized peak CCK (post-prandial)
k_GB_empty  = 0.30,  // GB emptying rate constant (h^-1)

// ---- Gallstone Dynamics ----
CSI_thresh    = 1.05,  // CSI threshold for nucleation
k_nucleat     = 0.0005,// crystal nucleation rate constant (mL/h per CSI unit above threshold)
k_grow        = 0.012, // stone growth rate constant (mL/h)
k_dissol      = 0.0,   // stone dissolution rate constant (mL/h) -- unused in the original (dead parameter; superseded by the computed k_dissol_eff), preserved verbatim, not renamed
K_DISSOL_UDCA = 0.025, // UDCA dissolution effect coefficient (/mmol/L UDCA in bile) (was E_UDCA_dis)
Stone_vol0    = 0.0,   // initial stone volume (mL), 0=prevention, >0=dissolution
Stone_max     = 5.0,   // maximum stone volume (mL, for Emax model of growth)

// ---- Inflammatory Markers ----
IL6_base    = 2.0,   // baseline IL-6 (pg/mL)
CRP_base    = 0.5,   // baseline CRP (mg/L)
k_IL6_prod  = 0.05,  // IL-6 production rate per stone volume unit
k_IL6_elim  = 0.15,  // IL-6 elimination rate (h^-1)
k_CRP_prod  = 0.30,  // CRP production stimulated by IL-6
k_CRP_elim  = 0.035, // CRP elimination rate (h^-1)

// ---- Simulation Flags ----
WLOSS       = 0,     // weight loss intervention (1=yes)
k_WL        = 0.002  // weight loss rate (fraction/h, ~1.5 kg/week)

$CMT
// PK compartments - UDCA (bespoke enterohepatic-cycling chain, see refactor notes)
GUT_UDCA   // gut (absorption) compartment for UDCA [mg]
CENT_UDCA  // plasma (systemic) UDCA [mg]
HEP_UDCA   // hepatic UDCA [mg]
BILE_UDCA  // biliary-duct UDCA [mg]
GB_UDCA    // gallbladder UDCA [mg]

// PK compartments - Statin (archetype 3 variant: depot+central, no peripheral)
GUT_STAT   // gut statin [mg]
CENT_STAT  // plasma statin [mg]

// PK compartments - Ezetimibe (archetype 3 variant: depot+central, no peripheral)
GUT_EZET   // gut ezetimibe [mg]
CENT_EZET  // plasma ezetimibe [mg]

// Bile acid & cholesterol dynamics (disease network, untouched)
BA_pool      // total bile acid pool [g]
CHOL_h       // hepatic free cholesterol [mmol]
CHOL_bil     // biliary cholesterol [mmol/L equivalent]
PL_bil       // biliary phospholipid [mmol/L equivalent]

// Gallbladder & stone (disease network, untouched)
GB_vol       // gallbladder volume [mL]
Crystal_mass // cholesterol crystal mass [mg]
Stone_V      // gallstone volume [mL]

// Inflammatory markers (disease network, untouched)
IL6          // IL-6 [pg/mL]
CRP_plas     // CRP [mg/L]

$MAIN
// Derived PK volumes
double V1_UDCA_L = V1_UDCA * BWT;
double V1_STAT_L = V1_STAT * BWT;
double V1_EZET_L = V1_EZET * BWT;

// ---- Statin: exposed concentration + Hill effect (rename of E_STAT, gamma=1) ----
double C_STAT = CENT_STAT / V1_STAT_L;
double EFFECT_STAT = (EMAX_STAT * pow(C_STAT, GAMMA_STAT)) / (pow(EC50_STAT, GAMMA_STAT) + pow(C_STAT, GAMMA_STAT));

// ---- Ezetimibe: exposed concentration + Hill effect (rename of E_EZET, gamma=1) ----
// NOTE: EFFECT_EZET is computed here exactly as the original computed E_EZET, but --
// exactly as in the original -- it is never read by any dxdt_ line or disease term
// anywhere in this model. Preserved verbatim (dangling), not wired to anything new.
// See chol_refactor_notes.md and UPSTREAM_ISSUES.md #131.
double C_EZET = CENT_EZET / V1_EZET_L;
double EFFECT_EZET = (EMAX_EZET * pow(C_EZET, GAMMA_EZET)) / (pow(EC50_EZET, GAMMA_EZET) + pow(C_EZET, GAMMA_EZET));

// ---- UDCA: two concentration sites (bespoke enterohepatic chain) ----
// C_UDCA = biliary UDCA (umol/L) -- this is the site that drives every UDCA PD effect below,
// exactly as the original’s own $MAIN-computed C_UDCA_bile did (this is a rename only).
// C_UDCA_PLAS = plasma UDCA (mg/L) -- informational only; not read by any PD term in the
// original (mirrors the abdominal-aortic-aneurysm doxycycline precedent: tissue site drives
// PD, plasma is kept only as a non-exposed diagnostic). $TABLE below also recomputes both of
// these quantities fresh, under different (legacy) names, for reporting -- see the $TABLE
// block and chol_refactor_notes.md for why that second compute-site is preserved rather than
// collapsed into this one (they are not numerically interchangeable).
double C_UDCA      = (BILE_UDCA / 392.6) * 1000.0;  // umol/L, MW UDCA = 392.6
double C_UDCA_PLAS = CENT_UDCA / V1_UDCA_L;         // mg/L

// FXR activation fraction (by BA pool via enterohepatic-delivered BAs)
double BA_conc_portal = (BA_pool / 3.5) * 40.0;  // umol/L proxy
double FXR_act = BA_conc_portal / (KD_FXR + BA_conc_portal);

// Weight during weight loss (if enabled)
// SOLVERTIME -> TIME: build-compat fix, mrgsolve 2.0.1 does not expose SOLVERTIME
// (_ODETIME_) inside $MAIN; see chol_refactor_notes.md and UPSTREAM_ISSUES.md #131.
double BWT_t = BWT * (1.0 - WLOSS * k_WL * TIME);
if(BWT_t < BWT * 0.75) BWT_t = BWT * 0.75;

// Cholesterol synthesis rate (inhibited by statin and SREBP2 feedback)
double k_CHOL_syn_eff = k_CHOL_syn * (1.0 - EFFECT_STAT) * (1.0 + 0.3*(1.0 - CHOL_h/CHOL_h0));

// Biliary cholesterol secretion rate (reduced by FXR, UDCA effect in bile)
// EFFECT_UDCA_DISSOL is UDCA’s own biliary effect term (renamed from the inline
// E_UDCA_dis * C_UDCA_bile_norm product, common-subexpression-factored since the
// original computed the identical product twice within this same $MAIN block, once
// here and again for k_dissol_eff below).
double C_UDCA_norm = C_UDCA / 500.0;  // normalize to typical biliary UDCA
double EFFECT_UDCA_DISSOL = K_DISSOL_UDCA * C_UDCA_norm;
double E_UDCA_CSI = E_FXR_CHOL + EFFECT_UDCA_DISSOL;
if(E_UDCA_CSI > 0.70) E_UDCA_CSI = 0.70;
double k_CHOL_bil_eff = k_CHOL_bil * (1.0 - E_UDCA_CSI) * (1.0 - EFFECT_STAT * 0.3);

// Compute CSI (Cholesterol Saturation Index)
// CSI = CHOL_bil / (0.1875 * BA_bil + 0.1429 * PL_bil) from Admirand-Small diagram
double BA_bil = BA_pool / 0.10;  // rough [BA]_bile from pool size
double CSI = CHOL_bil / (0.1875 * BA_bil + 0.1429 * PL_bil + 1e-6);

// Stone dissolution rate (UDCA-enhanced)
double k_dissol_eff = EFFECT_UDCA_DISSOL * 0.5;

// Crystal nucleation (occurs above CSI threshold)
double delta_CSI = (CSI > CSI_thresh) ? (CSI - CSI_thresh) : 0.0;
double nucleat_rate = k_nucleat * delta_CSI * GB_vol;

// Stone growth (sigmoidal inhibition by bile capacity)
double growth_rate = k_grow * delta_CSI * Stone_V * (1.0 - Stone_V / Stone_max);

// IL-6 production from stone-induced inflammation
double k_IL6_stim = k_IL6_prod * Stone_V * (Stone_V > 0.1 ? 1.0 : 0.0);

// Initial conditions
// _INIT(<CMT>) -> <CMT>_0: build-compat fix, mrgsolve 2.0.1 does not accept the
// _INIT() macro idiom; see chol_refactor_notes.md and UPSTREAM_ISSUES.md #131.
double Stone_V_IC = Stone_vol0;
if(NEWIND <= 1) {
    Stone_V_0      = Stone_vol0;
    Crystal_mass_0 = Stone_vol0 * 100.0;  // mg per mL stone
    BA_pool_0      = BA_pool0;
    CHOL_h_0       = CHOL_h0;
    CHOL_bil_0     = CHOL_bil0;
    PL_bil_0       = PL_bil0;
    GB_vol_0       = GB_vol0;
    IL6_0          = IL6_base;
    CRP_plas_0     = CRP_base;
}

$ODE
// ---- UDCA PK (bespoke enterohepatic chain) ----
double dose_UDCA_per = DOSE_UDCA / 3.0;  // TID dosing handled via event table

// Gut -> Plasma (with first-pass hepatic extraction)
dxdt_GUT_UDCA  = -KA_UDCA * GUT_UDCA;
double UDCA_absorbed = KA_UDCA * GUT_UDCA;
double UDCA_systemic = UDCA_absorbed * (1.0 - FRAC_HEP_UDCA);
double UDCA_liver_in = UDCA_absorbed * FRAC_HEP_UDCA;

// Plasma (simplified 1-cpt for plasma)
dxdt_CENT_UDCA = UDCA_systemic - (CL_UDCA / V1_UDCA_L) * CENT_UDCA;

// Hepatic UDCA
double UDCA_bile_out = FRAC_BILE_UDCA * KEHC_UDCA * HEP_UDCA;
dxdt_HEP_UDCA  = UDCA_liver_in - KEHC_UDCA * HEP_UDCA;

// Biliary UDCA (in bile ducts + GB)
double UDCA_gb_in = 0.40 * UDCA_bile_out;  // fraction going to GB during fasting
dxdt_BILE_UDCA = UDCA_bile_out - 0.50 * KEHC_UDCA * BILE_UDCA;

// GB UDCA (concentration during fasting)
dxdt_GB_UDCA   = UDCA_gb_in - KEHC_UDCA * GB_UDCA;

// ---- Statin PK ----
dxdt_GUT_STAT  = -KA_STAT * GUT_STAT;
dxdt_CENT_STAT = KA_STAT * GUT_STAT * F_STAT - (CL_STAT / V1_STAT_L) * CENT_STAT;

// ---- Ezetimibe PK ----
dxdt_GUT_EZET  = -KA_EZET * GUT_EZET;
dxdt_CENT_EZET = KA_EZET * GUT_EZET * F_EZET - (CL_EZET / V1_EZET_L) * CENT_EZET;

// ---- Bile Acid Pool ----
// BA synthesis (suppressed by FXR feedback from UDCA/BAs)
// EFFECT_UDCA_BA is UDCA’s own effect term on BA synthesis (renamed from the inline
// K_BA_UDCA(was E_UDCA_BA) * C_UDCA_norm product; kept here in $ODE, exactly where the
// original computed this product, per the guide’s "keep a calculation in the block the
// original used it in" rule -- not hoisted into $MAIN.
double EFFECT_UDCA_BA = K_BA_UDCA * C_UDCA_norm;
double BA_syn_rate = (BA_synth0 / 24.0) * (1.0 - 0.50 * FXR_act) * (1.0 + EFFECT_UDCA_BA);
// BA fecal loss
double BA_fecal = kBA_fecal * BA_pool;
// Net BA pool dynamics
dxdt_BA_pool = BA_syn_rate - BA_fecal;

// ---- Hepatic Cholesterol ----
// Synthesis inhibited by statin, LDLR-mediated uptake adds
double CHOL_uptake = 0.20 * (1.0 + EFFECT_STAT * 0.8);  // LDLR upregulated by statin
dxdt_CHOL_h = k_CHOL_syn_eff + CHOL_uptake - k_CHOL_deg * CHOL_h - k_CHOL_bil_eff;

// ---- Biliary Cholesterol & Phospholipids ----
dxdt_CHOL_bil = k_CHOL_bil_eff - 0.08 * CHOL_bil;  // secretion minus removal/dilution
dxdt_PL_bil   = k_PL_bil - 0.06 * PL_bil;           // PL secretion vs removal

// ---- Gallbladder Volume ----
// Fasting: fills slowly; post-prandial: CCK -> empties
dxdt_GB_vol = k_GB_fill * (GB_vol0 - GB_vol) - k_GB_empty * CCK_peak * GB_vol;
if(GB_vol < GB_vol_min) dxdt_GB_vol = 0.0;

// ---- Crystal & Stone Dynamics ----
// Crystal nucleation and growth
dxdt_Crystal_mass = nucleat_rate * 50.0 + growth_rate * 80.0 - k_dissol_eff * Crystal_mass;
if(Crystal_mass < 0.0) dxdt_Crystal_mass = 0.0;

// Stone volume
dxdt_Stone_V = growth_rate - k_dissol_eff * Stone_V;
if(Stone_V < 0.0) dxdt_Stone_V = 0.0;

// ---- Inflammatory Markers ----
dxdt_IL6     = IL6_base * k_IL6_elim + k_IL6_stim - k_IL6_elim * IL6;
dxdt_CRP_plas = k_CRP_prod * IL6 / (IL6_base + 1.0) * CRP_base - k_CRP_elim * CRP_plas;

$TABLE
double CSI_out = CHOL_bil / (0.1875 * (BA_pool / 0.10) + 0.1429 * PL_bil + 1e-6);
// UDCA_plas_conc / UDCA_bile_conc_umol / STAT_plas_conc / EZET_plas_conc: a fresh,
// $TABLE-cadence recompute of the same physical quantities C_UDCA_PLAS/C_UDCA/C_STAT/C_EZET
// represent in $MAIN. NOT a harmless textual duplicate -- confirmed empirically (see
// chol_refactor_notes.md) that $MAIN-cadence and $TABLE-cadence values genuinely diverge
// throughout the run (not just at dose instants), because $MAIN is evaluated once per
// record using state from the start of that record while $TABLE reads the fully-integrated
// state at the exact report time. The $MAIN-cadence values (C_UDCA, C_UDCA_PLAS, C_STAT,
// C_EZET below) are what the original’s own PD equations actually consumed and are kept as
// the exposed, pluggable concentration per this fork’s convention; these $TABLE-cadence
// values are preserved verbatim alongside them as the original’s own fresh diagnostic
// report columns (byte-identical formulas/names to the original, only compartment names
// renamed) -- not collapsed into the $MAIN versions, since doing so would silently change
// what gets reported.
double UDCA_plas_conc = CENT_UDCA / (V1_UDCA * BWT);        // mg/L
double UDCA_bile_conc_umol = (BILE_UDCA / 392.6) * 1000.0;  // umol/L
double STAT_plas_conc = CENT_STAT / (V1_STAT * BWT);        // mg/L
double EZET_plas_conc = CENT_EZET / (V1_EZET * BWT);        // mg/L -- never captured in the original at all; now captured (see notes)
double Stone_mm = pow(Stone_V * 6.0 / 3.14159, 1.0/3.0) * 10.0; // approx diameter mm (sphere)
double BA_pool_g = BA_pool;
double CHOL_sat_pct = CSI_out * 100.0;

$CAPTURE CSI_out UDCA_plas_conc UDCA_bile_conc_umol STAT_plas_conc EZET_plas_conc
$CAPTURE Stone_mm BA_pool_g CHOL_sat_pct
$CAPTURE C_UDCA C_UDCA_PLAS C_STAT C_EZET EFFECT_STAT EFFECT_EZET EFFECT_UDCA_DISSOL EFFECT_UDCA_BA
'

## ---- Compile model ----------------------------------------------------------
mod <- mcode("cholelithiasis_qsp_refactored", code)

## ---- Helper: create dosing events -------------------------------------------
make_events <- function(dose_UDCA = 750, dose_STAT = 0, dose_EZET = 0,
                        dur_days = 365, freq_UDCA = "TID") {
  per_UDCA <- if(freq_UDCA == "TID") dose_UDCA / 3 else dose_UDCA
  int_UDCA <- if(freq_UDCA == "TID") 8 else 24   # hours between doses

  ev_list <- list()

  if(dose_UDCA > 0)
    ev_list$udca  <- ev(amt = per_UDCA, cmt = "GUT_UDCA",
                        ii = int_UDCA, addl = dur_days * (24 / int_UDCA) - 1)
  if(dose_STAT > 0)
    ev_list$stat  <- ev(amt = dose_STAT,   cmt = "GUT_STAT",
                        ii = 24, addl = dur_days - 1)
  if(dose_EZET > 0)
    ev_list$ezet  <- ev(amt = dose_EZET,   cmt = "GUT_EZET",
                        ii = 24, addl = dur_days - 1)

  if(length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = 1))
  Reduce(ev_seq, ev_list)
}

## ---- Treatment Scenarios -----------------------------------------------------
scenarios <- list(
  "Scenario 1: No Treatment (Natural History)" = list(
    DOSE_UDCA = 0, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 2: UDCA 750 mg/day (Standard)" = list(
    DOSE_UDCA = 750, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 3: UDCA 1050 mg/day (High Dose)" = list(
    DOSE_UDCA = 1050, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 4: UDCA + Simvastatin 40 mg" = list(
    DOSE_UDCA = 750, DOSE_STAT = 40, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 5: Ezetimibe 10 mg Prevention" = list(
    DOSE_UDCA = 0, DOSE_STAT = 0, DOSE_EZET = 10,
    Stone_vol0 = 0.0, WLOSS = 0
  ),
  "Scenario 6: Lifestyle (Weight Loss) + UDCA" = list(
    DOSE_UDCA = 750, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 1
  )
)

## ---- Run simulations ---------------------------------------------------------
sim_duration  <- 365     # days
sim_delta     <- 1       # hourly resolution? No, daily output: every 24h
sim_end_h     <- sim_duration * 24
sim_times     <- seq(0, sim_end_h, by = 24)  # daily output

run_scenario <- function(sc_name, params) {
  ev_data <- make_events(
    dose_UDCA = params$DOSE_UDCA,
    dose_STAT = params$DOSE_STAT,
    dose_EZET = params$DOSE_EZET,
    dur_days  = sim_duration
  )

  out <- mod %>%
    param(DOSE_UDCA = params$DOSE_UDCA,
          DOSE_STAT  = params$DOSE_STAT,
          DOSE_EZET  = params$DOSE_EZET,
          Stone_vol0 = params$Stone_vol0,
          WLOSS      = params$WLOSS) %>%
    mrgsim(events = ev_data,
           end    = sim_end_h,
           delta  = 24,
           obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(Scenario = sc_name,
           Day = time / 24)

  return(out)
}

results <- bind_rows(
  mapply(run_scenario, names(scenarios), scenarios, SIMPLIFY = FALSE)
)

## ---- Plot 1: Stone Volume Dissolution / Progression -------------------------
p1 <- results %>%
  ggplot(aes(x = Day, y = Stone_V, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.1) +
  labs(title    = "Gallstone Volume Over Time",
       subtitle = "Scenarios 1–4: Stone Volume (mL) — Dissolution vs. Growth",
       x = "Time (days)", y = "Stone Volume (mL)",
       color = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 3))

## ---- Plot 2: Cholesterol Saturation Index -----------------------------------
p2 <- results %>%
  ggplot(aes(x = Day, y = CHOL_sat_pct, color = Scenario)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 20, y = 102, label = "CSI = 1.0 (Saturation threshold)",
           color = "red", size = 3.5) +
  geom_line(linewidth = 1.0) +
  labs(title    = "Biliary Cholesterol Saturation Index (%)",
       subtitle = "CSI > 100% = lithogenic bile",
       x = "Time (days)", y = "CSI (%)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 3: Bile Acid Pool Dynamics ----------------------------------------
p3 <- results %>%
  ggplot(aes(x = Day, y = BA_pool_g, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Total Bile Acid Pool (g)",
       x = "Time (days)", y = "BA Pool (g)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 4: UDCA Biliary Concentration (Scenarios with UDCA) ---------------
p4 <- results %>%
  filter(grepl("UDCA", Scenario)) %>%
  ggplot(aes(x = Day, y = UDCA_bile_conc_umol, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title    = "UDCA Biliary Concentration",
       subtitle = "Steady-state biliary UDCA (µmol/L)",
       x = "Time (days)", y = "UDCA in Bile (µmol/L)") +
  theme_classic(base_size = 13)

## ---- Plot 5: Inflammatory Markers -------------------------------------------
p5 <- results %>%
  ggplot(aes(x = Day, y = CRP_plas, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Plasma CRP Over Time",
       x = "Time (days)", y = "CRP (mg/L)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 6: Hepatic Cholesterol Over Time ----------------------------------
p6 <- results %>%
  ggplot(aes(x = Day, y = CHOL_h, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Hepatic Cholesterol Pool (mmol)",
       x = "Time (days)", y = "Hepatic Cholesterol (mmol)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Summary Table: Key Endpoints at 6 months and 12 months ----------------
endpoint_summary <- results %>%
  filter(Day %in% c(0, 180, 365)) %>%
  group_by(Scenario, Day) %>%
  summarise(
    Stone_Vol_mL    = round(mean(Stone_V), 3),
    Stone_Diam_mm   = round(mean(Stone_mm), 1),
    CSI_pct         = round(mean(CHOL_sat_pct), 1),
    BA_pool_g       = round(mean(BA_pool_g), 2),
    UDCA_bile_uM    = round(mean(UDCA_bile_conc_umol), 0),
    CRP_mgL         = round(mean(CRP_plas), 2),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Day)

print(endpoint_summary)

## ---- Sensitivity Analysis: UDCA dose vs dissolution rate at Day 365 --------
udca_doses <- seq(250, 1500, by = 250)
dose_response <- lapply(udca_doses, function(d) {
  ev_d <- make_events(dose_UDCA = d, dur_days = 365)
  out  <- mod %>%
    param(DOSE_UDCA = d, Stone_vol0 = 0.5) %>%
    mrgsim(events = ev_d, end = 365*24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    filter(time == 365*24) %>%
    mutate(Dose_UDCA = d,
           Pct_dissolv = (1 - Stone_V / 0.5) * 100)
  out
}) %>% bind_rows()

p_dose_resp <- dose_response %>%
  ggplot(aes(x = Dose_UDCA, y = Pct_dissolv)) +
  geom_point(size = 3, color = "#1565C0") +
  geom_line(color = "#1565C0", linewidth = 1.1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(title    = "UDCA Dose-Response: Stone Dissolution at 12 Months",
       subtitle = "% stone dissolved from initial volume (0.5 mL)",
       x = "UDCA Daily Dose (mg/day)", y = "Stone Dissolution (%)") +
  theme_classic(base_size = 13)

cat("\n=== Cholelithiasis QSP Model — Simulation Complete ===\n")
cat("Scenarios simulated:", length(scenarios), "\n")
cat("Simulation duration:", sim_duration, "days\n")
cat("\nPlots generated: p1 (stone volume), p2 (CSI), p3 (BA pool),\n")
cat("  p4 (UDCA bile conc), p5 (CRP), p6 (hepatic CHOL), p_dose_resp\n")
cat("\nEndpoint summary table printed above.\n")
