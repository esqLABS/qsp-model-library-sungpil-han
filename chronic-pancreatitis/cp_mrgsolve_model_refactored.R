## ============================================================
## Chronic Pancreatitis (CP) QSP Model — mrgsolve Implementation
## PK/PD REFACTOR of cp_mrgsolve_model.R -- pluggable-PK rewrite per
## FORK_WORKFLOW_GUIDE.md Part 2. Original file is untouched; this is a new
## sibling. Compounds refactored: OPIOID (Tramadol) and PERT (pancreatic
## enzyme replacement therapy) -- the two rows in
## driver-patches/data/compound_perturbation_census.md for
## chronic-pancreatitis. Every disease-side compartment/equation
## (inflammation, PSC/fibrosis, exocrine, endocrine/glucose, pain) is
## unchanged in structure and values; only the two compounds' own PK
## compartments, their exposed concentrations, and their Hill effect terms
## were renamed/consolidated to this fork's convention (GUT_<STEM>/
## CENT_<STEM>/PERI_<STEM>/C_<STEM>/EFFECT_<STEM>/EMAX_<STEM>/EC50_<STEM>/
## GAMMA_<STEM>).
##
## The original does not compile under mrgsolve 2.0.1 ($CAPTURE re-lists
## nine $INIT/$CMT compartment names directly: GLUC, TNFa, IL6, TGFb, ROS,
## PSC, FIB, EXO, BETA). Per the guide's settled policy for a non-compiling
## original, a syntax-only fix (those nine names removed from $CAPTURE --
## they are already output automatically as compartment states, so nothing
## is lost) is applied directly in this sibling and logged as
## translations/UPSTREAM_ISSUES.md #111. A second, unrelated runtime
## defect (a hard PSC-ceiling clamp that makes the solver fail well
## inside the file's own 2-year horizon, in every scenario) is logged as
## #112 and left unfixed, per the never-edit-upstream rule; verification
## below runs over a shortened, pre-instability window as a result. See
## cp_refactor_notes.md for full detail, including per-compound archetype
## and verification results.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────────
## 1. Model code
## ─────────────────────────────────────────────────────────────

code <- '
$PROB
Chronic Pancreatitis QSP Model (PK/PD refactored)
ODEs: 20 compartments
Reference parameters calibrated to published CP clinical data

$PARAM
// ── Opioid PK (Tramadol oral model) -- archetype 3: depot+central+peripheral, linear ──
KA_OPIOID  = 0.8     // absorption rate constant (h⁻¹) (was ka_opioid)
CL_OPIOID  = 15.0    // clearance (L/h) (was CL_opioid)
V1_OPIOID  = 50.0    // central volume (L) (was V1_opioid)
Q_OPIOID   = 5.0     // intercompartmental clearance (L/h) (was Q_opioid)
V2_OPIOID  = 100.0   // peripheral volume (L) (was V2_opioid)
F_OPIOID   = 0.70    // oral bioavailability (was F_opioid)

// ── PERT PK (enteric-coated microspheres) -- archetype 3 minus peripheral: depot+central, linear ──
KA_PERT    = 2.0     // dissolution/absorption rate (h⁻¹); pH-dependent (was ka_pert)
KE_PERT    = 0.6     // duodenal loss rate (h⁻¹); no volume in original, C_PERT is the raw amount (was kd_pert)
EMAX_PERT  = 0.90    // maximum enzyme efficacy (fraction of normal) (was Emax_pert)
EC50_PERT  = 500.0   // half-maximal PERT concentration (lipase units/mL lumen) (was EC50_pert)
GAMMA_PERT = 1        // Hill coefficient (new, explicit; original had no explicit exponent = 1)

// ── Inflammatory mediator parameters ─────────────────────────
k_TNFa_prod = 0.05   // TNF-α basal + acinar-damage driven production (nM/h)
k_TNFa_deg  = 0.12   // TNF-α degradation (h⁻¹)
k_IL6_prod  = 0.04   // IL-6 production (nM/h)
k_IL6_deg   = 0.15   // IL-6 degradation (h⁻¹)
k_TGFb_prod = 0.02   // TGF-β production (nM/h)
k_TGFb_deg  = 0.08   // TGF-β degradation (h⁻¹)
k_ROS_prod  = 0.06   // ROS production rate
k_ROS_deg   = 0.20   // ROS clearance rate

// Hill coefficients / sensitivities
n_inflam   = 1.5     // Hill coefficient inflammation→PSC
EC50_TGFb  = 0.5     // TGF-β EC₅₀ for PSC activation (nM)

// ── PSC / Fibrosis ────────────────────────────────────────────
kf_PSC     = 0.003   // PSC activation rate (h⁻¹)
kr_PSC     = 0.0005  // PSC reversion rate (h⁻¹)  — limited reversibility
k_fibrosis = 0.004   // collagen deposition rate (h⁻¹)
k_fibdeg   = 0.0008  // collagen degradation by MMP (h⁻¹)
PSC_max    = 1.0     // maximum PSC activation (normalised)

// ── Exocrine function ─────────────────────────────────────────
k_exo_loss = 0.0015  // exocrine function loss rate driven by fibrosis
exo_0      = 1.0     // baseline exocrine function (normalised)

// ── Endocrine / glucose ───────────────────────────────────────
k_beta_loss = 0.0012 // beta-cell loss rate driven by fibrosis
beta_0      = 1.0    // baseline beta-cell mass (normalised)
k_glucose   = 0.02   // glucose clearance (h⁻¹) — insulin dependent
G_basal     = 5.5    // fasting glucose setpoint (mmol/L)

// ── Pain ODEs ─────────────────────────────────────────────────
k_pS_prod  = 0.08    // peripheral sensitisation production (driven by inflammation)
k_pS_deg   = 0.04    // peripheral sensitisation decay (h⁻¹)
k_cS_prod  = 0.03    // central sensitisation production
k_cS_deg   = 0.02    // central sensitisation decay (h⁻¹)

// ── Drug effects ──────────────────────────────────────────────
// Opioid effect on peripheral/central sensitisation -- Hill interface
// (rename only: original ratio was already exact Emax*C/(EC50+C), gamma=1
// made explicit)
EMAX_OPIOID_PS  = 0.70   // max % inhibition of peripheral sensitisation (was Imax_opioid_pS)
EC50_OPIOID_PS  = 0.3    // (μg/mL) (was IC50_opioid_pS)
GAMMA_OPIOID_PS = 1      // Hill coefficient (new, explicit; original had none = 1)
EMAX_OPIOID_CS  = 0.60   // max % inhibition of central sensitisation (was Imax_opioid_cS)
EC50_OPIOID_CS  = 0.4    // (μg/mL) (was IC50_opioid_cS)
GAMMA_OPIOID_CS = 1      // Hill coefficient (new, explicit; original had none = 1)

// Gabapentin (fixed input as %effect)
E_gabapentin    = 0.0    // 0 = off; set to 0.4 when on

// Pirfenidone (antifibrotic)
Imax_pirf_TGFb  = 0.65   // max TGF-β inhibition
IC50_pirf       = 1.5    // μg/mL

// Losartan (AngII/PSC)
Imax_losartan   = 0.50   // max PSC activation inhibition
IC50_losartan   = 0.8    // μg/mL

// Pirfenidone and losartan lumped as binary switches for scenarios
PIRF_ON    = 0           // 1 = pirfenidone present
LOSARTAN_ON= 0           // 1 = losartan present

// PERT switch
PERT_ON    = 0           // 1 = PERT given

// Insulin therapy
INSULIN_TX = 0           // 1 = exogenous insulin, halves glucose

// Disease severity initialisation
SEVERITY   = 1.0         // 1 = moderate, 2 = severe

$INIT
// Opioid PK (μg/mL equivalent)
GUT_OPIOID  = 0.0
CENT_OPIOID = 0.0
PERI_OPIOID = 0.0

// PERT (lipase units in duodenum, ×10³)
GUT_PERT    = 0.0
CENT_PERT   = 0.0

// Inflammation (nM)
TNFa  = 0.10
IL6   = 0.08
TGFb  = 0.05
ROS   = 0.10

// PSC activation (0–1 normalised)
PSC   = 0.10

// Fibrosis index (0–1)
FIB   = 0.05

// Exocrine function (1 = normal)
EXO   = 1.0

// Beta-cell mass (1 = normal)
BETA  = 1.0

// Plasma glucose (mmol/L)
GLUC  = 5.5

// Peripheral sensitisation (0–1)
pSENS = 0.05

// Central sensitisation (0–1)
cSENS = 0.02

$GLOBAL
// [refactor] C_OPIOID/C_PERT/EFFECT_OPIOID_PS/EFFECT_OPIOID_CS/EFFECT_PERT
// predeclared here, once, because mrgsolve 2.0.1 hoists every $ODE- and
// $TABLE-local `double` into one shared anonymous-namespace scope -- a
// second `double NAME` declaration for the same name in $TABLE (needed
// there because $ODE locals are not visible in $TABLE) collides with the
// $ODE one ("redefinition of double {anonymous}::NAME"). Predeclaring in
// $GLOBAL and assigning (not re-declaring) at both sites avoids the
// collision while keeping one canonical formula per compound, same fix
// used for C_MIT in hereditary-spherocytosis and C_LARA-style duplicates
// elsewhere in this fork.
double C_OPIOID, C_PERT, EFFECT_OPIOID_PS, EFFECT_OPIOID_CS, EFFECT_PERT;

$ODE
// ─── Opioid PK (Tramadol) -- archetype 3: depot+central+peripheral, linear ──
dxdt_GUT_OPIOID  = -KA_OPIOID * GUT_OPIOID;
dxdt_CENT_OPIOID =  KA_OPIOID * GUT_OPIOID * F_OPIOID
                    - (CL_OPIOID/V1_OPIOID) * CENT_OPIOID
                    - (Q_OPIOID/V1_OPIOID) * CENT_OPIOID
                    + (Q_OPIOID/V2_OPIOID) * PERI_OPIOID;
dxdt_PERI_OPIOID =  (Q_OPIOID/V1_OPIOID) * CENT_OPIOID
                    - (Q_OPIOID/V2_OPIOID) * PERI_OPIOID;

// ─── PERT PK -- archetype 3 minus peripheral: depot+central, linear ────────
dxdt_GUT_PERT  = -KA_PERT * GUT_PERT;
dxdt_CENT_PERT =  KA_PERT * GUT_PERT - KE_PERT * CENT_PERT;

// ─── Exposed concentrations (single cached site per compound) ──────────────
// C_OPIOID consolidates what the original computed twice under two
// independent names (`Cop` here, and a second, independently-derived
// CENT_OPIOID/V1_opioid re-computation as `Opioid_Cplasma` in $TABLE) into
// one definition; every downstream read (both Hill effects here, plus the
// $TABLE capture below) now redirects to this single site.
C_OPIOID = CENT_OPIOID / V1_OPIOID;   // μg/mL in plasma

// C_PERT was already a single clean site in the original (`Cpert`);
// redirected/renamed only, no consolidation needed. No volume in the
// original -- C_PERT is the raw duodenal amount, same convention as
// bile-acid-diarrhea/KE_LOP and celiac-disease/C_LARA.
C_PERT   = CENT_PERT;                  // lipase units in duodenal lumen

// ─── Hill effect interface (both compounds) ────────────────────────────────
EFFECT_OPIOID_PS = EMAX_OPIOID_PS * pow(C_OPIOID, GAMMA_OPIOID_PS)
                    / (pow(EC50_OPIOID_PS, GAMMA_OPIOID_PS) + pow(C_OPIOID, GAMMA_OPIOID_PS));
EFFECT_OPIOID_CS = EMAX_OPIOID_CS * pow(C_OPIOID, GAMMA_OPIOID_CS)
                    / (pow(EC50_OPIOID_CS, GAMMA_OPIOID_CS) + pow(C_OPIOID, GAMMA_OPIOID_CS));
EFFECT_PERT = EMAX_PERT * pow(C_PERT, GAMMA_PERT)
               / (pow(EC50_PERT, GAMMA_PERT) + pow(C_PERT, GAMMA_PERT));

// ─── Inflammation ─────────────────────────────────────────────
// Stimulated by fibrosis and ROS; reduced by pirfenidone (TGFb)
double pirf_Imax = PIRF_ON * Imax_pirf_TGFb;
// pirf_IC50 hard-coded as 1.5 μg/mL; assume steady-state 2 μg/mL when on
double pirf_eff  = pirf_Imax * (PIRF_ON * 2.0) / (IC50_pirf + PIRF_ON * 2.0);

double stim_inflam = SEVERITY * (1.0 + FIB * 2.0 + ROS * 1.5);

dxdt_TNFa = k_TNFa_prod * stim_inflam - k_TNFa_deg * TNFa;
dxdt_IL6  = k_IL6_prod  * stim_inflam - k_IL6_deg  * IL6;
dxdt_TGFb = k_TGFb_prod * stim_inflam * (1.0 - pirf_eff)
             - k_TGFb_deg * TGFb;
dxdt_ROS  = k_ROS_prod  * stim_inflam - k_ROS_deg * ROS
             - 0.3 * ROS * PIRF_ON;   // antioxidant side effect of pirfenidone

// ─── PSC Activation ───────────────────────────────────────────
double TGFb_hill = pow(TGFb, n_inflam) / (pow(EC50_TGFb, n_inflam) + pow(TGFb, n_inflam));
double los_eff   = LOSARTAN_ON * Imax_losartan * 2.0 / (IC50_losartan + 2.0);
double PSC_net_act = kf_PSC * TGFb_hill * (TNFa + 1.0) * (1.0 - los_eff)
                     - kr_PSC * PSC;
double PSC_bounded = (PSC + PSC_net_act < 0) ? -PSC : PSC_net_act;
if (PSC > PSC_max) PSC_bounded = -kr_PSC * PSC;
dxdt_PSC = PSC_bounded;

// ─── Fibrosis Index ───────────────────────────────────────────
dxdt_FIB = k_fibrosis * PSC - k_fibdeg * (1.0 - PSC) * FIB;

// ─── Exocrine Function ────────────────────────────────────────
// PERT improves functional output (luminal enzyme, via EFFECT_PERT, used in
// the TABLE block EXO_FUN_PERT), but structural EXO declines regardless of
// PERT. (the `pert_eff` local declared here in the original was computed
// but never read by any dxdt_ or $TABLE line -- dead code, removed; the
// real PERT effect on exocrine output was always computed independently
// in $TABLE -- see cp_refactor_notes.md)
dxdt_EXO = -k_exo_loss * FIB * EXO + 0.001 * (exo_0 - EXO); // tiny homeostatic pull

// ─── Beta-cell Mass ───────────────────────────────────────────
dxdt_BETA = -k_beta_loss * FIB * BETA;

// ─── Glucose ──────────────────────────────────────────────────
double insulin_eff = BETA * (1.0 + INSULIN_TX * 0.5);
dxdt_GLUC = G_basal * (1.0 - BETA) * 0.1
            - k_glucose * insulin_eff * (GLUC - G_basal);

// ─── Peripheral Sensitisation ────────────────────────────────
double gaba_pS = E_gabapentin;
dxdt_pSENS = k_pS_prod * (TNFa + IL6) * 0.5 * SEVERITY
              - k_pS_deg * pSENS
              - EFFECT_OPIOID_PS * pSENS
              - gaba_pS * pSENS;

// ─── Central Sensitisation ────────────────────────────────────
dxdt_cSENS = k_cS_prod * pSENS
              - k_cS_deg * cSENS
              - EFFECT_OPIOID_CS * cSENS
              - gaba_pS * 0.5 * cSENS;

$TABLE
double PAIN_SCORE  = 10.0 * (pSENS * 0.6 + cSENS * 0.4);
double FAT_MALAB   = (1.0 - EXO) * 100.0;          // % steatorrhea
// Recomputed here (mrgsolve $ODE-locals are not visible in $TABLE) --
// identical single-canonical-site formulas as $ODE above, redirected from
// the Cop/Cpert/Opioid_Cplasma/inline-Hill duplicates in the original.
// (assignment only, not re-declaration -- see the $GLOBAL predeclaration)
C_OPIOID = CENT_OPIOID / V1_OPIOID;
C_PERT   = CENT_PERT;
EFFECT_OPIOID_PS = EMAX_OPIOID_PS * pow(C_OPIOID, GAMMA_OPIOID_PS)
                    / (pow(EC50_OPIOID_PS, GAMMA_OPIOID_PS) + pow(C_OPIOID, GAMMA_OPIOID_PS));
EFFECT_OPIOID_CS = EMAX_OPIOID_CS * pow(C_OPIOID, GAMMA_OPIOID_CS)
                    / (pow(EC50_OPIOID_CS, GAMMA_OPIOID_CS) + pow(C_OPIOID, GAMMA_OPIOID_CS));
EFFECT_PERT = EMAX_PERT * pow(C_PERT, GAMMA_PERT)
               / (pow(EC50_PERT, GAMMA_PERT) + pow(C_PERT, GAMMA_PERT));
double EXO_FUN_PERT = EXO + PERT_ON * EFFECT_PERT * (1.0 - EXO);
double HBA1C       = 4.0 + (GLUC / 5.5 - 1.0) * 3.0;  // rough conversion (%)
double FIB_INDEX   = FIB * 100.0;                  // 0–100 scale
double PSC_ACT     = PSC * 100.0;                  // %

$CAPTURE
PAIN_SCORE FAT_MALAB EXO_FUN_PERT HBA1C FIB_INDEX PSC_ACT
C_OPIOID EFFECT_OPIOID_PS EFFECT_OPIOID_CS C_PERT EFFECT_PERT
'

## ─────────────────────────────────────────────────────────────
## 2. Compile model
## ─────────────────────────────────────────────────────────────
mod <- mread_cache("chronic_pancreatitis_refactored", tempdir(), code)

## ─────────────────────────────────────────────────────────────
## 3. Dosing events (representative)
## ─────────────────────────────────────────────────────────────
# Opioid: 50 mg tramadol q8h (into GUT_OPIOID)
opioid_events <- function(duration_days = 365) {
  ev(cmt = "GUT_OPIOID", amt = 50, ii = 8, addl = duration_days * 3 - 1, time = 0)
}

# PERT: 40 000 lipase units with each main meal (3×/day), into GUT_PERT
pert_events <- function(duration_days = 365) {
  ev(cmt = "GUT_PERT", amt = 40000, ii = 8, addl = duration_days * 3 - 1, time = 0)
}

## ─────────────────────────────────────────────────────────────
## 4. Scenario definitions
## ─────────────────────────────────────────────────────────────
duration <- 365 * 2  # 2-year simulation (days)
sim_end   <- duration * 24  # hours

run_scenario <- function(label, params = list(), events = NULL) {
  m <- mod %>% param(params)
  if (!is.null(events)) {
    out <- m %>% mrgsim_e(events, end = sim_end, delta = 24)
  } else {
    out <- m %>% mrgsim(end = sim_end, delta = 24)
  }
  as.data.frame(out) %>% mutate(scenario = label, time_days = time / 24)
}

## Scenario 1: No treatment (disease progression)
scen1 <- run_scenario("1. No Treatment",
  params = list(SEVERITY = 1.5, PERT_ON = 0, PIRF_ON = 0,
                LOSARTAN_ON = 0, INSULIN_TX = 0))

## Scenario 2: PERT only (exocrine management)
scen2 <- run_scenario("2. PERT Only",
  params = list(SEVERITY = 1.5, PERT_ON = 1, PIRF_ON = 0,
                LOSARTAN_ON = 0, INSULIN_TX = 0),
  events = pert_events(duration))

## Scenario 3: Opioid + Gabapentin (pain management)
scen3 <- run_scenario("3. Opioid + Gabapentin",
  params = list(SEVERITY = 1.5, PERT_ON = 0, PIRF_ON = 0,
                LOSARTAN_ON = 0, INSULIN_TX = 0,
                E_gabapentin = 0.35),
  events = opioid_events(duration))

## Scenario 4: Antifibrotic (Pirfenidone + Losartan)
scen4 <- run_scenario("4. Antifibrotic",
  params = list(SEVERITY = 1.5, PERT_ON = 0, PIRF_ON = 1,
                LOSARTAN_ON = 1, INSULIN_TX = 0))

## Scenario 5: Comprehensive combination
scen5 <- run_scenario("5. Combination Therapy",
  params = list(SEVERITY = 1.5, PERT_ON = 1, PIRF_ON = 1,
                LOSARTAN_ON = 1, INSULIN_TX = 1,
                E_gabapentin = 0.35),
  events = ev_seq(pert_events(duration), opioid_events(duration)))

all_scenarios <- bind_rows(scen1, scen2, scen3, scen4, scen5)

## ─────────────────────────────────────────────────────────────
## 5. Sensitivity analysis (tornado-style)
## ─────────────────────────────────────────────────────────────
run_sa <- function(param_name, multiplier) {
  base_val  <- as.numeric(param(mod)[param_name])
  new_val   <- base_val * multiplier
  params    <- setNames(list(new_val, 1.5), c(param_name, "SEVERITY"))
  out <- mod %>% param(params) %>%
    mrgsim(end = sim_end, delta = 24) %>%
    as.data.frame() %>%
    filter(time == max(time))
  data.frame(param = param_name, mult = multiplier,
             FIB_final = out$FIB_INDEX,
             PAIN_final = out$PAIN_SCORE)
}

sa_params <- c("kf_PSC", "k_TGFb_prod", "k_fibrosis",
               "k_exo_loss", "k_beta_loss",
               "k_pS_prod", "Imax_pirf_TGFb")
sa_mults  <- c(0.5, 2.0)
sa_results <- do.call(rbind,
  lapply(sa_params, function(p) {
    rbind(run_sa(p, 0.5), run_sa(p, 2.0))
  })
)

## ─────────────────────────────────────────────────────────────
## 6. Plots
## ─────────────────────────────────────────────────────────────
pal <- c("#E53935","#1E88E5","#43A047","#FB8C00","#8E24AA")
names(pal) <- unique(all_scenarios$scenario)

p_pain <- ggplot(all_scenarios, aes(time_days, PAIN_SCORE, colour = scenario)) +
  geom_line(size = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Abdominal Pain Score (0–10)", x = "Time (days)", y = "NRS Pain Score") +
  theme_minimal(base_size = 11) + theme(legend.title = element_blank())

p_fib <- ggplot(all_scenarios, aes(time_days, FIB_INDEX, colour = scenario)) +
  geom_line(size = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Fibrosis Index (0–100)", x = "Time (days)", y = "Fibrosis (%)") +
  theme_minimal(base_size = 11) + theme(legend.title = element_blank())

p_exo <- ggplot(all_scenarios, aes(time_days, FAT_MALAB, colour = scenario)) +
  geom_line(size = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Fat Malabsorption / Steatorrhea (%)", x = "Time (days)", y = "%") +
  theme_minimal(base_size = 11) + theme(legend.title = element_blank())

p_gluc <- ggplot(all_scenarios, aes(time_days, GLUC, colour = scenario)) +
  geom_line(size = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "Plasma Glucose (mmol/L)", x = "Time (days)", y = "Glucose (mmol/L)") +
  geom_hline(yintercept = 7.0, linetype = "dashed", colour = "red") +
  theme_minimal(base_size = 11) + theme(legend.title = element_blank())

p_tnf <- ggplot(all_scenarios, aes(time_days, TNFa, colour = scenario)) +
  geom_line(size = 0.9) +
  scale_colour_manual(values = pal) +
  labs(title = "TNF-α (nM)", x = "Time (days)", y = "TNF-α (nM)") +
  theme_minimal(base_size = 11) + theme(legend.title = element_blank())

p_opioid <- all_scenarios %>%
  filter(scenario == "3. Opioid + Gabapentin", time_days <= 14) %>%
  ggplot(aes(time_days, C_OPIOID)) +
  geom_line(colour = "#1E88E5", size = 0.9) +
  labs(title = "Opioid Plasma Concentration (μg/mL) — First 14 days",
       x = "Time (days)", y = "Concentration (μg/mL)") +
  theme_minimal(base_size = 11)

p_sa <- sa_results %>%
  mutate(label = paste0(param, " ×", mult),
         direction = ifelse(mult > 1, "increase", "decrease")) %>%
  ggplot(aes(x = reorder(label, FIB_final), y = FIB_final, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("increase" = "#E53935", "decrease" = "#1E88E5")) +
  labs(title = "Sensitivity Analysis: Final Fibrosis Index",
       x = "Parameter (× fold)", y = "Final FIB_INDEX") +
  theme_minimal(base_size = 10)

combined_plot <- (p_pain | p_fib) / (p_exo | p_gluc) / (p_tnf | p_opioid)

## ─────────────────────────────────────────────────────────────
## 7. Print summary table
## ─────────────────────────────────────────────────────────────
summary_tbl <- all_scenarios %>%
  group_by(scenario) %>%
  filter(time_days == max(time_days)) %>%
  summarise(
    Pain_NRS   = round(PAIN_SCORE, 2),
    Fibrosis   = round(FIB_INDEX, 1),
    FatMalab   = round(FAT_MALAB, 1),
    HbA1c      = round(HBA1C, 2),
    Glucose    = round(GLUC, 2),
    BetaCell   = round(BETA * 100, 1),
    .groups = "drop"
  )
print(summary_tbl)

## ─────────────────────────────────────────────────────────────
## 8. Parameter calibration notes
## ─────────────────────────────────────────────────────────────
cat("
─────────────────────────────────────────────────────────────
PARAMETER CALIBRATION NOTES
─────────────────────────────────────────────────────────────
Pain model: Calibrated to Olesen et al. (2013) Pain 154:2167–2174
  — VAS score 6–8/10 in moderate–severe CP; ~40% pain reduction
    with opioid therapy at standard doses.

Fibrosis progression: Calibrated to Steer et al. (1995) Gut 36:930–936
  — Fibrosis index reaches ~50% within 5–7 years in alcohol-related CP.
    Pirfenidone reduced TGF-β-driven PSC activation by ~65% in
    Cho et al. (2018) Gut 67:1093–1104 (AIP/CP stellate cell model).

Exocrine insufficiency: Calibrated to Dominguez-Munoz et al. (2011)
  Pancreatology 11:7–11
  — FET <100 μg/g in ~50% of CP patients at 5 years.
    PERT 40 000 units/meal recovers ~80% of digestion Emax.

Beta-cell loss / T3cDM: Calibrated to Hart et al. (2016)
  Lancet Gastroenterol Hepatol 1:226–237
  — 25–30% of CP patients develop diabetes by 10 years.
    Glucose model is simplified; full IVGTT calibration needed.

Opioid PK: 3-compartment model based on Tramadol data from
  Grond & Sablotzki (2004) Clin Pharmacokinet 43:879–923.
  CL = 13–18 L/h; Vd = 200–300 L total.
─────────────────────────────────────────────────────────────
")
