## ============================================================
## Essential Thrombocythemia (ET) — QSP Model (PK/PD-refactored)
## mrgsolve ODE-based PK/PD Simulation
##
## This is a fork-owned sibling of `et_mrgsolve_model.R` (never edit the
## original — see FORK_WORKFLOW_GUIDE.md Part 2). It re-derives the same
## disease biology unchanged, and rewrites the PK/PD of the four in-scope
## compounds (Anagrelide/ANA, Hydroxyurea/HU, Pegylated interferon/PIFN,
## Ruxolitinib/RUX — all four rows for this file in
## `driver-patches/data/compound_perturbation_census.md`) into the guide's
## named "pluggable PK" convention: GUT_<STEM>/CENT_<STEM>/PERI_<STEM>,
## C_<STEM> as the single exposed concentration, EFFECT_<STEM> as the named
## Hill interface. See `et_refactor_notes.md` for the full rationale,
## build-compat fixes applied, and verification result.
##
## Compartments (21):
##   HSC  – Hematopoietic stem cell pool (JAK2-mutant fraction)
##   MKP  – Megakaryocyte progenitors
##   MK   – Mature megakaryocytes
##   PLT  – Circulating platelets (×10⁹/L)
##   TPO  – Serum thrombopoietin (pg/mL)
##   JAK2 – JAK2 allele burden (fraction 0-1)
##   SPL  – Spleen size (cm below costal margin)
##   GUT_HU / CENT_HU / PERI_HU       – Hydroxyurea depot/central/peripheral (µg/mL-equivalent)
##   GUT_ANA / CENT_ANA / PERI_ANA    – Anagrelide depot/central/peripheral (ng/mL-equivalent)
##   GUT_RUX / CENT_RUX / PERI_RUX    – Ruxolitinib depot/central/peripheral (ng/mL-equivalent)
##   GUT_PIFN / CENT_PIFN / PERI_PIFN – Peg-IFN-α2a depot/central/peripheral (µg/mL-equivalent)
##   RISK_T  – Cumulative thrombosis risk (AU)
##   RISK_MF – Cumulative MF transformation risk (AU)
##
## Key References:
##   Harrison 2005 NEJM (HU vs ANA PT1 trial)
##   Gisslinger 2013 Blood (ANAHYDRET trial)
##   Verstovsek 2020 Leukemia (RUX in ET)
##   Kiladjian 2013 Haematologica (Peg-IFN in ET)
##   Barosi 2009 Leukemia (ELN response criteria)
##   Rumi 2014 Blood (CALR vs JAK2 prognosis)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ── Model code ──────────────────────────────────────────────
code <- '
$PROB
Essential Thrombocythemia QSP Model (PK/PD-refactored)
21-compartment ODE | 7 treatment scenarios

$PARAM
// ── Disease biology parameters ──
k_HSC_prod  = 0.01    // HSC basal production rate (1/day)
k_HSC_diff  = 0.008   // HSC differentiation rate (1/day)
k_MKP_diff  = 0.15    // MKP differentiation to MK (1/day)
k_MK_mature = 0.10    // MK maturation rate (1/day)
k_PLT_prod  = 8.0     // platelet production per MK unit (×10⁹/L/day)
k_PLT_destr = 0.115   // platelet destruction rate (1/day, t½≈8d)
k_TPO_prod  = 0.50    // TPO basal production (pg/mL/day)
k_TPO_elim  = 0.012   // TPO elimination (1/day)
k_TPO_abs   = 0.00008 // TPO absorption by platelets (per PLT/day)
k_JAK2_exp  = 0.003   // JAK2 clone expansion rate (1/day)
k_SPL_grow  = 0.008   // spleen growth rate (cm/day)
k_SPL_shrk  = 0.05    // spleen regression rate (1/day)
phi_JAK2    = 3.5     // JAK2 gain-of-function amplification on HSC prod

// ── Steady-state reference values ──
PLT_ss = 900           // baseline ET platelet (×10⁹/L)
PLT_norm = 250         // normal platelet count
MK_ss  = 1.0           // reference MK pool
TPO_ss = 85            // normal TPO (pg/mL)
HSC_ss = 1.0

// ── PK: Hydroxyurea (renamed ka_HU → KA_HU; CL/V1/Q/V2 already matched convention) ──
F_HU   = 0.80          // bioavailability
KA_HU  = 1.40          // absorption (1/h)
CL_HU  = 4.20          // clearance (L/h)
V1_HU  = 28.0          // central Vd (L)
Q_HU   = 1.80          // intercompartmental CL (L/h)
V2_HU  = 12.0          // peripheral Vd (L)
MW_HU  = 76.05         // mol. weight (g/mol) → µg/mL to µM conv (unused by ODE, kept verbatim)

// ── PK: Anagrelide (renamed ka_ANA → KA_ANA) ──
F_ANA  = 0.70
KA_ANA = 4.60          // 1/h
CL_ANA = 9.50          // L/h
V1_ANA = 18.0
Q_ANA  = 2.10
V2_ANA = 8.00

// ── PK: Ruxolitinib (renamed ka_RUX → KA_RUX) ──
F_RUX  = 0.95
KA_RUX = 2.30          // 1/h
CL_RUX = 22.0
V1_RUX = 72.0
Q_RUX  = 8.50
V2_RUX = 38.0

// ── PK: Peg-IFN-α2a (SC depot, weekly; stem renamed pIFN → PIFN throughout) ──
F_PIFN  = 0.84
KA_PIFN = 0.020        // 1/h (slow SC absorption)
CL_PIFN = 0.038        // L/h (hepatic clearance)
V1_PIFN = 8.00
Q_PIFN  = 0.30
V2_PIFN = 4.50

// ── Dosing amount parameters ──
// Build-compat fix (see et_refactor_notes.md, logged as UPSTREAM_ISSUES.md
// #61): the original file never declared DOSE_HU/DOSE_ANA/DOSE_RUX/
// DOSE_pIFN/DOSE_ASA in its own $PARAM block even though its $ODE reads
// them directly (only ever set at the R side via param()/mrgsim(param=)),
// so it fails to compile at all under mrgsolve 2.0.1
// ("’DOSE_HU’ was not declared in this scope", etc). Declared here with
// the same safe defaults the original’s own R wrapper already used
// (`mod2 <- mod %>% param(DOSE_HU = 0, ...)`). Purely a missing-declaration
// fix, no value invented beyond the original’s own default of 0.
DOSE_HU   = 0
DOSE_ANA  = 0
DOSE_RUX  = 0
DOSE_PIFN = 0
DOSE_ASA  = 0

// ── PD: Effect parameters (renamed Emax_* → EMAX_*, gam_* → GAMMA_*) ──
// Hydroxyurea (on MKP)
EC50_HU   = 3.50       // µg/mL for 50% effect on MKP
EMAX_HU   = 0.85       // max fractional MKP suppression
GAMMA_HU  = 1.20       // Hill coefficient

// Anagrelide (on MK maturation)
EC50_ANA  = 0.025      // µg/mL = 25 ng/mL (Ki PDE3 ~36 nM → ~16 ng/mL ANA)
EMAX_ANA  = 0.75
GAMMA_ANA = 1.50

// Ruxolitinib (on JAK2 signaling → MKP proliferation)
EC50_RUX  = 0.15       // µg/mL = 150 ng/mL (IC50 JAK2 ~3.3 nM × 404 g/mol ≈ 1.3 ng/mL, shifted by protein binding)
EMAX_RUX  = 0.80
GAMMA_RUX = 1.80

// Peg-IFN-α2a (on JAK2 clone suppression)
EC50_PIFN  = 0.008     // µg/mL = 8 ng/mL
EMAX_PIFN  = 0.70      // max JAK2 clone suppression rate
GAMMA_PIFN = 1.00

// ── Thrombosis & MF risk parameters ─────────────────────
lambda_T   = 0.0005    // baseline thrombosis hazard/day (≈1.8%/yr low-risk)
alpha_T    = 1.80      // platelet exponent for thrombosis risk
delta_JAK2 = 0.80      // JAK2 allele burden contribution to thrombosis
lambda_MF  = 0.00012   // baseline MF transformation hazard/day (≈0.4%/yr)
delta_ASXL = 1.50      // ASXL1 co-mutation MF multiplier (set via param)

// ── Co-variate flags ────────────────────────────────────
ASXL1_pos  = 0         // 1 if ASXL1 positive
TP53_pos   = 0
age_gt60   = 0         // 1 if age >60
prior_thrombo = 0      // 1 if prior thrombosis

$INIT
// $CMT block removed: this mrgsolve build (2.0.1) rejects a bare $CMT list
// plus a same-named $INIT block as "Duplicated model names" (same defect
// class as UPSTREAM_ISSUES.md #34/#36/#42/.../#51). $INIT alone declares
// every compartment, same fix as issue #51.
HSC    = 1.0
MKP    = 1.2
MK     = 1.5
PLT    = 900.0
TPO    = 60.0
JAK2   = 0.55           // 55% JAK2 V617F allele burden at baseline
SPL    = 2.5            // 2.5 cm splenomegaly at baseline
GUT_HU    = 0.0
CENT_HU   = 0.0
PERI_HU   = 0.0
GUT_ANA   = 0.0
CENT_ANA  = 0.0
PERI_ANA  = 0.0
GUT_RUX   = 0.0
CENT_RUX  = 0.0
PERI_RUX  = 0.0
GUT_PIFN  = 0.0
CENT_PIFN = 0.0
PERI_PIFN = 0.0
RISK_T  = 0.0
RISK_MF = 0.0

$MAIN
// Single administration at t=0, same mechanism the original used (a
// constant DOSE_<STEM> parameter driving a decaying absorption profile
// from time 0), now expressed as a real depot compartment initial amount
// instead of an analytic exp() term baked into the central-compartment
// ODE (see et_refactor_notes.md for the algebraic equivalence proof).
GUT_HU_0   = F_HU   * DOSE_HU;
GUT_ANA_0  = F_ANA  * DOSE_ANA;
GUT_RUX_0  = F_RUX  * DOSE_RUX;
GUT_PIFN_0 = F_PIFN * DOSE_PIFN;

$ODE
// ── Drug PK: Hydroxyurea — Archetype 3 (depot + central + peripheral) ──
// NOTE: the original’s own central/peripheral exchange is NOT the textbook
// mass-balance form the guide’s archetype-3 template shows (Q/V1 outflow,
// Q/V2 return). It instead applies the SAME Q/V1 coefficient to the
// (CENT-PERI) difference in both directions, plus an independent Q/V2
// loss term applied only to the peripheral compartment. Preserved exactly
// as written (see et_refactor_notes.md) — this is a rename/reorganize of
// the analytic single-dose absorption into a real depot, not a change to
// the original’s own (non-standard) exchange kinetics.
dxdt_GUT_HU  = -KA_HU * GUT_HU;
dxdt_CENT_HU =  KA_HU * GUT_HU - (CL_HU/V1_HU) * CENT_HU - (Q_HU/V1_HU) * (CENT_HU - PERI_HU);
dxdt_PERI_HU =  (Q_HU/V1_HU) * (CENT_HU - PERI_HU) - (Q_HU/V2_HU) * PERI_HU;

// ── Drug PK: Anagrelide — Archetype 3 (same non-standard exchange form) ──
dxdt_GUT_ANA  = -KA_ANA * GUT_ANA;
dxdt_CENT_ANA =  KA_ANA * GUT_ANA - (CL_ANA/V1_ANA) * CENT_ANA - (Q_ANA/V1_ANA) * (CENT_ANA - PERI_ANA);
dxdt_PERI_ANA =  (Q_ANA/V1_ANA) * (CENT_ANA - PERI_ANA) - (Q_ANA/V2_ANA) * PERI_ANA;

// ── Drug PK: Ruxolitinib — Archetype 3 (same non-standard exchange form) ──
dxdt_GUT_RUX  = -KA_RUX * GUT_RUX;
dxdt_CENT_RUX =  KA_RUX * GUT_RUX - (CL_RUX/V1_RUX) * CENT_RUX - (Q_RUX/V1_RUX) * (CENT_RUX - PERI_RUX);
dxdt_PERI_RUX =  (Q_RUX/V1_RUX) * (CENT_RUX - PERI_RUX) - (Q_RUX/V2_RUX) * PERI_RUX;

// ── Drug PK: Peg-IFN-α2a — Archetype 3 (same non-standard exchange form) ──
dxdt_GUT_PIFN  = -KA_PIFN * GUT_PIFN;
dxdt_CENT_PIFN =  KA_PIFN * GUT_PIFN - (CL_PIFN/V1_PIFN) * CENT_PIFN - (Q_PIFN/V1_PIFN) * (CENT_PIFN - PERI_PIFN);
dxdt_PERI_PIFN =  (Q_PIFN/V1_PIFN) * (CENT_PIFN - PERI_PIFN) - (Q_PIFN/V2_PIFN) * PERI_PIFN;

// ── The exposed concentration (single redirect site per compound) ──
// The original`s own central-compartment state is already concentration-
// valued (no separate mass/volume split anywhere in its CL/Q terms — see
// et_refactor_notes.md), so C_<STEM> is a direct alias of CENT_<STEM>,
// not CENT_<STEM>/V1_<STEM> (that would double-apply the volume already
// folded into the CL/V1, Q/V1, Q/V2 rate constants).
double C_HU   = CENT_HU;
double C_ANA  = CENT_ANA;
double C_RUX  = CENT_RUX;
double C_PIFN = CENT_PIFN;

// ── PD effect functions (named Hill interface; rename only, same shape/values as original) ──
double EFFECT_HU   = EMAX_HU   * pow(C_HU,   GAMMA_HU)   / (pow(EC50_HU,   GAMMA_HU)   + pow(C_HU,   GAMMA_HU));
double EFFECT_ANA  = EMAX_ANA  * pow(C_ANA,  GAMMA_ANA)  / (pow(EC50_ANA,  GAMMA_ANA)  + pow(C_ANA,  GAMMA_ANA));
double EFFECT_RUX  = EMAX_RUX  * pow(C_RUX,  GAMMA_RUX)  / (pow(EC50_RUX,  GAMMA_RUX)  + pow(C_RUX,  GAMMA_RUX));
double EFFECT_PIFN = EMAX_PIFN * pow(C_PIFN, GAMMA_PIFN) / (pow(EC50_PIFN, GAMMA_PIFN) + pow(C_PIFN, GAMMA_PIFN));

// Aspirin: fixed fraction COX-1 inhibition (binary flag DOSE_ASA>0) — out
// of scope (not one of this file`s four redirect-classified compounds),
// left completely as-is beyond the shared $PARAM build-compat fix above.
double E_ASA = (DOSE_ASA > 0) ? 0.80 : 0.0;   // 80% TXA2 suppression with low-dose ASA

// ── Disease biology ODEs (unchanged) ─────────────────────

// TPO feedback: low PLT → high TPO → drives MK proliferation
double TPO_effect = TPO / (TPO + TPO_ss) * 2.0;  // normalized 0-2

// HSC pool: JAK2 mutation drives clonal expansion
dxdt_HSC = k_HSC_prod * (1.0 + JAK2 * phi_JAK2) * HSC_ss
           - k_HSC_diff * HSC;

// MKP pool: driven by HSC differentiation, suppressed by HU and RUX
dxdt_MKP = k_HSC_diff * HSC * TPO_effect
           - k_MKP_diff * (1.0 + EFFECT_HU + EFFECT_RUX) * MKP;

// Mature MK: from MKP maturation; anagrelide blocks here
dxdt_MK  = k_MKP_diff * MKP
           - k_MK_mature * (1.0 + EFFECT_ANA * 0.6) * MK
           - k_MK_mature * 0.4 * MK;

// Circulating platelets: produced by MK, destroyed normally
double k_PLT_destr_eff = k_PLT_destr * (1.0 + SPL / 10.0);  // splenomegaly increases destruction
dxdt_PLT = k_PLT_prod * k_MK_mature * MK
           - k_PLT_destr_eff * PLT;

// TPO: liver production, consumed by platelets (inverse feedback)
dxdt_TPO = k_TPO_prod
           - k_TPO_elim * TPO
           - k_TPO_abs  * PLT * TPO;

// JAK2 allele burden: clonal expansion offset by pIFN suppression
dxdt_JAK2 = k_JAK2_exp * JAK2 * (1.0 - JAK2)
             - EFFECT_PIFN * 0.015 * JAK2;

// Spleen size: grows with thrombocytosis, shrinks with therapy
double SPL_drive = (PLT / PLT_norm - 1.0);
dxdt_SPL = k_SPL_grow * (SPL_drive > 0 ? SPL_drive : 0.0)
           - k_SPL_shrk * EFFECT_RUX * SPL;

// ── Risk accumulation ODEs ───────────────────────────────
// Thrombosis risk: driven by platelet burden and JAK2 allele load
double thromb_mod = (1.0 + age_gt60 * 0.5) * (1.0 + prior_thrombo * 0.8);
dxdt_RISK_T  = lambda_T * pow(PLT / PLT_norm, alpha_T)
               * (1.0 + delta_JAK2 * JAK2)
               * thromb_mod
               * (1.0 - E_ASA * 0.4);  // aspirin partial protection

// MF transformation risk: driven by JAK2 burden and co-mutations
double mf_mod = 1.0 + ASXL1_pos * delta_ASXL + TP53_pos * 2.0;
dxdt_RISK_MF = lambda_MF * JAK2 * mf_mod;

$TABLE
double PLT_count = PLT;
double JAK2_AB   = JAK2 * 100.0;   // convert to percent
double Spleen    = SPL;
double TPO_level = TPO;
double MK_pool   = MK;

// ELN response classification
double CHR_flag = (PLT <= 400.0) ? 1.0 : 0.0;
double PHR_flag = (PLT <= 600.0 || PLT <= 0.5 * 900.0) ? 1.0 : 0.0;
double CMR_flag = (JAK2 * 100.0 < 1.0) ? 1.0 : 0.0;

// Annual thrombosis hazard (instantaneous from risk accumulator)
double Hazard_thromb = lambda_T * pow(PLT / PLT_norm, alpha_T)
                       * (1.0 + delta_JAK2 * JAK2);
double Hazard_MF = lambda_MF * JAK2;

// Drug concentrations and Hill effects: unlike the original (whose $TABLE
// recomputed Cp_HU/Cp_ANA/Cp_RUX/Cp_pIFN under distinct names because
// $ODE and $TABLE do NOT share a scope in mrgsolve generally), this
// build`s $CAPTURE-driven auto-declaration means C_HU/EFFECT_HU/etc. are
// already visible here as the exact values $ODE last computed — no
// re-declaration needed (confirmed: redeclaring `double C_HU` here is a
// compile-time "redefinition" error under this build, since $CAPTURE
// forward-declares each name once for the whole translation unit).

$CAPTURE PLT_count JAK2_AB Spleen TPO_level MK_pool
         CHR_flag PHR_flag CMR_flag
         Hazard_thromb Hazard_MF
         C_HU C_ANA C_RUX C_PIFN
         EFFECT_HU EFFECT_ANA EFFECT_RUX EFFECT_PIFN
'

mod <- mcode("ET_QSP_refactored", code)

## ── Dosing regimens ─────────────────────────────────────────
## Daily doses converted to appropriate units; drug given as a single
## administration at t=0 via each compound`s GUT_<STEM> depot initial
## amount (F_<STEM>*DOSE_<STEM>, set in $MAIN), matching the original`s own
## single-decaying-dose approximation exactly (see et_refactor_notes.md).
## NOTE (inherited from the original, unchanged by this refactor): this is
## still an approximate single-dose PK approximation, not a repeated daily
## dosing loop — the original`s own header comment already flagged this
## as approximate ("for rigorous PK use depot compartment (done in
## production code)"); the refactor now uses a real depot compartment but
## preserves the same single-administration dosing mechanism, since
## changing *that* would not be a "rename/reorganize" but a behavioral
## change out of scope for this refactor.

mk_ev <- function(dose_hu   = 0,   # mg/day  → split q12h
                  dose_ana  = 0,   # mg/day  → split q12h
                  dose_rux  = 0,   # mg/day  (BID)
                  dose_pifn = 0,   # µg/week → once weekly SC
                  dose_asa  = 0,   # 1 = on, 0 = off
                  n_days    = 365) {

  events <- ev(time = 0, DOSE_HU   = dose_hu / 2 / 1e3,    # mg → g (use µg/mL units via Vd)
                         DOSE_ANA  = dose_ana / 2 / 1e3,
                         DOSE_RUX  = dose_rux / 2 / 1e3,
                         DOSE_PIFN = dose_pifn / 4 / 1e3,  # quarterly approximate
                         DOSE_ASA  = dose_asa)
  events
}

## ── Scenario definitions ────────────────────────────────────
scenarios <- list(
  list(name = "① No Treatment (Observation)",
       dose_hu = 0, dose_ana = 0, dose_rux = 0, dose_pifn = 0, dose_asa = 0,
       color = "#616161"),
  list(name = "② Low-dose Aspirin only",
       dose_hu = 0, dose_ana = 0, dose_rux = 0, dose_pifn = 0, dose_asa = 1,
       color = "#ff8f00"),
  list(name = "③ Hydroxyurea 500 mg/d + ASA",
       dose_hu = 500, dose_ana = 0, dose_rux = 0, dose_pifn = 0, dose_asa = 1,
       color = "#1565c0"),
  list(name = "④ Hydroxyurea 1500 mg/d + ASA",
       dose_hu = 1500, dose_ana = 0, dose_rux = 0, dose_pifn = 0, dose_asa = 1,
       color = "#0288d1"),
  list(name = "⑤ Anagrelide 2 mg/d + ASA",
       dose_hu = 0, dose_ana = 2, dose_rux = 0, dose_pifn = 0, dose_asa = 1,
       color = "#6a1b9a"),
  list(name = "⑥ Ruxolitinib 20 mg/d (BID)",
       dose_hu = 0, dose_ana = 0, dose_rux = 20, dose_pifn = 0, dose_asa = 1,
       color = "#00695c"),
  list(name = "⑦ Peg-IFN-α2a 90 µg/week + ASA",
       dose_hu = 0, dose_ana = 0, dose_rux = 0, dose_pifn = 90, dose_asa = 1,
       color = "#c62828")
)

## ── Run simulations ─────────────────────────────────────────
run_scenario <- function(sc, n_days = 730) {
  params <- c(DOSE_HU   = sc$dose_hu   / 2 / 1e3,
              DOSE_ANA  = sc$dose_ana  / 2 / 1e3,
              DOSE_RUX  = sc$dose_rux  / 2 / 1e3,
              DOSE_PIFN = sc$dose_pifn / 4 / 1e3,
              DOSE_ASA  = sc$dose_asa)

  mrgsim(mod,
         param = params,
         end   = n_days, delta = 1,
         obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(scenario = sc$name,
           color    = sc$color)
}

## Parameters for dose availability in ODE (global params not in $PARAM block)
## Here we approximate by setting them as $PARAM and overriding per run
mod2 <- mod %>%
  param(DOSE_HU = 0, DOSE_ANA = 0, DOSE_RUX = 0, DOSE_PIFN = 0, DOSE_ASA = 0)

results <- purrr::map_dfr(scenarios, function(sc) {
  params_vec <- c(DOSE_HU   = sc$dose_hu   / 1e3,
                  DOSE_ANA  = sc$dose_ana  / 1e3,
                  DOSE_RUX  = sc$dose_rux  / 1e3,
                  DOSE_PIFN = sc$dose_pifn / 1e3,
                  DOSE_ASA  = sc$dose_asa)

  mrgsim(mod2 %>% param(as.list(params_vec)),
         end = 730, delta = 1, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(scenario = sc$name, color = sc$color)
})

## ── Plots ────────────────────────────────────────────────────
pal <- setNames(sapply(scenarios, `[[`, "color"),
                sapply(scenarios, `[[`, "name"))

p_plt <- ggplot(results, aes(time, PLT_count, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(400, 600), linetype = "dashed", color = c("darkgreen","orange")) +
  annotate("text", x = 700, y = 405, label = "CHR threshold (400)", size = 3, hjust = 1) +
  annotate("text", x = 700, y = 605, label = "PHR threshold (600)", size = 3, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "Platelet Count Dynamics (×10⁹/L)",
       x = "Day", y = "Platelets (×10⁹/L)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_jak2 <- ggplot(results, aes(time, JAK2_AB, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "purple") +
  annotate("text", x = 700, y = 2, label = "CMR threshold (<1%)", size = 3, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "JAK2 V617F Allele Burden (%)",
       x = "Day", y = "JAK2 Allele Burden (%)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p_risk <- ggplot(results, aes(time, Hazard_thromb * 365 * 100, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Annual Thrombosis Hazard (%/year)",
       x = "Day", y = "Thrombosis Hazard (%/yr)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p_spleen <- ggplot(results, aes(time, Spleen, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Spleen Size (cm below costal margin)",
       x = "Day", y = "Spleen Size (cm)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

combo_plot <- (p_plt | p_jak2) / (p_risk | p_spleen) +
  plot_annotation(
    title    = "Essential Thrombocythemia — QSP Simulation (7 Scenarios, 2-year)",
    subtitle = "mrgsolve ODE | JAK2/CALR/MPL Mutations · Megakaryopoiesis · Thrombosis · Drug PK/PD",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 10, color = "grey40"))
  )

print(combo_plot)

## ── Response summary table ───────────────────────────────────
response_tbl <- results %>%
  filter(time %in% c(90, 180, 365, 730)) %>%
  group_by(scenario, time) %>%
  summarise(PLT     = round(mean(PLT_count), 1),
            JAK2_AB = round(mean(JAK2_AB),   1),
            CHR     = round(mean(CHR_flag) * 100, 0),
            Spleen  = round(mean(Spleen), 2),
            Thromb_Haz = round(mean(Hazard_thromb * 365 * 100), 2),
            .groups = "drop") %>%
  mutate(Day = paste0("Day ", time)) %>%
  select(Scenario = scenario, Day, PLT, `JAK2 AB(%)` = JAK2_AB,
         `CHR(%)` = CHR, `Spleen(cm)` = Spleen,
         `Thrombo Haz(%/yr)` = Thromb_Haz)

cat("\n=== ET QSP Simulation — Response Summary (refactored) ===\n\n")
print(as.data.frame(response_tbl), row.names = FALSE)

## ── Calibration notes ───────────────────────────────────────
cat("
=== Parameter Calibration Notes (unchanged from the original — see
et_refactor_notes.md for what changed structurally) ===

Hydroxyurea:
  - Harrison 2005 NEJM (PT1): HU 1000-2000mg → CHR ~60% at 1yr
  - Gisslinger 2013 (ANAHYDRET): HU non-inferior to ANA for CHR
  - PK: Tmax ~1h, t½ ~2-3h, Vd ~0.5 L/kg (Gwilt 1998)
  - EC50 adjusted to yield ~600 PLT at 1000mg dose after 90 days

Anagrelide:
  - Harrison 2005 NEJM (PT1): ANA 2mg/d → CHR ~36% at 1yr
  - More arterial events vs HU but fewer venous
  - PK: Tmax ~1h, t½ ~1.3h (Anagrelide Study Group 1997)
  - PDE3A Ki ~36 nM (Psaila 2012 JEM)
  - EC50_ANA tuned to 25-30 ng/mL effective concentration

Ruxolitinib (in ET):
  - Verstovsek 2020 Leukemia: RUX 10mg BID → PLT normalization ~60%
  - RESPONSE-2: spleen volume reduction ≥35% in 40% pts
  - JAK2 allele burden reduction modest (~15-20%)
  - PK: Tmax ~2h, t½ ~3h, Vd ~72L (Shi 2011)

Peg-IFN-α2a:
  - Kiladjian 2013 Haematologica: 45% CHR, JAK2 MR in 18%
  - Quintás-Cardama 2009 J Clin Oncol: 76% CHR (escalating dose)
  - Unique: preferential suppression of JAK2+ HSC clone
  - PK: Tmax ~80h SC, t½ ~80h (Mager 2010)

Thrombosis risk:
  - Low-risk ET: ~1-2%/yr thrombosis
  - High-risk ET: ~4-5%/yr
  - JAK2 V617F: 2x risk vs CALR (Rumi 2014 Blood)
  - Low-dose ASA: ~40% reduction in microvascular events
")
