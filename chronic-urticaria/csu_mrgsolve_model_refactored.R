## ============================================================
## Chronic Spontaneous Urticaria (CSU) — mrgsolve QSP Model (REFACTORED)
## IgE/FcεRI Pathway · Mast Cell Activation · Type-2 Inflammation
## 18 ODE compartments · 7 treatment scenarios
##
## Refactored sibling of csu_mrgsolve_model.R -- pluggable-PK rewrite per
## FORK_WORKFLOW_GUIDE.md Part 2. The original is never edited; see
## csu_refactor_notes.md for the archetype used per compound, the upstream
## build defect found and fixed here (syntax-only), and the verification
## result against the original.
##
## Confirmed compound identities (from code comments/dosing labels in the
## original, not the census file's generic classifier labels):
##   AH  = H1-antihistamine, cetirizine prototype (dosing scenario literally
##         named "Cetirizine 10 mg QD" in the original)
##   OMA = Omalizumab (anti-IgE)
##   DUP = Dupilumab (anti-IL-4Rα)
##   BTK = BTK inhibitor, remibrutinib prototype (comment-confirmed)
##
## Calibration references (unchanged from original):
##   GLACIAL (Omalizumab 300mg q4wk) — Kaplan 2013 JACI
##   ASTERIA I/II (Omalizumab) — Saini 2015 JACI; Maurer 2013 NEJM
##   LIBERTY-CSU CUPID A/B (Dupilumab) — Simpson 2023 NEJM
##   H1-antihistamine PK — Simons 2004 JACI
##   Omalizumab PopPK — Lowe 2009 J Allergy Clin Immunol
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## ----------------------------------------------------------
## 1. Model code block
## ----------------------------------------------------------

csu_code <- '
$PROB
  Chronic Spontaneous Urticaria (CSU) QSP Model -- REFACTORED (pluggable PK)
  IgE/FcεRI-Mast Cell Axis + Type-2 Cytokine Network
  18 ODE compartments

  Refactor notes (see csu_refactor_notes.md for full detail):
  - AH (cetirizine): Archetype 3 minus peripheral (GUT_AH + CENT_AH), linear.
  - OMA (omalizumab): Archetype 3 (GUT/CENT/PERI_OMA) plus a bespoke
    drug-target binding system against the disease’s own free-IgE pool
    (kept exactly as the original modeled it -- see notes for why this is
    bespoke rather than a clean Archetype 4).
  - DUP (dupilumab): Archetype 3 (GUT/CENT/PERI_DUP), linear.
  - BTK (remibrutinib): Archetype 3 minus peripheral (GUT_BTK + CENT_BTK).
  - Every compound’s effect is renamed to the EFFECT_<STEM> Hill interface
    with explicit EMAX_/EC50_/GAMMA_ parameters (GAMMA_BTK and GAMMA_DUP
    added as 1.0 since the original had no explicit Hill exponent for
    those two -- a rename, not a refit; all parameter values are copied
    unchanged from the original).

$PARAM @annotated
  // ---- H1-antihistamine (cetirizine) PK: Archetype 3 minus peripheral ----
  KA_AH    : 1.5    : /h   Cetirizine (AH) absorption rate constant
  CL_AH    : 4.2    : L/h  Cetirizine (AH) clearance
  V1_AH    : 50.0   : L    Cetirizine (AH) volume of distribution
  F_AH     : 0.73   : -    Cetirizine (AH) oral bioavailability
  EMAX_AH  : 0.88   : -    Cetirizine (AH) max inhibition of mast cell degranulation
  EC50_AH  : 0.12   : mg/L Cetirizine (AH) EC50 for H1R blockade
  GAMMA_AH : 1.5    : -    Cetirizine (AH) Hill coefficient (original n_AH)

  // ---- Omalizumab PK: Archetype 3 (2-cpt SC) ----
  KA_OMA   : 0.011  : /h   Omalizumab SC absorption rate constant (~64h tmax)
  CL_OMA   : 0.0049 : L/h  Omalizumab clearance
  V1_OMA   : 3.5    : L    Omalizumab central volume
  V2_OMA   : 3.1    : L    Omalizumab peripheral volume
  Q_OMA    : 0.003  : L/h  Omalizumab inter-compartment clearance
  F_OMA    : 0.62   : -    Omalizumab SC bioavailability
  KON_OMA  : 45.0   : /nM/h Omalizumab-IgE binding (on) rate
  KOFF_OMA : 0.0005 : /h   Omalizumab-IgE dissociation (off) rate

  // ---- Dupilumab PK: Archetype 3 (2-cpt SC) ----
  KA_DUP   : 0.0087 : /h   Dupilumab SC absorption rate constant
  CL_DUP   : 0.0071 : L/h  Dupilumab clearance
  V1_DUP   : 4.8    : L    Dupilumab central volume
  V2_DUP   : 2.9    : L    Dupilumab peripheral volume
  Q_DUP    : 0.0024 : L/h  Dupilumab inter-compartment clearance
  F_DUP    : 0.64   : -    Dupilumab SC bioavailability
  EMAX_DUP_IL4  : 0.96 : - Dupilumab max inhibition of IL-4 signalling
  EMAX_DUP_IL13 : 0.95 : - Dupilumab max inhibition of IL-13 signalling
  EC50_DUP : 0.008  : mg/L Dupilumab EC50 (IL-4Rα), shared by both effects
  GAMMA_DUP: 1.0    : -    Dupilumab Hill coefficient (original had none; rename default)

  // ---- BTK inhibitor (remibrutinib) PK: Archetype 3 minus peripheral ----
  KA_BTK   : 2.1    : /h   BTKi absorption rate constant
  CL_BTK   : 38.0   : L/h  BTKi clearance (high first-pass)
  V1_BTK   : 280.0  : L    BTKi volume of distribution
  F_BTK    : 0.36   : -    BTKi oral bioavailability
  EMAX_BTK : 0.92   : -    BTKi max inhibition of MC activation
  EC50_BTK : 0.045  : mg/L BTKi EC50
  GAMMA_BTK: 1.0    : -    BTKi Hill coefficient (original had none; rename default)

  // ---- IgE biology (disease-side; not compound-specific) ----
  ksyn_IgE : 0.0015 : nM/h Basal IgE synthesis rate
  kdeg_IgE : 0.0050 : /h   Free IgE degradation (t1/2 ~140h)
  IgE0     : 300.0  : nM   Baseline free IgE (elevated CSU ~300 IU/mL equiv)

  // ---- FcεRI / Mast Cell ----
  FcεRI_tot: 1.0    : rel  Total FcεRI expression on mast cells (normalised)
  karm_MC  : 0.08   : /h   Rate of FcεRI arming by IgE
  kdisarm  : 0.02   : /h   Spontaneous FcεRI disarming / receptor turnover
  kact_MC  : 0.15   : /h   Mast cell activation rate (armed MC + autoantigen)
  kdeact_MC: 0.12   : /h   Mast cell deactivation rate
  MC0      : 1.0    : rel  Baseline mast cell priming state
  kprime_IL33: 0.04 : /h/u IL-33 priming of mast cells

  // ---- Histamine PK (skin/plasma) ----
  krel_H   : 2.5    : /h   Histamine release rate from activated MC
  kdeg_Hs  : 0.85   : /h   Histamine degradation in skin
  kdeg_Hp  : 3.2    : /h   Histamine degradation in plasma
  ktrans_H : 0.15   : /h   Skin-to-plasma transfer
  Hist0    : 0.1    : nM   Baseline histamine

  // ---- IL-31 / IL-33 dynamics (disease-side; IL-4/IL-13 network params
  //      below are carried over unused, exactly as in the original -- see
  //      refactor notes) ----
  ksyn_IL4 : 0.008  : nM/h IL-4 synthesis (Th2/ILC2) -- unused, see notes
  kdeg_IL4 : 0.55   : /h   IL-4 degradation -- unused, see notes
  ksyn_IL13: 0.010  : nM/h IL-13 synthesis -- unused, see notes
  kdeg_IL13: 0.48   : /h   IL-13 degradation -- unused, see notes
  ksyn_IL31: 0.005  : nM/h IL-31 synthesis (itch mediator)
  kdeg_IL31: 0.62   : /h   IL-31 degradation
  ksyn_IL33: 0.006  : nM/h IL-33 synthesis
  kdeg_IL33: 0.70   : /h   IL-33 degradation

  // ---- Eosinophil dynamics (declared, unused in any ODE -- see notes) ----
  keo_in   : 0.003  : /h   Eosinophil tissue recruitment -- unused, see notes
  keo_out  : 0.025  : /h   Eosinophil clearance -- unused, see notes
  Eo0      : 1.0    : rel  Baseline skin eosinophil level -- unused, see notes

  // ---- Disease activity & UAS7 mapping ----
  UAS7_max : 42.0   : score Maximum UAS7 score
  UAS7_0   : 30.0   : score Baseline UAS7 score (moderate-severe CSU)
  kUAS_H   : 8.0    : -    UAS7 sensitivity to skin histamine -- unused, see notes
  kUAS_IL31: 5.0    : -    UAS7 sensitivity to IL-31 (itch component) -- unused, see notes
  IgE_norm : 300.0  : nM   Normalisation IgE (= IgE0)

$CMT @annotated
  // PK compartments (renamed to convention; same declaration order as the
  // original so 1-based compartment indices used by dosing events are
  // unchanged between the original and this file)
  GUT_AH    : Cetirizine (AH) GI depot (mg)
  CENT_AH   : Cetirizine (AH) plasma (mg/L equiv)
  GUT_OMA   : Omalizumab SC depot (mg)
  CENT_OMA  : Omalizumab central (mg)
  PERI_OMA  : Omalizumab peripheral (mg)
  GUT_DUP   : Dupilumab SC depot (mg)
  CENT_DUP  : Dupilumab central (mg)
  PERI_DUP  : Dupilumab peripheral (mg)
  GUT_BTK   : BTKi GI depot (mg)
  CENT_BTK  : BTKi plasma (mg/L)

  // PD compartments (unchanged names -- IgE_free/IgE_OMA are the disease’s
  // own free-IgE pool and its omalizumab-bound complex; kept as originally
  // named because IgE_free is also read directly by the FcεRI-arming
  // equation independent of any drug -- see refactor notes, "bespoke")
  IgE_free : Free IgE (nM)
  IgE_OMA  : IgE-Omalizumab complex (nM)
  MC_primed: Armed (IgE-loaded) mast cell index (rel)
  MC_act   : Activated mast cell index (rel)
  Hist_skin: Skin histamine (nM)
  Hist_plasm: Plasma histamine (nM)
  IL31_skin: Skin IL-31 (nM)
  IL33_skin: Skin IL-33 (nM)

$MAIN
  // Steady-state initial conditions (identical to original -- see refactor
  // notes for a finding regarding IgE_free_0 not actually being at the
  // model’s own dynamic equilibrium)
  IgE_free_0 = IgE0;
  MC_primed_0 = MC0;
  Hist_skin_0 = Hist0;
  Hist_plasm_0 = Hist0 * 0.1;
  IL31_skin_0 = ksyn_IL31 / kdeg_IL31;
  IL33_skin_0 = ksyn_IL33 / kdeg_IL33;

$GLOBAL
  // EFFECT_<STEM> declared here (bare, no initializer) so each can be
  // bare-assigned in $ODE and still be $CAPTURE-able. A $PARAM member
  // compiles as a read-only reference inside $ODE in this mrgsolve build,
  // so a value recomputed every timestep cannot be a $PARAM -- same
  // constraint already documented and validated in this corpus (see e.g.
  // copd/copd_refactor_notes.md, breast-cancer/bc_refactor_notes.md).
  double EFFECT_AH;
  double EFFECT_BTK;
  double EFFECT_DUP_IL4;
  double EFFECT_DUP_IL13;

$ODE
  // ----------------------------------------------------------------
  // PK: Cetirizine (AH), Archetype 3 minus peripheral (depot + central)
  C_AH = CENT_AH / V1_AH;   // mg/L
  dxdt_GUT_AH  = -KA_AH * GUT_AH;
  dxdt_CENT_AH =  KA_AH * F_AH * GUT_AH - (CL_AH / V1_AH) * CENT_AH;

  // ----------------------------------------------------------------
  // PK: Omalizumab, Archetype 3 (2-cpt SC), plus bespoke IgE-binding term
  C_OMA = CENT_OMA / V1_OMA;   // mg/L
  dxdt_GUT_OMA  = -KA_OMA * GUT_OMA;
  dxdt_CENT_OMA =  KA_OMA * F_OMA * GUT_OMA - (CL_OMA + Q_OMA) / V1_OMA * CENT_OMA
                   + Q_OMA / V2_OMA * PERI_OMA
                   - KON_OMA * C_OMA * IgE_free + KOFF_OMA * IgE_OMA;
  dxdt_PERI_OMA =  Q_OMA / V1_OMA * CENT_OMA - Q_OMA / V2_OMA * PERI_OMA;

  // ----------------------------------------------------------------
  // PK: Dupilumab, Archetype 3 (2-cpt SC)
  C_DUP = CENT_DUP / V1_DUP;   // mg/L
  dxdt_GUT_DUP  = -KA_DUP * GUT_DUP;
  dxdt_CENT_DUP =  KA_DUP * F_DUP * GUT_DUP - (CL_DUP + Q_DUP) / V1_DUP * CENT_DUP
                   + Q_DUP / V2_DUP * PERI_DUP;
  dxdt_PERI_DUP =  Q_DUP / V1_DUP * CENT_DUP - Q_DUP / V2_DUP * PERI_DUP;

  // ----------------------------------------------------------------
  // PK: BTKi (remibrutinib), Archetype 3 minus peripheral
  C_BTK = CENT_BTK / V1_BTK;   // mg/L
  dxdt_GUT_BTK    = -KA_BTK * GUT_BTK;
  dxdt_CENT_BTK   =  KA_BTK * F_BTK * GUT_BTK - (CL_BTK / V1_BTK) * CENT_BTK;

  // ----------------------------------------------------------------
  // Drug effects -- named Hill interface (rename, not a refit; parameter
  // values copied unchanged from the original’s kinh_*/EC50_*/n_AH)
  EFFECT_AH  = EMAX_AH  * pow(C_AH,  GAMMA_AH)  / (pow(EC50_AH,  GAMMA_AH)  + pow(C_AH,  GAMMA_AH));
  EFFECT_BTK = EMAX_BTK * pow(C_BTK, GAMMA_BTK) / (pow(EC50_BTK, GAMMA_BTK) + pow(C_BTK, GAMMA_BTK));
  // NOTE (preserved defect, see refactor notes + UPSTREAM_ISSUES.md): these
  // two are computed exactly as the original computed E_DUP_IL4/E_DUP_IL13,
  // but -- as in the original -- NEITHER is referenced by any dxdt_ equation
  // below. Dupilumab’s PK is fully modeled but, as coded, has zero effect on
  // disease dynamics. This is carried over unchanged (not fixed) because
  // wiring these in would be new pharmacology the original never had, not a
  // rename.
  EFFECT_DUP_IL4  = EMAX_DUP_IL4  * pow(C_DUP, GAMMA_DUP) / (pow(EC50_DUP, GAMMA_DUP) + pow(C_DUP, GAMMA_DUP));
  EFFECT_DUP_IL13 = EMAX_DUP_IL13 * pow(C_DUP, GAMMA_DUP) / (pow(EC50_DUP, GAMMA_DUP) + pow(C_DUP, GAMMA_DUP));

  // ----------------------------------------------------------------
  // IgE / Omalizumab binding (bespoke -- see refactor notes: the complex
  // IgE_OMA has no clearance term here, exactly as in the original)
  dxdt_IgE_free = ksyn_IgE - kdeg_IgE * IgE_free
                  - KON_OMA * C_OMA * IgE_free + KOFF_OMA * IgE_OMA;
  dxdt_IgE_OMA  = KON_OMA * C_OMA * IgE_free - KOFF_OMA * IgE_OMA;

  // ----------------------------------------------------------------
  // FcεRI arming / mast cell priming
  double fIgE = IgE_free / IgE_norm;   // normalised free IgE fraction
  dxdt_MC_primed = karm_MC * fIgE * (FcεRI_tot - MC_primed - MC_act)
                   - kdisarm * MC_primed
                   - kact_MC * (1.0 - EFFECT_AH) * (1.0 - EFFECT_BTK) * MC_primed
                   + kprime_IL33 * IL33_skin * (FcεRI_tot - MC_primed - MC_act);

  dxdt_MC_act    = kact_MC * (1.0 - EFFECT_AH) * (1.0 - EFFECT_BTK) * MC_primed
                   - kdeact_MC * MC_act;

  // ----------------------------------------------------------------
  // Histamine (skin and plasma)
  dxdt_Hist_skin  = krel_H * MC_act - kdeg_Hs * Hist_skin - ktrans_H * Hist_skin;
  dxdt_Hist_plasm = ktrans_H * Hist_skin - kdeg_Hp * Hist_plasm;

  // ----------------------------------------------------------------
  // Cytokines (type-2 network)
  dxdt_IL31_skin = ksyn_IL31 * (1.0 + 2.0 * MC_act)
                   - kdeg_IL31 * IL31_skin;

  dxdt_IL33_skin = ksyn_IL33 * (1.0 + 1.5 * MC_act)
                   - kdeg_IL33 * IL33_skin;

$TABLE
  // [discoverability] Single-site canonical `double C_<STEM> = <expr>;`
  // initializers -- identical formulas to the $ODE bare recomputes above.
  // mrgsolve does not carry an $ODE-scoped local into $TABLE, so this
  // recompute is what actually makes C_<STEM> available for reporting at
  // every $TABLE-evaluated row; it also gives each compound the literal
  // `double C_<STEM> = ...;` text that downstream tooling regexes for
  // (see FORK_WORKFLOW_GUIDE.md, "What makes a compound’s PK discoverable").
  // This is the fix for the "duplicate concentration site" pattern flagged
  // in the compound census: the original computed the same concentration
  // twice under two DIFFERENT names (e.g. C_AH in $ODE, CONC_AH in $TABLE);
  // here both sites use the one canonical name C_<STEM>.
  double C_AH  = CENT_AH / V1_AH;
  double C_OMA = CENT_OMA / V1_OMA;
  double C_DUP = CENT_DUP / V1_DUP;
  double C_BTK = CENT_BTK / V1_BTK;

  // IgE suppression (%) -- Omalizumab’s disease-facing readout. No
  // standalone EFFECT_OMA Hill term exists (bespoke, see refactor notes):
  // Omalizumab’s action is fully embedded in the IgE mass-balance ODEs
  // above (a state-dependent outcome of binding kinetics), not a separate
  // algebraic function of C_OMA the way EFFECT_AH/EFFECT_BTK are.
  double IgE_suppression = (1.0 - IgE_free / IgE_norm) * 100.0;

  // UAS7 surrogate (identical formula to the original)
  double MC_effect = MC_act / MC0;
  double H_effect  = Hist_skin / Hist0;
  double I_effect  = IL31_skin / (ksyn_IL31 / kdeg_IL31);
  double UAS7 = UAS7_0 * (0.5 * MC_effect * H_effect + 0.3 * I_effect + 0.2);
  if (UAS7 > UAS7_max) UAS7 = UAS7_max;
  if (UAS7 < 0.0)      UAS7 = 0.0;

  double WCU = (UAS7 <= 6.0) ? 1.0 : 0.0;
  double CR  = (UAS7 == 0.0) ? 1.0 : 0.0;

$CAPTURE @annotated
  C_AH            : Cetirizine (H1-antihistamine) plasma concentration, mg/L
  C_OMA           : Omalizumab plasma concentration, mg/L
  C_DUP           : Dupilumab plasma concentration, mg/L
  C_BTK           : BTKi (remibrutinib) plasma concentration, mg/L
  EFFECT_AH       : Cetirizine mast-cell-degranulation inhibition, Hill fraction
  EFFECT_BTK      : BTKi mast-cell-activation inhibition, Hill fraction
  EFFECT_DUP_IL4  : Dupilumab IL-4 signalling inhibition, Hill fraction -- not wired to any disease ODE, see notes
  EFFECT_DUP_IL13 : Dupilumab IL-13 signalling inhibition, Hill fraction -- not wired to any disease ODE, see notes
  IgE_suppression : Free IgE suppression vs baseline, percent
  UAS7            : UAS7 disease-activity surrogate score
  WCU             : Well-controlled urticaria flag, UAS7 <= 6
  CR              : Complete response flag, UAS7 = 0
'

## ----------------------------------------------------------
## 2. Compile model
## ----------------------------------------------------------

mod <- mcode("csu_qsp_refactored", csu_code)

## ----------------------------------------------------------
## 3. Dosing regimens (identical amounts/timing to the original; compartment
##    names updated to the refactored convention -- same 1-based compartment
##    indices as the original since $CMT declaration order is unchanged)
## ----------------------------------------------------------

# H1-antihistamine: cetirizine 10 mg QD oral
dose_AH_std <- ev(amt = 10, cmt = "GUT_AH", ii = 24, addl = 27)   # 4 weeks

# High-dose antihistamine: 40 mg/day (4x10 mg)
dose_AH_high <- ev(amt = 40, cmt = "GUT_AH", ii = 24, addl = 83)  # 12 weeks

# Omalizumab 300 mg q4wk SC
dose_OMA_300 <- ev(amt = 300, cmt = "GUT_OMA", ii = 4*168, addl = 5)  # 6 doses (~24wk)

# Omalizumab 150 mg q4wk SC (lower dose)
dose_OMA_150 <- ev(amt = 150, cmt = "GUT_OMA", ii = 4*168, addl = 5)

# Dupilumab 300 mg q2wk SC (after 600 mg loading)
dose_DUP_LD  <- ev(time = 0,    amt = 600, cmt = "GUT_DUP")
dose_DUP_MD  <- ev(time = 336, amt = 300, cmt = "GUT_DUP", ii = 2*168, addl = 10)
dose_DUP <- c(dose_DUP_LD, dose_DUP_MD)

# BTKi (remibrutinib prototype) 25 mg QD oral
dose_BTK <- ev(amt = 25, cmt = "GUT_BTK", ii = 24, addl = 167)  # 24 weeks

## ----------------------------------------------------------
## 4. Treatment scenarios
## ----------------------------------------------------------

scenarios <- list(
  list(id = 1, name = "No treatment",              ev = ev()),
  list(id = 2, name = "Cetirizine 10 mg QD",        ev = dose_AH_std),
  list(id = 3, name = "High-dose AH 40 mg/day",     ev = dose_AH_high),
  list(id = 4, name = "Omalizumab 300 mg q4wk",     ev = dose_OMA_300),
  list(id = 5, name = "Omalizumab 300 mg + AH",     ev = c(dose_OMA_300, dose_AH_std)),
  list(id = 6, name = "Dupilumab 300 mg q2wk",      ev = dose_DUP),
  list(id = 7, name = "BTKi 25 mg QD",              ev = dose_BTK)
)

## ----------------------------------------------------------
## 5. Simulation function
## ----------------------------------------------------------

sim_scenario <- function(sc, end_h = 24 * 168, delta = 24) {
  out <- mod %>%
    mrgsim(
      events = sc$ev,
      end    = end_h,
      delta  = delta,
      obsonly = TRUE
    ) %>%
    as.data.frame() %>%
    mutate(
      scenario = sc$name,
      time_wk  = time / 168
    )
  out
}

results <- map_dfr(scenarios, sim_scenario)

## ----------------------------------------------------------
## 6. Plots
## ----------------------------------------------------------

# UAS7 over time
p_uas7 <- results %>%
  ggplot(aes(x = time_wk, y = UAS7, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 6,  linetype = "dashed", colour = "darkgreen", linewidth = 0.7) +
  geom_hline(yintercept = 0,  linetype = "dotted",  colour = "steelblue", linewidth = 0.7) +
  annotate("text", x = 22, y = 7.5, label = "WCU threshold (UAS7 <= 6)",
           colour = "darkgreen", size = 3) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 42, 7)) +
  labs(title    = "CSU Disease Activity — UAS7 by Treatment Scenario (refactored)",
       subtitle = "Chronic Spontaneous Urticaria QSP Model",
       x        = "Time (weeks)",
       y        = "UAS7 Score",
       colour   = "Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

print(p_uas7)

# Concentration profiles
p_conc <- results %>%
  ggplot(aes(x = time_wk, y = C_OMA, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Omalizumab Plasma Concentration",
       x      = "Time (weeks)",
       y      = "C_OMA (mg/L)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_conc)

# Mast cell activation
p_mc <- results %>%
  ggplot(aes(x = time_wk, y = MC_act, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Mast Cell Activation Index",
       x      = "Time (weeks)",
       y      = "MC Activation (rel.)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_mc)

# Skin histamine
p_hist <- results %>%
  ggplot(aes(x = time_wk, y = Hist_skin, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Skin Histamine Concentration",
       x      = "Time (weeks)",
       y      = "Skin Histamine (nM)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_hist)

## ----------------------------------------------------------
## 7. Summary table at key time points
## ----------------------------------------------------------

key_wks <- c(4, 12, 24)

summary_tbl <- results %>%
  filter(round(time_wk, 1) %in% key_wks) %>%
  group_by(scenario, time_wk) %>%
  slice_tail(n = 1) %>%
  summarise(
    UAS7_mean       = round(mean(UAS7), 1),
    WCU_pct         = round(mean(WCU) * 100, 1),
    IgE_sup_pct     = round(mean(IgE_suppression), 1),
    MC_act_rel      = round(mean(MC_act), 3),
    .groups = "drop"
  )

print(summary_tbl)
