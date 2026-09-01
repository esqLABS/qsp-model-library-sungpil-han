## ============================================================
##  Waldenström's Macroglobulinemia (WM) — mrgsolve QSP Model
##  Filename : wm_mrgsolve_model_refactored.R
##  Author   : Claude Code Routine (CCR)
##  Date     : 2026-06-27 (original) / refactored 2026-08-30
## ============================================================
##
##  PK/PD REFACTOR (fork-only sibling of wm_mrgsolve_model.R)
##  ────────────────────────────────────────────────────────────
##  Per FORK_WORKFLOW_GUIDE.md Part 2, all four of this file's
##  compounds are refactored to the guide's naming convention and
##  a named Hill (EFFECT_<STEM>) interface: Ibrutinib (IBR),
##  Rituximab (RTX), Venetoclax (VEN), Zanubrutinib (ZAN). The
##  disease model (BTK_occ, NFkB, LPC, PC, IgM, Hgb, Visc, BCL2,
##  CD20, Apop, NK) and the untouched compounds (bortezomib,
##  bendamustine) are unchanged apart from reading the new
##  EFFECT_IBR/EFFECT_ZAN/EFFECT_RTX/EFFECT_VEN names in place of
##  the original's local IBR_effect/ZAN_effect/RTX_CD20_inh/
##  VEN_BCL2_inhib variables (same arithmetic, same order).
##
##  Two pre-existing build/runtime defects were found in the
##  original (present regardless of this refactor) and are fixed
##  syntax-only, disclosed, in this file per the guide's settled
##  policy for "when the original doesn't compile at all" — see
##  wm_refactor_notes.md and translations/UPSTREAM_ISSUES.md #72:
##   1. `BMInf = fmin(...); dxdt_BMInf = 0;` assigns directly to a
##      $CMT compartment, which mrgsolve 2.0.1 rejects at compile
##      time ("assignment of read-only reference"). BMInf is
##      demoted to a local double in $ODE (algebraic, as the
##      original's own comment already said), recomputed from
##      LPC/PC in $TABLE for BMInf_pct — same formula, same value.
##   2. `_F(GUT) = value;` (bioavailability) compiles but crashes
##      the solver at runtime (signal 6) for every scenario,
##      dosed or not. Replaced with the modern `F_<CMT> = value;`
##      idiom, confirmed to apply identically.
##
##  COMPARTMENTS (19 ODE compartments; BMInf demoted to algebraic, see above)
##  ────────────────────────────────────────────────────────────
##  Drug PK (7)
##    1. GUT_IBR  — Ibrutinib gut (oral depot)              [was IBR_gut]
##    2. CENT_IBR — Ibrutinib central (plasma Cp)            [was IBR_C]
##    3. GUT_ZAN  — Zanubrutinib gut                          [was ZAN_gut]
##    4. CENT_ZAN — Zanubrutinib central                      [was ZAN_C]
##    5. CENT_RTX — Rituximab central                         [was RTX_C]
##    6. GUT_VEN  — Venetoclax gut                             [was VEN_gut]
##    7. CENT_VEN — Venetoclax central                         [was VEN_C]
##  PD / Disease (12)
##    8.  BTK_occ  — BTK occupancy (covalent, fraction 0-1; shared IBR+ZAN target)
##    9.  NFkB     — NF-κB activity (AU, driven by MYD88)
##   10.  LPC      — Lymphoplasmacytic cells (× 10⁹ cells, BM)
##   11.  PC       — IgM-secreting plasma cells (× 10⁹)
##   12.  IgM      — Serum IgM concentration (g/L)
##   13.  Hgb      — Hemoglobin (g/dL)
##   14.  Visc     — Serum viscosity (cP)
##   15.  BCL2     — Effective BCL-2 anti-apoptotic activity (AU)
##   16.  CD20     — Surface CD20 (× baseline, rituximab target)
##   17.  Apop     — Apoptosis rate modifier (AU)
##   18.  Protsm   — Proteasome activity (AU, bortezomib target — untouched)
##   19.  NK       — NK cell count (×10⁶/L, for ADCC)
##
##  TREATMENT SCENARIOS (7) — unchanged from the original
##  ────────────────────────────────────────────────────────────
##   1. Watch & Wait (natural history, symptomatic WM)
##   2. Ibrutinib monotherapy (420 mg/day — iNNOVATOR)
##   3. Ibrutinib + Rituximab (iR — INNOVATE trial)
##   4. Zanubrutinib monotherapy (ASPEN trial)
##   5. Rituximab-Bendamustine (R-Benda — 1st line)
##   6. Bortezomib + Rituximab + Dexamethasone (BDR)
##   7. Venetoclax (salvage / BTK-resistant WM)
##
##  CALIBRATION DATA (major clinical trials) — unchanged, see original
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ── Model code ────────────────────────────────────────────────
code <- '
$PARAM
  // ── Ibrutinib (IBR) PK — archetype 3, no peripheral (depot + central, linear)
  IBR_dose = 0,    // dead in the original (never read in $ODE/$MAIN/$TABLE; dosing is via ev() in the R driver) — preserved unrenamed, see refactor notes
  KA_IBR   = 1.2,  // absorption rate constant (1/h)              [was ka_ibr]
  F_IBR    = 0.10, // bioavailability ~10 % (fed state ~10-15 %)  [was F_ibr]
  CL_IBR   = 73.0, // clearance (L/h) — CYP3A4 dominated          [was CL_ibr]
  V1_IBR   = 820,  // central volume of distribution (L)         [was V_ibr]
  MW_IBR   = 440,  // ibrutinib molecular weight (nM conversion)  [was MW_ibr]
  EC50_IBR = 0.5,  // ibrutinib EC50 at BTK (nM)                  [was EC50_ibr_btk]
  EMAX_IBR = 1,    // new: original ratio saturates at 1 (no explicit Emax) — see notes
  GAMMA_IBR= 1,    // new: original had no explicit Hill coefficient — see notes

  // ── Zanubrutinib (ZAN) PK — same archetype as ibrutinib
  ZAN_dose = 0,    // dead parameter, same as IBR_dose — see refactor notes
  KA_ZAN   = 1.5,  // [was ka_zan]
  F_ZAN    = 0.60, // [was F_zan]
  CL_ZAN   = 29.0, // [was CL_zan]
  V1_ZAN   = 520,  // [was V_zan]
  MW_ZAN   = 471,  // [was MW_zan]
  EC50_ZAN = 0.30, // zanubrutinib EC50 at BTK (nM)                [was EC50_zan_btk]
  EMAX_ZAN = 1,    // new, see notes
  GAMMA_ZAN= 1,    // new, see notes

  // ── Rituximab (RTX) PK — archetype 1 (single compartment, linear)
  // plus a bespoke, disclosed extra loss term carried over unchanged; see notes
  RTX_dose = 0,    // dead parameter, see refactor notes
  k_RTX_on = 0.04, // dead parameter (never read anywhere), see refactor notes
  CL_RTX   = 0.21, // clearance (L/h)
  V1_RTX   = 4.5,  // central volume (L)                          [was V_RTX]
  K12_RTX  = 0.02, // bespoke: extra elimination-like term with no destination
                   // compartment or return flux (no naming-convention slot fits
                   // this — see refactor notes)                  [was k12_RTX]
  k21_RTX  = 0.010,// dead parameter (declared, never used anywhere in $ODE)
  EC50_RTX = 0.05, // rituximab EC50 for CD20 depletion (mg/L)     [was EC50_RTX_CD20]
  EMAX_RTX = 1,    // new, see notes
  GAMMA_RTX= 1,    // new, see notes

  // ── Venetoclax (VEN) PK — same archetype as ibrutinib/zanubrutinib
  VEN_dose = 0,    // dead parameter, see refactor notes
  KA_VEN   = 0.40, // [was ka_ven]
  F_VEN    = 0.72, // with food                                    [was F_ven]
  CL_VEN   = 12.0, // [was CL_ven]
  V1_VEN   = 256,  // [was V_ven]
  MW_VEN   = 868,  // venetoclax molecular weight                  [was VEN_MW]
  EC50_VEN = 0.5,  // venetoclax EC50 for BCL-2 inhibition (µM)
  EMAX_VEN = 1,    // new, see notes
  GAMMA_VEN= 1,    // new, see notes

  // ── BTK pharmacodynamics (shared disease compartment, driven by IBR+ZAN) ─
  kout_BTK     = 0.008, // BTK resynthesis rate (1/h; protein T½ ~3 d)

  // ── NF-κB dynamics ────────────────────────────────────
  NFkB_base    = 1.0,
  MYD88_drive  = 0.60,
  kBTK_NFkB    = 0.35,
  EC50_BOR_NFkB= 50.0,
  kout_NFkB    = 0.50,

  // ── Tumor (LPC + PC) dynamics ─────────────────────────
  kel0_LPC   = 0.002,
  kprolif_LPC= 0.0045,
  NFkB_kLPC  = 0.30,
  LPC0       = 50,
  kconv_LPC  = 0.00015,
  kel0_PC    = 0.001,
  PC0        = 10,
  KMAX_BM    = 200,

  // ── IgM secretion & clearance ─────────────────────────
  ksec_IgM   = 0.06,
  kel_IgM    = 0.004,
  IgM0       = 25.0,

  // ── Hemoglobin ────────────────────────────────────────
  Hgb0       = 9.5,
  kprod_Hgb  = 0.030,
  kel_Hgb    = 0.0028,
  BMInf_Hgb  = 0.60,

  // ── Serum viscosity ───────────────────────────────────
  Visc_base  = 1.5,
  kIgM_visc  = 0.080,
  Visc_exp   = 1.4,

  // ── BCL-2 / apoptosis ─────────────────────────────────
  BCL2_base  = 1.0,
  kNFkB_BCL2 = 0.40,
  kout_BCL2  = 0.10,

  // ── CD20 / ADCC (disease side) ────────────────────────
  kout_CD20  = 0.004,
  kADCC      = 0.025,

  // ── Proteasome (bortezomib — untouched, out of this refactor’s scope) ─
  BOR_dose   = 0,
  BOR_Cp     = 0,
  EC50_BOR   = 10.0,
  kout_Prot  = 0.050,
  kprot_kill = 0.012,

  // ── NK cells ──────────────────────────────────────────
  NK0        = 100,
  kprod_NK   = 0.008,
  kel_NK     = 0.006,

  // ── Bendamustine (untouched, out of this refactor’s scope) ───
  BENDA_kLPC = 0.010,
  BENDA_flag = 0

$CMT
  GUT_IBR CENT_IBR
  GUT_ZAN CENT_ZAN
  CENT_RTX
  GUT_VEN CENT_VEN
  BTK_occ
  NFkB
  LPC PC
  IgM
  Hgb
  Visc
  BCL2
  CD20
  Apop
  Protsm
  NK

$MAIN
  // Initial conditions
  // build-compat fix: "_F(CMT) = value;" is not valid/safe mrgsolve 2.0.1
  // syntax here -- it compiles (POST /model_manifest succeeds) but crashes
  // the solver at runtime (POST /run_simulation, signal 6 / "munmap_chunk():
  // invalid pointer") for every scenario, dosed or not. Replaced with the
  // modern "F_<CMT> = value;" idiom, which applies identically (confirmed:
  // steady-state CENT_IBR/ZAN/VEN scale exactly by F_IBR/F_ZAN/F_VEN as
  // expected) -- see refactor notes / UPSTREAM_ISSUES.md #72
  F_GUT_IBR = F_IBR;
  F_GUT_ZAN = F_ZAN;
  F_GUT_VEN = F_VEN;

  LPC_0    = LPC0;
  PC_0     = PC0;
  IgM_0    = IgM0;
  Hgb_0    = Hgb0;
  Visc_0   = Visc_base + kIgM_visc * pow(IgM0, Visc_exp);
  BCL2_0   = BCL2_base;
  CD20_0   = 1.0;
  Apop_0   = 0.0;
  Protsm_0 = 100.0; // 100 AU = fully active
  NK_0     = NK0;
  BTK_occ_0= 0.0;
  NFkB_0   = NFkB_base + MYD88_drive;

$ODE
  // ── Drug PK ───────────────────────────────────────────────

  // Ibrutinib (IBR): archetype 3 without a peripheral compartment
  double C_IBR_NGML = CENT_IBR / V1_IBR * 1000.0;       // ng/mL diagnostic (matches original’s $CAPTURE unit)
  double C_IBR      = C_IBR_NGML / MW_IBR * 1000.0;     // nM — the exposed concentration (PD reads this)
  dxdt_GUT_IBR  = -KA_IBR * GUT_IBR;
  dxdt_CENT_IBR =  KA_IBR * GUT_IBR - (CL_IBR / V1_IBR) * CENT_IBR;

  // Zanubrutinib (ZAN): same archetype as ibrutinib
  double C_ZAN_NGML = CENT_ZAN / V1_ZAN * 1000.0;
  double C_ZAN      = C_ZAN_NGML / MW_ZAN * 1000.0;     // nM
  dxdt_GUT_ZAN  = -KA_ZAN * GUT_ZAN;
  dxdt_CENT_ZAN =  KA_ZAN * GUT_ZAN - (CL_ZAN / V1_ZAN) * CENT_ZAN;

  // Rituximab (RTX): archetype 1 (single compartment, linear elimination)
  // plus a bespoke, disclosed extra loss term (K12_RTX) carried over
  // unchanged from the original — see refactor notes
  double C_RTX = CENT_RTX / V1_RTX; // mg/L — the exposed concentration
  dxdt_CENT_RTX = -(CL_RTX / V1_RTX) * CENT_RTX - K12_RTX * CENT_RTX;

  // Venetoclax (VEN): same archetype as ibrutinib/zanubrutinib
  double C_VEN = (CENT_VEN / V1_VEN * 1000.0) / MW_VEN * 1000.0; // µM
  dxdt_GUT_VEN  = -KA_VEN * GUT_VEN;
  dxdt_CENT_VEN =  KA_VEN * GUT_VEN - (CL_VEN / V1_VEN) * CENT_VEN;

  // ── BTK Occupancy (shared disease compartment driven by IBR + ZAN) ──
  // IBR and ZAN compete for BTK; fraction occupied grows then
  // declines as new BTK protein is synthesised
  double EFFECT_IBR = EMAX_IBR * pow(C_IBR, GAMMA_IBR) / (pow(EC50_IBR, GAMMA_IBR) + pow(C_IBR, GAMMA_IBR));
  double EFFECT_ZAN = EMAX_ZAN * pow(C_ZAN, GAMMA_ZAN) / (pow(EC50_ZAN, GAMMA_ZAN) + pow(C_ZAN, GAMMA_ZAN));
  double BTK_input  = fmax(EFFECT_IBR, EFFECT_ZAN); // dominant inhibitor, same as original
  dxdt_BTK_occ = BTK_input * (1.0 - BTK_occ) - kout_BTK * BTK_occ;

  // ── NF-κB activity ────────────────────────────────────────
  // MYD88-L265P = constitutive drive; BTK modulates amplitude
  // Bortezomib stabilises IκBα → reduces NF-κB
  double BTK_contrib = kBTK_NFkB * (1.0 - BTK_occ);
  double BOR_inhib_NF= EC50_BOR_NFkB / (EC50_BOR_NFkB + (100.0 - Protsm));
  double NFkB_input  = NFkB_base + MYD88_drive + BTK_contrib;
  dxdt_NFkB = NFkB_input * BOR_inhib_NF - kout_NFkB * NFkB;

  // ── BCL-2 (NF-κB → BCL2; Venetoclax inhibits) ───────────
  double EFFECT_VEN = EMAX_VEN * pow(C_VEN, GAMMA_VEN) / (pow(EC50_VEN, GAMMA_VEN) + pow(C_VEN, GAMMA_VEN));
  dxdt_BCL2 = kNFkB_BCL2 * NFkB * (1.0 - EFFECT_VEN) - kout_BCL2 * BCL2;

  // ── Apoptosis composite ───────────────────────────────────
  // Increases with BTK inhibition, BCL-2 inhibition, proteasome inhibition
  double apo_drug = BTK_occ * 0.4
                  + EFFECT_VEN * 0.5
                  + (100.0 - Protsm) / 100.0 * kprot_kill / kel0_LPC;
  dxdt_Apop = apo_drug - 0.5 * Apop;

  // ── CD20 (Rituximab depletes surface CD20) ────────────────
  double EFFECT_RTX = EMAX_RTX * pow(C_RTX, GAMMA_RTX) / (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX));
  dxdt_CD20 = kout_CD20 * (1.0 - CD20) - EFFECT_RTX * CD20;

  // ── NK cells (ADCC capacity) ──────────────────────────────
  dxdt_NK = kprod_NK * NK0 - kel_NK * NK;

  // ── Proteasome activity (bortezomib inhibits) ─────────────
  // BOR_Cp is set externally per cycle event; simplified
  double BOR_inhib = BOR_Cp / (BOR_Cp + EC50_BOR);
  dxdt_Protsm = kout_Prot * (100.0 - Protsm) - BOR_inhib * Protsm;

  // ── BM Infiltration fraction (build-compat fix: demoted from a $CMT
  // state to a local double — mrgsolve 2.0.1 rejects direct assignment to
  // a compartment name; see refactor notes / UPSTREAM_ISSUES.md #72) ───
  double total_tumor = LPC + PC;
  double BMInf = fmin(total_tumor / KMAX_BM, 1.0);

  // ── LPC (lymphoplasmacytic cells) dynamics ────────────────
  // Proliferation driven by NF-κB; apoptosis enhanced by drugs
  double LPC_prolif = kprolif_LPC * (1.0 + NFkB_kLPC * NFkB) *
                      LPC * (1.0 - total_tumor / KMAX_BM);
  double LPC_apop   = (kel0_LPC + apo_drug) * LPC;
  double LPC_ADCC   = kADCC * NK * CD20 * EFFECT_RTX * LPC / (LPC + 1.0);
  double LPC_BENDA  = BENDA_flag * BENDA_kLPC * LPC;
  double LPC_conv   = kconv_LPC * LPC; // LPC → PC differentiation
  dxdt_LPC = LPC_prolif - LPC_apop - LPC_ADCC - LPC_BENDA - LPC_conv;

  // ── Plasma cell dynamics ──────────────────────────────────
  double PC_apop  = (kel0_PC + apo_drug * 0.6) * PC;
  double PC_ADCC  = kADCC * 0.5 * NK * CD20 * EFFECT_RTX * PC / (PC + 1.0);
  double PC_BENDA = BENDA_flag * BENDA_kLPC * 0.7 * PC;
  dxdt_PC = LPC_conv - PC_apop - PC_ADCC - PC_BENDA;

  // ── Serum IgM ─────────────────────────────────────────────
  double IgM_prod = ksec_IgM * PC;
  double IgM_elim = kel_IgM * IgM;
  dxdt_IgM = IgM_prod - IgM_elim;

  // ── Hemoglobin ───────────────────────────────────────────
  double suppression = 1.0 - BMInf_Hgb * BMInf;
  double Hgb_prod = kprod_Hgb * Hgb0 * suppression;
  dxdt_Hgb = Hgb_prod - kel_Hgb * Hgb;

  // ── Serum Viscosity ───────────────────────────────────────
  double Visc_new = Visc_base + kIgM_visc * pow(fmax(IgM, 0.0), Visc_exp);
  dxdt_Visc = (Visc_new - Visc) * 0.5; // lag toward equilibrium

$TABLE
  double IgM_gL     = IgM;
  double Hgb_gdL    = Hgb;
  // build-compat: BMInf is no longer a $CMT state (see $ODE); recomputed
  // here from LPC/PC directly, identical formula, same numeric result
  double BMInf_pct  = fmin((LPC + PC) / KMAX_BM, 1.0) * 100.0;
  double Visc_cP    = Visc;
  double BTK_pct    = BTK_occ * 100.0;
  double NFkB_AU    = NFkB;
  double BCL2_AU    = BCL2;
  // aliases of the $ODE locals (not recomputed fresh from state here) --
  // matches the original’s own $TABLE pattern exactly (it aliased
  // Cp_ibr_ngmL/Cp_zan_ngmL/Cp_RTX/Cp_ven_uM the same way), including the
  // original’s own quirk of reporting the pre-dose $ODE-local value on the
  // rare row where an observation time coincides exactly with a dose time
  // (confirmed needed for an exact match -- see refactor notes)
  double Cp_IBR     = C_IBR_NGML;
  double Cp_ZAN     = C_ZAN_NGML;
  double Cp_RTX_mgl = C_RTX;
  double Cp_VEN_uM  = C_VEN;
  double LPC_cells  = LPC;
  double PC_cells   = PC;
  double NK_cells   = NK;
  double Protsm_AU  = Protsm;
  double CD20_frac  = CD20;
  double Apop_AU    = Apop;
  // Hyperviscosity flag (>3.5 cP = symptomatic threshold)
  double HVS_flag   = (Visc > 3.5) ? 1.0 : 0.0;
  // Response category (based on IgM% change from baseline)
  double IgM_change_pct = (IgM0 > 0) ? (IgM - IgM0) / IgM0 * 100.0 : 0.0;

$CAPTURE
  Cp_IBR Cp_ZAN Cp_RTX_mgl Cp_VEN_uM
  C_IBR C_ZAN C_RTX C_VEN
  EFFECT_IBR EFFECT_ZAN EFFECT_RTX EFFECT_VEN
  BTK_pct NFkB_AU BCL2_AU
  LPC_cells PC_cells IgM_gL Hgb_gdL
  BMInf_pct Visc_cP NK_cells Protsm_AU CD20_frac Apop_AU
  HVS_flag IgM_change_pct
'

mod <- mcode("WM_QSP_refactored", code)

## ============================================================
##  Helper: build dosing event tables
## ============================================================

# Ibrutinib 420 mg/day continuously
ibr_ev <- function(days = 730) {
  ev(cmt = "GUT_IBR",
     amt = 420 * 0.10,  # dose × F already factored in F_GUT_IBR
     ii  = 24, addl = days - 1, time = 0)
}

# Zanubrutinib 160 mg BID
zan_ev <- function(days = 730) {
  ev(cmt = "GUT_ZAN",
     amt = 160 * 0.60,
     ii  = 12, addl = days * 2 - 1, time = 0)
}

# Rituximab 375 mg/m² (assuming BSA 1.8 m²) q4 weeks × 6 then q8 weeks
rtx_ev <- function(n_cycles = 6, start = 0) {
  times <- start + seq(0, by = 28 * 24, length.out = n_cycles)
  ev(cmt = "CENT_RTX",
     amt = 375 * 1.8,    # total dose in mg
     time = times)
}

## ============================================================
##  Simulation scenarios
## ============================================================

run_scenario <- function(scenario, days = 730) {
  end_h <- days * 24
  delta <- 24   # hourly output → daily

  base_param <- list(BENDA_flag = 0, BOR_Cp = 0)

  e <- switch(scenario,
    "Watch_Wait" = ev(time = 0, amt = 0, cmt = "GUT_IBR"),

    "Ibrutinib"  = ibr_ev(days),

    "iR"         = ev_c(ibr_ev(days), rtx_ev(6, 0)),

    "Zanubrutinib" = zan_ev(days),

    "R_Benda" = {
      # R-Bendamustine: 6 cycles × 28d; Benda 90mg/m² d1-2
      rtx_part  <- rtx_ev(6, 0)
      # Simplified: toggle BENDA_flag via model param changes per cycle
      rtx_part
    },

    "BDR" = {
      rtx_ev(6, 0)
    },

    "Venetoclax" = {
      # Ramp: wk1=20 mg, wk2=50 mg, wk3=100 mg, wk4=200 mg, wk5+=400 mg
      ramp_amt <- c(20, 50, 100, 200, 400) * 0.72
      ramp_wk  <- c(0, 7, 14, 21, 28)
      ev(cmt  = "GUT_VEN",
         amt  = rep(ramp_amt, times = c(7, 7, 7, 7, days - 28)),
         time = unlist(mapply(function(a, s) s * 24 + seq(0, by = 24, length.out = ifelse(a == tail(ramp_amt,1), days-28, 7)),
                              ramp_amt, ramp_wk, SIMPLIFY = FALSE)))
    }
  )

  # Extra scenario-specific parameter modifications
  extra_p <- list()
  if (scenario == "R_Benda") {
    extra_p <- list(BENDA_flag = 0)  # simplified; toggle per cycle in full sim
  }
  if (scenario == "BDR") {
    extra_p <- list(BOR_Cp = 200)  # representative exposure during cycle
  }

  mrgsim(mod,
    events = e,
    param  = extra_p,
    end    = end_h,
    delta  = delta,
    obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(scenario = scenario,
           day      = time / 24)
}

scenarios <- c("Watch_Wait", "Ibrutinib", "iR", "Zanubrutinib",
                "R_Benda", "BDR", "Venetoclax")

results <- bind_rows(lapply(scenarios, function(s) {
  tryCatch(run_scenario(s), error = function(e) NULL)
}))

## ============================================================
##  Plot function
## ============================================================

scenario_colors <- c(
  "Watch_Wait"   = "#7F8C8D",
  "Ibrutinib"    = "#2980B9",
  "iR"           = "#1A5276",
  "Zanubrutinib" = "#8E44AD",
  "R_Benda"      = "#27AE60",
  "BDR"          = "#D35400",
  "Venetoclax"   = "#C0392B"
)
scenario_labels <- c(
  "Watch_Wait"   = "Watch & Wait",
  "Ibrutinib"    = "Ibrutinib 420 mg/d",
  "iR"           = "Ibrutinib + Rituximab",
  "Zanubrutinib" = "Zanubrutinib 160 mg BID",
  "R_Benda"      = "R-Bendamustine",
  "BDR"          = "Bortezomib+Rituximab+Dex",
  "Venetoclax"   = "Venetoclax (salvage)"
)

plot_panel <- function(var, ylab, title, yint = NULL) {
  p <- ggplot(results, aes(x = day, y = .data[[var]],
                            color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1.1, alpha = 0.9) +
    scale_color_manual(values = scenario_colors, labels = scenario_labels) +
    scale_linetype_manual(values = c("solid","dashed","dotdash",
                                     "longdash","twodash","dotted","solid"),
                          labels = scenario_labels) +
    labs(x = "Day", y = ylab, title = title,
         color = "Scenario", linetype = "Scenario") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"))
  if (!is.null(yint))
    p <- p + geom_hline(yintercept = yint, linetype = "dotted",
                         color = "#E74C3C", linewidth = 0.7)
  p
}

p_IgM  <- plot_panel("IgM_gL",      "IgM (g/L)",       "Serum IgM over Time", yint = 7)
p_Hgb  <- plot_panel("Hgb_gdL",     "Hemoglobin (g/dL)","Hemoglobin over Time", yint = 10)
p_LPC  <- plot_panel("LPC_cells",   "LPC (×10⁹)",       "Lymphoplasmacytic Cells")
p_Visc <- plot_panel("Visc_cP",     "Viscosity (cP)",   "Serum Viscosity", yint = 3.5)
p_BTK  <- plot_panel("BTK_pct",     "BTK Occupancy (%)", "BTK Target Occupancy")
p_BMInf<- plot_panel("BMInf_pct",   "BM Infiltration (%)","Bone Marrow Infiltration")

combined <- (p_IgM | p_Hgb) / (p_LPC | p_Visc) / (p_BTK | p_BMInf) +
  plot_annotation(
    title   = "Waldenström's Macroglobulinemia — QSP Simulation",
    subtitle = "7 Treatment Scenarios · 730-day time horizon",
    theme   = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 11))
  )

print(combined)

## ============================================================
##  Response rate summary at 12 months
## ============================================================
response_summary <- results %>%
  filter(day >= 360 & day <= 362) %>%
  group_by(scenario) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    IgM_reduction_pct = (IgM0_val - IgM_gL) / IgM0_val * 100,
    Response = case_when(
      IgM_reduction_pct >= 90 ~ "VGPR",
      IgM_reduction_pct >= 50 ~ "PR",
      IgM_reduction_pct >= 25 ~ "MR",
      IgM_reduction_pct >= 0  ~ "SD",
      TRUE                    ~ "PD"
    ),
    Hgb_normalized = Hgb_gdL >= 11.0
  ) %>%
  select(scenario, IgM_gL, IgM_reduction_pct, Hgb_gdL,
         BMInf_pct, Visc_cP, Response, Hgb_normalized)

# Replace IgM0_val reference with constant
response_summary <- results %>%
  filter(day >= 360 & day <= 362) %>%
  group_by(scenario) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    IgM0_const = 25.0,
    IgM_reduction_pct = (IgM0_const - IgM_gL) / IgM0_const * 100,
    Response = case_when(
      IgM_reduction_pct >= 90 ~ "VGPR (≥90% IgM reduction)",
      IgM_reduction_pct >= 50 ~ "PR  (≥50% IgM reduction)",
      IgM_reduction_pct >= 25 ~ "MR  (≥25% IgM reduction)",
      IgM_reduction_pct >= 0  ~ "SD",
      TRUE                    ~ "PD"
    )
  ) %>%
  select(Scenario = scenario,
         IgM_gL, IgM_reduction_pct, Hgb_gdL, BMInf_pct, Visc_cP, Response)

cat("\n=== Response Summary at 12 months (360 d) ===\n")
print(response_summary, n = Inf)

## ============================================================
##  Clinical trial calibration targets
## ============================================================
cat("\n=== Clinical Trial Calibration Targets ===\n")
tibble::tribble(
  ~Trial,         ~Regimen,          ~ORR,   ~VGPR_CR, ~mPFS,
  "iNNOVATOR",    "Ibrutinib mono",  "91.5%","30.4%",  "69% at 2y",
  "INNOVATE",     "Ibru+Rituximab",  "92%",  "43%",    "82% at 30mo",
  "ASPEN",        "Zanubrutinib",    "93.7%","28.4%",  "84% at 18mo",
  "Rummel 2013",  "R-Benda",         "96%",  "44%",    "69 mo",
  "Dimopoulos13", "BDR",             "83%",  "22%",    "43 mo",
  "Castillo 2018","Venetoclax",      "84%",  "36%",    "Not reached"
) %>% print(n = Inf)
