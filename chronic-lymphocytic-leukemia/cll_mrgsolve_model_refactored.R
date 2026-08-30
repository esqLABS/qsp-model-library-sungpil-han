# ============================================================================
# Chronic Lymphocytic Leukemia (CLL) — QSP mrgsolve Model (PK/PD REFACTORED)
# ============================================================================
# Refactored sibling of cll_mrgsolve_model.R per FORK_WORKFLOW_GUIDE.md Part 2.
# The original is left completely untouched; this file only reorganizes and
# renames each of the three compounds' PK/PD into the guide's pluggable
# naming convention (GUT_/CENT_/PERI_/REC_FREE_/COMPLEX_/C_/EFFECT_<STEM>).
# No disease biology, parameter value, or dosing schedule was changed.
#
# Compartments: 18 ODE states (same 18 as the original, renamed where the
#               naming convention applies)
# Scenarios   : 6 treatment regimens (unchanged from the original)
# Calibration : CLL14 (Fischer 2019 NEJM), RESONATE-2 (Burger 2015 NEJM),
#               MURANO (Seymour 2018 NEJM), SEQUOIA (Shadman 2023 NEJM)
#
# See cll_refactor_notes.md for: per-compound archetype classification,
# the renaming map, the pre-existing mrgsolve-2.0.1 build defects found and
# fixed here (syntax-only, disclosed, logged as UPSTREAM_ISSUES.md #45), and
# the qspserver mrgsolve_api verification results.
# ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ── Model Definition ─────────────────────────────────────────────────────────
cll_model <- '
$PARAM @annotated
// ── Ibrutinib (IB) PK: depot + central, single compartment, linear elimination ──
KA_IB   : 0.50   : Ibrutinib absorption rate (h-1) [was Ka_IB]
V1_IB   : 10000  : Central volume of distribution (L) [was Vd_IB]
CL_IB   : 980    : Clearance (L/h); t1/2~7h
F_IB    : 0.25   : Oral bioavailability (fasted)

// ── Ibrutinib (IB) — BTK covalent occupancy (irreversible, PD-only: no feedback into ibrutinib own PK) ──
KINACT_IB : 0.10   : Max BTK inactivation rate (h-1, covalent pseudo-1st order) [was kinact_BTK]
KD_IB     : 1.5    : Half-maximal ibrutinib concentration for BTK inactivation (nM) [was Ki_BTK]
KDEG_IB   : 0.010  : BTK protein turnover (h-1; t1/2~2.9 d, de novo synth) [was kdeg_BTK]
RTOT_IB   : 100    : Total normalized BTK pool (%) [new named param; was the literal 100.0 in original dxdt_BTK_FREE synth term and $INIT]
EMAX_IB   : 0.70   : Max CLL proliferation/survival inhibition by BTKi [was Emax_BTK]
EC50_IB   : 50.0   : EC50 BTK occupancy for CLL PD effect (%) [was EC50_BTK]
GAMMA_IB  : 1      : Hill coefficient [new; original had none, i.e. implicit 1]

// ── Venetoclax (VEN) PK: depot + central + peripheral, linear (archetype 3) ──
KA_VEN  : 0.30   : Venetoclax absorption rate (h-1)
V1_VEN  : 250    : Central volume (L)
V2_VEN  : 500    : Peripheral volume (L)
CL_VEN  : 65     : Clearance (L/h); t1/2~26h
Q_VEN   : 10     : Inter-compartment CL (L/h)
F_VEN   : 0.50   : Bioavailability (with fat meal; fasted ~35%)

// ── Venetoclax (VEN) — BCL-2 occupancy (quasi-equilibrium relaxation, PD-only: no feedback into venetoclax own PK) ──
KD_VEN   : 0.01   : BCL-2 apparent Kd for venetoclax (nM; ~0.01 nM measured) [was Ki_BCL2]
RTOT_VEN : 100    : Total BCL-2 (normalised units, =100 at baseline) [was BCL2_tot]
KOUT_VEN : 0.10   : Rate to reach quasi-SS occupancy (h-1) [was kout_BCL2]
EMAX_VEN : 0.90   : Max CLL apoptosis induction by venetoclax [was Emax_BCL2]
EC50_VEN : 35.0   : EC50 BCL-2 occupancy (%) [was EC50_BCL2]
GAMMA_VEN: 1      : Hill coefficient [new; original had none, i.e. implicit 1]

// ── Obinutuzumab (OBI) PK/TMDD: no depot (IV), central+peripheral, CD20 target binding feeds back into OBI own central-compartment mass balance (genuine TMDD, QSS-relaxation form) ──
V1_OBI   : 3.4    : Central volume IgG (L)
V2_OBI   : 3.0    : Peripheral volume (L)
CL_OBI   : 0.020  : Linear clearance (L/h)
Q_OBI    : 0.150  : Inter-comp CL (L/h)
KD_OBI   : 0.001  : Obi-CD20 apparent Kd (mg/L; empirical TMDD) [was Kd_CD20]
KINT_OBI : 0.0003 : CD20-complex internalization rate (h-1; Type II low) [was kint_CD20]
KSYN_OBI : 0.050  : CD20 synthesis rate constant (h-1) [was ksyn_CD20]
KDEG_OBI : 0.005  : CD20 basal degradation (h-1) [was kdeg_CD20]
RTOT_OBI : 100    : Baseline CD20 level (normalised) [was CD20_0]
EMAX_OBI : 0.75   : Max CLL kill by anti-CD20 [was Emax_CD20]
EC50_OBI : 50.0   : EC50 CD20 occupancy (%) [was EC50_CD20]
GAMMA_OBI: 1      : Hill coefficient [new; original had none, i.e. implicit 1]

// ── Disease Model ─────────────────────────────────────────────────────────
kprol_CLL : 0.0030  : CLL net growth rate (h-1; doubling ~9-10 months)
Kmax_ALC  : 300.0   : Carrying capacity ALC (x1e9/L)
ALC_BASELINE : 50.0 : Baseline ALC (x1e9/L); typical symptomatic patient [renamed from ALC_0 -- that name collided with mrgsolve auto-generated <CMT>_0 init symbol for compartment ALC once $INIT was moved into $MAIN to fix a build defect; see refactor notes]
BM_0      : 70.0    : Baseline BM infiltration (%)
LN_0      : 60.0    : Baseline LN burden (relative %)

// ── BTKi redistribution (lymphocytosis) ─────────────────────────────────
kegress   : 0.012  : BM/LN -> PB egress rate under BTKi (h-1)
egress_thr: 20.0   : BTK occupancy threshold to trigger egress (%)

// ── MCL-1 resistance dynamics ────────────────────────────────────────────
kin_MCL1  : 0.008  : MCL-1 upregulation rate (h-1) under venetoclax
kout_MCL1 : 0.050  : MCL-1 normalisation rate (h-1)
MCL1_max  : 4.0    : Max MCL-1 fold-upregulation

// ── NK cell activation ────────────────────────────────────────────────────
kin_NK    : 0.005  : NK activation rate (h-1) per unit CD20 occupancy
kout_NK   : 0.020  : NK activation decay (h-1)
NK_max    : 3.0    : Max NK fold activation

// ── Dose flags (1=active, 0=off); pre-existing dead parameters -- never read in $ODE/$MAIN/$TABLE in the original, preserved as-is ──
use_IB    : 0   : Give ibrutinib (1=yes)
use_VEN   : 0   : Give venetoclax (1=yes)
use_OBI   : 0   : Give obinutuzumab (1=yes)

$CMT @annotated
GUT_IB        : Ibrutinib gut depot (mg) [was DEPOT_IB]
CENT_IB       : Ibrutinib central amount (mg)
REC_FREE_IB   : Free BTK (% normalized, 100=fully free) [was BTK_FREE]
COMPLEX_IB    : BTK-ibrutinib covalent complex (%) [was BTK_OCC]
GUT_VEN       : Venetoclax gut depot (mg) [was DEPOT_VEN]
CENT_VEN      : Venetoclax central (mg)
PERI_VEN      : Venetoclax peripheral (mg)
REC_FREE_VEN  : Free BCL-2 (normalized units) [was BCL2_FREE]
COMPLEX_VEN   : Venetoclax-BCL2 complex (normalized) [was BCL2_OCC]
CENT_OBI      : Obinutuzumab central (mg)
PERI_OBI      : Obinutuzumab peripheral (mg)
REC_FREE_OBI  : Free CD20 receptor (normalized) [was CD20_FREE]
COMPLEX_OBI   : Obi-CD20 complex (normalized) [was CD20_OCC]
ALC           : Peripheral blood CLL (x1e9/L)
BM_CLL        : Bone marrow CLL infiltration (%)
LN_CLL        : Lymph node CLL burden (%)
MCL1_ADAPT    : MCL-1 adaptive upregulation (fold)
NK_ACT        : NK cell activation (fold)

$MAIN
// Build-defect workaround (syntax-only, non-numeric): the original set initial
// values via a separate $INIT block, which mrgsolve 2.0.1 refuses to combine
// with an annotated $CMT block ("Duplicated model names"). Moved into $MAIN
// using mrgsolve <CMT>_0 idiom; every value below is copied verbatim from
// the original $INIT block. See refactor notes / UPSTREAM_ISSUES.md #45.
GUT_IB_0       = 0;
CENT_IB_0      = 0;
REC_FREE_IB_0  = RTOT_IB;
COMPLEX_IB_0   = 0;
GUT_VEN_0      = 0;
CENT_VEN_0     = 0;
PERI_VEN_0     = 0;
REC_FREE_VEN_0 = RTOT_VEN;
COMPLEX_VEN_0  = 0;
CENT_OBI_0     = 0;
PERI_OBI_0     = 0;
REC_FREE_OBI_0 = RTOT_OBI;
COMPLEX_OBI_0  = 0;
ALC_0          = 50.0;
BM_CLL_0       = 70.0;
LN_CLL_0       = 60.0;
MCL1_ADAPT_0   = 1.0;
NK_ACT_0       = 1.0;

$ODE
// ── Ibrutinib (IB) PK ────────────────────────────────────────────────────
dxdt_GUT_IB   = -KA_IB * GUT_IB;
// CENT_IB in mg; C_IB is the exposed concentration (nM)
dxdt_CENT_IB  =  KA_IB * F_IB * GUT_IB - (CL_IB / V1_IB) * CENT_IB;
double C_IB   = (CENT_IB / V1_IB) * (1000.0 / 440.5); // mg/L -> nM (MW=440.5)

// ── Ibrutinib (IB) — BTK covalent occupancy ─────────────────────────────
double k_inact_IB = KINACT_IB * C_IB / (KD_IB + C_IB);
// de novo synth restores free BTK; covalent complex decays at same kdeg
dxdt_REC_FREE_IB = KDEG_IB * RTOT_IB - KDEG_IB * REC_FREE_IB - k_inact_IB * REC_FREE_IB;
dxdt_COMPLEX_IB  = k_inact_IB * REC_FREE_IB - KDEG_IB * COMPLEX_IB;
double EFFECT_IB = EMAX_IB * pow(COMPLEX_IB, GAMMA_IB) / (pow(EC50_IB, GAMMA_IB) + pow(COMPLEX_IB, GAMMA_IB));

// ── Venetoclax (VEN) PK ──────────────────────────────────────────────────
dxdt_GUT_VEN  = -KA_VEN * GUT_VEN;
dxdt_CENT_VEN =  KA_VEN * F_VEN * GUT_VEN
                 - (CL_VEN + Q_VEN) / V1_VEN * CENT_VEN
                 + Q_VEN / V2_VEN * PERI_VEN;
dxdt_PERI_VEN =  Q_VEN / V1_VEN * CENT_VEN - Q_VEN / V2_VEN * PERI_VEN;
double C_VEN  = (CENT_VEN / V1_VEN) * (1000.0 / 868.4); // MW=868.4 g/mol

// ── Venetoclax (VEN) — BCL-2 occupancy (quasi-steady state approach) ─────
double OCC_SS_VEN = RTOT_VEN * C_VEN / (KD_VEN + C_VEN);
dxdt_REC_FREE_VEN = KOUT_VEN * (RTOT_VEN - OCC_SS_VEN - REC_FREE_VEN);
dxdt_COMPLEX_VEN  = KOUT_VEN * (OCC_SS_VEN - COMPLEX_VEN);
double OCC_FRAC_VEN = (REC_FREE_VEN + COMPLEX_VEN > 0.001) ?
                       COMPLEX_VEN / (REC_FREE_VEN + COMPLEX_VEN) * 100.0 : 0;
double EFFECT_VEN = EMAX_VEN * pow(OCC_FRAC_VEN, GAMMA_VEN) / (pow(EC50_VEN, GAMMA_VEN) + pow(OCC_FRAC_VEN, GAMMA_VEN));

// ── Obinutuzumab (OBI) PK / TMDD ─────────────────────────────────────────
double C_OBI = CENT_OBI / V1_OBI;
double KON_OBI = 0.01; // apparent kon (L/mg/h); dead -- declared, never used downstream (pre-existing in original as k_on_CD20, preserved verbatim)
double OCC_SS_OBI = RTOT_OBI * C_OBI / (KD_OBI + C_OBI);
dxdt_CENT_OBI = -(CL_OBI + Q_OBI) / V1_OBI * CENT_OBI
                 + Q_OBI / V2_OBI * PERI_OBI
                 - KINT_OBI * (OCC_SS_OBI - COMPLEX_OBI) * V1_OBI * 0.01;
dxdt_PERI_OBI =  Q_OBI / V1_OBI * CENT_OBI - Q_OBI / V2_OBI * PERI_OBI;
dxdt_REC_FREE_OBI = KSYN_OBI * RTOT_OBI - KDEG_OBI * REC_FREE_OBI
                    - KINT_OBI * (OCC_SS_OBI - COMPLEX_OBI) * 0.5;
dxdt_COMPLEX_OBI  = KINT_OBI * (OCC_SS_OBI - COMPLEX_OBI) * 0.5
                    - KDEG_OBI * COMPLEX_OBI;
double OCC_FRAC_OBI = (REC_FREE_OBI + COMPLEX_OBI > 0.001) ?
                       COMPLEX_OBI / (REC_FREE_OBI + COMPLEX_OBI) * 100.0 : 0;
double EFFECT_OBI = EMAX_OBI * pow(OCC_FRAC_OBI, GAMMA_OBI) / (pow(EC50_OBI, GAMMA_OBI) + pow(OCC_FRAC_OBI, GAMMA_OBI));

// ── Drug-effect consumption (disease-level combination; arithmetic unchanged from original) ──
double E_BTKi  = EFFECT_IB;
double E_BCL2i = EFFECT_VEN / MCL1_ADAPT;
double E_CD20  = EFFECT_OBI * NK_ACT;

// Composite kill rates per compartment
double kill_ALC = (E_BTKi * 0.50 + E_BCL2i * 0.80 + E_CD20 * 0.50) * ALC;
double kill_BM  = (E_BTKi * 0.30 + E_BCL2i * 0.90 + E_CD20 * 0.40) * BM_CLL;
double kill_LN  = (E_BTKi * 0.60 + E_BCL2i * 0.70 + E_CD20 * 0.70) * LN_CLL;

// BTKi redistribution: CLL cells egress from BM/LN to PB
double do_egress  = (COMPLEX_IB > egress_thr) ? 1.0 : 0.0;
double egress_BM  = do_egress * kegress * BM_CLL;
double egress_LN  = do_egress * kegress * LN_CLL * 0.6;

// ── Disease Compartment ODEs ─────────────────────────────────────────────
double ALC_pos = (ALC > 0) ? ALC : 0;
double ALC_growth = kprol_CLL * ALC_pos * (1.0 - ALC_pos / Kmax_ALC);
dxdt_ALC = ALC_growth - kill_ALC + egress_BM + egress_LN;
if(ALC < 0.001) dxdt_ALC = 0;

double BM_pos = (BM_CLL > 0) ? BM_CLL : 0;
dxdt_BM_CLL = kprol_CLL * 0.8 * BM_pos * (1.0 - BM_pos / 100.0) - kill_BM - egress_BM;
if(BM_CLL < 0.001) dxdt_BM_CLL = 0;

double LN_pos = (LN_CLL > 0) ? LN_CLL : 0;
dxdt_LN_CLL = kprol_CLL * 1.2 * LN_pos * (1.0 - LN_pos / 100.0) - kill_LN - egress_LN;
if(LN_CLL < 0.001) dxdt_LN_CLL = 0;

// ── MCL-1 Adaptive Resistance ─────────────────────────────────────────────
double stim_MCL1 = OCC_FRAC_VEN / 100.0;
dxdt_MCL1_ADAPT = kin_MCL1 * stim_MCL1 * (MCL1_max - MCL1_ADAPT)
                  - kout_MCL1 * (MCL1_ADAPT - 1.0);
if(MCL1_ADAPT < 1.0) dxdt_MCL1_ADAPT = 0;

// ── NK Cell Activation ────────────────────────────────────────────────────
dxdt_NK_ACT = kin_NK * (OCC_FRAC_OBI / 100.0) * (NK_max - NK_ACT)
              - kout_NK * (NK_ACT - 1.0);
if(NK_ACT < 1.0) dxdt_NK_ACT = 0;

$TABLE
// ── Concentrations ────────────────────────────────────────────────────────
double C_IB_ngmL   = (CENT_IB  / V1_IB) * 1000.0;
double C_VEN_ngmL  = (CENT_VEN / V1_VEN) * 1000.0;
double C_OBI_ugmL  = CENT_OBI / V1_OBI * 1000.0;

// ── Occupancies ───────────────────────────────────────────────────────────
double BTK_OCC_out  = COMPLEX_IB;
double BCL2_OCC_out = (REC_FREE_VEN + COMPLEX_VEN > 0.001) ?
                       COMPLEX_VEN / (REC_FREE_VEN + COMPLEX_VEN) * 100.0 : 0;
double CD20_OCC_out = (REC_FREE_OBI + COMPLEX_OBI > 0.001) ?
                       COMPLEX_OBI / (REC_FREE_OBI + COMPLEX_OBI) * 100.0 : 0;

// ── Response flags (simplified IWCLL 2018) ───────────────────────────────
double ALC_pch  = (ALC_BASELINE > 0) ? (ALC - ALC_BASELINE) / ALC_BASELINE * 100.0 : 0;
int CR_flag  = (ALC < 4.0 && BM_CLL < 30.0 && LN_CLL < 20.0) ? 1 : 0;
int PR_flag  = (ALC_pch < -50.0 && !CR_flag) ? 1 : 0;
int PD_flag  = (ALC_pch >  50.0 && ALC > ALC_BASELINE * 1.5) ? 1 : 0;
int MRD_neg  = (ALC < 0.1 && BM_CLL < 5.0) ? 1 : 0;

double BURDEN = (ALC / Kmax_ALC * 100.0 + BM_CLL + LN_CLL) / 3.0;

$CAPTURE @annotated
BTK_OCC_out  : BTK occupancy (%)
BCL2_OCC_out : BCL-2 occupancy (%)
CD20_OCC_out : CD20 occupancy (%)
C_IB_ngmL    : Ibrutinib (ng/mL)
C_VEN_ngmL   : Venetoclax (ng/mL)
C_OBI_ugmL   : Obinutuzumab (ug/mL)
BURDEN       : Composite tumor burden (%)
CR_flag      : CR achieved (1=yes)
PR_flag      : PR achieved (1=yes)
MRD_neg      : MRD-undetectable (1=yes)
ALC_pch      : ALC % change from baseline
C_IB         : Ibrutinib exposed concentration (nM) -- the single pluggable PK covariate
C_VEN        : Venetoclax exposed concentration (nM) -- the single pluggable PK covariate
C_OBI        : Obinutuzumab exposed concentration (mg/L) -- the single pluggable PK covariate
EFFECT_IB    : Ibrutinib Hill effect on CLL (0-1; pure drug term, exact rename of original E_BTKi)
EFFECT_VEN   : Venetoclax Hill effect on CLL (0-1; pure drug term before MCL1-resistance division, exact rename of original E_BCL2i numerator)
EFFECT_OBI   : Obinutuzumab Hill effect on CLL (0-1; pure drug term before NK-activation multiplier, exact rename of original E_CD20 own ratio)
'

# ── Compile Model ─────────────────────────────────────────────────────────────
mod <- mcode("CLL_QSP_refactored", cll_model)

# ── Helper: build dosing events (renamed compartments only; schedule unchanged) ──
dose_events <- function(scenario = 1, end_days = 730) {
  ev_list <- list()

  # Scenario 1: Ibrutinib monotherapy 420 mg QD
  if (scenario %in% c(1)) {
    ev_list[["IB"]] <- ev(amt = 420, cmt = "GUT_IB",
                          time = 0, ii = 24, addl = end_days - 1)
  }

  # Scenario 2: Venetoclax monotherapy (5-week ramp → 400 mg QD)
  if (scenario %in% c(2, 4, 5, 6)) {
    ramp <- data.frame(
      amt  = c(20, 50, 100, 200, 400),
      addl = c(6,   6,   6,   6,   end_days * 7 - 28) / 1,
      time = c(0, 168, 336, 504, 672)  # hours: wk0,1,2,3,4+
    )
    for (i in seq_len(nrow(ramp))) {
      ev_list[[paste0("VEN_ramp", i)]] <- ev(
        amt  = ramp$amt[i], cmt = "GUT_VEN",
        time = ramp$time[i], ii = 24,
        addl = ceiling((ramp$time[min(i + 1, 5)] - ramp$time[i]) / 24) - 1
      )
    }
    ev_list[["VEN_main"]] <- ev(amt = 400, cmt = "GUT_VEN",
                                time = 672, ii = 24,
                                addl = end_days * 24 - 1)
  }

  # Scenario 3: Obinutuzumab (cycles 1-6 q28d; cycle 1 split: D1=100mg, D2=900mg, D15=1000mg)
  if (scenario %in% c(3, 4, 6)) {
    obi_times <- c(0, 24, 336, 672, 1344, 2016, 2688, 3360) # h: D1,D2,D15,C2-C6
    obi_amts  <- c(100, 900, 1000, 1000, 1000, 1000, 1000, 1000)
    for (i in seq_along(obi_times)) {
      ev_list[[paste0("OBI_", i)]] <- ev(amt = obi_amts[i], cmt = "CENT_OBI",
                                         time = obi_times[i])
    }
  }

  # Scenario 5: Ibrutinib + Venetoclax (MRD-guided fixed-duration)
  if (scenario == 5) {
    ev_list[["IB5"]] <- ev(amt = 420, cmt = "GUT_IB",
                           time = 0, ii = 24, addl = end_days - 1)
  }

  # Scenario 6: Ven + Obi (CLL14 regimen: Obi C1-C6 + Ven 12 months)
  if (scenario == 6) {
    ev_list[["IB6"]] <- NULL  # no ibrutinib -- pre-existing in the original: scenario 6
    # is labeled "Triplet IB+VEN+OBI" below but never actually doses ibrutinib;
    # its dosing is byte-identical to scenario 4. Observed, not fixed -- see
    # refactor notes (this is an R-driver-script inconsistency, not a DSL defect).
  }

  do.call(c, ev_list)
}

# ── Scenario Definitions ──────────────────────────────────────────────────────
scenarios <- list(
  "1_ibrutinib_mono"      = list(use_IB=1, use_VEN=0, use_OBI=0, label="Ibrutinib 420mg QD"),
  "2_venetoclax_mono"     = list(use_IB=0, use_VEN=1, use_OBI=0, label="Venetoclax 400mg QD"),
  "3_obinutuzumab_mono"   = list(use_IB=0, use_VEN=0, use_OBI=1, label="Obinutuzumab x6 cycles"),
  "4_ven_obi_cll14"       = list(use_IB=0, use_VEN=1, use_OBI=1, label="VEN+OBI (CLL14)"),
  "5_ib_ven_combo"        = list(use_IB=1, use_VEN=1, use_OBI=0, label="Ibrutinib+Venetoclax"),
  "6_triplet"             = list(use_IB=1, use_VEN=1, use_OBI=1, label="Triplet IB+VEN+OBI")
)

END_DAYS <- 730  # 2 years simulation

# ── Run all scenarios ─────────────────────────────────────────────────────────
run_scenario <- function(scen_id, params, end_days = END_DAYS) {
  e <- dose_events(scen_id, end_days)
  p <- list(use_IB  = params$use_IB,
            use_VEN = params$use_VEN,
            use_OBI = params$use_OBI,
            ALC_BASELINE = 50, BM_0 = 70, LN_0 = 60)
  out <- mrgsim(mod, ev = e, param = p,
                start = 0, end = end_days * 24, delta = 12) %>%
    as_tibble() %>%
    mutate(scenario = params$label,
           time_days = time / 24)
  out
}

# Map scenario number to param list
scen_nums <- c(1, 2, 3, 4, 5, 6)
names(scen_nums) <- names(scenarios)

results <- bind_rows(lapply(seq_along(scenarios), function(i) {
  tryCatch(
    run_scenario(scen_nums[i], scenarios[[i]]),
    error = function(e) { message("Scenario ", i, " error: ", e$message); NULL }
  )
}))

# ── Plot 1: ALC over time ─────────────────────────────────────────────────────
p1 <- ggplot(results, aes(time_days, ALC, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 4, linetype = "dashed", color = "grey50") +
  annotate("text", x = 5, y = 5, label = "IWCLL CR threshold (4×10⁹/L)",
           hjust = 0, size = 3.5, color = "grey50") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "CLL — Absolute Lymphocyte Count (ALC) by Treatment",
       x = "Time (days)", y = "ALC (×10⁹/L)", color = "Scenario") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# ── Plot 2: BTK and BCL-2 occupancy ──────────────────────────────────────────
p2 <- results %>%
  filter(scenario %in% c("Ibrutinib 420mg QD", "VEN+OBI (CLL14)",
                          "Ibrutinib+Venetoclax", "Triplet IB+VEN+OBI")) %>%
  select(time_days, scenario, BTK_OCC_out, BCL2_OCC_out) %>%
  pivot_longer(c(BTK_OCC_out, BCL2_OCC_out),
               names_to = "target", values_to = "occupancy") %>%
  mutate(target = recode(target,
    BTK_OCC_out  = "BTK Occupancy (%)",
    BCL2_OCC_out = "BCL-2 Occupancy (%)")) %>%
  ggplot(aes(time_days, occupancy, color = scenario, linetype = target)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 95, linetype = "dotted", color = "navy") +
  annotate("text", x = 5, y = 96, label = "95% BTK target", hjust = 0, size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Target Occupancy: BTK vs BCL-2",
       x = "Time (days)", y = "Occupancy (%)",
       color = "Scenario", linetype = "Target") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

# ── Plot 3: Tumor burden composite ───────────────────────────────────────────
p3 <- ggplot(results, aes(time_days, BURDEN, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Composite Tumor Burden (ALC + BM + LN average)",
       x = "Time (days)", y = "Composite burden (%)", color = "Scenario") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

# ── Plot 4: MCL-1 resistance & NK activation ──────────────────────────────────
p4 <- results %>%
  select(time_days, scenario, MCL1_ADAPT, NK_ACT) %>%
  pivot_longer(c(MCL1_ADAPT, NK_ACT)) %>%
  ggplot(aes(time_days, value, color = scenario, linetype = name)) +
  geom_line(linewidth = 0.9) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Adaptive Resistance (MCL-1) and NK Activation",
       x = "Time (days)", y = "Fold change",
       color = "Scenario", linetype = "Variable") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")

# ── Summary response table ────────────────────────────────────────────────────
response_table <- results %>%
  group_by(scenario) %>%
  summarise(
    ALC_nadir      = round(min(ALC, na.rm = TRUE), 2),
    ALC_nadir_day  = time_days[which.min(ALC)],
    CR_achieved    = ifelse(any(CR_flag == 1), "Yes", "No"),
    MRD_neg_pct    = round(mean(MRD_neg, na.rm = TRUE) * 100, 1),
    BM_final       = round(last(BM_CLL), 1),
    LN_final       = round(last(LN_CLL), 1),
    .groups = "drop"
  )

cat("\n=== CLL Treatment Response Summary (2-year simulation) ===\n")
print(response_table)

# ── Print plots ───────────────────────────────────────────────────────────────
print(p1)
print(p2)
print(p3)
print(p4)

# ── Clinical trial calibration notes ─────────────────────────────────────────
cat("
=== Calibration Reference Points ===
RESONATE-2 (ibrutinib vs chlorambucil, Burger 2015 NEJM):
  Ibrutinib: ORR 86%, 2yr PFS 74%; typical ALC nadir ~6-12 months
  Model target: ALC decline to <10 x1e9/L by 6-12 months with BTKi

CLL14 (venetoclax+obinutuzumab, Fischer 2019 NEJM):
  VEN+OBI: 2yr PFS 88.2% vs 64.1% (chlorambucil+obi)
  MRD-negative rate: 76.4% (blood), 57.0% (marrow)
  Model target: ALC <4 and BM <30 by 12 months in most VEN+OBI patients

MURANO (venetoclax+rituximab, Seymour 2018 NEJM):
  VEN+R: 2yr PFS 84.9% vs 36.3% (bendamustine+rituximab)
  uMRD blood: 83% at EOT
  Model target: Deep response with combination VEN regimens

SEQUOIA (zanubrutinib vs chlorambucil, Shadman 2023 NEJM):
  zanubrutinib: 2yr PFS 85.5%, ORR 94.6%
  Less AF vs ibrutinib (2.5% vs 10.1% in ALPINE trial)

Parameter calibration: kprol_CLL=0.003 h-1 yields LDT ~9.7 months
  consistent with intermediate-risk newly diagnosed CLL.
  Emax and EC50 values calibrated to achieve ~86% ORR for ibrutinib
  and ~76% MRD-neg for VEN+OBI at 12 months.
")
