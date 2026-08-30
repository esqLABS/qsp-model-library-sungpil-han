# =============================================================================
# Rheumatoid Arthritis (RA) QSP Model — mrgsolve
# 류마티스 관절염 정량적 시스템 약리학 모델
# =============================================================================
# Based on: Yuraszeck et al. (2010), Rullmann et al. (2005),
#           Miossec & Kolls (2012), van Vollenhoven et al. (2014)
# Drug PK: Certara PopPK reports, FDA clinical pharmacology reviews
# =============================================================================
#
# PK/PD REFACTOR (see FORK_WORKFLOW_GUIDE.md Part 2 and
# rheumatoid-arthritis/ra_refactor_notes.md for full detail):
#   Only tocilizumab's PK + IL-6R binding block, and its interface into the
#   disease equations, was touched. Adalimumab, methotrexate, and
#   baricitinib are byte-for-byte identical to ra_mrgsolve_model.R.
#   - SC_TCZ/IL6R_FREE/DRUG_RL/OCC_TCZ/F_SC/KON/KOFF/KSYN_R/KDEG_R/RTOT0
#     renamed to convention: GUT_TCZ/REC_FREE_TCZ/COMPLEX_TCZ/EFFECT_TCZ/
#     F_TCZ/KON_TCZ/KOFF_TCZ/KSYN_TCZ/KDEG_TCZ/RTOT_TCZ.
#   - This is archetype 4 (TMDD); the receptor-binding ODEs are unchanged
#     pharmacology, renamed only.
#   - EFFECT_TCZ is a RENAME of OCC_TCZ, not a fit: OCC_TCZ = DRUG_RL/RTOT is
#     already exactly Hill-shaped (gamma=1) in C_TCZ at quasi-steady-state,
#     confirmed analytically and by nls() (R^2=1.0000, see notes). EFFECT_TCZ
#     still computes from the exact ODE state ratio, so this is a pure
#     rename with zero added approximation.
#   - $CAPTURE also had CRP/INFLAM/VDH/CART/STAT3_P/FLS_ACT (and OCC_TCZ,
#     now EFFECT_TCZ) de-duplicated against $CMT: mrgsolve 2.0.1 refuses to
#     build a model whose $CAPTURE repeats a compartment name, and the
#     ORIGINAL file has this defect too (confirmed independently; logged to
#     translations/UPSTREAM_ISSUES.md). Compartment values still appear in
#     every output automatically, so this changes nothing about what is
#     reported -- it only removes an illegal duplicate declaration required
#     to make either file buildable under the currently installed mrgsolve.
#
# Verified against the original (rerun via the qspserver mrgsolve_api
# service) on the original file's own TCZ IV, TCZ SC, and TCZ+MTX combo
# dosing scenarios: every $CAPTURE'd output matched the original EXACTLY
# (max abs diff 0.0) across the full 1-year daily time grid.
# =============================================================================

library(mrgsolve)
library(tidyverse)
library(patchwork)

# Model code
code <- '
$PROB
RA QSP Model: IL-6 Receptor Occupancy, CRP, DAS28, Bone Erosion
Drugs: Tocilizumab (IV/SC), Methotrexate, Adalimumab, Baricitinib

$PARAM @annotated
// === TOCILIZUMAB (TCZ) PK — archetype 4, TMDD (refactored: renamed to convention) ===
CL_TCZ   : 0.52   : L/h   Tocilizumab linear clearance
V1_TCZ   : 3.72   : L     Central volume
V2_TCZ   : 2.35   : L     Peripheral volume
Q_TCZ    : 0.56   : L/h   Intercompartmental CL
F_TCZ    : 0.80   :       SC bioavailability (was F_SC)
KA_TCZ   : 0.108  : 1/h   SC absorption rate constant
KON_TCZ  : 0.0033 : L/nmol/h  Drug-IL6R association (was KON)
KOFF_TCZ : 0.019  : 1/h   Drug-IL6R dissociation (was KOFF)
KSYN_TCZ : 0.1    : nmol/L/h  IL-6R synthesis (was KSYN_R)
KDEG_TCZ : 0.007  : 1/h   IL-6R degradation (was KDEG_R)
// RTOT_TCZ (was RTOT0): unused in the ODE itself in the original -- see refactor notes
RTOT_TCZ : 14.3   : nmol/L    Baseline membrane IL-6R
// Hill-interface parameters for EFFECT_TCZ. Not a fit: the original
// occupancy ratio COMPLEX_TCZ/(REC_FREE_TCZ+COMPLEX_TCZ) is an EXACT Hill
// function of C_TCZ at quasi-steady-state (binding equilibrates in <1h,
// fast vs. the PK/dosing timescale): Occ(C) = C/(C+Kd), Kd=(KOFF+KDEG)/KON.
// Confirmed analytically and by nls() sampling C_TCZ over the scenarios own
// concentration range (~1e-3-2000 nmol/L, 18 log-spaced points held to
// steady state): Emax=1.000000, EC50=7.878788, gamma=1.000000, R^2=1.0000
// (max abs residual ~3e-10). EFFECT_TCZ below still uses the exact ODE
// state ratio (zero extra approximation); these are exposed only so a
// future external-covariate substitution has the closed-form ready to use.
EMAX_TCZ : 1      :       TCZ occupancy-effect Emax (derived, see above)
EC50_TCZ : 7.878788 : nmol/L  TCZ occupancy-effect EC50 = (KOFF_TCZ+KDEG_TCZ)/KON_TCZ
GAMMA_TCZ: 1      :       TCZ occupancy-effect Hill coefficient (derived)

// === ADALIMUMAB PK ===
CL_ADA   : 0.0114 : L/h   Adalimumab clearance
V1_ADA   : 2.84   : L     Central volume
KA_ADA   : 0.0073 : 1/h   SC absorption rate

// === METHOTREXATE ===
KA_MTX   : 0.3    : 1/h   MTX oral absorption rate
CL_MTX   : 5.4    : L/h   MTX clearance
V_MTX    : 28     : L     MTX volume of distribution
KPOLY    : 0.02   : 1/h   Polyglutamate formation rate
KDEPOLY  : 0.005  : 1/h   Polyglutamate hydrolysis

// === BARICITINIB PK ===
KA_BARI  : 0.69   : 1/h   Baricitinib oral absorption
CL_BARI  : 12.1   : L/h   Baricitinib clearance
V_BARI   : 76     : L     Baricitinib volume
IC50_JAK : 5.9    : ng/mL JAK1/2 IC50

// === IL-6 & SIGNALING DYNAMICS ===
KOUT_IL6 : 0.0289 : 1/h   IL-6 elimination
KIN_IL6  : 0.0722 : nmol/L/h  Basal IL-6 synthesis (RA state)
KOUT_STAT3: 0.1   : 1/h   STAT3 dephosphorylation
IL6_EC50 : 1.5    : nmol/L    IL-6 EC50 for STAT3 activation
STAT3MAX : 100    :       STAT3 phosphorylation Emax

// === CRP DYNAMICS ===
KOUT_CRP : 0.0289 : 1/h   CRP elimination (t1/2 ~24h)
KIN_CRP0 : 2.89   : mg/L/h    Baseline CRP synthesis rate
CRP_EC50 : 0.5    :       STAT3 EC50 for CRP production
CRP_BASE : 15     : mg/L  RA baseline CRP

// === INFLAMMATION / SYNOVIUM ===
K_INFLAM_IN  : 0.01 :     Inflammation stimulus rate
K_INFLAM_OUT : 0.001 :    Inflammation natural decay
INFLAM_BASE  : 10  :      Baseline inflammation score (RA)
K_STAT3_INFLAM: 0.2 :     STAT3 drive on inflammation
K_TNF_INFLAM : 0.15 :     TNF drive on inflammation
INFLAM_EMAX  : 20  :      Emax for drug effect on inflammation

// === FLS ACTIVATION ===
KFLS_IN  : 0.001  : 1/h   FLS activation rate
KFLS_OUT : 0.0005 : 1/h   FLS resolution rate

// === BONE EROSION (Sharp/vdH score) ===
KVDH_PROG: 0.00015 :      Joint damage progression rate per unit inflammation
KVDH_INHIB: 0.7   :       Drug effect on progression (max inhibition)
VDH_BASE : 10     :       Baseline Sharp/vdH score at treatment start

// === CARTILAGE (JSN) ===
K_CART_DESTR: 0.00008 :   Cartilage loss rate per inflammation unit
CART_INIT: 100    :       Initial cartilage integrity (%)

// === DISEASE ACTIVITY ===
DAS28_COEFF: 1.0  :       DAS28 scaling

$CMT @annotated
// Tocilizumab PK (archetype 4, TMDD; renamed to convention)
GUT_TCZ      : SC depot tocilizumab [nmol] (was SC_TCZ)
CENT_TCZ     : Central tocilizumab [nmol]
PERI_TCZ     : Peripheral tocilizumab [nmol]
REC_FREE_TCZ : Free membrane IL-6R [nmol/L] (was IL6R_FREE)
COMPLEX_TCZ  : TCZ-IL6R complex [nmol/L] (was DRUG_RL)

// Adalimumab PK
SC_ADA    : SC depot adalimumab [nmol]
CENT_ADA  : Central adalimumab [nmol]

// MTX
GUT_MTX   : MTX gut lumen [mg]
CENT_MTX  : Central MTX [mg]
POLY_MTX  : MTX polyglutamates [relative units]

// Baricitinib
GUT_BARI  : Baricitinib gut [mg]
CENT_BARI : Central baricitinib [mg]

// Signaling/PD
IL6       : Free IL-6 [nmol/L]
STAT3_P   : Phospho-STAT3 [relative units]
CRP       : C-reactive protein [mg/L]
INFLAM    : Synovial inflammation index [0-20]
FLS_ACT   : FLS activation level

// Structural damage
VDH       : Sharp/vdH total score
CART      : Cartilage integrity [%]

$MAIN
double C_TCZ  = CENT_TCZ / V1_TCZ;   // nmol/L
double C_ADA  = CENT_ADA / V1_ADA;
double C_MTX  = CENT_MTX / V_MTX;
double C_BARI = CENT_BARI / V_BARI;

// Tocilizumab Hill interface: fractional TCZ-IL6R occupancy.
// This is the compound effect on the disease system in the original file
// (used directly, undampened, as the fractional block on IL-6 signaling) --
// renamed (COMPLEX_TCZ / current total receptor), not refit; see the
// $PARAM block above and the refactor notes for why this qualifies as
// already Hill-shaped. RTOT_now_TCZ is the dynamic free+complex total
// (distinct from the RTOT_TCZ parameter, the baseline/steady-state value).
double RTOT_now_TCZ = REC_FREE_TCZ + COMPLEX_TCZ;
double EFFECT_TCZ   = COMPLEX_TCZ / (RTOT_now_TCZ + 1e-6);

// JAK inhibition by baricitinib
double JAK_INH = C_BARI / (C_BARI + IC50_JAK);

// IL-6 signal blocked by TCZ (receptor occupancy) + JAK inhibition
double IL6_sig = (IL6 / (IL6 + IL6_EC50)) * (1 - EFFECT_TCZ) * (1 - 0.8*JAK_INH);

// MTX anti-inflammatory effect
double MTX_EFF_INFLAM = 0.35 * POLY_MTX / (POLY_MTX + 0.5);

// TNF inhibition by adalimumab (simplified: free drug = serum concentration)
double TNF_INH_ADA = C_ADA / (C_ADA + 0.035);  // IC50 ~35 pmol/L -> using nmol/L scale

// Combined anti-inflammatory effect (only the TCZ term redirected to EFFECT_TCZ;
// ADA/MTX/BARI terms and the combination structure are untouched)
double DRUG_ANTI_INFLAM = 1 - (1-EFFECT_TCZ)*(1-TNF_INH_ADA)*(1-MTX_EFF_INFLAM)*(1-0.6*JAK_INH);

$ODE
// --- Tocilizumab PK (archetype 4, TMDD) ---
dxdt_GUT_TCZ  = -KA_TCZ * GUT_TCZ * F_TCZ;
dxdt_CENT_TCZ = KA_TCZ*GUT_TCZ*F_TCZ - (CL_TCZ/V1_TCZ)*CENT_TCZ
                - (Q_TCZ/V1_TCZ)*CENT_TCZ + (Q_TCZ/V2_TCZ)*PERI_TCZ
                - V1_TCZ*(KON_TCZ * C_TCZ * REC_FREE_TCZ - KOFF_TCZ * COMPLEX_TCZ);
dxdt_PERI_TCZ = (Q_TCZ/V1_TCZ)*CENT_TCZ - (Q_TCZ/V2_TCZ)*PERI_TCZ;

// IL-6R target binding (TMDD)
dxdt_REC_FREE_TCZ = KSYN_TCZ - KDEG_TCZ*REC_FREE_TCZ - KON_TCZ*C_TCZ*REC_FREE_TCZ + KOFF_TCZ*COMPLEX_TCZ;
dxdt_COMPLEX_TCZ  = KON_TCZ*C_TCZ*REC_FREE_TCZ - KOFF_TCZ*COMPLEX_TCZ - KDEG_TCZ*COMPLEX_TCZ;

// --- Adalimumab PK ---
dxdt_SC_ADA   = -KA_ADA * SC_ADA;
dxdt_CENT_ADA = KA_ADA*SC_ADA - CL_ADA*C_ADA;

// --- Methotrexate PK ---
dxdt_GUT_MTX  = -KA_MTX * GUT_MTX;
dxdt_CENT_MTX = KA_MTX*GUT_MTX - CL_MTX*C_MTX;
dxdt_POLY_MTX = KPOLY*C_MTX - KDEPOLY*POLY_MTX;

// --- Baricitinib PK ---
dxdt_GUT_BARI  = -KA_BARI * GUT_BARI;
dxdt_CENT_BARI = KA_BARI*GUT_BARI - (CL_BARI/V_BARI)*CENT_BARI;

// --- IL-6 dynamics ---
double IL6_prod = KIN_IL6 * (1 + 0.5*INFLAM/INFLAM_BASE);
dxdt_IL6 = IL6_prod - KOUT_IL6*IL6;

// --- STAT3 phosphorylation ---
double STAT3_drive = STAT3MAX * IL6_sig;
dxdt_STAT3_P = STAT3_drive - KOUT_STAT3*STAT3_P;

// --- CRP dynamics (driven by STAT3) ---
double CRP_drive = KIN_CRP0 * STAT3_P / (STAT3_P + CRP_EC50*STAT3MAX);
dxdt_CRP = CRP_drive - KOUT_CRP*CRP;

// --- Inflammation index (synovial) ---
double INFLAM_drive = K_INFLAM_IN * (K_STAT3_INFLAM*STAT3_P/STAT3MAX
                                    + K_TNF_INFLAM*(1-TNF_INH_ADA))
                      * INFLAM_BASE;
double INFLAM_decay = (K_INFLAM_OUT + K_INFLAM_IN*DRUG_ANTI_INFLAM)*INFLAM;
dxdt_INFLAM = INFLAM_drive - INFLAM_decay;

// --- FLS activation ---
dxdt_FLS_ACT = KFLS_IN*INFLAM - KFLS_OUT*FLS_ACT*(1 + 2*DRUG_ANTI_INFLAM);

// --- Bone erosion (Sharp/vdH) ---
double PROG_RATE = KVDH_PROG * INFLAM * (1 - KVDH_INHIB*DRUG_ANTI_INFLAM);
dxdt_VDH = PROG_RATE;

// --- Cartilage loss ---
dxdt_CART = -K_CART_DESTR * INFLAM * (1 - DRUG_ANTI_INFLAM) * (CART/100.0);

$TABLE
// DAS28-CRP (simplified model)
// DAS28 = 0.56*sqrt(TJC) + 0.28*sqrt(SJC) + 0.36*ln(CRP+1) + 0.014*GH
// Proxy: TJC~14*(INFLAM/20), SJC~10*(INFLAM/20), GH proportional to VAS
double TJC28_est = 14.0 * (INFLAM / 20.0);
double SJC28_est = 10.0 * (INFLAM / 20.0);
double GH_est = 50.0 * (INFLAM / 20.0);
double DAS28_est = 0.56*sqrt(TJC28_est+0.01) + 0.28*sqrt(SJC28_est+0.01)
                 + 0.36*log(CRP+1) + 0.014*GH_est;

// Clinical response criteria
double ACR20 = (DAS28_est < 5.1) ? 1 : 0;
double ACR50 = (DAS28_est < 3.8) ? 1 : 0;
double ACR70 = (DAS28_est < 3.2) ? 1 : 0;
double REMISSION = (DAS28_est < 2.6) ? 1 : 0;
double LDA = (DAS28_est < 3.2) ? 1 : 0;

// HAQ-DI proxy (0-3 scale)
double HAQ_DI_est = 1.5 * (INFLAM/20.0) + 0.02*VDH;

// Receptor occupancy %
double RO_pct = 100 * EFFECT_TCZ;

$CAPTURE DAS28_est HAQ_DI_est RO_pct
         C_TCZ C_ADA C_MTX C_BARI ACR20 ACR50 ACR70 REMISSION
         TNF_INH_ADA EFFECT_TCZ
'

# Build model
mod <- mcode("RA_QSP", code)

# =============================================================================
# Initial conditions (RA patient at baseline)
# =============================================================================
init_RA <- init(mod,
  REC_FREE_TCZ = 14.3,
  IL6       = 2.5,
  STAT3_P   = 50,
  CRP       = 15,
  INFLAM    = 10,
  FLS_ACT   = 5,
  VDH       = 10,
  CART      = 85
)

# =============================================================================
# Simulation Scenarios
# =============================================================================

# Scenario 1: No Treatment (disease progression)
e_none <- ev(time=0, amt=0, cmt=1)
out_none <- mod %>% init_RA %>%
  mrgsim(ev=e_none, end=8760, delta=24) %>%  # 1 year
  as_tibble() %>% mutate(scenario="1. No Treatment")

# Scenario 2: MTX monotherapy 20mg/week oral
e_mtx <- ev(time=seq(0, 8736, by=168), amt=20, cmt="GUT_MTX")
out_mtx <- mod %>% init_RA %>%
  mrgsim(ev=e_mtx, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="2. MTX 20mg/wk")

# Scenario 3: Tocilizumab IV 8 mg/kg q4w (assume 70kg -> 560mg -> ~3733 nmol)
# MW TCZ ~148 kDa, 560mg/148000 g/mol * 1e6 nmol/mol = 3784 nmol
e_tcz_iv <- ev(time=seq(0, 8736, by=672), amt=3784, cmt="CENT_TCZ")
out_tcz_iv <- mod %>% init_RA %>%
  mrgsim(ev=e_tcz_iv, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="3. TCZ IV 8mg/kg q4w")

# Scenario 4: TCZ SC 162mg q2w (162mg / 148kDa = 1095 nmol)
e_tcz_sc <- ev(time=seq(0, 8736, by=336), amt=1095, cmt="GUT_TCZ")
out_tcz_sc <- mod %>% init_RA %>%
  mrgsim(ev=e_tcz_sc, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="4. TCZ SC 162mg q2w")

# Scenario 5: TCZ IV + MTX combination
e_combo2 <- ev(time=seq(0,8736,672), amt=3784, cmt="CENT_TCZ") +
            ev(time=seq(0,8736,168), amt=20, cmt="GUT_MTX")
out_combo <- mod %>% init_RA %>%
  mrgsim(ev=e_combo2, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="5. TCZ IV + MTX")

# Scenario 6: Baricitinib 4mg QD (oral)
e_bari <- ev(time=seq(0, 8759, by=24), amt=4, cmt="GUT_BARI")
out_bari <- mod %>% init_RA %>%
  mrgsim(ev=e_bari, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="6. Baricitinib 4mg QD")

# Scenario 7: Adalimumab 40mg SC q2w (40mg / 148kDa = 270 nmol, assuming same MW scale)
e_ada <- ev(time=seq(0, 8736, by=336), amt=270, cmt="SC_ADA")
out_ada <- mod %>% init_RA %>%
  mrgsim(ev=e_ada, end=8760, delta=24) %>%
  as_tibble() %>% mutate(scenario="7. Adalimumab 40mg q2w")

# Combine all scenarios
all_sims <- bind_rows(out_none, out_mtx, out_tcz_iv, out_tcz_sc,
                       out_combo, out_bari, out_ada)

# Convert time to weeks
all_sims <- all_sims %>% mutate(week = time/168)

# =============================================================================
# Visualization
# =============================================================================
scenario_colors <- c(
  "1. No Treatment"       = "#E74C3C",
  "2. MTX 20mg/wk"        = "#E67E22",
  "3. TCZ IV 8mg/kg q4w"  = "#3498DB",
  "4. TCZ SC 162mg q2w"   = "#85C1E9",
  "5. TCZ IV + MTX"       = "#1A5276",
  "6. Baricitinib 4mg QD" = "#8E44AD",
  "7. Adalimumab 40mg q2w"= "#27AE60"
)

p1 <- ggplot(all_sims, aes(week, DAS28_est, color=scenario)) +
  geom_line(linewidth=1) +
  geom_hline(yintercept=2.6, linetype="dashed", color="darkgreen", alpha=0.7) +
  geom_hline(yintercept=3.2, linetype="dashed", color="gold3", alpha=0.7) +
  scale_color_manual(values=scenario_colors) +
  labs(title="DAS28-CRP over Time", x="Week", y="DAS28-CRP", color="") +
  theme_minimal(base_size=12) +
  annotate("text", x=max(all_sims$week)*0.9, y=2.4, label="Remission", color="darkgreen", size=3) +
  annotate("text", x=max(all_sims$week)*0.9, y=3.0, label="LDA", color="gold3", size=3)

p2 <- ggplot(all_sims, aes(week, CRP, color=scenario)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=scenario_colors) +
  labs(title="CRP (mg/L) over Time", x="Week", y="CRP (mg/L)", color="") +
  theme_minimal(base_size=12)

p3 <- ggplot(all_sims, aes(week, VDH, color=scenario)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Sharp/vdH Score (Bone Erosion)", x="Week", y="Sharp/vdH Score", color="") +
  theme_minimal(base_size=12)

p4 <- ggplot(all_sims, aes(week, INFLAM, color=scenario)) +
  geom_line(linewidth=1) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Synovial Inflammation Index", x="Week", y="Inflammation Score", color="") +
  theme_minimal(base_size=12)

combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Rheumatoid Arthritis QSP Model — Treatment Scenario Comparison",
    subtitle = "mrgsolve ODE simulation | 52-week follow-up",
    theme = theme(plot.title = element_text(face="bold", size=14))
  ) +
  plot_layout(guides="collect") &
  theme(legend.position="bottom")

# print(combined_plot)

# =============================================================================
# 52-week summary table
# =============================================================================
week52 <- all_sims %>%
  filter(abs(week - 52) < 0.5) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, DAS28_est, CRP, INFLAM, VDH, CART, HAQ_DI_est, RO_pct) %>%
  rename(
    Scenario    = scenario,
    DAS28_52wk  = DAS28_est,
    CRP_52wk    = CRP,
    Inflamm     = INFLAM,
    SharpvdH    = VDH,
    Cartilage   = CART,
    HAQ_DI      = HAQ_DI_est,
    TCZ_RO_pct  = RO_pct
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

print(week52)

# =============================================================================
# Dose-response analysis for TCZ IV
# =============================================================================
doses_TCZ_mg <- c(2, 4, 6, 8, 10)  # mg/kg (70kg patient)
doses_TCZ_nmol <- doses_TCZ_mg * 70 / 148 * 1000  # nmol

dose_resp <- map_dfr(seq_along(doses_TCZ_nmol), function(i) {
  e_dr <- ev(time=seq(0,8736,672), amt=doses_TCZ_nmol[i], cmt="CENT_TCZ")
  mod %>% init_RA %>%
    mrgsim(ev=e_dr, end=8760, delta=24) %>%
    as_tibble() %>%
    filter(abs(time/168 - 24) < 0.5) %>%
    slice(1) %>%
    mutate(dose_mgkg = doses_TCZ_mg[i], dose_nmol = doses_TCZ_nmol[i])
})

cat("\n=== TCZ Dose-Response at Week 24 ===\n")
print(dose_resp %>% select(dose_mgkg, DAS28_est, CRP, EFFECT_TCZ, REMISSION))
