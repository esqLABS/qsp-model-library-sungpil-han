# =====================================================================
# X-Linked Hypophosphatemia (XLH) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : PHEX (Xp22.11) loss-of-function → osteocyte FGF23 overproduction
#            → FGF23-FGFR1c/alphaKlotho renal signaling → NPT2a/NPT2c
#            internalization → renal phosphate wasting (TmP/GFR down) +
#            suppressed 1-alpha-hydroxylase (low/inappropriately-normal
#            calcitriol) → chronic hypophosphatemia → defective growth
#            plate/osteoid mineralization (pediatric rickets / adult
#            osteomalacia), enthesopathy, myopathy, dental abscesses.
#   PK/PD  : Burosumab (fully human anti-FGF23 IgG1 mAb, SC Q2W pediatric /
#            Q4W adult) directly neutralizes circulating FGF23, restoring
#            NPT2a/c-mediated renal Pi reabsorption AND de-repressing
#            1-alpha-hydroxylase. Conventional therapy (oral phosphate
#            salts + active vitamin D / calcitriol-alfacalcidol) bypasses
#            the FGF23 axis entirely via direct substrate replacement,
#            producing transient serum Pi spikes and hypercalciuria that
#            carry nephrocalcinosis / secondary-tertiary hyperparathyroidism
#            risk with chronic use.
#   Outputs: Rickets Severity Score (RSS), annualized growth velocity (AGV),
#            height Z-score, serum phosphate, TmP/GFR, serum 1,25(OH)2D,
#            PTH, bone-specific ALP (BSAP), 6-minute walk test (6MWT, adult),
#            WOMAC pain score, nephrocalcinosis risk accumulator, urine
#            calcium/creatinine.
#   References (calibration): Carpenter et al. NEJM 2018 (PMID 29791829;
#            phase 2 pediatric burosumab dose-finding 0.3-2.0 mg/kg Q2W/Q4W),
#            Imel et al. Lancet 2019 (PMID 31104833; CL303 phase 3 RCT,
#            burosumab vs continued conventional therapy, RSS/growth
#            superiority), Whyte et al. Lancet Diabetes Endocrinol 2019
#            (PMID 31104830; phase 2, ages 1-4 yr), Insogna et al. JBMR
#            2018 (PMID 29947083; AXLES1 adult phase 3 placebo-controlled
#            RCT), Carpenter et al. JCI 2014 (PMID 24569459; KRN23/burosumab
#            first-in-human single-dose PK/PD), Carpenter et al. JBMR 2011
#            (PMID 21538511; clinician's guide / conventional-therapy
#            natural history), Haffner et al. Nat Rev Nephrol 2019
#            (PMID 31068690; international consensus diagnosis/management
#            guideline), Skrinar et al. J Endocr Soc 2019 (PMID 31259293;
#            lifelong disease-burden survey).
#
# =====================================================================
# PK/PD REFACTOR (see FORK_WORKFLOW_GUIDE.md Part 2 and
# xlh_refactor_notes.md for full detail):
#   All three compounds in this file were touched (per
#   driver-patches/data/compound_perturbation_census.md, all three rows
#   classified "Redirect concentration (clean single site)"). There is no
#   disease-model code outside these three compounds' own PK/PD blocks, so
#   nothing else in the file changed except the two build-compat fixes below.
#
#   - Burosumab: archetype 3 (depot + central + peripheral, linear).
#     BURO_DEPOT/BURO_CENT/BURO_PERIPH/BURO_CP renamed to
#     GUT_BURO/CENT_BURO/PERI_BURO/C_BURO. No TMDD/receptor-binding ODEs are
#     present anywhere in the original despite burosumab being a mAb -- FGF23
#     neutralization (FGF23_NEUT, renamed EFFECT_BURO) is a plain algebraic
#     Hill function of C_BURO, not a receptor-occupancy state. EC50_BURO_FGF23/
#     HILL_BURO renamed EC50_BURO/GAMMA_BURO; EMAX_BURO=1 added (implicit in
#     the original, made explicit for the named interface). Rename only, not
#     a refit -- the original was already exactly this shape.
#   - Oral phosphate and oral calcitriol: bespoke variant of the archetypes
#     (depot -> a state that is *itself* the concentration signal PD reads,
#     no separate mass-then-divide-by-volume central compartment) -- kept
#     exactly as the original modeled it. PHOSORAL_GUT/PHOSORAL_SIG renamed
#     GUT_PHOSORAL/C_PHOSORAL; CALC_GUT/CALC_CENT renamed GUT_CALC/C_CALC.
#     PHOS_ORALBOOST/CALC_EXOG_BOOST (already plain Hill ratios, gamma=1
#     implicit) renamed EFFECT_PHOSORAL/EFFECT_CALC; EMAX_PHOS_ORALBOOST/
#     EMAX_CALC_BOOST renamed EMAX_PHOSORAL/EMAX_CALC; GAMMA_PHOSORAL=1 and
#     GAMMA_CALC=1 added explicitly. Rename only, not a refit.
#   - Two pre-existing mrgsolve-2.0.1 build defects were found while
#     verifying (present in the untouched original; not introduced by this
#     refactor; logged to translations/UPSTREAM_ISSUES.md, not fixed
#     upstream). Both required a syntax-only, non-numeric fix directly in
#     this file per the guide's settled policy for this situation:
#       1. `$CAPTURE` listed 14 names that are also `$CMT` compartment
#          names -- removed from `$CAPTURE` (compartments always appear in
#          mrgsolve output regardless of `$CAPTURE` membership).
#       2. `$MAIN`'s `if (NEWIND <= 1) { CMT = value; ... }` idiom assigns
#          directly to compartment names, which are read-only references
#          under `$PLUGIN autodec` in this mrgsolve build -- switched to the
#          `<cmt>_0` initial-value idiom (`NPT2_0 = NPT2_BASE;`, etc.).
#     Neither fix changes any number; see verification result below and
#     xlh_refactor_notes.md for full disclosure.
#
#   Verified via the qspserver mrgsolve_api service (`POST /model_manifest`,
#   `POST /run_simulation`) against the original file's own "Untreated",
#   "Conventional_Ped_PhosCalcitriol" (exercises both oral phosphate and
#   oral calcitriol), and "Burosumab_Adult_1p0mgkg_Q4W" dosing scenarios,
#   1-year daily grid: every shared `$CAPTURE`d output matched the original
#   EXACTLY (max abs diff 0.0) in all three scenarios.
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

xlh_code <- '
$PROB
# X-Linked Hypophosphatemia (XLH) QSP model
# 19 ODE compartments: 7 drug PK + 12 disease/PD/clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Burosumab PK (2-cpt, SC, mAb) — Carpenter 2014 JCI / Carpenter 2018 NEJM
// Archetype 3 (depot + central + peripheral, linear). No TMDD/receptor-
// binding ODEs in the original -- FGF23 neutralization below is a plain
// algebraic Hill function of concentration, not a receptor-occupancy state.
// ============================================
KA_BURO  : 0.017  : Burosumab SC absorption rate (1/h)      // Tmax ~4-7 d
CL_BURO  : 0.0086 : Burosumab linear clearance (L/h)         // t1/2 ~19 d
V1_BURO  : 3.2    : Burosumab central Vd (L)
V2_BURO  : 2.6    : Burosumab peripheral Vd (L)
Q_BURO   : 0.010  : Burosumab inter-compartmental clearance (L/h)
F_BURO   : 0.70   : Burosumab SC bioavailability

// ============================================
// Oral phosphate PK (transient absorption spike) — conventional therapy
// Bespoke: depot + a state that is itself the concentration/exposure
// signal PD reads (no separate mass/volume central compartment) -- kept
// exactly as modeled in the original.
// ============================================
KA_PHOSORAL : 0.90 : Oral phosphate absorption rate (1/h)
KE_PHOSORAL : 0.60 : Oral phosphate "excess" clearance/equilibration rate (1/h)
V_PHOSORAL  : 12.0 : Apparent Vd for oral phosphate excess signal (L)
F_PHOSORAL  : 0.60 : Oral phosphate net fractional absorption

// ============================================
// Oral calcitriol/alfacalcidol PK — conventional therapy
// Bespoke: same concentration-state pattern as oral phosphate above.
// ============================================
KA_CALC  : 0.35  : Oral calcitriol absorption rate (1/h)
KE_CALC  : 0.029 : Oral calcitriol elimination rate (1/h)   // t1/2 ~24 h
V_CALC   : 40    : Apparent Vd calcitriol (L)
F_CALC   : 0.60  : Oral bioavailability

// ============================================
// Disease baseline (untreated XLH natural history)
// ============================================
FGF23_BASE   : 800   : Baseline intact FGF23 (RU/mL; markedly elevated vs ~50 normal)
NPT2_BASE    : 0.40  : Baseline normalized NPT2a/c surface activity (0-1; XLH suppressed)
TMPGFR0      : 1.8   : Baseline TmP/GFR (mg/dL; markedly low vs 2.8-4.2 normal)
PHOS0        : 2.2   : Baseline serum phosphate (mg/dL; hypophosphatemic)
CALCITRIOL0  : 22    : Baseline serum 1,25(OH)2D (pg/mL; inappropriately low/normal)
PTH0         : 65    : Baseline intact PTH (pg/mL; upper-normal/secondary HPT)
BSAP0        : 45    : Baseline bone-specific ALP (µg/L; elevated turnover)
RSS0         : 5.5   : Baseline Rickets Severity Score (Thacher 0-10 scale, pediatric)
AGV_BASE_XLH : 4.3   : Untreated XLH annualized growth velocity (cm/yr)
HEIGHTZ0_XLH : -1.8  : Baseline height Z-score
SIXMWT0      : 380   : Baseline 6-minute walk distance (m, adult)
WOMAC0       : 45    : Baseline WOMAC pain/stiffness score (0-100)
NEPHROCALC0  : 0.05  : Baseline nephrocalcinosis risk index (0-1)
UCACR0       : 0.15  : Baseline urine Ca/Cr ratio (mg/mg)

// ============================================
// PD parameters — Burosumab (FGF23 neutralization)
// Hill interface: rename, not a fit -- the original's FGF23_NEUT was
// already exactly Emax*C^gamma/(EC50^gamma+C^gamma) shaped.
// ============================================
EC50_BURO : 400  : Burosumab conc. for half-max FGF23 neutralization (ng/mL)
GAMMA_BURO       : 1.0  : Hill coefficient
EMAX_BURO        : 1.0  : Max fractional FGF23 neutralization (implicit in the original; named for the interface)
EMAX_NPT2_RESCUE: 0.55 : Max fractional NPT2 surface-activity rescue at full FGF23 neutralization
K_NPT2          : 0.03 : Rate constant, NPT2 activity response (1/h)
K_TMPGFR        : 0.02 : Rate constant, TmP/GFR response (1/h)
TMPGFR_MAX      : 3.4  : Ceiling TmP/GFR under full rescue (mg/dL)
K_PHOS          : 0.08 : Rate constant, serum phosphate equilibration (1/h)
K_CALCITRIOL    : 0.015: Rate constant, endogenous calcitriol response (1/h)
CALCITRIOL_MAX  : 48   : Ceiling endogenous calcitriol under full FGF23 blockade (pg/mL)

// ============================================
// PD parameters — Oral phosphate & calcitriol (conventional therapy, direct substrate)
// Hill interface: rename, not a fit -- both boosts were already plain
// Emax*X/(EC50+X) ratios (gamma=1 implicit) in the original.
// ============================================
EMAX_PHOSORAL : 1.8  : Max transient serum Pi rise per oral phosphate dose (mg/dL)
EC50_PHOSORAL        : 8.0  : Oral phosphate exposure for half-max Pi boost
GAMMA_PHOSORAL       : 1.0  : Hill coefficient (implicit in the original; named for the interface)
EMAX_CALC      : 20   : Max serum calcitriol-equivalent rise from exogenous calcitriol (pg/mL)
EC50_CALC            : 15   : Calcitriol conc. for half-max effect (ng/mL)
GAMMA_CALC           : 1.0  : Hill coefficient (implicit in the original; named for the interface)

// ============================================
// PTH feedback
// ============================================
K_PTH        : 0.05  : Rate constant, PTH response to Ca/Pi/calcitriol (1/h)
PTH_CA_GAIN  : 30    : PTH rise per unit deficit in effective serum calcium signal
KOUT_BSAP    : 0.006 : Rate constant, BSAP decline with mineralization improvement (1/h)
BSAP_MIN     : 18    : Floor BSAP with sustained mineralization rescue

// ============================================
// Growth / clinical translation
// ============================================
K_RSS        : 0.010 : Rate constant, RSS improvement (1/h, slow structural)
RSS_MIN      : 0.5   : Floor RSS under sustained rescue
GAIN_AGV_XLH : 2.6   : Max AGV increment from full mineralization rescue (cm/yr)
K_HEIGHTZ    : 0.10  : Rate constant, height Z-score integration
K_SIXMWT     : 0.008 : Rate constant, 6MWT improvement (1/h)
SIXMWT_GAIN  : 65    : Max 6MWT distance gain (m)
K_WOMAC      : 0.012 : Rate constant, WOMAC pain improvement (1/h)
WOMAC_GAIN   : 22    : Max WOMAC score reduction

// ============================================
// Safety / nephrocalcinosis risk
// ============================================
K_UCACR      : 0.05  : Rate constant, urine Ca/Cr response (1/h)
UCACR_CONV_GAIN : 0.35 : Urine Ca/Cr rise from chronic conventional therapy
K_NEPHROCALC : 0.0008: Nephrocalcinosis risk accumulation rate (per h, driven by UCACR excess)
HYPERPI_THRESH : 5.5 : Serum Pi threshold above which overcorrection risk accrues (mg/dL)

// ============================================
// Adherence
// ============================================
ADHERENCE_BURO : 1.0 : Fraction of scheduled burosumab doses taken
ADHERENCE_CONV : 1.0 : Fraction of scheduled conventional-therapy doses taken (GI tolerability)

$CMT @annotated
GUT_BURO     : Burosumab SC depot (mg)                          // was BURO_DEPOT
CENT_BURO    : Burosumab central compartment (mg)                // was BURO_CENT
PERI_BURO    : Burosumab peripheral compartment (mg)             // was BURO_PERIPH
GUT_PHOSORAL : Oral phosphate gut/absorption depot (mg)          // was PHOSORAL_GUT
C_PHOSORAL   : Oral phosphate transient exposure signal (a.u.)   // was PHOSORAL_SIG
GUT_CALC     : Oral calcitriol gut depot (µg)                    // was CALC_GUT
C_CALC       : Oral calcitriol central concentration (ng/mL equiv) // was CALC_CENT
NPT2         : Normalized NPT2a/c renal surface activity (0-1)
TMPGFR       : TmP/GFR (mg/dL)
PHOS         : Serum phosphate (mg/dL)
CALCITRIOL   : Serum 1,25(OH)2D (pg/mL)
PTH          : Intact PTH (pg/mL)
BSAP         : Bone-specific alkaline phosphatase (µg/L)
RSS          : Rickets Severity Score (0-10)
HEIGHTZ_XLH  : Height Z-score
SIXMWT       : 6-minute walk distance (m)
WOMAC        : WOMAC pain/stiffness score (0-100)
UCACR        : Urine Ca/Cr ratio (mg/mg)
NEPHROCALC   : Nephrocalcinosis risk index (0-1)

$MAIN
F_GUT_BURO       = F_BURO * ADHERENCE_BURO;
F_GUT_PHOSORAL   = F_PHOSORAL * ADHERENCE_CONV;
F_GUT_CALC       = F_CALC * ADHERENCE_CONV;

// Build-compat fix (syntax-only, non-numeric; see UPSTREAM_ISSUES.md and
// xlh_refactor_notes.md): the original assigned directly to compartment
// names here (`NPT2 = NPT2_BASE;`, etc.), which mrgsolve 2.0.1 treats as a
// read-only reference under `$PLUGIN autodec`. Switched to the `<cmt>_0`
// initial-value idiom -- identical effect, compiles under this build.
if (NEWIND <= 1) {
  NPT2_0       = NPT2_BASE;
  TMPGFR_0     = TMPGFR0;
  PHOS_0       = PHOS0;
  CALCITRIOL_0 = CALCITRIOL0;
  PTH_0        = PTH0;
  BSAP_0       = BSAP0;
  RSS_0        = RSS0;
  HEIGHTZ_XLH_0= HEIGHTZ0_XLH;
  SIXMWT_0     = SIXMWT0;
  WOMAC_0      = WOMAC0;
  UCACR_0      = UCACR0;
  NEPHROCALC_0 = NEPHROCALC0;
}

$ODE
// ---- PK: Burosumab (2-cpt, SC) [archetype 3: depot + central + peripheral, linear] ----
dxdt_GUT_BURO  = -KA_BURO * GUT_BURO;
dxdt_CENT_BURO =  KA_BURO * GUT_BURO - CL_BURO/V1_BURO*CENT_BURO
                    - (Q_BURO/V1_BURO)*CENT_BURO + (Q_BURO/V2_BURO)*PERI_BURO;
dxdt_PERI_BURO =  (Q_BURO/V1_BURO)*CENT_BURO - (Q_BURO/V2_BURO)*PERI_BURO;
double C_BURO = CENT_BURO / V1_BURO * 1000.0;   // ng/mL equivalent -- the exposed concentration

// ---- PK: Oral phosphate (transient absorption spike; bespoke -- depot + a directly
//      concentration-valued "central" state, no separate mass/volume split; kept as-is,
//      structure unchanged) ----
dxdt_GUT_PHOSORAL = -KA_PHOSORAL * GUT_PHOSORAL;
dxdt_C_PHOSORAL   =  KA_PHOSORAL * GUT_PHOSORAL / V_PHOSORAL - KE_PHOSORAL * C_PHOSORAL;

// ---- PK: Oral calcitriol (bespoke -- same concentration-state pattern as oral phosphate) ----
dxdt_GUT_CALC = -KA_CALC * GUT_CALC;
dxdt_C_CALC   =  KA_CALC * GUT_CALC / V_CALC - KE_CALC * C_CALC;

// ---- PD: FGF23 neutralization by burosumab (Hill interface; already this shape in the
//      original -- EMAX_BURO=1 made explicit, rename not a refit) ----
double EFFECT_BURO = EMAX_BURO * pow(C_BURO, GAMMA_BURO) / (pow(EC50_BURO, GAMMA_BURO) + pow(C_BURO, GAMMA_BURO));

// ---- NPT2a/c surface activity: rescued by FGF23 neutralization ----
double NPT2_TARGET = NPT2_BASE + EMAX_NPT2_RESCUE * EFFECT_BURO;
if (NPT2_TARGET > 1.0) NPT2_TARGET = 1.0;
dxdt_NPT2 = K_NPT2 * (NPT2_TARGET - NPT2);

// ---- TmP/GFR: tracks NPT2 activity ----
double TMPGFR_TARGET = TMPGFR0 + (TMPGFR_MAX - TMPGFR0) * (NPT2 - NPT2_BASE) / (1.0 - NPT2_BASE + 1e-6);
dxdt_TMPGFR = K_TMPGFR * (TMPGFR_TARGET - TMPGFR);

// ---- Serum phosphate: driven by TmP/GFR rescue (burosumab, sustained) + transient oral phosphate spike (conventional) ----
double EFFECT_PHOSORAL = EMAX_PHOSORAL * pow(C_PHOSORAL, GAMMA_PHOSORAL) / (pow(EC50_PHOSORAL, GAMMA_PHOSORAL) + pow(C_PHOSORAL, GAMMA_PHOSORAL));
double PHOS_TARGET = PHOS0 + (TMPGFR - TMPGFR0) * 1.15 + EFFECT_PHOSORAL;
dxdt_PHOS = K_PHOS * (PHOS_TARGET - PHOS);

// ---- Endogenous calcitriol: de-repressed by FGF23 neutralization; exogenous calcitriol adds directly ----
double CALCITRIOL_TARGET = CALCITRIOL0 + (CALCITRIOL_MAX - CALCITRIOL0) * EFFECT_BURO;
double EFFECT_CALC = EMAX_CALC * pow(C_CALC, GAMMA_CALC) / (pow(EC50_CALC, GAMMA_CALC) + pow(C_CALC, GAMMA_CALC));
dxdt_CALCITRIOL = K_CALCITRIOL * (CALCITRIOL_TARGET - CALCITRIOL) + 0.02*EFFECT_CALC;

// ---- PTH: rises with low serum phosphate/calcitriol, falls as both normalize ----
double PTH_TARGET = PTH0 + PTH_CA_GAIN * (-(PHOS - PHOS0)/4.0) - 0.4*(CALCITRIOL - CALCITRIOL0);
if (PTH_TARGET < 10) PTH_TARGET = 10;
dxdt_PTH = K_PTH * (PTH_TARGET - PTH);

// ---- BSAP: declines as mineralization improves (chronic phosphate/calcitriol sufficiency) ----
double BSAP_TARGET = BSAP_MIN + (BSAP0 - BSAP_MIN) * (1.0 - (PHOS - PHOS0)/(TMPGFR_MAX - TMPGFR0 + 1e-6));
if (BSAP_TARGET < BSAP_MIN) BSAP_TARGET = BSAP_MIN;
dxdt_BSAP = KOUT_BSAP * (BSAP_TARGET - BSAP);

// ---- RSS: slow structural improvement tracking chronic phosphate sufficiency ----
double MINERAL_INDEX = (PHOS - PHOS0) / (4.5 - PHOS0 + 1e-6);
if (MINERAL_INDEX < 0) MINERAL_INDEX = 0;
if (MINERAL_INDEX > 1) MINERAL_INDEX = 1;
double RSS_TARGET = RSS0 - (RSS0 - RSS_MIN) * MINERAL_INDEX;
dxdt_RSS = K_RSS * (RSS_TARGET - RSS);

// ---- Growth: AGV/height Z-score respond to mineralization index (pediatric) ----
double AGV_CALC_XLH = AGV_BASE_XLH + GAIN_AGV_XLH * MINERAL_INDEX;
dxdt_HEIGHTZ_XLH = K_HEIGHTZ * ( (AGV_CALC_XLH - AGV_BASE_XLH)/3.0 - (HEIGHTZ_XLH - HEIGHTZ0_XLH)*0.03 );

// ---- 6MWT / WOMAC: functional & pain endpoints track mineralization + myopathy improvement (adult) ----
double SIXMWT_TARGET = SIXMWT0 + SIXMWT_GAIN * MINERAL_INDEX;
dxdt_SIXMWT = K_SIXMWT * (SIXMWT_TARGET - SIXMWT);
double WOMAC_TARGET = WOMAC0 - WOMAC_GAIN * MINERAL_INDEX;
if (WOMAC_TARGET < 5) WOMAC_TARGET = 5;
dxdt_WOMAC = K_WOMAC * (WOMAC_TARGET - WOMAC);

// ---- Nephrocalcinosis risk: driven by hypercalciuria from conventional therapy (oral Ca load via calcitriol) + phosphate overcorrection ----
double UCACR_TARGET = UCACR0 + UCACR_CONV_GAIN * EFFECT_CALC / EMAX_CALC;
dxdt_UCACR = K_UCACR * (UCACR_TARGET - UCACR);
double HYPERPI_EXCESS = PHOS > HYPERPI_THRESH ? (PHOS - HYPERPI_THRESH) : 0.0;
dxdt_NEPHROCALC = K_NEPHROCALC * ( UCACR + HYPERPI_EXCESS*0.5 ) * (1.0 - NEPHROCALC);

// Build-compat fix (syntax-only, non-numeric; see UPSTREAM_ISSUES.md and
// xlh_refactor_notes.md): the original listed 14 names in $CAPTURE that are
// also $CMT compartment names (NPT2, TMPGFR, PHOS, CALCITRIOL, PTH, BSAP,
// RSS, HEIGHTZ_XLH, SIXMWT, WOMAC, UCACR, NEPHROCALC, and -- after renaming
// -- C_PHOSORAL, C_CALC). mrgsolve 2.0.1 refuses to build any model whose
// $CAPTURE repeats a compartment name. Removed here; every compartment
// state is included in mrgsolve output regardless of $CAPTURE membership,
// so this changes nothing about what is reported.
$CAPTURE C_BURO EFFECT_BURO AGV_CALC_XLH EFFECT_PHOSORAL EFFECT_CALC
'

xlh_mod <- mcode("xlh_qsp", xlh_code)

# =====================================================================
# Treatment scenarios (10) — dosing via event tables
#   Pediatric simulations: 5-year-old XLH patient, WT ~18 kg, 52-week horizon
#   Adult simulations: WT ~70 kg, 24-week (AXLES1-comparable) or 52-week horizon
# =====================================================================
WT_PED   <- 18  # kg, representative 5-yr-old XLH body weight
WT_ADULT <- 70  # kg, representative adult body weight

make_ev <- function(amt, ii, addl, cmt) {
  ev(time = 0, amt = amt, ii = ii, addl = addl, cmt = cmt)
}

scenarios <- list(
  "1_Untreated_NaturalHistory"        = NULL,
  "2_Conventional_Ped_PhosCalcitriol" = c(
                                          make_ev(10 * WT_PED, 6, 4*365, "GUT_PHOSORAL"),  # ~10 mg/kg QID oral phosphate
                                          make_ev(30 * WT_PED/1000, 24, 364, "GUT_CALC")    # ~30 ng/kg/day calcitriol, QD depot equivalent
                                        ),
  "3_Burosumab_Ped_0p8mgkg_Q2W"       = make_ev(0.8 * WT_PED, 336, 25, "GUT_BURO"),
  "4_Burosumab_Ped_2p0mgkg_Q2W"       = make_ev(2.0 * WT_PED, 336, 25, "GUT_BURO"),
  "5_Burosumab_Adult_1p0mgkg_Q4W"     = make_ev(1.0 * WT_ADULT, 672, 12, "GUT_BURO"),
  "6_Switch_Conventional_to_Burosumab"= make_ev(0.8 * WT_PED, 336, 25, "GUT_BURO"),   # CL303-style switch: modeled as burosumab from time 0 (post-switch phase)
  "7_Conventional_PoorAdherence_GI"   = make_ev(10 * WT_PED, 6, 4*365, "GUT_PHOSORAL"),
  "8_Burosumab_Supratherapeutic_Overcorrection" = make_ev(3.5 * WT_PED, 336, 25, "GUT_BURO"),
  "9_Conventional_LongTerm_TertiaryHPT" = make_ev(30 * WT_PED/1000, 24, 4*364, "GUT_CALC"),
  "10_Burosumab_Adult_LowAdherence_60pct" = make_ev(1.0 * WT_ADULT, 672, 12, "GUT_BURO")
)

run_scenario <- function(name, ev_obj, adherence_buro = 1.0, adherence_conv = 1.0, end = 8760) {
  m <- xlh_mod %>% param(ADHERENCE_BURO = adherence_buro, ADHERENCE_CONV = adherence_conv)
  if (!is.null(ev_obj)) {
    out <- m %>% ev(ev_obj) %>% mrgsim(end = end, delta = 24) %>% as_tibble()
  } else {
    out <- m %>% mrgsim(end = end, delta = 24) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Example run (uncomment to execute):
# results <- bind_rows(
#   run_scenario("1_Untreated_NaturalHistory", NULL),
#   run_scenario("2_Conventional_Ped_PhosCalcitriol", scenarios[["2_Conventional_Ped_PhosCalcitriol"]]),
#   run_scenario("3_Burosumab_Ped_0p8mgkg_Q2W", scenarios[["3_Burosumab_Ped_0p8mgkg_Q2W"]]),
#   run_scenario("4_Burosumab_Ped_2p0mgkg_Q2W", scenarios[["4_Burosumab_Ped_2p0mgkg_Q2W"]]),
#   run_scenario("5_Burosumab_Adult_1p0mgkg_Q4W", scenarios[["5_Burosumab_Adult_1p0mgkg_Q4W"]]),
#   run_scenario("7_Conventional_PoorAdherence_GI", scenarios[["7_Conventional_PoorAdherence_GI"]], adherence_conv = 0.6),
#   run_scenario("8_Burosumab_Supratherapeutic_Overcorrection", scenarios[["8_Burosumab_Supratherapeutic_Overcorrection"]]),
#   run_scenario("9_Conventional_LongTerm_TertiaryHPT", scenarios[["9_Conventional_LongTerm_TertiaryHPT"]]),
#   run_scenario("10_Burosumab_Adult_LowAdherence_60pct", scenarios[["10_Burosumab_Adult_LowAdherence_60pct"]], adherence_buro = 0.6)
# )
#
# ggplot(results, aes(time/24, RSS, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Day", y = "Rickets Severity Score", title = "XLH: RSS trajectory by scenario")

# =====================================================================
# Calibration notes:
#  - Untreated FGF23_BASE (~800 RU/mL), PHOS0 (2.2 mg/dL), TMPGFR0 (1.8 mg/dL)
#    reflect typical baseline XLH biochemistry vs. age-matched normal ranges
#    (Carpenter 2011 JBMR clinician's guide; Carpenter 2018 NEJM baseline).
#  - EMAX_NPT2_RESCUE / EC50_BURO calibrated so that burosumab 0.8-2.0
#    mg/kg Q2W raises TmP/GFR and serum Pi into the low-normal range within
#    ~4-12 weeks, consistent with Carpenter 2018 NEJM dose-finding and the
#    Imel 2019 Lancet (CL303) phase-3 RCT superiority over conventional
#    therapy on RSS/RGI-C and growth at week 40/64.
#  - Insogna 2018 JBMR (AXLES1) adult 1.0 mg/kg Q4W week-24 primary analysis
#    anchors SIXMWT_GAIN / WOMAC_GAIN magnitude and time course.
#  - Conventional therapy (oral phosphate + calcitriol) produces only a
#    transient serum Pi spike (C_PHOSORAL fast decay) rather than sustained
#    TmP/GFR correction, reflecting its bypass of the FGF23-NPT2 axis --
#    consistent with Imel 2019 Lancet showing inferior/incomplete rickets
#    healing vs burosumab despite chronic multi-dose-daily dosing.
#  - UCACR_CONV_GAIN / K_NEPHROCALC calibrated so chronic conventional
#    therapy (scenario 9) accumulates nephrocalcinosis risk over years,
#    consistent with reported 20-80% pediatric prevalence under conventional
#    therapy vs. low incremental risk with burosumab monotherapy.
#  - Burosumab 2-cpt PK (KA/CL/V1/V2/Q) approximates the ~19-day terminal
#    half-life and dose-proportional exposure reported in Carpenter 2014 JCI
#    (KRN23 first-in-human) and subsequent population PK analyses.
# =====================================================================
