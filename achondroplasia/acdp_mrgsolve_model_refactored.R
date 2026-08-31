# =====================================================================
# Achondroplasia (ACH) — mrgsolve QSP Model (REFACTORED, pluggable PK)
#   Refactor of `acdp_mrgsolve_model.R` per FORK_WORKFLOW_GUIDE.md Part 2.
#   Original untouched — see `acdp_refactor_notes.md` for full rationale.
#
#   Scope: all three of this file's census rows in
#   `driver-patches/data/compound_perturbation_census.md` are genuine,
#   externally-dosed drugs (each has its own depot->central PK compartment
#   pair and is actually dosed by an `ev()` scenario somewhere in the
#   original file) — Vosoritide (VOS), TransCon CNP/navepegritide (TCNP),
#   Infigratinib (INFIG). "Released free-CNP (TCNP)" is NOT endogenous CNP:
#   it is the pharmacologically active moiety released from a dosed,
#   sustained-release SC prodrug depot (scenario 6 doses `TCNP_DEPOT`
#   directly) — see acdp_refactor_notes.md for the full check. The census
#   has been updated accordingly (all three rows marked as refactored
#   drugs, not "delete PK compartment").
#
#   Build-compatibility fixes (disclosed, syntax-only, non-numeric — see
#   notes and UPSTREAM_ISSUES.md #124): the ORIGINAL does not compile
#   under mrgsolve 2.0.1, for two independent reasons, both fixed here:
#     1. `$CAPTURE` re-listed fifteen `$CMT` compartment names directly
#        (`PERK CGMP_SIG CHONDRO HEIGHT_CM HEIGHTZ FMAREA SPCANALZ AHI
#        OTITIS BMIZ MAP_BP HR PHOS VOS_CP TCNP_CP INFIG_CP`), which
#        mrgsolve 2.0.1 rejects. Dropped (compartments are exposed as
#        output columns automatically); replaced with this refactor's own
#        `C_<STEM>`/`EFFECT_<STEM>` capture entries.
#     2. `$MAIN`'s `if (NEWIND <= 1) { ... }` block set thirteen disease
#        compartments' initial conditions with bare compartment-name
#        assignment (`PERK = PERK_BASE;` etc.) instead of the standard
#        `<cmt>_0 =` idiom. Once (1) is fixed, this fails to compile
#        ("assignment of read-only reference"). Switched to `PERK_0 =
#        PERK_BASE;` etc. — same values, same semantics.
#   Both fixes are confirmed syntax-only/non-numeric by the verification
#   below (exact match, max abs diff 0.0, against a patched-original
#   scratch copy carrying only these same two fixes).
#
#   A third, pre-existing finding (not a build defect, not fixed, disclosed
#   only — UPSTREAM_ISSUES.md #125): the CNP-axis Hill term uses a
#   non-integer exponent (GAMMA_VOS = 1.5, was HILL_VOS) applied to
#   CNP_TOTAL = C_VOS + C_TCNP. Once Vosoritide (fast, ~15-min half-life)
#   decays within a dosing interval to a concentration many orders of
#   magnitude below EC50_VOS, ordinary floating-point/adaptive-solver
#   roundoff can push the true near-zero value slightly *negative*;
#   `pow(negative, 1.5)` is `NaN` in C++, and the NaN then poisons the
#   whole model forever (CGMP_SIG -> PERK -> every downstream output).
#   Confirmed identical in the untouched original (same formula, same
#   trigger) — see notes. Every Vosoritide-dosing scenario in this file
#   hits it at simulated hour 8, identically in both models.
#
#   PK/PD  : Vosoritide (CNP analog, SC QD) · TransCon CNP/navepegritide
#            (sustained-release CNP prodrug, SC QW) · Infigratinib
#            (FGFR1-3 TKI, PO QD) · Growth hormone (off-label, historical,
#            flag-only — no PK compartment in the original, untouched)
#   Outputs: Annualized growth velocity (AGV), cumulative height, height
#            Z-score, foramen magnum area, spinal canal Z-score, OSA-AHI,
#            otitis media rate, BMI-Z, hemodynamic/safety signals, serum
#            phosphate (FGFR1-off-target).
#   References (calibration): Savarirayan et al. NEJM 2019 (PMID 31269546;
#            phase 2 dose-finding, 2.5/7.5/15/30 µg/kg), Savarirayan et al.
#            Lancet 2020 (PMID 32891212; phase 3, 52-wk placebo-controlled,
#            ΔAGV +1.57 cm/yr), Savarirayan et al. 2021 phase 3 extension
#            (PMID 34341520), Ascendis Pharma ApproaCH phase 3 (TransCon
#            CNP/navepegritide; topline Sep-2024, LS-mean ΔAGV +1.49 cm/yr;
#            FDA-approved Feb-2026 as YUVIWEL for age ≥2 yr, EMA decision
#            expected Q4-2026), BridgeBio/QED PROPEL 2 phase 2 (infigratinib;
#            NEJM 2025, PMID 39555818) and PROPEL 3 phase 3 (NEJM 2026;
#            best-in-class AGV + first significant body-proportionality
#            improvement; well tolerated, mild/transient asymptomatic
#            hyperphosphatemia ~4%, no discontinuations; NDA planned Q3-2026,
#            not yet approved), Horton 1978 (PMID 690757; ACH-specific
#            growth curves), Hunter 1998 (foramen magnum/cervicomedullary
#            natural history), White 2020 (AAP ACH health supervision
#            guideline).
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

acdp_code <- '
$PROB
# Achondroplasia (ACH) QSP model -- PK/PD refactor (pluggable PK)
# 20 ODE compartments: 6 drug PK + 13 disease/PD/clinical
# Renamed to this fork PK/PD naming convention: GUT_<STEM>/CENT_<STEM>,
# CL_<STEM>/V1_<STEM>/KA_<STEM>, C_<STEM> (exposed concentration),
# EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>/EFFECT_<STEM> (Hill interface).
# All three census-row compounds (Vosoritide=VOS, TransCon CNP=TCNP,
# Infigratinib=INFIG) are genuinely externally-dosed drugs -- see
# acdp_refactor_notes.md for the per-compound dosing-route check.

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Vosoritide PK (no depot->central linear elim; bespoke: state IS conc.) -- Savarirayan 2019/2020
// ============================================
KA_VOS  : 6.0   : Vosoritide SC absorption rate (1/h) // Tmax ~10-15 min
CL_VOS  : 2.77  : Vosoritide elimination rate constant (1/h)   // t1/2 ~15 min (NPR-C/NEP); was KE_VOS
V1_VOS  : 15.0  : Vosoritide apparent Vd (L)
F_VOS   : 0.70  : Vosoritide SC bioavailability

// ============================================
// TransCon CNP (navepegritide) PK -- sustained-release prodrug depot, SC QW
// ============================================
KA_TCNP : 0.018 : Prodrug release rate from SC depot (1/h) // sustained over ~1 wk; was KREL_TCNP -- release, not absorption, but same depot->central structural role
CL_TCNP : 2.5    : Released free-CNP elimination rate constant (1/h)  // fast, like native CNP; was KE_TCNP
V1_TCNP : 15.0   : Free-CNP apparent Vd (L)

// ============================================
// Infigratinib PK (depot->central, linear elim) -- QED/BridgeBio PROPEL
// ============================================
KA_INFIG  : 0.35  : Oral absorption rate (1/h)
CL_INFIG  : 7.7   : Apparent clearance (L/h)
V1_INFIG  : 480   : Apparent Vd (L)                    // t1/2 ~20-24 h
F_INFIG   : 0.60  : Oral bioavailability

// ============================================
// Growth hormone (off-label, historical use) -- flag only, no PK compartment in the original
// ============================================
GH_EFFECT : 0.06  : Fractional AGV boost from exogenous GH (small, historical)

// ============================================
// Disease baseline (untreated ACH natural history)
// ============================================
PERK_BASE   : 1.00  : Baseline normalized pERK activity (ACH, hyperactive)
PERK_NORMAL : 0.35  : Reference pERK activity (non-ACH)
AGV_BASE    : 3.90  : Untreated ACH annualized growth velocity (cm/yr) // Horton curves / placebo arm
HEIGHT0     : 85.0  : Starting height at model entry (cm; ~age 5)
HEIGHTZ0    : -5.0  : Starting height Z-score (ACH-specific curve)
FMAREA0     : 280   : Baseline foramen magnum area (mm^2; ACH reduced vs ~450 normal)
SPCANALZ0   : -2.5  : Baseline spinal canal diameter Z-score
AHI0        : 4.0   : Baseline OSA-AHI (events/h; mixed obstructive+central)
OTITIS0     : 4.5   : Baseline otitis media episodes/yr
BMIZ0       : 0.8   : Baseline BMI Z-score (ACH-specific; obesity-prone)
MAP0        : 82    : Baseline mean arterial pressure (mmHg)
HR0         : 90    : Baseline heart rate (bpm)
PHOS0       : 4.2   : Baseline serum phosphate (mg/dL)

// ============================================
// Hill interface -- CNP/NPR-B axis (Vosoritide + TransCon CNP; shared receptor pool, see notes)
// ============================================
EC50_VOS    : 8.0   : CNP-class conc. for half-max cGMP signal (ng/mL)
GAMMA_VOS   : 1.5   : Hill coefficient, CNP-NPR-B // was HILL_VOS
EMAX_VOS    : 1.0   : Hill ceiling [math-implied: original ratio saturates at 1; EMAX_PERKINH below applies the real downstream ceiling, unchanged]
EMAX_PERKINH: 0.42  : Max fractional pERK inhibition via PKGII-RAF1 (15 ug/kg plateau)
KOUT_CGMP   : 3.0   : cGMP signal decay rate (1/h)

// ============================================
// Hill interface -- Infigratinib (direct FGFR3 kinase block)
// ============================================
EC50_INFIG  : 25.0  : Infigratinib conc. for half-max FGFR3 inhibition (ng/mL)
GAMMA_INFIG : 1.2   : Hill coefficient // was HILL_INFIG
EMAX_INFIG  : 0.55  : Max fractional pERK inhibition (direct kinase block) // was EMAX_INFIG_PERK
K_OFFTARGET_FGFR1: 0.30 : Fractional off-target FGFR1 inhibition (growth-plate/renal-phosphate)

// ============================================
// Growth-plate -> clinical translation (calibrated to trial delta-AGV)
// ============================================
GAIN_AGV      : 4.20  : Max achievable AGV increment from full pERK normalization (cm/yr)
KOUT_GROWTH   : 0.30  : Rate constant, height integration smoothing (1/yr equiv, converted /h)
K_FM          : 0.015 : Rate constant, foramen-magnum-area response to chronic pERK rescue (1/h, slow)
FMAREA_MAX    : 420   : Ceiling foramen magnum area under sustained rescue (mm^2)
K_SPCANAL     : 0.010 : Rate constant, spinal canal Z-score slow response (1/h)
SPCANALZ_MAX  : -1.0  : Ceiling spinal canal Z-score under sustained rescue
K_OTITIS      : 0.05  : Rate constant, otitis frequency decline with age/growth (1/h)
K_BMI         : 0.02  : Rate constant, BMI-Z response to mobility improvement (1/h)
BMIZ_MIN      : 0.2   : Floor BMI-Z with improved mobility

// ============================================
// Hemodynamic / safety PD (CNP class vasodilatory effect)
// ============================================
EMAX_MAP_DROP : 8.0   : Max transient MAP drop with CNP-class agonism (mmHg)
EC50_MAP      : 10.0  : Conc. for half-max MAP effect (ng/mL)
KOUT_MAP      : 1.5   : MAP recovery rate (1/h)
EMAX_HR_RISE  : 12.0  : Max reflex tachycardia (bpm)
KOUT_HR       : 1.2   : HR recovery rate (1/h)
EMAX_PHOS_RISE: 1.3   : Max serum phosphate rise from off-target FGFR1 block (mg/dL)
KOUT_PHOS     : 0.10  : Phosphate equilibration rate (1/h)

// ============================================
// Adherence / dosing fidelity
// ============================================
ADHERENCE   : 1.0   : Fraction of scheduled vosoritide doses actually taken (1.0 = full)
GH_ON       : 0     : Flag (0/1), exogenous growth hormone co-administration

$CMT @annotated
GUT_VOS     : Vosoritide SC depot (ug)               // was VOS_DEPOT
CENT_VOS    : Vosoritide central concentration (ng/mL equiv, amount/V) // was VOS_CP
GUT_TCNP    : TransCon CNP SC prodrug depot (ug)      // was TCNP_DEPOT
CENT_TCNP   : Free CNP moiety central concentration (ng/mL equiv) // was TCNP_CP
GUT_INFIG   : Infigratinib gut compartment (mg)       // was INFIG_GUT
CENT_INFIG  : Infigratinib central concentration (ng/mL equiv, amount/V) // was INFIG_CP
PERK        : Normalized pERK (MAPK) activity
CGMP_SIG    : cGMP/PKGII counter-regulatory signal (a.u.)
CHONDRO     : Growth-plate chondrocyte proliferation index (a.u., 0-1)
HEIGHT_CM   : Cumulative height (cm)
HEIGHTZ     : Height Z-score (ACH-specific curve)
FMAREA      : Foramen magnum area (mm^2)
SPCANALZ    : Spinal canal diameter Z-score
AHI         : OSA apnea-hypopnea index (events/h)
OTITIS      : Otitis media episode rate (episodes/yr)
BMIZ        : BMI Z-score
MAP_BP      : Mean arterial pressure (mmHg)
HR          : Heart rate (bpm)
PHOS        : Serum phosphate (mg/dL)

$MAIN
F_GUT_VOS   = F_VOS * ADHERENCE;
F_GUT_INFIG = F_INFIG;

// Build-defect workaround (mrgsolve 2.0.1): the original assigns these
// initial conditions with bare compartment-name references
// (PERK = PERK_BASE; etc.), which this mrgsolve build rejects with
// "assignment of read-only reference" once the also-broken $CAPTURE
// duplication (see notes) is fixed. Switched to the standard <cmt>_0
// idiom -- same values, same semantics, confirmed via /run_simulation to
// initialize identically to the intended targets. See acdp_refactor_notes.md
// and translations/UPSTREAM_ISSUES.md #124.
if (NEWIND <= 1) {
  PERK_0      = PERK_BASE;
  CGMP_SIG_0  = 0;
  CHONDRO_0   = 1 - (PERK_BASE - PERK_NORMAL);
  HEIGHT_CM_0 = HEIGHT0;
  HEIGHTZ_0   = HEIGHTZ0;
  FMAREA_0    = FMAREA0;
  SPCANALZ_0  = SPCANALZ0;
  AHI_0       = AHI0;
  OTITIS_0    = OTITIS0;
  BMIZ_0      = BMIZ0;
  MAP_BP_0    = MAP0;
  HR_0        = HR0;
  PHOS_0      = PHOS0;
}

$ODE
// ---- PK ----
dxdt_GUT_VOS  = -KA_VOS * GUT_VOS;
dxdt_CENT_VOS =  KA_VOS * GUT_VOS / V1_VOS - CL_VOS * CENT_VOS;

dxdt_GUT_TCNP  = -KA_TCNP * GUT_TCNP;
dxdt_CENT_TCNP =  KA_TCNP * GUT_TCNP / V1_TCNP - CL_TCNP * CENT_TCNP;

dxdt_GUT_INFIG  = -KA_INFIG * GUT_INFIG;
dxdt_CENT_INFIG =  KA_INFIG * GUT_INFIG / V1_INFIG - (CL_INFIG / V1_INFIG) * CENT_INFIG;

// ---- Exposed concentrations (bespoke: each central compartment state
// is already the concentration, not an amount -- see notes) ----
double C_VOS   = CENT_VOS;
double C_TCNP  = CENT_TCNP;
double C_INFIG = CENT_INFIG;

// ---- CNP/NPR-B -> cGMP -> PKGII-RAF1 counter-signal ----
// Bespoke: Vosoritide and TransCon-released free CNP are both direct
// NPR-B agonists (analog / native ligand for the same receptor pool), so
// the original sums their concentration before the single receptor-Hill
// nonlinearity, rather than computing two independent Hill effects. Kept
// as one shared named interface (EFFECT_VOS), TCNP contributing no
// independent EFFECT_TCNP -- see acdp_refactor_notes.md for the same
// reasoning as the acne-vulgaris ISO/OXO precedent.
double CNP_TOTAL = C_VOS + C_TCNP;
double EFFECT_VOS = EMAX_VOS * pow(CNP_TOTAL, GAMMA_VOS) / (pow(EC50_VOS, GAMMA_VOS) + pow(CNP_TOTAL, GAMMA_VOS));
dxdt_CGMP_SIG = KOUT_CGMP * (EFFECT_VOS - CGMP_SIG);

// ---- pERK: baseline hyperactivation, reduced by CNP-axis signal AND/OR direct FGFR3 block ----
double EFFECT_INFIG = EMAX_INFIG * pow(C_INFIG, GAMMA_INFIG) / (pow(EC50_INFIG, GAMMA_INFIG) + pow(C_INFIG, GAMMA_INFIG));
double CNP_INHIB    = EMAX_PERKINH * CGMP_SIG;
double TOTAL_INHIB  = 1 - (1 - CNP_INHIB) * (1 - EFFECT_INFIG);   // combine non-additively (no monotherapy combo in practice)
dxdt_PERK = 2.0 * ( PERK_BASE * (1 - TOTAL_INHIB) - PERK );

// ---- Chondrocyte proliferation rescue (inversely tied to pERK) ----
double CHONDRO_TARGET = 1 - (PERK - PERK_NORMAL);
if (CHONDRO_TARGET < 0) CHONDRO_TARGET = 0;
if (CHONDRO_TARGET > 1) CHONDRO_TARGET = 1;
dxdt_CHONDRO = 1.5 * (CHONDRO_TARGET - CHONDRO);

// ---- Growth translation: instantaneous AGV (cm/yr) drives cumulative height ----
// AGV rises monotonically with CHONDRO rescue relative to the untreated baseline CHONDRO0
double CHONDRO0  = 1 - (PERK_BASE - PERK_NORMAL);
double AGV_CALC  = AGV_BASE + GAIN_AGV * (CHONDRO - CHONDRO0) + AGV_BASE * GH_EFFECT * (GH_ON > 0 ? 1.0 : 0.0);
dxdt_HEIGHT_CM = AGV_CALC / 8760.0;   // cm/yr -> cm/h
dxdt_HEIGHTZ   = 0.15 * ( (AGV_CALC - AGV_BASE) / 4.0 - (HEIGHTZ - HEIGHTZ0)*0.02 );

// ---- Skull base / spine: slow structural response to sustained chondrocyte rescue ----
double FM_TARGET = FMAREA0 + (FMAREA_MAX - FMAREA0) * (CHONDRO - CHONDRO0) / (1 - CHONDRO0 + 1e-6);
dxdt_FMAREA = K_FM * (FM_TARGET - FMAREA);

double SPC_TARGET = SPCANALZ0 + (SPCANALZ_MAX - SPCANALZ0) * (CHONDRO - CHONDRO0) / (1 - CHONDRO0 + 1e-6);
dxdt_SPCANALZ = K_SPCANAL * (SPC_TARGET - SPCANALZ);

// ---- OSA-AHI: mild improvement with growth-plate rescue (airway/midface), otherwise stable ----
double AHI_TARGET = AHI0 - 1.0 * (CHONDRO - CHONDRO0)/(1-CHONDRO0+1e-6);
if (AHI_TARGET < 1.0) AHI_TARGET = 1.0;
dxdt_AHI = 0.02 * (AHI_TARGET - AHI);

// ---- Otitis media: declines with age (Eustachian maturation), minor drug modulation ----
dxdt_OTITIS = -K_OTITIS * (OTITIS - 1.0);

// ---- BMI-Z: modulated by mobility improvement proxy (spinal/limb rescue) ----
double BMIZ_TARGET = BMIZ0 - (BMIZ0 - BMIZ_MIN) * (CHONDRO - CHONDRO0)/(1-CHONDRO0+1e-6);
dxdt_BMIZ = K_BMI * (BMIZ_TARGET - BMIZ);

// ---- Hemodynamic safety: transient CNP-class vasodilation / reflex tachycardia ----
// (a distinct physiological readout from the EFFECT_VOS receptor-Hill above;
// uses the shared CNP_TOTAL ligand pool directly, as the original did)
double MAP_DROP = EMAX_MAP_DROP * CNP_TOTAL / (EC50_MAP + CNP_TOTAL);
dxdt_MAP_BP = KOUT_MAP * ( (MAP0 - MAP_DROP) - MAP_BP );
double HR_RISE = EMAX_HR_RISE * (MAP0 - MAP_BP) / (EMAX_MAP_DROP + 1e-6);
dxdt_HR = KOUT_HR * ( (HR0 + HR_RISE) - HR );

// ---- Off-target FGFR1 (infigratinib class): serum phosphate rise ----
double PHOS_RISE = EMAX_PHOS_RISE * K_OFFTARGET_FGFR1 * EFFECT_INFIG / EMAX_INFIG;
dxdt_PHOS = KOUT_PHOS * ( (PHOS0 + PHOS_RISE) - PHOS );

$CAPTURE AGV_CALC C_VOS C_TCNP C_INFIG EFFECT_VOS EFFECT_INFIG CNP_TOTAL
'

acdp_mod <- mcode("acdp_qsp_refactored", acdp_code)

# =====================================================================
# Treatment scenarios (10) — dosing via event tables
#   All simulations: 5-year-old ACH patient, 52-week (1 yr) horizon
#   unless noted; weight assumed ~15 kg for µg/kg / mg/kg conversions.
#   Compartment targets renamed only (same doses/timing as the original):
#     VOS_DEPOT -> GUT_VOS, TCNP_DEPOT -> GUT_TCNP, INFIG_GUT -> GUT_INFIG
# =====================================================================
WT <- 15  # kg, representative 5-yr-old ACH body weight

make_ev <- function(amt, ii, addl, cmt, tinf = 0) {
  ev(time = 0, amt = amt, ii = ii, addl = addl, cmt = cmt)
}

scenarios <- list(
  "1_Untreated_NaturalHistory" = NULL,
  "2_Vosoritide_15ugkg_QD"     = make_ev(15 * WT, 24, 364, "GUT_VOS"),
  "3_Vosoritide_2p5ugkg_QD"    = make_ev(2.5 * WT, 24, 364, "GUT_VOS"),
  "4_Vosoritide_7p5ugkg_QD"    = make_ev(7.5 * WT, 24, 364, "GUT_VOS"),
  "5_Vosoritide_30ugkg_QD"     = make_ev(30 * WT, 24, 364, "GUT_VOS"),
  "6_TransConCNP_QW"           = make_ev(100 * WT, 168, 51, "GUT_TCNP"),
  "7_Infigratinib_PO_QD"       = make_ev(0.5 * WT, 24, 364, "GUT_INFIG"),
  "8_GrowthHormone_offlabel"   = NULL,   # modeled via GH_ON flag / param override, no PK cmt
  "9_Vosoritide_plus_FMDsurgery" = make_ev(15 * WT, 24, 364, "GUT_VOS"),
  "10_Vosoritide_PoorAdherence_60pct" = make_ev(15 * WT, 24, 364, "GUT_VOS")
)

run_scenario <- function(name, ev, gh_on = 0, adherence = 1.0) {
  m <- acdp_mod %>% param(GH_ON = gh_on, ADHERENCE = adherence)
  if (!is.null(ev)) {
    if (grepl("PoorAdherence", name)) {
      # simulate 60% adherence by dropping 40% of scheduled doses stochastically
      set.seed(42)
      full <- ev
      out <- m %>% ev(full) %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
    } else {
      out <- m %>% ev(ev) %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
    }
  } else {
    out <- m %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Example run (uncomment to execute):
# results <- bind_rows(
#   run_scenario("1_Untreated_NaturalHistory", NULL),
#   run_scenario("2_Vosoritide_15ugkg_QD", scenarios[["2_Vosoritide_15ugkg_QD"]]),
#   run_scenario("3_Vosoritide_2p5ugkg_QD", scenarios[["3_Vosoritide_2p5ugkg_QD"]]),
#   run_scenario("4_Vosoritide_7p5ugkg_QD", scenarios[["4_Vosoritide_7p5ugkg_QD"]]),
#   run_scenario("5_Vosoritide_30ugkg_QD", scenarios[["5_Vosoritide_30ugkg_QD"]]),
#   run_scenario("6_TransConCNP_QW", scenarios[["6_TransConCNP_QW"]]),
#   run_scenario("7_Infigratinib_PO_QD", scenarios[["7_Infigratinib_PO_QD"]]),
#   run_scenario("8_GrowthHormone_offlabel", NULL, gh_on = 1),
#   run_scenario("9_Vosoritide_plus_FMDsurgery", scenarios[["9_Vosoritide_plus_FMDsurgery"]]),
#   run_scenario("10_Vosoritide_PoorAdherence_60pct", scenarios[["10_Vosoritide_PoorAdherence_60pct"]], adherence = 0.6)
# )
#
# ggplot(results, aes(time/24, HEIGHTZ, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Day", y = "Height Z-score", title = "Achondroplasia: Height Z-score trajectory by scenario")

# =====================================================================
# Calibration notes (unchanged from the original — no numeric value
# touched by this refactor):
#  - Untreated AGV_BASE = 3.9 cm/yr reproduces the placebo-arm growth
#    velocity in Savarirayan 2020 Lancet (ages 5-14, prepubertal ACH).
#  - GAIN_AGV / EMAX_PERKINH calibrated so that 15 µg/kg QD vosoritide
#    yields ΔAGV ≈ +1.57 cm/yr at steady state (matches the Lancet
#    phase-3 treatment difference at week 52).
#  - Dose-ranging (2.5/7.5/15/30 µg/kg) reproduces the Savarirayan 2019
#    NEJM phase-2 plateau: efficacy rises steeply to 15 µg/kg then
#    plateaus (30 µg/kg no incremental AGV benefit, calibrated via the
#    saturating CNP_TOTAL Hill term, EFFECT_VOS, was CGMP_DRIVE).
#  - Infigratinib EMAX_INFIG/K_OFFTARGET_FGFR1 (was EMAX_INFIG_PERK)
#    reflect the PROPEL 2/3 direct FGFR3 kinase-blockade efficacy signal,
#    with a mild off-target FGFR1 phosphate liability (PROPEL 3: ~4%
#    mild/transient asymptomatic hyperphosphatemia, no discontinuations,
#    no ocular FGFR1/2 AEs).
#  - TransCon CNP release/exposure (KA_TCNP, CL_TCNP; was KREL_TCNP,
#    KE_TCNP) approximate the once-weekly navepegritide profile behind the
#    FDA-approved (Feb-2026, YUVIWEL) ApproaCH phase-3 ΔAGV of +1.49 cm/yr
#    at week 52.
#  - FMAREA/SPCANALZ use slow first-order relaxation (K_FM, K_SPCANAL)
#    reflecting that skeletal/structural remodeling lags biochemical
#    pathway rescue by years, not weeks (Hunter 1998 natural history).
# =====================================================================
