## ============================================================
## Osteoarthritis (OA) QSP Model — mrgsolve Implementation
## ============================================================
## Disease: Primary knee osteoarthritis (KOA)
## Scope   : Chondrocyte biology · ECM catabolism · Synovitis ·
##           Subchondral bone remodelling · Pain signalling ·
##           Drug PK/PD (NSAID · IA-CS · HA · Sprifermin · Tanezumab)
## Version : 1.0  (2026-06-23)
##
## Key References (parameter calibration):
##   Thijssen E et al. Osteoarthritis Cartilage 2015 (MMP-13 dynamics)
##   Hunter DJ et al. NEJM 2019 (OA pathophysiology review)
##   Lohmander LS et al. Ann Rheum Dis 2005 (uCTX-II biomarker)
##   Berenbaum F et al. Nat Rev Rheumatol 2013 (OA inflammation)
##   Hochberg MC et al. Osteoarthritis Cartilage 2017 (JSW loss)
##   Kloppenburg M et al. Ann Rheum Dis 2020 (OARSI guidelines)
##   Dahlberg L et al. Lancet 2016 (sprifermin FORWARD RCT)
##   Schnitzer TJ et al. NEJM 2019 (tanezumab Phase-3)
##   Conaghan PG et al. Arthritis Rheum 2016 (celecoxib PK/PD)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─────────────────────────────────────────────────────────────
## 1. MODEL CODE STRING
## ─────────────────────────────────────────────────────────────

oa_code <- '
$PROB
Osteoarthritis QSP Model
- 20 ODE compartments
- Drug PK: NSAID (celecoxib), IA corticosteroid, HA, Sprifermin, Tanezumab
- Disease PD: IL-1b, MMP-13, ADAMTS-5, Collagen II, Aggrecan,
              Chondrocyte, Synovitis, Osteoclast, JSW, PGE2, Pain, uCTX-II, COMP

$PARAM
// ── Baseline/Normal States ──────────────────────────────────
// [PK/PD refactor build-compat] ColII_0/Chondro_0/JSW_0 renamed to
// ColII_BASE/Chondro_BASE/JSW_BASE (values unchanged) -- mrgsolve 2.0.1
// auto-generates a "<compartment>_0" init symbol for every $INIT
// compartment, so a $PARAM entry of that exact name for a same-named
// compartment (ColII, Chondro, JSW all exist as $INIT compartments) is a
// hard compile-time redeclaration. Not NSAID/TANZ-scoped; disclosed in
// oa_refactor_notes.md and logged in translations/UPSTREAM_ISSUES.md.
ColII_BASE   = 100     // Collagen II baseline (% normal)
Agg_0        = 100     // Aggrecan baseline (% normal)
Chondro_BASE = 100     // Chondrocyte count baseline (% normal)
JSW_BASE     = 5.5     // Joint space width baseline (mm), typical KOA at KL2
OC_0         = 1.0     // Osteoclast activity baseline (normalized)
OB_0         = 1.0     // Osteoblast activity baseline (normalized)
Syn_0        = 0.5     // Synovitis score baseline (0-3)

// ── Disease Progression Parameters ──────────────────────────
k_IL1b_syn   = 0.10    // IL-1b synthesis rate (pg/mL/day) from synovitis
k_IL1b_deg   = 0.35    // IL-1b elimination rate (/day), t1/2 ~2d
k_TNFa_syn   = 0.08    // TNF-a synthesis (pg/mL/day)
k_TNFa_deg   = 0.40    // TNF-a elimination (/day)
k_MMP13_syn  = 0.15    // MMP-13 synthesis (ng/mL/day), IL1b/TNFa driven
k_MMP13_deg  = 0.30    // MMP-13 elimination (/day)
k_ADAM5_syn  = 0.12    // ADAMTS-5 synthesis (/day)
k_ADAM5_deg  = 0.28    // ADAMTS-5 elimination (/day)
K_IL1_mmp    = 2.0     // EC50 for IL-1b→MMP-13 (pg/mL)
Emax_IL1_mmp = 3.0     // Emax for IL-1b→MMP-13 fold induction
K_IL1_adam5  = 2.5     // EC50 for IL-1b→ADAMTS-5 (pg/mL)
k_ColII_syn  = 0.008   // Collagen II synthesis (/day), Chondro-dependent
k_ColII_deg  = 0.12    // Collagen II MMP-13-dependent degradation
K_mmp_col    = 2.0     // EC50 for MMP-13→ColII deg (ng/mL)
k_Agg_syn    = 0.010   // Aggrecan synthesis (/day)
k_Agg_deg    = 0.15    // Aggrecan ADAMTS-5-dependent deg (/day per normalized unit)
K_adam5_agg  = 1.5     // EC50 for ADAMTS-5→Aggrecan (normalized)
k_Chondro_prof= 0.003  // Chondrocyte proliferation rate (/day)
k_Chondro_apt= 0.006   // Chondrocyte apoptosis rate (/day) baseline + IL1b-driven
Emax_IL1_apt = 2.0     // Emax for IL-1b-driven apoptosis
K_IL1_apt    = 5.0     // EC50 for IL-1b→apoptosis (pg/mL)
k_Syn_syn    = 0.05    // Synovitis genesis (/day), ECM-loss driven
k_Syn_deg    = 0.10    // Synovitis resolution (/day)
Emax_ecm_syn = 1.5     // ECM loss → synovitis amplification
K_ecm_syn    = 30.0    // EC50 for ECM loss → synovitis (% loss)
k_OC_form    = 0.08    // Osteoclast formation rate (RANKL/OPG driven)
k_OC_deg     = 0.15    // Osteoclast apoptosis rate (/day)
Emax_IL1_OC  = 1.8     // IL-1b-driven RANKL → OC increase Emax
K_IL1_OC     = 3.0     // EC50 (pg/mL)
k_OB_form    = 0.06    // Osteoblast formation rate
k_OB_deg     = 0.10    // Osteoblast elimination
k_JSW_loss   = 0.0012  // JSW loss rate (mm/day), MMP-13 and ADAMTS-5 driven
k_JSW_repair = 0.0001  // Very slow JSW passive repair
k_PGE2_syn   = 0.20    // PGE2 synthesis (ng/mL/day), COX-2/IL-1b driven
k_PGE2_deg   = 0.60    // PGE2 elimination (/day), t1/2~28h
k_Pain_syn   = 0.08    // Pain generation rate (VAS units/day)
k_Pain_deg   = 0.12    // Pain resolution rate (/day)
Pain_baseline= 10.0    // Residual/structural pain floor (VAS)
Emax_PGE2_pain=60.0    // Emax PGE2→pain contribution (VAS units)
K_PGE2_pain  = 0.5     // EC50 PGE2→pain (ng/mL)
Emax_Syn_pain= 20.0    // Synovitis → pain Emax
K_Syn_pain   = 1.0     // EC50 synovitis
k_CTXII_syn  = 0.05    // uCTX-II synthesis (nmol/mmolCr/day), ColII deg product
k_CTXII_deg  = 0.25    // uCTX-II elimination (/day)
k_COMP_rel   = 0.04    // COMP release rate (μg/mL/day)
k_COMP_deg   = 0.20    // COMP elimination (/day)

// ── NSAID PK/PD (Celecoxib 200 mg) ──────────────────────────
// [PK/PD refactor] renamed to the guide's naming convention (values
// unchanged). Archetype 3 minus peripheral-via-Q/V2 (depot+central+a
// second, partition-linked tissue compartment) -- see refactor notes.
KA_NSAID     = 1.20    // Absorption rate constant (/h)
F_NSAID      = 0.74    // Oral bioavailability (fraction)
V1_NSAID     = 455     // Volume of distribution (L), population mean
CL_NSAID     = 27.7    // Clearance (L/h), t1/2≈11.5h
KP_NSAID     = 0.40    // Plasma-to-joint partition (synovial fluid:plasma ~0.4) -- bespoke, no convention slot
EC50_NSAID   = 0.042   // EC50 celecoxib for COX-2 (μg/mL) (was IC50_cox2)
EMAX_NSAID   = 0.95    // Max COX-2 inhibition fraction (was Emax_cox2)
GAMMA_NSAID  = 1.2     // Hill coefficient (was hill_cox2)
// Naproxen alternative: ka=1.1/h, F=0.99, Vd=11L, CL=0.55L/h, IC50=0.28μg/mL

// ── IA Corticosteroid PK/PD (Triamcinolone acetonide 40 mg) ─
k_IACS_abs   = 0.008   // Absorption from joint (slow, t1/2≈3.5d) (/h)
k_IACS_el    = 0.25    // Plasma elimination (/h), t1/2~2.7h
Vd_IACS      = 99      // Plasma Vd (L/70kg)
EC50_GR      = 0.05    // EC50 for GR-mediated NF-kB inhibition (μg/mL joint)
Emax_GR_NF   = 0.75    // Max NF-kB inhibition by GCS

// ── IA Hyaluronic Acid PK ────────────────────────────────────
k_HA_deg     = 0.041   // HA degradation in joint (/h), t1/2≈17h
Emax_HA_pain = 0.35    // HA viscosupplement max pain reduction fraction
K_HA_pain    = 0.8     // EC50 HA joint conc (mg/mL)
Emax_HA_IL1  = 0.25    // HA → IL-1b reduction Emax

// ── Sprifermin (FGF-18) PK/PD ───────────────────────────────
k_Sprif_deg  = 0.044   // Elimination from joint (/h), t1/2≈16 min converted to /h = ln2/(16/60h)
Emax_FGF_col = 0.40    // Sprifermin max ColII synthesis increase fraction (FORWARD RCT)
K_FGF18      = 0.002   // EC50 (μg/mL joint)

// ── Tanezumab (anti-NGF mAb) PK/PD ──────────────────────────
// [PK/PD refactor] renamed to the guide's naming convention (values
// unchanged). Checked for TMDD: this file models tanezumab with plain
// linear depot+central PK and a direct Emax*C/(EC50+C) pain-reduction
// ratio -- no receptor-binding ODE state anywhere, so this is archetype 3
// minus peripheral (depot+central, linear), not archetype 4/TMDD, despite
// tanezumab being an anti-NGF mAb in reality. GAMMA_TANZ=1: the original
// had no explicit Hill exponent (a plain first-order ratio).
KA_TANZ      = 0.0167  // SC absorption (/h), t1/2_abs≈41h → F=0.73
F_TANZ       = 0.73    // SC bioavailability
V1_TANZ      = 7.8     // Vd (L/70kg) for IgG
CL_TANZ      = 0.0133  // CL (/h), t1/2≈23d
EMAX_TANZ    = 45.0    // Max pain reduction (VAS units, tanezumab anti-NGF) (was Emax_tanz_pain)
EC50_TANZ    = 0.001   // EC50 for tanezumab plasma→NGF binding (was K_tanz_NGF)
GAMMA_TANZ   = 1       // Hill coefficient -- original had no explicit exponent

// ── KL Grade / Disease Severity Modifier ─────────────────────
KL_grade     = 2       // Kellgren-Lawrence grade at baseline (1-4)
// KL modifies disease load: KL_mult = 0.5 + 0.25*KL_grade
KL_mult      = 1.0     // Runtime calculated from KL_grade
// Age effect
Age_yr       = 62      // Patient age (years)
k_age_mod    = 0.005   // Age-dependent disease rate modifier (/yr above 45)

// ── Input Flags ──────────────────────────────────────────────
NSAID_flag   = 0       // 1 = NSAID dosing active
IACS_flag    = 0       // 1 = IA corticosteroid given
HA_flag      = 0       // 1 = IA HA viscosupplementation
Sprif_flag   = 0       // 1 = Sprifermin dosing
TANZ_flag    = 0       // 1 = Tanezumab dosing (was Tanz_flag)

$INIT
// ── Initial States ───────────────────────────────────────────
// Drug PK compartments (all zero at t=0)
// [PK/PD refactor] NSAID/TANZ compartments renamed to the guide's
// convention (GUT_/CENT_/PERI_<STEM>); IACS/HA/Sprifermin untouched.
GUT_NSAID      = 0
CENT_NSAID     = 0
PERI_NSAID     = 0
A_IACS_joint   = 0
A_IACS_plasma  = 0
A_HA_joint     = 0
A_Sprif_joint  = 0
GUT_TANZ       = 0
CENT_TANZ      = 0

// ── Disease State Initial Conditions ─────────────────────────
IL1b    = 5.0     // pg/mL (elevated vs normal ~1 pg/mL in mild KOA)
TNFa    = 3.0     // pg/mL
MMP13   = 8.0     // ng/mL (elevated in OA synovial fluid)
ADAM5   = 2.0     // normalized ADAMTS-5 activity
ColII   = 70.0    // % normal (KL2 = moderate loss)
Aggrecan= 65.0    // % normal
Chondro = 75.0    // % normal chondrocyte count
Synovitis=1.0     // 0–3 scale (mild-moderate)
OC_act  = 1.3     // osteoclast activity (slightly elevated)
OB_act  = 1.1     // osteoblast activity
JSW     = 3.8     // mm (mild joint space narrowing from 5.5 baseline)
PGE2_jt = 1.5     // ng/mL joint PGE2
VASPain = 45.0    // VAS pain 0-100 (moderate)
uCTXII  = 0.40    // nmol/mmolCr urine CTX-II
COMP_s  = 12.0    // μg/mL serum COMP

$ODE
// ─────────────────────────────────────────────────────────────
// Helper variables
// ─────────────────────────────────────────────────────────────
double age_mod  = 1.0 + k_age_mod * (Age_yr - 45.0);
double KL_mod   = 0.5 + 0.25 * KL_grade;

// [PK/PD refactor] NSAID concentration -- normalized to one cached site per
// tissue (was computed a second time, independently, in $TABLE under the
// names C_nsaid_pl/C_nsaid_jt; that $TABLE copy is now redirected to read
// these renamed compartments instead of re-deriving from old names -- see
// $TABLE below). CP_NSAID (plasma) is informational only, no PD term reads
// it, exactly as in the original (C_NSAID_plasma was computed but never
// used downstream there either). C_NSAID (joint) is the single PD-facing
// exposed concentration -- the only one COX-2 inhibition reads.
double CP_NSAID = CENT_NSAID / V1_NSAID;                 // plasma conc (μg/mL) -- informational only
double C_NSAID  = PERI_NSAID / (V1_NSAID * 0.01);        // joint conc (μg/mL), ~1% of total in joint -- PD-facing

// COX-2 inhibition fraction (Hill equation). EFFECT_NSAID renamed from
// COX2_inh -- pure rename, already this exact Emax/EC50/Hill shape.
double EFFECT_NSAID = EMAX_NSAID * pow(C_NSAID, GAMMA_NSAID) /
                  (pow(EC50_NSAID, GAMMA_NSAID) + pow(C_NSAID, GAMMA_NSAID));
double COX2_eff = 1.0 - EFFECT_NSAID * NSAID_flag;   // residual COX-2 activity

// IA corticosteroid joint concentration (μg/mL)
double C_IACS = A_IACS_joint;
// GR-mediated NF-kB inhibition
double GR_inh = (Emax_GR_NF * C_IACS / (EC50_GR + C_IACS)) * IACS_flag;
double NF_kB_eff = 1.0 - GR_inh;   // multiplier on NF-kB driven processes

// HA joint concentration (mg/mL)
double C_HA  = A_HA_joint;
double HA_pain_inh = Emax_HA_pain  * C_HA / (K_HA_pain  + C_HA) * HA_flag;
double HA_IL1_inh  = Emax_HA_IL1   * C_HA / (K_HA_pain  + C_HA) * HA_flag;

// Sprifermin joint concentration (μg/mL)
double C_Sprif = A_Sprif_joint;
double FGF18_eff = Emax_FGF_col * C_Sprif / (K_FGF18 + C_Sprif) * Sprif_flag;

// [PK/PD refactor] Tanezumab plasma concentration -- normalized to one
// cached site (was computed a second time, independently, in $TABLE under
// the name C_Tanz_pl; that $TABLE copy is now redirected to read this
// renamed compartment -- see $TABLE below). No separate tissue site for
// tanezumab in this file, so a single C_TANZ is exposed.
double C_TANZ = CENT_TANZ / V1_TANZ;
// EFFECT_TANZ renamed from Tanz_pain_red -- tanezumab sequesters NGF ->
// proportional pain reduction (VAS units). Pure rename, already
// Emax*C/(EC50+C) shape (GAMMA_TANZ=1, original had no Hill exponent).
double EFFECT_TANZ = EMAX_TANZ * pow(C_TANZ, GAMMA_TANZ) / (pow(EC50_TANZ, GAMMA_TANZ) + pow(C_TANZ, GAMMA_TANZ)) * TANZ_flag;

// ECM loss (0–100, where 0=normal, 100=complete loss)
double ECM_loss = 100.0 - 0.5*(ColII + Aggrecan);
double ECM_loss_pos = (ECM_loss < 0) ? 0.0 : ECM_loss;

// ─────────────────────────────────────────────────────────────
// DRUG PK ODEs
// ─────────────────────────────────────────────────────────────
// NSAID
dxdt_GUT_NSAID    = -KA_NSAID * GUT_NSAID;
dxdt_CENT_NSAID   = F_NSAID * KA_NSAID * GUT_NSAID * NSAID_flag
                      - (CL_NSAID / V1_NSAID) * CENT_NSAID
                      - KP_NSAID * (CL_NSAID / V1_NSAID) * CENT_NSAID;
dxdt_PERI_NSAID   = KP_NSAID * (CL_NSAID / V1_NSAID) * CENT_NSAID
                      - (CL_NSAID / V1_NSAID) * PERI_NSAID;

// IA Corticosteroid
dxdt_A_IACS_joint   = -k_IACS_abs * A_IACS_joint;
dxdt_A_IACS_plasma  = k_IACS_abs * A_IACS_joint - k_IACS_el * A_IACS_plasma;

// IA Hyaluronic Acid
dxdt_A_HA_joint     = -k_HA_deg * A_HA_joint;

// Sprifermin (IA FGF-18)
dxdt_A_Sprif_joint  = -k_Sprif_deg * A_Sprif_joint;

// Tanezumab (SC mAb)
dxdt_GUT_TANZ    = -KA_TANZ * GUT_TANZ;
dxdt_CENT_TANZ   = F_TANZ * KA_TANZ * GUT_TANZ - (CL_TANZ / V1_TANZ) * CENT_TANZ;

// ─────────────────────────────────────────────────────────────
// DISEASE PD ODEs
// ─────────────────────────────────────────────────────────────

// IL-1β (pg/mL) — master catabolic cytokine
// Synthesis driven by synovitis, age/KL severity; inhibited by GCS and HA
double IL1b_pos = (IL1b < 0.01) ? 0.01 : IL1b;
double IL1b_syn_rate = k_IL1b_syn * Synovitis * age_mod * KL_mod
                       * NF_kB_eff * (1.0 - HA_IL1_inh);
dxdt_IL1b = IL1b_syn_rate - k_IL1b_deg * IL1b;

// TNF-α (pg/mL) — co-driver with IL-1b
double TNFa_syn_rate = k_TNFa_syn * Synovitis * NF_kB_eff;
dxdt_TNFa = TNFa_syn_rate - k_TNFa_deg * TNFa;

// MMP-13 (ng/mL) — dominant collagenase in OA cartilage
// Induced by IL-1b, TNFa; inhibited by TIMP/GCS
double IL1b_TNFa_stim = (IL1b / (K_IL1_mmp + IL1b)) * Emax_IL1_mmp
                        + TNFa / (1.0 + TNFa);
double MMP13_syn_rate = k_MMP13_syn * IL1b_TNFa_stim * NF_kB_eff;
dxdt_MMP13 = MMP13_syn_rate - k_MMP13_deg * MMP13;

// ADAMTS-5 (normalized activity) — dominant aggrecanase
double ADAM5_syn_rate = k_ADAM5_syn * (IL1b / (K_IL1_adam5 + IL1b))
                        * NF_kB_eff;
dxdt_ADAM5 = ADAM5_syn_rate - k_ADAM5_deg * ADAM5;

// Collagen Type II (% normal, 0–100)
// Synthesis: chondrocyte-dependent, FGF18 anabolic boost
// Degradation: MMP-13 dependent (Hill), floor at 0
double ColII_pos = (ColII < 0.5) ? 0.5 : ColII;
double col_syn = k_ColII_syn * (Chondro / Chondro_BASE) * (1.0 + FGF18_eff);
double col_deg = k_ColII_deg * MMP13 / (K_mmp_col + MMP13) * ColII_pos * age_mod;
dxdt_ColII = col_syn - col_deg;

// Aggrecan (% normal, 0–100)
double Agg_pos = (Aggrecan < 0.5) ? 0.5 : Aggrecan;
double agg_syn = k_Agg_syn * (Chondro / Chondro_BASE) * (1.0 + 0.5*FGF18_eff);
double agg_deg = k_Agg_deg * ADAM5 / (K_adam5_agg + ADAM5) * Agg_pos * age_mod;
dxdt_Aggrecan = agg_syn - agg_deg;

// Chondrocyte population (% normal)
double Chondro_pos = (Chondro < 1.0) ? 1.0 : Chondro;
double chondro_prolif = k_Chondro_prof * Chondro_pos;
double chondro_apopt  = (k_Chondro_apt + k_Chondro_apt * Emax_IL1_apt *
                         IL1b / (K_IL1_apt + IL1b)) * Chondro_pos * age_mod;
dxdt_Chondro = chondro_prolif - chondro_apopt;

// Synovitis score (0–3)
// Driven by ECM damage (FN fragments, DAMP release), amplified by macrophage
// Resolved by anti-inflammatory mediators; GCS reduces
double Syn_drive = k_Syn_syn * (1.0 + Emax_ecm_syn * ECM_loss_pos /
                                (K_ecm_syn + ECM_loss_pos)) * age_mod;
double Syn_res   = k_Syn_deg * Synovitis * NF_kB_eff;
dxdt_Synovitis = Syn_drive - Syn_res;

// Osteoclast activity (normalized)
// RANKL/OPG imbalance driven by IL-1b, TNFa
double OC_drive = k_OC_form * (1.0 + Emax_IL1_OC * IL1b / (K_IL1_OC + IL1b));
dxdt_OC_act = OC_drive - k_OC_deg * OC_act;

// Osteoblast activity (normalized)
// Slightly activated by TGF-b from OC; GCS suppresses new bone
double OB_drive = k_OB_form * (1.0 + 0.2 * OC_act) * NF_kB_eff;
dxdt_OB_act = OB_drive - k_OB_deg * OB_act;

// Joint Space Width (mm)
// Loss driven by ECM degradation and OC activity; very slow passive
double JSW_pos = (JSW < 0.5) ? 0.5 : JSW;
double jsw_loss = k_JSW_loss * (MMP13 / (K_mmp_col + MMP13) +
                                ADAM5 / (K_adam5_agg + ADAM5)) * JSW_pos * age_mod;
double jsw_gain = k_JSW_repair * (1.0 + FGF18_eff);
dxdt_JSW = jsw_gain - jsw_loss;

// PGE2 in joint (ng/mL)
// Driven by COX-2 activity (which is inhibited by NSAID)
double PGE2_syn  = k_PGE2_syn * IL1b / (2.0 + IL1b) * COX2_eff;
dxdt_PGE2_jt = PGE2_syn - k_PGE2_deg * PGE2_jt;

// VAS Pain Score (0–100)
// PGE2 (inflammatory pain), synovitis, structural pain, central sensitization
// Reduced by NSAID (via PGE2), HA (viscosupplementation), tanezumab (NGF)
double PGE2_pain = Emax_PGE2_pain * PGE2_jt / (K_PGE2_pain + PGE2_jt);
double Syn_pain  = Emax_Syn_pain  * Synovitis / (K_Syn_pain + Synovitis);
double Struct_pain = Pain_baseline * (1.0 - JSW / JSW_BASE); // structural pain
double Pain_target = PGE2_pain + Syn_pain + Struct_pain
                     - HA_pain_inh * 20.0
                     - EFFECT_TANZ;
double Pain_target_clamp = (Pain_target < 0.0) ? 0.0 :
                           (Pain_target > 100.0) ? 100.0 : Pain_target;
dxdt_VASPain = k_Pain_syn * (Pain_target_clamp - VASPain);

// uCTX-II (nmol/mmolCr) — urinary type II collagen degradation biomarker
// Reflects active collagen catabolism; released when MMP-13 cleaves ColII
double CTX2_syn  = k_CTXII_syn * col_deg * 10.0;  // scaled to nM range
dxdt_uCTXII = CTX2_syn - k_CTXII_deg * uCTXII;

// Serum COMP (μg/mL) — reflects cartilage stress/damage
double COMP_drive = k_COMP_rel * (1.0 + Synovitis + MMP13 / 5.0);
dxdt_COMP_s = COMP_drive - k_COMP_deg * COMP_s;

$TABLE
// ── Derived Outputs ─────────────────────────────────────────
// [PK/PD refactor] reads $ODE's EFFECT_NSAID directly (cross-block read of
// an already-declared $ODE local, not a redeclaration) -- the exact
// pattern the original already used for this line (COX2_inh, computed in
// $ODE, read here unchanged); only the identifier is renamed.
double COX2_pct_inh = EFFECT_NSAID * NSAID_flag * 100.0;

// KOOS (0–100 scale, higher = better)
// Components: pain, symptoms, ADL, sport, QoL
// Approximate from VAS pain and JSW
double KOOS_pain = 100.0 - VASPain * 0.60;
double KOOS_func = 50.0 + (JSW / JSW_BASE) * 30.0 - (Synovitis / 3.0) * 10.0;
double KOOS_est  = 0.6 * KOOS_pain + 0.4 * KOOS_func;
KOOS_est = (KOOS_est < 0) ? 0 : (KOOS_est > 100) ? 100 : KOOS_est;

// WOMAC pain subscale (0–20, lower = better)
double WOMAC_pain_sub = VASPain / 5.0;

// Cartilage volume estimate (% baseline)
double CartVol_pct = 0.5 * (ColII / ColII_BASE) + 0.5 * (Aggrecan / Agg_0);

// ECM degradation score
double ECM_degradation = 100.0 - CartVol_pct * 100.0;

// [PK/PD refactor] $TABLE compiles into the same scope as $ODE in
// mrgsolve (reusing an $ODE-declared name here is a hard "redefinition"
// compile error, confirmed empirically -- see refactor notes), so these
// three lines must keep their own, $TABLE-local names (distinct from
// $ODE's C_NSAID/CP_NSAID/C_TANZ). Normalized: redirected their RHS to
// the renamed canonical compartments/parameters instead of the old
// A_NSAID_plasma/A_NSAID_joint/A_Tanz_plasma names -- this is now the
// single site downstream code should read for the $TABLE-scope value.
double C_nsaid_pl = CENT_NSAID / V1_NSAID;
double C_nsaid_jt = PERI_NSAID / (V1_NSAID * 0.01);
double C_IACS_jt  = A_IACS_joint;
double C_HA_jt    = A_HA_joint;
double C_Sprif_jt = A_Sprif_joint;
double C_Tanz_pl  = CENT_TANZ / V1_TANZ;

// GR inhibition percentage
double GR_inhibition_pct = ((Emax_GR_NF * C_IACS_jt / (EC50_GR + C_IACS_jt))
                             * IACS_flag) * 100.0;

$CAPTURE
// [PK/PD refactor build-compat] the 15 disease-state names removed here
// (IL1b, TNFa, MMP13, ADAM5, ColII, Aggrecan, Chondro, Synovitis, OC_act,
// OB_act, JSW, PGE2_jt, VASPain, uCTXII, COMP_s) are $INIT compartments,
// already emitted automatically as output columns; mrgsolve 2.0.1 rejects
// re-listing a compartment in $CAPTURE. Not NSAID/TANZ-scoped; disclosed
// in oa_refactor_notes.md and logged in translations/UPSTREAM_ISSUES.md.
// All 15 are still present in every simulation output under their same
// names (mrgsolve emits every compartment state unconditionally).
COX2_pct_inh C_nsaid_pl C_nsaid_jt C_IACS_jt C_HA_jt C_Sprif_jt C_Tanz_pl
GR_inhibition_pct KOOS_est WOMAC_pain_sub CartVol_pct ECM_degradation
'

## ─────────────────────────────────────────────────────────────
## 2. COMPILE MODEL
## ─────────────────────────────────────────────────────────────
mod <- mcode("oa_qsp_refactored", oa_code)

## ─────────────────────────────────────────────────────────────
## 3. DOSING REGIMENS
## ─────────────────────────────────────────────────────────────

# Simulation period: 2 years (730 days)
sim_end <- 730

# NSAID: Celecoxib 200 mg BID (every 12h), continuous
# Dose = 200 mg; Vd = 455 L → 200 mg dose input
nsaid_regimen <- ev(amt = 200, ii = 12, addl = sim_end * 2, cmt = "GUT_NSAID")

# IA Corticosteroid: Triamcinolone acetonide 40 mg Q3 months (4 injections/year)
# 40 mg dose into IACS joint compartment
iacs_regimen <- ev(time = c(0, 91, 182, 273, 365, 456, 547, 638),
                   amt = 40, cmt = "A_IACS_joint")

# IA Hyaluronic Acid: Synvisc 16 mg × 3 weekly (2 courses, at 0 and 6m)
ha_regimen <- ev(time = c(0, 7, 14, 182, 189, 196),
                 amt = 16, cmt = "A_HA_joint")

# Sprifermin (FGF-18): 30 μg IA every 12 weeks × 5 (FORWARD RCT inspired)
# Dose in μg directly into joint
sprif_regimen <- ev(time = c(0, 84, 168, 252, 336, 420, 504, 588, 672),
                    amt = 0.030, cmt = "A_Sprif_joint")

# Tanezumab: 2.5 mg SC every 8 weeks (anti-NGF mAb, Phase-3 inspired)
tanz_regimen <- ev(time = c(0, 56, 112, 168, 224, 280, 336, 392, 448,
                             504, 560, 616, 672, 728),
                   amt = 2.5, cmt = "GUT_TANZ")

## ─────────────────────────────────────────────────────────────
## 4. TREATMENT SCENARIOS
## ─────────────────────────────────────────────────────────────

run_scenario <- function(mod, regimen, params_override = list(), label = "Scenario") {
  p_mod <- do.call(param, c(list(mod), params_override))
  out   <- mrgsim(p_mod, ev = regimen, end = sim_end, delta = 1,
                  carry_out = "evid,amt,cmt")
  as_tibble(out) %>% mutate(Scenario = label)
}

# Scenario 1: Disease progression, no treatment
scen1 <- run_scenario(mod, ev(),
  params_override = list(NSAID_flag=0, IACS_flag=0, HA_flag=0, Sprif_flag=0, TANZ_flag=0),
  label = "1. No Treatment (Natural History)")

# Scenario 2: NSAID continuous (Celecoxib 200 mg BID)
scen2 <- run_scenario(mod, nsaid_regimen,
  params_override = list(NSAID_flag=1, IACS_flag=0, HA_flag=0, Sprif_flag=0, TANZ_flag=0),
  label = "2. Celecoxib 200mg BID")

# Scenario 3: IA Corticosteroid (Triamcinolone 40mg Q3mo)
scen3 <- run_scenario(mod, iacs_regimen,
  params_override = list(NSAID_flag=0, IACS_flag=1, HA_flag=0, Sprif_flag=0, TANZ_flag=0),
  label = "3. IA Triamcinolone 40mg Q3mo")

# Scenario 4: IA Hyaluronic Acid (viscosupplementation)
scen4 <- run_scenario(mod, ha_regimen,
  params_override = list(NSAID_flag=0, IACS_flag=0, HA_flag=1, Sprif_flag=0, TANZ_flag=0),
  label = "4. IA Hyaluronic Acid 3x")

# Scenario 5: Sprifermin (DMOAD, cartilage structural modification)
scen5 <- run_scenario(mod, sprif_regimen,
  params_override = list(NSAID_flag=0, IACS_flag=0, HA_flag=0, Sprif_flag=1, TANZ_flag=0),
  label = "5. Sprifermin 30μg IA Q12w")

# Scenario 6: Tanezumab (anti-NGF, pain-focused)
scen6 <- run_scenario(mod, tanz_regimen,
  params_override = list(NSAID_flag=0, IACS_flag=0, HA_flag=0, Sprif_flag=0, TANZ_flag=1),
  label = "6. Tanezumab 2.5mg SC Q8w")

# Scenario 7: Combination — NSAID + IA Corticosteroid
combo_regimen7 <- ev(nsaid_regimen, iacs_regimen)
scen7 <- run_scenario(mod, combo_regimen7,
  params_override = list(NSAID_flag=1, IACS_flag=1, HA_flag=0, Sprif_flag=0, TANZ_flag=0),
  label = "7. Celecoxib + IA Triamcinolone")

# Scenario 8: Sprifermin + NSAID (structural + symptomatic)
combo_regimen8 <- ev(nsaid_regimen, sprif_regimen)
scen8 <- run_scenario(mod, combo_regimen8,
  params_override = list(NSAID_flag=1, IACS_flag=0, HA_flag=0, Sprif_flag=1, TANZ_flag=0),
  label = "8. Sprifermin + Celecoxib")

# Combine all scenarios
all_scen <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6, scen7, scen8)

## ─────────────────────────────────────────────────────────────
## 5. SENSITIVITY ANALYSIS — Parameter Uncertainty
## ─────────────────────────────────────────────────────────────

# Vary key parameters ±30% for no-treatment scenario
sens_params <- c("k_MMP13_syn", "k_ColII_deg", "k_Agg_deg", "k_JSW_loss", "k_Pain_syn")
sens_results <- list()

for (p in sens_params) {
  for (mult in c(0.7, 1.0, 1.3)) {
    pmod <- param(mod, setNames(list(param(mod)[[p]] * mult), p))
    out  <- mrgsim(pmod, ev = ev(),
                   param = list(NSAID_flag=0, IACS_flag=0, HA_flag=0, Sprif_flag=0, TANZ_flag=0),
                   end = sim_end, delta = 7) %>%
            as_tibble() %>%
            mutate(SensParam = p, SensMult = as.character(mult))
    sens_results[[paste0(p, "_", mult)]] <- out
  }
}
sens_df <- bind_rows(sens_results)

## ─────────────────────────────────────────────────────────────
## 6. PLOTS
## ─────────────────────────────────────────────────────────────

theme_oa <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "#E3F2FD"),
        plot.title = element_text(face = "bold", size = 12))

scen_colors <- c("1. No Treatment (Natural History)" = "#B71C1C",
                 "2. Celecoxib 200mg BID"            = "#1565C0",
                 "3. IA Triamcinolone 40mg Q3mo"     = "#2E7D32",
                 "4. IA Hyaluronic Acid 3x"          = "#6A1B9A",
                 "5. Sprifermin 30μg IA Q12w"        = "#E65100",
                 "6. Tanezumab 2.5mg SC Q8w"        = "#00838F",
                 "7. Celecoxib + IA Triamcinolone"   = "#F9A825",
                 "8. Sprifermin + Celecoxib"         = "#AD1457")

# Plot 1: VAS Pain over time
p1 <- ggplot(all_scen, aes(time, VASPain, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "OA QSP Model: VAS Pain Score over 2 Years",
       x = "Time (days)", y = "VAS Pain (0–100)",
       color = "Treatment Scenario") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_oa
print(p1)

# Plot 2: Joint Space Width (structural outcome)
p2 <- ggplot(all_scen, aes(time, JSW, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Joint Space Width (mm) — Structural Outcome",
       x = "Time (days)", y = "JSW (mm)") +
  geom_hline(yintercept = 2.0, linetype = "dashed", color = "gray50",
             linewidth = 0.5) +
  annotate("text", x = 700, y = 2.1, label = "TKR threshold ~2mm",
           size = 3, color = "gray50") +
  theme_oa
print(p2)

# Plot 3: Collagen II and Aggrecan
p3 <- all_scen %>%
  pivot_longer(c(ColII, Aggrecan), names_to = "ECM", values_to = "Value") %>%
  ggplot(aes(time, Value, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ECM, scales = "fixed",
             labeller = labeller(ECM = c(ColII = "Collagen II (%)", Aggrecan = "Aggrecan (%)"))) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Cartilage ECM Components over Time",
       x = "Time (days)", y = "% Baseline") +
  scale_y_continuous(limits = c(0, 110)) +
  theme_oa
print(p3)

# Plot 4: Inflammatory Biomarkers
p4 <- all_scen %>%
  pivot_longer(c(IL1b, TNFa, MMP13), names_to = "Marker", values_to = "Value") %>%
  ggplot(aes(time, Value, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Marker, scales = "free_y",
             labeller = labeller(Marker = c(IL1b = "IL-1β (pg/mL)",
                                            TNFa = "TNF-α (pg/mL)",
                                            MMP13 = "MMP-13 (ng/mL)"))) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Inflammatory Mediators over Time",
       x = "Time (days)", y = "Concentration") +
  theme_oa
print(p4)

# Plot 5: Circulating Biomarkers (uCTX-II, COMP)
p5 <- all_scen %>%
  pivot_longer(c(uCTXII, COMP_s), names_to = "Biomarker", values_to = "Value") %>%
  ggplot(aes(time, Value, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Biomarker, scales = "free_y",
             labeller = labeller(Biomarker = c(uCTXII = "uCTX-II (nmol/mmolCr)",
                                               COMP_s = "Serum COMP (μg/mL)"))) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Circulating OA Biomarkers",
       x = "Time (days)", y = "Concentration") +
  theme_oa
print(p5)

# Plot 6: KOOS Estimated Score
p6 <- ggplot(all_scen, aes(time, KOOS_est, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Estimated KOOS Score (0–100, higher = better)",
       x = "Time (days)", y = "KOOS (estimated)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_oa
print(p6)

# Plot 7: Sensitivity analysis — JSW at 2 years
sens_end <- sens_df %>%
  filter(time == max(time)) %>%
  group_by(SensParam, SensMult) %>%
  summarise(JSW_end = mean(JSW), Pain_end = mean(VASPain), .groups = "drop")

p7 <- ggplot(sens_end, aes(SensParam, JSW_end, fill = SensMult)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("0.7" = "#4CAF50", "1.0" = "#2196F3", "1.3" = "#F44336"),
                    labels = c("−30%", "Baseline", "+30%")) +
  labs(title = "Sensitivity Analysis: JSW at 2 Years",
       x = "Parameter", y = "JSW (mm)", fill = "Parameter Change") +
  coord_flip() +
  theme_oa
print(p7)

# Summary table at 1 year (day 365)
summary_1yr <- all_scen %>%
  filter(abs(time - 365) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    VAS_Pain     = round(mean(VASPain), 1),
    JSW_mm       = round(mean(JSW), 2),
    KOOS         = round(mean(KOOS_est), 1),
    ColII_pct    = round(mean(ColII), 1),
    Aggrecan_pct = round(mean(Aggrecan), 1),
    MMP13_ngmL   = round(mean(MMP13), 2),
    uCTXII       = round(mean(uCTXII), 3),
    COMP_ug_mL   = round(mean(COMP_s), 2),
    .groups = "drop"
  )

cat("\n═══════════ OA QSP: 1-Year Outcome Summary ═══════════\n")
print(summary_1yr, n = 10)

cat("\n═══════════ Parameter Reference ═══════════\n")
cat("KL Grade at baseline:", param(mod)$KL_grade, "\n")
cat("Patient Age:", param(mod)$Age_yr, "years\n")
cat("JSW_BASE (normal reference):", param(mod)$JSW_BASE, "mm\n")
cat("Simulation period: 730 days (2 years)\n")
