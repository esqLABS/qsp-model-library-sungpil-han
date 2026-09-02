$PROB
Endometriosis QSP ODE Model — REFACTORED (pluggable PK)
24-compartment PK/PD model
E2 dual-feedback HPO + local E2/aromatase/PGE2 loop + pain/BMD
Same math as endo_mrgsolve_model.R, reorganized to the fork-wide PK
naming convention (GUT_/DEPOT_/CENT_/CL_/V1_/C_/EFFECT_<STEM>).

$PARAM @annotated
// ---- Leuprolide (LEUP), GnRH agonist depot — bespoke dual-absorption
//      archetype: a rapid GUT path plus a slow-release SC DEPOT path,
//      both feeding one CENT compartment, no peripheral compartment.
//      Renamed from Gut_leup/Depot_leup/Central_leup/ka_leup/ka_dep/
//      Vd_leup/CL_leup/F_leup/dose_dep. ----
KA_LEUP       : 0.8    : Rapid-phase gut absorption rate constant (/h)
KDEP_LEUP     : 0.005  : Slow-release SC depot release rate constant (/h)
V1_LEUP       : 20.0   : Central volume of distribution (L)
CL_LEUP       : 8.0    : Clearance (L/h)
F_LEUP        : 0.95   : Bioavailability (depot SC)
DOSE_DEP_LEUP : 3750   : Leuprolide depot dose reference (mcg / 28d) (UNUSED in model body — same as the original, actual dose amount is set directly by ev() amt= at the R level)

// ---- Elagolix (ELA) PK — archetype 3 minus peripheral (gut + central).
//      Renamed from Gut_ela/Central_ela/ka_ela/Vd_ela/CL_ela/F_ela. ----
KA_ELA : 1.2   : Absorption rate constant (/h)
V1_ELA : 90.0  : Volume of distribution (L)
CL_ELA : 12.0  : Clearance (L/h)
F_ELA  : 0.57  : Oral bioavailability

// ---- Dienogest (DIE) PK — archetype 3 minus peripheral (gut + central).
//      Renamed from Gut_die/Central_die/ka_die/Vd_die/CL_die/F_die. ----
KA_DIE : 1.5   : Absorption rate constant (/h)
V1_DIE : 47.0  : Volume of distribution (L)
CL_DIE : 6.4   : Clearance (L/h)
F_DIE  : 0.91  : Oral bioavailability

// ---- Letrozole (LET) PK — archetype 3 minus peripheral (gut + central).
//      Renamed from Gut_let/Central_let/ka_let/Vd_let/CL_let/F_let. ----
KA_LET : 1.8   : Absorption rate constant (/h)
V1_LET : 103.0 : Vd = 1.87 L/kg x 55 kg (L)
CL_LET : 2.1   : Clearance (L/h)
F_LET  : 0.99  : Oral bioavailability

// ---- Norethindrone acetate (NETA, add-back) PK — archetype 3 minus
//      peripheral (gut + central), simplified.
//      Renamed from Gut_neta/Central_neta/ka_neta/Vd_neta/CL_neta/F_neta. ----
KA_NETA : 1.3  : Absorption rate (/h)
V1_NETA : 60.0 : Vd (L)
CL_NETA : 7.5  : Clearance (L/h)
F_NETA  : 0.65 : Bioavailability

// ---- HPO axis parameters (unchanged — physiological feedback constants,
//      not part of any one compound stem) ----
GnRH_base    : 100.0  : Baseline GnRH signal (arbitrary units)
k_FSH_prod   : 0.5    : FSH production rate (IU/L/h)
k_FSH_deg    : 0.08   : FSH degradation rate (/h)
FSH_base     : 5.0    : Baseline FSH (IU/L)
k_LH_prod    : 0.4    : LH production rate (IU/L/h)
k_LH_deg     : 0.15   : LH degradation rate (/h)
LH_base      : 5.0    : Baseline LH (IU/L)
k_E2_prod    : 0.6    : E2 production rate (pg/mL/h)
k_E2_deg     : 0.3    : E2 degradation rate (/h)
E2_base      : 100.0  : Baseline E2 (pg/mL, mid-follicular)
k_P4_prod    : 0.4    : P4 production rate (ng/mL/h)
k_P4_deg     : 0.2    : P4 degradation rate (/h)
P4_base      : 2.0    : Baseline P4 (ng/mL, follicular phase)
IC50_E2_FSH  : 150.0  : IC50 of E2 on FSH secretion (pg/mL)
IC50_E2_LH   : 200.0  : IC50 of E2 on LH secretion (pg/mL)
hill_E2      : 2.0    : Hill coefficient for E2 feedback

// ---- Desensitization (GnRH agonist) — unchanged; consumed by EFFECT_LEUP
//      and EFFECT_ELA below inside the shared GnRHR_occ ODE state ----
k_desens_on  : 0.001  : Rate of GnRH-R down-regulation (/h)
k_desens_off : 0.0001 : Rate of GnRH-R recovery (/h)
GnRHR_base   : 1.0    : Baseline GnRH-R activity (normalized)

// ---- Lesion dynamics (unchanged) ----
k_lesion_grow : 0.0015 : Lesion growth rate (/h)
k_lesion_die  : 0.0008 : Lesion regression/apoptosis rate (/h)
Lesion_base   : 2.0    : Baseline lesion volume (cm3)
E2_EC50_lesion: 50.0   : E2_local EC50 for lesion proliferation (pg/mL equiv)
E2_hill_lesion: 1.5    : Hill for E2 on lesion

// ---- Inflammation: IL-6 (unchanged) ----
k_IL6_prod   : 0.05   : IL-6 production rate (pg/mL/h)
k_IL6_deg    : 0.5    : IL-6 degradation rate (/h)
IL6_base     : 50.0   : Baseline IL-6 (pg/mL)
IL6_lesion_k : 0.02   : IL-6 induction by lesion (per cm3/h)

// ---- Inflammation: PGE2 (unchanged) ----
k_PGE2_prod  : 0.1    : PGE2 production rate (pg/mL/h)
k_PGE2_deg   : 0.8    : PGE2 degradation rate (/h)
PGE2_base    : 200.0  : Baseline PGE2 (pg/mL)
PGE2_E2_k    : 0.05   : PGE2 induction by E2_local (per pg/mL/h)

// ---- NGF & pain (unchanged) ----
k_NGF_prod   : 0.02   : NGF production rate (pg/mL/h)
k_NGF_deg    : 0.15   : NGF degradation rate (/h)
NGF_base     : 100.0  : Baseline NGF (pg/mL)
NGF_IL6_k    : 0.001  : NGF induction by IL-6 (per pg/mL/h)
k_pain_on    : 0.1    : Pain activation rate (NRS/h per NGF above baseline)
k_pain_off   : 0.05   : Pain resolution rate (/h)
Pain_base    : 4.0    : Baseline pain NRS (0-10)

// ---- BMD (unchanged) ----
k_BMD_loss    : 0.0000015 : Rate of BMD loss per unit E2-deficiency (/h)
k_BMD_restore : 0.00001   : Rate of BMD restoration (/h)
BMD_base      : 1.0        : Baseline BMD lumbar (g/cm2)
E2_ref_BMD    : 60.0       : Reference E2 for zero net BMD change (pg/mL)

// ---- Aromatase at lesion (unchanged) ----
k_arom_prod   : 0.3    : Aromatase induction rate
k_arom_deg    : 0.2    : Aromatase degradation rate (/h)
Arom_base     : 1.0    : Baseline aromatase activity (normalized)
E2_local_base : 300.0  : Baseline E2 at lesion (pmol/L)
k_E2local_prod: 0.25   : E2 local production rate
k_E2local_deg : 0.18   : E2 local degradation rate (/h)
PGE2_arom_k   : 0.0005 : PGE2 -> aromatase induction
E2_arom_k     : 0.0003 : E2_local -> aromatase positive feedback

// ---- Drug Hill/effect parameters, renamed to the fork-wide convention.
//      Values copied verbatim from the original (rename, not a refit).
//      Elagolix used pow(C_ela, 1.0) explicitly -> GAMMA_ELA = 1.
//      None of Letrozole/Dienogest had an explicit Hill exponent either
//      -> GAMMA_LET = GAMMA_DIE = 1. ----
EMAX_ELA : 0.95 : Elagolix max GnRH-R blockade (was Emax_GnRHR_ela)
EC50_ELA : 3.0  : Elagolix EC50 for GnRH-R antagonism, ng/mL (was IC50_GnRHR_ela)
GAMMA_ELA: 1.0  : Elagolix Hill coefficient (explicit pow(.,1.0) in original)

EMAX_LET : 0.99 : Letrozole max aromatase inhibition (was Emax_aro_let)
EC50_LET : 1.0  : Letrozole EC50 for aromatase inhibition, ng/mL (was IC50_aro_let)
GAMMA_LET: 1.0  : Letrozole Hill coefficient (none in original; =1)

EMAX_DIE : 0.80 : Dienogest max antiproliferative effect (was Emax_prol_die)
EC50_DIE : 5.0  : Dienogest EC50 for antiproliferative effect, ng/mL (was IC50_prol_die)
GAMMA_DIE: 1.0  : Dienogest Hill coefficient (none in original; =1)

// NETA reused dienogest’s own potency constants in the original
// (0.5 * Emax_prol_die / (IC50_prol_die + C_neta)) rather than declaring
// its own. EMAX_NETA = 0.5 x 0.80 = 0.40 folds that 0.5 factor in so
// EFFECT_NETA = EMAX_NETA*C_NETA/(EC50_NETA+C_NETA) reproduces the exact
// same number; EC50_NETA reuses the identical 5.0 value. See notes.
EMAX_NETA : 0.40 : NETA max add-back bone-protective effect (= 0.5 x original Emax_prol_die; promotion of a hardcoded 0.5 factor into a named param, not a refit)
EC50_NETA : 5.0  : NETA EC50, ng/mL (reuses the original’s IC50_prol_die value directly)
GAMMA_NETA: 1.0  : NETA Hill coefficient (none in original; =1)

// ---- NSAID/COX-2 constants — declared in the ORIGINAL file but never
//      referenced by any $ODE/$TABLE expression (comment there says PGE2
//      inhibition by NSAIDs is "handled externally via flag, not modeled
//      here as explicit ODE"). No NSAID compound is in this refactor’s
//      scope. Preserved as-is (UNUSED), same as the original. ----
IC50_pge2_cox : 0.5   : NSAID IC50 for COX-2/PGE2 (ng/mL) (UNUSED in model body)
Emax_pge2_cox : 0.85  : Emax PGE2 suppression (NSAID/celecoxib) (UNUSED in model body)

// ---- Flags (0 = off, 1 = on) — unchanged names; not part of the PK
//      naming convention’s scope (scenario toggles, not compound stems) ----
use_leup     : 0  : Leuprolide depot flag
use_ela_low  : 0  : Elagolix 150 mg/d flag
use_ela_high : 0  : Elagolix 200 mg BID flag
use_die      : 0  : Dienogest flag
use_let      : 0  : Letrozole flag
use_neta     : 0  : Norethindrone acetate add-back flag
use_ocp      : 0  : Combined OCP flag

$CMT @annotated
Gut_leup      : Leuprolide rapid-phase gut compartment (mcg)
Depot_leup    : Leuprolide SC depot, slow release (mcg)
Central_leup  : Leuprolide central plasma (mcg)
Gut_ela       : Elagolix gut (mg)
Central_ela   : Elagolix central plasma (mg)
Gut_die       : Dienogest gut (mg)
Central_die   : Dienogest central plasma (mg)
Gut_let       : Letrozole gut (mg)
Central_let   : Letrozole central plasma (mg)
Gut_neta      : NETA gut (mg)
Central_neta  : NETA central plasma (mg)
GnRHR_occ     : GnRH receptor activity (normalized, 0-1 active)
FSH_plasma    : FSH in plasma (IU/L)
LH_plasma     : LH in plasma (IU/L)
E2_plasma     : Estradiol plasma (pg/mL)
P4_plasma     : Progesterone plasma (ng/mL)
Lesion        : Ectopic lesion volume (cm3)
IL6_peritoneal: Peritoneal IL-6 (pg/mL)
PGE2_local    : Local PGE2 at lesion (pg/mL)
NGF_lesion    : NGF at lesion (pg/mL)
Pain_score    : Pain NRS (0-10)
BMD_lumbar    : Lumbar spine BMD (g/cm2)
Aromatase_act : Aromatase activity at lesion (normalized)
E2_local      : E2 at ectopic lesion (pmol/L)

$MAIN
// --- Canonical, single-formula exposed concentrations. This is the
//     normalized redirect point: the original computed this same formula
//     a second time under a different name (C_leup_out etc.) purely for
//     $TABLE capture — here both sites share the identical formula and
//     canonical C_<STEM> name. mrgsolve hoists this "double C_<STEM> = ..."
//     declaration into storage shared with $ODE/$TABLE, so $TABLE below
//     does not redeclare it — it does a bare reassignment (no "double") to
//     refresh it from the live, post-solve state at each reporting time,
//     matching what the original’s separate $TABLE computation produced. ---
double C_LEUP = Central_leup / V1_LEUP;        // mcg/L = ng/mL
double C_ELA  = Central_ela  / V1_ELA * 1000;  // mg/L -> ng/mL (x1000)
double C_DIE  = Central_die  / V1_DIE * 1000;  // ng/mL
double C_LET  = Central_let  / V1_LET * 1000;  // ng/mL
double C_NETA = Central_neta / V1_NETA * 1000; // ng/mL

// --- Elagolix (ELA): competitive GnRH-R antagonist. Renamed from
//     ela_inh; a rename of an already Hill-shaped ratio, not a refit. ---
double EFFECT_ELA = (use_ela_low + use_ela_high) > 0 ?
    EMAX_ELA * pow(C_ELA, GAMMA_ELA) / (pow(EC50_ELA, GAMMA_ELA) + pow(C_ELA, GAMMA_ELA)) : 0.0;

// --- Leuprolide (LEUP): GnRH agonist. BESPOKE — not a Hill ratio. The
//     original models the well-known agonist flare/desensitization
//     pharmacology as an explicit ODE state (GnRHR_occ, below in $ODE),
//     not as an algebraic concentration-response curve. EFFECT_LEUP is a
//     gated pass-through of C_LEUP (renamed from leup_active) that drives
//     the desensitization term of that shared ODE; the actual dynamic
//     (flare vs. suppress) lives in the GnRHR_occ ODE itself, unchanged
//     from the original. See endo_refactor_notes.md. ---
double EFFECT_LEUP = (use_leup > 0 && C_LEUP > 0.01) ? C_LEUP : 0.0;

// OCP suppression effect on E2 (not one of this refactor’s 5 compounds —
// a regimen flag, not an independently-dosed drug PK — left unchanged) ---
double ocp_supp = (use_ocp > 0) ? 0.70 : 0.0;

// --- Letrozole (LET): aromatase inhibitor. Renamed from ai_inh. ---
double EFFECT_LET = (use_let > 0) ?
    EMAX_LET * pow(C_LET, GAMMA_LET) / (pow(EC50_LET, GAMMA_LET) + pow(C_LET, GAMMA_LET)) : 0.0;

// --- Dienogest (DIE): antiproliferative progestin. Renamed from
//     die_antiprol. ---
double EFFECT_DIE = (use_die > 0) ?
    EMAX_DIE * pow(C_DIE, GAMMA_DIE) / (pow(EC50_DIE, GAMMA_DIE) + pow(C_DIE, GAMMA_DIE)) : 0.0;

// --- Norethindrone acetate (NETA): add-back bone-protective progestin.
//     Renamed from neta_bmd_prot; EMAX_NETA already folds in the
//     original’s hardcoded 0.5 scale factor (see $PARAM note), so this
//     reproduces the exact same number as the original’s
//     "0.5 * Emax_prol_die * C_neta / (IC50_prol_die + C_neta)". ---
double EFFECT_NETA = (use_neta > 0) ?
    EMAX_NETA * pow(C_NETA, GAMMA_NETA) / (pow(EC50_NETA, GAMMA_NETA) + pow(C_NETA, GAMMA_NETA)) : 0.0;

// OCP dienogest-component antiproliferative (not one of this refactor’s
// 5 compounds — left unchanged) ---
double ocp_antiprol = (use_ocp > 0) ? 0.5 : 0.0;

$ODE

// =============================================================
// PK ODEs
// =============================================================

// --- Leuprolide ---
// Depot SC -> plasma (slow release)
dxdt_Depot_leup = -KDEP_LEUP * Depot_leup;
// Gut (rapid phase)
dxdt_Gut_leup   = -KA_LEUP * Gut_leup;
// Central plasma
dxdt_Central_leup = F_LEUP * KA_LEUP * Gut_leup
                  + F_LEUP * KDEP_LEUP * Depot_leup
                  - (CL_LEUP / V1_LEUP) * Central_leup;

// --- Elagolix ---
dxdt_Gut_ela     = -KA_ELA * Gut_ela;
dxdt_Central_ela = F_ELA * KA_ELA * Gut_ela - (CL_ELA / V1_ELA) * Central_ela;

// --- Dienogest ---
dxdt_Gut_die     = -KA_DIE * Gut_die;
dxdt_Central_die = F_DIE * KA_DIE * Gut_die - (CL_DIE / V1_DIE) * Central_die;

// --- Letrozole ---
dxdt_Gut_let     = -KA_LET * Gut_let;
dxdt_Central_let = F_LET * KA_LET * Gut_let - (CL_LET / V1_LET) * Central_let;

// --- NETA ---
dxdt_Gut_neta     = -KA_NETA * Gut_neta;
dxdt_Central_neta = F_NETA * KA_NETA * Gut_neta - (CL_NETA / V1_NETA) * Central_neta;

// =============================================================
// HPO Axis ODEs
// =============================================================

// --- GnRH-R occupancy (normalized activity 0-1) ---
// Increases due to baseline GnRH drive; decreases due to desensitization
// (leuprolide, via EFFECT_LEUP) or antagonist blockade (elagolix, via
// EFFECT_ELA). In natural state, GnRHR_occ ~ GnRHR_base = 1. This ODE
// structure — the actual flare/suppress pharmacology of a GnRH agonist —
// is kept exactly as the original built it; only the driving terms are
// renamed.
double gnrhr_drive = k_desens_on * EFFECT_LEUP;   // agonist drives desensitization
double gnrhr_recov = k_desens_off * (GnRHR_base - GnRHR_occ); // recovery toward baseline
double gnrhr_antag = k_desens_on * EFFECT_ELA * GnRHR_occ * 10;  // antagonist blocks
double gnrhr_ocp   = k_desens_on * ocp_supp * GnRHR_occ * 2;  // OCP indirect suppression

dxdt_GnRHR_occ = gnrhr_recov - gnrhr_drive - gnrhr_antag - gnrhr_ocp;
// Clamp: ensure GnRHR_occ remains in [0.01, 1.0] via ceiling/floor
// (mrgsolve does not natively clamp; we use max(0.01, ...) in TABLE)

// --- FSH ---
// Produced proportional to GnRH-R activity, inhibited by E2 (neg feedback)
double E2_eff = (E2_plasma > 0) ? E2_plasma : 0.01;
double E2_FSH_inh = pow(E2_eff, hill_E2) /
                    (pow(IC50_E2_FSH, hill_E2) + pow(E2_eff, hill_E2));

dxdt_FSH_plasma = k_FSH_prod * GnRHR_occ * (1.0 - E2_FSH_inh)
                - k_FSH_deg * (FSH_plasma - FSH_base);

// --- LH ---
double E2_LH_inh = pow(E2_eff, hill_E2) /
                   (pow(IC50_E2_LH, hill_E2) + pow(E2_eff, hill_E2));

dxdt_LH_plasma = k_LH_prod * GnRHR_occ * (1.0 - E2_LH_inh)
               - k_LH_deg * (LH_plasma - LH_base);

// --- Estradiol (E2) plasma ---
// Produced by FSH/LH-driven granulosa/theca cells
// OCP and GnRH suppression reduce E2 production
double FSH_stim_E2 = (FSH_plasma > 0) ? FSH_plasma / FSH_base : 1.0;
double LH_stim_E2  = (LH_plasma  > 0) ? LH_plasma  / LH_base  : 1.0;

dxdt_E2_plasma = k_E2_prod * FSH_stim_E2 * LH_stim_E2 * (1.0 - ocp_supp)
               - k_E2_deg * E2_plasma;

// --- Progesterone (P4) plasma ---
// Simplified: P4 driven by LH (luteal phase approximation, time-averaged)
dxdt_P4_plasma = k_P4_prod * LH_stim_E2 * (1.0 - ocp_supp)
               - k_P4_deg * P4_plasma;

// =============================================================
// Disease ODEs
// =============================================================

// --- Aromatase activity at lesion ---
// Induced by: PGE2 (via cAMP/SF-1), E2_local (positive feedback)
// Inhibited by: letrozole (via EFFECT_LET)
double PGE2_stim_arom = PGE2_arom_k * PGE2_local;
double E2loc_stim_arom = E2_arom_k * E2_local;
double arom_AI_inh = EFFECT_LET;  // letrozole

dxdt_Aromatase_act = k_arom_prod * (1.0 + PGE2_stim_arom + E2loc_stim_arom)
                   - k_arom_deg * Aromatase_act
                   - arom_AI_inh * k_arom_deg * Aromatase_act;

// --- E2 local at lesion ---
// Produced by aromatase from circulating androgens
// Augmented by systemic E2
// Inhibited by letrozole (via EFFECT_LET)
double E2_sys_contrib = 0.001 * E2_plasma;  // systemic -> local (pmol/L scale)
dxdt_E2_local = k_E2local_prod * Aromatase_act * (1.0 - EFFECT_LET) + E2_sys_contrib
              - k_E2local_deg * E2_local;

// --- IL-6 peritoneal ---
// Induced by: lesion volume (lesion cells produce IL-6)
// Degraded naturally
double IL6_induction = IL6_lesion_k * Lesion;
double die_IL6_supp = (use_die > 0 || use_ocp > 0) ?
    0.3 * (EFFECT_DIE + ocp_antiprol) : 0.0;

dxdt_IL6_peritoneal = k_IL6_prod + IL6_induction
                    - k_IL6_deg * IL6_peritoneal
                    - die_IL6_supp * IL6_peritoneal;

// --- PGE2 local ---
// Induced by: COX-2 (driven by E2_local, IL-6, NFkB)
// Degraded naturally
// Inhibited by NSAIDs (handled externally via flag, not modeled here as
// explicit ODE — same as the original; no NSAID compound in this refactor)
double E2_loc_pge2 = PGE2_E2_k * E2_local;
double IL6_pge2    = 0.01 * IL6_peritoneal;

dxdt_PGE2_local = k_PGE2_prod * (1.0 + E2_loc_pge2 + IL6_pge2)
                - k_PGE2_deg * PGE2_local;

// --- NGF at lesion ---
// Key pain mediator; induced by IL-6 and local inflammation
double NGF_IL6_stim = NGF_IL6_k * IL6_peritoneal;

dxdt_NGF_lesion = k_NGF_prod * (1.0 + NGF_IL6_stim)
                - k_NGF_deg * NGF_lesion;

// --- Pain score (NRS 0-10) ---
// Driven by NGF and PGE2 (peripheral + central sensitization)
// Modeled as dynamic process with activation and offset
double NGF_above_base = (NGF_lesion > NGF_base) ? (NGF_lesion - NGF_base) / NGF_base : 0.0;
double PGE2_pain_stim = 0.005 * (PGE2_local > PGE2_base ?
    PGE2_local - PGE2_base : 0.0);
double pain_stim = k_pain_on * (NGF_above_base + PGE2_pain_stim);

dxdt_Pain_score = pain_stim * (10.0 - Pain_score)
                - k_pain_off * (Pain_score - Pain_base);

// --- Lesion volume ---
// Growth driven by E2_local (via ER -> proliferation), IL-6, VEGF (implicit)
// Regression driven by apoptosis, immune clearance, progestins
// E2_local in pmol/L -> convert to pg/mL equiv (1 pmol/L ~ 0.272 pg/mL for E2)
double E2_local_pgmL = E2_local * 0.272;
double E2_grow_stim = pow(E2_local_pgmL, E2_hill_lesion) /
    (pow(E2_EC50_lesion, E2_hill_lesion) + pow(E2_local_pgmL, E2_hill_lesion));

// Progestin antiproliferative effect (dienogest + OCP-dienogest-component
// + NETA add-back, at the reduced 0.3 weight the original used)
double total_progestin_eff = EFFECT_DIE + ocp_antiprol
    + ((use_neta > 0) ? 0.3 * EFFECT_NETA : 0.0);
total_progestin_eff = (total_progestin_eff > 0.9) ? 0.9 : total_progestin_eff;

// Net lesion dynamics
dxdt_Lesion = k_lesion_grow * E2_grow_stim * (1.0 - total_progestin_eff)
            + 0.0001 * IL6_peritoneal  // IL-6 promotes lesion growth
            - k_lesion_die * Lesion;

// --- BMD lumbar ---
// Bone loss accelerated when E2 < reference (e.g., post-GnRHa)
// Add-back therapy (NETA, via EFFECT_NETA) partially protects
double E2_deficit = (E2_ref_BMD - E2_plasma > 0) ? E2_ref_BMD - E2_plasma : 0.0;
double bmd_loss_rate = k_BMD_loss * E2_deficit * (1.0 - EFFECT_NETA * 0.7);
double bmd_restore_rate = k_BMD_restore * E2_plasma / E2_ref_BMD;

dxdt_BMD_lumbar = bmd_restore_rate * (BMD_base - BMD_lumbar) - bmd_loss_rate;

$TABLE
// Canonical exposed concentrations. mrgsolve hoists a $MAIN "double X = ..."
// declaration into file-scope storage shared with $ODE/$TABLE, so C_LEUP
// etc. already exist here — redeclaring them with "double" again is a
// duplicate-declaration compile error (confirmed live). The correct,
// established pattern (see driver-patches/HANDOFF.md’s discoverability-
// audit note) is a bare reassignment: this recomputes the identical
// formula from the live, post-solve state at each reporting time — the
// same normalization the $MAIN copy above documents, just restated at the
// point $CAPTURE needs it, matching what the original’s own separate
// $TABLE recomputation (C_leup_out etc.) produced.
C_LEUP = Central_leup / V1_LEUP;
C_ELA  = Central_ela  / V1_ELA * 1000;
C_DIE  = Central_die  / V1_DIE * 1000;
C_LET  = Central_let  / V1_LET * 1000;
C_NETA = Central_neta / V1_NETA * 1000;

// Named Hill/effect terms, recomputed here with the same formula as
// $MAIN, purely for capture/discovery (feeds no ODE from this block).
// Same hoisting rule applies: bare reassignment, no "double".
EFFECT_ELA = (use_ela_low + use_ela_high) > 0 ?
    EMAX_ELA * pow(C_ELA, GAMMA_ELA) / (pow(EC50_ELA, GAMMA_ELA) + pow(C_ELA, GAMMA_ELA)) : 0.0;
EFFECT_LEUP = (use_leup > 0 && C_LEUP > 0.01) ? C_LEUP : 0.0;
EFFECT_LET = (use_let > 0) ?
    EMAX_LET * pow(C_LET, GAMMA_LET) / (pow(EC50_LET, GAMMA_LET) + pow(C_LET, GAMMA_LET)) : 0.0;
EFFECT_DIE = (use_die > 0) ?
    EMAX_DIE * pow(C_DIE, GAMMA_DIE) / (pow(EC50_DIE, GAMMA_DIE) + pow(C_DIE, GAMMA_DIE)) : 0.0;
EFFECT_NETA = (use_neta > 0) ?
    EMAX_NETA * pow(C_NETA, GAMMA_NETA) / (pow(EC50_NETA, GAMMA_NETA) + pow(C_NETA, GAMMA_NETA)) : 0.0;

// Clamp GnRHR_occ to [0.01, 1.0]
double GnRHR = GnRHR_occ < 0.01 ? 0.01 : (GnRHR_occ > 1.0 ? 1.0 : GnRHR_occ);

// E2 suppression from baseline (%)
double E2_pct_suppress = 100.0 * (1.0 - E2_plasma / E2_base);

// Lesion change from baseline (%)
double Lesion_pct_change = 100.0 * (Lesion - Lesion_base) / Lesion_base;

// Pain reduction from baseline
double Pain_delta = Pain_base - Pain_score;

// BMD percent change from baseline
double BMD_pct = 100.0 * (BMD_lumbar - BMD_base) / BMD_base;

$CAPTURE C_LEUP C_ELA C_DIE C_LET C_NETA
         EFFECT_LEUP EFFECT_ELA EFFECT_DIE EFFECT_LET EFFECT_NETA
         GnRHR E2_pct_suppress Lesion_pct_change Pain_delta BMD_pct

$SIGMA
0.04  // proportional residual error (CV = 20%, var = 0.04 on log scale)

$OMEGA @block
0.09 0.04 0.09   // IIV: ka, Vd, CL (approximate, for future pop-PK expansion)

