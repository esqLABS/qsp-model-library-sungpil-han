## =============================================================================
## Anti-NMDA Receptor Encephalitis (AIE) — QSP mrgsolve ODE Model
## PK/PD REFACTOR (this file)
## =============================================================================
## Disease:   Autoimmune Encephalitis (Anti-NMDAR Encephalitis)
## Reference: Dalmau et al. 2008 Ann Neurol; Titulaer et al. 2013 Lancet Neurol
##            Hughes et al. 2010 NEJM; Bost et al. 2021 J Neuroinflammation
##
## This file renames five compounds to the fork's PK/PD refactor convention
## (FORK_WORKFLOW_GUIDE.md, Part 2) so each exposes a single named
## concentration (C_<STEM>) and a single named Hill effect (EFFECT_<STEM>):
##   - IVIG  (intravenous immunoglobulin)      -- Archetype 2 (2-cmt, no depot)
##   - MP    (methylprednisolone)              -- Archetype 2 (2-cmt, no depot)
##   - RTX   (rituximab, anti-CD20)             -- Archetype 2 (2-cmt, no depot)
##   - TCZ   (tocilizumab, anti-IL-6R)          -- Archetype 2 (2-cmt, no depot)
##   - CPX   (4-OH-cyclophosphamide, active metabolite) -- Archetype 1 (1-cmt)
## See aie_refactor_notes.md for the full account, including:
##   (a) a pre-existing, unrelated build defect (`$SET` missing commas,
##       `$CMT`+`$INIT` jointly redeclaring all compartments, `$CAPTURE`
##       directly re-listing compartment names) fixed here with syntax-only
##       changes and logged as UPSTREAM_ISSUES.md #138;
##   (b) a genuine, in-scope PK bug discovered in the original: TCZ's
##       peripheral compartment was never given its own state -- it was
##       written to reuse the CPX_ACT compartment ("reusing CPX_ACT slot for
##       TCZ2 peripheral"), which silently broke BOTH compounds' kinetics
##       (CPX's own declared CL_CPX/Vc_CPX were dead code; TCZ's peripheral
##       elimination used the wrong volume). Fixed here by giving each
##       compound its own dedicated compartment, using each compound's own
##       already-declared parameter values (no new numbers introduced) --
##       exactly the "clean single site" redirect these two compounds were
##       classified for. This changes CPX- and TCZ-driven trajectories from
##       the original in scenarios 5 and 6; see refactor notes for the
##       disclosed, expected verification deviation this causes;
##   (c) IVIG's own duplicate/inconsistent concentration-reporting site
##       (`$MAIN`'s Hill-consumed `Cp_IVIG_mgl` was x1000 of the separately
##       recomputed, non-multiplied `$TABLE` `Cp_IVIG_out`) normalized to the
##       single canonical `C_IVIG` value PD actually reads (the Hill-consumed
##       one); this changes only the *reported* IVIG concentration column,
##       not any disease-side dynamics (the Hill formula itself is untouched).
##
## Compartments: 25 ODEs (16 disease + 9 PK; was 22 disease+PK, +3 net from
##               giving TCZ a real peripheral compartment)
## Scenarios: 6 treatment strategies (unchanged from original)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# MODEL CODE
# ─────────────────────────────────────────────────────────────────────────────
aie_model_code <- '
$PROB
Anti-NMDA Receptor Encephalitis (AIE) QSP Model v1.0 -- PK/PD refactored
- Immune: GCB -> Plasmablast/LLPC -> Anti-NMDAR IgG (serum & CSF)
- CNS: BBB integrity -> Microglia -> NMDAR surface density
- Neurotransmitter: Glutamate/GABA E/I imbalance
- Clinical: mRS-proxy, Seizures, Cognition, Psychiatry
- Drug PK/PD: IVIG (2-CMT), Methylprednisolone (2-CMT),
              Rituximab (2-CMT), Tocilizumab (2-CMT, dedicated peripheral),
              4-OH-CPX (1-CMT, dedicated)

$SET delta=0.25, end=365, start=0

$INIT
GCB      = 100.0  // Germinal Center B cells (arbitrary units, 100=normal)
PB       = 50.0   // Plasmablasts (short-lived)
LLPC     = 200.0  // Long-Lived Plasma Cells (key Ab source)
MB       = 100.0  // Memory B Cells (CD27+)
AB_SERUM = 0.01   // Serum Anti-NMDAR IgG (relative titer, 1=disease peak); near-zero at true disease onset
AB_CSF   = 0.001  // CSF Anti-NMDAR IgG (relative)

BBB      = 1.0    // BBB integrity index (1=intact, 0=fully disrupted)
MG       = 1.0    // Microglia activation (1=resting, up=activated)
NMDAR    = 1.0    // Surface NMDA-R density (1=normal, down=disease)
GLU      = 1.0    // Synaptic/extrasynaptic glutamate (relative, 1=normal)
IL6_CNS  = 1.0    // CNS IL-6 (relative, 1=baseline)
GFAP     = 1.0    // Reactive astrocyte marker (relative)

CRS      = 0.0    // Clinical severity score (0-10, mRS proxy)
SZ       = 0.0    // Seizure frequency (events/week)
COG      = 1.0    // Cognitive index (1=normal -> 0=severely impaired)
PSY      = 0.0    // Psychiatric symptom score (0-10)

CENT_IVIG = 0.0   // IVIG central (mg) [was IVIG1]
PERI_IVIG = 0.0   // IVIG peripheral (mg) [was IVIG2]
CENT_MP   = 0.0   // Methylprednisolone central (mg) [was MP1]
PERI_MP   = 0.0   // Methylprednisolone peripheral (mg) [was MP2]
CENT_RTX  = 0.0   // Rituximab central (mg) [was RTX1]
PERI_RTX  = 0.0   // Rituximab peripheral (mg) [was RTX2]
CENT_TCZ  = 0.0   // Tocilizumab central (mg) [was TCZ1]
PERI_TCZ  = 0.0   // Tocilizumab peripheral (mg) -- NEW dedicated compartment;
                  // the original never had one of its own (see notes)
CENT_CPX  = 0.0   // 4-OH-Cyclophosphamide active, central (mg) [was CPX_ACT,
                  // now a dedicated compartment -- no longer double-booked as
                  // TCZ2 peripheral]

$PARAM
// ── Disease Natural History ─────────────────────────────────────────────────
// (unchanged -- no compound PK/PD role; not in scope for this refactor)
// GCB dynamics (triggered by antigen stimulus at t=0)
k_GCB_stim  = 0.08    // Antigen-driven GCB expansion rate (/d)
k_GCB_death = 0.015   // GCB basal death/exit
k_GCB_MB    = 0.003   // Memory B -> GCB reactivation

// Plasmablast
k_PB_from_GCB = 0.06  // GCB -> PB differentiation
k_PB_death    = 0.12  // PB half-life ~6 d

// Long-Lived Plasma Cell (slow turnover, key sustained Ab source)
k_LLPC_in     = 0.01  // PB -> LLPC seeding rate
k_LLPC_death  = 0.0006 // LLPC t1/2 ~3 years -> delta ~ 0.0006/d

// Memory B Cell
k_MB_from_GCB = 0.03
k_MB_death    = 0.003

// Antibody dynamics
k_Ab_prod      = 0.0003  // Ab production per LLPC+PB unit (relative/d)
k_Ab_clear_ser = 0.030   // IgG serum clearance (t1/2 = 23 d, ln2/23)
k_Ab_transport = 0.08    // Serum->CSF transfer (BBB-dependent, /d)
k_Ab_clear_CSF = 0.12    // CSF Ab clearance (faster, /d)

// BBB dynamics
k_BBB_repair   = 0.04   // BBB self-repair rate (/d)
k_BBB_dmg_MG   = 0.10   // Microglia-mediated BBB damage
k_BBB_dmg_IL6  = 0.04   // IL-6-mediated BBB damage
BBB_min        = 0.15   // Minimum BBB integrity in severe disease

// Neuroinflammation
k_MG_act       = 0.5    // CSF Ab -> microglia activation (per unit Ab/d)
k_MG_res       = 0.08   // Microglia resolution rate
k_IL6_MG       = 0.25   // MG -> IL-6 secretion
k_IL6_clear    = 0.18   // IL-6 clearance
k_GFAP_MG      = 0.12   // MG -> astrocyte reactivity
k_GFAP_res     = 0.09

// NMDAR dynamics
k_NMDAR_base   = 0.020  // Basal NMDAR synthesis (proportional to deficit)
k_NMDAR_intern = 0.30   // Ab-induced NMDAR internalization (/unit Ab/d)
k_NMDAR_recov  = 0.012  // NMDAR recovery once Ab cleared
NMDAR_min      = 0.08   // Minimum NMDAR density (near-complete loss)

// Glutamate E/I imbalance
k_GLU_exc      = 0.25   // GLU excess when NMDAR down (interneuron hypofunction)
k_GLU_clear    = 0.60   // Glutamate uptake/clearance rate
GLU_max        = 4.0    // Maximum relative glutamate level

// Clinical endpoints
k_CRS_NMDAR    = 3.0    // NMDAR loss -> CRS contribution
k_CRS_GLU      = 1.2    // Excess GLU -> CRS
k_CRS_recover  = 0.015  // CRS spontaneous partial recovery
CRS_max        = 10.0

k_SZ_GLU       = 0.8    // Glu excess -> seizure frequency
SZ_thresh      = 1.6    // GLU threshold for clinical seizures
k_SZ_res       = 0.10

k_COG_loss     = 0.25   // NMDAR loss -> cognitive impairment rate
k_COG_recov    = 0.008  // Cognitive recovery (slow, NMDAR-dependent)

k_PSY_DA       = 0.90   // Dopaminergic disinhibition -> psychosis
k_PSY_res      = 0.06   // Psychiatric symptom resolution

// ── Drug PK Parameters (refactored naming; values unchanged from original) ──
// IVIG: 2-compartment, no depot (Rojas et al. 2015; Wang 2019) -- Archetype 2
CL_IVIG  = 0.210  // L/d (~8.75 mL/h), FcRn-dependent
V1_IVIG  = 3.7    // L central volume [was Vc_IVIG]
V2_IVIG  = 25.0   // L peripheral [was Vp_IVIG]
Q_IVIG   = 1.2    // L/d inter-compartment clearance

// Methylprednisolone (MP): 2-compartment, no depot (Mollmann 1992; Derendorf 1993) -- Archetype 2
CL_MP    = 24.0   // L/h x 24 = 576 L/d
V1_MP    = 28.0   // L (0.4 L/kg x 70 kg) [was Vc_MP]
V2_MP    = 56.0   // L (0.8 L/kg x 70 kg) [was Vp_MP]
Q_MP     = 72.0   // L/d

// Rituximab (RTX): 2-compartment, no depot (Berinstein 1998; Maloney 1997) -- Archetype 2
CL_RTX   = 0.336  // L/d (~14 mL/h); t1/2 ~21 d
V1_RTX   = 3.4    // L [was Vc_RTX]
V2_RTX   = 4.4    // L [was Vp_RTX]
Q_RTX    = 0.43   // L/d

// Tocilizumab (TCZ): 2-compartment, no depot (Nishimoto 2009; Gibiansky 2012) -- Archetype 2
// FIXED (in-scope): the original never gave TCZ its own peripheral state --
// it silently reused the CPX_ACT compartment for it (a genuine PK bug, see
// refactor notes). CL_TCZ/V1_TCZ/V2_TCZ/Q_TCZ are the original values,
// unchanged; V2_TCZ (Vp_TCZ) was declared but never actually used correctly
// by the original’s dxdt_TCZ1/dxdt_CPX_ACT pair -- now properly wired.
CL_TCZ   = 0.55   // L/d (linear + nonlinear components simplified)
V1_TCZ   = 3.5    // L [was Vc_TCZ]
V2_TCZ   = 2.9    // L [was Vp_TCZ -- declared in original but not used correctly]
Q_TCZ    = 0.48   // L/d

// 4-OH-Cyclophosphamide active (CPX): simplified 1-compartment -- Archetype 1
// FIXED (in-scope): the original declared CL_CPX/Vc_CPX for this compound but
// its own dxdt_CPX_ACT never actually used them -- CPX_ACT was overwritten by
// TCZ’s peripheral-elimination formula instead (dead parameters, see notes).
// Values unchanged from the original’s own declaration; now actually wired.
CL_CPX   = 96.0   // L/d (4-OH-CPX ~4 h t1/2)
V1_CPX   = 30.0   // L [was Vc_CPX]

// ── Drug PD Parameters (refactored naming; values unchanged from original) ──
// RTX -> B cell depletion (Emax, Hill) -- rename, not a refit
EC50_RTX  = 8.0    // mcg/mL for 50% B cell depletion
EMAX_RTX  = 0.98   // Max fraction depleted [was Emax_RTX]
GAMMA_RTX = 2.0    // Hill coefficient [was gamma_RTX]

// IVIG -> accelerated endogenous IgG catabolism via FcRn saturation -- rename, not a refit
EC50_IVIG  = 12.0   // mg/mL-equivalent IgG serum concentration [was EC50_IVIG_FcRn]
EMAX_IVIG  = 4.5    // Fold-increase in Ab clearance rate (max) [was Emax_IVIG_cat]
GAMMA_IVIG = 1.0    // New -- original had no explicit Hill exponent (plain ratio)

// Methylprednisolone -> anti-inflammatory (IL-6, TNF suppression) + BBB stabilization
// Two distinct pharmacological actions in the original, both read the same
// concentration/EC50 but have their own Emax -- kept separate per the guide
// ("never collapse several drugs into one shared Hill term"); here it is one
// compound with two genuinely different effect endpoints, not two drugs.
EC50_MP      = 0.25   // mcg/mL (in V1)
EMAX_MP_ANTI = 0.80   // Max anti-inflammatory fraction [was Emax_MP_anti]
EMAX_MP_BBB  = 0.55   // Max BBB stabilization fraction [was Emax_MP_BBB]
GAMMA_MP     = 1.0    // New -- original had no explicit Hill exponent (plain ratio, shared by both endpoints)

// Tocilizumab -> IL-6 blockade (sigmoid Emax) -- rename, not a refit
EC50_TCZ  = 2.5    // mcg/mL (in V1)
EMAX_TCZ  = 0.95   // Max IL-6 signaling inhibition [was Emax_TCZ_IL6]
GAMMA_TCZ = 1.0    // New -- original had no explicit Hill exponent (plain ratio)

// 4-OH-CPX -> lymphocyte/plasma cell killing -- rename, not a refit
EC50_CPX  = 500.0  // ng/mL (V1)
EMAX_CPX  = 0.92   // Max fraction killed [was Emax_CPX_kill]
GAMMA_CPX = 1.0    // New -- original had no explicit Hill exponent (plain ratio)

// Treatment flag parameters (0=off, 1=on; set via events) -- unrelated to any
// compound-PK naming role; dead/unused in the original ($MAIN/$ODE/$TABLE
// never reference it), preserved as-is
DOSE_IVIG_GIVEN = 0     // Flag: IVIG course given (for FcRn saturation model)

$MAIN
// ── Exposed concentrations (single canonical site per compound) ────────────
// Bare assignment (no `double`): each of these is also independently
// recomputed and $CAPTURE-d in $TABLE below (see that block for why), which
// auto-declares a model member for the name -- a second `double` declaration
// here would collide with that auto-declared member under mrgsolve 2.0.1
// ("redefinition of ... member"), confirmed empirically the same way this
// fork’s autoimmune-polyendocrinopathy refactor already documented for its
// own C_CSA/C_JAKI. Formulas are an exact carry-over of the original’s own
// Cp_RTX_mcg / Cp_TCZ_mcg / Cp_MP_mcg / Cp_CPX_ng / Cp_IVIG_mgl (same two-step
// division-then-multiply arithmetic, just written on one line and renamed).
C_RTX  = (CENT_RTX / V1_RTX) * 1000;    // mcg/mL
C_TCZ  = (CENT_TCZ / V1_TCZ) * 1000;    // mcg/mL
C_MP   = (CENT_MP  / V1_MP)  * 1000;    // mcg/mL
C_CPX  = (CENT_CPX / V1_CPX) * 1e6;     // ng/mL
C_IVIG = (CENT_IVIG / V1_IVIG) * 1000;  // mg/mL-equivalent (matches the value the
                                         // original’s Hill term actually consumed,
                                         // Cp_IVIG_mgl -- see refactor notes on the
                                         // duplicate/inconsistent reporting site)

// ── Drug Effect Calculations (Emax/Hill) -- same formulas as the original, ──
// ── renamed. These ARE genuine `double` declarations (unlike C_<STEM> above):
// EFFECT_<STEM> is not itself $CAPTURE-d under this exact name anywhere (only
// a separately-named, independently-recomputed EFFECT_<STEM>_out is, in
// $TABLE below), so there is no auto-declared member to collide with here --
// this is the single declaration site for each of these six names, used
// directly by the disease $ODE equations below (`$MAIN` and `$ODE` share one
// compiled scope; only `$TABLE` is a separate scope for this build).────────
double EFFECT_RTX = EMAX_RTX * pow(C_RTX, GAMMA_RTX) /
             (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX));

double EFFECT_IVIG = EMAX_IVIG * pow(C_IVIG, GAMMA_IVIG) /
              (pow(EC50_IVIG, GAMMA_IVIG) + pow(C_IVIG, GAMMA_IVIG));

double EFFECT_MP_ANTI = EMAX_MP_ANTI * pow(C_MP, GAMMA_MP) /
                 (pow(EC50_MP, GAMMA_MP) + pow(C_MP, GAMMA_MP));
double EFFECT_MP_BBB  = EMAX_MP_BBB  * pow(C_MP, GAMMA_MP) /
                 (pow(EC50_MP, GAMMA_MP) + pow(C_MP, GAMMA_MP));

double EFFECT_TCZ = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ) /
             (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));

double EFFECT_CPX = EMAX_CPX * pow(C_CPX, GAMMA_CPX) /
             (pow(EC50_CPX, GAMMA_CPX) + pow(C_CPX, GAMMA_CPX));

// ── Dopaminergic Disinhibition (NMDAR down -> PV interneuron loss -> DA up) ─
double NMDAR_safe = (NMDAR < 0.01) ? 0.01 : NMDAR;
double DA_disinhibition = (1.0 - NMDAR_safe) * 0.8 + (GLU - 1.0) * 0.2;
if(DA_disinhibition < 0) DA_disinhibition = 0;

// ── Effective BBB factor for Ab transport ────────────────────────────────────
double BBB_open = 1.0 - BBB;  // 0=intact, 1=fully disrupted
if(BBB_open < 0) BBB_open = 0;

$ODE
// ═══════════════════════════════════════════════════════════════════════════
// IMMUNE COMPARTMENTS
// ═══════════════════════════════════════════════════════════════════════════

// GCB: antigen-driven expansion (logistic), memory reactivation,
//      depleted by RTX & CPX
dxdt_GCB = k_GCB_stim * GCB * (1.0 - GCB / 500.0)
           + k_GCB_MB  * MB
           - k_GCB_death * GCB
           - EFFECT_RTX * k_GCB_death * GCB * 5.0   // RTX -> ADCC/CDC of CD20+ GCB
           - EFFECT_CPX * GCB * 0.3;                 // CPX kills proliferating GCB

// Plasmablast: from GCB, short-lived, partially CD20- (less RTX sensitive)
dxdt_PB = k_PB_from_GCB * GCB
          - k_PB_death * PB
          - EFFECT_CPX * PB * 0.5;

// LLPC: seeded from PB, very slow turnover (~years)
//       CPX partially effective; bortezomib more effective (modeled via CPX slot)
dxdt_LLPC = k_LLPC_in  * PB
            - k_LLPC_death * LLPC
            - EFFECT_CPX * LLPC * 0.15;     // CPX partial LLPC effect

// Memory B: from GCB, longer-lived, heavily depleted by RTX
dxdt_MB = k_MB_from_GCB * GCB
          - k_MB_death   * MB
          - EFFECT_RTX * MB * 0.9;          // RTX very effective vs memory B

// ═══════════════════════════════════════════════════════════════════════════
// ANTIBODY DYNAMICS
// ═══════════════════════════════════════════════════════════════════════════

// Serum IgG: produced by LLPC + PB
//            cleared at baseline rate + IVIG-accelerated catabolism
double Ab_prod     = k_Ab_prod * (LLPC + 2.0 * PB);
double Ab_clear    = k_Ab_clear_ser * (1.0 + EFFECT_IVIG);  // IVIG speeds up clearance

dxdt_AB_SERUM = Ab_prod - Ab_clear * AB_SERUM;

// CSF IgG: transported from serum through disrupted BBB
//          (intrathecal synthesis also occurs but simplified here)
dxdt_AB_CSF = k_Ab_transport * AB_SERUM * (BBB_open + 0.05)
              - k_Ab_clear_CSF * AB_CSF;

// ═══════════════════════════════════════════════════════════════════════════
// BBB INTEGRITY
// ═══════════════════════════════════════════════════════════════════════════

// BBB repairs towards 1; damaged by microglia activation and IL-6
// MP stabilizes BBB; TCZ reduces IL-6 damage
double IL6_active = IL6_CNS * (1.0 - EFFECT_TCZ);
double MG_active  = MG > 1.0 ? (MG - 1.0) : 0.0;

dxdt_BBB = k_BBB_repair  * (1.0 - BBB) * (1.0 + EFFECT_MP_BBB)
           - k_BBB_dmg_MG  * MG_active * BBB
           - k_BBB_dmg_IL6 * (IL6_active - 1.0) * BBB;
// Floor constraint (handled by clamping in TABLE)

// ═══════════════════════════════════════════════════════════════════════════
// NEUROINFLAMMATION
// ═══════════════════════════════════════════════════════════════════════════

// Microglia: activated by CSF Ab; resolved by MP anti-inflammatory effect
dxdt_MG = k_MG_act * AB_CSF * (MG > 0 ? 1.0 : 0.0)
          - k_MG_res * MG_active * (1.0 + EFFECT_MP_ANTI * 2.0);

// CNS IL-6: from activated microglia; blocked by TCZ
dxdt_IL6_CNS = k_IL6_MG * MG_active
               - k_IL6_clear * (IL6_CNS - 1.0)
               - EFFECT_TCZ * k_IL6_MG * MG_active;

// Astrocyte reactivity (GFAP proxy)
dxdt_GFAP = k_GFAP_MG * MG_active
             - k_GFAP_res * (GFAP - 1.0);

// ═══════════════════════════════════════════════════════════════════════════
// NMDAR SURFACE DENSITY
// ═══════════════════════════════════════════════════════════════════════════

// NMDAR: baseline synthesis toward 1.0
//        internalization proportional to CSF Ab x surface NMDAR
//        recovery when Ab cleared (slow, half-life ~7-14 d without Ab)
double NMDAR_safe2 = NMDAR > NMDAR_min ? NMDAR : NMDAR_min;

dxdt_NMDAR = k_NMDAR_base  * (1.0 - NMDAR_safe2)
             - k_NMDAR_intern * AB_CSF * NMDAR_safe2
             + k_NMDAR_recov  * (1.0 - NMDAR_safe2) * (AB_CSF < 0.05 ? 1.0 : 0.2);

// ═══════════════════════════════════════════════════════════════════════════
// GLUTAMATE / E:I IMBALANCE
// ═══════════════════════════════════════════════════════════════════════════

// GLU excess: NMDAR down -> PV interneuron hypofunction -> disinhibition
// GLU returns toward 1 when NMDAR recovers
double NMDAR_loss = 1.0 - NMDAR;
if(NMDAR_loss < 0) NMDAR_loss = 0;

dxdt_GLU = k_GLU_exc * NMDAR_loss
           - k_GLU_clear * (GLU - 1.0);
// Floor at 1 (no reduction below baseline in this model)

// ═══════════════════════════════════════════════════════════════════════════
// CLINICAL ENDPOINTS
// ═══════════════════════════════════════════════════════════════════════════

// CRS (mRS proxy): driven by NMDAR loss + glutamate excess
double disease_load = k_CRS_NMDAR * (1.0 - NMDAR) + k_CRS_GLU * (GLU - 1.0)
                    + 0.5 * (MG - 1.0) + 0.3 * (AB_CSF / 0.1);
dxdt_CRS = disease_load - k_CRS_recover * CRS;
// Cap handled in TABLE

// Seizure frequency: excess glutamate above seizure threshold
double GLU_over = GLU > SZ_thresh ? (GLU - SZ_thresh) : 0.0;
dxdt_SZ = k_SZ_GLU * GLU_over - k_SZ_res * SZ;

// Cognitive index: progressive loss proportional to NMDAR deficit x time
// Partial recovery when NMDAR > 0.7
double COG_safe = COG > 0.0 ? COG : 0.0;
dxdt_COG = -k_COG_loss * NMDAR_loss * COG_safe
           + k_COG_recov * (1.0 - COG_safe) * (NMDAR > 0.70 ? 1.0 : 0.0);

// Psychiatric symptoms: dopaminergic disinhibition
dxdt_PSY = k_PSY_DA * DA_disinhibition - k_PSY_res * PSY;

// ═══════════════════════════════════════════════════════════════════════════
// DRUG PK — IVIG (Archetype 2: 2-compartment, no depot)
// ═══════════════════════════════════════════════════════════════════════════
dxdt_CENT_IVIG = -(CL_IVIG + Q_IVIG)/V1_IVIG * CENT_IVIG + Q_IVIG/V2_IVIG * PERI_IVIG;
dxdt_PERI_IVIG =  Q_IVIG/V1_IVIG * CENT_IVIG - Q_IVIG/V2_IVIG * PERI_IVIG;

// ─────────────────────────────────────────────────────────────────────────────
// DRUG PK — METHYLPREDNISOLONE (Archetype 2: 2-compartment, no depot)
// ─────────────────────────────────────────────────────────────────────────────
dxdt_CENT_MP = -(CL_MP + Q_MP)/V1_MP * CENT_MP + Q_MP/V2_MP * PERI_MP;
dxdt_PERI_MP =  Q_MP/V1_MP * CENT_MP - Q_MP/V2_MP * PERI_MP;

// ─────────────────────────────────────────────────────────────────────────────
// DRUG PK — RITUXIMAB (Archetype 2: 2-compartment, no depot)
// ─────────────────────────────────────────────────────────────────────────────
dxdt_CENT_RTX = -(CL_RTX + Q_RTX)/V1_RTX * CENT_RTX + Q_RTX/V2_RTX * PERI_RTX;
dxdt_PERI_RTX =  Q_RTX/V1_RTX * CENT_RTX - Q_RTX/V2_RTX * PERI_RTX;

// ─────────────────────────────────────────────────────────────────────────────
// DRUG PK — TOCILIZUMAB (Archetype 2: 2-compartment, no depot)
// FIXED (in-scope, see header/notes): dedicated PERI_TCZ compartment instead
// of the original’s silent reuse of CPX_ACT. Standard CL/Q/V archetype-2
// formula, using the original’s own CL_TCZ/V1_TCZ/V2_TCZ/Q_TCZ values.
// ─────────────────────────────────────────────────────────────────────────────
dxdt_CENT_TCZ = -(CL_TCZ + Q_TCZ)/V1_TCZ * CENT_TCZ + Q_TCZ/V2_TCZ * PERI_TCZ;
dxdt_PERI_TCZ =  Q_TCZ/V1_TCZ * CENT_TCZ - Q_TCZ/V2_TCZ * PERI_TCZ;

// ─────────────────────────────────────────────────────────────────────────────
// DRUG PK — 4-OH-CYCLOPHOSPHAMIDE (Archetype 1: single compartment)
// FIXED (in-scope, see header/notes): dedicated compartment, no longer
// overwritten by TCZ’s peripheral-elimination formula. This restores the
// original’s own commented-out, never-wired intent:
//   "# dxdt_CPX_ACT = -CL_CPX/Vc_CPX * CPX_ACT  [handled via scenario events]"
// using CL_CPX/V1_CPX exactly as the original declared them.
// ─────────────────────────────────────────────────────────────────────────────
dxdt_CENT_CPX = -CL_CPX/V1_CPX * CENT_CPX;

$TABLE
// ── Clamping (unchanged from original) ───────────────────────────────────────
double BBB_c    = BBB    < BBB_min ? BBB_min : BBB;
double NMDAR_c  = NMDAR  < NMDAR_min ? NMDAR_min : NMDAR;
double GLU_c    = GLU    < 1.0 ? 1.0 : (GLU > GLU_max ? GLU_max : GLU);
double COG_c    = COG    < 0.0 ? 0.0 : (COG > 1.0 ? 1.0 : COG);
double PSY_c    = PSY    < 0.0 ? 0.0 : PSY;
double SZ_c     = SZ     < 0.0 ? 0.0 : SZ;
double CRS_c    = CRS    < 0.0 ? 0.0 : (CRS > CRS_max ? CRS_max : CRS);
double MG_c     = MG     < 1.0 ? 1.0 : MG;

// ── Derived clinical outputs (unchanged from original) ───────────────────────
double NMDAR_pct     = NMDAR_c * 100;
double Ab_norm       = AB_SERUM;   // Relative titer (1=peak)
double mRS_est       = CRS_c > 6.0 ? 6.0 : CRS_c;  // Cap at mRS 6
double response_flag = NMDAR_c > 0.7 ? 1.0 : 0.0;  // Functional recovery

// ── Discoverability fix: single, contiguous, literal `double C_<STEM> = ...;`
// per compound, recomputed independently from state (this file’s own
// $MAIN/$ODE locals are not visible in $TABLE -- confirmed empirically for
// this corpus, see e.g. duchenne-muscular-dystrophy’s refactor notes), and
// matching this file’s own original convention of recomputing every $TABLE
// quantity from state rather than reusing a $MAIN local (the original did
// exactly this for Cp_RTX_out/Cp_TCZ_out/Cp_MP_out/Cp_IVIG_out). This is the
// single canonical `double` declaration for each of these names; $MAIN’s
// touches are bare (see $MAIN comment) to avoid colliding with the member
// mrgsolve auto-declares from the $CAPTURE listing below.
double C_RTX = (CENT_RTX / V1_RTX) * 1000;
double C_TCZ = (CENT_TCZ / V1_TCZ) * 1000;
double C_MP = (CENT_MP / V1_MP) * 1000;
double C_CPX = (CENT_CPX / V1_CPX) * 1e6;
double C_IVIG = (CENT_IVIG / V1_IVIG) * 1000;

double EFFECT_RTX_out = EMAX_RTX * pow(C_RTX, GAMMA_RTX) /
                        (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX));
double EFFECT_IVIG_out = EMAX_IVIG * pow(C_IVIG, GAMMA_IVIG) /
                         (pow(EC50_IVIG, GAMMA_IVIG) + pow(C_IVIG, GAMMA_IVIG));
double EFFECT_MP_ANTI_out = EMAX_MP_ANTI * pow(C_MP, GAMMA_MP) /
                            (pow(EC50_MP, GAMMA_MP) + pow(C_MP, GAMMA_MP));
double EFFECT_MP_BBB_out  = EMAX_MP_BBB  * pow(C_MP, GAMMA_MP) /
                            (pow(EC50_MP, GAMMA_MP) + pow(C_MP, GAMMA_MP));
double EFFECT_TCZ_out = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ) /
                        (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));
double EFFECT_CPX_out = EMAX_CPX * pow(C_CPX, GAMMA_CPX) /
                        (pow(EC50_CPX, GAMMA_CPX) + pow(C_CPX, GAMMA_CPX));

$CAPTURE
NMDAR_pct Ab_norm mRS_est CRS_c SZ_c COG_c PSY_c
BBB_c MG_c
C_IVIG C_MP C_RTX C_TCZ C_CPX
EFFECT_RTX_out EFFECT_IVIG_out EFFECT_MP_ANTI_out EFFECT_MP_BBB_out EFFECT_TCZ_out EFFECT_CPX_out
response_flag
'

# ─────────────────────────────────────────────────────────────────────────────
# Compile model
# ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("aie_qsp_refactored", aie_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# EVENT BUILDER — helper to construct dosing events (cmt= targets renamed to
# match the refactored compartments; amounts/timing unchanged from original)
# ─────────────────────────────────────────────────────────────────────────────

# IVIG: 2 g/kg = 140 g total (70 kg patient), administered over 5 days
#       140,000 mg total / 5 doses = 28,000 mg/dose IV
make_IVIG <- function(start_day) {
  ev(amt=28000, cmt="CENT_IVIG", time=start_day, ii=1, addl=4)  # 5 daily doses
}

# Methylprednisolone: 1 g/d IV x 5 d (1000 mg/dose)
make_MP <- function(start_day) {
  ev(amt=1000, cmt="CENT_MP", time=start_day, ii=1, addl=4)
}

# Rituximab: 375 mg/m2 x 4 weekly doses (1.73 m2 BSA = 648 mg ~ 650 mg)
make_RTX_weekly <- function(start_day) {
  ev(amt=650, cmt="CENT_RTX", time=start_day, ii=7, addl=3)
}

# Rituximab: 1000 mg x 2 (2 weeks apart) — alternative regimen
make_RTX_biweekly <- function(start_day) {
  ev(amt=1000, cmt="CENT_RTX", time=start_day, ii=14, addl=1)
}

# Cyclophosphamide: 750 mg/m2 = 1300 mg q28d (monthly) for 6 cycles
# ~25% -> 4-OH-CPX active metabolite -> CENT_CPX dose = 325 mg
make_CPX <- function(start_day) {
  ev(amt=325, cmt="CENT_CPX", time=start_day, ii=28, addl=5)
}

# Tocilizumab: 8 mg/kg q4wk = 560 mg q28d for 6 months
make_TCZ <- function(start_day) {
  ev(amt=560, cmt="CENT_TCZ", time=start_day, ii=28, addl=5)
}

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIOS (unchanged from original — same amounts, same timing)
# ─────────────────────────────────────────────────────────────────────────────
# Initiation at day 14 post-symptom onset (typical diagnosis delay)
tx_start <- 14

scenarios <- list(
  # 1. Natural history — no treatment
  "01_NatHist" = ev(time=0, amt=0, cmt="CENT_IVIG"),

  # 2. First-Line: IVIG + IV Methylprednisolone (standard of care)
  "02_IVIG_MP" = ev(make_IVIG(tx_start), make_MP(tx_start)),

  # 3. First-Line + Plasmapheresis (PE simulated as 80% IgG removal event)
  #    PE modeled as direct reduction via $INIT manipulation at tx_start
  "03_IVIG_MP_PE" = ev(make_IVIG(tx_start), make_MP(tx_start),
                       # PE: 5 exchanges, start day 14, every other day
                       ev(amt=0, cmt="CENT_IVIG", time=tx_start)   # placeholder PE
                       ),

  # 4. Second-Line: Rituximab (weekly x 4) — given at day 30 (after 1L failure)
  "04_Rituximab" = ev(make_IVIG(tx_start), make_MP(tx_start),
                      make_RTX_weekly(30)),

  # 5. Second-Line: Cyclophosphamide (monthly x 6)
  "05_Cyclophosphamide" = ev(make_IVIG(tx_start), make_MP(tx_start),
                             make_CPX(30)),

  # 6. Refractory: Tocilizumab (IL-6 blockade) — day 60 failure of 2L
  "06_Tocilizumab" = ev(make_IVIG(tx_start), make_MP(tx_start),
                        make_RTX_weekly(30),
                        make_TCZ(60))
)

# ─────────────────────────────────────────────────────────────────────────────
# RUN ALL SCENARIOS
# ─────────────────────────────────────────────────────────────────────────────
run_scenario <- function(scen_name, events) {
  out <- mod %>%
    mrgsim(events=events, end=365, delta=0.5) %>%
    as.data.frame() %>%
    mutate(scenario = scen_name)
  return(out)
}

results <- bind_rows(
  run_scenario("1. No Treatment",         scenarios[["01_NatHist"]]),
  run_scenario("2. IVIG + MP",            scenarios[["02_IVIG_MP"]]),
  run_scenario("3. IVIG + MP + PE",       scenarios[["03_IVIG_MP_PE"]]),
  run_scenario("4. + Rituximab",          scenarios[["04_Rituximab"]]),
  run_scenario("5. + Cyclophosphamide",   scenarios[["05_Cyclophosphamide"]]),
  run_scenario("6. + Tocilizumab",        scenarios[["06_Tocilizumab"]])
)

results$scenario <- factor(results$scenario, levels = unique(results$scenario))

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE AT KEY TIMEPOINTS
# ─────────────────────────────────────────────────────────────────────────────
summary_table <- results %>%
  filter(time %in% c(0, 14, 30, 60, 90, 180, 365)) %>%
  select(scenario, time, NMDAR_pct, mRS_est, SZ_c, COG_c, PSY_c, CRS_c) %>%
  mutate(
    NMDAR_pct = round(NMDAR_pct, 1),
    mRS_est   = round(mRS_est, 1),
    SZ_c      = round(SZ_c, 2),
    COG_c     = round(COG_c, 3),
    PSY_c     = round(PSY_c, 2),
    CRS_c     = round(CRS_c, 2)
  ) %>%
  rename(
    "Scenario"         = scenario,
    "Day"              = time,
    "NMDAR (%)"        = NMDAR_pct,
    "mRS (est)"        = mRS_est,
    "Seizures/wk"      = SZ_c,
    "Cognitive Index"  = COG_c,
    "Psych Score"      = PSY_c,
    "CRS"              = CRS_c
  )
print(summary_table)

# ─────────────────────────────────────────────────────────────────────────────
# CLINICAL TRIAL CALIBRATION NOTES
# ─────────────────────────────────────────────────────────────────────────────
# Parameters calibrated against key clinical evidence:
#
# Titulaer et al. 2013 (Lancet Neurol, n=577):
#   - 81% patients improved with first-line therapy within 4 wk
#   - 97% recovered with or without relapses
#   - Model: NMDAR recovery >70% by wk8 under IVIG+MP -> response_flag=1
#
# Dalmau et al. 2008 (Ann Neurol):
#   - CSF NMDAR Ab critical; serum alone insufficient
#   - k_Ab_transport calibrated so CSF Ab peaks 7-14d after serum peak
#
# Nosadini et al. 2015 (J Neurol Neurosurg Psychiatry):
#   - Rituximab (second-line): 79% improvement in refractory cases
#   - RTX GCB depletion >99% within 1 wk; PB partial
#
# Tatencloux et al. 2015 (Eur J Paed Neurol):
#   - Cyclophosphamide effective in refractory pediatric AIE (small series)
#   - k_LLPC_death enhanced by CPX by ~15%
#
# Lee et al. 2016 (J Neuroimmunol):
#   - Tocilizumab (IL-6 inhibitor): anecdotal/case series evidence
#   - TCZ reduces BBB permeability via IL-6 blockade
#   - VEGF/MMP9 pathway downstream
#
# IVIG PK calibration (Rojas et al. 2015, Neurology):
#   - t1/2 IVIG ~21 d, but FcRn saturation reduces endogenous IgG t1/2 to ~5-8d
#   - EMAX_IVIG = 4.5-fold increase in IgG clearance rate at full dose

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────
sc_colors <- c(
  "1. No Treatment"       = "#E53935",
  "2. IVIG + MP"          = "#1E88E5",
  "3. IVIG + MP + PE"     = "#43A047",
  "4. + Rituximab"        = "#FB8C00",
  "5. + Cyclophosphamide" = "#8E24AA",
  "6. + Tocilizumab"      = "#00ACC1"
)

# Plot 1: NMDAR surface density
p1 <- ggplot(results, aes(time, NMDAR_pct, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  geom_hline(yintercept=70, linetype="dashed", color="grey50", linewidth=0.5) +
  annotate("text", x=350, y=72, label="Recovery threshold (70%)", size=2.8, hjust=1) +
  labs(title="Surface NMDA-R Density", x="Day", y="NMDAR Surface (%)", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Plot 2: Clinical severity (mRS proxy)
p2 <- ggplot(results, aes(time, mRS_est, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  scale_y_continuous(limits=c(0,6), breaks=0:6) +
  labs(title="Clinical Severity (mRS estimate)", x="Day", y="mRS (0-6)", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Plot 3: Serum Anti-NMDAR IgG
p3 <- ggplot(results, aes(time, Ab_norm, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  labs(title="Serum Anti-NMDAR IgG (Relative Titer)", x="Day", y="Relative IgG Titer", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Plot 4: Cognitive index
p4 <- ggplot(results, aes(time, COG_c, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  labs(title="Cognitive Index (1=Normal)", x="Day", y="Cognitive Index", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Plot 5: BBB integrity
p5 <- ggplot(results, aes(time, BBB_c, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  labs(title="BBB Integrity (1=Intact)", x="Day", y="BBB Integrity Index", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Plot 6: Seizure frequency
p6 <- ggplot(results, aes(time, SZ_c, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_color_manual(values=sc_colors) +
  labs(title="Seizure Frequency (events/week)", x="Day", y="Seizures/week", color="Scenario") +
  theme_bw(base_size=9) + theme(legend.position="right")

# Combined plot
combined <- (p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title = "Anti-NMDA Receptor Encephalitis: QSP Model — 6 Treatment Scenarios",
    subtitle = "Dalmau & Graus 2018 | Titulaer et al. 2013",
    theme = theme(plot.title = element_text(face="bold", size=11))
  )
print(combined)

# ─────────────────────────────────────────────────────────────────────────────
# DRUG PK PLOT
# ─────────────────────────────────────────────────────────────────────────────
pk_data <- results %>%
  filter(scenario %in% c("2. IVIG + MP", "4. + Rituximab", "6. + Tocilizumab")) %>%
  select(scenario, time, C_IVIG, C_MP, C_RTX, C_TCZ) %>%
  pivot_longer(cols=c(C_IVIG, C_MP, C_RTX, C_TCZ), names_to="Drug", values_to="Concentration")

pk_data$Drug <- factor(pk_data$Drug,
  labels=c("IVIG (mg/mL)", "MP (mcg/mL)", "RTX (mcg/mL)", "TCZ (mcg/mL)"),
  levels=c("C_IVIG","C_MP","C_RTX","C_TCZ"))

p_pk <- ggplot(pk_data, aes(time, Concentration, color=scenario, linetype=Drug)) +
  geom_line(linewidth=0.8) +
  facet_wrap(~Drug, scales="free_y", ncol=2) +
  labs(title="Drug PK Profiles — Selected Scenarios",
       x="Day", y="Concentration", color="Scenario", linetype="Drug") +
  theme_bw(base_size=9)
print(p_pk)

cat("\nAIE QSP model (refactored) complete.\n")
cat("  Compartments: 25 ODEs (6 immune, 6 CNS, 4 clinical, 9 PK)\n")
cat("  Scenarios: 6 treatment strategies (Natural -> Tocilizumab)\n")
cat("  Key biomarkers: NMDAR%, mRS, CSF Ab, BBB, Cognition, Seizures\n")
