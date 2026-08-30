# ============================================================
# Reactive Arthritis (ReA) — mrgsolve QSP Model
# ============================================================
# Disease: Infection-triggered (Chlamydia/GI pathogens), HLA-B27-associated
# Compartments: 26 total (16 disease + 10 drug PK)
# Treatment scenarios: 5 (natural history → TNF-inhibitor)
# Key references:
#   Carter JD 2010 (Ann Intern Med, PMID 20957185) — chlamydia ReA
#   Kvien TK 1994 (Arthritis Rheum) — ReA outcomes
#   Braun J 2014 (Arthritis Res Ther) — anti-TNF in ReA
# ============================================================
#
# ============================================================
# PK/PD REFACTOR (see rea_refactor_notes.md for full detail)
# ============================================================
# Fork-owned sibling of rea_mrgsolve_model.R (never edit the original --
# see FORK_WORKFLOW_GUIDE.md). Scope: the file's five compounds with rows
# in driver-patches/data/compound_perturbation_census.md -- NSAID, SSZ
# (Sulfasalazine), MTX (Methotrexate), TNFi (renamed TNFI, TNF inhibitor),
# IL17i (renamed IL17I, IL-17 inhibitor) -- each renamed to the guide's
# GUT_/CENT_/KA_/F_/EC50_/EMAX_/GAMMA_/C_/EFFECT_<STEM> convention as
# Archetype 3 minus peripheral (depot + central, linear). All disease/immune
# equations, all 5 scenarios, and all plotting/summary code are otherwise
# byte-identical to the original except for the renamed compartment strings
# noted inline below.
#
# Two pre-existing, compound-unrelated build defects were fixed
# syntax-only, disclosed in full in rea_refactor_notes.md and logged as a
# new translations/UPSTREAM_ISSUES.md entry: (1) $CMT+$INIT jointly
# redeclaring every compartment ("Duplicated model names"), fixed by
# deleting $INIT and moving its values into $MAIN's "<CMT>_0" idiom;
# (2) $ODE's dxdt_CHRON_PATH reading a bare `chlamydia` identifier never
# declared anywhere in the original (not even in this R wrapper), fixed by
# adding `chlamydia : 0` to $PARAM -- matching the value the original's own
# R wrapper already implicitly used everywhere. The checked-in
# rea_mrgsolve_model.R still carries both defects and does not build.

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================================
# MODEL CODE
# ============================================================
rea_code <- '
$PROB
Reactive Arthritis (ReA) QSP Model
- 16 disease/immune compartments (pathogen, innate, adaptive, cytokines, joint, endpoints)
- 10 drug PK compartments (NSAID, SSZ, MTX, TNF-inhibitor, IL-17i)
- Calibrated to: Carter 2010 (chlamydia ReA), Kvien 1994 (outcomes),
                 Braun 2014 (etanercept ReA), Hammer 1995 (SSZ ReA)
Time unit: hours; concentrations in pg/mL (cytokines) or µg/mL (drugs)

Refactor note: PK/PD refactor per FORK_WORKFLOW_GUIDE.md Part 2: NSAID, SSZ, MTX,
TNFi (renamed TNFI), IL17i (renamed IL17I) each renamed to the guide’s
GUT_/CENT_/KA_/F_/EC50_/EMAX_/GAMMA_/C_/EFFECT_<STEM> convention. Every one
is Archetype 3 minus peripheral (depot + central, linear, no back-flux),
using a bespoke KE_<STEM> elimination rate constant rather than CL_<STEM>/
V1_<STEM> because none of the five compounds’ central compartments in the
original are divided by a volume before use (NSAID_C/SSZ_C/MTX_C/TNFi_C/
IL17i_C are themselves the concentration state) -- same precedent as
KE_LOP/KE_OND/KE_RIF in bile-acid-diarrhea/bam_mrgsolve_model_refactored.R.
IL17i and TNFi are antibody biologics but this original models neither as
TMDD: both are a plain first-order depot->central PK feeding a plain
Emax*C/(EC50+C) ratio, no free-receptor/complex compartments anywhere in
the file -- see refactor notes.

Build-fix note: pre-existing upstream defects, both independent of every one of
these five compounds’ own PK (full detail + exact error text in
rea_refactor_notes.md and translations/UPSTREAM_ISSUES.md):
(1) the original declares a standalone $CMT block *and* a separate $INIT
    block assigning starting values to the same compartment names --
    mrgsolve 2.0.1 treats $INIT as its own compartment-declaring block, so
    using both together redeclares every compartment twice ("Duplicated
    model names"). Fixed here, syntax-only, by deleting $INIT and moving
    its assignments to the modern "<CMT>_0" idiom in $MAIN (values
    unchanged).
(2) $ODE’s dxdt_CHRON_PATH reads a bare identifier `chlamydia` that is
    never declared anywhere in $PARAM, $CMT, or the surrounding R wrapper
    (’chlamydia’ was not declared in this scope) -- a missing declaration,
    not a name collision. Fixed here, syntax-only, by adding
    `chlamydia : 0` to $PARAM, matching the value the original’s own R
    wrapper already implicitly used everywhere (it never sets `chlamydia`
    via param()/mrgsim(param=...) in any of its 5 scenarios, so 0 is not an
    invented value -- it is the original’s own default-by-omission).
Both fixes are syntax-only, non-numeric, and unrelated to any scoped
compound’s PK/PD. The checked-in original still has both defects exactly
as written and does not build.

$PARAM @annotated
// ----- Pathogen -----
k_path_decay  : 0.008  : Pathogen natural decay rate (hr-1)
k_path_innate : 0.05   : Innate-mediated pathogen killing (hr-1 per unit INNATE)
k_abx_gi      : 0      : GI antibiotic killing switch (0=off, use event to set)
k_abx_chlamyd : 0      : Anti-chlamydial killing (doxycycline/azithromycin)
k_chron_decay : 0.001  : Chronic chlamydial persistence decay rate (hr-1)
chlamydia     : 0      : build-fix -- never declared in original; 0 matches its own implicit default (see $PROB)

// ----- Innate immune -----
k_innate_act  : 0.3    : Innate immune activation rate by pathogen
k_innate_decay: 0.05   : Innate immune resolution rate (hr-1)
k_neut_prod   : 0.2    : Neutrophil production by INNATE & IL17 (hr-1)
k_neut_decay  : 0.03   : Neutrophil clearance (hr-1)
k_macro_act   : 0.15   : Macrophage activation rate (hr-1)
k_macro_res   : 0.03   : Macrophage resolution rate (hr-1)

// ----- T cell differentiation -----
k_Th1_diff    : 0.08   : Th1 differentiation rate (hr-1)
k_Th17_diff   : 0.07   : Th17 differentiation rate (hr-1)
k_Treg_prod   : 0.03   : Regulatory T cell production (hr-1)
k_Th1_decay   : 0.02   : Th1 decay (hr-1)
k_Th17_decay  : 0.025  : Th17 decay (hr-1)
k_Treg_decay  : 0.01   : Treg decay (hr-1)
HLAB27_eff    : 1.6    : HLA-B27 fold increase in Th17 differentiation
USE_HLAB27    : 1      : 1=HLA-B27 positive, 0=negative

// ----- Cytokines -----
k_TNF_prod    : 4.0    : TNF-alpha production (pg/mL/hr per unit MACRO+TH1)
k_TNF_decay   : 0.4    : TNF-alpha clearance (hr-1, t1/2~1.5h)
k_IL17_prod   : 2.5    : IL-17A production by TH17 (pg/mL/hr)
k_IL17_decay  : 0.25   : IL-17A clearance (hr-1)
k_IL6_prod    : 3.5    : IL-6 production (pg/mL/hr per unit MACRO)
k_IL6_decay   : 0.35   : IL-6 clearance (hr-1, t1/2~2h)
k_IL10_prod   : 1.2    : IL-10 production by TREG (pg/mL/hr)
k_IL10_decay  : 0.2    : IL-10 clearance (hr-1)
k_IFNg_prod   : 1.5    : IFN-gamma production by TH1 (pg/mL/hr)
k_IFNg_decay  : 0.3    : IFN-gamma clearance (hr-1)

// ----- Synovial/joint -----
k_synov_act   : 0.12   : Synovitis activation (per unit cytokine drive)
k_synov_res   : 0.015  : Spontaneous synovitis resolution (hr-1)
k_cart_dmg    : 0.00003: Cartilage damage rate (hr-1 per unit SYNOV)
k_pain_act    : 0.25   : Pain activation coefficient
k_pain_res    : 0.08   : Pain resolution rate (hr-1)
k_CRP_prod    : 1.8    : CRP production per IL6 (mg/L per pg/mL per hr)
k_CRP_decay   : 0.036  : CRP clearance (hr-1, t1/2~19h)
k_jcount_rate : 0.004  : Joint count equilibration rate (hr-1)

// ----- Baseline steady-state values -----
TNF_base      : 3.0    : Baseline TNF-alpha (pg/mL)
IL17_base     : 1.5    : Baseline IL-17A (pg/mL)
IL6_base      : 2.0    : Baseline IL-6 (pg/mL)
CRP_base      : 1.0    : Baseline CRP (mg/L)
PAIN_base     : 5.0    : Baseline VAS pain (0-100)

// ----- NSAID PK/PD (naproxen 500 mg BID) -----
// refactor: ka_NSAID->KA_NSAID, ke_NSAID->KE_NSAID (bespoke rate const,
// no CL/V1 split), F_NSAID unchanged, Vd_NSAID->V1_NSAID (still dead --
// never read anywhere in $ODE, exactly as in the original), IC50_NSAID->
// EC50_NSAID, Emax_NSAID->EMAX_NSAID, GAMMA_NSAID=1 added (original had no
// Hill exponent, plain ratio; pow(x,1)=x makes this a rename not a fit).
KA_NSAID      : 1.2    : NSAID absorption rate (hr-1)
KE_NSAID      : 0.33   : NSAID elimination rate (hr-1, t1/2=2.1h naproxen; bespoke rate const, no volume in original)
F_NSAID       : 0.95   : NSAID oral bioavailability
V1_NSAID      : 0.12   : (was Vd_NSAID) NSAID apparent volume (L/kg) -- dead parameter, unused in $ODE, same as original
EC50_NSAID    : 1.8    : NSAID COX-2 EC50 (µg/mL)
EMAX_NSAID    : 0.85   : Max COX-2 inhibition / PGE2 reduction
GAMMA_NSAID   : 1      : Hill exponent (original had none; rename not a fit)

// ----- Sulfasalazine PK/PD (2 g/day → maintenance 3 g/day) -----
KA_SSZ        : 0.25   : SSZ absorption (hr-1, slow due to colonic cleavage)
KE_SSZ        : 0.033  : SSZ elimination (hr-1, t1/2~21h sulfapyridine; bespoke rate const, no volume in original)
F_SSZ         : 0.10   : SSZ bioavailability (parent compound <10%)
EC50_SSZ      : 10.0   : SSZ EC50 for TNF/IL-6 suppression (µg/mL sulfapyridine)
EMAX_SSZ      : 0.55   : SSZ max cytokine suppression
GAMMA_SSZ     : 1      : Hill exponent (original had none; rename not a fit)

// ----- Methotrexate PK/PD (15 mg/wk SC) -----
KA_MTX        : 0.7    : MTX absorption (hr-1, SC faster than oral)
KE_MTX        : 0.065  : MTX elimination (hr-1, t1/2~10h; bespoke rate const, no volume in original)
F_MTX         : 0.78   : MTX SC bioavailability
EC50_MTX      : 0.45   : MTX EC50 for T cell suppression (µmol/L)
EMAX_MTX      : 0.72   : MTX max lymphocyte suppression
GAMMA_MTX     : 1      : Hill exponent (original had none; rename not a fit)

// ----- TNF inhibitor PK/PD (etanercept 50 mg/wk SC) -----
// refactor: stem TNFi -> TNFI. Checked for TMDD: none present -- plain
// first-order PK feeding a plain Emax*C/(EC50+C) ratio, no free-receptor/
// complex compartments anywhere in the file. See refactor notes.
KA_TNFI       : 0.009  : ETN SC absorption (hr-1, Tmax~72h)
KE_TNFI       : 0.0055 : ETN elimination (hr-1, t1/2~125h; bespoke rate const, no volume in original)
F_TNFI        : 0.76   : ETN SC bioavailability
EC50_TNFI     : 1.2    : TNFi EC50 for TNF neutralization (µg/mL)
EMAX_TNFI     : 0.95   : Max TNF neutralization
GAMMA_TNFI    : 1      : Hill exponent (original had none; rename not a fit)

// ----- IL-17 inhibitor PK/PD (secukinumab 300 mg SC monthly) -----
// refactor: stem IL17i -> IL17I. Checked for TMDD: none present -- same
// plain first-order PK / plain Emax*C/(EC50+C) ratio as TNFI above.
KA_IL17I      : 0.006  : SECU SC absorption (hr-1, Tmax~6d)
KE_IL17I      : 0.0025 : SECU elimination (hr-1, t1/2~27d; bespoke rate const, no volume in original)
F_IL17I       : 0.73   : SECU SC bioavailability
EC50_IL17I    : 0.7    : IL17i EC50 for IL-17A neutralization (µg/mL)
EMAX_IL17I    : 0.93   : Max IL-17A neutralization
GAMMA_IL17I   : 1      : Hill exponent (original had none; rename not a fit)

$CMT @annotated
PATH           : Pathogen load (normalized, 0-1)
CHRON_PATH     : Chronic synovial pathogen persistence (Chlamydia)
INNATE         : Innate immune activation score
NEUT           : Neutrophil activation
MACRO          : Macrophage activation
TH1            : Th1 cell count (normalized)
TH17           : Th17 cell count (normalized)
TREG           : Regulatory T cell count (normalized)
TNF_c          : TNF-alpha (pg/mL)
IL17_c         : IL-17A (pg/mL)
IL6_c          : IL-6 (pg/mL)
IL10_c         : IL-10 (pg/mL)
IFNg_c         : IFN-gamma (pg/mL)
SYNOV          : Synovial inflammation score (0-10)
CARTDMG        : Cumulative cartilage damage (0-1)
PAIN           : VAS pain score (0-100)
CRP_c          : CRP (mg/L)
JCOUNT         : Swollen joint count (0-68)
GUT_NSAID      : (was NSAID_DEPOT) NSAID gut depot (mg equivalent)
CENT_NSAID     : (was NSAID_C) NSAID plasma concentration -- the exposed concentration (µg/mL)
GUT_SSZ        : (was SSZ_DEPOT) SSZ gut depot (mg equivalent)
CENT_SSZ       : (was SSZ_C) SSZ plasma -- the exposed concentration (µg/mL)
GUT_MTX        : (was MTX_DEPOT) MTX SC/oral depot (mg equivalent)
CENT_MTX       : (was MTX_C) MTX plasma -- the exposed concentration (µmol/L)
GUT_TNFI       : (was TNFi_DEPOT) TNFi SC depot (mg)
CENT_TNFI      : (was TNFi_C) TNFi plasma -- the exposed concentration (µg/mL)
GUT_IL17I      : (was IL17i_DEPOT) IL-17i SC depot (mg)
CENT_IL17I     : (was IL17i_C) IL-17i plasma -- the exposed concentration (µg/mL)

$MAIN
// build-fix: see $PROB note (1): $INIT deleted (redeclared every
// compartment against $CMT), values moved here unchanged using the
// modern "<CMT>_0" idiom.
PATH_0       = 1.0;
CHRON_PATH_0 = 0.0;
INNATE_0     = 0.0;
NEUT_0       = 0.5;
MACRO_0      = 0.5;
TH1_0        = 0.3;
TH17_0       = 0.2;
TREG_0       = 0.2;
TNF_c_0      = 3.0;
IL17_c_0     = 1.5;
IL6_c_0      = 2.0;
IL10_c_0     = 2.0;
IFNg_c_0     = 1.0;
SYNOV_0      = 0.2;
CARTDMG_0    = 0.0;
PAIN_0       = 5.0;
CRP_c_0      = 1.0;
JCOUNT_0     = 0.5;
GUT_NSAID_0  = 0;
CENT_NSAID_0 = 0;
GUT_SSZ_0    = 0;
CENT_SSZ_0   = 0;
GUT_MTX_0    = 0;
CENT_MTX_0   = 0;
GUT_TNFI_0   = 0;
CENT_TNFI_0  = 0;
GUT_IL17I_0  = 0;
CENT_IL17I_0 = 0;

$ODE
// ============================================================
// Drug effect functions (Hill equation, Emax model)
// refactor: each compound’s exposed concentration (C_<STEM>) is a direct
// alias of its central compartment -- no volume divides it in the
// original (NSAID_C/SSZ_C/MTX_C/TNFi_C/IL17i_C were themselves the
// concentration state), same as C_LOP/C_OND/C_RIF in
// bile-acid-diarrhea/bam_mrgsolve_model_refactored.R. EFFECT_<STEM> is the
// named Hill interface; GAMMA_<STEM>=1 makes pow(x,1)=x, so this is a
// rename of the original’s E_NSAID/E_SSZ/E_MTX/E_TNFi/E_IL17i ratios, not
// a refit -- same "+1e-12" guard the original used is kept verbatim.
double C_NSAID = CENT_NSAID;
double EFFECT_NSAID = EMAX_NSAID * pow(C_NSAID, GAMMA_NSAID)
                      / (pow(EC50_NSAID, GAMMA_NSAID) + pow(C_NSAID, GAMMA_NSAID) + 1e-12);
double C_SSZ = CENT_SSZ;
double EFFECT_SSZ = EMAX_SSZ * pow(C_SSZ, GAMMA_SSZ)
                    / (pow(EC50_SSZ, GAMMA_SSZ) + pow(C_SSZ, GAMMA_SSZ) + 1e-12);
double C_MTX = CENT_MTX;
double EFFECT_MTX = EMAX_MTX * pow(C_MTX, GAMMA_MTX)
                    / (pow(EC50_MTX, GAMMA_MTX) + pow(C_MTX, GAMMA_MTX) + 1e-12);
double C_TNFI = CENT_TNFI;
double EFFECT_TNFI = EMAX_TNFI * pow(C_TNFI, GAMMA_TNFI)
                     / (pow(EC50_TNFI, GAMMA_TNFI) + pow(C_TNFI, GAMMA_TNFI) + 1e-12);
double C_IL17I = CENT_IL17I;
double EFFECT_IL17I = EMAX_IL17I * pow(C_IL17I, GAMMA_IL17I)
                      / (pow(EC50_IL17I, GAMMA_IL17I) + pow(C_IL17I, GAMMA_IL17I) + 1e-12);

// HLA-B27 modulator
double HB27 = 1.0 + (HLAB27_eff - 1.0) * USE_HLAB27;

// IL-10 suppression term (regulatory brake)
double IL10_supp = 1.0 / (1.0 + IL10_c / 6.0);

// ============================================================
// 1. Pathogen dynamics
// ============================================================
dxdt_PATH = -(k_path_decay + k_path_innate * INNATE + k_abx_gi) * PATH;
// Chlamydia: small fraction seeds persistent synovial compartment
dxdt_CHRON_PATH = 0.002 * chlamydia - k_chron_decay * CHRON_PATH - k_abx_chlamyd * CHRON_PATH;

// ============================================================
// 2. Innate immune
// ============================================================
dxdt_INNATE = k_innate_act * (PATH + CHRON_PATH * 0.5) - k_innate_decay * INNATE;
dxdt_NEUT   = k_neut_prod * (INNATE + IL17_c / 20.0) - k_neut_decay * NEUT;
dxdt_MACRO  = k_macro_act * (INNATE + TH1 * 0.4 + IFNg_c / 20.0)
              - k_macro_res * MACRO;

// ============================================================
// 3. T cell differentiation
// ============================================================
double Th1_in  = k_Th1_diff  * MACRO * IL10_supp * (1.0 - EFFECT_MTX);
double Th17_in = k_Th17_diff * MACRO * HB27 * IL10_supp * (1.0 - EFFECT_MTX * 0.8);
double Treg_in = k_Treg_prod * (1.0 + TH1 * 0.2);

dxdt_TH1  = Th1_in  - k_Th1_decay  * TH1;
dxdt_TH17 = Th17_in - k_Th17_decay * TH17;
dxdt_TREG = Treg_in - k_Treg_decay * TREG;

// ============================================================
// 4. Cytokine dynamics
// ============================================================
double TNF_prod  = k_TNF_prod  * (MACRO + TH1 * 0.7)
                   * (1.0 - EFFECT_TNFI) * (1.0 - EFFECT_SSZ * 0.5) * (1.0 - EFFECT_MTX * 0.6);
double IL17_prod = k_IL17_prod * TH17 * (1.0 - EFFECT_IL17I);
double IL6_prod  = k_IL6_prod  * MACRO * (1.0 - EFFECT_SSZ * 0.45) * (1.0 - EFFECT_MTX * 0.4);
double IL10_prod = k_IL10_prod * TREG;
double IFNg_prod = k_IFNg_prod * TH1;

dxdt_TNF_c  = TNF_prod  - k_TNF_decay  * TNF_c;
dxdt_IL17_c = IL17_prod - k_IL17_decay * IL17_c;
dxdt_IL6_c  = IL6_prod  - k_IL6_decay  * IL6_c;
dxdt_IL10_c = IL10_prod - k_IL10_decay * IL10_c;
dxdt_IFNg_c = IFNg_prod - k_IFNg_decay * IFNg_c;

// ============================================================
// 5. Synovial inflammation & joint damage
// ============================================================
double synov_drive = k_synov_act * (TNF_c / 5.0 + IL17_c / 3.0 + IL6_c / 5.0
                                    + CHRON_PATH * 2.0);
dxdt_SYNOV = synov_drive - k_synov_res * SYNOV;

// Irreversible cartilage damage
dxdt_CARTDMG = k_cart_dmg * SYNOV * (1.0 - CARTDMG);

// Pain: driven by PGE2-proxy (SYNOV, TNF, IL6), reduced by NSAIDs
double pain_drive = k_pain_act * (SYNOV * 8.0 + TNF_c / 4.0 + IL6_c / 6.0);
dxdt_PAIN = pain_drive - k_pain_res * PAIN - EFFECT_NSAID * k_pain_res * PAIN;

// CRP: IL-6 → hepatocyte → CRP (production, t1/2~19h)
dxdt_CRP_c = k_CRP_prod * IL6_c - k_CRP_decay * CRP_c;

// Joint count: slow equilibration to synovitis-based target
double jcount_target = 12.0 * SYNOV / (SYNOV + 1.5);
dxdt_JCOUNT = k_jcount_rate * (jcount_target - JCOUNT);

// ============================================================
// 6. Drug PK (first-order absorption + elimination)
// refactor: pure identifier rename per compound (Archetype 3 minus
// peripheral); no structural change from the original.
// ============================================================
dxdt_GUT_NSAID  = -KA_NSAID * GUT_NSAID;
dxdt_CENT_NSAID =  KA_NSAID * F_NSAID * GUT_NSAID - KE_NSAID * CENT_NSAID;

dxdt_GUT_SSZ    = -KA_SSZ * GUT_SSZ;
dxdt_CENT_SSZ   =  KA_SSZ * F_SSZ * GUT_SSZ - KE_SSZ * CENT_SSZ;

dxdt_GUT_MTX    = -KA_MTX * GUT_MTX;
dxdt_CENT_MTX   =  KA_MTX * F_MTX * GUT_MTX - KE_MTX * CENT_MTX;

dxdt_GUT_TNFI   = -KA_TNFI * GUT_TNFI;
dxdt_CENT_TNFI  =  KA_TNFI * F_TNFI * GUT_TNFI - KE_TNFI * CENT_TNFI;

dxdt_GUT_IL17I  = -KA_IL17I * GUT_IL17I;
dxdt_CENT_IL17I =  KA_IL17I * F_IL17I * GUT_IL17I - KE_IL17I * CENT_IL17I;

$TABLE
// Capture key outputs for plotting
// refactor: identifier rename only (renamed compartments referenced below);
// same original output names preserved for backward comparability.
capture VAS_pain    = PAIN;
capture swollen_jt  = JCOUNT;
capture CRP_obs     = CRP_c;
capture TNF_obs     = TNF_c;
capture IL17_obs    = IL17_c;
capture IL6_obs     = IL6_c;
capture IL10_obs    = IL10_c;
capture IFNg_obs    = IFNg_c;
capture SYNOV_obs   = SYNOV;
capture CARTDMG_pct = CARTDMG * 100.0;
capture PATH_obs    = PATH;
capture CHRON_obs   = CHRON_PATH;
capture TH1_obs     = TH1;
capture TH17_obs    = TH17;
capture TREG_obs    = TREG;
capture Th17_Treg   = TH17 / (TREG + 0.01);
capture NSAID_obs   = CENT_NSAID;
capture SSZ_obs     = CENT_SSZ;
capture MTX_obs     = CENT_MTX;
capture TNFi_obs    = CENT_TNFI;
capture IL17i_obs   = CENT_IL17I;

$CAPTURE @annotated
C_NSAID      : NSAID exposed concentration (µg/mL) -- direct alias of CENT_NSAID, no volume in original
EFFECT_NSAID : NSAID fractional COX-2/PGE2 inhibition (Hill interface, GAMMA=1 rename of original ratio)
C_SSZ        : Sulfasalazine exposed concentration (µg/mL) -- direct alias of CENT_SSZ, no volume in original
EFFECT_SSZ   : Sulfasalazine fractional TNF/IL-6 suppression (Hill interface, GAMMA=1 rename of original ratio)
C_MTX        : Methotrexate exposed concentration (µmol/L) -- direct alias of CENT_MTX, no volume in original
EFFECT_MTX   : Methotrexate fractional T-cell/cytokine suppression (Hill interface, GAMMA=1 rename of original ratio)
C_TNFI       : TNF inhibitor exposed concentration (µg/mL) -- direct alias of CENT_TNFI, no volume in original; not TMDD
EFFECT_TNFI  : TNF inhibitor fractional TNF neutralization (Hill interface, GAMMA=1 rename of original ratio)
C_IL17I      : IL-17 inhibitor exposed concentration (µg/mL) -- direct alias of CENT_IL17I, no volume in original; not TMDD
EFFECT_IL17I : IL-17 inhibitor fractional IL-17A neutralization (Hill interface, GAMMA=1 rename of original ratio)
'

mod <- mcode("ReactiveArthritis_QSP_refactored", rea_code)

# ============================================================
# SCENARIO DEFINITIONS
# ============================================================
# All scenarios: ReA triggered by infection at t=0 (event: PATH=1.0 already at INIT)
# Simulation duration: 365 days = 8760 hours

sim_time <- seq(0, 8760, by = 6)  # 6-hr steps

# ---- Helper: build NSAID dosing events ----
# Naproxen 500 mg BID (every 12h) — dose in mg converted to µg/mL equivalent via depot
# In model, GUT_NSAID gets dose_mg; PK handles absorption
# refactor: cmt string updated for the renamed depot (was "NSAID_DEPOT")
make_nsaid_events <- function(start_hr, dur_hr, dose_mg = 500, interval_hr = 12) {
  times <- seq(start_hr, start_hr + dur_hr - interval_hr, by = interval_hr)
  ev(cmt = "GUT_NSAID", time = times, amt = dose_mg, addl = 0)
}

# refactor: cmt string updated for the renamed depot (was "SSZ_DEPOT")
make_ssz_events <- function(start_hr, dur_hr, dose_mg = 1000, interval_hr = 12) {
  times <- seq(start_hr, start_hr + dur_hr - interval_hr, by = interval_hr)
  ev(cmt = "GUT_SSZ", time = times, amt = dose_mg)
}

# refactor: cmt string updated for the renamed depot (was "MTX_DEPOT")
make_mtx_events <- function(start_hr, start_dose = 15, dur_wk = 48) {
  times <- seq(start_hr, start_hr + (dur_wk - 1) * 168, by = 168)
  ev(cmt = "GUT_MTX", time = times, amt = start_dose)
}

# refactor: cmt string updated for the renamed depot (was "TNFi_DEPOT")
make_tnfi_events <- function(start_hr, dose_mg = 50, interval_hr = 168, n_doses = 52) {
  times <- seq(start_hr, start_hr + (n_doses - 1) * interval_hr, by = interval_hr)
  ev(cmt = "GUT_TNFI", time = times, amt = dose_mg)
}

# refactor: cmt string updated for the renamed depot (was "IL17i_DEPOT")
make_il17i_events <- function(start_hr, dose_mg = 300, interval_hr = 720, n_doses = 12) {
  # Loading: 300 mg at wk 0,1,2,3,4, then monthly
  load_times <- seq(start_hr, start_hr + 4 * 168, by = 168)
  maint_times <- seq(start_hr + 4 * 168 + interval_hr, start_hr + n_doses * interval_hr,
                     by = interval_hr)
  times <- c(load_times, maint_times)
  ev(cmt = "GUT_IL17I", time = times, amt = dose_mg)
}

# ============================================================
# SCENARIO 1: No Treatment (Natural History)
# ============================================================
# HLA-B27 positive (USE_HLAB27=1), Chlamydia GU trigger
s1_params <- list(USE_HLAB27 = 1, k_abx_gi = 0, k_abx_chlamyd = 0)
s1_ev <- ev(time = 0, cmt = "PATH", amt = 0, rate = 0)  # no treatment events

out_s1 <- mod %>%
  param(s1_params) %>%
  mrgsim(ev = s1_ev, end = 8760, delta = 6, carry_out = "time") %>%
  as.data.frame() %>%
  mutate(scenario = "S1: No Treatment", day = time / 24)

cat("Scenario 1 (No Treatment) complete. Final CRP:", round(tail(out_s1$CRP_obs, 1), 1),
    "mg/L; VAS pain:", round(tail(out_s1$VAS_pain, 1), 1), "\n")

# ============================================================
# SCENARIO 2: NSAIDs Only (Naproxen 500 mg BID for 12 weeks)
# ============================================================
s2_params <- list(USE_HLAB27 = 1, k_abx_gi = 0, k_abx_chlamyd = 0)
# NSAID for 12 weeks = 2016 hours
s2_ev <- make_nsaid_events(start_hr = 0, dur_hr = 2016)

out_s2 <- mod %>%
  param(s2_params) %>%
  mrgsim(ev = s2_ev, end = 8760, delta = 6, carry_out = "time") %>%
  as.data.frame() %>%
  mutate(scenario = "S2: NSAIDs (Naproxen 500mg BID)", day = time / 24)

cat("Scenario 2 (NSAID only) complete. 12-wk CRP:",
    round(out_s2$CRP_obs[out_s2$day == 84][1], 1), "mg/L\n")

# ============================================================
# SCENARIO 3: Antibiotics + NSAIDs (Chlamydia-triggered ReA)
# ============================================================
# Doxycycline 100 mg BID × 3 months = k_abx_chlamyd switch; plus NSAIDs
s3_params_list <- lapply(1:nrow(data.frame(time = sim_time)), function(i) {
  data.frame(time = sim_time[i])
})
s3_params <- list(USE_HLAB27 = 1, k_abx_gi = 0.5, k_abx_chlamyd = 0.3)

# Antibiotics for 12 weeks, then stop (parameter switch via tgrid — simplify with constant)
# For simplicity: fixed antibiotic effect during first 2016 hr then 0
# Implement via two-phase simulation
out_s3 <- mod %>%
  param(list(USE_HLAB27 = 1, k_abx_gi = 0.5, k_abx_chlamyd = 0.3)) %>%
  mrgsim(ev    = make_nsaid_events(0, 2016),
         end   = 2016, delta = 6, carry_out = "time") %>%
  as.data.frame()

s3_state <- as.list(tail(out_s3, 1)[, !names(out_s3) %in%
            c("ID","time","day","scenario","VAS_pain","swollen_jt","CRP_obs",
              "TNF_obs","IL17_obs","IL6_obs","IL10_obs","IFNg_obs","SYNOV_obs",
              "CARTDMG_pct","PATH_obs","CHRON_obs","TH1_obs","TH17_obs","TREG_obs",
              "Th17_Treg","NSAID_obs","SSZ_obs","MTX_obs","TNFi_obs","IL17i_obs")])
# Continue without antibiotics
# refactor: init() column list updated for the 10 renamed PK compartments
# (was NSAID_DEPOT/NSAID_C/SSZ_DEPOT/SSZ_C/MTX_DEPOT/MTX_C/TNFi_DEPOT/
# TNFi_C/IL17i_DEPOT/IL17i_C)
out_s3b <- mod %>%
  param(list(USE_HLAB27 = 1, k_abx_gi = 0, k_abx_chlamyd = 0)) %>%
  init(as.list(tail(out_s3, 1)[, c("PATH","CHRON_PATH","INNATE","NEUT","MACRO",
       "TH1","TH17","TREG","TNF_c","IL17_c","IL6_c","IL10_c","IFNg_c","SYNOV",
       "CARTDMG","PAIN","CRP_c","JCOUNT",
       "GUT_NSAID","CENT_NSAID","GUT_SSZ","CENT_SSZ","GUT_MTX","CENT_MTX",
       "GUT_TNFI","CENT_TNFI","GUT_IL17I","CENT_IL17I")])) %>%
  mrgsim(end = 8760 - 2016, delta = 6, carry_out = "time") %>%
  as.data.frame() %>%
  mutate(time = time + 2016)

out_s3 <- bind_rows(out_s3, out_s3b) %>%
  mutate(scenario = "S3: Antibiotics (Doxycycline 12wk) + NSAIDs", day = time / 24)

cat("Scenario 3 (ABX + NSAID). 24-wk CRP:",
    round(out_s3$CRP_obs[which.min(abs(out_s3$day - 168))], 1), "mg/L\n")

# ============================================================
# SCENARIO 4: NSAIDs + Sulfasalazine (Chronic/persistent ReA)
# ============================================================
# SSZ 500 mg BID × 2 wk → 1000 mg BID × 2 wk → 1500 mg BID maintenance (titration)
s4_params <- list(USE_HLAB27 = 1, k_abx_gi = 0, k_abx_chlamyd = 0)

# refactor: cmt strings updated for the renamed SSZ depot (was "SSZ_DEPOT")
ssz_ev_phase1 <- ev(cmt = "GUT_SSZ", time = seq(0, 330, by = 12), amt = 500)
ssz_ev_phase2 <- ev(cmt = "GUT_SSZ", time = seq(336, 666, by = 12), amt = 1000)
ssz_ev_phase3 <- ev(cmt = "GUT_SSZ", time = seq(672, 8760 - 12, by = 12), amt = 1500)
nsaid_s4      <- make_nsaid_events(0, 2016)

s4_ev <- c(ssz_ev_phase1, ssz_ev_phase2, ssz_ev_phase3, nsaid_s4)

out_s4 <- mod %>%
  param(s4_params) %>%
  mrgsim(ev = s4_ev, end = 8760, delta = 6, carry_out = "time") %>%
  as.data.frame() %>%
  mutate(scenario = "S4: NSAIDs + Sulfasalazine (2-3 g/day)", day = time / 24)

cat("Scenario 4 (NSAIDs + SSZ). 24-wk CRP:",
    round(out_s4$CRP_obs[which.min(abs(out_s4$day - 168))], 1),
    "mg/L; Joint count:", round(out_s4$swollen_jt[which.min(abs(out_s4$day - 168))], 1), "\n")

# ============================================================
# SCENARIO 5: TNF Inhibitor (Etanercept, refractory/HLA-B27+ axial ReA)
# ============================================================
# Start NSAID + SSZ; add ETN 50 mg/wk at week 12 (after SSZ trial)
s5_params <- list(USE_HLAB27 = 1, k_abx_gi = 0, k_abx_chlamyd = 0)

ssz_s5  <- make_ssz_events(0, 8760, dose_mg = 1000)
nsaid_s5 <- make_nsaid_events(0, 2016)
etn_ev   <- make_tnfi_events(start_hr = 2016, dose_mg = 50, interval_hr = 168, n_doses = 52)

s5_ev <- c(ssz_s5, nsaid_s5, etn_ev)

out_s5 <- mod %>%
  param(s5_params) %>%
  mrgsim(ev = s5_ev, end = 8760, delta = 6, carry_out = "time") %>%
  as.data.frame() %>%
  mutate(scenario = "S5: NSAIDs + SSZ → ETN (TNF-i, wk 12+)", day = time / 24)

cat("Scenario 5 (ETN). 52-wk CRP:",
    round(out_s5$CRP_obs[which.min(abs(out_s5$day - 364))], 1),
    "mg/L; VAS pain:", round(out_s5$VAS_pain[which.min(abs(out_s5$day - 364))], 1), "\n")

# ============================================================
# COMBINED DATA
# ============================================================
all_scenarios <- bind_rows(out_s1, out_s2, out_s3, out_s4, out_s5)

scenario_colors <- c(
  "S1: No Treatment"                              = "#8B0000",
  "S2: NSAIDs (Naproxen 500mg BID)"              = "#FF8C00",
  "S3: Antibiotics (Doxycycline 12wk) + NSAIDs"  = "#228B22",
  "S4: NSAIDs + Sulfasalazine (2-3 g/day)"       = "#4169E1",
  "S5: NSAIDs + SSZ → ETN (TNF-i, wk 12+)"      = "#9400D3"
)

# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

plot_cytokines <- function(data, ymax_TNF = 80, ymax_IL17 = 50, ymax_IL6 = 60) {
  p1 <- ggplot(data, aes(x = day, y = TNF_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    ylim(0, ymax_TNF) +
    labs(title = "TNF-α (pg/mL)", x = "Day", y = "pg/mL", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p2 <- ggplot(data, aes(x = day, y = IL17_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    ylim(0, ymax_IL17) +
    labs(title = "IL-17A (pg/mL)", x = "Day", y = "pg/mL", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p3 <- ggplot(data, aes(x = day, y = IL6_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    ylim(0, ymax_IL6) +
    labs(title = "IL-6 (pg/mL)", x = "Day", y = "pg/mL", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  (p1 | p2 | p3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

plot_clinical_endpoints <- function(data) {
  p1 <- ggplot(data, aes(x = day, y = VAS_pain, color = scenario)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = scenario_colors) +
    geom_hline(yintercept = c(30, 70), linetype = "dashed", color = "gray50") +
    annotate("text", x = 350, y = 30, label = "Mild (30)", size = 3, color = "gray40") +
    annotate("text", x = 350, y = 70, label = "Severe (70)", size = 3, color = "gray40") +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    ylim(0, 100) +
    labs(title = "VAS Pain Score (0–100)", x = "Day", y = "mm", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p2 <- ggplot(data, aes(x = day, y = CRP_obs, color = scenario)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = scenario_colors) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray50") +
    annotate("text", x = 350, y = 5.5, label = "ULN 5 mg/L", size = 3, color = "gray40") +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    labs(title = "CRP (mg/L)", x = "Day", y = "mg/L", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p3 <- ggplot(data, aes(x = day, y = swollen_jt, color = scenario)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = scenario_colors) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    scale_x_continuous(breaks = seq(0, 365, 52)) +
    ylim(0, 12) +
    labs(title = "Swollen Joint Count (0–28)", x = "Day", y = "Count", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  (p1 | p2 | p3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

plot_immune_cells <- function(data) {
  p1 <- ggplot(data, aes(x = day, y = TH17_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Th17 Cells (normalized)", x = "Day", y = "AU", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p2 <- ggplot(data, aes(x = day, y = TH1_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Th1 Cells (normalized)", x = "Day", y = "AU", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  p3 <- ggplot(data, aes(x = day, y = TREG_obs, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = scenario_colors) +
    labs(title = "Treg Cells (normalized)", x = "Day", y = "AU", color = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  (p1 | p2 | p3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

plot_drug_pk <- function(data) {
  s2_data <- filter(data, scenario == "S2: NSAIDs (Naproxen 500mg BID)")
  s4_data <- filter(data, scenario == "S4: NSAIDs + Sulfasalazine (2-3 g/day)")
  s5_data <- filter(data, scenario == "S5: NSAIDs + SSZ → ETN (TNF-i, wk 12+)")

  p1 <- ggplot(filter(s2_data, day <= 14), aes(x = day, y = NSAID_obs)) +
    geom_line(color = "#FF8C00", linewidth = 1.2) +
    geom_hline(yintercept = 1.8, linetype = "dashed", color = "red") +
    annotate("text", x = 12, y = 2.1, label = "IC50 COX-2", size = 3, color = "red") +
    labs(title = "NSAID Plasma (µg/mL)\nFirst 2 weeks", x = "Day", y = "µg/mL") +
    theme_bw(base_size = 11)

  p2 <- ggplot(filter(s4_data, day <= 90), aes(x = day, y = SSZ_obs)) +
    geom_line(color = "#4169E1", linewidth = 1.2) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
    annotate("text", x = 75, y = 11.5, label = "IC50", size = 3, color = "red") +
    labs(title = "Sulfasalazine Plasma (µg/mL)\nFirst 3 months", x = "Day", y = "µg/mL") +
    theme_bw(base_size = 11)

  p3 <- ggplot(filter(s5_data, day >= 84 & day <= 280), aes(x = day, y = TNFi_obs)) +
    geom_line(color = "#9400D3", linewidth = 1.2) +
    geom_hline(yintercept = 1.2, linetype = "dashed", color = "red") +
    annotate("text", x = 270, y = 1.5, label = "IC50", size = 3, color = "red") +
    labs(title = "Etanercept Plasma (µg/mL)\nWeeks 12–40", x = "Day", y = "µg/mL") +
    theme_bw(base_size = 11)

  (p1 | p2 | p3)
}

# ============================================================
# OUTCOME SUMMARY TABLE
# ============================================================
summary_timepoints <- c(84, 168, 364)  # 12 wk, 24 wk, 52 wk

outcome_table <- all_scenarios %>%
  filter(day %in% sapply(summary_timepoints, function(tp) {
    all_scenarios$day[which.min(abs(all_scenarios$day - tp))]
  })) %>%
  group_by(scenario, day) %>%
  summarise(
    `VAS Pain`    = round(mean(VAS_pain), 1),
    `CRP (mg/L)`  = round(mean(CRP_obs), 1),
    `Swollen Jts` = round(mean(swollen_jt), 1),
    `IL-17A (pg/mL)` = round(mean(IL17_obs), 1),
    `TNF-α (pg/mL)`  = round(mean(TNF_obs), 1),
    `Synovitis (AU)` = round(mean(SYNOV_obs), 2),
    `Cart Damage (%)` = round(mean(CARTDMG_pct), 2),
    .groups = "drop"
  ) %>%
  mutate(Week = round(day / 7))

print(outcome_table)

# ============================================================
# CALIBRATION CHECK (Kvien 1994 — reactive arthritis outcomes)
# ============================================================
# Reference: 80 patients, 6-month follow-up
# Untreated / NSAID only: ~50% resolution at 6 months (CRP normalization)
# Target at day 168 (24 wk): S1 CRP should still be elevated (>5 mg/L in ~50%)
cat("\n=== Calibration Check (Kvien 1994 targets) ===\n")
cat("S1 (No Tx) CRP at 24 wk:",
    round(out_s1$CRP_obs[which.min(abs(out_s1$day - 168))], 1),
    "mg/L [Target: >5 mg/L in ~50% untreated]\n")
cat("S2 (NSAID) CRP at 12 wk:",
    round(out_s2$CRP_obs[which.min(abs(out_s2$day - 84))], 1),
    "mg/L [Target: symptomatic relief but not full remission]\n")
cat("S5 (ETN) CRP at 24 wk:",
    round(out_s5$CRP_obs[which.min(abs(out_s5$day - 168))], 1),
    "mg/L [Target: <5 mg/L, Braun 2014 response ~70-80%]\n")

# ============================================================
# EXAMPLE PLOTS
# ============================================================
if (interactive()) {
  print(plot_cytokines(all_scenarios))
  print(plot_clinical_endpoints(all_scenarios))
  print(plot_immune_cells(all_scenarios))
  print(plot_drug_pk(all_scenarios))
}

cat("\n=== Reactive Arthritis QSP Model Ready ===\n")
cat("Use plot_cytokines(), plot_clinical_endpoints(),\n")
cat("    plot_immune_cells(), plot_drug_pk() to visualize.\n")
cat("Outcome summary stored in 'outcome_table'.\n")
