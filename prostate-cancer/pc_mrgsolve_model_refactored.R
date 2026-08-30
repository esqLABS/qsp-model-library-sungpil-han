################################################################################
# Prostate Cancer QSP Model (mrgsolve) -- PK/PD REFACTORED
# ============================================================
# Author  : QSP Disease Model Library (CCR)
# Date    : 2026-06-23 (original); refactored per FORK_WORKFLOW_GUIDE.md Part 2
#
# SCOPE
# -----
#   1. HPG Axis: GnRH -> LH -> Testosterone -> DHT
#   2. Androgen Receptor (AR) Signaling: AR activation, PSA production
#   3. Tumor Cell Kinetics: Proliferating + Quiescent cells, PSA
#   4. PI3K/AKT pathway (PTEN loss -> AKT -> mTOR)
#   5. Bone Metastasis: RANKL-mediated osteoclast/osteoblast dynamics, BMD
#   6. Drug PK/PD:
#      a) GnRH Agonists  (Leuprolide monthly depot)
#      b) GnRH Antagonist (Degarelix, Relugolix oral)
#      c) AR Pathway Inhibitors (Enzalutamide, Abiraterone)
#      d) Docetaxel chemotherapy
#      e) PARP inhibitor (Olaparib)
#      f) Bone agents (Denosumab)
#   7. Clinical Biomarkers: PSA, Testosterone, BMD, rPFS
#
# [refactor] Per FORK_WORKFLOW_GUIDE.md Part 2 (pluggable PK, named Hill
# interface). This file has EIGHT real, independently-dosed compounds with
# their own PK, but the census (driver-patches/data/compound_perturbation_
# census.md) had only ONE row for this file, "AR Signaling (DEG)" -- itself
# mislabeled (see pc_refactor_notes.md: "DEG" is Degarelix, not a
# degradation-rate parameter, but the census's process-description display
# name followed the same wrong-label-right-stem pattern documented in
# neonatal-hyperbilirubinemia's refactor). All eight were identified by
# reading the code and renamed to the fork's naming convention (C_<STEM>
# exposed concentration, EFFECT_<STEM> named disease effect):
#   Leuprolide (LEUP), Degarelix (DEG), Relugolix (REL), Enzalutamide (ENZ),
#   Abiraterone (ABI), Docetaxel (DOC), Olaparib (OLA), Denosumab (DEN).
# The census was corrected accordingly (8 rows now, one per compound).
#
# Renamed (values unchanged from the original except where noted):
#   Leup_depot/Leup_c/Flare_eff      -> GUT_LEUP/CENT_LEUP/FLARE_LEUP
#   kLeup_rel/kLeup_elim/V_Leup      -> KA_LEUP/KE_LEUP/V1_LEUP
#   GnRH_flare/flare_decay           -> FLARE_MULT_LEUP/FLARE_DECAY_LEUP
#   (hardcoded 0.97, 2.0)            -> EMAX_LEUP, EC50_LEUP (new named params)
#   Deg_sc/Deg_c                     -> GUT_DEG/CENT_DEG
#   kDeg_abs/kDeg_elim               -> KA_DEG/KE_DEG
#   Deg_Emax/Deg_EC50                -> EMAX_DEG/EC50_DEG
#   Rel_gut/Rel_c                    -> GUT_REL/CENT_REL
#   kRel_abs/F_Rel/kRel_elim/V_Rel   -> KA_REL/F_REL/KE_REL/V1_REL
#   (hardcoded 0.98)/Rel_EC50        -> EMAX_REL (new)/EC50_REL
#   Enz_gut/Enz_c                    -> GUT_ENZ/CENT_ENZ
#   kEnz_abs/F_Enz/kEnz_elim/V_Enz   -> KA_ENZ/F_ENZ/KE_ENZ/V1_ENZ
#   Enz_Emax/Enz_EC50                -> EMAX_ENZ/EC50_ENZ
#   Abi_gut/Abi_c                    -> GUT_ABI/CENT_ABI
#   kAbi_abs/F_Abi/kAbi_elim/V_Abi   -> KA_ABI/F_ABI/KE_ABI/V1_ABI
#   Abi_Emax/Abi_EC50/Abi_MW         -> EMAX_ABI/EC50_ABI/MW_ABI
#   Doc_c/Doc_p                      -> CENT_DOC/PERI_DOC (bespoke, see notes)
#   kDoc_elim1/k12/k21/elim2         -> KE1_DOC/K12_DOC/K21_DOC/KE2_DOC
#   Doc_Emax/Doc_EC50                -> EMAX_DOC/EC50_DOC
#   Ola_gut/Ola_c                    -> GUT_OLA/CENT_OLA
#   kOla_abs/F_Ola/kOla_elim/V_Ola   -> KA_OLA/F_OLA/KE_OLA/V1_OLA
#   Ola_Emax/Ola_EC50                -> EMAX_OLA/EC50_OLA
#   Den_sc/Den_c                     -> GUT_DEN/CENT_DEN
#   kDen_abs/kDen_elim/V_Den         -> KA_DEN/KE_DEN/V1_DEN
#   Den_Emax/Den_EC50                -> EMAX_DEN/EC50_DEN
#   Leup_suppress/Deg_blockade/Rel_blockade/Enz_AR_inh/Abi_T_inh/Doc_kill/
#   Ola_kill(minus HRR_def)/Den_RANKL_inh
#     -> EFFECT_LEUP/EFFECT_DEG/EFFECT_REL/EFFECT_ENZ/EFFECT_ABI/EFFECT_DOC/
#        EFFECT_OLA/EFFECT_DEN
# New (not in original; make explicit a Hill shape the original already had
# implicitly as a plain ratio): GAMMA_<STEM>=1 for all 8 compounds;
# EMAX_LEUP=0.97, EC50_LEUP=2.0, EMAX_REL=0.98 (were hardcoded literals).
#
# Bespoke deviations (disclosed, not archetype-forcing):
# - DOC keeps the original's dual-elimination 2-compartment shape (central
#   AND peripheral each eliminate) -- does not fit archetype 2's CL/Q/V
#   template (elimination only from central), so renamed in place instead
#   (KE1_DOC/K12_DOC/K21_DOC/KE2_DOC) rather than forced into CL/Q/V.
# - DEGARELIX'S OWN VOLUME IS WRONG IN THE ORIGINAL: V_Deg=1000.0 is
#   declared but never used; the original's own dxdt_Deg_c divides by
#   V_Den (Denosumab's volume, 3.0 L) instead. Preserved *numerically*
#   (V1_DEG=3.0, the value actually in effect) but decoupled *structurally*
#   from Denosumab's own V1_DEN (each compound now has its own independent
#   parameter, per the "not interleaved" plumbing requirement). See notes
#   and translations/UPSTREAM_ISSUES.md #67.
#
# Pre-existing upstream build defects (unrelated to any compound's own
# archetype choice), fixed syntax-only, disclosed, logged as
# translations/UPSTREAM_ISSUES.md #67:
#  1. 16 of 33 $INIT @annotated lines were missing the required description
#     field ("Leup_c : 0.0" with no third field) -- does not compile.
#  2. $CMT + $INIT jointly redeclared all 33 compartments (same family of
#     defect as #61/#63/#64/#65/#66) -- does not compile. Fixed by dropping
#     $INIT entirely and setting every compartment's initial value via
#     `<CMT>_0 = value;` in $MAIN instead (same values).
#  3. $CAPTURE duplicated all 18 of its own entries against $CMT -- does
#     not compile. Fixed by dropping the duplicated names ($CAPTURE now
#     lists only the 16 new C_<STEM>/EFFECT_<STEM> derived quantities;
#     mrgsolve reports every compartment's own state via
#     /model_manifest's outputPaths regardless of $CAPTURE).
#  4. Leuprolide's own GnRH-agonist-effect guard referenced an undeclared
#     `self.trt_leup` -- does not compile. Handled as an in-scope Leuprolide
#     design decision (not a generic defect): the guard is dropped and the
#     unconditional flare/suppression formula is used directly, which is
#     (a) what the file's own PK never actually got to run given defects
#     1-3 above, (b) what the surrounding comment describes ("flare then
#     desensitize"), and (c) numerically neutral whenever Leuprolide is
#     absent (Leup_c=Flare_eff=0 already gives the guarded branch's value,
#     1.0), confirmed by this file's own untreated-scenario verification.
# The checked-in original (pc_mrgsolve_model.R) is untouched and still
# carries all four defects exactly as written.
#
# All disease-side content (HPG axis, AR signaling, tumor kinetics,
# PI3K/AKT, bone metastasis) is untouched apart from the renamed drug
# variables it reads. Full rationale, per-compound archetype, and the
# floating-point-scale verification results (7 of the original's own
# scenarios + 3 constructed single-dose checks for the 3 compounds no
# shipped scenario ever doses) are in pc_refactor_notes.md.
#
# KEY CLINICAL PARAMETERS CALIBRATED TO:
#   - Testosterone nadir: <50 ng/dL (castrate) within 4 wk of ADT initiation
#   - PSA response: ~90% decline from baseline at 3 months with ADT
#   - Enzalutamide OS benefit: ~4 months in mCRPC (AFFIRM trial)
#   - Abiraterone PSA response rate: ~29% in mCRPC (COU-AA-301)
#   - Docetaxel: 3.0-month OS benefit in mCRPC (TAX 327)
#   - BRCA2-mutant: Olaparib ORR 33% (PROfound trial)
#   - Radium-223: 3.6-month OS benefit (ALSYMPCA trial)
#
# TREATMENT SCENARIOS (run with event tables)
#   1. Untreated / Natural History
#   2. ADT Alone (Leuprolide 7.5 mg IM monthly)
#   3. ADT + Enzalutamide (ARPI doublet)
#   4. ADT + Abiraterone (CYP17A1 inhibitor)
#   5. Docetaxel 75 mg/m^2 q3w x 6 cycles
#   6. Olaparib (BRCA2-mutant mCRPC)
#   7. Sequential ADT -> ARPI -> Docetaxel
#
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ==============================================================================
# 1. MODEL CODE (mrgsolve C++ ODE block)
# ==============================================================================

code <- '
$PROB
Prostate Cancer QSP Model (refactored: pluggable PK, named Hill interface)
Compartments: HPG axis, AR signaling, tumor cells, bone metastasis, drug PK
Refactor: renamed all 8 real compounds PK to the forks naming convention
(C_<STEM>, EFFECT_<STEM>, ...). See pc_refactor_notes.md for detail.

$PARAM @annotated
// ---- HPG Axis ----
kLH_prod   : 8.0   : LH basal production rate (IU/L/day)
kLH_deg    : 0.96  : LH degradation rate (/day, t1/2~17h)
GnRH_base  : 1.0   : Baseline GnRH tone (normalized)
kT_prod    : 0.25  : Testosterone production rate (nmol/L/day per LH)
kT_deg     : 0.35  : Testosterone degradation (/day, t1/2~2h)
T_adrenal  : 0.05  : Adrenal androgen contribution to T (nmol/L/day)
f5alpha    : 0.10  : Fraction T -> DHT by 5-alpha-reductase
kDHT_deg   : 0.50  : DHT degradation rate (/day)
T_baseline : 15.0  : Baseline testosterone (nmol/L, ~432 ng/dL)
LH_base    : 5.0   : Baseline LH (IU/L)

// ---- AR Signaling ----
kAR_synth  : 0.05  : AR protein synthesis rate (nmol/cell/day)
kAR_deg    : 0.05  : AR protein basal degradation (/day)
kon_AR     : 2.0   : DHT-AR binding on-rate (1/nmol/day)
koff_AR    : 0.5   : DHT-AR dissociation rate (/day)
k_nuc      : 1.5   : AR-DHT nuclear translocation rate (/day)
k_nuc_off  : 0.3   : AR nuclear export rate (/day)
kAR_nuc_deg: 0.2   : Nuclear AR degradation (/day)
PSA_kprod  : 0.002 : PSA production rate per nuclear AR-cell (ng/mL/nmol/10^9cells/day)
kPSA_deg   : 0.10  : PSA degradation rate (/day, t1/2~7 days)
AR0        : 1.0   : Initial AR protein (normalized units)

// ---- Tumor Cell Kinetics ----
k_prolif   : 0.06  : Tumor proliferation rate (/day, doubling ~12d)
k_death    : 0.01  : Basal tumor cell death rate (/day)
k_quiesce  : 0.02  : Entry into quiescence (/day)
k_unquiesce: 0.015 : Exit from quiescence (/day)
k_death_q  : 0.005 : Quiescent cell death (/day)
TC0        : 1.0   : Initial tumor cell burden (normalized = 1)
TC_cap     : 1000.0: Tumor cell carrying capacity (normalized)
AR_prolif_EC50: 0.5: AR nuclear occupancy for half-max proliferation (normalized)

// ---- PI3K/AKT Pathway (PTEN-loss phenotype) ----
PTEN_loss  : 0.7   : Fraction of PTEN loss (0=normal, 1=complete loss)
kAKT_base  : 0.3   : Baseline AKT activity
kAKT_max   : 1.0   : Max AKT activity (PTEN-null)
k_AKT_AR   : 0.3   : AKT enhancement of AR activity
k_AKT_BCL2 : 0.2   : AKT enhancement of BCL2 (anti-apoptotic)

// ---- Bone Metastasis ----
kOC_form   : 0.05  : Osteoclast formation rate
kOC_deg    : 0.15  : Osteoclast degradation (/day)
kOB_form   : 0.04  : Osteoblast formation rate
kOB_deg    : 0.12  : Osteoblast degradation (/day)
kRANKL     : 0.8   : RANKL-driven osteoclast activation
kOPG       : 0.4   : OPG inhibition of RANKL
kBMD_form  : 0.003 : BMD formation rate by OB (/day)
kBMD_resorb: 0.005 : BMD resorption rate by OC (/day)
BMD0       : 1.0   : Baseline BMD (T-score normalized)
k_bonehom  : 0.02  : Tumor bone-homing rate
BSI0       : 0.0   : Initial bone scan index

// ---- Leuprolide (LEUP) PK: GnRH agonist, 7.5 mg monthly IM depot ----
KA_LEUP    : 0.033 : Depot release rate (/day, ~21d sustained) [was kLeup_rel]
KE_LEUP    : 0.693 : Elimination rate constant (/day, t1/2~1h iv, depot sustained) [was kLeup_elim]
V1_LEUP    : 40.0  : Volume of distribution (L) [was V_Leup]
FLARE_MULT_LEUP : 3.0 : Initial GnRH flare multiplier (first 7 days) [was GnRH_flare]
FLARE_DECAY_LEUP: 0.5 : Flare decay rate (/day) [was flare_decay]
EMAX_LEUP  : 0.97  : Max GnRH suppression by leuprolide [was hardcoded 0.97]
EC50_LEUP  : 2.0   : EC50 for GnRH suppression (ng/mL) [was hardcoded 2.0]
GAMMA_LEUP : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Degarelix (DEG) PK: GnRH antagonist, 240 mg SC loading ----
KA_DEG     : 0.15  : SC absorption rate (/day) [was kDeg_abs]
KE_DEG     : 0.023 : Elimination rate constant (/day, t1/2~28d) [was kDeg_elim]
// NOTE: the originals own V_Deg=1000.0 is declared but never referenced in
// $ODE; the originals Deg_c equation actually divides by V_Den (Denosumabs
// volume, 3.0 L) instead -- an apparent copy-paste error. V1_DEG=3.0 below
// reproduces the originals *actual* numeric behaviour, now under
// Degarelixs own name (decoupled from Denosumabs parameter). See notes.
V1_DEG     : 3.0   : Effective volume of distribution (L) [see NOTE above; originals V_Deg was dead]
EMAX_DEG   : 0.98  : Max GnRH-R blockade by degarelix [was Deg_Emax]
EC50_DEG   : 0.001 : EC50 for GnRH-R blockade (ug/mL) [was Deg_EC50]
GAMMA_DEG  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Relugolix (REL) PK: GnRH antagonist, 120 mg QD oral ----
KA_REL     : 1.4   : Oral absorption rate (/day) [was kRel_abs]
F_REL      : 0.12  : Oral bioavailability [was F_Rel]
KE_REL     : 1.0   : Elimination rate constant (/day, t1/2~16.5h) [was kRel_elim]
V1_REL     : 2800.0: Volume of distribution (L) [was V_Rel]
EMAX_REL   : 0.98  : Max GnRH-R blockade by relugolix [was hardcoded 0.98]
EC50_REL   : 0.005 : EC50 for GnRH-R blockade (ng/mL) [was Rel_EC50]
GAMMA_REL  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Enzalutamide (ENZ) PK: AR pathway inhibitor, 160 mg QD oral ----
KA_ENZ     : 1.5   : Absorption rate (/day) [was kEnz_abs]
F_ENZ      : 0.84  : Oral bioavailability [was F_Enz]
KE_ENZ     : 0.114 : Elimination rate constant (/day, t1/2~5.8d) [was kEnz_elim]
V1_ENZ     : 110.0 : Volume of distribution (L/kg x 70 kg) [was V_Enz]
EMAX_ENZ   : 0.95  : Max AR inhibition by enzalutamide [was Enz_Emax]
EC50_ENZ   : 3.0   : EC50 for AR blockade (uM) [was Enz_EC50]
GAMMA_ENZ  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Abiraterone (ABI) PK: CYP17A1 inhibitor, 1000 mg QD + prednisone ----
KA_ABI     : 0.8   : Absorption rate (/day) [was kAbi_abs]
F_ABI      : 0.10  : Oral bioavailability (fasted) [was F_Abi]
KE_ABI     : 1.7   : Elimination rate constant (/day, t1/2~10h) [was kAbi_elim]
V1_ABI     : 19669.0: Volume of distribution (L) [was V_Abi]
MW_ABI     : 391.6 : Abiraterone molecular weight [was Abi_MW]
EMAX_ABI   : 0.95  : Max CYP17A1 inhibition [was Abi_Emax]
EC50_ABI   : 0.05  : EC50 for CYP17A1 inhibition (uM) [was Abi_EC50]
GAMMA_ABI  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Docetaxel (DOC) PK: cytotoxic chemo, 75 mg/m^2 IV q3w (bespoke 2-cmt) ----
KE1_DOC    : 3.94  : Central-compartment elimination rate (/day) [was kDoc_elim1]
K12_DOC    : 1.5   : Distribution to peripheral compartment (/day) [was kDoc_k12]
K21_DOC    : 0.8   : Return from peripheral (/day) [was kDoc_k21]
KE2_DOC    : 0.231 : Peripheral-compartment elimination rate (/day, t1/2~3d) [was kDoc_elim2]
EMAX_DOC   : 0.90  : Max tumor cell kill by docetaxel [was Doc_Emax]
EC50_DOC   : 0.05  : EC50 for cytotoxicity (uM) [was Doc_EC50]
GAMMA_DOC  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Olaparib (OLA) PK: PARP inhibitor, 300 mg BID oral ----
KA_OLA     : 1.4   : Absorption rate (/day) [was kOla_abs]
F_OLA      : 0.66  : Bioavailability [was F_Ola]
KE_OLA     : 1.6   : Elimination rate constant (/day, t1/2~11h) [was kOla_elim]
V1_OLA     : 167.0 : Volume of distribution (L) [was V_Ola]
EMAX_OLA   : 0.80  : Max kill in HRR-deficient (BRCA2-mut) [was Ola_Emax]
EC50_OLA   : 0.1   : EC50 for PARP inhibition (uM) [was Ola_EC50]
GAMMA_OLA  : 1.0   : Hill exponent [no explicit exponent in original]
HRR_def    : 0.0   : HRR deficiency status (0=proficient, 1=deficient); gates Olaparibs effect downstream, not part of its own Hill term

// ---- Denosumab (DEN) PK: RANKL antibody, 120 mg SC q4w ----
KA_DEN     : 0.062 : SC absorption rate (/day, t1/2 of absorption ~8d) [was kDen_abs]
KE_DEN     : 0.023 : Elimination rate constant (/day, t1/2~28d) [was kDen_elim]
V1_DEN     : 3.0   : Volume of distribution (L) [was V_Den]
EMAX_DEN   : 0.95  : Max RANKL inhibition [was Den_Emax]
EC50_DEN   : 0.5   : EC50 for RANKL neutralization (ug/mL) [was Den_EC50]
GAMMA_DEN  : 1.0   : Hill exponent [no explicit exponent in original]

// ---- Disease Progression Parameters ----
k_CRPC     : 0.003 : Rate of acquiring CRPC resistance (/day)
AR_v7_time : 365.0 : Time to ARv7 emergence (days, in CRPC)
k_ARv7     : 0.002 : Rate of ARv7 emergence under ARPI pressure

$CMT @annotated
// HPG Axis
LH        : Luteinizing hormone (IU/L)
T         : Testosterone (nmol/L)
DHT       : Dihydrotestosterone (nmol/L)

// AR Signaling
AR_free   : Free AR protein (nmol/cell, normalized)
AR_DHT    : AR-DHT cytoplasmic complex
AR_nuc    : Nuclear AR-DHT complex
PSA       : Serum PSA (ng/mL)

// Tumor Cell Kinetics
TC_p      : Proliferating tumor cells (normalized)
TC_q      : Quiescent tumor cells (normalized)
CRPC_frac : Fraction of castration-resistant cells (0-1)
ARv7_frac : Fraction of ARv7-positive cells (0-1)

// PI3K/AKT
AKT_act   : Active AKT (normalized 0-1)

// Bone Metastasis
OC        : Osteoclasts (normalized)
OB        : Osteoblasts (normalized)
BMD       : Bone mineral density (normalized T-score)
BoneMets  : Bone metastasis burden (normalized)

// Drug PK - LEUP (Leuprolide)
GUT_LEUP  : Leuprolide depot (mg)
CENT_LEUP : Leuprolide central (mg, amount)
FLARE_LEUP: GnRH flare effect

// Drug PK - DEG (Degarelix)
GUT_DEG   : Degarelix SC depot (mg)
CENT_DEG  : Degarelix central (mg, amount)

// Drug PK - REL (Relugolix)
GUT_REL   : Relugolix gut (mg)
CENT_REL  : Relugolix central (mg, amount)

// Drug PK - ENZ (Enzalutamide)
GUT_ENZ   : Enzalutamide gut (mg)
CENT_ENZ  : Enzalutamide central (mg, amount)

// Drug PK - ABI (Abiraterone)
GUT_ABI   : Abiraterone gut (mg)
CENT_ABI  : Abiraterone central (mg, amount)

// Drug PK - DOC (Docetaxel)
CENT_DOC  : Docetaxel central conc (uM) [concentration-state, dosed directly]
PERI_DOC  : Docetaxel peripheral conc (uM)

// Drug PK - OLA (Olaparib)
GUT_OLA   : Olaparib gut (mg)
CENT_OLA  : Olaparib central (mg, amount)

// Drug PK - DEN (Denosumab)
GUT_DEN   : Denosumab SC depot (mg)
CENT_DEN  : Denosumab central (mg, amount)

$MAIN
// ---- Initial conditions (replaces the originals separate $INIT block,
// which jointly redeclared these 33 compartments alongside $CMT -- a
// pre-existing mrgsolve 2.0.1 build defect; see pc_refactor_notes.md) ----
LH_0 = 5.0;
T_0 = 15.0;
DHT_0 = 1.5;
AR_free_0 = 1.0;
AR_DHT_0 = 0.5;
AR_nuc_0 = 0.3;
PSA_0 = 4.0;
TC_p_0 = 1.0;
TC_q_0 = 0.2;
CRPC_frac_0 = 0.01;
ARv7_frac_0 = 0.0;
AKT_act_0 = 0.0;
OC_0 = 1.0;
OB_0 = 1.0;
BMD_0 = 1.0;
BoneMets_0 = 0.0;
GUT_LEUP_0 = 0.0;
CENT_LEUP_0 = 0.0;
FLARE_LEUP_0 = 0.0;
GUT_DEG_0 = 0.0;
CENT_DEG_0 = 0.0;
GUT_REL_0 = 0.0;
CENT_REL_0 = 0.0;
GUT_ENZ_0 = 0.0;
CENT_ENZ_0 = 0.0;
GUT_ABI_0 = 0.0;
CENT_ABI_0 = 0.0;
CENT_DOC_0 = 0.0;
PERI_DOC_0 = 0.0;
GUT_OLA_0 = 0.0;
CENT_OLA_0 = 0.0;
GUT_DEN_0 = 0.0;
CENT_DEN_0 = 0.0;

// ==============================
// COMPOUND CONCENTRATIONS & NAMED HILL EFFECTS
// (computed here in $MAIN, matching the *originals own* placement of the
// analogous Leup_suppress/Enz_AR_inh/Doc_kill/etc quantities -- mrgsolve
// evaluates $MAIN once per output record, before that records own $ODE
// integration runs, so these values are piecewise-constant across each
// output interval, updated once per interval, exactly as in the original.
// This was verified deliberately: moving this block into $ODE instead
// (continuous re-evaluation) changes the integrated trajectory measurably
// -- confirmed by a direct comparison run, see pc_refactor_notes.md --
// because the originals own dynamics genuinely depend on this once-per-
// interval update cadence, not merely on a display/capture timing nuance.
// Reproducing the originals placement (not "improving" it) is therefore
// the numerically-faithful choice per this guides verification mandate.
// ==============================

// LEUP: depot + central, linear (archetype 3 minus peripheral)
double C_LEUP = CENT_LEUP / V1_LEUP * 1000.0;
double EFFECT_LEUP = EMAX_LEUP * C_LEUP / (EC50_LEUP + C_LEUP);
// Flare term is a separate PD sub-effect tied to its own state (FLARE_LEUP),
// not a function of C_LEUP -- kept as a disclosed bespoke combination, not
// folded into EFFECT_LEUP itself. The originals own guard here referenced
// an undeclared `self.trt_leup`, which does not compile; removed (see notes).
double GnRH_agonist_effect = (1.0 + FLARE_MULT_LEUP * FLARE_LEUP) * (1.0 - EFFECT_LEUP);

// DEG: depot + central, linear (archetype 3 minus peripheral)
double C_DEG = CENT_DEG / V1_DEG;
double EFFECT_DEG = EMAX_DEG * C_DEG / (EC50_DEG + C_DEG);

// REL: depot + central, linear (archetype 3 minus peripheral)
double C_REL = CENT_REL / V1_REL * 1e6;
double EFFECT_REL = EMAX_REL * C_REL / (EC50_REL + C_REL);

double GnRH_total = GnRH_base * (1.0 - EFFECT_DEG) * (1.0 - EFFECT_REL) * GnRH_agonist_effect;

// ENZ: depot + central, linear (archetype 3 minus peripheral)
double C_ENZ = CENT_ENZ / V1_ENZ;
double EFFECT_ENZ = EMAX_ENZ * C_ENZ / (EC50_ENZ + C_ENZ);

// ABI: depot + central, linear (archetype 3 minus peripheral)
double C_ABI = CENT_ABI / V1_ABI * 1e6 / MW_ABI;
double EFFECT_ABI = EMAX_ABI * C_ABI / (EC50_ABI + C_ABI);

// DOC: bespoke 2-compartment (elimination from BOTH central and peripheral;
// does not fit archetype 2s CL/Q/V shape, which eliminates only centrally).
// Dosed directly into CENT_DOC (concentration-state, not an amount depot) --
// preserved as-is per the originals own dosing convention.
double C_DOC = CENT_DOC;
double EFFECT_DOC = EMAX_DOC * C_DOC / (EC50_DOC + C_DOC);

// OLA: depot + central, linear (archetype 3 minus peripheral)
double C_OLA = CENT_OLA / V1_OLA;
double EFFECT_OLA = EMAX_OLA * C_OLA / (EC50_OLA + C_OLA);

// DEN: depot + central, linear (archetype 3 minus peripheral)
double C_DEN = CENT_DEN / V1_DEN;
double EFFECT_DEN = EMAX_DEN * C_DEN / (EC50_DEN + C_DEN);

// ---- AKT Activity (PTEN-loss model) ----
double AKT_ss = kAKT_base + (kAKT_max - kAKT_base) * PTEN_loss * TC_p / (TC_p + 0.5);

// ---- AR Nuclear Occupancy (effective, accounting for CRPC mechanisms) ----
double AR_nuc_eff = AR_nuc * (1.0 - EFFECT_ENZ) + ARv7_frac * 0.5;
double AR_nuc_norm = AR_nuc_eff / (AR_nuc_eff + 0.3);

// ---- Tumor Cell Proliferation Rate ----
double prolif_AR  = AR_nuc_norm;
double prolif_AKT = k_AKT_AR * AKT_act;
double prolif_eff = fmax(0.0, prolif_AR + prolif_AKT);
double k_prolif_eff = k_prolif * prolif_eff;

// Apoptosis enhanced by drug treatment (Olaparib gated by HRR_def, as original)
double k_death_eff = k_death * (1.0 + EFFECT_DOC + HRR_def * EFFECT_OLA)
                     * (1.0 + k_AKT_BCL2 * (1.0 - AKT_act));

// ---- RANKL Signaling for Bone ----
double RANKL_eff = kRANKL * BoneMets * (1.0 - EFFECT_DEN) /
                   (1.0 + kOPG * OB);

$ODE
// ==============================
// HPG AXIS
// ==============================
double LH_prod = kLH_prod * GnRH_total;
double LH_deg  = kLH_deg * LH;
dxdt_LH = LH_prod - LH_deg;

double T_prod = kT_prod * LH * (1.0 - EFFECT_ABI) + T_adrenal;
double T_deg  = kT_deg * T;
double neg_fb = T / (T + T_baseline);   // negative feedback on GnRH
dxdt_T = T_prod - T_deg;

double DHT_synth = f5alpha * T;
double DHT_deg   = kDHT_deg * DHT;
dxdt_DHT = DHT_synth - DHT_deg;

// ==============================
// AR SIGNALING
// ==============================
double AR_free_synth = kAR_synth;
double AR_free_deg   = kAR_deg * AR_free;
double AR_bind       = kon_AR * DHT * AR_free * (1.0 - EFFECT_ENZ);
double AR_unbind     = koff_AR * AR_DHT;
dxdt_AR_free = AR_free_synth - AR_free_deg - AR_bind + AR_unbind;

double AR_nuc_in  = k_nuc * AR_DHT;
double AR_nuc_out = k_nuc_off * AR_nuc + kAR_nuc_deg * AR_nuc;
dxdt_AR_DHT = AR_bind - AR_unbind - AR_nuc_in;
dxdt_AR_nuc = AR_nuc_in - AR_nuc_out;

// PSA production proportional to nuclear AR x total tumor cells
double TC_total   = TC_p + TC_q;
double PSA_prod   = PSA_kprod * AR_nuc_eff * TC_total;
double PSA_elim   = kPSA_deg * PSA;
dxdt_PSA = PSA_prod - PSA_elim;

// ==============================
// AKT ACTIVITY (quasi-equilibrium)
// ==============================
dxdt_AKT_act = 5.0 * (AKT_ss - AKT_act);  // fast equilibration

// ==============================
// TUMOR CELL KINETICS
// ==============================
double TC_total_now = TC_p + TC_q;
double logistic     = 1.0 - TC_total_now / TC_cap;

// CRPC subpopulation can proliferate despite ADT
double CRPC_prolif = CRPC_frac * k_prolif * 0.8;  // AR-independent proliferation
double net_prolif  = (k_prolif_eff + CRPC_prolif) * logistic;

dxdt_TC_p = net_prolif * TC_p
             - k_death_eff * TC_p
             - k_quiesce * TC_p
             + k_unquiesce * TC_q;

dxdt_TC_q = k_quiesce * TC_p
             - k_unquiesce * TC_q
             - k_death_q * TC_q * (1.0 + EFFECT_DOC * 0.5);

// CRPC fraction growth: accelerated by ARPI pressure
double CRPC_growth = k_CRPC * (1.0 - CRPC_frac)
                     * (1.0 + EFFECT_ENZ * 2.0 + EFFECT_ABI * 1.5);
dxdt_CRPC_frac = CRPC_growth;

// ARv7 emergence: accelerated by ARPI pressure in CRPC context
double ARv7_growth = CRPC_frac * k_ARv7 * (1.0 - ARv7_frac)
                     * (1.0 + EFFECT_ENZ * 3.0);
dxdt_ARv7_frac = ARv7_growth;

// ==============================
// BONE METASTASIS
// ==============================
// Bone homing driven by CXCL12/CXCR4
dxdt_BoneMets = k_bonehom * TC_p * (1.0 - BoneMets / 10.0);

// Osteoclast dynamics (RANKL-driven by bone mets)
dxdt_OC = kOC_form * (1.0 + RANKL_eff) - kOC_deg * OC;

// Osteoblast dynamics (ET-1, Wnt from tumor cells in sclerotic mets)
double ET1_eff = BoneMets * 0.5;  // endothelin-1 osteoblast stimulation
dxdt_OB = kOB_form * (1.0 + ET1_eff) - kOB_deg * OB;

// BMD dynamics
dxdt_BMD = kBMD_form * OB - kBMD_resorb * OC;

// ==============================
// DRUG PK - LEUP (Leuprolide, GnRH agonist)
// ==============================
dxdt_GUT_LEUP   = -KA_LEUP * GUT_LEUP;
dxdt_CENT_LEUP  =  KA_LEUP * GUT_LEUP - KE_LEUP * CENT_LEUP;
dxdt_FLARE_LEUP = -FLARE_DECAY_LEUP * FLARE_LEUP;

// ==============================
// DRUG PK - DEG (Degarelix, GnRH antagonist)
// ==============================
dxdt_GUT_DEG  = -KA_DEG * GUT_DEG;
dxdt_CENT_DEG =  KA_DEG * GUT_DEG - KE_DEG * CENT_DEG;

// ==============================
// DRUG PK - REL (Relugolix, oral GnRH antagonist)
// ==============================
dxdt_GUT_REL  = -KA_REL * GUT_REL;
dxdt_CENT_REL =  KA_REL * F_REL * GUT_REL - KE_REL * CENT_REL;

// ==============================
// DRUG PK - ENZ (Enzalutamide, oral ARPI)
// ==============================
dxdt_GUT_ENZ  = -KA_ENZ * GUT_ENZ;
dxdt_CENT_ENZ =  KA_ENZ * F_ENZ * GUT_ENZ - KE_ENZ * CENT_ENZ;

// ==============================
// DRUG PK - ABI (Abiraterone, oral CYP17A1 inhibitor)
// ==============================
dxdt_GUT_ABI  = -KA_ABI * GUT_ABI;
dxdt_CENT_ABI =  KA_ABI * F_ABI * GUT_ABI - KE_ABI * CENT_ABI;

// ==============================
// DRUG PK - DOC (Docetaxel, 2-compartment IV, bespoke dual elimination)
// ==============================
dxdt_CENT_DOC = -(KE1_DOC + K12_DOC) * CENT_DOC + K21_DOC * PERI_DOC;
dxdt_PERI_DOC =  K12_DOC * CENT_DOC - (K21_DOC + KE2_DOC) * PERI_DOC;

// ==============================
// DRUG PK - OLA (Olaparib, oral PARP inhibitor)
// ==============================
dxdt_GUT_OLA  = -KA_OLA * GUT_OLA;
dxdt_CENT_OLA =  KA_OLA * F_OLA * GUT_OLA - KE_OLA * CENT_OLA;

// ==============================
// DRUG PK - DEN (Denosumab, SC RANKL antibody)
// ==============================
dxdt_GUT_DEN  = -KA_DEN * GUT_DEN;
dxdt_CENT_DEN =  KA_DEN * GUT_DEN - KE_DEN * CENT_DEN;

$CAPTURE @annotated
C_LEUP     : Leuprolide plasma conc (ng/mL) [was Leup_c]
C_DEG      : Degarelix plasma conc (ug/mL) [was Deg_c]
C_REL      : Relugolix plasma conc (ng/mL) [was Rel_c]
C_ENZ      : Enzalutamide plasma conc (uM, nominal) [was Enz_c]
C_ABI      : Abiraterone plasma conc (uM) [was Abi_c]
C_DOC      : Docetaxel central conc (uM) [was Doc_c]
C_OLA      : Olaparib plasma conc (uM, nominal) [was Ola_c]
C_DEN      : Denosumab plasma conc (ug/mL) [was Den_c]
EFFECT_LEUP: Leuprolide suppression fraction (Hill ratio component)
EFFECT_DEG : Degarelix GnRH-R blockade fraction
EFFECT_REL : Relugolix GnRH-R blockade fraction
EFFECT_ENZ : Enzalutamide AR-blockade fraction
EFFECT_ABI : Abiraterone CYP17A1-inhibition fraction
EFFECT_DOC : Docetaxel tumor-kill fraction
EFFECT_OLA : Olaparib PARP-inhibition fraction (pre HRR_def gating)
EFFECT_DEN : Denosumab RANKL-inhibition fraction
'

# ==============================================================================
# 2. COMPILE MODEL
# ==============================================================================

mod <- mcode("ProstateCancer_QSP_refactored", code)
cat("Model compiled successfully.\n")
cat("Compartments:", length(init(mod)), "\n")
cat("Parameters:", length(param(mod)), "\n")

# ==============================================================================
# 3. DEFINE TREATMENT SCENARIOS
# ==============================================================================

# Time frame: 3 years (1095 days)
end_time <- 1095
delta    <- 1  # daily output

# Helper: convert mg to uM for IV bolus
mg_to_uM_docetaxel <- function(mg, V_L = 6.0, MW = 861.9) {
  (mg / MW / V_L) * 1e6  # uM
}

# --- Scenario 1: Untreated (natural history) ---
sc1_events <- ev()  # no treatment
sc1_name   <- "Untreated"

# --- Scenario 2: ADT Alone (Leuprolide 7.5 mg IM monthly) ---
# [renamed cmt: Leup_depot -> GUT_LEUP]
sc2_events <- ev(time = seq(0, 1080, by = 28),  # monthly
                 cmt  = "GUT_LEUP",
                 amt  = 7.5,
                 evid = 1) %>%
  mutate(Flare_eff = ifelse(time == 0, 1.0, 0.0)) %>%
  filter(TRUE)

# Simplified: just dose leuprolide depot
sc2_events <- ev(
  data.frame(
    time = seq(0, 1080, by = 28),
    cmt  = "GUT_LEUP",
    amt  = 7.5,
    evid = 1
  )
)
sc2_name <- "ADT (Leuprolide)"

# --- Scenario 3: ADT + Enzalutamide (160 mg QD) ---
# Start ADT at day 0, add enzalutamide at day 0 (upfront combination)
# [renamed cmt: Enz_gut -> GUT_ENZ]
enz_daily <- ev(
  data.frame(
    time = seq(0, 1094),
    cmt  = "GUT_ENZ",
    amt  = 160,
    evid = 1
  )
)
sc3_events <- c(sc2_events, enz_daily)
sc3_name   <- "ADT + Enzalutamide"

# --- Scenario 4: ADT + Abiraterone (1000 mg QD) ---
# [renamed cmt: Abi_gut -> GUT_ABI]
abi_daily <- ev(
  data.frame(
    time = seq(0, 1094),
    cmt  = "GUT_ABI",
    amt  = 1000,
    evid = 1
  )
)
sc4_events <- c(sc2_events, abi_daily)
sc4_name   <- "ADT + Abiraterone"

# --- Scenario 5: Docetaxel 75 mg/m^2 q3w x 6 cycles ---
# Start at day 0 (chemo-naive mCRPC)
# [renamed cmt: Doc_c -> CENT_DOC]
doc_cycles <- ev(
  data.frame(
    time = seq(0, 5 * 21, by = 21),  # 6 cycles
    cmt  = "CENT_DOC",
    amt  = mg_to_uM_docetaxel(135),  # 75 mg/m^2 x 1.8 m^2
    evid = 1
  )
)
sc5_events <- c(sc2_events, doc_cycles)
sc5_name   <- "ADT + Docetaxel (x6)"

# --- Scenario 6: Olaparib (HRR-deficient mCRPC, 300 mg BID) ---
# [renamed cmt: Ola_gut -> GUT_OLA]
ola_bid <- ev(
  data.frame(
    time = c(outer(seq(0, 1094), c(0, 0.5), "+")),
    cmt  = "GUT_OLA",
    amt  = 300,
    evid = 1
  )
)
sc6_events <- c(sc2_events, ola_bid)
sc6_name   <- "ADT + Olaparib (HRR-def)"

# --- Scenario 7: Sequential ADT -> ARPI -> Docetaxel ---
# ADT: day 0-365, add enzalutamide at day 180 (CRPC transition),
# add docetaxel at day 450 (docetaxel-switch after ARPI failure)
# [renamed cmt: Enz_gut -> GUT_ENZ, Doc_c -> CENT_DOC]
enz_from180 <- ev(
  data.frame(
    time = seq(180, 449),
    cmt  = "GUT_ENZ",
    amt  = 160,
    evid = 1
  )
)
doc_from450 <- ev(
  data.frame(
    time = seq(450, 450 + 5 * 21, by = 21),
    cmt  = "CENT_DOC",
    amt  = mg_to_uM_docetaxel(135),
    evid = 1
  )
)
sc7_events <- c(sc2_events, enz_from180, doc_from450)
sc7_name   <- "Sequential ADT->ARPI->Docetaxel"

# ==============================================================================
# 4. RUN SIMULATIONS
# ==============================================================================

run_scenario <- function(events, name, pars = list()) {
  sim_mod <- mod
  if (length(pars) > 0) sim_mod <- param(sim_mod, .x = pars)

  out <- mrgsim(sim_mod,
                events = events,
                end    = end_time,
                delta  = delta,
                carry_out = "evid") %>%
    as_tibble() %>%
    mutate(Scenario = name)
  return(out)
}

# Run all scenarios
results_list <- list(
  run_scenario(sc1_events, sc1_name),
  run_scenario(sc2_events, sc2_name),
  run_scenario(sc3_events, sc3_name),
  run_scenario(sc4_events, sc4_name),
  run_scenario(sc5_events, sc5_name),
  run_scenario(sc6_events, sc6_name, pars = list(HRR_def = 1.0)),
  run_scenario(sc7_events, sc7_name)
)

results_all <- bind_rows(results_list)

# ==============================================================================
# 5. DEFINE CLINICAL ENDPOINTS
# ==============================================================================

clinical_endpoints <- results_all %>%
  group_by(Scenario) %>%
  summarise(
    # PSA endpoints
    PSA_baseline    = PSA[time == 0][1],
    PSA_nadir       = min(PSA, na.rm = TRUE),
    PSA_nadir_time  = time[which.min(PSA)],
    PSA50_response  = any(PSA < (PSA[time == 0][1] * 0.5), na.rm = TRUE),
    PSA_doubling    = {
      psa_vals <- PSA[time > 180 & PSA > 0]
      t_vals   <- time[time > 180 & PSA > 0]
      if (length(psa_vals) > 5) {
        fit <- lm(log(psa_vals) ~ t_vals)
        log(2) / max(coef(fit)[2], 1e-6)
      } else { NA_real_ }
    },
    # Testosterone endpoints
    T_nadir         = min(T, na.rm = TRUE),
    T_castrate      = any(T < 1.73, na.rm = TRUE),  # 50 ng/dL = 1.73 nmol/L
    # Tumor endpoints
    TC_max          = max(TC_p + TC_q, na.rm = TRUE),
    TC_final        = (TC_p + TC_q)[time == max(time)][1],
    # Bone endpoints
    BMD_final       = BMD[time == max(time)][1],
    BoneMets_final  = BoneMets[time == max(time)][1],
    # Resistance
    CRPC_frac_final = CRPC_frac[time == max(time)][1],
    ARv7_frac_final = ARv7_frac[time == max(time)][1],
    .groups = "drop"
  )

cat("\n=== CLINICAL ENDPOINTS SUMMARY ===\n")
print(clinical_endpoints %>% select(Scenario, PSA_nadir, PSA_nadir_time,
                                     PSA50_response, T_castrate,
                                     TC_final, BMD_final, CRPC_frac_final))

# ==============================================================================
# 6. SENSITIVITY ANALYSIS
# ==============================================================================

# PSA response sensitivity to PTEN loss and HRR deficiency
sa_params <- expand.grid(
  PTEN_loss = seq(0, 1, by = 0.25),
  HRR_def   = c(0, 1)
)

sa_results <- mapply(function(pten, hrr) {
  run_scenario(sc6_events,
               paste0("PTEN=", pten, " HRR=", hrr),
               pars = list(PTEN_loss = pten, HRR_def = hrr)) %>%
    filter(time %in% c(0, 90, 180, 365)) %>%
    mutate(PTEN_loss = pten, HRR_def = hrr)
}, sa_params$PTEN_loss, sa_params$HRR_def, SIMPLIFY = FALSE) %>%
  bind_rows()

# ==============================================================================
# 7. PLOT RESULTS
# ==============================================================================

theme_qsp <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#2c3e50"),
    strip.text       = element_text(color = "white", face = "bold"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pal7 <- c("#e74c3c","#3498db","#2ecc71","#f39c12",
          "#9b59b6","#1abc9c","#e67e22")

# Plot 1: PSA over time by scenario
p1 <- ggplot(results_all, aes(x = time / 30.4, y = PSA,
                               color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 4.0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "gray70") +
  scale_y_log10(limits = c(0.01, 1000),
                breaks  = c(0.1, 1, 10, 100),
                labels  = c("0.1", "1", "10", "100")) +
  scale_color_manual(values = pal7) +
  labs(title = "PSA Dynamics Under Different Treatment Scenarios",
       subtitle = "Dashed: 4 ng/mL upper normal; Dotted: 0.2 ng/mL (deep response)",
       x = "Time (months)", y = "PSA (ng/mL, log scale)",
       color = "Treatment") +
  theme_qsp
print(p1)

# Plot 2: Testosterone over time
p2 <- results_all %>%
  filter(Scenario %in% c("Untreated", "ADT (Leuprolide)",
                          "ADT + Enzalutamide", "ADT + Abiraterone")) %>%
  ggplot(aes(x = time / 30.4, y = T * 28.84,  # nmol/L -> ng/dL
             color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  annotate("text", x = 1, y = 55, label = "Castrate threshold (50 ng/dL)",
           hjust = 0, color = "red", size = 3) +
  scale_color_manual(values = pal7[1:4]) +
  labs(title = "Testosterone Suppression: GnRH Agents",
       x = "Time (months)", y = "Testosterone (ng/dL)",
       color = "Treatment") +
  theme_qsp
print(p2)

# Plot 3: Tumor Burden (TC_p + TC_q)
p3 <- results_all %>%
  mutate(TC_total = TC_p + TC_q) %>%
  ggplot(aes(x = time / 30.4, y = TC_total,
             color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = pal7) +
  labs(title = "Total Tumor Cell Burden Over Time",
       x = "Time (months)", y = "Tumor Cell Burden (normalized)",
       color = "Treatment") +
  theme_qsp
print(p3)

# Plot 4: Resistance Mechanisms
p4 <- results_all %>%
  filter(Scenario %in% c("ADT (Leuprolide)",
                          "ADT + Enzalutamide",
                          "Sequential ADT->ARPI->Docetaxel")) %>%
  pivot_longer(cols = c(CRPC_frac, ARv7_frac),
               names_to = "Mechanism", values_to = "Fraction") %>%
  mutate(Mechanism = recode(Mechanism,
    CRPC_frac = "CRPC Subpopulation",
    ARv7_frac = "ARv7-positive Cells")) %>%
  ggplot(aes(x = time / 30.4, y = Fraction * 100,
             color = Scenario, linetype = Mechanism)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = pal7[c(2, 3, 7)]) +
  labs(title = "Emergence of Resistance Mechanisms",
       subtitle = "CRPC = AR-independent growth; ARv7 = ligand-independent AR splice variant",
       x = "Time (months)", y = "Resistant Cell Fraction (%)",
       color = "Treatment", linetype = "Mechanism") +
  theme_qsp
print(p4)

# Plot 5: Bone Metastasis & BMD
p5 <- results_all %>%
  filter(Scenario %in% c("Untreated", "ADT (Leuprolide)",
                          "ADT + Enzalutamide")) %>%
  pivot_longer(cols = c(BoneMets, BMD),
               names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = recode(Variable,
    BoneMets = "Bone Metastasis Burden",
    BMD      = "Bone Mineral Density (normalized)")) %>%
  ggplot(aes(x = time / 30.4, y = Value,
             color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = pal7[1:3]) +
  labs(title = "Bone Metastasis and BMD Dynamics",
       x = "Time (months)", y = "Value",
       color = "Treatment") +
  theme_qsp
print(p5)

# Plot 6: AR Signaling (Nuclear AR, PSA)
p6 <- results_all %>%
  filter(Scenario %in% c("Untreated", "ADT (Leuprolide)",
                          "ADT + Enzalutamide",
                          "ADT + Abiraterone")) %>%
  pivot_longer(cols = c(AR_nuc, PSA),
               names_to = "Variable", values_to = "Value") %>%
  mutate(Variable = recode(Variable,
    AR_nuc = "Nuclear AR (normalized)",
    PSA    = "PSA (ng/mL)")) %>%
  ggplot(aes(x = time / 30.4, y = Value,
             color = Scenario, group = Scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = pal7[1:4]) +
  labs(title = "AR Signaling Suppression by ADT and ARPI",
       x = "Time (months)", y = "Value",
       color = "Treatment") +
  theme_qsp
print(p6)

# Plot 7: Sensitivity Analysis (PSA response at 3 months)
p7 <- sa_results %>%
  filter(time == 90) %>%
  ggplot(aes(x = PTEN_loss, y = PSA,
             color = factor(HRR_def),
             group = interaction(PTEN_loss, HRR_def))) +
  geom_point(size = 3) +
  geom_line(aes(group = factor(HRR_def)), linewidth = 0.8) +
  scale_color_manual(values = c("#2ecc71", "#e74c3c"),
                     labels = c("HRR Proficient", "HRR Deficient")) +
  labs(title = "Sensitivity: PSA at 3 months vs PTEN Loss Severity",
       subtitle = "Under ADT + Olaparib scenario",
       x = "PTEN Loss Fraction (0=Normal, 1=Complete Loss)",
       y = "PSA (ng/mL)",
       color = "HRR Status") +
  theme_qsp
print(p7)

# ==============================================================================
# 8. SUMMARY REPORT
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PROSTATE CANCER QSP MODEL (REFACTORED) - SIMULATION REPORT\n")
cat(strrep("=", 70), "\n\n")
cat("Model scope  : HPG axis + AR signaling + Tumor kinetics + Bone mets\n")
cat("Scenarios    : 7 treatment regimens\n")
cat("Time horizon : 3 years (1095 days)\n\n")
cat("Key clinical calibration targets:\n")
cat("  - Castrate testosterone (<50 ng/dL): ADT achieves within 4 weeks\n")
cat("  - PSA >=50% decline: Expected in ADT-sensitive disease\n")
cat("  - CRPC emergence: ~12-24 months under ADT alone\n")
cat("  - ARv7: Emerges under prolonged ARPI pressure\n")
cat("  - BMD decline: Progressive under ADT without bone agents\n\n")
cat("Resistance mechanisms modeled:\n")
cat("  - CRPC_frac: AR-independent (AR bypass, PI3K/AKT)\n")
cat("  - ARv7_frac: Ligand-independent AR splice variant\n")
cat("  - PTEN_loss: Activates AKT -> compensates for AR inhibition\n\n")
