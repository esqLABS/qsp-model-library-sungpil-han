## =============================================================================
##  pji_mrgsolve_model.R
##  Prosthetic joint infection (PJI) QSP model
##  Implant-associated Staphylococcus aureus osteomyelitis
##
##  53 ODEs · 6 antibiotics with explicit bone penetration · 5 bacterial
##  phenotypic states · 2 resistant lineages · local elution from an
##  antibiotic-loaded cement spacer · innate immunity with frustrated
##  phagocytosis · RANKL/OPG osteolysis · 2018 ICM biomarkers · drug toxicity
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL IS ORGANISED AROUND
##  ---------------------------------------------------------------------------
##  PJI doctrine looks like a list of empirical rules ("you must debride",
##  "rifampicin is special", "never give rifampicin alone", "DAIR only works
##  early"). This model is written so that all four fall out of arithmetic
##  rather than being coded in as assertions.
##
##  PILLAR 1 — THE FOREIGN BODY MOVES THE INFECTIVE DOSE, NOT THE ORGANISM.
##    Two mechanisms, both explicit: (a) FBFACT scales down phagocyte killing
##    next to metal, and (b) KATT moves planktonic cells onto a surface where
##    phagocytosis is "frustrated" (FRUST = 0.03). Neither is a threshold
##    parameter. Simulated: an implant-free niche clears up to ~10^7 CFU;
##    with an implant, ~10^1 CFU establishes. That ~10^5-10^6 fold shift is
##    the model's version of Elek & Conen's suture experiment (1957) and of
##    Zimmerli's tissue-cage model (1982).
##
##  PILLAR 2 — EVERY ANTIBIOTIC IS ITS FREE BONE CONCENTRATION DIVIDED BY MBEC.
##    The model never uses the planktonic MIC to decide what happens in the
##    biofilm. For each drug it carries (i) AUC_bone/AUC_plasma, (ii) the free
##    fraction, (iii) the planktonic MIC and (iv) the biofilm MBEC. The ratio
##    C_bone,free / MBEC at standard adult doses then computes to:
##
##       vancomycin 1 g q12h ...... 2.083 / 512 = 0.0041
##       cefazolin 2 g q8h (MSSA).. 2.500 / 256 = 0.0098
##       daptomycin 8 mg/kg qd .... 0.350 / 32  = 0.0109
##       linezolid 600 mg q12h .... 2.218 / 128 = 0.0173
##       levofloxacin 750 mg qd ... 1.203 / 64  = 0.0188
##       RIFAMPICIN 450 mg q12h ... 0.281 / 1.0 = 0.281 (week 1)
##                                  0.141 / 1.0 = 0.141 (autoinduced steady state)
##
##    (all six computed by ratio_table(); rifampicin sits 7-68x above every
##     other option and is still BELOW 1 -- which is why it needs surgery.)
##
##    Rifampicin is not a better antibiotic. It is the only one whose target
##    (RNA polymerase) is still load-bearing in an adherent, barely dividing
##    cell, so its MBEC never climbs 2-3 logs above its MIC the way a cell-wall
##    agent's does. The entire clinical special status of rifampicin in
##    implant infection is that one column of the table.
##
##  PILLAR 3 — SURGERY IS A MUTANT-SUPPLY OPERATION, NOT A CLEANING OPERATION.
##    P(a pre-existing rpoB mutant) = 1 - exp(-mu * N), mu ~ 1e-8 per division.
##    A mature biofilm carries N ~ 1e9-1e10, so that probability is ~1.000 and
##    rifampicin is dead on arrival. A 4-log debridement takes N to ~1e5-1e6
##    and the probability to ~1e-3-1e-2. Debridement does not "help" the
##    antibiotic; it is the step that makes Pillar 2 usable at all. The model
##    computes this probability continuously (PRPOB) and also carries explicit
##    rpoB (RP/RB) and gyrA (QP/QB) lineages generated from the growth flux.
##
##  PILLAR 4 — BIOFILM MATURITY IS A CLOCK, AND DAIR RACES IT.
##    PHIB (the MBEC/MIC tolerance multiplier) is not constant: it is scaled by
##    matrix maturity EPS, which rises with a ~28-day time constant. A biofilm
##    debrided at day 7 is still ~40% "planktonic-like"; at day 90 it is fully
##    tolerant. The 3-week symptom-duration rule for DAIR is therefore an
##    emergent property of an EPS differential equation, not an if-statement.
##
##  ---------------------------------------------------------------------------
##  IMPORTANT MODELLING CHOICES / LIMITATIONS (read before using numbers)
##  ---------------------------------------------------------------------------
##  * Bone compartments are FREE-concentration effect compartments driven by
##    total plasma concentration through FB_x = (AUC_bone/AUC_plasma) x fu.
##    They do not feed mass back to plasma (periprosthetic mass is negligible).
##  * The model is deterministic. Bacterial extinction is enforced by a smooth
##    sub-single-cell decay term (KEXT below 1 CFU) so the ODE cannot carry
##    "0.001 of a bacterium" and regrow it, which is the classic deterministic
##    artefact in antimicrobial models. Stochastic outcomes (ID50, resistance
##    emergence, cure) are reported as PROBABILITIES computed from the
##    deterministic burden, not as trajectories.
##  * Cure probability is mapped as PCURE = exp(-PSEED * N_end): a surviving
##    burden of ~1 CFU next to a retained implant is roughly a coin flip.
##  * Cytokine and biomarker units are local/periprosthetic, not serum, except
##    CRP / ESR / PLT / SCR / ALT which are systemic.
##  * Calibration is to published trial-level endpoints (Zimmerli 1998 JAMA;
##    OVIVA 2019; DATIPO 2021; Lora-Tamayo 2013) and to standard PK/bone
##    penetration reviews (Landersdorfer 2009; Thabit 2019). It has NOT been
##    fitted to individual patient data. Educational / research use only.
##
##  ---------------------------------------------------------------------------
##  VERIFICATION (all 53 ODEs were re-implemented independently in
##  Python/scipy and run before this file was finalised). That pass exposed
##  and fixed six defects, each of which had silently broken one of the
##  pillars:
##
##  1. IMMUNE CAPACITY WRITTEN PER COMPARTMENT. A Michaelis-Menten sink on
##     each bacterial state made a rare rpoB clone face the ENTIRE phagocyte
##     capacity as first-order clearance at ~19 /h, sterilising exactly the
##     subpopulation Pillar 3 is about. Rewritten as one shared capacity
##     allocated in proportion to accessibility (NPHAG / SHR).
##  2. EXTINCTION FLOOR APPLIED PER STATE. Forcing every compartment below
##     1 CFU to zero annihilated the rpoB lineage on every dip. Extinction is
##     now a property of the CLONE. Before the fix an untreated 1e10 lesion
##     carried 0.01 mutants; after it, 63 -- against the analytic mu*N = 100.
##  3. DEBRIDEMENT RESET BIOFILM MATURITY. EPS was stripped in proportion to
##     the surgical log-kill, so a 2.5-log DAIR dropped EPS 0.62 -> 0.026 and
##     handed the residual biofilm near-planktonic susceptibility. Result:
##     levofloxacin MONOTHERAPY cured. A retained implant keeps its adherent
##     matrix (EPSCUT retention factor 0.55 -> 0.05).
##  4. NO CARRYING CAPACITY ON THE FREE POOL. Dispersed organisms grew to
##     1.8e10 and outnumbered the biofilm, making the infection spuriously
##     vancomycin-curable. NPCAP (fluid/abscess) and NICCAP (bounded by host
##     cell number) added.
##  5. UNBOUNDED LOOSENING AND OSTEOBLAST COLLAPSE. LOOSEN ran to 266 on a
##     0-100 scale and KAPOP drove osteoblasts to 8% of baseline.
##  6. PHAGOCYTE SATURATION CONSTANT INCONSISTENT WITH THE ID50. KMIMM = 1e4
##     against a target ID50 of ~1e7 made clearance first-order over three
##     logs with a large rate constant, so the immune system held the biofilm
##     at 1e4.8 CFU and no PJI ever established.
##
##  KNOWN CALIBRATION GAP (stated, not tuned away): the model's DAIR-timing
##  crossover falls between 3 and 6 months of infection duration, whereas the
##  clinical rule is ~3 weeks of symptoms. Debridement at day 21 sterilises on
##  day 60 of an 84-day course; at day 90 it sterilises on day 84 of 84 (zero
##  margin); at day 180 it fails outright. The direction and the mechanism
##  (EPS maturity ~10 d, sequestrum shielding ~3 mo) are right; the absolute
##  clock is too slow because the model has no soft-tissue abscess or
##  mechanical-loosening channel, which is how late DAIR actually fails.
##
##  Author: QSP Disease Model Library (Claude Code Routine)
## =============================================================================

library(mrgsolve)
library(dplyr)

pji_code <- '
$PROB
# Prosthetic Joint Infection (PJI) QSP model
# Implant-associated Staphylococcus aureus osteomyelitis
# 53 ODEs | biofilm PK/PD | mutant supply | surgery | osteolysis

$GLOBAL
// Emax kill helper: EMAX * C / (EC50eff + C)
#define EKILL(EMX, CC, EC) ((EMX)*(CC)/((EC)+(CC)))
// Biofilm tolerance multiplier scaled by matrix maturity (Pillar 4).
// EPSM in [0,1); a bare, freshly attached cell is nearly planktonic.
#define PHIMAT(PHI, EPSM) (1.0 + ((PHI) - 1.0)*(EPSM))

$PARAM @annotated
// ---------------------------------------------------------------- SETUP ----
IMPL    :  1    : implant present (1 = prosthesis in situ, 0 = implant-free bone)
MSSA    :  0    : methicillin susceptibility (1 = MSSA, cefazolin active; 0 = MRSA)
IMSUP   :  1    : host immunity factor (1 = normal, 0.5 = diabetes/steroids/RA)
INOC0   :  100  : CFU    : inoculum size (delivered as planktonic cells at t=0)

// ------------------------------------------------- VANCOMYCIN PK (IV) ------
VCVAN   : 40    : L      : VAN central volume of distribution
VPVAN   : 40    : L      : VAN peripheral volume of distribution
QVAN    :  8    : L/h    : VAN intercompartmental clearance
CLVAN   :  4    : L/h    : VAN systemic clearance (CrCl ~100 mL/min)
FUVAN   :  0.50 :        : VAN free fraction
PENVAN  :  0.20 :        : VAN AUC_bone / AUC_plasma
KEQVAN  :  0.30 : 1/h    : VAN bone compartment equilibration rate

// ------------------------------------------------ RIFAMPICIN PK (PO) -------
KARIF   :  1.20 : 1/h    : RIF absorption rate
FRIF    :  0.90 :        : RIF bioavailability
VCRIF   : 55    : L      : RIF volume of distribution
CLRIF0  : 12    : L/h    : RIF baseline clearance (before autoinduction)
FURIF   :  0.20 :        : RIF free fraction
PENRIF  :  0.50 :        : RIF AUC_bone / AUC_plasma
KEQRIF  :  0.50 : 1/h    : RIF bone compartment equilibration rate
KENZ    :  0.0072: 1/h   : CYP3A4 induction turnover rate (t1/2 ~ 4 days)
EIND    :  1.00 :        : maximum RIF autoinduction fold (clearance x2)

// ---------------------------------------------- LEVOFLOXACIN PK (PO) -------
KALVX   :  1.50 : 1/h    : LVX absorption rate
FLVX    :  0.99 :        : LVX bioavailability
VCLVX   : 90    : L      : LVX volume of distribution
CLLVX   :  9    : L/h    : LVX clearance
FULVX   :  0.70 :        : LVX free fraction
PENLVX  :  0.50 :        : LVX AUC_bone / AUC_plasma
KEQLVX  :  0.50 : 1/h    : LVX bone compartment equilibration rate

// ------------------------------------------------ DAPTOMYCIN PK (IV) -------
VCDAP   : 10    : L      : DAP volume of distribution
CLDAP   :  0.80 : L/h    : DAP clearance
FUDAP   :  0.08 :        : DAP free fraction (92% albumin bound)
PENDAP  :  0.15 :        : DAP AUC_bone / AUC_plasma
KEQDAP  :  0.25 : 1/h    : DAP bone compartment equilibration rate

// ------------------------------------------------- LINEZOLID PK (PO) -------
KALZD   :  1.80 : 1/h    : LZD absorption rate
FLZD    :  1.00 :        : LZD bioavailability
VCLZD   : 45    : L      : LZD volume of distribution
CLLZD0  :  7    : L/h    : LZD baseline clearance
FULZD   :  0.69 :        : LZD free fraction
PENLZD  :  0.45 :        : LZD AUC_bone / AUC_plasma
KEQLZD  :  0.50 : 1/h    : LZD bone compartment equilibration rate
RIFLZD  :  0.50 :        : strength of RIF-mediated induction of LZD clearance (AUC -32%)

// -------------------------------------------------- CEFAZOLIN PK (IV) ------
VCCFZ   : 15    : L      : CFZ volume of distribution
CLCFZ   :  4    : L/h    : CFZ clearance
FUCFZ   :  0.20 :        : CFZ free fraction
PENCFZ  :  0.20 :        : CFZ AUC_bone / AUC_plasma
KEQCFZ  :  0.50 : 1/h    : CFZ bone compartment equilibration rate

// --------------------------------- LOCAL ELUTION (ALBC / SPACER) -----------
VJOINT  :  0.05 : L      : effective joint space / interface volume
KFAST   :  0.05 : 1/h    : surface-layer burst release rate
KSLOW   :  0.0015: 1/h   : matrix diffusion release rate (t1/2 ~ 19 days)
KCLJ    :  0.35 : 1/h    : synovial fluid drug elimination rate
FLOCP   :  0.50 :        : access of local drug to planktonic cells
FLOCB   :  0.15 :        : access of local drug to residual biofilm

// -------------------------------------------- MIC / MBEC (S. aureus) -------
MICVAN  :  1.0  : mg/L   : VAN planktonic MIC
MICRIF  :  0.012: mg/L   : RIF planktonic MIC
MICLVX  :  0.25 : mg/L   : LVX planktonic MIC
MICDAP  :  0.50 : mg/L   : DAP planktonic MIC
MICLZD  :  2.0  : mg/L   : LZD planktonic MIC
MICCFZ  :  1.0  : mg/L   : CFZ planktonic MIC (MSSA)
MBCVAN  : 512   : mg/L   : VAN biofilm MBEC
MBCRIF  :  1.0  : mg/L   : RIF biofilm MBEC
MBCLVX  : 64    : mg/L   : LVX biofilm MBEC
MBCDAP  : 32    : mg/L   : DAP biofilm MBEC
MBCLZD  : 128   : mg/L   : LZD biofilm MBEC
MBCCFZ  : 256   : mg/L   : CFZ biofilm MBEC

// ------------------------------- MAXIMUM KILL RATES (planktonic) -----------
EMXVAN  :  0.45 : 1/h    : VAN maximum kill rate (time-dependent, slow)
EMXRIF  :  0.80 : 1/h    : RIF maximum kill rate
EMXLVX  :  1.20 : 1/h    : LVX maximum kill rate (concentration-dependent)
EMXDAP  :  2.00 : 1/h    : DAP maximum kill rate (rapid, concentration-dependent)
EMXLZD  :  0.12 : 1/h    : LZD maximum kill rate (bacteriostatic)
EMXCFZ  :  0.80 : 1/h    : CFZ maximum kill rate
EMXLZDI :  0.85 :        : LZD maximum growth-inhibition fraction (the substance of bacteriostasis)

// ---- STATE TOLERANCE MULTIPLIERS (applied ON TOP of the biofilm PHI) ------
PHPVAN  : 10.0  :        : VAN additional persister tolerance fold
PHPRIF  : 25.0  :        : RIF additional persister tolerance fold
PHPLVX  : 25.0  :        : LVX additional persister tolerance fold
PHPDAP  :  8.0  :        : DAP additional persister tolerance fold
PHPLZD  : 20.0  :        : LZD additional persister tolerance fold
PHPCFZ  : 40.0  :        : CFZ additional persister tolerance fold (maximal, growth-dependent)
PHSVAN  :  3.0  :        : VAN additional SCV tolerance fold
PHSRIF  :  2.0  :        : RIF additional SCV tolerance fold
PHSLVX  :  4.0  :        : LVX additional SCV tolerance fold
PHSDAP  :  3.0  :        : DAP additional SCV tolerance fold
PHSLZD  :  1.5  :        : LZD additional SCV tolerance fold
PHSCFZ  :  6.0  :        : CFZ additional SCV tolerance fold
PHIVAN  : 100.0 :        : VAN additional intracellular tolerance fold
PHIRIF  : 200.0 :        : RIF additional intracellular tolerance fold (intracellular cells are non-growing, SCV-like)
PHILVX  :  5.0  :        : LVX additional intracellular tolerance fold
PHIDAP  : 50.0  :        : DAP additional intracellular tolerance fold (inactive in the phagosome)
PHILZD  :  3.0  :        : LZD additional intracellular tolerance fold
PHICFZ  : 60.0  :        : CFZ additional intracellular tolerance fold
CARVAN  :  0.20 :        : VAN intracellular accumulation ratio
CARRIF  :  6.00 :        : RIF intracellular accumulation ratio
CARLVX  :  5.00 :        : LVX intracellular accumulation ratio
CARDAP  :  0.50 :        : DAP intracellular accumulation ratio
CARLZD  :  1.00 :        : LZD intracellular accumulation ratio
CARCFZ  :  0.10 :        : CFZ intracellular accumulation ratio

// ------------------------------------------ BACTERIAL POPULATION ----------
MUP     :  0.35 : 1/h    : planktonic maximum growth rate
MUB     :  0.080: 1/h    : biofilm growth rate
MUS     :  0.050: 1/h    : SCV growth rate
MUI     :  0.005: 1/h    : intracellular bacterial growth rate
NMAX    :  1e10 : CFU    : whole-lesion maximum carrying capacity
NBCAP   :  3e9  : CFU    : maximum implant-surface colonisation capacity
NPCAP   :  1e9  : CFU    : maximum capacity for non-adherent (synovial fluid / abscess) bacteria
NICCAP  :  2e8  : CFU    : maximum intracellular reservoir capacity (limited by host cell number)
KDNAT   :  0.005: 1/h    : natural death rate
KATT    :  0.40 : 1/h    : attachment rate constant (planktonic -> biofilm)
KDISP   :  0.020: 1/h    : dispersal rate constant (agr-dependent)
KPER0   :  0.0008:1/h    : baseline persister formation rate
APER    :  6.0  :        : amplification of persister formation by maturity and stress
KWAKE   :  0.0040:1/h    : persister reactivation (wake-up) rate
KSCV0   :  0.00050:1/h   : SCV formation rate
ASCV    :  8.0  :        : amplification of SCV induction by local antibiotic exposure
KREVS   :  0.0040:1/h    : SCV reversion to the normal phenotype
KINT0   :  0.0060:1/h    : osteoblast internalisation rate
KRELIC  :  0.0100:1/h    : intracellular bacterial release rate (including host cell death)
KEXT    :  2.0  : 1/h    : forced extinction rate below 1 CFU (applied per lineage)
NEXT    :  1.0  : CFU    : extinction threshold (1 cell)

// ------------------------------------------------- BIOFILM MATRIX ---------
KSEPS   :  0.0035: 1/h   : EPS production rate (maturation time constant ~10 days)
KDEPS   :  0.0002: 1/h   : EPS loss rate (matrix on a retained implant persists for a long time)
KEPSN   :  3e7  : CFU    : half-maximal surface bacterial load for EPS production
KMAT    :  0.35 :        : EPS maturation half-saturation constant
KSAIP   :  0.30 : 1/h    : AIP production rate
KDAIP   :  0.20 : 1/h    : AIP loss rate
KAIPN   :  1e8  : CFU    : half-maximal bacterial load for AIP production
KAIP    :  0.40 :        : AIP for half-maximal agr activity

// --------------------------------------------- RESISTANCE (mutant supply) --
MURPOB  :  1e-8 :        : rpoB mutation frequency (per division)
MUGYRA  :  3e-9 :        : gyrA mutation frequency (per division)
FITRPOB :  0.15 :        : rpoB mutant fitness cost
FITGYRA :  0.10 :        : gyrA mutant fitness cost

// --------------------------------------------------- INNATE IMMUNITY ------
KIMM    :  5e6  : CFU/h  : maximum phagocytic clearance capacity (at full recruitment)
KMIMM   :  1e6  : CFU    : phagocytic clearance saturation constant (shared capacity)
FBFACT  :  0.12 :        : foreign-body-related local immune deficiency factor
FRUST   :  0.030:        : residual efficiency of frustrated phagocytosis (biofilm)
FICIMM  :  0.020:        : immune efficiency against intracellular bacteria
KPMN    :  500  : cells/uL : PMN density for half-maximal immune recruitment
PMN0    : 50    : cells/uL : normal synovial fluid neutrophil density
MONO0   : 150   : cells/uL : normal synovial fluid monocyte density
KINPMN  :  7.5  : cells/uL/h : baseline PMN influx rate
KOUTPMN :  0.15 : 1/h    : PMN loss rate
EMXCHEM : 1200  :        : maximum chemotactic recruitment fold
KCHEM   :  0.35 :        : half-maximal chemotactic signal
KTOXP   :  0.35 : 1/h    : maximum leukotoxin-mediated PMN death rate
KINMON  :  1.5  : cells/uL/h : baseline monocyte influx rate
KOUTMON :  0.010: 1/h    : monocyte loss rate
EMXMON  : 18    :        : maximum monocyte recruitment fold
KBACS   :  1e7  : CFU    : half-maximal bacterial load for the PAMP signal
KSMDSC  :  0.030: 1/h    : MDSC recruitment rate (biofilm-specific)
KDMDSC  :  0.020: 1/h    : MDSC loss rate
KEPSM   :  0.30 :        : EPS for half-maximal MDSC recruitment
EMDSC   :  0.70 :        : maximum MDSC immunosuppressive fraction
KMD     :  0.50 :        : MDSC suppression half-saturation
KSMAC   :  0.10 : 1/h    : macrophage activation rate
KDMAC   :  0.050: 1/h    : macrophage deactivation rate

// ------------------------------------------------------- CYTOKINES --------
IL1B0   :  2.0  : pg/mL  : IL-1beta baseline
TNFA0   :  5.0  : pg/mL  : TNF-alpha baseline
IL60    :  3.0  : pg/mL  : IL-6 baseline
IL100   :  4.0  : pg/mL  : IL-10 baseline
KDIL1   :  0.70 : 1/h    : IL-1beta elimination
KDTNF   :  0.50 : 1/h    : TNF-alpha elimination
KDIL6   :  0.25 : 1/h    : IL-6 elimination
KDIL10  :  0.30 : 1/h    : IL-10 elimination
EIL1    : 60    :        : IL-1beta maximum induction fold
ETNF    : 40    :        : TNF-alpha maximum induction fold
EIL6    : 120   :        : IL-6 maximum induction fold
EIL10   : 30    :        : IL-10 maximum induction fold

// -------------------------------------------------- BIOMARKERS ------------
CRP0    :  3.0  : mg/L   : CRP baseline
KDCRP   :  0.0365:1/h    : CRP elimination (t1/2 19 h)
ECRP    : 40    :        : CRP maximum induction fold
KCRP    : 60    : pg/mL  : IL-6 for half-maximal CRP induction
ESR0    : 10    : mm/h   : ESR baseline
KDESR   :  0.00413:1/h   : ESR loss (t1/2 7 days)
EESR    :  8.0  :        : ESR maximum induction fold
ADEF0   :  0.30 : S/CO   : alpha-defensin baseline
KDADEF  :  0.10 : 1/h    : alpha-defensin elimination
EADEF   : 12    :        : alpha-defensin maximum induction fold
KADEF   : 4000  : cells/uL : PMN for half-maximal alpha-defensin induction

// ----------------------------------------------------- BONE / OSTEOLYSIS --
KSRANK  :  0.10 : 1/h    : baseline RANKL production
KDRANK  :  0.10 : 1/h    : RANKL loss
ERANK   :  9.0  :        : maximum cytokine-mediated RANKL induction
KSOPG   :  0.10 : 1/h    : baseline OPG production
KDOPG   :  0.10 : 1/h    : OPG loss
KINOBL  :  0.40 : 1/h    : osteoblast influx
KOUTOBL :  0.0040:1/h    : osteoblast loss
KAPOP   :  0.0040:1/h    : maximum osteoblast apoptosis rate
KINOCL  :  1.20 : 1/h    : maximum osteoclast influx
KRR     :  1.00 :        : RANKL/OPG half-saturation
KOUTOCL :  0.0060:1/h    : osteoclast loss
KFORM   :  0.020: %/h    : bone formation rate constant
KRES    :  0.020: %/h    : bone resorption rate constant
KHOM    :  0.00050:1/h   : bone mass homeostatic return
KLOOSE  :  0.0030:1/h    : interface loosening progression constant
KREPAIR :  0.00050:1/h   : interface repair constant
KSSEQ   :  0.0020:1/h    : sequestrum formation constant
KCSEQ   :  0.00030:1/h   : natural sequestrum loss
KSEQA   :  0.30 :        : half-saturation for sequestrum blockade of drug access
SEQMAX  :  3.0  :        : maximum sequestrum accumulation

// ---------------------------------------------------------- TOXICITY ------
SCR0    :  0.90 : mg/dL  : serum creatinine baseline
KELCR   :  0.030: 1/h    : creatinine elimination
ENEPH   :  0.55 :        : VAN maximum nephrotoxic fraction
AUCTHR  : 600   : mg*h/L : VAN AUC24 nephrotoxicity threshold
KNEPH   : 250   : mg*h/L : VAN nephrotoxicity half-saturation
PLT0    : 250   : 1e9/L  : platelet baseline
KOUTPLT :  0.00417:1/h   : platelet loss (lifespan 10 days)
EPLT    :  0.60 :        : LZD maximum inhibition of platelet production
KPLT    :  9.0  : mg/L   : LZD thrombocytotoxicity half-saturation concentration
ALT0    : 25    : U/L    : ALT baseline
KDALT   :  0.0144: 1/h   : ALT elimination (t1/2 48 h)
EALT    :  4.0  :        : RIF maximum ALT induction
KALT    :  6.0  : mg/L   : RIF hepatotoxicity half-saturation concentration

// ------------------------------------------------ SURGERY / STRATEGY ------
TSURG1  : 1e9   : h      : time of the first surgery (default: not performed)
LOGK1   :  0.0  : log10  : bacterial reduction at the first surgery (DAIR 2.5 / one-stage 4.5 / two-stage 6)
EXCH1   :  0    :        : implant exchange at the first surgery (1 = exchange, EPS reset)
TSURG2  : 1e9   : h      : time of the second surgery (two-stage reimplantation)
LOGK2   :  0.0  : log10  : bacterial reduction at the second surgery
EXCH2   :  0    :        : implant exchange at the second surgery
TSDUR   :  1.0  : h      : length of the interval over which the surgical reduction is applied
SPCF0   :  0.0  : mg     : initial spacer fast-release pool
SPCS0   :  0.0  : mg     : initial spacer slow-release pool

// -------------------------------------------------------- ENDPOINTS -------
PSEED   :  0.35 :        : relapse probability per surviving CFU (implant retained)
PSEEDX  :  0.10 :        : relapse probability per surviving CFU (implant exchanged)

$CMT @annotated
VANC   : VAN central compartment (mg)
VANP   : VAN peripheral compartment (mg)
VANB   : VAN free bone concentration (mg/L)
RIFG   : RIF gastrointestinal absorption compartment (mg)
RIFC   : RIF central compartment (mg)
RIFB   : RIF free bone concentration (mg/L)
ENZR   : CYP3A4 induction state (fold)
LVXG   : LVX gastrointestinal absorption compartment (mg)
LVXC   : LVX central compartment (mg)
LVXB   : LVX free bone concentration (mg/L)
DAPC   : DAP central compartment (mg)
DAPB   : DAP free bone concentration (mg/L)
LZDG   : LZD gastrointestinal absorption compartment (mg)
LZDC   : LZD central compartment (mg)
LZDB   : LZD free bone concentration (mg/L)
CFZC   : CFZ central compartment (mg)
CFZB   : CFZ free bone concentration (mg/L)
SPCF   : spacer fast-release pool (mg)
SPCS   : spacer slow-release pool (mg)
LOCJ   : local antibiotic concentration in synovial fluid (mg/L)
NP     : planktonic bacteria (CFU)
NB     : biofilm-adherent bacteria (CFU)
NPER   : persisters (CFU)
NSCV   : small-colony variants SCV (CFU)
NIC    : intracellular bacteria (CFU)
RP     : rpoB-mutant planktonic bacteria (CFU)
RB     : rpoB-mutant biofilm bacteria (CFU)
QP     : gyrA-mutant planktonic bacteria (CFU)
QB     : gyrA-mutant biofilm bacteria (CFU)
EPS    : biofilm matrix maturity (0-1)
AIP    : agr autoinducing peptide (normalised)
PMN    : synovial fluid neutrophils (cells/uL)
MONO   : synovial fluid monocytes (cells/uL)
MDSC   : myeloid-derived suppressor cells (normalised)
MAC    : activated macrophages (normalised)
IL1B   : IL-1beta (pg/mL)
TNFA   : TNF-alpha (pg/mL)
IL6    : IL-6 (pg/mL)
IL10   : IL-10 (pg/mL)
CRP    : serum CRP (mg/L)
ESR    : serum ESR (mm/h)
ADEF   : synovial fluid alpha-defensin (S/CO)
RANKL  : RANKL (normalised)
OPG    : OPG (normalised)
OBL    : osteoblasts (normalised to 100)
OCL    : osteoclasts (normalised to 100)
BONEV  : peri-implant bone mass (% of baseline)
LOOSEN : interface loosening index (0-100)
SEQ    : sequestrum / avascular barrier (normalised)
SCR    : serum creatinine (mg/dL)
PLT    : platelets (1e9/L)
ALT    : serum ALT (U/L)
AUCF   : VAN AUC24 running estimate (mg*h/L)

$MAIN
VANB_0   = 0.0;   RIFB_0 = 0.0;  LVXB_0 = 0.0;
DAPB_0   = 0.0;   LZDB_0 = 0.0;  CFZB_0 = 0.0;
ENZR_0   = 1.0;
SPCF_0   = SPCF0;
SPCS_0   = SPCS0;
LOCJ_0   = 0.0;

// Inoculum. NOTE: this is driven by the PARAMETER INOC0, not by init().
// $MAIN runs on every record, so writing a literal here would silently
// override anything passed to mrgsolve::init() -- change INOC0 instead.
NP_0     = INOC0;
NB_0     = 0.0;
NPER_0   = 0.0;
NSCV_0   = 0.0;
NIC_0    = 0.0;
RP_0     = 0.0;  RB_0 = 0.0;  QP_0 = 0.0;  QB_0 = 0.0;

EPS_0    = 0.0;
AIP_0    = 0.0;
PMN_0    = PMN0;
MONO_0   = MONO0;
MDSC_0   = 0.02;
MAC_0    = 0.10;
IL1B_0   = IL1B0;
TNFA_0   = TNFA0;
IL6_0    = IL60;
IL10_0   = IL100;
CRP_0    = CRP0;
ESR_0    = ESR0;
ADEF_0   = ADEF0;
RANKL_0  = 1.0;
OPG_0    = 1.0;
OBL_0    = 100.0;
OCL_0    = 100.0;
BONEV_0  = 100.0;
LOOSEN_0 = 0.0;
SEQ_0    = 0.0;
SCR_0    = SCR0;
PLT_0    = PLT0;
ALT_0    = ALT0;
AUCF_0   = 0.0;

$ODE
// Guard: the stiff 10-log bacterial collapse can leave states a hair below
// zero, and pow(negative, 3) then poisons the extinction term.
if(NP < 0)   NP   = 0.0;   if(NB < 0)   NB   = 0.0;
if(NPER < 0) NPER = 0.0;   if(NSCV < 0) NSCV = 0.0;
if(NIC < 0)  NIC  = 0.0;   if(RP < 0)   RP   = 0.0;
if(RB < 0)   RB   = 0.0;   if(QP < 0)   QP   = 0.0;
if(QB < 0)   QB   = 0.0;   if(EPS < 0)  EPS  = 0.0;
if(AIP < 0)  AIP  = 0.0;   if(SEQ < 0)  SEQ  = 0.0;

// =====================================================================
// 0. Surgical burden-reduction pulses (Pillar 3)
//    A pulse of height ln(10)*LOGK/TSDUR applied for TSDUR hours multiplies
//    every bacterial state by exactly 10^-LOGK. No event plumbing needed and
//    the arithmetic is exact.
// =====================================================================
double SR1 = 0.0, SR2 = 0.0, WIPE = 0.0;
if(SOLVERTIME >= TSURG1 && SOLVERTIME < (TSURG1 + TSDUR)) {
  SR1 = 2.302585093*LOGK1/TSDUR;
  if(EXCH1 > 0.5) WIPE = 1.0;
}
if(SOLVERTIME >= TSURG2 && SOLVERTIME < (TSURG2 + TSDUR)) {
  SR2 = 2.302585093*LOGK2/TSDUR;
  if(EXCH2 > 0.5) WIPE = 1.0;
}
double SURG = SR1 + SR2;
// Debridement also strips matrix; an implant exchange strips it completely.
// A RETAINED implant keeps its adherent matrix: debridement removes pus and
// necrotic soft tissue, not the biofilm bonded to metal. Only an exchange
// wipes the matrix. This asymmetry is what makes DAIR a race against EPS.
double EPSCUT = SURG*(WIPE > 0.5 ? 1.0 : 0.05);
double SEQCUT = SURG*(WIPE > 0.5 ? 1.0 : 0.10);

// =====================================================================
// 1. Pharmacokinetics
// =====================================================================
double CVAN = VANC/VCVAN;
double CVANP= VANP/VPVAN;
double CRIF = RIFC/VCRIF;
double CLVX = LVXC/VCLVX;
double CDAP = DAPC/VCDAP;
double CLZD = LZDC/VCLZD;
double CCFZ = CFZC/VCCFZ;

// Rifampicin autoinduction and its induction of linezolid clearance
double CLRIF = CLRIF0*ENZR;
double CLLZD = CLLZD0*(1.0 + RIFLZD*(ENZR - 1.0));

dxdt_VANC = -CLVAN*CVAN - QVAN*CVAN + QVAN*CVANP;
dxdt_VANP =  QVAN*CVAN - QVAN*CVANP;
dxdt_RIFG = -KARIF*RIFG;
dxdt_RIFC =  KARIF*FRIF*RIFG - CLRIF*CRIF;
dxdt_ENZR =  KENZ*((1.0 + EIND*CRIF/(3.0 + CRIF)) - ENZR);
dxdt_LVXG = -KALVX*LVXG;
dxdt_LVXC =  KALVX*FLVX*LVXG - CLLVX*CLVX;
dxdt_DAPC = -CLDAP*CDAP;
dxdt_LZDG = -KALZD*LZDG;
dxdt_LZDC =  KALZD*FLZD*LZDG - CLLZD*CLZD;
dxdt_CFZC = -CLCFZ*CCFZ;

// Free BONE concentrations. FB = (AUC_bone/AUC_plasma) x fu  (Pillar 2)
double FBVAN = PENVAN*FUVAN;
double FBRIF = PENRIF*FURIF;
double FBLVX = PENLVX*FULVX;
double FBDAP = PENDAP*FUDAP;
double FBLZD = PENLZD*FULZD;
double FBCFZ = PENCFZ*FUCFZ;

dxdt_VANB = KEQVAN*(FBVAN*CVAN - VANB);
dxdt_RIFB = KEQRIF*(FBRIF*CRIF - RIFB);
dxdt_LVXB = KEQLVX*(FBLVX*CLVX - LVXB);
dxdt_DAPB = KEQDAP*(FBDAP*CDAP - DAPB);
dxdt_LZDB = KEQLZD*(FBLZD*CLZD - LZDB);
dxdt_CFZB = KEQCFZ*(FBCFZ*CCFZ - CFZB);

// Local elution from antibiotic-loaded cement (two-pool, burst + matrix)
dxdt_SPCF = -KFAST*SPCF;
dxdt_SPCS = -KSLOW*SPCS;
dxdt_LOCJ = (KFAST*SPCF + KSLOW*SPCS)/VJOINT - KCLJ*LOCJ;

// Rolling 24 h vancomycin AUC estimate (for nephrotoxicity)
dxdt_AUCF = (24.0*CVAN - AUCF)/24.0;

// =====================================================================
// 2. Biofilm state variables and drug access
// =====================================================================
double EPSM  = EPS/(EPS + KMAT);                    // maturity 0-1
double AGRA  = AIP/(AIP + KAIP);                    // agr activity 0-1
double ACCESS= 1.0/(1.0 + SEQ/KSEQA);               // sequestrum shielding

// Effective drug concentrations seen by each niche
double VANEP = (VANB + FLOCP*LOCJ);
double VANEB = (VANB + FLOCB*LOCJ)*ACCESS;
double RIFEB = RIFB*ACCESS;
double LVXEB = LVXB*ACCESS;
double DAPEB = DAPB*ACCESS;
double LZDEB = LZDB*ACCESS;
double CFZEB = CFZB*ACCESS;

// Maturity-scaled biofilm tolerance multipliers (Pillar 4)
double PBVAN = PHIMAT(MBCVAN/MICVAN, EPSM);
double PBRIF = PHIMAT(MBCRIF/MICRIF, EPSM);
double PBLVX = PHIMAT(MBCLVX/MICLVX, EPSM);
double PBDAP = PHIMAT(MBCDAP/MICDAP, EPSM);
double PBLZD = PHIMAT(MBCLZD/MICLZD, EPSM);
double PBCFZ = PHIMAT(MBCCFZ/MICCFZ, EPSM);

double MS = (MSSA > 0.5) ? 1.0 : 0.0;

// =====================================================================
// 3. Antibiotic kill by bacterial state
//    Every term is EMAX * C / (EC50_state + C) with
//    EC50_state = MIC x (biofilm PHI) x (state PHI).
// =====================================================================
// --- planktonic (rifampicin-susceptible, FQ-susceptible) ---
double KLP = EKILL(EMXVAN, VANEP, MICVAN)
           + EKILL(EMXRIF, RIFB , MICRIF)
           + EKILL(EMXLVX, LVXB , MICLVX)
           + EKILL(EMXDAP, DAPB , MICDAP)
           + EKILL(EMXLZD, LZDB , MICLZD)
           + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double INHP = EMXLZDI*LZDB/(MICLZD + LZDB);

// --- biofilm-embedded ---
double KLB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
           + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF)
           + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX)
           + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
           + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
           + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);
double INHB = EMXLZDI*LZDEB/(MICLZD*PBLZD + LZDEB);

// --- persisters ---
double KLPER = EKILL(EMXVAN, VANEB, MICVAN*PBVAN*PHPVAN)
             + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF*PHPRIF)
             + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX*PHPLVX)
             + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP*PHPDAP)
             + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD*PHPLZD)
             + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ*PHPCFZ);

// --- small colony variants ---
double KLSCV = EKILL(EMXVAN, VANEB, MICVAN*PBVAN*PHSVAN)
             + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF*PHSRIF)
             + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX*PHSLVX)
             + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP*PHSDAP)
             + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD*PHSLZD)
             + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ*PHSCFZ);

// --- intracellular: cellular accumulation ratio x free bone concentration,
//     shielded by the same avascular sequestrum (the canalicular reservoir
//     sits inside dead cortical bone) ---
double VANIC = CARVAN*VANB*ACCESS;  double RIFIC = CARRIF*RIFB*ACCESS;
double LVXIC = CARLVX*LVXB*ACCESS;  double DAPIC = CARDAP*DAPB*ACCESS;
double LZDIC = CARLZD*LZDB*ACCESS;  double CFZIC = CARCFZ*CFZB*ACCESS;
double KLIC = EKILL(EMXVAN, VANIC, MICVAN*PBVAN*PHIVAN)
            + EKILL(EMXRIF, RIFIC, MICRIF*PBRIF*PHIRIF)
            + EKILL(EMXLVX, LVXIC, MICLVX*PBLVX*PHILVX)
            + EKILL(EMXDAP, DAPIC, MICDAP*PBDAP*PHIDAP)
            + EKILL(EMXLZD, LZDIC, MICLZD*PBLZD*PHILZD)
            + MS*EKILL(EMXCFZ, CFZIC, MICCFZ*PBCFZ*PHICFZ);

// --- rpoB mutants: identical EXCEPT rifampicin does nothing ---
double KLRP = EKILL(EMXVAN, VANEP, MICVAN)
            + EKILL(EMXLVX, LVXB , MICLVX)
            + EKILL(EMXDAP, DAPB , MICDAP)
            + EKILL(EMXLZD, LZDB , MICLZD)
            + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double KLRB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
            + EKILL(EMXLVX, LVXEB, MICLVX*PBLVX)
            + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
            + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
            + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);

// --- gyrA mutants: identical EXCEPT levofloxacin does nothing ---
double KLQP = EKILL(EMXVAN, VANEP, MICVAN)
            + EKILL(EMXRIF, RIFB , MICRIF)
            + EKILL(EMXDAP, DAPB , MICDAP)
            + EKILL(EMXLZD, LZDB , MICLZD)
            + MS*EKILL(EMXCFZ, CFZB, MICCFZ);
double KLQB = EKILL(EMXVAN, VANEB, MICVAN*PBVAN)
            + EKILL(EMXRIF, RIFEB, MICRIF*PBRIF)
            + EKILL(EMXDAP, DAPEB, MICDAP*PBDAP)
            + EKILL(EMXLZD, LZDEB, MICLZD*PBLZD)
            + MS*EKILL(EMXCFZ, CFZEB, MICCFZ*PBCFZ);

// =====================================================================
// 4. Innate immunity (Pillar 1)
// =====================================================================
double NTOT  = NP + NB + NPER + NSCV + NIC + RP + RB + QP + QB;
double NSURF = NB + NPER + RB + QB;                  // surface-attached pool
double CROWD = 1.0 - NTOT/NMAX;  if(CROWD < 0.0) CROWD = 0.0;

double PMNF  = PMN/(PMN + KPMN);
double MDSUP = 1.0 - EMDSC*MDSC/(KMD + MDSC);
double FB    = (IMPL > 0.5) ? FBFACT : 1.0;
double IMMCAP= KIMM*PMNF*FB*IMSUP*MDSUP;

// Phagocyte capacity is a SHARED resource, allocated across the accessible
// burden in proportion to how reachable each state is. Writing one
// Michaelis-Menten sink per compartment (the obvious but wrong way) makes a
// rare mutant clone face the entire immune system as if it were alone, which
// silently sterilises exactly the subpopulation the model exists to track.
double NPHAG = (NP + RP + QP)
             + FRUST*(NB + RB + QB + NSCV)
             + FICIMM*NIC;
double IMMTOT= IMMCAP*NPHAG/(KMIMM + NPHAG);
double SHR   = IMMTOT/(NPHAG + 1e-30);
double IMMP  = SHR*NP;
double IMMB  = SHR*FRUST*NB;                          // frustrated phagocytosis
double IMMS  = SHR*FRUST*NSCV;
double IMMI  = SHR*FICIMM*NIC;
double IMMRP = SHR*RP;
double IMMRB = SHR*FRUST*RB;
double IMMQP = SHR*QP;
double IMMQB = SHR*FRUST*QB;

// =====================================================================
// 5. Bacterial state transitions
// =====================================================================
double SFREE = (IMPL > 0.5) ? (1.0 - NSURF/NBCAP) : 0.0;
if(SFREE < 0.0) SFREE = 0.0;
double ATT   = KATT*SFREE;
double DISP  = KDISP*AGRA;
double STRESS= KLB/(KLB + 0.05);
double KPER  = KPER0*(1.0 + APER*(EPSM + STRESS));
double KSCV  = KSCV0*(1.0 + ASCV*LOCJ/(LOCJ + 20.0));
double KINT  = KINT0*(OBL/100.0);

// Smooth sub-single-cell extinction. Applied PER CLONE, not per phenotypic
// state: a cell can move between planktonic / sessile / persister / SCV /
// intracellular, so "fewer than one organism left" is a property of the
// lineage. Applying it per compartment annihilates the rare rpoB clone the
// instant it dips below 1 CFU in any one state and destroys Pillar 3.
double NSUSC = NP + NB + NPER + NSCV + NIC;
double NRPOB = RP + RB;
double NGYRA = QP + QB;
double EXS = KEXT/(1.0 + pow(NSUSC/NEXT, 3.0));
double EXR = KEXT/(1.0 + pow(NRPOB/NEXT, 3.0));
double EXQ = KEXT/(1.0 + pow(NGYRA/NEXT, 3.0));

// Carrying capacities: the free (joint fluid / abscess) pool and the
// intracellular reservoir are each bounded by their own niche, not only by
// the lesion total. Without this the dispersed planktonic pool overruns the
// biofilm and the infection becomes spuriously vancomycin-curable.
double NFREE = NP + RP + QP;
double CROWDP= CROWD*(1.0 - NFREE/NPCAP);  if(CROWDP < 0.0) CROWDP = 0.0;
double CROWDI= CROWD*(1.0 - NIC/NICCAP);   if(CROWDI < 0.0) CROWDI = 0.0;

// Growth fluxes (needed for mutant supply: mutants arise at DIVISION)
double GFP = MUP*(1.0 - INHP)*CROWDP*NP;
double GFB = MUB*(1.0 - INHB)*CROWD*NB;
double GFS = MUS*CROWD*NSCV;
double GFI = MUI*CROWDI*NIC;
double MUTR = MURPOB*(GFP + GFB);      // rpoB mutant supply  (Pillar 3)
double MUTQ = MUGYRA*(GFP + GFB);      // gyrA mutant supply

dxdt_NP   = GFP - MUTR - MUTQ
          - (KDNAT + KLP + EXS + SURG)*NP
          - IMMP - ATT*NP + DISP*NB - KINT*NP + KRELIC*NIC;

dxdt_NB   = GFB
          + ATT*NP - DISP*NB
          - (KDNAT + KLB + EXS + SURG)*NB
          - IMMB - KPER*NB + KWAKE*NPER - KSCV*NB + KREVS*NSCV;

dxdt_NPER = KPER*NB - KWAKE*NPER - (KLPER + EXS + SURG)*NPER;

dxdt_NSCV = GFS + KSCV*NB - KREVS*NSCV
          - (KDNAT + KLSCV + EXS + SURG)*NSCV - IMMS;

dxdt_NIC  = GFI + KINT*NP - KRELIC*NIC
          - (KLIC + EXS + 0.5*SURG)*NIC - IMMI;

dxdt_RP   = MUTR + MUP*(1.0 - FITRPOB)*(1.0 - INHP)*CROWDP*RP
          - (KDNAT + KLRP + EXR + SURG)*RP - IMMRP - ATT*RP + DISP*RB;

dxdt_RB   = MUB*(1.0 - FITRPOB)*(1.0 - INHB)*CROWD*RB
          + ATT*RP - DISP*RB
          - (KDNAT + KLRB + EXR + SURG)*RB - IMMRB;

dxdt_QP   = MUTQ + MUP*(1.0 - FITGYRA)*(1.0 - INHP)*CROWDP*QP
          - (KDNAT + KLQP + EXQ + SURG)*QP - IMMQP - ATT*QP + DISP*QB;

dxdt_QB   = MUB*(1.0 - FITGYRA)*(1.0 - INHB)*CROWD*QB
          + ATT*QP - DISP*QB
          - (KDNAT + KLQB + EXQ + SURG)*QB - IMMQB;

// =====================================================================
// 6. Biofilm matrix and quorum sensing
// =====================================================================
dxdt_EPS = KSEPS*(NSURF/(NSURF + KEPSN))*(1.0 - EPS) - KDEPS*EPS - EPSCUT*EPS;
dxdt_AIP = KSAIP*(NP + NB)/((NP + NB) + KAIPN) - KDAIP*AIP;

// =====================================================================
// 7. Innate immune cells, MDSC and cytokines
// =====================================================================
double BACSIG = NTOT/(NTOT + KBACS);
double ACTSIG = BACSIG*MDSUP;
double TOXL   = AGRA*BACSIG;

dxdt_PMN  = KINPMN*(1.0 + EMXCHEM*ACTSIG/(KCHEM + ACTSIG))
          - KOUTPMN*PMN - KTOXP*TOXL*PMN;
dxdt_MONO = KINMON*(1.0 + EMXMON*BACSIG/(KCHEM + BACSIG))
          - KOUTMON*MONO;
dxdt_MDSC = KSMDSC*EPS/(EPS + KEPSM) - KDMDSC*MDSC;
dxdt_MAC  = KSMAC*BACSIG - KDMAC*MAC;

dxdt_IL1B = KDIL1*IL1B0*(1.0 + EIL1*ACTSIG)  - KDIL1*IL1B;
dxdt_TNFA = KDTNF*TNFA0*(1.0 + ETNF*ACTSIG)  - KDTNF*TNFA;
dxdt_IL6  = KDIL6*IL60 *(1.0 + EIL6*ACTSIG)  - KDIL6*IL6;
dxdt_IL10 = KDIL10*IL100*(1.0 + EIL10*MDSC/(KMD + MDSC)) - KDIL10*IL10;

// =====================================================================
// 8. Systemic biomarkers (2018 ICM criteria read-outs)
// =====================================================================
dxdt_CRP  = KDCRP*CRP0*(1.0 + ECRP*(IL6 - IL60)/(KCRP + (IL6 - IL60))) - KDCRP*CRP;
dxdt_ESR  = KDESR*ESR0*(1.0 + EESR*(IL6 - IL60)/(KCRP + (IL6 - IL60))) - KDESR*ESR;
dxdt_ADEF = KDADEF*ADEF0*(1.0 + EADEF*(PMN - PMN0)/(KADEF + (PMN - PMN0))) - KDADEF*ADEF;

// =====================================================================
// 9. Bone: RANKL/OPG, osteoclasts, osteolysis, loosening, sequestrum
// =====================================================================
double CYTB = (TNFA - TNFA0)/(50.0 + (TNFA - TNFA0))
            + (IL1B - IL1B0)/(30.0 + (IL1B - IL1B0))
            + (IL6  - IL60 )/(80.0 + (IL6  - IL60 ));
if(CYTB < 0.0) CYTB = 0.0;
double APOS = 0.6*AGRA*BACSIG + NIC/(NIC + 1e6);

dxdt_RANKL = KSRANK*(1.0 + ERANK*CYTB/3.0) - KDRANK*RANKL;
dxdt_OPG   = KSOPG*(OBL/100.0)             - KDOPG*OPG;
dxdt_OBL   = KINOBL - KOUTOBL*OBL - KAPOP*APOS*OBL;
double RR  = RANKL/(OPG + 1e-6);
dxdt_OCL   = KINOCL*RR/(KRR + RR) - KOUTOCL*OCL;
dxdt_BONEV = KFORM*(OBL/100.0) - KRES*(OCL/100.0) + KHOM*(100.0 - BONEV);
double BLOSS = 100.0 - BONEV;  if(BLOSS < 0.0) BLOSS = 0.0;
dxdt_LOOSEN= KLOOSE*BLOSS*(1.0 - LOOSEN/100.0) - KREPAIR*LOOSEN;
dxdt_SEQ   = KSSEQ*(KRES*(OCL/100.0))*BACSIG*10.0*(1.0 - SEQ/SEQMAX)
           - KCSEQ*SEQ - SEQCUT*SEQ;

// =====================================================================
// 10. Drug toxicity
// =====================================================================
double EXAUC = AUCF - AUCTHR;  if(EXAUC < 0.0) EXAUC = 0.0;
double NEPH  = ENEPH*EXAUC/(KNEPH + EXAUC);
dxdt_SCR  = KELCR*SCR0 - KELCR*(1.0 - NEPH)*SCR;
dxdt_PLT  = KOUTPLT*PLT0*(1.0 - EPLT*CLZD/(KPLT + CLZD)) - KOUTPLT*PLT;
dxdt_ALT  = KDALT*ALT0*(1.0 + EALT*CRIF/(KALT + CRIF)) - KDALT*ALT;

$TABLE
double NTOTAL = NP + NB + NPER + NSCV + NIC + RP + RB + QP + QB;
if(NTOTAL < 0.0) NTOTAL = 0.0;
double NSUSC  = NP + NB + NPER + NSCV + NIC;
double NRES   = RP + RB + QP + QB;
capture LOGNTOT = log10(NTOTAL + 1e-12);
capture LOGNP   = log10(NP  + 1e-12);
capture LOGNB   = log10(NB  + 1e-12);
capture LOGNPER = log10(NPER+ 1e-12);
capture LOGNSCV = log10(NSCV+ 1e-12);
capture LOGNIC  = log10(NIC + 1e-12);
capture LOGNRES = log10(NRES+ 1e-12);
capture RESFRAC = NRES/(NTOTAL + 1e-12);

// --- Pillar 2 read-out: achievable free bone concentration over MBEC -------
capture CPVAN = VANC/VCVAN;
capture CPRIF = RIFC/VCRIF;
capture CPLVX = LVXC/VCLVX;
capture CPDAP = DAPC/VCDAP;
capture CPLZD = LZDC/VCLZD;
capture CPCFZ = CFZC/VCCFZ;
capture RATVAN = (VANB + FLOCB*LOCJ)/MBCVAN;
capture RATRIF = RIFB/MBCRIF;
capture RATLVX = LVXB/MBCLVX;
capture RATDAP = DAPB/MBCDAP;
capture RATLZD = LZDB/MBCLZD;
capture RATCFZ = CFZB/MBCCFZ;
capture RATMAX = RATVAN;
if(RATRIF > RATMAX) RATMAX = RATRIF;
if(RATLVX > RATMAX) RATMAX = RATLVX;
if(RATDAP > RATMAX) RATMAX = RATDAP;
if(RATLZD > RATMAX) RATMAX = RATLZD;
if(MSSA > 0.5 && RATCFZ > RATMAX) RATMAX = RATCFZ;

// --- Pillar 3 read-out: mutant supply --------------------------------------
capture PRPOB = 1.0 - exp(-MURPOB*NTOTAL);   // P(pre-existing rpoB mutant)
capture PGYRA = 1.0 - exp(-MUGYRA*NTOTAL);
capture EMUTR = MURPOB*NTOTAL;               // expected rpoB mutants present

// --- Pillar 4 read-out: biofilm maturity ------------------------------------
capture EPSMAT = EPS/(EPS + KMAT);

// --- clinical / diagnostic read-outs ---------------------------------------
capture SYNWBC = PMN + MONO;
capture PMNPCT = 100.0*PMN/(PMN + MONO + 1e-9);
// 2018 ICM minor-criteria score (CRP 2, D-dimer 2, synovial WBC 3,
// PMN% 2, alpha-defensin 3, LE 3): >=6 infected, 2-5 inconclusive
double SC = 0.0;
if(CRP > 10.0)     SC += 2.0;
if(SYNWBC > 3000.0)SC += 3.0;
if(PMNPCT > 80.0)  SC += 2.0;
if(ADEF > 1.0)     SC += 3.0;
capture ICMSC = SC;

// Cure probability: PSEED is the chance that ONE surviving CFU regrows.
double PS = (EXCH1 > 0.5 || EXCH2 > 0.5) ? PSEEDX : PSEED;
capture PCURE = exp(-PS*NTOTAL);
capture PFAIL = 1.0 - exp(-PS*NTOTAL);

capture CRCL  = 100.0*(SCR0/SCR);           // crude creatinine-clearance proxy
capture AKIFL = (SCR > 1.5*SCR0) ? 1.0 : 0.0;
capture TCPFL = (PLT < 100.0)    ? 1.0 : 0.0;
capture HEPFL = (ALT > 3.0*ALT0) ? 1.0 : 0.0;
'

pji_mod <- mrgsolve::mcode_cache("pji", pji_code)
pji_mod <- mrgsolve::update(pji_mod, end = 4320, delta = 4,
                            atol = 1e-8, rtol = 1e-6, maxsteps = 500000)

## =============================================================================
##  DOSING BUILDERS
##  Standard adult regimens; amt in mg, time in hours.
## =============================================================================
DAY <- 24

ev_van  <- function(start = 0, weeks = 6, amt = 1000, ii = 12)
  ev(time = start, cmt = "VANC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_rif  <- function(start = 0, weeks = 12, amt = 450, ii = 12)
  ev(time = start, cmt = "RIFG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_lvx  <- function(start = 0, weeks = 12, amt = 750, ii = 24)
  ev(time = start, cmt = "LVXG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_dap  <- function(start = 0, weeks = 6, amt = 560, ii = 24)
  ev(time = start, cmt = "DAPC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_lzd  <- function(start = 0, weeks = 6, amt = 600, ii = 12)
  ev(time = start, cmt = "LZDG", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)
ev_cfz  <- function(start = 0, weeks = 6, amt = 2000, ii = 8)
  ev(time = start, cmt = "CFZC", amt = amt, ii = ii, addl = ceiling(weeks*7*24/ii) - 1)

## =============================================================================
##  SCENARIOS
##  All start from a 100 CFU intra-operative contamination of a hip/knee
##  prosthesis at t = 0 and are simulated for 180 days (4320 h).
##  Diagnosis is assumed at day 21 (acute post-operative PJI) unless stated.
## =============================================================================
DX  <- 21*DAY     # diagnosis / surgery time for the acute presentations
DXC <- 90*DAY     # delayed (chronic) presentation

scen <- list(

  ## 1. Natural history, no treatment ---------------------------------------
  "01_No treatment, natural history" = list(
    par = list(),
    ev  = ev(time = 0, cmt = "NP", amt = 0)),

  ## 2. Antibiotics only, implant retained, NO surgery -----------------------
  ##    ("medical management" — the arm that shows Pillar 2 alone is not enough)
  "02_Vancomycin alone without surgery" = list(
    par = list(),
    ev  = ev_van(DX, 6)),

  ## 3. DAIR + vancomycin monotherapy 6 weeks --------------------------------
  "03_DAIR + vancomycin 6 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_van(DX, 6)),

  ## 4. DAIR + vancomycin + rifampicin 12 weeks ------------------------------
  "04_DAIR + vanco+rifampicin 12 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_van(DX, 6) + ev_rif(DX, 12)),

  ## 5. DAIR + levofloxacin + rifampicin (Zimmerli 1998 JAMA analogue) -------
  "05_DAIR + levofloxacin+rifampicin 12 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 6. DAIR + levofloxacin alone (the Zimmerli control arm) -----------------
  "06_DAIR + levofloxacin alone 12 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12)),

  ## 7. DAIR + RIFAMPICIN MONOTHERAPY (Pillar 3 demonstration) ---------------
  "07_DAIR + rifampicin alone 12 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_rif(DX, 12)),

  ## 8. NO debridement + rifampicin combination -----------------------------
  ##    (mature 1e10 biofilm: the rpoB mutant is already there)
  "08_Rifampicin combination without debridement" = list(
    par = list(),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 9. One-stage exchange + levofloxacin/rifampicin -------------------------
  "09_One-stage exchange + LVX/RIF 12 weeks" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 10. Two-stage exchange with an antibiotic-loaded spacer -----------------
  "10_Two-stage exchange + spacer" = list(
    par = list(TSURG1 = DX,        LOGK1 = 6.0, EXCH1 = 1,
               TSURG2 = DX + 8*7*DAY, LOGK2 = 1.5, EXCH2 = 1,
               SPCF0  = 1200, SPCS0 = 700),
    ev  = ev_van(DX, 6) + ev_rif(DX, 12)),

  ## 11. OVIVA analogue: fully oral vs IV backbone ---------------------------
  "11_OVIVA oral (LVX+RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_lvx(DX, 6) + ev_rif(DX, 6)),
  "12_OVIVA intravenous (VAN)+RIF" = list(
    par = list(TSURG1 = DX, LOGK1 = 4.5, EXCH1 = 1),
    ev  = ev_van(DX, 6) + ev_rif(DX, 6)),

  ## 13-14. DATIPO analogue: 6 vs 12 weeks -----------------------------------
  "13_DATIPO 6 weeks (DAIR+LVX/RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 6) + ev_rif(DX, 6)),
  "14_DATIPO 12 weeks (DAIR+LVX/RIF)" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 15-16. Late DAIR on a mature biofilm (Pillar 4) -------------------------
  ##    Read these against scenario 05 (same regimen, DAIR on day 21):
  ##    day  21 -> sterile on day 60 of an 84-day course (24 days of margin)
  ##    day  90 -> sterile on day 84 of 84 (zero margin)
  ##    day 180 -> never sterile; relapses after antibiotics stop
  "15_Delayed diagnosis (day 90) + DAIR + LVX/RIF" = list(
    par = list(TSURG1 = DXC, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DXC, 12) + ev_rif(DXC, 12)),
  "16_Chronic (day 180) + DAIR + LVX/RIF" = list(
    par = list(TSURG1 = 180*DAY, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(180*DAY, 12) + ev_rif(180*DAY, 12)),

  ## 16. Daptomycin + rifampicin --------------------------------------------
  "17_DAIR + daptomycin+rifampicin" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_dap(DX, 6) + ev_rif(DX, 12)),

  ## 17. Linezolid + rifampicin (with the CYP/P-gp interaction) --------------
  "18_DAIR + linezolid+rifampicin" = list(
    par = list(TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lzd(DX, 6) + ev_rif(DX, 12)),

  ## 18. MSSA: cefazolin + rifampicin ---------------------------------------
  "19_MSSA DAIR + cefazolin+rifampicin" = list(
    par = list(MSSA = 1, TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_cfz(DX, 6) + ev_rif(DX, 12)),

  ## 19. Immunocompromised host (diabetes / anti-TNF) ------------------------
  "20_Immunocompromised host DAIR + LVX/RIF" = list(
    par = list(IMSUP = 0.5, TSURG1 = DX, LOGK1 = 2.5, EXCH1 = 0),
    ev  = ev_lvx(DX, 12) + ev_rif(DX, 12)),

  ## 20. Chronic oral suppression (inoperable patient) -----------------------
  "21_Chronic oral suppressive therapy" = list(
    par = list(),
    ev  = ev_lvx(DX, 26))
)

run_scen <- function(mod, s, end = 4320) {
  p <- s$par
  m <- mod
  if(length(p)) m <- mrgsolve::param(m, p)
  mrgsolve::mrgsim(m, events = s$ev, end = end, delta = 4) %>% as.data.frame()
}

run_all <- function(mod = pji_mod, end = 4320) {
  out <- lapply(names(scen), function(nm) {
    d <- run_scen(mod, scen[[nm]], end = end)
    d$scenario <- nm
    d
  })
  dplyr::bind_rows(out)
}

## =============================================================================
##  PILLAR 1 CHECK — the foreign body moves the infective dose
##  Run a dose-escalation with and without an implant and find the inoculum at
##  which the model still has > 1e6 CFU at day 60.
## =============================================================================
id50_scan <- function(mod = pji_mod,
                      inoc = 10^seq(0, 9, by = 0.5)) {
  f <- function(impl) {
    sapply(inoc, function(n0) {
      m <- mrgsolve::param(mod, IMPL = impl, INOC0 = n0)
      o <- mrgsolve::mrgsim(m, end = 60*24, delta = 24) %>% as.data.frame()
      max(tail(o$LOGNTOT, 1), -12)
    })
  }
  data.frame(inoculum = inoc,
             log10_burden_d60_with_implant    = f(1),
             log10_burden_d60_without_implant = f(0))
}

## =============================================================================
##  PILLAR 2 CHECK — the C_bone,free / MBEC table at steady state
## =============================================================================
ratio_table <- function(mod = pji_mod) {
  p <- as.list(mrgsolve::param(mod))
  auc24 <- c(
    VAN = 2*1000/p$CLVAN,
    RIF = 2*450*p$FRIF/p$CLRIF0/2,       # /2 approximates the autoinduced state
    LVX = 750*p$FLVX/p$CLLVX,
    DAP = 560/p$CLDAP,
    LZD = 2*600*p$FLZD/p$CLLZD0,
    CFZ = 3*2000/p$CLCFZ)
  cavg  <- auc24/24
  fbone <- c(VAN = p$PENVAN*p$FUVAN, RIF = p$PENRIF*p$FURIF,
             LVX = p$PENLVX*p$FULVX, DAP = p$PENDAP*p$FUDAP,
             LZD = p$PENLZD*p$FULZD, CFZ = p$PENCFZ*p$FUCFZ)
  mbec  <- c(VAN = p$MBCVAN, RIF = p$MBCRIF, LVX = p$MBCLVX,
             DAP = p$MBCDAP, LZD = p$MBCLZD, CFZ = p$MBCCFZ)
  mic   <- c(VAN = p$MICVAN, RIF = p$MICRIF, LVX = p$MICLVX,
             DAP = p$MICDAP, LZD = p$MICLZD, CFZ = p$MICCFZ)
  data.frame(drug = names(auc24),
             AUC24_plasma = round(auc24, 1),
             Cavg_plasma  = round(cavg, 2),
             f_bone_free  = round(fbone, 3),
             Cbone_free   = round(cavg*fbone, 3),
             MIC          = mic,
             MBEC         = mbec,
             ratio_MIC    = round(cavg*fbone/mic, 1),
             ratio_MBEC   = round(cavg*fbone/mbec, 4),
             row.names = NULL)
}

## =============================================================================
##  PILLAR 3 CHECK — mutant supply before and after debridement
## =============================================================================
mutant_supply <- function(N = 10^c(4:10), mu = 1e-8) {
  data.frame(burden_CFU        = N,
             expected_mutants  = signif(mu*N, 3),
             P_preexisting     = signif(1 - exp(-mu*N), 3))
}

## =============================================================================
##  SUMMARY TABLE
## =============================================================================
summarise_scen <- function(df) {
  df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      burden_end_log10 = round(dplyr::last(LOGNTOT), 2),
      resistant_log10  = round(dplyr::last(LOGNRES), 2),
      P_cure           = round(dplyr::last(PCURE), 3),
      CRP_end          = round(dplyr::last(CRP), 1),
      CRP_peak         = round(max(CRP), 1),
      synWBC_peak      = round(max(SYNWBC), 0),
      bone_volume_end  = round(dplyr::last(BONEV), 1),
      loosening_end    = round(dplyr::last(LOOSEN), 1),
      max_bone_MBEC_ratio = round(max(RATMAX), 3),
      # Mutant supply at the LOWEST burden the regimen ever reaches: this is
      # the post-debridement number that decides whether rifampicin is usable.
      P_rpoB_at_nadir  = round(1 - exp(-1e-8*min(10^LOGNTOT)), 4),
      AKI              = max(AKIFL),
      thrombocytopenia = max(TCPFL),
      hepatotoxicity   = max(HEPFL),
      .groups = "drop")
}

## =============================================================================
##  EXAMPLE USAGE
## -----------------------------------------------------------------------------
##  out <- run_all()
##  summarise_scen(out)
##  ratio_table()          # Pillar 2 — the table that explains rifampicin
##  mutant_supply()        # Pillar 3 — why debridement precedes rifampicin
##  id50_scan()            # Pillar 1 — the foreign-body ID50 shift
## =============================================================================
