## ============================================================================
##  Chronic osteomyelitis · implant-associated bone and joint infection
##  Chronic Osteomyelitis & Implant-Associated Osteoarticular Infection
##  QSP model for mrgsolve  —  47 ODE states · 2 systemic drugs + 1 local depot
##                             · 5 bacterial subpopulations · 11 scenarios
##
##  ---------------------------------------------------------------------------
##  ORGANISING THESIS (why this model exists)
##  ---------------------------------------------------------------------------
##  In chronic osteomyelitis the drug and the bacteria do not live in the same
##  compartment. And treatment failure arises because two independent penalties are **multiplied**.
##
##  (1) The spatial penalty — sequestrum and implant biofilm have no perfusion.
##      The drug enters by diffusion alone, and on the way in it is lost to binding, degradation and bacterial uptake.
##      This is a reaction-diffusion problem, and the steady-state penetration efficiency closes with the Thiele modulus:
##
##          phi  = L * sqrt(k_loss / D_eff)
##          eta  = tanh(phi) / phi
##          C_seq = eta * C_bone
##
##      phi=1 -> eta 0.76 · phi=3 -> 0.33 · phi=6 -> 0.17 · phi=10 -> 0.10.
##      This is a 'concentration' problem, so dose and local delivery can push it. (EC50 shift)
##
##  (2) The phenotypic penalty — biofilm bacteria are not 'slow bacteria'. The kill rate of
##      beta-lactams, glycopeptides and aminoglycosides is proportional to the growth rate, so
##      at mu -> 0 the kill -> 0, and this is independent of concentration:
##
##          KILL_j = C_j/(C_j+EC50_j) * [ Emax_gd,j*(mu/mu_max) + Emax_gi,j ]
##
##      To a drug with Emax_gi = 0 the persister is not 'resistant' but 'invisible'.
##      No dose whatsoever can push this. (Emax collapse)
##
##  So the observed MBEC/MIC = 100~1000 is not a potency shift; it is mostly
##  the shadow of Emax collapse. This model separates the two and parameterises each.
##  (SCENARIO 7 below demonstrates that separation as an experiment: local vancomycin at
##   1000 mg/L — 1000-fold the MIC — still fails to sterilise the biofilm.)
##
##  Only the rifamycins retain their Emax in non-growing bacteria (Emax_gi is large). That is why
##  combining rifampicin is the rule in staphylococcal implant infection, and at the same time why
##  using it alone is absolutely contraindicated:  P(rpoB mutant present) = 1 - exp(-f_mut * B_total),
##  f_mut=1e-8 · B_total=1e9 -> P ~ 1.0. Rule and contraindication come from the same equation.
##
##  So the three interventions mathematically touch **different objects**:
##
##      surgical debridement/removal = change of initial condition (basin change) : lowers P0
##      rifampicin combination      = change of kill ceiling (RSTER change)       : lifts it above 1
##      duration of therapy         = change of time to reach                     : and nothing beyond that
##
##      RSTER = sum_j KILL_j(biofilm) / (mu_bf + persister supply)
##      RSTER <= 1  ->  a stable non-zero bacterial equilibrium. Extend the duration
##                      for ever and sterilisation never comes. (= why it becomes 'suppressive therapy')
##      RSTER >  1  ->  sterilisation is possible. All that remains is time, and that time
##                      depends on P0 only logarithmically, T_ster ~ ln(P0)/k_eff.
##
##  Finally, the pathophysiology **builds its own pharmacokinetic barrier**. This is the most
##  important closed loop in the model:
##
##      infection -> RANKL/OPG -> bone resorption + intramedullary pressure -> Haversian vessel thrombosis
##           -> PERF falls -> (a) bone interstitial drug exposure falls  (b) neutrophil access falls
##           -> bone necrosis -> SEQ up -> L up -> phi up -> eta down
##           -> C_seq down -> killing down -> infection persists  ->  (back to the start)
##
##  When SEQ grows 2 -> 15 cm3, L becomes 0.18 -> 0.32 cm and the eta of vancomycin
##  halves from 0.28 -> 0.16. Which means the same drug works less well in a patient who
##  presents late not because 'the organism got stronger' but because 'the geometry changed'.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION NOTES (basis of the parameters · key clinical trials)
##  ---------------------------------------------------------------------------
##  * bone:serum AUC ratio (ARBS/BRBS) — vancomycin 0.10-0.20, beta-lactams 0.15-0.30,
##    fluoroquinolones 0.30-0.50, linezolid 0.40-0.50, rifampicin 0.30-0.60,
##    clindamycin 0.40-0.80 (Landersdorfer 2009 Clin Pharmacokinet review;
##    Thabit 2019 Expert Opin Pharmacother). ARBS is defined as total-bone:total-plasma,
##    and the PD calculation uses the free concentration, i.e. multiplied by the tissue unbound fraction AFU.
##  * intracellular:extracellular ratio (ARIC/BRIC) — rifampicin 5-10, FQ 4-8, clindamycin 5-10,
##    vancomycin/beta-lactams <0.1 (Carryn 2003 Infect Dis Clin North Am;
##    Valour 2015 Antimicrob Agents Chemother — S. aureus inside osteoblasts).
##  * MBEC/MIC 100-1000 (Ceri 1999 J Clin Microbiol Calgary device). This model
##    reproduces that observation by raising EC50 only 3-fold (FBFEC 3.0, the EPS-binding share) and
##    leaving the rest to mu/mu_max gating — this is the testable claim of the model.
##  * persister awakening rate KWAKE 0.08/d — matched to the observation that reactivation of
##    staphylococcal biofilm persisters runs on a scale of days to weeks, and to the clinical time
##    axis of 'relapse 6-24 months after stopping' (Conlon 2013 Nature; Lewis 2010 Annu Rev Microbiol).
##  * rpoB mutation frequency 1e-8/division (Aubry-Damon 1998 AAC; O'Neill 2001 JAC).
##  * Zimmerli 1998 JAMA (RCT) — DAIR + ciprofloxacin ± rifampicin for 3-6 months:
##    treatment success 12/12 (100%) with rifampicin vs 7/12 (58%) without. In this model
##    the contrast of SCENARIO 3 (no rifampicin) vs 4 (with rifampicin) reproduces that result.
##  * OVIVA 2019 NEJM (Li HK, n=1054) — oral non-inferior to IV over the first 6 weeks (treatment
##    failure 13.2% vs 14.6%). In this model the route changes only F and CL, so it follows
##    structurally that a drug of adequate bioavailability gives the same outcome (SCENARIO 5).
##  * DATIPO 2021 NEJM (Bernard L) — 6 vs 12 weeks in prosthetic joint infection: 6 weeks was
##    NOT non-inferior (failure 18.1% vs 9.4%). Evidence that duration still matters.
##    The SCENARIO 9 duration sweep in this model reproduces that asymmetry (with adequate
##    debridement 6 weeks suffices; with the implant retained even 12 weeks is not enough).
##  * vancomycin AUC24 target 400-600 mg*h/L, AKI risk rises above >600
##    (Rybak 2020 consensus guideline).
##  * linezolid beyond >28 days: thrombocytopenia and mitochondrial toxicity climb steeply.
##  * rifampicin induction of CYP3A4/2C9/P-gp — linezolid AUC -30%, moxifloxacin -30%
##    (Egle 2005 Clin Pharmacol Ther; Weiner 2007 CID). Induction t1/2 ~ 2-3 d.
##  * in diabetic foot bone infection, the perfusion expressed by ABI/TcPO2 is the strongest predictor of outcome
##    (Lipsky 2012/2023 IWGDF guidance).
##
##  ---------------------------------------------------------------------------
##  Disclaimer
##  ---------------------------------------------------------------------------
##  A semi-quantitative QSP model for teaching and research. Do not use it for actual
##  prescribing, clinical decisions or regulatory submission. The parameters are illustrative approximations picked from the range of the published literature.
## ============================================================================

library(mrgsolve)
suppressMessages({
  has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
})

## ============================================================================
##  MODEL CODE
##  (Comments inside the model code are written in ASCII only, so that the generated
##   C++ compiles safely. The Korean prose is kept on the R-comment side of this file.)
## ============================================================================

com_code <- '
$PROB
Chronic osteomyelitis / implant-associated bone infection QSP model.
47 ODEs. Two systemic drugs (A = growth-dependent backbone, B = rifamycin-like
growth-INdependent partner) plus one local (cement/bead) depot of drug A.
Five bacterial subpopulations with distinct drug accessibility.
Central claim: failure = (spatial penalty eta) x (phenotypic penalty mu/mu_max).

$GLOBAL
#include <cmath>
#define SAFE 1.0e-12

// clamp helpers ------------------------------------------------------------
double posv(double x) { return (x > 0.0) ? x : 0.0; }
// A continuum ODE happily regrows a population from 1e-25 CFU. Real bacteria
// come in integers, so growth is switched off smoothly below ~1 CFU. Without
// this the deterministic solution always relapses, even after sterilisation.
double extf(double n) { return posv(n) / (posv(n) + 1.0); }
double capv(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}
// saturating (Hill, n=1) function
double hillf(double c, double ec50) {
  double cc = posv(c);
  return cc / (cc + ec50 + SAFE);
}
// Thiele-modulus penetration efficiency into an avascular slab.
//   phi = L*sqrt(k_loss/D_eff),  eta = tanh(phi)/phi
// This is the SPATIAL penalty. It is a concentration problem: raising the
// rim concentration (local delivery) raises the core concentration linearly.
double etaf(double L, double kloss, double deff) {
  double phi = posv(L) * sqrt(posv(kloss) / (posv(deff) + SAFE));
  if (phi < 1.0e-6) return 1.0;
  return tanh(phi) / phi;
}
// Killing law. THIS IS THE MODEL.
//   emgd = growth-DEPENDENT Emax  (beta-lactams, glycopeptides, aminoglycosides)
//   emgi = growth-INdependent Emax (rifamycins; partly daptomycin)
// murel = mu/mu_max in the pool being killed. For persisters murel = 0, so a
// drug with emgi = 0 has ZERO kill at ANY concentration. That is an efficacy
// ceiling, not a potency shift, and no dose can move it.
double killf(double conc, double ec50, double emgd, double emgi, double murel) {
  return hillf(conc, ec50) * (emgd * capv(murel, 0.0, 1.0) + emgi);
}

$PARAM @annotated
// ---- host / covariates ---------------------------------------------------
WT       :  70   : body weight (kg)
HOSTF    : 1.0   : host immune competence multiplier (diabetes/steroid -> <1)
PADSEV   : 0.0   : peripheral arterial disease severity 0-1 (global perfusion loss)
EGFR0    :  90   : baseline eGFR (mL/min/1.73m2)

// ---- drug A : growth-dependent backbone (defaults = VANCOMYCIN IV) -------
AV1      :  35   : central volume (L)
AV2      :  40   : peripheral volume (L)
AQ       :  60   : intercompartmental clearance (L/d)
ACL      :  97   : total clearance at EGFR0 (L/d)
AFRREN   : 0.90  : fraction of ACL that is renal
AKA      :  24   : oral absorption rate (1/d)
AFBIO    : 1.00  : oral bioavailability
AFU      : 0.50  : unbound fraction (plasma and tissue assumed equal)
ARBS     : 0.15  : total bone : total plasma AUC ratio
ARIC     : 0.05  : intracellular : extracellular partition
AKPB     : 6.0   : plasma<->bone interstitium equilibration (1/d)
AKLB     : 6.0   : local loss rate from bone interstitium (1/d)
AKDS     : 3.0   : bone<->sequestrum/biofilm equilibration (1/d)
AKIC     : 4.0   : bone<->intracellular equilibration (1/d)
ADEFF    : 0.004 : effective diffusivity in necrotic bone / EPS (cm2/d)
AKLS     : 2.0   : local loss rate inside sequestrum/biofilm (1/d)
AEC50    : 1.0   : free-drug EC50 for kill, ~MIC (mg/L)
AEMGD    : 4.0   : growth-DEPENDENT Emax (1/d)
AEMGI    : 0.03  : growth-INdependent Emax (1/d)  <- vancomycin ~ 0
AINDS    : 0.00  : sensitivity of hepatic CL_A to induction (0 = not a victim)
APLTOX   : 0.0   : platelet toxicity flag (linezolid = 1)
ANEPHTOX : 1.0   : nephrotoxicity flag (vancomycin = 1)

// ---- drug B : growth-INdependent partner (defaults = RIFAMPICIN PO) ------
BV1      :  45   : central volume (L)
BV2      :  35   : peripheral volume (L)
BQ       :  40   : intercompartmental clearance (L/d)
BCL      : 250   : total clearance, uninduced (L/d)
BFRREN   : 0.10  : renal fraction of BCL
BKA      :  20   : oral absorption rate (1/d)
BFBIO    : 0.70  : oral bioavailability
BFU      : 0.20  : unbound fraction
BRBS     : 0.40  : total bone : total plasma AUC ratio
BRIC     : 8.0   : intracellular : extracellular partition
BKPB     : 8.0   : plasma<->bone equilibration (1/d)
BKLB     : 6.0   : loss from bone interstitium (1/d)
BKDS     : 4.0   : bone<->sequestrum equilibration (1/d)
BKIC     : 6.0   : bone<->intracellular equilibration (1/d)
BDEFF    : 0.010 : effective diffusivity (cm2/d)
BKLS     : 2.0   : loss inside sequestrum (1/d)
BEC50    : 0.03  : free EC50 ~ MIC (mg/L)
BEMGD    : 6.0   : growth-dependent Emax (1/d)
BEMGI    : 1.50  : growth-INdependent Emax (1/d)  <- the whole point of rifampicin
BINDS    : 0.35  : autoinduction sensitivity of BCL
BRESOK   : 0.0   : residual activity of drug B against rpoB mutants (0 = none)
BHEPTOX  : 1.0   : hepatotoxicity flag

// ---- induction (rifampicin -> CYP3A4/2C9/P-gp) ---------------------------
KINDOUT  : 0.277 : induction turnover (1/d)  t1/2 = 2.5 d
EMAXIND  : 1.0   : maximal induction signal
EC50IND  : 1.0   : plasma conc of drug B for half-maximal induction (mg/L)

// ---- local antibiotic depot (PMMA cement / calcium sulfate beads) --------
KELUT    : 0.10  : first-order elution from the solid depot (1/d)
VLOC     : 0.05  : wound-fluid volume the depot elutes into (L)
KLOCOUT  : 4.0   : washout of wound fluid (1/d)
KLOCSEQ  : 2.0   : wound fluid -> biofilm/sequestrum equilibration (1/d)
KLOCSYS  : 0.15  : fraction of eluted drug reaching systemic circulation (1/d)

// ---- bacteria ------------------------------------------------------------
MUMAX    : 6.0   : maximal specific growth rate (1/d)
FMUPL    : 0.35  : planktonic mu as fraction of MUMAX
FMUBFA   : 0.12  : biofilm-rim mu as fraction of MUMAX
FMUIC    : 0.005 : intracellular mu as fraction of MUMAX (SCVs barely replicate)
BMAXPL   : 1.0e9 : planktonic carrying capacity (CFU)
BICMAX   : 1.0e7 : intracellular carrying capacity (CFU), set by host cell number
BFDENS   : 1.0e8 : biofilm carrying capacity per cm2 of colonisable surface
KATT     : 0.50  : planktonic -> biofilm attachment (1/d)
KDET     : 0.05  : biofilm -> planktonic detachment (1/d)
KDORM    : 0.004 : biofilm-active -> persister (1/d)
KWAKE    : 0.08  : persister -> biofilm-active AWAKENING (1/d)  <- sets T_ster
KINTL    : 0.02  : planktonic -> intracellular internalisation (1/d)
KREL     : 0.05  : intracellular -> planktonic release (1/d)
FMUT     : 1.0e-8: rpoB mutation frequency per division
FITCOST  : 0.10  : fitness cost of the rpoB mutation
FBFEC    : 3.0   : EC50 multiplier inside the biofilm (EPS binding only)
KEPS     : 1.0e-8: EPS production per biofilm CFU per day (arb. mass units)
KEPSDEG  : 0.05  : EPS turnover (1/d)
EPS50    : 1.0   : EPS mass for half-maximal added diffusion length

// ---- geometry of the sanctuary ------------------------------------------
LBF0     : 0.060 : thickness of a MATURE biofilm (cm)
KLEPS    : 0.010 : added diffusion length from EPS (cm)
KLSEQ    : 0.120 : added diffusion length per cm of sequestrum radius (cm)
AIMP0    : 30.0  : implant surface area at baseline (cm2)

// ---- immune -------------------------------------------------------------
KNEUMAX  : 3.0   : maximal neutrophil-mediated kill of planktonic cells (1/d)
FIMBF    : 0.06  : relative immune access to a MATURE biofilm
FIMBFX   : 0.85  : relative immune access to an IMMATURE (freshly debrided) biofilm
FIMBFPX  : 0.30  : relative immune access to persisters in an immature biofilm
EPSIM50  : 0.10  : EPS mass at which immune evasion is half-restored
AIMP50   : 5.0   : implant area giving half-maximal LOCAL GRANULOCYTE DEFECT (cm2)
FIMIC    : 0.00  : relative immune access to intracellular bacteria (none)
KICCLR   : 0.12  : immune clearance of INFECTED HOST CELLS, relative to KNEUMAX
KBPAMP   : 1.0e7 : bacterial load for half-maximal innate stimulation (CFU)
KINNEU   : 0.5   : neutrophil influx baseline (1/d)
KOUTNEU  : 0.5   : neutrophil egress/death (1/d)
ENEU     : 4.0   : maximal fold recruitment of neutrophils
KINMAC   : 0.2   : macrophage influx (1/d)
KOUTMAC  : 0.2   : macrophage turnover (1/d)
EMAC     : 2.0   : maximal fold macrophage expansion
KDIL1    : 4.0   : IL-1b elimination (1/d)
EIL1     : 12.0  : maximal IL-1b induction
KDTNF    : 6.0   : TNF elimination (1/d)
ETNF     : 8.0   : maximal TNF induction
KDIL6    : 4.0   : IL-6 elimination (1/d)
EIL6     : 20.0  : maximal IL-6 induction
WIL1IL6  : 0.30  : IL-1b contribution to IL-6

// ---- acute-phase biomarkers --------------------------------------------
KDCRP    : 0.876 : CRP elimination (1/d)  t1/2 = 19 h
CRPMAX   : 200   : maximal CRP (mg/L)
CRPBASE  : 2.0   : CRP floor (mg/L)
EC50CRP  : 40.0  : IL-6 for half-maximal CRP
KDESR    : 0.0693: ESR turnover (1/d)  t1/2 = 10 d  <- deliberately slow
ESRMAX   : 110   : maximal ESR (mm/h)
ESRBASE  : 8.0   : ESR floor (mm/h)
EC50ESR  : 12.0  : IL-6 for half-maximal ESR

// ---- bone remodelling / sequestrum -------------------------------------
OB0      : 100   : osteoblast reference
KOUTOB   : 0.05  : osteoblast turnover (1/d)
KAPOP    : 0.50  : S. aureus-induced osteoblast apoptosis (1/d)
KAPO50   : 1.0e6 : bacterial load for half-maximal apoptosis (CFU)
KSRANKL  : 2.0   : RANKL production (1/d)
KDRANKL  : 2.0   : RANKL turnover (1/d)
WIL1R    : 0.25  : IL-1b -> RANKL
WTNFR    : 0.15  : TNF -> RANKL
KSOPG    : 1.0   : OPG production scaled to osteoblast pool (1/d)
KDOPG    : 1.0   : OPG turnover (1/d)
OC0      : 10.0  : osteoclast reference
RR0      : 1.33  : healthy RANKL/OPG ratio
EOC      : 2.5   : maximal fold osteoclast expansion
KOC50    : 20.0  : excess RANKL/OPG ratio for half-maximal expansion
KOUTOC   : 0.05  : osteoclast turnover (1/d)
BM0      : 100   : viable cortical bone reference
BMMAX    : 100   : maximal bone mass
KFORM    : 0.10  : bone formation (1/d)
KRES     : 0.10  : osteoclastic resorption (1/d)
KNECR    : 0.010 : ischaemic necrosis rate constant (1/d)
PERFCRIT : 0.50  : perfusion below which bone becomes necrotic
KSEQF    : 0.004 : sequestrum formation per unit necrotic bone (cm3/d)
KSEQP    : 0.010 : sequestrum formation driven by pus/pressure (cm3/d)
KSEQR    : 0.020 : sequestrum resorption, perfusion-gated (1/d)
SEQ50    : 2.0   : sequestrum volume for half-maximal effects (cm3)

// ---- perfusion (the vicious cycle lives here) --------------------------
KPERF    : 0.15  : perfusion relaxation rate (1/d)
WDPUS    : 0.30  : pus/pressure weight in perfusion damage
WDSEQ    : 0.60  : sequestrum weight in perfusion damage
WDOC     : 0.30  : osteoclastic/inflammatory weight in perfusion damage
PERFMIN  : 0.02  : floor on perfusion
KPUS     : 2.0   : pus formation (mL/d at full stimulus)
KDRAIN   : 0.10  : spontaneous pus drainage (1/d)
PUS50    : 3.0   : pus volume for half-maximal effects (mL)

// ---- surgery (basin change) -------------------------------------------
TSURG    : 1.0e6 : time of index surgery (d); large = no surgery
SURGDUR  : 0.25  : duration over which the surgical effect is applied (d)
DEBLOG   : 0.0   : log10 reduction of biofilm+persisters achieved
DEBSEQF  : 0.0   : fraction of sequestrum excised (0-1)
DEBPUSF  : 0.0   : fraction of pus drained (0-1)
DEBICF   : 0.0   : fraction of intracellular reservoir removed with the bone
IMPREM   : 0.0   : 1 = implant removed at TSURG (A_imp -> 0)
FLAPB    : 0.0   : perfusion added by muscle/free flap (0-0.35)
TREVASC  : 1.0e6 : time of revascularisation (d)
REVASCB  : 0.0   : PADSEV reduction achieved by revascularisation (0-1)

// ---- therapy bookkeeping ---------------------------------------------
TABX0    : 0.0   : antibiotic start (d)
TABX1    : 1.0e6 : antibiotic stop (d)

// ---- toxicity ---------------------------------------------------------
KREPEG   : 0.10  : eGFR recovery rate (1/d)
EMAXAKI  : 0.55  : maximal fractional eGFR loss
AUCTOX50 : 33.0  : rolling exposure of drug A for half-maximal AKI (mg*d/L)
PLT0     : 250   : baseline platelets (10^9/L)
KOUTPLT  : 0.10  : platelet turnover (1/d)
EMAXPLT  : 0.60  : maximal fractional platelet suppression
PLTEC50  : 8.0   : rolling exposure of drug A for half-maximal suppression
TXD50    : 28.0  : days of therapy for half-maximal platelet effect
ALT0     : 25    : baseline ALT (U/L)
KOUTALT  : 0.15  : ALT turnover (1/d)
EMAXALT  : 3.0   : maximal fold ALT rise
ALTEC50  : 4.0   : rolling exposure of drug B for half-maximal ALT rise

// ---- pain / relapse ---------------------------------------------------
KPAIN    : 0.50  : pain relaxation rate (1/d)
WPPUS    : 0.35  : pus weight in pain
WPIL1    : 0.35  : IL-1b weight in pain
WPSEQ    : 0.30  : sequestrum weight in pain
IL150    : 10.0  : IL-1b for half-maximal pain
PREGROW  : 0.35  : probability that ONE surviving CFU re-establishes infection
KHAZ     : 0.02  : reservoir-burden accumulation constant

$CMT @annotated
// -- drug A -------------------------------------------------------------
GUTA  : drug A oral depot (mg)
CENA  : drug A central (mg)
PERA  : drug A peripheral (mg)
CBA   : drug A bone interstitium, total conc (mg/L)
CSA   : drug A sequestrum/biofilm, total conc (mg/L)
CIA   : drug A intracellular, total conc (mg/L)
// -- drug B -------------------------------------------------------------
GUTB  : drug B oral depot (mg)
CENB  : drug B central (mg)
PERB  : drug B peripheral (mg)
CBB   : drug B bone interstitium, total conc (mg/L)
CSB   : drug B sequestrum/biofilm, total conc (mg/L)
CIB   : drug B intracellular, total conc (mg/L)
// -- local depot --------------------------------------------------------
ALOC  : solid local antibiotic depot, drug A (mg)
CLOC  : wound-fluid concentration of drug A (mg/L)
// -- induction ----------------------------------------------------------
EIND  : enzyme/transporter induction state (0-1)
// -- bacteria -----------------------------------------------------------
BPL   : planktonic bacteria in vascularised bone (CFU)
BBFA  : biofilm, metabolically active rim (CFU)
BBFP  : biofilm persisters / dormant core (CFU)
BIC   : intracellular bacteria + SCV in osteoblasts/osteocytes (CFU)
BRES  : rpoB (rifampicin-resistant) subpopulation (CFU)
EPSM  : biofilm EPS matrix mass (arb.)
// -- immune -------------------------------------------------------------
NEU   : neutrophils in bone (relative)
MACR  : macrophages in bone (relative)
IL1   : IL-1beta (relative)
IL6   : IL-6 (relative)
TNFA  : TNF-alpha (relative)
// -- biomarkers ---------------------------------------------------------
CRP   : C-reactive protein (mg/L)
ESR   : erythrocyte sedimentation rate (mm/h)
// -- bone / structure ---------------------------------------------------
OB    : osteoblasts (relative)
OC    : osteoclasts (relative)
RKL   : RANKL (relative)
OPG   : osteoprotegerin (relative)
BM    : viable cortical bone mass (relative)
SEQ   : sequestrum / necrotic bone volume (cm3)
PERF  : perfusion fraction of the infected bone (0-1)
PUS   : abscess volume (mL)
AIMP  : colonisable implant surface area (cm2)
// -- toxicity -----------------------------------------------------------
EGFRC : eGFR (mL/min/1.73m2)
PLT   : platelets (10^9/L)
ALT   : ALT (U/L)
XA    : rolling 1-day exposure of drug A (mg*d/L)
XB    : rolling 1-day exposure of drug B (mg*d/L)
// -- clinical -----------------------------------------------------------
PAIN  : pain NRS 0-10
TXD   : cumulative days of antibiotic therapy (d)
HAZ   : cumulative reservoir-burden index
AUCA  : cumulative AUC of drug A (mg*d/L)
AUCB  : cumulative AUC of drug B (mg*d/L)

$MAIN
F_GUTA = AFBIO;
F_GUTB = BFBIO;

$ODE
// =========================================================================
// 0. surgical / adjuvant switches (applied over a short window at TSURG)
// =========================================================================
double insurg = (SOLVERTIME >= TSURG && SOLVERTIME < TSURG + SURGDUR) ? 1.0 : 0.0;
double kdeb   = insurg * (DEBLOG * log(10.0) / SURGDUR);                 // 1/d
double kseqx  = insurg * (-log(1.0 - capv(DEBSEQF, 0.0, 0.999)) / SURGDUR);
double kpusx  = insurg * (-log(1.0 - capv(DEBPUSF, 0.0, 0.999)) / SURGDUR);
double kicx   = insurg * (-log(1.0 - capv(DEBICF , 0.0, 0.999)) / SURGDUR);
double kimpx  = insurg * IMPREM * (log(1.0e4) / SURGDUR);
double flapon = (SOLVERTIME >= TSURG) ? FLAPB : 0.0;
double padnow = PADSEV * ((SOLVERTIME >= TREVASC) ? (1.0 - REVASCB) : 1.0);
double ontx   = (SOLVERTIME >= TABX0 && SOLVERTIME < TABX1) ? 1.0 : 0.0;

// =========================================================================
// 1. geometry of the sanctuary -> the SPATIAL penalty
//    L grows with EPS and with sequestrum radius, so the disease lengthens
//    its own diffusion path. This is the first half of the vicious cycle.
// =========================================================================
// Diffusion length is NOT a constant floor plus geometry. A freshly debrided,
// EPS-poor film is thin, so debridement improves penetration as well as
// reducing cell number. LBF0 is therefore the thickness of a MATURE film.
double leff  = LBF0 * hillf(posv(EPSM), EPS50)
             + KLEPS * hillf(posv(EPSM), EPS50)
             + KLSEQ * pow(posv(SEQ), 1.0/3.0);
double etaA  = etaf(leff, AKLS, ADEFF);
double etaB  = etaf(leff, BKLS, BDEFF);

// =========================================================================
// 2. perfusion -> gates BOTH drug delivery AND neutrophil access
// =========================================================================
double pnorm = capv(PERF, PERFMIN, 1.0);

// =========================================================================
// 3. PK of drug A (growth-dependent backbone)
// =========================================================================
double CA    = posv(CENA) / AV1;
double clA   = ACL * (AFRREN * capv(EGFRC / EGFR0, 0.05, 2.0)
                    + (1.0 - AFRREN) * (1.0 + AINDS * EIND));
double abtgt = ARBS * (1.0 + AKLB / AKPB);     // makes bone:plasma == ARBS at PERF = 1

dxdt_GUTA = -AKA * GUTA;
dxdt_CENA =  AKA * GUTA - (clA / AV1) * CENA
             - AQ * (CENA / AV1 - PERA / AV2)
             + KLOCSYS * CLOC * VLOC;
dxdt_PERA =  AQ * (CENA / AV1 - PERA / AV2);
dxdt_CBA  =  AKPB * pnorm * abtgt * CA - (AKPB * pnorm + AKLB) * CBA;
dxdt_CSA  =  AKDS * (etaA * CBA - CSA)
             + KLOCSEQ * (etaA * CLOC - CSA);   // local depot bypasses perfusion
dxdt_CIA  =  AKIC * (ARIC * CBA - CIA);

// =========================================================================
// 4. PK of drug B (rifamycin-like partner) + enzyme induction
//    Rifampicin GIVES growth-independent kill and TAKES partner exposure.
// =========================================================================
double CB    = posv(CENB) / BV1;
double clB   = BCL * (BFRREN + (1.0 - BFRREN) * (1.0 + BINDS * EIND));
double bbtgt = BRBS * (1.0 + BKLB / BKPB);

dxdt_GUTB = -BKA * GUTB;
dxdt_CENB =  BKA * GUTB - (clB / BV1) * CENB - BQ * (CENB / BV1 - PERB / BV2);
dxdt_PERB =  BQ * (CENB / BV1 - PERB / BV2);
dxdt_CBB  =  BKPB * pnorm * bbtgt * CB - (BKPB * pnorm + BKLB) * CBB;
dxdt_CSB  =  BKDS * (etaB * CBB - CSB);
dxdt_CIB  =  BKIC * (BRIC * CBB - CIB);

dxdt_EIND = KINDOUT * (EMAXIND * hillf(CB, EC50IND) - EIND);

// =========================================================================
// 5. local depot
// =========================================================================
dxdt_ALOC = -KELUT * ALOC;
dxdt_CLOC =  KELUT * posv(ALOC) / VLOC - (KLOCOUT + KLOCSYS) * CLOC;

// =========================================================================
// 6. free concentrations seen by each bacterial pool
//    plankton  -> bone interstitium     biofilm -> sequestrum core
//    intracell -> intracellular         rpoB    -> bone interstitium
//    (the resistant mutant must GROW to matter, and growing means living in
//     the nutrient-rich, better-perfused rim -> it sees the BEST exposure.
//     That is exactly why a partner drug can hold it down.)
// =========================================================================
double fA_pl = AFU * posv(CBA);
double fA_bf = AFU * posv(CSA);
double fA_ic = AFU * posv(CIA);
double fB_pl = BFU * posv(CBB);
double fB_bf = BFU * posv(CSB);
double fB_ic = BFU * posv(CIB);

// =========================================================================
// 7. growth rates per pool -> the PHENOTYPIC penalty
// =========================================================================
double btot  = posv(BPL) + posv(BBFA) + posv(BBFP) + posv(BIC) + posv(BRES);
double aseq  = 6.0 * pow(posv(SEQ), 2.0/3.0);        // cm2 of necrotic surface
double bfcap = BFDENS * (posv(AIMP) + aseq) + 1.0;
double bbf   = posv(BBFA) + posv(BBFP) + posv(BRES);

// The planktonic niche is consumed by planktonic cells only. Biofilm cells sit
// on a surface and do not compete for it. Charging biofilm against BMAXPL
// silently freezes the resistant mutant (which lives planktonically) and hides
// rifampicin monotherapy failure - the single most important thing this model
// has to say about drug choice.
double bplk  = posv(BPL) + posv(BRES);
double mupl  = MUMAX * FMUPL  * posv(1.0 - bplk / BMAXPL);
double fmat   = 1.0 - hillf(posv(EPSM), EPSIM50);   // 1 = immature, 0 = mature
double fmubfe = FMUBFA + (FMUPL - FMUBFA) * fmat;
double mubfa = MUMAX * fmubfe * posv(1.0 - bbf  / bfcap);
// the intracellular niche is bounded by the host cells that host it: as
// osteoblasts die the niche shrinks. Without this cap BIC grows without limit
// once therapy stops, which is both unphysiological and numerically fatal.
double icmax = BICMAX * (0.2 + 0.8 * capv(posv(BM) / BM0, 0.0, 1.0));
double muic  = MUMAX * FMUIC  * posv(1.0 - posv(BIC) / icmax);
double mures = mupl * (1.0 - FITCOST);

double rpl   = mupl  / MUMAX;
double rbfa  = mubfa / MUMAX;
double ric   = muic  / MUMAX;
double rres  = mures / MUMAX;

// =========================================================================
// 8. killing  (the two penalties multiply here)
// =========================================================================
double ecAbf = AEC50 * FBFEC;      // modest EPS-binding potency shift only
double ecBbf = BEC50 * FBFEC;

double kA_pl  = killf(fA_pl, AEC50, AEMGD, AEMGI, rpl );
double kA_bfa = killf(fA_bf, ecAbf, AEMGD, AEMGI, rbfa);
double kA_bfp = killf(fA_bf, ecAbf, AEMGD, AEMGI, 0.0 );
double kA_ic  = killf(fA_ic, AEC50, AEMGD, AEMGI, ric );
double kA_res = killf(fA_pl, AEC50, AEMGD, AEMGI, rres);

double kB_pl  = killf(fB_pl, BEC50, BEMGD, BEMGI, rpl );
double kB_bfa = killf(fB_bf, ecBbf, BEMGD, BEMGI, rbfa);
double kB_bfp = killf(fB_bf, ecBbf, BEMGD, BEMGI, 0.0 );
double kB_ic  = killf(fB_ic, BEC50, BEMGD, BEMGI, ric );
double kB_res = BRESOK * killf(fB_pl, BEC50, BEMGD, BEMGI, rres);

// immune killing is perfusion-gated: losing perfusion loses the host defence
// at the same time as it loses the drug. Two hits, one cause.
double pamp  = hillf(btot, KBPAMP);
double kimm  = KNEUMAX * (posv(NEU) / (KINNEU / KOUTNEU)) * pnorm * HOSTF;

// Immune ACCESS to the biofilm is not a constant. Two things gate it:
//  (a) biofilm MATURITY - a freshly debrided, EPS-poor biofilm is largely
//      accessible; a mature one is not. This is how debridement cures: it
//      resets maturity, not just cell number.
//  (b) the FOREIGN BODY. An implant surface causes a local granulocyte defect
//      (Zimmerli 1984 J Clin Invest, PMID 6323536; tissue-cage model described
//      in Zimmerli 1982 J Infect Dis, PMID 7119479) that does NOT go away
//      when the biofilm is cleaned off. This single term is why DAIR without
//      rifampicin fails while the same debridement with the implant REMOVED
//      succeeds - the mathematics of Zimmerli 1998 JAMA.
double fimb  = (1.0 - hillf(posv(EPSM), EPSIM50))
             * (1.0 - hillf(posv(AIMP), AIMP50));
double fimbfe = FIMBF + (FIMBFX - FIMBF) * fimb;
double fimbpe = FIMBFPX * fimb;

// =========================================================================
// 9. bacterial subpopulations
// =========================================================================
double surfav = capv((posv(AIMP) + aseq) / (AIMP0 + 6.0 * pow(SEQ50, 2.0/3.0)),
                     0.0, 3.0);

// every state is read through posv(): the pools fall many logs under therapy
// and a tiny negative numerical excursion fed back through its OWN growth term
// would otherwise diverge instead of decaying.
double bpl = posv(BPL);
double bfa = posv(BBFA);
double bfp = posv(BBFP);
double bic = posv(BIC);
double brs = posv(BRES);
double epm = posv(EPSM);

dxdt_BPL  =  mupl * bpl * extf(bpl)
           - KATT * surfav * bpl + KDET * bfa
           - KINTL * bpl + KREL * bic
           - (kA_pl + kB_pl + kimm) * bpl
           - FMUT * mupl * bpl
           - kdeb * bpl;

dxdt_BBFA =  mubfa * bfa * extf(bfa)
           + KATT * surfav * bpl - KDET * bfa
           - KDORM * bfa + KWAKE * bfp
           - (kA_bfa + kB_bfa + kimm * fimbfe) * bfa
           - FMUT * mubfa * bfa
           - kdeb * bfa;

dxdt_BBFP =  KDORM * bfa - KWAKE * bfp
           - (kA_bfp + kB_bfp + kimm * fimbpe) * bfp
           - kdeb * bfp;

// The intracellular niche is a PERSISTENCE niche, not a growth niche: SCVs
// barely replicate (FMUIC), they drain outwards (KREL) into a space where drug
// and neutrophils work, and the infected host cell itself is killed by the
// immune system (KICCLR). Neutrophils cannot reach the bacterium (FIMIC = 0)
// but they can remove the cell holding it - and that removal is perfusion-
// gated like everything else here, so a poorly perfused focus keeps it.
dxdt_BIC  =  muic * bic * extf(bic) + KINTL * bpl - KREL * bic
           - (kA_ic + kB_ic + kimm * (FIMIC + KICCLR)) * bic
           - kicx * bic;

dxdt_BRES =  FMUT * (mupl * bpl + mubfa * bfa)
           + mures * brs * extf(brs)
           - (kA_res + kB_res + kimm) * brs
           - kdeb * brs;

dxdt_EPSM =  KEPS * bfa - KEPSDEG * epm - kdeb * epm;

// =========================================================================
// 10. innate immune response
// =========================================================================
double neuss = (KINNEU / KOUTNEU) * (1.0 + ENEU * pamp * HOSTF);
double macss = (KINMAC / KOUTMAC) * (1.0 + EMAC * pamp);
dxdt_NEU  = KOUTNEU * (neuss - NEU);
dxdt_MACR = KOUTMAC * (macss - MACR);

double macrel = posv(MACR) / (KINMAC / KOUTMAC);
dxdt_IL1  = KDIL1 * (1.0 + EIL1 * pamp * macrel - IL1);
dxdt_TNFA = KDTNF * (1.0 + ETNF * pamp * macrel - TNFA);
dxdt_IL6  = KDIL6 * (1.0 + EIL6 * pamp * macrel + WIL1IL6 * posv(IL1) - IL6);

// =========================================================================
// 11. acute-phase biomarkers — same event, two different clocks
//     CRP t1/2 19 h (fast, useful) vs ESR t1/2 10 d (slow, misleading)
// =========================================================================
dxdt_CRP = KDCRP * (CRPBASE + CRPMAX * hillf(IL6, EC50CRP) - CRP);
dxdt_ESR = KDESR * (ESRBASE + ESRMAX * hillf(IL6, EC50ESR) - ESR);

// =========================================================================
// 12. bone remodelling, necrosis, sequestrum
// =========================================================================
double bload = posv(BPL) + posv(BIC);
dxdt_OB  = KOUTOB * (OB0 - OB) - KAPOP * OB * hillf(bload, KAPO50);
dxdt_RKL = KSRANKL * (1.0 + WIL1R * posv(IL1) + WTNFR * posv(TNFA)) - KDRANKL * RKL;
dxdt_OPG = KSOPG * (posv(OB) / OB0) - KDOPG * OPG;

double rr  = posv(RKL) / (posv(OPG) + 0.05);
// OC is driven by the RANKL/OPG ratio RELATIVE to its healthy value, and is
// bounded: an unbounded ratio (OPG -> 0 as osteoblasts die) would otherwise
// send OC to thousands and make the bone equation numerically explosive.
dxdt_OC  = KOUTOC * (OC0 * (1.0 + EOC * hillf(posv(rr / RR0 - 1.0), KOC50)) - OC);

double ischae = posv(1.0 - pnorm / PERFCRIT);
dxdt_BM  = KFORM * (posv(OB) / OB0) * posv(1.0 - BM / BMMAX)
         - KRES * (posv(OC) / OC0) * (posv(BM) / BM0)
         - KNECR * ischae * posv(BM);
dxdt_SEQ = KSEQF * posv(BM) * ischae
         + KSEQP * hillf(PUS, PUS50) * (posv(BM) / BM0)
         - KSEQR * posv(SEQ) * pnorm
         - kseqx * posv(SEQ);

// =========================================================================
// 13. perfusion — the hinge of the whole model
// =========================================================================
double dmg = WDPUS * hillf(PUS, PUS50)
           + WDSEQ * hillf(SEQ, SEQ50)
           + WDOC  * hillf(posv(posv(OC) / OC0 - 1.0), 1.5);
double perfss = capv((1.0 - capv(dmg, 0.0, 0.95)) * (1.0 - capv(padnow, 0.0, 0.9))
                     + flapon, PERFMIN, 1.0);
dxdt_PERF = KPERF * (perfss - PERF);

dxdt_PUS  = KPUS * (posv(NEU) / (KINNEU / KOUTNEU)) * pamp
          - KDRAIN * posv(PUS) - kpusx * posv(PUS);
dxdt_AIMP = -kimpx * posv(AIMP);

// =========================================================================
// 14. exposure accumulators and toxicity
// =========================================================================
dxdt_XA = CA - XA;                    // ~ rolling 1-day AUC of drug A
dxdt_XB = CB - XB;
dxdt_AUCA = CA;
dxdt_AUCB = CB;

double egfrt = EGFR0 * (1.0 - ANEPHTOX * EMAXAKI * hillf(XA, AUCTOX50));
dxdt_EGFRC = KREPEG * (egfrt - EGFRC);

double pltt = PLT0 * (1.0 - APLTOX * EMAXPLT * hillf(XA, PLTEC50) * hillf(TXD, TXD50));
dxdt_PLT   = KOUTPLT * (pltt - PLT);

double altt = ALT0 * (1.0 + BHEPTOX * EMAXALT * hillf(XB, ALTEC50));
dxdt_ALT   = KOUTALT * (altt - ALT);

dxdt_TXD = ontx;

// =========================================================================
// 15. pain and reservoir burden
// =========================================================================
double painss = 10.0 * (WPPUS * hillf(PUS, PUS50)
                      + WPIL1 * hillf(IL1, IL150)
                      + WPSEQ * hillf(SEQ, SEQ50)) / (WPPUS + WPIL1 + WPSEQ);
dxdt_PAIN = KPAIN * (painss - PAIN);
dxdt_HAZ  = KHAZ * log10(1.0 + posv(BBFP) + posv(BIC) + posv(BRES));

$TABLE
// ---- observable concentrations ------------------------------------------
double CAOBS = posv(CENA) / AV1;
double CBOBS = posv(CENB) / BV1;
double LEFFO = (LBF0 + KLEPS) * hillf(posv(EPSM), EPS50)
             + KLSEQ * pow(posv(SEQ), 1.0/3.0);
double PHIA  = LEFFO * sqrt(AKLS / (ADEFF + 1.0e-12));
double PHIB  = LEFFO * sqrt(BKLS / (BDEFF + 1.0e-12));
double ETAA  = etaf(LEFFO, AKLS, ADEFF);
double ETAB  = etaf(LEFFO, BKLS, BDEFF);
double BSRA  = (CAOBS > 1.0e-6) ? (posv(CBA) / CAOBS) : 0.0;   // apparent bone:serum

// ---- bacterial burden on log10 scale -----------------------------------
double BTOT  = posv(BPL) + posv(BBFA) + posv(BBFP) + posv(BIC) + posv(BRES);
double LBTOT = log10(BTOT  + 1.0e-6);
double LBPL  = log10(posv(BPL)  + 1.0e-6);
double LBFA  = log10(posv(BBFA) + 1.0e-6);
double LBFP  = log10(posv(BBFP) + 1.0e-6);
double LBIC  = log10(posv(BIC)  + 1.0e-6);
double LBRES = log10(posv(BRES) + 1.0e-6);

// ---- THE reservoir: what relapse is made of ----------------------------
double RESV  = posv(BBFP) + posv(BIC) + posv(BRES) + posv(BBFA);
double PRLPS = 1.0 - exp(-PREGROW * RESV);          // P(>=1 CFU re-establishes)

// ---- RSTER : is sterilisation even POSSIBLE with this regimen? ----------
//      RSTER = achievable biofilm kill / (biofilm growth + persister feed)
//      <=1 -> a stable non-zero equilibrium in bone. Duration is irrelevant.
// RSTER must be judged against the INTRINSIC biofilm growth rate, not the
// density-limited one. At carrying capacity the logistic term goes to zero, so
// a density-limited denominator makes RSTER blow up precisely where nothing is
// being killed. The question RSTER answers is "as this population falls away
// from capacity, does kill still outrun replication?" - that needs mu_intrinsic.
double FMATO  = 1.0 - hillf(posv(EPSM), EPSIM50);
double MUBF   = MUMAX * (FMUBFA + (FMUPL - FMUBFA) * FMATO);
double RBF    = MUBF / MUMAX;
double KILLBF = killf(AFU * posv(CSA), AEC50 * FBFEC, AEMGD, AEMGI, RBF)
              + killf(BFU * posv(CSB), BEC50 * FBFEC, BEMGD, BEMGI, RBF);
double KILLPS = killf(AFU * posv(CSA), AEC50 * FBFEC, AEMGD, AEMGI, 0.0)
              + killf(BFU * posv(CSB), BEC50 * FBFEC, BEMGD, BEMGI, 0.0);
double PNORMO = capv(PERF, PERFMIN, 1.0);
double FIMBO  = (1.0 - hillf(posv(EPSM), EPSIM50))
              * (1.0 - hillf(posv(AIMP), AIMP50));
double KIMMO  = KNEUMAX * (posv(NEU) / (KINNEU / KOUTNEU)) * PNORMO * HOSTF;
double KIMMBF = KIMMO * (FIMBF + (FIMBFX - FIMBF) * FIMBO);
double KIMMPS = KIMMO * FIMBFPX * FIMBO;
double PFEED  = KWAKE * posv(BBFP) / (posv(BBFA) + 1.0);
double RSTER  = (KILLBF + KIMMBF) / (MUBF + PFEED + 1.0e-9);

// ---- T_ster : how long, IF RSTER > 1 -----------------------------------
//      persisters leave only by being killed (KILLPS) or by waking (KWAKE)
double KEFF   = KILLPS + KIMMPS + KWAKE;
double TSTER  = (posv(BBFP) > 1.0) ? (log(posv(BBFP)) / (KEFF + 1.0e-9)) : 0.0;

// ---- misc clinical -----------------------------------------------------
double AUC24A = 24.0 * posv(XA);      // mg*h/L, comparable to the vanco target
double AUC24B = 24.0 * posv(XB);
double BMLOSS = 100.0 * (1.0 - posv(BM) / BM0);

$CAPTURE @annotated
CAOBS  : drug A plasma concentration (mg/L)
CBOBS  : drug B plasma concentration (mg/L)
LEFFO  : effective diffusion length (cm)
PHIA   : Thiele modulus, drug A
PHIB   : Thiele modulus, drug B
ETAA   : penetration efficiency, drug A
ETAB   : penetration efficiency, drug B
BSRA   : apparent bone:serum ratio, drug A
LBTOT  : log10 total bacterial burden
LBPL   : log10 planktonic
LBFA   : log10 biofilm-active
LBFP   : log10 persisters
LBIC   : log10 intracellular
LBRES  : log10 rpoB-resistant
RESV   : reservoir size (CFU)
PRLPS  : probability of relapse from the reservoir
MUBF   : intrinsic biofilm growth rate (1/d)
RSTER  : sterilisation ratio (>1 required)
KILLBF : total kill rate on biofilm rim (1/d)
KILLPS : total DRUG kill rate on persisters (1/d)
KIMMBF : immune kill rate on the biofilm rim (1/d)
KIMMPS : immune kill rate on persisters (1/d)
TSTER  : time to deplete the persister pool (d)
AUC24A : drug A AUC24 (mg*h/L)
AUC24B : drug B AUC24 (mg*h/L)
BMLOSS : cortical bone loss (%)
'

mod <- mcode("com_qsp", com_code, atol = 1e-8, rtol = 1e-8, maxsteps = 500000)

## ============================================================================
##  Initial conditions — day 0 is taken as the "moment of presentation" with chronic osteomyelitis.
##  (This is not a simulation of the onset of infection. It starts from a state in which
##   sequestrum, biofilm, abscess and hypoperfusion are already established.)
## ============================================================================

com_init <- c(
  GUTA = 0, CENA = 0, PERA = 0, CBA = 0, CSA = 0, CIA = 0,
  GUTB = 0, CENB = 0, PERB = 0, CBB = 0, CSB = 0, CIB = 0,
  ALOC = 0, CLOC = 0, EIND = 0,
  BPL  = 1e6,    # planktonic bacteria
  BBFA = 1e8,    # biofilm active rim
  BBFP = 5e6,    # persisters (about 5%)
  BIC  = 1e5,    # intracellular · SCV
  BRES = 10,     # rpoB-resistant subpopulation — not 0 even before treatment.
                 # At mutation-selection balance the frequency is ~ FMUT/FITCOST = 1e-7,
                 # so with 1e8 organisms the expectation is 10. Rifampicin alone fails
                 # not because resistance 'arises' but because it is already 'there'.
  EPSM = 20.0,   # equilibrium value of KEPS*BBFA/KEPSDEG
  NEU  = 3.0, MACR = 2.0, IL1 = 20, IL6 = 30, TNFA = 12,
  CRP  = 85, ESR = 80,
  OB   = 20, OC = 24, RKL = 8, OPG = 0.2,
  BM   = 100, SEQ = 2.0, PERF = 0.37, PUS = 3.0, AIMP = 30,
  EGFRC = 90, PLT = 250, ALT = 25, XA = 0, XB = 0,
  PAIN = 6, TXD = 0, HAZ = 0, AUCA = 0, AUCB = 0
)
mod <- update(mod, init = as.list(com_init))

CMTI <- setNames(seq_along(names(init(mod))), names(init(mod)))

## ============================================================================
##  Drug library
##  "The axis along which a drug is chosen in osteomyelitis is not the MIC but the product
##       (bone penetration) x (intracellular reach) x (kill of non-growing organisms)
##   of those three." — the last three columns of the table below are that axis.
## ============================================================================

drug_A_library <- list(

  vancomycin = list(
    AV1 = 35, AV2 = 40, AQ = 60, ACL = 97, AFRREN = 0.90, AKA = 24, AFBIO = 0.00,
    AFU = 0.50, ARBS = 0.15, ARIC = 0.05,
    ADEFF = 0.004, AKLS = 2.0,
    AEC50 = 1.0, AEMGD = 4.0, AEMGI = 0.03,
    AINDS = 0.00, APLTOX = 0.0, ANEPHTOX = 1.0),

  nafcillin = list(     # standard MSSA backbone (cefazolin treated almost identically)
    AV1 = 25, AV2 = 20, AQ = 40, ACL = 480, AFRREN = 0.30, AKA = 24, AFBIO = 0.00,
    AFU = 0.10, ARBS = 0.20, ARIC = 0.05,
    ADEFF = 0.010, AKLS = 2.0,
    AEC50 = 0.5, AEMGD = 6.0, AEMGI = 0.02,
    AINDS = 0.00, APLTOX = 0.0, ANEPHTOX = 0.2),

  levofloxacin = list(  # the rifampicin partner of the Zimmerli regimen · oral F ~ 1
    AV1 = 90, AV2 = 40, AQ = 40, ACL = 192, AFRREN = 0.80, AKA = 30, AFBIO = 0.99,
    AFU = 0.70, ARBS = 0.40, ARIC = 5.0,
    ADEFF = 0.020, AKLS = 1.5,
    AEC50 = 0.5, AEMGD = 9.0, AEMGI = 0.10,
    AINDS = 0.30, APLTOX = 0.0, ANEPHTOX = 0.0),

  linezolid = list(     # bacteriostatic — Emax_gd never exceeds mu
    AV1 = 45, AV2 = 25, AQ = 40, ACL = 168, AFRREN = 0.30, AKA = 40, AFBIO = 1.00,
    AFU = 0.70, ARBS = 0.45, ARIC = 5.0,
    ADEFF = 0.020, AKLS = 1.5,
    AEC50 = 2.0, AEMGD = 2.0, AEMGI = 0.08,
    AINDS = 0.35, APLTOX = 1.0, ANEPHTOX = 0.0),

  daptomycin = list(
    AV1 = 9,  AV2 = 6,  AQ = 10, ACL = 12, AFRREN = 0.80, AKA = 24, AFBIO = 0.00,
    AFU = 0.08, ARBS = 0.15, ARIC = 0.10,
    ADEFF = 0.003, AKLS = 3.0,
    AEC50 = 0.5, AEMGD = 10.0, AEMGI = 0.25,
    AINDS = 0.00, APLTOX = 0.0, ANEPHTOX = 0.3),

  clindamycin = list(
    AV1 = 60, AV2 = 30, AQ = 40, ACL = 300, AFRREN = 0.10, AKA = 30, AFBIO = 0.90,
    AFU = 0.20, ARBS = 0.60, ARIC = 8.0,
    ADEFF = 0.015, AKLS = 2.0,
    AEC50 = 0.25, AEMGD = 2.5, AEMGI = 0.08,
    AINDS = 0.35, APLTOX = 0.0, ANEPHTOX = 0.0)
)

drug_B_library <- list(

  rifampicin = list(    # the only drug with a large Emax_gi. Hence combination is the rule,
                        # monotherapy the contraindication — both from the same equation.
    BV1 = 45, BV2 = 35, BQ = 40, BCL = 250, BFRREN = 0.10, BKA = 20, BFBIO = 0.70,
    BFU = 0.20, BRBS = 0.40, BRIC = 8.0,
    BDEFF = 0.010, BKLS = 2.0,
    BEC50 = 0.03, BEMGD = 6.0, BEMGI = 1.50,
    BINDS = 0.35, BRESOK = 0.0, BHEPTOX = 1.0),

  none = list(          # no partner — makes Emax_gi 0
    BV1 = 45, BV2 = 35, BQ = 40, BCL = 250, BFRREN = 0.10, BKA = 20, BFBIO = 0.70,
    BFU = 0.20, BRBS = 0.40, BRIC = 8.0,
    BDEFF = 0.010, BKLS = 2.0,
    BEC50 = 0.03, BEMGD = 0.0, BEMGI = 0.00,
    BINDS = 0.00, BRESOK = 0.0, BHEPTOX = 0.0)
)

set_drugs <- function(m, A = "vancomycin", B = "none") {
  p <- c(drug_A_library[[A]], drug_B_library[[B]])
  do.call(param, c(list(m), p))
}

## ============================================================================
##  Dosing-event helpers
## ============================================================================

ev_iv <- function(amt, ii, days, cmt, start = 0, inf_h = 1) {
  n <- max(1, floor(days / ii))
  data.frame(ID = 1, time = start, amt = amt, ii = ii, addl = n - 1,
             cmt = CMTI[[cmt]], rate = amt / (inf_h / 24), evid = 1)
}
ev_po <- function(amt, ii, days, cmt, start = 0) {
  n <- max(1, floor(days / ii))
  data.frame(ID = 1, time = start, amt = amt, ii = ii, addl = n - 1,
             cmt = CMTI[[cmt]], rate = 0, evid = 1)
}
ev_bolus <- function(amt, cmt, time = 0) {
  data.frame(ID = 1, time = time, amt = amt, ii = 0, addl = 0,
             cmt = CMTI[[cmt]], rate = 0, evid = 1)
}
ev_none <- function() {
  data.frame(ID = 1, time = 0, amt = 0, ii = 0, addl = 0,
             cmt = 1, rate = 0, evid = 1)[0, ]
}
ev_bind <- function(...) {
  d <- do.call(rbind, Filter(function(x) nrow(x) > 0, list(...)))
  if (is.null(d) || nrow(d) == 0) return(ev_none())
  d[order(d$time), ]
}

run_scn <- function(m, events, end = 400, delta = 0.25) {
  ## If a dosing time beyond the observation window is left in the data set, lsoda tries to
  ## integrate out to that time and stalls. Trim them beforehand.
  if (nrow(events) > 0) events <- events[events$time <= end, , drop = FALSE]
  if (nrow(events) == 0) {
    mrgsim(m, end = end, delta = delta, atol = 1e-8, rtol = 1e-8,
           maxsteps = 2000000, hmax = 0.5)
  } else {
    mrgsim_d(m, data = events, end = end, delta = delta, atol = 1e-8,
             rtol = 1e-8, maxsteps = 2000000, hmax = 0.5)
  }
}

## ============================================================================
##  Scenario definitions
##  Each scenario is a combination of (drug · surgery · duration), and the point is that these
##  three axes each touch a different mathematical object: RSTER · P0 · T.
## ============================================================================

## --- S1. Natural history (no treatment) ------------------------------------
##  Expect: bacterial burden at equilibrium, SEQ growing steadily and PERF falling steadily.
##  A self-sustaining state whose vicious-circle gain exceeds 1 — it does not heal by itself.
scn01 <- function() {
  m <- set_drugs(mod, "vancomycin", "none")
  m <- param(m, TSURG = 1e6, TABX0 = 0, TABX1 = 0)
  list(name = "S1 natural history (untreated)", mod = m, ev = ev_none(), end = 400)
}

## --- S2. IV vancomycin 6 weeks, no debridement -----------------------------
##  Expect: CRP normalises and planktonic bacteria disappear (clinically it "looks better").
##  But RSTER << 1, so biofilm and persisters remain almost intact and it relapses after
##  stopping. The purpose of this scenario is to show that "the biomarker recovering" and
##  "the disease recovering" are different state variables.
scn02 <- function() {
  m <- set_drugs(mod, "vancomycin", "none")
  m <- param(m, TSURG = 1e6, TABX0 = 0, TABX1 = 42)
  e <- ev_iv(1250, 1/2, 42, "CENA", inf_h = 1.5)   # 1.25 g q12h x 6 weeks
  list(name = "S2 vancomycin 6 weeks, no debridement", mod = m, ev = e, end = 400)
}

## --- S3. Radical debridement + IV vancomycin 6 weeks (no rifampicin) -------
##  Expect: debridement lowers P0 by 5.5 log, so the residual reservoir shrinks and the
##  remaining persisters disappear only through awakening (KWAKE). Borderline cure at 6 weeks.
scn03 <- function() {
  m <- set_drugs(mod, "vancomycin", "none")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = 5.5, DEBSEQF = 0.95,
             DEBPUSF = 0.95, DEBICF = 0.90, IMPREM = 1, FLAPB = 0.20,
             TABX0 = 0, TABX1 = 44)
  e <- ev_iv(1250, 1/2, 44, "CENA", inf_h = 1.5)
  list(name = "S3 debridement + vancomycin 6 weeks (no rifampicin)", mod = m, ev = e, end = 400)
}

## --- S4. DAIR (implant retained) + levofloxacin + rifampicin 12 weeks ------
##  The Zimmerli 1998 JAMA regimen. The implant stays, so A_imp is maintained and the
##  debridement effect is only 3 log. It succeeds nonetheless for exactly one reason —
##  rifampicin lifts RSTER above 1 and kills persisters directly.
scn04 <- function() {
  m <- set_drugs(mod, "levofloxacin", "rifampicin")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = 3.0, DEBSEQF = 0.60,
             DEBPUSF = 0.95, DEBICF = 0.50, IMPREM = 0, FLAPB = 0.10,
             TABX0 = 0, TABX1 = 86)
  e <- ev_bind(
    ev_po(750, 1,   86, "GUTA"),      # levofloxacin 750 mg qd
    ev_po(450, 1/2, 86, "GUTB")       # rifampicin 450 mg bid
  )
  list(name = "S4 DAIR + levofloxacin/rifampicin 12 weeks", mod = m, ev = e, end = 400)
}

## --- S4b. DAIR + levofloxacin alone 12 weeks (rifampicin removed) ----------
##  The only difference from S4 is rifampicin. It is the control arm of the Zimmerli RCT.
scn04b <- function() {
  m <- set_drugs(mod, "levofloxacin", "none")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = 3.0, DEBSEQF = 0.60,
             DEBPUSF = 0.95, DEBICF = 0.50, IMPREM = 0, FLAPB = 0.10,
             TABX0 = 0, TABX1 = 86)
  e <- ev_po(750, 1, 86, "GUTA")
  list(name = "S4b DAIR + levofloxacin alone 12 weeks", mod = m, ev = e, end = 400)
}

## --- S5. Two-stage exchange + oral switch 6 weeks (OVIVA type) -------------
##  The implant is removed, so A_imp -> 0. Oral linezolid + rifampicin.
##  The route (IV/oral) changes only F and CL, so equal outcomes wherever bioavailability
##  is adequate are, in this model, not a theorem but a structure.
scn05 <- function() {
  m <- set_drugs(mod, "linezolid", "rifampicin")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = 5.0, DEBSEQF = 0.90,
             DEBPUSF = 0.95, DEBICF = 0.80, IMPREM = 1, FLAPB = 0.15,
             TABX0 = 0, TABX1 = 44)
  e <- ev_bind(
    ev_po(600, 1/2, 44, "GUTA"),      # linezolid 600 mg bid
    ev_po(450, 1/2, 44, "GUTB")       # rifampicin 450 mg bid
  )
  list(name = "S5 two-stage exchange + oral linezolid/rifampicin 6 weeks", mod = m, ev = e, end = 400)
}

## --- S6. Rifampicin monotherapy (absolutely contraindicated) ---------------
##  Expect: a steep initial kill -> the rpoB-resistant subpopulation takes over at 10-20 days.
##  P(mutant present) = 1 - exp(-1e-8 * 1e9) ~ 1.0, so this is not probability but
##  arithmetic.
scn06 <- function(debride = FALSE) {
  m <- set_drugs(mod, "vancomycin", "rifampicin")
  m <- param(m, AEMGD = 0, AEMGI = 0,        # backbone switched off = rifampicin alone
             TSURG = if (debride) 2 else 1e6, SURGDUR = 0.25,
             DEBLOG = if (debride) 1.5 else 0, DEBSEQF = if (debride) 0.60 else 0,
             DEBPUSF = if (debride) 0.95 else 0, DEBICF = if (debride) 0.50 else 0,
             IMPREM = 0, TABX0 = 0, TABX1 = 86)
  e <- ev_po(450, 1/2, 86, "GUTB")
  list(name = paste0("S6 rifampicin alone", if (debride) " + limited debridement" else " (no debridement)"),
       mod = m, ev = e, end = 250)
}

## --- S7. Local antibiotic beads alone (no debridement) ---------------------
##  ★ the sharpest result of this model ★
##  Even when the wound-fluid vancomycin concentration reaches ~1000 mg/L (1000-fold the MIC),
##  the biofilm is not sterilised. Because the ceiling is Emax, not concentration.
scn07 <- function() {
  m <- set_drugs(mod, "vancomycin", "none")
  m <- param(m, TSURG = 1e6, TABX0 = 0, TABX1 = 60)
  e <- ev_bolus(2000, "ALOC", time = 0)     # vancomycin 2 g in the cement
  list(name = "S7 local beads alone (fails despite the high concentration)", mod = m, ev = e, end = 300)
}

## --- S8. Diabetic foot bone infection (PAD) ± revascularisation ------------
##  PADSEV cuts PERF systemically, so drug exposure and neutrophil access fall at the same
##  time. In this model revascularisation sits in the same place as "antibiotic therapy".
scn08 <- function(revasc = FALSE) {
  m <- set_drugs(mod, "levofloxacin", "rifampicin")
  m <- param(m, PADSEV = 0.55, HOSTF = 0.65,
             TSURG = 2, SURGDUR = 0.25, DEBLOG = 4.0, DEBSEQF = 0.80,
             DEBPUSF = 0.90, DEBICF = 0.60, IMPREM = 0, FLAPB = 0.05,
             TREVASC = if (revasc) 3 else 1e6, REVASCB = if (revasc) 0.75 else 0,
             TABX0 = 0, TABX1 = 44)
  e <- ev_bind(ev_po(750, 1, 44, "GUTA"), ev_po(450, 1/2, 44, "GUTB"))
  list(name = paste0("S8 diabetic foot bone infection ", if (revasc) "+ revascularisation" else "(PAD uncorrected)"),
       mod = m, ev = e, end = 400)
}

## --- S9. Duration sweep (2/4/6/12/26 weeks) --------------------------------
##  Duration changes only the "time to reach". In a regimen with RSTER > 1 there is a point
##  at which the curve flattens; in a regimen with RSTER < 1 it flattens at no duration at all.
scn09 <- function(weeks, with_rif = TRUE, deblog = 3.0, impkeep = TRUE) {
  d <- weeks * 7
  m <- set_drugs(mod, "levofloxacin", if (with_rif) "rifampicin" else "none")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = deblog, DEBSEQF = 0.70,
             DEBPUSF = 0.95, DEBICF = 0.60, IMPREM = if (impkeep) 0 else 1,
             FLAPB = 0.10, TABX0 = 0, TABX1 = d + 2)
  e <- ev_bind(
    ev_po(750, 1, d, "GUTA"),
    if (with_rif) ev_po(450, 1/2, d, "GUTB") else ev_none()
  )
  list(name = sprintf("S9 %2d wk %s", weeks, if (with_rif) "+rifampicin" else "alone"),
       mod = m, ev = e, end = max(400, d + 250))
}

## --- S10. Long-term suppression (unremovable implant) ----------------------
##  A strategy that holds the equilibrium below threshold even when RSTER < 1. The model shows
##  explicitly that this is not "cure" but a "shift of the equilibrium".
scn10 <- function() {
  m <- set_drugs(mod, "clindamycin", "none")
  m <- param(m, TSURG = 2, SURGDUR = 0.25, DEBLOG = 2.0, DEBSEQF = 0.40,
             DEBPUSF = 0.90, DEBICF = 0.30, IMPREM = 0,
             TABX0 = 0, TABX1 = 730)
  e <- ev_po(600, 1/3, 730, "GUTA")     # clindamycin 600 mg tid, indefinitely
  list(name = "S10 long-term oral suppressive therapy", mod = m, ev = e, end = 900)
}

## --- S11. Timing of debridement (early vs chronic) -------------------------
##  Same drug, same duration, different timing. If the geometry (L) the disease built for
##  itself has already grown, the same regimen works less well.
scn11 <- function(tsurg = 2) {
  m <- set_drugs(mod, "levofloxacin", "rifampicin")
  m <- param(m, TSURG = tsurg, SURGDUR = 0.25, DEBLOG = 4.5, DEBSEQF = 0.85,
             DEBPUSF = 0.95, DEBICF = 0.70, IMPREM = 0, FLAPB = 0.10,
             TABX0 = tsurg, TABX1 = tsurg + 86)
  e <- ev_bind(
    ev_po(750, 1,   86, "GUTA", start = tsurg),
    ev_po(450, 1/2, 86, "GUTB", start = tsurg)
  )
  list(name = sprintf("S11 debridement performed on day %d", round(tsurg)),
       mod = m, ev = e, end = max(400, tsurg + 320))
}

## ============================================================================
##  Execution and summary
## ============================================================================

summarise_scn <- function(s, verbose = TRUE) {
  out <- as.data.frame(run_scn(s$mod, s$ev, end = s$end))
  txend <- as.numeric(param(s$mod)$TABX1)
  if (!is.finite(txend) || txend > s$end) txend <- s$end
  i_tx  <- which.min(abs(out$time - txend))
  i_end <- nrow(out)
  ## RSTER must not be read at "the end of treatment" — once sterilised, numerator and
  ## denominator are both meaningless. Read it at the moment that decides the outcome, i.e.
  ## early on, once the post-operative drug concentration has reached steady state (day 7 by default).
  i_d7  <- which.min(abs(out$time - min(7, txend)))

  res <- data.frame(
    scenario   = s$name,
    RSTER_d7   = round(out$RSTER[i_d7], 3),
    KILLBF_d7  = round(out$KILLBF[i_d7], 3),
    KILLPS_d7  = round(out$KILLPS[i_d7], 3),
    KIMMPS_d7  = round(out$KIMMPS[i_d7], 3),
    ETAA_d7    = round(out$ETAA[i_d7], 3),
    TSTER_d7   = round(out$TSTER[i_d7], 1),
    logB_tx    = round(out$LBTOT[i_tx], 2),
    logPersist = round(out$LBFP[i_tx], 2),
    logRes_tx  = round(out$LBRES[i_tx], 2),
    reservoir  = signif(out$RESV[i_tx], 3),
    P_relapse  = round(out$PRLPS[i_tx], 3),
    logB_end   = round(out$LBTOT[i_end], 2),
    CRP_tx     = round(out$CRP[i_tx], 1),
    ESR_tx     = round(out$ESR[i_tx], 1),
    SEQ_end    = round(out$SEQ[i_end], 2),
    PERF_end   = round(out$PERF[i_end], 3),
    boneloss   = round(out$BMLOSS[i_end], 1),
    eGFR_min   = round(min(out$EGFRC), 1),
    PLT_min    = round(min(out$PLT), 0),
    ALT_max    = round(max(out$ALT), 0),
    stringsAsFactors = FALSE
  )
  if (verbose) {
    cat("\n=== ", s$name, " ===\n", sep = "")
    print(t(res[, -1]))
  }
  invisible(list(sim = out, summary = res))
}

run_all <- function() {
  scns <- list(
    scn01(), scn02(), scn03(), scn04(), scn04b(), scn05(),
    scn06(FALSE), scn06(TRUE), scn07(), scn08(FALSE), scn08(TRUE), scn10(),
    scn11(2), scn11(120)
  )
  res <- lapply(scns, summarise_scn, verbose = FALSE)
  tab <- do.call(rbind, lapply(res, function(x) x$summary))
  cat("\n############ scenario summary ############\n")
  print(tab, row.names = FALSE)

  cat("\n############ duration sweep (S9) ############\n")
  sw <- list()
  for (w in c(2, 4, 6, 12, 26)) {
    for (rif in c(TRUE, FALSE)) {
      s <- scn09(w, with_rif = rif, deblog = 3.0, impkeep = TRUE)
      sw[[length(sw) + 1]] <- summarise_scn(s, verbose = FALSE)$summary
    }
  }
  swt <- do.call(rbind, sw)
  print(swt[, c("scenario", "RSTER_d7", "KILLPS_d7", "logB_tx", "logPersist",
                "reservoir", "P_relapse", "logB_end")], row.names = FALSE)

  invisible(list(scenarios = tab, sweep = swt, raw = res))
}

## ============================================================================
##  Structural diagnostics — claims of the model that can be checked without simulating
## ============================================================================

diag_penalties <- function() {
  cat("\n--- spatial penalty: Thiele modulus and penetration efficiency ---\n")
  Ls <- c(0.05, 0.10, 0.18, 0.25, 0.32, 0.50)
  for (drug in c("vancomycin", "levofloxacin", "daptomycin")) {
    p <- drug_A_library[[drug]]
    phi <- Ls * sqrt(p$AKLS / p$ADEFF)
    eta <- ifelse(phi < 1e-6, 1, tanh(phi) / phi)
    cat(sprintf("%-14s D_eff=%.3f k=%.1f : ", drug, p$ADEFF, p$AKLS))
    cat(paste(sprintf("L=%.2f eta=%.3f", Ls, eta), collapse = " | "), "\n")
  }
  cat("\n  reading: SEQ 2 -> 15 cm3 takes L from 0.18 -> 0.32 cm.\n")
  cat("        the eta of vancomycin halves across that range —\n")
  cat("        the reason the same drug works less well in a late patient is geometry.\n")

  cat("\n--- phenotypic penalty: kill rate by pool at the same concentration (1/d) ---\n")
  cat(sprintf("%-14s %8s %8s %8s %8s\n", "drug", "plank", "bf-rim", "persist", "gi/gd"))
  for (drug in names(drug_A_library)) {
    p <- drug_A_library[[drug]]
    h <- 0.8                       # assume a near-saturating concentration (removes the concentration problem)
    kpl <- h * (p$AEMGD * 0.35 + p$AEMGI)
    kbf <- h * (p$AEMGD * 0.12 + p$AEMGI)
    kps <- h * (p$AEMGD * 0.00 + p$AEMGI)
    cat(sprintf("%-14s %8.3f %8.3f %8.3f %8.3f\n", drug, kpl, kbf, kps,
                p$AEMGI / p$AEMGD))
  }
  p <- drug_B_library$rifampicin
  cat(sprintf("%-14s %8.3f %8.3f %8.3f %8.3f\n", "rifampicin",
              0.8 * (p$BEMGD * 0.35 + p$BEMGI),
              0.8 * (p$BEMGD * 0.12 + p$BEMGI),
              0.8 * (p$BEMGD * 0.00 + p$BEMGI),
              p$BEMGI / p$BEMGD))
  cat("\n  reading: the last column (gi/gd) is the drug ranking for this disease. Even at\n")
  cat("        saturating concentration the persister column is set by Emax_gi alone — dose-independent.\n")
}

diag_resistance <- function() {
  cat("\n--- rpoB resistance: arithmetic, not probability ---\n")
  for (B in c(1e6, 1e7, 1e8, 1e9, 1e10)) {
    cat(sprintf("  B_total=%.0e  P(mutant present) = %.4f\n", B, 1 - exp(-1e-8 * B)))
  }
  cat("  From the moment it passes 1e8, rifampicin alone is not a 'risk' but a 'booked failure'.\n")
}

diag_tster <- function() {
  cat("\n--- T_ster = ln(P0)/(KILLPS + KWAKE) : how debridement buys time ---\n")
  kwake <- 0.08
  for (kill in c(0.00, 0.01, 0.20, 1.11)) {
    cat(sprintf("  KILLPS=%.2f/d (%s)\n", kill,
        c("0.00" = "helpless against persisters", "0.01" = "vancomycin",
          "0.20" = "FQ alone", "1.11" = "with rifampicin")[sprintf("%.2f", kill)]))
    for (P0 in c(1e8, 1e5, 1e2, 1e1)) {
      cat(sprintf("      P0=%.0e -> T_ster = %6.0f d\n",
                  P0, log(P0) / (kill + kwake)))
    }
  }
  cat("  Read across for the effect of debridement (P0), down for that of rifampicin (KILLPS).\n")
  cat("  They do not substitute for each other — one acts logarithmically, the other linearly.\n")
}

## ---------------------------------------------------------------------------
##  Batch run (option-gated so that it does not fire when a Shiny app sources this file)
##
##      Rscript -e 'options(com.autorun = TRUE); source("com_mrgsolve_model.R")'
## ---------------------------------------------------------------------------
if (isTRUE(getOption("com.autorun", FALSE))) {
  diag_penalties()
  diag_resistance()
  diag_tster()
  invisible(run_all())
}
