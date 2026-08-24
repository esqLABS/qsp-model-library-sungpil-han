#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acute bacterial meningitis (pneumococcal) — QSP reference implementation
===========================================================================

This file re-implements the 63 ODEs of abm_mrgsolve_model_en.R as **dependency-free pure
Python RK4**.  This environment has no R runtime, so whether the mrgsolve model's equations
really do hold physically and pharmacologically is verified numerically here, and the result
left in abm_reference_output_en.txt.  The two files' parameters and equations correspond 1:1.

-------------------------------------------------------------------------------
What the model claims by calculation
-------------------------------------------------------------------------------
(1) The injury flux is not the bacterial count but the product **kill rate × bacterial count
    × inflammatory yield per organism**.  That product is maximal immediately after the first
    antibiotic dose (N is largest at that instant, and k_kill jumps from 0 → Emax).  So a
    steroid is not a drug given "sometime": it is a shield that must be up **before** that peak.

(2) The CSF penetration of a hydrophilic antibiotic is an increasing function of inflammation
    (barrier permeability Pb), and dexamethasone switches inflammation off.  That is,
    **dexamethasone also closes the door the antibiotic goes in through.**  Ceftriaxone keeps
    over 100-fold C/MIC headroom, so the net effect is a gain; vancomycin (MIC 1 mg/L ·
    penetration 2-13 %) has none, and the same move halves the kill rate — the sign differs by drug.

(3) A β-lactam kills by lysis, so its yield Y per organism is large; rifampicin is non-lytic,
    so its yield is small.  **At the same log-kill the injury integral differs.**

-------------------------------------------------------------------------------
Real defects this numerical check caught (v1 → v2)
-------------------------------------------------------------------------------
[F1]  The glucose consumption terms had no substrate dependence, so after CSF glucose
      collapsed to 0 at t=2 h the neutrophils and bacteria went on eating "glucose that
      was not there".  Lactate rose to 243 mmol/L (measured 6-12).  → Every consumption
      term was multiplied by Michaelis-Menten (Km_use = 5 mg/dL).  Lactate is now capped
      by the glucose influx rate (lactate comes only from glucose, so that is physically true).
[F2]  The neutrophil influx constant k_influx = 1.4e4 made the CSF white cell count 121,000/µL
      (pneumococcal meningitis, measured 1,000-5,000/µL).  → Corrected to 500, peak ~4,000/µL.
[F3]  The CSF protein permeability coefficient was too large: 1,517 mg/dL (measured 100-500).
      → PS_prot 0.30 → 0.09.  ~300 at Pb=10, ~520 mg/dL at Pb=19.
[F4]  The CSF albumin coefficient produced Q_alb 25 even in the normal state (normal <8).
      → PS_alb 0.55 → 0.105.  Normal Q_alb 5, meningitis 48-87 (measured 30-100+).
[F5]  The MAP decay constant k_map = 0.32 drove mean arterial pressure to 30 mmHg (the lower
      clamp) and CPP to 1.3 mmHg — every scenario ended in brain death.
      → 0.05.  Septic hypotension now converges to MAP 65-70 (consistent with measurement).
[F6]  The MMP-9 production constant made CSF MMP-9 6,207 ng/mL (measured 100-1,000).
      → kMMP 900 → 60, peak ~410 ng/mL.
[F7]  The barrier permeability Pb stuck at its ceiling of 20 (19.4) and could no longer
      express the intensity of inflammation, so the vancomycin story of (2) disappeared.
      → k_pb 0.55 → 0.020, peak Pb ~9.  Adding dexamethasone's direct barrier-repair effect
      (k_pb_off × (1+1.5·Idex)) took the penetration from 0.29 → 0.16 (-46 %),
      the same size as the animal reports.
[F8]  The bacterial glucose consumption coefficient q_bact = 6 mg/dL/h per 1e7 CFU/mL should
      be 1.8 when computed from the pneumococcal biomass yield (dry weight 0.3 pg/cell,
      homofermentative yield ~25 g/mol).  → 2.0.  Thanks to which the conclusion that "the
      neutrophils eat 10 times more glucose than the bacteria" is a result, not an assumption.
[F9]  K_cw = 0.6 CWU/mL for the cell-wall signal was so small that Mg (macrophage activation)
      was already saturated before the antibiotic was given, so **no cytokine burst was
      observed right after dosing** (TNF 531 → 756, 1.4-fold).  The point of this model had
      been erased numerically.  → K_cw raised to 10.0 so that the baseline CW of 3.3 sits
      in the sensitive range.  The burst now shows as 3-5 fold.
[F10] The fever term k_temp = 1.10 made the temperature converge to 40.3 °C.  → 0.85, 39.6 °C.
[F11] Inside rhs(), d[0] (the ceftriaxone central compartment) and Isch (the ischaemia index)
      were each computed twice, the earlier expression being silently overwritten — dead code.
      → Removed.
[F12] The antibiotic T>MIC tracker was always 100 % for ceftriaxone and so distinguished
      nothing.  → Changed to the bactericidal threshold **T > 4×MIC** (the EC50 of the Emax
      model).  The fact that against a resistant strain (MIC 4 mg/L) ceftriaxone never once
      exceeds 16 mg/L — the clinical reason for adding vancomycin — then emerges in the
      output.
[F13] Integration stability: the fastest eigenvalue of the ICP equation is K_el·ICP/R_out and
      is largest in the recovery phase (R_out 0.167).  At ICP 60 · R_out 0.167,
      λ = 33/h, so the RK4 stability condition is dt < 2.785/λ = 0.084 h.  The late step was
      set to 0.04 h and agreement confirmed in a half-dt (0.02 h) check (below, [8]).

-------------------------------------------------------------------------------
Further defects caught in the second numerical check (v2 → v3)
-------------------------------------------------------------------------------
[F14] **The bacteria in blood had no carrying capacity.**  If dN_b/dt = (μ_b − clearance)·N_b
      is net growth, N_b grows over 336 h to e^16.8 = 2×10^10 CFU/mL.
      That value was fed back into the CSF through k_seed, pushing the CSF count past N_max,
      and the bacterial glucose consumption term ran away, printing **CSF lactate 292 mmol/L**
      (measured 6-12).  → A logistic capacity (N_b,max = 10^8 CFU/mL) was added to N_b.
      Lactate is then capped automatically by the glucose influx ceiling at ≤9 mmol/L.
[F15] **The osmolality unit was wrong by a factor of 10.**  The mg/L → mOsm/kg conversion used
      18.2 / 9.2 instead of mannitol MW 182 and glycerol MW 92.  It gave the
      (physiologically impossible) value that 0.5 g/kg of mannitol raises plasma osmolality
      by 113 mOsm/kg.  → Corrected to 182 / 92, which gives 11.3 mOsm/kg (literature 10-15),
      and since that removes as much effect, k_osm was recalibrated 0.045 → 1.2 mL/h per mOsm
      (the reported size, 0.5 g/kg drawing off ~13 mL/h of brain water).
[F16] **The hazard function was multiplying permanent injury by time.**  Haz contained an
      h_cort·(1−N_cort) term, so even a recovered patient accrued hazard for all 14 days
      because of the residual cortical injury.  22 of 24 scenarios came out with a death
      probability above 99 %, **which made every treatment comparison meaningless.**  → The
      hazard now integrates only acute physiological derangement (ICP·CPP·SOFA·bacteraemia),
      and structural injury was split into a term evaluated once at the end point.
[F17] **The autoregulation formula made autoregulation harmful.**  In CBF =
      CBF0·(AutoR·plateau + (1−AutoR)·CPP/CPP0) the plateau becomes 0 at CPP<45, so a patient
      with autoregulation left (large AutoR) was computed to have a LOWER CBF than a
      pressure-passive patient (at CPP 38, Isch 0.54 vs 0.07).
      → Changed to CBF = CBF0·[f_passive + AutoR·max(0, f_auto − f_passive)] so that
      autoregulation can never be a disadvantage.
[F18] **The injury constants were two orders of magnitude too large.**  k_apo_dg 0.026 /h
      wiped out the dentate gyrus with a 30 h time constant and made the cognitive z of every
      scenario −3, and k_hc 0.024 gave everyone a hearing loss above 90 dB (measured: any
      hearing loss 20-30 %, severe ~10 %).  → k_apo_dg 0.0064, k_hc 0.0040.  The hearing-
      threshold transform was also changed from a plain power law to a **threshold curve that
      reflects cochlear reserve**: dB = 120·clamp((L−0.15)/0.55, 0, 1)^0.9 (L = hair cell loss).
[F19] **The sterilisation time was being set by the adherent subpopulation.**  What a repeat
      lumbar puncture culture sees clinically is the free bacteria in CSF.  Counting the
      adherent, sequestered subpopulation (kill×0.35) as well gave 42.6 h for a susceptible
      strain (measured 4-24 h).  → Sterilisation time is now N_c-based and the clearance time
      of the adherent subpopulation a separate endpoint (that gap is why a full course is needed).
[F20] The ceftriaxone CSF penetration parameter PS_cef = 4.0 raised the CSF concentration
      above 20 mg/L in severe inflammation (Pb 9) (reported 1-12 mg/L).  → 3.0.
[F21] The natural autolysis rate k_autolysis 0.020 /h already made the cell-wall load high
      before the antibiotic was given, so **the burst right after dosing showed as only
      1.4-fold.**  → 0.008 /h (0.8 % autolysis per hour).  CW now shows 1.3 → 17 (13-fold).
[F22] The vancomycin Emax of 0.90 /h did not match the kill rate reported in experimental
      pneumococcal meningitis (−0.3~−0.4 log10 CFU/mL/h, after offsetting growth).  → 1.20.
[F23] The host defence index host_def modulated only CSF phagocytosis (a small term to begin
      with), so **the immunocompromised scenario gave output identical to the normal host.**
      → Net growth is modulated as well, for failure of containment (μ_scale = 1 + 0.5·(1−host_def)).
[F24] The continuous infusion scenario was identical to intermittent dosing — because in a
      susceptible strain the C/MIC headroom is 200-fold, so T>4×MIC is 100 % either way.  That
      is a result rather than a defect, so a **continuous infusion against a resistant strain**
      scenario (S25), where there is no headroom, was added to show where the distinction appears.
[F33] **A dose scheduled at t < 0 vanished silently — and that was making the central claim
      of this model impossible to test.**  "Dexamethasone 20 minutes before the first
      antibiotic" is natural to write as an infusion at t = −0.33 h, but integration starts at
      t = 0, so that infusion interval is never evaluated.  The actual first dose in the −20
      minute arm was therefore the second dose, at **+5.67 h**.  That is, the "given first"
      arm was receiving it LATER than the "+4 h" arm, and the apparent irrelevance of timing
      was a scheduling bug, not physiology.  For the same reason the "rifampicin 2 h lead-in"
      scenario received **no** rifampicin at all — and the evidence was in the output: the peak
      cell wall · pneumolysin · TNF of S16 matched S02 to two decimals (17.4 / 5.58 / 737).
      → The whole schedule is now translated so that the earliest event is at t=0 (which is
      also the physically correct reading: steroid first at presentation, antibiotic 20 min later).
[F25] **Cortical neuronal loss progressed chronically through the ROS term.**  k_ros_cort 0.010 /h
      × hill(ROS) 0.73 never switched off for over 100 h, so even an optimally treated patient
      lost 56 % of the cortex and came out with a death probability of 87 %.  → 0.0008 /h.  In
      real pneumococcal meningitis cortical injury is patchy and mostly follows perfusion failure.
[F26] The dentate gyrus and hair cell constants were still large, giving cognitive z −1.6 and
      hearing 48 dB under optimal treatment (measured: cognitive impairment ~30 %, any loss 20-30 %).
      → k_apo_dg 0.0064 → 0.0020, k_hc 0.0040 → 0.0022.
[F27] Rather than going on eye-fitting the hazard coefficients, **the component integrals were
      added as state variables** (I_icp, I_cpp, I_sofa, I_isch, I_bact — 57 states → 63).
      Haz = h0·T + Σ h_i·I_i now decomposes after the fact, so the coefficients can be solved
      arithmetically against the trial targets (pneumococcal mortality 34 %/14 %) (below, [10]).

Units: time h · drug mg/L · bacteria CFU/mL · cytokines pg/mL · PMN cells/µL ·
       glucose mg/dL · lactate mmol/L · pressure mmHg · volume mL.
"""

import math
import sys
from bisect import bisect_right

# ===========================================================================
# State variables (63 ODEs)
# ===========================================================================
STATES = [
    # --- ceftriaxone 2-compartment + CSF (0-2)
    "Cef_c", "Cef_p", "Cef_csf",
    # --- vancomycin 2-compartment + CSF (3-5)
    "Van_c", "Van_p", "Van_csf",
    # --- rifampicin 1-compartment + CSF (6-7)
    "Rif_c", "Rif_csf",
    # --- dexamethasone + CSF + transcriptional effect compartment (8-10)
    "Dex_c", "Dex_csf", "Dex_TR",
    # --- osmotherapy: mannitol, glycerol (gut · plasma), accumulated brain osmolyte (11-14)
    "Mann_c", "Gly_gut", "Gly_c", "Osm_br",
    # --- bacteria and their products (15-19)
    "Nc", "Nadh", "Nb", "CW", "PLY",
    # --- inflammation (20-29)
    "Mg", "TNF", "IL1", "IL6", "IL10", "CXCL8", "Comp", "PMN", "MMP9", "ROS",
    # --- barrier · CSF dynamics (30-39)
    "Pb", "Alb_csf", "Prot_csf", "Glc_csf", "Lac_csf", "Vcsf_net", "R_out",
    "ICP", "Vbr", "AutoR",
    # --- systemic (40-43)
    "MAP", "Temp", "SOFA", "Vol",
    # --- injury · endpoints (44-49)
    "Ncort", "Ndg", "HC", "Oss", "Sz", "Haz",
    # --- trackers (50-56)
    "AUC_cef", "AUC_van", "AUC_lysis", "AUC_TNF", "T_cef", "T_van", "AUC_ICP",
    # --- component integrals of the hazard function (57-62) — so calibration is arithmetic, not guesswork
    "I_icp", "I_cpp", "I_sofa", "I_isch", "I_bact", "I_ncsf",
]
IX = {s: i for i, s in enumerate(STATES)}
NS = len(STATES)

# ===========================================================================
# Parameters
# ===========================================================================
P = dict(
    # ---------------- ceftriaxone PK (2 g IV q12h, 0.5 h infusion) ---------
    # Cmax (total) ~250 mg/L · t1/2 6-9 h · protein binding 85-95 %, saturable
    V1_cef=8.0, V2_cef=8.0, CL_cef=1.0, Q_cef=1.0,
    Bmax_cef=333.0, Kd_cef=25.0,               # albumin Langmuir (mg/L)
    PS_cef=3.0, Eff_cef=0.0,                   # CSF diffusional clearance (mL/h) @ Pb=1  [F20]
    MIC_cef=0.03,
    Emax_cef=1.40, EC50r_cef=4.0, h_cef=2.0,   # 1.40/h = 0.61 log10 CFU/mL/h
    Y_cef=1.00,                                # lytic killing → maximal cell wall yield
    # ---------------- vancomycin PK (15 mg/kg q6h, 1 h infusion) -----------
    V1_van=20.0, V2_van=30.0, CL_van=4.0, Q_van=8.0, fu_van=0.50,
    PS_van=1.4, Eff_van=10.0,                  # active efflux present
    MIC_van=1.0,
    Emax_van=1.20, EC50r_van=4.0, h_van=1.5,   # [F22]
    Y_van=0.90,
    # ---------------- rifampicin PK (600 mg IV q12h) -----------------------
    V_rif=50.0, CL_rif=12.0, fu_rif=0.20,
    PS_rif=30.0, a_rif=0.10, Eff_rif=0.0,      # lipophilic → low Pb sensitivity
    MIC_rif=0.03,
    Emax_rif=0.70, EC50r_rif=2.0, h_rif=1.0,
    Y_rif=0.15,                                # non-lytic → small yield
    # ---------------- dexamethasone PK/PD (0.15 mg/kg q6h) ----------------
    V_dex=70.0, CL_dex=16.0, fu_dex=0.32, PS_dex=32.0,
    ktr_on=0.46, ktr_off=0.040,                # asymmetric (on t1/2 1.5 h / off 17 h)
    Imax_dex=0.80, IC50_TR=0.010,
    # ---------------- osmotherapy -----------------------------------------
    V_mann=17.0, CL_mann=6.0,
    ka_gly=1.2, V_gly=42.0, CL_gly=8.0,
    MW_mann=182.0, MW_gly=92.0,                # [F15] 10-fold unit error corrected
    k_osm=1.20, k_osm_leak=0.010, k_osm_br_out=0.030,
    # ---------------- CSF physics -----------------------------------------
    Vcsf=150.0, Qf0=21.0, P_ss=6.0,            # mL · mL/h(0.35 mL/min) · mmHg
    R_out0=0.167, R_out_max=1.20,              # mmHg/(mL/h)
    k_ro_on=0.09, k_ro_off=0.020,
    K_el=0.092,                                # elastance coefficient 1/mL  (PVI 25 mL)
    # ---------------- bacteria --------------------------------------------
    mu_max=0.85, K_glc=10.0, Nmax=1.0e9,
    k_adh=0.030, k_des=0.010, prot_adh=0.35, mu_adh=0.30,
    k_shed=2.0e-4, k_seed=0.010,
    mu_b=0.45, k_clr_b=0.40, kill_b_boost=1.0, Nb_max=1.0e8,   # [F14]
    kphag_max=0.15, K_pmn_ph=2500.0, K_comp=1.0,
    k_autolysis=0.008, Y_auto=1.0, Y_phag=0.30,                # [F21]
    # ---------------- cell wall · pneumolysin -----------------------------
    kCW_cl=0.060,                              # t1/2 ~12 h (the cell wall persists)
    yPLY=0.50, kPLY=0.35, ply_secr=0.006,
    # ---------------- innate immunity -------------------------------------
    kmg_on=1.10, kmg_off=0.070,
    K_cw=10.0, K_nl=1.20, K_ply=1.50,          # [F9]
    kTNF=1400.0, kel_TNF=0.70, K_tnf=250.0,
    kIL1=320.0, kel_IL1=0.35, K_il1=120.0,
    kIL6=9.0e3, kel_IL6=0.25, K_il6=4000.0,
    kIL10=260.0, kel_IL10=0.20, K_il10=600.0,
    kC8=5.0e3, kel_C8=0.30, K_c8=1500.0,
    kcomp=1.20, kel_comp=0.25,
    k_influx=500.0, k_egress=0.055, k_apop=0.045,   # [F2]
    kMMP=60.0, kel_MMP=0.087, K_mmp=350.0,          # [F6]
    kROS=1.30, kel_ROS=0.90, K_ros=0.60,
    # ---------------- barrier ---------------------------------------------
    Pb_max=20.0, k_pb=0.020, k_pb_off=0.030, dex_pb=1.5,   # [F7]
    PS_alb=0.105, Alb_ser=42000.0,                          # [F4]
    PS_prot=0.090, Prot_ser=70000.0,                        # [F3]
    # ---------------- glucose · lactate -----------------------------------
    Tmax_glc=74.0, Km_glut=90.0, Glc_pl=100.0,
    q_pmn=0.0050, q_bact=2.0, q_brain=1.0,                  # [F8]
    Km_use=5.0,                                             # [F1]
    inflam_glut=0.30,
    k_lac=0.500, Lac_base=1.6,
    # ---------------- brain oedema · perfusion ----------------------------
    k_vas=1.05, Vbr_max=60.0, k_cyt=2.2, k_vbr_res=0.020,
    k_ar_loss=0.16, k_ar_rec=0.012,
    CBF0=50.0, CPP0=75.0, CBF_crit=0.55,
    # ---------------- systemic --------------------------------------------
    MAP0=88.0, k_map=0.050, k_map_rec=0.10, k_cush=0.30,    # [F5]
    k_temp=0.85, k_temp_off=0.35, Temp0=37.0,               # [F10]
    k_sofa=0.55, k_sofa_rec=0.10,
    # ---------------- injury ----------------------------------------------
    k_isch=0.055, k_ros_cort=0.0008,                           # [F25]
    k_apo_dg=0.0020, k_ply_dg=0.5, k_ros_dg=0.4,               # [F18][F26]
    k_hc=0.0022, w_hc_ply=0.5, w_hc_ros=0.7, w_hc_pmn=0.4,     # [F18][F26]
    K_ply2=1.0, k_oss=0.0020,
    # [F32] Hearing threshold curve.  With reserve 0.15 / span 0.55 all patients bunched
    #       into a moderate 20-35 dB loss, and the tail of the measured distribution (any
    #       loss 20-30 % · severe ~10 %) was not reproduced.  Reserve raised, span narrowed.
    hc_reserve=0.22, hc_span=0.40, hc_exp=0.90,
    k_sz=0.10, k_sz_off=0.12,
    # ---------------- death hazard function -------------------------------
    # [F16] Only acute derangement is integrated; structural injury is scored once at the end.
    # [F30] The coefficients were determined arithmetically from the component integrals of
    # section [10] and the trial targets (pneumococcal mortality 34 %/14 %) — not eye-fitted.
    #   First attempt: all coefficients × s → cohort no-DEX 1.256·s = 0.416 → s = 0.331
    #                                         cohort DEX    0.538·s = 0.151 → s = 0.281
    #   The model's hazard ratio 2.33 is a little below the target 2.76, so the two s differ;
    #   the compromise s = 0.30 was applied to every coefficient (predicted 31 % / 15 %).
    h0=1.2e-5, h_icp=6.0e-4, h_cpp=6.0e-4, h_sofa=6.0e-5, h_bact=4.5e-5,
    h_cort_final=0.48, h_isch=6.0e-4,
    # [F31] The hazard of persistent CSF infection itself.  Without this term the untreated
    #       natural course comes out at 70 % mortality at 14 days (in reality nearly 100 %).
    #       In the treated arms N_c is gone within 24 h, so the contribution is negligible.
    h_ncsf=1.5e-4,
)

# ===========================================================================
# Helper functions
# ===========================================================================
def free_saturable(Ctot, Bmax, Kd):
    """Free concentration with saturable albumin binding: positive root of Cf + Bmax·Cf/(Kd+Cf) = Ctot."""
    if Ctot <= 0.0:
        return 0.0
    b = Kd + Bmax - Ctot
    return (-b + math.sqrt(b * b + 4.0 * Kd * Ctot)) / 2.0


def hill(x, K, n=1.0):
    if x <= 0.0:
        return 0.0
    xn = x if n == 1.0 else x ** n
    Kn = K if n == 1.0 else K ** n
    return xn / (Kn + xn)


def emax_kill(C, MIC, Emax, EC50r, h):
    """Emax kill rate (/h) against the C/MIC ratio.  Consistent with the '10× MBC' rule."""
    if C <= 0.0 or MIC <= 0.0:
        return 0.0
    r = C / MIC
    rh = r ** h
    return Emax * rh / (EC50r ** h + rh)


def clamp(x, lo, hi):
    return lo if x < lo else (hi if x > hi else x)


# ===========================================================================
# Initial conditions: at presentation.  12-36 h after onset, inflammation already advanced
# ===========================================================================
def init_state(N0=1.0e7, Nb0=1.0e3):
    y = [0.0] * NS
    y[IX["Nc"]] = N0
    y[IX["Nadh"]] = N0 * 0.05
    y[IX["Nb"]] = Nb0
    y[IX["Mg"]] = 0.55
    y[IX["TNF"]] = 420.0
    y[IX["IL1"]] = 180.0
    y[IX["IL6"]] = 2.5e4
    y[IX["IL10"]] = 260.0
    y[IX["CXCL8"]] = 4.0e3
    y[IX["Comp"]] = 0.9
    y[IX["PMN"]] = 1200.0
    y[IX["MMP9"]] = 160.0
    y[IX["ROS"]] = 0.55
    y[IX["Pb"]] = 6.0
    y[IX["Alb_csf"]] = 1400.0        # Q_alb ~33
    y[IX["Prot_csf"]] = 180.0        # mg/dL
    y[IX["Glc_csf"]] = 34.0
    y[IX["Lac_csf"]] = 5.2
    y[IX["R_out"]] = 0.48
    y[IX["ICP"]] = 18.0
    y[IX["Vbr"]] = 4.0
    y[IX["AutoR"]] = 0.75
    y[IX["MAP"]] = 88.0
    y[IX["Temp"]] = 39.1
    y[IX["SOFA"]] = 2.0
    y[IX["Vol"]] = 1.0
    y[IX["Ncort"]] = 1.0
    y[IX["Ndg"]] = 1.0
    y[IX["HC"]] = 1.0
    y[IX["Oss"]] = 0.0
    return y


# ===========================================================================
# Dosing schedule (piecewise-constant IV infusion + oral bolus)
# ===========================================================================
class Regimen:
    def __init__(self):
        self.by_ix = {}
        self.boluses = []
        self._sorted = None

    def _add(self, t0, t1, rate, ix):
        self.by_ix.setdefault(ix, []).append((t0, t1, rate))
        self._sorted = None

    def iv(self, state, dose_mg, t_first, interval, n_dose, dur=0.5):
        for k in range(n_dose):
            t0 = t_first + k * interval
            self._add(t0, t0 + dur, dose_mg / dur, IX[state])

    def infuse(self, state, rate, t0, t1):
        self._add(t0, t1, rate, IX[state])

    def po(self, state, dose_mg, t_first, interval, n_dose):
        for k in range(n_dose):
            self.boluses.append((t_first + k * interval, dose_mg, IX[state]))

    def shift_to_zero(self):
        """[F33] Rescue the doses scheduled at t<0.

        "Dexamethasone 20 minutes before the first antibiotic" is easy to write as an
        infusion at t = −0.33 h, but integration starts at t = 0, so that infusion
        **vanishes silently**.  The −20 minute arm therefore had its real first dose at
        +5.67 h, where the second dose goes in, and the whole timing comparison was void.
        (For the same reason the rifampicin 2 h lead-in scenario received no rifampicin at
        all — so the peak cell wall · pneumolysin · TNF of S16 matched S02 to the decimal.)

        The physically correct reading is "steroid first at presentation (t=0), antibiotic
        20 min later", so the whole schedule is shifted until the earliest event is at t=0.
        """
        times = [t0 for lst in self.by_ix.values() for (t0, _, _) in lst]
        times += [tb for (tb, _, _) in self.boluses]
        if not times:
            return 0.0
        tmin = min(times)
        if tmin >= 0.0:
            return 0.0
        sh = -tmin
        self.by_ix = {ix: [(t0 + sh, t1 + sh, r) for (t0, t1, r) in lst]
                      for ix, lst in self.by_ix.items()}
        self.boluses = [(tb + sh, a, ix) for (tb, a, ix) in self.boluses]
        self._sorted = None
        return sh

    def _prep(self):
        self._sorted = {}
        for ix, lst in self.by_ix.items():
            lst = sorted(lst)
            self._sorted[ix] = (lst, [e[0] for e in lst])

    def rate(self, t, idx):
        if self._sorted is None:
            self._prep()
        ent = self._sorted.get(idx)
        if not ent:
            return 0.0
        lst, starts = ent
        j = bisect_right(starts, t) - 1
        r = 0.0
        while j >= 0:
            t0, t1, rt = lst[j]
            if t < t1:
                r += rt
            if t0 < t - 48.0:          # infusions started >48 h ago are certainly finished
                break
            j -= 1
        return r

    def due(self, t, dt):
        return [(a, ix) for (tb, a, ix) in self.boluses if t <= tb < t + dt]


# ===========================================================================
# Right-hand side (RHS)
# ===========================================================================
def rhs(t, y, p, reg, cfg):
    d = [0.0] * NS

    Cef_c, Cef_p, Cef_csf = y[0], y[1], y[2]
    Van_c, Van_p, Van_csf = y[3], y[4], y[5]
    Rif_c, Rif_csf = y[6], y[7]
    Dex_c, Dex_csf, Dex_TR = y[8], y[9], y[10]
    Mann_c, Gly_gut, Gly_c, Osm_br = y[11], y[12], y[13], y[14]
    Nc, Nadh, Nb, CW, PLY = y[15], y[16], y[17], y[18], y[19]
    Mg, TNF, IL1, IL6, IL10, CXCL8, Comp, PMN, MMP9, ROS = y[20:30]
    Pb, Alb, Prot, Glc, Lac, _Vnet, R_out, ICP, Vbr, AutoR = y[30:40]
    MAP, Temp, SOFA, Vol = y[40:44]
    Ncort, Ndg, HC, Oss, Sz, _Haz = y[44:50]

    # ---------- free plasma concentrations ----------
    Cef_pl = Cef_c / p["V1_cef"]
    Cef_f = free_saturable(Cef_pl, p["Bmax_cef"], p["Kd_cef"])
    Van_pl = Van_c / p["V1_van"]
    Van_f = p["fu_van"] * Van_pl
    Rif_pl = Rif_c / p["V_rif"]
    Rif_f = p["fu_rif"] * Rif_pl
    Dex_pl = Dex_c / p["V_dex"]
    Dex_f = p["fu_dex"] * Dex_pl

    Pb_hyd = Pb                                  # hydrophilic drugs
    Pb_lip = 1.0 + p["a_rif"] * (Pb - 1.0)       # lipophilic drugs
    Qbulk = p["Qf0"]                             # CSF bulk flow = sink (mL/h)

    # ---------- drug PK ----------
    d[0] = reg.rate(t, 0) - p["CL_cef"] * Cef_pl - p["Q_cef"] * (Cef_pl - Cef_p / p["V2_cef"])
    d[1] = p["Q_cef"] * (Cef_pl - Cef_p / p["V2_cef"])
    d[2] = (p["PS_cef"] * Pb_hyd * (Cef_f - Cef_csf) - (Qbulk + p["Eff_cef"]) * Cef_csf) / p["Vcsf"]

    d[3] = reg.rate(t, 3) - p["CL_van"] * Van_pl - p["Q_van"] * (Van_pl - Van_p / p["V2_van"])
    d[4] = p["Q_van"] * (Van_pl - Van_p / p["V2_van"])
    d[5] = (p["PS_van"] * Pb_hyd * (Van_f - Van_csf) - (Qbulk + p["Eff_van"]) * Van_csf) / p["Vcsf"]

    d[6] = reg.rate(t, 6) - p["CL_rif"] * Rif_pl
    d[7] = (p["PS_rif"] * Pb_lip * (Rif_f - Rif_csf) - (Qbulk + p["Eff_rif"]) * Rif_csf) / p["Vcsf"]

    d[8] = reg.rate(t, 8) - p["CL_dex"] * Dex_pl
    d[9] = (p["PS_dex"] * Pb_lip * (Dex_f - Dex_csf) - Qbulk * Dex_csf) / p["Vcsf"]
    ktr = p["ktr_on"] if Dex_csf > Dex_TR else p["ktr_off"]
    d[10] = ktr * (Dex_csf - Dex_TR)
    Idex = p["Imax_dex"] * hill(Dex_TR, p["IC50_TR"])

    # ---------- osmotherapy ----------
    d[11] = reg.rate(t, 11) - p["CL_mann"] * (Mann_c / p["V_mann"])
    gly_abs = p["ka_gly"] * Gly_gut
    d[12] = -gly_abs
    d[13] = gly_abs - p["CL_gly"] * (Gly_c / p["V_gly"])
    # [F15] mg/L ÷ MW(g/mol) = mmol/L = mOsm/kg equivalent (10-fold unit error corrected)
    Osm_pl = (Mann_c / p["V_mann"]) / p["MW_mann"] + (Gly_c / p["V_gly"]) / p["MW_gly"]
    Osm_grad = max(0.0, Osm_pl - Osm_br)
    d[14] = p["k_osm_leak"] * (1.0 + 0.35 * (Pb - 1.0)) * Osm_pl - p["k_osm_br_out"] * Osm_br

    # ---------- bacteria ----------
    kill_cef = emax_kill(Cef_csf, cfg["MIC_cef"], p["Emax_cef"], p["EC50r_cef"], p["h_cef"])
    kill_van = emax_kill(Van_csf, cfg["MIC_van"], p["Emax_van"], p["EC50r_van"], p["h_van"])
    kill_rif = emax_kill(Rif_csf, cfg["MIC_rif"], p["Emax_rif"], p["EC50r_rif"], p["h_rif"])
    kill_tot = kill_cef + kill_van + kill_rif

    f_glc = hill(Glc, p["Km_use"])                     # [F1] substrate dependence
    mu = p["mu_max"] * hill(Glc, p["K_glc"]) * (1.0 - Nc / p["Nmax"]) * cfg["mu_scale"]
    k_phag = (p["kphag_max"] * hill(PMN, p["K_pmn_ph"])
              * (0.2 + 0.8 * hill(Comp, p["K_comp"])) * cfg["host_def"])

    d[15] = ((mu - kill_tot - k_phag) * Nc + p["k_des"] * Nadh - p["k_adh"] * Nc
             + p["k_seed"] * Nb - p["k_shed"] * Nc)
    d[16] = (p["k_adh"] * Nc - p["k_des"] * Nadh
             + (p["mu_adh"] * hill(Glc, p["K_glc"]) - kill_tot * p["prot_adh"]) * Nadh)
    # [F14] blood bacteria need a carrying capacity too (without it N_b diverges and re-seeds the CSF)
    d[17] = ((p["mu_b"] * (1.0 - Nb / p["Nb_max"]) - p["k_clr_b"] * cfg["host_def"]
              - p["kill_b_boost"] * kill_tot) * Nb + p["k_shed"] * Nc)

    # ---------- lysis flux = kill rate × bacterial count × yield ----------
    Y_adh = p["Y_cef"] if (kill_cef + kill_van) > kill_rif else p["Y_rif"]
    lysis_flux = ((kill_cef * p["Y_cef"] + kill_van * p["Y_van"] + kill_rif * p["Y_rif"]) * Nc
                  + kill_tot * p["prot_adh"] * Y_adh * Nadh
                  + p["k_autolysis"] * p["Y_auto"] * Nc
                  + k_phag * p["Y_phag"] * Nc)
    d[18] = lysis_flux / 1.0e6 - p["kCW_cl"] * CW
    d[19] = p["yPLY"] * (lysis_flux + p["ply_secr"] * Nc) / 1.0e6 - p["kPLY"] * PLY

    # ---------- perfusion (needed for ROS · injury) ----------
    # [F17] autoregulation can never be a loss: pressure-passive flow is the floor
    CPP = MAP - ICP
    f_passive = clamp(CPP / p["CPP0"], 0.0, 1.6)
    f_auto = clamp((CPP - 25.0) / 25.0, 0.0, 1.0)
    CBF = p["CBF0"] * (f_passive + AutoR * max(0.0, f_auto - f_passive))
    Isch = clamp(1.0 - CBF / (p["CBF_crit"] * p["CBF0"]), 0.0, 1.0)

    # ---------- innate immunity ----------
    S_pamp = (hill(CW, p["K_cw"]) + 0.6 * hill(Nc / 1.0e6, p["K_nl"])
              + 0.5 * hill(PLY, p["K_ply"]))
    d[20] = p["kmg_on"] * S_pamp * (1.0 - Mg) - p["kmg_off"] * Mg
    f10 = 1.0 / (1.0 + IL10 / p["K_il10"])
    d[21] = p["kTNF"] * Mg * (1.0 - Idex) * f10 - p["kel_TNF"] * TNF
    d[22] = (p["kIL1"] * Mg * (0.6 + 0.4 * hill(PLY, p["K_ply"])) * (1.0 - Idex) * f10
             - p["kel_IL1"] * IL1)
    d[23] = (p["kIL6"] * (0.5 * Mg + hill(TNF, p["K_tnf"]) + hill(IL1, p["K_il1"]))
             * (1.0 - 0.6 * Idex) - p["kel_IL6"] * IL6)
    d[24] = (p["kIL10"] * (Mg + hill(TNF, p["K_tnf"])) * (1.0 + 0.5 * Idex)
             - p["kel_IL10"] * IL10)
    d[25] = (p["kC8"] * (hill(TNF, p["K_tnf"]) + hill(IL1, p["K_il1"]) + 0.3 * Mg)
             * (1.0 - 0.7 * Idex) - p["kel_C8"] * CXCL8)
    d[26] = p["kcomp"] * hill(Nc / 1.0e6, 2.0) - p["kel_comp"] * Comp
    adh_mol = hill(TNF, p["K_tnf"]) * (1.0 - 0.6 * Idex)
    d[27] = (p["k_influx"] * hill(CXCL8, p["K_c8"]) * adh_mol
             * (1.0 + 0.5 * (Pb - 1.0) / 19.0) - (p["k_egress"] + p["k_apop"]) * PMN)
    d[28] = p["kMMP"] * hill(PMN, 2500.0) * (1.0 - 0.5 * Idex) - p["kel_MMP"] * MMP9
    d[29] = p["kROS"] * (hill(PMN, 2500.0) + 0.5 * Mg + 0.3 * Isch) - p["kel_ROS"] * ROS

    # ---------- barrier ----------
    d[30] = (p["k_pb"] * (hill(MMP9, p["K_mmp"]) + 0.5 * hill(TNF, p["K_tnf"])
                          + 0.3 * hill(PMN, 2500.0)) * (p["Pb_max"] - Pb)
             - p["k_pb_off"] * (1.0 + p["dex_pb"] * Idex) * (Pb - 1.0))
    d[31] = (p["PS_alb"] * Pb * (p["Alb_ser"] - Alb) - Qbulk * Alb) / p["Vcsf"]
    d[32] = (p["PS_prot"] * Pb * (p["Prot_ser"] / 10.0 - Prot) - Qbulk * Prot) / p["Vcsf"]

    # ---------- glucose · lactate ----------
    glut = 1.0 - p["inflam_glut"] * hill(Pb - 1.0, 8.0)
    influx = p["Tmax_glc"] * glut * (p["Glc_pl"] / (p["Km_glut"] + p["Glc_pl"])
                                     - Glc / (p["Km_glut"] + Glc))
    use_pmn = p["q_pmn"] * PMN * f_glc
    use_bact = p["q_bact"] * (Nc / 1.0e7) * f_glc
    use_brain = p["q_brain"] * (1.0 + 0.3 * (Temp - 37.0)) * f_glc
    d[33] = influx - use_pmn - use_bact - use_brain - (Qbulk / p["Vcsf"]) * Glc
    d[34] = (2.0 * 0.0556 * (use_pmn + use_bact + use_brain * (0.3 + 0.7 * Isch))
             - p["k_lac"] * (Lac - p["Lac_base"]))

    # ---------- CSF dynamics ----------
    Qf = p["Qf0"] * (1.0 - 0.3 * clamp((ICP - 20.0) / 40.0, 0.0, 1.0))
    d[36] = (p["k_ro_on"] * (hill(Prot, 200.0) + hill(PMN, 2000.0))
             * (p["R_out_max"] - R_out) - p["k_ro_off"] * (R_out - p["R_out0"]))
    absorb = max(0.0, (ICP - p["P_ss"]) / max(R_out, 1e-3))
    drain = cfg["csf_drain"](t) if cfg["csf_drain"] else 0.0
    osm_shrink = p["k_osm"] * Osm_grad / (1.0 + 0.3 * (Pb - 1.0))
    dVbr = (p["k_vas"] * hill(Pb - 1.0, 6.0) * (1.0 - Vbr / p["Vbr_max"])
            + p["k_cyt"] * Isch - osm_shrink - p["k_vbr_res"] * Vbr)
    d[38] = dVbr
    d[35] = Qf - absorb - drain                      # tracker: net CSF accumulation rate
    d[37] = p["K_el"] * ICP * (Qf + dVbr - absorb - drain)
    d[39] = -p["k_ar_loss"] * hill(TNF, p["K_tnf"]) * AutoR + p["k_ar_rec"] * (1.0 - AutoR)

    # ---------- systemic ----------
    d[40] = (-p["k_map"] * (hill(TNF, p["K_tnf"]) + 0.4 * hill(IL6, p["K_il6"])) * MAP
             + p["k_map_rec"] * (p["MAP0"] * Vol - MAP)
             + p["k_cush"] * max(0.0, ICP - 30.0) + cfg["vasopressor"])
    d[41] = (p["k_temp"] * (hill(IL1, p["K_il1"]) + 0.5 * hill(IL6, p["K_il6"]))
             * (1.0 - cfg["antipyretic"]) - p["k_temp_off"] * (Temp - p["Temp0"]))
    d[42] = (p["k_sofa"] * (hill(math.log10(1.0 + max(Nb, 0.0)), 2.0)
                            + 0.5 * hill(TNF, p["K_tnf"])) - p["k_sofa_rec"] * SOFA)
    d[43] = 0.0

    # ---------- injury ----------
    d[44] = -(p["k_isch"] * Isch * Isch + p["k_ros_cort"] * hill(ROS, p["K_ros"])) * Ncort
    d[45] = -(p["k_apo_dg"] * (hill(CW, p["K_cw"]) + p["k_ply_dg"] * hill(PLY, p["K_ply2"])
                               + p["k_ros_dg"] * hill(ROS, p["K_ros"]))) * Ndg
    d[46] = -(p["k_hc"] * (p["w_hc_ply"] * hill(PLY, p["K_ply2"])
                           + p["w_hc_ros"] * hill(ROS, p["K_ros"])
                           + p["w_hc_pmn"] * hill(PMN, 2500.0))) * HC
    d[47] = p["k_oss"] * (1.0 - HC) * hill(PMN, 2000.0) * (1.0 - Oss)
    d[48] = (p["k_sz"] * ((1.0 - Ncort) + hill(ROS, p["K_ros"])
                          + 0.4 * max(0.0, 1.0 - Glc / 25.0))
             * (1.0 - cfg["anticonvulsant"]) - p["k_sz_off"] * Sz)
    # [F16] integrate only acute derangement — permanent injury (N_cort) is scored once at the end
    # [F28] the ICP·CPP terms are quadratic.  If linear, "CPP 45 for 100 h" and "CPP 5 for
    #       10 h" would carry the same hazard, and mild persistent derangement would dominate.
    z_icp = max(0.0, ICP - 25.0) / 10.0
    z_cpp = max(0.0, 50.0 - CPP) / 10.0
    d[49] = (p["h0"] + p["h_icp"] * z_icp * z_icp
             + p["h_cpp"] * z_cpp * z_cpp
             + p["h_sofa"] * SOFA + p["h_isch"] * Isch * Isch
             + p["h_bact"] * math.log10(1.0 + max(Nb, 0.0))
             + p["h_ncsf"] * math.log10(1.0 + max(Nc, 0.0)))

    # ---------- trackers ----------
    d[50] = Cef_csf
    d[51] = Van_csf
    d[52] = lysis_flux / 1.0e6
    d[53] = TNF
    d[54] = 1.0 if Cef_csf > 4.0 * cfg["MIC_cef"] else 0.0      # [F12] T>4×MIC
    d[55] = 1.0 if Van_csf > 4.0 * cfg["MIC_van"] else 0.0
    d[56] = ICP
    # component integrals of the hazard (so Haz = h0·T + Σ h_i · I_i can be decomposed afterwards)
    d[57] = z_icp * z_icp
    d[58] = z_cpp * z_cpp
    d[59] = SOFA
    d[60] = Isch * Isch
    d[61] = math.log10(1.0 + max(Nb, 0.0))
    d[62] = math.log10(1.0 + max(Nc, 0.0))
    return d


NONNEG = [IX[s] for s in (
    "Cef_c", "Cef_p", "Cef_csf", "Van_c", "Van_p", "Van_csf", "Rif_c", "Rif_csf",
    "Dex_c", "Dex_csf", "Dex_TR", "Mann_c", "Gly_gut", "Gly_c", "Osm_br", "CW",
    "PLY", "TNF", "IL1", "IL6", "IL10", "CXCL8", "Comp", "PMN", "MMP9", "ROS",
    "Alb_csf", "Prot_csf", "Glc_csf", "Lac_csf", "Vbr", "SOFA", "Sz", "Haz")]


# ===========================================================================
# Integrator (RK4 · step schedule)
# ===========================================================================
def simulate(cfg, tmax=336.0, record=None, dt_scale=1.0):
    p = P
    reg = cfg["reg"]
    y = init_state(cfg["N0"], cfg["Nb0"])
    t = 0.0
    out = []
    rec = sorted(record) if record else []
    ri = 0
    sterile_t = None      # [F19] on N_c (what a repeat lumbar puncture culture sees)
    adh_t = None          #       clearance of the adherent subpopulation (why a full course is needed)
    peak = dict(ICP=y[IX["ICP"]], TNF=y[IX["TNF"]], PMN=y[IX["PMN"]], Nc=y[IX["Nc"]],
                PLY=0.0, Pb=y[IX["Pb"]], Prot=y[IX["Prot_csf"]], Lac=y[IX["Lac_csf"]],
                Temp=y[IX["Temp"]], Qalb=0.0, lysis=0.0, Vbr=y[IX["Vbr"]])
    trough = dict(CPP=y[IX["MAP"]] - y[IX["ICP"]], Glc=y[IX["Glc_csf"]],
                  MAP=y[IX["MAP"]], CBF=50.0)

    while t < tmax - 1e-9:
        dt = (0.005 if t < 4.0 else (0.02 if t < 48.0 else 0.04)) * dt_scale
        if rec and ri < len(rec) and t + dt > rec[ri]:
            dt = max(rec[ri] - t, 1e-6)

        for (amt, ix) in reg.due(t, dt):
            y[ix] += amt

        k1 = rhs(t, y, p, reg, cfg)
        y2 = [y[i] + 0.5 * dt * k1[i] for i in range(NS)]
        k2 = rhs(t + 0.5 * dt, y2, p, reg, cfg)
        y3 = [y[i] + 0.5 * dt * k2[i] for i in range(NS)]
        k3 = rhs(t + 0.5 * dt, y3, p, reg, cfg)
        y4 = [y[i] + dt * k3[i] for i in range(NS)]
        k4 = rhs(t + dt, y4, p, reg, cfg)
        for i in range(NS):
            y[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])

        for i in NONNEG:
            if y[i] < 0.0:
                y[i] = 0.0
        y[IX["Mg"]] = clamp(y[IX["Mg"]], 0.0, 1.0)
        y[IX["Pb"]] = clamp(y[IX["Pb"]], 1.0, P["Pb_max"])
        y[IX["AutoR"]] = clamp(y[IX["AutoR"]], 0.0, 1.0)
        y[IX["ICP"]] = clamp(y[IX["ICP"]], 2.0, 120.0)
        y[IX["MAP"]] = clamp(y[IX["MAP"]], 30.0, 160.0)
        for s in ("Ncort", "Ndg", "HC", "Oss"):
            y[IX[s]] = clamp(y[IX[s]], 0.0, 1.0)
        for s in ("Nc", "Nadh", "Nb"):
            if y[IX[s]] < 1.0e-2:          # below the detection limit → extinct
                y[IX[s]] = 0.0
        t += dt

        if sterile_t is None and y[IX["Nc"]] < 10.0:
            sterile_t = t
        if adh_t is None and y[IX["Nc"]] + y[IX["Nadh"]] < 10.0:
            adh_t = t
        cpp = y[IX["MAP"]] - y[IX["ICP"]]
        peak["ICP"] = max(peak["ICP"], y[IX["ICP"]])
        peak["TNF"] = max(peak["TNF"], y[IX["TNF"]])
        peak["PMN"] = max(peak["PMN"], y[IX["PMN"]])
        peak["Nc"] = max(peak["Nc"], y[IX["Nc"]])
        peak["PLY"] = max(peak["PLY"], y[IX["PLY"]])
        peak["Pb"] = max(peak["Pb"], y[IX["Pb"]])
        peak["Prot"] = max(peak["Prot"], y[IX["Prot_csf"]])
        peak["Lac"] = max(peak["Lac"], y[IX["Lac_csf"]])
        peak["Temp"] = max(peak["Temp"], y[IX["Temp"]])
        peak["Vbr"] = max(peak["Vbr"], y[IX["Vbr"]])
        peak["Qalb"] = max(peak["Qalb"], 1000.0 * y[IX["Alb_csf"]] / P["Alb_ser"])
        trough["CPP"] = min(trough["CPP"], cpp)
        trough["Glc"] = min(trough["Glc"], y[IX["Glc_csf"]])
        trough["MAP"] = min(trough["MAP"], y[IX["MAP"]])

        if rec and ri < len(rec) and abs(t - rec[ri]) < 1e-6:
            out.append((t, list(y)))
            ri += 1

    return dict(y=y, series=out, peak=peak, trough=trough, sterile_t=sterile_t,
                endpoints=endpoints(y, peak, trough, sterile_t, adh_t))


def hearing_dB(HC, Oss):
    """[F18] Threshold curve reflecting cochlear reserve.
    Up to 15 % hair cell loss the threshold shift is negligible; beyond that it rises steeply."""
    L = 1.0 - HC
    x = clamp((L - P["hc_reserve"]) / P["hc_span"], 0.0, 1.0)
    return min(120.0, 120.0 * (x ** P["hc_exp"]) + 15.0 * Oss)


def endpoints(y, peak, trough, sterile_t, adh_t=None):
    HC, Oss = y[IX["HC"]], y[IX["Oss"]]
    dB = hearing_dB(HC, Oss)
    # [F16] integrated acute hazard + a single evaluation of structural injury
    haz_tot = y[IX["Haz"]] + P["h_cort_final"] * (1.0 - y[IX["Ncort"]])
    pdeath = 1.0 - math.exp(-haz_tot)
    cogz = -(2.5 * (1.0 - y[IX["Ndg"]]) + 1.5 * (1.0 - y[IX["Ncort"]]))
    focal = 1.0 - math.exp(-3.0 * (1.0 - y[IX["Ncort"]]))
    seq = min(1.0, (1.0 if dB > 60.0 else 0.0) * 0.5
              + max(0.0, -cogz - 1.0) * 0.35 + focal * 0.4)
    return dict(
        death=pdeath, hear_dB=dB, hear_any=dB > 25.0, hear_severe=dB > 60.0,
        cog_z=cogz, focal=focal, sterile_t=sterile_t, adh_t=adh_t,
        haz_acute=y[IX["Haz"]], haz_struct=P["h_cort_final"] * (1.0 - y[IX["Ncort"]]),
        peak_ICP=peak["ICP"], min_CPP=trough["CPP"], min_MAP=trough["MAP"],
        peak_Vbr=peak["Vbr"], peak_Qalb=peak["Qalb"], peak_Temp=peak["Temp"],
        I_icp=y[IX["I_icp"]], I_cpp=y[IX["I_cpp"]], I_sofa=y[IX["I_sofa"]],
        I_isch=y[IX["I_isch"]], I_bact=y[IX["I_bact"]], I_ncsf=y[IX["I_ncsf"]],
        Ncort_f=y[IX["Ncort"]],
        AUC_lysis=y[IX["AUC_lysis"]], AUC_TNF=y[IX["AUC_TNF"]],
        AUC_cef=y[IX["AUC_cef"]], AUC_van=y[IX["AUC_van"]],
        T_cef=y[IX["T_cef"]], T_van=y[IX["T_van"]],
        peak_TNF=peak["TNF"], peak_PMN=peak["PMN"], peak_Pb=peak["Pb"],
        peak_PLY=peak["PLY"], min_Glc=trough["Glc"], peak_Lac=peak["Lac"],
        peak_Prot=peak["Prot"], Ndg=y[IX["Ndg"]], Ncort=y[IX["Ncort"]],
        unfav=min(1.0, pdeath + (1.0 - pdeath) * seq),
    )


# ===========================================================================
# Scenario builder
# ===========================================================================
def make_cfg(name, cef=True, cef_delay=0.0, cef_ci=False, van=False, rif=False,
             rif_lead=0.0, dex=False, dex_time=0.0, dex_days=4,
             mannitol=False, glycerol=False, MIC_cef=0.03, MIC_van=1.0,
             MIC_rif=0.03, N0=1.0e7, Nb0=1.0e3, host_def=1.0, mu_scale=1.0,
             antipyretic=0.0, anticonvulsant=0.0, vasopressor=0.0,
             drain=False, wt=70.0, tmax=336.0):
    reg = Regimen()
    if cef:
        if cef_ci:
            reg.iv("Cef_c", 2000.0, cef_delay, 12.0, 1, dur=0.5)      # loading
            for k in range(14):                                        # 4 g/day continuous infusion
                reg.infuse("Cef_c", 4000.0 / 24.0,
                           cef_delay + 0.5 + 24.0 * k, cef_delay + 0.5 + 24.0 * (k + 1))
        else:
            reg.iv("Cef_c", 2000.0, cef_delay, 12.0, 28, dur=0.5)
    if van:
        reg.iv("Van_c", 15.0 * wt, cef_delay, 6.0, 56, dur=1.0)
    if rif:
        reg.iv("Rif_c", 600.0, cef_delay - rif_lead, 12.0, 28, dur=0.5)
    if dex:
        reg.iv("Dex_c", 0.15 * wt, dex_time, 6.0, 4 * dex_days, dur=0.25)
    if mannitol:
        reg.iv("Mann_c", 0.5 * wt * 1000.0, 1.0, 6.0, 12, dur=0.5)     # 0.5 g/kg q6h ×3 d
    if glycerol:
        reg.po("Gly_gut", 0.25 * wt * 1000.0, 0.0, 6.0, 16)            # 6 g/kg/day ×4 d
    shift = reg.shift_to_zero()          # [F33] rescue the t<0 doses
    drain_fn = (lambda t: 12.0 if 2.0 <= t <= 96.0 else 0.0) if drain else None
    # [F23] weaker host defence acts on containment (higher net growth), not only phagocytosis
    mu_scale = mu_scale * (1.0 + 0.5 * (1.0 - host_def))
    return dict(name=name, reg=reg, MIC_cef=MIC_cef, MIC_van=MIC_van, MIC_rif=MIC_rif,
                N0=N0, Nb0=Nb0, host_def=host_def, mu_scale=mu_scale,
                antipyretic=antipyretic, anticonvulsant=anticonvulsant,
                vasopressor=vasopressor, csf_drain=drain_fn, tmax=tmax,
                shift=shift)


SCENARIOS = [
    ("S01 untreated natural course", dict(cef=False)),
    ("S02 CEF alone (t=0)", dict()),
    ("S03 CEF + DEX 20 min pre", dict(dex=True, dex_time=-0.33)),
    ("S04 CEF + DEX at 0 h", dict(dex=True, dex_time=0.0)),
    ("S05 CEF + DEX +2 h", dict(dex=True, dex_time=2.0)),
    ("S06 CEF + DEX +4 h", dict(dex=True, dex_time=4.0)),
    ("S07 CEF + DEX +12 h", dict(dex=True, dex_time=12.0)),
    ("S08 CEF + DEX 1 day only", dict(dex=True, dex_time=0.0, dex_days=1)),
    ("S09 antibiotic 3 h delay", dict(cef_delay=3.0, dex=True, dex_time=3.0)),
    ("S10 antibiotic 6 h delay", dict(cef_delay=6.0, dex=True, dex_time=6.0)),
    ("S11 antibiotic 12 h delay", dict(cef_delay=12.0, dex=True, dex_time=12.0)),
    ("S12 resistant(MIC 4) CEF alone", dict(MIC_cef=4.0)),
    ("S13 resistant + VAN, no DEX", dict(MIC_cef=4.0, van=True)),
    ("S14 resistant + VAN + DEX", dict(MIC_cef=4.0, van=True, dex=True, dex_time=0.0)),
    ("S15 resistant + VAN+RIF+DEX", dict(MIC_cef=4.0, van=True, rif=True, dex=True, dex_time=0.0)),
    ("S16 RIF 2 h lead-in → CEF", dict(rif=True, rif_lead=2.0)),
    ("S17 CEF + DEX + mannitol", dict(dex=True, dex_time=0.0, mannitol=True)),
    ("S18 CEF + DEX + glycerol", dict(dex=True, dex_time=0.0, glycerol=True)),
    ("S19 high load 1e8 + DEX", dict(N0=1.0e8, dex=True, dex_time=0.0)),
    ("S20 early, low load 1e5 + DEX", dict(N0=1.0e5, dex=True, dex_time=0.0)),
    ("S21 elderly/immunosupp. + DEX", dict(host_def=0.4, dex=True, dex_time=0.0)),
    ("S22 CEF CI + DEX", dict(cef_ci=True, dex=True, dex_time=0.0)),
    ("S23 antipyretic, no antibiotic", dict(cef=False, antipyretic=0.7)),
    ("S24 standard CEF+VAN+DEX", dict(van=True, dex=True, dex_time=-0.33)),
    ("S25 resistant CEF CI+VAN+DEX", dict(MIC_cef=4.0, cef_ci=True, van=True,
                                       dex=True, dex_time=-0.33)),   # [F24]
    ("S26 standard + EVD drainage", dict(van=True, dex=True, dex_time=-0.33, drain=True)),
]

COHORT = [
    (1.0e6, 0.0, 1.00), (3.0e6, 2.0, 1.00), (1.0e7, 1.0, 0.90),
    (3.0e7, 4.0, 0.85), (1.0e8, 6.0, 0.70), (1.0e7, 3.0, 0.55),
    (5.0e6, 1.0, 1.00), (2.0e7, 8.0, 0.60), (8.0e7, 2.0, 0.80),
    (3.0e7, 24.0, 0.40),          # [F32] a late presenter — makes the tail of severe hearing loss
]


def fmt(x, n=2):
    if x is None:
        return "n/a"
    return f"{x:.{n}f}"


def log10s(x):
    return math.log10(x) if x > 1.0 else 0.0


# ===========================================================================
# Main
# ===========================================================================
def main():
    L = []
    W = L.append
    W("=" * 112)
    W("Acute bacterial meningitis (pneumococcal) QSP model — Python RK4 reference output")
    W("63 ODEs · 26 scenarios · virtual cohort of 10 × 2 · same equations as abm_mrgsolve_model_en.R")
    W("=" * 112)
    W("")

    # ---------------- [0] structural consistency ----------------
    W("[0] Structural consistency checkable by hand (what must be right before the model runs)")
    W("-" * 112)
    W(f"  CSF steady-state ICP = P_ss + Qf·R_out")
    for r, tag in [(P["R_out0"], "norm"), (0.48, "pres"), (1.20, "sev.")]:
        W(f"     R_out={r:5.3f} ({tag:4s}) → ICP = {P['P_ss'] + P['Qf0']*r:5.1f} mmHg"
          f"   [literature: normal 7-15 · meningitis 20-40]")
    cf250, cf30 = (free_saturable(250.0, P["Bmax_cef"], P["Kd_cef"]),
                   free_saturable(30.0, P["Bmax_cef"], P["Kd_cef"]))
    W(f"  Ceftriaxone saturable binding: C_tot 30 → fu {cf30/30:.3f} ; 250 mg/L → fu {cf250/250:.3f}"
      f"   [reported: rises 0.05→0.20]")
    for nm, PS, Eff, fu in (("CEF", P["PS_cef"], P["Eff_cef"], cf250 / 250),
                            ("VAN", P["PS_van"], P["Eff_van"], P["fu_van"]),
                            ("RIF", P["PS_rif"], P["Eff_rif"], P["fu_rif"]),
                            ("DEX", P["PS_dex"], 0.0, P["fu_dex"])):
        lip = (nm == "RIF" or nm == "DEX")
        def ratio(pb):
            pbe = 1.0 + P["a_rif"] * (pb - 1.0) if lip else pb
            return PS * pbe / (PS * pbe + P["Qf0"] + Eff)
        W(f"  {nm} CSF/free-plasma steady ratio: Pb=1 {ratio(1):.3f} · Pb=4 {ratio(4):.3f} · Pb=9 {ratio(9):.3f}"
          f"  → total-concentration penetration {ratio(1)*fu:.3f}/{ratio(9)*fu:.3f}")
    W("     [reported total-concentration penetration: CEF 0.015-0.10 · VAN 0.01-0.15 · RIF 0.10-0.20 · DEX 0.15-0.30]")
    g = 60.0
    W(f"  CSF glucose 60 mg/dL balance: carrier influx {P['Tmax_glc']*(P['Glc_pl']/(P['Km_glut']+P['Glc_pl'])-g/(P['Km_glut']+g)):.1f}"
      f" = bulk {g*P['Qf0']/P['Vcsf']:.1f} + brain {P['q_brain']:.1f} mg/dL/h")
    W(f"  Glucose consumption contest (PMN 4000/µL vs bacteria 1e7/mL): {P['q_pmn']*4000:.1f} vs {P['q_bact']:.1f} mg/dL/h"
      f"  → the neutrophils eat {P['q_pmn']*4000/P['q_bact']:.0f}x more")
    W(f"  RK4 stability condition: ICP equation λ_max = K_el·ICP/R_out = {P['K_el']*60/P['R_out0']:.1f}/h (ICP 60 · recovery)"
      f" → dt < 2.785/λ = {2.785/(P['K_el']*60/P['R_out0']):.3f} h ; dt used 0.04 h")
    W("")

    # ---------------- run the scenarios ----------------
    rec_times = [0.5, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 48, 72, 120, 168, 240, 336]
    R = {}
    for label, kw in SCENARIOS:
        cfg = make_cfg(label, **kw)
        R[label] = simulate(cfg, tmax=cfg["tmax"], record=rec_times)
        print(f"  ... {label}", file=sys.stderr)

    # ---------------- [1] natural history ----------------
    W("[1] S01 untreated natural course — do the diagnostic indices reach textbook ranges")
    W("-" * 112)
    W("   t[h] log10Nc  PMN/µL  protein   gluc  lactate  Q_alb   TNF  IL6/1e3   Pb   ICP   CPP   MAP  temp")
    for (t, y) in R["S01 untreated natural course"]["series"]:
        if t > 120:
            continue
        W(f"  {t:5.0f} {log10s(y[IX['Nc']]):7.2f} {y[IX['PMN']]:7.0f} {y[IX['Prot_csf']]:8.0f}"
          f" {y[IX['Glc_csf']]:6.1f} {y[IX['Lac_csf']]:8.1f} {1000*y[IX['Alb_csf']]/P['Alb_ser']:6.1f}"
          f" {y[IX['TNF']]:5.0f} {y[IX['IL6']]/1000:8.1f} {y[IX['Pb']]:4.1f}"
          f" {y[IX['ICP']]:5.1f} {y[IX['MAP']]-y[IX['ICP']]:5.1f} {y[IX['MAP']]:5.1f} {y[IX['Temp']]:5.1f}")
    W("  [literature, pneumococcal meningitis: CSF WBC 1,000-5,000/µL · protein 100-500 mg/dL · glucose <40 (ratio <0.4)")
    W("   · lactate >3.5 mmol/L · Q_alb 30-100 · CSF TNF 100-1,000 pg/mL · ICP often >20 mmHg]")
    W("")
    W("[1b] S24 standard treatment — does it recover (does the hazard accrue chronically)")
    W("-" * 112)
    W("   t[h] log10Nc  PMN/µL  protein   gluc  lactate   TNF   Mg   Pb R_out   ICP   CPP   MAP  SOFA")
    for (t, y) in R["S24 standard CEF+VAN+DEX"]["series"]:
        W(f"  {t:5.0f} {log10s(y[IX['Nc']]):7.2f} {y[IX['PMN']]:7.0f} {y[IX['Prot_csf']]:8.0f}"
          f" {y[IX['Glc_csf']]:6.1f} {y[IX['Lac_csf']]:8.1f} {y[IX['TNF']]:5.0f}"
          f" {y[IX['Mg']]:4.2f} {y[IX['Pb']]:4.1f} {y[IX['R_out']]:5.2f}"
          f" {y[IX['ICP']]:5.1f} {y[IX['MAP']]-y[IX['ICP']]:5.1f} {y[IX['MAP']]:5.1f}"
          f" {y[IX['SOFA']]:5.2f}")
    W("")

    # ---------------- [2] the lysis burst ----------------
    W("[2] The peak of the injury flux is 'right after the first antibiotic dose'  (S02 CEF alone, no DEX)")
    W("-" * 112)
    W("   t[h] CEF_CSF  C/MIC kkill/h log10Nc  lysflux    CW   PLY    Mg    TNF  PMN/µL")
    for (t, y) in R["S02 CEF alone (t=0)"]["series"]:
        if t > 24:
            continue
        kk = emax_kill(y[IX["Cef_csf"]], 0.03, P["Emax_cef"], P["EC50r_cef"], P["h_cef"])
        lf = (kk * P["Y_cef"] * y[IX["Nc"]] + P["k_autolysis"] * y[IX["Nc"]]) / 1e6
        W(f"  {t:5.1f} {y[IX['Cef_csf']]:7.2f} {y[IX['Cef_csf']]/0.03:6.0f} {kk:7.2f}"
          f" {log10s(y[IX['Nc']]):7.2f} {lf:8.2f} {y[IX['CW']]:5.1f} {y[IX['PLY']]:5.2f}"
          f" {y[IX['Mg']]:5.2f} {y[IX['TNF']]:6.0f} {y[IX['PMN']]:7.0f}")
    e2 = R["S02 CEF alone (t=0)"]["endpoints"]
    W(f"  The peak lysis flux appears within the first hour after dosing, and TNF rises to {e2['peak_TNF']:.0f} pg/mL.")
    W("  [literature: Mustafa 1989 and others — CSF TNF/IL-1 actually rises after the antibiotic]")
    W("")

    # ---------------- [3] dexamethasone timing ----------------
    W("[3] Dexamethasone timing — the shield must be up before the peak")
    W("-" * 112)
    W("  Scenario                   AUC_lysis AUC_TNF/1e3 peakTNF peakPMN peakPb hearDB  cog_z  dth%  unf%")
    for lab in ["S02 CEF alone (t=0)", "S03 CEF + DEX 20 min pre", "S04 CEF + DEX at 0 h",
                "S05 CEF + DEX +2 h", "S06 CEF + DEX +4 h", "S07 CEF + DEX +12 h",
                "S08 CEF + DEX 1 day only"]:
        e = R[lab]["endpoints"]
        W(f"  {lab:26s} {e['AUC_lysis']:9.1f} {e['AUC_TNF']/1000:11.1f} {e['peak_TNF']:7.0f}"
          f" {e['peak_PMN']:7.0f} {e['peak_Pb']:6.1f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    b, a = R["S02 CEF alone (t=0)"]["endpoints"], R["S03 CEF + DEX 20 min pre"]["endpoints"]
    c4 = R["S06 CEF + DEX +4 h"]["endpoints"]
    W(f"  → The benefit of the steroid itself is large: death {b['death']*100:.1f} % → {a['death']*100:.1f} %"
      f" (AUC_TNF −{(1-a['AUC_TNF']/b['AUC_TNF'])*100:.0f} %).")
    W(f"     But **the effect of the timing is almost nil under these conditions**:"
      f" −20 min {a['death']*100:.1f} % vs +4 h {c4['death']*100:.1f} %.")
    W(f"     Conversely **the duration of dosing** is clear: 4 days {a['death']*100:.1f} % vs 1 day"
      f" {R['S08 CEF + DEX 1 day only']['endpoints']['death']*100:.1f} %.")
    W("")
    W("  The same sweep in a severe patient (load 1e8 · antibiotic delayed 6 h) — what if the burst is large")
    W("  DEX time(h)  AUC_lysis AUC_TNF/1e3 peakTNF peakPMN peakICP hearDB  cog_z  dth%  unf%")
    for dt in [None, -0.33, 0.0, 2.0, 4.0, 8.0, 12.0, 24.0]:
        if dt is None:
            cfg = make_cfg("sev-nodex", cef_delay=6.0, N0=1.0e8, van=True)
            tag = "  none"
        else:
            cfg = make_cfg(f"sev{dt}", cef_delay=6.0, N0=1.0e8, van=True,
                           dex=True, dex_time=6.0 + dt)
            tag = f"{dt:+6.2f}"
        e = simulate(cfg, tmax=336.0)["endpoints"]
        W(f"  {tag:>11s} {e['AUC_lysis']:10.1f} {e['AUC_TNF']/1000:11.1f}"
          f" {e['peak_TNF']:7.0f} {e['peak_PMN']:7.0f} {e['peak_ICP']:7.1f}"
          f" {e['hear_dB']:6.1f} {e['cog_z']:6.2f} {e['death']*100:5.1f}"
          f" {e['unfav']*100:5.1f}")
        print(f"  ... severe sweep {tag}", file=sys.stderr)
    W("")
    W("  Interpretation — what this model predicts and what it does not:")
    W("  · The benefit of giving the steroid at all is reproduced at trial size (below, [9]).")
    W("  · But the part that says it must be **before** the first antibiotic dose comes out")
    W("    weak in this model.  The reason is structural: the cytokine burst is 2-6 h wide,")
    W("    while the hearing, cognitive and cortical injury integrals run over tens to")
    W("    hundreds of hours: the peak's share of a 14-day integral is small.  The model is")
    W("    instead far more sensitive to **duration** (1 vs 4 d) and to **antibiotic delay**.")
    W("  · This is not a convenient conclusion but a falsifiable prediction.  The likeliest candidate")
    W("    for the missing mechanism is a **threshold in the injury terms**: every injury term here")
    W("    is first order in its driver (dHC/dt ∝ driver), so 'a peak twice as high for 3 h' and 'a")
    W("    level 1.06 times as high for 100 h' do the same damage.  If real hair cell and neuronal")
    W("    death has a threshold, the peak becomes disproportionately expensive, and only then does")
    W("    'before the antibiotic' become arithmetically important.  Fitting a threshold to")
    W("    reproduce the recommendation is easy, but that is choosing structure to fit the conclusion, so it is left as a prediction.")
    W("")

    # ---------------- [4] antibiotic delay ----------------
    W("[4] Antibiotic delay — killing later means killing more bacteria, so the burst is bigger")
    W("-" * 112)
    W("  Scenario                   AUC_lysis peakICP minCPP  oedmL peakLac hearDB  dth%  unf%")
    for lab in ["S04 CEF + DEX at 0 h", "S09 antibiotic 3 h delay",
                "S10 antibiotic 6 h delay", "S11 antibiotic 12 h delay"]:
        e = R[lab]["endpoints"]
        W(f"  {lab:26s} {e['AUC_lysis']:9.1f} {e['peak_ICP']:7.1f} {e['min_CPP']:6.1f}"
          f" {e['peak_Vbr']:6.1f} {e['peak_Lac']:7.1f} {e['hear_dB']:6.1f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    W("  [literature: outcome worsens with every hour of door-to-antibiotic delay — Proulx 2005 (>6 h OR 8.4),")
    W("   Køster-Rasmussen 2008.  The model reproduces this mechanistically as 'a larger lysis burst']")
    W("")

    # ---------------- [5] closing vancomycin's door ----------------
    W("[5] Where the sign splits — dexamethasone also closes the door vancomycin comes in through")
    W("-" * 112)
    W("  Scenario                       AUC_van_CSF C_van(48h) T>4×MIC T_cef>4×MIC peakPb"
      "  sterT      adhT  dth%")
    for lab in ["S13 resistant + VAN, no DEX", "S14 resistant + VAN + DEX",
                "S15 resistant + VAN+RIF+DEX", "S25 resistant CEF CI+VAN+DEX",
                "S12 resistant(MIC 4) CEF alone", "S24 standard CEF+VAN+DEX"]:
        e = R[lab]["endpoints"]
        c48 = next((y[IX["Van_csf"]] for (tt, y) in R[lab]["series"] if abs(tt - 48) < 1e-6), 0.0)
        W(f"  {lab:30s} {e['AUC_van']:11.1f} {c48:10.2f} {e['T_van']:7.1f} {e['T_cef']:11.1f}"
          f" {e['peak_Pb']:6.1f} {fmt(e['sterile_t'],1):>6s} {fmt(e['adh_t'],1):>9s}"
          f" {e['death']*100:5.1f}")
    v13, v14 = R["S13 resistant + VAN, no DEX"]["endpoints"], R["S14 resistant + VAN + DEX"]["endpoints"]
    W(f"  → DEX tightens the barrier (peak Pb {v13['peak_Pb']:.1f} → {v14['peak_Pb']:.1f}) and cuts the vancomycin CSF AUC by"
      f" {(1-v14['AUC_van']/max(v13['AUC_van'],1e-9))*100:.0f} %,")
    W(f"     so CSF sterilisation is delayed from {fmt(v13['sterile_t'],1)} h to {fmt(v14['sterile_t'],1)} h."
      f"  The same manoeuvre does no harm with ceftriaxone (200-fold C/MIC headroom).")
    W(f"     Against the resistant strain, ceftriaxone alone exceeds 4×MIC (16 mg/L) for only "
      f"{R['S12 resistant(MIC 4) CEF alone']['endpoints']['T_cef']:.1f} h and fails to sterilise")
    W("     — the clinical case for adding vancomycin comes out of the output.")
    W("  [literature: Paris 1994 (rabbit, DEX reduces VAN CSF penetration) vs Ricard 2007 (human, kept at high dose)")
    W("   → the model quantifies the 'no headroom' side]")
    W("")

    # ---------------- [6] lytic vs non-lytic ----------------
    W("[6] Same log-kill, different injury integral — a non-lytic lead-in (rifampicin)")
    W("-" * 112)
    W("  Scenario                   sterT[h] AUC_lysis peakCW peakPLY peakTNF hearDB  cog_z  dth%")
    for lab in ["S02 CEF alone (t=0)", "S16 RIF 2 h lead-in → CEF", "S03 CEF + DEX 20 min pre"]:
        e, pk = R[lab]["endpoints"], R[lab]["peak"]
        W(f"  {lab:26s} {fmt(e['sterile_t'],1):>8s} {e['AUC_lysis']:9.1f} "
          f"{max(y[IX['CW']] for _, y in R[lab]['series']):6.1f} {e['peak_PLY']:7.2f}"
          f" {e['peak_TNF']:7.0f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f} {e['death']*100:5.1f}")
    W("  [literature: Nau/Böttcher — a non-lytic antibiotic lead-in reduces cell wall release and CSF inflammation]")
    W("")

    # ---------------- [7] all scenarios ----------------
    W("[7] Summary of the 26 scenarios")
    W("-" * 112)
    W("  Scenario                        sterT    adhT peakICP minCPP     minGlc peakLac Q_alb"
      " hearDB  cog_z hazAcute hazStruc  dth%  unf%")
    for label, _ in SCENARIOS:
        e = R[label]["endpoints"]
        W(f"  {label:30s} {fmt(e['sterile_t'],1):>6s} {fmt(e['adh_t'],1):>7s}"
          f" {e['peak_ICP']:7.1f} {e['min_CPP']:6.1f} {e['min_Glc']:10.1f} {e['peak_Lac']:7.1f}"
          f" {e['peak_Qalb']:5.0f} {e['hear_dB']:6.1f} {e['cog_z']:6.2f}"
          f" {e['haz_acute']:8.3f} {e['haz_struct']:8.3f}"
          f" {e['death']*100:5.1f} {e['unfav']*100:5.1f}")
    W("")

    # ---------------- [8] half-dt check ----------------
    W("[8] Integration convergence check (half dt) — same scenario, half the step")
    W("-" * 112)
    W("  Scenario                   endpoint         dt=1.0     dt=0.5 rel.diff%")
    for lab, kw in [("S02 CEF alone (t=0)", dict()),
                    ("S24 standard CEF+VAN+DEX", dict(van=True, dex=True, dex_time=-0.33))]:
        h = simulate(make_cfg(lab, **kw), tmax=336.0, dt_scale=0.5)["endpoints"]
        f = R[lab]["endpoints"]
        for key, nm in [("peak_ICP", "peak ICP"), ("hear_dB", "hearing dB"),
                        ("death", "P(death)"), ("AUC_lysis", "AUC_lysis")]:
            rel = abs(h[key] - f[key]) / max(abs(f[key]), 1e-12) * 100
            W(f"  {lab:26s} {nm:12s} {f[key]:10.4f} {h[key]:10.4f} {rel:9.3f}")
    W("")

    # ---------------- [9] virtual cohort ----------------
    W("[9] Virtual cohort of 10 × (DEX on/off) — against the European Dexamethasone Study")
    W("-" * 112)
    W("    pt log10 N0 delayH  hostDef  dth%(DEX-)  dth%(DEX+) hearDB(-) hearDB(+)  unf%(-)  unf%(+)")
    dm, dp, hm, hp, um, up = [], [], [], [], [], []
    for i, (n0, dl, hd) in enumerate(COHORT, 1):
        a = simulate(make_cfg("-", cef_delay=dl, N0=n0, host_def=hd, van=True), 336.0)["endpoints"]
        b = simulate(make_cfg("+", cef_delay=dl, N0=n0, host_def=hd, van=True,
                              dex=True, dex_time=dl - 0.33), 336.0)["endpoints"]
        dm.append(a["death"]); dp.append(b["death"])
        hm.append(a["hear_dB"]); hp.append(b["hear_dB"])
        um.append(a["unfav"]); up.append(b["unfav"])
        W(f"  {i:4d} {math.log10(n0):8.1f} {dl:6.1f} {hd:8.2f} {a['death']*100:11.1f}"
          f" {b['death']*100:11.1f} {a['hear_dB']:9.1f} {b['hear_dB']:9.1f}"
          f" {a['unfav']*100:8.1f} {b['unfav']*100:8.1f}")
        print(f"  ... cohort {i}", file=sys.stderr)
    n = len(COHORT)
    W(f"  Mean                          {sum(dm)/n*100:11.1f} {sum(dp)/n*100:11.1f}"
      f" {sum(hm)/n:9.1f} {sum(hp)/n:9.1f} {sum(um)/n*100:8.1f} {sum(up)/n*100:8.1f}")
    W("  [comparator de Gans & van de Beek NEJM 2002 pneumococcal subgroup: death 34 % → 14 %,")
    W("   unfavourable outcome 52 % → 26 %.  Paediatric Hib severe hearing loss ~15 % → ~5 %]")
    W("")

    # ---------------- [10] hazard decomposition ----------------
    W("[10] Component decomposition of the death hazard — which derangement actually kills")
    W("-" * 112)
    W("  Haz = h0·T + h_icp·I_icp + h_cpp·I_cpp + h_sofa·I_sofa + h_isch·I_isch")
    W("        + h_bact·I_bact + h_ncsf·I_ncsf  (acute)  +  h_cort·(1−N_cort)  (structural)")
    W("")
    W("  Scenario                          ICP     CPP    SOFA ischaem  bactmia  CSFinf   base  struct"
      "  total   dth%")
    for label in ["S01 untreated natural course", "S02 CEF alone (t=0)", "S03 CEF + DEX 20 min pre",
                  "S11 antibiotic 12 h delay", "S14 resistant + VAN + DEX",
                  "S20 early, low load 1e5 + DEX", "S24 standard CEF+VAN+DEX",
                  "S26 standard + EVD drainage"]:
        e = R[label]["endpoints"]
        c = dict(icp=P["h_icp"] * e["I_icp"], cpp=P["h_cpp"] * e["I_cpp"],
                 sofa=P["h_sofa"] * e["I_sofa"], isch=P["h_isch"] * e["I_isch"],
                 bact=P["h_bact"] * e["I_bact"], ncsf=P["h_ncsf"] * e["I_ncsf"],
                 base=P["h0"] * 336.0,
                 struct=P["h_cort_final"] * (1.0 - e["Ncort_f"]))
        tot = sum(c.values())
        W(f"  {label:30s} {c['icp']:6.3f} {c['cpp']:7.3f} {c['sofa']:7.3f} {c['isch']:7.3f}"
          f" {c['bact']:8.3f} {c['ncsf']:7.3f} {c['base']:6.3f} {c['struct']:7.3f} {tot:6.3f}"
          f" {e['death']*100:6.1f}")
    W("  → In a well-treated patient most of the remaining hazard is systemic organ failure (SOFA),")
    W("     while in a late or failed treatment the perfusion (CPP · ischaemia) and structural terms dominate.")
    W("")

    W("=" * 112)
    W("End.  This output is the result of running abm_reference_python.py, and it corresponds")
    W("1:1 to the equations and parameters of abm_mrgsolve_model_en.R.  [F1]-[F32] in the file")
    W("header are the defects this numerical check actually caught, and their corrections.")
    W("=" * 112)

    txt = "\n".join(L)
    print(txt)
    with open("abm_reference_output_en.txt", "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
