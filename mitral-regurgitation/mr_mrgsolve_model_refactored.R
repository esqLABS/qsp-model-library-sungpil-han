## =============================================================================
##  Mitral Regurgitation -- Quantitative Systems Pharmacology model
##  mrgsolve implementation, 43 ODEs
## =============================================================================
##
##  ORGANISING IDEA
##  ---------------
##  A regurgitant orifice produces ONE quantity that echocardiography reports --
##  the regurgitant volume RVol -- and that quantity is meaningless until it has
##  been divided by something.  Mitral regurgitation is therefore not one
##  disease but a family of diseases indexed by the DENOMINATOR:
##
##    RVol / C_LA(operating)  -> PRESSURE.  Sets congestion.  A virgin atrium
##       sits high on a steep pressure-volume curve, so 60 mL it has never seen
##       produces a giant v wave and floods the lung; a chronically dilated
##       atrium has shifted that curve right and flattened it, so the identical
##       60 mL is nearly silent.  Same lesion, opposite presentation.
##
##    RVol / SV_total  -> REGURGITANT FRACTION.  The grade.  Note the trap: in a
##       low-output patient a modest RVol is a large fraction, so the reported
##       severity rises as the ventricle fails.
##
##    RVol / LVEDV  -> PROPORTIONALITY.  Decides whether the VALVE or the
##       VENTRICLE is the disease, and therefore whether closing the valve can
##       possibly help.  This is the COAPT / MITRA-FR axis.
##
##    the afterload the leak REMOVES  -> CONTRACTILITY IS HIDDEN.  LVEF is
##       computed against a stroke volume that includes the leak, which is
##       ejected into a low-pressure sink.  EF therefore overstates the
##       ventricle: the model derives a post-operative EF of 0.50 at a
##       PRE-operative EF of 0.627, which is why the guideline threshold sits at
##       0.60 and not 0.50, and why EF FALLS after a successful repair.
##
##    EROA_true / k_PISA  -> THE MEASUREMENT HAS A DENOMINATOR TOO.  PISA assumes
##       a hemispheric flow convergence zone.  A degenerative orifice is round
##       and PISA is close to right; a functional orifice is a crescent along the
##       coaptation line and PISA overstates it.  The 0.4 cm2 severity threshold
##       was calibrated on the round one and is applied to the crescent.
##
##  STRUCTURE
##  ---------
##   Fast (algebraic, solved inside $ODE at every derivative evaluation):
##     closed-loop circulation -- time-varying-elastance LV and RV, systemic and
##     pulmonary Windkessels, exponential atrial pressure-volume law, and the
##     regurgitant path as a PARALLEL low-impedance drain competing with the
##     aorta.  Two nested monotone bisections: end-systolic pressure (inner) and
##     mean left atrial pressure closing total blood volume (outer).
##   Slow (the 43 ODE states):
##     1  LV chamber dilation, mass, collagen, contractility
##     2  LA dilation, fibrosis, atrial fibrillation burden
##     3  mitral annulus, leaflet area (rate-limited, PARTIAL adaptation)
##     4  degenerative orifice progression
##     5  pulmonary vascular resistance with an IRREVERSIBLE component
##     6  RV contractility and dilation
##     7  sympathetic, angiotensin, aldosterone, blood volume, eGFR, BNP
##     8  lung lymphatic adaptation of the oedema threshold
##     9  drug PK: furosemide, sacubitril/valsartan, beta-blocker, MRA, SGLT2i,
##        nitroprusside, dobutamine
##    10 device / surgical state: TEER, annuloplasty ring
##    11 endpoint accumulators: HF hospitalisation count, death hazard, symptoms
##
##  TIME UNIT: DAYS throughout.  The circulation equilibrates in seconds, so it
##  is treated as quasi-steady rather than integrated; only quantities with time
##  constants of hours or longer are ODE states.
##
##  CHRONIC DOSING enters as a continuous input RATE (mg/day) rather than as
##  discrete boluses.  The endpoints here have time constants of months, so only
##  average exposure matters, and this keeps the right-hand side smooth.  Use
##  rate_fur, rate_bb, rate_sac, rate_val, rate_mra, rate_sg.
##
##  CALIBRATION -- what is fitted and what is not
##  ---------------------------------------------
##   * The healthy subject is an EXACT fixed point of all 43 states.  The wall
##     stress, filling pressure and cardiac index setpoints are computed FROM the
##     healthy baseline rather than asserted, so no drift is possible.
##   * Baseline haemodynamics are resting physiology, not fitted: LVEDV 114 mL,
##     LVEF 0.620, CO 4.94 L/min, MAP 95, mean LAP 8.1, mean PAP 14.1, CVP 4.8,
##     LA volume 53 mL, PVR 1.2 WU, end-diastolic wall stress 12.7 mmHg,
##     end-systolic 81.8 mmHg, wall thickness 0.94 cm.
##   * ONE unstressed volume (V0sv) was solved to place mean LAP at 8.1 mmHg, and
##     four waveform-shape constants (c_ea, c_rv, cflow, c_s) fold the
##     within-beat pressure and velocity profiles into beat averages.
##   * TWO numbers in the whole model are fitted to OUTCOME data, and both come
##     from the COAPT CONTROL arm: h0_hfh and h0_death.  The hazard slopes are
##     set a priori from published hazard ratios (~1.22 per 5 mmHg of filling
##     pressure), NOT fitted.
##   * The two trial-matched virtual patients are built by forward simulation
##     from an index infarct; three knobs are fitted to three reported numbers
##     (LVEDV, LVEF, PISA-EROA) and everything else is a consequence.
##
##  KNOWN MISSES -- stated here rather than buried
##  ---------------------------------------------
##   1. THE CENTRAL MISS.  The model predicts a real benefit from edge-to-edge
##      repair in a MITRA-FR-like patient (HF admission 0.276 -> 0.147 at 12
##      months) where the trial found none.  It reproduces the ORDERING and its
##      mechanism -- the same device buys the COAPT patient 8.5% more forward
##      output and the MITRA-FR patient 2.2%, a 3.8-fold difference arising
##      purely from the orifice-to-ventricle ratio -- but the trials show a
##      change of SIGN, not a 3.8-fold difference.  An additive prognostic term
##      cannot repair this (additive terms cancel in a hazard ratio), and even
##      with NO orifice reduction the model sits below the observed MITRA-FR
##      event rate.  Either that device arm achieved much less than is assumed
##      here, or a large part of that cohort's risk lived outside any valve
##      model.  No attempt has been made to hide this by re-fitting.
##   2. Chronic primary MR in this model becomes congested faster than the
##      classic asymptomatic severe phenotype, because atrial remodelling is
##      slower than the ventricular feedback loop in the present calibration.
##      Wedge pressures in the multi-year primary-MR runs should be read as
##      qualitative.
##   3. The virtual MITRA-FR patient reaches an LVEDV of 237 mL against a
##      reported 256 mL (7.5% short).  In this model's geometry a ventricle that
##      large tethers its leaflets enough to open a wider orifice than the 0.31
##      cm2 reported, and the residual is the price of not forcing it.
##   4. Chamber geometry is spherical.  Real ventricles are ellipsoidal and real
##      functional orifices are crescentic; the sphere is used for wall stress
##      and papillary displacement, and the crescent enters only through k_PISA.
##   5. The regurgitant orifice is treated as constant within systole.  It is
##      demonstrably biphasic in secondary MR; only its systolic average is
##      represented.
##   6. Heart rate enters through beat timing and filling but there is no
##      explicit atrial-ventricular optimisation, so the AF penalty (a fixed 8%
##      loss of filling) is a lumped approximation.
##   7. The device is a step change in orifice area plus an added mitral
##      resistance.  Leaflet grasping geometry, single versus double clip and
##      device-specific durability are not represented.
##
##  VERIFICATION
##  ------------
##  Every equation in this file was independently re-implemented in Python/scipy
##  (mr_python_reference.py) and the reference output is committed as
##  mr_reference_output.txt.  That cross-check found and fixed EIGHT defects,
##  listed in README.md.  NOTE HONESTLY: no R toolchain was available in the
##  environment where this file was written, so the R code has NOT been executed.
##  The Python implementation is the executed reference; this file mirrors it
##  equation for equation and its numerical agreement has not been demonstrated.
##
##  USAGE
##  -----
##    library(mrgsolve); library(dplyr)
##    mod <- mread("mr_mrgsolve_model.R")
##
##    # healthy baseline: an exact fixed point
##    mod %>% mrgsim(end = 3650, delta = 10) %>% plot(EF + CI + Ppcw ~ time)
##
##    # acute papillary muscle rupture onto a virgin atrium
##    mod %>% param(EROApri_init = 0.60) %>% init(EROApri = 0.60) %>%
##        mrgsim(end = 30, delta = 0.25) %>% plot(RVol + vwave + Ppcw + CI ~ time)
##
##    # a COAPT-like patient, with and without edge-to-edge repair
##    coapt <- mod %>% init(Ees = 0.7028, V0d = 45.04, Aleaf = 10.42)
## =============================================================================

$PROB
# Mitral regurgitation: 43-ODE QSP model
- One regurgitant volume, five denominators.
- Closed-loop elastance circulation solved algebraically inside the ODE block,
  coupled to growth, remodelling, neurohormonal and device dynamics.
- Author: QSP disease model library
- Units: TIME IN DAYS.  Pressures mmHg, volumes mL, areas cm2,
  resistances mmHg.s/mL, compliances mL/mmHg.

$PARAM @annotated
// ---------------- subject -------------------------------------------------
BSA       : 1.90   : Body surface area (m2)
HR0       : 70.0   : Intrinsic heart rate (1/min)

// ---------------- LV chamber mechanics ------------------------------------
Ees0      : 2.90   : LV end-systolic elastance, healthy (mmHg/mL)
V0s_off   : 22.0   : ESPVR intercept offset, V0s = V0d - V0s_off (mL)
Aed       : 1.00   : EDPVR pressure scale (mmHg)
Bed0      : 0.0260 : EDPVR stiffness, healthy (1/mL)
V0d_h     : 30.0   : EDPVR volume intercept, healthy (mL)
V0d_max   : 240.0  : Ceiling on LV chamber dilation (mL)
LVmass_h  : 150.0  : LV mass, healthy (g)

// ---------------- waveform-shape constants (the only shape fits) ----------
c_ea      : 1.077  : Ea = c_ea * Rsys / T, reconciles Ea=Pes/SV with Rsys (-)
c_rv      : 2.20   : RV afterload elastance scale, mean PA -> end-systolic (-)
cflow     : 0.80   : Regurgitant waveform factor, <v> vs sqrt(<dP>) (-)
c_s       : 0.90   : Mean systolic LV pressure as a fraction of Pes (-)
t_iv      : 0.070  : Regurgitant time beyond ejection (s)
f_v       : 0.50   : Fraction of the LA v wave transmitted to the capillary (-)
xi_v      : 0.50   : Fraction of RVol retained in the LA during systole (-)

// ---------------- resistances and compliances -----------------------------
Rsys0     : 1.154  : Systemic vascular resistance, healthy (mmHg.s/mL)
Rpul0     : 0.0720 : Pulmonary vascular resistance, healthy = 1.2 WU
Rmv0      : 0.0020 : Native mitral resistance (mmHg.s/mL)
Csa       : 1.60   : Systemic arterial compliance (mL/mmHg)
Csv0      : 110.0  : Systemic venous compliance (mL/mmHg)
Cpa       : 5.00   : Pulmonary arterial compliance (mL/mmHg)
Cpu       : 14.0   : Pulmonary venous + capillary compliance (mL/mmHg)
V0sa      : 600.0  : Unstressed systemic arterial volume (mL)
V0sv      : 3060.646 : Unstressed systemic venous volume, SOLVED for LAP 8.1 (mL)
V0pa      : 50.0   : Unstressed pulmonary arterial volume (mL)
V0pu      : 200.0  : Unstressed pulmonary venous volume (mL)
Vtot_h    : 5000.0 : Total blood volume, healthy (mL)

// ---------------- left atrium: DENOMINATOR ONE ----------------------------
A_la      : 3.00   : LA pressure-volume scale (mmHg)
B_la0     : 0.0625 : LA pressure-volume stiffness, healthy (1/mL)
kB_la     : 0.565  : LA dilation softens the curve (exponent)
V0la_h    : 32.0   : LA unstressed volume, healthy (mL)
V0la_max  : 220.0  : Ceiling on LA dilation (mL)
s_fibla   : 1.20   : LA fibrosis stiffening at maximal fibrosis (-)

// ---------------- RV ------------------------------------------------------
Ees_rv0   : 0.500  : RV end-systolic elastance, healthy (mmHg/mL)
V0s_rv    : 20.0   : RV ESPVR volume intercept (mL)
Aed_rv    : 1.00   : RV EDPVR pressure scale (mmHg)
Bed_rv    : 0.0170 : RV EDPVR stiffness (1/mL)
V0d_rv_h  : 30.0   : RV EDPVR volume intercept, healthy (mL)

// ---------------- mitral valve geometry -----------------------------------
Ann_h     : 6.50   : Mitral annular area, healthy (cm2)
Ann_max   : 15.0   : Ceiling on annular area (cm2)
Aleaf_h   : 10.08  : Leaflet area available for coaptation, healthy (cm2)
CD_h      : 0.55   : Coaptation depth / tenting, healthy (cm)
kt        : 1.30   : Tethering -> required coaptation area (1/cm)
k_o       : 0.300  : Unmet coaptation area -> EROA (cm2/cm2)
kSI       : 0.55   : Dilation -> sphericity (-)
SI_h      : 0.62   : Sphericity index, healthy (-)
kCD       : 0.85   : Sphericity exponent in tenting (-)
k_pisa_sec : 0.45  : PISA inflation for a fully crescentic orifice (-)

// ---------------- remodelling gains ---------------------------------------
k_dil     : 0.00300 : ED wall stress -> chamber dilation (mL/(mmHg.day))
k_dil_rev : 0.45    : Regression is this fraction as fast, hysteresis (-)
k_hyp     : 0.01200 : Mass relaxation toward its target (1/day)
k_fib     : 0.00160 : Collagen deposition (1/day)
Fib_h     : 0.055   : Collagen volume fraction, healthy (-)
Fib_max   : 0.42    : Maximal collagen volume fraction (-)
s_fib     : 1.60    : Collagen -> EDPVR stiffness at Fib_max (-)
e_ang_fib : 0.55    : Angiotensin drive on fibrosis (-)
e_ald_fib : 0.45    : Aldosterone drive on fibrosis (-)
e_str_fib : 0.60    : Systolic stress drive on fibrosis (-)
k_eesloss : 0.000105 : Contractility loss per unit damage drive (1/day)
k_eesrec  : 0.000230 : Contractility recovery (1/day)
Ees_min   : 0.35    : Contractility floor (mmHg/mL)

// ---------------- LA and atrial fibrillation ------------------------------
k_la      : 0.00800 : LA pressure -> LA dilation (mL/(mmHg.day))
la_dead   : 3.00    : LA remodelling deadband (mmHg)
k_fibla   : 0.00125 : LA fibrosis rate (1/day)
k_af      : 0.00135 : AF substrate -> AF burden (1/day)
k_af_rev  : 0.00060 : AF regression (1/day)

// ---------------- annulus and leaflets ------------------------------------
k_ann     : 0.00135 : Dilation drive -> annular area (cm2/day)
k_ann_rev : 0.00050 : Annular regression (1/day)
k_ann_af  : 0.28    : AF contribution to annular dilation (-)
phi_leaf  : 0.35    : Fraction of excess coaptation demand leaflets recover (-)
k_leaf    : 0.00420 : Leaflet growth rate constant (1/day)
leaf_cap  : 0.00340 : Hard cap on leaflet growth (cm2/day)
Aleaf_max_f : 1.20  : Leaflets grow at most 20% above native area (-)
k_pri_prog : 0.0    : Degenerative orifice progression rate (1/day)

// ---------------- pulmonary vasculature and RV ----------------------------
k_pvr     : 4.20e-5 : LAP above 18 -> PVR rise ((mmHg.s/mL)/(mmHg.day))
k_pvr_rev : 0.00120 : PVR regression (1/day)
f_pvr_fix : 0.34    : Fraction of the PVR rise that NEVER regresses (-)
Rpul_max  : 0.62    : Ceiling on PVR (mmHg.s/mL)
k_rvloss  : 0.000135 : RV contractility loss per unit PA afterload drive (1/day)
k_rvrec   : 0.000210 : RV contractility recovery (1/day)
k_rvdil   : 0.00230 : RV dilation (mL/(mmHg.day))

// ---------------- neurohormonal -------------------------------------------
tau_ne    : 0.30   : Sympathetic time constant (day)
tau_ang   : 0.18   : Angiotensin time constant (day)
tau_ald   : 0.30   : Aldosterone time constant (day)
k_ne_co   : 0.85   : Baroreflex gain on cardiac index (-)
k_ang_co  : 1.10   : Baroreflex gain on angiotensin (-)
k_ald_ang : 1.00   : Angiotensin -> aldosterone gain (-)

// ---------------- volume and kidney ---------------------------------------
k_na      : 120.0  : Sodium/water retention per unit RAAS drive (mL/day)
k_na_esc  : 0.0550 : Pressure natriuresis / escape (1/day)
dVtot_max : 2200.0 : Largest sustainable intravascular overload (mL)
Vtot_min_f : 0.90  : Natriuresis fades out at this fraction of Vtot_h (-)
eGFR_h    : 72.0   : eGFR, healthy (mL/min/1.73m2)
k_gfr_perf : 0.62  : Forward perfusion term on eGFR (-)
k_gfr_cong : 0.85  : Venous congestion term on eGFR (-)
tau_gfr   : 9.0    : eGFR time constant (day)
tau_bnp   : 0.42   : BNP time constant (day)
k_bnp     : 1.35   : BNP gain (-)
Pcrit_h   : 20.0   : Alveolar oedema threshold, unadapted (mmHg)
k_lymph   : 0.00460 : Lymphatic adaptation rate (1/day)
Pcrit_max : 32.0   : Ceiling on the oedema threshold (mmHg)

// ---------------- endpoints -----------------------------------------------
h0_hfh    : 1.669505e-4 : HF hospitalisation baseline hazard, FITTED on COAPT control (1/day)
h0_death  : 1.139320e-4 : Death baseline hazard, FITTED on COAPT control (1/day)
a1_hfh    : 0.20   : Wedge slope for HF hospitalisation, a priori (-)
a2_hfh    : 0.20   : Cardiac index slope for HF hospitalisation, a priori (-)
b1_d      : 0.15   : Wedge slope for death, a priori (-)
b2_d      : 0.18   : Cardiac index slope for death, a priori (-)
b3_d      : 0.15   : CVP slope for death, a priori (-)
c1_ved    : 0.50   : LVEDV index slope, a priori; additive so it cannot change any HR (-)

// ---------------- drug PK -------------------------------------------------
KA_FUR    : 31.2   : Furosemide absorption (1/day)
CL_FUR    : 216.0  : Furosemide clearance (L/day)
V1_FUR     : 14.0   : Furosemide volume (L)
F_FUR     : 0.55   : Furosemide bioavailability (-)
EMAX_FUR  : 400.0  : Maximal furosemide natriuresis (mL/day)
EC50_FUR  : 1.10   : Furosemide EC50 (mg/L)
KBRAKE_FUR   : 0.55   : Diuretic braking rate (1/day)
KA_SAC    : 26.4   : Sacubitril absorption (1/day)
CL_SAC    : 8.40   : Sacubitrilat clearance (L/day)
V1_SAC     : 7.00   : Sacubitrilat volume (L)
F_SAC     : 0.60   : Sacubitril bioavailability (-)
EC50_SAC  : 0.90   : Sacubitrilat EC50 (mg/L)
KA_VAL    : 24.0   : Valsartan absorption (1/day)
CL_VAL    : 52.8   : Valsartan clearance (L/day)
V1_VAL     : 17.0   : Valsartan volume (L)
F_VAL     : 0.23   : Valsartan bioavailability (-)
EC50_VAL  : 1.60   : Valsartan EC50 (mg/L)
KA_BB     : 10.8   : Beta-blocker absorption (1/day)
CL_BB     : 1512.0 : Beta-blocker clearance (L/day)
V1_BB      : 290.0  : Beta-blocker volume (L)
F_BB      : 0.40   : Beta-blocker bioavailability (-)
EC50_BB   : 0.045  : Beta-blocker EC50 (mg/L)
KA_MRA    : 24.0   : Spironolactone absorption (1/day)
CL_MRA    : 45.6   : Canrenone clearance (L/day)
V1_MRA     : 60.0   : Canrenone volume (L)
F_MRA     : 0.70   : Spironolactone bioavailability (-)
EC50_MRA  : 0.055  : Canrenone EC50 (mg/L)
KA_SG     : 31.2   : SGLT2 inhibitor absorption (1/day)
CL_SG     : 312.0  : SGLT2 inhibitor clearance (L/day)
V1_SG      : 118.0  : SGLT2 inhibitor volume (L)
F_SG      : 0.78   : SGLT2 inhibitor bioavailability (-)
EC50_SG   : 0.030  : SGLT2 inhibitor EC50 (mg/L)
KE_SNP    : 60.0   : Nitroprusside effect-site rate (1/day)
EC50_SNP  : 1.00   : Nitroprusside EC50 (effect units)
KE_DOB    : 140.0  : Dobutamine effect-site rate (1/day)
EC50_DOB  : 1.00   : Dobutamine EC50 (effect units)

// ---------------- drug maximal effects ------------------------------------
EMAX_BB_HR  : 0.285 : Beta-blocker maximal heart rate reduction (-)
EMAX_BB_EES : 0.150 : Beta-blocker acute negative inotropy (-)
EMAX_BB_REC : 1.55  : Beta-blocker chronic contractility recovery gain (-)
EMAX_VAL_RSYS : 0.180 : Valsartan maximal SVR reduction (-)
EMAX_VAL_ANG  : 0.700 : Valsartan maximal AT1 blockade (-)
EMAX_SAC_CSV  : 0.250 : Neprilysin inhibition venodilation (-)
EMAX_SAC_NA   : 0.320 : Neprilysin inhibition natriuresis (-)
EMAX_SAC_FIB  : 0.260 : Neprilysin inhibition antifibrotic (-)
EMAX_MRA_FIB  : 0.520 : MRA antifibrotic (-)
EMAX_MRA_NA   : 0.150 : MRA natriuresis (-)
EMAX_SG_VOL   : 0.065 : SGLT2 inhibition volume setpoint reduction (-)
EMAX_SG_FIB   : 0.130 : SGLT2 inhibition antifibrotic (-)
EMAX_SNP_RSYS : 0.420 : Nitroprusside maximal SVR reduction (-)
EMAX_SNP_CSV  : 0.300 : Nitroprusside venodilation (-)
EMAX_DOB_EES  : 0.900 : Dobutamine maximal inotropy (-)
EMAX_DOB_HR   : 0.220 : Dobutamine chronotropy (-)

// ---------------- Hill coefficients: renamed, not fitted (all originally --
// ------------ plain concentration/EC50 ratios, i.e. implicit gamma = 1) --
GAMMA_BB  : 1.0    : Beta-blocker Hill coefficient (-)
GAMMA_VAL : 1.0    : Valsartan Hill coefficient (-)
GAMMA_SAC : 1.0    : Sacubitrilat Hill coefficient (-)
GAMMA_MRA : 1.0    : Canrenone Hill coefficient (-)
GAMMA_SG  : 1.0    : SGLT2 inhibitor Hill coefficient (-)
GAMMA_SNP : 1.0    : Nitroprusside Hill coefficient (-)
GAMMA_DOB : 1.0    : Dobutamine Hill coefficient (-)
GAMMA_FUR : 1.0    : Furosemide Hill coefficient (-)

// ---------------- interventions and infusions -----------------------------
rate_fur  : 0.0    : Furosemide input rate (mg/day)
rate_bb   : 0.0    : Beta-blocker input rate (mg/day)
rate_sac  : 0.0    : Sacubitril input rate (mg/day)
rate_val  : 0.0    : Valsartan input rate (mg/day)
rate_mra  : 0.0    : Spironolactone input rate (mg/day)
rate_sg   : 0.0    : SGLT2 inhibitor input rate (mg/day)
snp_inf   : 0.0    : Nitroprusside infusion (effect units)
dob_inf   : 0.0    : Dobutamine infusion (effect units)
crt       : 0.0    : Cardiac resynchronisation, 0 or 1 (-)
HR_fix    : -1     : Clamp heart rate; negative disables (1/min)
f_Rsys    : 1.0    : External scaling of systemic resistance (-)
Rmv_add   : 0.0    : External addition to mitral resistance (mmHg.s/mL)
dV0d_force : -1    : Forced chamber dilation rate; negative disables (mL/day)
V0d_stop  : 1e9    : Stop forced dilation at this V0d (mL)
EROApri_init : 0.0    : Convenience copy of the initial degenerative orifice (cm2)

$CMT @annotated
V0d      : LV EDPVR volume intercept, chamber dilation (mL)
LVmass   : LV mass (g)
Fib      : LV collagen volume fraction (-)
Ees      : LV end-systolic elastance (mmHg/mL)
V0la     : LA unstressed volume (mL)
Fibla    : LA collagen volume fraction (-)
AFb      : Atrial fibrillation burden (0-1)
Ann      : Mitral annular area (cm2)
Aleaf    : Leaflet area available for coaptation (cm2)
EROApri  : Degenerative regurgitant orifice (cm2)
Rpul     : Pulmonary vascular resistance (mmHg.s/mL)
Rpulfix  : Irreversible component of PVR (mmHg.s/mL)
Eesrv    : RV end-systolic elastance (mmHg/mL)
V0drv    : RV EDPVR volume intercept (mL)
NE       : Sympathetic drive index (1 = normal)
Ang      : Angiotensin II index (1 = normal)
Ald      : Aldosterone index (1 = normal)
Vtot     : Total blood volume (mL)
eGFR     : eGFR (mL/min/1.73m2)
BNP      : BNP index (1 = normal)
Pcrit    : Alveolar oedema threshold (mmHg)
GUT_FUR   : Furosemide depot (mg)
CENT_FUR   : Furosemide central (mg)
BRAKE_FUR    : Diuretic braking index (-)
GUT_SAC   : Sacubitril depot (mg)
CENT_SAC   : Sacubitrilat central (mg)
GUT_VAL   : Valsartan depot (mg)
CENT_VAL   : Valsartan central (mg)
GUT_BB    : Beta-blocker depot (mg)
CENT_BB    : Beta-blocker central (mg)
GUT_MRA   : Spironolactone depot (mg)
CENT_MRA   : Canrenone central (mg)
GUT_SG    : SGLT2 inhibitor depot (mg)
CENT_SG    : SGLT2 inhibitor central (mg)
CENT_SNP   : Nitroprusside effect site (-)
CENT_DOB   : Dobutamine effect site (-)
BBDUR_BB    : Cumulative beta-blockade exposure (day)
HFH      : Cumulative expected HF hospitalisations (-)
CumHz    : Cumulative death hazard (-)
NYHAi    : Symptom index (1-4.4)
TEERrmv  : Mitral resistance added by device (mmHg.s/mL)
TEERf    : Fraction of orifice abolished by device (-)
RingF    : Annuloplasty ring flag, 1 = annulus clamped (-)

$MAIN
_F(1) = 1.0;   // depots receive continuous rates, set in $ODE

// ---- initial conditions: the healthy subject, an exact fixed point --------
V0d_0      = V0d_h;
LVmass_0   = LVmass_h;
Fib_0      = Fib_h;
Ees_0      = Ees0;
V0la_0     = V0la_h;
Fibla_0    = Fib_h;
AFb_0      = 0.0;
Ann_0      = Ann_h;
Aleaf_0    = Aleaf_h;
EROApri_0  = EROApri_init;
Rpul_0     = Rpul0;
Rpulfix_0  = 0.0;
Eesrv_0    = Ees_rv0;
V0drv_0    = V0d_rv_h;
NE_0       = 1.0;
Ang_0      = 1.0;
Ald_0      = 1.0;
Vtot_0     = Vtot_h;
eGFR_0     = eGFR_h;
BNP_0      = 1.0;
Pcrit_0    = Pcrit_h;
NYHAi_0    = 1.0;

$GLOBAL
#define SQ(x) ((x)*(x))

// Collagen turnover balanced analytically so that Fib_h is an EXACT fixed
// point rather than an approximate one.
#define KFIBDEG (k_fib   * (Fib_max - Fib_h) / Fib_h)
#define KFIBLADEG (k_fibla * (Fib_max - Fib_h) / Fib_h)

// Healthy setpoints.  These are the values the algebraic circulation returns at
// the healthy state and are hard-coded here so that the growth laws have an
// exact zero.  They were computed FROM the baseline, not asserted: see
// make_setpoints() in mr_python_reference.py.
#define SIG_ED_SET 12.65531
#define SIG_ES_SET 81.83330
#define PPCW_SET    8.13000
#define PLA_SET     8.13000
#define PPA_SET    14.06430
#define PSV_SET     4.76080
#define CI_SET      2.60277
#define VED_SET   113.96100

namespace {
  double sphere_r(double V) {
    if (V < 1.0) V = 1.0;
    return pow(3.0 * V / (4.0 * M_PI), 1.0 / 3.0);
  }
  double wall_h(double Vcav, double Vwall) {
    double r = sphere_r(Vcav);
    if (Vwall < 1.0) Vwall = 1.0;
    double o3 = r * r * r + 3.0 * Vwall / (4.0 * M_PI);
    return pow(o3, 1.0 / 3.0) - r;
  }
}

// Every quantity the algebraic circulation produces, shared between $ODE and
// $TABLE so the beat is solved once per evaluation.
double g_HR, g_T, g_Tr, g_Ea, g_Rsys, g_Csv, g_Cla, g_Bla, g_Bed;
double g_EROA, g_EROApisa, g_CD, g_Areq, g_resv, g_Rmv;
double g_Pla, g_Ppcw, g_vwave, g_Vla, g_Ved, g_Ves, g_Ped, g_Pes;
double g_SVf, g_RVol, g_SVtot, g_EF, g_EFfwd, g_RF, g_CO, g_CI;
double g_MAP, g_Ppa, g_Psv, g_Vrved, g_sig_ed, g_sig_es, g_h_ed, g_dPmv;
double g_LVEDVi, g_PVRWU, g_VTImr;

$ODE
// =========================================================================
//  FAST BLOCK.  The circulation is quasi-steady on a daily time scale, so it
//  is solved algebraically here rather than integrated.  Two nested monotone
//  bisections: end-systolic pressure (inner) and mean LA pressure closing
//  total blood volume (outer).  Both residuals are strictly monotone, so each
//  root is unique and bisection cannot converge to the wrong one.
// =========================================================================

// ---- drug concentrations: one exposed C_<STEM> per compound ---------------
double C_BB  = CENT_BB  / V1_BB;
double C_VAL = CENT_VAL / V1_VAL;
double C_SAC = CENT_SAC / V1_SAC;
double C_MRA = CENT_MRA / V1_MRA;
double C_SG  = CENT_SG  / V1_SG;
double C_FUR = CENT_FUR / V1_FUR;
double C_SNP = CENT_SNP;   // effect-site compartment is itself the exposed concentration
double C_DOB = CENT_DOB;   // effect-site compartment is itself the exposed concentration

// ---- drug effects: one named EFFECT_<STEM>[_x] per compound's own action --
// OCC_<STEM> is the compound's raw (Emax=1) Hill fraction, shared across a
// compound's several named effects where the original scaled the same
// concentration ratio by more than one Emax (BB, SAC, MRA, SG).
double OCC_BB  = pow(C_BB,  GAMMA_BB)  / (pow(EC50_BB,  GAMMA_BB) + pow(C_BB,  GAMMA_BB));
double OCC_SAC = pow(C_SAC, GAMMA_SAC) / (pow(EC50_SAC, GAMMA_SAC) + pow(C_SAC, GAMMA_SAC));
double OCC_MRA = pow(C_MRA, GAMMA_MRA) / (pow(EC50_MRA, GAMMA_MRA) + pow(C_MRA, GAMMA_MRA));
double OCC_SG  = pow(C_SG,  GAMMA_SG)  / (pow(EC50_SG,  GAMMA_SG) + pow(C_SG,  GAMMA_SG));

double EFFECT_BB_HR   = EMAX_BB_HR   * OCC_BB;
double EFFECT_BB_EES  = EMAX_BB_EES  * OCC_BB;
double EFFECT_VAL_RSYS = EMAX_VAL_RSYS * pow(C_VAL, GAMMA_VAL) / (pow(EC50_VAL, GAMMA_VAL) + pow(C_VAL, GAMMA_VAL));
double EFFECT_SAC_CSV = EMAX_SAC_CSV * OCC_SAC;
double EFFECT_SAC_NA  = EMAX_SAC_NA  * OCC_SAC;
double EFFECT_SAC_FIB = EMAX_SAC_FIB * OCC_SAC;
double EFFECT_MRA_FIB = EMAX_MRA_FIB * OCC_MRA;
double EFFECT_MRA_NA  = EMAX_MRA_NA  * OCC_MRA;
double EFFECT_SG_VOL  = EMAX_SG_VOL  * OCC_SG;
double EFFECT_SG_FIB  = EMAX_SG_FIB  * OCC_SG;
double EFFECT_SNP_RSYS = EMAX_SNP_RSYS * pow(C_SNP, GAMMA_SNP) / (pow(EC50_SNP, GAMMA_SNP) + pow(C_SNP, GAMMA_SNP));
double EFFECT_SNP_CSV  = EMAX_SNP_CSV  * pow(C_SNP, GAMMA_SNP) / (pow(EC50_SNP, GAMMA_SNP) + pow(C_SNP, GAMMA_SNP));
double EFFECT_DOB_EES = EMAX_DOB_EES * pow(C_DOB, GAMMA_DOB) / (pow(EC50_DOB, GAMMA_DOB) + pow(C_DOB, GAMMA_DOB));
double EFFECT_DOB_HR  = EMAX_DOB_HR  * pow(C_DOB, GAMMA_DOB) / (pow(EC50_DOB, GAMMA_DOB) + pow(C_DOB, GAMMA_DOB));

// ---- beat timing ---------------------------------------------------------
double HR = HR0 * (1.0 + 0.30 * (NE - 1.0)) * (1.0 - EFFECT_BB_HR)
            * (1.0 + EFFECT_DOB_HR) * (1.0 + 0.16 * AFb);
if (HR < 38.0)  HR = 38.0;
if (HR > 175.0) HR = 175.0;
if (HR_fix > 0) HR = HR_fix;
double T = 60.0 / HR;
double LVET = 0.42 - 0.0016 * HR;
if (LVET < 0.150) LVET = 0.150;
double Tr = LVET + t_iv;
double Tdias = T - LVET;  if (Tdias < 0.10) Tdias = 0.10;

// ---- effective parameters after drugs ------------------------------------
double Rsys = Rsys0 * (1.0 + 0.34 * (Ang - 1.0) + 0.20 * (NE - 1.0))
              * (1.0 - EFFECT_VAL_RSYS)
              * (1.0 - EFFECT_SNP_RSYS) * f_Rsys;
if (Rsys < 0.25) Rsys = 0.25;
double Csv = Csv0 * (1.0 + EFFECT_SAC_CSV + EFFECT_SNP_CSV);
double EesE = Ees * (1.0 - EFFECT_BB_EES) * (1.0 + EFFECT_DOB_EES)
              * (1.0 + 0.10 * (NE - 1.0));
if (EesE < 0.10) EesE = 0.10;
double Ea  = c_ea * Rsys / T;
double V0s = V0d - V0s_off;  if (V0s < 2.0) V0s = 2.0;
double Bed = Bed0 * (1.0 + s_fib * (Fib - Fib_h) / (Fib_max - Fib_h));
if (Bed < 0.004) Bed = 0.004;
double Vwall = LVmass / 1.05;

// ---- LEFT ATRIUM: denominator one ----------------------------------------
// Exponential pressure-volume law.  The OPERATING (incremental) compliance is
// what divides the regurgitant volume, and it is small in a virgin atrium high
// on a steep curve, large in a dilated one whose curve has shifted right and
// flattened, and small again once atrial fibrosis stiffens it back.
double Bla = B_la0 * pow(V0la_h / (V0la > 1.0 ? V0la : 1.0), kB_la)
             * (1.0 + s_fibla * (Fibla - Fib_h) / (Fib_max - Fib_h));
if (Bla < 1e-4) Bla = 1e-4;
if (Bla > 0.30) Bla = 0.30;

// ---- VALVE GEOMETRY: the orifice as a geometric DEFICIT ------------------
double SI = SI_h + kSI * (V0d / V0d_h - 1.0);  if (SI > 0.95) SI = 0.95;
double r_ratio = sphere_r(V0d + 70.0) / sphere_r(V0d_h + 70.0);
double CD = CD_h * r_ratio * pow(SI / SI_h, kCD);
double Areq_l = Ann * (1.0 + kt * (CD > CD_h ? (CD - CD_h) : 0.0));
double resv = Aleaf - Areq_l;                    // COAPTATION RESERVE
double EROAsec = k_o * (resv < 0.0 ? -resv : 0.0);
double EROAt = (EROApri + EROAsec) * (1.0 - TEERf);
if (EROAt < 0.0) EROAt = 0.0;
double cap = 0.55 * Ann;
EROAt = cap * tanh(EROAt / cap);               // an orifice cannot exceed its annulus
EROAt *= (1.0 - 0.22 * crt);                   // resynchronised closing force
double frac_sec = EROAsec / ((EROApri + EROAsec) > 1e-9 ? (EROApri + EROAsec) : 1e-9);
double kpisa = 1.0 + k_pisa_sec * frac_sec;

double Rmv = Rmv0 + TEERrmv + Rmv_add;
double Kv  = 50.0 * cflow;

// ---- nested bisection ----------------------------------------------------
double Pla = 0.2;
double Pla_hi = 140.0;
double Ved=0;
double Pes_l=0;
double SVf=0;
double RVol_l=0;
double SVtot_l=0;
double Ped_l=0;
double Q=0;
double MAP_l=0;
double Ppa=0;
double Psv=0;
double Vrved=0;

for (int outer = 0; outer < 80; ++outer) {
  double Pmid = 0.5 * (Pla + Pla_hi);
  // ---- inner: solve the beat at this atrial pressure ---------------------
  double SVt_ = 60.0;
  double Vd_ = 100.0;
  double Pe_ = 100.0;
  double SVf_ = 60.0;
  double RV_ = 0.0;
  double Pd_ = 8.0;
  for (int fp = 0; fp < 5; ++fp) {
    double SVt_prev = SVt_;
    double dPmv_l = Rmv * SVt_ / Tdias;
    Pd_ = Pmid - dPmv_l;  if (Pd_ < 0.3) Pd_ = 0.3;
    Vd_ = V0d + log1p(Pd_ / Aed) / Bed;
    Vd_ *= (1.0 - 0.08 * AFb);                 // loss of atrial transport
    // f(Pes_l) is strictly DECREASING: both the forward and the regurgitant
    // path drain the ventricle harder as pressure rises.
    double lo = 0.05;
    double hi = 500.0;
    for (int k = 0; k < 60; ++k) {
      double m = 0.5 * (lo + hi);
      double dP = c_s * m - Pmid;  if (dP < 0.5) dP = 0.5;
      double f = EesE * (Vd_ - m / Ea - EROAt * Tr * Kv * sqrt(dP) - V0s) - m;
      if (f > 0.0) lo = m; else hi = m;
    }
    Pe_ = 0.5 * (lo + hi);
    SVf_ = Pe_ / Ea;
    double dP2 = c_s * Pe_ - Pmid;  if (dP2 < 0.5) dP2 = 0.5;
    RV_ = EROAt * Tr * Kv * sqrt(dP2);
    SVt_ = SVf_ + RV_;
    if (fabs(SVt_ - SVt_prev) < 1e-9) break;
  }
  // ---- volume conservation residual: strictly INCREASING in Pla ----------
  double Q_ = SVf_ / T;
  double MAP_ = Rsys * Q_;
  double Ppa_ = Pmid + Rpul * Q_;
  double Ea_rv = c_rv * Rpul / T;
  double eesrv = (Eesrv > 0.02 ? Eesrv : 0.02);
  double Vrv_ = V0s_rv + (SVf_ * (eesrv + Ea_rv) + Pmid) / eesrv;
  double xp = Bed_rv * (Vrv_ - V0drv > 0.0 ? Vrv_ - V0drv : 0.0);
  if (xp > 6.0) xp = 6.0;
  double Psv_ = Aed_rv * (exp(xp) - 1.0);
  if (Psv_ < 0.4)  Psv_ = 0.4;
  if (Psv_ > 40.0) Psv_ = 40.0;
  double Vla_ = V0la + log1p(Pmid / A_la) / Bla;
  double Vb = V0sa + Csa * MAP_ + V0sv + Csv * Psv_ + V0pa + Cpa * Ppa_
              + V0pu + Cpu * Pmid + (Vd_ - 0.5 * SVt_) + Vla_
              + (Vrv_ - 0.5 * SVf_);
  if (Vb - Vtot < 0.0) Pla = Pmid; else Pla_hi = Pmid;
  Ved=Vd_; Pes_l=Pe_; SVf=SVf_; RVol_l=RV_; SVtot_l=SVt_; Ped_l=Pd_;
  Q=Q_; MAP_l=MAP_; Ppa=Ppa_; Psv=Psv_; Vrved=Vrv_;
}
Pla = 0.5 * (Pla + Pla_hi);

// ---- read-outs -----------------------------------------------------------
double Vla   = V0la + log1p(Pla / A_la) / Bla;
double Claop = 1.0 / (Bla * (Pla + A_la));            // DENOMINATOR ONE
double vwave_l = 90.0 * tanh(xi_v * RVol_l / Claop / 90.0);
double Ppcw_l  = Pla + f_v * vwave_l;
double Ves   = Ved - SVtot_l;
double h_ed  = wall_h(Ved, Vwall);
double h_es  = wall_h(Ves > 5.0 ? Ves : 5.0, Vwall);
double sig_ed = Ped_l * sphere_r(Ved) / (2.0 * h_ed);
double sig_es = Pes_l * sphere_r(Ves > 5.0 ? Ves : 5.0) / (2.0 * h_es);
double CI_l = Q * 60.0 / 1000.0 / BSA;

g_HR=HR; g_T=T; g_Tr=Tr; g_Ea=Ea; g_Rsys=Rsys; g_Csv=Csv; g_Cla=Claop;
g_Bla=Bla; g_Bed=Bed; g_EROA=EROAt; g_EROApisa=EROAt*kpisa; g_CD=CD;
g_Areq=Areq_l; g_resv=resv; g_Rmv=Rmv; g_Pla=Pla; g_Ppcw=Ppcw_l; g_vwave=vwave_l;
g_Vla=Vla; g_Ved=Ved; g_Ves=Ves; g_Ped=Ped_l; g_Pes=Pes_l; g_SVf=SVf; g_RVol=RVol_l;
g_SVtot=SVtot_l; g_EF=SVtot_l/(Ved>1.0?Ved:1.0); g_EFfwd=SVf/(Ved>1.0?Ved:1.0);
g_RF=RVol_l/(SVtot_l>1e-6?SVtot_l:1e-6); g_CO=Q*60.0/1000.0; g_CI=CI_l; g_MAP=MAP_l;
g_Ppa=Ppa; g_Psv=Psv; g_Vrved=Vrved; g_sig_ed=sig_ed; g_sig_es=sig_es;
g_h_ed=h_ed; g_dPmv=Rmv*SVtot_l/Tdias; g_LVEDVi=Ved/BSA; g_PVRWU=Rpul/0.06;
g_VTImr = RVol_l / (EROAt > 1e-9 ? EROAt : 1e-9);

// =========================================================================
//  SLOW BLOCK.  The 43 ODEs.
// =========================================================================

// ---- LV structure --------------------------------------------------------
double dsig = sig_ed - SIG_ED_SET;
double gain = (dsig > 0.0) ? k_dil : k_dil * k_dil_rev;
double room = (dsig > 0.0) ? (1.0 - V0d / V0d_max) : 1.0;
if (room < 0.0) room = 0.0;
dxdt_V0d = gain * dsig * room;
if (dV0d_force >= 0.0) dxdt_V0d = (V0d < V0d_stop) ? dV0d_force : 0.0;

// Growth is directional.  Volume overload builds an ECCENTRIC ventricle, mass
// tracking cavity size.  The concentric term is driven by SYSTOLIC stress --
// and mitral regurgitation UNLOADS systole, so that term is NEGATIVE here.
// This is exactly why chronic MR is a dilating, thin-walled, low-fibrosis
// phenotype, and why its hypertrophy never protects the ventricle.
double conc = sig_es / SIG_ES_SET - 1.0;  if (conc < 0.0) conc = 0.0;
double M_target = LVmass_h * pow(Ved / VED_SET, 0.90) * (1.0 + 0.55 * conc);
dxdt_LVmass = k_hyp * (M_target - LVmass);

double fibd = 1.0 + e_ang_fib * (Ang - 1.0) + e_ald_fib * (Ald - 1.0)
              + e_str_fib * (sig_es / SIG_ES_SET - 1.0);
if (fibd < 0.0) fibd = 0.0;
fibd *= (1.0 - EFFECT_MRA_FIB) * (1.0 - EFFECT_SAC_FIB)
        * (1.0 - EFFECT_SG_FIB);
dxdt_Fib = k_fib * fibd * (Fib_max - Fib) - KFIBDEG * Fib;

double dmg = (sig_es / SIG_ES_SET - 1.0 > 0.0 ? sig_es / SIG_ES_SET - 1.0 : 0.0)
             + 0.85 * (Fib / Fib_h - 1.0 > 0.0 ? Fib / Fib_h - 1.0 : 0.0)
             + 0.60 * (NE - 1.0 > 0.0 ? NE - 1.0 : 0.0);
double bbf = BBDUR_BB / 60.0;  if (bbf > 1.0) bbf = 1.0;
double EFFECT_BB_REC = EMAX_BB_REC * bbf * OCC_BB;
double rec = EFFECT_BB_REC
             + 0.55 * (1.0 - sig_es / SIG_ES_SET > 0.0 ? 1.0 - sig_es / SIG_ES_SET : 0.0);
double flr = (Ees - Ees_min) / (Ees0 - Ees_min);  if (flr < 0.0) flr = 0.0;
dxdt_Ees = -k_eesloss * dmg * Ees * flr + k_eesrec * rec * (Ees0 - Ees);
dxdt_BBDUR_BB = OCC_BB;

// ---- LA ------------------------------------------------------------------
double la_ex = Ppcw_l - PPCW_SET - la_dead;
if (la_ex > 0.0) {
  double lr = 1.0 - V0la / V0la_max;  if (lr < 0.0) lr = 0.0;
  dxdt_V0la = k_la * la_ex * lr;
} else {
  // reverse remodelling returns the atrium toward its OWN size and no further
  double def = (V0la - V0la_h) / V0la_h;  if (def < 0.0) def = 0.0;
  dxdt_V0la = -0.40 * k_la * (-la_ex) * def;
}
double flad = 1.0 + 0.90 * (Ppcw_l / PPCW_SET - 1.0) + 0.60 * (Ald - 1.0);
if (flad < 0.0) flad = 0.0;
flad *= (1.0 - EFFECT_MRA_FIB);
dxdt_Fibla = k_fibla * flad * (Fib_max - Fibla) - KFIBLADEG * Fibla;

double afs = (Vla / 110.0 - 1.0 > 0.0 ? Vla / 110.0 - 1.0 : 0.0)
             + 2.2 * ((Fibla - Fib_h) / (Fib_max - Fib_h) - 0.25 > 0.0
                      ? (Fibla - Fib_h) / (Fib_max - Fib_h) - 0.25 : 0.0);
dxdt_AFb = k_af * afs * (1.0 - AFb) - k_af_rev * AFb;

// ---- annulus and leaflets ------------------------------------------------
double annd = (V0la / V0la_h - 1.0 > 0.0 ? V0la / V0la_h - 1.0 : 0.0)
              + 0.85 * (V0d / V0d_h - 1.0 > 0.0 ? V0d / V0d_h - 1.0 : 0.0)
              + k_ann_af * AFb;
if (RingF > 0.5) {
  dxdt_Ann = 0.0;                      // a ring clamps the annulus
} else {
  double ar = 1.0 - Ann / Ann_max;  if (ar < 0.0) ar = 0.0;
  dxdt_Ann = k_ann * annd * ar - k_ann_rev * (Ann - Ann_h);
}

// Leaflet plasticity.  Adaptation is RATE-LIMITED and PARTIAL: measured leaflet
// growth recovers only phi_leaf of the excess coaptation demand.  The residual
// deficit at adaptive equilibrium is (1 - phi_leaf) of the excess; the deficit
// BEFORE adaptation has had time to occur is the whole of it, which is why the
// SPEED of ventricular dilation matters and not only its final size.  The
// target never falls below native area: leaflets do not atrophy.
double demand = Aleaf_h + phi_leaf * (Areq_l - Aleaf_h > 0.0 ? Areq_l - Aleaf_h : 0.0);
double gapl = demand - Aleaf;
double Amax = Aleaf_max_f * Aleaf_h;
double raw = leaf_cap * tanh(k_leaf * gapl / leaf_cap);
if (raw > 0.0) {
  double hr_ = 1.0 - Aleaf / Amax;  if (hr_ < 0.0) hr_ = 0.0;
  raw *= hr_;
} else {
  raw *= 0.30;
}
dxdt_Aleaf = raw;
dxdt_EROApri = k_pri_prog * EROApri;

// ---- pulmonary vasculature and RV ---------------------------------------
double over = Pla - 18.0;  if (over < 0.0) over = 0.0;
dxdt_Rpulfix = f_pvr_fix * k_pvr * over;
double pr = 1.0 - Rpul / Rpul_max;  if (pr < 0.0) pr = 0.0;
double pexc = Rpul - Rpul0 - Rpulfix;  if (pexc < 0.0) pexc = 0.0;
dxdt_Rpul = k_pvr * over * pr - k_pvr_rev * pexc;

double rvl = Ppa / PPA_SET - 1.0;  if (rvl < 0.0) rvl = 0.0;
double rvf = (Eesrv - 0.10) / (Ees_rv0 - 0.10);  if (rvf < 0.0) rvf = 0.0;
double rvu = 1.0 - Ppa / PPA_SET;  if (rvu < 0.0) rvu = 0.0;
dxdt_Eesrv = -k_rvloss * rvl * Eesrv * rvf + k_rvrec * rvu * (Ees_rv0 - Eesrv);
double rvd = (V0drv - V0d_rv_h) / V0d_rv_h;  if (rvd < 0.0) rvd = 0.0;
dxdt_V0drv = k_rvdil * ((Psv - PSV_SET > 0.0 ? Psv - PSV_SET : 0.0)
                        - 0.45 * (PSV_SET - Psv > 0.0 ? PSV_SET - Psv : 0.0) * rvd);

// ---- neurohormonal ------------------------------------------------------
double baro = CI_SET / (CI_l > 0.4 ? CI_l : 0.4) - 1.0;
if (baro < 0.0) baro = 0.0;
if (baro > 1.5) baro = 1.5;
dxdt_NE = (1.0 + k_ne_co * baro - NE) / tau_ne;
double EFFECT_FUR_RAAS = C_FUR / EC50_FUR;  if (EFFECT_FUR_RAAS > 1.5) EFFECT_FUR_RAAS = 1.5;
double Angt = (1.0 + k_ang_co * baro) * (1.0 + 0.45 * EFFECT_FUR_RAAS);
dxdt_Ang = (Angt - Ang) / tau_ang;
dxdt_Ald = (1.0 + k_ald_ang * (Ang - 1.0) - Ald) / tau_ald;

// ---- volume, kidney, biomarkers -----------------------------------------
double EFFECT_FUR = EMAX_FUR * pow(C_FUR, GAMMA_FUR)
               / (pow(EC50_FUR * (1.0 + BRAKE_FUR), GAMMA_FUR) + pow(C_FUR, GAMMA_FUR));
double natri = EFFECT_FUR
               * (eGFR / eGFR_h);
natri += 260.0 * EFFECT_SAC_NA + 260.0 * EFFECT_MRA_NA;
// Natriuresis cannot continue into hypovolaemia: as intravascular volume
// approaches its floor the interstitium stops refilling the vasculature.
double Vmin = Vtot_min_f * Vtot_h;
double sat = (Vtot - Vmin) / (Vtot_h - Vmin);
if (sat < 0.0) sat = 0.0;  if (sat > 1.0) sat = 1.0;
natri *= sat;
dxdt_BRAKE_FUR = KBRAKE_FUR * ((EFFECT_FUR_RAAS > 3.0 ? 3.0 : EFFECT_FUR_RAAS) - BRAKE_FUR);
double raasd = Ald - 1.0 + 0.55 * (Ang - 1.0);  if (raasd < 0.0) raasd = 0.0;
double vroom = 1.0 - (Vtot - Vtot_h) / dVtot_max;  if (vroom < 0.0) vroom = 0.0;
double Vtar = Vtot_h * (1.0 - EFFECT_SG_VOL);
dxdt_Vtot = k_na * raasd * vroom - k_na_esc * (Vtot - Vtar) - natri;

double gp = 1.0 - CI_l / CI_SET;  if (gp < 0.0) gp = 0.0;
double gc = Psv - PSV_SET;      if (gc < 0.0) gc = 0.0;
double gfrt = eGFR_h * (1.0 - k_gfr_perf * gp - k_gfr_cong * gc / 30.0);
if (gfrt < 6.0) gfrt = 6.0;
dxdt_eGFR = (gfrt - eGFR) / tau_gfr;

double b1 = sig_ed / SIG_ED_SET - 1.0;  if (b1 < 0.0) b1 = 0.0;
double b2 = Ppcw_l / PPCW_SET - 1.0;      if (b2 < 0.0) b2 = 0.0;
dxdt_BNP = (1.0 + k_bnp * (b1 + 0.9 * b2) - BNP) / tau_bnp;

double lym = Ppcw_l - Pcrit + 4.0;  if (lym < 0.0) lym = 0.0;
dxdt_Pcrit = k_lymph * lym * (1.0 - Pcrit / Pcrit_max)
             - 0.0025 * (Pcrit - Pcrit_h);

// ---- endpoints ----------------------------------------------------------
// The hazard accumulators do not feed back into any other state, so the
// baseline hazards enter strictly LINEARLY.  Clamp the EXPONENT against
// overflow rather than the rate itself: capping the rate would destroy that
// linearity and would silently make different arms look identical.
double argh = a1_hfh * (Ppcw_l - 16.0) / 5.0 + a2_hfh * (2.40 - CI_l) / 0.50
              + c1_ved * (g_LVEDVi - 75.0) / 25.0;
if (argh > 18.0) argh = 18.0;
dxdt_HFH = h0_hfh * exp(argh);
double argd = b1_d * (Ppcw_l - 16.0) / 5.0 + b2_d * (2.40 - CI_l) / 0.50
              + b3_d * (Psv - 8.0) / 5.0 + c1_ved * (g_LVEDVi - 75.0) / 25.0;
if (argd > 18.0) argd = 18.0;
dxdt_CumHz = h0_death * exp(argd);

double ny = 1.0 + 1.45 * (Ppcw_l - 12.0 > 0.0 ? Ppcw_l - 12.0 : 0.0) / 8.0
            + 1.25 * (2.60 - CI_l > 0.0 ? 2.60 - CI_l : 0.0) / 0.60;
if (ny > 4.4) ny = 4.4;
dxdt_NYHAi = (ny - NYHAi) / 3.0;

// ---- drug PK: continuous input rates ------------------------------------
dxdt_GUT_FUR = F_FUR * rate_fur - KA_FUR * GUT_FUR;
dxdt_CENT_FUR = KA_FUR * GUT_FUR - CL_FUR / V1_FUR * CENT_FUR;
dxdt_GUT_SAC = F_SAC * rate_sac - KA_SAC * GUT_SAC;
dxdt_CENT_SAC = KA_SAC * GUT_SAC - CL_SAC / V1_SAC * CENT_SAC;
dxdt_GUT_VAL = F_VAL * rate_val - KA_VAL * GUT_VAL;
dxdt_CENT_VAL = KA_VAL * GUT_VAL - CL_VAL / V1_VAL * CENT_VAL;
dxdt_GUT_BB  = F_BB  * rate_bb  - KA_BB  * GUT_BB;
dxdt_CENT_BB  = KA_BB * GUT_BB - CL_BB / V1_BB * CENT_BB;
dxdt_GUT_MRA = F_MRA * rate_mra - KA_MRA * GUT_MRA;
dxdt_CENT_MRA = KA_MRA * GUT_MRA - CL_MRA / V1_MRA * CENT_MRA;
dxdt_GUT_SG  = F_SG  * rate_sg  - KA_SG  * GUT_SG;
dxdt_CENT_SG  = KA_SG * GUT_SG - CL_SG / V1_SG * CENT_SG;
dxdt_CENT_SNP = KE_SNP * (snp_inf - CENT_SNP);
dxdt_CENT_DOB = KE_DOB * (dob_inf - CENT_DOB);

// device and ring states are event-driven, not dynamic
dxdt_TEERrmv = 0.0;
dxdt_TEERf   = 0.0;
dxdt_RingF   = 0.0;

$TABLE
double HRo    = g_HR;
double EROA   = g_EROA;
double EROApisa = g_EROApisa;
double RVol   = g_RVol;
double RF     = g_RF;
double SVtot  = g_SVtot;
double SVfwd  = g_SVf;
double EDV    = g_Ved;
double ESV    = g_Ves;
double EF     = g_EF;
double EFfwd  = g_EFfwd;
double CO     = g_CO;
double CI     = g_CI;
double MAP    = g_MAP;
double Pes    = g_Pes;
double Ped    = g_Ped;
double LAP    = g_Pla;
double vwave  = g_vwave;
double ClaOp  = g_Cla;
double Ppcw   = g_Ppcw;
double LAvol  = g_Vla;
double PAP    = g_Ppa;
double CVP    = g_Psv;
double PVRWU  = g_PVRWU;
double sigED  = g_sig_ed;
double sigES  = g_sig_es;
double hED    = g_h_ed;
double dPmv   = g_dPmv;
double VTImr  = g_VTImr;
double CoaptResv = g_resv;
double Areq   = g_Areq;
double CDepth = g_CD;
double LVEDVi = g_LVEDVi;
double Oedema = (g_Ppcw > Pcrit) ? 1.0 : 0.0;
double SurvP  = exp(-CumHz);
double PropRV = g_RVol / (g_Ved > 1.0 ? g_Ved : 1.0);

$CAPTURE @annotated
HRo    : Heart rate (1/min)
EROA   : True effective regurgitant orifice area (cm2)
EROApisa : PISA-reported orifice area, inflated for a crescent (cm2)
RVol   : Regurgitant volume (mL)
RF     : Regurgitant fraction (-)
SVtot  : Total stroke volume (mL)
SVfwd  : Forward stroke volume (mL)
EDV    : LV end-diastolic volume (mL)
ESV    : LV end-systolic volume (mL)
EF     : LVEF, computed against a stroke volume containing the leak (-)
EFfwd  : Forward ejection fraction (-)
CO     : Forward cardiac output (L/min)
CI     : Cardiac index (L/min/m2)
MAP    : Mean arterial pressure (mmHg)
Pes    : LV end-systolic pressure (mmHg)
Ped    : LV end-diastolic pressure (mmHg)
LAP    : Mean left atrial pressure (mmHg)
vwave  : Left atrial v wave amplitude (mmHg)
ClaOp  : DENOMINATOR ONE, operating LA compliance (mL/mmHg)
Ppcw   : Effective pulmonary capillary wedge pressure (mmHg)
LAvol  : Left atrial volume (mL)
PAP    : Mean pulmonary artery pressure (mmHg)
CVP    : Central venous pressure (mmHg)
PVRWU  : Pulmonary vascular resistance (Wood units)
sigED  : End-diastolic wall stress (mmHg)
sigES  : End-systolic wall stress (mmHg)
hED    : End-diastolic wall thickness (cm)
dPmv   : Mean diastolic transmitral gradient (mmHg)
VTImr  : Regurgitant velocity-time integral (cm)
CoaptResv : Coaptation reserve, MR exists only when negative (cm2)
Areq   : Coaptation area demanded by the geometry (cm2)
CDepth : Coaptation depth / tenting (cm)
LVEDVi : LV end-diastolic volume index (mL/m2)
Oedema : Alveolar oedema flag (-)
SurvP  : Survival probability (-)
PropRV : Regurgitant volume per unit LVEDV, proportionality (-)

// ---- pluggable-PK interface: one C_<STEM> and named EFFECT_<STEM>[_x] -----
// ---- per compound (FORK_WORKFLOW_GUIDE.md Part 2) -------------------------
C_BB   : Beta-blocker plasma concentration (mg/L)
C_VAL  : Valsartan plasma concentration (mg/L)
C_SAC  : Sacubitrilat plasma concentration (mg/L)
C_MRA  : Canrenone (spironolactone active metabolite) plasma concentration (mg/L)
C_SG   : SGLT2 inhibitor plasma concentration (mg/L)
C_FUR  : Furosemide plasma concentration (mg/L)
C_SNP  : Nitroprusside effect-site concentration (effect units)
C_DOB  : Dobutamine effect-site concentration (effect units)
EFFECT_BB_HR   : Beta-blocker heart-rate-reduction effect (-)
EFFECT_BB_EES  : Beta-blocker acute negative-inotropy effect (-)
EFFECT_BB_REC  : Beta-blocker chronic contractility-recovery effect (-)
EFFECT_VAL_RSYS : Valsartan systemic-resistance-reduction effect (-)
EFFECT_SAC_CSV : Sacubitril venodilation effect (-)
EFFECT_SAC_NA  : Sacubitril natriuresis effect (-)
EFFECT_SAC_FIB : Sacubitril antifibrotic effect (-)
EFFECT_MRA_FIB : MRA antifibrotic effect (-)
EFFECT_MRA_NA  : MRA natriuresis effect (-)
EFFECT_SG_VOL  : SGLT2i volume-setpoint-reduction effect (-)
EFFECT_SG_FIB  : SGLT2i antifibrotic effect (-)
EFFECT_SNP_RSYS : Nitroprusside systemic-resistance-reduction effect (-)
EFFECT_SNP_CSV : Nitroprusside venodilation effect (-)
EFFECT_DOB_EES : Dobutamine inotropy effect (-)
EFFECT_DOB_HR  : Dobutamine chronotropy effect (-)
EFFECT_FUR     : Furosemide natriuresis effect (mL/day)
EFFECT_FUR_RAAS : Furosemide RAAS-activation driver (-, clamped C_FUR/EC50_FUR)

## =============================================================================
##  SCENARIOS
##  ---------------------------------------------------------------------------
##  library(mrgsolve); library(dplyr)
##  mod <- mread("mr_mrgsolve_model.R")
##  GDMT <- list(rate_fur = 40, rate_bb = 100, rate_sac = 194,
##               rate_val = 206, rate_mra = 25, rate_sg = 10)
##
##  1. HEALTHY BASELINE -- an exact fixed point of all 43 states
##     mod %>% mrgsim(end = 3650, delta = 10) %>% plot(EF + CI + Ppcw + EDV ~ time)
##
##  2. ACUTE severe MR: papillary muscle rupture onto a virgin atrium.
##     Denominator one at its smallest -- expect a giant v wave (~36 mmHg), a
##     wedge of ~35 mmHg, cardiac index ~1.8 and an ejection fraction that looks
##     SUPRANORMAL (0.77) while forward EF is 0.34.
##     mod %>% init(EROApri = 0.60) %>% mrgsim(end = 14, delta = 0.1) %>%
##         plot(RVol + vwave + ClaOp + Ppcw + EF + EFfwd ~ time)
##
##  3. The SAME regurgitant volume in a chronically adapted atrium.
##     Run scenario 2's orifice for years first, or start from a dilated atrium:
##     mod %>% init(EROApri = 0.60, V0la = 190, Fibla = 0.055) %>%
##         mrgsim(end = 30, delta = 0.5) %>% plot(vwave + ClaOp + Ppcw ~ time)
##
##  4. MILD MR IS STABLE.  Ten years, no therapy, and the coaptation reserve
##     never crosses zero.
##     mod %>% init(EROApri = 0.10) %>% mrgsim(end = 3650, delta = 10) %>%
##         plot(EROA + RF + EDV + CoaptResv + Ppcw ~ time)
##
##  5. THE VORTEX.  A moderate orifice that eventually opens a functional one on
##     top of itself, once leaflet growth is exhausted.
##     mod %>% init(EROApri = 0.30) %>% mrgsim(end = 3650, delta = 10) %>%
##         plot(EROA + CoaptResv + Aleaf + Areq + EDV + RF ~ time)
##
##  6. DEGENERATIVE PROGRESSION with watchful waiting.
##     mod %>% init(EROApri = 0.06) %>% param(k_pri_prog = 0.00069) %>%
##         mrgsim(end = 4380, delta = 10) %>% plot(EROA + EF + EDV + LAvol ~ time)
##
##  7. THE OPERATIVE THRESHOLD.  Take a compensated severe primary MR at day
##     400, then abolish the orifice with everything else unchanged.  EF falls
##     0.619 -> 0.490 with contractility untouched, while cardiac index RISES
##     2.11 -> 2.80.  Post-operative EF crosses 0.50 when pre-operative EF is
##     0.627 -- which is where the guideline's 60% comes from.
##     pre  <- mod %>% init(EROApri = 0.42) %>% mrgsim(end = 400)
##     # then re-run from the day-400 state with EROApri = 0
##
##  8. COAPT-like patient (disproportionate MR), medical therapy alone.
##     coapt <- mod %>% init(coapt_init) %>% param(GDMT)
##     coapt %>% mrgsim(end = 730, delta = 5) %>% plot(HFH + SurvP + Ppcw + CI ~ time)
##
##  9. COAPT-like patient WITH edge-to-edge repair.  Expect EROA 0.280 -> 0.090,
##     wedge -5.96 mmHg, cardiac index +8.5%, mitral gradient 3.2 mmHg.
##     coapt %>% init(TEERf = 0.68, TEERrmv = 0.030) %>% mrgsim(end = 730, delta = 5)
##
## 10. MITRA-FR-like patient (proportionate MR) with the SAME device.  The
##     identical procedure buys only +2.2% cardiac index because the ventricle,
##     not the valve, is the disease.  See KNOWN MISS 1: the model still
##     predicts more benefit here than the trial found.
##     mfr_init <- list(Ees = 1.9503, V0d = 141.1263, LVmass = 319.87,
##                      Fib = 0.0376, V0la = 33.026, Fibla = 0.0499,
##                      Ann = 8.2989, Aleaf = 12.5540, Vtot = 4907.68)
##     mfr <- mod %>% init(mfr_init) %>% param(GDMT)
##     mfr %>% init(TEERf = 0.68, TEERrmv = 0.030) %>% mrgsim(end = 365, delta = 5)
##
## 11. DEVICE EFFICACY SWEEP.  In COAPT the endpoint improves monotonically with
##     TEERf; in MITRA-FR a perfect valve barely moves it.
##     lapply(c(0.4,0.55,0.68,0.8,0.9,1.0), function(f)
##       mfr %>% init(TEERf = f, TEERrmv = 0.030) %>% mrgsim(end = 365))
##
## 12. TIMING: the same operation at years 0 to 6.  The irreversible component
##     of pulmonary vascular resistance is a ratchet, so an identical procedure
##     buys progressively less.
##     coapt %>% mrgsim(end = 1460) # then intervene from successive states
##
## 13. HEART RATE, both ways.  Slower means MORE regurgitation per beat and LESS
##     per minute; clamp the rate and read off the optimum.
##     lapply(seq(45, 110, 5), function(h)
##       mod %>% init(EROApri = 0.42) %>% param(HR_fix = h) %>% mrgsim(end = 1))
##
## 14. NITROPRUSSIDE in acute MR: afterload reduction is anti-regurgitant with
##     the orifice UNCHANGED, because the leak and the aorta are in parallel.
##     mod %>% init(EROApri = 0.60) %>% param(snp_inf = 4) %>%
##         mrgsim(end = 2, delta = 0.05) %>% plot(RF + SVfwd + RVol + Ppcw ~ time)
##
## 15. SACUBITRIL/VALSARTAN for 180 days, effect decomposed into an
##     instantaneous impedance component and a slow geometric one by re-running
##     the day-180 state with and without drug present.
##
## 16. ANNULOPLASTY RING vs TEER vs REPLACEMENT, three years.  The ring clamps
##     the annulus but leaves the tethering arm of the loop intact, and the
##     regurgitation comes back.
##     ring <- coapt %>% init(RingF = 1, Ann = 5.50) %>% mrgsim(end = 1095)
##     mvr  <- coapt %>% init(TEERf = 0.97, TEERrmv = 0.012) %>% mrgsim(end = 1095)
##
## 17. ATRIAL FUNCTIONAL MR: a normal-sized ventricle, an atrium and annulus
##     dilated by atrial fibrillation.
##     mod %>% init(AFb = 0.9, V0la = 190, Ann = 9.5) %>% mrgsim(end = 1825, delta = 10)
## =============================================================================
