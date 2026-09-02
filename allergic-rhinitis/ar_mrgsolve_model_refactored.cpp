$PROB
  Allergic Rhinitis QSP Model -- PK/PD refactor (pluggable PK)
  IgE/Mast Cell . Th2 Cytokines . Eosinophils . PK/PD
  Cetirizine (CETI) | Fluticasone propionate (FP) | Montelukast (MLKT) | Omalizumab (OMA)
  Renamed to this fork’s naming convention: GUT_<STEM>/CENT_<STEM>/PERI_<STEM>,
  C_<STEM> (exposed concentration), EMAX_<STEM>/EC50_<STEM>/GAMMA_<STEM>/EFFECT_<STEM>
  (Hill interface). See ar_refactor_notes.md for per-compound archetype detail.

$GLOBAL
// C_FP/EFFECT_FP are $GLOBAL macros, not $ODE-scope doubles, following the
// clostridioides-difficile-infection precedent (see cdi_mrgsolve_model_
// refactored.R and cdi_refactor_notes.md): Fluticasone is dosed directly
// into CENT_FP (no depot), so an $ODE-scope double reads a stale PRE-dose
// value on the duplicate report row mrgsolve/qspserver emits at the exact
// instant of a dose -- confirmed live via qspserver mrgsolve_api (see
// ar_refactor_notes.md for the before/after comparison). A macro is
// re-expanded textually wherever referenced (including in $TABLE at ITS
// own execution time), so it always reads the current, already-dosed
// state and eliminates the artifact entirely rather than shrinking it.
// This affects only the diagnostic/reporting duplicate-row value, never
// the actual ODE state trajectory.
#define C_FP      (CENT_FP / V1_FP)
#define EFFECT_FP (EMAX_FP * pow(C_FP, GAMMA_FP) / (pow(EC50_FP, GAMMA_FP) + pow(C_FP, GAMMA_FP)))

$PARAM @annotated
  // --- Allergen ---
  ALLERGEN_SS  : 1.0   : Steady-state nasal allergen load (AU/mL, normalized)
  K_ALLERGEN   : 0.1   : Allergen clearance rate (1/h)

  // --- IgE dynamics ---
  KSY_IGE      : 0.002 : IgE synthesis rate (AU/h, baseline)
  KDEG_IGE     : 0.005 : Free IgE degradation rate (1/h, t1/2~138h)
  K_BIND_MAST  : 0.1   : IgE binding to mast cell FceRI (1/h)
  K_OFF_MAST   : 0.001 : IgE dissociation from FceRI (1/h)
  MAST_TOTAL   : 1.0   : Total mast cell FceRI capacity (AU, normalized)

  // --- Mast cell activation & degranulation ---
  EC50_CROSS   : 0.5   : Allergen EC50 for crosslinking IgE-Ag (AU)
  HILL_CROSS   : 2.0   : Hill coefficient for crosslinking
  KDEG_HIST    : 2.0   : Histamine degradation (1/h, DAO/HNMT)
  KDEG_LT      : 1.0   : CysLT degradation (1/h)
  KHIST_PROD   : 5.0   : Max histamine production per unit mast activation
  KLT_PROD     : 2.0   : Max CysLT production per unit mast activation
  KMAST_REC    : 0.05  : Mast cell recharging rate after degranulation (1/h)

  // --- Th2 cytokines (IL-4, IL-5, IL-13) ---
  KSY_TH2      : 0.05  : Th2 cell activation rate driven by IgE/allergen
  KDEG_TH2     : 0.1   : Th2 cell decay (1/h)
  KSY_IL4      : 1.0   : IL-4 synthesis per Th2 unit (pg/mL/h)
  KSY_IL5      : 0.8   : IL-5 synthesis per Th2 unit
  KSY_IL13     : 1.2   : IL-13 synthesis per Th2 unit
  KDEG_IL4     : 0.5   : IL-4 clearance (1/h)
  KDEG_IL5     : 0.3   : IL-5 clearance (1/h)
  KDEG_IL13    : 0.4   : IL-13 clearance (1/h)
  TH2_BASE     : 0.5   : Baseline Th2 activity (atopic, normalized)

  // --- Eosinophil dynamics ---
  EOS_BLOOD_0  : 300   : Baseline blood eosinophil (cells/uL)
  KEO_PROD     : 0.05  : Eosinophil production stimulated by IL-5 (cells/uL/h)
  KEO_SURV     : 0.02  : IL-5-mediated eosinophil survival boost (1/h)
  KEO_DEATH    : 0.08  : Eosinophil apoptosis rate (1/h, t1/2~8.7h)
  KEO_MIGRATE  : 0.01  : Blood->nasal tissue migration (1/h) per chemokine unit
  KCHEMOKINE   : 2.0   : Eotaxin/CCL11 driven chemotaxis coefficient
  KEO_TISSUE0  : 10.0  : Baseline nasal tissue eosinophil (cells/uL tissue)
  KEO_TIS_DEATH: 0.05  : Nasal tissue eosinophil apoptosis (1/h)

  // --- Symptom model ---
  HIST_EC50    : 1.0   : Histamine EC50 for symptom generation
  LT_EC50      : 0.5   : CysLT EC50 for congestion
  EOS_EC50     : 50.0  : Nasal eosinophil EC50 for late symptoms
  SNEEZE_MAX   : 3.0   : Maximum sneezing score
  RHINO_MAX    : 3.0   : Maximum rhinorrhea score
  CONG_MAX     : 3.0   : Maximum congestion score
  PRUR_MAX     : 3.0   : Maximum pruritus score

  // =====================================================
  // PK -- Cetirizine (CETI): oral, depot->central, linear elim
  // Archetype 3 minus peripheral. Tmax ~1h, t1/2 ~10h, Vd 0.56 L/kg
  // =====================================================
  KA_CETI      : 0.9   : Cetirizine absorption rate (1/h)
  CL_CETI      : 7.0   : Cetirizine clearance (L/h)
  V1_CETI      : 70.0  : Cetirizine central volume (L)                 // was VD_CETI
  F_CETI       : 0.70  : Oral bioavailability cetirizine
  EC50_CETI    : 15.0  : Cetirizine EC50 for H1R occupancy (ng/mL)      // was H1_IC50_CETI
  GAMMA_CETI   : 1.0   : Hill coefficient, H1R occupancy                // was H1_HILL
  EMAX_CETI    : 1.0   : Hill ceiling [math-implied: original H1_RO ratio saturates at 1]

  // =====================================================
  // PK -- Fluticasone propionate (FP): intranasal local mucosal
  // Bespoke single-compartment (no depot; see refactor notes for the
  // preserved-not-fixed constant zero-order input term). Systemic exposure <2%.
  // =====================================================
  KA_FP        : 1.5   : FP nasal mucosal absorption rate constant (1/h)
  CL_FP        : 2.0   : FP local elimination rate constant from nasal mucosa (1/h)  // was CL_FP_LOCAL
  V1_FP        : 1.0   : FP mucosal volume (L, local compartment)       // was VD_FP_LOCAL
  BIOCONV_FP   : 200.0 : Intranasal dose bioconversion constant (200 ug ~ 450 nM local) // was FP_DOSE_BIOCONV
  EC50_FP      : 0.5   : FP EC50 for GR-mediated cytokine suppression (nM) // was GR_IC50_FP
  GAMMA_FP     : 1.2   : Hill coefficient, GR occupancy                  // was GR_HILL
  EMAX_FP      : 1.0   : Hill ceiling [math-implied: original GR_OCC ratio saturates at 1]

  // =====================================================
  // PK -- Montelukast (MLKT): oral, depot->central, linear elim
  // Archetype 3 minus peripheral. Tmax 3-4h, t1/2 5.5h, BA ~64%
  // =====================================================
  KA_MLKT      : 0.5   : Montelukast absorption rate (1/h)
  CL_MLKT      : 45.0  : Montelukast clearance (L/h)
  V1_MLKT      : 10.0  : Montelukast central volume (L)                 // was VD_MLKT
  F_MLKT       : 0.64  : Oral bioavailability montelukast
  EC50_MLKT    : 2.0   : Montelukast EC50, CysLT1 blockade (ng/mL)       // was CYSLTR1_IC50
  GAMMA_MLKT   : 1.0   : Hill coefficient, CysLT1 blockade [no explicit Hill term in original -> gamma=1]
  EMAX_MLKT    : 1.0   : Hill ceiling [math-implied: original CYSLTR1_INH ratio saturates at 1]

  // =====================================================
  // PK -- Omalizumab (OMA): SC, depot->central->peripheral (Archetype 3)
  // + bespoke TMDD-style IgE capture (binds the disease model’s own
  // shared IGE_FREE state, not an independent compound-exclusive pool)
  // t1/2 ~26d, ka ~0.004/h, Vss ~78 mL/kg
  // =====================================================
  KA_OMA       : 0.004 : Omalizumab SC absorption rate (1/h)
  CL_OMA       : 0.14  : Omalizumab central clearance (mL/h)
  V1_OMA       : 3140  : Omalizumab central volume (mL, ~45 mL/kg x 70kg)  // was VC_OMA
  V2_OMA       : 2360  : Omalizumab peripheral volume (mL)                 // was VP_OMA
  Q_OMA        : 0.40  : Intercompartmental clearance (mL/h)
  KON_OMA      : 1.0   : Omalizumab-IgE association rate                  // was KON_IGE
  KOFF_OMA     : 0.0001: Omalizumab-IgE dissociation (KD ~0.1 nM)         // was KOFF_IGE

$INIT @annotated
  // Allergen
  AG       : 0.0   : Nasal allergen (AU)

  // IgE compartments
  IGE_FREE : 50.0  : Free IgE (IU/mL equivalent)
  IGE_MAST : 0.5   : Mast cell-bound IgE (normalized)

  // Mast cell state
  MAST_ACT : 0.0   : Mast cell activation level (0-1)
  MAST_CHG : 1.0   : Mast cell granule charge (0-1)

  // Mediators
  HISTAMINE : 0.0  : Nasal histamine (normalized AU)
  CYS_LT   : 0.0   : Cysteinyl leukotrienes (AU)

  // Th2 / cytokines
  TH2      : 0.5   : Th2 cell activity (normalized, atopic baseline)
  IL4      : 5.0   : IL-4 (pg/mL)
  IL5      : 3.0   : IL-5 (pg/mL)
  IL13     : 8.0   : IL-13 (pg/mL)

  // Eosinophils
  EOS_B    : 300.0 : Blood eosinophil (cells/uL)
  EOS_N    : 10.0  : Nasal tissue eosinophil (cells/uL tissue)

  // Cetirizine PK (depot + central)
  GUT_CETI  : 0.0  : Cetirizine depot (mg)                 // was CETI_D
  CENT_CETI : 0.0  : Cetirizine central (mg)                // was CETI_C

  // Fluticasone PK (nasal local, single compartment)
  CENT_FP  : 0.0   : Fluticasone local nasal conc (nM-normalized) // was FP_LOC

  // Montelukast PK
  GUT_MLKT  : 0.0  : Montelukast depot (mg)                 // was MLKT_D
  CENT_MLKT : 0.0  : Montelukast central (mg)                // was MLKT_C

  // Omalizumab PK (depot + 2-comp + IgE complex)
  GUT_OMA   : 0.0  : Omalizumab SC depot (mg)                // was OMA_D
  CENT_OMA  : 0.0  : Omalizumab central (mg/mL x volume -> mg) // was OMA_C
  PERI_OMA  : 0.0  : Omalizumab peripheral (mg)               // was OMA_P
  COMPLEX_OMA : 0.0 : Omalizumab-IgE complex (IU-equivalents) // was OMA_IGE

$MAIN
  // Allergen challenge: set AG from input or steady-state
  // Dose event (compartment 1) delivers allergen pulse

$ODE
  // ==============================================================
  // 1. Allergen kinetics
  // ==============================================================
  double Ag = AG;
  dxdt_AG = -K_ALLERGEN * AG;  // cleared by mucociliary; replenished via dose event

  // ==============================================================
  // 2. IgE dynamics + Omalizumab (OMA) TMDD-style capture
  //    (bespoke: OMA binds the disease model’s own shared IGE_FREE
  //    state, not an independent compound-exclusive receptor pool --
  //    see refactor notes)
  // ==============================================================
  double degen_free = KDEG_IGE * IGE_FREE;
  double bind_rate  = K_BIND_MAST * IGE_FREE * (MAST_TOTAL - IGE_MAST);
  double unbind_rate= K_OFF_MAST  * IGE_MAST;

  // Omalizumab: exposed concentration and IgE-capture term
  double C_OMA = CENT_OMA / V1_OMA;  // mg/mL -- exposed concentration (canonical site)
  double ige_capture = KON_OMA * C_OMA * IGE_FREE - KOFF_OMA * COMPLEX_OMA;

  dxdt_IGE_FREE = KSY_IGE * 200.0 - degen_free - bind_rate + unbind_rate - ige_capture;
  dxdt_IGE_MAST = bind_rate - unbind_rate;
  dxdt_COMPLEX_OMA = ige_capture;

  // ==============================================================
  // 3. Mast cell crosslinking & activation
  // ==============================================================
  double crosslink_frac = pow(Ag, HILL_CROSS) / (pow(EC50_CROSS, HILL_CROSS) + pow(Ag, HILL_CROSS));
  double mast_trigger   = crosslink_frac * IGE_MAST * MAST_CHG;  // depends on bound IgE and granule charge

  dxdt_MAST_ACT = mast_trigger - 0.5 * MAST_ACT;
  dxdt_MAST_CHG = KMAST_REC * (1.0 - MAST_CHG) - mast_trigger * MAST_CHG;

  // ==============================================================
  // 4. Histamine & CysLT release/degradation
  // ==============================================================
  dxdt_HISTAMINE = KHIST_PROD * MAST_ACT - KDEG_HIST * HISTAMINE;
  dxdt_CYS_LT   = KLT_PROD   * MAST_ACT - KDEG_LT   * CYS_LT;

  // ==============================================================
  // 5. Th2 / cytokine cascade
  // ==============================================================
  double th2_drive = TH2_BASE + KSY_TH2 * Ag * IGE_FREE / (1.0 + IGE_FREE);
  dxdt_TH2  = th2_drive - KDEG_TH2 * TH2;
  dxdt_IL4  = KSY_IL4  * TH2 - KDEG_IL4  * IL4;
  dxdt_IL5  = KSY_IL5  * TH2 - KDEG_IL5  * IL5;
  dxdt_IL13 = KSY_IL13 * TH2 - KDEG_IL13 * IL13;

  // ==============================================================
  // 6. Eosinophil dynamics
  // ==============================================================
  double eos_prod    = KEO_PROD * IL5 + KEO_SURV * EOS_B * IL5 / (1.0 + IL5);
  double eos_death_b = KEO_DEATH * EOS_B;
  double chemokine   = KCHEMOKINE * IL5 * IL13;
  double eos_migrate = KEO_MIGRATE * chemokine * EOS_B;
  double eos_death_n = KEO_TIS_DEATH * EOS_N;

  dxdt_EOS_B = eos_prod - eos_death_b - eos_migrate;
  dxdt_EOS_N = eos_migrate - eos_death_n;

  // ==============================================================
  // 7. Drug PK
  // ==============================================================

  // 7a. Cetirizine (CETI) -- Archetype 3 minus peripheral
  dxdt_GUT_CETI  = -KA_CETI * GUT_CETI;
  dxdt_CENT_CETI =  KA_CETI * F_CETI * GUT_CETI - CL_CETI / V1_CETI * CENT_CETI;

  // 7b. Fluticasone propionate (FP) -- bespoke, single compartment.
  // Preserved verbatim from the original, INCLUDING its constant
  // dose-independent "zero-order input" term (KA_FP*BIOCONV_FP/V1_FP is a
  // fixed number that does not depend on any depot state) -- this is a
  // pre-existing behavioral quirk of the original, not introduced by this
  // refactor; disclosed, not fixed, in ar_refactor_notes.md and
  // UPSTREAM_ISSUES.md.
  dxdt_CENT_FP = KA_FP * BIOCONV_FP / V1_FP - CL_FP * CENT_FP;
  // CENT_FP additionally replenished directly by dose events (no separate depot)

  // 7c. Montelukast (MLKT) -- Archetype 3 minus peripheral
  dxdt_GUT_MLKT  = -KA_MLKT * GUT_MLKT;
  dxdt_CENT_MLKT =  KA_MLKT * F_MLKT * GUT_MLKT - CL_MLKT / V1_MLKT * CENT_MLKT;

  // 7d. Omalizumab (OMA) -- Archetype 3 (depot->central->peripheral); IgE
  // capture (ige_capture, C_OMA) computed above in block 2, matching where
  // the original computed it (needed live, every substep, by dxdt_IGE_FREE).
  dxdt_GUT_OMA  = -KA_OMA * GUT_OMA;
  dxdt_CENT_OMA =  KA_OMA * GUT_OMA - (CL_OMA/V1_OMA + Q_OMA/V1_OMA) * CENT_OMA + Q_OMA/V2_OMA * PERI_OMA;
  dxdt_PERI_OMA =  Q_OMA/V1_OMA * CENT_OMA - Q_OMA/V2_OMA * PERI_OMA;

  // ==============================================================
  // 8. Drug pharmacodynamic effects (inhibitory) that feed back into
  //    the live ODE system (must stay in $ODE, matching the original --
  //    see "keep a calculation in the block the original used it in").
  //    Cetirizine’s and Montelukast’s Hill terms do NOT feed back into
  //    any dxdt_ here (same as the original, where the analogous
  //    ODE-block h1ro/cysltr1_inh locals were computed but never used by
  //    any dxdt_ line) -- both are single-sourced instead, computed once
  //    in $TABLE where their only real consumers (the symptom-score
  //    outputs) already lived. See ar_refactor_notes.md.
  // ==============================================================

  // Fluticasone: GR occupancy -> suppress cytokines (IL-4, IL-5, IL-13) & eosinophil.
  // C_FP/EFFECT_FP are $GLOBAL macros (see top of file), not doubles declared
  // here -- avoids the dose-instant stale-value reporting artifact while
  // still being available at $ODE’s live-integration time exactly as before
  // (nM; V1_FP = 1, so numerically identical to CENT_FP directly, matching
  // the original’s undivided use).

  // Apply FP effect on cytokines (transrepression)
  dxdt_IL4  -= EFFECT_FP * 0.6 * IL4;
  dxdt_IL5  -= EFFECT_FP * 0.7 * IL5;
  dxdt_IL13 -= EFFECT_FP * 0.6 * IL13;
  dxdt_EOS_N -= EFFECT_FP * 0.5 * EOS_N;

$TABLE
  // ==============================================================
  // Secondary PK variables -- exposed concentrations (canonical single site
  // per compound) and reporting-only unit conversions
  // ==============================================================
  double C_CETI      = CENT_CETI / V1_CETI * 1000.0;    // ng/mL (mg/L -> ng/mL)
  double C_MLKT       = CENT_MLKT / V1_MLKT * 1000.0;    // ng/mL
  double OMA_CP_REPORT = CENT_OMA / (V1_OMA/1000.0);     // ug/mL (reporting only; was OMA_CP)

  // Cetirizine: H1R occupancy Hill interface. Rename-only from the
  // original’s plain ratio (H1_RO = C^gamma/(EC50^gamma+C^gamma)*100);
  // EMAX_CETI/GAMMA_CETI/EC50_CETI carry the original’s values unchanged.
  double EFFECT_CETI = EMAX_CETI * pow(C_CETI, GAMMA_CETI) / (pow(EC50_CETI, GAMMA_CETI) + pow(C_CETI, GAMMA_CETI));

  // Montelukast: CysLT1 blockade Hill interface. Rename-only from the
  // original’s plain ratio (CYSLTR1_INH = C/(EC50+C)*100, implicit gamma=1).
  double EFFECT_MLKT = EMAX_MLKT * pow(C_MLKT, GAMMA_MLKT) / (pow(EC50_MLKT, GAMMA_MLKT) + pow(C_MLKT, GAMMA_MLKT));

  // Percentage-scale reporting variables (same names/values the original exposed)
  double H1_RO        = EFFECT_CETI * 100.0;   // %
  double GR_OCC_FP     = EFFECT_FP    * 100.0;   // %
  double CYSLTR1_INH  = EFFECT_MLKT  * 100.0;   // %
  double FP_LOCAL_NM  = C_FP;                    // nM (unchanged name/value)
  double CETI_CP      = C_CETI;                  // ng/mL (unchanged name/value)
  double MLKT_CP       = C_MLKT;                  // ng/mL (unchanged name/value)
  double OMA_CP        = OMA_CP_REPORT;           // ug/mL (unchanged name/value)

  // Free IgE reduction (omalizumab; vs baseline 50) -- descriptive-only
  // readout, not an algebraic Hill gate consumed elsewhere: omalizumab’s
  // actual PD action is already fully expressed by the TMDD-style
  // ige_capture ODE term above (see refactor notes).
  double IGE_REDUCTION_PCT = (1.0 - IGE_FREE/50.0) * 100.0;
  // Fractional bound-IgE occupancy, the TMDD-archetype analogue of
  // COMPLEX_<STEM>/RTOT_<STEM>: omalizumab has no fixed RTOT (total IgE is
  // itself a dynamic disease state, not a compound-exclusive constant
  // pool), so this uses the instantaneous free+bound total in its place.
  double EFFECT_OMA = COMPLEX_OMA / (IGE_FREE + COMPLEX_OMA + 1.0e-9);

  // ==============================================================
  // Drug-modified mediator outputs
  // ==============================================================
  double H1_EFF_HIST = HISTAMINE * (1.0 - EFFECT_CETI);      // Histamine x (1-H1RO)
  double LT_EFF      = CYS_LT   * (1.0 - EFFECT_MLKT);       // LT x (1-CysLT1 block)

  // ==============================================================
  // Symptom scores (0-3 scale, Emax model)
  // ==============================================================
  // Sneezing: mainly histamine + sensory nerve
  double SNEEZE = SNEEZE_MAX * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST);

  // Rhinorrhea: histamine + LTs + glandular
  double RHINORRHEA = RHINO_MAX * (0.6 * H1_EFF_HIST + 0.4 * LT_EFF) /
                      (HIST_EC50 + 0.6 * H1_EFF_HIST + 0.4 * LT_EFF);

  // Congestion: mainly LTs + PGD2 (approximated by CYS_LT) + eosinophil
  double EOS_CONG = EOS_N / (EOS_EC50 + EOS_N);
  double CONGESTION = CONG_MAX * (0.5 * LT_EFF / (LT_EC50 + LT_EFF) + 0.3 * EOS_CONG + 0.2);

  // Pruritus: histamine + Th2 cytokines (IL-31 proxy -> IL-13)
  double IL13_NORM = IL13 / (IL13 + 8.0);
  double PRURITUS = PRUR_MAX * (0.7 * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST) + 0.3 * IL13_NORM);

  // TNSS (0-12)
  double TNSS = SNEEZE + RHINORRHEA + CONGESTION + PRURITUS;

  // ==============================================================
  // Biomarker endpoints
  // ==============================================================
  double FREE_IGE_IU  = IGE_FREE;          // IU/mL
  double BLOOD_EOS_UL = EOS_B;             // cells/uL
  double NASAL_EOS_UL = EOS_N;             // cells/uL tissue
  double TRYPTASE_MCG = MAST_ACT * 15.0;  // serum tryptase proxy (ng/mL)

$CAPTURE @annotated
  CETI_CP     : Cetirizine plasma concentration, unchanged-name reporting alias of C_CETI (ng/mL)
  MLKT_CP      : Montelukast plasma concentration, unchanged-name reporting alias of C_MLKT (ng/mL)
  OMA_CP       : Omalizumab plasma concentration, unchanged-name reporting alias (ug/mL)
  FP_LOCAL_NM  : Fluticasone local nasal concentration, unchanged-name reporting alias of C_FP (nM)
  H1_RO        : H1 receptor occupancy, unchanged-name reporting alias of EFFECT_CETI*100 (%)
  GR_OCC_FP     : GR occupancy, unchanged-name reporting alias of EFFECT_FP*100 (%)
  CYSLTR1_INH  : CysLT1 receptor inhibition, unchanged-name reporting alias of EFFECT_MLKT*100 (%)
  SNEEZE       : Sneezing score (0-3)
  RHINORRHEA   : Rhinorrhea score (0-3)
  CONGESTION   : Congestion score (0-3)
  PRURITUS     : Pruritus score (0-3)
  TNSS         : Total Nasal Symptom Score (0-12)
  FREE_IGE_IU  : Free IgE (IU/mL)
  BLOOD_EOS_UL : Blood eosinophil (cells/uL)
  NASAL_EOS_UL : Nasal tissue eosinophil (cells/uL tissue)
  TRYPTASE_MCG : Serum tryptase proxy (ng/mL)
  H1_EFF_HIST  : Histamine effective after H1 blockade (normalized AU)
  LT_EFF       : CysLT effective after CysLT1 blockade (AU)
  IGE_REDUCTION_PCT : Free IgE reduction vs baseline (%)
  C_CETI       : Cetirizine exposed concentration, canonical single site (ng/mL)
  C_FP         : Fluticasone exposed concentration, canonical single site (nM)
  C_MLKT       : Montelukast exposed concentration, canonical single site (ng/mL)
  C_OMA        : Omalizumab exposed concentration, canonical single site (mg/mL)
  EFFECT_CETI  : Cetirizine Hill effect, H1R occupancy fraction (0-1)
  EFFECT_FP    : Fluticasone Hill effect, GR occupancy fraction (0-1)
  EFFECT_MLKT  : Montelukast Hill effect, CysLT1 blockade fraction (0-1)
  EFFECT_OMA   : Omalizumab bespoke effect, fractional bound-IgE occupancy (0-1)
  // COMPLEX_OMA is a $CMT compartment -- exposed as an output column
  // automatically; not re-listed here (mrgsolve 2.0.1 rejects a
  // compartment name appearing in $CAPTURE, same defect class as the
  // original’s bare-compartment capture list -- see refactor notes).

