# Chagas disease / chronic Chagas cardiomyopathy — QSP model

**Chagas disease / chronic Chagas cardiomyopathy (CCC)**

| | |
|---|---|
| Mechanistic map | [`chg_qsp_model.dot`](../../../chagas-disease/chg_qsp_model.dot) · [SVG](../../../chagas-disease/chg_qsp_model.svg) · [PNG](../../../chagas-disease/chg_qsp_model.png) — 231 nodes, 397 edges, 30 clusters |
| ODE model | [`chg_mrgsolve_model.R`](../../../chagas-disease/chg_mrgsolve_model.R) — 68 ODEs, 28 scenarios, validation harness |
| Dashboard | [`chg_shiny_app.R`](../../../chagas-disease/chg_shiny_app.R) — 14 tabs |
| References | [`chg_references.md`](../../../chagas-disease/chg_references.md) — 114 PubMed-verified citations |

---

## 1. The question this model exists to answer

Four randomised clinical trials conducted in Chagas disease produced four
results that appear to contradict one another.

| Trial | Result |
|---|---|
| **STOP-CHAGAS** | Posaconazole cleared blood PCR **better** than benznidazole (93.3% vs 89.7% at 30 days). Yet sustained clearance at 180 days was 13.3% vs 86.7%, and at 360 days 16% vs 96%. |
| **CHAGASAZOL** | The same reversal was independently reproduced. |
| **BENDITA** | Benznidazole for **2 weeks** was effectively equivalent to **8 weeks** (83% vs 89%). 150 mg was equivalent to 300 mg, and once-weekly dosing was equivalent to daily dosing. |
| **BENEFIT** | In 2,854 patients with established cardiomyopathy it doubled the PCR-clearance rate (66.2% vs 33.5%), yet clinical outcomes did not change at all (HR 0.93, 95% CI 0.81–1.07). |

These are not contradictions. **Each of the four results is a measurement on a
different time-scale**, and this model builds those three time-scales as
three physically separate objects, so that all four results come out as
**outputs** rather than exceptions.

```
Time-scale 1  killing                units: hours~days      → determines blood qPCR during treatment
Time-scale 2  dormant reservoir      units: weeks~months    → determines relapse and cure
Time-scale 3  myocardial remodelling units: years~decades   → determines clinical events
```

Each trial measured **the wrong time-scale** for its own endpoint.

---

## 2. The model's five structural axes

### AXIS 1 — Static and cidal are two different properties, and qPCR measures the one that does not produce cure

Circulating trypanosomes are not a reservoir but a **flux** — the output of
intracellular amastigote replication. Anything that stops replication empties
the blood.

- **Azoles (CYP51 / ergosterol)**: ergosterol is needed when a *dividing*
  cell builds membrane. Blocking it arrests growth but does not kill. →
  `EMAXST_A = 0.99`, `EMAXCR_A = 0.025/day`
- **Nitroimidazoles (TcNTR-1)**: glyoxal and nitroso reactive metabolites
  form DNA adducts even in cells that are **not dividing**. → `EMAXST_N =
  0.90`, `EMAXCR_N = 2.0/day`

**Emergent result**: posaconazole wins on the endpoint the trial measured
during treatment, and loses on the endpoint that matters afterwards. This was
not coded in — it comes out of the ratio of the two parameters.

### AXIS 2 — Because immune pressure is antigen-dependent, a static drug suppresses the parasite while disarming the host at the same time

```c
double KIMM = KIMM0 * TH1 * (1.0 - STAT);
```

To kill an infected cell, that cell must **present** amastigote peptide, and
for that the amastigote must synthesise and turn over protein — in other
words, it must replicate. A drug that stops replication does not hand the
parasite over to the immune system; it **hides** it. Without the `(1 -
STAT)` factor, in this model simply blocking replication would sterilise the
patient within weeks, and posaconazole would end up curing Chagas disease.
The second, slower arm is that TH1 is antigen-driven, and the antigen
collapses once replication stops.

### AXIS 3 — Cure is a threshold crossing, not an integral of exposure

There exists a small population of non-replicating, dormant amastigotes
(`KDORM`/`KWAKE`, roughly 0.3% of the burden at steady state) that are drug
**tolerant** rather than drug resistant. Sterile cure is modelled as a
stochastic extinction event, `P(cure) = exp(-P_EST × N_min)`. The number of
surviving organisms falls exponentially, and that exponential drops below the
extinction scale **within the first two weeks**, so cure probability
**saturates**.

**Emergent result**: 2 weeks ≈ 4 weeks ≈ 8 weeks. Every additional week adds
nothing to the cure numerator and everything to the two toxicity risks —
which is why the discontinuation rate in BENDITA's 2-week arm was **0%**,
while the 8-week arm's was not.

### AXIS 4 — The injury rate is the sum of a parasite term and autocatalytic terms, and their ratio is a function of time, not of parasite burden

```
dCMYO/dt = -(INJ_BASE + INJ_PAR + INJ_AUT + INJ_ISC + INJ_STR) × CMYO
```

Only `INJ_PAR` involves the parasite. Each of the other three closes its own
loop by itself.

| Term | Loop |
|---|---|
| `INJ_AUT` | Injury → myocardial myosin release → autoantibodies → injury |
| `INJ_ISC` | Endothelin → microvascular constriction → ischaemia → endothelial injury → endothelin |
| `INJ_STR` | Cardiomyocyte loss → wall thinning → wall stress → cardiomyocyte loss |

Both autocatalytic channels enter **squared**. This is not curve-fitting: a
heart with 5% of the final autoantibody titre and 5% of the final
microvascular abnormality is not self-destructing by 5%. Both channels need a
critical mass of damaged tissue to feed on and grow.

At every time point the model computes an **internal counterfactual**,
`PAF(t) = 1 - INJ(parasite=0)/INJ`.

**The emergent result, and the point of this model**: PAF starts at 0.86 and
falls to 0.18 by 25 years. **There is no age term anywhere in the model.**
BENEFIT's null result then becomes a prediction — the trial enrolled patients
after the very quantity it was trying to move had already left.

### AXIS 5 — Patients die of two different things, and the two compete

Arrhythmic death is driven not by the **amount** of fibrosis but by its
**spatial heterogeneity** (`SCARH`) and by sympathetic denervation
(`SYMPD`). Pump-failure death is driven by EF. Because the model does not add
hazards but integrates a genuine cumulative-incidence function, `dCIF_i/dt =
h_i × S`, suppressing one mode of death necessarily unmasks the other.

---

## 3. Validation

Running `Rscript chg_mrgsolve_model.R` prints the table below. **Across 21
comparable anchors, the median |log₁₀ ratio| = 0.028** (a median ratio of
about 1.06).

| anchor | quantity | observed | predicted | ratio |
|---|---|---|---|---|
| A6 | BENDITA placebo, sustained clearance | 0.03 | 0.030 | 1.00 |
| A6 | BENDITA benznidazole 300 mg × 8 weeks | 0.89 | 0.900 | 1.01 |
| A6 | BENDITA benznidazole 300 mg × 4 weeks | 0.89 | 0.897 | 1.01 |
| A6 | BENDITA benznidazole 300 mg × **2 weeks** | 0.83 | 0.835 | 1.01 |
| A6 | BENDITA benznidazole 150 mg × 4 weeks | 0.83 | 0.893 | 1.08 |
| A6 | BENDITA benznidazole 300 mg once weekly × 8 weeks | 0.83 | 0.899 | 1.08 |
| A8 | STOP-CHAGAS benznidazole, 30-day PCR negative | 0.897 | 0.900 | 1.00 |
| A8 | STOP-CHAGAS posaconazole, 30-day PCR negative | 0.933 | 0.930 | 1.00 |
| A9 | STOP-CHAGAS benznidazole, 180-day sustained | 0.867 | 0.900 | 1.04 |
| A9 | STOP-CHAGAS **posaconazole, 180-day sustained** | 0.133 | 0.100 | 0.75 |
| A10 | STOP-CHAGAS posaconazole, 360-day PCR positive | 0.84 | 1.00 | 1.19 |
| A4 | BENEFIT composite endpoint, placebo arm, 5.4 years | 0.291 | 0.299 | 1.03 |
| A4 | BENEFIT composite endpoint, benznidazole arm | 0.275 | 0.293 | 1.07 |
| A4 | BENEFIT hazard ratio | 0.93 | 0.977 | 1.05 |
| A1 | BENEFIT PCR clearance, placebo arm | 0.335 | 0.335 | 1.00 |
| A5 | BENEFIT baseline PCR-positivity rate | 0.605 | 0.665 | 1.10 |
| A13 | LGE mass in established CCC (g) | 15.2 | 11.5 | 0.76 |
| A14 | Annual all-cause mortality, severe CCC | 0.079 | 0.089 | 1.12 |
| A14 | Annual sudden-death mortality | 0.026 | 0.020 | 0.78 |
| A14 | Annual pump-failure mortality | 0.035 | 0.033 | 0.93 |
| A14 | Annual stroke mortality | 0.004 | 0.0023 | 0.57 |

**The most important validation fact**: the two calibration constants
carrying the greatest weight — the nitro class's dormant-form killing rate
(`EMAXCD_N`) and the per-organism re-establishment probability (`P_EST`) —
were fitted **only to BENDITA's duration series**, and then, left untouched,
predicted STOP-CHAGAS. STOP-CHAGAS is a different trial, a different drug
class, and a different endpoint.

### Reported failures (failures, reported not repaired)

1. **Overpredicts the BENDITA 150 mg and once-weekly arms by 8%** (0.89 vs an
   observed 0.83). Cure in the model is determined purely by exposure, and
   the lower adherence actually seen in the low-dose/intermittent-dosing arms
   is not in the model.
2. **Predicts 0.100 for STOP-CHAGAS posaconazole at 180 days (observed
   0.133)**, and predicts a 360-day PCR-positivity rate of 1.00 (observed
   0.84). Both values are the same defect in the same direction — the
   model's azole fails **more completely** than what was observed. The
   observed residual success sits in the reinfection/misclassification band,
   and given that the placebo arm itself was 10%, posaconazole is
   indistinguishable from placebo.
3. **Underpredicts LGE mass by 24%** (11.5 g vs 15.2 g). This model's
   collagen is a single homogeneous fraction, whereas actual LGE visualises
   localised clumps of replacement fibrosis.
4. **Underpredicts stroke mortality by 43%.** Atrial fibrillation does not
   emerge in this model and must be switched on via the `AFIB` parameter;
   the default patient does not have it switched on. This is a genuine
   structural gap.
5. **Predicts a BENEFIT HR of 0.977 (observed point estimate 0.93).** The
   model's prediction lies within the observed 95% CI (0.81–1.07) but sits
   closer to 1 than the observed point estimate does.

### Defects found and fixed during development (defects found during integration)

- **Integrator blow-up**: early builds produced NaN at year 20 of the
  simulation. The cause was an unbounded wall-stress feedback — dilation →
  wall thinning → rising stress → faster dilation was a positive feedback
  with no ceiling. Fixed by saturating the stress drive and capping EDV at
  450 mL.
- **Moving baseline**: when the immune initial values were left as rough
  guesses, the parasite compartment slowly drifted out of immune control, so
  every treatment effect was measured against a moving baseline. Now every
  immune initial value is the **analytical steady state** at the baseline
  antigen level, and `KTH` is set so that `TH1_ss = 1`.
- **A static drug turning curative**: before the `(1 - STAT)` factor was
  introduced, merely blocking replication let standing immune pressure
  sterilise the patient within weeks, so the model was claiming posaconazole
  cured Chagas disease.
- **A deterministic model denying cure**: 0.26 surviving organisms recover to
  300,000 within 7 weeks, so in a deterministic ODE every patient relapses.
  The `STERILE` switch follows the branch of the extinction lottery in which
  the last organism died, and `p_cure()` weights that branch.
- **Coarse sampling wrecking the cure probability**: `p_cure()` reads the
  trajectory's **trough**. Running with `delta = 91.3` misses the trough, so
  the cure probability comes out randomly wrong. Always use `delta = 1` when
  assessing cure — the validation harness does.
- **Premature PAF collapse**: when the autocatalytic channels were entered
  linearly, PAF fell to 0.33 within 2 years of infection, contradicting the
  one thing the paediatric-treatment literature states clearly. Switching to
  squared dependence gave PAF 0.86 → 0.68 (2 years) → 0.18 (25 years).

---

## 4. The single most useful number this model produces

Printed by the final block of `Rscript chg_mrgsolve_model.R`. `f` = the
fraction of patients who are sterilised.

| Time of treatment | PAF | f | Maximum achievable HR | Required sample size (29% event rate) |
|---|---|---|---|---|
| Infection year 2 (age 32, asymptomatic phase) | 0.680 | 1.00 | 0.320 | **83** |
| Infection year 2 (age 32, asymptomatic phase) | 0.680 | 0.20 | 0.864 | 5,064 |
| Infection year 25 (age 55, established CCC) | 0.182 | 1.00 | 0.818 | 2,687 |
| Infection year 25 (age 55, established CCC) | 0.182 | 0.20 | 0.964 | **78,852** |

BENEFIT randomised 2,854 patients and reported an HR of 0.93 (0.81–1.07).
**The model's reading**: at that disease stage, with that sterilisation
fraction, this trial was about 30-fold too small to detect the effect its
own hypothesis predicted — and, enrolling **the same drug, the same
disease**, 20 years earlier, it would have been adequately powered.

### The effect of treatment timing (same drug · same dose · same susceptibility, only the timing differs)

| Scenario | Time of treatment | EF at age 80 (%) | 50-year cumulative mortality |
|---|---|---|---|
| S19 | Infection year 2 | 48.4 | 0.365 |
| S20 | Infection year 10 | 41.7 | 0.535 |
| S21 | Infection year 18 | 37.0 | 0.696 |
| S22 | Infection year 25 (BENEFIT) | 33.9 | 0.804 |
| S23 | No treatment | 25.4 | 0.947 |

### Demonstrating competing risks (20 years from age 55, SUSC = 1)

| | Sudden-death CIF | Pump-failure CIF | Overall mortality |
|---|---|---|---|
| No treatment | 0.276 | 0.674 | 0.999 |
| Amiodarone | 0.213 | 0.735 | 0.999 |
| ICD | **0.091** | **0.847** | 0.999 |

ICD reduces sudden death by 67%, but that death does not disappear — it
**reappears later as pump failure**. A model that simply adds up hazards
cannot produce this table.

---

## 5. List of compartments (68 ODEs)

| Group | Compartments |
|---|---|
| PK (16) | `GUT_B` `CEN_B` `PER_B` (benznidazole) · `GUT_A` `CEN_A` `PER_A` (azoles: posaconazole/ravuconazole) · `GUT_N` `CEN_N` (nifurtimox) · `GUT_F` `CEN_F` (fexinidazole) · `GUT_M` `CEN_M` `PER_M` `DEA` (amiodarone + desethylamiodarone) · `GUT_C` `CEN_C` (carvedilol) |
| Parasite (5) | `PBLD` `PRH` `PRX` `PDH` `PDX` |
| Antigen · antibody (2) | `ANTIG` `ABIG` |
| Immune (7) | `MPHI` `TH1` `TREG` `TNFA` `IL10` `TGFB` `AAB` |
| Myocardium (10) | `CMYO` `MFB` `COL` `SCARH` `APEX` `MVD` `ET1` `EDV` `ESV` `HYPM` |
| Conduction · autonomic (4) | `SYMPD` `PSYMD` `COND` `SANF` |
| Neurohormonal (4) | `ANGII` `ALDO` `NE` `BNP` |
| Digestive (3) | `ENSN` `ESOD` `COLD` |
| Clinical events (8) | `HTOT` `CIFSCD` `CIFHF` `CIFSTK` `CIFCMP` `CIFPPM` `CIFNHF` `CIFVT` |
| Safety (5) | `BZNCUM` `TAUD` `HRASH` `HNEU` `AMITIS` |
| Counterfactual (2) | `DMGCUM` `DMGNOP` |

## 6. Scenarios (28)

| Series | Scenarios |
|---|---|
| Natural history | S01 low susceptibility (no progression) · S02 moderate (BENEFIT phenotype) · S03 high susceptibility (severe CCC) · S04 digestive form |
| BENDITA | S05 placebo · S06 300 mg×8 weeks · S07 300 mg×4 weeks · S08 300 mg×2 weeks · S09 150 mg×4 weeks · S10 300 mg once weekly×8 weeks |
| STOP-CHAGAS | S11 benznidazole · S12 posaconazole · S13 combination · **S14 structural prediction: posaconazole for 12 months** · S15 fosravuconazole |
| Other drugs | S16 nifurtimox · S17 fexinidazole short course · S18 hypothetical sterilising agent (10× dormant-form killing) |
| Treatment timing | S19 infection year 2 · S20 year 10 · S21 year 18 · S22 year 25 (BENEFIT) · S23 placebo control |
| Cardiac treatment | S24 heart-failure therapy alone · S25 amiodarone · S26 ICD · S27 anticoagulation |
| Strain | S28 TcI strain (same benznidazole) |

**S14 is this model's most falsifiable prediction**: even 12 months of
posaconazole leaves cure probability at zero. A drug that cannot reach the
reservoir cannot cure, no matter the duration.

---

## 7. Known limitations

1. **Deterministic per patient.** Cure is reported as a probability, and
   population-level incidence comes from an explicit virtual population over
   the susceptibility covariate `SUSC` — not from ETAs within a single run.
2. **No transmission dynamics.** Inoculum is an initial condition, not a
   force of infection. Reinfection is not modelled, which is one of the
   reasons the model cannot explain the residual PCR-positivity seen in
   trials run in endemic areas.
3. **No drug–drug interactions.** Posaconazole's CYP3A4 inhibition combined
   with amiodarone matters in practice but is not in the model.
4. **Atrial fibrillation does not emerge.** It must be switched on via the
   `AFIB` parameter. This is the direct cause of the underprediction of
   stroke mortality.
5. **The digestive form is coarse** (3 states). It was included only because
   it shares its enteric-neuron-loss mechanism with the cardiac form; its
   megacolon predictions should not be trusted quantitatively.
6. **No reactivation.** Reactivation under immunosuppression (transplant,
   HIV) is present in the map but not in the ODEs.
7. `p_cure()` reads the trough of the trajectory, so it **must be run with
   `delta = 1`**.

---

## 8. How to run

```bash
# Render the mechanistic map
dot -Tsvg chg_qsp_model.dot -o chg_qsp_model.svg
dot -Tpng -Gdpi=150 chg_qsp_model.dot -o chg_qsp_model.png

# Compile the model + validation table + power calculation
Rscript chg_mrgsolve_model.R

# Dashboard (requires shiny, ggplot2, dplyr, tidyr)
Rscript -e 'shiny::runApp("chg_shiny_app.R")'
```

Environment used for validation: R 4.3.3, mrgsolve 2.0.1, Graphviz 2.43.0.

---

## 9. Disclaimer

This is a QSP model for educational and research purposes. It was built from
published literature and clinical trial data but has not been independently
validated or certified, and must not be used for actual clinical
decision-making, prescribing, or regulatory submission. The parameters are
illustrative approximations, and separate fitting and validation against
real patient data are required.
