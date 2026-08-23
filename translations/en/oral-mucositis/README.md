# Oral Mucositis — QSP Model
### Chemotherapy- and Radiotherapy-Induced Oral Mucositis · Quantitative Systems Pharmacology

<a href="om_qsp_model.svg"><img src="om_qsp_model.png" width="820" alt="Oral mucositis QSP mechanistic map"></a>

---

## The one idea

Oral mucositis is not a toxicity that a drug "has". It is **a failure of the balance
sheet of a regenerating epithelium**.

The oral mucosa is a conveyor belt.

```
S  →  P1 → P2 → P3  →  D  →  shed
│      transit-          differentiated
clonogenic  amplifying   barrier
basal cell  compartment  (post-mitotic)
```

An ulcer exists precisely and only **while D is below the threshold at which the
barrier keeps its continuity**.

And every cytotoxic insult — alkylator, antimetabolite, or photon — strikes the
**proliferative end (S, P)** of this belt and leaves the **differentiated end (D)**
untouched. This one asymmetry makes the whole disease.

- The insult is invisible for **one barrier lifetime**. Then the barrier collapses,
  because nothing has arrived to replace what was shed. The latent period is not a
  signalling riddle but a matter of `1/k_shed`.
- The ulcer persists **until the belt has refilled**.

This model writes that belt as 50 differential equations, runs it, and establishes by
calculation where its own claims hold up and where they collapse.

---

## What held and what collapsed (results, including the refutations)

Every number came out of actually running `om_python_reference.py`, and the full text
is in [`om_reference_output.txt`](../../../oral-mucositis/om_reference_output.txt).

### 1. The central claim survived only by half

The original claim was "**the time of onset is owned by the damage term and the
duration by the regeneration term**". The damage term (`sens`) and the regeneration
term (`lamS`) were each swung over a 2-fold range and tested as a 2×2.

| Chemoradiation 70 Gy + cisplatin | Variation in time of onset | Variation in duration |
|---|---|---|
| 2-fold change in the damage term | **193.7%** | 139.8% |
| 2-fold change in the regeneration term | 74.9% | 148.2% |
| **Separation ratio** | **2.59** | **1.06** |

The damage term moves the time of onset 2.6 times more — the first half of the claim
lives. For the duration the two terms move it **almost identically** — the second half
dies. In retrospect that is obvious. The duration is **the time taken to climb back
out of a depth**, and what sets the depth is the damage term. The defensible statement
is the weaker one.

> The time of onset belongs to the damage term.
> The duration is shared, and **only the regeneration term can shorten it without
> reducing the cytotoxic dose**.

Where this weaker statement earns its keep is section 8. It is exactly this
distinction that separates the drugs one can give a patient whose chemotherapy dose is
already fixed from the drugs one cannot.

### 2. Cryotherapy is **not** a delivery problem — it is a temperature problem

The mechanism the guidelines cite is that "ice reduces local blood flow and thereby
lowers drug delivery". The mucosa was put into the model as a physiologically sized
**perfusion-limited compartment** (equilibration time constant about 2 minutes), and a
6-hour cooling was decomposed into its two limbs.

| 6-hour cooling | Mucosal AUC | Duration of severe disease |
|---|---|---|
| No cooling | 0.5192 | 6.00 days |
| Both limbs | 0.5173 | **0.00 days** |
| Blood-flow limb only (Q10 = 1) | 0.5173 | **5.90 days** |
| Temperature limb only (fQ = 1) | 0.5192 | **0.00 days** |

Cutting blood flow by 78% drops the mucosal AUC by only **0.4%** and has almost no
clinical effect (6.00 → 5.90 days). A tissue that comes to equilibrium with plasma in
two minutes cannot be protected by cutting off its blood flow — blood flow changes the
**rate** of equilibration, not the exposure. What carries the whole of the effect is
the **Q10** on the alkylation chemistry and on cell cycle progression.

### 3. Which regimen ice works in is **derived**

Ice acts only while it is in the mouth. Therefore

```
benefit = (1 − f_cryo) × (fraction of the damage AUC falling inside the ice window)
        = (1 − f_cryo) × (1 − 2^(−T_ice / t½))
```

The same **30-minute** cooling was applied to three regimens differing only in the
length of exposure.

| Regimen | Reduction in mucosal AUC | Peak ulcer area (cooled / uncooled) |
|---|---|---|
| 5-FU bolus 425 mg/m² d1-5 (t½ 12 min) | 18.4% | 0.000 / 0.164 |
| Melphalan 200 mg/m² (t½ 1.2 h) | 4.4% | 0.538 / 0.677 |
| 5-FU 96-hour continuous infusion | **0.1%** | 0.961 / 0.962 |

The pattern by which MASCC/ISOO recommend cryotherapy for bolus 5-FU and high-dose
melphalan and do not recommend it for continuous-infusion fluoropyrimidines is
**reproduced** from the overlap of half-life and cooling window alone. No parameter
was adjusted to any cryotherapy trial.

### 4. Palifermin — the prediction has the right direction and too small a magnitude, and one sign came out inverted

Exactly **one** number was spent in stage 2: the 9-day grade 3/4 duration in the
placebo arm of Spielberger 2004. The palifermin arm is a prediction.

| | Placebo | Palifermin | Observed |
|---|---|---|---|
| WHO≥3 duration | 9.00 days | **6.50 days** | 9 → 3 days |
| Peak WHO grade | 4 | **3** | 4 → 3/4 |
| Days of opioid analgesic use | 11.11 | **8.81** | 11 → 7 |
| Time of onset | 11.10 | 11.40 | — |

The direction is right and the time of onset barely moves (just the shape a drug on
the regeneration axis ought to have). But **the magnitude is underpredicted** (−28%
against the observed −67%). That is because `eKGF` was left as a structural parameter
and not fitted to the palifermin arm; had it been fitted, this would no longer be a
prediction. It is reported as a miss.

**And the author's prior expectation was refuted with its sign reversed as well.**
Both the map and the model header explained the label's requirement for a 24-hour
separation as "KGF enlarges the proliferative pool and thereby offers the cytotoxic a
bigger target". For alkylator and photon regimens this is **structurally impossible**
— neither carries a cell-cycle-specific term and the killing is first-order in the
pool, so enlarging the pool only raises the absolute number killed proportionally,
while **the fraction killed does not change**. Once the penalty disappears, all that
remains is the 2-day half-life of the KGF effect pool, and therefore the **closer** to
the cytotoxic it is given, the monotonically **better**.

| Interval, last pre-dose ↔ alkylator | −1 day | 0 days | +1 day | +3 days |
|---|---|---|---|---|
| Duration of severe disease | **3.60 days** | 3.80 | 4.50 | 5.70 |
| Against placebo | −60.0% | −57.8% | −50.0% | −36.7% |

The model **cannot** reproduce the label's separation requirement from mucosal
dynamics alone. It must not be read as evidence supporting that requirement.

**Where the penalty is real is with cell-cycle-active agents.** In bolus 5-FU,
switching the single coupling parameter (`fcyc_KGF`) off to 0 lets palifermin
**abolish the ulcer completely** (duration 4.7 days → **0.00 days**). The coupling
eats almost the whole of the benefit. And a five-day dosing course **cannot be
separated** from a growth factor whose effect pool has a 2-day half-life — this is the
model's explanation for palifermin being used in transplant conditioning and not in
ordinary chemotherapy. Conditioning is an **impulse**, so it can be scheduled; a
multi-day course cannot.

### 5. Fractionation: the mucosal price and the tumour benefit compete on the same time axis

| Schedule | Onset | Duration | Tumour BED (Gy10) |
|---|---|---|---|
| 70 Gy / 35 fx / 7 weeks (standard) | 17.5 | 44.0 | 71.5 |
| 81.6 Gy / 68 fx b.i.d. / 7 weeks | **12.9** | **52.4** | **79.3** |
| 70 Gy / 35 fx / 6 weeks (6 fractions a week) | 12.7 | 46.9 | 76.1 |
| 54 Gy / 36 fx t.i.d. / 12 days (CHART-type) | **6.7** | 28.0 | 62.1 |

A 7-day unplanned interruption recovers 6 days of mucosa (44.0 → 38.0 days,
interruption d21-28) and shaves 4.7 Gy10 off the tumour BED.

### 6. The WHO scale is **not** a low-resolution ulcer-area meter — it is a different instrument

The same intervention (6-hour cooling) was measured on both scales while only the
baseline severity was varied.

| sens | Effect on the area scale | Effect on the WHO scale |
|---|---|---|
| 0.70 | −97.8% | 0.0% |
| 1.00 | −97.0% | −100.0% |
| 1.30 | −93.9% | −100.0% |
| 1.70 | −84.9% | −97.9% |

This section was written to hunt for a **ceiling effect**, and what it found was a
**floor effect**. An ordinal endpoint cannot express anything between "severe" and
"not severe", so an intervention that brings a patient just below the threshold and an
intervention that abolishes the disease are recorded **identically** as −100%. A
continuous endpoint separates the two and preserves the clinically meaningful fact
that the effect attenuates with severity (−97.8% → −84.9%). A trial powered on the
WHO grade is blind to the **magnitude** of the effect it is measuring.

---

## Defects the calculation exposed

This model exposed six real defects; all were fixed and all remain recorded in
comments in the code.

1. **A 70 Gy prescription delivered only 20 Gy.** A 2 Gy fraction is delivered over a
   10-minute window, and an adaptive integrator handed a 60-day interval simply steps
   straight over that window. The delivered dose was governed by **the integrator's
   step placement** rather than by the prescription (mucosal BED 24 against the
   correct 84). → Introduced `breakpoints()`, which makes every event boundary an
   explicit integration-interval boundary.
2. **State was not carried between intervals.** If `t_eval` is given, `sol.y[:,-1]` is
   **the last requested time point**, not the end of the interval. At every short
   window the state froze at its starting value, and a 70 Gy course delivered
   **0 Gy**. → Always request the end of the interval as well.
3. **The tissue became permanently sterile.** Because `dS/dt` is proportional to S,
   **S = 0 is an absorbing state**. The first calibration attempt at hitting the onset
   target annihilated the clonogenic pool, and severe mucositis **persisted for 39
   days and never healed**. → The missing biology was **a quiescent,
   label-retaining stem cell reservoir**. Adding it removed the absorbing state and
   gave the regeneration axis a second time constant (`kact`).
4. **The nesting of the calibration was the wrong way round.** The intent was to fix
   the potency from the time of onset, but under an impulse insult the time of onset is
   almost **flat** in potency (6.70 → 6.45 days for a 33% change in potency). The
   duration carries the whole of the dose-response (5.55 → 8.60 days). → The
   assignment was inverted, as the measured sensitivities dictated. This mistake is
   what exposed defect 3.
5. **The hyperfractionation was not hyperfractionation.** Expressing
   `fx_per_week = 10` as one fraction a day spreads 68 fractions over **68 days**
   rather than 34. Since the overall treatment time is the very heart of an
   altered-fractionation comparison, that arm was answering a different question from
   the one on its label (the tumour BED of the hyperfractionated arm came out
   **lower** than standard). → Added real b.i.d./t.i.d. support.
6. **The population simulation quietly reported 0.** The simulation-horizon argument
   leaked into the schedule generator, every worker died with a `TypeError`, and the
   incidence was reported as **0.000** rather than as an exception. → Pull the
   argument out first, and count and report the failures.

And one further, **conceptual** error: stage 1 fitted the deterministic patient to a
38-day **duration** of severe mucositis, but that number is defined only for patients
in whom severe mucositis actually occurred. The stage 1 patient is therefore **the
median affected patient, not the median enrolled patient**, and a population centred
on that patient has an incidence close to 100% rather than 34–43%.
→ A **location** parameter (`sens_med = 0.628`) was put on the population and fitted
to the radiotherapy-alone incidence.

---

## Structural checks — all pass

| Check | Result |
|---|---|
| TNF → NF-κB positive feedback loop gain | 0.067 – 0.292, **all < 1** (no divergence) |
| Accuracy of prescribed dose delivery | 70 Gy/35 fx → mucosal BED **84.0000** (theoretical 84.0000) |
| Integrator convergence (rtol 1e-7 against 1e-9) | difference in peak ulcer area **4.0e-06** |
| Is the drug-free steady state a genuine attractor | barrier drift over 60 days **1.3e-06**, ulcer area **6.1e-06** |

---

## What was bought and at what price (what was fitted)

**Eight published numbers were spent on eight parameters.** Because the fit is exactly
determined, the agreement in section 1 is **arithmetic, not evidence**. The evidence
begins with section 2.

| Stage | Parameter | Number spent | Value |
|---|---|---|---|
| 1 | `greg` | HDM severe onset 6.5 days | 5.849 |
| 1 | `pot_mel` | HDM severe duration 6.0 days | 23.22 |
| 1 | `rad_pot` | Chemoradiation onset 19.0 days | 0.5051 |
| 1 | `lamS` | Chemoradiation duration 38.0 days | 0.06047 |
| 2 | `cy_equiv` | Spielberger placebo arm duration 9 days | 50.14 mg/m² |
| 3 | `sens_med` | Trotti radiotherapy-alone incidence 34% | 0.6283 |
| 3 | `pot_cis` | Trotti chemoradiation incidence 43% | 15.63 |
| 4 | `pot_5fu` | Bolus 5-FU ulcerative stomatitis 5 days | 1.839 |

**Where the fit does not close, that is reported as it stands**: against a target of
6.5 days, the model's HDM time of onset is 7.50 days. The other three land on target
at 6.00 / 18.90 / 38.20.

`pot_mtx` is **constrained by nothing**. Methotrexate is carried for completeness
(GVHD prophylaxis in allogeneic transplantation) but appears in none of the 17
scenarios, so its potency is a placeholder, and it is marked as one.

A caution on the identifiability of `pot_cis`: the incidence rises from 0.33 to 1.00
over a narrow interval, so this parameter is pinned only to within that band. Any
claim that depends on the **magnitude** of the cisplatin contribution is weak to that
extent.

---

## Files

| File | Contents |
|---|---|
| [`om_qsp_model.dot`](om_qsp_model.dot) | Mechanistic map source — **165 nodes · 18 clusters · 237 edges** |
| [`om_qsp_model.svg`](om_qsp_model.svg) / [`.png`](om_qsp_model.png) | Rendering (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`om_python_reference.py`](../../../oral-mucositis/om_python_reference.py) | **The reference implementation that actually runs** — 50 ODEs, interval-wise integrator |
| [`om_calibrate.py`](../../../oral-mucositis/om_calibrate.py) | Four-stage calibration (nested bisection + parallel virtual population) |
| [`om_analysis.py`](../../../oral-mucositis/om_analysis.py) | The experiments of 12 sections — produces every number above |
| [`om_reference_output.txt`](../../../oral-mucositis/om_reference_output.txt) | The full text of that run's output |
| [`calib.log`](../../../oral-mucositis/calib.log) | Calibration log (including the failed attempts) |
| [`om_calibration.json`](../../../oral-mucositis/om_calibration.json) | The eight fitted values — both R and Python read them from here |
| [`om_mrgsolve_model.R`](om_mrgsolve_model.R) | mrgsolve model — 50 ODEs, 17 scenarios |
| [`om_mrgsolve_template.R`](om_mrgsolve_template.R) · [`mkmrgsolve.py`](../../../oral-mucositis/mkmrgsolve.py) | Generation template and generator (prevents R/Python parameter drift) |
| [`om_shiny_app.R`](om_shiny_app.R) | Shiny dashboard — **15 tabs** |
| [`om_references.md`](om_references.md) | **115 papers**, every one resolved live against PubMed |
| [`mkrefs.py`](mkrefs.py) | The generator for that reference list |

### Running

```bash
python3 om_calibrate.py 120     # calibration (tens of minutes)
python3 om_analysis.py          # all 12 sections → om_reference_output.txt
python3 mkmrgsolve.py           # generate om_mrgsolve_model.R
python3 mkrefs.py               # re-resolve the reference list (network required)
dot -Tsvg om_qsp_model.dot -o om_qsp_model.svg
dot -Tpng -Gdpi=150 om_qsp_model.dot -o om_qsp_model.png
```

```r
library(mrgsolve); library(dplyr)
mod <- mread("om_mrgsolve_model.R")
out <- mod %>% data_set(scen_hdm_cryo(6)) %>% mrgsim(end = 45, delta = 0.05)
plot(out, ULC + WHO + VAS + ANC ~ time)
shiny::runApp("om_shiny_app.R")
```

---

## The 17 scenarios

| # | Scenario | Onset | Duration | Peak area | Pain AUC |
|---|---|---|---|---|---|
| 1 | HDM 200 mg/m², no prophylaxis | 7.50 | 6.00 | 0.677 | 50.9 |
| 2 | HDM 140 mg/m² (dose reduced) | — | 0.00 | 0.163 | 24.5 |
| 3 | HDM 200 + cryotherapy 30 min | 7.90 | 4.40 | 0.538 | 43.6 |
| 4 | HDM 200 + cryotherapy 6 h | — | 0.00 | 0.020 | 4.4 |
| 5 | HDM 200 + palifermin (separated) | 8.90 | 1.80 | 0.385 | 36.3 |
| 6 | HDM 200 + palifermin (concurrent) | — | 0.00 | 0.295 | 30.2 |
| 7 | HDM 200 + photobiomodulation (daily) | — | 0.00 | 0.090 | 10.8 |
| 8 | HDM 200 + glutamine | — | 0.00 | 0.323 | 28.6 |
| 9 | HDM 200 + cryotherapy + palifermin | — | 0.00 | 0.017 | 2.1 |
| 10 | TBI-VP16-Cy conditioning (placebo) | 11.10 | 9.00 | 0.836 | 59.5 |
| 11 | TBI-VP16-Cy + palifermin | 11.40 | 6.50 | 0.675 | 46.3 |
| 12 | Head and neck 70 Gy/35 fx + cisplatin | 17.50 | 44.00 | 0.996 | 267.2 |
| 13 | Head and neck 70 Gy/35 fx, radiotherapy alone | 18.90 | 38.10 | 0.987 | 226.7 |
| 14 | Head and neck hyperfractionated 81.6 Gy/68 fx b.i.d. | 12.90 | 52.40 | 0.997 | 303.6 |
| 15 | Head and neck chemoradiation + benzydamine | 19.20 | 38.90 | 0.992 | 235.7 |
| 16 | 5-FU bolus 425 mg/m² d1-5 | — | 0.00 | 0.164 | 23.8 |
| 17 | 5-FU 96-hour continuous infusion 4000 mg/m² | 5.80 | 13.40 | 0.962 | 83.3 |

("—" = severe (WHO≥3) mucositis not reached. The `0.00` for numbers
2·4·6·7·8·9·16 is **a floor effect of the instrument** and is not a claim of 100%
prevention. See sections 6 and 8.)

---

## Limitations

- **This was written in an environment without an R toolchain.**
  `om_mrgsolve_model.R` and `om_shiny_app.R` mirror the executed Python reference
  implementation equation by equation, but **have not themselves been run**. Treat the
  first `mread()` as a syntax check.
- The palifermin effect size is **underpredicted** (−28% against the observed −67%).
  Fitting `eKGF` to the palifermin arm makes it go away, but then it is no longer a
  prediction.
- Because the deterministic patient sits just above the WHO≥3 threshold (peak area
  0.68, threshold 0.35), **floor effects** are pervasive on the ordinal endpoint.
  Incidence-type claims must be made from the virtual population and not from the
  table.
- The temperature limb of cryotherapy assumes `Q10 = 2.5` and `ΔT = 13 °C`
  **structurally**. Since nothing was fitted to any cryotherapy trial, the pattern in
  section 3 is a prediction, but the absolute size of the effect is directly
  proportional to these two numbers.
- The oral microbiota, pain, and myelosuppression axes are each present only as
  minimal structure. The corresponding sections of `om_references.md` give the basis
  for gauging the price of that simplification.
- The cell-cycle dependence of radiation was not put inside the LQ formalism (the
  killing is first-order in the pool). This is the **structural premise** of the
  conclusion in section 4 that alkylator and photon regimens carry no
  target-enlargement penalty.

---

## Disclaimer

This is a qualitative and semi-quantitative QSP model for educational and research
purposes. It was assembled on the basis of the published literature and clinical trial
data but has not been independently verified or certified, and **must not be used
directly for actual clinical decision-making, prescribing, or regulatory submission.**
