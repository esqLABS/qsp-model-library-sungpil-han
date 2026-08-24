# Distal Renal Tubular Acidosis (dRTA) — QSP Model

> A 56-ODE model that writes the H⁺ pump of the α-intercalated cell as a
> **saturating actuator**, distributes the acid load across **three sinks with
> different time constants**, and describes alkali therapy as a problem of
> **delivery rate**.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map | [`drta_qsp_model.dot`](drta_qsp_model.dot) · [SVG](drta_qsp_model.svg) · [PNG](drta_qsp_model.png) — **171 nodes / 267 edges / 19 clusters** |
| ⚙️ mrgsolve model | [`drta_mrgsolve_model.R`](drta_mrgsolve_model.R) — **56 ODEs / 28 scenarios** |
| 📊 Shiny dashboard | [`drta_shiny_app.R`](drta_shiny_app.R) — **10 tabs** |
| 📚 References | [`drta_references.md`](drta_references.md) — **463 items / 24 sections** (every one looked up through the NCBI E-utilities) |
| 🧪 Validation log | [`drta_validation_log.md`](drta_validation_log.md) — the raw tool output for every number below |
| 🐍 Reference implementation | [`drta_python_reference.py`](drta_python_reference.py) · [`drta_python_report.py`](drta_python_report.py) — a dependency-free RK4 harness that actually integrates and calibrates the 56 equations |

<a href="drta_qsp_model.svg"><img src="drta_qsp_model.png" width="100%" alt="dRTA QSP map"></a>

---

## 1. What this model does differently

dRTA is usually described as "a disease in which the distal tubule cannot excrete
acid, so the blood bicarbonate falls". Yet there is one thing in the clinical data
that this account does not explain. It is a pair of numbers measured at the same
visit in the 6-year follow-up cohort of ADV7103 (Sibnayal®) (B22CS, n=30).

| Measurement | Value | Interpretation |
|---|---|---|
| plasma HCO₃⁻ | **22.0 ± 3.2 mmol/L** | "very nearly normal" |
| lumbar spine BMD z-score | **−1.1 ± 1.0** | "clearly abnormal" |

The same patient, the same visit, the same lesion. This model treats the
discrepancy not as noise or measurement error but as **a structural fact the
model is obliged to reproduce**.

The core structure is this. The renal acid excretion control system is written as
an **integral controller plus a saturating actuator**. The error
e = (HCO₃⁻ set point − measured value) drives V-ATPase trafficking (the integral
term, τ≈days) and cellular pH sensing (the proportional term, τ≈hours), and the
resulting output VH is **clamped to [0, 1]**. dRTA is not "a disease in which the
pathway has gone" but **a state in which this loop is saturated**.

Once the actuator is against the rail, the remaining acid load has to keep being
taken up somewhere. The model divides that sink into three, each with a different
time constant.

| | Sink | Time constant | Measured? |
|---|---|---|---|
| SINK 1 | ECF bicarbonate | minutes | **yes — this is the only one that is measured** |
| SINK 2 | intracellular/protein non-bicarbonate buffer system | hours | no (net flux at steady state = 0) |
| SINK 3 | **the carbonate of bone** | **years** | no |

At steady state the net flux through SINK 2 is 0, so the unmet acid load is taken
up entirely by SINK 3. The plasma HCO₃⁻ then becomes **a value set not by the size
of the lesion but by the dose-response curve of bone**. That is, **plasma
bicarbonate is not a flux but a ratio.** The 22.0 and the −1.1 above are one acid
load seen through two windows with different gear ratios, and the model produces
both at once.

The bone term is **rectified** in pH — `max(0, pH_ref − pH)`. Rectification is
convex, so by Jensen's inequality **an HCO₃⁻ profile with large fluctuations
consumes more bone even at the same mean.** It is by this route that the dosing
schedule reaches bone mineral density.

The same logic applies to potassium. Plasma K⁺ is also a ratio — about 1 mmol/L
per 300 mmol of total body deficit. And since distal Na⁺ reabsorption, if it is
not electrically offset by H⁺ secretion, must be offset by K⁺ secretion, **a
lesion of the H⁺ pump is in itself a K⁺-wasting lesion**. This is the mechanistic
reason why sodium-salt alkali worsens the very hypokalaemia it is meant to treat,
and why alkali in dRTA has to be a potassium salt.

**The urine pH is not a fitted state variable.** It is solved by bisection at
every derivative evaluation as the root of the proton balance within the tubule,
`J_H(pH) = delivered HCO₃⁻ + titratable acid + trapped NH₄⁺ + protonated
citrate`. NH₄⁺ trapping, titratable acid, bicarbonaturia and brushite
supersaturation therefore move together rather than drifting out of step with one
another.

---

## 2. Validation — what the model passes

Independent observations that were not built into the parameters.

| Validation item | Literature value | Model |
|---|---|---|
| normal adult plasma HCO₃⁻ / pH / pCO₂ | 24–26 / 7.40 / 40 | **24.7 / 7.403 / 40.8** |
| normal adult NAE (Western diet) | 43–60 mEq/day | **67 mEq/day** |
| maximal NAE (under acid loading) | 300–450 mEq/day | **254 mEq/day** (dietary acid ×4) |
| normal free-flow urine pH | 5.5–6.5 | **6.32** (adult) · 6.35 (child) · 5.94 (infant) |
| minimum normal urine pH after NH₄Cl loading | **< 5.45** | **4.61** ✔ |
| in dRTA urine pH cannot fall below 5.5 | > 5.5 | **6.92** (complete form, child) ✔ |
| FE HCO₃⁻ (distal type) | < 3 % | **< 1 %** ✔ |
| slope of urinary Ca against NAE (Lemann 1999) | **0.035 mmol/mEq** | **0.0357 mmol/mEq** ✔ |
| hypokalaemia in complete dRTA | 2.5–3.5 | **3.10** (child) · 2.60 (severe) ✔ |
| hyperchloraemia in complete dRTA | 108–115 | **111.6** ✔ |
| hypercalciuria in complete dRTA | > 4 mg/kg/day | **4.8 mg/kg/day** ✔ |
| hypocitraturia in complete dRTA | < 1.7 mmol/day | **0.48 mmol/day** ✔ |

The brushite supersaturation `SS` is a **relative index**. Its absolute
normalisation has not been rigorously calibrated, so healthy individuals span
0.8–1.6 and do not sit clearly below 1.0. Every conclusion about supersaturation
in this document is therefore a **relative comparison between otherwise identical
conditions**, and such a comparison is invariant to the normalisation.

The Lemann slope was **left as a validation target rather than made a structural
parameter**. Calciuria is described separately along three routes — (i) the
dietary acid load, (ii) the calcium buffered out of bone, and (iii) acid-mediated
inhibition of distal TRPV5 — and the composite slope that comes out of a dietary
acid titration is 0.0357.

### The differential diagnosis is derived from the model

An acute NH₄Cl load of 1.87 mEq/kg in an adult (the 0.1 g/kg protocol):

| Lesion | HCO₃⁻ before loading | Urine pH before loading | Minimum urine pH | Acidification |
|---|---|---|---|---|
| normal | 24.8 | 6.32 | **4.61** | ✔ |
| **incomplete dRTA** | **24.8 (normal!)** | 6.29 | **6.20** | ✘ |
| pure gradient defect | 24.8 (normal) | 6.31 | 5.81 | ✘ |
| pure rate defect | 17.1 | 6.90 | 6.88 | ✘ |
| complete dRTA | 16.8 | 6.91 | 6.89 | ✘ |

**Incomplete dRTA is not assumed, it is derived.** The clinical definition
(plasma bicarbonate normal, acidification fails) falls straight out of the
equations without a parameter of its own.

### And the incomplete form is not a gentle slope but a **cliff**

When the model was being built, the expectation was that incomplete dRTA would be
"a state that quietly gnaws at bone beneath a normal bicarbonate". The simulation
did not say that. On a normal diet the patient with the incomplete form is
**indistinguishable from a healthy individual on every routine test** — it is just
that only 5 % of the actuator reserve is left. Raise the dietary acid load by only
30 % and that 5 % is gone.

| Dietary acid (× normal) | 0.7 | 1.0 | **1.3** | 1.6 | 2.0 |
|---|---|---|---|---|---|
| plasma HCO₃⁻ | 24.68 | 24.67 | **20.79** | 17.21 | 14.50 |
| remaining reserve | 29.7 % | 5.4 % | **0 %** | 0 % | 0 % |
| urinary citrate (mmol/day) | 2.51 | 2.50 | **0.73** | 0.72 | 0.72 |
| bone base flux (mEq/day) | 0.3 | 0.5 | **10.8** | 29.3 | 55.2 |
| *(healthy comparison)* | | 0.5 | | **0.9** | |

Going from a dietary acid of 1.0 to 1.3, the bone flux jumps **22-fold** and the
urinary citrate collapses to a third. Across the same dietary increase a healthy
individual retains 55 % of the reserve and the bone flux stays at 0.9 mEq/day. The
transition is abrupt because this is not a gentle loss of function but
**actuator saturation**.

**Incomplete dRTA is therefore not "a mild disease" but "a severe disease that
looks normal and is one dietary step away".** A testable prediction: raising the
dietary acid load in a patient with stones and incomplete dRTA should change the
urinary citrate and the bone turnover markers **in a step**, and there should be
no such step in a healthy individual.

---

## 3. Two of the author's hypotheses were rejected by the simulation

Two of the hypotheses made when the model was being built were **wrong, and the
simulation caught it.** They are left here rather than deleted.

### Rejection ①: "rate matching" — void as regards plasma bicarbonate

The hypothesis was this. Base arriving faster than the rate of endogenous acid
production exceeds the proximal reabsorption threshold and is thrown away as
bicarbonaturia, so a prolonged-release formulation yields a higher plasma
bicarbonate at the same daily dose.

The release rate was swept over a **16-fold range** (ADV7103 BID, 1.0 mEq/kg/day,
complete form, child):

| Bicarbonate release rate (/h) | 0.10 | 0.215 | 0.40 | 0.70 | 1.10 | 1.60 |
|---|---|---|---|---|---|---|
| plasma HCO₃⁻ (mmol/L) | 23.37 | 23.34 | 23.32 | 23.31 | 23.30 | 23.30 |

**Across the whole range the difference is 0.07 mmol/L — clinically zero.** The
direction comes out as predicted (slower is higher) but there is no magnitude.
The reason is also plain from the calculation: the patient is treated to an HCO₃⁻
of 20–23 and the proximal threshold is about 25.7, so **the threshold is never
exceeded in the first place.** The 30 % of "wasted alkali" that the model computes
is not dose overshoot but **the baseline bicarbonaturia that always leaves in
alkaline urine** — a property of the disease, not of the schedule.

### Rejection ②: "escape from the NaDC1 Tm" — void as regards citraturia too

The hypothesis was this. Since NaDC1 is Tm-limited, giving citrate as a bolus
escapes reabsorption and so increases citraturia further. The two endpoints would
therefore want exactly opposite kinetics.

| Citrate release rate (/h) | 0.10 | 0.20 | 0.40 | 0.70 | 1.10 | 1.60 |
|---|---|---|---|---|---|---|
| urinary citrate (mmol/day) | 1.65 | 1.64 | 1.63 | 1.63 | 1.63 | 1.63 |

**A difference of 0.02 mmol/day — again zero.** Because hepatic first-pass
extraction is 85 %, there is little citrate reaching the circulation to be
filtered in the first place, and that amount cannot push the filtered load above
the Tm.

### And yet the release rate **moves stone risk by 27 %**

In the same sweep, a cell that had been assumed to be nothing moved. The plasma
bicarbonate and the urine pH are **identical to two decimal places** and yet the
brushite supersaturation differs.

| Bicarbonate release rate (/h) | 0.10 | 0.215 | 0.40 | 0.70 | 1.10 | 1.60 |
|---|---|---|---|---|---|---|
| urine pH | 6.72 | 6.72 | 6.72 | 6.72 | 6.72 | 6.72 |
| urinary Ca/citrate | 0.159 | 0.169 | 0.179 | 0.187 | 0.192 | **0.194** |
| **brushite supersaturation** | **2.21** | 2.40 | 2.56 | 2.69 | 2.77 | **2.80** |

Slow release lowers the supersaturation by **27 %.** The mechanism is
rectification. The hypocalciuric action of alkali follows the **instantaneous**
acid load and urinary calcium is rectified at zero, so a bolus suppresses
calciuria for only two hours whereas a prolonged-release form suppresses it all
day. Jensen's inequality operates here in the same way as it does in bone.

**The conclusion is therefore inverted: prolonged-release alkali is not an
acid-base drug but a stone drug.** And this agrees exactly with the pattern B21CS
actually reported — the mean difference in bicarbonate was modest
(non-inferiority p<0.0001, superiority p=0.0008), while what improved greatly was
the **urinary Ca/citrate ratio** (falling below the stone-risk threshold in 56 %
of previous non-responders to SoC). Had the original hypothesis been right, the
bicarbonate should have opened up a wide gap, and it did not.

---

## 4. The 43 % → 90 % is a matter of exposure, not of pharmacokinetics

B21CS (PMID 32712761) reported that switching the same 37 patients from standard
therapy (median three times daily) to ADV7103 twice daily raised the plasma
bicarbonate response rate from **43 % to 90 %**. How much of that is explained by
changing nothing but the schedule at the same daily mEq/kg was decomposed here
(28 virtual patients, fixed at 0.90 mEq/kg/day, with adherence fixed as an
**explicit input**).

| Condition | Response rate | HCO₃⁻ | Trough | Urinary citrate | Supersaturation | Bone flux |
|---|---|---|---|---|---|---|
| SoC IR three times daily, adherence 0.92 | 39.3 % | 20.10 | 19.45 | 1.08 | 4.34 | 10.9 |
| ADV7103 twice daily, adherence 0.92 | 46.4 % | 20.80 | 20.17 | 1.34 | 3.39 | **8.2** |
| **effect of the schedule alone** | **+7.1 pp** | +0.70 | | | −22 % | **−25 %** |

**The pharmacokinetics of the schedule accounts for only 7 of the 47 points.** To
generate the rest, exposure has to move. SoC adherence was swept (the ADV arm
fixed at 0.92):

| SoC adherence | 0.60 | 0.70 | 0.80 | 0.92 |
|---|---|---|---|---|
| SoC response rate | 10.7 % | 17.9 % | 32.1 % | 39.3 % |
| difference vs ADV (0.92) | **+35.7 pp** | +28.6 pp | +14.3 pp | +7.1 pp |

To reproduce the roughly 47 points observed, **the actual SoC adherence would
have to be near 0.6**. This is consistent with what the trial itself reported —
palatability +25 mm VAS, gastrointestinal discomfort −14.2 mm, dosing frequency
3→2 per day, and the accounts in the qualitative research of parents giving up
work because of dosing at school and dosing overnight.

**A testable prediction: standard therapy given under directly observed therapy
(DOT) should perform almost as well as ADV7103 on plasma bicarbonate.** If that
prediction is wrong, this conclusion of the model is wrong.

Meanwhile the bone flux falls by **25 %** although the HCO₃⁻ difference is only
0.70 mmol/L. The convexity of the rectified bone term has amplified a small
difference in means into a large difference in flux, and this is the route by
which the lumbar spine BMD z improved from −1.1 to −0.8 in B22CS.

---

## 5. Model specification

### The 56 state variables

- **drug/gut (11)** — immediate-release bicarbonate and citrate; the
  immediate-release citrate granules and the prolonged-release bicarbonate
  granules of ADV7103 held **separately**; the cation of the salt; plasma
  citrate; thiazide; colecalciferol and 25-OH-D
- **acid-base core (5)** — ECF HCO₃⁻, non-bicarbonate buffer base, PaCO₂,
  proximal tubular cell pH, α-intercalated cell pH
- **renal actuator/adaptation (6)** — VH (the saturating actuator), pendrin,
  ammoniagenic capacity, NaDC1, fraction of functioning nephrons, interstitial
  fibrosis
- **electrolytes/volume (8)** — K⁺, total body K deficit, Cl⁻, ECF volume,
  aldosterone, ionised Ca²⁺, PTH, phosphate
- **bone (8)** — rapidly exchangeable base pool, bone mineral content,
  osteoclasts and osteoblasts, bone ALP, osteomalacia index, BMD z-score,
  cumulative base withdrawn
- **urine/stones (8)** — urinary Ca, citrate, pH and volume; brushite
  supersaturation; nephrocalcinosis; stone burden
- **systemic/clinical (10)** — IGF-1, height z-score, muscle strength, adherence,
  gastrointestinal irritation, hearing threshold shift, time below threshold,
  integrated acid exposure, wasted alkali, administered alkali

### The 28 scenarios

5 natural history (healthy child / untreated complete form in a child · severe
infant · adult / incomplete form) ·
5 diagnostic tests (NH₄Cl loading × 4 lesions, dietary acid titration) ·
10 schedule experiments (KHCO₃ vs NaHCO₃ vs K-citrate, TID/QID/BID/OD, ADV7103,
counterfactual with the prolonged-release component removed, counterfactual with
the immediate-release citrate removed, matched mean with differing trough) ·
3 clinical trial reproductions (the B21CS switch, B22CS at 6 years, Guittet urine
pH in healthy individuals) ·
5 comorbidity and co-medication (thiazide co-administration, low-acid diet,
acquired Sjögren's form + immunosuppression, delayed diagnosis + stage 3 CKD,
hypokalaemic paralytic crisis)

### The 10 tabs of the Shiny app

Patient profile · acid balance and the three sinks · alkali PK · urine chemistry ·
stone/nephrocalcinosis risk · **schedule laboratory** · **opposed-kinetics
sweep** · bone, growth and CKD · diagnostic tests · scenario library

At the head of each tab, the claim that tab is meant to test (or refute) is
stated explicitly. Tabs 6 and 7 were built so that the user can reproduce the
rejections of sections 3 and 4 above directly.

---

## 6. Running it

```bash
# mechanistic map
dot -Tsvg  drta_qsp_model.dot -o drta_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize="115,30!" drta_qsp_model.dot -o drta_qsp_model.png
```

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "shiny"))
library(mrgsolve)
mod <- mread("drta_mrgsolve_model.R")

# a single scenario
out <- sim_scenario(mod, "S15_ADV7103_BID")
plot(out, HCO3_e + urine_pH + UCit_day + SS_inst ~ time/24)

# summary table of the 28 scenarios
run_all_scenarios(mod)

# reproduction of the B21CS virtual cohort (schedule vs exposure decomposition)
res <- run_B21CS(mod, n = 40)
c(SoC = res$responder_soc, ADV = res$responder_adv)

# dashboard
shiny::runApp("drta_shiny_app.R")
```

---

## 7. Honesty note — limitations and known discrepancies

The environment in which this repository runs **had no R runtime.** All 56
equations were therefore first integrated in a dependency-free Python RK4
reference implementation, and calibrated there. In the course of that, **13
genuine defects** came to light, and at each place one was fixed a `BUG FIX #n`
comment records what was wrong and why. Representative ones:

1. the **sign of the buffer flux was reversed** — the buffer system became
   positive feedback and HCO₃⁻ diverged to 5×10⁹ mmol/L within 14 hours
2. nephron loss was driven by the **gap between target and current** fibrosis
   (which is 0 at steady state), so established nephrocalcinosis became harmless
3. respiratory compensation **applied Winter's formula at the normal point** — a
   healthy person sat at a pCO₂ of 45–47 and a pH of 7.36, slowly dissolving
   their own bone
4. the endogenous citrate correction was subtracted wrongly, feeding **every
   healthy individual a phantom acid of −57 mEq/day**, with the result that the
   V-ATPase controller was pinned at the ceiling in the baseline state
5. pure integral control (kI=0.22/h) traversed the whole range twice a day, so
   the controller degenerated into a **bang-bang oscillator**
6. `setup()` **overwrote the parameters with hard-coded constants**, silently
   discarding caller overrides, and double-scaled if called twice
7. the **daily integral of the meal-linked diurnal shape was 2 h rather than
   24 h** — the total NEAP entered at about 35 % of what was intended, so every
   patient looked far less acidotic than they really were
8. ammoniagenic capacity was driven **by the plasma HCO₃⁻ error alone** — and
   because the kidney is compensating, that error is small, so the ammonia arm
   could not be mobilised persistently and a high-protein diet alone left the
   controller railed
9. renal K⁺ secretion was **linear** in plasma K⁺ — plasma K⁺ was effectively
   defenceless, so K-citrate produced 6.6 mmol/L and Na alkali 1.4 mmol/L
10. the gastrointestinal irritation drive was normalised wrongly, so that **any
    immediate-release bolus saturated at the ceiling** — tolerability became an
    on/off switch rather than dose-dependent

### What could not be matched

- **The normalisation rate of the urinary Ca/citrate.** B21CS reported that 56 %
  of previous non-responders to SoC came below the stone-risk threshold, but in
  the virtual cohort the model barely produced that transition at all (0/17). The
  Ca/citrate of the severe patients is too far from 0.33.
- **The absolute level of the response rate.** A perfectly titrated cohort is
  100 % controlled in the model, so producing the observed 43 % baseline required
  real-world under-titration to be put in explicitly (`TIT_SHIFT`, `TIT_CV` — the
  **only** two values fitted to the trial, and fitted **only** to the 43 %
  baseline).
- **The paediatric alkali requirement.** The model's 10-year-old with the
  complete form of dRTA requires about 1.0 mEq/kg/day, whereas the clinical
  recommendation at this age is 1–3 mEq/kg/day. It only reaches the lower bound.
- **The BMD z-score breaks down as an extrapolation at large bone fluxes.**
  Because `kz_loss` was calibrated over the therapeutic range (fluxes of
  1–8 mEq/day), in the range where the flux is 30–55 mEq/day, as in the untreated
  severe form, the 1-year z-score becomes an unphysical value of anywhere from
  −10 to −20. It can be trusted within the treatment scenarios and cannot be
  used for untreated long-term prediction.
- **The absolute scale of the supersaturation is not calibrated.** Healthy
  individuals span an SS of 0.8–1.6 and do not sit clearly below the 1.0
  threshold. Only relative comparisons are valid.
- **The infant urine pH of 5.94** is somewhat lower than in reality. It is
  because the infant NEAP/kg is high at 2.0 mEq/kg/day and drives the actuator
  hard.
- **Hearing, enamel dysplasia and haemolysis** are included only as simple
  indices running from genotype to phenotype, and have no dynamics.
- `LES` and `LES_grad` were **not calibrated independently per genotype.** Their
  values match the qualitative ordering of severity in the literature.
- Fitted constants are marked `[FITTED]` and literature-based ones by their PMID.
  The calibration anchor table at the end of `drta_references.md` places the
  observed and modelled values side by side.

### How to falsify this model

1. If the plasma bicarbonate response rate does not come close to that of
   ADV7103 when standard therapy is given **under direct observation** → the
   exposure hypothesis of section 4 is wrong.
2. If there is no difference in **24-hour urinary calcium** between the
   immediate-release and prolonged-release forms at the same daily dose → the
   rectification mechanism of section 3 is wrong.
3. If, when the dietary acid load is raised in a patient with stones and
   incomplete dRTA, the urinary citrate and the bone turnover markers change
   **gently rather than in a step** → the actuator saturation structure of
   section 2 is wrong.

---

## ⚠️ Disclaimer

This is a QSP model for educational and research purposes. It has not undergone
independent validation and must not be used for clinical decision-making,
prescribing, or regulatory submission. The parameters are approximations for the
sake of explanation.
