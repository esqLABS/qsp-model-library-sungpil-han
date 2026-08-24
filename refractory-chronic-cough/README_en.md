# Refractory Chronic Cough — QSP Model

**Cough Hypersensitivity Syndrome**

> Chronic cough is not a symptom of another disease — it is **a disorder of the cough reflex
> arc itself.**

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`rcc_qsp_model_en.dot`](rcc_qsp_model_en.dot) · [SVG](rcc_qsp_model_en.svg) · [PNG](rcc_qsp_model_en.png) |
| ⚙️ mrgsolve model | [`rcc_mrgsolve_model.R`](rcc_mrgsolve_model.R) — 69 compartments · 241 parameters · 25 scenarios · 14 diagnostics |
| 📊 Shiny dashboard | [`rcc_shiny_app.R`](rcc_shiny_app.R) — 10 tabs |
| 📚 References | [`rcc_references.md`](rcc_references.md) — 120 PubMed citations |

---

## 1. Why this model (The premise)

For 50 years, chronic cough has been treated as **a problem of finding a cause** — find and
treat asthma, reflux, or post-nasal drip. This paradigm fails in 20–40% of specialist-clinic
patients, and the way it fails is highly suggestive: even after the "cause" is treated and its
objective evidence has disappeared, **the cough remains.**

The modern reinterpretation (Chung, Morice, McGarvey, Mazzone, Undem) views refractory chronic
cough as **a neuropathic disorder of the vagal reflex arc** — the airway's analogue of
neuropathic pain. Two features correspond exactly to pain.

| Cough | Pain correlate |
|---|---|
| **Hypertussia** — exaggerated coughing to a provoking stimulus | hyperalgesia |
| **Allotussia** — coughing triggered by talking, laughing, cold air, or perfume | allodynia |
| **Laryngeal paraesthesia** — a "tickle" | spontaneous pain |

The third is the key. This "tickle" arises from **ectopic spontaneous discharge** that occurs
even when there is nothing wrong with the airway. A patient can cough with an airway that is
completely normal.

---

## 2. Model structure — three layers in series

Every drug acts on exactly **one layer**, and that layer determines which endpoints it can
move.

```
LAYER 1  Peripheral gain   epithelial injury → pannexin-1 → ATP → vagal C-fibre P2X3
         (airway)          NGF/TrkA upregulates P2X3 expression itself over weeks
                           ▶ the layer P2X3 antagonists act on

LAYER 2  Central gain      afferent bombardment → nTS NMDA wind-up, substance P/NK-1,
         (brainstem)       microglial BDNF, loss of GABA/glycine inhibition
                           ▶ the layer that outlives the trigger
                           ▶ the layer gabapentin, morphine, and amitriptyline act on

LAYER 3  Cortical control  the drive first becomes a conscious "urge to cough", and
         (person)          only becomes motor coughing if prefrontal inhibition is insufficient
                           ▶ the layer speech therapy acts on — and so does the placebo response
```

$$\text{COUGH} = f\big(\underbrace{\text{afferent drive} \times \text{central gain}}_{\text{TUSS}} - \underbrace{\text{threshold}}_{\text{descending inhibition} \times \text{cortical inhibition}}\big) \times \text{arousal gate}$$

---

## 3. The central pharmacological claim of this model

### Efficacy and dysgeusia are **a single occupancy curve**, and the difference is **a single number**

Four P2X3 antagonists have entered phase 2b/3 trials for chronic cough. There is one efficacy
mechanism (blocking the vagal afferent **homotrimeric P2X3**) and one adverse-effect mechanism
(blocking **heterotrimeric P2X2/3** on taste-bud type II/III cells). Because ATP is **the
neurotransmitter** from taste bud to gustatory nerve for all five basic tastes, the taste system
has **no receptor reserve.**

| Drug | Selectivity (P2X2/3 : P2X3) | Efficacy | Taste adverse events |
|---|---|---|---|
| Gefapixant | ~6-fold | present (weak) | **58–69%** |
| Sivopixant | ~40-fold | unconfirmed | ~6–13% |
| Eliapixant | ~70-fold | present (weak) | ~10–21% · **discontinued for DILI** |
| Camlipixant | ~1500-fold | present (larger) | ~6% |

**How the model tests this claim** — fixing a single gefapixant 45 mg BID dose and sweeping
selectivity alone across four orders of magnitude (diagnostic 7):

| Selectivity | Cough reduction vs placebo | Taste adverse events | 12-week retention |
|---:|---:|---:|---:|
| 1-fold | −18.1% | **70.8%** | 79% |
| 6-fold (gefapixant) | −18.6% | **52.2%** | 89% |
| 70-fold (eliapixant) | −19.2% | 3.6% | 100% |
| 1500-fold (camlipixant) | −19.2% | **3.0%** | 100% |

**Efficacy stays essentially flat while the taste adverse events collapse.** This is why 6-fold
selectivity looks "somewhat selective" but is nowhere near enough.

### Class ceiling

The P2X3 pathway carries only **part** of the cough drive. Even with complete target
occupancy, the TRPV1/TRPA1/mechanoreceptor arm remains untouched. Fitting the COUGH-1 phase 3
result brings the P2X3-dependent fraction in a typical patient down to **about 5%**
(`WP2X = 0.047`). This is why this drug class never shows morphine-level cough suppression.

---

## 4. The placebo arm is not a number, it is **a mechanism**

Chronic cough is the area in respiratory medicine with **the largest placebo response in an
objective physiological measure** (machine-counted 24-hour cough frequency falls by 30–45%).
Treating the placebo arm as a mere number overestimates every drug effect by more than 2-fold.
Three mechanisms all operate in the placebo arm.

1. **Patients are enrolled at their worst** — `FLARE` decays with a half-life of ~41 days
   (regression to the mean, expressed as a mechanism)
2. **Trial participation itself trains inhibition** — diaries, monitoring, and contact with
   clinicians raise `CORT` with a time constant of ~3 weeks (incidentally achieving what speech
   therapy does deliberately)
3. **The natural resolution of the trigger**

⚠️ **A screening-to-randomisation run-in was essential.** Injecting the flare at the moment of
randomisation makes the placebo arm **worsen by 64%** over 12 weeks — the opposite of real
cough trials. Patients are screened when symptoms are worst and randomised **weeks later**,
with baseline already measured while the flare is declining. The model takes day 0 of the
trial to be the state after 28 simulated days of flare. With this correction, the placebo arm
becomes **−32%** at 12 weeks.

---

## 5. Validation — a head-to-head comparison against published trials

All values are **placebo-adjusted, at the corresponding population and time point.**

| Trial | Endpoint | Observed | Model |
|---|---|---:|---:|
| **COUGH-1 (12 weeks)** | Gefapixant 45 mg BID, 24h cough vs placebo | **−18.5%** | **−18.6%** |
| COUGH-1 | Gefapixant 45 mg, taste adverse events | 58% | 52% |
| PAGANINI (12 weeks) | Eliapixant 75 mg BID | −17.6% | −22.9% |
| PAGANINI | Eliapixant 150 mg, taste adverse events | ~21% | 26% |
| Sivopixant 2b (12 weeks) | 150 mg QD | −18.6% (ns) | −22.9% |
| **Ryan 2012 (10 weeks)** | Gabapentin 1800 mg/d, LCQ | +1.80 | **+1.99** |
| Ryan 2012 | Gabapentin, cough VAS (mm) | −12.1 | **−11.3** |
| Ryan 2012 | Gabapentin, cough frequency | −27% | **−24.6%** |
| Morice 2007 (4 weeks) | Morphine SR 10 mg BID, LCQ | +3.2 | +2.88 |
| Vertigan 2016 (4 weeks) | Pregabalin + speech therapy, LCQ | +2.5 | **+2.23** |
| PSALTI (4 weeks) | Speech therapy, LCQ | +1.53 | +1.12 |
| **CANAL (IPF cough)** | Nalbuphine ER, weekly cough vs placebo | −53% | **−53.6%** |
| **ICS, non-eosinophilic RCC** | 24h cough vs placebo | **no effect** | **−0.0%** |
| ICS, eosinophilic cough | 24h cough vs placebo | large effect | −18.8% |
| **PPI, non-acid reflux** | 24h cough vs placebo | **no effect** | −5.6% |
| PPI, confirmed acid reflux | 24h cough vs placebo | moderate | −10.3% |
| ACEi withdrawal | cough at 4 weeks (vs baseline) | mostly resolved | 40% |

### The **negative controls** are a stronger validation than the positive ones

The lack of effect of ICS and PPI **was not dictated to the model.** `PPI` suppresses only the
`ACIDR` state and has **no term at all** acting on `NACID`, and `ICS` acts only through the
eosinophil set point (`ICS_NFK = 0`). Both null results fall directly out of the structure.

### What does not fit — recorded without hiding it

| Trial | Observed | Model | Problem |
|---|---:|---:|---|
| COUGH-2 (24 weeks) gefapixant 45 | −14.6% | −25.8% | **too strong** |
| SOOTHE (4 weeks) camlipixant 50 | ~−34% | −10.5% | **too weak** |
| COUGH-1 gefapixant **15 mg** | −1.4% (ns) | −13.0% | cannot reproduce the dose cliff |

The first two are in **opposite directions.** The model's drug effect **builds too slowly and
then keeps growing**, whereas the real drug class reaches most of its effect by ~4 weeks and
then plateaus or declines. Sweeping `WP2X` and the wind-up decay rate moves both numbers
**together**, so fixing one worsens the other — this is a **missing mechanism (tolerance or
receptor escape)**, not a parameter misspecification. **Do not trust P2X3 predictions beyond
12 weeks.**

The 15 mg cliff cannot be reproduced by any occupancy-based model (15 mg BID achieves
substantial P2X3 occupancy). This is not a defect of this implementation but **an unsolved
problem in the field**, and it has been left unfitted.

---

## 6. The model **disproved** three of its own design assumptions

The value of a QSP model lies in revealing what is wrong, not just what is right. The three
assumptions below were all written as claims in the header, and all were retracted by the
model's own diagnostics.

### ① Bistability and hysteresis — **absent**

It was built expecting the vicious circle to produce two stable states. Sweeping loop gain
**upward** from health and **downward** from established disease, the two overlapped **within
0.3 coughs/h.** This model has no two attractors — only **a steep, monotonic sensitivity
curve.** (An earlier version claimed the bifurcation parameter was *trigger magnitude*, and
that too was disproved by a sweep — the trigger disappears within weeks and plays no role in
sustaining what follows.)

The resolution of post-infectious cough (1.8 coughs/h) versus refractory persistence
(21 coughs/h) **is reproduced** — not by crossing a fold, but by moving along a continuous
curve made by loop gain × peripheral sensitivity (oestrogen status, age).

### ② Efficacy ≠ effectiveness in magnitude — **almost no difference**

Dysgeusia → discontinuation → reduced exposure → dilution of the ITT estimate. This gap was
expected to be a large part of the reason gefapixant's effect is small. **It is not**: ITT
−18.6% vs per-protocol −19.2%, a difference of 0.6 percentage points. Almost all of the small
effect is because **the P2X3-dependent fraction of the drive is small**, not because of
dropout.

### ③ Benefit of early treatment — **absent**

Expectation: treating before central sensitisation sets in should work better. Result: −21.2%
at 10 weeks of illness vs −21.4% at 5 years. The reason shows up in the columns printed
alongside — wind-up and microglial activation **already reach steady state by 10 weeks**
(0.31 vs 0.29). If early clinical treatment really is better, the mechanism for that is
**absent from this model.**

In addition, diagnostic 12 showed that the **direction** of the objective/subjective
dissociation was the opposite of what was expected (it is peripheral drugs, not cortical
intervention, that give the larger count:urge ratio), and this was corrected in the header
along with the reason.

---

## 7. Prediction: challenge-test dissociation

Capsaicin opens TRPV1 directly and **bypasses P2X3.** Inhaled ATP opens P2X3. Because the model
computes both thresholds from the same state vector, this dissociation is **predicted, not
dictated.**

| Treatment | Change in capsaicin C5 | Change in ATP threshold | 24h cough |
|---|---:|---:|---:|
| Placebo | 1.00× | 0.99× | 0.83× |
| Gefapixant 45 | **1.01×** | **1.13×** | 0.75× |
| Camlipixant 50 | **1.01×** | **1.20×** | 0.74× |
| Morphine 10 | 1.04× | 1.03× | **0.48×** |

This is why the capsaicin challenge is an unsuitable pharmacodynamic biomarker for this drug
class. ⚠️ Note, however, that the model's ATP shift (1.1–1.2-fold) is far smaller than the
near one-order-of-magnitude log shift reported experimentally — **read only the direction and
the dissociation.**

---

## 8. Disease is **generated, not assigned**

No virtual patient is initialised in a diseased state. All start healthy, receive a trigger,
and are simulated for `PRERUN` days. The cough at the moment of enrolment is **the cough the
model produced.**

365-day drift of a healthy subject (`INSULT0 = 0`): 24h cough frequency **0.37%**, central
gain **0.084%**, loop gain **8×10⁻¹¹%**, capsaicin C5 **0.046%**.

| Phenotype | 24h cough | Waking cough | LCQ | VAS | Central gain | Capsaicin C5 |
|---|---:|---:|---:|---:|---:|---:|
| Healthy (female) | 1.9 | 2.7 | 20.5 | 1 | 1.00 | 26.0 |
| Healthy (male) | 1.4 | 1.9 | 20.6 | 1 | 1.00 | 61.9 |
| **Refractory chronic cough** | **21.1** | **29.4** | **11.0** | **55** | 1.44 | 17.0 |
| Post-infectious (resolved) | 1.8 | 2.6 | 20.6 | 1 | 0.98 | 26.3 |
| Eosinophilic (NAEB/CVA) | 16.6 | 23.1 | 12.4 | 45 | 1.20 | 18.0 |
| Acid reflux | 13.9 | 19.4 | 13.3 | 39 | 1.25 | 19.7 |
| Non-acid reflux | 13.7 | 19.0 | 13.4 | 38 | 1.25 | 19.7 |
| ACEi cough | 11.7 | 16.3 | 14.3 | 33 | 1.20 | 14.9 |
| Cough with IPF | 32.7 | 45.5 | 7.9 | 79 | 1.56 | 45.6 |
| Enrolled trial population | 21.8 | 30.4 | 10.8 | 56 | 1.45 | 17.0 |

It is deliberate that the acid/non-acid reflux pair is set to start from **nearly identical
cough** (13.9 vs 13.7) — so that the PPI comparison reflects only the **mechanistic**
difference, not a baseline difference.

---

## 9. Circadian structure is also **generated**

Cough almost ceasing during sleep is because the cortical urge pathway shuts off, and the
morning rise follows overnight mucus accumulation. Neither is imposed — both emerge.

- Waking : 24-hour cough frequency ratio = **1.39**
- Sleep nadir = **1.7 coughs/h** (vs ~30 while awake)
- Post-waking peak = 4 hours after waking

⚠️ The real peak is **immediately** after waking, whereas the model's is 4 hours later — the
mucus-clearance kinetics are too slow.

---

## 10. Usage

```r
# Required packages: mrgsolve, ggplot2, shiny
source("rcc_mrgsolve_model.R")

report()                      # all 14 diagnostics (about 3-5 minutes)

check_baseline()              # 1. healthy-subject drift (~0%)
check_phenotypes()            # 2. phenotype baselines produced by the natural history
check_bistability()           # 3. sensitivity sweep — no bistability (negative result)
check_window()                # 7. ★ therapeutic window — selectivity sweep
check_eff_vs_effectiveness()  # 8. ITT vs per-protocol (negative result)
check_challenges()            # 9. capsaicin vs ATP dissociation
validate_trials()             # 10. head-to-head comparison against published trials
check_early_late()            # 11. benefit of early treatment (negative result)
check_sensitivity()           # 12. sensitivity of the headline results
check_pro_dissociation()      # 13. objective/subjective dissociation (direction corrected)
check_circadian()             # 14. circadian structure

# Individual scenario
ic <- settle(pat_trial)                                   # natural-history pre-run
d  <- daily(run(pat_trial, list(), ev_gefapixant(45),
                days = 84, ic = ic))

# Dashboard
shiny::runApp("rcc_shiny_app.R")
```

Map rendering:
```bash
dot -Tsvg rcc_qsp_model_en.dot -o rcc_qsp_model_en.svg
dot -Tpng -Gdpi=150 rcc_qsp_model_en.dot -o rcc_qsp_model_en.png
```
> `newrank=true` is required. There is a cycle between clusters (e.g. baclofen → TLESR), which
> makes graphviz's default ranking fail with `trouble in init_rank`.

---

## 11. Two important implementation pitfalls

**① `$MAIN` runs on every record.** Leaving the initial-condition block as is **silently
overwrites** the state vector passed via `init()`. The entire 730-day natural-history pre-run
gets discarded and a healthy subject is simulated instead, with no error raised. Every `*_0`
assignment is guarded by the `SETIC` parameter.

**② The self-injury term must be saturated.** Writing cough → mechanical epithelial injury as
a proportional term makes a patient coughing 30 times an hour strip the airway ten times over.
Self-injury without saturation makes a model in which the upper stable state does not exist at
all. Every self-injury term takes the form `SAT(x, k)`, with a `CTOL` dead-band around the
healthy operating point.

---

## 12. Limitations

- This is a **typical-subject model.** Cough counts follow a log-normal distribution across
  patients, and this model predicts only the geometric-mean trajectory — it **does not predict
  the distribution of responders.** Between-patient variability is supplied as parameter sets,
  not an OMEGA block.
- Adherence is applied as a **population-average exposure multiplier.** This is a valid
  approximation for the ITT geometric mean and an invalid one for any individual patient.
- The selectivity ratios, IC50s, and unbound fractions of the four P2X3 antagonists were
  collected from preclinical reports run under different conditions and are **the most
  uncertain parameters in the model.** The qualitative result on the therapeutic window is
  robust; the predicted taste-adverse-event percentages are not.
- Eliapixant hepatotoxicity is idiosyncratic and is included only **as a flag, not a
  dose-response.**
- Speech therapy reproduces the LCQ effect (+1.12 vs +1.53) but fails to reproduce the
  objective cough-frequency effect (−16.8% vs −41%).
- The list of parameters without experimental support is set out in the final table of
  [`rcc_references.md`](rcc_references.md).

---

## ⚠️ Disclaimer

A qualitative/semi-quantitative QSP model for education and research purposes. It has not been
independently validated or certified, and **must not be used for actual clinical
decision-making, prescribing, or regulatory submission.**
