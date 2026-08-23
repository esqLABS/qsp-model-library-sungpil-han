# Refeeding Syndrome (RFS) QSP Model

> **The disease is a flux mismatch, not a low number on a blood test.**
> Λ<sub>P</sub> = J<sub>demand</sub> / J<sub>supply</sub>. The clinician sets the
> numerator through the glucose infusion rate; the patient's starvation history
> sets the denominator. Neither of them can see the problem in the admission
> panel, because the compartment the laboratory measures holds **0.06 %** of the
> phosphorus that is actually at risk.

<p align="center">
  <a href="../../../refeeding-syndrome/rfs_qsp_model.svg"><img src="../../../refeeding-syndrome/rfs_qsp_model.png" width="900" alt="Refeeding syndrome QSP mechanistic map"></a>
</p>

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map (147 nodes · 17 clusters · 229 edges) | [`rfs_qsp_model.dot`](../../../refeeding-syndrome/rfs_qsp_model.dot) · [`.svg`](../../../refeeding-syndrome/rfs_qsp_model.svg) · [`.png`](../../../refeeding-syndrome/rfs_qsp_model.png) |
| ⚙️ mrgsolve ODE model (48 ODEs · 172 params · 16 scenarios) | [`rfs_mrgsolve_model.R`](../../../refeeding-syndrome/rfs_mrgsolve_model.R) |
| 📊 Shiny dashboard (10 tabs) | [`rfs_shiny_app.R`](../../../refeeding-syndrome/rfs_shiny_app.R) |
| 📚 References (106 refs · 96 PubMed links) | [`rfs_references.md`](../../../refeeding-syndrome/rfs_references.md) |
| 🐍 Independent validation implementation (Python/scipy) | [`rfs_reference_model.py`](../../../refeeding-syndrome/rfs_reference_model.py) → [`rfs_reference_output.txt`](../../../refeeding-syndrome/rfs_reference_output.txt) |

---

## 1. The thesis

Almost every clinical guideline describes refeeding syndrome as
**"the disease where serum phosphate falls."** This model rewrites it as a
**flux mismatch**.

```
Λ_P  =  J_demand / J_supply
```

* **J_demand** — the phosphate that must enter the cell because glycolytic
  flux pulls up the intracellular organic-phosphate set point (G6P,
  F1,6-bisP, ATP, creatine phosphate, 2,3-DPG), plus the phosphate deposited
  into newly synthesised lean mass. **Set by the clinician, through the
  glucose infusion rate.**
* **J_supply** — absorbed dietary phosphate + prescribed intravenous
  phosphate + net efflux from bone + what the kidney spares by shutting down
  excretion. **Set by the patient's history.**

This ratio is dangerous for three arithmetic reasons, and all three are
stoichiometry, not fitted parameters.

**(1) The compartment that is measured is not the compartment at risk.**

| Compartment | Amount | Note |
|------|-----|------|
| Extracellular inorganic phosphate | **14 mmol** | ← **this is the only thing measured** |
| Intracellular phosphate | 3,216 mmol | 67 mmol/kg lean body mass |
| Bone phosphate | 19,550 mmol | 85 % of the total |
| **Total body phosphorus** | **22,780 mmol** | |

What the lab sees is **0.063 %** of the total. A daily refeeding demand of
60 mmol equals **4.2 times the entire measured compartment per day**.

**(2) The kidney's reserve is a one-time spend.** Urinary phosphate
excretion is a threshold function,
`GFR × max(0, Pser − TmP/GFR)`. The most the kidney can contribute is "stop
excreting", and that amount is about 14 mmol/d and is **used up within 24
hours.** After that, no matter how far phosphate falls, the kidney has
nothing left to give.

**(3) Bone is slow, and magnesium can switch it off.** Net efflux from bone
is about 10 mmol/d and is PTH-dependent, but hypomagnesaemia suppresses PTH
secretion. That is, the last remaining endogenous supply line is blocked
**by the second electrolyte that the same syndrome is simultaneously
depleting.**

### The two clocks

Thiamine and phosphate both take **about three weeks to deplete**, but
**their refilling time constants differ by two orders of magnitude.** At
supraphysiological plasma concentrations, thiamine bypasses the saturable
ThTR-1/ThTR-2 transporters and enters tissue by **passive diffusion**,
whereas phosphate must be **pumped into 40 kg of cells** against saturable
transporters.

---

## 2. What was fitted, and what was not

**Only 3 numbers were fitted to refeeding syndrome itself.**

| Parameter | Value | Role |
|----------|-----|------|
| `SINS` | 0.750 | how much glycolytic flux raises the cellular phosphate set point |
| `KFILL` | 0.0145 /h | the rate at which the cell approaches that set point (t½ 48 h) |
| `HAZ` | 1.00 | overall scale of the mortality hazard function |

Everything else came from **normal physiology** (total-body P/K/Mg and their
bone/intracellular/extracellular fractions, 67 mmol P and 72 mmol K per kg
lean body mass, serum reference ranges, TmP/GFR, GFR of 6 L/h, the urinary
P 29 · K 70 · Mg 4.9 mmol/d that balance a normal diet, Cunningham REE,
Forbes fat/lean partitioning, an insulin half-life of 5 minutes) or from
**pharmacological/nutritional data unrelated to refeeding** (P 22 · K 38 ·
Mg 6 mmol per 1000 kcal of enteral formula, saturable oral thiamine
absorption, the 0.6–1.0 mmol/L fall in serum potassium seen when
insulin-glucose is given for hyperkalaemia, the Ca×P solubility product,
QTc sensitivity).

### An emergent check

Thiamine body stores (26.5 mg) and biological half-life (14 days) were taken
from tracer studies. **The daily requirement was not.** From those two
values alone, the model computes a requirement of **1.42 mg/d**, and the
published RDA is **1.1–1.4 mg/d**.

The healthy steady state is a **numerically exact steady state** across all
44 physiological state variables (max |dy/dt| = 6.2 × 10⁻¹⁵ /h).

---

## 3. Five real defects that independent validation found

Every equation was independently re-implemented in Python/scipy
([`rfs_reference_model.py`](../../../refeeding-syndrome/rfs_reference_model.py)). In the
process, **five real defects** surfaced, and all five were fixed.

1. **The cellular phosphate set point had been written against a "well-fed
   control" baseline.** As a result, refeeding a starved patient at
   30 kcal/kg/d never crossed the set point, and **phosphate did not fall at
   all** — the disease mechanism itself was not firing.
2. The set point was then driven by **insulin**, which has **the wrong
   sign**. When hypophosphataemia suppresses glucose disposal, blood glucose
   and insulin **rise**, so demand grows at exactly the moment the cell can
   least meet it. It is now driven by **glycolytic flux**.
3. Renal potassium excretion had **no depletion adaptation.** During
   starvation, intracellular potassium fell to 23 % of normal — a value
   incompatible with survival.
4. **PaCO₂ appeared as an unbounded ratio inside an exponential hazard
   function.** Two scenarios produced 100 % mortality, which was purely a
   numerical artefact.
5. **75 % of insulin secretion sat in a glucose-independent basal term**, so
   insulin — and phosphate demand as a whole — barely responded to feeding.

---

## 4. Key findings (verified numbers)

### 4.1 Admission serum phosphate carries almost no information

Six patients whose starvation duration spans a **4-fold** range were given
**identical**, unprophylaxed refeeding:

| History | Admission P | Admission K | Intracellular K | P nadir | Mortality |
|------|--------|--------|----------|----------|--------|
| 30 d @ 50 % | 1.03 | 4.08 | 82.7 % | 0.68 | 0.41 % |
| 60 d @ 32 % | 1.02 | 3.60 | 75.6 % | 0.45 | 1.07 % |
| 120 d @ 22 % | 1.02 | 2.62 | 71.8 % | 0.71 | **1.46 %** |

Admission phosphate is essentially unmoved, **1.02–1.03** (within measurement
error, dead centre of the reference range). Mortality across the same six
patients differs by **3.6-fold**.

> **A counter-to-expectation result, reported as found:** the phosphate
> nadir is **not monotonic** in starvation duration. It is deepest at
> 60–75 days and becomes **shallower again** by 120 days. The most severely
> depleted patients have the lowest lean mass and insulin sensitivity, so
> they cannot generate the glycolytic flux that creates the demand in the
> first place. Yet mortality keeps climbing — driven by potassium,
> magnesium, and thiamine. **In the highest-risk patients, a shallow
> phosphate nadir is not reassurance — it is a sign that serum phosphate is
> the wrong indicator.**

### 4.2 Demand tracks glucose, not calories

Total energy was **fixed at 25 kcal/kg/d** while only the carbohydrate
fraction was varied:

| CHO % | GIR (mg/kg/min) | Λ_P | P nadir | Mortality |
|-------|-----------------|-----|----------|--------|
| 20 | 0.87 | 0.00 | 1.02 | 0.56 % |
| 40 | 1.74 | 0.89 | 0.96 | 0.60 % |
| **50** | **2.17** | **1.54** | **0.64** | 0.71 % |
| 60 | 2.60 | 2.18 | 0.42 | 1.29 % |
| 100 | 4.34 | 4.12 | 0.13 | **6.44 %** |

Every row delivers the same calories. **"How many calories" is the wrong
question; "how many mg/kg/min of glucose" is the right one.** The nadir
breaks precisely where Λ_P crosses 1 (CHO 40–50 %).

### 4.3 The feed itself is a phosphate prescription

**The same 20 kcal/kg/d, the same patient**, with only the route and
composition changed:

| Prescription | P nadir | GIR | Lactate | Mortality |
|------|----------|-----|------|--------|
| IV glucose alone, 100 % CHO | **0.07** | 3.47 | 10.5 | 5.75 % |
| Enteral formula, 100 % CHO | 0.19 | 3.47 | 6.9 | 3.73 % |
| Enteral formula, 40 % CHO | **1.00** | 1.39 | 4.0 | 0.64 % |

A complete enteral formula carries P 22 mmol · K 38 mmol · Mg 6 mmol per
1000 kcal, plus thiamine, **all together.** A bag of IV dextrose carries
**nothing at all** unless someone adds it.

### 4.4 A 2×2: which prophylaxis actually does the work

Aggressive refeeding (30 kcal/kg/d from hour 0) × electrolyte
supplementation × thiamine:

| | P nadir | Lactate | Mortality | Wernicke |
|---|---|---|---|---|
| Neither | 0.45 | 5.0 | 1.07 % | 4 % |
| Electrolytes only | **1.02** | 4.5 | 0.54 % | 4 % |
| Thiamine only | 0.45 | **3.9** | 0.51 % | **0 %** |
| Both | **1.02** | **3.9** | **0.36 %** | **0 %** |

The two prophylactic measures are **neither interchangeable nor
overlapping.** Electrolytes fix the phosphate nadir and do nothing for
Wernicke's. Thiamine fixes Wernicke's and lactate and does **nothing for
phosphate.**

### 4.5 NICE 2006 vs ASPEN 2020

| | NICE-type | ASPEN-type |
|---|---|---|
| Starting kcal/kg/d | 10 | 15 |
| P nadir | 1.01 | 1.02 |
| Energy delivered over 14 days | 17,641 | **19,967** |
| Mortality | 0.36 % | 0.36 % |

**Provided supplementation is co-prescribed**, faster advancement delivers
**13 % more energy** at the same risk. As a control, **slow advancement
without supplementation** is worse than fast advancement with
supplementation — P nadir 0.49, mortality 0.82 %, Wernicke's 7 %. This
model's reading of the NICE/ASPEN debate is that **supplementation, not the
calorie ceiling, was the active ingredient.**

### 4.6 Magnesium: how it is given, and the door that opens

The same 0.4 mmol/kg/d given as a 2-hour bolus vs a 24-hour continuous
infusion:

* Bolus — Mg **0.72** mmol/L on day 3
* Continuous infusion — Mg **0.88** mmol/L on day 3

Renal magnesium handling is a **threshold function**, so the bolus's peak is
not retained and is excreted. The prescription records the same daily dose,
but **the patient does not receive it.**

---

## 5. Results that ran against expectation (reported, not tuned away)

### 5.1 The prediction about oral thiamine **failed**

This section was originally written to show that "oral thiamine cannot
reverse acute deficiency." **The model disagreed, and that is reported as
is.**

| Regimen | TK 3 h | TK 12 h | TK 24 h | TK 3 d | Wernicke |
|------|--------|---------|---------|--------|----------|
| None | 0.01 | 0.02 | 0.04 | 0.26 | 10 % |
| Oral 300 mg/d | 0.23 | 0.97 | 1.00 | 1.00 | 0 % |
| Oral 1500 mg/d | 0.28 | 0.98 | 1.00 | 1.00 | 0 % |
| IV 100 mg/d | **0.70** | 1.00 | 1.00 | 1.00 | 0 % |
| Oral 300, vomiting | 0.01 | 0.04 | **0.14** | 0.69 | **5 %** |
| IV 100, vomiting | **0.70** | 1.00 | 1.00 | 1.00 | 0 % |

The saturable transporter was **not the rate-limiting step.** Its ceiling is
about 48 mg/d of absorption, while the daily requirement is only 1.4 mg —
even saturated absorption delivers roughly **30 times** the requirement.
That raising the oral dose from 300 to 1500 mg/d changes nothing was
expected, but that **switching to IV made no difference either (at the
24-hour mark)** was not.

The IV route earns its keep in the last two rows. When absorption is
impaired (vomiting, ileus, alcohol-related mucosal injury), the oral arm
collapses and the IV arm does not. And at the 3-hour mark, IV is still much
faster — 0.70 vs 0.23. The model's verdict is that **"IV because oral can't
deliver enough" is the wrong reason, and "IV because this patient's gut
can't be trusted" is the right one.** Since the population in question is
exactly the latter, **the clinical recommendation survives, but its
justification is replaced.**

### 5.2 "Thiamine before glucose" buys nothing on the timescale usually taught

| Thiamine timing | Lactate | Wernicke | Mortality |
|------------------|------|----------|--------|
| 1 h **before** glucose | 3.9 | 0 % | 0.51 % |
| **Simultaneous** with glucose | 3.9 | 0 % | 0.51 % |
| 6 h late | 4.5 | 1 % | 0.55 % |
| 24 h late | 5.0 | 2 % | 0.72 % |
| Never given | 5.0 | 4 % | 1.07 % |

Giving it 1 hour early and giving it simultaneously are **indistinguishable.**
What the model supports is the weaker, and probably the true, claim —
**thiamine must be given on day 0.**

Phosphate runs the opposite way. Starting on day 0 abolishes the nadir
(0.97), but even a 24-hour delay does not (0.48), because by then the cell
has already drawn what it needs from a 14-mmol pool. Phosphate given late
still normalises the day-7 value — that is, **the lab recovers fast and the
patient recovers slowly.**

### 5.3 The U-shaped curve for phosphate supplementation is real but lies outside the clinical range

Benefit **saturates at 0.5 mmol/kg/d** (the guideline-recommended range,
which the model never sees exceeded in practice). The Ca×P solubility
product is only crossed at **12 mmol/kg/d — 20 times the recommended
dose**, and only there does ionised calcium fall and mortality rise again.

Within this model the asymmetry is one-sided: **under-replacement costs
0.50 percentage points of mortality, and over-replacement carries no
measurable cost until the dose becomes absurd.** However, the model has
neither soft-tissue calcification, nor phosphate-infusion arrhythmia, nor
variable renal function — all three are real reasons clinicians limit the
rate. So this curve's flat top should be read as
**"this model cannot see the harm," not "there is no harm."**

---

## 6. Model structure

| Block | State variables | Key point |
|------|----------|------|
| Feeding · gut | `A_CHO A_FAT A_PRO PGUT KGUT MGGUT THGUT` | whether micronutrients arrive together depends on the route |
| Glucose · insulin | `GLU INS X GCG GLY FFA BHB` | the incretin term (35 %) distinguishes enteral from IV |
| Phosphate | `PE PI PB` | set-point tracking + threshold-type renal excretion + PTH-dependent bone efflux |
| Potassium | `KE KI` | insulin-driven uptake + K-driven aldosterone + ROMK downregulation |
| Magnesium | `MGE MGI` | threshold-type renal excretion → a bolus is wasted |
| Calcium axis | `CAE CAB PTH CTD FGF` | Mg → PTH → bone phosphate supply line |
| Thiamine | `THP THT` | saturable transporter + passive-diffusion pathway |
| Metabolism | `LAC` | TPP-deficient PDH gate → type B lactic acidosis |
| Organ energetics | `ATPM PCRM DPG ATPD` | requires **both** Pi and TPP |
| Body composition · heart | `FM FFM AT LVM NAB ALD` | atrophied left ventricle + insulin-driven sodium retention |
| Outcomes | `H_ARR H_CF H_RF H_WE SURV AUCP CUMP CUMK CUMKCAL` | three pathways: arrhythmia, heart failure, respiratory failure |

**16 scenarios**: NICE-type / ASPEN-type / aggressive+no prophylaxis /
aggressive+electrolytes / aggressive+thiamine / aggressive+both /
slow-advance+no-prophylaxis / IV-glucose-alone / enteral-isocaloric /
mixed-fuel / alcohol±IV-thiamine±oral-thiamine / Mg bolus vs infusion /
K-without-Mg supplementation.

---

## 7. How to reproduce

```bash
# Mechanistic map
dot -Tsvg rfs_qsp_model.dot -o rfs_qsp_model.svg
dot -Tpng -Gdpi=150 rfs_qsp_model.dot -o rfs_qsp_model.png

# Independent validation implementation (numpy + scipy)
python3 rfs_reference_model.py > rfs_reference_output.txt

# mrgsolve model
Rscript rfs_mrgsolve_model.R

# Shiny dashboard
R -e 'shiny::runApp("rfs_shiny_app.R")'
```

---

## 8. Limitations (stated, not buried)

* **`rfs_mrgsolve_model.R` was written in a container that cannot run R, so
  it has not been compile-verified.** The equations themselves were
  verified via the independent Python implementation (all 48 ODEs confirmed
  to have a `dxdt_` assignment, parentheses/braces confirmed balanced), and
  the two files come from the same equation list. Even so, the possibility
  of a remaining C++ compilation error cannot be excluded.
* In the worst-case scenario (IV glucose with no phosphate at all, no
  supplementation), the nadir reaches as low as 0.07 mmol/L. Reported
  real-world nadirs are 0.10–0.30. The model has none of the actual
  clinician's rescue interventions, and only two brakes — transporter
  saturation and glycolytic block — so it overestimates at the extremes.
* Mortality is a hazard model with a **single fitted scale**. The
  arrhythmia/heart-failure/respiratory-failure branches are structurally
  separated but are not each independently validated numbers.
* There is no sepsis, multi-organ failure, delirium, refeeding-related
  fatty liver, or any vitamin deficiency other than thiamine. Real
  refeeding deaths frequently involve these.
* The starvation phase assumes a constant intake fraction. Real histories
  are episodic, and a binge-restrict pattern would deplete thiamine faster
  than this model does.
* Bone is a single mixed compartment. Because chronic-starvation
  osteopenia is not represented, the model likely **overestimates bone
  buffering capacity** in the most chronically depleted patients.
* Intracellular phosphate is a single soft-tissue pool. Erythrocyte
  2,3-DPG, myocardium, and diaphragm each read off serum phosphate rather
  than having their own compartment.
* The model does not reproduce **magnesium-refractory hypokalaemia** as
  strongly as the textbook account (day-3 values of 3.77 vs 3.86 mmol/L). It
  may be underestimating a phenomenon that is commonly reported in practice.
* Renal function is fixed. The acute kidney injury common in this
  population would shift every threshold in the model.

---

> ⚠️ **This is a QSP model for education and research purposes. It must
> not be used for clinical decision-making, prescribing, or regulatory
> submission.**
