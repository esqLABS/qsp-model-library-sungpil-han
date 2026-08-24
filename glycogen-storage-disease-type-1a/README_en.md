# Glycogen storage disease type Ia (von Gierke disease) QSP model
### Glycogen Storage Disease Type Ia · von Gierke Disease · G6PC1 deficiency

<p align="center">
  <a href="gsd1a_qsp_model.svg">
    <img src="gsd1a_qsp_model.png" width="880" alt="GSD Ia QSP mechanistic map">
  </a>
</p>

| File | Contents |
|---|---|
| [`gsd1a_qsp_model.dot`](gsd1a_qsp_model.dot) · [SVG](gsd1a_qsp_model.svg) · [PNG](gsd1a_qsp_model.png) | Mechanistic map — 167 nodes, 17 clusters |
| [`gsd1a_mrgsolve_model.R`](gsd1a_mrgsolve_model.R) | mrgsolve ODE model (47 compartments) + 14 scenario functions |
| [`gsd1a_shiny_app_en.R`](gsd1a_shiny_app_en.R) | Shiny dashboard (11 tabs) |
| [`gsd1a_references_en.md`](gsd1a_references_en.md) | 100 references, including a component correspondence table |
| [`gsd1a_reference_model.py`](gsd1a_reference_model.py) | Independent Python/scipy re-implementation — for verification |
| [`gsd1a_reference_output.txt`](gsd1a_reference_output.txt) · [`gsd1a_scenario_results.json`](gsd1a_scenario_results.json) | The computed output of every number below |

---

## What this model says, in one sentence

**Glucose-6-phosphatase (G6Pase) is the single final step through which glycogenolysis
and gluconeogenesis <ins>both</ins> have to pass. Lose it and hepatic glucose
production does not fall — the liver switches from being a <ins>source of glucose to
being a sink</ins>** — because glucokinase keeps turning in the forward direction and
there is no enzyme to reverse it.

Everything else is overflow. So this model is written not as four sub-models but as
**one mass balance on hepatic G6P**, and the four classical syndromes fall out of the
branching ratios by themselves.

```
                       ┌─ free glucose      ✖ blocked (0 %)
    liver G6P ─────────┼─ glycogen          50 %  →  hepatomegaly
    (one branch point) ├─ lactate           46 %  →  lactic acidosis
                       └─ pentose phosphate  4 %  →  hyperuricaemia
                            └ acetyl-CoA → lipid → hypertriglyceridaemia
```

These are the actual branching ratios the model computes at steady state in a treated
GSD Ia patient — read off the flux decomposition (`fluxes()`), not assumed values.

---

## How it was verified

All 47 ODEs were implemented independently **twice** — once inside mrgsolve's C++, once
in Python/scipy. It is the cheapest way to catch the kinds of error QSP models commit
most readily (dimensionally wrong equations, parameters that make a clinically
catastrophic drug look beneficial). Every number below was obtained by integration, and
at the bottom of this document **the seven defects caught in the process and the one
refuted hypothesis** are recorded.

---

## Result 1 — the cover time of cornstarch is set by <ins>the release rate, not the total amount</ins>

The commonsense model of cornstarch is a reservoir: eat D grams, consume it at the
deficit rate, and it lasts D/deficit. That prediction (`t_reservoir_h`) is **wrong**, and
it is more wrong the larger the body.

The reason is that a first-order release depot supplies at `k_dis·A(t)` and this value
decays. Cover ends **the moment the release rate drops below the deficit rate**, with
starch still left in the gut:

```
    release(t) = k_dis · D · F · exp(−k_dis·t)  >  deficit
    ⟹  t_cover = (1/k_dis) · ln( k_dis · D · F / deficit )      ← logarithmic in dose
```

The time gained by doubling the dose is therefore **fixed at ln2/k_dis** — independently
of the starting dose and of body weight. The simulation confirms this:

| Patient | Requirement (mg/kg/min) | Residual EGP | Reservoir prediction | **Actual** | On doubling the dose | **Increment** | Increment on switching to extended-release |
|---|---|---|---|---|---|---|---|
| 6 months, 7 kg | 6.38 | 0.94 | 4.16 h | **3.10 h** | 4.74 h | **+1.64 h** | +0.04 h |
| 1 year, 10 kg | 6.09 | 0.89 | 4.35 h | **3.22 h** | 4.84 h | **+1.62 h** | +0.06 h |
| 5 years, 18 kg | 4.35 | 0.63 | 5.98 h | **4.06 h** | 5.68 h | **+1.62 h** | +0.44 h |
| 14 years, 50 kg | 2.72 | 0.39 | 8.60 h | **5.18 h** | 6.78 h | **+1.60 h** | +0.90 h |
| 30 years, 70 kg | 2.16 | 0.31 | 10.81 h | **5.92 h** | 7.54 h | **+1.62 h** | +1.18 h |

Across five patients whose weights differ 10-fold the increment is effectively identical
at **1.60–1.64 hours**, and it agrees with the hand-solved prediction
`ln2/0.45 = 1.54 hours`. This is why "just give more starch" does not work, and why the
night is bridged by **slower starch rather than more starch**. As the last column shows,
however, the benefit of extended release depends strongly on body size — effectively zero
in infants (+0.04 h) and meaningful only in adults (+1.18 h). The price of a lower initial
release rate is one the infant's high deficit rate cannot absorb, which is consistent with
extended-release starch not being recommended in infants.

---

## Result 2 — blood glucose level and neuroglycopenia are <ins>separated by lactate</ins>

The brain oxidises lactate. At a blood lactate of 7 mmol/L, lactate supplies about 55–60%
of cerebral oxidative metabolism (van Hall 2009, Boumezbeur 2010). The model sets out the
cerebral fuel budget explicitly and back-calculates **the blood glucose giving the same
cerebral ATP flux** (the isofuel glucose):

| Blood lactate (mmol/L) | 0.8 | 1.5 | 2.5 | 4.0 | 6.0 | 8.0 | 10.0 |
|---|---|---|---|---|---|---|---|
| Isofuel glucose (mmol/L) | **2.80** | 2.41 | 1.99 | 1.56 | **1.19** | 0.95 | 0.77 |
| Lactate fraction of cerebral fuel | 15 % | 25 % | 35 % | 45 % | 54 % | 60 % | 64 % |

The 2.8 mmol/L at which a normal child convulses corresponds, at a lactate of 6 mmol/L, to
1.2 mmol/L. This is why a GSD Ia child playing quite happily at a blood glucose of 1.5–2.0
is not an observational error.

### And here is where the clinical trap appears

**Normalise the lactate without raising the glucose and the threshold comes back.** The
model computes this — holding glucose fixed at 1.9 mmol/L:

| | lactate 6.5 mmol/L | lactate 1.6 mmol/L |
|---|---|---|
| Fuel adequacy index | 0.881 | **0.800** |
| Symptoms? | no | **yes** |

Convulsions can appear at a blood glucose the patient tolerated for years. This is
quantitative grounds for why judging that control has been "improved" from the biochemical
markers alone is dangerous.

---

## Result 3 — the counter-regulatory response is not merely futile but <ins>harmful</ins>

Glucagon challenge test (a 20-fold rise in glucagon in the model, 60 minutes):

| | Change in glucose | Change in lactate |
|---|---|---|
| Normal control | **+0.34 mmol/L** | +0.04 mmol/L |
| GSD Ia | **+0.00 mmol/L** | **+5.29 mmol/L** |

Glucagon mobilises glycogen **down a blocked route**. Not one µmol of glucose moves and the
entire amount comes out as lactate. The lactate cost per 1 mmol/L of glucose gained is not
"very large" but **undefined** — there is no numerator to divide by. This is the structure
by which hypoglycaemia turns itself into lactic acidosis, and it is why emergency glucagon
injection in GSD I is not merely pointless but harmful.

---

## Result 4 — for restored activity the two markers have <ins>different functional forms</ins>

When gene therapy restores activity a (14 years, 50 kg):

| a | Fasting tolerance | Daily starch requirement | Lactate | Urate |
|---|---|---|---|---|
| 0 | 5.18 h | 178 g | 3.64 | 7.11 |
| 0.02 | 5.74 h | 164 g | 2.94 | 6.27 |
| 0.05 | 6.90 h | 149 g | 2.34 | 5.53 |
| **0.08** | **8.58 h** | 138 g | 1.98 | 5.05 |
| 0.12 | 12.10 h | 128 g | 1.68 | 4.59 |
| 0.25 | 22.42 h | 111 g | 1.21 | 3.74 |

The starch requirement falls **linearly** in a (178 → 111 g), while fasting tolerance has
the form `1/(deficit − a·Vmax)` and therefore diverges **hyperbolically** (5.2 → 22.4 h).
Because the two curves have different shapes:

- **the critical activity for an unbroken 8-hour night is a\* ≈ 0.08 (8 %)**, and
- at that point the reduction in starch is only **22.5 %**.

Read the other way round it matters more: **below a\* the patient's night is unchanged even
as the biochemical markers and the starch dose improve.** This is exactly where the fact
that DTX401's primary endpoint is the percentage reduction in starch intake meshes with the
model, and the model specifies numerically from where improvement in a marker may be read
as improvement in a life.

---

## Result 5 — the durability of AAV in children is determined by <ins>hepatic growth, not silencing</ins>

AAV episomes do not replicate. A growing liver dilutes them. The same delivered dose
(a₀ = 0.22) given at different ages:

| Age at dosing | Remaining hepatic growth | Peak activity | Activity at 10 years | **Retention** |
|---|---|---|---|---|
| 2 years | ×4.31 | 0.210 | 0.035 | **16.5 %** |
| 6 years | ×3.10 | 0.213 | 0.047 | **22.2 %** |
| 12 years | ×1.70 | 0.217 | 0.082 | **37.8 %** |
| 18 years | ×1.12 | 0.219 | 0.118 | **53.9 %** |
| 30 years | ×1.00 | 0.219 | 0.130 | **59.4 %** |

**Peak activity is the same whatever age it is given at.** All that differs is retention,
and that is set by the remaining hepatic growth factor, not by promoter silencing. And
because anti-capsid neutralising antibodies make re-dosing impossible, **the age at dosing
becomes a decision more important than the dose and irreversible**. This is where the
logical requirement for a redosable LNP-mRNA (mRNA-3745) as an alternative to AAV in
children who have not finished growing comes from — mRNA is not diluted.

---

## Result 6 — urate is a quotient, so the two levers <ins>multiply rather than add</ins>

```
    urate = production × (1 − E_allopurinol)
            ─────────────────────────────────────────
            clearance ÷ (1 + lactate/Ki + ketone/Ki′)
```

> **The author's prior hypothesis was refuted.** Because lactate enters the denominator
> while also contributing to the numerator via the pentose phosphate pathway, dietary
> control was expected to overwhelm allopurinol. The simulation answered otherwise — the
> single-agent effects are almost identical.

| | Mean lactate | Urate | vs baseline |
|---|---|---|---|
| Poorly controlled baseline | 4.48 | 14.47 mg/dL | — |
| Allopurinol 300 mg alone | 4.48 | 6.78 | **−7.68** |
| Intensified diet alone | 3.91 | 6.56 | **−7.90** |
| Both | 3.91 | 3.02 | **−11.44** |

What was confirmed is the structural prediction:

- additive prediction **−15.58** → off by 4.14 mg/dL from the actual
- multiplicative prediction **−11.38** → agrees with the actual **−11.44** to within
  0.06 mg/dL

**Add allopurinol to a patient whose diet has just been intensified and it falls far less
than the sum of the two effects, and reading that shortfall as non-adherence is wrong.**

---

## Result 7 — the hypoketotic signature <ins>emerges without being instructed</ins> (a falsifiable confirmation)

Hypoketotic hypoglycaemia is the bedside discriminator separating GSD I from GSD 0/III/VI.
The model has no term saying "lower the ketones". Because ChREBP is switched on **by G6P**
and not by insulin, malonyl-CoA stays high even during fasting, and CPT-1 is simply
inhibited. After a 6-hour fast:

| | Glucose | 3-OHB | Free fatty acids | Lactate | malonyl-CoA |
|---|---|---|---|---|---|
| Normal control | 3.50 | 0.357 | 0.30 | 0.58 | 0.55 |
| GSD Ia | **1.43** | **0.271** | **0.60** | **8.89** | **2.40** |

GSD Ia has ketones that are **lower** despite a glucose more than 2 mmol/L lower and
lipolysis twice as high. The ratio of ketones to lipolysis is 1.19 in the normal control
against **0.45** in GSD Ia — a 2.6-fold difference. This is the model's strongest
validation item. Had this value come out the other way round, the entire lipid module would
have had to be discarded.

---

## Result 8 — 30 years: a borderline difference in markers is <ins>not a borderline difference in outcome</ins>

An 8-year-old, 26 kg patient, on two regimens with similar total daily carbohydrate but
different **distribution** (continuous vs three meals a day plus a long night):

| | Good control | Poor control | Ratio |
|---|---|---|---|
| Mean lactate | 3.87 | 4.60 | ×1.19 |
| Glucose nadir | 4.55 | 1.33 | — |
| Urate | 6.66 | 14.03 | ×2.1 |
| **Hepatic adenoma burden (30 years)** | **0.023** | **0.563** | **×24.5** |
| Cumulative HCC risk | 0.011 | 0.681 | ×65 |
| UACR (30 years) | 200 | 373 mg/g | ×1.9 |
| Height SDS | −1.21 | −2.99 | — |
| Lumbar spine BMD Z | −1.34 | −2.25 | — |

**A 1.19-fold difference in mean lactate becomes a 24.5-fold difference in adenoma
burden.** Because the metabolic control index enters with an exponent of 2.6 and the
adenoma term is itself self-amplifying. A biochemical difference one would call
"borderline" in clinic is not borderline in the outcome it predicts.

Adding an ACE inhibitor reduces the 30-year UACR from 200 to 69 mg/g, a **65.5 %**
reduction (the same direction and magnitude as the observation of Martens 2009).

---

## Result 9 — GSD Ib: a drug that looks like a switch

In SLC37A4 (G6PT) deficiency the catalytic subunit is intact but the substrate cannot enter
the endoplasmic reticulum. On top of that a lesion specific to neutrophils arises — when
1,5-anhydroglucitol is phosphorylated by hexokinase, **only G6PT can export that phosphate
ester**, and in Ib there is no exit.

| | Baseline | 4 weeks | 12 weeks | 12 months |
|---|---|---|---|---|
| 1,5-AG (µg/mL) | 14.3 | 4.7 | 4.7 | 4.7 |
| Neutrophil 1,5-AG6P | 1.54 | 0.49 | 0.49 | 0.49 |
| **ANC (×10⁹/L)** | **0.79** | **3.54** | **3.54** | **3.54** |

Because marrow suppression enters with a Hill coefficient of 3.2, the ANC response looks
like a **threshold** rather than a dose-response — once 1,5-AG6P falls past the critical
value, recovery is almost complete. This is consistent with the clinical observations of
Wortmann 2020 NEJM and Grünert 2022 (n = 112).

---

## Result 10 — head-to-head comparison of overnight regimens (5 years, 18 kg, a 9-hour night)

| Regimen | Glucose nadir | Time <3.9 | Time <3.0 | Time below fuel index | Mean lactate | Peak lactate |
|---|---|---|---|---|---|---|
| Uncooked starch 1.6 g/kg q4h | 3.32 | 1.54 h | 0.00 h | 0 h | 2.72 | 3.45 |
| Uncooked starch 1.6 g/kg q6h | 2.51 | 3.10 h | 1.54 h | 0 h | 3.94 | 6.99 |
| Uncooked starch 2.4 g/kg once | 1.98 | 4.87 h | 3.39 h | 0 h | 4.06 | 7.99 |
| Extended-release 2.0 g/kg once | 2.41 | 4.95 h | 2.51 h | 0 h | 3.88 | 7.49 |
| Continuous infusion 7 mg/kg/min | **4.60** | **0.00 h** | 0.00 h | 0 h | **2.21** | 3.05 |
| Pump stops at hour 5 during infusion | **1.42** | 3.31 h | **3.07 h** | 0 h | 3.24 | 7.13 |

Two things must be read together.

1. **The time below the fuel index is zero on every regimen.** Even on regimens where the
   glucose went below 2.0 — because the lactate rose correspondingly and fed the brain
   instead. Result 2 reappears here.
2. **The regimen that is safest while running is the most vulnerable when it stops.**
   Continuous infusion keeps insulin sustained so the liver mobilises nothing, and hence
   the glucose nadir after a stop is 1.42, the lowest of any regimen. Safety and
   vulnerability come from the same mechanism.
   (The intercurrent-illness + vomiting scenario has the same structure: 24 hours without
   rescue gives a glucose of 0.24, lactate 7.03, HCO₃⁻ 16.4 / starting IV glucose at
   8 mg/kg/min at 12 hours gives glucose 8.37, lactate 3.15, HCO₃⁻ 21.4.)

---

## Result 11 — overtreatment is not free

A normal liver exports surplus glucose again. A GSD Ia liver cannot — all the surplus is
trapped as G6P and leaves as glycogen, lactate and fat. An 8-year-old, 26 kg, varying the
delivery ratio relative to requirement:

| Delivery ratio | Glucose | Lactate | Insulin | Glycogen | Urate | Liver volume |
|---|---|---|---|---|---|---|
| 0.60 | 3.58 | 4.07 | 69 | 51 mg/g | 7.75 | ×1.33 |
| 0.75 | 4.28 | 4.05 | 97 | 85 | 6.94 | ×1.49 |
| **0.90** | **4.86** | **4.37** | 122 | 109 | **6.80** | ×1.76 |
| 1.05 | 5.38 | 4.94 | 144 | 126 | 7.04 | ×1.97 |
| 1.25 | 6.01 | 5.76 | 170 | 143 | 7.66 | ×2.17 |
| 1.45 | 6.62 | 6.66 | 193 | 155 | 8.38 | ×2.31 |

Urate traces a **U** (minimum at 0.90). Lactate and liver volume increase monotonically
through the excess range. This is the mechanistic explanation for why giving more
carbohydrate in pursuit of euglycaemia is explicitly warned against in the guidelines, and
it shows that the optimum is not "wherever glucose is highest".

---

## Model structure

**47 ODEs.** Intestinal starch delivery (3) · systemic glucose (1) · hepatic carbon (3) ·
lactate / acid-base (3) · lipid (5) · hormones (6) · purine / urate (5) · kidney (5) ·
long-term liver (2) · growth / bone (2) · exposure integrals (3) · gene / mRNA therapy (4) ·
GSD Ib (4) · adjunctive drugs (1)

**Drug PK/PD.** Uncooked cornstarch and extended-release waxy maize starch (first-order
release, `k_dis` 0.45 vs 0.28/h) · continuous nasogastric infusion ·
allopurinol → oxypurinol (t½ 22 h) → xanthine oxidase inhibition · ACE inhibitor ·
statin / fenofibrate · potassium citrate · empagliflozin (SGLT2 → urinary 1,5-AG
excretion) · G-CSF · AAV8-G6PC · LNP-mRNA · IV glucose

**Fourteen scenarios** (the `SCENARIOS` list) — natural history, uncooked starch q4h/q6h,
extended release, continuous infusion, pump stop, intercurrent illness ± rescue, urate,
lipids, kidney, gene therapy, mRNA, empagliflozin, overtreatment, 30-year course, activity
dose-response.

**Eleven Shiny tabs** — patient profile · dietary PK · glucose and the safe time window ·
cerebral fuel budget · the G6P branch point · metabolic biomarkers · liver and long-term
prognosis · kidney · gene / mRNA therapy · GSD Ib · scenario comparison.

### Running

```r
source("gsd1a_mrgsolve_model.R")
scn_uccs(PATIENTS$child_5y, dose_gkg = 1.6, interval = 4) |> plot()
shiny::runApp("gsd1a_shiny_app_en.R")
```

```bash
python3 gsd1a_reference_model.py       # regenerates every number above
dot -Tsvg gsd1a_qsp_model.dot -o gsd1a_qsp_model.svg
dot -Tpng -Gdpi=150 gsd1a_qsp_model.dot -o gsd1a_qsp_model.png
```

---

## Defects caught by the second implementation

Left here for documentation. All of them were in the first implementation and came to light
in the course of comparison with the Python re-implementation.

1. **The glycogenolysis flux had been multiplied by 60** — while making mmol/h, minutes were
   multiplied in again on a constant that was already in hours. Glycogen collapsed to
   0.3 mg/g at steady state, i.e. the model had abolished glycogen in a "glycogen storage
   disease".
2. **The plasma triglyceride volume of distribution was wrong by 10-fold** (`Vpl×10` instead
   of `Vpl×100`). The normal control's triglycerides came out at 13 mg/dL.
3. **Glycolytic pyruvate was being double counted** — the same carbon was exported as lactate
   and also sent to lipogenesis. Lactate and triglycerides inflated together. Resolved by
   branching with `fLacOut`.
4. **Glycogen synthase had no insulin arm** — with only G6P allosteric activation, a normal
   liver, in which G6P is low, could not store glycogen at all. Rewritten as the product of
   two arms, covalent (insulin) × allosteric (G6P).
5. **The lactate crisis had no brake.** The term by which acidosis inhibits PFK-1 was
   missing, so lactate went past 20 mmol/L up to the ceiling. This term is not a
   convenience; it is the only negative feedback that actually ends the crisis.
6. **The phosphorylase operating range was too narrow.** The Hill constants had been matched
   to resting hormone concentrations, so the enzyme was already near its ceiling, and the
   glucagon challenge test raised glucose by only 0.01 mmol/L even in the normal control.
   Result 3 itself did not hold.
7. **Three long-term states diverged.** Adenoma burden reached 57,041 at 30 years (fifty
   thousand times the liver volume) and height SDS −8.05. Each was rewritten with a logistic
   ceiling and with "relaxation towards a target + a gate at the time of skeletal
   maturation" respectively. A number a human cannot have is not a result just because the
   model emitted it.
8. **1,5-anhydroglucitol, a concentration, had been made proportional to body weight.** The
   whole GSD Ib arm died, so the ANC was always 0.

And **one refuted hypothesis**: in Result 6 the prior expectation that "lactate control
overwhelms allopurinol" was wrong. The single-agent effects were almost identical. What was
right instead was the **subadditivity** prediction that follows from the multiplicative
structure, accurate to 0.06 mg/dL.

---

## Limitations

Set down honestly.

- **Residual EGP is lower than in the literature.** The model gives about 15% of
  requirement (0.89 at 1 year, 0.31 mg/kg/min in adults). Isotope studies (Tsalikian 1984,
  Powell 1981) suggest 25–35% of normal. Raising residual EGP raises the fasting glucose
  nadir to 2.5–3.0 mmol/L, which is far too high for untreated GSD Ia. The two observations
  could not be satisfied at once, and the choice was made to match the fasting trajectory.
  This tension remains unresolved.
- **At chronic steady state the relationship between liver volume and control status is the
  opposite of the clinical one.** In the model, poor control (pulsatile delivery) lowers mean
  glycogen and so the liver is actually smaller (×1.43 vs ×1.83). In reality the liver of a
  poorly controlled patient is larger. This appears to be a limitation of approximating
  daytime over-intake as sinusoidal delivery.
- **The metabolic control index (MCI) is a construct, not a measured quantity.** The exponent
  of 2.6 was matched to epidemiological observations of adenoma incidence, not measured
  independently. The 24.5-fold factor in Result 8 is sensitive to this exponent.
- **The normal control group is a comparator, not a precise model.** Its 14-hour fasting
  glucose of 2.9–3.2 mmol/L is lower than in real healthy children.
- The kidney, growth and bone modules are semi-quantitative. Only the direction and rough
  magnitude of the time course can be trusted.
- G6PC3 deficiency, GSD Ic/Id, pregnancy and the course after liver transplantation are not
  in the model.

---

## ⚠️ Disclaimer

This model is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was constructed from the public literature and clinical trial data but has
not been independently validated or certified, and **must not be used directly for real
clinical decision-making, prescribing, or regulatory submission.** The parameters and
assumptions are illustrative approximations, and fitting and validation against real
patient data are separately required.
