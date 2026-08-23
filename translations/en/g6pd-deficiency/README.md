# G6PD Deficiency — QSP Model

**Glucose-6-Phosphate Dehydrogenase Deficiency** · drug- and fava-bean-induced acute
haemolysis

<a href="g6pd_qsp_model.svg"><img src="g6pd_qsp_model.png" width="820" alt="G6PD QSP mechanistic map"></a>

---

## In one line

G6PD deficiency is usually written down as a **number** — "activity below 10% of
normal". This model claims that the number is the wrong object, and that this disease
is a **function**. The mature erythrocyte has neither nucleus nor ribosome and cannot
make new enzyme. The enzyme is loaded **once**, in the erythroblast, and after that it
only declines:

> ### E(a) = E₀ · exp(−ln2 · a / τ)
> **The variant is τ, the slope of this curve.** Not the intercept.

| Variant | E₀ | τ (days) | Activity left in a 115-day-old cell | What the assay reports |
|---|---|---|---|---|
| B (normal) | 1.00 | 62 | **27.6%** | 100% |
| A− (Class III) | 0.55 | 13 | **0.12%** | 16.5% |
| Mediterranean (Class II) | 0.05 | 12 | **0.0064%** | 2.1% |

The third column is the most important number in this repository. And the fourth
column — what the laboratory actually reports — is the **age-weighted mean** of that
curve. A mean hides a slope.

---

## Why this explains the whole disease

Given an oxidant load OX, a cell survives only if its own NADPH regeneration capacity
exceeds it. Because E(a) is a clean exponential, **the critical age solves in closed
form**:

> ### a\* = (τ / ln2) · ln( E₀ · V<sub>max</sub> / OX<sub>eff</sub> )

Note that **there is no fitted haemolysis rate constant** in this expression. All that
is in it is the genotype (E₀, τ), the biochemical capacity (V<sub>max</sub>) and the
exposure (OX).

Haemolysis is not a rate constant multiplying a pool but **a blade passing through the
age histogram**. And everything else follows from that.

---

## Nine results the model generates

Every number is a simulation result. They were computed separately with
`g6pd_mrgsolve_model.R` (mrgsolve) and `g6pd_reference_check.py` (pure Python, an
independent implementation), and the two agree to 2–3 significant figures.

### ① The same drug, the same dose, the same "severe deficiency" label — a different disease

| | Baseline Hb | Nadir Hb | Fall | Day of nadir | Day 60 | a\* | Erythrocytes at risk | Reticulocytes | Haemoglobinuria |
|---|---|---|---|---|---|---|---|---|---|
| **Normal** | 15.00 | 15.00 | 0.0% | — | 15.00 | 120 d | 0% | 1.0% | none |
| **A−** | 15.00 | **11.77** | −21.5% | 5.7 d | **14.73** | 84 d | 30% | 3.3% | none |
| **Mediterranean** | 14.06 | **8.74** | −37.9% | 4.1 d | **12.55** | 36 d | 70% | 8.9% | 0.25 d |

This is the result of giving primaquine 30 mg/day **for the whole 60 days without
stopping**.

A− bottoms out and then **comes back close to baseline while still taking the drug.**
This is exactly the curve Dern observed in 1954, and the model has not a single
parameter that knows about the phenomenon. Because the slope is steep, a\* falls
*inside* the distribution, only the oldest 30% die, every surviving cell is young, and
the reticulocytes newly released are full of enzyme.

In Mediterranean the curve is flat and low. The resistant population it can leave
behind is thin, the new cells the marrow supplies are deficient too and so reach a\*
again within weeks, and a substantial fraction exceed the remaining buffering capacity
so far that they burst **inside the vessel rather than in the spleen**. What is reached
is not recovery but a low plateau propped up by very nearly maximal erythropoiesis.

> Two patients, one drug, one dose, both "severe deficiency". The difference between
> "continue and observe" and "stop and transfuse" is **a single decay constant**.

### ② The half-life is the toxicology

| | Nadir Hb | Fall | a\* | Erythrocytes at risk | Haemoglobinuria |
|---|---|---|---|---|---|
| Tafenoquine 300 mg **single dose** (t½ 15 days) | 9.71 | **−35.3%** | 66 d | 45% | 1.1 d |
| Primaquine 15 mg × 14 days (t½ 7 hours, **30% less total dose**) | 12.96 | −13.6% | 97 d | 19% | none |

A reticulocyte rescue call takes 4–7 days to arrive. A drug that switches off quickly
strikes at a rate the system can answer, and if it cannot answer, **stopping actually
works.** A drug that switches off slowly outruns the rescue, and **there is nothing to
stop.** This is the mechanistic content of why the FDA requires quantitative G6PD
testing (above 70% of normal) only for tafenoquine, and permits weekly primaquine in
mild deficiency without testing.

### ③ Once-weekly dosing is not a dose reduction but a different mechanism

Applying the WHO-recommended regimen (primaquine 45 mg **once weekly** × 8 weeks) to
an A− patient: nadir 12.45 (**−17.0%**), and after that merely a **sawtooth** that
never gets any deeper. Each pulse shaves off only the thin old layer, and the 7-day
interval gives the reticulocytes time to fill that place.

One honest clue the model volunteers: **the first pulse itself costs very nearly as
much as daily dosing.** What the weekly regimen buys is not a gentle start but **the
absence of accumulation**.

### ④ The methylene blue contraindication is not a fact to memorise but a result of the calculation

**One** NADPH pool feeds two consumers — glutathione reductase (oxidative defence)
and the NADPH-methaemoglobin reductase that converts methylene blue into leucomethylene
blue, the species that actually does the reducing. Raise methaemoglobin with dapsone
100 mg/day and then give the antidote on day 30:

| Time after the dose | → | 0 h | 2.4 h | 12 h | 24 h |
|---|---|---|---|---|---|
| **Normal** MetHb % | | 7.2 | **2.2** | 4.4 | 5.5 |
| **A−** MetHb % | | 7.9 | 6.6 | **8.3** | 7.2 |
| **A−** oxidant load | | 0.18 | **0.68** | 0.34 | 0.19 |

In the normal subject it collapses within 2.4 hours. In A− it barely moves, and **at
the same time the oxidant load jumps 3.8-fold** — because the dye that was not reduced
is itself a fresh oxidant. No parameter says "methylene blue is contraindicated in
G6PD deficiency". It simply computes that way.

### ⑤ The oxidant dose of rasburicase is set by the tumour, not by the prescription

Urate oxidase runs `urate + O₂ + H₂O → allantoin + H₂O₂` **mole for mole**. So the
peroxide load is stoichiometrically proportional to the urate pool that is destroyed.
The same 0.2 mg/kg:

| | Nadir Hb | Fall | Peak MetHb |
|---|---|---|---|
| Normal, urate 20 mg/dL | 15.00 | 0.0% | 1.0% |
| A−, urate 20 mg/dL | 10.32 | −31.2% | 5.0% |
| Mediterranean, urate 20 mg/dL | 7.09 | **−49.6%** | 5.8% |
| Mediterranean, urate **8** mg/dL | 9.86 | −29.9% | 4.2% |
| Mediterranean, **allopurinol** (no uricase) | 14.05 | −0.1% | 1.3% |

Two axes have to be read separately. **The genotype axis**: the same prescription,
three outcomes. **The tumour axis**: the same prescription and the same genotype, and
the nadir moves by 20 percentage points. Which is why, when the G6PD status is unknown,
the rational manoeuvre is not a dose reduction but **allopurinol, which never makes the
peroxide at all**.

### ⑥ The diagnostic assay lies at precisely the moment it is needed

The quantitative assay reports the **age-weighted mean over the whole** circulating
erythrocyte population. But the haemolysis has just erased the cells that were old and
enzyme-free, and their place has been taken by reticulocytes carrying **1.5 times** the
activity of a young cell.

A genuine A− patient: **16.5%** before exposure → **29.5%** during haemolysis. Test at
the nadir and the patient looks almost twice as normal as they are. A normal result
during an acute episode does not exclude deficiency. **The retest belongs three months
later** — one turn of the erythrocyte lifespan.

Give tafenoquine 300 mg to a heterozygous woman (60% normal cells) and the whole-blood
assay reads **60.9%**. That is below the FDA's 70% threshold, but it is a value many
laboratories will pass as normal, and it is far above the 30% line. And yet her
haemoglobin falls by **22.9%**. Because 40% of her erythrocytes have a\* = 19 days.

### ⑦ Neonatal kernicterus is multiplication, not addition

A 3 kg term infant, from 12 hours of age to day 12:

| | Peak total bilirubin (mg/dL) | Peak **free** bilirubin (nM) |
|---|---|---|
| Neither | 4.62 | 3.0 |
| G6PD deficiency only | 15.57 | 20.7 |
| UGT1A1 (TA)7/7 only | 8.58 | 6.9 |
| **Both** | **23.15** | **104.4** |

On the total-bilirubin scale it looks almost additive (predicted 19.5 against 23.15
observed). But on the **free bilirubin** scale — the species that actually crosses the
blood-brain barrier — it is **4.2 times** the additive prediction (predicted 24.6
against 104.4 observed). Because albumin binding saturates.

> **The non-linearity lives at exactly the site where the damage happens.** A model
> that tracks only total bilirubin misses it. G6PD deficiency is a leading cause of
> kernicterus worldwide not because the deficiency is especially powerful in itself
> but because it **multiplies** with a common second defect.

### ⑧ Self-limitation is a property neither of the drug nor of the variant

Repeating scenario 1 during a parvovirus B19 aplastic crisis (marrow reserve 15%):

| | Baseline | Nadir | Fall | Day 60 |
|---|---|---|---|---|
| A−, normal marrow | 15.00 | 11.77 | −21.5% | **14.73** |
| A−, aplastic crisis | 10.65 | 8.24 | −22.7% | **10.21** |

**The depth of the nadir is very nearly the same** — a\* knows nothing about the
marrow. What has gone is the recovery. Self-limitation was **what the reticulocyte arm
was answering** all along.

### ⑨ The enzyme that makes the oxidant also makes the cure

Primaquine 30 mg/day in an A− patient who is a CYP2D6 poor metaboliser: haemoglobin
15.00 unchanged, a\* 120 days, erythrocytes at risk 0%. **And yet the assay still
reads 16.5%** — the genotype is unchanged and only the bioactivation has gone. Because
the very enzyme that turns primaquine into an oxidant turns primaquine into a
radical-cure agent, this patient escapes the haemolysis and **fails the radical cure at
the same time.**

---

## Files

| File | Contents |
|---|---|
| [`g6pd_qsp_model.dot`](g6pd_qsp_model.dot) | Source of the mechanistic map — **176 mechanistic nodes · 234 edges · 17 clusters** (plus 64 invisible layout-only edges) |
| [`g6pd_qsp_model.svg`](g6pd_qsp_model.svg) / [`.png`](g6pd_qsp_model.png) | The renders (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`g6pd_mrgsolve_model.R`](g6pd_mrgsolve_model.R) | **48-ODE** mrgsolve model + variant/patient generators + **22 scenarios** |
| [`g6pd_reference_check.py`](../../../g6pd-deficiency/g6pd_reference_check.py) | An independent re-implementation (pure-Python RK4). Cross-validates the numbers above |
| [`g6pd_reference_output.txt`](../../../g6pd-deficiency/g6pd_reference_output.txt) | The output of that run, verbatim |
| [`g6pd_shiny_app.R`](g6pd_shiny_app.R) | 8-tab interactive dashboard |
| [`g6pd_references.md`](g6pd_references.md) | **77 references**, 15 sections, PubMed links |

### Running it

```bash
# redraw the mechanistic map
dot -Tsvg g6pd_qsp_model.dot -o g6pd_qsp_model.svg
dot -Tpng -Gdpi=150 g6pd_qsp_model.dot -o g6pd_qsp_model.png

# the 22 scenarios (R + mrgsolve)
Rscript -e 'source("g6pd_mrgsolve_model.R"); print(g6pd_age_table()); print(g6pd_all()$summary)'

# cross-validate the same thing with no dependencies at all (Python is enough, about 6 min)
python3 g6pd_reference_check.py --curves

# the interactive dashboard
Rscript -e 'shiny::runApp("g6pd_shiny_app.R")'
```

---

## Model structure

### The mechanistic map (17 clusters)

1. Genotype · WHO variant classification (where τ is set)
2. The pentose phosphate pathway — the erythrocyte's only source of NADPH
3. The glutathione redox network
4. **Erythrocyte age structure — the thesis of this model**
5. Sources of oxidative stress (every trigger converges on a single OX)
6. The haemoglobin oxidation chain (metHb → hemichrome → Heinz bodies)
7. Membrane damage · band 3 clustering · eryptosis
8. Extravascular haemolysis (spleen, slow, the A− route)
9. Intravascular haemolysis (free Hb, haptoglobin, AKI — the fava bean/rasburicase route)
10. The erythropoietic rescue call (a 4–7 day delay)
11. Methaemoglobin and the methylene blue trap (the second NADPH consumer)
12. Bilirubin · neonatal kernicterus
13. Pharmacokinetics (the half-life is the toxicology)
14. Clinical measures
15. Diagnosis (where the assay lies)
16. Management
17. Legend

### The ODE model (48 compartments)

| Compartment group | Number | Contents |
|---|---|---|
| Drug PK | 12 | primaquine (+ oxidant equivalents) · tafenoquine · dapsone · rasburicase · urate · methylene blue · divicine |
| **Erythrocyte age bins** | **12** | R1–R12, each 10 days wide, a 120-day Erlang chain |
| **Damage bins** | **12** | Z1–Z12, the hemichrome/Heinz load per age bin |
| Erythropoiesis | 4 | three marrow transit compartments + circulating reticulocytes |
| Plasma measures | 8 | metHb · bilirubin · free Hb · haptoglobin · tubular injury · creatinine · cumulative haemolysis · cumulative filtered Hb |

Two central structural choices:

- **Glutathione was reduced to a quasi-steady-state algebraic expression.** GSH
  turnover is on a timescale of minutes and the disease on a timescale of days, so this
  is a safe reduction, and it is what allows a\* to come out in closed form.
- **Intravascular haemolysis is governed by a ratio and not by a difference.**
  `ω = OX / (V_max·E(a) + OX_ref)` — a cell whose buffering capacity is effectively
  zero bursts in situ, before the spleen can deal with it. This is what generates the
  fact that Mediterranean has haemoglobinuria and A− does not.

---

## Verification

Two independent implementations — R/mrgsolve (LSODA) and pure Python (fixed-step
RK4) — were written from the same equations in different languages, with different
integrators and different state orderings, and they agree to 2–3 significant figures
across all 22 scenarios. That is not agreement by chance.

In the parameter calibration, **the only freely fitted value is the oxidant scale
factor**, and even that was anchored to a single curve, Dern 1954's A− primaquine
curve. Everything else is a literature central value. The following three came out
right **without being fitted**:

- The assay value in the untreated state falls into each variant's published activity
  class (B 100% · A− 16.5% · Mediterranean 2.1% · heterozygous woman 60.9%).
- Class I variants show a mild chronic non-spherocytic haemolytic anaemia at baseline.
- Mediterranean is slightly low at baseline, Hb 14.06, even with no drug (the oldest
  bin is already being shaved off by the baseline oxidant load alone).

---

## Honest limitations

1. **No fitting was done.** The values were fixed at literature central values without
   formal estimation or validation against individual patient data. This is a
   moderate-fidelity model for teaching and hypothesis generation.
2. **The age structure is a 12-compartment Erlang approximation.** It is broader than
   the real lifespan distribution (CV about 29%), which makes a\* a narrow band rather
   than a blade, and age differences finer than 10 days cannot be resolved.
3. **The total bilirubin of adult haemolytic jaundice comes out somewhat high.** In
   scenario 1 it is 7.2 mg/dL, roughly 1.5–2 times the real 3–5 mg/dL. Because the
   model does not carry the increase in hepatic conjugation reserve with load.
4. **The rebound after recovery is excessive.** In several scenarios the haemoglobin at
   day 60 exceeds baseline by more than 10% (for example 16.76 in 4a). EPO suppression
   was put in only linearly, and this is larger than the real rebound.
5. **The Class I (CNSHA) phenotype is milder than in reality.** It does not reproduce
   chronic transfusion dependence.
6. **The heterozygous woman is run twice and mixed** (`mix_mosaic()`). That is exact
   for the erythrocyte populations, but the plasma measures are a mass-weighted
   recombination and therefore an approximation.
7. **Nothing can be said about the minute-scale dynamics of glutathione**
   (the quasi-steady-state reduction).
8. In the extreme scenarios (for example favism −65.6%, rasburicase −49.6%) the model
   knows nothing of transfusion or death, so it simply proceeds arithmetically. In that
   range only the direction and the ordering of the absolute values should be read.

---

## ⚠️ Disclaimer

This is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was built from the public literature but has not been independently
validated or certified, and **must not be used directly for real clinical
decision-making, prescribing, or regulatory submission.**

For the supporting literature see [`g6pd_references.md`](g6pd_references.md) (77 references, PubMed links).
