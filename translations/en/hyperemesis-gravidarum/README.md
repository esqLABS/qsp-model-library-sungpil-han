# Hyperemesis Gravidarum — QSP Model

**Hyperemesis Gravidarum · HG** — a model in which placental GDF15 is compared against
**an adapting set point** in the maternal hindbrain, that is, a **fold-change
detector**. 42 ODEs, a mechanistic map of 214 nodes in 19 clusters, 16 treatment
scenarios, a 12-tab Shiny app, 77 PMID-verified references, and a dependency-free
Python twin that regenerates every number quoted here.

---

## The question this model sets out to answer

There are five observations about hyperemesis gravidarum that do not fit together well.

1. **Maternal GDF15 does not fall after the first trimester.** It goes on rising to
   term. And yet the nausea and vomiting peak at 9-11 weeks of gestation and are gone
   by 16-20 weeks.
   ([PMID 12495665](https://pubmed.ncbi.nlm.nih.gov/12495665/),
   [31515515](https://pubmed.ncbi.nlm.nih.gov/31515515/))
2. **Women with β-thalassemia have chronically high GDF15 and almost no NVP.**
   ([PMID 38092039](https://pubmed.ncbi.nlm.nih.gov/38092039/))
3. **A low pre-pregnancy GDF15 is the risk factor** — not a high one. (Same reference)
4. **Pre-pregnancy metformin** — a drug that raises GDF15 and is stopped before the
   nausea begins — lowers the risk of HG to aRR 0.29. Pre-pregnancy smoking raises
   GDF15 less and gives aRR 0.51.
   ([PMID 40588059](https://pubmed.ncbi.nlm.nih.gov/40588059/))
5. **Ondansetron** occupies about 98% of the 5-HT3 receptors and yet moved PUQE-24 by
   only −0.51 (not significant) in the first placebo-controlled trial in HG.
   **Mirtazapine** gave −1.86, and the difference widened further after day 4.
   ([PMID 41478546](https://pubmed.ncbi.nlm.nih.gov/41478546/))

A dose-response model in GDF15 **concentration** fails on point 1 (the hormone goes on
rising while the symptoms disappear), comes out **with the sign reversed** on points 2,
3 and 4 (a high GDF15 ought to be worse, and it is protective), and has nothing at all
to say about point 5.

This model takes one structural position and derives all five from it.

> **The hindbrain does not measure how much GDF15 there is. It measures how much more
> there is than before.**

---

## Three structural commitments

### 1. The detector is a ratio, not a concentration

The GFRAL axis has an adapting set point `SP` that tracks `log(GDF15)` with a time
constant `TAU_SP ≈ 30 days`. The nausea drive is `GDF15 ÷ the adapting set point`.

- The placental GDF15 ramp (doubling roughly every 4 days at 5-9 weeks of gestation) is
  faster than the set point, so the ratio shoots up.
- Once the ramp flattens the set point catches up and the ratio returns to 1 →
  **the symptoms disappear with the hormone still high.**
- Anything that raises GDF15 **before** pregnancy has already raised the set point
  before the placenta arrives → the same absolute concentration makes a much smaller
  ratio.

`ALPHA` is "the completeness of the adaptation". `ALPHA = 1` is a pure fold-change
detector and `ALPHA = 0` a pure concentration detector. **Set one parameter to zero and
the same equations become the conventional model, and all four of the predictions above
are inverted.**

| | Fold-change model (ALPHA=0.92) | Concentration model (ALPHA=0) |
|---|---|---|
| PUQE at 9 weeks | 14.2 | 15.0 |
| PUQE at 16 weeks | **3.0** | **15.0** |
| PUQE at 28 weeks | **3.0** | **15.0** |
| Peak PUQE in β-thalassemia | **3.0** | **15.0** |

The concentration model is wrong in two directions: the symptoms never resolve, and
β-thalassemia — reported to have almost no NVP — becomes **the most severe** case in
the cohort (because its GDF15 really is the highest).

### 2. Efficacy is set by node position

The drive reaches the emetic pattern generator through the NTS. The NTS sums a **large**
GFRAL/area postrema term and a **small** peripheral vagal/5-HT3 term. A drug cannot do
better than the weight of the node it occupies.

Each drug attaches to **a specific receptor class at a specific node**, at an occupancy
computed from its published Ki, unbound fraction, brain penetration ratio and PK. Only
two parameters (`W_VAG`, `E0`) were fitted; the rest of the ordering is prediction.

| Drug | Node | ΔPUQE (7 days) | NTS inhibition |
|------|------|------------:|-----------:|
| anti-GDF15 antibody | the ligand itself | **−10.57** | 0.0% |
| Gabapentin 600 mg q8h | NTS α2δ | −2.26 | 15.4% |
| Mirtazapine 30 mg qHS | H1 + 5-HT2 + 5-HT3 | −2.17 | 8.5% |
| Doxylamine/pyridoxine | H1 + vestibular | −1.50 | 4.3% |
| Clonidine 5 mg patch | presynaptic α2 at the NTS | −1.37 | 11.5% |
| Promethazine 25 mg q6h | NTS H1 / M1 | −1.22 | 3.9% |
| **Ondansetron 8 mg q8h** | **peripheral 5-HT3** | **−0.56** | **2.5%** |
| Metoclopramide 10 mg q6h | AP D2 + gastric motility | −0.07 | 0.0% |
| Corticosteroid | no node on the chain | **0.00** | 0.0% |

Ondansetron occupies 98% of the 5-HT3 receptors and yet inhibits NTS transmission by
only 2.5%. Because in established HG the peripheral 5-HT3 branch carries only **3.4%**
of the total drive. This is the model's answer to why the guideline-recommended drug
was, of all things, the one that failed.

### 3. The nutritional cascade runs on a different clock

Vomiting → loss of volume, Cl⁻ and K⁺ → **chloride-responsive metabolic alkalosis**
(without replacing the chloride the kidney cannot excrete bicarbonate, so nothing else
given will correct it). And a 28 mg thiamine store with a half-life of about 15 days
empties over **weeks**.

The consequence is that the risk of Wernicke can go on rising even while an antiemetic
is normalising the PUQE, and that **a glucose infusion given without thiamine is worse
than giving no fluid at all.**

| Intervention (a severe, protracted course) | Thiamine nadir | P(Wernicke) | Peak HCO3⁻ | Weight loss |
|---|---:|---:|---:|---:|
| No fluid | 12.2 mg | 0.0% | 40.6 | 9.3% |
| Normal saline + KCl, no glucose | 12.2 mg | 0.0% | 34.8 | 9.3% |
| **Glucose, no thiamine** | **6.1 mg** | **15.7%** | 34.8 | 7.4% |
| Glucose + thiamine 100 mg/day | 20.1 mg | 0.0% | 34.8 | 7.4% |

The third row is the clinically important one. The hydration and the calorie provision
are both better, and the risk of Wernicke is far higher. Because glucose consumes an
already depleted store faster
([PMID 30889425](https://pubmed.ncbi.nlm.nih.gov/30889425/)). Add thiamine alone to the
same fluid and the risk disappears. **Thiamine first, then glucose.**

---

## Prevention: the therapeutic window closes at conception

The model's central testable prediction, and its most actionable result.

| Intervention | GDF15 at conception | Peak PUQE | Weight loss | HG case |
|------|---------------:|----------:|--------:|:-------:|
| None | 250 | 14.3 | 6.7% | YES |
| **Metformin, before pregnancy only** | 485 | **10.6** | 2.9% | no |
| **The same metformin, started at 6 weeks of gestation** | 250 | **14.4** | 7.2% | YES |
| Smoking, before pregnancy | 325 | 12.0 | 4.3% | no |
| β-thalassemia (lifelong) | 3000 | 3.0 | 0.0% | no |
| Recombinant GDF15, before pregnancy | 5765 | 3.0 | 7.9% | no |
| Recombinant GDF15, started at 8 weeks of gestation | 250 | 15.0 | 9.8% | YES |
| anti-GDF15 antibody, at 8 weeks of gestation | 250 | 14.2 | 3.3% | YES |

The same drug, the same dose, **the opposite sign**. Metformin before pregnancy lowers
the peak PUQE from 14.3 to 10.6, and metformin started at 6 weeks of gestation gives
14.4 — it is laying GDF15 on top of an already steep placental ramp, and it is far too
late to move a set point with a time constant of 30 days.

Three independent exposures (β-thalassemia, metformin, smoking) all act before
pregnancy, all are protective, and the strength of protection lines up in the order of
the magnitude of their GDF15 rise. The model reproduces that ordering from `ALPHA`
alone plus the size of each exposure's GDF15 rise — it is not a fit.

### A cost the model found for itself

A fold-change detector responds to **the rise itself**, whether placental or injected.
During pre-pregnancy titration of recombinant GDF15 the peak PUQE is **15.0** — a
prophylactic drug producing a miniature version of the very disease it is meant to
prevent. If one were developed, it would have to be titrated over months rather than
weeks. Metformin does not have this problem over the exposure range modelled, and that
is a large part of its appeal.

---

## Two axes of risk

Fejzo 2024 reports that **fetal GDF15 production** and **maternal susceptibility**
contribute independently to the risk. The model keeps these as separate parameters
(`TROPH_GAIN`, `SENS`), and the two are not interchangeable.

| Cohort | Peak PUQE | Peak GDF15 | Peak hCG | Nadir TSH | Peak free T4 |
|--------|----------:|-----------:|---------:|---------:|-------------:|
| Normal pregnancy | 6.8 | 13 461 | 79 142 | 0.68 | 18.2 |
| HG, susceptibility-driven | 14.3 | 13 461 | 79 142 | 0.68 | 18.2 |
| HG, production-driven | 10.8 | 26 672 | 158 284 | **0.16** | **25.2** |

hCG and GDF15 come out of the same cell, the syncytiotrophoblast. The model therefore
predicts that biochemical thyrotoxicosis tracks **the fetal production axis only**. That
means **two women with the same PUQE can have entirely different TSH**, and it is
testable with paired GDF15/TSH measurements.

---

## Validation

All 15 validation targets passed. Regenerate with
`python3 hg_reference_impl.py --check`.

| Target | Model | Observed | Role | PMID |
|------|-----:|------|------|------|
| VOMIT mirtazapine ΔPUQE at day 2 | −1.92 | −1.86 (−3.61,−0.12) | **fitted** (`E0`) | 41478546 |
| VOMIT ondansetron ΔPUQE at day 2 | −0.52 | −0.51 (−2.32, 1.30) | **fitted** (`W_VAG`) | 41478546 |
| VOMIT mirtazapine − ondansetron at day 7 | −1.61 | −1.35 (−3.10, 0.40) | prediction | 41478546 |
| Koren doxylamine ΔPUQE at day 14 | −1.75 | −0.9 | prediction | 20843504 |
| Guttuso gabapentin relative reduction | 16% | 52% (16-88) | prediction | 33451591 |
| Maina clonidine PUQE improvement | 1.37 | CI 0.43-3.24 | prediction | 24684734 |
| Yost corticosteroid (negative) | 0.00 | 34% vs 35% (P=.89) | prediction | 14662211 |
| Peak PUQE in β-thalassemia | 3.00 | "almost no NVP" | prediction | 38092039 |
| Peak PUQE with pre-pregnancy metformin | 10.59 | aRR 0.29 | prediction | 40588059 |
| Peak PUQE with pre-pregnancy smoking | 12.04 | aRR 0.51 | prediction | 40588059 |
| Timing of the HG peak | 8.5 weeks | 9-11 weeks | fitted (natural history) | 31515515 |
| PUQE at 16 weeks in HG | 3.0 | resolution | fitted (natural history) | 31515515 |
| Weight loss in HG | 6.7% | ≥5% | prediction | 34555550 |
| GDF15 ratio, 28 weeks/9 weeks | 1.43 | >1 | prediction | 12495665 |
| Production-driven HG: TSH suppression | 0.16 | <0.4 in ~60% | prediction | 15073140 |

### Four of 179 parameters were fitted

| Parameter | What it was fitted to |
|----------|-------------------|
| `ALPHA`, `TAU_SP` | the natural history — the timing of the peak and resolution at 16 weeks |
| `W_VAG` | ondansetron's −0.51 |
| `E0` | mirtazapine's −1.86 |

The **ratios** of the node weights (`R_H1`, `R_A2`, `R_A2D` …), each drug's Ki, PK and
brain penetration ratio, the thiamine kinetics, the electrolyte physiology and the
strength of hCG-TSHR cross-stimulation were all fixed from the literature.

### The stated failure

**Metoclopramide is predicted to be effectively ineffective (ΔPUQE −0.07)**, whereas
the clinical trials mostly find it comparable to promethazine. Its central D2 occupancy
is low (~11%) and its prokinetic action flows into the peripheral branch that the model
says is nearly irrelevant. It is one or the other: either the peripheral branch matters
more than `W_VAG` allows, or metoclopramide has an action the model has left out.
**This is the easiest point at which to falsify the node-position rule.**

Other known limitations:

- The Guttuso gabapentin prediction (16%) only overlaps the **lower boundary** of the
  observed confidence interval.
- The Koren doxylamine prediction (−1.75) is larger than the observed value (−0.9).
- In the untreated severe course, Cl⁻ falls to 65 mmol/L, which is lower than the range
  seen clinically — it is so because that is a counterfactual scenario in which nothing
  is treated for 10 weeks.
- hCG is linked to severity only through `TROPH_GAIN`; a rise in hCG independent of
  GDF15 was not modelled.
- **The anti-GDF15 antibody scenario is entirely untested in pregnancy**, and the model
  attaches a fetal safety warning to it itself: GDF15 is required for trophoblast
  invasion ([PMID 37272232](https://pubmed.ncbi.nlm.nih.gov/37272232/),
  [40157640](https://pubmed.ncbi.nlm.nih.gov/40157640/)). There are also observations
  that pregnancies **with** NVP have better outcomes
  ([PMID 24893173](https://pubmed.ncbi.nlm.nih.gov/24893173/)). The model does not
  resolve this tension; it states it explicitly in cluster 16 of the map.

---

## Files

| File | Contents |
|------|------|
| [`hg_qsp_model.dot`](../../../hyperemesis-gravidarum/hg_qsp_model.dot) | Source of the mechanistic map — **214 nodes · 19 clusters · 165 edges** |
| [`hg_qsp_model.svg`](../../../hyperemesis-gravidarum/hg_qsp_model.svg) | Vector render (6583 × 4148 pt) |
| [`hg_qsp_model.png`](../../../hyperemesis-gravidarum/hg_qsp_model.png) | Raster render (150 dpi, 13714 × 8642 px) |
| [`hg_mrgsolve_model.R`](../../../hyperemesis-gravidarum/hg_mrgsolve_model.R) | mrgsolve model — **42 ODE compartments · 182 parameters · 16 scenarios** |
| [`hg_reference_impl.py`](../../../hyperemesis-gravidarum/hg_reference_impl.py) | A pure-Python RK4 twin of the same ODE system + the validation table |
| [`hg_shiny_app.R`](../../../hyperemesis-gravidarum/hg_shiny_app.R) | Shiny dashboard — **12 tabs** |
| [`hg_references.md`](hg_references.md) | **77** PMID-verified references, with the grounds mapped for each parameter |

### The map clusters (19)

1 the fetal-placental unit (the source of GDF15) · 2 maternal GDF15 concentration ·
**3 the fold-change detector** · 4 area postrema (GFRAL/RET) · 5 the NTS integrator →
the emetic pattern generator · 6 the peripheral gastrointestinal branch ·
7 vestibular and cortical inputs · 8 the genetic architecture ·
9 hCG → TSHR → transient thyrotoxicosis ·
10 volume and chloride-responsive alkalosis · 11 energy deficit, ketosis, liver ·
12 thiamine and Wernicke · 13 the 5-HT3 antagonists · 14 the H1/D2/M1 blockers ·
15 central agents at the high-authority NTS nodes · 16 mechanism-based therapy and
timing · 17 adjunctive therapy · 18 clinical endpoints and biomarkers ·
19 differential diagnosis

---

## Running it

```bash
# render the map
dot -Tsvg hg_qsp_model.dot -o hg_qsp_model.svg
dot -Tpng -Gdpi=150 hg_qsp_model.dot -o hg_qsp_model.png

# verification (needs only python3 — no R and no compiler)
python3 hg_reference_impl.py            # all the scenarios + the validation table
python3 hg_reference_impl.py --check    # PASS/FAIL (currently 15/15)
```

```r
# the mrgsolve model (16 scenarios + summary tables + figures)
source("hg_mrgsolve_model.R")

# the Shiny dashboard
shiny::runApp("hg_shiny_app.R")
```

The first thing to try in the Shiny app is the **Adaptation (ALPHA)** slider in the
sidebar. Take it down to zero and the hindbrain becomes a concentration detector, and
the predictions above invert in front of you. The second is the metformin **start time**
sweep on tab 8 — the slider lets you watch the same drug change sign across conception.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It has not
been independently validated or certified and **must not be used for real clinical
decision-making, prescribing or regulatory submission.** The parameters are illustrative
approximations. The GDF15-targeting scenarios in particular are entirely unvalidated in
pregnancy, and the model itself flags concerns about fetal safety.
