# Tardive Dyskinesia (TD) — QSP Model

> A quantitative systems pharmacology model that treats the involuntary
> movements appearing after long-term antipsychotic use **not as a list of
> adverse effects but as a hysteresis with two memories**.
> One is **postsynaptic D2 supersensitivity (RUP)**, which reverses completely
> with a time constant of months; the other is **structural damage (SDAM,
> bistable)**, which sustains itself; and the two levers available in the clinic —
> reducing the causative drug and inhibiting VMAT2 — sit on **opposite sides of
> the same synapse**. Because of that separation, this model computes "when can
> discontinuation still reverse it", "why does discontinuation make the symptoms
> worse at first" and "why is the drug effect certain while the disease is
> unchanged" all from a single structure.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT, 20 clusters · 200 nodes · 272 edges) | [`td_qsp_model.dot`](td_qsp_model.dot) |
| 🖼️ Map (SVG) | [`td_qsp_model.svg`](td_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi) | [`td_qsp_model.png`](td_qsp_model.png) |
| ⚙️ mrgsolve ODE model (40 ODEs · 10 scenarios · 11 analysis functions) | [`td_mrgsolve_model.R`](td_mrgsolve_model.R) |
| 📊 Shiny dashboard (9 tabs) | [`td_shiny_app.R`](td_shiny_app.R) |
| 🔬 Independent implementation for verification (numpy/scipy) | [`td_reference_check.py`](td_reference_check.py) |
| 📚 References (72 items, PubMed links) | [`td_references.md`](td_references.md) |

---

## 1. The disease in one paragraph

Tardive dyskinesia is an involuntary hyperkinetic syndrome that appears after
several months or more of exposure to dopamine D2 receptor blockers
(antipsychotics, and also non-psychiatric drugs such as metoclopramide).
Oro-lingual-masticatory movements are the commonest, and choreiform movements of
the limbs and trunk, tardive dystonia, tardive akathisia and respiratory
dyskinesia may be superimposed. The incidence is 5-6% per year with
first-generation antipsychotics and 3-4% per year with second-generation agents,
rising to as much as 25-30% per year in those aged 55 and over. The mechanisms
discussed include (i) postsynaptic supersensitivity to chronic D2 blockade
(increased receptor density, an increased D2High fraction, increased
Gi/β-arrestin gain), (ii) oxidative damage from autoxidation of cytosolic
dopamine, MAO-derived H2O2 and inhibition of mitochondrial complex I, (iii) the
structural change of loss of striatal cholinergic and parvalbumin interneurons
together with a deficit of GABAergic transmission, and (iv) the resulting
imbalance between the direct and indirect pathways and pallido-thalamic
disinhibition. Treatment takes two routes: reducing the causative drug or
switching to an agent with low occupancy such as clozapine, and reducing the
presynaptic dopamine supply with a VMAT2 inhibitor (valbenazine,
deutetrabenazine). The clinical difficulty is that these two demand different
prices, that dose reduction **worsens** the symptoms over the first several
weeks, and that in some patients discontinuation does not bring improvement at
all.
**The aim of this model is to explain those three difficulties in terms of
differences in the state of the system rather than in the size of a drug
effect.**

## 2. The model's structural claim

```
(A) Reversible memory: postsynaptic D2 supersensitivity RUP
    dRUP/dt = KIN_R · Hill(OCC; 0.70, n=4) · RISKMOD · (1+0.5·E_ACH)
              − KOUT_R · RUP / (1 + 0.8·SDAM)
    → driven only above the occupancy threshold (~70%). Once the drive is gone
      it disappears completely with τ = 300 d. With this alone, TD would always
      have to be reversible.

(B) Quasi-irreversible memory: structural damage SDAM  ← the core of this model
    dSDAM/dt = KIN_S · Hill(0.55·ROS + 0.45·RUP; 0.85, n=8) · (1−SDAM)
               − KOUT_S · SDAM
    and, decisively, SDAM feeds back onto its own driver:
        ros_drive += 2.2 · Hill(SDAM; 0.50, n=4)     (mitochondrial self-amplification)
    Because this feedback is saturating (negligible below SDAM<0.3, maximal
    above >0.7), the SDAM equation is bistable:
        lower stable point SDAM = 0        (repair wins)
        unstable threshold  SDAM ≈ 0.50-0.55
        upper stable point  SDAM ≈ 0.85-0.90 (self-sustaining even off drug)

(C) The observed dyskinesia is not read directly off these two states. It is the
    output of the basal ganglia loop.
    D2STIM = DA_SYN · (1 + 2.2·RUP) · (1 − 0.55·OCC) · (1+0.5·E_ACH)/(1+0.35·E_AMA)
    EXC    = max(D2STIM − 1, 0) + 1.2·SDAM
    IND → GPe → STN → GPi/SNr → THAL → AIMS = 26 · Hill(THAL−1; 0.55, 2)
```

The two terms in (C) generate the entire clinical content of this model.

* **`(1 − 0.55·OCC)`** — ongoing blockade **masks** the symptom it has itself
  created. Lift the blockade and AIMS rises **before** it falls.
* **`+1.2·SDAM`** — a term **independent** of occupancy. Once the latch has
  closed, there is a floor on AIMS below which no manipulation of the causative
  drug can bring it.

And the two levers sit on opposite sides of the same synapse. **Lowering the D2
blockade** removes the drive to RUP but unmasks the symptoms and pays in
psychosis. **Lowering the dopamine supply** (VMAT2 inhibition) removes the
symptoms but touches neither RUP nor SDAM, and pays in parkinsonism, depression
and sedation.

## 3. Computed results

Every figure below is **the output of these equations** and not a literature
value. All of it is reproduced by
[`td_reference_check.py`](td_reference_check.py), an independent implementation
of the same equations in numpy/scipy. The index patient = 58 years old,
risperidone-equivalent 8 mg/day, FGA-level oxidative burden (RISK_FGA 1.6), D2
occupancy 0.86, AIMS (day 730) = 10.59.

### 3.1 The bistability really does exist (`TD_latch_bistability`)

Scanning the sign of dSDAM/dt in the complete absence of drug gives three roots.

| | SDAM | Meaning |
|---|---|---|
| lower stable point | 0 | repair wins — complete recovery |
| unstable threshold | 0.50 – 0.55 | go above this and there is no coming back |
| upper stable point | 0.85 – 0.90 | self-sustaining without the causative drug = persistent TD |

### 3.2 Cumulative dose is not an exposure metric (`TD_exposure_threshold`)

The same drug is given split differently across dose × duration, and followed for
up to 6 years after discontinuation.

| Dose (risp-eq mg/day) | Duration (days) | D2 occupancy | Occupancy-days | Drive-days | SDAM at discontinuation | AIMS 6 years after discontinuation | Outcome |
|---|---|---|---|---|---|---|---|
| 16 | 300 | 0.931 | 280 | 227 | 0.399 | **4.87** | persistent |
| 8 | 460 | 0.856 | 394 | 318 | 0.552 | 4.91 | persistent |
| 6 | 560 | 0.811 | 453 | 359 | 0.615 | 4.92 | persistent |
| 4 | 700 | 0.729 | 509 | 376 | 0.653 | 4.91 | persistent |
| 3 | 900 | 0.660 | 591 | 393 | 0.673 | 4.90 | persistent |
| 2 | 1500 | 0.543 | 813 | 397 | 0.212 | **0.18** | recovery |
| 1.5 | 1825 | 0.451 | 824 | 269 | 0.018 | **0.00** | recovery |

280 occupancy-days (occupancy 0.93) leaves persistent TD, while **2.9 times as
many**, 813 occupancy-days (occupancy 0.54), recovers without a trace. What
determines the outcome is not the total amount of exposure but **the time spent
above the plasticity threshold**.

### 3.3 The point of no return is a step 5 days wide (`TD_reversibility_window`)

| Time of discontinuation (days) | SDAM at discontinuation | AIMS at discontinuation | Peak AIMS after discontinuation | Day peak reached | AIMS at +2 years | AIMS at +6 years | Outcome |
|---|---|---|---|---|---|---|---|
| 90 | 0.042 | 0.37 | 2.08 | +33 | 0.12 | 0.01 | recovery |
| 180 | 0.177 | 2.78 | 5.11 | +32 | 0.79 | 0.16 | recovery |
| 270 | 0.318 | 5.18 | 7.42 | +33 | 2.14 | 0.54 | recovery |
| **285** | 0.341 | 5.53 | 7.73 | +33 | 2.95 | **0.90** | recovery |
| **290** | 0.348 | 5.64 | 7.83 | +33 | 4.43 | **4.82** | persistent |
| 365 | 0.449 | 7.10 | 9.08 | +33 | 6.09 | 4.88 | persistent |
| 730 | 0.733 | 10.59 | 12.02 | +30 | 7.39 | 4.96 | persistent |
| 1460 | 0.872 | 12.36 | 13.53 | +29 | 8.34 | 5.04 | persistent |

t\* = day 285–290 (about 9.5 months). What is notable is that **SDAM at that
moment is 0.34, well below the self-sustaining threshold of 0.50–0.55**. Because
RUP goes on driving the damage for months after the drug is stopped, the state at
the time of discontinuation alone does not tell you whether that patient has
already passed the window. **What the model says: whether it is reversible is
determined not by the present state but by the trajectory to come.**

### 3.4 The withdrawal paradox and its crossover time (`TD_withdrawal_crossover`)

The dose is changed at day 730 and compared against the trajectory in which it is
maintained.

| Strategy | d730 | d744 | d760 | d820 | d1095 | d1825 | d2555 | Peak | Crossover time |
|---|---|---|---|---|---|---|---|---|---|
| maintain 8 mg | 10.59 | 10.66 | 10.74 | 11.01 | 11.82 | 12.61 | 12.80 | 12.80 | — |
| reduce to 6 mg | 10.59 | 10.79 | 10.92 | 11.13 | 11.81 | 12.49 | 12.66 | 12.66 | +318 days |
| reduce to 4 mg | 10.59 | 10.97 | 11.20 | 11.29 | 11.66 | 12.08 | 12.20 | 12.20 | +237 days |
| reduce to 2 mg | 10.59 | 11.26 | 11.64 | 11.43 | 10.79 | 10.09 | 9.88 | 11.64 | +154 days |
| complete discontinuation | 10.59 | 11.60 | 12.02 | 11.55 | 9.44 | 6.16 | 5.14 | 12.02 | **+135 days** |
| discontinuation + clozapine | 10.59 | 11.48 | 11.95 | 11.52 | 9.55 | 6.60 | 5.70 | 11.96 | +134 days |

The larger the reduction, the larger the initial worsening (+1.43 AIMS, peaking
on day 30), the sooner the crossover, and the larger the eventual benefit.
**Reduce slowly and the worsening is small, but the crossover is late and the
benefit small too** — these are two ends of the same curve, not two different
strategies. And the fact that every reduction looks like a loss for 4-5 months
explains why a clinician who has started a reduction ends up **reversing that
reduction in exactly that interval**.

### 3.5 VMAT2 inhibitors: dose-response and CYP2D6 (`TD_vmat2_dose_response`, `TD_cyp2d6_panel`)

At the 6-week time point (day 772), with the antipsychotic maintained unchanged.

| Drug / dose | Cave (ng/mL) | VMAT2 occupancy | Synaptic DA | AIMS | ΔAIMS | % | PARK | DEPR | QTc (ms) |
|---|---|---|---|---|---|---|---|---|---|
| valbenazine 40 mg | 12.0 | 0.316 | 1.124 | 8.70 | −1.89 | −17.9 | 0.449 | 0.190 | 0.96 |
| valbenazine 80 mg | 24.0 | 0.480 | 0.892 | 7.08 | **−3.51** | −33.2 | 0.579 | 0.261 | 1.92 |
| valbenazine 120 mg | 36.0 | 0.581 | 0.740 | 5.82 | −4.78 | −45.1 | 0.677 | 0.304 | 2.88 |
| deutetrabenazine 24 mg | 10.0 | 0.345 | 1.085 | 8.44 | −2.15 | −20.3 | 0.468 | 0.203 | 2.50 |
| deutetrabenazine 36 mg | 15.0 | 0.441 | 0.949 | 7.50 | **−3.09** | −29.2 | 0.544 | 0.245 | 3.75 |
| deutetrabenazine 48 mg | 20.0 | 0.513 | 0.843 | 6.68 | −3.91 | −36.9 | 0.612 | 0.277 | 5.00 |

The −3.51 (−33%) of valbenazine 80 mg and the −3.09 (−29%) of deutetrabenazine
36 mg are the same order of magnitude as the effects reported in KINECT-3 and in
ARM-TD/AIM-TD (literature: about −3.2 and −3.0 respectively). **What matters is
the last columns of that same table**: at every dose, SDAM at the y5 time point
is identically 0.885, and RUP is essentially unchanged.

CYP2D6 moves efficacy and toxicity **together along a single exposure axis** — at
a fixed valbenazine dose of 80 mg:

| Phenotype (CL multiplier) | Cave | VMAT2 occupancy | ΔAIMS % | PARK | QTc (ms) |
|---|---|---|---|---|---|
| UM (1.6×) | 15.0 | 0.366 | −22.1 | 0.482 | 1.20 |
| normal (1.0×) | 24.0 | 0.480 | −33.2 | 0.579 | 1.92 |
| IM (0.7×) | 34.3 | 0.569 | −43.5 | 0.664 | 2.74 |
| PM (0.5×) | 48.0 | 0.649 | −54.2 | 0.740 | 3.84 |
| PM + CYP3A4 inhibition (0.35×) | 68.5 | 0.725 | −62.6 | 0.807 | 5.48 |

The model's reading, in other words, is that the label's dose cap for PMs exists
not because "there is no effect" but **because benefit and harm grow together on
the same axis**.

### 3.6 Suppression and modification are different verbs (`TD_suppression_vs_modification`)

Treatment for 2 years (y2→y4), then stopped at y4.

| Intervention | AIMS y2 | AIMS y4 | RUP y4 | SDAM y4 | AIMS y4+8 weeks | AIMS y6 |
|---|---|---|---|---|---|---|
| no intervention | 10.59 | 12.36 | 1.445 | 0.872 | 12.41 | 12.74 |
| valbenazine 80 (stopped at y4) | 10.59 | **8.88** | 1.413 | 0.872 | **12.33** | 12.71 |
| deutetrabenazine 36 (stopped at y4) | 10.59 | 9.30 | 1.418 | 0.872 | 12.35 | 12.72 |
| switch to clozapine (y2) | 10.59 | **7.70** | **0.295** | 0.871 | **7.48** | **6.00** |
| clozapine + valbenazine | 10.59 | 4.69 | 0.284 | 0.872 | 7.39 | 5.97 |

The VMAT2 inhibitor lowers AIMS by 3.5 points over two years, but 8 weeks after
stopping it is indistinguishable from the no-intervention arm (12.33 vs 12.41).
There is **no** rebound overshoot — only the effect disappears. Switching to
clozapine, by contrast, brings RUP down from 1.45 to 0.30, and the benefit
persists even after discontinuation (y6: 6.00 vs 12.74). In this model the only
disease-modifying lever is the switch to clozapine.

### 3.7 Combination: a bridge early, sub-additive late (`TD_combination_interaction`)

A 2×2 of switch-to-clozapine × valbenazine. The effect is the change in AIMS
relative to day 730.

| Time point | neither | valbenazine | clozapine | both | E(valbenazine) | E(clozapine) | Interaction |
|---|---|---|---|---|---|---|---|
| day 772 (6 weeks) | +0.21 | −3.51 | **+1.32** | **−2.35** | −3.72 | +1.11 | +0.05 |
| day 1095 | +1.23 | −2.33 | −1.04 | −4.60 | −3.56 | −2.27 | −0.00 |
| day 1825 | +2.02 | −1.42 | −3.99 | −5.82 | −3.44 | −6.02 | +1.62 |
| day 2555 | +2.21 | −1.21 | −4.90 | −5.78 | −3.42 | −7.11 | **+2.54** |

At the 6-week time point the switch to clozapine on its own **worsens** AIMS
(+1.32) while the combination improves it (−2.35). That is, one of the clinical
roles of a VMAT2 inhibitor is computed out here: **a bridge that covers the
withdrawal-emergent worsening during the switch**. Conversely, at the 7-year time
point the interaction term is +2.54, i.e. sub-additive — both are pushing against
the same structural floor and so eat into each other's room. A single combination
is both complementary and redundant, depending on when you look.

### 3.8 The frontier of the opposed levers (`TD_opposed_levers`)

At the day 1095 time point. The same reduction in AIMS is paid for with
different things.

| Strategy from y2 | D2 occupancy | VMAT2 occupancy | AIMS | PSYCH | PARK | DEPR | ADHER | RUP |
|---|---|---|---|---|---|---|---|---|
| maintain 8 mg | 0.861 | 0 | 11.82 | 0.262 | 0.325 | 0.115 | 0.776 | 1.315 |
| reduce to 4 mg | 0.736 | 0 | 11.66 | 0.326 | 0.187 | 0.083 | 0.816 | 1.146 |
| reduce to 2 mg | 0.555 | 0 | 10.79 | 0.496 | 0 | 0.038 | 0.859 | 0.839 |
| complete discontinuation | 0 | 0 | 9.44 | **1.000** | 0 | 0 | 0.769 | 0.520 |
| switch to clozapine 350 | 0.279 | 0 | 9.55 | **0.391** | 0 | 0.003 | **0.894** | **0.548** |
| valbenazine 40 | 0.853 | 0.316 | 9.84 | 0.265 | 0.392 | 0.203 | 0.736 | 1.307 |
| valbenazine 80 | 0.842 | 0.480 | 8.26 | 0.270 | **0.524** | 0.284 | 0.686 | 1.295 |
| valbenazine 120 | 0.834 | 0.581 | 7.00 | 0.274 | 0.623 | 0.333 | 0.655 | 1.286 |
| 4 mg + valbenazine 80 | 0.704 | 0.480 | 7.99 | 0.349 | 0.368 | 0.249 | 0.720 | 1.101 |

On AIMS alone valbenazine 120 is the best; once psychosis is included, complete
discontinuation is the worst; and **the only thing that improves all three axes
(AIMS · psychosis · motor adverse effects) at once is the switch to clozapine.**
One more thing computed out here: because of the loop by which adverse effects
erode adherence (ADHER 0.776 → 0.686), a high dose of a VMAT2 inhibitor also
makes exposure to the causative drug irregular.

### 3.9 Host and adjunctive treatment (`TD_risk_scan`)

The same prescription (8 mg risp-eq, 5 years), different patients.

| Host | AIMS y1 | AIMS y3 | AIMS y5 | RUP y5 | SDAM y5 | Day of latch crossing |
|---|---|---|---|---|---|---|
| 25 years old, SGA | 4.23 | 10.37 | 11.43 | 1.175 | 0.880 | 585 |
| 58 years old, SGA | 6.58 | 11.74 | 12.60 | 1.512 | 0.884 | 458 |
| 58 years old, FGA (reference) | 7.10 | 11.82 | 12.61 | 1.518 | 0.885 | **400** |
| 75 years old, FGA | 8.35 | 12.73 | 13.46 | 1.833 | 0.886 | 366 |
| 58 years old, FGA, diabetes | 8.63 | 12.92 | 13.64 | 1.911 | 0.886 | 357 |
| 58 years old, FGA, oestrogen-replete | 5.47 | 10.43 | 11.28 | 1.131 | 0.885 | 424 |
| 58 years old, FGA, high-risk genotype | 8.68 | 13.08 | 13.81 | 1.987 | 0.886 | 382 |
| 58 years old, FGA + benztropine 2 mg | **10.69** | 14.40 | 15.01 | 1.964 | 0.886 | 382 |
| 58 years old, FGA + Ginkgo biloba EGb761 | **4.64** | 11.03 | 12.43 | 1.468 | 0.870 | **780** |

Co-administration of an anticholinergic raises the 1-year AIMS from 7.10 to 10.69
(it is not unmasking, it raises the D2 stimulation gain itself). Ginkgo biloba
(antioxidant) is computed to be **the only adjunctive therapy that delays the
latch itself** (latch at day 400 → 780). At the 5-year point, however, the
difference narrows to 12.61 versus 12.43 — **the model's conclusion is "delay",
not "prevention".**

### 3.10 Results that run against the hypothesis the model set out to express

This model was originally written in the expectation that "suppressing the
symptoms with a VMAT2 inhibitor creates room to reduce the causative drug in the
meantime, so the combination will be synergistic". The computed result **depends
on when you look**: at 6 weeks it is complementary (a bridge), but at 7 years the
interaction term is +2.54, i.e. **sub-additive**. Another thing that ran contrary
to expectation is that **the size of the withdrawal worsening does not determine
which strategy is better** — complete discontinuation has both the highest peak
AIMS (12.02) and the best final outcome (5.14), whereas a small dose reduction
looks the safest while its 5-year benefit is nearly zero (12.49 vs 12.61). And a
third: antioxidant adjunctive therapy delays the latch by almost a factor of two
and yet improves the 5-year AIMS by only 1.4% — in this model, **an intervention
that buys time and an intervention that changes the outcome are different
things.**

## 4. Parameter calibration anchors

| Target | Literature value | Model parameter / output |
|---|---|---|
| risperidone active moiety PK | CL/F ~5 L/h, t½ ~20 h, 4 mg → ~33 ng/mL | `CL_AP` 120 L/day, `V2_AP` 100 L |
| striatal D2 occupancy | EC50 ~12 ng/mL (4 mg 73%, 8 mg 86%) | `EC50_D2` 12, `HILL_D2` 1.25 |
| therapeutic / EPS threshold | 65% / 78-80% | `OCC50_P` 0.58, `OCC50_PK` 0.78 |
| clozapine | 350 mg → ~390 ng/mL, D2 20-40% | `CL_CLZ` 900 L/day → OCC 0.28 |
| valbenazine exposure | NBI-98782 Cave ~24 ng/mL (~2× in PMs) | `FM_VAL` 0.30, `CL_NBI` 1000, `CYP2D6` |
| deutetrabenazine exposure | (α+β)-HTBZ ~15 ng/mL, QTc ~4 ms | `FM_DTB` 0.50, `QT_HTB` 0.25 |
| KINECT-3 (80 mg, 6 weeks) | AIMS −3.2 (baseline ~10) | model −3.51 (baseline 10.59) |
| ARM-TD / AIM-TD (36 mg) | AIMS ~−3.0 | model −3.09 |
| annual incidence | FGA 5-6%/year, SGA 3-4%/year, >55 years 25-30%/year | `RISK_FGA`, RISKMOD = 1+0.015(age−40) |
| withdrawal-emergent dyskinesia | peaks within 4-6 weeks of dose reduction | model peak = day 30 after discontinuation |
| persistence after discontinuation | about one third recover, the rest persist for years | SDAM bistability (threshold 0.52), `KOUT_S` 3e-4 |
| Ginkgo biloba EGb761 240 mg | AIMS improvement (RCT) | `ANTIOX_MAX` 0.85 → latch day 400 → 780 |

## 5. How to run

```bash
# 1) render the mechanistic map
dot -Tsvg td_qsp_model.dot -o td_qsp_model.svg
dot -Tpng -Gdpi=150 td_qsp_model.dot -o td_qsp_model.png

# 2) mrgsolve model — runs all 11 analyses and prints the tables
Rscript td_mrgsolve_model.R

# 3) Shiny dashboard (9 tabs)
R -e 'shiny::runApp("td_shiny_app.R", port = 8080)'

# 4) reproduce the same figures without R / mrgsolve (only numpy + scipy needed)
python3 td_reference_check.py          # everything
python3 td_reference_check.py --quick  # scenarios only
```

`td_reference_check.py` is an independent implementation that transcribes the
`$ODE` block of `td_mrgsolve_model.R` as it stands, and it regenerates every
table in section 3 above. The two implementations give identical results because
(i) states and parameters correspond 1:1, (ii) both treat dosing as a **zero-order
daily input**, and (iii) both treat a change of prescription as a smooth step of
the form `0.5(1+tanh((t−t0)/1))`, so the solver never meets a discontinuity.

## 6. Model composition

**PK / exposure (19 states)**: antipsychotic oral, LAI depot, central,
peripheral and effect site; clozapine 2-compartment; valbenazine → NBI-98782
(including peripheral); deutetrabenazine → (α+β)-HTBZ (including peripheral);
amantadine; benztropine; Ginkgo biloba effect compartment; botulinum local effect
compartment; cumulative occupancy-time.

**Disease / circuit / endpoints (21 states)**: cytosolic, vesicular and synaptic
dopamine; depolarisation block; D2 supersensitivity RUP; oxidative stress ROS;
structural damage SDAM; indirect-pathway MSN, GPe, STN, GPi and thalamus; AIMS;
psychosis burden; parkinsonism; mood burden; QTc; adherence; and the multiplier
for the clinician's dose-escalation policy.

**10 scenarios**: natural history · FGA-level high occupancy · 50% dose reduction
· complete discontinuation · switch to clozapine · valbenazine ·
deutetrabenazine · benztropine co-administration · Ginkgo biloba + amantadine ·
clinician dose-escalation policy.

**11 analysis functions**: one for each table in section 3 above.

## 7. Limitations

* AIMS, psychosis, parkinsonism and mood burden are **normalised burden
  indices**, and PARK/DEPR/PSYCH are not calibrated in the units of an actual
  scale such as SAS, PANSS or an incidence rate. Only their **relative changes
  and ordering** can be interpreted.
* SDAM is an aggregate variable for "structural damage" and is not calibrated to
  the rate of loss of any particular cell population (cholinergic interneurons,
  parvalbumin FSIs). The bistability is a **structural hypothesis** and not a
  fact verified in human data.
* The antipsychotics have been lumped into a single risperidone equivalent. Fast
  dissociation (clozapine, quetiapine) is approximated only through occupancy and
  `EC50_D2_CLZ`, and D2 partial agonists (aripiprazole) are on the map but not in
  the ODEs.
* Dosing is a zero-order daily input. Abrupt swings in occupancy caused by real
  missed doses (a common trigger for withdrawal-emergent dyskinesia) are
  represented only as an average decrement in the `ADHER` state.
* This model is a **single-individual simulator** and has no stochastic structure
  for population variability (IIV) or for predicting incidence. The "incidence of
  5%/year" was used only for parameter calibration.
* The effect sizes of Ginkgo biloba, vitamin E and amantadine lie in an area of
  weak evidence, and the model's "delay" conclusion is no more than a structural
  prediction built on that weak evidence.

## 8. ⚠️ Disclaimer

This model is a **qualitative to semi-quantitative QSP model for educational and
research purposes**. It was constructed from the published literature and
clinical trial data but has not been independently validated or certified, and
**must not be used directly for actual clinical decision-making, prescribing, or
regulatory submission.** In particular, reducing or discontinuing an
antipsychotic carries the serious risk of relapse of psychosis, and no output of
this model can be used as a basis for that judgement.
