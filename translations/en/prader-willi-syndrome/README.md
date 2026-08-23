# Prader-Willi Syndrome (PWS) — QSP Model

**Prader-Willi Syndrome · one convertase, five branches, a bistable appetite**

> PWS is not written here as a **list of hormone deficiencies**. The lesion is written
> as **a single scalar** (the activity of the prohormone convertase PC1/3), and each of
> the five branches diverging from it is given exactly one **dimensionless escape
> ratio**. Then which branch fails first, which branch accumulates precursor, and
> **which measurement is in principle blind to which information** all fall out of
> algebra rather than of fitting. At that point a good part of the PWS literature
> becomes arithmetic rather than assertion.

| | |
|---|---|
| **Mechanistic map** | [`pws_qsp_model.svg`](pws_qsp_model.svg) · [`pws_qsp_model.dot`](pws_qsp_model.dot) — 223 nodes / 316 edges / 20 clusters |
| **mrgsolve model** | [`pws_mrgsolve_model.R`](pws_mrgsolve_model.R) — **64 ODEs**, 15 scenarios, 10 drugs |
| **Shiny app** | [`pws_shiny_app.R`](../../../prader-willi-syndrome/pws_shiny_app.R) — 10 tabs (mechanism drawn alongside the endpoints) |
| **References** | [`pws_references.md`](pws_references.md) — **357 papers**, every PMID looked up in PubMed |
| **Reference implementation** | [`pws_reference_model.py`](../../../prader-willi-syndrome/pws_reference_model.py) · [output](../../../prader-willi-syndrome/pws_reference_output.txt) — dependency-free Python RK4 |
| **Calibration and derived-quantity analysis** | [`pws_calibration.py`](../../../prader-willi-syndrome/pws_calibration.py) · [output](../../../prader-willi-syndrome/pws_calibration_output.txt) — all in closed form |

<a href="pws_qsp_model.svg"><img src="pws_qsp_model.png" width="820" alt="PWS QSP map"></a>

---

## 1. What is written differently here

### ① One enzyme, five substrates — and a single escape ratio sets the order

Loss of paternal expression at 15q11-q13 converges through SNORD116, NDN, and MAGEL2
onto NHLH2, and the outcome is a single thing: **reduced PC1/3 (PCSK1) activity** in
the hypothalamus and the islets. The model writes that as the scalar `PC13` and hangs
five precursors on it. Each branch has exactly one **escape ratio**:

```
eps_i = d_i · Km_i / kcat_i          (precursor escape ÷ processing)
```

Solving the steady state in the sub-saturating regime gives

```
L_i = 1 − (1+eps_i)/(1+eps_i/PC13)     product loss           [Eq. A]
R_i = (1+eps_i)/(PC13+eps_i)           precursor accumulation  [Eq. B]
```

At `PC13 = 0.40`, with no per-branch fitted values:

| Branch | eps | Product loss | Precursor |
|---|---|---|---|
| Pro-oxytocin → oxytocin | 6.00 | **−56.2%** | 1.09× |
| POMC → α-MSH | 1.50 | **−47.4%** | 1.32× |
| proGHRH → GHRH | 1.00 | −42.9% | 1.43× |
| Proinsulin → insulin | 0.15 | −16.4% | **2.09×** |
| Proghrelin → acyl-ghrelin | 0.10 | −12.0% | **2.20×** |

The two columns **run in opposite directions.** `L` grows with eps and `R` shrinks
with eps, so **the branch that loses product is precisely the branch that does not
accumulate precursor**. This is not a choice but a compulsion.

### ② So the two easy measurements are each exactly blind to the other

Divide Eq. B by (1 − Eq. A) and eps **cancels**:

```
precursor / product = (1+eps/PC13)/(PC13+eps) = 1 / PC13
```

All five branches give exactly 2.5000000000 (PC13 = 0.4). And since steady-state mass
conservation gives `processing + escape = S`:

| What you measure | What it actually measures | What it is blind to |
|---|---|---|
| Precursor + product (cross-reacting antibody) | the synthesis rate `S` | **exactly to PC1/3** |
| Precursor / product (the ratio) | `1/PC13` | **exactly to eps** — the same value in every branch |
| Product alone | `S/(1+eps/PC13)` | the only quantity that carries branch information |

That is, **the mutual disagreement of the plasma and CSF oxytocin literature in PWS is
not noise but a prediction**. An antibody that recognises precursor and product
together measures the synthesis rate, and that value is independent of the lesion.
And a trial that uses such a measurement as a target-engagement biomarker is looking
at **a quantity that cannot respond**. Meanwhile the proinsulin:insulin ratio and the
pro-oxytocin:oxytocin ratio must give **the same number** in the same patient, and
that number grades the convertase rather than the phenotype.

### ③ Satiety is a harmonic mean, and the melanocortin arm is in series, not in parallel

```
SAT = 1 / Σ (w_i / x_i)
```

In a harmonic mean **the weakest arm is rate-limiting**. And α-MSH is not an arm of
its own — because the MC4R satiety signal passes by way of the PVN oxytocin neurones,
α-MSH enters the oxytocin arm as **a saturating input gain**:

```
relay(MC) = (MC/(K+MC))·(K+1),   K = 0.20  →  relay(∞) = 1.20  (ceiling)
```

The block in PWS is in pro-oxytocin processing, that is, **below MC4R**. So MC4R
agonism is pressed down from above by exactly the quantity PWS lacks (oxytocin
availability). At the values for a 12-year-old with PWS:

| Intervention | SAT | ΔSAT |
|---|---|---|
| Baseline | 0.4562 | — |
| **All four non-oxytocin arms normalised** | 0.4830 | +0.0268 |
| MC4R arm driven to saturation | 0.5379 | +0.0817 |
| OXTR arm ×2.1 (carbetocin level) | 0.7139 | +0.2577 |

Fixing **all** four arms only reaches 0.483. A harmonic mean is not rescued by healthy
terms. And setmelanotide already reaches most of the relay ceiling at 3 mg, so
**there is almost nothing left for a stronger MC4R agonist** — this is the structural
statement of why the same molecule works in POMC, LEPR, and PCSK1 deficiency (lesion
above MC4R) and does not work in PWS (below).

The oxytocin receptor sits on **both sides of the synapse**. Postsynaptically it
raises the satiety arm; presynaptically it inhibits AgRP. The whole of the difference
between carbetocin and setmelanotide comes out of this asymmetry, and there is no
per-drug fitted parameter.

### ④ Food seeking is a bistable state — food security does not lower the drive, it restores a state

```
dS/dt = k_on·DRVe·(1−S) + G(REINF)·S⁴/(K⁴+S⁴)·(1−S) − k_off·S·SAT
```

Sweeping CUE (external food cues and access) changes the fixed-point structure:

| CUE | Fixed points |
|---|---|
| 0.15 | 0.0361 stable · 0.2630 **saddle** · 0.7305 stable |
| 0.28 (standard management) | 0.0663 stable · 0.2489 **saddle** · 0.7332 stable |
| 0.45 | 0.1090 stable · 0.2239 **saddle** · 0.7366 stable |
| 0.60 (**above 0.58**) | **a single stable point at 0.7395** |

At `CUE* = 0.5785` a **saddle-node** occurs and the low state is **annihilated.**
Above it, no dose can hold on to a state the vector field has erased. The control
group is monostable (0.0504) at any CUE.

So food security is **not on the same axis** as the other interventions. It does not
change the depth of a state, it changes **which states exist**. And this is why both
DESTINY-PWS and CARE-PWS required a stable food-secure environment as an enrolment
condition, and why their effect sizes must not be extrapolated to a family without
such an environment.

### ⑤ Growth hormone arrives at the airway twice — the window is not a fit

The fast clock: IGF-1 → lymphoid and adenotonsillar proliferation (τ ≈ 20 d,
regressing with a 100 d adaptation).
The slow clock: lean mass → respiratory muscle strength (τ = 75 d) + loss of fat mass
(months).

Starting GH at age 6:

| Week | LYMPH | RMS | %fat | AHI (ordinary airway) | AHI (vulnerable airway) |
|---|---|---|---|---|---|
| 0 | 1.262 | 0.794 | 50.3 | 11.18 | 11.18 |
| 4 | 1.668 | 0.843 | 48.7 | **12.80** | **21.17** |
| 8 | 1.612 | 0.884 | 47.0 | 11.79 | 19.63 |
| 52 | 1.445 | 1.068 | 31.8 | 4.93 | 7.82 |
| 104 | 1.384 | 1.175 | 17.8 | 2.90 | 4.68 |

That the peak sits at 4-8 weeks is not something fitted. It is interference between
`τ_lymphoid < τ_muscle < τ_fat`, and in closed form too
`t_peak = ln((A·k_f)/(B·k_s))/(k_f − k_s)` gives 3-8 weeks — the quantity least
sensitive to the amplitudes. And airway vulnerability multiplies **only the
obstructive term**, so it enlarges the window and does not enlarge the plateau. That
is why the same drug is dangerous for eight weeks and protective thereafter. This is
the mechanistic content of the recommendations for "polysomnography before starting
GH and 6-8 weeks after starting" and "do not start during a respiratory infection".

### ⑥ Intake is tied to **what is served** — so weight gain comes before hyperphagia

The caregiver serves age-normative portions. Because expenditure in PWS is low (lean
mass ↓ × activity ↓, 0.72× the predicted value), **normal intake becomes surplus**.
Without any appetite parameter changing, the Miller stages are reproduced in order:

| Age | Stage | Weight | %fat | EI/requirement | Weight SDS | HQ-CT | SEEK |
|---|---|---|---|---|---|---|---|
| 0.5 | 1a | 3.9 | 30.8 | **0.484** | −3.79 | 0.1 | 0.005 |
| 1.0 | 1b | 6.0 | 37.6 | 0.629 | −3.11 | 0.1 | 0.006 |
| 2.0 | 1b | 12.7 | 55.7 | 0.729 | **−0.07** | 0.2 | 0.011 |
| 3.0 | **2a** | 16.6 | 57.4 | 0.767 | **+1.00** | 0.5 | 0.019 |
| 6.0 | 2b | 23.7 | 53.0 | 0.855 | +1.32 | 4.0 | 0.093 |
| 8.0 | **3** | 33.4 | 50.9 | 1.000 | +2.51 | **14.6** | **0.698** |

In stage 1a the oral-motor constraint (`SUCK`) keeps EI/requirement below 1, and
HQ-CT only rises **five years after the weight SDS has turned positive** in stage 2a.
The hyperphagia gate is itself a developmental clock.

---

## 2. The drug panel — biomarkers and endpoints are nearly orthogonal

Age 12, food security, background growth hormone. At 8-13 weeks, placebo-corrected
ΔHQ-CT:

| Arm | ΔHQ-CT | ΔAG | ΔAG:UAG | Δ%fat | ΔHbA1c |
|---|---|---|---|---|---|
| Carbetocin 3.2 mg TID | **−1.88** | +2.0% | −0.0% | +0.04 | +0.00 |
| Carbetocin 9.6 mg TID | **−1.36** | +1.5% | −0.0% | +0.24 | +0.00 |
| DCCR 5.1 mg/kg | **−1.53** | +5.4% | +0.0% | +0.77 | **+0.38** |
| Setmelanotide 3 mg | −0.30 | +0.8% | −0.0% | +0.78 | +0.01 |
| Livoletide 60 µg/kg | −0.21 | +0.2% | **−21.9%** | +1.13 | +0.02 |
| Octreotide LAR 30 mg | −0.50 | **−72.4%** | +0.0% | +1.36 | +0.96 |
| Semaglutide 2.4 mg QW | −2.77 | +4.2% | −0.0% | **−5.35** | −0.16 |
| Metformin 1500 mg | +0.16 | +3.5% | −0.0% | +0.44 | −0.14 |

**The ghrelin columns and the HQ-CT column are nearly orthogonal.** The drug that
moves the ghrelin biomarker most moves the endpoint least. In the model this is a
single number — GHS-R1a occupancy is **already saturated** at PWS concentrations
(K = 300 pg/mL, PWS ≈ 870):

| AG (pg/mL) | Occupancy | Arm/control | Elasticity |
|---|---|---|---|
| 350 (control) | 0.5385 | 1.0000 | 0.462 |
| 870 (PWS, 2.65× control) | 0.7436 | **1.3809** | 0.256 |
| 1800 | 0.8571 | 1.5918 | 0.143 |

All that a 2.65-fold hyperghrelinaemia buys is **+38.1%** on the receptor arm, and
abolishing it completely cannot return more than **3.8%** of the drive. Octreotide is
the extreme case — it corrects the biomarker best, moves the endpoint least, and once
background GH is taken away it collapses an axis PWS is already short of (ΔIGF-1
−83.6% against −61.7% for stopping GH alone).

### Carbetocin's inverted dose-response is analytical

```
E(C) = E1·C/(K1+C) − E2·C/(K2+C)          OXTR action − V1a cross-activity
E'(C) = 0  →  C* = [K1·√(E2K2) − K2·√(E1K1)] / [√(E1K1) − √(E2K2)]
             = 1.2121 ng/mL,   E(C*) = 0.4924,   E(∞) = 0.0500
```

| mg TID | Cavg (ng/mL) | E(C) | E/E(C*) | ΔHQ-CT |
|---|---|---|---|---|
| 0.8 | 0.516 | 0.4233 | 0.8597 | −1.69 |
| 1.6 | 1.031 | 0.4897 | 0.9946 | −1.93 |
| **3.2** | 2.063 | 0.4649 | 0.9443 | **−1.89** |
| 6.4 | 4.125 | 0.3693 | 0.7501 | −1.59 |
| **9.6** | 6.187 | 0.3031 | 0.6156 | **−1.36** |
| 14.4 | 9.281 | 0.2416 | 0.4906 | −1.12 |

The 90% window is C ∈ [0.594, 2.503] ng/mL, that is **0.92-3.88 mg TID**. 3.2 mg is
inside it, 9.6 mg is on the descending limb (62% of the maximum effect). The model
reproduces the dose ordering of CARE-PWS **from V1a cross-activity alone, with no
per-dose parameter**, and predicts that the useful range is **closed on both sides**
— that a dose-finding study going below 3.2 mg would find the effect **maintained**
is a testable prediction.

### DCCR's therapeutic index is fixed in one stroke by two EC50s

Efficacy (AgRP KATP) and hyperglycaemia (β-cell KATP) are **two tissues carrying the
same channel**. The ratio is monotone and bounded: 0.8381 as `C→0`, 0.5714 as `C→∞`
(both ends closed by the Emax and EC50 ratios). No dose escapes this ratio.
Only a tissue-selective opener can change it, and that is the chemical target the
model puts quantitatively.

**A subgroup effect may be a branch rather than noise:** the same 5.1 mg/kg gives
−1.53 in the moderate stratum and −1.74 in the severe stratum (baseline HQ-CT 9.96
against 10.82, SEEK 0.066 against 0.155). This is because in the severe stratum the
patient sits closer to the upper branch, so lowering the drive can carry them across
the separatrix. If that is so, the effect will **replicate**.

---

## 3. The latch, and the development target the model proposes

A patient who has spent ages 10-12 with free access sits at SEEK = 0.7065. Two years
after standard management is restored:

| | SEEK | HQ-CT | %fat | REINF |
|---|---|---|---|---|
| Month 0 | 0.7065 | 15.72 | 50.8 | 0.500 |
| Month 6 | 0.6903 | 17.32 | 43.4 | 0.642 |
| Month 24 | 0.6922 | 17.44 | **38.4** | 0.667 |

**Standard management removes the fat but does not remove the state.** And nothing in
the panel gets the patient off the upper branch (at 6 months: carbetocin 3.2 mg
0.6239, DCCR 0.6923, semaglutide 0.6694, strict environment 0.6880, strict
environment + carbetocin 0.6204 — all upper). So the model does not claim that "some
dose will do it" but instead **solves for the value required**:

```
Current satiety                        SAT   = 0.4653
Upper saddle-node of a latched patient  SAT_c = 0.7696   (1.65× the achieved value)
Oxytocin arm gain required                      2.40×
Carbetocin's analytical ceiling                 1.49×
Shortfall                                       1.61×
```

That is, this is **not a potency problem but an efficacy problem**, and carbetocin is
already at its own optimum. And every HQ-CT reduction the panel produces happens
**inside the latched state** — a trial powered on mean HQ-CT can succeed with every
patient still latched. Because SEEK is bistable, the model predicts that the HQ-CT
response will be **a bimodal distribution** rather than a shift.

(One amusing bit of directionality: **the latched patient is the easier target.**
Being latched under management means being chronically unrewarded, and that lowers the
self-sensitisation gain propping the state up — a patient on the lower branch needs
3.45× where a latched patient needs 2.40× (`pws_calibration_output.txt` §6b). The
patients most likely to be moved off a branch by an OXTR agonist are the stratum with
the highest baseline HQ-CT in a food-secure environment, and that is the same stratum
the DCCR subgroup analysis pointed at.)

---

## 4. Natural history and the metabolic paradox

At age 25:

| | Weight | BMI | %fat | Height | Height SDS | HQ-CT | AHI | TEE/requirement |
|---|---|---|---|---|---|---|---|---|
| Non-PWS control | 68.9 | 22.3 | 20.7 | 175.6 | −0.19 | 4.2 | 1.6 | 0.993 |
| PWS free access | 122.9 | 55.1 | 68.1 | 149.3 | −3.75 | 14.5 | 22.7 | 1.063 |
| PWS food security | 59.3 | 26.0 | 29.0 | 150.9 | −3.54 | 12.0 | 3.5 | 0.714 |
| PWS security + GH (from age 1) | 62.4 | 21.6 | **19.0** | **169.9** | **−0.96** | 11.4 | 3.0 | 0.963 |

Growth hormone moves height, lean mass, and %fat, and **does not move HQ-CT**
(at age 14, untreated 10.2 against GH 9.8). The body-composition axis and the satiety
axis share the lesion but do not share the pathway. In scoliosis the two effects
almost cancel — Cobb 12.5° against 12.7° at age 14, while height is +18 cm. This is
because GH worsens it by raising the growth velocity and improves it by raising muscle
tone, and that is exactly what the trials reported.

**The metabolic paradox** comes out of one parameter, not two. In a fat-mass-matched
comparison:

| | BMI | %fat | Insulin | HOMA-IR | Adiponectin | HbA1c |
|---|---|---|---|---|---|---|
| PWS, free access | 55.1 | 68.1 | **13.1** | 3.29 | **6.2** | 5.67 |
| Non-PWS, hyperphagic | 50.8 | 66.4 | 15.2 | 3.70 | 3.8 | 5.56 |

The same convertase that starves the satiety arm also throttles proinsulin
processing. So PWS arrives at the same fat mass with **less insulin and more
adiponectin**.

**Puberty** too is a two-directional effect of the same deletion. MKRN3 is a brake, so
losing it **brings the gate forward**, while the hypothalamic lesion **lowers the
amplitude** — a puberty that begins early, is weak, and does not progress
(testosterone at age 18, control 600 against PWS 132 ng/dL). The late closure of the
bone age is why even late GH buys height, and why sex hormone replacement closes that
window as it takes effect.

**The trade-off in management** is explicit too. Unrealised drive is exactly the
frustration term, so:

| Environment | %fat | HQ-CT | Behaviour score | Frustration |
|---|---|---|---|---|
| Free access | 53.6 | 15.57 | 38.0 | 0.183 |
| Normal portions, no titration | 32.8 | 17.91 | **42.2** | **0.330** |
| Standard PWS management | 34.5 | 9.80 | 35.4 | 0.074 |
| Strict environment | 32.7 | 9.46 | 35.1 | 0.063 |

Carbetocin is **the only arm in the panel that lowers both at once** — because
oxytocin enters the drive term and the behaviour term through the same receptor.

---

## 5. Numerical integrity

Relative change after a 12-year run when dt is halved from 0.125 to 0.0625 d:

| Output | dt=0.125 | dt=0.0625 | Relative difference |
|---|---|---|---|
| %fat | 26.88183 | 26.86764 | 5.3e-04 |
| HQ-CT | 9.95914 | 9.95957 | 4.3e-05 |
| IGF-1 | 416.24969 | 416.33118 | 2.0e-04 |
| AHI | 3.33654 | 3.33304 | **1.1e-03** |
| Height | 152.98526 | 152.98948 | 2.8e-05 |
| SAT | 0.45622 | 0.45621 | 2.7e-05 |

Maximum drift 1.1e-03. Energy is conserved exactly by construction — because lean
mass follows a developmental target and **fat is the energy buffer**.

---

## 6. Files

| File | Contents |
|---|---|
| [`pws_qsp_model.dot`](pws_qsp_model.dot) | Mechanistic map source (20 clusters, with reading-order comments) |
| [`pws_qsp_model.svg`](pws_qsp_model.svg) · [`.png`](pws_qsp_model.png) | Rendering (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`pws_mrgsolve_model.R`](pws_mrgsolve_model.R) | 64-ODE mrgsolve model + 15 scenario drivers |
| [`pws_shiny_app.R`](../../../prader-willi-syndrome/pws_shiny_app.R) | 10-tab dashboard (the mechanism tabs come before the endpoint tabs) |
| [`pws_reference_model.py`](../../../prader-willi-syndrome/pws_reference_model.py) | Dependency-free Python RK4 reference implementation |
| [`pws_reference_output.txt`](../../../prader-willi-syndrome/pws_reference_output.txt) | The source of every number in the document above |
| [`pws_calibration.py`](../../../prader-willi-syndrome/pws_calibration.py) | Closed-form derived quantities · branch structure · sensitivities |
| [`pws_calibration_output.txt`](../../../prader-willi-syndrome/pws_calibration_output.txt) | Its output |
| [`pws_references.md`](pws_references.md) | 357 papers, every PMID looked up in PubMed |

Reproducing:

```bash
python3 pws_reference_model.py > pws_reference_output.txt   # everything (~a few minutes)
python3 pws_reference_model.py --fast                       # quickly, with a coarse dt
python3 pws_calibration.py     > pws_calibration_output.txt
python3 pws_calibration.py --sens                           # including local sensitivities
dot -Tsvg pws_qsp_model.dot -o pws_qsp_model.svg
dot -Tpng -Gdpi=150 pws_qsp_model.dot -o pws_qsp_model.png
```

Because this was written in an environment without an R runtime, the mrgsolve file is
a **verbatim port of equations that were first run and verified in Python**. The six
real defects the Python runs exposed are marked at the point of correction with
`[DEFECT n]` (a missing normalisation that latched the control group into
hyperphagia, a Cunningham intercept that starved the neonate, a height of 8 m from
writing IGF-1 as a power, double-counting of lean mass, the error of using gastric
contents as a satiety signal on a daily scale, and the error of charging the activity
cost against total fat).

---

## 7. Limitations

- **Everything below the daily timescale has been averaged out.** The real half-lives
  of ghrelin, insulin, GH, and LH are on the order of minutes, and the model handles
  them as daily-scale pools. The 30-minute postprandial dynamics are not this model's
  question.
- **Drug concentrations are apparent daily-scale exposures**, fitted to the 24-hour
  AUC. They do not reproduce a true Cmax. Somatropin's `FPOTGH = 0.63` is a
  pulsatility correction.
- **Arousal and cognition are not written as state variables.** Daytime
  hypersomnolence is represented only through the airway and hypoxia pathways.
- **Mortality is not quantified.** Only the pathway is left on the mechanistic map, as
  a dashed edge.
- **The lean-mass effect of testosterone enters only through IGF-1, muscle tone, and
  the gonadal factor.**
- **The evidence level for the GLP-1 arm is low.** In PWS it is mostly observational,
  so it must not be read with the same weight as the carbetocin and DCCR arms.
- **Only three lines were calibrated** (`KGHRD`, `ECBMAX`/`EV1AMX`, `EDZMAX`). None of
  the remaining structural results was fitted to any clinical number — that boundary
  is tabulated in §0 of
  [`pws_references.md`](pws_references.md).

---

## ⚠️ Disclaimer

This model is a **qualitative and semi-quantitative QSP model for educational and
research purposes**. It was assembled on the basis of the published literature and
clinical trial data but has not been independently verified or certified, and **must
not be used directly for actual clinical decision-making, prescribing, or regulatory
submission.** The parameters and assumptions are approximations for the sake of
exposition, and fitting and validation against real patient data is separately
required.
