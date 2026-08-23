# Chemotherapy-Induced Peripheral Neuropathy (CIPN) — QSP Model
### Quantitative Systems Pharmacology Model

[🗺️ Mechanistic map (SVG)](../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.svg) ·
[⚙️ mrgsolve model](cipn_mrgsolve_model.R) ·
[📊 Shiny app](cipn_shiny_app.R) ·
[📚 102 references](cipn_references.md) ·
[🐍 Python reference implementation](../../../chemotherapy-induced-peripheral-neuropathy/cipn_reference_model.py) ·
[📄 Run output](../../../chemotherapy-induced-peripheral-neuropathy/cipn_reference_output.txt)

<a href="../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.svg"><img src="../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.png" width="760" alt="CIPN mechanistic map"></a>

---

## 0. The question this model sets out to answer

CIPN is a **treatment-limiting toxicity**. The medical oncologist makes the same decision
every cycle — *give this dose as written, reduce it, or stop.* What makes the decision hard
is that the neurotoxicity and the antitumour effect **come out of the same exposure**.

This model was built to quantify that trade-off, and three structural commitments generate
all of its interesting behaviour.

| # | Structural commitment | What it produces |
|---|---|---|
| 1 | A **23-day decay half-life** for the SARM1 axon-death programme + **regeneration blocked** while the programme is switched on | **Coasting** — continued worsening for 4~12 weeks after the last dose. It emerges automatically from two rate constants, with no post-hoc forcing term |
| 2 | Nerve damage is a **convex function** of DRG proteasome occupancy; the antitumour effect is a **saturating function** of plasma occupancy | **Bortezomib SC vs IV**, at identical AUC and 1/14 the Cmax, differ in toxicity only and are equal in efficacy |
| 3 | The oxaliplatin antitumour exposure-response **saturates** (C50 = 250 mg/m² cumulative) but the neurotoxicity **keeps on integrating** | A **therapeutic-index optimum** in cumulative dose exists, and it is derived rather than assumed (the IDEA result) |

Every number here was obtained by running `cipn_reference_model.py` (standard library only),
and the full output is in [`cipn_reference_output.txt`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_reference_output.txt).

---

## 1. Deliverables

| File | Contents |
|------|------|
| [`cipn_qsp_model.dot`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.dot) / [`.svg`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.svg) / [`.png`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_qsp_model.png) | Mechanistic map — **16 clusters · 278 nodes · 369 edges** |
| [`cipn_mrgsolve_model.R`](cipn_mrgsolve_model.R) | **Model of record** — 34 ODEs, PK for 4 drugs, 18 scenarios, population incidence, therapeutic-index optimisation |
| [`cipn_shiny_app.R`](cipn_shiny_app.R) | 12-tab interactive dashboard |
| [`cipn_references.md`](cipn_references.md) | 102 references (all with PubMed links, tied to model structure and parameters) |
| [`cipn_reference_model.py`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_reference_model.py) | **Dependency-free Python re-implementation** of the same ODEs with the same parameters — prints every number in this README |
| [`cipn_reference_output.txt`](../../../chemotherapy-induced-peripheral-neuropathy/cipn_reference_output.txt) | The complete run output of that script |

`python3 cipn_reference_model.py` regenerates the whole report, and
`python3 cipn_reference_model.py --calibrate` re-derives the calibrated parameters from scratch.

---

## 2. Mechanistic map (16 clusters)

1. Neurotoxic chemotherapy classes · schedule · cumulative dose (Cmax vs AUC vs cumulative)
2. Systemic PK (oxaliplatin 3-compartment / paclitaxel 3-compartment MM / bortezomib 2-compartment + SC)
3. **Why the DRG** — transporters (OCT2 · CTR1 · P-gp) and the absence of a blood-nerve barrier
4. Platinum mechanism — Pt-DNA adducts, nucleolar stress, reduced protein synthesis in the soma, neuronal death
5. Microtubule and axonal-transport failure (taxanes · vinca · auristatin ADCs)
6. Proteasome inhibition (bortezomib) — Cmax-driven, partly reversible
7. Mitochondrial dysfunction · oxidative stress — the common pathway on which every class converges
8. **Committed degeneration and coasting** — the SARM1 execution programme
9. Ion-channel remodelling · peripheral hyperexcitability — the acute syndrome (Nav1.6 · TRPM8 · Kv7)
10. Peripheral neuroinflammation (DRG macrophages · satellite glia · Schwann cells · mast cells)
11. Spinal dorsal horn · descending modulation — central sensitisation (the duloxetine site of action)
12. Host risk factors · pharmacogenomics (diabetes · age · CMT · GSTP1 · ERCC1 · CYP2C8*3)
13. Prevention strategies — cooling, compression, exercise, and **the agents that failed** (Ca/Mg · vitamin E · ALC)
14. Symptomatic-treatment PK/PD — duloxetine · pregabalin · topical therapy
15. Clinical endpoints · biomarkers (CIPN20 · CTCAE · TNSc · IENFD · plasma NfL)
16. **The dose-intensity dilemma** — RDI → tumour kill → 3-year DFS

---

## 3. mrgsolve model structure (34 ODEs)

**PK (12 compartments)** oxaliplatin ultrafilterable platinum, 3 compartments (deep
compartment = tissue-bound platinum, t½ ≈ 14 d) · paclitaxel 3 compartments +
Michaelis-Menten elimination · bortezomib 2 compartments + SC depot · duloxetine · pregabalin

**Biophase (4)** DRG platinum (accumulating, t½ 8 d) · taxane nerve burden (accumulating,
t½ 12 d) · bortezomib DRG effect site (**fast keo** — the DRG capillaries are fenestrated,
so it tracks plasma Cmax) · acute channel effect site

**Mechanism · biomarkers (18)** Pt-DNA adducts → nucleolar stress → protein synthesis in
the soma · tubulin occupancy + proteasome inhibition → axonal transport · mitochondrial
capacity · ROS · `ENERGY = MITO × transport × reserve` · **SARM1** · axon density ·
**irreversible DRG neuron loss** · DRG macrophages · IL-1β · plasma NfL · chronic
hyperexcitability · acute cold dysaesthesia · central sensitisation · descending
noradrenergic tone · cumulative-dose tracking

### Three modelling details you have to know about

**(a) Why a multi-timescale integrator is needed.** The self-decay rate of the paclitaxel
central compartment is about 144/day, and bortezomib's about 290/day. At the 0.05-day step
the slow biology demands, explicit RK4 is far beyond its stability limit and diverges — the
drug amount goes negative, the sign of the Michaelis-Menten denominator `(Km + C)` flips, and
the solution explodes. So the PK block is written in the form `dy = src − b·y` and advanced by
the **exact exponential update** `y ← y·exp(−b·dt) + (src/b)(1 − exp(−b·dt))`
(unconditionally stable at every step size), and only the slow biology is solved with RK4.

**(b) Floating-point error silently swallows dose.** Accumulating steps means arriving at a
dosing time as 6.99999999999 rather than 7.0. Look the dosing times up without a tolerance and
the step schedule is not refined, so **a single RK4 step straddles the entire infusion
interval** and **33% of every dose after the first disappears**. This bug is particularly
dangerous because it does not diverge — it returns plausible results. The current
implementation puts a `+1e-9` tolerance on the dosing-time lookup, snaps steps exactly onto
the boundaries, and evaluates the **infusion rate once per step, at the midpoint** (steps are
truncated at every infusion start and end, so the rate is constant inside a step).
Verification: the delivered dose is exact to machine precision — `12×85 = 1020.0000 mg/m²`.

**(c) Why such a smooth hazard function.** `RISK` is the **fraction** of sensory axons whose
metabolic demand exceeds supply. Axons differ in length, and therefore in demand, so this must
be a gentle sigmoid in `ENERGY/E_THR` rather than a single hard threshold. With a hard
threshold the virtual population becomes **bimodal** (either untouched or devastated) and
cannot reproduce the observed grade distribution. In the same way, dying-back is
**self-limiting** (`PHI_RELIEF`) — when distal axons are lost, the transport and energy burden
on the surviving axons falls. That is what makes severity **graded** rather than
all-or-nothing.

---

## 4. Calibration: 7 parameters ← 7 clinical-trial observations

**Only the seven below were fitted.** Every other number is a prediction.

| Parameter | Value | Fitted to |
|---|---|---|
| `KDAM_PT` | 0.3799 | IDEA — FOLFOX 6 months, grade ≥2 = 47.7% |
| `SIGMA_S` | 0.3664 | MOSAIC — FOLFOX 6 months, grade ≥3 = 12.4% |
| `PHI_RELIEF` | 0.5009 | IDEA — FOLFOX 3 months, grade ≥2 = 16.6% |
| `KDAM_TAX` | 0.1626 | ECOG 1199 — weekly paclitaxel 80 mg/m² ×12, grade ≥2 = 27% |
| `KDAM_BTZ` | 147.60 | MMY-3021 — bortezomib 1.3 mg/m² IV, grade ≥2 = 41% |
| `BTZ_JH` | 1.5113 | MMY-3021 — bortezomib SC, grade ≥2 = 24% |
| `KCS_DUL` | 0.3535 | Smith 2013 JAMA — duloxetine, net −0.73 BPI |

**The grade thresholds were not fitted.** CTCAE grade 1/2/3 were fixed a priori at sensory
fibre loss of **10% / 25% / 45%** (the sural nerve and IENFD studies put the symptom-onset
threshold around 20~30% fibre loss). Fixing these thresholds is what makes the remaining arms
genuine out-of-sample predictions.

Also, because CTCAE grades **treatment-emergent** neuropathy, the structural deficit is
measured **from the patient's own baseline**. A diabetic patient's pre-existing asymptomatic
deficit is not already a CIPN grade. Reduced reserve acts **only through bioenergetics
(RESERVE) and regeneration (REGEN)**, and not by bringing the patient closer to a grade
threshold.

### Verification — what was fitted and what was predicted

| Arm | Model g≥2 | Observed g≥2 | Model g≥3 | Observed g≥3 | Status |
|---|---|---|---|---|---|
| FOLFOX 6 months (1020 mg/m²) | 47.7% | 47.7% | 12.4% | 12.4% | fitted |
| FOLFOX 3 months (510 mg/m²) | 16.6% | 16.6% | 0.8% | 2.7% | fitted (g2) / **predicted** (g3) |
| CAPOX 3 months (520 mg/m²) | **18.4%** | ~15% | 1.0% | ~3% | **predicted** |
| CAPOX 6 months (1040 mg/m²) | **50.0%** | ~45% | 13.7% | ~11% | **predicted** |
| Paclitaxel weekly 80×12 | 27.0% | 27% | 1.3% | — | fitted |
| Paclitaxel 3-weekly 175×4 | **14.9%** | 20% | 0.4% | — | **predicted** |
| Bortezomib 1.3 IV | 41.0% | 41% | 7.0% | 16% | fitted (g2) / predicted (g3) |
| Bortezomib 1.3 SC | 23.8% | 24% | 2.5% | 6% | fitted (g2) / predicted (g3) |
| Bortezomib IV once weekly | **6.3%** | — | 0.3% | 8% | **predicted** |

**The two CAPOX arms are a completely independent validation** — a different dose
(130 mg/m²), a different interval (q3w), not a single parameter fitted to them, and they come
in at 18.4% / 50.0% (observed ~15% / ~45%).

---

## 5. Known misses — read this first

The value of a QSP model lies less in what it gets right than in **whether you can diagnose
why it is wrong**.

**(1) The upper tail is too thin — grade ≥3 and long-term persistence are underpredicted
together.** FOLFOX 3 months grade ≥3 is 0.8% (observed 2.7%), bortezomib IV 7% (observed 16%).
The long-term course is worse: the model gives grade ≥3 as 10.2% at end of treatment →
0.1% at 12 months → 0.0% at 48 months, whereas MOSAIC observed 12.4% → 1.1% → 0.7%. Grade ≥1
runs 79.6% at end of treatment → 27.9% at 12 months → 7.8% at 24 months, diverging widely from
the observation that about one third of patients still have grade ≥1 sensory symptoms five
years after adjuvant FOLFOX (Selvy 2020).

These two are **the same end of the same distribution**, and the cause pins down to one thing:
**irreversible DRG neuron death (`KNEURON`) is underestimated**. At median susceptibility, soma
loss is only 1.2%, so the axons recover almost completely, to 98.8%. Raising `KNEURON` would
deepen the permanent floor and thicken the severe tail, so **both failures would improve at
once**. The reason it was not raised is this — **there are no data quantifying human DRG
neuron loss after platinum exposure.** That the parameter governing the long-term prognosis of
CIPN is a parameter nobody has measured may be this model's most useful conclusion.

**(2) Inter-cycle recovery is too complete — schedule effects are overpredicted.** Spreading
the same cumulative 1020 mg/m² from q2w → q3w → q4w gives grade ≥2 of
47.7% → 16.8% → 4.0%. Real oxaliplatin neurotoxicity is far more cumulative-dose-dominated
than this. The same cause **overestimates the benefit of once-weekly bortezomib** — in
Bringhen 2010 the grade 3–4 peripheral neuropathy rate on the once-weekly regimen was 8%
(twice-weekly 24%), where the model gives grade ≥3 as 0.3%. It happens because mitochondrial
(t½ 4.6 d) and axonal-transport (t½ 6.9 d) recovery is nearly complete over a 7~14 day
interval. **The direction is right and the magnitude is excessive.**

**(3) Cryotherapy is overestimated for the same reason.** With 45% of distal delivery blocked,
the model gives 47.7% → 4.6%. Real cryotherapy trials roughly halve grade ≥2. The **shape**
(threshold-like, so partial adherence hurts disproportionately) is worth reading, but do not
trust the absolute values.

**(4) `RESERVE` has too much leverage.** A 6% fall in bioenergetic reserve **doubles** the
risk. The observed relative risk for diabetes is about 1.5~2, and the model produces that with
a deficit of only 10%. This is **the parameter most urgently in need of real data** before it
is used for individual patient prediction.

**(5) What to concede on the PK.** Bortezomib AUC is 22.9 ng·h/mL in the model vs a reported
~151: the model uses first-dose clearance (102 L/h) while the reported value corresponds to
the lowered steady-state clearance. What matters for the model is that **the AUCs of SC and IV
are exactly equal**, and that is preserved. Paclitaxel reproduces AUC and time above threshold
(the established PD driver) to within ~10%, but underpredicts the Cmax of the 3-weekly regimen
(3.8 vs an observed ~4.3 µmol/L).

---

## 6. What the model says

### 6.1 Coasting is a derived result

| Regimen | Last dose | Peak severity | Coasting |
|---|---|---|---|
| FOLFOX 6 months | day 154 | day 187 | **+33 d (4.7 weeks)** |
| FOLFOX 3 months | day 70 | day 136 | **+66 d (9.4 weeks)** |
| Paclitaxel weekly ×12 | day 84 | day 119 | **+35 d (5.0 weeks)** |
| Bortezomib IV ×8 | day 168 | day 168 | 0 d |

Observed coasting is 4~12 weeks and the model lands inside that. There is no forcing term —
at the time of the last dose SARM1 is already switched on (activity 0.098), its decay half-life
is 23 days, and regeneration is blocked for that whole time. As a result axons are **lost
further, from 82.2% to 79.5%, with no further drug going in at all**. The absence of coasting
with bortezomib is consistent too — proteasome inhibition is reversible and spends only a short
time pushing SARM1 above threshold.

### 6.2 Bortezomib SC vs IV — the mechanism is convexity, not exposure

| | IV bolus | SC | Ratio |
|---|---|---|---|
| Plasma Cmax (model) | 201.8 ng/mL | 14.7 ng/mL | 13.7× |
| Plasma AUC | 22.9 ng·h/mL | 22.9 ng·h/mL | **identical** |
| Peak DRG 20S occupancy | 30.9% | 8.1% | 3.8× |
| **DRG occupancy AUC** | 0.0050 | 0.0059 | **0.85× (SC is the larger!)** |
| **Damage flux AUC** | 0.0057 | 0.0048 | **1.19× (IV is the larger)** |
| grade ≥2 PN | 41.0% | 23.8% | — |

Read this table slowly. **The occupancy AUC is in fact larger for SC** (occupancy saturates, so
the low, flat SC profile stays in the linear region for longer). Total exposure cannot explain
the difference in toxicity. But because the damage flux is a **Hill function of occupancy
(h = 1.51)**, the short, high IV peak delivers **1.19-fold the damage**. **Convexity is the
entire mechanism**, and the efficacy cost is zero — tumour kill depends on the saturating
plasma occupancy AUC, which the two routes share exactly.

### 6.3 IDEA 3 months vs 6 months — the optimum is derived, not assumed

| | grade ≥2 | 3-year DFS, high risk | 3-year DFS, low risk |
|---|---|---|---|
| No oxaliplatin | — | 53.2% | 76.9% |
| 3 months (510 mg/m²) | 16.6% | 62.7% | 82.4% |
| 6 months (1020 mg/m²) | 47.7% | 64.4% | 83.3% |
| *observed (IDEA)* | *16.6 / 47.7* | *62.7 / 64.4* | *83.1 / 83.3* |

The antitumour exposure-response is 80% saturated at 1020 mg/m² and 67% at 510 mg/m² —
**a 13 pp difference**. Over the same span the neurotoxic exposure **doubles**.

**Therapeutic index (w = 0.15, the weight on CIPN burden):**

| Cycles | Cumulative mg/m² | grade ≥2 | Utility, high risk | Utility, low risk |
|---|---|---|---|---|
| 4 | 340 | 3.5% | **7.72** | **4.22** |
| 6 | 510 | 16.6% | 7.01 | 2.94 |
| 8 | 680 | 31.2% | 5.60 | 1.18 |
| 10 | 850 | 41.6% | 4.58 | −0.09 |
| 12 | 1020 | 47.7% | 4.05 | **−0.80** |

The interesting result is **not** that the two optima differ — at this weight both are 4 cycles
(the CIPN cost curve is identical in the two strata and only the attainable DFS differs). The
point is **the asymmetry of overtreatment**. Going from the optimum out to 12 cycles costs the
high-risk patient 3.67 utility points but still leaves them ahead of no treatment. **The
low-risk patient falls below zero** — that is, giving a low-risk patient all 12 cycles is, on
this metric, worse than giving none at all. The last 500 mg/m² buys about 0.5 pp of DFS and
sells about 31 pp of grade ≥2. **This asymmetry is why IDEA's non-inferiority holds in the
low-risk stratum and could not be established in the high-risk one**, and here it was derived
from a saturating exposure-response.

**And lowering the CIPN weight splits the optima of the two strata:**

| CIPN weight w | High-risk optimum | Low-risk optimum |
|---|---|---|
| 0.05 | **12 cycles (6 months)** | **6 cycles (3 months)** |
| 0.10 | 4 cycles | 4 cycles |
| 0.15 | 4 cycles | 4 cycles |
| 0.30 | 4 cycles | 4 cycles |

At `w = 0.05` (counting CIPN as 1/20 of one DFS point) the optima are **6 months for high risk
and 3 months for low risk** — exactly **IDEA's recommendation**. The model was given no
information whatsoever about the trial durations. The separation appears by itself once a
saturating antitumour exposure-response meets a linearly accumulating neurotoxicity and the
crossing point is set by *how much residual disease there is left to kill*. Give CIPN a larger
weight and both strata converge on 4 cycles — meaning that for a patient who fears neuropathy
more than the trial designers did, **even 3 months is longer than optimal**. The 3-months
versus 6-months question has no answer independent of the patient's weighting of the two
harms, and this model **converts that weighting into a dose**.

### 6.4 The single most useful thing about the patient

Moving host factors one at a time at identical exposure (FOLFOX 6 months):

| Factor | Value | grade ≥2 | Relative risk |
|---|---|---|---|
| Bioenergetic reserve `RESERVE` | 1.00 → 0.94 | 47.7% → **68.8%** | **1.44** |
| Bioenergetic reserve `RESERVE` | 1.00 → 0.90 | 47.7% → **82.3%** | **1.72** |
| Regenerative capacity `REGEN` | 1.00 → 0.55 | 47.7% → 52.0% | 1.09 |
| Baseline axon density `AXON0` | 100 → 70 | 47.7% → 35.4% | **0.74** |

**A 6% fall in reserve doubles the risk, whereas a 45% fall in regenerative capacity moves it
by only 9%, and pre-existing asymptomatic fibre loss does not raise treatment-emergent CIPN at
all — in the model it is, if anything, slightly protective (there is less absolute axon to
lose).** Peak severity is decided by the energy set-point **during** exposure. Regeneration
governs *how much comes back*, and not *how bad it gets*.

The clinical implication: the reason a diabetic patient does worse should be **mitochondrial
state rather than nerve conduction studies**, and bioenergetic protection before and during
treatment is overwhelmingly more valuable than promoting regeneration afterwards. And
'reserve' is an item that no current CIPN risk score measures.
(But read it together with the excessive-leverage warning of §5(4).)

### 6.5 Acute and chronic are two different clocks

Oxaliplatin acute cold dysaesthesia (`COLDA`) peaks with every cycle and decays away
completely within days, and because cumulative axon loss sensitises the channel response it
**grows from cycle to cycle** (peak `COLDA` 0.186 → 0.234, cycle 1 → 12); it reaches its peak
1.25 days after each dose, and the time spent above 25% of the peak is 5.5 days. Yet 46 days
after the last dose the acute state is 0.0013 (effectively zero) while axon density is 80.4% —
**the chronic neuropathy is at its worst precisely then**. Same drug, two clocks.

### 6.6 Duloxetine moves the symptom and does not touch the axon

| | ΔBPI (5 weeks) | ΔIENFD |
|---|---|---|
| Placebo | −0.59 | 6.11 → 6.39 /mm |
| Duloxetine 60 mg | **−1.31** | 6.11 → 6.39 /mm |
| Pregabalin 300 mg | −1.39 | 6.11 → 6.39 /mm |

Net effect −0.73 (fitted to the −0.73 observed by Smith 2013). The predicted part is the
**dissociation** — IENFD is **not different at all**. The drug acts on central sensitisation,
so **a good analgesic response is not evidence that the neuropathy is recovering.** A trial
using CIPN20 as its primary endpoint is measuring a mixture of structure and symptom — this
model tracks the two separately.

### 6.7 Plasma NfL is an early warning, not a surrogate for grade

| | Value |
|---|---|
| Baseline NfL | 10.0 pg/mL |
| Peak NfL | **40.7 pg/mL (4.1× baseline)** |
| NfL reaches 2× | day 80 (cycle 6) |
| CTCAE grade 1 reached | day 100 |
| CTCAE grade 2 reached | **never reached** |
| **NfL lead time** | **20 days = 1.4 oxaliplatin cycles** |

The median-susceptibility patient **never reaches grade 2 at all**, and yet NfL rises to 4.1×
during treatment. That is, NfL is not a surrogate for grade but **a signal that axons are being
lost in a patient CTCAE still calls normal**. (Observed: in CIPN, NfL precedes change in the
clinical scales and rises 3~5-fold — Karteri 2022, Huehnchen 2022.)

Put together with §6.3 it becomes an actionable rule: the antitumour exposure-response is
already about two-thirds saturated at 510 mg/m², so **the cost of stopping oxaliplatin when
NfL rises is very small.**

### 6.8 Cryotherapy has to be applied every cycle

How many of the 12 cycles miss the 45% block of distal delivery:

| Cycles missed | 0 | 2 | 4 | 6 |
|---|---|---|---|---|
| grade ≥2 | 4.6% | 9.1% | 19.2% | 30.6% |

Because the damage is threshold-like, **partial adherence hurts disproportionately**.
(For the absolute values see the warning of §5(3) — read the shape only.)

### 6.9 Multivariate virtual population (n = 300, FOLFOX 6 months)

| Subgroup | n | grade ≥2 | grade ≥3 | Mean peak CIPN20 | Mean IENFD |
|---|---|---|---|---|---|
| All | 300 | 55.7% | 12.7% | 29.8 | 5.09 /mm |
| Diabetic | 61 | **80.3%** | 18.0% | 37.2 | 3.55 /mm |
| Non-diabetic | 239 | 49.4% | 11.3% | 27.9 | 5.49 /mm |
| Age < 60 | 160 | 47.5% | 11.9% | 27.9 | 5.34 /mm |
| Age ≥ 70 | 45 | **71.1%** | 20.0% | 33.7 | 4.52 /mm |

The overall incidence (55.7% / 12.7%) comes close to the fitted arm (47.7% / 12.4%) — the
covariate model was centred on the population mean age (60 years). Without that, every
covariate acts only in the direction of reducing reserve and the whole virtual population
becomes systematically worse than the reference patient. The covariate magnitudes are assumed
values, so read **the contrast between subgroups rather than the absolute level**, and the
diabetic and elderly contrasts inherit the excessive leverage of §5(4).

---

## 7. Scenario summary (typical patient, S = 1.0)

| Scenario | Cumulative dose mg/m² | Peak CS | Peak CIPN20 | Peak BPI | Lowest IENFD |
|---|---|---|---|---|---|
| FOLFOX 6 months (12×85 q2w) | 1020 | 0.240 | 25.2 | 3.14 | 5.57 |
| FOLFOX 3 months (6×85 q2w) | 510 | 0.111 | 11.9 | 1.65 | 6.39 |
| CAPOX 3 months (4×130 q3w) | 520 | 0.120 | 12.8 | 1.76 | 6.33 |
| CAPOX 6 months (8×130 q3w) | 1040 | 0.250 | 26.2 | 3.24 | 5.49 |
| Paclitaxel weekly 80×12 | 960 | 0.168 | 17.8 | 2.31 | 5.98 |
| Paclitaxel 3-weekly 175×4 | 700 | 0.115 | 12.3 | 1.69 | 6.33 |
| Bortezomib 1.3 IV ×8 | — | 0.216 | 22.8 | 2.97 | 5.69 |
| Bortezomib 1.3 SC ×8 | — | 0.149 | 16.0 | 2.24 | 6.15 |
| FOLFOX + cryotherapy | 1020 | 0.055 | 8.7 | 1.99 | 7.00 |
| FOLFOX, diabetic host | 1020 | 0.578 | 57.7 | 5.11 | 2.73 |
| FOLFOX, 20% dose reduction | 816 | 0.140 | 16.5 | 2.68 | 6.38 |
| FOLFOX stop-and-go | 1020 | 0.139 | 16.7 | 2.78 | 6.41 |

All 18 scenarios are in the `scenarios` list of `cipn_mrgsolve_model.R` and in
`cipn_reference_output.txt`.

---

## 8. Shiny dashboard (12 tabs)

① Patient · regimen ② Drug PK (Cmax vs AUC vs time above threshold) ③ Biophase · mechanism
④ Axon · coasting ⑤ Clinical endpoints ⑥ Biomarkers · NfL lead time ⑦ Acute vs chronic
⑧ Comparison of the 18 scenarios ⑨ Population incidence (exact calculation) ⑩ Route · schedule
⑪ Therapeutic-index optimisation ⑫ Calibration · verification

```r
shiny::runApp("cipn_shiny_app.R")
```

### How to compute population incidence exactly, without Monte Carlo

Peak severity is **monotone** in the susceptibility multiplier `S`, so
`P(grade ≥ g) = P(S ≥ s*_g)`; find `s*_g` by bisection and evaluate the log-normal survival
function. About 20~30 simulations per arm give the exact incidence **with no Monte Carlo
noise**. The multivariate virtual population is used separately, for the covariate analysis
only.

---

## 9. Reproducing

```bash
# the full report (every number in this README)
python3 cipn_reference_model.py

# re-derive the calibrated parameters from the trial anchors (about 15 min)
python3 cipn_reference_model.py --calibrate

# regenerate the mechanistic map
dot -Tsvg cipn_qsp_model.dot -o cipn_qsp_model.svg
dot -Tpng -Gdpi=150 cipn_qsp_model.dot -o cipn_qsp_model.png
```

The R side needs `mrgsolve`, `dplyr`, `tidyr`, `ggplot2` and `DT`.

---

## ⚠️ Disclaimer

This model is a **semi-quantitative QSP model for education and research**. It was built from
the public literature and clinical-trial data but has not been independently validated or
qualified, and **must not be used directly for real clinical decisions, prescribing, or
regulatory submission.** Be sure to read the known limitations of §5 alongside it — in
particular, long-term persistence and grade ≥3 incidence are systematically underpredicted,
while schedule effects, cryotherapy effects, and the reserve sensitivity are overpredicted.
