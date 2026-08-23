# Essential tremor — QSP model

**Essential Tremor · ET** — a quantitative systems pharmacology model based on an oscillator formulation

<a href="et_qsp_model.svg"><img src="et_qsp_model.png" width="560" alt="ET QSP map"></a>

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`et_qsp_model.dot`](et_qsp_model.dot) · [SVG](et_qsp_model.svg) · [PNG](et_qsp_model.png) | 201 nodes · 272 edges · 18 clusters |
| mrgsolve ODE model | [`et_mrgsolve_model.R`](et_mrgsolve_model.R) | 48 ODEs · 25 scenarios |
| Independent verification implementation | [`et_verify.py`](../../../essential-tremor/et_verify.py) · [output](../../../essential-tremor/et_verify_output.txt) | Pure Python RK4, no dependencies |
| Shiny dashboard | [`et_shiny_app.R`](../../../essential-tremor/et_shiny_app.R) | 10 tabs |
| References | [`et_references.md`](../../../essential-tremor/et_references.md) · [verification script](et_reference_check.py) | 126 papers, every PMID verified |

---

## 1. The single idea of this model

Tremor is **not the level of some substance but a limit cycle**. So the state variable of
this model is not "tremor" but the **amplitude envelope** of an oscillator just past a
supercritical Hopf bifurcation.

```
G_total = G0·(1+PROG) · [ w_C·Φ_C + w_P·Φ_P ]
Φ_C     = ( a_O·φ_olive + a_R·φ_cblthal ) · φ_thal · φ_ctx
Φ_P     = φ_spindle · φ_nmj
μ       = G_total − 1                            ← the bifurcation parameter
dr/dt   = (1/τ_A)·( r·(μ − β·r²) + ε ) / (1+|μ|)
r*      = √(μ/β)   (μ>0)      ·   ≈0   (μ<0)
```

Three things follow from these five lines **as arithmetic rather than as assumptions**.

**(1) Amplitude is set by gain, frequency by delay.** Every drug in this file enters only
the `φ_*` terms (gain) and touches `τ_loop` (delay) not at all. Consequently **every drug
changes only the amplitude and leaves the frequency alone** — as is actually observed.
Frequency instead falls as ageing and Purkinje loss lengthen the loop (model: 5.01 Hz at
age 60 → the 4.5 Hz range 20 years later).

**(2) Because of the √ law, the same drug abolishes tremor in mild disease and barely
touches it in severe disease.** Since the curve has a vertical tangent at G=1, a patient
near the threshold is pushed below it and the tremor disappears, while a patient far from
it moves a little along the flat part of the curve.

**(3) The ceiling on treatment is set by topology, not by potency.** The inferior olive is
**one of two parallel branches** of the central loop (weight `a_O`), whereas the Vim
thalamus is a **serial factor**. So even a perfect Cav3.1 blocker cannot bring Φ_C below
`a_R`, whereas a Vim lesion multiplies the whole central term by 0.08.

A fourth consequence is about **measurement** rather than biology. Because clinical scales
are the **logarithm** of amplitude (`grade = 2 + 2·log₁₀(A_cm)`, Elble), a 50% reduction
on accelerometry is always 0.60 points per item wherever the baseline sits.

---

## 2. The principal results the model derived (all values from running `et_verify.py`)

### 2.1 The same drug, the same occupancy, four patients — deriving the responder/non-responder split

Propranolol LA 160 mg, 24 weeks. β₂ occupancy is **identical** at 0.939 in all four
patients.

| Phenotype | G0 | Baseline amplitude | On treatment | Change in amplitude | Baseline TETRAS | Change in TETRAS |
|---|---|---|---|---|---|---|
| Mild | 1.15 | 0.977 cm | 0.015 cm | **−98.4 %** | 19.9 | **−19.65** |
| Moderate | 1.60 | 1.917 cm | 1.247 cm | −35.0 % | 26.0 | −3.87 |
| Severe | 6.00 | 5.504 cm | 4.729 cm | −14.1 % | 42.2 | −1.54 |
| Very severe | 12.0 | 8.160 cm | 7.122 cm | −12.7 % | 46.9 | −1.36 |

The content of this table is that the responder/non-responder dichotomy repeatedly
reported in β-blocker trials is **not two kinds of biology but a question of whether one
equation crosses a threshold**. A testable prediction: the response rate should be
inversely correlated with baseline severity.

### 2.2 The 'discrepancy' between accelerometry and TETRAS is not a discrepancy but a logarithm

On propranolol 160 mg, **accelerometry −35.0 %**; the same result is **TETRAS −3.87
points** (−14.9 % of a 26.0 baseline), and **−0.41 points** on a single upper-limb item.
To gain one point the amplitude has to fall **3.16-fold**.

This is exactly the pattern in the literature whereby drug trials report 4–6 points on
TETRAS/FTM while accelerometry reports 40–60% reductions, and **it was derived from the
log transformation rather than fitted**. Reporting "% change in TETRAS" is therefore a
scale error.

### 2.3 β blockade: the site is peripheral, and that is why β₁ selectivity fails

| Regimen | β₂ occupancy | Change in amplitude | Change in TETRAS | Heart rate |
|---|---|---|---|---|
| Propranolol LA 60 mg | 0.833 | −24.6 % | −2.56 | 57.7 |
| Propranolol LA 120 mg | 0.917 | −32.2 % | −3.51 | 54.3 |
| Propranolol LA 160 mg | 0.939 | −35.0 % | −3.87 | 53.2 |
| Propranolol LA 240 mg | 0.961 | −38.4 % | −4.36 | 52.0 |
| Propranolol LA 320 mg | 0.972 | −40.7 % | −4.70 | 51.4 |
| **Atenolol 100 mg** (β₁ selective) | 0.353 | **−4.3 %** | −0.41 | 51.5 |
| **Nadolol 120 mg** (non-selective · peripherally restricted) | 0.991 | **−38.2 %** | −4.30 | 50.2 |

The three drugs lower heart rate to much the same degree (51–58 bpm) while the tremor
effect diverges from −4 % to −38 %. Because the benefit comes **only from β₂ occupancy**.
And nadolol, which barely enters the brain, produces effectively the same result as
propranolol — by the model's decomposition, **93 % of propranolol's effect at 160 mg is
peripheral β₂ and the central non-β component is only −2.4 percentage points**. What is
gained by raising the dose to 240–320 mg is that small central component too (β₂ occupancy
is already saturated at 0.917 at 120 mg), and fatigue and central adverse effects attach
over the same range.

### 2.4 The contraindication is the efficacy

`OCCB2` appears in **exactly two places** in the model: spindle gain and the airway. Hence
the same number.

| Drug | β₂ occupancy | Change in amplitude | FEV₁ (with comorbid asthma) |
|---|---|---|---|
| Propranolol 160 | 0.939 | −35.0 % | 2.18 L (**−32.0 %**) |
| Nadolol 120 | 0.991 | −38.2 % | 2.12 L (−33.7 %) |
| Atenolol 100 | 0.353 | −4.3 % | 2.81 L (−12.2 %) |

Buy tremor benefit through β₂ blockade and the airway effect attaches **at an almost fixed
ratio**. The asthma contraindication is not a fact to be memorised separately; it is the
same term as the efficacy.

### 2.5 Primidone — the active molecule is the parent

Primidone t½ ≈ 10 h, phenobarbital t½ ≈ 100 h. The model was given **only two potencies
and two half-lives**, and the time course is an output.

| Dose (q8h) | Change in amplitude | P_parent | P_phenobarbital | Parent's contribution | Somnolence |
|---|---|---|---|---|---|
| 250 mg/day | −18.0 % | 0.206 | 0.054 | **79.2 %** | 18.9 |
| 500 mg/day | −25.5 % | 0.310 | 0.099 | 75.8 % | 28.2 |
| 750 mg/day | −31.6 % | 0.371 | 0.136 | 73.2 % | 33.7 |

At 12 hours on the first day of 250 mg/day, most of the steady-state effect is already
present (−13.8 % against −15.0 %), and at that point phenobarbital contributes **1 %** of
the total GABA-A potentiation.
And **reaching the same potentiation as primidone 250 mg/day with phenobarbital alone
requires a plasma level of 66.9 mg/L** — the therapeutic range is 10–40 mg/L. This is
where the reason phenobarbital is useless as a tremor drug is derived.

### 2.6 Combination therapy is **supra-additive** — a result that refuted the author's prior expectation

The author expected the combination to be **sub-additive** (each drug shaves the gain and
the amplitude is a square root). `et_verify.py` refuted this: because **√ is concave**,
pushing μ towards 0 moves into the region where the mapping is steepest, and the second
drug earns **more** than it does alone.

```
propranolol alone −35.0 %  |  primidone alone −18.0 %  |  simple sum −52.9 %
actual combination        −70.8 %  →  17.8 percentage points greater than the sum
gain: G 1.612 → 1.260 (propranolol) / 1.412 (primidone) / 1.065 (combination)
```

This is consistent in direction with propranolol + primidone being the clinical standard
and occasionally producing dramatic responses. It is presented as **a testable
prediction**, and recorded not as a result that came out as expected but as one where the
expectation was wrong.

### 2.7 Ethanol: suppression and rebound are not additional assumptions but consequences of asymmetric dynamics

Acute tolerance rises fast (Mellanby, τ_on 1 h) and decays slowly (τ_off 5 h). Ethanol
itself is gone in 4–5 hours. So **the adaptation necessarily outlasts the drug**, and that
is the rebound.

2 standard drinks (28 g): nadir **−34.1 %** (2.0 h) → peak rebound **+17.9 %** (8 h) →
+1.0 % at 24 h

| Standard drinks | Nadir | Time of nadir | Rebound | Intoxication index |
|---|---|---|---|---|
| 1 | −20.4 % | 1.50 h | +10.6 % | 20.1 |
| 2 | −34.1 % | 2.00 h | +19.5 % | 36.5 |
| 3 | −42.5 % | 2.26 h | +23.2 % | 47.8 |
| 4 | −47.9 % | 2.76 h | +25.0 % | 55.8 |

**The self-medication trap**: on 3 drinks a day for 90 days, the morning tremor with no
alcohol on board rises from +18.1 % on day 1 to **+22.1 %** on day 89, and the evening
effect is eroded from +1.1 % to +8.3 %. This is a mechanistic account of the raised risk
of alcohol use disorder in ET patients, and it constitutes quantitative grounds for
telling a patient self-treating with alcohol rather than a prescription that **the tremor
itself gets worse**.

**1-octanol** has a different efficiency per unit of intoxication because of the
difference in cerebellar/cortical distribution: 2 drinks of ethanol give a gain/intoxication
ratio of **0.95**, and octanol 8 mg/kg gives **4.81** — about 5 times the tremor benefit
for the same intoxication. (The absolute effect of octanol is smaller: −15.1 % against
−34.8 %.)

### 2.8 T-type calcium blockers — the ceiling is made by topology (this model's headline)

| Setting | φ_olive | Φ_C | G | Change in amplitude |
|---|---|---|---|---|
| Human ET, 100 mg/day | 0.536 | 0.838 | 1.450 | −14.4 % |
| Human ET, **perfect** Cav3 blockade | 0.011 | 0.654 | 1.273 | **−33.2 %** (the ceiling) |
| **Harmaline rat** (a_O = 1), the same 100 mg/day | 0.536 | 0.536 | 1.235 | **−60.5 %** |

The same drug, the same channel occupancy; what changed is **one parameter (`a_O`)**.
Harmaline oscillation is by definition purely olivary (a_O=1, serial) so blockade shaves
the total gain directly, whereas in humans the olive is one of two parallel branches, so
even blocking it perfectly leaves `a_R` behind.

An `a_O` sweep (the maximum reachable under perfect blockade):

```
a_O = 0.20 → −17.3 %   0.35 → −33.0 %   0.50 → −54.0 %
a_O = 0.65 → −98.2 %   0.80 → −98.4 %   1.00 → −98.4 %
```

The threshold lies between a_O ≈ 0.50 and 0.65. **The very fact that T-type blockers were
not curative in humans therefore gives an upper bound of a_O < 0.62** — a failed clinical
trial measured a parameter. It explains why a target that was curative in the animal model
stalled at a moderate effect in humans by **circuit topology** rather than by potency or
exposure, and it is falsifiable.

### 2.9 The reason surgery beats drugs is topology (and volume produces two outcomes)

| Lesion volume | Effective lesioned fraction | G | Change in amplitude | Ataxia index |
|---|---|---|---|---|
| 20 mm³ | 0.116 | 1.503 | −9.4 % | 7.1 |
| 40 mm³ | 0.427 | 1.275 | −32.9 % | 13.3 |
| 60 mm³ | 0.672 | 1.040 | −74.4 % | 18.8 |
| **90 mm³** | 0.850 | 0.870 | **−98.4 %** | 25.7 |
| 120 mm³ | 0.921 | 0.802 | −98.4 % | 31.6 |
| 250 mm³ | 0.986 | 0.739 | −98.4 % | 49.0 |
| 400 mm³ | 0.996 | 0.730 | −98.4 % | 60.6 |

Efficacy half-saturates at V50 = 45 mm³ while ataxia keeps rising to V50 = 260 mm³.
**Above 90 mm³ the additional benefit is zero and the additional ataxia keeps
accumulating** — the volume therapeutic window of MRgFUS was derived rather than chosen.

DBS frequency is not a dial but a **switch** (Hill 4, f₅₀ 80 Hz):

```
10 Hz  +17.6 %   30 Hz  +9.7 %   50 Hz  −3.9 %   80 Hz  −25.5 %
100 Hz −43.5 %   130 Hz −63.5 %  185 Hz −82.7 %  250 Hz −94.0 %
```

Even **the worsening of tremor** at low frequency (entrainment) is reproduced, and the
clinical ">100 Hz" rule comes out of a single Hill coefficient.

**Habituation is a matter of rewiring capacity, not of elapsed time.** Five years after a
120 mm³ lesion:

| Rewiring capacity | 1 year | 3 years | 5 years |
|---|---|---|---|
| 0.45 (typical) | −98.4 % | −48.5 % | −29.0 % |
| 0.75 | −75.5 % | −26.7 % | −8.2 % |
| 0.95 | −61.2 % | −15.3 % | +3.4 % |
| 1.00 | −58.3 % | −12.6 % | **+6.1 %** |

The moment φ_thal crosses the threshold (0.375) the tremor returns **relatively abruptly**
— it is a bistable transition rather than a gradual worsening. This is consistent with the
pattern in which DBS tolerance is reported in a minority and often suddenly.

### 2.10 Botulinum toxin — the therapeutic window is made by precision, not by dose

**Separate SNAP-25 pools** are placed in the tremor-predominant muscles (wrist flexors and
extensors) and in the grip muscles (flexor digitorum profundus), and the only thing
separating them is the spillover fraction `f_spill`.

| Injection | SNAP_T | SNAP_G | Change in amplitude | Grip strength | QUEST |
|---|---|---|---|---|---|
| 100 U, ultrasound/EMG **guided** (f=0.15) | 0.455 | 0.870 | **−44.2 %** | **92.6 %** | 27.6 |
| 100 U, **unguided** (f=0.45) | 0.600 | 0.658 | −25.2 % | **68.6 %** | 37.1 |
| 50 U, unguided (f=0.45) | 0.774 | 0.811 | −10.0 % | 87.5 % | 34.2 |
| 150 U, guided (f=0.15) | 0.307 | 0.811 | −69.2 % | 87.5 % | 22.6 |

**Halving the dose merely reduces the tremor effect and the grip weakness together; it
does not open the window** (−10.0 % / grip 87.5 %). What opens the window is `f_spill`: a
guided injection reduces tremor **more** while losing **less** grip. This is where the
reason the early trials, designed with fixed doses and no guidance, failed on grip is
derived, and it is a clinically actionable conclusion.

### 2.11 The model reproduces the differential-diagnostic test itself

The mass loading test. Two diseases come out of the same equations.

| Condition | μ | f₀ (mechanical resonance) | Observed dominant frequency | Amplitude |
|---|---|---|---|---|
| ET, unloaded | +0.606 | 8.50 Hz | 5.01 Hz | 1.908 cm |
| ET, +500 g | +0.606 | 4.95 Hz | **4.85 Hz** | 2.154 cm |
| Enhanced physiological tremor (EPT), unloaded | **−0.443** | 8.50 Hz | 8.48 Hz | 0.513 cm |
| EPT, +500 g | −0.443 | 4.95 Hz | **4.95 Hz** | 0.513 cm |
| EPT + propranolol 160 | −0.587 | 8.50 Hz | 8.19 Hz | 0.031 cm |

**Loading moves the ET peak by −0.16 Hz and the EPT peak by −3.54 Hz.** Because ET is a
limit cycle with μ>0 whose frequency comes from the central delay (1/τ_loop), whereas EPT
has μ<0 and no limit cycle, being noise passed through the mechanical resonance, so its
frequency simply is f₀. This is exactly the discriminator used in the clinic and the
laboratory, and it was **derived from the structure**. That the same propranolol
essentially abolishes EPT (0.513 → 0.031 cm) while reducing ET by only a third is the same
reason.

### 2.12 Head tremor appears when the cervical effector crosses G=1

| Phenotype | μ_cervical | Head amplitude | Upper-limb amplitude |
|---|---|---|---|
| HDG 0.55, G0 1.6 | **−0.117** | 0.001° (absent) | 1.908 cm |
| HDG 0.85, G0 1.6 | +0.365 | 0.891° | 1.908 cm |
| HDG 0.85, G0 2.4 | +1.048 | 1.510° | 2.909 cm |
| HDG 0.85, G0 6.0 | +4.119 | 2.994° | 5.492 cm |

Head tremor is not a separate disease but **the event of the cervical effector crossing its
own threshold as G0 rises** — consistent with the observation that head tremor is a marker
of patients with long disease duration and severe disease. And the same propranolol gives
**−21.0 % in the arm and −9.7 % in the head**: because the β₂ spindle mechanism carries
less weight in the cervical loop (w_P,neck = 0.15), and this too is a derived result.

### 2.13 Progression: constant mechanistic progression produces a decelerating amplitude

`dA/dt = (dG/dt)/(2√μ)` is largest when μ is small.

| Year | G | μ | Amplitude | TETRAS | Annual rate of amplitude increase |
|---|---|---|---|---|---|
| 0 | 1.606 | 0.606 | 1.908 cm | 25.9 | — |
| 1 | 1.687 | 0.687 | 2.027 cm | 26.5 | **+6.3 %/year** |
| 2 | 1.783 | 0.783 | 2.160 cm | 27.0 | +6.5 %/year |
| 5 | 2.094 | 1.094 | 2.536 cm | 30.1 | +5.8 %/year |
| 10 | 2.653 | 1.653 | 3.086 cm | 35.2 | +4.3 %/year |
| 20 | 3.912 | 2.912 | 4.019 cm | 39.1 | **+2.5 %/year** |

While the gain rises linearly, **the proportional rate of worsening slows 2.5-fold**. This
is the same direction as patients' accounts of worsening quickly in the first few years
after onset and then flattening out for a long time.

### 2.14 An explicit hypothesis: can early suppression be disease-modifying (default OFF)

`KEXC` is the term by which oscillation amplitude feeds back onto cerebellar damage. **Its
default is 0, and it is reported as a hypothesis rather than a result.**

| At 10 years | PROG | G | Amplitude | TETRAS |
|---|---|---|---|---|
| KEXC 0, untreated | 0.500 | 2.653 | 3.086 cm | 35.2 |
| KEXC 0.3, untreated | 0.690 | 3.092 | 3.446 cm | 36.9 |
| KEXC 0.3, propranolol for 10 years | 0.608 | 2.322 | 2.746 cm | 34.2 |

If the hypothesis is true, ten years of suppression reverses 43 % of the progression.
**This is an open question, and the model does not answer it; it only says how large the
answer would have to be to be detectable.**

---

## 3. Verification

`et_verify.py` is **not** a convenience wrapper around the R model but a second
implementation rewritten from scratch so as to disagree if either one is wrong (pure
Python, hand-written RK4, no numpy or scipy).

### 3.1 Is the slow-envelope reduction legitimate?

The full two-dimensional oscillator was integrated at 5.5 Hz and its peak amplitude
compared with `r*`.

| μ | Peak of the full oscillator | r* (envelope) | Relative error |
|---|---|---|---|
| 0.15 | 0.3949 | 0.3873 | 1.97 % |
| 0.60 | 0.7784 | 0.7746 | 0.50 % |
| 1.40 | 1.1857 | 1.1832 | 0.21 % |
| 5.00 | 2.2374 | 2.2361 | **0.06 %** |

### 3.2 The seven real defects this verification found and fixed

They are marked `[FIXED — defect N]` at the relevant places in the file.

| # | Defect | Symptom | Fix |
|---|---|---|---|
| 1 | The envelope relaxation time was `τ_A/(2μ)` — it **shrinks** as severity rises | RK4 returned **NaN** at G0=12 | Divided by `(1+|μ|)`. The fixed point `r*` is unchanged; the relaxation time is bounded above by τ_A/2 |
| 2 | `NZ = 1+KCAT(RAG−1)` | **Physiological tremor of −0.73** (negative) on propranolol | Only the β₂-dependent fraction of physiological tremor can be removed (`FNZ`), floor 0.30 |
| 3 | `B2REG` multiplied the **agonist concentration** | Atenolol **worsened** tremor while on treatment (+2.5 %) | Being a receptor-number effect, moved to the β₂-mediated **gain term** → appears only in withdrawal |
| 4 | `REB = 1+KREBF·ADAPTF`, a single symmetric τ | Destroyed the acute suppression by itself (φ_cbl>1 at the peak) and the **rebound was +0.1 %** (effectively absent) | The rebound is the **non-opposed part** of the adaptation (`ADAPTF − KAF·P_RAW`), asymmetric with τ_on 1 h / τ_off 5 h |
| 5 | `R_UL²` in the progression driving term | **OverflowError** when the envelope diverged | Bounded above |
| 6 | `KD_AG = 30 nM` → a resting spindle β₂ occupancy of 0.0066 carrying 60 % of the peripheral gain | RAG had 150-fold headroom → **a hyperthyroid "EPT" patient at μ=1.41**, i.e. possessing a limit cycle (in flat contradiction with the definition of EPT) | `KD_AG=0.35 nM`, resting occupancy 0.36, RAG capped at 2.75 |
| 7 | A fixed damping **coefficient** → ζ falls as 1/√J | Resonance sharpened at +500 g and the ET amplitude **doubled** (1.91→3.70 cm) | Fixed the damping **ratio** instead, reflecting co-contraction under load |

In addition, this verification **refuted the author's prior hypothesis** (§2.6: the
combination is supra-additive, not sub-additive). And it caught a defect on the runner side
as well: the large step needed for multi-year runs (dt=1 h) is 5.7 times the envelope
relaxation time (0.175 h), so **the 20-year progression curve came out non-monotonic**
(2.16 cm at 2 years, 2.09 cm at 5 years — while μ increases monotonically). It was fixed
to pin the envelope at quasi-steady state for dt ≥ 0.2 h. mrgsolve is unaffected because
LSODA reduces its own step.

### 3.3 The relationships between the scales were verified too

Reading the trough alone **underestimates** a once-daily drug (propranolol 160 mg: trough
−31 % against a dosing-interval mean of −37 %). Because a trial's assessment time is not
the trough, both this file and the R model report **the final dosing-interval mean**.

---

## 4. Where the model disagrees with the literature (reported as it stands, not reconciled)

Items where it was judged better to write them down than to hide them.

1. **Primidone comes out weak.** −31.6 % at 750 mg/day, whereas the strongest reports
   suggest reductions of 50–60 %. The relative ordering against propranolol is right
   (roughly equivalent).
2. **Baseline TETRAS comes out 3–6 points high.** The mild patient is at 19.9 points,
   whereas the mild-to-moderate baseline in real trials is roughly 15–25. The per-item log
   anchor (1 cm ≈ 2 points) is straight from the literature, so the difference comes from
   the reconstruction approach of deriving all 12 items from amplitude.
3. **The ET amplitude increases by +13 % under load** (the 2.5.11 table). In the clinic it
   is reported to change hardly at all, or to fall slightly. Fixing the damping ratio
   removed the 2-fold error but the sign remains. It is left as a falsifiable prediction.
4. **Atenolol comes out almost ineffective at −4.3 %.** Some trials report modest efficacy
   for atenolol. The model depends entirely on the selectivity assumption of a β₂ Ki of
   1000 nM.
5. **The magnitude of the ethanol rebound (+18 to 25 %) has not been validated
   quantitatively.** The direction and time course (a 4–10 hour window) agree with reports
   but the magnitude rests on `KREBF` alone.
6. **There is no placebo effect in the model.** The placebo response in ET trials is not
   small, so every change here is an absolute value against no treatment rather than
   against placebo.
7. **`w_P = 0.40` (the peripheral loop carrying 40 % of the gain) is a strong claim.** It
   was back-calculated from nadolol's efficacy and propranolol's effect size, but it is not
   an independently measured value.

---

## 5. Files and how to run them

```bash
# render the map
dot -Tsvg et_qsp_model.dot -o et_qsp_model.svg
dot -Tpng -Gdpi=150 et_qsp_model.dot -o et_qsp_model.png

# independent verification (no dependencies; run without --quick for a finer step)
python3 et_verify.py --quick

# re-verify the references (NCBI E-utilities)
python3 et_reference_check.py --verify
python3 et_reference_check.py --harvest
python3 et_reference_check.py --emit
```

```r
# the mrgsolve model + the 25 scenarios
source("et_mrgsolve_model.R")

# Shiny dashboard (10 tabs)
shiny::runApp("et_shiny_app.R")
```

### Compartment composition (48 ODEs)

| Group | Number | Contents |
|---|---|---|
| PK · effect sites | 26 | propranolol (2 compartments + brain) · atenolol · nadolol · primidone/phenobarbital/PEMA (+2 brain) · topiramate · gabapentin · ethanol (MM elimination + brain) · 1-octanol · T-type blocker (+brain) · botulinum 2 depots |
| Target states | 2 | 2 pools of functional SNAP-25 by muscle |
| Adaptation · induction | 5 | acute/chronic GABA-A opposing adaptation, β receptor upregulation, enzyme induction, tolerance to somnolence |
| Disease | 4 | Purkinje integrity, dentate disinhibition, thalamic rewiring, gain progression |
| **Oscillators** | **3** | Upper-limb · cervical · laryngeal amplitude envelopes |
| Organ systems | 6 | FEV₁, HCO₃⁻, body weight, bone mineral density, ALT, cognition |

---

## 6. Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It was
constructed starting from the public literature but has not been fitted or validated
against patient data. **It must not be used for clinical decision-making, prescribing, or
regulatory submission.** The parameters are illustrative approximations, and every number
in §2 is **an output of this model**, not a literature value.
