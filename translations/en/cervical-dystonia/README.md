# Cervical Dystonia — QSP Model

> **A dose-response curve that flattens is not a drug that has saturated.
> It is a disease the needle can reach that has run out.**
>
> This plateau is usually read as *pharmacological saturation*, and therefore
> ignored. This model reads it as a **measurement**.

In cervical dystonia, the dystonic torque that chemodenervation can remove is only
the part the injectate **actually reached**. Writing that share as φ,

```
        φ  =  ρ  ·  Σ w_m                    (sum over the injected muscles)
              ↑        ↑
     fraction of      which muscles
     each muscle      were on the
     the needle       list
     reached
```

and for **any** botulinum toxin — at any dose, product, potency, dilution, or
interval — the dystonic torque load `L` satisfies

```
        L(t)  ≥  L_min  =  (1 − φ) · D_cen
```

**There is not a single drug parameter on the right-hand side.** No dose, no
potency, no serotype, no dilution, no interval. φ is a property of the injection
*plan* and of the patient's *anatomy*.

And φ is observable. If the side effects keep rising while only the efficacy curve
flattens, then it is not that the toxin has used something up but that the needle
has used up the tissue it can reach. Reading the plateau that way and fitting it
(A1 stage 1):

```
        ρ  =  0.546          Σw (standard 4 muscles) = 0.65   →   φ = 0.355
                             Σw (extended 6 muscles) = 0.87   →   φ = 0.475
```

> **A standard superficial injection never touches about 64 % of the dystonic
> drive.** That number was not assumed from anatomy; it was recovered from the
> *shape* of the published dose-response curve.

---

## Four levers, with values attached (A2)

The currency is the **mean TWSTRS improvement score over a steady-state 12-week
cycle**. It is the value the patient actually lives with, and the only measure on
which a treatment that wears off and a permanent one can be compared honestly.

| Lever | Where the ceiling moves to | Headroom | Share of the whole design space |
|---|---|---|---|
| Currently achieved — 240 U q12wk | 14.62 points | — | — |
| **1. More or better toxin** (from the current φ up to the ceiling) | 18.98 | **+4.37** | 16.0 % |
| **2. Needle placement accuracy** ρ 0.55 → 1.00 | 34.64 | **+15.66** | **57.6 %** |
| **3. Extending the target list** (adding semispinalis + obliquus capitis inferior) | 41.57 | **+6.93** | 25.5 % |
| **4. Blocking nerve-terminal sprouting** | 14.87 | **+0.25** | 0.9 % |

```
    dose axis      ...   4.37 points  —  and there it hits the ceiling
    geometry axis  ...  22.59 points  —  5.2× the dose axis
```

240 U q12wk **already extracts 77.0 % of everything obtainable from a perfect,
permanent, side-effect-free block of the same motor units.** Raising the dose
60-fold reaches 99.8 % — that is, *almost all of the small thing*, at the price of
a dysphagia probability of 13 % → 21 %.

### Three conclusions that follow

**(i) The dose axis is the worst-value axis, and it carries all the harm.**
It is not merely inefficient; it **has a ceiling**, and that ceiling is the
smallest of the four.

**(ii) Placement accuracy and the target list are not substitutes for dose but
preconditions for it.**
Applied with 240 U held fixed they gain little — list extension +2.95, perfect
placement +12.85, both +17.11. That is because they do not *climb towards* the
ceiling, they *move* it. Dose climbs; φ moves. **A trial that improves imaging
guidance only while holding total dose fixed is testing the wrong combination**,
and this is a concrete prediction as to why the effect in such trials reads as
marginal.

**(iii) That ceiling is itself clinically unreachable.**
A perfect block puts the probability of neck weakness / head drop at 10 % under
the standard plan and at **93.5 %** under the extended plan with perfect placement.
The reason is that the posterior muscles generating the dystonic torque are the
very muscles holding the head up. **The disease's antagonist is posture's prime
mover.** This is not a dosing problem but a geometric conflict, and no improvement
in toxin resolves it.

---

## Deliverables

| File | Contents |
|---|---|
| [`cdys_qsp_model.dot`](../../../cervical-dystonia/cdys_qsp_model.dot) · [SVG](../../../cervical-dystonia/cdys_qsp_model.svg) · [PNG](../../../cervical-dystonia/cdys_qsp_model.png) | Mechanistic map — 17 clusters / 167 nodes / 219 edges. Therapy nodes are coloured by **which term of the model they touch, not by drug class** |
| [`cdys_mrgsolve_model.R`](../../../cervical-dystonia/cdys_mrgsolve_model.R) | 70-ODE mrgsolve model (muscle 8×5 · diffusion 8 · immune 6 · central 4 · clinical 3 · oral drug PK 6 · accounting 3), 15 scenarios, a full set of read-out functions |
| [`cdys_shiny_app.R`](../../../cervical-dystonia/cdys_shiny_app.R) | 10-tab dashboard. **Tab 2 is the point of this app** — it draws the ceiling and lets you try to cross it with the dose slider |
| [`cdys_references.md`](cdys_references.md) | 87 references, **every one actually resolved against PubMed**. Each section states which part of the model it supports |
| [`cdys_resolve_refs.py`](../../../cervical-dystonia/cdys_resolve_refs.py) | The script that generated that reference list. A device for not writing PMIDs from memory |
| [`cdys_reference_check.py`](../../../cervical-dystonia/cdys_reference_check.py) | **An independent numpy/scipy port of the same equation system.** Every number in this README was computed here |
| [`cdys_reference_output.txt`](../../../cervical-dystonia/cdys_reference_output.txt) | The full output of that script (A0–A13), verbatim |

This repository's build environment has no R toolchain. Rather than publish an ODE
model that has never been integrated, the same system was ported a second time and
integrated with scipy LSODA, and every number below is output from
`python3 cdys_reference_check.py`. If the two ports disagree, one of them is wrong.

```bash
python3 cdys_reference_check.py            # everything (A0–A13)
python3 cdys_reference_check.py --only A2  # the core result only
python3 cdys_reference_check.py --list
```

---

## Operator classification

Every therapy is classified by **which term it touches**. This classification
matters because each class **hits a different ceiling**, and they do not substitute
for one another.

| Operator class | Term touched | Examples | Ceiling |
|---|---|---|---|
| **CHEMODENERVATION** | T_m of the muscles reached | BoNT/A · B, dose, potency, long-acting toxin | blocked by φ (+4.37) |
| **GEOMETRY** | ρ, Σw → **φ itself** | needle EMG / ultrasound guidance, muscle selection, selective peripheral denervation | **moves** the ceiling (+22.6) |
| **CENTRAL DRIVE** | D_cen, G_cen | GPi-DBS, trihexyphenidyl, baclofen | strong, but mostly surgical |
| **GATING** | RecInh, SurrInh | clonazepam, GABAergic drugs | adjunctive |
| **AFFERENT** | aff (spindle afferent) | BoNT's intrafusal blockade, sensorimotor retraining | the only peripheral route to the centre |
| **NOCICEPTIVE** | Pain | analgesics, BoNT's direct antinociceptive action | see A13(4) |
| **IMMUNE ESCAPE** | Nab gate | low-protein products, serotype switching, interval extension | see A7 · A8 |

---

## Two clocks — the drug is long gone before the effect reaches its peak (A4)

Sternocleidomastoid, 50 U:

| Species | Peak | Post-peak half-life |
|---|---|---|
| `A` free toxin | 11.58 U (day 0.25) | **0.25 d** |
| `B` internalising | 25.77 U (day 0.50) | 1.50 d |
| `C` **active light chain** | 29.43 U (day 5.75) | **30.0 d**  ← calibrated parameter |
| `S` intact SNAP-25 | minimum 0.0338 (day 6.0) | resynthesis half-life 5.0 d |

**On day 2 only 0.004 % of the injected free toxin is left, and the clinical
effect has not even reached its peak yet.** Duration of action is **not** a
pharmacokinetic property of the injectate. Four sequential and progressively
slower processes set it: free-toxin residence (0.25 d) → endosomal translocation
(1.2 d) → **active light-chain persistence (28 d, calibrated)** → SNAP-25
resynthesis (5.0 d).

### The safety factor makes E(S) a threshold function

| S | r(S) release capacity | E(S) transmission efficacy |
|---|---|---|
| 1.00 | 0.992 | 1.000 |
| 0.50 | 0.940 | 0.998 |
| 0.35 | 0.771 | 0.985 |
| 0.20 | 0.500 | 0.889 |
| 0.10 | 0.111 | 0.395 |
| 0.05 | 0.015 | 0.059 |

**Destroying half the substrate produces no measurable weakness at all.** The
clinical effect only begins after S falls below about 0.20, because a normal
terminal releases roughly three times what is required (the neuromuscular junction
safety factor). This threshold is the cause of the two results below.

---

## Dose → duration is logarithmic, dose → diffusion is linear (A3)

| Dose (U) | S_min | Nadir ΔTWSTRS | Duration (d) | Pharyngeal light-chain peak | Swallowing burden (deficit·d) | Duration/burden |
|---|---|---|---|---|---|---|
| 60 | 0.123 | −3.61 | 0 | 0.70 | 0.0 | — |
| 120 | 0.065 | −7.59 | 67 | 1.41 | 0.0 | — |
| **240** | 0.034 | **−10.88** | **109** | 2.87 | 3.7 | 29.2 |
| 480 | 0.018 | −13.00 | 150 | 5.98 | 20.4 | 7.4 |
| 960 | 0.009 | −14.33 | 190 | 12.95 | 42.4 | 4.5 |
| 1920 | 0.005 | −15.19 | 225 | 30.22 | 66.8 | 3.4 |
| 3840 | 0.004 | −15.70 | 250 | 73.25 | 92.3 | 2.7 |

```
    duration          vs ln(dose) :  53.7 d per e-fold = 37.2 d per doubling   (R² = 0.993)
    pharyngeal LC     vs dose     :  4.58 U per 240 U                          (R² = 0.991)
```

**A logarithm cannot outrun a straight line.** And the pharyngeal light chain is in
fact **supralinear**: at 16 times the dose the pharyngeal load is not 16-fold but
**25.5-fold**. Once terminal binding in the injected muscle saturates, the
remaining toxin is free to diffuse. **The harm axis steepens at exactly the point
where the efficacy axis flattens.**

> Note: `P(dysphagia)` itself saturates at about 21 %. Probability is the wrong
> currency for this comparison — a bounded scale hides an unbounded exposure. Had
> only the probability been reported, the therapeutic index would have appeared to
> *recover* at high dose. It does not recover.

### Dilution is a nearly free safety lever (A6)

240 U fixed, volume alone varied:

| Dilution | Volume | Nadir ΔTWSTRS | Duration (d) | Swallowing deficit | P(dysphagia) |
|---|---|---|---|---|---|
| 200 U/mL | 1.2 mL | −11.33 | 116 | 0.055 | **4.5 %** |
| 50 U/mL (reference) | 4.8 mL | −10.88 | 109 | 0.581 | 11.0 % |
| 12.5 U/mL | 19.2 mL | −9.26 | 87 | 0.958 | **19.8 %** |

Across a 16-fold dilution range efficacy is almost flat while dysphagia risk is
not. In this model **concentration is the cheapest safety lever available at the
bedside**, and unlike dose reduction it costs almost nothing in efficacy.

---

## Sprouting is **not** the reason the effect wears off (A5)

The received explanation for loss of effect is nerve-terminal sprouting. Sprouting
itself is real, its timing is roughly right, and this model reproduces it. What is
being tested is the **causal claim**.

Sprouts grow from the same axon and fill vesicles from the **same cytosolic SNARE
pool** the light chain has been destroying. If so, sprout release must also be
gated by S. If so, sprouting cannot restore transmission while the toxin is still
working, and by the time enough toxin has gone for sprouts to function, the parent
terminal functions too. Then sprouting predicts almost nothing about duration.

| k_sp | Nadir ΔTWSTRS | Duration (d) | Benefit-time (pt·d/168 d) |
|---|---|---|---|
| **0×** (complete block) | −11.02 | **112** | 1104 |
| 1× (calibrated) | −10.88 | **109** | 1079 |
| 4× | −10.82 | 107 | 1068 |

**Abolishing sprouting entirely lengthens duration by 3 days and increases
benefit-time by 2.4 %.** For comparison, quadrupling the dose increases
benefit-time by **66.2 %**. Sprouting's share of the recovery rises to a maximum
of **9.3 %** on day 42 and then falls again; it is never dominant.

**The clock that ends it is the convolution of the light chain's own decay (28 d)
with SNAP-25 resynthesis (5.0 d).** A12 arrives at the same conclusion by an
entirely different route, sensitivity analysis: the elasticity of duration is
−0.84 for k_LC while **k_sp · k_rg is ±0.04, the smallest of any kinetic parameter
in the model.** A5 gets there by removing sprouting, A12 by perturbing it.

> **A falsifiable prediction.** In an agent that blocks sprouting — or in a patient
> group whose sprouting is impaired — the difference in duration should be
> **days, not weeks**. If a real anti-sprouting intervention greatly extends
> duration, then the SNARE-gating assumption above is wrong and sprouts are filling
> vesicles from a pool the toxin cannot reach.

---

## The central ratchet — why the trough falls cycle by cycle (A9)

BoNT blocks not only extrafusal terminals but **intrafusal (γ) terminals as well**.
It therefore lowers the abnormal spindle afferent drive, and that feeds the
plastic state `D_cen` (half-life 87 days). **The peripheral paralysis is cyclic;
the central effect is a ratchet.**

| Cycle | D_cen (end) | Trough TWSTRS | Nadir TWSTRS | S_min (SCM) |
|---|---|---|---|---|
| 1 | 0.918 | 42.96 * | 32.06 | 0.0338 |
| 2 | 0.874 | 36.44 | 29.04 | 0.0300 |
| 4 | 0.841 | 33.08 | 27.05 | 0.0295 |
| 7 | 0.830 | 32.26 | 26.45 | 0.0295 |
| 10 | 0.829 | **32.15** | 26.38 | 0.0295 |

\* The "trough" of cycle 1 is the pre-treatment baseline (t=0).

**The trough TWSTRS falls by −10.81 points over 10 cycles while the maximum
SNAP-25 cleavage does not change at all. The improvement is not in the muscle.**
This is the model's account of what clinicians call 'cumulative benefit', and it
yields a sharp prediction: **the trough should keep improving even through skipped
cycles, and should only relapse slowly after discontinuation.**

Stopping after 6 cycles (last injection day 420):

| day | TWSTRS | D_cen |
|---|---|---|
| 504 | 32.26 | 0.832 |
| 672 | 39.77 | 0.911 |
| 840 | 42.06 | 0.975 |
| 1200 | 42.88 | 0.999 |

It relapses towards the baseline of 42.94 but **never reaches it**, because D_cen
unwinds with an 87-day half-life.

---

## The 12-week rule is not explained by antibody risk (A7, A8)

Conventional clinical practice is **"do not inject more often than every 12 weeks —
short intervals induce neutralising antibodies, and potency lost to antibodies does
not come back (an absorbing state)"**. That is a quantitative claim, so it can be
tested with an immunogenicity submodel fitted to the **observed** antibody
incidence.

Searching for the interval that maximises cumulative 5-year benefit-time (all
products given as A-equivalent 240 U):

| Interval | inco benefit (pt·d) | ona benefit | rima-B benefit | ona residual potency | inco swallowing burden |
|---|---|---|---|---|---|
| **6 weeks (42 d)** | **32223** | **31819** | **29539** | 82.1 % | 567.4 |
| 8 weeks (56 d) | 30535 | 30093 | 25977 | 89.2 % | 305.0 |
| 12 weeks (84 d) | 25460 | 25131 | 19823 | 95.1 % | 129.5 |
| 26 weeks (182 d) | 14962 | 14906 | 10796 | 98.6 % | 39.1 |

**There is no interior optimum.** On the antibody axis alone, for every product
**shorter is simply better** — because at the observed antibody rates the median
patient never loses meaningful potency. The stated rationale for the 12-week
practice does not survive its own arithmetic.

An interior optimum appears **only far out in the tail** of the antibody-propensity
distribution:

| k_b multiple | Population percentile | Optimal interval |
|---|---|---|
| 1× | 50.0 th | 42 d |
| 3× | 88.9 th | 42 d |
| 10× | 99.5 th | **140 d** |
| 30× | 99.99 th | **240 d** |

So the 12-week rule is **optimal only for patients at or beyond roughly the 99th
percentile**. A rule calibrated to the tail is being applied to everyone. That may
still be the right policy — an absorbing loss is worth insuring against, and an
individual's propensity cannot be measured in advance — but it is **a decision
about the variance, not about the average patient**, and it should be stated that
way.

**What actually constrains the interval** is the last column. Halving the interval
doubles the number of injections (2.0×, 21 → 43) but multiplies cumulative
swallowing exposure by **4.4** (129.5 → 567.4) while benefit rises only +27 %. The
harm axis is **supralinear in injection frequency too**, just as it was in dose
(A3). This is where the 12-week rule belongs. Not with antibodies.

### Secondary non-response and serotype structure (A8)

Virtual cohort n = 20 000 (log-normal variate σ = 0.9 on k_b), 5 years:

| Product | Interval | 5-year secondary non-response | Median Nab | 99th-percentile Nab |
|---|---|---|---|---|
| incobotulinumtoxinA | 12 weeks | **0.00 %** | 0.0145 | 0.126 |
| onabotulinumtoxinA | 12 weeks | 5.10 % | 0.158 | 2.94 |
| onabotulinumtoxinA | 8 weeks | **12.88 %** | 0.269 | 9.27 |

Reported values are about 1–3 % for current onabotulinumtoxinA and about 0–1.1 %
for incobotulinumtoxinA, so **the model over-predicts the antibody rate of modern
formulations by roughly a factor of two** (A13(5)). The direction and the ordering —
that an 11-fold difference in protein load essentially abolishes the event — are
right, but the magnitude is not.

**Serotype structure.** BoNT-B cleaves VAMP rather than SNAP-25 and shares no
neutralising epitopes with A. Following a high responder (99.7 th percentile) over
28 cycles:

| Cycle | Nab_A | A residual potency | ΔTWSTRS if A continued | If switched to B at cycle 12 |
|---|---|---|---|---|
| 4 | 1.13 | 45.0 % | −12.93 | −12.93 |
| 12 | 5.16 | 6.7 % | −3.78 | −3.78 |
| **14** | 6.43 | 4.8 % | −3.05 | **−13.51** |
| 22 | 12.81 | 1.7 % | −1.41 | −7.06 |
| 26 | 16.91 | 1.1 % | −0.98 | −4.79 |

Cumulative 6.4-year benefit 9667 → **16270 pt·d (+68.3 %)**. But after the switch
Nab_B climbs to 3.53 and by cycle 26 it is back down to −4.79 — **the B pool
immunises independently, so the rescue is finite.** B also comes with an autonomic
cost: at equipotent dose, the salivary compartment deficit is 0.794 versus 0.100
(probability of dry mouth 15.5 % versus 4.8 %).

---

## Operator comparison (A10)

The same standard 240 U plan for 3 cycles plus one intervention. Values at cycle 3.

| Intervention | Operator | Nadir TWSTRS | Trough | ΔTWSTRS (nadir) |
|---|---|---|---|---|
| BoNT-A 240 U alone | CHEMODENERVATION | 27.70 | 34.05 | −15.24 |
| + baclofen 60 mg/d | CENTRAL DRIVE | 26.43 | 32.79 | −16.51 |
| + clonazepam 2 mg/d | GATING | 26.18 | 32.58 | −16.76 |
| + trihexyphenidyl 20 mg/d | CENTRAL DRIVE | 24.56 | 30.92 | −18.38 |
| + selective peripheral denervation 40 % | **CEILING (φ)** | 21.82 | 28.73 | **−21.12** |
| + GPi-DBS (from day 0) | **CENTRAL DRIVE** | 19.15 | 26.37 | **−23.79** |

The two strongest operators are the two that are not oral drugs: **raising φ and
lowering central drive.** No licensed oral drug in this model reaches either of
them at a useful potency — the same structural gap the influenza model found for
infected-cell killing.

---

## Calibration (A1) — seven parameters, seven reported anchors

| Parameter | Value | What it was fitted to |
|---|---|---|
| `k_cl` | 0.13613 /d | ↘ 4-week ΔTWSTRS = −10.5 at 240 U |
| `k_LC` | 0.02465 /d (t½ **28.1 d**) | ↘ 12-week ΔTWSTRS = −6.5 at 240 U |
| **`ρ`** | **0.5463** | ↘ 480 U nadir = −13.0 (**the dose-response plateau**) |
| `dys_d50`, `dys_k` | 1.716, 0.543 | dysphagia 11 % at 240 U / 5.5 % at 120 U |
| `k_b` | 8.2 × 10⁻⁵ | onaBoNT-A q12wk 5-year neutralising antibodies ≈ 2 % |
| `nw_d50` | 0.6536 | neck weakness ≈ 9 % at 240 U |

A note on `k_LC`: independent estimates of intraneuronal BoNT/A light-chain
persistence span weeks to months. The fitted 28.1 days falls inside that range,
but **it was not constrained to do so.**

### Predicted time course

| Dose | wk 1 | wk 2 | **[wk 4]** | wk 6 | wk 8 | **[wk 12]** | wk 16 | Nadir | Nadir day |
|---|---|---|---|---|---|---|---|---|---|
| Placebo | +0.00 | −0.00 | −0.01 | −0.01 | −0.01 | −0.01 | −0.01 | −0.01 | — |
| 60 U | −1.63 | −3.33 | −3.48 | −2.97 | −2.58 | −2.29 | −2.24 | −3.61 | 21 d |
| 120 U | −3.86 | −6.57 | −7.58 | −6.73 | −5.39 | −3.58 | −3.00 | −7.59 | 26 d |
| **240 U** | −5.28 | −8.29 | **−10.50** | −10.83 | −9.92 | **−6.50** | −4.34 | −10.88 | 38 d |
| 480 U | −5.86 | −8.80 | −11.45 | −12.67 | −12.99 | −11.04 | −7.20 | −13.00 | 54 d |

Only the two bracketed cells and the 480 U nadir were used in the fit. Everything
else is prediction.

**Held out of the fit and used for validation:** (a) 60 and 120 U should respond
shallowly and sub-proportionally → **met**; (b) above 480 U dysphagia should keep
rising after efficacy has stopped → **met**; (c) the peak effect should fall
between weeks 2 and 4 → **not met** (38 days); (d) the duration the patient
perceives should be 10–12 weeks → **not met** (109 days). (c) and (d) are reported
in A13(2).

---

## The uncomfortable top of the sensitivity analysis (A12)

When each parameter is perturbed by +25 %, the two that move the **duration**
read-out most are `L50` (elasticity +5.54) and `kp_off` (+3.60). Neither is toxin
biology — they are the position of the clinical severity curve and the pain time
constant.

They dominate because "duration" here is defined as **the time at which a fixed
threshold (MCID) is re-crossed**. Anything that shifts the whole curve up or down
moves that crossing point a great deal.

> **A substantial part of the 'duration of effect' that trials report is therefore
> not a property of how long the drug acts, but a property of where the threshold
> sits relative to the response.**
> This is the same measurement problem A13(2) runs into from the other side.

Among mechanistic parameters alone the ordering is clean and consistent with A4 and
A5: `k_LC` (−0.84) → `S50` · `k_cl` · `k_syn` (±0.5) → `k_sp` · `k_rg` (±0.04). The
owners of the three read-outs (maximum effect, duration, diffusion risk) are almost
disjoint, which is why one dose knob cannot trade between them (A3).

---

## Five discrepancies — where this model does not fit (A13)

### (1) The **factorisation** of φ = ρ·Σw **is not identified** — the most important limitation

The plateau identifies the **product φ**. It does **not** identify ρ and Σw
separately. The split between "placement accuracy (+15.66)" and "target list
(+6.93)" in the table above rests on the assumed muscle torque weights `w_m`, and
**that is the weakest link in the whole construction.**

What is identified: that **the total non-dose headroom (22.59 points) is 5.2 times
the dose headroom (4.37 points)**. That conclusion does not depend on the split.
What does depend on the split is the advice about which of the two improvements
those 22.59 points should be invested in.

### (2) The effect-time curve is too square

The model reaches 10 % of maximum effect at 1.5 days, 50 % at 8.5 days, 90 % at 22
days, and stays above the MCID (4.5 points) for **109 days** (15.6 weeks). Reported
values are a median time to first effect of about 7 days and a patient-perceived
loss at 10–12 weeks. **So the model rises too fast and falls too slowly** — while
sitting exactly on the two fitted TWSTRS time points.

These are not two errors but one: **fitting a rating scale at two visit time points
does not reproduce the onset and offset the patient perceives.** Perceived benefit
therefore cannot be an affine function of the TWSTRS total. Forcing both would
require a separate perceptual read-out with its own threshold and hysteresis, and
there is none here.

### (3) One Hill function cannot capture both ends of the dose-response

The model reproduces 240 U > 120 U with 480 U only marginally better than 240 U.
What it **cannot** reproduce with a single `L50` is the near-linearity over 60–180 U
that some series report. That is **evidence about the *shape* of `sev_ss(L)`, not
evidence about potency**, and a model fitted with `k_cl` alone would have quietly
absorbed the shape error into a toxin parameter.

### (4) Pain improves more than the torque reduction explains

The model's nadir pain change is about −2.0/20, whereas the reported TWSTRS-pain
change at 240 U is −2.5 ~ −3.5; closing that gap requires a **direct
antinociceptive term** — that is, BoNT must be doing more than releasing the
muscle; it must be **blocking CGRP and substance P release** at nociceptive
terminals. This is a falsifiable claim rather than a fit, which is why
`pain_direct` is left at a default of 0 so that the gap remains visible.

### (5) One `k_b` cannot hold both formulation eras

The calibration lands the median patient's Nab exactly on target, but because the
memory B-cell recall term (1 + β·B) is non-linear, the cohort's tail becomes
heavier than the log-normal assumption implies. As a result the n = 20 000 cohort
gives a secondary non-response rate of **5.10 %** for current onabotulinumtoxinA
q12wk, whereas the reported value is about 1–3 %. **The model over-predicts modern
formulations by roughly a factor of two.**

The ordering is right — an 11-fold difference in protein load essentially abolishes
the event for incobotulinumtoxinA (0.00 %), and shortening the interval to 8 weeks
raises it to 12.88 %. Beyond that it is wrong in two directions.

| | Model | Reported |
|---|---|---|
| Modern 5 ng/100 U | 5.10 % | 1–3 % |
| Historical 25 ng/100 U | 39.46 % | ~9.5 % |
| Ratio for a 5-fold increase in load | **7.7×** | about 4× |

**(a) The absolute value is too high** — the median patient's Nab is structurally
on target, but the recall term (1 + β·B) makes the cohort tail heavier than a
log-normal `k_b` alone would predict.
**(b) It is too steep in load** — a 5-fold protein load raises the event rate
7.7-fold, whereas the historical comparison gives about 4-fold.

(b) matters more: **a pure antigen-dose model with a single saturation constant
cannot be this steep and hit the absolute value at the same time.** Either antigen
processing saturates more strongly than `Ka_ag` allows, or individual variability
lies on several axes (HLA type, prior exposure) rather than on `k_b` alone. Neither
was put in, so this line is a gap rather than a fit, and **it is the only one of
A1's seven anchors that the model misses badly**.

### What is not in the model

**Primary non-response from mis-selected targets**, as distinct from a low φ;
treating needle EMG accuracy as a distribution rather than a binary; phasic versus
tonic dystonia; head tremor; the physiology of the sensory trick (geste
antagoniste); and the anticholinergic cognitive burden that in reality constrains
trihexyphenidyl far more than efficacy does. The abo-/rima- unit conversions are
clinical rules of thumb and **must not be read as claims of potency equivalence.**

---

## State vector (70 ODEs)

| Block | States | Count |
|---|---|---|
| Muscles m = 1…8 | `A_m` free toxin · `B_m` internalised · `C_m` active light chain · `S_m` intact SNARE · `Q_m` sprout capacity | 40 |
| Swallowing compartment (pharyngeal constrictors) | `A_sw`, `B_sw`, `C_sw`, `S_sw` | 4 |
| Autonomic compartment (salivary glands · ganglia) | `A_au`, `B_au`, `C_au`, `S_au` | 4 |
| Humoral immunity, serotypes A / B | `Ag`, `Bmem`, `Nab` × 2 | 6 |
| Central · spinal | `D_cen`, `RecInh`, `SurrInh`, `Cbll` | 4 |
| Clinical | `Sev`, `Pain`, `Disab` | 3 |
| Oral adjunct PK | trihexyphenidyl · baclofen · clonazepam (depot + central) | 6 |
| Accounting | `AUCben`, `AUCdys`, `CumU` | 3 |
| | | **70** |

### The muscle table — where the factorisation of φ lives

| m | Muscle | w_m | Relative mass | Pharyngeal proximity | Posterior | Standard | Extended |
|---|---|---|---|---|---|---|---|
| 1 | sternocleidomastoid (contralateral) | 0.20 | 1.00 | 0.30 | | ✓ | ✓ |
| 2 | splenius capitis (ipsilateral) | 0.24 | 1.60 | 0.08 | ✓ | ✓ | ✓ |
| 3 | trapezius | 0.13 | 2.20 | 0.04 | ✓ | ✓ | ✓ |
| 4 | levator scapulae | 0.08 | 0.80 | 0.05 | ✓ | ✓ | ✓ |
| 5 | semispinalis (**deep**) | 0.12 | 1.80 | 0.12 | ✓ | | ✓ |
| 6 | scalene group | 0.06 | 0.90 | 0.28 | | | |
| 7 | obliquus capitis inferior (**deep**) | 0.10 | 0.50 | 0.18 | ✓ | | ✓ |
| 8 | longus colli (effectively unreachable) | 0.07 | 0.70 | 0.60 | | | |

---

## Products

| Product | Serotype | Unit conversion | Protein load (ng/100 label-U) | Light-chain decay |
|---|---|---|---|---|
| incobotulinumtoxinA | A | 1.00 | **0.44** | 1.0× |
| onabotulinumtoxinA | A | 1.00 | **5.00** | 1.0× |
| abobotulinumtoxinA | A | 0.34 | 0.87 | 1.0× |
| daxibotulinumtoxinA | A | 1.00 | 0.50 | **0.50×** (phenomenological) |
| rimabotulinumtoxinB | B | 0.03 | 0.10 | 1.45× |

The antigen is **the protein load, not the units**: at the same nominal dose, ona
versus inco differs 11-fold. Serotypes A and B share no neutralising epitopes, so
switching to B revives potency — but **only finitely**, because the B pool
immunises independently.

### A long-acting toxin does not move the ceiling (A11)

The three conventional serotype-A products are **identical** after unit conversion
— in this model they differ only in unit scale and protein load, and protein load
acts on the antibody pool, not on the muscle.

daxibotulinumtoxinA is **deeper and longer**. In this model depth and duration are
not independent, because both are set by light-chain persistence. What matters is
the **exchange rate**:

| | incobotulinumtoxinA | daxibotulinumtoxinA | Change |
|---|---|---|---|
| Nadir ΔTWSTRS | −10.88 | −12.72 | **+16.9 %** |
| Duration | 109 d | 263 d | **+141 %** |

It buys **8.4 times** more duration than depth. This is the same asymmetry as the
reported values — about 24 weeks versus about 12 weeks in cervical dystonia, with
broadly similar maximum TWSTRS change — and the model's ratio (2.4×) comes close to
the reported ratio (2×) **without having been fitted to it.**

**But the ceiling did not move.** Both products still fall short of the asymptotic
ceiling at φ = 0.355 (−13.24). A long-acting toxin is a **convenience** improvement
— the same control in fewer visits — and only a marginal efficacy improvement. It
buys benefit **per injection**, not the benefit **per cycle** that A2's geometry
levers buy.

---

## Reproduce

```bash
# mechanistic map
dot -Tsvg cdys_qsp_model.dot -o cdys_qsp_model.svg
dot -Tpng -Gdpi=150 cdys_qsp_model.dot -o cdys_qsp_model.png

# every number in this README (without R)
python3 cdys_reference_check.py > cdys_reference_output.txt

# re-verify the references (actually queried against PubMed)
python3 cdys_resolve_refs.py
```

```r
# the mrgsolve port
library(mrgsolve)
mod <- mread("cdys_mrgsolve_model.R")
env <- mrgsolve::env_get(mod)
out <- env$sim_scenario(mod, "S03_std_240U_q12wk")
plot(out, TWSTRS_TOTAL ~ time)
env$run_all_scenarios(mod)         # summary of the 15 scenarios

# Shiny dashboard (tab 2 is the point)
shiny::runApp("cdys_shiny_app.R")
```

---

## ⚠️ Disclaimer

This is a **qualitative to semi-quantitative QSP model for education and research**.
It was assembled from the public literature and clinical trial data but has not
been independently validated or certified, and **must not be used directly for real
clinical decision-making, prescribing, or regulatory submission.** In particular,
the value ρ = 0.55 above is a *model-based inference* about the shape of published
dose-response curves, not a measurement that any individual patient's injection
reaches 55 % of the muscle.
