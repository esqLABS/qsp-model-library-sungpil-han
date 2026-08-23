# Generalized Anxiety Disorder — QSP Model

> **One corticolimbic gain, four factors, four clocks — and a fifth clock inside the ruler.**

<p align="center">
  <a href="../../../generalized-anxiety-disorder/gad_qsp_model.svg"><img src="../../../generalized-anxiety-disorder/gad_qsp_model.png" width="820" alt="GAD QSP mechanistic map"></a>
</p>

---

## 0. The claim in one line

This model writes anxiety as **a single dimensionless corticolimbic gain**

```
        E_amy × S_glu
Φ  =  ───────────────────────────
      C_pfc·(1 + e_ex·EXPECT) × I_gaba
```

It is a **ratio** with two factors in the numerator and two in the denominator, and
each of the four factors is moved by a different drug class **with a different time
constant**. And all of this is read through a **saturating scale**, HAM-A, which
lays its own **fifth clock** (expectancy effect + regression to the mean) on top.

| Factor | Meaning | The drugs that move it | Time constant (as computed by the model) |
|---|---|---|---|
| `I_gaba` | GABA-A inhibitory efficiency (denominator) | benzodiazepines | **0.5 d** (effect t50 2.5 d) |
| `S_glu` | presynaptic glutamate release probability (numerator) | pregabalin (α2δ-1) | **2.0 d** (effect t50 4.5 d) |
| `C_pfc` | prefrontal top-down regulatory capacity (denominator) | SSRIs/SNRIs, CBT | **40 d** (effect t50 35 d) |
| `E_amy` | amygdala/BNST excitatory drive (numerator) | CBT, chronic 5-HT | a few weeks |
| `EXPECT`+`FLUCT` | expectancy effect + regression to the mean at enrolment | (common to every arm) | **the fastest** — 57% of the 8-week placebo effect by week 1 |

Three results come out of this algebra **by calculation and not by assumption** —
and three priors held in advance are **refuted and reported as they stand**.

---

## 1. The files in this directory

| File | Contents |
|---|---|
| [`gad_qsp_model.dot`](../../../generalized-anxiety-disorder/gad_qsp_model.dot) · [`.svg`](../../../generalized-anxiety-disorder/gad_qsp_model.svg) · [`.png`](../../../generalized-anxiety-disorder/gad_qsp_model.png) | Mechanistic map — **206 nodes · 19 clusters · 301 edges** |
| [`gad_python_reference.py`](../../../generalized-anxiety-disorder/gad_python_reference.py) | **The runnable reference implementation** (49 ODEs). Every number in the README, in mrgsolve and in Shiny comes from here |
| [`gad_calibrate.py`](../../../generalized-anxiety-disorder/gad_calibrate.py) | Calibration. **Five parameters, six numbers.** Records even which attempt failed and why |
| [`gad_analysis.py`](../../../generalized-anxiety-disorder/gad_analysis.py) → [`gad_reference_output.txt`](../../../generalized-anxiety-disorder/gad_reference_output.txt) | The complete run of all 12 scenario blocks (472 lines) |
| [`gad_mrgsolve_model.R`](../../../generalized-anxiety-disorder/gad_mrgsolve_model.R) | mrgsolve model (50 compartments). The reference implementation carried across equation by equation |
| [`mkmrgsolve.py`](../../../generalized-anxiety-disorder/mkmrgsolve.py) | The generator for the R file above — makes it impossible for the calibrated values and the initial conditions to diverge |
| [`gad_shiny_app.R`](../../../generalized-anxiety-disorder/gad_shiny_app.R) | Shiny dashboard (**12 tabs**) |
| [`gad_references.md`](../../../generalized-anxiety-disorder/gad_references.md) · [`mkrefs.py`](../../../generalized-anxiety-disorder/mkrefs.py) | **118 references / 113 unique PMIDs**, all looked up live on PubMed |
| `gad_calibration.json` · `gad_attractor_cache.npz` | The calibration results and the drug-free attractor cache |
| `calib.log` · `calib_stage1_joint.log` | The calibration logs (including the failed first attempt) |

Rendering and running:

```bash
dot -Tsvg gad_qsp_model.dot -o gad_qsp_model.svg
dot -Tpng -Gdpi=150 gad_qsp_model.dot -o gad_qsp_model.png
python3 gad_calibrate.py      # re-calibrate the five parameters
python3 mkmrgsolve.py         # regenerate gad_mrgsolve_model.R
python3 gad_analysis.py       # regenerate gad_reference_output.txt
python3 mkrefs.py             # regenerate gad_references.md
Rscript -e 'shiny::runApp("gad_shiny_app.R")'
```

---

## 2. What was fitted and what was predicted

**Only five parameters were fitted.**

| Parameter | Fitted to |
|---|---|
| `fluct0` = 3.85, `kfl` = 2.00 /d | The **shape** of the placebo-arm curve — Khan 2011 (PMID 21694613) week 1 −5.94 / week 8 −11.10, Rickels 2005 (PMID 16143734) week 4 −8.40 |
| `k5ht_eff` = 0.126 | The escitalopram–placebo difference **−2.45** (Slee 2019 NMA, PMID 30712879) |
| `emax_pgb` = 0.094 | The pregabalin–placebo difference **−2.79** (Slee 2019 NMA) |
| `emax_a2` = 0.571 | The alprazolam 1.5 mg–placebo difference **−2.50** (Rickels 2005) |

`dvisit` = 0.30 and `e_ex` = 0.55 were **not fitted but fixed structurally**. The
reason is in §6.

**The result of the fit (six numbers):**

```
  placebo week 1        :   -5.75   target   -5.94
  placebo week 4        :   -9.60   target   -8.40      <-- miss
  placebo week 8        :  -10.09   target  -11.10      <-- miss
  escitalopram 10 delta :   -2.45   target   -2.45
  pregabalin 300 delta  :   -2.79   target   -2.79
  alprazolam 1.5 delta  :   -2.50   target   -2.50
```

The placebo-arm residual is structural. The placebo curves of the two trials are
almost **straight**, at −5.94 → −8.40 → −11.10, whereas this model's placebo
response **plateaus** around week 4. Put honestly, this model cannot produce a
placebo curve that keeps falling all the way to week 8.

**Everything else is out-of-sample.** That is true of every number in §3–§5 below.

---

## 3. Out-of-sample scorecard

### 3.1 Slee 2019 Lancet network meta-analysis (PMID 30712879) — placebo-subtracted difference at 8 weeks

| Drug | Model | Reported | |
|---|---|---|---|
| escitalopram 10 mg | −2.45 | −2.45 | ← fitted |
| pregabalin 300 mg | −2.79 | −2.79 | ← fitted |
| venlafaxine ER 150 mg | **−4.24** | −2.69 | over-predicted |
| duloxetine 60 mg | **−2.53** | −3.13 | under-predicted |
| quetiapine XR 150 mg | **−0.57** | −3.60 | **a big miss** |
| buspirone 45 mg | **−0.27** | (effective, but the sample is small) | under-predicted |

That the two SNRIs are wrong in opposite directions is meaningful — the NET
occupancy curve (EC50 260 ng/mL vs 46 ng/mL) is too generous to venlafaxine and too
stingy to duloxetine.

### 3.2 Khan 2011 quetiapine XR (PMID 21694613) — this model's largest failure

| Arm | Week 1 model | Week 1 reported | Week 8 model | Week 8 reported |
|---|---|---|---|---|
| Placebo | −5.75 | −5.94 | −10.09 | −11.10 |
| Quetiapine 50 | −6.24 | −7.47 | −10.61 | −13.31 |
| Quetiapine 150 | −6.18 | −8.19 | −10.66 | −13.54 |
| Quetiapine 300 | −6.10 | −7.23 | −10.68 | −11.87 |

The model gives quetiapine almost no effect. In this model quetiapine's only
anxiolytic route is NET inhibition by norquetiapine (occupancy 0.295 at 150 mg),
whereas the real drug works through 5-HT2A/2C antagonism, α1 blockade and H1
blockade, and this model **does not have** those routes. Nor does it reproduce the
reported inverted-U dose-response (300 mg < 150 mg) — even adding
adverse-event-driven dropout in the virtual population of §9, it stays monotonic at
−10.16 (300) vs −9.88 (150). This is not a parameter problem but a **structural
deficit**, and it is reported as it stands.

### 3.3 Allgulander 2006 randomised discontinuation trial (PMID 16316482) — **the decisive test of the model**

Simulating the design exactly as it was written — 12 weeks open-label →
randomisation of the responders → 76 weeks of follow-up — with nothing but the
acute-phase calibration:

| | Model | Reported |
|---|---|---|
| Open-label responder fraction | 65% | 76% |
| Relapse rate at 76 weeks, escitalopram continued | **17%** | **19%** |
| Relapse rate at 76 weeks, switched to placebo | **66%** | **56%** |
| Hazard ratio (placebo/drug) | **3.83** | **4.04** |

The acute-phase calibration never saw any of these numbers. Stop the drug and
`C_pfc` falls back down along its own time constant (30 d), and that is what pulls
the two arms apart — Φ_n at 76 weeks 1.14 (continued) against 1.41 (placebo).

### 3.4 Psychic/somatic subscale dissociation (Rickels 2005)

Rickels reported it as a pattern of significance: pregabalin significant on
**both** psychic and somatic, alprazolam significant on psychic only (somatic
p = 0.21). The model reproduces that **direction** — the somatic/psychic ratio of
the placebo-subtracted difference is 1.33 for pregabalin and 0.89 for alprazolam.
Because `S_glu` feeds the muscle-tension item directly, whereas `I_gaba` reaches the
somatic items only through Φ.

---

## 4. Three priors the model refuted

By the rule of this repository, when an expectation held before the model was built
collapses under the calculation, it is not quietly deleted but written down as it
stands.

**(1) "The SSRI dose-response is flat because of the SERT occupancy hyperbola" —
no.** The occupancy certainly does flatten (10 mg 0.744 → 20 mg 0.853, +15%
relative). And yet the HAM-A difference goes −2.45 → −4.23, **73% larger**. Because
as the reuptake term approaches the floor, extracellular 5-HT rises steeply
(1.98 → 2.62). Baldwin 2006 (PMID 16946363) reported that 20 mg is not
significantly better than 10 mg, whereas the model predicts an advantage of 1.8
points. The clinical flatness has to come from something that is not in this model
(most plausibly dropout at 20 mg).

**(2) "Combination is sub-additive because of scale saturation" — the opposite.**
Φ is **exactly multiplicative** in all three pairs (errors −0.4%, −1.2%, −0.7%).
That is the strong result of this section. On HAM-A, though, it is 107%, 115%,
112% — **super-additive**. The HAM-A link function is **convex** while the effect is
small and only turns concave as it approaches the floor, and the 3–4 point effects
of GAD trials sit in the convex region. So the clinical prediction is inverted:
**combination therapy in GAD ought to exceed the sum of the monotherapy
differences.**

**(3) "Because Φ is a ratio, the drug–placebo difference grows with severity" — it
cancels out.** Across mild/moderate/severe tertiles the difference is
−1.68 / −1.76 / −1.46: **flat and non-monotonic**. Severe patients sit on a steeper
part of the scale, but **their placebo response is larger to the same degree**
(−12.34 against −8.21). The argument that power is gained by selecting severe
patients does not follow in this model.

---

## 5. What survives

**(a) The occupancy-effect gap comes out of the calculation.** SERT is half blocked
in **0.5 days** and HAM-A moves half way in **35 days** — **70-fold**. This is not
fitted; it is two slow stages connected in series:

| Stage | t50 |
|---|---|
| SERT occupancy | 0.5 d |
| 5-HT1A autoreceptor desensitisation (the gate) | 10.5 d |
| Extracellular 5-HT | 14.5 d |
| BDNF | 23.5 d |
| C_pfc | 40.0 d |
| HAM-A | 35.0 d |

At the same occupancy (0.574 → 0.744) the acute rise in 5-HT is **12.9%** and the
chronic rise **98.6%**. Because the gate has opened. And that is why the onset curve
is not exponential but **sigmoid**.

**(b) The fifth clock does not shrink the difference. It buries it.** At week 1 the
placebo arm moves −5.75 points, and the largest drug increment laid on top of that
is −1.75 points. A clinician observing an individual patient at week 1 is **looking
at the fifth clock**, whatever was prescribed. The drug's clock is visible only in
the *difference*.

**(c) Assay sensitivity is derived.** Sweeping `e_ex` alone (not one line of the
pharmacology changes):

| e_ex | placebo ΔHAM-A | esc10 diff | pgb300 diff | lzp3 diff |
|---|---|---|---|---|
| 0.00 | −3.85 | −1.57 | −2.03 | −1.75 |
| 0.55 | −10.09 | −2.45 | −2.79 | −2.91 |
| 1.00 | −16.95 | −3.10 | −3.20 | −3.17 |
| 1.50 | −20.03 | **−0.04** | **−0.21** | **−0.09** |

At a site where the expectancy effect is very large, **the measurable effect
disappears even though the pharmacology is identical**. Because expectancy raises
`C_pfc`, `C_pfc` sits in the **denominator** of Φ, and HAM-A is a saturating
function of Φ. A common cause of death of failed multicentre trials comes out here
as a result rather than as an assumption.

**(d) Benzodiazepines: one subunit makes two clocks.**
On continuous lorazepam 3 mg/d

| | Day 1 | Day 14 | Day 168 |
|---|---|---|---|
| α1 pool `R_a1` (sedation) | 0.895 | **0.380** | 0.377 |
| α2/3 pool `R_a2` (anxiolysis) | 0.999 | 0.972 | **0.892** |

Sedation loses 42% within two weeks while the anxiolytic pool is still 0.892 at
day 168. This is the subunit-level explanation of the clinical observation that
"tolerance develops to sedation but hardly at all to the anxiolytic effect".

**(e) The benzodiazepine bridge takes back exactly what it lent.** Laying three
weeks of lorazepam on top of escitalopram earns −1.76 points at week 1, but as
`DEPEND` unwinds over the taper it converges to −0.22 by week 6 and −0.01 by week
12.

---

## 6. Two identifiability failures, met and recorded

The failed attempts have been left in these files rather than erased. See the
docstring of [`gad_calibrate.py`](../../../generalized-anxiety-disorder/gad_calibrate.py)
and [`calib_stage1_joint.log`](../../../generalized-anxiety-disorder/calib_stage1_joint.log).

1. **The first joint fit missed placebo week 1 structurally.** Because `kfl` (the
   decay rate of regression to the mean at enrolment) had been pinned by hand at
   0.075/d. The expectancy term has to travel through WORRY · SLEEPD · AUTON and
   therefore **cannot in principle be fast**; the only thing that can be fast is
   `FLUCT`, and that was precisely what was tied down. Release `kfl` and week 1 came
   right.

2. **"Real down-regulation" and "selection bias" in the placebo response cannot be
   separated by three summary numbers.** Fitting `fluct0`, `kfl` and `dvisit`
   together, the optimiser pushed `dvisit` to its lower bound (0.02) and dumped the
   entire placebo response into the additive `FLUCT` term (9.72 points, τ 5.8 d). On
   those three numbers alone it is a perfect fit, but doing so **quietly switches
   off the biological expectancy route on which the assay-sensitivity argument of
   §5(c) depends.** So `dvisit` and `e_ex` are not fitted but fixed structurally,
   and swept explicitly in §5(c) instead. **Every claim in this repository that
   depends on the size of the expectancy effect is a claim about a structural
   choice, and has been marked as such.**

`kfl` ended up pinned at its upper bound (2.0/d, τ 0.35 d). That is, the fit wants
the increment present at enrolment to disappear **effectively instantly** — meaning
that it is **baseline-score inflation** rather than a real fluctuation in symptoms.
That too is written down as it stands.

---

## 7. Defects found and fixed during construction

Every one of these would have stayed plausible-looking had the Python reference
implementation never actually been run:

1. **The HPA axis diverged.** The loop gain of CRH → cortisol → GR down-regulation →
   weakened feedback → CRH exceeded 1. Bounded by making excess cortisol enter every
   downstream effect as `x/(1+0.5x)`.
2. **The gain of the whole disease loop also exceeded 1.** Amygdala → sleep
   deprivation → amygdala diverged. Every downstream consequence of Φ was made to
   read not `phi_n` but the saturating excess term
   `z = (Φn−1)/(1+(Φn−1)/3)`, and a logistic ceiling was put on `E_amy`
   (EAMAX = 3) and a ceiling on `SLEEPD` (4).
3. **The CBT windows overwrote one another.** With more than one CBT arm in a single
   run, the later one was switching the earlier one off — caught because the "CBT
   alone" arm came out **numerically identical** to placebo. Fixed by making them
   accumulate.
4. **The drug effect was fitted to the absolute change.** The result was that
   alprazolam's total change (−10.9) was matched exactly while the
   placebo-subtracted difference was wrong by more than a factor of two. The drug
   parameters had absorbed all the error of the placebo model. Changed to fit the
   **difference**.
5. **Pregabalin had two routes to the score.** The direct symptom-relief terms
   (`b_pgb`, `h_pgb`) had been set large, so the calibration pushed `emax_pgb` down
   to 0.029 — that is, a state in which `S_glu` does nothing at all and yet the drug
   works. The claim that "S_glu is pregabalin's factor" is then empty. The direct
   route was reduced and the structural choice made explicit.
6. **The scale of disease severity was out by a third.** The *stable* HAM-A at
   DIS = 1 was 25.9, whereas the 24–27 that trials report is the score *at
   enrolment*. The enrolment criterion selects the peaks of the fluctuation. Fixed
   by dividing `kdis_*` by 2.5, giving a stable score of 19.4 and an enrolment score
   of 23.1. An error that could not have been seen without the virtual population.

---

## 8. Numerical verification

- **The integrator**: a hand-written RK4 (dt = 30 min, vectorised over subjects)
  checked against `scipy` LSODA (rtol 1e-9) — after 7 days of escitalopram 10 mg QD
  + pregabalin 300 mg BID, the relative error on HAM-A is **1.3e-06**, and
  1.1e-03 on the worst state variable (a decaying absorption compartment).
- **Step size**: halving dt (to 1/96 d) moves HAM-A by only **2.0e-05**.
- **The attractor**: the drug-free equilibrium point is identical to four decimal
  places at dt = 1/8, 1/12 and 1/24 (the equilibrium is a fixed point of the RK4
  map, so it is independent of the step).
- **The mrgsolve port**: the initial conditions in `$MAIN` are a fit of the
  reference implementation's attractor to a cubic in DIS, with a maximum relative
  error of **5.7e-03** (SLEEPD). Because the generator ([`mkmrgsolve.py`](../../../generalized-anxiety-disorder/mkmrgsolve.py))
  writes the calibrated values and the initial conditions in together, the two
  implementations cannot diverge.

**One honest limitation:** there is no R toolchain in this environment, so
[`gad_mrgsolve_model.R`](../../../generalized-anxiety-disorder/gad_mrgsolve_model.R) and
[`gad_shiny_app.R`](../../../generalized-anxiety-disorder/gad_shiny_app.R) **have not been run.** The R files were
written by checking them equation by equation against the Python reference
implementation, but they have not been verified in themselves. Also, mrgsolve has no
way to express the discrete expectancy update at a visit as a dose, so the same
mapping is approximated with a fast-decaying `VISCUE` compartment
(kvis = 48/d) — stated explicitly in the header.

---

## 9. Model structure

**49 ODEs** (50 compartments in mrgsolve, including `VISCUE`):

| Block | State variables |
|---|---|
| PK (20) | escitalopram, two compartments + effect site; venlafaxine + ODV; duloxetine; pregabalin; benzodiazepine + effect site; buspirone + 1-PP; quetiapine + norquetiapine |
| Neurotransmission (4) | `SHT`, `AUTO` (5-HT1A autoreceptor), `NE`, `A2AUTO` |
| Plasticity · circuitry (4) | `BDNF`, `CPFC`, `EAMY`, `TRAF` (α2δ transport) |
| GABA-A (3) | `RA1` (sedation, τ 2.6 d), `RA2` (anxiolysis, τ 42 d), `DEPEND` (adaptation, τ 20 d) |
| HPA (4) | `CRH`, `ACTH`, `CORT`, `GR` |
| Symptom layer (4) | `SNS`, `AUTON`, `SLEEPD`, `WORRY` |
| Trial machinery (2) | `EXPECT`, `FLUCT` |
| Adverse effects (6) | `RNAU`, `RDIZZ`, `RH1`, `SEXD`, `WT`, `RACT` |
| Comorbidity · risk (2) | `MADRSS`, `CUMHAZ` |

**Constants fixed by measurement** (not fitted):

| Constant | Value | Source |
|---|---|---|
| escitalopram SERT EC50 | 5 ng/mL (occupancy 0.80 at 10 mg) | PET/TDM (PMID 38287888, 27557550) |
| venlafaxine NET EC50 | 260 ng/mL → occupancy 0.21–0.44 at 75–225 mg | Arakawa 2019, measured 8–61% (PMID 30649319) |
| benzodiazepine binding-site EC50 | **96 ng/mL** lorazepam equivalent → occupancy **~20%** at anxiolytic doses | Atack 2007 [¹¹C]flumazenil (PMID 17164474) |
| pregabalin t½ / renal clearance | 6.3 h / proportional to CrCl | PMID 17940637, 12638396 |
| α1 = sedation, α2/3 = anxiolysis, α5 = tolerance to sedation | structural | knock-in studies (PMID 15282283) |

**Seventeen scenarios** are in [`gad_reference_output.txt`](../../../generalized-anxiety-disorder/gad_reference_output.txt):
drug-free natural history · comparison of the four clocks · Rickels 4 weeks, 5 arms ·
Khan 8 weeks, 4 arms · dose-response sweeps (escitalopram at 4 doses, venlafaxine at
3 doses × 3 CYP2D6 phenotypes, duloxetine at 2 doses, buspirone) ·
placebo/expectancy sweep · benzodiazepine after 12 weeks, abrupt discontinuation vs
4-week taper vs continuation · three combination pairs · the benzodiazepine bridge ·
CBT alone/in combination · virtual population of 6 arms (n = 200, dropout · LOCF) ·
Allgulander randomised discontinuation (n = 260 → 169×2, 76 weeks) · renal
function/CYP2D6/adherence · verification of the integrator.

---

## 10. Shiny dashboard (12 tabs)

Patient · prescription → PK/occupancy → **the four factors** → **the four clocks** →
clinical endpoints → trial simulator → scenario comparison → dose-response →
benzodiazepines (tolerance · rebound) → adverse effects → **the assay-sensitivity
sweep** → relapse prevention. Each tab is built so as to reconstruct interactively
the calculations of §3–§5 above.

---

## 11. References

[`gad_references.md`](../../../generalized-anxiety-disorder/gad_references.md) — **118 entries / 113 unique PMIDs**, in 20 sections.
Every title, journal, year, author and PMID was looked up live by [`mkrefs.py`](../../../generalized-anxiety-disorder/mkrefs.py) through the
NCBI `esearch` + `esummary`, and none of it was written from memory. Each entry also
records **what the model took from that paper (its intent)**, so that where the
search picked up a paper other than the one intended it is visible at a glance (they
are not all right — that is the purpose of this annotation).

> ⚠️ The default sort of the E-utilities is not relevance but **date**. Leave out
> `sort=relevance` and "anxiety and cardiovascular risk" returns a 2026 alopecia
> meta-analysis. That is exactly what happened on the first build, and it is left as
> a comment in `mkrefs.py` together with a record of the accident.

---

## ⚠️ Disclaimer

This is a QSP model for educational and research purposes. It was built from the
public literature and clinical-trial data but has not been independently validated
or certified, and **must not be used directly for real clinical decision-making,
prescribing, or regulatory submission.** In particular, do not quote a number from
this model without reading the limitations set out in §3.2 (quetiapine), §4 (the
refuted priors), §6 (non-identifiability) and §8 (the R code has not been run).
