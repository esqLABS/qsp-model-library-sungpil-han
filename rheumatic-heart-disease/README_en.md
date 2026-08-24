# Rheumatic Heart Disease (RHD) QSP Model

**Acute rheumatic fever → rheumatic mitral stenosis · RHD**

Describes the entire course from acute rheumatic fever (ARF) to rheumatic
mitral stenosis in a single quantitative systems pharmacology (QSP) model.
The starting point of this model is the following sentence:

> **RHD is managed by two numbers, and neither of them is the number that actually matters.**

---

## 1. What This Model Claims (The three claims)

### 1-1. Adherence is not protection — the two numbers reverse the ranking

The number secondary prevention programmes audit is **injections received /
injections scheduled**. The quantity that actually protects the patient is
**the time during which the serum penicillin G concentration exceeds
0.02 µg/mL**.

Benzathine penicillin G (BPG) is a depot showing flip-flop kinetics. The
serum concentration tracks **release** (t½ ≈ 9 days), not elimination
(t½ 30 minutes), and its amplitude is inversely proportional to body size.
So the period covered by a single 1.2 MU injection varies with body weight:

| Body weight | Days above 0.02 µg/mL | 28-day interval | 21-day interval |
|---|---|---|---|
| 40 kg | 23.4 d | Falls short | Meets |
| 55 kg | 20.6 d | Falls short | Close |
| 70 kg | 18.1 d | Falls short | Falls short |
| 100 kg | 14.8 d | Falls short | Falls short |

A conclusion follows from this. **At adult body weights, a single injection
never covers a 28-day interval, under any scenario.** And the audited number
and the protecting number can rank two patients in **opposite** order:

| Patient | Doses received | Time actually protected |
|---|---|---|
| 55 kg, 3-week interval, **80%** administered | 80% | **79.0%** |
| 95 kg, 4-week interval, **100%** administered | 100% | **60.0%** |

The patient with 100% adherence ends up less protected than the patient with
80% adherence. This is not a figure of speech but arithmetic, and it explains
why the programme metric correlates poorly with recurrence rate.

### 1-2. The valve is destroyed by two arms, and they answer to different masters

```
Immune arm      d(MVA)/dt = -2√(πA) · KFI · VIT              ← what penicillin blocks
Autonomous arm  d(MVA)/dt = -2√(πA) · KFS · (MVG/4) · BRAKE  ← what it cannot block
```

The autonomous arm is **shear stress**. By a simplified Bernoulli relation,
the square of the mean transvalvular velocity is **exactly the mean gradient
(MVG) divided by 4**. And since MVG is proportional to (flow/area)², this arm
**accelerates itself** as the valve narrows.

The model calculates the point where the two arms cross (a computed result,
not an assumption):

| MVA (cm²) | Shear arm (cm²/yr) | Immune arm (cm²/yr) | Dominant |
|---|---|---|---|
| 4.5 | 0.030 | 0.036 | Immune |
| 3.5 | 0.040 | 0.046 | Immune |
| 3.0 | 0.048 | 0.049 | Immune |
| **2.75** | — | — | **Crossover point** |
| 2.5 | 0.058 | 0.051 | Shear |
| 1.5 | 0.092 | 0.048 | Shear |
| 0.8 | 0.131 | 0.040 | Shear |

**Below an MVA of 2.75 cm², the valve's own shear loss exceeds even the
expected immune loss of a patient given no prophylaxis at all.** In other
words, most of the benefit prophylaxis confers on the *valve* lies in the
range **before** stenosis is even measurable — that is latent RHD, the
population echocardiographic screening detects, and exactly the population
the GOAL trial randomised.

At the same time, prophylaxis's effect on preventing *recurrence* is large
across the entire range. A 25-year simulation:

| Regimen | Time protected | Recurrent ARF | Final MVA |
|---|---|---|---|
| No prophylaxis | 0% | 0.85 episodes | 3.28 cm² |
| 4-week interval, full dose | 78.3% | 0.17 episodes | 3.63 cm² |
| 3-week interval, full dose | 92.9% | 0.05 episodes | 3.69 cm² |

**Recurrence falls by 94%, but only 0.41 cm² of valve area is preserved.**
One drug, two endpoints, entirely different effect sizes. A trial powered
for one cannot answer for the other. This model reports the two claims
separately.

### 1-3. The Gorlin block — the patient decompensates while the valve stays unchanged

```
Mean gradient MVG = ( cardiac output / (37.7 · effective area · diastolic filling time) )²
```

In the table below, **valve area is fixed at 1.5 cm² in every row.** Only
demand and rhythm change:

| State | Heart rate | Mean gradient | Left atrial pressure | NYHA |
|---|---|---|---|---|
| Rest, sinus rhythm | 72 | 6.4 mmHg | 12.4 | 1.00 |
| **Pregnancy (CO ×1.5)** | 106 | **23.0 mmHg** | **29.0** | **2.41** |
| New-onset atrial fibrillation | 107 | 11.9 mmHg | 17.9 | 1.00 |
| + metoprolol 100 mg | 99 | 10.8 mmHg | 16.8 | 1.00 |

Pregnancy raises the gradient 3.6-fold without touching the valve, and
pushes left atrial pressure into the pulmonary oedema range. This is why
mitral stenosis "presents" in the second–third trimester of pregnancy.

Reading the same equation in reverse yields an **optimal heart rate**.
Lowering heart rate increases total diastolic filling period (DFP) per
minute, but stroke volume is capped by ventricular capacity, so an optimum
exists, and its value falls as the valve narrows — again a computed result,
not an assumption:

| MVA | Optimal heart rate | Maximum cardiac output at that rate |
|---|---|---|
| 2.0 cm² | 108 bpm | 11.9 L/min |
| 1.5 cm² | 89 bpm | 9.7 L/min |
| 1.0 cm² | 66 bpm | 7.2 L/min |
| 0.8 cm² | 55 bpm | 6.1 L/min |

And one more thing: **in a patient who has already reached the left atrial
pressure ceiling (LAPMAX), rate control cannot lower the gradient.** The
gradient is fixed at LAPMAX − LVEDP, so the diastole that the beta-blocker
buys is spent on forward flow instead. The patient improves, yet **only the
measured gradient fails to move** (cardiac output 4.46 → 4.81 L/min). Below
the ceiling (MVA 1.4, at rest), the same drug lowers the gradient instead.

---

## 2. Incidental Findings

**Steroids are a matter of duration, not dose.** The half-life of
cross-reactive antibody is 60 days, yet the 6-week course used in the trial
covers only 38% of the antibody-time integral. Suppressed valvulitis mostly
just happens later:

| Prednisolone course | Valve area preserved |
|---|---|
| 2 weeks | 5% |
| **6 weeks (duration used in the trial)** | **14%** |
| 12 weeks | 27% |
| 26 weeks | 42% |

The model reproduces the Cochrane null result (PMID 26017576) while
**simultaneously stating what would need to be different.** A 26-week course
has never been trialled.

**The "nine-day window" is a statement about carriage, not about the drug.**
The classic Denny/Wannamaker finding is that treatment within 9 days of
pharyngitis onset still prevents ARF. In this model, a host who clears the
organism quickly is already clear by day 8, so late treatment reduces
antigen exposure by only 6.8%. But in a slow-clearing (carrier-type) host,
the same day-8 treatment reduces it by 76%:

| Natural clearance time | Antigen reduction from day-8 treatment |
|---|---|
| 12 days | 6.8% |
| 21 days | 56.4% |
| 37 days | 75.9% |
| 60 days | 85.3% |

In other words, the size of the window is determined entirely by **how long
an untreated host carries the organism.** This is a testable statement.

---

## 3. Files

| File | Contents |
|---|---|
| [`rhd_qsp_model.dot`](rhd_qsp_model.dot) | Mechanistic map source — **161 nodes, 17 clusters** |
| [`rhd_qsp_model.svg`](rhd_qsp_model.svg) · [`.png`](rhd_qsp_model.png) | Rendered map (150 dpi) |
| [`rhd_mrgsolve_model.R`](rhd_mrgsolve_model.R) | **39 ODEs**, 16 scenarios, mrgsolve |
| [`rhd_shiny_app.R`](rhd_shiny_app.R) | **11-tab** interactive dashboard |
| [`rhd_references.md`](rhd_references.md) | **123 PubMed references** (all verified via API) |
| [`rhd_verify_python.py`](rhd_verify_python.py) | Independent Python/scipy reimplementation + 47 anchors |
| [`rhd_verification_output.txt`](rhd_verification_output.txt) | Verification run output (47/47 passed) |

### The Map's 17 Clusters

Exposure · pharyngitis / host susceptibility / molecular mimicry / acute
rheumatic fever (Jones) / valvulitis / valve remodelling (ratchet) /
**Gorlin block** / left atrium · rhythm · thrombus / pulmonary vasculature ·
right ventricle / systemic haemodynamics · symptoms / BPG depot PK / oral
penicillin · primary prevention / anti-inflammatory therapy / rate control ·
diuresis / anticoagulation / mechanical intervention / endpoints ·
programme metrics

### 39 State Variables (summary)

BPG immediate/sustained-release depot · oral absorption · penicillin central
compartment · **tonsillar effect site** · GAS burden · mucosal immunity ·
antigen · ASO · immune memory · cross-reactive antibody · valvulitis ·
valve oedema · **MVA** · valve calcification · acute/chronic MR · aortic
valve involvement · left atrial volume · AF burden · pulmonary vascular
resistance · right ventricular function · pulmonary congestion · volume
status · CRP · aspirin · prednisolone · beta-blocker · digoxin · warfarin
PK · prothrombin complex · cumulative ARF · cumulative emboli ·
**cumulative time protected (TPROT)**

---

## 4. Verification

All 39 ODEs were **independently reimplemented** in Python/scipy and run
against 47 published/derived anchors (**47/47 passed**). The **7 real
defects** found and fixed in the process are recorded below rather than
silently fixed — knowing what kind of errors hide in a model like this is
as useful as the model itself.

1. **Immune-memory feedback divergence.** Writing the
   `AG → MEM → recurrence risk → AG` loop without saturation gives a loop
   gain of 2.9, which diverges. Fixed by making MEM a capped pool
   (`MEMMAX`). Before this was found, antigen grew without bound.
2. **180-fold immune-cascade scale error.** With the original parameters,
   background exposure alone collapsed MVA from 4.5 → 3.45 cm² in 60 days.
   The entire cascade gain was recalibrated to match a per-episode valve
   loss of 0.35 cm².
3. **An ODE artefact in which GAS never goes extinct.** A continuous
   variable never reaches zero, so once mucosal immunity decayed, a single
   pharyngitis episode recurred forever (∫VIT accumulated to 997 by day
   400). Fixed by introducing an extinction floor (`GEXT`).
4. **Oral penicillin had no effect for lack of an effect compartment.**
   Driving killing off the 30-minute serum half-life leaves almost no
   effective time at a 12-hour dosing interval. A tonsillar effect-site
   compartment (`KEO`, t½ 3 hours) was added. Because the programme's
   0.02 µg/mL threshold is a **plasma** value, TPROT continues to be
   computed from plasma.
5. **Infinite Jacobian in the AFB^0.7 term.** Its derivative is infinite at
   0, and every simulation starts exactly there, so the integrator stalled
   (terminated OOM). Replaced with a saturating Hill form.
6. **Right ventricular function index going negative.** An additive damage
   term sent RVF to −0.43 at a mean pulmonary artery pressure of 53 mmHg.
   Fixed by making the damage term proportional to RVF so it decays toward
   zero.
7. **Heart rate did not respond to demand.** Modelling pregnancy purely
   through cardiac-output demand leaves left atrial pressure short of the
   pulmonary oedema range. Introducing `KHRDEM` and switching from a static
   calculation to a **dynamic** simulation exposed the feedback loop
   congestion → sympathetic drive → heart rate → shortened diastole →
   rising gradient — the calculation itself shows why the sympathetic
   reflex is actually harmful in mitral stenosis.

### Calibrated vs Predicted

**Calibrated (anchored):** BPG depot concentration curve (Kaplan 1989),
penicillin body-size scaling (Hand 2019, Neely 2014), progression rate of
established MS (Sagie 1996, Gordon 1992), the Gorlin constant (1951),
gradient-severity grading (Baumgartner 2009), reduction in recurrence with
BPG (Manyemba 2002), embolic rate in rheumatic AF, warfarin dose-INR,
immediate and 5-year valve area post-PMBV (Iung 1999, Palacios 2002),
reactive pulmonary hypertension in severe MS.

**Predicted (not fitted):** the adherence/protection rank reversal · the
2.75 cm² crossover point · the area-dependent optimal heart rate · the
steroid duration-response · that at the left atrial pressure ceiling rate
control buys cardiac output rather than gradient · that the size of the
nine-day window is determined by carriage duration.

### The Most Exposed Claim

The shear arm is calibrated to a **cohort-average** progression rate of
0.09 cm²/yr at an MVA of 1.5 cm², and then accelerates as the valve narrows
further. But **Sagie 1996 reported an association in the opposite
direction — progression was slower in narrower valves.** As a result, this
model's 10-year post-valvuloplasty valve area (0.81 cm²) is more pessimistic
than the literature's "about 40% restenosis at 10 years."

Sensitivity analysis shows this conclusion is nearly insensitive to the
shear arm's *area dependence* (the 10-year value stays 0.81 even when the
brake coefficient is changed from 1.6 → 20), because the gradient saturates
at LAPMAX − LVEDP as the valve narrows. **So if this trajectory is wrong,
the error lies in the *magnitude* of KFS — that is, in reading the cohort
average as if every valve progressed that way.** A restenosis cohort with
serial echocardiography could settle this directly.

### What It Cannot Explain

The result from INVICTUS (2022, PMID 36036525) that vitamin K antagonists
were superior to rivaroxaban is **encoded only empirically.** This model has
no mechanism to explain why, and does not pretend to.

---

## 5. Usage

```bash
# render the map
dot -Tsvg rhd_qsp_model.dot -o rhd_qsp_model.svg
dot -Tpng -Gdpi=150 rhd_qsp_model.dot -o rhd_qsp_model.png

# verification (requires numpy, scipy) — 47 anchors
python3 rhd_verify_python.py
```

```r
# mrgsolve model and 16 scenarios
source("rhd_mrgsolve_model.R")

# Shiny dashboard (11 tabs)
shiny::runApp("rhd_shiny_app.R")
```

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes.
It is based on public literature but has not been independently validated
or certified, and **must not be used for clinical decision-making,
prescribing, or regulatory submission.** The parameters are illustrative
approximations and require separate fitting and validation against real
patient data.
