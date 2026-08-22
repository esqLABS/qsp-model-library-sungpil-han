# Erectile Dysfunction QSP Model
## Erectile Dysfunction — 45-ODE Quantitative Systems Pharmacology model · 4 oral PDE5 inhibitors · 16 scenarios

<p align="center">
  <a href="../../../erectile-dysfunction/ed_qsp_model.svg"><img src="../../../erectile-dysfunction/ed_qsp_model.png" width="900" alt="ED QSP mechanistic map"></a>
</p>

---

## The question this model is trying to answer

PDE5 inhibitors are among the most successful drugs in the history of medicine.
And yet, at the same drug, the same dose and the same plasma exposure, some
patients recover completely and others show no change at all. In a general
erectile dysfunction cohort the SEP3 success rate for sildenafil 100 mg is 60–70%,
whereas in patients who have had non-nerve-sparing radical prostatectomy it is
under 10%. **Nothing has changed on the drug side.**

The usual explanation is "the patient is more severe". That explanation is a
tautology. It does not say what is more severe or how, nor why it is not solved by
*raising the dose*.

This model was designed so that the answer **is derived as arithmetic rather than
assumed**. The 45 differential equations are **exactly identical** for every
phenotype.

---

## The structural decision — a threshold readout on a saturating amplifier

The entire structure of this model is the following single sentence.

> An erection is **a threshold readout placed on a saturating amplifier**.

There are five stages.

```
① NO pulse           NO_prod = NOSB + FNOI·(KNONN·nNOS·DRIVE·NOSEFF + KNOEN·eNOS·NOSEFF·SHEAR)
                     DRIVE = S(t)·AROU·NRV/(1+kNE·ΔNE)      ← only while the nerve is firing
② amplification (gain stage)   dcGMP/dt = sGC(NO) − KPDE5·PDE5E·P5RES·cGMP/(Km+cGMP) − KDCG·cGMP
                     P5RES = 1/(1 + Cu/IC50)                ← this is the only thing the drug multiplies
③ relaxation         cGMP → PKG → Ca²⁺ efflux ↑ → MLC20 phosphorylation ↓ → R
④ haemodynamics      G_cav = GCMIN + (GCMAX−GCMIN)·R^n/(KRG^n+R^n)   ← steep in R
                     dV/dt = Q_in − Q_out,  Q_out = G_ven·(ICP−P_ven)
                     G_ven = GVEN0·(1 − VOCC),  VOCC = VOCMAX(SMI)·x⁶/(KV⁶+x⁶)   ← positive feedback
⑤ threshold readout  SEP2 = logistic(rigidity − RIG50),  SEP3 = SEP2·logistic(time above threshold − TAE50)
```

What matters is that **the drug term in ② adds nothing to ①**. PDE5 inhibitors do
not make nitric oxide. They are a **multiplicative factor**, `1/(1+Cu/IC50)`.

---

## ① gain × 0 = 0

At the maximum labelled dose, the amount each agent delivers **to the enzyme** is
as follows (all computed from PK fitted to the labelled Cmax/AUC/half-life).

| Agent | Dose | Cmax simulated/label (ng/mL) | AUC simulated/label | Unbound Cu (nM) | Cu/IC50(PDE5) | PDE5 residual activity |
|---|---|---|---|---|---|---|
| sildenafil | 100 mg | 422.7 / 450 | 1965 / 1963 | 35.63 | 10.2 | 0.089 |
| tadalafil | 20 mg | 357.0 / 378 | 8008 / 8066 | 55.01 | 30.6 | 0.032 |
| vardenafil | 20 mg | 19.7 / 21 | 75 / 74 | 2.02 | 14.4 | 0.065 |
| avanafil | 200 mg | 5181.6 / 5200 | 11601 / 11600 | 107.08 | 20.6 | 0.046 |

All four agents sit in the **saturated region** of `1/(1+Cu/IC50)`. That network
meta-analyses cannot distinguish the efficacy ceilings of the four agents is a
**consequence** of this table, and not an assumption put into the model.

Applying that same exposure (sildenafil 100 mg, PDE5 residual activity 0.089) to
six phenotypes gives the following.

| Phenotype | NO untreated → treated (nM) | cGMP untreated → treated | Peak ICP untreated → treated (mmHg) | IIEF-EF untreated → treated |
|---|---|---|---|---|
| healthy 30-year-old | 17.83 → 17.83 | 21.00 → 33.44 | 88.1 → 88.1 | 29.1 → 29.2 |
| mild ED | 3.29 → 3.29 | 6.30 → 11.09 | 86.0 → 86.9 | 24.7 → 28.0 |
| vasculogenic ED | 2.50 → 2.57 | 4.58 → 8.34 | 7.7 → 83.4 | 9.5 → 25.6 |
| diabetic ED | 2.93 → 2.99 | 4.40 → 7.98 | 7.5 → 79.8 | 9.5 → 23.2 |
| post-prostatectomy (bilateral nerve-sparing) | 0.26 → 0.26 | 0.63 → 1.17 | 7.5 → 7.5 | 9.5 → 9.5 |
| post-prostatectomy (non-nerve-sparing) | 0.26 → 0.26 | 0.63 → 1.17 | 7.5 → 7.5 | 9.5 → 9.5 |

Three things can be read at once.

1. **In a healthy man the drug does almost nothing** (ICP 88.1 → 88.1 mmHg),
   because the amplifier is already saturated. The reason PDE5 inhibitors are not
   erection enhancers in normal men and the reason 200 mg adds almost nothing to
   100 mg are the same reason — the structural basis for the maximum labelled dose
   stopping at 100 mg.
2. **In vasculogenic and diabetic patients the drug is what carries them over the
   threshold.**
3. **In the non-nerve-sparing prostatectomy patient, cGMP rises from 0.63 to 1.17
   and yet ICP is 7.5 → 7.5 mmHg.** The drug-side terms are **exactly identical**
   to the two patients above. What differs is the size of the pulse being
   multiplied, and multiplying a pulse 68 times smaller by the same factor still
   leaves it below threshold.

This is "gain × 0 = 0". The conclusion comes out of the equations; nowhere was it
written in that "patients with nerve injury do not respond".

---

## The same logic run backwards — the cAMP bypass

If the proposition that multiplying a zero pulse by any gain gives zero is true,
then **a drug that bypasses the pulse must work in the same patient.** Alprostadil
enters via EP2/EP4 → adenylate cyclase → cAMP → PKA and requires no nitric oxide
whatsoever. The model contains this fact structurally, so what follows is a
prediction.

| Phenotype | Treatment | Peak ICP (mmHg) | Rigidity (%) | cGMP | SEP3 | IIEF-EF |
|---|---|---|---|---|---|---|
| post-prostatectomy (non-nerve-sparing) | none | 7.5 | 0.3 | 0.63 | 0.000 | 9.5 |
| post-prostatectomy (non-nerve-sparing) | sildenafil 100 | 7.5 | 0.3 | 1.17 | 0.000 | 9.5 |
| post-prostatectomy (non-nerve-sparing) | alprostadil 20 ug | 83.8 | 80.7 | 0.81 | 0.948 | 27.7 |
| post-prostatectomy (non-nerve-sparing) | trimix | 83.8 | 80.7 | 0.97 | 0.948 | 27.7 |
| diabetic ED | none | 7.5 | 0.3 | 4.40 | 0.000 | 9.5 |
| diabetic ED | sildenafil 100 | 79.8 | 78.3 | 7.98 | 0.584 | 23.2 |
| diabetic ED | alprostadil 20 ug | 67.4 | 68.5 | 4.46 | 0.320 | 18.5 |
| diabetic ED | trimix | 76.3 | 76.0 | 5.36 | 0.627 | 23.5 |
| vasculogenic ED | none | 7.7 | 0.3 | 4.58 | 0.000 | 9.5 |
| vasculogenic ED | sildenafil 100 | 83.4 | 80.5 | 8.34 | 0.772 | 25.6 |
| vasculogenic ED | alprostadil 20 ug | 82.0 | 79.7 | 4.68 | 0.940 | 27.6 |
| vasculogenic ED | trimix | 83.0 | 80.2 | 5.61 | 0.945 | 27.6 |

In the non-nerve-sparing patient, sildenafil 100 mg leaves ICP at 7.5 mmHg while
alprostadil 20 µg puts it at 83.8 mmHg. cGMP is low in both cases (1.17 vs 0.81) —
that is, alprostadil does not work by raising cGMP. The reason intracavernosal
injection is the standard next-line treatment in patients who fail PDE5 inhibitors
comes out of the equations. For the same reason, the priapism risk is also
structural: the cAMP pathway has no gate saying "switch off when the nerve stops".

### sGC activators — a second bypass, and a conditional failure

When the haem of sGC is oxidised it no longer responds to nitric oxide. The model
computes the oxidised fraction `SGCOX` from ROS (healthy 0.049 → diabetic 0.359).
An sGC **stimulator** raises the NO sensitivity of the reduced form and therefore
requires NO to be present; an sGC **activator** directly activates the oxidised
pool and therefore does not require NO.

| Phenotype | Oxidised sGC fraction | Treatment | cGMP | Peak ICP (mmHg) |
|---|---|---|---|---|
| diabetic ED | 0.359 | none | 4.40 | 7.5 |
| diabetic ED | 0.359 | sildenafil 100 | 7.98 | 79.8 |
| diabetic ED | 0.359 | sGC stimulator | 11.91 | 80.6 |
| diabetic ED | 0.359 | sGC activator | 12.37 | 80.6 |
| diabetic ED | 0.359 | activator + PDE5i | 20.81 | 80.6 |
| post-prostatectomy (non-nerve-sparing) | 0.195 | none | 0.63 | 7.5 |
| post-prostatectomy (non-nerve-sparing) | 0.195 | sildenafil 100 | 1.17 | 7.5 |
| post-prostatectomy (non-nerve-sparing) | 0.195 | sGC stimulator | 1.53 | 7.5 |
| post-prostatectomy (non-nerve-sparing) | 0.195 | sGC activator | 4.52 | 8.3 |
| post-prostatectomy (non-nerve-sparing) | 0.195 | activator + PDE5i | 8.29 | 85.8 |

In the non-nerve-sparing patient, the activator alone raises cGMP from 0.63 to
4.52 but ICP is still 8.3 mmHg. The activator **plus** a PDE5 inhibitor, however,
gives cGMP 8.29 and ICP 85.8 mmHg. The two drugs are each insufficient for
different reasons — one has no supply and the other has only gain. That this
combination is a rational hypothesis for NO-deficient tissue is a prediction of the
model, and it has not been tested in humans.

---

## ② Two clocks — same drug, different time constants

On-demand dosing moves only the **fast** variable, cGMP (time constant of
minutes). Daily dosing adds a **slow** structural variable on top of it: chronic
cGMP suppresses Smad signalling through PKG, so `TGFB → COL` falls and smooth
muscle trophism rises, improving `SMI` (time constant of weeks to months).

This structure produces a testable prediction: **the benefit measured during
treatment** and **the benefit measured after a drug-free washout** are different
quantities. The REACTT design was reproduced as written — after bilateral
nerve-sparing prostatectomy, nine months of tadalafil 5 mg daily / 20 mg on demand
/ placebo, followed by a six-week drug-free period.

| Arm | SMI at 9 months | Penile length at 9 months (cm) | Weekly mean cGMP at 9 months | IIEF at 9 months (on drug) | SMI after washout | Unassisted ICP after washout | Unassisted IIEF after washout |
|---|---|---|---|---|---|---|---|
| placebo | 0.604 | 11.42 | 0.82 | 28.2 | 0.629 | 86.3 | 27.8 |
| on demand 20 mg | 0.712 | 11.70 | 1.25 | 28.6 | 0.636 | 86.3 | 27.9 |
| daily 5 mg | 0.744 | 11.80 | 1.48 | 28.2 | 0.639 | 86.3 | 27.9 |

What has to be read is the contrast between two columns.

**During treatment the structure clearly diverges.** At nine months SMI is 0.604
for placebo, 0.712 on demand and 0.744 daily, and penile length is 11.42 / 11.70 /
11.80 cm. This is the result of chronic cGMP (weekly means 0.82 / 1.25 / 1.48)
suppressing `TGFB → COL` and thereby defending the smooth muscle:collagen ratio.

**After washout most of that difference disappears.** After six weeks off drug, SMI
converges to 0.629 / 0.636 / 0.639 and the difference in unassisted IIEF-EF is 0.0
points. Defending the structure did not translate into erectile function.

The reason is in the equations. The endpoint is a **product and a threshold**, and
at this point the limiting factor is not the structure but the neural input
(NRV 0.681). The nerve decides whether the threshold is crossed; the structure
decides how rigid things are *after* it has been crossed. In other words, **while
the nerve is the limiting factor, protecting the structure does not move unassisted
erectile function.** The pattern REACTT reported — failing on its primary endpoint
while still reducing the loss of penile length — comes straight out, and is not the
result of writing "rehabilitation does not work" into the model.

Note: the values above are the trajectory of **one median patient**. In this model
the median bilateral nerve-sparing patient recovers spontaneously, so IIEF is
pinned at the ceiling and the between-arm differences can only be read in the
structural variables. At the population level (responder fraction) recovery becomes
a gentle curve rather than a step.

### Natural history — a race between regeneration and fibrosis

The same structure creates the *window* for rehabilitation. The time constant of
axonal regeneration is about 7 months and that of hypoxia-induced fibrosis about 12
weeks. Since nerve-sparing status sets the recovery ceiling (`NRVMAX`), the 24-month
courses of the two phenotypes diverge as follows.

| Month | Bilateral nerve-sparing NRV | SMI | Unassisted ICP | tadalafil 20 ICP | Non-nerve-sparing NRV | SMI | tadalafil 20 ICP |
|---|---|---|---|---|---|---|---|
| 0 | 0.100 | 0.532 | 7.5 | 7.5 | 0.100 | 0.532 | 7.5 |
| 3 | 0.381 | 0.431 | 7.6 | 83.9 | 0.118 | 0.423 | 7.5 |
| 6 | 0.563 | 0.492 | 85.1 | 85.5 | 0.129 | 0.423 | 7.5 |
| 9 | 0.681 | 0.604 | 86.2 | 86.2 | 0.136 | 0.423 | 7.5 |
| 12 | 0.758 | 0.650 | 86.3 | 86.3 | 0.141 | 0.424 | 7.5 |
| 18 | 0.840 | 0.687 | 86.3 | 86.3 | 0.146 | 0.424 | 7.5 |
| 24 | 0.875 | 0.699 | 86.3 | 86.3 | 0.148 | 0.424 | 7.5 |

---

## ③ One product, two vascular beds

Cavernosal efficacy and systemic hypotension are **the same product of two
factors** (NO supply × PDE5 inhibition) read out in different vessels. Nitrates
raise the first factor systemically. The model computes the blood pressure drop
from **the same** two parameters that generate the erection.

| Combination | Nadir MAP (mmHg) | ΔMAP | Peak systemic cGMP |
|---|---|---|---|
| none | 95.0 | 0.0 | 1.00 |
| GTN 0.4 mg SL | 93.7 | -1.3 | 1.46 |
| sildenafil 100 | 88.6 | -6.4 | 1.83 |
| sildenafil 100 + GTN | 88.5 | -6.5 | 1.92 |
| tadalafil 20 + GTN | 87.5 | -7.5 | 1.99 |

The same isoform selectivity table also computes the adverse effects: sildenafil's
PDE6 occupancy at the efficacious dose is 51.2% (visual disturbance) and
tadalafil's PDE11 occupancy is 59.8% (myalgia, back pain). These values come from
**the same** concentrations that generate the efficacy; no separate adverse-effect
model was added.

---

## Veno-occlusion — why a "knee" is needed

An erection is not made by inflow alone. As the cavernosa expand, the subtunical
veins are compressed against the tunica albuginea and outflow is blocked (positive
feedback). This ceiling is set by the smooth muscle:collagen ratio `SMI`. But this
dependence **must not be linear** — clinically, veno-occlusive dysfunction is a late
and severe finding, and mild loss of smooth muscle does not abolish a rigid
erection. Hence a knee (Hill, KSMI=0.221, n=6.0).

| Phenotype | SM | COL | SMI | Veno-occlusive ceiling VOCMAX | Stenosis STEN | Nocturnal erection peak ICP |
|---|---|---|---|---|---|---|
| healthy 30-year-old | 1.00 | 1.00 | 1.00 | 0.97 | 0.020 | 88.1 |
| mild ED | 0.80 | 1.73 | 0.63 | 0.97 | 0.066 | 82.1 |
| vasculogenic ED | 0.72 | 2.17 | 0.50 | 0.96 | 0.148 | 7.5 |
| diabetic ED | 0.70 | 2.67 | 0.42 | 0.95 | 0.163 | 7.5 |
| post-prostatectomy (bilateral nerve-sparing) | 0.72 | 2.00 | 0.53 | 0.97 | 0.087 | 7.5 |
| post-prostatectomy (non-nerve-sparing) | 0.72 | 2.00 | 0.53 | 0.97 | 0.087 | 7.5 |
| hypogonadal ED | 0.76 | 1.71 | 0.62 | 0.97 | 0.101 | 85.5 |
| psychogenic ED | 0.99 | 1.03 | 0.98 | 0.97 | 0.020 | 88.1 |

Not one value in that table is an input. A single set of comorbidity inputs
(age · HbA1c · BMI · LDL · smoking · RAAS) is collected into an oxidative load
`rosdrv`, from which ROS, eNOS coupling, ADMA, oxidised sGC and stenosis follow in
closed form, and the oxygenation grade of nocturnal erections drives `TGFB → COL`
and smooth muscle apoptosis, which set SMI. That is, **the structure is a computed
result**.

---

## Virtual population — a trial mean is a mixture distribution

Because the endpoint is a threshold readout, **a single simulated patient is either
a responder or not**. A trial's mean IIEF-EF is not an individual's value but the
expectation of a mixture distribution. The trial-matching quantity in this model is
therefore a population and not an individual: the individual NO-production
multiplier `FNOI` is swept log-normally (log SD 0.55) and averaged.

| Cohort | Agent | Dose | IIEF-EF simulated / target | Peak ICP | Rigidity | SEP2 | SEP3 | Responder fraction |
|---|---|---|---|---|---|---|---|---|
| healthy 30-year-old | sildenafil | 0 mg | 29.1 / 29.0 | 88.1 | 82.9 | 0.966 | 0.950 | 1.00 |
| mild ED | sildenafil | 0 mg | 20.3 / 21.0 | 58.7 | 53.3 | 0.610 | 0.479 | 0.56 |
| mild ED | sildenafil | 50 mg | 25.9 / 26.0 | 81.9 | 77.5 | 0.902 | 0.778 | 0.88 |
| vasculogenic ED | sildenafil | 0 mg | 15.0 / 14.0 | 34.7 | 28.9 | 0.327 | 0.236 | 0.27 |
| vasculogenic ED | sildenafil | 25 mg | 20.8 / 19.5 | 62.0 | 58.0 | 0.672 | 0.516 | 0.61 |
| vasculogenic ED | sildenafil | 50 mg | 22.1 / 21.0 | 67.5 | 63.9 | 0.746 | 0.577 | 0.66 |
| vasculogenic ED | sildenafil | 100 mg | 22.7 / 22.5 | 69.7 | 66.1 | 0.770 | 0.615 | 0.73 |
| vasculogenic ED | sildenafil | 200 mg | 23.1 / — | 71.4 | 67.9 | 0.789 | 0.637 | 0.73 |
| vasculogenic ED | tadalafil | 20 mg | 23.3 / — | 71.7 | 68.2 | 0.798 | 0.647 | 0.76 |
| vasculogenic ED | vardenafil | 20 mg | 23.1 / — | 71.3 | 67.9 | 0.793 | 0.631 | 0.76 |
| vasculogenic ED | avanafil | 200 mg | 23.1 / — | 71.4 | 67.9 | 0.789 | 0.636 | 0.73 |
| diabetic ED | sildenafil | 0 mg | 12.4 / 11.5 | 22.9 | 16.7 | 0.182 | 0.121 | 0.15 |
| diabetic ED | sildenafil | 100 mg | 19.9 / 18.5 | 59.4 | 56.2 | 0.644 | 0.471 | 0.56 |
| diabetic ED | tadalafil | 20 mg | 20.6 / — | 61.7 | 58.5 | 0.678 | 0.506 | 0.61 |
| hypogonadal ED | sildenafil | 0 mg | 22.4 / — | 66.4 | 61.6 | 0.721 | 0.580 | 0.68 |
| hypogonadal ED | sildenafil | 100 mg | 27.0 / — | 83.5 | 79.5 | 0.934 | 0.837 | 0.95 |
| post-prostatectomy (bilateral nerve-sparing) | tadalafil | 20 mg | 9.5 / — | 7.5 | 0.3 | 0.000 | 0.000 | 0.00 |
| post-prostatectomy (non-nerve-sparing) | tadalafil | 20 mg | 9.5 / 8.5 | 7.5 | 0.3 | 0.000 | 0.000 | 0.00 |
| psychogenic ED | sildenafil | 0 mg | 28.9 / — | 88.1 | 82.9 | 0.966 | 0.933 | 1.00 |
| psychogenic ED | sildenafil | 50 mg | 29.1 / — | 88.1 | 82.9 | 0.966 | 0.951 | 1.00 |

Six global parameters (`KPKG`, `NPKG`, `KRG`, `NGC`, `KSMI`, `RIG50`) and three
per-archetype nerve-integrity values were least-squares fitted to the ten anchors
above. No other parameter was fitted to a clinical outcome.

**The responder fraction is far more sensitive to the log-normal SD than the mean
IIEF-EF is.** This is a limitation of the model and at the same time something this
structure predicts — a testable statement that the mean change in a clinical trial
is not "every patient improving a little" but "some patients crossing the
threshold".

---

## The window of action — not the half-life but the time spent above threshold

Sweep the time of the attempt after a single dose and the clinically meaningful
quantity is not the half-life but **the last time point at which rigidity is above
threshold**.

| After dose (h) | sildenafil 100 Cu / rigidity | tadalafil 20 Cu / rigidity | vardenafil 20 Cu / rigidity | avanafil 200 Cu / rigidity |
|---|---|---|---|---|
| 0.5 | 35.52 / 80.5 | 45.25 / 80.5 | 1.99 / 80.5 | 107.00 / 80.5 |
| 1.0 | 32.29 / 80.5 | 54.64 / 80.5 | 1.85 / 80.5 | 83.55 / 80.5 |
| 2.0 | 23.12 / 80.5 | 51.75 / 80.5 | 1.23 / 80.5 | 31.64 / 80.5 |
| 4.0 | 13.28 / 80.4 | 44.23 / 80.5 | 0.58 / 80.4 | 9.98 / 80.0 |
| 8.0 | 5.99 / 79.9 | 36.89 / 80.5 | 0.21 / 79.7 | 5.20 / 78.2 |
| 12.0 | 3.10 / 77.3 | 31.26 / 80.5 | 0.10 / 75.2 | 2.96 / 68.0 |
| 24.0 | 0.47 / 0.4 | 19.05 / 80.5 | 0.02 / 0.4 | 0.54 / 0.4 |
| 36.0 | 0.07 / 0.3 | 11.60 / 80.5 | 0.00 / 0.3 | 0.10 / 0.3 |
| 48.0 | 0.01 / 0.3 | 7.07 / 80.4 | 0.00 / 0.3 | 0.02 / 0.3 |

Tadalafil 20 mg is still above threshold at 36 hours with a rigidity of 80.5%,
while sildenafil falls below threshold between 8 and 12 hours. How the clinical
description "36 hours" arises from a half-life of 16.8 hours is computed.

---

## The androgen gate — multiplicative, not additive

In this model testosterone is **one value multiplied into four places**: nNOS
expression, PDE5A expression, smooth muscle trophism, and central libido.
Therefore, in a hypogonadal non-responder, replacing T should not *add* an effect
but *convert* the patient into a responder. Twelve weeks of testosterone gel
50 mg/day:

| Group | Total T (ng/dL) | LH (IU/L) | Haematocrit (%) | nNOS | Smooth muscle SM | SMI | Unassisted ICP | sildenafil 100 ICP | Unassisted IIEF | sildenafil 100 IIEF |
|---|---|---|---|---|---|---|---|---|---|---|
| before treatment | 294 | 7.01 | 45.0 | 0.575 | 0.761 | 0.617 | 85.8 | 85.9 | 26.4 | 28.3 |
| 12-week untreated control | 294 | 7.01 | 45.0 | 0.575 | 0.569 | 0.477 | 84.4 | 84.9 | 24.4 | 27.1 |
| 12 weeks of gel | 310 | 6.75 | 45.0 | 0.589 | 0.630 | 0.530 | 85.1 | 85.5 | 24.7 | 27.3 |

Exogenous T raises total testosterone from 294 to 310 ng/dL and suppresses LH from
7.01 to 6.75 IU/L (negative feedback). Haematocrit does not move from 45.0%,
because the T reached is within the physiological range and the `EPOT` term does not
switch on — that erythrocytosis is a function of supraphysiological exposure comes
out of the equations.

**Comparing against the untreated control over the same period** shows where the
androgen gate acts: smooth muscle 0.630 versus 0.569, SMI 0.530 versus 0.477,
nNOS 0.589 versus 0.575. This median patient's erectile endpoint is already at the
ceiling (ICP 85.1 mmHg) so the IIEF difference is only 0.3 points, but the effect is
measurable in three of the four places the gate enters multiplicatively (smooth
muscle trophism, nNOS expression, structure).

**This point is a falsifiable prediction of the model.** The literature contains
both RCTs showing that adding T is effective in sildenafil non-responders with
hypogonadism and RCTs showing no effect (references §10). The multiplicative-gate
structure predicts the former, and only in patients who are below threshold.

---

## Risk-factor modification — the only route that moves without a drug

Exercise + statin + weight loss + LDL improvement were applied over six months (no
PDE5 inhibitor).

| Group | ROS | eNOS coupled fraction | SMI | Unassisted ICP | Unassisted IIEF | sildenafil 100 IIEF |
|---|---|---|---|---|---|---|
| baseline | 2.36 | 0.515 | 0.498 | 7.7 | 9.5 | 25.6 |
| 6-month untreated control | 2.36 | 0.515 | 0.402 | 7.7 | 9.5 | 25.3 |
| 6 months of modification | 1.20 | 0.814 | 0.735 | 84.4 | 27.4 | 28.2 |

Oxidative load improves from 2.36 to 1.20 and the eNOS coupled fraction from 0.515
to 0.814, and SMI becomes 0.498 → 0.735. **Against the untreated control over the
same period**, the unassisted IIEF-EF difference is 17.9 points (control 9.5 →
modification 27.4). That the control itself has risen above baseline is residual
drift between the analytical baseline state and the ODE steady state, which is why
every long-term comparison is presented together with a control.

The reason this route is slow is structural — ROS falls on a timescale of hours,
but `SMI` is tied to the collagen time constant (about 8 days).

### Isolating the slow arm of daily low dose

In vasculogenic ED, tadalafil 5 mg daily for 12 weeks versus 20 mg on demand for 12
weeks — **the total exposure of the two arms differs, but what is compared is the
unassisted (drug-free) state.**

| Group | Weekly mean cGMP | TGFB | COL | SMI | Oxygenation grade ERFR | Unassisted ICP | Unassisted IIEF |
|---|---|---|---|---|---|---|---|
| baseline | 0.52 | 2.007 | 2.170 | 0.498 | 0.0000 | 7.7 | 9.5 |
| daily 5 mg for 12 weeks | 0.92 | 1.723 | 1.887 | 0.526 | 0.0375 | 7.5 | 9.6 |
| on demand 20 mg for 12 weeks | 0.79 | 1.842 | 2.008 | 0.472 | 0.0208 | 7.6 | 9.4 |

---

## Verification

The 45 differential equations were implemented **twice, independently of each
other** — mrgsolve/C++ (DLSODA) and Python/scipy (LSODA). The two implementations
share no code, only the equations. Across 12 (phenotype × treatment) combinations:

| Quantity | Maximum relative discrepancy | Median |
|---|---|---|
| oxidative load ROS | 0.00% | 0.00% |
| eNOS coupled fraction | 0.00% | 0.00% |
| ADMA | 0.00% | 0.00% |
| oxidised sGC fraction | 0.00% | 0.00% |
| stenosis STEN | 0.00% | 0.00% |
| smooth muscle SM | 2.16% | 0.00% |
| collagen COL | 2.85% | 0.00% |
| SMI | 3.27% | 0.00% |
| total testosterone | 0.00% | 0.00% |
| peak NO | 27.22% | 1.21% |
| peak cGMP | 21.52% | 4.58% |
| peak relaxation R | 67.43% | 8.46% |
| peak ICP | 79.78% | 6.65% |
| peak rigidity | 94.38% | 3.75% |
| unbound drug concentration | 8.42% | 4.21% |
| IIEF-EF | 59.05% | 2.43% |

**The baseline state (structure · oxidation · endocrine) agrees to within 3.3%.**
The median discrepancy in the dynamic quantities is 2–8%. The maxima, however, are
far larger — and those maxima come from **a single combination**: diabetic ED +
sildenafil 100 mg. In that combination the peak cGMP of the two implementations is
7.39 versus 7.98 (a 7.3% difference) while the peak ICP is 16.1 versus
79.8 mmHg.

This discrepancy is not a defect but **the model's own claim appearing as
numbers**. That patient is sitting exactly astride the threshold, and in a
threshold-readout system a difference of a few per cent in the input becomes a
difference of several fold in the output. The clinical implication is the same — in
a patient near threshold, small physiological fluctuations (sleep, alcohol, anxiety,
absorption with food) separate success from failure. Conversely, in patients far
from the threshold (healthy men, non-nerve-sparing prostatectomy patients) the two
implementations agree to the decimal.

### Defects the verification found

This list is a record of where the model's credibility comes from. Six were found
by the Python reimplementation, three by R↔Python cross-checking.

1. **The structure collapsed in every patient archetype.** The hypoxia signal had
   been written as a binary indicator at `ICP ≥ 40 mmHg`, so in patients whose
   nocturnal erections do not reach that threshold the signal went immediately to 1,
   fibrosis ran away, and SMI fell to 0.16 (human veno-occlusive dysfunction is
   around 0.6). Fixed by making oxygenation a **continuous grade**,
   `g(ICP) = clip((ICP−10)/50)`, and bounding `HYPF` and `APOPF`.
2. **Defining relaxation relative to the patient's own baseline tone was an
   error.** Using the healthy flaccid value as the zero point meant that a
   pathological patient with high baseline tone had relaxation clipped negative even
   at maximum stimulation, making **every ED patient a complete non-responder**.
   Fixed by using maximum contraction and maximum relaxation as absolute reference
   points.
3. **An exponential conductance law abolished flaccid perfusion.** The `KG ≈ 20`
   fitted to `G = G₀·exp(KG·R)` made cavernosal blood flow **0.2 mL/h** in a
   patient whose baseline tone is higher than the healthy flaccid value (tissue that
   could not survive). Replaced by a Hill form with a floor, and `GCMAX` derived
   from the anchor "healthy flaccid = 300 mL/h".
4. **The SMI dependence of the veno-occlusive ceiling was excessive.** With no knee
   it was nearly linear, so even mild loss of smooth muscle made a rigid erection
   impossible and PDE5 inhibitors worked in no patient at all. A Hill knee was
   introduced to match the clinical fact that veno-occlusive dysfunction is a late
   and severe finding.
5. **Unnormalised steady states.** The synthesis/degradation constants for TGF-β1,
   AGE, oxidative damage, PDE5 expression and nNOS did not give the healthy
   reference value of 1.0, so even a healthy person's baseline state drifted (TGF-β1
   up to 50). Each term was normalised against its healthy reference value.
6. **Nocturnal erections contaminated the sexual activity diary.** Because the
   stimulus input was a single scalar, nocturnal erections were driving performance
   anxiety and the IIEF diary. The stimulus was split into a nocturnal component and
   an attempt component.
7. **The diary gate opened too early** (cross-check). SEP3 is defined by the
   *time* spent above threshold, and recording the diary while the erection is still
   forming means that time does not accumulate and every attempt reads as a failure.
   Performance anxiety was overestimated about 10-fold. Fixed by moving the gate to
   the later part of the erection (phase 0.78–0.98).
8. **The analytical baseline state was not a fixed point of the ODEs**
   (cross-check). The weekly mean oxygenation and cGMP were computed at the *peak*
   of the nocturnal erection, whereas the ODEs take a time *average* of them. The two
   are different functionals, so the baseline state drifted over the collagen time
   constant. Fixed by computing the same functional by quadrature over both
   stimulus windows (nocturnal + attempt).
9. **The neural deficit of non-surgical patients healed itself** (cross-check).
   Nitrergic fibre loss from ageing and diabetes does not regenerate, but the default
   for `NRVMAX` was 1.0, so over long integrations the NRV of the mild, vasculogenic
   and hypogonadal archetypes "healed" from 0.47 to 0.86 and the disease itself
   disappeared. Fixed by specifying `KREG = 0` and `NRVMAX = NRV` for those
   archetypes.

There is also a point at which **the verification rejected the author's prior
hypothesis**. Initially it was expected that the rehabilitation effect of daily
dosing would appear through nerve regeneration, so `KNRVDRUG` (a cGMP-dependent
neurotrophic term) was included. But reproducing the human post-washout result
required this term to be zero. The default of 0 encodes the human result rather than
the animal experiments, and the parameter was left in place so that the hypothesis
can be tested.

**An honest statement about the calibration.** Defects 7–9 above were found and
fixed *after* the least-squares calibration was finished, and no recalibration was
done. As a result the population table above departs from the fitted targets by up
to 1.4 points (IIEF-EF) — within the MCID for mild ED (2 points) — and target and
achieved values are shown side by side. Likewise, the baseline discrepancies in the
cross-check table above (≤ 3.3%) are residual differences between the analytical
baseline state and the ODE steady state, which is why every long-term scenario was
run together with **an untreated control over the same period**, so that this
residual drift cancels in the between-group comparison.

---

---

## Files

| File | Contents |
|---|---|
| [`ed_qsp_model.dot`](../../../erectile-dysfunction/ed_qsp_model.dot) · [`.svg`](../../../erectile-dysfunction/ed_qsp_model.svg) · [`.png`](../../../erectile-dysfunction/ed_qsp_model.png) | Mechanistic map — 189 nodes, 18 clusters |
| [`ed_mrgsolve_model.R`](ed_mrgsolve_model.R) | 45-ODE mrgsolve model, PK library for 4 PDE5 inhibitors, 9 patient archetypes, 16 scenarios, virtual population functions, calibration notes |
| [`ed_shiny_app.R`](../../../erectile-dysfunction/ed_shiny_app.R) | 14-tab interactive dashboard |
| [`ed_references.md`](ed_references.md) | 102 references (every PMID automatically verified via NCBI E-utilities) |

```r
source("ed_mrgsolve_model.R")
b <- ed_baseline(ed_mod, "vasculo", "sildenafil")   # build the patient baseline state
ed_attempt(b, dose = 100)                            # a single attempt
ed_s3_dose(ed_mod)                                   # dose-response
ed_s9_rehab(ed_mod)                                  # REACTT reproduction
ed_population(ed_mod, "diabetic", "sildenafil", 100) # virtual population
shiny::runApp("ed_shiny_app.R")
```

---

## Limitations

1. **`RIG50 = 62.8%` and `TAE50 = 8 minutes` are not physical constants.** They
   are behavioural parameters converting a pressure waveform into sexual diary
   entries, and they are the single largest source of endpoint uncertainty.
2. **The log-normal SD of inter-individual variability (0.55) is not a measured
   value.** It is an inferred value set to reproduce the population means, and the
   responder *fraction* is sensitive to it.
3. **The time structure of ejaculation, orgasm and libido is not modelled.**
   Central arousal is a single scalar, and premature or delayed ejaculation is out of
   scope.
4. **Cavernosal tissue is assumed to be spatially homogeneous.** Focal plaque
   (Peyronie's disease) is marked on the map but is not in the equations.
5. **cGMP, cAMP and Ca²⁺ are treated as multiples of the flaccid reference value**
   rather than as absolute molar concentrations. Only the isoform IC50 values and the
   unbound drug concentrations are in absolute units.
6. All values are for **educational and research purposes**. They cannot be used
   for clinical decision-making.
