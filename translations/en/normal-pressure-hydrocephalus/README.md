# Idiopathic Normal Pressure Hydrocephalus (iNPH) — Quantitative Systems Pharmacology Model

[![Mechanistic Map](../../../normal-pressure-hydrocephalus/inph_qsp_model.png)](../../../normal-pressure-hydrocephalus/inph_qsp_model.svg)

> **Summary in one sentence.** iNPH is not a pressure disease but **a resistance
> disease**, and that resistance can only be observed through **a lossy transducer
> that depends on compliance**. And what sets the safety margin of treatment is not
> biology but **hydrostatics**.

---

## 0. The thesis

The name "normal pressure hydrocephalus" is a clinical observation, but it is at
the same time **the definition of a modelling problem**. The Davson equation makes
mean intracranial pressure an **affine function** of outflow resistance.

```
ICP = Pss + If · Rout              (Davson)
C   = 1 / (E1 · (P − P0))          (Marmarou)
AMP = ΔVp / C                      (pulse amplitude)
```

The intercept `Pss` (dural venous sinus pressure) **carries no information at all
about the disease.** Therefore

```
d ln ICP / d ln Rout = If·Rout / (Pss + If·Rout)  <  1     always
```

In this model's patient that value is **0.506** — **49% of the Rout signal is
thrown away.** The pulse amplitude has no such intercept. Because
`AMP = ΔVp·E1·(P − P0)`, it **inherits** the pressure rise and **multiplies** it by
the rise in elastance. So the same disease gives

| | Healthy 75-year-old | iNPH | Ratio |
|---|---|---|---|
| **Rout** [mmHg/(mL/min)] | 9.00 | 19.00 | **2.11×** |
| Mean ICP (supine) [mmHg] | 9.65 | **13.15** | 1.36× |
| Compliance C (supine) [mL/mmHg] | 0.550 | 0.297 | 0.54× |
| **ICP pulse amplitude AMP** [mmHg] | 2.00 | **4.81** | **2.41×** |
| Pulsatile transmantle gradient [mmHg] | 0.052 | 0.268 | 5.17× |
| Periventricular white matter perfusion [mL/100g/min] | 20.9 | 15.9 | 0.76× |
| Gait speed [m/s] | 1.05 | 0.61 | — |
| MMSE | 28.9 | 24.9 | — |

**13.15 mmHg is inside the normal range of 7–15 mmHg.** The resistance became
2.11-fold while the mean pressure became only 1.36-fold, and the same disease
raised the pulse amplitude 2.41-fold, past the literature threshold of "AMP >
4 mmHg predicts shunt response". **The disease barely appears in the very variable
its name points at.**

---

## 1. Deliverables

| File | Contents |
|---|---|
| [`inph_qsp_model.dot`](../../../normal-pressure-hydrocephalus/inph_qsp_model.dot) · [`.svg`](../../../normal-pressure-hydrocephalus/inph_qsp_model.svg) · [`.png`](../../../normal-pressure-hydrocephalus/inph_qsp_model.png) | Mechanistic map — **20 clusters, 255 nodes** |
| [`inph_reference_model.py`](../../../normal-pressure-hydrocephalus/inph_reference_model.py) | **The reference implementation.** Pure standard-library Python, 45 ODEs + 1 closed-form fast state. Every number below is computed here |
| [`inph_model_report.txt`](../../../normal-pressure-hydrocephalus/inph_model_report.txt) | The **committed output** of that script. The source of every figure in this README |
| [`inph_mrgsolve_model.R`](inph_mrgsolve_model.R) | mrgsolve implementation (46 compartments), 21 scenarios, titration-map and delay-sweep functions |
| [`inph_shiny_app.R`](inph_shiny_app.R) | 13-tab interactive dashboard |
| [`inph_references.md`](inph_references.md) | **156 PubMed links, all with author, year and journal verified against the E-utilities API.** §19 conversely states which parts are *not* supported by the literature |

### Model structure (45 ODEs + a closed-form hydraulics block)

| Block | Compartments | State variables |
|---|---|---|
| Craniospinal hydraulics · structure | 10 | `Vv` `Vsas` `Vsdh` `Rout` `E1` `Vplast` `Wpv` `AQ` `Cart` `Autoreg` |
| White matter · glia | 5 | `WMint` `WMperm` `Myel` `Astro` `Micro` |
| Amyloid · tau · CSF biomarkers | 7 | `Ab_isf` `Ab_plq` `Tau_isf` `A_ab` `A_pt` `A_nfl` `A_lrg` |
| Drug PK | 11 | acetazolamide (gut · central · **saturable erythrocyte binding**), loop diuretic, solifenacin, donepezil, melatonin |
| Systemic safety | 3 | `HCO3` `Kser` `AZesc` (pharmacodynamic escape) |
| Device | 2 | `Occl` `Vdr` |
| Clinical endpoints | 7 | `Gslow` `Gfast` `Cslow` `Cfast` `Urin` `Headx` `SDHhaz` |
| Acute perturbation | 1 | `dVac` (lumbar puncture volume deficit) |

**Mean intracranial pressure is not a compartment.** It is obtained **in closed
form** at every derivative evaluation from the CSF flow balance (which is linear in
P once the valve state is determined). Craniospinal mechanics equilibrate on a
timescale of minutes while the disease progresses over years, so solving the fast
block algebraically is what keeps this system a **non-stiff 45-state system**
rather than a stiff 50-state one.

**Posture is a weight, not an oscillation.** Shunt hydraulics differ by tens of
cmH₂O between supine and upright. Instead of resolving the 24-hour cycle, the
hydraulics block is evaluated **twice** per derivative call (supine and upright) and
averaged with the weight `f_up` (the fraction of the day spent upright). The
hydrostatic physics is kept exact while the integrator stays on daily steps.

---

## 2. Siphoning — the complication rate is a hydrostatics problem, not a biology problem

A shunt is **a pressure-controlled resistance connected in parallel with the broken
absorption pathway**. When upright, the fluid column from ventricle to peritoneum
**adds directly** to the driving pressure. No biological parameter in this model
touches that term.

| Configuration | ICP supine | ICP upright | Daily mean ICP | Shunt flow [mL/day] | Relative to production |
|---|---|---|---|---|---|
| untreated | 13.15 | 5.15 | 8.35 | 0 | 0.00 |
| valve 10 cmH₂O, no protection | 11.28 | **−15.65** | **−4.88** | **1002** | **1.99×** |
| valve 10 + membrane ASD 60% | 11.28 | 1.90 | 5.65 | 205 | 0.41× |
| valve 10 + **gravitational 30 cmH₂O** | 11.28 | 3.85 | 6.82 | 116 | 0.23× |
| valve 20 cmH₂O, no protection | 13.15 | −9.15 | −0.23 | 650 | 1.29× |
| **ETV** (no distal column at all) | 12.75 | 4.75 | 7.95 | 0 | 0.00 |

Decomposed by posture: an unprotected 10 cmH₂O valve **drains 142 mL/day (28% of
production) when supine and holds ICP at 11.28 mmHg — a device that works well.**
But **when upright it drains 1576 mL/day, 3.13 times production, and pulls ICP down
to −15.65 mmHg.** The valve is not mis-specified. It is being asked to deal with
**a 45 cmH₂O hydrostatic column that it cannot sense.**

A gravitational assist is **an opening-pressure term that switches on with
posture**, and it contains nothing one would call pharmacology. And that is the
whole of the solution.

---

## 3. ICP after shunting is a property of the device, not of the disease

The flow balance with the valve open is

```
If = (P − Pss)/Rout + (P + h − Popen − Pdist)/Rsh
P  = [ If + Pss/Rout + (Popen + Pdist − h)/Rsh ] / (1/Rout + 1/Rsh)
```

Since `Rsh ≪ Rout`, every Rout term receives a low weight. Holding the valve fixed
and sweeping only disease severity:

| Rout | ICP without shunt | ICP with shunt | dP/dRout (none) | dP/dRout (shunt) |
|---|---|---|---|---|
| 9.0 | 9.65 | 9.65 | — | — |
| 13.0 | 11.05 | 11.03 | 0.3500 | 0.3457 |
| 19.0 | 13.15 | 11.28 | 0.3500 | **0.0405** |
| 23.0 | 14.55 | 11.37 | 0.3500 | 0.0246 |
| 26.0 | 15.60 | 11.43 | 0.3500 | 0.0186 |

Without a shunt the slope is exactly `If` = 0.350 mmHg per unit Rout **at every
severity**. Over the range where the valve is open (Rout ≥ 19) the mean slope is
**0.0279**, so the device absorbs **92.0% of the residual influence of the
disease.** That is, a shunt **replaces the patient's outflow resistance with the
manufacturer's outflow resistance.**

And look at the first two rows. **At Rout 9–13 the pressure is identical with and
without a shunt — because the valve does not open at all when supine.** A 10 cmH₂O
valve plus 5 cmH₂O of intraperitoneal pressure is 11.0 mmHg, which is above normal
supine ICP. **The same hardware is non-functional in a mild patient and dominant in
a severe one.** This is the model deriving why iNPH shunts have to use low opening
pressures, and at the same time it means that a postoperative ICP measurement tells
you more about the valve than about the patient.

---

## 4. The valve titration map — what sets the width of the window is posture

24 months, `f_up = 0.60` (14.4 hours a day upright). Utility = gait gain −
1.1 × subdural incidence − 0.03 × headache index.

**Unprotected**

| Popen [cmH₂O] | Daily mean ICP | Upright ICP | ΔGait [m/s] | Hygroma [mL] | Subdural [%] | Utility |
|---|---|---|---|---|---|---|
| 2 | 1.58 | −9.14 | **−0.336** | 91.8 | 53.6 | −1.229 |
| 10 | 1.76 | −8.95 | −0.093 | 51.9 | 28.2 | −0.702 |
| 20 | 1.93 | −7.09 | +0.163 | 14.7 | 0.0 | −0.080 |
| 28 | 3.20 | −3.75 | +0.187 | 0.0 | 0.0 | +0.045 |
| **36** | 6.29 | 1.39 | **+0.130** | 0.0 | 0.0 | **+0.130** |

**Gravitational assist 30 cmH₂O**

| Popen [cmH₂O] | Daily mean ICP | Upright ICP | ΔGait [m/s] | Hygroma [mL] | Subdural [%] | Utility |
|---|---|---|---|---|---|---|
| 2 | 1.99 | −0.98 | +0.265 | 1.9 | 0.0 | +0.205 |
| **4** | 3.06 | 0.09 | **+0.244** | 0.0 | 0.0 | **+0.217** |
| 8 | 5.64 | 2.66 | +0.142 | 0.0 | 0.0 | +0.142 |
| 12 | 8.19 | 5.22 | −0.035 | 0.0 | 0.0 | −0.035 |
| 14–36 | 8.83 | 5.63 | −0.060 | 0.0 | 0.0 | −0.060 |

**A gravitational assist does not merely lower the complication rate at a fixed
setting — it moves the optimum itself.** Without protection the best setting is
36 cmH₂O (+0.130 m/s); with a gravitational unit the optimum drops to 4 cmH₂O and
delivers +0.244 m/s. **What the device buys is 32 cmH₂O of usable titration range,
and the value of that range is +0.115 m/s of gait.**

Look at what happens at either end. The model was not instructed to do either of
these things.
- **Set it too high** and the valve does not open at all when supine, and the
  patient simply follows the natural history (the flat tail of the gravitational arm
  reproduces the 24-month deterioration of S1 exactly).
- **Set it too low without gravitational protection** and gait becomes **worse than
  not operating at all.** Because the mass effect of the subdural hygroma and the
  low-pressure headache cancel out a hydraulic improvement that genuinely happened.

The unprotected arm is **pushed into a high, weakly effective corner** because it
has no other safe option. The gravitational arm is not.

---

## 5. Why the tap test misses responders — and a physics error found here

The textbook treatment decays the pressure effect of a tap with
`τ = Rout·C ≈ 6 min`. That is **a linearisation of the flow balance**, and for a
40 mL perturbation it is off by an order of magnitude.

```
linearised: recovery rate = dVac/(Rout·C) = 40 mL / 5.6 min ≈ 7 mL/min
actual:     below sinus pressure, absorption = 0 and shunt = 0
            → recovery is limited to the choroid plexus secretion rate of 0.35 mL/min
            → refilling 40 mL takes at least 114 min
```

While it carried this error, the model produced the conclusion — contrary to the
entire literature — that **"the tap test can never be positive"**. Fixing it with
the exact balance:

| After tap [h] | ICP supine | AMP | Gait | ΔGait |
|---|---|---|---|---|
| 0.00 | −2.00 | 0.48 | 0.614 | +0.000 |
| 0.50 | −1.98 | 0.48 | 0.654 | +0.040 |
| 1.00 | −1.81 | 0.48 | 0.669 | +0.054 |
| 3.00 | 13.11 | 4.80 | 0.670 | +0.056 |
| 12.00 | 13.12 | 4.80 | 0.651 | +0.037 |
| 24.00 | 13.12 | 4.80 | 0.636 | +0.021 |
| 48.00 | 13.12 | 4.80 | 0.622 | +0.007 |

The pressure stays at −2.0 mmHg for the first 114 minutes (the tap pushed the
pressure onto the flat foot of the P–V curve, the region where compliance is
saturated). Then at 3 hours it **springs back** to 13.1 mmHg — almost all of the
recovery happens in the last few mL, because the curve is exponential. **The gait
peak is +0.060 m/s at 1.6 hours and is still +0.021 m/s at 24 hours. Mean ICP
returns to normal long before the patient does — so the effect variable cannot be
mean pressure.**

| | Gain [m/s] | Relative to shunt |
|---|---|---|
| shunt at 12 months | **+0.192** | 100% |
| tap test peak | +0.060 | 31% |
| 72-hour continuous drainage (ELD) | +0.031 | 16% |

The **time window** over which each perturbation holds the patient above the
decision threshold: at a threshold of +0.02 m/s the tap gives 25.1 hours and the
ELD 81.5 hours. At a threshold of +0.05 the tap gives 4.7 hours and the ELD 0.

**The central claim:** sensitivity here is not a statement about biology but **a
statement about threshold crossing**. At a 0.10 m/s threshold the tap test reads
**negative**, and yet that patient gains +0.192 m/s from surgery. In a model with
no responder/non-responder covariate at all, **a false negative is manufactured
purely by the measurement.**

### A discrepancy reported rather than calibrated away (a miss)

This model rates a single 40 mL tap (+0.060) **above 72 hours of continuous
drainage (+0.031)**, which is **the reverse** of the sensitivity ordering in the
literature (tap 26–61%, ELD 50–100%). The reason is structural and worth stating.
The model's fast arm responds to the **depth** of the reduction in pulse amplitude,
and a 40 mL bolus (AMP 4.81 → 0.48) is far deeper than 10 mL/h of continuous
drainage. Meanwhile the slow arm has `tau_gait = 34 days`, so it converts 72 hours
of sustained decompression into almost no measurable gait change. In other words,
**the model reproduces depth and fails to reproduce duration.** The real advantage
of continuous drainage is probably duration — three days of continuous improvement
gives far more opportunities to cross a threshold. Fixing this would need either a
faster tissue-level arm or an explicit repeated-measures model, and neither can be
calibrated with the data available, so **neither was added.**

---

## 6. Early versus late surgery — the reversible indices converge

One patient, one set of hardware (4 cmH₂O + gravitational 30), operated at 6, 12,
24, 36 and 60 months after symptom onset. **Every row is read at 36 months after
surgery, so the post-operative exposure is identical.**

| Delay [months] | Evans | AMP | Wpv | WMint | **WMperm** | **Gait** | MMSE |
|---|---|---|---|---|---|---|---|
| 6 | 0.396 | 2.13 | 0.54 | 0.818 | **0.144** | **0.912** | 27.8 |
| 12 | 0.400 | 2.13 | 0.54 | 0.803 | 0.159 | 0.894 | 27.7 |
| 24 | 0.407 | 2.14 | 0.55 | 0.771 | 0.191 | 0.856 | 27.4 |
| 36 | 0.415 | 2.14 | 0.55 | 0.739 | 0.224 | 0.816 | 27.1 |
| 60 | 0.431 | 2.14 | 0.56 | 0.673 | **0.292** | **0.736** | 26.4 |

Between a 6-month and a 60-month delay, **the pulse amplitude ends up within
0.01 mmHg (0.4%) and the interstitial water within 0.03 (4.8%) — the hydraulics
recover either way.** What does not converge is the irreversible pool: **WMperm
0.144 vs 0.292 (2.03-fold)**, and with it **gait 0.912 vs 0.736 m/s — 0.176 m/s as
the price of waiting** — and MMSE 27.8 vs 26.4.

**The indices the clinician can see on the postoperative images are precisely the
indices that forget the delay.**

---

## 7. Acetazolamide — a ceiling set by arithmetic, a dose that buys only acidosis

Because `ICP = Pss + If·Rout`, the maximum pressure that can be removed by
suppressing CSF production is `Emax · If · Rout = 0.50 × 0.35 × 19.0 = 3.32 mmHg`,
and **`Pss` can never be touched.** The sustainable ceiling after escape is
**0.73 mmHg**. A shunt does not act on this axis — it replaces `Rout` itself.

| Dose (BID) | Free concentration [mg/L] | E_az day 3 [%] | E_az day 90 [%] | ICP supine | ΔGait day 90 | HCO₃ | K |
|---|---|---|---|---|---|---|---|
| none | 0.000 | 0.0 | 0.0 | 13.21 | +0.001 | 24.0 | 4.20 |
| 125 mg | 0.187 | 38.9 | 12.7 | 12.29 | +0.019 | 22.1 | 4.00 |
| 250 mg | 0.358 | 41.1 | 12.1 | 12.37 | +0.018 | **20.9** | 3.88 |
| 500 mg | 0.701 | 42.5 | 11.7 | 12.41 | +0.017 | **19.5** | 3.72 |
| 1000 mg | 1.391 | 43.2 | 11.5 | 12.43 | +0.017 | **18.0** | 3.57 |

**Three things happen at once and only one of them is dose-dependent.**

1. **Choroid plexus occupancy saturates.** The free-concentration EC50 is
   0.028 mg/L, so 125 mg BID is already high on the Hill curve. Four further
   doublings from there barely increase the day-3 effect (38.9% → 43.2%).
2. **The effect escapes.** Compare the day-3 and day-90 columns. The acute
   reduction in CSF production has mostly gone by three months (`tau_esc` = 18 days,
   `f_esc` = 0.78). **Acetazolamide behaves like a temporary tap test, not like a
   shunt.**
3. **The acidosis does not escape, and it is dose-dependent.** Because the systemic
   acid-base EC50 (1.20 mg/L) is about 43 times the choroid plexus EC50.

Put together: above 125–250 mg BID **the therapeutic index deteriorates
monotonically with essentially no hydraulic gain.** It is a ceiling on one axis and
a slope on another; it is not a titration. As the dose rises, **the escape grows
too**, so the day-90 effect actually *falls* slightly.

### Non-linear PK — saturable erythrocyte carbonic anhydrase binding

| Dose (BID) | Central [mg] | Erythrocyte-bound [mg] | Erythrocyte fraction | Total concentration [mg/L] | Free concentration [mg/L] |
|---|---|---|---|---|---|
| 62.5 mg | 28.15 | 88.60 | **0.759** | 2.011 | 0.1005 |
| 125 mg | 52.46 | 100.13 | 0.656 | 3.747 | 0.1874 |
| 250 mg | 100.20 | 107.97 | 0.519 | 7.157 | 0.3579 |
| 500 mg | 196.22 | 113.15 | 0.366 | 14.016 | 0.7008 |
| 1000 mg | 389.55 | 116.30 | **0.230** | 27.825 | 1.3912 |

The drug has to **fill a deep erythrocyte pool first.** At low doses 76% of it is
in there, and once the pool saturates that falls to 23% — **the reason the apparent
half-life lengthens with dose lies in distribution, not elimination.**

---

## 8. CSF biomarkers are measured in a space that is diluted and drains badly

A CSF concentration is the influx flux divided by the clearance.

```
c = J(brain→CSF) / [ (If + Qsh)/Vcsf + kdeg ]
```

In iNPH **the numerator (glymphatic efflux) falls, and the turnover term in the
denominator falls too** (because Vcsf has grown). Since the two terms push c in
**opposite directions**, **the fact that Aβ42 and p-tau are low together is a
statement about which term won, and not direct evidence of Alzheimer pathology.**

| | Healthy 75-year-old | iNPH | Ratio |
|---|---|---|---|
| CSF volume [mL] | 135.0 | 157.0 | 1.16× |
| CSF turnover [/day] | 3.733 | 3.210 | 0.86× |
| AQP4 polarisation | 0.850 | 0.450 | 0.53× |
| Cortical plaque burden | 0.116 | 0.228 | 1.97× |
| **CSF Aβ42** [pg/mL] | 798 | **457** | **0.57×** |
| **CSF p-tau** [pg/mL] | 45.1 | **31.3** | **0.69×** |
| CSF NfL [pg/mL] | 814 | 1418 | 1.74× |
| CSF LRG [a.u.] | 1.00 | 2.73 | 2.73× |

**Aβ42 is low and p-tau is low as well** — exactly the opposite combination to
Alzheimer's disease, where p-tau *rises* — and this is the signature of the iNPH
CSF profile. The model **reproduces this combination as the arithmetic of dilution
and reduced efflux, and to do so it put in no iNPH-specific term at all.**

### The post-shunt trajectory — flagged as a prediction

| Month | CSF turnover | AQP4 polarisation | Shunt flow | Aβ42 | p-tau | NfL |
|---|---|---|---|---|---|---|
| 0 | 5.831 | 0.450 | 412 | 457 | 31.3 | 1387 |
| 1 | 5.911 | 0.458 | 411 | **274** | 18.6 | 685 |
| 3 | 5.993 | 0.490 | 411 | 291 | 19.6 | 659 |
| 6 | 6.049 | 0.538 | 411 | 315 | 21.0 | 629 |
| 12 | 6.077 | 0.615 | 410 | 352 | 22.9 | 583 |
| 24 | 6.067 | 0.698 | 408 | **387** | 24.2 | 540 |

The trajectory is **not monotonic.** The moment the shunt goes in, clearance rises
first and drops Aβ42 from 457 to 274 (dilution); then over the following months, as
glymphatic efflux recovers, it climbs back to 387. The magnitudes of the two effects
are **within 1.18-fold** of each other. This model therefore **does not claim that
the direction is robust.** It claims that the sign is determined by a single ratio
that has never been measured (`kflux_ab × ΔAQ` versus `ΔQsh/Vcsf`), and it claims
that **a trial result reported as "no significant change" is entirely compatible
with large changes in both underlying fluxes.**

---

## 9. Sensitivity — the model refuted the premise it was built on

Elasticity of the 24-month gait gain, `d ln(ΔGait) / d ln(parameter)`, ±10% central
difference, at the titrated shunt setting (4 + gravitational 30). Baseline gain
+0.2443 m/s.

| Rank | Parameter | Value | Elasticity | What it is |
|---|---|---|---|---|
| 1 | **`kW_out`** | 1.300 | **+2.273** | drainage rate of periventricular interstitial water |
| 2 | **`Hcol_cm`** | 45 | **+1.653** | ventricle-to-peritoneum hydrostatic column = **the patient's height** |
| 3 | **`kAQ_rec`** | 0.008 | **+0.978** | rate of recovery of AQP4 polarisation |
| 4 | `CBF_crit` | 18 | −0.379 | perfusion threshold for white matter injury |
| 5 | `k_rep` | 0.0165 | +0.256 | rate of white matter functional recovery |
| 6 | `Rout_init` | 19 | −0.138 | the patient's CSF outflow resistance |
| 7 | `Rsh` | 2.5 | −0.116 | shunt hydraulic resistance |
| 8 | `k_perm` | 1.0e-4 | −0.096 | rate of irreversible conversion |
| … | | | | |
| **14** | **`kappa_tm`** | 0.060 | **−0.049** | fraction of the ICP pulsation transmitted across the mantle |

**The structural centre of this model is the pulsatile transmantle gradient, so
`kappa_tm` was expected to dominate. It did not — it is 10th of 14.** Both the
expectation and its refutation have been left in the report.

The reason `kappa_tm` drops out is a saturation the model was not designed to
create. The transmantle driving pressure is comfortably above the creep yield
threshold `dP_yield` at every severity tested, so perturbing `kappa_tm` moves
ventricular creep but does not move **what actually limits white matter recovery —
the periventricular water and the perfusion it costs.** **`kappa_tm` determines the
mechanism; `kW_out` and `kAQ_rec` determine the answer.**

Two implications worth stating plainly.

1. Measuring `kappa_tm` better sharpens the story but changes the predictions
   hardly at all. Conversely, **measuring the interstitial water drainage axis
   changes the predictions a great deal** — that is, the uncertainty this model
   cites most is not the uncertainty that matters most.
2. **`Hcol_cm` is second — above every pharmacological and pulsatile parameter.**
   In a model containing five drugs, **the patient's height matters more to gait at
   two years than any of the drugs do.** Because height is what determines which
   siphon the valve has to be chosen against. This is the same conclusion §2 reached
   from hydraulics alone, but here it was reached by a completely independent route.

> For the full elasticity table see §11 of [`inph_model_report.txt`](../../../normal-pressure-hydrocephalus/inph_model_report.txt).

---

## 10. Model calibration — what was fitted, what is wrong, and what is a guess

### Parameters fitted to literature values (only nine)

| Parameter | Value | Target |
|---|---|---|
| `Rout_init` | 19 | iNPH infusion-test outflow resistance 15–25 vs normal < 13 mmHg/(mL/min) |
| `E1_init` / `Vp_init` | 0.222 / 1.35 | iNPH AMP ≈ 4–5, normal ≈ 2 mmHg (response-prediction threshold 4 mmHg) |
| `kE1_rec` | 0.010 | 35–40% reduction in pulse amplitude over months after a working shunt |
| `a_fast` | 0.075 | maximum tap-test gait gain inside the positive-decision range (5–10%) |
| `k_creep` | 26 | untreated Evans index drifting about 0.01 per year (not running away) |
| `k_sdh` | 1.55 | hygroma of an unprotected valve in the tens-of-mL range where it becomes symptomatic |
| `kflux_ab` / `kflux_tau` | 98000 / 26400 | healthy control CSF Aβ42 ≈ 800, p-tau ≈ 45 pg/mL |

### Discrepancies reported rather than hidden by calibration

1. **Response rate.** This model makes almost every patient with a raised `Rout` a
   responder, because the only non-response mechanisms are `WMperm` and the
   APOE/plaque terms. Real cohorts improve in 60–80% at one year, and
   placebo-controlled (valve-closed) trials are more sceptical still.
   → **It overestimates the value of shunting in an unselected population.** S6 and
   S21 must be read as the **best case**.
2. **Ventricular size.** It reduces the Evans index by 0.01–0.02 after a successful
   shunt. Much of the literature reports no significant change, so it **still couples
   size too tightly to pressure.**
3. **Acetazolamide.** It gives 39–43% acutely and 12% chronically. The clinical
   evidence is weaker than that. S12 and S13 are **upper bounds**.
4. **The tap versus continuous drainage ordering is inverted.** Quantified in §5 and
   not calibrated away.
5. **The frequency band of the pulsatility.** Recent literature reports that brain
   deformation correlates better with vasomotion, which is slower than the cardiac
   pulse. This model hangs `kappa_tm` on the cardiac component (`AMP`). **The model
   does not adjudicate which is right, and there are no data with which to
   adjudicate.**
6. **A single CSF compartment.** There are reports that ventricular and lumbar CSF
   biomarkers differ, but this model reads the tap test (lumbar) and the shunt
   (ventricular) at the same concentration.
7. **The healthy control MMSE** sits at 28–29 rather than 29–30, because the
   amyloid block gives a 75-year-old a non-zero plaque burden.

### Uncalibrated structural guesses

| Parameter | Value | Status |
|---|---|---|
| `kW_out` | 1.300 | **the parameter the answer actually hangs on (1st in §9), and it is an order-of-magnitude estimate** |
| `kAQ_rec` | 0.008 | rate of AQP4 repolarisation after decompression. No human measurement (3rd in §9) |
| `kappa_tm` | 0.060 | never measured in a human ventricle. **Expected to be the most influential and came 14th in §9** |
| `DESH` | 1.55 | morphology multiplier. Pure estimate |
| `k_perm` | 1.0e-4 | **this one parameter determines the result of §6 (early versus late)** |
| `Hcol_cm` | 45 | varies with the patient's height. The dominant term of §2 |
| `k_haz` | 1.5e-5 | **so the subdural incidences (%) in §4 are an ordinal scale and not calibrated incidences** |
| `etv_dRout` | 0.06 | **a structural argument** that the resistance lies downstream of the stoma, not a measurement |
| `az_tau_esc` / `az_f_esc` | 18 days / 0.78 | the mechanism is accepted but has never been measured at this rate in iNPH |
| `tau_on` / `tau_off` | 29 min / 21.6 h | the asymmetry of the fast arm. A modelling choice, but **falsifiable** — it predicts a fast onset and slow decay of tap-test improvement, which is the pattern reported clinically |

### Defects found and fixed by the model's own checks during development

| # | Defect | Consequence |
|---|---|---|
| a | The acute volume perturbation was handled as a **linear offset `dVac/C` plus 3 fixed-point iterations** | at a 40 mL tap the iteration **failed to converge**, oscillating between −0.2 and −22.5 mmHg before settling on the −30 clamp. Replaced by the Marmarou **integral** form `P = P0 + (Pbase−P0)·exp(E1·dV)` (explicit and monotonic) |
| b | **Tap-test refilling was 20-fold too fast — the most serious defect found.** `τ = Rout·C` (the linearisation) implies 7 mL/min whereas the choroid plexus secretion rate is 0.35 mL/min | pressure returned to normal in 30 minutes → the model produced the wrong conclusion that **"the tap test can never be positive"**. With the exact balance, absorption is zero below sinus pressure, so recovery is **production-limited** and 40 mL takes 114 min |
| c | A symmetric first-order filter (τ = 13 h) integrates almost nothing from a 2-hour perturbation | the fast arm was made **asymmetric** (charging 29 min / discharging 21.6 h) |
| d | Ventricular recoil had **no plasticity floor** | shunted patients' ventricles returned to a normal 25 mL within a year — the opposite of the entire literature. `Vplast` added |
| e | **There was no valve-closed branch** | at high opening pressures it computed a **negative shunt flow** (reflux from peritoneum to ventricle), and the high-pressure arm of §4 looked therapeutic |
| f | The mean transmantle gradient depended only on `Rout` | it **did not respond to a shunt at all**, so periventricular water was never drained. Scaled by `par_frac`, the fraction still absorbed through the high-resistance parenchymal route |
| g | The **production rates for reactive astrocytosis and microglial activation were an order of magnitude larger than their sinks** | the untreated equilibrium was `Astro` = 20.7 on a 0–3 scale → quietly drove AQP4 polarisation to zero |
| h | The AQP4 loss term acted even on **baseline** astrocytosis | **glymphatic function in the healthy control drifted over five years** and gait deteriorated from 1.05 to 0.90 as a result. Fixed to act only on the excess |
| i | `Rout` progression was unconditional | **the healthy control's outflow resistance drifted to `Rout_max`.** Gated on the abnormality already being present |
| j | Acetazolamide elimination was **applied to the free pool but scaled on total drug** | half-life of 70 days, concentrations 100-fold too high. `CL_u = CL_total/fu` |
| k | Acetazolamide had **no pharmacodynamic escape** | **the carbonic anhydrase inhibitor came out better than a shunt** — the exact opposite of the clinical record |

---

## 11. The 21 scenarios

At 24 months (S10 and S11 at 48 months, S18 at 36 months). ΔGait is the change
relative to each scenario's own baseline value.

| # | Scenario | Gait | ΔGait | MMSE | GS | AMP | Hygroma | Subdural % |
|---|---|---|---|---|---|---|---|---|
| S1 | untreated iNPH | 0.554 | −0.060 | 24.4 | 7 | 5.60 | 0.0 | 0.0 |
| S2 | healthy 75-year-old control | 1.082 | +0.031 | 29.1 | 0 | 2.03 | 0.0 | 0.0 |
| S3 | tap test 40 mL | 0.619 | +0.004 | 25.0 | 7 | 4.80 | 0.6 | 0.0 |
| S4 | continuous lumbar drainage 10 mL/h × 72 h | 0.646 | +0.032 | 25.2 | 7 | 3.76 | 0.0 | 0.0 |
| S5 | fixed valve 10, no siphon protection | 0.535 | **−0.093** | 24.2 | 6 | 5.32 | **51.9** | **28.2** |
| S6 | **valve 4 + gravitational 30 (titrated setting)** | **0.871** | **+0.244** | 27.2 | 3 | 2.26 | 0.0 | 0.0 |
| S6b | valve 10 + gravitational 30 (underdrainage) | 0.669 | +0.051 | 25.4 | 6 | 3.96 | 0.0 | 0.0 |
| S7 | programmable, high pressure 16 cmH₂O | 0.554 | −0.060 | 24.4 | 7 | 5.60 | 0.0 | 0.0 |
| S8 | low pressure 4 cmH₂O, **no gravitational unit** | 0.353 | **−0.275** | 22.7 | 8 | 5.27 | **81.8** | **48.3** |
| S9 | stepwise downward titration 16 → 8 | 0.776 | +0.161 | 26.4 | 4 | 2.87 | 0.0 | 0.0 |
| S10 | early shunt (6-month delay), 48 months | 0.914 | **+0.299** | 27.9 | 2 | 2.14 | 0.0 | 0.0 |
| S11 | late shunt (36-month delay), 48 months | 0.758 | **+0.143** | 26.2 | 4 | 2.21 | 0.0 | 0.0 |
| S12 | acetazolamide 250 mg BID | 0.592 | −0.022 | 24.7 | 7 | 4.60 | 0.0 | 0.0 |
| S13 | acetazolamide 500 mg BID | 0.590 | −0.024 | 24.7 | 7 | 4.64 | 0.0 | 0.0 |
| S14 | shunt + acetazolamide 250 BID | 0.873 | +0.246 | 27.2 | 3 | 2.23 | 0.0 | 0.0 |
| S15 | shunt + comorbid AD (APOE ε4) | 0.853 | +0.226 | 26.6 | 3 | 2.26 | 0.0 | 0.0 |
| S16 | unprotected valve + very active (`f_up` 0.80) | 0.300 | **−0.330** | 22.5 | 8 | 6.33 | **93.3** | **54.3** |
| S17 | shunt + solifenacin 5 mg | 0.871 | +0.244 | **25.2** | 3 | 2.26 | 0.0 | 0.0 |
| S18 | shunt occlusion at 18 months + revision | 0.868 | +0.241 | 27.3 | 3 | 2.21 | 0.0 | 0.0 |
| S19 | shunt + melatonin 2 mg | 0.959 | **+0.332** | 28.3 | 1 | 2.14 | 0.0 | 0.0 |
| S20 | ETV in communicating iNPH | 0.563 | −0.053 | 24.4 | 7 | 5.44 | 0.0 | 0.0 |
| S21 | **the full bundle** (gravitational 4, melatonin, donepezil) | **0.959** | **+0.332** | **29.9** | **1** | 2.14 | 0.0 | 0.0 |

A few contrasts worth reading.

- **S5 vs S8 vs S6.** In the same situation without a gravitational unit, lowering
  the opening pressure from 10 to 4 grows the hygroma from 51.9 to 81.8 mL and makes
  gait worse still, −0.093 → −0.275. **The same 4 cmH₂O used together with a
  gravitational unit gives +0.244.** The setting itself is neither good nor bad —
  only its combination with siphon protection has meaning.
- **S7 = S1.** A 16 cmH₂O valve is identical to no treatment to three decimal
  places, because the valve does not open when supine. There was an operation and
  there was no treatment.
- **S20 ≈ S1.** In communicating iNPH, ETV at −0.053 is essentially
  indistinguishable from no treatment (−0.060). Because the resistance is downstream
  of the stoma — but in exchange, **there is no distal catheter, so the siphoning risk
  is also zero.**
- **The value of S17 is not in the gait but in the MMSE.** Solifenacin improves the
  urinary indices but its central antimuscarinic action **costs 2.0 MMSE points**,
  27.2 → 25.2. It does not touch gait or the hydraulics at all.
- **S19 vs S21.** Melatonin raises gait from +0.244 to +0.332 through the glymphatic
  axis. Adding donepezil on top leaves gait unchanged (0.959) and raises only the
  MMSE, 28.3 → 29.9 — **the two drugs are on different axes and neither replaces the
  valve.**
- **S12 and S13 are both still negative.** Acetazolamide is better than no
  treatment (−0.022 vs −0.060) but does not stop the progression of the disease, and
  500 mg BID is **worse** than 250 mg BID (−0.024 vs −0.022) — because the escape is
  larger.

---

## 12. How to reproduce

```bash
# the reference implementation — no dependencies, prints every number in this README
python3 inph_reference_model.py             # the full report (about 20 min)
python3 inph_reference_model.py --brief     # the key numbers only

# re-render the mechanistic map
dot -Tsvg inph_qsp_model.dot -o inph_qsp_model.svg
dot -Tpng -Gdpi=150 inph_qsp_model.dot -o inph_qsp_model.png
```

```r
# the mrgsolve model (R)
source("inph_mrgsolve_model.R")
all <- run_all()
print(summarise_scenarios(all), n = 25, width = Inf)
print(titration_map())     # the titration map of §4
print(delay_sweep())       # the two clocks of §6

# Shiny dashboard
shiny::runApp("inph_shiny_app.R")
```

**The Python implementation is the verified implementation.** The R code and the
Python code use identical parameter names and identical equations, but R is not
installed in this container, so **the R/Shiny code has not been verified by
execution.** The two intended differences between the implementations are stated in
the header of `inph_mrgsolve_model.R` (the integration scheme for `dVac` and the
handling of the upper bound on `E1`).

---

## ⚠️ Disclaimer

This model is a **QSP model for educational and research purposes**. It was built
on the published literature but has not been independently verified or certified,
and **must not be used for actual clinical decision-making, prescribing or
regulatory submission.** In particular, **the valve titration map of §4 contains an
uncalibrated hazard coefficient (`k_haz`) and therefore cannot form the basis of an
actual valve setting.**
