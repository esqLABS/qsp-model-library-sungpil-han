# Severe Traumatic Brain Injury — QSP Model
### Severe TBI · intracranial pressure, perfusion and metabolism · 47-ODE QSP model

<a href="../../../traumatic-brain-injury/tbi_qsp_model.svg"><img src="../../../traumatic-brain-injury/tbi_qsp_model.png" width="640" alt="sTBI QSP mechanistic map"></a>

---

## Why this model is shaped differently

Almost every model in this repository is a **cascade**. An antigen or a mutation
moves a cytokine, the cytokine moves an effector, the effector moves a clinical
score, and a drug cuts one arrow.

Severe traumatic brain injury is not that kind of structure. Severe TBI is a
**constraint problem**, and the constraint is that the skull is a closed box.

```
    Monro-Kellie:   V_blood + V_CSF + V_brain + V_lesion = constant
```

The variable we treat, intracranial pressure (ICP), is therefore **not a state
variable with a biology of its own.** It is the **residual** of a volume balance,
read off an exponential pressure-volume curve. Five conclusions follow immediately
from this single structural fact, and this model is built to **compute** them rather
than to assert them.

The governing equation of this model is one line.

```
dICP/dt · [ PVI/(ICP·ln10) + Ca − dVv/dICP ]
      = dCa/dt·(MAP−ICP) + Ca·dMAP/dt + (If − Io − Q_evd)
        + Jw_normal + Jw_injured + dV_hem/dt
```

Everything a drug can do is to change one term on the right-hand side.

---

## Files

| File | Contents |
|---|---|
| [`tbi_qsp_model.dot`](../../../traumatic-brain-injury/tbi_qsp_model.dot) | Mechanistic map — **236 nodes · 335 edges · 20 clusters** |
| [`tbi_qsp_model.svg`](../../../traumatic-brain-injury/tbi_qsp_model.svg) / [`.png`](../../../traumatic-brain-injury/tbi_qsp_model.png) | Rendered output (150 dpi) |
| [`tbi_reference_model.py`](../../../traumatic-brain-injury/tbi_reference_model.py) | **Reference implementation.** Pure standard library, 47 ODEs, fixed-step RK4 |
| [`tbi_reference_output.txt`](../../../traumatic-brain-injury/tbi_reference_output.txt) | The complete stdout of the file above — **the sole source of every number below** |
| [`tbi_mrgsolve_model.R`](../../../traumatic-brain-injury/tbi_mrgsolve_model.R) | mrgsolve model (a 1:1 port of the Python implementation) + 12 scenarios |
| [`tbi_shiny_app.R`](../../../traumatic-brain-injury/tbi_shiny_app.R) | Shiny dashboard, 11 tabs |
| [`tbi_references.md`](../../../traumatic-brain-injury/tbi_references.md) | 132 references + parameter-provenance table + **a table of where the model is wrong** |

```bash
python3 tbi_reference_model.py           # regenerates every number below
dot -Tsvg tbi_qsp_model.dot -o tbi_qsp_model.svg
Rscript -e 'shiny::runApp("tbi_shiny_app.R")'
```

> **Provenance rule.** Every number in this README was obtained by running
> `tbi_reference_model.py` and is contained verbatim in
> `tbi_reference_output.txt`. Not one number was written from memory. Where the
> Python implementation and the mrgsolve implementation disagree, **the Python side
> is the specification.**

---

## Verification comes first (R0: the healthy brain must sit still)

The reference point (MAP 88, ICP 10, PaCO2 40, 37 °C, no drugs) was not hand-tuned
but **back-solved**: CSF formation = absorption, net capillary water movement = 0,
renal sodium and volume balance = 0, autoregulatory set-point = the actual flow. As
a result the healthy brain does not move when integrated for 12 hours.

| Quantity | t=0 | t=12 h | 12-hour drift |
|---|---:|---:|---:|
| ICP (mmHg) | 10.0000 | 9.9962 | **−0.0038** |
| CBF (mL/100g/min) | 49.5289 | 49.6205 | +0.0916 |
| CMRO2 (mL O₂/100g/min) | 3.2990 | 3.3170 | +0.0179 |
| PbtO2 (mmHg) | 22.6663 | 22.4697 | −0.1966 |
| SjvO2 (%) | 58.9125 | 58.7726 | −0.1399 |
| Lactate/pyruvate ratio | 15.7530 | 15.1971 | −0.5559 |
| Plasma Na (mmol/L) | 140.0000 | 140.0000 | **0.0000** |
| Plasma osmolality (mOsm/kg) | 293.0000 | 293.0000 | **0.0000** |
| Excess brain water (mL) | 0.0000 | 0.3138 | +0.3138 |
| ECF glutamate (µM) | 2.0000 | 2.0040 | +0.0040 |
| ECF K⁺ (mmol/L) | 3.0000 | 3.0000 | +0.00002 |
| Infarct fraction | 0.0000 | 0.0000 | **0.0000** |

The small residual drift (brain water +0.31 mL, CMRO2 +0.018) comes from the
temperature rising very slowly, and it is reported as it is rather than hidden.

**Step-size convergence:** from dt = 0.0025 → 0.02 min the ICP, infarct fraction,
brain water, CSF volume and PbtO2 at the 24-hour time point **agree to five decimal
places**. All results were computed with dt = 0.02.

**Segmented integration = continuous integration:** computations that take the state
out and put it back in mid-course (the reserve probe, treatment branches) are handed
the absolute time `t0`. Omit that and terms with an explicit time dependence (the
haematoma expansion window, the fibrinolytic burst) restart from the beginning every
time, **manufacturing a much sicker patient.** This error actually occurred during
development, and after the fix the 24/48/72-hour ICP values of
`10.09 / 13.235 / 11.85` were confirmed to agree between the two approaches to three
decimal places.

---

## Five computed conclusions

### 1 · ICP tells you almost nothing about the room that is left

Because `dP/dV ∝ P`, the same 10 mL is invisible at ICP 8 and lethal at ICP 25.
This is usually taught with a diagram; here it is **measured**: the computation was
branched at several points in the course, a mass was grown inside the cranial vault
at 0.5 mL/min, and the volume that went in before ICP reached 30 was recorded.

| Time (h) | ICP (mmHg) | CSF volume (mL) | Excess brain water (mL) | **Reserve left to ICP 30 (mL)** |
|---:|---:|---:|---:|---:|
| 0 | 10.00 | 110.7 | 0.0 | **76.5** |
| 8 | 9.32 | 94.0 | −1.7 | **68.5** |
| 12 | 9.31 | 92.9 | −0.9 | **67.0** |
| 18 | 9.75 | 85.9 | 7.9 | **55.5** |
| 21 | 10.28 | 76.5 | 18.8 | **44.5** |
| 24 | 12.03 | 65.3 | 32.9 | **32.5** |
| 27 | 16.40 | 57.9 | 46.2 | **22.0** |
| 30 | 20.97 | 54.2 | 55.6 | **15.0** |

**Over 21 hours the monitored ICP moved 0.28 mmHg and 42% of the reserve
disappeared.** Out to 30 hours, ICP moves 10.97 mmHg and 80% of the reserve is gone.
The two columns are not measuring the same thing, and only one of them is on the
monitor screen.

And the **ledger of which compartment paid** is computed too. At t=0, of the 76.5 mL
of reserve, 58.0 mL is CSF displacement — the cheapest buffer in the head and the
first to empty — and once venous blood has contributed 12.6 mL the veins are flat.
After that there is nothing but elastic storage against an exponential curve. It is
not oedema that eats the reserve — at the 12-hour point the excess brain water is
still −0.9 mL.

### 2 · There is no plateau wave in the model. There is only a feedback loop

There is no oscillator in this code, no waveform generator, no periodic forcing
term. There is one loop.

```
CPP falls → autoregulatory vasodilatation → cerebral blood volume rises → ICP rises → CPP falls further
```

The same 12 mmHg / 3 min fall in arterial pressure was applied at several points in
the same patient.

| Time (h) | Baseline ICP | Peak ICP | Rise | ICP at +30 min | Ca max | Verdict |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 10.00 | 10.87 | 0.87 | 9.64 | 0.261 | no wave |
| 12 | 9.31 | 10.52 | 1.21 | 10.07 | 0.290 | small transient |
| 21 | 10.28 | 11.44 | 1.16 | 10.84 | 0.296 | small transient |
| 24 | 12.03 | 14.50 | **2.48** | 13.14 | 0.303 | large transient |
| 27 | 16.40 | 19.54 | **3.14** | 17.83 | **0.305** | **large transient** |
| 30 | 20.97 | 23.36 | 2.39 | 21.99 | **0.305** | large transient |

Early, when compliance is high, the same stimulus moves 0.87 mmHg and comes straight
back. Once compliance has fallen, the same stimulus moves 3.14 mmHg and has still
not come back 30 minutes later. The `Ca max` column shows the mechanism — the
autoregulator is sitting against its dilatation ceiling (0.305).

And if the loop is the cause, deleting the loop must delete the wave. Removing
autoregulation drops the rise at the 15-hour point from **7.60 → 1.21 mmHg** and at
the 18-hour point from **8.93 → 2.62 mmHg**. That is the cleanest evidence that this
is not a numerical artefact of a stiff ODE.

### 3 · Rosner and Lund are both right, and right for the same patient at different times

A vasopressor does two things to intracranial volume, in opposite directions.
Neither was put in separately; both come out of terms that had to be there anyway
(the compliance sigmoid and the Starling equation).

* **Rosner's arm** — when flow exceeds demand the autoregulator **constricts** and
  arterial blood leaves the box. It requires a living autoregulator, and it is
  **fast** (tone time constant 20 s).
* **Lund's arm** — capillary hydrostatic pressure rises, and where the barrier is
  broken there is nothing to stand against that pressure, so water comes **into** the
  box. It requires a broken barrier, and it is **slow**.

The fast gain and the slow cost do not average. **They cross.**

Noradrenaline +0.20 µg/kg/min, autoregulation intact:

| Elapsed (min) | ΔICP (mmHg) | Δarterial blood (mL) | Δbrain water (mL) | Capillary Pc | Which arm wins |
|---:|---:|---:|---:|---:|---|
| 1 | **−0.58** | −0.72 | +0.01 | 39.6 | Rosner (ICP falls) |
| 5 | **−1.36** | −1.87 | +0.08 | 39.2 | Rosner |
| 20 | −0.93 | −1.92 | +0.33 | 40.2 | Rosner |
| 60 | −0.27 | −1.94 | +0.82 | 42.9 | Rosner |
| 90 | +0.17 | −1.96 | +1.13 | 45.2 | — crossover — |
| 120 | **+0.69** | −1.94 | +1.43 | 47.5 | Lund (ICP rises) |

In a patient without autoregulation the sign is opposite from the first minute
(+0.94, arterial blood +1.08 mL, i.e. passive distension).

The **time** of the crossover is stubborn: 101 minutes if the barrier is closed, and
it still crosses even when the barrier is wide open. So it is not only Lund's
filtration that ends the early gain — most of it is **because the CSF system is an
integrator**. Take blood volume out to lower the pressure and CSF absorption
immediately slows, CSF refills, and the pressure returns to the value its own
hydraulics specify. A purely vascular manoeuvre is **buying a transient against a
compartment that refills**.

The **size** of the late effect is where the barrier lives. The arterial-blood column
delivers almost the same −1.9 mL in every row, so the differences between rows come
entirely from the water column.

Either way the clinical implication is the same, and it is uncomfortable. The bedside
test everyone performs — raise the vasopressor and watch the ICP for a few minutes —
samples from the **only window in which this trade looks favourable**. Rosner
measured that window and Lund measured the twelve hours after it. Averaging over
patients or over observation periods in a regime where the sign changes is the surest
way to manufacture a null trial.

### 4 · Osmotherapy does not treat the lesion. It taxes the healthy brain

The reflection coefficient σ of the blood-brain barrier is the whole of the
pharmacology of osmotherapy, and it is **local**. Where the barrier is intact
σ ≈ 1 and 1 mOsm/kg is a driving force of 19.3 mmHg. Where it is broken σ collapses
to 0, the osmoles cross along with the water, and the gradient does nothing at all.

So the model was made to do an accounting that is impossible at the bedside: of the
brain water that left, **how much came from the injured region?**

The same bolus (3% NaCl 250 mL) was given **varying only the state of the barrier**.

| Barrier openness | ΔICP (mmHg) | Total water (mL) | From normal region | From injured region | Injured-region share | Effective σ_injured |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | −8.24 | −14.63 | −8.74 | −5.89 | **40.3 %** | 0.970 |
| 0.20 | −8.31 | −14.31 | −9.06 | −5.26 | 36.7 % | 0.792 |
| 0.40 | −8.33 | −13.56 | −9.21 | −4.36 | 32.1 % | 0.613 |
| 0.60 | −8.16 | −12.33 | −9.20 | −3.13 | 25.4 % | 0.435 |
| 0.80 | −7.20 | −9.85 | −8.38 | −1.47 | 14.9 % | 0.256 |
| 1.00 | −7.77 | −10.48 | −9.63 | −0.85 | **8.1 %** | 0.078 |

The water contribution of the injured region falls **40.3 % → 8.1 %** as the barrier
opens, tracking the reflection coefficient σ in the last column exactly. Once the
barrier is fully open the osmoles cross with the water and the gradient does nothing
there. The total water removed also falls, 14.6 → 10.5 mL.

A statement sharper than the textbook's comes out of this. **Osmotherapy acts on
tissue that still has a barrier.** Its effect on the lesion therefore disappears over
the days after injury — on exactly the timetable on which inflammatory barrier
opening arrives — whereas its effect on the healthy brain does not disappear.
Osmotherapy late in the course is almost entirely **a manoeuvre performed on normal
tissue**, taking water out of normal brain to make room for the lesion to expand.
The ceiling of this therapy is therefore set not by how swollen the lesion is but by
**how much water the healthy brain can safely lose**.

Tolerance and rebound, likewise, are not separate clinical maxims but **one term read
twice**. Brain osmolyte concentration is not a constant but a state variable, so each
dose faces a target that has moved towards it (tolerance); and the very broken
barrier that makes mannitol powerless lets mannitol into the lesion, so that once it
has disappeared from the plasma it **inverts the gradient** (rebound).

### 5 · Hyperventilation is a loan, and the model computes the interest

Two things come out of the single line `PaCO2_eff = PaCO2 × 24 / HCO3_csf`: because
choroid plexus bicarbonate follows with a 6-hour time constant the effect
**disappears**, and returning a normal PaCO2 against an adapted bicarbonate is
itself a **hypercapnic stimulus**.

PaCO2 held at 30 for 24 hours, then returned to 38:

| Time | HCO3_csf | PaCO2_eff | ΔICP vs control | PbtO2 | CBF penumbra |
|---:|---:|---:|---:|---:|---:|
| 1.1 h | 23.2 | 31.0 | −1.16 | 9.9 | 20.2 |
| 8 h | 21.0 | 34.3 | −1.45 | 5.9 | 10.1 |
| 24 h | 20.1 | 35.9 | −1.67 | 5.2 | 6.1 |
| 26 h (after return) | 20.5 | **44.4** | **+1.88** | 5.2 | 7.1 |

As the bicarbonate follows 24.0 → 20.1, `PaCO2_eff` returns towards normal,
31.0 → 35.9 — because the vessels do not see PaCO2, they see perivascular pH. And
returning PaCO2 to 38 after 24 hours makes `PaCO2_eff` **44.4** against the adapted
bicarbonate, so it acts as a hypercapnic stimulus — ICP ends up **+1.88 mmHg higher
than the control.**

Use the same manoeuvre for only **30 minutes** rather than 24 hours and the tissue
hypoxic debt per 1 mmHg·h of ICP gain is an entirely different number. Same drug,
same dose, and the cost per unit of gain differs by orders of magnitude — because
the gain decays with the bicarbonate while the oxygen debt accumulates for as long as
the vessels are constricted. This is the quantitative content of the recommendation
for "short, targeted hyperventilation only".

---

## Pricing every ICP therapy in its own currency

There are only four volumes in the cranial vault, so there are only four things any
therapy can do. At one moment in one patient (the instant at which ICP has just
crossed 22 during tier-0 management) the computation was branched twelve ways and
each therapy was applied for 1 hour. **The point is not the ICP column but the
columns to its right.**

Branch point: t = 30.7 h, ICP 22.19, CPP 64.6, PbtO2 7.04, CSF 53.9 mL,
noradrenaline 0.060 µg/kg/min.

| Therapy (1 hour) | ΔICP | CPP | MAP | PbtO2 | CBF penumbra | Na⁺ | Osmolar gap | Urine output |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Control (no intervention) | 0.00 | 62.2 | 87.3 | 6.44 | 15.8 | 140.0 | 0.0 | 1.6 |
| EVD drainage at 12 mmHg | **−4.85** | 67.1 | 87.3 | 7.22 | 16.4 | 140.0 | 0.0 | 1.6 |
| PaCO2 38 → 33 | −4.24 | 66.5 | 87.3 | 6.94 | **13.9 ↓** | 140.0 | 0.0 | 1.6 |
| PaCO2 38 → 28 | −7.12 | 69.3 | 87.3 | 7.29 | **11.9 ↓↓** | 140.0 | 0.0 | 1.6 |
| 3% NaCl 250 mL | −8.50 | 70.7 | 87.3 | 7.84 | 17.0 | **144.9** | 0.0 | 1.9 |
| 23.4% NaCl 30 mL | **−10.66** | 72.9 | 87.3 | 8.11 | 17.3 | **145.7** | 0.0 | 1.6 |
| Mannitol 0.5 g/kg | −8.71 | 70.9 | 87.3 | 7.88 | 17.1 | 140.1 | **9.2** | **4.3** |
| Propofol 3 → 6 | −1.15 | 60.1 | **84.1 ↓** | 6.68 | 15.2 | 140.0 | 0.0 | 1.6 |
| Thiopental 5 mg/kg/h | −1.45 | 60.3 | **83.9 ↓** | 6.75 | 15.0 | 140.0 | 0.0 | 1.6 |
| Hypothermia 33.5 °C | −4.98 | 67.2 | 87.3 | 7.61 | **13.2 ↓** | 140.0 | 0.0 | 1.6 |
| Noradrenaline +0.15 | −0.22 | 75.9 | 100.7 | 6.53 | 16.3 | 140.0 | 0.0 | 1.6 |
| Decompressive craniectomy | −1.42 | 63.6 | 87.3 | 6.68 | 16.0 | 140.0 | 0.0 | 1.6 |

CSF drainage buys ICP without spending perfusion, sodium or arterial pressure — the
only free row, and small and exhaustible in exchange. The two hyperventilation rows
buy ICP while **selling penumbral perfusion, 15.8 → 13.9 → 11.9.** Hypertonic saline
pays in sodium (140.0 → 145.7), mannitol in osmolar gap (9.2) and urine output
(1.6 → 4.3 mL/min), and the sedative and the barbiturate pay in arterial pressure
(87.3 → 83.9). That is why tier escalation usually does not arrive alone but paired
with a vasopressor.

Two rows deserve particular attention. **Hypothermia** lowers ICP by −4.98 while
**raising** PbtO2 from 6.44 → 7.61 — because metabolic suppression directly improves
the supply-demand balance when perfusion is short. The price is penumbral flow, which
falls 15.8 → 13.2. And **noradrenaline** at this moment moves ICP by only −0.22 while
raising CPP from 62.2 → 75.9. The 1-hour observation window lies before the crossover
of R4 (about 100 minutes), so Rosner's arm is still winning — but it is already
nearly spent.

---

## Why a global monitor cannot see dying tissue

Jugular bulb oxygen saturation is a flow-weighted global average, so it is dominated
by the healthy majority. The brain parenchymal oxygen probe is placed in the small
minority of compartments that are actually dying, behind a much higher local
resistance. In this model the two are **different numbers computed from different
perfusions**. Over a 48-hour course:

| Condition | Time |
|---|---:|
| ICP < 22 mmHg **while** PbtO2 < 20 mmHg | **41.7 h** |
| SjvO2 > 55 % **while** PbtO2 < 20 mmHg | **42.2 h** |

This is time that a protocol watching ICP alone cannot see, and far from being rare
it is most of the trajectory. This is the entire logic of BOOST-2/3.

Let us also be clear about what this does **not** say: it does not say that treating
PbtO2 improves outcome, only that ICP and PbtO2 carry different information. The
reason BOOST-3 exists is precisely that the first of those statements does not imply
the second.

---

## Running tiered therapy as a closed loop: tier 0 alone is worse than no therapy

The SIBICC tier algorithm was implemented as **feedback** rather than as events on a
clock. The same patient (severity 0.62, 42 years old) managed for 72 hours under
several protocols:

| Management | Peak ICP | Hours ICP>22 | ICP burden (mmHg·h) | PbtO2 burden | Final infarct | P(poor) | P(death) |
|---|---:|---:|---:|---:|---:|---:|---:|
| No ICP therapy | 27.9 | 47.2 | 340.5 | 728.2 | 0.198 | 0.927 | 0.531 |
| **tier 0 alone** (sedation + CPP 65) | **37.6** | 47.2 | **688.6** | 770.6 | 0.198 | 0.954 | **0.785** |
| tiers 0–1 (EVD + hypertonic saline) | 22.2 | **0.5** | **7.0** | 657.6 | 0.194 | 0.645 | 0.282 |
| tiers 0–1, mannitol instead | 22.8 | 1.2 | 6.8 | 657.5 | 0.194 | 0.644 | 0.282 |
| tiers 0–2 (+hyperventilation) | 22.2 | 0.7 | 7.9 | 653.4 | 0.196 | 0.652 | 0.288 |
| tiers 0–3 (+barbiturate/craniectomy) | 22.2 | 0.7 | 7.9 | 653.4 | 0.196 | 0.652 | 0.288 |
| ICP threshold 18 (aggressive) | 18.5 | 0.0 | 0.0 | 639.8 | 0.197 | 0.616 | 0.262 |
| ICP threshold 25 (permissive) | 25.2 | 9.7 | 42.2 | 669.9 | 0.195 | 0.728 | 0.356 |
| Autoregulation lost | 22.0 | 0.0 | 18.1 | 669.5 | 0.200 | 0.688 | 0.318 |
| Age 68, identical injury | 22.2 | 0.7 | 7.9 | 653.4 | 0.196 | **0.823** | **0.534** |

Read three things, in order.

**First, the ladder works, and almost all of it works on one rung.** Going from
sedation alone to sedation + drainage + osmotherapy drops the ICP burden from
688.6 → 7.0 mmHg·h. Tiers 2 and 3 **add nothing** in this patient — because tier 1
is already holding it at the threshold, so the escalation criterion is never met
again. That is the quantitative expression of why those are upper tiers and not
routine prescriptions.

**Second, the protocol does not follow a timetable; it argues with the disease.** The
log shows escalation → osmotic bolus → stabilisation and de-escalation 6 hours later
→ oedema progresses and it escalates again, repeatedly, throughout the course.
Nowhere in the code is it written that it should do that.

**Third, and this is the uncomfortable one: tier 0 alone is worse than giving no ICP
therapy at all** (ICP burden 688.6 vs 340.5, probability of death 0.785 vs 0.531).
It is not a bug and not an artefact of the outcome model. Tier 0 here means sedation
plus a noradrenaline loop defending CPP 65, and R4 has already shown what a
vasopressor does in a brain with a leaking barrier over hours rather than minutes.
**Defending perfusion pressure alone, with no means of taking volume out, is running
Lund's arm with nothing to oppose it.** Drainage and osmotherapy are not merely
additional treatments — they are **what makes the vasopressor safe**.

---

## Decompressive craniectomy: why DECRA and RESCUEicp diverged

A young patient (28 years old) with a severe injury (severity 0.92) and almost no
reserve, managed with tiers 0–2:

| Craniectomy ICP threshold | Peak ICP | ICP burden | Final infarct | P(poor) | P(death) | Time performed |
|---|---:|---:|---:|---:|---:|---|
| Not performed | 29.8 | 158.0 | **0.301** | 0.835 | 0.462 | — |
| 16 (very early) | 23.5 | **59.1** | **0.301** | 0.779 | **0.383** | t = 1.9 h |
| 20 (DECRA-like) | 25.7 | 114.0 | **0.302** | 0.818 | 0.435 | t = 2.6 h |
| 25 (intermediate) | 27.6 | 177.2 | **0.301** | 0.842 | 0.473 | t = 25.7 h |

**The ICP burden column and the final infarct column say different things, and that
disagreement is the answer.** The earliest threshold cuts the pressure-time burden
substantially, 158.0 → 59.1 mmHg·h. Yet the final non-surviving tissue fraction is
**the same to three decimal places in every row (0.301).** Enlarge the box and the
**pressure** the brain experiences changes; how much of it survives does **not**.

This is RESCUEicp's result stated as a mechanism — craniectomy is an excellent
pressure therapy and a poor tissue therapy, so it converts death from intracranial
hypertension into survival **with the amount of destroyed brain unchanged**. And to
DECRA's finding that early surgery is harmful, this model says the following: **the
intracranial physiology of early surgery here is plainly favourable. If the trial says
otherwise, the harm cannot be intracranial**, and must live in the things this model
does not have — the surgery itself, the syndrome of the trephined, subdural hygroma,
reoperation, and the fact that randomisation also operates on patients who never
needed surgery in the first place. The very **fact that a mechanistic model fails to
reproduce the harm** localises where the harm lives.

Whether that survival was worth having is something no computation can adjudicate.

---

## Reserve is a function of age and injury severity

Loss of the basal cisterns is a Marshall/Rotterdam CT criterion and the single
strongest CT predictor of raised intracranial pressure. Cerebral atrophy is why an
eighty-year-old tolerates the subdural haematoma that would take a twenty-year-old's
life. **The two are the same quantity** — how much CSF there was to give up — so they
enter the model in one place.

| Injury severity | Age | Initial CSF (mL) | Peak ICP | Hours ICP>22 | Final infarct fraction |
|---:|---:|---:|---:|---:|---:|
| 0.62 | 25 | 88.8 | **31.0** | **51.0 h** | 0.199 |
| 0.62 | 42 | 103.3 | 27.9 | 47.5 h | 0.198 |
| 0.62 | 70 | 127.1 | **20.2** | **0.0 h** | 0.197 |

It is the **same injury**. The infarct fractions are nearly the same too. And yet the
twenty-five-year-old spends more than two days in intracranial hypertension while the
seventy-year-old never once crosses the threshold. That is the result of a single
line.

---

## ⚠️ Where this model is wrong

This section has **not been deleted**, for the sake of honesty.

**Therapeutic hypothermia.** Mechanistically it ought to work. CMRO2 falls with a
Q10 of 2.3, the required flow falls with it, the arterioles constrict, blood leaves
the box and ICP comes down. And, separately, the supply-demand ratio improves so more
of the penumbra ought to survive. Eurotherm3235 and POLAR confirmed the ICP effect
and did not confirm an outcome benefit. Eurotherm reported **harm**.

This model reproduces the ICP effect and the rewarming rebound but **does not
reproduce the worse outcome.** The reason is structural: this model has no pneumonia,
no immunosuppression, no coagulopathy, no shivering and no arrhythmia, and hypothermia
is delivered instantly and for free. Those are precisely the mechanisms through which
the trials see the benefit disappear.

To put it plainly: **a QSP model containing only the mechanisms by which a therapy is
believed to work has no choice but to predict that the therapy works.** This is the
single most important failure mode of mechanistic modelling in drug development, and
the honest response is not to tune the existing parameters until the answer looks
right but to **name the missing compartment**.

Beyond that:

- The **outcome logistic** (R19) is fitted to three calibration anchors, and this is
  stated in the source. It should be read as an **ordering**, not as individual
  probabilities.
- **DECRA/RESCUEicp**: with no surgical complications, no syndrome of the trephined
  and no reoperation, the net harm of early craniectomy is not reproduced. It can show
  why the ICP arms of the two trials came out as they did, but it cannot adjudicate
  whether that survival was worth having — nothing computable can.
- It is a **single-patient deterministic model** with no inter-individual variability.
- The two-compartment penumbra/core tissue model is a crude approximation of real
  regional heterogeneity.

---

## Model structure (47 ODEs)

| Group | State variables |
|---|---|
| Craniospinal mechanics · haemodynamics | `Pic` `x_aut` `MAP` `V_csf` |
| CO₂ / perivascular pH adaptation | `HCO3` |
| Two-compartment brain water · osmoles | `V_int` `V_inj` `Osm_int` `Osm_inj` `Mann_inj` |
| Systemic osmotic state | `Na_ecf` `V_ecf` `Mann_c` `Mann_p` |
| Mass lesion · coagulopathy | `V_hem` `Fib` |
| Excitotoxic · metabolic cascade | `Glu` `K_ec` `Ca_i` `MitoD` `ROS` `Lac` `Glc_br` |
| Tissue fate | `F_core` `F_pen` |
| Neuroinflammation · barrier | `Micro` `Cyto` `Neut` `MMP9` `BBB_mech` `BBB_infl` |
| Circulating biomarkers | `GFAP` `UCHL1` `NfL` `S100B` |
| Temperature | `Temp` |
| Drug PK | propofol ×4, thiopental ×2, noradrenaline, TXA ×2 |
| Cumulative burden | `D_icp` `D_cpp` |

**Closed-loop protocol.** The SIBICC tier algorithm is implemented as **feedback** on
the simulated patient rather than as events on a clock. The tier rises when ICP
exceeds the threshold for 5 minutes and falls when it has been controlled for
6 hours. An osmotic bolus is given only when ICP exceeds the threshold **and** at the
same time sodium and osmolality permit it. Noradrenaline is a PI loop on CPP.
Nothing other than the hysteresis timers reads the wall clock.

---

## References

- Full results and commentary: [`tbi_reference_output.txt`](../../../traumatic-brain-injury/tbi_reference_output.txt)
- 132 references and the parameter-provenance table: [`tbi_references.md`](../../../traumatic-brain-injury/tbi_references.md)
- Structural basis of the haemodynamic core: Ursino & Lodi, *J Appl Physiol*
  1997;82:1256-69
  ([PMID 9104864](https://pubmed.ncbi.nlm.nih.gov/9104864/))
- Tiered therapy: Hawryluk et al., SIBICC, *Intensive Care Med* 2019
  ([PMID 31659383](https://pubmed.ncbi.nlm.nih.gov/31659383/))

> **Disclaimer.** This is a qualitative to semi-quantitative QSP model for
> educational and research purposes. It has not been independently validated or
> certified and must not be used for actual clinical decision-making, prescribing or
> regulatory submission.
