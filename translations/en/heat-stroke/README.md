# Heat Stroke QSP Model

**Exertional heat stroke (EHS) · Classic / non-exertional heat stroke (NEHS)**

The starting point of this model is a single sentence.

> **Heat stroke is the event in which the heat-balance equation loses its fixed
> point. Everything that happens after that is a clock and a dose.**

There are three regimes, and **every intervention that appears in this file acts on
exactly one of them.** This is where the reason that heat stroke is treated so
badly comes from — almost everything that has been tried acts on regime 3.

```
regime 1  COMPENSABLE     A steady state exists. Core temperature plateaus.
                          This is heat strain, not heat stroke.

regime 2  UNCOMPENSABLE   No steady state exists.
                          dTc/dt = (H_prod − Q_loss)/C > 0
                          The clock runs at the rate the equation computes exactly.

regime 3  COMMITTED       The accumulated thermal dose has latched a bistable
                          inflammatory switch. Cooling now returns the temperature
                          but does not return the patient.
```

---

## 1. What this model computes (eight numbers)

Every number below came out of running `hs_analysis.py`, and the raw output is in
[`hs_verification_output.txt`](../../../heat-stroke/hs_verification_output.txt).

### 1-1. The critical wet-bulb temperature is a **result**, not an assumption — and 35 °C is the resting value

Define compensability not as a simulation problem but as the algebraic question
**"does the heat-balance equation have a root"**, and the following comes out.

| Metabolic rate | Critical wet-bulb temperature (30 % RH) | (60 % RH) | (90 % RH) |
|---|---|---|---|
| Rest (84 W) | 34.9 °C | 35.8 °C | **36.4 °C** |
| Light work (200 W) | 31.8 | 32.7 | 33.4 |
| Moderate (400 W) | 28.6 | 29.6 | 30.3 |
| Hard work (600 W) | 25.2 | 26.3 | 27.0 |
| Very hard work (900 W) | **19.7** | 20.8 | 21.6 |

**The 35 °C survivability limit of the climate literature was never entered as a
parameter, and it comes straight out of the equation.** And so does the fact that
it is the **resting limit** — someone who is working has already lost the fixed
point at a wet-bulb temperature of 20 °C. The measurements of Vecellio (2022) and
Wolf (2023) are the second and third rows of this table.

### 1-2. Heat acclimatisation **cannot move** this boundary (and the model says why)

Across all 15 conditions the mean shift is **−0.27 °C** (range −0.51 ~ −0.03), that
is, no meaningful shift, and the sign is not consistently favourable either. The
reason is visible in the solution itself:

```
At every point on the boundary   E_actual == E_max,environment
Sweating capacity               unacclimatised 668 W  /  acclimatised 1233 W
```

**Neither is ever reached.** The limiting factor is the air, not the sweat gland,
and therefore doubling the output of the gland buys nothing.

This is not a hedge but a **falsifiable structural claim**: heat acclimatisation —
and every sweat-side countermeasure — has to act on the **transient and on
cardiovascular reserve** rather than on which environments are survivable. §1-8
verifies that consequence.

### 1-3. The same 40 °C threshold gives the two diseases warning times that differ **19-fold**

| | Reaches 38.5 °C | Reaches 40.0 °C | Reaches 42.0 °C | 40→42 rate of rise | **Warning time** |
|---|---|---|---|---|---|
| EHS 900 W, 35 °C/80 % | 13 min | 26 min | 46 min | 0.100 °C/min | **20 min** |
| EHS 1100 W, football equipment | 9 min | 15 min | 23 min | 0.250 °C/min | **8 min** |
| NEHS elderly + anticholinergic, 40 °C | 17 h | 24 h | 30 h | 0.0053 °C/min | **375 min** |
| NEHS elderly + anticholinergic, 43 °C | 6 h | 10 h | 15 h | 0.0075 °C/min | 268 min |

The same threshold, the same dose law, **two clocks running at rates that differ
19-fold.** The warning the thermometer gives is inversely proportional to the heat
imbalance. This is the **first half** of the reason classic heat stroke kills more
often than exertional heat stroke — not intrinsic lethality but the dose already
paid before arrival.

**A single anticholinergic flips whether the fixed point exists.** In the same
40 °C / 55 % RH room, the elderly person without an anticholinergic crawls along
above the boundary and reaches only 38.1 °C at 72 hours, while the one with it
loses the fixed point and runs away to 40.1 °C at 24 hours.

### 1-4. Time and temperature do not exchange 1:1

Sapareto–Dewey: `dD/dt = R^(43−Tc)`, R = 0.25 (below 43 °C), 0.50 (above).

| Core | 40.0 | 41.0 | 42.0 | 42.5 | 43.0 | 43.5 | 44.0 |
|---|---|---|---|---|---|---|---|
| CEM43 / min | 0.016 | 0.062 | 0.250 | 0.500 | 1.000 | 1.414 | 2.000 |

Every degree below 43 °C **divides the dose rate by four.** Transplanting this dose
law, established in oncological hyperthermia, into heat stroke is the core of this
model, and it is what puts "which cooler" and "when do you start" on **one and the
same axis**.

### 1-5. Cooling modality is a **potency** and delay is an **exposure time** — the exchange rate is about 5 minutes

Total CEM43 paid from collapse at 42.0 °C down to 38.6 °C:

| Modality (fitted UA) | Delay 0 min | 10 min | 30 min | 60 min | Time to cool |
|---|---|---|---|---|---|
| Ice-water immersion 2 °C (23.5 W/K) | **0.90** | 2.97 | 6.22 | 9.52 | 17 min |
| Cold-water immersion 14 °C (27.7) | 1.12 | 3.16 | 6.35 | 9.60 | 22 min |
| Tarp cooling 10 °C (20.1) | 1.35 | 3.36 | 6.49 | 9.69 | 27 min |
| Cold shower 20 °C (15.9) | 1.85 | 3.78 | 6.81 | 9.90 | 40 min |
| Evaporative + convective (16.2) | 2.32 | 4.19 | 7.11 | 10.09 | 51 min |
| Endovascular catheter (4.6) | 3.18 | 4.91 | 7.64 | 10.42 | 68 min |
| Ice packs alone (1.6) | 6.27 | 7.56 | 9.58 | 11.66 | 136 min |
| Passive (0) | 15.67 | 15.67 | 15.67 | 15.67 | 364 min |

**Read a row across and delay overwhelms modality. Read a column down and modality
is worth a factor of a few.** Look at how the difference between best and worst
narrows to 9.52 against 11.66 in the 60-minute column — once you are late, it
scarcely matters what you cool with.

The exchange rate is computed:

| Collapse temperature | CEM43 saved by evaporative → ice-water immersion | Cost of 1 minute of delay | **A better cooler = how many minutes of delay?** |
|---|---|---|---|
| 41.0 °C | 0.38 | 0.062 | **6.0 min** |
| 42.0 °C | 1.42 | 0.250 | **5.7 min** |
| 43.0 °C | 4.57 | 1.000 | **4.6 min** |

**The exchange rate is pinned at around 5 minutes almost independently of peak
temperature, while the absolute risk grows exponentially.** This is precisely the
quantitative content of "cool first, transport second": modality is worth a factor
of a few, delay has no ceiling.

### 1-6. The **dose** that commits is a property of the patient; the affordable **delay** is a property of the cooler

A saddle-node bistability with HMGB1 as the commitment variable:

```
dH/dt = drive + A·H³/(H³+K³) − c_eff·H

Fixed points:  OFF = 0      unstable = 10.91 ng/mL      ON = 50.3 ng/mL
Saddle-node:   c_eff = 0.00498/min  (current c_eff = 0.00307/min)
```

Find, for each modality by bisection, the delay at which the switch latches, and
then read off the dose corresponding to that delay:

| Modality | Delay at commitment | CEM43 paid by then |
|---|---|---|
| Ice-water immersion 2 °C | **44.7 min** | 8.03 |
| Cold-water immersion 14 °C | 43.4 min | 7.99 |
| Tarp cooling | 42.9 min | 8.05 |
| Cold shower | 40.5 min | 8.06 |
| Evaporative + convective | 38.3 min | 8.09 |
| Endovascular catheter | 34.1 min | 8.09 |
| Ice packs alone | **17.2 min** | 8.36 |

**The dose is effectively a constant at 8.1 ± 0.1 CEM43 (spread 0.36), while the
delay ranges from 17 to 45 minutes.** The patient sets the dose; the cooler sets the
delay that can be afforded. These two numbers are the entire reason for cooling in
the field.

(This value of 8.1 is calibrated against the epidemiology: cooling within 30
minutes does not commit and 45 minutes does — a value fitted to the fact that all
274 people immersed immediately in the Falmouth cohort survived, and to the
observation that delayed cases do badly.)

### 1-7. Antipyretics cannot work, and **not because the drug is weak**

Fever raises the set-point and the body defends it. An antipyretic lowers the
set-point and the defence is released. In heat stroke the **set-point is still
37 °C** and the effectors are already saturated with respect to it.

Impose a 0.5 °C shift in set-point and — **the absolute effect is nearly the same in
all three regimes** (of the order of a thousandth of a degree per minute). What
differs by more than an order of magnitude is **what it is a fraction of**: in the
compensable region it is **83 %** of the whole rate of change, in the uncompensable
region **1.0 %**, and not even the sign is reliably favourable there.

The cost is billed separately. In a liver where CYP2E1 is induced and glutathione
is already being consumed:

| | Peak ALT | Trough GSH | Peak core temperature |
|---|---|---|---|
| No paracetamol | 225 U/L | 1.00 | 41.93 |
| 1 g IV | 229 | 0.90 | **41.93** |
| 4 g / 24 h | 240 | 0.75 | **41.93** |

The peak temperature is unchanged to the second decimal place. The liver is not.
Dantrolene has the same structure — RyR1 is **the mechanism of malignant
hyperthermia, not the mechanism of this**, so it reduces muscle heat production by
only 12 %, and the model reproduces exactly the null observed in Bouchama's RCT.

### 1-8. Prevention moves the boundary; rescue only shortens the clock

| Intervention | Time to reach 42 °C | CEM43 at 90 min |
|---|---|---|
| Reference: unacclimatised · no fluid · 900 W | 46 min | 53.96 |
| Heat acclimatisation (10–14 days) | 46 min | 49.58 |
| Fluid intake 0.7 L/h | 47 min | 52.34 |
| Work rate 900 → 600 W | 100 min | 3.49 |
| **Work rate 900 → 400 W** | **never** | **0.33** |
| Wind speed 1.5 → 4 m/s | 65 min | 16.01 |
| Shade (120 W of solar load removed) | 59 min | 24.90 |
| **Humidity 80 % → 40 %** | **109 min** | **1.21** |
| Acclimatisation + fluid + shade | 59 min | 24.79 |

**Only an intervention that restores the fixed point produces `never`.** All the
rest merely buy minutes. And, as §1-2 foretold, acclimatisation and drinking do
almost nothing here either — because in this environment the limiting factor is
the air.

---

## 2. The 16 scenarios

The biology is held fixed and only what the name says is changed. Every difference
among 01–08 is therefore attributable **to the cooling strategy alone**. All arms
reach the same peak temperature (42.0 °C, the collapse trigger) and all are
encephalopathic at that moment — that is not an outcome but the diagnostic
criterion. What separates the arms is **the dose paid afterwards** and **whether the
switch latched**.

| Scenario | CEM43 | HMGB1 | Committed | IL-6 | ALT | CK | Cr | PLT | DIC | GCS 24h | GCS 7d |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 No cooling until hospital (60 min) | 12.42 | **50.3** | **YES** | 1016 | **910** | 89396 | **2.53** | 102 | 5 | 8.2 | **9.6** |
| 02 Field ice-water immersion (<5 min) | **3.70** | 0.0 | no | 834 | 144 | 42916 | 1.26 | 163 | 3 | 11.2 | **15.0** |
| 03 Field CWI 20 min | 6.78 | 0.0 | no | 892 | 225 | 59876 | 1.33 | 156 | 3 | 11.0 | 15.0 |
| 04 Tarp cooling 10 min | 5.31 | 0.0 | no | 865 | 184 | 52685 | 1.30 | 159 | 3 | 11.1 | 15.0 |
| 05 Evaporative cooling 10 min | 6.15 | 0.0 | no | 883 | 207 | 57766 | 1.32 | 156 | 3 | 11.0 | 15.0 |
| 06 Ice packs alone 10 min | 9.67 | 0.0 | no | 956 | 353 | 76125 | 1.46 | 143 | 3 | 10.5 | 15.0 |
| 07 4 °C saline 2 L only | 7.50 | 0.0 | no | 687 | 231 | 61716 | **1.06** | 178 | 3 | 11.7 | 15.0 |
| 08 CWI 20 min + cooled fluids | 6.33 | 0.0 | no | 635 | 193 | 56159 | **1.05** | 183 | 3 | 11.9 | 15.0 |
| 09 Heat-acclimatised runner, CWI 20 min | 6.77 | 0.0 | no | 959 | 204 | **40115** | 1.54 | 138 | 3 | 10.4 | 14.2 |
| 10 Drinking during play, CWI 20 min | 6.78 | 0.0 | no | 891 | 225 | 59874 | 1.16 | 162 | 3 | 11.2 | 15.0 |
| 11 + paracetamol 1 g | 6.78 | 0.0 | no | 894 | 229 | 59866 | 1.33 | 155 | 3 | 10.9 | 15.0 |
| 12 + ibuprofen 800 mg | 6.77 | 0.0 | no | 900 | 230 | 59825 | **1.39** | 152 | 3 | 10.9 | 15.0 |
| 13 + dantrolene 2.5 mg/kg | 6.78 | 0.0 | no | 892 | 225 | 59876 | 1.33 | 156 | 3 | 11.0 | 15.0 |
| 14 + hydrocortisone 200 mg | 6.78 | 0.0 | no | **649** | 225 | 59876 | 1.32 | 158 | 3 | 11.0 | 15.0 |
| 15 60-min delay + thrombomodulin alfa | 12.42 | **0.0** | **no** | 969 | **324** | 89396 | **1.45** | 148 | 3 | 10.7 | **15.0** |
| 16 Classic NEHS, found late | **21.47** | **50.3** | **YES** | 921 | 910 | 87744 | 2.21 | 114 | 5 | 9.0 | **9.6** |

A few ways to read it:

* **01 vs 02.** The same patient, the same collapse temperature. Only the timing of
  cooling differs, 60 minutes against 4. ALT 910 against 144, GCS at day 7 9.6
  against 15.0. A single delay separates permanent neurological sequelae from
  complete recovery.
* **11–14 are barely distinguishable from 03.** That is the point. Neither
  antipyretics nor dantrolene nor steroids change the outcome in this model — except
  that ibuprofen raises creatinine from 1.33 to 1.39 (blockade of renal
  prostaglandins), and hydrocortisone brings IL-6 down from 892 to 649 without that
  leading to any outcome at all.
* **The creatinine in 07 and 08** is the lowest at 1.06 and 1.05. The real benefit of
  cold fluids is not cooling but **volume resuscitation**.
* **15** is the only pharmacological success in this table. With **exactly the same
  thermal history** as 01 (CEM43 12.42, CK 89396), the switch does not latch and GCS
  at day 7 is 15.0.
* **16** has CEM43 21.47 — it arrives having **already** paid a larger dose than any
  EHS arm.

---

## 3. The only drug that touches the switch

In thrombomodulin alfa (ART-123) the lectin-like domain degrades HMGB1. In this
model that is a term **added to the clearance** of the commitment variable, and the
saddle-node calculation gives qualitatively different results depending on dose:

| Dose | c_eff | Result |
|---|---|---|
| None | 0.00307/min | Unstable 10.91, ON 50.3 ng/mL |
| Half | 0.00427/min | Unstable 14.66, **ON 33.0 ng/mL (bistability retained)** |
| Standard 380 U/kg/day | 0.00547/min | **MONOSTABLE — the ON state does not exist** |
| Double | 0.00707/min | MONOSTABLE |

Moving the start time of rTM within the 60-minute delay (the exposure that commits):

| rTM start | Standard dose: committed? | GCS day 7 | Half dose: committed? | GCS day 7 |
|---|---|---|---|---|
| 30 min | no | 15.0 | no | 15.0 |
| 60 min | no | 15.0 | no | 15.0 |
| 120 min | no | 15.0 | **YES** | 9.7 |
| 240 min | no | 15.0 | YES | 9.7 |
| 8 h | no | 15.0 | YES | 9.7 |
| 24 h | no | 15.0 | YES | 9.7 |
| 48 h | no | 15.0 | YES | 9.7 |
| None | YES | 9.6 | YES | 9.6 |

The two blocks have to be read together. **At the standard dose there is no time
window at all** — once c_eff crosses the saddle-node the ON state *ceases to exist*,
so even a switch already latched is reversed. **At half dose the ON state remains at
33 ng/mL and the window reappears** (between 60 and 120 minutes).

That is, in this model **the therapeutic window is not a pharmacokinetic property of
the drug but a property of which side of the saddle-node the dose falls on.** This
is falsifiable, and it is the most exposed therapeutic claim in this model — the
sepsis-DIC trial of the same drug (SCARLET) was null.

---

## 4. Files

| File | Contents |
|---|---|
| [`hs_qsp_model.dot`](../../../heat-stroke/hs_qsp_model.dot) | Mechanistic map source — **128 nodes, 13 clusters** |
| [`hs_qsp_model.svg`](../../../heat-stroke/hs_qsp_model.svg) / [`.png`](../../../heat-stroke/hs_qsp_model.png) | Rendered output (150 dpi) |
| [`hs_mrgsolve_model.R`](../../../heat-stroke/hs_mrgsolve_model.R) | **50-ODE** mrgsolve model + 16 scenarios |
| [`hs_shiny_app.R`](../../../heat-stroke/hs_shiny_app.R) | **10-tab** Shiny dashboard |
| [`hs_references.md`](hs_references.md) | **110** PubMed references, with "how the model uses it" stated section by section |
| [`hs_core.py`](../../../heat-stroke/hs_core.py) | Independent Python/scipy implementation — the executable source of truth |
| [`hs_calibrate.py`](../../../heat-stroke/hs_calibrate.py) | Numerical fit of the cooling conductance UA to published cooling rates |
| [`hs_analysis.py`](../../../heat-stroke/hs_analysis.py) | The script that generates every number above |
| [`hs_verification_output.txt`](../../../heat-stroke/hs_verification_output.txt) | The raw output of that run |
| [`hs_smoke.py`](../../../heat-stroke/hs_smoke.py) | Physiological anchor checks |

---

## 5. Model structure

**50 ODEs.** Thermal nodes 3 (core · muscle · skin) · thermal dose 2 (CEM43 raw /
protected) · HSP70 · water · plasma volume · sweat-gland fatigue 3 · barrier and
endotoxin 2 · cytokines 5 · HMGB1 and NETs 2 · endothelium and coagulation 7 ·
organ injury 13 (including the injured-myocyte pool) · drug PK 11 · cortisol ·
cumulative fluid · 1 for bookkeeping.

**Calibration anchors** (all reproduced in `hs_smoke.py` / `hs_calibrate.py`):

| Anchor | Target | Model |
|---|---|---|
| Ice-water immersion cooling rate | 0.22 °C/min | 0.220 (fitted) |
| Cold-water immersion 14 °C | 0.17 | 0.170 (fitted) |
| Passive cooling | 0.02–0.03 | 0.016 |
| Core fall from 2 L of 4 °C fluid | 1.2–1.7 °C | 1.4–1.7 (incremental) |
| EHS rate of rise (900 W, 35 °C/80 %) | 0.10–0.20 °C/min | 0.100 |
| NEHS rate of rise (40 °C room, 12–24 h window) | 0.1–0.5 °C/h | 0.20 °C/h |
| Resting critical wet-bulb temperature | ~35 °C | 35.7 |
| Muscle–core difference during exercise | 1–2 °C | 0.66–0.9 |
| Time of CK peak | 24–48 h | 25.5 h |
| Time of myoglobin peak | a few hours | 5.2 h |
| Time of AST peak | 24–48 h | 24.3 h |
| Time of ALT peak | 48–72 h | 45.6 h |

---

## 6. Defects found by verification (ten of them)

R/mrgsolve is absent from this build environment, so the Python/scipy
implementation is **the executable source of truth** and the R file is a
transcription from the same equation sheet. Building a runnable implementation and
holding it against the physiological anchors exposed the defects below, all of
which have been fixed.

1. **The cooler's conductance was being bypassed.** The device UA was applied to the
   skin node, but the core–skin conductance created by skin blood flow (up to
   520 W/K) overwhelmed the UA (≈24), so the device became effectively irrelevant.
   The consequence was that **ice water (0.481 °C/min) came out slower than 14 °C
   water (0.495)** — the sign was inverted. The cause was making cold-induced
   vasoconstriction complete, whereas in reality the central hyperthermic drive
   overwhelms the local cold signal (Proulx 2003). Fixed by adding a floor
   proportional to the central drive.
2. **Passive cooling came out at 0.071 °C/min** (measured 0.02). A transient of the
   initial condition, in which the shell was not in equilibrium, was being misread
   as a cooling rate. Fixed by adding a step that solves muscle and skin to steady
   state with the core held fixed.
3. **Two litres of 4 °C fluid dropped the core by 5.08 °C** (measured 1.2–1.7).
   Because the enthalpy sink was applied to the core node alone. Fluid mixes with
   blood and distributes in proportion to perfusion, so this was fixed with a
   perfusion-weighted distribution.
4. **The HMGB1 switch was not bistable.** The autocatalytic term had been written
   `A·H·hill(H)`, which makes production supralinear at every H, so the OFF state
   became unstable and the switch latched **unconditionally**. Saturating production
   `A·hill(H)` against linear clearance is what gives two stable states.
5. **An uncompensable exposure integrated up to 59.6 °C and then went NaN.** Having
   no fixed point *is* this disease, so it must be handled explicitly rather than
   integrated through. Added a roll-off of the Q₁₀ term above 42 °C (protein
   denaturation) plus a termination event.
6. **The cytokines were about 20-fold too large.** TNF settled around 2000 pg/mL
   (measured 50–300), dragging the whole coagulation cascade with it, and the
   platelets of every arm collapsed to 23–29. The production constants for TNF,
   IL-6, IL-1β and IL-10 were recalibrated.
7. **HSP70 was induced up to 58-fold** (measured 2–10-fold). This divided the thermal
   dose by 24 and the commitment drive of the NEHS arms disappeared. Correcting the
   induction rate produced the right physics as a by-product — a 45-minute
   exertional exposure develops almost no thermotolerance, while a multi-hour
   classic exposure develops a good deal.
8. **The myoglobin cast burden accumulated to ~3000** (a variable used as
   `min(1, CAST)`). GFR was pinned at its floor of 0.05 in every arm and creatinine
   came out at 12–17 mg/dL. The coefficient was corrected by three orders of
   magnitude.
9. **Latching the switch changed the outcomes almost not at all** (ALT 130 against
   124). HMGB1's only downstream connection was by way of the cytokines, and with
   that this model's central claim — "cooling returns the temperature but does not
   return the patient" — is hollow. Direct terms from HMGB1 to liver, kidney and
   central nervous system were added. Now the committed and uncommitted arms separate
   at ALT 910 against 144 and GCS at day 7 9.6 against 15.0.

10. **CK peaked 1.7 hours after collapse** (clinically 24–48 hours). Because CK had
    been written to be released simultaneously with the muscle thermal dose, whereas
    injured myocytes keep leaking over a period of hours. Adding a single injured
    myocyte transit pool moved the CK peak to 25.5 hours and, as a by-product,
    spontaneously reproduced the clinical observation that myoglobin (half-life
    90 minutes) peaks before CK (half-life 36 hours) — 5.2 hours against 25.5 hours.

**Verification also refuted one of the author's hypotheses.** The initial
expectation was that a single UA per cooling modality would predict the
water-temperature dependence, but the model predicts the ratio of cooling rates for
2 °C against 14 °C water as **1.55** whereas the measured meta-analytic value is
**1.16–1.29**. The likely causes are that the conductance is linear and that
shivering-thermogenesis offset and incomplete immersion are not modelled. It is
flagged as **the most exposed thermodynamic prediction** of this model; the
scenarios use individually fitted UA values per modality, so the cooling rates
reported there do agree with measurement.

---

## 7. Limitations

* **A single deterministic patient.** There is no between-subject variability.
  Genetic predisposition (RYR1 variants), recent fever, sleep deprivation and
  obesity enter only as parameters.
* **The resting core equilibrium is 36.8 °C** (normal 37.0). Because active cold
  defence was not included, which biases the model towards a slightly later onset
  of NEHS.
* **The commitment dose of 8.1 CEM43 is a value calibrated against the
  epidemiology**, not an independently measured one. The calibration of absolute
  HMGB1 concentrations rests on small cohorts.
* **The cooling conductance UA is an empirical lumped parameter**, not a
  first-principles surface coefficient. It genuinely varies with body surface area,
  body fat and the extent of immersion.
* **The prediction that rTM reverses an already-latched switch** is a direct
  consequence of the saddle-node calculation, but the clinical evidence consists
  only of retrospective registry studies, and the sepsis-DIC trial of the same drug
  was null.
* Death is not modelled. The integration is terminated at 44 °C and the outcome is
  reported only through organ markers and whether commitment occurred.

---

## 8. Running it

```bash
# the mechanistic map
dot -Tsvg hs_qsp_model.dot -o hs_qsp_model.svg
dot -Tpng -Gdpi=150 hs_qsp_model.dot -o hs_qsp_model.png

# regenerate every number (numpy + scipy required, about 20 min)
python3 hs_analysis.py > hs_verification_output.txt

# refit the cooling conductance
python3 hs_calibrate.py

# physiological anchor checks
python3 hs_smoke.py

# the R side (mrgsolve required)
Rscript -e 'source("hs_mrgsolve_model.R")'
shiny::runApp("hs_shiny_app.R")
```

---

> ⚠️ **Disclaimer.** This is a semi-quantitative QSP model for educational and
> research purposes. It has not been independently validated or certified and must
> not be used for clinical decisions, prescribing, or regulatory submission. Real
> heat stroke is an emergency, and immediate cooling in the field is — as it happens,
> also the conclusion of this model — established practice guidance.
