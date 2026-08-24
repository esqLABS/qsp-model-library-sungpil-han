# Malignant Hyperthermia — QSP Model

> **One driving variable (myoplasmic Ca²⁺), two stores whose sizes differ a hundredfold.**
> The whole of this model is derived from that ratio.

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map (134 nodes · 17 clusters) | [`mh_qsp_model.dot`](mh_qsp_model.dot) · [`.svg`](mh_qsp_model.svg) · [`.png`](mh_qsp_model.png) |
| ⚙️ mrgsolve ODE model (39 compartments · 18 scenarios) | [`mh_mrgsolve_model.R`](mh_mrgsolve_model.R) |
| 📊 Shiny dashboard (11 tabs) | [`mh_shiny_app.R`](mh_shiny_app.R) |
| 📚 References (112) | [`mh_references.md`](mh_references.md) |

---

## 1. The thesis

Malignant hyperthermia is usually learnt as "a list of signs". This model claims that
**the order of that list is not convention but arithmetic on two numbers**.

There is **only one** driving variable — free myoplasmic Ca²⁺ (set by the RyR1 open
probability). Everything the anaesthetist can see is that single driving variable
**integrated into a store**, and the two stores that matter differ by more than two
orders of magnitude in the time they take to fill.

```
CO₂   store   K_CO2B = 60 mL / mmHg
heat  store   C_CORE = 158 000 J / °C
```

Dividing each by **the same excess flux produced by the same event** (in the model's
step experiment, at the 5-minute mark, CO₂ 403 mL/min and heat 160 W):

```
0.149 min / mmHg  (EtCO₂)          16.5 min / °C  (core temperature)
```

That is, **1 °C of temperature consumes as much time as 111 mmHg of EtCO₂**. With
ventilation fixed, the plateau ceiling of EtCO₂ is 149 mmHg — so **the capnograph's
scale is entirely used up before the thermometer moves one degree.** The name
"malignant *hyperthermia*" points at **the thing that happens last**.

### Verification: the step experiment (anaesthetic pre-equilibrated in muscle, trigger pulled at t = 0)

| EtCO₂ | Time reached | | Core temperature | Time reached |
|-------|-----------|---|----------|-----------|
| +10 mmHg | **8.2 min** | | +0.5 °C | 9.8 min |
| +20 mmHg | **18.5 min** | | +1.0 °C | 20.5 min |
| +40 mmHg | 42.8 min | | +1.5 °C | **31.0 min** |
| | | | 38.5 °C | **40.8 min** |

---

## 2. The second axis — the term that starts it and the term that sustains it are different

A volatile anaesthetic opens RyR1. But Ca²⁺ overload plus mitochondrial Ca²⁺
accumulation generates ROS/RNS, and that **oxidises RyR1** (the state variable
`SENS`). Oxidised RyR1 behaves as though 1.3 activity units of anaesthetic were still
present, and its reversal half-life is about 116 minutes — far slower than the rate at
which the anaesthetic washes out.

The consequence of that is the entire history of this disease before 1979:

| Intervention | Outcome |
|------|------|
| No intervention | non-survivable physiology at 140 min |
| Vaporiser shut off + 10 L/min O₂ flush | **non-survivable at 152 min** |
| + activated charcoal filter | **non-survivable at 156 min** |

**Remove the trigger and only the START term disappears; the SUSTAIN term keeps
running.** Dantrolene blocks **the channel itself**, regardless of what opened it.
This is why dantrolene dropped mortality from ~70–80% to a few per cent in an era
when "turn off the anaesthetic" had already been standard for decades.

---

## 3. The third axis — every intervention other than dantrolene treats the "readout"

The model's treatment arms (recognition at 55 min, MHS_high, sevoflurane 2% +
suxamethonium) quantify this. The result is uncomfortable.

| Treatment arm | Peak core temperature | Peak EtCO₂ | CK | Peak K⁺ | Lowest pH | Outcome |
|--------|--------------|-----------|-----|---------|---------|------|
| No intervention | 41.5 °C | 121 | 22 531 | 8.1 | 6.85 | death at 140 min |
| Vaporiser shut off + flush only | 42.0 °C | 123 | 25 611 | 8.2 | 6.85 | death at 152 min |
| + activated charcoal filter | 42.1 °C | 124 | 26 183 | 8.3 | 6.85 | death at 156 min |
| **Hyperventilation ×3 alone** | 42.5 °C | **50.7 (normal)** | 31 894 | 6.8 | 7.21 | **death at 164 min** |
| **Aggressive cooling alone** | **38.1 °C (normal)** | 121 | 15 769 | 7.6 | 6.85 | **death at 141 min** |
| Dantrolene 10 mg/kg alone | 41.3 °C | 65 | 10 497 | 5.5 | 7.17 | survival (GFR 0.71) |
| Bundle − dantrolene | 38.1 °C | 64 | **46 525** | 6.9 | 7.18 | survival (GFR 0.33) |
| **Full bundle + dantrolene** | 37.4 °C | 51 | **4 460** | 5.0 | 7.26 | survival (GFR 0.75) |
| TIVA prophylaxis (no volatile) | 37.0 °C | 40 | 86 | 4.1 | 7.37 | no event at all |

Look at the two "masking" arms. **Hyperventilation alone makes the capnograph normal
and cooling alone makes the thermometer normal, and each of them kills the patient.**
This is not a straw man — both hyperventilation and cooling are legitimate components
of the MHAUS bundle. The point is that **neither is a treatment, and each destroys the
signal that would have prompted treatment**.

And running the identical bundle with only the drug removed isolates the drug's
contribution: **CK 46 525 → 4 460 (tenfold), GFR 0.33 → 0.75.** That is the value of
the only agent that acts on Po.

---

## 4. The fourth axis — the steady-state set point is in the sarcolemma, not in the SR

At steady state the SR release flux and the SERCA reuptake flux **cancel** (otherwise
the SR would empty within seconds). So what determines the sustained myoplasmic Ca²⁺
concentration? **The sarcolemmal balance** — STIM1/Orai1 store-operated entry (gain
`KSOCE`) against PMCA/NCX efflux.

The model therefore predicts that **the trigger is the SR release channel but the
maintenance phase of MH is a disease of extracellular Ca²⁺ entry**, and it does in fact
reproduce the in-vitro observation that removing extracellular Ca²⁺ abolishes the
halothane contracture:

| Ca²⁺ in the bath | 100% | 50% | 20% | 0% |
|-----------|------|-----|-----|-----|
| 2% halothane contracture (g) | 5.63 | 1.30 | 0.17 | **0.00** |

---

## 5. Other results that were derived rather than assumed

**MAC is not a scale for RyR1.** MAC is a measure of anaesthetic potency defined in
the central nervous system, and the volume % corresponding to 1 MAC differs eightfold,
from 0.75% for halothane to 6.0% for desflurane. RyR1 reads molecular concentration
and not MAC, so the time of onset differs by agent (1 MAC, no suxamethonium, minutes
to a sustained rise of EtCO₂ +8 mmHg):

| Genotype | Halothane | Enflurane | Isoflurane | Sevoflurane | Desflurane |
|--------|--------|----------|-----------|-----------|-----------|
| MHN | — | — | — | — | — (no event) |
| RYR1 low penetrance | 87 | 96 | 114 | 136 | 197 |
| CACNA1S | 64 | 70 | 80 | 91 | 117 |
| RYR1 high penetrance | **40** | 42 | 47 | 52 | **61** |
| RYR1-RM | 28 | 29 | 32 | 34 | 38 |

**The weakest trigger produces the latest and most easily missed presentation.**

**Washing the machine is not washing the patient.** Circuit concentration (ppm) after
the vaporiser is shut off:

| Minute | 1 | 5 | 30 | 90 | 180 | 480 |
|---|---|---|----|----|-----|-----|
| 10 L/min flush | 9 572 | 6 819 | 994 | 96 | 39 | 8 |
| Activated charcoal filter | 52 | 9 | 2 | 1 | 0 | 0 |

Flushing alone **cannot get below the anaesthetic the patient is exhaling** — that is
its asymptote. A charcoal filter also adsorbs the anaesthetic in the expirate, so it
has no such floor. Even so, as the treatment arms in the table above show, **the
patient's own muscle store is untouched**, so the clinical outcome barely moves.

**The arithmetic of cooling.** Whole-body heat capacity 232 kJ/°C. One litre of 4 °C
normal saline = 138 kJ = **0.60 °C** (the MHAUS recommendation of 3 L → 1.79 °C).
Surface cooling at 200 W = 3.1 °C/h — almost exactly matching the 249 W of heat
production in fulminant MH. That is, **cooling only buys time; it does not end the
event.** Cooling's only genuinely mechanistic contribution is through the Q10 = 2.6
term (a colder channel opens less), and that is a **second-order effect** on the
driving variable.

**A non-depolarising neuromuscular blocker has precisely zero effect on the
rigidity.** Rocuronium acts at the nAChR, and the Ca²⁺ of MH comes from RyR1 —
**downstream**. Which is why masseter spasm persisting after complete blockade is
diagnostic rather than paradoxical. Once ATP is exhausted the rigidity becomes
**rigor**, and rigor consumes no further ATP and cannot be reversed.

**RQ exceeds 1.** When bicarbonate buffers lactate, CO₂ that did not come from O₂ is
released, so VCO₂ rises more than VO₂ does.

**The contraction throttles its own blood flow.** Metabolic hyperaemia opens muscle
perfusion fivefold, but the contracture compresses the capillaries and closes 40% of
them. Their product sets the aerobic ceiling and drives compartment syndrome.

**Recrudescence is not bad luck but the consequence of underdosing.** A single 2.5
mg/kg dose with no maintenance → `SENS` 0.73 at 4 hours, and as the effect-site
concentration falls below 2 mg/L the event **reignites** (EtCO₂ rebounds by +33 mmHg,
CK 88 413 at 30 hours). Give two doses and `SENS` is 0.23 with **no rebound at all**.

**The price of delay.** Outcome by the delay from recognition to dosing:

| Delay (min) | 0 | 15 | 30 | 45 | 60 | ≥90 |
|-----------|---|----|----|----|----|-----|
| CK | 4 292 | 5 964 | 7 656 | 10 112 | 13 378 | — |
| Peak K⁺ | 5.0 | 5.8 | 6.2 | 6.7 | 7.2 | — |
| GFR | 0.70 | 0.61 | 0.55 | 0.47 | 0.39 | — |
| Outcome | survival | survival | survival | survival | survival | **death at 140 min** |

Beyond about 75 minutes the drug **arrives after the point of no return**. The curve
does not worsen — it **ends**.

**The difference between Ryanodex and Dantrium is not the molecule but the delay.**
Dantrium/Revonto come as 20 mg vials, each of which has to be dissolved in 60 mL of
sterile water, so 2.5 mg/kg for a 70 kg adult needs 9 vials, 540 mL and about 18
minutes of shaking. Ryanodex is 250 mg/5 mL, about 2 minutes. **The only thing that
differs in the model is the input delay**, and the value of those 16 minutes is CK
4 460 → 5 775 (+29%) and GFR 0.75 → 0.69.

**The contracture test (IVCT/CHCT) measures the very parameter the genotype moves.**
Which is why a functional test is still better than sequencing. European MH Group
criteria (≥0.2 g at halothane ≤2% or caffeine ≤2 mM):

| Genotype | 2% halothane (g) | 2 mM caffeine (g) | Verdict |
|--------|--------------|----------------|------|
| MHN | 0.00 | 0.00 | **MHN** ✓ |
| RYR1 low penetrance | 2.31 | 1.11 | MHS ✓ |
| CACNA1S | 3.90 | 2.52 | MHS ✓ |
| RYR1 high penetrance | 5.63 | 5.44 | MHS ✓ |
| RYR1-RM | 6.69 | 6.66 | MHS ✓ |

---

## 6. What was fitted and what was derived

**Only five things were fitted.**

1. `EC50_VOL`, `EC50_CAF` and `PO_BASE` by genotype — so that the simulated
   contracture test classifies MHN as MHN and every MHS variant as MHS under the EMHG
   criteria.
2. `KSOCE` by genotype — so that sustained myoplasmic Ca²⁺ reaches about 1.05 µM.
3. `K_INJ` — so that CK is in the thousands with immediate treatment and above 20 000
   without treatment.
4. Dantrolene PK (Flewellen 1983).
5. The relative RyR1 potency per volume % of the five volatile anaesthetics — so that
   the order of onset at 1 MAC is halothane < enflurane < isoflurane < sevoflurane <
   desflurane (Wedel 1993).

**Everything else is falsifiable output**: the 111 mmHg/°C ratio and the order of the
signs, the failure of removing the anaesthetic on its own, the two masking arms,
dantrolene's tenfold separation in CK, the recrudescence that appears only at the low
dose, RQ > 1, the extracellular Ca²⁺ dependence of the contracture, the arithmetic of
cooling, the futility of non-depolarising blockade, and compartment syndrome.

---

## 7. Where the model departs from the textbooks (not concealed)

**(1) "A rise of 1–2 °C every 5 minutes" is thermodynamically implausible.** That rate
requires a **net** heat accumulation of 500–1 000 W. The model's fulminant event peaks
at 3.5 °C/h, and at that point whole-body VO₂ is 733 mL/min (3.7 times baseline) and
heat production is 249 W. Allowing for the aerobic ceiling (cardiac output ~10 L/min ×
200 mL O₂/L × 0.85 extraction = 1.7 L O₂/min ≈ 570 W), the textbook figure is
possible only **transiently, only in the core compartment, and only when very nearly
the whole musculature is recruited**. The model states this discrepancy rather than
hiding it.

**(2) The contribution of suxamethonium is weak.** In the model, suxamethonium
produces an immediate and self-limiting masseter Ca²⁺ transient but brings the onset of
the sustained event forward by only **about one minute**. That is far weaker than the
epidemiological association in which suxamethonium is involved in the majority of
fulminant MH events, and it is **the most weakly grounded structural choice in this
model**.

**(3) The delay-complication curve is too steep.** A virtual population (genotype mix
45/35/20%, log-normal variability in EC50_VOL, KSOCE, muscle mass, ATPase and K_INJ,
recognition time 55 ± 8 min, 60 subjects per delay):

| Delay (min) | 0 | 15 | 30 | 45 | 60 | 90 |
|-----------|---|----|----|----|----|-----|
| Complication rate | 8% | 18% | 37% | 33% | 50% | 58% |
| Mortality | 2% | 3% | 17% | 8% | 18% | 18% |

The logistic slope gives an **odds ratio of 2.38 per 30 minutes**, against the **1.61**
observed by Larach 2010. The direction and monotonicity are right, but the slope is
excessive. Possible reasons why the real registry data are shallower are (a) the
inclusion of many mild events that never reach a complication regardless of delay,
(b) reverse causation, in which the more severe cases are treated sooner, and (c) the
uncertainty of the recorded time of the "first sign".

**(4) The timing of the recrudescence is late.** The model produces recrudescence only
in the low-dose arm, but the rebound appears late in the 30-hour window. The observed
median in Burkman 2007 is about 13 hours. The reversal half-life of `SENS` (116 min) is
tuned to the acute timescale of hours and decays too quickly to produce a late
recrudescence.

---

## 8. Model structure

**39 ODE compartments** (70 kg adult, time unit = minutes)

| Block | Compartments |
|------|------|
| Volatile anaesthetic PK (6) | circuit · alveolar · vessel-rich group · muscle (well-perfused/bulk) · fat |
| Suxamethonium (1) | plasma amount |
| Dantrolene (4) | central · peripheral 1 · peripheral 2 · effect site |
| Calcium (3) | myoplasm · SR · mitochondria |
| Energy (3) | high-energy phosphate (ATP+PCr) · glycogen · muscle lactate |
| Gas / acid-base (4) | muscle PCO₂ · arterial PCO₂ · whole-body lactate · bicarbonate |
| Heat (2) | core · shell |
| Muscle injury (5) | sarcolemmal integrity · CK · myoglobin · K⁺ · phosphate |
| Kidney (3) | tubular casts · GFR fraction · creatinine |
| Coagulation (4) | thermal dose · fibrinogen · platelets · D-dimer |
| Monitors (4) | heart rate · cumulative VO₂ · cerebral injury index · **SENS (RyR1 oxidation)** |

**Eleven Shiny dashboard tabs**: ① patient · event ② trigger PK ③ **the two stores**
④ Ca²⁺ · energy ⑤ clinical endpoints ⑥ **masking** ⑦ scenario comparison ⑧ time to
dosing ⑨ contracture test ⑩ biomarkers · organs ⑪ model card

**Eighteen mrgsolve scenarios**: MHN control · no treatment · vaporiser shut off alone
· activated charcoal filter · TIVA prophylaxis · hyperventilation alone · cooling
alone · dantrolene alone · bundle − drug · Ryanodex · Dantrium ·
30 min delay · 60 min delay · low-dose recrudescence · low-penetrance variant ·
desflurane 1 MAC · low extracellular Ca²⁺ · antioxidant (NAC)

---

## 9. How to run it

```r
# the ODE model and all the analyses
source("mh_mrgsolve_model.R")
run_all()                    # summary of the 18 scenarios
capacitance_experiment()     # the time constants of the two stores
ivct_panel()                 # the simulated contracture test
ivct_ca_dependence()         # extracellular Ca²⁺ dependence
onset_panel()                # time of onset by genotype × agent
delay_curve()                # the price of a delay in dosing
washout()                    # circuit flush vs activated charcoal
cooling_arithmetic()         # the arithmetic of cooling
virtual_population()         # complication rates in the virtual population

# the dashboard
shiny::runApp("mh_shiny_app.R")

# rendering the map
dot -Tsvg mh_qsp_model.dot -o mh_qsp_model.svg
dot -Tpng -Gdpi=150 mh_qsp_model.dot -o mh_qsp_model.png
```

**How it was verified.** All 39 ODEs were independently re-implemented in
Python/scipy (LSODA) and checked numerically first, then carried over to mrgsolve.
Five defects came to light in the course of that and were corrected:
(i) the non-susceptible (MHN) control drifted into MH by itself after 4 hours — the
RyR1 Hill coefficient and the MHN EC50 were too low; (ii) divergences such as CK
2 500 000 U/L and K⁺ 300 mmol/L — an error in which the myoplasmic Ca²⁺ plateau was
being pinned at 3 µM by the SOCE/PMCA balance; (iii) oxygen supply at steady state was
always sufficient so that no lactic acidosis ever developed — the mitochondrial Ca²⁺
overload-driven uncoupling and the capillary compression by the contracture were
missing; (iv) divergence of the renal model (creatinine 79 mg/dL) — the GFR recovery
term was being divided by the cast burden; (v) a residual drift in which the MHS
patient rose to 38.2 °C even on **a non-triggering anaesthetic** — the force–pCa curve
was too shallow, so a resting Ca²⁺ of 0.27 µM produced 8% muscle activation
(corrected from Hill 3 to 5).

The verification also **refuted one of the author's prior hypotheses**: the initial
model treated cooling as purely "treating the readout", but it emerged that through the
Q10 = 2.6 term cooling really does have a second-order effect on the driving variable.
The account in §3 above has been revised accordingly.

---

## 10. Emergency information and disclaimer

- **MHAUS 24-hour hotline (North America):** 1-800-644-9737 (outside the USA
  +1-209-417-3722) · <https://www.mhaus.org/>
- **European Malignant Hyperthermia Group:** <https://www.emhg.org/>

> ⚠️ This is a **qualitative / semi-quantitative QSP model for educational and
> research purposes**. It has not been validated against patient-level data and
> **must not be used for any clinical decision-making.** If malignant hyperthermia is
> genuinely suspected, follow the MHAUS/EMHG protocol and contact the hotline in your
> country immediately.
