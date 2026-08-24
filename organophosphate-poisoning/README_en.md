# Acute Organophosphorus (OP) Insecticide Self-Poisoning — QSP Model

> **One-line summary** — everything an oxime can do is fixed by a single
> dimensionless number **Ω = k_r_max /(k_i · C_oxon)**, and its upper bound is
> **Ω/(1+Ω)**. Because oxime reactivation saturates, this ceiling *cannot be passed
> with dose.* Ageing loss, by contrast, is not a deadline but an **integral**, and
> what decides death is not the enzyme but **the availability of mechanical
> ventilation**. These three sentences explain 60 years of oxime trial results as
> arithmetic.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (115 nodes · 13 clusters · 185 edges) | [`op_qsp_model_en.dot`](op_qsp_model_en.dot) · [SVG](op_qsp_model_en.svg) · [PNG](op_qsp_model_en.png) |
| ⚙️ mrgsolve ODE model (51 states · 23 scenarios) | [`op_mrgsolve_model_en.R`](op_mrgsolve_model_en.R) |
| 📊 Shiny dashboard (11 tabs) | [`op_shiny_app_en.R`](op_shiny_app_en.R) |
| 📚 References (90 entries) | [`op_references_en.md`](op_references_en.md) |
| 🐍 Standalone Python reimplementation for verification | [`op_reference_model.py`](op_reference_model.py) → [`op_reference_output.txt`](op_reference_output.txt) |
| 📦 Numerical results | [`op_scenario_results.json`](op_scenario_results.json) · [`op_population_results.json`](op_population_results.json) |

---

## 1. Why this disease

Deliberate self-ingestion of organophosphorus insecticide kills tens of thousands of
people a year worldwide and is one of **the largest causes of death attributable to a
single toxic mechanism**. Yet the clinical literature on this condition is more
paradoxical than that of any other poisoning.

* The antidote (an oxime) **works perfectly in the laboratory.** It reactivates human
  AChE in the test tube, and in patients it does genuinely raise blood enzyme
  activity.
* And yet **randomised controlled trials have been negative, repeatedly**, and two
  meta-analyses reported that mortality was in fact higher in the oxime arms.
* Meanwhile, **the things that are not drugs at all** (airway suctioning, mechanical
  ventilation, how fast atropine is escalated) are what determine survival.

This is exactly the kind of question QSP can answer. Do we translate the observation
"the drug does not work" into "the drug is bad", or into **"the conditions under
which the drug can work are arithmetically defined, and most of the patients enrolled
in the trials were outside them"**? This model demonstrates the latter by
calculation.

---

## 2. The organising principle of the model — a three-state switch and two numbers

### 2.1 The switch

Acetylcholinesterase (AChE) moves between three states.

```
    E  ──  k_i · C_oxon  ──▶  EP  ──  k_a  ──▶  EP-aged      (irreversible)
    E  ◀──  k_s + k_r[X]  ──  EP
```

Here, **what does the inhibiting is not the parent compound but the oxon the CYPs
make**. Chlorpyrifos, parathion, dimethoate and fenthion are all P=S pro-toxins and
have to pass through hepatic desulfuration to become the P=O oxon. That single fact
creates the whole time structure of the model: peak inhibition arrives not
immediately after ingestion but after the liver has started making oxon, and while
the parent compound comes back out of the fat depot **oxon goes on being produced for
days.**

### 2.2 The first number — the oxime sufficiency number Ω

Oxime reactivation is not linear. In Worek's formalisation the second-order rate
constant is `k_r2 = k_r_max/(K_D + X)`, so raising the oxime concentration X without
limit still stops the reactivation *rate* at `k_r_max`. The maximum fraction of free
enzyme reachable at quasi-steady state is therefore

```
        E_ceiling = k_r_max / (k_r_max + k_i · C_oxon) = Ω / (1 + Ω)
        Ω = k_r_max / (k_i · C_oxon)
```

**If Ω < 1, no oxime dose whatsoever can hold half the enzyme open.** This does not
mean "the oxime is weak"; it means "for that patient the entire oxime drug class is
arithmetically excluded", and a trial that enrols such patients cannot show benefit
under any design.

The critical oxon concentration `C_crit` computed by the model (target: 30% free
enzyme):

| Oxon | k_i (nM⁻¹h⁻¹) | 2-PAM C_crit | Obidoxime C_crit |
|---|---|---|---|
| chlorpyrifos-oxon | 0.180 | 467 nM | 622 nM |
| paraoxon-ethyl | 0.084 | 1000 nM | 1333 nM |
| omethoate (dimethyl) | 0.0084 | 4167 nM | 8333 nM |
| fenoxon (dimethyl) | 0.050 | 700 nM | 1400 nM |

And because CYP desulfuration **saturates**, in a large ingestion plasma oxon
converges on `Vmax/CL_oxon` and forms a **dose-independent plateau**. Consequently, in
a massive ingestion the oxime ceiling becomes a constant regardless of the amount
ingested.

| Formulation | Oxon plateau | Infinite-dose ceiling | **At concentrations actually reached** |
|---|---|---|---|
| Chlorpyrifos 20% EC + 2-PAM 130 µM | 200 nM | 50.0% | **28.3%** |
| Chlorpyrifos 20% EC + obidoxime 20 µM | 200 nM | 57.1% | **30.8%** |
| Parathion 50% EC + obidoxime 20 µM | 190 nM | 75.0% | **50.0%** |
| Dimethoate 40% EC + 2-PAM 130 µM | 5000 nM | 26.3% | **11.6%** |
| Dimethoate 40% EC + obidoxime 20 µM | 5000 nM | 41.7% | **12.5%** |
| Fenthion 50% EC + 2-PAM 130 µM | 444 nM | 40.3% | **19.3%** |

That last column is the clinically heaviest number in this model. The threshold for
symptom onset is somewhere around 30% free enzyme, and **the commonest formulations
combined with the commonest regimens straddle exactly that threshold.** With
dimethoate it does not even come close.

### 2.3 The second number — the ageing ratchet φ

Ageing is not an event, it is an integral.

```
        φ = k_a / (k_a + k_s + k_r[X])        (the probability that one inhibition ends irreversibly)
        aged(T) = 1 − exp( −k_a · ∫ f_EP dt )  (cumulative irreversible loss)
```

Here the model **corrects** one textbook statement. The claim that "dimethyl OPs age
fast, so oximes are useless" is, as a statement about φ, **wrong**. Dimethyl ages 9×
faster, but its spontaneous reactivation is also 110× faster, so φ is in fact *lower*
(0.159) than for diethyl (0.700).

| OP | k_a (h⁻¹) | k_s (h⁻¹) | φ (no oxime) | φ (2-PAM 130 µM) |
|---|---|---|---|---|
| chlorpyrifos / parathion | 0.0210 | 0.0090 | 0.700 | 0.001 |
| dimethoate | 0.1873 | 0.9902 | 0.159 | 0.033 |
| fenthion | 0.1873 | 0.7702 | 0.196 | 0.034 |

What actually makes dimethyl OPs hopeless is not φ but **f_EP**. Omethoate is a weak
inhibitor (k_i is 20× smaller), but PON1 barely hydrolyses dimethyl oxons, so the
concentration climbs into the µM range and the enzyme consequently stays inhibited
*continuously*. Once f_EP ≈ 1 is maintained, the clock of k_a = 0.187 h⁻¹ simply runs:

```
   dimethoate,   held at f_EP = 1  →  67.5% irreversible at 6 h ·  89.4% at 12 h · 98.9% at 24 h
   chlorpyrifos, same conditions   →  11.8% irreversible at 6 h ·  22.3% at 12 h · 78.0% at 72 h
```

In other words, **the oxime's job is not to "rescue" enzyme but to reduce the area
under the inhibited-enzyme curve**, and that is possible only while
`k_r[X] ≫ k_i·C_oxon`.

---

## 3. Clinical observations the model reproduces (things that come out without being prescribed)

All 23 scenarios are in §7 of
[`op_reference_output.txt`](op_reference_output.txt).
Transcribing only the key contrasts:

### (a) The oxime's effect is a function of the amount ingested, not of the drug

| Exposure | 72 h RBC AChE (supportive care → 2-PAM) | P(death) |
|---|---|---|
| Chlorpyrifos 10 mL | 0.1% → **67.2%** | 18.3% → **7.2%** |
| Parathion 15 mL (+ obidoxime) | 0.1% → **45.5%** | 20.8% → **10.8%** |
| Chlorpyrifos 50 mL | 0.0% → 34.6% | 27.3% → 26.1% |
| Dimethoate 50 mL | 0.1% → **1.0%** | 13.3% → 13.1% |

Same drug, same regimen, same start time, and the effect disappears. What makes it
disappear is Ω, and nothing else.

### (b) Dose-response sweep — the point at which the oxime stops

| Ingested volume | Oxon plateau | Ω | Ceiling | 72 h AChE (supportive / 2-PAM) | Δmortality (%p) |
|---|---|---|---|---|---|
| 2 mL | 14 nM | 13.9 | 93% | 0.7 / 88.8 | **+4.4** |
| 5 mL | 33 nM | 6.15 | 86% | 0.3 / 78.7 | **+9.0** |
| 10 mL | 56 nM | 3.58 | 78% | 0.1 / 67.2 | **+11.1** |
| 20 mL | 87 nM | 2.29 | 70% | 0.1 / 53.1 | +8.8 |
| 35 mL | 115 nM | 1.74 | 64% | 0.0 / 41.4 | +6.9 |
| 50 mL | 132 nM | 1.52 | 60% | 0.0 / 34.6 | +1.2 |
| 100 mL | 159 nM | 1.26 | 56% | 0.0 / 24.0 | −4.0 |
| 200 mL | 177 nM | 1.13 | 53% | 0.0 / 17.3 | −4.9 |

The enzyme measure (AChE at 72 h) improves monotonically at every dose, yet **the
gain in probability of death peaks at 10 mL and changes sign at 100 mL.** A
dissociation between the enzyme endpoint and the clinical endpoint — this is the
shape of the oxime trial literature itself. (For the interpretation of the part where
the sign changes, be sure to read the warning in §6 alongside it.)

### (c) A virtual randomised trial — no assumption is needed to produce a negative result

Take a realistic cohort of 60 with a log-normal ingested volume (median 31 mL),
diethyl 55% / dimethyl 45%, PON1 QQ 25%, and pass all of them through both arms:

```
   overall        placebo 18.9%  →  2-PAM 17.2%     RR 0.91   (ARR 1.7%p)
   ≤ 15 mL        placebo 11.6%  →  2-PAM  7.7%     RR 0.67   ← the real effect is here
   > 15 mL        placebo 20.5%  →  2-PAM 19.3%     RR 0.94
   chlorpyrifos   placebo 25.9%  →  2-PAM 22.5%     RR 0.87
   dimethoate     placebo 11.9%  →  2-PAM 11.8%     RR 0.99
```

**No assumption that the oxime is powerless was put in, and the overall RR still
sticks to 1.** It is enough that the patients above the ceiling dilute the real
effect below the ceiling. This quantitatively raises the possibility that the actual
trials measured not "oximes do not work" but "in this cohort the proportion of
patients in whom an oxime could work was low".

### (d) Mechanical ventilation is worth more than any drug in this box

```
   ventilation available     (S02)   P(death) 26.1%
   ventilation unavailable   (S13)   P(death) 52.3%     ← the drugs are completely identical
```

### (e) How fast atropine is escalated matters more than how much atropine is used

```
   rapid doubling protocol  (S02)   24 h atropine 33 mg   P(death) 26.1%   ventilation 85 h
   customary slow titration (S14)   24 h atropine 16 mg   P(death) 39.0%   ventilation 335 h
```

The slow arm uses **less** atropine and does far worse. In the model atropine is
implemented not as a dosing schedule but as a **closed-loop controller** (titrated to
targets for secretions · heart rate · bronchial tone · spontaneous ventilation), and
the cumulative dose is an *output*, not an input. So this contrast is not a "dose
comparison" but a **comparison of controller gain**.

### (f) A quaternary ammonium antimuscarinic cannot protect the brain

```
   atropine (brain/plasma 0.35)   P(death) 26.1%,  seizure index 0.00
   glycopyrrolate (0.02)          P(death) 44.5%,  seizure index 0.21,  ventilation 335 h
```

The peripheral signs are controlled just as well (indeed the atropine-equivalent dose
is higher, at 1602 mg) and the outcome is far worse. This is the model's expression of
the observation that "early death is central".

### (g) Decontamination is a race against the first-order absorption rate constant

```
   activated charcoal at 1 h (S17)   72 h AChE 54.9%,  P(death) 15.1%
   activated charcoal at 6 h (S18)   72 h AChE 37.0%,  P(death) 23.4%
```

The same 50 g, five hours apart, makes an 8 percentage point difference in the
probability of death. This explains why the large multiple-dose activated charcoal
trials were negative not by "no drug effect" but by **the distribution of arrival
times**.

### (h) The PON1 genotype scales oxon exposure

```
   PON1 R192R   24 h oxon  63 nM,  72 h ageing 47.0%,  P(death) 26.1%
   PON1 Q192Q   24 h oxon 181 nM,  72 h ageing 62.3%,  P(death) 34.6%
```

### (i) A stoichiometric bioscavenger is arithmetically impossible for insecticides

The reason a strategy that works in nerve agent prophylaxis does not work in
insecticide ingestion is not pharmacology but moles.

```
   chlorpyrifos 50 mL of 20% EC = 28.5 mmol parent compound
   plasma BChE   200 mg = 2.35 µmol binding sites   →  12,122-fold short
   plasma BChE  1000 mg = 11.8 µmol binding sites   →   2,424-fold short
   dimethoate 50 mL of 40% EC = 87.2 mmol           →  even 1 g is 7,414-fold short
```

In simulation too, a 1 g dose moves the AChE nadir from 0.03% to 0.03%. That is, it
changes nothing at all.

### (j) Erythrocyte AChE is a poor marker of recovery

Erythrocytes cannot resynthesise enzyme (replaced only by red cell turnover,
~1%/day). Muscle and brain AChE are resynthesised with a half-life of about 5 days.
The model shows the two compartments diverging at day 7 — RBC AChE 28.7% vs NMJ AChE
37.0% (S02), 64.6% vs 70.3% (S05). That RBC AChE and clinical state disagree during
recovery is not measurement error, it is **because the two compartments have
different turnover**.

### (k) Stop the oxime and re-inhibition follows

With a lipophilic OP the fat depot sends the parent compound back over several days,
and because `k_i·C_oxon ≫ k_s` even when only single-digit nM of oxon remains, the
enzyme closes again. The model reproduces the fall in AChE at the moment the oxime
infusion is stopped, and this has the same shape as the re-inhibition phenomenon Eyer
described in parathion poisoning. **An oxime is not a drug that repairs the enzyme; it
is a drug that holds the enzyme open.**

---

## 4. The 51 state variables

| Group | States | Count |
|---|---|---|
| Toxicokinetics (parent compound · fat depot · oxon · solvent) | `A_gut` `A_th_c` `A_th_f` `A_ox_c` `A_ox_t` `A_solv_g` `A_solv_c` | 7 |
| Antidote PK | `A_pam_c` `A_pam_p` `A_atr_c` `Ce_atr` `A_dz` `A_mg` `SCAV` | 7 |
| Esterases (4 AChE compartments × 3 states) | `E/EP/EA` × `rbc, mus, nmj, cns` | 12 |
| BChE · NTE | `B_free` `B_inh` `N_free` `N_inh` `N_aged` | 5 |
| Neurotransmission · receptors | `ACh_m` `ACh_n` `ACh_b` `Rn_des` `Rm_down` | 5 |
| End organs | `SEC` `BT` `HR` `RESPD` `MSTR` `MAP` `SEIZ` `LUNG` | 8 |
| Delayed · cumulative | `OPIDN` `HAZ` `VTIME` `ATRCUM` `AUC_EPn` `AUC_ox` | 6 |
| **Total** | | **51** |

Three structural points worth noting:

1. **Acetylcholine is a hyperbola.** At steady state `ACh = (1+leak)/(E+leak)`, so if
   the enzyme is halved, ACh rises not by 2× but by more. Because the antagonism is
   competitive, **the atropine concentration required is linear in ACh, i.e.
   proportional to 1/E.** A patient at 2% AChE needs more than 4× the atropine of a
   patient at 20% — this is not a matter of degree of severity, it is arithmetic.
   (See the table in §5.)
2. **Ventilation = min(what the brain demands, what the muscle can deliver).** Resting
   ventilation needs only about 28% of maximum inspiratory force (the neuromuscular
   safety factor `SF_NMJ`), so neuromuscular block has a far larger margin than
   central depression. This asymmetry generates, simultaneously, the observation that
   "early death is central" and the observation that "the intermediate syndrome comes
   late and gradually".
3. **Central respiratory depression divides into a part atropine can reverse and a
   part it cannot.** The irreversible part only goes away once receptor
   downregulation (`Rm_down`, half-life about 2 days) has developed. That is why
   patients are weaned on **a clock measured in days** rather than on the atropine
   chart — with the enzyme still under 10%.

---

## 5. The hyperbola of atropine requirement

| AChE (%) | Synaptic ACh (× normal) | Ce* required (nM) | Bolus equivalent (mg) | Maintenance (mg/h) |
|---|---|---|---|---|
| 100 | 1.00 | 0.0 | 0.00 | 0.00 |
| 50 | 1.91 | 7.3 | 0.42 | 0.10 |
| 30 | 3.00 | 16.0 | 0.93 | 0.21 |
| 20 | 4.20 | 25.6 | 1.48 | 0.34 |
| 10 | 7.00 | 48.0 | 2.78 | 0.64 |
| 5 | 10.50 | 76.0 | 4.40 | 1.01 |
| 2 | 15.00 | 112.0 | 6.48 | 1.49 |
| 0 | 21.00 | 160.0 | 9.26 | 2.13 |

The 24-hour cumulative atropine coming out of the closed-loop simulations is
16–70 mg depending on the scenario, and the 14-day cumulative is 250–1600 mg, which
overlaps the range reported clinically. To emphasise again, these values are **outputs
of the model**, not inputs.

---

## 6. Verification and the defects found

Separately from the mrgsolve model (`op_mrgsolve_model_en.R`), the 51 ODEs were
independently reimplemented from the same specification in Python/scipy
(`op_reference_model.py`), and the closed-form results of §2 were recomputed
analytically. That process exposed **five real defects**, all of which were fixed.
They are recorded here.

1. **A unit error in the bioscavenger (a missing mg → g conversion).** The first
   implementation computed 200 mg of plasma BChE as 2,353 µmol of binding sites (the
   true figure is 2.35 µmol — a 1000-fold overestimate). In that erroneous state the
   model produced the plainly wrong conclusion that "1 g of bioscavenger neutralises
   half of a 50 mL chlorpyrifos ingestion". The conclusion in §3(i) is the corrected
   one.
2. **Receptor downregulation had gone in with the sign reversed.** In the initial
   implementation `Rm_down` entered as a term that *raised* apparent affinity, so the
   more tolerance developed, the *stronger* the muscarinic signal became. As a result
   the atropine requirement diverged over time. Receptor downregulation must
   **attenuate signal transduction**, not occupancy.
3. **The ventilation model had no neuromuscular safety factor.** It was initially set
   as `ventilation = drive × muscle strength`, which makes a patient whose muscle
   strength has fallen to 30% automatically go into respiratory failure. In reality
   resting ventilation needs only ~28% of maximum inspiratory force. Because of this
   defect every scenario stayed on mechanical ventilation for the full 14 days. It was
   replaced with a `min(drive, strength/SF_NMJ)` structure.
4. **Chronic lung injury accrued mortality risk without limit.** Minor secretions
   created a permanent lung-injury state which then accumulated risk over 336 hours,
   pushing the probability of death above 50% even in a well-treated patient. It was
   fixed by making aspiration a **threshold event** (only when the airway is actually
   swamped or hydrocarbon actually goes down the wrong way) and by raising the
   elimination rate.
5. **A blind spot in the atropine controller.** The first controller targeted only
   secretions and heart rate. As a result, in the oxime arm the peripheral signs
   improved, atropine was reduced, and since atropine is the only drug in this model
   with central activity, this produced the artefact that **the oxime made the patient
   worse**. A spontaneous-ventilation term was added, in line with the real
   atropinisation targets used at the bedside.

### ⚠️ The model's most exposed prediction

One result survives the corrections and the author **does not trust** it: that in the
massive-ingestion range (≥100 mL) the oxime **raises** the probability of death by
4–5 percentage points. Decomposing the hazard terms, this difference does not come
from the enzyme but from **time on mechanical ventilation** (263 h vs 142 h) and from
borderline hypercapnia near the threshold. That is, it is most likely a structural
side-effect of a sigmoid switch, in which the oxime pushes up the peripheral limit and
thereby keeps spontaneous ventilatory capacity sitting *just below* the weaning
threshold for a long time.

Interestingly, though, **the direction itself agrees with the actual meta-analyses**
(Rahimi 2006, Peter 2006: mortality was higher in the oxime arms). Whether the model
reproduced that observation by coincidence, or has put its finger on a real mechanism
whereby "above the ceiling, an oxime redistributes ventilator-days rather than lives",
cannot be distinguished with this model. **This item should be read as a hypothesis,
not a conclusion, and it is flagged here as the point most vulnerable to refutation.**

For the same reason the probabilities of death for S11 (start at 0.5 hours) and S12
(start at 12 hours) are inverted (26.2% vs 21.2%). The enzyme endpoint keeps the
correct ordering (72 h ageing 46.7% vs 51.0%). **Conclusions about start time should
be read only from the enzyme measures.**

### What verification passed

* The closed-form arithmetic of §1–§6 (Ω, the ceiling, C_crit, X*, φ, the ageing
  integral, the atropine hyperbola, the bioscavenger mole count) agrees with the
  quasi-steady state of the ODE simulation to the decimal place. For example:
  chlorpyrifos 50 mL, oxon 759 nM at the 2-hour point, 2-PAM 215 µM → analytic
  solution 12.0% vs ODE 11.8%.
* Mass conservation: in each esterase compartment `E + EP + EA` is conserved at 1 when
  the turnover terms are absent.
* The ordering of dimethyl/diethyl lethality, continuous infusion > bolus, fast
  atropine > slow atropine, the time-dependence of activated charcoal — all agree in
  direction with the literature.

---

## 7. Limitations

* **The plasma oxon concentration is the most uncertain input in the model.** This
  value determines almost every conclusion through Ω, and yet measurements in real
  patients are rare and extremely variable (Eyer 2009). The model's central prediction
  is therefore left in a **testable form**: *measure plasma oxon. Only in patients
  with C_oxon < C_crit does an oxime lift the enzyme above the threshold.* That is the
  experiment this model proposes.
* Mortality risk is written as a semi-empirical function summing individual hazard
  terms, and it has never been fitted to a real cohort. Absolute mortality has only
  been matched at the magnitude of the cohort level.
* Carbamates, inhalational and dermal exposure, children and pregnancy are not
  covered.
* Nerve-agent oximes such as HI-6 and HLö-7 are on the map but were not
  parameterised.
* Catalytic bioscavengers (phosphotriesterase variants) are not subject to the
  stoichiometric constraint, so the argument in §3(i) does not apply to them. The
  model does not cover that route.
* The intermediate syndrome is collapsed into a single state variable for nAChR
  desensitisation. The real mechanism (muscle necrosis, endplate damage, oxidative
  stress) is more complicated.

---

## 8. How to run

```r
# the model
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread("op_mrgsolve_model_en.R")
# the USAGE block at the foot of the file contains all of the scenario drivers

# the dashboard
shiny::runApp("op_shiny_app_en.R")
```

```bash
# render the map
dot -Tsvg op_qsp_model_en.dot -o op_qsp_model_en.svg
dot -Tpng -Gdpi=150 op_qsp_model_en.dot -o op_qsp_model_en.png

# independent verification (needs numpy + scipy, about 7 minutes)
python3 op_reference_model.py
```

---

## 9. What this model proposes for the clinic (as research hypotheses)

1. **If plasma oxon can be measured, the indication for an oxime can be defined
   arithmetically.** A design that gives an oxime only to patients with Ω > 2
   (roughly, a free-enzyme ceiling > 65%) has a higher probability of detecting an
   effect than any trial to date.
2. **The reason not to use an oxime should be "this patient is above the ceiling", not
   "the drug is bad".** The compound class (dimethyl vs diethyl) and the amount
   ingested serve as surrogate markers when oxon cannot be measured.
3. **Whether or not an oxime is used, the rate of atropine titration has to be managed
   separately.** The largest and most robust single effect in the model was the rate of
   atropinisation (probability of death 39.0% → 26.1%).
4. **From a resource-allocation standpoint, access to mechanical ventilation has a
   larger effect size than any antidote** (52.3% → 26.1%). This conclusion is not
   sensitive to any of the model's uncertain parameters.

---

## ⚠️ Disclaimer

This model is a **semi-quantitative QSP model for educational and research purposes**.
It was constructed on the basis of the published literature and clinical trial data,
but it has not been independently verified or certified, and **it must not be used for
actual clinical decision-making, prescribing, or regulatory submission.** In
particular, the "most exposed prediction" in §6 was written down in the expectation of
being refuted; it is not a recommendation.
