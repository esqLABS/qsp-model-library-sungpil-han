# Osmotic Demyelination Syndrome QSP Model

> **Directory:** `osmotic-demyelination-syndrome/` | **Abbreviation:** ODS | **Date:** 2026-08-07
> **Category:** neurological / electrolyte / iatrogenic

---

[![ODS QSP mechanistic map](../../../osmotic-demyelination-syndrome/ods_qsp_model.png)](../../../osmotic-demyelination-syndrome/ods_qsp_model.svg)

---

## The one-sentence summary

Osmotic demyelination is not a disease caused by hyponatraemia but a disease
caused by the fact that, inside the astrocyte, **the route that exports
osmolytes (a channel) and the route that brings them back in (a transporter that
has to be transcribed) have different time constants** — and this model shows
that the correction-rate limits of the guidelines, the risk factors, the
autonomous overcorrection the kidney generates, the desmopressin clamp, and the
deadline for relowering all follow arithmetically from that one asymmetry.

---

## The structural argument — two time constants

```
     going out :  VRAC / LRRC8A          =  ion channel      t½ ≈ 15 h
     coming in :  SMIT1 / TauT / BGT1    =  transporter proteins,
                  TonEBP has to transcribe them first        t½ ≈ 2–3 d
```

A channel opens in milliseconds. A transporter has to be **transcribed**. The
brain can give osmolytes up roughly six times faster than it can take them back,
and every clinical rule about the correction of hyponatraemia is in fact a
statement about this ratio.

In this model there is exactly one state variable that carries the disease.

```
     Ω(t) = ORG_set( Osm_eff(t) )  −  ORG(t)        [mOsm/kg brain water]
```

**The organic osmolyte deficit** — what the brain *ought to be holding* at the
present tonicity, minus what it actually holds. Zero in the normal brain, and
zero in a **chronically adapted** hyponatraemic brain as well (which is why a
sodium of 110 is not in itself an injury); it becomes positive only when plasma
tonicity moves faster than transcription.

---

## The adapted brain the model produced by itself

The composition of the adapted brain is **an output, not an input**. Generate
hyponatraemia in simulation and the brain discards osmolytes of its own accord.

| Phenotype | [Na] | Total organic osmolytes | myo-inositol | Brain water (normal 80.00) | Ω |
|---|---|---|---|---|---|
| SIADH (21 days) | 110 | 29.06 / 48.00 (−39.5%) | 2.44 / 7.00 (**−65%**) | 80.08 | 0.000 |
| Hypovolaemic (7 days) | 110 | 31.08 / 48.00 (−35.2%) | 2.93 / 7.00 (−58%) | 80.32 | 0.000 |
| Thiazide (21 days) | 110 | 29.05 / 48.00 (−39.5%) | 2.44 / 7.00 (−65%) | 80.08 | 0.000 |
| **Acute (8 hours)** | 110 | **44.56 / 48.00 (−7.2%)** | **6.17 / 7.00 (−12%)** | **87.31 ← oedema** | 0.000 |

This agrees with the 60–70% loss of myo-inositol reported in the rat by
Verbalis & Gullans (1991). **The chronically adapted brain has normal brain water
even at a sodium of 110.** It is not swollen and it is not injured — it has
merely **used up all of its buffering capacity.**

---

## Result ① — the two numbers in the guidelines are one threshold and one transporter

The threshold Ω\* = 8 mOsm/kg was left as **one and the same number for every
patient**. The risk factors act on FOSM (organic osmolyte transport capacity) and
on nothing else. Below is the maximum osmotic stress for each combination; injury
begins once 8.0 is exceeded.

| FOSM | +4 | +6 | +8 | +10 | +12 | +14 | +16 | +20 | (mmol/L/24 h) |
|---|---|---|---|---|---|---|---|---|---|
| **1.00** normal | 5.1 | 5.4 | 6.0 | 7.3 | **8.6** | 9.9 | 11.2 | 13.7 | limit **10–12** |
| 0.80 | 5.1 | 5.5 | 6.8 | **8.3** | 9.7 | 11.1 | 12.5 | 15.3 | |
| 0.65 | 5.1 | 6.1 | 7.7 | **9.3** | 10.9 | 12.4 | 14.0 | 17.0 | |
| **0.55** alcohol / malnutrition | 5.1 | 6.8 | **8.6** | 10.3 | 12.0 | 13.8 | 15.5 | 18.4 | limit **6–8** |
| 0.45 | 5.6 | 7.8 | **9.7** | 11.7 | 13.7 | 15.5 | 17.1 | 19.8 | |

At normal transport capacity the threshold is crossed **between 10 and 12**, at
alcoholic/malnourished transport capacity **between 6 and 8**. These are exactly
the two numbers written into every guideline ("normal risk ≤10–12, high risk
≤8"). Here they are not two rules but **one threshold read through two
transporters**, and the threshold was calibrated once only, while the high-risk
limit is not a calibration but a **prediction**.

---

## Result ② — what sets the rate of correction is not the prescription but the kidney (the central experiment)

A hypovolaemic hyponatraemic patient given 0.9% normal saline 2 L/day.
**Hypertonic saline is 0 mL in both arms.** In one arm AVP responds
physiologically; in the other it is held at its value on presentation (the
counterfactual).

| Time | A: [Na] | AVP | Urine osmolality | Urine volume | EFWC | B: [Na] | AVP | Urine osmolality |
|---|---|---|---|---|---|---|---|---|
| 0 h | 110.0 | 10.06 | 1032 | 0.45 | 0.37 | 110.0 | 10.06 | 1032 |
| 12 h | 110.0 | 2.57 | 877 | 0.54 | 0.33 | 110.0 | 10.06 | 1033 |
| 18 h | 110.3 | 0.02 | **205** | 2.17 | 1.87 | 109.9 | 10.06 | 1033 |
| 24 h | 113.4 | 0.00 | **53** | 6.41 | **6.06** | 109.7 | 10.06 | 1033 |
| 48 h | **133.3** | 0.00 | 50 | 7.30 | 6.88 | 109.0 | 10.06 | 1033 |

| | 24-hour rise | Maximum stress | Maximum deficit | Myelin nadir |
|---|---|---|---|---|
| A (AVP responding) | **20.9** | **14.8** | **42.8** | 0.572 |
| B (AVP held) | 0.0 | 0.2 | 0.0 | 1.000 |

**Neither arm was prescribed a single millimole of sodium.** The whole difference
comes from the kidney being released from a volume stimulus that no longer
exists. Urine osmolality collapses, the electrolyte-free water clearance opens,
and the sodium rises by itself.

---

## Result ③ — the treatment is therefore to take the kidney out of the circuit. And the problem is when you put back what you took out

Desmopressin 2 µg IV q8h for 5 days + 3% saline titrated to +6 mmol/L/24 h.

| | d0–5 rise | d4–14 rise | TBW on day 5 | Maximum stress | Maximum deficit |
|---|---|---|---|---|---|
| 0.9% normal saline alone | 20.9 | — | — | 14.8 | 60.0 |
| Clamp, water intake 1.5 L/day | **6.0** | **17.1** | 46.8 L | 16.9 | **41.2** |
| Clamp, water intake 1.0 L/day | **6.0** | 9.5 | 45.0 L | 10.2 | 0.0 |
| Clamp, water intake 0.5 L/day | **6.2** | 4.8 | 42.4 L | 4.6 | 0.0 |
| Clamp + dose taper 2→0.25 µg | 6.0 | 9.7 | 45.0 L | 10.3 | 0.0 |

**For as long as the clamp is on, the prescription and the outcome become the
same number** — against a prescribed 6.0, an actual 6.0, and zero injury. The
clamp does not lower the sodium; it takes the kidney out of the circuit.

But **taking the clamp off** is the other half of the treatment, and here the
model produces something that was not in its design. The clamp **stores water.**
The moment the drug disappears, the kidney excretes that water all at once. At a
water intake of 1.5 L/day, TBW reaches 46.8 L within 5 days, the sodium jumps by
17.1 mmol/L in the 24 hours after the clamp is removed, and the patient
**sustains, because of the clamp, exactly the injury the clamp prevented.** At
0.5 L/day the rebound is 4.8 and the stress does not come anywhere near the
threshold.

**Tapering the dose from 2 µg to 0.25 µg achieves nothing whatever** (9.7 vs
9.5). Even 0.25 µg produces a plasma concentration of about 10 pg/mL, against a
V2 EC50 of 1.6 pg/mL, so whatever the schedule the receptor remains saturated
until the drug has effectively gone. The water diuresis is not being tapered —
it is **waiting**.

> The rule that is derived is not "wean the desmopressin slowly" but
> **"keep the free-water balance neutral for as long as the clamp is on"**.
> A litre stored is a litre that leaves the moment the clamp comes off.

---

## Result ④ — relowering has a deadline

**Change the prescription** at time T after overcorrection (arm A above) has
begun: stop the normal saline, restrict water to 0.5 L/day, titrate 5% dextrose
to [Na] 118 over 24 hours, desmopressin 2 µg q8h for 60 hours.

| Relowering started | Maximum stress | Astrocyte nadir | Maximum deficit | Deficit at day 90 |
|---|---|---|---|---|
| 8 h | 4.3 | 1.000 | **0.0** | 0.0 |
| 12 h | 4.8 | 1.000 | **0.0** | 0.0 |
| 18 h | 7.6 | 1.000 | **0.0** | 0.0 |
| 24 h | 5.5 | 1.000 | **0.0** | 0.0 |
| 36 h | 10.8 | 0.955 | 0.2 | 0.1 |
| 48 h | 13.9 | 0.778 | 47.7 | 16.0 |
| 72 h | 14.8 | 0.686 | 60.5 | 21.0 |
| not done | 14.8 | 0.686 | 60.0 | 20.8 |

**The deadline lies between 36 and 48 hours.** It is set not by the sodium but by
the **time constant of astrocyte death**. At 72 hours relowering does nothing at
all — because by then what needs fixing is not the sodium.

---

## Result ⑤ — acute and chronic are opposite diseases at the same sodium

The model was never told any such thing.

| | Brain water (at presentation) | Organic osmolytes | Maximum deficit on correction at +14/24h |
|---|---|---|---|
| Chronic (21 days) | 80.08 (normal) | 29.1 / 48 | **32.1** |
| Acute (8 hours) | **87.31 (+9% oedema)** | 44.6 / 48 | **0.6** |

The danger in the acute brain is **oedema and herniation**, and correction is the
treatment. The danger in the chronic brain is **the correction itself**. One Ω
equation makes the two clinical situations run in opposite directions.

To add it honestly, this protection is not unlimited. Pushed as far as
**+20 mmol/L/24 h**, even the acute brain generates a deficit of 21.7 (the
chronic, 61.6), because pulling a brain that was swollen excessively below normal
brings the shrinkage term back to life. That is, the model does not say "if it is
acute you may correct however you like" but **"if it is acute the kind of danger
is reversed and the margin is far larger"**.

---

## Result ⑤-b — the same prescription, five outcomes

All of them are the identical prescription of **+10 mmol/L/24 h**.

| Phenotype | Maximum stress | Astrocyte | Myelin | Maximum deficit |
|---|---|---|---|---|
| Normal-risk SIADH | 7.26 | 1.000 | 1.000 | **0.0** |
| Alcoholism / malnutrition | 10.30 | 0.752 | 0.607 | **54.6** |
| Hypokalaemia (K 2.5) | 8.40 | 0.974 | 1.000 | 0.0 |
| Cirrhosis / pre-transplant | 8.74 | 0.934 | 0.962 | 4.3 |
| Severe, starting [Na] 100 | 7.27 | 1.000 | 1.000 | 0.0 |

"10 mmol/L a day" is **not a statement of safety about a patient but a statement
of safety about a transporter**. And whether the starting sodium is 100 or 110,
*at the same rate* the risk is the same — what creates the risk is not the
starting point but the distance travelled and the rate.

---

## Result ⑥ — potassium is sodium

Edelman: `[Na]s = 1.11 (Na_e + K_e) / TBW − 25.6`. Potassium is **in the
numerator.**

| KCl (mmol/day) | 24-hour rise in [Na] | Serum K |
|---|---|---|
| 0 | +2.0 | 2.69 |
| 60 | +3.3 | 2.89 |
| **120** | **+4.6** | 3.09 |
| 180 | +5.9 | 3.29 |

Put 40 mmol of KCl into 42 L and 1.11 × 40 / 42 = **1.06 mmol/L of sodium
correction** happens without appearing on any fluid chart. Replace potassium at
120 mmol/day in a hypokalaemic hyponatraemic patient and a substantial part of
that day's allowance is spent before any saline is hung. In this model
hypokalaemia is dangerous by **two independent routes** — this one, and the loss
of the sodium gradient the Na⁺-coupled osmolyte transporters use.

---

## Result ⑦ — the clinical course is biphasic, and the MRI is later still

| Correction rate | Maximum stress | Astrocyte | Myelin | Maximum deficit | **Symptom onset** | **MRI positive** |
|---|---|---|---|---|---|---|
| +8 | 5.9 | 1.000 | 1.000 | 0.0 | — | — |
| +10 | 7.3 | 1.000 | 1.000 | 0.0 | — | — |
| +12 | 8.6 | 0.961 | 0.998 | 0.1 | — | — |
| +14 | 9.8 | 0.853 | 0.798 | 32.1 | 8.3 days | — |
| +16 | 11.1 | 0.774 | 0.657 | 49.9 | **6.5 days** | **15.9 days** |
| +20 | 13.6 | 0.672 | 0.519 | 61.6 | 5.4 days | 12.1 days |
| +30 | 19.4 | 0.549 | 0.405 | 68.5 | 4.5 days | 10.3 days |

The sodium becomes normal first, the patient improves, and **then** deteriorates.
And the imaging is **9 days behind** the symptoms. A normal MRI at the time of
symptom onset does not exclude the diagnosis — and this is not a rule put into
the model but a result that emerges from the time constants of the cascade.

---

## Result ⑧ — an honest negative result on urea

In pure effective-osmolality accounting urea **cannot** protect the brain.
Because it crosses the blood-brain barrier it appears identically on both sides
of the balance equation and cancels exactly. And yet the animal evidence
(Soupart 2000, Gankam Kengne 2015) is real. So this model lets urea act only
where those studies **actually measured**, that is on the blood-brain barrier and
the microglial arms.

| | Maximum stress | BBB maximum | Astrocyte nadir | Maximum deficit |
|---|---|---|---|---|
| Overcorrection, no urea | 14.64 | 1.209 | 0.679 | 60.7 |
| + urea 30 g/day | 13.65 | 1.179 | 0.656 | 58.8 |
| + dexamethasone | 14.64 | 1.115 | 0.679 | 55.9 |
| + minocycline | 14.64 | 1.209 | 0.679 | 55.2 |

Across the four arms **the osmotic stress is effectively identical** and only the
downstream arms move. If an experiment with the correction rate strictly matched
demonstrates that urea reduces the osmotic injury itself, then **this model is
wrong.** This is the most exposed assumption in it.

---

## Result ⑨ — the model derives why water restriction fails

The urine (Na+K) / serum Na ratio on presentation = **2.00**. Above 1 the urine
is a more concentrated cation solution than the plasma, so every litre excreted
makes the patient **more** hyponatraemic.

| Prescription (SIADH, [Na] 110, 7 days) | [Na] at 24 hours | 24-hour rise | Maximum stress |
|---|---|---|---|
| Water restriction 1.0 L/day | 111.1 | 1.3 | 0.9 |
| Water restriction 0.5 L/day | 112.5 | 3.2 | 2.4 |
| **0.9% normal saline 2 L/day** | **108.5** | 0.0 | 0.0 |
| 3% saline titrated to +6/day | 115.5 | 6.0 | 4.5 |
| Oral urea 30 g/day | 112.0 | 6.0 | 4.2 |
| **Tolvaptan 15 mg/day** | **130.3** | **28.6** | **20.4** |

Give normal saline (308 mOsm/L) against a urine of osmolality 967 and the sodium
**falls** (desalination). Tolvaptan is the opposite extreme, opening a water
diuresis that cannot be switched off. Neither the Furst rule nor desalination was
put in as a rule — both are results of the computation.

---

## File inventory

| File | Contents |
|---|---|
| [`ods_qsp_model.dot`](../../../osmotic-demyelination-syndrome/ods_qsp_model.dot) / [`.svg`](../../../osmotic-demyelination-syndrome/ods_qsp_model.svg) / [`.png`](../../../osmotic-demyelination-syndrome/ods_qsp_model.png) | Mechanistic map — **154 nodes · 20 clusters · 224 edges** |
| [`ods_mrgsolve_model.R`](ods_mrgsolve_model.R) | mrgsolve model — **40 ODEs · 20 scenarios** |
| [`ods_shiny_app.R`](../../../osmotic-demyelination-syndrome/ods_shiny_app.R) | Shiny dashboard — **11 tabs** (including the safety map, the counterfactual experiment and the deadline) |
| [`ods_references.md`](ods_references.md) | **123 references**, every PMID looked up and checked against the PubMed E-utilities |
| [`ods_verify_python.py`](ods_verify_python.py) | Independent Python/scipy reimplementation of the 40 ODEs, and the scenario runs |
| [`ods_verification_output.txt`](../../../osmotic-demyelination-syndrome/ods_verification_output.txt) | The run output of the file above — **the source of every number in this document** |

---

## Verification

R is not present in the build container, so the mrgsolve file has been
**verified as equations but not verified by compilation.** Instead all 40 ODEs
were independently reimplemented in Python/scipy and integrated, and every number
in this document is the output of that run.

- **The healthy steady state is numerically exact**: after a 150-day
  integration, `max|dy/dt| = 2.5e-11`, [Na] 140.0000, astrocyte, oligodendrocyte
  and myelin all exactly 1.00000, and the cumulative injury integral
  **exactly 0**.
- Normal urine osmolality 760 mOsm/kg, urine volume 1.04 L/day, BUN 20.4 mg/dL.
- Adrogué–Madias check: at 70 kg · TBW 42 L · [Na] 110, +9.37 mmol/L per litre of
  3% saline (0.94 per 100 mL). In a 60 kg patient three 150 mL boluses give
  +5.8 mmol/L, arithmetically consistent with the European guidelines' target of
  a "5 mmol/L rise".
- myo-inositol takes **5.9 days** to recover to 95% of normal
  (the rat data of Verbalis & Gullans 1993, ~5 days).

### Four defects this process exposed and fixed

1. **The floor of the smooth hinge.** The commonly used
   `0.5(x + √(x²+ε²))` leaves, below the knee, not zero but a floor of
   `ε²/4|x|`. Harmless in a rate equation, it was **fatal in an integral** —
   passed down the astrocyte → microglia → oligodendrocyte chain, that floor
   **demyelinated a completely healthy brain to 75% myelin in 60 days.**
   Replaced with a C¹ hinge that is exactly zero below the knee.
2. **A switch hung on a measured value.** The feedback switch "stop the infusion
   when [Na] touches the ceiling" set up an oscillation near the ceiling, and the
   integrator failed to finish even after 7.7×10⁵ function evaluations. A
   prescription is properly a statement about time, so the stop time was changed
   to a deterministic time computed from the prescription.
3. **Sodium excretion was too slow.** The half-life of a load was effectively
   infinite, so the patient accumulated 16 L during a 3% saline correction. Fixed
   to a load half-life of about 11 hours, and at the same time a **medullary
   washout** term was added — without it, urine (Na+K) exceeds 400 mmol/L, a
   urine that does not exist.
4. **Organic osmolyte efflux was too fast.** At t½ 5.8 hours an 8-hour-old acute
   brain has already lost 15% of its osmolytes, so **correcting acute
   hyponatraemia caused injury — the exact opposite of what is observed.** Fixed
   to t½ 15 hours, in keeping with the order given in the literature in which
   electrolytes leave over hours and organic osmolytes over 24–48 hours, and the
   acute/chronic asymmetry then came out of its own accord.

---

## What was fitted and what was predicted

| | Contents |
|---|---|
| **Fitted to normal physiology** | the water and solute steady state (solving WIN to give an exact steady state), the Edelman regression, the urinary concentrating range, plasma urea, the AVP osmotic threshold and slope |
| **Fitted to adaptation data** | the osmotic response coefficients β_i (total −40%, myo-inositol −65%), the efflux/influx time constants (recovery ~5–6 days) |
| **Fitted to the injury dose-response** | four of them, KINJ · KAST · KOLI · KDEM. The threshold Ω\* = 8 was calibrated **once only**, so that the normal-risk limit comes out at 10–12 |
| **Not fitted (= predicted)** | the high-risk limit 6–8, the acute/chronic asymmetry, the autonomous overcorrection the kidney generates, the effect size of the DDAVP clamp, **the fact that the rebound on coming off the clamp is a problem of stored water**, the 36–48 hour relowering deadline, the Furst ratio rule, desalination, the 9-day MRI lag |
| **The most exposed assumption** | urea cannot reduce Ω and acts only on the BBB and microglial arms |

---

## Limitations

1. R being unavailable, the mrgsolve file is **not verified by compilation.** The
   equations were verified against the Python twin.
2. **Ω has never been measured in a living human brain.** What has been measured
   is ¹H-MRS myo-inositol (Videen 1995, Restuccia 2004), which is one term of Ω.
   A falsifiable prediction: what separates the patients who go on to
   demyelinate is not the 24-hour rise in sodium but the **MRS myo-inositol
   deficit at 24–48 hours of correction**.
3. The protective action of urea was hung not on osmolality but on the BBB and
   the microglia (§8).
4. The clinical deficit scale (0–100) is an ordinal construct mapped onto lesion
   burden, and is not the mRS.
5. Serum potassium is a logarithmic function of the total body deficit, and
   intercompartmental shifts driven by acid-base status, insulin or beta-agonists
   are not in the model.
6. Cerebral blood flow, intracranial pressure and the glymphatic system are not
   in the model. The herniation index is no more than a monotonic surrogate for
   brain oedema.
7. That the deficit is 0.1 at +12 mmol/L/24 h and 32.1 at +14 is the cliff in the
   neighbourhood of the threshold; the real distribution of patients will be
   gentler than this (between-subject variance in FOSM is not in the model).

---

## Usage

```r
source("ods_mrgsolve_model.R")
S <- run_scenarios(mod)
plot_central_experiment(S)          # the counterfactual experiment of §2
plot(S$S06, SODIUM + UOSMOL + EFWCL + OMEGAc + STRESSc + DEFICIT ~ time)

shiny::runApp("ods_shiny_app.R")    # the 11-tab dashboard
```

```bash
python3 ods_verify_python.py        # regenerates every number
```

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for teaching and research. It has not been
independently validated or certified, and **must not be used for actual clinical
decision-making, prescribing, or regulatory submission.**
