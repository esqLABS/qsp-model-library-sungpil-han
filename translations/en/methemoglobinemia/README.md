# Methaemoglobinaemia — QSP Model
## Quantitative Systems Pharmacology model of methaemoglobinaemia

<p align="center">
  <a href="mhb_qsp_model.svg">
    <img src="mhb_qsp_model.png" width="900" alt="Methaemoglobinaemia QSP mechanistic map">
  </a><br>
  <sub><a href="mhb_qsp_model.svg">View full-resolution SVG</a> · 153 nodes · 209 edges · 19 clusters (18 mechanisms + 1 legend)</sub>
</p>

---

## The claim this model is built to make

> **Methaemoglobinaemia is not a disease of saturation but a disease of delivery,
> and the two numbers available at the bedside are wrong in opposite directions.**
>
> **And the antidote and the poison are the same molecule. The only thing that
> separates the two is a single electron supply (NADPH), and the most important
> comorbidity in this disease — G6PD deficiency — is exactly what cuts off that
> supply.**

As with the other models in this repository, what matters is not "what was put
in" but **"what came out despite not being put in."**

The following are **absent** from this model's list of 139 parameters:

- A severity scale, or a threshold such as "danger above some % MetHb"
- Methylene blue's **7 mg/kg cumulative ceiling**
- A **contraindication switch for G6PD deficiency**
- A **rebound** term for dapsone poisoning
- The pulse oximeter's **85% floor** constant
- A rule that sulfhaemoglobin and HbM do not respond to methylene blue

And yet every result below comes out anyway. What was put in is only five
**structures**.

| What was put into the model (structure) | What the model computed (result) |
|---|---|
| A catalytic co-oxidation cycle whose cofactor is glutathione | Why 250 mg of benzocaine can produce 31% MetHb — stoichiometry alone caps it at 0.5% |
| **One** finite NADPH supply with two consumers | The whole G6PD gradient: working → slowing → failing → **failing while haemolysing** |
| A saturable productive branch plus a non-saturable futile branch | Why the methylene blue dose–response bends back, and why the ceiling is **cumulative, not per-dose** |
| Ferric fraction → P50 and Hill n (Darling–Roughton) | Why the equivalent anaemia for "MetHb 30%" is not 10.5 but **8.31 g/dL** |
| Two extinction coefficients per pigment plus one linear calibration equation | The 85% floor, and the fact that **the point where sensitivity collapses is exactly where the patient is dying** |
| Two clocks running at different speeds | The rebound in dapsone poisoning — the model was never told that rebound exists |

Every figure was actually computed by
[`mhb_reference_check.py`](../../../methemoglobinemia/mhb_reference_check.py), and
the full output is in
[`mhb_reference_output.txt`](../../../methemoglobinemia/mhb_reference_output.txt).
**All 42 anchors pass**, and the 139 parameters in the R file and the Python file
were cross-checked with **zero mismatches**.

---

## 1. The headline: what a percentage actually costs

When a patient is told "your methaemoglobin is 30%," the arithmetic actually
done at the bedside is `15 × 0.70 = 10.5 g/dL`. **That arithmetic is wrong, and
this model computes how wrong.**

A ferric (Fe³⁺) subunit does two things to its neighbouring subunits. It
abolishes cooperativity (n → 1), and it shifts the curve **leftward** (Darling &
Roughton 1942). In other words, the remaining haemoglobin **will not let go** of
the oxygen it is still carrying.

**Equivalent Hb (EQUIV Hb)** = the concentration of structurally normal
haemoglobin that would produce the **same tissue PO₂** at the same VO₂ and the
same cardiac output.

| MetHb % | Naïve arithmetic Hb×(1−f) | **Equivalent Hb** | Hb neutralised | % of remaining capacity | PvO₂ | SpO₂ |
|---:|---:|---:|---:|---:|---:|---:|
| 15 | 12.75 | **11.41** | 1.34 | 10.5% | 34.3 | 90.0 |
| 20 | 12.00 | **10.32** | 1.68 | 14.0% | 32.5 | 88.8 |
| **30** | **10.50** | **8.31** | **2.19** | **20.8%** | 28.6 | 87.4 |
| 40 | 9.00 | **6.57** | 2.43 | 27.0% | 24.3 | 86.6 |
| 50 | 7.50 | **5.17** | 2.33 | 31.1% | **19.1** | 86.1 |

Look at the fifth column. The share neutralised by the leftward shift grows
**with severity** — from 10% at 15%, to 21% at 30%, to 31% at 50%. The error is
largest exactly where the decision matters most.

And at MetHb 50%, PvO₂ is 19.1 mmHg — **below the anaerobic threshold**. Nobody
told the model that "50% is dangerous"; it comes out that way regardless.

---

## 2. The same percentage is not the same disease

%MetHb is a **ratio**. Tissue oxygen delivery is a **product**. When the model
holds both, an anaemic patient's vulnerability stops being a warning label and
becomes arithmetic.

| Hb | MetHb 20% → PvO₂ | MetHb 30% → PvO₂ | SpO₂ (20% / 30%) |
|---:|---:|---:|:--|
| 17.0 | 34.8 | 31.0 | 88.8 / 87.4 |
| 15.0 | 32.5 | 28.6 | 88.8 / 87.4 |
| 12.0 | 28.4 | 24.5 | 88.8 / 87.4 |
| 9.0 | 23.1 | **18.9** ⚠ | 88.8 / 87.4 |
| 7.0 | **18.0** ⚠ | **13.5** ⚠ | 88.8 / 87.4 |

**Look at the right-hand column. The pulse oximeter reads exactly the same
number in every row.** The SpO₂ difference between Hb 7 and Hb 17 is 0.0 points.
The tissue PO₂ difference is 16.8 mmHg. No bedside measurement distinguishes
these two patients.

The textbook table of "15% asymptomatic / 30% dyspnoea / 50% acidosis"
**implicitly assumes Hb 15 g/dL.** For a patient with Hb 8, it is simply the
wrong table.

---

## 3. The 85% floor is a property of the device, not the blood

A pulse oximeter computes exactly one thing, `R = A₆₆₀ / A₉₄₀`, and feeds that
value into a **linear** calibration equation, `SpO₂ = 110 − 25R`. Any pigment
that absorbs similarly at 660 nm and 940 nm pulls R toward 1, and evaluating the
line at R = 1 gives `110 − 25 = 85`.

| MetHb % | R | SpO₂ | Co-oximetry SaO₂ | Gap | **−dSpO₂/d(%MetHb)** |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.346 | 100.0 | 96.8 | 3.2 | — |
| 5 | 0.608 | 94.8 | 91.9 | 2.9 | **0.90** |
| 15 | 0.800 | 90.0 | 81.9 | 8.1 | 0.35 |
| 30 | 0.902 | 87.4 | 67.1 | 20.4 | 0.14 |
| 50 | 0.956 | 86.1 | 47.3 | 38.8 | 0.05 |
| 70 | 0.982 | 85.4 | 27.8 | 57.6 | 0.03 |
| 85 | 0.993 | 85.2 | 13.6 | **71.6** | **0.02** |

The last column is the point of this section. **The device's sensitivity
collapses as the disease worsens.** Over 0–5% it moves 0.9 points per percentage
point; over 60–85% it moves 0.02. The device effectively stops responding
exactly in the range where the patient is dying, and it does so smoothly, with
no alarm.

As a by-product, the visible cyanosis threshold also comes out of the model:
**MetHb 9.1%** (clinical teaching says "10–15%"). All that was needed was to
write cyanosis as the sum of three pigments, each given a threshold
(deoxyhaemoglobin 5 g/dL, MetHb 1.5, sulfHb 0.5).

---

## 4. Methylene blue — one molecule, two directions

The model has no maximum-dose parameter. What it has is a single **fork**.
Leucomethylene blue (LMB) can donate two electrons either to methaemoglobin or
to oxygen. The first branch is proportional to `MHB/(Km+MHB)`, so it
**saturates and switches off once the job is done.** The second branch does not
saturate, and **keeps running for as long as the drug is present.**

| MB (mg/kg) | MetHb @2h | Hb @48h | Heinz bodies | ROS exposure (AUC) |
|---:|---:|---:|---:|---:|
| 0.5 | 9.82 | 14.98 | 0.002 | 168 |
| 1.0 | 5.05 | 14.94 | 0.011 | 308 |
| 2.0 | 2.86 | 14.68 | 0.060 | 559 |
| **3.0** | **2.50** ← lowest | 14.30 | 0.133 | 783 |
| 5.0 | 2.75 | 13.53 | 0.276 | 1177 |
| 7.0 | 3.44 | **12.91** | 0.392 | 1516 |
| 10.0 | 4.84 | **12.25** | 0.511 | 1947 |

The methaemoglobin column bottoms out and then **climbs back up**, while the
haemoglobin and Heinz-body columns keep getting worse. **This is why the
clinical ceiling is a "cumulative dose," not a "single dose":** the harm is the
**time integral** of the futile branch running for as long as the drug is
present, while the benefit switches off the instant the substrate is gone. A
constraint on a time integral is inherently a cumulative constraint.

---

## 5. G6PD — why the antidote fails, and then does harm

The model contains no rule saying "methylene blue is contraindicated in G6PD
deficiency." What it contains is a single finite NADPH supply with two
consumers: the flavin reductase that makes leucomethylene blue, and the
glutathione reductase that keeps the membrane intact.

| G6PD activity | MetHb @2h | Lowest GSH (µM) | Heinz bodies | Hb @72h | Verdict |
|---:|---:|---:|---:|---:|:--|
| 1.00 | 2.85 | 1897 | 0.06 | 14.65 | Works |
| 0.40 | 16.00 | 570 | 0.03 | 14.83 | Works |
| 0.25 | 22.37 | 120 | 0.23 | 13.83 | Works slowly |
| 0.15 (A⁻) | 22.82 | 50 | 0.44 | 12.86 | Fails |
| 0.08 | 18.26 | 22 | 0.63 | 12.06 | **Fails + haemolyses** |
| 0.02 (Mediterranean) | 12.85 | 5 | 0.83 | **11.10** | **Fails + haemolyses** |

What emerges is not a binary rule but a **gradient**. And the GSH column shows
why: at low activity, methylene blue is not merely ineffective — it **competes
with glutathione for the electrons that were holding the cell together.**

One more thing — moving down the table, MetHb @2h **drops back down** from
22.8% (12.85%). This is not an error but a prediction: when G6PD is extremely
low, even the oxidant's **co-oxidation cycle itself** cannot regenerate
glutathione and stalls. In other words, a patient with severe G6PD deficiency,
on oxidant exposure, **makes less methaemoglobin and haemolyses instead** — and
that is exactly the clinical picture actually observed.

---

## 6. Dapsone — two clocks, and an inevitable rebound

| Dapsone | Steady-state MetHb |
|---|---|
| 50 mg/d | 3.07% |
| 100 mg/d | 5.69% |
| 200 mg/d | 12.39% |
| 300 mg/d | 21.61% |
| 100 mg/d **+ cimetidine** | 3.57% (**a 37% reduction**) |

Dapsone 2 g overdose, untreated: peak MetHb **58.7%** (at 26.8 hours), **above
20% for 80 hours.**

Give a single dose of methylene blue 2 mg/kg at 6 hours, and:

> 27.2% → 2.7% (at 7.9 hours) → **rebounds back to 25.8%** (at 52.2 hours)

The model was never told that dapsone poisoning rebounds. It rebounds because
two clocks run at different speeds. **Clock 1 (the cause)**: dapsone has a
half-life of roughly 30 hours plus enterohepatic recirculation, so the oxidant
influx lasts for **days**. **Clock 2 (the antidote)**: methylene blue's effect
on the blood is over within **one hour**.

And that leads to a therapeutic argument — the question of **which term is
rate-limiting**:

| Strategy | MetHb AUC (%·h) | Hb @144h | Cumulative MB |
|---|---:|---:|---:|
| Methylene blue ×5 (up to the ceiling) | 1052 | 12.64 | 7 mg/kg |
| Methylene blue ×2 + cimetidine | **913** | **13.64** | 4 mg/kg |

Methylene blue acts on the **elimination term** and is exhausted within an
hour. Cimetidine acts on the **generation term** and keeps acting. In a
poisoning whose cause lasts for days, a short-acting drug must be re-dosed
until it hits the cumulative ceiling, but the strategy that reduces the cause
carries no such constraint.

---

## 7. Chronic 25% is fine; acute 25% is not

Patients with congenital CYB5R3 deficiency are visibly blue, run
methaemoglobin of 15–30%, and are **perfectly well.** The model has no
"chronic" switch. What it has is 2,3-BPG and bone marrow, both of which respond
to tissue hypoxia. And **what 2,3-BPG does to P50 is exactly the opposite
direction of what methaemoglobin does.**

| | MetHb | 2,3-BPG | P50eff | Hb | **PvO₂** | SpO₂ |
|---|---:|---:|---:|---:|---:|---:|
| Chronic (250 days) | 22.5% | 1.09× | 26.4 | 15.08 | **35.5** | 88.3 |
| Same % but acute | 22.5% | 1.00× | 23.8 | 15.00 | **31.6** | 88.3 |

**Same percentage, same oximeter reading, a 3.9 mmHg difference in the tissue.**
This compensation is neither "tolerance" nor "the brain adapting" — it is a
measurable rightward shift offsetting part of a measurable leftward shift.

For the same reason, methylene blue in congenital methaemoglobinaemia is a
**cosmetic procedure.** The colour changes and tissue PO₂ barely moves.

In the **HbM variant** (FHBM = 0.25), giving methylene blue does not move
MetHb at all. The model has no rule stating "HbM does not respond" — the HbM
ferric haem simply never enters the substrate term of the reduction equation in
the first place. Sulfhaemoglobin is handled the same way.

---

## 8. Bypassing haemoglobin altogether — Boerema's arithmetic

The three classes of treatment attack **different terms**:

| Term | Treatment | Character |
|---|---|---|
| **Generation** | Stopping the causative drug · cimetidine · activated charcoal | Slow but **continuous** |
| **Elimination** | Methylene blue · ascorbic acid · riboflavin | Fast but has a ceiling and needs NADPH |
| **Delivery** | Hyperbaric oxygen · exchange transfusion | The only thing that works when 1 and 2 cannot |

| Situation | Dissolved O₂ | PvO₂ | |
|---|---:|---:|:--|
| MetHb 70%, room air | 0.29 mL/dL | 5.4 mmHg | Anaerobic |
| MetHb 70%, 100% O₂ (1 ATA) | 1.92 | 16.7 | Still anaerobic |
| MetHb 70%, 100% O₂ @ 2.8 ATA | 5.40 | 131.5 | Aerobic |
| **MetHb 95%, 100% O₂ @ 2.8 ATA** | **5.40** | **131.6** | **Aerobic** |

The arithmetic is Boerema's. Resting extraction is `VO₂/(CO×10) = 5.0 mL/dL`,
and oxygen at 2.8 ATA dissolves `0.003 × 1800 = 5.4 mL/dL` in plasma. **5.4 >
5.0.** In other words, hyperbaric oxygen covers resting metabolism even with
essentially no functioning haemoglobin — and that is precisely the situation in
which methylene blue no longer has any substrate left to act on.

Ascorbic acid is not "a slow methylene blue." It is a non-enzymatic bimolecular
reaction, and no dose the kidney can tolerate makes it competitive. It is used
**when methylene blue is contraindicated**, not in place of it when it is not.

---

## 9. The pharmacology of the antidote itself

Methylene blue is a potent MAO-A inhibitor (IC50 ≈ 5 µM).

| | Peak synaptic serotonin (normal = 1) |
|---|---:|
| Methylene blue alone | 1.29 |
| SSRI alone | 1.62 |
| **Methylene blue + SSRI** | **2.62** |

---

## 10. A falsified hypothesis — left in rather than deleted

The autocatalytic term for nitrite (`KNIT2`) was added **expecting it to
produce the clinically known "lag then surge" shape. It failed to.**

| | Peak MetHb at 5 g | Time from 5% → 30% |
|---|---:|---:|
| With autocatalytic term | 62.1% | 0.41 h |
| Autocatalytic term removed | 61.6% | 0.41 h |

The reason is itself a result: at poisoning doses, this reaction is **not
rate-limited by its own chemical speed.** On the way up, **absorption** is
rate-limiting; at the peak, **stoichiometry** is. Sodium nitrite 5 g is 72.5
mmol against roughly 186 mmol of haem, and unlike arylhydroxylamines, **nitrite
is consumed as it reacts.** That is why, for nitrite patients, "how much was
ingested" predicts outcome far better than it does for dapsone patients.

The autocatalysis is real chemistry. It simply is not the rate-limiting step in
vivo. The term has been left in the model, and this section has been left in
the output. **Deleting a failed hypothesis is how a model stops being
falsifiable.**

---

## 11. The most fragile assumption (not hidden)

This model's headline result depends **entirely on `ALPHAM`** (the magnitude of
the Darling–Roughton leftward shift). Yet the **quantitative human data** on
exactly how much methaemoglobin lowers the P50 of the remaining ferrous
subunits are **sparse.**

| Leftward shift | P50eff (MetHb 30%) | Equivalent Hb at MetHb 30% |
|---|---:|---:|
| 0 (none) | 26.8 | **10.39** — the naïve arithmetic (10.50) becomes essentially correct, and this model's claim disappears |
| Half | 24.8 | 9.31 |
| **Model value (α = 0.50)** | 22.8 | **8.31** |
| 1.5× | 20.8 | 7.40 |

**This model is falsifiable, and the experiment it needs is not a new one:**
measure co-oximetry %MetHb and directly measured P50 together in the same
samples, and look at the slope. If the slope is zero, this model's central
claim is wrong.

Other limitations honestly noted:

1. `E660M` and `E940M` were set **equal**. What the model actually predicts is
   the **shape** of the floor (an asymptotic floor, collapsing sensitivity),
   not its exact height.
2. `EB5` is an **in-vivo flux fraction**, not an assay activity fraction.
   Congenital type I is usually reported at 10–20% residual activity, but
   reproducing the observed MetHb of 15–30% requires an in-vivo flux of around
   3.5%. This discrepancy is real and is not hidden.
3. The haemolysis term is a **two-compartment population approximation**, not a
   cell-age-structure model. It reproduces the direction and rough magnitude,
   but should not be used to predict transfusion thresholds.
4. With no R available in the container, the mrgsolve file underwent
   **equation verification, not compilation verification.**

---

## Files

| File | Contents |
|---|---|
| [`mhb_qsp_model.dot`](mhb_qsp_model.dot) · [`.svg`](mhb_qsp_model.svg) · [`.png`](mhb_qsp_model.png) | Mechanistic map — 153 nodes, 209 edges, 19 clusters |
| [`mhb_mrgsolve_model.R`](mhb_mrgsolve_model.R) | mrgsolve model — 42 ODEs, 139 parameters, 23 scenarios |
| [`mhb_shiny_app.R`](mhb_shiny_app.R) | Shiny dashboard — 13 tabs |
| [`mhb_reference_check.py`](../../../methemoglobinemia/mhb_reference_check.py) | Independent Python/scipy re-implementation and verification |
| [`mhb_reference_output.txt`](../../../methemoglobinemia/mhb_reference_output.txt) | Full output of the script above (the source of every figure in this README) |
| [`mhb_references.md`](mhb_references.md) | 171 references — every PMID actually looked up via the PubMed E-utilities |

### Reproduce

```bash
python3 mhb_reference_check.py --all                      # 42/42 anchors pass
python3 mhb_reference_check.py --params mhb_mrgsolve_model.R   # 139/139 parameters match
dot -Tsvg mhb_qsp_model.dot -o mhb_qsp_model.svg
dot -Tpng -Gdpi=150 mhb_qsp_model.dot -o mhb_qsp_model.png
Rscript -e 'shiny::runApp("mhb_shiny_app.R")'
```

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational and
research purposes.** It was built from published literature but has not been
independently validated or certified, and **must not be used directly for
clinical decision-making, prescribing, or regulatory submission.** In
particular, the dose of methylene blue, its cumulative ceiling, and whether to
use it in G6PD deficiency must always follow current clinical guidelines and
toxicological consultation.
