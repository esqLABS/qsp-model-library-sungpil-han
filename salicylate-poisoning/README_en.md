# Salicylate (Aspirin) Poisoning — Quantitative Systems Pharmacology Model

> Between the number the lab reports and the number that kills the patient sit
> **two multiplications that nobody measures**.
> This model computes those two multipliers explicitly and checks whether the
> clinical paradoxes of salicylate poisoning follow arithmetically from them.

---

## 1. The thesis

Salicylate poisoning is classified — and dialysis decided — by a single
number, **total plasma salicylate concentration**. But what actually kills the
patient is **salicylate inside the brain cell**, and between the two sit two
unmeasured multiplications.

```
C_brain  =  C_total  ×  fu(C_total, pH)  ×  f_n(pH_plasma) / f_n(pH_brain)
            ─────────   ────────────────    ─────────────────────────────
            the only     ① albumin           ② pH partitioning of a weak
            value the    saturation          acid (pKa 3.0) — brain
            lab reports  150 mg/L → 8.5%      intracellular pH is buffered,
                         800 mg/L → 43%       plasma pH is not
```

② is the core of this model, and it needs physiology that ordinary PK models
leave out. **CO₂ crosses the cell membrane instantly and is absorbed inside
the cell by a non-bicarbonate buffer system (β ≈ 35 mmol/L per pH).** So
raising PaCO₂ from 25 → 60 mmHg gives:

| | Change |
|---|---|
| Plasma pH | falls **0.380** |
| Brain intracellular pH | falls **0.143** (moves 2.7-fold less) |
| Brain:plasma partition coefficient | rises **1.73-fold** |

That is, **the gradient itself changes.** Acute respiratory acidosis and
metabolic acidosis are never equivalent for the brain.

---

## 2. Three results the model derives (derived, not assumed)

### (1) The ventilator is a narrow-therapeutic-index drug — this model's central experiment

Eight hours after a 30 g ingestion, **a single switch** is flipped. Arm A is
intubated at a conventional minute ventilation; arm B is intubated matched to
the patient's own minute ventilation (2.30-fold normal). Everything else is
identical.

| Time | Plasma A | Plasma B | Brain A | Brain B | Brain A/B |
|---|---|---|---|---|---|
| 8 h (just before the switch) | 813.8 | 813.8 | 143.3 | 143.3 | 1.00 |
| 9 h | **581.6** | 811.6 | **162.3** | 144.5 | 1.12 |
| 12 h | 577.5 | 796.4 | 158.1 | 141.3 | 1.12 |
| 24 h | 547.1 | 692.0 | 128.7 | 106.6 | 1.21 |
| 36 h | 496.7 | 569.6 | 98.6 | 71.9 | 1.37 |

Within one hour, pH goes 7.47 → 7.11, PaCO₂ 25 → 57 mmHg, brain:plasma ratio
0.178 → 0.279 (1.57-fold). **The lab value falls 28% while the brain
concentration rises.** Cumulative 36-hour CNS exposure is 1.31-fold.

A milder form of the same lesion is **co-ingested opioid** — suppressing
respiratory drive by only 45% is enough to raise brain concentration
(scenario S13).

### (2) Sodium bicarbonate is the only agent that touches both multipliers at once

| Intervention | Urine pH (12 h) | Renal clearance | Plasma at 24 h | Brain at 24 h |
|---|---|---|---|---|
| Supportive care only (fluids) | 6.10 | 3.5 mL/min | 670 | 109 |
| Equal-volume saline diuresis | 6.10 | 16.6 | 409 | 43 |
| **NaHCO₃ + K** | **7.54** | **38.5** | **191** | **7.7** |

- **The plasma-pH arm**: brain concentration falls within 40 minutes of
  dosing, before the urine has even responded.
- **The urine-pH arm**: at urine pH 5, 97% of filtered drug is reabsorbed by
  non-ionic diffusion, but at urine pH 7.5–8 that falls below 3%. An 11-fold
  difference in clearance.
- **Acetazolamide** alkalinises the urine while **acidifying the blood** — the
  worst possible combination available in this model.

### (3) The rate-limiting reagent of this treatment is potassium

Urine alkalinises only once plasma bicarbonate crosses the **renal tubular
bicarbonate threshold**, and hypokalaemia raises that threshold (+5.5 mmol/L
per 1 mmol/L of K deficit).

| | K repleted | K not repleted |
|---|---|---|
| Time to reach urine pH 7.5 | 11.7 h | **21.9 h** |
| Plasma HCO₃ needed at that point | — | **+2.9 mmol/L more** |
| Renal clearance at 12 h | 38.5 mL/min | 17.9 |
| Brain concentration at 48 h | reference | **2.97-fold** |

Withholding potassium does not make alkalinisation *impossible — it delays it
by 10 hours*, at the cost of tolerating more severe alkalaemia to get there —
the model draws this distinction precisely as "delay and cost," not "failure."

---

## 3. Same number, different patient (why a single level means nothing)

| Situation | Plasma mg/L | fu | Arterial pH | Brain pHi | Brain mg/L | Brain/plasma |
|---|---|---|---|---|---|---|
| Acute 30 g, 6 h | 809 | 0.43 | 7.51 | 7.11 | 133.0 | 0.164 |
| Acute 30 g, 30 h | 611 | 0.35 | 7.41 | 7.04 | 92.6 | 0.151 |
| Intubated at conventional ventilation | 582 | 0.40 | 7.11 | 6.96 | **162.3** | **0.279** |
| **Chronic, elderly, day 14** | **216** | 0.33 | 6.75 | 6.92 | **109.2** | **0.505** |

A chronic patient's plasma concentration of 216 mg/L (21.6 mg/dL) is commonly
read as "mild–moderate." But because it is at equilibrium, albumin is low,
and acidaemia is present, brain concentration is 109 mg/L — matching an acute
massive-ingestion patient. An acute patient who first reaches the same plasma
concentration (0.9 hours post-ingestion) has a brain concentration of only
3.0 mg/L — **a 36-fold difference in the brain:plasma ratio.**

This is **why the Done nomogram fails**. The model does not assume the
nomogram; it treats it as *output* and scores it (Verification report,
Section 8): in 4 of 6 situations, the nomogram's grade and the model's
computed brain-concentration grade disagree.

---

## 4. Deliverables

| File | Contents |
|---|---|
| `sal_qsp_model.dot` / `.svg` / `.png` | Mechanistic map — **125 nodes, 16 clusters**. Exposure · absorption · protein binding · pH partitioning · hepatic metabolism · renal handling · uncoupling · acid-base · CNS · other organs · electrolytes · treatment · iatrogenic harm · measurement · clinical endpoints |
| `sal_mrgsolve_model_en.R` | mrgsolve model with **28 ODEs, 14 scenarios**. Includes a stepwise parameter-change runner |
| `sal_shiny_app_en.R` | **11-tab** dashboard (patient · PK · the two multipliers · acid-base · kidney · endpoints · scenario comparison · ventilator experiment · nomogram · EXTRIP · lab panel) |
| `sal_verify_python.py` | **Independent Python/scipy re-implementation** of the 28 ODEs |
| `sal_verification_output.txt` | Verification report — **54 of 54** literature-based anchors pass |
| `sal_references_en.md` | **83 references**, every PMID confirmed against PubMed |

---

## 5. Model structure (28 ODEs)

| Group | State variables |
|---|---|
| GI tract | `AST` stomach · `AGUT` small intestine · `ACONC` bezoar · `AAC` activated charcoal |
| Systemic | `AASA` unhydrolysed aspirin · `ACENT` central · `APER` tissue · `ACNS` brain |
| Metabolism·excretion | `ASU` salicyluric acid · `APG`/`AAG` glucuronide conjugates · `AGA` gentisic acid · `AUR` urinary excretion · `AHD` dialysis removal · `GLY` glycine pool |
| Acid-base | `HCO3` · `PACO2` · `BBB` brain buffer base · `BBT` muscle buffer base · `LAC` · `KET` · `FATIG` respiratory-muscle fatigue |
| Electrolytes·volume | `KBAL` potassium balance · `VECF` extracellular fluid |
| Toxicodynamics | `TCORE` core temperature · `GLUB` brain glucose · `LUNG` lung injury · `CNSI` cumulative CNS exposure |

The three key algebraic equations are implemented with the **identical
algorithm** in mrgsolve's `$GLOBAL` and in Python.

1. `freeSal()` — Newton solution for two-site albumin binding
2. `phIntra()` — Newton solution for intracellular pH from buffer base and
   non-bicarbonate buffering capacity
3. `phUrine()` — urine pH from urinary bicarbonate concentration

---

## 6. What was fitted and what was predicted

**Fitted — matched only to normal or therapeutic-dose data**

- Two-site albumin binding (fu 8.5% @150, 29% @500, 43% @800 mg/L)
- Levy's dose-limited metabolism (salicyluric acid Vmax/Km, phenolic
  glucuronide) → single 650 mg dose t½ 3.6 h, steady state 202 mg/L at
  3.9 g/day
- Tubular non-ionic back-diffusion constant (unchanged fraction excreted in
  acidic urine at therapeutic doses, 1–6%)
- Intracellular buffering capacity (³¹P-MRS observation that brain pHi is
  defended against acute PaCO₂ changes)

**Predicted — nothing here was fitted to it**

- t½ of 37 h and apparent Vd rising 0.11 → 0.32 L/kg in overdose
- Respiratory alkalosis first, then a mixed disorder (pH 7.51, PaCO₂ 24 mmHg
  at 6 hours)
- The intubation catastrophe, and the accompanying phenomenon in which
  **plasma concentration actually falls**
- The 10-hour potassium delay
- Hypoglycorrhachia — brain glucose falling to 0.43 mmol/L while blood glucose
  is normal
- The Done nomogram's misclassification pattern

---

## 7. Defects the verification found

The Python re-implementation was not decorative — it **caught four genuine
errors**, all of which were fixed.

1. **Brain intracellular pH was being computed from a fixed intracellular
   bicarbonate.** That let acute hypercapnia drop brain pH by exactly as much
   as plasma pH, **erasing the partitioning effect entirely.** Only after
   switching to a buffer-base-plus-buffering-capacity approach did this
   model's central experiment hold together.
2. **Single-site albumin binding** saturated to a free fraction of 64% by
   800 mg/L, leaving no room for acidaemia to displace further. Refitted to
   two binding sites.
3. **No maintenance fluid intake term** meant every long-duration simulation
   ended in a spurious contraction acidosis with GFR pinned to its floor
   (making every chronic scenario meaningless).
4. **Untitrated bicarbonate infusion** drove plasma HCO₃ up to 45 mmol/L,
   completely masking the potassium failure mode. Adding the clinical rule of
   stopping infusion at an arterial pH ceiling (7.60) let the potassium
   dependency show through properly.

---

## 8. Limitations (stated, not hidden)

- **Brain concentration, the core variable of this model, has never been
  measured in a living human.** It was computed from animal partitioning
  experiments (Hill 1971 · PMID 5089754) and the pH physics of a weak acid,
  not fitted to clinical outcomes.
- The falsifiable prediction is **direction and magnitude**: immediately after
  a PaCO₂ change, plasma (−28%) and brain (+12%) move in opposite directions,
  and by 24 hours the brain is +21%.
- The clearance benefit of urinary alkalinisation over supportive care
  (11-fold) is sensitive to the comparator (urine pH · urine output of the
  supportive-care patients). The "10–20-fold" figures in the literature are
  usually comparisons against dehydrated · acidic-urine patients.
- Renal tubular secretion (OAT1/3) is present on the map but is absorbed into
  filtration · back-diffusion in the ODEs.
- Paediatric parameters (beyond weight scaling), pregnancy, and clearance
  differences between dialysis membrane types are not addressed.

---

## 9. Running it

```bash
# Render the mechanistic map
dot -Tsvg sal_qsp_model.dot -o sal_qsp_model.svg
dot -Tpng -Gdpi=150 sal_qsp_model.dot -o sal_qsp_model.png

# Independent verification (Python)
pip install numpy scipy
python3 sal_verify_python.py            # 54/54 pass report

# mrgsolve model
R -e "source('sal_mrgsolve_model_en.R'); print(summarise_scenario(run_scenario('S07')))"

# Shiny dashboard
R -e "shiny::runApp('sal_shiny_app_en.R', port = 8080)"
```

---

## 10. Scenarios (14 scenarios)

| ID | Contents |
|---|---|
| S01 | Therapeutic-dose aspirin 650 mg q4h × 5 days |
| S02 | Acute 30 g, supportive care (fluids) only |
| S03 | Acute 30 g, no fluids — renal clearance collapses |
| S04 | Acute 30 g + saline diuresis 250 mL/h |
| S05 | Acute 30 g + NaHCO₃ + potassium (standard treatment) |
| S06 | Acute 30 g + NaHCO₃, no potassium |
| **S07** | Acute 30 g, intubated at 8 hours **with conventional PaCO₂** |
| **S08** | Acute 30 g, intubated at 8 hours **matched to minute ventilation** |
| S09 | Acute 30 g + NaHCO₃, haemodialysis at 6 hours |
| S10 | Acute 30 g + NaHCO₃, dialysis delayed to 18 hours |
| S11 | Enteric-coated 30 g + 45% bezoar, multi-dose activated charcoal |
| S12 | Chronic salicylism: 1.3 g q6h, age 80, GFR 60, albumin 28 |
| S13 | Acute 30 g + co-ingested opioid (45% suppression of respiratory drive) |
| S14 | Acute 30 g + NaHCO₃ + glucose (blood glucose 11 mmol/L) |

---

## 11. References

`sal_references_en.md` — 83 entries, grouped into 18 sections. Every PMID was
confirmed by lookup against PubMed E-utilities, and the 2 that could not be
confirmed are separated out in Section 17.
`[CAL]` marks evidence used to *fit* a parameter, `[VAL]` marks evidence used
to *validate* it without fitting.

Five key references:

- Palmer BF, Clegg DJ. Salicylate Toxicity. *N Engl J Med* 2020 — [PMID 32579814](https://pubmed.ncbi.nlm.nih.gov/32579814/)
- Hill JB. Experimental salicylate poisoning: effects of altering blood pH on tissue and plasma salicylate. *Pediatrics* 1971 — [PMID 5089754](https://pubmed.ncbi.nlm.nih.gov/5089754/)
- Levy G. Pharmacokinetics of salicylate elimination in man. *J Pharm Sci* 1965 — [PMID 5862532](https://pubmed.ncbi.nlm.nih.gov/5862532/)
- Stolbach AI, et al. Mechanical ventilation was associated with acidemia in salicylate-poisoned patients. *Acad Emerg Med* 2008 — [PMID 18821862](https://pubmed.ncbi.nlm.nih.gov/18821862/)
- Juurlink DN, et al. EXTRIP recommendations for salicylate poisoning. *Ann Emerg Med* 2015 — [PMID 25986310](https://pubmed.ncbi.nlm.nih.gov/25986310/)

---

> ⚠️ **Disclaimer**: a QSP model for educational · research purposes. It has
> not been independently validated · certified and must not be used for
> actual clinical decision-making, prescribing, or regulatory submission.
