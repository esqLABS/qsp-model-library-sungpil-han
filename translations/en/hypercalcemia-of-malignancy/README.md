# Hypercalcaemia of Malignancy — QSP Model

> **In one sentence** — plasma calcium is not a stored quantity but **the residual of
> four fluxes**, and those fluxes turn on **three clocks** with entirely different
> timescales (kidney = hours, gut = a day, bone = days). The drain that is the kidney
> has a **ceiling**, and the hypercalcaemia itself drags that ceiling down — so the
> crisis is not a large number but a **bifurcation**, and what fluid does is not to
> remove calcium but **to put the ceiling back up**.

<p align="center">
  <a href="../../../hypercalcemia-of-malignancy/mah_qsp_model.svg"><img src="../../../hypercalcemia-of-malignancy/mah_qsp_model.png" width="820" alt="Hypercalcaemia of malignancy QSP mechanistic map"></a><br>
  <sub>140 nodes · 18 subgraph clusters · 222 edges — click the <a href="../../../hypercalcemia-of-malignancy/mah_qsp_model.svg">SVG</a> for a larger view</sub>
</p>

---

## 1. Three structural theses

### THESIS 1 — calcium is a **residual**, and the four terms sit on **different clocks**

```
    dCa_ECF/dt  =  J_res   −  J_form  +  J_gut   −  U_Ca
                     |          |          |         |
      clock 3 (bone, 2–5 d) ───┘          |         └── clock 1 (kidney, hours)
            clock 2 (gut, one day) ───────┘
```

The ECF calcium pool is only **36 mmol**, while the glomerulus filters **259 mmol** a
day — the equivalent of turning it over completely seven times a day. Every treatment
moves **exactly one** of these four terms, and **there is no drug that is both fast
and fundamental.**

| Treatment | Term it moves | Onset | Duration |
|---|---|---|---|
| Normal saline | clock 1 (kidney) | hours | only while it is running |
| Calcitonin | clocks 1 + 3 | hours | **48 hours** (receptor down-regulation) |
| Zoledronic acid | clock 3 (bone) | 2–4 days | about 4 weeks |
| Denosumab | clock 3 (bone) | 4–10 days | about 4 weeks, independent of renal function |
| Glucocorticoid | clock 2 (+3) | 2–4 days | only in the calcitriol mechanism |
| Antineoplastic therapy | **removes the input itself** | weeks | the only one that keeps improving |

The guideline combination of calcitonin plus a bisphosphonate is **not additive
chemistry but a relay between two clocks**. The model reproduces that handover
interval exactly (note C):

| Time after treatment | Saline | Saline+ZA | Saline+CT+ZA | **Calcitonin's share** |
|---|---|---|---|---|
| 6 h | 3.350 | 3.307 | 3.120 | **0.188** |
| 12 h | 3.218 | 3.114 | 2.946 | 0.167 |
| 48 h | 3.228 | 2.898 | 2.831 | 0.067 |
| 96 h | 3.215 | 2.728 | 2.703 | 0.025 |
| 7 d | 3.233 | 2.641 | 2.634 | **0.007** |

Over the same interval, calcitonin receptor availability `R_CT` falls 1.00 → 0.51
(12 h) → 0.26 (24 h) → **0.11 (48 h)** and then recovers to 0.45 by day 7. **The whole
of the benefit is spent in the first 48 hours, and those 48 hours are precisely the
interval in which the bisphosphonate is not yet working.** This complementarity was
not fitted — `CT_KOUT` was set from the duration of the calcitonin response and
`ZOL_KUP` from the time to onset of the bisphosphonate, independently of each other.

### THESIS 2 — the drain has a **ceiling**, and the hypercalcaemia **drags that ceiling down**

Hold total calcium fixed, equilibrate everything else (volume, AQP2, symptoms, GFR,
bone, parathyroid), then read the steady-state urinary calcium: the curve **does not
saturate, it turns over and comes back down** (note A):

| Ca held at | Urinary Ca | GFR | ECF volume | FE_Ca | Urine output |
|---|---|---|---|---|---|
| 2.40 | 5.20 | 180.0 | 1.000 | 2.01% | 2.01 L/day |
| 2.80 | 16.80 | 179.1 | 0.995 | 5.60% | 1.82 |
| 3.00 | 22.40 | 177.7 | 0.988 | 7.05% | 1.52 |
| **3.20** | **27.16** | 175.2 | 0.976 | 8.17% | 1.15 |
| 3.40 | 13.04 | 120.1 | 0.838 | 5.43% | 0.84 |
| 3.60 | 2.22 | 60.8 | 0.733 | 1.75% | 0.60 |
| 3.80 | 0.00 | 28.2 | 0.653 | 0 | 0.49 |

**Between 3.2 and 3.6 the kidney loses 92% of its excretory capacity.** And the cause
of that is the very calcium it was trying to excrete (CaSR → AQP2 down-regulation →
nephrogenic diabetes insipidus, and nausea/vomiting → ECF contraction → enhanced
proximal tubular Na⁺/Ca²⁺ reabsorption + a fall in GFR → a reduced filtered load →
calcium rises).

Sweeping the unregulated calcium input (`J_EXO`) upwards and finding by bisection the
point at which the steady state disappears yields a **saddle-node bifurcation**
(note B):

| Condition | Bifurcation point (mmol/day) | Ca at that point | What it fails on |
|---|---|---|---|
| **Baseline (no supportive care)** | **23.07** | 3.008 | circulatory collapse |
| + saline 1.2 L/day | 40.90 | 3.549 | renal excretory capacity |
| + saline 2.4 L/day | 51.43 | 4.964 | renal excretory capacity |
| Distal Tm +60% (the renal action of PTHrP alone) | **10.92** | 3.020 | circulatory collapse |
| Distal Tm +60% + saline 2.4 L/day | 39.75 | 4.972 | renal excretory capacity |
| Thiazide | 20.04 | 3.011 | circulatory collapse |
| AQP2 suppression removed | 25.62 | 3.073 | renal excretory capacity |
| Nausea/vomiting arm removed | 31.84 | 3.202 | renal excretory capacity |
| Residual nephrons 60% / 30% | 12.89 / 7.25 | 2.777 / 2.651 | circulatory collapse |

Four things can be read from this table.

1. **Saline at only 1.2 L/day nearly doubles the tolerable input, and changes the
   failure mode itself** — from circulatory collapse to a simple exceedance of renal
   excretory capacity. That sentence is why fluid is first-line treatment.
2. **The renal action of PTHrP on its own** (without adding a single extra osteoclast
   anywhere in the model) drops the bifurcation point from 23.07 to 10.92, **less than
   half.** A humoral tumour tips the same patient over with less than half the bone
   resorption of an osteolytic tumour.
3. A thiazide shaves **13%** off the bifurcation point — the quantitative version of
   "avoid it".
4. Deleting the concentrating defect moves the bifurcation point by 11%, and
   **deleting the nausea/vomiting arm moves it by 38%.** That means the dominant
   element in this loop is not the renal side but **the gastrointestinal side**, which
   is a testable prediction, and it differs from what the conventional account
   emphasises.

### THESIS 3 — the four mechanisms are not four names for the same disease. **They load different arms.**

Immediately before treatment (day 12), **with total calcium identical**, the
physiology is entirely different (note D):

| | Humoral (PTHrP) | Local osteolytic | Ectopic PTH | Calcitriol |
|---|---|---|---|---|
| Total calcium | 3.550 | 3.550 | 3.547 | 2.835 |
| Ionised calcium | 1.716 | 1.715 | 1.713 | 1.404 |
| PTH (pmol/L) | 1.20 | 1.20 | **10.99** | 2.44 |
| PTHrP (pmol/L) | **14.28** | 0.50 | 0.50 | 0.50 |
| 1,25D (pmol/L) | 89 | 64 | 113 | **415** |
| Phosphate (mmol/L) | **0.86** | 1.38 | 0.91 | 1.19 |
| **Bone resorption J_res** | **20.5** | **34.3** | 18.6 | 13.3 |
| Intestinal absorption J_gut | 0.26 | −0.58 | 1.09 | **12.73** |
| **Urinary calcium U_Ca** | **5.9** | **17.0** | 5.2 | 16.9 |
| FE_Ca | **2.2%** | **6.5%** | 2.1% | 5.6% |
| Distal calcium Tm | **43.3** | **26.3** | 41.4 | 30.4 |
| CTX (ng/mL) | 0.897 | 1.502 | 0.817 | 0.582 |

Read the bone resorption row and the urinary excretion row alongside each other. **The
osteolytic patient has to dissolve 67% more bone than the humoral patient to reach the
same calcium** — because PTHrP opens the tap and closes the drain at the same time
(distal Tm 26.3 → 43.3, +65%). One ligand does two things in the same direction.

Look too at what the intestinal arm is doing in the other three columns: 0.26,
**−0.58**, 1.09 mmol/day. Anorexia and suppressed 1,25D have already closed that arm,
and in myeloma it is frankly negative. "Stop the calcium supplements" is almost always
right and almost always insufficient.

**The treatment results follow with no further assumptions:**

| Scenario | Nadir Ca | Days held below Ca 2.70 | CTX nadir |
|---|---|---|---|
| Humoral + saline + ZA (s05) | 2.635 | 15.0 | 0.153 |
| Osteolytic + saline + ZA (s11) | **2.402** | **25.9** | 0.302 |
| Calcitriol + ZA (s12) | 2.599 | 28.5 | **0.109** |
| Calcitriol + prednisolone (s13) | **2.411** | 22.4 | 0.382 |

The last two rows are the crux. In the calcitriol mechanism, prednisolone reaches a
**deeper and faster** nadir than zoledronic acid (2.411 vs 2.599, 0.85 days to
normalisation vs 1.50 days) while suppressing bone resorption **three times less**
(CTX 0.382 vs 0.109). The two drugs are working in different arms, and **the model can
say which arm is carrying the disease.**

---

## 2. One agreement that is not a coincidence (an internal consistency check)

It was not fitted. The net influx at day 12 (J_res − J_form + J_gut) is **13.18** for
the humoral case and **24.23** mmol/day for the osteolytic one — nearly a twofold
difference for the same plasma calcium. And the bifurcation points of note B,
computed entirely separately, are **10.92** for the humoral renal phenotype and
**23.07** for a kidney PTHrP has not touched. **Both mechanisms sit just past their own
bifurcation point, by a similar relative margin.** That is why two such different
physiologies arrive at the same number. The two calculations share parameters, but
neither was adjusted to fit the other.

---

## 3. The self-delivering drug

It is not the bloodstream that delivers zoledronic acid to the osteoclast. It binds to
hydroxyapatite and is then endocytosed by the osteoclast **while resorption is going
on**, so the rate of uptake is **proportional to the very flux the drug is trying to
abolish**:

```
    uptake = ZOL_KUP × (drug on the bone surface) × (J_res / J_res_normal)
```

Three observations come out of that one term at once (notes E/F):

| J_res at dosing | Peak intracellular concentration | Day of the peak | CTX nadir/baseline |
|---|---|---|---|
| 8.0 (normal turnover) | 0.418 | 8.40 | 0.190 |
| 20.5 (humoral) | 0.610 | 6.15 | 0.170 |
| 34.3 (osteolytic) | 0.725 | **4.85** | 0.201 |

**The sicker the bone, the faster and the more of the drug loads itself.** And when
resorption stops, delivery stops too; what is left on the surface is buried in the
matrix (`ZOL_BURY`) and the intracellular drug disappears with a 15-day half-life —
which is why a single dose lasts four weeks. The same term is also why the
dose-response flattens at the top:

| Dose | Nadir Ca | Days normocalcaemic out of 30 |
|---|---|---|
| 1 mg | 2.763 | 0.0 |
| 2 mg | 2.686 | 5.6 |
| **4 mg** | **2.635** | **15.0** |
| 8 mg | 2.604 | 23.2 |
| 16 mg | 2.587 | 26.4 |

Going from 4 to 8 mg buys 8.2 days; going from 8 to 16 mg buys only 3.2. What creates
the ceiling is not target occupancy but **the delivery step**.

---

## 4. Renal failure inverts the choice of drug (note H)

Starting from 30% residual nephrons, the same tumour presents not at 3.55 but at
**4.40 mmol/L** (because the bifurcation point has come down from 23.1 to 7.25).
Against that background:

| | Nadir Ca | Days normocalcaemic | Outcome |
|---|---|---|---|
| Zoledronic acid (s09) | 2.967 | 0.0 | reaches the collapse boundary on day 46 |
| Denosumab (s10) | **2.621** | **7.7** | no collapse |

Denosumab has **no renal step anywhere in its disposition route**, and it does not
require bone resorption in order to deliver itself. It is the only comparison in this
file in which the two antiresorptives part company decisively, and the model reaches
that conclusion from **the structure of the disposition** rather than from an assumed
difference in potency.

---

## 5. Loop diuretics, when the question is put properly (note I)

It is not "furosemide is harmful". The story is more specific than that.

| Scenario | Nadir Ca | Days normocalcaemic | Outcome |
|---|---|---|---|
| s05 adequate saline + ZA | 2.635 | 15.0 | — |
| s15 adequate saline + ZA + **furosemide** | **2.254** | **19.2** | benefit |
| s23 inadequate saline (1.2 L/day) + ZA | 2.642 | 14.0 | — |
| s16 inadequate saline + ZA + **furosemide** | 3.453 | 0.0 | **collapse boundary on day 2.6** |

Furosemide does in fact increase calcium excretion, and the model says so. At the same
time, if the volume underneath is not replenished, it can turn a recoverable illness
into **circulatory collapse within three days.** **The variable that decides which of
the two happens is not the diuretic.** Suki's original 1970 report was done with
thorough volume repletion, and LeGrand's 2008 systematic review criticised the practice
where it is not — and the model places both papers in the same table without
contradiction.

---

## 6. Where corrected calcium fails (note G)

"The true equivalent" = the total calcium that would give the same ionised calcium at
an albumin of 40 g/L.

| Albumin | Total calcium | Ionised | Payne correction | True equivalent | **Payne error** |
|---|---|---|---|---|---|
| 35 | 3.55 | 1.879 | 3.65 | 3.76 | −0.11 |
| 30 | 3.55 | 1.997 | 3.75 | 3.99 | −0.24 |
| 25 | 3.55 | 2.130 | 3.85 | 4.26 | **−0.41** |
| 22 | 3.55 | 2.219 | 3.91 | 4.44 | **−0.53** |
| 25 | 2.40 | 1.440 | 2.70 | 2.88 | −0.18 |

The error is not a constant offset: **it grows with both the hypoalbuminaemia and the
calcium value.** The patient in whom this correction is least trustworthy is exactly
the patient it was created for. Scenario 18 is the dynamic version — give the same
tumour to a patient with an albumin of 22 g/L and **the total calcium is lower
(3.42 vs 3.55) while the ionised calcium is much higher (2.06 vs 1.72)**.

pH acts on the same scale: at total calcium 3.00 and albumin 40, the ionised calcium is
1.563 at pH 7.30 and 1.437 at pH 7.50. So **the alkalosis of vomiting is a genuine
negative feedback arm inside the vicious circle** — it destroys volume while protecting
the ionised calcium.

---

## 7. Deliverables

| File | Contents |
|---|---|
| [`mah_qsp_model.dot`](../../../hypercalcemia-of-malignancy/mah_qsp_model.dot) · [`.svg`](../../../hypercalcemia-of-malignancy/mah_qsp_model.svg) · [`.png`](../../../hypercalcemia-of-malignancy/mah_qsp_model.png) | Mechanistic map — 140 nodes, 18 clusters, 222 edges |
| [`mah_mrgsolve_model.R`](mah_mrgsolve_model.R) | 50-ODE mrgsolve model, 205 parameters, 23 scenarios, calibration notes A–J |
| [`mah_shiny_app.R`](../../../hypercalcemia-of-malignancy/mah_shiny_app.R) | 10-tab Shiny dashboard (including the ceiling curve and the bifurcation calculation) |
| [`mah_references.md`](mah_references.md) | 105 references, every PMID checked with the NCBI E-utilities |

### The scenario list

| # | Scenario | What it shows |
|---|---|---|
| 1 | Healthy control (90 days) | an exact steady state (every variable drifting under 0.2%) |
| 2 | Humoral, untreated | the natural course → the collapse boundary |
| 3 | + saline alone | how far one fast arm on its own can get |
| 4 | + saline + calcitonin | the first 48 hours, and the tachyphylaxis |
| 5 | + saline + zoledronic acid 4 mg | the standard |
| 6 | + calcitonin + zoledronic acid | the guideline triple therapy — the relay interval |
| 7 | Zoledronic acid alone (no saline) | why the slow arm alone is not enough |
| 8 | + saline + denosumab | deeper and for longer |
| 9 / 10 | ZA vs denosumab in renal failure (30%) | the point at which the choice of drug inverts |
| 11 | Osteolytic + ZA (calcium-matched) | the key comparison for THESIS 3 |
| 12 / 13 | ZA vs prednisolone in lymphoma | which arm is carrying the disease |
| 14 | Ectopic PTH + cinacalcet + ZA | the PTH-dependent mechanism |
| 15 / 16 / 23 | Furosemide × degree of volume repletion | the resolution of the diuretic controversy |
| 17 | Haemodialysis (low-calcium dialysate) | the only means of removing calcium faster than bone |
| 18 | Albumin 22 g/L | the divergence between total and ionised calcium |
| 19 / 20 | Immobilisation · thiazide | quantifying the precipitating factors |
| 21 | + effective antineoplastic therapy | removing the input rather than buffering it |
| 22 | Zoledronic acid q28d × 2 | the redosing interval |

---

## 8. Validation

- Every ODE was independently implemented in **Python/scipy (LSODA, rtol 1e-7)** for
  the calibration.
- The `$ODE` block of the R file was then **mechanically extracted**, compiled with
  `g++ -Wall` (passing with no warnings, and all 50 derivatives confirmed assigned) and
  integrated with **fourth-order Runge-Kutta at a step of 0.002 days**.
- Across four scenarios × 60 days the two implementations agree **to within 0.0005% on
  every state variable and every derived quantity**.
- The healthy control is **an exact steady state by construction** — `KL_BASE`,
  `KO_BASE`, `PTG_SYN`, `K1A`, `D25_IN`, `MG_FE`, `ORAL_H2O` and `GUT_SEC` were each
  **solved analytically** to make the corresponding derivative zero, not fitted. Over
  90 days the drift is +0.004% in total calcium, +0.10% in PTH, −0.18% in 1,25D,
  −0.03% in osteoclast number, −0.07% in free RANKL and −0.03% in CTX.
- The reference values reproduced: total calcium 2.400 · ionised 1.1995 · filtered load
  259 mmol/day · urinary calcium 5.21 mmol/day · FE_Ca 2.01% · GFR 180 L/day · urine
  output 2.02 L/day · bone turnover 8.0 mmol/day in both directions ·
  net intestinal absorption 5.22 mmol/day · segmental reabsorption proximal 65% /
  TAL 25% / distal 8%.

---

## 9. What is not in this model (note J — honestly)

There is **no** infection, thrombosis, hypercoagulability, analgesia, geometry of bone
metastases, skeletal-related events (SRE), osteonecrosis of the jaw, acute-phase
reaction, atrial fibrillation, digoxin interaction or mortality model. The "collapse
boundary" (an ECF deficit above 38% or a total calcium above 5.0 mmol/L) is **the edge
of the calibrated domain, not a prediction of death.** The tumour is a single logistic
compartment with no resistance, no heterogeneity and no metastatic dynamics, and the
antineoplastic therapy of scenario 21 should be read as "the input has been removed"
rather than as the simulation of a particular regimen. Bone is a single remodelling
compartment and does not distinguish trabecular from cortical bone, so **it has nothing
to say about fracture risk.** Every parameter is a population point estimate with no
individual variability — this model reproduces the mechanisms and their order, not
individual patients.

The **modelling choices** that are not supported by the literature are set out
separately in [`mah_references.md` §17](mah_references.md).
Of those, the one that is at once most uncertain and most influential on the results is
the **vomiting arm (`VOM_MAX`)** — by note B the bifurcation point is far more sensitive
to this arm (38%) than to the concentrating defect (11%), which makes it both the
parameter most worth measuring and the parameter the model is least confident about.

---

## 10. Running it

```r
# the model
source("mah_mrgsolve_model.R")
out <- run_scenario("s05")                 # humoral + saline + zoledronic acid
plot(mah, out, "Ca_tot,iCa,J_res,U_Ca,CTX,vol")

# the dashboard
shiny::runApp("mah_shiny_app.R")
```

```bash
# redrawing the map
dot -Tsvg mah_qsp_model.dot -o mah_qsp_model.svg
dot -Tpng -Gdpi=150 -Gsize="80,50" mah_qsp_model.dot -o mah_qsp_model.png
```

---

> ⚠️ This is a **QSP model for educational and research purposes**. It was built from
> the public literature but has not been independently validated or certified, and must
> not be used directly for real clinical decision-making, prescribing or regulatory
> submission.
