# Neonatal Hyperbilirubinaemia (Neonatal Jaundice) · Kernicterus — QSP Model

> **The value that is measured is not the value that does the damage.**
> And what joins the two is a saturable binding isotherm.

This directory holds a quantitative systems pharmacology (QSP) model of neonatal
hyperbilirubinaemia. It consists of 34 ODEs, a 137-node mechanistic map, 10
treatment scenarios, an 11-tab Shiny app, and 97 references verified against
PubMed.

| File | Contents |
|------|------|
| [`nhb_qsp_model.dot`](nhb_qsp_model.dot) | Mechanistic map source (15 clusters · 137 nodes · 179 edges) |
| [`nhb_qsp_model.svg`](nhb_qsp_model.svg) · [`nhb_qsp_model.png`](nhb_qsp_model.png) | Rendered map |
| [`nhb_mrgsolve_model.R`](nhb_mrgsolve_model.R) | mrgsolve 34-ODE model + 10 scenarios + 8 analysis functions |
| [`nhb_shiny_app.R`](../../../neonatal-hyperbilirubinemia/nhb_shiny_app.R) | 11-tab interactive dashboard |
| [`nhb_reference_check.py`](../../../neonatal-hyperbilirubinemia/nhb_reference_check.py) | **Independent implementation** (pure python RK4). The source of every number here |
| [`nhb_reference_output.txt`](../../../neonatal-hyperbilirubinemia/nhb_reference_output.txt) | Run output of the file above (A1-A14) |
| [`nhb_references.md`](nhb_references.md) | 97 references (all bibliographically verified through the PubMed API) |

```bash
dot -Tsvg nhb_qsp_model.dot -o nhb_qsp_model.svg
dot -Tpng -Gdpi=150 nhb_qsp_model.dot -o nhb_qsp_model.png
python3 nhb_reference_check.py > nhb_reference_output.txt   # ~100 s
Rscript -e 'shiny::runApp("nhb_shiny_app.R")'
```

---

## Why this disease was modelled this way

Neonatal jaundice is common (60-80 % of term infants are visibly jaundiced) and
mostly resolves on its own. Even so, kernicterus still occurs, and worldwide G6PD
deficiency is its largest single cause. This combination — "common and benign, but
rarely devastating" — has a particular character as a modelling target. **The risk
does not lie in the quantity we measure (total serum bilirubin, TSB) but in the way
that quantity is translated into another quantity (free bilirubin, Bf).**

This model does not put "bilirubin goes up and phototherapy brings it down" at its
centre. It puts the following isotherm at its centre.

```
TSB = Bf + A·K1·Bf/(1+K1·Bf) + A·n2·K2·Bf/(1+K2·Bf)
```

A is the molar albumin concentration, K1 the binding constant of the high-affinity
site (≈4.5×10⁷ M⁻¹) and K2 that of the low-affinity site (≈1×10⁶ M⁻¹). The
left-hand side is the number the laboratory reports, and every threshold in the
guidelines is written against it. Bf is the **only** chemical species that crosses
the blood-brain barrier.

Four things follow from this single line, and none of them was put into the model
as a rule — all four are **computed results**.

---

## Theme 1 · The AAP "neurotoxicity risk factor" correction is computed

Solving the isotherm backwards at Bf = 30 nM gives, in different babies, **the TSB
that means the same risk** (`nhb_reference_output.txt` A1b):

| Condition | Albumin (g/dL) | TSB at Bf=30 nM | Difference | B/A molar ratio |
|------|------|------|------|------|
| Term, pH 7.40, no displacing drug | 3.5 | 18.58 | (reference) | **0.604** |
| Albumin 3.0 (AAP risk factor) | 3.0 | 15.93 | **−2.65** | **0.604** |
| Albumin 2.5 | 2.5 | 13.27 | −5.31 | **0.604** |
| Acidosis pH 7.15 (K1 × 0.70) | 3.5 | 15.59 | −2.99 | 0.506 |
| Ceftriaxone-type displacement | 3.5 | 15.59 | −2.99 | 0.506 |
| Preterm 30 weeks (K1 × 0.75) | 2.8 | 12.93 | −5.65 | 0.525 |
| Sepsis + acidosis + albumin 2.6 | 2.6 | 10.59 | −7.99 | 0.463 |

The AAP 2022 guideline lowers the phototherapy threshold by **about 2 mg/dL** when
a neurotoxicity risk factor is present. From albumin 3.0 g/dL alone the model
yields **−2.65 mg/dL**. It has never been shown the guideline. The value comes out
of the binding stoichiometry.

### And why the B/A ratio is a good index, and where it fails

At a fixed Bf,

```
B/A = K1·Bf/(1+K1·Bf) + n2·K2·Bf/(1+K2·Bf)
```

**there is no albumin term in this expression**. The first three rows of the table
above differ in albumin and differ in TSB, yet B/A is exactly the same 0.604. That
is the mathematical content of the claim that "B/A is better than TSB".

But K1 does appear in the expression. In the acidosis and displacement rows, at the
same Bf = 30 nM, B/A falls to 0.506. **B/A is robust to problems of binding
capacity and helpless against problems of binding affinity.** That is, in the baby
with sepsis, acidosis and ceftriaxone — precisely the baby at greatest risk — B/A
is wrong in the reassuring direction. No amount of diligent albumin measurement
recovers this error.

---

## Theme 2 · TSB is the integral of a difference of two fluxes, and ETCOc identifies it

Because `d(burden)/dt = production − removal`, a baby with double production and a
baby with half the removal can trace **the same TSB curve**. Yet the two babies
respond to entirely different drugs.

HO-1 releases exactly one molecule of CO per molecule of bilirubin. Carboxyhaemoglobin
is therefore a direct read-out of the production flux (A3, A8):

| State | ETCOc (ppm) | Hb change |
|------|------|------|
| Physiological jaundice | 1.38 | none |
| ABO alloimmunisation | 2.48 | slight |
| Rh(D) alloimmunisation | 3.37 | 13.5 → 14.0 (d7, after treatment) |
| G6PD oxidative crisis (24 h) | 9.95 | 17.0 → 15.8 |
| UGT1A1*6/*6 (removal lesion) | **1.38** | none |

The last row is the point. The `*6/*6` baby reaches a TSB of 12.8 while ETCOc is
entirely normal. A production lesion and a removal lesion are separable at the
bedside, and the model even tells you which drug can work. Stannsoporfin and IVIG
enter only the production term, phenobarbital and gene therapy only the removal
term, and **phototherapy is in neither of them** (Theme 4).

---

## Theme 3 · The assay counts photoisomers as bilirubin

Phototherapy is not one reaction but two.

1. **Reversible** (4Z,15Z) ⇄ (4Z,15E) configurational isomerisation — fast, and it
   reaches a photostationary state.
2. **Irreversible** cyclisation of the E-isomer → lumirubin — this is what actually
   carries the excretion, and it does not require conjugation (t½ ≈ 1 h).

Both photoisomers are measured as "bilirubin" by the usual assays. The result (A5):

| Time | Reported TSB | Native (4Z,15Z) | Lumirubin | E-isomer | Photoisomer % |
|------|------|------|------|------|------|
| 24 h (just before the lamps go on) | 7.83 | 7.81 | 0.000 | 0.012 | 0.2 |
| 30 h | 8.93 | 6.82 | 0.409 | 1.703 | 23.6 |
| 48 h | 7.28 | 5.55 | 0.351 | 1.377 | 23.7 |
| 72 h (lamps off) | 6.06 | 4.64 | 0.286 | 1.130 | 23.4 |
| **74 h** | **5.72** | 4.58 | 0.071 | 1.066 | 19.9 |
| 84 h | 6.80 | 6.01 | 0.000 | 0.792 | 11.6 |
| 96 h | 8.07 | 7.51 | 0.000 | 0.555 | 6.9 |

This is an ABO-alloimmunised newborn given intensive phototherapy from 24 to 72
hours. There are three points.

- During phototherapy, **26 % of the reported TSB is photoisomer**. The
  photostationary-state E-isomer is 20.6 %, in agreement with the reported
  literature values (20-25 %).
- **The reported TSB rises in the first few hours after the lamps go on**
  (7.83 → 8.93). The native pigment is already coming down. The whole of the rise
  is the photostationary state being established.
- When the lamps go off, TSB **falls first** (5.72 at 74 h), because lumirubin
  washes out with a t½ of 1 h. Only after that does it climb again through thermal
  E→Z back-conversion and continuing production (+2.00 mg/dL over 24 hours). **Part
  of the "rebound" is not newly made pigment but photochemical bookkeeping.**

---

## Theme 4 · Phototherapy is the only route that bypasses UGT1A1, and the reason it fails as the child grows is not the usual explanation

Lumirubin does not require conjugation. That is why phototherapy works even in
Crigler-Najjar type I, where UGT1A1 activity is zero. And in the same model the
phenobarbital arm **overlaps completely** with the untreated arm (A11a):

| Condition | d3 | d7 | d14 |
|------|------|------|------|
| CN-I, no treatment | 9.65 | 17.31 | 25.71 |
| CN-I, continuous intensive phototherapy | 3.79 | 3.46 | 3.27 |
| CN-I, phototherapy + phenobarbital | 3.79 | 3.46 | **3.27** |
| CN-II, no treatment | 9.55 | 16.68 | 23.40 |
| CN-II, phenobarbital | 9.48 | 15.98 | **20.28** |

This is not a rule. Induction enters as `GENO × ont(t) × (1 + E·C/(EC50+C))`, so
when `GENO = 0` the induction term vanishes along with the enzyme. A gene product
that does not exist cannot be induced.

### And this is where the usual explanation collapses

Every photoreaction rate in this model is `photon delivery × concentration within
the irradiated layer`, and never `× amount`. The photoreaction happens in the
volume the light reaches, so it must be proportional to the irradiated **area** and
the local **concentration**. (During development there was a version that wrote
this term against amount, and in that case photon delivery scaled as BSA×W and the
conclusion below was quietly inverted. It was caught by comparing the derivatives
of the two implementations.)

The photoremoval rate per kg therefore falls as `BSA/W ~ W^-0.30`. This is exactly
the term usually offered as the reason CN-I patients lose phototherapy as they
grow. Holding body size fixed, quasi-steady-state TSB was computed with three
nested models (A11b, 12 h/day intensive phototherapy):

| Age | W (kg) | BSA (cm²) | BSA/W | Production (mg/kg/d) | TSB [G] | TSB [G+P] | TSB [G+P+O] |
|------|------|------|------|------|------|------|------|
| 0 | 3.40 | 2209 | 649.8 | 7.68 | 4.80 | 4.80 | 4.87 |
| 0.25 y | 6.00 | 3253 | 542.1 | 6.53 | 5.46 | 4.65 | 4.88 |
| 0.5 y | 7.80 | 3892 | 498.9 | 5.99 | 5.80 | 4.52 | 4.93 |
| 1 y | 9.60 | 4552 | 474.0 | 5.38 | 6.03 | 4.22 | 4.93 |
| 2 y | 12.20 | 5473 | 448.6 | 4.76 | 6.28 | 3.89 | 5.15 |
| 4 y | 16.30 | 6844 | 419.9 | 4.22 | 6.60 | 3.63 | 5.79 |
| 6 y | 20.50 | 8105 | 395.4 | 3.99 | 6.90 | 3.59 | 6.45 |
| 8 y | 25.30 | 9417 | 372.2 | 3.84 | 7.21 | 3.60 | 6.95 |
| 10 y | 31.90 | 11008 | 345.1 | 3.69 | 7.61 | 3.66 | 7.31 |
| 14 y | 50.00 | 14990 | 299.8 | 3.53 | 8.43 | 3.88 | 7.98 |
| 18 y | 66.00 | 17937 | 271.8 | 3.43 | **9.05** | **4.04** | **8.36** |

- **[G] geometry alone**: the limiting value rises 1.89-fold. The
  surface-area-to-mass argument is real.
- **[G+P] plus the fall in production per kg**: the limiting value **flattens**
  (0.84-fold). BSA/W shrinks to 0.42-fold from birth to age 18, but bilirubin
  production per kg also shrinks to 3.80/8.50 = 0.45-fold, so **the two terms very
  nearly exactly cancel.**
- **[G+P+O] plus skin optics**: the loss returns (1.72-fold). Applying that ratio
  to the neonatal limiting value of ~15 mg/dL observed clinically gives ~26 mg/dL
  in adolescence.

So **"the lamps stop working because the surface-area-to-mass ratio shrinks" cannot
on its own account for the clinical course.** What is left is the optical term (the
skin thickens and gains pigment, so the pigment pool within the penetration depth
of blue light shrinks) and the exposure time actually achievable. This distinction
changes behaviour. What must be defended is not surface area but **delivered
irradiance and hours per day**, and the real answer is to replace the enzyme before
the optical term wins.

### Enzyme replacement — 10 % of adult activity is enough

A simulation of a single dose of AAV8-hUGT1A1 on day 60 with the lamps withdrawn on
day 65 (A11c):

| Day | TSB | Transgene (% of adult UGT1A1) | Bf (nM) |
|------|------|------|------|
| 60 (dose) | 3.95 | 0.00 | 2.7 |
| 65 (lamps removed) | 3.57 | 2.93 | 2.4 |
| 70 | 7.79 | 5.00 | 7.0 |
| 80 | 10.36 | 7.50 | 10.3 |
| 100 | 9.31 | 9.38 | 8.9 |
| 140 | 9.31 | 9.96 | 8.9 |

**About 10 % of adult activity converts CN-I into a CN-II phenotype** — exactly
what the first-in-human GNT0003 results showed. The model also says why so little
enzyme is enough: the pathway is not saturated, so activity and removal rate are
still in the linear region.

---

## Other things the model computes

### Exchange transfusion is not "fast phototherapy" (A6, A7)

170 mL/kg is exchanged over three hours. Nowhere in the model is there a sentence
saying "it removes 50 %". The amount removed is the integral `∫ Q_ET·C_p dt` and is
limited by the rate at which the extravascular pool refills the plasma.

| Time | TSB | Bf (nM) | Albumin | Hb | Antibody load |
|------|------|------|------|------|------|
| 30 h (just before exchange) | 9.87 | 9.8 | 2.98 | 13.2 | 0.287 |
| 33 h (exchange ends) | 5.92 | 3.3 | 3.78 | 13.1 | 0.062 |
| 39 h (+6 h) | 6.09 | 3.6 | 3.78 | 13.1 | 0.062 |
| 81 h | 3.63 | 2.1 | 3.75 | 13.8 | 0.058 |

**TSB comes down 40 % while Bf comes down 66 %.** It is because the donor plasma
resets albumin (2.98 → 3.78 g/dL) and washes maternal antibody out from 0.287 to
0.062. Exchange transfusion acts **simultaneously** on all three terms of the free
bilirubin expression — burden, binding capacity, and production rate. Phototherapy
acts on only one.

How the treatment ladder separates in severe Rh disease detected late (A6,
readmission at 48 hours):

| Treatment | TSB 48 h | 60 h | 72 h | Peak Bf | Time above exchange threshold | P(kernicterus) |
|------|------|------|------|------|------|------|
| No treatment | 15.77 | 18.24 | 20.45 | 135.9 | 172 h | 0.996 |
| Early phototherapy at 18 h | 9.11 | 8.67 | 8.27 | 10.3 | 0 | 0.002 |
| Phototherapy at 48 h | 15.78 | 14.88 | 11.75 | 30.8 | 0 | 0.002 |
| + IVIG | 15.78 | 13.56 | 9.85 | 30.8 | 0 | 0.002 |
| + IVIG + stannsoporfin | 15.78 | 12.22 | 7.60 | 30.8 | 0 | 0.002 |
| + IVIG + exchange transfusion | 15.78 | 10.15 | 7.03 | 30.8 | 0 | 0.002 |

Note that intensive phototherapy started at 18 hours alone solves everything. The
treatment ladder acquires meaning **only when detection is late**.

### G6PD deficiency — the quantitative reason a normal-looking haemoglobin is no reassurance (A8)

Destruction of 1 g/dL of haemoglobin is 0.85 g/kg of Hb, that is 29 mg/kg of
bilirubin. Spread through a distribution space of 2.5 dL/kg that becomes **a TSB of
11.6 mg/dL**. A 3 g/dL fall in Hb carries about 35 mg/dL of pigment with it.

| Arm | Peak TSB | Peak Bf | Hb nadir | ETCOc |
|------|------|------|------|------|
| No trigger | 8.54 | 8.3 | 17.00 | 1.38 |
| 12 h trigger, no treatment | 19.73 | 35.6 | 16.43 | 8.87 |
| 12 h trigger + intensive phototherapy | 19.21 | 27.3 | 16.43 | 8.87 |
| 24 h trigger, no treatment | 28.63 | 115.1 | 15.81 | 9.95 |
| 24 h trigger + intensive phototherapy | 21.83 | 37.6 | 15.81 | 9.95 |
| 24 h trigger + phototherapy + exchange transfusion | 20.81 | 33.1 | 14.03 | 9.75 |

Hb falls only from 17.0 to 15.8, yet TSB reaches 28.6. "The Hb is fine" is no
grounds for reassurance.

### Genotype × feeding — prolonged jaundice is an interaction (A9)

| Genotype | Feeding | Peak TSB | t_peak (h) | d7 | d14 | d28 |
|------|------|------|------|------|------|------|
| Wild type | Breast milk | 8.54 | 108 | 7.29 | 2.89 | 1.84 |
| Wild type | Formula | 8.01 | 96 | 6.23 | 2.52 | 1.74 |
| Wild type | Breast milk, suboptimal intake | 9.58 | 126 | 9.00 | 3.64 | 1.89 |
| *28 heterozygous | Breast milk | 10.30 | 144 | 10.10 | 4.89 | 2.29 |
| *28 heterozygous | Breast milk, suboptimal intake | 11.54 | 162 | 11.53 | 6.87 | 2.49 |
| Gilbert (*28/*28) | Breast milk | 13.53 | 210 | 13.16 | 11.37 | 3.94 |
| Gilbert | Formula | 12.62 | 192 | 12.50 | 9.48 | 3.32 |
| Gilbert | Breast milk, suboptimal intake | **14.77** | 228 | 14.06 | **13.37** | 5.27 |
| *6/*6 | Breast milk | 12.76 | 192 | 12.62 | 9.83 | 3.36 |
| *6/*6 | Breast milk, suboptimal intake | 14.04 | 216 | 13.62 | 12.04 | 4.25 |

Suboptimal intake raises the peak without touching UGT1A1. Because the
enterohepatic term `KREAB·BGU·(1−occ)` is a flux competing with removal,
**feeding support has the same units as drug effect.** That is why AAP 2022 calls
this "suboptimal intake hyperbilirubinaemia", and in the model the contribution of
breast milk itself (β-glucuronidase) is late and gentle while the contribution of
intake volume is early and large.

### Treating enterohepatic blockade pharmacologically (A10)

| Intervention | Peak TSB | d7 TSB | Cumulative excretion (mg) |
|------|------|------|------|
| None (breast milk, normal intake) | 8.54 | 7.29 | 122 |
| Oral agar 250 mg/kg/day | 7.64 | 5.60 | 138 |
| UDCA 10 mg/kg q12h | 7.91 | 5.99 | 137 |
| Formula supplementation | 8.01 | 6.23 | 132 |
| Phenobarbital 5 mg/kg/day | 6.73 | 3.25 | 159 |

### Phototherapy dose-response — area beats irradiance (A4)

| Irradiance (µW/cm²/nm) | Exposed BSA | TSB 24 h | 48 h | 72 h | Change over 24 h |
|------|------|------|------|------|------|
| 0 | — | 7.81 | 11.75 | 14.81 | +50.4 % |
| 8 | 0.35 | 7.82 | 11.13 | 12.40 | +42.4 % |
| 15 | 0.35 | 7.82 | 10.49 | 11.00 | +34.2 % |
| 30 | 0.35 | 7.82 | 9.58 | 9.32 | +22.6 % |
| 15 | 0.80 | 7.82 | 8.69 | 7.90 | +11.1 % |
| **30** | **0.80** | 7.83 | 7.28 | 6.07 | **−7.0 %** |
| 50 | 0.80 | 7.83 | 6.39 | 5.14 | −18.3 % |
| 30 | 1.00 | 7.83 | 6.56 | 5.30 | −16.2 % |

Widening the exposed area from 0.35 to 0.80 (+22.6 → −7.0 %) is a great deal
better than raising irradiance from 30 to 50 (−7.0 → −18.3 %). The saturation was
not assumed; it is optical (`I50 = 45 µW/cm²/nm`).

### Same TSB, different disease — the iso-TSB experiment (A12)

TSB is **held fixed** at 18 mg/dL and only the phenotype is changed.

| Phenotype | Bf (nM) | Brain bilirubin (nM-eq) | 24 h injury accumulation rate\* |
|------|------|------|------|
| Term, albumin 3.4, pH 7.40 | 29.8 | 29.8 | 0.0000 |
| Term, albumin 2.8 | 48.3 | 48.3 | 0.1592 |
| Late preterm 35 weeks (binding ×0.85, BBB ×1.6) | 51.6 | 82.5 | 0.5702 |
| Preterm 30 weeks (binding ×0.75, BBB ×2.2) | 104.8 | 230.5 | 2.3464 |
| Term + sepsis/acidosis pH 7.15, albumin 2.8 | 68.9 | 137.9 | 1.2346 |
| Term + ceftriaxone (displacement ×0.70) | 48.9 | 48.9 | 0.1668 |

\*Unconstrained accumulation rate. The ODE itself saturates at INJ = 1.

**At the same 18 mg/dL the 30-week preterm accumulates injury and the term infant
accumulates none at all.** Nowhere in the model is there a sentence saying "preterm
infants are more vulnerable". The binding isotherm and the barrier permeability say
so.

### Threshold-based closed-loop phototherapy (A13)

In the scenario library, phototherapy switches on when that baby's own AAP
threshold is crossed and off when it falls 2 mg/dL below it. **Phototherapy hours
are therefore an output of the model.**

| Scenario | Peak TSB | t_pk | Peak Bf | Phototherapy h | Above exchange threshold | ETCOc | Hb d7 | Injury |
|------|------|------|------|------|------|------|------|------|
| S1 physiological term, breast milk | 8.54 | 108 | 8.3 | 0 | 0 | 1.38 | 17.4 | 0.000 |
| S2 suboptimal intake, delayed stooling | 9.33 | 100 | 9.4 | 0 | 0 | 1.38 | 16.5 | 0.000 |
| S3 ABO alloimmunisation (DAT+) | 18.43 | 140 | 29.7 | 0 | 0 | 2.48 | 16.4 | 0.000 |
| S4 Rh disease + IVIG | 20.50 | 96 | 40.7 | 10 | 0 | 3.37 | 14.0 | 0.000 |
| S5 Rh disease + IVIG + exchange transfusion | 11.98 | 112 | 12.8 | 0 | 0 | 3.37 | 14.7 | 0.000 |
| S6 G6PD 24 h oxidative crisis | 23.31 | 84 | 46.0 | 22 | **14** | 9.95 | 16.4 | 0.001 |
| S7 late preterm 35 weeks + risk factors | 16.04 | 130 | 32.4 | **40** | 0 | 3.22 | 16.4 | 0.042 |
| S8 UGT1A1*6/*6 prolonged jaundice | 12.76 | 194 | 14.8 | 0 | 0 | 1.38 | 17.4 | 0.000 |
| S9 Crigler-Najjar II + phenobarbital | 20.28 | 336 | 35.6 | 0 | 0 | 1.38 | 17.4 | 0.001 |
| S10 ABO + stannsoporfin | 11.98 | 154 | 13.4 | 0 | 0 | 2.46 | 16.4 | 0.000 |

Two things stand out.

- **S7** needs the most phototherapy (40 hours). Its haemolytic burden is moderate,
  but its threshold is the lowest and its free bilirubin at a given TSB is the
  highest.
- **S9 is a useful failure.** The Crigler-Najjar II baby sits **just below** the
  threshold curve at 20 mg/dL. The AAP curves were built on the premise of a
  transient jaundice with a removal rate rising behind it, and a conjugation defect
  does not keep that premise. A threshold-based controller **precisely undertreats
  the baby who is not going to come down on its own.**

---

## The 34 state variables

| # | State | Meaning | Unit |
|---|------|------|------|
| 1 | `W` | Body weight (tracking WHO reference, including the physiological loss nadir) | kg |
| 2 | `HB` | Haemoglobin | g/dL |
| 3 | `RET` | Reticulocytes (marrow output signal) | % |
| 4 | `ABMAT` | Maternal alloantibody load bound to neonatal red cells | – |
| 5 | `HO1` | Haem oxygenase-1 activity | – |
| 6 | `COHB` | Carboxyhaemoglobin → ETCOc | % |
| 7 | `BLEX` | Haemoglobin in extravascular sequestered blood (cephalohaematoma) | g |
| 8 | `BP` | Plasma unconjugated bilirubin | mg |
| 9 | `BEX` | Extravascular (skin, interstitial) unconjugated bilirubin — **the photon target** | mg |
| 10 | `BLIV` | Hepatocyte unconjugated bilirubin (ligandin pool) | mg |
| 11 | `BCON` | Hepatocyte conjugated bilirubin | mg |
| 12 | `BGC` | Intestinal luminal conjugated bilirubin | mg |
| 13 | `BGU` | Intestinal luminal unconjugated bilirubin (reabsorbable) | mg |
| 14 | `BSTL` | Cumulative faecal and renal excretion | mg |
| 15 | `LUMI` | Plasma lumirubin | mg |
| 16 | `EZ` | Plasma (4Z,15E) configurational photoisomer | mg |
| 17 | `ALB` | Serum albumin | g/dL |
| 18 | `UGT` | UGT1A1 activity (relative to adult) | – |
| 19 | `TGX` | AAV8-hUGT1A1 transgene activity | – |
| 20 | `BBR` | Basal ganglia bilirubin | nM-eq |
| 21 | `INJ` | Cumulative neural injury index | – |
| 22 | `ABRD` | Auditory brainstem response deficit (**reversible**) | – |
| 23-24 | `ASNMP`, `CSNMP` | Stannsoporfin IM depot · plasma | mg · µg/mL |
| 25-26 | `APB`, `CPB` | Phenobarbital gut depot · plasma | mg · mg/L |
| 27-28 | `IGGC`, `IGGP` | IVIG central · peripheral | g |
| 29-30 | `AUDCA`, `CUDCA` | UDCA gut depot · plasma | mg · µmol/L |
| 31 | `GBIND` | Intestinal luminal binder (agar/charcoal/zinc) | mg |
| 32 | `AUCX` | AUC above the AAP phototherapy threshold | mg·h/dL |
| 33 | `PTH` | Cumulative phototherapy exposure | h |
| 34 | `ETV` | Cumulative exchange transfusion volume | mL/kg |

---

## Calibration

| Target | Literature value | Model value | Source section |
|------|------|------|------|
| Bilirubin production | 6-10 mg/kg/day (twice the adult value) | 7.68 | A14a |
| Distribution space | larger than plasma, mostly extravascular | 2.50 dL/kg, 80 % extravascular | A14b |
| Physiological peak TSB (breast milk) | 8-10 mg/dL, 96-120 h | 8.54, 108 h | A2 |
| d7 / d14 TSB | 5-8 / 2-5 | 7.29 / 2.89 | A2 |
| Free bilirubin (TSB 15, albumin 3.5) | 15-25 nM | 19.6 | A1a |
| ETCOc normal / haemolytic | 1.3-1.7 / 2-4 ppm | 1.38 / 2.5-3.4 | A3 |
| 24 h fall under intensive phototherapy | 30-40 % | 39 % | A4 |
| Photostationary-state E-isomer | 20-25 % of TSB | 20.6 % | A5 |
| Immediate TSB fall from exchange transfusion | ~50 % | 40 % | A7 |
| CN-II untreated plateau value | 6-25 mg/dL | 23.2 | A11a |
| UGT1A1 maturation (adult by 14 weeks) | adult level | 99.9 % | A14d |
| RK4 integration error (dt 0.20 → 0.01) | — | 8.5403 → 8.5404 | A14e |

---

## Verification

There is no R runtime in this environment. The equations and parameters were
therefore developed and calibrated in an **independent pure python RK4
implementation** (`nhb_reference_check.py`), and every number in this README and in
the header of the R file comes from that run output
(`nhb_reference_output.txt`). If the two files disagree, the one that was actually
executed is the python side.

The transcription into the R model was verified as follows.

1. The C++ of `$GLOBAL`/`$ODE`/`$TABLE` was extracted and compiled standalone with
   `g++ -Wall` — no warnings.
2. The derivatives of the two implementations were compared at a deliberately
   awkward state: 31.5 hours after birth, phototherapy running, a double-volume
   exchange transfusion in progress, the G6PD oxidative trigger active, and IVIG,
   stannsoporfin, phenobarbital, UDCA and an intestinal binder all present in the
   body.
3. **All 34 derivatives and the derived outputs (TSB, Bf, B/A, TcB, ETCOc) agreed to
   within a relative error of 1e-8.** The single exception was 1.3e-9 for brain
   bilirubin, which comes from the free-bilirubin bisection iteration running 40
   times in C++ and 60 times in python.

This establishes that the two implementations are the same model. Whether mrgsolve
handles the R plumbing (`$PARAM @annotated` parsing, event scheduling, `param()`
updates) without trouble is, however, **not verified.**

---

## Known biases and limitations (read these)

1. **A double-volume exchange transfusion removes 40-45 % of the total body
   burden.** The classical teaching figure is ~25 %. The cause is treating the
   extravascular compartment as **a single homogeneous compartment**. Because the
   entire extravascular pool refills the plasma with one fast rate constant, the
   fraction of the burden that photons and exchange transfusion can reach is
   overestimated.
2. **For the same reason the phototherapy limiting value in CN-I is optimistic**
   (4.9-8.4 against a reported 15-25 mg/dL). These two biases are **linked.**
   `PSXKG` (plasma↔interstitium permeability) also determines phototherapy
   efficacy, so fitting one breaks the other. Phototherapy calibration was given
   priority. The conclusion of Theme 4 concerns the **ratio** between columns and is
   therefore insensitive to this offset.
3. **`TAUUGT` is not UGT1A1 protein maturation alone.** It bundles the maturation of
   OATP1B1 and of ligandin into a single time constant, calibrated against the TSB
   trajectory. It is not a measurement of UGT1A1 protein amount.
4. **The AAP 2022 thresholds are an analytic approximation** (a gestational-age
   plateau plus a 36-hour rise constant). They are not the tabulated values and
   cannot be used in care.
5. **The direction of the skin optics terms (`FOPTMIN`, `TAUOPT`) is established but
   the values are illustrative.** Theme 4 rests on the **cancellation** of two
   well-measured terms (BSA/W and production per kg), not on the values of these
   two parameters.
6. **The CO model is linear in production.** It gives ETCOc of 9-10 ppm in a violent
   G6PD crisis, which is higher than the 2-5 ppm of the reported cohorts. That is
   because the reported cohorts were not that extreme, and the linearity itself is
   stoichiometrically correct.
7. **Closed-loop phototherapy is implemented differently in python (hysteresis) and
   in R (logistic smoothing).** This is to keep a discontinuous switch out of the
   differential equation solver, and the phototherapy-hours numbers may differ
   slightly between the two implementations.
8. **No outcomes**: long-term neurodevelopmental outcome, readmission cost,
   maternal-infant bonding, and retinal and DNA oxidative damage are on the map only
   and not in the ODEs.

---

## ⚠️ Disclaimer

This model is a **qualitative and semi-quantitative QSP model for education and
research**. It was built from the published literature and clinical data, but it has
not been independently verified or certified, and **it must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.** In
particular, the AAP thresholds in this model are approximations, and actual care must
follow the original guideline and the judgement of the responsible clinicians.
