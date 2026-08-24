# Bartter Syndrome & Gitelman Syndrome QSP Model
### Bartter & Gitelman Syndrome — Salt-Losing Tubulopathies · BGS

> **One lesion, two locations.** This model does not treat Bartter syndrome and
> Gitelman syndrome as two separate diseases. Its purpose is to impose only a
> **single structural difference — placing the same NaCl-transport defect above
> the macula densa (TAL) or below it (DCT)** — and verify whether every clinical
> finding that distinguishes the two conditions **falls out as a computed
> result.**

---

## 1. Deliverables

| File | Content | Scale |
|------|------|------|
| [`bgs_qsp_model.dot`](bgs_qsp_model.dot) | Mechanistic map (Graphviz source) | **229 nodes · 296 edges · 15 clusters** |
| [`bgs_qsp_model.svg`](bgs_qsp_model.svg) | Vector rendering (zoomable, searchable) | 338 KB |
| [`bgs_qsp_model.png`](bgs_qsp_model.png) | Raster rendering (150 dpi) | 9.8 MB |
| [`bgs_mrgsolve_model.R`](bgs_mrgsolve_model.R) | mrgsolve ODE model | **38 ODEs · 17 scenarios · 10 genotype presets** |
| [`bgs_shiny_app.R`](bgs_shiny_app.R) | Shiny dashboard | **10 tabs** |
| [`bgs_references_en.md`](bgs_references_en.md) | References | **71 PubMed links · 10 sections** |

Reproducing the renders:

```bash
dot -Tsvg bgs_qsp_model.dot -o bgs_qsp_model.svg
dot -Tpng -Gdpi=150 bgs_qsp_model.dot -o bgs_qsp_model.png
```

Running it:

```r
source("bgs_mrgsolve_model.R")
print(bgs_vpop())                       # full genotype sweep -> phenotype table
bgs_summary(bgs_run_all(c("S3","S5","S6"), days = 365))
shiny::runApp("bgs_shiny_app.R")        # dashboard
```

---

## 2. Disease Background

Bartter syndrome and Gitelman syndrome are **autosomal recessive salt-losing
tubulopathies** caused by impaired NaCl reabsorption in the renal tubule. Their
shared presentation is **hypokalaemic, hypochloraemic metabolic alkalosis +
hyperreninaemic hyperaldosteronism + normal blood pressure** — a combination
that is itself the signature of "the kidney is losing salt."

| | Bartter | Gitelman |
|---|---|---|
| Lesion site | Thick ascending limb of Henle's loop (**TAL**) | Distal convoluted tubule (**DCT**) |
| Relative to macula densa | **Upstream / part of the lesion** | **Downstream** |
| Gene | `SLC12A1`(I) · `KCNJ1`(II) · `CLCNKB`(III) · `BSND`(IVa) · `CLCNKA+KB`(IVb) · `MAGED2`(V) | `SLC12A3` |
| Onset | Antenatal–neonatal (I·II·IVa·V) / childhood (III) | Adolescence–adulthood |
| Urinary PGE₂ | **↑↑ 2–10-fold** (hyperprostaglandin E syndrome) | **Normal** |
| Urinary Ca | **Hypercalciuria** → nephrocalcinosis → CKD | **Hypocalciuria** |
| Serum Mg | Normal to mildly reduced (antenatal form usually normal) | **↓↓ Markedly reduced** |
| Urinary concentrating ability | Severely impaired (polyuria, polyhydramnios) | Largely preserved |
| Hearing loss | Present in type IV (barttin is also expressed in the inner ear) | Absent |
| COX inhibitors | **Effective (disease-modifying)** | **Largely ineffective** |

Gitelman syndrome is the most common inherited tubulopathy, with a prevalence
of roughly 1–10/40,000 (carrier frequency ~1%); Bartter syndrome affects
roughly 1/1,000,000.

---

## 3. The Model's Core Idea — Location Creates the Phenotype

The model contains **no genotype-specific switch whatsoever.** Only two
capacity parameters move.

```
FTAL  = NaCl-transport capacity of the thick ascending limb (0-1)
FDCT  = NCC-transport capacity of the distal convoluted tubule (0-1)
```

### 3.1 The macula densa senses NaCl through its own NKCC2

The macula densa sits **at the end of the TAL.** This creates a familiar
paradox. In Bartter syndrome, TAL reabsorption is broken, so the luminal NaCl
that **reaches** the macula densa is actually **higher**. So why do renin and
COX-2 rise?

The answer: the macula densa cell must also take up NaCl through **its own
apical NKCC2** in order to sense it, and **that same NKCC2 carries the same
lesion.** Even though the luminal concentration is high, the cell reads it as
"no salt." This is exactly the logic by which furosemide raises luminal NaCl
while still stimulating renin.

In the model this is a single line:

```cpp
LUMREL  = LDCT / LDCTR;                                  // NaCl delivered (high)
MDSENSE = ATAL * pow(LUMREL, EMD);                       // NaCl sensed (low)
```

Because `ATAL` is multiplied in, **a low FTAL lowers the sensed signal, and
FDCT never appears in `MDSENSE` at all.** The anatomical order is the order of
the equation. On tab 6 of the Shiny app, **the divergence between these two
lines is Bartter syndrome itself.**

### 3.2 What follows from this (results not coded into the model)

| Result | Derivation path |
|---|---|
| **Urinary PGE₂ rises only with a TAL lesion** | `MDSENSE down -> COX-2 induction -> mPGES-1 -> PGE2` |
| **COX inhibition is effective only in Bartter syndrome** | Gitelman never engages that loop at all, so identical target occupancy produces no effect |
| **PGE₂ positive feedback (vicious cycle)** | `PGE2 -> further NKCC2 inhibition -> MDSENSE falls further` |
| **The sign of urinary calcium flips** | TAL: loss of K⁺ recycling collapses the lumen-positive potential (+8 mV), abolishing the paracellular driving force for Ca²⁺ (hypercalciuria) **vs.** DCT lesion: TAL intact + volume contraction increases proximal tubular Ca reabsorption (hypocalciuria) |
| **Serum Mg: Gitelman < type III < antenatal Bartter** | In the antenatal form, the surge in distal Mg delivery is recovered by a hypertrophied, TRPM6-upregulated DCT (`KMGLOAD`); in Gitelman, the DCT itself is atrophic and cannot recover it |
| **Type III overlaps with Gitelman** | ClC-Kb is expressed in both TAL and DCT — the only genotype that lowers FTAL and FDCT simultaneously |
| **K⁺ does not rise until Mg is repleted first** | Hypomagnesaemia releases the intracellular Mg²⁺ block on ROMK (`KMGROMK`) |
| **Alkalosis corrects with KCl, not KHCO₃** | Chloride depletion raises the renal bicarbonate threshold (`KCLTH`) |
| **The thiazide-loading response is flat in Gitelman** | The NCC that would be blocked is already inactive |

---

## 4. mrgsolve Model Structure (38 ODEs)

| Group | Compartments | Compartment list |
|---|---|---|
| Pharmacokinetics | 13 | Indometacin (gut · central) · celecoxib · amiloride · spironolactone -> canrenone · enalapril -> enalaprilat · oral KCl · oral Mg · oral NaCl |
| Fluid & electrolytes | 9 | ECF volume · Na · extracellular K · **intracellular K** · Cl · HCO₃ · extracellular Mg · **bone/intracellular Mg store** · filtered Ca |
| Signalling | 6 | Renin (PRA) · Ang II · aldosterone · COX-2 · PGE₂ · vasopressin |
| Tubular adaptation | 5 | NCC phosphorylation · ENaC expression · TRPV5 · TRPM6 · **DCT mass** |
| Long-term outcomes | 5 | Nephrocalcinosis grade · functional GFR · IGF-1 · renal SDS · gastric mucosal integrity |

### 4.1 Nephron-segment chain (algebraic)

```
Filtered Na  ->  PT (Ang II-dependent)  ->  TAL (FTAL x PGE2 inhibition x CaSR)
             ->  [macula densa sensing]  ->  DCT (FDCT x NCC phosphorylation)
             ->  ASDN (saturable, ENaC/aldosterone, amiloride block)  ->  urinary Na
```

Distal delivery (`LASDN`) is the coupling variable that **drives K⁺ secretion
and H⁺ secretion simultaneously**, and it is the single cause of the
hypokalaemic metabolic alkalosis shared by both syndromes.

### 4.2 Chloride balance closed by urinary electroneutrality

Rather than fitting chloride as a separate parameter, it is closed with the
urinary electroneutrality equation.

```
U_Cl = U_Na + U_K + 0.55*J_H+ - U_HCO3 - U_(unmeasured anions)
```

This lets "why does alkalosis improve when NaCl is given" fall out without a
separate assumption.

### 4.3 Calibration targets (60 kg adult, eGFR 100)

| Phenotype | K | Mg | HCO₃ | uCa/Cr | uPGE₂ | PRA | Urine output |
|---|---|---|---|---|---|---|---|
| Normal | 4.2 | 0.85 | 25 | 0.30-0.40 | 1.0 | 1.0 | 1.5 L/d |
| Gitelman | 2.8-3.1 | **0.48-0.55** | 30-32 | **< 0.07** | 1.0-1.3 | 2-3 | 1.9 L/d |
| Bartter III | 2.7-3.1 | 0.50-0.62 | 31-34 | 0.25-0.45 | **3-4** | 4-5 | 3.0 L/d |
| Bartter I/II | 2.5-3.0 | **0.72-0.85** | 30-34 | **0.7-1.1** | **6-9** | 8-10 | > 6 L/d |

Units: mmol/L · mmol/mmol · x upper limit of normal.
Underlying references: [`bgs_references_en.md`](bgs_references_en.md)
§1 (Bettinelli 1992, Peters 2002, Konrad 2021, Blanchard 2017).

---

## 5. Scenarios (17 + genotype sweep)

| # | Scenario | What it shows |
|---|---|---|
| S1 | Normal control | Baseline physiology |
| S2 | Gitelman, untreated | Natural history |
| S3 | Bartter III, untreated | Natural history |
| S4 | Antenatal Bartter I, untreated (child) | Polyuria · growth delay · nephrocalcinosis |
| S5 | Bartter III + indometacin 2 mg/kg/d | **Effective** |
| S6 | Gitelman + indometacin 2 mg/kg/d | **Same drug, no effect** <- the key contrast |
| S7 | Bartter I + celecoxib | Trade-off between COX-2 selectivity's renal efficacy and GI safety |
| S8 / S9 | Gitelman + amiloride / spironolactone | Comparing two K⁺-sparing pathways |
| S10 | Gitelman + KCl + Mg | **K does not rise without Mg** |
| S11 | Bartter III combination therapy (indometacin + amiloride + KCl + Mg + NaCl) | Maximal treatment |
| S12 | Bartter I + indometacin + rhGH | Growth endpoints |
| S13 | Thiazide-loading test | **Flat** in Gitelman |
| S14 | Furosemide-loading test | **Flat** in Bartter I/II |
| S15 | NSAID use during gastroenteritis | PG-dependent loss of GFR -> **AKI** |
| S16 | Mg once daily vs. split into 4 doses | Absorption saturation, osmotic diarrhoea |
| S17 | Combination therapy + 60% adherence | Pill burden is the true rate-limiting step |
| S18 | Genotype-library sweep | Two-dimensional phenotype space |

The pivotal experiment is **S5 vs. S6.** Both scenarios use the identical
drug, identical dose, and identical PK/PD parameters, and the trajectory of
`INH2` (renal COX-2 inhibition fraction) is effectively the same. The only
difference is whether the circuit the target sits on is actually turning in
that patient. This prediction matches the Blanchard 2015 JASN randomised
crossover trial (a small ΔK from indometacin in Gitelman syndrome, with a high
discontinuation rate).

---

## 6. Mechanistic Map (15 clusters)

[![BGS QSP map](bgs_qsp_model.png)](bgs_qsp_model.svg)

1. Genetic lesions — Bartter I-V · Gitelman · phenocopies (CaSR, KCNJ10, HNF1B, CLDN16/19, TRPM6, anti-NCC autoantibodies)
2. **Thick ascending limb** — NKCC2 · ROMK K⁺ recycling · ClC-Kb/barttin · lumen-positive potential · claudin-16/19
3. **Distal convoluted tubule** — NCC · WNK-SPAK-KLHL3 network · Kir4.1/5.1 · TRPM6 · TRPV5
4. Proximal tubule and the aldosterone-responsive distal nephron (ENaC · ROMK · BK · H⁺-ATPase · pendrin · AQP2)
5. **Macula densa / JGA** — where location becomes phenotype
6. RAAS · volume regulation · vasopressin
7. Systemic electrolytes and acid-base — including the Mg-ROMK link
8. Neuromuscular and cardiac complications (tetany · QTc · torsade de pointes · chondrocalcinosis · type IV hearing loss)
9. Growth failure and endocrine sequelae (GH/IGF-1 resistance)
10. Structural renal progression (nephrocalcinosis -> interstitial fibrosis -> CKD)
11. Antenatal/neonatal module (polyhydramnios · preterm birth · spontaneous resolution of MAGED2)
12. Drug PK/PD
13. Iatrogenic risk (NSAID GI/AKI risk · hyperkalaemia · Mg-induced diarrhoea · **non-adherence**)
14. Diagnosis, endpoints, and differential diagnosis
15. Legend

---

## 7. Clinical Implications (What the Model Says)

1. **Target occupancy is not efficacy.** In Gitelman syndrome, indometacin
   inhibits renal COX-2 to the same degree as in Bartter syndrome, yet
   produces no clinical benefit. QSP's most basic lesson — look at the
   **circuit**, not the target — is reproduced here quantitatively.
2. **Magnesium comes before potassium.** As long as hypomagnesaemia keeps
   ROMK disinhibited, increasing oral KCl simply runs off into the urine (S10
   vs. KCl alone).
3. **It must be the chloride salt.** Potassium bicarbonate or citrate does not
   correct the alkalosis.
4. **The benefit of NSAIDs is something to titrate, not maximise.** The same
   PGE₂ that drives the disease also protects the gastric mucosa and supports
   glomerular filtration (S15).
5. **Pill burden is the true rate-limiting step.** In S17, 60% adherence
   cancels out most of the biochemical benefit of combination therapy — in a
   disease treated with 6-20 tablets a day for life, that is a design problem,
   not a pharmacology problem.
6. **In children the primary endpoint is growth.** Renal-height SDS
   trajectory, not serum K, is the integrated readout of the TAL lesion, and
   the one axis that PGE₂ inhibition can reverse.

---

## 8. Limitations

- Blood pressure is treated only as an output variable; no haemodynamic model
  is included.
- PTH, calcitriol, and bone metabolism are collapsed into strong homeostatic
  terms.
- Changes in nephron number with renal development and growth, and the fetal
  circulation, are not modelled (the antenatal form is approximated with
  paediatric parameters).
- The spontaneous resolution of type V (MAGED2) is represented only as a
  preset; HIF-1alpha oxygen dependence itself is not solved as a differential
  equation.
- Parameters are a **semi-quantitative** calibration to population summary
  statistics from the public literature and have not been fitted or validated
  against individual patient data.

---

## Disclaimer

This model is a QSP model for **educational and research purposes**. It has
not been independently validated or certified and **must not be used for
real clinical decision-making, prescribing, or regulatory submission.**

---

*Part of the QSP disease model library · [back to repository root](../README_en.md)*
