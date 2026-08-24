# Marfan Syndrome — QSP Model

**Directory:** `marfan-syndrome/`
**Abbreviation:** MFS
**Added:** 2026-06-25
**Category:** Connective tissue disease / Cardiovascular genetic disease

---

## Overview

Marfan syndrome (MFS) is an autosomal dominant connective tissue disease
caused by germline mutations in the **FBN1 gene (Chr 15q21.1)**.
Global prevalence is 1/5,000-1/10,000, occurring equally in men and women.

The defective fibrillin-1 protein disrupts the microfibril network in the
extracellular matrix (ECM), leading to:

- **Impaired TGF-β sequestration** -> increased free TGF-β1/2 -> SMAD2/3 +
  ERK/MAPK overactivation
- **ECM degradation** (MMP-2/9 upregulation) -> aortic medial degeneration
- **Aortic root dilation** -> aortic regurgitation, dissection risk
- **Skeletal features:** tall stature, arachnodactyly, scoliosis, pectus
  excavatum/carinatum
- **Ocular features:** lens dislocation (60-70%), myopia, retinal detachment

---

## Mechanistic Map

[![MFS QSP Mechanistic Map](mfs_qsp_model.png)](mfs_qsp_model.svg)

> Click to view the full-resolution SVG.

### 14 Clusters

| Cluster | Content |
|----------|------|
| ① Genetic & Molecular Foundation | FBN1 mutation, fibrillin-1 protein, LTBP, microfibrils |
| ② TGF-β Signalling | SMAD2/3, SMAD4, SMAD7, canonical/non-canonical pathways |
| ③ MAPK/ERK Pathway | RAS-RAF-MEK-ERK, p38, JNK, PI3K/AKT |
| ④ ECM Remodelling | MMP-2/9/13, TIMP-1/2, elastin degradation, collagen fragmentation |
| ⑤ Vascular SMC Biology | VSMC phenotype switching, apoptosis, NOX4, NF-κB |
| ⑥ Aortic Pathology | Aortic root, sinus of Valsalva, STJ, ascending aorta, Laplace stress |
| ⑦ Cardiac Manifestations | Aortic regurgitation, mitral valve prolapse, LV dilation, dissection |
| ⑧ Haemodynamic Parameters | HR, SBP, dP/dt_max, PWV, pulse wave velocity |
| ⑨ Skeletal System | Tall stature, arachnodactyly, pectus excavatum/carinatum, scoliosis |
| ⑩ Ocular System | Lens dislocation, myopia, retinal detachment, glaucoma |
| ⑪ Other Systemic Features | Dural ectasia, pneumothorax, hernia, sleep apnoea |
| ⑫ Drug PK | Atenolol two-compartment, losartan/EXP-3174 PK |
| ⑬ Drug PD | β1-blockade, AT1R blockade, TGF-β inhibition mechanisms |
| ⑭ Clinical Endpoints & Biomarkers | Z-score, annual growth rate, AR grade, plasma TGF-β |

**Total nodes:** 130+ nodes | **Clusters:** 14

---

## mrgsolve ODE Model (Compartmental Model)

**File:** `mfs_mrgsolve_model.R`

### 20 Compartments

| Compartment | Description |
|------|------|
| DEPOT_ATN | Atenolol gut absorption compartment |
| C1_ATN | Atenolol central compartment |
| C2_ATN | Atenolol peripheral compartment |
| DEPOT_LOS | Losartan gut absorption compartment |
| C1_LOS | Losartan central compartment |
| C_EXP3174 | EXP-3174 active metabolite |
| TGFb | Plasma free TGF-β1 [ng/mL] |
| pSMAD | Phosphorylated SMAD2/3 (fold) |
| pERK | Phosphorylated ERK1/2 (fold) |
| MMP | Circulating MMP activity [U/mL] |
| Ao_Diam | Aortic root diameter [mm] |
| AR_Grade | Aortic regurgitation grade (0-4) |
| HR | Heart rate [bpm] |
| SBP | Systolic blood pressure [mmHg] |
| dPdt | dP/dt_max [mmHg/s] |
| NT_proBNP | NT-proBNP [pg/mL] |
| LVEDD | Left ventricular end-diastolic diameter [mm] |
| TGFb_plasma_obs | Observed plasma TGF-β1 |
| Systemic_score | Ghent systemic score (0-20) |

### Six Treatment Scenarios

| Scenario | Treatment | Clinical basis |
|----------|------|----------|
| 1. Untreated | — | Natural history (Salim 1994; Rossig 2019) |
| 2. Atenolol 50 mg QD | β1-blockade | PHN RCT — Lacro et al. NEJM 2014 |
| 3. Atenolol 100 mg QD | β1-blockade (high dose) | Shores et al. NEJM 1994 |
| 4. Losartan 50 mg QD | ARB | PHN RCT — Lacro et al. NEJM 2014 |
| 5. Losartan 100 mg QD | ARB (high dose) | COMPARE — Radonic et al. EHJ 2010 |
| 6. Atenolol + losartan | Combination therapy | AIMS — Forteza et al. JACC 2016 |

### Key Clinical Trial Calibration Data

| Trial | Drug | Duration | Key result |
|----------|------|------|----------|
| Shores 1994 NEJM | Propranolol | 10 years | Atenolol arm vs. untreated: 50% reduction in aortic growth rate |
| Lacro 2014 NEJM (PHN) | Atenolol vs. losartan | 3 years | Similar between arms (aortic root Z-score down ~0.12) |
| Radonic 2010 EHJ (COMPARE) | Losartan 100 mg | 2 years | No significant difference in aortic growth rate vs. atenolol |
| Forteza 2016 JACC (AIMS) | Irbesartan vs. atenolol | 3 years | Equivalent aortic growth rate (0.48 vs. 0.52 mm/yr) |
| Brooke 2008 NEJM | Losartan (paediatric) | ~2.3 years | Markedly reduced growth rate in the losartan arm (historical control) |
| Habashi 2006 Science | Losartan (Fbn1+/- mice) | — | Complete suppression of aortic dilation and TGF-β signalling |

---

## Shiny Dashboard (Interactive Dashboard)

**File:** `mfs_shiny_app.R`

### Seven Tabs

| Tab | Content |
|----|------|
| ① Patient Profile | Ghent criteria table, pathophysiology overview, mechanistic map preview |
| ② Drug PK | Atenolol/losartan/EXP-3174 concentration-time curves, PK summary |
| ③ TGF-β / Mol. PD | TGF-β1, p-SMAD2/3, p-ERK1/2, MMP kinetics |
| ④ Cardiovascular Endpoints | Aortic root diameter, Z-score, AR grade, HR/dP/dt, LVEDD, SBP |
| ⑤ Scenario Comparison | Side-by-side comparison of the 6 treatment scenarios, 5-year summary table |
| ⑥ Biomarkers | NT-proBNP, TGF-β, Ghent score, annual growth rate, threshold reference table |
| ⑦ Surgical Decision Support | Surgical threshold visualisation, ESC/AHA guidelines, time to threshold by treatment |

---

## References

**File:** `mfs_references_en.md` — 50 total

| Section | Citations |
|------|---------|
| 1. Genetics & Molecular Pathogenesis | 8 |
| 2. TGF-β Signalling | 6 |
| 3. Aortic Pathology & Natural History | 6 |
| 4. Clinical Trials — β-Blockers | 4 |
| 5. Clinical Trials — ARBs | 6 |
| 6. Pharmacokinetics | 5 |
| 7. MMP / ECM | 4 |
| 8. Ophthalmology | 2 |
| 9. Skeletal & Systemic | 3 |
| 10. Surgery | 3 |
| 11. QSP / Modelling | 3 |

---

## Usage

```bash
# 1) Render the mechanistic map (requires Graphviz)
dot -Tsvg mfs_qsp_model.dot -o mfs_qsp_model.svg
dot -Tpng -Gdpi=150 mfs_qsp_model.dot -o mfs_qsp_model.png
```

```r
# 2) Run the mrgsolve ODE model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
source("mfs_mrgsolve_model.R")
# -> runs 6 scenarios x 5-year simulation, prints a summary table

# 3) Run the Shiny dashboard
install.packages(c("shiny", "shinydashboard", "DT"))
shiny::runApp("mfs_shiny_app.R")
```

---

## File Structure

```
marfan-syndrome/
├── README.md                  ← this file
├── mfs_qsp_model.dot          ← Graphviz mechanistic map source
├── mfs_qsp_model.svg          ← SVG rendering
├── mfs_qsp_model.png          ← PNG rendering (150 dpi)
├── mfs_mrgsolve_model.R       ← mrgsolve ODE QSP model
├── mfs_shiny_app.R            ← Shiny interactive dashboard
└── mfs_references_en.md          ← References (50 total)
```
