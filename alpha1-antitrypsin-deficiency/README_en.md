# Alpha-1 Antitrypsin Deficiency (AATD) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-AATD%20%2F%20Alpha--1%20Antitrypsin%20Deficiency-red)]()
[![Category](https://img.shields.io/badge/Category-Rare%20Genetic%20%2F%20Pulmonary%20%2F%20Hepatic-orange)]()
[![Nodes](https://img.shields.io/badge/Mechanistic%20Map-130%2B%20nodes%20%7C%2010%20clusters-blue)]()
[![ODE](https://img.shields.io/badge/ODE%20Model-20%20compartments%20%7C%206%20scenarios-green)]()
[![Shiny](https://img.shields.io/badge/Shiny%20App-6%20tabs-purple)]()
[![Refs](https://img.shields.io/badge/References-54%20PubMed-brightgreen)]()

---

## Disease Overview

**Alpha-1 antitrypsin deficiency (AATD)** is a hereditary metabolic
disease caused by mutations in the *SERPINA1* gene, with a distinctive
dual pathophysiology that affects the liver and lungs simultaneously. The
most severe phenotype, the **ZZ genotype (PIZZZ)**, occurs in about
1 in 3,000–5,000 people worldwide.

### The Core Dual Pathophysiology

| Pathway | Mechanism | Result |
|------|------|------|
| **Liver (toxic gain of function)** | Glu342Lys mutation → accumulation of loop-sheet polymers of Z-AAT protein in the endoplasmic reticulum → ER stress → hepatocyte injury/fibrosis | Cirrhosis, hepatocellular carcinoma (HCC risk ↑10-40 fold) |
| **Lung (loss of function)** | Serum AAT deficiency (below the 11 µM protective threshold) → unchecked activity of alveolar neutrophil elastase (NE) → elastin degradation → panlobular emphysema | Rapid FEV1 decline (2-3x normal), severe COPD |

---

## Model Components

| Deliverable | File | Specification |
|--------|------|------|
| Mechanistic map | [`aatd_qsp_model.dot`](aatd_qsp_model.dot) / [`.svg`](aatd_qsp_model.svg) / [`.png`](aatd_qsp_model.png) | **130+ nodes, 10 clusters** |
| mrgsolve ODE model | [`aatd_mrgsolve_model.R`](aatd_mrgsolve_model.R) | **20 compartments, 6 treatment scenarios** |
| Shiny dashboard | [`aatd_shiny_app.R`](aatd_shiny_app.R) | **6-tab interactive app** |
| References | [`aatd_references.md`](aatd_references.md) | **54 PubMed citations, 13 sections** |

---

## Mechanistic Map

[![AATD QSP Mechanistic Map](aatd_qsp_model.png)](aatd_qsp_model.svg)

### 10 Subgraph Clusters

| # | Cluster | Key elements |
|---|---------|----------|
| 1 | Genetics & genotype | SERPINA1, M/Z/S/Null alleles, MM/MZ/SZ/ZZ/SS genotypes |
| 2 | ER protein quality control | Endoplasmic reticulum folding (calnexin, BiP/GRP78), Z-AAT polymers, UPR (IRE1α/PERK/ATF6), ERAD, autophagy |
| 3 | Hepatic pathology | Hepatocyte injury, ER stress, hepatic stellate cell (HSC) activation, TGF-β1, collagen deposition, Metavir F0-F4 fibrosis, cirrhosis, HCC |
| 4 | AAT biology & distribution | M-AAT/Z-AAT secretion, serum/ELF distribution, inhibition of NE/PR3/cathepsin G, anti-inflammatory function |
| 5 | Protease-antiprotease balance | Pulmonary PMN recruitment, NE release/activation, MPO/ROS, MMP-12, MMP-9, SLPI, elafin, TIMP-3, AAT oxidation |
| 6 | Pulmonary pathology | Degradation of elastin/fibronectin/collagen/laminin, ECM destruction, panlobular emphysema, FEV1, FVC, DLCO, RV/TLC, bronchiectasis |
| 7 | Pulmonary inflammatory cascade | IL-8, LTB4, TNF-α, IL-1β, IL-17, IFN-γ, IL-10, NF-κB, AP-1, complement C5a, CD8/Th17/NK cells |
| 8 | Drug PK/PD | Augmentation therapy (Prolastin-C/Zemaira/Aralast/Glassia), siRNA (fazirsiran/belcesiran), correctors (VX-864/GSK3117391), NE inhibitors (alvelestat/lonodelestat/brensocatib), gene therapy (rAAV/CRISPR), LABA/LAMA/ICS |
| 9 | Treatment PD effects | Serum AAT↑, ELF AAT↑, pulmonary NE↓, Z-AAT polymer↓, slowed FEV1 decline, preserved CT lung density, reduced exacerbations |
| 10 | Clinical endpoints | SGRQ, CAT, mMRC, 6MWD, serum AAT, ELF AAT, CT emphysema index, hepatic biomarkers, mortality |

---

## mrgsolve ODE Model

### Compartment Structure (20)

| Group | Compartments | Description |
|------|------|------|
| **Hepatic Z-AAT** | `ZAAT_ER`, `ZAAT_Poly`, `HSC_act`, `Liver_coll`, `Liver_fib` | Z-AAT in the ER, polymers, hepatic stellate cells, collagen, fibrosis |
| **Serum AAT** | `AAT_C1`, `AAT_C2` | Serum AAT 2-compartment PK (V1=3.76 L/kg, t½=4.5 days) |
| **Pulmonary pathway** | `PMN_lung`, `IL8_lung`, `NE_free`, `MMP12_lung`, `Elastin`, `FEV1_pct` | PMN, IL-8, free NE, MMP-12, elastin, FEV1 |
| **Drug PK** | `AUG_C1`, `AUG_C2`, `NEi_A`, `NEi_C`, `siRNA_Eff`, `Gene_Eff` | Augmentation therapy 2-compartment, NE inhibitor, siRNA effect compartment, gene therapy |

### 6 Treatment Scenarios

| Scenario | Treatment | Supporting trial |
|---------|------|-------------|
| **S1** No treatment | — | Natural history (Tanash 2010, Janus 1985) |
| **S2** AAT augmentation therapy | Prolastin-C 60 mg/kg IV weekly | RAPID trial (Chapman 2015, *Lancet*) |
| **S3** Fazirsiran (siRNA) | GalNAc-siRNA 200 mg SQ q12wk | ARO-AAT phase 2 (Strnad 2022, *NEJM*) |
| **S4** NE inhibitor | Alvelestat 60 mg BID PO | Phase 2 (McElvaney 2020, *AJRCCM*) |
| **S5** Gene therapy | rAAV-SERPINA1 single dose | Flotte 2004, Mueller 2008 |
| **S6** Augmentation + NE inhibitor | Combination therapy | Simulation-based prediction |

### Key Parameter Calibration Basis

| Parameter | Value | Source |
|---------|-----|------|
| ZZ serum AAT | 6-8 mg/dL (~2-7 µM) | Crystal 1990 |
| Protective threshold | >11 µM (57 mg/dL) | RAPID trial 2015 |
| ELF/serum ratio | 0.10 | Hubbard 1991 |
| AAT half-life | 4.5 days | Prolastin-C PK label |
| Annual FEV1 decline (ZZ) | 50-200 mL/yr | Dirksen 1999/2009 |
| siRNA mRNA knockdown | ~88% | Strnad 2022 |
| Alvelestat NE inhibition | ~75% | McElvaney 2020 |

---

## Shiny App Tab Structure

| Tab | Contents |
|----|------|
| 1. Patient profile | Genotype, smoking status, FEV1, CT emphysema index, risk summary, drug targets |
| 2. Drug PK / AAT levels | Serum AAT over time, ELF AAT, ELF NE inhibition (%), PK parameter table |
| 3. Pulmonary PD | Free NE activity, elastin content (%), MMP-12, pulmonary PMN burden |
| 4. Clinical endpoints | FEV1 (% predicted), CT emphysema index, SGRQ quality of life, annual exacerbation risk |
| 5. Scenario comparison | Graphs comparing FEV1/serum AAT/Z-polymer across the 6 treatment arms, 5-year outcome summary table |
| 6. Biomarkers | Hepatic Z-AAT polymer, Metavir fibrosis, hepatic stellate cell activation, desmosine (an elastin-degradation marker) |

---

## Quick Start

```r
# 1. Render the mechanistic map
# dot -Tsvg aatd_qsp_model.dot -o aatd_qsp_model.svg
# dot -Tpng -Gdpi=150 aatd_qsp_model.dot -o aatd_qsp_model.png

# 2. mrgsolve simulation
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))
source("aatd_mrgsolve_model.R")   # automatically runs and visualises the 6 scenarios

# 3. Shiny dashboard
install.packages(c("shiny", "shinydashboard", "plotly", "DT", "shinycssloaders"))
shiny::runApp("aatd_shiny_app.R")
```

---

## Key Clinical Features Summary

- **Prevalence**: ZZ genotype about 1/3,000-5,000 (European ancestry > Asian ancestry)
- **Diagnostic delay**: average 5-7 years (typical of a rare disease)
- **Pulmonary phenotype**: 80% of ZZ adults, lower-lobe-predominant panlobular emphysema
- **Hepatic phenotype**: neonatal cholestasis, adult cirrhosis/HCC
- **Effect of smoking**: accelerates FEV1 decline 2-3 fold
- **ERS/ATS guideline**: AAT augmentation therapy recommended when FEV1 <80% and ZZ is confirmed

---

## References

54 PubMed citations — see [aatd_references.md](aatd_references.md)
Key clinical trials: RAPID trial (2015, *Lancet*), ARO-AAT phase 2 (2022, *NEJM*), alvelestat phase 2 (2020, *AJRCCM*)

---

*Model generated by: Claude Code Routine (CCR) | 2026-06-24*
