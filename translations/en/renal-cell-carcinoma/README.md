# Renal Cell Carcinoma (ccRCC) — QSP Model

[![Nodes](https://img.shields.io/badge/Nodes-125%2B-blue)](../../../renal-cell-carcinoma/rcc_qsp_model.dot)
[![ODE Compartments](https://img.shields.io/badge/ODE%20Compartments-18-green)](../../../renal-cell-carcinoma/rcc_mrgsolve_model.R)
[![Regimens](https://img.shields.io/badge/Regimens-7-orange)](../../../renal-cell-carcinoma/rcc_mrgsolve_model.R)
[![References](https://img.shields.io/badge/References-60-red)](../../../renal-cell-carcinoma/rcc_references.md)

## Overview

Clear cell renal cell carcinoma (ccRCC) accounts for about 75% of kidney cancers, and its core pathogenic mechanism is **VHL tumour suppressor gene mutation** (60–90%) causing pVHL loss → HIF-1α/2α accumulation → VEGF · CA9 · GLUT1 overexpression. This model is a Quantitative Systems Pharmacology (QSP) framework integrating the VHL/HIF/VEGF oxygen-sensing pathway, PI3K/AKT/mTOR, MAPK/RAS, the tumour immune microenvironment (TME), and drug PK/PD.

## Pathogenesis Summary

| Pathway | Key Event |
|------|------------|
| VHL/HIF oxygen sensing | VHL mutation → pVHL loss → failure of HIF-2α ubiquitination → HRE transcriptional activation |
| Angiogenesis | HIF-2α → VEGF-A↑ → VEGFR2 phosphorylation → PI3K/MAPK → endothelial cell proliferation and migration |
| mTOR feedback | PI3K/AKT → mTORC1 activation → increased HIF translation (positive feedback) |
| Immune evasion | VEGF → MDSC recruitment, Treg induction, PD-L1 expression → suppression of CD8+ T cell function |
| Metabolic reprogramming | HIF-1α → GLUT1 · LDHA↑ → Warburg effect |

## Model Files

| File | Content |
|------|------|
| [`rcc_qsp_model.dot`](../../../renal-cell-carcinoma/rcc_qsp_model.dot) | Graphviz mechanistic map (125+ nodes, 11 clusters) |
| [`rcc_qsp_model.svg`](../../../renal-cell-carcinoma/rcc_qsp_model.svg) | SVG vector map |
| [`rcc_qsp_model.png`](../../../renal-cell-carcinoma/rcc_qsp_model.png) | PNG thumbnail (150 dpi) |
| [`rcc_mrgsolve_model.R`](../../../renal-cell-carcinoma/rcc_mrgsolve_model.R) | mrgsolve ODE model (18 compartments, 7 treatment scenarios) |
| [`rcc_shiny_app.R`](../../../renal-cell-carcinoma/rcc_shiny_app.R) | 7-tab Shiny interactive dashboard |
| [`rcc_references.md`](../../../renal-cell-carcinoma/rcc_references.md) | 60 PubMed references (14 sections) |

## Mechanistic Map

[![ccRCC QSP Model](../../../renal-cell-carcinoma/rcc_qsp_model.png)](../../../renal-cell-carcinoma/rcc_qsp_model.svg)

*Click to view the enlarged SVG*

## ODE Compartments (18)

| Group | Compartment | Description |
|------|------|------|
| Sunitinib PK | DEPOT_SUN, CENT_SUN, PERI_SUN, MET_SUN | 2-compartment PK + SU12662 active metabolite |
| Nivolumab TMDD | CENT_NIV, PERI_NIV, PD1_FREE, PD1_BOUND | Target-mediated pharmacokinetics (PD-1 binding) |
| Belzutifan | DEPOT_BEZ, CENT_BEZ | HIF-2α inhibitor PK |
| VHL/HIF/VEGF | pVHL, HIF2A, VEGF, VEGFR2_ACT | Oxygen-sensing pathway |
| mTOR | mTOR_ACT | mTORC1 activity (AU) |
| Tumor (Simeoni) | TUM_W1, TUM_W2, TUM_W3, TUM_VOL | 3-transit TGI model |
| Immune TME | CD8_T, TREG, MDSC | Tumour immune microenvironment |

## Treatment Scenarios

| # | Regimen | Clinical Basis | Key Parameters |
|---|------|----------|--------------|
| 1 | Untreated | Control | — |
| 2 | Sunitinib 50 mg (4/2 schedule) | NEJM 2007 (Motzer) | CL=51.8 L/h, SU12662 active metabolite |
| 3 | Nivolumab + Ipilimumab | CheckMate 214 (Motzer 2018) | TMDD: kon=0.32 nM⁻¹h⁻¹ |
| 4 | Pembrolizumab + Axitinib | KEYNOTE-426 (Rini 2019) | Dual VEGFR+PD-1 blockade |
| 5 | Cabozantinib + Nivolumab | CheckMate 9ER (Choueiri 2021) | MET/VEGFR2 + PD-1 |
| 6 | Cabozantinib 60 mg | METEOR (Choueiri 2015) | IC50=0.006 µM |
| 7 | Everolimus 10 mg | RECORD-1 (Motzer 2008) | mTOR IC50=0.15 nM |
| 8 | Belzutifan 120 mg | LITESPARK-005 (Choueiri 2023) | HIF-2α IC50=0.018 µM |

## Shiny App Tabs

| Tab | Content |
|----|------|
| 1. Patient Profile | IMDC risk, initial tumour volume, simulation duration settings |
| 2. PK | Sunitinib/SU12662/Nivolumab/Belzutifan concentration-time curves |
| 3. VHL/HIF/VEGF Pathway | HIF-2α · VEGF · VEGFR2 · mTOR dynamics, HIF-2α SS vs pVHL curve |
| 4. Tumor Dynamics | Simeoni TGI tumour volume curve, waterfall chart (BOR), TGI statistics |
| 5. Immune TME | CD8 · Treg · MDSC dynamics, PD-1 occupancy, MDSC vs VEGF scatter plot |
| 6. Scenario Comparison | Simultaneous comparison of 8 regimens (tumour volume + endpoint table) |
| 7. Biomarker Dashboard | Biomarker Z-score heatmap (week 12 snapshot) |

## How to Run

### mrgsolve Model

```r
library(mrgsolve)
source("rcc_mrgsolve_model.R")
# 위 스크립트가 자동으로 컴파일 → 시뮬레이션 → 플롯 생성
```

### Shiny App

```r
library(shiny)
shiny::runApp("rcc_shiny_app.R")
```

### Graphviz Rendering

```bash
dot -Tsvg rcc_qsp_model.dot -o rcc_qsp_model.svg
dot -Tpng -Gdpi=150 rcc_qsp_model.dot -o rcc_qsp_model.png
```

## Clinical Calibration Summary

| Clinical Trial | Regimen | mPFS (Actual) | Model Optimisation |
|----------|------|------------|------------|
| CheckMate 214 | Nivo+Ipi | 11.6 mo | PD-1 occ, CD8 kill |
| KEYNOTE-426 | Pembro+Axitinib | 15.1 mo | Axitinib IC50, dual blockade |
| CheckMate 9ER | Cabo+Nivo | 16.6 mo | Cabo MET/VEGFR IC50 |
| CLEAR | Len+Pembro | 23.9 mo | (lenvatinib not separately modelled) |
| METEOR | Cabozantinib | 7.4 mo | Cabo monotherapy dose-response |
| RECORD-1 | Everolimus | 4.0 mo | mTOR IC50 calibration |
| LITESPARK-005 | Belzutifan | ORR 22% | HIF-2α inhibition IC50 |
