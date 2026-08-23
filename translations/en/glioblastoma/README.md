# Glioblastoma Multiforme (GBM) — QSP Model

[![Disease](https://img.shields.io/badge/Disease-Glioblastoma-red)](../../../glioblastoma/) [![Category](https://img.shields.io/badge/Category-CNS%20Oncology-blueviolet)](../../../glioblastoma/) [![Drugs](https://img.shields.io/badge/Drugs-TMZ%20%C2%B7%20Bevacizumab%20%C2%B7%20Pembrolizumab%20%C2%B7%20TTF-blue)](../../../glioblastoma/) [![Status](https://img.shields.io/badge/Status-Complete-brightgreen)](../../../glioblastoma/)

## Overview

**Glioblastoma Multiforme (GBM)** is a WHO Grade IV malignant glioma, the most common and prognostically poorest primary brain tumour in adults. Median survival is approximately 14–22 months (depending on MGMT methylation status), and the 5-year survival rate is below 5%.

This model provides a comprehensive QSP framework spanning the major treatment strategies, including the **Stupp protocol** (TMZ + radiotherapy → adjuvant TMZ), **bevacizumab**, **pembrolizumab**, and **TTF**.

---

## Mechanistic Map

[![GBM Mechanistic Map](../../../glioblastoma/gbm_qsp_model.png)](../../../glioblastoma/gbm_qsp_model.svg)

*Click the image to view the full-resolution SVG.*

| Component | Content |
|------|------|
| **Total nodes** | 220+ |
| **Subgraph clusters** | 13 |
| **Key pathways** | EGFR→PI3K/AKT/mTOR, RAS-MAPK, p53-RB, TMZ-MGMT-MMR, VEGF-angiogenesis, immune checkpoints, GSC stem-cell maintenance |

### Cluster list

| # | Cluster | Key components |
|---|---------|--------------|
| 1 | Genomic alterations | IDH1/2, MGMT, EGFR amp, PTEN loss, TP53, TERT, CDKN2A deletion |
| 2 | RTK/RAS-MAPK | EGFR, PDGFR, MET, GRB2-SOS, RAS-RAF-MEK-ERK |
| 3 | PI3K/AKT/mTOR | PI3K, PIP3, PTEN, AKT, TSC1/2, mTORC1/2, S6K1, 4E-BP1 |
| 4 | Cell cycle & apoptosis | CDK4/6-CycD, RB-E2F, p16/p21, BCL-2 family, caspase cascade |
| 5 | DNA damage & repair | O6-MeG, MGMT, MMR, DSB, ATM/ATR-CHK1/2, NHEJ/HR |
| 6 | Tumour microenvironment | TAM-M1/M2, MDSCs, CAFs, TGF-β, IL-6/STAT3, IDO1, ARG1 |
| 7 | Angiogenesis | HIF-1α, VEGF-A/B/C, VEGFR1/2, ANG1/2, FGF2, neovascularisation |
| 8 | Immune response & checkpoints | CD8+ CTL, Treg, PD-1/PD-L1, CTLA-4, TIM-3, LAG-3, TIGIT |
| 9 | Glioblastoma stem cells | CD133+, SOX2, Notch, Wnt/β-cat, SHH/GLI, BMI1, EZH2 |
| 10 | Drug PK | TMZ 2-cmt+BBB, bevacizumab 2-cmt+VEGF binding, pembrolizumab |
| 11 | Blood-brain barrier | BBB tight junction, P-gp, BCRP, Kp,brain, blood-tumour barrier (BTB) |
| 12 | Radiotherapy | LQ model (α/β=10 Gy), OER, TMZ-RT synergy, SRS/SBRT |
| 13 | Clinical endpoints | Tumour volume, RANO criteria, PFS/OS, KPS, MGMT/IDH prognostic factors |

---

## mrgsolve ODE Model

### Compartment design (18 compartments)

| Category | Compartment | Description |
|----------|------|------|
| **TMZ PK** | Gut, Cp_tmz, Cp2_tmz, Cbrain, O6MeG | GI absorption → 2-compartment PK → BBB penetration → O6-MeG lesion |
| **BEV PK** | BEV_Cp, BEV_Cp2, VEGF_free, BEV_VEGF | 2-compartment PK + VEGF binding/dissociation kinetics |
| **Anti-PD1** | APD1_Cp, PD1_occ | 1-compartment PK + PD-1 receptor occupancy |
| **Tumour cells** | Ts, Tr, GSC | Sensitive/resistant/stem-cell subpopulations |
| **Immune** | CD8_eff, Treg_c, TAM_M2 | CD8+ T cells, regulatory T cells, M2-TAM |
| **Angiogenesis** | NV | VEGF-driven neovascularisation index |

### Key equations

**Tumour growth (Gompertz model)**
```
dTs/dt = Ts × [−kg × ln(N/K) − kill_TMZ − kill_RT − kill_CD8 − k_resist]
```

**O6-MeG lesion (TMZ mechanism)**
```
dO6MeG/dt = k_O6 × C_brain − kMGMT × O6MeG − k_O6deg × O6MeG
```
*MGMT methylated: kMGMT = 0.05 h⁻¹ (low repair)*
*MGMT unmethylated: kMGMT = 0.40 h⁻¹ (high repair)*

**RT effect (Linear-Quadratic model)**
```
SF = exp(−αD − βD²)  [α=0.30 Gy⁻¹, β=0.030 Gy⁻², α/β=10 Gy]
```

### Parameter calibration basis

| Parameter | Value | Source |
|----------|---|------|
| TMZ CL | 11.4 L/h | Ostermann 2004 Clin Cancer Res |
| TMZ V₁ | 22.5 L | Baker 2003 Clin Cancer Res |
| TMZ Kp,brain | 0.28 | Ostermann 2004 |
| BEV CL | 0.207 L/day | Lu 2008 Cancer Chemother Pharmacol |
| GBM α/β | 10 Gy | Joiner & van der Kogel 2009 |
| kg (Gompertz) | 0.003 day⁻¹ | Calibrated to Stupp 2005 |

### Treatment Scenarios (7)

| # | Scenario | Supporting trial | Key mechanism |
|---|---------|-------------|---------|
| S1 | Untreated control | — | Natural Gompertz course |
| S2 | Stupp (MGMT methylated) | Stupp 2005 NEJM, Hegi 2005 NEJM | TMZ+RT → adjuvant TMZ; O6-MeG accumulation |
| S3 | Stupp (MGMT unmethylated) | Hegi 2005 NEJM | High MGMT expression → rapid O6-MeG repair → resistance |
| S4 | Stupp + bevacizumab | AVAGLIO 2014, RTOG0825 | VEGF neutralisation → angiogenesis suppression |
| S5 | Stupp + TTF | EF-14 2017 JAMA | Electric field → disrupts mitotic cells |
| S6 | Pembrolizumab + TMZ | Keynote-028, Reardon 2020 | PD-1 blockade → CD8 reactivation |
| S7 | Bevacizumab alone (salvage) | Friedman 2009 JCO | VEGF neutralisation (recurrent GBM) |

---

## Shiny Dashboard

### Tab layout (7 tabs)

| Tab | Content |
|---|------|
| **1. Patient Profile** | Age, KPS, MGMT/IDH status, extent of resection, treatment selection |
| **2. Drug PK** | TMZ plasma/brain PK, bevacizumab PK, pembrolizumab PK |
| **3. DNA Damage & MGMT** | O6-MeG lesion dynamics, MGMT methylation comparison |
| **4. Tumor Dynamics** | Ts/Tr/GSC cell populations, tumour volume trends |
| **5. Scenario Comparison** | Simultaneous comparison of 7 treatment strategies |
| **6. TME & Biomarkers** | CD8/Treg/TAM immune microenvironment, VEGF, PD-1 occupancy |
| **7. Clinical Endpoints** | Tumour diameter (RANO), treatment kill-rate components |

### How to Run

```r
# Install R packages
install.packages(c("shiny", "mrgsolve", "dplyr", "ggplot2",
                   "tidyr", "shinydashboard", "DT", "gridExtra"))

# Run the app
shiny::runApp("gbm_shiny_app.R")
```

---

## Key Clinical Prognostic Factors

| Factor | Favourable prognosis | Unfavourable prognosis | Median OS difference |
|------|----------|---------|------------|
| MGMT methylation | Methylated (+) | Unmethylated (−) | 21.7 vs 12.6 months (Hegi 2005) |
| IDH mutation | IDH-mutant | IDH-WT | ~31 vs ~15 months |
| Extent of resection | GTR | Biopsy | ~3–4-month difference |
| KPS | ≥70 | <70 | Independent prognostic factor |
| Age | <50 years | >70 years | Age-independent prognostic factor |

---

## File List

| File | Description |
|------|------|
| [`gbm_qsp_model.dot`](../../../glioblastoma/gbm_qsp_model.dot) | Graphviz mechanistic map source (220+ nodes, 13 clusters) |
| [`gbm_qsp_model.svg`](../../../glioblastoma/gbm_qsp_model.svg) | Vector-format map (zoomable in browser) |
| [`gbm_qsp_model.png`](../../../glioblastoma/gbm_qsp_model.png) | Raster-format map (150 dpi) |
| [`gbm_mrgsolve_model.R`](../../../glioblastoma/gbm_mrgsolve_model.R) | mrgsolve ODE QSP model (18 compartments, 7 scenarios) |
| [`gbm_shiny_app.R`](../../../glioblastoma/gbm_shiny_app.R) | Shiny interactive dashboard (7 tabs) |
| [`gbm_references.md`](gbm_references.md) | 60+ references (13 sections) |

---

## Pathophysiology Summary

```
EGFR amp / PTEN loss / IDH1 WT
        ↓
PI3K/AKT/mTOR hyperactivation → cell proliferation
RAS/MAPK hyperactivation → invasion
        ↓
Tumour growth + angiogenesis (VEGF-A↑ → HIF-1α)
        ↓
Immunosuppression (M2-TAM + Treg + PD-L1↑)
        ↓
TMZ resistance (MGMT unmethylated) + GSC recurrence
```

**Treatment strategy**:
- **First-line**: TMZ+RT (radiotherapy) (Stupp protocol) → adjuvant TMZ × 6 cycles
- **MGMT methylated**: excellent TMZ response (HR 0.45, Hegi 2005)
- **At recurrence**: bevacizumab ± irinotecan, pembrolizumab, lomustine
- **Newer approach**: TTF (Optune®) → median OS of 20.9 months in the EF-14 study

---

*Model generated: Claude Code Routine (CCR), 2026-06-26*
*Reference guidelines: EANO 2021, NCCN CNS Tumors v2.2024*
