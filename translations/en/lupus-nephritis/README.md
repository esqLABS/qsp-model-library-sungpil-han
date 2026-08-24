# Lupus Nephritis — QSP Model

> **ISN/RPS Class III/IV/V · Glomerulonephritis · SLE-driven kidney disease**

[![Mechanistic Map](../../../lupus-nephritis/ln_qsp_model.png)](../../../lupus-nephritis/ln_qsp_model.svg)

---

## Overview

Lupus nephritis (LN) is one of the most serious organ complications of
**systemic lupus erythematosus (SLE)**, occurring in roughly 40-60% of SLE
patients. It presents a complex pathophysiology in which immune complexes
(IgG/dsDNA/C1q) deposit in the glomerulus, leading to complement
activation, podocyte injury, proteinuria, and declining renal function.

Lupus Nephritis is a major cause of morbidity and mortality in SLE, requiring careful mechanistic modeling of the immune-renal axis to design effective treatment strategies.

---

## Key Pathophysiology

| Pathway | Key elements |
|------|----------|
| **Type I IFN signalling** | pDC -> IFN-α -> BAFF up -> B-cell hyperactivation |
| **B-cell hyperactivation** | Tfh help -> GC -> long-lived plasma cells -> anti-dsDNA IgG |
| **Immune complexes** | anti-dsDNA + C1q -> IC -> glomerular mesangial deposition |
| **Complement activation** | C1q -> C4 -> C3 -> C5a (neutrophil chemoattraction) + C5b-9 (podocyte attack) |
| **Podocyte injury** | Slit diaphragm (nephrin) breakdown -> proteinuria (UPCR >= 0.5) |
| **Fibrosis** | TGF-β1 -> myofibroblasts -> interstitial fibrosis -> declining eGFR |

---

## Model Files

| File | Description |
|------|------|
| [`ln_qsp_model.dot`](../../../lupus-nephritis/ln_qsp_model.dot) | Graphviz mechanistic map source (15 subgraph clusters) |
| [`ln_qsp_model.svg`](../../../lupus-nephritis/ln_qsp_model.svg) | Vector graphic map |
| [`ln_qsp_model.png`](../../../lupus-nephritis/ln_qsp_model.png) | Raster map (150 dpi) |
| [`ln_mrgsolve_model.R`](../../../lupus-nephritis/ln_mrgsolve_model.R) | mrgsolve ODE model + 5 treatment scenarios |
| [`ln_shiny_app.R`](../../../lupus-nephritis/ln_shiny_app.R) | Shiny dashboard (6 tabs) |
| [`ln_references.md`](../../../lupus-nephritis/ln_references.md) | 45 PubMed references |

---

## Model Specifications

### Mechanistic Map
- **100+ nodes**, 15 subgraph clusters
- Coverage: innate immunity (NETosis, pDC, TLR, cGAS-STING) · adaptive
  immunity (B/T cells) · complement · glomerulus · tubulointerstitium ·
  drug PK/PD (6 agents)

### mrgsolve ODE Model (20 compartments)
```
PK compartments (9):
  MPA_gut, MPA_plasma, MPAG_gut        <- MMF/MPA 2-cmt + EHC
  HCQ_blood, HCQ_tissue                <- HCQ 2-cmt
  VCS_plasma                           <- Voclosporin
  BEL_central, BEL_periph              <- Belimumab 2-cmt IV
  ANI_central, ANI_periph              <- Anifrolumab 2-cmt IV

Immune compartments (5):
  B_naive, B_GC, Plasma_cell, Tfh, Treg

Antibody/complement (3):
  Anti_dsDNA, C3, C4

Renal compartments (3):
  Podocyte_inj, Proteinuria (UPCR), eGFR
```

### Treatment Scenarios (5)
| # | Scenario | Basis |
|---|---------|------|
| 1 | MMF 3g/day + HCQ (SoC) | ACR guidelines 2012 |
| 2 | MMF + HCQ + voclosporin 23.7mg BID | AURORA 1 (Rovin 2021) |
| 3 | MMF + HCQ + belimumab 10mg/kg q4w | BLISS-LN (Furie 2020) |
| 4 | MMF + HCQ + anifrolumab 300mg q4w (high IFN signature) | TULIP-LN (Jayne 2022) |
| 5 | CYC induction (Euro-Lupus) -> MMF maintenance | Houssiau 2002 |

### Shiny Dashboard (6 Tabs)
| Tab | Content |
|----|------|
| ① Patient Profile | LN class, baseline UPCR/eGFR, IFN signature, BW |
| ② Drug PK | Concentration-time curves for MPA, HCQ, voclosporin, belimumab, anifrolumab |
| ③ Immune Biomarkers | B/T cell compartments, anti-dsDNA, C3/C4, BAFF/IC load |
| ④ Renal Function | Time course of eGFR, UPCR, podocyte injury index |
| ⑤ Clinical Endpoints | Achievement of CRR/PRR, SLEDAI renal score, waterfall plot |
| ⑥ Scenario Comparison | Direct comparison of two regimens |

---

## Clinical Endpoints

| Endpoint | Definition |
|-----------|------|
| **Complete Renal Response (CRR)** | UPCR < 0.5 g/g + eGFR >= 60 mL/min/1.73m² |
| **Partial Renal Response (PRR)** | UPCR < 1.0 g/g + eGFR >= 60 mL/min/1.73m² |
| **SLEDAI Renal** | Sum of haematuria + proteinuria + pyuria + casts (maximum 16 points) |

---

## Quick Start

```r
# Install mrgsolve
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
library(mrgsolve)
library(dplyr)

# Run the model
source("ln_mrgsolve_model.R")

# Run the Shiny app
library(shiny)
shiny::runApp("ln_shiny_app.R")
```

---

## Key Parameters

| Parameter | Value | Source |
|---------|-----|------|
| MMF oral bioavailability (F) | 0.94 | Djebli 2012 |
| HCQ volume of distribution | 800 L/kg | Tett 1993 |
| Voclosporin CL/F | 31.6 L/h | Ito 2022 |
| Belimumab CL | 2.8 mL/h/kg | Dall'Era 2010 |
| Plasma cell half-life | ~230 days | Literature estimate |
| IC_load -> podocyte injury | k = 0.03/day | AURORA calibration |

---

## Key Clinical Trials in the References

- **AURORA 1/2** (Rovin et al., Lancet 2021): demonstrated voclosporin's superiority — CRR 41% vs. 23%
- **BLISS-LN** (Furie et al., NEJM 2020): belimumab superiority — CRR 43% vs. 32%
- **ACCESS** (Appel et al., JASN 2009): MMF vs. CYC — confirmed non-inferiority
- **Euro-Lupus** (Houssiau et al., 2002): non-inferiority of low-dose CYC

---

*This model was generated by Claude Code Routine (CCR) on 2026-06-25.*
