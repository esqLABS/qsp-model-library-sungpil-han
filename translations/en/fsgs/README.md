# Focal Segmental Glomerulosclerosis (FSGS) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Renal/Urological

[![FSGS QSP Model](../../../fsgs/fsgs_qsp_model.png)](../../../fsgs/fsgs_qsp_model.svg)

## Overview

Focal segmental glomerulosclerosis (FSGS) is a major cause of nephrotic syndrome characterised by podocyte injury, accounting for approximately 35~50% of adult nephrotic syndrome. In primary FSGS, a circulating permeability factor (hypothesised to include suPAR, CLCF1, and others) stimulates podocyte surface receptors, causing foot process effacement and proteinuria. Untreated, the disease progresses to end-stage renal failure within a few years, and RAAS inhibition, steroids, calcineurin inhibitors, and rituximab are the principal treatment options. Sparsentan (DUPLEX trial) has recently been approved as a new targeted therapy.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Circulating permeability factor | CLCF1/suPAR → activates podocyte β3 integrin | Foot process effacement, proteinuria |
| Podocyte injury | Loss of nephrin, podocin, and synaptopodin → breakdown of the glomerular filtration barrier | Proteinuria ≥ 3.5 g/day |
| Glomerulosclerosis | TGF-β → mesangial cell activation → collagen deposition | Segmental glomerulosclerosis, reduced GFR |
| Complement pathway | C3b/MAC formation → podocyte cell death | Accelerated glomerular injury |
| RAAS overactivation | Ang II → raised intraglomerular pressure → worsened proteinuria | Accelerated decline in renal function |
| Chronic renal inflammation | Proteinuria → peritubular inflammation → interstitial fibrosis | Progressive decline in eGFR |
| Secondary FSGS | Obesity, solitary kidney, hypertension → glomerular hypertrophy and sclerosis | Proteinuria, chronic renal failure |

## Drug Targets

- **RAAS inhibitors (ACEi/ARB)**: reduce intraglomerular pressure, suppress proteinuria — the foundational treatment for all FSGS
- **Prednisolone (high-dose steroids)**: protects podocytes, reduces proteinuria — first line in primary FSGS
- **Calcineurin inhibitors (tacrolimus, ciclosporin)**: direct podocyte protection plus immunosuppression — for steroid-resistant disease
- **Rituximab**: B-cell depletion and direct podocyte stabilisation (via the sphingomyelin phosphodiesterase pathway)
- **Sparsentan**: a dual antagonist (AT1R + endothelin ETA receptor) — the DUPLEX trial

## Model Files

| File | Description |
|------|------|
| [fsgs_qsp_model.dot](../../../fsgs/fsgs_qsp_model.dot) | Graphviz mechanistic map source (approximately 143 nodes / 12 clusters) |
| [fsgs_qsp_model.svg](../../../fsgs/fsgs_qsp_model.svg) | SVG vector image (scalable) |
| [fsgs_qsp_model.png](../../../fsgs/fsgs_qsp_model.png) | PNG image (150 dpi) |
| [fsgs_mrgsolve_model.R](../../../fsgs/fsgs_mrgsolve_model.R) | mrgsolve ODE model (approximately 23 compartments / 6 treatment scenarios) |
| [fsgs_shiny_app.R](../../../fsgs/fsgs_shiny_app.R) | Shiny dashboard |
| [fsgs_references.md](../../../fsgs/fsgs_references.md) | References (approximately 67 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: disease-state compartments (circulating permeability factor, podocyte fraction, foot process effacement index, proteinuria, eGFR, sclerosis fraction, TGF-β, complement activation, glomerular inflammation) + drug PK compartments (prednisolone, tacrolimus, rituximab, sparsentan PK)
- **Key treatment scenarios**: ① natural course (no treatment), ② prednisolone monotherapy, ③ prednisolone + tacrolimus (steroid-resistant protocol), ④ prednisolone + tacrolimus + rituximab (refractory FSGS), ⑤ sparsentan monotherapy (based on the DUPLEX trial), ⑥ full combination therapy
- **Calibration/evidence**: parameters referenced from the DUPLEX trial (sparsentan), the Cattran FSGS international registry data, the D'Agati classification, and the FSGS natural history literature

## Shiny Dashboard

Comprises 6 tabs: ① Patient Profile (sets baseline proteinuria, eGFR, suPAR, and FSGS subtype), ② PK tab (blood concentration of each drug, tacrolimus trough monitoring), ③ Podocyte/renal PD tab (podocyte fraction, proteinuria, and eGFR trends), ④ Clinical Endpoints (complete/partial remission, time to 50% eGFR decline), ⑤ Scenario Comparison (simultaneous comparison of 6 treatment strategies), ⑥ Biomarkers (TGF-β, complement, sclerosis fraction, and tubular injury index trends).

## Usage

```r
library(mrgsolve)
mod <- mread("fsgs_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("fsgs_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg fsgs_qsp_model.dot -o fsgs_qsp_model.svg
```

## References

For detailed citations, see [fsgs_references.md](../../../fsgs/fsgs_references.md) (approximately 67 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
