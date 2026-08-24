# Polymyositis (PM) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Autoimmune/Rheumatic

[![PM QSP Model](pm_qsp_model.png)](pm_qsp_model.svg)

## Overview
Polymyositis (PM) is an autoimmune myositis in which CD8+ cytotoxic T cells directly attack muscle fibres, occurring in approximately 1~5 per 100,000 population and more common in women. MHC-I overexpression and CD8+ T-cell perforin/granzyme B-mediated invasion of muscle fibres are the principal mechanisms of tissue damage, and myositis-specific autoantibodies such as anti-Jo-1 (an antisynthetase antibody) serve as diagnostic and prognostic markers. Proximal muscle weakness, elevated serum CK, electromyographic abnormalities, and CD8 T-cell infiltration on muscle biopsy form the diagnostic criteria. High-dose steroids are the first-line treatment, with steroid-sparing immunosuppressants (MTX, AZA) added, and in some cases rituximab or a JAK inhibitor.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| MHC-I expression pathway | IFN-β/α stimulation → STAT1 activation → MHC-I overexpression | Enables CD8 T-cell recognition of muscle fibres |
| CD8 cytotoxic pathway | Perforin/granzyme B, TRAIL/FasL | Muscle fibre necrosis, CK release |
| Th1 inflammatory pathway | IFN-γ, TNF-α, IL-6 | Amplified intramuscular inflammation |
| B-cell/autoantibody pathway | MSAs such as anti-Jo-1/anti-Mi-2, plasma cells | Antisynthetase syndrome, ILD risk |
| Impaired muscle regeneration pathway | Chronic inflammation → satellite cell depletion | Irreversible muscle weakness |
| ILD-associated pathway | IFN, TGF-β, IL-6 → pulmonary fibrosis | Interstitial lung disease in anti-Jo-1-positive patients |

## Drug Targets
- **Prednisolone**: broad anti-inflammatory/immunosuppressive action — inhibits NF-κB/AP-1, first-line therapy
- **Methotrexate (MTX)**: antagonises folate metabolism, accumulates polyglutamated metabolites — a steroid-sparing agent
- **Azathioprine (AZA)/6-TGN**: suppresses purine synthesis, maintains remission
- **Rituximab (RTX)**: an anti-CD20 monoclonal antibody → depletes B cells/plasma cells, for refractory MSA-positive (e.g. anti-Jo-1) PM
- **JAK inhibitors**: inhibit JAK1/2 → block IFN signalling → suppress MHC-I expression

## Model Files
| File | Description |
|------|------|
| [pm_qsp_model.dot](pm_qsp_model.dot) | Graphviz mechanistic map source (approximately 193 nodes / 12 clusters) |
| [pm_qsp_model.svg](pm_qsp_model.svg) | SVG vector image (scalable) |
| [pm_qsp_model.png](pm_qsp_model.png) | PNG image (150 dpi) |
| [pm_mrgsolve_model.R](pm_mrgsolve_model.R) | mrgsolve ODE model (approximately 28 compartments / 6 treatment scenarios) |
| [pm_shiny_app.R](pm_shiny_app.R) | Shiny dashboard |
| [pm_references.md](pm_references.md) | References (approximately 53 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: PK compartments for prednisolone, MTX, 6-TGN, RTX (2 compartments + CD20 binding), IVIG, and JAK inhibitor, plus PD compartments for naive/effector CD8, CD4 Th1, B cells, plasma cells, autoantibody, IFN-γ, TNF-α, IL-6, MHC-I, muscle inflammation, CK, and MMT8 (muscle strength)
- **Key treatment scenarios**: ① no treatment, ② prednisolone monotherapy, ③ prednisolone+MTX, ④ prednisolone+AZA, ⑤ rituximab+prednisolone, ⑥ JAK inhibitor+prednisolone
- **Calibration/evidence**: MMT8 muscle strength score and CK normalisation time calibrated from the RIM trial (rituximab in PM/DM) and IMACS international network clinical data

## Shiny Dashboard
Comprises tabs for patient profile (MSA status, muscle strength score, concurrent ILD), drug PK (including RTX CD20 saturation), immune PD (CD8, B cells, cytokines), muscle clinical endpoints (MMT8, CK), treatment scenario comparison (remission rate, steroid dose), and biomarkers (IFN-γ, autoantibodies).

## Usage
```r
library(mrgsolve)
mod <- mread("pm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("pm_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg pm_qsp_model.dot -o pm_qsp_model.svg
```

## References
For detailed citations, see [pm_references.md](pm_references.md) (approximately 53 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
