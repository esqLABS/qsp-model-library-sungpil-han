# Relapsing Polychondritis (RPC) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Autoimmune/Rheumatic

[![RPC QSP Model](../../../relapsing-polychondritis/rpc_qsp_model.png)](../../../relapsing-polychondritis/rpc_qsp_model.svg)

## Overview

Relapsing polychondritis (RPC) is a rare autoimmune disease that targets the cartilage and elastic tissue of the ears, nose, trachea, bronchi, joints, and eyes, with an annual incidence of approximately 3.5 per million. Autoantibodies against type II, IX, and XI collagen and matrilin-1, together with CD4+ T-cell responses, destroy the cartilage matrix, and airway cartilage involvement can be life-threatening through tracheal stenosis. Corticosteroids are the core treatment in the acute phase, with immunosuppressants such as dapsone, methotrexate, and azathioprine used to suppress relapse, and anti-TNF, anti-IL-6, and anti-CD20 biologics used in refractory cases.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Collagen autoantigen recognition | Type II/IX/XI collagen antigens → autoantibody (IgG) + CD4+ Th1 activation | Initiation of immune-mediated cartilage damage |
| Th1/Th17 infiltration | IFN-γ, IL-17 → recruitment of perichondrial macrophages and neutrophils | Cartilage matrix degradation, acute inflammation |
| MMP-mediated cartilage degradation | TNF, IL-1β → secretion of MMP-1/3/13, ADAMTS-5 | Loss of cartilage matrix (proteoglycan, collagen) |
| Airway cartilage involvement | Subglottic and tracheal cartilage inflammation → fibrosis, softening | Airway stenosis, risk of respiratory failure |
| Inner-ear/ocular vascular inflammation | Inner-ear vasculitis → sensorineural hearing loss; scleritis, uveitis | Sensory organ damage |
| Complement activation | Autoantibodies → C1q, C3a, C5a activation → neutrophil influx into cartilage | Accelerated cartilage destruction |
| Fibrosis/remodelling | TGF-β → fibroblast activation → airway fibrosis | Irreversible airway deformation |

## Drug Targets

- **Corticosteroids (prednisolone)**: NF-κB inhibition → broad anti-inflammatory action; first-line therapy for acute flares, high dose essential in airway involvement
- **Dapsone**: inhibits neutrophil migration and MPO; suppresses mild relapses, prevents cutaneous and auricular cartilage flares
- **Methotrexate**: folate antagonism → suppressed T-cell proliferation; used as a steroid-sparing agent
- **Azathioprine/mycophenolate**: purine synthesis inhibition; maintenance immunosuppression
- **Anti-TNF agents (infliximab, etanercept)**: TNF blockade → reduced MMP and IL-1β; for refractory airway and joint involvement
- **Anti-IL-6R (tocilizumab)**: blocks IL-6 signalling; normalises systemic inflammatory markers
- **Anti-CD20 (rituximab)**: B-cell depletion → reduced autoantibodies; for severe relapsing cases

## Model Files

| File | Description |
|------|------|
| [rpc_qsp_model.dot](../../../relapsing-polychondritis/rpc_qsp_model.dot) | Graphviz mechanistic map source (approximately 100+ nodes / 10 clusters) |
| [rpc_qsp_model.svg](../../../relapsing-polychondritis/rpc_qsp_model.svg) | SVG vector image (scalable) |
| [rpc_qsp_model.png](../../../relapsing-polychondritis/rpc_qsp_model.png) | PNG image (150 dpi) |
| [rpc_mrgsolve_model.R](../../../relapsing-polychondritis/rpc_mrgsolve_model.R) | mrgsolve ODE model (approximately 21 compartments / approximately 41 scenarios) |
| [rpc_shiny_app.R](../../../relapsing-polychondritis/rpc_shiny_app.R) | Shiny dashboard |
| [rpc_references.md](../../../relapsing-polychondritis/rpc_references.md) | References (approximately 63 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: PK compartments for steroids/dapsone/methotrexate/biologics + autoantibody titre dynamics, Th1 activation index, MMP production compartment, cartilage thickness damage compartment, and an airway diameter stenosis prediction module
- **Key treatment scenarios**: untreated natural course, prednisolone monotherapy, dapsone, prednisolone + methotrexate, prednisolone + azathioprine, infliximab, tocilizumab, rituximab, combined steroid + biologic, and many others
- **Calibration/evidence**: based on clinical cohort data related to the McAdam diagnostic criteria, the Mathian et al. tocilizumab case series, and the Moulis et al. French national cohort literature

## Shiny Dashboard

Comprises 6 tabs: (1) **Patient Profile** — set the site of involvement (ear, nose, airway, joint, eye) and autoantibody status; (2) **PK Profile** — blood concentrations of steroids/immunosuppressants/biologics; (3) **Key PD Measures** — dynamics of autoantibody, MMP, and IL-1β suppression; (4) **Clinical Endpoints** — auricular flare frequency, airway diameter, hearing, ocular disease activity; (5) **Scenario Comparison** — comparison of flare frequency and long-term damage across treatment strategies; (6) **Biomarkers** — CRP, ESR, type II collagen antibody, pulmonary function (FEV1).

## Usage

```r
library(mrgsolve)
mod <- mread("rpc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("rpc_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg rpc_qsp_model.dot -o rpc_qsp_model.svg
```

## References

For detailed citations, see [rpc_references.md](../../../relapsing-polychondritis/rpc_references.md) (approximately 63 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
