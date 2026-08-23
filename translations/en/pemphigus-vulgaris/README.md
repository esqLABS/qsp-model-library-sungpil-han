# Pemphigus Vulgaris (PV) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Dermatology

[![PV QSP Model](../../../pemphigus-vulgaris/pv_qsp_model.png)](../../../pemphigus-vulgaris/pv_qsp_model.svg)

## Overview

Pemphigus vulgaris (PV) is a rare autoimmune blistering disease in which autoantibodies (mainly IgG4) against the epidermal intercellular adhesion proteins desmoglein-3 (Dsg3), and in some cases desmoglein-1 (Dsg1), disrupt keratinocyte cell-cell adhesion, causing intraepidermal blistering and mucosal erosions. The estimated worldwide annual incidence is approximately 1–5 per million, and the disease is genetically more common in Ashkenazi Jewish and Mediterranean populations. Left untreated, extensive skin loss and sepsis carry a high mortality, so steroids and rituximab are the core therapies. Studies expanding the indication of efgartigimod (anti-FcRn) are currently underway.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| TFH-GC B-cell axis | TFH (IL-21) → GC B cells → memory B cells/plasma cells → anti-Dsg3 IgG | Sustained antibody production, relapse |
| Anti-Dsg3 IgG4 binding | IgG4 antibody → destroys Dsg3 target, activates p38-MAPK → desmosome breakdown | Intraepidermal blistering (acantholysis) |
| Complement activation | Dsg3-IgG1 complex → C1q → C3b → amplified inflammation | Accelerated epidermal damage |
| Impaired Treg function | CD4+Foxp3+ Treg↓ → failure to suppress TFH/GC B cells | Breakdown of self-tolerance |
| FcRn-IgG recycling | FcRn prevents IgG degradation → prolongs anti-Dsg3 antibody half-life | Sustained high antibody titre |
| Keratinocyte signalling | Loss of Dsg3 → EGFR/p38 → cell shrinkage, apoptosis | Blister formation, expanding erosion |
| Long-lived plasma cells (LLPC) | Bone marrow LLPC → produces antibody resistant to steroids/rituximab | Relapse, treatment resistance |

## Drug Targets

- **High-dose prednisolone**: broad immunosuppression; suppresses TFH/GC B cells, expands Treg → first-line therapy
- **Rituximab**: anti-CD20 B-cell depletion; superior CR rate versus steroids alone in the RITUX3 trial
- **Mycophenolate mofetil (MMF)**: inhibits purine synthesis → blocks B-cell/T-cell proliferation; used for maintenance immunosuppression
- **Efgartigimod**: anti-FcRn → rapidly lowers circulating IgG (including anti-Dsg3 antibody) levels; BALLAD study
- **IV methylprednisolone pulse**: induces rapid remission in severe acute disease; followed by tapering after high-dose steroids

## Model Files

| File | Description |
|------|------|
| [pv_qsp_model.dot](../../../pemphigus-vulgaris/pv_qsp_model.dot) | Graphviz mechanistic map source (approximately 159 nodes / 12 clusters) |
| [pv_qsp_model.svg](../../../pemphigus-vulgaris/pv_qsp_model.svg) | SVG vector image (scalable) |
| [pv_qsp_model.png](../../../pemphigus-vulgaris/pv_qsp_model.png) | PNG image (150 dpi) |
| [pv_mrgsolve_model.R](../../../pemphigus-vulgaris/pv_mrgsolve_model.R) | mrgsolve ODE model (approximately 23 compartments / 6 treatment scenarios) |
| [pv_shiny_app.R](../../../pemphigus-vulgaris/pv_shiny_app.R) | Shiny dashboard |
| [pv_references.md](../../../pemphigus-vulgaris/pv_references.md) | References (approximately 40 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: prednisolone (2 oral compartments + peripheral), rituximab IV (2 compartments), MMF oral (1 compartment), IVIg (1 compartment), efgartigimod (1 compartment + FcRn); immature B cells, GC B cells, memory B cells, SLPC, LLPC, TFH, Treg, anti-Dsg3 IgG, Dsg3 expression, blister/PDAI score, complement activation, and cumulative corticosteroid-related bone loss
- **Key treatment scenarios**: ① high-dose CS monotherapy (existing standard) ② rituximab + low-dose CS (RITUX3) ③ MMF + moderate CS ④ rituximab + MMF combination ⑤ IV methylprednisolone pulse + rituximab (severe disease) ⑥ efgartigimod + low-dose CS
- **Calibration/evidence**: parameters referenced from the RITUX3 (rituximab + low-dose CS vs. high-dose CS) and BALLAD (efgartigimod) clinical trial data

## Shiny Dashboard

Comprises 6 tabs: ① **Patient Profile** (sets baseline anti-Dsg3 titre, PDAI score, and mucosal involvement), ② **PK** (plasma drug concentration and FcRn occupancy, rituximab CD20 occupancy), ③ **Key PD Metrics** (anti-Dsg3 IgG, GC B-cell, and LLPC trends), ④ **Clinical Endpoints** (PDAI, CR/PR rate, time to relapse), ⑤ **Scenario Comparison** (direct comparison of 6 treatment strategies), ⑥ **Biomarkers** (anti-Dsg3 IgG4, IgG1, and cumulative steroid toxicity trends).

## Usage

```r
library(mrgsolve)
mod <- mread("pv_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("pv_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg pv_qsp_model.dot -o pv_qsp_model.svg
```

## References

For detailed citations, see [pv_references.md](../../../pemphigus-vulgaris/pv_references.md) (approximately 40 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
