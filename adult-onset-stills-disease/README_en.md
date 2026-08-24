# Adult-Onset Still's Disease (AOSD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Autoimmune/Rheumatic

[![AOSD QSP Model](aosd_qsp_model.png)](aosd_qsp_model.svg)

## Overview

Adult-onset Still's disease (AOSD) is a systemic autoinflammatory disease of unknown cause, a rare condition with an incidence of approximately 1.5–2.2 per million population. Its core pathogenic mechanism is excess IL-1β/IL-18 production driven by NLRP3 inflammasome activation, leading to an IL-6/TNF-α/IFN-γ cytokine storm, with systemic macrophage activation as a hallmark feature. Clinically, the triad of high fever (≥ 39°C, occurring 1–2 times daily), a transient salmon-coloured skin rash, and arthritis is characteristic, together with marked hyperferritinaemia (often above 10,000 ng/mL). Macrophage activation syndrome (MAS), a serious, life-threatening complication, occurs in approximately 10–15% of cases. IL-1 blockers (anakinra, canakinumab) and the IL-6 blocker (tocilizumab) are the standard biologic therapies.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| NLRP3 inflammasome | Pathogen/danger signal → caspase-1 → pro-IL-1β cleavage | Excess IL-1β production, fever |
| IL-18 axis | Macrophage IL-18 secretion → NK-cell/CD8⁺ T-cell activation | Induces MAS, increases IFN-γ |
| IL-6 pathway | STAT3 activation → acute-phase protein synthesis | Sharp rise in CRP/ferritin, anaemia |
| Macrophage activation | M1 polarisation, enhanced phagocytosis | Hepatosplenomegaly, cytopenia (MAS) |
| TNF-α signalling | NF-κB → expression of multiple inflammatory genes | Systemic inflammation, joint damage |
| Ferritin excess | Excess iron-storage protein → direct immune dysregulation | Serves as a diagnostic biomarker |
| NK-cell dysfunction | IL-18-induced NK-cell exhaustion | Loss of cytotoxic immune surveillance |

## Drug Targets

- **Anakinra (IL-1Ra)**: an IL-1 receptor antagonist; excellent early response in the acute phase, short half-life requiring daily SC dosing
- **Canakinumab (anti-IL-1β antibody)**: SC dosing every 4 weeks; effective in AOSD and MAS prevention
- **Tocilizumab (anti-IL-6R antibody)**: blocks IL-6 signalling; normalises CRP/ferritin, improves arthritis
- **Corticosteroids**: first-line therapy; effective for acute suppression but limited by long-term side effects and limited MAS-prevention effect
- **Tofacitinib (JAK inhibitor)**: growing evidence for blocking the IFN-γ pathway in refractory AOSD and MAS

## Model Files

| File | Description |
|------|------|
| [aosd_qsp_model.dot](aosd_qsp_model.dot) | Graphviz mechanistic map source (approximately 161 nodes / 13 clusters) |
| [aosd_qsp_model.svg](aosd_qsp_model.svg) | SVG vector image (scalable) |
| [aosd_qsp_model.png](aosd_qsp_model.png) | PNG image (150 dpi) |
| [aosd_mrgsolve_model.R](aosd_mrgsolve_model.R) | mrgsolve ODE model (approximately 21 compartments / multiple treatment scenarios) |
| [aosd_shiny_app.R](aosd_shiny_app.R) | Shiny dashboard |
| [aosd_references.md](aosd_references.md) | References (approximately 67 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: anakinra (SC, 2 compartments) + canakinumab (SC, 2 compartments) + tocilizumab (IV, 2 compartments) + corticosteroid and tofacitinib PK compartments + PD compartments for IL-1β, IL-6, IL-18, IFN-γ, TNF-α, ferritin, CRP, activated macrophages, NK cells, AOSD activity, and MAS risk
- **Key treatment scenarios**: ① no treatment, ② corticosteroid monotherapy, ③ anakinra + steroid, ④ canakinumab, ⑤ tocilizumab, ⑥ tofacitinib (refractory disease), ⑦ combination (steroid + biologic)
- **Calibration/evidence**: parameters referenced from the SOBI ANAKIN clinical study (anakinra), the Yamaguchi classification criteria, and the Fautrel criteria

## Shiny Dashboard

The dashboard comprises a patient profile tab (fever pattern, skin rash, baseline ferritin), a drug PK profile tab, an IL-1/IL-6/IL-18/IFN-γ cytokine dynamics tab, a clinical endpoints tab (fever, arthritis score, ferritin), a treatment scenario comparison tab, and an MAS risk biomarker tab.

## Usage

```r
library(mrgsolve)
mod <- mread("aosd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("aosd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg aosd_qsp_model.dot -o aosd_qsp_model.svg
```

## References

For detailed citations, see [aosd_references.md](aosd_references.md) (approximately 67 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
