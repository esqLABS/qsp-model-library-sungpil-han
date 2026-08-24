# Psoriatic Arthritis (PsA) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Autoimmune/Rheumatic

[![PsA QSP Model](psa_qsp_model.png)](psa_qsp_model.svg)

## Overview

Psoriatic arthritis (PsA) is a chronic inflammatory arthritis that
accompanies or precedes psoriatic skin lesions, occurring in about
25-30% of psoriasis patients, with a worldwide prevalence of 0.1-0.25%.
Its core pathogenic mechanism is the TNF-α and IL-23/IL-17 axis, which
simultaneously mediates enthesitis, synovial arthritis, dactylitis, spinal
inflammation (spondylitis), and nail psoriasis. The IL-23p19, IL-17A, TNF,
and JAK pathways are the main treatment targets, and biologics (anti-TNF,
IL-17 inhibitors, IL-12/23 and IL-23p19 inhibitors) together with JAK
inhibitors (tofacitinib, upadacitinib, filgotinib) make up the standard of
care.

## Key Pathophysiological Pathways

| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| IL-23/Th17 axis | DC/macrophage IL-23 → Th17 differentiation → secretion of IL-17A, IL-17F, IL-22 | Cutaneous psoriatic plaques, joint inflammation |
| TNF-α pathway | Synoviocyte/macrophage TNF → NF-κB → inflammatory cytokine cascade | Synovial thickening, bone erosion, elevated CRP |
| Enthesitis pathway | Mechanical stress + IL-17/TNF → activation of enthesis fibroblasts | Enthesis pain, imbalance of bone formation and resorption |
| JAK-STAT signalling | IL-6/IL-2/IFN-γ → JAK1/3-STAT3/STAT1 activation | Synovial proliferation, immune cell differentiation |
| RANKL/OPG pathway | TNF, IL-17 → osteoclast activation → bone erosion | Progressive radiographic bone damage |
| IL-12/Th1 axis | IL-12 → Th1 → IFN-γ → M1 macrophage polarisation | Chronic granulomatous inflammation, systemic symptoms |
| Keratinocyte hyperactivation | IL-17, IL-22 → keratinocyte proliferation and CXCL secretion | Formation of cutaneous psoriatic plaques |

## Key Drug Targets

- **Anti-TNF agents (adalimumab, etanercept, certolizumab)**: direct neutralisation of TNF-α; suppresses joint, skin, and radiographic progression
- **IL-17A inhibitors (secukinumab, ixekizumab)**: block IL-17A; excellent for skin, enthesitis, and axial lesions
- **IL-23p19 inhibitors (guselkumab, risankizumab, tildrakizumab)**: selective IL-23 inhibition; induces durable remission
- **IL-12/23p40 inhibitor (ustekinumab)**: dual blockade of IL-12 and IL-23; controls both skin and joint lesions simultaneously
- **JAK inhibitors (tofacitinib, upadacitinib)**: suppresses intracellular JAK-STAT signalling; broad cytokine blockade
- **DMARDs (methotrexate, sulfasalazine, leflunomide)**: mainly for peripheral joints; also partially effective for skin lesions

## Model Files

| File | Description |
|------|------|
| [psa_qsp_model.dot](psa_qsp_model.dot) | Graphviz mechanistic map source (about 100+ nodes / 15 clusters) |
| [psa_qsp_model.svg](psa_qsp_model.svg) | SVG vector image (zoomable) |
| [psa_qsp_model.png](psa_qsp_model.png) | PNG image (150 dpi) |
| [psa_mrgsolve_model.R](psa_mrgsolve_model.R) | mrgsolve ODE model (about 25 compartments / about 21 scenarios) |
| [psa_shiny_app.R](psa_shiny_app.R) | Shiny dashboard |
| [psa_references_en.md](psa_references_en.md) | References (about 41, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: 2-compartment PK for anti-TNF/anti-IL-17/anti-IL-23/JAK inhibitor drugs plus a TNF/IL-17A/IL-23 production-elimination PD module, a synovial inflammation index (SII), a PASI score compartment, and a joint-damage-progression (mTSS) compartment
- **Key treatment scenarios**: untreated natural history, methotrexate, adalimumab, secukinumab, guselkumab, ustekinumab, tofacitinib, adalimumab + methotrexate combination, and others
- **Calibration/basis**: parameters calibrated against data from the FUTURE trials (secukinumab), DISCOVER (guselkumab), PSUMMIT (ustekinumab), OPAL BROADEN (tofacitinib), and ADEPT (adalimumab)

## Shiny Dashboard

Structured into 6 tabs: (1) **Patient profile** — sets joint subtype,
skin involvement, comorbidities; (2) **PK profile** — blood concentration
over time for biologics/small molecules; (3) **Key PD measures** — IL-17A,
TNF, IL-23 inhibition rate; (4) **Clinical endpoints** — ACR20/50/70, PASI
75/90/100, mTSS change; (5) **Scenario comparison** — comparison of
monotherapy and combination treatments; (6) **Biomarkers** — CRP, ESR,
enthesis ultrasound score.

## Usage

```r
library(mrgsolve)
mod <- mread("psa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("psa_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg psa_qsp_model.dot -o psa_qsp_model.svg
```

## References

For full citations, see [psa_references_en.md](psa_references_en.md) (about 41 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
