# Pseudogout (CPPD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Endocrine/Metabolic

[![CPPD QSP Model](../../../pseudogout/cppd_qsp_model.png)](../../../pseudogout/cppd_qsp_model.svg)

## Overview

Pseudogout is a crystal-induced arthritis caused by the deposition of calcium pyrophosphate dihydrate (CPPD) crystals in articular cartilage and synovium. It typically triggers sudden inflammatory attacks in large joints such as the knee and wrist, and its prevalence rises markedly after age 60. The core pathogenic mechanism is recognition of CPPD crystals by the NLRP3 inflammasome, which drives massive release of IL-1β and IL-18; colchicine, NSAIDs, and steroids are the principal treatment targets. Secondary CPPD is associated with metabolic abnormalities such as hyperparathyroidism, hypomagnesaemia, and haemochromatosis.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| CPPD crystal formation | Excess inorganic pyrophosphate (PPi) production via ANKH/ENPP1, transglutaminase 2 activation | Cartilage calcification, radiographic chondrocalcinosis |
| NLRP3 inflammasome activation | CPPD crystals → K+ efflux → assembly of the NLRP3/ASC/caspase-1 complex | Maturation and secretion of IL-1β/IL-18 |
| Neutrophil infiltration | IL-1β-induced CXCL8/C5a → massive neutrophil recruitment into the joint | Acute joint inflammation, pain |
| COX-2/prostaglandin pathway | Increased PGE2 production by macrophages and synoviocytes | Vasodilation, pain amplification |
| Matrix metalloproteinases (MMPs) | IL-1β/TNF → MMP-3/MMP-13 expression | Cartilage degradation, joint damage |
| Metabolic-abnormality-associated pathway | PTH excess → raised extracellular Ca²⁺/PPi | Promotes secondary CPPD |
| Chronic low-grade inflammation | M1 macrophages/TGF-β → fibrosis, chronic arthropathy | Degenerative joint damage |

## Drug Targets

- **Colchicine**: inhibits tubulin polymerisation → blocks neutrophil chemotaxis and NLRP3 activation; used for attack prevention and acute treatment
- **NSAIDs (naproxen, indomethacin)**: inhibit COX-1/2 → reduce PGE2; first-line treatment for acute attacks
- **Corticosteroids (prednisolone, intra-articular triamcinolone injection)**: broad anti-inflammatory action via NF-κB inhibition; an alternative when NSAIDs are contraindicated
- **IL-1 blockers (anakinra, canakinumab)**: directly neutralise IL-1β; used for recurrent or refractory attacks
- **Methotrexate/hydroxychloroquine**: suppress chronic recurrence (off-label)

## Model Files

| File | Description |
|------|------|
| [cppd_qsp_model.dot](../../../pseudogout/cppd_qsp_model.dot) | Graphviz mechanistic map source (100+ nodes / 12 clusters) |
| [cppd_qsp_model.svg](../../../pseudogout/cppd_qsp_model.svg) | SVG vector image (scalable) |
| [cppd_qsp_model.png](../../../pseudogout/cppd_qsp_model.png) | PNG image (150 dpi) |
| [cppd_mrgsolve_model.R](../../../pseudogout/cppd_mrgsolve_model.R) | mrgsolve ODE model (approximately 21 compartments / approximately 27 scenarios) |
| [cppd_shiny_app.R](../../../pseudogout/cppd_shiny_app.R) | Shiny dashboard |
| [cppd_references.md](cppd_references.md) | References (approximately 62 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: 2-compartment PK for colchicine/NSAID/steroid/IL-1Ra + a PD module running from NLRP3 activation → IL-1β production → neutrophil recruitment → joint inflammation index (JII), including a cartilage PPi accumulation compartment
- **Key treatment scenarios**: untreated natural course, colchicine monotherapy, NSAID monotherapy, intra-articular steroid injection, subcutaneous anakinra, single-dose canakinumab, colchicine + NSAID combination, prophylactic long-term low-dose colchicine, and others
- **Calibration/evidence**: parameters set from Dalbeth et al. *Lancet* 2021 (CPPD epidemiology), Ea et al. *Arthritis Rheum* NLRP3 mechanistic studies, CRESCENT trial data (canakinumab), and colchicine PK literature (Terkeltaub et al.)

## Shiny Dashboard

Comprises 6 tabs: (1) **Patient Profile** — sets age, comorbidities, and triggering factors; (2) **PK Profile** — blood/tissue concentration time courses for each drug; (3) **Inflammation Markers** — IL-1β, neutrophil, and CRP dynamics; (4) **Clinical Endpoints** — simulated joint pain score and attack frequency; (5) **Scenario Comparison** — comparison between monotherapy and combination treatment; (6) **Biomarkers** — serum PPi, uric acid, and synovial fluid crystal index.

## Usage

```r
library(mrgsolve)
mod <- mread("cppd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("cppd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg cppd_qsp_model.dot -o cppd_qsp_model.svg
```

## References

For detailed citations, see [cppd_references.md](cppd_references.md) (approximately 62 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
