# Ankylosing Spondylitis (AS) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Autoimmune/Rheumatic

[![AS QSP Model](as_qsp_model.png)](as_qsp_model.svg)

## Overview

Ankylosing spondylitis (AS) is a chronic inflammatory arthritis mainly
affecting the spine and sacroiliac joints, with an HLA-B27 positivity
rate of about 85-95% and a population prevalence of 0.1-0.5%. The core of
its pathogenesis is the HLA-B27-related unfolded protein response (UPR)
and IL-23/IL-17-axis enthesitis. Gut microbial dysbiosis and intestinal
mucosal inflammation promote excess IL-23 production, which activates
Th17 cells in entheseal tissue to secrete IL-17A. Chronic inflammation
leads to a dual outcome: on one hand bone erosion (damage), and on the
other, ectopic bone formation (new bone, syndesmophytes) driven by
TNF-induced activation of Wnt/BMP signalling. TNF inhibitors and IL-17A
inhibitors are effective at improving BASDAI/ASDAS, and JAK inhibitors are
also used.

## Key Pathophysiological Pathways

| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| HLA-B27/UPR | Misfolded HLA-B27 heavy-chain dimerisation → ER stress → IL-23 | Activation of the Th17 axis |
| IL-23/IL-17 axis | IL-23 → RORγt⁺ T cells → secretion of IL-17A/F | Enthesitis and spinal inflammation |
| TNF-α signalling | TNF → NF-κB → activation of osteoclasts and chondrocytes | Bone erosion, disc damage |
| New bone formation | Wnt/DKK1 balance → BMP → fibroblast ossification | Syndesmophytes, spinal ankylosis |
| Gut-spine axis | Gut microbial dysbiosis → increased mucosal IL-23 | Triggers spondyloarthritis |
| RANKL/OPG | Inflammation → RANKL↑, OPG↓ → osteoclast proliferation | Bone erosion, osteoporosis |
| CRP/disease activity | Systemic inflammation → elevated ASDAS/BASDAI | Reduced quality of life, functional loss |

## Key Drug Targets

- **TNF inhibitors (adalimumab, etanercept, certolizumab, etc.)**: about 40-50% achieve BASDAI50 (50% improvement in BASDAI); excellent suppression of imaging-detected inflammation, incomplete suppression of new bone formation
- **IL-17A inhibitors (secukinumab, ixekizumab)**: at least as effective as TNF inhibitors in AS; caution needed with concomitant Crohn disease
- **JAK inhibitors (tofacitinib, upadacitinib)**: oral administration, effective in patients refractory to TNF/IL-17 blockade
- **IL-23 inhibitors (risankizumab, guselkumab)**: effective in the peripheral-arthritis-predominant subtype; efficacy in axial spondylitis still under study
- **NSAIDs**: first-line symptomatic treatment; long-term use may suppress new bone formation

## Model Files

| File | Description |
|------|------|
| [as_qsp_model.dot](as_qsp_model.dot) | Graphviz mechanistic map source (about 191 nodes / 10 clusters) |
| [as_qsp_model.svg](as_qsp_model.svg) | SVG vector image (zoomable) |
| [as_qsp_model.png](as_qsp_model.png) | PNG image (150 dpi) |
| [as_mrgsolve_model.R](as_mrgsolve_model.R) | mrgsolve ODE model (about 22 compartments / multiple treatment scenarios) |
| [as_shiny_app.R](as_shiny_app.R) | Shiny dashboard |
| [as_references.md](as_references.md) | References (about 56, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: adalimumab/etanercept (SC PK) + secukinumab (SC 2-compartment) + tofacitinib/upadacitinib (oral) + NSAID PK compartments plus TNF, IL-17A, IL-23, IL-6, CRP, RANKL, OPG, osteoclast, bone-erosion, new-bone/mSASSS, and a composite BASDAI disease-activity PD compartment
- **Key treatment scenarios**: (1) no treatment, (2) NSAID alone, (3) adalimumab, (4) etanercept, (5) secukinumab, (6) tofacitinib, (7) switch to IL-17 blockade after TNF failure
- **Calibration/basis**: referenced against the MEASURE trials (secukinumab), ATLAS/ABILITY adalimumab/certolizumab data, and parameters for a model of mSASSS progression

## Shiny Dashboard

Structured into a patient profile tab (HLA-B27, baseline BASDAI, imaging
findings), a drug PK and TNF/IL-17 inhibition tab, a disease activity
(BASDAI/ASDAS) trajectory tab, a structural damage (mSASSS, new bone) tab,
a treatment scenario comparison tab, and a biomarker (CRP, RANKL, mSASSS)
tab.

## Usage

```r
library(mrgsolve)
mod <- mread("as_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("as_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg as_qsp_model.dot -o as_qsp_model.svg
```

## References

For full citations, see [as_references.md](as_references.md) (about 56 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
