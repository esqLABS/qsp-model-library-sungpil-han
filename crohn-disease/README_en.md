# Crohn's Disease (CD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Gastroenterology/Hepatobiliary

[![CD QSP Model](cd_qsp_model.png)](cd_qsp_model.svg)

## Overview

Crohn's disease is a transmural chronic inflammatory bowel disease that can occur anywhere in the digestive tract, including the small and large intestine, affecting more than 3 million people worldwide and occurring predominantly in young adults. An excessive Th1/Th17-driven immune response, acting through the TNF-α, IL-12/23, and IL-17 cytokine network, causes increased intestinal permeability and the formation of transmural granulomas. Anti-TNF biologics (infliximab, adalimumab) form the foundation of biological therapy, with anti-IL-12/23 (ustekinumab), anti-integrin (vedolizumab), and JAK inhibitors (upadacitinib) established as additional treatment options.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Th1 overactivation | IL-12 → Th1 differentiation → IFN-γ, TNF-α secretion | Transmural granulomatous inflammation |
| Th17 pathway | IL-23 → Th17 → IL-17, IL-22 secretion | Barrier damage, mucosal immune dysregulation |
| TNF-α-centred inflammation | TNF-α → NF-κB → cytokine cascade | Ulceration, stricture, fistula formation |
| Impaired barrier function | Reduced tight-junction proteins (ZO-1, claudins) | Bacterial translocation, endotoxaemia |
| Bone marrow migration/neutrophil infiltration | IL-8, MIP-1α → neutrophil recruitment | Mucosal damage, raised faecal calprotectin |
| Mesenteric adipose tissue inflammation | Creeping fat → adipokine/cytokine secretion | Determines disease location, promotes stricture |
| Bone metabolism abnormality | Steroids plus chronic inflammation → reduced BMD | Osteoporosis risk |

## Drug Targets

- **Anti-TNF antibodies (infliximab, adalimumab, certolizumab)**: neutralise TNF-α → inhibit NF-κB, promote mucosal healing
- **Ustekinumab (anti-IL-12/23)**: blocks the p40 subunit → suppresses Th1/Th17 differentiation (UNIFI/CERTIFI)
- **Vedolizumab (anti-α4β7 integrin)**: blocks gut-selective lymphocyte homing (GEMINI trial)
- **Upadacitinib (JAK1 inhibitor)**: suppresses IL-6/IL-12/IL-23 signalling → pan-cytokine suppression (U-EXCEL)
- **Azathioprine/6-MP (thiopurines)**: inhibit purine synthesis → suppress proliferating lymphocytes; combined with biologics to suppress antibody formation

## Model Files

| File | Description |
|------|------|
| [cd_qsp_model.dot](cd_qsp_model.dot) | Graphviz mechanistic map source (approximately 392 nodes / 10 clusters) |
| [cd_qsp_model.svg](cd_qsp_model.svg) | SVG vector image (scalable) |
| [cd_qsp_model.png](cd_qsp_model.png) | PNG image (150 dpi) |
| [cd_mrgsolve_model.R](cd_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / 8 treatment scenarios) |
| [cd_shiny_app.R](cd_shiny_app.R) | Shiny dashboard |
| [cd_references.md](cd_references.md) | References (approximately 61 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: PK compartments for infliximab, adalimumab, ustekinumab, vedolizumab, thioguanine, prednisolone, and upadacitinib (10 compartments) + PD compartments for TNF-α, IL-12/23, IL-17, Th17, Th1, Treg, neutrophils, mucosal inflammation index (MI), CRP, faecal calprotectin (FC), bone mineral density, and haemoglobin (12 compartments)
- **Key treatment scenarios**: (1) no treatment, (2) steroid monotherapy, (3) infliximab induction-maintenance, (4) adalimumab SC, (5) ustekinumab IV→SC, (6) vedolizumab IV, (7) azathioprine + infliximab, (8) oral upadacitinib
- **Calibration/evidence**: infliximab PK based on Ng CM et al. (Clin Pharmacokinet 2010); clinical response rates based on the ACCENT I/II, CHARM, and UNIFI trial data

## Shiny Dashboard

Comprises 6 tabs: (1) Patient Profile — sets lesion location, disease severity, and immunogenicity risk; (2) PK tab — plasma concentration time series for biologics and the small molecule; (3) Key PD Metrics — TNF-α, IL-12/23, Th1/Th17/Treg, CRP; (4) Clinical Endpoints — CDAI, faecal calprotectin, mucosal healing; (5) Scenario Comparison — 1-year outcomes across 8 treatment strategies; (6) Biomarkers — CRP, faecal calprotectin, haemoglobin, bone mineral density, and drug concentration TDM

## Usage

```r
library(mrgsolve)
mod <- mread("cd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("cd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg cd_qsp_model.dot -o cd_qsp_model.svg
```

## References

For detailed citations, see [cd_references.md](cd_references.md) (approximately 61 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
