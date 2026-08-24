# Primary Biliary Cholangitis (PBC) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Gastroenterology/Hepatobiliary

[![PBC QSP Model](pbc_qsp_model.png)](pbc_qsp_model.svg)

## Overview
Primary biliary cholangitis (PBC) is a chronic cholestatic liver disease in which autoimmune attack on the small intrahepatic bile ducts drives progressive cholestasis and liver fibrosis. It occurs overwhelmingly in middle-aged women (female:male = 10:1), with a prevalence of approximately 40 per 100,000. Antimitochondrial antibodies (AMA, particularly anti-PDC-E2) are positive in over 95% of cases, and CD4/CD8 T-cell and NK-cell-mediated injury to cholangiocytes is the core mechanism. UDCA is the standard first-line treatment, with FXR agonists (obeticholic acid) and PPAR agonists (elafibranor, seladelpar, bezafibrate) added in patients with an inadequate response to UDCA.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Autoimmune biliary injury pathway | AMA → PDC-E2 recognition, CD4/CD8/NK-cell infiltration | Destruction of small bile ducts, cholestasis |
| Bile acid toxicity pathway | Accumulation of hydrophobic bile acids → cell death, oxidative stress | Hepatocyte and cholangiocyte injury |
| FXR-FGF19 pathway | FXR activation → FGF19 release → suppressed bile acid synthesis | Reduced bile acid load |
| PPAR anti-inflammatory/antifibrotic pathway | PPARα/δ activation → reduced bile acid toxicity, suppressed inflammation | Reduced ALP/GGT, suppressed fibrosis |
| Hepatic fibrosis pathway | TGF-β stimulation of hepatic stellate cells (HSC) → collagen deposition | Cirrhosis, portal hypertension |
| Enterohepatic bile acid circulation | Intestinal bile acid reabsorption, enterohepatic circulation | Regulation of the bile acid pool |

## Drug Targets
- **UDCA (ursodeoxycholic acid)**: substitutes a hydrophilic bile acid, provides cytoprotection and immunomodulation — the first-line standard treatment
- **Obeticholic acid (OCA, FXR agonist)**: activates FXR → increases FGF19 → suppresses bile acid synthesis (POISE trial)
- **Elafibranor (PPARα/δ)**: reduces bile acid toxicity and inflammation (ELATIVE phase 3)
- **Seladelpar (PPARδ)**: regulates bile acid reabsorption, anti-inflammatory effect (RESPONSE phase 3)
- **Bezafibrate (PPARα)**: suppresses bile acid synthesis, normalises ALP (BEZURSO trial)

## Model Files
| File | Description |
|------|------|
| [pbc_qsp_model.dot](pbc_qsp_model.dot) | Graphviz mechanistic map source (approximately 181 nodes / 10 clusters) |
| [pbc_qsp_model.svg](pbc_qsp_model.svg) | SVG vector image (scalable) |
| [pbc_qsp_model.png](pbc_qsp_model.png) | PNG image (150 dpi) |
| [pbc_mrgsolve_model.R](pbc_mrgsolve_model.R) | mrgsolve ODE model (approximately 21 compartments / 7 treatment scenarios) |
| [pbc_shiny_app.R](pbc_shiny_app.R) | Shiny dashboard |
| [pbc_references.md](pbc_references.md) | References (approximately 57 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: gut/central PK compartments for UDCA, OCA, elafibranor, seladelpar, and bezafibrate, plus PD compartments for AMA, Th1 cells, cholangiocyte injury, bile acid toxicity, FGF19, ALP, bilirubin, GGT, fibrosis, and IgM
- **Key treatment scenarios**: ① no treatment, ② UDCA monotherapy (standard), ③ UDCA+OCA (POISE regimen), ④ UDCA+elafibranor (ELATIVE), ⑤ UDCA+seladelpar (RESPONSE), ⑥ UDCA+bezafibrate (BEZURSO), ⑦ triple therapy (UDCA+ELF+SEL)
- **Calibration/evidence**: calibrated with reference to the ALP response rates and GLOBE score improvements from the POISE (OCA), ELATIVE (elafibranor), RESPONSE (seladelpar), and BEZURSO (bezafibrate) phase 3 trials

## Shiny Dashboard
Comprises tabs for patient profile (UDCA responsiveness, PBC-40 symptom score, liver fibrosis stage), drug PK profile, bile acid/FXR-PPAR PD, liver biochemistry clinical endpoints (ALP, bilirubin, GGT), treatment scenario comparison (ALP normalisation rate, GLOBE score), and autoimmune/fibrosis biomarkers (AMA, IgM, Fibroscan).

## Usage
```r
library(mrgsolve)
mod <- mread("pbc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("pbc_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg pbc_qsp_model.dot -o pbc_qsp_model.svg
```

## References
For detailed citations, see [pbc_references.md](pbc_references.md) (approximately 57 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
