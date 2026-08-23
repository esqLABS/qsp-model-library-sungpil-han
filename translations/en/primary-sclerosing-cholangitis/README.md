# Primary Sclerosing Cholangitis (PSC) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology/Hepatobiliary

[![PSC QSP Model](../../../primary-sclerosing-cholangitis/psc_qsp_model.png)](../../../primary-sclerosing-cholangitis/psc_qsp_model.svg)

## Overview
Primary sclerosing cholangitis (PSC) is a chronic cholestatic liver disease characterised by progressive multifocal stricturing and fibrosis of the intrahepatic and extrahepatic bile ducts, with a prevalence of approximately 5~10 per 100,000 population. Inflammatory bowel disease (IBD, mainly ulcerative colitis) is present in 70~80% of cases, and median survival after diagnosis is 12~15 years. The core mechanism runs from gut dysbiosis-driven translocation of LPS via the gut-liver axis → bile duct epithelial injury and TH17-immune imbalance → LOXL2-mediated hepatic stellate cell activation → periductal fibrosis and stricturing. The annual risk of cholangiocarcinoma (CCA) is markedly elevated at 1.5~2%, making regular surveillance essential. There is currently no effective medical therapy; symptomatic management and, ultimately, liver transplantation remain the only curative option.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Gut-Liver Axis | Impaired gut barrier function → portal translocation of LPS → TLR4 activation | Triggers bile duct epithelial inflammation |
| Bile duct epithelial injury pathway | TH17/IL-17A, TNF-α, biliary epithelial senescence | Cholangiocyte loss |
| FXR-bile acid regulatory pathway | Reduced FXR activity → excess bile acid synthesis → toxic accumulation | Biliary toxicity, cholestasis |
| Hepatic stellate cell activation pathway | TGF-β/LOXL2 → HSC activation, collagen crosslinking | Periductal fibrosis and stricturing |
| Biliary fibrosis progression pathway | Col1a1 accumulation, increased Fibroscan stiffness | Portal hypertension, cirrhosis |
| Cholangiocarcinoma risk pathway | Chronic cholestasis, DNA damage, inflammation | CCA risk of 1.5~2% per year |

## Drug Targets
- **UDCA**: improves bile acid hydrophilicity and provides cytoprotection — has not been shown to improve survival in PSC, and may be harmful at high doses, so guideline opinions are divided
- **NorUDCA (24-norursodeoxycholic acid)**: a more hydrophilic bile acid substitute for the biliary tree, in PSC-specific phase 2/3 trials
- **Obeticholic acid (FXR agonist)**: FXR activation → suppresses bile acid synthesis and has antifibrotic activity — a PSC phase 3 trial is ongoing
- **Bezafibrate (PPARα)**: reduces bile acid toxicity and has anti-inflammatory activity — under clinical exploration in PSC
- **Antibiotics (metronidazole, vancomycin)**: modulate the gut microbiota — small studies ongoing in PSC
- **Liver transplantation**: for end-stage PSC and high cholangiocarcinoma risk — the curative treatment

## Model Files
| File | Description |
|------|------|
| [psc_qsp_model.dot](../../../primary-sclerosing-cholangitis/psc_qsp_model.dot) | Graphviz mechanistic map source (approximately 213 nodes / 11 clusters) |
| [psc_qsp_model.svg](../../../primary-sclerosing-cholangitis/psc_qsp_model.svg) | SVG vector image (scalable) |
| [psc_qsp_model.png](../../../primary-sclerosing-cholangitis/psc_qsp_model.png) | PNG image (150 dpi) |
| [psc_mrgsolve_model.R](../../../primary-sclerosing-cholangitis/psc_mrgsolve_model.R) | mrgsolve ODE model (approximately 27 compartments / 5 treatment scenarios) |
| [psc_shiny_app.R](../../../primary-sclerosing-cholangitis/psc_shiny_app.R) | Shiny dashboard |
| [psc_references.md](../../../primary-sclerosing-cholangitis/psc_references.md) | References (approximately 43 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: separate biliary/plasma PK compartments for UDCA, OCA, NorUDCA, and bezafibrate, plus PD compartments for LPS, gut barrier, FXR, bile acid pool, biliary hydrophilicity index, IL-17A, TNF-α, IL-6, Treg/IL-10, bile duct epithelial health, biliary senescence, hepatic stellate cells, Col1a1, LOXL2, ALP, bilirubin, Fibroscan, portal pressure, and cholangiocarcinoma risk
- **Key treatment scenarios**: ① untreated natural course, ② UDCA 15 mg/kg/day, ③ OCA 10 mg/day, ④ UDCA+OCA combination, ⑤ bezafibrate 400 mg/day
- **Calibration/evidence**: ALP and fibrosis progression rate parameters calibrated from the PSC natural history cohort (Boonstra 2013) and data from the PRIMROSE (NorUDCA) and AESOP (OCA for PSC) clinical trials

## Shiny Dashboard
Comprises tabs for patient profile (IBD status, bile duct stricture location, Amsterdam PSC risk score), bile acid/FXR PK/PD, biliary inflammation/fibrosis PD, liver function clinical endpoints (ALP, bilirubin, Fibroscan), treatment scenario comparison (ALP response, fibrosis suppression), and cholangiocarcinoma risk/biomarkers (LOXL2, IL-17A, Col1a1).

## Usage
```r
library(mrgsolve)
mod <- mread("psc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("psc_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg psc_qsp_model.dot -o psc_qsp_model.svg
```

## References
For detailed citations, see [psc_references.md](../../../primary-sclerosing-cholangitis/psc_references.md) (approximately 43 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
