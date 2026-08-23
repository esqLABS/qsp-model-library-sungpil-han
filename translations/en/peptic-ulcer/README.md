# Peptic Ulcer Disease (PUD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology/Hepatobiliary

[![PUD QSP Model](../../../peptic-ulcer/pud_qsp_model.png)](../../../peptic-ulcer/pud_qsp_model.svg)

## Overview
Peptic ulcer disease is a chronic ulcer occurring in the gastric or duodenal mucosa, with a worldwide prevalence of approximately 10%, causing substantial hospitalisation and bleeding complications each year. The core pathogenic mechanism is an imbalance between mucosal defence and attack driven by *Helicobacter pylori* infection and NSAID use: H. pylori's CagA and VacA toxins cause mucosal damage, while NSAIDs reduce prostaglandin production via COX-1 inhibition. The principal drug targets are the gastric acid secretion pump (H⁺/K⁺-ATPase), the H. pylori organism, and mucosal protective mechanisms.

## Key Pathways
| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| H. pylori toxin pathway | CagA T4SS, VacA, oxidative stress | Mucosal epithelial injury, inflammation |
| Gastric acid secretion pathway | H⁺/K⁺-ATPase (proton pump), gastrin/histamine | Raised intragastric acidity |
| Mucosal defence pathway | Mucus, prostaglandin E2, bicarbonate | Maintenance of the mucosal protective barrier |
| NSAID pathway | COX-1 inhibition → reduced PGE2 | Reduced mucosal blood flow, ulcer induction |
| Inflammatory amplification pathway | IL-1β, TNF-α, NF-κB activation | Neutrophil infiltration, tissue damage |
| Ulcer healing pathway | EGF, TGF-β, cell regeneration | Reduced ulcer area |

## Drug Targets
- **Proton pump inhibitors (PPIs)**: omeprazole, esomeprazole — irreversibly inhibit H⁺/K⁺-ATPase, block gastric acid secretion
- **H2 receptor antagonists (H2RAs)**: ranitidine, famotidine — block the histamine H2 receptor, suppress nocturnal acid secretion
- **H. pylori eradication antibiotics**: amoxicillin (AMX), clarithromycin (CLR) — triple/quadruple therapy
- **Cytoprotective agent**: misoprostol — a PGE1 analogue providing mucosal protection
- **Combined antibacterial regimen**: metronidazole, bismuth — quadruple therapy for resistant strains

## Model Files
| File | Description |
|------|------|
| [pud_qsp_model.dot](../../../peptic-ulcer/pud_qsp_model.dot) | Graphviz mechanistic map source (approximately 174 nodes / 9 clusters) |
| [pud_qsp_model.svg](../../../peptic-ulcer/pud_qsp_model.svg) | SVG vector image (scalable) |
| [pud_qsp_model.png](../../../peptic-ulcer/pud_qsp_model.png) | PNG image (150 dpi) |
| [pud_mrgsolve_model.R](pud_mrgsolve_model.R) | mrgsolve ODE model (approximately 18 compartments / 5 treatment scenarios) |
| [pud_shiny_app.R](../../../peptic-ulcer/pud_shiny_app.R) | Shiny dashboard |
| [pud_references.md](../../../peptic-ulcer/pud_references.md) | References (approximately 55 articles, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: separate GI/central/peripheral PK compartments for PPI, H2RA, AMX, CLR, and NSAID, plus H. pylori burden (log10), fractional proton pump activity, intragastric pH, mucosal layer (mucin), prostaglandin, inflammation score, and ulcer area
- **Key treatment scenarios**: ① PPI BID monotherapy, ② H2RA BID monotherapy, ③ triple therapy (PPI+AMX+CLR), ④ NSAID monotherapy, ⑤ NSAID + PPI combination
- **Calibration/evidence**: parameters calibrated with reference to eradication rates and ulcer healing rates from the Maastricht V guidelines and major H. pylori eradication trials (MLST data)

## Shiny Dashboard
An interactive dashboard comprising tabs for patient profile (H. pylori status, NSAID use, risk factors), drug PK simulation, intragastric pH and degree of acid suppression, ulcer area change and healing rate, treatment scenario comparison (eradication rate, recurrence rate), and key biomarkers (CRP, PG, HP burden).

## Usage
```r
library(mrgsolve)
mod <- mread("pud_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("pud_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg pud_qsp_model.dot -o pud_qsp_model.svg
```

## References
For detailed citations, see [pud_references.md](../../../peptic-ulcer/pud_references.md) (approximately 55 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
