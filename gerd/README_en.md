# Gastroesophageal Reflux Disease (GERD) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Gastroenterology/Hepatobiliary

[![GERD QSP Model](gerd_qsp_model.png)](gerd_qsp_model.svg)

## Overview

Gastroesophageal reflux disease (GERD) is a condition in which gastric acid and contents reflux into the oesophagus, causing mucosal damage and symptoms (heartburn, regurgitation), with a prevalence of approximately 10~20% among Western adults. The core pathophysiology comprises transient lower oesophageal sphincter relaxations (TLESRs) or reduced LES function, decreased oesophageal acid clearance, and impaired oesophageal mucosal defence mechanisms. The disease spans a broad spectrum from erosive reflux disease (ERD) to Barrett's oesophagus, and PPIs (proton pump inhibitors) and P-CABs (potassium-competitive acid blockers), which potently suppress gastric acid secretion, are the principal therapies.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| LES dysfunction | Increased TLESR frequency, reduced basal LES pressure | Increased frequency of acid reflux |
| Gastric acid secretion | H+/K+-ATPase → maintains gastric pH 1~2 | Oesophageal damage upon reflux |
| Oesophageal acid exposure | Refluxed acid → oesophageal mucosal epithelial injury | Elevated DeMeester score |
| Mucosal defence mechanisms | Reduced mucus, bicarbonate, and epithelial proliferation | Erosion and ulcer formation |
| Inflammatory cascade | Acid/pepsin → IL-8, TNF-α → neutrophil infiltration | Mucosal inflammation, erosion |
| Barrett's oesophagus progression | Chronic acid exposure → intestinal metaplasia | Increased risk of oesophageal adenocarcinoma |
| Oesophageal motility disorder | Weakened oesophageal peristalsis → delayed acid clearance | Increased nocturnal reflux |

## Drug Targets

- **PPIs** (omeprazole, esomeprazole, pantoprazole): irreversibly inhibit H+/K+-ATPase → maintain gastric pH ≥ 4, the standard treatment
- **P-CABs** (vonoprazan, tegoprazan): potassium-competitive reversible H+/K+-ATPase inhibition → faster onset, superior nocturnal acid suppression
- **H2 receptor antagonists** (famotidine, ranitidine): block histamine H2R → adjunctive acid suppression
- **Prokinetics** (metoclopramide, dopamine antagonists): promote gastric emptying, raise LES pressure
- **Alginates/antacids**: physical acid neutralisation and a barrier against acid reflux

## Model Files

| File | Description |
|------|------|
| [gerd_qsp_model.dot](gerd_qsp_model.dot) | Graphviz mechanistic map source (approximately 176 nodes / 11 clusters) |
| [gerd_qsp_model.svg](gerd_qsp_model.svg) | SVG vector image (scalable) |
| [gerd_qsp_model.png](gerd_qsp_model.png) | PNG image (150 dpi) |
| [gerd_mrgsolve_model.R](gerd_mrgsolve_model.R) | mrgsolve ODE model (approximately 23 compartments / 6 treatment scenarios) |
| [gerd_shiny_app.R](gerd_shiny_app.R) | Shiny dashboard |
| [gerd_references.md](gerd_references.md) | References (approximately 41 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK compartments (separate GI/central compartments for PPI, H2RA, and P-CAB) + PD compartments (H+/K+-ATPase activity, intragastric pH, gastric acid secretion rate, LES pressure/TLESR frequency, oesophageal acid exposure time, mucosal damage index, Barrett's progression risk, inflammation index)
- **Key treatment scenarios**: ① no treatment, ② standard-dose PPI (omeprazole 20 mg QD), ③ high-dose PPI (esomeprazole 40 mg QD), ④ P-CAB (vonoprazan 20 mg QD), ⑤ H2RA (famotidine 20 mg BID), ⑥ PPI + H2RA combination (suppressing nocturnal breakthrough symptoms)
- **Calibration/evidence**: parameters referenced from the Metz PPI PD model, the vonoprazan PHALCON-EE trial, and normal values from oesophageal pH-impedance monitoring

## Shiny Dashboard

Comprises 6 tabs: ① Patient Profile (sets baseline acid secretion, LES function, erosion grade, H. pylori status), ② PK tab (PPI/P-CAB/H2RA blood concentration, CYP2C19 metaboliser type), ③ Gastric acid/oesophageal PD tab (intragastric pH, acid exposure time, TLESR frequency trends), ④ Clinical Endpoints (symptom remission rate, erosion healing rate, Barrett's risk), ⑤ Scenario Comparison (simultaneous comparison of 6 treatment strategies), ⑥ Biomarkers (H+/K+-ATPase activity, mucosal damage index, inflammation trends).

## Usage

```r
library(mrgsolve)
mod <- mread("gerd_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("gerd_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg gerd_qsp_model.dot -o gerd_qsp_model.svg
```

## References

For detailed citations, see [gerd_references.md](gerd_references.md) (approximately 41 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
