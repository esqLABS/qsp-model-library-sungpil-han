# Chronic Gastritis (CGAST) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology & Hepatobiliary

[![CGAST QSP Model](../../../chronic-gastritis/cgast_qsp_model.png)](../../../chronic-gastritis/cgast_qsp_model.svg)

## Overview

*Helicobacter pylori* infection is the most common cause of chronic
gastritis, estimated to infect about 44% of the world's population.
Following the Correa cascade, disease progresses from chronic gastritis
to atrophic gastritis, intestinal metaplasia, dysplasia, and finally
gastric cancer, with a cytokine cascade of IL-8, IL-1β, TNF-α, and others
via the NF-κB pathway driving mucosal damage. PPI-based triple and
quadruple eradication therapy is the key treatment strategy: it eradicates
H. pylori, resolves mucosal inflammation, and interrupts the Correa
cascade.

## Key Pathophysiological Pathways

| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| H. pylori infection and NF-κB activation | CagA/VacA toxins → NF-κB → IL-8 secretion | Neutrophil infiltration, mucosal damage |
| Th1/Treg imbalance | IFN-γ ↑, IL-10 reduced | Persistent chronic inflammation, loss of gastric glands |
| Abnormal acid regulation | Gastrin ↑ (body-predominant gastritis), increased acid secretion | Risk of peptic ulcer |
| Oxidative stress | ROS generation, damage to the mucus mucosal barrier | Cellular damage, promotes atrophy |
| Correa cascade — atrophy | Loss of parietal and chief cells, reduced PG I/II ratio | Reduced acid secretion, bacterial overgrowth |
| Progression of intestinal metaplasia | CDX2 expression, replacement by intestinal-type epithelium | Risk of malignant transformation |
| Recovery after eradication | Resolution of inflammation, partial reversal of atrophy (in early stages) | Reduced cancer risk |

## Key Drug Targets

- **PPIs (proton pump inhibitors)**: omeprazole, esomeprazole — inhibit H⁺/K⁺-ATPase, block acid secretion → support antibiotic efficacy
- **Vonoprazan (P-CAB)**: a potassium-competitive acid blocker, faster and more potent acid suppression than PPIs
- **Amoxicillin**: inhibits cell-wall synthesis, primary bactericidal action against H. pylori
- **Clarithromycin**: inhibits the 50S ribosomal subunit, the core of standard triple therapy
- **Metronidazole**: antianaerobic activity, an alternative for clarithromycin resistance
- **Bismuth (BSS)**: mucosal protection, anti-H. pylori activity → a component of quadruple therapy

## Model Files

| File | Description |
|------|------|
| [cgast_qsp_model.dot](../../../chronic-gastritis/cgast_qsp_model.dot) | Graphviz mechanistic map source (about 439 nodes / 10 clusters) |
| [cgast_qsp_model.svg](../../../chronic-gastritis/cgast_qsp_model.svg) | SVG vector image (zoomable) |
| [cgast_qsp_model.png](../../../chronic-gastritis/cgast_qsp_model.png) | PNG image (150 dpi) |
| [cgast_mrgsolve_model.R](../../../chronic-gastritis/cgast_mrgsolve_model.R) | mrgsolve ODE model (about 22 compartments / 7 treatment scenarios) |
| [cgast_shiny_app.R](../../../chronic-gastritis/cgast_shiny_app.R) | Shiny dashboard |
| [cgast_references.md](cgast_references.md) | References (about 61, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: PPI, amoxicillin, clarithromycin, and metronidazole PK compartments plus H. pylori burden, NF-κB, IL-8, IL-1β, TNF-α, IFN-γ, IL-10, neutrophils, Th1, Treg, gastrin, gastric acid, mucosa, an atrophy score, an intestinal-metaplasia score, and symptom PD compartments
- **Key treatment scenarios**: (1) untreated natural history, (2) PPI alone, (3) standard triple therapy (PPI+AMX+CLR), (4) bismuth quadruple therapy, (5) metronidazole quadruple therapy, (6) vonoprazan-based triple therapy, (7) 5-year follow-up after eradication
- **Calibration/basis**: eradication rates for each scenario are based on data from major guidelines such as Malfertheiner et al. (Gut 2017)

## Shiny Dashboard

Structured into 6 or more tabs: (1) patient profile — sets infection
stage and antibiotic resistance; (2) PK tab — PPI and antibiotic plasma
concentrations; (3) key PD measures — H. pylori burden, cytokines
(IL-8, IL-1β, TNF-α), gastrin and gastric acid; (4) clinical endpoints —
atrophy/intestinal-metaplasia score, symptom score; (5) scenario
comparison — outcomes across the 7 eradication regimens; (6) biomarkers —
PGI/PGII ratio (an atrophy marker), CagA status

## Usage

```r
library(mrgsolve)
mod <- mread("cgast_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("cgast_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg cgast_qsp_model.dot -o cgast_qsp_model.svg
```

## References

For full citations, see [cgast_references.md](cgast_references.md) (about 61 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
