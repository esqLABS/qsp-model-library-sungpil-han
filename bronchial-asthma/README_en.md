# Bronchial Asthma (BA) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Respiratory

[![BA QSP Model](ba_qsp_model.png)](ba_qsp_model.svg)

## Overview
Bronchial asthma is the most common chronic airway disease, affecting
about 300 million people worldwide, characterised by reversible airway
obstruction, airway hyperresponsiveness, and chronic airway inflammation
triggered by exposure to allergens, viruses, and pollutants. Th2/
eosinophil-centred type 2 inflammation (IL-4, IL-5, IL-13, TSLP) is at the
core of allergic asthma, and inhaled corticosteroids (ICS) plus LABA are
the first-line standard of care. For severe eosinophilic asthma,
biologics targeting IL-5 (mepolizumab, benralizumab), IL-4Rα (dupilumab),
and TSLP (tezepelumab) are used.

## Key Pathophysiological Pathways
| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| TSLP-Th2 activation | Epithelial cells → TSLP → DC → Th2 differentiation | IgE production, onset of allergic inflammation |
| IL-5-mediated eosinophil production/activation | IL-5 → bone marrow eosinophil differentiation and migration into blood | Blood and tissue eosinophilia |
| IL-13-mediated airway remodelling | IL-13 → goblet cell hyperplasia, mucus hypersecretion, AHR | FEV1 decline, chronic airway obstruction |
| IgE-mast cell axis | FcεRI → mast cell degranulation → histamine, leukotrienes | Early allergic reaction, bronchospasm |
| Airway smooth muscle contraction | β2-AR signalling (cAMP), leukotriene receptors | Reversible bronchoconstriction |
| Airway remodelling | TGF-β, basement-membrane thickening, smooth-muscle hypertrophy | Irreversible airway obstruction |

## Key Drug Targets
- **ICS (inhaled corticosteroids)**: GR-mediated suppression of airway inflammation broadly — basic therapy at every step of asthma
- **LABA (long-acting β2 agonists)**: β2-AR → cAMP → airway smooth-muscle relaxation — combined with ICS
- **Mepolizumab (anti-IL-5)**: neutralises IL-5 → reduces eosinophils — severe eosinophilic asthma (MENSA/SIRIUS)
- **Benralizumab (anti-IL-5Rα)**: blocks IL-5Rα → depletes eosinophils via ADCC — dosed q4w then q8w
- **Dupilumab (anti-IL-4Rα)**: dual blockade of IL-4/IL-13 — type 2-high asthma (QUEST/VENTURE)
- **Tezepelumab (anti-TSLP)**: blocks TSLP upstream → suppresses the entire type 2 cascade (NAVIGATOR)
- **Omalizumab (anti-IgE)**: neutralises free IgE, downregulates FcεRI — severe allergic asthma

## Model Files
| File | Description |
|------|------|
| [ba_qsp_model.dot](ba_qsp_model.dot) | Graphviz mechanistic map source (about 215 nodes / 16 clusters) |
| [ba_qsp_model.svg](ba_qsp_model.svg) | SVG vector image (zoomable) |
| [ba_qsp_model.png](ba_qsp_model.png) | PNG image (150 dpi) |
| [ba_mrgsolve_model_en.R](ba_mrgsolve_model_en.R) | mrgsolve ODE model (about 28 compartments / 5 treatment scenarios) |
| [ba_shiny_app.R](ba_shiny_app.R) | Shiny dashboard |
| [ba_references_en.md](ba_references_en.md) | References (about 41, PubMed links) |

## mrgsolve Model (ODE Model)
- **Compartment structure**: biologic PK (2-compartment SC PK each for mepolizumab, benralizumab, dupilumab, tezepelumab, plus TMDD for omalizumab), ICS/LABA PK (pulmonary and systemic), and immune PD (TSLP, IL-5, IL-13, blood eosinophils, tissue eosinophils, ASM tone, mucus, FEV1)
- **Key treatment scenarios**: (1) ICS/LABA alone, (2) + mepolizumab, (3) + benralizumab, (4) + dupilumab, (5) + tezepelumab
- **Calibration/basis**: referenced against FEV1 and eosinophil-reduction data from the MENSA (mepolizumab), SIROCCO (benralizumab), QUEST (dupilumab), and NAVIGATOR (tezepelumab) clinical trials

## Shiny Dashboard
Structured into tabs for patient profile (eosinophil count, IgE, asthma
severity, subtype selection), pharmacokinetics (plasma concentration for
each biologic), airway inflammation PD (TSLP, IL-5, IL-13, eosinophils),
clinical endpoints (FEV1, exacerbation rate, symptom score), treatment
scenario comparison (overlay of the 5 biologics), and biomarker prediction
(drug response by biomarker).

## Usage
```r
library(mrgsolve)
mod <- mread("ba_mrgsolve_model_en.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("ba_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg ba_qsp_model.dot -o ba_qsp_model.svg
```

## References
For full citations, see [ba_references_en.md](ba_references_en.md) (about 41 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
