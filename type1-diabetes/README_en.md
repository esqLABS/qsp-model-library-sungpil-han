# Type 1 Diabetes (T1DM) — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Endocrine/Metabolic

[![T1DM QSP Model](t1dm_qsp_model.png)](t1dm_qsp_model.svg)

## Overview

Type 1 diabetes (T1DM) is an autoimmune disease in which autoreactive
CD8+ T cells (CTLs) selectively destroy pancreatic beta cells, causing
absolute insulin deficiency. Worldwide prevalence is about 8.5 million
(2022), with about 180,000 children and adolescents newly diagnosed each
year. Beta-cell loss progresses slowly over years before diagnosis, and
by the time of clinical presentation, 70-90% of beta-cell mass has
already been lost. Intensive insulin therapy is the foundation of
glycaemic control, and the CD3 antibody teplizumab is the first
immunomodulatory therapy approved for delaying disease onset.

## Key Pathophysiological Pathways

| Pathway | Key molecules/mechanism | Clinical result |
|------|----------------|-----------|
| Autoimmune beta-cell attack | GAD, IA-2, ZnT8 autoantibodies; MHC-I-mediated beta-cell lysis by CD8+ CTLs | Reduced beta-cell mass, loss of C-peptide |
| Treg dysfunction | Reduced Foxp3+ Treg activity, breakdown of immune tolerance | Failure to suppress CTLs, persistent inflammation |
| Insulin deficiency | Beta-cell loss → cessation of basal and postprandial insulin secretion | Loss of glycaemic control, risk of ketoacidosis |
| Glucose-insulin kinetics | Dysregulated hepatic glucose production (EGP), reduced peripheral glucose utilisation | Fasting and postprandial hyperglycaemia |
| Impaired glucagon counter-regulation | Absence of insulin suppression → glucagon oversecretion, increased EGP | Poor recovery after hypoglycaemia |
| HbA1c kinetics | Mean glucose → accumulation of glycated haemoglobin | Risk of chronic complications (retinopathy, nephropathy, neuropathy, vascular disease) |

## Key Drug Targets

- **Basal insulin** — insulin degludec, glargine: fasting glucose control, prevention of nocturnal hypoglycaemia
- **Mealtime (bolus) insulin** — rapid-acting insulin aspart, lispro: suppresses postprandial hyperglycaemia
- **Insulin pump (CSII) + hybrid closed loop (HCL)** — automated insulin delivery linked to continuous glucose monitoring (CGM)
- **CD3 antibody — teplizumab**: reprogrammes autoreactive effector T cells, delays disease onset (TrialNet TN10)
- **GAD vaccine / anti-CD20 (rituximab)** — investigational: blocks the beta-cell autoimmune pathway

## Model Files

| File | Description |
|------|------|
| [t1dm_qsp_model.dot](t1dm_qsp_model.dot) | Graphviz mechanistic map source (about 236 nodes / 12 clusters) |
| [t1dm_qsp_model.svg](t1dm_qsp_model.svg) | SVG vector image (zoomable) |
| [t1dm_qsp_model.png](t1dm_qsp_model.png) | PNG image (150 dpi) |
| [t1dm_mrgsolve_model.R](t1dm_mrgsolve_model.R) | mrgsolve ODE model (about 20 compartments / 6 treatment scenarios) |
| [t1dm_shiny_app.R](t1dm_shiny_app.R) | Shiny dashboard |
| [t1dm_references.md](t1dm_references.md) | References (about 52, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: 3 compartments for beta cells/CTL/Treg (immune), 1 compartment for C-peptide, 2 compartments for glucose kinetics (Gp, Gt), 5 compartments for insulin PK (SC1, SC2, central, peripheral, remote action), 1 compartment for glucagon, 1 for HbA1c, 1 for teplizumab PK, 2 for CGM/APC, 3 for meal glucose (Qsto1, Qsto2, Qgut), and 1 for mean glucose — about 20 compartments in total
- **Key treatment scenarios**: (1) untreated T1DM, (2) MDI (basal + bolus insulin), (3) CSII (insulin pump), (4) hybrid closed loop (HCL), (5) teplizumab (stage 2 onset prevention), (6) teplizumab + MDI combination
- **Calibration/basis**: glucose-insulin kinetics based on the Bergman minimal model; the 2-year delay effect from the TrialNet TN10 trial (teplizumab); the Dalla Man meal-absorption model; and the HbA1c-complication relationship from the DCCT study

## Shiny Dashboard

Structured into 6 or more tabs: a patient profile tab (sets body weight,
age, residual beta-cell function, baseline HbA1c), PK visualisation
(insulin and teplizumab blood concentration-time curves), glycaemic PD
measures (fasting/postprandial glucose, simulated CGM), clinical
endpoints (HbA1c, hypoglycaemia frequency, preservation of C-peptide),
treatment scenario comparison (beta-cell mass and HbA1c trajectory), and a
biomarker panel (CTL, Treg, C-peptide, glucagon).

## Usage

```r
library(mrgsolve)
mod <- mread("t1dm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("t1dm_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg t1dm_qsp_model.dot -o t1dm_qsp_model.svg
```

## References

For full citations, see [t1dm_references.md](t1dm_references.md) (about 52 entries).

---
*This model is a qualitative/semi-quantitative QSP model for educational and research purposes and must not be used directly for clinical decision-making.*
