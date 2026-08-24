# Type 2 Diabetes (T2DM) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Endocrine/Metabolic

[![T2DM QSP Model](t2dm_qsp_model.png)](t2dm_qsp_model.svg)

## Overview

Type 2 diabetes (T2DM) is the most common metabolic disease, arising from the combined effects of peripheral insulin resistance and pancreatic beta-cell dysfunction, and affects approximately 530 million people worldwide (IDF 2021). Blood glucose is initially maintained by compensatory insulin hypersecretion, but as beta-cell exhaustion progresses this eventually leads to overt hyperglycaemia. Excess free fatty acids (FFA) and ectopic fat accumulation are the key mechanisms of hepatic and peripheral insulin resistance, and reduced GLP-1 secretion together with DPP-4 overactivity further worsens glycaemic control. Diverse drug classes including metformin, GLP-1 receptor agonists, SGLT-2 inhibitors, and insulin are used through complementary mechanisms.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Peripheral insulin resistance | IRS-1 serine phosphorylation, reduced GLUT4 translocation, FFA/ceramide toxicity | Reduced skeletal muscle glucose uptake |
| Hepatic insulin resistance | Hepatic fat accumulation → suppressed insulin signalling → increased EGP | Fasting hyperglycaemia, excess hepatic glucose production |
| Beta-cell dysfunction | Glucotoxicity, lipotoxicity, oxidative stress, amyloid deposition → beta-cell loss | Progressive decline in insulin secretion |
| Reduced GLP-1/incretin effect | Reduced intestinal L-cell GLP-1 secretion, increased DPP-4 degradation | Reduced postprandial insulin response, failure of glucagon suppression |
| Increased renal glucose reabsorption | SGLT2 overexpression, raised renal glucose threshold | Reduced urinary glucose loss → contributes to hyperglycaemia |
| Vascular/renal damage complications | Chronic hyperglycaemia → AGE accumulation, reduced glomerular filtration rate, albuminuria | Decreased eGFR, increased UACR, cardiovascular risk |

## Drug Targets

- **Metformin** — inhibits hepatic EGP (via AMPK activation), weight neutral: first-line therapy
- **GLP-1 receptor agonists** — semaglutide: promotes insulin secretion, suppresses glucagon, reduces weight, cardiovascular protection (SUSTAIN, PIONEER)
- **SGLT-2 inhibitors** — empagliflozin: blocks renal glucose reabsorption, protects against heart failure and renal disease (EMPA-REG OUTCOME)
- **DPP-4 inhibitors** — sitagliptin: inhibits GLP-1 degradation, lowers glucose, weight neutral
- **Sulfonylureas** — glimepiride: promotes beta-cell insulin secretion (hypoglycaemia risk)
- **Insulin** — degludec: basal glycaemic control, essential in advanced beta-cell failure

## Model Files

| File | Description |
|------|------|
| [t2dm_qsp_model.dot](t2dm_qsp_model.dot) | Graphviz mechanistic map source (approximately 331 nodes / 11 clusters) |
| [t2dm_qsp_model.svg](t2dm_qsp_model.svg) | SVG vector image (scalable) |
| [t2dm_qsp_model.png](t2dm_qsp_model.png) | PNG image (150 dpi) |
| [t2dm_mrgsolve_model.R](t2dm_mrgsolve_model.R) | mrgsolve ODE model (approximately 27 compartments / 7 treatment scenarios) |
| [t2dm_shiny_app.R](t2dm_shiny_app.R) | Shiny dashboard |
| [t2dm_references.md](t2dm_references.md) | References (approximately 40 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK (metformin 3 compartments, empagliflozin 1 compartment, semaglutide SC 2 compartments, DPP-4 inhibitor 1 compartment, sulfonylurea 1 compartment, insulin SC 2 compartments, pioglitazone 1 compartment) + glucose-insulin dynamics (Gp, Gt, Ip, X-action) + endocrine (glucagon, GLP-1, beta-cell mass) + insulin resistance (IR_H, IR_P, FFA) + body weight + clinical endpoints (HbA1c, eGFR, UACR), for a total of approximately 27 compartments
- **Key treatment scenarios**: (1) no treatment, (2) metformin monotherapy, (3) metformin + empagliflozin, (4) metformin + semaglutide, (5) DPP-4 inhibitor, (6) insulin, (7) triple combination (Met + Empa + Sema)
- **Calibration/evidence**: parameters referenced from EMPA-REG OUTCOME (cardiorenal protection with empagliflozin), SUSTAIN-6 (MACE reduction with semaglutide), UKPDS (long-term benefits of metformin), and the UKPDS beta-cell function decline model

## Shiny Dashboard

The dashboard comprises six or more tabs: patient profile (setting body weight, baseline HbA1c, eGFR, and insulin resistance level), PK visualisation (blood concentration-time curves for each drug), glycaemic PD metrics (fasting glucose and HbA1c time series), clinical endpoints (weight change, eGFR, UACR), treatment scenario comparison (HbA1c trajectories for 7 regimens), and a biomarker panel (GLP-1, glucagon, beta-cell mass).

## Usage

```r
library(mrgsolve)
mod <- mread("t2dm_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("t2dm_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg t2dm_qsp_model.dot -o t2dm_qsp_model.svg
```

## References

For detailed citations, see [t2dm_references.md](t2dm_references.md) (approximately 40 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
