# Hashimoto's Thyroiditis (HT) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Endocrine/Metabolic

[![HT QSP Model](ht_qsp_model.png)](ht_qsp_model.svg)

## Overview

Hashimoto's thyroiditis (chronic autoimmune thyroiditis) is the most common autoimmune thyroid disease, in which thyroid destruction by anti-TPO (thyroid peroxidase) and anti-Tg (thyroglobulin) autoantibodies and autoreactive T cells causes progressive hypofunction. Global prevalence is approximately 1~2%, and it occurs 7~10 times more frequently in women. Genetic factors such as HLA-DR3/DR5, CTLA-4, and PTPN22 act together with environmental factors including iodine excess, smoking, and stress. As autoimmune destruction progresses, T4 production falls and TSH rises, leading to overt hypothyroidism, for which levothyroxine (LT4) replacement is the standard treatment.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| T-cell autoimmunity | Th1 predominance → IFN-γ/TNF-α → thyrocyte apoptosis (Fas/FasL) | Thyrocyte destruction |
| B cells/autoantibodies | Anti-TPO IgG → complement activation/ADCC → further tissue damage | Elevated anti-TPO/Tg |
| HPT axis dysregulation | Thyroid destruction → T4↓ → increased pituitary TSH secretion → goitre | Elevated TSH, goitre |
| Antioxidant deficiency | Selenium deficiency → reduced GPx activity → ROS accumulation → thyrocyte damage | Persistent inflammation, sustained anti-TPO |
| Treg dysfunction | Reduced FOXP3+ Treg → disrupted Th17/Th1 balance | Sustained autoimmunity |
| Hypothyroidism | Loss of thyroid parenchyma → insufficient T4/T3 production | Fatigue, bradycardia, oedema, hyperlipidaemia |

## Drug Targets

- **Levothyroxine (LT4)**: T4 replacement → peripheral conversion to T3; normalisation of TSH is the treatment goal
- **Selenium**: supports selenoproteins (GPx, deiodinases) → reduces anti-TPO (evidence from the SELENOIT trial)
- **Liothyronine (LiT3)**: direct T3 replacement; combination therapy for patients with persistent symptoms on LT4 monotherapy
- **High-dose LT4 TSH-suppression therapy**: keeps TSH low with the aim of shrinking a goitre
- **Glucocorticoids**: used transiently in acute painful thyroiditis (subacute-onset Hashimoto's)
- **Methimazole (a thionamide)**: used short term during the initial hyperthyroid phase of Hashitoxicosis

## Model Files

| File | Description |
|------|------|
| [ht_qsp_model.dot](ht_qsp_model.dot) | Graphviz mechanistic map source (approximately 100+ nodes / 10 clusters) |
| [ht_qsp_model.svg](ht_qsp_model.svg) | SVG vector image (scalable) |
| [ht_qsp_model.png](ht_qsp_model.png) | PNG image (150 dpi) |
| [ht_mrgsolve_model.R](ht_mrgsolve_model.R) | mrgsolve ODE model (approximately 19 compartments / 7 treatment scenarios) |
| [ht_shiny_app.R](ht_shiny_app.R) | Shiny dashboard |
| [ht_references_en.md](ht_references_en.md) | References (approximately 44 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: includes the HPT axis (TRH, TSH), intrathyroidal T4/T3 synthesis, plasma T4/T3/fT4/fT3/rT3, tissue T4/T3, immune compartments (Th1, Treg, B cells, anti-TPO antibody, anti-Tg antibody, thyroid damage index), and drug PK (LT4 gut/central/peripheral, LiT3, selenium) compartments
- **Key treatment scenarios**: ① untreated baseline (progressive hypofunction), ② LT4 100 mcg/day, ③ selenium 200 mcg/day, ④ LT4 + selenium combination, ⑤ LT4 + LiT3 combination therapy, ⑥ high-dose LT4 (TSH suppression), ⑦ LT4 75 mcg + selenium 200 mcg
- **Calibration/evidence**: TSH/fT4 normalisation curves qualitatively calibrated from the SELENOIT trial (Ventura 2017), selenium data from Gartner et al. (JCEM 2002), and T4+T3 combination therapy data from Celi et al. (JCEM 2011)

## Shiny Dashboard

The dashboard comprises a patient profile tab (baseline TSH, anti-TPO titre, degree of thyroid damage, selenium deficiency status), thyroid hormone PK/PD dynamics, immune markers (anti-TPO/anti-Tg antibodies, T-cell balance), clinical endpoints (TSH normalisation, fT4 level), comparison of 7 treatment scenarios, and a biomarker tab (TSH, anti-TPO, selenium).

## Usage

```r
library(mrgsolve)
mod <- mread("ht_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("ht_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg ht_qsp_model.dot -o ht_qsp_model.svg
```

## References

For detailed citations, see [ht_references_en.md](ht_references_en.md) (approximately 44 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
