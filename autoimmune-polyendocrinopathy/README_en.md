# Autoimmune Polyendocrinopathy Syndrome (APS / APECED)
## Autoimmune Polyendocrinopathy Syndrome Type 1 — QSP Model

> **QSP Disease Model Library** · A Quantitative Systems Pharmacology (QSP) model automatically generated via Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Endocrine/Metabolic

[![APS QSP Model](aps_qsp_model.png)](aps_qsp_model.svg)

---

## Overview

Autoimmune Polyendocrinopathy-Candidiasis-Ectodermal Dystrophy (APECED) is a rare monogenic disease in which mutation of the **AIRE** (autoimmune regulator) gene destroys central immune tolerance, causing multiple endocrine organs to come under autoimmune attack. Also called APS Type 1, its three cardinal components are chronic mucocutaneous candidiasis (CMC), hypoparathyroidism, and Addison's disease.

- **Prevalence**: 1/25,000 in Finland, 1/14,400 in Sardinia, 1/9,000 among Iranian Jews
- **Age of onset**: mainly in childhood (mean age at first sign of CMC: 5 years)
- **Inheritance**: autosomal recessive (AIRE gene, 21q22.3), >300 mutation types identified
- **Diagnostic criteria**: at least 2 of the 3 major symptoms, or a family history of APECED

---

## Pathophysiology Summary

| Stage | Mechanism | Result |
|------|------|------|
| Loss of AIRE function | Inability of thymic medullary epithelial cells (mTEC) to express tissue-specific antigens (TSA) | Failure of negative selection of autoreactive T cells |
| Breakdown of central tolerance | Release of autoreactive CD4+/CD8+ T cells into the periphery | Multi-organ autoimmune attack |
| Autoantibody generation | Autoreactive T-cell help → B cells → IgG autoantibodies | Anti-21-OH (adrenal), anti-NALP5 (parathyroid), anti-GAD65 (pancreas), anti-TPO (thyroid) |
| Anti-IFN-α antibodies | APS1-specific pathology (neutralising anti-IFN-α/IFN-ω antibodies) | Increased CMC risk, vulnerability to viral infection |
| Organ destruction | Direct CTL damage + complement-mediated destruction | Loss of adrenal cortex, parathyroid, beta-cell, and thyroid function |

---

## Model Architecture

### Mechanistic Map

- **Clusters**: 11 subgraphs
- **Number of nodes**: 140+ (immune system, 4 target organs, drug PK/PD, clinical endpoints)
- **Key pathway**: AIRE → thymic negative selection → autoreactive T cells → autoantibodies → organ destruction → clinical phenotype

| Cluster | Content |
|---------|------|
| 1. AIRE & thymus | AIRE gene/protein, mTEC, TSA expression, negative selection, Treg generation |
| 2. Peripheral tolerance | Treg, CTLA-4, PD-1, regulatory B cells, anergy |
| 3. Autoantibodies | B-cell activation, GC response, plasma cells, 10 major autoantibodies |
| 4. Cytokines | IFN-γ, IL-12/23, IL-17, TNF-α, JAK-STAT, NF-κB |
| 5. Adrenal pathology | Adrenal cortex destruction, cortisol, ACTH, HPA axis, aldosterone |
| 6. Parathyroid pathology | PTG destruction, PTH, Ca²⁺, vitamin D activation, neuromuscular excitability |
| 7. Pancreatic pathology | Beta-cell destruction, insulin, blood glucose, HbA1c |
| 8. Thyroid pathology | Thyroid destruction, TSH, FT4, metabolic effects |
| 9. Hormone replacement PK | Hydrocortisone, fludrocortisone, calcitriol, levothyroxine |
| 10. Immunosuppressant PK | Cyclosporine A, abatacept, rituximab, tofacitinib, IVIG |
| 11. Clinical endpoints | Organ function scores, autoantibody titres, QoL, APS composite score |

---

### mrgsolve ODE Model (20 Compartments)

| # | Compartment | Unit | Normal value |
|---|------|------|--------|
| 1 | AIRE_func | 0-1 scale | 1.0 (normal) |
| 2 | AutoT_pool | cells/µL | 2.0 |
| 3 | Treg_pool | cells/µL | 15.0 |
| 4 | AutoAb_adren | U/mL (anti-21-OH) | 1.0 |
| 5 | AutoAb_PTG | U/mL (anti-NALP5) | 1.0 |
| 6 | AutoAb_beta | U/mL (anti-GAD65) | 1.0 |
| 7 | AutoAb_thy | U/mL (anti-TPO) | 1.0 |
| 8 | Adrenal_fn | % (0-100) | 100 |
| 9 | Cortisol_c | µg/dL | 12.0 |
| 10 | PTG_fn | % (0-100) | 100 |
| 11 | PTH_plasma | pg/mL | 40.0 |
| 12 | Ca_serum | mg/dL | 9.4 |
| 13 | Beta_mass | % (0-100) | 100 |
| 14 | Insulin_p | pmol/L | 60.0 |
| 15 | Glucose_p | mg/dL | 90.0 |
| 16 | Thyroid_fn | % (0-100) | 100 |
| 17 | TSH_plasma | mIU/L | 2.0 |
| 18 | FT4_plasma | ng/dL | 1.2 |
| 19-23 | Drug_CsA/Aba/RTX/JAKi/HC | ng/mL or µg/dL | 0 |

---

### Treatment Scenarios (7)

| # | Scenario | AIRE severity | Treatment |
|---|---------|------------|------|
| 1 | Natural course (severe) | 90% loss of function | None |
| 2 | HRT alone | 90% loss of function | HC 20mg/day |
| 3 | HRT + CsA | 90% loss of function | HC + cyclosporine A 3.5mg/kg/day |
| 4 | HRT + abatacept | 90% loss of function | HC + Abatacept 125mg SC/week |
| 5 | HRT + rituximab | 90% loss of function | HC + RTX 375mg/m² q6mo |
| 6 | HRT + tofacitinib | 90% loss of function | HC + JAKi 10mg/day |
| 7 | Early intervention (mild) | 30% loss of function | HC 15mg/day + early initiation |

---

## Key Drug PK/PD Parameters

| Drug | Mechanism | Key PK | Clinical effect |
|------|------|---------|---------|
| Hydrocortisone 20mg/day | GR agonist | F=0.96, t½=1.5h, CL=1.1L/min | Normalisation of AM cortisol (8-20 µg/dL) |
| Fludrocortisone 100µg/day | MR agonist | F=0.99, t½=3.5h | Normalisation of Na+/K+ balance |
| Calcitriol 0.5µg/day | VDR agonist | t½=5-8h | Normalisation of Ca²⁺ (8.5-10.5 mg/dL) |
| Levothyroxine 75-100µg/day | Thyroid hormone replacement | F=0.75-0.80, t½=9 days | Normalisation of TSH (0.4-4.0 mIU/L) |
| Cyclosporine A 3.5mg/kg/day | Calcineurin inhibition | F=0.30 (variable), t½=8-24h | ~50% AutoT suppression at Cp=150ng/mL |
| Abatacept 125mg SC/week | CD28/B7 blockade | F=0.79, t½=13 days | Improved Treg/AutoT ratio |
| Rituximab 375mg/m² q6mo | Anti-CD20 (B-cell depletion) | IV, t½=21 days | ~85% reduction in autoantibodies |
| Tofacitinib 10mg/day | JAK1/3 inhibition | F=0.74, t½=3h | 70% blockade of IFN-γ/IL-17 signalling |

---

## Clinical Endpoints

| Endpoint | Normal range | APS1 target |
|------|----------|----------|
| AM cortisol (08:00) | 8-20 µg/dL | ≥8 µg/dL on HRT |
| Serum Ca²⁺ (corrected) | 8.5-10.5 mg/dL | 8.0-9.0 mg/dL (safe range) |
| Serum PTH | 15-65 pg/mL | Unmeasurable (loss of function) |
| Fasting blood glucose | 70-100 mg/dL | <130 mg/dL (T1DM) |
| HbA1c | <5.7% | <7.0% (T1DM target) |
| TSH | 0.4-4.0 mIU/L | 0.5-2.5 mIU/L (on LT4 treatment) |
| Free T4 | 0.8-1.8 ng/dL | 1.0-1.5 ng/dL |
| Anti-21-OH Ab | <1 U/mL | Monitoring (>1 → adrenal risk) |
| Organ function score | 100% | Target: maintain >50% |

---

## Annual Screening Panel (APS1 Lifetime Monitoring)

Because new components can appear throughout the lifetime of APS1 patients, annual autoantibody screening is essential.

| Antibody | Target organ | Threshold | Action |
|------|---------|--------|------|
| Anti-21-hydroxylase IgG | Adrenal | >1 U/mL | AM cortisol, ACTH stimulation test |
| Anti-NALP5 IgG | Parathyroid | >10 U/mL | Serum Ca, PTH measurement |
| Anti-GAD65 IgG | Pancreatic β cells | >5 U/mL | Fasting glucose, OGTT |
| Anti-TPO/TG IgG | Thyroid | >34 IU/mL | TSH measurement |
| Anti-IFN-α neutralising antibody | Systemic antiviral defence | Positive | Confirmatory marker for APS1 |
| Anti-gastric parietal cell/intrinsic factor Ab | Gastrointestinal | Positive | Vitamin B12 level |
| Anti-17α-OH IgG | Gonads | Positive | LH/FSH/oestradiol |

---

## Files

| File | Description |
|------|------|
| [aps_qsp_model.dot](aps_qsp_model.dot) | Graphviz mechanistic map source (140+ nodes, 11 clusters) |
| [aps_qsp_model.svg](aps_qsp_model.svg) | SVG vector image (scalable, interactive) |
| [aps_qsp_model.png](aps_qsp_model.png) | PNG raster image (150 dpi) |
| [aps_mrgsolve_model.R](aps_mrgsolve_model.R) | mrgsolve ODE model (20 compartments, 7 scenarios) |
| [aps_shiny_app.R](aps_shiny_app.R) | Shiny dashboard (7 tabs: patient/immune/PK/endocrine/endpoints/scenarios/biomarkers) |
| [aps_references.md](aps_references.md) | 60 references (including PubMed links) |

---

## Key References

- Husebye ES et al. *N Engl J Med.* 2018;378:1132–1141
- Perheentupa J. *J Clin Endocrinol Metab.* 2006;91:2843–2850
- Anderson MS et al. *Science.* 2002;298:1395–1401
- Alimohammadi M et al. *N Engl J Med.* 2008;358:1018–1028
- Landegren N et al. *Sci Rep.* 2016;6:20104
