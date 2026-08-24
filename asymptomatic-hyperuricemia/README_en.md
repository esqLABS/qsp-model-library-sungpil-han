# Asymptomatic Hyperuricemia (AHY)
## QSP Disease Model — Quantitative Systems Pharmacology Model

---

## Overview

**Asymptomatic hyperuricemia (AHY)** refers to a state in which serum uric acid
concentration exceeds ≥7.0 mg/dL in men or ≥6.0 mg/dL in women, but without clinical
symptoms such as gout flares, tophi, or uric acid nephrolithiasis.

It shows a prevalence of 15-25% among adults worldwide, higher among Asians, and acts
over the long term as a risk factor for gout, chronic kidney disease (CKD),
hypertension, and cardiovascular disease.

---

## Mechanistic Map

[![AHY QSP Model](ahy_qsp_model_en.png)](ahy_qsp_model_en.svg)

> Click to open the full SVG vector image (100+ nodes, 15 subgraph clusters).

---

## Pathophysiology Summary

### Key uric acid metabolism pathways

| Pathway | Key enzyme/transporter | Clinical significance |
|------|-----------------|----------|
| Purine biosynthesis (de novo) | PRPP synthetase, enzyme complex | Overproduction → SUA↑ |
| Purine catabolism | Xanthine oxidase (XO/XDH) | Target of allopurinol and febuxostat |
| Renal reabsorption | URAT1 (SLC22A12), GLUT9 (SLC2A9) | Reduced uric acid excretion (underexcretor type) |
| Renal secretion | OAT1, OAT3, ABCG2, NPT1/4 | Reduced secretion → SUA↑ |
| Intestinal secretion | ABCG2 (intestinal) | Q141K variant → reduced secretion |
| MSU crystals | Saturation point 6.8 mg/dL | Dynamic balance of nucleation, growth, and dissolution |

### Long-term effects of hyperuricemia

```
SUA↑ (≥7.0mg/dL)
├─ Direct toxicity → endothelial damage → NO↓ → BP↑ → increased cardiovascular risk
├─ XO byproducts (H₂O₂, O₂⁻) → oxidative stress → increased atherosclerosis
├─ MSU crystal formation (>6.8mg/dL) → NLRP3 → IL-1β → gout flare
├─ Renal tubular damage → interstitial fibrosis → decreased eGFR → worsening CKD
└─ Increased insulin resistance → metabolic syndrome → NAFLD, T2DM
```

---

## Main Pathogenic Subtypes

| Type | Frequency | Cause | Characteristics |
|------|------|------|------|
| Underexcretor type | ~85% | URAT1↑, ABCG2↓, GFR↓ | Urinary UA <600mg/day |
| Overproducer type | ~10% | XO hyperactivity, HGPRT deficiency | Urinary UA >800mg/day |
| Mixed type | ~5% | Combined causes | Combined diet + genetics |

---

## Key Genetic Variants

| Gene | Variant | Effect | Frequency (Asian) |
|--------|------|------|------------|
| ABCG2 | Q141K (rs2231142) | 50%↓ intestinal + renal secretion | 34% |
| SLC22A12 | rs11602903 | Altered URAT1 reabsorption function | 5-20% |
| SLC2A9 | rs16890979 | Increased GLUT9 reabsorption | 10-25% |
| XDH | Multiple | Increased XO activity | Rare |
| HLA-B*58:01 | — | Allopurinol toxicity risk | 6-8% (Korea/China) |

---

## Drug PK/PD Parameters

### XO inhibitors

| Parameter | Allopurinol 300mg | Febuxostat 80mg |
|---------|----------------|-----------------|
| Bioavailability (F) | 90% | 84% |
| Peak concentration | 2-3 μg/mL | 3.4 μg/mL |
| Half-life (t½) | 1-2h (oxypurinol: 18-30h) | 5-8h |
| XO IC50 | Oxypurinol: 8.0mg/L | 0.001mg/L (non-purine) |
| SUA-lowering effect | 30-40% | 40-53% |
| Renal excretion | Main route | Dual route (hepatic/renal) |

### Uricosuric agents

| Parameter | Lesinurad 200mg | Dotinurad 4mg |
|---------|----------------|--------------|
| Bioavailability | ≈100% | ≈100% |
| t½ | 5h | 14-17h |
| Target transporter | URAT1+OAT4 | URAT1-selective |
| SUA reduction | Additional 15-30% with combination | 15-20% alone |

### Pegloticase (biologic agent)

| Parameter | Value |
|---------|---|
| Route of administration | IV 8mg q2weeks |
| t½ | 6-14 days (PEGylated) |
| Mechanism | UA → allantoin (10-fold↑ solubility) |
| SUA reduction | >80% (within 24 hours) |
| Failure rate | 40-50% (anti-PEG antibodies) |

---

## mrgsolve ODE Model Structure

### 19 Compartments

| # | Compartment | Description |
|---|------|------|
| 1 | UA_plasma | Serum uric acid (mg/dL) |
| 2 | UA_tissue | Tissue uric acid pool (mg) |
| 3 | XO_free | Xanthine oxidase activity (normalized) |
| 4 | Oxypurinol | Oxypurinol blood concentration (mg/L) |
| 5 | Febuxostat_C | Febuxostat blood concentration (mg/L) |
| 6 | Uricosuric_C | Uricosuric agent blood concentration (mg/L) |
| 7 | URAT1_free | Active URAT1 fraction |
| 8 | UrinaryUA | Urinary uric acid (mg/day) |
| 9 | MSU_depot | MSU crystal deposit amount (mg) |
| 10 | Endothelial_fn | Endothelial function (0-1) |
| 11 | NO_level | Nitric oxide level (rel.) |
| 12 | BP | Mean arterial pressure (mmHg) |
| 13 | GFR | Glomerular filtration rate (mL/min/1.73m²) |
| 14 | IL1beta | IL-1β (pg/mL) |
| 15 | CRP | hs-CRP (mg/L) |
| 16 | InsulinResist | HOMA-IR |
| 17 | CV_risk_score | Cumulative cardiovascular risk score |
| 18 | Tophus_vol | Tophus volume (mm³) |
| 19 | ABCG2_frac | Functional ABCG2 fraction |

---

## Treatment Scenarios

| Scenario | Description | Expected SUA change (2 years) |
|---------|------|-------------------|
| 1. Untreated AHY | SUA=7.5 baseline | →8.2mg/dL (gradual increase) |
| 2. Allopurinol 300mg | 30-40% XO inhibition | →5.2mg/dL (target achieved) |
| 3. Febuxostat 80mg | 40-53% XO inhibition | →4.5mg/dL (target achieved) |
| 4. Allopurinol + uricosuric combination | Dual mechanism | →4.0mg/dL (aggressive treatment) |
| 5. High-fructose diet + alcohol | Lifestyle risk | →9.5mg/dL (sharp rise) |
| 6. ABCG2 Q141K + febuxostat 120 | Genetic variant + maximum dose | →4.8mg/dL |
| 7. Aggressive target SUA<6 | Allopurinol 600mg + lifestyle correction | →5.0mg/dL |

---

## Shiny App Tab Structure

| Tab | Content |
|----|------|
| 1. Patient Profile | Disease overview, initial risk assessment, SUA time series |
| 2. Drug PK | Blood concentration-time curves, XO inhibition rate |
| 3. Uric Acid Dynamics | SUA, urinary UA, MSU crystal deposition, gout flare risk |
| 4. Cardiovascular/Renal Effects | eGFR, blood pressure, cardiovascular risk, endothelial function/NO |
| 5. Scenario Comparison | Simultaneous comparison of 7 treatment scenarios, 2-year outcome table |
| 6. Biomarkers | IL-1β, hs-CRP, HOMA-IR, SUA-CRP scatter plot |

---

## Model Files

| File | Description |
|------|------|
| [ahy_qsp_model_en.dot](ahy_qsp_model_en.dot) | Graphviz mechanistic map source (15 clusters, 100+ nodes) |
| [ahy_qsp_model_en.svg](ahy_qsp_model_en.svg) | SVG vector image (zoomable, high resolution) |
| [ahy_qsp_model_en.png](ahy_qsp_model_en.png) | PNG raster image (150 dpi) |
| [ahy_mrgsolve_model_en.R](ahy_mrgsolve_model_en.R) | mrgsolve ODE model (19 compartments, 7 scenarios, dose-response analysis) |
| [ahy_shiny_app_en.R](ahy_shiny_app_en.R) | Shiny dashboard (6-tab interactive simulator) |
| [ahy_references_en.md](ahy_references_en.md) | 45 references (PubMed links, classified by section) |

---

## Clinical Treatment Guidelines

| Society | Treatment recommendation (AHY) |
|------|----------------|
| ACR 2020 | SUA ≥9mg/dL + comorbidities → consider treatment (weak recommendation) |
| EULAR 2016 | ULT is generally not recommended in asymptomatic AHY |
| Chinese Rheumatology Association 2023 | SUA ≥8mg/dL + cardiovascular/renal risk → consider treatment |
| Japanese Society of Gout | Persistent SUA ≥8mg/dL → treatment recommended |

> **Key controversy**: Randomized controlled trial evidence on whether urate-lowering
> therapy improves cardiovascular and renal outcomes in AHY remains insufficient.

---

*Date created: 2026-06-20 | Claude Code Routine (CCR) | Disease category: Chronic disease/Metabolic*
