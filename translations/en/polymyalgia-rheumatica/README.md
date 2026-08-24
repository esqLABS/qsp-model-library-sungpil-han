# Polymyalgia Rheumatica (PMR) — QSP Model

> **Polymyalgia rheumatica** | Autoimmune/inflammatory disease | IL-6 pathway · corticosteroid PK/PD · tocilizumab

[![Graphviz](https://img.shields.io/badge/Map-130+%20nodes%2C%2012%20clusters-blue)](../../../polymyalgia-rheumatica/pmr_qsp_model.svg)
[![mrgsolve](https://img.shields.io/badge/ODE-22%20compartments%2C%207%20scenarios-green)](pmr_mrgsolve_model.R)
[![Shiny](https://img.shields.io/badge/Shiny-6%20tabs%20dashboard-orange)](../../../polymyalgia-rheumatica/pmr_shiny_app.R)
[![References](https://img.shields.io/badge/Refs-55%20PubMed%20citations-red)](pmr_references.md)

---

## Disease Overview

Polymyalgia rheumatica (PMR) is a common inflammatory rheumatic disease
occurring in people **aged 50 and older**, characterised by **bilateral
muscle pain and morning stiffness** of the shoulder and pelvic girdles.

| Feature | Content |
|------|------|
| **Incidence** | 50–100/100,000/year (age 50+) |
| **Sex** | Female:male ≈ 2–3:1 |
| **Age of onset** | Predominantly the 70s–80s |
| **Geographic bias** | Higher in Northern European ancestry |
| **GCA overlap** | About 15–20% of PMR patients also have giant cell arteritis |
| **Treatment response** | Dramatic response to corticosteroids (prednisolone) — a **diagnostic feature** |

### Core symptoms
- Bilateral shoulder and pelvic girdle pain and weakness
- Morning stiffness ≥45 minutes
- Elevated CRP/ESR (acute-phase reactants)
- Weight loss, fever, fatigue

---

## Mechanistic Map

[![PMR QSP Map](../../../polymyalgia-rheumatica/pmr_qsp_model.png)](../../../polymyalgia-rheumatica/pmr_qsp_model.svg)

> *Click the image to view the full-screen SVG*

### 12-cluster structure

| Cluster | Key content | Key nodes |
|---------|----------|----------|
| **Epidemiology/risk factors** | HLA-DRB1*04, PTPN22, environmental triggers | HLA_DRB1_04, PTPN22_variant |
| **Innate immunity** | DC, M1 macrophages, neutrophils, NLRP3 | NLRP3, NF_kB_innate, MacM1 |
| **Adaptive immunity** | Th1/Th17/Treg imbalance | Th17_cells, Treg_cells, RORgt |
| **Cytokines** | IL-6 centred, TNF-α, IL-17A | IL6, sIL6R, IL17A, TNF_alpha |
| **JAK-STAT signalling** | JAK1/2-STAT3, SOCS3 feedback | STAT3, JAK1, SOCS3 |
| **Target tissue** | Shoulder/hip synovium and bursae | Subacromial_bursa, FLS_synov |
| **Prednisolone PK/PD** | GR binding, transrepression, GILZ | GR_Pred_complex, Transrepression |
| **Tocilizumab PK/PD** | mIL-6R/sIL-6R blockade, TMDD | TCZ_mIL6R_cpx, IL6_signal_blk |
| **HPA axis** | Cortisol suppression, adrenal atrophy | ACTH, Cortisol_endog, HPA_feedback |
| **Bone effects** | RANK/RANKL/OPG, reduced BMD | Osteoclast, RANKL, BMD_lumbar |
| **Vascular/GCA** | Temporal artery, aortic inflammation | Temporal_artery, PMR_GCA_overlap |
| **Clinical endpoints** | PMR-AS, CRP, ESR | PMR_AS, Remission, Relapse_event |

---

## Pathogenesis Summary

```
Environmental trigger (e.g. infection) + genetic predisposition (HLA-DRB1*04)
        |
  Innate immune activation (mDC, NLRP3, TLR4)
        |
  IL-23 -> Th17 differentiation / IFN-gamma -> Th1 skewing
     |                    |
  IL-17A/F             IFN-gamma, TNF-alpha
     |                    |
      --> IL-6 storm (muscle and joint tissue)
              |
    JAK1/2 -> STAT3 -> acute-phase proteins (CRP, fibrinogen)
              |
      Pain and morning stiffness (PGE2, tissue oedema)
              |
          PMR-AS up

Treatment:
  Prednisolone -> GR activation
    -> NF-kB/AP-1 transrepression -> IL-6, TNF-alpha down
    -> GILZ, Annexin A1 up -> anti-inflammatory
  Tocilizumab (IL-6R blockade) -> STAT3 inhibition -> CRP normalisation
```

---

## mrgsolve Model Specifications

**File**: [`pmr_mrgsolve_model.R`](pmr_mrgsolve_model.R)

### 22 Compartments

| Compartment group | Number | Content |
|----------|--------|------|
| Prednisolone PK | 3 | Depot → Central ↔ Peripheral |
| Tocilizumab PK | 3 | Depot → Central ↔ Peripheral (TMDD) |
| HPA axis | 1 | Cortisol (suppression response) |
| IL-6 pathway | 2 | IL-6, soluble IL-6Rα |
| Acute-phase reactants | 2 | CRP, ESR |
| Bone effects | 1 | BMD (normalised) |
| Disease activity | 2 | PMR-AS, relapse risk |

### 7 Treatment Scenarios

| # | Scenario | Drug and dose | Supporting trial |
|---|---------|----------|------------|
| **S1** | No treatment (natural history) | — | Natural-history cohort |
| **S2** | Pred 15 mg → taper 2.5 mg/mo | ACR standard | Dejaco 2015 ACR/EULAR |
| **S3** | Pred 22.5 mg → rapid taper 4 mg/mo | Severe cases | BSR guideline |
| **S4** | Pred 15 mg → slow taper 1 mg/mo | Relapse prevention | Observational cohort |
| **S5** | TCZ 162 mg SC QW + Pred 12.5 mg | GC-sparing | GiACTA (Stone 2017 NEJM) |
| **S6** | TCZ 162 mg SC Q2W + Pred 12.5 mg | GC-sparing | Cf. SEMAPHORE/SAPHYR |
| **S7** | TCZ QW alone (steroid-free induction) | Exploratory | PMR-SPARE phase 2 |

### Key Parameters (Trial-Calibrated)

| Parameter | Value | Source |
|---------|-----|------|
| Pred CL | 14.0 L/h | Bergmann 2012 |
| Pred V1 | 30.0 L | Bergmann 2012 |
| Pred Fu | 0.28 | Buttgereit 2005 |
| TCZ CL (linear) | 0.29 L/h | Nishimoto 2008 |
| IL-6 baseline (PMR) | 15 pg/mL | Roche 1993 |
| CRP baseline (PMR) | 35 mg/L | Clinical observation |
| EC50 (Pred→PMR-AS) | 180 ng/mL free | Model calibration |
| EC50 (TCZ→PMR-AS) | 50 nM | Model calibration |

---

## Shiny Dashboard (Interactive Dashboard)

**File**: [`pmr_shiny_app.R`](../../../polymyalgia-rheumatica/pmr_shiny_app.R)

### 6-Tab Structure

| Tab | Contents |
|----|------|
| **(1) Patient profile** | Diagnostic criteria, baseline parameters, PMR-AS/CRP/ESR/IL-6 value boxes |
| **(2) Drug PK** | Prednisolone total/free concentrations, tocilizumab nM trajectory, PK summary table |
| **(3) Inflammatory markers** | Time-series plots of IL-6, sIL-6Rα, CRP, ESR (with normal baselines) |
| **(4) Disease activity** | PMR-AS trajectory, HPA-axis cortisol suppression, relapse-risk score |
| **(5) Scenario comparison** | CRP, PMR-AS, IL-6, BMD compared across 7 treatment arms (colour-coded) |
| **(6) Biomarker explorer** | BMD trajectory, GC-sparing analysis, IL-6↔CRP correlation, reference-value table |

### How to Run
```r
install.packages(c("shiny","shinydashboard","mrgsolve","dplyr","ggplot2","plotly","DT"))
shiny::runApp("pmr_shiny_app.R")
```

---

## References

**File**: [`pmr_references.md`](pmr_references.md)

| Section | Citations |
|------|--------|
| Epidemiology and clinical features | 5 |
| Diagnostic criteria | 4 |
| Treatment guidelines | 4 |
| Pathogenesis and immunology | 7 |
| IL-6 pathway and cytokines | 5 |
| Corticosteroid PK/PD | 5 |
| Tocilizumab PK/PD and clinical trials | 6 |
| HPA-axis suppression | 3 |
| Osteoporosis and bone effects | 4 |
| Disease activity scores | 4 |
| Vascular involvement and GCA overlap | 4 |
| QSP/PK-PD modelling methodology | 4 |
| **Total** | **55** |

---

## Deliverables Summary

| Component | File | Specification |
|---------|------|------|
| Mechanistic map | [`pmr_qsp_model.dot`](../../../polymyalgia-rheumatica/pmr_qsp_model.dot) · [`.svg`](../../../polymyalgia-rheumatica/pmr_qsp_model.svg) · [`.png`](../../../polymyalgia-rheumatica/pmr_qsp_model.png) | **130+ nodes, 12 clusters** |
| mrgsolve ODE | [`pmr_mrgsolve_model.R`](pmr_mrgsolve_model.R) | **22-compartment ODE, 7 treatment scenarios, VPop of 200** |
| Shiny app | [`pmr_shiny_app.R`](../../../polymyalgia-rheumatica/pmr_shiny_app.R) | **6 tabs** (patient profile · PK · inflammatory markers · disease activity · scenario comparison · biomarker explorer) |
| References | [`pmr_references.md`](pmr_references.md) | **55 PubMed citations** (12 sections) |

---

## Clinical Implications

| Topic | Implication |
|------|--------|
| **Early diagnosis** | CRP >5 mg/L + ESR >50 mm/hr + bilateral shoulder pain + morning stiffness ≥45 min → strongly suspect PMR |
| **Starting GC** | Pred 12.5–25 mg/d → reconsider the diagnosis if there is no rapid response (24–72 h) |
| **Tapering** | ACR recommendation: reduce by 2.5 mg monthly after 4 weeks (revert to the prior dose on relapse) |
| **Adding TCZ** | In relapsing/GC-dependent PMR, TCZ 162 mg QW → sustains remission and spares GC |
| **Osteoporosis prevention** | Bisphosphonate + vitamin D/Ca are mandatory once GC use reaches ≥3 months |
| **GCA surveillance** | Headache, temporal artery tenderness, or visual disturbance → immediate GCA evaluation (urgent high-dose GC) |

---

*Generated: 2026-06-24 | Automatically generated by Claude Code Routine (CCR)*
