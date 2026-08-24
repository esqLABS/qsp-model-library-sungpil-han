# Hemophilia A — QSP Model

**Directory:** `hemophilia-a/` | **Abbreviation:** HA | **Model version:** 1.0.0 | **Written:** 2026-06-23

---

## Disease Overview

Hemophilia A (HA) is an X-linked recessive bleeding disorder caused by deficiency of coagulation factor VIII (FVIII). It is classified as severe (FVIII <1 IU/dL), moderate (1–5 IU/dL), or mild (>5 IU/dL), and occurs at a frequency of about 1 in 10,000 males worldwide. Repeated severe intra-articular bleeding (haemarthrosis) progresses to irreversible haemophilic arthropathy, causing extreme deterioration in quality of life.

Core pathophysiology:
- **Intrinsic-pathway defect**: FVIII acts as a cofactor for FIXa in forming the tenase complex; deficiency reduces FXa generation → insufficient thrombin generation
- **Unstable fibrin clot**: absence of a thrombin burst reduces fibrin polymer strength
- **Inhibitors**: about 30% of severe HA patients develop anti-FVIII IgG antibodies → requiring bypassing therapy

---

## Core Treatment Paradigms

| Strategy | Drug | Mechanism | Key clinical trial |
|------|------|------|-------------|
| **FVIII replacement** | Standard half-life (SHL) FVIII | Direct FVIII replacement | Manco-Johnson 2007 NEJM |
| **Extended half-life** | Fc-fusion/PEGylated FVIII | ↑ t½ via FcRn recycling | A-LONG (Mahlangu 2014) |
| **Non-factor prophylaxis** | **Emicizumab** (Hemlibra) | FIXa–FX bispecific antibody, mimics FVIII | HAVEN 1/3/4 (Oldenburg 2017) |
| **Antithrombin inhibition** | **Fitusiran** | siRNA knockdown of AT mRNA → ↑ thrombin | ATLAS-INH (Young 2023) |
| **TFPI inhibition** | **Marstacimab** | Anti-TFPI antibody → amplifies the extrinsic pathway | BASIS study |
| **Gene therapy** | Valoctocogene roxaparvovec | AAV5-vector FVIII expression | HOPE-B (Ozelo 2022) |

---

## Mechanistic Map

[![Hemophilia A QSP Map](ha_qsp_model.png)](ha_qsp_model.svg)

> Click to view the high-resolution SVG file.

**Composition:** 167 nodes · 10 subgraph clusters

| Cluster | Key content | Treatment target |
|---------|---------|---------|
| **Vascular injury · platelets** | TF exposure · vWF · GPIb · GPIIbIIIa · TXA2 · ADP | Initial plug from platelet activation |
| **Extrinsic pathway** | TF · FVIIa · TF/FVIIa complex · TFPI inhibition | Marstacimab (TFPI neutralisation) |
| **FVIII biology · intrinsic pathway** | F8 gene · vWF protection · FVIIIa · Xase complex | FVIII replacement therapy · EHL |
| **Common pathway · thrombin generation** | Prothrombinase · thrombin burst · FXIIIa · fibrin | TGA/ETP clinical indices |
| **Natural anticoagulant mechanisms** | AT · protein C/S · EPCR · thrombomodulin · tPA | Fitusiran (AT reduction) |
| **Inhibitor immunology** | TH cells · B cells · anti-FVIII IgG · BU titre | Immune tolerance induction (ITI) · emicizumab |
| **Drug PK/PD** | SHL/EHL FVIII · emicizumab · fitusiran · desmopressin | Exposure-response of each drug |
| **Bleeding phenotype** | Haemarthrosis · intramuscular · intracranial · GI bleeding · ABR | Prophylaxis target ABR <3 |
| **Haemophilic arthropathy** | Synovial iron deposition · ROS · cartilage destruction · Pettersson score | Joint protection through prophylaxis |
| **Clinical endpoints** | FVIII trough ≥1% · ETP · ABR · EQ-5D · joint score | Zero-bleed phenotype |

---

## mrgsolve ODE Model Specifications

**File:** [`ha_mrgsolve_model.R`](ha_mrgsolve_model.R)

### Compartments — 16

| # | Compartment | Unit | Description |
|---|------|------|------|
| 1 | `FVIII_C` | IU/dL | FVIII central compartment |
| 2 | `FVIII_P` | IU/dL | FVIII peripheral compartment |
| 3 | `EMIC_SC` | mg | Emicizumab SC depot |
| 4 | `EMIC_C` | mg/L | Emicizumab central compartment |
| 5 | `EMIC_P` | mg/L | Emicizumab peripheral compartment |
| 6 | `FITU_SC` | mg | Fitusiran SC depot |
| 7 | `FITU_C` | mg/L | Fitusiran central compartment |
| 8 | `AT_mRNA` | rel. | Antithrombin mRNA (baseline = 1) |
| 9 | `AT_prot` | rel. | Antithrombin protein (baseline = 1) |
| 10 | `Inhibitor` | BU/mL | FVIII inhibitor titre |
| 11 | `Thrombin_ETP` | norm. | Thrombin generation potential |
| 12 | `CumBleeds` | count | Cumulative bleed count |
| 13 | `JointScore` | 0-100 | Pettersson joint score |
| 14 | `QoL` | 0-1 | Quality of life (EQ-5D) |
| 15 | `Synovitis` | 0-1 | Synovial inflammation index |
| 16 | `FVIII_eff` | IU/dL | Effective FVIII activity (FVIII + emicizumab) |

### Treatment Scenarios — 7

| Scenario | Regimen | Administration | Clinical calibration |
|---------|------|---------|---------|
| 1 | No prophylaxis | On demand | ABR ~30 (untreated severe HA) |
| 2 | SHL-FVIII prophylaxis | 25 IU/kg 3×/week IV | Manco-Johnson 2007; ABR ~3–4 |
| 3 | EHL-FVIII prophylaxis | 50 IU/kg every 3-4 days IV | A-LONG 2014; ABR ~2–3 |
| 4 | Emicizumab Q1W | 1.5 mg/kg SC Q1W (3 mg/kg loading ×4) | HAVEN 3 2018; ABR 1.5 |
| 5 | Emicizumab Q4W | 6 mg/kg SC Q4W (after loading) | HAVEN 4 2019; ABR 2.4 |
| 6 | Fitusiran Q1M | 80 mg SC Q1M | ATLAS-INH 2023; ABR ~0 |
| 7 | FVIII + emicizumab | Combination | Surgery/high-risk periods |

---

## Shiny Dashboard (Interactive Dashboard)

**File:** [`ha_shiny_app.R`](ha_shiny_app.R) | **Tabs: 6**

| Tab | Key features |
|----|---------|
| **1. Patient Profile** | Weight · FVIII severity · inhibitor status · treatment selection; value-box summary |
| **2. FVIII PK** | FVIII concentration-time curve (linear/log); 1%/15% trough reference lines; PK summary statistics |
| **3. PD Core Metrics** | ETP thrombin-generation index; tracking of AT protein reduction by fitusiran |
| **4. Bleed Risk & ABR** | Instantaneous ABR over time; combined joint-score + QoL plot; value box |
| **5. Scenario Comparison** | Simultaneous comparison of 6 scenarios; long-term ABR + joint-score trends |
| **6. Biomarkers** | FVIII activity vs. ETP scatter plot; inhibitor-titre dynamics; clinical-outcome summary table |

```r
# How to run
install.packages(c("shiny", "bslib", "plotly", "dplyr", "tidyr", "ggplot2"))
shiny::runApp("hemophilia-a/ha_shiny_app.R")
```

---

## References

**File:** [`ha_references.md`](ha_references.md)

55 PubMed citations — key sections:
- WFH guidelines · epidemiology (5)
- FVIII biology · molecular pathology (5)
- Thrombin generation (5)
- FVIII pharmacokinetics (5)
- SHL/EHL FVIII clinical trials (4)
- Emicizumab HAVEN 1/3/4 (6)
- Fitusiran ATLAS-INH (4)
- Marstacimab and others (3)
- Inhibitor development (4)
- Haemophilic arthropathy (4)
- QoL · patient-reported outcomes (3)
- Gene therapy (2)
- QSP/PK-PD modelling (5)

---

## Summary of Key Clinical Parameters

| Parameter | Severe HA (untreated) | SHL-FVIII prophylaxis | Emicizumab Q1W |
|---------|----------------|--------------|-------------|
| ABR | ~30/year | ~3–4/year | ~1.5/year |
| FVIII trough | <1 IU/dL | 1–5 IU/dL | ~15 IU/dL equivalent |
| ETP (vs. normal) | ~15% | ~40–60% | ~70–80% |
| QoL (EQ-5D) | 0.55–0.65 | 0.75–0.85 | 0.85–0.92 |
| 10-year joint score | >50 | 15–25 | 8–15 |

---

## Deliverables Summary

| Deliverable | File | Content |
|--------|------|------|
| 🗺️ Mechanistic map | [`ha_qsp_model.dot/.svg/.png`](ha_qsp_model.svg) | **167 nodes, 10 clusters** (vascular injury · extrinsic pathway · FVIII/intrinsic pathway · common pathway · anticoagulant mechanisms · inhibitor immunology · drug PK/PD · bleeding phenotype · arthropathy · clinical endpoints) |
| ⚙️ mrgsolve ODE | [`ha_mrgsolve_model.R`](ha_mrgsolve_model.R) | **16-compartment ODE** (FVIII 2 compartments · emicizumab 3 compartments · fitusiran 2 compartments · AT mRNA/protein · inhibitor · ETP · CumBleeds · joint score · QoL · synovitis · FVIII_eff), **7 treatment scenarios** |
| 📊 Shiny app | [`ha_shiny_app.R`](ha_shiny_app.R) | **6 tabs** (patient profile · FVIII PK · PD core metrics · bleed risk · scenario comparison · biomarkers), bslib darkly, plotly, built-in ODE simulator |
| 📚 References | [`ha_references.md`](ha_references.md) | **55 PubMed citations** (HAVEN 1/3/4 · ATLAS-INH · A-LONG · HOPE-B · Manco-Johnson 2007, etc.) |
