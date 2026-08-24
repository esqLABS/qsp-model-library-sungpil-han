# Hemophilia B (Christmas Disease) — QSP Model

**Directory:** `hemophilia-b/` | **Abbreviation:** HB | **Model version:** 1.0.0 | **Date:** 2026-07-01

---

## Disease Overview

Hemophilia B is an X-linked recessive bleeding disorder caused by deficiency
of coagulation factor IX (FIX). It is about 4-5-fold rarer than hemophilia A
(roughly 1 in 25,000-30,000 male births). It is classified as severe
(FIX <1 IU/dL), moderate (1-5 IU/dL), or mild (>5 IU/dL), and repeated
intra-articular bleeding (hemarthrosis) can progress to irreversible
haemophilic arthropathy.

Core pathophysiology:
- **F9 gene defect**: missense/nonsense/large-deletion variants in the *F9*
  gene at Xq27.1 impair hepatocyte FIX synthesis · secretion; the Hemophilia
  B Leiden variant is notable for spontaneous improvement after puberty
- **Vitamin K-dependent gamma-carboxylation**: the FIX Gla domain must be
  carboxylated via the GGCX/VKORC1 pathway to enable phospholipid-surface
  binding and coagulant activity; warfarin blocks this pathway
- **Intrinsic tenase deficiency**: FIXa fails to assemble the tenase complex
  with FVIIIa cofactor on the Ca2+/phospholipid surface, so FX activation and
  thrombin generation are reduced
- **Inhibitors**: anti-FIX IgG antibodies develop in about 1-3% of severe HB
  (lower than in HA), with higher risk in null-mutation genotypes, and carry
  the HB-specific risk of **nephrotic syndrome · anaphylaxis** as immune
  tolerance induction (ITI) complications
- **FIX-Padua (R338L) variant**: a naturally occurring hyperactive variant
  (~8-fold comparative activity) used as the transgene in current AAV gene
  therapy vectors

---

## Core Treatment Paradigms

| Strategy | Drug | Mechanism | Key trial |
|------|------|------|-------------|
| **Standard half-life (SHL) FIX replacement** | Nonacog alfa (BeneFIX) · Rixubis | Direct FIX replacement, IV | Registry population PK |
| **Extended half-life (EHL) — Fc fusion** | Eftrenonacog alfa (Alprolix) | FcRn recycling → t½ ~82h | B-LONG (Powell 2013 NEJM) |
| **Extended half-life (EHL) — albumin fusion** | Albutrepenonacog alfa (Idelvion) | Albumin fusion → t½ ~102-104h | PROLONG-9FP (Santagostino 2016) |
| **Extended half-life (EHL) — glycoPEGylation** | Nonacog beta pegol (Refixia) | GlycoPEGylation → t½ ~93h | pathfinder2 (Collins 2014) |
| **TFPI inhibition (non-replacement)** | **Concizumab** | Anti-TFPI antibody → amplifies the extrinsic pathway | explorer7/8 (Shapiro/Chowdary) |
| **TFPI inhibition (non-replacement)** | **Marstacimab** | Anti-TFPI antibody, weekly SC | BASIS (Pipe 2023 NEJM) |
| **Antithrombin inhibition (non-replacement)** | **Fitusiran** | siRNA knockdown of AT mRNA → thrombin↑ | ATLAS-A/B, ATLAS-INH (Young 2023) |
| **Gene therapy** | Etranacogene dezaparvovec (Hemgenix) | AAV5 vector, expresses FIX-Padua | HOPE-B (Pipe 2023 NEJM) |
| **Gene therapy** | Fidanacogene elaparvovec (Beqvez) | AAV-Rh74var vector, expresses FIX-Padua | BENEGENE-2 (Cuker 2024) |

---

## Mechanistic Map

[![Hemophilia B QSP Map](hb_qsp_model.png)](hb_qsp_model.svg)

> Click to view the high-resolution SVG file.

**Composition:** 125 nodes · 13 subgraph clusters

| Cluster | Key content | Therapeutic target |
|---------|---------|---------|
| **F9 gene · hepatic synthesis · vitamin K cycle** | F9 variant types · GGCX · VKORC1 · Gla carboxylation · FIX-Padua | Gene therapy transgene design |
| **Intrinsic pathway · tenase complex** | FXII/FXI/FIXa · FVIIIa cofactor · tenase · FX activation | FIX replacement therapy target |
| **Common pathway · thrombin generation · fibrin** | Prothrombinase · thrombin burst · FXIIIa · TAFI · plasmin | ETP clinical metric |
| **Bleeding phenotype · haemophilic arthropathy** | Spontaneous hemarthrosis · synovial iron deposition · MMP · target joint | Prophylaxis goal of ABR <3 |
| **FIX replacement PK (SHL/EHL)** | BeneFIX/Alprolix/Idelvion/Refixia · FcRn recycling | Trough-maintenance strategy |
| **Non-replacement rebalancing therapy** | Concizumab · marstacimab (TFPI) · fitusiran (AT) | Factor-independent hemostasis |
| **AAV gene therapy** | Hepatocyte transduction · episomal expression · capsid immune response | One-time curative-intent approach |
| **Inhibitor immunology** | Anti-FIX IgG · Bethesda titer · ITI · nephrotic syndrome · anaphylaxis | HB-specific ITI risk management |
| **Pharmacodynamics · laboratory biomarkers** | one-stage/chromogenic assay · TGA · reagent discrepancy | EHL assay-error correction |
| **Clinical endpoints · PRO** | ABR/AJBR · HJHS · Haem-A-QoL · surgical hemostasis | Zero-bleed phenotype |
| **Lifecycle · special-population correction factors** | Neonatal physiologically low FIX · pediatric Vd · hepatic/renal function | Dose individualisation |
| **Long-term comorbidity · history-based risk** | Chronic pain · orthopedic surgery · transfusion history · mortality | Quality of life · life expectancy |
| **Adjunctive hemostatic measures** | Tranexamic acid · desmopressin (ineffective) · physical therapy | Local/adjunctive management |

---

## mrgsolve ODE Model Specifications

**File:** [`hb_mrgsolve_model.R`](hb_mrgsolve_model.R)

### Compartments — 24

| # | Compartment | Unit | Description |
|---|------|------|------|
| 1 | `FIX_C` | IU/dL | FIX central compartment (SHL) |
| 2 | `FIX_P` | IU/dL | FIX peripheral compartment (SHL) |
| 3 | `FIXe_C` | IU/dL | FIX central compartment (EHL: Fc/albumin/glycoPEGylation option) |
| 4 | `FIXe_P` | IU/dL | FIX peripheral compartment (EHL) |
| 5 | `AAV_Vector` | rel. | Circulating AAV vector genomes |
| 6 | `Transduced_Hep` | 0-1 | Fraction of transduced hepatocytes |
| 7 | `Capsid_Antigen` | rel. | Transient capsid antigen-presentation pool (drives immune response only, ~5-week half-life) |
| 8 | `Transgene_Expr` | IU/dL | Endogenous FIX-Padua transgene expression |
| 9 | `Capsid_Immune` | 0-1 | Capsid-specific immune activation |
| 10 | `ALT_level` | fold | Fold rise in ALT (hepatotoxicity marker) |
| 11 | `CONC_SC` | mg | Concizumab SC depot |
| 12 | `CONC_C` | ng/mL | Concizumab central compartment |
| 13 | `MARS_SC` | mg | Marstacimab SC depot |
| 14 | `MARS_C` | ng/mL | Marstacimab central compartment |
| 15 | `FITU_SC` | mg | Fitusiran SC depot |
| 16 | `FITU_C` | mg/L | Fitusiran central compartment |
| 17 | `AT_mRNA` | rel. | Antithrombin mRNA (baseline = 1) |
| 18 | `AT_prot` | rel. | Antithrombin protein (baseline = 1) |
| 19 | `Inhibitor` | BU/mL | Anti-FIX inhibitor titer |
| 20 | `Thrombin_ETP` | norm. | Thrombin generation potential |
| 21 | `CumBleeds` | count | Cumulative bleed count |
| 22 | `JointScore` | 0-100 | Haemophilic arthropathy score |
| 23 | `Synovitis` | 0-1 | Synovial inflammation index |
| 24 | `QoL` | 0-1 | Quality of life (utility derived from Haem-A-QoL) |

### Treatment Scenarios — 10

| Scenario | Regimen | Route of administration | Clinical calibration |
|---------|------|---------|---------|
| 1 | No prophylaxis | On-demand | ABR ~28 (untreated severe HB) |
| 2 | SHL-rFIX prophylaxis | 40 IU/kg 2×/week IV | Registry population PK |
| 3 | EHL-rFIX-Fc prophylaxis | 50 IU/kg Q7-10 days IV | B-LONG 2013; t½ ~82h |
| 4 | EHL-rFIX-albumin prophylaxis | 75 IU/kg Q14 days IV | PROLONG-9FP 2016; t½ ~102-104h |
| 5 | GlycoPEGylated-rFIX prophylaxis | 40 IU/kg Q7 days IV | pathfinder2 2014; t½ ~93h |
| 6 | Concizumab daily SC | Loading 210mg + maintenance 15mg/day SC | explorer7/8 |
| 7 | Marstacimab weekly SC | Loading 300mg + maintenance 150mg/week SC | BASIS 2023 |
| 8 | Fitusiran monthly SC | 50 mg/month SC | ATLAS-A/B 2023 |
| 9 | AAV gene therapy, single dose | Single IV vector dose + steroid-response management | HOPE-B 2023 |
| 10 | Inhibitor-positive patient | ITI protocol + on-demand bypassing-agent management | Null-mutation high-risk scenario |

---

## Shiny Dashboard

**File:** [`hb_shiny_app.R`](hb_shiny_app.R) | **Tabs: 8**

| Tab | Key features |
|----|---------|
| **1. Patient Profile** | Weight · FIX severity · genotype (null mutation) · inhibitor status · treatment selection; Value Box summary |
| **2. FIX PK** | SHL/EHL FIX concentration-time curves (linear/log); 1%·15% trough reference lines |
| **3. Non-Factor Rebalancing** | Concizumab/marstacimab concentrations; tracking of AT protein reduction by fitusiran |
| **4. Gene Therapy** | Time course of AAV transgene expression; capsid immune response and ALT hepatotoxicity kinetics |
| **5. PD Core Metrics** | ETP thrombin generation metric; total effective FIX-equivalent (replacement + transgene + TFPI) |
| **6. Bleed Risk & Arthropathy** | Instantaneous ABR over time; combined joint-score + QoL plot; Value Box |
| **7. Scenario Comparison** | Simultaneous comparison of 9 treatment options; long-term ABR + joint-score trends |
| **8. Inhibitor & Biomarkers** | FIX activity vs. ETP scatter plot; inhibitor titer kinetics; clinical outcome summary table |

```r
# How to run
install.packages(c("shiny", "bslib", "plotly", "dplyr", "tidyr", "ggplot2"))
shiny::runApp("hemophilia-b/hb_shiny_app.R")
```

---

## References

**File:** [`hb_references.md`](hb_references.md)

71 verified PubMed citations — key sections:
- Disease overview · epidemiology (WFH guidelines)
- F9 gene · molecular genetics · FIX biology (including FIX-Padua)
- Coagulation cascade · thrombin generation
- FIX pharmacokinetics — SHL products
- FIX pharmacokinetics — EHL products (B-LONG · PROLONG-9FP · pathfinder2)
- Non-replacement rebalancing therapy (concizumab · marstacimab · fitusiran)
- AAV gene therapy (HOPE-B · BENEGENE-2)
- Inhibitors · immune complications (nephrotic syndrome · anaphylaxis)
- Haemophilic arthropathy
- Clinical endpoints · QoL · guidelines

---

## Summary of Key Clinical Parameters

| Parameter | Severe HB (untreated) | SHL-rFIX prophylaxis | EHL-rFIX-Fc prophylaxis | AAV gene therapy |
|---------|----------------|--------------|-----------------|--------------|
| ABR | ~28/year | ~13-16/year | ~7-11/year | ~4-6/year |
| FIX trough/expression | <1 IU/dL | 1-5 IU/dL | 5-10 IU/dL | ~15-30 IU/dL (sustained) |
| ETP (vs. normal) | ~15% | ~30-40% | ~40-55% | ~60-75% |
| QoL (utility) | 0.55-0.65 | 0.70-0.80 | 0.75-0.85 | 0.80-0.90 |

---

## Deliverables Summary

| Deliverable | File | Contents |
|--------|------|------|
| 🗺️ Mechanistic map | [`hb_qsp_model.dot/.svg/.png`](hb_qsp_model.svg) | **125 nodes, 13 clusters** (F9 gene · intrinsic pathway · common pathway · bleeding phenotype · FIX replacement PK · non-replacement rebalancing · AAV gene therapy · inhibitor immunology · biomarkers · clinical endpoints · lifecycle · long-term comorbidity · adjunctive hemostasis) |
| ⚙️ mrgsolve ODE | [`hb_mrgsolve_model.R`](hb_mrgsolve_model.R) | **24-compartment ODE** (FIX SHL/EHL 4 compartments · AAV gene therapy 5 compartments · concizumab/marstacimab 4 compartments · fitusiran/AT 4 compartments · inhibitor · ETP · CumBleeds · joint score · QoL · Synovitis), **10 treatment scenarios** |
| 📊 Shiny app | [`hb_shiny_app.R`](hb_shiny_app.R) | **8 tabs** (patient profile · FIX PK · non-replacement rebalancing · gene therapy · PD core metrics · bleed risk/arthropathy · scenario comparison · inhibitor/biomarkers), bslib darkly, plotly, built-in ODE simulator |
| 📚 References | [`hb_references.md`](hb_references.md) | **71 PubMed citations** (HOPE-B · B-LONG · PROLONG-9FP · pathfinder2 · BASIS · ATLAS-A/B, etc.) |
