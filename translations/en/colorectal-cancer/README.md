# Colorectal Cancer (CRC) — QSP Model

> **Category**: Oncology · Gastroenterology | **Abbreviation**: CRC | **Version**: 1.0 | **Date**: 2026-06-23

---

## Disease Overview

Colorectal cancer (CRC) is a malignant tumour arising in the colon or
rectum, and is the third most common and second most lethal cancer
worldwide (GLOBOCAN 2020: about 1.93 million new patients and 930,000
deaths per year). Most CRC follows the **Vogelstein sequence**
(APC→KRAS→SMAD4→TP53), which progresses from normal epithelium to adenoma
to invasive carcinoma, while about 15% progresses rapidly via the MSI-H
pathway.

---

## Key Pathophysiological Pathways

| Pathway | Key molecules | Clinical significance |
|------|-----------|-----------|
| **Wnt/β-catenin** | APC(mut), β-catenin, TCF/LEF | APC mutation in ~80%, c-Myc/cyclin D1 overexpression |
| **RAS/MAPK** | KRAS(40%), NRAS(5%), BRAF V600E(10%) | Predicts response to anti-EGFR therapy |
| **PI3K/AKT** | PIK3CA(20%), PTEN loss | mTOR activation, resistance mechanism |
| **TP53** | R175H/R248W/R273H hotspots | Cell-cycle arrest, impaired apoptosis |
| **TGF-β/EMT** | SMAD4(35%), ZEB1/TWIST | Metastasis, poor prognosis |
| **MSI/MMR** | MLH1/MSH2/MSH6/PMS2 loss | High tumour mutational burden → sensitivity to immunotherapy |

---

## Consensus Molecular Subtypes (CMS)

| CMS | Key characteristics | Treatment strategy |
|-----|-----------|-----------|
| **CMS1** (MSI Immune) | MSI-H, BRAF V600E, high immune activity | Anti-PD-1/PD-L1 (pembrolizumab) |
| **CMS2** (Canonical) | High CIN, WNT/MYC activation | FOLFOX + anti-EGFR (RAS-WT) |
| **CMS3** (Metabolic) | KRAS mutation, abnormal lipid metabolism | FOLFIRI/FOLFOX ± BEV |
| **CMS4** (Mesenchymal) | TGF-β activity, EMT, stromal infiltration | Anti-angiogenic ± anti-EGFR |

---

## Model Architecture

### Mechanistic Map

[![CRC QSP Map](../../../colorectal-cancer/crc_qsp_model.png)](../../../colorectal-cancer/crc_qsp_model.svg)

- **Node count**: 130+ (12 subgraph clusters)
- **Clusters**: Wnt · RAS/MAPK · PI3K/AKT · TP53/Apoptosis · Cell Cycle · TGF-β/EMT · TME · Immune Checkpoints · Angiogenesis · MSI/MMR · Drug PK/PD · Clinical Endpoints

### mrgsolve ODE Model

**20 ODE compartments:**

| Compartment group | Variables | Modelled content |
|-----------|------|-------------|
| 5-FU PK | FU1, FU2, FU_ic | 2-compartment + intracellular activation, DPD metabolism |
| Oxaliplatin | OX1, OX_DNA | Platinum-DNA adduct formation/repair |
| Irinotecan | IRI1, SN38, SN38G | Reflects UGT1A1*28 polymorphism |
| Bevacizumab | BEV1, BEV_VEGF | TMDD model, VEGF neutralisation |
| Cetuximab | CTX1 | EGFR occupancy |
| Pembrolizumab | PEM1 | PD-1 occupancy |
| Tumour | Ts, Tr | Logistic growth of sensitive/resistant cells |
| Biomarkers | CEA, ctDNA | Linked to tumour burden |
| Immune | CD8eff | T-cell activity (PD-1 blockade) |
| VEGF, EGFR_occ, PD1_occ | 3 | Target-molecule kinetics |

**Treatment Scenarios (7):**

| # | Scenario | Supporting trial |
|---|---------|--------------|
| 1 | No treatment (natural history) | — |
| 2 | FOLFOX6 | MOSAIC (André 2004) |
| 3 | FOLFIRI | GERCOR (Tournigand 2004) |
| 4 | FOLFOX + bevacizumab | NO16966 (Saltz 2008) |
| 5 | FOLFIRI + cetuximab (RAS-WT) | CRYSTAL (Van Cutsem 2009) |
| 6 | FOLFIRI + bevacizumab | TRIBE (Falcone 2013) |
| 7 | Pembrolizumab (MSI-H) | KEYNOTE-177 (André 2020) |

### Shiny App Structure (6 Tabs)

| Tab | Contents |
|----|------|
| Patient profile | Age, BSA, RAS/BRAF/MSI status, CMS subtype, treatment selection |
| PK profile | Concentration-time for 5-FU/oxaliplatin/irinotecan/bevacizumab/cetuximab/pembrolizumab |
| Key PD measures | RECIST change in tumour diameter, sensitive/resistant cells, CD8 T-cell activity |
| Clinical endpoints | Best overall response, estimated PFS, RECIST classification |
| Scenario comparison | Simultaneous comparison graph and summary table across multiple regimens |
| Biomarkers | CEA, ctDNA, free VEGF, EGFR/PD-1 occupancy |

---

## Drug Summary

| Drug | Mechanism | Indication |
|------|------|--------|
| 5-FU/LV | TS inhibition → dTTP depletion | First-line for all CRC |
| Oxaliplatin | Platinum-DNA adducts | FOLFOX backbone |
| Irinotecan | Top1 poison (SN-38) | FOLFIRI backbone |
| Bevacizumab | Anti-VEGF-A (TMDD) | First-line metastatic CRC |
| Cetuximab | Anti-EGFR (RAS-WT) | CRYSTAL/OPUS |
| Pembrolizumab | Anti-PD-1 | First-line MSI-H/dMMR CRC |
| Regorafenib | Multi-kinase TKI | Third line and beyond |

---

## Reference Files

| File | Description |
|------|------|
| [`crc_qsp_model.dot`](../../../colorectal-cancer/crc_qsp_model.dot) | Graphviz mechanistic map source |
| [`crc_qsp_model.svg`](../../../colorectal-cancer/crc_qsp_model.svg) | Vector graphic (zoomable) |
| [`crc_qsp_model.png`](../../../colorectal-cancer/crc_qsp_model.png) | 150 dpi raster image |
| [`crc_mrgsolve_model.R`](../../../colorectal-cancer/crc_mrgsolve_model.R) | mrgsolve ODE model (20 compartments, 7 scenarios) |
| [`crc_shiny_app.R`](../../../colorectal-cancer/crc_shiny_app.R) | Shiny interactive dashboard (6 tabs) |
| [`crc_references.md`](../../../colorectal-cancer/crc_references.md) | 40 references (grouped by section) |
