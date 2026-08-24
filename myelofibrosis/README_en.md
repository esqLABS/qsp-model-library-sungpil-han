# Myelofibrosis (MF) QSP Model

> **Directory**: `myelofibrosis/` | **Date**: 2026-06-23 | **Category**: Oncology/Haematology

[![MF QSP Map](mf_qsp_model.png)](mf_qsp_model.svg)

---

## Disease Overview

**Myelofibrosis (MF)** is a chronic myeloproliferative neoplasm (MPN) arising from haematopoietic stem cells (HSC), characterised by collagen/reticulin fibrosis within the bone marrow, ineffective haematopoiesis, abnormal megakaryocyte proliferation, and **extramedullary haematopoiesis (EMH)** centred on the spleen and liver. It is classified as Primary MF (PMF) and MF transformed from polycythaemia vera/essential thrombocythaemia (Post-PV/ET MF).

- **Incidence**: approximately 3,000~4,000 new cases per year in the United States, with a median age at onset of 67 years
- **5-year survival**: ~35% for intermediate-2/high-risk groups, ~70% for low-risk groups (by DIPSS Plus)
- **AML transformation rate**: 5~10% per year, markedly increased in the presence of high-risk molecular abnormalities (ASXL1, EZH2, SRSF2, IDH1/2, TP53)

### Molecular Pathogenesis

| Mutation | Frequency | Mechanism |
|---------|-----|------|
| **JAK2 V617F** | ~55–60% | Constitutive JAK2 activation → STAT3/5 hyperphosphorylation → unrestrained proliferation |
| **CALR exon 9 (Type 1/2)** | ~25–30% | Aberrant activation of the MPL receptor → enhanced JAK/STAT signalling |
| **MPL W515L/K** | ~5–10% | TPO receptor mutation → enhanced JAK/STAT signalling |
| **Triple-negative** | ~10% | Other, unidentified mechanisms |
| **HMR mutations** | Various | ASXL1 (30–40%), EZH2 (5%), SRSF2 (10%), IDH1/2 (5%) — high-risk |

---

## Model Architecture

### Four Deliverables

| Deliverable | File | Specification |
|--------|------|----------|
| 🗺️ **Mechanistic map** | [`mf_qsp_model.dot`](mf_qsp_model.dot) / [`.svg`](mf_qsp_model.svg) / [`.png`](mf_qsp_model.png) | **204 nodes, 12 clusters, 224 edges** |
| ⚙️ **mrgsolve ODE model** | [`mf_mrgsolve_model.R`](mf_mrgsolve_model.R) | **23-compartment ODE**, 6 treatment scenarios |
| 📊 **Shiny dashboard** | [`mf_shiny_app.R`](mf_shiny_app.R) | **6-tab** interactive app |
| 📚 **References** | [`mf_references.md`](mf_references.md) | **36** PubMed citations |

---

### Mechanistic Map — 12 Clusters

| # | Cluster | Key components |
|---|---------|-------------|
| 1 | **Genetic/mutational drivers** | JAK2V617F, CALR T1/T2, MPL W515L/K, ASXL1, EZH2, SRSF2, IDH1/2, TP53, TET2, DNMT3A, DIPSS/MIPSS70 risk scores |
| 2 | **JAK/STAT signalling** | JAK1/2/TYK2, STAT1/3/5a/5b, pSTAT3/pSTAT5, SOCS1/3, SHP2, BCL2/MCL1/CCND1/MYC/PIM1 |
| 3 | **HSC niche** | LT-HSC → ST-HSC → MPP → CMP/GMP/MEP → MkP/BFU-E; CXCL12/CXCR4, SCF/c-KIT, TPO, EPO |
| 4 | **Bone marrow microenvironment/fibrosis** | MSC, osteoblasts/osteoclasts, abnormal megakaryocytes, TGF-β1/PDGF/bFGF/CTGF, collagen I/III, MF grade 0–3, osteosclerosis |
| 5 | **Cytokine storm** | IL-1β/6/8/10/12/13, TNF-α, IFN-γ, CXCL10, CCL2, NF-κB, mTOR, PI3K, RAS/ERK, NLRP3 |
| 6 | **Extramedullary haematopoiesis** | Circulating CD34+, splenic/hepatic EMH, splenic/hepatic volume, marrow failure |
| 7 | **Haematological outcomes** | Hgb dynamics, RBC production/destruction, platelets, WBC/neutrophils, circulating blasts, transfusion dependence |
| 8 | **Thrombosis/vasculature** | Platelet activation, TXA2, thrombin generation, TF/FXa cascade, PAI-1, endothelial dysfunction, DVT/PE/BCS/PVT |
| 9 | **Drug PK** | Ruxolitinib 2-compartment (ka/CL/V1/V2/Q/F1), fedratinib, pacritinib, momelotinib PK |
| 10 | **Drug PD** | JAK1/2 inhibition (Emax/IC50), pSTAT3/5 suppression, SVR35/TSS50 output, anaemia/thrombocytopenia side effects |
| 11 | **Clinical endpoints** | SVR35, TSS50, bone marrow histological response, CHR, OS, PFS, AML transformation, IWG-MRT/ELN 2023 criteria |
| L | **Legend** | Guide to node shapes and colours |

---

### mrgsolve ODE Model — 23 Compartments

```
Ruxolitinib PK   : DEPOT_RUX, CENT_RUX, PERI_RUX             (3 compartments)
Fedratinib PK    : DEPOT_FED, CENT_FED                        (2 compartments)
Pacritinib PK    : DEPOT_PAC, CENT_PAC                        (2 compartments)
BET inhibitor PK : CENT_BET  (pelabresib)                     (1 compartment)
JAK/STAT PD      : pSTAT3, pSTAT5                             (2 compartments)
Clonal dynamics  : NHSC (malignant HSC), NHSC_N (normal HSC)  (2 compartments)
Erythroid lineage: PROG_E, RET, RBC                           (3 compartments)
Megakaryocytic lineage: MEG_P, PLT                            (2 compartments)
Macroscopic measures: SPLEEN, FIBROSIS                        (2 compartments)
Cytokines        : IL6, TNF                                   (2 compartments)
Symptom score    : TSS                                        (1 compartment)
```

**Clinical Trials Used for Calibration**

| Scenario | Trial | Key outcome |
|---------|---------|---------|
| Ruxolitinib 20 mg BID | COMFORT-I (Verstovsek 2012) | SVR35 41.9%, TSS50 45.9% |
| Ruxolitinib 15 mg BID | COMFORT-I low-platelet cohort | SVR35 ~28% |
| Fedratinib 400 mg QD | JAKARTA (Pardanani 2015) | SVR35 36%, TSS50 36% |
| Pacritinib 200 mg BID | PERSIST-2 (Mesa 2017) | SVR35 18%, platelets <50×10⁹/L |
| Ruxolitinib + Pelabresib | MANIFEST-2 (Pemmaraju 2024) | SVR35 66% vs 35% |
| Untreated | Natural history | Spleen increases by ~10% per year |

---

### Shiny Dashboard — 6 Tabs

| Tab | Content |
|----|------|
| **① Patient Profile** | Age/sex/diagnosis, DIPSS Plus risk, mutations (JAK2/CALR/MPL), baseline spleen volume/Hgb/platelets, symptom score |
| **② Pharmacokinetics (PK)** | Drug selection, dose/interval input, blood concentration-time curve (Cp vs time), Cmax/Cmin/AUC/t½ table |
| **③ PD Biomarkers** | pSTAT3/pSTAT5 inhibition profile, JAK2 V617F VAF trend, IL-6/TNF-α cytokine dynamics |
| **④ Clinical Endpoints** | Spleen volume (SVR35 threshold), Hgb change, platelet count, TSS score, bone marrow fibrosis grade |
| **⑤ Treatment Comparison** | Multi-drug checkboxes, SVR35/TSS50/Hgb butterfly plot, waterfall spleen response |
| **⑥ Biomarker Dynamics** | JAK2 VAF–spleen response correlation, cytokine heatmap, KM curve for AML transformation risk |

---

## How to Run

### 1. Render the Mechanistic Map

```bash
# Graphviz must be installed
dot -Tsvg myelofibrosis/mf_qsp_model.dot -o myelofibrosis/mf_qsp_model.svg
dot -Tpng -Gdpi=150 myelofibrosis/mf_qsp_model.dot -o myelofibrosis/mf_qsp_model.png
```

### 2. Run the mrgsolve Model

```r
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr"))
library(mrgsolve)
source("myelofibrosis/mf_mrgsolve_model.R")
# Simulation and plotting run automatically
```

### 3. Run the Shiny App

```r
install.packages(c("shiny", "shinydashboard", "ggplot2", "plotly", "DT"))
shiny::runApp("myelofibrosis/mf_shiny_app.R")
```

---

## Key Clinical Trials Summary

| Trial | Drug | Population | SVR35 | TSS50 | Median OS | PMID |
|-------|------|------|-------|-------|---------|------|
| **COMFORT-I** | Ruxolitinib 20mg BID vs placebo | Int-2/High PMF | 41.9% vs 0.7% | 45.9% vs 5.3% | — | [22375971](https://pubmed.ncbi.nlm.nih.gov/22375971/) |
| **COMFORT-II** | Ruxolitinib vs BAT | Int-2/High PMF | 28% vs 0% | — | NR vs 27.9mo | [22375970](https://pubmed.ncbi.nlm.nih.gov/22375970/) |
| **JAKARTA** | Fedratinib 400mg QD | Int-2/High, Rux-naive | 36% vs 1% | 36% vs 6% | — | [26003172](https://pubmed.ncbi.nlm.nih.gov/26003172/) |
| **PERSIST-2** | Pacritinib 200mg BID | PLT <100×10⁹/L | 18% vs 3% | — | — | [29049469](https://pubmed.ncbi.nlm.nih.gov/29049469/) |
| **SIMPLIFY-1** | Momelotinib 200mg QD | Rux-naive | 26.5% vs 29% | 28.4% vs 42.2% | — | [28930484](https://pubmed.ncbi.nlm.nih.gov/28930484/) |
| **MANIFEST-2** | Pelabresib + Ruxolitinib | Rux-naive | 65.9% vs 35.2% | 52.3% vs 37.5% | NR | [39504566](https://pubmed.ncbi.nlm.nih.gov/39504566/) |

---

## References

See [`mf_references.md`](mf_references.md) for the full reference list (36 PubMed citations).

Key citations:
- Verstovsek S et al. *N Engl J Med* 2012;366:799–807 (COMFORT-I) · PMID [22375971](https://pubmed.ncbi.nlm.nih.gov/22375971/)
- Pardanani A et al. *J Clin Oncol* 2015;33:2771–2779 (JAKARTA) · PMID [26003172](https://pubmed.ncbi.nlm.nih.gov/26003172/)
- Klampfl T et al. *N Engl J Med* 2013;369:2379–2390 (CALR mutations) · PMID [24325356](https://pubmed.ncbi.nlm.nih.gov/24325356/)
- James C et al. *Nature* 2005;434:1144–1148 (JAK2 V617F discovery) · PMID [15793561](https://pubmed.ncbi.nlm.nih.gov/15793561/)

---

*This model was automatically generated by Claude Code Routine (CCR) on 2026-06-23.*  
*It was produced for educational and research purposes and must not be used directly for clinical decision-making.*
