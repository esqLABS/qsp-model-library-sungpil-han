# Non-Small Cell Lung Cancer (NSCLC) — QSP Model

> A quantitative systems pharmacology (QSP) model of **non-small cell lung cancer
> (NSCLC)**  
> Driver mutations · Targeted therapies · Immune checkpoint inhibitors · Chemotherapy · Resistance mechanisms

[![NSCLC QSP Model](nsclc_qsp_model.png)](nsclc_qsp_model.svg)

---

## Disease Overview

Non-small cell lung cancer (NSCLC) accounts for about 85% of lung cancers and is the
leading cause of cancer death worldwide.

**Key Pathophysiology:**
- **Driver mutations:** EGFR (exon19 del/L858R/T790M), KRAS G12C, ALK/ROS1 fusion, BRAF V600E, MET exon14, RET, NTRK
- **Oncogenic signaling:** RAS-RAF-MEK-ERK, PI3K-AKT-mTOR, JAK-STAT3 → unrestrained proliferation
- **Immune evasion:** PD-L1 overexpression → suppression of CD8+ T-cell function; Treg/MDSC accumulation; IDO1/TGF-β
- **Tumor microenvironment (TME):** HIF1α → VEGF → tumor angiogenesis; EMT → metastasis
- **Drug resistance:** T790M/C797S (EGFR), bypass signaling (MET amp, HER2 amp), SCLC transformation

---

## Model Structure

### Mechanistic Map
| Item | Value |
|------|-----|
| Node count | 231 |
| Cluster count | 10 |
| Edge count | 301 |
| Render engine | sfdp |

**10 clusters:**
1. **Oncogene drivers** — EGFR/KRAS/ALK/ROS1/BRAF/MET/HER2/RET/NTRK/NRG1/PIK3CA/STK11/KEAP1/TP53
2. **Oncogenic signaling** — RAS/RAF/MEK/ERK, PI3K/AKT/mTOR, JAK/STAT3, SRC/FAK
3. **Tumor biology** — HIF1α-VEGF-angiogenesis, EMT, CDK4/6-RB1 cell cycle, BCL2/BAX apoptosis, TERT
4. **Immune checkpoint** — PD-L1/PD-1, CTLA4/CD28/B7, CD8+/NK/Treg/TAM, IDO1, LAG3/TIM3/TIGIT
5. **Drug PK/PD** — Osimertinib, Alectinib, Sotorasib, Pembrolizumab, Atezolizumab, Cisplatin, Pemetrexed
6. **Resistance mechanisms** — T790M, C797S, MET/HER2 amp, AXL/FGFR bypass, SCLC transformation
7. **Biomarkers** — PD-L1 TPS, TMB, ctDNA, ALK FISH, CEA, CYFRA 21-1, OS/PFS/ORR
8. **Tumor dynamics** — Sensitive/resistant cell pools, Gompertz growth, RECIST response categories
9. **Toxicity** — irAE (pneumonitis/colitis/nephritis), TKI rash/ILD, chemo myelosuppression
10. **Clinical endpoints** — OS, PFS, ORR, DCR, RECIST CR/PR/SD/PD, Stage IIIA-IVB, CNS/bone mets

### mrgsolve ODE Model
| Item | Value |
|------|-----|
| Number of ODE compartments | 19 |
| Drug PK compartments | 9 (Osimertinib 2 + Alectinib 2 + Sotorasib 2 + Pembrolizumab 1 + Cisplatin 1 + Pemetrexed 1) |
| Disease PD compartments | 10 (tumor sensitive/resistant/total, CD8+ Teff, PD1 occupancy, Treg, CEA, ctDNA, ANC, PD-L1) |
| Treatment scenarios | 7 |
| Simulation duration | 24 months |

**7 treatment scenarios:**
| Scenario | Indication | Supporting trial |
|----------|--------|-------------|
| S1 | No treatment (natural progression) | — |
| S2 | Osimertinib 80mg/day | FLAURA (PFS 18.9 mo) |
| S3 | Alectinib 600mg BID | ALEX (PFS 34.8 mo) |
| S4 | Carboplatin + Pemetrexed Q3W | ECOG 1594 |
| S5 | Pembrolizumab 200mg Q3W (PD-L1≥50%) | KEYNOTE-024 (PFS 10.3 mo) |
| S6 | Pembrolizumab + Carboplatin + Pemetrexed | KEYNOTE-189 (PFS 9.0 mo) |
| S7 | Sotorasib 960mg/day (KRAS G12C+) | CodeBreaK100 (ORR 37%) |

### Shiny Dashboard (6 tabs)
| Tab | Content |
|----|------|
| 1. Patient Profile | Age/ECOG/Stage/histology/molecular profile/PD-L1/TMB; treatment indication matrix |
| 2. Pharmacokinetics (PK) | Multi-cycle concentration-time curves per drug; Cmax/AUC/t½ table; renal/hepatic dose adjustment |
| 3. Tumor Response | Sensitive/resistant cell trajectories; RECIST waterfall; spider plot |
| 4. Biomarkers | CEA, ctDNA allele fraction, ANC nadir, PFS KM curve; resistance prediction |
| 5. Scenario Comparison | Multi-way comparison of 7 treatments; ORR/DCR/PFS table; forest plot + clinical trial benchmark |
| 6. Toxicity/Safety | Hematologic toxicity simulation; irAE probability; organ-specific traffic-light indicators; dose modification recommendations |

---

## Files

| File | Description |
|------|------|
| [`nsclc_qsp_model.dot`](nsclc_qsp_model.dot) | Graphviz DOT mechanistic map source |
| [`nsclc_qsp_model.svg`](nsclc_qsp_model.svg) | Vector graphic |
| [`nsclc_qsp_model.png`](nsclc_qsp_model.png) | Raster image |
| [`nsclc_mrgsolve_model.R`](nsclc_mrgsolve_model.R) | mrgsolve ODE model + simulation code |
| [`nsclc_shiny_app.R`](nsclc_shiny_app.R) | Shiny interactive dashboard |
| [`nsclc_references.md`](nsclc_references.md) | 78 references (classified by section) |

---

## How to Run

```bash
# Render the mechanistic map
sfdp -Tsvg nsclc_qsp_model.dot -o nsclc_qsp_model.svg
sfdp -Tpng -Gdpi=96 nsclc_qsp_model.dot -o nsclc_qsp_model.png
```

```r
# mrgsolve simulation
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "survminer"))
library(mrgsolve)
source("nsclc_mrgsolve_model.R")

# Shiny dashboard
install.packages(c("shiny", "shinydashboard", "ggplot2", "dplyr", "plotly", "DT"))
shiny::runApp("nsclc_shiny_app.R")
```

---

## Key Parameters

| Parameter | Value | Unit | Basis |
|----------|-----|------|------|
| Baseline tumor growth rate (kg) | 0.002 | /day | NSCLC TGI literature |
| Osimertinib t½ | 48 | h | FLAURA PK |
| Alectinib t½ | 33 | h | ALEX PK |
| Pembrolizumab t½ | 27 | days | KEYNOTE biologic PK |
| Resistance emergence rate (kr) | 0.0001 | /day | Clinical resistance-emergence model |
| Osimertinib Emax | 0.92 | — | Calibrated to FLAURA median PFS |
| PD-L1 TPS threshold | 50% | % | KEYNOTE-024 selection criterion |

---

## Key References

- Mok TS et al. (2009) IPASS: gefitinib in EGFR+ NSCLC. *NEJM*. PMID: 19692680
- Soria JC et al. (2018) FLAURA: osimertinib vs 1st-gen TKI. *NEJM*. PMID: 29151359
- Peters S et al. (2017) ALEX: alectinib vs crizotinib. *NEJM*. PMID: 28586279
- Reck M et al. (2016) KEYNOTE-024: pembrolizumab vs chemo (PD-L1≥50%). *NEJM*. PMID: 27718847
- Gandhi L et al. (2018) KEYNOTE-189: pem+chemo. *NEJM*. PMID: 29658856
- Skoulidis F et al. (2021) CodeBreaK100: sotorasib. *NEJM*. PMID: 34096690

See [`nsclc_references.md`](nsclc_references.md) for the full list of 78 references.

---

*Generated by Claude Code Routine (CCR) — 2026-06-23*
