# Pheochromocytoma/Paraganglioma (PPGL) — QSP Model

> **Directory:** `pheochromocytoma/` | **Abbreviation:** PPGL | **Date:** 2026-06-25

[![PPGL QSP Mechanistic Map](ppgl_qsp_model.png)](ppgl_qsp_model.svg)

---

## Disease Overview

Pheochromocytoma is a catecholamine-secreting neuroendocrine tumor arising from
**chromaffin cells of the adrenal medulla**; when it arises from extra-adrenal
sympathetic ganglia it is classified as **paraganglioma**. The two tumor types are
collectively referred to as **PPGL**.

| Feature | Content |
|------|------|
| **Prevalence** | 2–8 per 100,000 population; 5% of incidentally discovered adrenal tumors |
| **Heritability** | 40% of all cases are germline variants (SDH, VHL, RET, NF1, MAX, etc.) |
| **Malignancy rate** | 10–17% overall; up to 40–80% with SDHB variants |
| **Main symptoms** | Headache-sweating-palpitation triad (WHO criteria); hypotension·hypertension·crisis episodes |
| **Diagnostic biomarkers** | Plasma free metanephrines (NMN/MN), 24h urinary catecholamines, chromogranin-A |
| **Surgical cure rate** | Localized PPGL: 80–95% biochemical cure after R0 resection |

### Key pathophysiology

1. **Catecholamine oversecretion**: excess NE/EPI secreted via the TH (tyrosine
   hydroxylase)-dependent biosynthetic pathway → vasoconstriction through α₁ receptors
   → causes **hypertensive crisis**
2. **Genetic clustering**: SDHB/C/D variants (pseudohypoxia cluster) vs. RET/NF1/TMEM127
   variants (kinase signaling cluster) — differing clinical phenotype and malignancy
3. **Cardiovascular toxicity**: chronic catecholamine excess → catecholamine
   cardiomyopathy, takotsubo-like cardiomyopathy, arrhythmia

---

## Mechanistic Map

| File | Format |
|------|------|
| `ppgl_qsp_model.dot` | Graphviz source |
| `ppgl_qsp_model.svg` | Vector (high resolution) |
| `ppgl_qsp_model.png` | Raster (150 dpi) |

### 11 clusters (130+ nodes)

| Cluster | Key components |
|---------|------------|
| ① Genetic/molecular drivers | SDHB/C/D/A · VHL · RET · NF1 · MAX · EPAS1/HIF-2α · pseudohypoxia · kinase cluster · malignancy |
| ② Tumor biology | Chromaffin cells · tumor proliferation · apoptosis · HIF-1α · VEGF-A · angiogenesis · CgA · NSE · metastasis |
| ③ Catecholamine biosynthesis | Tyrosine→TH→DOPA→AADC→dopamine→DBH→NE→PNMT→EPI; COMT·MAO metabolism |
| ④ Storage/secretion | Chromaffin granules · VMAT2 · Ca²⁺ influx · exocytosis · SNARE · NET reuptake |
| ⑤ Adrenergic receptor signaling | α₁/α₂/β₁/β₂/β₃-AR · Gq/Gi/Gs · PLC→IP3/DAG→PKC · cAMP→PKA · MAPK |
| ⑥ Cardiovascular effects | SBP · DBP · MAP · HR · CO · SVR · LVEF · catecholamine cardiomyopathy · hypertensive crisis · arrhythmia |
| ⑦ Metabolic effects | Hepatic glycogenolysis · glucagon↑ · insulin suppression · hyperglycemia · FFA · BAT thermogenesis · BMR↑ · weight loss |
| ⑧ α-blocker PK | Phenoxybenzamine (irreversible) 2-compartment · doxazosin (selective α₁) 1-compartment PK |
| ⑨ Systemic therapy PK | Sunitinib 2-compartment · ¹³¹I-MIBG · ¹⁷⁷Lu-DOTATATE · CVD chemotherapy |
| ⑩ Pharmacodynamics/biomarkers | α₁/β receptor occupancy · TH inhibition · VEGFR inhibition · plasma NMN/MN · CgA · RECIST response |
| ⑪ Surgical/perioperative management | Preoperative alpha/beta blockade · metyrosine · IV fluids · laparoscopic/open surgery · intraoperative crisis · postoperative hypotension |

---

## mrgsolve ODE Model

**File:** `ppgl_mrgsolve_model.R`

### ODE Compartments (20)

| # | Compartment | Biological meaning |
|---|------|------------|
| 1–3 | `PHE_gut`, `PHE_C`, `PHE_P` | Phenoxybenzamine GI · central · peripheral compartments |
| 4 | `DOX_C` | Doxazosin central compartment |
| 5 | `MET_C` | Metyrosine central compartment |
| 6 | `BB_C` | Beta-blocker (propranolol) central |
| 7–8 | `SUNIT_C`, `SUNIT_P` | Sunitinib central · peripheral (malignant PPGL) |
| 9 | `TH_act` | TH enzyme activity (inhibited by metyrosine) |
| 10 | `NE_store` | Chromaffin granule NE storage pool |
| 11 | `NE_plasma` | Plasma NE |
| 12 | `EPI_plasma` | Plasma EPI |
| 13 | `TUMvol` | Tumor volume (mL) |
| 14 | `VEGF_tum` | Plasma VEGF (pg/mL) |
| 15 | `SBP` | Systolic blood pressure (mmHg) |
| 16 | `DBP` | Diastolic blood pressure (mmHg) |
| 17 | `HR` | Heart rate (bpm) |
| 18 | `GLU` | Plasma glucose (mmol/L) |
| 19 | `FFA` | Free fatty acids (mmol/L) |
| 20 | `CgA_plasma` | Plasma chromogranin-A (ng/mL) |

### Treatment Scenarios (6)

| Scenario | Regimen | Clinical basis |
|---------|------|---------|
| S0 | No treatment (natural course) | — |
| S1 | **Phenoxybenzamine** 60 mg/d × 14 days → surgery | Kinney 2002 J Cardiothorac Vasc Anesth |
| S2 | **Doxazosin** 16 mg/d × 14 days → surgery | Shao 2016 World J Surg (meta-analysis) |
| S3 | **PHE + metyrosine 2 g/d + propranolol** triple combination | Steinsapir 1997 Arch Intern Med |
| S4 | **Sunitinib** 37.5 mg/d (metastatic malignant PPGL) | Niemeijer 2014 J Clin Endocrinol Metab |
| S5 | **Metyrosine alone** (inoperable, symptom control) | Engelman 1968 NEJM |

### Clinical Calibration

| Clinical trial/study | Calibration target |
|-------------|---------|
| Lentschener 2009 Hypertension | Equivalence of phenoxybenzamine vs. doxazosin in perioperative blood-pressure control |
| Steinsapir 1997 | Metyrosine TH inhibition 40–80%, reduced catecholamine synthesis |
| Niemeijer 2014 | Sunitinib in malignant PPGL: ORR 25%, tumor stabilization |
| Averbuch 1988 | CVD chemotherapy ORR 37%, clinical response rate 79% |

---

## Shiny Dashboard

**File:** `ppgl_shiny_app.R`

| Tab | Key content |
|----|---------|
| 1. Patient Profile | PPGL type · tumor volume · baseline NE/EPI · SBP/HR · CgA · genetic variants |
| 2. Drug PK | PHE/DOX/MET/BB/sunitinib blood-concentration simulation |
| 3. Catecholamines | Plasma NE · EPI · NMN · MN dynamics; TH activity; CgA tumor marker |
| 4. Cardiovascular | SBP · DBP · HR · α-blockade ratio · metabolic effects (glucose · FFA) |
| 5. Tumor (malignant) | Tumor volume · VEGF · RECIST response assessment; treatment-option comparison table |
| 6. Scenario Comparison | Parallel comparison of 6 scenarios; summary table of key metrics |

---

## Usage

```bash
# Render the mechanistic map
dot -Tsvg ppgl_qsp_model.dot -o ppgl_qsp_model.svg
dot -Tpng -Gdpi=150 ppgl_qsp_model.dot -o ppgl_qsp_model.png
```

```r
# Run the mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2"))
source("ppgl_mrgsolve_model.R")

# Run the Shiny dashboard
install.packages(c("shiny","shinydashboard","tidyr"))
shiny::runApp("ppgl_shiny_app.R")
```

---

## References

**File:** `ppgl_references.md` — 45 PubMed citations (12 sections)

| Section | Count |
|------|------|
| Landmark reviews · guidelines | 5 |
| Epidemiology · genetics | 4 |
| Catecholamine biosynthesis · metabolism | 5 |
| Preoperative management · alpha blockade | 6 |
| Metyrosine PK/PD | 3 |
| Systemic therapy for malignant PPGL | 6 |
| Cardiovascular · hemodynamic effects | 4 |
| Biochemical diagnosis · biomarkers | 3 |
| Imaging · localization | 2 |
| Molecular pathophysiology | 3 |
| Sunitinib PK modeling | 2 |
| QSP/PK-PD modeling | 2 |
