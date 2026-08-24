# Hereditary Hemochromatosis (HH) QSP Model

> **Directory:** `hereditary-hemochromatosis/` | **Abbreviation:** HH | **Date:** 2026-06-25

[![HH QSP Mechanistic Map](../../../hereditary-hemochromatosis/hh_qsp_model.png)](../../../hereditary-hemochromatosis/hh_qsp_model.svg)

---

## Disease Overview

**Hereditary hemochromatosis (HH)** is an autosomal recessive inherited metabolic
disease in which **HFE gene mutations (C282Y/H63D)** reduce hepcidin production,
resulting in excessively increased intestinal iron absorption and consequent iron
accumulation in the liver, heart, pancreas, joints, and pituitary gland.

| Feature | Content |
|------|------|
| **Prevalence** | 1/200–400 in people of northern European descent |
| **Main gene** | HFE (chr 6p21.3) |
| **Key mutations** | C282Y (~85%), H63D (~15%), S65C (rare) |
| **Penetrance** | 15–50% (biochemical); ~10% (clinical expression) |
| **Total body iron** | Normal 3–5 g → HH 20–40 g |
| **Diagnostic criteria** | TSAT >45% + elevated serum ferritin → HFE genotyping |
| **Treatment** | Phlebotomy — standard; iron chelators (DFO/DFX/DFP) |

---

## Mechanistic Map Clusters (11, 130+ nodes)

| Cluster | Key components |
|---------|------------|
| ① Genetic/molecular drivers | HFE C282Y/H63D, HFE-TfR1/TfR2 complex, HJV, TMPRSS6, the four HH types |
| ② Hepcidin regulation (BMP-SMAD) | BMP2/4/6/9, BMPR-I/II, SMAD1/5/8, SMAD4, IL-6 STAT3, ERFE, hepcidin deficiency |
| ③ Iron absorption/transport | DMT1, DCYTB, ferroportin (FPN), hephaestin, Tf-TBI, NTBI, LPI, TSAT |
| ④ Cellular iron storage | Ferritin (H/L chain), hemosiderin, IRP1/2, IRE, LIP, mitochondrial iron, ferritinophagy |
| ⑤ Erythropoiesis/circulation | EPO, BFU-E/CFU-E, reticulocytes, mature RBCs, RES macrophages, ERFE |
| ⑥ Liver damage | TfR1/ZIP14 hepatocyte uptake, LIC, ROS/Fenton reaction, HSC activation, fibrosis, cirrhosis, HCC |
| ⑦ Extrahepatic organ damage | Cardiac T2* MRI, cardiomyopathy/arrhythmia, pancreatic β-cells, bronze diabetes, chondrocalcinosis, hypopituitarism |
| ⑧ Chelation therapy PK | DFO (2-compartment SC), DFX (oral once daily), DFP (oral TID), hepcidin agonists (experimental) |
| ⑨ Phlebotomy and clinical management | Induction/maintenance phlebotomy, dietary restriction, target ferritin 50–100 ng/mL |
| ⑩ Clinical diagnosis | TSAT, ferritin, HFE genotype, MRI R2*, FibroScan, echocardiography, HbA1c |
| ⑪ Treatment PD/outcomes | Body iron burden, iron removal per agent, fibrosis reversal, reduced HCC risk |

---

## ODE Compartments (20)

| # | Compartment | Unit | Biological meaning |
|---|------|------|------------|
| 1 | `TBI` | mg | Plasma transferrin-bound iron |
| 2 | `NTBI_pool` | mg | Non-transferrin-bound iron (toxic) |
| 3 | `FPN` | AU | Relative ferroportin activity |
| 4 | `HEPC` | ng/mL | Plasma hepcidin |
| 5 | `FERRITIN` | ng/mL | Serum ferritin (biomarker) |
| 6 | `LIVER_Fe` | mg | Total liver iron accumulation |
| 7 | `HEART_Fe` | mg | Cardiac iron accumulation |
| 8 | `PANCR_Fe` | mg | Pancreatic iron accumulation |
| 9 | `RBC_Fe` | mg | Hemoglobin-bound iron in red blood cells |
| 10 | `MACRO_Fe` | mg | RES macrophage iron (recycling) |
| 11 | `LIV_FIB` | 0–4 | Liver fibrosis score (Metavir) |
| 12 | `BCELL` | 0–1 | Pancreatic β-cell function fraction |
| 13 | `EF` | % | Cardiac ejection fraction |
| 14 | `DFO_C` | mg | Deferoxamine plasma concentration |
| 15 | `DFO_FE` | mg | Ferrioxamine complex |
| 16 | `DFX_C` | mg | Deferasirox plasma |
| 17 | `DFX_FE` | mg | DFX-Fe complex |
| 18 | `DFP_C` | mg | Deferiprone plasma |
| 19 | `DFP_FE` | mg | DFP-Fe complex |
| 20 | `HBA1C` | % | HbA1c (marker of bronze diabetes) |

---

## Treatment Scenarios (6)

| Scenario | Regimen | Clinical calibration basis |
|---------|------|-------------|
| S1 | **No treatment**, natural course of C282Y homozygotes | Allen 2008 NEJM; 28% of men have iron overload |
| S2 | **Phlebotomy** 500 mL q2wk (induction) → q3months (maintenance) | Niederau 1996 Gastroenterology; AASLD 2011 |
| S3 | **Deferoxamine (DFO)** 40 mg/kg/day SC | Hoffbrand 2003 Blood; ~70 mg Fe removed/24h |
| S4 | **Deferasirox (DFX)** 20 mg/kg/day oral | Cappellini 2014 Blood; ~30 mg/kg/day balance |
| S5 | **Deferiprone (DFP)** 75 mg/kg/day TID (cardiac-focused) | Anderson 2002 Lancet; T2* +27% vs. DFO +13% |
| S6 | **DFP+DFO combination** (shuttle chelation) | Tanner 2007 Circulation; most effective for cardiac iron |

---

## Clinical Calibration

| Metric | Calibration target | Basis |
|------|-----------|------|
| Untreated ferritin rise | +50–100 ng/mL/year | Bacon 2011 Hepatology |
| Phlebotomy ferritin decrease | ~50 ng/mL/session | Phatak 2010 Am J Hematol |
| LIC normalization | 2–3 years of regular phlebotomy | Brissot 2018 Nat Rev Dis Primers |
| DFP cardiac T2* improvement | +27% at 1 year | Anderson 2002 Lancet |
| DFP+DFO ferritin decrease | −1019 ng/mL/year | Tanner 2007 Circulation |
| HCC risk reduction | 30× → ~5× with phlebotomy | Elmberg 2003 Gastroenterology |

---

## Shiny App Tab Structure (6)

| Tab | Content |
|----|------|
| ① Patient Profile | Genotype, age, baseline iron status, comparison with normal ranges |
| ② Drug PK | Chelator PK profiles and parameter comparison table |
| ③ Iron Kinetics | TBI · NTBI · liver iron · RBC iron dynamics; ferritin · TSAT trends |
| ④ Organ Iron Load | LIC/fibrosis, cardiac T2*/EF, pancreatic β-cell/HbA1c, hepcidin |
| ⑤ Scenario Compare | Comparison of ferritin/LIC/T2*/fibrosis across 6 scenarios |
| ⑥ Biomarkers | Yearly biomarker table, HbA1c trend, iron balance, risk radar chart |

---

## Usage

```bash
# 1. Render the mechanistic map
sfdp -Tsvg hh_qsp_model.dot -o hh_qsp_model.svg
sfdp -Tpng -Gdpi=150 hh_qsp_model.dot -o hh_qsp_model.png
# or convert SVG → PNG
rsvg-convert -d 150 -p 150 hh_qsp_model.svg -o hh_qsp_model.png
```

```r
# 2. Run the mrgsolve model
install.packages(c("mrgsolve","dplyr","ggplot2","gridExtra"))
source("hh_mrgsolve_model.R")

# 3. Run the Shiny dashboard
install.packages(c("shiny","shinydashboard","plotly","DT"))
shiny::runApp("hh_shiny_app.R")
```

---

## File List

| File | Description |
|------|------|
| `hh_qsp_model.dot` | Graphviz mechanistic map source (130+ nodes, 11 clusters) |
| `hh_qsp_model.svg` | Vector-format mechanistic map |
| `hh_qsp_model.png` | Raster-format mechanistic map (150 dpi) |
| `hh_mrgsolve_model.R` | mrgsolve ODE model (20 compartments, 6 scenarios) |
| `hh_shiny_app.R` | Shiny dashboard (6 tabs) |
| `hh_references.md` | References, 55 articles (15 sections) |
| `README.md` | This file |

---

*Generated by Claude Code Routine — 2026-06-25*
