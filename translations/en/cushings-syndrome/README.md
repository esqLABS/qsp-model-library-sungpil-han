# Cushing's Syndrome — QSP Model

> **Directory:** `cushings-syndrome/` | **Abbreviation:** CS | **Date:** 2026-06-24
> **Category:** Endocrine Disease

[![CS QSP mechanistic map](cs_qsp_model.png)](../../../cushings-syndrome/cs_qsp_model.svg)

---

## Disease Overview

**Cushing's syndrome (CS)** is an endocrine disease caused by excess endogenous or exogenous glucocorticoids.

| Item | Content |
|------|------|
| Prevalence | 10-15 per million per year (endogenous) |
| Most common cause | Cushing's disease (pituitary ACTH adenoma) — 70% |
| Other causes | Ectopic ACTH (10%), adrenal adenoma (15%), adrenal carcinoma (5%) |
| Standard treatment | Transsphenoidal surgery (pituitary adenoma) / adrenalectomy (adrenal cause) |
| Second-line treatment | Pasireotide, ketoconazole, metyrapone, osilodrostat, mifepristone |
| Core biochemical markers | UFC >150 μg/24h, LNSC >10 nmol/L, 1mg DST >1.8 μg/dL |

---

## Core Mechanisms (13 clusters, 140+ nodes)

| Cluster | Key mechanisms |
|---------|-----------|
| 1. Hypothalamus | CLOCK/BMAL1 circadian rhythm, CRH/AVP secretion, somatostatin/dopamine inhibition |
| 2. Anterior pituitary | CRHR1 → Gs → PKA → CREB → POMC → ACTH; USP8 mutation (~50% of Cushing's disease); CDK4/6-Rb proliferation |
| 3. Adrenocortical biosynthesis | MC2R → cAMP → PKA → StAR → CYP11A1 → ... → CYP11B1 → cortisol |
| 4. GR signalling | GR-α/β, HSP90, FKBP51/52, nuclear translocation, GRE/nGRE, AP-1/NF-κB tethering repression, GILZ, MKP1, SGK1 |
| 5. Metabolic complications | ↑ PEPCK/G6Pase (liver), ↓ GLUT4 (periphery), insulin resistance, ↑ visceral fat, muscle wasting, osteoporosis |
| 6. Cardiovascular/renal | RAAS hyperactivation (Ang II→AT1R→aldosterone), sodium retention, hypertension, VTE, dyslipidaemia |
| 7. Immune suppression | NF-κB/AP-1 inhibition → ↓ cytokines, lymphopenia, ↓ NK cells, ↑ infection risk |
| 8. CNS effects | Hippocampal atrophy, ↓ BDNF, NMDA excitotoxicity, ↓ serotonin, depression and cognitive decline |
| 9. Pasireotide PK/PD | 2-compartment SC; SSTR5>SSTR2; Gi→↓cAMP→ACTH inhibition (max 65%) |
| 10. Steroidogenesis inhibitors | Ketoconazole (CYP17A1+11B1), metyrapone (CYP11B1), osilodrostat (CYP11B1/B2) |
| 11. GR antagonist · cabergoline | Mifepristone (competitive GR antagonism), cabergoline (D2R→ACTH inhibition), mitotane (adrenolytic) |
| 12. Clinical diagnosis/assessment | UFC, LNSC, 1mg DST, high-dose DST, CRH stimulation test, IPSS |
| 13. Aetiological classification | Cushing's disease, ectopic ACTH, adrenal adenoma/carcinoma, PBMAH, McCune-Albright, cyclical Cushing's |

---

## mrgsolve ODE Model (21 compartments)

| Module | Compartments | Key dynamics |
|------|------|------------|
| HPA axis | CRH, ACTH_PIT, ACTH_PL | Circadian CRH synthesis + tumour ACTH secretion + GR negative feedback (Hill n=2) |
| Adrenal steroid | F_ADR, F_PL | Michaelis-Menten (ACTH→cortisol); Bliss combination drug inhibition |
| GR dynamics | GR_FREE, GR_BOUND, GR_NUC | 2-step (cytoplasm→nucleus) binding/dissociation ODE |
| Metabolism | GLUCOSE, INSULIN, VAT, MUSCLE, BMD, BP | GR-driven metabolic changes; insulin feedback |
| Clinical output | UFC_ACC, LNSC | Proportional secretion model |
| Pasireotide PK | A_PAS_C, A_PAS_P | 2-compartment; SC absorption; Emax=0.62 |
| Ketoconazole PK | A_KETO | 1-compartment; dual CYP17A1+11B1 inhibition; Emax=0.72 |
| Metyrapone PK | A_METY | 1-compartment; selective CYP11B1 inhibition; Emax=0.82 |
| Osilodrostat PK | A_OSILO | 1-compartment; high potency (EC50=0.15 μg/mL); Emax=0.85 |
| Mifepristone PK | A_MIFE | 1-compartment; high Vd (115 L); GR occupancy Emax=0.82 |

---

## Clinical Calibration Data for Treatment Scenarios

| Scenario | Drug | Clinical trial | Key result |
|---------|------|----------|----------|
| Natural course (no treatment) | — | Historical cohort | UFC >500 μg/24h, accumulating complications |
| Pasireotide 0.6 mg BID | Pasireotide | PASPORT-CUSHINGS (Colao 2012 NEJM) | UFC normalisation 22-24% (6 months) |
| Ketoconazole 400 mg BID | Ketoconazole | Castinetti 2014 Eur J Endocrinol | UFC normalisation 49%, caution for hepatotoxicity |
| Osilodrostat 5 mg BID | Osilodrostat | LINC 3/4 (Pivonello 2020 Lancet DE) | UFC normalisation 86% (maintenance phase) |
| Mifepristone 600 mg QD | Mifepristone | SEISMIC (Fleseriu 2012 JCEM) | Glycaemic/BP clinical response 87% (GR antagonism) |
| Post-surgical remission | — | Meta-analysis | Recurrence rate ~20% (within 5 years) |

---

## QSP Model Component Files

| Component | File | Specification |
|---------|------|-----|
| 🗺️ Mechanistic map (DOT) | [`cs_qsp_model.dot`](cs_qsp_model.dot) | **140+ nodes, 13 clusters** |
| 🖼️ SVG vector image | [`cs_qsp_model.svg`](cs_qsp_model.svg) | Zoomable, high resolution |
| 🖼️ PNG raster image | [`cs_qsp_model.png`](cs_qsp_model.png) | 150 dpi |
| ⚙️ mrgsolve ODE model | [`cs_mrgsolve_model.R`](cs_mrgsolve_model.R) | **21-compartment ODE, 6 treatment scenarios** |
| 📊 Shiny app | [`cs_shiny_app.R`](cs_shiny_app.R) | **8-tab dashboard** |
| 📚 References | [`cs_references.md`](cs_references.md) | **55 PubMed citations (10 sections)** |

---

## Shiny App Tab Layout (8 tabs)

| Tab | Content |
|----|------|
| 1. Patient profile | Baseline cortisol/ACTH/glucose settings, aetiology selection, diagnostic criteria |
| 2. HPA axis/PK dynamics | Circadian rhythm, ACTH dynamics, cortisol time series, drug PK |
| 3. Steroid biosynthesis | Enzyme-inhibition dynamics, pathway overview, cortisol synthesis inhibition |
| 4. Clinical indices | UFC, LNSC, remission-criteria assessment, dexamethasone suppression test |
| 5. Scenario comparison | Side-by-side comparison of 5 treatments (cortisol/ACTH/glucose/summary table) |
| 6. Biomarker panel | GR nuclear activation, HPA feedback, dexamethasone test simulation |
| 7. Metabolic complications | Glucose/insulin, body composition (VAT+muscle), BMD, blood pressure |
| 8. Virtual population analysis | N=10-500 virtual patients, response rate, UFC distribution |

---

## Diagnostic Biomarker Reference Values

| Biomarker | Normal | Cushing's syndrome | Unit |
|-----------|------|------------|------|
| UFC 24h | < 50 | > 150 (often >500) | μg/24h |
| LNSC (midnight salivary) | < 4 | > 10 | nmol/L |
| 1mg DST cortisol | < 1.8 | > 1.8 (not suppressed) | μg/dL |
| Plasma ACTH | 10–46 | ↑ (pituitary/ectopic) / ↓ (adrenal) | pg/mL |
| Fasting glucose | < 5.6 | 5.6–11.1+ | mmol/L |
| BMD T-score | > -1.0 | < -1.0 (often < -2.5) | |
| Systolic blood pressure | < 140 | > 140 (often 150-180) | mmHg |

---

## References (Key 5)

1. Colao et al. (2012) *N Engl J Med* — pasireotide phase III
2. Pivonello et al. (2020) *Lancet Diabetes Endocrinol* — osilodrostat LINC 3
3. Nieman et al. (2018) *J Clin Endocrinol Metab* — treatment guidelines
4. Fleseriu et al. (2012) *J Clin Endocrinol Metab* — mifepristone SEISMIC
5. Miller & Auchus (2011) *Endocr Rev* — steroidogenesis mechanisms

---

*Model created: 2026-06-24 | Claude Code Routine | QSP Library [pipetcpt/qsp](https://github.com/pipetcpt/qsp)*
