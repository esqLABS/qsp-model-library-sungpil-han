# Fabry Disease QSP Model (Fabry Disease Quantitative Systems Pharmacology Model)

> **Directory:** `fabry-disease/` | **Abbreviation:** FBR | **Date:** 2026-06-24
> **Category:** Rare · Lysosomal Storage Disease (X-linked lysosomal storage disease)

[![FBR QSP mechanistic map](fbr_qsp_model.png)](../../../fabry-disease/fbr_qsp_model.svg)

---

## Disease Overview

**Fabry disease** is an X-linked lysosomal storage disease caused by mutations in the *GLA* gene (Xq22.1). Deficiency or reduced function of the enzyme α-galactosidase A (α-Gal A) causes the glycosphingolipid **globotriaosylceramide (Gb3/GL-3)** and its deacylated derivative **lyso-Gb3** to accumulate progressively throughout the body, including the kidneys, heart, brain, skin, and peripheral nervous system.

| Characteristic | Detail |
|------|------|
| Prevalence | 1:40,000–1:117,000 (classic form) / 1:3,000–1:10,000 (including late-onset forms) |
| Inheritance | X-linked — hemizygous males severely affected, heterozygous females variably affected |
| Major phenotypes | Classic form (residual enzyme <1%), late-onset cardiac/renal forms (residual 1–30%) |
| Diagnostic delay | Average 10–20 years (symptom onset to diagnosis) |
| Causes of death | Renal failure, cardiac events (sudden death, heart failure), stroke |

---

## Core Mechanisms (14 clusters)

| Cluster | Key Content |
|---------|----------|
| 1. Genetic basis | GLA gene mutations (Xq22.1): missense ~60%, nonsense/frameshift, splicing; amenable mutations (~40%, migalastat-eligible) |
| 2. α-Gal A enzyme biology | ER synthesis → Golgi M6P phosphorylation → M6P receptor → lysosomal delivery → Gb3 hydrolysis at pH 4.5–5.0 |
| 3. Glycosphingolipid metabolism | Ceramide → GCS → GlcCer → LacCer → A4GalT → Gb3 biosynthesis; Gb3 deacylation → lyso-Gb3 (toxic signalling molecule) |
| 4. Renal pathology | Podocyte Gb3 accumulation → foot-process effacement → proteinuria (increased UPCR) → TGF-β fibrosis → FSGS → declining eGFR → ESRD |
| 5. Cardiac pathology | Cardiomyocyte Gb3 → left ventricular hypertrophy (increased LVMi) → diastolic dysfunction → myocardial fibrosis (LGE) → arrhythmia → sudden death |
| 6. Neurological pathology | DRG Gb3 → small-fibre neuropathy → neuropathic pain (BPI-SF) + CNS vascular endothelium → white-matter lesions → TIA/stroke |
| 7. Other organs | Angiokeratoma, anhidrosis, corneal verticillata, gastrointestinal dysmotility, sensorineural hearing loss |
| 8. Inflammatory cascade | Lyso-Gb3 → TLR4/NF-κB → IL-6/TNF-α → NLRP3 inflammasome → decreased eNOS → endothelial activation |
| 9. ERT PK/PD | Agalsidase beta (1 mg/kg Q2W, t½ ~45 min), alfa (0.2 mg/kg Q2W), pegunigalsidase alfa (1 mg/kg Q4W, t½ ~80 h) → M6P → lysosomal Gb3 degradation (Emax ~80%) |
| 10. Chaperone therapy (migalastat) | 150 mg PO QOD; stabilises misfolded α-Gal A protein; amenable mutations only; ATTRACT trial (non-inferior vs ERT) |
| 11. Substrate reduction therapy (SRT) | Lucerastat 1000 mg TID (GCS inhibition, decreased Gb3 precursor); MODIFY trial (BPI-SF −1.5 points); venglustat (CNS-penetrant) |
| 12. Biomarkers | Plasma lyso-Gb3 (μg/L, most sensitive), urinary Gb3 (nmol/mg Cr), α-Gal A activity (nmol/h/mg), DBS newborn screening |
| 13. Clinical endpoints | eGFR slope, UPCR, LVMi, BPI-SF pain, EQ-5D QoL, MSSI severity |
| 14. Natural history | Classic males (childhood onset) vs. late-onset forms vs. female carriers; Fabry Registry data; diagnostic delay averaging 10–20 years |

---

## mrgsolve ODE Model (22 compartments)

| Module | Compartments | Key Dynamics |
|------|------|------------|
| Agalsidase beta PK | A_AGAB_C, A_AGAB_P, A_AGAB_LYS | 2-compartment + lysosomal delivery; CL=0.42 L/h; t½ ~45 min |
| Agalsidase alfa PK | A_AGAA_C, A_AGAA_P, A_AGAA_LYS | 2-compartment; CL=0.55 L/h |
| Migalastat PK | A_MIG_GUT, A_MIG_C | 1-compartment oral; ka=0.82/h, F=75%, t½~3.5h |
| Lucerastat PK | A_LUC_GUT, A_LUC_C | 1-compartment oral; IC50_GCS=0.18 μg/mL, Emax=42% |
| α-Gal A enzyme | E_GalA | ERT (Emax 70 nmol/h/mg) + migalastat (Emax 6) + baseline residual activity |
| Glycosphingolipid | GB3_PLM, GB3_KID, GB3_HRT, LGB3_PLM | Synthesis minus enzyme-dependent degradation; SRT reduces upstream supply |
| Inflammation | INFLAM | lyso-Gb3-driven production k_in minus clearance k_out ODE |
| Renal function | eGFR, UPCR | Decline dependent on cumulative Gb3_KID; protected by ERT |
| Cardiac function | LVMi | Hypertrophy driven by Gb3_HRT; reversed by ERT |
| Neuropathic pain | PAIN | BPI-SF; decreases with enzyme activity |

---

## Clinical Evidence for the 6 Treatment Scenarios

| Scenario | Treatment | Clinical Trial | Key Result |
|---------|--------|---------|---------|
| S1: Natural history | None | Mehta 2009 Eur J Clin Invest | eGFR −3–12/yr, LVMi increases annually, lyso-Gb3 30–80 μg/L |
| S2: Agalsidase beta | 1 mg/kg IV Q2W | FABRY-001 (Eng 2001 NEJM); Banikazemi 2007 AIM | 61% reduction in composite renal/cardiac/cerebral events |
| S3: Agalsidase alfa | 0.2 mg/kg IV Q2W | Schiffmann 2001 Ann Intern Med | Improved neuropathic pain, stabilised renal function |
| S4: Migalastat | 150 mg PO QOD | ATTRACT (Germain 2016 NEJM); Hughes 2017 Lancet | Non-inferior to ERT; eGFR slope −0.3 mL/min/yr |
| S5: Pegunigalsidase alfa | 1 mg/kg IV Q4W | BRIGHT (Schiffmann 2021 JAMA Intern Med) | Stable eGFR, lyso-Gb3 −50%, Q4W convenience |
| S6: ERT + lucerastat | Agalsidase beta + 1000 mg TID | MODIFY (Lenders 2022 Lancet DE) | Additional Gb3 reduction, BPI-SF pain improvement of −1.5 points |

---

## QSP Model Files

| Deliverable | File | Specification |
|--------|------|------|
| 🗺️ Mechanistic map | [`fbr_qsp_model.dot`](fbr_qsp_model.dot) · [`fbr_qsp_model.svg`](fbr_qsp_model.svg) · [`fbr_qsp_model.png`](fbr_qsp_model.png) | **138 nodes, 14 clusters** |
| ⚙️ mrgsolve ODE | [`fbr_mrgsolve_model.R`](fbr_mrgsolve_model.R) | **22-compartment ODE**, **6 treatment scenarios** |
| 📊 Shiny app | [`fbr_shiny_app.R`](fbr_shiny_app.R) | **8 tabs** (patient profile · PK/enzyme · Gb3 dynamics · renal · cardiac · scenario comparison · biomarkers · virtual population) |
| 📚 References | [`fbr_references.md`](fbr_references.md) | **60 PubMed citations** (14 sections) |

---

## Usage

```bash
# 1. Render the mechanistic map (Graphviz)
dot -Tsvg fbr_qsp_model.dot -o fbr_qsp_model.svg
dot -Tpng -Gdpi=150 fbr_qsp_model.dot -o fbr_qsp_model.png
```

```r
# 2. Run the mrgsolve ODE model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "tidyr", "patchwork"))
source("fbr_mrgsolve_model.R")

# 3. Launch the Shiny dashboard
shiny::runApp("fbr_shiny_app.R")
```

---

## Key Parameter Summary

| Drug | Regimen | t½ | Mechanism of Action | Clinical Efficacy |
|------|------|-----|---------|----------|
| Agalsidase beta (Fabrazyme) | 1 mg/kg IV Q2W | ~45 min | M6P → lysosomal Gb3 degradation (Emax ~80%) | 61% reduction in composite events |
| Agalsidase alfa (Replagal) | 0.2 mg/kg IV Q2W | ~45–110 min | M6P → lysosomal Gb3 degradation (Emax ~70%) | Improved neuropathic pain and renal function |
| Pegunigalsidase alfa (Elfabrio) | 1 mg/kg IV Q4W | ~80 hours | PEGylation-extended t½ (Emax ~85%) | Stable eGFR, Q4W convenience |
| Migalastat (Galafold) | 150 mg PO QOD | ~3.5 hours | α-Gal A chaperone (EC50 ~0.25 μg/mL) | Non-inferior to ERT in amenable mutations |
| Lucerastat (combination) | 1000 mg PO TID | ~8 hours | GCS inhibition, IC50 ~0.18 μg/mL (Emax 42%) | Decreased Gb3, BPI-SF −1.5 points |

---

*Claude Code Routine (CCR) — automatically generated QSP model | 2026-06-24*
