# Interstitial Cystitis / Bladder Pain Syndrome (IC/BPS)

> **Directory:** `interstitial-cystitis/` | **Abbreviation:** IC/BPS | **Date:** 2026-06-26  
> **Prevalence:** occurs in 2–7% of women and 0.5% of men; an estimated 3–8 million affected in the US  
> **Key keywords:** GAG layer · mast cells · neurogenic inflammation · central sensitization · TRPV1 · bladder pain

---

## Disease overview

**Interstitial cystitis / bladder pain syndrome (IC/BPS)** is a chronic urological
disease characterized by bladder pain, urgency, and frequency without a urinary tract
infection or other obvious cause. It has a serious impact on quality of life and is a
heterogeneous disease in which multiple mechanisms act in combination.

### Main subtypes

| Feature | **Hunner type (inflammatory)** | **Non-Hunner type (majority)** |
|------|----------------------|----------------------|
| Proportion | 5–15% | 85–95% |
| Cystoscopic findings | Hunner lesion (mucosal lesion) | Petechial hemorrhage, glomerulations |
| Main pathology | Lymphocyte/plasma cell infiltration, IgG4+ cells | GAG layer deficiency, mast cell infiltration |
| Inflammatory markers | IL-6, IFN-γ very high | Moderate |
| Recommended treatment | CsA, fulguration | PPS, BoNTA, intravesical instillation therapy |

---

## Mechanistic pathway summary

```
Trigger (stress · infection · autoimmunity · genetics · allergens)
  │
  ▼
GAG layer deficiency (heparan sulfate · chondroitin sulfate · hyaluronic acid deficiency)
  │  → urothelial permeability ↑ → urinary K+ leak → nerve depolarization
  │
  ├─ Mast cell activation → histamine · tryptase · PGE2 · TNF-α
  │    └─ PAR2 activation → urothelial barrier↓ (positive feedback)
  │
  ├─ C-fiber (TRPV1+, P2X3+) sensitization
  │    └─ Retrograde SP · CGRP release → neurogenic inflammation → further mast cell activation
  │
  ├─ NGF ↑ → TrkA → C-fiber proliferation + TRPV1↑
  │
  ├─ Spinal sensitization (wind-up · NMDA · BDNF)
  │    └─ Central sensitization (ACC · insula · prefrontal cortex) → pain amplification
  │
  └─ Bladder wall remodeling
       ├─ TGF-β1 → collagen deposition → fibrosis
       └─ Bladder capacity↓ → frequency · OLS score↑
```

---

## Mechanistic Map

[![IC/BPS QSP Mechanistic Map](ic_bps_qsp_model.png)](ic_bps_qsp_model.svg)

*Click to go to the high-resolution SVG map (140+ nodes, 9 clusters)*

---

## mrgsolve ODE model (`ic_bps_mrgsolve_model.R`)

### Compartment structure (22 ODEs)

| Category | Compartment | Description |
|------|------|------|
| **Drug PK (8)** | PPS_GUT, PPS_CENT | Pentosan polysulfate (F=6%) |
| | HYD_GUT, HYD_CENT | Hydroxyzine (F=80%) |
| | CSA_GUT, CSA_CENT | Cyclosporine A (F=35%) |
| | AMI_GUT, AMI_CENT | Amitriptyline (F=50%) |
| **Disease PD (14)** | GAG | GAG layer integrity (0–1) |
| | PERM | Urothelial permeability (0–1) |
| | MC | Mast cell activity index |
| | HIST | Histamine level |
| | SP | Substance P |
| | NGF | Nerve growth factor |
| | IL6 | Interleukin-6 |
| | TNF | TNF-α |
| | C_FIBER | C-fiber sensitization index |
| | SPINAL | Spinal sensitization index |
| | CENTRAL | Central sensitization index |
| | CAP | Functional bladder capacity (mL) |
| | PAIN | VAS pain score (0–10) |
| | OLS | O'Leary-Sant score (0–20) |

### 7 treatment scenarios

| # | Regimen | Dose | Supporting trial |
|---|------|------|--------------|
| S1 | No treatment (natural course) | — | Natural-course observation |
| S2 | **Oral PPS (Elmiron)** | 100mg TID | Nickel 2005 *Urology*; Hanno 2003 *Urology* |
| S3 | **Hydroxyzine** | 25mg QD | Sant 2003 *J Urol*; Theoharides 1991 |
| S4 | **Intravesical DMSO** | 50% 50mL, q2wk×6 doses | Sant 1987 *Urology*; Fowler 1981 |
| S5 | **Cyclosporine A** | 3mg/kg/day (Hunner type) | Sairanen 2005 *J Urol*; Forrest 2012 |
| S6 | **Botulinum toxin A** | 100U intravesical injection | Gottsch 2011 *J Urol*; Kuo 2010 |
| S7 | **Triple combination therapy** | PPS + hydroxyzine + amitriptyline | Foster 2010 *J Urol* |

### Basis for key parameters

- **GAG synthesis/degradation constants**: Parsons 2007 *Urology* — the GAG deficiency mechanism
- **Mast cell activation rate**: Peeker 2000 *J Urol* — mast cell density in IC/BPS
- **C-fiber sensitization**: Nazif 2007 *Urology* — neural upregulation in IC/BPS
- **PPS bioavailability 6%**: Nickel 2005 — PPS PK data
- **CsA half-life 24h**: Forrest 2012 — calcineurin inhibitor in IC/BPS

---

## Shiny dashboard (`ic_bps_shiny_app.R`)

### Tab structure (8 tabs)

| Tab | Content |
|----|------|
| 1. **Overview** | Disease introduction, QSP structure, subtype comparison table, value boxes |
| 2. **Patient Profile** | Baseline biomarker radar chart, UPOINT domain scores, subtype classification |
| 3. **PK — Drug Levels** | Drug concentration-time profiles, PK parameter table |
| 4. **PD — Biomarkers** | Individual PD biomarker dynamics, heatmap (baseline vs. post-treatment) |
| 5. **Clinical Endpoints** | VAS pain, OLS score, bladder capacity, voiding frequency |
| 6. **Scenario Comparison** | Simultaneous comparison of 7 scenarios, summary table |
| 7. **Subtype Explorer** | Hunner vs. non-Hunner immune profiles, CsA response comparison |
| 8. **Sensitivity Analysis** | 1-way sensitivity analysis, tornado diagram |

### How to run it

```r
# Install R packages
install.packages(c("shiny", "shinydashboard", "mrgsolve", "dplyr",
                   "ggplot2", "tidyr", "DT", "plotly", "patchwork"))

# Run the Shiny app
shiny::runApp("ic_bps_shiny_app.R")
```

---

## References summary (`ic_bps_references.md`)

**62** PubMed citations in total, classified into 17 sections:

1. Clinical guidelines and disease definition (4)
2. Epidemiology and prevalence (3)
3. GAG layer deficiency and the urothelial barrier (5)
4. Mast cell pathophysiology (4)
5. Neurogenic inflammation and sensory pathways (6)
6. Immunopathology and subtypes (5)
7. Central sensitization and pain mechanisms (4)
8. Diagnosis and outcome measures (4)
9. Treatment: PPS (3)
10. Treatment: intravesical DMSO (3)
11. Treatment: hydroxyzine, amitriptyline, multimodal (4)
12. Treatment: cyclosporine A (3)
13. Treatment: botulinum toxin A (3)
14. Treatment: neuromodulation (2)
15. New/investigational treatments (4)
16. Comorbidities and systemic associations (3)
17. QSP modeling methodology (2)

---

## File list

| File | Description |
|------|------|
| [`ic_bps_qsp_model.dot`](ic_bps_qsp_model.dot) | Graphviz mechanistic map (140+ nodes, 9 clusters) |
| [`ic_bps_qsp_model.svg`](ic_bps_qsp_model.svg) | Vector-format map (high resolution) |
| [`ic_bps_qsp_model.png`](ic_bps_qsp_model.png) | Raster-format map (150 dpi) |
| [`ic_bps_mrgsolve_model.R`](ic_bps_mrgsolve_model.R) | mrgsolve ODE QSP model (22 compartments, 7 scenarios, VP n=200) |
| [`ic_bps_shiny_app.R`](ic_bps_shiny_app.R) | Shiny interactive dashboard (8 tabs) |
| [`ic_bps_references.md`](ic_bps_references.md) | 62 references (17 sections) |
| [`README.md`](README.md) | This document |

---

## Clinical application notes

### Evidence-based drug selection guide

```
Confirmed IC/BPS diagnosis
  │
  ├─ Cystoscopy → Hunner lesion confirmed?
  │    ├─ YES (Hunner type)
  │    │    → Fulguration/laser ablation (first line)
  │    │    → CsA 3mg/kg/day (relapse prevention)
  │    │    → Steroids + immunosuppression
  │    │
  │    └─ NO (non-Hunner type)
  │         → Step 1: behavioral therapy + dietary management + physical therapy
  │         → Step 2: PPS 100mg TID + hydroxyzine 25–50mg QD
  │         → Step 3: intravesical instillation (DMSO/heparin/lidocaine)
  │         → Step 4: amitriptyline 25–75mg QD (central sensitization)
  │         → Step 5: BoNTA 100U intravesical injection (q6 months)
  │         → Step 6: neuromodulation (sacral nerve stimulation)
  └─
```

### Predictors of treatment response (QSP model)

| Factor | Good-response group | Poor-response group |
|------|------------|------------|
| GAG layer index | > 50% | < 30% |
| Mast cell activity | Moderate | Very high |
| Central sensitization index | Low | High (with widespread pain) |
| Subtype | Non-Hunner (PPS) | Hunner (CsA required) |
| Bladder capacity | > 200 mL | < 150 mL (BoNTA preferred) |

---

*Model version 1.0 · Claude Code Routine · 2026-06-26*
