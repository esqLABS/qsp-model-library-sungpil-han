# Allergic Rhinitis (AR) — QSP Model

[![AR Mechanistic Map](ar_qsp_model.png)](ar_qsp_model.svg)

**Category:** Allergy / Immunology
**Directory:** `allergic-rhinitis/`
**ARIA classification:** Intermittent/persistent x mild/moderate-severe

---

## Pathophysiology

Allergic rhinitis (AR) is an **IgE-mediated type I hypersensitivity
reaction**, affecting roughly 10-40% of the world's population, and is a
leading chronic allergic disease with tens of billions of dollars in annual
healthcare costs and lost productivity. It arises from repeated exposure to
inhaled allergens such as pollen, house dust mites, animal dander, and
mould spores, and is linked to asthma, allergic conjunctivitis, chronic
sinusitis, and the atopic march.

| Core pathological pathway | Detailed mechanism |
|-------------|------------|
| **Sensitisation phase** | Allergen -> capture by dendritic cells (DC) -> MHC II presentation -> Th2 differentiation (GATA3+) -> B-cell IgE class switching |
| **IgE-mast cell axis** | Free IgE -> FcεRI binding -> mast cell sensitisation; on re-exposure, IgE-antigen crosslinking -> degranulation |
| **Early-phase response (0-60 min)** | Release of histamine, PGD2, LTC4/D4, PAF -> H1R/CysLT1R activation -> sneezing, rhinorrhoea, itching |
| **Late-phase response (4-24 h)** | Th2 (IL-4/5/13)/ILC2 -> eosinophil influx (eotaxin/CCL11) -> nasal congestion, sustained inflammation |
| **Epithelial alarmin pathway** | Allergen proteases -> barrier disruption -> TSLP/IL-33/IL-25 secretion -> ILC2 activation |
| **Neurogenic inflammation** | H1R -> trigeminal nerve activation -> axon reflex -> SP/CGRP release -> NANC vasodilation/secretion |
| **Chronic remodelling** | IL-13 -> goblet cell hyperplasia; MBP/ECP (eosinophil granule proteins) -> epithelial damage -> reduced olfaction |

---

## Model Architecture

### Mechanistic Map (`ar_qsp_model.dot/.svg/.png`)

12 subgraph clusters, **100+ nodes**:

| Cluster | Key content |
|---------|---------|
| 1. Allergen exposure | Pollen/house dust mite/animal dander/mould, nasal mucosal contact, protease activity |
| 2. Epithelial barrier/innate alarmins | Epithelial TSLP/IL-33/IL-25 secretion, ILC2 activation, NLRP3 inflammasome |
| 3. APC sensitisation/Th2 differentiation | Plasmacytoid/myeloid DC, MHC II presentation, OX40L, Th2/Th1/Treg differentiation |
| 4. B cells/IgE class switching | IgE CSR, plasma cells, free IgE, ImmunoCAP-specific IgE |
| 5. Mast cell/basophil activation | FcεRI, IgE-antigen crosslinking, degranulation, Syk/PLCγ signalling |
| 6. Early-phase inflammatory mediators | Histamine, tryptase, heparin, PGD2, TXA2, LTC4/D4/E4, PAF, SP, CGRP, VIP |
| 7. Receptor activation | H1/H2/H4R, DP1/CRTH2, CysLT1/2R, PAFR, TRKR, α1/α2R |
| 8. Late-phase response/eosinophils | IL-4/5/13/9/31, eotaxin-1/3/MCP-4/RANTES, MBP/ECP/EPX/EDN |
| 9. Nasal physiology/symptoms | Vasodilation, glands, nerves, congestion, TNSS (0-12)/RQLQ |
| 10. PK: oral H1 antihistamines | Cetirizine, loratadine, fexofenadine, bilastine (absorption -> plasma -> H1RO) |
| 11. PK: intranasal corticosteroids | Fluticasone, mometasone, budesonide, triamcinolone (GR translocation -> transcriptional repression) |
| 12. Biologics/LTRA/immunotherapy | Montelukast, omalizumab, dupilumab, mepolizumab, SCIT/SLIT, oxymetazoline |

---

### mrgsolve ODE Model (`ar_mrgsolve_model.R`)

**25 state variables (ODE compartments):**

| Compartment | Variable | Description |
|------|------|------|
| Allergen | `AG` | Nasal allergen burden (AU) |
| IgE | `IGE_FREE`, `IGE_MAST` | Free IgE, mast-cell-bound IgE |
| Mast cells | `MAST_ACT`, `MAST_CHG` | Activation, granule loading status |
| Mediators | `HISTAMINE`, `CYS_LT` | Histamine, cysteinyl leukotrienes |
| Th2/cytokines | `TH2`, `IL4`, `IL5`, `IL13` | Th2 activity, pg/mL |
| Eosinophils | `EOS_B`, `EOS_N` | Blood (cells/μL), nasal tissue |
| Cetirizine PK | `CETI_D`, `CETI_C` | Gut, central compartments |
| Fluticasone PK | `FP_LOC` | Local nasal concentration (nM) |
| Montelukast PK | `MLKT_D`, `MLKT_C` | Gut, central compartments |
| Omalizumab PK | `OMA_D`, `OMA_C`, `OMA_P`, `OMA_IGE` | SC -> two-compartment/IgE complex |

**Clinical trial calibration parameters:**
- Cetirizine at Css gives H1 receptor occupancy >=80% (Yanai 1995, *Allergy*)
- FP 200 μg/d -> 35-40% reduction in TNSS vs. placebo (Meltzer 2005, *JACI*)
- Montelukast 10 mg/d -> 20-25% reduction in TNSS (Philip 2002, *Clin Exp Allergy*)
- Omalizumab -> >95% reduction in free IgE (Fahy 1997, *AJRCCM*)
- Nasal eosinophilia reduced by ~50% with FP treatment (Holgate 2003, *Allergy*)

---

### Treatment Scenarios (7)

| No. | Scenario | Key mechanism |
|------|---------|-------------|
| 1 | Natural history (allergen alone) | Untreated control |
| 2 | Cetirizine 10 mg/d | H1R blockade -> suppresses histamine effect |
| 3 | Fluticasone (FP) 200 μg/d | GR activation -> suppresses IL-4/5/13 and eosinophils |
| 4 | Montelukast 10 mg/d | CysLT1R blockade -> suppresses congestion and vascular permeability |
| 5 | Cetirizine + fluticasone (combination) | Simultaneous H1R + GR blockade |
| 6 | Omalizumab 300 mg q4w | Free IgE capture -> FcεRI downregulation |
| 7 | Triple combination (cetirizine + FP + MLKT) | Multi-pathway blockade |

---

### Shiny App (`ar_shiny_app.R`)

**6 tabs:**

| Tab | Content |
|----|------|
| 1. Patient profile | ARIA classification, baseline IgE/eosinophils/IL-5/TNSS, pathophysiology summary |
| 2. Drug PK | Cetirizine/montelukast/omalizumab plasma concentration time series, H1RO/GR/CysLT1 occupancy |
| 3. Mediators/biomarkers | Histamine, CysLT, IL-4/5/13, eosinophils (blood + nasal), free IgE, mast cell activity |
| 4. Symptoms/TNSS | TNSS (0-12) time series, individual symptom scores (sneezing/rhinorrhoea/congestion/itching, each 0-3) |
| 5. Scenario comparison | Overlaid TNSS graphs for 7 treatment scenarios + 12-week outcome summary table |
| 6. Mechanistic map | Description of the 12-cluster model structure |

---

## Key Clinical Trial Evidence

| Drug | Trial | Key result |
|------|---------|---------|
| Cetirizine | Yanai 1995 *Allergy* | PET study: Css H1-RO >=80% |
| Fluticasone (FP) | Meltzer 2005 *JACI* | TNSS -38% vs. placebo (8 weeks, n=615) |
| Montelukast | Philip 2002 *Clin Exp Allergy* | Seasonal AR TNSS -22% |
| Omalizumab | Meltzer 2000 *JACI* | Improved nasal symptom score in perennial AR |
| Dupilumab | Bachert 2016 *NEJM* | Nasal polyposis + AR NPS -2.06 |
| SCIT (pool) | Durham 1999 *NEJM* | Effect sustained beyond 3 years, specific IgG4 up |

---

## Files

| File | Description |
|------|------|
| `ar_qsp_model.dot` | Graphviz DOT mechanistic map (12 clusters, 100+ nodes) |
| `ar_qsp_model.svg` | SVG vector image |
| `ar_qsp_model.png` | PNG raster image (150 dpi) |
| `ar_mrgsolve_model.R` | mrgsolve ODE QSP model (25 compartments, 7 scenarios) |
| `ar_shiny_app.R` | Shiny dashboard (6 tabs) |
| `ar_references.md` | 60 references (epidemiology, pathophysiology, pharmacology, clinical trials, PK/PD) |
| `README.md` | This document |

---

## Disclaimer

This model is a qualitative/semi-quantitative QSP model for **education and
research purposes**. It has not been independently validated or certified
and must not be used directly for real clinical decision-making.
