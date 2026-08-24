# Schizophrenia QSP Model

**Category**: Neuropsychiatric Disorder
**Pathogenesis**: Dopamine hyperactivity (mesolimbic) / hypoactivity (mesocortical) + NMDA hypofunction + GABAergic PV interneuron deficit + serotonergic abnormalities + neuroinflammation
**Date**: 2026-06-20

---

## 1. Mechanistic Map

[![Schizophrenia QSP Model](../../../schizophrenia/sch_qsp_model.png)](../../../schizophrenia/sch_qsp_model.svg)

> Click to view the full-resolution SVG.

### Key Clusters (10)

| # | Cluster | Key content |
|---|----------|-----------|
| ① | Neurodevelopmental risk factors | DISC1, NRG1, DTNBP1, COMT Val158Met, C4A; environmental stress, cannabis |
| ② | Dopamine pathways | Mesolimbic (↑), mesocortical (↓), nigrostriatal (EPS), tuberoinfundibular (PRL); D1/D2/D3 receptors |
| ③ | Glutamate/NMDA | NMDA hypofunction → disinhibition of PV interneurons → cortical disinhibition → ↑ subcortical DA |
| ④ | Serotonin | DRN, 5-HT2A/2C/1A; SGA blockade of 5-HT2A is key to restoring mesocortical DA |
| ⑤ | GABAergic interneurons | PV+ cell deficit (GAD67 ↓ 25-50%), loss of gamma oscillations → cognitive deficits |
| ⑥ | Neuroinflammation/oxidative stress | ↑ IL-6, IL-1β, TNF-α; C4A complement → excess synaptic pruning; ↓ BDNF |
| ⑦ | Antipsychotic PK | HAL, RIS/PALI, CLZ, ARI — 2-compartment models, including metabolites |
| ⑧ | Pharmacodynamics/receptor occupancy | D2 occupancy 65-80% = therapeutic window; 5-HT2A occupancy >80% = SGA-driven negative-symptom improvement |
| ⑨ | Clinical endpoints | PANSS positive/negative/general, RBANS cognition, relapse rate, functional recovery |
| ⑩ | Adverse effects | EPS, tardive dyskinesia, metabolic syndrome, hyperprolactinaemia, QTc prolongation, agranulocytosis |

---

## 2. mrgsolve ODE Model (`sch_mrgsolve_model.R`)

### Compartments (22)

| Compartment type | Compartment names |
|-----------|-----------|
| **HAL PK** | GUT_HAL, CENT_HAL, PERI_HAL |
| **RIS/PALI PK** | GUT_RIS, CENT_RIS, PERI_RIS, CENT_PALI |
| **CLZ PK** | CENT_CLZ |
| **ARI/dARI PK** | GUT_ARI, CENT_ARI, CENT_dARI |
| **Dopamine PD** | DA_MESOLIM, DA_MESOCORT, DA_NIGROSTR |
| **Other PD** | PRL_CMPT, PV_ACT |
| **Clinical** | PANSS_POS, PANSS_NEG, PANSS_GEN |
| **Biomarkers** | BDNF_CMPT, IL6_CMPT, EPS_RISK |

### Drug PK Parameters (key)

| Drug | F (%) | CL (L/h) | Vc (L) | t½ | Brain Kp |
|------|--------|-----------|--------|-----|----------|
| Haloperidol (FGA) | 65 | 15 | 20 | 18-24h | 12 |
| Risperidone (SGA) | 74 | 25 | 30 | 3h (→21h PALI) | 7 |
| Clozapine (TRS) | 55 | 30 | 50 | 12h | 6 |
| Aripiprazole (partial D2) | 87 | 3.6 | 245 | 75h | ~15 |

### Treatment Scenarios (7)

| # | Scenario | Represents |
|---|----------|------|
| 1 | No treatment (natural course) | Comparator baseline |
| 2 | Haloperidol 10 mg/d (standard FGA) | 1st generation |
| 3 | Haloperidol 5 mg/d (low dose) | Low-dose FGA |
| 4 | Risperidone 4 mg/d (SGA) | 2nd-generation standard |
| 5 | Clozapine 300 mg/d (TRS) | Treatment-resistant |
| 6 | Aripiprazole 15 mg/d (partial D2) | Partial agonist |
| 7 | Risperidone 2 mg/d (low-dose SGA) | Low-dose SGA |

### Key Clinical Calibration Sources
- **CATIE 2005** (Lieberman JA, NEJM): comparative efficacy of antipsychotics
- **EUFEST 2008** (Kahn RS, Lancet): first-episode schizophrenia drug comparison
- **Kapur 2000** (AJP): PET study of the D2 occupancy therapeutic window
- **Nordstrom 1995** (AJP): PET study of clozapine's low D2 occupancy

---

## 3. Shiny Dashboard (`sch_shiny_app.R`)

**6-tab layout**:

| Tab | Content |
|----|------|
| ① Patient profile | Demographics, baseline PANSS settings, treatment selection, pathogenesis summary |
| ② PK profile | Plasma concentration-time curves, brain concentration, approach to steady state |
| ③ D2/5-HT2A occupancy | Receptor occupancy vs. time, therapeutic-window visualisation, drug fingerprint |
| ④ PANSS clinical endpoints | PANSS total/subscale scores, dopamine-pathway dynamics, response rate |
| ⑤ Scenario comparison | Simultaneous comparison of the 7 treatments (PANSS, D2 occupancy, EPS) |
| ⑥ Biomarkers | Prolactin, EPS risk index, BDNF, IL-6, PV interneuron activity |

**How to run**:
```r
library(shiny); runApp("sch_shiny_app.R")
```

---

## 4. References (`sch_references.md`)

45 PubMed references (11 sections):
- Epidemiology and clinical overview
- The dopamine hypothesis
- The glutamate/NMDA receptor hypothesis
- GABAergic interneurons and circuit dysfunction
- The serotonergic system
- Neuroinflammation and oxidative stress
- Genetics and neurodevelopment
- Antipsychotic pharmacokinetics
- D2/5-HT2A occupancy (PET studies)
- Clinical trials and treatment efficacy
- QSP/computational modelling

---

## 5. File List

```
schizophrenia/
├── sch_qsp_model.dot        # Graphviz mechanistic map (10 clusters, 160+ nodes)
├── sch_qsp_model.svg        # SVG vector image
├── sch_qsp_model.png        # PNG image (150 dpi)
├── sch_mrgsolve_model.R     # mrgsolve ODE model (22 compartments, 7 scenarios)
├── sch_shiny_app.R          # Shiny dashboard (6 tabs)
├── sch_references.md        # 45 references
└── README.md                # this file
```

---

## 6. Key Summary of Schizophrenia Pathophysiology

### The Dual Dopamine Hypothesis

```
Mesolimbic pathway:       VTA → NAc    (↑ in SCZ) → positive symptoms (hallucinations, delusions)
Mesocortical pathway:     VTA → DLPFC  (↓ in SCZ) → negative symptoms + cognitive deficits
Nigrostriatal pathway:    SNc → striatum (D2 >80% block) → EPS
Tuberoinfundibular-pituitary: hypothalamus→pituitary (D2 block) → hyperprolactinaemia
```

### Differences Between SGA and FGA

| | FGA (e.g. haloperidol) | SGA (e.g. risperidone) |
|-|------------------------|------------------------|
| D2 occupancy | High (~78%) | High (~75%) |
| 5-HT2A occupancy | Low (~30%) | High (~96%) |
| Positive symptoms | Good | Good |
| Negative symptoms | Limited | Improved (5-HT2A → ↑ mesocortical DA) |
| EPS | High | Low (5-HT2A → restores nigrostriatal DA) |
| Metabolic effects | Low | Moderate-high (5-HT2C, H1, M1 blockade) |

### The Distinctiveness of Aripiprazole
- **Partial D2 agonist**: acts as an antagonist under DA hyperactivity; acts as an agonist under DA hypoactivity
- **Dopamine stabilisation**: balances mesolimbic (↓) with mesocortical (↑)
- **5-HT1A partial agonist**: restores PV interneuron function + ↑ mesocortical DA
