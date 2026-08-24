# Waldenström's Macroglobulinemia (WM) — QSP Model

[![WM](wm_qsp_model.png)](wm_qsp_model.svg)

**Category:** Hematologic Oncology · **Directory:** [`waldenstrom-macroglobulinemia/`](.)

---

## Pathophysiology

Waldenström's macroglobulinemia (WM) is a rare B-cell neoplasm classified as
**lymphoplasmacytic lymphoma (LPL)** in the bone marrow (BM) and lymph nodes. The
diagnostic criteria are ① ≥10% lymphoplasmacytic infiltration in the bone marrow and
② detection of a monoclonal IgM paraprotein, and it accounts for about 1-2% of all
non-Hodgkin lymphomas.

### Key pathogenic mechanisms

| Pathway | Detailed mechanism |
|------|------------|
| **MYD88 L265P mutation** (~95%) | Spontaneous MyD88 dimerization → IRAK4/IRAK1 → TRAF6 → IKKβ → constitutive NF-κB activation; a direct MYD88→HCK→BTK axis |
| **CXCR4 WHIM mutation** (~35%) | Reduced receptor internalization → sustained CXCL12/SDF-1 signaling → PI3K/AKT → increased BM retention; associated with ibrutinib resistance |
| **BCR/BTK signaling** | Antigen→BCR→LYN/SYK→BTK→PLCγ2→PKCβ/DAG/IP3→CBM complex→IKK→NF-κB |
| **PI3K/AKT/mTOR** | BCR and CXCR4 → PI3K-δ/γ → PIP3 → AKT(pT308) → mTORC1 → protein synthesis/proliferation |
| **BM microenvironment** | Stromal cell CXCL12/BAFF/APRIL → tumor cell survival; mast cell IL-6 paracrine signaling; Treg suppression of NK cells |
| **Excess IgM production** | LPC → plasma cell differentiation (BLIMP1/IRF4) → pentameric IgM secretion → hyperviscosity syndrome, cryoglobulinemia, anti-MAG neuropathy |

---

## Model Specifications

| Component | Content |
|----------|------|
| **Mechanistic map** | 12 subgraph clusters, 120+ nodes (Drug PK · BTK/BCR · MYD88/NF-κB · PI3K/AKT/mTOR · CXCR4/BM · LPC differentiation · TME · apoptosis · IgM complications · drug PD · clinical outcomes · genomics) |
| **mrgsolve ODE** | 20 compartments (Ibrutinib/Zanubrutinib/Rituximab/Venetoclax PK + BTK occupancy · NF-κB · LPC · PC · IgM · Hgb · viscosity · BM infiltration · BCL-2 · CD20 · apoptosis · proteasome · NK) |
| **Treatment scenarios** | 7 (Watch & Wait · Ibrutinib · Ibrutinib+Rituximab · Zanubrutinib · R-Bendamustine · BDR · Venetoclax) |
| **Shiny app** | 7 tabs (Patient Profile/IPSSWM · Drug PK · Key PD · Clinical Endpoints · Scenario Comparison · Biomarkers · About) |
| **References** | 62 (diagnosis · epidemiology · molecular biology · clinical trials · PK/PD · complications · new drugs) |

---

## Clinical Trial Calibration Data

| Trial | Regimen | ORR | VGPR/CR | PFS |
|------|------|-----|---------|-----|
| iNNOVATOR (Treon 2015) | Ibrutinib 420 mg/d | 91.5% | 30.4% | 69% at 2 years |
| INNOVATE (Dimopoulos 2018) | Ibrutinib + Rituximab | 92% | 43% | 82% at 30 months |
| ASPEN (Tam 2020) | Zanubrutinib 160 mg BID | 93.7% | 28.4% | 84% at 18 months |
| Rummel 2013 | R-Bendamustine | 96% | 44% | mPFS 69 mo |
| Dimopoulos 2013 | BDR | 83% | 22% | mPFS 43 mo |
| Castillo 2018 | Venetoclax | 84% | 36% | Not reached |

---

## Key Drug Mechanisms

| Drug | Target | Mechanism |
|------|------|------|
| **Ibrutinib** | BTK C481 (covalent) | Blocks the BCR/MYD88→BTK pathway → NF-κB ↓, inhibits proliferation |
| **Zanubrutinib** | BTK C481 (covalent) | More selective than ibrutinib; minimizes EGFR/ITK-related side effects |
| **Rituximab** | CD20 | ADCC · CDC · direct apoptosis; caution for IgM flare |
| **Bortezomib** | 26S Proteasome β5 | Stabilizes IκBα → NF-κB ↓; UPR → ER stress-induced apoptosis |
| **Venetoclax** | BCL-2 BH3 (Ki <1 nM) | BAX/BAK release → MOMP → caspase cascade |
| **Bendamustine** | DNA double strands | Alkylation plus purine-analogue effect |

---

## Files

| File | Description |
|------|------|
| [`wm_qsp_model.dot`](wm_qsp_model.dot) | Graphviz mechanistic map source |
| [`wm_qsp_model.svg`](wm_qsp_model.svg) | Vector map (interactive view) |
| [`wm_qsp_model.png`](wm_qsp_model.png) | Raster map (150 dpi) |
| [`wm_mrgsolve_model.R`](wm_mrgsolve_model.R) | mrgsolve ODE model + simulation of 7 scenarios |
| [`wm_shiny_app.R`](wm_shiny_app.R) | 7-tab interactive Shiny dashboard |
| [`wm_references_en.md`](wm_references_en.md) | 62 references (classified by section) |

---

## Usage

```r
# mrgsolve model
library(mrgsolve)
source("wm_mrgsolve_model.R")   # compile the model + automatically run 7 scenarios

# Shiny dashboard
shiny::runApp("waldenstrom-macroglobulinemia/wm_shiny_app.R")
```

```bash
# Graphviz rendering
dot -Tsvg wm_qsp_model.dot -o wm_qsp_model.svg
dot -Tpng -Gdpi=150 wm_qsp_model.dot -o wm_qsp_model.png
```

---

*Date created: 2026-06-27 · Claude Code Routine (CCR)*
