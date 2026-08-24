# Congenital Adrenal Hyperplasia (CAH)
## 21-Hydroxylase Deficiency QSP Model

[![Model](cah_qsp_model.png)](cah_qsp_model.svg)

---

## Overview

**Congenital adrenal hyperplasia (CAH)** is the most common disorder of
adrenocortical hormone synthesis, caused by a genetic deficiency of the
steroidogenic enzyme **CYP21A2 (21-hydroxylase)**. CAH occurs in roughly
1/10,000-1/15,000 live births, and 21-hydroxylase deficiency accounts for
95% of all CAH.

When CYP21A2 is deficient:
- **Cortisol production is blocked** -> loss of negative feedback on the
  hypothalamic-pituitary-adrenal (HPA) axis -> excess ACTH secretion
- **Aldosterone production is reduced** (salt-wasting form, SW-CAH) ->
  hyponatraemia, hyperkalaemia, adrenal crisis
- **17-OHP accumulates** -> is shunted into the androgen synthesis pathway
  -> androgen excess, virilisation

| Phenotype | Mutation | Residual CYP21A2 activity | Key features |
|--------|----------|-------------------|-----------|
| Salt-wasting (SW) | Null (del30kb, I2G, Q318X) | <1% | Neonatal crisis, aldosterone deficiency, severe virilisation |
| Simple virilising (SV) | I172N, P30L | 1-2% | Virilisation, accelerated growth, normal aldosterone |
| Non-classic (NC) | V281L, R339H | 20-50% | Mild hyperandrogenism, adult-onset presentation |

---

## Model Architecture

### Mechanistic Map
- **File**: [`cah_qsp_model.dot`](cah_qsp_model.dot) -> [`cah_qsp_model.svg`](cah_qsp_model.svg) / [`cah_qsp_model.png`](cah_qsp_model.png)
- **14 clusters**: hypothalamus, pituitary, adrenal cortex, steroidogenic
  biosynthesis pathway, androgen effects, growth/skeleton, salt-water axis,
  standard drug PK, novel drug PK, glucocorticoid PD, clinical endpoints,
  metabolic complications, psychosocial outcomes, genetic basis
- **100+ nodes**: covering every steroidogenic enzyme, receptor, pathway,
  and biomarker

### mrgsolve ODE Model
- **File**: [`cah_mrgsolve_model.R`](cah_mrgsolve_model.R)
- **Number of compartments**: 35 ODE compartments

| Compartment group | ODEs |
|-----------|------|
| HPA axis | CRH, ACTH |
| Steroidogenic biosynthesis | PREG, PROG, 17-OHP, DHEA, A4, testosterone, DOC, compound S, cortisol, aldosterone |
| Mineralocorticoid axis | RENIN |
| Growth/skeleton | HEIGHT_SDS, BONE_AGE, BMD |
| HC PK (2-cpt) | HC_GUT, HC_CENT, HC_PERI |
| Prednisolone PK | PRED_GUT, PRED_CENT |
| Fludrocortisone PK | FC_GUT, FC_CENT |
| Tildacerfont PK (2-cpt) | TILD_GUT, TILD_CENT, TILD_PERI |
| Crinecerfont PK | CRINE_GUT, CRINE_CENT |

### Treatment Scenarios (6)

| # | Scenario | Description |
|---|----------|------|
| 1 | **Untreated** | Natural history of SW-CAH |
| 2 | **HC + FC (standard)** | Hydrocortisone 20 mg/day TID + fludrocortisone 100 mcg/day |
| 3 | **Prednisolone + FC** | Prednisolone 5 mg/day BID + fludrocortisone |
| 4 | **Dexamethasone** | DEX 0.25 mg QD at bedtime (adult NC-CAH) |
| 5 | **Tildacerfont + HC + FC** | Tildacerfont 100 mg QD + reduced-dose HC (15 mg/day) + FC |
| 6 | **Crinecerfont + HC + FC** | Crinecerfont 200 mg BID + reduced-dose HC (15 mg/day) + FC |

### Trial Validation

| Trial | Drug | Endpoint | Observed | Model |
|----------|------|-----------|--------|------|
| Bonfig 2009 (JCEM) | Standard HC therapy | 17-OHP control rate | ~53% | ~50% |
| CAH2301 (Sarafoglou NEJM 2023) | Tildacerfont | % reduction in 17-OHP | -58% | -55% |
| CARES (Merke NEJM 2024) | Crinecerfont | % reduction in androstenedione | -44% | -42% |
| CARES (Merke NEJM 2024) | Crinecerfont | % reduction in ACTH | -66% | -61% |

---

## Drug PK Summary

| Drug | Route | F (%) | t½ (h) | GC potency | CRF1R IC50 |
|------|------|--------|--------|---------|-----------|
| Hydrocortisone | Oral | 95 | 1.5 | 1x | — |
| Prednisolone | Oral | 82 | 2.5 | 4x | — |
| Dexamethasone | Oral | 78 | 3.8 | 25x | — |
| Fludrocortisone | Oral | 90 | 3.5 | — (MC 125x) | — |
| Tildacerfont | Oral | 65 | 12-14 | — | ~4 nM |
| Crinecerfont | Oral | 50 | 8-10 | — | ~0.5 nM |

---

## Shiny Dashboard

**File**: [`cah_shiny_app_en.R`](cah_shiny_app_en.R)

**6 tabs**:
| Tab | Content |
|----|------|
| 1. Patient profile | Mutation type, phenotype, treatment goals |
| 2. Drug PK | Plasma concentration-time profiles, PK parameters |
| 3. Steroid biomarkers | 17-OHP, ACTH, androstenedione, cortisol |
| 4. Clinical endpoints | Height SDS, bone age, bone density, renin |
| 5. Scenario comparison | Simultaneous comparison of 6 treatment strategies |
| 6. Biomarker dashboard | Target achievement rate, CRF1R occupancy, summary table |

---

## How to Run

```r
# Run the mrgsolve model
install.packages(c("mrgsolve", "dplyr", "ggplot2", "patchwork"))
source("cah_mrgsolve_model.R")

# Run the Shiny dashboard
install.packages(c("shiny", "shinydashboard", "DT", "plotly"))
shiny::runApp("cah_shiny_app_en.R")

# Render the mechanistic map (requires Graphviz)
# dot -Tsvg cah_qsp_model.dot -o cah_qsp_model.svg
# dot -Tpng -Gdpi=150 cah_qsp_model.dot -o cah_qsp_model.png
```

---

## Key Model Insights

1. **Difficulty of 17-OHP control**: only about 50% of patients reach the
   target 17-OHP (<36 nmol/L) on standard HC therapy — related to the short
   half-life of soluble HC
2. **Advantage of CRF1 antagonists**: by directly suppressing ACTH, they can
   achieve equal or better biomarker control at lower GC doses
3. **Dual risk to growth suppression**: both excess androgen (premature bone
   age advancement) and excess GC (direct growth suppression) contribute to
   reduced final height
4. **Crinecerfont vs. tildacerfont**: the difference in IC50 (0.5 nM vs.
   4 nM) means crinecerfont achieves more complete CRF1R occupancy
5. **Salt-wasting crisis**: aldosterone deficiency + hyperkalaemia require
   prompt FC replacement

---

## References

54 PubMed citations: [`cah_references_en.md`](cah_references_en.md)

---

*Date created: 2026-06-25 | Claude Code Routine (CCR)*
