# Alopecia Areata (AA) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Dermatology

[![AA QSP Model](../../../alopecia-areata/aa_qsp_model.png)](../../../alopecia-areata/aa_qsp_model.svg)

## Overview

Alopecia areata (AA) is a T-cell-mediated autoimmune disease targeting the hair follicle, a common form of hair loss affecting approximately 2% of the world's population at least once in their lifetime. Its core pathogenic mechanism is the collapse of the "immune privilege" that the hair follicle normally maintains: once the mechanism suppressing MHC-I expression is lost, CD8⁺ NKG2D⁺ T cells and NK cells attack follicular cells. IFN-γ and IL-15 amplify this immune attack via JAK1/2-STAT1/STAT5 signalling, and further T-cell recruitment via CXCL10 sustains the process. Severity ranges from patchy scalp alopecia areata through alopecia totalis (the entire scalp) to alopecia universalis (the entire body). JAK inhibitors (baricitinib, ritlecitinib, deucravacitinib) have recently been approved and introduced into clinical practice.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Collapse of immune privilege | Loss of MHC-I suppression → exposure of follicular antigens | Initiation of CD8⁺ T-cell attack |
| IL-15/NKG2D axis | IL-15 → increased NKG2D expression on NK cells and CD8⁺ T cells | Follicular cytotoxicity |
| JAK-STAT signalling | IFN-γ/IL-15 → JAK1/2 → STAT1/3/5 | Inflammatory gene expression, sustained immune attack |
| IFN-γ/CXCL10 feedback | CXCL10 secretion → further recruitment of CXCR3⁺ T cells | Amplification cycle of inflammation |
| Treg dysfunction | Reduced/impaired Treg → failure of immune suppression | Loss of self-tolerance |
| Hair follicle cycle | Immune attack → premature termination of anagen | Hair shedding, sparse hair |
| JAK3/TYK2-dependent pathway | IL-4/IL-13 (with concurrent atopy) → Th2 component | Important in atopic AA |

## Drug Targets

- **Baricitinib (JAK1/2 inhibitor)**: significant improvement in SALT score for moderate-to-severe AA in the BRAVE-AA1/2 trials; FDA approved in 2022
- **Ritlecitinib (JAK3/TEC inhibitor)**: ALLEGRO trial; indicated for age 12 and above
- **Ruxolitinib (JAK1/2 inhibitor)**: available as both topical and oral formulations; downregulates IFN-γ/CXCL10
- **Dupilumab (anti-IL-4Rα)**: some efficacy in AA with concurrent atopy; blocks the Th2 pathway
- **Corticosteroids**: topical, intralesional, and systemic; transient effect, limited for long-term treatment

## Model Files

| File | Description |
|------|------|
| [aa_qsp_model.dot](../../../alopecia-areata/aa_qsp_model.dot) | Graphviz mechanistic map source (approximately 164 nodes / 11 clusters) |
| [aa_qsp_model.svg](../../../alopecia-areata/aa_qsp_model.svg) | SVG vector image (scalable) |
| [aa_qsp_model.png](../../../alopecia-areata/aa_qsp_model.png) | PNG image (150 dpi) |
| [aa_mrgsolve_model.R](aa_mrgsolve_model.R) | mrgsolve ODE model (approximately 22 compartments / multiple treatment scenarios) |
| [aa_shiny_app.R](../../../alopecia-areata/aa_shiny_app.R) | Shiny dashboard |
| [aa_references.md](aa_references.md) | References (approximately 55 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: oral PK for baricitinib/ruxolitinib (2 compartments) + dupilumab SC PK + JAK3-binding compartment + PD compartments for NKG2DL, NK cells, naive CD8, effector CD8, Treg, IFN-γ, IL-15, CXCL10, pSTAT1, pSTAT5, immune privilege index, anagen hair, hair density, and systemic inflammation
- **Key treatment scenarios**: ① no treatment (natural course), ② baricitinib 2 mg/day, ③ baricitinib 4 mg/day, ④ ruxolitinib, ⑤ dupilumab (concurrent atopy), ⑥ steroid + JAK inhibitor combination
- **Calibration/evidence**: parameters referenced from BRAVE-AA1/2 (baricitinib) and the foundational immune parameters of Xing et al. Nat Med 2014

## Shiny Dashboard

The dashboard comprises a patient profile tab (SALT score, extent of hair loss, presence of concurrent atopy), a JAK inhibitor PK and target occupancy tab, an immune cell dynamics tab (CD8, Treg, NK), a hair density and SALT score change tab, a treatment scenario comparison tab, and a cytokine biomarker tab (IFN-γ, CXCL10).

## Usage

```r
library(mrgsolve)
mod <- mread("aa_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("aa_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg aa_qsp_model.dot -o aa_qsp_model.svg
```

## References

For detailed citations, see [aa_references.md](aa_references.md) (approximately 55 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
