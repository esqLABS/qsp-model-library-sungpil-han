# Vitiligo (VIT) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Dermatology

[![VIT QSP Model](../../../vitiligo/vit_qsp_model.png)](../../../vitiligo/vit_qsp_model.svg)

## Overview

Vitiligo is an autoimmune skin disease in which CD8+ T-cell-mediated melanocyte destruction produces depigmented patches on the skin, with a worldwide prevalence of approximately 0.5~2%. Autoreactive CD8+ T cells recognise melanocyte-specific antigens (Melan-A, PMEL17, and others) and cause cytotoxic destruction, and a self-amplifying feedback loop of IFN-γ → JAK1/2-STAT1 → CXCL10 secretion → CD8+ T-cell skin homing is central to lesion expansion. The JAK1/2 inhibitor ruxolitinib (topical) cream was the first to receive FDA approval and demonstrated efficacy for repigmentation, and afamelanotide (an MC1R agonist) and NB-UVB phototherapy are used in combination.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Oxidative stress/melanocyte vulnerability | H₂O₂ accumulation, induced NKG2D ligand expression → NK/NKT cell activation | Initial melanocyte damage and antigen exposure |
| IFN-γ — JAK-STAT1 axis | CD8+ T-cell/NK-cell IFN-γ → JAK1/2 → pSTAT1 → CXCL9/CXCL10 secretion | Amplifies CXCR3+ CD8+ T-cell skin homing |
| CD8+ CTL melanocyte killing | CXCL10 gradient → CD8+ T-cell lesion homing → perforin/granzyme B cytotoxicity | Reduced melanocyte density, expansion of depigmented patches |
| Treg dysfunction | Reduced Foxp3+ Treg skin infiltration → failure to suppress the IFN-γ/CXCL10 axis | Sustained/expanding autoimmune response |
| Reduced MITF signalling | IFN-γ → MITF suppression → reduced melanin synthesis, reduced melanocyte survival | Cessation of melanin production, skin depigmentation |
| Hair follicle melanocyte reservoir | NB-UVB → reactivation/migration of hair-follicle-stored melanocytes | Repigmentation beginning around the hair follicle |

## Drug Targets

- **JAK1/2 inhibitor (topical) — ruxolitinib 1.5% cream**: suppresses pSTAT1 → reduces CXCL10 → blocks CD8+ T-cell homing, repigmentation (TRuE-V1/2 trials)
- **JAK1/2 inhibitor (oral) — ruxolitinib 10 mg BID**: suppresses the systemic IFN-γ/CXCL10 axis, considered for widespread vitiligo
- **MC1R agonist — afamelanotide**: an α-MSH analogue that promotes proliferation and migration of hair follicle melanocytes, accelerating repigmentation when combined with NB-UVB
- **NB-UVB phototherapy**: immune modulation (increases Treg) plus mobilisation of hair-follicle-stored melanocytes
- **Topical calcineurin inhibitor — tacrolimus**: blocks T-cell activation, second-line treatment for facial and flexural vitiligo

## Model Files

| File | Description |
|------|------|
| [vit_qsp_model.dot](../../../vitiligo/vit_qsp_model.dot) | Graphviz mechanistic map source (approximately 262 nodes / 11 clusters) |
| [vit_qsp_model.svg](../../../vitiligo/vit_qsp_model.svg) | SVG vector image (scalable) |
| [vit_qsp_model.png](../../../vitiligo/vit_qsp_model.png) | PNG image (150 dpi) |
| [vit_mrgsolve_model.R](../../../vitiligo/vit_mrgsolve_model.R) | mrgsolve ODE model (approximately 21 compartments / 5 treatment scenarios) |
| [vit_shiny_app.R](../../../vitiligo/vit_shiny_app.R) | Shiny dashboard |
| [vit_references.md](vit_references.md) | References (approximately 32 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK (ruxolitinib oral absorption/plasma/skin, 3 compartments; afamelanotide SC/plasma, 2 compartments) + immune compartments (NKG2D ligand, NKG2D activation, CD8+ CTL, Treg, IFN-γ, CXCL10, pSTAT1) + melanocyte, MITF, melanin synthesis, and hair follicle reservoir compartments + composite metrics (inflammation index, VASI, cumulative repigmentation), for a total of approximately 21 compartments
- **Key treatment scenarios**: (1) placebo (natural progression), (2) ruxolitinib cream BID, (3) ruxolitinib cream QD, (4) oral ruxolitinib 10 mg BID, (5) afamelanotide + NB-UVB combination
- **Calibration/evidence**: parameters referenced from TRuE-V1/2 (24-week F-VASI response to topical ruxolitinib), RECAL (afamelanotide + NB-UVB), and VASI score epidemiological data

## Shiny Dashboard

The dashboard comprises six or more tabs: patient profile (setting body weight, baseline VASI, lesion distribution, and Fitzpatrick skin type), PK visualisation (ruxolitinib skin/plasma concentration-time curves), immune PD metrics (IFN-γ/CXCL10/pSTAT1/CD8+ T-cell time series), clinical endpoints (VASI change, repigmentation area ratio), treatment scenario comparison (long-term melanocyte recovery across 5 regimens), and a biomarker panel (melanin content, hair follicle reservoir, Treg level).

## Usage

```r
library(mrgsolve)
mod <- mread("vit_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("vit_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg vit_qsp_model.dot -o vit_qsp_model.svg
```

## References

For detailed citations, see [vit_references.md](vit_references.md) (approximately 32 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
