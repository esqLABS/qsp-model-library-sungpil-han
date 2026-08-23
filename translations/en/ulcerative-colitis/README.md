# Ulcerative Colitis (UC) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Gastroenterology/Hepatobiliary

[![UC QSP Model](../../../ulcerative-colitis/uc_qsp_model.png)](../../../ulcerative-colitis/uc_qsp_model.svg)

## Overview

Ulcerative colitis (UC) is a chronic relapsing inflammatory bowel disease confined to the colonic mucosa, beginning in the rectum and extending continuously into the proximal colon. Prevalence in Western populations is approximately 100~300 per 100,000, and it is also steadily increasing in Asia and Korea. Its core pathogenic mechanism can be summarised as impaired barrier function → exposure to gut microbial antigens → Th2/Th17-driven mucosal inflammation → excess secretion of TNF-α, IL-17, and IL-13. Therapies with diverse mechanisms, including 5-ASA, biologics (anti-TNF, anti-integrin, anti-IL-12/23), JAK inhibitors, and S1P modulators, are used in a stepwise fashion.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Impaired barrier function | Reduced claudin/occludin expression, reduced mucin layer | Mucosal penetration by gut bacterial antigens |
| Th2-driven inflammation | TSLP/IL-33 → ILC2 → excess IL-13/IL-5 secretion | Mucosal epithelial damage, eosinophil infiltration |
| TNF-α pathway | Macrophages/Th1 → TNF-α → NF-κB → mucosal ulceration | Bloody/mucous stool, abdominal pain, raised Mayo score |
| IL-17 pathway | Th17 → IL-17A → neutrophil recruitment, CXCL8 induction | Neutrophilic infiltration of the intestinal mucosa, worsening colitis |
| Regulatory T-cell dysfunction | Weakened Treg→IL-10 pathway, breakdown of immune tolerance | Chronic relapsing course |
| Gut homing pathway | α4β7-MAdCAM-1 axis, CCR9-CCL25 → lymphocyte homing to the colon | A targetable pathway (vedolizumab) |

## Drug Targets

- **Anti-TNF-α — infliximab**: neutralises TNF-α, induces mucosal healing (ACT-1/2 trials)
- **Anti-α4β7 integrin — vedolizumab**: blocks gut-selective lymphocyte homing (GEMINI 1/2)
- **JAK inhibitor — tofacitinib**: inhibits JAK1/3, blocks cytokine signalling (OCTAVE Induction/Sustain)
- **Anti-IL-12/23 — ustekinumab**: blocks the IL-12/IL-23 p40 subunit, suppresses Th1/Th17 (UNIFI)
- **S1P receptor modulator — ozanimod**: sequesters lymphocytes in secondary lymphoid organs, blocks gut homing (TRUE NORTH)

## Model Files

| File | Description |
|------|------|
| [uc_qsp_model.dot](../../../ulcerative-colitis/uc_qsp_model.dot) | Graphviz mechanistic map source (approximately 422 nodes / 12 clusters) |
| [uc_qsp_model.svg](../../../ulcerative-colitis/uc_qsp_model.svg) | SVG vector image (scalable) |
| [uc_qsp_model.png](../../../ulcerative-colitis/uc_qsp_model.png) | PNG image (150 dpi) |
| [uc_mrgsolve_model.R](../../../ulcerative-colitis/uc_mrgsolve_model.R) | mrgsolve ODE model (approximately 35 compartments / 6 treatment scenarios) |
| [uc_shiny_app.R](../../../ulcerative-colitis/uc_shiny_app.R) | Shiny dashboard |
| [uc_references.md](../../../ulcerative-colitis/uc_references.md) | References (approximately 47 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK (infliximab, 2 compartments + TMDD; vedolizumab, 1 compartment + TMDD; tofacitinib oral, 2 compartments; ozanimod plus active metabolite, 3 compartments; ustekinumab SC, 2 compartments) + cytokines (TNFα, IL-17, IL-13, IL-10) + immune cells (Th2, Th17, Treg, neutrophils) + disease markers (Mayo score, mucosal healing index, CRP, faecal calprotectin), for a total of approximately 35 compartments
- **Key treatment scenarios**: (1) placebo, (2) infliximab 5 mg/kg IV induction/maintenance, (3) vedolizumab SC, (4) tofacitinib 10 mg BID induction → 5 mg BID maintenance, (5) ustekinumab IV induction → SC maintenance, (6) ozanimod 0.92 mg/day
- **Calibration/evidence**: parameters referenced from ACT-1/2 (infliximab mucosal healing rate), GEMINI (vedolizumab clinical remission rate), OCTAVE (tofacitinib Mayo score response), UNIFI (ustekinumab), and TRUE NORTH (ozanimod) clinical trial data

## Shiny Dashboard

The dashboard comprises six or more tabs: patient profile (setting body weight, baseline Mayo score, CRP, and faecal calprotectin), PK visualisation (blood concentration-time curves for each biologic and small molecule), mucosal PD metrics (Mayo score and mucosal healing index time series), clinical endpoints (comparison of clinical remission and endoscopic response rates), treatment scenario comparison (long-term efficacy across 6 regimens), and a biomarker panel (cytokine and immune cell profiles).

## Usage

```r
library(mrgsolve)
mod <- mread("uc_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("uc_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg uc_qsp_model.dot -o uc_qsp_model.svg
```

## References

For detailed citations, see [uc_references.md](../../../ulcerative-colitis/uc_references.md) (approximately 47 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
