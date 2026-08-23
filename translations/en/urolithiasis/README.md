# Chronic Recurrent Urolithiasis (URI) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README.md](../README.md) · Category: Renal/Urological

[![URI QSP Model](../../../urolithiasis/uri_qsp_model.png)](../../../urolithiasis/uri_qsp_model.svg)

## Overview

Chronic recurrent urolithiasis arises through the pathway of urinary supersaturation of calcium, oxalate, uric acid, and citrate → crystal nucleation → attachment at the base of a Randall plaque → stone growth, and is a common urological disease with a lifetime recurrence rate reaching 50%. Worldwide prevalence is approximately 10~15% and rising in developed countries. Causes range from genetic conditions such as primary hyperoxaluria type 1 (PH1) to metabolic abnormalities such as hypercalciuria and hyperuricosuria. Increased fluid intake and a low-salt, low-protein diet are foundational, and drugs such as thiazide diuretics, potassium citrate, allopurinol, and tamsulosin are selected according to stone type.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Urinary supersaturation | Excess excretion of calcium, oxalate, phosphate, and uric acid → increased ion activity product | Nucleation of CaOx, CaP, and uric acid crystals |
| Randall plaque | Apatite deposition in the renal papillary basement membrane → an anchor site for calcium oxalate stone attachment | Increased risk of recurrent stone growth |
| Reduced citrate inhibition | Hypocitraturia → deficiency of a crystallisation inhibitor | Promotes stone growth, increases recurrence risk |
| PH1 hyperoxaluria | AGXT mutation → hepatic oxalate overproduction → severe renal oxalate deposition | Early renal failure, renal and systemic oxalosis |
| Renal inflammation/injury | Crystal deposition → oxidative stress → tubular epithelial injury → inflammation | Reduced eGFR, progression of chronic renal failure |
| Urinary pH regulation | pH ≤ 5.5 → uric acid crystallisation; pH ≥ 7 → risk of CaP precipitation | Target pH differs by stone type |

## Drug Targets

- **Thiazide diuretics — HCTZ 25 mg/day**: promotes distal tubular calcium reabsorption → reduces urinary calcium (hypercalciuric type)
- **Potassium citrate — 60 mEq/day**: increases urinary citrate plus urinary alkalinisation → inhibits CaOx and uric acid stones
- **Allopurinol — 300 mg/day**: inhibits XO → reduces uric acid production → prevents uric acid stones and hyperuricosuric CaOx stones
- **Lumasiran** — an siRNA-based AGXT-complementary therapy: inhibits hepatic HAO1 in PH1 → reduces oxalate production (ILLUMINATE-A/B)
- **Tamsulosin — an α1-blocker**: relaxes ureteral smooth muscle → promotes spontaneous stone passage

## Model Files

| File | Description |
|------|------|
| [uri_qsp_model.dot](../../../urolithiasis/uri_qsp_model.dot) | Graphviz mechanistic map source (approximately 431 nodes / 10 clusters) |
| [uri_qsp_model.svg](../../../urolithiasis/uri_qsp_model.svg) | SVG vector image (scalable) |
| [uri_qsp_model.png](../../../urolithiasis/uri_qsp_model.png) | PNG image (150 dpi) |
| [uri_mrgsolve_model.R](../../../urolithiasis/uri_mrgsolve_model.R) | mrgsolve ODE model (approximately 19 compartments / 6 treatment scenarios) |
| [uri_shiny_app.R](../../../urolithiasis/uri_shiny_app.R) | Shiny dashboard |
| [uri_references.md](../../../urolithiasis/uri_references.md) | References (approximately 38 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: drug PK (HCTZ, 3 compartments; allopurinol + oxypurinol, 3 compartments; potassium citrate, 2 compartments; tamsulosin, 2 compartments) + physiological compartments (plasma Ca, oxalate, uric acid, citrate: 4 compartments) + stone size, renal inflammation, and eGFR (3 compartments), for a total of approximately 19 compartments
- **Key treatment scenarios**: (1) no treatment (CaOx stone former), (2) HCTZ 25 mg/day, (3) potassium citrate 60 mEq/day, (4) allopurinol 300 mg/day, (5) PH1 + lumasiran, (6) lifestyle modification plus combination therapy
- **Calibration/evidence**: parameters referenced from ILLUMINATE-A/B (lumasiran in PH1), 24-hour urine analysis data from the AUA stone guidelines, and thiazide meta-analyses (effect on reducing urinary calcium)

## Shiny Dashboard

The dashboard comprises six or more tabs: patient profile (setting body weight, fluid intake, stone type, baseline eGFR, and presence of a genetic cause), PK visualisation (blood concentration-time curves for each drug), urine chemistry PD metrics (urinary calcium/oxalate/citrate and urine pH time series), clinical endpoints (stone size trajectory, eGFR change, recurrence rate), treatment scenario comparison (long-term stone growth suppression across 6 regimens), and a biomarker panel (plasma calcium, uric acid, citrate, and renal inflammation index).

## Usage

```r
library(mrgsolve)
mod <- mread("uri_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("uri_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg uri_qsp_model.dot -o uri_qsp_model.svg
```

## References

For detailed citations, see [uri_references.md](../../../urolithiasis/uri_references.md) (approximately 38 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
