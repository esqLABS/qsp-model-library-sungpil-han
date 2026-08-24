# Autoimmune Encephalitis (AIE) — QSP Model

> **QSP Disease Model Library** · A quantitative systems pharmacology (QSP) model automatically generated via the Claude Code Routine.
> Parent library → [../README_en.md](../README_en.md) · Category: Neurology

[![AIE QSP Model](aie_qsp_model.png)](aie_qsp_model.svg)

## Overview

Autoimmune encephalitis (AIE) is a group of autoimmune disorders in which antibodies directed against neuronal cell-surface antigens (chiefly the NMDA receptor, LGI1, CASPR2, and GABA-B receptor) directly impair brain function. Anti-NMDAR encephalitis is the most common form, occurring predominantly in young women and sometimes accompanied by an ovarian teratoma. In the pathogenesis, autoantibodies produced by B cells cross the blood-brain barrier (BBB) and bind the GluN1 subunit of the NMDA receptor in the hippocampus and cerebral cortex, driving receptor internalisation and downregulation and thereby causing synaptic transmission failure. Clinically, the disease follows a characteristic progression from neuropsychiatric symptoms (psychosis, cognitive impairment) → seizures → movement disorder → decreased consciousness. When the response to first-line immunotherapy (steroids, IVIG, plasma exchange) is inadequate, second-line therapy (rituximab, cyclophosphamide) is initiated.

## Key Pathways

| Pathway | Key molecule/mechanism | Clinical outcome |
|------|----------------|-----------|
| Antibody production | Germinal-centre B cells → plasma cells → anti-NMDAR IgG | Elevated CSF antibody titre |
| BBB permeation | Inflammation → tight-junction damage → IgG influx into CSF | Increased CNS antibody concentration |
| NMDAR internalisation | Antibody-GluN1 binding → clathrin-mediated internalisation | Decreased synaptic NMDAR density |
| Microglial activation | NMDAR loss → microglial activation → IL-6/TNF-α release | Persistent neuroinflammation |
| Dopamine imbalance | NMDAR↓ → dopaminergic disinhibition → dopamine excess | Psychotic symptoms |
| Excitotoxicity | Excess glutamate → Ca²⁺ influx | Cognitive impairment |
| Seizure threshold | NMDAR↓ → disrupted cortical excitation/inhibition balance | Epileptic seizures |

## Drug Targets

- **Methylprednisolone**: first-line immunotherapy; stabilises the BBB and suppresses inflammation; high-dose pulse followed by oral tapering
- **IVIG (intravenous immunoglobulin)**: neutralises antibodies and saturates Fc receptors; infused over 3–5 days; can be combined with steroids
- **Plasmapheresis**: rapidly removes pathogenic antibodies; a preferred option in life-threatening cases
- **Rituximab (anti-CD20)**: depletes B cells → suppresses antibody production; the standard second-line therapy and relapse-prevention agent
- **Tocilizumab (anti-IL-6R)**: growing evidence for suppressing IL-6-mediated BBB damage in refractory or relapsing AIE
- **Tumour removal**: where a teratoma is present, tumour resection is central to treatment as it removes the immune-stimulating source

## Model Files

| File | Description |
|------|------|
| [aie_qsp_model.dot](aie_qsp_model.dot) | Graphviz mechanistic map source (approximately 183 nodes / 12 clusters) |
| [aie_qsp_model.svg](aie_qsp_model.svg) | SVG vector image (scalable) |
| [aie_qsp_model.png](aie_qsp_model.png) | PNG image (150 dpi) |
| [aie_mrgsolve_model.R](aie_mrgsolve_model.R) | mrgsolve ODE model (approximately 25 compartments / multiple treatment scenarios) |
| [aie_shiny_app.R](aie_shiny_app.R) | Shiny dashboard |
| [aie_references.md](aie_references.md) | References (approximately 65 articles, PubMed links) |

## mrgsolve Model (ODE Model)

- **Compartment structure**: B cells (germinal-centre, plasma, long-lived plasma, and memory cells) + serum/CSF anti-NMDAR antibody + BBB permeability + microglia/IL-6/GFAP + NMDAR expression, glutamate, cognitive function, psychosis, and seizure PD compartments + IVIG (2 compartments) + methylprednisolone (2 compartments) + rituximab (2 compartments) + tocilizumab (2 compartments) PK compartments
- **Key treatment scenarios**: ① no treatment (natural course), ② steroid monotherapy, ③ IVIG monotherapy, ④ steroid + IVIG combination (first-line standard), ⑤ plasma exchange + steroid, ⑥ second-line therapy (rituximab), ⑦ tocilizumab (refractory disease)
- **Calibration/evidence**: parameters referenced from the Titulaer et al. Lancet Neurol 2013 outcome cohort, the Dalmau et al. original description, and rituximab case series

## Shiny Dashboard

The dashboard comprises a patient profile tab (antibody titre, presence of a teratoma, initial severity), an immunotherapy PK and B-cell dynamics tab, an NMDAR expression and CSF antibody concentration tab, a clinical endpoints tab (mRS, cognition, seizures), a treatment strategy comparison tab (first-line vs. second-line immunotherapy), and a neuroinflammation biomarker tab (GFAP, IL-6, antibody titre).

## Usage

```r
library(mrgsolve)
mod <- mread("aie_mrgsolve_model.R")
out <- mrgsim(mod, end = 365)
plot(out)
# Shiny dashboard:
shiny::runApp("aie_shiny_app.R")
```
```bash
# Render the mechanistic map
dot -Tsvg aie_qsp_model.dot -o aie_qsp_model.svg
```

## References

For detailed citations, see [aie_references.md](aie_references.md) (approximately 65 articles).

---
*This model is a qualitative/semi-quantitative QSP model intended for education and research, and must not be used directly for clinical decision-making.*
