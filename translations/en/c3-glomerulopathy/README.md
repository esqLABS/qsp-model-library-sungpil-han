# C3 Glomerulopathy (C3G: DDD & C3GN) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking genetic (CFH/
> CFI/CFB/C3/MCP variants, CFHR5 duplication) and acquired (C3 Nephritic
> Factor [C3NeF], anti-Factor H autoantibody) drivers of uncontrolled
> alternative-pathway (AP) C3-convertase (C3bBb) amplification to glomerular
> C3b deposition (mesangial/subendothelial in C3 glomerulonephritis [C3GN],
> intramembranous ribbon-like dense deposits in Dense Deposit Disease [DDD]),
> mesangial proliferation, podocyte injury, proteinuria, TGF-β-driven
> interstitial fibrosis, nephron loss, eGFR decline, ESKD and kidney
> transplant recurrence — coupled to upstream alternative-pathway inhibitor
> (iptacopan [Factor B], pegcetacoplan [C3/C3b], danicopan [Factor D]) and
> terminal-pathway inhibitor (eculizumab, ravulizumab [C5]) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`c3g_qsp_model.dot`](c3g_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`c3g_qsp_model.svg`](c3g_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`c3g_qsp_model.png`](c3g_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`c3g_mrgsolve_model.R`](../../../c3-glomerulopathy/c3g_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`c3g_shiny_app.R`](c3g_shiny_app.R) |
| 📚 References             | [`c3g_references.md`](c3g_references.md) |

---

## 1. Disease in One Paragraph

C3 Glomerulopathy (C3G) is a group of rare kidney diseases in which genetic
or acquired dysregulation of the complement **alternative pathway (AP)**
excessively stabilises the C3 convertase (C3bBb), leading to glomerular
deposition dominated by C3 degradation products. It splits into two
subtypes: **Dense Deposit Disease (DDD)**, histologically characterised by
ribbon-like intramembranous electron-dense deposits, and **C3
glomerulonephritis (C3GN)**, characterised by mesangial/subendothelial
deposits. The core pathogenesis is loss of regulation from **loss-of-function
variants in Factor H (CFH)/Factor I (CFI)**, convertase stabilisation from
**gain-of-function variants in Factor B (CFB)/C3**, or — most commonly —
acquired dysregulation by **C3 Nephritic Factor (C3NeF)**, an autoantibody
that binds the C3bBb convertase and blocks its Factor H-mediated decay. This
stabilised convertase continuously cleaves C3 in the fluid phase and on the
glomerular surface, depositing C3b and driving downstream C5-convertase
formation and terminal-pathway activation via the membrane attack complex
(C5b-9), which in turn causes mesangial proliferation, podocyte injury, and
proteinuria. The cumulative pathological glomerular injury accelerates
TGF-β-mediated tubulointerstitial fibrosis and nephron loss, so that roughly
half of patients reach end-stage kidney disease (ESKD) within 10 years, and
histological recurrence occurs in up to 50-90% after kidney
transplantation. Characteristic extrarenal manifestations can include ocular
**drusen** and **acquired partial lipodystrophy**. Standard therapy (RAAS
blockers, immunosuppressants) has had limited efficacy, but the oral
**Factor B inhibitor iptacopan** and the subcutaneous **C3/C3b inhibitor
pegcetacoplan**, which act directly on the upstream alternative pathway, have
recently produced significant reductions in proteinuria in phase 3 trials. A
mechanistically important point of differentiation is that **C5 inhibitors
(eculizumab, ravulizumab)**, which block only the terminal pathway, leave
upstream C3 deposition ongoing and so show heterogeneous responses.

## 2. Mechanistic Map Clusters (16 clusters, 120 nodes)

1. Genetic/acquired predisposition (CFH/CFI/CFB/C3/MCP variants, CFHR5
   duplication, monoclonal gammopathy, infection triggers, genetic panel
   testing)
2. Autoantibody-mediated mechanisms (C3NeF, anti-Factor H/B/C3b
   autoantibodies, plasmablast clones, persistent vs. transient titre)
3. Alternative-pathway amplification loop (C3 tick-over -> Factor B/D ->
   C3bBb convertase -> properdin stabilisation -> C3b deposition/C3a)
4. Complement regulator network (Factor H/I, MCP/CD46, DAF/CD55, CR1, net
   regulatory reserve)
5. Terminal pathway/membrane attack complex (C5 convertase -> C5a -> C5b-9
   MAC -> lytic/non-lytic signalling, sC5b-9)
6. Glomerular deposition histopathology (mesangial C3 deposition, MPGN
   pattern, C3GN vs. DDD EM classification, crescents, GBM double contours)
7. Glomerular cell injury (mesangial proliferation, endothelial activation,
   podocyte C5aR1/C3aR signalling, foot-process effacement, leukocyte
   infiltration)
8. Tubulointerstitial injury and fibrotic progression (RAAS activation,
   TGF-β1/CTGF, myofibroblasts, IF/TA, glomerulosclerosis, nephron-loss
   feedback)
9. Extrarenal clinical manifestations (drusen/retinopathy, acquired partial
   lipodystrophy, infection susceptibility)
10. Kidney transplantation and recurrence (persistent systemic AP
    dysregulation, allograft recurrence, graft failure)
11. Iptacopan PK/PD (oral Factor B inhibitor)
12. Pegcetacoplan PK/PD (C3/C3b inhibitor)
13. C5 inhibitor PK/PD (eculizumab/ravulizumab)
14. Danicopan PK/PD (oral Factor D inhibitor, combination therapy)
15. Biomarkers (serum C3/C4, Factor B/H/I, sC5b-9, C3NeF titre, AP function
    assays)
16. Clinical endpoints (complete/partial remission, proteinuria, eGFR
    trajectory, ESKD, transplant survival, histological activity)

## 3. mrgsolve Model (21 ODE Compartments)

* **Drug PK (10 compartments)** — iptacopan gut depot/plasma (2),
  pegcetacoplan subcutaneous depot/plasma (2), eculizumab central/peripheral
  (2, two-compartment monoclonal antibody), ravulizumab central/peripheral
  (2, FcRn recycling), danicopan gut depot/plasma (2).
* **Disease/PD (10 compartments)** — alternative-pathway convertase activity
  index (AP_ACTIVITY), serum C3 (C3_LEVEL), soluble terminal complex
  (SC5B9), glomerular deposition burden (GLOM_DEPOSIT), mesangial activity
  index (MESANGIAL), surviving podocyte fraction (PODO_FRAC), fibrosis
  index (FIBROSIS), residual nephron fraction (NEPHRON_FRAC), proteinuria
  (UPCR), eGFR.
* **Time tracking (1 compartment)** — FX_YEARS.
* Core design: **upstream (Factor B/D, C3) inhibitors** (iptacopan,
  pegcetacoplan, danicopan) lower AP_ACTIVITY itself, reducing both
  glomerular deposition (GLOM_DEPOSIT) and the terminal complex (SC5B9),
  whereas **terminal-pathway (C5) inhibitors** (eculizumab/ravulizumab) were
  implemented to block only SC5B9 via a c5_block term, leaving the upstream
  C3b deposition pathway intact — this is a deliberate mechanistic
  representation of why C5 inhibitors show heterogeneous and limited
  responses in real-world C3G treatment. Danicopan is designed as an
  add-on therapy on a C5-inhibitor background that partially suppresses
  residual (breakthrough) AP activity.

### 10 Scenarios

| # | Scenario | Calibration basis |
|---|---|---|
| 1 | Natural history - C3GN phenotype (SEVERITY=1.0) | Servais 2012 Kidney Int |
| 2 | Natural history - DDD phenotype, high-titre C3NeF (SEVERITY=1.35) | Smith 2019 Nat Rev Nephrol |
| 3 | Iptacopan 200 mg BID | Bomback 2025 Lancet (APPEAR-C3G) |
| 4 | Pegcetacoplan 1080 mg SC twice weekly | Dixon 2023 Kidney Int Rep |
| 5 | Eculizumab 900 mg IV q2w (off-label) | Bomback 2012 CJASN |
| 6 | Ravulizumab weight-based IV q8w | Lee 2019 Blood (PK bridging) |
| 7 | Danicopan 150 mg TID + eculizumab combination | Risitano 2021 Br J Haematol (extrapolated from PNH EVH) |
| 8 | CFH loss-of-function variant + iptacopan | Gale 2010 Lancet (CFHR5 nephropathy) |
| 9 | Post-transplant recurrence + pegcetacoplan prophylaxis | Zand 2014 J Am Soc Nephrol |
| 10 | Iptacopan discontinued after 2 years -> relapse | Extrapolated from the APPEAR-C3G extension-phase trial design |

## 4. Shiny Dashboard (8 Tabs)

1. **Patient profile** — phenotype (DDD/C3GN severity), C3NeF titre, CFH
   variant, age, simulation duration.
2. **PK** — plasma concentrations by drug (iptacopan/pegcetacoplan/
   eculizumab/ravulizumab/danicopan).
3. **Key PD metrics (complement system)** — AP activity index, serum C3, and
   sC5b-9.
4. **Clinical endpoints** — eGFR trajectory (with ESKD threshold line),
   UPCR, time to ESKD.
5. **Glomerular histology** — deposition burden and mesangial activity,
   podocyte fraction and fibrosis.
6. **Scenario comparison** — overlaid comparison across multiple scenarios
   with a summary table.
7. **Biomarkers** — serum C3, sC5b-9, AP activity, deposition burden,
   mesangial index.
8. **References** — the full reference list.

## 5. How to Run

```bash
# 1) Render the mechanistic map
dot -Tsvg c3g_qsp_model.dot -o c3g_qsp_model.svg
dot -Tpng -Gdpi=150 c3g_qsp_model.dot -o c3g_qsp_model.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("c3g_mrgsolve_model.R") %>% param(SEVERITY = 1.0, C3NEF_TITER = 1.0)
e_ipta <- ev(amt = 200, cmt = "IPTA_GUT", time = 0, ii = 12, addl = 2*365*2-1)
out <- mod %>% ev(e_ipta) %>% mrgsim(end = 24*365*10, delta = 24)  # 10-year follow-up
plot(out, c("AP_ACTIVITY", "C3_LEVEL", "UPCR", "EGFR"))

# 3) Run the Shiny dashboard
shiny::runApp("c3g_shiny_app.R")
```

## 6. Key Clinical Calibration Basis

| Endpoint | Comparator | Basis |
|---|---|---|
| Natural-history 10-year ESKD rate (~50%) | Registry cohort | Servais 2012 Kidney Int |
| Iptacopan 6-month proteinuria reduction (35.1% vs. placebo) | Randomised controlled phase 3 | Bomback 2025 Lancet (APPEAR-C3G) |
| Pegcetacoplan 48-week proteinuria reduction >=50% | Single-arm phase 2 | Dixon 2023 Kidney Int Rep |
| Eculizumab response heterogeneity (benefit only in sC5b-9 responders) | Open-label proof of concept | Bomback 2012 CJASN |
| Ravulizumab C5 saturation and PK bridging | PNH programme | Lee 2019 Blood |
| Danicopan add-on suppression of breakthrough activity (extrapolated from PNH EVH) | Randomised controlled phase 3 (ALPHA) | Risitano 2021 Br J Haematol |
| Post-transplant recurrence rate (up to 50-90%) | Retrospective cohort | Zand 2014 J Am Soc Nephrol |
| CFHR5 nephropathy autosomal-dominant founder gene duplication | Family cohort | Gale 2010 Lancet |

## 7. Model Validation Status

This container does not have an initial R/mrgsolve runtime environment
installed (no `Rscript`), so the mrgsolve model has been completed only to
the stage of **literature-based parameter design and self-review of the
code** (all 21 compartments confirmed to correspond 1:1 across `$CMT`/`$ODE`/
`$INIT`, and the logic deliberately separating the site of action of upstream
AP inhibitors from that of terminal (C5) inhibitors was reviewed) — the
values have not been verified by an actual compile-and-integrate run. The
`.dot` file was actually rendered with Graphviz `dot` (installed via
`apt-get install graphviz`) to generate and check the SVG/PNG (120 nodes, 16
clusters). The 32 PubMed references were individually verified by cross-
checking author/year/journal information against the actual PubMed and
journal pages via WebSearch (4 items not indexed with a PMID, such as
conference abstracts, are instead marked with a journal webpage link). It is
recommended to run the model as described in "How to Run" above, in an
environment with mrgsolve/R available, to confirm the numerical integration
results.
