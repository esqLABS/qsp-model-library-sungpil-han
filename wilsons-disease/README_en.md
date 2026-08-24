# Wilson's Disease (WD) QSP Model

> **Directory**: `wilsons-disease/` | **Abbreviation**: WD | **Date**: 2026-06-24

[![WD QSP mechanistic map](wd_qsp_model.png)](wd_qsp_model.svg)

---

## Disease Overview

**Wilson's disease (WD)** is a disorder of copper metabolism caused by
autosomal recessive mutations in the *ATP7B* gene. ATP7B is a P-type
Cu-ATPase that exports copper into bile from hepatocytes and loads copper
onto ceruloplasmin; when it loses function, copper accumulates in the
liver, brain, cornea, and kidney, causing progressive hepatic and
neurological damage.

| Feature | Value |
|------|-----|
| Prevalence | 1/30,000 (carrier frequency 1/90) |
| Age of onset | 5–35 years (reported up to 45) |
| Inheritance | Autosomal recessive (AR) |
| Gene location | 13q14.3 (ATP7B, 21 exons) |
| Common mutations | p.His1069Gln (Europe, 35%), p.Arg778Leu (Asia, 20%) |
| Leipzig diagnostic score | ≥4 = confirmed diagnosis |

---

## Pathophysiology Summary

```
Copper intake -> GI absorption (CTR1/ATP7A) -> portal vein -> hepatocyte uptake
                                              |
                                    ATP7B function x
                                              |
                              +-- biliary excretion down --+
                              |                            |
                    hepatic Cu accumulation (>250 ug/g)  ceruloplasmin down
                              |
                    Cu excess -> Fenton reaction -> ROS -> oxidative stress
                              |
                    MT saturation -> NCBC increase -> systemic distribution
                              |
                  +-----------+------------+
                  |           |            |
              brain Cu    renal toxicity  corneal Cu
          (basal ganglia   Fanconi        KF ring
           preference)
                  |
          neuropsychiatric symptoms (UWDRS up)
```

---

## Comparison of Treatment Drugs

| Drug | Mechanism | Main indication | NCBC reduction | Main adverse effects |
|------|------|-----------|---------|-----------|
| **D-Penicillamine** | Copper chelation → urinary excretion | Hepatic WD, first-line | ~60–70% | Nephrotoxicity, SLE-like syndrome, neurological worsening (50%) |
| **Trientine** | Copper chelation (2nd-line to DPA) | On DPA adverse effects | ~50–60% | Gastrointestinal discomfort (mild) |
| **Zinc Acetate** | Induces intestinal MT → blocks Cu absorption | Maintenance therapy, pregnancy, children | ~30–40% | Gastrointestinal irritation, Fe interaction |
| **ALXN1840 (TTM)** | Forms a TTM-Cu-albumin tripartite complex | Neurological WD (first-line), hepatic | **~98%** | Mild (ATLAS trial) |

---

## QSP Model Files

| Component | File | Specification |
|---------|------|-----|
| Mechanistic map | [`wd_qsp_model.dot`](wd_qsp_model.dot) | **119 nodes, 11 clusters** |
| SVG | [`wd_qsp_model.svg`](wd_qsp_model.svg) | Vector format |
| PNG | [`wd_qsp_model.png`](wd_qsp_model.png) | 150 DPI |
| mrgsolve ODE | [`wd_mrgsolve_model.R`](wd_mrgsolve_model.R) | **24-compartment ODE**, **8 treatment scenarios** |
| Shiny app | [`wd_shiny_app_en.R`](wd_shiny_app_en.R) | **8-tab** interactive dashboard |
| References | [`wd_references_en.md`](wd_references_en.md) | **60 references** (13 sections) |

---

## Mechanistic Map Clusters (11)

| # | Cluster | Key nodes |
|---|---------|---------|
| 1 | GI copper absorption | CTR1, DMT1, MT_gut, ATP7A, portal vein |
| 2 | Hepatocyte copper chaperones | ATOX1, CCS, COX17, SOD1, MT_hepatic |
| 3 | ATP7B & biliary excretion | ATP7B_WT/Mutant, TGN, bile vesicles, ceruloplasmin |
| 4 | Systemic copper & biomarkers | Cp_serum, NCBC, 24h urinary Cu, Leipzig score |
| 5 | Hepatic pathology & fibrosis | Cu-ROS, Kupffer cells, HSC, TGF-β, fibrosis, ALT/AST |
| 6 | Brain copper & neuropsychiatry | Basal ganglia, dopaminergic neurons, tremor/dystonia, UWDRS |
| 7 | Multi-organ toxicity | KF ring, Fanconi syndrome, haemolytic anaemia |
| 8 | DPA PK/PD | Absorption → chelation → urinary excretion, the paradox of neurological worsening |
| 9 | Zinc PK/PD | Induction of intestinal MT, blockade of Cu absorption |
| 10 | Trientine & ALXN1840 PK/PD | Tripartite complex, faecal excretion |
| 11 | Clinical outcomes | Leipzig score, treatment response, liver transplant indications |

---

## mrgsolve ODE Model Compartments (24)

```
[Drug PK -- 8 compartments]
GUT_DPA -> CENT_DPA           D-Penicillamine, 2-compartment PK
GUT_ZN  -> CENT_ZN            Zinc, 2-compartment PK
GUT_TRI -> CENT_TRI           Trientine, 2-compartment PK
GUT_TTM -> CENT_TTM           ALXN1840, 2-compartment PK

[Copper kinetics -- 7 compartments]
CU_GI                         GI copper (absorption pool)
CU_HEP                        Hepatic copper (the core accumulation site)
MT_HEP                        Hepatic metallothionein-bound copper
CU_NCBC                       NCBC (free copper, the toxic fraction)
CP_SERUM                      Serum ceruloplasmin
CU_URINE                      Urinary copper excretion rate
CU_BRAIN / CU_KIDNEY / CU_CORNEA  Organ copper

[Hepatic pathology -- 3 compartments]
ROS_HEP                       Hepatic ROS index
ALT_SERUM                     Serum ALT
FIBROSIS                      Fibrosis score (Metavir F0-F4)

[Neurological -- 1 compartment]
NEURODEGENERATION              Neurodegeneration index
```

---

## Treatment Scenarios (8)

| Scenario | Drug | Clinical basis | Key result |
|---------|------|---------|---------|
| S1 | Untreated WD | Natural history | Hepatic Cu ↑↑, fibrosis progression, neurodegeneration |
| S2 | DPA 500 mg TID | Walshe 1956 | NCBC ↓60%, possible early neurological worsening |
| S3 | Zinc 50 mg TID | Brewer 1998 | Absorption blocked, NCBC ↓30–40% |
| S4 | Trientine 500 mg TID | Weiss 2013 | NCBC ↓50%, fewer adverse effects than DPA |
| S5 | ALXN1840 15 mg QD | ATLAS 2022 | **NCBC ↓98%** |
| S6 | DPA → zinc switch (after 1 year) | AASLD guideline | Initial chelation followed by a maintenance strategy |
| S7 | ALXN1840 + trientine combination | Hypothetical exploration | Maximal NCBC reduction |
| S8 | Normal WT control | Reference | Normal Cp, NCBC<10, normal ALT |

---

## Shiny App Tab Structure (8 Tabs)

1. **Patient profile**: Leipzig score calculator, ATP7B mutation distribution, genotype-phenotype
2. **Drug PK**: simulated blood concentrations of DPA/zinc/trientine/ALXN1840
3. **Copper kinetics**: time series of serum copper, NCBC, ceruloplasmin, organ distribution
4. **Hepatic outcomes**: ALT, fibrosis, ROS trajectory, liver-transplant risk index
5. **Neurological/ophthalmic outcomes**: brain copper, UWDRS, KF ring trajectory, DPA's paradoxical worsening
6. **Scenario comparison**: side-by-side comparison of the 8 scenarios against the ATLAS trial
7. **Biomarker explorer**: correlations among copper biomarkers, diagnostic panel
8. **Model information**: parameter sources, clinical trial summary

---

## Key References

- Członkowska A, et al. *Nat Rev Dis Primers* 2018 (**comprehensive disease review**)
- Bandmann O, et al. *Lancet Neurol* 2015 (**neurological mechanisms**)
- Schilsky ML, et al. *NEJM Evid* 2022 (**ATLAS trial — ALXN1840**)
- EASL. *J Hepatol* 2012 (**diagnostic and treatment guideline**)
- Brewer GJ, et al. *J Lab Clin Med* 1998 (**15-year zinc follow-up**)
