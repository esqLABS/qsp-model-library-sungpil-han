# Lymphangioleiomyomatosis (LAM) — QSP Model

[![LAM QSP Model](lam_qsp_model.png)](lam_qsp_model.svg)

**Category**: Rare Lung Disease / Tumor Suppressor Gene Disorder
**Abbreviation**: LAM
**Directory**: [`lymphangioleiomyomatosis/`](.)

---

## Disease Overview

**Lymphangioleiomyomatosis (LAM)** is a rare, progressive cystic lung disease that occurs mainly in
women of reproductive age.
Loss-of-function mutations in the TSC1 or TSC2 gene hyperactivate mTORC1, causing LAM cells
(smooth-muscle-like neoplastic cells) to infiltrate and destroy the lung parenchyma, forming
thin-walled cysts in both lungs.

| Item | Detail |
|------|------|
| Prevalence | 1-9 per 100,000 women (approx. 3,000-6,000 worldwide) |
| Typical age of onset | Women of reproductive age (20s-40s) |
| Genetic background | TSC2 mutation (sporadic, 80%), TSC1/2 germline mutation (TSC-LAM) |
| Core pathway | TSC1/2 → Rheb-GTP → mTORC1 hyperactivation |
| Diagnostic biomarker | Serum VEGF-D >800 pg/mL (sensitivity ~73%, specificity ~100%) |
| Major comorbidities | Renal angiomyolipoma (AML, 50-60%), pneumothorax, chylothorax |
| First-line treatment | Sirolimus 2mg/day (ERS 2022, ATS 2017) |
| Prognosis | ~120 mL/yr annual FEV1 decline (untreated); stabilises with treatment |

---

## QSP Model Core Mechanisms (Key Mechanisms)

### 1. The TSC1/TSC2 Two-Hit Model
- **TSC-LAM**: a TSC1 or TSC2 germline mutation plus somatic LOH → complete loss of the TSC complex
- **Sporadic LAM**: a de novo somatic TSC2 mutation → loss of function of tuberin (the TSC2 protein)
- The TSC complex (hamartin-tuberin) is Rheb's GTPase-activating protein (GAP) — when it is lost,
  Rheb-GTP rises sharply

### 2. The mTOR Hyperactivation Pathway
```
Loss of TSC2 function
    ↓
Rheb-GTP up (~2.5x normal)
    ↓
mTORC1 hyperactivation (~4x normal)
    ↓
S6K1-pT389 up / 4E-BP1-P up
    ↓
Increased LAM cell proliferation, survival, invasion
```

### 3. LAM Cell Biology
- **Origin**: uterine/pelvic LAM cells → hematogenous/lymphatic spread → lung colonisation
- **Phenotype**: HMB-45+, α-SMA+, desmin+, ERα+, PR+ (smooth-muscle plus melanocytic features)
- **Invasion**: secretion of MMP-2, MMP-9, MMP-13 → ECM degradation → peribronchial cyst formation
- **VEGF-D**: serum VEGF-D secretion → VEGFR-3 → lymphangiogenesis → chylothorax/retroperitoneal
  lymphangioleiomyoma

### 4. Estrogen Regulation
- ERα expression → E2 activates the non-genomic PI3K/Akt pathway → Akt phosphorylation-mediated
  inhibition of TSC2 → further upregulation of mTORC1
- Progression slows after menopause, and disease tends to worsen during pregnancy → confirming
  estrogen dependence

### 5. Drug PK/PD (mTOR Inhibitors)

| Drug | Mechanism | Target Blood Concentration | Key Clinical Evidence |
|------|------|----------------|----------------|
| **Sirolimus** 2mg/day | FKBP12 binding → allosteric mTORC1 inhibition | 5-15 ng/mL (trough) | MILES (NEJM 2011) |
| **Everolimus** 10mg/day | FKBP12 binding → mTORC1 inhibition | 5-10 ng/mL | EXIST-2 (Lancet 2016) |

---

## ODE Model Structure (18 Compartments)

| No. | Compartment | Description |
|------|------|------|
| 1-3 | SIRO_GUT/C/P | Sirolimus PK (oral, 2-compartment) |
| 4-6 | EVER_GUT/C/P | Everolimus PK |
| 7 | RHEB_GTP | Rheb-GTP fraction |
| 8 | MTORC1 | mTORC1 activity (normalised) |
| 9 | S6K1_P | S6K1-pT389 phosphorylation |
| 10 | EBPP1 | 4E-BP1 phosphorylation |
| 11 | LAM_CELLS | LAM cell burden |
| 12 | VEGFD | Serum VEGF-D (pg/mL) |
| 13 | MMP_ACT | MMP activity (normalised) |
| 14 | ESTROGEN | Estrogen level |
| 15 | CYST_VOL | Pulmonary cyst volume (%) |
| 16 | FEV1_PCT | FEV1 (% predicted) |
| 17 | DLCO_PCT | DLCO (% predicted) |
| 18 | AML_VOL | Renal AML volume (mL) |

---

## Treatment Scenarios (5)

1. **Untreated** — natural history (~120 mL/yr FEV1 decline)
2. **Sirolimus 2mg/day** — calibrated to the MILES trial
3. **Everolimus 10mg/day** — calibrated to EXIST-2 (50% AML reduction)
4. **Sirolimus, discontinued after 12 months** — simulates FEV1 re-decline after MILES-style
   discontinuation
5. **Everolimus + GnRH agonist** — combined mTOR inhibition and estrogen blockade

---

## Clinical Calibration

| Clinical Trial | Key Result | Model Implementation |
|---------|---------|---------|
| **MILES** (McCormack, NEJM 2011) | Sirolimus: FEV1 +153mL vs. control, VEGF-D down ~30% | FEV1 stabilisation, VEGF-D reduction PD |
| **Bissler et al.** (NEJM 2008) | Sirolimus: AML -47%, regrows after discontinuation | AML ODE compartment |
| **Johnson et al.** (NEJM 2010) | Natural history: FEV1 decline ~117 mL/yr | kFEV1_decline parameter |
| **EXIST-2** (Kingswood, Lancet 2016) | Everolimus: AML reduction >50% | kAML_shrink parameter |
| **Young et al.** (Ann Intern Med 2011) | VEGF-D >800 pg/mL: sensitivity 73%, specificity 100% | Initial value VEGFD_LAM=1500 |

---

## Files

| File | Description |
|------|------|
| [`lam_qsp_model.dot`](lam_qsp_model.dot) | Graphviz mechanistic map source |
| [`lam_qsp_model.svg`](lam_qsp_model.svg) | Vector mechanistic map |
| [`lam_qsp_model.png`](lam_qsp_model.png) | Raster mechanistic map (150 dpi) |
| [`lam_mrgsolve_model.R`](lam_mrgsolve_model.R) | 18-compartment ODE PK/PD model |
| [`lam_shiny_app.R`](lam_shiny_app.R) | 6-tab interactive Shiny dashboard |
| [`lam_references.md`](lam_references.md) | 50 PubMed references |

---

## Shiny App Tab Structure (6 Tabs)

| Tab | Key Contents |
|----|---------|
| 1. Patient Profile | Patient profile, severity gauge, disease overview |
| 2. Drug PK | Concentration-time curves, steady-state PK, PK parameter table |
| 3. mTOR Pathway (PD) | mTORC1 inhibition rate, S6K1-pT389, 4E-BP1, Rheb-GTP |
| 4. Clinical Endpoints | FEV1, DLCO, cyst volume, estimated 6-minute walk distance |
| 5. Scenario Comparison | Comparison of 5 treatment scenarios, 12/24-month outcome table |
| 6. Biomarker Dashboard | VEGF-D, S6K1, MMP, AML volume panel |

---

## Mechanistic Map Clusters (14 Subgraphs, 120+ Nodes)

1. Genetic Basis (TSC1/TSC2)
2. LAM Cell Origin & Phenotype
3. Upstream mTOR Signaling Inputs (PI3K/Akt, AMPK, HIF-1α)
4. TSC/Rheb/mTOR Core Axis
5. mTORC1 Downstream Substrates (S6K1, 4E-BP1, autophagy, TFEB)
6. LAM Cell Biology (proliferation, migration, invasion, MMP)
7. Lung Pathology (cysts, ECM degradation, airflow obstruction)
8. Lymphatic & Extra-Pulmonary Involvement (VEGF-D, AML, chylothorax)
9. Hormonal Regulation (E2, ERα, GnRH)
10. Drug PK (Sirolimus & Everolimus)
11. Drug PD (mTOR inhibition, PD biomarkers)
12. Adverse Effects (stomatitis, pneumonitis, hyperlipidemia)
13. Disease Biomarkers (VEGF-D, S6K1, CT score, spirometry)
14. Clinical Endpoints & Management (FEV1, transplant, guidelines)

---

## Reference Summary

- McCormack FX et al. **MILES Trial**. *N Engl J Med* 2011;364:1595-1606.
- Kingswood JC et al. **EXIST-2**. *Lancet* 2016;387:1629-1638.
- Young LR et al. **VEGF-D Biomarker**. *Ann Intern Med* 2011;154:743-751.
- Johnson SR et al. **Natural History**. *N Engl J Med* 2010;363:950-959.
- Gupta N et al. **ATS/JRS Guidelines**. *Am J Respir Crit Care Med* 2017;196:1337-1348.

Full reference list, 50 entries: [`lam_references.md`](lam_references.md)
