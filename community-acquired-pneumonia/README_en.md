# Community-Acquired Pneumonia (CAP) — QSP Model

<p align="center">
  <a href="cap_qsp_model_en.svg"><img src="cap_qsp_model_en.png" width="720" alt="CAP QSP mechanistic map"></a>
</p>

A Quantitative Systems Pharmacology (QSP) model of community-acquired
pneumonia, built around *Streptococcus pneumoniae* as the archetype. It ties
bacterial burden, **PAMP release that depends on how the bacteria are
killed**, host inflammatory amplification, alveolar-capillary barrier
destruction and intrapulmonary shunt, systemic organ failure, and the PK/PD of
ceftriaxone · amoxicillin · azithromycin · levofloxacin · hydrocortisone into
a single system of differential equations.

| File | Contents |
|---|---|
| [`cap_qsp_model_en.dot`](cap_qsp_model_en.dot) | Mechanistic map source (14 clusters, 150+ nodes, 270+ edges) |
| [`cap_qsp_model_en.svg`](cap_qsp_model_en.svg) / [`cap_qsp_model_en.png`](cap_qsp_model_en.png) | Rendered map (Graphviz `dot`, PNG 150 dpi) |
| [`cap_mrgsolve_model.R`](cap_mrgsolve_model.R) | mrgsolve model — **38 ODE compartments**, 9 predefined scenarios |
| [`cap_shiny_app.R`](cap_shiny_app.R) | Shiny dashboard — 8 tabs |
| [`cap_references_en.md`](cap_references_en.md) | 71 PubMed-verified references (every PMID cross-checked) |

---

## What this model is trying to say (The three structural commitments)

The textbook CAP model ends at "the antibiotic kills the bacteria, so the
patient improves." That single picture cannot simultaneously satisfy the six
facts the literature demands. So three structural choices were made.

### ① fT>MIC and fAUC/MIC are not inputs but **integrated outputs**

Nowhere in the model is it written that "beta-lactams are time-dependent,
quinolones are concentration-dependent." Both drugs use the **same** saturating
Emax function — against free ELF concentration, with EC50 proportional to the
strain's MIC. The only difference is the Hill coefficient (1.5 for
ceftriaxone, 2.6 for levofloxacin). `TAM` (fT>MIC) and `AUCF` (fAUC/MIC) are
**cumulative compartments** that integrate that curve over time, so both
metrics are observations, not parameters.

### ② Bactericidal killing is a **source** of PAMPs, not a sink

Beta-lactams inhibit PBPs and trigger the LytA autolysin, bursting the
bacteria. This dumps cell-wall material (LTA · peptidoglycan) and pneumolysin
into the alveolus all at once. Non-lytic killing (macrolides, and to a lesser
degree quinolones) releases far less.

```
dPAMP/dt = YLIVE·(growth of live bacteria) + YHOST·(phagocytic killing)
         + (YLYS_BL·k_βlactam + YLYS_FQ·k_quinolone + YLYS_MAC·k_macrolide)·BE
         − KPAMP·PAMP
```

`YLYS_BL = 1.00`, `YLYS_FQ = 0.35`, `YLYS_MAC = 0.12`.
This single term produces **the inflammatory surge and transient clinical
worsening seen 6-24 hours after the first dose**, and simultaneously explains
the **macrolide paradox** — why a survival benefit from adding a macrolide is
observed in severe pneumococcal CAP even when the strain is macrolide
**resistant**: because the benefit flows not through the bactericidal pathway
but through the toxin/cytokine pathway (`MAC_PLY`, `MAC_ANTI`, `MAC_MIG`).
Setting `MIC_AZI = 64` drives `kazi ≈ 0`, but the three immunomodulatory terms
remain untouched.

This structure is falsifiable: setting `YLYS_BL = 0` erases both the initial
surge and most of the macrolide benefit.

### ③ Corticosteroids have two arms, and **the sign of the net effect is conditional**

A single effect compartment, `HCE`, splits into two directions.

| Arm | Term | Direction |
|---|---|---|
| NF-κB transrepression | `GC_ANTI` (`IMAX_GC` 0.72) | cytokine amplifier · barrier damage · shunt · vasoplegia ↓ (**benefit**) |
| Opsonophagocytosis suppression | `GC_PHAG` (`IMAX_GCP` 0.45) | host bactericidal capacity ↓ (**harm**) |

So the net effect depends on whether the antibiotic has already cleared the
bacteria enough that host-mediated killing has become **unnecessary**. With a
susceptible strain and an effective beta-lactam, steroids are protective (the
CAPE-COD direction); with a resistant strain, or viral pneumonia where the
antibiotic has no target to hit, **the same dose is harmful**.

Scenarios 6 and 7 hold the steroid dose, the host, and every other parameter
identical, and vary only the MIC. Yet the sign of `Mortality_pct` flips.
Setting `IMAX_GCP = 0` makes this claim disappear entirely, so it too is
falsifiable.

---

## Mechanistic Map (14 clusters)

1. Host susceptibility · portal of entry — age, comorbidity, aspiration, mucociliary clearance, secretory IgA, vaccination
2. Pathogen · virulence factors — capsule, pneumolysin, LytA, PspA/PspC, NanA, atypical organisms, viruses
3. Alveolar epithelium · barrier — AT1/AT2, tight junctions, surfactant, ENaC/Na-K-ATPase, hyaline membrane
4. Innate immune recognition — TLR2/4/9, NOD2, NLRP3, cGAS-STING, MyD88-NF-κB, gasdermin-D
5. Cellular effectors — alveolar macrophages, neutrophils, NETs, MPO/elastase, monocyte M1/M2, efferocytosis
6. Cytokine network — TNF-α, IL-1β, IL-6, IL-8, CXCL1/2, CCL2, IL-17A, IL-10, TGF-β
7. Acute-phase response · biomarkers — CRP, PCT, MR-proADM, lymphopenia
8. Gas exchange · respiratory mechanics — V/Q mismatch, shunt, HPV, dead space, compliance, PaO₂/FiO₂
9. Systemic spread · organ failure — bacteremia, endothelial activation, vasoplegia, lactate, AKI, DIC, myocardial injury, SOFA
10. Antibiotic PK — absorption, protein binding, ELF penetration, lung-tissue accumulation, renal/biliary clearance
11. Antibiotic PD · resistance — PBP2x/2b, ermB/mefA, GyrA/ParC, efflux pumps, persisters, **mode of lysis**
12. Adjunctive therapy — steroids (GR · transcriptional repression · immunosuppression), oxygen/HFNC/NIV, fluids/vasopressors, vaccination
13. Resolution · repair — efferocytosis, SPM, AT2→AT1, matrix remodeling, organizing pneumonia
14. Clinical outcomes — resolution of fever, clinical stability (Halm), treatment failure, ICU transfer, length of stay, 28-day mortality, relapse

---

## mrgsolve Model Structure (38 ODEs)

| Compartment group | Count | State variables |
|---|---|---|
| Pharmacokinetics | 11 | `CEFC` `CEFP` `AMXD` `AMXC` `AZID` `AZIC` `AZIT` `LVXD` `LVXC` `HCC` `HCE` |
| Pathogen | 5 | `BE` alveolar extracellular bacteria · `BP` persister/intracellular bacteria · `BB` bacteremia · `PAMP` · `VIR` |
| Immune | 9 | `AM` `NB` `NL` `TNF` `IL1` `IL6` `IL8` `IL10` `SPM` |
| Biomarkers | 3 | `CRP` `PCT` `LYM` |
| Lung | 4 | `PERM` `EDEM` `SURF` `FIBP` |
| Systemic | 3 | `MAPD` `LACT` `AKID` |
| Accumulators | 3 | `TAM` (fT>MIC) · `AUCF` (fAUC) · `CHZ` (cumulative mortality hazard) |

The time unit is **hours** — resolving exposure changes within a dosing
interval and the post-dose PAMP surge is the whole point of the model. The
default simulation window is 336 h (14 days).

Output metrics: `PaO2_FiO2`, `SpO2_pct`, `Shunt_fraction`, `Temperature`,
`Resp_rate`, `Heart_rate`, `MAP_mmHg`, `Lactate`, `SOFA_score`, `CURB65`,
`Clinically_stable` (implemented exactly per the Halm criteria),
`Mortality_pct = 100·(1 − e^(−CHZ))`.

### 9 scenarios

| # | Scenario | What it shows |
|---|---|---|
| 1 | Untreated natural history | Host-only clearance → a 7-9 day "crisis," high cumulative risk |
| 2 | Amoxicillin 1 g q8h (mild, outpatient) | Where an oral beta-lactam alone is sufficient |
| 3 | Ceftriaxone 2 g q24h (ward) | Reference monotherapy, fT>MIC ~100% |
| 4 | Ceftriaxone + azithromycin | Bacterial curve almost identical, but PAMP·IL-8·shunt are lower |
| 5 | Levofloxacin 750 mg qd | Same bacteriological endpoint, a different PK/PD route (fAUC/MIC ~70) |
| 6 | Severe CAP + hydrocortisone (CAPE-COD) | Steroid in the **protective** direction |
| 7 | Same steroid + a penicillin-non-susceptible strain | Steroid in the **harmful** direction — sign reversal |
| 8 | Ceftriaxone, 5-day short course | Safety of stopping after stabilisation; relapse is governed by `BP` |
| 9 | Antibiotic delayed 18 hours (severe) | Risk accrued during the delay is not recovered afterward |

### Sensitivity dials worth turning

| Parameter | Manipulation | Result |
|---|---|---|
| `YLYS_BL` | 1.0 → 0.2 | The post-dose inflammatory surge disappears, most of the macrolide benefit is lost |
| `IMAX_GCP` | 0.45 → 0 | Removes the harmful steroid arm → scenario 7 is no longer harmful |
| `MIC_AZI` | 0.12 → 64 | ermB resistance — bactericidal action dies, immunomodulation survives |
| `KPERS` | 0.010 → 0 | Relapse after short-course therapy disappears entirely |
| `CRCL` | 90 → 25 | Levofloxacin accumulates; ceftriaxone, cleared biliarily, is nearly unchanged |
| `VIRAL` | 0 → 1 | PCT stays low while CRP rises — the two markers dissociate |

---

## Usage

```r
# 1) Simulate the model
library(mrgsolve); library(dplyr); library(ggplot2)
mod <- mread_cache("cap_mrgsolve_model.R")

# Ceftriaxone 2 g IV q24h × 7 days + azithromycin 500 mg × 5 days, first dose at hour 6
ev_abx <- ev(amt = 2000, cmt = "CEFC", time = 6, ii = 24, addl = 6) +
          ev(amt = 500,  cmt = "AZID", time = 6, ii = 24, addl = 4)

out <- mod |> mrgsim(events = ev_abx, end = 336, delta = 0.25) |> as_tibble()
ggplot(out, aes(time/24, Bacteria_log10)) + geom_line()

# 2) Batch-run all 9 scenarios — see CAP_simulate_scenarios() in the comments at
#    the bottom of the model file

# 3) Shiny dashboard
shiny::runApp("cap_shiny_app.R")
```

Re-rendering the map:

```bash
dot -Tsvg cap_qsp_model_en.dot -o cap_qsp_model_en.svg
dot -Tpng -Gdpi=150 cap_qsp_model_en.dot -o cap_qsp_model_en.png
```

---

## Calibration Anchors

- `KPHAG_AM`/`KPHAG_N`/`KG` were tuned so the untreated natural history
  reproduces the 7-9 day "crisis" of pre-antibiotic-era pneumococcal pneumonia
- Ceftriaxone 2 g q24h: Cmax ~150 mg/L, trough ~10 mg/L, fu ~0.10, ELF/plasma
  free fraction ~0.4 → fT>MIC ~100% against an MIC of 0.25 mg/L
- Levofloxacin 750 mg qd: AUC24 ~100 mg·h/L, fu 0.70, ELF/plasma ~1.2 →
  fAUC/MIC ~70
- Azithromycin: serum Cmax ~0.4 mg/L, yet ELF 1-3 mg/L (the reason the deep
  lung compartment is modeled separately)
- Halm 1998: median time to clinical stability of ~3 days with effective
  treatment
- CAPE-COD (Dequin 2023): 28-day mortality in severe CAP of 6.2% vs 11.9%

---

## ⚠️ Disclaimer

This model is a **qualitative · semi-quantitative QSP model for educational
and research purposes**. It was built from published literature but has not
been fitted · validated against patient-level data, and **must not be used
for clinical decision-making · prescribing · regulatory submission.**
