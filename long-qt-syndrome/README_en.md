# Congenital Long QT Syndrome (LQTS) — QSP Model

> A cardiac channelopathy causing delayed ventricular repolarisation (QTc
> prolongation). Depending on the causative gene — KCNQ1 (IKs, LQT1), KCNH2/hERG
> (IKr, LQT2), SCN5A (late-INa, LQT3) — patients show different triggers
> (exercise/auditory stimuli/sleep) and different beta-blocker responsiveness, and
> can progress via early afterdepolarisation (EAD) and repolarisation dispersion to
> Torsades de Pointes (TdP) -> ventricular fibrillation (VF) -> sudden cardiac
> death (SCD). Acquired (drug-induced) LQTS and the CiPA (Comprehensive in vitro
> Proarrhythmia Assay) multi-ion-channel-block paradigm are also modelled.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`lqts_qsp_model_en.dot`](lqts_qsp_model_en.dot) |
| 🖼️ Map (SVG)             | [`lqts_qsp_model_en.svg`](lqts_qsp_model_en.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`lqts_qsp_model_en.png`](lqts_qsp_model_en.png) |
| ⚙️ mrgsolve ODE model     | [`lqts_mrgsolve_model_en.R`](lqts_mrgsolve_model_en.R) |
| 📊 Shiny dashboard        | [`lqts_shiny_app.R`](lqts_shiny_app.R) |
| 📚 References             | [`lqts_references.md`](lqts_references.md) |

---

## 1. Background & Epidemiology

The prevalence of congenital LQTS is estimated at about 1 in 2,000 (Schwartz 2009
*Circulation*). The great majority follow autosomal dominant inheritance
(Romano-Ward syndrome), while the autosomal recessive form caused by biallelic
KCNQ1/KCNE1 mutations (Jervell-Lange-Nielsen syndrome) is accompanied by congenital
sensorineural deafness. By causative gene, cases are classified as LQT1 (KCNQ1,
~35-40%), LQT2 (KCNH2/hERG, ~25-30%), and LQT3 (SCN5A, ~5-10%), with the remainder
being rare subtypes such as ANK2 (LQT4), CACNA1C (LQT8, Timothy syndrome),
CALM1-3 (LQT14-16), TRDN, and AKAP9. According to data from international LQTS
registries, the incidence of cardiac events (syncope, cardiac arrest, sudden death)
in untreated high-risk groups varies greatly by genotype, QTc duration, and sex
(Goldenberg 2008 *Circulation*; Sauer 2007 *JACC*), with QTc > 500 ms, female sex in
LQT2, and a history of syncope being the leading high-risk factors. Drug-induced
(acquired) LQTS arises from a range of non-cardiac drugs that block the hERG
channel (antibiotics, antipsychotics, antihistamines, etc.), and the CiPA
initiative assesses arrhythmia risk using a qNet-like index that reflects the
balance of multiple ion channels (hERG/ICaL/late-INa) rather than hERG block alone.

## 2. Pathophysiology

* **Genotype-specific channel defects**: LQT1 is a loss-of-function of KCNQ1
  (+ the KCNE1 beta-subunit) that reduces IKs (the slow delayed-rectifier K+
  current); in particular, the "safety-valve" mechanism that normally augments IKs
  under sympathetic activation fails to engage, so QTc paradoxically prolongs
  further during exercise or swimming. LQT2 arises from trafficking or pore defects
  in KCNH2 (hERG) that reduce IKr (the rapid delayed-rectifier K+ current), and is
  vulnerable to auditory stimuli, emotion, and the postpartum period. LQT3
  involves incomplete inactivation of SCN5A, increasing persistent/late INa and
  prolonging phase 2 of the action potential, with risk rising at rest, during
  sleep, and with bradycardia.
* **Action potential and repolarisation reserve**: reduced IKs/IKr or increased
  late-INa all reduce repolarisation reserve and prolong action potential duration
  (APD). When the APD gradient between epicardium, M cells, and endocardium widens,
  transmural dispersion of repolarisation (TDR) increases, and reactivation of the
  ICaL window current during the prolonged phase 2/3 produces early
  afterdepolarisations (EADs). If an EAD propagates successfully, a functional
  re-entrant circuit (phase-2 reentry) can form and progress to Torsades de
  Pointes.
* **Autonomic triggers**: sympathetic activation normally increases both ICaL and
  IKs together in the healthy heart, keeping QT relatively stable, but in LQT1 the
  mutant KCNQ1 channel fails to respond to PKA phosphorylation, so this balance is
  broken. LQT2 is vulnerable to auditory/emotional stimuli, and LQT3 to bradycardia,
  sleep, and other low-sympathetic-tone states.
* **Acquired/drug-induced LQTS (CiPA)**: drugs that block the hERG channel produce
  a phenotype similar to congenital LQT2, and hypokalaemia amplifies this effect
  through extracellular-K+-dependent reduction of IKr conductance. CiPA's qNet
  concept aims to predict actual TdP risk more accurately by integrating the
  degree of simultaneous ICaL and late-INa block rather than hERG block alone.

## 3. Model Structure

### 3.1 Mechanistic map — 14 clusters, 129 nodes

1. Genetics — ion-channel genotype (KCNQ1/KCNH2/SCN5A/rare genes)
2. Ion-channel dysfunction (molecular level: trafficking, gating, inactivation defect)
3. Ventricular action potential and ionic currents (Phase 0-4, INa/Ito/ICaL/IKs/IKr/IK1/INaK/INaCa)
4. Transmural heterogeneity and repolarisation dispersion (epicardium/M cell/endocardium/Purkinje, TDR, Tpeak-Tend)
5. EAD/arrhythmic substrate and triggering (ICaL reactivation, Ca2+ overload, CaMKII, re-entry)
6. Autonomic triggers (sympathetic activation, trigger specificity by exercise/auditory/sleep)
7. TdP -> VF -> SCD cascade
8. CiPA/hERG drug-block pharmacology (acquired LQTS, qNet index)
9. Beta-blocker PK/PD (propranolol, nadolol)
10. Mexiletine PK/PD (late-INa block, LQT3-targeted therapy)
11. Potassium/electrolyte regulation (oral K+, spironolactone, hypokalaemia risk)
12. Procedural/device interventions (LCSD, ICD)
13. Clinical diagnosis (QTc correction formulas, T-wave morphology, exercise stress testing, Schwartz score)
14. Clinical endpoints and risk stratification (QTc, TdP probability, SCD risk, registries)

### 3.2 mrgsolve model — 23 ODE compartments

* **Drug PK (14 compartments)** — propranolol (2-compartment: gut/central/peripheral),
  nadolol (2-compartment, renal clearance), mexiletine (1-compartment,
  CYP2D6[major] + CYP1A2[minor] metabolism), oral KCl (2-compartment),
  spironolactone (2-compartment), a hypothetical QT-prolonging drug X (2-compartment,
  a CiPA-style hERG blocker).
* **Disease PD (9 compartments)** — serum K+, genotype-specific channel-conductance
  indices (GKs/GKr/GNa-late), a sympathetic-drive index (SYMP_DRIVE), a lumped
  repolarisation-reserve/QTc surrogate (QTC), an EAD probability/substrate index
  (EAD_SUBSTRATE), cumulative TdP risk hazard (TDP_HAZARD), and cumulative expected
  TdP event count (TDP_EVENTS, the Poisson-intensity integral of a time-to-event
  survival analysis).
* QTc is computed with a simplified lumped turnover model: genotype-specific
  IKs/IKr deficits and late-INa excess are combined by weighted sum to compute a
  target QTc, which converges via first-order turnover (a mechanism-based reduced
  model, not the full 40-state O'Hara-Rudy action-potential model).

### 3.3 10 treatment/simulation scenarios

| # | Scenario | Calibration basis |
|---|---|---|
| 1 | Untreated LQT1 + exercise trigger | Priori 2004 NEJM, Ackerman 1999 Mayo Clin Proc |
| 2 | LQT1 + propranolol 2 mg/kg/day | Vincent 2009 Circulation |
| 3 | LQT1 + nadolol 1 mg/kg/day | Vincent 2009 Circulation (nadolol comparator arm) |
| 4 | Untreated LQT2 + auditory/emotional trigger, female | Schwartz 2001 Circulation, Priori 2003 Circulation |
| 5 | LQT2 + beta-blocker + K+/spironolactone | Schwartz 2001 Circulation (genotype-specific effects) |
| 6 | Untreated LQT3 + sleep/bradycardia trigger | Schwartz 1995 Circulation, Zareba 1998 NEJM |
| 7 | LQT3 + mexiletine | Moss 2000/2008 series, Ruan 2007 Circulation |
| 8 | LQT3 + mexiletine + propranolol combination | Schwartz 1995 Circulation (genotype-specific response) |
| 9 | Arbitrary genotype (LQT2) + addition of QT-prolonging drug X (acquired-on-congenital overlap) | CiPA (Colatsky 2016, Vicente 2019) |
| 10 | High-risk LQT2 after cardiac arrest: rescue with LCSD + beta-blocker + ICD | Schwartz 2004 Circulation, Zareba 2003 |

## 4. Shiny Dashboard (8 tabs)

1. **Patient/genotype profile** — set genotype, trigger, sex, syncope history.
2. **Drug PK** — plasma concentrations of propranolol/nadolol/mexiletine/QT-prolonging drug X.
3. **Ion channel/PD** — GKs/GKr/GNa-late conductance indices, sympathetic drive, qNet-like index.
4. **QTc & ECG surrogate** — QTc trajectory over time, risk threshold lines.
5. **TdP/arrhythmia risk** — EAD substrate index, cumulative risk hazard, cumulative TdP event probability.
6. **Clinical endpoints** — summary table of syncope/TdP/SCD surrogate endpoints.
7. **Scenario comparison** — simultaneous comparison of the 8 scenarios (QTc, TdP probability).
8. **Biomarkers/risk stratification** — Schwartz risk-score calculator, risk biomarker table.

## 5. How to Run

```bash
# 1) Render the mechanistic map
dot -Tsvg lqts_qsp_model_en.dot -o lqts_qsp_model_en.svg
dot -Tpng -Gdpi=150 lqts_qsp_model_en.dot -o lqts_qsp_model_en.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","shinydashboard","DT"))
library(mrgsolve)
mod <- mread("lqts_mrgsolve_model_en.R")
# See the LQTS_simulate_scenarios() helper in the comments at the bottom of the file
out <- mod %>% param(GENOTYPE=1, TRIGGER=1) %>%
  ev(amt=40, cmt="GUT_PROP", ii=8/24, addl=270) %>%
  mrgsim(end=2160, delta=6)
plot(out, "QTc_ms,TdP_event_probability,GKs_idx")

# 3) Run the Shiny dashboard
shiny::runApp("lqts_shiny_app.R")
```

## 6. Reference Summary

`lqts_references.md` contains 64 references organised into 10 sections (genetics,
ion-channel electrophysiology, epidemiology/registries, genotype-phenotype/
triggers, beta-blockers, mexiletine/LQT3, ICD/LCSD, CiPA/drug-induced LQTS,
diagnosis/risk stratification, reviews). Key references: Priori 2004 JAMA
(genotype-specific beta-blocker response), Vincent 2009 Circulation (propranolol
vs. nadolol), the Moss/Ruan mexiletine series (LQT3-targeted therapy), Schwartz
2001 Circulation (genotype-trigger specificity), and Colatsky 2016/Vicente 2019
(the CiPA qNet concept).

## 7. Limitations

* For research/education/hypothesis-generation only; must not be used for clinical
  decision-making.
* The QTc/repolarisation-reserve model is a reduced lumped turnover model, not a
  full multi-compartment cardiomyocyte action-potential model (e.g. O'Hara-Rudy,
  ten Tusscher). Real EAD occurrence is stochastic and depends on cell-to-cell
  coupling and tissue-level wave propagation, but this model simplifies it to a
  single deterministic index, EAD_SUBSTRATE. The transition to TdP is approximated
  as a Poisson-intensity integral (cumulative hazard), which reflects the
  qualitative direction of event rates observed in actual clinical trials/
  registries (relative risk reduction by genotype and by treatment) but whose
  absolute quantitative values have not been validated.
* Parameters such as mexiletine's late-INa block EC50/Emax and the hERG-IC50 of
  QT-prolonging drug X are approximations inferred from representative literature
  values and have not been fitted to individual patient data.
* This is a typical-value model that does not include inter-individual variability
  (IIV).
* This container has no R/mrgsolve runtime installed, so the mrgsolve code was
  completed through literature-based design and self-review of the code
  (compartments/dimensions/boundary-value checks), but its numbers were not
  verified by actual compilation and integration. The `.dot` file was actually
  rendered with Graphviz `dot`, and SVG/PNG generation was confirmed.
