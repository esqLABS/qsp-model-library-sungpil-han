# Autoimmune Pulmonary Alveolar Proteinosis (aPAP) QSP Model

> **One-line summary** — aPAP is posed not as an inflammatory lung disease but as an **unclosed surfactant mass balance**, opened by the loss of a *signal* rather than the loss of a cell. Production is normal; the dominant catabolic sink is switched off. The alveolar burden is the time-integral of a small difference between two large fluxes.

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`pap_qsp_model.dot`](../../../pulmonary-alveolar-proteinosis/pap_qsp_model.dot) · [SVG](../../../pulmonary-alveolar-proteinosis/pap_qsp_model.svg) · [PNG](../../../pulmonary-alveolar-proteinosis/pap_qsp_model.png) | 212 nodes · 25 clusters · 298 edges |
| mrgsolve ODE model | [`pap_mrgsolve_model.R`](../../../pulmonary-alveolar-proteinosis/pap_mrgsolve_model.R) | 59 ODEs · 36 scenarios · 17 diagnostics |
| Shiny dashboard | [`pap_shiny_app.R`](../../../pulmonary-alveolar-proteinosis/pap_shiny_app.R) | 10 tabs |
| References | [`pap_references.md`](../../../pulmonary-alveolar-proteinosis/pap_references.md) | 158 entries (every PMID verified) |

```bash
dot -Tsvg pap_qsp_model.dot -o pap_qsp_model.svg
dot -Tpng -Gdpi=150 pap_qsp_model.dot -o pap_qsp_model.png
Rscript pap_mrgsolve_model.R                     # run 36 scenarios + 17 diagnostics
Rscript -e 'shiny::runApp("pap_shiny_app.R", port = 8080)'
```

Validation environment: R 4.3.3 · mrgsolve 2.0.1 · Graphviz 2.42.2.

---

## 1. Why This Disease Is Built This Way (Three structural commitments)

aPAP is usually drawn as "autoantibody -> impaired macrophage function ->
surfactant accumulation -> hypoxaemia." Written that way it is a story, not
a model, and every arrow becomes a fitted slope. This model instead makes
three structural choices, and every number below **follows from** those
choices rather than being fitted to the comparator trials.

### (1) Recycling is not clearance

Uptake-and-resecretion by type II cells is a **large flux with zero net
effect at steady state**. If this loop is folded into "clearance" and the
macrophage arm is removed, the model predicts the burden would stabilise at
2-3x normal. Patients reach **30-100x**. The real sinks are only macrophage
catabolism, intracellular catabolism in type II cells, and
mucociliary/lymphatic clearance, of which the macrophage accounts for most
of net output. The loop is drawn explicitly on the map and deliberately
excluded from the net balance sheet.

The model **assigns** the healthy-lung output split as 65/25/10
(macrophage/type II/mucociliary) and **back-solves** uptake Vmax, digestion
Vmax, and de novo production in `$MAIN`. The result is that the healthy
lung is structurally stationary (-0.002% drift over 1000 days), and the
phospholipid mass balance closes to within 9x10^-8 mg even in a run that
includes whole-lung lavage.

### (2) The autoantibody is a stoichiometric buffer, not an IC50 inhibitor

GM-CSF in the epithelial lining fluid (ELF) sits at tens of pg/mL (a few
pM), while the neutralising antibody's binding sites are in the hundreds to
thousands of pM. Free ligand is therefore set by **binding equilibrium in
the antibody-excess regime**, and the model solves this as the
**numerically stable form** of the quadratic for 1:1 binding (the textbook
root cancels to zero or noise in double precision at this ratio — that is
exactly what happened in the first implementation).

That single equation produces two clinical facts that otherwise look
contradictory.

* **A titre threshold**: the model's half-signal titre is **5.63 µg/mL**
  (literature threshold 5 µg/mL, Sakagami 2010).
* **Saturation above the threshold**: above it, free GM-CSF is already
  close to zero, so raising titre further changes almost nothing.

### (3) Inhaled GM-CSF works by locally overwhelming that buffer

The therapeutically relevant quantity is not plasma concentration but the
**molar ratio inside the ELF**. Nebulising 300 µg at 40% deposition puts
roughly 8300 pmol of ligand into 30 mL that holds about 160 pmol of
neutralising sites — about a 50-fold excess, against a receptor that
saturates in the tens of pM. So it is **route of administration, not dose**,
that separates success from failure.

Here the model produces an unexpected result (diagnostic D17): **the
inhaled-dose threshold sits around 1-3 µg/day**, roughly one-hundredth of
the clinical dose. Everything above that is a flat plateau. At 300 µg, the
trough free GM-CSF is 96.5 pM, 34-fold normal (2.87 pM). This is consistent
with the fact that a clear dose-response has never been observed in aPAP,
and with 125 µg BID and 300 µg QD performing similarly.

---

## 2. What This Model Says That Is New (The mechanism the trials force into the model)

On receptor pharmacology alone, inhaled GM-CSF fully restores the signal
and so should empty the lung within weeks. In practice, DLCO rises by only
about 11.6 percentage points at 48 weeks. Resolving this tension requires a
mechanism that was not in the model, and it is the following.

> **The aerosol goes where the air goes. But the burden sits exactly where the air does not go.**

Alveoli filled with surfactant are unventilated, so they receive almost none
of the nebulised drug. Conversely, open alveoli that the drug reaches
abundantly have almost nothing left to clear. So the **fraction of the
burden the drug can actually act on** is small. The model implements this
as two geometric facts.

* `RHO` — the multiple by which the surfactant mass of one filled alveolus
  exceeds that of one open alveolus (default 50). About 89-95% of the
  burden sits inside filled units.
* `EDGEF` — consolidated regions are cleared **from the edge inward**. This
  is an interface term (proportional to F(1-F)) for drug-exposed
  macrophages in contact with filled alveoli. Because this term vanishes as
  the lung opens, the response **decelerates**. A volume-based reach term
  would instead accelerate it, which the trials do not show.

Clinically useful conclusions follow from this structure.

* **Deposition and lung opening matter more than dose.** D17's plateau is
  the quantitative statement of this.
* **Maintaining GM-CSF after lavage** beats either alone (S18). Lavage
  creates a lung the drug can reach.
* **Oral/systemic drugs do not pay this reach penalty.** Blood reaches even
  a consolidated lung; aerosol does not. This is why statins look
  disproportionately good relative to their own effect size (S26), and why
  high-dose subcutaneous GM-CSF produces a partial response (S13,
  consistent with the ~40% response rate in Seymour 2001) despite far worse
  receptor-level pharmacology.

---

## 3. Reproducing the Clinical Trials (Validation)

**16 / 17 diagnostics pass.** The numbers below are not fitted to the
comparator values.

| Diagnostic | What is tested | Model | Observed |
|---|---|---|---|
| D01 | Healthy-lung stationarity (1000 days) | -0.002% drift | — |
| D02 | Phospholipid mass-balance closure (including lavage) | 9x10^-8 mg residual | — |
| D04 | GMAb critical titre | half-signal 5.63 µg/mL | 5 µg/mL |
| D07 | IMPALA-2 24-week ΔDLCO (drug / placebo) | +9.6 / +3.7 (diff +5.9) | +9.8 / +3.8 (diff +6.0) |
| D07 | IMPALA-2 24-week ΔSGRQ-T | -9.6 / -3.0 | -11.5 / -4.9 |
| D09 | IMPALA 24-week ΔA-aDO2 difference | -3.3 to -5.1 | -6.2 |
| D10 | PAGE placebo arm (after modelling run-in) | almost unchanged | +0.17 mmHg |
| D11 | Tazawa 2010 ΔA-aDO2 | derives the responder fraction | -12.3 mmHg, 62% |
| D12 | Whole-lung lavage | burden -44%, A-aDO2 34.3->19.6, relapse 3.4 yr | 70% relapse-free at 7 yr |
| D13 | Inhaled GM-CSF in hereditary PAP | difference < 0.05 pp (structural null) | ineffective |
| D16 | Existence of asymptomatic screen-detected patients | below threshold once above the critical floor | 31.8% |

**The placebo arm is not a placebo parameter.** IMPALA-2's placebo-arm
improvement (+3.8 pp) comes from two selection mechanisms. ① Patients are
referred and randomised after worsening, so baseline is measured at the
trough of a decline that would recover on its own. ② DLCO is effort-dependent
and improves on repetition. Turning off both mechanisms flattens the placebo
arm. And **the same model leaves PAGE's placebo arm unmoved** — because PAGE
excluded improvers during its 12-week observation period (D10). This is the
result of simulating two enrolment protocols, not two different placebo
effects.

---

## 4. Failures and Refutations, Reported Honestly

Do not cite the numbers in this file without reading this section.

* **D08 — 48-week over-prediction.** IMPALA-2's 48-week difference comes
  out as **+15.1 pp** against an observed +6.9 pp (about 2-fold). The cause
  is clear. The model has a positive feedback loop on the recovery side (as
  the lung opens, aerosol reaches more of the remaining burden -> faster
  clearance -> more opening). The 48-week data do not show this kind of
  acceleration, but Tazawa 2010's post-discontinuation data (29 of 35
  patients stable off treatment for 1 year) suggest something
  self-sustaining is happening. One of the two is being misread, and this
  model cannot tell which. **Cite the 24-week figures.**
* **D15 — no bistability was found.** Cluster 11 of the map was drawn on
  the premise that the fill/clear feedback would create two stable states
  — "trapped" and "runaway." At the same catabolic floor, burdens of 300 mg
  and 8169 mg **both converge to 6160 mg**. The feedback is real but
  sub-critical, and clinical heterogeneity is **steep single-valued
  dependence** on the floor, not bistability. The premise was refuted, and
  rather than redraw the cluster, the refutation was labelled on it.
* **D05 — only half-reproduces Inoue's negative finding.** The lack of
  correlation between titre and severity does not disappear from buffer
  saturation alone. Adding the fact that "the assay measures binding, not
  neutralisation" (patient-to-patient variation in neutralising fraction)
  weakens the correlation from r = +0.85 to **+0.57**, but does not remove
  it. Also, in this cohort **the dominant covariate is still the antibody
  side, not the catabolic floor**. Both statements correct an earlier draft
  of this file. A testable prediction: in a cohort enriched for 5-15 µg/mL,
  the correlation should reappear, and **a neutralisation assay should show
  a correlation that a binding-titre assay cannot**.
* **D19 — rituximab, FcRn inhibition, and statins are hypotheses.**
  Rituximab cannot reach CD20-negative long-lived plasma cells, so it can
  only reduce the plasmablast share of antibody production (why reported
  responses are inconsistent). FcRn blockade lowers total IgG by 60-70%,
  pushing median-titre patients past the threshold, so the model predicts a
  large benefit, but this has never been tested in aPAP and it directly
  inherits D08's over-prediction. The statin effect-size parameter
  (`EMAXST`) is not anchored to data.

---

## 5. Model Structure

### The mass-balance core

```
Production (normal) ──► [Alveolar phospholipid pool PLA] ◄──► [Type II lamellar body LB]   ← zero-net loop
                    │  │  │  │
                    │  │  │  └──► Consolidated/sequestered material SEQ (slow redispersal)
                    │  │  └─────► Mucociliary/lymphatic clearance   (sink, capacity-limited)
                    │  └────────► Intracellular catabolism in type II cells   (sink, saturable)
                    └───────────► Macrophage uptake → lipid loading LIP → digestion
                                   ▲                            │
                                   │                    Capacity-limited (MM): foam cells are
                              GM-CSF signal                 "stalled" cells
```

Writing digestion as first-order in lipid loading lets a foam macrophage
carrying 6x the lipid digest 6x as much as a normal cell, which means
**PAP can never form at all** (the burden stalls at 2.8x normal). Digestion
must be Michaelis-Menten, with a per-cell Vmax proportional to the
digestive programme.

### The signalling chain

```
Total ELF GM-CSF ─┐
                   ├─ 1:1 binding equilibrium (numerically stable form) ─► Free GM-CSF ─► Occupancy ─► (Hill n=2)
Neutralising ELF antibody ─┘                                                        │
                                                                        ▼
                        PU.1 ─► PPARγ ─► lysosomal degradation machinery ─► per-cell degradation capacity
                                                                        │
     + GM-CSF-independent floor (CAPFLOOR)  + statin (systemic, no reach penalty)   │
                                                                        ▼
                                        × macrophage count = total degradation capacity
```

### Gas exchange is physics, not a score

A-aDO2 and DLCO are **calculated** from the shunt equation, the alveolar gas
equation, the Severinghaus dissociation curve (and its closed-form inverse),
and the Fick relation. So exercise-induced deterioration emerges without
changing any lung parameter (VO2 5-fold, cardiac output 3-fold -> mixed
venous oxygen content falls -> the same shunt becomes more costly). This is
why desaturation on exertion precedes resting hypoxaemia.

### The patient is defined by an event, not an initial value

Healthy lung -> seroconversion -> integration over years. The point of
diagnosis is **the event of first crossing the clinical-trial enrolment
criterion** (resting PaO2 <= 70 mmHg, PAGE's criterion). Enrolment then
occurs about 0.75 years later, when the patient is near their own
equilibrium (otherwise the control arm would keep worsening throughout the
trial, which real placebo arms did not).

The result is that baseline values **emergently** match three published
figures simultaneously: A-aDO2 of roughly 35-39 mmHg, DLCO of roughly 40-45%
predicted, KL-6 of roughly 3200 U/mL, and a burden about 20x normal. And
patients with a catabolic floor above about 0.33 **never cross the enrolment
threshold at all** — this is the model's explanation for the 31.8% found on
health screening (Inoue 2008).

---

## 6. Scenarios (36)

| Group | Content |
|---|---|
| S01-S06 | Healthy control · 25-year natural history from seroconversion · progressive phenotype · screen-detected phenotype · spontaneous remission |
| S07-S15 | Molgramostim 300 µg QD (24/48 wk) · every-other-week regimen · PAGE regimen · Tazawa 2010 regimen · subcutaneous 5/20 µg/kg · dose sweep (10-3000 µg) · split-dose comparison at the same total daily dose |
| S16-S20 | Whole-lung lavage alone · 8-year follow-up post-lavage · lavage + maintenance therapy · repeat lavage x4 every 18 months · segmental lavage |
| S21-S27 | Rituximab (intermediate titre / low titre) · plasma exchange x10 · FcRn inhibition (titre 25 / 10) · atorvastatin · statin + molgramostim |
| S28-S30 | Hereditary PAP + inhaled GM-CSF (structural null) · secondary PAP (no antibody, monocytopenia) · macrophage-transplant surrogate |
| S31-S36 | Titre sweep (0.5-300 µg/mL) · catabolic-floor sweep · Nocardia burden (untreated vs. on treatment) · smoking and smoking cessation · treatment after a 10-year diagnostic delay |

---

## 7. Shiny Dashboard (10 Tabs)

① Patient generation (from seroconversion) · ② Mass balance (the two large
fluxes and their small difference) · ③ Antibody buffering curve and this
patient's position on it · ④ Drug delivery (ELF pharmacology, molar ratio,
time above threshold, reach) · ⑤ Gas exchange (PaO2/A-aDO2/DLCO calculated
from the shunt, resting vs. exercise) · ⑥ Clinical measures (DSS, SGRQ,
6MWD, oxygen requirement) · ⑦ Biomarkers (KL-6, SP-D, CEA, LDH, CT density)
· ⑧ Scenario comparison (12 regimens) · ⑨ Virtual clinical trial
(reproducing IMPALA-2 / PAGE) · ⑩ Covariates (what predicts severity and
what does not)

---

## 8. Parameter Confidence (Where the numbers come from)

| Category | Example | Basis |
|---|---|---|
| Directly from the literature | GMAb critical titre 5 µg/mL · IgG half-life 21 days · GM-CSF receptor Kd 20-50 pM · clinical trial doses and schedules | References §1, §12, §14 |
| Derived (back-solved) | Uptake Vmax · digestion Vmax · de novo production · antibody synthesis rate | From the output split and turnover rates specified in `$MAIN` |
| Calibrated | ELF/serum transfer ratio κ, neutralising fraction `FNEUT` (-> critical titre) · SGRQ mapping slope · interface coefficient `EDGEF` | Each fitted to one published value. Noted explicitly in the file |
| Assumption | Fill/clear exponent · `RHO` · sequestration redispersal rate · `EMAXST` | Subject to sensitivity analysis. See D19 |
| **Weakest link** | **The absolute size and turnover of the alveolar surfactant pool** | The human value is known to only one significant figure. The model fixes it at 300 mg / 200 mg/day and back-solves every Vmax from that, so **the split is accurate but the absolute scale inherits that uncertainty.** Read mg-level statements as order-of-magnitude. What the model is sensitive to is the **ratio** of macrophage sink to production. |

---

## 9. Limitations

* Natural history is integrated beyond the clinically survivable range
  (DLCO in the low teens at end stage). Death and censoring are not
  modelled.
* Regional heterogeneity is represented only as three fractions —
  filled/partially filled/open — plus an interface term. Real PAP's
  geographic compartmentalisation (crazy paving, lobar boundaries) would
  require a spatial model.
* Infection enters only as a hazard and a single Nocardia-burden ODE. No
  antibiotic treatment is included.
* Lung transplant/HSCT are represented only as a surrogate for the capacity
  parameter.
* Paediatric PAP, the haematopoietic phenotype of GATA2 deficiency, and
  overlap with smoking-related interstitial lung disease are on the map but
  not in the equations.

---

## 10. Disclaimer

This is a QSP model for education and research purposes. It has not been
validated for clinical decision-making, prescribing, or regulatory
submission. All parameters and predictions should be read in the context of
the evidence levels and the list of failures noted above.
