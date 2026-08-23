# Common variable immunodeficiency (CVID) — QSP model

**Primary antibody deficiency**

> CVID is **two diseases** sharing one laboratory definition,
> and the treatment we give treats **only one** of them.

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`cvid_qsp_model.dot`](cvid_qsp_model.dot) · [SVG](cvid_qsp_model.svg) · [PNG](cvid_qsp_model.png) — 234 nodes · 20 clusters |
| ⚙️ mrgsolve model | [`cvid_mrgsolve_model.R`](cvid_mrgsolve_model.R) — 74 ODEs · 224 parameters · 30 scenarios · 34 diagnostics (34 PASS) |
| 📊 Shiny dashboard | [`cvid_shiny_app.R`](cvid_shiny_app.R) — 10 tabs |
| 📚 References | [`cvid_references.md`](../../../common-variable-immunodeficiency/cvid_references.md) — 124 PubMed citations (all confirmed by lookup through the E-utilities) |

---

## 1. The premise

Every CVID patient satisfies the same criteria. IgG below −2SD for age, reduced IgA
and/or IgM, impaired vaccine responses, no other explanation. Out of this single
definition come two almost unrelated diseases.

| | ARM 1 · antibody deficiency | ARM 2 · immune dysregulation |
|---|---|---|
| **Presentation** | Recurrent upper and lower respiratory infection (encapsulated bacteria) | Autoimmune cytopenias · GLILD · enteropathy · granulomas · polyclonal lymphoproliferation · lymphoma |
| **Its essence** | Not having antibody | A disorder of T cells and tolerance |
| **Frequency** | Effectively everyone | About 20–30% |
| **Predictability** | Predicted **quantitatively** from the serum IgG concentration | Genotype- and phenotype-dependent |
| **Immunoglobulin replacement** | **Largely solves it** | **Essentially no effect** |
| **Contribution to mortality** | Low after replacement | **This is what kills people** (Resnick 2012: RR ≈ 11-fold) |

A CVID model that models antibody alone has modelled **only the part that is already
solved.** So this model runs the two arms **in parallel** from a common genetic and
cellular root, has each drug act on exactly one arm, and reproduces **from the
structure** the fact that maximal-dose IgG replacement barely moves GLILD — a
structure, not an assertion.

### And a third layer that neither arm can reverse

```
ARM 1 (antibody deficiency, replaceable) ─┐
                                          ├─→  irreversible structural damage  ─→  clinical endpoints
ARM 2 (immune dysregulation, not)        ─┘     (bronchiectasis · interstitial fibrosis · NRH)
```

Bronchiectasis is the archetype, and **an irreversible state variable sits inside a
positive feedback loop** — it is a ratchet.

```
colonisation → IL-8 → neutrophil elastase → airway wall destruction
      → bronchiectasis (irreversible) → reduced mucociliary clearance
      → mucus stasis → more colonisation ──┘
```

So the most important number in CVID is not the dose but the **diagnostic delay**
(historically 4–7 years). A large part of this model's reason for existing is to put a
price on that delay.

---

## 2. Five design decisions that define the model

### (1) The exposure-response is a Hill function, not a straight line

The meta-regression of 17 IVIG studies by Orange 2010 is the quantitative backbone of
this field: **each 100 mg/dL rise in IgG trough reduces the pneumonia incidence by
about 27%**, with 0.113 events/patient-year at a trough of 500 mg/dL. That slope is
robust.

The corollary usually quoted alongside it — **zero pneumonia at a trough of 1400
mg/dL** — is an artefact of *fitting a straight line to a rate*. And put it into a
model and the model will say "push to 1400 and forget about it".

This model instead uses

```
rate = RMAX / (1 + (OPSONIN / C50)^HILL)     C50 = 200 mg/dL, HILL = 2.4, RMAX = 0.7745/year
```

and the parameters are values **solved out of the anchors, not asserted**:

| Test | Model | Target |
|---|---|---|
| d(ln rate)/dC at C = 700 | **−29.0% / 100 mg/dL** | −27% (Orange) |
| Incidence at C = 500 | **0.113/year** | 0.113/year (Orange) |
| Downward extrapolation to C = 250 / 100 | **0.40 / 0.72 /year** | The untreated and severe range |
| At a trough of 900 | **0.028/year** | 0.02–0.05 observed on replacement |
| Healthy control | **0.019/year** | Adult community-acquired pneumonia incidence |

The linear form cannot be extrapolated either below or above it. The Hill form improves
without reaching sterility and creates **a knee of diminishing returns around 800–1000
mg/dL**.

### (2) Risk is driven by the instantaneous concentration, not the trough — and that gives a result for free

It is the mechanistically honest choice. And that choice gives one result for free.
Because `rate(C)` is **convex** over the therapeutic range, **Jensen's inequality**
means that a fluctuating profile with *the same mean concentration* has a **higher mean
event rate** than a flat one.

That is, weekly SCIG being superior to 4-weekly IVIG at identical AUC is not
immunology but **the shape of the curve**. The model measures the size of this penalty
(diagnostic D6):

| Regimen | Mean IgG | Trough | Fluctuation | Pneumonia/year | Excess |
|---|---|---|---|---|---|
| IVIG 500 q4wk | 959 | 731 | 897 | 0.0275 | **+9.6%** |
| IVIG 375 q3wk | 959 | 781 | 668 | 0.0259 | +3.2% |
| SCIG 171 weekly (×1.37) | 944 | 895 | 90 | 0.0254 | +1.2% |
| SCIG 86 twice weekly | 947 | 926 | 35 | 0.0251 | 0.0% |
| fSCIG 538 q4wk | 939 | 728 | 501 | 0.0279 | +11.2% |

IVIG has a pneumonia rate **9.6% higher despite** a mean **1.8% higher**. A pure
convexity penalty. The model was not built to produce this result; it fell out.

### (3) The dose-response of pneumonia and of sinusitis are different

Serum IgG penetrates poorly to the airway surface, and secretory IgA is **not replaced
by any product**. So the mucosal infection rate has **a floor independent of IgG**,
which is why patients on replacement still get 2–3 episodes of sinusitis and
bronchitis a year.

| IVIG dose | Trough | Pneumonia (relative) | Sinusitis (relative) |
|---|---|---|---|
| 300 mg/kg | 590 | 1.000 | 1.000 |
| 400 | 713 | 0.591 | 0.891 |
| 500 | 833 | 0.391 | 0.833 |
| 600 | 949 | 0.278 | 0.799 |
| 800 | 1176 | 0.161 | 0.754 |
| 1000 | 1395 | **0.104** | **0.726** |

**Clinical implication:** dose escalation for recurrent **pneumonia** is worthwhile
(−90%). Dose escalation for persistent **sinusitis** is largely pointless (−27%, even
at a 3.3-fold escalation). What acts on that floor is not dose but azithromycin
(model: 1.52/year vs a limit of 2.31/year at any dose).

### (4) The GLILD equations contain no IgG term — and the diagnostic refuted the author's expectation

IgG does not appear **directly** in `dxdt_LYMPHAGG` or `dxdt_GRAN`. Diagnostic D7 was
originally written to require **exactly zero** across an IVIG sweep from 300 to 1000
mg/kg.

**The model refuted that expectation.** A residual sensitivity of 8.3% remains. Traced
back, it is the **indirect** path
`IgG → colonisation → chronic activation → CD21low B → lymphoid aggregates`, and that
path is biologically real (chronic infection really does drive lymphoproliferation).
So it is **reported rather than removed**. The claim the model can defend is the weaker
but more useful one:

> Across the whole clinical dose range, replacement therapy moves GLILD by **8%**.
> In the same model, rituximab + azathioprine moves it by **80%**. **About a 10-fold
> difference.**

### (5) Irreversible states have no removal term

`BE` (bronchiectasis), `FIB` (interstitial fibrosis), `FEV1_IRR` and `NRH` are
monotonically non-decreasing **by construction** — by the form of the equations, not by
a choice of parameters. Diagnostic D8 verifies monotonicity numerically across all 30
scenarios (30/30).

---

## 3. The price of diagnostic delay (the headline result)

The same IVIG 500 mg/kg q4wk. A different start time. At 20 years:

| Diagnostic delay | IgG trough | Bronchiectasis (Reiff 0–18) | FEV1 (%pred) | Cumulative pneumonias | QoL |
|---|---|---|---|---|---|
| 1 year | 731 | 1.4 | **93.2** | 0.9 | 0.843 |
| 4 years | 769 | 6.7 | **73.8** | 2.4 | 0.788 |
| 7 years | 810 | 11.5 | **53.1** | 4.5 | 0.733 |
| 15 years | 729 | 15.7 | **44.4** | 12.5 | 0.703 |

**The trough is the same in all four cases.** The whole of the difference is
irreversible. The FEV1 gap between a 1-year and a 15-year delay is **48.8 %pred**, and
it is not recovered at any dose.

> The most valuable intervention in CVID is therefore not a better immunoglobulin but
> a **single earlier serum IgG measurement**.

**Note:** the optimal combined regimen (early diagnosis + SCIG to a trough of ~1180 +
azithromycin prophylaxis) reduces pneumonia from 0.0224 to 0.0123/year and sinusitis
from 2.68 to 1.52/year against the standard regimen. Yet the whole of that effect
combined is **smaller than the loss from one year of delay** — in the same model, one
year of delay already leaves 0.96 Reiff points.

---

## 4. Four trade-offs

### ① Route and interval vs trough
Table (2) above. The same milligrams, a different trough. It was built to **fall out**
of the PK (depot · lymphatic absorption · FcRn recycling · saturable catabolism are all
explicit nodes).

### ② Immunosuppression vs infection
Every ARM 2 treatment raises net immunosuppression, and that comes straight back as
ARM 1 susceptibility. But there is an **asymmetry** specific to CVID:

| | Normal host | CVID (on replacement) |
|---|---|---|
| Total B cells after rituximab | **−98.1%** | **−98.1%** |
| Serum IgG after rituximab | −11.1% | **−2.7%** |

The main humoral price of rituximab in other diseases is secondary
hypogammaglobulinaemia. **In CVID that price has already been paid and is already being
replaced.** Rituximab is therefore relatively cheap in CVID, and the real price is not
humoral but **cellular**. (Incidentally: because the B-cell sink disappears, soluble
BAFF rises **2.38-fold** — observed 2–5-fold.)

### ③ Splenectomy vs sepsis
Four strategies in refractory ITP, at 20 years:

| Strategy | Platelets | Cumulative invasive infections | Net immunosuppression |
|---|---|---|---|
| High-dose IVIG + prednisone | 76 | — | 0.200 |
| Rituximab | 54 | — | 0.001 |
| **Splenectomy** | 82 | **0.207** | 0.000 |
| **Eltrombopag** | 72 | **0.083** | 0.000 |

Splenectomy raises platelets well (it removes the main **site of clearance** of
antibody-coated platelets), but the autoantibody itself is unchanged. And in asplenia
in a patient with no opsonising antibody of their own, the encapsulated-organism risk
is **multiplicative, not additive** — invasive infections in the model are 2.5 times
those on eltrombopag.

### ④ A fixed trough target vs an individual "biological IgG level"
Bonagura's point is that the dose which stops infections differs from patient to
patient and can be far higher than any population target. The model contains both a
fixed-target strategy and a treat-to-effect strategy.

**Protein-losing enteropathy is the extreme case.** The model **solves** the dose:

| | Trough | IgG clearance |
|---|---|---|
| Standard patient, 500 mg/kg q4wk | 811 | 1.57 dL/day |
| PLE patient, the same 500 mg/kg | **378** | 2.98 dL/day (1.89-fold) |
| PLE patient, **the solved 1426 mg/kg** | 823 | — |

Note: the required **dose ratio of 2.85-fold** exceeds the **clearance ratio of
1.89-fold**. Because FcRn-mediated catabolism is itself saturable, pushing the
concentration back up raises clearance **again**. It is not "just give twice as much".

---

## 5. Model structure (74 ODEs)

| Compartment group | Number | Contents |
|---|---|---|
| Immunoglobulin replacement PK | 3 | SC depot · central · peripheral (interstitial) |
| B-cell development · class-switch block | 9 | transitional → naive → GC → switched memory / IgM memory / CD21low / plasma cells |
| BAFF / APRIL / sBCMA | 3 | receptor-weighted sink |
| T-cell help · regulation · cytokines | 10 | naive and memory CD4 · cTfh · Treg · activated CD8 · TEMRA · IFN-γ · IL-21 · IL-6 · CXCL13 |
| Innate immunity · mucosal barrier · microbial translocation | 6 | barrier · LPS · sCD14 · monocytes · ISG · airway neutrophils |
| Antibody breadth (donor pool) | 1 | |
| Pathogens · cumulative infection | 6 | colonisation · Pseudomonas · GI pathogens · pneumonia/sinusitis/invasive counters |
| Airway inflammation and the ratchet | 6 | inflammation · mucus · **BE (irreversible)** · **FEV1_IRR (irreversible)** · FEV1_REV · exacerbations |
| GLILD · fibrosis | 4 | lymphoid aggregates · granulomas · **FIB (irreversible)** · respiratory failure |
| Lymphoproliferation · spleen · liver | 3 | spleen · lymph nodes · **NRH (irreversible)** |
| Autoimmune cytopenias | 3 | autoantibody · platelets · haemoglobin |
| Enteropathy · albumin | 2 | (→ feeds back into IgG clearance) |
| Malignancy risk | 1 | |
| Mortality hazard · QoL · cumulative immunosuppression | 4 | infectious / non-infectious hazards separated |
| Drug PK · effect compartments | 13 | rituximab · abatacept · sirolimus · leniolisib · prednisolone · azathioprine · eltrombopag · azithromycin · JAK |

### The initial conditions are solved, not written by hand

- **Healthy baseline**: a 40-year pre-run with `HEALTHY=1` → reference physiology. The
  subsequent 20-year drift is diagnostic D1 (IgG drift **0.0076%**).
- **CVID baseline at the time of symptom onset**: because the class-switch block
  precedes symptoms, a 25-year pre-run with the block switched on and damage
  accumulation switched off (`DAMAGEON=0`) → the CVID steady state of the immune
  compartments. Then only the damage states are initialised to zero (t=0 = **the first
  symptom**, the ratchet not yet turning). The long-lived IgG plasma cell pool is
  seeded with its **legacy value** from before the block, so untreated serum IgG starts
  at 250 mg/dL and **falls** over decades — this is a prediction of the model, and it
  agrees with the observation that untreated hypogammaglobulinaemia deepens.

---

## 6. Validation (34 diagnostics, 34 PASS)

Every one of them is a test that **can fail**. They were run under real mrgsolve 2.0.1.

| ID | Test | Result |
|---|---|---|
| D1 | 20-year drift of healthy baseline IgG | 0.0076% |
| D2 | 5-year stability of the CVID onset state | 0.357% |
| D3 | Terminal IgG half-life after a single IV dose | **36.5 days** (target 30–40) |
| D4 | SC/IV AUC ratio → EU correction factor | F=0.746, **1.341** (target 1.37) |
| D5 | Exposure-response slope at C=700 | **−29.0%/100 mg/dL** (Orange −27%) |
| D5b–d | Anchor point · downward extrapolation · healthy control floor | 0.113 / 0.40 · 0.72 / 0.019 /year |
| D6 | Agreement of mean IgG under AUC matching | 959 vs 946 |
| D6b | **Excess risk arising from the profile shape alone** | **+7.7%** |
| D6c | Fluctuation of IVIG 400/500 mg/kg | 654 / 898 mg/dL |
| D7 | IgG dose sensitivity of GLILD (structural zero) | 8.3% (no direct term) |
| D7b/c | Rituximab+AZA vs IgG dose | −79.7% vs 8.3% → **10.1-fold** |
| D8 | BE · FIB monotonically non-decreasing (all scenarios) | **30/30** |
| D9/9b | Monotone loss of FEV1 · BE with delay | 48.8 %pred, ceiling not reached |
| D10/10b | Rituximab asymmetry (identical depletion, different IgG price) | −98.1%/−98.1%, −11.1%/−2.7% |
| D11/11b | BAFF: after rituximab / untreated CVID | 2.38-fold / 1.70-fold |
| D12–12d | PLE trough collapse · solved dose · FcRn overshoot | 378, 1426 mg/kg, 2.85 vs 1.89 |
| D13/13b | Platelets on the 4 ITP strategies / infection price of splenectomy | 76 · 54 · 82 · 72 / 0.207 vs 0.083 |
| D14/14b | Mortality hazard ratio for non-infectious complications / non-infectious share | **10.5-fold** (Resnick ~11) / 100% |
| D15/15b | Leniolisib index lesion / naive B % | **−21%** (observed −39%) / rising |
| D16–16c | Divergence of the pneumonia vs sinusitis dose-response / azithromycin | −82% vs −19% / 1.52 vs 2.31 |

### What the model refuted about itself is reported as it stands

1. **D7 (structural zero):** exactly zero was expected, and an indirect sensitivity of
   8.3% remained. The path is real, so it is reported rather than removed.
2. **D15 (leniolisib):** removing an arbitrary lymph-node shrinkage term and leaving
   only PI3Kδ inhibition **underpredicts** at −21% against the observed −39%. It was
   not made to match.
3. **An error found during writing:** of five PMIDs written from memory in the draft,
   four pointed at entirely different papers (lymphoedema · mantle cell lymphoma ·
   asthma prescribing · multiple myeloma). All were replaced with the results of NCBI
   E-utilities lookups, and all 124 references were confirmed the same way.

---

## 7. How to run it

```bash
# rendering the mechanistic map
dot -Tsvg cvid_qsp_model.dot -o cvid_qsp_model.svg
dot -Tpng -Gdpi=150 cvid_qsp_model.dot -o cvid_qsp_model.png

# model + 30 scenarios + 34 diagnostics (a real mrgsolve run, takes several minutes)
Rscript cvid_mrgsolve_model.R

# interactive dashboard (10 tabs)
Rscript -e 'shiny::runApp("cvid_shiny_app.R")'
```

Packages required: `mrgsolve` (≥2.0), `dplyr`, `tidyr`, `ggplot2`, `shiny`, `DT`, and
Graphviz (`dot`).

---

## 8. Shiny dashboard (10 tabs)

| Tab | Contents |
|---|---|
| 1 · Patient profile | Immune phenotype, 20-year outcomes, EUROclass/Freiburg coordinates |
| 2 · IgG PK | Whole time course + a zoom on the steady-state profile, trough/fluctuation compared by route |
| 3 · Exposure-response | Hill form against linear form, local slope, **visualisation of the Jensen penalty** |
| 4 · Infection | Pneumonia · sinusitis · invasive infection, **the dose-response divergence** |
| 5 · Irreversible damage | Each stage of the ratchet, delay comparison, lung function |
| 6 · ARM 2 | GLILD · fibrosis · spleen · lymph nodes, cytopenias, enteropathy |
| 7 · The price of immunosuppression | Net immunosuppression → susceptibility → infection, **the rituximab asymmetry** |
| 8 · Scenario comparison | 5 comparison sets × 10 endpoints, summary table |
| 9 · Biomarkers | Opsonin vs mucosal defence (= the IgA gap), B-cell composition, BAFF and the sink |
| 10 · Validation & diagnostics | Table against the literature anchors, results of the 34 diagnostics |

---

## 9. Limitations

- **`DYSGENO` is a single aggregate variable.** As Chapel 2008 showed, the real
  phenotypes (autoimmunity / lymphoproliferation / enteropathy / lymphoma) are
  considerably exclusive of one another. The model approximates this by raising only
  one arm in each scenario family, but it does not generate phenotype exclusivity
  itself as a mechanism.
- **There is no population variability.** It is a single representative patient. The
  distribution of the "individual biological IgG level" and a population-level
  comparison of fixed-target against treat-to-effect strategies can only be handled
  properly by adding `$OMEGA`.
- **The mortality hazard constants are only weakly pinned.** The non-infectious hazard
  ratio of 10.5-fold matched Resnick's roughly 11-fold, but the absolute survival
  curves need separate calibration against cohort data.
- **The untreated IgG decline trajectory in CVID** is a prediction of the model and has
  not been validated quantitatively (the direction agrees with observation).
- **The relation between bronchiectasis score and FEV1** (1.7 %pred per Reiff point) is
  an approximation taken from a cross-sectional correlation.
- Infections are modelled as a **continuous hazard rate**. Discrete events and their
  attendant acute trajectories (admission, courses of IV antibiotics) are averaged out.

---

## 10. Disclaimer

This model is a **quantitative systems pharmacology model for educational and research
purposes**. It was calibrated to aggregate data from the public literature but has not
been validated or certified against individual patient data, and **must not be used
directly for real clinical decision-making, prescribing, or regulatory submission.**
