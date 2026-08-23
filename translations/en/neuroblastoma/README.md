# High-Risk Neuroblastoma — QSP Model

> **Summary in one sentence.** The **efficacy and the dose-limiting toxicity** of
> anti-GD2 antibody are carried by **two different Fc effector arms** (ADCC is
> first-order in bound IgG per cell, CDC second-order), and those two arms sit in
> **different compartments whose permeabilities differ 200-fold**. From these two
> facts alone, **as arithmetic rather than as assertion**, follow: why this drug
> eradicates marrow minimal residual disease but cannot eradicate bulk solid
> tumour, why raising the dose buys no benefit, and why the K322A Fc mutation is
> the only thing that moves the therapeutic index.

> **In one sentence.** Anti-GD2 antibody's efficacy and its dose-limiting toxicity
> are carried by **two different Fc effector arms** — ADCC linear in bound IgG per
> cell, CDC quadratic — sitting in **compartments whose permeabilities span
> 200-fold**. From those two facts alone, and not from any assertion, follow: why
> the drug clears marrow minimal residual disease but not bulk soft tissue, why
> dose escalation buys pain and not efficacy, and why the K322A Fc mutation is the
> only lever that moves the therapeutic index.

---

## Contents

| File | Contents |
|------|------|
| [`nb_qsp_model.dot`](../../../neuroblastoma/nb_qsp_model.dot) | Mechanistic map source — **162 nodes · 247 edges · 14 clusters** |
| [`nb_qsp_model.svg`](../../../neuroblastoma/nb_qsp_model.svg) / [`.png`](../../../neuroblastoma/nb_qsp_model.png) | Rendered map (`dot -Tpng -Gdpi=150`) |
| [`nb_mrgsolve_model.R`](../../../neuroblastoma/nb_mrgsolve_model.R) | **46-ODE** mrgsolve model + regimen builder + 9 scenarios + virtual cohort |
| [`nb_shiny_app.R`](nb_shiny_app.R) | **9-tab** Shiny dashboard |
| [`nb_references.md`](nb_references.md) | **79 references** (every PMID checked against the PubMed API) |

Rendering / running:

```bash
dot -Tsvg nb_qsp_model.dot -o nb_qsp_model.svg
dot -Tpng -Gdpi=150 nb_qsp_model.dot -o nb_qsp_model.png
Rscript -e 'source("nb_mrgsolve_model.R"); print(nb_scen_dose())'
Rscript -e 'shiny::runApp("nb_shiny_app.R")'
```

---

## 1. The two structural decisions

There are effectively only two things this model *chose*; everything else is a
consequence of them.

### Decision 1 — the PD driver is not occupancy but **the number of bound IgG per cell**

GD2 is expressed at about 8×10⁶ copies per cell. Converted to moles, that is
**13.3 nmol of antigen per gram of tumour**. One course of dinutuximab
(17.5 mg/m²/day × 4 days, BSA 0.8 m²) is **373 nmol**. So it is not that the
antigen greatly outnumbers the antibody — that is not the problem. What is
decisive is that **capillary permeability is overwhelmingly small relative to the
amount of antigen**.

The permeability-surface area product reported in the literature for monoclonal
antibodies in solid tumours is about 0.3–3 µL/(h·g) per gram, i.e. a range of
0.007–0.07 L/(d·kg). Taking the midpoint, 0.020 L/(d·kg), the antibody that can
be delivered per gram per day at a plasma concentration of 80 nM is only
**1.6 × 10⁻³ nmol**. Filling 13.3 nmol of antigen would take 8,000 days.

**Conclusion: in a solid tumour, antigen occupancy cannot exceed 1% at any
clinical dose.** Run the model and occupancy at clinical doses is of order 10⁻³.
Occupancy therefore cannot be the PD driver. ADCC responds to **how many** IgG
molecules are attached to the cell surface, and that threshold is of order
thousands to tens of thousands. The model therefore uses

```
BPC (bound IgG per cell) = bound antibody (nmol) × 6.022e14 / cell number
ADCC driving force = BPC / (BPC + ADCC50),   ADCC50 = 2 × 10⁴ /cell
```

This reformulation was forced during model construction, when **arithmetic done by
hand contradicted the occupancy-based expression originally intended** (see §6,
the list of defects).

### Decision 2 — three target compartments, three permeabilities, two Hill coefficients

There are three places the antibody has to reach, and they do not share a
permeability.

| Compartment | PS/g (L/(d·kg)) | Basis |
|------|------|------|
| Solid tumour | **0.020** | midpoint of the literature range + an interstitial fluid pressure (IFP) penalty that grows with volume |
| Marrow | **0.300** (×15) | sinusoidal vessels, no barrier — **assumption** |
| Dorsal root ganglion (DRG) | **4.000** (×200) | fenestrated capillaries, **no blood-nerve barrier** — **assumption** |

And the two Fc effector functions **have different Hill coefficients because their
stoichiometry differs**:

- **ADCC — Hill n = 1.** One Fc binding one FcγR suffices.
- **CDC — Hill n = 2.** C1q has to bridge **two adjacent** surface-bound Fc
  molecules (references 31–32). CDC is therefore **supralinear** in surface
  density.

Put these two decisions in, run the model, and at clinical dose the bound IgG per
cell is

| | Solid tumour | Marrow | DRG/nerve |
|---|---|---|---|
| BPC (molecules/cell) | **1.4 × 10⁴** | **3.0 × 10⁶** | **2.0 × 10⁵** |
| Relevant threshold | ADCC50 = 2×10⁴ | ADCC50 = 2×10⁴ | CDC50 = 5×10⁴ |
| Relative to threshold | **0.7× — below threshold** | **150× — saturated** | **4× — saturated** |

That is, **the therapeutic arm (marrow) is already saturated, the toxic arm
(nerve) is already saturated, and only the solid-tumour arm is stuck below its own
threshold.** This arrangement generates the entire clinical behaviour of the drug.

---

## 2. Nine derived axes

Every number is what the scenario functions in `nb_mrgsolve_model.R` actually
printed. Cross-checked against a mirror implementation in Python.

### Axis 1 — same drug, same dose, two different dose-response geometries

Dinutuximab was swept 1000-fold, from 0.175 to 175 mg/m²/day (5 cycles of
immunotherapy).

| Dose (mg/m²/day) | Cmax (nM) | BPC tumour | BPC marrow | ADCC AUC solid | ADCC AUC marrow | Pain AUC |
|---|---|---|---|---|---|---|
| 0.175 | 0.9 | 1.05e3 | 1.36e5 | 0.062 | 11.47 | **5.2** |
| 0.875 | 4.5 | 6.8e2 | 5.95e5 | 0.315 | 15.68 | 82.0 |
| 1.75 | 9.0 | 1.39e3 | 1.03e6 | 0.628 | 16.76 | 170.8 |
| 8.75 | 44.8 | 7.07e3 | 2.48e6 | 2.724 | 18.11 | 346.9 |
| **17.5 (clinical)** | **89.7** | **1.42e4** | **3.01e6** | **4.570** | **18.35** | **388.1** |
| 35 | 179.5 | 2.84e4 | 3.37e6 | 6.899 | 18.49 | 417.9 |
| 175 | 897.7 | 1.42e5 | 3.73e6 | 11.907 | 18.60 | 453.1 |

**Raise the dose 10-fold** above the clinical dose: solid-tumour ADCC goes up
**2.6-fold** (and 19-fold from 1.75 → 175, entirely linear), marrow ADCC
**1.014-fold**, pain **1.17-fold**. Conversely, **lower it 10-fold**: marrow ADCC
retains 91% (16.76/18.35) while pain exposure falls **56%** (170.8 vs 388.1). Go
all the way down to 0.175 mg/m² and 62% of marrow ADCC is retained while only
**1.3%** of the pain remains.

> **Derived conclusion.** The optimal dose depends on **where** the lesion is. If
> marrow MRD is the only problem, the clinical dose is already far up the
> diminishing-returns arm and additional dose buys nothing but pain. If solid
> lesions remain, escalation is right — but what forbids it is not the tumour, it
> is **the pain sitting on the other Fc arm**. Which is how the story moves on to
> K322A (§Axis 3).

### Axis 2 — delivery (permeability) and dose are perfectly interchangeable, but pain is not

| Tumour permeability | BPC tumour | ADCC AUC | log kill |
|---|---|---|---|
| ×0.10 | 1.40e3 | 0.630 | 0.00 |
| ×1 (reference) | 1.42e4 | 4.570 | 1.62 |
| ×2 | 2.83e4 | 6.892 | 2.62 |
| ×10 | 1.40e5 | 11.869 | 4.65 |
| ×25 | 3.40e5 | 13.440 | 5.17 |

Permeability ×2 is **identical to three decimal places** with dose ×2
(35 mg/m²: BPC 2.836e4, ADCC 6.899, logkill 2.62). Because the two quantities
appear in the model only as a product.

**Fate of the infused antibody** (1,867 nmol total over 5 cycles):

| | Bound (nmol) | % of dose |
|---|---|---|
| bound to the solid tumour deposit (2×10⁶ cells) | 8.0e-5 | 0.0000043% |
| bound to the marrow deposit (1×10⁶ cells) | 2.0e-3 | 0.00011% |
| bound to dorsal root ganglion (2 g of tissue) | **2.09** | **0.112%** |

> **Derived conclusion.** Improving delivery 25-fold is worth exactly as much as
> raising the dose 25-fold, and **it does not touch the pain compartment**. The
> development direction this model points at is not escalation but delivery
> (vascular normalisation, IFP reduction, local routes of administration). And
> 99.89% of what is infused reaches no target at all and is eliminated — the
> absolute molar comparison above is dominated by the mass of the target tissue
> (2 g of nerve vs 0.002 g of deposit), so it must be read together with the
> **per-cell** comparison (the table in §1). The two views say different things:
> absolute moles say where the drug is consumed, the per-cell numbers say where
> the effect comes from.

### Axis 3 — the K322A Fc mutation: the only lever that moves the therapeutic index

| | Pain AUC | Peak pain | Complement nadir (CPL) | ADCC AUC solid | ADCC AUC marrow |
|---|---|---|---|---|---|
| ch14.18 (dinutuximab) | 388.1 | 2.53 | 0.35 | 4.570 | 18.350 |
| hu14.18K322A (C1q ×0.1) | **86.7** | **0.59** | 0.84 | **4.570** | **18.350** |
| pain matched by dose instead: 0.914 mg/m²/day (**19.2-fold reduction**) | 86.7 | — | — | **0.329** | 15.756 |

- The Fc mutation reduces pain **4.48-fold** while holding ADCC **identical to
  three decimal places**.
- Obtaining the same pain **by dose reduction** requires a 19.2-fold cut, and the
  price is a **93% loss** of solid-tumour ADCC and a **14% loss** of marrow ADCC.

> **A non-obvious emergent result.** C1q binding was reduced 10-fold but pain fell
> only 4.48-fold. The reason is feedback: with CDC weakened, complement pool
> consumption falls, so available complement recovers from 0.35 → 0.84, and that
> recovery gives back more than half of the attenuation. This was not put in; it
> came out.
>
> **Clinical implication.** If pure marrow MRD is the only target, dose reduction
> is worth almost as much as Fc engineering (86% of marrow ADCC retained). If
> there is a solid component, **dose reduction is a disaster and only the Fc
> mutation works.** The two strategies are not the same thing, and what is being
> treated is the fork in the road.

### Axis 4 — why this is an MRD drug for **after** consolidation

| Tumour burden | BPC tumour | Peak ADCC (1/day) | log kill over the immunotherapy window |
|---|---|---|---|
| 0.001–0.1 g | 1.26e4 | 0.049 | **0.62** |
| 1 g | 1.24e4 | 0.048 | 0.60 |
| 7 g | 1.07e4 | 0.042 | 0.46 |
| 26 g (half-maximal point) | ~6.3e3 | — | — |
| 60 g | 2.99e3 | 0.015 | **0.00** |
| 300 g | 4.85e2 | 0.003 | **0.00** |

Below 1 g, BPC is **entirely independent of burden** (permeability, antigen and
interstitial volume are all proportional to mass, so they cancel). The burden
dependence comes only from the IFP term, and the **half-maximal point is about
26 g**; above 60 g the log kill is 0.

> **Derived conclusion.** The same drug at the same dose on the same schedule
> kills 0.62 log in MRD and kills **nothing** in a 60 g mass. The sequence
> induction chemotherapy → surgery → high-dose consolidation → **and only then**
> immunotherapy is not a convention, it is what this curve demands.

### Axis 5 — why IL-2 failed and GM-CSF works (SIOPEN HR-NBL1 reproduced)

| | Peak effector | NK | Treg | Exhaustion | ANC | ADCC AUC solid/marrow | log kill |
|---|---|---|---|---|---|---|---|
| Antibody alone | 20.2 | 200 | 40 | 0.00 | 4.0 | 3.742 / 15.202 | 1.249 |
| + GM-CSF | **43.7** | 213 | 40 | 0.00 | **9.2** | **4.629 / 18.574** | **1.646** |
| + IL-2 | 20.3 | 200 | **47** | 0.03 | 4.0 | 3.686 / 14.987 | 1.224 |
| + GM-CSF + IL-2 (COG) | 43.7 | 213 | 47 | 0.03 | 9.2 | 4.570 / 18.350 | 1.620 |

- **GM-CSF: marrow ADCC +22%** (15.20 → 18.57). The mechanism is cell **number** —
  ANC rises from 4.0 to 9.2, and in this model **ANC *is* the granulocyte ADCC
  effector**.
- **IL-2: marrow ADCC −1.4%** (15.20 → 14.99). The IL-2 concentration reached at
  clinical doses is far below the NK expansion EC50 (40) but above the Treg EC50
  (5) — because CD25 is the high-affinity receptor. So only Treg rises, 40 → 47,
  and the net effect is **slightly negative**.

> **An honest distinction.** The **sign** of the IL-2 null result comes from having
> put in the EC50 ordering (Treg < NK), so it is partly built into the structure.
> The **magnitude** (−1.4%) and the +22% for GM-CSF are outputs. This arrangement
> reproduces the "no benefit + more toxicity" that HR-NBL1 observed.

### Axis 6 — isotretinoin: exposure that destroys itself, and an effect that disappears if the sequence is wrong

**Exposure** (160 mg/m²/day split BID × 14 days):

| | Day 1 Cmax | Day 4 | Day 8 | Day 14 | AUC |
|---|---|---|---|---|---|
| capsule swallowed whole (F 1.00) | 2.54 | **3.67** | 3.07 | 2.96 µM | 43.9 |
| capsule opened, mixed with food (F 0.60) | 1.54 | 2.33 | 1.97 | **1.91 µM** | 28.1 |

CYP26A1 autoinduction makes exposure fall on its own after a day-4 peak. And
**opening the capsule and mixing it into food drops the day-14 peak concentration
to 1.91 µM, below the reported 2 µM target** — from changing not the drug but
**the way it is given**.

**Sequence** (same patient, same total dose):

| | Burden after induction | Nadir | Maximum differentiated fraction | Day of relapse |
|---|---|---|---|---|
| No isotretinoin | 0.00158 | 7.53e-6 | 0.00 | **774** |
| given after consolidation (standard) | 0.00158 | 3.60e-7 | 0.48 | **896** |
| given **concurrently** with induction | **0.0179** | 2.81e-6 | **0.99** | **808** |

> **Derived conclusion.** Isotretinoin buys 122 days when the sequence is right and
> only 34 days when given concurrently with induction — **72% of the effect is lost
> to a scheduling error.** The mechanism is visible: given concurrently, **99%** of
> the tumour moves into the differentiated (post-mitotic) pool, and that pool is
> invisible to S/M-phase-dependent cytotoxics. The post-induction burden is 11.3
> times worse. ADCC has no proliferation dependence, so it does not have the same
> problem, which is why concurrent use with immunotherapy is fine. **Those two
> sentences are the administration sequence of the actual protocol.**

### Axis 7 — ¹³¹I-MIBG: the dose stolen by carrier and by the patient's home medicines

18 mCi/kg (13,320 MBq), 20 kg, 25 g tumour:

| | Whole-body dose | mGy/MBq | **Tumour dose** | ANC nadir | Thyroid |
|---|---|---|---|---|---|
| carrier-free (NCA) | 2.67 Gy | 0.200 | **10.6 Gy** | 0.34 | 228 MBq |
| carrier-added ×10 | 2.67 Gy | 0.200 | 9.1 Gy (−14%) | 0.34 | 228 |
| carrier-added ×100 | 2.67 Gy | 0.200 | **4.2 Gy (−60%)** | 0.34 | 229 |
| NCA + labetalol not stopped | 2.67 Gy | 0.200 | **3.2 Gy (−70%)** | 0.34 | 229 |
| NCA + tricyclic not stopped | 2.67 Gy | 0.200 | 5.2 Gy (−51%) | 0.34 | 228 |
| NCA + KI thyroid blockade | 2.67 Gy | 0.200 | 10.6 Gy | 0.34 | **22.8** |

> **The key point.** Carrier addition and NET-blocking drugs lower **only the
> tumour dose** and do not change the whole-body dose (= marrow toxicity) by **a
> single digit**. It is a pure loss of therapeutic index. One tablet of labetalol
> that was not stopped takes away 70% of the tumour dose.

**Dosimetry-based prescribing** (2 Gy whole-body dose target):

| Biological clearance rate | Effective half-life | Activity required | mCi/kg | Tumour dose |
|---|---|---|---|---|
| 0.10 /day | 3.7 days | 6,729 MBq | **9.1** | 4.8 Gy |
| 0.19 /day (reference) | 2.5 days | 9,978 MBq | 13.5 | 7.5 Gy |
| 0.35 /day | 1.6 days | 15,755 MBq | **21.3** | 13.1 Gy |

> **Derived conclusion.** A fixed 18 mCi/kg prescription gives a slow clearer
> **twice** the target whole-body dose and gives a fast clearer **too little** —
> the case for dosimetry-based prescribing comes out as numbers. But there is a
> counterpart: fixing the whole-body dose makes **the tumour dose vary 2.7-fold**
> (4.8 → 13.1 Gy). In uniformising toxicity, dosimetry **creates variability in
> efficacy.** This was a result the model did not anticipate.

### Axis 8 — ALK inhibition: not a mutation problem but a **free fraction** problem

| | Free Cmax | Free trough | IC50 | Trough/IC50 | Mean inhibition |
|---|---|---|---|---|---|
| Crizotinib 280 mg/m² BID / R1275Q | 0.127 | 0.119 µM | 0.30 | **0.40** | 16.1% |
| Crizotinib / F1174L | 0.127 | 0.119 | 2.50 | **0.05** | 2.6% |
| **Lorlatinib** 95 mg/m² QD / R1275Q | 0.706 | 0.484 | 0.05 | **9.7** | **50.8%** |
| Lorlatinib / F1174L | 0.706 | 0.484 | 0.05 | 9.7 | 50.8% |
| Crizotinib at **double dose** / F1174L | 0.254 | 0.237 | 2.50 | 0.09 | 5.0% |

Even against a **sensitive** mutation, crizotinib's free trough is only 0.4 times
the IC50, so it inhibits just 16%. Lorlatinib, at 9.7 times, inhibits 51% (92% of
its own maximum). That 24-fold difference is made **roughly half and half by
protein binding (3.8-fold) and affinity (6-fold)** — not by one or the other.
Doubling crizotinib does not rescue F1174L (2.6% → 5.0%).

> That crizotinib produced a response in only 1 of 11 ALK-mutated neuroblastomas in
> ADVL0912 is predicted by this arithmetic.

### Axis 9 — sodium thiosulfate: the ear-versus-tumour trade-off as numbers

| | Cumulative ototoxicity | Burden after induction | Day of relapse |
|---|---|---|---|
| No STS | 0.060 | 0.00158 | **896** |
| STS, no tumour protection | **0.015** | 0.00158 | 896 |
| STS, 10% tumour protection | 0.015 | 0.00286 | 892 (−4 days) |
| STS, 20% tumour protection | 0.015 | 0.00585 | 858 (**−38 days**) |
| STS, 30% tumour protection | 0.015 | 0.0146 | 818 (**−78 days**) |

> **Derived conclusion.** A 75% reduction in ototoxicity remains a good deal **as
> long as tumour protection does not exceed 20%**. At 20% you pay 5 weeks, at 30%
> 11 weeks, of relapse-free time. This table quantifies what the concern raised by
> ACCL0431 means — except that `FTUMSTS` is not a measured value but **a scanned
> assumption**, and the experiment that would fix it (measuring intracellular
> platinum accumulation in disseminated lesions with and without STS) is the only
> unknown on this axis.

### Bonus axis — FcγR genotype: toxicity unchanged, efficacy differing by 1.6 log

| Genotype | ADCC AUC solid/marrow | log kill | Pain AUC |
|---|---|---|---|
| F/F (low affinity, 1.00) | 4.570 / 18.350 | **1.620** | 388.1 |
| V/F (1.35) | 6.179 / 24.772 | 2.331 | **388.1** |
| V/V (high affinity, 1.80) | 8.235 / 33.030 | **3.233** | **388.1** |

Because pain is complement-mediated, it is **entirely independent** of FcγR
genotype (388.1 in all three groups). A patient homozygous for the low-affinity
allele therefore **pays the toxicity of the full dose while receiving a fortieth
of the killing**. Genotype moves the therapeutic index 40-fold without changing
the dose.

---

## 3. Calibration — what was fitted and what came out

| Target | Target (literature) | Model output | Calibrated / predicted |
|---|---|---|---|
| Dinutuximab Cmax | ~11.5 µg/mL | **11.88 µg/mL** | calibrated (`V1AB`) |
| Dinutuximab terminal t½ | ~10 days | **10.4 days** | calibrated (`CLAB`) |
| Isotretinoin peak concentration | 2–4 µM | **2.5–3.7 µM** | calibrated (`VRA`,`CLRA`) |
| CYP26A1 autoinduction | exposure falls within a course | day 4 3.67 → day 14 2.96 µM | prediction |
| ¹³¹I-MIBG whole-body dose | 0.20 mGy/MBq | **0.200** | **calibrated** (`SWB`) — not a prediction |
| MIBG tumour dose | median ~15–30 Gy | **10.6 Gy** (25 g lesion) | calibrated (`VMAXNET`) |
| ANC nadir (induction) | grade 4 (<0.5) | **0.50** | calibrated (`SLOPE_CT`) |
| Platelet nadir | grade 3–4 | **59** | calibrated |
| EBRT 21.6 Gy surviving fraction | derived from LQ | 0.022 (α 0.15, α/β 10) | derived |
| **ANBL0032 control arm 2-year EFS** | **46%** | **40.0%** (n=20) | reproduced |
| **ANBL0032 immunotherapy 2-year EFS** | **66%** | **90.0%** (n=20) | **overestimate — see §4** |
| **HR-NBL1: added benefit of IL-2** | **none** | **−1.4%** | reproduced |
| **ADVL0912 crizotinib low activity** | 1 of 11 | free trough = 0.4 × IC50 | reproduced |
| K322A pain reduction | reduction reported | 4.48-fold reduction, ADCC unchanged | prediction |

Virtual cohort (n = 20, log-normal variability on parameters + FcγR genotype +
MYCN 45%):

| Arm | 1-year EFS | 2-year EFS | 3-year EFS |
|---|---|---|---|
| Isotretinoin alone | 95.0% | **40.0%** | 30.0% |
| Full COG (Ab + GM-CSF + IL-2) | 100% | **90.0%** | 60.0% |

At n = 20 the 95% confidence intervals on 2-year EFS are roughly 19–61% and
77–100% respectively. The control arm overlaps ANBL0032 (46%); the immunotherapy
arm does not.

---

## 4. Reported failures — what was left unfixed

What the model cannot do is not hidden. Each item carries a **structural
diagnosis**.

**F1. It overestimates the benefit of immunotherapy (90% vs 66% observed).** The
control arm is right but the immunotherapy arm is 24 percentage points too high.
**Diagnosis:** this model has only two mechanisms of immunotherapy failure —
antigen-low (MES) escape and poor delivery. Acquired resistance, ADCC-refractory
subclones and loss of exposure through HACA are all absent. So once ADCC exceeds
the regrowth rate, extinction **necessarily** follows. More than half of the
"relapse despite immunotherapy" seen in real practice has causes this model cannot
express.

**F2. The time course of pain differs from the clinical description.** The model
predicts pain reaching half its maximum 28 hours after the first infusion and its
maximum at 69 hours, decaying slowly thereafter (by day: 0.72 → 2.29 → 2.53 →
2.49 → 2.46). Clinically it is described as **worst at the first infusion**.
**Diagnosis:** the only attenuating mechanism in this model is complement
consumption (CPL 0.85 → 0.35), and that is not fast enough to make day 1 the peak.
Nociceptive sensitisation (wind-up) and an anticipatory pain component are absent
from the model.

**F3. Infusion time has no effect at all on pain — contradicting clinical
practice.** The pain AUC for a 0.75-hour and a 24-hour infusion are 388.2 versus
388.2 and the peaks 2.53 versus 2.53, i.e. **identical**. That is because the
nerve compartment reaches equilibrium in about 2 hours and CDC is already
saturated. Yet in reality the tolerability profiles of long infusions
(dinutuximab beta) and short ones (naxitamab) are reported to differ. **Both a
diagnosis and this model's most useful output:** if infusion rate really does
matter clinically, then one of three things is wrong —
(i) nerve GD2 density is higher than the assumed 0.5×10⁶/cell, so there is no
saturation,
(ii) `PSG_NRV` is lower than 4.0, so equilibration is slower,
(iii) pain responds not to equilibrium occupancy but to the **rate of complement
activation**.
**One experiment separates the three: is the pain score proportional to Cmax or to
AUC?** The model answers AUC, so if a Cmax dependence is observed then (iii) is
right and the CDC term in this model has to become a flux expression rather than an
equilibrium one.

**F4. It fails to reproduce the effect of consolidation intensity (and even
inverts the ordering).**

| | Nadir burden | Day of relapse |
|---|---|---|
| single ASCT, CEM | 3.60e-7 | **896** |
| single ASCT, BuMel | 4.35e-7 | **896** |
| tandem ASCT | 2.22e-7 | **880** |

HR-NBL1 found BuMel superior to CEM (3-year EFS 50% vs 38%) and ANBL0532 found
tandem superior to single (3-year EFS 61.6% vs 48.4%). The model produces **no
difference**, and for tandem **marginally worse**. **Diagnosis:** in this model the
time of relapse is governed not by the depth of consolidation but by **the nadir
and regrowth rate that immunotherapy sets**. The difference in nadir is less than
2-fold and that is worth only 16 days of relapse time. In other words, the model
**is missing what consolidation actually achieves** — the most plausible candidate
being eradication of chemotherapy-resistant subclones, which a single quiescent
pool (TQ) does not represent. Fitting the intensity effect as a parameter without
putting that in would be concealment rather than calibration, so it was left as it
is.

**F5. The MIBG tumour dose is below the reported median (10.6 Gy vs ~15–30 Gy).**
The value is for a 25 g lesion and dose per gram falls with mass, so the direction
is right; but `VMAXNET` is a value **calibrated** to hit the target dose and does
not rest on an independent measurement of NET expression. The carrier and
NET-blockade effects in Axis 7 are therefore trustworthy as **ratios** and not
trustworthy as **absolute Gy**.

**F6. Antigen-low (MES) escape has almost no effect.** Raising the MES fraction
from 0 all the way to **0.60** takes marrow ADCC AUC from 18.365 to 18.028
(**−1.8%**) and log kill from 1.620 to 1.615, i.e. essentially unchanged.

| MES fraction | 0.00 | 0.05 | 0.20 | 0.40 | 0.60 |
|---|---|---|---|---|---|
| ADCC AUC marrow | 18.365 | 18.350 | 18.296 | 18.194 | 18.028 |
| log kill | 1.620 | 1.620 | 1.619 | 1.618 | 1.615 |

**Diagnosis:** the marrow arm is saturated at 150 times its threshold, so cutting
mean density by more than half does not move the driving force. The model
therefore predicts that **antigen loss only starts to matter once density falls
close to the threshold** — that is, that the relationship between GD2 density and
response should be a **step** rather than a gentle slope. That is testable, but
because the antigen-low relapse seen clinically is not represented in this model,
it is one of the direct causes of F1 (overestimated immunotherapy effect). Real
escape is more likely to be **selection of a completely negative subclone** than a
fall in mean density, and that cannot be expressed by the structure of this model,
which has no discrete subpopulations.

**F7. Things not implemented.** Increased clearance from anti-drug antibody
(HACA) (present as a node on the map but with no ODE), eflornithine (DFMO)
maintenance, GD2 CAR-T, checkpoint inhibition, haploidentical NK infusion, the
kinetics of reversible ADRN⇄MES switching (represented only as a fixed fraction),
risk of second malignancy, second-look surgery. Morphine for the pain is a comment,
not a state variable.

---

## 5. Measurements that would overturn this model's conclusions (falsification)

It is made explicit that the model's two structural decisions rest on
**assumptions**.

1. **The permeability ordering (solid < marrow < nerve).** Only the solid-tumour
   value has a literature basis; marrow ×15 and nerve ×200 rest on nothing but the
   qualitative argument that "both are vascular beds without a barrier".
   **Decisive experiment:** measure %ID/g of radiolabelled anti-GD2 in marrow,
   nerve and solid lesions in the same subject. If the ordering inverts, every
   conclusion in §1–2 collapses.
2. **That pain is complement-mediated.** The basis is the **indirect** evidence of
   K322A (remove C1q and pain falls), not a direct measurement of complement at
   the nerve. **Decisive experiment:** the Cmax-versus-AUC experiment of F3, and
   whether pain falls when C5 inhibition (eculizumab or similar) is added.
3. **A nerve GD2 density of 0.5×10⁶/cell and 2 g of GD2+ nerve tissue.** These are
   effectively pure assumptions, and these two numbers decide whether the nerve arm
   is saturated. F3 and F6 hang on them.
4. **ADCC50 = 2×10⁴, CDC50 = 5×10⁴ molecules/cell.** The **ordering** of the
   thresholds (CDC requires the higher density) is a general principle in the
   antibody literature, but these are not values measured for GD2.
5. **`FTUMSTS`** (tumour protection by STS). The whole of Axis 9 sits on this
   unknown.

---

## 6. Defects found by running it

**D1. The occupancy-based PD expression was arithmetically impossible.** The first
implementation had ADCC ∝ occupancy and CDC ∝ occupancy². Working it out by hand:
against 13.3 nmol of antigen per gram, permeability delivers only 1.6×10⁻³ nmol
per day, so occupancy is permanently of order 10⁻³, and with that expression ADCC
is zero at any dose. The entire model was reformulated in terms of **bound IgG per
cell**. This is what became the model's central claim.

**D2. Bound antibody did not disappear along with the cells that died.** `BNDT`
had no cell-loss term, so as burden fell the bound IgG per cell diverged — at a
burden of 1,400 cells, BPC reached 6.9×10⁵/cell (17% of the antigen itself) and
ADCC became self-reinforcing. A term `min(g,0)·BNDT` removing bound antibody in
proportion to the fractional cell loss was added.

**D3. As the burden approached zero the tumour antibody subsystem became singular
and the integrator stalled.** `PS_T`, interstitial volume and antigen amount are
all proportional to cell number, so all three go to zero together. An extinction
floor of one cell (`NEXT = 1e-9`) was introduced.

**D4. The hard extinction cutoff and `min(gT,0)` made the right-hand side
non-smooth, and LSODA failed to converge in 1 of 10 cohort members**
(`Repeated convergence failures`). They were replaced by smooth equivalents — a
Hill-1 switch `N/(N+NEXT)` for growth and `min(x,0) ≈ ½(x − √(x²+ε))` for the
minimum. Results were preserved across the replacement (days of relapse 528/896
unchanged).

**D5. The NET Vmax for MIBG was 3,400-fold too large.** The entire infused MIBG
mass was trapped in the tumour and the tumour dose came out at **2×10⁶ Gy**.
Worse, in that state uptake is limited by mass rather than by the transporter, so
**the carrier competition effect itself disappeared** and the difference between
NCA and carrier-added ×100 looked like only 1.66-fold. Once Vmax was calibrated to
the reported tumour dose, that difference emerged as 2.5-fold (10.6 vs 4.2 Gy). A
wrong parameter was **erasing even a qualitative conclusion.**

**D6. The ALK inhibitor PK treated mg as µmol.** Free concentrations came out
100-fold too low (0.005 µM) and every inhibition was under 1%. Crizotinib
280 mg/m² BID is 497 µmol per dose and apparent clearance is ~100 L/h in adults
(722 L/d by allometry to 20 kg). After the fix, the free trough of 0.119 µM
matched the reported value and the conclusions of Axis 8 held.

**D7. The dosing schedule lookup was an O(n) linear scan, repeated at every
right-hand side evaluation.** One 900-day protocol took about 30 seconds. It was
precompiled into per-state step functions, bringing that under 1 second.

---

## 7. Structure

**46 ODEs**, units: days · 10⁹ cells · nmol/nM · µmol/µM · MBq · Gy · 10⁹/L

| Block | State variables |
|---|---|
| Tumour (5) | `TP` proliferating · `TQ` quiescent (G0) · `TD` differentiated (post-mitotic) · `TM` marrow MRD · `VTU` volume |
| Antibody (8) | `ABC` `ABP` systemic · `ABT`/`BNDT` solid tumour · `ABM`/`BNDM` marrow · `ABN`/`BNDN` nerve |
| Immune (5) | `NKB` `NKT` `TREG` `NKEXH` `CPL` (complement pool) |
| Myelosuppression (6) | `PROL` `TR1` `TR2` `TR3` `ANC` `PLT` (Friberg structure) |
| Retinoid (3) | `RAG` `RAA` `FIND` (CYP26A1 induction) |
| Cytotoxics (2) | `CTC` `CTP` |
| MIBG (7) | `AWB` `MBC` `MASS` (moles including cold carrier) `MBT` `DWB` `DTU` `THY` |
| ALK (2) | `AKG` `AKC` |
| Cytokines (4) | `IL2D` `IL2C` `GMD` `GMC` |
| Endpoints (4) | `PAIN` `HVA` `OTO` `PLATAUC` |

The four tumour pools have **four different drug sensitivities** — and this is what
forces the order of administration:

| | Cytotoxics | ADCC | MIBG β⁻ / EBRT | Retinoid |
|---|---|---|---|---|
| `TP` proliferating | **all of it** | yes | yes | → moves to `TD` |
| `TQ` quiescent | ×0.12 | yes | **yes** (proliferation-independent) | — |
| `TD` differentiated | **×0.02** | **yes** | yes | — |
| `TM` marrow | ×0.55 | **×1.25 (the strongest)** | yes | — |

Why `MASS` is a separate state: ¹³¹I decays but the carrier molecules do not, so
the **specific activity in blood changes with time**. Without this term, carrier
competition cannot be expressed.

---

## 8. Usage

```r
source("nb_mrgsolve_model.R")

# standard COG high-risk regimen versus isotretinoin alone
a <- nb_run(e = nb_regimen(immuno = FALSE), end = 1200)
b <- nb_run(e = nb_regimen(immuno = TRUE),  end = 1200)
nb_relapse(a); nb_relapse(b)

# the nine axes
nb_scen_dose()      # dose-response geometry of the two Fc arms
nb_scen_perm()      # permeability versus dose
nb_scen_k322a()     # the Fc mutation
nb_scen_burden()    # why it is an MRD drug
nb_scen_cytokine()  # IL-2 failure, GM-CSF success
nb_scen_retinoid()  # autoinduction · bioavailability · schedule
nb_scen_mibg()      # carrier competition and dosimetry
nb_scen_alk()       # free fraction
nb_scen_sts()       # ear versus tumour

# virtual cohort EFS
co <- nb_cohort(200)
mean(nb_efs(co, immuno = TRUE) > 730)
```

The Shiny dashboard offers the same axes interactively across nine tabs:
① patient · regimen ② pharmacokinetics ③ the two Fc arms ④ delivery barriers
⑤ tumour burden · PD ⑥ clinical endpoints ⑦ scenario comparison ⑧ MIBG dosimetry
⑨ virtual cohort (EFS)

```r
shiny::runApp("nb_shiny_app.R")
```

---

## ⚠️ Disclaimer

This model is for **research and educational purposes**. It must not be used for
clinical care, dosing decisions or protocol design. A substantial number of the
parameters are **assumptions** rather than literature-based, and §5 lists them
together with the measurement that would falsify each. The list of failures in §4
is a set of clinical facts the model **cannot reproduce**, and among them F1 and F4
are limitations severe enough to rule out any use of this model for prognostic
prediction.

This model is for **research and educational purposes only** and must not be used
for clinical care, dosing decisions, or protocol design. Many parameters are
assumptions rather than measured values (§5, with the falsifying measurement named
for each), and §4 lists clinical facts the model fails to reproduce — F1 and F4 in
particular rule out any prognostic use.

79 references: [`nb_references.md`](nb_references.md)
