# Bronchopulmonary dysplasia (BPD) — QSP model

> **In one line.** A QSP model that formalises BPD not as "a disease in which injury
> accumulates and the lung is wrecked" but as a **scheduling problem**.
> Alveolarisation is a **scheduled developmental programme** running on the
> postmenstrual age (PMA) clock, and injury does not destroy alveoli — it **closes
> the growth gate G and halts the programme**. Because the integral of the window
> `S_dev(PMA)` is finite, **a day spent halted becomes irrecoverable surface area**.
> The central prediction of this model is therefore a single one — **timing beats
> potency.**

---

## 0. Deliverables

| File | Contents |
|------|------|
| [`bpd_qsp_model.dot`](bpd_qsp_model.dot) · [SVG](bpd_qsp_model.svg) · [PNG](bpd_qsp_model.png) | Mechanistic map — 227 nodes / 16 clusters / 382 edges |
| [`bpd_mrgsolve_model.R`](bpd_mrgsolve_model.R) | mrgsolve model — 37 ODEs, 19 scenarios, PMA-maturing neonatal PK |
| [`bpd_shiny_app.R`](bpd_shiny_app.R) | Shiny dashboard — 10 tabs |
| [`bpd_references.md`](bpd_references.md) | 127 references (all verified by actual lookup through the PubMed E-utilities) |
| [`bpd_reference_impl.py`](bpd_reference_impl.py) | **A dependency-free Python re-implementation** — the code that actually computes every number below |

Every number appearing below is output from `python3 bpd_reference_impl.py`. None of them
were copied from memory or from the literature, and where the README and the script output
disagree, **the script is right and the README is out of date**.

---

## 1. Why a "scheduling model" rather than an "injury model"

Most BPD models in the literature are injury models. Oxygen and pressure damage the lung,
the damage accumulates, and the lung is left scarred. That framework predicts **"give a
more potent anti-injury drug at any time and things improve"**. Forty years of clinical
trials say otherwise. DART helped extubation but did not reduce BPD at 36 weeks, and
PREMILOC, stimulating the same receptor far more weakly, started on day 1 and improved
survival without BPD. Late high-dose dexamethasone reduced BPD but increased cerebral palsy
and was effectively abandoned.

This model uses the following structure instead.

```
    dALV/dt  =  k · S_dev(PMA(t)) · G(t) · (1 − ALV/ALVmax)  −  (elastolytic loss)

    G(t)  =  [ VEGF^wv · cGMP^wn · IGF1^wi · RA^wr · EPC^we ]
             ─────────────────────────────────────────────── · (1 − E_GR,antiprolif)
                     (1 + a·TGFβ1) · (1 + b·IL-1β)

    dLOSTW/dt  =  k · S_dev(PMA(t)) · (1 − G(t)) · (1 − ALV/ALVmax)
```

Three things follow immediately from this structure.

1. **Timing beats potency.** Because `S_dev` is decaying, a drug that opens the gate 40%
   on day 1 buys more surface area than one that opens it 100% on day 21. **The crossover
   day is computable.**
2. **Treatments fall into three classes by "which term they touch", and it is that class,
   not the potency, that sets the ceiling.** (a) exposure reduction (b) transmission
   inhibition (c) programme support. A drug that touches none of the three terms moves the
   numbers for a few hours and changes nothing at 36 weeks. The designated control of this
   model is furosemide, and it was put in deliberately.
3. **The readout is not the level of support but the lost window (LOSTW).** Two babies can
   have the same FiO2 at 36 weeks with entirely different LOSTW, and what predicts the
   adult FEV1 plateau is LOSTW, not FiO2.

---

## 2. Model structure (37 ODEs)

| Layer | Compartments |
|----|------|
| **PK (11)** | Caffeine depot · central (PMA-maturing CL), dexamethasone, hydrocortisone, budesonide (lung · systemic), hepatic retinol store · plasma retinol, azithromycin, sildenafil, cumulative dexamethasone-equivalent ledger |
| **Inflammation (5)** | NF-κB, IL-1β, neutrophils/elastase, TGF-β1, ROS |
| **Growth gate (5)** | VEGF, NO bioavailability, IGF-1, retinoate signalling, endothelial progenitor cells |
| **Structure (5)** | Alveolar surface area, microvascular density, septal thickening, airway smooth muscle, lung water |
| **Support (3)** | Support intensity, time fraction of invasive ventilation, Ureaplasma load |
| **Pulmonary vasculature (2)** | Vascular remodelling, right ventricular hypertrophy |
| **Accumulators (6)** | Body weight, **LOSTW**, neurodevelopmental cost, ventilator days, oxygen days, mortality hazard |

Three key design decisions:

- **One receptor, three ligands.** Dexamethasone, hydrocortisone and budesonide are all
  converted to a **dexamethasone-equivalent concentration (mg/L)** and put through a single
  occupancy equation. Three different Emax values then branch off it — NF-κB inhibition
  (the benefit), compliance/extubation (what DART measured), and **an antiproliferative
  effect on septation itself (the hidden cost)**. Without that last term the model predicts
  "steroids, the earlier and the more the better", which is precisely the prediction the
  trials refuted.
- **Body weight is a state variable.** It sets the mg/kg dose, scales every clearance, and
  poor growth returns to the gate via IGF-1. This is how the loop is represented in which
  a steroid buys lung while suppressing growth, which lowers IGF-1, which closes again the
  gate the steroid was opening.
- **The oxygen requirement is split into two terms.** A *diffusion/reserve* term (surface
  area deficit, with a reserve range below `DTHR` in which FiO2 does not move at all) and a
  *shunt* term (atelectasis · oedema · fibrosis, which has no reserve and raises FiO2
  immediately). This separation is the quantitative basis of the claim that "FiO2 at 36
  weeks is a poor marker of developmental deficit".

---

## 3. Computed results

Every table and number below is output from `python3 bpd_reference_impl.py`. The reference
infant is deliberately in bad circumstances — **25 weeks' gestation, histological
chorioamnionitis, delivery-room intubation, a significant PDA for 14 days, no adjunctive
therapy.** Because this is the population the landmark trials actually recruited.

### A0. PK consistency — maturation of caffeine clearance (reproducing the literature anchors)

| PMA (weeks) | CL (L/h/kg) | t½ (h) | Css, citrate 5 mg/kg/d (mg/L) |
|---|---|---|---|
| 26 | 0.0065 | 90.3 | 16.0 |
| 28 | 0.0087 | 67.9 | 12.0 |
| 32 | 0.0143 | 41.1 | 7.3 |
| 44 | 0.0430 | 13.7 | 2.4 |
| 60 | 0.0958 | 6.1 | 1.1 |

It reproduces the literature anchors exactly (t½ ≈ 90–110 h in the extremely preterm infant
→ ≈ 5–6 h at term-corrected age; trough 10–20 mg/L on citrate 5 mg/kg/d). Peak concentration
after a 20 mg/kg citrate load, **14.8 mg/L**.

**One receptor, two ligands.** DART dexamethasone peak concentration 0.0416 mg/L → pulmonary
GR occupancy **0.912**. PREMILOC hydrocortisone peak concentration 0.294 mg/L → pulmonary GR
occupancy **0.746**. That is, **DART has the higher occupancy.** The reason it nonetheless
loses on the 36-week endpoint is that it arrives late.

### A1. The window — the programme still remaining at birth

Total integral of S_dev = 24.724 week-units; integral to 36 weeks PMA = 11.124.

| Gestational age | Consumed in utero | Remaining ex utero | Of which before 36 weeks (of the total) | Alveolar surface area at birth |
|---|---|---|---|---|
| 24 | 4.2% | 95.8% | **40.8%** | 0.123 |
| 25 | 6.0% | 94.0% | 39.0% | 0.168 |
| 28 | 14.0% | 86.0% | 31.0% | 0.354 |
| 32 | 29.0% | 71.0% | **16.0%** | 0.683 |

A 24-week infant must complete **40.8% of the alveolarisation programme inside the NICU
before the 36-week clock stops**, whereas a 32-week infant has only 16.0% left — a
**2.6-fold** difference. The gestational-age gradient comes out of this geometry rather than
from a fitted risk coefficient.

### A2. Timing vs potency — the crossover day = **14 days**

| Start day | 40% effect | 100% effect | Does 100% beat "40% on day 1"? |
|---|---|---|---|
| 1 | 0.7732 | 0.8036 | yes |
| 7 | 0.7582 | 0.7891 | yes |
| 10 | 0.7491 | 0.7802 | yes |
| **14** | 0.7352 | **0.7659** | **NO** |
| 21 | 0.7124 | 0.7413 | NO |
| 35 | 0.6737 | 0.6953 | NO |

Untreated, 0.6168. **The 40% drug started on day 1 (0.7732) is better than the perfect drug
started on day 14 (0.7659).** In potency it is a 2.5-fold difference, and 13 days of delay
erases it.

### A3. Is caffeine a lung drug or a ventilator drug?

| Arm | ALV36 | Ventilator days | grade | LOSTW |
|---|---|---|---|---|
| No caffeine | 0.6168 | 46.0 | 2 | 0.4203 |
| Caffeine, both mechanisms | 0.6729 | 29.7 | 1 | 0.3605 |
| Caffeine, direct A2A only | 0.6427 | 42.0 | 2 | 0.3921 |

Of the total benefit of +0.0560, **the direct A2A mechanism contributes +0.0259 (46%)** and
**exposure reduction +0.0301 (54%)**. Ventilation is shortened by **16.4 days**. That is, in
this model the majority of the caffeine benefit is not "a drug that is good for the lung" but
**"a drug that gets the baby off the ventilator sooner"**. A caffeine analogue without the
respiratory-centre effect should lose half the benefit — a falsifiable prediction.

Virtual population (a 90-point deterministic grid, 25–29 weeks' gestation, chorioamnionitis
35%):

| | Moderate–severe (grade 2–3) BPD | Death |
|---|---|---|
| Placebo | 11.9% | 9.5% |
| Caffeine | 5.3% | 7.2% |

**The model's limitation is stated as it is:** "requires oxygen at 36 weeks" (any BPD)
**saturates at 100%** in this population and so carries no information, and consequently the
treatment effect on that item is structurally exactly zero. This is because the baseline
practice pattern is deliberately bad. CAP's control rate was 47.2%, so **the model
overestimates the mild end.** The comparable value is the effect on grade 2–3 (−6.6
percentage points), and presenting the saturated column as if it were a "negative result"
would be the dishonest choice.

### A4. Steroids — the same receptor, three schedules, three conclusions

| Arm | ALV36 | grade | Ventilator days | Cumulative DEXEQ (mg/kg) | NDI cost | Net utility |
|---|---|---|---|---|---|---|
| No steroid | 0.6168 | 2 | 46.0 | 0.000 | 0.263 | −0.118 |
| PREMILOC HC d1–10 | **0.7048** | 1 | 20.3 | 0.340 | 0.440 | −0.198 |
| DART dex d14–23 | 0.6552 | 1 | 31.7 | 0.890 | 0.738 | −0.332 |
| Early high-dose dex d1–7 | 0.6765 | 1 | 26.3 | **3.002** | **2.796** | **−1.258** |
| Intratracheal budesonide | 0.6842 | 1 | 25.4 | 0.100 | 0.180 | −0.081 |

- **DART's extubation effect is real:** the time fraction of invasive ventilation goes from
  0.835 on day 14 to 0.175 on day 24. And yet it buys less 36-week surface area than
  PREMILOC (0.6552 against 0.7048). *Through the same receptor equation*, the arm with the
  higher occupancy loses.
- Population-level "survival without BPD": untreated 0.0% · **PREMILOC 31.5%** · DART 3.7% ·
  early high-dose dex 7.1% · budesonide 5.7%.
- **Why early high-dose dexamethasone was abandoned comes out as arithmetic:** it does buy
  lung (0.6765 > 0.6168), but the cumulative dexamethasone equivalent is 3.0 mg/kg, most of
  it given in the first week of life, so the NDI cost rises to 2.80 and the net utility of
  −1.258 is the worst of all arms.

### A5. Furosemide — the numbers that move and the endpoint that does not

| Day of life | FiO2 (reference) | FiO2 (furosemide) | Difference |
|---|---|---|---|
| 10 (dosing starts) | 0.473 | 0.473 | 0.000 |
| 12 | 0.488 | 0.469 | **−0.020** |
| 14 | 0.499 | 0.480 | −0.020 |
| 28 | 0.476 | 0.473 | −0.002 |

Alveolar surface area at 36 weeks 0.6168 → 0.6195 (**+0.44%**), grade 2 → 2.
**It improves oxygenation within 48 hours and leaves the developmental deficit essentially
untouched.** Because it reaches none of the three terms. It is this model's designated
"placebo-shaped drug" and was included on purpose.

### A6. Azithromycin — not a BPD drug but a subgroup drug

| | ALV36 | grade | Ventilator days |
|---|---|---|---|
| Ureaplasma +, not treated | 0.5109 | 3 | 61.7 |
| Ureaplasma +, treated | 0.6163 | 2 | 46.2 |
| Ureaplasma −, not treated | 0.6168 | 2 | 46.0 |
| Ureaplasma −, treated | 0.6168 | 2 | 46.0 |

The benefit in colonised infants is **+0.1054** and in non-colonised infants **+0.0000**.
An unselected trial dilutes the first number with the second. This is a **clinical trial
design conclusion** produced by a mechanistic model.

### A7. BPD-PH — and the point at which the model disagrees with the trial record

PVR = (1/CAP_rel^1.2) × (1 + tone) × (1 + remodelling).

| | PVR36 | Alveolar surface area | Microvascular density |
|---|---|---|---|
| Reference | 1.845 | 0.6168 | 0.7138 |
| Sildenafil (from day 28) | 1.602 (−13.2%) | 0.6803 (**+10.3%**) | 0.7987 (+11.9%) |
| Early bundle (no PDE5i) | 1.498 (−18.8%) | 0.7825 (+26.9%) | 0.8259 (+15.7%) |

⚠️ **Here the model disagrees with the clinical trials, and that is not hidden.**
The intended claim was that "a PDE5 inhibitor moves only the tone term and cannot touch
structure". **The numbers above do not say that.** Because cGMP is a weighted member of the
growth gate G (WN = 0.2), raising it also raises alveolarisation, and the model credits
sildenafil from day 28 with +10.3% of surface area. NO CLD and the preterm iNO trials do not
support this. So one of two things is true and the model does not yet know which — (a) WN is
too large (NO/cGMP permits angiogenesis but is not the rate-limiting step for septation in
an injured preterm lung), or (b) iNO and oral PDE5 inhibition do not in fact raise cGMP in
the compartment that matters (that is, the PD is right and the exposure assumption is
wrong). A sensitivity analysis setting `WN = 0.05` is the obvious way to test this.

### A8. SpO2 targets — the trade-off exists only "when there is hyperoxia to be saved"

**Severely affected infant (25 weeks, chorioamnionitis, intubated, caffeine):**

| Target | FiO2 at 36 weeks | ALV36 | Survival probability | PVR36 |
|---|---|---|---|---|
| 85–89% (low) | 0.272 | **0.6842** | 0.5998 | 3.525 |
| 91–95% (high) | 0.314 | 0.6729 | **0.8705** | 1.719 |

**Mildly affected infant (28 weeks, no chorioamnionitis, LISA + caffeine):**

| Target | FiO2 at 36 weeks | ALV36 | Survival probability |
|---|---|---|---|
| 85–89% (low) | 0.217 | 0.8244 | 0.7441 |
| 91–95% (high) | 0.221 | **0.8331** | **0.9416** |

In the severely affected infant **the lung optimum and the survival optimum diverge** —
exactly the shape NeOProM · SUPPORT · BOOST-II reported (a between-band FiO2 spread of
0.041). In the mildly affected infant, however, the trade-off **disappears**: already near
room air, there is no hyperoxia for the low band to save (an FiO2 spread of 0.004), only the
hypoxia cost remains, and the high band wins in every column. That is, this model predicts
that **a permissive saturation policy has to be targeted at infants who actually have an
oxygen requirement and is a net cost in those who do not**. The trials randomised the policy
but did not randomise the policy × severity interaction, so they could not have generated
this prediction. It is the most testable claim in this file.

### A9. The escalating loop — it amplifies but does not bifurcate

The antenatal hit (KANTE) is swept twice: with the loop live, and with the support-mediated
injury terms switched off (KVILI = KOX = 0).

| KANTE | ALV36 (loop on) | ALV36 (loop off) | FiO2 | grade |
|---|---|---|---|---|
| 0.0 | 0.7240 | 0.8087 | 0.251 | 1 |
| 1.2 | 0.6496 | 0.7874 | 0.342 | 2 |
| 2.4 | 0.6133 | 0.7751 | 0.400 | 2 |
| 3.0 | 0.6022 | 0.7703 | 0.416 | 2 |

Surface area lost across the whole sweep: 0.1218 with the loop, 0.0385 without →
**an amplification factor of 3.17.** The largest single-step drop is only 0.0436, so **the
response is smooth and monotonic.** In this parameterisation the escalating loop amplifies
but **does not bifurcate (the loop gain is subcritical).** That is, "evolving BPD" is not a
switch but a gradual slide; there is no particular day on which the baby is lost, and no
threshold to wait for before intervening. This is a real result, so it is reported rather
than hidden.

### A10. Full scenario table (36 weeks PMA)

| Scenario | ALV36 | % of ideal | LOSTW | FiO2 | gr | Vent days | PVR | NDI | Net utility |
|---|---|---|---|---|---|---|---|---|---|
| S00 Reference (no adjunctive therapy) | 0.6168 | 61.7 | 0.4203 | 0.394 | 2 | 46.0 | 1.845 | 0.263 | −0.118 |
| S01 Caffeine (CAP dose) | 0.6729 | 67.3 | 0.3605 | 0.314 | 1 | 29.7 | 1.719 | 0.172 | −0.078 |
| S02 Caffeine, respiratory-centre mechanism OFF | 0.6427 | 64.3 | 0.3921 | 0.363 | 2 | 42.0 | 1.808 | 0.247 | −0.111 |
| S03 LISA / nCPAP | 0.7089 | 70.9 | 0.3245 | 0.247 | 1 | 8.4 | 1.619 | 0.112 | −0.051 |
| S04 VTV + permissive hypercapnia | 0.6647 | 66.5 | 0.3705 | 0.318 | 1 | 37.6 | 1.727 | 0.079 | −0.035 |
| S05 Early PDA closure (day 3) | 0.6458 | 64.6 | 0.3897 | 0.348 | 2 | 37.9 | 1.763 | 0.230 | −0.104 |
| S06 PREMILOC HC d1–10 | 0.7048 | 70.5 | 0.3298 | 0.263 | 1 | 20.3 | 1.633 | 0.440 | −0.198 |
| S07 DART dex d14–23 | 0.6552 | 65.5 | 0.3836 | 0.324 | 1 | 31.7 | 1.753 | 0.738 | −0.332 |
| S08 Early high-dose dex d1–7 | 0.6765 | 67.7 | 0.3595 | 0.299 | 1 | 26.3 | 1.698 | 2.796 | −1.258 |
| S09 Intratracheal budesonide | 0.6842 | 68.4 | 0.3508 | 0.290 | 1 | 25.4 | 1.680 | 0.180 | −0.081 |
| S10 Azithromycin, Ureaplasma + | 0.6163 | 61.6 | 0.4209 | 0.395 | 2 | 46.2 | 1.847 | 0.264 | −0.119 |
| S11 Azithromycin, Ureaplasma − | 0.6168 | 61.7 | 0.4203 | 0.394 | 2 | 46.0 | 1.845 | 0.263 | −0.118 |
| S12 Vitamin A (Tyson) | 0.6264 | 62.6 | 0.4099 | 0.382 | 2 | 44.4 | 1.830 | 0.256 | −0.115 |
| S13 Furosemide alone | 0.6195 | 62.0 | 0.4175 | 0.390 | 2 | 45.2 | 1.837 | 0.260 | −0.117 |
| S14 Sildenafil (from day 28) | 0.6803 | 68.0 | 0.3536 | 0.318 | 1 | 39.0 | 1.602 | 0.235 | −0.106 |
| **S15 Early bundle (day 1)** | **0.8099** | **81.0** | **0.2128** | 0.225 | **0** | 1.0 | 1.456 | 0.334 | **+0.782** |
| S16 Late bundle (same drugs, day 21) | 0.7170 | 71.7 | 0.3174 | 0.265 | 1 | 26.2 | 1.652 | 0.684 | −0.308 |

**Two things are the point of this table.**

1. **S15 against S16 — the same drugs, 20 days apart.** Surface area relative to ideal 81.0%
   against 71.7% (**a 9.3 percentage point difference**), grade 0 against 1, net utility
   +0.782 against −0.308. Of all 19 arms, **the early bundle is the only one with a positive
   net utility.**
2. **An irreducible floor exists.** Even the best realistic bundle **fails to recover 19.0%
   of the ideal.** A virtual infant with mechanical and oxidative injury set **entirely to
   zero** still only reaches **90.9%** of the ideal. Because part of the deficit occurred
   before birth and part of the window is already closed. A model in which the best bundle
   recovers 100% would be a model that is lying.

---

## 4. Limitations, stated honestly

- **This is not a model fitted to individual patient data.** Twelve coefficients were
  calibrated by coordinate descent so as to reproduce the **direction and approximate
  magnitude** of the landmark trials (CAP, PREMILOC, DART, Tyson, NeOProM and others). The
  calibration targets and the residuals are all visible if you run
  `bpd_reference_impl.py`.
- **It does not agree completely with the network meta-analysis.** This model's
  "the earlier the better" timing narrative fits the sequence PREMILOC (day 1, positive) →
  SToP-BPD (days 7–14, negative) → Watterberg 2022 (days 14–28, negative) well, but the
  Cochrane network meta-analysis rates mid-to-late dexamethasone as rather effective on BPD
  endpoints. The timing claim is therefore **a mechanistic hypothesis not contradicted by
  the trial record**, and not a fact the trials established.
- **The width of the rising limb of `S_dev` (`SDW_L = 6` weeks) is a modelling choice.** The
  onset of human secondary septation is itself under debate. Narrowing this value steepens
  the gestational-age gradient and widening it flattens it. Experiment A1 has exposed this
  assumption, so do not trust it — change it and see.
- **The neurodevelopmental cost weight `W_NDI` is a value judgement, not a measurement.** The
  model can compute how much lung a steroid buys and how much dexamethasone equivalent it
  spends, but it cannot compute the exchange rate between one ventilator day and one point of
  developmental index. Move the slider on tab 10 of the Shiny app and the ranking inverts,
  and **that instability is the honest answer**.
- **The sildenafil PK is extrapolated from term-infant population PK.** The conclusions of A7
  depend on the *structure* of the PVR equation and not on this PK, but do not use the
  absolute concentrations anywhere.
- R was not installed in this session's environment. Consequently **`bpd_mrgsolve_model.R`
  has not been verified by execution.** Instead `bpd_reference_impl.py`, which
  re-implements the same equations using only the standard library, was run and verified, and
  all the numbers above came from there. The R file and the Python file were written so that
  `$PARAM` / `P0` and `$ODE` / `rhs()` correspond one to one.

---

## 5. Usage

```bash
# all experiments (A0 to A10)
python3 bpd_reference_impl.py

# individual experiments
python3 bpd_reference_impl.py A2      # timing vs potency
python3 bpd_reference_impl.py A5      # furosemide

# re-rendering the mechanistic map
dot -Tsvg bpd_qsp_model.dot -o bpd_qsp_model.svg
dot -Tpng -Gdpi=150 bpd_qsp_model.dot -o bpd_qsp_model.png
```

```r
# mrgsolve
source("bpd_mrgsolve_model.R")
BPD_simulate_scenarios()          # endpoint table for the 19 scenarios
BPD_window_remaining()            # fraction of the programme remaining at birth
BPD_caffeine_decomposition()      # decomposition of the caffeine benefit
BPD_caffeine_pk_check()           # PK consistency check
BPD_timing_vs_potency()           # crossover-day sweep

# Shiny
shiny::runApp("bpd_shiny_app.R")
```

---

## ⚠️ Disclaimer

This is a mechanistic model for educational and research purposes. It has not been fitted or
validated against individual patient data, has not been prospectively qualified, and has not
been reviewed by any regulatory agency. In particular, the neonatal doses quoted here are
**there to explain the model parameters** and must not be used for prescribing. Nothing in
this directory should be used in the care of a real infant.
