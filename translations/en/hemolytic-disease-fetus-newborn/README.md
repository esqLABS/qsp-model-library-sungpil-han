# Hemolytic Disease of the Fetus and Newborn (HDFN)

**QSP model — mechanistic map · 45-ODE mrgsolve model · 12-tab Shiny app · 96 PubMed references**

---

## The One Sentence

The rate at which fetal red cells are destroyed is a **product of three
factors**, and every treatment for this disease owns **exactly one** of them.

```
destruction  =  A(t)            x   f_ag(t)              x   M(t)
                antibody              antigen-positive        reticuloendothelial
                actually delivered    fraction of              clearance
                to the fetus          circulating haemoglobin  capacity
```

| Factor | What it is | Who moves it |
|------|----------|-----------------|
| **A** | Maternal titre x placental conveyor capacity (e^0.268 per week, ~210-fold from 16 to 36 weeks) | Nipocalimab (FcRn blockade), plasma exchange, high-dose IVIG (FcRn competition), and upstream of all of them, anti-D prophylaxis |
| **f_ag** | Fraction of circulating red cells that carry the antigen | **Intrauterine transfusion (IUT)** — donor red cells are antigen-negative, so transfusion does not add haemoglobin so much as **dilute the substrate** |
| **M** | Clearance capacity of splenic and hepatic macrophages, and its saturation | High-dose IVIG (FcγR blockade), splenic volume |

f_ag is the crux of this map. Because of this factor — which nobody names —
the rate of post-transfusion haemoglobin decline flattens with each
successive procedure (model: 0.29 -> 0.18 -> 0.17 g/dL/day, averaged over 10
surveillance seeds). It is why a transfusion's effect outlasts the
haemoglobin it delivered.

---

## Two Clocks

1. **The conveyor clock — e^0.268 per week.** Placental FcRn transport
   capacity increases roughly 210-fold between 16 and 36 weeks. After being
   diluted into a fetus that is itself growing, what the fetus actually
   experiences — the fetal:maternal antibody concentration ratio — becomes
   about 17-fold higher at term than at 19.5 weeks. The same maternal titre
   becomes an entirely different exposure depending on gestational age.
   *Early-onset severe HDFN is therefore not a different disease* — it is
   simply a titre high enough to be dangerous even while the conveyor is
   still weak. This is also why FcRn blockers must start at 14-16 weeks:
   not to treat the disease, but to intercept the conveyor before it can
   carry a dangerous dose.

2. **The clearance clock — a discontinuity at delivery.** In utero,
   bilirubin is cleared through the placenta — that is, through the
   *maternal liver* — with a half-life of about 1 hour. So a fetus with 8
   g/dL of haemolysis is anaemic but not jaundiced. At delivery that route
   is severed, and UGT1A1 activity at roughly 6% of adult level takes its
   place. **The haemolysis is unchanged, but the disease changes** — because
   one number in the denominator has changed.

3. And one slow state variable: maternal IgG in the neonate has a half-life
   of 3 weeks. Haemolysis outlives delivery by 6 weeks — the point at which
   a marrow suppressed by transfusion is least able to answer.

---

## What Was Fitted and What Was Predicted

**Only eight numbers were fitted.** All are published summary statistics, and
they were fitted in four mutually independent stages (`hdfn_calibrate.py`,
logged in [`calib.log`](../../../hemolytic-disease-fetus-newborn/calib.log)).

| Stage | Parameters | Target |
|------|----------|------|
| 1 | `v0`, `g_pl` | Fetal:maternal total IgG ratio 0.075 (19.5 wk), 1.25 (39 wk) — Malek 1996 |
| 2 | `kops` | First IUT at 26 weeks at anti-D 15 IU/mL — Nishie 2012 cohort mean 26.1 wk |
| 3 | `ksens`, `fmh_ante` | 16% sensitisation with no prophylaxis, 1.6% with 300 µg postpartum only |
| 4 | `ugt_birth`, `ugt_t50` | Physiological jaundice peak of 8 mg/dL at ~4 days in **healthy** term neonates |
| — | `emh_thresh` | Positioned so overt ascites appears at 5-6 g/dL |

`Kres`, `dev`, `nip_pl_pen`, `do2_alpha` (=1), and `visc_k` (=0) are fixed
**structurally** rather than fitted, and reported instead as sensitivities.

What the fitter never saw — i.e., predictions:

### 1. The **slowing** of the post-transfusion decline rate is a derived result

| Interval | Model | Reported |
|------|------|--------|
| After 1st IUT | 0.291 g/dL/day | **0.40 (SD 0.25)** — Nishie 2012 |
| After 2nd IUT | 0.183 g/dL/day | (not measured) |
| After 3rd IUT | 0.171 g/dL/day | (not measured) |

**The ordering is predicted, and the magnitude is under-predicted by about a
third** (roughly 0.5 SD). Reported without hiding it, decomposed into three
terms for the first interval.

| Term | Contribution |
|----|----------|
| Destruction of antigen-positive cells (4.55 g Hb) | +0.345 g/dL/day |
| Ageing of donor red cells (1.41 g Hb) | +0.107 g/dL/day |
| Volume change (114 -> 101 mL) | **-0.144** g/dL/day |
| Total | 0.308 g/dL/day |

The sign of the third term is the opposite of what one would expect: because
the volume of destroyed red cells shrinks faster than plasma volume
increases, fetoplacental blood volume actually **decreases**, so dilution
does not help the decline — it **fights** it by 0.14 g/dL/day. Shortening
donor red-cell lifespan from 70 to 45 to 30 days raises the decline rate
from 0.290 to 0.322 to 0.359, closing only about half of the remaining gap.
So donor-cell survival alone does not explain it, and at least one more term
is missing — candidates are bystander loss of donor cells, sequestration
within 24 hours of the procedure, or a true interval shorter than the model
captures. The model keeps the adult value of 70 days and reports the gap.

### 2. The MCA-PSV threshold of 1.5 MoM is not an empirical cutoff

If cerebral oxygen delivery is defended (blood flow proportional to 1/Hb),
then velocity x haemoglobin is constant, so

```
PSV MoM = 1 / Hb MoM      ->      1.50 MoM <=> Hb 0.667 MoM
```

That is, **the cutoff is the arithmetic mirror image of the definition of
moderate anaemia (0.65 MoM) itself**, not a product of an ROC analysis.

Sensitivity and false-positive rate are purely properties of **measurement
error**, so the model can be asked what measurement error would reproduce
the reported figures. Using the Mari 2000 cohort composition
(non-anaemic 41, mild 35, moderate 4, severe 31):

| Measurement CV | Sensitivity | False-positive rate |
|---------|--------|----------|
| 10% | 98% | 5% |
| 15% | 97% | 10% |
| 20% | 94% | 14% |
| Reported | 100% | **12%** |

Sensitivity is reproduced at any CV (severe fetuses sit well above the
cutoff). The informative number is the **false-positive rate**, and it pins
Doppler reproducibility CV at 15-20% — a testable statement about the
examiner and the technique, not about the disease.

Adding a separate viscosity term breaks this result: raising `visc_k` from
0 to 0.45 to 0.90 to 1.80 shifts the Hb corresponding to 1.5 MoM from 0.667
to 0.703 to 0.733 to 0.780 MoM, and pushes the false-positive rate from 10
to 15 to 20 to 31%. Only `visc_k = 0` matches the observed 12%. The physical
reason is the same one: the increase in blood flow produced by lower
viscosity is itself the **mechanism** of oxygen-delivery defence, so
counting it separately is double counting.

### 3. The hydrops threshold is also a derived result

If oxygen delivery is defended and maximal cardiac-output reserve is
2.2-fold, then the point at which the fetus can no longer defend delivery
is `Hb MoM = 1/2.2 = 0.45`, i.e. roughly 5.5 g/dL. In the model, overt
ascites first appears at **0.49 MoM (6.28 g/dL, GA 29.9 weeks)** — consistent
with Nicolaides' clinical rule of roughly 5 g/dL or below. Umbilical venous
pressure rises with anaemia and then falls again at the most extreme values
(the non-monotonic pattern reported by Ville 1994), which also falls out of
the model unforced.

### 4. The 72-hour window for anti-D prophylaxis is not an immunological time constant

**An uncoated fetal red cell is just a red cell** — it survives in the
maternal circulation with a half-life of about 70 days. So a 3-day delay
loses only a few per cent of the antigen-exposure integral. That is why the
window is measured in days rather than minutes, and why it closes at around
2 weeks.

The dosing rule is **stoichiometry**, not pharmacokinetics: 20 µg of anti-D
(=100 IU/mL) per mL of fetal red cells. So 300 µg = 1500 IU covers 15 mL of
fetal red cells, i.e. 30 mL of fetal whole blood. The model's 300 µg curve
**bends exactly at 30 mL of fetal whole blood.**

| Regimen | Population sensitisation rate |
|------|-------------|
| No prophylaxis | 16.00% *(fitting target)* |
| 300 µg postpartum (within 24 h) | 1.60% *(fitting target)* |
| 100 µg postpartum (within 24 h) | 2.34% |
| 300 µg postpartum, 72 h delay | 2.17% |
| 300 µg postpartum, 7-day delay | 3.20% |
| 300 µg postpartum, 14-day delay | 4.75% |
| 28-week antenatal only | 2.28% |
| **28-week antenatal + 300 µg postpartum** | **0.26%** — prediction, observed 0.1-0.4% |
| 28-week antenatal + 600 µg postpartum | 0.22% |

This prediction holds because the residual risk is not an immunological
failure but an **arithmetic** one: the tail of the fetomaternal-haemorrhage
distribution that a fixed 300 µg dose cannot cover.

### 5. Simulating the UNITY trial as published

Rather than assuming Moise 2024's (NEJM) enrolment criterion, it was
**reproduced inside the model** — a virtual mother is enrolled only *if,
untreated, she would have needed a transfusion before 24 weeks*. A
log-normal prior is placed on anti-D titre, mothers are screened by that
criterion, and then dosed with nipocalimab 30 mg/kg/week.

| | Model | UNITY |
|---|------|-------|
| Live birth at 32+ weeks, no IUT | 11/30 = **37%** | 7/13 = **54%** (95% CI 25-81) |
| Hydrops | 1/30 | 0/13 |

The median anti-D among enrolled virtual mothers was 52 IU/mL (IQR 33-63),
and the screening pass rate was 38%. The prediction falls inside the
reported 95% confidence interval (25-81%). This result depends on **a
single distributional assumption** — the spread (GSD) of the prior — and its
sensitivity is reported alongside it in a table. The sharpest inference this
simulation produces concerns the drug itself:

> If placental syncytiotrophoblast FcRn saw the same nipocalimab
> concentration as vascular endothelial FcRn, transfer would be blocked by
> more than 99% and **no fetus in the UNITY cohort would have needed a
> transfusion.** Six of 13 did. So placental FcRn is far harder to saturate
> than endothelial FcRn (the model's assumption: 15% interstitial
> penetration).

### 6. Anti-K is not anti-D — the same anaemia, with the bilirubin missing

Anti-Kell destroys **red-cell progenitors**, not circulating red cells
(Vaughan 1998). In the model this is a single parameter, `kell_kill`, raised
from zero, and the consequence is that entirely different secondary markers
follow from the same haemoglobin. The clinical corollary is a warning the
model *produces*: **in Kell sensitisation, amniotic ΔOD450 underestimates
anaemia by 3.2-fold** (anti-D 0.0398 vs. anti-K 0.0124 under identical
conditions) — because the anaemia is created by cells that never reach the
circulation, so no haem is released. A surveillance strategy based on
bilirubin is structurally wrong in Kell disease; one based on MCA-PSV is
not.

### 7. Hydrops begins at 5-6 g/dL

In the untreated 15 IU/mL case, overt ascites first appears at **GA 29.9
weeks, Hb 6.28 g/dL (0.49 MoM)**. Albumin is suppressed by EMH, umbilical
venous pressure rises from 4.5 to 7.2 mmHg, and filtration increases from
17 to 71 mL/day while lymphatic return stalls at 17-59. A healthy fetus, by
definition, does not drift (lymphatic capacity = baseline filtration).

### 8. The neonatal period: low bilirubin, and a nadir six weeks later

A fetus that received four transfusions is born with a bilirubin **lower
than its disease severity would suggest**, because most of its circulating
haemoglobin is donor haemoglobin that nothing is attacking — the same f_ag
that flattened the in-utero decline rate also flattens the jaundice. And the
haemoglobin nadir arrives not in the first week but at **roughly day 39 of
life** (peak bilirubin 16.9 mg/dL, zero exchange transfusions, four top-up
transfusions, iron load 187 mg): maternal IgG is still present with its
3-week half-life, and a marrow suppressed by months of transfusion has no
reserve left. This is why these infants need top-up transfusions for 2-3
months.

---

## Defects Found by Running the Model (Not by Reading It)

This model exposed **eleven** genuine defects **through execution**, all
recorded in section 16 of `hdfn_reference_output.txt`. Because a defect that
produces plausible-looking numbers is the most dangerous kind, each entry
also records what it did.

1. Writing the circulation as "fixed blood volume + a decaying excess"
   made a transfusion look right immediately after the procedure but
   haemoconcentrated the fetus to 20 g/dL a day later. It was rewritten as
   **plasma volume + red-cell volume**, which exposed a second error the
   first had been hiding — packed red cells were diluting the fetal
   antibody, when in reality they do not.
2. Erythropoiesis responded only to anaemia and not to plethora, so a
   fetus made plethoric by transfusion kept producing red cells at maximum
   rate. This offset most of the post-transfusion decline and made the
   decline rate insensitive to the destruction rate. EPO was fixed to
   respond to Hb MoM **bidirectionally**.
3. Production was written as an absolute rate, which does not suit a
   fetus whose mass quadruples. It is now a **multiple of demand** (ageing
   plus growth in blood volume), so `k_prod = 1` is not a fitted number.
4. A baseline imbalance in the Starling balance caused a healthy fetus to
   accumulate ascites from 18 weeks onward, and the mortality-hazard
   function responded to it. Lymphatic capacity is now **defined as
   baseline filtration**.
5. With no long-lived plasma cells, every sensitised mother's titre
   collapsed during pregnancy (15 -> 0.06 IU/mL by 33 weeks).
6. Sensitisation was driven by the **clearance flux** of fetal cells,
   which makes a 72-hour window impossible. It is now driven by the
   exposure **integral**, from which the 72-hour window is a derived
   result.
7. Bilirubin was given the total-body-water volume of distribution (0.75
   L/kg). Because it is albumin-bound, it should be 0.20 L/kg. The
   arithmetic of physiological jaundice does not close with the wrong
   volume.
8. UGT1A1 maturation was written as a first-order process with a 6-day
   time constant, which reached 40% of adult activity by day 3 of life.
9. **Exchange transfusion added donor cells without removing any.** Donor
   mass grew with each procedure, so did its ageing, so did bilirubin
   production, which re-triggered the threshold — a positive feedback loop
   that reached total bilirubin of 1.2x10^7 mg/dL and 69 "exchanges" in one
   neonate. A double-volume exchange is volume-neutral and removes roughly
   85% of circulating content.
10. **The "fixed 2-week interval IUT" protocol never actually ran.** No
    code set the first trigger, so that arm silently reported zero
    transfusions and became numerically identical to no treatment.
11. **In the IVIG mechanism-decomposition analysis, the "FcγR-only" arm
    was implemented by zeroing the IVIG pool**, which removes both
    mechanisms at once. The arm looked ineffective and the narrative drew
    the opposite conclusion from the truth. The model now has an
    `ivig_compete` switch.

Also recorded, two silent bugs found not in the model but in the **analysis
scripts**: a bare percent sign inside a format string ("40% of adult") was
read as a `%o` conversion and killed the script in the final block of a
40-minute run, and the section-4 viscosity sensitivity was printing a
constant rather than actually solving for the shifted threshold.
`hdfn_section_rerun.py` exists because of this.

**Author expectations disproved by computation** (reported rather than
deleted):

- "Transfusion works by adding haemoglobin." -> It works at least as much
  by **removing the substrate**. Section 2 separates the two effects.
- "IVIG's job is done by FcγR blockade." -> **No.** After the switch was
  fixed, nearly all of the effect on delivered antibody turned out to be
  FcRn competition (5.57 -> 1.20 IU/mL), with FcγR blockade contributing
  almost nothing — because at this antibody burden the reticuloendothelial
  system is not the rate-limiting step.
- "An FcRn blocker that lowers maternal IgG by 70% should abolish the
  disease." -> It actually does, provided placental FcRn is not much
  harder to saturate than endothelial FcRn. UNITY's failure rate implies
  exactly that.
- "MCA-PSV needs a viscosity term." -> Adding one moves the 1.5 MoM
  threshold away from the definition of moderate anaemia and pushes the
  false-positive rate past the observed value (section 4).

---

## Files

| File | Content |
|------|------|
| [`hdfn_qsp_model.dot`](../../../hemolytic-disease-fetus-newborn/hdfn_qsp_model.dot) · [`.svg`](../../../hemolytic-disease-fetus-newborn/hdfn_qsp_model.svg) · [`.png`](../../../hemolytic-disease-fetus-newborn/hdfn_qsp_model.png) | Mechanistic map — **201 nodes / 18 clusters**. Boxes = state variables (ODEs), ovals = mediators/fluxes, hexagons = drugs/procedures, diamonds = clinical endpoints, octagons = **derived thresholds or arithmetic identities** |
| [`hdfn_mrgsolve_model.R`](../../../hemolytic-disease-fetus-newborn/hdfn_mrgsolve_model.R) | **45-ODE** mrgsolve model + an R driver that implements the obstetric protocols as a discrete-time controller + **16 scenarios** |
| [`hdfn_shiny_app.R`](../../../hemolytic-disease-fetus-newborn/hdfn_shiny_app.R) | **12-tab** dashboard. Tab 4 shows A x f_ag x M as separate trajectories, and every treatment control in the sidebar is labelled with which factor it moves |
| [`hdfn_references.md`](../../../hemolytic-disease-fetus-newborn/hdfn_references.md) | **96** references. Every one resolved live against NCBI esearch/esummary, with the *intent* — what the model took from that paper — recorded alongside each entry so a mismatch between the retrieved paper and the intent is visible at a glance |
| [`hdfn_python_reference.py`](../../../hemolytic-disease-fetus-newborn/hdfn_python_reference.py) | **The executed reference implementation.** Integrates the 45 ODEs with scipy. Where the mrgsolve file and this file disagree, this file is correct |
| [`hdfn_calibrate.py`](../../../hemolytic-disease-fetus-newborn/hdfn_calibrate.py) · [`calib.log`](../../../hemolytic-disease-fetus-newborn/calib.log) · [`hdfn_calibration.json`](../../../hemolytic-disease-fetus-newborn/hdfn_calibration.json) | The four-stage calibration and its full log |
| [`hdfn_analysis.py`](../../../hemolytic-disease-fetus-newborn/hdfn_analysis.py) · [`hdfn_reference_output.txt`](../../../hemolytic-disease-fetus-newborn/hdfn_reference_output.txt) | The script that produces every number above, and its output |
| [`mkrefs.py`](../../../hemolytic-disease-fetus-newborn/mkrefs.py) · [`refs_raw.json`](../../../hemolytic-disease-fetus-newborn/refs_raw.json) · [`refs_meta.json`](../../../hemolytic-disease-fetus-newborn/refs_meta.json) | The reference-resolution script and its raw payload — every citation can be re-derived |

## 45 Compartments (ODEs)

**Maternal (13)** Total IgG · anti-D IgG1 · anti-D IgG3 · memory B cells ·
plasma cells · uncoated fetal red cells · coated fetal red cells · passive
anti-D · IM depot · sensitisation signal · nipocalimab central/peripheral ·
IVIG-derived IgG

**Fetal antibody (4)** Total IgG · anti-D IgG1 · anti-D IgG3 · nipocalimab

**Circulation/red cells (8)** Estimated fetal weight · plasma volume ·
autologous antigen-positive Hb · donor antigen-negative Hb · reticulocytes ·
opsonised Hb · progenitors · extramedullary haematopoiesis

**Homeostasis (7)** EPO · reticuloendothelial capacity · albumin · ascites ·
cardiac decompensation index · lactate · nucleated red cells

**Bilirubin (7)** Plasma · tissue · enteric · photoisomers · amniotic fluid ·
UGT1A1 · cumulative neurotoxic exposure

**Records/risk (6)** Cumulative hazard function · iron load · postnatal
weight · cumulative antibody delivered to the fetus · cumulative destruction
· cumulative ageing

---

## Honest Limitations

- **The post-transfusion decline rate is under-predicted by roughly a
  third** (0.291 vs. 0.40 g/dL/day, SD 0.25). Even shortening donor
  red-cell lifespan to 30 days only reaches 0.359, so donor-cell survival
  alone does not explain it. The model keeps the adult value of 70 days
  and reports the gap.
- **It does not reproduce the UK "low-risk" band.** Even in the <4 IU/mL
  band, 44% still need a transfusion and 88% reach neonatal exchange
  transfusion in the virtual cohort (n=100). Because destruction is linear
  in titre while the conveyor is exponential in time, a 4-fold lower titre
  does not make the disease 4-fold milder — it **delays it** by ln(4)/0.268
  = 5.2 weeks. A threshold non-linearity at the antibody-to-phagocytosis
  step appears to be missing (section 14).
- MCA-PSV false-positive rate is under-predicted at a measurement CV of
  10% (5% vs. 12%). Reproducing 12% requires a CV of 15-20% (section 4).
- Maternal total IgG falls only about 12% during pregnancy, not the 30-40%
  Malek measured. Dilution and the placental sink are present, but **there
  is no term for reduced maternal IgG synthesis** — this is probably the
  missing term.
- The UNITY simulation is only as good as **one distributional
  assumption** — the spread of the anti-D titre prior. A sweep is reported
  alongside it showing how widely the answer moves when that assumption
  changes.
- Neonatal IVIG has no dosing route: the FcγR term is driven only by
  **maternal** IVIG.
- Amniotic ΔOD450 is modelled as a single mixed compartment. It does not
  quantitatively reproduce the Liley zones, only the anti-D vs. anti-K
  contrast.
- **Because the environment has no R toolchain**, `hdfn_mrgsolve_model.R`
  and `hdfn_shiny_app.R` mirror the executed Python reference
  implementation equation-by-equation but **have not themselves been run.**
- All parameters are literature-based approximations, and this model is
  for education and research purposes. It must not be used for real
  clinical decision-making.
