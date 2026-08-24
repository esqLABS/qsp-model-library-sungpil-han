# Congenital Hyperinsulinism (CHI)

> **The disease is written as a single number: the fraction `g` of
> K<sub>ATP</sub> conductance surviving in the beta cell.**
> And every drug is classified by **how it depends** on that one number —
> diazoxide **multiplies** `g`, octreotide **adds** to the voltage divider,
> and glucagon · glucose · ersodetug · surgery are **independent of** `g`.
> So the clinical fact that "recessive K<sub>ATP</sub>-CHI does not respond
> to diazoxide" is, in this model, **arithmetic, not an assumption**.

<p align="center">
  <a href="chi_qsp_model.svg"><img src="chi_qsp_model.png" width="900" alt="CHI QSP mechanistic map"></a>
</p>

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`chi_qsp_model.dot`](chi_qsp_model.dot) · [SVG](chi_qsp_model.svg) · [PNG](chi_qsp_model.png) | 150 nodes · 14 clusters · 199 edges |
| mrgsolve ODE model | [`chi_mrgsolve_model.R`](chi_mrgsolve_model.R) | **36 ODEs** · 10 genotypes · 10 scenarios |
| Shiny app | [`chi_shiny_app_en.R`](chi_shiny_app_en.R) | 10 tabs |
| References | [`chi_references_en.md`](chi_references_en.md) | 100 (99 PubMed links) |

---

## 1. Why a Single Number Suffices

The beta cell's resting membrane potential is a **conductance divider**.

```
              (G_KATP + G_GIRK)·E_K  +  g_leak·E_leak
      V_m  =  ────────────────────────────────────────
                  G_KATP + G_GIRK + g_leak

              G_KATP = g · g_max · P_open(ATP/ADP)
```

`g` enters at exactly one place. But because it is a divider, **three
classes of drugs fall out on their own.**

| Drug | Where it acts | `g` dependence | As `g→0` |
|---|---|---|---|
| **Diazoxide** | `P_open` (shifts the channel's ATP/ADP setpoint R₅₀ rightward) | **Multiplies** `g` | **Exactly 0** — there are no channels to open |
| **Octreotide** | SSTR2 → Gi → `G_GIRK` (a **second channel**, unaffected by mutation) | **Adds** to the divider | Still works (59 %) |
| **Glucagon / glucose** | Substrate supply | **Independent** | Works |
| **Ersodetug** | The insulin **receptor** (downstream) | **Independent** | Works |
| **Surgery** | Beta-cell **mass B** (a different symbol) | **Independent** | Works, but is irreversible |

Diazoxide and octreotide are both called "K-channel drugs," yet the same
equation produces opposite conclusions for them. This is the central claim
of this model.

### Verified Results

| `g` | Untreated IV glucose | **Diazoxide 15 mg/kg/d** | **Octreotide 30 µg/kg/d** |
|---:|---:|---:|---:|
| 0.30 | 6.8 mg/kg/min | **100 %** reduction | 100 % reduction |
| 0.20 | 10.1 | **49 %** | 100 % |
| 0.10 | 11.1 | **4.9 %** | 100 % |
| 0.05 | 11.2 | **1.1 %** | 71 % |
| 0.02 | 11.2 | **0.6 %** | **59 %** |

It is clearer as a dose titration (mean glucose, mg/dL, 0 → 15 mg/kg/day):

| Genotype | 0 | 5 | 10 | 15 |
|---|---:|---:|---:|---:|
| Normal beta cell | 99.8 | 122.2 | 134.5 | 142.1 |
| K<sub>ATP</sub> dominant (`g`=0.60) | 80.4 | 93.5 | **100.3** | 104.5 |
| GDH-HI (GLUD1) | 76.3 | 91.9 | **100.2** | 105.5 |
| K<sub>ATP</sub> recessive (`g`=0.02) | 31.8 | 31.9 | 31.9 | **31.9** |

The last row is **completely flat.** This reproduces, **without
calibration**, the clinically observed split that "biallelic recessive
K<sub>ATP</sub>-CHI does not respond to diazoxide, while GDH-HI/HNF4A/
dominant forms almost always do."

---

## 2. What Was Calibrated, and So What Is a Prediction

**Nine** numbers were used, and **eight of them are normal neonatal
physiology**.

| | What was calibrated | Target | Achieved |
|---|---|---|---|
| N1 | Basal insulin at a glucose of 75 mg/dL | 5 µU/mL | 5.0 |
| N2–N5 | The **shape** of the normal glucose→insulin dose-response (secretion ratio vs 75 mg/dL) | 0.25× at 45 · 1.8× at 90 · 12× at 150 · 24× at 250 | 0.254 · 2.04 · 12.76 · 23.25 |
| N6 | Neonatal brain glucose consumption / BBB GLUT1 K<sub>m</sub> / brain requirement | 4.2 mg/kg/min · 40 mg/dL · 4.0 | As given |
| N7 | Term-infant hepatic glycogen | ≈2000 mg/kg | 1885 |
| N8 | Insulin IC₅₀ **sequence** | Lipolysis 12 < ketogenesis 15 < glycogenolysis 30 < gluconeogenesis 45 µU/mL | As given |
| **C1** | **The single number taken from CHI itself**: `gGIRK`=1.2 — octreotide halves the glucose requirement of severe diffuse CHI | ≈50 % | 64 % |

Diazoxide's potency (`kdzx`=0.60) was also fitted to the **normal** beta
cell: the condition being that 10 mg/kg/day should not abolish secretion
but should produce mild iatrogenic hyperglycaemia (134 mg/dL).

### So Everything Below Is a Prediction

| | Prediction | Model | Literature/clinical |
|---|---|---|---|
| **P1** | Total glucose supply in severe recessive diffuse CHI | IV 11.2 + enteral 4.7 = **15.9 mg/kg/min** | 15–20 mg/kg/min |
| **P2** | Genotype split in diazoxide response | **0.1 %** at `g`=0.02, 100 % at `g`≥0.3 | Recessive non-response / response otherwise |
| **P3** | Octreotide retains effect even at `g`=0 | **59 %** (vs diazoxide's 0.6 % at the same point) | Partial effect, cannot substitute for surgery |
| **P4** | **Severity saturation** | `g`=0.10 → 11.07, `g`=0.02 → 11.21 mg/kg/min | No severity gradient between null and severe hypomorphs |
| **P5** | Loss of the brain-fuel safety margin | **43.7 mg/dL** at zero ketones, **25.3 mg/dL** at BOHB 2 mM → **an 18.4 mg/dL loss** | CHI target of 70, vs ketotic hypoglycaemia tolerated into the 40s |
| **P6** | Glucagon stimulation test | CHI **+62 mg/dL** (glycogen 2599) vs normal 20 h fast **+16** (glycogen 2) | A rise >30 mg/dL is diagnostic |
| **P7** | A **narrow window** for resection extent | Residual 0.5 = still glucose-dependent, 0.3–0.2 = euglycaemic, ≤0.1 = **diabetic** | The real dilemma of persistent hypoglycaemia vs surgical diabetes |
| **P8** | GDH-HI signature | **−37** after a leucine load vs dominant K<sub>ATP</sub>'s −16 mg/dL; ammonia at 150 µmol/L, **150.0 → 150.0** with diazoxide | Protein sensitivity + diazoxide-unresponsive hyperammonaemia |
| **P9** | **Why** nifedipine fails | 5.8 % — the mechanism is right, but the beta-cell EC₅₀ of ≈1200 ng/mL is about 8-fold the plasma Cmax | No clinical benefit reproduced |
| **P10** | The quantitative requirement for a focal lesion | (lesion fraction × local density) must be ≈ **1.0** to match diffuse-disease severity. Severity depends **only on that product** (0.10×5 = 4.23, 0.05×10 = 4.29) | The high-density, enlarged beta-cell histology of focal adenomatous hyperplasia |

---

## 3. The Second Axis: Insulin Deletes the Backup Fuel Before It Takes the Glucose

This is the part of the disease most often described but almost never
quantified.

Because of the insulin IC₅₀ sequence (N8), **any insulin concentration
high enough to cause hypoglycaemia has already switched off ketogenesis.**
There is also a second, independent reason: the liver's fasting switch
**arms only once glycogen is depleted**, and in CHI, insulin protects
glycogen, so the switch never arms in the first place.

So writing brain fuel as a **sum** — glucose + ketones + lactate against
requirement — the CHI-specific glucose target falls out arithmetically.

| Plasma BOHB | Glucose at which brain fuel is exhausted | Fraction supplied by ketones |
|---:|---:|---:|
| 0.0 mM (CHI) | **43.7 mg/dL** | 0 % |
| 0.5 | 36.1 | 5.8 % |
| 1.0 | 31.3 | 10.7 % |
| 2.0 (ketotic hypoglycaemia) | **25.3 mg/dL** | 18.8 % |
| 4.0 | 19.0 | 30.0 % |

**Hyperinsulinism deletes an 18.4 mg/dL safety margin before glucose has
even fallen by 1 mg/dL.** Requiring a CHI child to stay at 70 mg/dL while
a child with ketotic hypoglycaemia tolerates the 40s is not arbitrary
conservatism — it is this subtraction.

It was decisive here that the brain MCT's K<sub>m</sub> (≈6 mM) is the
**rate-limiting step**. Setting K<sub>m</sub> to 1.5 mM calculates a
tolerable glucose of 16 mg/dL for ketotic hypoglycaemia, which disagrees
with the clinic, and that discrepancy **forced** the choice of 6 mM — a
structural constraint, not a calibration.

---

## 4. Surgery: Reduces `B`, Leaves `g` Alone

Surgery does not touch `g`. It changes only the mass `B`. So **a single
equation** produces both failure modes of CHI surgery.

| Residual beta-cell mass | IV glucose | Mean glucose | Interpretation |
|---:|---:|---:|---|
| 1.00 | 11.2 | 69.9 | Pre-surgery |
| 0.50 | 3.4 | 69.9 | **Still glucose-dependent** (insufficient resection) |
| 0.30 | 0.0 | 77.9 | Euglycaemic |
| 0.20 | 0.0 | 97.0 | Euglycaemic |
| 0.10 | 0.0 | **176.9** | **Diabetic** |
| 0.02 | 0.0 | **341.5** | **Diabetic** |

The window is only a residual mass of 0.20–0.30. And **growth re-reads the
same residual mass**: a 2% residual that was adequate at 3.5 kg becomes
0.70% of normal at 10 kg, 0.35% at 20 kg, and 0.18% at 40 kg. Even if
nothing new happens to the pancreas, the fact that about half of children
who undergo subtotal pancreatectomy become diabetic by adolescence is
**because the denominator grows**.

Focal disease is different. Lesionectomy sets `w_ab = 0`, leaving a
normal pancreas with `g`=1 — a cure.

---

## 5. Model Structure (36 ODEs)

| System | State variables |
|---|---|
| Systemic glucose | `GLU` blood glucose · `GLUi` interstitial fluid (CGM, 10-minute lag) |
| Insulin axis | `INS` · `X` receptor-bound insulin · `CPEP` C-peptide |
| Beta cell | `CABn`/`CABa` Ca signalling in the normal/abnormal populations · `CAMP` · `BMASS` mass |
| Liver | `GLY` glycogen (**a state variable** — this is where the glucagon test's discriminating power comes from) |
| Fat · ketones | `FFA` · `BOHB` |
| Counter-regulation | `GCG` · `EPI` · `CORT` |
| Other metabolism | `LAC` · `NH3` · `AA` leucine · `GGUT` gut glucose |
| Pharmacokinetics | Diazoxide 3 compartments · octreotide 2 + `ROCT` desensitisation · glucagon · ersodetug 2 compartments · sirolimus 2 · nifedipine · exendin(9-39) |
| Outcome integrals | `AUCHYPO` · `TFUEL` cumulative time of brain-fuel deficit · `DEV` neurodevelopmental deficit |
| Controller | `GIRi` — the integral term auto-titrating IV glucose (**the glucose requirement becomes the output**) |

It matters that the glucose requirement is not a parameter but **the
output of a closed-loop controller** — because that number is exactly the
metric used clinically to gauge CHI severity.

---

## 6. Verification — and the Defects It Caught

All 36 ODEs were **independently reimplemented in pure Python RK4** for
cross-validation. In the process, **five real defects** surfaced and were
fixed, and that is the only reason the numbers above can be trusted.

1. **Basal amino-acid drive was too large** (`kAA`=6). Normal basal
   insulin came out at 19–22 µU/mL instead of 5. Fixed by a grid search
   refitting (`R50`, `nR`, `gKmax`, `kAA`) to the *shape* of the normal
   dose-response.
2. **Glycogen synthesis ran on basal insulin alone.** As a result,
   "fasting" was not really fasting, glycogen never depleted, and both the
   normal neonate's fasting tolerance and ketone rise disappeared. Fixed
   by switching synthesis to be driven by **portal glucose influx (Ra)**
   — without this fix, P6 does not hold.
3. **Unit mismatch for nifedipine and sirolimus.** Concentrations were in
   mg/L while EC₅₀ was written in ng/mL, a 1000-fold discrepancy that made
   both drugs appear to have no effect at all.
4. **The glucagon dose was 1000-fold too small** (a missing µg/kg →
   pg/mL conversion). Because of this, the glucagon stimulation test came
   out at +1.7 mg/dL. There was also no term for glucagon/PKA **pushing
   insulin's glycogenolysis IC₅₀ rightward** (PKA vs PP1 competition), so
   no pharmacological dose of glucagon could overcome insulin's block.
5. **Diazoxide was written as "forcing a fixed fraction of channels
   open."** As a result, 2.5 mg/kg/day abolished secretion in the normal
   beta cell, raising glucose to 260 mg/dL. Rewriting it, per the
   literature (Shyng 1997), as **an MgADP-dependent shift of the
   channel's ATP/ADP setpoint** produced mild hyperglycaemia in the
   normal beta cell, and **simultaneously** preserved the structure of
   multiplying `g`, so P2 survived.

Verification also revealed that **the closed loop quietly cancels out the
rise in the glucagon test.** The actual test cannot be reproduced unless
the glucose infusion is held fixed — a trap in the experimental design,
recorded in the file.

---

## 7. Tensions Stated Honestly (nothing hidden)

- **T1 — the most exposed number in this model.** The predicted plasma
  insulin for untreated severe diffuse CHI is **≈83 µU/mL**, higher than
  reported critical-sample values (10–50 µU/mL). This is not a free
  choice but **forced arithmetic**: there is no way to explain a glucose
  requirement of 16 mg/kg/min at neonatal insulin sensitivity with
  20 µU/mL. One of three things is true — (i) a single critical sample
  underestimates 24-hour exposure, (ii) neonatal insulin sensitivity is
  higher than the model assumes, or (iii) the model's `Uidmax`/`KI` are
  wrong. **This is a falsifiable point.**
- **T2** — normal fasting glucose stays at ≈72 mg/dL at 18–24 h (higher
  than the measured 55–65). This is because gluconeogenesis has no
  substrate ceiling. It does not affect the CHI branch (which is
  insulin-driven), and it makes the model **conservative** for normal
  fasting hypoglycaemia.
- **T3** — 75% at ersodetug 9 mg/kg is larger than the published early
  clinical effect size. Because the allosteric insulin-receptor
  antibody's E<sub>max</sub>/EC₅₀ have not been publicly established, this
  branch should be read only **structurally** (as evidence that a
  `g`-independent mechanism exists), not quantitatively.
- **T4** — sirolimus (19%) was written as a mass/secretion effect because
  its mechanism is genuinely unresolved. The number carries no more
  weight than that assumption.
- **T5** — the fact that hypothalamic K<sub>ATP</sub> carries the same
  mutation (impairing hypoglycaemia awareness itself) is included on the
  map only structurally and **not quantified**. Counter-regulatory
  failure appears only as a consequence of α-cells being suppressed by
  local insulin.
- **T6** — setting the GCK-activating form to `KGshift`=0.45 produces a
  fairly severe phenotype (61.8 mg/dL untreated). Because the actual
  severity of GCK-HI varies greatly by mutation, this value is only a
  single illustrative example.

---

## 8. Usage

```bash
# the map
dot -Tsvg chi_qsp_model.dot -o chi_qsp_model.svg
dot -Tpng -Gdpi=150 chi_qsp_model.dot -o chi_qsp_model.png
```

```r
# model + 10 scenarios
source("chi_mrgsolve_model.R")
print(s2)   # the disease axis: a sweep over g
print(s3)   # g-dependence of diazoxide vs octreotide  <- the central result
print(s7)   # resection extent

# interactive explorer (10 tabs)
shiny::runApp("chi_shiny_app_en.R")
```

Required packages: `mrgsolve`, `dplyr`, `tidyr`, `ggplot2`, `shiny`, `DT`.

---

## 9. Clinical Background Summary

Congenital hyperinsulinism is **the most common cause of persistent
neonatal hypoglycaemia** (about 1/40,000–1/50,000 live births, up to
1/2,500 in consanguineous populations), and carries an unusually high
risk of neurological damage because inappropriate insulin secretion
causes hypoglycaemia while **simultaneously blocking the alternative fuel
(ketones)** — neurodevelopmental abnormalities are reported in 25–50%
despite treatment.

Diagnosis relies on a critical sample at the time of hypoglycaemia
(detectable insulin, elevated C-peptide, suppressed BOHB and free fatty
acids), the glucagon stimulation test, and distinguishing focal from
diffuse disease with genetic testing and ¹⁸F-DOPA PET. The treatment
ladder runs glucose → diazoxide → octreotide/lanreotide → glucagon →
(focal) lesionectomy or (diffuse, refractory) subtotal pancreatectomy,
with ersodetug, which targets the insulin receptor, being developed as a
new axis.

This model's contribution is to **rearrange that ladder as dependence on
a single parameter.** Detailed justification is in
[`chi_references_en.md`](chi_references_en.md), which states, reference by
reference, "which term of this model depends on that reference."

---

> ⚠️ **Disclaimer** — This is a qualitative/semi-quantitative QSP model for
> educational and research purposes. It was built from public literature
> but has not been independently validated or certified, and **must not be
> used for actual clinical decision-making, prescribing, or regulatory
> submission.** The parameters are illustrative approximations, and the
> tensions in §7 (T1–T6) remain unresolved.
