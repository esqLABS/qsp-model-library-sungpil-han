# Intrahepatic Cholestasis of Pregnancy (ICP) — QSP Model

> **One-line summary.** This model does not treat ICP as a question of "what is
> the total bile acid". It **resolves the bile acids into five species**, moves
> them through six maternal and three fetal compartments, and only at the end
> **reconstructs, via an observation equation**, the laboratory values used in the
> clinic. Doing it that way turns the two paradoxes of this field from assertions
> into arithmetic — **why stillbirth risk suddenly turns upward at
> 100 µmol/L**, and **why a drug that lowers total bile acids by a third failed to
> change perinatal outcomes in a randomised trial**.

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map (191 nodes · 182 edges · 20 clusters) | [`icp_qsp_model.dot`](../../../intrahepatic-cholestasis-of-pregnancy/icp_qsp_model.dot) · [`.svg`](../../../intrahepatic-cholestasis-of-pregnancy/icp_qsp_model.svg) · [`.png`](../../../intrahepatic-cholestasis-of-pregnancy/icp_qsp_model.png) |
| ⚙️ mrgsolve ODE model (78 ODEs, 35 scenarios) | [`icp_mrgsolve_model.R`](icp_mrgsolve_model.R) |
| 📊 Shiny dashboard (10 tabs) | [`icp_shiny_app.R`](icp_shiny_app.R) |
| 📚 References (202 items, every PMID looked up) | [`icp_references.md`](icp_references.md) |
| 🐍 Python RK4 reference implementation | [`icp_reference_model.py`](icp_reference_model.py) · [output](icp_reference_output.txt) |
| 🔬 Fitting, ablation and trial-reproduction analysis | [`icp_calibration.py`](../../../intrahepatic-cholestasis-of-pregnancy/icp_calibration.py) · [output](../../../intrahepatic-cholestasis-of-pregnancy/icp_calibration_output.txt) |

---

## 1. Why it has to be computed species by species

In the clinic the severity of ICP is stratified by a single number, the **total
bile acid (TBA)** (≥10 diagnostic · ≥40 severe · ≥100 most severe). But that
number is measured by the 3α-hydroxysteroid dehydrogenase enzymatic assay, so it
**counts the administered UDCA and its metabolite hyocholic acid along with
everything else.** And what actually harms the fetus is not the total but the
hydrophobic species (LCA · DCA · CDCA).

This model therefore carries five species separately: CA · CDCA · DCA · LCA ·
UDCA.

```
maternal:  hepatocyte → bile/gallbladder → small intestine → colon → plasma   (+ renal excretion, 6α-hydroxylation sink)
fetal:     plasma/extracellular fluid ↔ fetal liver → fetal bile → meconium,  + amniotic fluid
link:      placenta  ── saturable active efflux J_F2M  +  bidirectional diffusion J_PASS
```

The three different "bile acids" that are computed:

| Variable | Meaning | Normal pregnancy | Severe ICP | Most severe ICP |
|------|-----|----------|----------|-----------|
| `TBA` | the clinical laboratory value (includes UDCA and HCA) | 5.3 | 58.2 | 137.7 |
| `TBA_ENDO` | the endogenous fraction | 4.1 | 55.1 | 133.7 |
| `FCL` | fetal hydrophobic load Σ wᵢ·C_fetal | 1.15 | 6.31 | 17.20 |

The genetic predisposition vectors (`GBSEP`, `GMDR3`, `GSULT`) of the three
strata were **not specified in advance**; they were found in exactly the way the
literature defines the strata — that is, by the TBA band each vector generates.

---

## 2. The timing of onset and the postpartum resolution were not fitted

The lesion consists of two inhibition terms. Estradiol-17β-glucuronide inhibits
BSEP competitively (cis), and progesterone sulfate exerts trans-inhibition and
FXR antagonism. The precursors of both metabolites are made by the placenta, and
placental production rises exponentially until term and **becomes 0 at
delivery.**

The following therefore falls out as it stands, with no change of parameters:

| Scenario | Crossing of the diagnostic threshold (10) | Peak TBA | Return to <10 after delivery |
|----------|-------------------|---------|-----------------|
| Normal pregnancy | never crossed | 5.3 | 1 day |
| ICP, <40 band | **34.4 weeks** | 15.8 | 2 days |
| ICP, 40–99 band | 22.1 weeks | 58.2 | 4 days |
| ICP, ≥100 band | 20.0 weeks | 137.7 | 7 days |
| Twins, <40 genotype | **26.7 weeks** | 43.6 | 3 days |

The mild genotype first crosses the threshold at 34 weeks (third-trimester
onset) and normalises 2 days after delivery, and **multiplying that same genotype
by the single twin factor (sex steroid load ×1.55) brings the onset forward by
8 weeks and raises the peak 2.8-fold**. This is a mechanistic reading of why ICP
is far commoner in twin pregnancy, and it points the same way as the case reports
in which selectively reducing one twin reverses the disease.

---

## 3. The 100 µmol/L threshold sits in the fetal myocardium, not in the placenta

This is the central result of the model, and **the place where the model
corrected the author's expectation**.

The map was drawn on the hypothesis that "saturation of the placental transporter
creates the threshold". So, against the same stillbirth data (Ovadia 2019:
0.13% / 0.28% / 3.44% at TBA <40 / 40–99 / ≥100), the fit was performed **three
times, changing nothing but which variable the hazard function is written in**.
The smaller the exponent required, the more of the non-linearity that variable
already carries.

| Variable the hazard function was written in | <40 | 40–99 | ≥100 | Fold | **Required exponent** |
|---|---|---|---|---|---|
| maternal total bile acid (the laboratory value) | 15.8 | 58.2 | 137.7 | 8.7× | **2.95** |
| fetal hydrophobic load FCL | 1.94 | 6.31 | 17.20 | 8.9× | **2.55** |
| myocardial arrhythmia index ARRI | 0.011 | 0.082 | 0.414 | **36.8×** | **1.55** |

Even after passing all the way through placental transport, the fold change moves
barely at all, from 8.7 to 8.9. **Transport alone does not account for the
threshold.** Yet across the same three pregnancies the arrhythmia index spreads
36.8-fold, and then an exponent of 1.55 suffices.

The ablation experiments confirm this (section I of `icp_calibration.py`):

| Ablation | ARRI fold (≥100 / <40) | Required exponent |
|------|----------------------|----------|
| (reference) | 36.8 | 1.55 |
| placental transporter desaturated (V_P ×20) | 45.9 | 1.55 |
| species-specific diffusive permeabilities made uniform | 40.1 | 1.50 |
| cytotoxicity weights made uniform | 39.2 | 1.50 |
| connexin-43 uncoupling reduced **in magnitude only** | 36.8 | 1.55 |
| connexin-43 uncoupling **linearised** | **12.7** | **2.15** |
| calcium overload term removed | 24.0 | 1.80 |

What destroys the threshold is the **threshold shape** of the connexin-43
uncoupling. Reducing only its magnitude leaves the fold change unchanged
(multiplying a Hill function by a constant leaves the shape the same), whereas
making it linear collapses it. That is a sharper claim, and a claim it is harder
to satisfy by accident.

**A clinically testable prediction:** risk should track the **hydrophobic
fraction of the fetal pool**, not the maternal total bile acid. Species-resolved
quantification of bile acids in cord blood at delivery tests this directly, and
the maternal laboratory value cannot do so even in principle.

As a by-product, **why it is acute and why it is not caught by fetal movement
monitoring** also falls out: the risk is carried not by a "progressive
deterioration in fetal wellbeing" but by a step change in intercellular coupling
in the myocardium.

Meanwhile, saturation of placental transport does have a job of its own — it sets
the **absolute level of fetal exposure**, and it makes the cord:maternal ratio
fall as far as 0.30 in mild disease and then rise back above 0.40 in the ≥100
band (the point at which transporter occupancy exceeds 1).

---

## 4. The arithmetic by which PITCHES was negative

In the PITCHES randomised trial, UDCA lowered total bile acids yet failed to
change the composite perinatal endpoint (24.8% vs 27.8%, not significant).
Running the model against the trial's population composition:

| Endpoint | Placebo | UDCA | Relative change | Value reported in the trial |
|-----------|------|------|----------|-----------|
| total bile acids (µmol/L) | 35.6 | 23.3 | **−34.5%** | decreased |
| fetal hydrophobic load FCL | 3.87 | 2.40 | −38.0% | not measurable |
| ALT (U/L) | 124.8 | 76.0 | −39.1% | decreased |
| pruritus VAS (0–10) | 4.61 | 4.27 | **−0.34** | **−0.7 cm** |
| stillbirth % | 0.20 | 0.20 | −0.8% | 2 cases vs 1 case |
| spontaneous preterm birth <37 weeks % | 9.75 | 8.77 | −10.1% | — |
| neonatal intensive care unit ≥4 h % | 11.63 | 10.45 | −10.1% | — |
| **composite %** | **20.39** | **18.45** | **−9.5%** | **27.8 → 24.8 (−10.8%)** |

Of the three components of the composite endpoint, **two are not functions of
bile acids but functions of gestational age and of the decision to deliver**. The
model is built that way too (delivery before 37 weeks is a delivery decision, and
the probability of neonatal intensive care unit admission is a function of GA).
Stillbirth accounts for only **0.97%** of the composite in the placebo arm, so
abolishing stillbirth entirely moves the composite by **0.20 percentage points**.

Converted into power:

| Sample size per arm required for 80% power on the component |
|---|
| composite endpoint — **6,262** |
| stillbirth alone — **more than 10 million** |
| neonatal intensive care unit alone — 12,993 |

With 604 patients it **could not have detected its own mechanism.** This result is
therefore not in contradiction with the stillbirth signal reported by the
individual patient data meta-analysis — the two studies were measuring different
things.

---

## 5. The bile acid axis and the pruritus axis respond to different drugs

The result of adding each prescription to a pregnancy in the 40–99 band (section
F of `icp_calibration.py`):

| Prescription | TBA | ΔTBA | FCL | ΔFCL | Autotaxin | VAS | ΔVAS | ALT | INR |
|------|-----|------|-----|------|---------|-----|------|-----|-----|
| none | 58.2 | — | 6.31 | — | 3.50 | 5.6 | — | 167 | 1.10 |
| UDCA 1000 mg/day | 35.9 | −38% | 3.51 | −44% | 3.47 | 5.2 | **−0.4** | 103 | 1.00 |
| UDCA 1500 mg/day | 31.3 | −46% | 2.96 | −53% | 3.46 | 5.1 | −0.5 | 88 | 1.00 |
| **rifampicin 600 mg/day** | 44.6 | **−23%** | 3.48 | −45% | **1.70** | **1.0** | **−4.6** | 71 | **1.28** |
| UDCA + rifampicin | 45.7 | −21% | 3.25 | −48% | 1.70 | 1.0 | −4.6 | 60 | 1.00 |
| colestyramine 16 g/day | 32.5 | −44% | 3.32 | −47% | 3.46 | 5.1 | −0.4 | 119 | **1.13** |
| SAMe 1000 mg/day IV | 33.0 | −43% | 3.51 | −44% | 3.47 | 5.2 | −0.4 | 106 | 1.06 |
| IBAT inhibitor (hypothetical) | **28.9** | **−50%** | 3.09 | −51% | 3.46 | 5.1 | −0.5 | 116 | 1.05 |
| naltrexone 50 mg/day | 58.2 | +0% | 6.31 | +0% | 3.50 | 3.7 | −1.9 | 167 | 1.10 |
| antihistamine 12 mg/day | 58.2 | +0% | 6.31 | +0% | 3.50 | 4.3 | −1.3 | 167 | 1.10 |

**Read the ΔFCL column and the ΔVAS column together.** The drugs that bring the
bile acid pool down (UDCA · colestyramine · IBAT inhibitor) all leave the
pruritus very nearly where it was, and the drug that abolishes the pruritus
(rifampicin) is only middling on bile acids. This rule was not put into the
model — it comes out of a structure in which autotaxin is **driven by
progesterone sulfate and estradiol, with the bile acid term small and
saturating**. That structure was chosen so as to satisfy four observations at
once:

1. circulating autotaxin does not fall with UDCA, but does fall with rifampicin
   and with nasobiliary drainage
2. the degree of pruritus correlates poorly with total bile acids
3. in a substantial proportion of patients the pruritus precedes the biochemical
   abnormality by several weeks
4. asymptomatic hypercholanaemia exists

A model that drives pruritus from bile acids reproduces none of these.

**The practical implication:** a trial with pruritus as its primary endpoint
cannot rank bile-acid-lowering agents, and a trial with bile acids as its primary
endpoint cannot rank antipruritic agents. That is roughly the state of the ICP
literature at present. Since a trial directly comparing UDCA with rifampicin
(TURRIFIC) has been conducted, this prediction of the model — that rifampicin
leads by a wide margin on pruritus and trails on bile acids — is directly
testable.

The `INR` column also deserves a look. Colestyramine and rifampicin lower the
intraluminal bile acid concentration (by binding and by pool contraction
respectively) and so reduce micelle formation, meaning they **trade vitamin K for
bile acids**. UDCA and the IBAT inhibitor do not pay that price. The model
computes the exchange rate.

### Two incidental results (things that simply fell out of the 35-scenario table)

**The ABCB4 phenotype separates itself out.** Scenario 23, with MDR3 lowered to
30%, has a term total bile acid of **14.3 µmol/L** (very nearly normal) yet an
ALT of **272 U/L**, the highest of all the scenarios. It does not cross the
diagnostic threshold until 35.7 weeks. No such rule was put into the model; it
comes out of a structure in which the biliary bile acid:phospholipid ratio drives
the injury — MDR3 deficiency is not a disease that raises bile acids but a
disease that **raises the canalicular detergency ratio**. This points the same
way as the clinical picture in patients with ABCB4 variants, in whom GGT is
elevated and bile acids are comparatively low.

**Delivery at 36 weeks fails the PITCHES composite endpoint by definition.** The
composite endpoint of scenarios 27, 28 and 35 (most severe, delivery at 36 weeks)
is **100%**. Since "delivery before 37 weeks" is part of the definition of the
composite, following exactly the management the guidelines recommend at
≥100 µmol/L makes that pregnancy meet the endpoint automatically. In other words,
this composite endpoint **counts the recommended treatment as a failure.** It is
the other face of the same story as the arithmetic in section 4 about why PITCHES
could not give a useful answer.

---

## 6. Timing of delivery — where it agrees with the guidelines and where it does not

The delivery decision is written not as a threshold but as **the minimisation of
an integral**: expected loss = cumulative stillbirth risk + `WMORB` × (neonatal
intensive care unit + respiratory distress). `WMORB` (the utility ratio of
neonatal morbidity to stillbirth) is the one judgement value in this model, and it
is reported rather than hidden.

At `WMORB` = 0.055 (milli-expected-loss; only the location of the minimum is
meaningful):

| Stratum | 35 weeks | 36 | 37 | 38 | 39 | 40 | **Optimum** |
|------|------|----|----|----|----|----|---------|
| TBA <40 | 42.94 | 26.31 | 15.62 | 9.79 | 7.09 | 6.08 | **40 weeks** |
| TBA ~26 | 43.14 | 26.58 | 15.96 | 10.22 | 7.64 | 6.76 | **40 weeks** |
| TBA 40–99 | 44.16 | 27.90 | 17.69 | 12.46 | 10.51 | 10.43 | **40 weeks** |
| TBA ≥100 | 59.24 | 47.01 | **41.58** | 41.87 | 46.13 | 52.94 | **37 weeks** |
| TBA ≥100 + UDCA | 51.97 | 36.79 | 27.85 | 24.11 | **23.86** | 25.71 | **39 weeks** |

For the **≥100 band** the answer is 37 weeks, pointing in the same direction as
the guidelines (36–37 weeks). And with treatment the optimum **is pushed out to
39 weeks** — in this model the value of UDCA lies not in a biochemical endpoint
but in **about two weeks of gestation**. That is a testable claim, and as far as
we know an untested one.

**The 40–99 band does not reproduce the guidelines.** The model computes 39–40
weeks where the guidelines recommend 37–38. The reason is arithmetic: the excess
stillbirth risk in that band is 0.28% against 0.13%, i.e. **0.15 percentage
points**, and no reasonable weighting makes it outweigh the neonatal cost of
delivering 2–3 weeks early. For 37 weeks to be optimal, `WMORB` would have to be
< 0.01 (neonatal intensive care unit admission worth less than one hundredth of a
stillbirth), and at that weighting the <40 band gets pulled forward along with
it, giving a conclusion no guideline recommends.

**Do not read this as a victory for the model.** Two interpretations are both
available and the model alone cannot adjudicate between them. (a) The basis for
the recommendation does not lie in the stillbirth numbers. (b) The model has
missed something real in that band. The strongest candidate on the (b) side is
**within-individual variability** — someone measured at 70 µmol/L may spend days
above 100 between visits, and this model runs smooth trajectories. A version
driven by serial measurements is the way to test this, and it is the most useful
direction of extension we have identified.

---

## 7. Mechanistic dissection of UDCA — and the lithocholic acid problem

Switching off each mechanism of UDCA one at a time (cells = maternal TBA / fetal
hydrophobic load):

| Ablation | 500 mg/day | 1000 mg/day | 2000 mg/day |
|------|-----------|------------|------------|
| reference | 43.5 / 4.43 | 35.9 / 3.51 | 28.4 / 2.60 |
| effect of BSEP insertion into the plasma membrane removed | 48.0 / 4.92 | 40.9 / 4.01 | 32.8 / 2.98 |
| intestinal 7β-dehydroxylation (→LCA) removed | 42.3 / 4.37 | 33.7 / 3.42 | 24.5 / 2.47 |
| hepatic LCA sulfation removed | 51.9 / **6.22** | 48.8 / **6.21** | 47.7 / **6.71** |

How to read it: **BSEP insertion into the plasma membrane is the therapeutic
mechanism** (switch it off and a substantial part of the effect disappears). The
route by which the intestinal flora convert UDCA into lithocholic acid
(cytotoxicity weight 1.00, the highest of the five species) is, over this dose
range, **a real but gentle burden** — switching it off improves the fetal load
from 3.51 to 3.42 at 1000 mg and from 2.60 to 2.47 at 2000 mg, and the gap widens
with dose. That is because the LCA-generating route is first-order in dose while
the hepatic sulfation that handles it saturates.

The last row is the most important. Switching off hepatic LCA sulfation
**inverts the dose-response** — the fetal hydrophobic load goes 6.22 → 6.21 →
6.71, so that above 1000 mg escalating the dose stops helping and starts harming.
That is to say, the entire dose safety of UDCA hangs on the spare capacity of
SULT2A1. This is the model's reading of why "UDCA is well tolerated in ordinary
ICP and becomes a genuine worry in the most severe cholestasis, where that enzyme
system is already fighting a far larger substrate load", and it is a testable
prediction: when high-dose UDCA is used in the most severe ICP, the **lithocholic
acid fraction** in the circulation and in cord blood is what must be watched, and
the total bile acid does not show this risk.

> This section is the point at which the model actually told the author about a
> problem during development. In the version that lacked the LCA detoxification
> route, UDCA made the fetal pool 43% more hydrophobic, so **the drug came out
> harmful**, and that was the signal that the SULT2A1 sulfation route had to be
> put in.

---

## 8. Nine actual defects that surfaced during development

With no R runtime available, every equation was first run in dependency-free
Python RK4 (`icp_reference_model.py`), and in the course of that nine actual
defects came to light. Each is marked `[DEFECT n]` in the mrgsolve file.

| # | Symptom | Cause and fix |
|---|------|-----------|
| 1 | maternal plasma bile acids diverged to **negative** values | the plasma volume had been written as the pure plasma volume (4 L), making it a compartment with a 7-minute half-life. Bile acids are over 95% albumin-bound and the bound pool equilibrates with the interstitium and within body water, so the apparent volume of distribution (18 L) is the right one |
| 2 | circulating UDCA 60+ µmol/L (never observed), or UDCA **displacing endogenous bile acids and raising them** | caused by putting UDCA inside the shared BSEP denominator. Canalicular efflux of UDCA is mainly via MRP2 and does not compete |
| 3 | the genotype→phenotype mapping **exceeded the clinical range by an order of magnitude** (term TBA 200–860) | BSEP had been placed near Vmax (denB ~4.4), so a 1.3-fold reduction in transport capacity was amplified into a 6-fold hepatocellular load. It was moved away from saturation |
| 4 | the basolateral escape valve was already 76% engaged in normal pregnancy | because it had been written as a term proportional to FXR. There was no headroom left to express the 5- to 20-fold induction seen in human cholestatic liver. It was changed to a Hill function of FXR occupancy |
| 5 | UDCA made the fetal pool 43% more hydrophobic, so **the drug came out harmful** | there was no route for handling the lithocholic acid made by the intestinal flora. SULT2A1 sulfation was added |
| 6 | UDCA accumulated once the cholestasis was blocked, reaching 100+ µmol/L in the circulation | renal excretion of UDCA glucuronide is far better than that of the endogenous conjugates |
| 7 | cord blood bile acids **diverged without bound** to 700 µmol/L | the placental diffusion term had been written as one-way, maternal→fetal, as in some published models. Once the transporter saturates, that leaves the fetal compartment with no concentration-independent escape route. Diffusion is bidirectional |
| 8 | UDCA lowered pruritus by **3.3 cm** on a 10 cm scale | because autotaxin was driven by maternal bile acids (KA_CH = 17). The observed value in the trial is 0.7 cm, and circulating autotaxin does not fall with UDCA |
| 9 | ALT 61 U/L and 86% meconium staining in normal pregnancy | the injury and meconium hazard functions had no threshold and no gestational-age gating. The stillbirth hazard function was also accumulating before 24 weeks, so the background risk was absorbing four weeks' worth of the pre-viability period |

---

## 9. A complete list of the fitted constants

Only **8** constants in this model were fitted to data, and three of them came
from a single paper.

| Constant | What it was fitted to |
|------|-----------------|
| `HSB0` = 1.089e-5 /day | the 0.13% stillbirth rate in the <40 band of Ovadia 2019 |
| `HSBSC` = 3.190e-3 /day | the 3.44% in the ≥100 band of the same paper |
| `HN` = 1.55 | the 0.28% in the 40–99 band of the same paper (residual 0.45%) |
| `HPT0`, `HM0`, `HMSC` | the background rates of spontaneous preterm birth and meconium in normal pregnancy |
| `VSCALE` = 4.81 | VAS cm per unit of the central pruritus state |
| `WMORB` = 0.055 | the utility ratio of neonatal morbidity to stillbirth — **a judgement value, adjustable in tab 8** |

Everything else is either a transport, binding or turnover constant (within the
literature range, linked item by item in `icp_references.md`), or a scale set so
that normal pregnancy comes out at normal values. Fitting three data points with
three constants is not in itself evidence — the evidence is **the size of the
exponent that comes out**, and section 3 is the measurement of that exponent
against three different driving variables.

---

## 10. What the model misses and where it conflicts with the literature

The full table is in the last section of `icp_references.md`. The main points:

- **The timing of delivery in the 40–99 band** does not reproduce the guidelines
  (section 6).
- **`V_P` (maximal placental transport capacity)** is a modelling choice fitted to
  cord blood concentrations in normal pregnancy, because no directly measured
  human value could be found.
- The *ordering* of the **`WTOX` cytotoxicity weights** has a basis in
  cardiomyocyte experiments, but the *spacing* is a modelling choice.
- **`HLP` = 4** is not the cooperativity of a single binding event but the
  apparent value of the multi-stage LPAR–TRPV1–spinal GRPR amplification. Do not
  read it mechanistically.
- **Direct fetal cardiac protection by UDCA** is reported in the literature, but
  the model expresses it only as `WTOX[UDCA]=0.02` and gives it no separate term.
  There is **a possibility that the fetal protection afforded by UDCA is
  underestimated**.
- **The thrombotic risk of ICP**, **the long-term neurodevelopment of the
  children**, and **inter-individual differences in pruritus by MRGPRX4
  genotype** are all in the literature but outside the scope of the model.
- **Within-individual variability of maternal bile acids** is not represented.
  This is the most important direction of extension.
- **Vitamin K status exceeds 1** (around 1.4 in the UDCA scenarios). That UDCA
  enlarges the intraluminal bile acid pool and thereby improves micelle formation
  is a real direction in itself, but writing the absorption term as first-order
  produces the meaningless value of "a vitamin K status 44% better than normal".
  Read the direction only, not the magnitude.
- **The magnitude of the worsening of vitamin K by rifampicin** (VITK 0.16, INR
  1.28) is in the right direction but is quite likely exaggerated.

---

## 11. How to run it

```bash
# render the map
dot -Tsvg icp_qsp_model.dot -o icp_qsp_model.svg
dot -Tpng -Gdpi=150 icp_qsp_model.dot -o icp_qsp_model.png

# Python reference implementation (all scenarios, without R; --quick for just 6)
python3 icp_reference_model.py
python3 icp_calibration.py          # the whole fitting, ablation and trial reproduction
```

```r
# mrgsolve
library(mrgsolve); library(dplyr)
mod <- mread("icp_mrgsolve_model", ".")
out <- mod %>% param(GBSEP = 0.65, GMDR3 = 0.50, GSULT = 2.40) %>%
  mrgsim(end = 154, delta = 1) %>% as_tibble()
# The 35 scenarios are in the comments at the bottom of icp_mrgsolve_model.R.

# Shiny (10 tabs)
shiny::runApp("icp_shiny_app.R")
```

---

## ⚠️ Disclaimer

This is a qualitative to semi-quantitative QSP model for educational and research
purposes. It was constructed from the published literature but has not been
independently validated or certified, and **must not be used for actual clinical
decision-making, prescribing, or regulatory submission.** In particular, the
delivery-timing calculations in section 6 conflict in part with the guidelines,
and that discrepancy must not be read as evidence that the model is the side that
is right.
