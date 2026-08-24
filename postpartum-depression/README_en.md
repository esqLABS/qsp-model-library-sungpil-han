# Postpartum Depression QSP Model

**Peripartum-onset major depressive disorder · Quantitative Systems Pharmacology**

<a href="ppd_qsp_model_en.svg"><img src="ppd_qsp_model_en.png" width="760" alt="PPD QSP mechanistic map"></a>

| File | Contents |
|------|------|
| [`ppd_qsp_model_en.dot`](ppd_qsp_model_en.dot) · [`.svg`](ppd_qsp_model_en.svg) · [`.png`](ppd_qsp_model_en.png) | Mechanistic map — 14 clusters (+ legend), 207 nodes, 294 edges |
| [`ppd_mrgsolve_model_en.R`](ppd_mrgsolve_model_en.R) | mrgsolve ODE model — 38 compartments, 10 scenarios |
| [`ppd_shiny_app_en.R`](ppd_shiny_app_en.R) | Shiny dashboard — 9 tabs |
| [`ppd_reference_check.py`](ppd_reference_check.py) | An independent implementation of the same equations (pure standard-library RK4) |
| [`ppd_reference_output.txt`](ppd_reference_output.txt) | The **actual run output** of that script — the source of every number in this README |
| [`ppd_references_en.md`](ppd_references_en.md) | 109 references, every PMID looked up and verified directly on PubMed |

> **There is no R runtime in this repository environment.** Every number below was
> therefore obtained by **actually running** `ppd_reference_check.py` (the same 38
> differential equations and the same parameters, re-implemented in pure Python RK4),
> and its raw output is `ppd_reference_output.txt`. If the R file gives different
> values, it is the R file that is wrong.

---

## 1. The single idea of this model

Postpartum depression has not been modelled here as "allopregnanolone falls and
brexanolone puts it back". This model is organised around one structural fact: that
this is **a product of two factors with different time constants**.

```
G_tonic(t) = R_δ(t) × [ 1 + Emax · PAM(t)^h / (EC50^h + PAM(t)^h) ]
             ~~~~~~     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             receptor arm                 ligand arm
             t½ ≈ 120 h                   t½ ≈ 12 h
```

During pregnancy the receptor arm is down-regulated in accordance with a
**homeostatic set point**, so the product stays exactly at the non-pregnant set point.

```
R_δ,ss(PAM) = G_target / (1 + potentiation(PAM))
```

So **every steady state of this system has identical tonic inhibition** (Table 1 of
the run output):

| State | PAM (nM) | potentiation | R_δ set point | **G_tonic** |
|------|---------:|-------------:|-----------:|------------:|
| Non-pregnant | 6.0 | 0.037 | 0.964 | **1.000** |
| Mid-pregnancy | 60.0 | 0.687 | 0.593 | **1.000** |
| Term | 160.0 | 1.498 | 0.400 | **1.000** |
| Postpartum day 3 | 4.8 | 0.027 | 0.973 | **1.000** |

That is, **postpartum depression is not a steady state of this system.** No
steady-state analysis whatever can generate it. It is **a transient between two
identical steady states**. When PAM has already fallen to 4.8 nM while R_δ is still at
its term value of 0.400, the product is 0.411, which is a **59 % loss of tonic
inhibition** relative to the set point — a figure neither factor shows on its own.

---

## 2. Results the model **computed** (as opposed to assumed)

### Result 1 — the window of onset is derived from the two time constants alone

Give it only the two numbers, a ligand half-life of 12 h and a receptor half-life of
120 h, and the duration of the deficit follows:

| Deficit in tonic inhibition | Duration |
|---|---|
| > 25 % | 149 h = **6.2 days** |
| > 10 % | 308 h = **12.8 days** |
| > 5 % | 428 h = **17.8 days** |

The code has never been told "4–6 weeks postpartum". Because **the width of the window
is set by KR and KR has never been measured in humans** (assumption A8), that width
should be read only as an order-of-magnitude estimate.

### Result 2 — an identical hormonal trajectory, and one vulnerability parameter separates "the blues" from "an episode"

The steroid collapse is exactly the same in every row; only V differs.

| Phenotype | V | E/I threshold | Peak EPDS | HAM-D at day 42 |
|---|---:|---:|---:|---:|
| Resilient (good support) | 0.80 | 3.29 | 2.3 | 4.0 |
| Population mean | 1.00 | 2.63 | **7.1** | 7.0 |
| At risk (history of depression) | 1.40 | 1.88 | 16.1 | 23.9 |
| High risk (history of PPD) | 1.70 | 1.55 | 19.7 | **27.9** |
| Very high risk | 2.00 | 1.32 | 20.6 | 29.6 |

The woman with mean V goes through a self-limiting nadir that stays **below** the EPDS
screening cut-off (12/13) and recovers by herself — this is what this model calls
**the baby blues**. The high-risk woman **latches** into a persistent episode. It is a
formal version of Bloch's 2000 experiment (PMID 10831472), in which an identical
steroid withdrawal provoked symptoms only in women with a history of PPD.

### Result 3 — that the 60 → 90 µg/kg/h dose-response is flat is a prediction of the Hill function

| Infusion rate µg/kg/h | Css ng/mL | Brain nM-eq | potentiation |
|---:|---:|---:|---:|
| 30 | 30.0 | 188 | 1.632 |
| 60 | 60.0 | 377 | 2.081 |
| 90 | 90.0 | 565 | 2.244 |
| 120 | 120.0 | 754 | 2.323 |

Raise the exposure by 50 % and potentiation rises by only **7.8 %**. Site 2, on the
other hand — direct channel activation, which is what produces the sedation
(PMID 17108970) — has no saturation in this range. So "efficacy is flat while only the
sedation increases" is a **prediction**, and it is exactly what phase 3 study 1
reported as BRX60 −19.5 / BRX90 −17.7 (PMID 30177236).

The same arithmetic says one thing more: the steady-state exposure at 90 µg/kg/h is
283 nM, **3.5 times** the term allopregnanolone the model uses (80 nM). The infusion
does not "restore pregnancy levels" — it **exceeds** them severalfold, which is why
the dose-limiting toxicity is sedation and not efficacy.

### Result 4 — comparison with the clinical trial endpoints

Enrolment at postpartum day 21, V = 1.70. Change in HAM-D17:

| Arm | Baseline | 60 h | Day 3 | Day 15 | Day 28 | Day 45 |
|---|---:|---:|---:|---:|---:|---:|
| Placebo (inpatient, 60 h) | 28.8 | **−15.6** | −15.3 | −13.8 | −14.5 | −15.1 |
| Brexanolone 60 µg/kg/h | 28.8 | **−19.2** | −19.8 | −13.7 | −14.3 | −15.1 |
| Brexanolone 90 µg/kg/h | 28.8 | −19.4 | −20.1 | −13.8 | −14.3 | −15.1 |
| Placebo (outpatient) | 28.8 | −5.8 | −6.6 | **−13.1** | −14.3 | −15.1 |
| Zuranolone 50 mg × 14 days | 28.8 | −6.4 | −7.4 | **−15.5** | −13.4 | −14.7 |
| Zuranolone 30 mg × 14 days | 28.8 | −6.2 | −7.1 | −14.5 | −13.5 | −14.8 |
| Zuranolone 50 mg × 3 days | 28.8 | −6.4 | −7.4 | −12.8 | −14.1 | −15.0 |
| Sertraline 50→100 mg | 28.8 | −5.8 | −6.6 | −13.9 | −17.5 | −21.8 |
| Esketamine 0.25 mg/kg once | 28.8 | −5.8 | −6.6 | −13.1 | −14.4 | −15.1 |
| Sleep protection alone | 28.8 | −5.8 | −6.6 | −13.3 | −14.7 | −15.5 |
| CBT/IPT alone | 28.8 | −6.0 | −6.8 | −14.3 | −16.0 | −17.0 |

**Observed values (all verified directly on PubMed):**

| Trial | Endpoint | Observed |
|---|---|---|
| Kanes 2017 phase 2 (PMID 28619476) | 60 h | brexanolone −21.0 / placebo −8.8 |
| Meltzer-Brody 2018 study 1 (PMID 30177236) | 60 h | BRX60 −19.5 / BRX90 −17.7 / placebo −14.0 |
| Meltzer-Brody 2018 study 2 | 60 h | BRX90 −14.6 / placebo −12.1 |
| ROBIN 2021 (PMID 34190962) | Day 15 | zuranolone 30 mg −17.8 / placebo −13.6 |
| SKYLARK 2023 (PMID 37491938) | Day 15 | zuranolone 50 mg −15.6 / placebo −11.6 |

### Result 5 — half of the enormous placebo response is mechanism (and half is not)

The placebo arms of these trials improved by 8.8–14.0 points, which is **larger than
the drug-placebo difference in every one of the trials.** In the model, part of that is
not noise but a state variable.

- **Reproduced:** the 60-hour inpatient placebo response is larger than the 15-day
  outpatient placebo response **at every enrolment time point**. There is only one
  reason — 60 hours of continuous nursing care protects the mother's sleep, and sleep
  is an active input in this model. The brexanolone trials were inpatient and their
  placebo arm recorded −14.0 at 60 hours; the zuranolone trials were outpatient and
  took 15 days to reach −11.6. The model was never taught that ordering.
- **Not reproduced:** after calibration the placebo response is almost **flat** with
  respect to the time of enrolment (about 1 point from day 7 to day 180), even though
  the δ receptor pool moves substantially, from 0.66 to 0.97. The fitted care term
  overwhelms the mechanistic one. The attractive hypothesis that "receptor recovery
  explains a substantial part of the placebo response" **does not hold** at these
  weights.

### Result 6 — the speed of onset distinguishes what efficacy at week 6 cannot

| Arm | Day 1 | Day 3 | Day 7 | Day 15 |
|---|---:|---:|---:|---:|
| Brexanolone 60 (GABA_A) | −9.9 | **−19.8** | −15.8 | −13.7 |
| Zuranolone 50 (GABA_A) | −2.9 | −7.4 | −12.9 | −15.5 |
| Esketamine (glutamate) | −2.8 | −6.6 | −10.5 | −13.1 |
| Sertraline (monoamine) | −2.8 | **−6.6** | −10.7 | −13.9 |

Sertraline has not been made weak. SERT occupancy is immediate (78 % at day 15) and
its final effect is competitive. The delay comes from **two slow states placed in
series** — 5-HT₁A autoreceptor desensitisation (t½ 10 days) and structural plasticity
(t½ 10 days).

### Result 7 — maternal exposure and infant exposure are different questions

| Regimen | Maternal Cmax ng/mL | Milk ng/mL | RID % | Infant plasma nM |
|---|---:|---:|---:|---:|
| Brexanolone 90 (60 h) | 88.1 | 132.2 | **0.65** | 0.50 |
| Zuranolone 50 mg/day (14 days) | 246.3 | 369.5 | 5.47 | 1.15 |
| Sertraline 100 mg/day | 54.5 | 98.2 | 0.80 | 0.33 |

The infant plasma concentrations are **two orders of magnitude below** the model's own
potentiation EC50 (120 nM). That is, the risk of maternal sedation and the risk of
infant sedation are not the same question, and this model separates them explicitly
rather than inferring the second from the maternal dose. Because measured brexanolone
milk data exist (PMID 35869362), the RID column is **falsifiable**.

---

## 3. Where this model **failed** — and what that failure tells us

Recorded honestly. The model reproduces (a) the severity at enrolment, (b) the 60-hour
brexanolone response, (c) the 15-day zuranolone response and (d) the magnitude of both
placebo responses. But **it fails to reproduce the maintenance of the drug-placebo
separation at 28–45 days** — every arm converges on a HAM-D of about 14. Once the
receptor window has closed, the late attractor of this system is a partially recovered
state, and no drug parameter can change the **position** of that attractor, only the
**speed** at which it is reached.

So the two obvious repairs were implemented as switches and **run**, rather than argued
about:

| Setting | Baseline | brx 60 h | brx day 30 | zur day 15 |
|---|---:|---:|---:|---:|
| Reference (both off) | 28.8 | −19.2 | −14.4 | −15.5 |
| KR_BOOST = 1.0 (the drug accelerates receptor recovery) | 28.7 | −17.9 | −14.2 | −14.0 |
| KR_BOOST = 3.0 | 28.7 | −16.7 | −14.1 | −13.3 |
| W_PAM_BDNF = 0.15 (giving the drug a plasticity arm) | 25.5 | −17.2 | −11.7 | −14.6 |
| W_PAM_BDNF = 0.35 | 23.9 | −16.2 | −10.7 | −16.6 |

**Both make the fit worse, and for the same reason.**

> Brexanolone **is allopregnanolone itself** (assumption A2 — a fact, not a modelling
> choice). So the drug and the endogenous ligand enter through **the same potentiation
> term**. Any parameter that gives the drug extra credit **gives pregnancy exactly the
> same extra credit.** The baseline column of the table above is precisely that
> (28.8 → 25.5 → 23.9).

KR_BOOST fails more directly still. **In a homeostatic set-point model a PAM cannot
accelerate the recovery of the receptor pool** — because a PAM lowers the very set
point the pool is chasing. Raise the speed of the chase and it merely gets to a lower
target faster.

The conclusion: the missing mechanism of persistence **cannot lie on the
allopregnanolone-potentiation axis.** It has to be something that the drug does and
pregnancy does not. One candidate that a trial could distinguish: a lasting change in
the sleep-behaviour loop that 60 hours of "being functional" buys (which in this model
would be expressed as a permanent reduction in A_SYMP or WAKE rather than as a receptor
effect).

Incidentally, this is also why **no prediction about the duration of dosing comes out
either.** A 3-day zuranolone course is inferior at day 15 (it has already stopped) but
indistinguishable at day 45. This model supports only the weak claim that "the course
has to continue if the acute gain is to be maintained"; it cannot say "a short course
relapses".

### The uncomfortable fact the sensitivity analysis states

The ±30 % univariate sensitivity ranking for the 15-day zuranolone endpoint:

| Rank | Parameter | max |Δ| (HAM-D) | Character |
|---|---|---:|---|
| 1 | THR0 | 2.63 | the E/I threshold — the disease boundary the model asserts ✓ |
| 2 | **K_NSP** | 2.19 | **a fitted non-specific care term (not mechanism)** |
| 3 | KDEC_SD | 1.76 | repayment of sleep debt |
| 4 | W_SD | 1.57 | sleep debt → excitatory load |
| 5–8 | EMAX_PAM · KOFF · EC50_PAM · ZUR_EQ | 0.79–0.96 | the drug arm (saturated, hence insensitive) |
| 10 | KR | 0.62 | sets the **width** of the window but the day-15 value less so |

That the second place is a non-mechanistic fitted parameter is a warning label: a
substantial part of this endpoint is calibration rather than mechanism. Conversely,
that the potency parameters of the drug itself sit mid-table is consistent with the
saturation of Result 3.

---

## 4. Model structure

**38 ODEs.** Placental and neurosteroid synthesis (PLAC, P4, DHP, ALLOP, ALLOB, E2,
PCRH) → the HPA axis (HCRH, ACTH, CORT, GRFN) → **receptor plasticity (RD, RG,
KCC2)** → monoamine and plasticity (MAOA, FIVEHT, AUTO, BDNF, SYN) → inflammation and
kynurenine (INFL, KYNR) → sleep (SLP, SDEBT) → excitation/inhibition (EXC) →
symptoms (SYMP) → bonding (BONDI) → drug PK (BRX1/2, ZURA/1/2, SERA/C, ESKC, AMPAS) →
psychotherapy (CBTP) → lactational transfer (MILKD, INFP).

**Three central couplings:**

1. **The multiplicative coupling** `G_tonic = R_δ × (1 + potentiation)` — this is what
   creates the transient.
2. **The homeostatic set point** `R_δ,ss = G_target/(1+potentiation)` — this is what
   makes every steady state identical, and thereby excludes the disease from any
   steady state.
3. **The closed vicious circle** symptoms → insomnia → sleep debt → excitatory load →
   symptoms. Together with the threshold, this loop creates two stable states
   (depressed/recovered), and that is why a temporary intervention can flip the state.
   The observation that it is a change in sleep rather than in hormones that predicts
   the timing of a PPD relapse (PMID 20708275) is the grounds for placing this loop as
   the proximal driver.

**Important:** this model **operates near a bifurcation** (assumption A7). That is
deliberate — it is the way it explains why a modest intervention can flip the outcome
and why the placebo trajectories are so gradual and so variable. But it also means that
the late endpoints are intrinsically sensitive to THR0 and KR. The sensitivity table
quantifies this rather than hiding it.

---

## 5. Assumptions — to be read before quoting any number

| # | Assumption | Character |
|---|------|------|
| A1 | Neurosteroids are handled as **total** plasma/brain concentrations, and the EC50 is in the same units → protein binding of over 99 % and brain partitioning are **absorbed into** the EC50. The EC50 here is not comparable with a patch-clamp EC50. | structural simplification |
| A2 | **Brexanolone = allopregnanolone.** Partition coefficient and potency are assigned by identity, not fitted. As §3 shows, this is not a convenience but a **constraint**. | fact |
| A3 | ZUR_EQ is the **only** parameter fitted to a zuranolone endpoint. | calibration |
| A4 | Term ALLO ≈ 80 nM, postpartum nadir ≈ 2 nM. There is a several-fold disagreement between assays (PMID 11238543), and the 2025 individual-patient-data meta-analysis (PMID 39511449) reports that **absolute concentrations do not separate cases from controls** — which is precisely the premise of this model. What was fitted is therefore the **ratio** of about 40-fold. | literature-based |
| A5 | K_CARE and K_NSP are **non-mechanistic** terms for expectancy and structured clinical contact. Fitted once on the placebo arms and then held fixed in every active arm. The model does not explain the expectancy effect. | an explicit fudge |
| A6 | The reference enrolment time is postpartum day 21. Real trials enrolled up to 6–12 months. | a scenario choice |
| A7 | The model operates near a bifurcation → the day-45 endpoints are intrinsically sensitive. | structural |
| A8 | **KR (recovery of the δ pool, t½ ≈ 5 days) has never been measured in humans.** Inferred from rodents (PMID 9789080). The entire width of the window rests on this one number. | a free parameter |

---

## 6. How to reproduce

```bash
# render the mechanistic map
dot -Tsvg  ppd_qsp_model_en.dot -o ppd_qsp_model_en.svg
dot -Tpng -Gdpi=150 ppd_qsp_model_en.dot -o ppd_qsp_model_en.png

# regenerate every number in this README (no R needed, about 15 minutes)
python3 ppd_reference_check.py > ppd_reference_output.txt

# the mrgsolve model (R required)
Rscript ppd_mrgsolve_model_en.R

# the interactive dashboard
Rscript -e 'shiny::runApp("ppd_shiny_app_en.R")'
```

---

## ⚠️ Disclaimer

**This is a QSP model for educational and research purposes.** It was built from the
public literature and the published endpoints of four clinical trials, but it has not
been validated against individual patient data, and as recorded in §3 there are
observations it fails to reproduce. **It cannot be used for clinical decision-making,
prescribing or regulatory submission.** Postpartum psychosis (1–2 per 1000) and
suicidal ideation are separate emergencies and lie outside the scope of this model — if
either is suspected, immediate specialist assessment is required.
