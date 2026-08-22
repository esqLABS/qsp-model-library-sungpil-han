# Osteosarcoma — QSP Model

> **Cure is a Poisson bet on lesions nobody can see, and it is paid for out of
> an exposure budget with three hard organ ceilings — one of which is the organ
> that clears the drug.**

<p align="center">
  <a href="../../../osteosarcoma/osa_qsp_model.svg">
    <img src="../../../osteosarcoma/osa_qsp_model.png" width="900" alt="Osteosarcoma QSP mechanistic map">
  </a>
</p>

| File | Contents |
|------|------|
| [`osa_qsp_model.dot`](../../../osteosarcoma/osa_qsp_model.dot) · [`.svg`](../../../osteosarcoma/osa_qsp_model.svg) · [`.png`](../../../osteosarcoma/osa_qsp_model.png) | Mechanistic map — 193 nodes, 15 clusters, 258 edges |
| [`osa_mrgsolve_model.R`](osa_mrgsolve_model.R) | **The model proper** — 55 ODE compartments, 141 parameters, 20 scenarios, protocol gate layer, virtual population layer |
| [`osa_shiny_app.R`](../../../osteosarcoma/osa_shiny_app.R) | Interactive dashboard (8 tabs) |
| [`osa_references.md`](../../../osteosarcoma/osa_references.md) | 64 references + calibration comparison table |
| [`osa_reference_model.py`](../../../osteosarcoma/osa_reference_model.py) | Python/scipy reimplementation of the same equations (for numerical verification in a container with no R) |
| [`osa_reference_output.txt`](../../../osteosarcoma/osa_reference_output.txt) | The full text of the actual run output (509 lines) |
| [`osa_scenario_results.json`](../../../osteosarcoma/osa_scenario_results.json) · [`osa_analysis_results.json`](../../../osteosarcoma/osa_analysis_results.json) | Machine-readable results |

Every number below is a value that `osa_reference_model.py` actually integrated
and printed, and it can be checked exactly as it stands in
`osa_reference_output.txt`. These are calculations, not claims.

---

## 1. Why this disease needs a model

The 5-year survival of osteosarcoma has **barely moved** since the early 1980s
(Mirabello 2009, SEER). The large randomised trials run in the interval have
been remarkably consistently negative:

- **EURAMOS-1 MAPIE** (618 poor responders randomised): EFS **HR 0.98**
  (95% CI 0.78–1.23)
- **EURAMOS-1 MAPifn** (good responders): no benefit
- **OS2006** (zoledronic acid, 318 patients): no EFS benefit, **local control
  actually worse**
- **EOI intensification trial**: **histological response improved, survival did
  not**
- **SARC028** (pembrolizumab): 5% response rate in osteosarcoma

Taken one at a time, each reads as "that drug did not work". The point of this
model is that they are **the same result coming out of the same structure**. And
that structure can be written in two sentences.

---

## 2. The organising thesis

### 2.1 Cure is a Poisson bet on lesions nobody can see

Let λ₀ be the number of occult micrometastatic lesions already present at
diagnosis, n₀ the number of clonogenic cells in one lesion, and K the log kill
delivered by treatment. Then

```
P(cure) = exp( −λ₀ · [ 1 − exp( −n₀ · e^(−K) ) ] )
           \______/     \________________________/
            how many        whether this one
            lesions         lesion died
```

**λ₀ is not a free parameter.** The cure rate with amputation alone is
exp(−λ₀), and in the literature that value is 0.15–0.20 (Link 1986,
Eilber 1987). λ₀ is therefore fixed at 1.80, and the model gives 0.168. From
that one anchor the rest follows — to reach the survival of 0.60 that MAP
actually achieves, K has to be **6.6 log10**, and the model delivers 6.61.

The reason this form matters is that K sits inside a **double exponential**. The
value of adding a little exposure changes sharply around the operating point,
and averaged over a whole population it almost vanishes. Why forty years of
intensification attempts produced a flat curve comes out of this (§3.G).

### 2.2 Methotrexate destroys the organ that excretes it

86% of the dose leaves by the kidney. The intratubular concentration is **the
excretion rate divided by the urine flow**, and solubility depends
exponentially on pH:

```
C_tub = CL_ren(KID) · C_plasma / UF(KID)
S(pH) = 0.86 · 10^(0.682·(pH − 5))  mM        [0.39 / 1.55 / 9.04 mg/mL at pH 5 / 6 / 7]
```

Whatever exceeds S crystallises and blocks the tubule. Obstruction lowers
**both** `CL_ren` and `UF`, and `C_tub` therefore **rises**. This is not a risk
factor but a **positive feedback loop**, and so what exists is not a risk curve
but a **critical point**.

---

## 3. Seven results — all of them calculations, none of them claims

### A. It is not a guideline, it is a bifurcation point

Integrating a single 12 g/m² course while varying only the urine pH:

| Threshold | Derived urine pH |
|---|---|
| tubular fluid begins to be supersaturated | **7.30** |
| delayed excretion (C48 > 1 µM) | **7.15** |
| AKI (loss of 25% or more of GFR) | **6.95** |

The protocol wording "maintain urine pH ≥ 7.0 and hydrate at 3 L/m²/day" is not
a convention. It is **the point at which the sign of this loop changes**. And
pH 7.0 lies **0.30 units inside** the supersaturation boundary — this is why
delayed excretion keeps occurring even under complete protocol compliance, and
why the response to a rising creatinine is *more bicarbonate, not less
methotrexate*.

### B. One unit of urine pH is worth 4.81× the hydration

`C_tub` is proportional to dose/flow and `S` is proportional to
10^(0.682·pH). Any factor f applied to `C_tub` therefore has the value of
**log₁₀(f)/0.682 pH units**:

| Intervention | pH equivalent |
|---|---|
| 2× hydration | 0.44 pH |
| 4× hydration | 0.88 pH |
| MTX dose 1/2 | 0.44 pH |
| MTX dose 1/4 | 0.88 pH |
| **urine pH 1.00 unit** | **4.81× on C_tub = 4.81× hydration = a 4.81-fold dose reduction** |

Confirmed by simulation (pH fixed at 6.6, single 12 g/m²):

| Hydration (L/m²/day) | C48 (µM) | Final eGFR | Delayed excretion |
|---|---|---|---|
| 1.5 | 14.67 | 55.7% | yes |
| 3.0 (standard) | 9.16 | 62.9% | yes |
| 6.0 | 1.86 | 81.2% | yes |
| 12.0 | 0.56 | 99.8% | **no** |
| **hydration left at 3.0, pH alone raised to 7.3** | **0.56** | **99.5%** | **no** |

**0.7 pH units does the same work as 4× hydration.** Alkalinisation is not
adjunctive supportive care. Arithmetically it is the largest lever on the safety
of this regimen, larger than the dose itself.

### C. The worse the kidney, the larger the loop gain

Same course, same pH 6.6, varying only the starting GFR:

| Starting eGFR | Supersaturation ratio C/S | C48 (µM) | AUC (mM·h) | Final eGFR |
|---|---|---|---|---|
| 100% | 3.45 | 9.16 | 8.58 | 62.9% |
| 90% | 3.55 | 11.88 | 9.30 | 60.9% |
| 80% | 3.66 | 15.49 | 10.14 | 59.0% |
| 70% | 3.77 | 20.24 | 11.12 | 57.2% |
| 60% | 3.90 | 26.55 | 12.28 | 55.4% |
| 50% | 4.03 | 34.92 | 13.66 | 53.7% |

Halving the starting GFR makes **the AUC of the same dose 1.59× and C48 3.8×**,
and the kidney ends at 53.7% instead of 62.9%. (The absolute loss is in fact
smaller the worse the starting point — there is less to lose. The damage
converges towards a common floor, and the point is that this floor is exactly
where 12 courses take you.) That is, 12 courses of methotrexate and 4 of
cisplatin are **not additively but multiplicatively nephrotoxic, through one
shared state variable**, and the nth course rewrites the pharmacokinetics of
courses n+1 … 12.

### D. Rescue therapy is diagnostically gated so that it arrives after the exposure it would have prevented is over

Glucarpidase is triggered by a **measured concentration**, and that
concentration is measured at 24 / 48 / 72 hours. But the exposure integral is
loaded at the front:

(urine pH 6.4, single 12 g/m²)

| Time of administration | C48 (µM) | C72 (µM) | AUC (mM·h) | AUC removed | Final eGFR |
|---|---|---|---|---|---|
| none | 12.07 | 1.899 | 8.92 | — | 59.1% |
| 12 h | 1.51 | 0.251 | 6.34 | 29% | 63.2% |
| 24 h | 1.02 | 0.367 | 7.99 | 10% | 60.5% |
| 36 h | 0.20 | 0.343 | 8.55 | 4% | 59.6% |
| **48 h (standard trigger)** | 12.07 | 0.177 | 8.78 | **2%** | 59.3% |
| 72 h | 12.07 | 1.899 | 8.90 | 0% | 59.2% |

A drug arriving at the standard trigger time (48–72 hours) can remove **2% of
what is left** and cannot touch **the 98% that has already arrived**. And at no
time point does the tubular injury integral change — the damage was already
paid for in the first 24 hours, before any monitoring threshold had been
crossed. **There is no post hoc rescue of the tubule.** The only real lever is
the one applied *before* the infusion starts (§B).

### E. The Huvos grade is a prognostic factor and manipulable at the same time

Running MAP while varying only the intrinsic resistance RES₀:

| RES₀ | Huvos % | Call | log10 kill | Cure rate | HR (vs 0.02) |
|---|---|---|---|---|---|
| 0.00 | 95.5 | good | 9.01 | 0.884 | 0.24 |
| 0.02 | 90.3 | good | 6.61 | 0.602 | 1.00 |
| 0.05 | 86.1 | poor | 5.66 | 0.408 | 1.77 |
| 0.10 | 82.0 | poor | 4.91 | 0.263 | 2.62 |
| 0.15 | 79.3 | poor | 4.44 | 0.207 | 3.10 |
| 0.20 | 77.1 | poor | 4.09 | 0.185 | 3.32 |
| 0.30 | 73.5 | poor | 3.56 | 0.172 | 3.46 |
| 0.45 | 68.8 | poor | 2.93 | 0.170 | 3.49 |

**r(Huvos %, micrometastatic log kill) = 0.97.** That is why the Huvos grade is
the most powerful prognostic factor in this disease.

But — Huvos is an integral over the **primary lesion**, and there the drug
penetration penalty (mineralised matrix) and the matrix-derived TGF-β growth
boost both apply. Cure is an integral over the **lung**, where neither applies.
**Move the same total drug from post-operative to pre-operative and Huvos moves
by +3.5 points while the micrometastatic log kill moves by +0.00.**

That is, the metric can improve without the patient improving. **A metric that
is prognostic and manipulable at once is the worst kind of metric on which to
base a treatment decision**, and response-adapted intensification does exactly
that.

### F. MAPIE's negative result is the arithmetic of its own design

Every kill term carries a **(1 − RES)** factor. A "poor responder" is by
definition a patient in whom that factor is small. And MAPIE bought its five
cycles of ifosfamide/etoposide **by giving up one doxorubicin and one cisplatin
out of a budget already against the wall**.

| | log10 kill | Cure rate | TRM | Survival | HR | Schedule adherence | eGFR nadir | ANC nadir |
|---|---|---|---|---|---|---|---|---|
| Good responder (RES₀ 0.02) — MAP | 6.61 | 0.602 | 1.0% | 0.596 | 1.00 | 100% | 92% | 0.39 |
| Good responder (RES₀ 0.02) — MAPIE | 6.93 | 0.659 | 3.0% | 0.639 | **0.86** | 76% | 62% | 0.06 |
| **Poor responder (RES₀ 0.20) — MAP** | 4.09 | 0.185 | 1.0% | 0.183 | 1.00 | 100% | 92% | 0.39 |
| **Poor responder (RES₀ 0.20) — MAPIE** | 4.30 | 0.197 | 3.0% | 0.191 | **0.98** | 76% | 62% | 0.06 |

**In the population in which the trial never randomised intensification (the
good responders) the HR is 0.86; in the population it did randomise the HR is
0.98.** The value EURAMOS-1 observed is **0.98** (0.78–1.23). The negative
result is not a surprise about the drugs but the arithmetic of the design.

### G. Three mechanisms behind the flat curve — none of them is the drug

**(i) Undetectability.** The between-patient CV of the achieved log kill is of
the order of 25% (SLC19A1/ABCC2 genotype, MTX clearance CV ~30%, MTHFR, the
dose intensity actually delivered after the gates, intrinsic sensitivity —
**none of which has ever been randomised**). What a realistic intensification
buys against that spread:

| Δ log10 kill | Increase in cure rate | HR | Patients per arm for 80% power |
|---|---|---|---|
| +0.10 | +0.019 | 0.94 | 10,827 |
| **+0.20** | **+0.037** | **0.88** | **2,764** |
| +0.30 | +0.054 | 0.83 | 1,256 |
| +0.50 | +0.087 | 0.73 | 475 |

MAPIE's actual net change is +0.2 log10. EURAMOS-1 randomised 618 patients per
arm. **That trial was not underpowered by accident; it was underpowered because
of the shape of the dose-response surface it was mounted on.**

How decisive the spread is can be seen like this — applying the same +0.5 log at
different CVs:

| CV (SD, nats) | Baseline cure rate | After +0.5 log | Absolute gain |
|---|---|---|---|
| 0.02 (0.30) | 0.675 | 0.871 | **+19.6 %p** |
| 0.10 (1.52) | 0.631 | 0.792 | +16.0 %p |
| 0.25 (3.81) | 0.602 | 0.689 | +8.7 %p |
| 0.65 (9.90) | 0.588 | 0.623 | **+3.5 %p** |

The same drug, the same half log. The only thing that changed is the spread.

**(ii) An indifference curve.** This one is worse. Halving the hydration
**increases exposure through delayed excretion and genuinely buys +0.4 log10 of
additional kill**. And the price paid for it is 12% treatment-related
mortality:

| | cure | TRM | **Survival** |
|---|---|---|---|
| Standard MAP | 0.602 | 1.0% | **0.596** |
| Half hydration | 0.674 | 11.7% | **0.595** |
| Urine pH 6.0 | 0.592 | 38.1% | **0.366** |

**The regimen sits on an indifference curve.** An intervention that moves along
the curve reads as negative however large it is. So one has to read **the
"survival" column, not the "kill" column** — several arms buy kill by wrecking
the kidney, the kill is real, and they still lose.

**Only one intervention moves the wall itself.** Dexrazoxane raises the cardiac
ceiling, and as a result doxorubicin 600 mg/m² + dexrazoxane gives the best
performance of all 20 arms (survival **0.649**, HR **0.83**). Everything else
redistributes the budget.

### (iii) Not the mean and not the tail — the shoulder

P(cure) **saturates at both ends**. The marginal value of 1 additional log of
kill is therefore 0 in a patient who will already be cured, and 0 in a patient
who cannot be saved by any achievable dose. Spending the same total additional
kill (a budget equivalent to +0.25 log for everyone) on different parts of the
distribution gives an entirely different answer:

| Who receives it | Per patient | Cure rate | Gain | HR |
|---|---|---|---|---|
| Everyone (uniform intensification) | +0.25 log | 0.649 | +0.048 | 0.85 |
| Bottom quartile (the worst responders) | +1.00 log | 0.632 | +0.030 | 0.90 |
| **Second quartile (those who narrowly missed)** | +1.00 log | **0.717** | **+0.116** | **0.65** |
| Third quartile | +1.00 log | 0.632 | +0.031 | 0.90 |
| Top quartile (already cured) | +1.00 log | 0.604 | +0.002 | 0.99 |
| Middle half (25th–75th percentile) | +0.50 log | 0.688 | +0.087 | 0.74 |

**The same total drug is worth 2.4× as much given to the second quartile as
sprayed over everyone, and given to the worst quartile it is worse than uniform
intensification** (+0.030 vs +0.048). A patient who is 4 log short is not saved
by one more log.

Locating the marginal value by decile:

| Decile | Achieved log10 kill | Marginal value |
|---|---|---|
| 1 | 0.74 – 4.50 | 0.000 |
| 2 | 4.50 – 5.23 | 0.002 |
| 3 | 5.23 – 5.75 | 0.129 |
| 4 | 5.75 – 6.18 | 0.499 |
| **5** | **6.18 – 6.60** | **0.592** ← maximum |
| 6 | 6.60 – 7.04 | 0.412 |
| 7 | 7.04 – 7.47 | 0.186 |
| 8 | 7.47 – 7.98 | 0.072 |
| 9 | 7.98 – 8.67 | 0.022 |
| 10 | 8.67 – 12.52 | 0.003 |

**The patients worth intensifying are those already close to the threshold who
just failed to cross it.**

And placing this alongside §E closes the argument. Because the necrosis fraction
correlates with the achieved log kill at r = 0.97, **this metric can identify
that band** — the "narrowly missed" patients at roughly 85–89% necrosis. Not
the 40%-necrosis patients. EURAMOS-1 used the metric accurately as a prognostic
tool and then spent the intensification budget on **the group that metric had
identified as the least salvageable**. **The band the arithmetic points to has
never once been randomised.**

---

## 4. Model structure

### 4.1 State variables (55 ODE compartments)

| Block | Compartments |
|---|---|
| MTX PK + renal loop | `A1` `A2` `A3` `PPT` `KID` |
| Folate · rescue therapy | `MPGT` `MPGM` `LV` `RFT` `RFM` `GLUC` |
| Doxorubicin · heart | `DOXC` `DOXP` `DOXOL` `CMV` `CUMDOX` `DEXR` |
| Cisplatin · tubule · cochlea | `CISC` `ADDT` `ADDK` `ADDC` `HL` `MG` |
| Ifosfamide · etoposide | `IFOA` `ETOC` |
| Tumour | `PRIM` `NEC` `RES` `LMET` `KI_MTX` `KI_DOX` `KI_CIS` `KI_IE` `KI_IMM` |
| Bone remodelling vicious cycle | `OCL` `TGFM` `RKL` `DENO` `ZOL` |
| Vasculature · hypoxia | `VASC` `HYP` |
| Immunity | `CTL` `M2` `MIFA` |
| Myelosuppression (Friberg) | `PROL` `TR1` `TR2` `TR3` `CIRC` |
| Mucosa · exposure · biomarkers · mortality | `MUC` `EMTX` `EDOX` `ECIS` `EIFO` `ALP` `TRM` |

Two key design decisions:

- **`TGFM` multiplies the growth term** (`grow = KG·(1 + 0.45·TGFM)·gompertz`),
  not the survival term. During chemotherapy the kill terms are one to two
  orders of magnitude larger than the growth term, so suppressing the vicious
  cycle moves osteolysis, ALP and imaging, and does not move survival. The
  zoledronic acid and denosumab results come out of here.
- **`PEN` (the penetration penalty) applies to the primary lesion only.**
  Micrometastases are fully perfused. This is why Huvos and cure are *different*
  integrals of the same parameters.

### 4.2 The protocol gate layer (go / no-go)

A real protocol does not give a scheduled course to a patient who has not
recovered. This layer is the mechanism by which **intensification spends the
budget instead of adding to it**:

| Agent | ANC minimum | Mucositis | Kidney | Other |
|---|---|---|---|---|
| MTX | ≥ 0.75 | < grade 3 | ≥ 0.60 | hold if KID < 0.55, resume above 0.70 |
| DOX | ≥ 1.00 | < grade 3 | — | LVEF ≥ 50% |
| CIS | ≥ 1.00 | < grade 3 | ≥ 0.60 | |
| IE | ≥ 1.00 | < grade 3 | ≥ 0.60 | |

A deferred course returns two weeks later at **−25% dose**, and after that it is
missed altogether. Result: schedule adherence MAP 100% vs MAPIE 76% —
reproducing, without assuming it, the same **25-point gap** as EURAMOS-1's
76% / 51%.

### 4.3 Treatment-related mortality (TRM)

Without this the model rewards nephrotoxicity (§G-ii). The infection term is
the **product** of the depth of neutropenia and the mucositis grade — mucosal
damage with neutrophils present is survived, neutropenia with intact mucosa is
survived, and the two together are not.

```
dTRM/dt = KTRM_BASE + KTRM_INF·exp(−ANC/0.35)·(MUC/4) + KTRM_REN·g_ren + KTRM_CARD·g_card
overall survival = P(cure) × exp(−TRM)
```

---

## 5. The 19 scenarios (summary)

Full text in [`osa_reference_output.txt`](../../../osteosarcoma/osa_reference_output.txt) §1. (The arm that moves IE pre-operatively is run at a resistance of RES₀ = 0.20 and therefore sits separately in §3.E. The R model includes it as a 20th named scenario.)

| Scenario | Huvos% | log10 | cure | TRM% | **Survival** | HR |
|---|---|---|---|---|---|---|
| Surgery alone (historical control) | 1.2 | 0.00 | 0.170 | 0.9 | **0.168** | 3.44 |
| **MAP (EURAMOS-1 arm A)** | 90.3 | 6.61 | 0.602 | 1.0 | **0.596** | 1.00 |
| MAPIE (intensification in poor responders) | 90.3 | 6.93 | 0.659 | 3.0 | **0.639** | 0.86 |
| Urine pH 6.0 | 89.6 | 6.56 | 0.592 | 38.1 | **0.366** | 1.94 |
| Urine pH 7.5 + full hydration | 90.3 | 6.61 | 0.602 | 1.0 | **0.596** | 1.00 |
| Half hydration | 93.3 | 7.02 | 0.674 | 11.7 | **0.595** | 1.00 |
| Urine pH 6.6 | 88.4 | 6.49 | 0.578 | 22.6 | **0.447** | 1.55 |
| Urine pH 6.6 + glucarpidase 48h | 88.2 | 6.47 | 0.575 | 22.4 | **0.446** | 1.56 |
| MAP + dexrazoxane | 90.3 | 6.61 | 0.602 | 1.0 | **0.596** | 1.00 |
| Doxorubicin 600 (unprotected) | 89.6 | 6.91 | 0.655 | 3.2 | **0.635** | 0.88 |
| **Doxorubicin 600 + dexrazoxane** | 89.6 | 6.91 | 0.655 | 1.0 | **0.649** | **0.83** |
| MAP + zoledronic acid (OS2006) | 90.4 | 6.61 | 0.601 | 1.0 | **0.595** | 1.00 |
| MAP + denosumab | 90.4 | 6.61 | 0.601 | 1.0 | **0.595** | 1.00 |
| MAP + mifamurtide (INT-0133) | 90.3 | 6.69 | 0.616 | 1.0 | **0.610** | 0.95 |
| AP alone (MTX omitted) | 61.1 | 4.81 | 0.248 | 0.9 | **0.246** | 2.71 |
| MAP, metastatic at diagnosis | 90.3 | 6.61 | 0.325 | 1.0 | **0.322** | 2.19 |
| MAP, intrinsically resistant tumour (RES₀ 0.25) | 75.2 | 3.81 | 0.176 | 1.0 | **0.174** | 3.38 |
| MAP + regorafenib maintenance | 89.2 | 6.60 | 0.598 | 1.0 | **0.593** | 1.01 |
| MAP, MTX 15 g/m² | 91.8 | 6.80 | 0.635 | 1.0 | **0.629** | 0.90 |

---

## 6. Calibration

| Target | Literature value | Model value |
|---|---|---|
| HDMTX C24 / C48 / C72 | < 10 / < 1 / < 0.1 µM | 8.8 / 0.56 / 0.042 |
| MTX solubility at pH 5 / 6 / 7 | 0.39 / 1.55 / 9.04 mg/mL | identical (log-linear fit) |
| target urine pH | ≥ 7.0 (protocol) | **derived value 7.30** |
| cure rate with surgery alone | 0.15–0.20 | 0.168 |
| MAP 5-year EFS | 0.54–0.59 | 0.596 |
| survival in metastatic disease | 0.20–0.30 | 0.322 |
| cumulative doxorubicin / cisplatin | 450 / 480 mg/m² | 450 / 480 mg/m² |
| LVEF after MAP | 62 → high 50s | 62.0 → 58.5 |
| high-frequency hearing loss after MAP | grade 1–2 | 12.2 dB |
| MAP ANC nadir | grade 4 | 0.39 ×10⁹/L |
| MAPIE EFS HR (poor responders) | 0.98 (0.78–1.23) | **0.98** |
| protocol completion rate MAP / MAPIE | 76% / 51% | 100% / 76% (25-point gap) |
| treatment-related mortality | ~1% | 1.0% (MAP) / 3.0% (MAPIE) |
| prognostic strength of Huvos | strongest prognostic factor | r = 0.97 |
| zoledronic acid EFS | no benefit | HR 1.00 |

---

## 7. How to run

```r
# The model proper (R + mrgsolve)
Rscript osa_mrgsolve_model.R

# Interactive dashboard (8 tabs)
shiny::runApp("osa_shiny_app.R")
```

```bash
# Numerical reference implementation (needs only numpy/scipy; regenerates every result above)
python3 osa_reference_model.py       # -> osa_reference_output.txt + JSON

# Render the mechanistic map
dot -Tsvg osa_qsp_model.dot -o osa_qsp_model.svg
dot -Tpng -Gdpi=150 osa_qsp_model.dot -o osa_qsp_model.png
```

Why `osa_reference_model.py` exists: the build container of this repository has
no R. The Python file uses **the same equations and the same parameter names**,
and guarantees that every number quoted in the README and in the map is a run
result and not a claim.

---

## 8. Where this model could be wrong (limitations)

- `KPPT`, `KDIS`, `PPT50`, `KINJ_PPT`, `KINJ_TUB` are **assumed parameters with
  no direct measurement**. That the derived critical urine pH falls near the
  guideline value is a *constraint* on those assumptions, not an independent
  validation. The **shape** of the result, however — the existence of a
  bifurcation point and the 4.81-fold pH/hydration exchange rate — comes out of
  the solubility law itself and was not fitted.
- A single scalar `RES` represents cross-resistance across all four cytotoxics.
  This is deliberately conservative: it makes intensification look **worse**
  than a per-agent resistance model would, and it is part of why the MAPIE
  prediction matches the trial result.
- `LUNG_IMM = 4` (the weighting of immune kill in the lung) is an assumption,
  and it is the only route through which mifamurtide acts. The predicted
  mifamurtide benefit (+1.4 %p) is therefore a *consequence* of that assumption
  and *not evidence for it*.
- `DEXR` and `MIFA` are not plasma drug but **effect-compartment surrogates**
  (iron chelation / TOP2B depletion and macrophage activation persist far longer
  than the parent molecule).
- The absolute TRM of the MAP arm (1.0%) matches observation, but the absolute
  level of schedule adherence (100% vs the reported 76%) does not. The model has
  no non-haematological or social reasons for discontinuation. **What is
  reproduced is the 25-point gap between the two arms, and that is the part the
  argument needs.**
- Local therapy (surgical margins, radiotherapy) is simplified to a single event
  at week 11 that removes the primary lesion. Local recurrence dynamics are not
  modelled.

---

## 9. Disclaimer

This is a qualitative to semi-quantitative QSP model for teaching and research
purposes. It was assembled from the public literature and clinical trial data
but has not been independently verified or certified, and **must not be used for
real clinical decision-making, prescribing, or regulatory submission.**
