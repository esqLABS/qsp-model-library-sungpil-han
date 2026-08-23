# Radiation-induced lung injury (RILI) — QSP model

**Radiation-Induced Lung Injury · Radiation Pneumonitis & Pulmonary Fibrosis**

> A radiotherapy plan is **not a dose but a distribution**. This model divides the
> lung into six bins of the dose-volume histogram (DVH), puts 10 state variables in
> each bin, and then writes the injury as **two loops with different time
> constants** — a subcritical inflammatory loop on a scale of weeks (self-resolving,
> steroid-responsive) and a **bistable** fibrotic loop on a scale of months (once
> crossed it does not come back, steroid-unresponsive). Out of these two separations,
> a substantial part of the radiation lung-toxicity literature becomes **arithmetic
> rather than assertion**.

| | |
|---|---|
| **Mechanistic map** | [`rili_qsp_model.svg`](../../../radiation-induced-lung-injury/rili_qsp_model.svg) · [`rili_qsp_model.dot`](../../../radiation-induced-lung-injury/rili_qsp_model.dot) — 244 nodes / 385 edges / 18 clusters |
| **mrgsolve model** | [`rili_mrgsolve_model.R`](rili_mrgsolve_model.R) — **79 ODEs** (6 DVH bins × 10 states + 19 global), 40 scenarios |
| **Shiny app** | [`rili_shiny_app.R`](../../../radiation-induced-lung-injury/rili_shiny_app.R) — 10 tabs |
| **References** | [`rili_references.md`](../../../radiation-induced-lung-injury/rili_references.md) — **236 papers**, every PMID verified by lookup |
| **Reference implementation** | [`rili_reference_model.py`](../../../radiation-induced-lung-injury/rili_reference_model.py) · [output](../../../radiation-induced-lung-injury/rili_reference_output.txt) — dependency-free Python RK4 |
| **Calibration · derived-quantity analysis** | [`rili_calibration.py`](../../../radiation-induced-lung-injury/rili_calibration.py) · [output](../../../radiation-induced-lung-injury/rili_calibration_output.txt) |

<a href="../../../radiation-induced-lung-injury/rili_qsp_model.svg"><img src="../../../radiation-induced-lung-injury/rili_qsp_model.png" width="820" alt="RILI QSP map"></a>

---

## 1. What is written differently here

### ① The plan is a vector — mean lung dose is only one moment of the histogram

If the lung is taken as a single compartment and described by the single number of
mean lung dose (MLD), then two plans with the same MLD are **by definition the same
patient**. In this model the plan is `(D_RX, NFX, v₁…v₆)`, every state variable
exists per bin, and the organ-level outputs (DLCO · FVC · CTCAE grade · fibrosis
score) are volume-weighted integrals over the bins. In each bin the dose passes
through a linear-quadratic (LQ) transformation together with the fraction size:

```
BED_b = D_b · (1 + d_b/(α/β)),    d_b = D_b/NFX,    RB_b = BED_b / T_course
```

Because lung has α/β = 3 Gy and tumour 10 Gy, **the fractionation scheme moves the
two tissues by different factors**. This is the only point at which fractionation
enters this file.

### ② The latent period is created by cell turnover, not by dose

Type II alveolar cells (AT2) that have taken lethal radiation damage **do not die
immediately.** They die when they next attempt to divide. So the model uses the
three stages `AT2 → DOOM → death`, and surfactant production is proportional to
**(AT2 + DOOM), not to surviving AT2**. That is to say a lethally damaged cell keeps
working until it dies, and as a consequence the alveolar lining keeps thinning for
weeks after the last fraction has finished. This is precisely why late-responding
tissue injury is **delayed**.

### ③ The two loops

**The fast loop** (weeks): death flux → DAMP → NF-κB → cytokines → oedema. It is
designed to be **subcritical**, satisfying a stability condition:

```
GIMM·(KAMP + KDAMP·KDCYT/KMCYT) + KSYSIN·KSOUT/KCLS  <  KCE
no drug   0.0717 < 0.150        durvalumab   0.0910 < 0.150
```

So pneumonitis resolves without treatment. Under the initial parameters this gain
was 0.42, exceeding 0.15, and as a result any plan at all — down to SBRT with an MLD
of 4.3 Gy — set the whole lung alight and ended in grade 4 ([DEFECT 3]).

**The slow loop** (months): TGF-β1 → myofibroblast → collagen → stiffness → TGF-β1.
It is **bistable**. The three fixed points follow algebraically from the parameters
and the simulation reproduces them:

| Fixed point | Active TGF-β1 | Myofibroblasts | Total collagen | Stability |
|---|---|---|---|---|
| healthy | 0.0401 | 0.00057 | **1.0000** | stable |
| separatrix | 0.2115 | 0.3563 | **1.5864** | **unstable** |
| fibrotic | 0.7767 | 0.9692 | **2.9400** | stable |

Bistability holds because the mechanotransduction term is written as **a threshold,
not a proportion**. Latent TGF-β1 is released by force transmitted through
integrins, and that force is only transmitted once the matrix has stiffened enough
to resist myofibroblast contraction. So the term is
`KTGFM·CX³/(KMCOL³+CX³)`. Written instead as `km·CX/(1+CX)` it becomes almost linear
at low collagen, and the loop gain **destroys** the healthy fixed point — the
unirradiated lung has nowhere to sit ([DEFECT 2]).

### ④ The two loops respond to different drugs and appear in different endpoints

This is the central clinical claim of the model. **A drug that reaches one loop
cannot be detected by an endpoint that reads the other loop.**

---

## 2. Results

Every number is output from `python3 rili_reference_model.py`, and the whole of it
is in [`rili_reference_output.txt`](../../../radiation-induced-lung-injury/rili_reference_output.txt).

### 2.1 Calibration — the QUANTEC dose-response

There are **only two** scale constants not identified by mechanism (`PNI50`,
`PNISL`) and they were fitted to the QUANTEC dose-response. The rest are either
rate constants from the literature or values fixed algebraically.

| Mean lung dose | Model | QUANTEC target |
|---|---|---|
| 13 Gy | 8.9 % | 10 % |
| 20 Gy | 21.2 % | 20 % |
| 24 Gy | 29.3 % | 30 % |

Natural history by plan (no drug):

| Plan | MLD | V20 | V40 | NTCP grade≥2 | Time of peak | Volume converted to fibrosis |
|---|---|---|---|---|---|---|
| SBRT 54/3 | 4.30 | 4.4 % | 1.1 % | **2.1 %** | 58.8 d | 6.5 % |
| hypofractionated 60/8 | 8.64 | 13.7 % | 5.4 % | 6.3 % | 56.8 d | 10.0 % |
| proton 60/30 | 12.22 | 21.8 % | 9.8 % | 8.0 % | 39.0 d | 4.5 % |
| IMRT 60/30 | 17.75 | 34.4 % | 16.7 % | 16.9 % | 40.0 d | 8.0 % |
| 3D-CRT 60/30 | 16.72 | 31.3 % | 18.7 % | 14.4 % | 40.0 d | 10.0 % |
| IMRT 74/37 | 24.94 | 44.0 % | 27.5 % | **26.7 %** | 36.2 d | 21.0 % |

Every peak falls 36–59 days (5.2–8.4 weeks) **after** the last fraction, inside the
observed 4–12 week window.

### 2.2 Dose sets the magnitude, turnover sets the timing

Shaking MLD 6-fold from 4.7 to 28.0 Gy:

| MLD (Gy) | 4.70 | 7.00 | 11.61 | 17.75 | 24.65 | 27.98 |
|---|---|---|---|---|---|---|
| Peak pneumonitis index | 0.192 | 0.277 | 0.447 | 0.671 | 0.922 | **1.043** |
| Time of peak (days post-RT) | 39 | 39 | 40 | 40 | 40 | **40** |

**The magnitude moves 5.4-fold while the timing moves one day.** This is how the
claim that the latent period is a property of AT2 turnover rather than of dose
becomes testable. (In the initial implementation the peak sat 12 days after the end
of RT and shaking `KMIT` 3.3-fold moved it only from 13 to 8 days — because `KMIT`
was not the governing constant. The endothelial and epithelial leak sources had been
bundled into a single coefficient, so the two terms were cancelling each other
([DEFECT 9]).)

### 2.3 ★ Same mean dose, opposite ranking

Two plans matched on MLD with **only the shape** of the histogram changed:

| | low-dose bath | focal hot spot |
|---|---|---|
| Mean lung dose | 15.53 Gy | 15.35 Gy |
| V5 | **92.9 %** | 39.8 % |
| V40 | 4.9 % | **18.8 %** |
| **Pneumonitis NTCP** | **15.1 %** | 11.7 % |
| **Volume converted to fibrosis** | 2.0 % | **12.0 %** |
| DLCO nadir | 31.4 | 45.1 |

**The two endpoints rank the two plans in opposite orders.** Pneumonitis is a
volume-weighted *sum* over a parallel organ, so it is worse in the plan with
low dose spread widely; fibrosis is a local *threshold*, so it is 6-fold worse in
the plan with the large high-dose volume. With mean lung dose alone this distinction
is impossible in principle.

### 2.4 Fibrosis is a threshold, not a curve

Uniform whole-lung irradiation, 2 Gy/fraction, collagen excess at day 730:

| Total dose | 36 Gy | 40 Gy | **42 Gy** | 44 Gy | 60 Gy |
|---|---|---|---|---|---|
| Collagen excess | 0.053 | **0.091** | **2.851** | 2.851 | 2.851 |

It **switches completely** between 40 Gy and 42 Gy. A 5 % increase in dose moves the
response from 0.09 to 2.85. This threshold is not a fitted value; it follows from the
critical activation of the slow loop, 0.0207, and from `KTGFR`, and it becomes a
mechanistic explanation for why radiation fibrosis in the clinic is **sharply
delimited by isodose lines** — whereas pneumonitis is not.

### 2.5 Steroids bend the acute curve and leave the late curve untouched

Prednisolone 60 mg/d, started on day 56:

| | +0 d | +14 d | +28 d | +56 d | Volume converted to fibrosis |
|---|---|---|---|---|---|
| no treatment | 0.6004 | 0.6586 | 0.6711 | 0.6123 | **8.0 %** |
| prednisolone 30 mg/d | 0.6004 | 0.3672 | 0.2863 | 0.2266 | **8.0 %** |
| prednisolone 60 mg/d | 0.6004 | 0.3649 | 0.2849 | **0.2246** | **8.0 %** |

After 56 days the index is 63 % lower but **the volume converted to fibrosis is
identical to the first decimal place.** That 30 mg and 60 mg are indistinguishable
is also a prediction — both are already on the plateau at about 125 times the IC50.
Started **after** the peak (day 120), neither the peak nor the NTCP moves at all
(0.6715 / 16.9 %, identical to no treatment).

### 2.6 Antifibrotics cannot be detected by a pneumonitis endpoint

| | Peak pneumonitis index | NTCP | Volume converted to fibrosis |
|---|---|---|---|
| no drug | 0.6715 | 16.9 % | 8.0 % |
| pirfenidone 0–365 d | **0.6715** | **16.9 %** | **0.0 %** |
| nintedanib 0–365 d | **0.6715** | **16.9 %** | **0.0 %** |
| lisinopril 0–730 d | **0.6715** | **16.9 %** | **0.0 %** |

All three drugs completely prevent fibrosis while **not changing the pneumonitis
index to the fourth decimal place.** A trial with pneumonitis as its primary
endpoint cannot detect any of these drugs. Conversely, a trial reading DLCO at 12
months cannot detect steroids.

### 2.7 ★ The window for antifibrotics closes within a week of the end of RT

Changing only the pirfenidone start day (RT ends on day 42, observed to day 900):

| Start day | 0 | 14 | 28 | **42** | **49** | 56 | 70 | none |
|---|---|---|---|---|---|---|---|---|
| Fibrosis score (900 d) | 0.007 | 0.007 | 0.008 | **0.013** | **0.235** | 0.235 | 0.235 | 0.236 |
| Volume converted to fibrosis | 0 % | 0 % | 0 % | **0 %** | **8.0 %** | 8.0 % | 8.0 % | 8.0 % |

**The window closes between day 42 and day 49, and it closes completely rather than
gradually.** Because the slow loop is bistable and LOX-crosslinked collagen is not a
substrate for MMP, an antifibrotic is not a dose-response drug but **a race against
the separatrix**. Given concurrently it prevents fibrosis completely; given two
months later it does nothing. This is a testable prediction, and it is the exact
opposite of the way such drugs are usually trialled (started after fibrosis has been
diagnosed).

In the implementation without a crosslinked pool, pirfenidone started a year after
the switch had flipped abolished fibrosis entirely — that is, the model claimed to
**reverse established radiation fibrosis**. No such thing happens ([DEFECT 11]).
Mature matrix is a mechanical **memory**.

### 2.8 Amifostine's effect is a property of the schedule, not of the drug

The plasma half-life of WR-1065 is 8 minutes. The protection factor is
`PFAMI·exp(−ln2·Δt/8min)`:

| Interval before irradiation | 0 min | 15 min | 30 min | 60 min | control |
|---|---|---|---|---|---|
| NTCP grade≥2 | **3.3 %** | 13.1 % | 15.9 % | 16.8 % | 16.9 % |

Delayed by 30–60 minutes the effect essentially disappears. Avasopasem (half-life 27
minutes) retains 12.5 % even at a 30-minute interval. This is a mechanistic
hypothesis for the conflicting results of the amifostine trials, and **it is stated
as a model-level hypothesis, not an observation.**

### 2.9 Durvalumab raises the gain rather than adding a term — so the excess risk tracks MLD

| | chemoRT alone | + durvalumab | Absolute increase |
|---|---|---|---|
| MLD 12.2 Gy | 8.0 % | 9.2 % | **+1.2 points** |
| MLD 17.8 Gy | 16.9 % | 19.2 % | **+2.3 points** |
| MLD 24.9 Gy | 26.7 % | 30.2 % | **+3.5 points** |

Writing checkpoint blockade as **the gain GIMM** of the fast loop (0.0717 → 0.0910,
effective time constant 12.8 → 17.0 days) makes the excess risk grow with MLD rather
than being a constant offset. PACIFIC reported any-grade 24.8 → 33.9 % and grade 3–4
2.6 → 3.4 %, and the model's grade≥2 excess falls between them. The prediction that
goes beyond the trial is that **the interaction is not spread evenly but concentrated
in the high-MLD subgroup**.

### 2.10 Same injury, different grade

The same plan, differing only in baseline reserve:

| Baseline DLCO | 100 % | 85 % | 60 % | 45 % |
|---|---|---|---|---|
| NTCP grade≥2 | 14.3 % | 16.9 % | 23.7 % | **30.6 %** |
| DLCO nadir | 44.1 | 37.5 | 26.5 | 19.8 |
| Volume converted to fibrosis | **8.0 %** | **8.0 %** | **8.0 %** | **8.0 %** |

The injury is identical in all four cases and **only the grade** moves. Because grade
is a threshold on injury divided by reserve.

### 2.11 That SBRT wins on both axes at once is geometry, not radiobiology

| Plan | Fraction size | Tumour BED₁₀ | MLD | TCP | NTCP | **UCP** |
|---|---|---|---|---|---|---|
| IMRT 60/30 | 2.0 Gy | 72.0 | 17.75 | 42.5 % | 16.9 % | 35.4 % |
| proton 60/30 | 2.0 Gy | 72.0 | 12.22 | 42.5 % | 8.0 % | 39.2 % |
| IMRT 74/37 | 2.0 Gy | 88.8 | 24.94 | 99.5 % | 26.7 % | 72.9 % |
| hypofractionated 60/8 | 7.5 Gy | 105.0 | 8.64 | 100 % | 6.3 % | 93.7 % |
| SBRT 54/3 | 18.0 Gy | 151.2 | 4.30 | 100 % | 2.1 % | **97.9 %** |

α/β of 3 against 10 moves the ratio **against the lung** as fraction size grows (in
SBRT's highest-dose bin the lung BED₃ is 353.8 Gy₃, 3.7 times IMRT's 95.1). SBRT
nonetheless wins because the **volume** irradiated falls — MLD 17.8 → 4.3 Gy. The two
contributions enter the model at different places, so they can be measured
separately.

### 2.12 The volume-effect exponent the model produces — it disagrees with the literature

**300 random DVHs** were generated at 60 Gy/30 fx, and the regression residual was
used to measure which gEUD exponent `a = 1/n` collapses the model's own risk onto a
single curve:

| a | 0.50 | **0.70** | 1.10 | 1.70 | 2.90 | 5.90 |
|---|---|---|---|---|---|---|
| n = 1/a | 2.00 | **1.43** | 0.91 | 0.59 | 0.34 | 0.17 |
| Residual SD | 0.0315 | **0.0158** | 0.0386 | 0.0828 | 0.1361 | 0.1889 |

The model produces **n = 1.43**. The literature value (Seppenwoolde 2003, NSCLC) is
**0.99**. That is, the lung in this model is **more** of a mean-dose organ than a
"mean-dose organ", so widely spread low dose matters more than the mean would
suggest. This is reported as it stands rather than tuned away. The direction itself
is consistent with the result in §2.3 (bath plans are worse for pneumonitis) and
points the same way as the observation that V5 emerges as a predictor in trimodality
series. The magnitude, however, disagrees with the literature by 44 %.

---

## 3. Predictions that disagree with the literature · limitations

Stated explicitly rather than tuned away.

1. **SBRT pneumonitis is underpredicted.** For peripheral SBRT the model gives
   grade≥2 of 2.1 % (baseline DLCO 85 %) or 3.5 % (inoperable COPD, reserve 55 %).
   Reported values are around 9–10 %. This is a direct consequence of a structure
   that uses MLD as the governing variable, and in real SBRT it appears that the
   steep dose fall-off, central location and pre-existing lung disease are not
   summarised by MLD.
2. **Volume-effect exponent n = 1.43 against the literature's 0.99** (§2.12).
3. **The amifostine time-window interpretation is a hypothesis.** The model shows
   that temporal coincidence can explain the discrepancy between trials, but whether
   that discrepancy was actually due to the dose-to-irradiation interval has not
   been measured.
4. **The 12-month DLCO decline is at the top of the literature range.** For IMRT
   60/30 it gives −19.8 %, where the reported range is roughly −10 to −20 %.
5. **Chemotherapy radiosensitisation is not an explicit mechanism.** Concurrent
   chemotherapy enters only by raising `KDAM`.
6. **Dose delivery has been smoothed.** Fractions are written as a continuous BED
   delivery rate rather than a daily bolus. Total BED and the fractionation
   dependence are exactly preserved, but interfraction sublethal damage repair is
   not represented.
7. **The CTCAE grade cut-points are a convention.** The thresholds turning the
   continuous index into grades were matched to the model's output distribution, not
   to an observed distribution.
8. **Discretisation into six bins.** Increasing the number of bins may change the
   derived value of n.
9. **`PNI50` and `PNISL` are fitted values** — they do not come from mechanism.
10. **Perfusion weighting of the low-dose bins is a user input.** The model does not
    know how perfusion is actually distributed in emphysema.

---

## 4. Validation

Because this was written in an environment with no R runtime, **every equation of the
mrgsolve file was first implemented and run in dependency-free Python (RK4, fixed
step)**, and [`rili_mrgsolve_model.R`](rili_mrgsolve_model.R) is a literal port of that
result. **Twelve real defects** came to light in the process, and each is left as a
`[DEFECT n]` comment at the place it was fixed.

Those that changed the structure:

| # | What was wrong | Symptom |
|---|---|---|
| 1 | Slow-loop gain of 2.74 at the healthy fixed point | Unirradiated lung fibrosed on its own |
| 2 | Mechanotransduction written as a proportion (`km·CX/(1+CX)`) | **Only one** fixed point — the fibrotic state alone exists |
| 3 | The bystander-killing loop was self-sustaining | Everything down to SBRT 4.3 Gy came out grade 4 |
| 4 | The oedema state saturated at ~7 | Every bin above 20 Gy identical — the volume effect cannot be asked about |
| 5 | An artificial saturation was imposed on the index | The volume effect n grew by construction (smuggling in the assumption) |
| 6 | Collagen time constant of 50 days | Conventional fractionation gave no fibrosis and only SBRT fibrosed — the exact opposite of the clinic |
| 7 | Durvalumab Emax of 0.90 | The fast loop exploded, DLCO 0.56 %pred in every group |
| 8 | Tumour repopulation applied unconditionally | 0.85 cells repopulated as e²³, TCP 0 % throughout |
| 9 | Two leak sources bundled into one coefficient | The peak sat 12 days after the end of RT (observed 4–12 weeks) |
| 10 | Surfactant proportional to surviving AT2 only | Recovery the instant the beam stopped — the latent period disappeared |
| 11 | No crosslinked pool | Antifibrotics reversed established fibrosis |
| 12 | Steroids reached only via the cytokine pathway | Symptomatic improvement stalled at 16 % |

The unirradiated control run (beam off) returns to baseline exactly even after 730
days: collagen drifts by 1.2 × 10⁻⁵ and DLCO is 84.998 %pred.

---

## 5. Model structure

### State variables (10 per bin × 6 bins = 60)

| State | Meaning |
|---|---|
`AT2` | Surviving type II alveolar cells
`DOOM` | Lethally damaged AT2 — still functioning, waiting to divide
`EC` | Microvascular endothelial integrity
`SURF` | Surfactant pool, produced from `(AT2+DOOM)`
`PERM` | Alveolar-capillary leak / oedema
`CYT` | Local pro-inflammatory cytokine pool (fast loop)
`TGFB` | Active TGF-β1 (slow loop)
`MFB` | Myofibroblast density
`COL` | Soluble collagen
`XCOL` | Crosslinked collagen — mechanical memory, not an MMP substrate

### Global (19)

`CYTS` (systemic cytokines — the only route to out-of-field pneumonitis) ·
`BEDL₁…₆` (cumulative lung BED₃) ·
`BEDT` · `TUMLN` · prednisolone (2) · pirfenidone (2) · nintedanib (2) ·
durvalumab (2) · lisinopril (2)

### Drugs and the loop each reaches

| Drug | Target | Loop reached |
|---|---|---|
prednisolone | GR → NF-κB transrepression + ENaC fluid clearance | **fast loop only**
amifostine / avasopasem | thiol radical scavenging / SOD mimetic | **initial damage only** (requires temporal coincidence)
pirfenidone | TGF-β signalling + collagen synthesis | **slow loop only**
nintedanib | PDGFR/FGFR/VEGFR | **slow loop only**
lisinopril | latent TGF-β activation | **slow loop only**
durvalumab | PD-L1 → fast-loop gain ↑ | **fast loop (unfavourably)**

---

## 6. Running

```bash
# reference implementation (runs without R, no dependencies)
python3 rili_reference_model.py       # 40 scenarios
python3 rili_calibration.py           # calibration + derived volume-effect exponent

# rendering the mechanistic map
dot -Tsvg rili_qsp_model.dot -o rili_qsp_model.svg
dot -Tpng -Gdpi=96 rili_qsp_model.dot -o rili_qsp_model.png
```

```r
# mrgsolve
library(mrgsolve); library(dplyr)
mod <- mread("rili_mrgsolve_model", ".")
mod %>% param(DRX = 60, NFX = 30,
              V1 = .40, V2 = .20, V3 = .12,
              V4 = .11, V5 = .09, V6 = .08) %>%
  mrgsim(end = 730, delta = 1) %>% as_tibble()

# Shiny dashboard (10 tabs)
shiny::runApp("rili_shiny_app.R")
```

---

## ⚠️ Disclaimer

This is a **qualitative / semi-quantitative QSP model** for educational and research
purposes. It was constructed from the public literature but has not been
independently validated or certified, and **must not be used for the evaluation of
real radiotherapy plans, for clinical decision-making, for prescribing, or for
regulatory submission.** In particular, the NTCP values and dose thresholds given
here cannot be used to approve a real patient's plan. The parameters and assumptions
are illustrative approximations and the limitations in §3 must be read alongside
them.
