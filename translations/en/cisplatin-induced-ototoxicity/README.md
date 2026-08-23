# Cisplatin-Induced Ototoxicity (CIO) — QSP Model

<a href="cio_qsp_model.svg"><img src="cio_qsp_model.png" width="640" alt="CIO QSP mechanistic map"></a>

> **One-line summary.** Cochlear platinum is not a concentration but a
> **reservoir**, and getting from plasma to the cochlea passes through a
> cascade of **plasma → stria vascularis → perilymph → hair cells**. The
> half-life of free plasma platinum is 22 minutes, yet perilymph platinum
> does not peak until 6.9 hours — so sodium thiosulfate given **6 hours**
> after cisplatin cannot touch 99.996% of the systemic exposure that has
> already occurred, while it is still **upstream of 79%** of the platinum
> the cochlea has yet to take up. The "6-hour rule" that clinical trials
> discovered is not a pharmacological trick but **a single time lag**, and
> this model derives it arithmetically.

| Deliverable | File |
|---|---|
| Mechanistic map (140 nodes · 174 edges · 15 clusters) | [`cio_qsp_model.dot`](cio_qsp_model.dot) · [SVG](cio_qsp_model.svg) · [PNG](cio_qsp_model.png) |
| mrgsolve ODE model (**73 ODEs**, 20 scenarios) | [`cio_mrgsolve_model.R`](cio_mrgsolve_model.R) |
| Shiny dashboard (**10 tabs**) | [`cio_shiny_app.R`](cio_shiny_app.R) |
| Independent Python verification implementation | [`cio_reference_model.py`](../../../cisplatin-induced-ototoxicity/cio_reference_model.py) |
| Verification run output (full) | [`cio_reference_output.txt`](../../../cisplatin-induced-ototoxicity/cio_reference_output.txt) |
| References (**158**, verified directly via PubMed) | [`cio_references.md`](cio_references.md) |

---

## 1. Four Things This Model Writes Differently (The four structural choices)

This model does not fit ototoxicity to a "dose-response curve." It makes
four structural choices and checks whether what is observed clinically
**can be derived** from them.

### (1) Cochlear Platinum Is a Reservoir, Not a Concentration, and That Reservoir Splits in Two

The cochlea has no lymphatic drainage. Inside the cell, platinum splits
into a **reversible, redox-active labile pool** (`PT`, half-life 60 days)
and a **covalently fixed bound pool** (`PTB`, half-life 2 years). What
mass spectrometry still detects decades later is the second reservoir,
yet **only the first one causes damage.**

This was written as a testable claim, not a convenience. Sweeping
`TRET_D` **240-fold** from 1 month to 20 years moves the measured bound
platinum by the same factor, 0.0001 → 0.7804, yet **the audiogram does
not move by a single dB to two decimal places** (verification output
Part 5b). The model's prediction: *residual platinum is a marker of past
exposure, not an agent of ongoing damage, and an intervention that
chelates it years later changes nothing.*

### (2) The Cochlea Is Reached Through a Cascade, Not Directly — and That Is the Treatment Window

`free plasma Pt → stria vascularis → perilymph → hair cells`

| Compartment | Cmax | Tmax |
|---|---|---|
| Free plasma platinum | 13.04 µM | 1.0 h (t½ **22.2 minutes**) |
| Stria vascularis platinum | 0.987 | 1.9 h |
| Perilymph platinum | 0.416 µM | **6.9 h** (t½ 4.9 h) |

Fraction of each integral that **remains** at time *t*:

| t (h) | Free plasma platinum | Cochlear uptake | Tumour adducts (well-perfused) | Tumour adducts (poorly perfused) |
|---:|---:|---:|---:|---:|
| 2 | 6.98 % | 98.4 % | 31.07 % | 53.56 % |
| 4 | 0.166 % | 90.6 % | 2.46 % | 14.54 % |
| **6** | **0.0039 %** | **79.0 %** | **0.15 %** | **3.72 %** |
| 8 | 0.00009 % | 66.7 % | 0.01 % | 0.94 % |
| 12 | ~0 % | 44.4 % | 0.00 % | 0.06 % |

Just read the 6-hour row. Systemic rescue therapy is **downstream** of
99.9961% of the systemic exposure that has already happened and 99.85% of
adduct formation in well-perfused tumour, yet **upstream** of 79% of the
platinum the cochlea has yet to take up. This asymmetry is the entirety
of the SIOPEL-6 protocol.

The reason thiosulfate can bypass this route is also structural: platinum
must cross the stria vascularis through a transporter, while thiosulfate,
being a small hydrophilic anion, goes directly to perilymph. The model
states this asymmetry explicitly as an assumption, and shows that without
it, 6-hour rescue would be impossible.

### (3) Because Glutathione Depletion Saturates, Each Band Has Its Own Critical Platinum Load

Writing glutathione depletion by the oxidant flux as saturating
(`KCON·FLUX·G/(KMG+G)`) yields, in closed form, the **critical labile
platinum load** at which the steady-state balance collapses:

```
PT*_j = KSYN · GMAX_j / (KCON · KROS)
```

The glutathione value at this threshold also falls out in closed form
with no free parameter: `G*² + KMG·G* − GMAX·KMG = 0`.

**Since the two thresholds are computed independently from different
state variables, their agreement is itself the test:**

| Band | PT* | G* | Predicted cycle for PT>PT* | Observed cycle for GSH<G* |
|---|---:|---:|---:|---:|
| 0.25k | 1.195 | 0.1667 | Never crossed | Never crossed |
| 1k | 0.888 | 0.1407 | Never crossed | Never crossed |
| 2k | 0.746 | 0.1272 | 10.06 | 10.10 |
| 4k | 0.620 | 0.1142 | 4.05 | 4.08 |
| 8k | 0.514 | 0.1021 | 2.05 | 2.10 |
| 12.5k | 0.454 | 0.0948 | 2.02 | 2.05 |

Every band agrees to within 0.05 cycles, and even the bands that never
cross agree on **which** band that is by both criteria. So ototoxicity is
not a programmed curve but **a threshold that sweeps toward the apex as
cumulative dose rises**.

### (4) Tonotopic Vulnerability Is the Product of Two Gradients, and Which One Shapes the Audiogram Is Measured

Bands are placed at their actual position on the Greenwood function
`f = 165.4(10^{2.1x} − 0.88)`. Vulnerability = `PT*_j / uptake_j` ∝
(reserve gradient) × (inverse of the uptake gradient): a basal/apical
uptake ratio of 1.98, an apical/basal reserve ratio of 2.63, multiplying
to **5.20**.

At an adult dose of 600 mg/m², the 8 kHz − 1 kHz audiogram gradient =
**41.6 dB**. Removing each gradient in turn:

| | 8 kHz − 1 kHz gradient | Fraction lost |
|---|---:|---:|
| Full model | 41.6 dB | — |
| Reserve gradient removed (`BGSH=0`) | 1.5 dB | **96.3 %** |
| Uptake gradient removed (`BUPT=0`) | 18.0 dB | 56.8 % |
| Both removed | 0.0 dB | 100 % |

**The shape of the audiogram is created by where defence is weaker, not
by where more platinum enters.** The uptake gradient alone cannot produce
the observed high-frequency-sloping audiogram.

---

## 2. Verified Results (Results, all from `cio_reference_output.txt`)

### 2.1 Thiosulfate Delay Sweep — SIOPEL-6 and ACCL0431 Fit into One Picture

Paediatric 480 mg/m² background. Same drug, same dose, same schedule —
only the **delay time** changes.

| Delay (h) | 8 kHz (dB) | Otoprotection | Tumour log-kill loss (well-perfused) | Tumour log-kill loss (poorly perfused) |
|---:|---:|---:|---:|---:|
| 0 | 3.1 | 88.6 % | **64.2 %** | **69.0 %** |
| 1 | 11.7 | 65.6 % | 22.0 % | 35.3 % |
| 2 | 13.4 | 60.8 % | 6.3 % | 17.9 % |
| 4 | 14.0 | 59.2 % | 0.4 % | 14.5 % |
| **6** | **15.4** | **55.6 %** | **0.0 %** | **1.2 %** |
| 8 | 17.3 | 50.2 % | 0.0 % | 0.3 % |
| 24 | 29.8 | 10.5 % | 0.0 % | 0.0 % |

- **55.6% otoprotection at 6 hours** is the same magnitude as SIOPEL-6's
  relative risk reduction (63% → 33%).
- **Simultaneous administration burns 64% of efficacy.** The reason
  thiosulfate must be given late is not pharmacodynamics but this single
  line.
- **In poorly perfused metastatic tumour, the same 6 hours is no longer
  free** (1.2% vs 0.0%, a 24-fold difference). And pulling it forward to
  4 hours jumps that to 14.5% — a testable mechanistic hypothesis for the
  divergent results of ACCL0431 and SIOPEL-6. (Stated explicitly as a
  model-level hypothesis, not an observation.)

### 2.2 Cumulative Dose Is the Exposure Variable, and Grading Is a Staircase on Top of It

| Cycle | Cumulative (mg/m²) | 1k | 2k | 4k | 8k | Brock | SIOP |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 300 | 0.6 | 1.1 | 2.6 | 6.6 | 0 | 0 |
| 4 | 400 | 1.2 | 2.6 | 6.6 | 16.8 | 0 | 1 |
| 5 | 500 | 2.2 | 5.2 | 13.7 | 31.0 | 0 | 1 |
| **6** | **600** | 3.7 | 9.2 | 23.4 | **45.3** | **1** | **2** |
| 8 | 800 | 8.4 | 21.0 | 44.9 | 65.9 | 2 | 3 |
| 10 | 1000 | 15.1 | 35.2 | 61.4 | 76.4 | 2 | 3 |

The clinical intuition that 400–500 mg/m² is the risk inflection point
emerges as an **output, not a parameter**.

### 2.3 Fractionated Dosing Buys Almost Nothing (a prediction that conflicts with intuition)

Same 600 mg/m²:

| Schedule | PTA (dB) | 8 kHz (dB) |
|---|---:|---:|
| 100 mg/m² every 3 weeks × 6 | 9.54 | 45.31 |
| 20 mg/m² for 5 days × 6 | 9.04 | 44.56 |
| 50 mg/m² every 10.5 days × 12 | 8.86 | 43.38 |

This is because perilymph completely smooths out the plasma peak. The
model shows that any benefit from fractionation can only come from
superlinearity in the death function, and that sweeping the Hill
coefficient from 1.0 → 3.0 caps the benefit at **0.50 dB**. If
fractionation is clinically meaningfully beneficial, the death function
must be far steeper than the one this model uses — a falsification
condition stated in testable form.

### 2.4 Post-Treatment Progression Is a Property of the Labile Pool

Untreated follow-up after adult 600 mg/m²:

| Month | PTA (dB) | 4 kHz | 8 kHz | Brock |
|---:|---:|---:|---:|---:|
| 0 | 9.54 | 23.43 | 45.31 | 1 |
| 3 | 13.44 | 32.50 | 57.44 | 1 |
| 6 | 13.69 | 33.08 | 58.13 | 1 |
| 24 | 13.33 | 32.69 | 58.03 | 1 |

Even after treatment ends, 8 kHz progresses a further **+12.7 dB** over
3–6 months before flattening out. Shortening `TLAB_D` to 1 week eliminates
the late progression, and lengthening it to 1 year runs it away to
+40 dB. In other words, **the magnitude of late progression is a direct
measurement of the labile pool's half-life.**

### 2.5 Nephrotoxicity Is Not Parallel Damage but an Amplifier of Ototoxicity

| | Final GFR | PTA (dB) | 8 kHz (dB) | Mean cochlear Pt |
|---|---:|---:|---:|---:|
| With renal function feedback | 79.9 | 9.54 | 45.31 | 0.5231 |
| GFR fixed | 100.0 | 8.91 | 43.62 | 0.5111 |
| Baseline CKD (GFR 45) | 35.6 | 12.08 | 51.34 | 0.5618 |

As GFR falls, free platinum AUC rises, and that raises the cochlear load
in the next cycle. Turning off the loop lowers 8 kHz by 1.7 dB — a small
but directionally and quantitatively stated prediction.

### 2.6 Route: The Round Window Creates a Delivery Gradient That Overlaps the Vulnerability Gradient

Paediatric 480 mg/m². Protection rate by band (%):

| | 0.25k | 1k | 2k | 4k | 8k | 12.5k |
|---|---:|---:|---:|---:|---:|---:|
| Systemic thiosulfate (6 h) | 39.9 | 48.9 | 54.7 | 58.0 | 52.7 | 43.7 |
| Intratympanic thiosulfate (161 µmol) | **20.5** | 34.8 | 45.0 | 53.3 | **53.3** | **48.0** |

The systemic route acts fairly uniformly across the audiogram, whereas
the intratympanic route delivers 20.5% at the apex and 53.3% at the base
— **the most delivery goes where damage is worst**. One cycle of
exposure:

| Route | Plasma AUC (µM·h) | **Tumour AUC** | Perilymph peak |
|---|---:|---:|---:|
| Systemic 20 g/m² | 6970 | **6970.2** | 27.4 µM |
| Intratympanic 161 µmol | 0.0000 | **0.0000** | **1190.4 µM** |

Perilymph concentration is 43-fold higher yet tumour exposure is
unmeasurable. The model's conclusion is that route is the only way to
eliminate the efficacy trade-off in §2.1 altogether.

### 2.7 Risk-Factor Decomposition Is Superadditive

Adult 600 mg/m² background, change in PTA:

| | Δ PTA (dB) | 8 kHz | Brock |
|---|---:|---:|---:|
| Age 1 | +6.78 | 56.6 | 1 |
| Aminoglycoside co-therapy | +9.92 | 62.5 | 2 |
| Noise 85 dBA | +5.82 | 57.1 | 1 |
| Furosemide | +17.67 | 69.3 | 2 |
| Baseline CKD (GFR 45) | +2.54 | 51.3 | 1 |
| **All five together** | **+55.95** | 81.2 | 4 |
| (sum of individual effects) | (+42.73) | | |

Because of the threshold structure, this is **1.31-fold superadditive**.
The clinical impression that risk grows "faster than additively" in
patients with multiple risk factors emerges from the structure.

### 2.8 Grading Systems Are an Output of the Model, Not an Input

Brock (40 dB) · SIOP Boston (20 dB) · ASHA · CTCAE v5.0 are all computed
from the per-band threshold shift. The same continuous trajectory steps
through the staircase at different points depending on the criterion —
showing that disagreement between grading systems is not measurement
error but **the geometry of quantisation**.

---

## 3. Verification

Built in an environment with no R/mrgsolve runtime, so **every equation
was first written and run as an independent Python implementation**
(`cio_reference_model.py`, scipy LSODA, PTA identical to five decimal
places at rtol 1e-6/1e-7/1e-8). The R file is a transcription of it, with
parameter names, state names, and order corresponding 1:1.

In the process, **three real defects** surfaced, each left as a comment
at the point where it was fixed.

1. **The integration driver was truncating state.** When `t_eval` is
   passed to `solve_ivp`, `sol.y[:,-1]` is the state at the **last
   requested output time**, not the end of the interval. State was
   silently truncated at every dosing interval, producing the physically
   impossible result that ototoxicity **decreased** when nephrotoxicity
   parameters were changed. The interval end is now always included in
   `t_eval` and verified with an `assert`.
2. **Intratympanic thiosulfate was applied to a fully mixed perilymph
   pool.** As a result, a single round-window bolus saturated the entire
   cochlea, making 40 µmol and 400 µmol indistinguishable, and erasing
   the very basal bias that is this route's reason for existing.
   Rewritten as a local term.
3. **Residual platinum was written as a single reservoir.** With only a
   single 2-year half-life pool, the death risk stayed at treatment-era
   levels forever, running progression away to **+46 dB** at 2 years. The
   labile/bound split emerged while fixing this defect, and the testable
   claim in §1(1) is a by-product of it.

Numerical soundness: no negative states, `OHC`/`IHC`/`SGN` never leave
[0,1], `EP ∈ [0, EP0]`, glutathione never negative, results unchanged
across 3 orders of magnitude of tolerance.

The R code was checked statically — the 73 `dxdt_` terms correspond
exactly 1:1 to the 73 `$CMT` states, `$MAIN`/`$ODE`/`$TABLE` brackets
balance (14/14, 312/312, 133/133), and the band index `TS7 = 8 kHz`
matches the Brock definition. The Shiny app was confirmed to have all 28
outputs present in both UI and server, and all 23 inputs declared.

---

## 4. Limitations the Model States About Itself (Stated, not tuned away)

1. **Furosemide co-therapy at +17.7 dB** matches the literature in
   direction, but its magnitude has never been validated against human
   data. It is the most exposed prediction in this model.
2. **The weekly low-dose regimen (S20)** predicts GFR 60.7 at a
   cumulative 360 mg/m² (vs 79.9 for 600 mg/m² every 3 weeks). This is
   because the tubular recovery half-life exceeds the 7-day interval, and
   needs separate validation.
3. **The ACCL0431 interpretation** is a model-level hypothesis, not an
   observation.
4. **This is a deterministic "representative patient" model.**
   Incidences such as the trial's 33% vs 63% are a product of
   inter-individual variability, and here pharmacogenomics
   (TPMT · COMT · ACYP2) were left as parameters without a variability
   model attached.
5. The perilymph and stria vascularis compartments are **tracer
   compartments** that do not feed back to plasma (the standard
   assumption that cochlear volume is negligibly small).

---

## 5. Usage

```r
# load the model and run the 20 scenarios
source("cio_mrgsolve_model.R")
print(cio_table(mod))

# thiosulfate delay sweep
sapply(c(0, 2, 4, 6, 8, 12), function(d) {
  o <- as.data.frame(cio_run(param(mod, AGE = 3),
                             cio_ev(6, 80, sts_cycles = 0:5, sts_delay = d)))
  c(delay = d, PTA = tail(o$PTA, 1), TS8k = tail(o$TS7, 1))
})

# dashboard
shiny::runApp("cio_shiny_app.R")
```

```bash
# re-run the independent Python verification (numpy + scipy)
python3 cio_reference_model.py
```

---

## 6. Disclaimer

This is a model for education and research. It must not be used for
clinical decision-making, dose determination, or the selection of
otoprotection protocols. The parameters were set based on the literature
but have not been formally fitted to individual clinical trial data, and
every number here is a model output under the stated assumptions.
