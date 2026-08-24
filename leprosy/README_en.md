# Leprosy (Hansen's Disease) QSP Model

> **Directory:** `leprosy/` | **Abbreviation:** LEP | **Date:** 2026-08-07
> **Category:** Infectious / Mycobacterial / Peripheral nerve

---

[![LEP QSP mechanistic map](lep_qsp_model.png)](lep_qsp_model.svg)

---

## One-sentence summary

Leprosy is tracked and staged by a single number — the **bacterial index (BI)** —
but the rate constant that governs that number is one no drug can touch: this
model separates the clock BI reads from the clock the drugs actually turn, and
shows that reaction · clofazimine loading · the steroid treatment window all
follow arithmetically from that separation.

---

## Structural thesis — two clocks

```
    BI = log10( live bacilli + dead bacilli, per gram )
    MI = 100 × live bacilli / (live bacilli + dead bacilli)      "morphological index"
```

A dead *M. leprae* is still an **acid-fast bacillus**. It disappears only once a
macrophage digests it, and that half-life is about 4–5 months. Hence:

| | What it reads | Who sets it | Timescale |
|---|---|---|---|
| **Clock 1 · killing** | d(live)/dt → **MI** | the drug | hours to weeks |
| **Clock 2 · clearance** | d(dead)/dt → **BI** | the host | t½ ≈ 137 days, **drug-independent** |

Four consequences follow from this separation.

### (1) BI cannot compare regimens

| Comparison | Difference in live-bacillus burden at day 21 | Difference in BI at year 1 |
|------|----------------------|---------------|
| MDT-MB (S02) vs dapsone monotherapy (S06) | **3.83 log10** | 0.24 log10 |
| MDT-MB (S02) vs 1/20 kill rate (S09) | **3.70 log10** | **0.11 log10** |

Even lowering the bactericidal rate constant 20-fold moves BI at one year
by only 0.11 log10 — within smear-reading error. The same two regimens are
**already separated by MI at day 21** (0.0014% vs 6.5%). A slowly falling BI
is neither evidence of treatment failure nor a criterion for stopping treatment.

### (2) Reaction (ENL) is driven by the *derivative* of antigen, not its *level*

Antigen is exposed to the immune system not when a bacillus **dies** but when
it is **digested**, and live bacilli actively suppress phagosome maturation.
So killing a bacterial population **releases, all at once, the breakdown of
the entire dead-bacillus pile that had accumulated up to that moment**:

* Antigen release rate: 1.56 AU/day pre-treatment → **3.48 AU/day (2.2-fold)**
  after MDT starts, then decaying with the dead-bacillus pool at t½ = 137 days
* Cumulative antigen load is **set by the diagnosis-time bacterial burden B₀** —
  changing the kill rate 20-fold shifts the cumulative value at 2,500 days
  by only 15%, and that residual is simply the extra bacterial growth the
  slower arm accrues
* So the drug **can move the ENL peak along the time axis but cannot remove
  its area**: a 20-fold difference in kill rate changes the ENL peak by
  **3.7-fold** (S02 44.8 vs S09 12.1) on the same BI trajectory.

The model reproduces the classical clinical spectrum without being asked to.
Because ENL is the product of `antigen release × antibody`, it appears **only
at the lepromatous (LL) pole** (LL peak 44.8 vs TT ~0), while Type 1 reaction
is the product of `recovered T-cell function × residual nerve antigen`, so it
peaks **at borderline disease** (BL is 4.3-fold LL). Two oppositely directed
gradients, multiplied together.

### (3) Clofazimine's anti-inflammatory action is a *loading* problem

Terminal half-life 70 days → the tissue depot reaches **25% of steady state
at 1 month, 83% at 6 months**. That is, it **arrives later** than the ENL peak
it was meant to blunt. Front-loading 300 mg/day in the first month raises the
day-30 depot **5.3-fold** and cuts cumulative ENL burden **by 49%** (S10).

> ⚠️ This 49% is **the single most exposed prediction in this model**. No
> clinical trial has ever run this comparison.

### (4) Nerve damage: a reversible pool leaks into an irreversible sink

```
    d(NFI_reversible)/dt = damage − recovery(k_rec, steroid-accelerated) − k_fix·NFI_reversible
    d(NFI_permanent)/dt = k_fix·NFI_reversible          k_fix = 1/180 day⁻¹
```

Steroids act **only on the reversible pool**. Two distinct consequences follow.

* **Duration**: the 20-week WHO tapering course leaves **2.4-fold** more
  permanent deficit than the 52-week course (S12 5.3 points vs S13 2.1 points)
  — the same direction in which the Cochrane review could not confirm a
  long-term benefit of the 20-week course
* **Promptness**: even a course of adequate duration is useless if started
  late. **Half of the recoverable deficit disappears within 117 days**, and by
  one year the delay is essentially equivalent to no treatment at all (S14
  10.6 points ≈ 99% of the untreated value of 13.1 points)

---

## By-product: the containment threshold

At low bacterial burden, cell-mediated killing is `KHOST × SPEC`, and the net
growth of a replicating clone is `MUMAX − KNAT`. The model therefore contains
a sharp immune set-point threshold:

```
    SPEC* = (MUMAX − KNAT) / KHOST = (0.05545 − 0.01043) / 2.00 = 0.0225
```

* `SPEC > 0.0225` — the host suppresses the residual clone → no relapse after
  completion of treatment
* `SPEC < 0.0225` — suppression fails → this is where relapse and
  **resistance selection** occur

This threshold is confirmed numerically by bisection search on the ODE system
(agreeing with the algebraic value to within 0.7%). Both the ten-year
expansion of a folP1-mutant clone to 8.2 log10/g under dapsone monotherapy
(S07) and MDT's suppression of that same clone in the same patient sit on
opposite sides of this threshold.

---

## Scenario results (17 scenarios, mrgsolve run values)

| # | Scenario | BI₀ | BI at 1 yr | ΔBI/yr | MI day 21 (%) | log live day 21 | ENL peak | ENL burden | T1R peak | Permanent nerve damage | Relapse index (%) |
|---|----------|-----|--------|--------|-------------|--------------|----------|----------|----------|--------------|-------------|
| S01 | Untreated lepromatous natural history | 5.89 | 5.90 | −0.01 | 19.7 | 5.18 | 4.0 | 71.3 | 0.0 | 2.4 | 100 |
| S02 | **WHO MDT-MB 12 months (reference)** | 5.89 | 5.09 | 0.80 | 0.0014 | 1.00 | 44.8 | 13.5 | 3.0 | 0.5 | 1.77 |
| S03 | WHO MDT-PB 6 months | 2.79 | 1.99 | 0.80 | 0.0002 | −2.96 | 0.0 | 0.0 | 0.5 | 0.0 | 0.00005 |
| S04 | Uniform MDT 6 months (in a multibacillary patient) | 5.89 | 5.09 | 0.80 | 0.0014 | 1.00 | 44.8 | 13.7 | 3.0 | 0.5 | **4.16** |
| S05 | MDT-MB 24 months (pre-1998) | 5.89 | 5.09 | 0.80 | 0.0014 | 1.00 | 44.8 | 13.5 | 3.0 | 0.5 | 0.29 |
| S06 | Dapsone monotherapy (historical) | 5.89 | 5.33 | 0.56 | **8.60** | **4.83** | 19.8 | 38.7 | 2.8 | 1.4 | — |
| S07 | Dapsone monotherapy + folP1 mutation, polar patient | 5.89 | 5.60 | 0.29 | 8.61 | 4.83 | 19.8 | 145.8 | 0.1 | 4.9 | 100 |
| S08 | Single-dose rifampicin contact prophylaxis (SDR-PEP) | −2.51 | −3.31 | 0.80 | 0.0019 | −7.29 | 0.0 | 0.0 | 0.0 | 0.0 | 0.006 |
| S09 | **Experiment: 1/20 kill rate** | 5.89 | 5.20 | **0.69** | **6.49** | **4.70** | **12.1** | 11.8 | 2.9 | 0.5 | 2.70 |
| S10 | **Experiment: 1-month clofazimine loading** | 5.89 | 5.09 | 0.80 | 0.0012 | 0.90 | **40.1** | **6.9** | 3.0 | 0.3 | 1.76 |
| S11 | MDT + ofloxacin/minocycline | 5.89 | 5.09 | 0.80 | 0.0004 | 0.41 | 44.8 | 13.5 | 3.0 | 0.5 | 1.65 |
| S12 | Type 1 reaction · 20-week WHO taper | 4.85 | 4.05 | 0.80 | 0.0003 | −0.76 | 0.0 | 0.0 | 9.5 | **5.3** | 0.004 |
| S13 | Type 1 reaction · 52-week course | 4.85 | 4.05 | 0.80 | 0.0003 | −0.76 | 0.0 | 0.0 | 8.5 | **2.1** | 0.003 |
| S14 | Type 1 reaction · same course, 180-day delay | 4.85 | 4.05 | 0.80 | 0.0003 | −0.76 | 0.0 | 0.0 | 12.9 | **10.6** | 0.0002 |
| S15 | Thalidomide 300 mg for ENL, 12 weeks | 5.89 | 5.09 | 0.80 | 0.0014 | 1.00 | 44.8 | **8.2** | 3.0 | 0.3 | 1.77 |
| S16 | MDT in a G6PD-deficient patient | 5.89 | 5.09 | 0.80 | 0.0014 | 1.00 | 44.8 | 13.5 | 3.0 | 0.5 | 1.77 |
| S17 | Rifampicin resistance · second-line ROM + clofazimine | 5.89 | 5.16 | 0.73 | 4.17 | 4.51 | 14.0 | 15.1 | 0.6 | 0.5 | 1.12 |

S16's safety axis: methemoglobin **18.4%** (normal 6.0%), haemoglobin nadir
**5.6 g/dL** (normal 10.7 g/dL). Clofazimine pigmentation reaches an index of
66 at 12 months and rises to 98 with the 24-month regimen (S05).

**Look at S02, S04, and S05 side by side.** The BI trajectory · MI · ENL are
all identical, and the only thing that diverges is the relapse index (0.29% →
1.77% → 4.16%). What treatment duration buys is not BI but the **residual
live-bacillus count**.

---

## Model structure

**38 ODE compartments**

| Group | Compartments |
|------|------|
| Rifampicin | `RIFG` `RIFC` `ENZ` (autoinduction) |
| Dapsone | `DAPG` `DAPC` `NOH` (hydroxylamine) `METHB` `HB` |
| Clofazimine | `CLOG` `CLO1` `CLO2` (tissue depot, t½ 70 days) |
| Prednisolone · thalidomide · second-line drugs | `PDNG` `PDNC` `THAG` `THAC` `ROMG` `ROMC` |
| Bacterial populations | `BG` (replicating) `BP` (persisters) `BR` (resistant) `BD` (**dead — governs BI**) |
| Antigen | `AG` `AGN` (intraneural) `AGC` (cumulative) |
| Humoral · ENL | `AB` `IC` `TNF` `NEU` |
| Cell-mediated immunity | `CMI` `TREG` |
| Reaction · nerve | `T1R` `NFIR` (reversible) `NFIP` (permanent) |
| Descriptive | `LES` `PIG` `HPA` `ALT` `ENLC` |

**Key equations**

```
k_clr  = KCLR0 · (1 + CLRBOOST·(1 − B_live/(B_live+KSUB)))   ← live bacilli block digestion
release = YB·(death flux) + YD·k_clr·B_dead                    ← YB=0.10, YD=0.90
k_host = KHOST · CMI · KSAT/(KSAT + B_live)                    ← granuloma containment saturates
CMI target = CMIMAX·SPEC/(1 + B_live/BSUP)·(1 − TREGSUP·TREG)  ← anergy is driven by *live bacilli*
ENL    = Hill₄(TNF, neutrophils),  immune complex = antigen × antibody
```

It matters that `CMI` anergy is pinned to **live-bacillus burden, not free
antigen**. Live bacilli collapse within days of effective bactericidal
treatment, but the antigen pool persists for years — which is why Type 1
reaction becomes a **treatment-onset phenomenon**, peaking at day 67 in the
model.

---

## Verification

Two independent checks were performed.

1. **`lep_verify_python.py`** — re-implements all 38 ODEs from scratch in
   Python/scipy and checks them against 54 published anchors. **54/54 pass**
   (`lep_verification_output.txt`).
2. **R ↔ Python cross-validation** — the mrgsolve and Python implementations
   agree to 4 significant figures over 361 days (e.g. `BLIVE` at day 359:
   2.859e-04 vs 2.859e-04).

Example anchors: rifampicin Cmax 7.9 mg/L·t½ 4.1 h, dapsone Css 1.63 mg/L·
methemoglobin 5.9%, clofazimine plasma Css 0.58 mg/L·tissue depot 25% at day
30, untreated LL BI 5.9·MI 19%, MDT year-1 BI decline 0.80 log10, dapsone
monotherapy reaching MI 0 at day 181, 4.0 log10 kill from a single 600 mg dose
of rifampicin, relapse index 1.8% after 12 months of MDT-MB.

### Defects found and fixed during tuning

* **An early version pinned anergy to free antigen** — since the antigen pool
  decays over years together with the dead bacilli, cell-mediated immunity
  rose gradually over 332 days, and Type 1 reaction became a flat plateau
  rather than an acute episode, which made the steroid-timing experiments
  meaningless. Switching to live-bacillus burden moved the peak to day 67.
* **A version with no threshold on nerve damage** — even the quiet tail at
  4/100 reaction activity kept eroding nerve continuously, so a 20-week
  steroid course made no difference across a 900-day integral. Adding a Hill
  threshold (T1R50 = 12, h = 5) let acute episodes dominate the damage.
* **A version with granuloma containment as a linear term independent of
  bacterial burden** — set strongly enough to prevent relapse in treated LL
  patients, it spontaneously cured untreated borderline patients; set weakly,
  every LL patient relapsed within 3 years. A saturating term
  `KSAT/(KSAT+B_live)` satisfied both requirements at once, and yielded the
  containment threshold SPEC* as a by-product.
* **A version in which persisters were killed too easily by host immunity** —
  the relapse hazard ratio between 6-month and 12-month regimens came out at
  19-fold (observed ~2-fold). Adding immune evasion for dormant persisters
  (PIMM = 0.02) brought it to 2.4-fold.

### What remains unvalidated

1. The disinhibition coefficient for dead-bacillus breakdown, `CLRBOOST = 2.0`
   — back-calculated from the observed rate of BI decline and never directly
   measured.
2. The 49% reduction in ENL from clofazimine loading — **an untested
   prediction**.
3. The containment threshold `SPEC* = 0.0225` — an internally derived value
   with no measured clinical counterpart.
4. Relapse is not read off the deterministic trajectory as a probability but
   computed as a Poisson escape probability from the **residual live-bacillus
   count** at the end of treatment, because a continuous ODE cannot represent
   the stochastic extinction of a population of a few thousand organisms.
5. The nerve-damage index (0–100) was shaped to match the WHO disability
   grade, not fitted to individual nerve-conduction studies.

---

## Files

| File | Contents |
|------|------|
| [`lep_qsp_model.dot`](lep_qsp_model.dot) | Mechanistic map source (123 nodes · 17 clusters) |
| [`lep_qsp_model.svg`](lep_qsp_model.svg) / [`.png`](lep_qsp_model.png) | Rendered output |
| [`lep_mrgsolve_model_en.R`](lep_mrgsolve_model_en.R) | 38-ODE mrgsolve model + 17 scenarios |
| [`lep_shiny_app_en.R`](lep_shiny_app_en.R) | 10-tab Shiny dashboard |
| [`lep_verify_python.py`](lep_verify_python.py) | Independent Python/scipy re-implementation + 54 anchors |
| [`lep_verification_output.txt`](lep_verification_output.txt) | Verification run output (54/54 pass) |
| [`lep_references_en.md`](lep_references_en.md) | 87 references (all PMIDs confirmed against PubMed) |

### Running it

```bash
# Mechanistic map
dot -Tsvg lep_qsp_model.dot -o lep_qsp_model.svg
dot -Tpng -Gdpi=150 lep_qsp_model.dot -o lep_qsp_model.png

# 17-scenario summary table
Rscript lep_mrgsolve_model_en.R

# Independent verification
python3 lep_verify_python.py

# Dashboard
Rscript -e 'shiny::runApp("lep_shiny_app_en.R")'
```

---

## ⚠️ Disclaimer

This model is a **QSP model for educational · research purposes**. It was
built from published literature but has not been independently validated ·
certified, and **must not be used for actual clinical decision-making,
prescribing, or regulatory submission.**
