# Lipodystrophy Syndromes — QSP Model

**Congenital generalised (CGL) · acquired generalised (AGL) · familial partial (FPLD) · acquired partial (APL)**

| Deliverable | File |
|---|---|
| Mechanistic map (192 nodes · 290 edges · 22 clusters) | [`lipo_qsp_model.dot`](../../../lipodystrophy-syndromes/lipo_qsp_model.dot) · [SVG](../../../lipodystrophy-syndromes/lipo_qsp_model.svg) · [PNG](../../../lipodystrophy-syndromes/lipo_qsp_model.png) |
| mrgsolve model (53 ODEs · 28 scenarios · 7 falsification tests) | [`lipo_mrgsolve_model.R`](../../../lipodystrophy-syndromes/lipo_mrgsolve_model.R) |
| Shiny dashboard (12 tabs) | [`lipo_shiny_app.R`](../../../lipodystrophy-syndromes/lipo_shiny_app.R) |
| Validation run output (actual execution log from mrgsolve 2.0.1) | [`lipo_validation_output.txt`](../../../lipodystrophy-syndromes/lipo_validation_output.txt) |
| References (121 PMIDs, all fully verified) | [`lipo_references.md`](../../../lipodystrophy-syndromes/lipo_references.md) |

---

## The Organising Thesis

Adipose tissue performs two unrelated jobs within a single tissue.

- **(a) A mechanical buffer** — the only compartment that can esterify and store large triglyceride influxes at low thermodynamic cost
- **(b) An endocrine organ** — the source of the leptin signal that reports the buffer's size to the hypothalamus

Lipodystrophy is **a single lesion** (loss of adipocyte mass/function), but its consequence is **two mutually independent deficits**.

```
Capacity deficit :  J_ov = max(0, J_in − ΣJ_st,i − OX)
Signal deficit   :  DEF  = 1 − S(L),   S(L) = L²/(L² + 4²)
```

And the signal deficit **raises J_in through overeating.** In other words the two deficits do not add — they **multiply.**
Unlike the other models in this repository, drugs here are not classified by "indication" but by **which of the two deficit terms they enter.**

| Drug | Term Affected |
|---|---|
| Metreleptin | Signal S → (and, as a consequence) J_in |
| Pioglitazone (PPARγ) | Capacity C_cap — but only while preadipocyte substrate remains |
| Volanesorsen (APOC3) · evinacumab (ANGPTL3) · fenofibrate · omega-3 | **Neither** — removes only the plasma TG that is a *byproduct* of the overflow |
| Insulin · metformin · SGLT2i | Glucose only |
| Very-low-fat diet | J_in directly |

### The GL/PL dissociation is not written into this file — it is derived

This model has **no** parameter saying "metreleptin works better in the generalised form."
There is only a single saturating transfer function `S(L)`, applied directly to total effective leptin.
Because S saturates, **the same dose produces the same receptor occupancy in every patient, but the
achievable change in S is capped at 1 − S(L_baseline)** — that is, capped by **a property of the patient.**

Actual run output (`inference1_dose_vs_deficit()`):

| Scenario | Baseline leptin | Deficit | **24h occupancy** | 24h signal | Time above EC50 | ΔHbA1c | ΔTG | ΔHFF |
|---|---|---|---|---|---|---|---|---|
| CGL + 0.06 mg/kg | 0.64 | 0.975 | **0.695** | 0.867 | 100% | **−2.40** | −70.8% | −86.0% |
| FPLD, low leptin + 0.06 | 5.40 | 0.354 | **0.761** | 0.895 | 100% | −0.37 | −15.9% | −17.3% |
| FPLD, normal leptin + 0.06 | 7.36 | 0.228 | **0.782** | 0.899 | 100% | −0.20 | −5.2% | −5.6% |
| FPLD, normal leptin + **0.13** (dose increase) | 7.36 | 0.228 | **0.855** | 0.906 | 100% | **−0.22** | −5.4% | −5.8% |

Occupancy climbs down the table while the effect shrinks **12-fold.** Target binding is 100% of the day in all four arms.
Doubling the dose to raise occupancy from 0.78 to 0.86 barely moves HbA1c, from −0.20 to −0.22.
**The constraint is the deficit, not the dose.**

---

## Four Structural Commitments

1. **One lesion, two deficits** — see the table above. Drugs are classified by which term they enter.
2. **Plasma TG is only stage 5 of the cascade, not the disease** — overflow is distributed in a fixed
   order: liver → muscle → plasma → pancreas → heart → kidney. So a drug that removes only plasma TG
   should barely move hepatic fat or HbA1c. This is not an assumption but a **prediction**, and
   scenarios 15/16 are its test.
3. **Only fibrosis is irreversible** — every other state reverts once its driver is removed (scenario
   10). Because FIB's regression rate is 1/67 of its formation rate, starting treatment early changes
   the *ceiling* of recovery, not its *speed*.
4. **Nowhere does the model assert positive-feedback bistability** — a withdrawal scenario was
   deliberately included so this claim could be checked. The steatosis → insulin resistance →
   lipogenesis loop exists, but its gain is below 1 (during development, setting this loop's gain too
   high produced 68% hepatic fat, so the coefficient was lowered).

---

## Phenotype = A Capacity Distribution, Not a Label

Phenotype is **not a categorical switch.** Only the capacity ceilings of each compartment (FC1/FC2/FCV),
expansion headroom (PREAD), and the secretion coefficient (LEPKM) are specified; leptin, hepatic fat,
TG, and HbA1c are **outputs obtained by running a 12-year burn-in.** Writing in the initial values by
hand would turn the phenotype from something *derived* into something merely *asserted*.

Run output (`baseline_table()`, each phenotype run to its own steady state):

| Phenotype | Fat (kg) | Capacity (kg) | Leptin | Deficit | Intake | J_in | J_ov | HFF% | TG | Glucose | HbA1c | Insulin | IR | ALT | FIB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Healthy control | 14.38 | 20.0 | 11.31 | 0.118 | 1895 | 157 | 27 | 4.5 | 116 | 99 | 5.43 | 10.9 | 1.41 | 22 | 0.00 |
| **CGL** (AGPAT2/BSCL2) | 0.62 | 0.64 | 0.64 | 0.975 | 2762 | 121 | 74 | **24.2** | **1060** | 185 | **8.73** | 78.8 | 8.64 | 90 | 3.08 |
| AGL (Lawrence) | 1.32 | 1.45 | 1.12 | 0.928 | 2608 | 115 | 68 | 14.5 | 850 | 169 | 8.10 | 53.4 | 6.32 | 55 | 2.36 |
| FPLD2, low leptin | 8.62 | 11.04 | 5.40 | 0.354 | 2087 | 121 | 39 | 5.9 | 234 | 119 | 6.17 | 22.3 | 2.80 | 27 | 0.00 |
| FPLD2, normal leptin | 10.97 | 14.54 | 7.36 | 0.228 | 1985 | 132 | 33 | 5.4 | 174 | 112 | 5.93 | 17.4 | 2.25 | 25 | 0.00 |
| APL (Barraquer-Simons) | 9.74 | 13.35 | 7.83 | 0.207 | 1968 | 129 | 33 | 3.6 | 139 | 106 | 5.67 | 13.5 | 1.76 | 20 | 0.00 |
| PLIN1 FPLD4 (elevated basal lipolysis) | 7.08 | 10.25 | 4.66 | 0.424 | 2143 | 133 | 43 | 7.9 | 336 | 128 | 6.53 | 28.5 | 3.52 | 33 | 0.08 |
| **Congenital leptin deficiency** (control) | 28.63 | 32.0 | 0.63 | 0.976 | 2591 | 186 | 50 | 9.0 | 153 | 141 | 7.04 | 38.3 | 4.59 | 37 | 0.65 |
| **Simple obesity** (control) | 29.84 | 36.6 | 23.94 | 0.095 | 2294 | 215 | 30 | 9.7 | 138 | 117 | 6.09 | 20.5 | 2.61 | 39 | 1.05 |

Two **control phenotypes** matter as much as the disease phenotypes.

- **Congenital leptin deficiency** — fat mass is *high* and signal is zero → a pure signal-deficit arm
- **Simple obesity** — fat mass is high and signal is also *high* (the transfer function saturates from above, plus the SOCS3 brake) → a pure capacity-excess arm

These two controls are what make the claim that "the two deficits are genuinely separable" testable.
And it is this same saturation, approached **from above**, that explains why leptin does not work in obesity.

---

## Seven Falsification Tests — Actual Run Output

| # | Claim | Test Method | Result |
|---|---|---|---|
| 1 | Effect = occupancy × deficit | Same mg/kg given to two phenotypes, plus a dose increase in high-leptin PL | **Confirmed.** Occupancy rises from 0.695 to 0.855 while ΔHbA1c goes from −2.40 to −0.22 |
| 2 | Leptin has an intake-independent component | Intake fixed at the untreated value (pair-feeding) | **Confirmed, but only partially.** **57%** of the hepatic-fat effect and **35%** of the HbA1c effect survive (ΔHFF −48.6%, ΔHbA1c −0.83, ΔTG −19%) |
| 3 | Plasma TG is a stage-5 marker | Removing only TG with an APOC3 ASO | **Confirmed.** TG −80%, hepatic fat −0.0% (HFF/TG ratio 0.000 vs. 1.185 for metreleptin) |
| 4 | PPARγ requires a preadipocyte substrate | **Identical exposure** (1290 ng/mL) given to CGL and FPLD | **Confirmed.** Capacity 0.64→0.64 vs. 11.0→12.1; ΔHbA1c −0.02 vs. −0.26 |
| 5 | A diet that cuts J_in by the same amount should reproduce the drug's effect on hepatic fat | Dietary fat fraction **titrated (uniroot)** to match J_in exactly | **Confirmed.** J_in 84.8 vs. 84.7 → ΔHFF −86% vs. −81%, but ΔHbA1c −2.40 vs. −1.53. The difference is exactly the leptin signal (drug arm S=0.502 vs. diet arm S=0.018 — the diet arm loses more fat and so **loses** signal) |
| 6 | Because only fibrosis is irreversible, the start time changes the ceiling | Starting in year 1 vs. year 10, followed for 20 years | **Confirmed.** At year 20, hepatic fat and HbA1c are identical (3.4%, 5.84), but FIB is 1.11 vs. 1.89 and eGFR is 65.9 vs. 44.6 |
| 7 | Neutralising antibodies eliminate the effect while preserving exposure | NEUT switch, contrasted with the withdrawal arm | **Confirmed.** At year 9 the **measured** leptin is 4.2 (indistinguishable from the normal arm's 4.0), yet effective free leptin is 0.79 vs. 4.03 and HbA1c is 8.13 vs. 5.84. It looks exactly like an adherence problem |

### Pancreatitis Is a Matter of *Time Above Threshold*, Not the Mean

| Scenario | Days with TG>1000 since start | Cumulative risk |
|---|---|---|
| CGL, untreated | 3647 (the whole period) | 4.600 |
| + metreleptin | 49 | 0.404 |
| + volanesorsen | 11 | 0.107 |
| + both | 4 | 0.008 |
| + fenofibrate/omega-3 | 9 | 0.429 |

Volanesorsen **moves hepatic fat not at all**, yet beats metreleptin on this endpoint. In other words,
which endpoint you look at determines which drug looks better.

### Adding Inter-Individual Variability Turns Stratification into Response Rates

| Phenotype | n | Mean ΔHbA1c | Responders (ΔHbA1c ≤ −1.0%) |
|---|---|---|---|
| CGL | 200 | −2.01 | 81.2% |
| FPLD2, low leptin | 200 | −0.42 | 11.0% |
| FPLD2, normal leptin | 200 | −0.21 | 4.0% |

---

## A Structure That Was Broken and Discarded During Development

The first draft had leptin **directly suppress hepatic VLDL secretion** (since this is observed in the
literature). That structure predicted that **hepatic fat would worsen by 71% under pair-fed
conditions.** Which makes sense — with substrate delivery unchanged, suppressing the liver's largest
disposal route can only enlarge the pool. Since the observed direction is the opposite, **the drop in
VLDL secretion seen with metreleptin must be a *consequence* of the shrinking hepatic triglyceride
pool** (substrate-driven VLDL, Adiels). So the direct suppression term was removed (`KLVLDL = 0`, left
in the code so it remains testable), and leptin's direct hepatic action was reduced to **increased
fatty-acid oxidation · decreased SREBP-1c/lipogenesis** alone. After this fix, test 2 flipped from −71%
to **+57% survival**, matching the human data of Brown 2018 (PMID 29723161). The result was not
adjusted to fit the parameters — **the structure was fixed, and the result followed.**

For the same reason, `$ODE` carries the traces of two further fixes, left as comments.

- Unless the recycled NEFA flux is entered on **both sides** of the buffer balance, the futile cycle
  looks like "free storage capacity" (the reason the adipose compartment emptied out on the first run).
- Gating the storage flux by insulin **sensitivity** alone empties the residual compartment. In reality
  adipocytes see a 4–5-fold higher insulin concentration, so it has to be written as
  **action = concentration × sensitivity**, which then reproduces the *hypertrophic* adipocytes of
  partial lipodystrophy. → The claim that **capacity**, not signal, is the constraint comes from this
  one line.
- `$TABLE` was recomputing overflow with a different gate than `$ODE`, so the reported J_ov diverged
  from the value that actually drove the differential equations (the two storage fluxes now match
  exactly, verifiable in the validation log as `JST == JREL` — an internal consistency check on
  compartment balance).

---

## Honest Calibration Gaps

| Item | Model | Literature | Note |
|---|---|---|---|
| FPLD2 hepatic fat fraction | 5.9% | 10–18% (PMID 28199729, etc.) | **Under-predicted.** The only terms that raise hepatic fat in the partial form are visceral lipolysis's portal delivery (FHVIS) and the subcutaneous-fat competition term (KFLIVA); raising both further pushes CGL above 40%. A single Michaelis disposal curve cannot simultaneously fit normal 3%, partial 14%, and generalised 28% |
| FPLD2 HbA1c | 6.17% | 7.5–8.5% | Same cause as above (lower hepatic fat → lower IR). FPLD2's clinical spectrum runs from IGT to overt diabetes, and this phenotype corresponds to a moderate case |
| Healthy-subject hepatic fat fraction | 4.5% | 2–4% | Slightly over-predicted |
| CGL plasma TG | 1060 mg/dL | 500–3000+ | Within range but on the high side of the median |
| Metreleptin ΔTG (CGL) | −70.8% | −40 to −60% | Somewhat too strong |
| Metreleptin ΔHFF (CGL) | −86.0% | −30 to −50% (cohort) / −86% (PMID 12021250) | Matches Petersen's direct measurement, but stronger than the cohort average |
| Hepatic fat under ANGPTL3 blockade | +slight increase | Not observed clinically | `KEVVLDL` is the **least well-constrained** coefficient in this file and has been set low, at 0.05. The model does not yet fully agree with the absence of excess steatosis in humans with ANGPTL3 loss of function |

Parameters marked `(struct)` are values where **the literature supports the existence of the
relationship but the coefficient itself has not been fitted** (e.g., the pancreatic partition fraction
of overflow, the insulin slope of proteinuria).

---

## Usage

```r
# Model + scenarios
setwd("lipodystrophy-syndromes")
source("lipo_mrgsolve_model.R")

baseline_table()                    # 9 phenotypes, each to its own steady state
summarise_all(run_all())            # 28 scenarios
inference1_dose_vs_deficit()        # occupancy x deficit
inference2_pairfed()                # fixed intake
inference3_decoupling()             # TG / hepatic-fat dissociation
inference4_substrate()              # PPARgamma substrate dependence
inference5_diet_equivalence()       # J_in-matched diet (titrated)
inference6_window()                 # treatment window
inference7_ada()                    # neutralising antibodies
pancreatitis_burden()               # time above threshold
responder_rates(n = 200)            # inter-individual variability -> response rate

# Dashboard (12 tabs)
shiny::runApp("lipo_shiny_app.R")
```

```bash
# Re-render the map
dot -Tsvg lipo_qsp_model.dot -o lipo_qsp_model.svg
dot -Tpng -Gdpi=150 lipo_qsp_model.dot -o lipo_qsp_model.png
```

Full validation log (actually run under mrgsolve 2.0.1): [`lipo_validation_output.txt`](../../../lipodystrophy-syndromes/lipo_validation_output.txt)

---

## Model Composition

**53 ODE compartments** (time unit = days; a 4-hour half-life and a 20-year fibrosis trajectory coexist in the same system)

- Drug PK/PD (23): metreleptin, 3 compartments + neutralising-antibody titre · pioglitazone ·
  volanesorsen, 3 compartments + APOC3 indirect response · fenofibrate · evinacumab + ANGPTL3 ·
  omega-3 · metformin · exogenous insulin · GLP-1 RA · SGLT2i
- Leptin axis and intake (4): endogenous leptin · SOCS3 brake · delayed signal · intake state
- Adipose compartments and capacity (5): gluteofemoral · upper truncal/facial · visceral · capacity
  multiplier · adaptive oxidation
- Ectopic pools filled by overflow (5): liver · muscle · pancreatic islet · myocardium · kidney
- Plasma lipids (3): TG · NEFA · adiponectin
- Glycaemic axis (5): IR delayed state · blood glucose · insulin · beta cells · HbA1c
- Organ damage (8): ALT · hepatic inflammation · fibrosis · proteinuria · eGFR · androgen index ·
  cumulative time above the TG threshold · cumulative pancreatitis risk

**28 scenarios** — every claim has a matched control. 7/8/9 are the same mg/kg plus a dose increase,
13/14 are identical pioglitazone exposure with different substrate, 6/11 are the same dose with and
without neutralising antibodies, 12/12b are a pair-fed pair, and 23a/23b differ only in start time.

The map, with **192 nodes · 290 edges · 22 clusters**, spells out the seven structural inferences and
their falsification methods in cluster 22, and its node text has been updated to match the model's
actual run output (for example, the high-dose-insulin node reads not "raised hepatic fat" but the
actual result, "HbA1c −1.4%, hepatic fat and TG unchanged").

---

## References

[`lipo_references.md`](../../../lipodystrophy-syndromes/lipo_references.md) — **all 121 PMIDs fully
verified.** Every identifier was looked up via NCBI E-utilities, and the title, first author, journal,
and year **as returned by PubMed** were recorded. Of the 8 identifiers written from memory during the
first draft, 4 pointed to entirely unrelated papers (a herpesvirus dihydrofolate reductase,
linkage-analysis software, and so on); rather than guess a correction, they were deleted. Two
annotations record points where the model and the literature **disagree.**

---

## ⚠️ Disclaimer

This is a QSP model for education and research. It was built from published literature but has not
been independently verified or certified, and must not be used for clinical decision-making,
prescribing, or regulatory submission. The parameters are illustrative approximations; fitting and
validation against real patient data would be required separately.
