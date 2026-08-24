# Chronic Hepatitis D (HDV) — QSP Model

> **HDV is not a "virus to be killed" but a two-substrate assembly line.**
> The genome is made by the **host's** RNA polymerase II (entirely independent of HBV
> replication — which is why nucleos(t)ide analogues have no effect on HDV at all), while
> the envelope is made by **HBV cccDNA**. So every drug is classified not by
> "how potent" but by **which flux it interrupts**.

<p align="center">
  <a href="hdv_qsp_model_en.svg"><img src="hdv_qsp_model_en.png" width="880" alt="HDV QSP mechanistic map"></a>
</p>

---

## 1. The Model's One Idea

There are three fluxes that maintain the infected-hepatocyte pool `Id`, and two of them are independent of NTCP.

| Flux | Equation | NTCP Dependence | Drug That Interrupts It | Contribution to Baseline Influx |
|---|---|---|---|---|
| **(E) Entry / re-infection** | `βd · Vd · Ib · (1 − OCC)` | Fully dependent | **Bulevirtide** | **52.9%** |
| **(C) Cell-to-cell spread** | `k_cc · Id · Ib · (1 − φ·OCC)` | Partially dependent (φ = 0.5) | Only partially | 13.2% |
| **(D) Division-mediated spread** | When an infected cell divides, both daughter cells are infected | **Independent** | None | 33.8% |

Summing (D) with the NTCP-independent part of (C) gives **40.4% of the baseline influx** — this is the
**floor** that an entry inhibitor cannot cross. The evidence is Giersch 2019 (*Gut*): HDV persists
through liver regeneration and is amplified by cell division.

And all three are multiplied by **envelopment (A).** Because HDV cannot leave a cell without an HBsAg
envelope, lonafarnib (FTase → blocks prenylation) and HBsAg-targeting siRNA interrupt the same flux at
different points.

---

## 2. Only 4 Parameters Were Fitted — Everything Else Is a Prediction

| Fitted | Anchor | Result |
|---|---|---|
| `r_oatp`, `Kd_ntcp` | **Two** fold-increases in total bile acids: 2 mg ×3.2, 10 mg ×13.0 | `Kd_NTCP = 0.701 nM` |
| `dth_Id_immune` | MYR301, 2 mg, 48 weeks, **71% virologic response rate** | `0.02636 /day` |
| `ALT_base` | MYR301, 2 mg, 48 weeks, **51% ALT normalisation rate** | `42.0 U/L` |

Every other parameter that could be fixed by the steady-state condition (`d/dt = 0`) was back-solved,
so **the untreated control is a steady state by construction** (drift ≤ 1.8% after a 2-year simulation).

### 2.1 Bile Acids Are Back-Solved into Target Occupancy — Occupancy Is Derived, Not Assumed

NTCP is both the HDV receptor and a bile-acid transporter. So at steady state

```
TBA_fold = (1 + r) / ((1 − occupancy) + r),   r = CL_OATP / CL_NTCP
```

holds, and the two published bile-acid fold-increases alone determine `Kd` and `r`.

| Dose | Css | **Derived occupancy** | Free NTCP | Bile acids |
|---|---|---|---|---|
| 2 mg | 1.72 nM | **0.710** | 29.0% | ×3.2 |
| 10 mg | 14.47 nM | **0.954** | 4.6% | ×13.0 |

**The approved 2 mg dose does not saturate NTCP.** The residual entry flux differs **6.3-fold** between
the two doses.

Independent check (not used in the fit): the back-solved `Kd = 0.70 nM` falls within the **reported in
vitro NTCP affinity (80 pM to a few nM)** of myrcludex B / bulevirtide. This value was obtained from
bile-acid pharmacology alone, without ever looking at binding-assay data.

---

## 3. Predicted Results (Not Used in the Fit)

### 3.1 Bulevirtide Has No First Phase — A Structural Consequence

An entry inhibitor cannot touch **already-infected cells**, so it cannot reduce secretion on day 1.
Interferon (blocks replication/production) and lonafarnib (blocks assembly) can. The model was never
taught this — it simply follows from which flux each drug interrupts.

| Regimen | d1 | d7 | d28 | d84 | d336 | Week-1 slope (log10/day) |
|---|---|---|---|---|---|---|
| BLV 2 mg | −0.00 | −0.03 | −0.18 | −0.60 | −2.70 | **−0.0050** |
| BLV 10 mg | −0.00 | −0.05 | −0.25 | −0.84 | −3.77 | −0.0071 |
| PegIFN alfa 180 µg | −0.12 | −0.46 | −0.68 | −1.25 | −4.07 | **−0.0661** |
| LNF 50 + RTV 100 BID | −0.19 | −0.62 | −0.78 | −1.18 | −3.17 | **−0.0889** |

### 3.2 Reproducing MYR301 (48 Weeks, Virtual Population n = 300)

| Regimen | Mean Δlog10 | RNA response | ALT normalisation | **Combined response** |
|---|---|---|---|---|
| Untreated (NUC only) | −0.00 | 0.0% | 0.7% | 0.0% |
| BLV 2 mg | −2.77 | 69.7% ◄anchor | 51.3% ◄anchor | **36.3%** |
| BLV 10 mg | −3.87 | 90.3% | 53.7% | **48.3%** |
| **Hypothetical complete entry blockade** | −4.07 | 92.7% | 54.0% | **50.3%** |
| *Observed (MYR301)* | | *4 / 71 / 76%* | *12 / 51 / 56%* | *2 / 45 / **48**%* |

**The central result.** Even though the residual entry flux falls 6.3-fold, the combined response moves
only from 36% to 48% (observed: 45% to 48%). This is because entry is only 52.9% of the influx, and the
rest lies where an entry inhibitor cannot reach.

**The ceiling for the entire drug class.** Even a hypothetical, non-toxic entry inhibitor at 100%
occupancy achieves only a **50.3%** combined response — 14.0 points above the approved dose. Relative
to no treatment, the entry-inhibition axis is already **72% exhausted** at 2 mg/day. **Dose escalation
is not the lever. Adding a mechanism is the lever.**

### 3.3 Only Interferon Produces a Sustained Response After Treatment Ends (a Prediction, Not a Fit)

Interferon does not reduce influx — it is the only current mechanism that **raises the death rate** of
infected cells and reverses exhaustion.

| Regimen | Δlog10 at end of treatment | RNA < LOD at +24 weeks | *Observed reference* |
|---|---|---|---|
| PegIFN alfa, 48 weeks | −4.10 | **26.7%** | *HIDIT-1 ~26–28%* |
| PegIFN alfa, 96 weeks | −8.55 | 93.3% | *HIDIT-2 ~31%* |
| PegIFN lambda, 48 weeks | −4.10 | 26.7% | *LIMT-1 36% (5/14)* |
| BLV 2 mg + PegIFN, 48 weeks | −5.45 | 55.0% | *MYR204 ~45% (10 mg combination)* |
| BLV 10 mg + PegIFN, 48 weeks | −5.91 | 63.3% | |
| BLV 2 mg, continuous 96 weeks | −6.30 | 61.3% (on treatment) | Relapses on discontinuation |

The synergy is not empirical but **structural**. Bulevirtide reduces the **influx** to the same pool,
and interferon raises the **outflux** — additive in the net growth rate, but because the endpoint is a
threshold on a decaying exponential, it appears super-additive.

### 3.4 The Lonafarnib Paradox — Blood Levels Fall While Intracellular Levels Rise

Blocking FTase stops assembly without stopping replication. The genome becomes trapped inside the cell.

| Week | Blood HDV RNA (log10) | Intracellular Rg (copies/cell) |
|---|---|---|
| 0 | 5.50 | 3000 |
| 1 | 4.88 | **4583 (+53%)** |
| 48 | 2.33 | 4581 |

**Blood HDV RNA measures the assembly flux, not the size of the reservoir.** Stopping at week 48
releases the trapped pool, producing a **+0.59 log rebound within 7 days** that then plateaus
(2.33 → 2.92 → 2.90 at day 84). In other words the rebound is **fast but partial** — release is
immediate once the assembly block is lifted, but the infected-cell pool itself has already shrunk.
Lonafarnib is the only current drug that drives these two quantities (blood vs. intracellular) apart.

### 3.5 ALT and HDV RNA Read Different Things

In the model, the ALT driving equation is `kill = (innate + CTL + IFN)·Id + κ·(entry flux)`.
Blocking entry removes a term **within days**, while `Id` (and hence blood RNA) declines **over
months**.

| Week | ALT (% of baseline) | Δlog10 HDV RNA |
|---|---|---|
| 4 | 68.9% | −0.18 |
| 12 | 50.5% | −0.60 |
| 48 | 38.3% | −2.70 |

By week 4, ALT has already fallen 31%, while HDV RNA has dropped only 0.18 log. This is the model's
explanation for why even virologic non-responders in MYR301 had ALT normalise, and it is **a
hypothesis** — the falsification condition is in section 5 below.

### 3.6 Bile Acids as a Therapeutic-Index Gauge (Efficacy Saturates, Cost Does Not)

| Dose | Occupancy | Free NTCP | Bile acids | Combined response |
|---|---|---|---|---|
| 0.5 mg | 0.305 | 69.5% | ×1.4 | 2.7% |
| 1 mg | 0.498 | 50.2% | ×1.9 | 14.0% |
| **2 mg** | **0.710** | **29.0%** | **×3.2** | **33.3%** |
| 5 mg | 0.897 | 10.3% | ×7.6 | 45.3% |
| 10 mg | 0.954 | 4.6% | ×13.0 | 46.0% |
| 20 mg | 0.978 | 2.2% | ×18.8 | 46.7% |

The approved 2 mg dose sits at the knee of the curve — **without saturating the target.** This means
the usual reason a dose-response curve flattens (receptor saturation) is not what's happening here.

### 3.7 What Determines the Floor (Mean Δlog10 at 48 Weeks, ±30% Perturbation)

| Parameter | −30% | Baseline | +30% | \|Impact\| | Role |
|---|---|---|---|---|---|
| **dth_Id_immune** | −1.94 | −2.69 | −3.50 | **0.783** | Killing of infected cells (differential death rate) |
| Kd_ntcp (target affinity) | −3.01 | −2.69 | −2.46 | 0.276 | Occupancy |
| frac_cc | −2.78 | −2.69 | −2.60 | 0.094 | Share of cell-to-cell spread |
| phi_ntcp | −2.60 | −2.69 | −2.78 | 0.091 | NTCP-dependent fraction of cell-to-cell spread |
| k_cure | −2.65 | −2.69 | −2.73 | 0.044 | Non-cytolytic clearance |
| **d_hep** | −2.69 | −2.69 | −2.69 | **0.000** | Hepatocyte turnover rate |

Sum of host parameters, 1.012, versus target affinity, 0.276 — host cell biology wins by **3.7-fold.**

**And one unexpected result.** `d_hep` (hepatocyte turnover rate) — the very parameter the
division-mediated floor is built on — has an impact of **exactly zero.** This is not numerical error;
it is algebra I did not do on purpose. Raising the background death rate raises the infected cells'
death rate and the division rate that refills their place **by exactly the same amount**, so they
cancel out precisely in the net dilution rate. What sets the floor is not how fast the liver
regenerates, but **how much faster** an infected hepatocyte dies than a neighbouring HBsAg-positive
cell (the differential, `dth_Id_immune`). This does not weaken the therapeutic argument — it sharpens
it. A drug that regenerates the whole liver faster does nothing; only something that makes infected
cells die **preferentially** (interferon, checkpoint release, therapeutic vaccination) can move the
floor.

---

## 4. Files

| File | Contents |
|---|---|
| [`hdv_qsp_model_en.dot`](hdv_qsp_model_en.dot) · [`.svg`](hdv_qsp_model_en.svg) · [`.png`](hdv_qsp_model_en.png) | Mechanistic map — **18 clusters, 157 nodes** |
| [`hdv_mrgsolve_model.R`](hdv_mrgsolve_model.R) | mrgsolve ODE model — **33 compartments, 14 scenarios** + bile-acid dose-response and floor-sensitivity helpers |
| [`hdv_reference_model.py`](hdv_reference_model.py) | A **dependency-free** Python reference implementation + virtual population (n = 300) + fitting routine |
| [`hdv_model_report.txt`](hdv_model_report.txt) | The computed results behind every figure above (A0–A14). Regenerate with `python3 hdv_reference_model.py` |
| [`hdv_shiny_app_en.R`](hdv_shiny_app_en.R) | Shiny dashboard — **9 tabs** (patient/floor · PK and occupancy · virology · endpoints · ALT-RNA dissociation · floor and ceiling · bile acids · scenario comparison · long-term outcomes) |
| [`hdv_references_en.md`](hdv_references_en.md) | **62** references, with PMIDs directly confirmed on PubMed marked separately from search links for unconfirmed entries |

### Reproduction

```bash
python3 hdv_reference_model.py                          # -> hdv_model_report.txt (about 25 minutes)
dot -Tsvg hdv_qsp_model_en.dot -o hdv_qsp_model_en.svg
dot -Tpng -Gdpi=150 hdv_qsp_model_en.dot -o hdv_qsp_model_en.png
```
```r
library(mrgsolve); mod <- mread("hdv_mrgsolve_model.R")
print(scenario_summary(mod)); print(blv_dose_response(mod)); print(floor_sensitivity(mod))
shiny::runApp("hdv_shiny_app_en.R")
```

---

## 5. Where This Model Is Wrong (Full List in Report Section A13)

Reported plainly, without smoothing it over.

1. **Dose separation — the clearest failure, and what was learned from it.** Tied to the bile-acid
   anchor, the model spreads 2 mg and 10 mg apart more than observed on the **virologic** endpoint
   (69.7% vs. 90.3%, observed 71% vs. 76%). This washes out in the combined response (36% vs. 48%,
   observed 45% vs. 48%), but much of that flattening comes from the **ALT ceiling**, not the floor.
   In other words, trying to explain both observations with one mechanism revealed that two mechanisms
   are needed. One of three possible interpretations is testable — in the same cells and the same
   drug, **HDV entry inhibition must saturate at a lower occupancy than bile-acid transport
   inhibition.**
2. **The ALT-RNA dissociation term (`κ·entry`) is a hypothesis.** Falsification condition: even in
   virologic non-responders, ALT should fall in proportion to (1 − occupancy). If, patient by patient,
   the ALT decrease tracks the HDV RNA decrease, this term is wrong.
3. **The floor parameters are not mutually identifiable.** Falsification/discrimination condition:
   serial liver-tissue HDAg staining during treatment. If division-mediated persistence dominates, the
   HDAg-positive **fraction** should stay nearly unchanged while only blood RNA falls; if cell-to-cell
   spread dominates, the fraction should fall too.
4. **The combined response for lonafarnib/siRNA is over-predicted.** D-LIVR reported a ~10% combined
   response for lonafarnib/ritonavir, while the model runs several-fold higher. The missing piece is
   probably the **cytotoxicity/immunogenicity of the trapped intracellular genome and unprenylated
   L-HDAg itself**, in which case blocking assembly should **add** a damage term rather than remove
   one.
5. **There is no adherence, dose reduction, or toxicity-driven discontinuation at all.** Lonafarnib's
   gastrointestinal toxicity and interferon's haematologic toxicity are the main reasons real-world
   regimens underperform their pharmacology, but every number in this file is a **per-protocol**
   number.
6. **Confirmed that the variance of the floor (`cv_floor`) across the population is not identifiable**
   and excluded it from the fit (the combined response is flat at 34–37% across the 0.10–1.60 range).
   What moves the response rate is the floor's **mean, not its variance.**
7. No spatial structure (cell-to-cell spread is treated as well-mixed mass action), a single HBsAg pool
   (HBsAg from integrated HBV DNA is not separated out), and no genotypes (everything is HDV-1).

---

## 6. Disease Background (One Paragraph)

HDV is a 1.7 kb circular negative-strand RNA virusoid that has no polymerase of its own, so it uses the
host's RNA Pol II for double rolling-circle replication, and it **must borrow HBV's HBsAg for its
envelope.** When ADAR1 edits the amber/W site (UAG→UGG), S-HDAg (p24) becomes L-HDAg (p27), 19 amino
acids longer, which gains a CXXQ prenylation motif (Cys²¹¹) that enables assembly. About 5% of
HBsAg-positive people (12 to 72 million) are infected, and among the chronic viral hepatitides, HDV has
the **fastest fibrosis progression** (the model's untreated progression rate is 0.149 Ishak units/year;
the literature reports 0.15–0.25, versus ~0.10 for HBV alone). The only approved treatment is
**bulevirtide** (an NTCP entry inhibitor, 2 mg SC once daily); peg-IFN alfa is used off-label, and
lonafarnib/ritonavir, peg-IFN lambda, and HBsAg-targeting siRNAs/NAPs are in development. **Only HBsAg
loss (a functional cure) is a true cure for HDV** — without an envelope, HDV cannot leave a cell.

---

> ⚠️ **Disclaimer.** This is a QSP model for education and research. It was built from published
> literature and clinical trial data but has not been independently verified, and must not be used for
> actual clinical decision-making, prescribing, or regulatory submission. `hdv_references_en.md` states
> the source of each parameter and flags the items marked "model hypothesis"; `hdv_model_report.txt`
> section A13 states the falsification condition for each hypothesis.
