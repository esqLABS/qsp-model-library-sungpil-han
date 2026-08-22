# Invasive Pulmonary Aspergillosis QSP Model
### Invasive Pulmonary Aspergillosis · IPA

<p align="center">
  <a href="../../../invasive-aspergillosis/ipa_qsp_model.svg">
    <img src="../../../invasive-aspergillosis/ipa_qsp_model.png" width="880" alt="Invasive pulmonary aspergillosis QSP mechanistic map">
  </a>
</p>

| File | Contents |
|---|---|
| [`ipa_qsp_model.dot`](../../../invasive-aspergillosis/ipa_qsp_model.dot) · [SVG](../../../invasive-aspergillosis/ipa_qsp_model.svg) · [PNG](../../../invasive-aspergillosis/ipa_qsp_model.png) | Mechanistic map — 171 nodes, 17 clusters, 259 edges |
| [`ipa_mrgsolve_model.R`](../../../invasive-aspergillosis/ipa_mrgsolve_model.R) | mrgsolve ODE model (53 compartments, 5 antifungals, 2 concomitant drugs) + 23 scenarios |
| [`ipa_shiny_app.R`](../../../invasive-aspergillosis/ipa_shiny_app.R) | Shiny dashboard (11 tabs) |
| [`ipa_references.md`](../../../invasive-aspergillosis/ipa_references.md) | 99 references — every PMID looked up and verified through NCBI E-utilities |
| [`ipa_reference_model.py`](../../../invasive-aspergillosis/ipa_reference_model.py) | Independent Python/scipy re-implementation — for verification |
| [`ipa_reference_output.txt`](../../../invasive-aspergillosis/ipa_reference_output.txt) · [`ipa_scenario_results.json`](../../../invasive-aspergillosis/ipa_scenario_results.json) | The computed output behind every number below |

> **Note**: this repository already contains an [allergic bronchopulmonary aspergillosis (ABPA)](../../../allergic-bronchopulmonary-aspergillosis/)
> model. Same fungus, but **the opposite disease**. ABPA is a Th2 hypersensitivity
> that arises because the immune system over-reacts; IPA is tissue invasion that
> arises because there is no immune system. In ABPA, steroids are the treatment;
> in IPA, steroids are the cause.

---

## What this model says, in one sentence

**The azole multiplies the growth term, and only neutrophils and the polyenes add
to the elimination term. This single arithmetic asymmetry accounts for very nearly
the whole of the pharmacotherapeutics of invasive aspergillosis.**

```
   dB/dt  =  k_grow · (1 − I_azole) · B   −   (k_host · N_eff + k_cidal) · B
             └─ where triazoles act ──┘       └─ where neutrophils and ────┘
                                                 polyenes act
```

The triazoles deplete ergosterol and so slow hyphal elongation. They do not lyse
hyphae. That is why a triazole **never appears in the elimination term**. The
consequences:

- The maximum value `Imax` of `I_azole` is structurally less than 1 (0.955 in this
  model).
- Therefore, with no neutrophils (`N_eff = 0`) and no polyene,
  `dB/dt ≥ k_grow·(1 − Imax)·B > 0` holds **always**.
- That is, **an azole alone cannot reduce the fungal burden of a host without
  neutrophils.** It can only make it slower.

Computing that lower bound as a number (model parameters; the calculation is in
`ipa_reference_model.py` §C):

| Quantity | Value |
|---|---|
| Growth lower bound at maximum azole effect | **0.00428 /h = 0.045 log10/day** |
| Doubling time in that state | **162 h (6.8 days)** |
| Minimum neutrophil function needed to make `dB/dt < 0` | **1.3% of maximal function** |

The last line governs the clinical behaviour of this disease. The azole is **a
bridge until the marrow comes back**, and what is needed is not normal
neutrophils but 1.3% of maximal function. If that 1.3% never arrives, the azole
fails at any dose.

---

## 1. Pharmacokinetic verification (against the label)

Drug was given without infection and steady-state exposure (day 12–13) computed.

| Regimen | Cmin (mg/L) | Cmax (mg/L) | AUC24 (mg·h/L) |
|---|---|---|---|
| Voriconazole 4 mg/kg q12h IV, **UM (\*17/\*17)** | 0.37 | 3.17 | 25.0 |
| Voriconazole 4 mg/kg q12h IV, RM (\*1/\*17) | 0.64 | 3.63 | 34.0 |
| Voriconazole 4 mg/kg q12h IV, **NM (\*1/\*1)** | 1.13 | 4.27 | 47.9 |
| Voriconazole 4 mg/kg q12h IV, IM (\*1/\*2) | 2.30 | 5.58 | 78.3 |
| Voriconazole 4 mg/kg q12h IV, **PM (\*2/\*2)** | 4.56 | 7.92 | 133.8 |
| Isavuconazole 200 mg q24h IV | 3.14 | 5.22 | 92.5 |
| Posaconazole DR tablet 300 mg q24h PO | 1.13 | 1.60 | 33.7 |
| Liposomal amphotericin B 3 mg/kg q24h | 13.15 | 73.25 | 552.4 |
| Liposomal amphotericin B 10 mg/kg q24h | 43.84 | 244.15 | 1841.4 |
| Anidulafungin 100 mg q24h | 3.38 | 6.46 | 110.7 |

**At an identical mg/kg dose, genotype alone spreads AUC24 from 25.0 to 133.8 —
a 5.4-fold range.** Weight-based dose calculation does not touch that variance at
all. This is why voriconazole needs TDM, and why TDM is *a device for reducing
variance rather than a device for raising mean exposure*.

### 1-1. And dose-exposure is not linear

The CYP2C19 route saturates (Michaelis–Menten) while only the CYP3A4/FMO3 route
is linear.

| Dose (mg/kg q12h) | AUC24 | Fold vs 3 mg/kg | Dose fold |
|---|---|---|---|
| 2.0 | 19.2 | 0.60 | 0.67 |
| 3.0 | 32.2 | 1.00 | 1.00 |
| 4.0 | 47.9 | 1.49 | 1.33 |
| 5.0 | 66.8 | 2.07 | 1.67 |
| 6.0 | 89.2 | 2.77 | 2.00 |
| 8.0 | 145.8 | 4.53 | 2.67 |

**Doubling the dose from 4 to 8 mg/kg makes exposure 3.04-fold.** Voriconazole's
therapeutic range is narrow not because the PD curve is steep but **because the
PK is non-linear**. This is not a conclusion of the model but a direct arithmetic
consequence of the structure put into it (saturable metabolism), and it runs in
the same direction as the supraproportionality reported by Purkins 2002.

---

## 2. Twenty-one treatment scenarios

Fungal inoculation on day 0, treatment started on day 3 (unless stated
otherwise), 12 weeks of observation. `logB` = log10 hyphal burden (CFUe),
`PERF` = minimum lesional perfusion index.

| # | Scenario | logB d14 | logB d42 | GM peak | PERF min | 6-week mortality % | 12-week mortality % |
|---|---|---|---|---|---|---|---|
| S01 | No treatment, persistent neutropenia | 10.48 | 10.48 | 2.21 | 0.14 | 96.7 | **100.0** |
| S02 | No treatment, neutrophil recovery d10 | 7.75 | 0.00 | 4.40 | 0.24 | 78.6 | 91.7 |
| S03 | **Voriconazole, persistent neutropenia** | 5.71 | **7.54** | 2.15 | 0.14 | 52.7 | **98.7** |
| S04 | Voriconazole NM, recovery d10 | 0.00 | 0.00 | 0.00 | 1.00 | 13.4 | **15.5** |
| S05 | Voriconazole UM (\*17/\*17), recovery d10 | 0.00 | 0.00 | 0.00 | 1.00 | 14.5 | 16.5 |
| S06 | Voriconazole PM (\*2/\*2), recovery d10 | 0.00 | 0.00 | 0.00 | 1.00 | 23.9 | **27.7** |
| S07 | UM + TDM dose increase (7.6 mg/kg) | 0.00 | 0.00 | 0.00 | 1.00 | 13.5 | 15.6 |
| S08 | Isavuconazole, recovery d10 | 0.00 | 0.00 | 0.00 | 1.00 | 13.5 | 15.6 |
| S09 | L-AmB 3 mg/kg, recovery d10 | 0.00 | 0.00 | 0.00 | 1.00 | 9.6 | **11.7** |
| S10 | **L-AmB 10 mg/kg, recovery d10** | 0.00 | 0.00 | 0.00 | 1.00 | 14.0 | **16.5** |
| S11 | Anidulafungin monotherapy, recovery d10 | 0.86 | 0.00 | 0.18 | 0.94 | 48.5 | 77.7 |
| S12 | Voriconazole + anidulafungin | 0.00 | 0.00 | 0.00 | 1.00 | 12.9 | 15.0 |
| S13 | TR34/L98H (VRC MIC 8) + voriconazole | 7.30 | 0.00 | 5.00 | 0.33 | 71.9 | 88.8 |
| S14 | TR34/L98H, switched to L-AmB on d10 | 3.41 | 0.00 | 4.81 | 0.44 | 70.9 | 88.2 |
| S15 | Steroid host (normal neutrophils) + voriconazole | 0.00 | 0.00 | 0.00 | 1.00 | 2.4 | 4.7 |
| S16 | Steroid host, taper from d14 | 0.00 | 0.00 | 0.00 | 1.00 | 2.4 | 4.7 |
| S17 | Voriconazole + G-CSF (from d3) | 0.00 | 0.00 | 0.00 | 1.00 | 12.9 | 15.0 |
| S18 | **Late start on d10, voriconazole** | 5.27 | 0.00 | 2.58 | 0.29 | 76.6 | **90.8** |
| S19 | Central nervous system involvement, voriconazole | 0.00 | 0.00 | 0.00 | 1.00 | 13.4 | 15.5 |
| S20 | Central nervous system involvement, L-AmB 5 mg/kg | 0.00 | 0.00 | 0.00 | 1.00 | 9.1 | 11.3 |
| S21 | Transplant, tacrolimus + voriconazole | 0.00 | 0.00 | 0.00 | 1.00 | 2.4 | 4.7 |

Three lines to read:

- **S03 (azole + persistent neutropenia)**: the burden falls to 10^5.71 by day 14
  and then **climbs back to 10^7.54 by day 42.** The drug is still being given and
  the fungus grows. This is the growth lower bound above, drawn as a graph. The
  azole has not failed; it has done exactly what an azole can do.
- **S09 vs S10 (L-AmB 3 vs 10 mg/kg)**: the high dose is **worse** (11.7% vs
  16.5%). The burden is completely eliminated in both arms, so there is nothing to
  gain on efficacy and only nephrotoxicity is added. This runs in the same
  direction as what the AmBiLoad trial actually observed.
- **S06 (PM genotype)**: fixed-dose PM is 27.7% against 15.5% for NM. The burden
  goes to 0 in both — **the whole difference is toxicity** (see §6).

---

## 3. Treatment delay: not a slope but a bifurcation

Voriconazole NM, neutrophil recovery on day 14. Only the day treatment starts is
varied.

| Start day | PERF min | 12-week mortality |
|---|---|---|
| 1 | 1.000 | 12.7 % |
| 2 | 1.000 | 15.2 % |
| 3 | 1.000 | 20.3 % |
| 4 | 0.998 | 25.0 % |
| 5 | 0.984 | **29.7 %** |
| 6 | 0.879 | **96.0 %** |
| 7 | 0.487 | 97.9 % |
| 8 | 0.212 | 98.9 % |
| 10 | 0.168 | 99.4 % |
| 14 | 0.168 | 99.5 % |

**Days 1–5 are a gentle slope (12.7 → 29.7%), and the single day from 5 to 6
moves 66.3 percentage points.** Compared with 17.0 points for the whole of the
preceding four days, this is not a dose-response but **a threshold**.

The position of that threshold is parameter-dependent, so do not carry the number
"6 days" into the clinic. What can be carried is the **shape**: the value of early
diagnosis is not linear, and past some day the same drug abruptly stops working.
This is the model's explanation of why a strategy of chasing the chest CT halo
sign to give pre-emptive therapy halves mortality, and at the same time a
falsifiable prediction.

---

## 4. A hypothesis of mine that my own model refuted

The map (§cluster 7) draws an orange positive feedback loop. **The hyphae destroy
the vessels, and those vessels deliver the drug that would kill the hyphae** — so
the reason late treatment fails would be that the drug cannot get there. That was
my hypothesis.

The way to check it is to delete it. Setting `DELIV_FLOOR = 1` means the drug
reaches the hyphae at the full ELF concentration no matter how thrombosed the
lesion is.

| Start day | Peak logB (feedback present) | Peak logB (ablated) | 12-week mortality (feedback) | 12-week mortality (ablated) |
|---|---|---|---|---|
| 3 | 5.72 | 5.72 | 20.3 | 20.3 |
| 5 | 7.56 | 7.56 | 29.7 | 29.7 |
| 6 | 8.49 | 8.48 | 96.0 | 96.0 |
| 8 | 10.12 | 10.07 | 98.9 | 98.9 |
| 14 | 10.48 | 10.48 | 99.5 | 99.5 |

**First failing start day: 6 with the feedback, 6 with it ablated. The cliff does
not move at all.**

The hypothesis was wrong. The threshold of treatment delay is **not a delivery
problem** but a question of whether the fungal burden at the moment a
growth-rate inhibitor takes over exceeds the size that drug can hold until the
marrow returns. The angioinvasion loop is real inside the model and it moves
perfusion, oxygenation and mortality risk, but **it is not the mechanism of the
cliff.** I have deleted that claim from the headers of the map and the R model
(the title of map cluster 7 and the `DELIV → BIFURC` edge are now labelled
"ablated: no effect").

---

## 5. Putting the drug lever and the marrow lever in the same units

12-week mortality %, same virtual patient.

| Neutrophil recovery day | No treatment | Voriconazole | Isavuconazole | L-AmB 3 | VRC + echinocandin |
|---|---|---|---|---|---|
| 5 | 13.1 | 9.7 | 9.7 | 9.0 | 9.6 |
| 10 | 91.7 | 15.5 | 15.6 | 11.7 | 15.0 |
| 14 | 94.9 | 20.3 | 20.5 | 13.8 | 19.3 |
| 21 | 97.4 | 29.2 | 29.2 | 17.4 | 26.8 |
| 28 | 98.6 | 38.2 | 37.9 | 20.8 | 34.2 |
| **No recovery** | 100.0 | **98.7** | **97.1** | **43.2** | **97.9** |

| Comparison | Spread (percentage points) |
|---|---|
| Difference among the four active regimens, recovery fixed at d10 | **3.9** |
| Difference among the four active regimens, recovery fixed at d14 | 6.7 |
| Difference among the four active regimens, no recovery | 55.5 |
| Difference among hosts, voriconazole fixed | **89.0** |
| Difference among hosts, L-AmB fixed | 34.2 |

**In the patients a clinical trial enrols — those whose marrow returns within two
weeks — the entire spread available to drug choice is 3.9–6.7 points, while the
spread the marrow produces at a fixed drug is 89.0 points.** This is the model's
explanation of why non-inferiority trials only ever show non-inferiority.

And the last row is this model's most testable prediction: **the polyene separates
from the azoles only in the patient whose marrow does not come back** (43.2% vs
97–99%). The reason is exactly the arithmetic of §0 — the polyene enters the
elimination term and the azole does not. If the trials that compared azoles
against polyenes did not enrol enough of this subgroup, the model claims those
trials missed a difference that really exists.

---

## 6. Galactomannan looks at flux, not at stock

The model makes GM release **proportional to flux** rather than **proportional to
burden**:

```
GM production = k_growth · (growth flux) + k_lysis · (death flux)
```

Growing hyphae shed GM, and dying hyphae spill what is left of theirs. A standing
mass of hyphae does not itself make GM. The result in the persistently
neutropenic host (the population in which GM is actually measured):

| Arm | GM d7 | GM d10 | GM d14 | GM d21 | logB d10 | logB d14 | logB d21 |
|---|---|---|---|---|---|---|---|
| No treatment | 0.06 | **2.14** | **0.22** | 0.00 | 10.46 | 10.48 | 10.48 |
| Voriconazole | 0.00 | 0.00 | 0.00 | 0.00 | 5.45 | 5.71 | 6.17 |
| L-AmB 3 mg/kg | 0.00 | 0.00 | 0.00 | 0.00 | 0.80 | 0.00 | 0.00 |
| Anidulafungin | 0.00 | 0.03 | **1.09** | 0.09 | 8.52 | 10.20 | 10.48 |

Three results come out of this.

1. **GM falls in the untreated patient.** 2.14 on day 10 → 0.22 on day 14. The
   burden has not come down; it has **reached carrying capacity** (10.46 → 10.48).
   The growth flux has gone to zero, so GM production has stopped. **A falling GM
   is not a treated infection.**
2. **The day-14 ranking inverts.** Anidulafungin GM 1.09 (burden 10.20) is
   **higher** than no-treatment GM 0.22 (burden 10.48). The higher GM belongs to
   the lower burden.
3. **On voriconazole GM never once exceeds 0.5 while the burden is 10^5.7** — a
   GM-negative living infection. This is the model's version of the "reduced GM
   sensitivity under mould-active therapy" reported by Marr 2005, and here it is
   not an assumption but is **derived** as a consequence of writing GM against
   flux.

---

## 7. What TDM buys is not efficacy but avoided toxicity

Voriconazole, neutrophil recovery d14, MIC 0.5.

| Strategy | AUC24 | Peak ALT (U/L) | Neurological/visual effect AUC | 12-week mortality |
|---|---|---|---|---|
| UM, fixed 4 mg/kg | 25.1 | 26 | 0.2 | 22.7 % |
| UM, TDM increase to 7.6 mg/kg | 62.6 | 67 | 180.1 | **20.5 %** |
| NM, fixed 4 mg/kg | 47.9 | 27 | 4.2 | 20.3 % |
| **PM, fixed 4 mg/kg** | 134.6 | **255** | **929.1** | **31.6 %** |
| **PM, TDM reduction to 1.8 mg/kg** | 47.0 | 31 | 2.4 | **20.1 %** |

- **What a TDM dose increase buys in a UM: 2.2 points.** Exposure was raised
  2.5-fold and this is all that comes back on efficacy. It is because the patient
  is already on the flat part of the efficacy curve.
- **What a TDM dose reduction buys in a PM: 11.5 points.** And all 11.5 come from
  toxicity — ALT 255 → 31, neurological effect AUC 929 → 2.4. The burden is 0 in
  both.

**So TDM is asymmetrically useful.** A little on the under-exposed side, five
times as much on the over-exposed side. CYP2C19 genotyping done in advance tells
you which way you are going to be wrong, and that is why CPIC recommends it.

---

## 8. The MIC ladder and the breakpoint

Voriconazole, neutrophil recovery d14. `fAUC24` = 20.1 mg·h/L (free fraction 0.42
applied).

| MIC (mg/L) | fAUC24/MIC | 12-week mortality |
|---|---|---|
| 0.125 | 161.0 | 19.6 % |
| 0.25 | 80.5 | 19.8 % |
| 0.5 | 40.2 | 20.3 % |
| 1.0 | 20.1 | 22.3 % |
| **2.0** | **10.1** | **29.0 %** |
| **4.0** | **5.0** | **91.7 %** |
| 8.0 | 2.5 | 94.3 % |

From MIC 0.125 to 1 the curve is essentially flat (19.6 → 22.3%). The cliff lies
**between MIC 2 and 4**, between `fAUC/MIC` 10 and 5. That is one to two doubling
dilutions above the EUCAST voriconazole clinical breakpoint (susceptible
≤1 mg/L) — **meaning the model predicts the breakpoint slightly too permissively,
and I record this as an approximation failure rather than a success.**

## 9. Resistance comes from the environment — the model says so itself

In-host resistance selection was left switched on at `k_mut = 2×10⁻⁸`/replication.
The resistant subpopulation fraction after six weeks:

| Arm | Resistant subpopulation fraction (day 42) | Total burden log10 |
|---|---|---|
| No treatment | 1.8×10⁻⁴⁶ | 0.00 |
| Voriconazole monotherapy | 1.8×10⁻⁶¹ | 0.00 |
| **Voriconazole, persistent neutropenia** | **3.4×10⁻⁷** | **7.54** |
| L-AmB 3 mg/kg | 5.7×10⁻¹²⁰ | 0.00 |
| Voriconazole + anidulafungin | 2.9×10⁻⁶¹ | 0.00 |

**In no arm does resistance reach a clinically meaningful fraction within six
weeks.** There is only one condition under which it could: the fungus has to keep
replicating under azole pressure, and that happens only in the patient whose
neutrophils never return. Even then it stops at 3.4×10⁻⁷.

The model's conclusion: **in-host selection is far too slow to be the main route
to the azole resistance observed clinically.** This agrees with the epidemiological
observation (Verweij, Snelders, Chowdhary, Meis) that the dominant route is
inhaling environmental TR34/L98H · TR46/Y121F strains, derived from agricultural
azole fungicides, from the outset — and the model did not assume that but showed
it by leaving the in-host route switched on and finding it insufficient.

And S13–S14 are the clinical price of that: against an MIC 8 strain voriconazole
gives 88.8%, and switching to L-AmB on day 10 still gives 88.2% — **switch late
and the switch is worth nothing.**

---

## 10. The azole-tacrolimus interaction

Before an azole is an antifungal it is a CYP3A4 inhibitor.

| Arm | Tacrolimus trough d14 (ng/mL) | Fold vs alone |
|---|---|---|
| Tacrolimus alone | 6.92 | 1.00 |
| + Voriconazole | 18.22 | **2.63** |
| + Isavuconazole | 10.00 | **1.45** |
| Tacrolimus 1/3 dose + voriconazole | 6.07 | 0.88 |

This is the same magnitude as the reported values (voriconazole about 3-fold,
isavuconazole about 1.4-fold), and **a reduction to one third restoring the target
concentration** is the actual practice of the transplant ward.

---

## 11. Central nervous system involvement — where the model contradicts clinical received wisdom

| Arm | Mean brain concentration (mg/L) | Kp | Brain burden day 42 |
|---|---|---|---|
| Voriconazole 4 mg/kg q12h | 0.414 | 0.50 | eliminated (<1 CFUe) |
| L-AmB 5 mg/kg | **0.751** | 0.045 | eliminated (<1 CFUe) |
| Isavuconazole 200 mg qd | 0.281 | 0.11 | eliminated (<1 CFUe) |
| Anidulafungin 100 mg qd | 0.039 | 0.02 | **10^8.69** |
| No treatment | 0.000 | — | **10^8.69** |

It is correct that the echinocandin fails (brain concentration 0.039 mg/L). But
**L-AmB comes out with a higher brain concentration than voriconazole** — because
even with an 11-fold lower partition coefficient its plasma concentration is
20-fold higher. This means **the model fails to reproduce the clinical received
wisdom that puts voriconazole first-line for central nervous system
aspergillosis.**

I know what mechanism is missing: liposomal amphotericin B is measured while
still trapped in the liposome and does not release free drug into the CSF. The
model handles total concentration only and does not apply a free fraction to the
brain compartment. **This is an error in the model and I record it as a point to
be fixed in the next revision.** Drug choice in central nervous system
aspergillosis must not be argued from this table.

---

## 12. Three hundred and twenty virtual patients: what predicts death

CYP2C19 genotype, timing of neutrophil recovery, treatment start day, MIC,
steroid exposure and inoculum were drawn from distributions and 320 patients
simulated (all on voriconazole).

**Overall predicted 12-week mortality 40.8%** — inside the observed IPA mortality
range.

| Neutrophil recovery day | n | Mortality |
|---|---|---|
| 5 | 33 | 14.4 % |
| 8 | 56 | 22.5 % |
| 10 | 58 | 29.5 % |
| 14 | 58 | 41.5 % |
| 21 | 56 | 51.7 % |
| 28 | 40 | 60.5 % |
| No recovery | 19 | **99.2 %** |

| CYP2C19 | n | Mortality | | Treatment start day | n | Mortality |
|---|---|---|---|---|---|---|
| UM (\*17/\*17) | 16 | 41.2 % | | 2 | 73 | 26.2 % |
| RM (\*1/\*17) | 95 | 39.3 % | | 3 | 83 | 29.3 % |
| NM (\*1/\*1) | 105 | 42.3 % | | 4 | 96 | 36.7 % |
| IM (\*1/\*2) | 86 | 37.6 % | | 6 | 43 | 69.8 % |
| PM (\*2/\*2) | 18 | **54.7 %** | | 9 | 25 | 87.7 % |

| MIC | n | Mortality | | Steroid | n | Mortality |
|---|---|---|---|---|---|---|
| 0.25 | 82 | 36.3 % | | Not used | 205 | 38.9 % |
| 0.5 | 150 | 37.4 % | | Used | 115 | 44.3 % |
| 1.0 | 59 | 38.8 % | | | | |
| 2.0 | 16 | 63.6 % | | | | |
| 8.0 | 13 | 89.8 % | | | | |

### Covariate ranking (spread in 12-week mortality, percentage points)

| Rank | Covariate | Spread |
|---|---|---|
| 1 | **Timing of neutrophil recovery** | **84.8** |
| 2 | Treatment start day | 61.5 |
| 3 | Voriconazole MIC | 53.5 |
| 4 | CYP2C19 genotype | 17.1 |
| 5 | Steroid exposure | 5.4 |

The two largest are **things that cannot be written on a prescription**. The third
(MIC) can be known from susceptibility testing, and the fourth (genotype) can be
narrowed with genotyping and TDM. That PM is the worst genotype (54.7%) is, as
seen in §7, **because of toxicity** — the efficacy is already sufficient and only
the exposure is tripled.

---

## 13. Verification: what came of writing the same 53 equations twice

`ipa_mrgsolve_model.R` (mrgsolve C++) and `ipa_reference_model.py`
(Python/scipy) were written independently of one another from the same equation
sheet. **The two implementations plus the regression tests caught 15 defects
during development.** Each defect is left as a comment on the line of the Python
source where it lived.

| # | Defect | How it showed itself |
|---|---|---|
| 1 | `sol.y` returned as a list over an interval with no point on the output grid, causing a crash | Running the steroid scenario |
| 2 | **Pure MM elimination made the PM genotype accumulate without limit** (AUC24 895, 7-fold the observed value) | Immediately, in the PK verification table |
| 3 | Posaconazole, L-AmB and echinocandin exposures inconsistent with the label | PK verification table |
| 4 | IL-6 proportional to burden without limit, giving **CRP 5405 mg/L** | Natural history scenario |
| 5 | Neutrophil recruitment driven by **lesion volume**, so that recruitment collapsed as the lesion shrank and the fungus relapsed to 10^8.6 in a patient with recovered neutrophils | d84 burden higher than d42 |
| 6 | A deterministic ODE regrowing from 10^-30 CFUe | Late relapse in eliminated arms |
| 7 | The azole EC50 coefficient rendered isavuconazole ineffective | The ISA arm identical to no treatment |
| 8 | Killing by recruited neutrophils alone left **a stable equilibrium at 10^4.3** — an immunocompetent host could not clear the fungus | The immunocompetent control failed to eliminate |
| 9 | **If there was no output point inside a dosing interval, the integration result was discarded and the previous state retained** — at dt≥2h a 1-hour infusion was deleted entirely | No steroid getting into the steroid arm |
| 10 | **When a dosing boundary fell between grid points, the state at the grid point before the boundary was carried into the next interval, so the drug mass of the intervening interval was lost** | Peak burden shifted by 0.34 log when a 0-rate dummy infusion was added, and identically so at rtol 1e-7 and 1e-11 |
| 11 | Perfusion loss linear in angioinvasion with no upper bound, so treatment on day 3 → PERF 0.94 and day 4 → 0.05, **a one-day switch** | Delay sweep |
| 12 | Drug delivery fully proportional to perfusion, so zero drug reached a thrombosed lesion | Same sweep |
| 13 | The azole EC50 sitting at 3.4-fold the achieved concentration, so a 3-fold fall in delivery dropped inhibition from 0.92 to 0.57 (a 4-day difference in marrow recovery moved mortality from 24 to 97%) | Recovery-day sweep |
| 14 | `k_grow` of 0.140/h took a 10^3 inoculum to carrying capacity in 5 days and the treatment window disappeared | Delay sweep |
| 15 | **No carrying capacity in the brain compartment, so the burden reached 10^24 CFUe**, and the log-based hazard term added ~2.6×10⁻³/h, **overestimating mortality in every failing arm** | Physically impossible values in section N |

Defect 10 was the most dangerous. Adding a **0-rate dummy infusion**, which
cannot change the physics, changed the answer by 0.34 log, and the value stayed
the same when tolerance was tightened from 1e-7 to 1e-11 — evidence of a
structural defect rather than numerical error. It is pinned down by a regression
test:

```
REGRESSION 1: a 0-rate dummy infusion must change nothing
  dt=0.5  plain 7.85449  null-breaks 7.85449  delta 2.09e-09
  dt=1.0  plain 7.85449  null-breaks 7.85449  delta 2.09e-09
  dt=2.0  plain 7.85449  null-breaks 7.85449  delta 2.09e-09
  dt=4.0  plain 7.85449  null-breaks 7.85449  delta 2.09e-09
REGRESSION 2: dose mass balance with the elimination routes switched off
  total in body after a 2-hour infusion of 6 mg/kg × 70 kg = 420.0000 mg (expected 420.0000)
  with dummy boundaries added                               = 420.0000 mg
```

### Numerical quality

| Check | Result |
|---|---|
| Occurrences of a negative state | 0 |
| Largest excursion of PERF outside [0,1] | 0.00e+00 |
| Survival function monotonically decreasing | True |
| Grid independence of logB(day 42), dt 1h vs 4h | difference 0.00e+00 |

---

## 14. Calibration targets, and what fitted well / what did not

| Target | Observed | Model | Verdict |
|---|---|---|---|
| Voriconazole 4 mg/kg IV NM trough | 1–3 mg/L | 1.13 | ✅ |
| Isavuconazole 200 mg qd AUC24 | about 100 mg·h/L | 92.5 | ✅ |
| L-AmB 3 mg/kg Cmax / AUC24 | about 83 / 555 | 73.3 / 552 | ✅ |
| Voriconazole dose-exposure supraproportionality | present | 3.04-fold from 4→8 mg/kg | ✅ |
| CYP2C19 exposure spread | 4–10-fold | 5.4-fold | ✅ |
| L-AmB 10 mg/kg has no benefit over 3 mg/kg and only more toxicity | AmBiLoad | 11.7% vs 16.5% | ✅ |
| Voriconazole-tacrolimus interaction | about 3-fold | 2.63-fold | ✅ |
| Isavuconazole-tacrolimus | about 1.4-fold | 1.45-fold | ✅ |
| Reduced GM sensitivity under mould-active therapy | Marr 2005 | a GM-negative 10^5.7 infection | ✅ |
| 12-week mortality of the virtual population | 30–60% | 40.8% | ✅ |
| EUCAST voriconazole breakpoint | 1 mg/L | cliff between 2 and 4 | ⚠️ permissive |
| Voriconazole > L-AmB in the central nervous system | clinical received wisdom | brain concentration comes out higher for L-AmB | ❌ not reproduced |
| Non-inferiority among the three azoles | SECURE and others | **encoded by the model, not predicted by it** | — |

The last row matters. The EC50 coefficients of the three triazoles were
**calibrated to be equipotent against wild-type strains at the licensed doses**.
That is, under wild-type, licensed-dose conditions a comparison between azoles is
an input to the model, not a prediction of it. What the model does **predict** is
how that equipotency breaks down as MIC, genotype, perfusion, site of infection
and timing of marrow recovery move away from the trial population.

---

## 15. Using the files

```bash
# Render the map
dot -Tsvg ipa_qsp_model.dot -o ipa_qsp_model.svg
dot -Tpng -Gdpi=150 ipa_qsp_model.dot -o ipa_qsp_model.png

# Regenerate every number above (about 30 minutes)
python3 ipa_reference_model.py

# mrgsolve model (23 scenarios + 9 analyses)
Rscript ipa_mrgsolve_model.R

# Dashboard
R -e 'shiny::runApp("ipa_shiny_app.R")'
```

## 16. Limitations of the model (these set the scope of use, so read them)

- **The comparison between azoles is an input, not an output** (§14).
- **Do not use this model for central nervous system drug choice** (§11). The free
  fraction of liposome-bound drug is not reflected in the brain compartment.
- **Echinocandin monotherapy** has no controlled trial as first-line therapy in
  neutropenia. That arm (S11) is extrapolation beyond the evidence.
- **The position of the treatment-delay threshold (6 days) is parameter-dependent.**
  What can be carried across is the shape — that a threshold exists — not the
  date.
- **This is a deterministic model.** Stochastic clearance in the neighbourhood of
  1 CFU is not represented, and only regrowth from fractional organisms is
  prevented, by a growth gate.
- **Immune reconstitution inflammatory syndrome (IRIS) on neutrophil recovery is
  not modelled.** The model treats returning neutrophils as purely beneficial,
  and the chest CT does not always look that way (air-crescent sign).
- **`PERF` could not be calibrated against human data.** The evidence is set out
  in the last section of [`ipa_references.md`](../../../invasive-aspergillosis/ipa_references.md),
  "Points where no literature exists".

---

This is a QSP model for education and research. **Do not use it for clinical
decision-making.**
