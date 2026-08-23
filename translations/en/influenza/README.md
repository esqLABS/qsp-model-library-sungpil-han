# Influenza A — QSP Model

> **All a drug can eliminate is virus that has not yet happened.**
> "Start antivirals within 48 hours" is usually taught as a rule about the
> *clock*. This model recasts it as a rule about **a quantity that is being
> depleted**.

The viral AUC that a drug started at time `t_rx` can reduce is bounded above
**exactly** by the following value.

```
        R(t_rx) = ∫_{t_rx}^{∞} ( log10 V(t) − LOD )⁺ dt      [residual AUC]
```

`R` is **a quantity determined solely by the untreated course.** No drug, no
potency, no mechanism enters into it. No molecule can cross this line.
Everything else in this repository is a calculation of (a) how fast `R`
declines and (b) how much of it each operator takes.

In the calibrated adult (`A1`):

| Dosing time | Target cells remaining | R(t) | % of total AUC |
|---|---|---|---|
| At symptom onset | 99.9 % | 18.61 | **97.9 %** |
| +12 h | 85.7 % | 16.99 | 89.4 % |
| **+24 h** (CAPSTONE-1 median enrolment time) | 0.1 % | 14.27 | **75.1 %** |
| +48 h (end of the licensed window) | 3.5 % | 8.65 | **45.5 %** |
| +96 h | 13.2 % | 1.96 | 10.3 % |

And at the median enrolment time, the calibrated baloxavir takes **65.6 %**
of what remains. So the **headroom** in the two directions is as follows.

```
    Potency headroom  (Emax → 1)               +0.09 log10·d
    Timing headroom   (dosed at symptom onset)  +4.34 log10·d      ← 48-fold
```

**Nearly all of the remaining benefit sits on an axis that a better molecule
cannot buy.** This is the quantitative content of "treat early." Late
treatment is not weak because the drug is weak — it is asked to act on a
quantity that has already mostly been spent.

---

## Operator Classification

Every treatment is classified by **which term** of the replication loop it
touches. This classification matters because each class leaves a
**differently shaped** signature on the time axis, and the classes are not
interchangeable (`A5`).

| Operator class | Term touched | Examples |
|---|---|---|
| **ENTRY** | β | mucosal IgA, neutralising monoclonal antibodies, (weakly) NAIs |
| **TRANSCRIPTION** | E → I | baloxavir (PA endonuclease, cap-snatching) |
| **PRODUCTION / RELEASE** | p | oseltamivir · zanamivir · peramivir, baloxavir, favipiravir |
| **MUTAGENESIS** | fraction of infectious progeny | favipiravir (lethal mutagenesis) |
| **VIRION CLEARANCE** | c | monoclonal antibodies, convalescent plasma |
| **TARGET PROTECTION** | T → R | type I·III interferon, the ISG antiviral state |
| **INFECTED-CELL DEATH** | δ | CD8 CTL, NK — **the strongest operator in the model, and no licensed anti-influenza drug works this way** |
| **IMMUNOPATHOLOGY DAMPING** | suppresses the terms above | corticosteroids (harmful direction here) |

---

## Deliverables

| File | Contents |
|---|---|
| [`flu_qsp_model.dot`](../../../influenza/flu_qsp_model.dot) · [SVG](../../../influenza/flu_qsp_model.svg) · [PNG](../../../influenza/flu_qsp_model.png) | Mechanistic map — 14 clusters / 145 nodes / 217 edges. Treatment nodes are coloured by the operator class above |
| [`flu_mrgsolve_model.R`](../../../influenza/flu_mrgsolve_model.R) | 50-ODE mrgsolve model (virus·epithelium 18 · immunity 9 · clinical 3 · drug PK 17 · accounting 3), 14 scenarios, 14 analysis functions |
| [`flu_shiny_app.R`](../../../influenza/flu_shiny_app.R) | 10-tab interactive dashboard (tab 2 is the point of this model) |
| [`flu_references.md`](flu_references.md) | **97 references, verified by PubMed lookup**, each stating which part of the model it supports, section by section |
| [`flu_reference_check.py`](../../../influenza/flu_reference_check.py) | **An independent numpy/scipy port of the same equation system.** Every number in this README was computed here, and it can be reproduced without R |
| [`flu_reference_output.txt`](../../../influenza/flu_reference_output.txt) | The full output of the script above (A0–A13), verbatim |

This repository's build environment has no R toolchain. Rather than publish
an ODE model that has never been integrated, the same system was ported
once more and integrated with scipy LSODA, and every number below is the
output of `python3 flu_reference_check.py`. If the two ports disagree, one
of them is wrong.

```bash
python3 flu_reference_check.py            # everything (A0–A13)
python3 flu_reference_check.py --only A1  # key result only
python3 flu_reference_check.py --list
```

---

## 0. Calibration Check (A0)

Natural course, no treatment, no drug:

| Item | Model | Target |
|---|---|---|
| Upper-airway peak titre | 6.62 log₁₀ TCID50/mL | 6.0–7.0 |
| Time to peak | 51.5 h p.i. | 48–72 h |
| Shedding duration (>LOD) | **6.3 d** | 4.0–5.5 d ← **long, reported as is** |
| Symptom onset | 24.2 h p.i. | 24–48 h |
| Peak symptom score | 12.4 / 21 | 10–15 |
| Peak temperature | 38.99 °C | 38.5–39.5 |
| Fever duration | 67 h | 48–96 h |
| Peak IL-6 | 60 pg/mL | 20–100 |
| **Minimum susceptible epithelium** | **0.1 % of T₀** | <20 % (target-cell limited) |
| Lower-airway peak titre | 3.87 log₁₀ | 2–4 (uncomplicated) |
| Peak CD8 expansion | 99 × naive | 20–100 × |

The most important row is the minimum susceptible epithelium. The fact that
the infection is **target-cell limited** is the condition that pins A1's
boundary.

---

## 1. Trial Ledger (A3) — Model vs Published Values

CAPSTONE-1 (Hayden 2018, NEJM 379:913, [PMID 30184455](https://pubmed.ncbi.nlm.nih.gov/30184455/)):

| Endpoint | Placebo (published / model) | Oseltamivir (published / model) | Baloxavir (published / model) |
|---|---|---|---|
| Time to alleviation of symptoms (TTAS) | 80.2 / **78.7** h | 53.8 / **63.5** h | 53.7 / **44.0** h |
| Cessation of viral shedding | 96.0 / **118.0** h | 72.0 / **85.2** h | 24.0 / **47.0** h |
| Day-2 titre change | −1.3 / **−1.32** log₁₀ | −2.8 / **−2.59** log₁₀ | −4.8 / **−4.76** log₁₀ |

**The virological columns fit well, and the symptom column is the weakest
part of this model.** A8 answers why — and that answer is the most testable
claim in this model.

All 14 scenarios (peramivir · favipiravir · antibody · combination · timing
variants · immunocompromised host · resistance profile · vaccinated host)
are in `flu_reference_output.txt`.

---

## 2. Resistance Is Competitive Release (A6)

The resistant subpopulation is **not seeded.** It arises from a mutation
term (μ = 2.5×10⁻⁵) acting on the wild type from the very first replication
cycle, and it grows when the drug **removes the competitor and hands over
target cells the wild type would otherwise have taken.**

Peak resistant-strain titre (log₁₀) from a dosing-time × potency grid:

| Dosing (after symptom onset) | T remaining | Emax 0 | 0.9 | 0.99 | 0.999 | 0.9999 | 0.999999 |
|---|---|---|---|---|---|---|---|
| **0 h** | 99.95 % | 2.58 | **6.19** | 6.19 | 6.19 | 6.19 | 6.19 |
| 12 h | 85.73 % | 2.58 | 2.15 | 2.75 | 2.88 | 2.90 | 2.90 |
| 24 h | 0.10 % | 2.58 | 2.53 | 2.53 | 2.53 | 2.53 | 2.53 |
| 48 h | 3.54 % | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 | 1.21 |
| 72 h | 7.10 % | −0.22 | −0.22 | −0.22 | −0.22 | −0.22 | −0.22 |

Two things that emerged without being put in:

1. **Resistance does not simply decline along the potency axis.** At Emax =
   0, the resistant strain loses to the more fit wild type. As the drug
   gets better, the resistant strain stops losing — because the competitor
   is removed and target cells are handed over. **Resistance selection is
   the result of the drug working, not of it failing.**
2. **The time axis has an easily missed precondition: release requires a
   field to release into.** Late dosing selects for almost no resistance —
   not because the drug is weak, but because the wild type has already
   consumed the epithelium and there is nothing left to hand over.
   **Resistance emergence and clinical benefit share the same
   precondition, and the same policy of "treat earlier" increases both.**

---

## 3. What the Symptom Endpoint Requires of the Model (A8)

CAPSTONE-1's virological endpoints are reproduced by viral kinetics alone.
**The symptom endpoint is not, and the reason lies in structure, not in
numbers.**

If the symptom score is driven purely by cytokine status (WVIR = 0), its
magnitude and timing are already set at the interferon peak, and that peak
has already passed **before** a drug given at the median enrolment time can
act. In that structure, **no antiviral, at any potency, can shorten the
illness.**

`WVIR` is the fraction of symptom drive that tracks **the titre at that
instant**:

| WVIR | Onset (h) | Placebo TTAS | Oseltamivir benefit | Baloxavir benefit |
|---|---|---|---|---|
| **0.00** | 42.8 | 55.5 h | **0.0 h** | **0.8 h** |
| 0.20 | 33.2 | 66.2 h | 4.0 h | 8.8 h |
| 0.40 | 27.0 | 74.2 h | 9.8 h | 21.0 h |
| **0.60** (calibrated value) | 24.2 | 78.7 h | 15.2 h | 34.7 h |
| 0.80 | 23.0 | 82.5 h | 22.8 h | 46.5 h |
| 1.00 | 22.0 | 86.2 h | 32.7 h | 54.2 h |

CAPSTONE-1 reported a 26.5-hour benefit for both drugs. **The published
symptom endpoint cannot be explained by a symptom model driven by cytokines
alone.** WVIR = 0.60 is the value paid to match that, and it is not a
fitting convenience but a **falsifiable claim** — the claim that a
substantial part of what the patient feels tracks *today's replication*,
not the magnitude of a cytokine wave that has already passed.

---

## 4. Operator Decomposition (A5)

Result of turning on each operator **alone at 95%** from t_rx (no PK,
matched comparison):

| Operator | TTAS (h) | Shedding (h) | day-2 (log) | AUC (log·d) | Epithelial loss |
|---|---|---|---|---|---|
| PRODUCTION / RELEASE (p) | 66.2 | 88.5 | −2.63 | 13.89 | 85.1 % |
| TRANSCRIPTION (E→I) | 73.5 | 114.5 | −1.52 | 18.23 | 85.7 % |
| ENTRY (β) | 78.2 | 117.5 | −1.34 | 18.92 | 84.1 % |
| VIRION CLEARANCE (c) | 65.8 | 88.0 | −2.70 | 13.24 | 84.2 % |
| **INFECTED-CELL DEATH (δ)** | **32.2** | **29.0** | **−4.16** | **9.56** | 72.4 % |
| TARGET PROTECTION (T→R) | 78.5 | 117.5 | −1.33 | 18.94 | 84.7 % |
| (none: placebo) | 78.7 | 118.0 | −1.32 | 19.00 | 85.6 % |

- ENTRY and TARGET PROTECTION barely move the 24-hour decline. **Cells
  already infected** at the moment the operator turns on keep producing
  regardless. These operators block only the *next* generation.
- PRODUCTION/RELEASE and VIRION CLEARANCE drop the titre by nearly the same
  magnitude. At quasi-steady state the titre is `p·I/c`, so these two are
  the numerator and denominator of the same fraction.
- **INFECTED-CELL DEATH is overwhelmingly the strongest.** Only this one
  *removes* the source; the rest merely modulate it. No licensed
  anti-influenza drug works this way — the immune system does — which is
  why the natural course is more decisive than any drug that acts on it.
- **TRANSCRIPTION is weaker than PRODUCTION at the same 95%.** This is not
  obvious: it turns away only cells not yet mature, whereas blocking
  production also silences cells that are already mature. Baloxavir beats
  NAIs in the ledger **because it achieves higher potency while acting on
  the weaker point**, not because its point of action is better.

---

## 5. Other Findings

- **Prophylactic dosing is the same molecule used on the other side of the
  peak (A10).** Dosed within 12 hours of infection, the model records no
  illness at all; dosed at 48 hours, it shortens the illness. What
  separates the two regimes is the sign of `(peak time − dosing time)`.
  What BLOCKSTONE ([PMID 32640124](https://pubmed.ncbi.nlm.nih.gov/32640124/))
  reported was also a reduction in **probability of illness** (1.9 % vs
  13.6 %), not duration.
- **Corticosteroids (A9).** The direction is right and the magnitude is
  small. Steroids lower the symptom score while increasing viral AUC
  (12.96 → 13.07) and shedding — the direction of the observational signal
  is reproduced by mechanism alone, without a single mortality term. But
  **the small magnitude is a clue about this calibration, not reassurance
  about the drug**: in this parameterisation, neither IFN-mediated target
  protection nor CD8 killing is the rate-limiting step in clearance, so
  suppressing them costs little. In a host where CD8 is rate-limiting
  (elderly · immunosuppressed — exactly the population the observational
  signal came from), the cost is far larger.
- **In the immunocompromised host, the window does not close (A12).** In
  the normal host, the drug competes with an immune response that would
  have cleared the virus anyway, so its marginal value vanishes over time
  (9.37 → 1.90 → 0.00 log·d at 24 / 96 / 168 h). Remove that response, and
  the decay slows (11.89 → 4.95 → 0.08).
- **Lower airway and bacterial superinfection (A11).** In the uncomplicated
  adult, bacterial burden rises about 4 log but does not cross the
  threshold **in any treatment arm.** Because the threshold itself is a
  stipulated, not calibrated, value, this module can only **rank**
  treatment arms and cannot say who develops pneumonia.

---

## 6. What This Model Does Not Reproduce

Points where the model disagrees with the literature are numbered and
reported, not hidden.

| # | Where | Content |
|---|---|---|
| **1** | A2 | With a Hill slope of 1 (Michaelis), no Emax can produce CAPSTONE-1's 24-hour titre decline (−4.8 log₁₀), since the floor of the residual production fraction is EC50/(C+EC50), which saturates near −3.3 log₁₀. The trial data is evidence about the **shape of the concentration-response curve**, not about potency. In addition, the calibrated in-vivo EC50 is 13-fold lower than the free-drug in-vitro value. |
| **2** | A12 | The prediction that late treatment retains value in the immunocompromised host. It agrees with clinical practice, but no randomised trial has tested it — **it is unfalsified, not validated.** |
| **3** | A8 | If symptoms track the virus at all, baloxavir should show a larger symptom benefit than oseltamivir. CAPSTONE-1 could not distinguish them (53.7 vs 53.8 hours) despite a 100-fold difference in day-2 titre. |
| **4** | A6 | Because the model is deterministic, a resistant lineage at the 10⁻⁴ level "exists" in every virtual patient. Actual I38T emergence is a stochastic establishment event that occurs in about 10 % of adults. The grid above should be read as **selection pressure**, not as a prediction of viral load. |
| **5** | A7 | **The order comes out reversed.** In this model, competitive release is dominated almost entirely by `T(t_rx)`, and **anything that slows the wild type leaves more field** — including pre-existing immunity. So pre-immune adults select the most resistance and children the least (resistant AUC 17.90 vs 3.19). This runs opposite to the observation that I38T is more common in children, and is the clearest sign that the well-mixed epithelium assumption is wrong for this question. The **immunocompromised host** (resistant AUC 29.33), however, is where both mechanisms point the same way, and the model gets it right. |

---

## 7. Structural Notes (mines easy to step on when editing the model)

- **A mixed-unit convention.** Cells are counted in absolute numbers, virus
  in titre per mL of wash fluid, and β absorbs the sample volume. As a
  result, (a) the absorption term `−β·T·V` has units of cells/day and
  cannot be subtracted from a titre (it is dropped by standard convention,
  absorbed into `c`), and (b) **β and p cannot be individually identified
  from titre data alone** — the growth rate fixes only the product
  `β·T₀·p`, and `p` is set by the *scale* of the titre. Neither should be
  cited on its own.
- **Blocking transcription does not freeze cells in eclipse.** Cells always
  leave eclipse via `KECL`, and the drug determines *where* that cell goes
  (productive vs abortively infected). Wiring it as blocking the exit
  instead creates an eclipse reservoir that keeps supplying productive
  cells for days after the drug starts working, and then **no potency at
  all** can reproduce baloxavir's 24-hour decline.
- **The extinction floor.** A deterministic ODE lets a lineage fall to
  10⁻²⁰ cells and then revive once target cells recover. Without a term
  that smoothly sends release to zero near 1 infected cell, the epithelial
  regeneration term manufactures a spurious second wave.

---

## 8. What This Model Is Not

- It is a deterministic, well-mixed model. Influenza replicates in a
  spatially structured mucosa, and stochastic extinction and the spatial
  separation of lineages matter **precisely for the resistance question
  this model asks** (discrepancies 4 and 5).
- Symptoms are a single-compartment composite score whose coupling
  coefficient to viral load (`WVIR`) is calibrated. This parameter is the
  **hardest to defend and most important** value in this model, because
  every symptom endpoint in every influenza trial hangs on it.
- The bacterial superinfection threshold is stipulated, not calibrated
  (A11).
- The PD of the monoclonal antibody, favipiravir, and peramivir are assumed
  values. Only oseltamivir and baloxavir are calibrated to clinical data.

---

> This is a QSP model for educational and research purposes. It has not
> been independently validated and must not be used for clinical
> decision-making, prescribing, or regulatory submission.
