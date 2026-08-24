# Cutaneous T-Cell Lymphoma (Mycosis Fungoides / Sézary Syndrome) — QSP Model
### Cutaneous T-Cell Lymphoma (Mycosis Fungoides / Sézary Syndrome) · Quantitative Systems Pharmacology

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map | [`ctcl_qsp_model.dot`](ctcl_qsp_model.dot) · [SVG](ctcl_qsp_model.svg) · [PNG](ctcl_qsp_model.png) — 187 nodes / 273 edges / 21 clusters |
| ⚙️ mrgsolve model | [`ctcl_mrgsolve_model.R`](../../../cutaneous-t-cell-lymphoma/ctcl_mrgsolve_model.R) — 67 ODEs, time unit = day, 24 scenarios + 6 structural experiments |
| 📊 Shiny dashboard | [`ctcl_shiny_app.R`](ctcl_shiny_app.R) — 13 tabs |
| 📚 References | [`ctcl_references.md`](ctcl_references.md) — 181 PubMed articles (checked in full via E-utilities) |

---

## The Organising Thesis

**One clone, three compartments, two host variables.**

```
                      ┌──────────────┐
        J_home ──────►│  SKIN        │◄── mSWAT reads only this box
     (CCR4·CLA)       │  T_RM phenotype │
                      └──────┬───────┘
                             │ J_egress = k·(1−CD69)·CCR7
                      ┌──────▼───────┐
                      │  BLOOD       │◄── B-score reads only this box
                      │  T_CM phenotype │
                      └──────┬───────┘
                      ┌──────▼───────┐
                      │  LYMPH NODE  │◄── N-score reads only this box
                      └──────────────┘

        GLOBAL response = AND of the three boxes above   ← Olsen 2011
```

And on the host side there are **two mutually independent** state variables.

```
   Variable A · how much was killed        NSK NSKR NTR NBL NBLR NLN NVS
   Variable B · who is left to keep watch  DCA E8SK E8BL NKB TREG
```

Every **ORR reads only variable A**, and every **TTNT/DOR is set by variable
B**. Cytotoxic chemotherapy moves the two **in opposite directions**, while
interferon and extracorporeal photopheresis move **only B**. That is how "a
drug with a high response rate that relapses quickly" and "a drug with a low
response rate that lasts a long time" both fall out of the same two lines of
algebra.

### Three constants that were not fitted

| Constant | Value | Source | What it does in this file |
|---|---|---|---|
| `ABCSK` | 0.157 | Antibody tissue-partition coefficient (skin) | Skin interstitial concentration of every antibody = 0.157 x plasma |
| `ABCLN` | 0.085 | Antibody tissue-partition coefficient (lymph node) | Lymph-node exposure |
| `NKRSK` | 0.12 | Lesional-skin/blood NK density ratio | The compartmental gap in ADCC |

Fixing these three turns the following from **an assertion into a computed
result**.

- Mogamulizumab's **skin CCR4 occupancy is already close to 1**. Even
  0.157 x the plasma concentration is far above the KD (0.02 µg/mL).
  **So what is scarce in the skin is not the drug but the effector cell.**
  Tripling the dose does not move the skin (structural experiment E1
  reproduces exactly this).
- If so, this drug's skin effect cannot be ADCC. Only two pathways remain in
  the file — (a) blocking the **CCR4 homing route** to cut off the blood as
  a *source*, and (b) **depleting the Treg, the cell that expresses CCR4
  most highly**, first, which releases CD8 surveillance. Both are edges that
  were already on the map, not ones added to fit the skin response rate.
- (b) **is the same event as mogamulizumab-associated rash (MAR)**. So the
  model predicts that rash and response go together.
- Antibody-drug conjugates do not have this gap. **Because MMAE does not
  require an effector cell**, brentuximab vedotin's only compartmental
  penalty is the 0.157 itself. For the same reason, this drug is relatively
  strong in the skin, whereas mogamulizumab is relatively strong in the
  blood.

### Stage is an equilibrium point, not an initial value

TNMB stage is a burden the patient has carried for years. Dropping it onto a
naive host and starting integration makes the "response" of the first few
weeks nothing but a transient. So the patient is built in four stages —
①equilibrate the host with the clone held fixed, ②release the clone and let
the skin/blood/lymph-node distribution settle on its own, then return **only
the total burden** to the staged value, ③fix it again to let lesion
morphology settle, ④back-solve the immune fitness `IMMF` at which the clone
stalls at that burden, then apply **the same surveillance deficit across all
stages**.

The result is that **the distinction between MF and SS becomes an output,
not an input.** Lowering only `CD69` from 0.85 to 0.10 in the same equations
moves the circulating clone from B0 to B2 (structural experiment E6).

---

## Compartments and ODE Structure (67 states)

| Group | State variables |
|---|---|
| Antibody PK | `MOG1` `MOG2` `BV1` `BV2` `MMAE1` `MMAE2` `ALEM` |
| Small-molecule PK | `ROM1` `ROM2` `VORD` `VOR` `BEXD` `BEX` `IFND` `IFN` `MTXD` `MTX` `GEM` |
| Local/physical-therapy exposure | `SDT` `TSEBC` `STER` `ABX` `ECPD` |
| Clones (sensitive/resistant/transformed) | `NSK` `NSKR` `NTR` `NBL` `NBLR` `NLN` `NVS` |
| Antigen expression (replicator dynamics) | `FCCR4` `CD30E` |
| Surveillance | `E8SK` `E8BL` `NKB` `TREG` `DCA` |
| Cytokines/biomarkers | `TH2` `IL31` `IL10` `IFNG` `TARC` `SIL2R` `LDH` |
| Lesion morphology | `APAT` `APLQ` `ATUM` |
| Barrier, microbiome, pruritus | `BARR` `SAUR` `SAG` `PRUR` `SENS` |
| Haematology/toxicity | `PROL` `TR1` `TR2` `TR3` `CIRC` `PLT` `PN` `TSH` `FT4` `TG` `MAR` `CD4N` |
| Risk integrals | `HZI` `HZD` `ADA` |

Lesion area has a **different time constant for each morphology** — 8 days
for patch, 25 for plaque, 12 for tumour. The reason a 4-week confirmation is
needed for response assessment comes out of the kinetics, not a rule.

---

## Why Resistance Is Split Into Two Numbers

The kill fraction left over in the drug-resistant subclones (`NSKR`, `NBLR`)
splits into two numbers.

```
   RESK = 0.12   only 12% of drug killing gets through   (8x resistance)
   RESI = 0.50   50% of immune killing still gets through (2x resistance)
```

This **4-fold asymmetry** is the only device in this model that separates
ORR from TTNT. A therapy that kills only with the drug barely touches the
resistant subpopulation, so its response is shallow and short-lived, while a
therapy that restores surveillance also reaches the resistant subpopulation,
so its response is slower but lasts longer. Total-skin electron beam (TSEB)
clearing the skin while **also burning out the skin-resident CD8** creates a
relapse for the same reason (structural experiment E5).

---

## Superantigen Feedback — a Pathway That Changes the Disease Without Killing the Clone

```
   Lesion density up -> barrier BARR down -> antimicrobial peptides down -> S. aureus SAUR up
        -> superantigen SAG up -> Vbeta-restricted clonal stimulation (proliferation x sag) and IL-31 up
        -> pruritus up -> scratching -> barrier down  (loop closes)
```

Antibiotics lower only `SAUR`. They do not touch the clone. Yet mSWAT and NRS
still fall (structural experiment E4). This is what anti-staphylococcal
treatment does in this file, and it is also, at the other end of the same
loop, why infection is the **leading cause of death** in this disease.

---

## Scenarios (24)

| # | Scenario | What this scenario asks |
|---|---|---|
| 01 | Untreated stage IB MF | Natural history |
| 02 | nbUVB + topical steroid | Depth limit of an epidermis-only treatment |
| 03 | Chlormethine gel | Surface alkylating agent |
| 04 | Low-dose TSEB, 12 Gy | Clears the skin and burns out surveillance |
| 05-07 | Bexarotene / interferon / combination | Cell-intrinsic vs. immunomodulatory |
| 08-09 | Vorinostat / romidepsin | HDAC inhibitors |
| 10-11 | Mogamulizumab (MF / Sézary) | The core of the compartmental gap |
| 12 | Vorinostat (Sézary) | MAVORIC control arm |
| 13 | ECP + interferon | A therapy that raises only the surveillance variable |
| 14 | Low-dose alemtuzumab | Clears the blood and clears the host too |
| 15-16 | Brentuximab (CD30 high / low) | CD30 heterogeneity |
| 17 | Romidepsin x2 -> brentuximab | Sequencing rationale via HDACi-induced CD30 |
| 18 | Gemcitabine x6 | The classic high-ORR, short-TTNT pattern |
| 19 | Low-dose methotrexate | |
| 20 | Large-cell-transformed MF + brentuximab | |
| 21 | Anti-staphylococcal antibiotics alone | A treatment that does not kill the clone |
| 22 | Mogamulizumab + ECP | Moving both variables at once |
| 23 | CCR4-low clone + mogamulizumab | Antigen loss |
| 24 | Mogamulizumab 3 mg/kg | Dose is not the limiting factor |

### Structural Experiments (E1-E6)

All run via `Sys.setenv(CTCL_RUN = "1"); source("ctcl_mrgsolve_model.R")`.

| Experiment | Question asked |
|---|---|
| **E1** | Is the skin gap an exposure gap or an effector-cell gap (3x dose vs. `NKRSK` -> 1) |
| **E2** | Same patient, four regimens — how do maximum response depth and TTNT diverge |
| **E3** | Do HDAC inhibitors raise CD30 to make a subsequent ADC stronger |
| **E4** | Do antibiotics change the disease without touching the clone |
| **E5** | Are TSEB's complete response and its subsequent relapse the same event |
| **E6** | Does CD69 alone turn MF into Sézary syndrome |

---

## Distinguishing Inputs from Outputs

**Inputs** — tissue residency (`CD69`, `CCR7`), `CLA`, `CCR4`, `CD30`,
drug-resistant fraction `FRES`, clonal drug sensitivity `SENSF`, NK density
`NKF`, immune fitness `IMMF`, transformed-subclone proliferation multiplier
`GTRF`, staged burden, body size, and dosing regimen.

**Outputs** — mSWAT and its patch/plaque/tumour decomposition, Sézary cell
count and B-score, N-score, GLOBAL response via the Olsen AND rule, PFS,
TTNT, pruritus NRS, serum TARC/sIL-2R/LDH, ANC/platelets, MMAE neuropathy,
bexarotene's central hypothyroidism and hypertriglyceridaemia, CD4 count,
infection risk and disease risk, and **whether the patient is MF or Sézary**.

---

## What Was Fitted and What Was Not

**Exactly three** potency parameters were fitted. Each was matched to one
published response rate, using a 90-patient virtual cohort (previously
treated MF/SS) simulated for 78 weeks and scored with the Olsen AND rule.

| Parameter | Value | Anchor | Model output |
|---|---|---|---|
| `EVOR` | 0.0394 | MAVORIC vorinostat **skin** ORR 15.8% | Fitted |
| `KADCC` | 0.0482 | MAVORIC mogamulizumab **blood** ORR 67.7% | 64.7% |
| `KBVKILL` | 0.5178 | ALCANZA brentuximab skin ORR (CD30>=10%) 56.3% | 56.9% |

The remaining potency terms (bexarotene, romidepsin, methotrexate,
gemcitabine, interferon, ECP, alemtuzumab, phototherapy, TSEB, chlormethine)
are **unfitted order-of-magnitude estimates**, and the response rates those
arms produce are predictions rather than calibrations. In particular,
**mogamulizumab's skin response is a prediction** — only the blood arm was
fitted.

## Limitations

- This is a semi-quantitative model for education and research. It must not
  be used for clinical decision-making, prescribing, or regulatory
  submission.
- **A known defect is recorded rather than hidden — the blood-compartment
  response rule is not yet reliable over long-term observation.** A patient
  with a baseline Sézary count of 250-1000/µL (B1) can drop below half of
  baseline simply through the slow settling of skin<->blood redistribution,
  so **running the untreated arm out to 78 weeks registers a formal blood
  response in roughly half of B1 patients.** This model's blood response
  rate should therefore be read only for B2 (>=1000/µL) patients, and only
  against a concurrent untreated arm. This behaviour does not appear in the
  skin, lymph-node, or GLOBAL scoring. **Because of this defect, this README
  does not include a virtual-cohort response-rate table.**
- **This container has no R installed, so the mrgsolve file has never been
  executed.** All 67 ODEs were developed, calibrated, and checked in an
  independent Python/scipy implementation and then transcribed; the R file
  has only been statically checked (every parameter, compartment, and
  CAPTURE name resolves, and every compartment has a derivative). Its first
  run should be treated as a **port test**, not a validated run.
- The virtual cohort's parameter distributions approximate literature
  values and are not fitted to real patient data.
- Daily oral agents (vorinostat, bexarotene, etc.) are represented as a
  zero-order daily input, because with a half-life of 2-7 hours and a PD
  time constant of several weeks, mean exposure dominates; this
  approximation is applied identically in both the R file and the Python
  prototype.
- Because large-cell transformation is a seed that grows exponentially in a
  deterministic ODE, the transformation probability itself is kept low and
  the aggressiveness of the transformed subclone is set explicitly via the
  patient covariate `GTRF`.
- Survival is represented only as hazard integrals (`HZD`, `HZI`), not as a
  competing-risks model.

## Disclaimer

This model is a **qualitative/semi-quantitative QSP model for education and
research purposes**, built from the public literature. It has not been
independently validated or certified and must not be used directly for real
clinical decision-making, prescribing, or regulatory submission.
