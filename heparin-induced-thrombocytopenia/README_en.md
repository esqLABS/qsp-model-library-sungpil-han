# Heparin-Induced Thrombocytopenia (HIT) — QSP Model
### Heparin-Induced Thrombocytopenia · including autoimmune HIT (aHIT), spontaneous HIT syndrome and VITT

<a href="hit_qsp_model_en.svg"><img src="hit_qsp_model_en.png" width="720" alt="HIT QSP mechanistic map"></a>

| Deliverable | File | Scale |
|---|---|---|
| 🗺️ Mechanistic map | [`hit_qsp_model_en.dot`](hit_qsp_model_en.dot) · [`.svg`](hit_qsp_model_en.svg) · [`.png`](hit_qsp_model_en.png) | **226 nodes · 22 clusters · 330+ edges** |
| ⚙️ mrgsolve ODE model | [`hit_mrgsolve_model.R`](hit_mrgsolve_model.R) | **59 ODEs · 31 scenarios · 14 diagnostic analyses** |
| 📊 Shiny dashboard | [`hit_shiny_app.R`](hit_shiny_app.R) | **10 tabs** |
| 📚 References | [`hit_references.md`](hit_references.md) | **159 PubMed citations** |

---

## 1. The organising idea

> **HIT is not a drug toxicity — it is an immune complex disease in which
> the drug forms one half of the antigen.**

Platelet factor 4 (PF4/CXCL4) is a 7.8 kDa tetramer with a net charge of
about +20, and heparin is a polyanion. At the molar ratio where charge is
nearly neutralised, the two do not simply bind — they polymerise into
**ultralarge complexes (ULC, >670 kDa).** Only the ULC is antigenic,
because only the ULC presents neoepitopes at a surface density that
multivalent IgG can bind, and only multivalent IgG can cross-link FcγRIIa.

This fact dictates a consequence that governs the entire model: **the
antigen is not a monotonic function of heparin dose.** So the model enters
the polyanion into the equations three times, in three different roles.

| Role | Model variable | Consequence |
|---|---|---|
| **① Substrate** | `HEPB` — the concentration of chains long enough (≥11–14 saccharide units) to bridge two PF4 tetramers | more chains means more antigen → **a monotonically increasing dose-response within a given drug** |
| **② Competitor** | `HEPT` — total chain concentration weighted by PF4-binding capacity. Short chains cannot bridge but still capture and steal PF4 | when `HEPT` ≫ PF4, each tetramer claims a separate chain and bridging fails → **heparin excess** |
| **③ PF4 mobiliser** | long chains displace PF4 from endothelial heparan sulfate (`KMOB`, `KHSPG`) | within minutes of a UFH bolus, plasma PF4 rises 10–20-fold → **the drug manufactures its own antigen partner** |

Because ① and ③ favour long chains while ② penalises high chain molarity,
**the clinical risk ordering is derived from chain chemistry, not
asserted.**

```
UFH therapeutic dose  ≫  LMWH therapeutic dose  >  prophylactic dose (LMWH/UFH)  ≫  danaparoid  ≫  fondaparinux (structural null)
```

And the single experimental fact any HIT model must reproduce — **platelets
activate at 0.1–0.3 U/mL heparin but the response is abolished at
100 U/mL** — is likewise not coded as a rule; it simply falls out of
evaluating term ② of the same equation at two points 250-fold apart on the
heparin axis. This is the confirmatory step of the serotonin release assay
(SRA).

---

## 2. Three serial layers and one self-closing loop

| Layer | Compartments | Contents |
|---|---|---|
| **LAYER 1 · Antigen** | `PF4P` `PF4EC` `PF4DNA` `ULC` | PF4 + polyanion → ULC. **Only stopping heparin** breaks this layer |
| **LAYER 2 · Immune** | `BCELL` `PBLAST` `LLPC` `MEMB` `IGGP` `IGGN` | ULC → B-cell clone → plasma cell → anti-PF4/heparin IgG |
| **LAYER 3 · Effector** | `IC` `PLT` `PLTA` `MP` `MOA` `TF` `NET` `THR` | IgG·ULC → FcγRIIa cross-linking → platelet activation and consumption, PS+ microparticles, monocyte tissue factor, a thrombin burst |

Every peculiar dynamic of Layer 2 arises **from the structure and is not
scripted in**:

- **a pre-existing clone** (`BCELL0`) → why IgG can be detected within
  4–5 days even on first exposure
- **plasma cells dominate, long-lived plasma cells are almost zero**
  (`FLL` = 0.0008) → why antibody is **transient**
- **`KRECALL` defaults to 0** → why re-exposure produces **no anamnestic
  boost**

### The autocatalytic loop

Activated platelets pour α-granule PF4 into the very plasma that is
activating them (`PLT_GRAN → PF4P → ULC → activation`). **Antigen creates
activation, and activation creates antigen.** This is why HIT is
explosive rather than gradual, and **why platelet transfusion is
anti-therapeutic** — it feeds the loop's substrate.

---

## 3. The central paradox, and why risk was given two pathways

> **A disease defined by a low platelet count, and the patient dies of
> thrombosis.**

The very FcγRIIa signal that tags platelets for splenic clearance also
releases procoagulant microparticles and induces monocyte tissue factor.
So the model routes `ACTSIG` to **both sinks**, giving thrombosis risk
**two additive pathways.**

| Pathway | Driver | What breaks it |
|---|---|---|
| `HZTHR` | thrombin | argatroban · bivalirudin · DOACs |
| `HZPLT` | platelet activation | **no anticoagulant breaks this one.** Only stopping heparin, IVIG, or plasma exchange |

This one structural choice splits drug classes across **two independent
axes rather than one**, and is why residual risk remains even in patients
treated with a DTI (model: 4.8% vs 51.0% untreated).

---

## 4. Two pharmacological traps — generated, not scripted as warnings

### (i) The warfarin / venous limb gangrene trap

Protein C has a half-life of 8 hours; prothrombin, 60–72 hours. The model
contains nothing but an ODE with those two rate constants — no other
device. Starting a vitamin K antagonist while thrombin generation is
still at its peak causes the model to reproduce the literature's
discrepancy on its own.

| After a 10 mg loading dose | Protein C | Prothrombin | Measured INR | Effective thrombin |
|---|---|---|---|---|
| 0 h | 100% | 100% | 1.00 | 2.1 nM |
| 6 h | 69.8% | 95.1% | 1.19 | 9.9 nM |
| 12 h | 51.6% | 90.3% | 1.37 | **10.7 nM** |
| 24 h | **37.5%** | **82.5%** | 1.57 | 8.8 nM |

This matches observation (protein C about 30–40% at 24 hours, prothrombin
about 80%), and the net effect is a **transient rise in thrombin.** And a
sweep of the day warfarin is started (D6) **generates** the guideline
rule as a monotonic dose-response rather than as text:

| Day warfarin started | Platelets at start | Probability of limb gangrene (45 days) |
|---|---|---|
| day 7 | 108 | **1.57%** |
| day 9 | 128 | 0.13% |
| day 11 | 166 | 0.00% |
| day 15 or later | 242+ | 0.00% |

### (ii) The argatroban–INR trap

Argatroban prolongs PT without lowering vitamin-K-dependent factors at
all. So the model computes INR from **two inputs** (true factor activity +
the DTI assay artefact) and reports a chromogenic factor X alongside it,
making the dissociation visible.

| Day 21 | Measured INR | True INR | Chromogenic factor X |
|---|---|---|---|
| with argatroban 1.20 µg/mL | **4.31** | **2.07** | **47.1%** |

Stopping argatroban on the strength of a measured INR of 4.31 would leave
the patient unprotected.

---

## 5. The iceberg is an output, not an assumption

Antibodies after exposure are common; disease is rare. The model splits
**antibody presence → platelet-activating capacity → clinical HIT** into
three sequential thresholds, and reports each level through the assay
that measures it. Diagnostic analysis D4 shows that what sets a patient's
position on the iceberg is **exposure intensity.**

| UFH infusion rate | Concentration | Max ULC | ELISA OD | ① seropositive | ② SRA positive | Platelet drop | ③ clinical HIT |
|---|---|---|---|---|---|---|---|
| 2 U/kg/h | 0.036 | 0.8 | 0.06 | — | — | 0.1% | — |
| 4 U/kg/h | 0.074 | 3.5 | 0.10 | — | — | 0.3% | — |
| **8 U/kg/h** | 0.157 | 25.6 | 1.24 | **✓** | **✓** | **47.9%** | **—** |
| 12 U/kg/h | 0.249 | 57.5 | 2.29 | ✓ | ✓ | 77.5% | **✓** |
| 18 U/kg/h | 0.404 | 115.0 | 2.40 | ✓ | ✓ | 78.4% | ✓ |

The 8 U/kg/h row is the waterline of the iceberg — seropositive and
functionally positive, but the platelet drop stalls at 47.9% and never
becomes clinical HIT. The cardiopulmonary bypass scenario gives the same
result: seroconversion at OD 0.82, but platelets fall only 21.4% —
consistent with the observation that seroconversion after cardiac surgery
is common while clinical HIT is 1–3%.

---

## 6. Every variant is a single coefficient

Autoimmune HIT, spontaneous HIT syndrome, and VITT all share this entire
model, and **exactly one thing changes**: `HIND` (the heparin-independence
coefficient) moves from 0 to 1. Their antibodies bind the heparin-binding
site of PF4 itself, clustering PF4 even without heparin present. Every
clinically strange feature of these variants then follows **with no
further modification.**

| `HIND` | ULC at 14 days after discontinuation | Nadir | SRA (no heparin) | SRA (100 U/mL) | Thrombosis (30 days) |
|---|---|---|---|---|---|
| 0.00 (classical HIT) | 0.00 | 99.9 | 0.0% | 0.0% | 4.8% |
| 0.15 | 0.08 | 82.9 | 79.4% | 79.6% | 6.3% |
| 0.30 | 0.14 | 73.5 | 88.5% | 88.5% | 8.0% |
| 0.60 (aHIT) | 1.80 | 63.6 | 92.2% | 92.2% | 56.4% |
| 0.90 (VITT) | 2.97 | 59.0 | 92.8% | 92.8% | 65.6% |

- **platelets keep falling even after heparin is stopped** — because there
  is nothing left to stop
- **the standard SRA loses its high-dose heparin inhibition** — the model
  reports 93% at low heparin, 93% at 100 U/mL, and 93% with no heparin at
  all. This is precisely why such sera need a PF4-enhanced assay, and why
  a classical SRA can come back negative
- **IVIG does what anticoagulation cannot** — receptor competition does
  not care where the antibody bound (aHIT: IVIG 7.7% vs argatroban alone
  56.4%)

---

## 7. The baseline state is solved, not fitted

Every baseline equilibrium is derived algebraically from the values
declared in `$MAIN` (`SPF4`, `KECD`, `KBS`, `KMK`, `STPO`, `KTF0`, `KVWS`,
`KTMS`, `KATS`, `KPT`, `KFBN`, `KDDS`, `KF12S`). So over a **200-day run
with no heparin, every state drifts by 0.0000%.** In this model, disease
is **generated by exposure, never assumed.**

---

## 8. Calibration anchors

| Item | Model | Observed |
|---|---|---|
| onset of platelet drop (therapeutic UFH) | 6.4 days | 5–10 days (median 6) |
| platelet nadir | 54 ×10⁹/L | median 55–60 |
| nadir below 20 is rare | 50–100 | a feature that distinguishes it from ITP/DIC |
| 30-day thrombosis, untreated | 51.0% | 20–50%+ |
| 30-day thrombosis, DTI-treated | 4.8% | about 6–14% |
| plasma PF4 on therapeutic UFH | 0.5 µg/mL | 0.05–0.2 (up to about 2 post-CPB) |
| protein C at 24 h after warfarin 10 mg | 37% | about 30–40% |
| prothrombin at 24 h after warfarin 10 mg | 83% | about 80% |
| argatroban steady state, 2 µg/kg/min | 1.20 µg/mL | about 1.2 |
| steady-state aPTT | 70 s | 60–80 (target 1.5–3×) |
| SRA, 0.1–0.3 U/mL heparin | 93% | >20% = positive |
| SRA, 100 U/mL (classical HIT) | 14% | <20% = confirmatory |
| SRA, 100 U/mL (VITT) | 93% | not inhibited |
| immunoassay seroreversion | 120–180 days | about 85–100 days |
| fondaparinux ULC formation | about 0 | no ULC forms (structural null) |

---

## 9. Negative and self-refuting results — reported, not removed

The mark of a good model is that it can refute its own design
assumptions. All of the following are things this model refuted about
itself, and they are reported rather than deleted.

1. **"About 50% residual thrombosis from stopping heparin alone" is not a
   property of discontinuation itself.** The model reproduces this figure
   only for **late diagnosis.** Recognised at day 5, residual risk is
   0.9%; at day 9, 18%; at day 14, 32%; at day 30, 51%. The model's claim
   is that "discontinuation alone" was never a single intervention to
   begin with — it is a family of interventions indexed by the day
   diagnosis was made, and the historical figure is a statement about how
   late HIT used to be recognised.

   | Day of recognition | IgG at discontinuation | Nadir | Peak thrombin | Thrombosis (30 days) |
   |---|---|---|---|---|
   | day 5 | 3.8 | 222.3 | 2.50 nM | **0.9%** |
   | day 7 | 23.1 | 99.8 | 9.35 nM | 9.1% |
   | day 9 | 100.1 | 58.4 | 12.79 nM | 18.0% |
   | day 14 | 532.1 | 54.1 | 14.70 nM | 32.1% |
   | day 30 | 1495.5 | 54.1 | 15.32 nM | **51.0%** |

2. **The received idea that a heparin flush sustains HIT is not
   supported.** This was originally the hypothesis built into the model as
   an explanation for item 1 above, and it failed. Subtherapeutic heparin
   concentrations barely generate ULC while still accelerating
   antithrombin, so the model places residual risk as flat or even
   slightly **lower** (9.1% at 0 U/day → 7.6% at 2000 U/day). The failed
   hypothesis is reported rather than removed.

3. **The pathogenic-antibody fraction `FPATH` does not create the
   iceberg.** Once therapeutic-dose exposure occurs, antibody is so
   abundant that even a 2% platelet-activating fraction is enough. In
   this model the iceberg is created by **exposure intensity** (§5), which
   is the stronger, more testable claim.

4. **Prophylactic UFH comes out at slightly lower risk than prophylactic
   LMWH** — the opposite direction from epidemiology. The model assigns
   all of UFH's advantage to chain length and concentration, so at the low
   concentrations reached by subcutaneous prophylactic dosing that
   advantage is exhausted. This is a genuine limitation, not a fitted
   result.

5. **Plasma exchange is nearly powerless here** — removing IgG does not
   remove the plasma cells producing it. The model says plasma exchange is
   a **bridge**, not a therapy.

6. **Rapid-onset reactivity outlasts the immunoassay-positive window**
   (D8). Even a few U/mL of residual platelet-activating IgG can drop
   counts within 48 hours on re-exposure to a large new dose of antigen.
   So the model treats the 100-day re-exposure rule as an
   **underestimate** of the true risk window — a testable claim, not a
   fitted one.

---

## 10. Using the files

```bash
# render the mechanistic map
dot -Tsvg hit_qsp_model_en.dot -o hit_qsp_model_en.svg
dot -Tpng -Gdpi=150 hit_qsp_model_en.dot -o hit_qsp_model_en.png

# run all 31 scenarios + 14 diagnostic analyses (requires mrgsolve)
Rscript hit_mrgsolve_model.R

# interactive dashboard (10 tabs)
Rscript -e 'shiny::runApp("hit_shiny_app.R")'
```

### mrgsolve model structure (59 ODEs, time unit = hours)

| Module | Compartments |
|---|---|
| Heparin/polyanion PK | `HEPD` `HEP` `PROT` `LMWD` `LMWC` `FOND` `FONC` `DANC` `SURG` |
| PF4 and antigen | `PF4P` `PF4EC` `PF4DNA` `ULC` `IC` |
| Immune | `BCELL` `PBLAST` `LLPC` `MEMB` `IGGP` `IGGN` |
| Platelets | `PLT` `PLTA` `MK1` `MK2` `TPO` `MP` |
| Leukocytes · endothelium | `MOA` `TF` `NET` `ECA` `VWF` `TM` |
| Coagulation | `FII` `FVII` `FIX` `FX` `PC` `PROS` `AT` `THR` `FBN` `DDIM` `FRG` |
| Clinical outcomes | `TEC` `NECR` `THRAUC` |
| Treatment PK | `ARGC` `ARGP` `BIVC` `RIVD` `RIVC` `IVGC` `IVGP` `WARD` `WARC` `RTXC` `RTXP` `PLEXA` `VITK` |

### Key tunable parameters

| Parameter | Default | Meaning |
|---|---|---|
| `HIND` | 0 | **heparin independence** — 0 for classical HIT, 0.3–0.7 for aHIT, 0.8–1.0 for VITT |
| `RSTAR` | 8 | the optimal PF4 tetramer : chain molar ratio |
| `LEN50` | 24 | chain length (saccharide units) at half-maximal bridging |
| `FCG` | 1.0 | FcγRIIa 131 genotype (HH 1.35 / HR 1.0 / RR 0.75) |
| `KRECALL` | **0** | memory recall — 0 matches observation |
| `FPATH` | 0.15 | platelet-activating IgG fraction |
| `HEPFN` | 1.0 | hepatic function (adjusts argatroban dosing) |
| `RENFN` | 1.0 | renal function (adjusts fondaparinux/danaparoid/DOAC dosing) |

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational
and research purposes.** It was built from published literature and
clinical trial data but has not been independently verified or certified,
and **must not be used directly for clinical decision-making,
prescribing, or regulatory submission.** HIT in particular is an
emergency requiring immediate expert judgement, and no output of this
model can substitute for an individual patient's anticoagulation
decisions.
