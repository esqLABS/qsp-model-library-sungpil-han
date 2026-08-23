# Infective Endocarditis (IE) — QSP Model

> **The MIC is a number measured on 10⁵ CFU/mL of free-floating,
> log-phase bacteria in broth. This lesion is 10⁹–10¹¹ CFU/g of
> stationary-phase bacteria sitting behind an avascular
> platelet–fibrin matrix 2–5 mm thick.**
>
> Every peculiarity of endocarditis therapy — the 4–6 week duration, the
> bactericidal-only rule, combination therapy, surgical indications, the
> ceiling on vancomycin exposure — is a **consequence** of this mismatch.
>
> **Infective endocarditis is a geometry problem wearing the clothes of a
> susceptibility problem.**

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (252 nodes · 22 clusters) | [`ie_qsp_model.dot`](../../../infective-endocarditis/ie_qsp_model.dot) · [SVG](../../../infective-endocarditis/ie_qsp_model.svg) · [PNG](../../../infective-endocarditis/ie_qsp_model.png) |
| ⚙️ mrgsolve ODE model (55 compartments · 208 parameters · 16 scenarios · 16 diagnostic functions) | [`ie_mrgsolve_model.R`](../../../infective-endocarditis/ie_mrgsolve_model.R) |
| 📊 Shiny dashboard (9 tabs) | [`ie_shiny_app.R`](../../../infective-endocarditis/ie_shiny_app.R) |
| 📚 References (149 items) | [`ie_references.md`](ie_references.md) |

<a href="../../../infective-endocarditis/ie_qsp_model.svg"><img src="../../../infective-endocarditis/ie_qsp_model.png" width="820" alt="IE QSP mechanistic map"></a>

---

## 1. The premise this model sets up differently

Almost every antimicrobial PK/PD model is built on **a single MIC**. But
the MIC is obtained by observing, for 18–20 hours, 10⁵ CFU/mL of
log-phase, planktonic bacteria stirred in nutrient broth. The endocarditis
lesion satisfies none of those conditions.

So instead of asking "is this organism susceptible to this drug," the
model builds in the following five facts as structure. Each one carries
the condition that it must **be reproduced, not fitted** to the clinical
trial that demonstrated it.

### (1) The ledger — logs needed vs logs delivered

Cure means driving the **entire** vegetation below 1 CFU. Starting from
10⁹–10¹⁰ means achieving **9–10 log₁₀ of sustained net killing.**
Bacteraemia can be cleared in 14 days, but this is why endocarditis needs
4–6 weeks. `IE_ledger()` prints a shell-by-shell statement — logs
required, logs delivered, and **at what depth the shortfall lies.**

```
              regimen   req_log10  delivered  surface   mid   core  sterile
           no therapy        9.33      -0.46    -0.44 -0.39  -0.64       no
     nafcillin 12 g/d        9.33       9.40    11.14  9.72   8.59      YES
      cefazolin 6 g/d        9.33      12.33    12.09 11.80  11.45      YES
   vancomycin AUC 500        9.33       5.51     7.65  6.06   4.67       no
   daptomycin 8 mg/kg        9.33      12.33    12.09 11.80  11.45      YES
```

Achieving 9.33 logs on average is not cure. It must be achieved **in every
shell.** Vancomycin delivers 7.65 logs at the surface but only 4.67 at the
core, and that gap is the failure.

### (2) The exponential — the penetration length λ

Drug concentration inside an avascular lesion falls as
`C(x) ≈ C_surface · exp(-x/λ)`. Here, `λ = √(D/k_bind)` is a property of
**the molecule and the matrix**, not of the dose. Doubling plasma AUC
simply raises this entire exponential curve, without changing its slope.

Core exposure measured at steady state (`IE_penetration()`, geometry
fixed):

| Drug | λ (mm) | Core AUC / free plasma AUC | Homogenate-equivalent |
|---|---|---|---|
| Rifampicin | 2.50 | **0.4712** | 0.764 |
| β-lactam | 1.10 | 0.0843 | 0.293 |
| Daptomycin | 0.50 | 0.0040 | 0.101 |
| Vancomycin | 0.45 | **0.0025** | 0.095 |
| Gentamicin | 0.30 | 0.0001 | 0.045 |

The last two columns must be read **together.** The "vegetation:serum
ratio of 0.2–0.4" reported in the literature is measured on **homogenised**
vegetation — a **volume-weighted average** dominated by the bulky outer
shell — and, on top of that, it is **total drug**, including what is bound
to fibrin. But that binding is exactly what shortens λ. **The standard
assay reassures you for exactly the wrong reason.**

### (3) A ceiling, not a maximum — vancomycin has no safe window

Core exposure **saturates** with respect to plasma exposure (because the
dose does not appear inside the exponential decay). AKI does not saturate
— instead it enters a positive-feedback loop: `AKI → CrCL ↓ → clearance ↓
→ exposure ↑ → AKI`.

Result of `IE_aucsweep()` (MRSA, 7 mm vegetation, 42 days):

| Target AUC24 | Core log-kill | Peak AKI index | Benefit per +100 | Harm per +100 | Benefit/harm |
|---|---|---|---|---|---|
| 400 | 3.67 | 0.28 | 1.74 | 0.13 | **13.4** |
| 500 | 5.31 | 0.48 | 1.64 | 0.20 | 8.2 |
| 600 | 6.81 | 0.79 | 1.50 | 0.31 | 4.8 |
| 700 | 8.26 | 1.28 | 1.45 | 0.49 | 3.0 |
| 800 | 9.77 | 2.17 | 1.51 | 0.89 | 1.7 |
| 1000 | 11.40 | 5.32 | 0.82 | 1.58 | **0.5** |

**The AUC24 that first achieves sterilisation is 800, and the AUC24 at
which KDIGO stage-1 AKI first appears is also 800.** For a lesion of this
size, there is no dial setting that protects the kidney while killing the
lesion. This is the mechanistic account of why vancomycin keeps being
displaced in staphylococcal endocarditis, and it means the 400–600 window
in the 2020 consensus guideline should be read not as "the dose that
reliably sterilises" but as "the range where additional exposure is still
worth it (benefit/harm > 3)."

### (4) The growth-rate trap — this is not resistance, it is **tolerance**

β-lactam killing requires an actively remodelling cell wall, so killing
rate is proportional to μ. In the core, μ ≈ 0. Here, **the MIC the lab
reports does not change, and tells you nothing.** Daptomycin (membrane
depolarisation) and rifampicin (transcription) do not require μ, which is
the sole mechanistic reason these two drugs are special inside a
vegetation, and why the bacteriostatic linezolid cannot handle a
10-log lesion on its own.

The model's `FSTAT` vector (fraction of activity retained against
non-dividing cells) is this entire argument in one place: β-lactam 0.02 ·
vancomycin 0.10 · daptomycin 0.55 · rifampicin 0.60 · gentamicin 0.05.

### (5) Resistance is not bad luck, it is arithmetic

rpoB mutations arise at a frequency of about 10⁻⁸ per replication. Growing
a lesion from a seed to 10⁹–10¹⁰ means **rifampicin-resistant organisms
already exist before the first dose is given.**

```
wild-type burden ...................... 1.892e+09 CFU
rpoB-mutant burden (EMERGENT) ......... 1.253e+02 CFU     ← 125 already present at diagnosis
naive expectation f_mut x N ........... 1.892e+01 CFU
```

Nothing was seeded into the model at t=0. The mutation term is
`MUTRATE × μ × B` (per-replication occurrence), so the mutants are a
**by-product** of the lesion's own growth, and the 6.6-fold excess over
the naive product (19 organisms) is Luria–Delbrück accumulation over the
growth period. So rifampicin monotherapy fails **deterministically, not
stochastically.**

---

## 2. What the model had to reproduce (validation targets T1–T10)

| # | Target | Literature | Model |
|---|---|---|---|
| T1 | MSSA + anti-staphylococcal penicillin, blood culture clearance | 2–3 days | **3.6 days** |
| T1 | MRSA + vancomycin (AUC 500) | 7–9 days | **6.2 days** |
| T1 | MRSA (MIC 2) + vancomycin | delayed | **10.2 days** |
| T1 | MRSA + daptomycin 8 mg/kg | 4–8 days | 2.3 days ⚠️ |
| T2 | Vancomycin vegetation:serum ratio | 0.2–0.4 (homogenate) | homogenate-equivalent 0.095, **core 0.0025** |
| T3 | Front-loading of embolic risk | ~65% | **89%** (pre-diagnosis + first 2 weeks of therapy) |
| T3 | Risk reduction after effective therapy | ~10-fold | week 1 vs week 5–6 **6.6-fold**, vs no therapy **48-fold** |
| T4 | Vegetation size → embolic risk | >10 mm, >15 mm | reproduced structurally (diameter^1.6 × vulnerability) |
| T5 | Adjunctive gentamicin (Cosgrove 2009) | **no benefit, nephrotoxic** | +0.18 log, AKI 0.00 → **0.63** |
| T6 | Adjunctive rifampicin (ARREST 2018) | **no benefit** | **0.00 log** in a bacteraemia-sized lesion |
| T7 | Aspirin (Chan 2003) | **no reduction in emboli, more bleeding** | P(embolism) 0.0137 → **0.0136** |
| T8 | Partial oral switch (POET 2019) | non-inferior | IV 10 days → oral: **sterile**, IV 10 days then stopped: recurrence at 10^9.8 |
| T9 | Early surgery (EASE 2012) | embolic events drove the composite endpoint | surgery at day 2 gives a **57% reduction** in post-diagnosis embolic probability; microbiological outcome nearly unchanged |
| T10 | Vancomycin AKI (Rybak 2020) | sharp rise above AUC ~600 | KDIGO stage 1 first appears at **AUC 800** |

**T5, T6, and T7 are negative trials.** A model that cannot reproduce a
negative result is a curve fit, not a mechanistic model. All three are
reproduced by **structural properties that were built in for other
reasons**:

- **Gentamicin** — uptake requires a proton motive force, so `OXY =
  AVL^1.5` multiplies it, and λ is the shortest of any drug in the file.
  It does nothing in the core. Nephrotoxicity, on the other hand, is
  driven by cumulative **plasma** AUC, which has nothing to do with
  geometry. **The benefit is geometric and the harm is not.**
- **Rifampicin** — ARREST randomised patients with *S. aureus*
  **bacteraemia**, of whom only about a tenth had endocarditis.
  Rifampicin's specific advantage is reaching diffusion-limited
  compartments, and there is nothing to give in a patient without one. In
  the model's row 1 (diagnosis at day 4, a 1.9 mm lesion), the benefit is
  exactly 0.00 log.
- **Aspirin** — acts on matrix accumulation. But in this model, the
  matrix-forming drive tracks the **actively growing bacteria** at the
  surface, and effective antibiotics kill those within 48 hours. That is,
  the antiplatelet agent subtracts **30% of an already near-zero rate.**
  Meanwhile, embolic risk = f(diameter) × **vulnerability**, and
  vulnerability tracks the same surface population, so it has already
  collapsed for reasons unrelated to platelets. There is no room left for
  the drug to do work. All that remains is the bleeding.

---

## 3. Where the model **disagrees** with current practice (recorded, not hidden)

`IE_negatives()`'s T6 table lays out three situations side by side.

```
                             setting  veg_mm_dx  vanco_end  vanco_rif_end  benefit_log10
 bacteraemia, minimal lesion (day 4)        1.9      -3.00          -3.00           0.00
            native-valve IE (day 21)        7.1       3.14          -3.00           6.14
        prosthetic-valve IE (day 21)        8.0       6.77          -3.00           9.77
```

Row 1 reproduces ARREST. But **row 2 is where the model's prediction
departs from the guidelines** — the guidelines restrict rifampicin to
prosthetic valves, while this model predicts a large benefit in
native-valve endocarditis too. For the guideline to be right and the
model wrong, at least one of the following must hold:

- in real vegetations, `LAM_RIF` is **much** smaller than 2.5 mm,
- the fitness cost of the rpoB mutant is smaller than `FIT_COST = 0.15`,
  letting it survive alongside the partner drug,
- the limiting factor was never the pharmacology at all, but rather
  **treatment discontinuation from drug interactions and hepatotoxicity**,
  which this model does not have.

`IE_lambda()` shows how small a move in λ flips the answer (λ_van 0.45 →
failure, 0.70 → sterile). This puts on display the fact that this is
simultaneously the most consequential and the hardest-to-measure parameter
in the model.

---

## 4. Honest limitations (Known Limitations)

1. **Daptomycin clears bacteraemia too fast** (2.3 days vs 4–8 days in the
   literature). This is a clear miss. The cause is that the model has no
   **metastatic foci** (vertebral osteomyelitis, splenic abscess, etc.),
   so the only source that keeps replenishing the bloodstream is a single
   intracellular reservoir. Nafcillin (3.6 days) and vancomycin (6.2 days)
   land within the target range, so the problem is concentrated in this
   one drug's `FCELL_DAP` value, not in the structure.
2. **The 3 shells are a coarse discretisation of a continuous diffusion
   field.** They capture the exponential decay and phenotype gradient, but
   the near-surface concentration profile is not accurate at short time
   scales.
3. **Rifampicin's co-drug induction terms (`FRIND_VAN`, `FRIND_DAP`) default
   to 0.** Both drugs are renally cleared and the reported interactions are
   inconsistent. They are exposed as parameters and can be switched on
   directly. The oral combination term (`FRIND_ORL`) is set to 1.0 because
   that interaction is real and large.
4. **Death is not a state variable.** It is a downstream function of heart
   failure, embolism, and persistent bacteraemia; turning it into a state
   would only add fitted parameters without adding mechanism.
5. **This is a single-lesion model.** Multivalve disease and CIED lead
   infection are not represented, and metastatic foci are a loss term, not
   a compartment.
6. **The front-loading of embolic risk is overstated** (89% vs ~65%). The
   direction and rough magnitude are right, but real patients continue to
   grow new lesions and re-embolise after diagnosis, which the model
   underweights.
7. **The vegetation's size at 6 weeks (7.3 → 2.3 mm) shrinks faster than
   clinical observation.** Real vegetations remain visible on echo for
   months even after sterilisation. `KLYS` should be read as the loss rate
   of the dense fibrin that acts as the drug-diffusion barrier, not as the
   shrinkage rate of the mass seen on echo.
8. **The oral-combination slot is a composite of a two-drug regimen**
   (`EMAX_ORL 0.65`, `FSTAT_ORL 0.30`, `LAM_ORL 1.6`), collapsing the
   actual regimens POET individualised by susceptibility into a single
   representative drug.

---

## 5. Model structure

**55 ODE compartments**, time in **hours**, bacterial states in
**absolute CFU counts** (because cure is a statement about an absolute
count).

| Group | Compartments |
|---|---|
| Antibiotic PK (7 drugs) | vancomycin 2-compartment · daptomycin 2-compartment · generic β-lactam 2-compartment · gentamicin 2-compartment · rifampicin (oral + auto-induced enzyme pool) · dalbavancin 2-compartment · oral combination 2-compartment + 3 exposure integrals |
| Drug within the vegetation | 7 drug-specific diffusion-delay compartments, each read out at 3 depths via `exp(-d/λ)` |
| Bacteria | **3 shells** (surface/mid/core) × **2 phenotypes** (growing/dormant) × **2 genotypes** (wild-type/rpoB mutant) = 12 |
| Bloodstream · reservoirs | blood CFU · intracellular reservoir (endothelial/phagocyte) |
| Lesion · host | matrix mass · perivalvular abscess · IL-6 · CRP · PCT · neutrophils |
| Haemodynamics | valve integrity (ratchet) · LVEDV · NT-proBNP · PR interval |
| Embolism | cumulative embolic risk · cumulative CNS embolic risk |
| Toxicity | AKI index · CK · ALT · cumulative gentamicin AUC |
| PD metric | integrated core β-lactam fT>MIC |

### The patient at diagnosis is **not an input**

`IE_seed()` simulates 21 days of natural history from a 10³ CFU seed. The
resulting state is the patient at the moment of diagnosis:

```
  vegetation diameter ............ 7.3 mm
  total viable burden ............ 10^9.33 CFU
  dormant fraction ............... 0.59
  pre-existing rpoB mutants ...... 112 CFU
  bacteraemia .................... 1.55 CFU/mL
  CRP / NT-proBNP ................ 83 mg/L / 500 pg/mL
  regurgitant fraction ........... 0.35
```

**None of these values were set by hand.** Vegetation size, dormant
fraction, and — crucially — the number of pre-existing mutants are all
consequences of how the lesion grew, not free parameters. Changing
`dx_day` produces a different disease:

```
 dx_day  veg_mm  log10_dx  dormant_frac  core_xMIC  logs_delivered  sterile
      6     2.1      8.06         0.414      2.252           11.06      YES
     10     2.7      8.35         0.454      2.012           11.35      YES
     14     3.8      8.74         0.502      1.620           11.74      YES
     18     5.7      9.11         0.555      1.152           11.77      YES
     25     9.3      9.53         0.637      0.622            7.30       no
```

A 2 mm vegetation has no protected core at all. Same organism, same drug,
same MIC, **a different disease.**

---

## 6. How to run it

```r
install.packages("mrgsolve")          # developed and verified on 2.0.1
source("ie_mrgsolve_model.R")
mod <- ie_model()

IE_report(mod)                        # runs every diagnostic in the file, in order
```

Individual diagnostic functions:

| Function | Contents |
|---|---|
| `IE_scenarios()` | master table of all 16 scenarios |
| `IE_ledger()` | **the ledger** — logs required/delivered per shell |
| `IE_penetration()` | core:plasma exposure ratio (geometry fixed) |
| `IE_clearance()` | time to blood-culture clearance (T1) |
| `IE_cefazolin()` | cefazolin vs nafcillin, how large an inoculum effect would need to be |
| `IE_rifmono()` | pre-existing rpoB mutants and the deterministic failure of monotherapy |
| `IE_aucsweep()` | vancomycin AUC — saturating benefit vs accelerating harm (T10) |
| `IE_duration()` | duration and recurrence |
| `IE_embolic()` | time distribution of embolic risk (T3) |
| `IE_surgery()` | sweep of surgical timing (T9) |
| `IE_negatives()` | **the three negative trials** (T5, T6, T7) |
| `IE_pve()` | prosthetic valve — where rifampicin is actually needed |
| `IE_infusion()` | continuous vs intermittent infusion |
| `IE_lambda()` | sensitivity to penetration length — the most honest table here |
| `IE_size()` | diagnostic delay → size → a different disease |
| `IE_poet()` | partial oral switch (T8) |

Shiny dashboard:

```bash
R -e 'shiny::runApp("ie_shiny_app.R", launch.browser = TRUE)'
```

9 tabs: ① the ledger ② patient · lesion ③ pharmacokinetics
④ **the exponential** ⑤ bacteriology ⑥ embolism · valve
⑦ toxicity ⑧ regimen comparison ⑨ exposure sweep.

---

## 7. Phenotype presets

| Preset | Settings |
|---|---|
| `MSSA_native` | nafcillin-like PK, MIC_BL 0.5, VIR 1.0 |
| `MRSA_native` | MIC_BL 64, MIC_VAN 1.0, MIC_DAP 0.5 |
| `MRSA_MIC2` | MIC_VAN 2.0 — the required AUC doubles to 800 and becomes unreachable |
| `MRSA_PVE` | λ_van 0.35 / λ_dap 0.40 (foreign-surface biofilm), higher abscess rate |
| `VGS_native` | ceftriaxone-like PK, MIC 0.03, low toxicity, low embolic potential |
| `EFAECALIS` | `EMAX_BL 0.25` — the cell-wall agent is **bacteriostatic**, ampicillin PK |
| `CONS_PVE` | slow growth, thick biofilm, low matrix dissolution rate |

The β-lactam slot holds only one drug at a time, so its PK must move
**together with** the dose (`IE_BLPK`: nafcillin / cefazolin /
ceftriaxone / ampicillin).

---

## ⚠️ Disclaimer

This model is for **educational, research, and hypothesis-generating
purposes only.** It must not be used for actual clinical decision-making,
prescribing, antibiotic selection, timing of surgery, or regulatory
submission. Treatment of infective endocarditis should follow current
guidelines (2023 ESC, 2015 AHA, 2023 Duke-ISCVID) and the judgement of a
multidisciplinary **endocarditis team** involving infectious disease,
cardiology, cardiac surgery, radiology, and clinical pharmacy.
