# Prosthetic Joint Infection (PJI) — QSP Model

> **Implant-associated *Staphylococcus aureus* osteomyelitis**
> A 171-node 18-cluster mechanistic map · 53-ODE mrgsolve model · 21
> treatment scenarios · 13-tab Shiny dashboard · 68 references, every PMID
> checked

![PJI QSP map](pji_qsp_model_en.png)

| File | Contents |
|------|------|
| [`pji_qsp_model_en.dot`](pji_qsp_model_en.dot) / [`.svg`](pji_qsp_model_en.svg) / [`.png`](pji_qsp_model_en.png) | Mechanistic map (171 nodes · 18 subgraphs) |
| [`pji_mrgsolve_model_en.R`](pji_mrgsolve_model_en.R) | 53-ODE mrgsolve model + 21 scenarios + verification functions |
| [`pji_shiny_app_en.R`](pji_shiny_app_en.R) | 13-tab interactive dashboard |
| [`pji_references_en.md`](pji_references_en.md) | 68 references (every PMID verified) |

---

## 1. The Question This Model Answers

Prosthetic joint infection treatment guidelines read, on the surface, like
a list of mutually unrelated rules of thumb.

- Debridement is mandatory
- Rifampicin is special
- But rifampicin must never be given alone
- DAIR (implant retention) only works within 3 weeks of symptom onset

This model does **not code these four in as rules**; it writes them so that
they emerge from three pieces of arithmetic. However the parameters are
perturbed, the rules themselves do not disappear, and if the arithmetic's
values change, the rules change with them.

---

## 2. The Three Arithmetic Pillars

### Pillar ① — The Foreign Body Changes the Inoculum Needed to Establish Infection, Not the Organism's Virulence

Only two independent mechanisms were put in. `FBFACT` (phagocytic killing
efficiency ×0.12 at the metal interface) and `KATT` (once planktonic
organisms move onto the surface, frustrated phagocytosis drops killing
efficiency to `FRUST`=0.03). No threshold was put in anywhere.

| Inoculum (CFU) | Implant present — day-60 log10 burden | No implant |
|---|---|---|
| 1 | Extinction | Extinction |
| 3 | **9.92 (infection established)** | Extinction |
| 10² – 10⁷ | 9.92 | Extinction |
| 3 × 10⁷ | 9.92 | **9.08 (established)** |

The same organism, the same host, the same inoculum of 100 CFU is an
infection next to metal and nothing at all otherwise. The no-foreign-body
threshold of ~10⁷ sits at the same point Elek & Conen measured in human
skin (10⁶–10⁷, ~10² with a silk suture present).

> **Honest caveat**: this model is deterministic and has no stochastic
> extinction. So "3 CFU with an implant" is not the measured ID50 (~10²)
> but its **lower bound**. What the model reproduces is the ≥10⁴-fold size
> of the shift, not the absolute value.

---

### Pillar ② — Every Antibiotic Is Summarised by a Single Number: "Free Bone Concentration ÷ MBEC"

This model does **not use the planktonic MIC** when deciding what happens
inside the biofilm. Each drug carries separately (i) AUC_bone/AUC_plasma,
(ii) the free fraction, (iii) the planktonic MIC, and (iv) the biofilm
MBEC, and at standard adult doses the ratios work out as follows.

| Regimen | AUC24 (mg·h/L) | Free bone concentration (mg/L) | MBEC | **C_bone,free / MBEC** |
|------|---:|---:|---:|---:|
| Vancomycin 1 g q12h | 500 | 2.083 | 512 | **0.0041** |
| Cefazolin 2 g q8h (MSSA) | 1500 | 2.500 | 256 | **0.0098** |
| Daptomycin 8 mg/kg qd | 700 | 0.350 | 32 | **0.0109** |
| Linezolid 600 mg q12h | 171 | 2.218 | 128 | **0.0173** |
| Levofloxacin 750 mg qd | 82 | 1.203 | 64 | **0.0188** |
| **Rifampicin 450 mg q12h (week 1)** | 68 | 0.281 | 1.0 | **0.281** |
| Rifampicin (autoinduced steady state) | 34 | 0.141 | 1.0 | **0.141** |

This is not because rifampicin is a better drug than the others. Its
target (RNA polymerase) remains essential even in adherent bacteria that
are barely dividing, so unlike a cell-wall synthesis inhibitor, its MBEC
does not climb 2–3 log above its MIC. The result is a 7–68-fold gap, and
**it still does not exceed 1.** There is no systemic route that crosses the
MBEC — the sole exception is the local concentration of antibiotic-loaded
cement (above 1000 mg/L for the first few days).

---

### Pillar ③ — Surgery Is Not Washing; It Is Manipulating the Supply of Mutants

P(pre-existing rpoB mutant) = 1 − e^(−μN), μ ≈ 10⁻⁸/division.

| Burden N | Expected mutants | P(pre-existing) |
|---|---:|---:|
| 10⁵ | 0.001 | 0.0010 |
| 10⁶ | 0.01 | 0.0100 |
| 10⁷ | 0.1 | 0.095 |
| 10⁸ | 1 | 0.632 |
| **10¹⁰ (mature biofilm)** | **100** | **≈ 1.000** |

The model reproduces this arithmetic **without being told to**: a lesion
left untreated at 10^9.92 accumulates **63 CFU** of rpoB mutants, against
an analytical μN of 100. And a scenario giving levofloxacin+rifampicin for
12 weeks without debridement ends with the resistant clone at **10^9.97
CFU** — occupying the entire lesion. Debridement does not "help" the
antibiotic; it is the precondition that makes Pillar ② usable at all.

---

### Pillar ④ — Biofilm Maturation Is a Clock, and DAIR Races It

Rather than holding the biofilm tolerance multiplier (MBEC/MIC) constant,
it is scaled by matrix maturity `EPS` (`PHIMAT`). EPS rises with a time
constant of about 10 days in response to surface bacterial burden.

| Post-infection | EPS | Maturity | Effective vancomycin EC50 |
|---|---:|---:|---:|
| Day 3 | 0.00 | 0.000 | 1 mg/L (same as planktonic) |
| Day 7 | 0.02 | 0.046 | 24 mg/L |
| Day 14 | 0.41 | 0.537 | 275 mg/L |
| Day 21 | 0.65 | 0.651 | 334 mg/L |
| Day 60 | 0.94 | 0.728 | 373 mg/L |
| Day 90 | 0.94 | 0.730 | 374 mg/L |

The biofilm in the first week of infection has almost the same
susceptibility as planktonic organisms, and by 3 weeks the effective EC50
has risen **more than 300-fold**, moving completely outside the 2.08 mg/L
that vancomycin can reach.

Changing only the timing of the same regimen (DAIR 2.5-log +
levofloxacin/rifampicin for 12 weeks):

| DAIR timing | Day sterility reached | Margin vs the 84-day course | Outcome |
|---|---|---|---|
| Day 21 of infection | Treatment day 60 | **24-day margin** | Cure |
| Day 90 of infection | Treatment day 84 | **0 days (the last day)** | Cutting it close |
| Day 180 of infection | Never reached | — | Relapse after antibiotics stop |

---

## 3. Reproducing Clinical Trials (Trial-level calibration checks)

Of the 21 scenarios simulated over 365 days, the ones that can be checked
against actual trials:

| Reference trial | Scenario | Model result | Actual |
|---|---|---|---|
| **Zimmerli 1998 JAMA** | 05 DAIR+LVX/RIF vs 06 DAIR+LVX alone | 05 cure / 06 failure (final 10^9.92) | RIF combination 100% vs alone 58% |
| **Rifampicin monotherapy contraindicated** | 07 DAIR+RIF alone | Cures deterministically, but **P(rpoB)=0.54 at the lowest burden** | Clinically contraindicated |
| **Debridement mandatory** | 08 LVX/RIF without debridement | Failure, resistant clone 10^9.97 | The principal risk factor for rifampicin resistance |
| **OVIVA 2019** | 11 oral LVX/RIF vs 12 IV VAN/RIF | Both arms cure | Oral non-inferior (13.2% vs 14.6% failure) |
| **DATIPO 2021** | 13 6 weeks vs 14 12 weeks (DAIR) | **6 weeks fails / 12 weeks cures** | 6 weeks is **not** non-inferior to 12 weeks |
| **Two-stage exchange superiority** | 09 one-stage vs 10 two-stage | Both cure, 10 has the lowest burden with P(rpoB)=0.000 | Two-stage has the lowest reinfection rate |

The rifampicin monotherapy result is reported as is, without softening:
after 2.5-log debridement the burden is 10^7.4, so the expected number of
rpoB mutants is 0.25, and the deterministic trajectory (= the majority of
patients) cures. **What the model gives as the basis for the
contraindication is not the trajectory but the 54% probability.**

---

## 4. Model Structure (53 ODEs)

| Block | State variables | Count |
|------|----------|------|
| PK — vancomycin · rifampicin · levofloxacin · daptomycin · linezolid · cefazolin | Central/peripheral/GI/**free bone concentration** compartments + CYP induction states | 17 |
| Local delivery | Spacer fast/slow pools + joint fluid concentration | 3 |
| Bacteria | Planktonic · biofilm · persisters · SCV · intracellular · rpoB clones (2) · gyrA clones (2) | 9 |
| Biofilm | EPS matrix · agr AIP | 2 |
| Immunity | Neutrophils · monocytes · MDSCs · macrophages · IL-1β · TNF-α · IL-6 · IL-10 | 8 |
| Biomarkers | CRP · ESR · α-defensin | 3 |
| Bone | RANKL · OPG · osteoblasts · osteoclasts · bone mass · loosening index · sequestrum | 7 |
| Toxicity | Creatinine · platelets · ALT · vancomycin AUC24 filter | 4 |

### Key Structural Decisions

- **The bone compartment is a free-concentration effect compartment**:
  `FB = (AUC_bone/AUC_plasma) × fu`. Mass is not returned to plasma
  (peri-implant mass is negligible).
- **Killing is `EMAX·C/(EC50_state + C)`**, with
  `EC50_state = MIC × biofilm PHI × state PHI`. There are five states —
  planktonic/biofilm/persister/SCV/intracellular — and the intracellular
  one multiplies the concentration by the intracellular accumulation ratio
  `CAR` (rifampicin 6, levofloxacin 5, vancomycin 0.2).
- **Resistant clones are generated from the division flux** — the number
  of divisions, not the burden, determines the mutant supply. rpoB clones
  neutralise only rifampicin, gyrA clones only fluoroquinolones.
- **Surgery is implemented without event plumbing** — a pulse of height
  ln(10)·LOGK/TSDUR applied for TSDUR time, which multiplies every
  bacterial compartment by exactly 10^−LOGK.
- **A retained implant keeps its mature matrix intact** (`EPSCUT`
  retention factor 0.05 vs 1.0 for exchange). This is where the reason
  DAIR races the clock comes from.
- **Rifampicin autoinduction** doubles clearance (ratio 0.28 → 0.14), and
  the same induction raises linezolid clearance by 50%, dropping its AUC
  by ~33% (Gandelman 2011's −32%).

---

## 5. Verification — 6 Defects Caught by an Independent Python/scipy Reimplementation

Before finalising the R model, all 53 ODEs were **independently
reimplemented** in Python/scipy and run. In the process, six defects
surfaced, each one quietly undermining one of the four pillars.

| # | Defect | Symptom | Fix |
|---|------|------|------|
| 1 | Immune capacity written as a Michaelis-Menten sink **per compartment** | The rare rpoB clone alone caught the entire phagocytic capacity at a first-order rate of ~19/h and went extinct → nullified Pillar ③ | A single shared capacity apportioned in proportion to accessibility (`NPHAG`/`SHR`) |
| 2 | The extinction threshold (<1 CFU) applied **per phenotypic state** | The clone was wiped out whenever it fell below 1 CFU in any one state. An untreated 10¹⁰ lesion had 0.01 mutants (vs an analytical value of 100) | Extinction changed to a property of the **clone** → 63, matching the analytical value |
| 3 | Debridement reset biofilm maturity | 2.5-log DAIR dropped EPS from 0.62 to 0.026, restoring the residual biofilm to planktonic-level susceptibility. **Levofloxacin monotherapy cured** | The adherent matrix on a retained implant cannot be removed (factor 0.55 → 0.05) |
| 4 | No carrying capacity on the non-adherent bacterial pool | Dispersed planktonic organisms grew to 1.8×10¹⁰, overwhelming the biofilm → the infection became curable by vancomycin | Added `NPCAP` (joint fluid · abscess) and `NICCAP` (host cell number limit) |
| 5 | Unbounded progression of the loosening index and osteoblasts | LOOSEN reached 266 on a 0–100 scale, osteoblasts fell to 8% of baseline | Added a saturation term, `KAPOP` 0.030 → 0.0040 |
| 6 | Phagocytic saturation constant mismatched the target ID50 | `KMIMM`=10⁴ against a target ID50 of 10⁷ → operated at a large first-order rate across a 3-log range, with immunity pinning the biofilm at 10^4.8 so that **PJI never established** | `KMIMM` 10⁶, `KPMN` 3000 → 500 |

### Remaining Calibration Gap (stated, not fixed)

The model's DAIR-timing threshold falls between **3–6 months** of
infection, while the clinical rule is **3 weeks** of symptoms. The
direction (later is worse) and the mechanism (EPS maturation ~10 days,
sequestrum shielding ~3 months) are correct, but the absolute clock is
slow. The reason is known — this model has no soft-tissue abscess or
mechanical loosening pathway, and in reality late DAIR fails for surgical
rather than microbiological reasons. Rather than tuning parameters to hit
3 weeks, the gap is recorded here.

---

## 6. Natural History (untreated, implant retained, 100 CFU inoculum)

| Timepoint | log10 burden | EPS | CRP (mg/L) | ESR | Joint fluid WBC (/μL) | PMN % | α-defensin | Bone mass % | Loosening |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Day 3 | 4.54 | 0.00 | 3.7 | 10 | 609 | 75 | 0.5 | 100.0 | 0.0 |
| Day 7 | 7.00 | 0.02 | 60.6 | 18 | 28,014 | 98 | 3.4 | 99.9 | 0.0 |
| Day 14 | 9.92 | 0.41 | 97.5 | 46 | 15,151 | 88 | 3.1 | 97.9 | 0.4 |
| Day 21 | 9.91 | 0.65 | 94.7 | 59 | 14,913 | 86 | 3.0 | 94.2 | 2.2 |
| Day 60 | 9.92 | 0.94 | 93.6 | 70 | 14,766 | 85 | 3.0 | 74.5 | 32.0 |
| Day 90 | 9.92 | 0.94 | 93.6 | 70 | 14,762 | 85 | 3.0 | 64.2 | 53.6 |
| Day 180 | 9.92 | 0.95 | 93.6 | 70 | 14,762 | 85 | 3.0 | 48.5 | 74.2 |

It transitions on its own from an acute-phase joint fluid WBC of
28,000/μL and PMN of 98% (acute thresholds >10,000 · >90%) to a
chronic-phase 14,800/μL and 85% (chronic thresholds >3,000 · >80%). This
decay arises from MDSC recruitment (`KSMDSC`, Heim 2014) and
leukotoxin-mediated neutrophil death, and was not coded in separately.

---

## 7. How to Run

```r
# model + 21 scenarios
source("pji_mrgsolve_model_en.R")

out <- run_all()            # 21 scenarios x 365 days
summarise_scen(out)         # summary of burden · resistance · cure probability · biomarkers · toxicity

ratio_table()               # Pillar ② — the table that explains rifampicin
mutant_supply()             # Pillar ③ — mutant supply before and after debridement
id50_scan()                 # Pillar ① — the ID50 shift caused by the foreign body

# dashboard
shiny::runApp("pji_shiny_app_en.R")
```

The Shiny app's 13 tabs: overview · PK (plasma vs bone) · **bone
concentration/MBEC ratio** · bacterial subpopulations · biofilm maturation
· local elution · **resistance/mutant supply** · host immunity ·
osteolysis · biomarkers · scenario comparison · toxicity · parameters.

---

## 8. Mechanistic Map (18 clusters · 171 nodes)

① Host and surgical risk factors and inoculum · ② implant surface ·
conditioning film · ③ adhesin MSCRAMMs · ④ biofilm matrix and structure ·
⑤ agr quorum sensing · dispersal · toxins · ⑥ bacterial subpopulation
state machine · ⑦ intracellular and osteocyte-lacunar reservoirs · ⑧
innate immunity and frustrated phagocytosis · ⑨ cytokines · systemic
inflammation · ⑩ bone remodelling · osteolysis · implant loosening · ⑪
antibiotic PK and bone penetration · ⑫ local delivery (ALBC) · ⑬ biofilm
pharmacodynamics (the MBEC barrier) · ⑭ resistance genetics and mutation
supply · ⑮ surgical strategy · ⑯ diagnosis · biomarkers (2018 ICM) · ⑰
drug toxicity · tolerability · ⑱ clinical endpoints

[View the full SVG](pji_qsp_model_en.svg)

---

## ⚠️ Disclaimer

This is a qualitative/semi-quantitative QSP model for educational and
research purposes. It was built from public literature and clinical trial
data, and all 53 ODEs were checked against an independent reimplementation,
but it has never been fitted to or validated against individual patient
data. **Do not use it for actual clinical decision-making, prescribing, or
regulatory submission.**
