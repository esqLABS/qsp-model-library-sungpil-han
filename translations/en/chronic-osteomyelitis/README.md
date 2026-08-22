# Chronic Osteomyelitis · Implant-Associated Bone Infection — QSP Model

*Quantitative Systems Pharmacology model · centred on Staphylococcus aureus*

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`com_qsp_model.dot`](com_qsp_model.dot) · [SVG](com_qsp_model.svg) · [PNG](com_qsp_model.png) | **158 nodes · 260 edges · 15 subgraphs** |
| mrgsolve model | [`com_mrgsolve_model.R`](com_mrgsolve_model.R) | **47 ODEs · 2 systemic drugs + 1 local depot · 5 bacterial subpopulations · 14 scenarios** |
| Shiny app | [`com_shiny_app.R`](../../../chronic-osteomyelitis/com_shiny_app.R) | **10 tabs** |
| References | [`com_references.md`](com_references.md) | **85 papers (every one confirmed by PubMed lookup)** |

Every number in this document is a value obtained by **actually running the R model
in this repository under mrgsolve 2.0.1**. None of them was estimated by hand.

---

## 1. Why this model exists

In chronic osteomyelitis the drug and the bacteria **do not live in the same
compartment.** And treatment failure arises from two mutually independent penalties
being **multiplied together**.

### ① The spatial penalty — dose can push through it

The sequestrum and the implant biofilm have no perfusion. The drug gets in by
diffusion alone, and while it is getting in it is lost to binding, degradation, and
bacterial uptake. This is a reaction-diffusion problem, and the steady-state
penetration efficiency closes in the Thiele modulus.

```
phi   = L · sqrt(k_loss / D_eff)
eta   = tanh(phi) / phi
C_seq = eta · C_bone
```

The values the model actually computes:

| Diffusion distance L (cm) | 0.05 | 0.10 | 0.18 | 0.25 | 0.32 | 0.50 |
|---|---|---|---|---|---|---|
| Vancomycin eta | 0.722 | 0.437 | 0.248 | 0.179 | 0.140 | 0.089 |
| Levofloxacin eta | 0.942 | 0.808 | 0.587 | 0.450 | 0.358 | 0.231 |
| Daptomycin eta | 0.581 | 0.315 | 0.176 | 0.126 | 0.099 | 0.063 |

This is a **concentration** problem. Raise the concentration at the margin 100-fold
and the deep concentration rises 100-fold with it.

And in the model `L` is not a parameter but **a value the lesion sets**
(`L = (LBF0+KLEPS)·EPS_maturity + KLSEQ·SEQ^(1/3)`):

| Sequestrum SEQ (cm³) | 0.5 | 2 | 5 | 10 | 15 | 25 | Immediately after debridement (0.1, EPS reset) |
|---|---|---|---|---|---|---|---|
| Diffusion distance L (cm) | 0.162 | **0.218** | 0.272 | 0.325 | **0.363** | 0.418 | **0.056** |
| Vancomycin eta | 0.276 | **0.205** | 0.164 | 0.138 | **0.123** | 0.107 | **0.680** |
| Levofloxacin eta | 0.632 | 0.506 | 0.417 | 0.353 | 0.317 | 0.276 | 0.917 |
| Rifampicin eta | 0.428 | 0.323 | 0.260 | 0.217 | 0.195 | 0.169 | 0.784 |

### ② The phenotypic penalty — no dose can push through it

Biofilm bacteria are not 'slow bacteria'. The kill rate of beta-lactams,
glycopeptides, and aminoglycosides is **proportional to growth rate**, so as mu → 0,
kill → 0, and this is independent of concentration.

```
KILL_j = C_j/(C_j + EC50_j) · [ Emax_gd,j · (mu/mu_max) + Emax_gi,j ]
```

Kill rate per pool (1/d) at saturating concentration (h = 0.8) — model-computed
values:

| Drug | Planktonic (mu/mumax=0.35) | Biofilm margin (0.12) | **Persister (0)** | Emax_gi / Emax_gd |
|---|---|---|---|---|
| Vancomycin | 1.144 | 0.408 | **0.024** | 0.007 |
| Nafcillin | 1.696 | 0.592 | **0.016** | 0.003 |
| Levofloxacin | 2.600 | 0.944 | **0.080** | 0.011 |
| Linezolid | 0.624 | 0.256 | **0.064** | 0.040 |
| Daptomycin | 3.000 | 1.160 | **0.200** | 0.025 |
| Clindamycin | 0.764 | 0.304 | **0.064** | 0.032 |
| **Rifampicin** | 2.880 | 1.776 | **1.200** | **0.250** |

The values in the persister column are set by **Emax_gi alone**, independently of
concentration. Rifampicin is **50 times** vancomycin. This one column is the sole
reason rifampicin is special in this disease, and at the same time it becomes the
testable claim that

> **the observed MBEC/MIC = 100~1000 is not a potency shift (a displacement of
> EC50) but is mostly the shadow of Emax collapse.**

The model grants the EC50 displacement due to EPS binding **only 3-fold**
(`FBFEC = 3.0`) and reproduces all the rest through mu gating.

### The three interventions touch mathematically different objects

```
RSTER = [ drug kill (biofilm) + immune kill (biofilm) ] / [ mu_bf + persister supply ]

RSTER ≤ 1  →  a stable non-zero bacterial equilibrium. Extending the duration never brings sterilisation.
RSTER > 1  →  sterilisation is possible. All that remains is time, and that time depends on P0 only logarithmically.

T_ster ≈ ln(P0) / (KILLPS + KIMMPS + KWAKE)
```

| Intervention | Mathematical object | Mode of action |
|---|---|---|
| **Surgical debridement / implant removal** | initial condition P0 · surface area A_imp | logarithmic scale |
| **Rifampicin combination** | kill ceiling RSTER | linear (does it cross 1) |
| **Duration of dosing** | time to arrival T | and nothing beyond that |

Spread `T_ster = ln(P0)/(KILLPS + KWAKE)` out into a table and it is immediately
visible that the two cannot substitute for each other (model-computed values,
KWAKE = 0.08/d):

| KILLPS (1/d) | P0 = 1e8 | 1e5 | 1e2 | 1e1 |
|---|---|---|---|---|
| 0.00 (a drug powerless against persisters) | 230 d | 144 d | 58 d | 29 d |
| 0.01 (vancomycin) | 205 d | 128 d | 51 d | 26 d |
| 0.20 (FQ alone) | 66 d | 41 d | 16 d | 8 d |
| **1.11 (rifampicin combination)** | **15 d** | **10 d** | **4 d** | **2 d** |

Read it across and it is the effect of debridement; read it down and it is the effect
of rifampicin. One acts logarithmically and the other linearly, so **doubling one
cannot stand in for the other.**

### And the pathophysiology builds the pharmacokinetic barrier itself

```
infection → RANKL/OPG → bone resorption + intramedullary pressure → thrombosis of Haversian vessels
          → PERF falls → (a) bone-interstitial drug exposure falls
                         (b) neutrophil access falls        ← one cause, two losses
          → bone necrosis → SEQ rises → L rises → phi rises → eta falls
          → C_seq falls → kill falls → infection persists → (back to the start)
```

This loop is what makes the difference in cure rate between early treatment and the
chronic phase. Read the table above again: as the sequestrum grows from 2 → 15 cm³
the diffusion distance L goes from 0.218 → 0.363 cm and vancomycin's eta **falls 40%,
from 0.205 → 0.123.** Conversely, clear the sequestrum out by debridement and reset
the EPS and the same drug's eta **rises 3.3-fold, from 0.205 → 0.680** — part of the
effect of debridement is not that it reduced the bacteria but that it **reduced the
diffusion distance**. Which means the reason the same drug works less well in a
patient who presents late is not that "the organism got stronger" but that **the
geometry changed**.

Because perfusion gates drug delivery and neutrophil access **simultaneously**
(`kimm ∝ PERF`), a loss of perfusion is an event that loses two things at once. That
is why in this model **revascularisation and free flaps sit in the same place as
'antibiotic therapy'** — S8 below shows this numerically.

---

## 2. Model structure

### 2.1 Mechanistic map — 15 clusters

`com_qsp_model.dot` (158 nodes · 260 edges)

1. Routes of infection · inoculum (haematogenous / contiguous spread / direct /
   diabetic foot; a foreign body lowers the critical inoculum)
2. **Bacterial subpopulations · state transitions** — planktonic / biofilm margin /
   persister / intracellular·SCV / rpoB-resistant
3. Biofilm architecture · MBEC (EPS, PIA/PNAG, eDNA, agr QS, sigB, toxin-antitoxin)
4. Host immunity (neutrophils, NETs, macrophages, IL-1β/IL-6/TNF, complement,
   staphylococcal immune evasion)
5. **Bone destruction · sequestrum formation** (osteoblast death → RANKL/OPG →
   osteoclasts → thrombosis → dead bone)
6. **Perfusion · diffusion barrier = the spatial penalty** (Thiele modulus,
   bone:serum ratio, PAD, revascularisation)
7. Systemic PK · DDI (oral/IV, renal and hepatic elimination, rifampicin CYP3A4/P-gp
   induction)
8. **Target-site PK — the three sanctuaries** (sequestrum / intracellular / abscess)
   + local depot
9. **PD · kill kinetics = the phenotypic penalty** (Emax_gd vs Emax_gi, MSW, RSTER,
   T_ster)
10. Antibiotic regimens (7 backbones + rifampicin + suppressive therapy + duration ·
    route)
11. Surgical intervention (radical debridement · DAIR · one-/two-stage exchange ·
    dead-space management · flap · amputation)
12. Toxicity (AKI · thrombocytopenia · hepatotoxicity · tendinopathy · C. difficile ·
    efficacy-toxicity utility)
13. Biomarkers · clinical endpoints (CRP vs ESR, culture and sonication, MRI, cure ·
    relapse · amputation)
14. Patient covariates · Cierny-Mader staging
15. **Summary of the core equations** (①~⑧)

Because the map carries genuine feedback loops between clusters, it cannot be laid
out by dot's local cluster ranking (`newrank = true` is required). That fact is
itself a structural feature of this disease.

### 2.2 mrgsolve model — 47 ODEs

| Compartment group | State variables | Count |
|---|---|---|
| Backbone (A) PK | `GUTA CENA PERA CBA CSA CIA` | 6 |
| Partner (B) PK | `GUTB CENB PERB CBB CSB CIB` | 6 |
| Local depot | `ALOC CLOC` | 2 |
| Enzyme induction | `EIND` | 1 |
| Bacteria · biofilm | `BPL BBFA BBFP BIC BRES EPSM` | 6 |
| Immunity | `NEU MACR IL1 IL6 TNFA` | 5 |
| Acute-phase response | `CRP ESR` | 2 |
| Bone · structure | `OB OC RKL OPG BM SEQ PERF PUS AIMP` | 9 |
| Toxicity · exposure | `EGFRC PLT ALT XA XB` | 5 |
| Clinical | `PAIN TXD HAZ AUCA AUCB` | 5 |
| | | **47** |

Plasma → bone interstitium → sequestrum → intracellular is chained together with a
**multiplicative penalty** at each step, and at the end of that product the kill
equation is gated again by growth rate. Each of the five bacterial subpopulations
**sees the concentration in a different compartment**:

| Subpopulation | Concentration it sees | mu/mu_max | Character |
|---|---|---|---|
| `BPL` planktonic | bone interstitium `CBA/CBB` | 0.35 | maximal drug susceptibility |
| `BBFA` biofilm margin | **sequestrum `CSA/CSB`** | 0.12 (when mature) | eta penalty + mu penalty |
| `BBFP` persister | **sequestrum `CSA/CSB`** | **0** | zero susceptibility without Emax_gi |
| `BIC` intracellular·SCV | **intracellular `CIA/CIB`** | 0.005 | only cell-penetrating drugs arrive |
| `BRES` rpoB-resistant | bone interstitium (to proliferate it has to be at the margin) | 0.32 | the partner drug suppresses it |

Three decisions set the character of this model.

- **The initial value of `BRES` is not 0 but 10 CFU.** At mutation-selection
  equilibrium the frequency is `FMUT/FITCOST = 1e-8/0.1 = 1e-7`, and with 1e8
  organisms the expectation is 10 organisms. **The reason rifampicin monotherapy
  fails is not that resistance 'arises' but that it is already 'there'.**
- **A foreign body creates a local granulocyte defect** (Zimmerli 1984,
  PMID 6323536). In the model the implant surface area `AIMP` directly cuts immune
  access (`AIMP50 = 5 cm²`), and this one term is what produces "the same
  debridement fails if you leave the implant in and succeeds if you take it out".
- **Immune access is a function of biofilm maturity.** A freshly debrided,
  EPS-poor film is almost accessible (`FIMBFX = 0.85`), and once mature, access is
  blocked (`FIMBF = 0.06`). Debridement does not only reduce the bacterial count, it
  **resets maturity** — this is the actual route by which debridement cures.

---

## 3. Verification — what came out of actually running the model

`Rscript -e 'options(com.autorun = TRUE); source("com_mrgsolve_model.R")'`

`RSTER` is read at the moment that determines the outcome (immediately after surgery,
at drug steady state, day 7). Read at end of therapy it is already past
sterilisation, so both numerator and denominator become meaningless.

### 3.1 The 14 standard scenarios

| # | Scenario | RSTER(7d) | eta(A) | Persister kill (drug+immune, /d) | End-of-therapy reservoir (CFU) | Relapse probability | Outcome |
|---|---|---|---|---|---|---|---|
| S1 | Natural history (untreated) | 0.28 | 0.21 | 0.000 + 0.001 | 1.1e8 | 1.00 | chronic progression |
| S2 | Vancomycin 6 weeks, **no debridement** | 0.39 | 0.19 | 0.001 + 0.000 | 5.4e9 | 1.00 | **failure** |
| S3 | Radical debridement + vancomycin 6 weeks | 1.14 | 0.59 | 0.005 + **0.765** | 9.1e-4 | 0.00 | **cure** |
| S4 | DAIR + levofloxacin/**rifampicin** 12 weeks | 1.07 | 0.74 | **0.615** + 0.075 | 2.3e-31 | 0.00 | **cure** |
| S4b | DAIR + levofloxacin **alone** 12 weeks | 0.47 | 0.72 | 0.012 + 0.065 | 2.9e9 | 1.00 | **failure** |
| S5 | Two-stage exchange + oral linezolid/rifampicin 6 weeks | 1.84 | 0.86 | 0.737 + 0.713 | 2.2e-31 | 0.00 | **cure** |
| S6 | **Rifampicin monotherapy**, no debridement | 0.97 | 0.19 | 0.316 + 0.001 | 5.9e9 | 1.00 | **resistant replacement** |
| S6b | Rifampicin monotherapy + limited debridement | 1.05 | 0.29 | 0.499 + 0.014 | 4.8e-7 | 0.00 | cure (stochastic region) |
| S7 | **Local beads alone**, no debridement | 0.95 | 0.19 | 0.026 + 0.001 | 6.3e9 | 1.00 | **failure** |
| S8 | Diabetic foot bone infection, PAD uncorrected | 0.77 | 0.73 | 0.501 + 0.031 | 8.0e8 | 1.00 | **failure** |
| S8b | Same regimen + **revascularisation** | 0.91 | 0.76 | 0.582 + 0.042 | 5.6e-8 | 0.00 | **cure** |
| S10 | Long-term oral suppressive therapy (removal impossible) | 0.40 | 0.47 | 0.009 + 0.019 | 6.7e9 | 1.00 | suppression (not cure) |
| S11 | Debridement on **day 2** | 1.20 | 0.83 | 0.714 + 0.096 | 1.2e-37 | 0.00 | **cure** |
| S11b | Same regimen, debridement on **day 120** | 0.36 | 0.47 | 0.000 + 0.000 | 1.8e7 | 1.00 | **failure** |

Against the clinical literature:

- **S4 vs S4b — reproduces Zimmerli 1998 JAMA (PMID 9605897).** The only difference
  is rifampicin. RSTER 1.07 vs 0.47, reservoir 2.3e-31 vs 2.9e9. The RCT result in
  which treatment success in the rifampicin-combination arm greatly exceeded the
  control arm comes, in this model, out of the single parameter `BEMGI`.
- **S3 / S5 — OVIVA (PMID 30699315).** Route changes only `F` and `CL` in the model,
  so the fact that oral and intravenous give the same outcome for a drug with
  adequate bioavailability is **derived structurally.** S5 is entirely oral and still
  cures.
- **S11 vs S11b — same drug, same duration, different timing.** eta falls
  0.83 → 0.47 and RSTER collapses 1.20 → 0.36. **The drug did not get worse; the
  geometry did.**
- **S8 vs S8b — perfusion is an antibiotic.** Same drug, same dose, same duration.
  Add revascularisation alone and the reservoir moves 16 log, 8.0e8 → 5.6e-8. The
  reason the IWGDF guidelines put perfusion assessment first comes out of this one
  line.
- **S2 is the clinically most dangerous scenario.** Give 6 weeks of vancomycin
  without debridement and at end of therapy CRP has risen to 130 mg/L and eGFR has
  fallen 90 → 65.6, while the bacteria remain at 5.4e9. **This is the case where you
  pay the toxicity and buy nothing.**

### 3.2 The arithmetic of rifampicin monotherapy — not a probability but a reservation

```
P(rpoB mutant present) = 1 − exp(−FMUT · B_total),  FMUT = 1e-8/division
```

| B_total | 1e6 | 1e7 | 1e8 | 1e9 | 1e10 |
|---|---|---|---|---|---|
| P(mutant present) | 0.010 | 0.095 | **0.632** | **1.000** | **1.000** |

The typical bacterial burden in chronic osteomyelitis is 1e8~1e9. Rifampicin
monotherapy is therefore not a 'risky choice' but a **reservation for failure**. In
the simulation S6 (rifampicin monotherapy without debridement) reduces the wild type
by 2 log while pushing the rpoB subpopulation up from **log 1.0 → log 8.65** —
replacement complete.

### 3.3 The duration sweep — the conclusion of this entire model

**Implant retained (DAIR, 3 log debridement):**

| Duration | Regimen | RSTER(7d) | Persister kill (drug) | End-of-therapy reservoir | Relapse probability |
|---|---|---|---|---|---|
| 2 weeks | Levo + rifampicin | 1.12 | 0.644 | 331 | 1.00 |
| **4 weeks** | Levo + rifampicin | 1.12 | 0.644 | **1.6e-6** | **0.00** |
| 6 weeks | Levo + rifampicin | 1.12 | 0.644 | 1.1e-12 | 0.00 |
| 12 weeks | Levo + rifampicin | 1.12 | 0.644 | 2.3e-32 | 0.00 |
| 26 weeks | Levo + rifampicin | 1.12 | 0.644 | 1.6e-81 | 0.00 |
| 2 weeks | Levo alone | 0.48 | 0.013 | 6.2e7 | 1.00 |
| 4 weeks | Levo alone | 0.48 | 0.013 | 1.4e8 | 1.00 |
| 6 weeks | Levo alone | 0.48 | 0.013 | 2.7e8 | 1.00 |
| 12 weeks | Levo alone | 0.48 | 0.013 | 1.4e9 | 1.00 |
| **26 weeks** | Levo alone | 0.48 | 0.013 | **4.7e9** | 1.00 |

> A regimen with RSTER > 1 has a **threshold** — 2 weeks is not enough, from 4 weeks
> the reservoir falls below 1 CFU, and beyond that extending it changes nothing.
>
> A regimen with RSTER < 1 **does not respond monotonically to duration** — going
> from 2 weeks to 26 weeks the reservoir **actually grows (6.2e7 → 4.7e9).** The
> duration is not insufficient, it is **irrelevant**. What is needed then is not a
> longer antibiotic course but **the operating theatre again**.

**Implant removed + radical debridement (5.5 log), no rifampicin:**

| Duration | End-of-therapy reservoir | Relapse probability |
|---|---|---|
| 2 weeks | 69.9 | 1.000 |
| 4 weeks | 0.129 | 0.044 |
| **6 weeks** | **1.2e-4** | **0.000** |
| 12 weeks | 2.1e-14 | 0.000 |

**The clinical standard of "6 weeks after radical debridement" is derived in the
model** — 4 weeks is borderline, from 6 weeks it is safe, 12 weeks adds no further
benefit. This table is the structural reason behind the results DATIPO
(PMID 34042388) and OVIVA (PMID 30699315) obtained on the 6-week side and the
duration-shortening side respectively.

### 3.4 The paradox of local antibiotics — the sharpest result in this model

S7: 2 g of vancomycin in cement, no debridement. The wound-fluid concentration really
is very high.

| Day | Wound fluid (mg/L) | Free concentration in sequestrum (mg/L) | × MIC | Biofilm kill (1/d) | Persister kill (1/d) | Biofilm log10 |
|---|---|---|---|---|---|---|
| 0.25 | 613 | 12.3 | 12× | 0.414 | 0.024 | 8.02 |
| **1** | **878** | **34.7** | **35×** | **0.474** | **0.028** | 8.02 |
| 7 | 491 | 19.0 | 19× | 0.445 | 0.026 | 7.90 |
| 28 | 60 | 1.8 | 1.8× | 0.194 | 0.011 | 9.28 |
| 56 | 3.7 | 0.09 | 0.1× | 0.016 | 0.001 | 9.77 |

Even on reaching **35 times** the MIC, biofilm kill at 0.474/d does not exceed the
biofilm growth rate of 0.72/d, and persister kill at 0.028/d is effectively zero. The
biofilm does not shrink.

> **Because the ceiling is not concentration but Emax.** The spatial penalty can be
> bypassed by local administration; the phenotypic penalty cannot. That
> antibiotic-impregnated cement is an **adjunct to** debridement and not a substitute
> for it is, here, a result of the equations rather than of clinical experience.

### 3.5 Biomarkers — the same event, different clocks

End of therapy in the scenarios that cured (S3·S4·S5·S8b·S11): **CRP 8.3 mg/L
(normalised) · ESR 19~22 mm/h (still raised).** In the scenarios that failed:
CRP 126~131 mg/L.

- CRP half-life 19 hours → fast and useful.
- ESR half-life about 10 days → **it stays high for weeks after the patient has
  recovered.** Unsuitable as an index for judging cure.
- And **CRP does not see the reservoir.** It tracks planktonic bacteria and
  cytokines, so it can normalise while the biofilm is still fully in place. The
  divergence of these two curves is the reason "the inflammatory markers have come
  down, so we can stop" is dangerous (Shiny tab ⑦-(b)).

### 3.6 Toxicity — efficacy that requires duration vs toxicity that duration creates

- Vancomycin: eGFR 90 → **65.6** (S2/S3). Includes the feedback in which falling
  renal function raises exposure again.
- Linezolid 6 weeks: platelets 250 → **215** (S5). Duration dependence implemented
  through `TXD50 = 28 d`.
- Rifampicin: ALT 25 → **51** (about 2-fold).

```
T* = argmax [ P_cure(T) − Σ w_k · Tox_k(T) ]
```

Efficacy **requires** duration and toxicity is **made by** duration. Debridement
pushes the `P_cure(T)` curve to the left and thereby reduces `T*` — **surgery is
toxicity reduction.**

---

## 4. Actual errors corrected during model development

The first time this model was run it gave clinically wrong answers. What was fixed is
left on the record.

| Symptom | Cause | Fix |
|---|---|---|
| Bacterial burden diverged to log 38 | the intracellular pool had no carrying capacity | `BICMAX`, proportional to host cell number |
| Osteoclasts exploded into the thousands | the `RANKL/OPG` ratio was unbounded (OPG → 0) | ratio relative to the healthy value `RR0` · capped Hill |
| Radical debridement + vancomycin **always** failed | intracellular SCV proliferated without limit and there was no clearance of infected cells | `FMUIC 0.05→0.005`, `KICCLR` added |
| The film regrew unchanged even after debridement | a fixed floor on the diffusion distance · immune access was a constant | `LBF0` made proportional to EPS maturity · `FIMBFX` |
| **DAIR + FQ alone came out as cure** (the opposite of Zimmerli's control arm) | the local granulocyte defect from the foreign body was missing | `AIMP50` — the term by which the implant cuts immune access |
| **No resistance emerged under rifampicin monotherapy** | biofilm bacteria ate up the planktonic carrying capacity so the mutant could not grow | planktonic carrying capacity consumed by the planktonic pool only |
| Regrowth from 1e-25 CFU after sterilisation | a continuum ODE knows nothing of integer bacteria | growth blocked below 1 CFU, `extf(n) = n/(n+1)` |
| RSTER diverged at carrying capacity | it was divided by the density-limited growth rate | divide by the **intrinsic** growth rate |

The last two matter methodologically. To handle sterilisation with a deterministic ODE
you must (1) prevent growth below 1 CFU and (2) hand everything below that over to
probability. This model hands it over via
`P_relapse = 1 − exp(−PREGROW · N)` (probability that a single organism causes
reinfection, 0.35).

---

## 5. How to run

```bash
# render the mechanistic map
dot -Tsvg com_qsp_model.dot -o com_qsp_model.svg
dot -Tpng -Gdpi=150 com_qsp_model.dot -o com_qsp_model.png

# structural diagnostics + the 14 scenarios + the full duration sweep
Rscript -e 'options(com.autorun = TRUE); source("com_mrgsolve_model.R")'

# an individual scenario
Rscript -e 'source("com_mrgsolve_model.R"); summarise_scn(scn04())'

# Shiny dashboard (10 tabs)
Rscript -e 'shiny::runApp("com_shiny_app.R", port = 8080)'
```

Dependencies: `mrgsolve` (>= 2.0), `shiny`, Graphviz `dot`.
Graphviz requires `fonts-nanum` in order to render Korean.

### Sensitivity of the main parameters

| Parameter | Default | What it governs |
|---|---|---|
| `DEBLOG` | 0~7 log | **The largest single determinant of outcome.** Initial condition P0 |
| `BEMGI` | 1.5 /d | Rifampicin's reason for existing. Set to 0 and no duration succeeds |
| `KWAKE` | 0.08 /d | Persister awakening → the denominator of `T_ster` |
| `AIMP` | 30 cm² | Biofilm carrying capacity + local granulocyte defect |
| `ADEFF`, `AKLS` | 0.004, 2.0 | `phi` → `eta`. The per-drug diffusion penalty |
| `PADSEV`, `FLAPB`, `REVASCB` | 0, 0.2, 0 | Perfusion → drug exposure and immune access simultaneously |
| `FMUT`, `FITCOST` | 1e-8, 0.10 | Pre-existing resistance frequency `FMUT/FITCOST = 1e-7` |
| `FBFEC` | 3.0 | The share of the EC50 displacement attributable to EPS binding (the remainder is mu gating) |

---

## 6. Limitations

- **A single lesion · a single species.** Polymicrobial infection (diabetic foot,
  pressure ulcer) and interspecies interaction are absent.
- **Space is a one-dimensional mean field.** `eta` is the steady-state solution for a
  homogeneous slab, so the irregular geometry of a real sequestrum, multiple lesions,
  and partial perfusion are all smeared into a single `L`.
- **When the reservoir is near 1 CFU the deterministic solution loses meaning.** That
  region is handed over to `P_relapse = 1 − exp(−p·N)`, but what is really needed is
  a stochastic-process simulation. S6b sits exactly in that region (reservoir
  4.8e-7 → deterministically a cure, but in reality a gamble).
- **The parameters are not fitted to an individual patient.** They are illustrative
  values picked from the range of the published literature, matched against clinical
  trial results only to the level of getting the qualitative ordering right.
- **Immunity is condensed onto neutrophils.** Adaptive immunity, antibody, and
  vaccine approaches are absent.
- Bone union and functional recovery are proxied only by `BM` and `PAIN`, which are
  some distance from real reconstructive outcomes.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It was
built on the basis of the published literature but has not been independently
verified or certified, and **must not be used for actual clinical decision-making,
prescribing, or regulatory submission.** The parameters and assumptions are
approximations for explanatory purposes, and fitting and validation against real
patient data are separately required.

The 85 supporting references are collected, with PubMed links, in
[`com_references.md`](com_references.md), and the
mapping from model equation → supporting reference is in the same file.
