# Acute Bacterial Meningitis (pneumococcal) — QSP Model

> **The moment at which the treatment becomes the cause of the injury, and the fact that the shield closes the door along with it.**
>
> Bacterial meningitis is usually taught as "kill the organism fast and give a
> steroid alongside". This model starts instead from **a single product**.
>
> ```
>      injury flux  =  k_kill(t) × N(t) × Y_lysis × (1 − E_dex(t))
> ```
>
> And it adds just one fact: **k_kill jumps from 0 to E_max at the first antibiotic
> dose, and N is at its maximum at precisely that moment.**
>
> So this product is **at its maximum in the first few hours of treatment.** The
> antibiotic saves the brain and, at that same moment, pours out the largest
> quantity of inflammatory cargo.

---

## The three sentences this model asserts

### ① The product — what is dangerous is not the number of organisms but the rate of killing

These are the values the simulation prints hour by hour (S02: ceftriaxone alone, no
steroid):

```
 t[h]  CEF_CSF  C/MIC  kkill/h  log10Nc  lysisflx   CW    PLY    Mg    TNF
  0.5     0.82     27     1.37     6.97     12.78   3.9   1.90  0.69   530
  1.0     2.57     86     1.40     6.79      8.74   9.3   4.14  0.82   627
  2.0     4.51    150     1.40     6.45      3.93  14.9   5.57  0.93   729
  3.0     5.24    175     1.40     6.10      1.77  17.0   5.17  0.95   728
  4.0     5.36    179     1.40     5.75      0.80  17.4   4.25  0.95   683
  8.0     3.99    133     1.40     4.39      0.03  15.0   1.36  0.93   515
 24.0     4.33    144     1.40     1.55      0.00   5.9   0.01  0.86   408
```

**The cell-wall load rises 13-fold, from 1.3 CWU/mL before dosing to 17.4 three to
four hours later.** Pneumolysin passes its peak at 2 hours. TNF goes from 530 to a
maximum of 737 pg/mL. The organisms disappear while the inflammation reaches its
maximum immediately afterwards — the "rise in CSF cytokines after antibiotic
administration" recorded in the literature is the arithmetical consequence of this
product.

### ② Delay is bad not because "there are more organisms" but because "there is more to kill"

The same antibiotic, with only the starting time changed:

| Antibiotic delay | AUC_lysis | Peak lactate | Time to sterilisation | Hearing loss dB | Mortality % |
|---|---|---|---|---|---|
| 0 h | **22.3** | 5.2 | 30.4 h | 0.1 | 11.3 |
| 3 h | 116.7 | 5.2 | 40.2 h | 0.2 | 16.6 |
| 6 h | 313.6 | 7.3 | 47.7 h | 8.9 | 20.5 |
| 12 h | **610.4** | 8.4 | 58.1 h | 19.2 | 25.2 |

**Twelve hours late and the total cargo poured out is 27 times larger.** Because N
has grown at μ ≈ 0.85 /h, and because all of that N has in the end to be killed.
"Give it fast" and "kill it slowly" do not conflict — the optimum is **as fast as
possible with the shield already raised**.

### ③ The two doors are the same door — and the sign differs from drug to drug

The blood-CSF barrier permeability `Pb` is opened by inflammation and closed by
dexamethasone. When it opens, hydrophilic antibiotics come in — and albumin,
neutrophils and water come in at the same time. So the steroid helps with one arm
and hinders with the other. **Which of the two wins is decided by that antibiotic's
C/MIC headroom.**

| Scenario | VAN CSF AUC | Peak Pb | **CSF sterilisation** | Mortality % |
|---|---|---|---|---|
| Susceptible organism · CEF+VAN+DEX (standard therapy) | 813 | 6.1 | **20.2 h** | 11.2 |
| Resistant organism (MIC 4) · CEF+VAN, no DEX | 1539 | 8.7 | **73.0 h** | 43.0 |
| Resistant organism · CEF+VAN**+DEX** | 1310 | 7.0 | **160.2 h** | 36.1 |
| Resistant organism · CEF+VAN+**RIF**+DEX | 789 | 6.1 | **39.9 h** | 12.0 |
| Resistant organism · CEF alone | — | 8.7 | **fails** | 78.6 |

How to read it:

- In a cephalosporin-resistant organism (MIC 4 mg/L), **when dexamethasone tightens
  the barrier, vancomycin CSF exposure falls by 15 % and sterilisation is delayed
  2.2-fold, from 73 h to 160 h.** The same manoeuvre does no harm at all to
  ceftriaxone against a susceptible organism — because there the C/MIC headroom is
  200-fold. **Same drug, same dose, opposite sign.**
- And that harm **disappears if a lipophilic drug is added**: rifampicin's
  penetration is almost independent of Pb (a = 0.1), so it gets in even when the
  steroid closes the door. Sterilisation 39.9 h, mortality 36.1 % → 12.0 %.
- In the resistant organism, ceftriaxone alone exceeds the bactericidal threshold
  (4×MIC = 16 mg/L) for only 75.8 h out of 336 and fails to sterilise. **The
  clinical rationale for combining vancomycin comes out as an output rather than as
  an assumption.**

---

## The same log-kill, a different injury integral

β-lactam killing is lytic, so the cargo yield per organism `Y` is at its maximum
(1.00); rifampicin is non-lytic, so it is 0.15. Give rifampicin two hours ahead:

| | Time to sterilisation | AUC_lysis | Peak CW | Peak PLY | Peak TNF | Mortality % |
|---|---|---|---|---|---|---|
| CEF alone | 28.5 h | 22.3 | 17.4 | 5.58 | 737 | 26.9 |
| **RIF 2 h ahead → CEF** | **17.5 h** | **16.6** | **13.4** | **4.75** | **699** | 25.6 |

**It sterilises faster and the cargo is 26 % smaller.** Because the amount killed is
the same but the yield is different.

---

## What this model says about steroids (and what it cannot say)

| | AUC_TNF | Peak PMN | Peak Pb | Hearing loss dB | Mortality % |
|---|---|---|---|---|---|
| No DEX | 39,800 | 3,364 | 8.4 | 0.3 | **26.9** |
| DEX 20 min before | 19,400 | 1,847 | 6.1 | 0.1 | **11.9** |
| DEX simultaneously | 18,700 | 1,849 | 6.1 | 0.1 | 11.3 |
| DEX +4 h | 20,100 | 2,327 | 6.4 | 0.1 | 11.7 |
| DEX +12 h | 22,100 | 2,960 | 6.9 | 0.1 | 12.9 |
| DEX for 1 day only (instead of 4) | 24,700 | 1,868 | 6.5 | 0.2 | **16.8** |

**The effect of giving it at all is large** (26.9 % → 11.9 %, AUC_TNF −51 %).
**The duration is clear too** (4 days 11.9 % vs 1 day 16.8 %).
And yet **the effect of the timing is almost nil in this model.** Rechecked in a
severely ill patient (inoculum 10⁸ · antibiotic delayed 6 h), it is 21.8 % at
−20 min vs 21.8 % at +4 h vs 24.1 % at +24 h.

The reason is structural. The cytokine burst is 2–6 hours wide, whereas the injury
integrals for hearing, cognition and cortex are on a scale of tens to hundreds of
hours, so the peak accounts for only a small share of a 14-day integral. **This is
not a convenient conclusion but a falsifiable prediction.** The most likely
candidate for the missing mechanism is **a threshold in the injury terms**: every
injury term in this model is first-order in its driver, so "a twofold higher peak
for 3 h" and "a 1.06-fold higher level for 100 h" produce the same injury. If real
hair-cell and neuronal death does have a threshold, the peak becomes
disproportionately expensive, and only then does "before the antibiotic" become
arithmetically important. Putting a threshold in and reproducing the recommendation
would be easy, but that is choosing the structure to suit the conclusion, so it is
not done here and is left as a prediction.

---

## Who actually eats the glucose

Low CSF glucose — hypoglycorrhachia — is usually explained by "the bacteria eating
the glucose". In this model the two consumption terms are each computed
independently:

- Bacteria: from a pneumococcal dry weight of 0.3 pg/cell and a homolactic
  fermentation yield of ≈ 25 g/mol → **2.0 mg/dL/h** per 10⁷ CFU/mL
- Neutrophils: from the reported glycolytic rate → **20 mg/dL/h** per 4,000 /µL

**The neutrophils eat ten times as much.** That is why the low glucose persists
after the CSF has become culture-negative (S24: sterilisation at 20.2 h, glucose
44.5 mg/dL at 48 h and 51.9 at 120 h), and why **CSF glucose cannot serve as an
early marker of sterilisation.** This is not an assumption but a result that follows
from the ratio of the two coefficients.

And because GLUT1 is a **bidirectional transporter**, the further CSF glucose falls
the larger the net influx becomes. That is why glucose never reaches 0 in a treated
patient.

---

## Recovery has to be computed too (S24, standard therapy)

To confirm that the hazard function is not accruing chronically, the recovery
trajectory has actually to return to normal:

```
 t[h] log10Nc  PMN/µL prot mg/dL   gluc  lactate  TNF   Mg   Pb R_out   ICP   CPP   MAP  SOFA
    0    6.97    1283      180   35.4      4.8   530 0.69  6.0  0.50  16.9  69.0  85.9  2.15
    4    5.51    2017      180   37.6      3.7   682 0.95  6.3  0.65  19.7  54.0  73.7  3.04
   12    2.37    2374      185   33.2      4.0   239 0.91  6.3  0.87  23.8  39.9  63.7  2.80
   24    0.00    1691      173   39.0      3.5   164 0.86  5.6  0.97  25.7  37.8  63.6  1.65
   48    0.00    1204      146   44.5      2.9   144 0.69  4.8  0.98  25.8  39.0  64.7  1.09
  120    0.00     536      109   51.9      2.2    50 0.07  3.6  0.92  24.7  46.5  71.2  0.53
  240    0.00       1       50   59.6      1.7     0 0.00  1.6  0.71  20.9  66.4  87.3  0.01
  336    0.00       0       31   60.0      1.7     0 0.00  1.0  0.58  18.1  69.9  88.0  0.00
```

The organisms are gone by 24 hours, the barrier closes from 6.0 to 1.0, ICP goes
25.7 → 18.1 mmHg, CPP 37.8 → 69.9 mmHg, and SOFA goes to 0.
**If this trajectory is not right, none of the endpoints below means anything.**

---

## What actually kills people (decomposition of the hazard)

```
Haz = h0·T + h_icp·I_icp + h_cpp·I_cpp + h_sofa·I_sofa + h_isch·I_isch
      + h_bact·I_bact + h_ncsf·I_ncsf   (acute)   +   h_cort·(1 − N_cort)   (structural)
```

| Scenario | ICP | CPP | SOFA | Ischaemia | Bacteraemia | CSF infection | Structural | Total | Mortality % |
|---|---|---|---|---|---|---|---|---|---|
| No treatment | 0.006 | 0.644 | 0.119 | 0.010 | 0.106 | 0.457 | 0.325 | 1.671 | 81.2 |
| CEF alone | 0.002 | 0.168 | 0.015 | 0.002 | 0.001 | 0.015 | 0.108 | 0.314 | 26.9 |
| CEF+DEX | 0.000 | 0.062 | 0.011 | 0.000 | 0.001 | 0.016 | 0.035 | 0.128 | 11.9 |
| 12 h delay | 0.001 | 0.149 | 0.021 | 0.000 | 0.005 | 0.042 | 0.068 | 0.290 | 25.2 |
| Resistant organism+VAN+DEX | 0.001 | 0.190 | 0.035 | 0.001 | 0.011 | 0.115 | 0.091 | 0.448 | 36.1 |
| Standard therapy | 0.000 | 0.060 | 0.010 | 0.000 | 0.001 | 0.011 | 0.035 | 0.121 | 11.2 |
| Standard therapy+EVD | 0.000 | 0.004 | 0.010 | 0.000 | 0.001 | 0.011 | 0.031 | 0.061 | 6.2 |

**Perfusion (CPP) is the largest term in every scenario, and structural cortical
injury comes next.** Meningitis kills through the brain, not through the blood.
Because this decomposition was built as state variables, the hazard coefficients did
not have to be matched by eye — they could be **solved arithmetically against the
clinical-trial targets** (below).

---

## Comparison with the clinical trial — a virtual cohort of 10 × (DEX with/without)

| Patient | Inoculum log₁₀ | Delay h | Host defence | Mortality %(DEX−) | Mortality %(DEX+) | Hearing dB(−) | Hearing dB(+) |
|---|---|---|---|---|---|---|---|
| 1 | 6.0 | 0 | 1.00 | 15.9 | 6.8 | 0.1 | 0.1 |
| 2 | 6.5 | 2 | 1.00 | 27.0 | 11.4 | 0.3 | 0.1 |
| 3 | 7.0 | 1 | 0.90 | 29.4 | 12.6 | 2.3 | 0.2 |
| 4 | 7.5 | 4 | 0.85 | 38.4 | 19.5 | 21.7 | 8.5 |
| 5 | 8.0 | 6 | 0.70 | 41.0 | 22.3 | 27.4 | 15.5 |
| 6 | 7.0 | 3 | 0.55 | 36.1 | 17.5 | 16.7 | 1.7 |
| 7 | 6.7 | 1 | 1.00 | 26.2 | 10.8 | 0.3 | 0.1 |
| 8 | 7.3 | 8 | 0.60 | 41.5 | 23.0 | 28.1 | 16.2 |
| 9 | 7.9 | 2 | 0.80 | 37.8 | 18.9 | 20.8 | 7.3 |
| 10 | 7.5 | 24 | 0.40 | 48.7 | 32.5 | 43.4 | 33.6 |
| **Mean** | | | | **34.2** | **17.5** | **16.1** | **8.3** |

| Measure | Model | Literature (pneumococcal, adults) |
|---|---|---|
| Mortality, no steroid | **34.2 %** | 34 % (de Gans & van de Beek, EDS 2002) |
| Mortality, with steroid | **17.5 %** | 14 % |
| Unfavourable outcome, without → with | **50.9 % → 26.1 %** | 52 % → 26 % |
| Any hearing loss (>25 dB) | **30 % → 10 %** | 20–30 % |

The hazard coefficients were determined from this comparison **arithmetically**:
because the component-wise integrals are held as state variables, a post hoc
decomposition `Haz = h0·T + Σ hᵢ·Iᵢ` is available, and solving 1.256·s = 0.416 for
no-DEX and 0.538·s = 0.151 for DEX against the cohort means gives s = 0.331 / 0.281.
The two values diverge because the model's hazard ratio (2.33) is slightly smaller
than the target (2.76), so the compromise value **s = 0.30** was applied to every
coefficient.

**What does not fit is written down as well**: the untreated natural history comes
out at 81.2 % mortality at 14 days, whereas in reality it is close to 100 %. This is
because the model has no state that absorbs brainstem death and herniation, and it
is an explicit limitation. Severe hearing loss (>60 dB) also appears in only one
member of the cohort (the patient who presented late), fewer than reported.

---

## Structural consistency that can be checked by hand

The things that have to come out right before the model is run — these are
identities, not fits:

| Check | Model | Literature |
|---|---|---|
| ICP = P_ss + Q_f·R_out (normal R_out) | 9.5 mmHg | normal 7–15 |
| The same equation, R_out ×7 | 31.2 mmHg | severe meningitis 25–40 |
| Ceftriaxone free fraction (30 → 250 mg/L) | 0.076 → 0.167 | saturable binding, 0.05 → 0.20 |
| CEF total-concentration CSF penetration (Pb 1 → 9) | 0.021 → 0.094 | 0.015–0.10 |
| VAN total-concentration CSF penetration (Pb 1 → 9) | 0.022 → 0.144 | 0.01–0.15 |
| RIF total-concentration CSF penetration | 0.118 → 0.144 | 0.10–0.20, independent of inflammation |
| DEX total-concentration CSF penetration | 0.193 → 0.235 | 0.15–0.30 |
| CSF glucose balance at 60 mg/dL | influx 9.3 = bulk 8.4 + brain 1.0 | — |
| RK4 stability condition (ICP eigenvalue 33.1/h) | dt < 0.084 h, 0.04 used | — |

In the integration-convergence check (dt halved), the relative differences in peak
ICP, hearing loss, probability of death and AUC_lysis are all below 0.54 %.

---

## Files

| File | Contents |
|---|---|
| [`abm_qsp_model_en.dot`](abm_qsp_model_en.dot) · [`.svg`](abm_qsp_model_en.svg) · [`.png`](abm_qsp_model_en.png) | Mechanistic map — 164 nodes · 21 clusters · 146 edges |
| [`abm_mrgsolve_model_en.R`](abm_mrgsolve_model_en.R) | 63-ODE mrgsolve model · 26 treatment scenarios · virtual cohort |
| [`abm_reference_python_en.py`](abm_reference_python_en.py) | A dependency-free independent Python RK4 implementation of the same 63 equations (for verification) |
| [`abm_reference_output_en.txt`](abm_reference_output_en.txt) | The full log of that run — the source of every number in this README |
| [`abm_shiny_app_en.R`](abm_shiny_app_en.R) | 12-tab interactive dashboard |
| [`abm_references_en.md`](abm_references_en.md) | 650 references (25 sections) — nothing but the PubMed API responses, transcribed |
| [`fetch_refs_en.py`](fetch_refs_en.py) | The script that generates that reference list (NCBI E-utilities) |

### The order in which to read the map (21 clusters)

1 host risk factors · 2 colonisation → invasion → entry into the CSF ·
3 virulence factors (the cargo to be released) ·
4 subpopulations and the killing term (where the product is made) · 5 innate immune
recognition · 6 cytokines · 7 neutrophils and proteolysis · 8 the blood-CSF barrier
(the two-way door) · 9 CSF dynamics and ICP ·
10 perfusion · autoregulation · ischaemia · 11 CSF metabolism · 12 neuronal injury
(two branches) · 13 antibiotic PK (the delivery stages) · 14 antibiotic PD
(headroom · lysis · the resistance cliff) ·
15 dexamethasone (the shield and the closing door) · 16 adjunctive therapy ·
17 cochlear injury · 18 cerebrovascular complications · 19 systemic sepsis ·
20 diagnosis and endpoints · 21 legend

---

## Model structure (63 ODEs)

| Compartments | State variables |
|---|---|
| 1–3 | ceftriaxone: central · peripheral · CSF |
| 4–6 | vancomycin: central · peripheral · CSF |
| 7–8 | rifampicin: central · CSF |
| 9–11 | dexamethasone: plasma · CSF · **the asymmetric transcriptional-effect compartment** |
| 12–15 | osmotherapy: mannitol · glycerol (gut · plasma) · accumulated intracerebral osmoles |
| 16–20 | bacteria: free in CSF · adherent/sequestered · in blood · cell wall · pneumolysin |
| 21–30 | macrophage activation · TNF · IL-1β · IL-6 · IL-10 · CXCL8 · complement · neutrophils · MMP-9 · ROS |
| 31–40 | barrier Pb · CSF albumin · protein · glucose · lactate · CSF accumulation · outflow resistance · **ICP** · brain water · autoregulation |
| 41–44 | MAP · temperature · SOFA · volume status |
| 45–50 | cortical neurons · hippocampal dentate gyrus · cochlear hair cells · labyrinthine ossification · seizure burden · acute hazard |
| 51–57 | exposure trackers (AUC_cef, AUC_van, AUC_lysis, AUC_TNF, T>4×MIC ×2, AUC_ICP) |
| 58–63 | component-wise integrals of the hazard (I_icp, I_cpp, I_sofa, I_isch, I_bact, I_ncsf) |

A few of the central mechanisms:

- **The asymmetric transcriptional-effect compartment**: glucocorticoid suppression
  switches on fast (t½ 1.5 h) and off slowly (t½ 17 h). This asymmetry is why a
  4-day course is sufficient, and at the same time it is the structure that
  (according to ①) ought to create a timing window.
- **The adherent/sequestered subpopulation**: it receives only 35 % of the killing
  effect, so it always disappears late. Sterilisation of free organisms in the CSF at
  20.2 h against clearance of the adherent subpopulation at 30.2 h — that gap is why
  a full course of treatment is needed.
- **Exponential intracranial compliance**: dP/dV = K·P (Marmarou PVI 25 mL). At
  steady state it reduces to the Davson relation ICP = P_ss + Q_f·R_out, and R_out
  was calibrated with that.
- **The efficacy of osmotherapy ∝ 1/Pb**: the barrier damage that made the oedema
  also lets the osmotic agent leak.

---

## The 33 defects caught by numerical verification

There was no R runtime in this environment, so the 63 equations were
**independently implemented in dependency-free Python RK4 and actually integrated**.
Thirty-three defects came out that would not otherwise have been found. The full
list is at [F1]–[F33] in the header of
[`abm_reference_python_en.py`](abm_reference_python_en.py), and each is left as a comment at the
place where it was fixed. The ones that changed the structure:

| # | Defect | Symptom |
|---|---|---|
| **F1** | The glucose consumption terms had no substrate dependence | Consumption continued after CSF glucose had reached 0, giving **lactate of 243 mmol/L** (measured 6–12). Lactate comes only from glucose |
| **F14** | The bacteria in blood had no carrying capacity | N_b diverged to **2×10¹⁰ CFU/mL** by day 14, and that then re-contaminated the CSF |
| **F15** | A molecular weight of 18.2 instead of 182 in the osmolality conversion | Mannitol 0.5 g/kg raised plasma osmolality by **113 mOsm/kg** (actual 11.3). A tenfold unit error |
| **F16** | The hazard function multiplied permanent damage by time | Even a recovered patient accrued hazard throughout the 14 days → **22 of the 26 scenarios had a probability of death above 99 %**, making every treatment comparison meaningless |
| **F17** | The autoregulation formula made autoregulation harmful | At CPP 38 a patient with autoregulation intact was **more** ischaemic than a pressure-passive one (0.54 vs 0.07) |
| **F25** | Cortical loss progressed chronically through the ROS term | Even an optimally treated patient lost 56 % of the cortex |
| **F28** | The ICP and CPP hazard terms were linear | "CPP 45 for 100 h" and "CPP 5 for 10 h" carried the same hazard → a mild but sustained excursion dominated the hazard function |
| **F29** | An extinction floor must not be written inside the derivative block | A discontinuous switch stalls the integrator. The mrgsolve version leaves the state untouched and decides by a reporting threshold alone |
| **F30** | The hazard coefficients were being matched by eye | The component-wise integrals were added as state variables and solved arithmetically against the clinical-trial targets (a single scale factor of 0.30) |
| **F33** | **Doses scheduled at t < 0 silently disappeared** | The actual first dose in the "dexamethasone 20 min before" arm was at **+5.67 h** — that is, the "given first" arm was receiving it later than the "+4 h" arm, and **the central claim of this model was in an untestable state.** The evidence was in the output: the peak cell wall, pneumolysin and TNF of the "rifampicin 2 h ahead" scenario matched the scenario without rifampicin to two decimal places (17.4 / 5.58 / 737) — because no rifampicin had been received at all |

F33 is the most important item on this list. Before it was fixed, "the timing of
dexamethasone does not matter" looked like a conclusion of the model, whereas it was
not physiology but a scheduling bug. After the fix the conclusion survived (with its
grounds completely changed), and it is left, as recorded above, as a **falsifiable
prediction**.

---

## Usage

```r
source("abm_mrgsolve_model_en.R")

structural_checks()          # the identities that can be checked by hand
claim1_burst()               # the peak of the product (the first 24 h)
claim2_timing()              # dexamethasone timing sweep
claim3_signflip()            # the four arms where the sign differs
claim4_delay()               # antibiotic delay sweep
hazard_decomposition()       # what actually kills
run_all_scenarios()          # the 26 scenarios
run_cohort()                 # virtual cohort, 10 patients × 2

shiny::runApp("abm_shiny_app_en.R")     # the 12-tab dashboard
```

In an environment without R, the reference implementation can be run just as it is:

```bash
python3 abm_reference_python_en.py      # no dependencies. regenerates abm_reference_output_en.txt
```

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It was
built from the public literature and clinical-trial data but has not been
independently validated or certified, and **must not be used for real clinical
decision-making, prescribing, or regulatory submission.** The limitations recorded
above (under-prediction of untreated mortality, the missing tail of severe hearing
loss, the absence of a threshold in the injury terms) are only those that are known;
there will be more that are not.
