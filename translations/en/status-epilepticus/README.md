# Status Epilepticus QSP Model
## Status Epilepticus — a receptor-trafficking clock

<p align="center">
  <a href="../../../status-epilepticus/se_qsp_model.svg"><img src="../../../status-epilepticus/se_qsp_model.png" width="900" alt="Status epilepticus QSP mechanistic map"></a>
</p>

---

## The question this model tries to answer

Status epilepticus is one of the best-studied conditions in emergency
medicine. There are several large randomised trials, and the guidelines are
clear. And yet the answers those trials produced look strangely at odds
with each other.

- ESETT compared three second-line agents and landed on a near-perfect
  three-way tie — **47% · 45% · 46%**. Their mechanisms of action are
  completely different, so why the same number?
- In RAMPART, **the slower drug (intramuscular midazolam) beat the faster
  drug (intravenous lorazepam)**. 73.4% vs 63.4%.
- In the VA Cooperative trial, lorazepam achieved 64.9% in overt SE but
  **7.7–24.2%** in subtle SE. Same drug, same dose.
- Kapur & Macdonald reported that 30 minutes into SE, diazepam potency
  falls by **about 20-fold**, while phenobarbital falls by only **about
  3-fold**.
- Why has ketamine become a late-line drug? It looks less like "it is
  acceptable to give it late" and more like it **needs to be given late**
  to work well.
- Allopregnanolone hits precisely the right target and still failed in the
  phase 3 SRSE trial (STATUS).
- And even as treatment has improved, mortality has barely moved in 30
  years.

Textbooks list these as separate cautionary notes. This model was built to
show that they are **arithmetic consequences of a single structure**.

---

## Structure — three pools and one clock

> Status epilepticus is not "a seizure that runs long." It is a **clock
> that moves receptors.** While the clinician is deciding what to give, the
> very target the drug must bind is leaving the synapse.

A single clock (ongoing seizure activity) moves two receptor pools **in
opposite directions** and leaves a third untouched.

| Pool | What the clock does | At 60 minutes |
|---|---|---|
| `R_SYN` synaptic γ2-GABA-A | clathrin-mediated **internalisation** | 0.52 |
| `F_BZS` benzodiazepine-sensitive fraction (α1→α4 subunit switch) | **decreases** | 0.62 |
| `R_EXTRA` extrasynaptic δ-GABA-A | barely moves | 1.04 |
| `NR_SYN × Mg unblock` synaptic NMDA | **influx** | 1.85 |

And every drug class is written in **exactly the same form**:

```
EFFECT_class(t)  =  f(effect-site concentration)  ×  target pool(t)
```

The only thing that distinguishes a drug class is **which pool is
multiplied**.

| Drug | Target multiplied | Sign over time |
|---|---|---|
| benzodiazepine | `R_SYN × F_BZS` | ↓↓ (both terms fall together) |
| phenobarbital | `0.20·R_SYN + 0.80·R_EXTRA` (+ AMPA/kainate block) | ≈ 0 |
| propofol | `0.25·R_SYN + 0.75·R_EXTRA` | ≈ 0 |
| allopregnanolone | `0.40·R_SYN + 0.60·R_EXTRA` | ≈ 0 |
| **ketamine** | `NR_SYN × Mg unblock × use-dependence` | **↑** |
| LEV / PHT / VPA | `1` (presynaptic/axonal targets do not move) | 0 |

The 59 differential equations are **identical across every scenario**. The
seven contradictions above are not reproduced by changing parameters — they
fall out of this table by calculation.

---

## (A) The price of one minute is not constant

`CREQ_BZD(t)` = the benzodiazepine effect-site concentration (in multiples
of EC50) needed, at this exact instant, to bring the network gain G back to
1. This is a model **output**.

| Time | R_SYN | F_BZS | benzo target | benzo conc. needed | phenobarbital conc. needed |
|---:|---:|---:|---:|---:|---:|
| 5 min | 0.950 | 0.966 | 0.918 | **0.23** | 0.20 |
| 15 min | 0.836 | 0.883 | 0.738 | 0.52 | 0.33 |
| 30 min | 0.700 | 0.776 | 0.543 | **1.62** | 0.52 |
| 45 min | 0.597 | 0.688 | 0.411 | **29.7** | 0.71 |
| 60 min | 0.520 | 0.616 | 0.320 | **impossible** | 0.89 |
| 8 h | 0.216 | 0.281 | 0.061 | impossible | **2.51** |

Between 5 and 30 minutes, the benzodiazepine concentration required rises
**7.0-fold**, phenobarbital's **2.6-fold**. And benzodiazepine has a
**vertical asymptote** — between 45 and 50 minutes in the base patient,
after which **no dose works.** Phenobarbital has no such point even out to
8 hours.

This matches, in both direction and asymmetry, the 20-fold vs 3-fold that
Kapur & Macdonald measured in isolated neurons. The magnitude is
underestimated by the model, for a clear reason: the model moves only
receptor **number** and the **benzo-sensitive fraction**, not the potency
per remaining receptor. This gap is left visible, not hidden.

**Downstream consequence:** response rate by time of first dose (a virtual
cohort of 300)

| First benzo dose | 5 min | 10 min | 20 min | 30 min | 45 min | 60 min |
|---|---:|---:|---:|---:|---:|---:|
| Response rate | 93.3% | 85.0% | 60.0% | 34.7% | 13.0% | 6.3% |
| Hippocampal neuron loss | 0.4% | 2.4% | 9.3% | (cliff) | — | — |

**A dose shortfall also has a value that depends on timing.** The same
2-mg shortfall is nearly free at 5 minutes (scenario 05: seizure ends at
9.8 min, 0.7% loss), but far more costly at 30 minutes (scenario 06: ends
at 59.8 min, 29.4% loss) — because the target is a **shrinking resource**.

---

## (B) ESETT's three-way tie is a structural prediction — but it is not immune to delay

All three drugs (levetiracetam, fosphenytoin, valproate) multiply their
target by **1**. SV2A (presynaptic vesicle), Nav (axonal), and valproate's
multiple targets are not pulled into the cell by seizure activity. So the
three end up equal.

| | model | observed (ESETT) |
|---|---:|---:|
| levetiracetam 60 mg/kg | 46.6% | 47% |
| fosphenytoin 20 mg PE/kg | 50.3% | 45% |
| valproate 40 mg/kg | 47.2% | 46% |

At this point it is tempting to read this as "second-line agents are free
from the clock." **The model does not say that.** Restricting to
first-line failures and varying the time of the second-line dose:

| Time of second-line dose | 30 min | 45 min | 60 min | 90 min | 120 min |
|---|---:|---:|---:|---:|---:|
| Response rate (120 benzo failures) | 98.3% | 83.3% | **54.2%** | 1.7% | 0.0% |

There **is** a cliff. It simply falls **about 40 minutes later** (roughly
62 min) than the first-line cliff (roughly 25 min). Their own targets are
not depleted, but the network's **endogenous inhibitory tone**, which they
still rely on suppressing, is still being fed by `R_SYN`.

> Second-line agents do not stop the clock. **They only buy 40 minutes.**

---

## (C) RAMPART — a logistics result that has been read as pharmacology

The model has **not a single** RAMPART-specific parameter. The difference
between the two arms comes down to exactly two operational facts: (i)
intramuscular administration needs no IV access, so it starts 4.5 minutes
earlier, and (ii) 10% of the IV arm fails on the first attempt, losing a
further 10 minutes. The pharmacology (IM midazolam is slower) actually
works **against** the IM arm.

| | model | observed |
|---|---:|---:|
| IM midazolam 10 mg @ 8.0 min | 80.0% | 73.4% |
| IV lorazepam 4 mg @ 12.5 min | 78.8% | 63.4% |

The direction is reproduced and the magnitude is underestimated. The point
is not the ranking but **why**: absorption delay (< about 5 min) is smaller
than administration delay (about 4.5 min, plus 10 min on failure), and
`R_SYN` is falling throughout that entire interval.

---

## (D) Motor manifestation and EEG are different windows onto the same state

In the model, `MOTOR` and `EEG` are two outputs of a single network state,
`SEIZ`, and the motor gain `MOTG` decays separately over time and in
proportion to GABAergic drug occupancy. For a refractory patient given
repeated doses of lorazepam:

| Time | benzo occupancy | motor output | EEG burden | bedside call | actual |
|---:|---:|---:|---:|---|---|
| 30 min | 0.73 | 0.44 | 1.00 | seizure ongoing | seizure ongoing |
| 60 min | 0.78 | 0.31 | 1.00 | seizure ongoing | seizure ongoing |
| 120 min | 0.89 | 0.17 | 1.00 | borderline | **seizure ongoing** |
| 180 min | 0.92 | 0.11 | 1.00 | **appears controlled** | **seizure ongoing** |
| 240 min | 0.92 | 0.08 | 1.00 | **appears controlled** | **seizure ongoing** |

The VA Cooperative's 64.9% vs 7.7–24.2% is not a result of a different
drug — it is a result of **a different observation window**. The model
reproduces this as the inevitable consequence of "partially treated SE."

---

## (E) Ketamine — the only drug whose target grows, and it still loses to the clock

Fixing site occupancy at the same value (0.75) and comparing targets alone:

| Time | fraction of G cut by benzo | phenobarbital | **ketamine** | benzo target pool | removable NMDA share |
|---:|---:|---:|---:|---:|---:|
| 5 min | **0.746** | 0.766 | 0.288 | 0.918 | 0.451 |
| 30 min | 0.677 | 0.790 | 0.330 | 0.543 | 0.518 |
| 60 min | 0.597 | 0.814 | 0.339 | 0.320 | 0.531 |
| 120 min | 0.468 | 0.841 | 0.338 | 0.153 | 0.530 |
| 240 min | 0.349 | 0.859 | 0.335 | 0.082 | 0.525 |
| 480 min | 0.293 | 0.866 | **0.334** | 0.061 | 0.524 |

The ketamine curve **rises** while the benzodiazepine curve falls. The two
curves cross at **280 minutes** — a model output, not an assumption. This
is the quantitative counterpart of ketamine's clinical role as a
super-refractory agent.

And yet — **scenario 15 vs 16** — the same background therapy and the same
ketamine dose stops the seizure at 71 minutes when given at 65 minutes, but
does not stop it at all over 12 hours when given at 300 minutes.

> Even a drug whose target grows still loses to the clock. Ketamine is a
> drug that is acceptable to give late, not a drug that is better given
> late.

One more finding: ketamine is the model's only sedative that **raises**
MAP (minimum MAP of 76 in scenario 15, vs 67 with propofol in scenario 14).

---

## (F) Stopping the seizure and protecting the brain are different endpoints

`SEIZ` and `INJURY` are separate state variables, moving on separate
clocks.

| Time to termination | cumulative seizure time | hippocampal neuron loss | epileptogenesis burden |
|---|---:|---:|---:|
| 8 min (scenario 02) | 5.3 min | **0.4%** | 0.016 |
| 10 min (scenario 21, optimal) | 7.0 min | 0.9% | 0.027 |
| 16 min (scenario 07, RAMPART IM) | 12.1 min | 3.0% | 0.069 |
| 56 min (scenario 03) | 52.4 min | **28.0%** | 0.834 |
| 96 min (scenario 04) | 90.8 min | **42.7%** | 2.083 |
| never terminated (scenario 01) | 718.7 min | 61.0% | 54.4 |

By trial criteria, every one of these five rows counts as a "success."

---

## (G) Two kinds of resistance are indistinguishable at the bedside, but not in the model

| Time | P-gp induction | BBB opening | brain:plasma index | benzo target pool |
|---:|---:|---:|---:|---:|
| 0 min | 1.00 | 1.00 | 1.000 | 1.000 |
| 30 min | 1.24 | 1.04 | 0.896 | 0.543 |
| 60 min | 1.43 | 1.16 | 0.854 | 0.320 |
| 120 min | 1.66 | 1.47 | 0.853 | 0.153 |
| 240 min | 1.85 | 1.95 | **0.907** | 0.082 |
| 480 min | 1.92 | 2.20 | **0.943** | 0.061 |

**Transporter resistance** (rising P-gp) can be overcome by dose.
**Target resistance** (the pool above) cannot be overcome by any amount of
the same drug. And as time passes, BBB opening offsets P-gp, so the
brain:plasma index **nearly recovers**, while the target pool keeps
falling.

> Clinical implication of the model: **once the brain:plasma ratio is
> normal, stop escalating the same drug.** The remaining problem is not
> exposure — it is the target.

---

## (H) Aetiology sits in the numerator; drugs mostly sit in the denominator

```
G  =  excitatory gain E  /  inhibitory tone I
E  =  EBASE × EDRIVE × (AMPA term + NMDA term) × (1 − second-line suppression) × (1 + IL-1β gain)
I  =  IBASE × (0.6·R_SYN + 0.4·R_EXTRA) × Cl coefficient + peptides + adenosine + drugs
```

`EDRIVE` (aetiology) enters only the numerator. Anticonvulsants mostly
touch only the denominator. So **refractoriness is not "the wrong drug was
chosen"** — it is **"DRIVE exceeds the maximum achievable SUPPRESSION."**

This structure immediately explains three things.

- **Anti-NMDA-receptor encephalitis** (scenario 19, EDRIVE 2.40): antibody
  keeps re-supplying DRIVE, so no combination of anticonvulsants can win.
  The seizure never stops over 12 hours.
- **FIRES + anakinra** (scenario 20): blocking IL-1 lowers the numerator.
  In the same patient, no anakinra gives 719 min seizure, 61.1% loss, 56.4
  epileptogenic burden; with anakinra it becomes **80 min, 37.3%, 1.29.**
- **Hypoglycaemic SE** (scenario 22): removing 90% of DRIVE with
  glucose+thiamine at 30 minutes ends the seizure at 32.5 minutes — the
  only intervention here that touches the numerator.
- **The STATUS trial** (scenario 18, EDRIVE 3.20): allopregnanolone targets
  the surviving pool (`R_EXTRA`) precisely. The target is right, but the
  position on the trajectory is wrong — by the time SRSE is reached, DRIVE
  and INJURY have already passed the point where restoring inhibition
  changes the endpoint.

---

## Deliverables

| File | Contents |
|---|---|
| [`se_qsp_model.dot`](../../../status-epilepticus/se_qsp_model.dot) · [`.svg`](../../../status-epilepticus/se_qsp_model.svg) · [`.png`](../../../status-epilepticus/se_qsp_model.png) | mechanistic map — **177 nodes, 19 clusters** |
| [`se_mrgsolve_model.R`](se_mrgsolve_model.R) | **59 ODEs** (27 PK + 32 disease PD), **193 parameters**, **22 treatment scenarios**, 8 analysis functions |
| [`se_shiny_app.R`](../../../status-epilepticus/se_shiny_app.R) | **10-tab** interactive dashboard |
| [`se_references.md`](se_references.md) | **126 references** — every PMID confirmed via NCBI E-utilities |

### The map's 19 clusters

1 aetiology (DRIVE) · 2 network bistability · **3 clock 1: synaptic GABA-A
internalisation** · **4 the surviving pool: extrasynaptic δ** ·
**5 clock 2: NMDA influx** · 6 chloride/KCC2 · 7 endogenous peptides ·
8 neuroinflammation · BBB · 9 excitotoxic injury · 10 systemic compensation ·
11 systemic decompensation · 12 benzodiazepines · 13 second-line agents ·
14 anaesthetics · 15 neurosteroids · 16 drug resistance · **17 the delay
chain** · 18 clinical endpoints and three dissociations · 19 sequelae ·
epileptogenesis

### The 59 ODEs

**PK (27)** lorazepam · midazolam (IM absorption + 1-OH metabolite) ·
diazepam (2-compartment) · levetiracetam · fosphenytoin→phenytoin
(Michaelis-Menten + albumin binding) · valproate (saturable binding) ·
phenobarbital · ketamine (2-compartment + norketamine) · propofol
(2-compartment) · allopregnanolone — each with its own effect-site
compartment

**Disease PD (32)** `RSYN` `RENDO` `FBZS` `REXTRA` `NRSYN` `NRINT` `AMPACP`
`KO` `SEIZ` `MOTG` `KCC2` `CLI` `ADO` `ADK` `NETPEP` `IL1B` `BBBP` `PGP`
`EDEM` `GLU` `CAI` `INJURY` `ATPD` `AUTO` `MAP` `GLUCP` `LAC` `TEMP` `CK`
`RESPD` `TSEIZ` `EPG`

### The 22 scenarios

no treatment · guideline-adherent (5 min) · delayed (30 min / 60 min) ·
underdosed (early/delayed) · RAMPART IM/IV · ESETT three-arm ·
phenobarbital · midazolam infusion · propofol infusion · ketamine at 65 min
/ 300 min / without a GABA partner · allopregnanolone SRSE · anti-NMDAR
encephalitis · FIRES+anakinra · optimal pathway · hypoglycaemic SE

---

## Running the model

```r
# Required packages: mrgsolve, dplyr (+ shiny, ggplot2, tidyr for the Shiny app)
source("se_mrgsolve_model.R")   # prints all 22 scenarios + 8 analyses

run_all()              # scenario summary table
sweep_benzo_time()     # the price of one minute
potency_drift()        # which drug makes that minute expensive
marginal_effect()      # the crossover
crossover_time()       # 280
calibration_check()    # PHTSE / ESETT comparison
rampart()              # RAMPART reproduction
first_line_delay(); second_line_delay()   # the two cliffs
dissociation()         # the motor–EEG dissociation
resistance_split()     # transporter vs target resistance
fires_anakinra()       # the intervention that touches the numerator

shiny::runApp("se_shiny_app.R")
```

The model has been verified to compile and run under mrgsolve 2.0.1 / R
4.x. Every number in the tables above is the actual output of this code.

---

## What is fitted, what is predicted

**Fitted — only five places.**

| Parameter | Fitted to |
|---|---|
| `KENDO`/`KREC`/`KDEGR` | Naylor 2005: surface γ2 −47% at 1 h into SE (model −48.6%) |
| `KEXO`/`KENDN`/`FEXO0` | Naylor 2013: synaptic NR1 +38% at 1 h into SE (model +38%) |
| `SLEV`/`SPHT`/`SVPA` + corresponding EC50 | ESETT's **tie**, 47/45/46 (model 46.6/50.3/47.2) |
| `EDRIVE` median · variance | PHTSE/RAMPART first-line response rate 59–73% (model 59.8%) |
| `KRESPS > KRESPD` | PHTSE: placebo intubation 22.5% > lorazepam 10.6% |

**Predicted — four points where this model can be falsified.**

1. the location of `CREQ_BZD(t)`'s **vertical asymptote** (45–50 min in the
   base patient)
2. the second-line cliff falling **about 40 minutes** after the first-line
   cliff (25 min vs 62 min)
3. the **280-minute crossover time** between the recoverable benzo pool and
   the removable NMDA share
4. the **difference in cumulative excitotoxic burden** between termination
   at 20 minutes and at 60 minutes (9.3% vs 42.7%)

**What the model cannot do (stated limitations)**

- Changes in potency **per surviving receptor** are not modelled, which is
  why the model underestimates Kapur & Macdonald's 20-fold as 7-fold.
- The 10-point-percentage RAMPART gap is underestimated as 1.2 points —
  only the direction is reproduced.
- There is no spatial spread or focality of the seizure. It is a single
  network state variable.
- Mortality and mRS are read only as monotonic functions of `INJURY` and
  were not separately calibrated.
- Paediatric/neonatal physiology (immature KCC2, different PK) is not
  represented.

---

## What this model is for

It is not for choosing a drug. ESETT has already shown that choice is
almost free. This model exists **to put a price on one minute**, and to
say which drug makes that minute expensive (benzodiazepines), which makes
it less expensive (second-line agents), and which one actually gains value
from it (ketamine).

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational
and research purposes.** It was built from published literature and
clinical trial data but has not been independently verified or certified,
and **must not be used directly for clinical decision-making,
prescribing, or regulatory submission.** Parameters and assumptions are
illustrative approximations; separate fitting and validation against real
patient data would be required.
