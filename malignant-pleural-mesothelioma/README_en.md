# Malignant Pleural Mesothelioma — QSP Model

> **One-line summary.** This model writes mesothelioma as a **sheet, not a ball**.
> Write the tumour burden as `V = A × h` (invaded pleural area × rind thickness) and
> the fact that mRECIST sees **only one factor** of that product becomes a computable
> result — even in the untreated state the imaging reads as **−14.7 %** (a response)
> while viable tumour cells have gone up **+44 %**.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (231 nodes · 18 clusters) | [`mpm_qsp_model.dot`](mpm_qsp_model.dot) · [SVG](mpm_qsp_model.svg) · [PNG](mpm_qsp_model.png) |
| ⚙️ mrgsolve model (51 ODEs · 24 scenarios) | [`mpm_mrgsolve_model.R`](mpm_mrgsolve_model.R) |
| 🐍 Reference implementation + verification (dependency-free Python RK4) | [`mpm_reference_model.py`](mpm_reference_model.py) · [output](mpm_calibration_output.txt) |
| 📊 Shiny dashboard (10 tabs) | [`mpm_shiny_app.R`](mpm_shiny_app.R) |
| 📚 References (188 items, actually queried against PubMed) | [`mpm_references.md`](mpm_references.md) · [fetch script](mpm_fetch_references.py) |

<p align="center">
  <a href="mpm_qsp_model.svg"><img src="mpm_qsp_model.png" width="760" alt="MPM QSP mechanistic map"></a>
</p>

---

## 1. Where this model departs from other tumour models

Almost every quantitative account of a solid tumour **assumes the tumour is a ball**.
Volume is the cube of a diameter, response is the change in that diameter, doubling
time is read off the volume, and drug arrives from the vessels surrounding the mass.

Pleural mesothelioma is not a ball. It is a **sheet (rind) covering up to ~1300 cm² of
the inner surface of a hemithorax at a thickness of 0.3–3 cm**. Writing the burden as

```
V(t) = A(t) × h(t)        A = invaded pleural area [cm²]
                          h = mean rind thickness [cm]
```

is not a matter of notation; it produces **four consequences**. Each is a particular
equation in this model, not a sentence.

### (a) The measurement is one factor of the product

modified RECIST is **the sum of six mutually perpendicular thicknesses**. It sees `h`
and cannot see `A`. And `h` is not viable tumour — the fibrous and necrotic stroma left
behind by dead cells is measured just the same, so

```
h_meas = (N + φ·M) / A          φ = 0.70
```

Two things follow. First, as the sheet **spreads it thins itself**. Second, the better
the treatment works, the larger the share of the measured thickness that is stroma
rather than tumour (**62 %** at cycle 6).

Model output (two columns for the same patient at the same time point):

| Scenario | Best mRECIST change | Change in viable cells at that point | Difference | Assessment |
|---|---:|---:|---:|:--:|
| Untreated (BSC) | **−14.7 %** | **+44.1 %** | +58.8 pp | SD |
| Cisplatin alone | **−32.4 %** | **+1.1 %** | +33.5 pp | **PR** |
| Cisplatin + pemetrexed | −59.1 % | −60.7 % | −1.6 pp | PR |

The imaging of an untreated patient reads as a "response", and cisplatin alone, which
in fact kills almost no tumour, **passes the partial-response criterion**. No arbitrary
assumption was put in; it follows from nothing more than the fact that `h_meas` is a
quantity divided by `A`.

### (b) Growth is a front, not a mass

A sheet advances at its **perimeter** (`P = 2√(πA)`). The proliferating pool is a shell
of fixed depth (`L_OX = 0.18 mm`) and is therefore **proportional to area**, whereas the
burden is proportional to **area × thickness**. As a result **the volume doubling time
lengthens by itself, with no resistance mechanism at all**.

| Day | Area cm² | Thickness cm | Volume cm³ | Doubling time d | f_exp | Effusion mL | SMRP nM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 450 | 0.80 | 360 | — | 0.095 | 1000 | 4.0 |
| 120 | 724 | 0.69 | 499 | 208 | 0.108 | 1320 | 6.6 |
| 365 | 1097 | 0.97 | 1069 | 263 | 0.076 | 935 | 14.7 |
| 730 | 1267 | 1.42 | 1803 | 641 | 0.052 | 491 | 24.5 |
| 1000 | 1292 | 1.63 | 2109 | **1196** | 0.045 | 318 | 28.5 |

### (c) Drug delivery arrives from one face only — and the front is thin

The rind is perfused from **one face only**, the one carrying the chest-wall and
visceral vessels; the face towards the pleural cavity is avascular. The mean relative
exposure of a slab of thickness `h` supplied from one side with penetration length `L`
is

```
fpen(L, h) = (L/h) · (1 − e^(−h/L))        1 if h≪L, L/h if h≫L
```

This one function carries the whole delivery-geometry argument.

| Thickness h | L_p (mm) | Systemic f | Intrapleural f | Union | **Sanctuary** | IgG f | T-cell f |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.20 cm | 0.92 | 0.407 | 0.150 | 0.556 | 0.444 | 0.122 | 0.182 |
| 0.80 cm | 0.76 | 0.095 | 0.037 | 0.132 | **0.868** | 0.025 | 0.038 |
| 2.40 cm | 0.70 | 0.029 | 0.013 | 0.042 | **0.958** | 0.008 | 0.012 |

That is, in a typical rind at the time of diagnosis **87 % of the cells sit in a
geometric sanctuary that no dose intensity reaches**. This is built into the model as
the reason a complete response is essentially never observed in mesothelioma. A 150 kDa
antibody and a T cell do worse.

But **the advancing margin is only 0.6 mm thick**, so there `f_sys = 0.691` and
`f_Tcell = 0.436`. **The front is reached and the mass is not.** So chemotherapy
**stops circumferential spread** long before it shrinks the mass — and what the hazard
is most sensitive to is not the mass but the spread (encasement).

Applying the same regimen to a thin and to a thick rind:

| Baseline state | f_exp(0) | Viable-cell nadir | Best mRECIST | Median survival |
|---|---:|---:|---:|---:|
| Thin 0.35 cm / 250 cm² | 0.257 | −93.7 % | −64.8 % | **15.3 mo** |
| Usual 0.80 cm / 450 cm² | 0.095 | −76.8 % | −59.1 % | 12.0 mo |
| Bulky 2.00 cm / 800 cm² | 0.034 | −51.1 % | −53.9 % | **4.4 mo** |

### (d) A sheet cannot produce a resection margin

"Macroscopic complete resection" is, for a sheet, **R1 by definition**. It collapses `h`
but does not touch `A`, and the microscopic sheet that remains regrows from the same
area in a wound flooded with IL-6 and TGF-β.

```
pre-op   day  69: area 461 cm²,  thickness 0.566 cm,  volume 261 cm³
post-op  day  72: area 461 cm²,  thickness 0.022 cm,  volume  10 cm³
              → 96.2 % of the volume removed,  0.0 % of the area removed
         day 300: area 783 cm²,  thickness 0.475 cm,  volume 372 cm³
         day 600: area 1164 cm², thickness 1.038 cm,  volume 1208 cm³
```

**MARS2 is an output of the model.** The model hazard ratio for extended
pleurectomy/decortication + chemotherapy against chemotherapy alone is **1.27** at 12
months and **1.35** at 24 months — the MARS2 reported value is **1.28**. Extrapleural
pneumonectomy gives **2.15** (MARS reported 1.90).

### Fifth — the pleural cavity is a compartment with a volume

The effusion dilutes intrapleural drug, and the rind **blocks the lymphatic stomata**
that drain the effusion. Formation rises with area, drainage falls with area, and the
space is finite, so **the effusion is not monotonic in tumour area**:

```
   day 0:  1000 mL  (A/S = 0.35)
 day 131:  1322 mL  (A/S = 0.57)   ← peak
 day 365:   935 mL  (A/S = 0.84)
day 1000:   318 mL  (A/S = 0.99)   ← the late hemithorax dries out
```

No rule was added in order to produce the "drying out" of the thorax at the end of the
disease. It is a consequence of the product.

---

## 2. Calibration to the clinical trials

**The control arms of the MPM trials differ from one another by more than four months**
(EMPHACIS cisplatin 9.3 months, CheckMate 743 chemotherapy 14.1 months, MARS2
chemotherapy 24.8 months). Absolute median survival is therefore not a common scale.
This model was **anchored on three absolute values and calibrated on within-trial
hazard ratios for everything else**. Every observed value was read directly from the
PubMed abstract of the paper concerned.

**Absolute anchors** (cohort = epithelioid 75 % / biphasic 13 % / sarcomatoid 12 %)

| Arm | Model | Observed | Source |
|---|---:|---:|---|
| Best supportive care | **7.7 mo** | 7.6 mo | MS01 ASC arm |
| Cisplatin alone | **9.0 mo** | 9.3 mo | EMPHACIS control arm |
| Cisplatin + pemetrexed | **12.0 mo** | 12.1 mo | EMPHACIS |

**Hazard ratios** (the model is read at 24 months)

| Comparison | Model HR 12mo | Model HR 24mo | Trial HR | Source |
|---|---:|---:|---:|---|
| cis+pem vs cis | 0.57 | **0.71** | 0.77 | EMPHACIS |
| + bevacizumab | 0.83 | **0.81** | 0.77 | MAPS |
| nivolumab+ipilimumab vs chemo | 1.52 | **0.78** | 0.74 | CheckMate 743 |
| pembrolizumab+chemo vs chemo | 0.88 | **0.82** | 0.79 | IND.227 / KEYNOTE-483 |
| second-line nivolumab vs placebo | 0.98 | **0.89** | 0.69 | CONFIRM |
| pegargiminase+chemo (non-epithelioid) | 0.68 | **0.74** | 0.71 | ATOMIC-Meso |
| extended P/D + chemo vs chemo | 1.27 | **1.35** | 1.28 | MARS2 |
| extrapleural pneumonectomy | 2.16 | **2.15** | 1.90 | MARS |

Six of the eight hazard ratios are within ±0.08 of the observed value.

---

## 3. Results that came out of the model

### 3.1 Histology — why chemotherapy and immunotherapy change places

`CHEMOS` (chemosensitivity) and `IMMINF`/`VISTA_S` (immune infiltration and
suppression) were **set from the pathology** and not fitted back from the outcome. The
**direction** of the crossover is an output.

| Histology | Chemo median survival | Nivo+ipi | HR 12mo | HR 24mo | CM743 HR |
|---|---:|---:|---:|---:|---:|
| Epithelioid | 13.7 mo | 9.8 mo | 1.98 | 0.78 | 0.86 |
| Biphasic | 8.8 mo | 8.1 mo | 0.82 | 0.43 | — |
| Sarcomatoid | 5.6 mo | 6.1 mo | 0.55 | **0.43** | **0.46** |

**CheckMate 743 is only reproduced once crossover in the control arm is put in.** About
44 % of the CM743 chemotherapy arm went on to receive subsequent systemic therapy. Mix
the "chemo → nivolumab" curve into the control arm in that proportion and:

```
 0 % crossover → HR 0.78        45 % crossover → HR 0.82   (trial 0.74)
sarcomatoid: no crossover 0.43 → 45 % crossover 0.46  (trial 0.46, an exact match)
```

**The epithelioid result sits on a threshold and is fragile.** Sweeping the maximum
immune kill rate over a 2.4-fold range:

| KKILL (/d) | Margin kill rate (/d) | Epithelioid HR 24mo |
|---:|---:|---:|
| 0.22 | 0.0296 | 1.11 |
| **0.33** | **0.0429** | **0.78** |
| 0.44 | 0.0552 | 0.52 |
| 0.53 | 0.0645 | 0.38 |

The quantity that matters is **the immune kill rate at the thin advancing margin
against that margin's own net proliferation rate (0.070 /d)**. Below it the front keeps
advancing and immunotherapy loses to six cycles of chemotherapy; close to it **the
front stops for as long as the drug is being given (two years, not five months)** and
immunotherapy wins on duration. The calibrated value is 0.043 /d — below the threshold,
on the steep part of the curve. **The epithelioid prediction of this model is therefore
badly conditioned**: a 30 % error in one rate constant moves the epithelioid HR from
1.1 to 0.5. This is stated rather than hidden.

### 3.2 Folate supplementation — toxicity comes down and efficacy does not

|  | Folate nM | Homocysteine μM | ANC nadir | Days ANC<0.5 | Delivered dose | Viable-cell nadir | Median survival |
|---|---:|---:|---:|---:|---:|---:|---:|
| Unsupplemented | 12.0 | 16.5 | **0.46** | **6 days** | 100 % | −77.6 % | 11.8 mo |
| Folate 400 µg + B12 | 43.4 | **7.6** | 0.61 | **0** | 100 % | −76.8 % | 12.0 mo |

The asymmetry is **the ratio of two rescue constants**. The folate term raises the
**marrow** EC50 by `1 + FOL/14 nM` and the **tumour** EC50 by only `1 + FOL/400 nM` —
**28.6-fold selectivity**. At 400 µg/d the marrow EC50 rises 4.21-fold and the tumour
EC50 1.11-fold. So grade 4 neutropenia disappears while tumour kill is unchanged.

### 3.3 Nephrotoxic positive feedback and a biomarker that lies

Cisplatin drives CrCl down, and because pemetrexed clearance is proportional to CrCl
the exposure goes up and the next nadir is deeper.

| Cycle | CrCl | pem AUC per cycle | ANC nadir | Viable cells | SMRP nM | SMRP with CrCl held fixed |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 95.0 | 9.20 | 0.61 | 259 | 4.02 | 4.02 |
| 3 | 88.4 | 9.94 | 1.25 | 175 | 2.93 | 2.74 |
| 6 | 82.2 | **10.75** | 1.21 | 82 | 1.48 | 1.32 |

Soluble mesothelin (SMRP) is a 40 kDa fragment filtered at the glomerulus. At day 120
viable tumour is at **−76.1 %** while SMRP has fallen only **−71.1 %**, and **+15.5 %
of that discrepancy is the kidney, not the tumour** (CrCl 95 → 79). The model's
quantitative recommendation is that **SMRP must be rescaled by a measured GFR before it
is read as a response biomarker**.

### 3.4 The antiangiogenic paradox

Blocking VEGF lowers the interstitial pressure (good for delivery) and prunes the
vessels (bad for delivery), and at the same time dries out the pleural cavity and slows
the front. The net effect is not monotonic in dose.

| Bevacizumab mg/kg | Free VEGF | Interstitial pressure mmHg | Vessel density | L_p mm | f_exp | Effusion mL | Median survival |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 128 | 9.9 | 0.72 | 0.781 | 0.184 | 482 | 12.0 |
| 2.5 | 26 | 5.5 | 0.30 | 0.617 | 0.148 | 224 | 12.8 |
| **7.5** | 11 | 4.7 | 0.18 | 0.513 | 0.122 | 185 | **13.0** |
| 15 | 5.9 | 4.4 | 0.14 | 0.461 | 0.109 | 172 | 13.0 |
| 40 | 2.3 | 4.1 | 0.11 | 0.416 | 0.098 | 163 | 13.0 |

**Penetration depth has already lost 34 % at 7.5 mg/kg** (0.781 → 0.513 mm) and the
survival gain plateaus there. The model predicts that half the MAPS dose is enough, and
the ground for that is not efficacy but **not losing any more delivery**.

### 3.5 Does the stroma a tumour makes become the barrier it meets — the author's hypothesis was rejected

Stroma disappears with a half-life of 87 days and the cycles are 21 days apart, so the
debris of cycle 1 is still in the diffusion path at cycle 4 and the stromal pool
**grows**, 144 → 196 cm³.

The author expected this to protect the later cycles. **It does not.** The exposed
fraction **rises** across the six cycles (0.095 → 0.165) and the kill per cycle
**accelerates** (−16.9 % → −26.5 %). This is because the rind thins faster than the
stroma accumulates. The stroma only **delays the improvement**. Isolating the size of
that with a `φ = 0` counterfactual (dead cells vanish instead of fibrosing):

```
with stroma (φ=0.70):  day 105 thickness 0.470 cm,  f_exp 0.165
no stroma (φ=0):       day 105 thickness 0.142 cm,  f_exp 0.470
→ the debris takes away 65 % of the exposed fraction that could have been recovered
```

So it is a real loss, but not the reversal that was expected. Where the stroma
**dominates** is not the kill but the **measurement** — at day 105, **62 %** of what the
scanner measures is stroma and not tumour.

### 3.7 The systemic and the intrapleural route cover different depths

The two routes advance from **opposite faces**, so their coverage does not overlap but
adds.

| Regimen | Viable-cell nadir | Maximum adhesion | Median survival |
|---|---:|---:|---:|
| Systemic chemotherapy alone | −76.8 % | 0.66 | 12.0 mo |
| Intrapleural alone | −82.1 % | 0.66 | 11.7 mo |
| Systemic + intrapleural | **−99.5 %** | 0.66 | **14.0 mo** |
| Talc pleurodesis → chemo | −76.9 % | 0.82 | 12.3 mo |
| Talc → chemo + intrapleural | −98.4 % | 0.82 | 13.9 mo |

Talc closes the very route it was meant to make usable — the effusion is better
controlled and yet "talc → both" is worse than "systemic + intrapleural".

**A caution about intrapleural therapy.** The model likes this route more than the
literature does, and the reason can be named: the instillate reaches the 0.6 mm
advancing margin from the free surface, so there `f = 0.24` while in the mass
`f = 0.021`. This depends almost linearly on an **assumed contact fraction of 0.55**
that nobody has measured. There is no trial comparing intrapleural with systemic
chemotherapy directly in MPM, so this is **a prediction, and a fragile prediction**.

### 3.6 Sequential and maintenance therapy (epithelioid)

| Regimen | Median survival | Days ANC<0.5 | Maximum irAE | Best mRECIST |
|---|---:|---:|---:|---:|
| Chemotherapy, 6 cycles | 12.0 mo | 0 | 0.00 | −59.1 % |
| Chemotherapy, 6 cycles + pemetrexed maintenance | **21.5 mo** | 0 | 0.00 | −74.2 % |
| Immunotherapy first, chemotherapy on progression | 9.1 mo | 0 | 0.65 | −32.1 % |
| Chemotherapy + immunotherapy added at day 140 | 12.5 mo | 0 | 0.65 | −59.2 % |
| Chemotherapy → second-line nivolumab | 12.1 mo | 0 | 0.35 | −59.1 % |
| Chemotherapy → second-line vinorelbine | 13.0 mo | 0 | 0.00 | −59.1 % |

The 21.5 months of maintenance pemetrexed is **a prediction of the model and it
conflicts with the literature** — the actual maintenance trials were negative or
marginal. The reason the gain is large in the model is clear: arrest of the front holds
only while the drug is being given, so keep giving it and it keeps stopping. This goes
on the list of untested predictions in §4.

---

## 4. What the model predicts but has not yet been measured

1. **mRECIST cannot see a 12-fold difference in drug access.** The author expected the
   partial-response rate to fall with baseline thickness. The model rejects that and
   says something sharper — hold the invaded area fixed and vary only the thickness, and
   **the measured response is flat while the actual kill is not**:

   | Baseline thickness | f_exp | Best mRECIST | Viable-cell nadir | Assessment |
   |---:|---:|---:|---:|:--:|
   | 0.25 cm | 0.333 | −63.1 % | **−96.2 %** | PR |
   | 0.50 cm | 0.160 | −62.2 % | −89.2 % | PR |
   | 0.80 cm | 0.095 | −59.1 % | −76.8 % | PR |
   | 1.20 cm | 0.061 | −59.0 % | −64.4 % | PR |
   | 1.80 cm | 0.039 | −61.6 % | −53.9 % | PR |
   | 2.50 cm | **0.028** | **−63.8 %** | **−47.5 %** | PR |

   The whole of a 10-fold range of thickness and a 12-fold range of exposed fraction is
   a partial response, and **mRECIST is flat at −59 ~ −64 % while the actual kill splits
   two-fold, from −96 % to −48 %**. This is because the measured thickness is dominated
   by loss of stroma and by spread of the sheet, and neither of those tracks the number
   of cells the drug reached. **The testable claim is the opposite of the expectation**:
   survival splits 15.3 against 4.4 months (§1(c)) while the depth of the mRECIST
   response should not stratify by baseline thickness.
2. **The fibrotic fraction of the measured thickness during treatment is not
   monotonic.** 28 % at diagnosis → 68.5 % at day 120 → 16 % at day 300 → 10 % at day
   500. Directly testable with an imaging–pathology correlation study.
3. **SMRP has to be rescaled by a measured GFR.** Over six cycles of cisplatin the size
   of the correction is +15.5 %.
4. **Intrapleural chemotherapy and talc pleurodesis are mutually exclusive because of
   geometry, not pharmacology.**
5. **Maintenance pemetrexed prolongs the arrest of the front** (model 21.5 vs 12.0
   months) — because this conflicts with the results of the existing maintenance trials,
   it is a prediction that can falsify the model.

---

## 5. Verification — and the five times the equations corrected the author

There is no R runtime in the environment in which this model was built. Rather than
simply hand over 51 ODEs that had never been run, **the whole system was written
twice** — once in mrgsolve C++, and once as the same equations in **dependency-free
Python RK4** in `mpm_reference_model.py`. Every number in this README and in the
repository table came out of running the Python version, and
`mpm_calibration_output.txt` is its verbatim output (114 individual simulations).

**Integration convergence.** With the two-stage step (dt = 0.02 d while the drug depots
are changing, 0.06 d thereafter), `V(400 d)` changes by **0.000 %** when dt is reduced
eight-fold.

**Five real defects that running the equations exposed** — each carries a `DEFECT #n`
comment where it was fixed:

| # | Defect | Symptom |
|---|---|---|
| 1 | No non-negativity floor | Immediately after an intrapleural instillation N went negative → thickness negative → `fpen()` negative → the sign of the rind flipped |
| 2 | No ceiling on the effusion | Formation rises with area while drainage falls, so it diverged to **1.3 million mL in a year**. The pleural cavity has a volume |
| 3 | The metastatic compartment grew without bound | It dominated the hazard in the second year → removing 89 % of the pleural tumour gave zero survival gain |
| 4 | A 10-minute infusion is shorter than the integration step | Whether an RK4 stage falls inside the infusion window decided the delivered dose. A **44 % difference** between dt 0.02 and 0.005, and the fitted EMAXP/EMAXC were compensating for the missing drug |
| 5 | Adhesion did not block the intrapleural route | Talc reduced the dilution volume and so came out **as if it potentiated intrapleural chemotherapy** — the opposite of what the map asserts |

**The simulation forced the author to withdraw two claims.**

- **§3.5** the hypothesis that the stroma would protect the later cycles — the exposed
  fraction rises instead and the kill per cycle accelerates. The stroma only delays
  (65 %), and where it dominates is the measurement.
- **§3.1** the early observation that the epithelioid HR was insensitive to the immune
  kill rate — that was an artefact of moving a compensating parameter along with it, and
  with the rest held fixed the HR moves a long way, 1.11 → 0.38. The real structure is
  **a threshold, not insensitivity**, and the calibration point sits on the steep part
  of that threshold.
- **§4.1** the expectation that the mRECIST partial-response rate would fall with
  baseline thickness — it came out flat (−59 ~ −64 %), and since that is the stronger
  claim the prediction was changed to it.

---

## 6. Where the model disagrees with the data (stated, not adjusted away)

1. **The histology × immunotherapy interaction is too flat.** CM743 split 0.86
   (epithelioid) / 0.46 (non-epithelioid). The model gives 0.83 / 0.46 after the
   crossover correction — sarcomatoid is exact but epithelioid still overestimates the
   benefit.
2. **The immunotherapy curves cross too late.** The model's HR is 1.52 at 12 months and
   0.78 at 24. The CM743 curves cross at about 8 months, so the model underestimates the
   early benefit and the median consequently falls on the wrong side of the crossing
   point.
3. **CONFIRM is the worst fit.** Second-line nivolumab, model HR 0.89 against an
   observed 0.69.
4. **Absolute median survival is low for the modern trial arms.** The cohort is anchored
   on EMPHACIS-era survival while MAPS, CM743 and MARS2 enrolled healthier patients. The
   only comparable thing is the hazard ratio, and the hazard ratio is what was fitted.
5. **The survival effect of supplementation is smaller than in the EMPHACIS subgroup.**
   The fully supplemented subset of the trial gave 13.3 against 10.0 months, whereas the
   model produces the difference in toxicity and only a small difference in survival
   (12.0 against 11.8). The EMPHACIS abstract itself says that supplementation "reduced
   toxicity without worsening survival", so the 13.3/10.0 contrast has era confounding
   mixed into it, and the model stands on the smaller side.
6. **The model likes intrapleural therapy more than the literature does.** The reason
   can be named — the instillate reaches the 0.6 mm advancing margin from the free
   surface. It is a fragile prediction that depends almost linearly on a **contact
   fraction (assumed 0.55)** that nobody has measured.

---

## 7. State variables (51 ODEs)

| Group | States |
|---|---|
| Tumour geometry (core) | `N` viable tumour · `M` necrotic/fibrous stroma · `A` invaded area · `Z` invasion depth · `MET` nodal/distant |
| Pleural cavity | `VEFF` effusion · `PSY` adhesion |
| Vasculature | `VEGF` · `RHOV` microvessel density |
| Inflammation · immunity | `IL6` · `TGFB` · `TEFF` · `TREG` · `TCLON` memory pool · `WOUND` surgical wound |
| Pemetrexed | `PEM_C` · `PEM_P` · `PEM_T` · `PEM_TP` tumour polyglutamates · `PEM_M` · `PEM_MP` marrow |
| Folate axis | `FOL` · `HCY` |
| Platinum | `CIS_F` free Pt · `CIS_B` bound Pt · `ADD` adducts · `ADDIP` intrapleural adducts |
| Biologicals | `BEV_C` · `BEV_P` · `NIV_C` · `NIV_P` · `IPI_C` |
| Arginine axis | `ADI_A` · `ADI_C` · `ADA` anti-drug antibody · `ARG` |
| Local | `IPD` drug in the pleural cavity |
| Marrow | `PROL` · `TR1` · `TR2` · `TR3` · `CIRC` (Friberg, MTT 110 h) |
| Kidney | `CRCL` · `CRCLSS` irreversible set-point |
| Other | `SMRP` · `LBM` · `CUMH` cumulative hazard · `AUCP` · `AUCC` · `VINE` · `IRAE` |

## 8. The 24 scenarios

Untreated · cisplatin alone · cis+pem (unsupplemented/supplemented) · carboplatin+pem ·
MAPS (+bevacizumab) · nivolumab+ipilimumab · pembrolizumab+chemo · chemo→second-line
nivolumab · ADI-PEG20+chemo · extended P/D+chemo (MARS2) · EPP+chemo+hemithoracic
radiotherapy · chemo+intrapleural cisplatin · talc pleurodesis→chemo ·
talc→chemo+intrapleural · chemo+PRMT5 inhibitor (MTAP-deleted) · pemetrexed maintenance ·
CrCl 52 mL/min · chemo→second-line vinorelbine · early (thin rind) · late (bulky) ·
intrapleural alone · immunotherapy first then chemotherapy · immunotherapy added after
chemotherapy

## 9. Running it

```bash
# map
dot -Tsvg mpm_qsp_model.dot -o mpm_qsp_model.svg
dot -Tpng -Gdpi=150 mpm_qsp_model.dot -o mpm_qsp_model.png

# reference implementation + full verification (~10 min, standard library only)
python3 mpm_reference_model.py > mpm_calibration_output.txt

# re-fetch the references (NCBI E-utilities)
python3 mpm_fetch_references.py

# R
Rscript -e 'source("mpm_mrgsolve_model.R"); print(run_all(), row.names = FALSE)'
shiny::runApp("mpm_shiny_app.R")
```

---

> ⚠️ This is a QSP model for educational and research purposes. It has not been
> independently validated and must not be used for clinical decision-making,
> prescribing, or regulatory submission. The parameters are approximations for the sake
> of illustration.
