# Alcohol Use Disorder (AUD) — QSP Model

**Quantitative systems pharmacology model**

<a href="../../../alcohol-use-disorder/aud_qsp_model.svg"><img src="../../../alcohol-use-disorder/aud_qsp_model.png" width="100%" alt="AUD QSP mechanistic map"></a>

*(click the map for the zoomable SVG)*

---

## 1. The premise

In every other model in this library **the dose is an input**. In AUD it is different.

> **The dose is an output.** The patient sets it, and what makes that decision is the
> pharmacology of the alcohol already drunk.

So this model was designed as a **closed loop** rather than a cascade. `drinks/day` is
not a covariate but a model variable, generated every evening from the state of the
reward, stress and control systems.

```
DRIVE ──► drinks/day ──► ethanol PK ──► MOR/dopamine reinforcement
  ▲                                        │
  │                                        ▼
  └── craving ◄── negative affect ◄── allostatic neuroadaptation
                                           │
             withdrawal revealed as BAC falls ──┘
```

The three consequences that follow from building it this way are the whole of this model.

1. **The disease is generated rather than assumed.** A healthy drinker and a severe AUD
   patient differ in **only two parameters**, environment (`ENVDRIVE`) and vulnerability
   (`VULN`). The escalation trajectory is emergent.
2. **Drugs do not "lower an endpoint".** Perturb one arm of the loop and the loop finds a
   new operating point. This is where the reason naltrexone and acamprosate move
   **different endpoints** comes from.
3. **Withdrawal is not a separate model.** It is what happens to *the same equations*
   when the ethanol input stops while the adaptation states are still high.

### Two pharmacological arms, explicitly separated

| | (A) The positive reinforcement arm | (B) The negative reinforcement arm |
|---|---|---|
| Mechanism | ethanol → β-endorphin → MOR → disinhibition of VTA GABA interneurons → phasic dopamine in the nucleus accumbens → **within-episode escalation (priming)** | chronic ethanol → GABA-A subunit reconfiguration + NR2B↑ + GLT-1↓ + dynorphin/CRF↑ → negative affect · insomnia · protracted withdrawal → **relief craving** |
| Drugs that cut it | naltrexone · 6-β-naltrexol · nalmefene · semaglutide · ondansetron · baclofen · topiramate (AMPA) | acamprosate · gabapentin · topiramate (GABA) · benzodiazepines (detoxification only) |
| Clinical signature | **DPDD (drinks per drinking day) · %HDD** | **%PDA (percentage of days abstinent) · time to first drink** |

A third arm is the upstream **aversive** arm — disulfiram → ALDH2 inactivation →
acetaldehyde → flushing. This arm is opened not by exposure but by **supervision**.

---

## 2. Deliverables

| File | Contents |
|------|------|
| [`aud_qsp_model.dot`](../../../alcohol-use-disorder/aud_qsp_model.dot) | Graphviz mechanistic map source — **21 modules · 248 nodes · 403 edges** |
| [`aud_qsp_model.svg`](../../../alcohol-use-disorder/aud_qsp_model.svg) / [`.png`](../../../alcohol-use-disorder/aud_qsp_model.png) | The rendered map (the PNG at 150 dpi) |
| [`aud_mrgsolve_model.R`](../../../alcohol-use-disorder/aud_mrgsolve_model.R) | mrgsolve model with **72 ODE compartments · 271 parameters · 21 scenarios** + validation report (actually run) |
| [`aud_shiny_app.R`](../../../alcohol-use-disorder/aud_shiny_app.R) | **Nine-tab** interactive dashboard (including a 13-regimen builder) |
| [`aud_references.md`](aud_references.md) | **132 papers** (21 sections) — every PMID looked up through the PubMed E-utilities and cross-checked against the title. The 15 that failed the check were excluded rather than guessed at, and are noted at the end of the list |

Reproduction:

```bash
dot -Tsvg aud_qsp_model.dot -o aud_qsp_model.svg
dot -Tpng -Gdpi=150 aud_qsp_model.dot -o aud_qsp_model.png
Rscript aud_mrgsolve_model.R          # runs the 21 scenarios + prints the validation report
Rscript -e 'shiny::runApp("aud_shiny_app.R")'
```

---

## 3. Model structure

### 72 ODE compartments

| Group | Number | States |
|---|---|---|
| Ethanol ADME · metabolism | 10 | `ETH_ST` `ETH_GUT` `ETH_C` `ACD` `ACT` `CYP2E1` `NADH` `LACT` `ETHAVG` `ETOHCUM` |
| Objective biomarkers | 4 | `PETH` `CDT` `GGT` `MCV` |
| Fast neurotransmission | 7 | `GABA_TONE` `GLU_NAC` `DA_NAC` `BEND` `NE_LC` `CRAVE` `DRIVE` |
| Slow neuroadaptation | 10 | `GABAA_SUB` `NMDA_UP` `GLT1` `DYN` `CRF_CEA` `NPY_CEA` `ALLO` `CUE` `HABIT` `PFC` |
| HPA axis | 4 | `CRH_HYP` `ACTH` `CORT` `GR_SENS` |
| Affect · sleep · withdrawal | 5 | `NEGAFF` `SLEEPD` `CIWA` `KINDLE` `AVEXP` |
| Liver · nutrition · cardiovascular | 9 | `STEAT` `LPS` `TNFA` `ALT` `AST` `FIB` `THIA` `SBP` `CMYO` |
| Drug PK | 23 | naltrexone 4 (+depot) · acamprosate 2 · disulfiram 3 (+`ALDH2A`) · topiramate 2 · gabapentin 2 · baclofen 2 · nalmefene 2 · ondansetron 2 · semaglutide 2 · benzodiazepine 2 |

### Three structural choices

**(1) Chronic exposure drives every slow adaptation.**
GABA-A subunit reconfiguration, NR2B upregulation, GLT-1 downregulation, dynorphin/CRF
induction and PFC erosion are driven not by the *instantaneous* ethanol concentration but
by `CHRON`, a sigmoid function of the 7-day moving average (`ETHAVG`). This is **the
point that separates someone drinking two drinks a day from someone escalating**. Drive
it from the instantaneous concentration and the two adaptation trajectories become
identical.

**(2) Incentive sensitisation has a threshold.**
The acquisition rate of `CUE` (cue reactivity) is proportional to **the square of the
excess above threshold** of the phasic dopamine response (`POSPART(REINF − 0.35)²`). This
is consistent with the observation that what drives incentive sensitisation is not total
exposure but *intermittent high-concentration peaks*, and it creates **two attractors** in
the model — stable controlled drinking and AUD.

**(3) The drinking window is a smooth Gaussian.**
Instead of a hard switch, a Gaussian intake rate centred at 20:00 with σ = 0.9 h is used.
The solution does not jump across a discontinuity, and as a result the drinks → peak BAC
relation matches the Widmark prediction (5 drinks → 0.075 g/dL, 3 drinks → 0.035 g/dL).

### Benzodiazepines integrate CIWA

The symptom-triggered protocol was implemented as a **continuous analogue** rather than as
events: `rate = BZDGAIN × max(0, CIWA − 8)`, capped at 7 days. So the total diazepam
equivalent dose emerges as **something the model computes** rather than as a prescription.

---

## 4. Scenarios (21)

| # | Scenario |
|---|---|
| S01 | Stable moderate drinker (no vulnerability, low-cue environment) |
| S02 | Untreated severe AUD (natural history) |
| S03 | Medical management / CBT alone (the placebo arm of COMBINE) |
| S04 | Acute cessation, no withdrawal treatment |
| S05 | Symptom-triggered benzodiazepine detoxification (CIWA-Ar ≥ 8) |
| S06 | Naltrexone 50 mg PO daily + MM |
| S07 | Extended-release naltrexone 380 mg IM q4wk + MM |
| S08 | Acamprosate 666 mg tid + MM |
| S09 | Naltrexone + acamprosate (the COMBINE combination cell) |
| S10 | Disulfiram 250 mg, **supervised** |
| S11 | Disulfiram 250 mg, **unsupervised** (~25% adherence) |
| S12 | Topiramate titrated to 300 mg/day |
| S13 | Gabapentin 600 mg tid (1800 mg/day) |
| S14 | Nalmefene 18 mg as needed (reduction goal) |
| S14R | Reduction-goal control (no drug, `ABSTG` 0.20) — nalmefene's matched comparator |
| S15 | Baclofen 60 mg tid (the high-dose protocol) |
| S16 | Semaglutide 0.25 → 1.0 mg SC weekly |
| S17 | Ondansetron 4 µg/kg bid, the 5-HTTLPR LL early-onset subtype |
| S18 | Naltrexone in an OPRM1 118G carrier (pharmacogenomics) |
| S19 | ALDH2*2 heterozygote (endogenous disulfiram), same environment · no drug |
| S20 | A kindled patient (5 prior detoxifications), abstinence with no drug |

---

## 5. Validation

The full report is reproduced by `Rscript aud_mrgsolve_model.R`. What follows is its
output.

### The baseline is a fixed point, not an assumption

In a 180-day simulation with no drug and no ethanol, the **maximum drift across the 24
states is 0.0000 %**. Because every homeostatic state was written in the form
`dxdt = k(X_base·(1+stimulus) − X)`, with stimulus at zero it is an algebraically exact
fixed point.

### The ethanol PK reproduces Widmark

A 75 kg man, single oral load of 0.8 g/kg:

| Item | Model | Literature expectation |
|---|---|---|
| Cmax BAC | **0.0848 g/dL** | 0.075–0.105 |
| Tmax | **1.30 h** | 0.5–1.2 |
| β-slope (descending limb) | **0.0137 g/dL/h** | 0.012–0.020 |
| Peak acetaldehyde (wild type) | **1.40 µM** | 1–3 |
| Peak acetaldehyde (ALDH2\*2 het) | **14.0 µM** | 5–20× wild type |

### The disease is generated (3 years, a cue-rich environment, a vulnerable host)

| Day | drinks/day | g/day | GABAAsub | NMDAup | GLT1 | CUE | HABIT | PFC | GGT | PEth |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 4.7 | 65 | 0.35 | 0.26 | 0.94 | 0.64 | 0.00 | 0.98 | 39 | 71 |
| 180 | 6.4 | 89 | 0.69 | 0.61 | 0.74 | 2.28 | 0.03 | 0.71 | 68 | 148 |
| 365 | 6.7 | 94 | 0.72 | 0.65 | 0.70 | 2.75 | 0.14 | 0.62 | 73 | 166 |
| 1095 | 6.8 | 96 | 0.73 | 0.66 | 0.69 | 2.90 | 0.25 | 0.60 | 74 | 172 |

The severe AUD state at 3 years: **6.8 drinks/day (96 g/day)**, PEth 172 ng/mL, GGT
74 U/L, CDT 1.86 %, MCV 96.2 fL, **AST/ALT 1.93** (105/54 U/L), steatosis 22 %,
SBP 130 mmHg, thiamine 0.59. All within the ranges observed at this level of consumption.

### The handover test — the loop is self-consistent

Cutting off the forced drinking history and handing over to the endogenous loop, then 90
days: **6.88 → 7.27 drinks/day (+5.6 %)**. That is, the loop sustains by itself the
consumption it generated (which means the parameters mesh with one another; a large
departure from this would mean the loop gain was wrong).

### Withdrawal is graded

| Scenario | Peak CIWA-Ar | Time of peak | Seizure risk | Diazepam equivalent |
|---|---|---|---|---|
| S04 untreated cessation | 11.5 | **42 h** | 2.1 % | — |
| S05 symptom-triggered benzodiazepine | 11.1 | — | 1.7 % | **19 mg** |
| S20 kindled (5 prior detoxifications) | **22.5** | **68 h** | **11.7 %** | — |

The peak arriving at 24–48 h, being later and more severe in the kindled patient, and the
seizure risk diverging more than 5-fold, all agree with the observed natural history.
Kindling is set up to accumulate *between* episodes (threshold CIWA ≥ 15, t½ ≈ 6.6 years)
so it does not saturate after a single detoxification.

### Drug exposures reproduce the literature values

| Agent | Model (steady state) | Reported |
|---|---|---|
| Naltrexone | 1.09 ng/mL (Cavg) | 1–10 (Cmax 5–10) |
| 6-β-naltrexol | **73.1 ng/mL** | 40–100 |
| MOR occupancy (oral 50 mg) | **97.9 %** | ≥ 90 % (PET) |
| 6-β-naltrexol on extended-release naltrexone | 3.9 ng/mL | 2–8 sustained |
| MOR occupancy on extended-release naltrexone | 94.3 % | ≥ 80 % throughout the month |
| Acamprosate | 300 ng/mL | 350–500 (2 g/day) |
| Disulfiram + metabolites | 0.55 mg/L | 0.5–2 |
| Topiramate | 6.17 mg/L | 5–8 (300 mg/day) |
| Gabapentin | 3.91 mg/L | 3–6 (1800 mg/day) |
| Baclofen | 512 ng/mL | 200–600 (high dose) |
| Semaglutide | 23.3 nM | 20–30 (1.0 mg weekly) |

That oral naltrexone produces 70 times more 6-β-naltrexol while the intramuscular
extended-release form does not is because first-pass conversion (`FMNOL` 0.55) and
systemic conversion (`FSYSNOL` 0.12) were kept **separate**. Without that distinction the
metabolite exposures of the two formulations cannot be matched simultaneously.

### Six-month clinical endpoints (the 5–26 week window)

| Scenario | g/day | DPDD | %HDD | %PDA | Craving |
|---|---|---|---|---|---|
| S01 stable moderate drinker | 25.8 | 2.79 | 15.9 | 33.9 | 0.05 |
| S02 untreated severe AUD | 95.7 | 6.95 | 63.3 | 1.6 | 5.87 |
| **S03 MM alone (the reference arm)** | **47.5** | **5.52** | **33.3** | **38.5** | **1.49** |
| S06 naltrexone PO | 37.2 | 4.98 | 26.4 | 46.6 | 0.94 |
| S07 extended-release naltrexone | 37.5 | 5.01 | 26.5 | 46.6 | 0.94 |
| S08 acamprosate | 43.4 | 5.91 | 30.0 | 47.6 | 0.95 |
| S09 naltrexone + acamprosate | 35.7 | 5.18 | 25.2 | 50.7 | 0.72 |
| S10 disulfiram (supervised) | 15.7 | 2.34 | 8.6 | 52.2 | 0.65 |
| S11 disulfiram (unsupervised) | 32.1 | 4.48 | 22.8 | 48.7 | 0.84 |
| S12 topiramate | 36.3 | 5.11 | 25.7 | 49.2 | 0.82 |
| S13 gabapentin | 43.2 | 5.99 | 29.8 | 48.5 | 0.90 |
| S14 nalmefene (as needed) | 45.7 | 4.60 | 32.4 | 29.1 | 1.19 |
| S14R reduction-goal control | 67.2 | 5.47 | 47.2 | 12.3 | 2.95 |
| S15 baclofen | 41.3 | 5.25 | 29.2 | 43.7 | 1.11 |
| S16 semaglutide | 39.7 | 5.09 | 28.1 | 44.2 | 1.03 |
| S17 ondansetron (LL) | 39.4 | 5.08 | 27.9 | 44.6 | 1.06 |
| S18 naltrexone (OPRM1 118G) | 33.4 | 4.63 | 23.7 | 48.6 | 0.84 |
| S19 ALDH2\*2 heterozygote | 33.6 | 3.50 | 22.9 | 31.3 | 0.21 |

### And this is the central result of the model — the dissociation

Change against MM alone (S03). The observed values are from the Jonas 2014 JAMA
meta-analysis and the registration trials.

| Agent | Δ%HDD (model) | Δ%HDD (literature) | ΔDPDD (model) | ΔDPDD (literature) | Δ%PDA (model) | Δ%PDA (literature) |
|---|---|---|---|---|---|---|
| Naltrexone 50 mg | **−7.0** | −4 ~ −6 | **−0.54** | −0.5 | **+8.2** | +4 |
| Extended-release naltrexone | −6.8 | (heavy-drinking events −25 %) | −0.51 | | +8.1 | |
| **Acamprosate** | **−3.3** | ~ 0 | **+0.40** | ~ 0 | **+9.1** | **+9** |
| Combination | −8.1 | | −0.34 | | +12.3 | |
| Topiramate 300 mg | −7.6 | −8 ~ −9 | −0.41 | −1.0 | +10.8 | +10 |
| Gabapentin 1800 mg | −3.5 | −9 | +0.47 | −1.0 | +10.0 | +13 |
| Nalmefene (vs S14R) | −14.8 | −7.7 (−2.3 d/month) | −0.87 | | +16.8 | |
| Disulfiram supervised | −24.7 | (open-label supervised SMD ~0.6) | −3.18 | | +13.7 | |
| Disulfiram unsupervised | −10.6 | ~ 0 (Fuller 1986) | −1.04 | | +10.3 | |
| Baclofen | −4.2 | mixed / small | −0.27 | | +5.3 | |
| Semaglutide | −5.2 | reduction in DPDD (phase 2) | −0.43 | | +5.7 | |
| Ondansetron (LL) | −5.4 | subtype-specific | −0.44 | | +6.1 | |
| Naltrexone (OPRM1 118G) | −9.6 | larger in G carriers | −0.88 | | +10.1 | |

**Naltrexone lowers DPDD by −0.54 and acamprosate raises it by +0.40.** And it does so
while both raise %PDA in the same direction. This inversion of sign is the dissociation
repeatedly observed in clinical trials, and in the model it emerges automatically because
the two drugs cut *different arms* of the loop — it was not tuned endpoint by endpoint.

### The statistical layer from deterministic trajectory to trial endpoint (made explicit)

The ODEs give only a single mean daily consumption. A real cohort does not drink the mean
every day; it mixes abstinent days and drinking days. So this distribution is **exposed
rather than hidden**:

```
probability that day i is abstinent   π = plogis(Z0 + Z1·ABSTG + Z2·ΔNEGAFF + Z3·ΔCRAVE)
on all other days                     drinks ~ zero-truncated NegBin(mu = DPDD, size = 2.4)
                                      DPDD  = μ_model / (1 − π)
```

`ΔNEGAFF` and `ΔCRAVE` are *fractional* reductions against a matched no-drug arm.
Drugs that cut the negative reinforcement arm (acamprosate · gabapentin) raise π, and
drugs that cut the positive reinforcement arm (naltrexone · nalmefene) lower μ.
The mapping contains not one drug-specific term.

---

## 6. What the model says

- **"How good is this drug" is the wrong question.** The right question is *which arm it
  cuts*. A patient trying to reduce heavy drinking days and a patient trying to maintain
  abstinence do not need the same drug.
- **Compare nalmefene with an abstinence-goal arm and it looks bad** (Δ%HDD +5.3).
  Compare it with a matched reduction-goal control (S14R) and it is −14.8. Measuring a drug
  whose label says "reduction" against an abstinence endpoint is a problem with the design,
  not with the drug — in the model this difference is the single parameter `ABSTG`.
- **Disulfiram's effect size comes from supervision, not from pharmacology.** The degree of
  ALDH2 inhibition in the supervised and unsupervised arms differs a bit more than 2-fold
  (0.16 against 0.44), whereas the aversive expectation (`AVEXP`) differs 5-fold, 2.19
  against 0.46. The reason people actually do not drink is not that they have experienced
  the flushing but that **they know they would**, and that knowledge is maintained only
  where a dose cannot be skipped.
- **The ALDH2\*2 heterozygote (S19) shows a reduction in consumption comparable to
  acamprosate with no drug at all** (33.6 g/day vs S03's 47.5). Same environment, same
  vulnerability, same loop — the only difference is acetaldehyde. The protective effect
  observed in East Asia is reproduced here. That the same acetaldehyde raises the risk of
  upper gastrointestinal cancer is also drawn on the map.
- **Detoxifying a patient who has been detoxified many times is not the same event.** At
  the same consumption the CIWA peak goes from 11.5 to 22.5 and the seizure risk from
  2.1 % to 11.7 %. The kindling index has a t½ ≈ 6.6 years and effectively does not come
  back.

---

## 7. Limitations

- **It is deterministic and semi-quantitative.** The parameters are order-of-magnitude
  estimates from the literature, not a fitted population model. The arm-level %HDD/%PDA
  come from the explicit statistical layer in §5 above, and were not obtained by counting
  drinking days in a single trajectory.
- **The escalation rate is compressed.** Real AUD develops over years to decades, whereas
  the model reaches equilibrium in 6–12 months. It can be stretched by lowering
  `KCUEON`/`KHABON`, but the time constants of the GABA-A and NMDA adaptations (days to
  weeks) are physiologically fixed, so slowing the whole thing to the real rate only
  increases the simulation cost.
- **Gabapentin's %HDD effect is underpredicted** (−3.5 against the literature's −9). In the
  model gabapentin acts only on the relief arm, so %PDA (+10.0, literature +13) matches but
  it does not lower heavy-drinking-day intensity enough. This may mean part of the effect in
  Mason 2014 is not explained by the relief arm, or it may be the consequence of a
  modelling choice not to give the α2δ action a consumption-suppressing component.
- **The effect sizes of nalmefene and supervised disulfiram are overpredicted** (about
  2-fold and 1.5-fold the literature respectively). Both drugs press consumption down
  directly through the `WPRIME` and `AVEXP` routes, because the loss of adherence and the
  dropout of real trials were not put into the model.
- **The model has two attractors.** The transition between controlled drinking and AUD is
  steep. The region in between (hazardous but stable drinking) is narrow in parameter
  space, and in reality many people are in that region.
- **Behaviour is a single scalar (`DRIVE`).** Situational specificity (does not drink at
  home but drinks in a bar), social reinforcement, price and availability, and treatment
  retention are all folded into the two parameters `ENVDRIVE` and `ABSTG`.
- **Hepatic fibrosis (`FIB`) barely moves in a six-month simulation.** Because its
  accumulation time constant is in years, it is decorative within this model's window.
- **`SEIZP` (seizure probability) and `%HDD` are post-hoc mappings, not ODE states.** Both
  are logistic functions and the coefficients are exposed in the code, so change them if
  you disagree.

---

## 8. Disclaimer

This is a QSP model for educational and research purposes. It was constructed from the
public literature and clinical trial data but has not been independently validated or
certified. **It must not be used directly for clinical decision-making, prescribing, or
regulatory submission.** In particular, alcohol withdrawal is a condition that can be
fatal if untreated, and this model's CIWA and seizure-risk outputs are not a substitute
for clinical assessment of any kind.

Full citation list: [`aud_references.md`](aud_references.md) (132 papers, PubMed-verified).
