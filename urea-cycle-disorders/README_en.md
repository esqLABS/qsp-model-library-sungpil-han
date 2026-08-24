# Urea Cycle Disorders (UCD) — QSP Model
### Ornithine transcarbamylase (OTC) deficiency as the index defect

<p align="center">
  <a href="ucd_qsp_model.svg">
    <img src="ucd_qsp_model.png" width="900" alt="UCD QSP mechanistic map">
  </a><br>
  <sub>Click to open a zoomable SVG · 183 nodes · 278 edges · 16 clusters</sub>
</p>

---

## 1. What This Model Sets Out to Explain

The urea cycle is the **only pathway by which the human body can discard
nitrogen in bulk.** Nitrogen from dietary protein and endogenous protein
breakdown reaches the periportal hepatocyte as ammonia, and is activated by
NAG through CPS1 → **OTC** → ASS1 → ASL → ARG1 to become urea, which is
excreted by the kidney. If any single step is blocked, the nitrogen comes
back.

The body's first response is not to let ammonia accumulate linearly but to
**buffer it as glutamine** (mainly via skeletal-muscle glutamine
synthetase, systemic flux ≈300 µmol/kg/h). As a result:

- **Plasma glutamine rises before ammonia does** — which is why it serves
  as an early-warning marker.
- **The disease is a cliff, not a slope** — the patient looks fine while
  the buffer holds, and once a nitrogen load exceeds residual
  urea-forming capacity for long enough, ammonia spikes abruptly.

The fatal organ is the brain, and the mechanism is **osmotic**. Astrocytes
are the only CNS cell carrying glutamine synthetase, so circulating
ammonia is trapped there as glutamine. Glutamine is an osmolyte, and
astrocytes defend their volume **by exporting myo-inositol.** But
**in chronic hyperammonaemia the myo-inositol reserve is already
depleted.** So a trigger that a normal brain would tolerate causes a
chronic patient to collapse catastrophically. This interaction is the
single most important non-intuitive behaviour in UCD, and this model
reproduces it explicitly as the `MINSB` compartment.

---

## 2. Deliverables

| File | Contents |
|------|------|
| [`ucd_qsp_model.dot`](ucd_qsp_model.dot) | Mechanistic map source (183 nodes · 278 edges · 16 clusters) |
| [`ucd_qsp_model.svg`](ucd_qsp_model.svg) | Vector rendering (zoomable, searchable) |
| [`ucd_qsp_model.png`](ucd_qsp_model.png) | Raster rendering (150 dpi) |
| [`ucd_mrgsolve_model.R`](ucd_mrgsolve_model.R) | 35-compartment ODE model + 14 scenarios + dose-conversion helpers |
| [`ucd_shiny_app_en.R`](ucd_shiny_app_en.R) | 8-tab interactive dashboard |
| [`ucd_references_en.md`](ucd_references_en.md) | 88 PubMed links + a parameter-to-literature cross-reference table |

---

## 3. Mechanistic Map — 16 Clusters

| # | Cluster | Key content |
|---|----------|-----------|
| 1 | Genetic substrate | NAGS · CPS1 · **OTC (Xp11.4)** · ASS1 · ASL · ARG1 · ORNT1 · citrin, allele grades, X-inactivation skewing in heterozygous females |
| 2 | Nitrogen load | Dietary protein, gut bacterial urease, intestinal glutaminase, muscle proteolysis, **catabolic triggers** (infection · fasting · surgery · postpartum · steroids · valproate · chemotherapy), anabolic rescue |
| 3 | Mitochondrial arm | NH4+ · HCO3- · AcCoA → **NAGS → NAG → CPS1** → carbamoyl phosphate → **OTC** → citrulline, ORNT1 transport |
| 4 | Cytosolic arm | ASS1 (aspartate = the 2nd nitrogen) → ASL → arginine → ARG1 → **urea**, fumarate re-entering the TCA cycle |
| 5 | Orotic acid overflow | Cytosolic CP efflux → CAD → DHODH → **orotic acid** — the diagnostic signal that distinguishes OTC from CPS1/NAGS |
| 6 | Glutamine trap | Periportal (high-capacity) / perivenous (high-affinity) compartmentalisation, muscle GS, renal glutaminase, deep reservoir |
| 7 | Ammonia kinetics | NH3/NH4+ partitioning (pKa 9.15), the **vicious cycle of respiratory alkalosis**, BBB diffusion, 100/200/360 thresholds |
| 8 | CNS pathophysiology | Astrocytic GS → brain glutamine → **myo-inositol depletion** → cell swelling → cerebral oedema → ICP → herniation; mPTP, NMDA excitotoxicity, EEG |
| 9 | Arginine | NO synthesis (ASL as a NOS-complex scaffold), creatine, growth; conversely accumulates in ARG1 deficiency |
| 10 | Phenylbutyrate axis | NaPBA / **GPB (lipase-dependent sustained release)** → PBA → β-oxidation → **PAA** → GLYATL1 + glutamine → **PAGN (2 nitrogens)** → OAT secretion; the PAA:PAGN ratio and neurotoxicity |
| 11 | Benzoate axis · activator | Benzoate + glycine → hippurate (**1 nitrogen**); **carglumic acid** (a NAG analogue, activates CPS1); L-citrulline (bypasses OTC), L-arginine; Na load · BCAA depletion |
| 12 | Extracorporeal removal · definitive therapy | Haemodialysis · CVVHDF · **post-dialysis rebound**, liver transplantation, AAV8-OTC (DTX301), OTC mRNA-LNP (ARCT-810), measurement of urea-forming flux |
| 13 | Renal handling · DDI | OAT1/3, probenecid, PAA accumulation in renal failure, hepatic failure, steroids |
| 14 | Hepatic · systemic involvement | Hepatomegaly, transaminase elevation, coagulopathy, trichorrhexis nodosa in ASL deficiency |
| 15 | Clinical endpoints | Crisis frequency, **coma duration**, neurocognition/IQ, growth, hospitalisation, survival, **natural protein tolerance** |
| 16 | Monitoring | Ammonia, glutamine, quantitative amino acids, orotic acid, urinary PAGN, PAA:PAGN, MRS, EEG |

---

## 4. mrgsolve Model (35 Compartments)

### 4.1 Compartment structure

- **PK (16)** — GPB intestinal pool · PBA gut/plasma · PAA · PAGN
  (plasma · cumulative urinary) · benzoate gut/plasma · hippurate
  (plasma · cumulative urinary) · citrulline gut/plasma · arginine
  gut/plasma · carglumic acid gut/plasma
- **Nitrogen balance (5)** — plasma ammonia · exchangeable glutamine ·
  **deep tissue glutamine reservoir** · cumulative urea nitrogen ·
  catabolic state
- **CNS (4)** — brain ammonia · brain glutamine (derived from circulating
  ammonia) · **brain myo-inositol** · excess brain water
- **Auxiliary pools (4)** — ornithine · glycine · BCAA · orotic acid
- **Gene therapy (2)** — vector template (AAV episome / LNP-mRNA) ·
  transgene-derived OTC activity
- **Integrators (4)** — AUC > 100 · cumulative time in the coma range ·
  neuronal injury index · excess PAA exposure

### 4.2 Self-calibrating baseline

Not a single initial value is hard-coded. In `$MAIN`, the pre-treatment
steady state is **solved numerically by bisection** from the patient
descriptors,

```
RA_N  =  A·C/(KM+C)  +  CL_NH4R·C  +  CL_RGLN·vgs(C)/(KGLNASE+CL_RGLN)
```

and from that root, glutamine, brain ammonia, brain glutamine,
myo-inositol, brain water, ornithine, orotic acid, glycine, and BCAA are
all derived in closed form and used as initial values. So **the untreated
arm is exactly at rest**, and every change seen in a simulation comes
purely from the mechanism.

Verification: healthy control (100% activity, protein 1.0 g/kg/d) →
ammonia 34.2 µmol/L, glutamine 601 µmol/L, brain glutamine 5.2 mmol/L,
orotic acid 0, protein tolerance 1.55 g/kg/d.

### 4.3 Nitrogen stoichiometry — the model's core

Glutamine carries **two** nitrogens. So glutamine synthesis subtracts 2 N
from the free-nitrogen pool, and glutaminase returns 2 N. Renal
glutaminase excretes 1 N as urinary NH4+ and returns 1 N. As a result:

- **1 molecule of PAGN = 2 nitrogens** permanently leave the body
- **1 molecule of hippurate = 1 nitrogen**

This is why phenylbutyrate is twice as nitrogen-efficient as benzoate on
a per-mole basis, and the model reproduces this structurally — not
through an arbitrary effect parameter.

---

## 5. 14 Predefined Scenarios and Results

Default patient: late-onset OTC male, 70 kg, residual activity 32%,
natural protein 0.50 g/kg/day (baseline ammonia 66 µmol/L, glutamine
879 µmol/L).

| # | Scenario | Duration | Final NH3 | Peak NH3 | Final Gln | Coma (h) | IQ | Protein tolerance |
|---|----------|------|---------|---------|---------|--------|-----|---------|
| 1 | Diet therapy alone | 90 d | 66.0 | 66.0 | 879 | 0 | 100 | 0.46 |
| 2 | Protein 0.45 + citrulline 170 mg/kg/d | 90 d | 54.8 | 57.9 | 778 | 0 | 100 | 0.46 |
| 3 | **NaPBA 20 g/d TID** + citrulline | 90 d | 41.5 | 66.0 | 597 | 0 | 100 | 0.59 |
| 4 | **GPB 17.5 mL/d TID** + citrulline | 90 d | 38.2 | 66.0 | 554 | 0 | 100 | 0.64 |
| 5 | Sodium benzoate 250 mg/kg/d | 90 d | 54.6 | 66.0 | 775 | 0 | 100 | 0.53 |
| 6 | GPB + benzoate combined | 90 d | 32.4 | 66.0 | 499 | 0 | 100 | 0.70 |
| 7 | Severe catabolic crisis (sepsis), no treatment escalation | 96 h | 265.3 | 268.2 | 2110 | 71.5 | 89.5 | 0.66 |
| 8 | Crisis + Ammonul + arginine + anabolic rescue (30 h) | 96 h | 25.6 | 95.8 | 440 | 0 | 100 | 0.63 |
| 9 | **Neonatal OTC (1%) — drug therapy only, no dialysis** | 168 h | 235.5 | 1181.7 | 1926 | 164.3 | 64.1 | 0.29 |
| 10 | **Same neonate + CVVHD (62 h)** | 168 h | 201.5 | 1181.7 | 1717 | 96.6 | 67.5 | 0.29 |
| 11 | NAGS deficiency + carglumic acid 150 mg/kg/d | 30 d | 30.5 | 117.2 | 557 | 0 | 100 | **1.35** |
| 12 | AAV8-OTC gene transfer (DTX301) + GPB | 90 d | 13.2 | 66.0 | 293 | 0 | 100 | **1.57** |
| 13 | GPB + renal and hepatic failure (PAA accumulation) | 21 d | 46.9 | 66.0 | 662 | 0 | 100 | 0.59 |
| 14 | Valproate initiation (OTC heterozygous female) | 10 d | 74.5 | 176.4 | 948 | 12.1 | 94.4 | 0.69 |

(NH3 and Gln in µmol/L, protein tolerance in g/kg/day. See
`UCD_summarise()` for the full column set.)

### Clinical observations the model reproduces

- **Scenarios 3 vs 4** — at the same PBA molar dose, GPB gives a lower
  PAA Cmax (85 → 61 µg/mL) and slightly better ammonia control. This
  matches the direction of the crossover-design result in Diaz 2013
  *Hepatology*, and inside the model the cause is a single **pancreatic
  lipase-dependent sustained-release step.**
- **Scenario 5** — benzoate accounts for only 15% of nitrogen disposal
  (PAGN accounts for 30%) — a direct consequence of the 1-nitrogen versus
  2-nitrogen stoichiometry.
- **Scenarios 9 vs 10** — in a neonatal crisis, dialysis drops ammonia
  from 1182 to 31 µmol/L within 4 hours, while **drug therapy alone can
  only bring it down to 236 over 108 hours.** And after dialysis ends,
  nitrogen redistributes from the deep glutamine reservoir, producing
  **rebound** (13 at 86 h → 181 µmol/L at 120 h). Even so, the final IQ
  gap between the two arms is only 3 points — the most sobering fact in
  this disease is that **it is not the treatment modality but the
  60-hour delay to diagnosis that determines prognosis.**
- **Scenario 11** — carglumic acid raises protein tolerance from 0.46 to
  1.35 g/kg/day in NAGS deficiency — the only UCD in which drug therapy
  cures the phenotype.
- **Scenario 13** — when renal and hepatic dysfunction overlap, PAA
  reaches 549 µg/mL and the PAA:PAGN ratio exceeds the threshold of 2.5
  (peak 19.3) — an early sign of conjugation saturation.

---

## 6. How to Run

```r
# 1) Build the model + run the scenarios
source("ucd_mrgsolve_model.R")
mod <- UCD_build()
sim <- UCD_simulate_scenarios(mod)          # all 14
UCD_summarise(sim)                          # summary table
UCD_plot(sim, "NH3", "Plasma ammonia (umol/L)", log = "y")

# 2) Individual patient simulation
mod |>
  mrgsolve::param(OTCACT = 0.15, PROT = 0.45, WT = 70) |>
  mrgsolve::mrgsim(events = ucd_gpb(17.5, n = 3, ii = 8, addl = ucd_addl(30, 8)),
                   end = 30*24, delta = 1, hmax = 0.25) |>
  plot(NH3 + GLNP + PAAUG ~ time)

# 3) Dashboard
shiny::runApp("ucd_shiny_app_en.R")
```

Required packages: `mrgsolve` (>= 1.0; developed and verified on 2.0.1);
the dashboard additionally needs `shiny`, `shinydashboard`, `dplyr`,
`tidyr`, `ggplot2`, `DT`.

### Dose-conversion helpers

```r
ucd_napba(20)                      # NaPBA 20 g/day TID
ucd_gpb(17.5)                      # GPB 17.5 mL/day TID (≈ same PBA molar dose as above)
ucd_benzoate(250, wt = 70)         # sodium benzoate 250 mg/kg/day
ucd_citrulline(170, wt = 70)       # L-citrulline 170 mg/kg/day
ucd_carglumic(150, wt = 70)        # carglumic acid 150 mg/kg/day
ucd_ammonul(wt = 70)               # Ammonul loading + 24-hour maintenance infusion
ucd_arginine_iv(400, wt = 70)      # L-arginine HCl 400 mg/kg IV
ucd_aav(200); ucd_mrna(900)        # AAV8-OTC / OTC mRNA-LNP
```

### Key parameters

| Parameter | Meaning | Default |
|----------|------|--------|
| `OTCACT` · `CPS1A` · `NAGSA` · `DISTALA` | Residual activity fraction at each locus | 0.32 / 1 / 1 / 1 |
| `PROT` · `PROT2`/`TPROT` · `PROT3`/`TPROT3` | Natural protein intake and its two change points | 0.50 g/kg/d |
| `CATAMP` · `CATT0` · `CATDUR` · `ANAB` | Catabolic trigger and anabolic rescue | 0 |
| `CLHD` · `HDT0` · `HDDUR` | Extracorporeal ammonia clearance rate and time window | 0 |
| `VPAI` · `TVPA` | Valproate's NAGS-inhibition fraction and timing | 0 |
| `ACTX` · `TTX` | Total activity after liver transplant and the transplant time | none |
| `WT` · `VBR` · `ALLOM` | Body weight · brain volume · allometric exponent (scalable down to neonates) | 70 / 1.35 / 0.75 |

---

## 7. Known Limitations

- The amino acid network is not solved explicitly and is instead reduced
  to a **single free-nitrogen pool.** Both of glutamine's nitrogens are
  subtracted from and credited back to that pool, which closes the
  balance consistently at the level of whole-body nitrogen but does not
  predict individual amino acid kinetics.
- `IQEST` is an **illustrative index tuned to reproduce the
  coma-duration–prognosis relationship**; it does not predict an
  individual's prognosis.
- Because placental clearance is not represented, the "normal ammonia"
  window immediately after birth is not captured. Scenarios 9/10 begin at
  12 hours of age.
- ARG1 deficiency (the arginine-accumulation type) and HHH syndrome are
  present on the map but are handled in the ODEs only through the
  `DISTALA` reduction.

---

> ⚠️ **Disclaimer**
> This model is a **semi-quantitative QSP model for educational and
> research purposes.** It was built on published literature but has not
> been independently validated or certified. **It must not be used
> directly for clinical decision-making, prescribing, or regulatory
> submission.**
