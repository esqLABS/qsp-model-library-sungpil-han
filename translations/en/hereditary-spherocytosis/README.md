# Hereditary spherocytosis QSP model

> **A red cell is two numbers.** Membrane area $S$ and volume $V$. That is the whole
> premise of this model. The two numbers set the minimum cylindrical diameter $D_c$, and
> $D_c$ is the only thing the wall of a splenic sinusoid can measure about a red cell.
> Unlike the other models in this repository, this one **calibrates only the normal red
> cell and predicts the disease.**

| Deliverable | File |
|--------|------|
| 🗺️ Mechanistic map (148 nodes · 18 clusters · 224 edges) | [`hsph_qsp_model.dot`](hsph_qsp_model.dot) · [SVG](hsph_qsp_model.svg) · [PNG](hsph_qsp_model.png) |
| 📐 Geometry kernel (closed-form, self-verifying) | [`hsph_geometry.py`](../../../hereditary-spherocytosis/hsph_geometry.py) |
| ⚙️ mrgsolve model (72 ODEs) | [`hsph_mrgsolve_model.R`](hsph_mrgsolve_model.R) |
| 🐍 The reference implementation that was actually run (the model of record) | [`hsph_python_reference.py`](../../../hereditary-spherocytosis/hsph_python_reference.py) |
| 🎯 Calibration ledger | [`hsph_calibrate.py`](../../../hereditary-spherocytosis/hsph_calibrate.py) · [`hsph_calibration.json`](../../../hereditary-spherocytosis/hsph_calibration.json) · [`calib.log`](../../../hereditary-spherocytosis/calib.log) |
| 🔬 Full run output of all 18 analysis sections | [`hsph_analysis.py`](../../../hereditary-spherocytosis/hsph_analysis.py) → [`hsph_reference_output.txt`](../../../hereditary-spherocytosis/hsph_reference_output.txt) |
| 📊 Shiny dashboard (14 tabs) | [`hsph_shiny_app.R`](hsph_shiny_app.R) |
| 📚 References (115 papers, live PubMed lookup) | [`hsph_references.md`](hsph_references.md) · [`mkrefs.py`](mkrefs.py) |

---

## 1. Why this disease is a geometry problem

The red cell membrane is **area-incompressible** (stretch it more than 3–4% and it
ruptures). So when a cell with area $S$ and volume $V$ is deformed into a cylinder capped
with hemispheres,

$$S = \pi D L,\qquad V = \frac{\pi D^2}{4}(L-D) + \frac{\pi D^3}{6}
\;\Longrightarrow\; V = \frac{SD}{4} - \frac{\pi D^3}{12}$$

The **smallest positive root** of this cubic is the diameter of the narrowest cylinder
that cell can pass through, i.e. the **minimum cylindrical diameter $D_c$**
(Canham & Burton 1968). Using the trigonometric solution it takes a closed form:

$$V_{\rm sph} = \frac{S^{3/2}}{6\sqrt{\pi}},\qquad s=\frac{V}{V_{\rm sph}},\qquad
\boxed{\,D_c = 2\sqrt{\frac{S}{\pi}}\,\cos\!\left(\frac{\arccos(-s)}{3}-\frac{2\pi}{3}\right)}$$

$V_{\rm sph}$ is the **maximum volume** that membrane can contain — that is, the critical
haemolytic volume. From this, **with no fitting**, the following follow.

| What is derived | Model | Literature |
|---|---|---|
| Normal red cell $D_c$ | 2.72 µm ($S{=}140,V{=}90$) — a range of 2.72–3.06 µm | ~2.8 µm |
| Normal 50% osmotic haemolysis | 0.418 %NaCl | 0.40–0.45 %NaCl |
| HS 50% osmotic haemolysis | 0.537 %NaCl ($S{=}118,V{=}86$) | 0.55–0.65 %NaCl |
| Reduction in EMA binding | the membrane area ratio itself | diagnostic criterion 16–21% |

**A point not in the textbooks.** Waugh (1992) reported that a 50-day-old red cell loses
10.5% of its area and 8.4% of its volume, and that "sphericity barely changes". Run the
same area loss with and without the accompanying volume loss and the change in $s$ is
$+0.047$ against $+0.105$ — **the same membrane loss makes the cell 2.21 times more
spherical when the volume does not go with it.** The reason HS vesicles cause disease is
not their number but **their cargo** (HS vesicles contain no haemoglobin).

---

## 2. The calibration ledger — seven numbers, and everything else is prediction

| Stage | Numbers consumed | Parameters obtained |
|---|---|---|
| **1. The normal red cell only** (nothing to do with HS) | mean lifespan 120 days · area −10.5% · volume −8.4% (Waugh) · geometric splenic clearance is 1% of the total in normals | `tau50` `kv_base` `kd_base` `k_ph` |
| **2. The reference moderate patient** | Hb 10.5 g/dL, **one number** | `fdef_mod` |
| **3. The opsonin arm** (Reliene 2002) | 140 IgG molecules/cell in a splenectomised band 3 patient · opsonin clearance is 50% of destruction pre-splenectomy | `K_cl` `k_ops` |

**Fixed structurally rather than fitted** (with the grounds stated):
`kv_def = 3.33·kv_base` · `kd_def = 6.67·kd_base` · `kd_cord = 41·kd_base` ·
`a_ent_def = 1.0` · `visc_k = ln2/3 = 0.231` (Chien's cytoplasmic viscosity–MCHC slope) ·
intravascular haemolysis restricted to a minor route (HS haemolysis is more than 95%
extravascular).

The value obtained at stage 2 is `fdef_mod = 0.420`, i.e. **58% spectrin content**. The
literature value for moderate HS is about 70%, so this model requires **a defect 12
percentage points larger than the literature** to produce textbook moderate disease. This
deviation was applied **identically** across the mild and severe ranges too, so the
remaining three genotypes are predictions rather than fits.

### What follows from that one number (all of it prediction)

| Measure | Model | Literature / criterion |
|---|---|---|
| Reticulocytes | 8.75 % | moderate >6% ✓ |
| Total bilirubin | 3.30 mg/dL | moderate >2 ✓ |
| Red cell lifespan | 17.5 days | HS 10–30 days ✓ |
| MCV | 85.1 fL | ✓ |
| MCHC | 35.26 g/dL | screening criterion ≥35.4 (almost matched) |
| 50% osmotic haemolysis | 0.571 %NaCl | HS 0.55–0.65 ✓ |
| EMA reduction | 20.9 % (against the normal population mean) | criterion 16–21 ✓ |
| Spleen volume | 385 mL | HS 300–800 ✓ |
| RDW | 11.7 % | screening criterion ≥14 ✗ (see the limitations below) |

---

## 3. Three results, and two failures

### (1) Normal red cells and HS red cells die in the same organ for different reasons

Decomposition of the clearance routes (run output §4):

| | Geometric (splenic cords) | Opsonin (IgG) | Senescence marker | Liver |
|---|---|---|---|---|
| Normal | 0.010 | 0.000 | **0.763** | 0.147 |
| Carrier | 0.373 | 0.000 | 0.252 | 0.349 |
| Mild HS | 0.672 | 0.000 | 0.044 | 0.280 |
| Moderate HS (ANK1) | **0.806** | 0.000 | 0.003 | 0.191 |
| Severe HS | **0.943** | 0.000 | 0.000 | 0.057 |
| Moderate HS (SLC4A1) | 0.599 | **0.266** | 0.001 | 0.134 |

That opsonin clearance is exactly zero in the spectrin/ankyrin genotypes is a structural
consequence: when band 3 leaves in the vesicles the surface density is conserved so no
clusters form, and IgG stays at the control level (45 molecules/cell). The opsonin arm
switches on **only in the band 3 genotype**.

Normal red cells die **because they are labelled** (senescence 76%), and HS red cells die
**because of their shape**. And the reason for that comes straight out of the geometry:
normal senescence loses area and volume at almost the same rate, so $D_c$ barely moves
(2.72 → 2.82 µm). A normal red cell cannot fail geometrically — which is why a molecular
label is needed.

### (2) Splenectomy fixes the anaemia and does not fix the cell — and the tests get **worse**

Moderate HS, before → after splenectomy: Hb 10.50 → 12.45, reticulocytes 8.75% → 3.83%,
lifespan 17.5 → 31.7 days, bilirubin 3.30 → 2.23, bilirubin production 932 → 642 mg/day.
But membrane area 110.7 → 108.4 µm², $D_c$ 3.316 → 3.756 µm, 50% haemolysis 0.571 → 0.655
%NaCl, MCHC 35.26 → 33.77. **Osmotic fragility does not improve; it worsens.** The reason
is structural: haemolysis depends on the **spleen-owned** term
$p_{\rm slow}(D_c)\times\pi_{\rm dest}$, while EMA and osmotic fragility depend on
$(S,V)$. After splenectomy the cell lives more than 100 days and keeps losing membrane, so
it becomes more spherical still — consistent with the observation that spherocytes persist
on the post-splenectomy blood film.

### (3) A parvovirus crisis is arithmetic

If production stops for $D$ days, $D/\text{lifespan}$ of the red cell mass disappears.
That alone explains why the same virus is harmless in a normal child and lethal in HS.

| | Baseline Hb | Lifespan | Arithmetic prediction $Hb\times 8/\text{lifespan}$ | 72-ODE model |
|---|---|---|---|---|
| Normal | 15.10 | 120 days | 1.01 | 0.72 |
| Mild HS | 13.16 | 42 days | 2.51 | 1.37 |
| Moderate HS | 10.50 | 17.5 days | 4.81 | 2.01 |

### Failure 1 — the amplifier hypothesis was not refuted but **unnecessary**

The model was designed around the idea that "the splenic cords do not merely eat
spherocytes, they **make** them" (area loss → $D_c$ rises → longer residence in the cords →
area loss). Switch that gain `cordamp` on and the phenotype worsens smoothly — but that is
**indistinguishable within the model** from raising `fdef`. And it was not needed either:
the phenotype of §2 already emerges at `cordamp = 1` (no amplification). So the reference
model switches the amplifier off.

**One real defect** came to light in this section and is recorded. The same parameter set
produced mean lifespans of 15.5, 18.6 and 30.8 days depending on the ODE integrator's
maximum step. The cause was not `cordamp` but **splenomegaly**. Because splenic blood flow
was written proportional to splenic mass, the gain of a second loop —
`hazard → erythrophagocytic load → splenic mass → blood flow → hazard` — exceeded 1, so the
spleen ran to its ceiling and the steady state was not unique. In a large spleen the flow
per gram falls — write flow $\propto$ mass$^{0.35}$ (capped at 1.6-fold) and the gain drops
below 1, and the model gives the same answer at max_step 4, 6 and 8 (§18).

### Failure 2 — genotype × splenectomy: the mechanism is reproduced and the direction is not

Reliene 2002 reported two things. ① band 3-deficient red cells carry up to 140 IgG per
cell, while spectrin/ankyrin-deficient cells and controls carry 60 or fewer. ② Splenectomy
helps spectrin/ankyrin deficiency **more**. In the model a single parameter separates the
two genotypes — **whether the shed vesicles carry band 3 out with them.**

① is reproduced: spectrin/ankyrin 45 (control level), band 3 151 (161 after splenectomy).
The mechanism is a single accounting rule. If band 3 leaves with the vesicle, band 3
**surface density** falls together with the area, so the **clusters** a low-affinity natural
antibody needs for bivalent binding are never formed at all.

② is not reproduced. In the model the band 3 arm benefits **more** from splenectomy (+3.48
against +1.95 g/dL). The cause is a single number that was assumed rather than measured —
the 72:28 spleen:liver split of opsonin clearance — and because the spectrin/ankyrin arm
has no opsonin clearance at all, this number acts only on the band 3 arm. So it can be
**back-calculated**: the liver would have to account for more than 97% of opsonin clearance
for the ordering of the Hb benefit to invert, and the ordering of the **lifespan** benefit
inverts at 65%. This is the quantitative implication of Reliene that the model reads out —
their result is not evidence about spherocytes but **evidence that clearance of heavily
IgG-coated red cells happens mostly in the liver**. The reference model leaves 28% as it is
and reports the failure.

### Failure 3 — bilirubin and Gilbert syndrome do not combine multiplicatively

`ugt_f = 0.287` is not a free parameter. It is the UGT1A1\*28 homozygous activity that
reproduces a Gilbert bilirubin of 2.2 mg/dL alone at a normal haemolytic load. Moderate
HS + Gilbert is then predicted at 11.2 mg/dL — twice the reported 4–7 mg/dL.
The model combines the two insults **multiplicatively** because conjugation is the only
processing step, and the data say they combine **sub-multiplicatively**. This is not a
tuning error but a structural statement: bilirubin handling must have a **load-inducible
component** (uptake, MRP2 excretion, alternative conjugates) that UGT1A1 kinetics alone
does not contain.

---

## 4. Defects found and fixed (things that came to light because the model was run)

1. **An order-of-magnitude error in the intracordal macrophage phagocytosis rate.** Writing
   `k_ph` as a per-day hazard made the normal red cell lifespan 0.44 days. Because a cell
   passes through the spleen 72 times a day, the destruction probability per pass has to be
   $\sim7\times10^{-4}$.
2. **Normal cells were inside the feedback loop.** At the initial `cordamp = 42` the loop
   gain exceeded 1 and the membrane of **normal** red cells collapsed over a few months.
3. **Volume loss was too slow.** `kd_base` was small, so normal senescent cells lost only
   area and $D_c$ moved into the filter. Fixed by making Waugh's −8.4% an explicit
   calibration target.
4. **The opsonin arm had no threshold.** A background IgG of 45 molecules/cell accounted for
   74% of clearance in normals. Reliene's control level (≤60) was introduced as a threshold.
5. **Dehydration protected the cell.** On $S/V$ geometry alone, losing volume makes $D_c$ go
   **down** — i.e. dehydration saves the cell. The clinic is the exact opposite. A
   cytoplasmic viscosity term had to be added, but attached to the **extraction step**
   rather than to the **residence time** (attach it to residence and dehydration amplifies
   membrane loss until even normal cells become supercritical).
6. **Intrasplenic cation leak is self-cancelling.** Put dehydration only in the cords and
   circulating MCHC cannot be raised to any value at all. A dehydrated cell is a cell about
   to be eaten, and what the blood count sees is the survivors in whom that **did not**
   happen. The high MCHC of HS is evidence that the leak is **constitutive** (a property of
   the membrane).
7. **The splenic mass–blood flow loop made the steady state non-unique** (§3 failure 1).
8. **Multivariate least squares kept converging to degenerate solutions.** Solutions with
   `kv_def → 0` that loaded all clearance onto `k_ph`, solutions with a viscosity
   coefficient 56 times the literature value, and so on. Because each target is monotonic in
   its own parameter, this was replaced by **one-dimensional bisection**.
9. **Intravascular haemolysis became an escape valve.** If the spleen did not remove the
   cell, the cell went all the way to a sphere and burst, with the result that 73% of HS
   haemolysis became intravascular (the real figure is <5%).
10. **The band 3 count at release did not track the area deficit.** Every genotype was born
    with a normal band 3 count laid on an already deficient membrane, so the surface density
    was high, and the spectrin/ankyrin arm, which should be at control level, came out at
    107 IgG/cell, **inverting the direction of the genotype × splenectomy prediction.**
11. **The integrator stepped over the transfusion window.** At `max_step = 2` the 0.5-day
    infusion window was never evaluated once, so the transfusion group quietly reported "no
    effect" (Hb 4.80 → 4.81).
12. **The reported spleen volume did not reflect the surgery.** The pre-splenectomy mass was
    displayed even after total splenectomy.

---

## 5. Honest limitations

- **There was no R toolchain in the environment.** `hsph_mrgsolve_model.R` and
  `hsph_shiny_app.R` **were not run.** The model of record is
  `hsph_python_reference.py`, which was run and calibrated, and the R files are an
  equation-for-equation transcription of it. The parameter values are injected
  mechanically from Python by `sync_r_params.py`.
- **The RDW falls far short** (11.7% against the screening criterion of ≥14%). In a
  structure where each of the 9 age cohorts has a single area and volume value there is no
  **within**-cohort variation, so the real population width cannot be generated. This is a
  structural limitation and cannot be fixed with parameters.
- **The circulating population always looks milder than the destroyed population.** The
  initial target of matching the EMA deficit at 16% against the circulating mean was
  **structurally incompatible** with a 20-day lifespan. This is not a bug but a result — a
  cell that deficient is already deep in the filter. The model reports the geometry at the
  moment of destruction (`A_end`) separately.
- **The severity mapping is off by one band.** At a spectrin content of 70% this model gives
  mild-to-moderate disease. So the genotype of the reference patient was fitted once (58%)
  and that deviation applied identically across all bands.
- **Stage 3 was calibrated at `fdef = 0.30` whereas the reference patient is 0.42.** As a
  result the predicted IgG of the band 3 arm exceeds the target of 140 by about 15%
  (151/161).
- **The neonatal section is the adult equations** run with neonatal UGT1A1 and EPO set
  points. There is neither a fetal-to-adult haemoglobin switch nor neonatal blood volume
  expansion. It is illustrative and not a neonatal model.
- **Mitapivat's** efficacy in HS is not pinned down in the literature. `EC50_m` comes from
  the published PK, and `Emax_atp`/`Emax_dpg` are placeholders taken from PK deficiency
  data.
- **The gallstone and iron sections** are semi-quantitative. The gallstone index is a
  dimensionless nucleation integral, and a single constant was matched so that untreated
  moderate HS reaches 50% at 30 years. The resulting 7.6% in normals (40 years) and 29.8%
  after splenectomy at age 6 (40 years) are predictions and are plausible, but HS + Gilbert
  at 92% (30 years) greatly overpredicts the reported roughly 5-fold excess risk — the same
  root as failure 3 in §3.
- **The HS + Gilbert combination in the neonatal section produces physically meaningless
  values** (peak bilirubin 82 mg/dL). Overlaying immature UGT1A1 and Gilbert as a product
  collapses the processing capacity, whereas real neonates have phototherapy, exchange
  transfusion and alternative excretion routes. This cell should be read only as showing at
  the extreme that "the model's bilirubin handling is multiplicative".
- **The virtual patient population (§17) is skewed mild.** Because the distribution was
  drawn from the literature's spectrin content ($fdef \sim N(0.28,0.09)$) without applying
  the 12 percentage point deviation found in §2. As a result 90% come out with Hb>11.
  Applying the deviation shifts the severity distribution to the right.

---

## 5b. Sensitivity — this is a geometry model

The top of the local sensitivity of Hb (each parameter ±20%, run output §17):

| Parameter | −20% | +20% | Spread |
|---|---|---|---|
| `A0` (membrane area at release) | 1.89 | 14.89 | **13.01** |
| `V0` (volume at release) | 13.58 | 2.12 | **11.46** |
| `D50` (filter position) | 9.08 | 12.53 | 3.45 |
| `w_esc` (residence time law) | 9.33 | 11.01 | 1.68 |
| `k_ph` (phagocytosis in the cords) | 11.18 | 9.84 | 1.34 |
| `cordamp` (the amplifier) | 10.68 | 10.32 | 0.36 |
| `tau50`, `k_sen`, `k_ops` | 10.50 | 10.50 | ~0 |

What dominates is **geometry** ($A_0$, $V_0$) and **the position of the filter**
($D_{50}$). The normal senescence constants (`tau50`, `k_sen`) have effectively no
influence on the Hb of moderate HS — because the cell does not live long enough to age.
The amplifier (`cordamp`) is nearly bottom of the list, which is consistent with the
conclusion of failure 1 in §3.

Transfusion introduces donor red cells with normal geometry, so it is simultaneously
haemoglobin replacement and **substrate dilution** (severe HS: Hb 3.06 → 8.73 g/dL, liver
iron 3.8 → 16.8 mg/g after 2 years).

---

## 6. How to reproduce

```bash
python3 hsph_geometry.py            # self-verification of the geometry kernel
python3 hsph_calibrate.py           # the 3-stage calibration → hsph_calibration.json
python3 hsph_analysis.py            # all 18 sections → hsph_reference_output.txt
python3 sync_r_params.py            # inject the Python parameters into the R $PARAM
python3 mkrefs.py --refresh         # re-query the references from PubMed
dot -Tsvg hsph_qsp_model.dot -o hsph_qsp_model.svg
dot -Tpng -Gdpi=150 hsph_qsp_model.dot -o hsph_qsp_model.png
```

On the R side:

```r
library(mrgsolve); library(dplyr)
mod <- mread("hsph_mrgsolve_model", ".")
mod %>% param(fdef = 0.42) %>% mrgsim(end = 1400, delta = 1) %>%
  plot(Hb + RETpct + MCHC + TBIL ~ time)
shiny::runApp("hsph_shiny_app.R")
```

---

## 7. Model structure (72 ODEs)

| Block | Number of states | Contents |
|---|---|---|
| 9 red cell age cohorts × 5 | 45 | cell number `N`, plus `N×area`, `N×volume`, `N×haemoglobin`, `N×band 3` (conserving the extensive quantities exactly under upstream transport) |
| Erythropoiesis | 5 | `EPO` `PROG` `ERB` `RETM` `RETB` |
| Transfused donor red cells | 1 | `NDON` (normal geometry — transfusion is also substrate dilution) |
| Spleen | 3 | `SPLV` `CORD` `MAC` |
| Haem metabolism · hepatobiliary | 7 | `HPT` `FHB` `BILU` `BILC` `BILE` `STONE` `LDH` |
| Iron | 4 | `HEPC` `FERR` `FELIV` `FESPL` |
| Mitapivat PK/PD | 5 | `MGUT` `MCEN` `MPER` `ATP` `DPG` |
| Crises · cofactors | 2 | `PARVO` `FOL` |

Age grid boundaries: 0 / 4 / 9 / 16 / 26 / 42 / 64 / 92 / 130 / 200 days.

The principal hazard terms:

$$h = \underbrace{f_{\rm pass}\,p_{\rm slow}(D_c)\,\bigl[1-e^{-k_{\rm ph}M\tau_c\eta}\bigr]}_{\text{geometric, spleen only}}
+ \underbrace{k_{\rm ops}\bigl(w_{\rm spl}\sigma M + w_{\rm liv}\bigr)\frac{x^m}{1+x^m}}_{\text{opsonin, spleen 72\% + liver 28\%}}
+ h_{\rm liv} + h_{\rm lys} + \underbrace{k_{\rm sen}(\tau/\tau_{50})^{m}}_{\text{senescence marker}}$$

where $\eta = e^{0.231(\mathrm{MCHC}-33)}$ is the cytoplasmic viscosity and acts **only on
the extraction step** (defect 5 in §4).

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It has not
been independently validated and must not be used for clinical decision-making,
prescribing, or regulatory submission.
