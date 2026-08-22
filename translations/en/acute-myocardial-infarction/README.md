# Acute Myocardial Infarction (STEMI) — QSP Model
## Acute Myocardial Infarction · Quantitative Systems Pharmacology

<p align="center">
  <a href="../../../acute-myocardial-infarction/ami_qsp_model.svg">
    <img src="../../../acute-myocardial-infarction/ami_qsp_model.png" width="900" alt="AMI QSP mechanistic map">
  </a><br>
  <sub><a href="../../../acute-myocardial-infarction/ami_qsp_model.svg">View the full-resolution SVG</a> · 246 nodes · 409 edges · 16 clusters (15 mechanistic + 1 legend)</sub>
</p>

---

## The claim this model is built to make

> **An infarct is not "the artery blocked and the muscle died".
> It is two races and one bifurcation.**

**Race 1 (minutes to hours).** The necrosis wavefront creeps from endocardium to
epicardium. Its speed is set by how far the residual collateral flow falls short of
the **basal metabolic demand** — basal demand, not resting demand. Tissue that can
pay the basal demand stops contracting and survives indefinitely (hibernation);
tissue that cannot, dies. This race is run against the time to reperfusion.

**Race 2 (the first few minutes after reperfusion).** Ischaemia leaves two things
behind. **Acid (H⁺)** holds the mitochondrial permeability transition pore (mPTP)
**shut**, and **succinate** becomes an oxidant burst the instant oxygen returns.
Restore the flow and the acid washes out and the succinate ignites — at the same
moment. Reperfusion is the treatment and simultaneously the second injury.

**The bifurcation (weeks to months).** Laplace wall stress → dilatation → still
greater wall stress is **positive feedback**. Concentric and eccentric hypertrophy
are the **negative feedback** that opposes it. And **scar cannot hypertrophy** —
compensation is the surviving myocardium's job. Whether a given infarct converges to
a stable ventricle or diverges into heart failure is therefore not a matter of dose
but **a race between two rates**.

As with the other models in this repository, what matters is not "what was put in"
but **"what came out although it was not put in"**. The following do **not exist** in
the parameter list:

- `golden_hour` — a golden-hour constant
- `critical_infarct_size` — a critical infarct size
- `no_reflow_switch` — a no-reflow switch
- `reperfusion_injury_fraction` — a reperfusion-injury fraction
- `therapeutic_window` — a cardioprotection window
- `severity` — a severity scale

Yet all of the results below come out.

| What was put into the model (structure) | What the model computed (result) |
|---|---|
| One transmural collateral gradient `GC` | the endocardium → epicardium necrosis wavefront, and a surviving epicardial rim |
| Comparing supply against **basal** demand | hibernating myocardium (akinetic but viable) arises by itself |
| Wavefront speed × time of reperfusion | the time–muscle curve is not a straight line but a **hyperbola** |
| Acid holds the pore shut, and acid washes out only with flow | the pH paradox, and a cardioprotection window of **minutes** |
| Succinate accumulation during ischaemia + reoxygenation | the oxidant burst becomes a **reperfusion event** (not an ischaemic one) |
| One term in which plasmin activates platelets | post-lysis reocclusion, and the reason P2Y₁₂ inhibitors exist |
| Patency rates of 60 % vs 95 % (observed in trials) | the lysis–PCI crossover point at **60–75 minutes** |
| Heart rate × blood pressure = demand | an intravenous beta-blocker reduces the infarct **only before** reperfusion |
| Marker release requires flow | reperfusion **brings forward** the marker peak (peak time 6.1 → 3.1 h) |
| Scar cannot hypertrophy (`1 − IS`) | the **separatrix** of remodelling appears as a computed result |
| ACEi and MRA merely block upstream of the receptor | GDMT moves not the infarct but the **separatrix** |

Every number was actually computed by
[`ami_reference_check.py`](../../../acute-myocardial-infarction/ami_reference_check.py),
and the full output is in
[`ami_reference_output.txt`](../../../acute-myocardial-infarction/ami_reference_output.txt).
**No number was written in by hand.**

---

## Files

| File | Contents |
|---|---|
| [`ami_qsp_model.dot`](../../../acute-myocardial-infarction/ami_qsp_model.dot) | Mechanistic map source — 16 clusters · 246 nodes · 409 edges |
| [`ami_qsp_model.svg`](../../../acute-myocardial-infarction/ami_qsp_model.svg) / [`.png`](../../../acute-myocardial-infarction/ami_qsp_model.png) | The rendered map (`dot -Tsvg` / `dot -Tpng -Gdpi=150`) |
| [`ami_mrgsolve_model.R`](../../../acute-myocardial-infarction/ami_mrgsolve_model.R) | mrgsolve model — **82 ODEs**, 5 transmural layers, 13 PK compartments, 7 reference scenarios |
| [`ami_reference_check.py`](../../../acute-myocardial-infarction/ami_reference_check.py) | numpy/scipy reimplementation — **the source of the numbers**. Regenerates every figure without R |
| [`ami_reference_output.txt`](../../../acute-myocardial-infarction/ami_reference_output.txt) | The full output of that script (16 experiments) |
| [`ami_shiny_app.R`](../../../acute-myocardial-infarction/ami_shiny_app.R) | Shiny dashboard — 11 tabs |
| [`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md) | 105 references + a table of **where the model disagrees with the literature** |

```bash
python3 ami_reference_check.py          # all experiments
python3 ami_reference_check.py 2 8      # experiments 2 and 8 only
dot -Tsvg ami_qsp_model.dot -o ami_qsp_model.svg
```

```r
library(mrgsolve)
mod <- mread("ami_mrgsolve_model.R")
mod %>% param(T_PCI = 1.5, T_ASA = 1, T_TIC = 1, T_HEP = 1,
              T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
        mrgsim(end = 24*180, delta = 0.05) %>% plot(IS + MVO_LV + EF + EDV ~ time)
```

---

## Model structure (82 ODEs)

The area at risk (AAR) is divided into **5 transmural layers**. The **only
difference between the layers is collateral flow** (`GC = 0.05 / 0.35 / 0.75 /
1.35 / 2.50`), and everything else about the wavefront is a consequence of that
gradient.

```
layer blocks (8 × 5 = 40)   E energy charge · G glycogen · H acid · C calcium ·
                            P mPTP opening · NI ischaemic necrosis · NR reperfusion necrosis · SUC succinate
global (29)                 ROS MVO THR PLN PLG FIB PAI DAMP NEU M1 M2 IL1 IL6 CRP
                            TGF MYOF COL MMP EDVS MASS THIN STUN NE ANG ALD BNP VOL
                            TNI CKMB
pharmacokinetics (13)       TNK×2 · ticagrelor×2 · metoprolol×2 · ramipril×2 ·
                            eplerenone · empagliflozin · colchicine · anti-IL-1β · CsA
```

`EDVS` is not the observed end-diastolic volume but the **structural (unloaded)
chamber size**. The observed EDV is `EDVS` displaced along the diastolic
pressure-volume curve according to the actual filling pressure, so that a single
variable holds both **day-1 dilatation** (incomplete emptying → raised filling
pressure) and **six-month dilatation** (growth of `EDVS` itself) together without
conflating them.

---

## Result 0 — an uninjured ventricle must not move

Before any remodelling result can be believed, **nothing must happen when there is
no infarct**. The wall-stress set point, contractility, and cardiac output are not
typed in but **derived** from the baseline geometry.

| | t = 0 | t = 365 days |
|---|---:|---:|
| EDV (mL) | 110.000 | 109.934 |
| ESV (mL) | 42.000 | 41.924 |
| EF (%) | 61.818 | 61.865 |
| Wall mass (g) | 150.000 | 150.180 |
| Wall stress | 9.240 | 9.223 |

**A one-year EDV drift of 0.06 %.** Every remodelling result below is read against
this baseline.

---

## Result 1 — the wavefront is not an assumption but a consequence of the collateral gradient

Complete occlusion, no reperfusion. Layer-by-layer collateral flow is
`0.005 / 0.035 / 0.075 / 0.135 / 0.250` and basal metabolic demand is `0.200`.
Nowhere is the model told which layer dies first.

| t (h) | L1 subendocardial | L2 | L3 midwall | L4 | L5 subepicardial | Infarct (%LV) | % of AAR |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 6.6 | 1.6 | 0.0 | 0.0 | 0.0 | 0.58 | 1.6 |
| 1.0 | 33.1 | 27.5 | 11.0 | 0.0 | 0.0 | 5.02 | 14.3 |
| 2.0 | 70.9 | 67.6 | 58.0 | 0.0 | 0.0 | 13.76 | 39.3 |
| 3.0 | 89.3 | 87.8 | 83.2 | 3.9 | 0.0 | 18.50 | 52.9 |
| 6.0 | 99.6 | 99.5 | 99.2 | 84.8 | 0.1 | 26.82 | 76.6 |
| 12 | 100.0 | 100.0 | 100.0 | 99.8 | 0.1 | 27.99 | 80.0 |
| 48 | 100.0 | 100.0 | 100.0 | 100.0 | **0.1** | 28.00 | 80.0 |

- **25 %** of the final infarct size is complete at **72 minutes**, **50 % at 123
  minutes**, **75 % at 219 minutes**, and **90 % at 291 minutes**.
- The 48-hour final infarct = **28.00 %LV = 80.0 % of the area at risk**.
- **Only 0.1 % of the subepicardial layer dies.** Because it can pay the basal
  demand. This is the "surviving rim" Reimer–Jennings measured in the dog, and it is
  not in the parameter list.

---

## Result 2 — the time–muscle curve is not a straight line but a hyperbola

Same model, same parameters. **Only `T_PCI` is moved.**

| Reperfusion (min) | Infarct (%LV) | Salvage index (%) | Ischaemic | Reperfusion | Reperfusion-injury fraction (%) | MVO (%LV) | EF at day 2 (%) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 15 | 0.04 | 99.9 | 0.00 | 0.04 | 98.3 | 2.84 | 54.5 |
| 30 | 0.99 | 97.2 | 0.48 | 0.50 | 51.2 | 2.32 | 53.2 |
| 45 | 5.05 | 85.6 | 2.65 | 2.40 | 47.6 | 3.03 | 52.3 |
| 60 | 8.90 | 74.6 | 5.06 | 3.85 | 43.2 | 3.55 | 50.5 |
| 90 | 13.71 | 60.8 | 9.73 | 3.99 | 29.1 | 4.13 | 48.2 |
| 120 | 16.17 | 53.8 | 12.87 | 3.29 | 20.4 | 4.46 | 47.6 |
| 180 | 19.18 | 45.2 | 16.24 | 2.94 | 15.3 | 4.99 | 47.1 |
| 240 | 21.81 | 37.7 | 18.52 | 3.29 | 15.1 | 5.64 | 46.6 |
| 360 | 26.51 | 24.2 | 22.84 | 3.67 | 13.9 | 7.05 | 45.7 |
| none | 28.00 | 20.0 | 24.14 | 3.87 | 13.8 | 20.25 | 45.1 |

**Muscle salvaged per 30 minutes of delay removed:**

| Interval | %LV / 30 min |
|---|---:|
| 30 → 45 min | **8.14** |
| 45 → 60 min | **7.74** |
| 60 → 75 min | 5.94 |
| 90 → 120 min | 2.46 |
| 240 → 300 min | 1.21 |
| 480 → 720 min | 0.03 |

Thirty minutes in the first hour is worth about **270 times** as much as thirty
minutes at the eighth hour. The model has never heard the phrase "golden hour" — the
wavefront has simply already passed.

**The reperfusion-injury fraction is non-monotone too**: it peaks at 43–51 % between
30 and 60 minutes and falls to 14 % beyond 240 minutes. Before the wavefront arrives
there is nothing to save, and after it has passed there is nothing left.

---

## Result 3 — the speed of the clock is set by collateral flow

| Collateral flow `COLL` | 30 min | 60 min | 90 min | 180 min | 360 min | none | Wavefront t50 (min) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.02 | 2.49 | 16.98 | 22.41 | 30.07 | 34.49 | 35.00 | 90 |
| 0.05 | 1.51 | 13.00 | 18.23 | 27.60 | 34.17 | 35.00 | 114 |
| 0.10 | 0.99 | 8.93 | 13.72 | 19.42 | 26.96 | 28.00 | 123 |
| 0.15 | 0.76 | 6.80 | 10.15 | 18.15 | 20.90 | 28.01 | 144 |
| 0.20 | 0.63 | 5.81 | 9.22 | 13.04 | 14.40 | 21.02 | 138 |
| 0.30 | 0.45 | 3.61 | 6.14 | 12.34 | 13.99 | 14.03 | 117 |

**At `COLL = 0.30` a six-hour delay (13.99 %LV) loses less muscle than a one-hour
delay at `COLL = 0.02` (16.98 %LV).** The "golden hour" is not a property of the
clock but **of that patient's collateral status**. At `COLL = 0.02` the surviving rim
disappears and 100 % of the area at risk dies.

---

## Result 4 — the pH paradox: the treatment takes the protection away

The pore is held shut by acid as `gate = 1/(1 + H/HP50)`, and acid is washed out
**only by flow**. Reperfusion at 90 minutes, subendocardial layer (L1), minutes
relative to the time of reperfusion:

| Relative to Trep (min) | H⁺ | Succinate | ROS | Ca²⁺ | Acid gate | mPTP | E | NR (% of layer) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −5 | 0.729 | 0.219 | 0.096 | 1.352 | 0.292 | 0.0395 | 0.000 | 1.82 |
| 0 | 0.688 | 0.199 | 0.131 | 1.476 | 0.304 | 0.0471 | 0.083 | 2.23 |
| +2 | 0.642 | 0.179 | 0.158 | 1.517 | 0.318 | 0.0519 | 0.199 | 2.39 |
| +10 | 0.441 | 0.101 | 0.227 | 1.443 | **0.405** | **0.0939** | 0.632 | 3.48 |
| +20 | 0.294 | 0.055 | 0.211 | 1.108 | **0.505** | **0.1430** | 0.850 | 5.51 |
| +40 | 0.143 | 0.021 | 0.124 | 0.507 | 0.678 | 0.1438 | 0.978 | 10.58 |
| +80 | 0.072 | 0.013 | 0.054 | 0.116 | 0.806 | 0.0441 | 0.999 | 15.45 |
| +160 | 0.060 | 0.012 | 0.050 | 0.038 | 0.834 | 0.0030 | 1.000 | 17.01 |

As the acid washes out the gate opens from 0.29 to 0.83, in the same window the
succinate burns and ROS peaks, and the pore opens **10–40 minutes after
reperfusion**. It is essentially over within 30 minutes. **There is no window
parameter in the model.**

---

## Result 5 — the cardioprotection window is minutes wide, and the prize moves

CsA does nothing but multiply the pore opening rate by `(1 − csa)`. The only thing
that changes between rows is **the timing**. Reference: PCI at 90 minutes, no CsA →
infarct 13.72 %LV (reperfusion component 3.99).

| CsA given | Infarct (%LV) | Reduction (%) |
|---|---:|---:|
| Trep − 30 min | 12.52 | **8.8** |
| Trep − 2 min | 12.73 | 7.2 |
| Trep + 0 min | 12.77 | 6.9 |
| Trep + 10 min | 13.08 | 4.7 |
| Trep + 30 min | 13.61 | 0.9 |
| Trep + 60 min | 13.72 | **0.0** |

And **the size of the prize depends on how long the ischaemia lasted**:

| Ischaemia (min) | No CsA | CsA | Reduction (%) | Reperfusion-injury fraction (%) |
|---:|---:|---:|---:|---:|
| 20 | 0.08 | 0.06 | **30.4** | 81.1 |
| 30 | 0.99 | 0.79 | 19.6 | 51.2 |
| 60 | 8.93 | 7.74 | 13.4 | 43.2 |
| 120 | 16.18 | 15.67 | 3.2 | 20.4 |
| 240 | 23.38 | 23.30 | **0.3** | 14.5 |

**The effect is non-monotone in ischaemic time.** A pragmatic trial that enrols a
wide band of ischaemic times dilutes its own effect size. This is **not** a claim
that ciclosporin works — CIRCUS and CYCLE were both neutral. It only shows that the
effect this model predicts is **exactly the kind a pragmatic trial is apt to lose**.

---

## Result 6 — the lysis versus PCI crossover point is a statement about patency *rates*

This is the most instructive experiment in the model. First hospital contact at
45 minutes.

**Part A — a single deterministic patient (artery opens in both arms).** The
crossover point comes at a PCI-related delay of **15–30 minutes**. Far shorter than
the guidelines' 60–120. It is not the model that is wrong but **the comparison**. A
single trajectory cannot express "in about 40 % of lysis patients the artery never
opens at all". And this is also where it becomes visible why partial patency is
almost as good as complete patency — once antegrade flow exceeds the basal demand
(0.20) the wavefront stops, so TIMI-2 flow alone arrests it. This is not a trick of
the model but physiology.

**Part B — the same comparison as a mixture of two subgroups.** The weights are not
model parameters but **the patency rates observed in trials** (about 60 % TIMI-3 at
90 minutes for a fibrinolytic with adjunctive therapy, about 95 % for primary PCI).

| | Infarct if successful | Infarct if failed | Expected value |
|---|---:|---:|---:|
| Lysis (60 % success) | 9.21 | 28.00 | **16.73 %LV** |
| PCI (95 % success) | see table | 28.00 | depends on delay |

| PCI-related delay (min) | E[infarct] PCI | E[infarct] lysis | Difference | Winner |
|---:|---:|---:|---:|:--|
| 0 | 6.48 | 16.73 | −10.25 | PCI |
| 30 | 12.70 | 16.73 | −4.03 | PCI |
| 60 | 15.72 | 16.73 | −1.01 | PCI |
| **75** | **16.77** | **16.73** | **+0.04** | lysis |
| 120 | 19.02 | 16.73 | +2.29 | lysis |

**A crossover at 60–75 minutes** — inside the guideline window (Nallamothu & Bates
60 min, Pinto 40–180 min). What matters is **where it came from**: not from flow
physiology but from **two patency rates**. The crossover point is a statement not
about **how well** each treatment works but about **how often** it works.

**Part C — reocclusion after lysis.** The raison d'être of antiplatelet therapy
comes out of one term (`KPLT_PLN`).

| Regimen | 90-min patency (%) | 24-hour patency (%) | Infarct (%LV) |
|---|---:|---:|---:|
| lysis alone | 48.2 | 6.3 | 28.03 |
| + aspirin | 55.9 | 11.8 | 28.03 |
| + aspirin + heparin | 61.3 | 24.4 | 35.00 |
| + aspirin + heparin + ticagrelor | 65.0 | **39.7** | **9.21** |

Because plasmin activates platelets, **the fibrinolytic partially undoes itself.**

---

## Result 7 — microvascular obstruction is not infarct size

MVO has its own equations (endothelial energy failure, leucocyte plugging, distal
embolisation, slow resolution) and **feeds back on flow**, so no-reflow sustains
itself. Nowhere is it forced to track infarct size.

| Scenario | Infarct (%LV) | MVO (%LV) | MVO/infarct | EDV at 180 d | EF at 180 d |
|---|---:|---:|---:|---:|---:|
| PCI 90 min, rich collaterals | 9.22 | 4.84 | 0.52 | 104.7 | 44.2 |
| PCI 90 min, poor collaterals | 21.85 | 6.22 | 0.28 | 146.2 | 42.2 |
| PCI 90 min, heavy embolisation | 13.91 | **13.33** | **0.96** | 116.4 | 44.3 |
| PCI 90 min, minimal embolisation | 13.73 | 3.17 | **0.23** | 114.0 | 43.9 |
| PCI 240 min | 23.38 | 8.09 | 0.35 | 172.3 | 37.7 |

**Two patients whose infarct sizes differ by only 0.18 %LV have MVO of 3.17 versus
13.33 %LV, a 4-fold difference.** Because the MVO/infarct ratio spreads from 0.23 to
0.96, MVO in this model is **not** a function of infarct size — that is the
experiment's first thesis.

> **The second thesis, however, the model fails to reproduce.** de Waha's 2017
> pooled CMR analysis reported that MVO predicts prognosis **independently of**
> infarct size. In this model a 4-fold difference in MVO changes the six-month EDV
> by only 2.5 mL (114.0 → 116.4). The reason is structural: because of `KMVO_R` (a
> resolution rate of 0.012/h, a half-life of about 2.4 days) the MVO disappears
> before it reaches the time scale of remodelling (weeks to months), and the model
> has no **permanent** sequelae pathway such as iron deposition from intramyocardial
> haemorrhage or microvascular rarefaction. Reported as a discrepancy rather than
> fitted away —
> [`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md) §14.

---

## Result 8 — the remodelling bifurcation: a critical infarct size comes out as a computation

Laplace: `σ = P·r/(2h)`. Dilatation raises `r` and so raises `σ` (positive
feedback); hypertrophy raises `h` and so lowers `σ` (negative feedback). And
**because scar cannot hypertrophy**, the brake is proportional to the surviving
fraction `(1 − IS)`. There is no GDMT in this table.

> **Read volumes only off the numbers on the stable branch.** Rows marked
> `DIVERGED` have reached the model's growth limit, and beyond that point what is
> being reported is **the verdict "this ventricle does not stabilise"**, not a
> prediction in millilitres. The model has no death and no other saturating
> mechanism.

| Reperfusion (min) | Infarct (%LV) | EDV 30 d | EDV 90 d | EDV 180 d | Growth (%) | EF 180 d | Verdict |
|---:|---:|---:|---:|---:|---:|---:|:--|
| 20 | 0.12 | 113.9 | 114.6 | 114.6 | 0.7 | 55.0 | **stable** |
| 40 | 3.71 | 116.8 | 119.8 | 121.1 | 3.6 | 52.0 | **stable** |
| 60 | 8.94 | 123.3 | 137.7 | 168.6 | 36.7 | 48.7 | diverging |
| 90 | 13.73 | 132.8 | 174.0 | — | — | 33.4 | **DIVERGED** (day 269) |
| 120 | 16.19 | 140.8 | 226.1 | — | — | 12.0 | **DIVERGED** (day 195) |
| 180 | 19.42 | 172.2 | — | — | — | — | **DIVERGED** (day 127) |
| none | 35.00 | — | — | — | — | — | **DIVERGED** (day 39) |

On the divergent branch **EF falls and BNP rises** — because once the hypertrophy
ceiling (`MASSMAX`) is reached the brake disappears and wall stress keeps climbing.
A ventricle that compensates and a ventricle that fails to compensate have
**qualitatively different trajectories**.

Sweeping the area at risk puts the **untreated separatrix at an infarct of about
4–5 %LV**: at AAR 0.05 (infarct 3.65 %LV) there is already 13.2 % growth, and from
AAR 0.10 (infarct 7.42 %LV) the growth ceiling is reached.

> **How this number should be read.** It is far smaller than the commonly cited
> critical infarct size of 18–20 %LV. But the cited figure comes from cohorts
> **receiving guideline therapy**. There is no GDMT in this table — that is, this is
> a pre-ACE-inhibitor-era ventricle. Result 9 tests exactly that difference, and
> **the treated separatrix moves to 18–20 %LV.**

---

## Result 9 — GDMT moves not the infarct but the separatrix

Guideline therapy starts at 24 hours — after the wavefront is over. It therefore
cannot change infarct size, and indeed does not. What changes is **whether the
mechanical loop diverges**. There is no anti-remodelling term in the model — only
receptor blockade upstream of the dilatation amplifier.

AAR 0.35 (infarct about 19.4 %LV), PCI at 3 hours, GDMT from 24 hours:

| Regimen | Infarct (%LV) | EDV 30 d | EDV 180 d | Growth (%) | EF 180 d | BNP | Ang II | Verdict |
|---|---:|---:|---:|---:|---:|---:|---:|:--|
| untreated | **19.42** | 172.2 | 834.9 | 384.7 | 29.2 | 252 | 2.71 | DIVERGED |
| ACEi alone | **19.43** | 120.9 | 171.8 | 42.1 | 42.1 | 1.84 | 1.42 | diverging |
| beta-blocker alone | **19.42** | 164.2 | 810.2 | 393.3 | 0.6 | 92.7 | 7.77 | DIVERGED |
| ACEi + BB | **19.42** | 131.2 | 175.7 | 34.0 | 43.1 | 1.78 | 1.31 | diverging |
| ACEi + BB + MRA | **19.42** | 127.2 | 151.5 | 19.1 | 44.0 | 1.30 | 1.31 | diverging |
| ACEi + BB + MRA + SGLT2i | **19.43** | 124.1 | 136.7 | **10.1** | **43.7** | 1.04 | 1.40 | **slowed** |

**Look at the infarct size column. It is identical to two decimal places
throughout.** A drug started at 24 hours cannot change a wavefront that is already
over. Yet the six-month EDV goes 835 → 137 mL and the EF goes 29 % → 44 %. And with
each drug added the growth rate falls **monotonically**, 42.1 → 34.0 → 19.1 →
10.1 %.

**And the separatrix really does move** (AAR sweep, PCI at 6 hours):

| AAR | Infarct (%LV) | Untreated growth (%) | Verdict | Full-GDMT growth (%) | Verdict |
|---:|---:|---:|:--|---:|:--|
| 0.08 | 5.90 | 37.1 | diverging | **−8.7** | **stable** |
| 0.12 | 8.94 | 96.5 | DIVERGED | **−7.0** | **stable** |
| 0.15 | 11.25 | 200.9 | DIVERGED | **−4.2** | **stable** |
| 0.20 | 15.13 | 305.1 | DIVERGED | 6.4 | slowed |
| 0.25 | 19.04 | 235.8 | DIVERGED | 18.9 | diverging |

> **The untreated separatrix at about 4–5 %LV → the treated separatrix at about
> 15–19 %LV.** The treated separatrix agrees with the critical infarct size cited in
> the literature — because the cited figure comes from treated cohorts. There is no
> anti-remodelling term in the model. There are only four constants: `KDIL_A` ·
> `KDIL_L` (the dilatation amplification by angiotensin and aldosterone) and
> `ACEI_EMAX` · `MRA_EMAX` (receptor blockade). The mechanism of SAVE, AIRE, and
> EPHESUS is contained in one table.

> **One place the model gets it wrong.** **Beta-blocker alone** **worsens**
> remodelling in this model (growth 393 %, EF 0.6 %). Because lowering heart rate
> lowers cardiac output, and that raises RAAS drive (Ang II 2.71 → 7.77). It is
> right as acute haemodynamics but the opposite of the clinical outcome
> (CAPRICORN, BHAT), and it is so because the model has none of the beta-blocker's
> antiarrhythmic, antioxidant, β1-receptor up-regulation-recovery, or direct
> reverse-remodelling effects. Combined with an ACEi the direction becomes right
> (ACEi alone 42.1 % → ACEi+BB 34.0 %), so the defect is confined to the monotherapy
> trajectory without RAAS blockade. Reported as it stands rather than smoothed away —
> [`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md) §14.

---

## Results 10–15 (summary)

The full numbers are in
[`ami_reference_output.txt`](../../../acute-myocardial-infarction/ami_reference_output.txt).

- **Experiment 10 · intravenous beta-blocker before reperfusion.** Metoprolol has
  **no** anti-necrotic term. Because perfusion occurs mostly in diastole, lowering
  heart rate lengthens diastole and so **increases collateral flow**, which delays
  the wavefront. The timing prediction is forced — given at 15 minutes the infarct
  falls by 2.0 %, after 45 minutes by 0.1 % or less, and after the artery is open by
  exactly 0. **The direction matches the split between METOCARD-CNIC (positive,
  before reperfusion) and EARLY-BAMI (neutral, later), but the magnitude is far
  smaller than the roughly 20 % reduction METOCARD-CNIC reported.** In this model
  the only route by which heart-rate control reaches the wavefront is collateral
  flow (basal metabolic demand does not depend on heart rate), and other mechanisms
  such as metoprolol's inhibition of neutrophil–platelet aggregation are absent.
  Reported as a discrepancy rather than fitted away — see
  [`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md) §14.
- **Experiment 11 · anti-inflammatory therapy — the model produces no effect at all
  (a negative result).** This experiment was designed to find the dose and timing
  optimum (a U-shaped response) of IL-1β blockade. **It was not found.** Sweeping
  the time of administration from 0.1 hours to 30 days and the strength of blockade
  from 0 to 0.995, the scar strength, thinning, six-month EDV, and EF are all
  identical to three decimal places. The reason is structural: the chain that builds
  the scar is M1 → M2 → TGF-β → myofibroblast → collagen and **does not pass through
  IL-1β**, and IL-1β contributes only 20 % of MMP activation. That is, there is no
  lever in this model for IL-1β blockade to grasp. **The clinical effects of CANTOS
  and COLCOT cannot be explained by this model.** The experiment is left in as a
  negative result rather than deleted — because it points to which mechanism needs
  to be added
  ([`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md) §14).
- **Experiment 12 · the washout phenomenon — only half of it is reproduced.**
  Because marker release is multiplied by flow, reperfusion **brings the marker peak
  forward**: peak time 6.1 hours with no reperfusion → 4.9 hours with PCI at
  240 minutes → 3.1 hours with PCI at 60 minutes. That is the time axis of the
  washout phenomenon and the model gets it right. **But it fails to reproduce the
  magnitude axis.** The peak troponin goes 39.7 → 50.2 → 53.5, preserving the same
  ordering as the actual infarct size, so the clinical observation that "the peak
  ranks patients in the wrong order" does not emerge. Moreover the AUC/infarct ratio
  also changes 2-fold, 135 → 68, so **the AUC too underestimates large infarcts.**
  Both defects appear to be because `WASH_FLOOR = 0.15` is too high, so that markers
  are released substantially even from occluded myocardium, and it was not adjusted
  to fit the data.
- **Experiment 13 · the price called bleeding.** Plasmin cannot tell coronary fibrin
  from intracranial fibrin. Halving the fibrinolytic dose roughly halves the
  bleeding indices too, but patency and muscle are lost along with them — the STREAM
  half-dose amendment in the elderly comes out of a single PK/PD chain. (The ICH
  index in the table is **a reporting figure**, not a validated risk model.)
- **Experiment 14 · structural sensitivity.** Mechanisms are cut one at a time. What
  carries the load is the **collateral gradient** (cut it and the wavefront loses its
  transmurality), the **basal-demand survival criterion**, and the **dilatation
  loop** (cut it and every infarct becomes harmless by six months). Removing the acid
  gate **shifts** necrosis from the reperfusion column to the ischaemic column
  without changing the total much — an internal-consistency check that the pore
  simply opens earlier.
- **Experiment 15 · the 7 reference scenarios.** These are the numbers the R model
  must reproduce.

| Scenario | Infarct (%LV) | MVO | Peak TnI | EF day 1 | EF day 180 | EDV day 180 | Scar |
|---|---:|---:|---:|---:|---:|---:|---:|
| S1 late presentation, artery not opened, untreated | 28.00 | 20.24 | 53.5 | 41.6 | 14.1 | *DIV | 0.654 |
| S2 primary PCI 90 min + DAPT, no GDMT | **13.72** | 5.36 | 47.7 | 49.4 | 33.4 | 326.1 | 0.740 |
| S3 primary PCI 90 min + full GDMT | **13.72** | 5.36 | 47.7 | 49.4 | **44.0** | **114.4** | 0.717 |
| S4 prehospital lysis 45 min + PCI 4 h + GDMT | 9.20 | 5.27 | 34.2 | 51.7 | 44.2 | 105.5 | 0.705 |
| S5 PCI 90 min + GDMT + IV BB + CsA + colchicine | **12.57** | 5.30 | 41.7 | 50.2 | 44.0 | 111.9 | 0.716 |
| S6 large AAR 0.45, PCI 4 h, full GDMT | 31.28 | 10.49 | 66.4 | 33.9 | 0.7 | *DIV | 0.708 |
| S7 small AAR 0.18, PCI 3 h, full GDMT | 9.85 | 3.46 | 22.3 | 56.1 | 47.2 | 98.2 | 0.671 |

  S2 and S3 have **the same infarct size** (13.72). The only difference is the drugs
  after 24 hours, and the six-month EDV goes 326 → 114 mL and the EF 33.4 → 44.0 %.
  S5 adds cardioprotection (an intravenous beta-blocker before reperfusion plus CsA
  just before reperfusion) and reduces the infarct 13.72 → **12.57 %LV** (8.4 %).
  `*DIV` marks a trajectory that reached the growth ceiling, and is a verdict rather
  than a volume in millilitres.

---

## Shiny app (11 tabs)

```r
install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
                   "DT","bslib","bsicons","scales"))
shiny::runApp("ami_shiny_app.R")
```

| Tab | Question it answers |
|---|---|
| 1 patient · lesion | who is this patient and what did the artery do |
| 2 pharmacokinetics | what is in the blood, and when |
| 3 necrosis wavefront | **race 1** — does the wavefront beat the clock |
| 4 reperfusion injury | **race 2** — acid washout versus oxidant burst |
| 5 microvascular obstruction | was the tissue reperfused, or only the artery |
| 6 inflammation · scar | is the scar built before the wall thins |
| 7 remodelling bifurcation | **the bifurcation** — convergence or divergence |
| 8 clinical endpoints | what would we have measured |
| 9 scenario comparison | same patient, different decisions |
| 10 biomarkers | troponin · CK-MB · CRP · BNP · fibrinogen |
| 11 sensitivity · structure | cut one mechanism and see what the model was standing on |

---

## What the model does **not** do (Deliberate limits)

This list exists to fix **what must not be asked of the model**.

- **It does not compute event rates.** The outputs go as far as infarct size, MVO,
  EF, EDV, BNP, and scar strength. The step from a surrogate to death or
  heart-failure hospitalisation is **not implemented.** Cluster 15 of the map is
  there to show the causal pathway.
- **There are no arrhythmia dynamics.** VF is an event, not a state variable, and it
  does not fit this model's time scale (hours to months).
- **There is no spatial structure of coronary anatomy.** The area at risk is a
  single scalar `AAR`, with no difference between LAD/RCA/LCx and no multivessel
  disease. There is no right ventricular infarction either.
- **The absolute position of the remodelling separatrix is not quantitatively
  calibrated** (the limitation of Result 8 above).
- **There are no comorbidities.** The covariates are only `AGE` · `PRECOND` ·
  `COLL`.
- **There is no population heterogeneity.** The mixture in experiment 6 is a
  weighting of two trajectories at **the reporting step**, not a stochastic process
  inside the ODEs.

The complete list, and the places where the model disagrees with the literature, are
in [`ami_references.md`](../../../acute-myocardial-infarction/ami_references.md)
§14–15.

---

## Disclaimer

A qualitative to semi-quantitative QSP model for education and research. It was
assembled on the basis of the public literature and clinical trial data but has not
been independently validated or certified, and **must not be used directly for real
clinical decision-making, prescribing, or regulatory submission.** The parameters and
assumptions are illustrative approximations, and fitting and validation against real
patient data would be needed separately.
