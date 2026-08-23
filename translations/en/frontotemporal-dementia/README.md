# Frontotemporal Dementia (FTD / FTLD) — QSP Model

> **One-sentence summary.** An anti-sortilin antibody restores plasma
> progranulin in GRN carriers perfectly to the normal range, yet clinical
> measures did not move. This model reproduces that failure as a
> **structural prediction** — because sortilin is simultaneously the
> **clearance pathway** that sets extracellular PGRN concentration and the
> **delivery pathway to the lysosomal pool** where PGRN actually acts,
> blocking it makes the ligand rise and the delivery route disappear **at
> the same time**.

---

## 1. Why FTD

Frontotemporal lobar degeneration (FTLD) is a leading cause of dementia
under age 60, and, unlike Alzheimer's disease, **the cholinergic system is
largely preserved**, so AChE inhibitors are ineffective and can even worsen
behavioural symptoms. Heritability is high, at roughly 30–40%, and three
Mendelian genes converge on **two distinct proteinopathies**:

| Gene | Mechanism | Pathology | Predicted onset age (observed) |
|---|---|---|---|
| **GRN** (17q21.31) | Progranulin haploinsufficiency | FTLD-TDP type A | 55.8 y (61.3) |
| **C9orf72** (9p21.2) | GGGGCC repeat expansion | FTLD-TDP type B + ALS | 54.3 y (58.2) |
| **MAPT** (17q21.31) | Tau mutation | FTLD-tau (**no** TDP-43) | 50.5 y (49.5) |
| Sporadic (~60–70%) | Heterogeneous | TDP-43 A–C / FET / tau | 59.2 y (~58) |

Observed values are the cohort mean from Moore 2020 *Lancet Neurol* (PMID
31810826, n=3403). **Onset age is an output of this model, not an input**
(§4 below).

---

## 2. Deliverables

| File | Contents |
|---|---|
| [`ftd_qsp_model.dot`](../../../frontotemporal-dementia/ftd_qsp_model.dot) | Mechanistic map source — **215 nodes · 313 edges · 14 clusters** |
| [`ftd_qsp_model.svg`](../../../frontotemporal-dementia/ftd_qsp_model.svg) | Vector rendering (zoom in to read) |
| [`ftd_qsp_model.png`](../../../frontotemporal-dementia/ftd_qsp_model.png) | Raster rendering (9706 × 5965) |
| [`ftd_mrgsolve_model.R`](../../../frontotemporal-dementia/ftd_mrgsolve_model.R) | **49-ODE** mrgsolve model + 13 scenarios + p-sensitivity analysis |
| [`ftd_reference_model.py`](../../../frontotemporal-dementia/ftd_reference_model.py) | Dependency-free RK4 reference implementation (same 49 states) — for validation |
| [`ftd_model_report.txt`](../../../frontotemporal-dementia/ftd_model_report.txt) | The **computed** output of the Python implementation above (the source of every figure below) |
| [`ftd_shiny_app.R`](../../../frontotemporal-dementia/ftd_shiny_app.R) | 9-tab interactive dashboard |
| [`ftd_references.md`](../../../frontotemporal-dementia/ftd_references.md) | 66 literature links — PMIDs clearly separated into **verified (✅ 44)** / **unverified (🔍 26)** |

Because R could not be run in this container, numerical validation of the
ODE system was carried out with the Python reference implementation. Every
number in `ftd_model_report.txt` is a **computed result, not a claim**.

---

## 3. The core of this model — sortilin's dual role

Progranulin is split into three compartments: **plasma / CSF / lysosome**.
Sortilin acts in two of them at once.

```
        f = sortilin's share of CSF PGRN clearance             (= 0.780)
        p = sortilin's share of lysosomal PGRN delivery         (= ???  never measured)
    theta = sortilin occupancy achieved by the drug

    CSF PGRN fold-rise       = 1 / (1 - f·theta)
    Lysosomal delivery ratio = (1 - p·theta) / (1 - f·theta)
```

This closed form agrees with the 49-state numerical integration (residual
≤0.02, report §C). So the following are **properties of the structure, not
coincidences of parameters**:

- The **sign of the effect is determined solely by sign(f − p)** — a
  property of the **biology**, not of the drug or the dose.
- **If p = f, the lysosomal pool does not change at any dose whatsoever.**
  The biomarker normalises perfectly while the delivery gain is
  mathematically zero.
- **If p > f, the drug raises the biomarker while lowering the target
  pool.**
- The ceiling as theta → 1 is **(1−p)/(1−f)**. **It cannot be crossed by
  dose** — raising occupancy makes the ligand rise and the delivery
  capacity disappear at the same rate.
- **Plasma PGRN always rises regardless of the value of p.** In other
  words, the licensing biomarker **cannot mathematically measure** the
  parameter that determines whether the drug works.

### 3.1 Validated figures (96 weeks, symptomatic GRN, latozinemab 60 mg/kg IV q4w)

| Compartment | Baseline | Post-dose | Fold |
|---|---|---|---|
| **Plasma PGRN** (licensing biomarker) | 35.0% | **94.5%** | 2.70× |
| CSF PGRN | 35.0% | 72.7% | 2.08× |
| **Lysosomal PGRN** (the actual site of action) | 35.0% | **49.1%** | 1.40× |

(% of wild-type for each compartment. Antibody CSF/plasma ratio **0.21%**,
peripheral occupancy 1.000, CNS occupancy 0.658.)

- 96-week CDR plus NACC-FTLD SB slowing: **0.65%** (−0.021 points out of
  +3.22 points of untreated progression). On a 0–24 scale, an effect under
  1% over 96 weeks is **undetectable by any realistically sized trial**.
  The model predicts **clinical failure sitting on top of clean target
  engagement**.
- At the same time the PD effects are real and measurable: CSF C1q −0.039,
  C3 −0.040, CSF BMP −19.5, lysosomal function +0.202. **Target engagement
  was not the problem.**

### 3.2 Sign reversal (a sweep over p)

| p | Lysosomal pool | Delivery ratio | Ceiling | **Plasma PGRN** | ΔCDR | Verdict |
|---|---|---|---|---|---|---|
| 0.00 | 73.5% | 2.100 | 4.545 | **94.5%** | −0.054 | benefit |
| 0.40 | 53.9% | 1.541 | 2.727 | **94.5%** | −0.027 | benefit |
| 0.50 | 49.1% | 1.402 | 2.273 | **94.5%** | −0.021 | benefit |
| 0.70 | 39.3% | 1.122 | 1.364 | **94.5%** | −0.008 | benefit |
| **0.78** | 35.4% | **1.011** | **1.000** | **94.5%** | −0.003 | **neutral (= f)** |
| 0.85 | 32.0% | 0.913 | 0.682 | **94.5%** | +0.001 | **harmful** |
| 0.95 | 27.1% | 0.773 | 0.227 | **94.5%** | +0.008 | **harmful** |

Look at the plasma PGRN column — **completely unchanged at 94.5% across
the entire sweep.** Whether the drug is helping, doing nothing, or making
the target pool worse, the biomarker looks identical. This is this model's
sharpest claim.

### 3.3 Gene therapy is not caught by this ceiling

AAV-GRN raises **production**, so the receptor is not cleared.

| Treatment | Plasma | CSF | **Lysosomal** | ΔCDR | ΔBMP |
|---|---|---|---|---|---|
| Untreated | 35.0% | 35.0% | 35.0% | — | — |
| latozinemab | **94.5%** | 72.7% | 49.1% | −0.021 | −19.5 |
| AAV-GRN (single dose) | 80.3% | 80.3% | **80.3%** | **−0.060** | −27.0 |

At a **lower plasma level**, gene therapy achieves 1.64× the lysosomal
recovery and 2.9× the predicted clinical effect. **The plasma PGRN
biomarker ranks these two treatments in the wrong order.**

→ The model's recommendation: ongoing progranulin-targeted trials should
add **lysosomal chemistry readouts (CSF BMP · glucosylsphingosine ·
cathepsin activity)** alongside plasma PGRN. This model simulates that
readout through the `BMPCSF` compartment.

---

## 4. Onset age is an output, not an input

A fixed genetic lesion alone cannot explain why a lifelong carrier stays
healthy for decades. So an explicit **age gate** was added (declining
autophagy/lysosomal capacity and microglial priming, `(age/62)^5`), and as
a result onset age became an output of the model.

| Gene | Predicted | Observed (Moore 2020) | Error |
|---|---|---|---|
| MAPT | 50.5 y | 49.5 y | **+1.0** |
| Sporadic | 59.2 y | ~58 y | **+1.2** |
| C9orf72 | 54.3 y | 58.2 y | −3.9 |
| GRN | 55.8 y | 61.3 y | **−5.5** |

The protective TMEM106B haplotype **delays** the predicted GRN onset from
55.8 to 57.1 y — consistent with the direction reported for rs1990622. The
model also wires TMEM106B into the lysosomal/ageing-vulnerability axis
rather than into PGRN concentration, which is consistent with a
7,071-person biofluid meta-analysis (PMID 38539243) reporting no
association between plasma PGRN and rs1990622. **This was not fitted — it
falls out of the wiring on its own.**

---

## 5. Other therapeutic axes (all validated figures)

| Treatment | Target engagement | Clinical effect | Interpretation |
|---|---|---|---|
| **C9orf72 ASO** | CSF poly-GP **−50.1%** (FOCUS-C9 reported ~50%) | ΔCDR **−0.035** | Target engagement is necessary but nowhere near sufficient. The same genotype also carries a C9orf72 **protein-haploinsufficiency** lysosomal penalty that the ASO cannot touch, and once started, the TDP-43 loop is self-templating. |
| **MAPT ASO** | Soluble tau **−45.7%**, mRNA −45.6% | Aggregated tau only **−17.4%**, ΔCDR −0.139 | Blocking synthesis only **starves** the self-templating deposits — it does not remove them. → Should be tested in presymptomatic/early cohorts. |
| **Trazodone** | 5-HT2A/SERT | ΔNPI **−5.8** (Lebert 2004: −6 ~ −10) | The only clearly positive behavioural RCT in FTD. |
| **SSRI** | SERT blockade | ΔNPI −6.2 | Reverses only the presynaptic serotonin deficit. Postsynaptic 5-HT2A loss cannot be reversed, so the effect stays partial. |
| **Donepezil** | AChE **blocked** (cholinergic tone 0.95 → 1.66) | **ΔCDR exactly 0.000**, ΔNPI **+3.2** | **The null result comes from the structure.** `ACH_INTEGRITY = 0.95` — because the cholinergic system is preserved in FTD, tone was never the rate-limiting variable for the cognitive index. The drug raises tone, but raises something CDR does not depend on, while only the agitation term is engaged, so behaviour worsens. A mechanistic contrast with AD. |

---

## 6. Natural-history validation

| Genotype | CDR slope/year | Plasma NfL | CSF NfL | Atrophy %/year | NPI |
|---|---|---|---|---|---|
| MAPT | 1.58 | 42.7 | 2484 | 2.16 | 22 |
| C9orf72 | 1.67 | 59.5 | 3402 | 2.20 | 25 |
| GRN | 1.76 | 57.5 | 3290 | 2.25 | 24 |
| Sporadic | 1.69 | 56.1 | 3214 | 2.18 | 23 |
| **Target** | **1.5–2.5** | **50–80** | **3000–5000** | **2–3** | **20–40** |

All six metrics fall within the literature range (only MAPT's NfL is
slightly low).

---

## 7. Where this model is wrong (reported without concealment)

In order of importance. All eight items are in report §J.

1. **No human measurement of `P_LYS_SORT` (p) exists.** The **entire
   sign** of the anti-sortilin prediction depends on this value. The
   default of 0.50 is below f_CSF=0.78, so it was set **on the side
   favourable to the drug**; set it adversarially and the drug is,
   frankly, harmful. What this model constrains is not the **value** of p
   but the **relationship** between p and f.

2. **There is no compensated latency period — this single defect creates
   two discrepancies.** Because damage accumulates monotonically from
   t=0, (a) GRN onset comes 5.5 years too early, and (b) the
   presymptomatic NfL lead time comes out at **~15 years** instead of the
   reported ~1–2 years (at a 1.5× threshold; it shrinks to 2.5 years at a
   4× threshold). This can only be fixed with **a different structure**,
   not a parameter — a compensatory regime is needed that holds
   biomarkers near normal until a threshold is crossed, then collapses.
   **This is the single most valuable improvement anyone could make to
   this model.**

3. `KSYN_NMD = 0.70` is a **fitted stand-in** used to reconcile plasma
   PGRN at ~1/3 of controls with an allelic dose of 1/2. Its consequence
   is real, though: plasma falls to 33% and CSF to 50%, so **a single
   gene-dosage scalar cannot fit both at once** — evidence that the two
   compartments genuinely have different clearance structures.

4. CNS antibody exposure is treated as a single homogeneous CSF
   compartment. Real cortical parenchymal distribution is regionally
   uneven, and **the salience network happens to be exactly where
   delivery matters most**. If cortical exposure is lower than the CSF
   surrogate, theta_CNS is an overestimate, and every anti-sortilin
   figure here is optimistic.

5. poly-GP is treated as a **marker** in equilibrium with repeat RNA, and
   poly-GR as **the toxic species**. If the causal DPR is actually
   poly-GA (proteasome damage), the ASO predictions change, and the CSF
   marker trials actually use ends up tracking the wrong arm of the
   pathway.

6. Two lumped regions cannot represent what actually determines the FTD
   phenotype — **which network fails first**. The svPPA scenario is a
   weight shift, not a mechanism.

7. Sporadic FTLD was given the same TDP-43 machinery as GRN, with no
   genetic driver, so its trajectory depends solely on the age gate. Real
   sporadic FTLD spans TDP-43 A–C, FET, and tau, and cannot be captured
   by a single parameterisation.

8. The age gate `(age/62)^5` is phenomenological. It does real work
   turning onset age into a prediction, but because the exponent was
   fitted to onset data, the good agreement in §4 is partly
   **calibration, not prediction**. What was not fitted is the
   **ordering of genotypes** and the **direction of the TMEM106B
   effect**.

---

## 8. Mechanistic map (14 clusters)

[![FTD QSP map](../../../frontotemporal-dementia/ftd_qsp_model.png)](../../../frontotemporal-dementia/ftd_qsp_model.svg)

1. Genetic architecture and risk modifiers (GRN · C9orf72 · MAPT · TBK1 · VCP · CHMP2B · TMEM106B)
2. Progranulin biology and sortilin trafficking ← **the drug-target axis**
3. Lysosomal / autophagy / proteostasis failure
4. TDP-43 proteinopathy and cryptic splicing (STMN2 · UNC13A · HDGFL2)
5. C9orf72 repeat toxicity — RNA foci · RAN translation · the three DPR pathways
6. Tau pathology (exon-10 splicing · 4R/3R · aggregation · propagation)
7. Neuroinflammation · microglia · **complement-mediated synaptic pruning**
8. Synapse loss · neuronal death · selective vulnerability (**von Economo neurons**)
9. Large-scale network degeneration and regional atrophy
10. Neurotransmitter systems — serotonergic loss versus **the preserved cholinergic system**
11. Disease-modifying drug PK/PD (anti-sortilin · AAV-GRN · C9 ASO · MAPT ASO)
12. Symptomatic pharmacology — **including the cholinergic null result**
13. Fluid and imaging biomarkers
14. Clinical phenotypes · endpoints · outcomes

---

## 9. Usage

```r
source("ftd_mrgsolve_model.R")

# Predict onset age (Moore 2020: MAPT ~50, C9 ~58, GRN ~61)
for (g in c("MAPT", "C9", "GRN", "sporadic")) {
  bi <- FTD_burnin(FTD_genotype(g))
  cat(sprintf("%-9s predicted onset age %.1f y\n", g, bi$onset_age))
}

# 13 scenarios
res <- FTD_simulate_scenarios()

# Sign reversal and invariance of the plasma biomarker
FTD_sweep_p()
```

```bash
# Recompute every figure with the reference implementation (no R needed)
python3 ftd_reference_model.py            # full report
python3 ftd_reference_model.py --brief    # headline figures only

# Re-render the map
dot -Tsvg ftd_qsp_model.dot -o ftd_qsp_model.svg
dot -Tpng -Gdpi=75 ftd_qsp_model.dot -o ftd_qsp_model.png
```

```r
shiny::runApp("ftd_shiny_app.R")   # 9-tab dashboard
```

---

## ⚠️ Disclaimer

This model is a **qualitative/semi-quantitative QSP model for educational
and research purposes.** It was built from published literature and
clinical trial data but has not been independently validated or
certified, and **must not be used directly for clinical decision-making,
prescribing, or regulatory submission.** In particular, `P_LYS_SORT` has
never been measured in humans, and it determines the sign of every
anti-sortilin-related conclusion in this model. The descriptions here of
commercial development programmes are modelling interpretations based on
public reporting, and are not themselves clinical conclusions.
