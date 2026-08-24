# Prolactin-Secreting Pituitary Tumour (Prolactinoma) — QSP Model
### Prolactinoma / lactotroph PitNET · drug-induced hyperprolactinaemia · stalk effect · macroprolactinaemia, in a single equation

<a href="prl_qsp_model.svg"><img src="prl_qsp_model.png" width="720" alt="Prolactinoma QSP mechanistic map"></a>

| Deliverable | File | Scale |
|---|---|---|
| 🗺️ Mechanistic map | [`prl_qsp_model.dot`](prl_qsp_model.dot) · [`.svg`](prl_qsp_model.svg) · [`.png`](prl_qsp_model.png) | **236 nodes · 21 clusters · 339 edges** |
| ⚙️ mrgsolve ODE model | [`prl_mrgsolve_model.R`](prl_mrgsolve_model.R) | **58 ODEs · 179 parameters · 34 scenarios · 15 diagnostic analyses** |
| 📊 Shiny dashboard | [`prl_shiny_app_en.R`](prl_shiny_app_en.R) | **11 tabs** |
| 📚 References | [`prl_references_en.md`](prl_references_en.md) | **150 PubMed citations (all looked up and verified via NCBI)** |

All simulations were run and verified on actual **mrgsolve 2.0.1 / R 4.3.3**.

---

## 1. The organising idea

> **Prolactinoma is the only pituitary tumour under inhibitory neuroendocrine control.
> That is why a receptor agonist is not merely a secretion suppressant but an
> antitumour agent.**
>
> **And among endocrine diseases, it is the only one in which the measured number
> and the signal the body sees can diverge, in either direction, by two orders of
> magnitude.**

Every other cell of the anterior lobe is *driven* by a releasing hormone, but the
lactotroph is uniquely *held down* by one — dopamine, descending through the portal
vessels from hypothalamic TIDA neurons. Remove dopamine and prolactin rises; add a
higher-affinity agonist and it falls. Unusually, **the same receptor also decides
whether that cell divides or survives.** This is why prolactinoma is the only
pituitary tumour worldwide for which medical therapy, not surgery, is first-line
treatment.

---

## 2. One occupancy equation, four directions

In the model, every dopaminergic ligand enters through **a single competitive
occupancy equation**. Each ligand is weighted by its own intrinsic efficacy, e.

```
SIGDRIVE = Σᵢ eᵢ (Cᵢ/Kᵢ) / (1 + Σⱼ Cⱼ/Kⱼ)      ×  D2R density
```

| Ligand | e | Outcome |
|---|---|---|
| Dopamine (portal) | 1.00 | Physiological tone |
| Cabergoline · quinagolide | 1.00 | SIGDRIVE ↑ → prolactin ↓ |
| Bromocriptine | 0.80 (partial) | **Structural ceiling** — cannot be overcome by dose |
| Aripiprazole | 0.25 (partial) | High affinity means it *exceeds* physiological tone even alone |
| Risperidone · paliperidone · amisulpride · haloperidol | 0.00 | Displaces dopamine → prolactin **rises** |

From this single line, the following are **derived rather than coded in as rules**:

- **Antipsychotic-induced hyperprolactinaemia**: an e=0 ligand displaces an e=1
  ligand → signal falls → prolactin rises (at risperidone-equivalent 6 mg/day,
  prolactin goes 54 → 73 ng/mL, continuing to climb over months because chronic D2
  blockade drives lactotroph hyperplasia)
- **Aripiprazole alone**: low efficacy but high affinity means it displaces
  dopamine (e=1) while still raising net signal, so prolactin falls **below
  normal** (model value 4.3 ng/mL) — the prolactin-lowering effect of
  aripiprazole reported clinically
- **Aripiprazole add-on**: adding it on top of risperidone replaces e=0 with
  e=0.25 → **signal rises → prolactin falls** (SIGDRIVE 0.0141 → 0.1503, a
  10.6-fold change). Not a new mechanism, just arithmetic
- **Stalk effect**: here it is the *delivery* term (portal dopamine) that is
  interrupted rather than the receptor term → only the normal cells are
  disinhibited, producing a **calculable ceiling**

---

## 3. One receptor, four branches, four clocks

D2 occupancy acts through four branches, each with its own IC50 and its own time
constant. The ordering is deliberate.

| Branch | Pathway | Time constant | Outcome |
|---|---|---|---|
| **① Secretory granule release** | Gi/o → GIRK K⁺ → hyperpolarisation → Ca²⁺ ↓ | Hours | The sharp first-day fall |
| **② Transcription** | cAMP/PKA ↓, ERK ↓ → Pit-1/ERα drive ↓ | Days | **Determines the chronic level** |
| **③ Cell volume** | Granule depletion + RER atrophy | About 2 months | **Reversible** shrinkage |
| **④ Proliferation/apoptosis** | PI3K/Akt·MAPK withdrawal, p27 ↑, TGF-β1 ↑ | Months to years | **Banked** shrinkage |

What the model *generates* from this structure:

- Prolactin falls to 66% of baseline within a day and to 22% by day 7, but **the
  chronic level is set by branch ②** — because of mass conservation at steady
  state (secretion = synthesis − degradation), not by assumption.
- Tumour shrinkage lags biochemical control by months.
- **Raising the dose in a macroadenoma whose prolactin is already normal but has
  not shrunk is rational** — because the four branches sit at different points on
  the same occupancy axis.
- Splitting shrinkage into **volume (reversible)** and **cell number
  (irreversible)** shows that an MRI that looks worse a month after stopping the
  drug reflects *cell re-swelling, not regrowth* (diagnostic analysis D14).

### An honest negative result — branches ① and ② are not as temporally separated as the map suggests

The actual output of D5 shows that both receptor branches **settle within 3
days**. What separates them is *role*, not *timing*. This weakens the map's
narrative, but it is reported as found.

---

## 4. The regulator is already saturated at the time of diagnosis

Prolactin drives its own inhibition through a short-loop feedback: PRL →
PRLR/JAK2/STAT5 in TIDA neurons → tyrosine hydroxylase ↑ → portal dopamine ↑. In
prolactinoma, **this loop is normal. What is broken is the receiver (reduced D2R
density).**

Two conclusions the model generates:

1. At the time of diagnosis, TIDA drive is already at **68%** of its maximum
   amplification (portal dopamine 5.0 → 11.8 nM). With no endogenous reserve left
   to recruit, **the only way to add more signal is a ligand with far higher
   affinity than dopamine itself** (cabergoline Ki 0.7 nM vs dopamine's apparent
   50 nM at pituitary D2).
2. **But this loop cannot restrain residual tumour.** See §7 below — the point at
   which this model refutes its own design hypothesis.

---

## 5. The measurement layer is part of the pathophysiology

The three classic management errors in this disease are not biological failures
at all — they are failures of **the mapping between biological prolactin and
reported prolactin**. So that mapping was implemented explicitly, as a separate
layer.

```
PRLIMM   = PRLB·(1+FDIM) + XMAC·PRLM                  ← the true analyte
PRLMEAS  = PRLIMM / (1 + (PRLIMM/KHOOK)^PHOOK)        ← two-site sandwich assay reading
PEGREC   = 100 · PRLB·(1+FDIM) / PRLIMM               ← PEG-precipitation recovery
```

| Pitfall | Direction | Model result |
|---|---|---|
| **Hook effect** | False **low** reading | A macroadenoma's true prolactin of 10,769 ng/mL (analyte 11,846) is **reported as 250 ng/mL**, and a 1:100 dilution recovers **9,998** |
| **Macroprolactinaemia** | False **high** reading | Reported 62.5 / biological 9.6 ng/mL, **PEG recovery 16.9%**, gonadal axis normal (ovulatory capacity 1.00, T-score 0.00) |
| **Glycosylation · 16 kDa cleavage** | Activity mismatch | Shown on the map as a separate fraction/pathway |

The hook is a curve, not a rule: almost linear up to an analyte of 2,500, peaking
around 3,000, and collapsing above that. **The confirmatory 1:100 dilution
evaluates the same curve 100-fold lower on the antigen axis** — structurally
identical to excess heparin abolishing the serotonin-release assay in this
repository's HIT model.

---

## 6. Safety is a selectivity ratio read at a different concentration

Efficacy is driven by the **pituitary biophase** (61.5-fold the plasma
concentration, with an equilibration half-life of ~4 days), while valvular risk
is driven by **plasma** 5-HT2B occupancy, because the valve has no such
reservoir. The two are therefore not proportional.

| Receptor | Outcome | Model value |
|---|---|---|
| D2 (pituitary biophase) | Prolactin control · tumour shrinkage | 6.83 nM at 1 mg/week = 9.8× Ki |
| **5-HT2B (plasma)** | Ergot valvulopathy | Prolactinoma at 1 mg/**week**: P(TR) **0.38%** vs Parkinson's disease at 3 mg/**day**: **25.8%** |
| D3 (mesolimbic) | Impulse-control disorder | Already **14.1%** at 1 mg/week — independent of the ergot scaffold |
| D2 (area postrema) | Nausea (including the tolerant state) | Day 1: 18.4 → day 84: 2.9 |

**The Parkinson's disease literature and the prolactinoma literature were never
actually in conflict — they simply were never on the same dose axis.** Sweeping
5-HT2B Ki across four orders of magnitude (1.2 → 1200 nM) leaves prolactin
control and tumour shrinkage **completely unchanged** while collapsing valvular
risk alone (D12). Quinagolide, by contrast, sits at the far end of that axis
because it lacks the ergot scaffold, yet lands **at the same point as
cabergoline on the D3/impulse-control axis** — which is why switching to it for
that purpose does not help.

---

## 7. Where this model refutes its own design hypothesis (D10)

The mechanistic map argues that "remission after stopping the drug occurs when
the residual mass is small enough for endogenous dopamine tone to restrain it."
The equations **reject this claim.**

```
The short loop is driven by prolactin. When the residual is small, prolactin is low,
so TIDA drive also returns to baseline and portal dopamine falls with it.
→ Endogenous tone is weakest exactly when the residual is smallest — the opposite of what is needed.
```

Measured: for a small residual, PRL 13.3 → portal DA 5.64 nM → SIGDRIVE 0.1014;
for a large residual, PRL 52.3 → portal DA 10.30 nM → SIGDRIVE 0.1708.

So this model **never produces a cure.** Wherever residual cells remain, relapse
always follows; a smaller residual only delays it (seeding 0.002 mL → 8.9 years,
3.2 mL → immediate). Time to relapse is **log-linear in the residual, with no
threshold.**

This reinterprets Dekkers 2010's pooled sustained-remission rate of ~21% **not as
a statement about cure, but as a statement about follow-up duration**. And it
leaves a testable prediction: **with 15 years of follow-up, biochemical relapse
should be near-universal, and any patient who does not relapse is evidence of a
mechanism absent from this model (senescence, infarction, immune-mediated
atrophy).** This is reported as found, not removed by adjusting parameters.

---

## 8. A generated phenotype (not imposed)

Every disease scenario first simulates **2 years of untreated, asymptomatic
natural history (PRESIM)** without any drug, and reads the state at that point
as "presentation". None of prolactin, amenorrhoea, bone density, or visual field
was entered by hand. Diagnostic analysis D1 reports a maximum 10-year drift of
**0.000060%** in the absence of a tumour — because every baseline equilibrium is
solved algebraically in `$MAIN`.

| Untreated natural history | PRL at presentation | PRL at 7 years | Volume | Visual field MD | T-score | Anovulation |
|---|---|---|---|---|---|---|
| Healthy control | 10.0 | 10.0 | 0.00 | -0.0 | -0.00 | 0% |
| Microadenoma 0.35 mL, D2R 0.50 | **43** | 131 | 0.41 → 0.83 | -0.0 | -1.11 | 83% |
| Macroadenoma 3.2 mL, D2R 0.35 | **784** | 3084 | 4.25 → 11.15 | -15.1 | -1.63 | 100% |

The published size-prolactin ordering (macroadenoma >250, microadenoma
50–250 ng/mL) falls out of **mass × receptor density**. Taking the normal
lactotroph pool as 0.15 mL, a 3.2 mL tumour is 21 times normal mass, multiplied
further by the deficient receiver.

---

## 9. Scenario summary excerpt (34 total)

| # | Scenario | Outcome |
|---|---|---|
| S4-S6 | Cabergoline 0.5 / 0.25×2 / 1 mg/week | Microadenoma 43 → 6.1 (77% shrinkage), macroadenoma 784 → 58 (73% shrinkage) |
| S5 | **Twice-weekly split vs once weekly** | 15% difference in biophase trough, **no difference** in prolactin → does not support splitting the dose for efficacy |
| S7 | Macroadenoma vs borderline compression | 72-hour visual field recovery is seen **only in the borderline case** — see §11 |
| S8-S10 | Bromocriptine / switch / quinagolide | SIGDRIVE 0.650 vs 0.908 vs 0.951 — the partial agonist's structural ceiling |
| S11 | Partial resistance vs true resistance | A 7-fold dose increase rescues an EC50 shift but not a reduced Emax |
| S12 | 10-year clonal escape | Resistant fraction 3.7% → 99.6%, generated by Darwinian selection |
| S15-S16 | Pregnancy | Drug withdrawal + ~300-fold E2 + placental lactogen. **Most of the increase is cell re-swelling** |
| S17 | Stalk effect | Portal DA 0.55 nM, reported PRL 85 — a calculable ceiling |
| S21 | **Macroprolactinaemia + cabergoline (structural null)** | Reported 62.5 → 10.5, **no change** in ovulatory capacity or T-score |
| S25-S26 | Surgery ± cabergoline | A 25% residual regrows on its own but not under dopaminergic control |
| S27-S28 | Temozolomide, MGMT 0.05 vs 0.80 | 78% shrinkage vs no response (indistinguishable from no alkylating agent) |
| S29-S30 | Cabergoline valvulopathy: weekly vs daily dosing | P(TR) 0.38% vs 25.8% |
| S31-S32 | Bone: what recovers and what is banked | About 62% of the 5-year loss is irreversible; sex-hormone replacement protects bone even without tumour control |

## 10. Diagnostic analyses (15 total)

`D1` Baseline steady-state (0.000060% drift) · `D2` Cabergoline PK (Cmax 68.6
pg/mL, t½ 107 h) · `D3` 22% of baseline at 14 days after a single dose · `D4`
Cabergoline vs bromocriptine dose-response (plateau 50.9 vs 80.5 ng/mL) · `D5`
Decomposition of the four branches · `D6` Hook curve · `D7` PEG recovery · `D8`
Stalk-effect ceiling (94 ng/mL with complete transection, 1/8 of the 3.2 mL
tumour) · `D9` Aripiprazole add-on paradox · `D10` Discontinuation
(self-refutation) · `D11` Valvulopathy dose axis · `D12` 5-HT2B/D2 selectivity
sweep · `D13` Geometry of pregnancy risk · `D14` Volume first, cell number later
· `D15` The permanent cost of diagnostic delay

---

## 11. Structural nulls and negative results (reported, not removed)

**Cases where the model should fail treatment, and does:**

- Macroprolactinaemia + cabergoline → only the number moves; the axis is
  unchanged
- True (reduced-Emax) resistance + 7-fold dose increase → no rescue
- MGMT-high tumour + temozolomide → no response
- Switching to quinagolide to avoid impulse-control disorder → no protection

**Negative and self-critical findings:**

1. **Does not produce remission** — discontinuation always ends in relapse,
   refuting the map's design hypothesis (§7).
2. **Increasing the dose beyond about 0.5-1 mg/week buys almost nothing
   biochemically** — the occupancy term is already saturated. The model argues
   **against** part of the practice it set out to explore (D4). (Shrinkage,
   however, sits on a less saturated branch, so dose escalation for that
   purpose remains justified.)
3. **Does not support twice-weekly split dosing** — except for
   tolerability/nausea purposes (S5).
4. **Does not reproduce 24-72-hour visual field recovery in macroadenomas.** The
   reason is geometric and hard to argue with: because volume scales with the
   cube of radius, a 3% volume reduction cannot move the apex of a 4 cm mass by
   even 1 mm. In the model, rapid recovery appears **only when the optic chiasm
   sits near a threshold**, which is likely the true substance of the case
   reports. No term was added to manufacture a fast response in macroadenomas
   (S7).
5. **Cabergoline drives prolactin below normal in most simulated patients** —
   not a fitted result but a **prediction of iatrogenic hypoprolactinaemia**.
6. **Population normalisation rates (Webster's 83% vs 59%) cannot be produced by
   a single-patient model.** What is reproduced is only the direction and the
   plateau difference; the rates themselves would require a population
   distribution of D2R density.
7. **A stalk-effect ceiling above 94 ng/mL cannot be produced by this
   structure.** If the values above 150 ng/mL reported in some series are real,
   they require hyperplasia larger than the model allows, or a co-existing
   lactotroph adenoma — this is stated as a **falsifiable claim, not a
   caveat**.

---

## 12. Limitations

- Cabergoline's **absolute bioavailability is unknown**, so `V2_CAB` and
  `CL_CAB` are apparent values, fitted to reproduce the published pg/mL range
  and the 63-109 h half-life. The biophase partition coefficient is not
  separately identifiable from plasma data alone.
- The hook constants `KHOOK` and `PHOOK` are **platform-specific**. The only
  general claim is the shape of the curve.
- **The menstrual cycle is not modelled.** Ovulatory capacity is a graded
  index, not a cycle simulator.
- A single "antipsychotic" compartment stands in for an entire class. Varying
  `EFF_AP` and `KI_AP` represents risperidone (e=0), aripiprazole (e=0.25), and
  prolactin-sparing agents (large `KI_AP`).
- Surgery is implemented as rapid first-order removal over a 0.5-day window,
  and surgical hypopituitarism, unlike compressive hypopituitarism, is added
  as a **permanent term**.

---

## 13. How to run

```bash
# Render the mechanistic map
dot -Tsvg prl_qsp_model.dot -o prl_qsp_model.svg
dot -Tpng -Gdpi=150 prl_qsp_model.dot -o prl_qsp_model.png

# 34 scenarios + 15 diagnostic analyses (about 10-15 minutes)
Rscript prl_mrgsolve_model.R

# Interactive dashboard (11 tabs)
Rscript -e 'shiny::runApp("prl_shiny_app_en.R", port = 8080)'
```

Required packages: `mrgsolve` (≥2.0), `shiny`, `ggplot2`, `dplyr`, `tidyr`, `DT`,
and Graphviz.

---

## 14. References

[`prl_references_en.md`](prl_references_en.md) — **150
references, every PMID looked up via NCBI E-utilities to confirm title,
journal, and year.** None was cited from memory. Each section states which
equation, parameter, or diagnostic analysis in the model that group of
references supports, and it also records candidates that failed the relevance
check and **parameters for which no verifiable citation could be found to
anchor them**.

## ⚠️ Disclaimer

This is a qualitative/semi-quantitative QSP model for educational and research
purposes. It has not been independently validated or certified and must not be
used directly for clinical decision-making, prescribing, or regulatory
submission.
