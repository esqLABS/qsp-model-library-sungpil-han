# Langerhans Cell Histiocytosis (LCH) — QSP Model

<p align="center">
  <a href="lch_qsp_model_en.svg"><img src="lch_qsp_model_en.png" width="720" alt="LCH QSP mechanistic map"></a>
</p>

A quantitative systems pharmacology (QSP) model of LCH — an inflammatory
myeloid neoplasm that begins when a bone-marrow dendritic-cell precursor
carrying a MAPK-pathway mutation settles in the wrong place. It combines,
into a single system of differential equations, the mutant precursor
reservoir, the ERK-driven senescence secretome, bone destruction,
irreversible hypothalamic-pituitary and cerebellar damage, and the PK/PD
of vemurafenib, dabrafenib+trametinib, vinblastine/prednisolone, and
cladribine+cytarabine.

| File | Contents |
|---|---|
| [`lch_qsp_model_en.dot`](lch_qsp_model_en.dot) | mechanistic map source (20 clusters, 213 nodes, 309 edges) |
| [`lch_qsp_model_en.svg`](lch_qsp_model_en.svg) / [`lch_qsp_model_en.png`](lch_qsp_model_en.png) | rendered map (Graphviz `dot`, PNG 150 dpi) |
| [`lch_mrgsolve_model.R`](lch_mrgsolve_model.R) | mrgsolve model — **62 ODEs**, 5 phenotypes, 12 scenarios |
| [`lch_python_twin.py`](lch_python_twin.py) | a dependency-free Python twin (RK4) + **51 quantitative checks passed** |
| [`lch_shiny_app_en.R`](lch_shiny_app_en.R) | Shiny dashboard — 11 tabs (including a falsification panel) |
| [`lch_references_en.md`](lch_references_en.md) | 136 references, all confirmed on PubMed |

```bash
python3 lch_python_twin.py     # 51/51 checks passed
Rscript -e 'shiny::runApp("lch_shiny_app_en.R")'
```

---

## What this model argues (Three structural commitments)

The textbook LCH model is "a clone proliferates and the drug kills it."
That single picture cannot satisfy six facts the literature demands at the
same time.

| # | Fact the literature demands |
|---|---|
| i | The same BRAF V600E produces a spontaneously resolving skull lesion in one child and fatal risk-organ disease in another |
| ii | Lesional Ki-67 is only 5–10%, yet the clinical picture is an inflammatory storm ([PMID 9588881](https://pubmed.ncbi.nlm.nih.gov/9588881/)) |
| iii | MAPK inhibitors produce a clinical response **within days**, while blood BRAF V600E cfDNA rarely disappears ([PMID 30718231](https://pubmed.ncbi.nlm.nih.gov/30718231/), [PMID 34383272](https://pubmed.ncbi.nlm.nih.gov/34383272/)) |
| iv | More than 75% reactivate within weeks to months of discontinuation ([PMID 31513482](https://pubmed.ncbi.nlm.nih.gov/31513482/), [PMID 28667012](https://pubmed.ncbi.nlm.nih.gov/28667012/)) |
| v | Cladribine works even though the target cells barely divide |
| vi | Central diabetes insipidus and neurodegeneration are irreversible, and correlate with **how long** disease was active rather than with **how severe** it was ([PMID 16047354](https://pubmed.ncbi.nlm.nih.gov/16047354/), [PMID 15049016](https://pubmed.ncbi.nlm.nih.gov/15049016/)) |

So three structural choices were made, each fitted with a **kill switch**
to keep it falsifiable.

### ① Phenotype is set by **cell-of-origin partitioning**, not growth rate

There is a single mutant precursor pool. What distinguishes the
phenotypes is entirely a set of **origin descriptors** — `FSR` (the
self-renewal capacity of the differentiation stage at which the driver
arose), `PRECM0` (its size), and the seeding partition
`THB/THS/THR/THP/THC/THL`.

**Not one of the lesion dynamics constants (`KPROL`, `KDL`, `KIMM`,
`KDP`, `KSR`, `LMAX`, `KSEED`, `KDC`) differs between phenotypes.**

| Phenotype | FSR | PRECM0 | θ(risk organ) | θ(pituitary) |
|---|---|---|---|---|
| SS-b, single-system bone | 0.10 | 0.020 | 0.04 | 0.06 |
| MS RO-negative | 0.72 | 0.12 | 0.06 | 0.14 |
| MS RO-positive | 1.00 | 0.45 | 0.34 | 0.08 |
| CNS-risk + pituitary | 0.85 | 0.25 | 0.10 | 0.30 |
| Adult pulmonary LCH | 0.24 | 0.06 | 0.02 | 0.02 |

This structure is a direct translation into equations of Berres 2014
([PMID 25646268](https://pubmed.ncbi.nlm.nih.gov/25646268/)) and
Xiao 2020's ([PMID 32750121](https://pubmed.ncbi.nlm.nih.gov/32750121/))
claim that "the differentiation stage at which the mutation arose
determines disease extent."

**Falsification:** giving the SS-b phenotype only the MS RO+ origin
descriptors — without touching a single rate constant — switches on
risk-organ burden, from 0.000 to 0.051.

### ② Lesions are sustained by **recruitment**, and MAPK inhibition is **cytostatic**

`KPROL = 0.022/day` was **deliberately** set below the minimum immune
clearance rate, `KDL + KIMM·(1−TREGMAX) = 0.0288/day`. So **no lesion can
sustain itself by local division alone**; it must keep being resupplied
from the precursor pool. This is the quantitative content of the low
Ki-67.

MAPK inhibitors lower only `CCND` (proliferation), `BCL` (survival
signal), and `SASP` (the secretory programme) — **they carry no killing
term at all** (`SL_MAPKI_KILL = 0`). Only nucleoside analogues,
vinblastine, and glucocorticoids carry a killing term. So the following
falls out **as a result, not an assumption**:

```
KILL_LES = E_ARAC·PFN + E_VBL·PFN + E_CLAD + EGR_KILL + SL_MAPKI_KILL·(IBRAF+IMEK)
                    └── S-phase-dependent   └── acts even on non-dividing cells (2-CdA)   └── fixed at 0
```

| Validated result | Value |
|---|---|
| week-1 DAS drop vs lesion-mass drop | −30% vs −17% (secretome falls first) |
| clinical-response/mass-reduction ratio (MAPKi vs 2-CdA/Ara-C) | 1.65 vs 0.80 |
| cfDNA at 12 months on MAPKi | 0.282 → 0.014 (**above** the 0.005 detection limit = still detectable) |
| reservoir at 12 months on MAPKi | 92% of baseline remains |
| after the same duration on 2-CdA/Ara-C | cfDNA below detection limit, reservoir extinct |
| DAS after discontinuation | 0.08 → 8.75 (relapse) |

Two auxiliary mechanisms create the plateau:

- **`FN_APO = 0.25`** — the bone-marrow niche supplies ERK-independent
  survival signals, so loss of BCL2A1 drives lesional cells into apoptosis
  far more strongly than reservoir cells. The reservoir does not
  disappear; it forms a **plateau.**
- **`PNICHE = 0.003`** — the minimum clone size at which self-renewal
  collapses (a deterministic stand-in for stochastic clonal extinction).
  MAPK inhibition never reaches it; **6 cycles** of 2-CdA/Ara-C do, but
  **3 cycles** do not.

**Prediction that falls out of this:** dose intensity has a **threshold**,
not a slope — at day 500 the reservoir is 2.2531 after 3 cycles and
0.000000 after 6.

**Falsification:** setting `SL_MAPKI_KILL = 0.25` abolishes relapse after
discontinuation (ratio 108.7 → 0.0) — directly contradicting literature
fact (iv).

### ③ Permanent sequelae are a function not of peak severity but of **the time-integral of active disease**

`AVPN` (AVP-secreting neurons), `ANTPIT`, `NEUR` (cerebellar/pontine
neurons), `BILF`, and `LUNGC` all move **monotonically** (decreasing or
increasing) with no recovery term. Central diabetes insipidus is declared
at `AVPN < 0.15`, clinical neurodegeneration at `NEUR < 0.70`.

```
dAVPN/dt = −KAVP·(LPIT/LREF)·√(IL1B/IL1B0)·AVPN        (cannot become positive)
dNEUR/dt = −KND·[(LCNS/LREF) + W_CDI_ND·1{CDI}·(LCNS/LREF + 0.15)]·NEUR
```

Two scenarios for the same patient, differing only in when treatment
starts (S9 vs S9b):

| | treated at day 14 | delayed 180 days |
|---|---|---|
| AVP neuron pool (d730) | 0.872 | **0.001** |
| Central diabetes insipidus | none | **permanent** |
| Neuron pool | 0.962 | 0.211 |
| Days of active disease `TTET` | 14 days | 185 days |

**And the model computes the crossover point.** Does "slow-but-deep"
cytotoxic therapy lose to "fast-but-shallow" cytostatic therapy on the CNS
endpoint? At current calibration, **there is no crossover if both start
at diagnosis** (deep arm AVPN 0.929 vs fast arm 0.878). The crossover
first appears once cytotoxic therapy is **delayed by 15 days.** The model
computes this boundary; it does not assume it.

---

## Model structure (62 ODEs)

| Block | # compartments | Compartments |
|---|---|---|
| Drug PK | 22 | vemurafenib (2-compartment + auto-induction) · dabrafenib (+ active metabolite) · trametinib (2-compartment) · cobimetinib · vinblastine (2-compartment) · prednisolone · cytarabine (+ intracellular ara-CTP) · cladribine (+ intracellular Cd-ATP) · 6-MP/MTX maintenance |
| Signalling | 5 | `ERK` `CCND` `BCL` `SASP` `GRE` |
| Cell pools | 10 | `PRECM` (marrow reservoir) `CIRC` `LBONE` `LSKIN` `LRO` `LPIT` `LCNS` `LLUNG` `TREG` `OCL` |
| Secretome | 6 | `IL1B` `TNFA` `IL6` `OSM` `MMP9` `RANKL` |
| Biomarkers | 4 | `CFDNA` `SCD163` `CRP` `FERR` |
| Irreversible-organ pools | 6 | `AVPN` `ANTPIT` `NEUR` `BILF` `LUNGC` `BVOL` |
| Myelosuppression (Friberg) | 4 | `PROL` `TR1` `TR2` `ANC` |
| Targeted-therapy toxicity | 2 | `SKTOX` (paradoxical activation) `LVEF` |
| Cumulative endpoints | 3 | `AUCERK` `TTET` `CUMDAS` |

Time is in **days**, concentrations in mg/L, lesion burden in burden units
(1 unit ≈ 10⁹ LCH cells). All PK parameters are referenced to 70 kg and
allometrically scaled to body weight (CL^0.75, V^1.0) — without this, a
12 kg infant's steady-state vemurafenib at 20 mg/kg/day is predicted at a
quarter of the observed range (20–60 mg/L).

### The genotype switch and paradoxical activation

At `GENO = 2` (MAP2K1/ARAF/driver-negative), the BRAF-inhibition term goes
to zero, and instead **paradoxical MAPK activation** switches on via CRAF
dimer transactivation. MEK inhibitors sit downstream and suppress both
cases.

| Validated result | Value |
|---|---|
| BRAF inhibitor in MAP2K1-driven disease | DAS 3.54 → 8.42 (no response, progression) |
| MEK inhibitor in the same genotype | d90 DAS 4.22 vs 8.42 |
| paradoxical skin toxicity when a MEK inhibitor is added | `SKTOX` 1.83 → 0.54 |

---

## 12 scenarios

| # | Scenario | Point |
|---|---|---|
| S1 | SS-b, observation | a self-renewal-poor origin → spontaneous resolution, re-ossification after bone lesions form |
| S2 | MS RO− LCH-III, vinblastine/prednisolone, 12 months | standard first-line therapy |
| S3 | MS RO+, 6-week non-response → 2-CdA/Ara-C rescue | how the 6-week response splits prognosis |
| S4 | MS RO+, 2-CdA/Ara-C 6 cycles, front-line | reservoir extinction → durable remission |
| S5 | vemurafenib, continuous | fast response, cfDNA plateau |
| S6 | vemurafenib, **discontinued** at 12 months | relapse + new permanent damage |
| S7 | dabrafenib + trametinib, continuous | deeper pERK suppression, paradoxical suppression |
| S8 | MAPKi 8-week bridge → consolidation with 2-CdA/Ara-C → full discontinuation | buying time with cytostasis, locking it in with killing |
| S9 / S9b | CNS-risk lesion, treatment delayed 180 days vs 14 days | direct test of structural choice ③ |
| S10 / S10b | adult pulmonary LCH, continued smoking vs smoking cessation at 3 months | nodules heal, cystic change persists |

Summary at 2 years (3 years for pulmonary LCH):

| Scenario | DAS | cfDNA | AVPN | NEUR | CDI | ND | TTET (days) |
|---|---|---|---|---|---|---|---|
| S1 observation | 0.00 | 0.00 | 0.999 | 1.000 | 0 | 0 | 0 |
| S2 LCH-III | 0.00 | 0.00 | 0.952 | 0.985 | 0 | 0 | 0 |
| S3 first-line failure → rescue | 0.00 | 0.00 | 0.987 | 0.993 | 0 | 0 | 0 |
| S4 2-CdA/Ara-C front-line | 0.00 | 0.00 | 0.986 | 0.988 | 0 | 0 | 1 |
| S5 vemurafenib, continuous | 0.02 | 0.00 | 0.969 | 0.981 | 0 | 0 | 4 |
| **S6 vemurafenib discontinued** | **8.75** | **0.93** | **0.004** | **0.270** | **1** | **1** | **316** |
| S7 dabrafenib+trametinib | 0.00 | 0.00 | 0.979 | 0.988 | 0 | 0 | 2 |
| S8 bridge→consolidation | 0.00 | 0.00 | 0.976 | 0.986 | 0 | 0 | 4 |
| **S9 diagnosis delayed 6 months** | 0.00 | 0.00 | **0.001** | **0.211** | **1** | **1** | **185** |
| S9b treated at 14 days | 0.00 | 0.00 | 0.872 | 0.962 | 0 | 0 | 14 |
| S10 pulmonary LCH, smoking cessation | 1.23 | 0.05 | 0.963 | 0.975 | 0 | 0 | 140 |
| S10b pulmonary LCH, continued smoking | 4.52 | 0.26 | 0.819 | 0.923 | 0 | 0 | 1083 |

The contrast between S6 and S9 is this model's point — **controlling the
disease is not enough on its own; when it was controlled determines the
damage left behind.** In S6, relapse after discontinuation produces fresh
permanent damage (AVPN 0.004, ND occurs).

---

## Verification

`lch_python_twin.py` reimplements **the same right-hand side and the same
parameter blocks** as the mrgsolve model in pure Python (RK4, no numpy
required), and runs 51 quantitative checks.

```
51 / 51 checks passed
```

Included checks:

- **PK plausibility** — paediatric vemurafenib Css 22.3 mg/L, trametinib
  trough 9.4 ng/mL, exposure reduction from auto-induction, combination
  pERK of 15.3% vs 18.6% for monotherapy
- **Natural history** — SS-b resolves spontaneously and re-ossifies with
  no diabetes insipidus; MS RO+ shows rising risk-organ burden and
  develops diabetes insipidus and biliary fibrosis
- **All of structural choices ①②③, each with its own kill switch**
- **Pharmacodynamic distinction** — cladribine killing is not
  cell-cycle-gated; limited CNS penetration makes CNS lesions respond less
  than systemic lesions (at d30, 0.62 of baseline remaining in the CNS vs
  0.10 in risk organs)
- **Toxicity** — grade-4 neutropenia on 2-CdA/Ara-C (nadir 0.05 × 10⁹/L)
  with recovery before the next cycle; vinblastine/prednisolone much
  milder (nadir 2.52); MEK-inhibitor LVEF nadir 54.8%
- **Numerical hygiene** — no negative states over a 2-year simulation;
  halving the integration step changes DAS by < 10⁻⁷

## What this model does not claim (Explicit limitations)

- At current calibration, **cytotoxic therapy is not slower than MAPK
  inhibitors** (5.5 days vs 21.0 days to DAS control). So starting both at
  the same time produces no crossover on the CNS endpoint. The crossover
  appears only once cytotoxic therapy is delayed by 15 days or more. The
  model **computes** this boundary; it does not assume it.
- The IL-17A/dendritic-cell fusion pathway, switchable via `IL17ON`, is
  **contested** in the literature
  ([PMID 18157139](https://pubmed.ncbi.nlm.nih.gov/18157139/) vs
  counter-argument). It defaults to on, but switching it off shows its
  contribution is limited to a drop in osteoclast output.
- `PNICHE` (the clonal-extinction threshold) is a deterministic stand-in
  for stochastic extinction. In reality it is a **probability** of
  relapse depending on residual clone size — a distribution, not a
  threshold.
- The lesion-burden unit is an arbitrary scale. cfDNA is likewise an
  internal model signal, not an actual VAF (%); interpret it clinically
  only as a relative value against the detection limit (`CF_LOD`).
- Immune clearance (`KIMM`) and Treg expansion are driven by a single term
  based on total systemic burden. Organ-specific local immune
  microenvironment differences are not represented.
- Radiotherapy, local steroid injection, denosumab, PD-1 blockade, and
  HSCT appear on the map but are not implemented as ODEs.

---

## References

All 136 references were confirmed on NCBI E-utilities for PMID existence
and title match — see
[`lch_references_en.md`](lch_references_en.md).

Key evidence:

- **Clonal origin** Willman 1994 ([PMID 8008029](https://pubmed.ncbi.nlm.nih.gov/8008029/)) · Badalian-Very 2010 ([PMID 20519626](https://pubmed.ncbi.nlm.nih.gov/20519626/))
- **Differentiation stage of origin** Berres 2014 ([PMID 25646268](https://pubmed.ncbi.nlm.nih.gov/25646268/)) · Xiao 2020 ([PMID 32750121](https://pubmed.ncbi.nlm.nih.gov/32750121/))
- **Senescence secretome** Bigenwald 2021 ([PMID 33958797](https://pubmed.ncbi.nlm.nih.gov/33958797/))
- **Intralesional trapping** Hogstad 2018 ([PMID 29263218](https://pubmed.ncbi.nlm.nih.gov/29263218/))
- **Low proliferative fraction** Brabencova 1998 ([PMID 9588881](https://pubmed.ncbi.nlm.nih.gov/9588881/))
- **80% pERK inhibition threshold** Bollag 2010 ([PMID 20823850](https://pubmed.ncbi.nlm.nih.gov/20823850/))
- **MAPKi response and residual clones** Donadieu 2019 ([PMID 31513482](https://pubmed.ncbi.nlm.nih.gov/31513482/)) · Eckstein 2019 ([PMID 30718231](https://pubmed.ncbi.nlm.nih.gov/30718231/))
- **Relapse after discontinuation** Cohen Aubart 2017 ([PMID 28667012](https://pubmed.ncbi.nlm.nih.gov/28667012/))
- **Rescue therapy** Donadieu 2015 ([PMID 26194764](https://pubmed.ncbi.nlm.nih.gov/26194764/))
- **Permanent sequelae** Grois 2006 ([PMID 16047354](https://pubmed.ncbi.nlm.nih.gov/16047354/)) · Haupt 2004 ([PMID 15049016](https://pubmed.ncbi.nlm.nih.gov/15049016/))
- **Neurodegeneration mechanism** Wilk 2023 ([PMID 38091952](https://pubmed.ncbi.nlm.nih.gov/38091952/)) · Vicario 2025 ([PMID 40081365](https://pubmed.ncbi.nlm.nih.gov/40081365/))

---

## Disclaimer

A computational model for educational and research purposes. It is not a
clinical practice guideline and must not be used for individual patient
treatment decisions. Parameters are representative values derived from
published literature and were not calibrated against individual patient
data.
