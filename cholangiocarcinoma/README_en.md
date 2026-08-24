# Cholangiocarcinoma (CCA) QSP Model

**Biliary tract cancer · intrahepatic / perihilar / distal cholangiocarcinoma**

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (230 nodes · 22 clusters) | [`cca_qsp_model.dot`](cca_qsp_model.dot) · [SVG](cca_qsp_model.svg) · [PNG](cca_qsp_model.png) |
| ⚙️ mrgsolve model (54 ODEs) | [`cca_mrgsolve_model.R`](cca_mrgsolve_model.R) |
| 🐍 Dependency-free reference twin | [`cca_reference_impl.py`](cca_reference_impl.py) |
| 📊 Shiny dashboard (13 tabs) | [`cca_shiny_app_en.R`](cca_shiny_app_en.R) |
| 📚 References (87 articles, PMID-verified) | [`cca_references.md`](cca_references.md) |

```bash
python3 cca_reference_impl.py --checks   # 32/32 structural self-checks
python3 cca_reference_impl.py 120        # 15 virtual-trial arms
dot -Tsvg cca_qsp_model.dot -o cca_qsp_model.svg
```

---

## Why this model exists

Cholangiocarcinoma is usually drawn like this — **risk factor → bile-duct epithelial mutation → tumour growth
→ chemotherapy shrinks it → the patient lives longer (or does not).** This picture cannot withstand the
following seven observations at once.

| | Observation | What a single-line cascade fails to explain |
|---|---|---|
| **A** | ABC-02 gemcitabine/cisplatin: response rate 26%, mPFS 8.0 months, mOS 11.7 months | Three-quarters of the patients who benefit never once reach a RECIST response. Response and survival are only loosely coupled |
| **B** | The most common terminal event is not the mass but the bile duct | Cholangitis and obstructive liver failure arrive well before tumour volume alone would be fatal |
| **C** | Biliary drainage has **no** antitumour action whatsoever | Yet without drainage no systemic therapy can be delivered at full dose (bilirubin holds both gemcitabine and cisplatin) |
| **D** | FGFR2-fusion tumours progress on pemigatinib through **polyclonal** kinase-domain mutations (V564F, N550K, …), and ctDNA moves before imaging does | If resistance were "induced," there would be no reason for it to be polyclonal |
| **E** | TOPAZ-1: median survival 11.5 → 12.8 months (barely moves), 24-month survival 10.4% → 24.9% (a large gap) | A model that tracks average tumour size **cannot open up a tail like that** |
| **F** | ClarIDHy ivosidenib: response rate ~2% yet a PFS hazard ratio of 0.37 | Benefit without shrinkage — cytostatic, not cytotoxic |
| **G** | FGFR2 fusions and IDH1 mutations occur essentially only in intrahepatic cholangiocarcinoma, while perihilar lesions fail through obstruction | Anatomical location itself changes the mode of failure |

So instead of a cascade, this model takes **three structural positions**.

---

### I. Delivered dose is not an input but an **output**

Every clinical trial protocol writes "gemcitabine 1000 mg/m², d1, d8," and every clinical trial reports a
relative dose intensity well below 100%. In this model, the prescription enters only as *intent*, and it is
multiplied by **four gates computed from the patient's state at the moment of dosing**.

```
Rule 1  Bilirubin       GATE = 1 / (1 + (BILI/2.6)^6)
Rule 2  Neutrophils     GATE = 1 / (1 + (1.0/ANC)^8)
Rule 3  Renal function  GATE = 1 / (1 + (45/CrCl)^8)     ← cisplatin only
Rule 4  Performance     GATE = 1 / (1 + (PS/2.6)^8)
```

Because a perihilar tumour causes obstruction, obstruction raises bilirubin, and bilirubin closes Rule 1, this
model contains a **closed positive-feedback loop** —
`growth → obstruction → dose withholding → growth`. On the map (`cca_qsp_model.dot`) this loop is literally
drawn as a bold red line labelled `LOOP 1/4 … 4/4`.

**Drainage kills no tumour cell anywhere in this file.** Yet it is the node with the largest leverage on the
map — because it is the only edge that breaks that loop. That is the point.

> **The falsifier is a single parameter.** With `GATE_ON = 0` all four rules disappear.
> Delaying drainage by 60 days should then have no effect on survival.

### II. Resistance is not induced, it is **selected**

The FGFR2 kinase-domain mutant clone `T_R` is already seeded at t = 0 by the Goldie–Coldman rule
`T_R(0) = T(0) · μ · ln(N)`. Nowhere in this file does a sensitive cell *convert* into a mutant cell. There is
no "time of resistance emergence" parameter.

**This position produced a negative result, and it is reported rather than hidden.** At a defensible seed
frequency (about 2.5 × 10⁻⁵), this clone **cannot** achieve clinical dominance within a 7-month PFS. The
arithmetic forbids it — `ln(0.2 / 2.5e-5) / 0.026` comes out to about 11 months. So the model assigns the
median PFS on an FGFR inhibitor to **drug-tolerant persister cells (`T_P`)**, and leaves the clone responsible
only for the events it can actually explain — loss of an established response after about a year, ctDNA rises
that precede imaging, and the cross-resistance asymmetry between the reversible and covalent inhibitors
(`RHO_PEM` 0.05 vs `RHO_FUT` 0.45). Consequently `MU_RES = 0` changes almost nothing within the first year,
and **the falsifier table below records exactly that.**

### III. Survival is **two competing risks** running on different clocks

```
h_tumour   slow clock : tumour volume/liver volume + hepatic reserve deficit + cachexia
h_biliary  fast clock : cholangitis burden + excess ALBI,  recurrent and largely preventable
S(t) = exp(−∫(h_tumour + h_biliary))
```

Population survival is the **average** of individual `S(t)` curves. That is how a single simulation produces
both a median and a 24-month tail at once. And it is exactly why immunotherapy — because it acts not on
everyone a little but only on the immune-engaged tumour fraction `PI_IMMUNE` — **barely moves the median while
opening up the tail.**

---

## What was fitted and what was not

Only 8 parameters were fitted to data: `K_GEM`, `K_CIS` (cytotoxic potency), `KSP0`, `ALPHA_P` (persister-cell
entry and residual killing), `KFGFR_KILL`, `KIMM`, and two hazard slopes, `HT_T` and `HB_CH`. Everything else
is a published PK parameter, a labelled dose, a physiological constant, or a structural choice fixed in
advance (Simeoni `PSI = 20`, 3 Friberg transit compartments, the Goldie–Coldman seeding rule).

---

## Virtual-Trial Results (n = 60/arm, `python3 cca_reference_impl.py 60`)

OS and PFS in months. **Bold = predicted, not used in fitting.**

| Scenario | mOS | mPFS | ORR | 12 mo | 24 mo | RDI | Actual trial |
|---|---|---|---|---|---|---|---|
| S1 Best supportive care + drainage | 5.8 | 1.7 | 0% | 6.6% | 0% | — | Historical BSC 4–6 months |
| S2 Gemcitabine/cisplatin | 11.0 | 7.4 | 31.7% | 45.5% | 4.0% | 94% | ABC-02 11.7 / 8.0 / 26.1% |
| S3 + Durvalumab | **13.1** | 7.6 | 45.0% | 54.0% | **20.3%** | 94% | TOPAZ-1 12.8 / — / 26.7% / 24-mo 24.9% |
| S4 Drainage delayed 60 days, no exchange | **10.2** | 6.0 | 28.3% | 43.7% | 18.4% | **68%** | (prediction) |
| S5 Plastic stent | **10.9** | 7.4 | 30.0% | 45.6% | 4.0% | **89%** | (prediction) |
| S6 Pemigatinib 14/7 | 10.9 | 6.2 | **40.0%** | 43.1% | 0% | — | FIGHT-202 21.1 / 6.9 / 35.5% |
| S7 Futibatinib continuous | 11.8 | 6.4 | 80.0% | 49.2% | 0.1% | — | FOENIX-CCA2 20.0 / 9.0 / 42% |
| S8 Ivosidenib | 8.1 | 4.4 | **0%** | 25.1% | 0% | — | ClarIDHy 10.3 / 2.7 / **2%** |
| S11 Pemigatinib continuous (prediction) | **11.6** | **6.4** | 65.0% | 47.6% | 0.1% | — | Never run |
| S9 R0 resection + adjuvant capecitabine | 31.2 | 19.5 (RFS) | — | 92.0% | 84.1% | — | BILCAP 53 / RFS 24.4 |
| S10 R0 resection, observation | 27.0 | 14.8 (RFS) | — | 91.9% | 72.9% | — | BILCAP 36 / RFS 17.5 |

### What came out without being fitted (predictions)

- **The TOPAZ-1 signature is reproduced.** Median survival versus control rises from 11.0 to 13.1 months
  (+2.1), and 24-month survival from 4.0% to **20.3%**. The actual trial found 11.5 → 12.8 (+1.3), 10.4% →
  24.9%. This pattern — a median that barely moves alongside a tail that opens up — **is a consequence of the
  mixture parameter, not a fitted value.** It cannot come out of an average-effect model in which everyone
  responds a little.
- **Delaying drainage by 60 days drops delivered dose from 94% to 68% and shortens median survival from 13.1
  to 10.2 months.** **Not one** antitumour parameter changed — only the timing of stent placement did.
- **Using a plastic instead of a metal stent** shortens the patency half-life from 239 to 90 days, dropping
  RDI from 94% to 89% and median survival from 11.0 to 10.9 months. Small, but in the right direction and
  arising purely from patency.
- **Ivosidenib comes out with a 0% response rate.** It only lowers the growth rate and is attached to a node
  with zero kill rate, which is exactly the shape of ClarIDHy's 2%.
- **S11 is an experiment that has never been run.** Removing pemigatinib's 1-week drug holiday raises the
  response rate from 40% to 65% and median survival from 10.9 to 11.6 months, because `T_S` and `T_P` both
  regrow at full rate during every off week — a testable prediction.

### Falsifiers (one parameter changed at a time, nothing else refitted)

| | Change | Prediction | Result |
|---|---|---|---|
| **F1** | `GATE_ON = 0` | The survival loss from delayed drainage should disappear | **Passes.** S4 goes from 10.2 to 13.1 months, RDI from 68% to 100%. The entire loss came from the gates |
| **F2** | `MU_RES = 0` | PFS on an FGFR inhibitor should lengthen | **Negative result — reported.** Identical to S6 to the decimal place (10.9 / 6.2 / 40.0%). The seeded clone does nothing at all within the observation window. See II above |
| **F3** | `PI_IMMUNE = 1` | Durvalumab should move the **median** | **Passes.** Median survival is pushed beyond the observation horizon, with 24-month survival at 76.9%. Since TOPAZ-1 says the median barely moves, the assumption that every tumour is immune-engaged is ruled out |
| **F4** | `FPEN_MIN = 1` | The gem/cis response rate should exceed the target | **Passes.** 31.7% → **86.7%**, median survival 11.0 → 17.3 months. Without the stromal penetration barrier, cholangiocarcinoma becomes a chemosensitive tumour |

### What misses (recorded without being asked)

1. **Response rates run high overall.** gem/cis 31.7% vs. actual 26.1%; futibatinib 80% vs. 42%. Too large a
   fraction of `T_S` dies before converting into a persister — especially with continuous oral dosing, because
   `e_sp` stays lower than with pulsed intravenous dosing.
2. **Survival on FGFR inhibitors is substantially underestimated.** S6 10.9 months vs. FIGHT-202's 21.1; S7
   11.8 vs. FOENIX-CCA2's 20.0. Patients in both of those trials were already a selected population enrolled
   with good performance status after first-line therapy, and the model does not represent that selection.
   `HT_T` being too steep for intrahepatic cholangiocarcinoma also contributes.
3. **Survival after adjuvant therapy is underestimated.** 31.2 / 27.0 months vs. BILCAP's 53 / 36. Recurrence-
   free survival (19.5 / 14.8 vs. 24.4 / 17.5) is much closer, meaning the post-recurrence trajectory is too
   fast — because a micro-lesion at the time of recurrence is built to grow at the same `LAM0` as the tumour
   at diagnosis.
4. **Misses 1 and 2 are likely the same defect.** Both come from the assumption that "the lesion remaining
   after treatment follows the same dynamics as the original tumour." Giving `T_P` its own, separate slow
   regrowth clock looks likely to improve both together, and that will be the first thing tried in the next
   revision.

---

## Model Structure

### Mechanistic map — 230 nodes, 22 clusters

0 falsification-target observations · 1 pathogenesis and chronic biliary injury · 2 cell of origin and
anatomical subtype · 3 genomic drivers · 4 the IDH1–2HG axis · 5 FGFR2 signalling and the phosphate axis ·
6 desmoplastic stroma and the drug-penetration barrier · 7 the tumour immune microenvironment ·
8 tumour cell populations (selected resistance) · 9 biliary obstruction · drainage · cholangitis ·
10 hepatic reserve · 11 gemcitabine PK/PD · 12 cisplatin PK/PD · 13 FGFR inhibitors · 14 IDH1 inhibitors ·
15 immune checkpoint inhibition · 16 other systemic therapy · 17 surgery and local therapy ·
18 **dose gates** · 19 biomarkers · 20 clinical endpoints and the two risks · 21 falsifiers

### 54 ODEs

| Compartment group | Count | Content |
|---|---|---|
| PK | 13 | Gemcitabine 2 compartments + intracellular dFdCTP · free platinum 2 compartments + Pt-DNA adducts · durvalumab 2 compartments · FGFR inhibitor absorption/central/covalent occupancy · ivosidenib · 5-FU |
| Tumour | 8 | `T_S` sensitive · `T_P` persister · `T_R` FGFR2-mutant clone · 3 damaged-cell transit compartments · metastatic burden · RECIST nadir tracker |
| Microenvironment | 5 | CAF · extracellular matrix · CD8 effector cells · suppressive compartment · IL-6 |
| Biliary/hepatic | 7 | Obstruction fraction · stent patency · bilirubin · ALP · albumin · functional hepatic reserve · cholangitis |
| Host/toxicity | 10 | 5 Friberg myelosuppression compartments · platelets · creatinine clearance · neuropathy · lean body mass · performance status |
| Biomarkers | 5 | CA 19-9 (confounded by cholestasis) · ctDNA · 2-HG · phosphate · immune-related adverse events |
| Survival/tracking | 6 | Cumulative tumour hazard · cumulative biliary hazard · S(t) · cumulative cisplatin · delivered dose · prescribed dose |

### Differences between the R model and its Python twin

`cca_reference_impl.py` is an exact twin of `cca_mrgsolve_model.R` — the same 54 states, the same names, the
same numbers. Every value in the table above came from the Python file, which has **no dependencies at all**
(neither numpy nor scipy). The two files differ at exactly two mechanical points.

- **Stent exchange.** Python re-stents when patency falls below 0.25 (a state-triggered event). mrgsolve has
  no state-triggered events, so the R driver schedules a fixed-interval exchange instead.
- **Integrator.** Because gemcitabine's plasma half-life is about 6 minutes, Python uses a fixed-step RK4 that
  shrinks the step to 0.002–0.01 day around infusions. mrgsolve uses LSODA and needs no such adjustment.

---

## Shiny Dashboard (13 tabs)

1 patient summary · 2 PK · **3 clones** · **4 dose gates** · 5 biliary obstruction · drainage ·
6 hepatic reserve · ALBI · 7 haematological toxicity · 8 immune · 9 biomarkers ·
10 RECIST · clinical endpoints · **11 the two risks** · 12 scenario comparison · **13 falsification tests**

On tab 4, moving only the drainage-delay slider shows delivered RDI fall without touching any antitumour
parameter. Switching on F1 should make that effect disappear — that is how to falsify this model.

---

## ⚠️ Disclaimer

This is a qualitative/semi-quantitative QSP model for education and research purposes. It was built from
public literature and clinical trial results but has not been independently validated or certified, and
**must not be used directly for clinical decision-making, prescribing, or regulatory submission.** Parameters
are illustrative approximations, and fitting and validation against real patient data would be needed
separately. The "What misses" section above is not a complete list.
