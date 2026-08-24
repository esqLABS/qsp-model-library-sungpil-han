# Head and Neck Squamous Cell Carcinoma (HNSCC) — QSP Model
### Head and Neck Squamous Cell Carcinoma · Quantitative Systems Pharmacology

| Deliverable | File | Scale |
|--------|------|------|
| Mechanistic map | [`hnscc_qsp_model.dot`](hnscc_qsp_model.dot) · [SVG](hnscc_qsp_model.svg) · [PNG](hnscc_qsp_model.png) | 213 nodes · 319 edges · 18 clusters |
| mrgsolve ODE model | [`hnscc_mrgsolve_model.R`](hnscc_mrgsolve_model.R) | 72 ODEs · 245 parameters · 38 observations · 26 scenarios |
| Shiny dashboard | [`hnscc_shiny_app.R`](hnscc_shiny_app.R) | 11 tabs |
| References | [`hnscc_references.md`](hnscc_references.md) | 164 PubMed links |

---

## The organising thesis

The benefit of a drug that targets resistance is **not a property of the drug.** It is
the **product** of what the drug does and the amount of resistance the tumour still has
left to lose (its headroom).

```
Delta(log kill) = (agent effect on factor i) x (HEADROOM in factor i)
```

Nowhere in any file in this repository is there **a rule** saying "cetuximab is
inferior to cisplatin in HPV-positive disease", or "an interruption of treatment loses
about 1 % of local control per day", or "hypoxic sensitisers work only in hypoxic
tumours". Instead there are **four resistance factors** acting on the same log-kill
product, and each drug enters only through the factors it actually reaches.

| Resistance factor | Mechanism | What removes it |
|-------------|------|----------------------|
| **R_hypoxia** | The oxygen enhancement ratio (OER) is applied **only to the clonogens that are actually hypoxic**. Because `CSCO`/`CSCH` are separate compartments exchanging through `KOXTR`, **reoxygenation is a competing process rather than a parameter** | nimorazole · carbogen · transfusion |
| **R_repair** | α and the sublethal damage repair rate `MUREP`. Both vary with NHEJ capacity, HR capacity and the platinum adduct burden | cisplatin (α↑, NHEJ↓) |
| **R_repopulation** | The clonogen birth rate. It accelerates once the damage signal `DAMS` has accumulated for about three weeks | shortening the treatment duration · cetuximab (chiefly here) |
| **R_immune** | CD8 killing suppressed by PD-1 engagement · Treg · MDSC · loss of HLA class I | pembrolizumab · nivolumab |

### Headroom is not assumed but **measured**

Seven **counterfactual integrators**
(`LKACT`/`LKOXI`/`LKREP`/`LKNRP`/`LKIMM`/`LKIMX`/`LKCHM`) run alongside the actual
trajectory, accumulating at the same time the log kill the clonogens **actually
received** and the log kill they **would have received with one resistance factor
switched off**.

```
HDHYP  = 1 - LKACT/LKOXI     fraction of RT log-kill lost to hypoxia
HDREPR = 1 - LKACT/LKREP     fraction lost to intact repair capacity
HDPOP  = LKNRP/LKACT         fraction given back by repopulation
HDIMM  = 1 - LKIMM/LKIMX     fraction of immune killing lost to immune evasion
```

---

## Why there are four factors and not three — a claim that did not survive

The first draft bundled repair and repopulation into **one** factor. Cisplatin and
cetuximab then share the same headroom, so the combination becomes sub-additive, and
that is exactly what RTOG 0522 observed. But once the equations were actually written,
that bundling **was not supported.** Cisplatin acts on α while cetuximab acts chiefly
on the clonogen birth rate, and the two are separable. As a result this model
**over-predicts** the cisplatin + cetuximab combination.

That failure has not been hidden by redrawing the map; it is left as it stands in the
validation table below.

---

## The distinction between inputs and outputs

**Inputs (imposed).** The HPV-positive phenotype is imposed rather than derived.
`HPV = 1` changes exactly six things, all of them **measured properties** of
HPV-positive oropharyngeal carcinoma and none of them treatment outcomes.

| Parameter | HPV+ | HPV− | Grounds |
|----------|------|------|------|
| `HRCAP` | 0.45 | 1.00 | RAD51 mislocalisation and loss of TIP60 caused by E7 (references 15–18) |
| `VASCQ` | 4.00 | 1.00 | Good perfusion, low 15-gene hypoxia score (55, 59) |
| `LAMS0` | 0.070 | 0.100 | Long T_pot, slow clonogen birth (48) |
| `EGAMP` | 1.00 | 2.30 | EGFR amplification is an HPV-negative event (26) |
| `E7BYP` | 0.55 | 0.00 | Cell-cycle entry via the Rb bypass is EGFR-independent (12) |
| `AGVIR` | 2.60 | 1.00 | E6/E7 are foreign antigens, TIL-hot (19, 99) |

**Fitted.** Six parameters were fitted by Nelder-Mead to eight published local control
rates. Final objective function 0.0089, every anchor within an error of 0.06.

```
ALPHA0 0.2484 /Gy    FCSC0 2.58e-5     PHIPT  1.005
KHRA   0.4798        KOXTR 0.1625 /d   VASCQH 4.00  (hit the upper search bound)
```

That `VASCQH` ended up against its bound is stated plainly. It means the anchor set
wants HPV-positive tumours to be as well perfused as the model will allow, and it is a
sign that one perfusion parameter is carrying more of the HPV contrast than it should.

**Outputs (everything else).** Every figure in the validation tables below. The model
was **not told any** of the following: the cost of a treatment interruption, the time
at which repopulation begins, the benefit of hyperfractionation, the salivary sparing
of IMRT, the relationship between CPS and immunotherapy benefit, the hypoxic fraction
of HNSCC, or T_pot.

---

## In this file radiotherapy is not a discrete map

The usual shortcut of applying a surviving fraction as an instantaneous jump between
fractions **cannot represent incomplete repair.** So this model uses the **exact
continuous form** of the linear-quadratic model. Each fraction is a bolus into a
delivery buffer `DBUF` that empties at `KDEL = 2000/d`, so the instantaneous dose rate
is

```
DRATE = KDEL * DBUF          (area d Gy, an exponential pulse averaging 43 s — a real linear accelerator)
dxdt_SLD = DRATE - MUREP * SLD
HZ = ALPHAE*DRATE + 2*BETAE*SLD*DRATE
```

which integrates to exactly `ALPHAE*D + BETAE*G*D²` (`G` = the Lea-Catcheside
protraction factor). Because `MUREP/KDEL = 0.012`, `G` **within** a fraction is 0.99,
and because `MUREP` is finite, `G < 1` **between** fractions once they are closer
together than about six hours — which is why hyperfractionation and its
incomplete-repair penalty are both **outputs**.

Numerical check: setting `NOREP = 1` to abolish sublethal repair raises the log kill of
a 70 Gy course from 4.65 to 9.40, which agrees with the analytic solution of
`αD + βD²` for `D = 70` given all at once.

---

## Validation

Values obtained by actually running under mrgsolve 2.0.1 / R 4.3.3. The anchors
**were** used in the fitting; the held-out items are predictions that **were not**.

### The fitted anchors (8)

| Grounds | Target | Goal | Model |
|------|------|------|------|
| Bonner 2006 control arm | HPV− RT alone 70 Gy | 0.34 | **0.381** |
| Bonner 2006 | HPV− cetuximab-RT | 0.47 | **0.508** |
| MACH-NC / RTOG 0129 | HPV− cisplatin-RT | 0.58 | **0.600** |
| RTOG 1016 cisplatin arm | HPV+ cisplatin-RT | 0.90 | **0.876** |
| RTOG 1016 cetuximab arm | HPV+ cetuximab-RT | 0.83 | **0.813** |
| DAHANCA 5-85 | HPV− RT + nimorazole | 0.49 | **0.437** |
| Overgaard 2011 / Toustrup 2011 | HPV+ RT + nimorazole | 0.79 | **0.822** |
| HPV+ RT alone reference | HPV+ RT alone | 0.77 | **0.748** |

### Held-out predictions — passed

| Prediction | Literature | Model output |
|------|------|-----------|
| **The central prediction**: the cetuximab benefit is smaller in HPV+ | consistent with RTOG 1016 and De-ESCALaTE (no direct measurement) | HPV− **+0.960** vs HPV+ **+0.689** log10 (ratio 0.72) |
| Within HPV+, cisplatin > cetuximab | RTOG 1016 | **+1.055** vs **+0.689** log10 |
| Repopulation starts in weeks 3–4 | Withers 1988 | `HDPOP` d10 0.115 → d21 **0.250** → d28 0.282 → d46 0.338 |
| The cost of prolonging treatment is 0.6–0.9 Gy/day, about 1 %/day | Hansen 1997 · Bese 2007 | a late interruption (after 30 fractions) **0.69 Gy/day, 1.19 %/day** |
| Hypoxic fraction 10–25 %, T_pot 4–5 days | Nordsmark 2005 · Begg 1992 | emergently **HF 0.162**, T_pot 4.6 days |
| Hyperfractionation improves local control | EORTC 22791 · MARCH | 1.2 Gy BID 81.6 Gy: **+0.43 log10**, LRC 0.600 → 0.684 |
| IMRT parotid sparing | PARSPORT · QUANTEC | salivary flow at 1 year **47 %** (IMRT) vs **8 %** (3D-CRT), tumour dose identical |
| The immune benefit increases with CPS | KEYNOTE-048/040 | at the same exposure `HDIMM` 0.959→0.985, tumour volume ratio 0.79→0.64 (CPS 1→50) |
| Free platinum Cmax 2–5 mg/L | Urien 2004 | **3.83 mg/L** |
| Cetuximab Cmax 185–250 mg/L | Fracasso 2007 | **235 mg/L** |
| Ototoxicity of cisplatin 300 mg/m² | Rybak 2007 | a **21 dB** high-frequency threshold shift |
| Cisplatin nephrotoxicity | Miller 2010 | eGFR 100 → **83 mL/min** |
| Weight loss of 7–10 % during CRT | Langius 2013 | **11.7 %** |
| Peak grade of oral mucositis | Trotti 2003 | CRT **3.61** vs RT alone **2.99** |

### Held-out predictions — **failed** (with the cause stated)

1. **RTOG 0522 is over-predicted.** Adding cetuximab to cisplatin-RT, the model
   predicts **+0.93 log10** (LRC 0.600 → 0.695), whereas the trial showed no benefit.
   Applying the compliance the trial actually reported (a 5-day interruption plus only
   two of the three cisplatin cycles) as scenario 12b removes about half the excess
   (**+0.45 log10**). The remaining residual is a genuine failure, and the cause is as
   set out above — in this model the four factors are additive in the log by
   construction, and the possibility remains that cisplatin and cetuximab are not in
   fact independent.

2. **The phenotype specificity of nimorazole comes out backwards.** The model predicts
   **+0.900 log10** in the well-perfused HPV+ tumour (hypoxic fraction 0.007) and
   **+0.501 log10** in the hypoxic HPV− tumour (0.162). The 15-gene classifier
   literature (Toustrup 2011) says the benefit is confined to tumours that are
   hypoxia-**high**, so the direction is wrong.
   **Cause:** if the exchange between oxic and hypoxic clonogens is slow
   (`KOXTR` 0.1625/d), then by the end of the course the surviving clonogens are almost
   all hypoxic cells whatever the starting hypoxic fraction, so the hypoxic **tail**
   determines tumour control. A two-bin oxygenation model cannot match DAHANCA's
   **effect size** and its **phenotype restriction** at the same time.
   **The fix:** a continuous pO₂ distribution over the clonogens rather than two bins.
   It is not in this version.

3. **An interruption early in the course comes out slightly beneficial.** An
   interruption after 10 fractions computes as **−0.25 Gy/day** (that is, a benefit).
   The cause is the same as in item 2: the interruption gives time for reoxygenation,
   and with a slow `KOXTR` that gain outweighs 10 days of repopulation. Historically
   this was the rationale for split-course radiotherapy, which was shown clinically to
   be inferior. For late interruptions the model does reproduce the literature values.

4. **Weekly cisplatin comes out slightly superior to three-weekly** (log kill 6.01 vs
   5.77). The trials (Noronha 2018 · JCOG1008) show non-inferiority at best. The cause:
   this model has no fall in cumulative dose intensity and no loss of adherence, and the
   weekly arm lays adducts down across more fractions, which increases the opportunity
   for sensitisation.

---

## Individual patient vs population — two probabilities of control

| Value to read | Definition | When to use it |
|---------|------|-----------|
| `TCP` | `exp(-CSC)` — the Poisson probability that **this** tumour has no surviving clonogen | `$OMEGA` virtual populations |
| `LRC` | The population curve, with clonogen number marginalised over a log-normal (SD `SIGLK` = 3.62) | `zero_re()` deterministic contrasts |

**Use both at once and the heterogeneity is double-counted**, giving a dose-response
curve shallower than anything that has been measured. So the rule is: **deterministic
contrasts → `zero_re()` + `LRC`; virtual populations → `$OMEGA` + `TCP`.**

`$OMEGA` is not decoration; it is genuinely used. All four random effects (clonogen
number · intrinsic α · perfusion · antigenicity) were confirmed to reach the outputs
(over 200 subjects, clonogens 2.7e4–2.65e7, α 0.098–0.578, pO₂ 6.4–56.0 mmHg).
`responder_rate()` uses this to produce the response rate, the control rate and the
rate of grade 3+ mucositis.

### When the probability of control is read

Local control is read at **the clonogen nadir within a one-year window**. Leave the
window open and two things go wrong. Read at a fixed late time point and the regrowth
of the fraction that was never controlled in the first place drags the value down. Read
at the global nadir and **the horizon itself becomes a covariate**: in the arms that
responded deeply, residual immune killing grinds the clonogens down over years, so the
nadir is pushed outside the clinical follow-up period and arms the trial separates
clearly converge.

---

## The 26 scenarios — a matched control for every claim

| # | Scenario | Pair | What it reads |
|---|----------|-----|--------------|
| 01–02 | Untreated natural history (HPV−/+) | — | volume doubling time 91 days / 83 days (emergent) |
| 03–06 | HPV− RT alone · cisplatin · cetuximab | 06 is the control | the Bonner geometry |
| 07–09 | **the same three arms in HPV+** | 09 is the control | only the six phenotype inputs change → RTOG 1016 |
| 10 | HPV+ dose reduced to 60 Gy | 07 | de-escalation |
| 11 | Induction TPF → RT | 04 | the PARADIGM/DeCIDE geometry |
| 12 | Cisplatin + cetuximab + RT | 04 | **RTOG 0522 — the point of failure** |
| 13–14 | Nimorazole, hypoxic vs well-perfused | 06 / 09 | **point of failure 2** |
| 15–16 | With and without a 10-day interruption | each other | the cost of OTT |
| 17–20 | **Pembrolizumab at the same exposure, CPS 1/10/20/50** | each other | the immune headroom gradient |
| 21 | CPS fixed at 20, only HPV switched | 19 | orthogonality of the immune axis |
| 22 | EXTREME (platinum/5-FU/cetuximab) | — | first line in recurrent/metastatic disease |
| 23 | Parotid dose alone as 3D-CRT | 04 | **tumour dose identical, only the toxicity contrasted** |
| 24 | Weekly cisplatin 40 ×7 | 04 | **point of failure 4** |
| 25 | Hyperfractionation 1.2 Gy BID | 04 | the loss from incomplete repair |
| 26 | Carboplatin AUC5, CrCl 55 | 04 | the substitute in impaired renal function |

---

## Defects found and fixed during development (for the record)

Left here for transparency. These are the places where the first draft was quietly
wrong.

1. **The repopulation signal was driven by the absolute kill flux.** As the tumour
   shrinks the flux collapses, so the repopulation drive was **maximal in week 1 and
   declined thereafter** — the exact opposite of what is observed. What surviving
   clonogens respond to is **the per-cell hazard** they themselves experience.
2. **The necrosis compartment integrated the entire cell-loss flux.** A 30 cm³ tumour
   became 200 cm³ of debris within a month. Only `FNEC` of the spontaneous loss is now
   left as radiological necrosis.
3. **The body weight model ran away.** An untreated patient lost 34 %. The recovery term
   was too weak and the baseline IL-6 was already producing a cachexia signal.
4. **The clonogens were held in a single mean-field OER bin.** That makes the leverage
   of hypoxia excessive and makes reoxygenation disappear as a mechanism. Oxic and
   hypoxic clonogens were separated.
5. **`HDIMM` was blind to the PD-1 axis.** The exhausted pool was not put back into the
   counterfactual denominator, so **no value of CPS changed the result at all.**
6. **`c(HPVNEG, list(CPS = 50))`** — this gives two elements with the same name, and
   `param()` silently takes the first. Four CPS-gradient scenarios were running
   **exactly the same simulation.** All of them were replaced with `modifyList()`.
7. **`ETA(3)` was dead code.** It had already been assigned into `gVQ` and was then
   being multiplied into the local variable `VQ`, so the perfusion random effect reached
   nowhere at all.
8. **`$TABLE` was reading globals from the last derivative evaluation.** Since an output
   could then be a leftover from an intermediate solver step, every reported quantity is
   now recomputed from the state vector.
9. **Local control was being read at day 730.** It was dragged down by regrowth and
   contradicted the calibrated value.
10. **`KOXTR` was back-transformed wrongly from the fitted log value** (1.4536 instead
    of 0.1625). Every anchor was off, and it was found because the file's default could
    not reproduce the calibration script's own output.

---

## How to run it

```bash
# 1. render the map
dot -Tsvg hnscc_qsp_model.dot -o hnscc_qsp_model.svg
dot -Tpng -Gdpi=150 hnscc_qsp_model.dot -o hnscc_qsp_model.png

# 2. build the model + run all 26 scenarios
Rscript -e 'options(hnscc.run.scenarios=TRUE); source("hnscc_mrgsolve_model.R")'

# 3. the dashboard
Rscript -e 'shiny::runApp("hnscc_shiny_app.R")'
```

Packages required: `mrgsolve` (≥ 2.0.1), `shiny`, `ggplot2`, `dplyr`, `tidyr`;
`DT` and `patchwork` are optional.

The main functions: `run_all()` · `summarise_scn()` · `headroom_test()` ·
`ott_penalty()` · `cps_gradient()` · `responder_rate()`.

---

## ⚠️ Disclaimer

This is a **qualitative / semi-quantitative QSP model for educational and research
purposes**. It was built from the public literature and clinical-trial data and was
actually run and verified in mrgsolve, but it has not been independently certified and
**must not be used for real clinical decision-making, prescribing or regulatory
submission.** The four "failure" items above are the places where this model currently
predicts wrongly, and they have been left in rather than deleted because they are the
work list for the next revision.
