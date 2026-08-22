# Allergic Bronchopulmonary Aspergillosis (ABPA) QSP Model
### Allergic Bronchopulmonary Aspergillosis · ABPA

<p align="center">
  <a href="../../../allergic-bronchopulmonary-aspergillosis/abpa_qsp_model.svg"><img src="../../../allergic-bronchopulmonary-aspergillosis/abpa_qsp_model.png" width="880" alt="ABPA QSP mechanistic map"></a>
</p>

**200-node · 18-cluster mechanistic map · 43-compartment mrgsolve ODE model · 18 treatment scenarios ·
10-tab Shiny dashboard · 97 references, every PMID verified through the PubMed API · dependency-free pure-Python
reference implementation (A0–A13)**

Every number appearing in this README was produced by running
`abpa_reference_implementation.py`, and its complete output is committed alongside as
`abpa_numerical_report.txt`. Where a claim and a number disagree, the disagreement has been
recorded rather than hidden.

---

## The question this model is built around

ABPA is usually described as "a hypersensitivity reaction to a mould", and its treatment as
"suppress the inflammation with steroids and reduce the fungus with antifungals". This model does
not use that description. It is built instead on a single structural claim.

> ***Aspergillus* populations are split across two compartments.**
>
> | | Location | Concentration the drug sees |
> |---|---|---|
> | `FLUM` | Airway lumen · adherent to mucosa | the plasma unbound concentration itself |
> | `FPLG` | **Inside the mucus plug (the sanctuary)** | unbound concentration × `f_pen` (≈ 0.10) |
>
> The two compartments are joined by trapping `k_in·PLUG·FLUM` and release `k_out(PLUG)·FPLG`.
> **The drug that kills the fungus cannot reach the plug, and the drug that dissolves the plug does not kill the fungus.**

Linearising this 2×2 system at low fungal burden gives the eradication condition in closed form.

```
        ⎡ g_l − k_in·PLUG − k_host − E        k_out          ⎤
    J = ⎢                                                    ⎥
        ⎣ k_in·PLUG                     g_p − k_out − f_pen·E ⎦

    Eradication condition:  trace(J) < 0  and  det(J) > 0
```

The only dose term entering the sanctuary row (the second diagonal element) is `f_pen·E`. So if
`g_p > k_out`, the luminal kill rate required is at least `(g_p − k_out)/f_pen`, and into that lower
bound **not one property of the antifungal enters except `f_pen`.** And `k_out` (the plug clearance
rate) is not an antifungal parameter but **a parameter of the type 2 inflammation (mucus) axis**.

That is the whole thesis of this model. Azoles and biologics are not two ways of doing one thing;
they stand in a relation where **one of them lowers the threshold the other has to clear**.

---

## What the numbers say

### ① The threshold is real, and itraconazole barely straddles it — A2

If you multiply through by `f_pen` and keep only the kill rate that reaches the sanctuary, the
required exposure separates as follows.
(Pre-treatment effective `k_out` = **0.2660 /d**, `g_p` = 0.35 /d)

| `f_pen` | Required luminal kill rate E\* (1/d) |
|---|---|
| 0.02 | 4.375 |
| 0.05 | 1.826 |
| **0.10 (model reference value)** | **0.957** |
| 0.20 | 0.505 |
| 1.00 (no sanctuary) | 0.112 |

At steady state on itraconazole 200 mg BID the parent is **1.279 mg/L** and the hydroxy metabolite
**2.202 mg/L** (metabolite/parent 1.72), giving an unbound equivalent of 0.00608 mg/L → **E_lum =
0.954 /d** (61.5% of the class ceiling Emax = 1.55). It **meets the required 0.957 in the third
decimal place.**

This is not the result of tuning; it is the result itself. It explains why response rates to azole
monotherapy are so erratic — **even though the model contains not one responder/non-responder
covariate.** Because the patient is sitting on the blade of the threshold.

### ② As severity rises the threshold passes the class ceiling — A2(d)

`k_out` is not a constant. Eosinophil peroxidase (EPX) mucin cross-linking and established
bronchiectasis both lower `k_out`, so as the disease progresses **E\* rises.**

`.` = itraconazole suffices · `o` = a stronger azole is needed · `X` = **no azole at any dose can do it**

```
  BRON\EPX      1.0      1.6      2.2      2.8      3.4      4.0
  0.0         0.26.    0.73.    1.11o    1.41o    1.64X    1.83X
  3.0         0.49.    0.96o    1.32o    1.59X    1.81X    1.98X
  6.0         0.70.    1.16o    1.49o    1.75X    1.94X    2.10X
  9.0         0.90.    1.33o    1.64X    1.88X    2.06X    2.21X
 12.0         1.07o    1.48o    1.77X    1.99X    2.16X    2.30X
 15.0         1.22o    1.61X    1.89X    2.09X    2.25X    2.38X
 18.0         1.36o    1.72X    1.98X    2.18X    2.33X    2.45X
```

**26 of the 42 cells (62%) exceed the Emax of the entire azole class.** Following the EPX = 1.6 line,
itraconazole becomes insufficient at BRON = **2.96**, and the class as a whole reaches its limit at
BRON = **13.61**. This is the arithmetical account of why late ABPA with established bronchiectasis
is refractory — not because the drugs are weak, but because the fungus is somewhere the drugs cannot
reach.

**The threshold is not fixed, however.** In A9 the state after eight years of untreated progression
(BRON 7.97) has E\* = 1.355, so itraconazole's 0.954 falls **SHORT**; but once treatment reduces the
fungal burden, the antigen → IL-13 → plug loop unwinds with it and E\* falls to 0.708. That is, this
is **a constraint at the moment of decision**, not a permanent sentence. That fact is reported as a
result too.

### ③ Total IgE *rises* on omalizumab — and the ABPA response criterion is written on that curve — A6

| | Baseline | Week 52 | Change |
|---|---|---|---|
| **Total IgE** (the measured quantity) | 1984 IU/mL | **6206 IU/mL** | **×3.13 rise** |
| **Free IgE** (the pharmacologically active species) | 1984 IU/mL | **77.3 IU/mL** | **−96.1%** |
| FcεRI receptor density | 1.000 | 0.730 | −27% |
| Effector activation | 0.984 | 0.662 | −32.7% |

Because omalizumab:IgE complexes are cleared more slowly than free IgE (`kelCX < kdegE`), total IgE
at steady state **must** rise. And yet the ABPA response criterion is "a 35–50% fall in total IgE".

Percentage change in total IgE at week 12 (the time point a trial's primary endpoint usually looks
at):

| Regimen | Total IgE change (week 12) | Free IgE suppression |
|---|---|---|
| Prednisolone, ISHAM taper | **−36%** (criterion met) | — |
| Itraconazole 200 mg BID | −29% | — |
| **Omalizumab 375 mg q2w** | **+213%** (criterion completely failed) | **−96%** |

So this criterion **scores as a failure** a drug that suppresses free IgE by 96%, and scores as a
success one that does not. This is not a matter of parameter choice — in A7, sweeping the complex
half-life 6.5-fold from 4 days to 26 days moves the total-IgE fold change from 1.58× to 9.69×, but
**at all seven values the criterion is not met, and free IgE suppression is always above 89%.**
**The direction is structural; the magnitude is unidentifiable.**

As for the discrepancy with the observed 2–5-fold, it is stated explicitly that the model cannot
adjudicate between three possibilities — (i) complexes are cleared faster than IgG, (ii) the assay
under-detects complexed IgE, (iii) IgE production itself falls (the side Lowe 2011 argues for). No
adjudication is made; it is left as it stands.

A second, quiet result: even at the classical target of free IgE 25 ng/mL, FcεRI occupancy is
**still about 57%** (receptor affinity ≈ 0.10 nM). What omalizumab lives off is not instantaneous
occupancy but the slow fall in **receptor density**, and this is why the clinical effect lags the
biomarker by weeks.

### ④ ABPA begins outside the omalizumab dosing table — and the functional form of that table is derived — A8

Omalizumab neutralises IgE stoichiometrically, so the requirement is not a concentration but a
**molar flux**.

```
375 mg q2w  →  20.41 nM/d supplied
IgE production needed to sustain a baseline of X IU/mL  =  0.00350 · X  nM/d
η = 0.5 (IgE neutralised per molecule of omalizumab)  →  flux crossover = baseline total IgE 2917 IU/mL
```

Integrating the 42-state model in earnest, the cliff appears **at exactly that point.**

| Baseline total IgE | Free IgE suppression (week 52) | Free IgE (ng/mL) | Effector reduction |
|---|---|---|---|
| 500 | 98.8% | 14.0 | 87.3% |
| 1000 | 98.5% | 36.6 | 71.0% |
| 1500 | 97.8% | 79.2 | 52.5% |
| 2000 | 96.1% | 185.7 | 32.7% |
| **3000** | **54.4%** | 3283 | **4.3%** |
| 5000 | 13.3% | 10399 | 0.4% |
| 12000 | 4.8% | 27415 | **0.0%** |

The dose required to achieve free IgE < 25 ng/mL is **very nearly exactly proportional** to baseline
IgE — **0.00685 ~ 0.00697 mg/(IU/mL·kg)** across a 20-fold range, a variation of under 1.8%:

| Baseline IgE | Required q2w dose | vs the table maximum of 375 mg | mg/(IU/mL·kg) | In the table? |
|---|---|---|---|---|
| 500 | 240 mg | ×0.64 | 0.00685 | yes |
| 1000 | 484 mg | ×1.29 | 0.00692 | **no** |
| 1500 | 728 mg | ×1.94 | 0.00694 | **no** |
| 2500 | 1217 mg | ×3.25 | 0.00695 | **no** |
| 5000 | 2438 mg | ×6.50 | 0.00697 | **no** |
| 10000 | 4881 mg | ×13.02 | 0.00697 | **no** |

The model has never seen the approved dosing table. Nevertheless it **derives the functional form of
that table (dose ∝ IgE × body weight) from the molar flux balance**, and places the crossover in the
same neighbourhood as the point where the table ends (total IgE 1500 IU/mL). The diagnostic criterion
for ABPA is total IgE > 1000 and typical baselines are 2000–5000 — **this disease begins where the
label ends.** This recasts "omalizumab non-response in ABPA" as a dosing problem with an arithmetical
answer, and predicts that response should correlate not with mg but with **mg/(IU/mL·kg)**.

### ⑤ Itraconazole's steroid-sparing effect depends on which steroid was used — A4·A5

Residual CYP3A4 activity at itraconazole steady state is **I = 0.2932**. From this it follows first
that **a single-site model is arithmetically impossible**: the ceiling on the fold change in a model
where only hepatic clearance is inhibited is `1/I` = **3.41**, whereas the observed
itraconazole–budesonide interaction is **4.2-fold**. The interaction must therefore occur at **two
sites in series, the gut wall and the liver**, and the fold changes multiply.

| Steroid | s_gut | s_hep | Model AUC fold | Observed | Error |
|---|---|---|---|---|---|
| Prednisolone | 0.00 | 0.15 | **1.119** | 1.24 | **−9.8%** (prediction, not a fit) |
| Methylprednisolone | 0.15 | 0.80 | 2.574 | 2.60 | −1.0% (fitted) |
| Budesonide (inhaled) | 0.70 | 0.70 | 3.917 | 4.20 | −6.7% (fitted) |

Two were fitted and **prednisolone was predicted**. It is under-predicted by 10%, and that fact is
recorded rather than erased by tuning.

And then the decomposition in A5 — how the 52-week FEV₁ gain from adding itraconazole to an
identical fixed steroid regimen divides between antifungal effect and drug interaction, measured on
**the exposure integral (AUC of GR effect) rather than the terminal trough**:

| Steroid background | AUC fold | Delivered GR exposure fold | ΔFEV₁ total | Antifungal | **DDI** | DDI share |
|---|---|---|---|---|---|---|
| Prednisolone 10 mg/d | 1.12 | 1.09 | 0.784 | 0.254 | **0.509** | 64.9% |
| Methylprednisolone 8 mg/d | 2.64 | 1.85 | 4.168 | 0.219 | **3.722** | 89.3% |
| Budesonide 1600 µg/d inhaled | 4.06 | 3.47 | 6.364 | 3.257 | **6.134** | 96.4% |

**I was wrong twice here and I report both.** At first I expected the interaction to account for
most of the effect on every background; then, seeing that the AUC rise for prednisolone is only
1.1-fold, I corrected myself to say that prednisolone-based trials are safe. **Neither survives its
own table.** Because two things change at once in that table — the DDI contribution grows with
CYP3A4 dependence (0.51 → 3.72 → 6.14), and the antifungal contribution is not constant either
(0.25 → 0.22 → 3.26). The latter is not because the fungus recognises steroids but is a **headroom
effect**: prednisolone 10 mg/d is already controlling the disease (FEV₁ 80.3, plug 1.08) so there is
almost nothing for the antifungal to earn, whereas inhaled budesonide alone barely controls it at all
(FEV₁ 69.9, plug 4.40) so the same pharmacology moves 6 points.
**The two effects push in the same direction, and efficacy endpoints alone cannot separate them.**

> **A falsifiable prediction:** an itraconazole steroid-sparing trial should report an effect
> **several times larger** on a methylprednisolone or inhaled-budesonide background than on a
> prednisolone background, and the difference is pharmacokinetic, not antifungal. **This comparison
> has never been performed.**

### ⑥ Where the drug interaction stops being a confounder and becomes the adverse event itself — A10

| Regimen | Budesonide (mg/L) | GR effect | **Morning cortisol** | BMD | HbA1c |
|---|---|---|---|---|---|
| Inhaled budesonide 1600 µg alone | 0.00000 | 0.0007 | 13.43 | 0.9967 | 5.49 |
| Itraconazole alone | 0.00000 | 0.0000 | 14.00 | 1.0000 | 5.40 |
| **Budesonide + itraconazole** | 0.00004 | 0.0216 | **11.96** | 0.9888 | 5.70 |
| Prednisolone 10 mg + itraconazole | — | 0.0811 | **8.43** | 0.9687 | 6.23 |

An inhaled steroid is chosen in order to avoid systemic exposure. Add a CYP3A4 inhibitor and the
first-pass step that made that true disappears. **Nothing in the disease model changed between rows
1 and 3** — the adrenal suppression is manufactured entirely in the PK block.

### ⑦ How long does emptying the sanctuary directly last — A11

| Regimen | Time for the sanctuary to refill to 50% of the untreated level |
|---|---|
| Bronchoscopic lavage alone | **56 days** |
| Lavage + itraconazole | not reached within 52 weeks |
| Lavage + itraconazole + dupilumab | not reached within 52 weeks |

Mechanical removal is the only intervention that takes `FPLG` to zero instantly, but the untreated
state is an **attractor**, so half of it is back in eight weeks. Emptying the sanctuary is not a
treatment; it is **the manoeuvre that puts the system into a state the drugs can reach**.

### ⑧ No Bliss synergy was found — reported as a negative result — A9

The first A9 ran a Bliss independence test in early disease, and because the two monotherapies
already reduce the fungal burden by more than 98%, every index sat at 1.0 and **the test itself had
no dynamic range.** Re-running it in the eight-year progressed state gives the same answer — **the
excess is slightly negative (sub-additive) on every index.** This model does not manufacture synergy
on fungal-burden endpoints, and that is how it is reported.

What the model actually says lies in two other columns. (i) **Composition**: even at a similar total
burden, the sanctuary share falls from 54.4% on itraconazole alone to **16.0%** on the combination —
the fungus that survives ends up somewhere the drug reaches. (ii) **Threshold**: dupilumab raises
`k_out` from 0.2236 to 0.3657 and so lowers E\* from 1.355 to 0.017.

### ⑨ The only irreversible endpoint — A12 (5 years)

| Regimen | Bronchiectasis /18 | Δ | FEV₁ | Exacerbations | Cumulative OCS (mg) | BMD | Net utility |
|---|---|---|---|---|---|---|---|
| No treatment | 5.934 | +2.934 | 62.28 | 6.41 | 0 | 1.000 | −11.13 |
| Intermittent prednisolone (rescue only) | 5.588 | +2.588 | 62.94 | 5.66 | 2 940 | 0.980 | −15.10 |
| Continuous itraconazole | 4.097 | +1.097 | 69.92 | 2.25 | 0 | 1.000 | −5.57 |
| Continuous dupilumab | 3.390 | +0.390 | 78.18 | 1.14 | 0 | 1.000 | −3.05 |
| **Dupilumab + itraconazole** | **3.382** | **+0.382** | 78.32 | 1.12 | 0 | 1.000 | **−3.02** |
| Prednisolone 10 mg maintenance | 3.307 | +0.307 | **79.78** | **0.90** | **18 250** | **0.750** | **−34.87** |

Prednisolone maintenance comes **first on every pulmonary index** and **overwhelmingly last on net
utility** (cumulative 18 250 mg, 25% loss of bone mineral density). Intermittent rescue therapy is
worse than no treatment at all — because it arrives late every time and pays the steroid price for
nothing. This is why one must not look at a single endpoint, and why the Shiny app never shows a
steroid curve without cortisol, BMD, HbA1c and cumulative dose beside it.

### ⑩ The answer this model gives is owned by a single never-measured number — A13

The +20% local sensitivities (elasticities `d ln y / d ln θ`) are reported split across **two E\*
columns**. `E*` is recomputed at the state the regimen actually reached (so a parameter can move it
**indirectly** by changing `PLUG` or `k_out`), while `E*|fixed` is computed at a fixed reference state
(PLUG 1.0, untreated `k_out`) and so isolates **the closed form itself**.

| Parameter | Sanctuary burden `FPLG` | `E*` (state allowed to move) | **`E*|fixed`** (closed form) |
|---|---|---|---|
| `g_p` (growth inside the plug) | **+482.8** | +15.118 | **+3.464** |
| `kout0` (plug clearance) | −2.821 | −5.000 | **−2.327** |
| `g_epx` (EPX cross-linking) | +2.392 | +2.798 | **+1.068** |
| `f_pen` (penetration) | −1.044 | −0.711 | **−0.775** |
| `k_host` (host killing) | −0.848 | −0.644 | −0.086 |
| **`Emax_af`** (antifungal potency) | −1.016 | −0.046 | **0.000** |
| **`EC50_af`** (antifungal efficacy) | +0.823 | +0.023 | **0.000** |

**The potency and the efficacy of the antifungal are exactly 0.000 in the closed form.** Not by
coincidence: it is the fact that no antifungal property other than `f_pen` enters the sanctuary row
of A2's Jacobian, appearing unchanged in a finite-difference sweep of the 42-state model. The small
residual in the `E*` column (−0.046) is **a state effect, not a potency effect** — killing more
fungus lowers the threshold along the antigen → IL-13 → plug → `k_out` path — and it is one to two
orders of magnitude smaller than that of `g_p`. **Everything that actually moves the threshold
belongs to the mucus axis.**

And the elasticity of the sanctuary burden with respect to `g_p` **varies 673-fold with the
regimen** — **0.72** on a well-controlled regimen (itraconazole + prednisolone), **482.8** on azole
monotherapy. Sensitivity does not merely grow near the eradication boundary; it **diverges.** Had the
parameter sweep been run only in the well-treated arms, `g_p` would have been reported as
unimportant.

> **`g_p` (the growth rate of *Aspergillus* inside the plug) has never been measured in the human
> airway.**
> Everything this model says about "can an azole clear ABPA" is a statement about `g_p`, and A2(e)
> maps where the answer flips: at `g_p` ≤ 0.20 it clears without any drug, at 0.29 itraconazole
> suffices, at 0.35 a stronger azole is needed, and **above 0.45 no azole is possible.** Anyone
> wishing to refute this model should measure that number rather than argue about the rest. The
> absence of a citation for `g_p` in the references is not an omission but a refusal to hide that
> fact.

---

## Files

| File | Contents |
|---|---|
| [`abpa_qsp_model.dot`](../../../allergic-bronchopulmonary-aspergillosis/abpa_qsp_model.dot) | 200-node · 18-cluster mechanistic map (Graphviz source) |
| [`abpa_qsp_model.svg`](../../../allergic-bronchopulmonary-aspergillosis/abpa_qsp_model.svg) · [`abpa_qsp_model.png`](../../../allergic-bronchopulmonary-aspergillosis/abpa_qsp_model.png) | Renders (`dot -Tpng -Gdpi=150`) |
| [`abpa_mrgsolve_model.R`](../../../allergic-bronchopulmonary-aspergillosis/abpa_mrgsolve_model.R) | 42-compartment mrgsolve model + scenario builder |
| [`abpa_reference_implementation.py`](../../../allergic-bronchopulmonary-aspergillosis/abpa_reference_implementation.py) | Dependency-free pure-Python reference implementation (A0–A13) |
| [`abpa_numerical_report.txt`](../../../allergic-bronchopulmonary-aspergillosis/abpa_numerical_report.txt) | The full output of the above script — the source of every number in this README |
| [`abpa_shiny_app.R`](../../../allergic-bronchopulmonary-aspergillosis/abpa_shiny_app.R) | 10-tab Shiny dashboard |
| [`abpa_references.md`](../../../allergic-bronchopulmonary-aspergillosis/abpa_references.md) | 97 references, every PMID verified against the API + annotations on their role in the model |

```bash
# render the map
dot -Tsvg abpa_qsp_model.dot -o abpa_qsp_model.svg
dot -Tpng -Gdpi=150 abpa_qsp_model.dot -o abpa_qsp_model.png

# numerical suite (standard library only, ~8 min)
python3 abpa_reference_implementation.py          # all of A0-A13
python3 abpa_reference_implementation.py A2 A8    # individual analyses

# mrgsolve / Shiny
Rscript -e 'library(mrgsolve); mod <- mread("abpa_mrgsolve_model", "."); print(mod)'
Rscript -e 'shiny::runApp("abpa_shiny_app.R")'
```

---

## ODE compartments (43 in mrgsolve · 42 in Python)

The **163 shared parameters of the two files agree to machine precision** (the R side adds only the
two counterfactual switches `DDI_OFF` and `AF_OFF`), and the compartment count differs by one because
the Python implementation handles the omalizumab SC depot not as a compartment but as an analytic
integral between events — since there is no feedback into the depot, the two are equivalent.

| Group | Compartments |
|---|---|
| Antigen · type 2 axis | `AG` `TH2` `IL13` `IL5` |
| Plasma cell · IgE TMDD | `PC` `TI` (total IgE) `TO` (total omalizumab) `OMAD` (R only) `FCER` |
| Eosinophils · granule chemistry | `EOSB` `EOSA` `EPX` |
| **Mucus plug** | `PLUG` |
| **Sanctuary split** | `FLUM` **`FPLG`** |
| Azole PK (non-linear) | `AITR` `ITRA` `OHIT` `AVOR` `VORI` |
| Inhaled antifungal | `AMB` |
| Steroid PK | `APRD` `PRED` `AMPD` `MPRD` `ABUD` `BUD` |
| HPA · toxicity | `CORT` `BMD` `HBA1C` `CUMO` |
| Biologic PK | `MEPD` `MEPO` `BEND` `BENR` `DUPD` `DUPI` `TEZD` `TEZE` |
| Structure · hazard | `BRON` `CHAZ` |
| Exposure integrals | `AUCP` `AUCCS` |

Each treatment is classified by **which term of the split it moves**:

| Treatment | Term moved |
|---|---|
| Itraconazole · voriconazole | the whole of `E_af` in `FLUM`, × `f_pen` ≈ 0.10 in `FPLG` |
| Inhaled liposomal amphotericin B | `E_af` at the mucosal surface; no systemic CYP3A4 interaction, but `f_pen` = 0.03 (worse) |
| Corticosteroids | `PLUG` falls quickly, but through inhibition of `k_host` **`FLUM` rises** |
| Omalizumab | FcεRI occupancy only; `PLUG` slightly, `FPLG` not at all |
| Mepolizumab · benralizumab | less EPX → less mucin cross-linking → `k_out` rises |
| Dupilumab | less MUC5AC drive → `PLUG` falls → `k_out` rises → **E\* falls** |
| Bronchoscopic lavage | `FPLG` and `PLUG` to ~0 instantly (but half back by day 56) |

---

## Calibration anchors

| Target | Anchor | Model value |
|---|---|---|
| Untreated phenotype | total IgE 2000 IU/mL · eosinophils 800/µL | 1984 · 797 |
| | FEV₁ ~68% pred · plug score 6/18 | 67.9 · 6.04 |
| Itraconazole PK | metabolite/parent 1.5–2-fold | 1.72 (parent 1.279, OH 2.202 mg/L) |
| Voriconazole PK | 200 mg BID → ~2 mg/L (NM) | mixed MM + linear clearance |
| CYP3A4 inhibition | methylprednisolone AUC ×2.6 (fitted) | 2.574 |
| | prednisolone AUC ×1.24 (**predicted**) | 1.119 (−9.8%) |
| | budesonide AUC ×4.2 (fitted) | 3.917 |
| Omalizumab | free IgE −90~99%, total IgE ×2–5 | −96.1%, ×3.13 |
| Steroid | prednisolone 10 mg/d maintenance → cortisol suppression | 14.0 → 8.9 µg/dL |
| Numerical verification | solutions agree across dt 0.005–0.04 | identical to four decimal places |
| | IgE flux balance residual | 3.3 × 10⁻⁴ nM/d |

---

## What the model gets wrong or cannot answer (stated failures and limits)

Every defect caught during the work, and every remaining limitation, is written down here.

1. **`g_p` has never been measured, and it owns the answer.** (⑩) The central claim of this model is
   falsifiable by a single measurement, and that measurement does not exist.
2. **The prednisolone interaction is under-predicted by 10%** (1.119 vs 1.24). It was not erased by
   fitting.
3. **The fold rise in total IgE is unidentifiable.** Only the direction is structural. The model's
   3.13-fold against an observed 2–5-fold is a calibrated value, and the model cannot adjudicate
   between the three explanations. (③)
4. **No Bliss synergy was found.** Dupilumab + azole is slightly sub-additive on every index of the
   fungal-burden endpoints. The threshold and the composition move, but that is not synergy. (⑧)
5. **The bronchiectasis integrator assumes causation.** Whether ABPA *causes* bronchiectasis or
   merely coexists with it is unsettled in the literature (references §⑪, Agarwal 2026). Every number
   in ⑨ inherits that assumption.
6. **The outcome layer has no mortality and no quality of life.** The net-utility weights are
   illustrative, not fitted.
7. **Azole resistance is on the map but not in the ODEs as dynamics.** `RESIST_SEL` exists only as a
   node; the process by which sub-MIC exposure inside the sanctuary selects for resistance is not
   modelled over time.
8. **Five real defects caught during development** (all fixed and re-run):
   a wrong 1000-fold unit factor in the oral input term · pure Michaelis–Menten voriconazole
   clearance diverging when the input rate exceeds Vmax · the omalizumab effector arm not being wired
   downstream, so that scenario 3 was identical to no treatment · A8 setting baseline IgE as an
   initial condition so that the model returned to its own set point and every row came out identical ·
   initial conditions not lying on the steady-state manifold, so untreated total IgE overshot to 3450
   and contaminated the 12-week readout. The last two genuinely changed the results tables, and the
   reason for each fix has been left in the code as a comment.

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for education and research. It was assembled from the public
literature but has not been independently validated or certified, and **must not be used for real
clinical decision-making, prescribing, or regulatory submission.** In particular, this model contains
**falsifiable claims** about the interpretation of clinical trials (④ · ⑤). Those are hypotheses, not
validated conclusions.
