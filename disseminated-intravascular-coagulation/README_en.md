# Disseminated Intravascular Coagulation (DIC) QSP Model
## Disseminated Intravascular Coagulation — Sepsis-induced coagulopathy · APL coagulopathy

<p align="center">
  <a href="dic_qsp_model.svg"><img src="dic_qsp_model.png" width="900" alt="DIC QSP mechanistic map"></a>
</p>

---

## The question this model sets out to answer

DIC is always described in the textbooks as having two faces. In sepsis,
microthrombi wreck the organs; in acute promyelocytic leukaemia (APL) the patient
dies of an intracerebral haemorrhage. The same name is attached to both, the exact
opposite thing happens, and so the treatment is opposite too — an anticoagulant on
one side, a haemostatic agent on the other.

But **why** does it become the opposite? The common explanation is a story about
temporal sequence: "as the consumptive coagulopathy progresses the factors are
eventually depleted and it tips over into bleeding". That explanation is wrong.
Septic DIC patients usually do not bleed even after several days, and the APL
patient bleeds from the first day.

This model is built so that the answer is **derived arithmetically rather than
assumed**. The 49 differential equations are **completely identical** for the two
phenotypes. All that differs is **four** parameters, and those four set **two
mutually independent clocks**.

---

## The structural decision — one consumption process, two independent clocks

The whole structure of this model is one sentence.

> DIC is **one** tissue-factor-driven process of thrombin generation and
> consumption, and the clinical phenotype is the result of reading that process
> through **two independently set clocks**.

Every coagulation protein is consumed at the **same** fractional consumption rate
`KCONS = KC0·THR/(KTHRC+THR)`. The consumption side is common. What diverges is the
**synthesis side** and the **clearance side**.

---

## CLOCK 1 — the **sign** of the acute-phase term

The hepatic synthesis term is multiplied by IL-6. But the sign differs from protein
to protein.

```
fibrinogen       synthesis × (1 + APOSF · APF)     ← positive acute-phase protein
antithrombin     synthesis ÷ (1 + ANEGAT · APF)    ← negative acute-phase protein
protein C        synthesis ÷ (1 + ANEGPC · APF)    ← negative acute-phase protein
free protein S   synthesis ÷ (1 + ANEGPS · APF)    ← via the rise in C4BP
```

A single `APF = ΔIL-6/(KIL6 + ΔIL-6)` acts in two directions.
Over a 48-hour simulation **the consumption rates of the two phenotypes are
effectively the same** (sepsis KCONS = 0.0184/h, APL 0.0204/h — APL is if anything
slightly higher). Multiply that same consumption rate by each one's acute-phase term
and solve for the steady state:

```
fibrinogen_ss = kd·FIB0·(1 + APOSF·APF) / (kd + KCONS·WFIB + fibrinogenolysis)

sepsis  IL-6 1400 pg/mL → APF 0.777 → synthesis ×3.33
        0.00693 × 300 × 3.33 / (0.00693 + 0.0184 + 0.0002) = 272 mg/dL
APL     IL-6  240 pg/mL → APF 0.370 → synthesis ×2.11
        0.00693 × 300 × 2.11 / (0.00693 + 0.0204 + 0.0212) =  91 mg/dL
```

**At almost the same consumption rate: 272 mg/dL and 91 mg/dL.**

The two clocks contribute about half each. Remove APL's fibrinogenolysis term
(0.0212/h, which is what CLOCK 2 created) and APL's steady state becomes
161 mg/dL. Conversely, leave APL's CLOCK 2 as it is and change only the acute-phase
term to the sepsis value and 91 → 143 mg/dL. That is, **272 → 161 is the difference
CLOCK 1 makes and 161 → 91 the difference CLOCK 2 makes.**

- The septic patient is in a state where **the brake is off but the substrate keeps
  being replenished**. → uninhibited thrombin + normal fibrinogen =
  **microthrombosis**, not bleeding.
- The APL patient has neither brake nor substrate. → **bleeding**.

**The reason septic DIC does not bleed is that fibrinogen is a positive acute-phase
protein.** That is not a rule put into the model but a result that came out of the
equations.

### The second consequence of CLOCK 1 — why the ISTH score misses it

The ISTH overt-DIC score awards fibrinogen one point **only when it is below
100 mg/dL**. But by the arithmetic above the septic patient's fibrinogen does
**not** fall that far. Looking at the 24-hour time point in the simulation:

| | Value | ISTH | SIC |
|---|---|---|---|
| Platelets | 135 ×10⁹/L | 0 points | 1 point |
| D-dimer | 4.3 mg/L | 2 points | (no such item) |
| PT prolongation | +2.4 s (INR 1.202) | 0 points | 1 point |
| Fibrinogen | 281 mg/dL | **0 points** | (no such item) |
| SOFA | 4.3 | (no such item) | 2 points |
| **Total** | | **2 points (non-overt)** | **4 points (SIC confirmed)** |

Same patient, same moment, and ISTH says "not overt DIC" while SIC already says
"coagulopathy present". Yet at this moment deposited fibrin `FDEP` is already
0.88 AU (68% of its 48-hour value) and SOFA is 4.3.

**The reason the SIC score deliberately left fibrinogen out comes out of the
equations.** The model was never told that "ISTH has low sensitivity". The two
scores were computed as **outputs** of the state vector, and that is how they came
out.

---

## CLOCK 2 — the **molar ratio** of PAI-1 : t-PA

t-PA and PAI-1 bind 1:1 and effectively irreversibly. What matters, therefore, is
not each concentration but the **molar ratio**. Because their molecular weights
differ (t-PA 68 kDa, PAI-1 45 kDa) this ratio is invisible in units of ng/mL.

```
t-PA   nM = ng/mL × 0.0147
PAI-1  nM = ng/mL × 0.0222 × (active fraction 0.60)

free t-PA = t-PA_nM × [ FLOC + (1−FLOC)/(1 + PAI-1_nM/KPAI) ] × annexin A2 amplification
```

Simulated 48-hour values:

| | t-PA | PAI-1 | **molar ratio** | **free t-PA** |
|---|---|---|---|---|
| Normal | 5.0 ng/mL | 20 ng/mL | 6.0 : 1 | 0.0336 nM |
| **Septic DIC** | 20.2 ng/mL | **564 ng/mL** | **42.1 : 1** | **0.0222 nM** |
| **APL** | 20.9 ng/mL | 53 ng/mL | 3.9 : 1 | **0.2882 nM** |

**In sepsis total t-PA rose four-fold while free t-PA fell to 0.66 times normal.**
This is the quantitative identity of "fibrinolytic shutdown". Measuring t-PA is
meaningless; **the ratio is what must be measured**. APL's free t-PA is **13 times**
that of sepsis.

### The consequence of CLOCK 2 — the sign of tranexamic acid flips by itself

TXA enters the equations in exactly **one place**, with one sign:

```
vPLN = ( KPG · free_tPA · PLG · FSURF  +  uPA term ) × TXAI ,    TXAI = 1/(1 + C_TXA/IC50)
```

This single term produces opposite clinical outcomes in the two phenotypes.

| | Sepsis (shutdown) | APL (hyperfibrinolysis) |
|---|---|---|
| Plasmin | 0.30 → 0.07 nM | 2.32 → 0.70 nM |
| Fibrinogen | 245 → 245 mg/dL (no change) | 101 → 138 mg/dL (**improved**) |
| Deposited fibrin FDEP | 1.29 → **1.59** (worse) | 0.19 → **0.51** (worse) |
| Bleeding index | 0.71 → 0.59 | 0.98 → 0.94 |
| 28-day mortality | 43.3% → **51.1%** | 32.8% → **36.2%** |

In sepsis it switches off the last of a fibrinolysis that was already close to zero,
so the stock grows. In APL fibrinogen recovers, but the fibrin stock grows at the
same time.

The model predicts that **TXA fails to produce a net benefit even in APL**. This
agrees with Avvisati's randomised trial (no benefit) and with the ELN
recommendation (do not use antifibrinolytics routinely, particularly given the
thrombotic risk when combined with ATRA). In the scenario that adds TXA on top of
ATRA, 28-day mortality **worsens** from 7.5% → 21.8%.

---

## Why the clinical trials read the way they do — three products

### Product (1) Heparin is a **catalyst**, and the AT term **saturates**

The textbook claim is "AT falls in DIC, therefore heparin does not work". Written
linearly (effect ∝ AT) the effect at AT 40% falls to 40%, which is far larger than
what is observed clinically.

The real chemistry is a catalytic reaction. One molecule of heparin mediates many
AT-thrombin reactions and saturates in AT in Michaelis-Menten fashion. So the model
writes it like this:

```
thrombin inactivation = KATII × [ AT                     ← spontaneous reaction, linear in AT
                                + (AMP−1) × ATcat ]      ← heparin catalysis, saturating in AT
        ATcat = AT/(KMAT + AT) × (KMAT + 1),   KMAT = 0.25
```

| AT activity | No heparin | UFH 0.5 IU/mL | Fold gain |
|---|---|---|---|
| 100% | 2.20 /h | 17.73 /h | 8.1× |
| 60% | 1.32 /h | 15.02 /h | 11.4× |
| **40%** | 0.88 /h | **12.83 /h** | 14.6× |
| **220%** | **4.84 /h** | 22.27 /h | 4.6× |

Two things come out here at once.

1. **At AT 40% the heparin ceiling is 72% of the normal-AT ceiling.** Not 68% —
   72%. It is a loss of **ceiling** that raising the dose does not recover, not a
   loss of potency. And its magnitude is far smaller than the 40% a linear model
   would claim — consistent with the clinical evidence saying "AT deficiency is not
   the sole cause of heparin resistance". (The model also carries a separate non-AT
   resistance term `KHEPN` due to PF4 and histones.)

2. **Even raising AT to 220% earns only +26% in the heparin arm.** Yet **in the
   uncatalysed spontaneous arm, 60% → 220% is a 3.7-fold change.**

That is, **antithrombin supplementation works through precisely the arm that
heparin renders irrelevant.** Give heparin alongside it and the incremental benefit
of AT supplementation cannot help but shrink. This is the arithmetical candidate for
why KyberSept was negative overall and positive **only in the subgroup not receiving
heparin**. In the simulation too, the increment from AT concentrate is −7.3%p
without heparin and −5.3%p with heparin, and the bleeding burden of the combination
arm is larger.

### Product (2) APC generation is a **product of four terms**

```
APC generation = KPCA × [thrombin] × TM × EPCR × [protein C]
```

Values from the 48-hour septic DIC simulation:

```
TM 0.318  ×  EPCR 0.559  ×  protein C 0.333  =  0.059
```

**5.9% of normal.** And the APC generated per unit of thrombin falls **16.6-fold**
relative to normal (0.636 → 0.038 nM APC per nM thrombin). Because three fractions
are multiplied together, not one of them is catastrophically low and yet the result
is catastrophic.

This is the pharmacological logic of giving **APC itself** (drotrecogin alfa)
instead of supplementing protein C — it bypasses all three terms. In the simulation
early drotrecogin alfa gives 43.3% → 37.4% (−5.9%p), almost the same as PROWESS's
−6.1%p.

### Product (3) Haemostatic failure is a product too, so **the worst term dominates**

```
haemostatic capacity = (1−h_platelet)(1−h_fibrinogen)(1−h_PT)(1−h_lysis)(1−h_drug)
bleeding index = 1 − haemostatic capacity
```

Comparing what to add for an APL patient already receiving ATRA (48 hours,
scenarios O/P/R):

| | Platelets | Fibrinogen | h_platelet | h_fibrinogen | h_lysis | **Bleeding index** | 28-day mortality |
|---|---|---|---|---|---|---|---|
| O · ATRA only | 33 | 108 | 0.664 | 0.416 | 0.784 | 0.966 | 7.5% |
| P · + fibrinogen concentrate | 32 | **218** | 0.675 | **0.103** | 0.840 | 0.961 | 9.2% |
| R · + platelet transfusion | **54** | 105 | **0.397** | 0.435 | 0.788 | **0.944** | **6.7%** |

Doubling fibrinogen from 108 → 218 mg/dL improves `h_fibrinogen` four-fold, from
0.416 → 0.103, **and yet the bleeding index barely moves, 0.966 → 0.961.** That is
because the other terms of the product (`h_platelet` 0.67, `h_lysis` 0.78) are
unchanged. Raise the platelets from 33 → 54, by contrast, and the bleeding index
moves further, 0.966 → 0.944, and mortality is the lowest of the three.

The clinical maxim that **the worst variable must be fixed first** comes out of the
multiplication.

Row P is, in addition, one of the model's clinically contentious predictions:
supplementing fibrinogen alone barely improves the bleeding, while amounting to
feeding substrate into a thrombin-driven process, so deposited fibrin rises from
0.19 → 0.30 and net mortality gets worse rather than better, 7.5% → 9.2%. This is
not a validated result and should be read as a **falsifiable model prediction**.

---

## Stock and rate are different state variables

The model splits fibrin into **two** states.

```
deposition rate   vDEP = KDEP · THR^1.5/(KTHD^1.5 + THR^1.5) · (FIB/FIB0)   ← rate
deposited stock   dFDEP/dt = vDEP − KLYS·PLN_eff·FDEP − KRES·FDEP           ← stock
organ damage      ∝ FDEP                                                    ← driven by the stock
```

Every anticoagulant acts on `vDEP`, and organ failure is decided by `FDEP`. And in
the shutdown state `PLN_eff ≈ 0`, so the stock has essentially only one exit,
`KRES` — a half-life of about 87 hours. The stock persists for days.

### But this structure **refuted** the author's prior hypothesis

I expected this structure to be able to explain PROWESS → PROWESS-SHOCK: "a drug
that only acts on the rate is useless if given after the stock has already
accumulated." Sweeping the dosing start time shows that it is not so.

| Dosing start | 0 h | 6 h | 12 h | 24 h | 30 h | 48 h | 72 h | 96 h |
|---|---|---|---|---|---|---|---|---|
| 28-day mortality | 37.8% | 37.4% | 36.9% | **36.1%** | **35.8%** | **35.9%** | 37.3% | 39.1% |

(No treatment 43.3%.) The benefit does not decline monotonically but is **maximal at
24–48 hours**, and only collapses beyond 72 hours. The maximum tracks the interval
over which fibrin deposition is fastest. Timing matters, but not in the direction of
"too late", and therefore this model does not explain PROWESS-SHOCK. The cause most
likely lies elsewhere (the recruited population, concomitant steroids, improvements
in standard care).

I leave this in rather than deleting it. The worth of a model lies in its being able
to refute the expectations of the person who built it.

---

## What the half-life ordering makes — and what it does not

| Protein | Half-life | Time constant τ under consumption | 48-hour simulation value |
|---|---|---|---|
| **protein C** | **6 h** | **7.5 h** | **33%** |
| Factor VII | 5 h | 6.7 h | 71% |
| Factor VIII | 10 h | 10.7 h | **123%** (acute-phase rise) |
| Factor V | 15 h | 14.1 h | 51% |
| Free protein S | 42 h | — | 39% |
| Factor X | 40 h | 28.8 h | 51% |
| Prothrombin | 65 h | 35.6 h | 48% |
| Antithrombin | 65 h | 49.8 h | 60% |
| Fibrinogen | 100 h | 39.6 h | 245 mg/dL (**normal**) |
| Platelets | 200 h | ~20 h | 78 ×10⁹/L |

protein C has τ = 7.5 hours and therefore **reaches its floor in about 22 hours**.
Antithrombin has a τ of 50 hours and fibrinogen 40 hours, and fibrinogen is on top
of that pushed up by acute-phase synthesis so that it never leaves the normal range.

**So the anticoagulant brake collapses long before the coagulation substrate is
depleted.** That is why the narrative of "consumption progresses and it tips over
into bleeding later" does not hold in sepsis. Interestingly, the fall in
antithrombin is mostly **not due to consumption**. Decomposing the AT loss pathways
at the 48-hour time point gives normal turnover 0.00642/h, capillary leak from
glycocalyx damage 0.00568/h (comparable in size to normal turnover) and elastase
degradation 0.00091/h, with synthesis suppressed by ×0.865 on top of that. The
stoichiometric consumption from thrombin inactivation is smaller than any of these.
**The only term heparin can reach is the smallest one of them** — and heparin, of
all things, makes that term larger.

---

## Files

| File | Contents |
|---|---|
| [`dic_qsp_model.dot`](dic_qsp_model.dot) | Mechanistic map source — 139 nodes · 17 clusters · 242 edges |
| [`dic_qsp_model.svg`](dic_qsp_model.svg) / [`.png`](dic_qsp_model.png) | Rendered map |
| [`dic_mrgsolve_model_en.R`](dic_mrgsolve_model_en.R) | 49-ODE mrgsolve model + 16 treatment scenarios + calibration record |
| [`dic_shiny_app_en.R`](dic_shiny_app_en.R) | 10-tab interactive dashboard |
| [`dic_references_en.md`](dic_references_en.md) | 80 references (78 PubMed links) + equation↔reference correspondence table |

### The 49 compartments

Trigger · inflammation 7 (`PATH BLAST IL6 TNF NET HIST TF`) · endothelium 3
(`TM EPCR GLX`) · coagulation factors 10
(`FII FV FVII FVIII FX AT PC PS TFPI FIB`) · active enzymes 3
(`THR FXA APC`) · platelets 2 (`PLT PACT`) · fibrin 2 (`SFM FDEP`) ·
fibrinolysis 6 (`TPA PAI1 PLG PLN A2AP TAFIA`) · markers 1 (`DD`) ·
organs 4 (`OKID OLIV OLUN OCNS`) · outcomes 2 (`BLEEDC CUMH`) ·
drug PK 9 (`HEPC ENXD ENXC RTM APCX TXAC ATRAD ATRAC ARGA`)

### The 16 treatment scenarios

| | Scenario | 48 h platelets | 48 h fibrinogen | 48 h FDEP | **28-day mortality** |
|---|---|---|---|---|---|
| A | Normal control | 249 | 299 | 0.00 | 0.5% |
| B | Same infection, DIC-resistant host | 131 | 346 | 0.77 | 28.4% |
| **C** | **Septic DIC · supportive care only** | **78** | **245** | **1.29** | **43.3%** |
| D | + UFH 18 U/kg/h | 132 | 347 | 0.79 | 28.3% |
| E | + enoxaparin 1 mg/kg q12h | 97 | 287 | 1.09 | 36.3% |
| F | + AT concentrate (no heparin) | 99 | 292 | 1.03 | 36.0% |
| G | + AT concentrate + UFH | 158 | 383 | 0.55 | 23.0% |
| H | + thrombomodulin alfa (ART-123) | 84 | 261 | 1.28 | 43.8% |
| I | + drotrecogin alfa, early (6 h) | 105 | 302 | 0.92 | 37.4% |
| J | + drotrecogin alfa, delayed (30 h) | 93 | 279 | 0.94 | 35.8% |
| L | + tranexamic acid | 76 | 245 | 1.59 | 51.1% |
| M | + argatroban (direct thrombin inhibition) | 116 | 322 | 0.91 | 32.2% |
| **N** | **APL DIC · supportive care only** | **31** | **101** | **0.19** | **32.8%** |
| O | APL + ATRA | 33 | 108 | 0.19 | 7.5% |
| Q | APL + ATRA + fibrinogen + TXA | 30 | 338 | 1.04 | 21.8% |
| R | APL + ATRA + platelets | 54 | 105 | 0.20 | 6.7% |

Look at C and N side by side. **APL's deposited fibrin (0.19) is 15% of sepsis's
(1.29).** Same equations, same consumption process, and yet the stock does not
accumulate — because free t-PA is 13 times higher. That is why the APL patient dies
of bleeding rather than of organ failure.

---

## Verification

All 49 ODEs of the mrgsolve model were **independently re-implemented in
Python/scipy (LSODA)** and integrated over 28 days. That process exposed the
following defects, which were fixed.

1. **The steady state was not a steady state.** In the initial implementation the
   fibrinogen of a "normal" patient receiving no stimulus at all collapsed from
   300 → 25 mg/dL within 24 hours. The cause was that the pathological driving terms
   were responding to the baseline values themselves (`fTNF = TNF/(K+TNF)` is
   already 0.077 at TNF=10). Only after every driving term had been rewritten
   against **the increment above baseline** did the normal state stand still.
   Currently: no state variable drifts by more than 1.2% over 28 days (the largest
   drift is platelets 250→247).

2. **APC was destroying the zymogen that gets measured.** It had been written so
   that APC directly consumed the FV/FVIII **zymogen** compartments, so that PT shot
   up to 32 seconds when drotrecogin alfa was given. Real APC inactivates the
   activated **FVa/FVIIIa cofactors**. Switching to a structure that multiplies the
   prothrombinase/tenase terms by `APCF = 1/(1 + APC/KIAPC)` removed the PT runaway
   and at the same time gave APC a substantive anticoagulant effect, which flipped
   the sign of the rTM scenario.

3. **The active TAFIa half-life was set 14-fold too long.** It had been left at
   `KDTAFI = 0.30/h` (t½ 2.3 hours), whereas real TAFIa has a t½ of about
   10 minutes. In addition, because ART-123 is the extracellular D1-D2-D3 fragment
   it activates protein C well but TAFI poorly (efficiency coefficient `FTAFIS`),
   and this had been set too high, at 0.25. Until both values were corrected the
   model predicted that **thrombomodulin alfa increases mortality** — the opposite
   sign to SCARLET. After correction it is effectively neutral (+0.5%p), which
   matches SCARLET's negative result.

4. **The heparin-AT relationship had been written as linear.** As a result AT
   concentrate + heparin gave the unrealistic result of lowering mortality from
   43% → 1.3%. It came down to a realistic magnitude once the chemistry — that
   heparin is a **catalyst** and that the AT term saturates — had been reflected
   (see "Product (1)" above).

5. **The bleeding-risk metric had no dynamic range at the top.**
   `BIDX = 1 − Π(1−hᵢ)` saturates at 1, so it could not express "bleeding has been
   greatly reduced" in APL. Solved by moving the risk onto `−ln Π(1−hᵢ)` (an
   additive scale of the deficit).

---

## Where this model is wrong (an honest record)

| Item | Clinical trial | Model | Verdict |
|---|---|---|---|
| drotrecogin alfa, early | PROWESS −6.1%p | −5.9%p | agrees |
| AT concentrate, no heparin | KyberSept subgroup −5.8%p | −7.3%p | agrees |
| AT concentrate + heparin | KyberSept 0%p, bleeding 22.0 vs 12.8% | increment −5.3%p, bleeding increased | direction only |
| thrombomodulin alfa | SCARLET −2.6%p (p=0.32) | +0.5%p | agrees (both neutral) |
| ATRA in APL | early haemorrhagic death reduced by 2/3~3/4 | 32.8% → 7.5% | agrees |
| antifibrinolytics in APL | no benefit, thrombosis when combined with ATRA | 7.5% → 21.8% (harm) | agrees |
| **Unfractionated heparin** | **at best a few %p on meta-analysis** | **−15.0%p** | **disagrees** |

**The UFH term is this model's most exposed prediction.** It should be read as a
falsifiable hypothesis and must not be cited as a result. The mechanisms most
likely to be missing are (i) heparin-induced bleeding into already damaged organs
and (ii) the difficulty of aPTT/anti-Xa monitoring in DIC, which means that real
practice does not reach the target exposure.

Other structural limitations:

- **The cardiovascular item of SOFA is absent.** Haemodynamics and vasopressors are
  not modelled, so the simulated SOFA is up to 4 points lower than the real one.
- The coagulation cascade is **lumped**. FIX, FXI, FXII and kallikrein are on the
  map but are compressed into a single thrombin feedback term in the ODEs.
- `FDEP` is in **arbitrary units** and is not resolved by organ. The differences
  between kidney/lung/liver are expressed only through damage-sensitivity constants.
- It is a **deterministic** model with no inter-individual variability (IIV) and no
  residual error. Mortality should be read as a cohort mean, not as an individual's
  probability.
- Heparin-induced thrombocytopenia (HIT) is on the map but not in the ODEs.

---

## Usage

```r
# 1) load the model and run every scenario
source("dic_mrgsolve_model_en.R")
res <- run_all()
summarise_dic(res) |> dplyr::filter(time == 48)
summarise_dic(res) |> dplyr::filter(time == 672) |> dplyr::select(scenario, MORT)

# 2) dosing start-time sweep (testing the PROWESS-SHOCK hypothesis)
apc_timing_sweep()

# 3) interactive dashboard
shiny::runApp("dic_shiny_app_en.R")

# 4) re-render the map
# dot -Tsvg dic_qsp_model.dot -o dic_qsp_model.svg
# dot -Tpng -Gdpi=150 -Gsize="80,40" dic_qsp_model.dot -o dic_qsp_model.png
```

To turn the two clocks yourself, only four parameters need to be changed:

```r
# continuous deformation from the sepsis phenotype into the APL phenotype
run_dic("sepsis", par = list(SIL6 = 0.10,   # switch CLOCK 1 off
                             SPAI = 0.15,   # flip CLOCK 2
                             KANX = 3.0,    # annexin A2 amplification
                             MSUP = 0.20))  # myelosuppression
```

---

## Disclaimer

This model is for **research and educational purposes**. It must not be used for
the diagnosis or treatment decisions of an actual patient. The simulated mortality
figures are consequences of the model structure and parameters, not clinical
predictions, and must be read together with the discrepancies listed in "Where this
model is wrong" above.
