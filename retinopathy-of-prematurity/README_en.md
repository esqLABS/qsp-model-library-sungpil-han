# Retinopathy of Prematurity (ROP) — QSP Model
### Retinopathy of Prematurity · Quantitative Systems Pharmacology

<a href="rop_qsp_model.svg"><img src="rop_qsp_model.png" width="880" alt="ROP QSP mechanistic map"></a>

| Deliverable | File | Scale |
|---|---|---|
| Mechanistic map | [`rop_qsp_model.dot`](rop_qsp_model.dot) · [SVG](rop_qsp_model.svg) · [PNG](rop_qsp_model.png) | 142 nodes · 222 edges · 13 clusters |
| mrgsolve model | [`rop_mrgsolve_model.R`](rop_mrgsolve_model.R) | 36 ODEs · 9 scenario functions · virtual population |
| Shiny app | [`rop_shiny_app.R`](rop_shiny_app.R) | 8 tabs (all 17 outputs pass `testServer`) |
| References | [`rop_references.md`](rop_references.md) | over 110 verified PMIDs |

---

## The organising thesis

ROP is usually described as "the disease of too much oxygen and then too
little". That phrasing conceals the **two quantities** that actually matter, and
conceals the fact that the two are set by different people.

```
SUPPLY = PCHOR                 The choroidal oxygen tension supplied to the
                               avascular retina. The ventilator sets it. The
                               choroid does not autoregulate, so it follows
                               PaO2 almost linearly.

DEMAND = PCRIT = MO2·h² / 2k   The choroidal oxygen tension "required" to
                               oxygenate the full thickness of an avascular
                               retina of thickness h and consumption rate MO2.
                               Maturation sets it. No drug can move it.
```

Write the two separately and the whole disease is derived from their collision.
`PINNER = max(0, PCHOR − PCRIT)` is the oxygen tension of the inner avascular
retina, the day on which `PCHOR = PCRIT` is the **phase 1 → phase 2 transition
point**, and the **sign of every oxygen intervention is the sign of
(t − t_transition)**.

The second structural decision is that **there are two VEGF pools, not one**.
The two pools read opposite sides of that collision.

| Pool | What it reads | What it sets | Under hyperoxia |
|---|---|---|---|
| **VFRONT** (front pool) | the **local** pO2 at the growing front (choroid − the front's own consumption) | vessel growth (KG gate) · endothelial survival (KSURV gate) | **collapse** → growth arrest + vessel loss |
| **VEGFR** (vitreous pool) | hypoxic intensity × avascular **area** | neovascularisation (KNV gate, Hill n=4) · plus disease | **explosion** after the transition |

With only one VEGF variable this model does not hold together, because that
single variable would have to be suppressed by hyperoxia in phase 1 (growth
arrest) and raised by hyperoxia in phase 2 (increased NV). **Splitting the pool
is the only structural device that makes the sign reversal representable.**

And third, why phase 1 is a *window* is **derived** in the model. The pO2 at the
front is `0.9·PCHOR − MO2·(0.6h)²/2k`, so as the retina matures and its own
consumption rate rises the front **shields itself** from systemic hyperoxia.
Only while it is immature, thin and low-consuming is the front defenceless. The
existence of the window of vulnerability is not assumed; the consumption term
creates it.

---

## Quantitative results (every number was produced by the code in this repository)

### AXIS 1 — The pulse oximeter is not an oxygen sensor, and its error is larger than the experimental contrast of the clinical trials

With Hill n = 2.7, when P50 moves from 19.0 mmHg (HbF-dominant) to 26.6 mmHg
(HbA after transfusion), **the same SpO2 corresponds to a PaO2 differing by a
factor of 1.40 across the whole range.** Put the other way round, PaO2 50 mmHg
reads as

- **SpO2 93.2%** in an untransfused preterm infant
- **SpO2 84.6%** after replacement with adult blood

— **an offset of 8.6 points**. The entire contrast the NeOProM trials prescribed
(85–89% versus 91–95%) is 6 points. That is, **the haemoglobin switch that
proceeds uncontrolled during a NICU stay moves retinal oxygen more than the
intervention those trials tested.**

Worse still is the resolution. dPaO2/dSpO2 is 1.26 mmHg per point at SpO2 87%,
4.43 at 95%, 8.89 at 97%, and **45.45 mmHg** at 99%. The instrument loses its
resolution entirely in exactly the range where the therapeutic threshold of
AXIS 2 lies.

| SpO2 (%) | PaO2 (HbF) | PaO2 (HbA) | Ratio | ΔPaO2 per point |
|---|---|---|---|---|
| 85 | 36.1 | 50.6 | 1.40 | 1.05 |
| 89 | 41.2 | 57.7 | 1.40 | 1.56 |
| 93 | 49.5 | 69.3 | 1.40 | 2.83 |
| 95 | 56.5 | 79.2 | 1.40 | 4.43 |
| 97 | 68.8 | 96.4 | 1.40 | 8.89 |
| 99 | 104.2 | 145.9 | 1.40 | 45.45 |

### AXIS 2 — A square-root law puts a ceiling on every oxygen-based phase 2 therapy. The ceiling is anatomical

The avascular retina is a flat slab supplied from one face only. The depth of
oxygenation is `L = √(2k·PCHOR/MO2)` — **doubling the depth requires
quadrupling the choroidal oxygen tension**. Because MO2 and h both rise with
maturation, PCRIT soars from 6.7 mmHg at PMA 28 weeks to 93.7 mmHg at 36 weeks.

| PMA (weeks) | h (µm) | PCRIT (mmHg) | Required PaO2 | Required SpO2 (HbF) | Required SpO2 (HbA) |
|---|---|---|---|---|---|
| 30 | 123.5 | 13.9 | 16.3 | 39.9 | 21.1 |
| 32 | 137.1 | 33.0 | 38.9 | 87.4 | 73.6 |
| 34 | 152.9 | 64.6 | 76.1 | 97.7 | 94.5 |
| **35** | **160.2** | **80.4** | **94.6** | **98.7** | **96.9** |
| 36 | 166.5 | 93.7 | 110.2 | 99.1 | 97.9 |

STOP-ROP prescribed 96–99%. In HbF-dominant blood that band spans PaO2
**61.7–104.2 mmHg** — it **straddles** the requirement (94.6 at PMA 35 weeks).
That is this model's explanation for STOP-ROP having been neither a success nor
a failure but a near miss (48% → 41%, adjusted OR 0.72, 95% CI 0.52–1.01): only
**part of the upper end** of the prescribed range clears the diffusion
threshold, and the pulse oximeter cannot tell you which infant is inside it.

### AXIS 3 — Reproducing NeOProM requires an **achieved 2-point** separation, not the prescribed 6 points

One parameter (KNV) was set so that treatment-requiring ROP in the
higher-target arm came to 15.6% (observed 14.9%). Nothing else about ROP was
fitted. Sweeping the achieved separation then gives:

| Achieved separation (SpO2 points) | 1.0 | **2.0** | 3.0 | 4.0 | 6.5 |
|---|---|---|---|---|---|
| ROP RR | 0.872 | **0.743** | 0.674 | 0.604 | 0.471 |
| Death RR | 1.063 | **1.159** | 1.280 | 1.411 | 1.715 |

The NeOProM observations are ROP RR **0.74** and death RR **1.17**. The model
lands on both at 2.0 points, and **at that same separation predicts an absolute
incidence of 11.6% in the lower arm (observed 10.9%)** — a number to which
nothing was fitted.

At the prescribed 6.5 points the model says RR 0.47. That is, **the mechanistic
effect of the policy is about three times the effect the trials measured, and
the difference is exposure, not biology.** An achieved separation of ~2 points
is the value that comes out when time within the target range is 50–60% and the
distributions overlap. This is not a claim about ROP but a **falsifiable claim
about achieved oximetry**.

### AXIS 4 — The two standard anti-VEGF doses are equimolar. Therefore whatever differs between them is not "dose"

```
bevacizumab 0.625 mg / 149 kDa = 4.195 nmol
ranibizumab 0.200 mg /  48 kDa = 4.167 nmol      ratio = 1.007
```

What differs is (1) valency — 2 versus 1, a binding-capacity ratio per dose of
**2.01**, (2) vitreous half-life — **9.82 days versus 7.19 days** in the adult
eye, and (3) the only large difference, the **Fc**. Binding capacity ×
half-life favours bevacizumab by **2.75-fold** in the eye, and FcRn recycling
makes it **41.7-fold** more favourable (= unfavourable) systemically — matching
the serum AUC ratio of **35-fold** that Avery et al. measured in adults. One
molecule is slightly better in the eye and far worse systemically, and **both
facts come out of the same Fc.**

### AXIS 5 — The dose was set on adult terms, but of the two volumes involved only one is not on adult terms

Ocular effect scales with dose/vitreous volume, systemic risk with dose/body
weight. A preterm vitreous of 1.1 mL against an adult 4.0 mL; a preterm 1 kg
against an adult 70 kg. So from the same "half the adult dose" one gets

- ocular concentration: **1.8×** the adult exposure
- systemic dose per body weight: **35×** the adult exposure

And inside the vitreous the drug is in **roughly 229,000-fold molar excess**
over its target (3813 nM of binding sites at 0.625 mg against VEGF
1500 pg/mL = 0.033 nM). So **efficacy is logarithmic in dose** (one half-life
is the price of each halving) while **systemic exposure is strictly linear**.
The case for dose reduction is this asymmetry, not caution.

### AXIS 6 — The PEDIG dose floor is a stoichiometric limit, not an affinity limit, and the model puts it where it was found

PEDIG de-escalated over a 312-fold range and saw its first failure only at the
very bottom: 0.031 mg 9/9, 0.016 mg 13/13, 0.008 mg 9/9, **0.004 mg 9/10,
0.002 mg 17/23**. Computing the three candidate constraints for a 1.1 mL eye:

| Dose (mg) | Binding sites (pmol) | Residual concentration at 4 weeks (nM, preterm scaling) | Residual concentration at 4 weeks (nM, adult t½) | PEDIG 4-week success |
|---|---|---|---|---|
| 0.031 | 416 | 1.763 | 26.2 | 9/9 |
| 0.016 | 215 | 0.910 | 13.5 | 13/13 |
| 0.008 | 107 | 0.455 | 6.8 | 9/9 |
| **0.004** | **53.7** | **0.227** | 3.4 | **9/10** |
| **0.002** | **26.9** | **0.114** | 1.7 | **17/23** |

- **Affinity**: the residual concentration at the floor is 2–4× the reported Kd
  of 0.058 nM — it bites exactly there.
- **Stoichiometry**: with a vitreous VEGF half-life of 1 hour, cumulative VEGF
  production over 4 weeks is **17.1 pmol** (2.9 pmol if 6 hours). The margin at
  0.004 mg is 3.1-fold and **falls to 1.6-fold at 0.002 mg** — at exactly the
  dose that failed.
- **Duration**: being logarithmic, each halving costs only one half-life.

Because affinity and stoichiometry bite together to **within 2-fold at the
observed floor, this model cannot distinguish them, and it says so.** The
measurement that would distinguish them is the vitreous VEGF elimination rate
in the preterm eye. Only the stoichiometric reading yields a falsifiable
prediction: **in APROP the floor should rise** — because binding capacity is
consumed by VEGF production and APROP makes more of it.

### AXIS 7 — The serum VEGF suppression measured after IVB is 10-fold weaker than the 1:1 binding prediction, and the gap is the complex

Sato et al. measured, in infants given a total IVB dose of 0.5 mg, serum
bevacizumab of **1214 ng/mL** and serum free VEGF of **269 pg/mL** (baseline
1628, i.e. 16.5%). Naive 1:1 binding computed from the measured drug
concentration and the published Kd predicts **1.4% of baseline** — 12-fold too
much suppression.

The cause of the gap is not affinity. It is that **VEGF bound to drug is
protected from elimination, so total serum VEGF rises and the bound fraction
stays high**. Setting complex elimination 17-fold slower than that of free VEGF
reproduces **252 pg/mL** at day 14 (observed 269).

This produces a directly testable prediction that no published study reports:
**after IVB, total (bound + free) serum VEGF should be roughly an order of
magnitude above baseline, and free VEGF should be below baseline over the same
period.**

### AXIS 8 — A low dose is not a safety compromise. Because of the VEGF window it is the better strategy in the eye as well

Normal vascularisation requires VEGF above KG; neovascularisation requires it
above KNV. A dose that drags the shared vitreous pool **below** KG halts normal
vessels as well and leaves **residual avascular retina** — the substrate for
late reactivation. In a severe reference infant:

| Dose | Time with VEGF < KG | Time of reactivation | Residual avascular fraction | Final spherical equivalent |
|---|---|---|---|---|
| bevacizumab 0.625 mg | 6.0 weeks | +8.1 weeks | 0.414 | −1.78 D |
| bevacizumab 0.125 mg | 4.4 weeks | +6.4 weeks | 0.342 | −1.68 D |
| bevacizumab 0.031 mg | 3.0 weeks | +5.0 weeks | 0.281 | −1.58 D |
| bevacizumab 0.004 mg | 1.4 weeks | +3.1 weeks | 0.203 | −1.47 D |
| ranibizumab 0.2 mg | 4.7 weeks | +6.3 weeks | 0.325 | −1.64 D |
| aflibercept 0.4 mg | 11.9 weeks | +14.0 weeks | 0.550 | −1.86 D |
| laser (85% ablation) | 0 | none | 0.018 | **−8.42 D** |

Each 4-fold dose reduction loses about 1.5 weeks of quiescence and gains about
0.06 of retina that goes on vascularising. A high dose does not merely carry
systemic risk — it **freezes the retina in an avascular state that will need
the drug again.** aflibercept is the extreme of this axis: quiet the longest,
and leaving the most avascular retina.

### AXIS 9 — The competing risk is arithmetic, and it is decisive

From NeOProM's own absolute rates, lowering the target range gives, **per
1000 infants**:

| Item | Value |
|---|---|
| Excess deaths | **+28.0** |
| Reduction in ROP treatments | **−40.0** |
| Excess severe NEC | **+23.0** |
| Unfavourable structural outcomes averted (ETROP 9.1%) | 3.64 |
| Unfavourable visual outcomes averted (ETROP 14.5%) | 5.80 |
| **Deaths per unfavourable structural outcome averted** | **7.7** |
| **Deaths per unfavourable visual outcome averted** | **4.8** |

Not one parameter of this model enters that calculation. It is trial data
divided.

---

## Validation

### What was fitted (calibration) — 7

| Parameter | What it was fitted to | Source |
|---|---|---|
| KNV, KPLV | treatment-requiring ROP 14.9% in the higher-target arm | NeOProM [PMID 29872859](https://pubmed.ncbi.nlm.nih.gov/29872859/) |
| achieved separation 2.0 points | ROP RR 0.74 | same |
| H0, AHYPOX | death 17.1%, RR 1.17 | same |
| FSYS, CLSYS | day-14 serum 1214 ng/mL, t½ 21 days, tmax 14 days | Sato [21930258](https://pubmed.ncbi.nlm.nih.gov/21930258/) · Kong [25613938](https://pubmed.ncbi.nlm.nih.gov/25613938/) |
| KDVSB | day-14 serum free VEGF 269 pg/mL | Sato |
| CLSYS (ranibizumab) | bev/ranib serum AUC ratio 35 | Avery [25001321](https://pubmed.ncbi.nlm.nih.gov/25001321/) |
| KMYOP | zone I laser spherical equivalent −8.44 D | Geloneck [25103848](https://pubmed.ncbi.nlm.nih.gov/25103848/) |

**What was not fitted**: the Krogh diffusion constant, retinal oxygen
consumption rate, retinal thickness, P50/Hill coefficient, vitreous half-life,
binding affinity, and the ETROP type 1/2 definitions — all taken straight from
published values.

### Predictions that passed (held-out)

| Prediction | Model | Observed |
|---|---|---|
| incidence of treatment-requiring ROP in the lower-target arm | 11.6% | 10.9% (NeOProM) |
| refraction after bevacizumab 0.625 mg | −1.78 D | −1.51 D (Geloneck zone I) |
| the required PaO2 at PMA 35 weeks falls inside the STOP-ROP 96–99% band | 94.6 mmHg against a band of 61.7–104.2 | OR 0.72 (0.52–1.01), a near miss |
| direction and magnitude of the bev/ranib difference in serum VEGF suppression | 12.7% versus ~100% (relative to baseline) | IVB ≫ IVR, significant at 2 · 4 · 8 weeks (Wu) |
| dose reduction abolishes systemic VEGF suppression | 0.625 mg 59.2 days versus 0.031 mg **0 days** | the hypothesis offered as the rationale for dose reduction |
| the threshold for needing ROP treatment corresponds to the zone II/III boundary | AVASC_transition ≈ 0.43 (= radius 0.75) | zone III disease is almost never treated |
| the closing of the phase 1 window of vulnerability | the consumption term shields the front (after PMA ~31 weeks) | consistent with the descriptions of Ashton and Patz |

### Predictions that failed — reported, not fixed

**F1. It does not reproduce the laser recurrence rate of BEAT-ROP.** Observed:
recurrence before 54 weeks PMA, bevacizumab 4% versus laser 22%. Across the
whole range of ablation fractions from 95% down to 45% the model gives laser 0%
(residual avascular retina 0.008 → 0.059, never enough to cross KNV again) and
anti-VEGF 100% — **the ordering is inverted.** The diagnosis is structural and
specific: laser failure is **spatially local** (a missed sector remains locally
hypoxic), and a single homogeneous avascular compartment cannot represent
sectors. Reproducing 22% would require a sectorised avascular compartment, and
this model has none. In addition, the model's near-universal late anti-VEGF
reactivation is closer to the long-term follow-up literature than to BEAT-ROP's
54-week primary endpoint.

**F2. It does not reproduce Mega Donna Mega** (AA:DHA supplementation halved
severe ROP from 33.3% → 15.8%, aRR 0.50). The model gives RR 1.00. The reason
is quantitative and there is a number to state: the only pathway the model
gives AA/DHA is oxidative-stress-mediated vessel loss, and at the calibrated
ROS level that pathway accounts for **at most 5.8% of the vessel loss rate**.
Abolishing it entirely cannot produce the observed halving. So **this model is
evidence that the AA/DHA effect is not simple antioxidant protection**, and the
actual pathway (membrane phospholipid composition at the growing front,
resolvin/lipoxin mediators, or a direct action on angiogenic tone) is not here.

**F3. It overestimates the incidence of "ROP of any degree".** 99% of the
simulated cohort reach at least stage 1, where the reported figure is ~60–70%.
No infant in the model completes vascularisation before the metabolic
transition — on the mild side of the cohort, angiogenesis is too slow relative
to maturation. The **treatment-requiring endpoint** (the one that was
calibrated, and the one every drug axis rests on) is right (15.6% versus
14.9%).

**F4. The STOP-ROP subgroup finding is neither confirmed nor refuted.** The
trial saw benefit only in the absence of plus disease (46% → 32%) and none in
its presence (52% → 57%). At the model's enrolment time only 25 of 621 infants
have plus disease, and in the calibrated cohort the "threshold" surrogate is
effectively unreachable, so the subgroup contrast cannot be evaluated. On the
type 1 endpoint the model gives 23.2% → 13.6% (RR 0.59), which is a **larger**
phase 2 oxygen benefit than STOP-ROP's OR 0.72 (upper limit 1.01).

**F5. In the model, serum free VEGF after ranibizumab transiently exceeds
baseline** (120% at day 14), because at low drug concentrations complex
accumulation is faster than free VEGF removal. Avery et al. measured a slight
**decrease** (nadir 14.4 against a baseline of 17 pg/mL). The direction of the
bev/ranib contrast is right; the absolute ranibizumab level overshoots.

**F6. The achieved-separation inference of AXIS 3 depends on the model's
phase 1 oxygen sensitivity**, and that sensitivity was calibrated to ROP
incidence, not to a direct measurement of front-pool VEGF. An independent
measurement of retinal VEGF against PaO2 in phase 1 would test it, but no such
measurement exists for the human preterm eye.

**F7. The preterm vitreous volume (1.1 mL) and the scaling exponent (2/3) of
vitreous half-life on eye volume are assumptions, not measurements.** The floor
in AXIS 6 moves with both. Using the unscaled adult half-life raises the 4-week
residual concentration 15-fold and the affinity constraint no longer bites (the
last column of the table above).

### Defects found and fixed during integration (for the record)

1. **With a linear HIF→VEGF mapping the three VEGF regimes do not separate.**
   Trying to fit phase 1 suppression (~250 pg/mL) · physiological growth (~400)
   · phase 2 neovascularisation (>1170) simultaneously with a linear mapping
   produced the contradiction that the amplification coefficient would have to
   be negative. Resolved by introducing a fourth-power exponent reflecting the
   cooperativity of HRE occupancy.
2. **If vitreous VEGF reads hypoxic intensity alone, oxygen policy barely moves
   severity.** Before the area term (AVASC/AVREF) was added, peak NV at SpO2
   87% versus 93% was 0.383 versus 0.390 — effectively identical. VEGF is a
   concentration, and the steady-state concentration in a well-mixed
   compartment is proportional to production (= area × rate per unit area), so
   area has to enter.
3. **If NV carrying capacity is fixed independently of avascular area, the area
   dependence is destroyed.** `(1 − NV/NVMAX)` was replaced by
   `(AVASC·NVMAX − NV)` so that the avascular area at the transition point
   propagates directly into peak NV.
4. **With a single VEGF pool, late vascularisation stops.** With only the front
   pool and no consumption-shielding term, VFRONT fell to 86–96 pg/mL in room
   air (PaO2 ~100), so even a term infant failed to complete vascularisation
   (VASC 0.77 at PMA 44 weeks). Adding the front's own consumption normalised
   VASC to 0.94–0.97 and at the same time derived the closing of the phase 1
   window of vulnerability.
5. **The ridge indicator switched stage 2 on from immediately after birth.** A
   saturating expression already gave a steady state of 0.46 at phase 1 VEGF, so
   every infant was at stage 2 in the first week of life. Replacing it with a
   steep Hill function of the same form as NV (KSHV=880, n=4) moved the
   appearance of stage 2 to PMA 33.4–35.7 weeks — consistent with the natural
   history.
6. **The STOP-ROP axis was initially enrolled at the wrong time.** The stage≥2
   trigger fired at PMA 31.4 weeks (because of defect 5) and raised oxygen
   **before** the transition, and the model produced a result of exactly the
   opposite sign (the supplemented arm worse, 33.5% versus 27.3%). Once defect 5
   was fixed and enrolment moved to PMA 33.7 weeks, the direction flipped. An
   incident that showed the sign reversal is not an artefact of the model but
   **a function of timing**.
7. **The plus disease threshold was low and switched on in phase 1.**
   Individuals with a high VEGF set point already exceeded PLUS 0.5 in phase 1,
   so Type 1 fired at PMA 28.8 weeks (observed: 34–37 weeks). Resolved by
   raising KPLV from 700 → 1300 and n from 3 → 4.
8. **A 1000-fold unit error in the stoichiometric calculation** reported the
   4-week VEGF demand as 0.02 pmol (true value 17.1 pmol). Before this error was
   fixed the stoichiometric constraint appeared never to bite, and the
   conclusion of AXIS 6 was the opposite.
9. **Passing a list to `param()` as a positional argument made mrgsolve fail**,
   and an empty list raised an error in the no-treatment scenario — replaced
   with `update(mod, param = list(...))`.
10. **In the Python prototype the front block sat ahead of the MO2/HRET
    definitions**, so the consumption term was being computed as 0 (PFRONT
    exactly 0.9·PCHOR). Fixing the block order was what finally made an SpO2
    response appear.

---

## The 36 state variables

| Group | State |
|---|---|
| Oxygen transport | `FHBF` |
| Hypoxic signalling | `HIF` `VEGFR` `VFRONT` `EPOR` `ANG2` |
| Vascular geometry | `VASC` `ABLA` |
| Pathology | `SHUNT` `NV` `PLUS` `VSLOW` `TGFB` `FIBRO` `DETACH` |
| Systemic | `IGF1` `WT` `INFL` `BPD` `ROS` |
| Drug PK | `DVIT` `DAQ` `DCEN` `DPER` |
| Pharmacodynamics | `VEGFS` `AUCDS` `AUCVS` |
| Cumulative / outcomes | `HAZD` `HYPBUR` `HYPRBUR` `MYOP` `VFLOSS` `AVBUR` `NVPEAK` `VADEF` `RETPR` |

Derived outputs include `PCRIT` `PINNER` `LPEN` `PFRONT` `P50` `ZONE` `STAGE`
`PLUSD` `TYPE1` `THRESH` `PDEATH` — the ICROP3 zone/stage and ETROP type 1
definitions are implemented exactly as published.

## The 9 scenario functions

| Function | What it contrasts |
|---|---|
| `scenario_oximetry()` | the same SpO2, two haemoglobins (AXIS 1) |
| `scenario_diffusion_ceiling()` | required oxygen tension and penetration depth by PMA (AXIS 2) |
| `scenario_oxygen_sweep()` | the phase 1 dose-response of SpO2 targets 85–99% (AXIS 3) |
| `scenario_supplemental_o2()` | supplemental oxygen at prethreshold (the STOP-ROP design) |
| `scenario_dose_ladder()` | the 9 doses of the PEDIG de-escalation ladder (AXIS 5·6·7) |
| `scenario_treatment()` | laser versus bev versus ranib versus aflibercept (AXIS 8) |
| `scenario_transfusion()` | at a fixed SpO2 target, transfusion moves P50 |
| `scenario_igf1()` | nutrition · IGF-1 · rhIGF-1 supplementation (F-validation) |
| `cohort_rop()` | the virtual population against NeOProM |

---

## How to run it

```bash
# 1. Render the map
dot -Tsvg rop_qsp_model.dot -o rop_qsp_model.svg
dot -Tpng -Gdpi=150 rop_qsp_model.dot -o rop_qsp_model.png

# 2. Build the model + run all 9 scenarios
Rscript -e 'source("rop_mrgsolve_model.R"); run_all_rop(cohort_n = 400)'

# 3. Dashboard
Rscript -e 'shiny::runApp("rop_shiny_app.R")'
```

Required packages: `mrgsolve`, `dplyr`, `shiny` (graphics use base R only).

---

## ⚠️ Disclaimer

This is a **semi-quantitative QSP model for teaching and research purposes**. It
was assembled from the public literature and clinical trial data, but it has not
been independently verified or certified, and **must not be used directly for
real clinical decision-making, prescribing, or regulatory submission.** The
failure list F1–F7 above in particular is not decoration — it states what the
model does not reproduce, and the model should not be used for questions that
fall under those items.
