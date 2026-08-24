# Retinal Vein Occlusion (CRVO / BRVO) Macular Oedema — QSP Model

**Anti-VEGF is a multiplier in front of a bracket it cannot change.**

---

## 0. The one structural claim

Retinal thickening in vein occlusion is a Starling filtration flux:

```
    Jv  =  Lp · S · [ (Pc − Pt)  −  σ · (πc − πt) ]
           \______/   \________/    \___________________/
         PERMEABILITY  PRESSURE           ONCOTIC
             ARM         ARM               ARM
```

Of this model's **27 ODEs, every term in which anti-VEGF appears sits in the first
factor alone (`Lp`, `σ`).** Not one of them is inside the bracket. And the `Pc`
that is inside the bracket is set by the occlusion:

```
    Pc = (Pa · Rv + Pv · Ra) / (Ra + Rv),      Rv = Rv0 · (1 + OCC) / (1 + COLL)
```

Four results follow from that asymmetry **as arithmetic rather than as assertion**.

---

## 1. There is a critical capillary pressure

Set `Jv = 0` with the permeability arm forced to its floor — that is, with a
hypothetical infinite dose of a perfect blocker — and solve for the capillary
pressure:

```
    Pc*  =  Pt + σmax · (πc − πt0)  =  12 + 0.92 × 21  =  31.32 mmHg
```

Above `Pc*` the bracket stays positive with `Lp` already basal, so **no dose of
any anti-VEGF agent can dry that macula.**

The model actually performs that experiment — infinite blockade of VEGF · PlGF ·
Ang-2 switched on at presentation, and whatever oedema survives is the pressure
arm by construction:

| Phenotype | Pc (mmHg) | Jv,floor (µm/d) | CST 6 months | CST 36 months |
|---|---|---|---|---|
| Mild occlusion | 20.9 | −9.09 | 250 | 250 |
| BRVO | 26.2 | −3.76 | 251 | 250 |
| Non-ischaemic CRVO | 28.9 | −2.14 | 254 | 250 |
| **Ischaemic CRVO** | **37.2** | **+5.09** | **291** | **286** |

Only ischaemic CRVO has `Jv,floor > 0` — that is, 36 µm of oedema remains even
with a perfect drug.

### The critical venous resistance

Solving `Pc = Pc*` for `Rv` gives `Rv* = (Pc* − Pv)·Ra / (Pa − Pc*)`, which
diverges as `Pa → Pc*`. **If the inflow arterial pressure is below `Pc*`, no
occlusion whatsoever can create a permanent oedema floor.**

| Pa (mmHg) | Rv* | vs normal (0.25) | dPc/dPa |
|---|---|---|---|
| 48.0 | 0.859 | 3.4× | 0.462 |
| 44.0 | 1.129 | 4.5× | 0.530 |
| **40.0** | **1.650** | **6.6×** | **0.623** |
| 36.0 | 3.060 | 12.2× | 0.754 |
| 34.0 | 5.343 | 21.4× | 0.842 |
| **32.0** | **21.059** | **84.2×** | **0.955** |
| 31.5 | 79.556 | 318× | 0.988 |
| 31.0 | Unreachable | — | — |

**Lower the inflow pressure by 8 mmHg, from 40 to 32 mmHg, and the venous
resistance needed before a permanent oedema floor can exist rises 12.8-fold, from
1.65 to 21.1.** This is non-linear because `dPc/dPa = Rv/(Ra+Rv)` approaches 1 as
the occlusion tightens — in a normal eye only 20% of the arterial pressure is
transmitted to the capillary, whereas in an obstructed eye it arrives almost
undamped.

Hypoxic arteriolar dilatation (a 35% fall in `Ra`) lowers that threshold
**further**:

| Ra | Pa = 40 | Pa = 36 | Pa = 32 |
|---|---|---|---|
| 1.00 (normal) | Rv* = 1.650 | 3.060 | 21.059 |
| 0.65 (hypoxic dilatation) | Rv* = 1.072 (4.3×) | 1.989 (8.0×) | 13.688 (54.8×) |

Chronicity lowers the threshold too. Once the barrier can no longer be fully
resealed (`TJmax = 0.70` → σ = 0.756, πt = 6.81):

```
    Pc*  :  31.32  →  25.76 mmHg
    Rv*  :   1.650 →   0.615   (2.68× lower)
```

**That is, an occlusion that could have been dried at presentation becomes an
occlusion that cannot be dried after two years of untreated oedema.** Chronicity
turns the eye into a floor eye.

### When the eye crosses the critical pressure

Comparing `Pc` day by day against `Pc* = 31.32 mmHg` in the untreated natural
history:

| Phenotype | Pc (day 42) | Day it falls below Pc* | Days above Pc* |
|---|---|---|---|
| Mild occlusion | 28.6 | Day 42 (from the outset) | 0 |
| BRVO | 35.3 | **Day 153** | 111 |
| Non-ischaemic CRVO | 33.7 | **Day 181** | 139 |
| Ischaemic CRVO | 37.9 | **Never** | 1054 |

BRVO and non-ischaemic CRVO drop below the critical pressure at around six months
as the collateral bypass matures; ischaemic CRVO never drops below it in three
years. **What divides the two clinical courses is this crossing time, not the
choice of agent.** In the model the difference between BRVO and CRVO is expressed
by a single `COLL_MAX` (2.20 vs 0.70) and a single non-recanalisable residual
`OCC_RES` (0.32 vs 0.55), and that alone produces the difference in prognosis.

---

## 2. The pressure arm amplifies whatever permeability remains

Because `Jv = Lp × bracket`, **the same residual rise in permeability makes more
water when the bracket is large.** How differently the gap between real
aflibercept treat-and-extend and infinite blockade — that is, exactly the same
pharmacological shortfall — appears in two eyes:

| | Real T&E (36 months) | Infinite blockade (36 months) | Difference |
|---|---|---|---|
| Non-ischaemic CRVO (small bracket) | 257 µm | 250 µm | **7.3 µm** |
| Ischaemic CRVO (large bracket) | 363 µm | 286 µm | **77.2 µm** |

**The same pharmacological shortfall is amplified 10.6-fold.** This is why "let us
try a different agent" fails in a high-pressure eye — the shortfall being
amplified is itself already small. What switching agent can recover is 7 µm; what
remains is the pressure arm.

---

## 3. Suppression duration and dry duration are different numbers

Suppression duration is **pharmacology**. With `Cret(t) = C0,eff·e^(−kt)` and a
suppression ratio `R = 1 + Cret/(α·KD)`:

```
    t_sup = (t½ / ln2) · ln[ C0,eff / ((R−1) · α · KD) ]
```

| Agent | Molar dose (nmol) | C0 vitreous (µM) | KD (pM) | t½ (d) | **t_sup (d)** |
|---|---|---|---|---|---|
| Ranibizumab 0.5 mg | 10.4 | 2.60 | 46.00 | 7.19 | **23.5** |
| Bevacizumab 1.25 mg | 8.4 | 2.10 | 58.00 | 9.82 | **25.8** |
| Brolucizumab 6 mg | 230.8 | 57.69 | 28.40 | 5.10 | **43.0** |
| Faricimab 6 mg | 40.3 | 10.07 | 30.00 | 7.50 | **43.8** |
| Aflibercept 2 mg | 17.4 | 4.35 | 0.49 | 9.10 | **96.2** |
| Aflibercept 8 mg | 69.6 | 17.39 | 0.49 | 9.10 | **114.4** |

### Why a 94-fold difference in affinity becomes a 4-fold difference in time

Decompose `ln(C0,eff / threshold)` into a **reservoir term** and an **affinity
term** and multiply by the half-life factor:

| Agent | ln(C0,eff) | −ln(threshold) | Sum | t½/ln2 | t_sup |
|---|---|---|---|---|---|
| Ranibizumab | 12.25 | −9.98 | 2.27 | 10.37 | 23.5 |
| Bevacizumab | 12.03 | −10.21 | 1.82 | 14.17 | 25.8 |
| Brolucizumab | 15.34 | −9.49 | 5.85 | 7.36 | 43.0 |
| Faricimab | 13.60 | −9.55 | 4.05 | 10.82 | 43.8 |
| Aflibercept | 12.76 | −5.44 | 7.32 | 13.13 | 96.2 |

Duration is **logarithmic** in dose and affinity and **linear** in half-life. The
intrinsic affinities of aflibercept and ranibizumab differ **94-fold**, yet their
suppression durations differ only **4.08-fold**. And raising the aflibercept dose
four-fold (2 → 8 mg) lengthens the duration from 96 to 114 days — only **19%**.

> This is why SCORE2 reported bevacizumab (KD 58 pM) as the equal of aflibercept
> (KD 0.49 pM) (+18.6 vs +18.9 letters, a 118-fold difference in affinity). A model
> in which effect is proportional to affinity is refuted by SCORE2; what survives
> is the logarithmic structure of `t_sup`.

### And dry duration is **the disease**

Tally "days of adequate VEGF suppression" and "days with CST < 310 µm" separately
within the same simulation:

| Scenario | Injections | Days suppressed | Days dry | Ratio |
|---|---|---|---|---|
| Ranibizumab PRN | 23 | 541 | 555 | 1.03 |
| Aflibercept T&E | 15 | 974 | 1033 | 1.06 |
| Faricimab T&E | 21 | 769 | 873 | 1.13 |
| **Ischaemic CRVO, aflibercept q4w** | **38** | **1053** | **4** | **0.00** |

**The ischaemic CRVO eye had VEGF adequately suppressed on 1053 of 1095 days, and
the macula was dry on four.** 38 injections, 3 years, VEGF held down throughout,
oedema present throughout. In the clinic this eye is classified as an "anti-VEGF
non-responder"; in the model it is not a non-responder but a **pressure-arm eye.**
Because `t_sup` is already satisfied, there is nothing to recover by switching
agent or shortening the interval.

**The diagnostic rule:** if it is suppressed and yet wet (ratio < 1), it is the
pressure arm.

---

## 4. Oedema is a state, vision is an integral

Visual loss is the sum of two terms:

- `L_ed = 36 · Wtot/(230 + Wtot)` — **a function of the present state.** It comes
  back when the macula dries.
- `L_pr = 42 · (1 − EZ)` — **the time integral of thickness × ischaemia.** It does
  not come back.

The price of delay is therefore billed **permanently**, whether or not the oedema
eventually resolves.

| Delay (months) | Starting BCVA | Peak BCVA | Mean BCVA (24–36 months) | EZ (36 months) | Permanent loss |
|---|---|---|---|---|---|
| 0 | 57.6 | 80.0 | 73.4 | 0.898 | — |
| 0.5 | 56.4 | 80.0 | 72.9 | 0.888 | −0.5 |
| 1 | 55.9 | 80.0 | 72.5 | 0.878 | −0.9 |
| 2 | 55.3 | 80.0 | 71.6 | 0.859 | −1.8 |
| 3 | 54.7 | 80.0 | 71.2 | 0.846 | −2.2 |
| 4 | 54.2 | 80.0 | 70.5 | 0.834 | −2.8 |
| **6** | 53.2 | 80.0 | **68.8** | 0.802 | **−4.6** |
| 9 | 51.7 | 80.0 | 66.8 | 0.758 | −6.6 |
| **12** | 50.2 | 80.0 | **64.9** | 0.722 | **−8.5** |

Worth noting: **peak BCVA is the same, 80.0, in every delay arm.** No arm "fails
to respond" — the oedema clears well in all of them. The difference is what is
left behind once the oedema has gone. This is the structure behind COPERNICUS,
where the immediate-treatment arm was at +13.0 letters at week 100 while the
deferred arm stayed at +1.5 letters.

---

## 5. Trial replication

Non-ischaemic CRVO, baseline = day of first injection (day 42). Baseline CST
613 µm, baseline BCVA 57.6 letters.

| Scenario | Injections | CST 6 months | ΔCST | BCVA 6 months | 12 months | 24 months | Observed target |
|---|---|---|---|---|---|---|---|
| No treatment | 0 | 499 | −114 | **−4.4** | −7.3 | −12.4 | CRUISE sham +0.8 |
| Ranibizumab 0.5 q4w×6 → PRN | 23 | 282 | −331 | **+12.3** | +4.0 | +5.3 | CRUISE +14.9 |
| Bevacizumab 1.25 q4w×6 → PRN | 22 | 284 | −329 | **+12.5** | +10.4 | +7.9 | SCORE2 +18.6 |
| **Aflibercept 2 mg q4w×6 → T&E** | 15 | 261 | −352 | **+17.8** | +12.1 | +17.4 | **COPERNICUS +17.3 / GALILEO +18.0** |
| Faricimab 6 mg q4w×6 → T&E | 21 | 262 | −351 | **+17.3** | +16.3 | +7.6 | COMINO +16.9 |
| Brolucizumab 6 mg q4w×6 → T&E | 26 | 266 | −348 | +16.8 | +15.8 | +15.2 | — |
| Aflibercept 8 mg q4w×3 → T&E | 12 | 261 | −352 | +17.5 | +17.7 | +14.4 | — |

BRVO (baseline CST 590 µm, BCVA 58.0):

| Scenario | Injections | CST 6 months | BCVA 6 months | 24 months | Observed target |
|---|---|---|---|---|---|
| No treatment | 0 | 378 | +5.3 | +11.2 | BVOS — spontaneous improvement occurs |
| Aflibercept q4w×6 → PRN | 7 | 255 | **+19.4** | +14.3 | **VIBRANT +17.0** |
| Ranibizumab q4w×6 → PRN | 8 | 314 | +13.9 | +13.3 | BRAVO +18.3 |

The aflibercept six-month result agrees with the two phase 3 trials to within
0.5 letters. Ranibizumab and bevacizumab come out 2–6 letters below the observed
values (see the limitations in §9 below).

---

## 6. Ischaemic CRVO, neovascular glaucoma, and rebound

| Scenario | Injections | NP (DA) | NVI | IOP | NVG | CST 36 months | BCVA 36 months |
|---|---|---|---|---|---|---|---|
| No treatment | 0 | 55.6 | 0.838 | 32.5 | **Yes** | 1130 | **0.0** |
| Aflibercept q4w continuously for 36 months | 38 | 38.1 | 0.002 | 15.0 | No | 363 | 24.0 |
| **Aflibercept q4w ×6, then stopped** | **6** | 54.3 | **0.837** | **32.5** | **Yes** | 1122 | **0.0** |
| Aflibercept ×6 + PRP (at 3 months), then stopped | 6 | 53.2 | 0.793 | 30.7 | Yes | 1091 | 0.0 |
| Aflibercept T&E | 38 | 38.1 | 0.002 | 15.0 | No | 363 | 24.0 |

Two things come out.

**(1) Anti-VEGF removes only the *signal* drive for neovascularisation.** In the
model `dNVI/dt ∝ f(NP) × g(VACT)`, so while drug is present the second factor is
near zero and NVI regresses, but the first factor — non-perfusion — remains as it
was. Stop, and NVI returns to 0.837 and neovascular glaucoma develops — an outcome
entirely different from the arm that was not stopped. **Stopping the drug while
non-perfusion remains is not ending treatment; it is switching the clock back on.**

**(2) In a pressure-arm eye, treat-and-extend collapses to the minimum interval.**
`ISCH_AFL_TAE` receives exactly the same 38 injections as `ISCH_AFL_CONT` — the eye
is wet at every visit, so the T&E rule keeps shortening the interval down to its
floor (4 weeks). In this eye T&E is not a means of economy but another name for
fixed q4w.

---

## 7. Steroids, blood pressure, and real-world care

### Steroids reach mediators anti-VEGF cannot — and there is a price

| Scenario | Implants | Injections | CST 6 months | BCVA 6 months | Cataract 36 months | IOP 36 months | BCVA 36 months |
|---|---|---|---|---|---|---|---|
| Dexamethasone q6mo, phakic | 6 | 0 | 487 | −6.9 | **0.828** | 19.6 | 29.4 |
| Dexamethasone q6mo, pseudophakic | 6 | 0 | 487 | −3.9 | 0.084 | 19.6 | **39.3** |
| Aflibercept T&E | 0 | 15 | 261 | +17.8 | 0.084 | 15.0 | **73.8** |
| Aflibercept T&E + dexamethasone ×2 | 2 | 15 | 266 | +17.3 | 0.084 | 15.0 | 74.5 |

**The same implant schedule is 9.9 letters worse in a phakic eye than in a
pseudophakic one** — and the whole difference is the lens (cataract 0.828 vs
0.084). Steroids suppress IL-6 (−80%), Ang-2 (−50%) and VEGF transcription (−70%)
simultaneously and raise junction-protein transcription directly (the model's
`E_DEX_TJ`), but with this model's parameters they are clearly inferior to
anti-VEGF alone, and the incremental benefit of combining them is small
(+0.7 letters).

### Blood pressure is one of the few arrows pointing left, into the pressure arm

| Scenario | Pc 36 months | CST 36 months | BCVA 36 months | Injections |
|---|---|---|---|---|
| Aflibercept T&E | 28.9 | 257 | 73.8 | 15 |
| **Aflibercept T&E + Pa −8 mmHg** | 24.7 | 250 | **77.4** | 15 |
| Pa −8 mmHg alone (no drug) | 24.7 | 296 | 61.5 | 0 |
| No treatment | 28.9 | 432 | 40.4 | 0 |
| Ischaemic CRVO, aflibercept T&E | 37.2 | 363 | 24.0 | 38 |
| **Ischaemic CRVO, aflibercept T&E + Pa −8** | 30.2 | **269** | **50.7** | **19** |

In a non-ischaemic eye, 8 mmHg is worth +3.5 letters. **In an ischaemic eye it is
worth +26.7 letters and the injection count halves, 38 → 19.** Touch the pressure
arm in an eye the pressure arm dominates and that much comes back. (The model's
`Pa` is retinal arteriolar inflow pressure and is not 1:1 with systemic MAP — what
clinical fall in blood pressure an 8 mmHg fall in inflow pressure corresponds to
is not something this model answers. See §9.)

### Under-treatment does not reduce the benefit; it reverses it

| Scenario | Injections | BCVA change (36 months) | CST 36 months | EZ | Letters per injection |
|---|---|---|---|---|---|
| Registry-like (6+3+2 = 11) | 11 | **+4.2** | 363 | 0.842 | +0.38 |
| Severe under-treatment (4+2+1 = 7) | 7 | **−4.8** | 390 | 0.696 | −0.68 |
| Protocol T&E | 15 | **+16.3** | 257 | 0.898 | +1.08 |

Cutting the injection count from 15 to 11 does not reduce the benefit
proportionally — it falls to a quarter, and at 7 injections it turns negative. The
reason is §4: oedema re-accumulates in every interval, and the ellipsoid zone
integrates **all** of those intervals.

24-month dose-response (aflibercept, cap on total injections):

| Cap | Actual | CST 24 months | BCVA gain | Letters per injection |
|---|---|---|---|---|
| 0 | 0 | 432 | −12.4 | — |
| 4 | 4 | 409 | −6.7 | −1.68 |
| 6 | 6 | 404 | −5.2 | −0.87 |
| 8 | 8 | 396 | −2.3 | −0.29 |
| 10 | 10 | 386 | +2.3 | +0.23 |
| **12** | **12** | **254** | **+17.4** | **+1.45** |
| 16 | 12 | 254 | +17.4 | +1.45 |
| 20 | 12 | 254 | +17.4 | +1.45 |

**The dose-response is not smooth; it has a threshold.** Two more injections, from
10 to 12, turn +2.3 letters into +17.4. Beyond 12 the T&E rule itself stops at 12,
so no more can be given — that is, **there is a place where "enough" sits, below
which almost all of it is wasted and above which none of it is needed.**

---

## 8. The occlusion-severity response surface and the virtual population

### Occlusion severity orders everything (aflibercept T&E, 24 months)

| OCC0 | Rv | Pc | Baseline CST | CST 24 months | Dry? | Injections | BCVA gain |
|---|---|---|---|---|---|---|---|
| 2.0 | 0.32 | 22.6 | 297 | 250 | Yes | 11 | +5.7 |
| 5.0 | 0.57 | 25.6 | 443 | 253 | Yes | 11 | +14.9 |
| 9.0 | 0.90 | 28.5 | 587 | 252 | Yes | 12 | +15.9 |
| 12.0 | 1.15 | 30.0 | 657 | 258 | Yes | 12 | **+16.7** |
| 16.0 | 1.48 | 31.6 | 719 | 268 | Yes | 12 | +14.6 |
| 22.0 | 1.98 | 33.3 | 779 | 282 | Yes | 14 | +9.7 |
| 30.0 | 2.65 | 35.0 | 825 | 302 | Yes | 21 | +3.9 |
| 45.0 | 3.90 | 36.7 | 872 | 346 | **No** | 25 | **−14.4** |

The gain is **non-monotonic**. It is maximal near OCC0 = 12 (+16.7) and falls away
on both sides — the mild eye because it has no vision to lose, the heavy eye
because it cannot get it back. And the injection count rises monotonically with
severity (11 → 25). **Severity raises the benefit up to a point, and past that
point it turns into cost.**

### Virtual population (n = 400, 24 months, aflibercept treat-and-extend)

Occlusion severity (`OCC0` log-normal, median 11, CV 0.55, truncated at 40), the
non-recanalisable residual, inflow pressure (normal, 40 ± 4.5 mmHg), the
collateral ceiling, photoreceptor vulnerability, the VEGF induction gain,
presentation delay (normal, 42 ± 25 days) and the number of injections actually
received (normal, 11 ± 4, truncated to 4–24) were each drawn at random.

| Metric | Value |
|---|---|
| Mean baseline CST · BCVA | 533 µm · 61.1 letters |
| Mean 24-month CST · change in BCVA | 361 µm · **+4.0 letters** |
| Mean number of injections | 9.6 |
| Dry macula (CST < 310 µm) | 65.5% |
| **Pressure-arm floor present (Jv,floor > 0)** | **19.2%** |
| Floor present and still wet | 14.8% |
| No floor and still wet (a genuine pharmacological shortfall) | 19.8% |
| ≥15 letters gained | 18.2% |
| ≥15 letters lost | 8.5% |
| Iris neovascularisation (NVI > 0.5) | 5.0% |

Stratified by whether a pressure-arm floor is present:

| Stratum | n | Pc | CST 24 months | BCVA change | Injections | EZ |
|---|---|---|---|---|---|---|
| **Floor present** | 77 | 34.8 | **650 µm** | **−13.2** | 10.5 | 0.617 |
| No floor | 323 | 26.0 | 292 µm | **+8.1** | 9.3 | 0.894 |

**The two strata received almost the same number of injections (10.5 vs 9.3).**
The difference is not dose but pressure. As correlation coefficients:

```
    Pearson r (24-month Pc      vs  24-month residual CST)  =  +0.816
    Pearson r (injection count  vs  24-month residual CST)  =  −0.076
```

**Capillary hydrostatic pressure explains residual thickness overwhelmingly better
than the number of doses given.** In this population the injection count explains
essentially none of the residual oedema (r² = 0.006 vs 0.666).

Stratified by injection count, the threshold reappears:

| Injections | n | BCVA change | CST | EZ |
|---|---|---|---|---|
| 0–4 | 30 | −1.4 | 406 | 0.778 |
| 5–6 | 30 | +0.3 | 386 | 0.782 |
| 7–8 | 62 | −1.1 | 395 | 0.781 |
| **9–12** | **256** | **+7.0** | **329** | **0.882** |

### A methodological warning: "gain from baseline" rewards late presentation

Stratify by presentation delay and, on the face of it, the patients who came late
gain more (0–25 days −5.2 letters, 45–70 days +7.2 letters). This is a
**measurement artefact** — `gain` is computed against the BCVA on the day of first
injection, so a patient who came early is compared against a baseline that has not
yet finished deteriorating. The columns genuinely sensitive to delay are 24-month
absolute BCVA and EZ:

| Presentation delay | n | Apparent gain | **24-month absolute BCVA** | **EZ** |
|---|---|---|---|---|
| 0–25 days | 86 | −5.2 | **64.5** | 0.853 |
| 25–45 days | 97 | +6.4 | **66.4** | 0.855 |
| 45–70 days | 154 | +7.2 | 65.9 | 0.842 |
| Beyond 70 days | 63 | +4.9 | **61.7** | **0.799** |

Read in absolute BCVA, presentation beyond 70 days is the worst. **That a trial
which enrolled patients with worse baselines reports a larger letter gain is
exactly the same arithmetic** — the size of the gain is in part a choice of
baseline.

---

## 9. What this model does not claim (limitations, stated plainly)

1. **`ALPHA = 120` is a calibration constant, not a mechanism.** The 120-fold gap
   between SPR-measured affinity and the apparent in vivo IC50 lumps competition
   with VEGFR2, local interstitial flux and retinal diffusion limitation into a
   single number. It is applied **identically** to every agent, so it cannot
   manufacture differences between agents, but it is not a measured quantity.
   Absolute suppression durations (96 days and so on) should be read only as
   orders of magnitude.
2. **Ranibizumab and bevacizumab come out 2–6 letters below the observed values.**
   The `t_sup` of both agents (23.5 and 25.8 days) is shorter than the monthly
   interval, so in the model some oedema returns at the end of every cycle. The
   real CRUISE and SCORE2 did not behave that way. That is, `ALPHA` may be
   penalising low-affinity agents excessively, which is a criticism pointing in the
   same direction as what SCORE2 says.
3. **`Pa` is not systemic blood pressure.** It is retinal arteriolar inflow
   pressure, and this model does not answer what clinical fall in blood pressure a
   change of 8 mmHg corresponds to. The blood-pressure result in §7 is a structural
   statement — "touch the pressure arm and this much comes back" — not a clinical
   recommendation that "this much antihypertensive buys this much improvement".
4. **`σ = σmax · TJ^0.55` is a modelling choice.** The functional coupling between
   the reflection coefficient and junction integrity was chosen so as to give
   qualitatively correct behaviour; it is not a measured relationship. It matters
   because it is precisely this relationship that sets `Pc*` — if this exponent is
   wrong, every threshold in §1 moves.
5. **The absolute thicknesses are calibrated values.** `KF`, `VMAX_PUMP` and
   `KM_PUMP` were tuned so that the CRVO baseline CST comes out near the trial
   baselines. What the model produces is **the difference between arms**; the
   absolute µm are calibration.
6. **The untreated natural history is more pessimistic than observed.** The model's
   untreated non-ischaemic CRVO is at −4.4 letters at six months whereas CRUISE
   sham was +0.8. Lowering `OCC_RES` (the non-recanalisable residual) would fit
   that, but then the model fails to reproduce RETAIN's report at 3–4 years that
   "half still require injections". This model chose to fit the long-term end.
7. **These are not individual patient predictions.** The parameters are
   population-level and the calibration targets are trial means. Nothing was fitted
   to any individual eye.
8. **`Pc` cannot be measured.** There is no clinical way of measuring capillary
   hydrostatic pressure. That is exactly why the diagnostic rule of §3 ("if it is
   suppressed and yet wet, it is the pressure arm") is the usable form — you do not
   need to know `Pc`, and you do know the drug concentration and the OCT.

---

## 10. Files

| File | Contents |
|---|---|
| [`rvo_qsp_model.dot`](rvo_qsp_model.dot) | Mechanistic map source — 25 clusters, 130+ nodes. The centre (cluster 13) is the Starling node, the left half is the pressure arm and the right half the permeability arm. Every anti-VEGF arrow is on the right. |
| [`rvo_qsp_model.svg`](rvo_qsp_model.svg) · [`rvo_qsp_model.png`](rvo_qsp_model.png) | Renderings (SVG / 150 dpi PNG) |
| [`rvo_mrgsolve_model.R`](rvo_mrgsolve_model.R) | mrgsolve model (27 ODEs) + agent library + phenotypes + PRN/T&E controllers + 31 scenarios + closed-form analysis + virtual population + one plot per claim (4) |
| [`rvo_reference_model.py`](rvo_reference_model.py) | **The independent reference implementation that was actually executed** (pure standard library, RK4). Every term of the R file re-implemented term by term. Every number in this README comes from here. |
| [`rvo_reference_output.txt`](rvo_reference_output.txt) | The full execution log of that file |
| [`rvo_scenario_results.json`](../../../retinal-vein-occlusion/rvo_scenario_results.json) | 31 scenarios + suppression durations + delay/dose/severity curves (machine-readable) |
| [`rvo_population_results.json`](../../../retinal-vein-occlusion/rvo_population_results.json) | Virtual population summary for 400 patients + 120 individual rows |
| [`rvo_shiny_app.R`](rvo_shiny_app.R) | 10-tab Shiny dashboard. Of the sliders, only `Pa` moves the pressure arm; every other treatment control moves the permeability arm — and a verdict banner states in words which of "pressure-arm eye / pharmacological shortfall / dry" applies. |
| [`rvo_references.md`](rvo_references.md) | 112 references, every PMID verified against PubMed E-utilities. Each entry states what it was used for in the model (parameter · structural assumption · calibration target). |

### How to run

```bash
# render the map
dot -Tsvg rvo_qsp_model.dot -o rvo_qsp_model.svg
dot -Tpng -Gdpi=150 rvo_qsp_model.dot -o rvo_qsp_model.png

# run the reference implementation (no R needed, pure Python; about 25 min)
python3 rvo_reference_model.py

# the mrgsolve model (R required)
Rscript -e 'source("rvo_mrgsolve_model.R"); res <- run_all(); print(analyse_critical_pressure())'

# the Shiny dashboard
Rscript -e 'shiny::runApp("rvo_shiny_app.R", port = 8080)'
```

### State vector (27 ODEs)

| # | State | Unit | Meaning |
|---|---|---|---|
| 1–2 | `DVIT` `ASYS` | nmol | Anti-VEGF vitreous reservoir / systemic |
| 3–6 | `VTONE` `PTONE` `ATONE` `IL6` | pM, pg/mL | VEGF-A · PlGF · Ang-2 · IL-6 **tone** (the value in the absence of drug) |
| 7 | `HIF` | 0–1 | Hypoxic signal |
| 8–9 | `OCC` `COLL` | — | Occlusion severity / collateral conductance |
| 10–11 | `NP` `MI` | DA, 0–3 | Non-perfused area / macular ischaemia index |
| 12 | `LEUK` | 0–1 | Leukostasis |
| 13–14 | `TJ` `CHRON` | 0–1 | Junction integrity / chronic remodelling |
| 15–16 | `W` `SRF` | µm | Intraretinal excess water / subretinal fluid |
| 17–19 | `EZ` `DRIL` `CUMED` | 0–1, µm·d | Ellipsoid zone · disorganisation of the inner layers · cumulative oedema exposure |
| 20–22 | `NVI` `IOP` `CAT` | 0–1, mmHg, 0–1 | Iris/angle neovascularisation · intraocular pressure · lens |
| 23–24 | `IMP` `CDEX` | ng, ng/mL | Dexamethasone implant remaining / vitreous concentration |
| 25–27 | `BCVAO` `ABLA` `TSUP` | letters, 0–1, d | Observed acuity (10-day lag) · PRP ablated fraction · cumulative days of suppression |

---

## ⚠️ Disclaimer

This is a qualitative and semi-quantitative QSP model for educational and research
purposes. It was assembled from the published literature and clinical trial data
but has not been independently validated or certified, and **must not be used for
real clinical decisions, prescribing, or regulatory submission.** The
blood-pressure results in §7 in particular are a structural statement about the
model, not a clinical recommendation.
