# Refactor notes — `igg4-related-disease/igg4rd_mrgsolve_model.R`

Scope: **all three modelled compounds** — Rituximab (RTX), Prednisone/prednisolone
(PRED), and Dupilumab (DUP). There is no fourth compound left untouched in this
file; every compartment, parameter, and effect equation in
`igg4rd_mrgsolve_model.R` belongs to one of these three PK blocks or to the
shared disease system (B cells, Tfh2/CTL4, cytokines, IgG4, fibrosis, IRI),
which is otherwise untouched.

## Archetype determination

Read the actual `dxdt_` lines rather than assuming from the header comment
("TMDD 2-compartment" / "SC, 2-compartment").

**Rituximab (RTX):** `CENT_RTX`/`PERI_RTX` (no depot, IV) with a genuine
receptor-binding pair `CD20_FREE`/`RTX_CD20`, each with its own ODE
(`ksyn_CD20 - kdeg_CD20*CD20_FREE - kon_r + koff_r` /
`kon_r - koff_r - kint_r`) — **Archetype 4 (TMDD)**, no depot variant.

**Dupilumab (DUP):** `SC_DUP` (SC depot) feeding `CENT_DUP` (no peripheral),
with its own receptor-binding pair `IL4RA_FREE`/`DUP_IL4RA` — **Archetype 4
(TMDD)**, depot-but-no-peripheral variant.

**Prednisone (PRED):** `GUT_PRED` (oral depot) feeding `CENT_PRED`, plain
linear elimination, no receptor binding at all — **Archetype 3 minus the
peripheral compartment** (depot + central, linear), the recurring corpus
variant already documented in `eosinophilic-esophagitis/eoe_refactor_notes.md`.

**"Compartment IS concentration" for all three centrals.** Every one of
`CENT_RTX`/`CENT_PRED`/`CENT_DUP` is read *undivided* by the original's own
`$TABLE` (`RTX_conc_ug_mL = CENT_RTX * 148000/1e6*1000`, `PRED_conc_nM =
CENT_PRED`, `DUP_conc_ug_mL = CENT_DUP * 146000/1e6*1000` — none divide by a
volume) and by the TMDD binding terms themselves (`kon_r = kon_RTX * CENT_RTX
* CD20_FREE` uses `CENT_RTX` directly, not `CENT_RTX/V1_RTX`). The dosing
helper `mg_to_nM()` pre-divides every dose by the compound's own volume
before injecting it as `amt`, confirming the state is intended as a
concentration throughout, not a mass amount — the same idiom documented for
tocilizumab in `sepsis/sep_refactor_notes.md` and for mepolizumab/cendakimab
in `eosinophilic-esophagitis/eoe_refactor_notes.md`. So `C_RTX = CENT_RTX`,
`C_PRED = CENT_PRED`, `C_DUP = CENT_DUP`, all undivided — `V1_RTX`/`V1_PRED`/
`V1_DUP` still appear (renamed only) inside each compound's own
clearance-rate-constant terms (`CL/V1`, `Q/V1`, `Q/V2`), exactly as in the
original.

## Renaming applied

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `CENT_RTX` | `CENT_RTX` (unchanged) | `SC_DUP` | `GUT_DUP` |
| `PERI_RTX` | `PERI_RTX` (unchanged) | `CENT_DUP` | `CENT_DUP` (unchanged) |
| `CD20_FREE` | `REC_FREE_RTX` | `IL4RA_FREE` | `REC_FREE_DUP` |
| `RTX_CD20` | `COMPLEX_RTX` | `DUP_IL4RA` | `COMPLEX_DUP` |
| `kon_RTX` | `KON_RTX` | `ka_DUP` | `KA_DUP` |
| `koff_RTX` | `KOFF_RTX` | `kon_DUP` | `KON_DUP` |
| `CD20_ss` | `RTOT_RTX` | `koff_DUP` | `KOFF_DUP` |
| `ksyn_CD20` | `KSYN_RTX` | `IL4RA_ss` | `RTOT_DUP` |
| `kdeg_CD20` | `KDEG_RTX` | `ksyn_IL4RA` | `KSYN_DUP` |
| `kint_RTX` | `KINT_RTX` | `kdeg_IL4RA` | `KDEG_DUP` |
| `Emax_RTX` | `EMAX_RTX` | `kint_DUP` | `KINT_DUP` |
| `ka_PRED` | `KA_PRED` | `Emax_DUP` | `EMAX_DUP` |
| `V_PRED` | `V1_PRED` | | |
| `Emax_PRED` | `EMAX_PRED` | | |

`CL_RTX`, `V1_RTX`, `Q_RTX`, `V2_RTX`, `CL_PRED`, `F_PRED`, `CL_DUP`,
`V1_DUP`, `EC50_PRED` already matched the convention and are unchanged.
Every renamed compartment's `dxdt_<name>` line was renamed identically
(`dxdt_CD20_FREE` → `dxdt_REC_FREE_RTX`, etc. — this needed an explicit
second substitution pass beyond a plain word-boundary rename, since `dxdt_`
is joined to the compartment name by an underscore with no regex word
boundary between them; caught and fixed before compiling). All parameter
*values* are copied verbatim from the original — nothing invented, nothing
dropped. `GUT_DUP`/`CENT_DUP`/`REC_FREE_DUP`/`COMPLEX_DUP` keep the same
declared order/position as `SC_DUP`/`CENT_DUP`/`IL4RA_FREE`/`DUP_IL4RA` did
(7th–10th compartment), so any numeric-index dosing is unaffected; the
outer R script's one **name**-based dosing call for this compartment
(`ev(..., cmt = "SC_DUP", ...)` in Scenario 6) was updated to
`cmt = "GUT_DUP"` — the only necessary edit to the post-DSL R script.

Added (new, not present in the original; naming convention's Hill-interface
slots): `GAMMA_RTX = 1`, `GAMMA_PRED = 1`, `EC50_DUP = 0.503333`,
`GAMMA_DUP = 1` (all derivations below).

## Hill interface

### Prednisone — rename, not a fit

`Eimmu_PRED = Emax_PRED * CENT_PRED / (EC50_PRED + CENT_PRED)` is already
the canonical shape. Renamed directly:
`EFFECT_PRED = EMAX_PRED * pow(C_PRED, GAMMA_PRED) / (pow(EC50_PRED,
GAMMA_PRED) + pow(C_PRED, GAMMA_PRED))`, `GAMMA_PRED = 1` (new; original had
no explicit Hill exponent). `EMAX_PRED = 0.85`, `EC50_PRED = 80.0` are the
original's own values, untouched.

### Rituximab and Dupilumab — TMDD occupancy, rename plus a documented closed form

Both originals compute the disease-facing effect as an **algebraic ratio of
the current ODE state** (`CD20_occ`/`IL4RA_occ`), not as output of a
separate downstream fit — the same situation documented in
`rheumatoid-arthritis/ra_refactor_notes.md` for tocilizumab. Per the guide,
this is a rename: `EFFECT_RTX = Ekill_RTX` (i.e. `EMAX_RTX * CD20_occ`,
`EMAX_RTX = 0.98` unchanged) and `EFFECT_DUP = IL4RA_occ` are computed as
plain aliases of the original's own arithmetic — zero added approximation.

**`EMAX_DUP` is a pre-existing dead parameter, not repurposed.** The
original declares `Emax_DUP = 0.90` ("max IL-4/IL-13 blockade") but never
reads it anywhere in `$ODE` — `IL4RA_occ` is used directly, unscaled, in
every downstream equation (`dxdt_GCB`, `dxdt_TFH2`, `dxdt_IL4`,
`dxdt_TGFB`). Renamed to `EMAX_DUP` and left just as unused as before (same
treatment as `RTOT_TCZ` in the RA precedent) — `EFFECT_DUP` does **not**
multiply by it, matching the original's own behaviour exactly.

**Analytic derivation: occupancy vs. concentration is exact Hill(gamma=1)
for both compounds, even though free-receptor and complex have *different*
decay constants (`KDEG_*` vs. `KINT_*`), unlike the RA/TCZ precedent's
single shared `KDEG_TCZ`.** At quasi-steady-state for fixed `C`:

```
REC_FREE = KSYN / (KDEG + a*C),   a = KON*KINT/(KOFF+KINT)
COMPLEX  = [KON*C/(KOFF+KINT)] * REC_FREE
Occ(C)   = COMPLEX/(REC_FREE+COMPLEX) = C / (C + Kd),   Kd = (KOFF+KINT)/KON
```

`KSYN`/`KDEG` (free-receptor turnover) drop out of the *fractional*
occupancy entirely — they only set `RTOT(C)`'s absolute size, not the
free:complex ratio. This gives:

- RTX: `Kd_RTX = (KOFF_RTX + KINT_RTX)/KON_RTX = (0.002+0.3)/0.55 =
  0.549091` nM
- DUP: `Kd_DUP = (KOFF_DUP + KINT_DUP)/KON_DUP = (0.001+0.15)/0.30 =
  0.503333` nM

**Numerical confirmation** (plain Python, `scipy.integrate.solve_ivp`,
`LSODA`, `rtol=1e-10`): integrated the isolated 2-state receptor-binding
system to steady state at 10 concentrations per compound, log-spaced from
below to far above each `Kd` (RTX: 0.01–2,252,252 nM, spanning the actual
Scenario-3/4/5 dose range; DUP: 0.001–456,621 nM, spanning Scenario 6's
range), and compared to the closed form:

| Compound | Max abs diff (occupancy) |
|---|---|
| RTX | 3.4e-15 |
| DUP | 3.0e-13 |

Effectively machine precision — R² = 1.0000000 in all but name. Binding
equilibrates in well under a day at every dose actually used (`kon_r` ~
1e5–1e6 nM/day at Scenario-3 doses vs. disease-compartment turnover rates of
0.005–0.4/day), so the quasi-steady-state approximation is well justified
relative to the dosing/output timescale, matching the same argument made in
the RA/TCZ precedent.

**`GAMMA_RTX = 1`, `EC50_DUP = 0.503333`, `GAMMA_DUP = 1` are exposed as new
`$PARAM` values documenting this closed form for future covariate
substitution — `EFFECT_RTX`/`EFFECT_DUP` themselves still compute from the
exact ODE state ratio, not the closed form, so this adds zero approximation
on top of an already-exact relationship** (same reasoning, same result, as
the RA precedent's `EMAX_TCZ`/`EC50_TCZ`/`GAMMA_TCZ`).

**Why there is no new `EC50_RTX`:** the original *already* declares
`EC50_RTX = 50.0` ("nM RTX (central) for 50% B cell kill") — but it is
**never read anywhere in the original's `$ODE`** (`Ekill_RTX = Emax_RTX *
CD20_occ` has no EC50 term at all; it is a dead parameter, same class as
`RTOT_TCZ` in the RA precedent). Its value (50) does not match the derived
`Kd_RTX` (0.549091) — these are two different numbers with two different
provenances (one an unused original guess, one a fresh analytic/numeric
derivation). Per "parameter values always come from the original — never
invent," the pre-existing `EC50_RTX = 50.0` was left completely untouched
rather than silently overwritten with the derived value (which would be an
undisclosed value change to a named parameter, even though it happens to
have zero effect on any simulated output either way, since it's dead in
both directions). The derived `Kd_RTX = 0.549091` is documented here and in
the header comment, not stored under any `$PARAM` name, to avoid two
different numbers sharing one name's "meaning."

## `$ODE`/`$TABLE` scope collision (mrgsolve 2.0.1), and how it was resolved

First compile attempt declared `C_RTX`/`EFFECT_RTX`/`C_PRED`/`EFFECT_PRED`/
`C_DUP`/`EFFECT_DUP` via `double X = ...;` in **both** `$ODE` (for the
disease equations to read) and `$TABLE` (recomputed fresh, since `$ODE` and
`$TABLE` locals are not visible to each other for read purposes — confirmed
by the fact the *original* file already recomputes `CD20_occ_pct`/
`IL4RA_occ_pct` fresh in `$TABLE` rather than reusing the `$ODE`-local
`CD20_occ`/`IL4RA_occ`). This failed:

```
68:10: error: redefinition of 'double {anonymous}::C_RTX'
37:10: note: 'double {anonymous}::C_RTX' previously declared here
```

mrgsolve 2.0.1 shares one C++ scope for any name that is also listed in
`$CAPTURE`, regardless of which block(s) declare it — so the *same* `double
X = ...;` pattern appearing once in `$ODE` and again in `$TABLE` is a
genuine redeclaration, not two independent locals. **Fixed the same way as
`eosinophilic-esophagitis/eoe_mrgsolve_model_refactored.R`'s
`C_MEPO_OUT`/`EFFECT_MEPO_OUT`:** the `$ODE`-local names (`C_RTX`,
`EFFECT_RTX`, `C_PRED`, `EFFECT_PRED`, `C_DUP`, `EFFECT_DUP`) are used
by the disease equations exactly as before (unchanged, zero numeric
effect); the `$TABLE`-recomputed duplicates carry an `_OUT` suffix
(`C_RTX_OUT`, `EFFECT_RTX_OUT`, `C_PRED_OUT`, `EFFECT_PRED_OUT`,
`C_DUP_OUT`, `EFFECT_DUP_OUT`) and are what is listed in `$CAPTURE` and
therefore what `/model_manifest`'s `outputPaths` and `/run_simulation`'s
`outputs` actually expose. Confirmed via `/model_manifest`: all six
`_OUT` names appear in `outputPaths`, and `EMAX_RTX`/`EC50_RTX`/
`GAMMA_RTX`/`EMAX_PRED`/`EC50_PRED`/`GAMMA_PRED`/`EMAX_DUP`/`EC50_DUP`/
`GAMMA_DUP` all appear in `parameters`.

## Upstream defects found and logged (not fixed upstream)

Two independent, pre-existing defects, both confirmed on the untouched
original via the qspserver `mrgsolve_api` container, both unrelated to any
of the three in-scope compounds' own PK/effect equations. Logged as
`translations/UPSTREAM_ISSUES.md` #52 and #53.

**#52 — build defect, fixed forward into the delivered file per the guide's
settled policy.** `double TGFB_eff = TGFB^2 / (TGFB_EC50_FIB^2 + TGFB^2);`
uses R's `^` power operator, invalid in the C++ mrgsolve compiles this DSL
block to (`invalid operands ... to binary 'operator^'`) — the model does
not compile at all, refactored or not, until this is fixed. Rewrote as
`pow(TGFB,2)` / `pow(TGFB_EC50_FIB,2)` **directly in
`igg4rd_mrgsolve_model_refactored.R`** (algebraically identical for all
real inputs; confirmed by the exact 0.0 max-abs-diff verification result
below, including on `Fibrosis_idx`/`ECM`, which depend on this exact
line). The tracked original is untouched and still does not compile as-is.

**#53 — a scenario-definition defect discovered while verifying, not a
build defect and not fixed anywhere (workaround used for verification
only).** Scenario 2 ("Prednisone 40mg/d taper") doses all 168 days with
`rate = -2`, which requires a modeled infusion-duration parameter
(`D_GUT_PRED`) that the original never defines anywhere. Reproducing this
exact dosing against the untouched (only build-compat-patched) original via
`/run_simulation` throws `[mrgsolve] modeled infusion duration D_CMT or Dn
must be positive when dosing record RATE is set to -2` — meaning the
original's own R script, run top to bottom as written, would itself crash
on this scenario. For verification purposes only, the same 168
time/amt/cmt triples were resubmitted with `rate = 0` (plain bolus) against
**both** models — this substitutes the administration-route representation
only; it says nothing about, and does not touch, PRED's own PK parameters
or effect equation (unchanged, pure rename, see above).

## Verification

Per the guide's mandatory protocol: ran all **six** of the original file's
own named dosing scenarios (`S1`–`S6`, exactly as constructed in the
original's post-DSL R code, dose amounts computed with the same
`mg_to_nM()`/taper formulas) through both the original (patched only for
the two defects above, in-memory scratch copies, never written to either
tracked file) and the refactored model, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`), full `end=730,
delta=1` daily grid as in the original (no shortening needed — none of the
six scenarios approached mrgsolve's default 20000-step budget).

| Scenario | Dosing | Time points | Result |
|---|---|---|---|
| S1: Untreated natural history | none (zero-dose event) | 732 | exact match, 0.0 |
| S2: Prednisone 40mg/d taper | 168 daily `GUT_PRED` doses (rate=0 substitution, see #53) | 899 | exact match, 0.0 |
| S3: RTX 1g IV D1+D15 (Khosroshahi) | 2× `CENT_RTX` bolus | 733 | exact match, 0.0 |
| S4: RTX 375mg/m² ×4 weekly | 4× `CENT_RTX` bolus | 735 | exact match, 0.0 |
| S5: RTX induction + 500mg q6m maintenance | 4× `CENT_RTX` bolus | 735 | exact match, 0.0 |
| S6: Dupilumab 300mg SC q2w | 51× `GUT_DUP` (was `SC_DUP`) bolus | 782 | exact match, 0.0 |

For every scenario, all **37 shared outputs** were compared point-by-point
across the full time grid: the 14 unchanged `$CAPTURE` names
(`RTX_conc_ug_mL`, `PRED_conc_nM`, `DUP_conc_ug_mL`, `Bcell_pct`,
`IgG4_mgdL`, `Fibrosis_idx`, `Activity_IRI`, `CD20_occ_pct`,
`IL4RA_occ_pct`, `TGFB_rel`, `PC_rel`, `TFH2_rel`, `CR_flag`, `PR_flag`)
plus the 23 raw compartment states, mapped through the renames above
(`CD20_FREE→REC_FREE_RTX`, `RTX_CD20→COMPLEX_RTX`, `SC_DUP→GUT_DUP`,
`IL4RA_FREE→REC_FREE_DUP`, `DUP_IL4RA→COMPLEX_DUP`, all others unchanged
names). **Result: exact match, max absolute difference 0.0 on every output,
every scenario** — the expected outcome for a pure structural
reorganization with no Hill-fitting (Archetypes 3/4, gamma=1 throughout,
`EFFECT_*` computed as literal aliases of the originals' own arithmetic).

Sanity-checked the dynamics are non-trivial, not degenerately zero on both
sides (Scenario 3, RTX 1g D1+D15): `IgG4_mgdL` falls from 450 to 127.99
mg/dL by day 84 (~71.6% reduction, directionally consistent with the
Khosroshahi calibration note in the model's own header, ~75–80% at 3
months — the original's own model, not something this refactor changed);
`Bcell_pct` nadirs at 0.3% of baseline (>99% depletion, consistent with the
Carruthers calibration note of <5/μL / ~95% depletion); `CD20_occ_pct`
saturates to ~99.9995% within 1 day of the first dose, confirming the
fast-binding assumption behind the quasi-steady-state Hill derivation
above. The new `_OUT` outputs behave sensibly in the refactored model:
`EFFECT_RTX_OUT` reads `0.98` (`= EMAX_RTX`) once CD20 is saturated, and
`C_RTX_OUT` matches `CENT_RTX` exactly at every timepoint, as designed.

## qspserver compatibility

Confirmed via `/model_manifest` on the finished `_refactored.R`'s extracted
DSL: `outputPaths` includes all 23 renamed/unchanged compartments, the 14
unchanged `$CAPTURE` names, and the 6 new `_OUT` names; `parameters`
includes every `$PARAM` entry, including the 4 new Hill-interface
parameters (`GAMMA_RTX`, `GAMMA_PRED`, `EC50_DUP`, `GAMMA_DUP`) and the
renamed originals (`EMAX_RTX`, `EC50_RTX`, `EMAX_PRED`, `EC50_PRED`,
`EMAX_DUP`). The extracted `.cpp` used for every verification call above
(not a delivered artifact — deleted after use) was confirmed byte-identical
to the quoted `code <- '...'` string re-extracted from the finished
`igg4rd_mrgsolve_model_refactored.R`.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`igg4-related-disease` rows for `DUP`, `PRED`, and `RTX`.
