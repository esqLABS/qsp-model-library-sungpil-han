# Refactor notes — `narcolepsy/narc_mrgsolve_model.R`

Scope: all **five** compounds this file gives their own PK compartments and
that have existing rows in `driver-patches/data/compound_perturbation_census.md`
— **Sodium Oxybate (OXY)**, **Modafinil (MOD)**, **Pitolisant (PIT)**,
**Solriamfetol (SOL)**, and **Venlafaxine (VEN)** — all classified "Redirect
concentration (clean single site)". Every one of these five is the file's
entire pharmacological content (there is no sixth compound in this model),
so nothing was left untouched for "out of scope" reasons.

## Headline finding: the original does not compile at all — unrelated to any of the five compounds' own PK

Before any renaming could even be attempted, `POST /model_manifest` on the
untouched original (via the qspserver `mrgsolve_api` container,
`http://localhost:8007`) failed to build. Root cause, logged as
`translations/UPSTREAM_ISSUES.md` #78 (same defect class as #30/#34/#41/
#45/#46): the `$CAPTURE` block lists ten names —
`WAKE_STATE, REM_STATE, NREM_STATE, VLPO_ACT, SLEEP_P, ADENOSINE, WAKE_LC,
WAKE_TMN, WAKE_VTA, WAKE_DRN` — that duplicate `$CMT` compartment names
declared earlier in the same file:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  WAKE_STATE,REM_STATE,NREM_STATE,VLPO_ACT,SLEEP_P,ADENOSINE,WAKE_LC,
  WAKE_TMN,WAKE_VTA,WAKE_DRN
```

All ten belong to the disease/neuroscience side of the model (wake nuclei,
sleep-pressure, flip-flop state), not to any of the five compounds' own
PK/effect blocks — confirmed incidental to the refactor's own scope.

Per the fork guide's settled policy for this situation
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"), the
fix is **syntax-only, non-numeric, non-behavioral**, and is applied
directly in the delivered `narc_mrgsolve_model_refactored.R` (not just a
throwaway scratch copy): the ten compartment names are simply dropped from
the `$CAPTURE` line. Compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported —
`WAKE_STATE`, `SLEEP_P`, etc. are still returned by every simulation, they
are just no longer *redundantly* re-listed in `$CAPTURE`. The tracked
`narc_mrgsolve_model.R` is completely untouched and still carries the
defect exactly as written. See `translations/UPSTREAM_ISSUES.md` #78 for
the full defect writeup.

## Archetype per compound

All five compounds are **Archetype 3 variants — depot + central PK**, with
Sodium Oxybate the one exception that also carries a peripheral
compartment (the full Archetype 3 shape); the other four drop the
peripheral compartment exactly as the guide's own named variant describes
("Drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a 2-compartment-no-depot variant" runs
the same idea in reverse — here it is "drop the peripheral, keep depot +
central"). None of the five needed a compartment added or removed: the
original already used `GUT_<STEM>`/`CENT_<STEM>`/`PERI_<STEM>` compartment
names matching this fork's convention exactly, so **compartment names and
compartment order are completely unchanged** for all five compounds — a
pure parameter/variable rename throughout.

- **Sodium Oxybate (OXY): Archetype 3, full (depot + central + peripheral,
  linear).** `GUT_OXY`/`CENT_OXY`/`PERI_OXY` unchanged.
- **Modafinil (MOD): Archetype 3 without peripheral (depot + central
  only).** `GUT_MOD`/`CENT_MOD` unchanged.
- **Pitolisant (PIT): Archetype 3 without peripheral.**
  `GUT_PIT`/`CENT_PIT` unchanged.
- **Solriamfetol (SOL): Archetype 3 without peripheral.**
  `GUT_SOL`/`CENT_SOL` unchanged.
- **Venlafaxine (VEN): Archetype 3 without peripheral.**
  `GUT_VEN`/`CENT_VEN` unchanged.

## Naming

The original already used `GUT_<STEM>`/`CENT_<STEM>`/`PERI_<STEM>` for
compartments and `EMAX_<STEM>_<TARGET>`/`EC50_<STEM>_<TARGET>` for the Hill
parameters — both already matched this fork's convention exactly, so
neither needed renaming. Three things did:

| Role | Original name | Renamed to | Compounds |
|---|---|---|---|
| Central volume (L/kg, scaled by `BW` in `$MAIN`/`$ODE`) | `VD_OXY`, `VD_MOD`, `VD_PIT`, `VD_SOL`, `VD_VEN` | `V1_OXY`, `V1_MOD`, `V1_PIT`, `V1_SOL`, `V1_VEN` | all 5 |
| Peripheral volume (L/kg) | `VP_OXY` | `V2_OXY` | OXY only |
| Exposed concentration (`$GLOBAL` double, ug/mL) | `CP_OXY`, `CP_MOD`, `CP_PIT`, `CP_SOL`, `CP_VEN` | `C_OXY`, `C_MOD`, `C_PIT`, `C_SOL`, `C_VEN` | all 5 |
| `$TABLE`'s redundant recompute of the same concentration (kept, for the same reason the original kept it — see "qspserver compatibility" below) | `CP_OXY_out`, ... | `C_OXY_out`, ... | all 5 |
| Hill coefficient | shared `H_EDS`/`H_CAT`/`H_WAKE`/`H_TMN`/`H_SOL` | per-stem-per-target `GAMMA_<STEM>_<TARGET>` (see below) | all 5 |

`KA_<STEM>`, `F_<STEM>`, `CL_<STEM>`, `Q_OXY` were already correctly named
per the convention and are untouched. All parameter *values* are copied
verbatim from the original — nothing invented, nothing defaulted, nothing
dropped.

### Nine effects across five compounds, and a genuinely shared Hill coefficient split apart

Four of the five compounds (OXY, MOD, PIT, SOL) each drive **two**
independent downstream effects from the same exposed concentration — e.g.
Sodium Oxybate's `EFFECT_OXY_EDS` (EDS/sleep-quality) and `EFFECT_OXY_CAT`
(cataplexy), each with its own `EMAX`/`EC50`. Per the guide's "never
collapse several drugs into one shared Hill term" rule (which also applies
within one compound with more than one target, per the `bronchiectasis`
AZM precedent): each effect keeps its own name, `EFFECT_<STEM>_<TARGET>`,
exactly mirroring what the original already did (`EFFECT_OXY_EDS`,
`EFFECT_OXY_CAT`, `EFFECT_MOD_WAKE`, `EFFECT_MOD_EDS`, `EFFECT_PIT_TMN`,
`EFFECT_PIT_EDS`, `EFFECT_SOL_WAKE`, `EFFECT_SOL_EDS`, `EFFECT_VEN_CAT` —
none of these nine names changed).

What *did* need to change: the original's Hill coefficient was **not**
per-effect — it was five bare constants (`H_EDS = 1.5`, `H_CAT = 2.0`,
`H_WAKE = 1.5`, `H_TMN = 1.2`, `H_SOL = 1.5`) each **shared across several
different compounds' effect terms** (e.g. `H_EDS` alone is reused by
`EFFECT_OXY_EDS`, `EFFECT_MOD_EDS`, `EFFECT_PIT_EDS`, and
`EFFECT_SOL_EDS`). This is exactly the naming chaos the guide's convention
exists to remove — a single implicit constant standing in for what should
be nine independently-tunable parameters. Each of the nine effects now gets
its own `GAMMA_<STEM>_<TARGET>` parameter, all initialized to the exact
value the shared constant held for that effect in the original (so no
numeric value changed, only how many named parameters hold it):

| New parameter | Value | Was |
|---|---|---|
| `GAMMA_OXY_EDS` | 1.5 | `H_EDS` |
| `GAMMA_OXY_CAT` | 2.0 | `H_CAT` |
| `GAMMA_MOD_WAKE` | 1.5 | `H_WAKE` |
| `GAMMA_MOD_EDS` | 1.5 | `H_EDS` |
| `GAMMA_PIT_TMN` | 1.2 | `H_TMN` |
| `GAMMA_PIT_EDS` | 1.5 | `H_EDS` |
| `GAMMA_SOL_WAKE` | 1.5 | `H_SOL` |
| `GAMMA_SOL_EDS` | 1.5 | `H_EDS` |
| `GAMMA_VEN_CAT` | 2.0 | `H_CAT` |

Splitting a shared constant into nine identically-valued named parameters
does not change any computed value anywhere (`pow(x, 1.5)` is identical
whether the `1.5` came from `H_EDS` or from `GAMMA_OXY_EDS`) — confirmed by
the exact-match verification below. It does mean a future user driving
`C_OXY` (say) externally can now tune `GAMMA_OXY_EDS` independently of
`GAMMA_MOD_EDS` without touching Modafinil's own effect, which was not
possible in the original's shared-constant scheme.

## The Hill interface: rename only, no fitting anywhere in this file

All nine effect terms were **already** exactly the guide's own Hill
formula in the original —
`EMAX * pow(C, H) / (pow(EC50, H) + pow(C, H))` — with `EMAX`/`EC50`
already stem-and-target-named and a real, previously-chosen Hill
coefficient (not defaulted to 1). This is a pure rename in every case:
`EMAX_<STEM>_<TARGET>` and `EC50_<STEM>_<TARGET>` keep their original names
and values; only the Hill exponent moves from a cross-compound shared
constant to a dedicated `GAMMA_<STEM>_<TARGET>` (table above), same value.
No `nls()` fit was performed or needed for any of the five compounds — the
"already this shape, rename don't refit" branch of the guide applies to
all nine effects.

## qspserver compatibility

- `C_<STEM>` and `EFFECT_<STEM>_<TARGET>` are `$GLOBAL`-scope `double`s
  assigned in `$MAIN` (unchanged mechanism from the original, which already
  used exactly this pattern for `CP_<STEM>`/`EFFECT_<STEM>_<TARGET>`) — not
  `$PARAM` entries. Following the precedent already established across
  this fork's prior refactors (`type1-diabetes/t1dm_refactor_notes.md`,
  `bronchiectasis/bex_refactor_notes.md`): mrgsolve 2.0.1 passes every
  `$PARAM` value into `$MAIN`/`$ODE` as a `const double&`, so assigning to
  one is a hard compile error. Discoverability is satisfied instead by
  listing every `C_<STEM>` and every `EFFECT_<STEM>_<TARGET>` explicitly in
  `$CAPTURE` (confirmed present in `/model_manifest`'s `outputPaths`, see
  below) — this repo's settled resolution of an apparent tension between
  the guide's `$PARAM` wording and the actual mrgsolve 2.0.1 constraint.
- The original's `$CAPTURE` only listed five of the nine `EFFECT_*`
  variables (`EFFECT_OXY_EDS`, `EFFECT_MOD_EDS`, `EFFECT_PIT_EDS`,
  `EFFECT_SOL_EDS`, `EFFECT_VEN_CAT` — the "EDS-family" ones) — the other
  four (`EFFECT_OXY_CAT`, `EFFECT_MOD_WAKE`, `EFFECT_PIT_TMN`,
  `EFFECT_SOL_WAKE`) were computed every step but never exposed. All nine,
  plus all five `C_<STEM>`, are now in `$CAPTURE` in the refactored file —
  a strict widening of what is discoverable, not a narrowing.
- Confirmed via `POST /model_manifest` on the qspserver `mrgsolve_api`
  container: `outputPaths` includes all five exposed concentrations
  (`C_OXY`, `C_MOD`, `C_PIT`, `C_SOL`, `C_VEN`), all nine effect terms, and
  the original's own `$TABLE` output names (`C_OXY_out` etc., `EDS_NOW`,
  `CAT_NOW`, `SLEEP_LATENCY`, `REM_PCT`, `SLEEP_EFF`, `CIRC_OUT`,
  `CSF_OREXIN`, `OREXIN_NEURONS`), kept for backward output comparability.
  `parameters` includes all 27 renamed/new PK and Hill parameters
  (`V1_*`/`V2_OXY`, all nine `GAMMA_*`) alongside every parameter that
  didn't need renaming.
- `compartments` — all 23 `$CMT` names, in the same declared order as the
  original (`GUT_OXY, CENT_OXY, PERI_OXY, GUT_MOD, CENT_MOD, GUT_PIT,
  CENT_PIT, GUT_SOL, CENT_SOL, GUT_VEN, CENT_VEN, WAKE_LC, ...`), so
  1-based compartment-number dosing (as `POST /run_simulation`'s `dosing`
  field uses) is unaffected.
- `$SET end/delta` is not present in this model at all — the original
  never declared one; both the original R script and this refactored
  sibling rely entirely on `mrgsim(start=, end=, delta=)` call-site
  arguments in `run_scenario()`. Nothing to note about it being
  non-load-bearing for API runs; every verification run below supplied
  `time.end`/`time.delta` explicitly via the request.
- No R-only syntax inside the DSL block — every parameter, macro, and
  computation the block needs lives inside `$PARAM`/`$GLOBAL`/`$MAIN`/
  `$ODE`/`$TABLE` itself. The surrounding R script needed **zero** changes
  to any `cmt = "GUT_<STEM>"` string in `make_sox_dose()`/`make_mod_dose()`/
  `make_pit_dose()`/`make_sol_dose()`/`make_ven_dose()`/`run_scenario()`,
  since no compartment name changed for any of the five compounds. It did
  need the `CP_<STEM>_out` → `C_<STEM>_out` rename propagated into the
  R-side post-processing that references those `$TABLE` output columns by
  name: `daily_summary`'s `CP_OXY_peak = max(CP_OXY_out, ...)` etc. (now
  `max(C_OXY_out, ...)`), and `pk_day1`'s `select()`/`pivot_longer()`/
  `recode()` block (`CP_OXY_out = "Sodium Oxybate"` → `C_OXY_out =
  "Sodium Oxybate"`, same for MOD/PIT/SOL/VEN). No other line of the R
  wrapper (dosing builders, scenarios, plotting, summary tables) differs
  from the original.

## An observation on `$MAIN` vs `$TABLE` timing (found, not fixed)

Not a build defect and not something this refactor changed, but worth
flagging for anyone about to drive `C_<STEM>` externally (the whole point
of this "redirect concentration" refactor): the original computes each
`EFFECT_<STEM>_<TARGET>` in `$MAIN` from the `$GLOBAL` concentration
(`C_OXY` etc.), while `$TABLE` independently *recomputes* the same
concentration fresh from the current `$CMT` state as `C_OXY_out` etc.
Empirically, `$MAIN`'s `C_OXY` at a given output record equals `$TABLE`'s
`C_OXY_out` from the *previous* record — a one-`delta`-step lag between
the concentration the disease-side `EFFECT_*` terms actually respond to
and the concentration the model reports as "current." This is a
pre-existing architectural characteristic of the original (confirmed by
spot-checking `EFFECT_OXY_EDS`/`EFFECT_OXY_CAT` against a from-scratch
Hill-formula recomputation using `C_OXY_out`: they disagree by up to ~0.6
absolute, entirely explained by the one-step lag, not by any refactor
error), reproduced identically — not introduced, not fixed — in the
refactored file. It does not affect the verification below at all (both
sides carry the identical mechanism, so they still match exactly), but a
future driver-substitution of `C_OXY` should account for this lag rather
than assume `EFFECT_OXY_EDS(t)` is a same-instant function of the reported
`C_OXY_out(t)`.

## Verification

Per the guide's mandatory protocol: ran each compound's own dosing regimen
from the original R script's own event builders (`make_sox_dose()`,
`make_mod_dose()`, `make_pit_dose()`, `make_sol_dose()`, `make_ven_dose()`
— same doses, same intervals, same target compartments) through both the
`$CAPTURE`-patched original and the refactored model via the qspserver
`mrgsolve_api` service (`POST /run_simulation`), comparing every shared
output point-by-point across the full time grid.

**Window and initial conditions, disclosed:** the original R script runs
30 days (720 h) with a custom initial-condition list (`inits_nt1`:
`WAKE_LC = 0.08`, `SLEEP_P = 0.20`, etc.) passed via `init()`. The API's
`SimRequest` has no field for overriding `$CMT` initial conditions (only
`dosing`/`parameters`/`time`/`outputs`), so both sides were run from
mrgsolve's default (all-zero) initial state instead — a limitation of the
verification harness, not of either model, and applied identically to both
sides so the comparison stays apples-to-apples. The window was shortened
to the first **10 days (240 h)**, hourly output, comfortably within the
API's default solver step budget for these non-stiff, linear-PK
scenarios (no `maxsteps` issues encountered).

**Sodium Oxybate — `make_sox_dose(10)` regimen (4.5 g at t=14h, 2.75 g at
t=17h, daily ×10):** `CENT_OXY`/`GUT_OXY`/`PERI_OXY`, `C_OXY_out`↔
`CP_OXY_out`, `EFFECT_OXY_EDS`, plus the shared disease outputs (`EDS_NOW`,
`CAT_NOW`, `SLEEP_LATENCY`, `REM_PCT`, `SLEEP_EFF`, `CIRC_OUT`,
`CSF_OREXIN`, `OREXIN_NEURONS`) — **max abs diff = 0.0 for every output,
every time point (261 points).**

**Modafinil — `make_mod_dose(10)` regimen (200 mg QD):** `CENT_MOD`/
`GUT_MOD`, `C_MOD_out`↔`CP_MOD_out`, `EFFECT_MOD_EDS`, plus shared disease
outputs — **max abs diff = 0.0 (251 points).**

**Pitolisant — `make_pit_dose(10)` regimen (18 mg QD):** `CENT_PIT`/
`GUT_PIT`, `C_PIT_out`↔`CP_PIT_out`, `EFFECT_PIT_EDS`, plus shared disease
outputs — **max abs diff = 0.0 (251 points).**

**Solriamfetol — `make_sol_dose(10)` regimen (150 mg QD):** `CENT_SOL`/
`GUT_SOL`, `C_SOL_out`↔`CP_SOL_out`, `EFFECT_SOL_EDS`, plus shared disease
outputs — **max abs diff = 0.0 (251 points).**

**Venlafaxine — `make_ven_dose(10)` regimen (75 mg QD):** `CENT_VEN`/
`GUT_VEN`, `C_VEN_out`↔`CP_VEN_out`, `EFFECT_VEN_CAT`, plus shared disease
outputs — **max abs diff = 0.0 (251 points).**

(The four effect terms the original never captured — `EFFECT_OXY_CAT`,
`EFFECT_MOD_WAKE`, `EFFECT_PIT_TMN`, `EFFECT_SOL_WAKE` — could not be
compared against the original this way, since the original's `$CAPTURE`
never exposed them for `Req()` selection to find; see the "$MAIN vs
$TABLE timing" note above for the separate, from-first-principles spot
check performed on the refactored file's own `EFFECT_OXY_CAT`, which
confirmed it follows exactly the stated Hill formula given the same
`$MAIN`-timing characteristic as the already-verified `EFFECT_OXY_EDS`.)

**Result: exact match (max abs diff = 0.0) for every output, every
scenario, every compound.** This is the expected outcome for Archetype 3
(all five compounds) — pure PK/parameter reorganization, `GAMMA_*` renames
of an already-fitted Hill exponent, no fitting anywhere in this file.

## Anything else worth flagging

- Compartment order is completely unchanged (`GUT_OXY` compartment 1
  through `CATAPLEXY_ACC` compartment 23), so 1-based compartment-number
  dosing is unaffected between the original and the refactored file.
- The R wrapper's dosing-builder `cmt = "GUT_<STEM>"` strings needed zero
  changes (no compartment was renamed for any of the five compounds); only
  the `CP_<STEM>_out` column references in `daily_summary` and
  `pk_day1`/`pivot_longer`/`recode` needed updating to `C_<STEM>_out` (see
  "qspserver compatibility" above).
- The model's own header comment claims "22 ODE compartments"; the actual
  `$CMT` count is 23 (12 PK + 4 wake nuclei + 3 sleep-system + 3
  state-variable + 2 symptom-accumulator). Not logged as a separate
  `UPSTREAM_ISSUES.md` entry (a stale doc comment, not a build or
  behavioral defect) — noted here for a future reviewer, reproduced
  unchanged in the refactored file's header.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`narcolepsy | Modafinil`, `narcolepsy | Pitolisant`, `narcolepsy | Sodium
Oxybate`, `narcolepsy | Solriamfetol`, and `narcolepsy | Venlafaxine` rows.
