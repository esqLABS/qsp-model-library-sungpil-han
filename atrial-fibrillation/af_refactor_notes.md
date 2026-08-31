# Refactor notes — `atrial-fibrillation/af_mrgsolve_model.R`

Compounds refactored: **Amiodarone (AMIO)**, **Apixaban (APIX)**, and
**Metoprolol (METRO)** — the three rows in
`driver-patches/data/compound_perturbation_census.md` classified
`atrial-fibrillation | Amiodarone`, `atrial-fibrillation | Apixaban`, and
`atrial-fibrillation | Metoprolol (METRO)`. These are the file's entire set
of exogenous compounds — no other compound exists in this model. Every
disease-side compartment and equation (AF burden, ERP/electrical
remodeling, LA size, fibrosis, the RAAS/ROS/SMAD2/3/IL-6 axis, NE, heart
rate, BNP, the FXa/thrombin coagulation cascade, CHA2DS2-VASc stroke risk,
QTc, IKr, and intracellular Ca²⁺) is untouched — confirmed by a
line-level diff of the delivered file against the original outside the six
edited regions (see "Delivered R file" below).

## Metabolite check (per the task's request to confirm from code, not guess)

**No desethylamiodarone (DEA), or any other active metabolite, is modelled
anywhere in this file.** Amiodarone's own PK is a plain 3-compartment
system: `GI_AMIO` (gut depot) → `C1_AMIO` (central plasma) ⇌ `C2_AMIO`
(fat/tissue depot, i.e. a *distribution* peripheral compartment, not a
biotransformation product). There is no `KMET`/`FM`-style formation-rate
parameter, no separate DEA compartment, and no comment anywhere in the
file suggesting a metabolite was intended and dropped. This differs from
the two possible precedents named in the task brief: `chagas-disease`
classifies its own amiodarone row as "Amiodarone and desethylamiodarone"
(that file *does* carry a DEA compartment, per its own census
classification — a genuinely different model, not yet refactored as of
this writing, so not literally available as a worked file to compare
against) and `aortic-dissection` (whose own census rows for this same
class of drug do not mention DEA at all — also not yet refactored). Not
finding a metabolite here is a fact about *this* file's amiodarone model,
independent of either.

## Archetype per compound

- **Amiodarone (AMIO)**: **archetype 3** (depot + central + peripheral,
  linear) with one disclosed deviation from the archetype's textbook
  CL/Q/V form (see below). `GUT_AMIO`/`CENT_AMIO`/`PERI_AMIO` (was
  `GI_AMIO`/`C1_AMIO`/`C2_AMIO`), `KA_AMIO` (was `ka_AMIO`), `F_AMIO`,
  `CL_AMIO`, `V1_AMIO`, `V2_AMIO` unchanged in value. Four separate,
  simultaneous PD effects (ERP prolongation, HR slowing, IKr block, ICaL
  block) — see Hill interface below.
- **Apixaban (APIX)**: **archetype 3 minus peripheral** (depot + central,
  linear, oral). `GUT_APIX`/`CENT_APIX` (was `GI_APIX`/`C1_APIX`),
  `KA_APIX` (was `ka_APIX`), `F_APIX`, `CL_APIX`, `V1_APIX` unchanged.
  One PD effect (FXa inhibition).
- **Metoprolol (METRO)**: **archetype 3 minus peripheral** (depot +
  central, linear, oral). `GUT_METRO`/`CENT_METRO` (was
  `GI_METRO`/`C1_METRO`), `KA_METRO` (was `ka_METRO`), `F_METRO`,
  `CL_METRO`, `V1_METRO` unchanged. One PD effect (HR reduction), which is
  then combined with amiodarone's own HR effect at the single point the
  disease equations use it (`HR_red_combined`, untouched formula, per the
  guide's "combine only at point of use" rule).

No compound needed TMDD (archetype 4) — every drug-target interaction in
this file is a plain saturating Hill ratio on free plasma/central
concentration, not an integrated receptor-binding ODE system.

### Disclosed deviation — amiodarone's `K12_AMIO`/`K21_AMIO` kept as literal, asymmetric micro-constants, not converted to a single `Q_AMIO`

The guide's archetype 3 template calls for the `k10`/`k12`/`k21`
micro-constant form to always be rewritten to `CL`/`Q`/`V`, since in a
mass-balanced two-compartment model `Q = k12·V1 = k21·V2` (one physical
inter-compartmental clearance, reachable equivalently from either side).
Checking that equality against this file's own calibrated values:

```
k12_AMIO * V1_AMIO = 0.020  * 40   = 0.8   /h·L  (i.e. an implied Q of 0.8 L/h)
k21_AMIO * V2_AMIO = 0.0002 * 4200 = 0.84  /h·L  (i.e. an implied Q of 0.84 L/h)
```

These two implied values differ by 5% — the original's `k12_AMIO`/
`k21_AMIO` are **not** mass-balanced, i.e. not actually derived from a
single physical `Q`. They are two independently-calibrated first-order
rate constants (the file's own comment: "`k12/k21` calibrated so
steady-state Cp ≈ 1–2.5 µg/mL" — a calibration target, not a mechanistic
derivation). Forcing them into a single `Q_AMIO/V1_AMIO`/`Q_AMIO/V2_AMIO`
pair, as the archetype's worked example shows, would silently pick one of
the two values over the other and change the calibrated depot-loading
trajectory — exactly the kind of "structural fix that changes the
simulated behavior" the guide says not to do for a pure-rename archetype.
Kept instead as **renamed-but-otherwise-literal** `K12_AMIO`/`K21_AMIO`
(central→peripheral and peripheral→central transfer rate constants,
values unchanged), applied to the amount-scale ODEs exactly as the
original did. `CL_AMIO`/`V1_AMIO` (elimination) were already in the
target `CL`/`V1` form and needed no change.

## PK/PD interface

- **Amiodarone**: `C_AMIO = CENT_AMIO / V1_AMIO` (µg/mL) is the single
  exposed concentration all four of amiodarone's own PD effects read.
  `PERI_AMIO / V2_AMIO` (the fat-depot concentration) remains an
  informational-only `$TABLE` output (`Cp_AMIO_depot`, unchanged name) —
  it was never read by any PD term in the original either, so it is not
  part of the pluggable interface, per the guide's "two [concentration
  sites] only when a genuinely different tissue site matters" — here the
  depot concentration is reported but never mechanistically acted on.
- **Apixaban**: `C_APIX = CENT_APIX / V1_APIX` (µg/mL, matching the units
  `EC50_APIX_FXA` was already calibrated in).
- **Metoprolol**: `C_METRO = CENT_METRO / V1_METRO * 1000.0` (ng/mL,
  matching `EC50_METRO_HR`'s own units).

### Hill interface — rename, not a fit, for all six effect terms

Every one of the six PD effect terms (four amiodarone, one metoprolol, one
apixaban) was already shaped exactly `Emax*C/(IC50+C)` (implicit Hill
coefficient of 1) in the original — a pure rename per the guide's rule,
with `GAMMA_<STEM>_<EFFECT> = 1` pulled out as an explicit new parameter
in every case (six new parameters total, all `=1`, matching the original's
implicit shape):

| Compound | Old name(s) | New name(s) |
|---|---|---|
| AMIO | `IC50_AMIO_ERP`/`Emax_AMIO_ERP` | `EC50_AMIO_ERP`/`EMAX_AMIO_ERP`/`GAMMA_AMIO_ERP` |
| AMIO | `IC50_AMIO_HR`/`Emax_AMIO_HR` | `EC50_AMIO_HR`/`EMAX_AMIO_HR`/`GAMMA_AMIO_HR` |
| AMIO | `IC50_AMIO_IKr`/`Emax_AMIO_IKr` | `EC50_AMIO_IKR`/`EMAX_AMIO_IKR`/`GAMMA_AMIO_IKR` |
| AMIO | `IC50_AMIO_Ca`/`Emax_AMIO_Ca` | `EC50_AMIO_CA`/`EMAX_AMIO_CA`/`GAMMA_AMIO_CA` |
| METRO | `IC50_METRO_HR`/`Emax_METRO_HR` | `EC50_METRO_HR`/`EMAX_METRO_HR`/`GAMMA_METRO_HR` |
| APIX | `IC50_APIX_FXa`/`Emax_APIX_FXa` | `EC50_APIX_FXA`/`EMAX_APIX_FXA`/`GAMMA_APIX_FXA` |

Amiodarone's four effect terms (`EFFECT_AMIO_ERP`, `EFFECT_AMIO_HR`,
`EFFECT_AMIO_IKR`, `EFFECT_AMIO_CA`) are kept as four independent, named
`EFFECT_AMIO_<...>` variables rather than one combined
`EFFECT_AMIO` — the original computes them as four mechanistically
distinct pharmacology (class III/IKr, negative chronotropy, direct IKr
block, ICaL block), each read by a different downstream disease equation
(`dERP_drug`, `HR_red_combined`, `IKr_drug_block`, `Ca_AMIO_blk`
respectively), so collapsing them into one term would lose which effect
drives which physiology. This is the same reasoning already applied in
this fork to takotsubo-syndrome's carvedilol (three simultaneous receptor
blocks kept as three separate `EFFECT_CAR_*` terms).

`HR_red_combined` (metoprolol + amiodarone's HR effects combined via a
Bloch/Greco independence model, `E1+E2-E1*E2`, with an 0.82 safety
ceiling) is left exactly as the original computed it — the one point the
disease equations combine two drugs' effects, per the guide's
"combine... only at the point of use" rule; its two inputs
(`EFFECT_AMIO_HR`, `EFFECT_METRO_HR`) are the independently-named,
independently-driveable terms.

## Duplicate concentration sites — all three compounds needed normalizing, not just the one the census flagged

The census table (a corpus-wide heuristic classifier, not ground truth —
see the guide's shared rules) tagged only Amiodarone as "normalize
duplicate concentration sites, then redirect", and Apixaban/Metoprolol as
already "clean single site". Reading the actual code found duplicate
concentration computations for **all three**:

- **Amiodarone**: `$ODE` computed `Cp_AMIO = C1_AMIO/V1_AMIO` (read by all
  four PD terms) and `$TABLE` separately recomputed `Cp_AMIO_out =
  C1_AMIO/V1_AMIO` (for `$CAPTURE`) — two different names, same formula,
  the "duplicate site" pattern the census's label describes.
- **Apixaban**: `$ODE` computed `Cp_APIX_ug = C1_APIX/V1_APIX` (µg/mL,
  read by the FXa effect) and `$TABLE` separately recomputed `Cp_APIX_ng =
  C1_APIX/V1_APIX*1000.0` (ng/mL, for `$CAPTURE`) — same underlying
  concentration, two names, two units.
- **Metoprolol**: `$ODE` computed `double Cp_METRO_ng = C1_METRO/V1_METRO
  *1000.0` (read by the HR effect) and `$TABLE` computed **the exact same
  name and formula again** — this is not merely a semantic duplicate but
  a literal duplicate `double` declaration, which is one of this file's
  three pre-existing mrgsolve 2.0.1 build defects (see "Upstream build
  status" below) — the census's classifier evidently treats identical
  formulas under identical names as "clean" rather than flagging the
  compile-time collision.

All three are normalized in the refactor: each compound now computes its
`C_<STEM>` exactly once (in `$ODE`, the block the original used for the
PD-reading computation — see below), and every `$TABLE` output that used
to recompute it independently (`Cp_AMIO_out`, `Cp_APIX_ng`,
`Cp_METRO_ng`, `AntiXa_pct`) now reads that single value instead.

### Kept in the block the original used it in

Every `C_<STEM>` and `EFFECT_<STEM>` remains computed in `$ODE`, exactly
where the original computed its own equivalent — none were moved to
`$MAIN` or vice versa. This matters per the guide's warning: these values
depend on live state (`CENT_AMIO` etc.), so `$ODE` (re-evaluated every
solver substep) vs `$MAIN` (evaluated once per interval) would produce a
different simulated trajectory, not just relocated code. Only the
`$TABLE`-side *redundant recomputations* (which fed nothing but display
output) were removed and replaced with a reference to the single `$ODE`-
computed value.

### `$GLOBAL`, not `$ODE`-local `double` — avoiding both the duplicate-declaration bug and the dose-instant artifact

All nine new/promoted symbols (`C_AMIO`, `C_APIX`, `C_METRO`,
`EFFECT_AMIO_ERP`, `EFFECT_AMIO_HR`, `EFFECT_AMIO_IKR`, `EFFECT_AMIO_CA`,
`EFFECT_APIX_FXA`, `EFFECT_METRO_HR`, plus `HR_red_combined`) are declared
once at `$GLOBAL` scope and only ever *assigned* (never re-declared with
`double`) inside `$ODE`. This is the same fix the guide credits to
`clostridioides-difficile-infection` for the dose-instant stale-read
`$CAPTURE` artifact, and it is also mechanically required here: this
file's own pre-existing `Cp_METRO_ng` bug (see above) is a live
demonstration of what happens when a `$CAPTURE`d symbol is `double`-
declared more than once in this mrgsolve build. `Cp_METRO_ng` itself keeps
a single `double` declaration in `$TABLE` (now `double Cp_METRO_ng =
C_METRO;`, i.e. one declaration instead of two) rather than being promoted
to `$GLOBAL`, since nothing else in the model needs to read it.

## Upstream build status — three pre-existing mrgsolve 2.0.1 defects, none related to the compounds' own PK

`POST /model_manifest` on the untouched original's own extracted DSL
returned HTTP 500 for three independent reasons, logged as
[`translations/UPSTREAM_ISSUES.md` #128](../translations/UPSTREAM_ISSUES.md):

1. **`$CAPTURE` re-lists seventeen `$CMT` compartment names directly**
   (`AF_BURDEN ERP LAsize Fibrosis FXa_act Thrombin HR_AF Stroke_risk QTc
   NE AngII ROS SMAD23 IL6 BNP Ca_i IKr`) — mrgsolve 2.0.1 rejects this.
   Fixed by dropping them (compartments are exposed as output columns
   automatically).
2. **`Cp_METRO_ng` is `double`-declared twice** (once in `$ODE`, once in
   `$TABLE`) — see above; resolved as a side effect of normalizing
   metoprolol's own duplicate concentration site.
3. **`$PARAM IL6_0` collides with mrgsolve's auto-generated
   `<compartment>_0` init symbol** for compartment `IL6` (same defect
   class as `hashimoto-thyroiditis`'s `Th1_0` etc., `#116`). `$MAIN`'s own
   `IL6_0 = IL6_0;` line was clearly meant to seed the `IL6` compartment's
   initial value from the baseline parameter, but both are spelled
   `IL6_0`, so it self-references, and the auto-generated compartment-init
   symbol is `const`, making the assignment a second, independent compile
   error. Fixed by renaming the parameter to `IL6_BASE` (and its `$MAIN`
   and `$ODE` `IL6_decay` references) — a syntax-only rename, the
   parameter's value (1.40) is untouched.

None of these three defects sits inside any of the three refactored
compounds' own PK/PD blocks — all are pre-existing, disease-model-side
(or file-wide `$CAPTURE`-block) build defects, logged rather than fixed
upstream, and carried into `af_mrgsolve_model_refactored.R` as disclosed,
syntax-only, non-numeric fixes (confirmed by the exact-match verification
below, run against a *patched-original* baseline carrying the identical
three fixes and nothing else).

## Verification

Via the qspserver `mrgsolve_api` service (`POST /model_manifest`,
`POST /run_simulation`, `http://localhost:8007`, requests spaced ≥2.2 s
apart), comparing a patched-original baseline (the untouched original's
DSL, with only the three syntax-only fixes above applied — nothing else
changed) against `af_mrgsolve_model_refactored.R`'s extracted DSL
(confirmed byte-identical between the delivered file's quoted string and
the scratch copy used for every API call below), using **the model's own
`build_events()` dosing scenarios**, same dosing, same `$MAIN`-driven
initial conditions (the API's own `/run_simulation` respects `$MAIN`'s
`NEWIND<=1` initialization block, confirmed by the first reported row of
every run matching the original's own `init_vals`), full **730-day / 2-year**
horizon at the file's own `dt=8h` output resolution — **no shortening
needed for any scenario**, none hit the API's default solver step-count
limit:

| Compound(s) exercised | Scenario reproduced | Window | Points | Max abs diff |
|---|---|---|---|---|
| METRO alone | `S2_rate_metro` (50 mg q12h) | 0–17520 h | 2192 | **0.0** |
| AMIO alone | `S3_rhythm_amio` (400 mg q12h ×28d load → 200 mg q24h maintenance) | 0–17520 h | 2193 | **0.0** |
| APIX alone | `S4_apix_only` (5 mg q12h) | 0–17520 h | 2192 | **0.0** |
| AMIO + APIX | `S6_amio_apix` (amiodarone load+maintenance + apixaban 5 mg q12h, i.e. the file's own "standard of care" scenario) | 0–17520 h | 2194 | **0.0** |

For every row, all 30 shared `$CAPTURE` outputs (`AF_BURDEN`, `ERP`,
`LAsize`, `Fibrosis`, `FXa_act`, `Thrombin`, `HR_AF`, `Stroke_risk`,
`QTc`, `NE`, `AngII`, `ROS`, `SMAD23`, `IL6`, `BNP`, `Ca_i`, `IKr`,
`Cp_AMIO_out`, `Cp_APIX_ng`, `Cp_METRO_ng`, `Cp_AMIO_depot`, `AntiXa_pct`,
`HR_reduct_pct`, `LA_volume_mL`, `CRP_proxy`, `NTproBNP_prox`, `QTc_warn`,
`QTc_alert`, `HR_controlled`, `HR_strict_ctl`) matched exactly at every
reported time point. This is the expected outcome per the guide's
tolerance rule — every one of this file's three compounds is a pure
structural reorganisation (rename, duplicate-site normalization,
Hill-shape rename, `K12`/`K21` kept literal) with no ODE-derived
Hill-fitting anywhere, so exact (0.0) agreement is the correct result, not
a loosened tolerance.

### `/model_manifest` parameter-value diff

All 57 of the original's own parameters map to the refactored model under
their new (or unchanged) name with **zero value mismatches** (automated
diff, `KA_AMIO`↔`ka_AMIO`, `EC50_AMIO_ERP`↔`IC50_AMIO_ERP`, etc., all 57
pairs exact). The refactored manifest carries exactly 6 new parameters
beyond that — `GAMMA_AMIO_ERP`, `GAMMA_AMIO_HR`, `GAMMA_AMIO_IKR`,
`GAMMA_AMIO_CA`, `GAMMA_METRO_HR`, `GAMMA_APIX_FXA` — all `=1`, the
explicit Hill coefficients pulled out per the guide's rename rule. `$CMT`
compartment order (and therefore every compound's 1-based dosing index)
is identical between original and refactored manifests at all 24
positions (e.g. position 1 is `GI_AMIO` in the original, `GUT_AMIO` in the
refactored; position 6 is `GI_METRO`/`GUT_METRO`).

All nine new interface symbols (`C_AMIO`, `C_APIX`, `C_METRO`,
`EFFECT_AMIO_ERP`, `EFFECT_AMIO_HR`, `EFFECT_AMIO_IKR`, `EFFECT_AMIO_CA`,
`EFFECT_APIX_FXA`, `EFFECT_METRO_HR`) are state-dependent quantities
recomputed every `$ODE` step, so — per this fork's established precedent
(`takotsubo-syndrome`, `clostridioides-difficile-infection`, and others)
and confirmed directly here (see "Upstream build status" point 3: a
`$PARAM` symbol is compiled as a `const double&`, so it cannot be
reassigned inside `$ODE` at all) — they cannot be `$PARAM` defaults; all
nine, plus `Cp_METRO_ng`, appear in `/model_manifest`'s `outputPaths`
instead (46 `outputPaths` total: 24 `$CMT` states + 22 `$CAPTURE`d
quantities, vs. the original's 24 `$CMT` + 30 `$CAPTURE`, the difference
being the 17 dropped compartment-name duplicates offset by the 9 new
interface symbols).

## Delivered R file

`af_mrgsolve_model_refactored.R` is the full original script (helper
functions, six scenario definitions, all twelve figures, summary tables)
with exactly six edits, confirmed by a line-level diff of the two files
outside the quoted DSL block:

1. `af_model_code` → `af_model_code_refactored` (the DSL variable itself).
2. The quoted DSL body swapped for the refactored, renamed DSL described
   above.
3. `mcode("AF_QSP_v2", af_model_code)` → `mcode("AF_QSP_v2_refactored",
   af_model_code_refactored)` — new registered model name, to avoid any
   compiled-model cache collision with the original.
4. `init_vals`'s seven PK compartment names updated to match the renamed
   `$CMT` block (`GI_AMIO`→`GUT_AMIO`, `C1_AMIO`→`CENT_AMIO`,
   `C2_AMIO`→`PERI_AMIO`, `GI_APIX`→`GUT_APIX`, `C1_APIX`→`CENT_APIX`,
   `GI_METRO`→`GUT_METRO`, `C1_METRO`→`CENT_METRO`) — all values still 0,
   matching the original.
5–6. Three `build_events()` dosing-target comments updated to name the
   new compartments (the numeric `cmt=` indices themselves are
   **unchanged** — compartment declaration order was preserved 1-for-1,
   so `cmt=1`/`4`/`6` still target the amiodarone/apixaban/metoprolol gut
   depots respectively).

No other line in the R script (scenario definitions, plotting, summary
tables, all of which read only `$CAPTURE`d output names that were kept
identical to the original) was touched.

## Anything else flagged

- No compound other than Amiodarone, Apixaban, and Metoprolol was
  touched — this is the file's complete compound set.
- All scratch/debug files used for API verification (`.cpp` DSL
  extractions, patched-original scratch copies, request/response JSON,
  the comparison/diff scripts) were created under the session scratchpad
  directory only and deleted after use — none were left in
  `atrial-fibrillation/` or anywhere else in the repo.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`atrial-fibrillation | Amiodarone`, `atrial-fibrillation | Apixaban`, and
`atrial-fibrillation | Metoprolol (METRO)` rows.
