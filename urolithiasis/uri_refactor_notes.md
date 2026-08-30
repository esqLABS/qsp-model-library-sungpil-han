# Refactor notes — `urolithiasis/uri_mrgsolve_model.R`

Scope: **HCTZ (thiazide), ALLO/OXP (allopurinol and its own active
metabolite oxypurinol), KCIT (potassium citrate), and TAM (tamsulosin)
only** — their PK blocks and their downstream effect equations. Lumasiran
(modeled only as an on/off flag applying a fixed 53% multiplier to
hepatic oxalate synthesis, `Luma_OX_inhib = Lumasiran * 0.53`, no PK
compartment of its own) is untouched — it is not one of the five
compounds this task covers and has no PK to refactor. All 7 disease-PD
compartments (`CA_PL`, `OX_PL`, `UA_PL`, `CIT_PL`, `STONE`, `INFL`,
`GFR_CUR`) and their kinetics are untouched except for a build-
compatibility fix (see below), which touches only the `$INIT`→`$MAIN`
initial-condition idiom and two dead-clamp deletions, not any dynamics.

## OXP is allopurinol's own active metabolite, not a second independently-dosed compound

The task brief flagged this as something to check rather than assume,
and the code confirms it cleanly, the same pattern already documented in
`sarcoidosis/sarc_refactor_notes.md` for PRED/PREDL:

- `OXP` has **no `GUT_`/depot compartment and no `ka_`/`F_` absorption
  parameters of its own** anywhere in the file, and no `Dose_Oxp`
  parameter exists.
- `OXP`'s only inflow term is a conversion computed directly from
  `ALLO2` (allopurinol's own central compartment):
  ```
  dxdt_ALLO2 = ka_Allo * ALLO1 - (CL_Allo/V_Allo) * ALLO2 - k_Allo_Oxp * ALLO2;
  dxdt_OXP   = k_Allo_Oxp * ALLO2 - (CL_Oxp/V_Oxp) * OXP;
  ```
  `k_Allo_Oxp` is literally named "Allo → Oxypurinol conversion rate" in
  its own `$PARAM` comment.
- The disease-effect parameters (`Emax_XO`/`EC50_XO_Oxp`/`Hill_XO`) are
  computed against `C_OXP` (derived from `OXP`), **never** against
  `C_ALLO` — allopurinol itself has no direct pharmacodynamic effect in
  this model; only its metabolite (oxypurinol) inhibits xanthine oxidase.
  `C_ALLO` itself is computed in `$MAIN` (`ALLO2 / V_Allo`) but never
  read anywhere else in the file — a dead, informational-only quantity,
  same treatment as `PRED`'s own plasma level in the sarcoidosis file.

**Handled as a genuine parent-metabolite chain, not two independently
dosed compounds:** `ALLO` keeps its own PK block (`GUT_ALLO`/`CENT_ALLO`)
and its own `C_ALLO`, which exists only informationally (no direct PD
read anywhere in the original) — `ALLO` has no `EFFECT_ALLO`, matching
the original having no direct allopurinol PD effect. `OXP` keeps its own
single-compartment PK (`CENT_OXP`, fed by the conversion out of
`CENT_ALLO`) and exposes the `C_OXP`/`EFFECT_OXP` pair the disease
equation (xanthine-oxidase inhibition feeding `UA_prod`) actually reads.
The census row for OXP (previously `Target/Pathway: ?`) is corrected
below to reflect the parent-metabolite relationship and the real target
(xanthine oxidase), not an independent compound with an unknown target.

## Archetype per compound

**HCTZ (thiazide) — Archetype 3** (depot + central + peripheral, linear
elimination), confirmed from the actual `dxdt_` lines:
```
dxdt_HCTZ1 = Rate_HCTZ - ka_HCTZ * HCTZ1;
dxdt_HCTZ2 = ka_HCTZ * HCTZ1 - (CL_HCTZ/V1_HCTZ)*HCTZ2 - (Q_HCTZ/V1_HCTZ)*HCTZ2 + (Q_HCTZ/V2_HCTZ)*HCTZ3;
dxdt_HCTZ3 = (Q_HCTZ/V1_HCTZ)*HCTZ2 - (Q_HCTZ/V2_HCTZ)*HCTZ3;
```
Already CL/Q/V-parameterized in the original (no micro-constant
rewrite needed). Renamed `HCTZ1→GUT_HCTZ`, `HCTZ2→CENT_HCTZ`,
`HCTZ3→PERI_HCTZ`, `ka_HCTZ→KA_HCTZ`, `Emax_HCTZ→EMAX_HCTZ`
(`CL_HCTZ`/`V1_HCTZ`/`V2_HCTZ`/`Q_HCTZ`/`F_HCTZ`/`EC50_HCTZ` already
matched the convention, unchanged).

**ALLO (allopurinol) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination), PK-only, no direct PD effect
(see above):
```
dxdt_ALLO1 = Rate_Allo - ka_Allo * ALLO1;
dxdt_ALLO2 = ka_Allo * ALLO1 - (CL_Allo/V_Allo) * ALLO2 - k_Allo_Oxp * ALLO2;
```
Renamed `ALLO1→GUT_ALLO`, `ALLO2→CENT_ALLO`, `ka_Allo→KA_ALLO`,
`CL_Allo→CL_ALLO`, `V_Allo→V1_ALLO`, `F_Allo→F_ALLO`,
`Dose_Allo→Dose_ALLO`, and the conversion rate `k_Allo_Oxp→KCONV_OXP`
(named for the metabolite it produces, following the `k_predl_conv`→
`FCONV_PREDL` precedent's naming logic).

**OXP (oxypurinol) — Archetype 1** (no depot, single compartment, linear
elimination), fed by the conversion inflow above:
```
dxdt_OXP = k_Allo_Oxp * ALLO2 - (CL_Oxp/V_Oxp) * OXP;
```
Renamed `OXP→CENT_OXP`, `CL_Oxp→CL_OXP`, `V_Oxp→V1_OXP`. `C_OXP`
(`= CENT_OXP/V1_OXP`) unchanged — already matched the `C_<STEM>`
convention in the original.

**KCIT (potassium citrate) — Archetype 3 minus the peripheral
compartment** (depot + central, linear elimination):
```
dxdt_KCIT1 = Rate_KCit - ka_KCit * KCIT1;
dxdt_KCIT2 = ka_KCit * KCIT1 - (CL_KCit/V_KCit) * KCIT2;
```
Renamed `KCIT1→GUT_KCIT`, `KCIT2→CENT_KCIT`, `ka_KCit→KA_KCIT`,
`CL_KCit→CL_KCIT`, `V_KCit→V1_KCIT`, `F_KCit→F_KCIT`,
`Emax_KCit→EMAX_KCIT`, `EC50_KCit→EC50_KCIT`, `Dose_KCit→Dose_KCIT`.
`C_KCIT` unchanged (already matched convention).

**TAM (tamsulosin) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination):
```
dxdt_TAM1 = Rate_Tam - ka_Tam * TAM1;
dxdt_TAM2 = ka_Tam * TAM1 - (CL_Tam/V_Tam) * TAM2;
```
Renamed `TAM1→GUT_TAM`, `TAM2→CENT_TAM`, `ka_Tam→KA_TAM`,
`CL_Tam→CL_TAM`, `V_Tam→V1_TAM`, `F_Tam→F_TAM`, `Emax_Tam→EMAX_TAM`,
`EC50_Tam→EC50_TAM`, `Dose_Tam→Dose_TAM`. `C_TAM` unchanged (already
matched convention). Kept the stem `TAM` (not `TAMS`, as used in the
unrelated `benign-prostatic-hyperplasia/bph_mrgsolve_model_refactored.R`
tamsulosin refactor) — the guide keys naming off the compound's
*existing* abbreviation stem in *this* file (`TAM1`/`TAM2`/`C_TAM`
throughout the original), not off a different file's independent choice
for the same real drug.

## Renaming summary (all four PK blocks)

| Original | Refactored |
|---|---|
| `HCTZ1` | `GUT_HCTZ` |
| `HCTZ2` | `CENT_HCTZ` |
| `HCTZ3` | `PERI_HCTZ` |
| `ka_HCTZ` | `KA_HCTZ` |
| `Emax_HCTZ` | `EMAX_HCTZ` |
| — | `GAMMA_HCTZ = 1.0` (new; original had bare `C/(EC50+C)`) |
| `HCTZ_Ca_effect` | `EFFECT_HCTZ` |
| `ALLO1` | `GUT_ALLO` |
| `ALLO2` | `CENT_ALLO` |
| `ka_Allo` | `KA_ALLO` |
| `CL_Allo` | `CL_ALLO` |
| `V_Allo` | `V1_ALLO` |
| `F_Allo` | `F_ALLO` |
| `Dose_Allo` | `Dose_ALLO` |
| `k_Allo_Oxp` | `KCONV_OXP` |
| `OXP` | `CENT_OXP` |
| `CL_Oxp` | `CL_OXP` |
| `V_Oxp` | `V1_OXP` |
| `Emax_XO` | `EMAX_OXP` |
| `EC50_XO_Oxp` | `EC50_OXP` |
| `Hill_XO` | `GAMMA_OXP` (value 1.2 preserved, not reset to 1) |
| `XO_inhib` | `EFFECT_OXP` |
| `KCIT1` | `GUT_KCIT` |
| `KCIT2` | `CENT_KCIT` |
| `ka_KCit` | `KA_KCIT` |
| `CL_KCit` | `CL_KCIT` |
| `V_KCit` | `V1_KCIT` |
| `F_KCit` | `F_KCIT` |
| `Emax_KCit` | `EMAX_KCIT` |
| `EC50_KCit` | `EC50_KCIT` |
| `Dose_KCit` | `Dose_KCIT` |
| — | `GAMMA_KCIT = 1.0` (new; original had bare `C/(EC50+C)`) |
| `KCit_Cit_effect` | `EFFECT_KCIT` |
| `TAM1` | `GUT_TAM` |
| `TAM2` | `CENT_TAM` |
| `ka_Tam` | `KA_TAM` |
| `CL_Tam` | `CL_TAM` |
| `V_Tam` | `V1_TAM` |
| `F_Tam` | `F_TAM` |
| `Emax_Tam` | `EMAX_TAM` |
| `EC50_Tam` | `EC50_TAM` |
| `Dose_Tam` | `Dose_TAM` |
| — | `GAMMA_TAM = 1.0` (new; original had bare `C/(EC50+C)`) |
| `Tam_pass` | `EFFECT_TAM` |

`C_HCTZ`, `C_ALLO`, `C_OXP`, `C_KCIT`, `C_TAM` are unchanged — the
original already named these exactly per the `C_<STEM>` convention.
All parameter *values* are copied verbatim from the original — nothing
invented except the three `GAMMA_*` = 1.0 additions above (each a rename
target the original's bare Emax-ratio shape mathematically already
implied) and `GAMMA_OXP`, which is a pure rename of the original's own
explicit `Hill_XO = 1.2`, not a new value.

## Hill interface: rename, not a fit — for all four compounds

Every one of HCTZ, KCIT, and TAM was already the exact plain-ratio shape
`Emax*C/(EC50+C)` (implicit `gamma=1`) in the original:
```
double HCTZ_Ca_effect  = Emax_HCTZ * C_HCTZ / (EC50_HCTZ + C_HCTZ);
double KCit_Cit_effect = Emax_KCit * C_KCIT / (EC50_KCit + C_KCIT);
double Tam_pass         = Emax_Tam  * C_TAM  / (EC50_Tam  + C_TAM);
```
Rewritten to the explicit `pow(C, GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))`
form with `GAMMA_*=1.0` — mathematically identical (`pow(x,1)=x`), a pure
rename with no curve fitting performed or needed.

OXP's effect (xanthine-oxidase inhibition) was already the exact Hill
shape too, with an explicit, non-unity exponent already present in the
original:
```
double XO_inhib = Emax_XO * pow(C_OXP, Hill_XO) / (pow(EC50_XO_Oxp, Hill_XO) + pow(C_OXP, Hill_XO));
```
Renamed term-for-term (`Emax_XO→EMAX_OXP`, `EC50_XO_Oxp→EC50_OXP`,
`Hill_XO→GAMMA_OXP`, value `1.2` preserved) — also a rename, not a fit.

No R² is reported for any of the four compounds because none applies —
every effect term was already in exactly the Hill shape the guide asks
for, so this is a "pull the parameters out and rename" exercise
throughout, matching Archetypes 1 and 3's "pure structural
reorganization" tolerance rule.

## `$CAPTURE` additions for discoverability (qspserver `/model_manifest`)

The original captured `C_HCTZ_mgL`/`C_OXP_mgL`/`C_KCIT_mEqL`/`C_TAM_mgL`
(kept unchanged, now just re-pointing at the renamed `$MAIN` doubles) but
never captured `C_ALLO` or any of the four `EFFECT_<STEM>` terms — so per
the guide's qspserver-compatibility requirement #4 ("every output needed
for verification or discovery must be in `$CAPTURE`"), five new captures
were added: `C_ALLO_OUT`, `EFFECT_HCTZ_OUT`, `EFFECT_OXP_OUT`,
`EFFECT_KCIT_OUT`, `EFFECT_TAM_OUT`. The `_OUT` suffix (rather than the
bare `C_ALLO`/`EFFECT_HCTZ` name) is required because those bare names
are already declared as `$MAIN`-local `double`s used directly by `$ODE`
— a `capture C_ALLO = C_ALLO;` self-reference is the
`redefinition of capture` collision already logged as its own defect
class (`UPSTREAM_ISSUES.md` #35/#43/#50/#59/#78/#80); giving the capture
a distinct name sidesteps it entirely, the same pattern used in
`sarcoidosis/sarc_refactor_notes.md`. Confirmed discoverable via
`/model_manifest`'s `outputPaths` (91 parameters total, including every
renamed `CL_`/`V1_`/`V2_`/`Q_`/`KA_`/`F_`/`EMAX_`/`EC50_`/`GAMMA_`/
`KCONV_OXP` entry for all four compounds — see Verification).

One subtlety worth flagging for whoever next touches this file: `C_HCTZ`,
`C_ALLO`, `C_OXP`, `C_KCIT`, `C_TAM`, and all four `EFFECT_<STEM>` terms
are computed in `$MAIN`, not `$ODE`/`$TABLE`. Because the original doses
continuously (a rate term inside `dxdt_`, never an `ev()`/dosing event —
this file's own "Edge case" per the guide), `$MAIN` runs only once per
simulation rather than being re-evaluated at every reported time point,
so every `$MAIN`-computed quantity captured in `$TABLE` is one reporting
interval stale relative to the compartment amount it was derived from
(confirmed: `C_HCTZ_mgL` at `t=24h` reads `0`, using `CENT_HCTZ`'s *prior*
value, not the `t=24h` value already advanced by the solver). This is a
pre-existing property of the original's own architecture, reproduced
identically (not introduced or fixed) by this refactor — the new
`_OUT`-suffixed captures re-expose the same already-lagged `$MAIN`
doubles rather than recomputing a "fresh," non-lagged version in
`$TABLE`, specifically so the refactor's own verification stays an exact
match rather than silently changing the timing behavior.

## Build-compat fix (mrgsolve 2.0.1) — logged as `UPSTREAM_ISSUES.md` #81

Confirmed via the qspserver `mrgsolve_api` container
(`http://localhost:8007`) that the **untouched original does not compile
at all**, in three layered defects, none specific to any of the four
refactored compounds:

1. `$OMEGA @block 0` / `$SIGMA @block 0` — mrgsolve 2.0.1 does not accept
   a trailing numeric argument on `@block`; it silently reads every
   number in the block as one flat diagonal sequence (10 numbers for
   `$OMEGA`, 3 for `$SIGMA`) instead of respecting the lower-triangular
   grouping the row layout implies, which then fails validation against
   the 4/3 declared `@labels` (`Length of labels (4) does not match the
   matrix rows (10)`).
2. `$CMT`+`$INIT` jointly redeclare all 17 compartments (same class as
   issues #34/#36/etc — `$INIT`'s `NAME = value` list is its own
   compartment declaration under mrgsolve 2.0.1, not a companion to
   `$CMT`).
3. Two `$ODE` lines write directly to compartment state
   (`if(STONE<0) STONE=0;`, `if(GFR_CUR<15) GFR_CUR=15;`) — mrgsolve
   2.0.1 passes every `$CMT` state as a `const double&`, so this is a
   hard compile error, same class as issue #36.

**Fix applied directly to the delivered
`uri_mrgsolve_model_refactored.R`** (not just a scratch copy), all
syntax-only and non-numeric: (a) `$OMEGA @block 0`→`$OMEGA @block` and
`$SIGMA @block 0`→`$SIGMA` (same values, now parsing as the intended
4x4/3x3 matrices matching the declared labels); (b) the `$INIT` block
deleted, its 17 assignments moved into `$MAIN` as `<CMT>_0 = value;`
under the refactor's renamed compartments (e.g. `GUT_HCTZ_0 = 0;`,
`CENT_OXP_0 = 0;`), same values, no compartment added or removed;
(c) the two illegal clamp lines deleted outright — confirmed a no-op
even where compilable, since only `dxdt_*` feeds mrgsolve's integrator
(same reasoning as `age-related-macular-degeneration/amd_refactor_notes.md`
for issue #36). The checked-in original (`uri_mrgsolve_model.R`) is
untouched and still carries all three defects exactly as written; an
identically syntax-fixed scratch copy of the original (same three
fixes, original compartment/param names, never committed) was built
in-memory purely to construct the verification comparison target. Full
disclosure and exact error traces are in `translations/UPSTREAM_ISSUES.md`
#81.

## qspserver compatibility checklist

- `model_content` is pure mrgsolve DSL text (no R wrapper) — confirmed by
  extracting the quoted `uri_model_code <- '...'` block byte-for-byte
  and building it standalone via `POST /model_manifest`.
- Every `KA_`/`CL_`/`V1_`/`V2_`/`Q_`/`F_`/`EMAX_`/`EC50_`/`GAMMA_`/
  `KCONV_OXP`/`Dose_` parameter for all four compounds lives in `$PARAM`
  (confirmed present in the manifest's `parameters` list, 91 total).
- `C_HCTZ_mgL`/`C_OXP_mgL`/`C_KCIT_mEqL`/`C_TAM_mgL`/`C_ALLO_OUT`/
  `EFFECT_HCTZ_OUT`/`EFFECT_OXP_OUT`/`EFFECT_KCIT_OUT`/`EFFECT_TAM_OUT`
  are all discoverable `outputPaths` (see above).
- `$SET end/delta` is untouched (not present in the original either — the
  R-side `run_sim()` wrapper supplies `end`/`delta` explicitly, and the
  original doses via a continuous-rate term in `dxdt_`, never `ev()`, so
  no dosing-event schedule is baked into the DSL to begin with); the
  verification below drives `time.start`/`time.end`/`time.delta` and
  every dose-rate parameter through the request, not through any
  model-internal default.
- No R-only syntax inside the DSL block — the extracted block compiled
  standalone with no surrounding R script.
- Solver step budget: the original's own 5-year/500-point default
  (`delta = 5*365*24/500 ≈ 87.6h`) is well under the API's default
  20000-step budget for every scenario tested; no shortening was needed.

## Verification

Per the guide's mandatory protocol: extracted the quoted DSL block from
both the build-compat-patched original (scratch copy only, never
committed) and the delivered `uri_mrgsolve_model_refactored.R`, confirmed
both build via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`, confirmed healthy throughout,
requests spaced ~2s apart), then ran the original file's own six named
dosing scenarios (`Dose_HCTZ`/`Dose_Allo`/`Dose_KCit` set directly as
`parameters`, matching the original's own continuous-rate dosing design
— no `ev()`/events needed) through `POST /run_simulation`, full 5-year
(43,800h) horizon at 500 points (`delta ≈ 87.6h`), identical to the
original's own `run_sim(duration_yr=5)` default:

| Scenario | Dosing (from the original's own R scenario blocks) | Result |
|---|---|---|
| 1. Untreated CaOx stone former | all doses 0 | exact match, max abs diff 0.0 |
| 2. HCTZ 25mg/day | `Dose_HCTZ=25` | exact match, max abs diff 0.0 |
| 3. K-Citrate 60 mEq/day | `Dose_KCit=60` | exact match, max abs diff 0.0 |
| 4. Allopurinol 300mg/day (MetSyn UA stone) | `Dose_Allo=300`, MetSyn=1 | exact match, max abs diff 0.0 |
| 6. Lifestyle + HCTZ + K-Citrate | `Dose_HCTZ=25, Dose_KCit=60` | exact match, max abs diff 0.0 |
| 7. MetSyn UA stone — combo Rx | `Dose_Allo=300, Dose_KCit=60`, MetSyn=1 | exact match, max abs diff 0.0 |

Every shared output was compared point-by-point across all 501 time
points: all 17 compartments (mapped through the renames given above) and
19 deterministic `$CAPTURE` outputs (`uCit_mgday`, `uVol_Lday`,
`urinepH`, `SS_CaOx_out`, `SS_CaP_out`, `SS_UA_out`, `StoneSize_mm`,
`GFR_out`, `Inflammation`, `C_HCTZ_mgL`, `C_OXP_mgL`, `C_KCIT_mEqL`,
`C_TAM_mgL`, `Hypercalciuria`, `Hyperoxaluria`, `Hyperuricosuria`,
`Hypocitraturia`, `LowUV`, `Stone_recur_risk`). **Result: exact match,
max abs diff = 0.0, for every one of these outputs in all six
scenarios** — the expected result for pure structural reorganization
with no Hill-refitting (every effect term was already an exact Hill
shape, see above).

**Excluded from the exact-match check, and why:** `uCa_mgday`,
`uOx_mgday`, `uUA_mgday` carry the original's own `$SIGMA`-drawn residual
error (`capture uCa_mgday = uCa_24 * (1.0 + EPS(1));`, etc.) and are
therefore genuinely stochastic from one API call to the next — confirmed
by calling the identical untouched-original model twice with identical
parameters and observing different `uCa_mgday` trajectories while
`StoneSize_mm` (and every other non-`EPS`-touched output) stayed
identical between the two calls. This is a pre-existing property of the
original's own design (residual observation error on three specific
urinary outputs), not something the refactor introduced or could
control; the *pre-noise* deterministic quantities these three are
derived from (`uCa_24`/`uOx_24`/`uUA_24`) are exercised indirectly via
the exactly-matching `Hypercalciuria`/`Hyperoxaluria`/`Hyperuricosuria`
threshold flags (which read `uCa_24` etc. directly, not the noised
capture), so the underlying deterministic math is still fully verified.

**Tamsulosin bespoke verification (none of the original's own six
scenarios doses it):** grepping the original's own scenario blocks
confirms `Dose_Tam` is never set to a nonzero value anywhere in the
file — tamsulosin's PK/PD is modeled but never actually exercised by any
scenario the author wrote or plotted. To verify it anyway, a bespoke run
was constructed using the same baseline diet parameters as Scenario 1
plus `Dose_Tam=400` (mg/day; an arbitrary nonzero value, not implying
this is a clinically realistic dose — tamsulosin's real-world dose is
0.4mg, but the original file's own dose parameter is unit-labeled
mg/day and was never populated with any specific value by the author, so
no "original" dose exists to reproduce) — this is not an invented dosing
*mechanism*, only a nonzero value for an existing, already-declared
`$PARAM` knob the original's own scenarios simply left at its `0`
default. Full 5-year horizon, same as above. Result: `GUT_TAM`/`CENT_TAM`
(mapped from `TAM1`/`TAM2`) and `C_TAM_mgL` matched exactly (max abs diff
0.0) between the build-fixed original and the refactored model.
`EFFECT_TAM_OUT` has no original counterpart to compare against (the
original never captured `Tam_pass`) — spot-checked instead by manual
recomputation: at the final time point, `C_TAM_mgL = 4.2857`,
`EMAX_TAM*C/(EC50_TAM+C) = 0.28*4.2857/(0.001+4.2857) = 0.279935`,
matching the model's own reported `EFFECT_TAM_OUT = 0.2799` (rounded to
the API's 4-decimal output precision).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
five `urolithiasis` rows (ALLO, HCTZ, KCIT, Oxypurinol renal (OXP),
Tamsulosin) — OXP's `Target/Pathway` corrected from `?` to "Xanthine
oxidase (as allopurinol's own active metabolite)" per the
parent-metabolite finding above; the other four rows' `Target/Pathway`
filled in from the actual code (NCC for HCTZ, none/PK-only for ALLO,
urinary citrate excretion for KCIT, alpha1A-adrenergic receptor /
ureteral smooth muscle relaxation for TAM).

## Anything else flagged

- Lumasiran is out of this refactor's scope (not one of the five named
  compounds, and it has no PK compartment of its own to refactor — it is
  a bare on/off flag applying a fixed 53% multiplier to hepatic oxalate
  synthesis). Untouched.
- The post-DSL R script required two small, disclosed edits to keep the
  file runnable standalone after the rename: Scenario 3/6's
  `Dose_KCit = 60` → `Dose_KCIT = 60`, and Scenario 4/7's
  `Dose_Allo = 300` → `Dose_ALLO = 300` (both simple `param()` list
  keys, not events or `cmt=` targets, since this file doses via
  continuous rate parameters rather than `ev()`). Everything else in the
  post-DSL script — `run_sim()`, all seven scenario definitions' other
  parameters, all four ggplot figures (none of their aesthetics reference
  a renamed compartment or capture — every column they plot,
  `StoneSize_mm`/`SS_CaOx_out`/`uCa_mgday`/`uCit_mgday`/`uOx_mgday`/
  `GFR_out`/`Inflammation`, is unchanged by this refactor), and the fluid-
  intake sensitivity analysis — is byte-identical to the original.
- `mcode()`'s model identifier was changed from `"urolithiasis_qsp"` to
  `"urolithiasis_qsp_refactored"` purely to avoid an `mrgsolve` `soloc`
  cache collision if both files are ever compiled in the same R session;
  non-behavioral.
