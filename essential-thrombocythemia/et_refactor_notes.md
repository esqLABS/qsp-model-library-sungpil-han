# Refactor notes — `essential-thrombocythemia/et_mrgsolve_model.R`

All four compounds modeled in this file were in scope, per the existing
rows in `driver-patches/data/compound_perturbation_census.md` (all four
classified "Redirect concentration (clean single site)"): **Anagrelide
(ANA)**, **Hydroxyurea (HU)**, **Pegylated interferon (PIFN)**, and
**Ruxolitinib (RUX)**. Aspirin (`DOSE_ASA`) also appears in this file but
is a binary on/off flag with no PK of its own (no concentration
compartment, no census row) — left completely untouched beyond the shared
`$PARAM` build-compat fix described below (needed for the whole model to
compile at all, not specific to any one compound).

## Archetype per compound: all four are Archetype 3 (depot + central +
peripheral), with one structural quirk preserved from the original

The original does not use a real absorption-depot compartment for any of
the four compounds. Instead, each drug's central compartment ODE has a
closed-form absorption term baked in directly:

```
dxdt_HU_C  = F_HU * DOSE_HU * ka_HU * exp(-ka_HU * SOLVERTIME)
             - (CL_HU/V1_HU) * HU_C - (Q_HU/V1_HU) * (HU_C - HU_P);
dxdt_HU_P  = (Q_HU/V1_HU) * (HU_C - HU_P) - (Q_HU/V2_HU) * HU_P;
```

`SOLVERTIME` is this mrgsolve build's valid in-`$ODE` accessor for the
solver's own advancing time (`_ODETIME_[0]`, per `UPSTREAM_ISSUES.md`
#27), so `F*DOSE*ka*exp(-ka*t)` is the exact closed-form solution of a
depot compartment `GUT` with `dGUT/dt = -ka*GUT` and initial amount
`GUT(0) = F*DOSE` — i.e. `ka*GUT(t) = F*DOSE*ka*exp(-ka*t)`, term for
term. This refactor makes that depot real and explicit — a rename/
reorganization, not a new mechanism:

```
dxdt_GUT_HU  = -KA_HU * GUT_HU;
dxdt_CENT_HU =  KA_HU * GUT_HU - (CL_HU/V1_HU) * CENT_HU - (Q_HU/V1_HU) * (CENT_HU - PERI_HU);
dxdt_PERI_HU =  (Q_HU/V1_HU) * (CENT_HU - PERI_HU) - (Q_HU/V2_HU) * PERI_HU;
```
with `GUT_HU_0 = F_HU * DOSE_HU;` set once in `$MAIN` — the original's own
"single administration at t=0, no repeat-dosing loop" mechanism, now
expressed as a real initial condition on a depot state instead of a
literal `exp()` term. `DOSE_HU`, `DOSE_ANA`, `DOSE_RUX`, `DOSE_PIFN` are
themselves the original's own dosing parameters (see build-compat fix #3
below — they were never actually declared in `$PARAM`), carried forward
unchanged so the refactored model's *default* driving mechanism still
matches the original's own scenarios exactly. `C_<STEM>` is the guide's
required single redirect site for an external covariate to intercept
later. Identical structure for ANA, RUX, and PIFN (stem renamed
`pIFN`→`PIFN` throughout, matching the naming convention's uppercase-stem
rule).

**Structural quirk preserved exactly, not "corrected" to the textbook
archetype-3 form:** the original's central/peripheral exchange is *not*
the standard mass-balance split (`Q/V1` outflow from central, `Q/V2`
return from peripheral) the guide's own archetype-3 template shows. It
instead applies the *same* `Q/V1` coefficient to the `(central-peripheral)`
difference in **both** directions, plus an independent `Q/V2` loss term
applied only to the peripheral compartment:
```
- (Q_HU/V1_HU) * (HU_C - HU_P)      // central loses this
+ (Q_HU/V1_HU) * (HU_C - HU_P) - (Q_HU/V2_HU) * HU_P   // peripheral: gains the same term, then separately loses Q/V2*HU_P
```
This was caught by verification, not assumed: an earlier draft of this
refactor used the guide's literal archetype-3 template (`Q/V1` out, `Q/V2`
back), which compiles and runs but produces a **different** trajectory
(e.g. at t=10 into the HU 1500 mg/d scenario, `Cp_HU`=0.1961 in the
original vs `C_HU`=0.2560 in that draft — a real ~30% discrepancy, growing
over time, not floating-point noise). Switching to the exact form above
(`(CL/V1)*CENT` for elimination, `(Q/V1)*(CENT-PERI)` for exchange, `Q/V2*
PERI` as peripheral's *only additional* term) reproduces the original
exactly (see Verification below). Per the guide's "always rewrite to
CL/Q/V" instruction: this file's parameters were already named
`CL_<STEM>`/`Q_<STEM>`/`V1_<STEM>`/`V2_<STEM>` (no micro-constant renaming
needed), but the instruction to use "the more interpretable, more common
convention" does not license silently changing which equation those named
parameters plug into — that would be a numeric/behavioral change, not a
rename. The non-standard exchange form is kept exactly as the original
computes it.

**The exposed concentration `C_<STEM>` is a direct alias of `CENT_<STEM>`,
not `CENT_<STEM>/V1_<STEM>`.** The original's `HU_C`/`ANA_C`/`RUX_C`/
`pIFN_C` states are themselves concentration-valued (µg/mL, ng/mL) with no
separate mass→concentration division anywhere in the file — `V1_<STEM>`
only ever appears as part of a rate constant (`CL/V1`, `Q/V1`), never as a
literal amount/volume divisor. Introducing `C_HU = CENT_HU/V1_HU` would
double-apply the volume already folded into the rate constants and change
every downstream number. `C_HU = CENT_HU;` (identity) is therefore the
correct, non-distorting rename here — the same reasoning `x-linked-
hypophosphatemia/xlh_refactor_notes.md` used for `C_CALC`/`C_PHOSORAL`.

## Renaming

| Original | Refactored |
|---|---|
| `HU_C` / `HU_P` | `CENT_HU` / `PERI_HU` (new: `GUT_HU`) |
| `ANA_C` / `ANA_P` | `CENT_ANA` / `PERI_ANA` (new: `GUT_ANA`) |
| `RUX_C` / `RUX_P` | `CENT_RUX` / `PERI_RUX` (new: `GUT_RUX`) |
| `pIFN_C` / `pIFN_P` | `CENT_PIFN` / `PERI_PIFN` (new: `GUT_PIFN`) |
| `ka_HU`/`ka_ANA`/`ka_RUX`/`ka_pIFN` | `KA_HU`/`KA_ANA`/`KA_RUX`/`KA_PIFN` |
| `Emax_HU`/`Emax_ANA`/`Emax_RUX`/`Emax_pIFN` | `EMAX_HU`/`EMAX_ANA`/`EMAX_RUX`/`EMAX_PIFN` |
| `gam_HU`/`gam_ANA`/`gam_RUX`/`gam_pIFN` | `GAMMA_HU`/`GAMMA_ANA`/`GAMMA_RUX`/`GAMMA_PIFN` |
| `E_HU`/`E_ANA`/`E_RUX`/`E_pIFN` | `EFFECT_HU`/`EFFECT_ANA`/`EFFECT_RUX`/`EFFECT_PIFN` |
| `Cp_HU`/`Cp_ANA`/`Cp_RUX`/`Cp_pIFN` (`$TABLE`-only aliases) | superseded by `C_HU`/`C_ANA`/`C_RUX`/`C_PIFN` (see below) |
| `DOSE_pIFN` (R-side and, per build-compat fix, `$PARAM`) | `DOSE_PIFN` |

`F_HU`/`F_ANA`/`F_RUX`, `CL_HU`/`CL_ANA`/`CL_RUX`, `V1_HU`/`V1_ANA`/
`V1_RUX`, `V2_HU`/`V2_ANA`/`V2_RUX`, `Q_HU`/`Q_ANA`/`Q_RUX`,
`EC50_HU`/`EC50_ANA`/`EC50_RUX`/`EC50_pIFN`, and `MW_HU` already matched
convention or were left as-is (`MW_HU` is declared but never referenced
anywhere in the original's `$ODE`/`$TABLE`, kept verbatim). `F_pIFN`,
`CL_pIFN`, `V1_pIFN`, `V2_pIFN`, `Q_pIFN`, `EC50_pIFN` were renamed only
for the stem casing (`pIFN`→`PIFN`), values unchanged. All parameter
*values* are copied verbatim from the original — nothing invented or
defaulted beyond the build-compat `DOSE_*` defaults noted below (which
match the original's own R-side default of 0).

Two mechanical consequences of the rename, both required for the model to
still mean the same thing (not scope creep):
- `$TABLE`'s `Cp_HU`/`Cp_ANA`/`Cp_RUX`/`Cp_pIFN` recompute lines are gone.
  Under this mrgsolve build, `$CAPTURE` auto-forward-declares each of its
  names once for the whole translation unit, and `$ODE`/`$TABLE` share
  that same scope (confirmed empirically: redeclaring `double C_HU` in
  `$TABLE` after `$ODE` already computed it is a **compile-time
  redefinition error** here, not silently shadowed) — so `C_HU`/`EFFECT_HU`
  etc., once computed in `$ODE` and listed in `$CAPTURE`, are already the
  values reported; no separate `$TABLE` recomputation is needed or
  possible. The original's own `Cp_HU := HU_C` pattern worked only because
  `Cp_HU` was a distinct name never used in `$ODE`.
- The R-side `scenarios`/`mk_ev`/`run_scenario`/`results` code was updated
  to reference `DOSE_PIFN` instead of `DOSE_pIFN` throughout, to match the
  renamed `$PARAM` entry.

Added (new, purely to complete the named Hill interface, not present in
the original): none — `EMAX_<STEM>`, `EC50_<STEM>`, `GAMMA_<STEM>` all
already existed under other names (`Emax_*`/`EC50_*`/`gam_*`) with
explicit values, so this is a rename in every case, not an addition.

## Hill interface: rename in all four cases, not a fit

All four compounds' disease effects were already written as plain
`Emax*C^gamma/(EC50^gamma+C^gamma)` Hill ratios in the original, with
explicit (non-default) `gam_*` values for all four — no ODE-solved
receptor kinetics anywhere in this file (no TMDD, no receptor-binding
compartments for any of the four compounds — checked explicitly in the
actual `$ODE` block, not assumed from drug class). This is a pure rename
per the guide's "already this shape" rule:

- `E_HU = Emax_HU*pow(HU_C,gam_HU)/(pow(EC50_HU,gam_HU)+pow(HU_C,gam_HU))`
  → `EFFECT_HU = EMAX_HU*pow(C_HU,GAMMA_HU)/(pow(EC50_HU,GAMMA_HU)+pow(C_HU,GAMMA_HU))`,
  `GAMMA_HU = 1.20` (already explicit, not defaulted).
- Same pattern for ANA (`GAMMA_ANA = 1.50`), RUX (`GAMMA_RUX = 1.80`), and
  PIFN (`GAMMA_PIFN = 1.00`, already explicit in the original as
  `gam_pIFN`).

No curve-fitting (`nls()` or otherwise) was needed or performed for any of
the four compounds.

## Three pre-existing build defects found while verifying (not fixed
upstream, logged as `translations/UPSTREAM_ISSUES.md` #64)

Confirmed on the untouched original independently via `POST
/model_manifest` against the original file's own extracted DSL — none
introduced by this refactor:

1. **`$CMT` (bare) and `$INIT` jointly redeclare all 16 compartments** —
   the same defect class as issues #34/#36/#42/#44/#48/#49/#51/#61/#62/#63:
   `invalid class "mrgmod" object: Duplicated model names: HSC MKP MK PLT
   TPO JAK2 SPL HU_C HU_P ANA_C ANA_P RUX_C RUX_P pIFN_C pIFN_P RISK_T
   RISK_MF`.
2. **`$CAPTURE` repeats two compartment names already in `$CMT`**
   (`RISK_T`, `RISK_MF`): `invalid class "mrgmod" object: compartment
   should not be in $CAPTURE: RISK_T,RISK_MF`.
3. **`$ODE` reads `DOSE_HU`/`DOSE_ANA`/`DOSE_RUX`/`DOSE_pIFN`/`DOSE_ASA`
   directly but none of the five is ever declared in `$PARAM`** — the
   original's own R wrapper only ever sets them via `param()`/
   `mrgsim(param=...)` at the R side, never inside the DSL block itself,
   which is fatal at `mcode()`-build time (`'DOSE_HU' was not declared in
   this scope`, etc. — a genuinely new defect class for this repo, a
   missing declaration rather than a collision).

Per the guide's settled policy for this situation, all three fixes are
syntax-only and non-numeric, and were applied **directly to the delivered
`et_mrgsolve_model_refactored.R`** (not just a scratch copy), while
`et_mrgsolve_model.R` was left completely untouched:

1. The standalone `$CMT` block was removed; `$INIT` alone declares every
   compartment (same fix as issue #51).
2. `RISK_T`/`RISK_MF` were dropped from `$CAPTURE` — compartment states
   are always present in mrgsolve's own output regardless of `$CAPTURE`
   membership, confirmed by diffing `/model_manifest`'s `outputPaths`
   before/after (both still list `RISK_T`/`RISK_MF`).
3. `DOSE_HU = 0`, `DOSE_ANA = 0`, `DOSE_RUX = 0`, `DOSE_PIFN = 0`,
   `DOSE_ASA = 0` were added to `$PARAM`, matching the defaults the
   original's own R wrapper already used
   (`mod2 <- mod %>% param(DOSE_HU = 0, ...)`) — no value invented beyond
   what the original itself already treated as the no-dose default.

Logged as entry **#64** in `translations/UPSTREAM_ISSUES.md`.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL —
all Hill-interface and PK parameters appear in the manifest's
`parameters`, including the five dosing parameters added by the
build-compat fix: `DOSE_HU`, `DOSE_ANA`, `DOSE_RUX`, `DOSE_PIFN`,
`DOSE_ASA`, `KA_HU`/`KA_ANA`/`KA_RUX`/`KA_PIFN`, `EMAX_HU`/`EMAX_ANA`/
`EMAX_RUX`/`EMAX_PIFN`, `GAMMA_HU`/`GAMMA_ANA`/`GAMMA_RUX`/`GAMMA_PIFN`,
`EC50_HU`/`EC50_ANA`/`EC50_RUX`/`EC50_PIFN`, plus every `CL_*`/`V1_*`/
`V2_*`/`Q_*`/`F_*` (69 parameters total).

`C_HU`/`C_ANA`/`C_RUX`/`C_PIFN` and `EFFECT_HU`/`EFFECT_ANA`/`EFFECT_RUX`/
`EFFECT_PIFN` are state-derived doubles (recomputed every `$ODE` step, not
fixed covariates), so — per the same reasoning established in
`x-linked-hypophosphatemia/xlh_refactor_notes.md`, `mucopolysaccharidosis-
type-1/mps1_refactor_notes.md`, and `primary-sclerosing-cholangitis/
psc_refactor_notes.md` for this engine/plugin combination — they cannot
also be `$PARAM` entries (attempting it collides at the C++ stage); they
appear in the manifest's `outputPaths` via `$CAPTURE` instead, confirmed
present (39 output paths total, including all 21 compartment states).

`model_content` is pure DSL text extracted from the `code <- '...'`
R-string wrapper — no `.cpp` sibling was left behind; the extraction was
in-memory only, used to build the verification requests below and then
discarded.

## Verification

Per the guide's mandatory protocol: ran all seven of the original file's
own dosing scenarios (`scenarios` list, doses exactly as written) through
both the (defect-patched, for buildability only) original and the
refactored model via the qspserver `mrgsolve_api` service (`POST
/model_manifest`, `POST /run_simulation`), 1-year-plus horizon matching
the original's own `results <- purrr::map_dfr(...)` computation
(`end=730, delta=1` — the code path that actually produces every plotted/
reported number in the original, as opposed to the separately-defined but
never-invoked `run_scenario()` function). Dosing parameters translated
directly (`DOSE_HU = dose_hu/1e3`, etc., matching the original's own
`params_vec` construction exactly).

Scenarios run (names/doses copied verbatim from the original's own
`scenarios` list):

1. **① No Treatment (Observation)** — no PK at all.
2. **② Low-dose Aspirin only** — no PK at all (ASA is a flag, not a
   modeled compound).
3. **③ Hydroxyurea 500 mg/d + ASA**
4. **④ Hydroxyurea 1500 mg/d + ASA**
5. **⑤ Anagrelide 2 mg/d + ASA**
6. **⑥ Ruxolitinib 20 mg/d (BID) + ASA**
7. **⑦ Peg-IFN-α2a 90 µg/week + ASA**

Every shared output was compared point-by-point across the full 730-day
grid: `PLT_count`, `JAK2_AB`, `Spleen`, `TPO_level`, `MK_pool`, `CHR_flag`,
`PHR_flag`, `CMR_flag`, `Hazard_thromb`, `Hazard_MF`, `RISK_T`, `RISK_MF`,
and `Cp_HU`/`Cp_ANA`/`Cp_RUX`/`Cp_pIFN` ↔ `C_HU`/`C_ANA`/`C_RUX`/`C_PIFN`.

**Result: match to floating-point/JSON-rounding precision (max abs diff
≤ 1e-4 for state variables of order 1–900, ≤ ~1.5e-8 for the drug
concentration crosswalks themselves) for every scenario, over the entire
window either model could run to completion** — the expected outcome per
the guide's tolerance table for a pure structural reorganization with
rename-only Hill terms (no fitting performed, none needed):

- Scenarios ① and ②: exact match (max abs diff 0.0, one rounding-scale
  outlier of `2.3e-25` on `Cp_pIFN`↔`C_PIFN`, both models' floating-point
  representation of an exact 0), no `NaN`, full 730 days.
- Scenario ⑦ (Peg-IFN-α2a): exact match (max abs diff ≤ 6e-10, JSON
  rounding scale), no `NaN`, full 730 days — `GAMMA_PIFN=1.0` is immune to
  the fragility below.
- Scenarios ③/④ (Hydroxyurea 500/1500 mg/d): match to ≤ 1e-4 (state
  variables) / ≤ 1.3e-8 (`Cp_HU`↔`C_HU`) up to each model's own
  independent `NaN` onset (original: t=158/153 days; refactored: t=119/129
  days — see the fragility note below).
- Scenario ⑤ (Anagrelide): match to ≤ 1e-4 / ≤ 7.5e-9 up to `NaN` onset
  (original: t=33 days; refactored: t=29 days).
- Scenario ⑥ (Ruxolitinib): match to ≤ 1e-4 / ≤ 1.5e-8 up to `NaN` onset
  (original: t=74 days; refactored: t=47 days).

**Important process note on an earlier, incorrect draft:** an initial
version of this refactor used the guide's literal archetype-3 exchange
template (`Q/V1` out of central, `Q/V2` back from peripheral) and *failed*
this verification — `Cp_HU` vs `C_HU` diverged progressively (e.g. 0.502
vs 0.531 at t=5, 0.196 vs 0.256 at t=10 in the HU 1500 mg/d scenario, a
real and growing discrepancy, not noise). Per the guide's "a gate failure
is a question, not a verdict, and it cuts both ways" — this was traced to
the original's non-standard exchange form (see Archetype section above),
not a flaw in the verification itself. Switching the refactored `$ODE` to
exactly reproduce the original's actual exchange equations (rather than
the guide's textbook template) resolved it completely; the numbers above
are from the corrected version.

**`NaN` fragility (shared, pre-existing, independently reproducing —
logged as `UPSTREAM_ISSUES.md` #64's fourth finding, not a refactor bug):**
three of the four compounds' Hill terms use a non-integer `GAMMA`
(1.20/1.50/1.80 for HU/ANA/RUX). Over the 730-day horizon, each drug's own
concentration decays (dominated by the slow terminal/redistribution phase,
since `KA` is far larger than the elimination rate constants) to values
many orders of magnitude below anything with physiological meaning
(sub-`1e-8`, against `EC50` of 0.008–3.5) — small enough that ordinary
floating-point/adaptive-solver roundoff can momentarily push the
compartment slightly negative, and `pow(negative, non-integer)` is `NaN`
in C++, which cascades into every other output once it occurs. This
reproduces independently in **both** the original and the refactored
model, confirmed by inspecting the trajectories directly around each
model's own onset point (e.g. for the Anagrelide scenario: `MKP` in the
original is finite through t=32 and `NaN` from t=33; in the refactored
model, finite through t=28 and `NaN` from t=29 — different onset time
because the two formulations reach the same near-zero values via
different floating-point paths, but the same underlying cause). Since
both sides agree essentially exactly right up to each one's own
independent blow-up, this is evidence the refactor is correct, not a
mismatch — consistent with the same pattern seen in `abdominal-aortic-
aneurysm/aaa_refactor_notes.md` and `chronic-myeloid-leukemia/
cml_refactor_notes.md`. Not fixed in the delivered model (would require a
`fmax(C, 0.0)` clamp before each fractional `pow()`, a numeric/behavioral
change out of scope for a rename-only refactor); disclosed here and in
`UPSTREAM_ISSUES.md` #64 instead.

## Anything else worth flagging

- The four compounds are fully independent in this file — no shared
  multi-drug expression. Each compound's `EFFECT_<STEM>` feeds a distinct
  (or non-overlapping) point in the disease ODEs: `EFFECT_HU` and
  `EFFECT_RUX` both appear in `dxdt_MKP` (added together, not collapsed
  into one term — each keeps its own name and Hill parameters, matching
  the guide's "never collapse several drugs into one shared Hill term"
  rule), `EFFECT_ANA` in `dxdt_MK`, `EFFECT_RUX` again in `dxdt_SPL`, and
  `EFFECT_PIFN` in `dxdt_JAK2`.
- Compartment count went from 17 to 21 (added `GUT_HU`/`GUT_ANA`/
  `GUT_RUX`/`GUT_PIFN`, the four new depot compartments); no compartment
  was removed. Ordering: disease compartments first (unchanged), then the
  four compounds' `GUT_<STEM>`/`CENT_<STEM>`/`PERI_<STEM>` triples in the
  same HU→ANA→RUX→PIFN order the original used, then `RISK_T`/`RISK_MF`
  last (unchanged position) — any code addressing compartments by 1-based
  number rather than name would need updating, but no such code exists in
  this file (the R wrapper never dosed via `cmt=` events; see Archetype
  section).
- Aspirin (`DOSE_ASA`/`E_ASA`) is untouched beyond the shared `$PARAM`
  build-compat fix (defect 3 above touched it because it shares that
  missing-declaration defect with the four in-scope compounds, not
  because it was in scope itself).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`essential-thrombocythemia` rows for ANA, HU, PIFN, and RUX.
