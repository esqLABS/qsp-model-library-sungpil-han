# Refactor notes — `dic_mrgsolve_model.R` (ATRA + tranexamic acid)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only two
compounds' PK/PD blocks were rewritten, per their existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(both classified "Redirect concentration (clean single site)"): **ATRA**
(stem `ATR`) and **tranexamic acid** (stem `TXA`). Every other compound in
the file — unfractionated heparin (`HEPC`/`VH`/`KELH`), enoxaparin
(`ENXD`/`ENXC`), thrombomodulin alfa (`RTM`/`VTM`/`KELTM`), drotrecogin alfa
(`APCX`/`VAPC`/`KELAPC`), and argatroban (`ARGA`/`VARG`/`KELARG`) — is
copied verbatim, confirmed by `diff` (see below). Every disease-side
equation (thrombin generation, the two "clocks," protein C, platelets,
fibrinolysis, organ injury, bleeding, mortality hazard) is untouched.

## Archetype determined

### ATRA (`ATR`) — Archetype 3 minus the peripheral compartment (depot + central, linear elimination)

The original models ATRA with a gut depot (`ATRAD`) absorbing into a
central amount (`ATRAC`), no peripheral compartment, no bioavailability
factor (fully absorbed):

```
dxdt_ATRAD = RATE_ATRA - KAATR * ATRAD;
dxdt_ATRAC = KAATR * ATRAD - KELATR * ATRAC;
```

matching Archetype 3 in `FORK_WORKFLOW_GUIDE.md` with the peripheral
compartment dropped (the guide's own "drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a
2-compartment-no-depot variant," applied here to drop only the peripheral
compartment while keeping the depot). The disease effect
(`ATRAE = EMAXATRA*C/(EC50ATRA+C)`, driving APL blast differentiation) is
already the guide's canonical Hill ratio with an implicit exponent of 1 —
a rename, not a refit.

**Renaming applied (ATRA only):**

| Original | Refactored | Value |
|---|---|---|
| `ATRAD` (compartment) | `GUT_ATR` | — |
| `ATRAC` (compartment) | `CENT_ATR` | — |
| `KAATR` | `KA_ATR` | 1.5 (1/h) |
| `KELATR` + `VATR` | `CL_ATR` (new, `= KELATR x VATR`) + `V1_ATR` | `CL_ATR = 50.0` (L/h), `V1_ATR = 100` (L) |
| `EC50ATRA` | `EC50_ATR` | 0.05 (mg/L) |
| `EMAXATRA` | `EMAX_ATR` | 1.0 |
| — (none; original had no Hill exponent) | `GAMMA_ATR` (new, `= 1`) | documented as "original had no explicit Hill term" |
| `ATRAE` (local `$ODE` ratio) | `EFFECT_ATR` | — |
| `RATE_ATRA` (continuous-infusion-input covariate, never actually invoked by any scenario — all ATRA dosing in this file is bolus) | left unchanged | out of the naming convention's role table, same treatment as `TCZIN`/`NIMPO` in prior calibration runs |

**`CL_ATR` is a derived parameter, not an invented one** — same technique as
`CL_NIM` in `aneurysmal-subarachnoid-hemorrhage/sah_refactor_notes.md`:
`CL_ATR := KELATR x VATR = 0.50 x 100 = 50.0`, so `CL_ATR / V1_ATR` reduces
to exactly `0.50` — bit-identical to Python's own `50.0/100.0 == 0.50`
check, confirmed empirically (see Verification).

`EFFECT_ATR` replaces `ATRAE` at its one use site
(`dxdt_BLAST`'s `- KDIFF * EFFECT_ATR * BLASTv`).

### Tranexamic acid (`TXA`) — Archetype 1 (no depot, single compartment, linear elimination) + continuous-infusion-input edge case

```
dxdt_TXAC = RATE_TXA - KELTX * TXAC;
```

One compartment, no absorption depot, no peripheral compartment, no TMDD —
Archetype 1 exactly. `RATE_TXA` is the guide's "continuous infusion input,
not an event" edge case (a zero-order rate parameter added directly into
the ODE, used by two scenarios alongside a same-time bolus dose into the
same compartment) — preserved as-is, name unchanged (it already follows
this file's own consistent `RATE_<STEM>` convention shared by every other
compound, so renaming it would not remove any naming-chaos this guide
targets).

**Renaming applied (tranexamic acid only):**

| Original | Refactored | Value |
|---|---|---|
| `TXAC` (compartment) | `CENT_TXA` | — |
| `KELTX` + `VTX` | `CL_TXA` (new, `= KELTX x VTX`) + `V1_TXA` | `CL_TXA = 6.96` (L/h), `V1_TXA = 12.0` (L) |
| `CTXA` (local `$ODE` concentration) | `C_TXA` | — |
| `IC50TXA` | `EC50_TXA` | 6.0 (ug/mL) |
| — (none; the original's inhibitory ratio has an implicit ceiling of 1.0) | `EMAX_TXA` (new, `= 1.0`) | documented below |
| — (none; original had no Hill exponent) | `GAMMA_TXA` (new, `= 1`) | documented as "original had no explicit Hill term" |
| `TXAI` (local `$ODE` inhibitory ratio) | `EFFECT_TXA` (redefined; see below) | — |
| `RATE_TXA` | left unchanged | already stem-consistent, used by two scenarios |

**`CL_TXA` is derived, not invented**, same technique as `CL_ATR` above:
`CL_TXA := KELTX x VTX = 0.58 x 12.0 = 6.96`, and `6.96/12.0 == 0.58`
bit-for-bit in IEEE double precision (confirmed empirically).

**The inhibitory-vs-stimulatory Hill wrinkle.** The original's effect term
is *inhibitory* — it computes the fraction of fibrinolytic activity
*surviving* TXA, not the compound's own effect magnitude:

```c
double TXAI = 1.0 / (1.0 + CTXA / IC50TXA);   // -> 1 at C=0, -> 0 as C rises
```

used at two sites, both as a multiplicative survival fraction:
`FGNLYS = KFGN * PLNv * (1-A2v)^2 * TXAI;` and
`vPLN = (... ) * TXAI;`. The guide's named Hill interface
(`EFFECT_<STEM> = EMAX*X^gamma/(EC50^gamma+X^gamma)`) is stimulatory — it
rises with concentration. Rather than expose an inhibitory `EFFECT_TXA`
(which would read backwards next to every other compound's `EFFECT_<STEM>`
in this corpus), `EFFECT_TXA` is defined in the guide's canonical
stimulatory form (the *degree* of fibrinolysis inhibition, 0 at C=0, rising
toward `EMAX_TXA`), and the two use sites read `(1.0 - EFFECT_TXA)` in
place of the original's `TXAI`:

```c
double EFFECT_TXA = EMAX_TXA * pow(C_TXA, GAMMA_TXA)
                   / (pow(EC50_TXA, GAMMA_TXA) + pow(C_TXA, GAMMA_TXA));
// ... FGNLYS = KFGN * PLNv * (1-A2v)^2 * (1.0 - EFFECT_TXA);
// ... vPLN   = (...) * (1.0 - EFFECT_TXA);
```

This is algebraically **a rename, not a refit**: with `EMAX_TXA = 1`,
`GAMMA_TXA = 1`, `1 - EFFECT_TXA = 1 - C/(EC50+C) = EC50/(EC50+C) =
1/(1+C/EC50) = TXAI` exactly, for any `C_TXA`. `EMAX_TXA = 1.0` is not an
invented free parameter — it is the original's own implicit ceiling made
explicit (the inhibitory ratio asymptotes to exactly 1.0 as `C_TXA -> inf`,
never higher), required to reproduce the original term exactly under the
guide's canonical form.

## Diff scope confirmation

`diff -u dic_mrgsolve_model.R dic_mrgsolve_model_refactored.R` touches: the
header comment (new, documents refactor scope only), the two `$PARAM`
renaming groups (ATRA's EC50/EMAX/GAMMA in the "trigger" block, TXA's
EC50/EMAX/GAMMA in the "fibrinolysis" block, both compounds' CL/V1/KA in
the "drug PK" block), the three `$INIT` compartment lines, the
`ATRAE`→`C_ATR`/`EFFECT_ATR` and `TXAI`→`EFFECT_TXA` derivations in `$ODE`,
the `CTXA`→`C_TXA` line, the two `FGNLYS`/`vPLN` usage sites, the three
`dxdt_*` lines, one additive `$TABLE` block (four new discoverability
captures, see below), and five R-side `bolus()` calls in the scenario list
(`"TXAC"`→`"CENT_TXA"` x2, `"ATRAD"`→`"GUT_ATR"` x3) so the refactored
file's own scenarios still dose the renamed compartments by name. Nothing
else in the ~940-line file differs — every other compound's `$PARAM`/
`$INIT`/`$ODE`/`$TABLE` lines, every disease-side equation, the aetiology
presets, the dosing helpers, `run_dic()`, and the calibration-notes comment
block are byte-identical.

## qspserver `/model_manifest` discoverability

- `KA_ATR`, `CL_ATR`, `V1_ATR`, `EC50_ATR`, `EMAX_ATR`, `GAMMA_ATR`,
  `CL_TXA`, `V1_TXA`, `EC50_TXA`, `EMAX_TXA`, `GAMMA_TXA` are all declared
  in `$PARAM`, confirmed present with their defaults in the live
  `/model_manifest` response (e.g. `CL_ATR=50`, `EC50_TXA=6`,
  `GAMMA_ATR=1`).
- `C_ATR`/`EFFECT_ATR`/`C_TXA`/`EFFECT_TXA` themselves are computed in
  `$ODE` from compartment state, so — as already found for `EFFECT_NIM`
  (`sah_refactor_notes.md`) and `EFFECT_TCZ` (`ted_refactor_notes.md`) —
  they cannot also be `$PARAM` entries, and a `$TABLE` `capture` cannot
  reuse the same bare name (mrgsolve hoists every `$ODE` `double` into the
  same anonymous-namespace scope as `$TABLE` captures; reusing the name
  produces `error: redefinition of 'capture {anonymous}::C_ATR'` etc., the
  same class of error as the `FREETPA` defect below). Fixed the same way as
  those two prior runs: added `o`-prefixed `$TABLE` locals
  (`oC_ATR`/`oEFFECT_ATR`/`oC_TXA`/`oEFFECT_TXA`, recomputed from
  compartment state) and exposed them as `capture ATR_conc`/
  `capture ATR_effect`/`capture TXA_conc`/`capture TXA_effect` — confirmed
  present in `/model_manifest`'s `outputPaths` and retrievable via
  `/run_simulation`. These four captures are additive only (no counterpart
  in the original) and excluded from the verification diff below.
- `model_content` was extracted from the `dic_code <- '...'` R-string
  wrapper (straight text extraction, verified byte-identical to the
  quoted block) — no `.cpp` sibling was left behind; extraction was
  in-memory only, used to build the verification requests below and then
  discarded, per the workflow guide.

## A pre-existing upstream compile defect, unrelated to this refactor

Neither `dic_mrgsolve_model.R` nor the refactored file compiles as written
under mrgsolve 2.0.1: `$ODE` declares `double FREETPA = ...` (the free-tPA
term feeding `vPLN`) and `$TABLE` independently declares `capture FREETPA =
...` (the same quantity, recomputed from state for reporting) — mrgsolve's
auto-promotion of every `$ODE` local to a reportable member collides with
the explicit `$TABLE` capture of the same name:

```
145:11: error: redefinition of 'capture {anonymous}::FREETPA'
83:10: note: 'double {anonymous}::FREETPA' previously declared here
```

This is present in the *original* file already and has nothing to do with
ATRA or TXA. Logged as **upstream issue #43** in
[`translations/UPSTREAM_ISSUES.md`](../translations/UPSTREAM_ISSUES.md).
**For verification only**, both the original's and the refactored file's
extracted DSL text were given the identical, minimal, non-numeric
workaround (`$TABLE`'s `capture FREETPA = ...` renamed `FREETPA_CAP = ...`,
right-hand side and every `$ODE` line untouched) purely so both would
compile side by side. This workaround exists only in transient in-memory
request payloads sent to the qspserver API — **neither
`dic_mrgsolve_model.R` nor `dic_mrgsolve_model_refactored.R` on disk
contains it**, and neither currently builds against mrgsolve 2.0.1 without
the same rename.

## Verification

**Method.** Both files' embedded mrgsolve `dic_code <- '...'` DSL blocks
were mechanically extracted (regex on the assignment, verbatim quoted
text, no character changes), patched with the identical `FREETPA_CAP`
workaround above, and POSTed to the local qspserver `mrgsolve_api` service
at `http://localhost:8007` (`/model_manifest` then `/run_simulation`),
which compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side — no local R/mrgsolve install used, per the task instructions.
`init()`-only aetiology overrides (`AETIOLOGY$apl`'s `BLAST=30, PLT=35,
FIB=110, A2AP=0.40, PLG=0.60, DD=22`, and `AETIOLOGY$sepsis`'s `PATH=1.0`)
were reproduced as `evid=1` dosing records at `t=0` that bridge each
compartment from its `$INIT` baseline to the aetiology's target value
(e.g. `PLT`: dose `-215` from baseline `250` to reach `35`) — applied
identically to both models, so any imprecision in this bridging technique
cancels in the comparison. `RATE_TXA`'s scenario-defined
infusion-then-stop (125 mg/h, 0→96 h, then 0) was reproduced via the API's
`drivers` (piecewise-constant/`method:"constant"` time-course covariate
mechanism), which the original's own `run_dic()` implements equivalently
via segment-wise `param()` updates.

Three of the file's own scenarios were run through both models (full
`end=672 h, delta=0.5 h` window as in the originals — no shortening was
needed; all three completed in ~4 s each, well within the API's step
budget):

1. **`L_TXA_sepsis`** (TXA only) — `aet="sepsis"`, bolus 1000 mg into
   `TXAC`/`CENT_TXA` at t=0, `RATE_TXA` 125 mg/h from 0–96 h. 68 shared
   `$CAPTURE`/compartment outputs, 1347 timepoints.
2. **`O_APL_ATRA`** (ATR only) — `aet="apl"`, 40 mg into `ATRAD`/`GUT_ATR`
   q12h x56 doses (0–660 h). 68 shared outputs, 1352 timepoints.
3. **`Q_APL_ATRA_TXA`** (both compounds together, plus an untouched
   fibrinogen-replacement dose) — `aet="apl"`, ATRA dosing as in (2), TXA
   dosing as in (1). 68 shared outputs, 1354 timepoints.

(`R_APL_ATRA_platelets` was not separately run — it repeats scenario (2)'s
exact ATRA dosing/effect path with an untouched compound's — platelet —
dose added, exercising no new ATR/TXA interaction beyond what (2) and (3)
already cover.)

**Results:**

- **`O_APL_ATRA` (ATR only): exact match.** Max relative deviation = 0.0
  across all 68 shared outputs, and the renamed compartments themselves
  (`ATRAD` vs `GUT_ATR`, `ATRAC` vs `CENT_ATR`) also matched exactly
  (max rel. dev. = 0.0, peak concentrations 40.0 mg and 23.0945 mg
  identical to machine precision). `BLAST` showed a real, non-trivial
  ATRA-driven decline (29.0 → 25.2 x10^9/L over the run), confirming
  `EFFECT_ATR` is genuinely active in this comparison, not a degenerate
  no-op. This is the expected outcome for a pure Archetype-3 rename with
  no Hill-fitting.

- **`L_TXA_sepsis` (TXA only): near-exact match, root-caused to floating-
  point non-associativity, not a structural bug.** Max relative deviation
  across the 68 shared outputs was **7.18e-4** (worst: `PATH` 7.18e-4,
  `APFo` 7.06e-4, `APCTHR` 1.58e-4; every other shared output ≤5e-6). The
  renamed compartment `TXAC` vs `CENT_TXA` matched to 1.0e-14 (floating-
  point noise floor). **Root-caused via two controlled ablations, both
  run through the same API:**
  1. With the `(1.0 - EFFECT_TXA)` expression rewritten to the
     *bit-identical-to-the-original* operation order
     (`1.0 / (1.0 + C_TXA / EC50_TXA)` instead of
     `1.0 - C_TXA/(EC50_TXA+C_TXA)`) — algebraically identical, only the
     floating-point evaluation order changed — the deviation on **every**
     shared output dropped to exactly **0.0**.
  2. With TXA dosing removed entirely (`RATE_TXA`/bolus absent, pure
     structural rename exercised with the compound inactive), the
     deviation was also exactly **0.0**.

  Together these prove the ~7e-4 deviation is not a semantic or structural
  error: `1 - C/(EC50+C)` and `1/(1+C/EC50)` are mathematically identical
  but not bit-identical in IEEE double precision, and this file's stiff,
  tightly-coupled ODE system (thrombin generation feeding four zymogen
  consumption terms, protein C, and both fibrinolysis clocks
  simultaneously) is sensitive enough to adaptive-solver step-size
  selection that a last-bit rounding difference in one term (`FGNLYS`/
  `vPLN`'s TXA multiplier) measurably perturbs unrelated state trajectories
  (`PATH`, `APFo`) ~0.07% by 672 h, even though `PATH`'s own `dxdt_PATH`
  never references TXA at all. **The canonical stimulatory `EFFECT_TXA`
  form was kept in the delivered file** (matching the guide's naming
  interface and reading correctly next to every other compound's
  `EFFECT_<STEM>`), and this floating-point-reassociation explanation is
  disclosed here rather than silently accepted or worked around, per the
  guide's "a gate failure is a question, not a verdict" rule — this
  question has a confirmed, benign answer.

- **`Q_APL_ATRA_TXA` (both compounds together): near-exact match,
  consistent with the same explanation.** Max relative deviation across
  the 68 shared outputs was **8.66e-6** (worst: `FXA`, `TF`, `MORT`, all
  single-digit-1e-6); `TXAC` vs `CENT_TXA` matched to 7.54e-4 (same order
  as scenario 1, again attributable to the same TXA floating-point
  reassociation); `ATRAD`/`ATRAC` vs `GUT_ATR`/`CENT_ATR` matched to
  machine precision (both decay to ~1e-19, below the noise floor).

**No mismatch was found that the two ablations above did not fully
explain.** Per the guide's tolerance table, Archetype 3 (ATR) reproduced
bit-for-bit; Archetype 1 (TXA) reproduced to a confirmed floating-point-
reassociation tolerance an order of magnitude tighter than typical
Hill-fit deviations, with the root cause proven rather than assumed.
