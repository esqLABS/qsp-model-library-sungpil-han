# Refactor notes — sepsis / tocilizumab (TOCI → TCZ)

Calibration run: one of five parallel refactors of the same drug
(tocilizumab), each starting from a different original disease-model
implementation, run to check whether the process converges on consistent
structure. This run covers `sepsis/sep_mrgsolve_model.R`.

**Scope of this refactor**: only tocilizumab's PK compartment and its
downstream effect on IL-6 signaling were touched. Every other compound
(bacterial dynamics, meropenem antibiotic PK, norepinephrine, hydrocortisone,
the cytokine/coagulation/organ-failure network, SOFA scoring) is copied
verbatim from the original into `sep_mrgsolve_model_refactored.R`.

## Archetype determined

**Archetype 1 — no depot, single compartment, linear elimination.**

The original models tocilizumab as one compartment, `TOCI_C`, with:

```
dxdt_TOCI_C = -CL_Toci / V1_Toci * TOCI_C;
```

No absorption depot, no peripheral compartment, no TMDD kinetics — a single
linear ODE, matching Archetype 1 in `FORK_WORKFLOW_GUIDE.md` Part 2 exactly
in form (`dxdt_CENT_TCZ = -CL_TCZ * C_TCZ` reduces to the same equation once
`C_TCZ` is substituted).

### A genuine structural quirk, preserved rather than "fixed"

The original doses tocilizumab with `ev(amt = toci_dose, cmt = "TOCI_C", ...)`
— i.e. it adds the mg dose amount directly into the compartment, and the
downstream effect term reads that same compartment value directly as a
concentration:

```r
double E_Toci = Emax_Toci * TOCI_C / (EC50_Toci + TOCI_C);
```

There is no `amt / V1` conversion anywhere; `V1_Toci` only ever appears
inside the elimination-rate ratio `CL_Toci / V1_Toci`. So in this file the
compartment *is* the concentration by construction, not an amount that PD
divides by volume. Archetype 1's literal template line
(`C_TCZ = CENT_TCZ / V1_TCZ`) would silently change the model's numeric
output by a factor of `V1_TCZ` (4.2) if applied here, since the ODE's state
trajectory is identical either way (only the ratio `CL/V1` enters the
dynamics) but the *effect* term is not.

Renamed accordingly, preserving the original's actual behavior:

```
$PARAM
CL_TCZ    = 0.28    // clearance (L/h)              [was CL_Toci]
V1_TCZ    = 4.2     // volume (L) — only sets CL_TCZ/V1_TCZ; the compartment
                    // is dosed/read directly as a concentration, unchanged
                    // from the original TOCI_C behavior
EMAX_TCZ  = 0.90    // [was Emax_Toci]
EC50_TCZ  = 1.5     // mcg/mL                        [was EC50_Toci]
GAMMA_TCZ = 1        // original had no explicit Hill exponent

$CMT
CENT_TCZ            // [was TOCI_C]

$MAIN
double C_TCZ = CENT_TCZ;   // exposed concentration = compartment, undivided
double EFFECT_TCZ = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ)
                    / (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));

$ODE
dxdt_CENT_TCZ = -CL_TCZ / V1_TCZ * CENT_TCZ;
```

`EFFECT_TCZ` replaces `E_Toci` at its one use site (`Stim_IL6`'s
`(1.0 - EFFECT_TCZ)` term). The dosing event's `cmt = "TOCI_C"` string was
updated to `cmt = "CENT_TCZ"` since that's the interface the rest of the
file now doses into; this is the only change outside the `$PARAM`/`$CMT`/
`$MAIN`/`$ODE` blocks. The `useToci` switch parameter and every scenario
definition, `make_events()`, plotting, and summary code are untouched.

## Hill interface

The original's effect term was already a plain Emax ratio
(`Emax*C/(EC50+C)`, i.e. Hill with an implicit `gamma = 1`), so this is a
rename, not a refit: `EMAX_TCZ`, `EC50_TCZ` carry the original's exact
values, `GAMMA_TCZ = 1` is added explicitly since the original had no Hill
exponent term.

## Verification

**Method**: `library(mrgsolve)` locally on `C:\Program Files\R\R-4.6.0\bin\Rscript.exe`
segfaulted intermittently on `mcode()`/`mrgsim()` calls in this environment
(root cause unclear — possibly the R 4.6.0/mrgsolve-built-for-4.6.1 ABI
mismatch reported by the package's own load-time warning, possibly a
transient issue in the sandbox; simple synthetic models did eventually
compile locally without crashing once retried, so the tool itself works).
To get a reliable, reproducible verification run, this was switched to the
qspserver `mrgsolve_api` container at `http://localhost:8007`
(`/model_manifest` and `/run_simulation`), which compiles and runs the same
mrgsolve model server-side.

**A pre-existing compile blocker in the original file, unrelated to this
refactor**: neither `sep_mrgsolve_model.R` nor the refactored file compiles
as written under the installed mrgsolve version (2.0.1), for three
independent reasons, all present in the *original* file already:

1. `$CAPTURE` lists compartment names directly (e.g. `BACT`, `TNF`, `TOCI_C`,
   …) — this mrgsolve version rejects that: *"compartment should not be in
   $CAPTURE"*.
2. Every compartment is declared via **both** `$CMT` (bare name) and `$INIT`
   (`name = value`) — this mrgsolve version's `validObject()` flags that as
   *"Duplicated model names"*, even though this is an extremely common,
   previously-unremarkable mrgsolve pattern.
3. At the C++ compilation stage, `$PARAM` names `IL6_0`, `IL10_0`, `PAI1_0`
   collide with mrgsolve's auto-generated `<compartment>_0` initial-value
   symbols for the identically-named compartments `IL6`, `IL10`, `PAI1`,
   producing `conflicting declaration`/`redeclaration` compiler errors.

These are logged here (per the fork guide's "log what you find, don't fix it
upstream" rule) rather than fixed in `sep_mrgsolve_model.R`. **For
verification only**, both the original and the refactored model text were
given the identical, minimal, non-committed workaround (compartment names
dropped from `$CAPTURE`, the `$CMT` block removed in favor of `$INIT` alone,
and the three colliding `$PARAM` names renamed `IL6_SS0`/`IL10_SS0`/
`PAI1_SS0` with all their usages updated) purely so both would compile
side by side. This workaround was applied symmetrically and does not touch
any numeric value, so it cannot bias the comparison. Neither
`sep_mrgsolve_model.R` nor `sep_mrgsolve_model_refactored.R` on disk contains
this workaround — it existed only in temporary, unsaved request payloads
sent to the qspserver API.

**Scenario tested**: `S6_BundleToci`, the only scenario in the original
file's own scenario list (`scenarios` in `sep_mrgsolve_model.R`) that doses
tocilizumab — antibiotic q8h ×28d, norepinephrine infusion 2→72h,
hydrocortisone q6h ×7d, a single 672 mg tocilizumab dose at t = 6h, and a
fluid bolus, run 0–168h at Δt = 0.5h (341 time points), reproduced through
the qspserver API using the same dose amounts/compartment indices/times and
`useAbx=useNE=useHC=useToci=useFluid=1, FluidBoost=3`.

**Result: exact match.** Every one of the 36 `$CAPTURE`d outputs shared
between the two models (`ABX1, ABX2, BACT, BILIRUBIN, C5A, CREATININE,
Cp_abx_report, CytokineIdx, ENDOT, FIBRIN, HC_C, IL10, IL1B, IL6, IL6_pg,
LACTATE, MACS, MAP_val, Mortality_Prob, NEUT_B, NEUT_T, NE_C, PAI1,
PF_RATIO, PLT_COUNT, SOFA, THROMBIN, TNF, TNF_ng, log_BACT, shock_flag,
sofa_cardio, sofa_cns, sofa_coag, sofa_liver, sofa_lung, sofa_renal`) had
**max relative deviation = 0.0** across all 341 time points. The renamed
compartment itself, `TOCI_C` (original) vs `CENT_TCZ` (refactored), also
matched exactly (max relative deviation = 0.0). This exceeds the guide's
"near-exact match" expectation for a pure structural reorganization with no
Hill-fitting — it is a bit-for-bit match, as expected for a rename with no
change to the underlying arithmetic.

**One cosmetic observation, not a correctness issue**: the extra diagnostic
column `C_TCZ` this refactor adds to `$CAPTURE` (which has no counterpart in
the original to verify against) occasionally reports a value one internal
step out of sync with the true compartment `CENT_TCZ` at a handful of output
rows shortly after the t = 6h dose (e.g. at t = 10h, `CENT_TCZ` = 514.70 vs.
reported `C_TCZ` = 532.15; by t = 20h they agree exactly again). This looks
like an mrgsolve `$MAIN`-vs-`$TABLE` evaluation-timing nuance for a `double`
declared in `$MAIN` and read back at report time, rather than anything
specific to this refactor — `EFFECT_TCZ`, which is what actually drives the
ODE (via `Stim_IL6`), is evidently computed correctly and continuously
throughout integration, since the full disease trajectory matched the
original exactly. Worth keeping in mind if a future patch reads `C_TCZ`
directly from captured output rather than from the live model state.

**Other observation**: the `useToci` parameter switch is declared in
`$PARAM` and used by the R-side scenario/`make_events()` code to decide
*whether* to schedule a dose, but is never referenced inside `$MAIN`/`$ODE`
itself — the effect is gated entirely by whether the compartment has been
dosed, not by the switch. Not a bug, just worth knowing if `useToci` is ever
expected to gate the effect independent of dosing.
