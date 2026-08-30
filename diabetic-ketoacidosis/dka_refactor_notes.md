# Refactor notes — `dka_mrgsolve_model.R` (insulin pharmacokinetics only)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2) and the
existing row in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
("Insulin pharmacokinetics", classified "Redirect concentration (clean
single site)"), only insulin's PK/effect block was rewritten. Insulin is the
**only** modeled compound in this file (this is a 42-ODE electrolyte/acid-base/
ketone model with a single drug, not a multi-drug file), so there is nothing
else to leave untouched in the pharmacological sense — but every other one of
the model's 42 compartments, its acid-base solve, ketone/lactate/renal/brain
blocks, and its entire R scenario-runner layer were left byte-for-byte
identical except for the handful of lines that reference the renamed insulin
identifiers by name (listed below). Confirmed by `diff`, see "Diff scope
confirmation" below.

## Archetype determination

**Bespoke — does not cleanly fit any of the guide's four archetypes.**

The insulin PK/PD block has four components, three of which are individually
close to a standard shape but which combine into something the guide's
worked examples don't cover:

1. **A subcutaneous absorption depot into a central pool**, structurally
   close to Archetype 3 (depot + central, no peripheral compartment) —
   *except* the central "compartment" (`INSP`, now `CENT_INS`) is written
   directly in **concentration units** (uU/mL) rather than as an amount
   divided by a volume:
   ```c
   dxdt_CENT_INS = ((Riv + abs_sc) * 1.0e6 - CL_INS * 1000.0 * CENT_INS) / (V1_INS * 1000.0);
   ```
   (`Riv` = IV infusion + endogenous secretion in U/h, `abs_sc` = absorption
   flux in U/h; the `1.0e6`/`1000.0` factors convert U/h and L into uU and
   mL so the state itself comes out in uU/mL.) This is algebraically
   equivalent to an amount-over-volume Archetype-3 central compartment, just
   pre-divided by volume at the point of writing the ODE — the census
   classifier's "Redirect concentration (clean single site)" label fits this
   exactly: there already is one, single, clean concentration variable, it
   just isn't produced by a `double C = A/V;` line the way Archetype 1-3's
   worked examples show it.
2. **Two downstream effect-site compartments, not one.** Peripheral actions
   (lipolysis suppression, glucose disposal stimulation, the potassium
   shift) read from a first-order effect site with tau 0.35 h (`INSEF`, now
   `EFF_PERI_INS`); hepatic glucose output suppression reads from a
   *separate* first-order effect site with tau 1.5 h (`INSES`, now
   `EFF_HEP_INS`). The file's own header (point 2 of the top-of-file
   commentary) states this separation, not any fitted rate, is what
   produces the model's characteristic treated-glucose-curve shape — this
   is exactly the guide's stated exception ("two [concentration variables]
   only when a genuinely different tissue site matters"), so both effect
   sites were kept, renamed but structurally untouched.
3. **Four disease-facing Hill/Emax terms, not one.** Insulin's four
   physiological actions modeled here (lipolysis suppression, hepatic
   glucose output suppression, glucose disposal stimulation, Na/K-ATPase
   potassium shift) are four independently-parameterized Emax/EC50 ratios,
   not one shared `EFFECT_<STEM>`. The file's own calibration notes at the
   bottom of the file call these "the four half-maximal constants
   [that] are the load-bearing parameters of the whole model" — collapsing
   them into a single `EFFECT_INS` would both misrepresent the model and
   contradict the guide's instruction not to flatten a mechanistically rich
   model for convenience. Given no naming-convention slot exists for "one
   compound, four disease axes", this was handled as a bespoke, documented
   extension: `EFFECT_INS_LIP`, `EFFECT_INS_HGP`, `EFFECT_INS_UP`,
   `EFFECT_INS_KSH`, each built from the shared `EMAX_INS_*`/`EC50_INS_*`/
   `GAMMA_INS_*` naming pattern, keyed by physiological axis instead of by
   compound alone.
4. **An endogenous-secretion / portal-privilege construction with no analog
   in the guide's archetypes.** Insulin enters `CENT_INS` from two sources
   fused into a single flux, `Riv = INS_IV + sec` (exogenous IV infusion
   plus glucose-driven endogenous secretion `sec`), and a separate,
   secretion-privileged "portal" quantity (`IeffP`) feeds the ketogenesis
   (CPT-1 gate) block only — modeling the fact that secreted insulin passes
   the liver first at high concentration, while exogenous insulin does not.
   This asymmetry (file header, point 3) is the mechanism that lets
   residual beta-cell function suppress ketosis while permitting severe
   hyperglycaemia (i.e. it is literally what distinguishes HHS from DKA in
   this model) and was preserved exactly — renamed (`INSEF`→`EFF_PERI_INS`
   inside `IeffP`'s definition) but not restructured.

Per the guide's "None of these fit" fallback: renamed to the convention,
kept isolated from the surrounding disease code (it already was — insulin's
PK/PD lives entirely in `$ODE` section 2, one contiguous block, plus the one
`IeffP` line inside section 4's ketogenesis block that reads insulin's
already-computed effect-site value), and documented here rather than forced
into Archetype 1/2/3's `double C_STEM = CMT/V;` shape or Archetype 4's TMDD
receptor-binding shape (there is no receptor-binding state here — the two
"effect sites" are first-order concentration-lag compartments, not a
bound-receptor pool).

## Naming applied

| Role | Original | Refactored |
|---|---|---|
| Absorption depot compartment | `INSSC` | `GUT_INS` |
| Central "compartment" (already stored as a concentration) | `INSP` | `CENT_INS` |
| Peripheral effect site (tau 0.35h) | `INSEF` | `EFF_PERI_INS` |
| Hepatic effect site (tau 1.5h) | `INSES` | `EFF_HEP_INS` |
| Central volume | `VINS` | `V1_INS` |
| Clearance | `CLINS` | `CL_INS` |
| Absorption rate (rapid-acting analogue, the only one actually used) | `KA_SC` | `KA_INS` |
| Absorption rate (regular insulin) | `KA_REG` | `KA_REG` (unchanged — see "dead parameter" below) |
| Peripheral effect-site time constant | `TINSE_F` | `TAU_PERI_INS` |
| Hepatic effect-site time constant | `TINSE_S` | `TAU_HEP_INS` |
| Exposed plasma concentration | *(none — `INSP` itself was read directly)* | `C_INS` (new `$ODE`-local double, `= CENT_INS`) |
| IR-scaled peripheral effect-site value | `Ieff` | `X_INS_PERI` |
| IR-scaled hepatic effect-site value | `IeffS` | `X_INS_HEP` |
| Lipolysis suppression | `EMAX_LIP`, `IC50_LIP` | `EMAX_INS_LIP`, `EC50_INS_LIP` (+ new `GAMMA_INS_LIP = 1.0`) |
| Hepatic glucose output suppression | `EMAX_HGP`, `IC50_HGP` | `EMAX_INS_HGP`, `EC50_INS_HGP` (+ new `GAMMA_INS_HGP = 1.0`) |
| Glucose disposal stimulation | *(none)*, `EC50_UP` | new `EMAX_INS_UP = 1.0`, `EC50_INS_UP` (+ new `GAMMA_INS_UP = 1.0`) |
| Na/K-ATPase potassium shift | *(none)*, `IC50_KSH` | new `EMAX_INS_KSH = 1.0`, `EC50_INS_KSH` (+ new `GAMMA_INS_KSH = 1.0`) |
| Disease-facing effect terms | `fLIP`, `fHGP`, `fUP`, `fKSH` (kept, now derived from the `EFFECT_INS_*` terms below) | `EFFECT_INS_LIP`, `EFFECT_INS_HGP`, `EFFECT_INS_UP`, `EFFECT_INS_KSH` |
| Portal privilege quantity (bespoke, ketogenesis-only) | `IeffP` | `IeffP` (name unchanged — internal quantity outside the naming convention's scope, same treatment the KTX calibration run gave `TCZIN`/`TCZON`; only its `INSEF` reference was updated to `EFF_PERI_INS`) |
| Basal plasma insulin, endogenous secretion params, portal factor | `INSB`, `SECMAX`, `KGSEC`, `PORTF` | unchanged (physiological constants outside the PK-naming table's scope) |

**The four Hill terms are a rename, not a refit.** The original's `fLIP`/
`fHGP` are already `1 - Emax*C/(C+IC50)` and `fUP`/`fKSH` are already
`C/(C+EC50)` — the Emax=1 special case of the same form. So:
```c
double EFFECT_INS_LIP = EMAX_INS_LIP * pow(X_INS_PERI, GAMMA_INS_LIP)
               / (pow(EC50_INS_LIP, GAMMA_INS_LIP) + pow(X_INS_PERI, GAMMA_INS_LIP));
double fLIP  = 1.0 - EFFECT_INS_LIP;
...
double EFFECT_INS_UP  = EMAX_INS_UP * pow(X_INS_PERI, GAMMA_INS_UP)
               / (pow(EC50_INS_UP, GAMMA_INS_UP) + pow(X_INS_PERI, GAMMA_INS_UP));
double fUP   = EFFECT_INS_UP;
```
with `GAMMA_INS_* = 1.0` throughout (the original had no explicit Hill
coefficient — same precedent the guide's own Archetype-4 worked example
uses for `GAMMA_TCZ`) and `EMAX_INS_UP = EMAX_INS_KSH = 1.0` (the original's
`C/(C+EC50)` form has an *implicit* Emax of 1; this makes it an explicit,
named parameter rather than a hidden constant — parallel to the `GAMMA=1`
precedent, not a new pharmacological claim). With these values, every
`EFFECT_INS_*` reduces algebraically to exactly the original expression for
any input, confirmed by the exact-match verification below.

**A genuinely dead parameter was found and left untouched.** `KA_REG`
("s.c. absorption, regular human insulin (1/h)") is declared in `$PARAM` but
never referenced anywhere in `$ODE` — the single `dxdt_INSSC`/`abs_sc` line
always uses `KA_SC` regardless of which insulin formulation a scenario might
conceptually intend. This is not a compile-blocking defect (an unused
parameter doesn't break the build) so it isn't in
`translations/UPSTREAM_ISSUES.md`, but it's worth flagging for whoever next
touches this file: either the model was meant to let scenarios choose an
absorption rate and that switch was never wired up, or `KA_REG` is
vestigial. Left as-is per the guide ("never invent or default a PK
parameter, and don't add one the original didn't have" — the converse also
applies: don't delete one either).

## Renaming outside the `$ODE`/`$PARAM`/`$CMT` block

Unlike the KTX calibration run (where no wrapper-level R code referenced the
tocilizumab compartment by name), this file's R scenario-runner layer *does*
set insulin's compartments by literal name in a few places, so the whole
sibling file — not just the quoted DSL string — had to stay internally
consistent for the file to still run standalone:

- `SIZE_PARAMS` (the body-weight-scaling vector used by `patient()`):
  `"VINS", "CLINS"` → `"V1_INS", "CL_INS"`.
- `run_protocol()`'s bolus line: `st$INSP <- st$INSP + bolus * 1e6 / (p$VINS * 1000)`
  → `st$CENT_INS <- st$CENT_INS + bolus * 1e6 / (p$V1_INS * 1000)`.
- Three `st$INSSC <- st$INSSC + ...` lines (in `sc8_transition()` and
  `sc9_subcutaneous()`) → `st$GUT_INS <- st$GUT_INS + ...`.
- The counter-regulatory-hormone block's glucagon set point
  (`double gcg_t = GCG0 * (1.0 + 2.2 * (1.0 - Ieff / (Ieff + 40.0)) ...`,
  section 3 of `$ODE`) reads insulin's own peripheral effect-site value
  directly — this is a consumer of insulin's exposed state, not another
  compound's code, so it was updated to `X_INS_PERI` along with the rest of
  the insulin block; it is the one place outside insulin's own two
  contiguous blocks (section 2, and the `IeffP` line in section 4) that
  needed a change.

## qspserver compatibility

- `C_INS` and each `EFFECT_INS_*` are `$ODE`-local `double`s, not `$PARAM`
  entries. This was tried the other way first (declaring them in `$PARAM`
  with a `= 0.0` default, per the guide's compatibility item 2) and it does
  not compile: mrgsolve 2.0.1 passes every `$PARAM` value into `$ODE` as a
  `const double&`, so assigning to it (`C_INS = CENT_INS;`) is a hard
  compile error ("assignment of read-only reference"). A repeat check of
  all five prior calibration runs' `_refactored.R` files
  (kidney-transplant-rejection, sepsis, thyroid-eye-disease,
  polymyalgia-rheumatica, rheumatoid-arthritis) confirms none of them put
  `C_TCZ`/`EFFECT_TCZ` in `$PARAM` either, for the same structural reason —
  this appears to be a hard mrgsolve constraint the guide's wording doesn't
  fully anticipate, not something this file did wrong.
- Discoverability is instead satisfied via `$CAPTURE`/`$TABLE`: `C_INS` and
  all four `EFFECT_INS_*` values are captured every step, under the names
  `C_INS_report`, `EFFECT_INS_LIP_report`, `EFFECT_INS_HGP_report`,
  `EFFECT_INS_UP_report`, `EFFECT_INS_KSH_report` — **not** the bare
  `C_INS`/`EFFECT_INS_*` names, deliberately. This file's `$TABLE` idiom is
  `capture NAME = EXPR;` (there is no separate `$CAPTURE @annotated` block
  in this model at all), and that macro declares a *new* `NAME` in the
  fused `$ODE`/`$TABLE` scope — `capture C_INS = C_INS;` would self-shadow
  the already-computed `$ODE` value and read uninitialised memory instead
  of it. Verified this reasoning is necessary, not merely cautious, by
  checking the C++ scoping rule directly (a variable's own name is in scope
  from its declarator onward, so its initializer expression already refers
  to itself). The `_report`-suffixed names are confirmed present in
  `/model_manifest`'s `outputPaths` (see verification below) and are listed
  in the updated census row.
- All 42 `$CMT` names (with the four insulin ones renamed) and every
  `$PARAM` (with the six new Hill parameters) are confirmed present and
  correctly defaulted via `/model_manifest` — see below.

## Diff scope confirmation

`diff -u dka_mrgsolve_model.R dka_mrgsolve_model_refactored.R` is 217 lines
(mostly additions — the new `GAMMA_INS_*`/`EMAX_INS_UP`/`EMAX_INS_KSH`
parameter lines and their explanatory comments), touching exactly: the
insulin `$PARAM` block, the four insulin `$CMT` lines, the four insulin
`$MAIN` init lines, `$ODE` section 2 in full plus the one `gcg_t` line in
section 3 and the one `IeffP` line in section 4, the three insulin
`$TABLE` capture lines (plus five new ones), `SIZE_PARAMS`, the bolus line
in `run_protocol()`, and the three `st$INSSC` lines in `sc8_transition()`/
`sc9_subcutaneous()`. Nothing else in the file differs: all 38 non-insulin
compartments, the entire acid-base/ketone/lactate/renal/water/brain/
respiratory `$ODE` blocks, every other `$TABLE` capture, `$MAIN`'s other 38
init lines, the fluid library, `patient()`'s non-insulin behavior, and every
scenario function's non-insulin logic are byte-for-byte identical.

## Verification

**Method.** Both the original's and the refactored sibling's embedded
mrgsolve DSL blocks were mechanically extracted (find the `code <- '...'`
assignment, pull the quoted text verbatim) and run through the qspserver
`mrgsolve_api` container (`http://localhost:8007`, confirmed healthy),
using `POST /model_manifest` to confirm both compile and expose the right
parameters/outputs, and `POST /run_simulation` to compare every shared
`$CAPTURE`d output across eleven scenarios built directly from dosing values
and parameter overrides already defined in the original file's own R code
(not invented ones):

1. **`sc0_untreated`** — `leadin(patient(), 48)`: default patient, 48 h,
   delta 0.25.
2. **`sc10_phenotypes`, all 8 cases** — each is a single `leadin()` call in
   the original file (`patient()` with the exact override list each
   phenotype uses: severe sepsis `ILL0=0.95`; child 30 kg, reproducing
   `patient()`'s own body-size scaling of `SECMAX`/`V1_INS` (`VINS`)/
   `CL_INS` (`CLINS`)/etc.; euglycaemic DKA `BETA=0.10, SGLT2=1,
   GLYCO0=60`; alcoholic ketoacidosis; HHS early/advanced; established AKI
   `GFRMAX=3.0`), each for that scenario's own duration (24-54 h).
3. **Constant IV insulin infusion** — `INS_IV = 7.0` (the file's own
   `ins = 0.10` U/kg/h dose at BW = 70, used by `sc1_standard`/
   `sc3_insulin_dose`'s low-dose arm) plus `BETA = 1.0` to also engage
   endogenous secretion, run continuously from t=0 for 24 h. This is the
   model's own native continuous-infusion dosing mechanism (`INS_IV` is a
   plain `$PARAM`, entered with no `ev()`/`data_set()` anywhere in the
   file — matching the guide's "continuous infusion input, not an event"
   edge case) driven for its full duration rather than through
   `run_protocol()`'s taper schedule (see "What could not be verified"
   below for why).
4. **SC bolus insulin dose** — 21 U (`0.3 * BW` at BW = 70, the exact
   opening-dose size `sc9_subcutaneous()` uses) via an NM-TRAN-style bolus
   dose into `$CMT` #19 (`GUT_INS`/`INSSC` in both models — the compartment
   order was not changed by the rename), with `ILL0 = 0.30` matching that
   scenario's own patient override, 20 h.

Every shared `$CAPTURE`d column (48 of the original's 48 non-compartment
outputs; the refactored sibling additionally captures the five new
`*_report` discoverability columns, which have no original counterpart to
diff against) was compared across the full output time grid for each
scenario.

**Result: exact match. Maximum absolute and relative deviation observed =
0.0 (bit-identical) across all 11 scenarios and every shared `$CAPTURE`d
column.** This is the expected outcome for a pure structural
reorganization/rename with no Hill-refitting, per the guide's tolerance
table.

**A pre-existing, unrelated build blocker was found and worked around only
for this verification, not fixed in either file** — logged as upstream
issue #37 in
[`translations/UPSTREAM_ISSUES.md`](../translations/UPSTREAM_ISSUES.md):
the pH-bisection block in `$ODE` declares `double lo = 5.60, hi = 8.30, pHm,
hres;` on one line, and mrgsolve 2.0.1's preprocessing drops the `double`
type from every name after the first, producing "was not declared in this
scope" for `hi`/`pHm`/`hres`. This has nothing to do with insulin (it's in
the acid-base block) and reproduces identically in the untouched original.
Verification here used in-memory-only copies of both models' DSL text with
that one line mechanically split into four separate `double` declarations
(same values, same order, no numeric change), applied identically to both
sides purely to make them buildable for comparison — **neither the tracked
original nor the delivered `_refactored.R` was changed**; both still
contain the single multi-declarator line exactly as written.

**A qspserver infrastructure flakiness was observed and is unrelated to
either model.** Two of the ~25 `/run_simulation` calls made during this
verification returned a transient `500` ("worker was killed by signal 6" /
"double free or corruption") on the *first* attempt and then succeeded with
a bit-identical result on an immediate retry with the identical payload.
This reproduced for both the original and the refactored DSL at different
points, under no obvious pattern (large-output S1 once, and once more on a
rerun of the same request), and is almost certainly a concurrency/memory
issue in the `mrgsolve_api` container itself (unrelated to model content,
since retrying the *exact same request* immediately succeeded and matched
every other passing run bit-for-bit) rather than a defect in either DSL.
Not logged in `UPSTREAM_ISSUES.md` since that file is for defects in the
disease-model corpus, not the shared qspserver service — flagged here so a
future verification run isn't surprised by an occasional transient 500.

**What could not be verified: `run_protocol()`'s live closed-loop feedback
itself.** `sc1_standard`/`sc2_levers`/`sc3_insulin_dose`/etc. all run
through `run_protocol()`, an R-side control loop that reads live glucose
and potassium from a short probe simulation every 0.25 h and adjusts
`RATE_FL`/`INS_IV`/`KCL`/fluid choice accordingly (dextrose-switch,
insulin-taper, and the ADA potassium-hold rule are all live feedback
branches, not a fixed time-based schedule). Reproducing this exactly via
the qspserver API would require chaining many short `/run_simulation` calls
with the full 42-compartment state carried forward between them. The API's
`SimRequest` schema was checked directly: its `events` field (state-variable
replace/add/multiply) is documented as deSolve-engine-only ("other engines
ignore the field") and confirmed inert for mrgsolve `model_content` by a
direct test (a `replace` event on `GLU`/`INSP` at t=0 had no effect on the
t=0 output row); mrgsolve models only honor the NM-TRAN-style `dosing` field
(additive amount into a numbered compartment) and constant `parameters`
overrides for the whole run — there is no way to inject an arbitrary
full-state initial condition or a piecewise-time-varying parameter through
this API for an mrgsolve model. Given that constraint, scenarios 3 and 4
above were built to exercise the same PK/PD block through the model's own
supported continuous-infusion and bolus-dosing mechanisms (using dose values
taken directly from the file's own scenarios) rather than attempting a
partial, easy-to-misjudge reimplementation of `run_protocol()`'s feedback
logic. Since the underlying ODE mathematics is confirmed bit-identical
across a wide range of static parameter/dosing conditions (natural history,
eight distinct phenotypes spanning the BETA/ILL0/SGLT2/ALCOHOL/GFRMAX/BW
space, continuous IV infusion, and SC bolus), and `run_protocol()`'s
decision logic itself is untouched R code that consumes only `GLUmM`/`KmM`
(both confirmed bit-identical outputs), there is no mechanism by which the
closed-loop scenarios could diverge between the two models — but this is an
inference from the verified building blocks, not a directly-executed check
of `run_protocol()` end-to-end, and is disclosed as such rather than
silently assumed.

## Anything else worth flagging

- This is the first calibration-style refactor in this fork where the
  refactored compound is the file's *only* modeled drug (all five prior
  worked examples were multi-compound files where most of the diff-scope
  discipline was about leaving 5-20 *other* compounds untouched). The
  "leave everything else alone" discipline here instead applied to the
  other 38 non-insulin compartments and the disease physiology built on
  top of them.
- The four-Hill-term, two-effect-site structure means this file's insulin
  PK/PD is the most mechanistically elaborate "single compound" case
  encountered in this fork's naming-convention work so far, despite having
  no TMDD/receptor-binding compartment at all — a reminder that "bespoke"
  and "Archetype 4" are not the same axis (mechanistic richness vs.
  receptor-binding-as-a-state are independent).
