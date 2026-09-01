# Refactor notes — `long-covid/pasc_mrgsolve_model.R`

## Scope of this pass — the census mislabeled two of this file's four rows, the other two were correct

The census (`driver-patches/data/compound_perturbation_census.md`) had
four rows for this file, all classified "Redirect concentration (clean
single site)":

| Compound (as listed) | Correct? |
|---|---|
| `MET` | Correct as listed — metformin |
| `SERT` | Correct as listed — sertraline |
| `Neuroinflammation / BBB (LDN)` | Mislabeled — LDN is low-dose naltrexone, not a process description |
| `Viral Kinetics (NIRM)` | Mislabeled — NIRM is nirmatrelvir, not a process description |

This is the same pattern already documented in
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md` and
`prostate-cancer/pc_refactor_notes.md`: the classifier's stem detection
was right (`LDN`, `NIRM` are genuinely this file's own PK-parameter
stems — `F_LDN`/`ka_LDN`/`CL_LDN`/`Vd_LDN`/`EC50_LDN`, and
`F_nirm`/`ka_nirm`/`CL_nirm`/`Vd_nirm`/`IC50_nirm`), but its display-name
derivation picked up the nearest `$PARAM` section comment
(`// ---- Neuroinflammation / BBB ----` for LDN, `// ---- Viral Kinetics
----` for NIRM) instead of the compound's own identity. Unlike
`prostate-cancer` (which also had 7 entirely *missing* compounds), this
file's census was complete at 4 rows — **no additional real compounds
were found** beyond these four. `use_anticoag` was checked and is *not*
a fifth compound: it is a bare 0/1 switch multiplying a fixed rate
constant (`kAnticoag`) directly in the `Fibrin` ODE, with no depot, no
central compartment, no concentration variable, and no dose amount
anywhere in the file — there is nothing to redirect (see "Anticoagulant
is out of scope" below).

## Real compound identities (per the task's questions)

- **NIRM = nirmatrelvir** (with ritonavir boosting implicit in the
  parameter comment "total clearance (boosted by RTV)"), oral antiviral,
  300 mg BID dosed for an extended 15-day course (Scenario 2; standard
  Paxlovid regimen is 5 days, this model's own comment/scenario explicitly
  tests the extended-course hypothesis referenced by the STOP-PASC trial
  cited in the file's own header comment).
- **MET = metformin**, oral antidiabetic/AMPK activator, 500 mg BID
  maintenance (Scenario 3), repurposed here per the COVID-OUT RCT
  rationale cited in the file's own header.
- **SERT = sertraline**, oral SSRI, 50 mg QD (Scenario 5).
- **LDN = low-dose naltrexone**, oral opioid-receptor/TLR4 antagonist,
  4.5 mg QD (Scenario 4) — the low-dose regimen used off-label for
  post-viral neuroinflammatory conditions, not naltrexone's on-label
  50 mg dose.

All four have their own gut depot + central concentration compartment
pair, and each is dosed independently via its own `use_<drug>` switch and
its own `ev()` in the R script (Scenarios 2-5), combined in Scenarios 6-7.

## Anticoagulant is out of scope (not a missed 5th compound)

`use_anticoag` (0/1) appears only in
`dxdt_Fibrin = ... - use_anticoag * kAnticoag * Fibrin`. There is no
`A_anticoag`/`C_anticoag`-style compartment, no PK parameter block for
it, and no dose amount is ever set for it anywhere in the R script (no
scenario turns it on: all seven `simulate_scenario()` calls pass `0` for
the `use_ac_` argument). It is a bare mechanistic on/off switch on a
fixed rate constant, not a dosed compound with its own exposure — nothing
in this file resembles a concentration to redirect. Left entirely
unchanged; not added to the census.

## Archetype per compound

All four are **Archetype 3 minus peripheral** (depot + central, linear,
no peripheral compartment) — a direct match, no forcing needed:

- `dxdt_GUT_<STEM> = -KA_<STEM> * GUT_<STEM>;`
- `dxdt_CENT_<STEM> = use_<drug> * (KA_<STEM>*GUT_<STEM>*F_<STEM>/V1_<STEM>) - KE_<STEM>*CENT_<STEM>;`

One structural note distinguishing this file from the guide's own
template: **the original's own central-compartment state is already a
concentration, not an amount.** The original's `C_nirm`/`C_met`/
`C_sert`/`C_LDN` are declared directly as `$CMT` states in "ug/L = ng/mL"
units and integrated as concentration-state ODEs (`dxdt_C_nirm =
use_nirm*(ka_nirm*A_nirm*F_nirm/Vd_nirm) - ke_nirm*C_nirm`), rather than
an amount-state divided by volume in a separate expression as the guide's
worked template shows. The refactor preserves this exactly: `CENT_<STEM>`
is the renamed concentration-state compartment (not an amount), and
`C_<STEM>` is a trivial identity of it (`C_NIRM = CENT_NIRM`), not a
`/V1_<STEM>` division — dividing a value that is already a concentration
by volume again would double-count the volume and change the numbers.
This is the same situation `prostate-cancer/pc_refactor_notes.md`
documents for Docetaxel's `Doc_c`/`Doc_p` (a bespoke case there, since
Docetaxel's own structure didn't otherwise fit an archetype); here all
four compounds share this same concentration-state convention and each
still cleanly fits Archetype 3 minus peripheral, so no bespoke handling
was needed.

`KE_<STEM>` (renamed from the original's `ke_nirm`/`ke_met`/`ke_sert`/
`ke_LDN`, computed in `$MAIN`) reuses the original's exact formula,
`(CL_<STEM>*1000)/(24.0*V1_<STEM>)`, unchanged. **This formula's own
units are worth flagging as a likely pre-existing modelling quirk, not
fixed here**: `CL_<STEM>` is declared in L/h and `V1_<STEM>` in L, so the
standard L/h -> /day conversion is `CL*24/V1`; the original instead
computes `CL*1000/(24*V1)`, i.e. divides by 24 rather than multiplying,
alongside an unexplained factor of 1000 with no stated unit rationale
(the four `CL_<STEM>`/`V1_<STEM>` value pairs don't correspond to a
mg-vs-ug or L-vs-mL correction that would explain it either). This
"defect" is not build-blocking and not touched — per "parameter values
always come from the original file" and "log what you find, don't fix it
upstream," the formula is carried forward exactly, renamed only. Not
logged as a numbered `UPSTREAM_ISSUES.md` entry because it doesn't block
compilation and isn't confirmed wrong (a nonstandard unit convention is
possible, just undocumented) — flagged here for a reviewer's awareness.

## The Hill interface

All four compounds' original effect terms were already plain
Emax-implicit Hill ratios (`C/(EC50+C)`, no explicit Emax multiplier, no
explicit exponent) — every `EFFECT_<STEM>` below is a **rename**, not a
fit, with `EMAX_<STEM> = 1.0` and `GAMMA_<STEM> = 1.0` added as new named
parameters per the guide ("pull out EMAX/EC50/GAMMA... GAMMA=1 if the
original had no explicit Hill coefficient"):

- **NIRM**: `EFFECT_NIRM = EMAX_NIRM*C_NIRM/(EC50_NIRM+C_NIRM)` — was
  `nirm_eff = C_nirm/(IC50_nirm+C_nirm)`. `EC50_NIRM` renamed from
  `IC50_nirm` (same value, 0.003 ug/mL); the original's own name called it
  an IC50 (biologically apt, since it inhibits viral replication), but
  this fork's naming convention uses `EC50_<STEM>` uniformly regardless
  of inhibitory/stimulatory direction (matching `EC50_LDN`, already named
  that way in the original for an equally inhibitory effect).
- **MET**: **one shared `EFFECT_MET`**, combined at each of its two
  points of use exactly as the original combined its own two identical
  ratios — `met_IL6 = 1.0 - EFFECT_MET` (IL-6 suppression) and
  `met_mito = EFFECT_MET * kMet_mito` (AMPK-mediated mito protection,
  `kMet_mito` kept as-is, a rate constant, not itself a Hill parameter).
  `EC50_MET` renamed from `EC50_met_IL6` (same value, 500 ng/mL) — the
  original's own name suggests it was intended only for the IL-6 pathway,
  but the code reuses the identical ratio (same `C_met`, same
  `EC50_met_IL6`) for the mito-protection term too, so one shared
  `EFFECT_MET` is a faithful rename of what the original actually
  computes, per the guide's "combine only at the point disease equations
  actually use them."
- **SERT**: `EFFECT_SERT = EMAX_SERT*C_SERT/(EC50_SERT+C_SERT)` — was
  `sert_5HT_eff = C_sert/(50.0+C_sert)`. `EC50_SERT = 50.0` is a **new
  named parameter**, promoted from a hardcoded literal the original wrote
  directly inline (its own comment already stated "EC50 ~50 ng/mL for
  SERT occupancy") — a promotion, not a fit, same pattern as
  `abdominal-aortic-aneurysm/Statin`'s NF-kB/ROS terms in
  `abdominal-aortic-aneurysm_refactor_notes.md`-equivalent precedent.
- **LDN**: `EFFECT_LDN = EMAX_LDN*C_LDN/(EC50_LDN+C_LDN)` — was
  `LDN_eff = C_LDN/(EC50_LDN+C_LDN)`; `EC50_LDN` already named exactly
  this in the original, unchanged (2.00 ng/mL).

**No Hill fitting was needed for any of the 4 compounds** — every effect
term was already a plain concentration ratio, consistent with the
exact-match verification below.

## A dead parameter, found and disclosed (not fixed)

`kLDN_MG = 1.50 /day` ("LDN EC50 for microglial suppression (PD)") is
declared in `$PARAM` but **never referenced anywhere** in `$ODE`,
`$MAIN`, or `$TABLE` — confirmed by grepping the whole original file for
`kLDN_MG`; it appears exactly once, in its own `$PARAM` declaration line.
LDN's actual microglial-suppression effect is computed entirely from
`EC50_LDN` (a separate, correctly-wired parameter) — `kLDN_MG` has zero
effect on any output regardless of its value. This looks like a
leftover/renamed-away parameter from an earlier draft (its own comment
describes exactly what `EC50_LDN` already does). Kept unchanged
(unrenamed, unused) in the refactored file, per "never invent or default
a PK parameter... don't add one the original didn't have" (equally: don't
remove one either) — disclosed here rather than silently dropped.

## Build defect: the original does not compile under mrgsolve 2.0.1

Three layered, file-wide defects unrelated to any compound's own PK —
full detail logged as
[`UPSTREAM_ISSUES.md` #68](../translations/UPSTREAM_ISSUES.md):

1. `$CMT` + `$INIT` jointly redeclare all 26 compartments (same family as
   #61/#63/#64/#65/#67).
2. `$CAPTURE C_nirm C_met C_sert C_LDN` duplicates 4 already-`$CMT` names
   (same family as #67's Defect 3).
3. `$TABLE` self-referentially re-declares `FSS`/`MoCA`/`SF36_PCS` via
   `capture NAME = NAME;` (same mechanism as `spinal-muscular-atrophy`'s
   #66).

**Fixes applied directly to the delivered
`pasc_mrgsolve_model_refactored.R`**, syntax-only and non-numeric: the
`$INIT` block is replaced by `<CMT>_0 = value;` assignments in `$MAIN`
(identical values); the standalone `$CAPTURE` line is dropped (the
refactor's own `$TABLE` captures give `C_NIRM`/`C_MET`/`C_SERT`/`C_LDN`
their own named entries anyway); and the three colliding locals are
renamed `FSS_calc`/`MoCA_calc`/`SF36_PCS_calc` (`SF36_PCS_calc`'s formula
updated to read `FSS_calc` in place of `FSS`), with each `capture` line
reading from the renamed local. None of this changes any numeric value —
proved by the exact-match verification below.

**One additional, non-blocking, unrelated finding:** `$CMT` declares
`V_GUT` (comment: "viral depot (gut absorption)") but no `dxdt_V_GUT`
line exists anywhere — mrgsolve issues a non-fatal audit warning
(`Missing differential equation(s): --| missing: dxdt_V_GUT`) and treats
it as a permanently-zero state. `V_GUT` is never dosed or read anywhere
else either — the viral kinetics submodel starts directly from
`V_PLASMA`. Inert dead state, not a behavioural bug; preserved unchanged.

## A naming subtlety specific to this file: local `$ODE` variable names cannot equal a `$TABLE` capture name

While building the refactored file, an initial draft declared `double
C_NIRM = CENT_NIRM;` and `double EFFECT_NIRM = ...;` directly inside
`$ODE` (for use in the PD equations there), then tried to also `capture
C_NIRM`/`capture EFFECT_NIRM` from `$TABLE`. This fails to compile with
exactly the same mechanism as Build Defect 3 above (and
`spinal-muscular-atrophy`'s #66): mrgsolve compiles `$PARAM`/`$CMT`/
`$MAIN`/`$ODE`/`$TABLE` doubles into **one shared member namespace**, not
block-scoped locals, so a `double C_NIRM` in `$ODE` and a `capture
C_NIRM = ...` in `$TABLE` collide as a redefinition of the same symbol.
Fixed by keeping `$ODE`'s own local variables lowercase
(`eff_nirm`/`eff_met`/`eff_sert`/`eff_ldn` — mirroring the original's own
local-variable style, `nirm_eff`/`met_IL6`/`sert_5HT_eff`/`LDN_eff`) and
reserving the exact `C_<STEM>`/`EFFECT_<STEM>` names for `$TABLE`'s
`capture` lines only, which recompute the identical formula from the
same current `$CMT` state at each output record (not a second,
independent calculation — same numbers, confirmed by the verification
below since `$TABLE` runs after `$ODE`'s integration for that record,
reading the same state `$ODE` used at that instant). This is disclosed
here as a genuine constraint of this mrgsolve build worth knowing before
attempting the same pattern in another file, not as an upstream defect
(nothing here was in the checked-in original, since the original never
attempted to capture a computed-effect variable at all, only the raw
`$CMT` concentrations).

## Verification

**Method.** Both the original's own model code (with the three build-
compat fixes above applied to an in-memory-only scratch copy, never to
the checked-in `pasc_mrgsolve_model.R`) and
`pasc_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
(`POST /model_manifest`, `POST /run_simulation`) at
`http://localhost:8007`, requests spaced ~2s apart, never more than one
in flight. The extracted refactored DSL was confirmed **byte-identical**
to the `_refactored.R`'s embedded quoted string except for one possessive
apostrophe in a comment (`"$ODE's own"` -> `"$ODE own"`), which cannot
appear inside an R single-quoted string literal and was rewritten without
changing meaning — same disclosed pattern as
`prostate-cancer/pc_refactor_notes.md`. Compartment order and 1-based
dosing indices are identical between the two models (nothing added,
removed, or reordered): `GUT_NIRM`=19, `CENT_NIRM`=20, `GUT_MET`=21,
`CENT_MET`=22, `GUT_SERT`=23, `CENT_SERT`=24, `GUT_LDN`=25, `CENT_LDN`=26
in both, confirmed from `/model_manifest`'s `outputPaths`.

**Scenarios run.** Both of the original's own most informative scenarios,
copied verbatim from the original's own R code (dose amounts, `ii`,
`addl`), full one-year horizon (`end=365, delta=1`, no shortening
needed — neither scenario approached the API's default `maxsteps`
budget):

- **S1** (untreated natural history): no dosing, all `use_*=0`.
- **S7** (full combination): all four compounds dosed at their own full
  regimens simultaneously — nirmatrelvir 300 mg BID x15d (`ii=0.5,
  addl=29`), metformin 500 mg BID x1yr (`ii=0.5, addl=728`), sertraline
  50 mg QD x1yr (`ii=1, addl=364`), LDN 4.5 mg QD x1yr (`ii=1,
  addl=364`), all four `use_*=1`.

For each scenario, all 12 shared disease-side outputs (`FSS`, `VO2max`,
`MoCA`, `POTS_HR`, `SF36_PCS`, `NfL_pg`, `mMRC`, `CRP_proxy`,
`Ddimer_out`, `Viral`, `Antigen`, `Reservoir`) and all 4 renamed PK
outputs (`C_nirm`/`C_NIRM`, `C_met`/`C_MET`, `C_sert`/`C_SERT`,
`C_LDN`/`C_LDN`) were compared point-by-point across the full time grid.

| Scenario | Dosing | n pts | Shared-output max abs diff | PK max abs diff |
|---|---|---|---|---|
| S1 Untreated | none | 366 | **0.0 exact** | **0.0 exact** (all 4 stay at 0) |
| S7 Full Combo | Nirm 300mg BID x15d + Met 500mg BID x1yr + Sert 50mg QD x1yr + LDN 4.5mg QD x1yr | 370 | **0.0 exact** | **0.0 exact** (peak `C_NIRM`~1209, `C_MET`~428, `C_SERT`~246, `C_LDN`~0.63 ng/mL, all matched bit-for-bit) |

**Result: bit-exact match on every output, every timepoint, both
scenarios.** This is a pure structural rename with no reordering of
floating-point operations that would introduce rounding drift (unlike,
e.g., `prostate-cancer`'s multi-year horizon, where reordered arithmetic
slowly compounded to a small but nonzero deviation over hundreds of days)
— here every renamed expression is evaluated in exactly the same order
and grouping as the original, so the solver sees bit-identical arithmetic
at every step. Fully consistent with the guide's tolerance for
Archetypes 1-3 ("pure structural reorganization... expect a near-exact
match").

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL
(and again against the exact string embedded in the delivered
`pasc_mrgsolve_model_refactored.R`, both return `200`): all `$PARAM`
entries are listed, including `F_NIRM, KA_NIRM, CL_NIRM, V1_NIRM,
EC50_NIRM, EMAX_NIRM, GAMMA_NIRM`; `F_MET, KA_MET, CL_MET, V1_MET,
EC50_MET, EMAX_MET, GAMMA_MET`; `F_SERT, KA_SERT, CL_SERT, V1_SERT,
EC50_SERT, EMAX_SERT, GAMMA_SERT`; `F_LDN, KA_LDN, CL_LDN, V1_LDN,
EC50_LDN, EMAX_LDN, GAMMA_LDN`. `C_<STEM>` and `EFFECT_<STEM>` (all 8)
are state-derived (`$TABLE`-computed from `$CMT` values, per the naming
subtlety above), so — per the same established reasoning in
`breast-cancer/bc_refactor_notes.md`,
`neonatal-hyperbilirubinemia/nhb_refactor_notes.md`, and
`prostate-cancer/pc_refactor_notes.md` (mrgsolve 2.0.1 compiles `$PARAM`
members as read-only references, so a value recomputed every timestep
from state cannot also live in `$PARAM`) — they cannot also be `$PARAM`
entries; all 8 appear in the manifest's `outputPaths` via `$CAPTURE`,
confirmed discoverable. All 26 compartments (`V_GUT` ... `CENT_LDN`) also
appear in `outputPaths` as ordinary compartments, in the same order/index
as the original.

No `.cpp` extraction file was left behind — extraction was in-memory
only, used to build the verification requests above and then discarded,
per the workflow guide.

## Anything else flagged

- No plot (`p1`-`p7`) in the R script references a renamed drug-specific
  compartment name except `p6`/`pk_data`, which read `C_nirm`/`C_met`/
  `C_sert`/`C_LDN` — updated to `C_NIRM`/`C_MET`/`C_SERT`/`C_LDN` in the
  refactored file, values unchanged (same columns, renamed only). All
  `ev()` calls (`ev_nirm`, `ev_met`, `ev_sert`, `ev_LDN`) updated from
  `cmt="A_<drug>"` to `cmt="GUT_<STEM>"`; every dose amount, `ii`, `addl`,
  and `time` is identical to the original.
- Every disease-side computation (viral kinetics, immune dysregulation,
  endothelial/coagulation, neuroinflammation, autonomic, mitochondrial)
  is untouched apart from reading the renamed drug-effect locals it
  already read before (`eff_nirm`/`met_IL6`/`met_mito`/`eff_sert`/
  `eff_ldn` in place of `nirm_eff`/`met_IL6`/`met_mito`/`sert_5HT_eff`/
  `LDN_eff`).
- Compartment order and 1-based indices are unchanged from the original
  for all 26 compartments, so no external code addressing this model by
  compartment number is affected.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
two mislabeled rows (`Neuroinflammation / BBB (LDN)`, `Viral Kinetics
(NIRM)`) were corrected to their real compound names (`LDN`, `NIRM`); the
two already-correct rows (`MET`, `SERT`) were left as classified. All
four remain "Redirect concentration (clean single site)" — each exposes
exactly one concentration variable (`C_<STEM>`) feeding exactly one
Hill-shaped effect term (`EFFECT_<STEM>`), confirmed by the exact-match
verification above.

## Discoverability fix

A corpus-wide discoverability audit found `C_NIRM`/`C_MET`/`C_SERT`/`C_LDN`
were not written as single contiguous `double C_<STEM> = <expr>;`
statements anywhere in the file. All four were exposed only via `capture
C_<STEM> = CENT_<STEM>;` lines inside `$TABLE` — a legitimate, working
pattern (not a bug: `CENT_<STEM>` is itself a concentration-state
compartment in this model, so the capture is a direct alias), but not
literal-text-discoverable by tooling that regexes for `double C_<STEM> =
...;`.

**Fix applied (Recipe B):** converted all four `capture C_<STEM> =
CENT_<STEM>;` lines in `$TABLE` to genuine `double C_<STEM> = CENT_<STEM>;`
declarations (identical formula/value — a rename of the declaration
mechanism, not a refit), and added an explicit `$CAPTURE C_NIRM C_MET
C_SERT C_LDN` block immediately after (this file previously had no
standalone `$CAPTURE` block — every other output, including
`EFFECT_NIRM`/`EFFECT_MET`/`EFFECT_SERT`/`EFFECT_LDN`, remains an inline
`capture NAME = expr;` line, untouched, unaffected by this fix). This
mirrors the fix already applied in
`autoimmune-polyendocrinopathy/aps_refactor_notes.md` for the same
`capture`-vs-`double` collision class.

**Verification:**
- `Rscript -e 'parse(...)'` succeeds with no errors.
- `grep`: `double C_<STEM> = ` and `EC50_<STEM>` both present for all four
  stems.
- qspserver `mrgsolve_api`: `/model_manifest` compiles the edited DSL and
  lists `C_NIRM`, `C_MET`, `C_SERT`, `C_LDN` in `outputPaths`.
  `/run_simulation` run against the file's own Scenario S2 (Nirmatrelvir
  300 mg BID×15d dosing, `use_nirm=1`, `end=365, delta=1`) shows
  `CENT_NIRM` and `C_NIRM` **exactly identical** (max abs diff = 0) between
  the pre-edit and post-edit DSL across all 367 reported rows, including
  the duplicate dose-instant row at t=0 — no dose-instant artifact at all
  for this compound, since `CENT_NIRM` is a state variable read directly,
  not a freshly-computed `$ODE`-local. A second run combining all four
  compounds' own dosing (Nirmatrelvir + Metformin + Sertraline + LDN
  simultaneously, matching the file's own Scenario S7 "Full Combo" doses)
  over a 60-day window confirmed `C_NIRM`/`C_MET`/`C_SERT`/`C_LDN` are also
  each exactly identical (max abs diff = 0) between pre- and post-edit.
