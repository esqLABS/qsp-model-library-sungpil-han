# Refactor notes — `t1dm_mrgsolve_model.R` (insulin + teplizumab)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2) and the
existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
("INS", classified "Redirect concentration (clean single site)"; "TEP",
classified "Redirect concentration inside #define macro"), only insulin's
and teplizumab's PK/effect blocks were rewritten, each independently named
per its own stem (`INS`, `TEP`) — not collapsed into a shared term. This is
a two-compound file with no other modeled drugs, so there is nothing else
in the pharmacological sense to leave untouched — but every other
compartment (beta-cell/immune, glucose/glucagon, C-peptide, HbA1c, CGM/APC,
meal absorption) and the entire R scenario-runner layer were left
byte-for-byte identical except for the handful of lines that reference the
renamed identifiers by name (listed below). Confirmed by `diff`, see "Diff
scope confirmation" below.

## Archetype determination

### INS (insulin) — bespoke

Insulin's PK/PD block is structurally close to several archetypes but does
not cleanly match any single one:

1. **A two-step sequential SC absorption depot** (`SC1` hexamer → `SC2`
   dimer/monomer, renamed `GUT1_INS`/`GUT2_INS`), feeding a **central
   compartment already written directly in concentration units** (mU/L,
   `Ic` → `CENT_INS`) rather than amount-over-volume — the same
   "already-a-concentration" pattern documented in
   `diabetic-ketoacidosis/dka_refactor_notes.md` for `CENT_INS` there. This
   part is close to Archetype 3 (depot + central + peripheral), except
   Archetype 3's worked example has one depot, not two.
2. **A genuine peripheral distribution compartment** (`Ip` → `PERI_INS`),
   exchanging with the central pool via micro rate constants `k12_ins`/
   `k21_ins` (renamed `K12_INS`/`K21_INS`). The guide's Archetype 2
   instructs "always rewrite to CL/Q/V" for a two-compartment linear PK
   model — **not done here**, deliberately: the original declares a
   "peripheral volume" `Vp_ins` (renamed `V2_INS`) but never uses it in
   `$ODE` at all (see "Dead parameters" below), and `k12_ins*Vc_ins`
   (0.04×3.5=0.14) does not equal `k21_ins*Vp_ins` (0.02×2.5=0.05) — so
   there is no consistent `Q = k12*V1 = k21*V2` to derive. Forcing a CL/Q/V
   reparameterization would have required inventing a volume relationship
   the original does not have, which the guide explicitly warns against
   ("don't add or remove compartments/parameters to make it fit"). The
   micro-rate-constant form was kept as-is, renamed only.
3. **A fused endogenous+exogenous+device-infusion flux into one central
   pool.** `CENT_INS` receives insulin from three sources added into the
   same `dxdt_CENT_INS` line: the SC absorption flux, glucose-driven
   endogenous beta-cell secretion (`secR`, computed in `$MAIN` from `Bm`
   and `Gp`), and the closed-loop-pump (APC) infusion term. This is
   structurally the same pattern as diabetic-ketoacidosis's `Riv = INS_IV +
   sec` fused flux (`dka_refactor_notes.md`, point 4) — preserved exactly,
   renamed only. `secR`/`secRate`/`IbasalEn` (the endogenous-secretion
   machinery) and the APC controller (`APC_mode`/`APC_ON`/`Kp_apc`/
   `Ki_apc`/`tgt_mid`/`APC_integral`) are physiological/device
   infrastructure outside insulin's own PK-naming-convention scope, left
   completely unrenamed — matching the DKA precedent's treatment of
   `INSB`/`SECMAX`/`KGSEC`/`PORTF`.
4. **Three disease-facing effect axes, not one**, all linear (non-Hill) in
   the same lagged effect compartment or in the raw central concentration:
   - `EFFECT_INS_EGP` — suppresses hepatic glucose output, reads the
     lagged effect compartment `EFF_INS` (was `X_ins`)
   - `EFFECT_INS_UP` — stimulates glucose uptake/disposal, also reads
     `EFF_INS`
   - `EFFECT_INS_GCG` — suppresses glucagon secretion, reads `C_INS`
     (the raw central concentration) **directly, not through the lagged
     effect compartment** — a genuine, distinct pharmacodynamic detail of
     the original (glucagon responds to immediate plasma insulin; glucose
     utilisation/EGP respond to the delayed "remote" action site) and
     preserved exactly, not flattened into a single shared term.

   None of these three is a saturating ratio in the original — each is a
   **plain linear multiple** of its input (`K_INS_EGP * EFF_INS`,
   `K_EFF_INS * EFF_INS`, `K_INS_GCG * (C_INS/15.0)`), with no `EC50`-style
   denominator anywhere. The guide's Hill interface offers two paths: rename
   an already-Hill-shaped term, or fit an Emax/Hill curve to ODE-derived
   kinetics. **Neither applies here** — the original term is neither
   already a saturating ratio nor does it emerge from otherwise-unwritable
   ODE-solved kinetics; it is simply linear and unbounded by construction.
   A Hill/Emax formula *always* saturates as its input grows (it can never
   equal an unbounded linear function for any finite `Emax`/`EC50`/`gamma`),
   so forcing one here would be a structural misrepresentation, not a
   rename — exactly the guide's warning ("a clean, non-standard structure
   beats a standard structure that's wrong"). Each term was therefore kept
   as a named, isolated, linear function of the compound's own concentration
   or effect-site state — satisfying the guide's actual requirement ("the
   compound's effect on disease is expressed as one named function of that
   concentration") without inventing saturation the original doesn't have.

Per the guide's "None of these fit" fallback: renamed to the convention,
isolated (insulin's PK/PD lives in one contiguous `$ODE` block, plus one
line each inside the EGP/Rd calculations in `$MAIN` and the glucagon block
in `$ODE`, all updated consistently), and documented here rather than
forced into Archetype 1/2/3's shape or an artificial Hill fit.

**A parameter reused for two distinct physiological roles.** `p2U`
(renamed `K_EFF_INS`) is used both as the effect-compartment's own
first-order rate constant (`dxdt_EFF_INS = K_EFF_INS * (C_INS - EFF_INS)`)
**and**, multiplied against `EFF_INS`, as `EFFECT_INS_UP`'s sensitivity
coefficient. This dual role exists in the original (same parameter, same
numeric value, both places) and is preserved exactly under the renamed
identifier — flagged here as a pre-existing modeling quirk, not something
introduced by the rename, and not a compile-blocking defect so not logged
in `UPSTREAM_ISSUES.md`.

### TEP (teplizumab) — Archetype 1

No depot, single compartment (`Ctep` → `CENT_TEP`), linear elimination
(`dxdt_Ctep = -(CL_tep/Vc_tep)*Ctep`) — a clean Archetype 1 match, with the
same "already-a-concentration" state-storage pattern as `CENT_INS` (no
`/V1_TEP` division needed in `C_TEP`'s definition; `V1_TEP` is used only
inside the elimination-rate term `CL_TEP/V1_TEP`, exactly mirroring
Archetype 1's own worked example `dxdt_CENT_TCZ = -CL_TCZ * C_TCZ;`).

**The census's classification does not describe the actual redirect
site.** The TEP row is labeled "Redirect concentration inside #define
macro". The file's `$GLOBAL` block does contain `#define TEPLIZ (Ctep >
0.001)`, and this macro does reference the compound's concentration — but
`TEPLIZ` is **never referenced anywhere else in the file** (confirmed by
`grep -n "TEPLIZ" t1dm_mrgsolve_model.R`, which matches only its own
declaration line). It is dead code, both before and after this rename. The
actual, and only, places teplizumab's concentration is read are two plain
`$ODE` lines (`dxdt_CTL`'s `- kCTL_tep * Ctep * CTL_pos` and `dxdt_Treg`'s
`+ kTreg_tep * Ctep`), which is exactly the "clean single site" shape the
*other* census classification label describes. Per the guide's instruction
that a gate/classifier mismatch is "a question, not a verdict" — this
looks like the classifier's static-text scan matched `Ctep`'s appearance
inside the `#define` line and attributed the redirect site to it, rather
than to the two live `$ODE` consumers. The refactor was built around the
actual (live) redirect site, `C_TEP`, defined once in `$ODE` and read by
both disease equations; the dead `TEPLIZ` macro was renamed for
consistency (`Ctep`→`CENT_TEP` inside it) but otherwise left exactly as
inert as it was.

**Teplizumab's two disease-facing effect terms are also linear, not
Hill-shaped**, for the same reason as INS's: `EFFECT_TEP_CTL = KCTL_TEP *
C_TEP` and `EFFECT_TEP_TREG = KTREG_TEP * C_TEP` are plain proportional
(unbounded) rate contributions in the original, with no EC50/saturation
term. Kept linear, renamed only — not fit to a Hill curve, for the same
structural reason given above for INS (a Hill formula cannot represent an
unbounded linear term without misrepresenting it). No TMDD/receptor-binding
state exists for teplizumab in this model (its mechanism — depleting
pathogenic CTLs while expanding Tregs — is captured entirely as a linear
rate contribution to the immune-cell ODEs, not as an explicit
antibody–receptor binding compartment), so Archetype 4 does not apply
either.

## Naming applied

| Role | Original | Refactored |
|---|---|---|
| **INS** absorption depot 1 (hexamer) | `SC1` | `GUT1_INS` |
| **INS** absorption depot 2 (dimer/monomer) | `SC2` | `GUT2_INS` |
| **INS** central compartment (already a concentration, mU/L) | `Ic` | `CENT_INS` |
| **INS** peripheral compartment | `Ip` | `PERI_INS` |
| **INS** lagged effect/remote compartment | `X_ins` | `EFF_INS` |
| **INS** absorption rate 1 | `ka1` | `KA1_INS` |
| **INS** absorption rate 2 | `ka2` | `KA2_INS` |
| **INS** clearance | `CL_ins` | `CL_INS` |
| **INS** central volume | `Vc_ins` | `V1_INS` |
| **INS** peripheral volume (declared, unused — see below) | `Vp_ins` | `V2_INS` |
| **INS** central→peripheral micro rate constant | `k12_ins` | `K12_INS` |
| **INS** peripheral→central micro rate constant | `k21_ins` | `K21_INS` |
| **INS** exposed concentration | *(none — `Ic` read directly)* | `C_INS` (new `$ODE`-local double, `= CENT_INS`) |
| **INS** effect-compartment rate constant / Rd sensitivity coefficient (dual role, see above) | `p2U` | `K_EFF_INS` |
| **INS** EGP-suppression coefficient | `kp1` | `K_INS_EGP` |
| **INS** glucagon-suppression coefficient | `kGcg_ins` | `K_INS_GCG` |
| **INS** disease-facing effect terms | *(inline, unnamed)* | `EFFECT_INS_EGP`, `EFFECT_INS_UP`, `EFFECT_INS_GCG` |
| **INS** intermediate absorption flux | `Ins_from_SC` | `Abs_INS` |
| **TEP** central compartment (already a concentration, μg/mL) | `Ctep` | `CENT_TEP` |
| **TEP** clearance | `CL_tep` | `CL_TEP` |
| **TEP** central "volume" | `Vc_tep` | `V1_TEP` |
| **TEP** exposed concentration | *(none — `Ctep` read directly)* | `C_TEP` (new `$ODE`-local double, `= CENT_TEP`) |
| **TEP** CTL-suppression coefficient | `kCTL_tep` | `KCTL_TEP` |
| **TEP** Treg-expansion coefficient | `kTreg_tep` | `KTREG_TEP` |
| **TEP** disease-facing effect terms | *(inline, unnamed)* | `EFFECT_TEP_CTL`, `EFFECT_TEP_TREG` |
| Dead macro referencing TEP's concentration (unused before and after) | `#define TEPLIZ (Ctep > 0.001)` | `#define TEPLIZ (CENT_TEP > 0.001)` |
| Endogenous secretion, APC controller, glucose/glucagon/beta-cell/immune/HbA1c/CGM/meal machinery | *(unchanged)* | *(unchanged — outside INS/TEP's own naming-convention scope)* |

**Both compounds' effect terms are renames, not refits.** The original's
`kp1*ins_eff`, `p2U*X_ins`, `kGcg_ins*(Ic/15.0)`, `kCTL_tep*Ctep`, and
`kTreg_tep*Ctep` are all already plain linear multiples with no `EC50`-style
denominator — pulling them out as named `EFFECT_<STEM>_<AXIS>` doubles with
the original's own coefficients changes zero numeric behaviour, confirmed
exactly by the verification below. No `EMAX`/`EC50`/`GAMMA` parameters were
invented for either compound, since neither original term has a saturating
shape to name (see "Archetype determination" above for why forcing one
would misrepresent the model, not just rename it).

**Dead parameters, found and left untouched (not upstream-logged — not
compile-blocking).**
- `V2_INS` (`Vp_ins`): declared "Peripheral volume (L)" but never
  referenced anywhere in `$ODE`/`$MAIN`/`$TABLE` — `PERI_INS`'s exchange
  with `CENT_INS` uses the micro rate constants `K12_INS`/`K21_INS`
  directly, with no volume term at all. Kept declared and unused, matching
  the guide's "don't delete a parameter the original didn't ask to
  delete" converse of "don't invent one."
- `ka_type`: declared "Current bolus absorption rate", set via
  `param(mod, ka_type = ...)` in every scenario function, but **never
  referenced anywhere in `$ODE`** — the model always uses the fixed
  `KA2_INS` (was `ka2`) for SC2→central absorption regardless of which
  "bolus type" a scenario nominally selects. Same dead-parameter shape as
  diabetic-ketoacidosis's `KA_REG` (`dka_refactor_notes.md`). Left as-is,
  name unchanged (not part of the active naming convention).
- `kp2` (declared "EGP suppression by portal insulin"): never referenced
  anywhere. Left unchanged (unrenamed) since it plays no role in any
  insulin-effect equation actually wired into the model — renaming a name
  that is never used would just add unnecessary diff noise.

**A pre-existing duplicate/dead-code assignment, found and preserved
exactly.** The original's `dxdt_Ic` is assigned **twice** in a row: the
first assignment includes a `* (-1.0)` sign-correction stub the author's
own comment flags as wrong ("APC infusion sign correction"), then a second,
differently-signed assignment immediately follows and silently overwrites
it (C++ semantics: the last assignment to a variable in a block wins, so
the first is dead code with zero effect on any simulated output). This is
**not** a compile blocker (both are syntactically valid statements), so it
is not in `UPSTREAM_ISSUES.md`, but it is worth flagging: the refactored
`dxdt_CENT_INS` preserves the identical double-assignment structure
(renamed) rather than "cleaning it up," per the guide's instruction that a
rename-only pass must not silently change behaviour or structure, even
where the structure contains harmless dead code.

**A pre-existing dimensional inconsistency in TEP's dosing, found and
preserved exactly (not upstream-logged — not compile-blocking).**
`tepliz_dose()` adds the raw dose amount (μg) directly into `CENT_TEP`,
which is nominally a **concentration** state (μg/mL) — i.e. the dose is
applied as if the distribution volume were exactly 1 mL, even though
`V1_TEP` (`Vc_tep` = 5 L) is declared and used elsewhere (only in the
elimination-rate term `CL_TEP/V1_TEP`). This is a genuine scientific
inconsistency in the original (the resulting post-dose `CENT_TEP` bump is
not actually `dose/V1_TEP`), but it is not a compile-time defect, and
"fixing" it would change simulated numeric behaviour — outside the scope
of a rename-only refactor. Preserved exactly; flagged here for whoever
next revisits this model's dosing.

## qspserver compatibility

- `C_INS`, `C_TEP`, and each `EFFECT_<STEM>_<AXIS>` are `$ODE`- or
  `$MAIN`-local `double`s, not `$PARAM` entries — mrgsolve 2.0.1 passes
  every `$PARAM` value into `$ODE` as a `const double&`, so assigning to
  one (`C_INS = CENT_INS;`) is a hard compile error, the same constraint
  documented in `dka_refactor_notes.md`'s "qspserver compatibility"
  section (and independently rediscovered here while testing an
  init-override workaround — see `UPSTREAM_ISSUES.md` #42).
- Discoverability is instead satisfied via a new `$CAPTURE` block. This
  file's `$TABLE` idiom is bare `double NAME = EXPR;` (no `capture` macro
  and no pre-existing `$CAPTURE` block at all — confirmed by
  `/model_manifest` on the untouched original: its `outputPaths` list only
  the 20 `$CMT` compartments, none of the disease's own `$TABLE` doubles
  like `Glucose_mgdL`). Since a `$TABLE`-scoped block cannot see `$ODE`'s
  local `C_INS`/`EFFECT_INS_*` variables directly (each mrgsolve block
  compiles into its own function), the refactored file recomputes each one
  from the underlying `$CMT` states (which are visible everywhere) under
  `_report`-suffixed names — `C_INS_report`, `EFFECT_INS_EGP_report`,
  `EFFECT_INS_UP_report`, `EFFECT_INS_GCG_report`, `C_TEP_report`,
  `EFFECT_TEP_CTL_report`, `EFFECT_TEP_TREG_report` — then lists all seven
  in a new `$CAPTURE` block, following the identical pattern used in
  `dka_refactor_notes.md`. Confirmed present in `/model_manifest`'s
  `outputPaths` (see verification below) and listed in the updated census
  row. (This refactor did not attempt to fix the pre-existing gap that
  none of the disease's *own* `$TABLE` outputs like `Glucose_mgdL` are
  discoverable via the API either — that is outside INS/TEP's scope and
  was left untouched, consistent with "leave everything else alone.")
- All 20 `$CMT` names (six renamed for INS/TEP) and every `$PARAM` entry
  (INS's 10 renamed + TEP's 4 renamed) are confirmed present and correctly
  defaulted via `/model_manifest` — see below.
- `$SET end/delta` is not present in this model at all (the original never
  declared one; both models rely entirely on `mrgsim(end=..., delta=...)`
  call-site arguments in the R scenario functions) — so there is nothing
  to note about it being non-load-bearing for API runs; every verification
  run below supplied `time.end`/`time.delta` explicitly via the request.

## Diff scope confirmation

`git diff --no-index type1-diabetes/t1dm_mrgsolve_model.R
type1-diabetes/t1dm_mrgsolve_model_refactored.R` touches exactly: the INS
and TEP `$PARAM` entries (renamed/relocated, `kGcg_ins` moved into the
Insulin PK section as `K_INS_GCG`), the six renamed `$CMT`/`$INIT` lines,
the `TEPLIZ` macro's one internal reference, the `$MAIN` EGP/Rd block (two
new named `EFFECT_INS_*` doubles replacing the anonymous `ins_eff`/inline
`p2U*X_ins`), the `$ODE` beta-cell/immune block's two TEP-effect lines
(plus the two new lines that compute `C_TEP`/`EFFECT_TEP_*` up front), the
entire insulin-PK `$ODE` block (renamed, `C_INS` added), the glucagon
block's one line, the teplizumab elimination line, the `$TABLE` block's two
renamed report lines plus seven new `_report`/`$CAPTURE` lines, and the R
helper/scenario functions' six compartment-name string literals plus their
`init(...)` overrides (`SC1`/`SC2`→`GUT1_INS`/`GUT2_INS`, `Ic`→`CENT_INS`,
`Ctep`→`CENT_TEP`). Nothing else — the beta-cell/immune ODEs, C-peptide,
meal absorption (Dalla Man), plasma/tissue glucose, HbA1c, CGM/APC control,
mean-glucose running average, all plotting code, CGM-metrics calculation,
and the sensitivity-analysis loop are byte-for-byte identical.

## Verification

**Method.** Both files' embedded mrgsolve DSL blocks were mechanically
extracted (`code <- '...'` → bare DSL text) and run through the qspserver
`mrgsolve_api` container (`http://localhost:8007`, confirmed healthy),
using `POST /model_manifest` to confirm both compile and expose the right
parameters/outputs, and `POST /run_simulation` to compare every shared
compartment output across four scenarios built directly from the original
file's own dosing constructors (`meal_event`/`bolus_event`/`tepliz_dose`)
and its scenario functions' own parameter overrides:

1. **V0 — undosed baseline.** Default parameters/initial conditions
   (matching both files' own `$INIT`), 1000 min, delta 5.
2. **V1 — insulin, MDI-style.** Reproduces `run_MDI()`'s own per-day
   building blocks exactly (3× daily 75 g meal + 8 U rapid bolus at
   07:00/12:00/18:30, plus 20 U basal at 22:00, into the same compartments
   `meal_event`/`bolus_event` use), extended across 2 days, with
   `APC_mode=0, ka_type=0.05` (`run_MDI`'s own overrides). Minute-scale
   time grid (0–2880 min, delta 5), matching the native unit of insulin's
   own rate constants.
3. **V2 — teplizumab, full 14-dose course.** Reproduces `tepliz_dose()`'s
   own dosing exactly: 14 IV doses of 51×70=3570 μg, `ii=1440` (1 day
   apart, `tepliz_dose`'s own `time = day*1440`), with
   `CTL0=0.7, kDest=0.0015, APC_mode=0` (`run_teplizumab()`'s own
   overrides). Run to 20000 min (~13.9 days, comfortably covering all 14
   scheduled doses), delta 60.
4. **V3 — combined insulin + teplizumab.** One insulin bolus + one meal at
   t=0 plus a 4-dose teplizumab course (`ii=1440, addl=3`) in the same run,
   `CTL0=0.8, kDest=0.002, APC_mode=0` (`run_tepliz_MDI()`'s own
   overrides), 4 days, delta 10 — confirms no cross-contamination between
   the two compounds' renamed identifiers.

Every shared `$CMT` compartment state (14 disease-side states common to
both files by name, plus the 6 renamed INS/TEP states compared
name-to-name — `SC1`↔`GUT1_INS`, `SC2`↔`GUT2_INS`, `Ic`↔`CENT_INS`,
`Ip`↔`PERI_INS`, `X_ins`↔`EFF_INS`, `Ctep`↔`CENT_TEP`) was compared across
the full output time grid for each scenario — compartment states are
always present in mrgsolve's `/run_simulation` output regardless of
`$CAPTURE`, so this comparison is independent of the `$TABLE`
discoverability gap noted above.

**Result: exact match. Maximum absolute and relative deviation observed =
0.0 (bit-identical) across all 4 scenarios and all 20 compared state
variables.** This is the expected outcome for a pure structural
reorganization/rename with no Hill-refitting, per the guide's tolerance
table.

**A qspserver infrastructure flakiness was observed once, unrelated to
either model.** One `/run_simulation` call during V3 returned a transient
`500` ("worker was killed by signal 6" / "corrupted size vs. prev_size")
on the first attempt, then succeeded with a bit-identical result on an
immediate retry with the identical payload — the same class of transient
container-level flakiness documented in `dka_refactor_notes.md`. Not
logged in `UPSTREAM_ISSUES.md` (that file is for defects in the
disease-model corpus, not the shared qspserver service).

**What could not be verified: the R scenario functions' own `init(...)`
state overrides.** Every scenario function in this file (`run_untreated`,
`run_MDI`, `run_CSII`, `run_HCL`, `run_teplizumab`, `run_tepliz_MDI`) calls
mrgsolve's `init(...)` to override starting values for `Bm`/`CTL`/`Treg`/
`Gp`/`Ic`(`CENT_INS`)/`Ctep`(`CENT_TEP`) before dosing (e.g.
`init(Bm = 0.5, CTL = 0.7, Treg = 0.4, Gp = 100, Ctep = 0)` in
`run_teplizumab()`). The qspserver `SimRequest` schema has no field for
arbitrary state-initial-condition overrides for an mrgsolve model — its
`events` field is documented deSolve-engine-only and confirmed inert for
mrgsolve `model_content` (same finding as `dka_refactor_notes.md`); a first
attempt to route around this by declaring the moved-`$INIT` values as
`$PARAM`-style `<CMT>_0` entries (which mrgsolve does support as an
initial-condition idiom in some configurations) failed to compile here,
because mrgsolve 2.0.1 already auto-generates a mutable `<CMT>_0` reference
for every compartment inside `$MAIN`, and a same-named `$PARAM` entry
collides with it (`conflicting declaration 'const double& Gt_0'` — see
`UPSTREAM_ISSUES.md` #42 for the full error and reasoning). All four
verification scenarios above therefore use each model's own **default**
`$INIT` starting values (identical in both files) combined with the
scenarios' own **parameter overrides and dosing** (both of which the API
does support), rather than their `init(...)` state overrides. Running from a range of
custom starting conditions was therefore not directly tested — but the
underlying ODE right-hand sides are confirmed byte-for-byte equivalent
(via V0's zero-dosing, zero-override baseline match, and V1–V3's
dosing-driven matches), so there is no mechanism by which a different
starting `CTL`/`Treg`/`Bm`/`Gp` value could cause the two models to
diverge — the renamed identifiers evaluate the identical arithmetic
regardless of the state they start from. Disclosed as an inference from
verified building blocks, not a directly-executed check of every named
scenario's exact starting point, consistent with the DKA precedent's
equivalent disclosure for its own closed-loop scenarios.

## Anything else worth flagging

- This is the second file in this fork's refactor history (after DKA)
  where the compound-effect terms are linear rather than Hill/Emax-shaped
  in the original — but unlike DKA (where the *shape* was already
  `Emax*C/(C+IC50)` and thus a trivial rename), here the shape is
  genuinely unbounded/proportional, making a "rename to EMAX/EC50/GAMMA"
  form actively wrong rather than merely unnecessary. Reported as such
  rather than silently wrapping a Hill formula around a linear term to
  match the guide's naming table more closely.
- Two independently-modeled compounds in one file, both needing full
  compound-specific naming (`INS`, `TEP`) but sharing no compartments or
  effect terms — verification scenario V3 specifically exercises both
  together to confirm the renames don't collide or cross-read each other's
  state.
