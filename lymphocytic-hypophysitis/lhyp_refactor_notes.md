# Refactor notes — `lymphocytic-hypophysitis/lhyp_mrgsolve_model.R`

**Scope of this pass.** Per `driver-patches/data/compound_perturbation_census.md`,
this file has exactly two rows, both classified **"Normalize duplicate
concentration sites, then redirect":**

- **Azathioprine (AZA)** — target: purine synthesis
- **Rituximab (RTX)** — target: CD20

Prednisolone is a third compound in this file (2-compartment oral PK, its
own `E_pred` Hill term feeding `ImmunoSupp` and several HPA/pituitary
terms) but has **no row of its own** in the census and was left completely
untouched: every `Pred_*`/`ka`/`Vc`/`Vp`/`CL_P`/`Q_P`/`F_pred`/`EC50_pred`/
`Emax_pred`/`E_pred`/`Cpred`/`Cpred_obs` name, value, and formula is
identical to the original. The entire disease side of the model (immune
cell dynamics, pituitary inflammation/function, HPA/HPT/HPG/GH-IGF1 axes,
prolactin, ADH) is likewise untouched except for the two forced,
mechanical substitutions where those equations read AZA's or RTX's effect
term by its old name.

## The assigned defect: duplicate concentration sites

Reading the original confirmed the census classification exactly. Both
AZA and RTX had their plasma concentration computed **independently, with
two different local names, in two different mrgsolve blocks that cannot
see each other's locals ($MAIN and $TABLE):**

- AZA: `double Caza = AZA_plasma / Vc_aza;` in `$MAIN` (feeds `E_aza`,
  which feeds `ImmunoSupp`) vs. `double Caza_obs = AZA_plasma / Vc_aza;`
  in `$TABLE` (reporting-only).
- RTX: `double Crtx = RTX_plasma / Vc_rtx;` in `$MAIN` (feeds `E_rtx`,
  used directly in `dxdt_Bn`/`dxdt_Bp`) vs.
  `double Crtx_obs = RTX_plasma / Vc_rtx;` in `$TABLE` (reporting-only).

Same formula, same numeric value at every timepoint (both blocks read the
same live state at the record they each evaluate), but two textually
separate definitions with no shared name — exactly the kind of thing that
silently drifts apart the next time someone edits one copy and not the
other, and exactly the case the guide's naming convention exists to
prevent ("that variable's *definition* is the single point where an
external covariate could later be substituted in" — currently there would
have to be two).

Normalized to **one name per compound, `C_AZA`/`C_RTX`**, predeclared once
in a new `$GLOBAL` block (the original had no `$GLOBAL` block at all) and
assigned identically at both of the original's own evaluation sites — the
same pattern already established in `hereditary-spherocytosis/
hsph_refactor_notes.md`, `cervical-cancer`, `breast-cancer`, `age-related-
macular-degeneration`, `membranous-nephropathy`, and `stable-angina`:
`$TABLE` cannot read a `$MAIN` local, and reading a `$GLOBAL` value another
block wrote earlier in the same interval risks a stale read (see the
guide's "keep a calculation in the block the original used it in" and
"dose-instant reporting artifact" sections), so rather than compute once
and share, the identical formula is kept at **both** of the original's own
call sites — now writing to one shared name instead of two independently
maintained ones. Concretely:

```c
$GLOBAL
double C_AZA, EFFECT_AZA, C_RTX, EFFECT_RTX;

$MAIN
C_AZA = CENT_AZA / V1_AZA;   // was: double Caza = AZA_plasma / Vc_aza;
C_RTX = CENT_RTX / V1_RTX;   // was: double Crtx = RTX_plasma / Vc_rtx;
EFFECT_AZA = EMAX_AZA * pow(C_AZA, GAMMA_AZA) / (pow(EC50_AZA, GAMMA_AZA) + pow(C_AZA, GAMMA_AZA));
EFFECT_RTX = EMAX_RTX * pow(C_RTX, GAMMA_RTX) / (pow(EC50_RTX, GAMMA_RTX) + pow(C_RTX, GAMMA_RTX));
...
$TABLE
C_AZA = CENT_AZA / V1_AZA;   // was: double Caza_obs = AZA_plasma / Vc_aza;
C_RTX = CENT_RTX / V1_RTX;   // was: double Crtx_obs = RTX_plasma / Vc_rtx;
```

Only the *concentration* was duplicated in the original (`Caza`/`Caza_obs`,
`Crtx`/`Crtx_obs`) — `E_aza`/`E_rtx` (now `EFFECT_AZA`/`EFFECT_RTX`) were
each computed exactly once, in `$MAIN` only, and were never re-derived in
`$TABLE`. So `EFFECT_AZA`/`EFFECT_RTX` are assigned once, in `$MAIN`,
matching where the original placed the equivalent calculation — not
reassigned in `$TABLE` — while `C_AZA`/`C_RTX` are reassigned in both
blocks, matching what the original actually did at each site. `Caza_obs`
and `Crtx_obs` are removed; `C_AZA`/`C_RTX` are captured directly instead
— a genuine single-site concentration *definition* (one formula, one
name) where the original had two independently maintained copies, even
though — as with every prior refactor that hit this same mrgsolve
constraint — the assignment statement itself still has to appear at each
of the original's two evaluation points.

## Archetype determination

**Azathioprine (AZA).** Two compartments in the original: `AZA_gut`
(oral depot) and `AZA_plasma` (central). Linear absorption (`ka_aza`) into
linear elimination (`CL_aza`/`Vc_aza`), no peripheral compartment.
**Archetype 3 minus peripheral (depot + central, linear).**

One quirk, preserved exactly: `F_aza = 0.50` (bioavailability) is declared
in `$PARAM` but **never referenced anywhere in `$ODE`/`$MAIN`/`$TABLE`**
(confirmed by grep — the only occurrence of `F_aza` in the whole DSL block
is its own `$PARAM` declaration). Bioavailability is instead applied
externally, by scaling the *dosed amount* in the R-side `aza_events()`
helper (`amt = daily_mg * 0.50`), not inside the model. `F_AZA` is kept as
a declared, renamed parameter (matching the naming convention's
expectation of an `F_<STEM>` entry) but is **still unused inside `$ODE`**,
exactly as in the original — adding it into the ODE would double-count
bioavailability against the R-side event scaling and change the model's
actual behavior, which the guide explicitly forbids ("never invent or
default a PK parameter" cuts both ways: don't silently "fix" an unused one
either). Disclosed here rather than silently normalized away.

**Rituximab (RTX).** Checked for TMDD per the task's own instruction,
since rituximab clinically undergoes target-mediated drug disposition via
CD20 binding on B cells. **The original does not model this**: `RTX_plasma`
is a single compartment with plain linear clearance
(`dxdt_RTX_plasma = -(CL_rtx/Vc_rtx) * RTX_plasma`) and no receptor,
complex, or binding-kinetics compartment exists anywhere in the file for
CD20 — B-cell depletion is represented purely as an algebraic Hill ratio
on plasma concentration (`E_rtx = Emax_rtx * Crtx / (EC50_rtx + Crtx)`)
multiplying the B-cell activation/plasma-cell-generation rates directly.
**Archetype 1 (no depot, single compartment, linear elimination).** No
TMDD structure existed to preserve; none was added (the guide's "none of
these fit"/"don't flatten" instruction is about not *removing*
mechanistic richness the original had — there was none here to keep).

## Renaming applied (values unchanged from the original)

| Original | Refactored | Value | Role |
|---|---|---|---|
| `AZA_gut` (cmt) | `GUT_AZA` | — | depot (oral) |
| `AZA_plasma` (cmt) | `CENT_AZA` | — | central |
| `ka_aza` | `KA_AZA` | 0.80 (1/h) | absorption rate |
| `Vc_aza` | `V1_AZA` | 30.0 (L) | central volume |
| `CL_aza` | `CL_AZA` | 12.0 (L/h) | clearance |
| `F_aza` | `F_AZA` | 0.50 | bioavailability — declared, **unused in `$ODE`** (see above), matching the original exactly |
| `EC50_aza` | `EC50_AZA` | 0.20 (mg/L equiv.) | Hill EC50 |
| `Emax_aza` | `EMAX_AZA` | 0.70 | Hill Emax |
| — (none) | `GAMMA_AZA` (new) | 1.0 | Hill exponent [implicit in the original] |
| `Caza` ($MAIN) / `Caza_obs` ($TABLE) | `C_AZA` | — | **the exposed concentration** (single, normalized site) |
| `E_aza` | `EFFECT_AZA` | — | rename, not a fit |
| `RTX_plasma` (cmt) | `CENT_RTX` | — | central (only compartment — no depot) |
| `Vc_rtx` | `V1_RTX` | 3.1 (L) | central volume |
| `CL_rtx` | `CL_RTX` | 0.015 (L/h) | clearance |
| `MW_rtx` | `MW_RTX` | 145000 (g/mol) | molecular weight (informational; never used in `$ODE`) |
| `EC50_rtx` | `EC50_RTX` | 0.05 (mg/L) | Hill EC50 |
| `Emax_rtx` | `EMAX_RTX` | 0.90 | Hill Emax |
| — (none) | `GAMMA_RTX` (new) | 1.0 | Hill exponent [implicit in the original] |
| `Crtx` ($MAIN) / `Crtx_obs` ($TABLE) | `C_RTX` | — | **the exposed concentration** (single, normalized site) |
| `E_rtx` | `EFFECT_RTX` | — | rename, not a fit |

`ImmunoSupp = 1.0 - (1.0-(1.0-E_pred))*(1.0-E_aza*0.5)` now reads
`(1.0-EFFECT_AZA*0.5)`; `dxdt_Bn`/`dxdt_Bp`'s `(1.0-E_rtx)` now read
`(1.0-EFFECT_RTX)`. No other equation changed. All parameter *values* are
copied verbatim from the original.

## Hill interface: rename, not a fit

Both `E_aza = Emax_aza*Caza/(EC50_aza+Caza)` and
`E_rtx = Emax_rtx*Crtx/(EC50_rtx+Crtx)` **are already** the canonical
`EMAX*C^gamma/(EC50^gamma+C^gamma)` shape with implicit `gamma=1`.
`EFFECT_AZA`/`EFFECT_RTX` are one-to-one renames with `GAMMA_AZA`/
`GAMMA_RTX = 1` added explicitly (`pow(C,1)` is mathematically identical
to `C`, not a behavioral change). No `nls()` fit was needed or performed
for either compound.

## When the original doesn't compile at all

The untouched original fails to compile under mrgsolve 2.0.1 with **two
independent, pre-existing defects**, both confirmed via
`POST /model_manifest` against the untouched original alone (no refactor
content involved):

**1. `$CMT` and `$INIT` jointly redeclare all 25 compartments.**
`$CMT` names and comments every compartment, then a separate `$INIT`
block re-declares the same 25 names with starting values. mrgsolve 2.0.1
treats `$INIT` as its own compartment-declaring block, not a companion to
`$CMT`:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: Pred_gut
  Pred_central Pred_periph AZA_gut AZA_plasma RTX_plasma Tn Te Tr Bn Bp APA
  PitInf PitFunc ACTH Cortisol TSH fT4 GH IGF1 FSH LH E2 PRL ADH
```

Logged as `translations/UPSTREAM_ISSUES.md` **#122** (same defect class
already logged for several other files in this corpus, e.g. #34, #114,
#117, #119).

**2. `$TABLE` ends in three bare, unheadered `capture NAME1 NAME2 ...`
lines** (lowercase, no preceding `$CAPTURE` marker), which mrgsolve 2.0.1
tries to compile as a call to an undeclared C++ function `capture(...)`
and rejects:

```
error: expected initializer before 'Caza_obs'
  653 |   capture Cpred_obs Caza_obs Crtx_obs PFS CA_ratio
```

Logged as `translations/UPSTREAM_ISSUES.md` **#123** (same defect class
already logged for `bronchiectasis/bex_mrgsolve_model.R`, issue #46).
None of the captured names duplicate a `$CMT` compartment name, so (unlike
issue #46) no names needed to be dropped — only the block marker itself.

Neither defect is inside AZA's or RTX's own PK/PD block, and neither was
fixed in the original, per the never-edit-upstream rule; the checked-in
`lhyp_mrgsolve_model.R` is completely untouched and still carries both
defects exactly as written.

**Fix applied directly to the delivered `lhyp_mrgsolve_model_refactored.R`**
(per the guide's settled policy — syntax-only, non-numeric, non-behavioral):
(a) the `$CMT` block is deleted; `$INIT` alone declares all 25
compartments, in the same order the original's `$CMT` block used —
declares no new compartment, changes no numeric value, and leaves
1-based dosing/output compartment indices unchanged; (b) each bare
`capture` line is given a real `$CAPTURE` header (`capture` → `$CAPTURE`,
three separate block headers, same names, same order — `Caza_obs`/
`Crtx_obs` are additionally absent from the refactored version of that
first line, but only because the AZA/RTX duplicate-site refactor itself
normalized them away, not because of this build-compat fix). Confirmed
via `POST /model_manifest` that the patched DSL compiles, and via
`POST /run_simulation` (below) that both fixes together reproduce a
patched-original whose output is byte-identical to the refactored file.

## Verification

**Method.** Both files' embedded `code <- '...'` DSL blocks were
mechanically extracted (find the quoted-string assignment, pull the
contents verbatim) and POSTed to the local qspserver `mrgsolve_api`
service at `http://localhost:8007` (`POST /model_manifest` then
`POST /run_simulation`), which compiles and runs each DSL block directly
with mrgsolve 2.0.1 server-side — no local R/mrgsolve install used. For
the "original" side of the comparison, the two disclosed build-compat
fixes above (#122, #123) were applied to a **scratch copy only** (never
to the checked-in `lhyp_mrgsolve_model.R`) so that a compiled baseline
exists to verify against at all. Requests were spaced ~2.2s apart and run
sequentially (never more than one in flight), respecting the service's
`max_concurrent_jobs: 2` limit and its documented history of crashing
under concurrent load.

One further build-compat issue surfaced only at the `/run_simulation`
step, not `/model_manifest`: submitting a multi-drug `dosing` array whose
entries are not sorted ascending by `time` (e.g. azathioprine's dosing
entry, constructed after several prednisolone-taper entries in R-source
order but starting *earlier* in simulated time) fails with `"the data set
is not sorted by time"`. This is a request-construction detail of how
this verification script built the `dosing` payload (the API's own
requirement, not a DSL defect) — fixed by sorting each scenario's combined
dosing list by `time` before submission; not logged as an upstream issue
since it is not a defect in either `.R` file.

**Scenarios run — all five of the file's own R-side scenario functions
(`pred_events`/`aza_events`/`rtx_events`), not invented ones, reproduced
via equivalent `dosing` entries (`ii`/`addl` reproducing each helper's
`seq(...)`-generated repeat-dose sequence exactly), full 2-year duration
(`end=17520h, delta=24h`, up to 737 points; no shortening needed — none
of the five approached the API's default solver-step budget):**

1. `S1_no_treatment` — no dosing (natural history)
2. `S2_pred_taper` — prednisolone-only 6-phase taper (matches R's `e_s2`)
3. `S3_pred_aza` — prednisolone taper + azathioprine 150 mg/day from day 30
   for 600 days (matches R's `e_s3`) — **exercises AZA**
4. `S4_rtx` — rituximab 1000 mg × 2 infusions, 180 days apart, + short
   prednisolone course (matches R's `e_s4`) — **exercises RTX**
5. `S5_pred_lowdose` — prednisolone 10 mg/day throughout (matches R's `e_s5`)

**Result: exact match, max abs diff 0.0**, across all 38 shared outputs
(compared by mapping the original's `AZA_gut`/`AZA_plasma`/`RTX_plasma`/
`Caza_obs`/`Crtx_obs` onto the refactored file's `GUT_AZA`/`CENT_AZA`/
`CENT_RTX`/`C_AZA`/`C_RTX`, plus every identically-named disease-state
compartment and derived `$CAPTURE` output — 33 identically-named +
5 renamed), at every timepoint, in **all five scenarios**. This is a pure
structural reorganization (rename + duplicate-site normalization, no
Hill-fitting), consistent with the guide's tolerance table for Archetypes
1–3: "expect a near-exact match... anything beyond floating-point-scale
deviation means a bug." Here the match is bit-exact (0.0), not merely
floating-point-scale.

`EFFECT_AZA`/`EFFECT_RTX` (new outputs, not present in the original at
all) were sanity-checked directly: in `S3_pred_aza`, `EFFECT_AZA` is
nonzero and `EFFECT_RTX` is exactly zero throughout (no RTX dosed); in
`S4_rtx`, `EFFECT_RTX` reaches 0.8999 (approaching `EMAX_RTX=0.90`, as
expected under sustained high exposure) while `EFFECT_AZA` is exactly
zero throughout (no AZA dosed) — confirming each Hill interface fires
only for its own compound's dosing, per the guide's naming-convention
item 4.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed,
build-compat-fixed DSL: `KA_AZA`, `V1_AZA`, `CL_AZA`, `F_AZA`, `EC50_AZA`,
`EMAX_AZA`, `GAMMA_AZA`, `V1_RTX`, `CL_RTX`, `MW_RTX`, `EC50_RTX`,
`EMAX_RTX`, `GAMMA_RTX` all appear in the manifest's `parameters` with
their original numeric defaults (83 parameters total, including every
untouched prednisolone/disease-model parameter). `C_AZA`/`EFFECT_AZA`/
`C_RTX`/`EFFECT_RTX` are state-derived (recomputed each interval from live
compartment state), so — per the same reasoning already established for
`C_LARA`/`EFFECT_LARA` in `celiac-disease/cd_refactor_notes.md`, `C_MIT`/
`EFFECT_MIT_*` in `hereditary-spherocytosis/hsph_refactor_notes.md`, and
~25 other refactors in this corpus that hit the same mrgsolve constraint
(`$PARAM` members are read-only inside `$ODE`/`$TABLE`) — they cannot also
be `$PARAM` entries; both are listed in `$CAPTURE` and confirmed present
in the manifest's `outputPaths`, discoverable as **outputs**. `GUT_AZA`/
`CENT_AZA`/`CENT_RTX` also appear in `outputPaths` as ordinary
compartments, at the same 1-based positions (4, 5, 6) the original's
`AZA_gut`/`AZA_plasma`/`RTX_plasma` occupied.

No `.cpp` extraction file, patched-scratch DSL, or verification script was
left behind — all were scratch-only, used to build the requests above and
then discarded.

## Anything else flagged

- The R-side `aza_events()`/`rtx_events()` helpers were updated only where
  they name the renamed compartment (`cmt = "AZA_gut"` → `"GUT_AZA"`,
  `cmt = "RTX_plasma"` → `"CENT_RTX"`), same dosing amounts, same timing,
  same five scenarios throughout. `pred_events()` and every prednisolone
  reference are untouched.
- Plot 6 (drug plasma concentrations) was updated to read `C_AZA`/`C_RTX`
  in place of the removed `Caza_obs`/`Crtx_obs`; its output (which
  scenario shows which trace) is unchanged.
- `lhyp_shiny_app.R` and `lhyp_qsp_model.dot`/`.svg`/`.png` are untouched
  — out of scope for this refactor (only the mrgsolve model file and its
  two required siblings are deliverables).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: both
`lymphocytic-hypophysitis` rows (Azathioprine, Rituximab (RTX)) filled in
with their exposed concentration/effect names and a notes summary per the
above.
