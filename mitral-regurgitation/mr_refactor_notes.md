# Refactor notes — `mr_mrgsolve_model.R` (furosemide, sacubitril/valsartan, beta-blocker, spironolactone/canrenone, SGLT2 inhibitor, nitroprusside, dobutamine)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), all eight
compounds this 43-ODE model treats were rewritten, per their existing rows
in [`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all classified "Redirect concentration (clean single site)"): **BB**
(beta-blocker, stem `BB`), **DOB** (dobutamine, stem `DOB`), **FUR**
(furosemide, stem `FUR`), **MRA** (spironolactone/canrenone, stem `MRA`),
**SAC** (sacubitril, stem `SAC`), **SG** (SGLT2 inhibitor, stem `SG`),
**SNP** (sodium nitroprusside, stem `SNP`), and **VAL** (valsartan, stem
`VAL`). These are the model's only compounds — nothing was left untouched.

## Archetype per compound

| Compound | Archetype | PK compartments | Notes |
|---|---|---|---|
| BB  | 3 minus peripheral (depot+central, linear) | `GUT_BB`/`CENT_BB` | 3 named effects sharing one occupancy fraction |
| VAL | 3 minus peripheral | `GUT_VAL`/`CENT_VAL` | 1 named effect + 1 vestigial Emax (see below) |
| SAC | 3 minus peripheral | `GUT_SAC`/`CENT_SAC` | 3 named effects sharing one occupancy fraction |
| MRA | 3 minus peripheral | `GUT_MRA`/`CENT_MRA` | 2 named effects sharing one occupancy fraction |
| SG  | 3 minus peripheral | `GUT_SG`/`CENT_SG` | 2 named effects sharing one occupancy fraction |
| FUR | 3 minus peripheral, plus a bespoke braking state | `GUT_FUR`/`CENT_FUR`/`BRAKE_FUR` | 2 named effects, one Hill-shaped, one a clamped linear ratio (see below) |
| SNP | bespoke: single effect-site compartment, continuous-infusion edge case | `CENT_SNP` | `C_SNP` is a direct alias (no volume division) |
| DOB | bespoke: single effect-site compartment, continuous-infusion edge case | `CENT_DOB` | `C_DOB` is a direct alias (no volume division); never dosed by any named scenario (see below) |

All eight are dosed as a **continuous zero-order input rate** added
directly inside `dxdt_GUT_<STEM>` (FUR/SAC/VAL/BB/MRA/SG, via
`rate_<stem>`) or directly into the effect-site compartment (SNP/DOB, via
`snp_inf`/`dob_inf`) — this is the guide's documented "continuous infusion
input, not an event" edge case, preserved exactly as the original does it.
The `rate_<stem>`/`snp_inf`/`dob_inf` dosing-input parameter names are
**not** part of the naming convention's table (they are external dosing
knobs, not PK/PD state) and were left unchanged.

## SAC/VAL: independent PK sharing no compartment

Sacubitril and valsartan are co-formulated clinically, but in this model
they have **fully independent PK**: separate depot/central compartments
(`GUT_SAC`/`CENT_SAC` vs. `GUT_VAL`/`CENT_VAL`), separate `KA`/`CL`/`V1`/`F`
parameters, and separate dosing-rate parameters (`rate_sac`, `rate_val`).
There is no shared dosing route or shared compartment to normalize — each
was refactored as its own independent archetype-3-minus-peripheral
compound, `C_SAC` and `C_VAL` each a clean single exposed concentration.

## FUR: two named effects, one Hill-shaped and one not

Furosemide's own block has, uniquely among the eight, **two structurally
different named effects** on the disease side, both functions of `C_FUR`
alone:

- `EFFECT_FUR` (natriuresis) — the original's `Emax_fur * Cfur / (EC50_fur
  * (1 + brake) + Cfur)` is a genuine Hill ratio, just with a
  *dynamically-braked* EC50 (the `(1 + brake)` term, `brake` being the
  model's own diuretic-braking state, now `BRAKE_FUR`). Renamed to
  `EMAX_FUR`/`EC50_FUR`/`GAMMA_FUR` (`GAMMA_FUR = 1`, rename not a fit) and
  written with the guide's `pow()` form; the `(1.0 + BRAKE_FUR)` multiplier
  on `EC50_FUR` is preserved unchanged, since it is FUR's own internal
  dynamic, not a second compound's interference.
- `EFFECT_FUR_RAAS` (RAAS-activation driver, feeds `Ang`'s target and
  `BRAKE_FUR`'s own dynamics) — the original's `furx = Cfur / EC50_fur`
  clamped at 1.5 is a **linear**, unsaturating ratio, not a Hill fraction
  (no `+ C` in the denominator), so it was renamed straight (`C_FUR /
  EC50_FUR`, same clamp) rather than forced into the `pow()` Hill template
  it structurally isn't.

`BRAKE_FUR`'s own `dxdt` uses `EFFECT_FUR_RAAS` (clamped at 3.0, same
redundant double-clamp as the original, since `EFFECT_FUR_RAAS` is already
≤ 1.5) — unchanged behaviour, renamed only.

## VAL: one vestigial, unused Emax preserved as-is

`Emax_val_ang` ("Valsartan maximal AT1 blockade") is declared in the
original's `$PARAM` but never read anywhere in `$ODE` — confirmed by
exact-name grep across the whole file. Renamed to `EMAX_VAL_ANG` and kept
in `$PARAM` (out of caution against deleting a parameter another part of
the fork's tooling might reference), but **no `EFFECT_VAL_ANG` was
invented** — there is nothing in the original to rename or redirect, and
inventing a new effect term the original never had would be scope creep,
not a rename. Only `EFFECT_VAL_RSYS` (SVR reduction) is real and named.

## DOB: no named scenario in the original

Dobutamine's effect-site compartment and both its named effects
(`EFFECT_DOB_EES`, inotropy; `EFFECT_DOB_HR`, chronotropy) are fully wired
into the disease side (`EesE`, `HR`), unlike VAL's vestigial parameter —
but **no scenario in the original file's own `## SCENARIOS` comment block
ever sets `dob_inf` away from its `0.0` default.** This is a documentation
gap in the original, not a code defect (dobutamine is real, clinically
appropriate acute inotropic support in decompensated MR, and its PK/PD is
fully implemented) — noted here rather than silently worked around. Since
there is no named scenario to reproduce verbatim, DOB was verified instead
via a direct `dob_inf` perturbation (see Verification below), structurally
identical to how the original's own Scenario 14 exercises `snp_inf`.

## SNP/DOB: bespoke effect-site archetype

Neither nitroprusside nor dobutamine has a depot/central PK split in the
original — each is a single first-order effect-site compartment driven
directly by a continuous infusion input already expressed in "effect
units" (`Ce_snp`/`Ce_dob`, now `CENT_SNP`/`CENT_DOB`), with no volume term
anywhere in its own dynamics (`dxdt_Ce_snp = ke_snp * (snp_inf - Ce_snp)`).
This does not fit archetypes 1–4 as written (no `CL`/`V1`, no depot) and
was handled as bespoke, per the guide's instruction for a genuine
non-archetype structure: renamed in place, `C_SNP`/`C_DOB` are **direct
aliases** of the renamed compartment (same pattern as `C_PHE`/`C_SNP`-style
already-concentration-valued states documented in this fork's
`pheochromocytoma`/`aortic-dissection` refactors), not divided by any
volume.

## Pre-existing build defect: the original does not compile under mrgsolve 2.0.1

`mr_mrgsolve_model.R` fails to compile as-written (confirmed via a live
`POST /model_manifest` call against the untouched original), for four
independent reasons, **none of which touch any of the eight compounds'
own PK/PD blocks** — all four are inside the shared haemodynamic/
circulation/remodelling machinery every compound's effect reads from or
writes into. Logged as
[`translations/UPSTREAM_ISSUES.md` #104](../translations/UPSTREAM_ISSUES.md):

1. Ten `$ODE`-local scratch variables (`Areq`, `dPmv`, `Pes`, `RVol`,
   `SVtot`, `Ped`, `MAP`, `vwave`, `Ppcw`, `CI`) shadow `$TABLE`'s own
   `double NAME = g_NAME;` declarations that populate the identically-named
   `$CAPTURE` outputs — mrgsolve compiles `$ODE` and `$TABLE` into the same
   scope, so this is a genuine C++ redefinition, not just a style clash.
2. `dP` is declared twice with `double` in sibling (not nested) `for`-loop
   bodies inside `$ODE`'s inner bisection — also a redefinition under this
   build, even though standard C++ block scoping would allow it.
3. The original's own `$PARAM EROApri_0` ("convenience copy of the initial
   orifice") collides with mrgsolve's auto-generated writable
   compartment-init symbol also named `EROApri_0` (auto-created because
   `$CMT` declares a compartment `EROApri`); `$MAIN`'s `EROApri_0 =
   EROApri_0;` is an unresolvable self-collision between the two meanings.
4. Four multi-variable `double a = x, b = y;` declarations (lines 569,
   570, 575, 584 of the original) lose their `double` keyword during this
   build, turning the whole statement into undeclared bare assignments.

**Fix applied directly to the delivered `_refactored.R`** (syntax-only,
non-numeric, per the guide's settled policy for a non-compiling original):
(1) renamed the ten shadowing `$ODE` locals with an `_l` suffix
(`Areq_l`, `dPmv_l`, `Pes_l`, `RVol_l`, `SVtot_l`, `Ped_l`, `MAP_l`,
`vwave_l`, `Ppcw_l`, `CI_l`) everywhere they are used internally — `$TABLE`
is completely untouched, so every `$CAPTURE` name and its meaning are
identical to the original; (2) renamed the second `dP` declaration and its
one downstream use to `dP2`; (3) renamed the colliding `$PARAM` to
`EROApri_init` and updated `$MAIN` to `EROApri_0 = EROApri_init;` (also
updated the one usage-example comment referencing the old parameter name,
for accuracy — a comment-only change); (4) split the four multi-variable
declarations into one `double X = value;` per variable, same initial
values. All four changes are mechanical renames/declaration-splits with
identical numerics — the verification below (exact match, max abs diff
0.0) is the proof that nothing behavioural changed.

## qspserver `/model_manifest` discoverability

Confirmed via a live `POST /model_manifest` call against the refactored
DSL: 184 parameters, 104 output paths. Every renamed PK parameter
(`KA_<STEM>`/`CL_<STEM>`/`V1_<STEM>`/`F_<STEM>` for the six depot+central
compounds, `KE_SNP`/`KE_DOB` for the two effect-site compounds) and every
renamed/new Hill parameter (`EMAX_<STEM>[_x]`/`EC50_<STEM>`/`GAMMA_<STEM>`)
is present in `parameters` with its original numeric default. `C_BB`,
`C_VAL`, `C_SAC`, `C_MRA`, `C_SG`, `C_FUR`, `C_SNP`, `C_DOB`, and all 16
`EFFECT_<STEM>[_x]` names are computed as plain `double` locals inside
`$ODE` (not `$PARAM` — a `$PARAM` declaration is read-only in mrgsolve
2.0.1 and `$ODE` cannot assign into it, the same pattern documented in
`ppgl_refactor_notes.md`) and are discoverable/retrievable instead via a
bare `$CAPTURE` addition; confirmed present in `outputPaths` and
retrievable through `/run_simulation`'s `outputs` selection (spot-checked
directly — see Verification).

## Verification

**Method.** Both DSLs — the original's own text with only the four
syntax-only build fixes above applied in-memory (never saved as a tracked
file; used solely as the "ground truth" reference since the true original
does not compile at all) and the delivered `_refactored.R`'s equivalent —
were POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`POST /model_manifest` then `POST
/run_simulation`), which compiles and runs each DSL block directly with
mrgsolve 2.0.1 server-side. Requests were spaced ~2s apart per the task's
concurrency note (`max_concurrent_jobs: 2`).

**Scenario A — all eight compounds simultaneously, from the healthy
baseline.** The original's own GDMT dosing list (`rate_fur=40, rate_bb=100,
rate_sac=194, rate_val=206, rate_mra=25, rate_sg=10`, from the file's own
`## SCENARIOS` comment) plus `snp_inf=4` (the exact value the original's
own Scenario 14 uses) plus `dob_inf=2` (DOB has no named scenario, see
above; chosen only to exercise its PK/PD alongside the other seven), run
from the default (healthy) initial state — `param()`-only, no `init()`
needed, over 100 days (`delta=5`, 21 points). Compared all 36 of the
original's own `$CAPTURE` outputs (unchanged names in both models):
**exact match, max absolute difference 0.0** at every output, every time
point.

**Scenario B — the original's own Scenario 14 verbatim** ("nitroprusside
in acute MR"): `init(EROApri=0.60)` (reproduced as a bolus dose of 0.60
into compartment index 10, `EROApri`'s position in `$CMT` — confirmed
identical in both models via `/model_manifest`'s `compartments` list) plus
`param(snp_inf=4)`, `end=2, delta=0.05` (42 points): **exact match, max
absolute difference 0.0** at every output, every time point. Values are
physiologically sensible and match the scenario's own narrative (e.g.
`vwave` rises from 0 toward the high-20s mmHg range as the acute orifice
opens onto the unadapted atrium).

**Discoverability spot-check.** A separate `/run_simulation` call
requesting `C_BB`, `C_FUR`, `C_SAC`, `C_VAL`, `C_MRA`, `C_SG`, `C_SNP`,
`C_DOB`, `EFFECT_BB_HR`, `EFFECT_FUR`, `EFFECT_SAC_FIB`,
`EFFECT_SNP_RSYS`, `EFFECT_DOB_EES` directly (refactored model only, GDMT +
snp_inf + dob_inf dosing) confirmed all of them are retrievable and
physiologically sensible (e.g. `C_SNP`/`C_DOB` track their infusion rates
exactly, `EFFECT_FUR` plateaus near 31 mL/day against its `EMAX_FUR=400`
ceiling, all Hill-shaped effects lie in [0,1)).

No tolerance was loosened anywhere in this pass — every comparison above
is exact (0.0), consistent with all eight compounds being pure structural
reorganizations (renames + archetype-conforming layout), never a
Hill-refit.
