# Refactor notes — `oa_mrgsolve_model_refactored.R`

Scope: **NSAID (celecoxib) and TANZ (tanezumab) only** — per the
`osteoarthritis | NSAID` and `osteoarthritis | TANZ` rows in
`driver-patches/data/compound_perturbation_census.md`, both classified
"Normalize duplicate concentration sites, then redirect." The file's other
three drugs (IA corticosteroid, IA hyaluronic acid, sprifermin) are not in
scope and are untouched (see diff-scope confirmation below), as is the
entire disease PD network (IL-1β, TNF-α, MMP-13, ADAMTS-5, Collagen II,
Aggrecan, chondrocyte, synovitis, osteoclast/osteoblast, JSW, PGE2, uCTX-II,
COMP) except for the two-token `EFFECT_TANZ`/`COX2_eff` reads it needed.

## What the duplicate-definition problem actually was

**NSAID.** Two concentrations (plasma and joint) were each computed once in
`$ODE` (`C_NSAID_plasma`, `C_NSAID_joint`) and then **independently
re-derived a second time in `$TABLE`**, under different names
(`C_nsaid_pl`, `C_nsaid_jt`), from the same raw compartments
(`A_NSAID_plasma`, `A_NSAID_joint`) divided by the same volume terms. The
`$ODE` copy fed the Hill/COX-2 term (`COX2_inh`); the `$TABLE` copy existed
only to be captured. Two names, two independent divisions, one physical
quantity, per tissue site.

**TANZ.** Same pattern, one site: `C_Tanz` computed once in `$ODE` (feeding
`Tanz_pain_red`) and re-derived independently in `$TABLE` as `C_Tanz_pl`
(`A_Tanz_plasma / Vd_tanz`, textually identical division, different name),
purely for capture.

## Redirect performed

**Within `$ODE`** (a single compiled scope): one cached `double C_NSAID =
PERI_NSAID / (V1_NSAID * 0.01)` (joint, PD-facing) and one cached `double
CP_NSAID = CENT_NSAID / V1_NSAID` (plasma, informational-only — exactly as
in the original, where `C_NSAID_plasma` was computed but never read
downstream either) are now the only sites that derive NSAID's
concentrations; `EFFECT_NSAID` reads `C_NSAID`. Likewise one cached `double
C_TANZ = CENT_TANZ / V1_TANZ` is the only site for tanezumab; `EFFECT_TANZ`
reads it.

**Between `$ODE` and `$TABLE`:** confirmed empirically (via the qspserver
`mrgsolve_api` container, same method documented in
`alkaptonuria/aku_refactor_notes.md` and
`kidney-transplant-rejection/ktx_refactor_notes.md`) that `$ODE` and
`$TABLE` share enough of mrgsolve's generated scope that **re-declaring**
the same name in both blocks is a hard compile error, but *reading* an
already-`$ODE`-declared local from `$TABLE` without redeclaring it is not —
this file's own original code already relied on exactly that (`$TABLE`'s
`COX2_pct_inh = COX2_inh * NSAID_flag * 100.0` reads `COX2_inh`, an `$ODE`
local, without ever redeclaring it). The refactor keeps that working
pattern: `$TABLE`'s `COX2_pct_inh` now reads `EFFECT_NSAID` (the renamed
`$ODE` local) directly, no second Hill computation.

For the two purely-diagnostic concentrations that `$TABLE` captures under
their own distinct, already-non-colliding names (`C_nsaid_pl`, `C_nsaid_jt`,
`C_Tanz_pl`), the safer, previously-verified-safe pattern from
`aku_refactor_notes.md` was followed instead of relying on further
cross-block reads: kept these three names as `$TABLE`-local (no rename
needed, they never collided with the new `$ODE` names `C_NSAID`/`CP_NSAID`/
`C_TANZ`), but **redirected their right-hand side** to the renamed
canonical compartments/parameters (`CENT_NSAID`/`V1_NSAID`,
`PERI_NSAID`/`V1_NSAID`, `CENT_TANZ`/`V1_TANZ` — was `A_NSAID_plasma`/
`Vd_nsaid`, `A_NSAID_joint`/`Vd_nsaid`, `A_Tanz_plasma`/`Vd_tanz`) instead
of the old names. The duplication across the `$ODE`/`$TABLE` boundary is an
unavoidable mrgsolve architecture fact for these two names (each is a
distinct evaluation, mathematically identical to its `$ODE` sibling by
construction), not something left inconsistent any more — both sides now
read from the same renamed compartments and parameters instead of two
independently-maintained old names.

## Naming convention applied

### NSAID — archetype 3 minus peripheral-via-Q/V2 (depot + central + a
second, partition-linked tissue compartment)

The original's joint compartment is not the guide's Q/V2 peripheral
archetype (no back-diffusion term `Q*(C1-C2)`); it is a one-way partition
off central clearance (`Kp_nsaid_jt * (CL/V1) * CENT`) that itself
eliminates at the same rate constant. Kept this bespoke structure exactly
as the original built it — only renamed, no compartments added or removed,
no algebra changed.

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `A_NSAID_gut` (cmt) | `GUT_NSAID` | `ka_nsaid` | `KA_NSAID` |
| `A_NSAID_plasma` (cmt) | `CENT_NSAID` | `F_nsaid` | `F_NSAID` |
| `A_NSAID_joint` (cmt) | `PERI_NSAID` | `Vd_nsaid` | `V1_NSAID` |
| `C_NSAID_plasma` (`$ODE` local) | `CP_NSAID` | `CL_nsaid` | `CL_NSAID` |
| `C_NSAID_joint` (`$ODE` local) | `C_NSAID` | `Kp_nsaid_jt` | `KP_NSAID` (bespoke, no convention slot) |
| `C_nsaid_pl`/`C_nsaid_jt` (`$TABLE` locals) | kept, RHS redirected | `IC50_cox2` | `EC50_NSAID` |
| `COX2_inh` (`$ODE` local) | `EFFECT_NSAID` | `Emax_cox2` | `EMAX_NSAID` |
| | | `hill_cox2` | `GAMMA_NSAID` |

`EFFECT_NSAID = EMAX_NSAID*C_NSAID^GAMMA_NSAID /
(EC50_NSAID^GAMMA_NSAID+C_NSAID^GAMMA_NSAID)` is a pure rename of
`COX2_inh` — already exactly this shape, `GAMMA_NSAID` is the original's
own fitted `hill_cox2 = 1.2`. `NSAID_flag`'s gating stays exactly where the
original had it (multiplying `EFFECT_NSAID` into `COX2_eff`, not baked into
`EFFECT_NSAID` itself) — preserved, not moved, per the guide's "keep a
calculation in the block/form the original used it in."

### TANZ — archetype 3 minus peripheral (depot + central, linear); checked
for TMDD, none present

Tanezumab is an anti-NGF monoclonal antibody, a molecule class that often
needs target-mediated drug disposition. **This file does not model TMDD**:
there is no receptor-binding ODE state anywhere for tanezumab — just a
depot→central linear PK (`dxdt_A_Tanz_depot`/`dxdt_A_Tanz_plasma`, no
peripheral compartment) and a direct `Emax*C/(EC50+C)` pain-reduction
ratio. Confirmed by reading the full `$ODE` block; no receptor pool, no
`kon`/`koff`, no free/bound antibody split. Archetype 3 minus peripheral,
same class as the NSAID/SSZ/TNFI compounds in
`reactive-arthritis/rea_refactor_notes.md`.

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `A_Tanz_depot` (cmt) | `GUT_TANZ` | `ka_tanz` | `KA_TANZ` |
| `A_Tanz_plasma` (cmt) | `CENT_TANZ` | `F_tanz` | `F_TANZ` |
| `C_Tanz` (`$ODE` local) | `C_TANZ` | `Vd_tanz` | `V1_TANZ` |
| `C_Tanz_pl` (`$TABLE` local) | kept, RHS redirected | `CL_tanz` | `CL_TANZ` |
| `Tanz_pain_red` (`$ODE` local) | `EFFECT_TANZ` | `Emax_tanz_pain` | `EMAX_TANZ` |
| `Tanz_flag` | `TANZ_flag` | `K_tanz_NGF` | `EC50_TANZ` |
| | | (none — implicit) | `GAMMA_TANZ = 1` |

`EFFECT_TANZ = EMAX_TANZ*C_TANZ^GAMMA_TANZ /
(EC50_TANZ^GAMMA_TANZ+C_TANZ^GAMMA_TANZ) * TANZ_flag` is a pure rename of
`Tanz_pain_red` — algebraically identical to the original's
`Emax*C/(EC50+C)` ratio with `GAMMA_TANZ=1` (the original had no explicit
Hill exponent). The flag multiplication stays inline in the effect term
itself, exactly as the original had it (unlike NSAID's flag, which the
original applied one step downstream) — each compound's own original
placement was preserved, not normalized to match the other.

All PK for both compounds is a pure rename — no value changed, no
compartment added, removed, or reordered (compartment 1-based positions are
unchanged: `GUT_NSAID`/`CENT_NSAID`/`PERI_NSAID` are still positions 1-3,
`GUT_TANZ`/`CENT_TANZ` still 8-9, exactly where `A_NSAID_gut` etc. and
`A_Tanz_depot`/`A_Tanz_plasma` were).

## Build-compatibility fixes (not NSAID/TANZ-scoped)

The untouched original does not compile under mrgsolve 2.0.1 at all, for
two reasons unrelated to either compound's own PK/PD block. Logged as
`translations/UPSTREAM_ISSUES.md` entry **#121**; syntax-only, non-numeric
fixes applied directly to the delivered `_refactored.R`/`.cpp`, per the
guide's settled policy:

1. **`$CAPTURE` re-lists all fifteen disease-state `$INIT` compartment
   names** (`IL1b TNFa MMP13 ADAM5 ColII Aggrecan Chondro Synovitis OC_act
   OB_act JSW PGE2_jt VASPain uCTXII COMP_s`), illegal under mrgsolve
   2.0.1. Fixed by dropping all fifteen from `$CAPTURE` — every one is
   still emitted automatically as a compartment output column, confirmed
   present under its same name in every verification run below.
2. **Three `$PARAM` names (`ColII_0`, `Chondro_0`, `JSW_0`) collide with
   mrgsolve's auto-generated `<compartment>_0` init symbol** for the
   identically-named compartments. Fixed by renaming to `ColII_BASE`/
   `Chondro_BASE`/`JSW_BASE` and updating every reference (`col_syn`,
   `agg_syn`, `Struct_pain`, `KOOS_func`, `CartVol_pct`, and the R
   wrapper's summary `cat()` line). `Agg_0` does not collide (the
   compartment is `Aggrecan`, not `Agg`) and was left unchanged.

Neither fix touches anything inside the NSAID or TANZ PK/PD blocks.

## qspserver compatibility

Confirmed via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`): the refactored DSL compiles and its
manifest lists `KA_NSAID, F_NSAID, V1_NSAID, CL_NSAID, KP_NSAID,
EC50_NSAID, EMAX_NSAID, GAMMA_NSAID` and `KA_TANZ, F_TANZ, V1_TANZ,
CL_TANZ, EMAX_TANZ, EC50_TANZ, GAMMA_TANZ` as discoverable `$PARAM`
entries (none of the old names remain), and its `outputPaths` list shows
`GUT_NSAID, CENT_NSAID, PERI_NSAID, GUT_TANZ, CENT_TANZ` in the same
1-based compartment positions the original's `A_NSAID_gut, A_NSAID_plasma,
A_NSAID_joint, A_Tanz_depot, A_Tanz_plasma` held — dosing into compartment
1 (NSAID) and compartment 8 (TANZ) is therefore unaffected. `C_NSAID`/
`CP_NSAID`/`C_TANZ`/`EFFECT_NSAID`/`EFFECT_TANZ` could not be added to
`$PARAM` themselves (they are computed `$ODE` locals, not overridable —
the same read-only-reference constraint documented in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md` and
`alkaptonuria/aku_refactor_notes.md`); their `$TABLE`/`$CAPTURE`-scope
values are discoverable via `outputPaths` as `C_nsaid_pl`, `C_nsaid_jt`,
`C_Tanz_pl`, and `COX2_pct_inh`.

A bare-DSL `.cpp` extraction of the `oa_code <- '...'` quoted block (a
straight text pull, confirmed byte-identical to the DSL embedded in
`_refactored.R`) was used as `model_content` for every `/model_manifest`
and `/run_simulation` call below; it was a scratch verification artifact
only, not committed to the repo, matching the rest of this batch's
practice.

## Verification

**Method.** Both DSL blocks (untouched-original-plus-build-compat-fix, and
the full refactored sibling) were extracted verbatim and POSTed to the
qspserver `mrgsolve_api` container's `/run_simulation`, spaced ~2.5 s
apart, reproducing the original file's own dosing scenarios exactly (same
`amt`/`ii`/`addl`/`time` as `nsaid_regimen`, `tanz_regimen`,
`iacs_regimen`, `sprif_regimen`, and each scenario's own `params_override`
flags), with a `730`-day/`delta=1` window matching the original's own
`sim_end`. No step-count shortening was needed — every run below completed
within the API's default solver budget at the file's full 2-year horizon.

**Result: exact match on every scenario tested.**

| Scenario (original's own) | Compound(s) exercised | Points compared | Max abs diff | Max rel diff |
|---|---|---|---|---|
| 2. Celecoxib 200 mg BID | NSAID alone | 732 | 0.0 | 0.0 |
| 6. Tanezumab 2.5 mg SC Q8W | TANZ alone | 745 | 0.0 | 0.0 |
| 7. Celecoxib + IA Triamcinolone | NSAID + untouched IACS | 740 | 0.0 | 0.0 |
| 8. Sprifermin + Celecoxib | NSAID + untouched Sprifermin | 741 | 0.0 | 0.0 |

Compared across all shared `$CAPTURE`d/compartment outputs: `IL1b, TNFa,
MMP13, ADAM5, ColII, Aggrecan, Chondro, Synovitis, OC_act, OB_act, JSW,
PGE2_jt, VASPain, uCTXII, COMP_s, COX2_pct_inh, C_nsaid_pl, C_nsaid_jt,
C_IACS_jt, C_HA_jt, C_Sprif_jt, C_Tanz_pl, GR_inhibition_pct, KOOS_est,
WOMAC_pain_sub, CartVol_pct, ECM_degradation`, and every drug-PK
compartment (`GUT_NSAID`/`A_NSAID_gut`, `CENT_NSAID`/`A_NSAID_plasma`,
`PERI_NSAID`/`A_NSAID_joint`, `GUT_TANZ`/`A_Tanz_depot`,
`CENT_TANZ`/`A_Tanz_plasma`, plus the three untouched IACS/HA/Sprifermin
compartments). This is the exact match expected for a pure structural
reorganization (archetype rename only, no Hill-fitting, no behavioral
change) per the guide's tolerance rule for archetypes 1–3.

## Anything else worth flagging

- `KL_mult` (declared `= 1.0 // Runtime calculated from KL_grade`) is dead
  in the original — never referenced in `$ODE`/`$TABLE` despite the
  comment implying it should modulate disease load. Not touched (outside
  NSAID/TANZ scope, and fixing it would change simulated behavior).
- The 15 disease-state compartments removed from `$CAPTURE` (build-compat
  fix #121) remain fully present in every simulation output automatically
  — confirmed identical values, same names, in every verification run
  above.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`osteoarthritis | NSAID` and `osteoarthritis | TANZ` rows.
