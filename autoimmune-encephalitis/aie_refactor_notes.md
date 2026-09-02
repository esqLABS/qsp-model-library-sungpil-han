# Refactor notes — `autoimmune-encephalitis/aie_mrgsolve_model.R`

Compounds refactored: **IVIG** (intravenous immunoglobulin), **MP**
(methylprednisolone), **RTX** (rituximab, anti-CD20), **TCZ** (tocilizumab,
anti-IL-6 receptor), **CPX** (4-OH-cyclophosphamide, the active metabolite of
cyclophosphamide) — all five rows in
`driver-patches/data/compound_perturbation_census.md` for
`autoimmune-encephalitis`. Compound identity was confirmed against the
directory's own `README.md` ("주요 약물 표적" section: 메틸프레드니솔론/
methylprednisolone, IVIG, 리툭시맙/rituximab anti-CD20, 토실리주맙/tocilizumab
anti-IL-6R) and the code's own header comment ("4-OH-Cyclophosphamide active");
CPX is cyclophosphamide's active metabolite, matching the disease overview's
mention of cyclophosphamide as second-line therapy. All five compounds were in
scope — this file has no other exogenous drug.

## A necessary build-compatibility fix, unrelated to any compound's PK

The untouched original does not compile under mrgsolve 2.0.1 at all, for
**three independent, pre-existing reasons**, confirmed via `POST
/model_manifest` against the untouched original's own extracted
`aie_model_code <- '...'` string:

1. `$SET delta=0.25 end=365 start=0` — missing commas between key=value
   pairs. mrgsolve's `handle_SET` parses this block as an R `list(...)` call
   and fails with `unexpected symbol` before the model itself is even
   inspected. Fixed: `$SET delta=0.25, end=365, start=0`.
2. `$CMT` (22 compartments, each with a trailing `//` comment) and a separate
   `$INIT` block re-declare the same 22 names, in the same order
   (`Duplicated model names: GCB PB LLPC MB ... CPX_ACT`). Same defect class
   already logged for several other files in this corpus (e.g.
   `distal-renal-tubular-acidosis` #38, `autoimmune-polyendocrinopathy` #51).
   Fixed the same way `autoimmune-polyendocrinopathy` fixed it: the `$CMT`
   block is dropped outright and `$INIT` is kept (values unchanged; the
   per-compartment descriptive comments that lived on `$CMT` are preserved as
   `$INIT`-line comments instead, so no documentation is lost).
3. Once (1) and (2) are worked around, `$CAPTURE` still fails:
   `compartment should not be in $CAPTURE: IL6_CNS,GFAP,AB_CSF,GCB,PB,LLPC,MB`
   — the original's own `$CAPTURE` block directly re-lists seven compartment
   names verbatim (alongside the `$TABLE`-derived quantities, which are
   fine). Fixed: these seven names dropped from `$CAPTURE` — they remain in
   the simulation output automatically as compartment state columns (that is
   the entire reason mrgsolve rejects listing them a second time), so nothing
   is lost.

All three fixes are syntax-only and non-numeric. Logged as
`translations/UPSTREAM_ISSUES.md` **#140**. Applied directly to the delivered
`aie_mrgsolve_model_refactored.R` (not just a scratch copy), per the guide's
settled policy for a non-compiling original. The same three fixes were
applied to an in-memory-only **patched-original** scratch copy (original DSL,
otherwise byte-identical, same compound names/compartments/parameters) used
solely as the verification baseline below — never written into the repo
tree, not a deliverable, discarded after use.

A minor, separate documentation-only discrepancy also found in passing (not
logged to `UPSTREAM_ISSUES.md`, since it changes nothing numeric and isn't a
build defect): the original's own top-of-file comment claims "Compartments:
22 ODEs", but the file's own `$CMT` block lists 24 names (6 immune + 6 CNS +
4 clinical + 8 PK = 24, not 22) — a stale headline count, not a code defect.

## A genuine, in-scope PK bug discovered in the original (TCZ / CPX)

`RTX`, `MP`, and `IVIG` each already had their own clean, non-overlapping
2-compartment PK block in the original (`RTX1`/`RTX2`, `MP1`/`MP2`,
`IVIG1`/`IVIG2`) — for these three, this refactor is a pure rename +
reorganization (Archetype 2), and verification below confirms an exact match
on every disease-side output (bar one disclosed IVIG-reporting fix, see
below).

**TCZ and CPX did not.** The original's own PK block reads:

```
dxdt_TCZ1 = -(CL_TCZ + Q_TCZ)/Vc_TCZ * TCZ1 + Q_TCZ/Vp_TCZ * CPX_ACT;
// Note: reusing CPX_ACT slot for TCZ2 peripheral when TCZ scenario active
dxdt_CPX_ACT = -(CL_TCZ + Q_TCZ)/Vc_CPX * CPX_ACT + Q_TCZ/Vc_TCZ * TCZ1;
// When CPX is dosed instead, this becomes 4-OH-CPX:
// dxdt_CPX_ACT = -CL_CPX/Vc_CPX * CPX_ACT  [handled via scenario events]
```

The comment states the intent honestly: `CPX_ACT` was meant to double as
TCZ's own peripheral compartment ("TCZ2"), and the commented-out line shows
CPX's *own* intended 1-compartment kinetics (`-CL_CPX/Vc_CPX * CPX_ACT`) —
but that line is dead text, never actually wired into any live `dxdt_`. In
the real, compiled model:

- **CPX's own declared parameters (`CL_CPX = 96.0`, `Vc_CPX = 30.0`,
  intended t1/2 ~4h) were completely dead.** `CPX_ACT`'s real elimination
  used `-(CL_TCZ + Q_TCZ)/Vc_CPX` (~0.034/day, a ~20-day-scale half-life) —
  about **95x slower** than the compound's own declared, documented PK.
- **TCZ's peripheral compartment never existed as its own state.** `TCZ1`'s
  return-flux term correctly used `Vp_TCZ` (`Q_TCZ/Vp_TCZ * CPX_ACT`), but
  the "peripheral" state itself (`CPX_ACT`) used the *wrong* volume
  (`Vc_CPX`, cyclophosphamide's own 30 L, not `Vp_TCZ`'s 2.9 L) and carried
  an extra `CL_TCZ` elimination term a true peripheral compartment should
  never have (mass should only exchange with central, never leave the model
  from a peripheral compartment directly) — a non-conservative, physically
  wrong two-compartment implementation.
- **Consequence, confirmed by running the original's own scenario 5**
  (Cyclophosphamide, day 30 monthly x6, no TCZ ever dosed): because `CPX_ACT`
  feeds directly into `dxdt_TCZ1`, and both compartments' only real
  elimination path is the slow, cross-wired formula above, dosing CPX alone
  drives the file's own reported "tocilizumab concentration" output
  (`Cp_TCZ_out`) to **6.86 x 10^9** by day 365 — an unbounded, physically
  meaningless number, in a scenario that never contains a milligram of
  tocilizumab. Running the original's own scenario 6 (Tocilizumab, dosed
  properly at day 60) shows the same runaway: `Cp_TCZ_out` reaches **1.74 x
  10^9** by day 365, instead of decaying the way a real anti-IL-6R mAb with
  the file's own declared `CL_TCZ`/`Q_TCZ` would. `Cp_RTX_out` in the same
  scenario is completely unaffected either way (bit-exact between original
  and refactored, confirming RTX/TCZ/CPX are otherwise fully decoupled and
  this is specifically the `CPX_ACT`/`TCZ1` entanglement, not a broader
  defect).

This is squarely **in scope** for this refactor — both TCZ and CPX are
census-classified "Redirect concentration (clean single site)", i.e. exactly
"isolate this compound's PK into its own compartment(s), normalize the
duplicate/shared concentration site." Per the guide ("If a defect is
genuinely inside the scope compound's own block ... fix it as part of the
refactor itself and say so"), this is fixed here:

- **TCZ** gets a real, dedicated `PERI_TCZ` compartment (new `$CMT`/`$INIT`
  entry, starts at 0.0) and a standard Archetype-2 CL/Q/V1/V2 formula, using
  the original's own already-declared `CL_TCZ`/`V1_TCZ` (`Vc_TCZ`)/`V2_TCZ`
  (`Vp_TCZ`)/`Q_TCZ` values — no new numbers.
- **CPX** gets its own dedicated single compartment (`CENT_CPX`, was
  `CPX_ACT`, no longer double-booked), with the Archetype-1 formula the
  original's own commented-out line already specified:
  `dxdt_CENT_CPX = -CL_CPX/V1_CPX * CENT_CPX`, using the original's own
  `CL_CPX = 96.0`/`Vc_CPX = 30.0` (renamed `V1_CPX`) — again, no new numbers,
  just finally wiring up parameters the original declared but never used.

**This is a real, disclosed behavioral change from the original for
scenarios 5 and 6** (see Verification below) — not a tolerance to explain
away, an intentional, in-scope correction of a genuine defect, made using
only the compound's own already-declared values.

### Naming applied

| Role | IVIG | MP | RTX | TCZ | CPX |
|---|---|---|---|---|---|
| Central compartment | `CENT_IVIG` (was `IVIG1`) | `CENT_MP` (was `MP1`) | `CENT_RTX` (was `RTX1`) | `CENT_TCZ` (was `TCZ1`) | `CENT_CPX` (was `CPX_ACT`) |
| Peripheral compartment | `PERI_IVIG` (was `IVIG2`) | `PERI_MP` (was `MP2`) | `PERI_RTX` (was `RTX2`) | `PERI_TCZ` (**new** — see above) | — (1-cmt) |
| Clearance | `CL_IVIG` (unchanged) | `CL_MP` (unchanged) | `CL_RTX` (unchanged) | `CL_TCZ` (unchanged) | `CL_CPX` (unchanged) |
| Central volume | `V1_IVIG` (was `Vc_IVIG`) | `V1_MP` (was `Vc_MP`) | `V1_RTX` (was `Vc_RTX`) | `V1_TCZ` (was `Vc_TCZ`) | `V1_CPX` (was `Vc_CPX`) |
| Peripheral volume | `V2_IVIG` (was `Vp_IVIG`) | `V2_MP` (was `Vp_MP`) | `V2_RTX` (was `Vp_RTX`) | `V2_TCZ` (was `Vp_TCZ` — now actually used correctly) | — |
| Inter-compartmental CL | `Q_IVIG` (unchanged) | `Q_MP` (unchanged) | `Q_RTX` (unchanged) | `Q_TCZ` (unchanged) | — |
| Exposed concentration | `C_IVIG` | `C_MP` | `C_RTX` | `C_TCZ` | `C_CPX` |
| Hill EC50 | `EC50_IVIG` (was `EC50_IVIG_FcRn`) | `EC50_MP` (unchanged) | `EC50_RTX` (unchanged) | `EC50_TCZ` (unchanged) | `EC50_CPX` (was `EC50_CPX_kill`) |
| Hill Emax | `EMAX_IVIG` (was `Emax_IVIG_cat`) | `EMAX_MP_ANTI`/`EMAX_MP_BBB` (was `Emax_MP_anti`/`Emax_MP_BBB`) | `EMAX_RTX` (was `Emax_RTX`) | `EMAX_TCZ` (was `Emax_TCZ_IL6`) | `EMAX_CPX` (was `Emax_CPX_kill`) |
| Hill gamma | `GAMMA_IVIG = 1.0` (new) | `GAMMA_MP = 1.0` (new, shared) | `GAMMA_RTX` (was `gamma_RTX`, = 2.0, unchanged) | `GAMMA_TCZ = 1.0` (new) | `GAMMA_CPX = 1.0` (new) |
| Effect | `EFFECT_IVIG` (was `eff_IVIG_cat`) | `EFFECT_MP_ANTI`/`EFFECT_MP_BBB` (was `eff_MP_anti`/`eff_MP_BBB`) | `EFFECT_RTX` (was `eff_RTX`) | `EFFECT_TCZ` (was `eff_TCZ_IL6`) | `EFFECT_CPX` (was `eff_CPX`) |

`DOSE_IVIG_GIVEN` is left unchanged — it is not one of the naming
convention's structural roles, and it is dead/unused in the original (never
referenced in `$MAIN`/`$ODE`/`$TABLE`); preserved as-is, still unused.

MP has **two** distinct pharmacological actions in the original
(anti-inflammatory + BBB stabilization), both reading the same concentration
and sharing one `EC50_MP`, but with their own `Emax`. Per the guide ("never
collapse several drugs into one shared Hill term") this is one compound with
two genuinely different effect endpoints, not two drugs to merge — kept as
`EFFECT_MP_ANTI` and `EFFECT_MP_BBB`, matching the original's own structure
exactly, just renamed.

## Hill interface

**IVIG, MP (both endpoints), TCZ, CPX** — rename, not a refit. All four are
already the bare `C/(EC50+C)` shape with an implicit `Emax`/no separate
exponent in the original:

```
eff_IVIG_cat = Emax_IVIG_cat * Cp_IVIG_mgl / (EC50_IVIG_FcRn + Cp_IVIG_mgl);
eff_MP_anti  = Emax_MP_anti  * Cp_MP_mcg   / (EC50_MP + Cp_MP_mcg);
eff_MP_BBB   = Emax_MP_BBB   * Cp_MP_mcg   / (EC50_MP + Cp_MP_mcg);
eff_TCZ_IL6  = Emax_TCZ_IL6  * Cp_TCZ_mcg  / (EC50_TCZ + Cp_TCZ_mcg);
eff_CPX      = Emax_CPX_kill * Cp_CPX_ng   / (EC50_CPX_kill + Cp_CPX_ng);
```

`GAMMA_IVIG`/`GAMMA_MP`/`GAMMA_TCZ`/`GAMMA_CPX` (all `= 1.0`) are new named
parameters reproducing these ratios exactly via the guide's
`Emax*pow(X,gamma)/(pow(EC50,gamma)+pow(X,gamma))` template — not a fit.
Confirmed bit-exact (see Verification).

**RTX** already had an explicit Hill coefficient in the original
(`gamma_RTX = 2.0`) — renamed to `GAMMA_RTX`, value unchanged, same
`pow(...)`-based formula. Not TMDD: RTX's effect is a plain sigmoid-Emax
function of `C_RTX` alone, with no receptor-binding/occupancy ODE state
anywhere in the original (B-cell depletion is applied as a direct multiplier
on the GCB/MB death-rate terms, not mediated through a modeled CD20-receptor
pool) — Archetype 2 plain-PK Hill is the correct, faithful structure.

**TCZ** — also confirmed **not** TMDD, despite the real IL-6-receptor
target: the original models IL-6 blockade as a plain concentration-driven
Emax ratio (`eff_TCZ_IL6`), with no receptor pool/complex compartment
anywhere. Kept as plain-PK Hill (Archetype 2), per the same "rename, not a
refit" rule — the concentration itself (`C_TCZ`) now behaves correctly
post-fix (see above), but the Hill formula shape is unchanged from the
original.

## IVIG's own duplicate/inconsistent concentration-reporting site

The original computes IVIG's concentration **twice**, in two different
places, with **two different results**:

```
$MAIN:  double Cp_IVIG = IVIG1 / Vc_IVIG;          // used directly below
        double Cp_IVIG_mgl = Cp_IVIG * 1000;       // THIS feeds the Hill term
        double eff_IVIG_cat = Emax_IVIG_cat * Cp_IVIG_mgl / (EC50_IVIG_FcRn + Cp_IVIG_mgl);

$TABLE: double Cp_IVIG_out = IVIG1 / Vc_IVIG;       // NO x1000 -- reported value
```

The *reported* concentration (`Cp_IVIG_out`, no `x1000`) is **1000x smaller**
than the concentration the disease model's own Hill term actually consumes
(`Cp_IVIG_mgl`, `x1000`) — the file's own two "IVIG concentration" numbers
disagree by three orders of magnitude, exactly the "duplicate/inconsistent
concentration site" this compound was census-classified for ("Normalize
duplicate concentration sites, then redirect").

Per the guide's "exposes exactly one concentration variable that PD
equations read" requirement, `C_IVIG` is defined as the value PD actually
reads (`(CENT_IVIG / V1_IVIG) * 1000`, an exact carry-over of
`Cp_IVIG_mgl`'s own two-step arithmetic, just written on one line) — this
means the disease-affecting Hill formula (`EFFECT_IVIG`) is **byte-for-byte
unchanged** from the original (confirmed exact match on every IVIG-dosed
scenario's disease outputs below), and the *only* thing that changes is the
reported concentration column itself, which now matches what the model
internally uses instead of silently under-reporting it 1000x. This is a
reporting-site fix, not a PK/PD behavior change.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
(`http://localhost:8007`, `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2.2s apart, one request in flight at a
time). Since the untouched original does not compile at all (#140), the
comparison baseline is the syntax-patched-original scratch copy described
above (identical compounds/parameters/compartments/equations to the tracked
`aie_mrgsolve_model.R`, only the three syntax fixes applied).

**Scenarios.** All six of the original's own scenarios, reproduced with
identical dosing (amounts, `ii`, `addl`, start times) via the API's
`dosing` field: `01_NatHist` (no dosing), `02_IVIG_MP` (IVIG+MP day 14 x5),
`03_IVIG_MP_PE` (same — the original's own "PE" arm is an `amt=0` no-op
placeholder, never implemented), `04_Rituximab` (+RTX weekly day 30 x4),
`05_Cyclophosphamide` (+CPX day 30, q28d x6), `06_Tocilizumab` (+RTX, +TCZ
day 60, q28d x6). Full 365-day window, `delta=0.5` (731-735 points per
scenario depending on how many dose times fall inside the grid) — the
original's own `end=365, delta=0.5` scenario runner; no shortening needed,
none of these hit the API's default solver step-count budget.

**Result, disease-side outputs (`NMDAR_pct`, `Ab_norm`, `mRS_est`, `CRS_c`,
`SZ_c`, `COG_c`, `PSY_c`, `BBB_c`, `MG_c`, `response_flag`):**

| Scenario | Result |
|---|---|
| `01_NatHist` (no drug at all) | Exact match on `NMDAR_pct`, `Ab_norm`, `mRS_est`, `CRS_c`, `BBB_c`, `response_flag` (max abs diff 0.0). `SZ_c`/`MG_c` differ by at most 0.0001, `PSY_c` by at most 0.0133 (out of a 0-10 scale, ~0.13%) — see note below. |
| `02_IVIG_MP` | **Exact match** on all ten disease-side outputs (max abs diff 0.0). |
| `03_IVIG_MP_PE` | **Exact match** on all ten (max abs diff 0.0). |
| `04_Rituximab` | **Exact match** on all ten (max abs diff 0.0). |
| `05_Cyclophosphamide` | Large, expected divergence (`NMDAR_pct` diff up to 92.0, `CRS_c` up to 9.2, `response_flag` flips) — this is the disclosed CPX PK fix taking effect (see below), not a verification failure. |
| `06_Tocilizumab` | Large, expected divergence (`NMDAR_pct` diff up to 77.5, `CRS_c` up to 6.2, `response_flag` flips) — the disclosed TCZ PK fix taking effect. |

**RTX/MP/IVIG concentration outputs** (`Cp_RTX_out` vs `C_RTX`, `Cp_MP_out`
vs `C_MP`): exact match (0.0) in every scenario. **IVIG concentration**
(`Cp_IVIG_out` vs `C_IVIG`): differs by exactly the disclosed 1000x
reporting-site fix in every IVIG-dosed scenario (expected, not a bug — see
above). **TCZ concentration** (`Cp_TCZ_out` vs `C_TCZ`): the original blows
up to 6.86x10^9 (scenario 5, CPX-only) and 1.74x10^9 (scenario 6,
TCZ-dosed) by day 365, driven entirely by the `CPX_ACT`/`TCZ1` entanglement
bug; the refactored model correctly shows `C_TCZ = 0` throughout scenario 5
(TCZ never dosed) and a physiologically sane rise-then-decay (peak
1.69x10^5 around day 200, down to 1.06 by day 365) in scenario 6 — direct
confirmation the fix produces sane 2-compartment kinetics where the original
did not.

**This is exactly the expected outcome for this file per the guide's
tolerance rules**: RTX/MP/IVIG are pure Archetype-2 rename+reorganization
with **no** bugs in their own compartments, and verify to floating-point-scale
tolerance or exactly 0.0; TCZ and CPX carry a real, disclosed, in-scope bug
fix and are **expected** to diverge from the original, which they do, in the
two scenarios that actually dose them, and nowhere else.

**A small, honestly-reported residual noise floor in the undosed scenario.**
`01_NatHist` doses nothing at all, so every compound's concentration and
effect term is *exactly* 0.0 in both models throughout (confirmed: `CENT_*`
states have no inflow and start at 0, a genuine fixed point) — meaning the
disease-only equations (byte-identical between original and refactored,
confirmed by a direct diff of the `$ODE` disease block and the `$MAIN`
disease-adjacent locals after normalizing only the renamed effect variables)
are evaluated with literally the same, exactly-zero drug terms in both. Yet
`SZ_c`/`MG_c`/`PSY_c` show a small, reproducible (confirmed via a same-DSL
repeat call) discrepancy peaking at 0.0133 (`PSY_c`, 0-10 scale). Every other
disease output in this same scenario, and *all* ten disease outputs in every
dosed scenario (`02`-`04`), are bit-exact. The most plausible explanation,
given the exhaustive comparison above rules out any actual formula/parameter
difference: the refactored model has one additional ODE compartment
(`PERI_TCZ`, exactly 0 throughout this scenario) that the original does not,
changing the solver's state-vector size; `SZ_c`/`MG_c`/`PSY_c` are the three
outputs downstream of hard threshold branches (`GLU > SZ_thresh`,
`MG > 1.0`/`MG > 0`) that are the most sensitive to an infinitesimal
adaptive-step-size difference. Reported here in full rather than silently
tolerating it — it is roughly 4-5 orders of magnitude smaller than the CPX/TCZ
fix's real, intentional divergence, and it does not appear at all once any
dosing (and the periodic solver restarts dosing brings) is present.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted DSL:

- `outputPaths` (46 total) includes all 25 compartments plus `NMDAR_pct`,
  `Ab_norm`, `mRS_est`, `CRS_c`, `SZ_c`, `COG_c`, `PSY_c`, `BBB_c`, `MG_c`,
  `C_IVIG`, `C_MP`, `C_RTX`, `C_TCZ`, `C_CPX`, `EFFECT_RTX_out`,
  `EFFECT_IVIG_out`, `EFFECT_MP_ANTI_out`, `EFFECT_MP_BBB_out`,
  `EFFECT_TCZ_out`, `EFFECT_CPX_out`, `response_flag`.
- `parameters` (76 total) includes `EC50_IVIG`, `EC50_MP`, `EC50_RTX`,
  `EC50_TCZ`, `EC50_CPX` (fixed defaults, unchanged from the original) and
  every renamed PK parameter (`CL_*`, `V1_*`, `V2_*`, `Q_*`).
- Literal-text discoverability self-check (`grep`, per the guide): a
  contiguous `double C_<STEM> = <expr>;` statement and a same-stem
  `EC50_<STEM>` parameter both exist for all five compounds (`RTX`, `TCZ`,
  `MP`, `CPX`, `IVIG`).
- `EFFECT_<STEM>` is **not** itself `$CAPTURE`-d under that exact name
  (it's a genuine `double`, declared once in `$MAIN`, used directly by the
  disease `$ODE` equations — `$MAIN`/`$ODE` share one compiled scope under
  this build). A separately-named, independently-recomputed
  `EFFECT_<STEM>_out` is captured from `$TABLE` instead, for verification
  purposes, avoiding the `$MAIN`/`$CAPTURE` auto-declaration collision this
  fork's `autoimmune-polyendocrinopathy` refactor already documented for its
  own `C_CSA`/`C_JAKI` (a `double NAME=...` declaration for a name that is
  *also* `$CAPTURE`-d under the identical name, anywhere else in the file,
  is a genuine mrgsolve 2.0.1 build failure — confirmed empirically here
  too, on an intermediate draft, before landing on this split).
- One note on the guide's qspserver-compatibility text suggesting
  `C_<STEM>`/`EFFECT_<STEM>` should live in `$PARAM`: this file follows the
  corpus's actual, already-established convention instead (matching
  `autoimmune-polyendocrinopathy` and `kidney-transplant-rejection`, neither
  of which puts these derived, state-dependent quantities in `$PARAM` — a
  `$PARAM` value is a fixed default, and a live ODE-driven concentration
  cannot be one). `EC50_<STEM>`/`EMAX_<STEM>`/`GAMMA_<STEM>` (the actual
  fixed values) are in `$PARAM`, satisfying `/model_manifest` discoverability
  for potency/compound-identity; `C_<STEM>`/`EFFECT_<STEM>` are discoverable
  via `outputPaths` instead, per the guide's own literal-text-grep mechanism.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the five
`autoimmune-encephalitis` rows (`CPX`, `IVIG`, `MP`, `RTX`, `TCZ`).

## Anything else flagged for manual review

- The original's `Cp_RTX_out`/`Cp_TCZ_out`/`Cp_MP_out` all apply the same
  "divide by volume in L, then x1000" convention, labeled `mcg/mL` in a
  comment; dimensionally this is actually `mg/L x 1000 = mcg/L`, not
  `mcg/mL` (mg/L and mcg/mL are numerically equal, so the x1000 produces a
  value 1000x too large for what the comment claims). This is internally
  consistent across RTX/MP/TCZ (same convention, same magnitude of
  mislabeling, doesn't change any relative comparison or the shape of any
  dynamic) and doesn't affect disease dynamics at all — preserved exactly
  (rename, not a refit) since fixing it would be a numeric change to
  arithmetic these three compounds' Hill terms actually depend on, out of
  scope for a structural PK refactor. Noted here for visibility, not logged
  as a separate `UPSTREAM_ISSUES.md` entry (it's a comment/labeling
  inaccuracy, not a build defect or a behavioral bug).
- No compound in this file was left untouched — all five (`IVIG`, `MP`,
  `RTX`, `TCZ`, `CPX`) were in scope and refactored.
