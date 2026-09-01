# Refactor notes — `stable-angina/sa_mrgsolve_model.R`

Five compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **BB** (bisoprolol,
beta-adrenergic receptor), **CCB** (amlodipine, calcium channel blocker),
**IVA** (ivabradine, I_f current), **NIT** (ISMN/nitrate, preload
reduction + coronary vasodilation), and **RAN** (ranolazine, late sodium
current / I_NaL). These are the only compounds this file models — nothing
else was touched. Every compartment already present in the original
(21 total) is still present, in the same order; none were added or
removed, so 1-based compartment numbers are unchanged between the
original and the refactored file.

## Archetypes

- **BB (bisoprolol)** — Archetype 3, depot + central + peripheral, linear
  elimination (`GUT_BB`/`CENT_BB`/`PERI_BB`, was `GUT_BB`/`CENTRAL_BB`/
  `PERIPH_BB`). Two independent, genuinely separate pharmacological
  actions already present in the original as plain ratios: HR reduction
  and SBP reduction. Named `EFFECT_BB_HR` and `EFFECT_BB_SBP` (see naming
  table below) — this is not "combining drugs into one shared Hill term"
  (which the guide forbids); it is one drug's two distinct, independently
  named effects, each reading the same single `C_BB`.
- **CCB (amlodipine)** — Archetype 3, depot + central + peripheral, linear
  elimination (`GUT_CCB`/`CENT_CCB`/`PERI_CCB`, was `GUT_CCB`/
  `CENTRAL_CCB`/`PERIPH_CCB`). Two effects: SBP reduction and coronary
  vasodilation (CBF increase). Named `EFFECT_CCB_SBP`/`EFFECT_CCB_CBF`.
- **RAN (ranolazine)** — Archetype 3 minus the peripheral compartment
  (depot + central, linear elimination; no peripheral compartment in the
  original, `GUT_RAN`/`CENT_RAN`, was `GUT_RAN`/`CENTRAL_RAN`). Two
  effects: I_NaL inhibition and exercise-capacity increase (CARISA). Named
  `EFFECT_RAN_INA`/`EFFECT_RAN_EX`.
- **IVA (ivabradine)** — Archetype 3 minus the peripheral compartment
  (depot + central, linear elimination; `GUT_IVA`/`CENT_IVA`, was
  `GUT_IVA`/`CENTRAL_IVA`). One effect (HR reduction via I_f inhibition):
  `EFFECT_IVA`.
- **NIT (ISMN nitrate)** — Archetype 3 minus the peripheral compartment
  (depot + central, linear elimination; `GUT_NIT`/`CENT_NIT`, was
  `GUT_NIT`/`CENTRAL_NIT`) **plus a bespoke tolerance state native to
  NIT's own block**, `TOL_NIT` (was `NIT_TOL`) — not covered by the
  guide's naming table, since it is neither a PK compartment nor a
  disease-side PD compartment: it is a first-order state that tracks
  continuous-vs-intermittent nitrate exposure and multiplicatively
  discounts both of NIT's own effect terms (the original's "nitrate
  tolerance" mechanism, kept exactly). Two effects: MVO2/preload reduction
  and coronary vasodilation. Named `EFFECT_NIT_MVO`/`EFFECT_NIT_CBF`, each
  multiplied by the tolerance discount factor (renamed `TOLEFF_NIT`, was
  the local `eff_NIT`) exactly as the original did.

No archetype-4 (TMDD) compound here — none of the five models receptor
binding as its own dynamic state; every effect is a plain concentration
ratio, confirmed by reading every `dxdt_*` line in the file.

## Naming (values unchanged from the original)

| Role | Original | Refactored |
|---|---|---|
| BB compartments | `GUT_BB`/`CENTRAL_BB`/`PERIPH_BB` | `GUT_BB`/`CENT_BB`/`PERI_BB` |
| BB PK params | `KA_BB`/`F_BB`/`CL_BB`/`V1_BB`/`Q_BB`/`V2_BB` | unchanged (already matched convention) |
| BB concentration | `Cp_BB` (`$ODE`-local, self-colliding capture) | `C_BB` |
| BB HR effect | `EC50_BB`/`EMAX_BB` | `EC50_BB_HR`/`EMAX_BB_HR` (+ new `GAMMA_BB_HR=1`) |
| BB SBP effect | `EC50_SBP_BB`/`EMAX_SBP_BB` | `EC50_BB_SBP`/`EMAX_BB_SBP` (+ new `GAMMA_BB_SBP=1`) |
| BB effect terms | `E_BB_HR`/`E_BB_SBP` | `EFFECT_BB_HR`/`EFFECT_BB_SBP` |
| CCB compartments | `GUT_CCB`/`CENTRAL_CCB`/`PERIPH_CCB` | `GUT_CCB`/`CENT_CCB`/`PERI_CCB` |
| CCB PK params | `KA_CCB`/`F_CCB`/`CL_CCB`/`V1_CCB`/`Q_CCB`/`V2_CCB` | unchanged |
| CCB concentration | `Cp_CCB` | `C_CCB` |
| CCB SBP effect | `EC50_CCB`/`EMAX_CCB` | `EC50_CCB_SBP`/`EMAX_CCB_SBP` (+ new `GAMMA_CCB_SBP=1`) |
| CCB CBF effect | `EC50CBF_CCB`/`EMAX_CBF_CCB` | `EC50_CCB_CBF`/`EMAX_CCB_CBF` (+ new `GAMMA_CCB_CBF=1`) |
| CCB effect terms | `E_CCB_SBP`/`E_CCB_CBF` | `EFFECT_CCB_SBP`/`EFFECT_CCB_CBF` |
| RAN compartments | `GUT_RAN`/`CENTRAL_RAN` | `GUT_RAN`/`CENT_RAN` |
| RAN PK params | `KA_RAN`/`F_RAN`/`CL_RAN`/`V1_RAN` | unchanged |
| RAN concentration | `Cp_RAN` | `C_RAN` |
| RAN I_NaL effect | `EC50_RAN`/`EMAX_RAN` | `EC50_RAN_INA`/`EMAX_RAN_INA` (+ new `GAMMA_RAN_INA=1`) |
| RAN exercise effect | `EC50EX_RAN`/`EMAX_EX_RAN` | `EC50_RAN_EX`/`EMAX_RAN_EX` (+ new `GAMMA_RAN_EX=1`) |
| RAN effect terms | `E_RAN_INa`/`E_RAN_EX` | `EFFECT_RAN_INA`/`EFFECT_RAN_EX` |
| IVA compartments | `GUT_IVA`/`CENTRAL_IVA` | `GUT_IVA`/`CENT_IVA` |
| IVA PK params | `KA_IVA`/`F_IVA`/`CL_IVA`/`V1_IVA` | unchanged |
| IVA concentration | `Cp_IVA` | `C_IVA` |
| IVA effect | `EC50_IVA`/`EMAX_IVA` | unchanged names (+ new `GAMMA_IVA=1`) |
| IVA effect term | `E_IVA_HR` | `EFFECT_IVA` |
| NIT compartments | `GUT_NIT`/`CENTRAL_NIT`/`NIT_TOL` | `GUT_NIT`/`CENT_NIT`/`TOL_NIT` |
| NIT PK params | `KA_NIT`/`F_NIT`/`CL_NIT`/`V1_NIT` | unchanged |
| NIT concentration | `Cp_NIT` | `C_NIT` |
| NIT MVO2 effect | `EC50_NIT`/`EMAX_NIT` | `EC50_NIT_MVO`/`EMAX_NIT_MVO` (+ new `GAMMA_NIT_MVO=1`) |
| NIT CBF effect | `EC50CBF_NIT`/`EMAX_CBF_NIT` | `EC50_NIT_CBF`/`EMAX_NIT_CBF` (+ new `GAMMA_NIT_CBF=1`) |
| NIT tolerance rates | `TOL_K1`/`TOL_K2` | `TOLK1_NIT`/`TOLK2_NIT` |
| NIT tolerance discount | `eff_NIT` (`$ODE`-local) | `TOLEFF_NIT` |
| NIT effect terms | `E_NIT_MVO`/`E_NIT_CBF` | `EFFECT_NIT_MVO`/`EFFECT_NIT_CBF` |

## Hill interface: renames, not fits

All nine effect terms (two each for BB/CCB/RAN/NIT, one for IVA) were
already written as a plain concentration ratio `Emax*C/(EC50+C)` in the
original — no ODE-solved receptor kinetics to approximate. Every one is a
rename, confirmed by the exact verification match below (not merely
"close"): `GAMMA_<STEM>_<EFFECT> = 1` was added explicitly for each
(matching this fork's established convention of writing the effect as
`EMAX*pow(C,GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))` even when `GAMMA=1`,
so the shape is visibly a named Hill interface rather than a bespoke
ratio) — algebraically identical to the original's plain ratio since
`pow(x,1) = x`. No `nls()` fitting was needed or performed for any of the
five compounds.

## Disease-side equations: unchanged shape, renamed inputs only

`HR_tgt`, `SBP_tgt`, `CBF_tgt`, `MVO2_mod`, `MVO2_eff`, `isch_tgt`, and
`ex_tgt` are byte-identical arithmetic to the original, with only the
renamed `EFFECT_*` terms substituted in place of the original's `E_*`
locals. Nothing about how the disease system combines multiple drugs'
effects was changed.

## `$PARAM` vs `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>`

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` (with a `= 0` default) for direct
`/model_manifest` discoverability. This was not done here, for the same
reason already documented in the AMD/membranous-nephropathy/sepsis/SAH/
AIP/breast-cancer refactors: mrgsolve 2.0.1 compiles `$PARAM` members as
**read-only references** inside `$ODE`, so a value that must be
recomputed every timestep from state (`C_BB`, `C_CCB`, `C_RAN`, `C_IVA`,
`C_NIT`, and all nine `EFFECT_*` terms) cannot also be declared in
`$PARAM`. Independently re-confirmed against the live `mrgsolve_api`
container with a minimal reproduction (a `$PARAM ... C_X = 0.0`
declaration reassigned in `$ODE` fails to build with `error: assignment
of read-only reference 'C_X'`). Instead, all fourteen are predeclared as
`double`s in `$GLOBAL` and listed bare in `$CAPTURE`: visible in every
simulation's output columns and in `/model_manifest`'s `outputPaths`,
confirmed directly —

```
"outputPaths": [..., "C_BB","C_CCB","C_RAN","C_IVA","C_NIT",
"EFFECT_BB_HR","EFFECT_BB_SBP","EFFECT_CCB_SBP","EFFECT_CCB_CBF",
"EFFECT_IVA","EFFECT_NIT_MVO","EFFECT_NIT_CBF","EFFECT_RAN_INA",
"EFFECT_RAN_EX", ...]
```

— just not in its `parameters` list.

## Pre-existing upstream build defect (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for a
reason unrelated to any of the five refactored compounds' own PK: `$ODE`
declares `double Cp_BB = CENTRAL_BB / V1_BB * 1000.0;` (and identically
for `Cp_CCB`/`Cp_RAN`/`Cp_IVA`/`Cp_NIT`), then `$TABLE` immediately writes
`capture Cp_BB = Cp_BB;` — a self-referential `capture NAME = NAME;` that
collides with the already-declared `double` of the same name:

```
68:11: error: redefinition of 'capture {anonymous}::Cp_BB'
   68 |   capture Cp_BB;
      |           ^~~~~
30:10: note: 'double {anonymous}::Cp_BB' previously declared here
   30 |   double Cp_BB;
      |          ^~~~~
```

(and identically four more times, for `Cp_CCB`/`Cp_RAN`/`Cp_IVA`/
`Cp_NIT`). Confirmed via `POST /model_manifest` on the untouched
original's own DSL, no changes involved. Logged as
[`UPSTREAM_ISSUES.md` #80](../translations/UPSTREAM_ISSUES.md). Per the
guide's settled policy for this situation, the fix is applied **directly
to the delivered `sa_mrgsolve_model_refactored.R`**, not just to a
scratch copy: rather than reproducing the same collision under the new
names (`capture C_BB = C_BB;` would fail identically), `C_BB`/`C_CCB`/
`C_RAN`/`C_IVA`/`C_NIT` and all nine `EFFECT_*` terms are predeclared in
`$GLOBAL`, assigned by plain statement in `$ODE`, and exposed via a bare
`$CAPTURE` list — the pattern already established in six prior refactors
in this fork (see the `$PARAM` vs `$GLOBAL` section above). This is
syntax-only and non-numeric. The checked-in original
(`sa_mrgsolve_model.R`) was left untouched and still carries this defect.

For the *verification* run below, a syntax-fixed scratch copy of the
**original** was also needed (since the untouched original does not
compile at all, there is nothing to compare against otherwise): the five
colliding `$ODE`-local doubles were renamed `Cp_BB_val`/`Cp_CCB_val`/
`Cp_RAN_val`/`Cp_IVA_val`/`Cp_NIT_val` (every downstream read updated to
match), while every `$TABLE` capture's own output name (`Cp_BB`, `E_BB`,
etc.) was left exactly as the original wrote it — this scratch copy
existed only in memory for the API calls below, was never written to
`sa_mrgsolve_model.R`, and is not part of the delivered file set.

## Verification

**Method.** Both the syntax-fixed scratch copy of the original (fix
described above, applied only for the purpose of getting it to compile at
all) and `sa_mrgsolve_model_refactored.R`'s embedded DSL were extracted as
bare mrgsolve DSL text and run through the qspserver `mrgsolve_api`
container (`POST /model_manifest`, `POST /run_simulation`) at
`http://localhost:8007`, spaced ~2s apart per request. The extraction was
confirmed byte-identical to the `_refactored.R`'s embedded quoted string
before running (checked, then discarded — no `.cpp` file is part of the
delivered set). Compartment numbers are identical (1-based) between the
original and the refactored file, since no compartment was added or
removed — only renamed — so no index remapping was needed for dosing.

All six of the original file's own named dosing scenarios were run,
exactly as defined in its own `make_regimen()`/`scenarios` list (dose
amounts, `ii`, `addl` all reproduced from the R script's own arithmetic;
`ii=24,addl=6` for QD dosing and `ii=12,addl=13` for BID/eccentric dosing
over the file's own `t_end=168`), over the file's own 1-week horizon
(`end=168`, `delta=0.5`, 337-340 points depending on dose timing):

1. **S1** — Untreated
2. **S2** — Bisoprolol 5 mg QD
3. **S3** — Bisoprolol + Amlodipine 5 mg
4. **S4** — Bisoprolol + Amlodipine + Ranolazine 1000 mg BID
5. **S5** — Bisoprolol + Ivabradine 5 mg BID
6. **S6** — Bisoprolol + Amlodipine + ISMN 40 mg (eccentric) + statin

Two additional analyses the original file's own R code also runs were
verified as well:

7. **Bisoprolol dose-response sweep** (1.25/2.5/5/10/20 mg QD, monotherapy,
   same 1-week horizon) — the original's own `dr_results` loop.
8. **2-year plaque progression, statin on/off** (`end=365*24*2=17520`,
   `delta=24`, 731 points, no dosing) — the original's own
   `long_no_statin`/`long_statin` runs.

Every shared `$CAPTURE`d output was compared point-by-point across the
full time grid for all eight runs: all 13 renamed compartments
(`GUT_BB`→`GUT_BB`, `CENTRAL_BB`→`CENT_BB`, `PERIPH_BB`→`PERI_BB`, ...,
`NIT_TOL`→`TOL_NIT`), the 8 disease-state compartments (unchanged names),
the 5 concentrations (`Cp_BB`→`C_BB`, etc.), the diagnostic aliases
(`HR_sim`/`SBP_sim`/`RPP_sim`/`CBF_sim`/`O2imbal`/`Isch`/`Angina`/
`ExCap`/`Plaque`/`NIT_tol`, unchanged names on both sides), and the three
effect diagnostics the original happened to capture (`E_BB`→
`EFFECT_BB_HR`, `E_IVA`→`EFFECT_IVA`, `E_RAN`→`EFFECT_RAN_INA`) — 37
paired outputs total.

**Result: exact match, max abs diff = 0.0 for every output, every one of
the 8 runs (6 named scenarios + dose-response sweep + 2-year plaque run),
no solver timeouts or step-count issues encountered at any horizon
(largest run: 731 points over a 2-year/17,520-hour window).** This is the
expected outcome per the guide's tolerance rule for pure structural
reorganization with no Hill-fitting — every renamed effect term computes
the identical arithmetic on the identical (renamed) inputs, and no
compartment was added, removed, or restructured, so there was no source
of numerical divergence to begin with.

## Collapsed redundant diagnostic captures (verified, not a numeric change)

The original's `$TABLE` separately recomputed `E_BB`, `E_IVA`, `E_RAN` —
each an exact duplicate of the already-computed `E_BB_HR`/`E_IVA_HR`/
`E_RAN_INa` internal locals (identical formula, identical inputs), purely
as extra diagnostic capture columns never read by the R-side scenario or
plotting code. In the refactored file these are not recomputed a second
time — `EFFECT_BB_HR`, `EFFECT_IVA`, and `EFFECT_RAN_INA` (their renamed
equivalents) are already `$GLOBAL`-declared and `$CAPTURE`d once, so the
duplicate `$TABLE` recomputation was dropped rather than tripled up under
new names. Confirmed harmless by the exact-match verification above
(`E_BB`↔`EFFECT_BB_HR`, `E_IVA`↔`EFFECT_IVA`, `E_RAN`↔`EFFECT_RAN_INA`
all matched to 0.0 across every scenario) — this changes what gets
computed twice, not what value ends up in any output column.

## Anything else worth flagging

- **`DOSE_BB`/`DOSE_CCB`/`DOSE_RAN`/`DOSE_IVA`/`DOSE_NIT` are declared in
  `$PARAM` but never read anywhere in `$ODE`/`$MAIN`/`$TABLE`.** All five
  dosing switches are dead parameters — the R script sets them per
  scenario (`DOSE_BB = ifelse(s$bb > 0, 1, 0)`, etc.) but nothing in the
  model ever consults them; actual dosing is driven entirely by the
  `ev()` events into the `GUT_*` compartments. Harmless (no build issue,
  no numeric effect — confirmed since the verification above never even
  passed these params for most runs and still matched exactly), so not
  logged as its own `UPSTREAM_ISSUES.md` entry, just disclosed here.
  Preserved unchanged (still declared, still unread) in the refactored
  file, since fixing dead-but-harmless parameters is out of scope for a
  rename-only refactor.
- The header comment's compartment count (originally "20 ODE
  compartments") was corrected to 21 in the refactored file's `$PROB`
  block, matching an actual count of the `$CMT` list (a pre-existing,
  purely cosmetic inaccuracy in the original's own comment, not a gated
  output).
- The R-side scenario/plotting code in `sa_mrgsolve_model_refactored.R`
  was updated only where it referenced a renamed identifier: the PK
  concentration plot's `pivot_longer(c(Cp_BB, Cp_CCB, ...))` now reads
  `c(C_BB, C_CCB, C_RAN, C_IVA, C_NIT)` and its scenario filter reads
  `C_BB`/`C_CCB`/`C_NIT`. All `ev()`/dosing calls (which address
  compartments by name, unaffected by any renaming since the same names
  are used consistently on both sides) and all six scenario definitions
  are otherwise identical to the original.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`stable-angina | BB`, `stable-angina | CCB`, `stable-angina | IVA`,
`stable-angina | NIT`, and `stable-angina | RAN` rows.

## Discoverability fix

A corpus-wide discoverability audit found `C_BB`, `C_CCB`, `C_RAN`,
`C_IVA`, and `C_NIT` were not written as single contiguous
`double C_<STEM> = <expr>;` statements anywhere in the file: all five
were `$GLOBAL`-declared as one bare forward declaration
(`double C_BB, C_CCB, C_RAN, C_IVA, C_NIT;`) with the actual values
assigned via bare reassignment (no `double`) once each in `$ODE`, and
never recomputed in `$TABLE` — a legitimate, working pattern, but not
literal-text-discoverable by tooling that regexes for
`double C_<STEM> = ...;`.

**Fix applied, per the standing recipe for this pattern:**
1. Removed the entire `double C_BB, C_CCB, C_RAN, C_IVA, C_NIT;` line
   from `$GLOBAL` (all five names in this batch were target compounds,
   so the whole line was removed; the adjacent `EFFECT_*` forward
   declares on the next line were untouched).
2. Since the file had no existing `$TABLE` recompute for any of the
   five (single-bare-assignment-site, per the audit), added five new
   lines to `$TABLE`, immediately before `$CAPTURE`, reusing the exact
   expression already used in `$ODE` for each:
   `double C_BB = CENT_BB / V1_BB * 1000.0;`,
   `double C_CCB = CENT_CCB / V1_CCB * 1000.0;`,
   `double C_RAN = CENT_RAN / V1_RAN * 1000.0;`,
   `double C_IVA = CENT_IVA / V1_IVA * 1000.0;`,
   `double C_NIT = CENT_NIT / V1_NIT * 1000.0;`.
   The `$ODE` bare reassignments and the bare `$CAPTURE` listing
   (`C_BB C_CCB C_RAN C_IVA C_NIT`) were left untouched — mrgsolve's
   DSL parser hoists any `double NAME = expr;` found anywhere in the
   block's text into a shared class-member declaration, so the new
   `$TABLE` declarations serve as the one declaration site both the
   `$ODE` bare reassignments and the bare `$CAPTURE` entries resolve
   against.

**Verification:**
- `Rscript -e 'parse("sa_mrgsolve_model_refactored.R")'` succeeds with
  no error (one straight apostrophe introduced in a first-draft comment
  for this fix broke the string body per this guide's own rule; fixed
  by switching to a curly apostrophe before re-parsing).
- `grep` confirms `double C_BB = `, `double C_CCB = `, `double C_RAN =
  `, `double C_IVA = `, `double C_NIT = ` and the matching
  `EC50_BB_HR`/`EC50_CCB_SBP`/`EC50_RAN_INA`/`EC50_IVA`/`EC50_NIT_MVO`
  all appear in the file.
- Extracted the DSL block (before and after this fix) and posted both
  to qspserver's `mrgsolve_api` `/model_manifest` (both compiled) and
  `/run_simulation` (both ran) against the file's own "S6: BB + ISMN 40
  mg + CCB" scenario (`make_regimen(bb_mg=5, ccb_mg=5, nit_mg=40,
  statin=1, t_end=168)`: bisoprolol 5 mg QD into `GUT_BB`, amlodipine 5
  mg QD into `GUT_CCB`, ISMN 40 mg BID into `GUT_NIT`, `STATIN_ON=1`;
  `delta = 0.5, end = 168`). Every `$CAPTURE`d column — all five
  `C_<STEM>` concentrations, all nine `EFFECT_<STEM>` terms, and every
  PD/clinical output — is numerically **identical** (max abs diff = 0
  at all 340 reported rows) between the pre-fix and post-fix DSL. This
  is a pure declaration-site move; no numeric behavior changed.
