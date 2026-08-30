# Refactor notes — `neurofibromatosis-type-1/nf1_mrgsolve_model.R`

Two compounds in this file, both classified "Redirect concentration (clean
single site)" in `driver-patches/data/compound_perturbation_census.md`:
**Mirdametinib** and **Selumetinib**. Both are MEK1/2 inhibitors, refactored
independently (each keeps its own `C_<STEM>`/`EFFECT_<STEM>` pair — never
collapsed into a shared Hill term). There is no third compound in this file.

## Archetype (both compounds)

**Archetype 3 minus the peripheral compartment** — depot + central, linear
elimination — confirmed against the actual `$ODE` lines, not assumed from
the file's own "3-cpt oral" header comment (there is no peripheral
compartment, no `Q`/`V2` term, for either drug anywhere in the model):

- Selumetinib: `SEL_GUT`/`SEL_CENT` (depot + central only) with
  `dxdt_SEL_CENT = KA_SEL*SEL_GUT - (CL_SEL/V_SEL)*SEL_CENT` — plain linear
  elimination, no receptor/target-binding term.
- Mirdametinib: same shape, `MIR_GUT`/`MIR_CENT`.

Renamed:

| Original | Refactored |
|---|---|
| `SEL_GUT` | `GUT_SEL` |
| `SEL_CENT` | `CENT_SEL` |
| `V_SEL` | `V1_SEL` |
| `SEL_CP` | `C_SEL` |
| `SEL_INHIB` | `EFFECT_SEL` |
| `MIR_GUT` | `GUT_MIR` |
| `MIR_CENT` | `CENT_MIR` |
| `V_MIR` | `V1_MIR` |
| `MIR_CP` | `C_MIR` |
| `MIR_INHIB` | `EFFECT_MIR` |

`KA_SEL`, `CL_SEL`, `F_SEL`, `EC50_SEL`, `KA_MIR`, `CL_MIR`, `F_MIR`,
`EC50_MIR` already matched the convention and are unchanged. All parameter
*values* are copied verbatim from the original — nothing invented, nothing
dropped. Added (new, not in the original, for the named Hill interface):
`EMAX_SEL=1`, `GAMMA_SEL=1`, `EMAX_MIR=1`, `GAMMA_MIR=1`.

### A subtlety: mrgsolve's `F_<cmt>` bioavailability magic variable

The original's `$MAIN` block contains:

```
F_SEL_GUT = F_SEL * ADHERENCE_SEL;
F_MIR_GUT = F_MIR * ADHERENCE_MIR;
```

This is not dead code and not a naming coincidence: mrgsolve recognizes a
`$MAIN` variable named `F_<compartment>` as the per-dose bioavailability
fraction for that compartment. Since the compartment is named `SEL_GUT`,
`F_SEL_GUT` is the exact name mrgsolve looks for — this is how
`ADHERENCE_SEL`/`ADHERENCE_MIR` actually reach the dosed amount (there is no
`F_SEL`/`F_MIR` multiplication inside the `$ODE` math itself). Renaming the
compartment to `GUT_SEL`/`GUT_MIR` therefore *requires* renaming this magic
variable identically, to `F_GUT_SEL`/`F_GUT_MIR` — missing this would
silently stop adherence/bioavailability from applying to either drug's
dosing, with no compile error to flag it. Done in the refactored file;
confirmed correct by the exact-match verification below (both original and
refactored ran with `ADHERENCE_SEL = ADHERENCE_MIR = 1.0`, matching
scenarios' default, so this path was exercised identically in both).

## Hill interface: rename, not a fit — for both compounds

The original's combined-inhibition block:

```
double SEL_INHIB = SEL_CP / (EC50_SEL + SEL_CP);
double MIR_INHIB = MIR_CP / (EC50_MIR + MIR_CP);
double TOTAL_INHIB = 1.0 - (1.0 - SEL_INHIB) * (1.0 - MIR_INHIB);
```

`SEL_INHIB` and `MIR_INHIB` are already exactly `Emax*C/(EC50+C)` with an
implicit `Emax=1`, `gamma=1` — a plain ratio, not something emerging from
ODE-solved receptor kinetics. Per the guide, this is a rename, not a fit:

```
double EFFECT_SEL = EMAX_SEL * pow(C_SEL, GAMMA_SEL) / (pow(EC50_SEL, GAMMA_SEL) + pow(C_SEL, GAMMA_SEL));
double EFFECT_MIR = EMAX_MIR * pow(C_MIR, GAMMA_MIR) / (pow(EC50_MIR, GAMMA_MIR) + pow(C_MIR, GAMMA_MIR));
double TOTAL_INHIB = 1.0 - (1.0 - EFFECT_SEL) * (1.0 - EFFECT_MIR);
```

With `EMAX_*=1`, `GAMMA_*=1` this is arithmetically identical to the
original line-for-line — confirmed by the exact-match verification below
(max abs diff 0.0, not merely small). `TOTAL_INHIB`'s Bliss-independence
combination structure itself is untouched; only its two inputs were
renamed, so each compound stays independently driveable (e.g. an external
covariate could substitute `C_MIR` alone without touching `EFFECT_SEL`'s
computation at all).

Sanity check from the verification runs: in the selumetinib-only scenario,
`EFFECT_MIR = 0` throughout (no mirdametinib dosed) and `EFFECT_SEL`
matches `TOTAL_INHIB` exactly at every timepoint (confirmed numerically,
e.g. both read `0.2306` at hour 168); the mirror holds for the
mirdametinib-only scenario.

## A pre-existing defect found while verifying (not fixed, logged upstream)

**`nf1_mrgsolve_model.R`, as checked in, does not compile under mrgsolve
2.0.1 at all** — confirmed via the qspserver `mrgsolve_api` service, on the
untouched original's own DSL, with two entangled defects (full detail,
including a minimal isolated reproduction, in
`translations/UPSTREAM_ISSUES.md` entry #41):

1. `$CAPTURE` lists ten names that duplicate `$CMT` compartment names
   (`RESIST`, `OPG_VOL`, `CNF_BURDEN`, `PAIN`, `QOL`, `VISION`, `LVEF`,
   `DERM_AE`, `CPK_AE`, `GROWTHZ`) — mrgsolve 2.0.1 refuses to build any
   model whose `$CAPTURE` repeats a compartment name (the same defect class
   seen in several other files in this corpus, e.g. `rheumatoid-arthritis`,
   `chronic-hypothyroidism`, `sepsis`).
2. **Once (1) is fixed in isolation, a second, distinct defect appears**:
   `$MAIN`'s initial-condition block uses bare compartment-name assignment
   (`PERK = PERK_BASE;`, `RESIST = 0;`, etc., 13 lines) instead of the
   standard mrgsolve `<cmt>_0` idiom. This bare form only compiles when the
   same name is *also* present in `$CAPTURE` — i.e. it silently depends on
   defect 1 being present. Confirmed with a minimal isolated 1-compartment
   reproduction: `X = XBASE;` in `$MAIN` compiles when `X` is duplicated
   into `$CAPTURE` (only failing later, at the `validObject` step — the
   same symptom as defect 1), but fails to compile with `error: assignment
   of read-only reference 'X'` once `X` is removed from `$CAPTURE`.
   Switching to the standard `X_0 = XBASE;` idiom compiles regardless of
   `$CAPTURE` membership, and was confirmed (via `/run_simulation` on the
   minimal model) to correctly initialize the state to the intended value.

Because defect 1 always fires first, the original file — as checked in —
has never been buildable under mrgsolve 2.0.1 by either route, so it is not
possible to observe from the original alone whether its `$MAIN` block's
bare-assignment idiom would otherwise have worked.

**Workaround applied** (per `FORK_WORKFLOW_GUIDE.md`'s instruction to fix
mechanically for a working verification harness, not upstream): the ten
duplicate names were removed from `$CAPTURE` (compartment states always
appear in mrgsolve's output regardless of `$CAPTURE`, so nothing is lost),
and all 13 `$MAIN` initial-condition lines were switched to the `<cmt>_0`
idiom (e.g. `PERK_0 = PERK_BASE;`). Neither change touches a parameter
value or a model dynamic. This pair of changes was applied to **both**: (a)
an in-memory-only scratch copy of the *original* file's DSL, used solely to
build a working comparison target, and (b) the delivered
`nf1_mrgsolve_model_refactored.R` itself — required so the delivered file
is actually buildable/usable through the qspserver `mrgsolve_api`, per the
guide's qspserver-compatibility section. `nf1_mrgsolve_model.R` itself was
never edited.

## qspserver `/model_manifest` discoverability

`POST /model_manifest` on the refactored DSL confirms every covariate this
guide's naming convention promises is present:

- `outputPaths` includes `C_SEL`, `C_MIR`, `EFFECT_SEL`, `EFFECT_MIR` (plus
  every renamed compartment and the other derived doubles) — all four are
  `$ODE`-computed, state-dependent doubles added directly to `$CAPTURE`
  (no name collision to work around here, since `C_SEL`/`C_MIR`/
  `EFFECT_SEL`/`EFFECT_MIR` are new names that don't duplicate any `$CMT` or
  `$PARAM` identifier — unlike some other files in this corpus that needed
  an `_out`-suffixed alias to avoid a `$PARAM`/`$ODE` collision).
- `parameters` lists `KA_SEL`, `CL_SEL`, `V1_SEL`, `F_SEL`, `EC50_SEL`,
  `EMAX_SEL`, `GAMMA_SEL` and the mirrored `_MIR` set, all with their
  original/derived default values, confirming `$PARAM` exposure of every
  renamed and newly-added parameter.

`C_SEL`/`C_MIR`/`EFFECT_SEL`/`EFFECT_MIR` themselves are **not** `$PARAM`
entries (state-dependent quantities recomputed from `CENT_SEL`/`CENT_MIR`
every step cannot be fixed `$PARAM` values — the same constraint documented
in `chronic-hypothyroidism/hypo_refactor_notes.md` and
`thyroid-eye-disease/ted_refactor_notes.md`); they are discoverable instead
via `$CAPTURE`/`outputPaths`, confirmed above.

## Verification

**Method.** Both `nf1_mrgsolve_model.R`'s and
`nf1_mrgsolve_model_refactored.R`'s embedded mrgsolve DSL blocks were
extracted mechanically (regex on the `nf1_code <- '...'` assignment, quoted
contents pulled verbatim, R's `\'` string-escaping undone) and POSTed to the
qspserver `mrgsolve_api` service (`http://localhost:8007`, confirmed
healthy) — `/model_manifest` first (both compile after the workaround
above), then `/run_simulation`.

Ran **2 of the original file's own scenarios** (dose amounts/timing copied
exactly from the original R script's own `scenarios` list, using
`BSA_PED = 1.10`, translated to the API's `ii`/`addl` dosing convention),
each exercising one compound at a time:

1. **Scenario 2** — Selumetinib 25 mg/m2 BID, pediatric, SPRINT-style
   (27.5 mg into `GUT_SEL`/`SEL_GUT`, q12h, addl=112).
2. **Scenario 4** — Mirdametinib 2 mg/m2 BID, pediatric, ReNeu-style
   (2.2 mg into `GUT_MIR`/`MIR_GUT`, q12h, addl=42).

Both ran the original file's own full 96-week/16128 h horizon at its own
`delta=24` h output grid — **no shortening needed**: both scenarios
returned in well under a second via the API, with no solver step-count
issue (the guide's known `maxsteps` caveat did not apply here).

Every `$CAPTURE`d output plus every raw compartment (all 19 states +
`C_SEL`/`C_MIR`/`TOTAL_INHIB`/`PERK_SUPPRESSION`/`PN_TOTAL`/
`PN_RESPONSE_PCT`, mapped name-for-name against the original's
`SEL_CP`/`MIR_CP`/etc.) was compared point-by-point across the full
674-point time grid for both scenarios.

**Result: exact match, max abs diff = 0.0 for every output, both
scenarios.** Expected for a pure structural reorganization (archetype 3
minus peripheral, both compounds) plus a zero-approximation Hill-term
rename, per the guide's tolerance rule for archetypes 1–3.

R syntax sanity check (not a substitute for the API verification above,
just confirming the new file parses as valid R): `Rscript -e
"parse('nf1_mrgsolve_model_refactored.R')"` succeeds, 10 top-level
expressions — identical count to `nf1_mrgsolve_model.R`.

## Anything else worth flagging

- No third compound in this file — Mirdametinib and Selumetinib are the
  only two PK/PD-modeled drugs in `nf1_mrgsolve_model.R` (the file's
  scenario 9, "Trametinib_offlabel_approx", re-parameterizes the
  mirdametinib PK/PD block rather than modeling trametinib independently —
  this R-level approximation is untouched, and its `V_MIR` reference in the
  commented-out example-run block was updated to `V1_MIR` to match the
  rename, with no other change).
- Compartment ordering is unchanged (`GUT_SEL` is still compartment 1,
  `CENT_SEL` 2, `GUT_MIR` 3, `CENT_MIR` 4, etc.) since renaming was done in
  place without reordering, and all ten disease/PD compartments after it
  keep their original names entirely (only the two drugs' own PK blocks and
  the two Hill-interface doubles were touched).
- Scenario `cmt=` strings in the R-level `scenarios` list (`"SEL_GUT"` /
  `"MIR_GUT"`) were updated to `"GUT_SEL"` / `"GUT_MIR"` to match the
  renamed compartments; dose amounts, `ii`, and `addl` values are otherwise
  byte-identical to the original.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`,
`neurofibromatosis-type-1 | Mirdametinib` and
`neurofibromatosis-type-1 | Selumetinib` rows.
