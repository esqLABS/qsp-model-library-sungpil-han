# Refactor notes — `ktx_mrgsolve_model.R` (calibration run: tocilizumab only)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only
tocilizumab's ("TCZ", modeled here as "anti-IL-6(R)", `CLTCZ`/`VTCZ`/`TCZ_C`)
PK block and its interface into the disease equations were rewritten. Every
other compound in the file (tacrolimus, cyclosporine, mycophenolate,
prednisolone/methylprednisolone, everolimus, rabbit ATG, basiliximab,
belatacept, rituximab, anti-CD38/felzartamab, IVIG, ganciclovir) was copied
verbatim — confirmed by diff, see below.

## Archetype determination

**Archetype 1 — no depot, single compartment, linear elimination**, combined
with the **"continuous infusion input, not an event" edge case.**

- PK: one compartment, `dxdt_TCZ_C = TCZIN - CLTCZ/VTCZ*TCZ_C`. No absorption
  compartment, no peripheral compartment, no TMDD (receptor binding is not
  modeled as its own dynamic state for this compound — IL-6 blockade is a
  plain concentration-driven ratio, `EMAXIL6 * CTCZ / (IC50IL6 + CTCZ)`).
  This matches Archetype 1 exactly.
- Dosing: `TCZIN` is **not** an `ev()`/`data_set()` event and is **not** a
  disguised `rate=`-driven infusion either. It is a zero-order input term
  computed directly inside `$ODE` from a sum of smooth rectangular window
  functions (`win(tt, start, dur)`, a logistic-product pulse of width `dur`
  days), one per monthly dose:
  ```
  double tczsum = 0.0;
  for(int k = 0; k < (int)TCZN; k++) tczsum += win(tt, TCZSTART + 28.0 * k, 1.0);
  double TCZIN = TCZON * 8.0 * WT * tczsum;
  ```
  This was confirmed by reading the scenario-running code at the bottom of
  the file (`run_scen()`, `mrgsim(end=TEND, delta=1, ...)`): there is no
  `ev()`, `data_set()`, or dosing-event table anywhere in the script. Every
  scenario is a pure `param()` call — the file's own header says as much
  ("WHY THERE ARE (ALMOST) NO DOSING EVENTS" / "Episodic IV agents ... enter
  through smooth rectangular windows"). So this is squarely the design
  guide's edge case, not an ordinary event-based bolus/infusion in disguise.
  **Preserved as-is** — `TCZIN`'s definition and the `win()` mechanism were
  not touched, and dosing is still driven purely by parameters
  (`TCZON`, `TCZSTART`, `TCZN`), never by an mrgsolve event.

## Renaming applied (tocilizumab only)

| Original | Refactored |
|---|---|
| `CLTCZ` | `CL_TCZ` |
| `VTCZ` | `V1_TCZ` |
| `TCZ_C` (compartment) | `CENT_TCZ` |
| `CTCZ` (local concentration var in `$ODE`) | `C_TCZ` |
| `EMAXIL6` | `EMAX_TCZ` |
| `IC50IL6` | `EC50_TCZ` |
| — (none; original had no Hill exponent) | `GAMMA_TCZ` (new, `= 1.0`, documented as "original had no explicit Hill term") |
| `IL6BLK` | `EFFECT_TCZ` |
| `dxdt_TCZ_C` | `dxdt_CENT_TCZ` |
| `TCZIN`, `TCZON`, `TCZSTART`, `TCZN` | unchanged (dosing driver/covariates, out of the naming convention's scope, and touching them would have broken the scenario list) |
| `CTCZr` (the `$TABLE`/`$CAPTURE` recomputation) | left named `CTCZr`, only its right-hand side updated to `CENT_TCZ / V1_TCZ`, so the `$CAPTURE`d column name is unchanged |

`EFFECT_TCZ` is written in the Hill form required by the guide:
```c
double EFFECT_TCZ = EMAX_TCZ * pow(C_TCZ, GAMMA_TCZ)
                    / (pow(EC50_TCZ, GAMMA_TCZ) + pow(C_TCZ, GAMMA_TCZ));
```
Since the original's effect term was already a plain Emax/EC50 ratio
(`EMAXIL6 * CTCZ / (IC50IL6 + CTCZ)`), this is a **rename, not a refit** —
`GAMMA_TCZ = 1.0` reproduces the original ratio exactly for any `C_TCZ`.

`EFFECT_TCZ` is used in exactly the two places `IL6BLK` was used
(both inside the B-cell/germinal-centre block, downstream of TCZ's PK):
`BHELP = TFH * (1.0 - 0.6 * EFFECT_TCZ) * (1.0 - 0.85 * COSBLKN)` and the
germinal-centre proliferation term in `dxdt_BGC`
(`... * (1.0 - 0.5 * EFFECT_TCZ) ...`). Only the `IL6BLK`/`EFFECT_TCZ` token
was touched on those two lines; every other term on them (`COSBLKN`, `MPAINH`,
`MTORINH`, `RTXKILL`, etc., belonging to other compounds) is untouched.

Note: `EMAXIL6`/`IL6BLK` never enters `NIS` (net immunosuppression) in either
`$ODE` or `$TABLE` — the original model deliberately keeps IL-6 blockade out
of the NIS scalar (it acts only on B-cell help / germinal-centre kinetics,
not on the T-cell/CNI/mTOR immunosuppression axis). This asymmetry was
preserved exactly, not "fixed" or normalized.

## Diff scope confirmation

`diff -u ktx_mrgsolve_model.R ktx_mrgsolve_model_refactored.R` touches exactly
9 hunks, all within: the `CLTCZ`/`VTCZ` param lines, the `EMAXIL6`/`IC50IL6`
param lines (+1 new `GAMMA_TCZ` line), the `TCZ_C` compartment line, the
`CTCZ`/`IL6BLK` derivation lines in `$ODE`, the `dxdt_TCZ_C` line, the two
`BHELP`/`dxdt_BGC` lines that reference `IL6BLK`, and the `CTCZr` line in
`$TABLE`. Nothing else in the ~1,013-line file differs (all other compounds'
`$PARAM`/`$CMT`/`$ODE`/`$TABLE`/`$CAPTURE` lines, all `$GLOBAL` helper
functions, and the entire scenario list and simulation-running code at the
bottom are byte-identical).

## Verification

**Method.** Both `ktx_mrgsolve_model.R`'s and `ktx_mrgsolve_model_refactored.R`'s
embedded mrgsolve code blocks were extracted, built with `mcode()`
(mrgsolve 2.0.1, `Rscript.exe` at `C:\Program Files\R\R-4.6.0\bin`), and run
through three parameter sets built from the *original file's own* scenario
list: `"01 Standard risk"` (baseline, TCZ off — sanity check that untouched
code paths are unaffected), the file's own `"22 dnDSA ABMR - anti-IL-6
(tocilizumab/clazakizumab)"` scenario (`ADHSTART=365, ADHFINAL=0.45,
ADHEND=730, TCZON=1, TCZSTART=730, TCZN=18`), and one additional stress case
turning TCZ on from day 0 for 26 monthly infusions (`TCZON=1, TCZSTART=0,
TCZN=26`) to exercise the dosing-input mechanism more thoroughly than the
file's own scenario alone (which only runs it for the last ~1.5 years of a
5-year simulation). Every `$CAPTURE`d output column was compared across the
full daily time grid (`delta=1`, `end=1825`, same solver tolerances as the
original: `atol=1e-8, rtol=1e-5`).

**Result: exact match. Maximum relative deviation observed = 0.0 (bit-identical)
across all three scenarios and every `$CAPTURE`d column.** This is the
expected outcome for Archetype 1 + the continuous-infusion edge case (pure
structural reorganization, no Hill-fitting) per the guide's tolerance table.

**A pre-existing, unrelated build blocker was found and worked around only
for this verification, not fixed in either file.** `$ODE` reuses the for-loop
iterator name `k` across four separate top-level loops (`belasum`, `tczsum`,
`fzbsum` x2, `plexsum`); mrgsolve's variable-hoisting preprocessor lifts every
`int k` declaration to the top of the generated C++ function without
deduplicating, so the generated source redeclares `int k` and fails to
compile under mrgsolve 2.0.1 / GCC 14.3.0 — for the untouched original file
just as much as for the refactored one, since the refactor copies these loops
verbatim. This has nothing to do with tocilizumab. It is logged as upstream
issue #27 in
[`translations/UPSTREAM_ISSUES.md`](../translations/UPSTREAM_ISSUES.md).
Verification here used in-memory-only copies of both models' code with the
five reused `k` names mechanically renamed to be unique (`kbe`, `ktc`, `kf1`,
`kf2`, `kpl` — rename only, no change to any bound, increment, or loop body),
applied identically to both models purely to make them buildable for
comparison. **Neither the checked-in original nor the delivered
`_refactored.R` was changed** — both still contain the original's `int k`
loops as written, and neither currently builds against mrgsolve 2.0.1 without
that same workaround. Since the workaround is identical text-for-text on both
sides of the comparison and touches only a loop counter's name (never TCZ,
never any dosing amount/timing), it cannot be masking a TCZ-specific
regression — the 0.0 deviation result stands on its own regardless of this
unrelated issue.

## Anything else worth flagging

- Surprising, but genuinely part of the original model rather than a defect:
  the `$TABLE` recomputation block never recomputes `IL6BLK`/`EFFECT_TCZ` at
  all (only `CTCZr`, the raw concentration, is captured) — the guide's
  worked archetypes don't show a compound whose Hill/effect term isn't
  itself part of `$CAPTURE`. This meant the equivalence check for TCZ's
  *effect* had to be indirect: verified through the trajectories of the
  downstream states it feeds (`BGC`, `PB`, `LLPC`, `DSA`, and everything the
  guide's own captured Banff/ABMR/`GS` outputs summarize), which is exactly
  what the 0.0-deviation scenario runs above cover — if `EFFECT_TCZ` had
  been mis-derived, `DSA`/`BANFF_CG`/`GS` at 5 years would have diverged.
- The other 415-file-corpus naming chaos the guide describes for tocilizumab
  (`TCZ_central`, `TCZC`, `A_TCZ1`, `C1_tcz`, `TOCZ_C`, ...) is corroborated
  here too: this file's own name was `TCZ_C`, a sixth spelling not in the
  guide's example list.
