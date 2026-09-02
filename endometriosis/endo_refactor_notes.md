# Endometriosis (endo) — PK/PD refactor notes

Sibling files: `endo_mrgsolve_model_refactored.R` (the rewritten model + driver
script) and `endo_mrgsolve_model_refactored.cpp` (bare mrgsolve DSL, extracted
verbatim from the quoted string in the `.R` file — confirmed byte-identical).
The original `endo_mrgsolve_model.R` is untouched.

## Compounds in scope (from `compound_perturbation_census.md`)

All five census rows for `endometriosis` were confirmed against the actual
code — all real, correctly identified drugs, standard endometriosis
pharmacotherapy, no census corrections needed:

| Compound | Original PK | Original effect term |
|---|---|---|
| Leuprolide | `Gut_leup`+`Depot_leup`+`Central_leup`, `ka_leup`/`ka_dep`/`Vd_leup`/`CL_leup`/`F_leup` | drives `GnRHR_occ` desensitization ODE (bespoke — see below) |
| Elagolix | `Gut_ela`+`Central_ela` | `ela_inh` = `Emax_GnRHR_ela * pow(C_ela,1.0) / (IC50_GnRHR_ela + C_ela)` |
| Dienogest | `Gut_die`+`Central_die` | `die_antiprol` = `Emax_prol_die * C_die / (IC50_prol_die + C_die)` |
| Letrozole | `Gut_let`+`Central_let` | `ai_inh` = `Emax_aro_let * C_let / (IC50_aro_let + C_let)` |
| Norethindrone acetate (NETA) | `Gut_neta`+`Central_neta` | `neta_bmd_prot` = `0.5 * Emax_prol_die * C_neta / (IC50_prol_die + C_neta)` (reused **dienogest's own** potency constants, scaled by 0.5 — not an independent potency) |

## What "duplicate concentration site" meant here, and how it was normalized

For all five drugs, the original computed the exposed concentration **twice
under two different names**: once in `$MAIN` (`C_leup`, `C_ela`, `C_die`,
`C_let`, `C_neta` — feeding the PD/Hill terms) and again in `$TABLE`
(`C_leup_out`, `C_ela_out`, `C_die_out`, `C_let_out` — feeding `$CAPTURE`;
NETA was never captured at all in the original). Both computations use the
identical formula, just under different names, at the two mrgsolve blocks
that need it (`$MAIN` locals are not visible inside `$TABLE`).

Normalization applied: give both sites the **same canonical name**,
`C_<STEM>`, and the same formula. This is not a behavior change — the
downstream ODE dynamics only ever read the `$MAIN` copy (never the `$TABLE`
one), exactly as in the original.

Implementation detail worth flagging for future refactors on this corpus:
the naive fix of writing `double C_LEUP = <expr>;` again inside `$TABLE`
**fails to compile** — confirmed live via the qspserver `mrgsolve_api`:

```
error: redefinition of 'double {anonymous}::C_LEUP'
```

mrgsolve hoists a `$MAIN`-declared `double X = ...` into file-scope storage
shared with `$ODE`/`$TABLE` (matches the cross-block value-sharing mechanism
`driver-patches/HANDOFF.md`'s discoverability-audit note describes for other
files in this corpus). The working fix, applied here, is a **bare
reassignment** in `$TABLE` (no `double` keyword) — `C_LEUP = <expr>;` — which
recomputes the same formula from the live, post-solve state at each
reporting time, matching what the original's separate `$TABLE` computation
produced. Same pattern applied to `EFFECT_ELA`, `EFFECT_LEUP`, `EFFECT_LET`,
`EFFECT_DIE`, `EFFECT_NETA`.

NETA's `C_NETA`/`EFFECT_NETA` are new captures (the original never exposed
NETA at all) — added per the naming convention's discoverability
requirement, not present in the original for direct before/after comparison,
but NETA's PK/effect equations are otherwise untouched math, and its
contribution is exercised (and indirectly validated) through Scenario 5's
`BMD_pct`/`Lesion_pct_change`, which do have an original counterpart and
match exactly (see Verification below).

## Archetype per compound

- **Leuprolide (LEUP)** — bespoke dual-absorption archetype: a rapid `Gut`
  path (`KA_LEUP`) plus a slow-release SC `Depot` path (`KDEP_LEUP`), both
  feeding one `Central` compartment, no peripheral compartment. Does not
  match archetypes 1–3 cleanly (two absorption inputs, not one), so it is
  disclosed as bespoke rather than forced into archetype 3's single-depot
  shape.
- **Elagolix, Dienogest, Letrozole, NETA (ELA/DIE/LET/NETA)** — archetype 3
  minus peripheral (Gut + Central only), linear elimination, matching the
  original exactly.

## The Hill interface — four rename-only, one bespoke

- **Elagolix, Letrozole, Dienogest** — the original's effect terms were
  already plain Hill ratios (`Emax*C/(IC50+C)`, with Elagolix's `pow(C,1.0)`
  making the (absent) Hill exponent explicit). Pulled out as
  `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>=1` named parameters, same values.
  A rename, not a refit.
- **NETA** — same shape, but the original computed it from **dienogest's**
  own `Emax_prol_die`/`IC50_prol_die` constants (times a hardcoded 0.5
  scale), not its own. `EMAX_NETA = 0.40` (= 0.5 x 0.80, the original's own
  arithmetic) and `EC50_NETA = 5.0` (identical value to `EC50_DIE`) are named
  params holding the exact same numbers the original's hardcoded expression
  produced — a promotion of a magic-number scale factor into a named param
  (same precedent as `abdominal-aortic-aneurysm`'s NF-kB/ROS promotion), not
  a refit.
- **Leuprolide — bespoke, disclosed.** The original does **not** model
  leuprolide's GnRH-agonist pharmacology as an algebraic Hill ratio. It
  models the real, well-known paradoxical flare/desensitization dynamic as
  an explicit ODE state, `GnRHR_occ` (GnRH-receptor activity, 0–1), driven
  by a desensitization term (`k_desens_on * leup_active`) competing against
  a first-order recovery term (`k_desens_off * (GnRHR_base - GnRHR_occ)`) —
  the same state that Elagolix's antagonism (`gnrhr_antag`) and the OCP flag
  also perturb. Per the guide's explicit instruction for this exact
  compound, that ODE structure is preserved unchanged (renamed variables
  only); `EFFECT_LEUP` is a disclosed bespoke term — a **gated
  concentration pass-through** (`C_LEUP` if `use_leup>0 && C_LEUP>0.01`,
  else 0), not a Hill ratio — that drives the shared `GnRHR_occ` ODE exactly
  as the original's `leup_active` did. No `EC50_LEUP`/`EMAX_LEUP` was
  invented: leuprolide's action here is a linear rate-constant term gated by
  a threshold, not a saturating dose-response curve, so it will not be
  discoverable by the `C_<STEM>`+`EC50_<STEM>` grep pattern — the same,
  correct outcome the guide describes for a construct with no meaningful
  EC50 (its own gene-therapy example). Disclosed rather than forced.

## Naming convention applied

`GUT_/DEPOT_/CENT_` intentionally kept as the original's own
`Gut_/Depot_/Central_` capitalization pattern for compartment names (not
uppercased to `GUT_LEUP` etc.) — matching the shared rule to preserve what
verification needs to line up 1:1 with the original's own compartment names
column-for-column; only the PK **parameters** (`KA_LEUP`, `KDEP_LEUP`,
`V1_LEUP`, `CL_LEUP`, `F_LEUP`, `KA_ELA`, `V1_ELA`, `CL_ELA`, `F_ELA`, ... one
set per stem) and the Hill/effect terms (`C_<STEM>`, `EMAX_<STEM>`,
`EC50_<STEM>`, `GAMMA_<STEM>`, `EFFECT_<STEM>`) were renamed to the fork-wide
convention. Disease-state compartments (`GnRHR_occ`, `FSH_plasma`,
`LH_plasma`, `E2_plasma`, `P4_plasma`, `Lesion`, `IL6_peritoneal`,
`PGE2_local`, `NGF_lesion`, `Pain_score`, `BMD_lumbar`, `Aromatase_act`,
`E2_local`) and the OCP regimen flags (`use_ocp`, `ocp_supp`,
`ocp_antiprol` — not one of the five compounds in scope, a combined-pill
regimen toggle, not an independently-dosed drug) are untouched, matching the
`abdominal-aortic-aneurysm` precedent of leaving non-PK physiological states
and out-of-scope flags alone.

## A note on the "$PARAM with =0 default" line in the guide

`FORK_WORKFLOW_GUIDE.md`'s "qspserver compatibility requirements" section
(point 2) reads literally as requiring every `C_<STEM>`/`EFFECT_<STEM>` to
be declared as a `$PARAM` with a `=0` default. That is **not** how the
actual, already-verified 114-file corpus does it (checked
`abdominal-aortic-aneurysm/aaa_mrgsolve_model_refactored.R`: `C_DOXY`,
`C_STAT`, `C_BB`, `EFFECT_*` are all plain `$ODE` doubles, never `$PARAM`),
and it directly contradicts the corpus-wide discoverability-audit fix
documented in `driver-patches/HANDOFF.md` (commit `2a93260`), which is the
actual, tested, audited contract: a single contiguous
`double C_<STEM> = <expr>;` (or, per this file's finding, a bare
reassignment sharing a `$MAIN`-hoisted declaration) plus a same-stem
`EC50_<STEM>` parameter anywhere in the file. This refactor follows the
audited, established practice (matches 114 prior files) rather than the
literal text of that one line, which appears to be aspirational/superseded
and unimplemented anywhere in the actual corpus.

## Upstream defects

None found. The original `endo_mrgsolve_model.R` compiles and runs cleanly
under mrgsolve 2.0.1 (confirmed via `/model_manifest` and `/run_simulation`
against the qspserver `mrgsolve_api`) — no syntax-only fix was needed for
either the original or the refactored sibling. No new `UPSTREAM_ISSUES.md`
entry required.

## Verification

Ran all 7 of the original file's own dosing scenarios (Scenario 0: no
treatment; 1: leuprolide depot 3.75 mg q28d; 2: elagolix 150 mg/d; 3:
elagolix 200 mg BID; 4: dienogest 2 mg/d; 5: letrozole 2.5 mg/d + NETA 5
mg/d add-back; 6: combined OCP cyclic dienogest) through both the original
and the refactored model via the qspserver `mrgsolve_api`
(`POST /model_manifest` then `POST /run_simulation`, requests spaced ~2 s
apart, one at a time), with identical dosing (matched NM-TRAN
`amt`/`cmt`/`ii`/`addl`/`time` records; `cmt` compartment indices verified
separately per model via each model's own `/model_manifest`
`outputPaths` order — the refactored file's `$CMT` block lists
`Depot_leup` before `Central_leup` where the original has the reverse
order, so `Depot_leup`'s dosing index is 2 in the refactored model vs. 3 in
the original; every other compartment index is identical between the two
files).

**Duration shortened to 1 year (8760 h, daily reporting, `delta=24`) from
the original's 5-year duration** — same precedent as other refactors in
this corpus needing to shorten a long-running scenario for solver/payload
size (e.g. `acne-vulgaris`'s 34-week scenario shortened to its own dosing
window). All five drugs' PK reach steady state and the HPO axis
re-equilibrates well within this window; the dosing regimen itself (route,
amount, interval) is identical to the original's own scenario definitions,
only total simulated duration differs.

**Result: exact match (max relative difference = 0.0) across every
comparable output, every scenario, every timepoint** (366–380 points
depending on scenario) — all 24 disease/PK state compartments
(`Gut_leup`/`Depot_leup`/`Central_leup`/... through `E2_local`) plus every
original `$CAPTURE`'d quantity (`C_leup_ng`↔`C_LEUP`, `C_ela_ng`↔`C_ELA`,
`C_die_ng`↔`C_DIE`, `C_let_ng`↔`C_LET`, and `GnRHR`, `E2_pct_suppress`,
`Lesion_pct_change`, `Pain_delta`, `BMD_pct` — all four of which share
identical names/formulas between the two files). This is exactly the
"near-exact match" the guide expects for pure structural reorganization
(archetypes 1–3, no Hill-fitting) — floating-point-scale deviation would
have been acceptable; 0.0 is what was actually observed, consistent with a
pure rename/reorganization introducing no numeric change.

NETA's new `C_NETA`/`EFFECT_NETA` captures have no original counterpart to
diff against directly (the original never captured NETA at all), but
Scenario 5 exercises NETA's PK+effect equations and its two downstream
disease outputs (`BMD_pct`, `Lesion_pct_change`) matched the original
exactly — a wrong NETA formula would have shown up there. `C_NETA`/
`EFFECT_NETA` were also spot-checked directly in the refactored output
(non-zero, plausible values once dosing starts, e.g. `C_NETA` rising to a
~3.1 ng/mL plateau, `EFFECT_NETA` to ~0.154) — consistent with NETA reaching
its expected low-dose add-back steady state.

Leuprolide's bespoke `EFFECT_LEUP`/`GnRHR_occ` dynamic was also spot-checked
directly: `GnRHR_occ` collapses to its 0.01 floor within the first reporting
interval under continuous agonist exposure (matching the original's
own aggressive `k_desens_on`-driven desensitization) — and this trajectory
is one of the exactly-matching captured columns (`GnRHR`) above, so the
bespoke ODE reproduces the original's flare/suppress dynamic bit-for-bit,
not just qualitatively.

`Rscript -e 'parse("endo_mrgsolve_model_refactored.R")'` succeeds with no
error. The `.cpp` extraction is confirmed byte-identical to the quoted
string in the `.R` file.

## Summary for the compound census

- Dienogest, Elagolix, Letrozole, Leuprolide, NETA: all confirmed real,
  correctly classified drugs; census method label ("Normalize duplicate
  concentration sites, then redirect") accurately describes what all five
  needed. No compounds the census missed; no misclassifications found.
