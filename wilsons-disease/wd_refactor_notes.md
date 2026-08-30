# Refactor notes — `wilsons-disease/wd_mrgsolve_model.R`

Scope: all four compounds carrying rows in
`driver-patches/data/compound_perturbation_census.md`: **D-Penicillamine
(DPA)**, **Zinc (ZN)**, **Trientine (TRI)**, **ALXN1840 / bis-choline
tetrathiomolybdate (TTM)**. All four are refactored — the original models
each purely as an externally-dosed drug (see the "Is Zinc really a drug
here?" section below for the check the task specifically asked for).

## The original compiles cleanly — no build-compat fix needed

Confirmed via `POST /model_manifest` on the untouched original's own DSL
(extracted verbatim from `wd_code <- '...'`): HTTP 200, all 37 outputs and
all `$PARAM @annotated` parameters listed with no error. Unlike ~20 other
files in this corpus, `wd_mrgsolve_model.R` has no `$CMT`+`$INIT`
redeclaration, no `$CAPTURE`/compartment collision, and no deprecated
idiom — it uses `$INIT @annotated` alone (no separate `$CMT` block) and
`$TABLE`-local `capture NAME = expr;` statements, both of which mrgsolve
2.0.1 accepts as written. **No entry was added to
`translations/UPSTREAM_ISSUES.md`** — there is no defect to log.

## Is Zinc really a drug here? (the task's specific check)

Checked explicitly, since zinc is genuinely dual-purpose biologically (an
endogenous trace mineral *and* a chelation therapy dosed as zinc
acetate/sulfate). In this model it is unambiguously the latter and only
the latter:

- Zinc has its own `GUT_ZN`/`CENT_ZN` depot+central PK compartments,
  its own absorption/clearance parameters (`ka_Zn`, `F_Zn`, `Vd_Zn`,
  `CL_Zn`), and its own `Zn_dose_mg` bookkeeping parameter — structurally
  identical in kind to DPA/TRI/TTM's own PK blocks.
- The original's own R driver doses it exactly like the other three
  drugs: `dose_ZN <- ev(amt = 50 * 0.15, cmt = "GUT_ZN", ii = 8, addl =
  ..., time = 0)` (Scenario 3, "Zinc Acetate 50 mg TID (Maintenance)").
- There is **no separate "endogenous/dietary zinc" state** anywhere in
  the model — no baseline zinc pool, no zinc homeostasis compartment, no
  zinc-deficiency dynamic. The only zinc-related state in the whole
  24-compartment system is this therapeutic PK pair.

So Zinc is a drivable external dose in this model, not a pure disease-state
variable — refactored as the fourth compound, same as the other three.
The census row is confirmed correct as classified ("Redirect concentration
(clean single site)"), not corrected.

## Archetypes

All four compounds share the same PK architecture in the original:
**Archetype 3 without a peripheral compartment** (depot + central only,
first-order absorption, first-order elimination) — `GUT_<STEM>` feeding
`CENT_<STEM>` via `ka_<STEM>`, `CENT_<STEM>` cleared via `CL_<STEM>` over
a fixed total volume (`Vd_<STEM>` per-kg × 70 kg). No compound has a
peripheral compartment or TMDD-style receptor binding in the original.

Convenient starting point: **the original's compartment names already
matched the guide's naming convention exactly** — `GUT_DPA`, `CENT_DPA`,
`GUT_ZN`, `CENT_ZN`, `GUT_TRI`, `CENT_TRI`, `GUT_TTM`, `CENT_TTM` — so
**no compartment was renamed**. Only the governing *parameters* and the
*local variable names* inside `$ODE` needed renaming to the convention:

| Original | Refactored | Role |
|---|---|---|
| `ka_DPA` | `KA_DPA` | absorption rate |
| `Vd_DPA` (0.5 L/kg × 70) | `V1_DPA` = 35.0 | central volume (pre-computed total, same number) |
| `CL_DPA` | `CL_DPA` | clearance (already matched) |
| `Kchel_DPA` | `KCHEL_DPA` | chelation potency (casing only, not a guide-named role) |
| `Cc_DPA` | `C_DPA` | **the exposed concentration** |
| `DPA_chel` | `EFFECT_DPA` | **the compound's effect on disease** |
| `ka_Zn` | `KA_ZN` | absorption rate |
| `Vd_Zn` (0.25 L/kg × 70) | `V1_ZN` = 17.5 | central volume |
| `CL_Zn` | `CL_ZN` | clearance |
| `IC50_Zn` | `EC50_ZN` | Hill EC50 |
| `Cc_ZN` | `C_ZN` | **the exposed concentration** |
| `Zn_eff` | `EFFECT_ZN` | **the compound's effect on disease** |
| `ka_TRI` | `KA_TRI` | absorption rate |
| `Vd_TRI` (0.4 L/kg × 70) | `V1_TRI` = 28.0 | central volume |
| `CL_TRI` | `CL_TRI` | clearance (already matched) |
| `Kchel_TRI` | `KCHEL_TRI` | chelation potency (casing only) |
| `Cc_TRI` | `C_TRI` | **the exposed concentration** |
| `TRI_chel` | `EFFECT_TRI` | **the compound's effect on disease** |
| `ka_TTM` | `KA_TTM` | absorption rate |
| `Vd_TTM` (1.0 L/kg × 70) | `V1_TTM` = 70.0 | central volume |
| `CL_TTM` | `CL_TTM` | clearance (already matched) |
| `Emax_TTM` | `EMAX_TTM` | Hill Emax |
| `EC50_TTM` | `EC50_TTM` | Hill EC50 (already matched) |
| `Cc_TTM` | `C_TTM` | **the exposed concentration** |
| `TTM_eff` | `EFFECT_TTM` | **the compound's effect on disease** |

`F_DPA`/`F_Zn→F_ZN`/`F_TRI`/`F_TTM` (bioavailability) and
`DPA_dose_mg`/`Zn_dose_mg→ZN_dose_mg`/`TRI_dose_mg`/`TTM_dose_mg` are
declared in `$PARAM` but **never referenced anywhere in the original's
`$ODE`** (confirmed by exact-name grep across the whole file — each
appears only in its own declaration line and in the R-side
`modifyList()` scenario-parameter lists, never inside the quoted DSL
block's `$ODE`). Bioavailability is instead applied by the R driver
**pre-scaling the administered dose amount** before building each `ev()`
(e.g. `amt = 500 * 0.55` for DPA), not via mrgsolve's native per-compartment
`F_<CMT>` mechanism. Preserved exactly as declared-but-unused, same
values, same casing convention applied to the `ZN` stem for consistency
with its now-uppercase PK parameters.

## The Hill interface — two different shapes, both preserved faithfully

### Zinc (ZN) and ALXN1840 (TTM): already exactly Hill-shaped — rename, not a fit

- `Zn_eff = Cc_ZN / (IC50_Zn + Cc_ZN)` → `EFFECT_ZN = EMAX_ZN * pow(C_ZN,
  GAMMA_ZN) / (pow(EC50_ZN, GAMMA_ZN) + pow(C_ZN, GAMMA_ZN))`. The
  original has **no explicit Emax term** — its formula's ceiling is
  implicitly 1 as `C_ZN -> Inf`. `EMAX_ZN = 1.0` is a **new** parameter,
  but it is *math-implied*, not invented: it names a value the original's
  own formula already produces at its limit, the same treatment the guide
  explicitly sanctions for `GAMMA_<STEM> = 1` when no Hill coefficient was
  declared. `GAMMA_ZN = 1.0` is likewise new (no explicit Hill coefficient
  in the original).
- `TTM_eff = Emax_TTM * Cc_TTM / (EC50_TTM + Cc_TTM)` → `EFFECT_TTM =
  EMAX_TTM * pow(C_TTM, GAMMA_TTM) / (pow(EC50_TTM, GAMMA_TTM) +
  pow(C_TTM, GAMMA_TTM))`. `EMAX_TTM`/`EC50_TTM` already existed in the
  original (renamed casing only for TTM's); `GAMMA_TTM = 1.0` is new (no
  explicit Hill coefficient in the original), same treatment as ZN.

### D-Penicillamine (DPA) and Trientine (TRI): genuinely linear, not Hill-shaped — preserved as linear, not force-fit

The original's DPA and Trientine "effect" terms —
`DPA_chel = Kchel_DPA * Cc_DPA` and `TRI_chel = Kchel_TRI * Cc_TRI` — are
**first-order (linear) in the drug's own concentration, with no
saturation on the drug side at all**. The only saturation anywhere near
these terms (e.g. `* CU_HEP / (CU_HEP + 10.0)` in `total_chelation`, or
`* CU_NCBC / (CU_NCBC + 5.0)` in `NCBC_chel_out`) is on the *substrate*
(hepatic or serum copper pool), not on drug concentration — that
saturation belongs to the disease side of the equation, not to DPA/TRI's
own PK-to-effect interface, and is untouched by this refactor.

Per the guide's Hill-interface section, forcing
`Emax*C^gamma/(EC50^gamma+C^gamma)` onto a term that is actually linear
and unbounded in `C` would either invent a false ceiling (misrepresenting
kinetics the original never modeled) or require inventing an `EC50` that
doesn't correspond to anything in the source — both of which the guide's
"don't force a fit" principle rules out. Instead, `EFFECT_DPA = KCHEL_DPA
* C_DPA` and `EFFECT_TRI = KCHEL_TRI * C_TRI` are exposed exactly as the
original computed them: a single named quantity per compound
(satisfying "the compound's effect on disease is expressed as one named
function of concentration"), just not Hill-shaped, because the original
isn't. Each is used at multiple downstream disease-equation sites — the
hepatic chelation flux, the NCBC chelation flux, the urinary copper
excretion rate, fibrosis regression, and (DPA only) neurological
recovery — combined only at each of those points, never pre-combined into
a shared cross-drug term (both `EFFECT_DPA` and `EFFECT_TRI` stay fully
independently driveable throughout).

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>` (the guide's "= 0 default" instruction)

Tried first, as the guide's qspserver-compatibility section literally
describes: declaring `C_DPA = 0` etc. directly in `$PARAM` and assigning
the recomputed value in `$ODE`. This does not compile under mrgsolve
2.0.1 — `$PARAM` values are injected into `$ODE` as **read-only
references**, so `double C_DPA = CENT_DPA / V1_DPA;` immediately after a
`$PARAM C_DPA = 0` declaration is a redeclaration/read-only-assignment
error. This is the same constraint already documented in
`pagets-disease/pbd_refactor_notes.md` and reported there as shared by
every prior model in this refactor batch (`rheumatoid-arthritis`,
`sepsis`, `thyroid-eye-disease`, `polymyalgia-rheumatica`,
`kidney-transplant-rejection`, and later scale-out batches) — this file
follows that same, actually-compiling precedent instead: `C_DPA`,
`EFFECT_DPA`, `C_ZN`, `EFFECT_ZN`, `C_TRI`, `EFFECT_TRI`, `C_TTM`,
`EFFECT_TTM` are plain `double` locals computed in `$ODE`, then listed in
a new `$CAPTURE` block. Confirmed via `POST /model_manifest` on the
refactored DSL: all eight appear in `outputPaths` (see Verification
below) — fully discoverable through the output side of the manifest, just
not the overridable `parameters` side, exactly as the pagets precedent
describes.

## Verification

Per the guide's mandatory protocol: extracted the bare DSL text from both
`wd_mrgsolve_model.R` (`wd_code <- '...'`, unmodified) and
`wd_mrgsolve_model_refactored.R` (same variable, refactored), and ran all
four of the original file's own single-compound dosing scenarios through
both, via the qspserver `mrgsolve_api` (`POST /model_manifest` once per
model, then `POST /run_simulation` per scenario/model pair), requests
spaced ~2.2s apart per the shared-service note:

1. **S2 — D-Penicillamine**: `dose_DPA` (`amt = 500*0.55 = 275` into
   `GUT_DPA`, `ii = 8`, TID).
2. **S3 — Zinc**: `dose_ZN` (`amt = 50*0.15 = 7.5` into `GUT_ZN`, `ii = 8`,
   TID).
3. **S4 — Trientine**: `dose_TRI` (`amt = 500*0.45 = 225` into `GUT_TRI`,
   `ii = 8`, TID).
4. **S5 — ALXN1840**: `dose_TTM` (`amt = 15*0.30 = 4.5` into `GUT_TTM`,
   `ii = 24`, QD — the ATLAS regimen).

**Shortened window (disclosed per the guide's solver-budget allowance):**
the original's own scenarios run for 5 years (`sim_hrs = 43800`) at daily
sampling. Verification instead used a 60-day window (`end = 1440`,
`delta = 24`, matching the original's own `by = 24` daily sampling
resolution) — long enough to fully exercise each compound's own
absorption/distribution/clearance to steady oscillation and to drive
every downstream disease-equation term it feeds (hepatic/NCBC chelation,
urinary excretion, fibrosis regression, neuro effects), not a
truncated/degenerate run. This is a concurrency/runtime consideration for
a shared service, not evidence of any solver failure — no scenario timed
out or errored at the full 5-year horizon either; the window was
shortened purely to keep the comparison lightweight, matching the same
allowance exercised in other refactor notes in this batch.

**Result: exact match, max abs diff = 0.0, every shared output, every
scenario**, compared point-by-point across the full 62-point time grid
(`Cc_DPA_out, Cc_ZN_out, Cc_TRI_out, Cc_TTM_out, Cu_hep_out, NCBC_out,
Cp_out, Cu_urine_out, Cu_brain_out, Cu_corn_out, ROS_out, ALT_out,
Fib_out, Neuro_out, Cu_serum_total, KF_index, UWDRS_approx`). This is the
expected outcome for a pure structural reorganization (all four
compounds are renames, not fits) — confirmed by the Python comparison
script (no tolerance loosening; genuinely `0.0` for every point).

`/model_manifest` on the refactored DSL additionally confirmed
`C_DPA`/`EFFECT_DPA`/`C_ZN`/`EFFECT_ZN`/`C_TRI`/`EFFECT_TRI`/`C_TTM`/
`EFFECT_TTM` all appear in `outputPaths`, and spot-checked non-degenerate,
correctly-shaped values in the S2/S3/S5 runs, e.g.:

- DPA (S2, t=48h): `C_DPA = 6.6931`, `EFFECT_DPA = 2.0079` =
  `KCHEL_DPA(0.3) * 6.6931` — exact linear relationship, as expected.
- ZN (S3, t=48h): `C_ZN = 0.7193`, `EFFECT_ZN = 0.2645` =
  `EMAX_ZN(1.0) * 0.7193 / (EC50_ZN(2.0) + 0.7193)` — exact Hill
  relationship, as expected.
- TTM (S5, t=72h): `C_TTM = 0.0915`, `EFFECT_TTM = 0.1516` =
  `EMAX_TTM(0.98) * 0.0915 / (EC50_TTM(0.5) + 0.0915)` — exact Hill
  relationship, as expected.
- `C_DPA`/`C_ZN`/`C_TRI`/`C_TTM` are exact pass-throughs of
  `Cc_DPA_out`/`Cc_ZN_out`/`Cc_TRI_out`/`Cc_TTM_out` at every timepoint
  (both now derive from the same expression), confirming the rename
  introduced no numeric drift.

No scratch/debug artifacts were left in the repository — the DSL
extraction and comparison scripts ran entirely from a session-local
scratchpad directory outside the repo tree.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, all
four `wilsons-disease` rows (ALXN1840/TTM, D-Penicillamine/DPA, Trientine,
Zinc/ZN).
