# Refactor notes — `graves-disease/gd_mrgsolve_model.R`

## Scope: three of the census's seven rows are genuine drugs; four are endogenous disease-state variables

This file carries seven rows in
`driver-patches/data/compound_perturbation_census.md`: **MMI**, **PROP**,
**PTU**, **T3**, **T4**, **TRAB**, **TSH**. Checked each against the actual
code (which compartments/parameters exist, whether anything has a dosing
route) rather than trusting the census's compound name alone, per the task's
specific flag:

| Row | Real drug? | Evidence |
|---|---|---|
| MMI (methimazole) | **Yes** | Own compartment `MMI` (dosed via `ev(cmt=1,...)` in scenarios 2, 4, 5, 6), own elimination parameter `kel_MMI`, own `IC50_MMI` |
| PROP (propranolol) | **Yes** | Own compartment `PROP_C` (dosed via `ev(cmt=5,...)` in scenario 4), own elimination `kel_PROP`, own `EC50_PROP_HR` |
| PTU (propylthiouracil) | **Yes** | Own compartment `PTU_C`, own elimination `kel_PTU`, own `IC50_PTU`/`IC50_PTU_D1` — has a full, ready-to-dose PK block even though none of the file's own six shipped scenarios happens to dose it (see below) |
| T3 | **No — endogenous.** Total serum T3, `$CMT T3_C`. No depot, no absorption/bioavailability parameter, no `ev()`/dataset dosing route anywhere in the file. Produced entirely by the model's own `T3_synth` term (`ksyn_T3 * synth_factor`) plus T4→T3 deiodination (`D1_conv`); consumed only by downstream disease equations (`fT3_C`, TSH feedback, bone, heart rate). This is the disease's own thyroid-hormone output, not something administered. |
| T4 | **No — endogenous.** Same reasoning as T3 (`$CMT T4_C`, no dosing route, produced by `T4_synth`). The file's own header comment claims a separate "C5 = Levothyroxine (exogenous T4) pool" and a `USE_LT4` flag exists, but there is **no such compartment anywhere in `$CMT`** — `USE_LT4` is declared in `$PARAM` and set by one scenario's `param()` call, but never referenced in `$ODE` at all (confirmed by exact-name grep), and that scenario's own `br_lt4` dosing event actually targets `cmt=5` (`PROP_C`, propranolol's compartment), not T4. So there is no working exogenous-T4 dosing pathway in this model at all; `T4_C` itself is purely the endogenous total-T4 pool. Logged as `UPSTREAM_ISSUES.md` #94 (below), not fixed here since LT4 is not one of this file's census rows. |
| TRAB (TRAb) | **No — endogenous.** `$CMT TRAb_C`, the disease's own stimulating autoantibody. No dosing route; produced solely by `kprod_TRAb * Bcell_C`. This is the autoimmune disease state itself (analogous to Wilsons-disease's copper pools or a disease's own cytokine level), not a drug. |
| TSH (thyroid-stimulating hormone) | **No — endogenous.** `$CMT TSH_C`, pituitary output of the HPT axis. No dosing route; produced by `TSH_prod = ksyn_TSH/fb_tot`, a pure negative-feedback term. |

**Conclusion, matching the task's own flagged suspicion exactly:** only MMI,
PROP, and PTU are genuine externally-dosed drugs. T3, T4, TRAb, and TSH are
disease-state variables the drugs act *through*, not drugs themselves — same
distinction as `wilsons-disease`'s copper pools or `neonatal-
hyperbilirubinemia`'s bilirubin pool (see those files' own refactor notes for
the same style of check). The census has been corrected accordingly (see
"Census" below) rather than force-fitting a PK block onto four variables that
have no dosing mechanism at all.

**PTU is a real drug with no exercising scenario.** Unlike T3/T4/TRAb/TSH,
PTU has a complete, structurally identical PK block to MMI (compartment,
elimination rate, two IC50s used in real downstream disease equations) — it
is simply never dosed by any of the six scenarios the R script ships. It is
refactored the same as MMI/PROP; verification for it necessarily uses one
additional, explicitly-disclosed non-original dosing scenario since the file
itself never exercises it (see Verification below).

## The original does not compile — build defects fixed, disclosed, logged

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `graves_model <- '...'`) returned HTTP 500:

```
invalid class "mrgmod" object: 1: reserved symbols in model names: F_MMI
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE:
MMI,PTU_C,RAI_THYR,PROP_C,TPO_inh,ThyMass,GO_act,BoneResor
```

Logged as `translations/UPSTREAM_ISSUES.md` **#93**. Two independent defects:

1. **`F_MMI` collides with mrgsolve's reserved `F_<compartment>` symbol**
   because a compartment is literally named `MMI`. `F_MMI` is declared but
   never referenced in `$ODE` anywhere (confirmed by exact-name grep), so
   this fixes itself as a side effect of this refactor's own naming-
   convention rename: the compartment becomes `CENT_MMI`, no longer
   colliding with the parameter name `F_MMI`. No separate change needed
   beyond the rename the guide already calls for.
2. **`$CAPTURE` re-lists eight `$CMT` compartment names directly**
   (`MMI PTU_C RAI_THYR PROP_C TPO_inh ThyMass GO_act BoneResor`), which
   mrgsolve 2.0.1 rejects — compartments are already exposed as output
   columns automatically, with no `$CAPTURE` entry needed. Fixed by simply
   dropping that line; nothing is lost from the output (all eight remain
   ordinary `$CMT` columns), and the line is replaced with the refactor's
   own new `C_<STEM>`/`EFFECT_<STEM>` capture entries (see qspserver
   discoverability below).

Both fixes are syntax-only and non-numeric, confirmed by the verification
below: a **patched-original** scratch copy (original DSL, unmodified except
for the same two syntax fixes — renaming the *unused* `F_MMI` parameter to
`Fbio_MMI_unused` rather than renaming any compartment, and dropping the
same `$CAPTURE` line) was used as the verification baseline, per the guide's
"apply the same syntax-only fix... directly into the `_refactored.R`" policy.
The patched-original scratch copy itself was not committed or left in the
repo — only used in-memory during verification, then discarded.

## A second, separate defect found while refactoring (logged, not fixed — out of scope)

While tracing every `ev(cmt=...)` dosing call to confirm compartment
identity (needed to build the verification requests below), found that two
of the file's six scenarios dose a flag-only intervention into a *different*
compound's own PK compartment — the same class of bug as `#90`
(huntingtons-disease branaplam-into-tominersen). Logged as
`UPSTREAM_ISSUES.md` **#94**:

- **Scenario 5** ("Block-and-Replace"): `br_lt4 <- ev(amt=0.10, ii=1.0,
  addl=730-1, cmt=5, time=90)` is commented "LT4 start at day 90", but
  `cmt=5` is `PROP_C`/`CENT_PROP` (propranolol), not a dedicated LT4 pool —
  this model has no LT4 compartment at all (see the T4 row above).
- **Scenario 6** ("Rituximab + MMI"): `rtx_events <- ev(amt=1, cmt=1,
  time=0) + ev(amt=1, cmt=1, time=14)` doses rituximab — a bare on/off flag
  in this model, `RTX_eff = (USE_RTX>0.5) ? RTX_Bcell : 0.0`, with no PK
  compartment of its own — into `cmt=1`, i.e. `MMI`/`CENT_MMI`.

Neither LT4 nor rituximab is one of this file's three census-row drugs, so
**not fixed here** — both scenarios are preserved byte-for-behavior-identical
in the refactored sibling (same `cmt=` targets, same amounts/timing). The
Scenario 6 verification below explicitly confirms this spurious `CENT_MMI`
contamination reproduces identically between the original and the refactored
DSL (not merely "not broken by the refactor" — literally the same numbers).

## Archetype: all three drugs are Archetype 1, with one shared bespoke deviation

MMI, PTU, and propranolol are all **Archetype 1 (no depot, single
compartment, linear elimination)** — no depot/absorption compartment exists
for any of the three; each is dosed by a single `ev()` bolus directly into
its own one-compartment PK state.

**Bespoke deviation, disclosed (same reasoning as
`neonatal-hyperbilirubinemia`'s `CENT_SNMP`/`CENT_PB`/`CENT_UDCA`):** the
guide's literal Archetype-1 template computes `C_<STEM> = CENT_<STEM> /
V1_<STEM>` (a volume division). This file's own three drug compartments
**never divide by their declared volumes at all** — `Vd_MMI`/`Vd_PTU`/
`Vd_PROP` are declared in `$PARAM` but never referenced anywhere in `$ODE`
(confirmed by exact-name grep, same as `ka_<STEM>`/`F_<STEM>` for all three
compounds — none of `ka_MMI`, `F_MMI`, `ka_PTU`, `F_PTU`, `ka_PROP`,
`F_PROP` is used either). The compartment itself already **is** the
concentration the original's own downstream Hill-shaped effect terms read
directly (`MMI_conc = MMI`, `PTU_conc = PTU_C`, `PROP_C` used bare). Dividing
by a volume here would not be a rename — it would change the numeric
trajectory of the drug's own decay curve (a bolus `amt=0.10` would become a
concentration of `0.10/V1` instead of `0.10`). So `C_MMI`/`C_PTU`/`C_PROP`
are defined as a plain **identity** of their compartment
(`double C_MMI = CENT_MMI;`), not a division, and `dxdt_CENT_<STEM> =
-CL_<STEM> * C_<STEM>` reproduces `dxdt_<orig> = -kel_<STEM> * <orig>`
exactly. `KA_<STEM>`, `V1_<STEM>` (from `Vd_<STEM>`), and `F_<STEM>` are kept
as declared-but-unused parameters, same values, same treatment as the
`wilsons-disease` precedent for genuinely dead parameters.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `MMI` (cmt) | `CENT_MMI` | -- | central, concentration-state (bespoke, see above) — also removes the `F_MMI` reserved-symbol collision |
| `ka_MMI` | `KA_MMI` | 8.0 | absorption rate [declared, unused in original — preserved as-is] |
| `kel_MMI` | `CL_MMI` | 4.0 | elimination rate constant |
| `Vd_MMI` | `V1_MMI` | 0.5 | volume [declared, unused in original — preserved as-is] |
| `F_MMI` | `F_MMI` | 0.93 | bioavailability [declared, unused in original — preserved as-is; no longer collides with the compartment name now that the compartment is `CENT_MMI`] |
| `IC50_MMI` | `EC50_MMI` | 0.15 | Hill EC50 (shared by both of MMI's downstream uses) |
| -- (none) | `EMAX_MMI` (new) | 1.0 | Hill ceiling [math-implied: MMI's own ratio saturates at 1 as `C→∞`; each downstream use (`TPO_inh_max`, `MMI_Bcell`) applies its own real ceiling on top — same treatment as `EMAX_ZN` in `wilsons-disease`] |
| -- (none) | `GAMMA_MMI` (new) | 1.0 | Hill coefficient [original had no explicit Hill exponent] |
| `PTU_C` (cmt) | `CENT_PTU` | -- | central, concentration-state (bespoke) |
| `ka_PTU` | `KA_PTU` | 5.0 | [declared, unused — preserved as-is] |
| `kel_PTU` | `CL_PTU` | 9.24 | elimination rate constant |
| `Vd_PTU` | `V1_PTU` | 0.5 | [declared, unused — preserved as-is] |
| `F_PTU` | `F_PTU` | 0.75 | [declared, unused — preserved as-is] |
| `IC50_PTU` | `EC50_PTU` | 1.0 | Hill EC50 — TPO-inhibition mechanism |
| -- (none) | `EMAX_PTU` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_PTU` (new) | 1.0 | Hill coefficient [new] |
| `IC50_PTU_D1` | `EC50_PTU_D1` | 5.0 | Hill EC50 — **separate** D1-deiodinase-inhibition mechanism (see Hill interface below) |
| -- (none) | `EMAX_PTU_D1` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_PTU_D1` (new) | 1.0 | Hill coefficient [new] |
| `PROP_C` (cmt) | `CENT_PROP` | -- | central, concentration-state (bespoke) |
| `ka_PROP` | `KA_PROP` | 6.0 | [declared, unused — preserved as-is] |
| `kel_PROP` | `CL_PROP` | 4.16 | elimination rate constant |
| `Vd_PROP` | `V1_PROP` | 3.9 | [declared, unused — preserved as-is] |
| `F_PROP` | `F_PROP` | 0.26 | [declared, unused — preserved as-is] |
| `EC50_PROP_HR` | `EC50_PROP` | 0.02 | Hill EC50 (renamed to drop the `_HR` suffix — the original already shared this one EC50 across *both* of propranolol's downstream uses, D1 inhibition and HR reduction, so `_HR` was misleadingly narrow) |
| -- (none) | `EMAX_PROP` (new) | 1.0 | Hill ceiling [math-implied] |
| -- (none) | `GAMMA_PROP` (new) | 1.0 | Hill coefficient [new] |

`T3_C`, `T4_C`, `TRAb_C`, `TSH_C`, and every parameter/equation touching only
those four (HPT-axis, B-cell/TRAb dynamics, thyroid-mass/TPO-synthesis
machinery) are **completely untouched** — not renamed, not restructured.

## The Hill interface — one drug (PTU) genuinely needs two, disclosed as bespoke

**MMI**: used at two disease sites in the original with the *same* `IC50_MMI`
but two different scaling factors (`TPO_inh_max` for TPO inhibition,
`MMI_Bcell` for B-cell suppression) — a single shared occupancy fraction,
scaled differently downstream. This fits the guide's interface cleanly as one
named term:

```
EFFECT_MMI = EMAX_MMI * pow(C_MMI, GAMMA_MMI) / (pow(EC50_MMI, GAMMA_MMI) + pow(C_MMI, GAMMA_MMI));
...
TPO_inh_MMI  = TPO_inh_max * EFFECT_MMI;   // was: TPO_inh_max * MMI_conc/(IC50_MMI+MMI_conc)
MMI_immu_eff = MMI_Bcell   * EFFECT_MMI;   // was: MMI_Bcell   * MMI_conc/(IC50_MMI+MMI_conc)
```

Exact rename (`EMAX_MMI=GAMMA_MMI=1` reproduces the original ratio exactly),
not a fit.

**Propranolol**: same shape — both downstream uses (D1 inhibition, HR
reduction) already shared the *same* `EC50_PROP_HR` in the original, just
scaled differently (`* 0.3` inline for D1, `* HR_max_stim` for HR). One
named term, exact rename:

```
EFFECT_PROP = EMAX_PROP * pow(C_PROP, GAMMA_PROP) / (pow(EC50_PROP, GAMMA_PROP) + pow(C_PROP, GAMMA_PROP));
D1_inh_PROP  = EFFECT_PROP;                          // scaled by *0.3 downstream, unchanged
PROP_HR_inh  = HR_max_stim * EFFECT_PROP;
```

**PTU: genuinely two different receptor affinities for the same drug — kept
as two named interfaces, not force-collapsed into one.** Unlike MMI and
propranolol, PTU's two downstream uses in the original use **two different
EC50s**: `IC50_PTU = 1.0` for TPO inhibition, `IC50_PTU_D1 = 5.0` for
D1-deiodinase inhibition (a real, distinct pharmacological affinity — PTU is
well known clinically for this extra-thyroidal D1-inhibiting activity that
MMI lacks). Collapsing these into one shared `EFFECT_PTU` term would either
invent a false shared EC50 (misrepresenting one of the two real affinities)
or silently drop one mechanism — both ruled out by the guide's "don't force
a bad fit" principle. Instead, per the guide's explicit allowance for a
structure that doesn't fit the single-`EFFECT_<STEM>` mold, PTU exposes
**two** named Hill interfaces from its one concentration:

```
EFFECT_PTU    = EMAX_PTU    * pow(C_PTU, GAMMA_PTU)    / (pow(EC50_PTU,    GAMMA_PTU)    + pow(C_PTU, GAMMA_PTU));    // TPO inhibition
EFFECT_PTU_D1 = EMAX_PTU_D1 * pow(C_PTU, GAMMA_PTU_D1) / (pow(EC50_PTU_D1, GAMMA_PTU_D1) + pow(C_PTU, GAMMA_PTU_D1)); // D1-deiodinase inhibition
...
TPO_inh_PTU = TPO_inh_max * EFFECT_PTU;   // was: TPO_inh_max * PTU_conc/(IC50_PTU+PTU_conc)
D1_inh_PTU  = EFFECT_PTU_D1;              // was: PTU_conc/(IC50_PTU_D1+PTU_conc), scaled by *0.5 downstream, unchanged
```

Both are exact renames (`EMAX=GAMMA=1` for each, reproducing each original
ratio exactly) — no `nls()` fit needed for any of the three drugs, since
every one of the original's effect terms was already a plain Hill ratio.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Tried declaring `C_MMI = 0` etc. directly in `$PARAM` (as the guide's
qspserver-compatibility section literally describes) and assigning the
recomputed value in `$ODE`; this does not compile under mrgsolve 2.0.1 for
the same reason already documented in `pagets-disease/pbd_refactor_notes.md`
and repeated in every refactor in this batch (`$PARAM` values are read-only
references inside `$ODE`, so a `double C_MMI = ...;` after `$PARAM C_MMI =
0` is a redeclaration/read-only-assignment error). Followed the same,
actually-compiling precedent: `C_MMI`, `EFFECT_MMI`, `C_PTU`, `EFFECT_PTU`,
`EFFECT_PTU_D1`, `C_PROP`, `EFFECT_PROP` are plain `double` locals computed
in `$ODE`, listed in `$CAPTURE` (replacing the original's now-invalid
compartment-duplicating `$CAPTURE` line — see defect #93 above). Confirmed
via `POST /model_manifest` on the refactored DSL: all seven appear in
`outputPaths` — discoverable through the output side of the manifest, not
the overridable-parameters side, exactly as the pagets precedent describes.

## Verification

**Method.** Extracted the bare DSL text from both `gd_mrgsolve_model.R`
(`graves_model <- '...'`, unmodified) and `gd_mrgsolve_model_refactored.R`
(same variable name, refactored) and ran them through the local qspserver
`mrgsolve_api` (`http://localhost:8007`), `POST /model_manifest` then
`POST /run_simulation`, requests spaced ~2.2s apart per the shared-service
note. Since the untouched original does not compile at all (defect #93), the
comparison baseline is a **patched-original** scratch copy — the original
DSL with only the same two syntax-only fixes applied (`F_MMI` →
`Fbio_MMI_unused`, i.e. the parameter renamed rather than the compartment,
to isolate the build fix from the refactor's own renaming; and the invalid
`$CAPTURE` line dropped) — confirmed compiling cleanly via `/model_manifest`
before use. This scratch copy was never written into the repo tree and is
not part of the deliverables.

**A second API constraint hit and worked around (not a defect, a request-
contract limitation):** the API's `/run_simulation` contract
(`app-mrgsolve-api/app/r_scripts/_runner.R`) has no `init`-override field —
only `parameters` ($PARAM overrides) and `dosing` (event-table bolus/infusion
records) are exposed. The original R scenarios set custom initial hormone
levels via mrgsolve's `init(TRAb_C=25, TSH_C=0.01, fT4_C=35, ...)`, which the
API cannot replicate directly. Confirmed this matters, not just cosmetic: a
first attempt using the model's bare `$MAIN` defaults (which leave
`TSH_C`/`T4_C`/`T3_C`/`fT4_C`/`fT3_C`/`rT3_C` all at 0) produced **all-NaN**
output from `t=0`, because `TSH_prod = ksyn_TSH / fb_tot` divides by zero
when `fT4_C = fT3_C = 0`. Worked around by adding `t=0` bolus dosing records
that ADD the exact delta from each compartment's own known default (from
`$MAIN`, or 0 where `$MAIN` sets nothing) up to the scenario's stated
`init()` value — e.g. `fT4_C` default 0 → bolus `amt=35` reaches the
scenario's `fT4_C=35`; `GO_act` default `TRAb0*kGO_TRAb/kGO_kel=133.3333`
→ bolus `amt=35-133.3333` reaches the scenario's `GO_act=35`. This is
disclosed as a verification-harness technique, not a change to either model.

**Scenarios run — the file's own, not invented, plus one disclosed exception
for PTU.**

1. **Scenario 2, "Methimazole 30mg/day"** (`USE_MMI=1`; `ev(amt=0.10,
   ii=0.33, addl=2189, cmt=1, time=0)`; 90-day window, `delta=1` — shortened
   from the original's 730-day/daily-sampling run per the guide's
   solver-budget allowance, long enough to reach MMI's own steady
   oscillation many times over given its `CL_MMI=4/day` elimination):
   **exact match, max abs diff 0.0**, on every shared deterministic output
   (`RAI_SERUM, RAI_THYR, TSH_C, T4_C, T3_C, fT4_C, fT3_C, rT3_C, TRAb_C,
   Bcell_C, BoneResor, HR_dev, GO_act, ThyMass, TPO_inh, HR_obs, BMD_loss`).
   `MMI` (original) vs `C_MMI` (refactored): exact match, `0.0522` at t=90 in
   both. `EFFECT_MMI = 0.2583 = 0.0522/(0.15+0.0522)` — correct Hill
   relationship.
2. **Scenario 4, "MMI + Propranolol"** (`USE_MMI=1, USE_PROP=1`; both
   `ev()`s, same 90-day window): **exact match, max abs diff 0.0** on every
   shared output. `MMI`/`C_MMI`: `0.0522` in both. `PROP_C`/`C_PROP`:
   `0.0884` in both. `EFFECT_PROP = 0.8156 = 0.0884/(0.02+0.0884)` — correct.
3. **Scenario 6, "Rituximab + MMI (refractory GO)"** (`USE_MMI=1, USE_RTX=1,
   USE_GCS=1`; `rtx_events` — two boluses into `cmt=1`, i.e. `MMI`/`CENT_MMI`,
   per the #94 contamination defect above, preserved unchanged): **exact
   match to floating-point-scale** — `TPO_inh` max abs diff `~1e-15`,
   `Bcell_C` max abs diff `~1e-35`, every other shared output `0.0`. `MMI`
   vs `C_MMI` both decay to `~-2.098e-28` by t=90 (numerically zero,
   floating-point noise around the true zero — two small `amt=1` boluses
   fully cleared by `CL_MMI=4/day` over 90 days), confirming the spurious
   rituximab-into-MMI contamination reproduces identically in both models,
   not merely "doesn't break" — the same numbers, bug included.
4. **PTU — no scenario in the original file doses it, so an additional,
   explicitly non-original scenario was used**, disclosed as such rather
   than silently substituted for a real one: `USE_PTU=1`, `ev(amt=2.0,
   ii=0.25, addl=359, cmt=2, time=0)` (PTU q6h for 90 days, arbitrary dose
   chosen only to produce a non-trivial concentration), same 90-day window.
   **Exact match, max abs diff 0.0**, on `TPO_inh, TSH_C, fT4_C, fT3_C`.
   `PTU_C` (original) vs `C_PTU` (refactored): exact match, `0.2204` at t=90
   in both. `EFFECT_PTU = 0.1806 = 0.2204/(1.0+0.2204)` (TPO-inhibition
   mechanism) and `EFFECT_PTU_D1 = 0.0422 = 0.2204/(5.0+0.2204)`
   (D1-deiodinase mechanism) — both correct, confirming the two-interface
   split reproduces the original's two distinct ratios exactly.

The four `_obs` outputs (`TSH_obs, fT4_obs, fT3_obs, TRAb_obs`) carry
`$SIGMA`/`EPS(1)` proportional residual error and are **not** compared for
exact match — they differ between any two independent stochastic draws of
the *same* model, seed not fixed by either the original or this refactor, so
a nonzero diff there is expected noise, not a mismatch. Their deterministic
inputs (`TSH_C, fT4_C, fT3_C, TRAb_C`) and the fully deterministic `HR_obs`,
`BMD_loss` all match exactly, which is the meaningful comparison.

This is the expected outcome for three Archetype-1 pure structural
reorganizations (renames plus one drug's Hill term split into two, with no
`nls()` fitting anywhere) — confirmed genuinely `0.0`/floating-point-scale
for every deterministic point, not a loosened tolerance.

`/model_manifest` on the refactored DSL additionally confirmed all of
`KA_MMI, CL_MMI, V1_MMI, F_MMI, EC50_MMI, EMAX_MMI, GAMMA_MMI, KA_PTU,
CL_PTU, V1_PTU, F_PTU, EC50_PTU, EMAX_PTU, GAMMA_PTU, EC50_PTU_D1,
EMAX_PTU_D1, GAMMA_PTU_D1, KA_PROP, CL_PROP, V1_PROP, F_PROP, EC50_PROP,
EMAX_PROP, GAMMA_PROP` appear in `parameters` with their original numeric
defaults, and `CENT_MMI, CENT_PTU, CENT_PROP, C_MMI, EFFECT_MMI, C_PTU,
EFFECT_PTU, EFFECT_PTU_D1, C_PROP, EFFECT_PROP` all appear in `outputPaths`.

No scratch/debug artifacts were left in the repository — DSL extraction,
the patched-original scratch copy, and the comparison scripts ran entirely
from a session-local scratchpad directory outside the repo tree, deleted
after use.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the three
genuine-drug rows (MMI, PROP, PTU) marked refactored/verified; the four
endogenous rows (T3, T4, TRAB, TSH) reclassified as "not a drug — endogenous
disease state" rather than force-fitted PK rows, per the task's specific
instruction.
