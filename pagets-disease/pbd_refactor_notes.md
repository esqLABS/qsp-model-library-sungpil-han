# Refactor notes — `pagets-disease/pbd_mrgsolve_model.R`

Scope: **Calcitonin (CTN)**, **Denosumab (DMB)**, and **Zoledronic Acid
(ZA)** — the three PBD-treatment compounds carrying rows in
`driver-patches/data/compound_perturbation_census.md`. **Alendronate
(ALN)** is out of scope and is byte-for-byte identical in behavior in
`pbd_mrgsolve_model_refactored.R`; see the "Shared-compartment edge case"
section below for the one mechanical exception this required.

## Archetypes

### Zoledronic Acid (ZA) — Archetype 2 (2-compartment linear, no depot) + a bespoke irreversible bone-trap sink

The original's `ZA_cen`/`ZA_per` pair is a standard central+peripheral
linear PK system with clearance/volumes already given as `CL_ZA`/`Vc_ZA`/
`Vp_ZA`/`Q_ZA` — Archetype 2, confirmed against the equations (no depot is
exercised by ZA's own dosing; ZA's own scenarios dose straight IV bolus
into the central compartment). Renamed to convention: `Vc_ZA`→`V1_ZA`,
`Vp_ZA`→`V2_ZA`, `ZA_cen`→`CENT_ZA`, `ZA_per`→`PERI_ZA`.

The original additionally tracks `ZA_bon`, an irreversible cumulative
bone-binding sink fed by `k_bon_ZA * ZA_cen` (a one-way trap, not a
reversible peripheral compartment) — this is real, well-documented
bisphosphonate pharmacology (irreversible incorporation into
hydroxyapatite) and is not one of the guide's four named PK roles.
Preserved exactly, renamed to `BON_ZA` (`k_bon_ZA`→`KBON_ZA`), per the
guide's "don't flatten a mechanistically rich PK model... because it's
inconvenient" and "if a compound's PK genuinely doesn't fit any of the
guide's four archetypes, don't force it" — this is Archetype 2 plus one
bespoke addition, not a flattening.

`C_ZA` is a straight pass-through of `CENT_ZA`: the original's `ZA_cen` is
already in concentration units (ng/mL comment in `$CMT`), not an amount
divided by volume — the dose itself is pre-divided by `Vc_ZA` in the R
driver before being injected (`dose_amt <- 5e6 / 18`). There is no
amount-space representation in the original to convert from, so `C_ZA =
CENT_ZA` (no division) is the faithful rename, not an invented shortcut.

### Calcitonin (CTN) — Archetype 3 without a peripheral compartment (depot + central only)

`CTN_abs`→`CENT_CTN`-style depot-to-central absorption with first-order
elimination (`CL_CTN`/`Vc_CTN`/`ka_CTN`), no peripheral compartment in the
original. Renamed to convention: `CTN_abs`→`GUT_CTN`, `CTN_cen`→
`CENT_CTN`, `Vc_CTN`→`V1_CTN`, `ka_CTN`→`KA_CTN`. Same concentration-space
design as ZA (`CTN_cen` is already pg/mL) — `C_CTN = CENT_CTN`, pass-through.

### Denosumab (DMB) — Archetype 3 (depot+central+peripheral, linear) + a bespoke mass-action RANKL-binding sink, NOT full TMDD

Checked explicitly per the task's instruction to look for TMDD/receptor-
binding kinetics, since DMB is an anti-RANKL antibody. The original's own
comment is accurate: `"RANKL neutralization (target-mediated): simplified
as linear clearance + RANKL binding"`. What is actually modeled:

- `DMB_abs`(amount, μg) → `DMB_cen`(conc, μg/mL) ↔ `DMB_per`(conc, μg/mL),
  linear `CL_DMB`/`Q_DMB`/`Vc_DMB`/`Vp_DMB`/`ka_DMB` — a completely
  standard Archetype 3 PK skeleton, renamed to `GUT_DMB`/`CENT_DMB`/
  `PERI_DMB`/`V1_DMB`/`V2_DMB`/`KA_DMB`.
- A single bimolecular loss term, `kbind_DMB * DMB_cen * RANKL_free`,
  appears symmetrically in `dxdt_DMB_cen` (draining DMB) and
  `dxdt_RANKL_free` (draining RANKL) — genuine mass-action target binding,
  not an algebraic occupancy shortcut, so it is not simply dropped.

This is **not** Archetype 4 (full TMDD): there is no complex/bound-drug
compartment, no `RTOT`, and no reverse (dissociation) reaction — `kdiss_DMB`
is declared in `$PARAM` but never referenced anywhere in `$ODE`/`$MAIN` in
the original (confirmed: single occurrence, the declaration itself).
Forcing an Archetype-4 shape (adding a `COMPLEX_DMB` compartment and a
`KOFF_DMB` term) would be adding pharmacology the original does not have,
which the guide explicitly prohibits ("don't add or remove compartments to
make it fit"). Renamed `kbind_DMB`→`KBIND_DMB`, `kdiss_DMB`→`KDISS_DMB`
(kept, disclosed unused). `C_DMB = CENT_DMB`, pass-through (same
concentration-space design as ZA/CTN).

## Hill interface

### ZA and CTN: rename, not a fit

Both `ZA_inhib_OC = Emax_ZA_OC * ZA_cen / (IC50_ZA_OC + ZA_cen)` and
`CTN_inhib_OC = Emax_CTN_OC * CTN_cen / (IC50_CTN_OC + CTN_cen)` are
*already* exactly the guide's Hill shape with `gamma = 1` — no ODE-derived
kinetics to fit, no receptor system to hold at steady state. Renamed
directly:

| Original | Refactored |
|---|---|
| `Emax_ZA_OC` | `EMAX_ZA` |
| `IC50_ZA_OC` | `EC50_ZA` |
| `Emax_CTN_OC` | `EMAX_CTN` |
| `IC50_CTN_OC` | `EC50_CTN` |
| `ZA_inhib_OC` | `EFFECT_ZA` |
| `CTN_inhib_OC` | `EFFECT_CTN` |

`GAMMA_ZA = 1` and `GAMMA_CTN = 1` are new (no explicit Hill coefficient in
the original), added per the guide's explicit allowance for this case.
Both `EMAX_ZA_OC`/`IC50_ZA_OC`/`Emax_CTN_OC`/`IC50_CTN_OC` were originally
declared in the "Osteoclast (OC) Dynamics" `$PARAM` section (not physically
next to the ZA/CTN PK parameters) — relocated to each compound's own
`$PARAM` block alongside `CL_<STEM>` etc., values unchanged, so the
Hill-interface parameters live with the rest of that compound's interface.

**Pre-existing dead code found while doing this (not fixed, disclosed):**
`OC_inhib_total` (a Bliss-independence combination of the two inhibition
fractions) is computed and then immediately overwritten — the very next
statement recomputes `dxdt_OC` from a completely different formula
(`eff_kdeg_OC = kdeg_OC * (1.0 + 3.0*ZA_inhib_OC + 1.5*CTN_inhib_OC)`,
additive with fixed 3x/1.5x weights, not Bliss independence) inside a
bare `{ }` block, and that second assignment is what survives. `OC_inhib_total`
and the first `dxdt_OC` assignment have **zero effect on simulated
dynamics** — confirmed by reading execution order (C++ statements execute
top-to-bottom; the last write to `dxdt_OC` wins). Preserved verbatim
(renamed `ZA_inhib_OC`/`CTN_inhib_OC` → `EFFECT_ZA`/`EFFECT_CTN` inside it,
since those are the same variables now under new names) rather than
deleted, since deleting dead code is not a pure rename. The actually-used
combination point is `eff_kdeg_OC`, which reads `EFFECT_ZA`/`EFFECT_CTN`
directly — this already satisfies "combine multiple drugs' effects only at
the point the disease equations actually use them," unchanged by the
refactor.

### DMB: bespoke, not forced into a Hill shape

The original never computes a bounded Emax-Hill effect for DMB — its
disease impact is entirely the mass-action term above, which is **linear
and unbounded** in `C_DMB` (no saturation, no Emax ceiling anywhere).
Fitting `Emax*C^gamma/(EC50^gamma+C^gamma)` to it would either force a
false ceiling (misrepresenting the original, which the guide's "don't
force a fit" principle rules out) or require re-deriving the whole
RANKL/OPG/OC steady-state system as a function of `C_DMB` alone (the
"ODE-derived kinetics" path) purely to manufacture a shape the original
never had. Instead, exposed faithfully in the original's own shape:

```
double EFFECT_DMB = KBIND_DMB * C_DMB;
double RANKL_bind_DMB = EFFECT_DMB * RANKL_free;
```

`EFFECT_DMB` is a **pseudo-first-order rate constant (1/h), not a 0–1
fraction** — this is a factoring of the exact same product the original
computed twice under two different expressions (`kbind_DMB * DMB_cen`
appeared once via `RANKL_bind_DMB` and once again inline in
`dxdt_DMB_cen`), given a name so it satisfies "the compound's effect on
disease is expressed as one named variable." Nothing about `dxdt_RANKL_free`
or `dxdt_DMB_cen` was changed beyond this renaming/factoring — confirmed
by the exact-match verification below.

## Shared-compartment edge case: `ZA_abs` (Alendronate, out of scope)

`ZA_abs` is a depot compartment whose outflow uses `ka_ALN` (Alendronate's
own absorption rate, not a ZA-owned parameter) and feeds directly into
`ZA_cen`/`CENT_ZA`. The original explicitly frames this as reusing ZA's
distribution compartments for a same-class oral bisphosphonate ("Oral
route uses ZA_abs depot (for alendronate analog scenario)"). Since
Alendronate carries no row in the compound census and is out of this
task's scope, `ZA_abs` was **left untouched by name** — it is not one of
ZA's own PK roles (its governing rate constant belongs to a different
compound entirely), so renaming it to `GUT_ZA` would misrepresent it as
ZA's own depot when the parameter driving it says otherwise.

The one unavoidable, purely mechanical consequence: `ZA_abs`'s single
downstream reference (the compartment it flows into) is now named
`CENT_ZA` instead of `ZA_cen`, so `dxdt_ZA_abs` is untouched but the target
of its inflow in `dxdt_CENT_ZA` reads the new name. The R driver's
Alendronate scenarios (`run_alendronate`, and the Alendronate portion of
`run_sequential`) still dose into `cmt = "ZA_abs"`, unchanged, with the
exact same `ka_ALN`/dose amounts. Only the ZA-IV portion of
`run_sequential`, and the standalone ZA scenarios (`run_za_iv`,
`run_za_supportive`), had their `cmt = "ZA_cen"` argument updated to
`"CENT_ZA"` — this is literally ZA's own dosing target, in scope.
Alendronate's numeric behavior is unchanged; only a compartment name
string it happens to flow through downstream was updated for internal
consistency.

## Build-compat fix (non-numeric, disclosed)

The original's `$CAPTURE` lists fifteen names that are already `$CMT`
compartments (`ZA_cen ZA_bon CTN_cen DMB_cen RANKL_free OPG_free OCpre OC
OBpre OB BMD bsALP NTX CTX_s Pain`). mrgsolve 2.0.1 rejects this outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: ZA_cen,ZA_bon,CTN_cen,DMB_cen,RANKL_free,OPG_free,OCpre,OC,OBpre,OB,BMD,bsALP,NTX,CTX_s,Pain
```

Confirmed on the untouched original alone via `POST /model_manifest`
(no changes involved). Per the guide's settled policy, the fifteen names
were removed from `$CAPTURE` (every compartment still appears in mrgsolve
output automatically regardless of `$CAPTURE` — confirmed by diffing
`/model_manifest` `outputPaths` before/after, still lists all fifteen)
and `C_ZA EFFECT_ZA C_CTN EFFECT_CTN C_DMB EFFECT_DMB` were added so the
new interface variables are discoverable. Applied **directly to the
delivered `pbd_mrgsolve_model_refactored.R`**, never to
`pbd_mrgsolve_model.R`. Logged as
`translations/UPSTREAM_ISSUES.md` **#56**.

## A note on the guide's "$PARAM with = 0" instruction for `C_<STEM>`/`EFFECT_<STEM>`

Tried first, as written: declaring `C_ZA = 0` etc. in `$PARAM` and then
assigning the recomputed value in `$ODE` (`C_ZA = CENT_ZA;`). This fails to
compile under mrgsolve 2.0.1 — `$PARAM` values are injected into
`$ODE`/`$MAIN` as **read-only references**, confirmed via the compiler
error surfaced through the qspserver `mrgsolve_api`:

```
234:6: error: assignment of read-only reference 'C_ZA'
```

A parameter that is only ever read cannot also be the thing a `double
C_ZA = CENT_ZA;`-style computation writes into. This is not a corner case
specific to this file — it is why every other model in this refactor batch
(`rheumatoid-arthritis`, `sepsis`, `thyroid-eye-disease`,
`polymyalgia-rheumatica`, `kidney-transplant-rejection`, and the later
scaled-out batches) declares `C_<STEM>`/`EFFECT_<STEM>` as plain `double`
locals in `$MAIN`/`$ODE` instead — none of the 29 prior deliverables in
this repo actually put a computed `C_<STEM>` in `$PARAM`, for the same
reason. This file follows that same, actually-compiling precedent:
`C_ZA`/`C_CTN`/`C_DMB`/`EFFECT_ZA`/`EFFECT_CTN`/`EFFECT_DMB` are `double`
locals, exposed via `$CAPTURE` — `/model_manifest`'s `outputPaths` field
lists them (confirmed below), so they remain fully discoverable, just via
the output side of the manifest rather than the overridable-`parameters`
side. This is worth a maintainer's attention if the guide's wording is
meant to be taken literally in a future pass.

## Pre-existing unused parameters (found, not fixed)

`F_ZA`, `F_CTN`, `F_DMB` (bioavailability) and `KDISS_DMB` (was
`kdiss_DMB`) are declared in `$PARAM` but never referenced anywhere in the
original's `$ODE`/`$MAIN` (each confirmed by exact-word grep: the
declaration line is the only occurrence). Dose scaling for CTN and DMB is
instead done entirely in the R driver's dosing amounts (`dose_ctn_pg <-
800`, `dose_dmb_ug <- 60000`), bypassing these parameters. Preserved as
declared-but-unused, per "don't invent, don't drop the original's
parameters" — not logged as an upstream defect since it doesn't affect
compilation or produce incorrect behavior, just dead weight.

## Verification

Per the guide's mandatory protocol: ran three of the original file's own
dosing scenarios through both the (`$CAPTURE`-patched-only, for
buildability) original and `pbd_mrgsolve_model_refactored.R`, via the
qspserver `mrgsolve_api` (`POST /model_manifest`, `POST /run_simulation`),
requests spaced ~2s apart per the shared-service note. `/model_manifest`
confirmed compilation of both, and confirmed `C_ZA`/`EFFECT_ZA`/`C_CTN`/
`EFFECT_CTN`/`C_DMB`/`EFFECT_DMB` are present in `outputPaths`.

Scenario windows were shortened from the original's own 730-day (2-year)
horizon per the guide's allowance (a 20-compartment ODE system run daily
for 2 years is unnecessarily heavy for a shared, concurrency-limited
service) — each window still fully exercises the relevant compound's
absorption/distribution/clearance and, for DMB, both scheduled doses:

1. **ZA IV** — `run_za_iv`'s own dose (5 mg → 277,777.78 ng/mL-equivalent
   bolus into `ZA_cen`/`CENT_ZA`, `cmt=2`), 180-day window, `delta=12`.
2. **CTN SC** — `run_calcitonin`'s own dosing (800 pg/mL-equivalent daily
   ×180 into `CTN_abs`/`GUT_CTN`, `cmt=5`, via `ii=24, addl=179`), 180-day
   window, `delta=12`.
3. **DMB SC** — `run_denosumab`'s own dosing (60,000 μg into
   `DMB_abs`/`GUT_DMB`, `cmt=7`, at day 0 **and** day 182), 184-day window
   (extended 2 days past the second dose) to exercise both doses,
   `delta=12`.

Every shared output (`ZA_abs, ZA_cen/CENT_ZA, ZA_per/PERI_ZA,
ZA_bon/BON_ZA, CTN_abs/GUT_CTN, CTN_cen/CENT_CTN, DMB_abs/GUT_DMB,
DMB_cen/CENT_DMB, DMB_per/PERI_DMB, RANKL_free, OPG_free, OCpre, OC,
OBpre, OB, BMD, bsALP, NTX, CTX_s, Pain, RANKL_OPG_R, OC_fold, OB_fold,
pct_bsALP, pct_NTX`) was compared point-by-point across the full time grid
for all three scenarios (compartment numbering unchanged, since renaming
was done in place with no reordering/insertion/deletion).

**Result: exact match, max abs diff = 0.0 for every shared output, every
scenario.** `RANKL_OPG_R` (`RANKL_free/OPG_free`) is `NaN` at the earliest
timepoints in *both* the original and refactored runs — the API has no
`init` override (only `parameters`, which map to `$PARAM` not initial
compartment values, and this model has no `$MAIN` init block), so all
compartments start at 0 and `0/0` is `NaN` identically in both; this is a
property of the comparison setup (zero initial conditions for a fair
apples-to-apples run), not a discrepancy — confirmed identical `NaN`
pattern in both series. This is the expected outcome for a pure rename
(ZA, CTN) and for a rename/factoring with zero added approximation (DMB) —
no curve-fit was performed for any of the three compounds, and the
verification bears that out.

Real, non-degenerate dynamics were confirmed alongside the diff (not just
"both are flat/zero"), e.g. refactored `CENT_ZA` = `[0, 277777.78,
62072.39, 36634.09, 22790.84, ...]`, `EFFECT_ZA` reaching its `EMAX_ZA =
0.9` ceiling almost immediately after the bolus, and `C_DMB`/`EFFECT_DMB`
rising through the DMB depot-absorption phase — all bit-identical to the
original's equivalent (renamed) series.

**A one-row reporting artifact, disclosed (does not affect the match
above):** `C_<STEM>`/`EFFECT_<STEM>`, computed as `$ODE`-local doubles,
read as the *pre-dose* value in the observation row that coincides with a
bolus event, catching up to the correct value one row later (e.g.
`CENT_ZA` at `t=12h` correctly shows `277777.78` post-bolus, but `C_ZA` at
the same row reads `0`, then `C_ZA` at `t=24h` correctly shows
`62072.39`, matching `CENT_ZA` exactly from that point on). This is the
same class of dosing/reporting timing artifact independently documented in
`rheumatoid-arthritis/ra_refactor_notes.md` (there for a `$MAIN`-local
double) — it reproduces here for an `$ODE`-local double too, so it is a
general mrgsolve dosing+capture interaction, not something introduced by
where these variables are declared. It affects only the newly-exposed
`C_<STEM>`/`EFFECT_<STEM>` diagnostics (which have no counterpart in the
original to mismatch against) — every raw compartment and `$TABLE`-derived
value compared above shows no such lag and matches exactly at every row.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
three `pagets-disease` rows (Calcitonin, Denosumab, Zoledronic Acid).
