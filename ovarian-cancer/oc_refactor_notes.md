# Refactor notes — `ovarian-cancer/oc_mrgsolve_model.R`

Five compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **Bevacizumab (BEV)**,
**Carboplatin (CAR)**, **Niraparib (NIRA)**, **Olaparib (OLA)**, and
**Paclitaxel (PAC)**. These are the only compounds modeled in this file —
nothing else was touched.

## Archetypes

- **CAR (carboplatin)** — **Archetype 2** (no depot, two-compartment,
  linear elimination), confirmed against the actual equations:
  `dxdt_CAR_C1 = -(CL/V1)*C1 - (Q/V1)*C1 + (Q/V2)*C2`, `dxdt_CAR_C2 =
  (Q/V1)*C1 - (Q/V2)*C2`, IV bolus dosing directly into the central
  compartment.
- **PAC (paclitaxel)** — **Archetype 2, extended with a second peripheral
  compartment** (no depot, three-compartment, linear elimination). The
  file's own header comment claims "Michaelis-Menten nonlinear" PK
  (Gianni 1995 JCO), but the actual `dxdt_PAC_C1` equation is pure linear
  elimination — constant `CL_PAC`, no `Km`/`Vmax` term anywhere. Logged as
  a documentation/implementation mismatch in
  [`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md), item 6.
  The linear three-compartment structure that *is* actually implemented
  was refactored faithfully (`CENT_PAC`/`PERI_PAC`/`PERI2_PAC`).
- **OLA (olaparib)** — **Archetype 3** (depot + central + peripheral,
  oral, linear): `dxdt_OLA_gut = -ka*gut`, `dxdt_OLA_C1 = F*ka*gut/V1 -
  (CL+Q)/V1*C1 + Q/V2*C2`, matching the guide's template almost exactly.
- **NIRA (niraparib)** — **bespoke**: the original declares `ka_NIRA`/
  `F_NIRA` (an absorption rate and bioavailability, mirroring olaparib's
  own depot pattern) but never declares a `NIRA_gut` depot compartment in
  `$CMT`, and `dose_niraparib()` doses directly into `NIRA_C1`. The term
  meant to represent absorption instead reads `NIRA_C1` itself:
  `dxdt_NIRA_C1 = (F_NIRA*ka_NIRA*NIRA_C1)/V1_NIRA - (CL_NIRA/V1_NIRA)*
  NIRA_C1 - (Q_NIRA/V1_NIRA)*NIRA_C1 + (Q_NIRA/V2_NIRA)*NIRA_C2` — a
  stray self-referencing term, not a genuine depot flux. This doesn't
  match any of the guide's four archetypes cleanly (it has the *shape* of
  Archetype 2 but with an erroneous extra term baked into the central
  compartment's own dynamics) — handled as bespoke: renamed to the
  convention and preserved *exactly*, stray term included, rather than
  "fixed" (see "A found defect, not fixed" below and
  [`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md), item 5).
- **BEV (bevacizumab)** — **Archetype 2** (no depot, two-compartment,
  linear, day-scale), confirmed against the actual equations —
  **no TMDD anywhere in this file**: plain linear PK, no receptor/complex
  compartment. VEGF neutralization is a first-order mass-action sink on a
  separate `VEGF` biomarker compartment (`kbind_BEV*C_BEV*VEGF`, renamed
  `VEGF_bind_BEV`, was `BEV_effect`), matching the *pattern* used for the
  same drug in `cervical-cancer/cc_mrgsolve_model.R` — but this file's own
  disease equations (`dxdt_TV`, `dxdt_CA125`, `dxdt_CD8T`) never read
  `VEGF` or any bevacizumab-derived quantity at all, unlike cc's file,
  which does have a real (if separate) anti-angiogenic kill term. See
  "No `EFFECT_BEV`" below.

## Renaming (values unchanged from the original)

| Original | Refactored |
|---|---|
| `CAR_C1` / `CAR_C2` | `CENT_CAR` / `PERI_CAR` |
| `PAC_C1` / `PAC_C2` / `PAC_C3` | `CENT_PAC` / `PERI_PAC` / `PERI2_PAC` |
| `OLA_gut` / `OLA_C1` / `OLA_C2` | `GUT_OLA` / `CENT_OLA` / `PERI_OLA` |
| `NIRA_C1` / `NIRA_C2` | `CENT_NIRA` / `PERI_NIRA` |
| `BEV_C1` / `BEV_C2` | `CENT_BEV` / `PERI_BEV` |
| `ka_OLA` / `ka_NIRA` | `KA_OLA` / `KA_NIRA` (case only, matching `KA_<STEM>`) |
| `BEV_effect` | `VEGF_bind_BEV` (rename only, identical formula) |
| `pac_eff` (`PAC_C1/(PAC_C1+100.0)`) | `EFFECT_PAC` (Hill form made explicit) |
| inline `if(OLA_C1>0.1) 0.8*OLA_C1/(OLA_C1+500)` | `EFFECT_OLA` (Hill form made explicit, guard preserved) |
| inline `if(NIRA_C1>0.1) 0.7*NIRA_C1/(NIRA_C1+2000)` | `EFFECT_NIRA` (Hill form made explicit, guard preserved) |
| `CA125_0` (dead param) | `CA125_BASE` (build-compat rename, see below) |
| `CD8T_0` (dead param) | `CD8T_BASE` (build-compat rename, see below) |

New (not in the original, for the named Hill interface; `EMAX`/`GAMMA`
make explicit a shape the original already had implicitly as a plain
ratio — **rename, not a refit**): `EMAX_PAC=1`, `EC50_PAC=100.0`,
`GAMMA_PAC=1`; `EMAX_OLA=0.8`, `EC50_OLA=500`, `GAMMA_OLA=1`;
`EMAX_NIRA=0.7`, `EC50_NIRA=2000`, `GAMMA_NIRA=1`.

`CL_CAR`/`V1_CAR`/`Q_CAR`/`V2_CAR`, `CL_PAC`/`V1_PAC`/`Q2_PAC`/`V2_PAC`/
`Q3_PAC`/`V3_PAC`, `F_OLA`/`CL_OLA`/`V1_OLA`/`Q_OLA`/`V2_OLA`,
`F_NIRA`/`CL_NIRA`/`V1_NIRA`/`Q_NIRA`/`V2_NIRA`, and
`CL_BEV`/`V1_BEV`/`Q_BEV`/`V2_BEV` already matched the naming convention
and are unchanged.

## Exposed concentration: `C_<STEM>` is a direct alias, not `CENT/V1`

Same finding as `cervical-cancer/cc_refactor_notes.md`: this original's
own dosing functions already divide the bolus dose by `V1` before adding
it to the central compartment (`dose_carboplatin`: `amt=dose_mg/V1`;
`dose_paclitaxel`: `amt=dose_mg*1000/V1`; `dose_bevacizumab`:
`amt=dose_mg/V1`), and the original's own `$TABLE` captured the central
compartment directly as "the concentration" (`capture CAR_Conc = CAR_C1;`
etc., no further `/V1`). Olaparib and niraparib are already
concentration-scaled by construction inside their own `dxdt_*_C1`
equations (the `/V1` division is folded into the absorption flux itself,
so the compartment's own units are already a concentration). So for all
five compounds here, `C_<STEM> = CENT_<STEM>` is a direct alias —
confirmed correct by the verification below (exact match against the
original's own `*_Conc` outputs at every timepoint).

## A found defect, not fixed: niraparib's missing depot

Per the shared "log what you find, don't fix it" rule, niraparib's
missing `NIRA_gut` compartment (see "Archetypes" above) is preserved
*exactly* in the refactored file — `dxdt_CENT_NIRA` still includes the
stray `(F_NIRA*KA_NIRA*CENT_NIRA)/V1_NIRA` term reading the central
compartment's own amount. Numerically this term is small
(`+0.73*0.36/537 ≈ +0.00049` per time unit against a combined
`-(CL_NIRA+Q_NIRA)/V1_NIRA ≈ -0.0392` elimination, about a 1.2% relative
offset) — not catastrophic, but a genuine authoring defect, not something
this refactor invents or corrects. Logged as
[`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md), item 5.

## No `EFFECT_BEV`

Unlike `cervical-cancer/cc_mrgsolve_model.R` (where an inline
`BEV_C1/(BEV_C1+10.0)` ratio inside `kill_bev` gives a real, if separate,
anti-angiogenic tumor-kill term to rename), **this file has no
bevacizumab-driven tumor-kill term at all** — `dxdt_TV`, `dxdt_CA125`,
and `dxdt_CD8T` never read `VEGF` or any bevacizumab-derived quantity.
`VEGF_bind_BEV` (the renamed `BEV_effect`) only feeds `dxdt_VEGF`, which
itself has no downstream consumer. So there is nothing to rename into an
`EFFECT_BEV` — inventing one would add pharmacology the original never
had, out of scope for a rename-only refactor. `C_BEV` is still fully
exposed (renamed, captured, discoverable via `/model_manifest`) for
future external driveability; `VEGF_bind_BEV` is exposed the same way.
Logged as [`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md),
item 8.

## Hill interface: renames, not fits

All disease-facing effect terms in this file were already plain
concentration ratios — no ODE-solved receptor-occupancy kinetics to
approximate, so every one of these is a rename, confirmed exact by the
verification below (no `nls()` fitting was needed or attempted for any of
the five compounds):

- **PAC**: `pac_eff = PAC_C1/(PAC_C1+100.0)` is `1*C^1/(EC50^1+C^1)` with
  implicit `Emax=1`, `gamma=1`. `EMAX_PAC=1`/`GAMMA_PAC=1` were added to
  make this explicit; `EFFECT_PAC` computes the identical ratio on
  `C_PAC`.
- **OLA**: `if(OLA_C1>0.1) parp_trap += 0.8*OLA_C1/(OLA_C1+500)` is
  `0.8*C^1/(EC50^1+C^1)` guarded by a threshold. `EMAX_OLA=0.8`/
  `GAMMA_OLA=1` were added; `EFFECT_OLA` computes the identical guarded
  ratio on `C_OLA` (the `>0.1` guard is preserved exactly — it is not
  merely a divide-by-zero safeguard, since the ratio is well-defined at
  `C_OLA=0` too; skipping it below the threshold changes the value by a
  few parts in 10⁴ in that narrow band, so it was kept, not "cleaned up").
- **NIRA**: same shape as OLA, `EMAX_NIRA=0.7`, `EC50_NIRA=2000`, guard
  preserved identically.
- **CAR**: no Hill/Emax term at all for carboplatin's own kill mechanism
  — `kill_Pt` is driven by `Pt_DNA` (a linear adduct-turnover state, not
  a saturating ratio), which is unchanged (renamed inputs only:
  `dxdt_Pt_DNA = k_adduct*C_CAR - k_repair*(1-0.6*eff_HRD)*Pt_DNA`).
- **BEV**: no Hill/Emax term exists to rename — see "No `EFFECT_BEV`"
  above.

Per the guide's "multiple drugs, one pathway" rule, `EFFECT_OLA` and
`EFFECT_NIRA` are kept fully separate through `$ODE`/`$TABLE` and combined
only at the point `dxdt_HRD` actually uses them: `dxdt_HRD = k_HRD_in *
(eff_HRD * (EFFECT_OLA + EFFECT_NIRA)) - k_HRD_out * HRD` — identical
arithmetic to the original's `parp_trap = (guarded OLA term) + (guarded
NIRA term); parp_trap *= eff_HRD; dxdt_HRD = k_HRD_in*parp_trap - ...`.

## `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>` (not `$PARAM`)

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` for direct `/model_manifest` discoverability.
Not done here, for the same reason already documented in the AMD,
membranous-nephropathy, breast-cancer, and cervical-cancer refactors:
mrgsolve 2.0.1 compiles `$PARAM` members as **read-only references**
inside `$ODE`, so a value recomputed every timestep from state cannot
also be declared in `$PARAM`. All nine (`C_CAR`, `C_PAC`, `C_OLA`,
`C_NIRA`, `C_BEV`, `EFFECT_PAC`, `EFFECT_OLA`, `EFFECT_NIRA`,
`VEGF_bind_BEV`) are predeclared as `double`s in `$GLOBAL` and listed in
`$CAPTURE`: visible in every simulation's output columns and in
`/model_manifest`'s `outputPaths`, just not in its `parameters` list. As
in the cervical-cancer refactor, each is *also* recomputed directly from
the live `$CMT` state inside `$TABLE` (redundantly with the identical
computation already needed inside `$ODE`), to avoid the same
simultaneous-dose/observation-row `$GLOBAL` staleness bug that refactor
found and fixed — confirmed closed here too by the exact-match
verification below (including at every dose/observation-coincident row).

## Pre-existing upstream build defects (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for four
layered, build-blocking reasons plus three non-blocking modeling
findings, all unrelated to the compound PK math itself except item 5
below (niraparib, kept unfixed per the log-don't-fix rule). Logged in
full as [`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md).
Summary:

1. **`$CMT`+`$INIT` jointly redeclare all 18 compartments** (`Duplicated
   model names`). Fixed: the `$INIT` block is dropped; its 18 assignments
   become `<CMT>_0 = value;` lines in `$MAIN` (identical values, renamed
   compartment identifiers).
2. **Two `$PARAM` names collide with mrgsolve's auto-generated `<CMT>_0`
   init symbols**: `CA125_0`, `CD8T_0` (both dead — never read anywhere
   except two already-broken `$MAIN` lines with a typo'd target,
   `CA125_0_`/`CD8T_0_`, that never reached mrgsolve's real init symbols
   even before the collision). Renamed to `CA125_BASE`/`CD8T_BASE`.
3. **A `$ODE`-local `double VEGF_free` collides with `$TABLE`'s `capture
   VEGF_free`** (`redefinition of {anonymous}::VEGF_free`). Fixed: the
   local alias is removed; its one use (feeding `VEGF_bind_BEV`) reads
   `VEGF` directly instead — an identity substitution, not a numeric
   change.
4. **Six `$ODE` lines clamp a compartment state directly**
   (`if(Pt_DNA<0) Pt_DNA=0;` and five more), which mrgsolve 2.0.1 rejects
   as `assignment of read-only reference`. Confirmed (same reasoning as
   the cervical-cancer refactor's #76 finding) that this was **always a
   behavioral no-op** — only `dxdt_*` feeds the integrator, so reassigning
   the bare state inside `$ODE` can never affect the next solver step or
   the reported trajectory. All six were deleted outright; confirmed
   dead, not a numeric change, by the exact-match verification below.

All four fixes are syntax-only and non-numeric, applied **directly to
the delivered `oc_mrgsolve_model_refactored.R`**, per the guide's settled
policy for a non-compiling original. The checked-in `oc_mrgsolve_model.R`
is left completely untouched and still carries all four defects exactly
as written. The same four fixes were applied identically to an
in-memory-only scratch copy of the original purely so it could build for
the comparison below.

5-7. See "Non-blocking modeling defects" 5-7 in
[`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md): niraparib's
missing depot (kept unfixed, see "A found defect, not fixed" above),
paclitaxel's "Michaelis-Menten" comment vs. its actual linear
implementation, and the olaparib/niraparib dosing-schedule time-unit
mismatch described next.

## The olaparib/niraparib maintenance-dosing horizon defect

Confirmed while building the verification requests: `dose_olaparib()`/
`dose_niraparib()` compute event times as `start_d*24`/`(start_d+dur_d)
*24` — an hour-scaled axis — while the rest of the model (chemo cycle
intervals, Gompertz growth rate, CA-125 turnover, the R script's own
`mrgsim(..., end=730, delta=1)` call labeled "2-year simulation (days)")
runs on a day-scaled axis, with no reconciliation between the two.
Scenario S4's `dose_olaparib(start_d=126, dur_d=604)` schedules its first
dose at `time=3024` and its last at `time=17508`. Run through the
qspserver `mrgsolve_api`, requesting `time: {end: 730, delta: 1}`
produces a returned time grid unioned with the event times that actually
extends to `t=17508` — about 24x the requested/intended horizon — meaning
no olaparib or niraparib dose ever lands inside the model's own intended
0–730-day window at all. This reproduces **identically** (bit-for-bit
matching outputs) on both the original (patched only for the four build
defects above) and the refactored file, since neither the DSL math nor
the R-side `dose_*()` argument *values* were changed — only `cmt=`
identifiers were renamed. Disclosed here and logged as
[`UPSTREAM_ISSUES.md` #79](../translations/UPSTREAM_ISSUES.md), item 7 —
not fixed, since correcting the schedule would be a real behavioral
change to the original's own scenario construction, out of scope for a
rename-only refactor.

Because this defect makes the literal S4/S5/S6 scenarios extend the
solve/output horizon far past 730 (up to `t≈17508`), the verification
below reports both the full extended-horizon comparison (as the API
actually returns it, given the original's own dosing) and confirms this
is where olaparib's and niraparib's own PK blocks *are* actually exercised
(non-zero `OLA_Conc`/`NIRA_Conc`) — despite landing outside the window
the file's own header and summary code assume the run stays within.

## Verification

**Method.** Both the (four-defects-patched, in-memory-only) original and
`oc_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
at `http://localhost:8007` (`POST /model_manifest`, `POST
/run_simulation`), spaced ~2s apart per request. `/model_manifest`
confirmed the refactored file compiles cleanly and lists every
`C_<STEM>`/`EFFECT_<STEM>`/`VEGF_bind_BEV` in `outputPaths` (61
parameters total, including the five new `EMAX`/`EC50`/`GAMMA` sets).

All six of the original file's own dosing scenarios were run, exactly as
defined in its own `dose_*()`/`ev_seq()` calls (dose amounts, `V1`-scaled
`amt`, `time` sequences all copied verbatim, translated to the API's
dosing-record convention — compartment indices are unchanged 1-based
positions, since the refactor renamed identifiers but did not add, drop,
or reorder any compartment; dosing records were time-sorted before
submission, a request-format requirement of the API's dosing path):

1. **S1** — Untreated (natural history)
2. **S2** — Carboplatin + Paclitaxel × 6 cycles (standard 1st line)
3. **S3** — Carbo+Pacli × 6 → Bevacizumab maintenance
4. **S4** — Carbo+Pacli × 6 → Olaparib maintenance (BRCA+, SOLO-1)
5. **S5** — Carbo+Pacli × 6 → Niraparib maintenance (all-comers, PRIMA)
6. **S6** — Carbo+Pacli+Bev × 6 → Olaparib+Bev maintenance (PAOLA-1)

All twelve shared `$CAPTURE`d outputs (`CAR_Conc`, `PAC_Conc`, `OLA_Conc`,
`NIRA_Conc`, `BEV_Conc`, `VEGF_free`, `TumorVol`, `CA125_lvl`,
`PtDNA_rel`, `HRD_dmg`, `CD8T_rel`, `TV_change`) were compared
point-by-point across the full time grid for all six scenarios.

**Result: exact match at floating-point scale, for every output, every
scenario.**

| Scenario | Grid points | Max abs diff (any output) |
|---|---|---|
| S1 (untreated) | 731 | 0.0 |
| S2 (Carbo+Pacli) | 743 | 0.0024 (CA125_lvl, order ~300 → relative ~8e-6) |
| S3 (+Bev maint.) | 765 | 0.0001 (CA125_lvl) |
| S4 (+Olaparib maint.) | 1951 (extended to t≈17508, see defect above) | 0.0024 (CA125_lvl) |
| S5 (+Niraparib maint.) | 1347 (extended to t≈17496) | 0.0024 (CA125_lvl) |
| S6 (PAOLA-1) | 1973 (extended to t≈17508) | 0.0001 (CA125_lvl) |

All differences are floating-point/solver-tolerance scale, consistent
with the original's and the refactor's mathematically identical
equations being evaluated in a different (but equivalent) additive order
(e.g. `-(CL+Q)/V1*C` vs. `-(CL/V1)*C - (Q/V1)*C`) through an adaptive
stiff solver over a long horizon — not evidence of a bug. This is the
expected outcome per the guide's tolerance rule for pure structural
reorganization with no Hill-fitting.

`OLA_Conc` and `NIRA_Conc` were confirmed non-zero and matching exactly
in S4/S5 (e.g. niraparib: `300` at first dose `t=3024`, rising to a
plateau of `594.2019` by `t≈17450`, identical on both sides) — i.e. the
Archetype 3 (olaparib) and bespoke (niraparib) PK blocks were genuinely
exercised by these scenarios, just far outside the model's own intended
0–730-day window (see the dosing-horizon defect above). No solver
step-count issues or timeouts were encountered at any scenario's full
horizon, including the extended ones.

## Anything else worth flagging

- The R-side scenario/plotting code was updated only where it referenced
  a renamed identifier: every `dose_*()` function's `cmt=` argument
  (`"CAR_C1"` → `"CENT_CAR"`, `"PAC_C1"` → `"CENT_PAC"`, `"OLA_gut"` →
  `"GUT_OLA"`, `"NIRA_C1"` → `"CENT_NIRA"`, `"BEV_C1"` → `"CENT_BEV"`).
  All dosing amounts, intervals, `ev_seq()` calls, six scenario
  definitions, and the `$TABLE` capture *names* (`CAR_Conc`, `PAC_Conc`,
  `OLA_Conc`, `NIRA_Conc`, `BEV_Conc`, `VEGF_free`, `TumorVol`,
  `CA125_lvl`, `PtDNA_rel`, `HRD_dmg`, `CD8T_rel`, `TV_change`) are
  otherwise unchanged from the original, so none of the downstream
  R-side summary/plotting code needed any further edits.
- `BRCAmut` is declared in `$PARAM` and set from the R side
  (`mod_S5`/`mod_S6` use `BRCAmut=0`) but is never read anywhere in
  `$MAIN`/`$ODE`/`$TABLE` — a pre-existing unused parameter, unrelated to
  any of the five refactored compounds, left exactly as-is.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`ovarian-cancer | Bevacizumab`, `ovarian-cancer | Carboplatin`,
`ovarian-cancer | Niraparib`, `ovarian-cancer | Olaparib`, and
`ovarian-cancer | Paclitaxel total CL` rows.
