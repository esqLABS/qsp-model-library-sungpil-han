# Refactor notes — `mn_mrgsolve_model.R` (rituximab, cyclophosphamide, tacrolimus)

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only three
compounds' PK blocks and their interfaces into the disease equations were
rewritten: **Rituximab (RTX)**, **Cyclophosphamide (CPX)**, and **Tacrolimus
(TAC)** — the three rows already present for this file in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md).
Avacopan (`AVA_*`, investigational C5aR inhibitor) and the ACEi/ARB RAAS
suppression flag (`ACEI_eff`) are the only other "drugs" in this file and
were left completely untouched — confirmed by diff, see below.

## Archetype determination

### Rituximab (RTX) — target CD20

**Archetype 2 (no depot, two compartments, linear elimination) augmented
with a TMDD-like binding sink**, closest to **Archetype 4** in spirit but
not a clean fit — see the deviations below.

- PK distribution: `CENT_RTX`/`PERI_RTX`, linear CL/Q/V kinetics — a plain
  Archetype 2 structure (`Vc_RTX`/`Vp_RTX` renamed `V1_RTX`/`V2_RTX`,
  `CL_RTX`/`Q_RTX` already matched the convention).
- **Binding sink, not classic TMDD.** The original also models an explicit
  drug–receptor complex (`RTX_CD20_bound`, renamed `COMPLEX_RTX`) with
  `kon_RTX`/`koff_RTX`/`kint_RTX` (renamed `KON_RTX`/`KOFF_RTX`/`KINT_RTX`),
  which is the guide's Archetype 4 shape on paper. Three things keep it from
  being a clean Archetype-4 match:
  1. **The "receptor" in the association flux is `CD20_B`** — the disease's
     own CD20+ B-cell compartment (also driving `Plasma_cells` production
     elsewhere in the model), not a drug-owned pool. Per the precedent set
     in the AMD refactor (VEGF_FREE/VEGF_BOUND kept unrenamed because they
     are shared disease state), **`CD20_B` is NOT renamed to
     `REC_FREE_RTX`**, and **no `RTOT_RTX` is introduced** — `CD20_B` is a
     dynamic compartment with its own turnover (`kprol_B`, `kdeg_B`), not a
     fixed total receptor pool, so there is no well-defined "occupancy
     fraction" to construct.
  2. **The binding sink never feeds back into `CD20_B`.** `ON_RTX`/`OFF_RTX`/
     `INT_RTX` only move mass between `CENT_RTX`, `PERI_RTX`, and
     `COMPLEX_RTX` — `dxdt_CD20_B` does not reference `COMPLEX_RTX` or
     `INT_RTX` at all, confirmed by reading the original's `$ODE` line by
     line. So `COMPLEX_RTX` is a real elimination pathway for the drug (a
     mass sink that depletes `CENT_RTX`), but it has **no pharmacodynamic
     consequence** in this model as written.
  3. **The actual B-cell-depleting effect is a separate, second, unrelated
     term**: `RTX_B_kill = keff_RTX * RTX_conc * CD20_B` (renamed
     `B_KILL_RTX = EFFECT_RTX * CD20_B`), driven directly by **free plasma
     concentration**, not by occupancy or by the complex compartment.

  This is preserved exactly as designed — the guide says "do not flatten a
  mechanistically rich PK model," so `COMPLEX_RTX` stays as a real dynamic
  compartment (it is not decorative to the drug's own PK — it does change
  `CENT_RTX`'s time-course) even though it is decorative to the *disease*
  side. Nothing was added or removed; only renamed.

- **`EFFECT_RTX` is NOT the guide's Emax/EC50/Hill template**, and this is
  the most important disclosure in this file. The Hill interface section
  says: "if the original's effect term is already [an Emax/EC50 ratio], this
  is a rename;" "if it emerges from ODE-solved kinetics, fit an
  approximation." RTX's effect term is neither — it is already a **single,
  closed-form, one-line algebraic function of `C_RTX`** (so nothing needs
  fitting), but that function is `KEFF_RTX * C_RTX` — **bilinear
  (mass-action), unsaturating**. As `C_RTX -> infinity`, the per-capita kill
  rate grows without bound; there is no `EC50_RTX` at which the effect is
  half-maximal, because the original never wrote one. Forcing an
  `EMAX_RTX`/`EC50_RTX`/`GAMMA_RTX` Hill fit onto this term would either (a)
  require picking an arbitrary saturation point the original doesn't have,
  silently changing high-concentration behavior, or (b) require a
  `GAMMA_RTX` fit against no real saturating data (since the original never
  saturates within the scenarios' concentration range). Per the guide's
  "None of these fit" escape hatch and its instruction not to force a bad
  fit, **`EFFECT_RTX = KEFF_RTX * C_RTX` is kept in its native functional
  form** — renamed and exposed exactly like an `EFFECT_<STEM>` term, but
  documented as a rate constant (1/day), not a fraction. `KEFF_RTX` is the
  renamed `keff_RTX`, unchanged in value.

### Cyclophosphamide (CPX) — no receptor target (alkylating cytotoxic)

**Bespoke — none of the guide's 4 PK archetypes fit.** All four archetypes
describe a single compound's own absorption/distribution (depot, central,
peripheral). CPX instead models a **prodrug -> active-metabolite conversion
chain**: `CPx_cent` (parent, renamed `CENT_CPX`) irreversibly converts to
`CPx_metab` (active 4-OH-cyclophosphamide, renamed `MET_CPX`), which is
itself eliminated by a separate first-order rate. This is a metabolic
transformation, not a distributional one, so it does not match Archetype 1
(single compartment), 2 (parent/peripheral), or 3 (depot/central) — none of
those describe a second compound-derived compartment produced *by* the
first. Kept as the original's two sequential compartments, renamed to the
fork's style (`CENT_CPX` = parent, `MET_CPX` = active metabolite; no
existing table entry covers "metabolite," so `MET_` was chosen by
extension of `GUT_`/`PERI_`/`CENT_`).

- `CL_CPx`/`Vd_CPx`/`km_CPx`/`kel_met` renamed `CL_CPX`/`V1_CPX`/`KM_CPX`/
  `KEL_CPX`. Values unchanged.
- **`C_CPX` exposes the metabolite concentration (`MET_CPX`), not the parent
  (`CENT_CPX`)** — because that is what the PD term (`EFFECT_CPX`, plasma
  cell kill) actually reads in the original (`CPx_met_conc`, never
  `CPx_cent`/its concentration). This matches the guide's "exactly one
  concentration variable that PD equations read" — the parent's own
  concentration is PK-internal only, never read by any PD equation, so it
  is not separately exposed.
- `EFFECT_CPX` **is** a rename-not-refit: the original's `Emax_CPx = CPx_Emax
  * CPx_met_conc / (CPx_EC50 + CPx_met_conc + 1e-9)` is already a plain
  Emax/EC50 ratio. `GAMMA_CPX = 1.0` is new (no explicit Hill coefficient in
  the original), reproducing the original ratio exactly.

### Tacrolimus (TAC) — calcineurin inhibitor, no receptor-binding kinetics

**Archetype 3 variant (depot + central, no peripheral)** — the guide
explicitly names this as "drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a
2-compartment-no-depot variant," and this file is the mirror case: depot +
central *without* a peripheral compartment.

- `TAC_gut`/`TAC_blood` renamed `GUT_TAC`/`CENT_TAC`; `Vd_TAC`/`ka_TAC`
  renamed `V1_TAC`/`KA_TAC`. Values unchanged.
- **`F_TAC` (bioavailability, 0.25) is declared in the original's `$PARAM`
  but never applied anywhere in `$ODE`** — `dxdt_TAC_blood` has no `F_TAC`
  term, so the original's own dosing is effectively full-bioavailability
  despite declaring `F_TAC = 0.25`. This looks like an unfinished parameter,
  not a build defect (it doesn't block compilation), so it is **not** logged
  in `UPSTREAM_ISSUES.md`. It is disclosed here instead: `F_TAC` is kept in
  `$PARAM`, unused, exactly as the original left it, so as not to introduce
  a numeric change under the guise of a rename.
- `EFFECT_TAC` **is** a rename-not-refit: the original's `Emax_TAC = TAC_Emax
  * TAC_conc / (TAC_EC50 + TAC_conc + 1e-9)` is already a plain Emax/EC50
  ratio. `GAMMA_TAC = 1.0` is new, reproducing the original exactly.

## Full naming table applied

| Original | Refactored | Compound |
|---|---|---|
| `Vc_RTX` | `V1_RTX` | RTX |
| `Vp_RTX` | `V2_RTX` | RTX |
| `Kd_RTX` | `KD_RTX` (unused in `$ODE`, same as original — see below) | RTX |
| `kon_RTX` | `KON_RTX` | RTX |
| `koff_RTX` | `KOFF_RTX` | RTX |
| `kint_RTX` | `KINT_RTX` | RTX |
| `keff_RTX` | `KEFF_RTX` | RTX |
| `RTX_cent` | `CENT_RTX` | RTX |
| `RTX_peri` | `PERI_RTX` | RTX |
| `RTX_CD20_bound` | `COMPLEX_RTX` | RTX |
| `RTX_conc` | `C_RTX` | RTX |
| (new; was inlined) | `EFFECT_RTX` | RTX |
| `Vd_TAC` | `V1_TAC` | TAC |
| `ka_TAC` | `KA_TAC` | TAC |
| `TAC_gut` | `GUT_TAC` | TAC |
| `TAC_blood` | `CENT_TAC` | TAC |
| `TAC_conc` | `C_TAC` | TAC |
| `TAC_EC50` | `EC50_TAC` | TAC |
| `TAC_Emax` | `EMAX_TAC` | TAC |
| (new; none in original) | `GAMMA_TAC = 1.0` | TAC |
| `Emax_TAC` | `EFFECT_TAC` | TAC |
| `CL_CPx` | `CL_CPX` | CPX |
| `Vd_CPx` | `V1_CPX` | CPX |
| `km_CPx` | `KM_CPX` | CPX |
| `kel_met` | `KEL_CPX` | CPX |
| `CPx_cent` | `CENT_CPX` | CPX |
| `CPx_metab` | `MET_CPX` | CPX |
| `CPx_met_conc` | `C_CPX` | CPX |
| `CPx_EC50` | `EC50_CPX` | CPX |
| `CPx_Emax` | `EMAX_CPX` | CPX |
| (new; none in original) | `GAMMA_CPX = 1.0` | CPX |
| `Emax_CPx` | `EFFECT_CPX` | CPX |

`CL_RTX` and `Q_RTX` already matched the convention and were left as-is
(names only — same values). `AVA_*`, `ACEI_eff`, `Bss`/`kprol_B`/`kdeg_B`,
`PC_ss`/`kprol_PC`/`kdeg_PC`, `Ab_ss`/`kdeg_Ab`, and every disease-side
parameter/compartment (`CD20_B`, `Plasma_cells`, `Anti_PLA2R1`,
`IgG_deposit`, `Complement_MAC`, `Podocyte_inj`, `GBM_thick`, `Proteinuria`,
`Serum_alb`, `eGFR`, `AngII`, `Aldosterone`) are untouched, since none of
them belong to the three scope compounds.

The R-side glue code (`make_events()`'s `ev(cmt = ...)` dosing targets, and
the PK-profile `ggplot`/`aes()` column references in Figure 3) was updated
to match the renamed compartments (`"RTX_cent"` -> `"CENT_RTX"`, `"TAC_gut"`
-> `"GUT_TAC"`, `"CPx_cent"` -> `"CENT_CPX"`, and the plotted
`RTX_cent`/`TAC_blood`/`CPx_metab` columns -> `CENT_RTX`/`CENT_TAC`/
`MET_CPX`) so the delivered file still runs standalone in R end-to-end
(compiles, simulates all 8 of the original's own scenarios, and produces
the same 4 figures), not just through the qspserver API.

## `KD_RTX` and `F_TAC`: unused parameters, preserved as unused

Two parameters are declared but never read anywhere in the original's
`$ODE`:

- `Kd_RTX` (comment: "koff = Kd*kon") — `koff_RTX` is set directly
  (`0.0001`), not derived from `Kd_RTX * kon_RTX` at runtime (their product
  is numerically consistent, `0.010 * 0.01 = 0.0001`, but that is a
  documentation coincidence, not a live dependency).
- `F_TAC` (bioavailability) — see the Tacrolimus section above.

Both are kept in `$PARAM`, renamed (`KD_RTX`, `F_TAC` unchanged), and left
just as unused as in the original. Applying them would change numeric
behavior and is out of scope for a rename-only refactor.

## `$PARAM` members cannot be reassigned in `$ODE` under mrgsolve 2.0.1

Per the guide's qspserver compatibility requirement #2, `C_<STEM>` and
`EFFECT_<STEM>` should ideally live in `$PARAM` (with a `= 0` default) so
they are discoverable via `/model_manifest` without reading the DSL text.
This was tested directly against the running `mrgsolve_api` container
before writing the final file:

```
$PARAM CL = 1.0, V1 = 2.0, C_X = 0.0
$CMT CENT
$ODE
C_X = CENT / V1;
dxdt_CENT = -CL*C_X;
```
fails to build with:
```
64:5: error: assignment of read-only reference 'C_X'
```

This confirms, independently, the same constraint already documented in the
AMD refactor (`age-related-macular-degeneration/amd_refactor_notes.md`):
mrgsolve 2.0.1 compiles `$PARAM` members as read-only references inside
`$ODE`, so a value that must be recomputed every timestep (as `C_RTX`,
`C_TAC`, `C_CPX`, `EFFECT_RTX`, `EFFECT_TAC`, `EFFECT_CPX` all must be, since
they depend on state) cannot also be declared in `$PARAM`. Per that
precedent, all six are declared in `$GLOBAL` (as `double`s) and added to
`$CAPTURE` instead — visible in every simulation's output columns and in
`/model_manifest`'s `outputPaths`, but not in its `parameters` list. This is
a known, disclosed limitation shared with the AMD refactor, not something
new introduced here.

## Pre-existing upstream build defect (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, unrelated to
any of the three refactored compounds: `$CAPTURE` lists eighteen names that
are already `$CMT` compartments (`RTX_cent RTX_peri RTX_CD20_bound
TAC_blood CPx_cent CPx_metab CD20_B Plasma_cells Anti_PLA2R1 IgG_deposit
Complement_MAC Podocyte_inj GBM_thick Proteinuria Serum_alb eGFR AngII
Aldosterone`), which mrgsolve 2.0.1 rejects outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: RTX_cent,RTX_peri,...
```

Confirmed via `POST /model_manifest` on the untouched original, no changes
involved. Logged as
[`UPSTREAM_ISSUES.md` #48](../translations/UPSTREAM_ISSUES.md). Per the
guide's settled policy for this situation ("When the original doesn't
compile at all"), the fix is applied **directly to the delivered
`mn_mrgsolve_model_refactored.R`**, not just to a scratch copy: the
eighteen compartment names were removed from `$CAPTURE`. This is
syntax-only and non-numeric — mrgsolve always reports every compartment's
state regardless of `$CAPTURE` (confirmed: `/model_manifest`'s
`outputPaths` for the refactored file still lists all the renamed
compartments, e.g. `CENT_RTX`, `PERI_RTX`, ..., even though none of them
appear in `$CAPTURE` any more). The checked-in original
(`mn_mrgsolve_model.R`) was left untouched and still carries this defect.

## Verification

**Method.** Both the original's own `mn_code` block (with only the
`$CAPTURE` build-compat fix above applied, so it would compile at all) and
`mn_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
(`POST /model_manifest`, `POST /run_simulation`) at `http://localhost:8007`.
Five scenarios were built from the original file's own `make_events()`/
`run_scenario()` R code (not invented doses), all to the original's own
`sim_duration = 730` days at `delta = 1`:

1. `no_treatment` — no dosing at all (sanity check that untouched
   disease-only code paths are unaffected).
2. `RTX_mono_MENTOR` — 1 g IV bolus-equivalent (`1000/3.1*1000/1000 ≈
   322.58` µg/mL) at day 0 and day 180.
3. `TAC_mono` — 1.75 mg BID (`0.05 mg/kg/day / 2`, 70 kg) from day 0,
   `ii=0.5, addl=1080` (the original's `seq(0, 540, by=0.5)`).
4. `Ponticelli` — 5.0 µg/mL-equivalent CPx daily for 30 days, three cycles
   (days 0–29, 60–89, 120–149), matching the original's `cpx_times`.
5. `RTX_TAC_combo` — TAC BID for the first 180 days (`ii=0.5, addl=358`)
   followed by a single 1 g RTX dose at day 180 (STARMEN-like sequential
   regimen).

All 23 shared `$CAPTURE`d/compartment outputs (`CENT_RTX`/`RTX_cent`,
`PERI_RTX`/`RTX_peri`, `COMPLEX_RTX`/`RTX_CD20_bound`, `GUT_TAC`/`TAC_gut`,
`CENT_TAC`/`TAC_blood`, `CENT_CPX`/`CPx_cent`, `MET_CPX`/`CPx_metab`, plus
all disease/clinical-endpoint states — `CD20_B`, `Plasma_cells`,
`Anti_PLA2R1`, `IgG_deposit`, `Complement_MAC`, `Podocyte_inj`, `GBM_thick`,
`Proteinuria`, `Serum_alb`, `eGFR`, `AngII`, `Aldosterone`, `CR_prot`,
`PR_prot`, `B_depl`, `AntiPLA2R_neg`) were compared pointwise across the
full daily time grid.

**Result: exact match. Maximum absolute and relative deviation = 0.0
(bit-identical) across all five scenarios and every compared output,
including every downstream disease/clinical-endpoint state fed by the three
compounds' PD terms.** This is the expected outcome per the guide's
tolerance table for pure structural reorganization (Archetype 2/3 + the
CPX bespoke chain), since every rename in this file — including
`EFFECT_RTX`'s native bilinear form — was a **rename, not a refit**: nothing
here required curve-fitting an ODE-derived quantity, so no non-zero
tolerance is expected or was observed.

`/model_manifest` was also called on the finished `_refactored.R`'s DSL to
confirm every `CL_/V1_/V2_/Q_/KD_/KON_/KOFF_/KINT_/KEFF_/KA_/F_/EC50_/EMAX_/
GAMMA_/KM_/KEL_` parameter this naming convention promises is present
(68 parameters total, all three compounds' full parameter sets confirmed
present by name).

## Diff scope confirmation

The only lines that differ between `mn_mrgsolve_model.R` and
`mn_mrgsolve_model_refactored.R` are: the three compounds' `$PARAM` lines
(plus 2 new `GAMMA_*` lines and inline `[refactor]` comments), the `$CMT`
compartment names for the 3 compounds' 7 compartments, the `$GLOBAL`
declaration line, the `$MAIN` init lines for those 7 compartments, the
`$ODE` lines that compute `C_RTX`/`C_TAC`/`C_CPX`/`EFFECT_RTX`/
`EFFECT_TAC`/`EFFECT_CPX` and the three compounds' `dxdt_` lines (plus the
two downstream lines, `B_KILL_RTX`/`PC_kill`/`Ab_prod`, that reference the
renamed effect terms), the `$CAPTURE` block (compartment names removed per
the build-compat fix; `C_<STEM>`/`EFFECT_<STEM>` added), and the R-side
`make_events()`/`ggplot` references to the renamed compartments described
above. Every disease-only compartment, every RAAS/avacopan/ACEi parameter
and equation, every plotting/summary block not touching the three
compounds' own PK/PD columns, and the entire scenario-running/figure-saving
code at the bottom of the file are otherwise byte-identical — confirmed by
diffing the two files' embedded DSL blocks (223 changed lines out of 432)
and manually reviewing every hunk.

## Discoverability fix

A corpus-wide audit found that `C_RTX`/`C_TAC`/`C_CPX` were each declared as
a bare `$GLOBAL` forward-declare (`double C_RTX, C_TAC, C_CPX;`) with the
actual value assigned via a bare reassignment (`C_RTX = CENT_RTX;`, no
`double`) inside `$ODE` — correct mrgsolve, but not a single contiguous
`double C_<STEM> = <expr>;` statement, so it is invisible to the driver-PK
dashboard's source-text pattern matching (per FORK_WORKFLOW_GUIDE.md Part 2,
"What makes a compound's PK discoverable by downstream tooling").

**Fix applied, all three compounds:**
- Removed `C_RTX, C_TAC, C_CPX` from the `$GLOBAL` line (nothing else on
  that line, so it now reads `double EFFECT_TAC, EFFECT_CPX, Emax_AVA;` /
  `double EFFECT_RTX;` only — the two `EFFECT_*` names are untouched,
  out of scope for this fix).
- Added, in `$TABLE` immediately before `$CAPTURE`:
  ```
  double C_RTX = CENT_RTX;
  double C_TAC = CENT_TAC;
  double C_CPX = MET_CPX;
  ```
  reusing verbatim the same expressions the original `$ODE` block already
  used (`C_RTX = CENT_RTX;`, `C_TAC = CENT_TAC;`, `C_CPX = MET_CPX;` — the
  active 4-OH-cyclophosphamide metabolite, not the parent `CENT_CPX`,
  exactly as before). The `$ODE`-side bare reassignments are unchanged;
  they now update the member mrgsolve auto-declares from this `$TABLE`
  initializer, the same mechanism already used successfully elsewhere in
  this fork (`clostridioides-difficile-infection` precedent).

**Note on `EC50_RTX`:** the discoverability contract also expects a
same-stem `EC50_<STEM>` parameter. RTX has none — by design, not an
oversight: as documented above and in `$ODE`'s own comment, `EFFECT_RTX` is
a bilinear mass-action kill-rate term (`KEFF_RTX * C_RTX`), not an Emax/EC50
Hill ratio, because the original never saturates as `C_RTX` grows. So
`double C_RTX = ` now appears but `EC50_RTX` legitimately does not; `C_TAC`
and `C_CPX` both pair correctly with `EC50_TAC`/`EC50_CPX`.

**Verification.** `Rscript -e 'parse(...)'` succeeds with zero errors. Both
the pre-fix (`git show HEAD:...`) and post-fix DSL were extracted and run
through the qspserver `mrgsolve_api` (`POST /model_manifest`,
`POST /run_simulation`, ~2s apart) with the file's own `RTX_mono_MENTOR`
scenario (RTX 1000/3.1 µg/mL-equivalent bolus at day 0 and day 180, CL/Q/V
default parameters, `end=730, delta=1`). `/model_manifest` confirms
`C_RTX`, `C_TAC`, `C_CPX` are now present in `outputPaths` alongside
`EC50_TAC`/`EC50_CPX` in the parameter list. `/run_simulation` on both DSL
strings returned identical 731-row time grids; every captured output
(`CENT_RTX`, `C_RTX`, `C_TAC`, `C_CPX`, `EFFECT_RTX`, `EFFECT_TAC`,
`EFFECT_CPX`, `Proteinuria`, `CD20_B`) matched with max absolute
difference **0** at every timepoint — no divergence, not even at the
dose-instant rows. Confirms the relocation changes nothing numeric.
