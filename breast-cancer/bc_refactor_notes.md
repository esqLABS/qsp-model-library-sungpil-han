# Refactor notes — `breast-cancer/bc_mrgsolve_model.R`

Four compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **Letrozole (LETRO)**,
**Olaparib (OLAP)**, **Palbociclib (PALBO)**, **Trastuzumab (TRAS)**. These
are the only compounds modeled in this file — nothing else was touched.

## Archetypes

- **LETRO** — depot + central, linear elimination, no peripheral
  compartment in the original (`GUT_LETRO`/`CENT_LETRO`). A variant of
  archetype 3 with the peripheral leg dropped (as the guide's own archetype
  3 already anticipates: "drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a
  2-compartment-no-depot variant" — the symmetric case, keep the depot and
  drop the peripheral leg, applies the same way).
- **OLAP** — same structure as LETRO: depot + central, linear, no
  peripheral compartment (`GUT_OLAP`/`CENT_OLAP`).
- **PALBO** — depot + central; the original also declared a third
  compartment, `PERI_PALBO`, but it is dead code (see below) and was
  dropped, leaving the same depot+central shape as LETRO/OLAP.
- **TRAS** — no depot, two-compartment, linear elimination
  (`CENT_TRAS`/`PERI_TRAS`, IV bolus/infusion dosing directly into
  `CENT_TRAS`). **Archetype 2, confirmed against the actual equations, not
  assumed from trastuzumab's real-world target-mediated clearance**: the
  original's PK block is `dxdt_CENT_TRAS = -k10*CENT_TRAS - k12*CENT_TRAS +
  k21*PERI_TRAS` / `dxdt_PERI_TRAS = k12*CENT_TRAS - k21*PERI_TRAS` — plain
  linear two-compartment disposition, with the day⁻¹-to-hour⁻¹ unit
  conversion carried over unchanged. There is no receptor-occupancy state
  variable anywhere in the file; HER2 blockade is a phenomenological Hill
  ratio on free-drug concentration (`Emax_HER2 * Cp_tras / (EC50_HER2 +
  Cp_tras)`), not target-mediated disposition. This is **not** archetype 4
  — nothing was flattened, because there was no TMDD structure to begin
  with.

## Renaming (values unchanged from the original)

| Original | Refactored |
|---|---|
| `ka_palbo` / `Vd_palbo` / `CL_palbo` | `KA_PALBO` / `V1_PALBO` / `CL_PALBO` |
| `Emax_CDK` / `EC50_CDK` / `hill_CDK` | `EMAX_PALBO` / `EC50_PALBO` / `GAMMA_PALBO` |
| `ka_letro` / `Vd_letro` / `CL_letro` | `KA_LETRO` / `V1_LETRO` / `CL_LETRO` |
| `KAI` | `EC50_LETRO` |
| `CL_tras` / `Vd1_tras` / `Vd2_tras` / `Q_tras` | `CL_TRAS` / `V1_TRAS` / `V2_TRAS` / `Q_TRAS` |
| `Emax_HER2` / `EC50_HER2` | `EMAX_TRAS` / `EC50_TRAS` |
| `ka_olap` / `Vd_olap` / `CL_olap` | `KA_OLAP` / `V1_OLAP` / `CL_OLAP` |
| `Cp_palbo` / `Cp_letro` / `Cp_tras` / `Cp_olap` | `C_PALBO` / `C_LETRO` / `C_TRAS` / `C_OLAP` |
| `CDK_inh` | `EFFECT_PALBO` |
| `HER2_block` | `EFFECT_TRAS` |
| `letro_inh` | `EFFECT_LETRO` |

New (not present in the original, for the named Hill interface): `EMAX_LETRO
= 1`, `GAMMA_LETRO = 1`, `GAMMA_TRAS = 1` — these make explicit a shape the
original already had implicitly as a plain ratio (see "Hill interface"
below). `EMAX_PALBO`/`EC50_PALBO`/`GAMMA_PALBO` are renames, not new values —
the original already had an explicit Hill term for palbociclib.

`GUT_PALBO`, `CENT_PALBO`, `GUT_LETRO`, `CENT_LETRO`, `CENT_TRAS`,
`PERI_TRAS`, `GUT_OLAP`, `CENT_OLAP` already matched the compartment-naming
convention and are unchanged.

## Dropped as dead code (verified to change nothing)

- **`PERI_PALBO`**: declared in the original's `$CMT` block but
  `dxdt_PERI_PALBO = 0` and it is never read anywhere else in the model
  (`Cp_palbo` is computed only from `CENT_PALBO / Vd_palbo`). It stays at
  its initial value of 0 for the entire simulation and contributes nothing.
  Removed in the refactor.
- **`Cp2_tras`**: a local `double` computed in the original
  (`PERI_TRAS / Vd2_tras`) but never subsequently read by anything (not by
  an ODE, not by `$TABLE`, not by `$CAPTURE`). Removed as dead code.

Both removals were confirmed, not assumed: the verification run below
reproduces every other output exactly, which would not be possible if
either dropped quantity had actually fed into anything.

## No `EFFECT_OLAP`

Unlike the other three compounds, **the original never connects `Cp_olap`
to any disease equation at all.** Olaparib gets its own PK block
(`GUT_OLAP`/`CENT_OLAP`, absorption, clearance) and `Cp_olap` is captured as
an output, but no `$ODE` line, `$MAIN` expression, or `$TABLE` term ever
reads it. The scenario that is supposed to show olaparib's pharmacology
(Scenario 6, "Olaparib (OlympiAD, BRCAm)") in fact produces its different
tumor trajectory purely through **static `$PARAM` overrides passed at call
time** (`kprol = 0.0012`, `EC50_CDK = 200.0`, renamed `EC50_PALBO` here) —
i.e. by hand-tuning the CDK4/6 and proliferation parameters for the
scenario, not by any mechanistic action of the olaparib molecule itself.
This looks like an authoring shortcut in the original rather than a
deliberate PARP-inhibition mechanism.

Per the guide ("inventing an effect the original doesn't have" is out of
scope for a rename), no `EFFECT_OLAP` was added. `C_OLAP` is still exposed
in `$GLOBAL`/`$CAPTURE` (satisfying the "exposes exactly one concentration
variable" requirement for future driverability), but there is no
compound-effect term to rename because none exists in the source. This is
flagged here as a finding about the original, not fixed — see the
corresponding note added to
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md).

## Hill interface: renames, not fits

All three compounds that do have a disease-facing effect term were already
written as a plain concentration ratio — no ODE-solved receptor kinetics to
approximate, so every one of these is a rename, confirmed by the exact
verification match below (not merely "close"):

- **PALBO**: `CDK_inh = Emax_CDK * Cp_palbo^hill_CDK / (EC50_CDK^hill_CDK +
  Cp_palbo^hill_CDK)` was already an explicit three-parameter Hill term
  (`gamma = 1.5`). `EFFECT_PALBO` uses the identical expression with the
  renamed parameters.
- **TRAS**: `HER2_block = Emax_HER2 * Cp_tras / (EC50_HER2 + Cp_tras)` is
  algebraically `Emax*C^1/(EC50^1+C^1)` — a Hill term with an implicit
  `gamma = 1`. `GAMMA_TRAS = 1` was added purely to make this explicit;
  `EFFECT_TRAS` computes the exact same ratio.
- **LETRO**: `letro_inh = Cp_letro / (Cp_letro + KAI)` is likewise
  `1*C^1/(EC50^1+C^1)` with implicit `Emax = 1`, `gamma = 1`. `EMAX_LETRO =
  1` and `GAMMA_LETRO = 1` were added to make this explicit; `EFFECT_LETRO`
  computes the exact same ratio and feeds into `dxdt_AROMATASE` with the
  same `* 10.0` scaling factor the original used, untouched.

No `nls()` fitting was needed or performed for any of the four compounds.

## `$PARAM` vs `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>`

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` (with a `= 0` default) for direct
`/model_manifest` discoverability. This was not done here, for the same
reason already documented in the AMD and membranous-nephropathy refactors:
mrgsolve 2.0.1 compiles `$PARAM` members as **read-only references** inside
`$ODE`, so a value that must be recomputed every timestep from state
(`C_PALBO`, `C_LETRO`, `C_TRAS`, `C_OLAP`, `EFFECT_PALBO`, `EFFECT_LETRO`,
`EFFECT_TRAS` — all of which read a compartment) cannot also be declared in
`$PARAM`. Instead, all seven are predeclared as `double`s in `$GLOBAL` and
listed in `$CAPTURE`: they are visible in every simulation's output columns
and in `/model_manifest`'s `outputPaths`, just not in its `parameters`
list. This was independently re-confirmed against the live `mrgsolve_api`
container with the same minimal reproduction used in the prior write-ups
(a `$PARAM ... C_X = 0.0` declaration reassigned in `$ODE` fails to build
with `error: assignment of read-only reference 'C_X'`).

## Pre-existing upstream build defect (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for a reason
unrelated to any of the four refactored compounds' own PK: `$CAPTURE` lists
twelve names that are already `$CMT` compartments (`TUMOR CSC ER_SIGNAL
CDK46_ACT HER2_SIGNAL PD_L1 CD8_EFF TREG E2_PLASMA AROMATASE Ki67 CA153`),
which mrgsolve 2.0.1 rejects outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: TUMOR,CSC,ER_SIGNAL,CDK46_ACT,HER2_SIGNAL,PD_L1,CD8_EFF,TREG,E2_PLASMA,AROMATASE,Ki67,CA153
```

Confirmed via `POST /model_manifest` on the untouched original's own DSL,
no changes involved. Logged as
[`UPSTREAM_ISSUES.md` #57](../translations/UPSTREAM_ISSUES.md). Per the
guide's settled policy for this situation, the fix is applied **directly to
the delivered `bc_mrgsolve_model_refactored.R`**, not just to a scratch
copy: the twelve compartment names were removed from `$CAPTURE`. This is
syntax-only and non-numeric — mrgsolve always reports every compartment's
state regardless of `$CAPTURE` (confirmed: `/model_manifest`'s
`outputPaths` for the refactored file still lists all twelve compartments
even though none of them appear in `$CAPTURE` any more). The checked-in
original (`bc_mrgsolve_model.R`) was left untouched and still carries this
defect.

## Verification

**Method.** Both the original's own model code (with only the `$CAPTURE`
build-compat fix above applied, so it would compile at all) and
`bc_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
(`POST /model_manifest`, `POST /run_simulation`) at `http://localhost:8007`,
spaced ~2s apart per request. The extracted refactored DSL was confirmed
byte-identical to the `_refactored.R`'s embedded quoted string before
running.

All six of the original file's own dosing scenarios were run, exactly as
defined in its own `ev()`/`run_scenario()` calls (dose amounts, `ii`,
`addl`, `time`, and `param_override` all copied verbatim, translated to the
API's dosing-record convention; compartment numbers adjusted only for the
1-based shift caused by dropping `PERI_PALBO` in the refactored model —
e.g. `GUT_LETRO` is compartment 4 in the original, compartment 3 in the
refactored file), over the model's own 1-year horizon (`end = 8760`,
`delta = 24`, daily output):

1. **Scenario 1** — Letrozole monotherapy, 2.5 mg daily oral
2. **Scenario 2** — Palbociclib (125 mg daily-equivalent) + Letrozole (PALOMA-2)
3. **Scenario 3** — "Ribociclib" (150 mg palbociclib-equivalent) + Letrozole (MONALEESA-2)
4. **Scenario 4** — Trastuzumab (560 mg load, 420 mg q3w) + Docetaxel, with
   the original's `kHER2`/`ER_SIGNAL_0`/`HER2_SIGNAL_0` overrides (CLEOPATRA)
5. **Scenario 5** — Pembrolizumab + Chemo, with the original's
   `kCD8_recruit`/`kCD8_kill`/`kprol`/`ER_SIGNAL_0`/`HER2_SIGNAL_0`
   overrides (KEYNOTE-522, TNBC)
6. **Scenario 6** — Olaparib (300 mg q12h), with the original's
   `kprol`/`EC50_CDK` overrides (renamed `EC50_PALBO` for the refactored
   run) (OlympiAD, BRCAm)

Every `$CAPTURE`d output was compared point-by-point across the full
367–368-point time grid for all six scenarios: `Cp_palbo`/`C_PALBO`,
`Cp_letro`/`C_LETRO`, `Cp_tras`/`C_TRAS`, `Cp_olap`/`C_OLAP`, `TUMOR`, `CSC`,
`ER_SIGNAL`, `CDK46_ACT`, `HER2_SIGNAL`, `PD_L1`, `CD8_EFF`, `TREG`,
`E2_PLASMA`, `AROMATASE`, `Ki67`, `CA153`, `TGR`, `response`.

**Result: exact match, max abs diff = 0.0 for every output, every
scenario, no solver timeouts or step-count issues encountered at the
model's own full 1-year horizon.** This is the expected outcome per the
guide's tolerance rule for pure structural reorganization with no
Hill-fitting — since every renamed effect term computes the identical
arithmetic on the identical (renamed) inputs, and the two dropped
quantities (`PERI_PALBO`, `Cp2_tras`) were confirmed dead, there was no
source of numerical divergence to begin with.

## Anything else worth flagging

- `kAI_aro`, `kPDL1_ind`, `EC50_PD1`, and `kKi67_on` are declared in the
  original's `$PARAM` block with descriptive comments implying they drive
  aromatase inhibition, PD-L1 induction, anti-PD1 effect, and Ki-67
  response respectively, but **none of them is referenced anywhere in
  `$ODE`** — the corresponding dynamics all use hardcoded magic-number
  rate constants instead (e.g. `dxdt_AROMATASE` uses `0.01` twice, not
  `kAI_aro`). These are outside the scope of this refactor (none belongs to
  LETRO/OLAP/PALBO/TRAS's own PK/PD block) and were left exactly as the
  original had them — unused, unrenamed, untouched.
- Compartment numbering is not preserved for anything after `CENT_PALBO`,
  because `PERI_PALBO` was removed (see above) — every compartment from
  `GUT_LETRO` onward shifts down by one position relative to the original.
  Any external code addressing this model's compartments by 1-based number
  rather than by name is affected; the R script's own `ev()` calls (which
  address by name) and the qspserver API calls used for verification (which
  were given the corrected numbers explicitly) are unaffected.
- The R-side scenario/plotting code in `bc_mrgsolve_model_refactored.R` was
  updated only where it referenced a renamed identifier: the `EC50_CDK`
  override in Scenario 6 (`EC50_PALBO` now) and the `Cp_palbo`/`Cp_letro`
  column references in Plot 5 (`C_PALBO`/`C_LETRO` now). All dosing
  `ev()` calls, all other `param_override` lists, and all six scenario
  definitions are otherwise identical to the original. The summary `cat()`
  block's compartment count was corrected from 22 (already inaccurate in
  the original, which has 21) to 20, reflecting the `PERI_PALBO` removal —
  a cosmetic print statement, not a gated/verified output.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`breast-cancer | LETRO`, `breast-cancer | OLAP`, `breast-cancer | PALBO`,
and `breast-cancer | TRAS` rows.
