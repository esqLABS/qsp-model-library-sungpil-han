# Refactor notes — `hypoparathyroidism/hypopt_mrgsolve_model.R`

Compounds refactored: **Encaleret (ENC)** and **Thiazide / hydrochlorothiazide
(THZ)** — the two rows in `driver-patches/data/compound_perturbation_census.md`
classified `hypoparathyroidism | Encaleret` and `hypoparathyroidism |
Thiazide (THZ)`, both "Redirect concentration inside #define macro". Every
other compound in this file (oral calcium, calcitriol, alfacalcidol,
cholecalciferol, the three exogenous-PTH molecules, IV calcium gluconate,
magnesium) and every disease-side compartment/equation is byte-for-byte
untouched.

## Why "inside #define macro" needed more than a rename

Both compounds' concentrations were computed as `$GLOBAL` C-preprocessor
macros (`#define CTHZ (THZC/VTHZ*1000.0)`, `#define CENC (ENCC/VENC*1000.0)`),
used both inside `$ODE` (feeding `FRRAW`/`FRCA`/`UCARATE`) and inside `$TABLE`
(feeding `FRCAPCT`/`FEXCA`). Per the guide's qspserver-compatibility
requirement 4, `C_THZ`/`EFFECT_THZ`/`C_ENC`/`EFFECT_ENC_*` must be
`$CAPTURE`d so they are discoverable outputs. That collides with keeping
them as `#define` macros: mrgsolve auto-declares a same-named local for
every `$CAPTURE`d identifier, and a `#define` of that exact identifier would
text-substitute mrgsolve's own generated declaration too (the preprocessor
does not know the difference between the user's code and mrgsolve's
generated code) — this is a stricter version of the collision already
documented for `distal-renal-tubular-acidosis` (`drta_refactor_notes.md`,
"no `double` on `$ODE`-scoped capture names"; that file's collision was
about an explicit `double NAME = ...` redeclaration, this one is about the
preprocessor mangling mrgsolve's own auto-generated declaration).

Resolution, following the `#define`-macro guidance in this task and the
drta precedent for the general mechanism: `C_THZ`, `EFFECT_THZ`, `C_ENC`,
`EFFECT_ENC_SETPT` and `EFFECT_ENC_RENAL` were pulled out of `$GLOBAL`
entirely and are no longer macros at all. They are:
- assigned as bare statements (no `double`) inside `$ODE`, relying on
  mrgsolve's own auto-declaration from being listed in `$CAPTURE`, and
- recomputed as ordinary `double` locals inside `$TABLE`, directly from
  current state — the same pattern this file already uses for every other
  `$TABLE` quantity (e.g. `double CATOTAL = CATOT;`), just written out in
  full instead of via a macro, since these five specific names can no
  longer be macros.

Both were confirmed to compile via `/model_manifest` (see Verification).
The two `$GLOBAL` macros that *consume* these values, `SETPTE` (CaSR
set-point shift) and `FRRAW` (raw fractional Ca reabsorption), are
unchanged in shape — they just read `EFFECT_ENC_SETPT` / `EFFECT_THZ` +
`EFFECT_ENC_RENAL` instead of the old `EMXENC*CENC/(EC50ENC+CENC)` /
`FRTHZM + FRENCM` inline expressions. No other macro, compartment, or
equation in the file was touched.

`$ODE` was reordered (not just edited in place) so the new thiazide and
encaleret blocks — each compound's own compartments, concentration, and
effect term(s), i.e. "its own clearly-delimited block" per the guide's
requirement 1 — sit immediately after the gut-lumen line and before
`dxdt_CAE`, which is the first equation that reads `EFFECT_THZ`/
`EFFECT_ENC_RENAL` (via `UCARATE`→`FRCA`→`FRRAW`). Unlike a `#define`
macro (pure text substitution, order-independent), these are now real
sequential assignments, so they must execute before their first use in the
same `$ODE` call — this is the one place order genuinely matters in this
rewrite. `dxdt_PTHE` (which needs `EFFECT_ENC_SETPT` via `SETPTE`) comes
later in the block, so it is unaffected either way.

## Archetype

**Both compounds: archetype 3, without a peripheral compartment** (depot +
central, linear elimination, first-order absorption with an mrgsolve
dosing-level bioavailability fraction) — a straight rename, not a
reparameterization; the original already used `KA`/`CL`/`V` (not
micro-constants), so no `k10→CL/V1`-style conversion was needed.

- **Thiazide**: `GUT_THZ`/`CENT_THZ` (was `THZDEP`/`THZC`), `KA_THZ`
  (`KATHZ`), `F_THZ` (`FTHZ`), `CL_THZ` (`CLTHZ`), `V1_THZ` (`VTHZ`), all
  values unchanged.
- **Encaleret**: `GUT_ENC`/`CENT_ENC` (was `ENCDEP`/`ENCC`), `KA_ENC`
  (`KAENC`), `F_ENC` (`FENC`), `CL_ENC` (`CLENC`), `V1_ENC` (`VENC`), all
  values unchanged.

Bioavailability is applied exactly as the original did — via mrgsolve's own
`F_<CMT>` dosing-bioavailability mechanism (`F_GUT_THZ = F_THZ;` /
`F_GUT_ENC = F_ENC;` in `$MAIN`, renamed from `F_THZDEP = FTHZ;` /
`F_ENCDEP = FENC;`), not by multiplying `F` into the absorption term inside
the ODE (the alternative the guide's archetype-3 example shows). Both
representations give identical central-compartment trajectories for a pure
depot→central→elimination chain with no other loss from the depot (only the
depot's own reported amount would differ by a factor of `F` between the two
representations); the original's own mechanism was kept so the refactor is
an exact rename with zero arithmetic change anywhere, not a second,
functionally-equivalent-but-differently-scaled reparameterization.

## Hill interface

### Thiazide — rename, not a fit

The original computes `FRTHZM = EMXTHZ*CTHZ/(EC50THZ+CTHZ)`, already
exactly `Emax·C/(EC50+C)` (implicit `Emax` from `EMXTHZ`, implicit `gamma=1`,
no separate ceiling constant). Per the guide's rename rule:

```
EMAX_THZ  = 0.012   // was EMXTHZ, unchanged
EC50_THZ  = 25.0    // was EC50THZ, unchanged (ng/mL)
GAMMA_THZ = 1.0      // no exponent in the original -> 1
EFFECT_THZ = EMAX_THZ*pow(C_THZ, GAMMA_THZ)
            / (pow(EC50_THZ, GAMMA_THZ) + pow(C_THZ, GAMMA_THZ));
```

`EFFECT_THZ` replaces the old `FRTHZM` term verbatim inside `FRRAW`. No fit
was needed or performed.

### Encaleret — rename, but with two named effect terms instead of one

Encaleret is pharmacologically a single calcilytic (CaSR negative allosteric
modulator) but the original genuinely models **two separate downstream
actions of the same concentration**, sharing one potency (`EC50ENC=100`)
but two different maximal effects:

1. `SETPTE = SETPT*(1.0 + EMXENC*CENC/(EC50ENC+CENC))` — rightward shift of
   the parathyroid CaSR set-point (`EMXENC=0.40`), consumed only by
   `SECFRAC`→`PTHSEC`→`dxdt_PTHE`.
2. `FRENCM = EMXENCK*CENC/(EC50ENC+CENC)` — direct increment to TAL
   fractional Ca reabsorption (`EMXENCK=0.022`), consumed only by `FRRAW`.

These are not "one drug on a combined multi-drug expression" (the case the
guide's "combine only where disease equations actually use them" rule
targets) — they are one drug's two real, independent physiological actions,
with different ceilings, feeding two unrelated equations. Collapsing them
into a single `EFFECT_ENC` would either drop one of the two real actions or
silently sum a dimensionless set-point-shift fraction with a
reabsorption-fraction increment, which are not the same quantity. Per the
guide's "None of these fit" fallback ("a clean, non-standard structure
beats a standard structure that's wrong"), this file uses **two** named
effect terms sharing one concentration and one potency:

```
EC50_ENC = 100.0  // was EC50ENC, unchanged (ng/mL), shared by both sites
GAMMA_ENC = 1.0    // no exponent in the original -> 1
EMAX_ENC_SETPT = 0.40   // was EMXENC, unchanged
EMAX_ENC_RENAL = 0.022  // was EMXENCK, unchanged
EFFECT_ENC_SETPT = EMAX_ENC_SETPT*pow(C_ENC, GAMMA_ENC)
                  / (pow(EC50_ENC, GAMMA_ENC) + pow(C_ENC, GAMMA_ENC));
EFFECT_ENC_RENAL = EMAX_ENC_RENAL*pow(C_ENC, GAMMA_ENC)
                  / (pow(EC50_ENC, GAMMA_ENC) + pow(C_ENC, GAMMA_ENC));
```

`SETPTE` now reads `EFFECT_ENC_SETPT`; `FRRAW` now reads
`EFFECT_ENC_RENAL` in place of the old `FRENCM`. Both are pure renames of
already-Hill-shaped ratios — no fit was needed or performed for either.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted DSL
(`http://localhost:8007`):

- `KA_THZ`, `F_THZ`, `V1_THZ`, `CL_THZ`, `EMAX_THZ`, `EC50_THZ`,
  `GAMMA_THZ`, `KA_ENC`, `F_ENC`, `V1_ENC`, `CL_ENC`, `EMAX_ENC_SETPT`,
  `EMAX_ENC_RENAL`, `EC50_ENC`, `GAMMA_ENC` all appear in the manifest's
  `parameters` (fixed values, unchanged from the original).
- `C_THZ`, `EFFECT_THZ`, `C_ENC`, `EFFECT_ENC_SETPT`, `EFFECT_ENC_RENAL`
  are state-dependent (derived from `CENT_THZ`/`CENT_ENC`), so — per the
  same reasoning `distal-renal-tubular-acidosis` used for `C_HCTZ` — they
  cannot be `$PARAM` defaults. All five appear in the manifest's
  `outputPaths` (via `$CAPTURE`), confirmed discoverable.
- No separate `.cpp` extraction step was needed beyond the mechanical one:
  `hypopt_mrgsolve_model.R` is the `hypopt_code <- '...'; mcode_cache(...)`
  wrapper pattern the guide's extraction step targets. The quoted DSL block
  was extracted verbatim from both files for every call below; no
  character of the DSL itself was altered by the extraction.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`), comparing the untouched original against
`hypopt_mrgsolve_model_refactored.R`.

**No upstream defect was hit.** Unlike several prior files in this fork,
`hypopt_mrgsolve_model.R` compiled cleanly on the first `/model_manifest`
call with no workaround, and stayed compiled through every scenario below
— nothing was logged to `translations/UPSTREAM_ISSUES.md`.

**`/model_manifest`**: both the original and refactored DSL compiled
successfully; the refactored manifest additionally lists all 15 renamed
parameters above in `parameters` and `GUT_THZ`, `CENT_THZ`, `GUT_ENC`,
`CENT_ENC`, `C_THZ`, `EFFECT_THZ`, `C_ENC`, `EFFECT_ENC_SETPT`,
`EFFECT_ENC_RENAL` in `outputPaths`.

**`/run_simulation`, thiazide**: reproduced the dosing composition of the
original's own scenario 4 (`"기존 치료 + 티아지드"`, `S[[4]]`) — post-surgical
HypoPT (`PTMASS=0.05, SETPT=1.20`), oral calcium 500 mg q24h into `GUTCA`
(cmt 1), calcitriol 0.25 ug q12h into `CTRDEP` (cmt 14), cholecalciferol
1000 IU q24h into `D3DEP` (cmt 18), and hydrochlorothiazide 25 mg q12h into
`THZDEP`/`GUT_THZ` (cmt 20, unchanged index — the rename does not move any
compartment's position in `$CMT`) — run for the scenario's own full
duration, **180 days** (`end=4320h`), at `delta=24h` (185 points), no
shortening needed: this model did not hit the API's default `maxsteps`
budget at this duration. A finer-resolution spot check (`delta=2h` over the
first 240h, 125 points) was also run and matched.

**`/run_simulation`, encaleret**: reproduced the dosing composition of the
original's own scenario 11 (`"ADH1 + 엔칼레렛"`, `S[[11]]`) — ADH1
(`PTMASS=1.00, SETPT=0.78, CASRGN=1.50, CASROFF=0.012`), encaleret 60 mg
q12h into `ENCDEP`/`GUT_ENC` (cmt 22) — run for the scenario's own full
duration, 180 days, at `delta=24h` (185 points), plus the same `delta=2h`/
240h spot check.

Both scenarios started from the model's own `$MAIN`-derived default initial
state (the healthy steady-state formulas, e.g. `CAE_0 = CAT0*VECF*10.0`)
with the scenario's disease parameters applied via the request's
`parameters` field, rather than from the R-side `hypopt_baseline()`
pre-equilibrated state (that helper runs the model forward 400 days with no
therapy first, entirely in R, which the API has no path to reproduce). This
means the compared trajectories are transiently relaxing toward the
diseased steady state over the run rather than starting from it — expected
and disclosed, and irrelevant to what this check verifies (that the
refactored DSL computes identically to the original given identical
inputs), which does not depend on where the run starts from.

**Result: exact match, not just near-exact.** Every compared output —
all 25 disease-side `$CAPTURE`/state quantities shared between both files
(`CATOTAL`, `CACORRD`, `CAIONMM`, `PISER`, `MGSER`, `PTHPG`, `PTHENDO`,
`D125PG`, `D25NG`, `TMPPERG`, `FRCAPCT`, `FEXCA`, `UCA24`, `UPI24`,
`CAABS`, `FABSPCT`, `CAPROD`, `GUTCA`, `CAE`, `CABONE`, `CAMIN`, `PIE`,
`PIBONE`, `MGE`, `PTHE`), plus each compound's own PK states under their
old/new names (`THZDEP`/`GUT_THZ`, `THZC`/`CENT_THZ` for the thiazide run;
`ENCDEP`/`GUT_ENC`, `ENCC`/`CENT_ENC` for the encaleret run) — matched with
**maximum absolute difference exactly 0.0** across every time point in
every run (180-day runs at `delta=24h` and the 240h `delta=2h` spot checks,
four run-pairs total). This is a stronger result than the API's own
4-decimal JSON rounding floor documented in other refactors in this fork
(e.g. `distal-renal-tubular-acidosis`) — here the numbers are bit-identical
as returned, consistent with this being a pure rename with no arithmetic
change anywhere (no `CL=k*V`-style reparameterization was needed, unlike
`drta`'s HCTZ). `EFFECT_THZ` ranged 0–0.009, `EFFECT_ENC_SETPT` ranged
0–0.146, and `EFFECT_ENC_RENAL` ranged 0–0.008 over the two dosed windows,
confirming the Hill arithmetic is active and not just trivially zero.

## Anything else flagged

- No upstream defect found or logged.
- No compound other than Encaleret and Thiazide was touched. Oral calcium,
  calcitriol, alfacalcidol, cholecalciferol, the three exogenous-PTH
  molecules (rhPTH(1-84)/teriparatide/TransCon), IV calcium gluconate, and
  magnesium repletion — their compartments, parameters, macros, and dosing
  helpers — are byte-identical to `hypopt_mrgsolve_model.R`.
- The R-side dosing helpers `hypopt_thiazide()` and `hypopt_encaleret()`
  were updated to dose into `"GUT_THZ"`/`"GUT_ENC"` (was `"THZDEP"`/
  `"ENCDEP"`) so the refactored sibling's own R scenario-driver code stays
  internally consistent with the renamed compartments; this is a rename of
  the dosing target only, not a behavioral change (same compartment,
  same 1-based index, same dose amounts/timing).

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`hypoparathyroidism | Encaleret` and `hypoparathyroidism | Thiazide (THZ)`
rows.
