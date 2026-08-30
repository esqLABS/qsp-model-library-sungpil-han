# Rosacea (ros) — PK/PD refactor notes

**Scope.** One compound only: **isotretinoin** (stem `ISO`), the compound
this file's row in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
classifies as *"Redirect concentration (clean single site)"*. The census's
`ISO` abbreviation alone doesn't disambiguate — the model was read to
confirm: this file's `ISO` compartments/parameters (`ISOG`/`ISOP`, `KAISO`,
`VISO`, `CLISO`, `ISOEC50`, `ISOEMAX`, `ISOAIC50`, `CISO`) are unambiguously
**isotretinoin** (oral retinoid, `dose_isotretinoin()` regimen builder,
docstring scenario "S12 isotretinoin 20 mg/day", calibration anchor
"isotretinoin 20 mg/day -71% [lesion reduction], published anchor -70 to
-90%"), not ivermectin (which is separately modelled in this same file under
the stem `IVM`, e.g. `IVMSK`/`IVMFO`).

Every other compound modelled in `ros_mrgsolve_model.R` — ivermectin (`IVM`),
metronidazole (`MTZ`), azelaic acid (`AZA`), brimonidine (`BRM`),
oxymetazoline (`OXY`), minocycline (`MIN`), doxycycline (`DOX`) — is
completely untouched in `ros_mrgsolve_model_refactored.R`: same compartment
names, same parameter names, same values, same ODEs.

## Archetype

**Archetype 3 minus peripheral** (depot + central, linear elimination, no
peripheral compartment, no bioavailability term in the original):

```
GUT_ISO  (renamed from ISOG)   -- oral depot, mg
CENT_ISO (renamed from ISOP)   -- central amount, mg
KA_ISO   (renamed from KAISO)  -- absorption rate, 1/day
V1_ISO   (renamed from VISO)   -- apparent volume, L
CL_ISO   (renamed from CLISO)  -- apparent clearance, L/day

dxdt_GUT_ISO  = -KA_ISO * GUT_ISO;
dxdt_CENT_ISO =  KA_ISO * GUT_ISO - (CL_ISO / V1_ISO) * CENT_ISO;
C_ISO = CENT_ISO / V1_ISO;   -- exposed concentration, mg/L
```

The original had no bioavailability parameter for isotretinoin (full-dose
absorption into `ISOP`/`CENT_ISO`) — none was invented; `F_ISO` simply does
not exist, matching the original exactly.

## Two effects, one concentration

Isotretinoin genuinely has **two independent disease-facing mechanisms** in
the original, both driven by the same single concentration:

1. **Sebosuppression** — `EISO = ISOEMAX * CISO / (ISOEC50 + CISO)`, feeding
   `dxdt_SEB` (reduces the sebaceous-gland target, shrinking the Demodex
   follicular carrying-capacity habitat indirectly).
2. **Innate/TLR2-KLK5 damping** — `FISOA = 1 / (1 + CISO / ISOAIC50)`, a
   multiplicative fraction-remaining term consumed in two places:
   `dxdt_KLK` (`KLKD * FAZAK * FDOXK * FISOA`) and the TLR2 drive `QTL`
   (`... * FIVMT * FISOA`).

Both are kept as **separate named `EFFECT_ISO_*` terms**, not collapsed into
one. The guide's caution against burying an effect "inside a combined
multi-drug expression" is about mixing *different drugs'* effects together
(the reason this corpus keeps `EFFECT_TCZ`, `EFFECT_ADA`, etc. distinct from
each other) — it is not a mandate to compress one compound's own distinct
mechanisms into a single term when the original itself does not. This same
file's doxycycline block has an even more pronounced version of the same
pattern (five separate `FDOX*` inhibition fractions: `FDOXK`, `FDOXM`,
`FDOXI`, `FDOXO`, `FABX`), left untouched here since doxycycline is out of
scope for this row.

`C_ISO = CENT_ISO / V1_ISO` is the single exposed concentration both
`EFFECT_ISO_SEB` and `EFFECT_ISO_INN` are computed from — satisfying "exactly
one concentration variable that PD equations read" even though there are two
downstream effect terms.

## Hill interface: renames, not fits

Both of isotretinoin's effect terms were already plain concentration ratios
in the original — no ODE-solved kinetics to approximate — so both are exact
renames, not fits:

- **`EFFECT_ISO_SEB`** = `EMAX_ISO_SEB * C_ISO^GAMMA_ISO_SEB /
  (EC50_ISO_SEB^GAMMA_ISO_SEB + C_ISO^GAMMA_ISO_SEB)`, with
  `EMAX_ISO_SEB = 0.85` (renamed from `ISOEMAX`), `EC50_ISO_SEB = 0.30`
  (renamed from `ISOEC50`), `GAMMA_ISO_SEB = 1` (new — the original had no
  explicit Hill coefficient; its own ratio is already `Emax*C^1/(EC50^1+C^1)`
  with an implicit `gamma = 1`). Algebraically identical to the original's
  `ISOEMAX * CISO / (ISOEC50 + CISO)`.
- **`EFFECT_ISO_INN`** = `EMAX_ISO_INN * C_ISO^GAMMA_ISO_INN /
  (EC50_ISO_INN^GAMMA_ISO_INN + C_ISO^GAMMA_ISO_INN)`, with
  `EC50_ISO_INN = 0.30` (renamed from `ISOAIC50`), `GAMMA_ISO_INN = 1` (new,
  implicit in the original), and **`EMAX_ISO_INN = 1.0`** — new, but *not* a
  free choice: it is forced by the original's own algebra, not invented
  independently. The original's `FISOA = 1/(1+CISO/ISOAIC50)` is exactly
  `1 - [1 * CISO/(CISO+ISOAIC50)]`, i.e. an implicit-`Emax=1` Hill
  fractional-*inhibition* term. Downstream code computes
  `FISOA = 1.0 - EFFECT_ISO_INN`, which is algebraically identical to the
  original's `FISOA` for every value of `C_ISO` (confirmed by the exact
  verification match below, not just by hand-algebra).

No `nls()` fitting was needed or performed.

## `$PARAM` vs `$GLOBAL` for `C_ISO`/`EFFECT_ISO_*`

Per the guide's qspserver compatibility requirement #2, these should ideally
live in `$PARAM` (with a `= 0` default) for direct `/model_manifest`
discoverability. This was not done, for the reason already documented in
several sibling refactors (e.g. `copd/copd_mrgsolve_model_refactored.R`,
`breast-cancer/bc_refactor_notes.md`): mrgsolve 2.0.1 compiles `$PARAM`
members as **read-only references** inside `$ODE`, so a value that must be
recomputed every timestep from state (`C_ISO`, `EFFECT_ISO_SEB`,
`EFFECT_ISO_INN` — all of which read `CENT_ISO`) cannot also be declared in
`$PARAM`. Instead, all three are predeclared as `double`s in `$GLOBAL` and
listed in `$CAPTURE`: confirmed present in `/model_manifest`'s `outputPaths`
for the refactored model (and absent from its `parameters` list, as
expected for a `$GLOBAL`-declared value) — see Verification below. The
existing `$TABLE` output `CISOO` (isotretinoin plasma concentration, mg/L)
is kept, unchanged in name and value, now defined as a plain alias
(`double CISOO = C_ISO;`) rather than recomputing `ISOP / VISO` a second
time.

## Build compatibility

**No pre-existing build defect was found.** `ros_mrgsolve_model.R`'s own DSL
compiled cleanly under mrgsolve 2.0.1 via the qspserver `mrgsolve_api`'s
`/model_manifest` endpoint, as-is, with no syntax changes. Nothing is logged
to `translations/UPSTREAM_ISSUES.md` for this file.

## Verification

**Method.** The quoted `ros_code <- '...'` DSL blocks were extracted
verbatim from both `ros_mrgsolve_model.R` (original) and
`ros_mrgsolve_model_refactored.R` (this file) as bare mrgsolve DSL text, and
run through the qspserver `mrgsolve_api` container (`POST /model_manifest`,
`POST /run_simulation`) at `http://localhost:8007`, requests spaced ~2s
apart. `$INIT`-declared compartment order is unchanged by the rename (12th
compartment is `ISOG`/`GUT_ISO`, 13th is `ISOP`/`CENT_ISO` in both files —
confirmed by listing both files' `$INIT @annotated` blocks in order).

Both models' own `/model_manifest` compiled successfully with no errors.
`/model_manifest` on the refactored model confirms: `KA_ISO`, `V1_ISO`,
`CL_ISO`, `EC50_ISO_SEB`, `EMAX_ISO_SEB`, `GAMMA_ISO_SEB`, `EC50_ISO_INN`,
`EMAX_ISO_INN`, `GAMMA_ISO_INN` all present in `parameters`; `GUT_ISO`,
`CENT_ISO`, `C_ISO`, `EFFECT_ISO_SEB`, `EFFECT_ISO_INN` all present in
`outputPaths`; `C_ISO`/`EFFECT_ISO_SEB`/`EFFECT_ISO_INN` correctly *absent*
from `parameters` (as expected for `$GLOBAL`-declared values, per the
section above).

**Scenario.** The original's own `S12` scenario dosing (`dose_isotretinoin`:
20 mg/day oral, `ii = 1`, `addl = 111` → 112 daily doses) was reproduced via
the API's `dosing` field (`amt = 20, cmt = 12` [1-based, `GUT_ISO`/`ISOG`],
`ii = 1, addl = 111`), over the model's own 16-week window (`end = 112,
delta = 1`, matching `$SET` and `S12`'s `days = 7*16`). Phenotype parameters
were set to the original's own `PPR-severe` row from `ros_phenotypes()`
(`SPROT=3.6, SNEUR=2.0, SMITE=16.0, SFIBR=2.0, ANDROG=1.1, ADBREV=0.8,
TRIGB=0.45, UVLOAD=0.35`), matching `S12`'s phenotype.

One simplification, disclosed: the run started from the model's plain
`$INIT` defaults rather than reproducing the R-side 10-year phenotype
burn-in (`ros_init_at()`). The stateless API's `SimRequest` has no
initial-condition-override field, so the burn-in (which only rewires
*non-ISO* state variables to their chronic phenotype steady state — `GUT_ISO`/
`CENT_ISO` start at 0 in the `$INIT` defaults regardless of burn-in, since no
isotretinoin is present during the burn-in itself) isn't reachable through
it. Since `EFFECT_ISO_SEB`/`EFFECT_ISO_INN`/`C_ISO` are pure functions of
`C_ISO` alone, and every other output's dependence on isotretinoin is only
*through* those two terms, an identical starting state for both models is
what makes this an apples-to-apples verification of the compound's own
refactor — it is not a claim to reproduce the full clinical `S12` phenotype
trajectory (which would require burn-in) but a full 16-week integration.

**Result.** All 36 shared `$CAPTURE`/state outputs compared
(`SEB KLK LL37 TLR2 IL1B MMP9 ROS MC NEU TH17 IL17 VEGF VDEN BARR DEMO CEA
PSA IGA ILC TELSC PHYGR FLFREQ STING DLQI OSDI ERYIDX ERYS1 ERYS2 TRIGO
CDOXO CISOO OCC_A2 TE_A2 TE_MMP TE_ABX TE_KILL`) over all 114 output rows:

- 33 of 36 outputs matched **exactly** (max abs diff `0.0`), including
  `SEB`/`KLK`/`TLR2` (the state variables isotretinoin's two effects feed
  into directly), `DEMO`/`BARR` (indirect downstream), and `CISOO` itself.
- 3 outputs (`CEA`, `DLQI`, `ERYIDX`) showed a **max abs diff of exactly
  `1.0e-4`**, at 1, 2, and 4 isolated time points respectively (out of 114),
  with every other time point exact. `1.0e-4` is the API's own JSON output
  rounding granularity (values serialize to 4 decimal places); this is
  floating-point-scale noise from `pow(C_ISO, 1.0)` vs. a plain
  multiplication (the refactor introduces `pow(x, GAMMA)` calls with
  `GAMMA = 1`, which some libm implementations do not evaluate bit-identical
  to `x`), occasionally tipping the last displayed digit at isolated
  timesteps — not a structural or numeric mismatch. Per the guide's
  tolerance for a pure structural reorganization (Archetypes 1-3), this
  counts as a match, not a bug: the deviation is at the serialization
  rounding floor itself, not accumulated drift.
- `GUT_ISO`/`ISOG` and `CENT_ISO`/`ISOP` (the renamed PK compartments
  themselves) matched **exactly** (max abs diff `0.0`) across all 114 rows.

No mismatch beyond this floating-point-scale noise was found. Nothing was
loosened or adjusted to force a pass.
