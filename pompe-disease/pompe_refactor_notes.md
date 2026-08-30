# Refactor notes — `pompe-disease/pompe_mrgsolve_model.R`

Four compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md` (all four classified
"Redirect concentration (clean single site)"): **Alglucosidase alfa
(ALGLU)**, **Avalglucosidase alfa (AVAL)**, **Cipaglucosidase alfa
(CIPA)**, **Miglustat (MIG)**. The file also models AAV9-hGAA gene therapy
(`AAV_X`) and a Rituximab-based ITI regimen (`RTX_C`) — both out of scope
(no census row for either) and left completely untouched.

## Checking the actual scenarios, per the task brief

The three enzyme-replacement therapies (ALGLU, AVAL, CIPA) are real-world
mutually-exclusive treatment arms — a patient receives one ERT, not several
at once. The original file's own `pompe_run()` scenario switch confirms
this structurally: `"alglu"` doses only `ALGLU_C`, `"aval"` doses only
`AVAL_C`, and `"cipa_mig"` doses `CIPA_C` **and** `MIG_A` together — there
is no scenario that doses two ERTs simultaneously, and no scenario doses
Miglustat without Cipaglucosidase alfa. This matches the task brief's
expectation exactly: Miglustat is modeled strictly as Cipaglucosidase
alfa's co-administered enzyme stabilizer (see `EFFECT_MIG`/`C_CIPA_STAB`
below), never dosed alone. Nothing about the four compounds' own PK/PD
blocks is mutually exclusive in the *code* (all three ERTs can be
dosed as easily as one — the mutual exclusivity is a clinical/scenario
convention, not a structural constraint), so the refactor treats all four
compounds independently and lets the scenario definitions (preserved
faithfully in the R wrapper) continue to reflect that real-world usage
pattern.

## Archetype per compound

- **Alglucosidase alfa (ALGLU) — Archetype 2** (no depot, two
  compartments, linear elimination), matched exactly, no deviation:
  `dxdt_CENT_ALGLU = -(CL_ALGLU+Q_ALGLU)/V1_ALGLU*CENT_ALGLU +
  Q_ALGLU/V2_ALGLU*PERI_ALGLU`, `dxdt_PERI_ALGLU = Q_ALGLU/V1_ALGLU*
  CENT_ALGLU - Q_ALGLU/V2_ALGLU*PERI_ALGLU`. IV dosing directly into the
  central compartment (`KA_ALGLU = 0`, unused, consistent with the
  original's own comment "Not used (IV infusion)").
- **Avalglucosidase alfa (AVAL) — Archetype 2**, identical structure to
  ALGLU, same shape.
- **Cipaglucosidase alfa (CIPA) — Archetype 2**, identical structure to
  ALGLU/AVAL. Its own PK block is unchanged by the fact that Miglustat
  modifies its *effective potency downstream* (see Hill interface below)
  — the PK compartments and elimination kinetics are untouched by
  Miglustat co-administration, only the disease-facing effect term is.
- **Miglustat (MIG) — Archetype 3 variant** (depot + central, linear, no
  peripheral compartment) — exactly the variant the guide's own archetype
  3 anticipates ("drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a
  2-compartment-no-depot variant" — the symmetric case, keep the depot,
  drop the peripheral leg, applies the same way): `dxdt_GUT_MIG =
  -KA_MIG*GUT_MIG`, `dxdt_CENT_MIG = KA_MIG*GUT_MIG*F_MIG -
  (CL_MIG/V1_MIG)*CENT_MIG`.

No TMDD/receptor-binding kinetics anywhere in this file for any of the four
compounds — every PK block is plain linear compartmental disposition, and
every disease-facing effect (see below) is a Michaelis–Menten/Hill ratio on
plasma concentration, not a receptor-occupancy state solved by its own ODEs.

## Renaming (values unchanged from the original)

| Original | Refactored |
|---|---|
| `ALGLU_C` / `ALGLU_P` | `CENT_ALGLU` / `PERI_ALGLU` |
| `AVAL_C` / `AVAL_P` | `CENT_AVAL` / `PERI_AVAL` |
| `CIPA_C` / `CIPA_P` | `CENT_CIPA` / `PERI_CIPA` |
| `MIG_A` / `MIG_C` | `GUT_MIG` / `CENT_MIG` |
| `V_MIG` | `V1_MIG` |
| `MIG_STAB` | `EMAX_MIG` |
| `MIG_EC50` | `EC50_MIG` |
| `Cp_alglu` / `Cp_aval` / `Cp_cipa` / `Cp_mig` | `C_ALGLU` / `C_AVAL` / `C_CIPA` / `C_MIG` |
| `Cp_cipa_stab` | `C_CIPA_STAB` |
| `mig_eff` | `EFFECT_MIG` |
| `uptake_alglu` / `uptake_aval` / `uptake_cipa` | `EFFECT_ALGLU` / `EFFECT_AVAL` / `EFFECT_CIPA` |

`KA_ALGLU`, `CL_ALGLU`, `V1_ALGLU`, `Q_ALGLU`, `V2_ALGLU`, `CL_AVAL`,
`V1_AVAL`, `Q_AVAL`, `V2_AVAL`, `CL_CIPA`, `V1_CIPA`, `Q_CIPA`, `V2_CIPA`,
`KA_MIG`, `CL_MIG`, `F_MIG` already matched the naming convention
(compound-specific PK parameters already keyed to the right stem in the
original) and are unchanged. All parameter *values* are copied verbatim —
nothing invented or defaulted.

New: `GAMMA_MIG = 1` (a genuine new `$PARAM`, since it's a fresh shape
constant not derived from any other existing parameter — same treatment as
`EMAX_LETRO`/`GAMMA_TRAS` in the breast-cancer refactor).

**Unused/dead parameters, left untouched:** `M6P_ALGLU`, `M6P_AVAL`,
`M6P_CIPA` are declared in `$PARAM` with descriptive comments ("Mannose-6-P
per mole alglucosidase", "~15-fold higher M6P content") but never
referenced anywhere in `$ODE` or `$TABLE` — the actual M6P-receptor uptake
potency differential between compounds is carried entirely by `RHO_AVAL`/
`RHO_CIPA` (which *are* used, see Hill interface below). This looks like
the author intended `M6P_ALGLU`/`M6P_AVAL`/`M6P_CIPA` to feed `RHO_AVAL`/
`RHO_CIPA` (`RHO_AVAL = 3.0` is in the right ballpark for `M6P_AVAL/
M6P_ALGLU = 37.5/2.5 = 15`, scaled down; `RHO_CIPA = 1.8` vs `M6P_CIPA/
M6P_ALGLU = 25/2.5 = 10` is not as close a match) but never wired the
connection up. Out of scope for a rename — these three parameters already
carry the correct stem and are simply dead code in the original; left
exactly as found.

## Hill interface: renames plus one disclosed derivation choice, no fitting

All four compounds' disease-facing effect terms were already written as
Michaelis–Menten/Hill-shaped ratios in the original (implicit `gamma = 1`
throughout) — no ODE-solved receptor kinetics to approximate anywhere, so
every one of these is a rename or an exact algebraic derivation, confirmed
by the exact verification match below (not merely "close"):

- **Miglustat**: `mig_eff = MIG_STAB * Cp_mig / (MIG_EC50 + Cp_mig)` →
  `EFFECT_MIG = EMAX_MIG * pow(C_MIG, GAMMA_MIG) / (pow(EC50_MIG,
  GAMMA_MIG) + pow(C_MIG, GAMMA_MIG))` with `GAMMA_MIG = 1` added
  explicitly. `pow(x, 1.0) == x` exactly (IEEE-754 guarantee), so this
  introduces no floating-point drift.
- **Alglucosidase alfa**: `uptake_alglu = VMAX_UPT * Cp_alglu / (KM_UPT +
  Cp_alglu) * (1 - ada_block)` → `EFFECT_ALGLU = EMAX_ALGLU *
  pow(C_ALGLU,1) / (pow(EC50_ALGLU,1) + pow(C_ALGLU,1)) * (1 -
  ada_block)`, with `EMAX_ALGLU = VMAX_UPT` and `EC50_ALGLU = KM_UPT`.
- **Avalglucosidase alfa**: `uptake_aval = VMAX_UPT * RHO_AVAL * Cp_aval /
  (KM_UPT/RHO_AVAL + Cp_aval) * (1 - 0.5*ada_block)` → `EFFECT_AVAL =
  EMAX_AVAL * C_AVAL / (EC50_AVAL + C_AVAL) * (1 - 0.5*ada_block)`, with
  `EMAX_AVAL = VMAX_UPT * RHO_AVAL` and `EC50_AVAL = KM_UPT / RHO_AVAL` —
  algebraically identical to the original's ratio.
- **Cipaglucosidase alfa**: `uptake_cipa = VMAX_UPT * RHO_CIPA *
  Cp_cipa_stab / (KM_UPT + Cp_cipa_stab) * (1 - 0.5*ada_block)` →
  `EFFECT_CIPA = EMAX_CIPA * C_CIPA_STAB / (EC50_CIPA + C_CIPA_STAB) * (1 -
  0.5*ada_block)`, with `EMAX_CIPA = VMAX_UPT * RHO_CIPA` and `EC50_CIPA =
  KM_UPT`. `C_CIPA_STAB = C_CIPA * (1 + EFFECT_MIG)` — the
  Miglustat-stabilized concentration Cipaglucosidase alfa's own effect term
  reads — is unchanged in structure from the original's `Cp_cipa_stab`.

No `nls()` fitting was needed or performed for any of the four compounds —
every renamed effect term computes the identical arithmetic on the
identical (renamed) inputs.

### Design decision: `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>` for
ALGLU/AVAL/CIPA are `$ODE`-computed doubles, not independent `$PARAM`
entries — disclosed, not a fitting shortcut

`VMAX_UPT`, `KM_UPT`, `RHO_AVAL`, and `RHO_CIPA` are shared, pre-existing
`$PARAM` entries that feed all three ERT compounds' uptake terms (`VMAX_UPT`
and `KM_UPT` directly parameterize ALGLU's term and are the base that
AVAL's/CIPA's terms scale via `RHO_AVAL`/`RHO_CIPA`). Rather than
hardcoding `EMAX_AVAL = 24.0` (`= VMAX_UPT * RHO_AVAL = 8.0 * 3.0`) and
`EC50_AVAL = 0.2333...` (`= KM_UPT / RHO_AVAL = 0.7 / 3.0`, a
non-terminating decimal) as independent `$PARAM` literals, they are
computed as `double`s in `$ODE` from the original shared parameters every
step:

```
double EMAX_ALGLU = VMAX_UPT;
double EC50_ALGLU = KM_UPT;
double GAMMA_ALGLU = 1.0;
double EMAX_AVAL  = VMAX_UPT * RHO_AVAL;
double EC50_AVAL  = KM_UPT / RHO_AVAL;
double GAMMA_AVAL = 1.0;
double EMAX_CIPA  = VMAX_UPT * RHO_CIPA;
double EC50_CIPA  = KM_UPT;
double GAMMA_CIPA = 1.0;
```

Two reasons, both disclosed here rather than silently chosen: (1) hardcoding
`EC50_AVAL` as a truncated decimal literal would introduce a genuine,
avoidable ~1e-10-scale numerical deviation relative to computing
`KM_UPT/RHO_AVAL` at double precision — exactly the kind of drift the
guide's tolerance rule for a pure structural reorganization warns against,
even though 1e-10 would likely still read as "floating-point-scale" to a
casual reviewer; and (2) a literal copy would silently desynchronize from
`VMAX_UPT`/`KM_UPT`/`RHO_AVAL`/`RHO_CIPA` if those are ever overridden via
the API's `parameters` field for a differently-parameterized run — the
`$ODE`-computed form always tracks whatever value those shared parameters
actually hold at run time, exactly like `C_<STEM>` already does. All nine
are exposed via `$CAPTURE` (not `$PARAM`) for the same reason `C_<STEM>`
and `EFFECT_<STEM>` are — see below. `GAMMA_MIG`, by contrast, is not
derived from any existing parameter (there was no explicit Hill
coefficient in the original at all for Miglustat), so it was added as a
genuine new `$PARAM` entry, no different from `EMAX_LETRO=1`/
`GAMMA_TRAS=1` in the breast-cancer refactor.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>` (and the derived
Hill constants above)

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` (with a `= 0` default) for direct
`/model_manifest` discoverability. This was not done, for the same reason
already documented in the AMD, membranous-nephropathy, breast-cancer, and
x-linked-hypophosphatemia refactors: mrgsolve 2.0.1 compiles `$PARAM`
members as read-only references inside `$ODE`, so a value recomputed every
timestep from state (`C_ALGLU`, `C_AVAL`, `C_CIPA`, `C_MIG`,
`EFFECT_ALGLU`, `EFFECT_AVAL`, `EFFECT_CIPA`, `EFFECT_MIG`, `C_CIPA_STAB`,
and the nine `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>` doubles for
ALGLU/AVAL/CIPA above) cannot also be declared in `$PARAM`. Instead, all of
these are `double`s computed in `$ODE` and listed in `$CAPTURE`: visible in
every simulation's output columns and in `/model_manifest`'s
`outputPaths`, confirmed present (49 total output paths, up from the
original's 22 compartments — 22 compartments + 27 captured quantities).
`GAMMA_MIG`, `EMAX_MIG` (renamed from `MIG_STAB`), and `EC50_MIG` (renamed
from `MIG_EC50`) *are* real, unmodified `$PARAM` entries — they are pure
constants, never reassigned, so the read-only restriction does not apply
to them; confirmed present in `/model_manifest`'s `parameters` list
alongside every PK parameter (`KA_ALGLU`/`CL_ALGLU`/`V1_ALGLU`/`Q_ALGLU`/
`V2_ALGLU`, `CL_AVAL`/`V1_AVAL`/`Q_AVAL`/`V2_AVAL`, `CL_CIPA`/`V1_CIPA`/
`Q_CIPA`/`V2_CIPA`, `KA_MIG`/`CL_MIG`/`V1_MIG`/`F_MIG`, and
`VMAX_UPT`/`KM_UPT`/`RHO_AVAL`/`RHO_CIPA`).

## Pre-existing upstream build defects (fixed syntax-only in the delivered
file, logged as `translations/UPSTREAM_ISSUES.md` #71)

**The original does not compile under mrgsolve 2.0.1 at all**, and even
patched to compile, one of its own `$TABLE` mechanisms turns out to be
silently non-functional. Two independent defects, confirmed via `POST
/model_manifest` and `POST /run_simulation` on the untouched original:

1. **11 `$PARAM @annotated` lines have no description field** (`Q_AVAL`,
   `V2_AVAL`, `CL_CIPA`, `V1_CIPA`, `Q_CIPA`, `V2_CIPA`, `M6P_CIPA`,
   `DIAPH_LOSS`, `DIAPH_GAIN`, `SMWT_MAX`, `SMWT_MIN`), which mrgsolve
   2.0.1's annotated-parameter parser rejects outright (`Error: improper
   annotation format`). Seven of these eleven (`Q_AVAL`/`V2_AVAL`,
   `CL_CIPA`/`V1_CIPA`/`Q_CIPA`/`V2_CIPA`/`M6P_CIPA`) belong to
   Avalglucosidase alfa's and Cipaglucosidase alfa's own PK blocks — in
   scope for this refactor, so fixed as part of it. The remaining four
   (`DIAPH_LOSS`, `DIAPH_GAIN`, `SMWT_MAX`, `SMWT_MIN`) are disease-side,
   out of scope, but the whole file cannot compile until every one of the
   eleven has a description, so they were fixed as a pure build-compat pass
   alongside the seven in-scope ones. A short, accurate description was
   added to each; no value changed.
2. **`$TABLE` ends in 13 bare `capture NAME;` lines** — one bare identifier
   per line, semicolon-terminated, no `$CAPTURE` header anywhere in the
   file. This compiles without error under this mrgsolve build, but is
   silently non-functional: none of the 13 named quantities (`Cp_alglu`,
   `Cp_aval`, `Cp_cipa`, `Cp_mig`, `Cp_rtx`, `SMWT`, `VENT_RISK`,
   `SF36_PCS`, `NTproBNP`, `EF_LV`, `CK`, `tissue_supply`, `ada_block`) is
   actually retrievable via `/run_simulation` — every one fails "not a
   compartment or captured item." Replaced with three `$CAPTURE` header
   lines (grouping the same 13 names, plus the newly-added
   `EFFECT_ALGLU`/`EFFECT_AVAL`/`EFFECT_CIPA`/`EFFECT_MIG`/`C_CIPA_STAB`
   and the nine derived `EMAX`/`EC50`/`GAMMA` doubles); same values, now
   actually exposed.

Both fixes are syntax-only and non-numeric, applied directly to the
delivered `pompe_mrgsolve_model_refactored.R` per the guide's settled
policy for a non-compiling original — never to the checked-in
`pompe_mrgsolve_model.R`, which is untouched and still carries both
defects exactly as written. Confirmed non-numeric by the verification
below: all four scenarios match the (identically defect-patched) original
exactly, max abs diff 0.0.

## Verification

**Method.** Both the original's own model code (with the two build-compat
fixes above applied, so it would compile and its `$TABLE` outputs would
actually be retrievable) and `pompe_mrgsolve_model_refactored.R`'s embedded
DSL were extracted as bare mrgsolve DSL text and run through the qspserver
`mrgsolve_api` container (`POST /model_manifest`, `POST /run_simulation`)
at `http://localhost:8007`, requests spaced ~2s apart. The extracted
refactored DSL was confirmed byte-identical to the `_refactored.R`'s
embedded quoted string before running (and re-confirmed after a later
cosmetic edit to a comment's issue-number cross-reference, which changed
nothing else).

Four scenarios were run, covering all four in-scope compounds and matching
the original file's own `pompe_run()` scenario definitions (dose amounts
and q2w/q24h intervals copied verbatim, translated to the API's
`dosing`/`ii`/`addl` convention; compartment numbers are 1-based and
unchanged by the rename since compartment declaration order is unchanged).
The verification window was shortened from the original scenario runner's
default 3 years to **180 days** (daily-equivalent output every 15 days) —
not because of any solver step-count limit encountered (none was), but
because it is ample to exercise the full q2w ERT dosing cadence (13 doses)
and the daily Miglustat dosing cadence (180 doses) while keeping the
verification light given the API's documented 2-concurrent-job limit:

1. **`no_tx`** — no dosing at all. Sanity check that all 22 compartments
   and every captured quantity hold at (or return to) baseline with all
   four compounds silent.
2. **`alglu`** — Alglucosidase alfa alone, 20 mg/kg (1400 mg for a 70 kg
   patient) q2w into `CENT_ALGLU` (cmt 1), 13 doses.
3. **`aval`** — Avalglucosidase alfa alone, same 20 mg/kg q2w schedule
   into `CENT_AVAL` (cmt 3).
4. **`cipa_mig`** — Cipaglucosidase alfa (20 mg/kg q2w into `CENT_CIPA`,
   cmt 5) **and** Miglustat (195 mg once daily into `GUT_MIG`, cmt 7)
   dosed together — the original's own combination regimen, and the only
   scenario that doses Miglustat at all.

Every shared output was compared point-by-point across the full time grid
for all four scenarios: `Cp_alglu`/`C_ALGLU`, `Cp_aval`/`C_AVAL`,
`Cp_cipa`/`C_CIPA`, `Cp_mig`/`C_MIG`, every PK compartment
(`ALGLU_C`/`CENT_ALGLU`, `ALGLU_P`/`PERI_ALGLU`, `AVAL_C`/`CENT_AVAL`,
`AVAL_P`/`PERI_AVAL`, `CIPA_C`/`CENT_CIPA`, `CIPA_P`/`PERI_CIPA`,
`MIG_A`/`GUT_MIG`, `MIG_C`/`CENT_MIG`), `GAA_M`, `GAA_C`, `GAA_D`,
`GLYC_M`, `GLYC_C`, `GLYC_D`, `HEX4`, `ADA_T`, `LVMI`, `MM_IDX`,
`DIAPH_F`, `FVC_UP`, `SMWT`, `VENT_RISK`, `SF36_PCS`, `NTproBNP`, `EF_LV`,
`CK`, `tissue_supply`, `ada_block`.

**Result: exact match, max abs diff = 0.0 for every output, every
scenario.** This is the expected outcome per the guide's tolerance rule
for pure structural reorganization with rename-only (and exact-derivation)
Hill terms — every renamed/derived quantity computes the identical
arithmetic on the identical (renamed) inputs, so there was no source of
numerical divergence to begin with.

**Sanity check (not required for verification, but confirms the new
outputs are wired correctly):** in scenario 4, `EMAX_CIPA` = 14.4 (=
`VMAX_UPT * RHO_CIPA` = 8.0 × 1.8) and `EC50_CIPA` = 0.7 (= `KM_UPT`) hold
constant throughout as expected for derived constants; `EFFECT_MIG`
equilibrates to ~0.178 as `C_MIG` approaches its daily-dosing steady state
of ~1.55 mg/L (`EMAX_MIG * 1.5515 / (EC50_MIG + 1.5515) = 0.35 * 1.5515 /
3.0515 ≈ 0.1779`, matching); `C_CIPA_STAB` tracks `C_CIPA * (1 +
EFFECT_MIG)` exactly at every sampled point.

`/model_manifest` on the refactored DSL confirms 49 output paths (22
compartments + 27 captured quantities, up from the original's 22
compartments with 0 actually-retrievable captured quantities, see defect 2
above) and 73 parameters, including every renamed/new PK and Hill-interface
parameter (`GAMMA_MIG`, `EMAX_MIG`, `EC50_MIG`, `V1_MIG`, and every
existing `KA_ALGLU`/`CL_ALGLU`/`V1_ALGLU`/`Q_ALGLU`/`V2_ALGLU`-style PK
parameter for all four compounds).

## Anything else worth flagging

- **AAV9-hGAA gene therapy (`AAV_X`/`AAV_DOSE`/`AAV_kexp`/`AAV_DECAY`/
  `AAV_GAIN`/`AAV_NAB`) and Rituximab-based ITI (`RTX_C`/`RTX_KIN`/
  `RTX_KOUT`/`RTX_ADA_K`/`RTX_EC50`) are completely untouched** — no
  census row exists for either, and neither compartment, parameter, or
  local variable belonging to them was renamed or touched. `Cp_rtx` (the
  local double `RTX_C / V1_ALGLU`, an approximation the original itself
  flags with `// approx`) is left with its original name for the same
  reason, and is still captured via `$CAPTURE` exactly as before.
- **`ag_drive` (antigen exposure driving anti-GAA ADA production) still
  sums all three ERTs' concentrations** (`(C_ALGLU + C_AVAL + C_CIPA) *
  ADA_AMP`, renamed inputs only). This is *not* a violation of "never
  collapse several drugs into one shared Hill term" — `ag_drive` is not a
  disease-effect Hill term at all, it is an immunogenicity-subsystem input
  (antigen exposure) feeding `ADA_T`'s own dynamics, structurally identical
  to the original and untouched beyond the necessary rename of its three
  inputs.
- **`ada_block` (shared ADA-neutralization multiplier) is applied with a
  different weight per compound** — full weight for `EFFECT_ALGLU`, half
  weight (`0.5*ada_block`) for `EFFECT_AVAL` and `EFFECT_CIPA` — exactly as
  the original had it. Each compound's `EFFECT_<STEM>` remains a distinct,
  independently-driveable named quantity; only the shared `ada_block` input
  (itself disease-state-derived, not drug-concentration-derived) is common
  across them, same as `ag_drive` above.
- **Miglustat's effect is combined with Cipaglucosidase alfa only at the
  point of actual use** (`C_CIPA_STAB = C_CIPA * (1 + EFFECT_MIG)`, feeding
  `EFFECT_CIPA`), not inside Cipaglucosidase alfa's own PK block —
  `C_CIPA` itself remains Cipaglucosidase alfa's own, undisturbed plasma
  concentration (`CENT_CIPA / V1_CIPA`), matching the guide's "combine them
  only at the point the disease equations actually use them" rule.
- Compartment ordering is fully preserved — no compartment was added,
  removed, or reordered (`CENT_ALGLU` is still compartment 1, `GUT_MIG` is
  still compartment 7, etc.), so any external code addressing this model's
  compartments by 1-based number rather than by name is unaffected.
- The R-side `pompe_run()` scenario runner's `cmt =` dosing targets were
  updated to the renamed compartment names (`"CENT_ALGLU"`, `"CENT_AVAL"`,
  `"CENT_CIPA"`, `"GUT_MIG"`); the `"aav_gt"` and `"alglu_iti"` scenarios'
  `"AAV_X"`/`"RTX_C"` targets are untouched (out of scope). All dose
  amounts, `q2w()` helper logic, and `param` overrides are otherwise
  identical to the original.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`pompe-disease | Alglucosidase alfa`, `pompe-disease | Avalglucosidase
alfa`, `pompe-disease | Cipaglucosidase alfa`, and `pompe-disease |
Miglustat` rows.
