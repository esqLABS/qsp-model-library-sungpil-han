# Refactor notes — `alzheimers-disease/ad_mrgsolve_model.R`

Three compounds in this file, all refactored independently per the fork
guide's classification (`driver-patches/data/compound_perturbation_census.md`,
all three rows tagged "Redirect concentration (clean single site)"):
**Donepezil (DON)**, **Lecanemab (LEC)**, and **Memantine (MEM)**. No other
compound exists in this file.

## Archetype — DON (donepezil)

**Archetype 3 (depot + central + peripheral, linear) plus one bespoke,
non-standard fourth "CNS tissue" compartment**, following the same
precedent as the extra tissue compartments in
`acute-intermittent-porphyria/aip_refactor_notes.md`. The original's
`DON_GUT → DON_CENT ⇄ DON_PERI` block is a textbook oral depot/central/
peripheral linear PK system, renamed straight to convention:
`GUT_DON`/`CENT_DON`/`PERI_DON`. The original central volume/peripheral
volume were named backwards from this guide's convention (`V2_DON` = the
*central* volume, `V3_DON` = the *peripheral* volume) — fully renamed to
`V1_DON`/`V2_DON` per the guide's table (central=`V1`, peripheral=`V2`),
not left as `V2`/`V3`.

Beyond the three-compartment linear system, the original adds a fourth
compartment, `DON_CNS` (renamed `CNS_DON`), governed by a first-order
relaxation ODE toward a fraction (`Kp_DON`≈0.18) of the *free* plasma
concentration: `dxdt_CNS_DON = 0.5*(KP_DON*CENT_DON/V1_DON - CNS_DON)`.
This is not a mammillary peripheral compartment (no back-flux into
`CENT_DON`) and not receptor-binding kinetics (archetype 4) — just a
pseudo-equilibrium brain-exposure surrogate, the same "extra role beyond
the guide's four-slot table" pattern as `LIV_GIV`/`LIV_HEM` in
`acute-intermittent-porphyria` and `SYS_VIT` in
`age-related-macular-degeneration`.

`CNS_DON` — not `CENT_DON` — is the *only* compartment the AChE-inhibition
term reads, so it is the exposed concentration: `C_DON = CNS_DON`. Like
`LIV_GIV`/`TOCI_C`, `CNS_DON` is already declared and integrated as a
concentration state (mg/L, via the relaxation ODE), not an amount divided
by a volume, so `C_DON = CNS_DON` is a straight alias, no division.

**Dead code preserved as-is (present in the original, not introduced by
this refactor):** `double CL2 = CL_DON; double V2 = V2_DON; double Q3 =
Q_DON; double V3 = V3_DON;` in the original `$MAIN` are computed and never
read again anywhere in `$MAIN`/`$ODE`/`$TABLE` — every actual PK equation
uses `CL_DON`/`Q_DON`/`V2_DON`/`V3_DON` directly. Kept in the refactored
file (renamed to reference the now-renamed parameters, `V2 = V1_DON;
V3 = V2_DON;`, otherwise they would reference undefined symbols) as inert,
functionally identical dead locals — this is a mechanical consequence of
renaming the parameters they alias, not a new edit to the original's own
(inert) logic.

## Archetype — MEM (memantine)

**Archetype 1 (no depot, single central compartment, linear elimination)
plus the same kind of bespoke second "CNS tissue" compartment as DON.**
The original declares `KA_MEM` (an absorption rate constant) in `$PARAM`,
but its own `$ODE` never uses it — `dxdt_MEM_CENT = -CL_MEM/V_MEM*MEM_CENT`
is the entire PK equation for the central compartment, with no
`dxdt_MEM_GUT`/depot compartment anywhere in `$CMT`, and the file's own
dosing (`mem_dose <- ev(amt=20, ..., cmt="MEM_CENT")`) doses `MEM_CENT`
directly. So this is archetype 1, not archetype 3 with a depot — `KA_MEM`
is a dead parameter in the original (same "declared but never read"
pattern as `V1_GIV`/`V2_GIV` in the acute-intermittent-porphyria
calibration run) and is carried through unchanged (still declared, still
unused) in the refactored file, per the "don't invent, don't drop the
original's parameters" rule.

`MEM_CENT` renamed to `CENT_MEM` (`V_MEM`→`V1_MEM` per the central-volume
convention). The bespoke `MEM_CNS`→`CNS_MEM` compartment follows the exact
same pseudo-equilibrium relaxation pattern as `CNS_DON`
(`dxdt_CNS_MEM = 0.3*(KP_MEM*CENT_MEM/V1_MEM - CNS_MEM)`), and is again the
only compartment the drug's disease-facing effect term reads, so
`C_MEM = CNS_MEM` (straight alias, same reasoning as `C_DON`).

## Archetype — LEC (lecanemab)

**Archetype 2 (no depot, two-compartment, linear — already written in the
original in CL/Q/V form, no micro-constant conversion needed) plus the same
kind of bespoke third "CNS tissue" compartment as DON/MEM.** `LEC_CENT ⇄
LEC_PERI` (`CL_LEC`, `V1_LEC`, `Q_LEC`, `V2_LEC`) is dosed by IV directly
into the central compartment — a completely ordinary, fully **linear**
2-compartment mAb PK system, renamed straight to convention
(`CENT_LEC`/`PERI_LEC`; `V1_LEC`/`V2_LEC` were already correctly assigned
to central/peripheral in the original, unlike DON's backwards `V2`/`V3`).

**Is this TMDD?** The task asked this explicitly because lecanemab is a
monoclonal antibody, where target-mediated drug disposition (nonlinear
clearance via a saturable receptor pool) is common. Checked directly: **no.
There is no free-receptor/drug-complex compartment pair anywhere in this
file** — nothing plays the `REC_FREE_TCZ`/`COMPLEX_TCZ` role from the
guide's archetype-4 template. What the original has instead is a third,
bespoke `LEC_CNS` (renamed `CNS_LEC`) compartment, fed by a fixed-rate BBB
influx (`KBBB_LEC*LEC_CENT`, no saturation) and cleared by a fixed
first-order rate (`KCNS_LEC*LEC_CNS`) *plus* one additional loss term,
`inh_proto*0.05*LEC_CNS` (renamed `EFFECT_LEC*0.05*CNS_LEC`), representing
protofibril-bound drug being consumed. `inh_proto`/`EFFECT_LEC` itself is a
plain algebraic saturation ratio (`Ccns_LEC/(KD_proto_LEC+Ccns_LEC)`), not
an ODE-integrated receptor-occupancy state — there is no dynamic free/bound
target pool that itself evolves and feeds back into a modeled synthesis/
degradation balance, which is what would make this genuinely archetype 4.
The self-referential loss term (the compound's own effect term feeding
back into its own CNS clearance) stays entirely inside LEC's own PK block —
it never reads a foreign disease-PD state — so it does not violate the "PK
lives in its own block" rule either. Net conclusion: **linear 2-compartment
plasma PK + a saturable-but-not-receptor-dynamic CNS loss term**, not TMDD;
kept as archetype 2 + bespoke CNS compartment, the same structural family
as DON and MEM rather than a fourth, different pattern.

`CNS_LEC` is the only compartment the protofibril-neutralization and
plaque-deposition disease terms read, so `C_LEC = CNS_LEC` (straight alias,
same reasoning as `C_DON`/`C_MEM`).

## The Hill interface — all three are renames, not fits

All three original effect terms are already exact Hill-shaped (or
Hill-with-implicit-gamma=1) ratios — no ODE-solved kinetics, receptor
occupancy, or curve-fitting involved for any of the three compounds:

- **DON**: `inh_AChE = Imax_DON * pow(Cb_DON,HILL_DON) /
  (pow(IC50_AChE_DON,HILL_DON)+pow(Cb_DON,HILL_DON))` with
  `Imax_DON = 1.0` hardcoded in the original's `$MAIN` (not previously a
  `$PARAM`). Renamed verbatim: `EMAX_DON = 1.0` (promoted to an explicit
  `$PARAM`, same value, same precedent as `EMAX_HEM`/`EMAX_GIV` in
  `acute-intermittent-porphyria/aip_refactor_notes.md`), `EC50_DON =
  IC50_AChE_DON = 6.7e-6`, `GAMMA_DON = HILL_DON = 1.2`. `EFFECT_DON` is a
  straight rename of `inh_AChE`.
- **MEM**: `inh_NMDA = Cb_MEM/(IC50_NMDA_MEM+Cb_MEM)` — no exponent applied
  in the original despite `HILL_MEM = 1.0` being declared in `$PARAM` (a
  dead parameter in the original, like `KA_MEM`). `EMAX_MEM = 1.0` added
  explicit (forced by the ratio's own asymptote as `C→∞`, same "algebra-
  forced, not invented" reasoning as `EMAX_HEM`). `EC50_MEM =
  IC50_NMDA_MEM = 0.8`, `GAMMA_MEM = HILL_MEM = 1.0`. Since `GAMMA_MEM=1`
  and `pow(x,1)=x` exactly, applying the exponent in the refactored
  `EFFECT_MEM` changes nothing numerically versus the original's
  exponent-free ratio — confirmed by the bit-exact verification below, not
  just asserted.
- **LEC**: `inh_proto = Ccns_LEC/(KD_proto_LEC+Ccns_LEC)` — same pattern as
  MEM, `HILL_LEC = 1.0` declared but never applied as an exponent in the
  original. `EMAX_LEC = 1.0` added explicit, `EC50_LEC = KD_proto_LEC =
  5.4e-5`, `GAMMA_LEC = HILL_LEC = 1.0`. Same `pow(x,1)=x` no-op reasoning
  as MEM.

Each `EFFECT_<STEM>` is used at exactly the same site(s) the original's
`inh_*` variable was used, nothing merged or combined:
`EFFECT_DON` → `kdeg_ACh_eff` (cholinergic pathway, one site); `EFFECT_MEM`
→ the synaptic-protection term in `dxdt_SYN` (one site); `EFFECT_LEC` →
three sites, all inside the amyloid cascade and LEC's own PK
(`dxdt_AB_PROTO`, `dxdt_AB_PLAQUE`, and `dxdt_CNS_LEC`'s own self-loss
term) — the same three sites `inh_proto` fed in the original, none of them
added or dropped.

## `C_<STEM>`/`EFFECT_<STEM>` are `$ODE`-local + `$CAPTURE` only, not `$PARAM`

Attempted putting all six (`C_DON`, `EFFECT_DON`, `C_MEM`, `EFFECT_MEM`,
`C_LEC`, `EFFECT_LEC`) in `$PARAM` per the qspserver-compatibility
section's literal request, and independently confirmed (not just cited)
the same mrgsolve 2.0.1 build incompatibility already documented in
`age-related-macular-degeneration/amd_refactor_notes.md` and
`acute-intermittent-porphyria/aip_refactor_notes.md`: adding `C_DON = 0` to
`$PARAM` while `$ODE` still does `double C_DON = CNS_DON;` fails to build
with

```
221:10: error: assignment of read-only reference 'C_DON'
  221 |   C_DON  = CNS_DON;
```

— a `$PARAM` member is passed into `$ODE` as a read-only reference in this
build, so a per-step-recomputed local and a settable `$PARAM` default are
mutually exclusive for the same name, exactly as found before. Given that
conflict, all six are declared as plain `$ODE`-local `double`s and exposed
via `$CAPTURE` under their bare names (no `_out` suffix needed — no name
collision exists in this file). Fully visible via `POST /run_simulation`'s
`outputs` selection (confirmed in `/model_manifest`'s `outputPaths`); **not**
visible via `/model_manifest`'s parameter listing (confirmed empirically —
none of the six appear under `parameters`). A future driver-patch redirect
of `C_DON`/`C_MEM`/`C_LEC` will need to replace the `double C_<STEM> =
CNS_<STEM>;` line itself, not override a `$PARAM` default.

## One pre-existing upstream defect found while verifying (not fixed in the checked-in original; logged as `translations/UPSTREAM_ISSUES.md` #45)

The untouched original does not compile under mrgsolve 2.0.1 at all — the
defect is unrelated to any of the three compounds' own PK, so per the
guide's settled policy it was fixed **directly in the delivered
`ad_mrgsolve_model_refactored.R`** (not just a scratch workaround), fully
disclosed here, and logged upstream:

**`$CAPTURE` repeats ten disease-state compartment names already declared
in `$CMT`** (`AB_MONO, AB_OLIGO, AB_PROTO, AB_PLAQUE, TAU_SOL, TAU_PHOS,
TAU_AGG, NEURO_INFLAM, ACH, SYN`):

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  AB_MONO,AB_OLIGO,AB_PROTO,AB_PLAQUE,TAU_SOL,TAU_PHOS,TAU_AGG,
  NEURO_INFLAM,ACH,SYN
```

**Fix applied, in the delivered `_refactored.R` only:** the ten names were
dropped from `$CAPTURE`. This changes nothing numerically or in what is
reported — mrgsolve compartments are already present in every simulation's
output regardless of `$CAPTURE` membership, `$CAPTURE` only controls
*additional* `$TABLE`-computed doubles. The verification run below
confirms all ten states are still fully present and bit-exact in the
refactored model's output under their own names. The untouched
`ad_mrgsolve_model.R` still carries the defect forward unfixed, per the
never-edit-upstream rule.

## Verification

**Method**: qspserver `mrgsolve_api` container, `POST /model_manifest` and
`POST /run_simulation`, `http://localhost:8007` — no local R/mrgsolve used,
per this fork's verification protocol. The bare DSL text was extracted from
both `code <- '...'` blocks (original + the build-fix described above) for
`model_content`; the extracted refactored text was confirmed byte-identical
to the `code <- '...'` block actually shipped in
`ad_mrgsolve_model_refactored.R`.

**Scenarios tested**: the original file's own named scenarios 2
("Donepezil_10mg"), 3 ("Memantine_20mg"), and 5 ("Lecanemab") from its own
`scenarios` list — one compound dosed at a time, avoiding the ambiguity of
`ev_seq()`'s sequential (not simultaneous) semantics in scenario 4's combo
regimen, which is irrelevant to verifying that each compound's own PK/PD
block still computes identically after renaming. Time grid matches the
original's own `simulate_scenario(duration_years = 3)` call exactly:
`end = 3*8760 = 26280`, `delta = 8` (3287 points); no shortening was
needed, every run completed well inside the API's default step budget.

- **Donepezil**: `amt=10, cmt=GUT_DON(1), ii=24, addl=729` (matches
  `don_dose`).
- **Memantine**: `amt=20, cmt=CENT_MEM(5), ii=24, addl=729` (matches
  `mem_dose`).
- **Lecanemab**: `amt=700, cmt=CENT_LEC(7), ii=336, addl=25`. The
  original's `rate = -2` (mrgsolve's "duration-parameter-driven infusion"
  code) has no matching `D_LEC_CENT`/`D_CENT_LEC` parameter defined
  anywhere in the file — confirmed by trying it literally through the API,
  which fails identically for both models with `modeled infusion duration
  D_CMT or Dn must be positive when dosing record RATE is set to -2`.
  Translated faithfully to the comment's own stated intent ("infusion over
  1h") with a literal zero-order rate, `rate = amt / 1h = 700`, applied
  identically to both models — same approach as
  `acute-intermittent-porphyria/aip_refactor_notes.md`'s hemin dosing.

For each scenario, every original compartment/`$CAPTURE` output was
compared point-by-point across the full 3287-point time grid against its
renamed counterpart in the refactored model:

- **Donepezil** (`DON_GUT/DON_CENT/DON_PERI/DON_CNS` ↔
  `GUT_DON/CENT_DON/PERI_DON/CNS_DON`, plus `Cp_Donepezil`,
  `Ccns_Donepezil`, `AChE_inhibition`, `MMSE`, `ADAS_Cog`, `CDR_SB`, `ACH`,
  `SYN`): **max abs diff = 0.0, max rel diff = 0.0** for every output at
  every non-`NaN` timepoint (see the NaN note below).
- **Memantine** (`MEM_CENT/MEM_CNS` ↔ `CENT_MEM/CNS_MEM`, plus
  `Cp_Memantine`, `NMDAR_occupancy`, `MMSE`, `ADAS_Cog`, `CDR_SB`, `SYN`):
  **max abs diff = 0.0, max rel diff = 0.0** across all 3287 points, no
  `NaN`s in this scenario.
- **Lecanemab** (`LEC_CENT/LEC_PERI/LEC_CNS` ↔ `CENT_LEC/PERI_LEC/CNS_LEC`,
  plus `Cp_Lecanemab`, `Proto_neutralized`, `AB_PROTO`, `AB_PLAQUE`, `MMSE`,
  `AmyloidPET_CL`, `CSF_pTau181`): **max abs diff = 0.0, max rel diff =
  0.0** across all 3287 points, no `NaN`s in this scenario.

This is a pure structural reorganization for all three compounds (every
effect term is an exact rename, `GAMMA_MEM`/`GAMMA_LEC=1` applied via
`pow(x,1)=x` with no numeric effect), so the guide's "near-exact match
expected" standard for archetypes 1–3 is satisfied — in practice the match
came out bit-for-bit exact, same as every prior rename-only refactor in
this fork.

**Pre-existing `NaN` blowup, identical in both models, not introduced by
this refactor:** in the Donepezil scenario only, `ACH` (and every
downstream state/output depending on it — `SYN`, `MMSE`, `COGNITION`, ...)
turns `NaN` at `time = 23528` h (~2.7 years into the 3-year run) with no
visibly divergent value beforehand, and stays `NaN` through the rest of the
horizon. Confirmed at the identical time index with byte-identical
pre-blowup values in both the original and the refactored model — a
pre-existing long-horizon solver artifact (see `UPSTREAM_ISSUES.md` #45 for
detail), not a difference introduced by the rename. All 333 `NaN` entries
occur at the same indices in both models for every affected output; every
one of the preceding 2954 points matches exactly (0.0/0.0). The Memantine
and Lecanemab scenarios, run to the same 3-year horizon, show no such
blowup.

**`$PARAM`-vs-`$CAPTURE` confirmation:** `/model_manifest` was called on
both models, and every parameter name/value pair was diffed
programmatically (not eyeballed). Original: 66 `$PARAM`s; refactored: 69
(the three added `EMAX_*`, since `GAMMA_*`/`EC50_*` replace `HILL_*`/
`IC50_*`/`KD_proto_*` 1-for-1 rather than adding new names). 55 names are
shared between the two files with byte-identical values (every disease
parameter — amyloid/tau/neuroinflammation/cholinergic/synaptic/cognitive/
disease-modifier — plus the compound parameters whose old name and new
name happen to coincide, e.g. `CL_DON`, `Q_LEC`). **One coincidental name
collision, disclosed for rigor:** `V2_DON` exists in *both* files but means
something different in each — in the original it is the *central* volume
(594 L); in the refactored file the guide's naming convention reassigns
`V2_<STEM>` to the *peripheral* role, so `V2_DON` there is the renamed
`V3_DON` (1180 L). This is expected (the guide's central/peripheral slot
numbering does not match this compound's original numbering) and does not
indicate any value carried over incorrectly — the diff script caught it
by comparing values, not just name presence, and every equation that reads
`V2_DON` in the refactored `$ODE`/`$MAIN`/`$TABLE` uses it exclusively in
the peripheral-volume role (see `$ODE` above), never mixed with the
central role. Every other renamed PK/Hill parameter for all three
compounds (`KA_DON`, `V1_DON`, `EMAX_DON`, `EC50_DON`, `GAMMA_DON`,
`V1_MEM`, `EMAX_MEM`, `EC50_MEM`, `GAMMA_MEM`, `V1_LEC`, `V2_LEC`,
`EMAX_LEC`, `EC50_LEC`, `GAMMA_LEC`, etc.) is listed with the original's
own values under its new name, none invented; confirmed `C_DON`/
`EFFECT_DON`/`C_MEM`/`EFFECT_MEM`/`C_LEC`/`EFFECT_LEC` do **not** appear in
the manifest's parameter list but do appear in `outputPaths`, per the
`$PARAM`-vs-`$CAPTURE` finding above.

**Diagnostic sanity check (no original counterpart to diff against):**
`C_DON` ranges 0 – 0.0057 mg/L with `EFFECT_DON` reaching ~99.97% AChE
inhibition at steady dosing (consistent with a sub-µg/L `EC50_DON` and
sustained daily 10 mg dosing); `C_MEM` ranges 0 – 0.269 mg/L with
`EFFECT_MEM` up to ~25% NMDAR occupancy (consistent with `EC50_MEM = 0.8`
mg/L against typical memantine steady-state trough exposure); `C_LEC`
ranges 0 – 1.03 mg/L with `EFFECT_LEC` reaching ~99.99% protofibril
neutralization during the biweekly-infusion maintenance phase (consistent
with a sub-µg/L `EC50_LEC` and repeated 700 mg IV dosing) — all three
ranges are physiologically sensible given the model's own parameters.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`alzheimers-disease | DON`, `alzheimers-disease | LEC`, and
`alzheimers-disease | MEM` rows: target compartments `CNS_DON`/`CNS_MEM`/
`CNS_LEC` respectively (all three CNS tissue compartments, not the plasma
central compartments); exposed variables `C_DON`/`EFFECT_DON` (mg/L;
fraction 0–1, AChE inhibition), `C_MEM`/`EFFECT_MEM` (mg/L; fraction 0–1,
NMDAR occupancy), `C_LEC`/`EFFECT_LEC` (mg/L; fraction 0–1, protofibril
neutralization).
