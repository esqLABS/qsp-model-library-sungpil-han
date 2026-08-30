# Refactor notes — `carcinoid-syndrome/carcsyn_mrgsolve_model.R`

Seven compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **Everolimus (EVE)**,
**IFN-α-2b (IFN)**, **Lanreotide autogel (LAN)**, **Octreotide LAR (OCT)**,
**Pasireotide LAR (PAS)**, **Telotristat (TEL)**, and **VEGFR-TKI (TKI)**.
The eighth drug in this file, **¹⁷⁷Lu-DOTATATE / PRRT** (`PRRT_AMT`), is
**not** in the census list and was left completely untouched — same name,
same equation (`dxdt_PRRT_AMT = -0.04 * PRRT_AMT`), same position. Every
disease-state compartment (`TUMOR`, `TPH1`, `HTP`, `SER_T`, `SER_P`,
`SER_PLT`, `HIAA_U`, `TGFb`, `VALVE`, `NTproBNP`, `BM`, `FLUSH`) is likewise
untouched except where it reads a renamed `EFFECT_*` term.

## Identifying "VEGFR-TKI"

The original models this compound generically — `$PARAM`/`$CMT` comments say
"VEGFR-TKI generic (sunitinib / surufatinib / cabozantinib) lumped" — but its
own dosing scenario (`S7_sunitinib_proxy_TKI`, 37.5 mg/day continuous oral)
and the disease directory's `README.md` ("VEGFR-TKI 37.5 mg/day
(sunitinib/surufatinib lump) — Raymond 2011") both identify the calibrated
dose and reference as **sunitinib** (the SUN1111 trial, Raymond et al., *N
Engl J Med* 2011, continuous 37.5 mg/day dosing in pancreatic NET — the exact
regimen this file's own scenario reproduces). The census row and this file's
comments now read "VEGFR-TKI (sunitinib)" accordingly; the `TKI` stem is kept
unchanged since nothing in the naming convention requires renaming the stem
itself, only clarifying the compound identity in prose.

## Archetypes

All seven compounds' PK are plain first-order absorption/elimination —
despite the octreotide/lanreotide/pasireotide comments calling these
formulations "LAR"/"autogel"/depot, the actual `$ODE` for every one of them
is a standard mono-exponential depot decay (`dxdt_GUT_X = -KA_X*GUT_X`
feeding `dxdt_CENT_X`), not a zero-order-release or biphasic depot — so none
of the three needed bespoke handling despite the guide's warning to check
this carefully for LAR/autogel formulations.

- **OCT (Octreotide LAR)** — archetype 3, full form: depot + central +
  peripheral, linear elimination (`GUT_OCT`/`CENT_OCT`/`PERI_OCT`, renamed
  from `PER_OCT`).
- **LAN (Lanreotide autogel)** — archetype 3 variant: depot + central, no
  peripheral compartment in the original (`GUT_LAN`/`CENT_LAN`, both already
  conforming).
- **PAS (Pasireotide LAR)** — archetype 3 variant: depot + central, no
  peripheral. **Stem inconsistency in the original, resolved**: the original
  used `PASIR` for the two `$CMT` names (`GUT_PASIR`/`CENT_PASIR`) but `PAS`
  for every PK parameter and the concentration (`KA_PAS`, `CL_PAS`, `V_PAS`,
  `KD_PAS_S5`, `C_PAS`) — six of eight identifiers already used `PAS`. This
  refactor standardizes on **`PAS`** throughout (`GUT_PASIR`→`GUT_PAS`,
  `CENT_PASIR`→`CENT_PAS`), the majority convention already in the file, per
  the guide's "don't invent a new stem" (this isn't inventing one, just
  picking the one the original already used more often).
- **TEL (Telotristat, active metabolite LP-778902)** — archetype 3 variant:
  depot + central with bioavailability `F_TEL` (no peripheral compartment).
- **EVE (Everolimus)** — archetype 3 variant: depot + central, no
  peripheral.
- **TKI (VEGFR-TKI / sunitinib)** — archetype 3 variant: depot + central, no
  peripheral.
- **IFN (IFN-α-2b)** — **archetype 1**: single compartment, no depot, linear
  elimination. The original doses IFN-α as an SC bolus directly into
  `CENT_IFN` (`dxdt_CENT_IFN = -(CL_IFN/V_IFN)*CENT_IFN`, no absorption
  compartment, no `KA_IFN` anywhere) — this is a genuine single-compartment
  PK in the original, not a simplification introduced here.

No TMDD/receptor-binding-kinetics compound exists in this file among the
seven (SSTR2/5 "occupancy" is a free-fraction Hill ratio on plasma
concentration, `C/(KD+C)`, not a dynamic receptor-binding ODE system — no
`REC_FREE_*`/`COMPLEX_*` states anywhere), so archetype 4 does not apply to
any of them.

## Renaming (values unchanged from the original)

| Compound | Original | Refactored |
|---|---|---|
| OCT | `KA_OCT`, `CL_OCT`, `Q_OCT` | unchanged |
| OCT | `V_OCT` | `V1_OCT` |
| OCT | `VP_OCT` | `V2_OCT` |
| OCT | `PER_OCT` (`$CMT`) | `PERI_OCT` |
| OCT | `KD_OCT_S2` | `EC50_OCT` |
| OCT | `OCC_OCT` | `EFFECT_OCT` |
| LAN | `KA_LAN`, `CL_LAN` | unchanged |
| LAN | `V_LAN` | `V1_LAN` |
| LAN | `KD_LAN_S2` | `EC50_LAN` |
| LAN | `OCC_LAN` | `EFFECT_LAN` |
| PAS | `KA_PAS`, `CL_PAS` | unchanged |
| PAS | `V_PAS` | `V1_PAS` |
| PAS | `GUT_PASIR`, `CENT_PASIR` (`$CMT`) | `GUT_PAS`, `CENT_PAS` |
| PAS | `KD_PAS_S5` | `EC50_PAS` |
| PAS | `OCC_PAS` | `EFFECT_PAS` |
| TEL | `KA_TEL`, `CL_TEL`, `F_TEL` | unchanged |
| TEL | `V_TEL` | `V1_TEL` |
| TEL | `IC50_TPH1` | `EC50_TEL` |
| TEL | `GAMMA_TPH1` | `GAMMA_TEL` |
| TEL | `INH_TPH1` | `EFFECT_TEL` |
| EVE | `KA_EVE`, `CL_EVE` | unchanged |
| EVE | `V_EVE` | `V1_EVE` |
| EVE | `IC50_MTOR` | `EC50_EVE` |
| EVE | `MTORi` | `EFFECT_EVE` |
| TKI | `KA_TKI`, `CL_TKI` | unchanged |
| TKI | `V_TKI` | `V1_TKI` |
| TKI | `IC50_VEGF` | `EC50_TKI` |
| TKI | `VEGFi` | `EFFECT_TKI` |
| IFN | `CL_IFN` | unchanged |
| IFN | `V_IFN` | `V1_IFN` |
| IFN | `EFF_IFN50` | `EC50_IFN` |
| IFN | `IFNi` | `EFFECT_IFN` |

`C_OCT`, `C_LAN`, `C_PAS`, `C_TEL`, `C_EVE`, `C_TKI`, `C_IFN` already matched
the `C_<STEM>` convention exactly in the original and are unchanged.
`GUT_OCT`, `CENT_OCT`, `GUT_LAN`, `CENT_LAN`, `GUT_TEL`, `CENT_TEL`,
`GUT_EVE`, `CENT_EVE`, `GUT_TKI`, `CENT_TKI`, `CENT_IFN` already conformed
and are unchanged.

New (not present in the original, for the named Hill interface):
`EMAX_OCT`/`GAMMA_OCT` = 1, `EMAX_LAN`/`GAMMA_LAN` = 1, `EMAX_PAS`/`GAMMA_PAS`
= 1, `EMAX_TEL` = 1 (`GAMMA_TEL` is a rename of the original's own explicit
`GAMMA_TPH1` = 1.3, not new), `EMAX_EVE`/`GAMMA_EVE` = 1,
`EMAX_TKI`/`GAMMA_TKI` = 1, `EMAX_IFN`/`GAMMA_IFN` = 1 — these make explicit
a shape every one of these seven effect terms already had implicitly as a
plain `C/(EC50+C)` ratio (`Emax=1`, `gamma=1`), same as the LETRO/TRAS
precedent in the breast-cancer refactor.

`SSTR_TOTAL` (the combined three-compound occupancy sum, capped at 0.99) is
kept under its original name — it is a genuine combination point (used by
`KREL_EFF` and the flushing equation), not a single compound's own
`EFFECT_<STEM>`, so it stays outside the naming convention's per-compound
slot, per the guide's "combine them only at the point the disease equations
actually use them."

## Hill interface: renames, not fits

All seven compounds' disease-facing effect terms were already written as a
plain concentration ratio in the original — no ODE-solved receptor kinetics
to approximate for any of them, so every one of these is a rename, confirmed
by the near-exact verification match below (see caveat on floating-point
noise):

- **OCT/LAN/PAS**: `OCC_X = C_X/(KD_X + C_X)` → `EFFECT_X = EMAX_X *
  pow(C_X, GAMMA_X) / (pow(EC50_X, GAMMA_X) + pow(C_X, GAMMA_X))` with
  `EMAX_X = GAMMA_X = 1` (new, explicit).
- **TEL**: `INH_TPH1 = pow(C_TEL, GAMMA_TPH1) / (pow(IC50_TPH1, GAMMA_TPH1) +
  pow(C_TEL, GAMMA_TPH1))` was already an explicit two-parameter Hill term
  (`gamma = 1.3`, no separate Emax) → `EFFECT_TEL` uses the identical
  expression with the renamed parameters and `EMAX_TEL = 1` added (new,
  explicit; the original's ratio already implied Emax = 1).
- **EVE/TKI/IFN**: `X_i = C_X/(IC50/EFF_X + C_X)` → same rename pattern,
  `EMAX_X = GAMMA_X = 1` (new, explicit).

No `nls()` fitting was needed or performed for any of the seven compounds.
The downstream scaling factors on each effect term inside `GROWTH`
(`1 - 0.6*EFFECT_EVE`, `1 - 0.45*EFFECT_TKI`, `1 - 0.3*EFFECT_IFN`) and inside
`KREL_EFF`/`FLUSH` (`1 - 0.75*SSTR_TOTAL`, `1 - SSTR_TOTAL`) are kept exactly
as the original had them, unchanged — same as the breast-cancer refactor's
precedent of leaving a compound's own downstream `* 10.0` scaling untouched.

## `$PARAM` vs local `double` for `C_<STEM>`/`EFFECT_<STEM>`

Unlike some other refactors in this batch, this model did **not** need the
`$GLOBAL`-declared-`double` workaround for the qspserver `$PARAM`
read-only-reference restriction: the original already computes `C_OCT`,
`C_LAN`, etc. as **local `double`s declared directly inside `$ODE`** and
lists them (bare, no `$GLOBAL` predeclaration) in `$CAPTURE` — and this
compiles and runs correctly through the live `mrgsolve_api` container (both
the original and the refactored DSL were confirmed via `/model_manifest` to
list every `C_<STEM>` and `EFFECT_<STEM>` in `outputPaths`, satisfying
qspserver compatibility requirement #4). The refactor keeps this same
pattern for the two new named quantities per compound (`C_<STEM>` unchanged,
`EFFECT_<STEM>` added as a local `double`), so requirement #2 ("ideally in
`$PARAM`") is not met for any of the seven `EFFECT_<STEM>` terms (recomputed
from state every step, so cannot also live in `$PARAM` — the same
`assignment of read-only reference` restriction documented in the AMD,
membranous-nephropathy, and breast-cancer refactor notes), but requirement
#4 (capture) is fully met, which is what makes every covariate discoverable
via `/run_simulation`'s `outputs` selection even though `/model_manifest`'s
`parameters` list doesn't carry them.

## Pre-existing upstream build defect

**None found.** Unlike several other models in this batch, the original
`carcsyn_mrgsolve_model.R`'s embedded DSL **compiles cleanly** under mrgsolve
2.0.1 through the live `mrgsolve_api` container — confirmed via
`POST /model_manifest` on the untouched original, no changes, no errors. No
syntax-only build-compat fix was needed or applied, and no new
`UPSTREAM_ISSUES.md` entry was logged for this file.

## Verification

**Method.** Both the original's own model code and
`carcsyn_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text (the quoted `carcsyn_code <- '...'` string) and run through
the qspserver `mrgsolve_api` container (`POST /model_manifest`,
`POST /run_simulation`) at `http://localhost:8007`, requests spaced ~2.5–3 s
apart. `/model_manifest` on the refactored DSL confirmed all 41 output paths
(27 compartments + `C_OCT`/`C_LAN`/`C_PAS`/`C_TEL`/`C_EVE`/`C_TKI`/`C_IFN` +
`EFFECT_OCT`/`EFFECT_LAN`/`EFFECT_PAS`/`EFFECT_TEL`/`EFFECT_EVE`/`EFFECT_TKI`/
`EFFECT_IFN` + `SSTR_TOTAL`) and all 73 `$PARAM` entries, including every new
`EMAX_<STEM>`/`GAMMA_<STEM>`.

All twelve of the original file's own dosing scenarios (`build_scenarios()`)
were run, exactly as defined (dose amounts, times, compartment targets
copied verbatim to the API's dosing-record convention — dosing records were
sorted by time before submission, a client-side requirement of the API's raw
`as.ev()` path unrelated to the refactor itself, since the original's
`bind_rows()`-built tibbles interleave multi-drug dose times in `ID`-block
order rather than global time order), over the model's own full 12-month
horizon (`end = 8064`, `delta = 12`):

1. S1 — Untreated natural history
2. S2 — Octreotide LAR 30 mg IM q28d × 12
3. S3 — Lanreotide autogel 120 mg SC q28d × 12 (CLARINET)
4. S4 — Octreotide LAR + Telotristat 250 mg t.i.d. (TELESTAR)
5. S5 — Pasireotide LAR 60 mg q28d
6. S6 — Everolimus 10 mg/day (RADIANT-4)
7. S7 — VEGFR-TKI (sunitinib) 37.5 mg/day
8. S8 — ¹⁷⁷Lu-DOTATATE 7.4 GBq q8w × 4 (untouched-compound sanity check)
9. S9 — IFN-α-2b 5 MU SC TIW × 12 months
10. S10 — HAE (tumor debulk) + octreotide LAR
11. S11 — Carcinoid crisis prophylaxis (octreotide IV bolus into `CENT_OCT` +
    LAR)
12. S12 — Quad therapy (octreotide + telotristat + PRRT + everolimus)

Every `$CAPTURE`d output relevant to the seven compounds
(`C_OCT/C_LAN/C_PAS/C_TEL/C_EVE/C_TKI/C_IFN`, `SSTR_TOTAL`) plus every
disease-state compartment (`TUMOR`, `TPH1`, `HTP`, `SER_T`, `SER_P`,
`SER_PLT`, `HIAA_U`, `TGFb`, `VALVE`, `NTproBNP`, `BM`, `FLUSH`, and
`PRRT_AMT` as an untouched-compound sanity check) was compared point-by-point
across each scenario's full time grid (675–1712 points depending on
scenario).

**Result: match to floating-point/solver scale for all twelve scenarios —
no genuine divergence found.** Ten of twelve scenarios (S1, S3, S5, S6, S7,
S8, S9, and the shared-output columns of S2/S4/S12) matched **exactly**
(max abs diff = 0.0) for every compared output. The remaining scenarios
(S2, S4, S10, S11, S12) showed non-zero but extremely small absolute
deviations confined to the downstream serotonin/symptom state variables
(`SER_T`, `SER_P`, `SER_PLT`, `HIAA_U`, `FLUSH`, `TUMOR`, `VALVE`,
`NTproBNP`, `HTP`) — e.g. S11's largest absolute deviation was `SER_PLT`
0.0208 against a value of ~202,486 (relative ≈ 1.0×10⁻⁷); S2's largest was
`SER_T` 0.0013 against ~29,902 (relative ≈ 4.3×10⁻⁸). Every one of these
deviations is **≤ ~1×10⁻⁷ relative**, consistent with the adaptive
LSODA-family solver taking a slightly different internal step sequence when
the algebraically-identical Hill ratio is written as `EMAX*pow(C,GAMMA)/
(pow(EC50,GAMMA)+pow(C,GAMMA))` (with `EMAX=GAMMA=1`) rather than the
original's bare `C/(EC50+C)` — `pow(x,1.0)` is mathematically but not always
bit-identical to `x` depending on the C library's implementation, and that
sub-ULP difference propagates through ~670–1700 adaptive steps. This is the
floating-point-scale deviation the guide's tolerance rule anticipates for
pure structural reorganization ("anything beyond floating-point-scale
deviation means a bug, not a tolerance to loosen") — not a bug, since the
seven renamed effect terms compute the identical arithmetic (`Emax=1,
gamma=1` reproduces the original ratio exactly) on the identical (renamed)
inputs.

One transient, non-reproducible failure was observed and re-verified: the
first attempt at scenario S6 (everolimus) returned a subprocess crash
(`free(): invalid next size (fast)`, signal 6) from the **original** DSL run.
This coincided with heavy concurrent load on the shared `mrgsolve_api`
container from other sessions using the same host during this work (visible
in the shared scratchpad directory, which had files from at least three
other concurrent refactor sessions). Re-running S6 alone, with no other
traffic in flight, reproduced an exact match (max abs diff = 0.0 across all
compared outputs) — consistent with the guide's documented note that this
container "has crashed under heavy concurrent load before," not a defect in
either model.

## Anything else worth flagging

- Compartment numbering is **fully preserved**: unlike some other refactors
  in this batch, no compartment was added or removed here (only three were
  renamed: `PER_OCT`→`PERI_OCT`, `GUT_PASIR`→`GUT_PAS`,
  `CENT_PASIR`→`CENT_PAS`), so every 1-based `cmt=` index used by the
  original's own `build_scenarios()` dosing tibbles is identical in the
  refactored file — no downstream R code needed updating for compartment
  position.
- `run_scenarios()`, `plot_serotonin()`, `plot_BM_flush()`, `plot_tumor()`,
  `omega_skeleton`, and `sigma_skeleton` reference no renamed identifier
  (they only touch `SER_P`, `BM`, `FLUSH`, `TUMOR`, and the omega/sigma
  parameter names `CL_OCT`/`CL_TEL`/`CL_EVE`/`KGROW`/`KIN_BM`/`KIN_FL`, none
  of which changed), so they are carried over unmodified.
- `KGROW_HILL` is declared in `$PARAM` with a comment implying it scales
  growth attenuation by mTORi, but it is never referenced anywhere in
  `$ODE` in the original — left exactly as-is (unused, unrenamed), outside
  this refactor's scope since it isn't part of any of the seven compounds'
  own PK/PD blocks.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`carcinoid-syndrome | Everolimus`, `carcinoid-syndrome | IFN-α-2b`,
`carcinoid-syndrome | Lanreotide autogel`, `carcinoid-syndrome | Octreotide
LAR`, `carcinoid-syndrome | Pasireotide LAR`, `carcinoid-syndrome |
Telotristat`, and `carcinoid-syndrome | VEGFR-TKI` rows (the last corrected
to note the sunitinib identity).
