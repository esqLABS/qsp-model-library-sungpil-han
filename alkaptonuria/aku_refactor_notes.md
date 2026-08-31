# Refactor notes — `aku_mrgsolve_model_refactored.R`

Scope: **nitisinone (NT) only** — per the `alkaptonuria | Nitisinone (NT)` row in
`driver-patches/data/compound_perturbation_census.md`, classified "Normalize
duplicate concentration sites, then redirect." Nitisinone is the file's only
drug (the disease's own metabolic pathway — Phe/Tyr/HPP/HPLA/HGA — is not a
"compound" in the census sense); every non-NT line is byte-for-byte identical
to the original (confirmed by `diff`, see below).

## What the duplicate-definition problem actually was

Nitisinone's plasma concentration (`NTCEN/V1`) was computed **four separate
times**, under **three different names**, with no single site downstream code
could be pointed at:

1. `$ODE`, cached once: `double CNT = NTCEN/V1; if(CNT<0.0) CNT=0.0;` — used
   by `KMAPP` (the HPD competitive-inhibition term) and by `dxdt_CUMNT`.
2. `$ODE`, `dxdt_NTCEN`'s own elimination/distribution line re-derived the
   *identical* expression **twice more, inline, uncached**:
   `-(CLNT*CYP3A4/V1)*NTCEN - Q*(NTCEN/V1 - NTPER/V2)` — two more evaluations
   of `NTCEN/V1` that never touched the `CNT` variable computed 30-odd lines
   above them.
3. `$TABLE`, a fourth, independent re-derivation under a fourth name:
   `double CNTo = NTCEN/V1; if(CNTo<0.0) CNTo=0.0;`, feeding `KMAPPo` →
   `JHPDo` → the captured trial endpoint `UHGA24`.

Four sites, three spellings (`CNT`, the two inline occurrences, `CNTo`), for
one physical quantity — exactly the naming/consolidation chaos the design
guide describes for tocilizumab across the corpus, reproduced *within a single
file* here.

## Redirect performed

**Within `$ODE`** (a single compiled scope, so true consolidation is possible
and was done): one cached `double C_NT = CENT_NT/V1_NT; if(C_NT<0.0) C_NT=0.0;`
is now the *only* site that derives the concentration; every downstream read
(`EFFECT_NT`/`KMAPP`, `dxdt_CENT_NT`, `dxdt_PERI_NT`, `dxdt_CUMNT`) redirects
to it. The two previously-inline, uncached re-derivations inside
`dxdt_NTCEN`/`dxdt_PERI_NT` are gone.

**Between `$ODE` and `$TABLE`:** these compile into the same translation unit
in mrgsolve — confirmed empirically against the qspserver `mrgsolve_api`
container, which rejected a first attempt at literally reusing the name
`C_NT` in both blocks with a hard compiler error:
```
93:10: error: redefinition of 'double {anonymous}::C_NT'
28:10: note: 'double {anonymous}::C_NT' previously declared here
```
So a single C++ variable cannot span both blocks; `$TABLE` must re-derive the
concentration itself. This is the same constraint already documented for this
exact "normalize duplicate concentration sites" classification in
`kidney-transplant-rejection/ktx_refactor_notes.md` (`CTCZr` kept distinct
from `C_TCZ` for the identical reason). Followed that precedent: kept
`$TABLE`'s local name as `CNTo` (was already distinct, so no rename needed
there) but **redirected its right-hand side** to the renamed, canonical
compartments/parameters (`CENT_NT/V1_NT`, was `NTCEN/V1`) instead of the
old names — the duplication across the block boundary is an unavoidable
mrgsolve architecture fact, not something left "inconsistent" any more.

## Naming convention applied (archetype 3: depot + central + peripheral, linear)

| Original | Refactored | Original | Refactored |
|---|---|---|---|
| `NTGUT` (cmt) | `GUT_NT` | `KA` | `KA_NT` |
| `NTCEN` (cmt) | `CENT_NT` | `V1` | `V1_NT` |
| `NTPER` (cmt) | `PERI_NT` | `V2` | `V2_NT` |
| `CLNT` | `CL_NT` | `Q` | `Q_NT` |
| `FBIO` | `F_NT` | `CYP3A4` | `CYP3A4_NT` |
| `CNT` (`$ODE` local) | `C_NT` | `HNT` | `GAMMA_NT` |
| `CNTo` (`$TABLE` local) | `CNTo` (kept, see above) | | |

All PK is a pure rename — no value changed, no compartment added or removed,
no algebra changed beyond substituting the cached `C_NT` for repeated inline
`NTCEN/V1`/`NTPER/V2` expressions (algebraically identical: `Q_NT*(C_NT -
PERI_NT/V2_NT)` is the same rearrangement of terms the original already used,
just reading the cached, clamped concentration instead of a raw inline
division). Bioavailability stays on mrgsolve's native `F_<CMT>` reserved-name
mechanism (`F_GUT_NT = F_NT;` in `$MAIN`, was `F_NTGUT = FBIO;`) rather than
being rewritten into the design guide's inline `KA*F*GUT` form — the two are
algebraically equivalent for the amount reaching `CENT_NT` (verified: for a
constant `F`, `GUT(t) = F·Dose·e^(-KA t)` either way), and the native
mechanism is what the original already used; switching it would change what
the `GUT` compartment itself reports without changing pharmacology, so it was
preserved as-is.

`HNT` → `GAMMA_NT`: this genuinely is the guide's "Hill coefficient" slot
(the original's own comment already called it "Hill exponent on HPD
inhibition"), so it was renamed for convention compliance. `KI_NT` was left
as `KI_NT`, not forced to `EC50_NT` — see below.

## The Hill interface: kept bespoke, not forced into Emax/EC50

Nitisinone's effect on the disease is **competitive inhibition of HPD**,
raising the enzyme's apparent Km:
`KMAPP = KMHPD*(1 + (C_NT/KI_NT)^GAMMA_NT)`. This is named as a single
variable, `EFFECT_NT = pow(C_NT/KI_NT, GAMMA_NT)`, used everywhere `KMAPP` is
computed (`$ODE` and, as `EFFECT_NTo` in `$TABLE`, for the identical
block-scoping reason above). It is **deliberately not** rewritten into the
guide's saturating `EMAX*X^GAMMA/(EC50^GAMMA+X^GAMMA)` template: that template
is bounded (approaches `EMAX` as `X→∞`), while competitive inhibition of an
enzyme's Km is structurally unbounded — `EFFECT_NT` keeps growing as `C_NT`
rises, which is exactly the mechanism the model's own header comment
describes at length ("residual flux ~ inversely proportional to dose"; "6-fold
difference in residual for a 5-fold dose ratio"). Forcing this into an
Emax/EC50 ratio would silently cap a mechanism the original leaves uncapped —
this is a rename, not a refit, and no curve-fitting was performed anywhere in
this file (`GAMMA_NT` is the original's own fitted `HNT = 1.45`, unchanged).

## Diff scope confirmation

`diff aku_mrgsolve_model.R aku_mrgsolve_model_refactored.R` touches only:
the `KA`/`V1`/`V2`/`Q`/`CLNT`/`FBIO`/`CYP3A4` param lines, the `KI_NT`/`HNT`
param lines, the `NTGUT`/`NTCEN`/`NTPER` `$CMT` lines, the `F_NTGUT` line in
`$MAIN`, the `CNT` definition and its three downstream reads in `$ODE`
(`KMAPP`/`EFFECT_NT`, `dxdt_NTGUT`/`dxdt_NTCEN`/`dxdt_NTPER`, `dxdt_CUMNT`),
the `CNTo`/`KMAPPo` lines and the new `CNTo`/`EFFECT_NTo` `$CAPTURE` entry in
`$TABLE`, the `aku_init()` compartment-name arguments, and the
`calibrate_aku()` fitted-parameter list's `HNT`→`GAMMA_NT` key. Every other
line — the entire Phe/Tyr/HPP/HPLA/HGA pathway, ochronosis, cartilage/spine/
valve/stone/pain/AKUSSI machinery, the 24-scenario list, and all calibration/
validation/analysis functions — is untouched.

## qspserver compatibility

Confirmed via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`): **the pristine, unmodified original
compiles as-is** (no build-compat fix was needed for this file — the "when
the original doesn't compile at all" section of the design guide does not
apply here). The refactored DSL also compiles; its manifest confirms
`KA_NT, V1_NT, V2_NT, Q_NT, CL_NT, F_NT, CYP3A4_NT, KI_NT, GAMMA_NT` are all
discoverable `$PARAM` entries (none of the old names remain), and its
`compartments` list shows `GUT_NT, CENT_NT, PERI_NT` in the same 1-based
positions the original's `NTGUT, NTCEN, NTPER` held (dosing into compartment
1 is therefore unaffected). `C_NT` itself could not be added to `$PARAM` (it
is computed, not overridable — assigning to a `$PARAM`-declared name is a
read-only-reference build error, the same constraint documented in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`); its `$ODE`-scope name is
`C_NT` and its `$TABLE`/`$CAPTURE`-scope name is the necessarily-distinct
`CNTo` (both discoverable via `outputPaths`), and `EFFECT_NT`/`EFFECT_NTo`
likewise.

## Verification

**Method.** Both DSL blocks were extracted verbatim from `code <- '...'`
(here `AKU_CODE <- '...'`) and POSTed to the qspserver `mrgsolve_api`
container's `/run_simulation`, spaced ~2.5 s apart. Since the model's own
`$MAIN` is deliberately empty of any `<CMT>_0` initial-condition assignment
(by design — see the file's own comment: "init() must win"), and the API's
`SimRequest` schema exposes no field that can set compartment initial values
for the mrgsolve engine (`parameters` only overrides `$PARAM`; `events` is
explicitly deSolve-engine-only per its own schema description, confirmed
inert for this model), every run necessarily starts from the DSL's own
all-zero `$CMT` defaults rather than the R wrapper's physiological
birth-state — irrelevant to verifying *this* refactor, since the original and
refactored DSLs were run from the identical (zero) state with identical
dosing, so any divergence between them would still isolate a genuine
NT-related regression.

Dosing was built directly as `DoseSpec` records (`cmt=1`, matching both
files' first-declared compartment), reproducing the file's own scenarios'
doses (`dose_mg * 1000/329.25` umol, `ii=1` day, `addl` days of repeat
dosing) — S02–S05-style (1/2/4/8 mg/day) and S06-style (10 mg/day) magnitudes.

**A pre-existing, model-wide (not nitisinone-specific) numerical fragility
was found, reproduced identically in the untouched original.** From an
all-zero start, every scenario tested — regardless of dose, and confirmed
also with **zero** drug present — degrades to `NaN` across every
`$CAPTURE`d output within 2–3 days of simulated time, at the *exact same*
time index in the original and the refactored model every time. This matches
the failure mode the file's own header comments describe at length (a
Michaelis-Menten denominator, `KMAPP + CHPP`, passing through zero on
integrator overshoot) — the model's own R wrapper (`sim_aku()`) works around
it with non-default `atol=1e-10, rtol=1e-8, maxsteps=2e6` and, essential here,
physiological (non-zero) initial conditions via `init()`. Confirmed this is
dose-magnitude-independent above a small threshold (0.1 mg and 1 mg and 10 mg
all fail at the same *style* of NaN onset, just at slightly different day
counts consistent with faster-rising `C_NT`/`EFFECT_NT` reaching the failure
condition sooner at higher doses); doses small enough to avoid it (≤0.01 mg)
instead hit a **separate, unrelated** `deSolve`/`lsoda` "h_ = 0" degenerate-
step error that also reproduces identically with **no dosing record at all**
— an artifact of the event table, not of nitisinone's magnitude. The
qspserver API exposes no `atol`/`rtol`/`maxsteps` override (confirmed absent
from its OpenAPI schema), so this model cannot presently be driven to a
multi-year horizon through this API regardless of which DSL (original or
refactored) is used. **Not logged as a new `UPSTREAM_ISSUES.md` entry**: this
is not a compile defect (both DSLs build and manifest cleanly) and the
model's own author already extensively documents the non-default tolerances
and non-zero initial state this file requires — the gap is the qspserver
API's fixed defaults, not a defect in the original model file.

**Result within the finite window both models could run:** exact match,
every time. Four paired runs (1 mg/28 d daily, 10 mg/4 y monthly-observed,
0.1 mg/20 d daily, plus the manifest check) all show **max abs diff = 0.0,
max rel diff = 0.0** across every one of the 37 shared `$CAPTURE`d/compartment
outputs (`CNTo`, `NTGUT`/`GUT_NT`, `NTCEN`/`CENT_NT`, `NTPER`/`PERI_NT`,
`CUMNT`, `UHGA24`, `CHGAo`, `CTYRo`, `CAKUSSI`, `PIGTOT`, and the rest), for
every time point up to the shared `NaN` onset, and the `NaN` onset itself
lands at the identical index in every paired comparison (e.g. index 3 of 30
for the 1 mg/day run, index 2 of 49 for the 10 mg/day run). This is the exact
match expected for archetype 3 (pure PK reorganization) plus a rename-only
Hill term (no fitting was performed anywhere in this file).

## Anything else worth flagging

- The `SRC_INS`/`KI_OAT`/`CLHGA0`/`CLHPPU`/`CLHPLAU`/`VCONJ`/`VRENTYR`
  calibration parameters and the whole Phe/Tyr/HGA pathway are completely
  untouched — they are not nitisinone's own PK/PD block.
- `KI_NT` was deliberately **not** renamed to `EC50_NT`: it is a true
  inhibition constant for an unbounded competitive-inhibition term, not the
  half-maximal concentration of a saturating ratio, so forcing the EC50 label
  onto it would misdescribe the mechanism it actually parameterizes.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`alkaptonuria | Nitisinone (NT)` row.
