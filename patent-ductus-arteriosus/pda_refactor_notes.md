# Refactor notes — `pda_mrgsolve_model.R` (ibuprofen, indomethacin, acetaminophen)

**Scope of this pass.** Per `FORK_WORKFLOW_GUIDE.md` Part 2 and the three
existing rows in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
(all three classified "Redirect concentration inside `#define` macro"),
the PK and ductal-COX effect of all three modeled drugs — **Acetaminophen
(APAP)**, **Ibuprofen (IBU)**, **Indomethacin (IND)** — were refactored.
These are the file's only three modeled compounds, so there is nothing
pharmacological left untouched to report; every other part of the model
(prostanoid dynamics, ductal contraction/geometry, wall-O2 bistability,
haemodynamics, lung mechanics, kidney/gut/platelet/bilirubin physiology,
outcome hazards, and the entire 16-scenario R runner) is byte-for-byte
identical to the original except for the handful of lines that reference
the renamed drug identifiers by name (see "Diff scope confirmation").

## Archetype determination

**Ibuprofen and indomethacin: Archetype 2 (no depot, two compartments,
linear) plus one bespoke ductal effect-site compartment each.** Both are
dosed directly into a central compartment (no absorption depot — IV/enteral
lysine dosing goes straight into the central pool), have a peripheral
compartment via `Q`/`V1`/`V2`, and — beyond the guide's four archetypes —
each also has a third, slow first-order **effect-site** compartment
(`IBUE`/`INDE`, mg/L, stored directly in concentration units) representing
slowly-reversible tight binding at the ductal COX channel, distinct from
the plasma/systemic concentration used at the renal/gut/platelet sites.
This is the same kind of extension the guide's own text anticipates
("[expose] two [concentration variables] only when a genuinely different
tissue site matters") and the same shape as the `dka` refactor's bespoke
`EFF_PERI_INS`/`EFF_HEP_INS` precedent (two effect sites, not one).
Clearance for both is postnatally time-varying (`CL_IBU0`/`CLMAT_IBU`,
`CL_IND0`/`CLMAT_IND` — a maturation function of postnatal age, not a
fixed archetype constant), preserved unchanged; the guide's `CL_<STEM>`
slot doesn't have a maturation-function variant, so this is kept as-is
rather than forced into a static-constant shape it doesn't have.

**Acetaminophen: Archetype 3 (depot + central + peripheral, linear).** A
clean textbook match — enteral depot (`APAPG`) with `KA_APAP`/`F_APAP` into
a central/peripheral two-compartment PK, no effect site (acetaminophen's
ductal action is instead peroxide-tone-gated at the concentration itself,
not through a slow-equilibrating site).

## The Hill interface: bespoke, and why

The original combines all three drugs' ductal effect inside one `#define`
macro (`PDA_ALG`), in two genuinely different ways:

1. **Ibuprofen and indomethacin compete for the same arachidonate channel**
   — the original's `occD = UIBUE/KI_IBU + UINDE/KI_IND` is a real,
   physical competitive-occupancy sum (Ariëns/Gaddum-style), and
   `ICHAN = occD/(1+occD)` is the resulting *joint* saturation. This is not
   two independent Hill curves later combined by Bliss-style
   `1-(1-e1)(1-e2)` — the two drugs share one denominator because they
   compete for one binding site. That coupling is real pharmacology, not
   an artifact of the macro, so it cannot be removed by any rename.
2. **Acetaminophen acts at a physically separate peroxidase site** —
   `IPEROX = UAPAP/(UAPAP+ic50D)` was already, on its own, exactly the
   guide's canonical Hill shape (Emax=1, gamma=1, EC50 = a
   peroxide-tone-scaled `ic50D`).

Given (1), forcing ibuprofen and indomethacin into two fully independent
`EMAX_<STEM>·Xᵞ/(EC50ᵞ+Xᵞ)` terms that combine only after being computed
would either be wrong (dropping the shared-site competition) or would
require reinventing the same combination internally anyway. The refactor
therefore keeps the physical coupling but **splits the original's single
`occD` sum into its two named addends before summing them** — a pure
algebraic rename, not a re-fit:

```c
double EFFECT_IBU  = U_IBU_EFF / KI_IBU;   // ibuprofen's own channel-occupancy contribution
double EFFECT_IND  = U_IND_EFF / KI_IND;   // indomethacin's own channel-occupancy contribution
double occD  = EFFECT_IBU + EFFECT_IND;    // <-- the ONE combination point; identical to the original
double ICHAN = occD / (1.0 + occD);
```
`EFFECT_IBU` and `EFFECT_IND` are each a genuine, single-compound function
of that compound's own effect-site concentration (`U_IBU_EFF`/`U_IND_EFF`,
in turn derived from `C_IBU`/`C_IND`) — satisfying "keep each compound's
`EFFECT_<STEM>` separate; combine them only at the point the disease
equations actually use them" as literally as the underlying physics
allows. Each addend individually is *not* itself bounded to \[0,1) the way
a textbook Hill term is (with the other drug's dose at zero it reduces
exactly to `Ci/(Ci+Ki)`, which *is* Hill-shaped; with both nonzero the
individual addends are pre-saturation odds, not fractional occupancies) —
this is disclosed here rather than dressed up with a cosmetic
`EMAX_IBU=1`/`GAMMA_IBU=1` that would misrepresent where the real Emax=1
saturation actually lives (the *combined* `ICHAN`, not each addend). This
is the bespoke handling the guide's "None of these fit" fallback allows,
and it is an exact rename: verification below confirms it reproduces the
original bit-for-bit.

For acetaminophen (2), the promotion is the guide's ordinary "rename, not
refit" case:
```c
double ic50D  = EC50_APAP * (1.0 + perox / KPEROX);
double EFFECT_APAP = EMAX_APAP * pow(U_APAP, GAMMA_APAP)
                    / (pow(ic50D, GAMMA_APAP) + pow(U_APAP, GAMMA_APAP));
```
with new explicit `EMAX_APAP = 1.0` and `GAMMA_APAP = 1.0` (the original's
plain `C/(C+IC50)` form had both implicitly; promoted to named parameters,
same precedent as the guide's own `GAMMA_TCZ=1`/the dka refactor's
`EMAX_INS_UP=1`) — with these values `EFFECT_APAP` reduces algebraically to
exactly the original `IPEROX` expression for any input.

Combination at the point the disease equations use it, unchanged in
structure from the original:
```c
double IPEROX = EFFECT_APAP;
double ICOXD  = 1.0 - (1.0 - ICHAN) * (1.0 - IPEROX);   // channel x peroxidase, same as original
```

**Scope decision — only the ductal (central) pathway got a named
`EFFECT_<STEM>`.** All three drugs also act at three *other* COX sites in
this model — renal (`ICOXK`), gut mucosal (`ICOXG`), and (ibuprofen/
indomethacin only) platelet (`ICOXPLT`) — each with its own site-specific
Ki/IC50 (`KI_IBU_K`, `KI_IBU_G`, `KI_IBU_PLT`, etc.) and used only for
safety/side-effect hazards (NEC, IVH bleeding-time, bilirubin
displacement), not the ductal-closure endpoint the census's
"COX(-1/-2)/peroxidase" target/pathway column names. These three sites
were **not** given their own `EFFECT_<STEM>` variables — they still use
the renamed `U_IBU`/`U_IND`/`U_APAP` (plasma-level, not effect-site,
concentrations) internally, structurally identical to the original. This
mirrors the AAA-Doxycycline precedent in the census (an informational,
non-exposed diagnostic concentration kept alongside the one that's
actually exposed) and keeps the naming-convention promise ("the compound's
effect on disease") scoped to the disease this file is actually about
(patent ductus arteriosus), not every downstream physiological consequence
the drugs happen to also have. Flagged here for whoever next wants to make
the renal/gut/platelet axes independently driveable too — the concentration
variables they read (`U_IBU`, `U_IND`, `U_APAP`) are already exposed and
renamed, so extending this would only require adding three more named
`EFFECT_<STEM>_<SITE>` terms, not restructuring anything.

## Naming applied

| Role | Original | Refactored |
|---|---|---|
| Central compartment, ibuprofen | `IBU1` | `CENT_IBU` |
| Peripheral compartment, ibuprofen | `IBU2` | `PERI_IBU` |
| Ductal effect site, ibuprofen (bespoke) | `IBUE` | `EFF_IBU` |
| Central compartment, indomethacin | `IND1` | `CENT_IND` |
| Peripheral compartment, indomethacin | `IND2` | `PERI_IND` |
| Ductal effect site, indomethacin (bespoke) | `INDE` | `EFF_IND` |
| Absorption depot, acetaminophen | `APAPG` | `GUT_APAP` |
| Central compartment, acetaminophen | `APAP1` | `CENT_APAP` |
| Peripheral compartment, acetaminophen | `APAP2` | `PERI_APAP` |
| Exposed plasma concentration | `CIBU`, `CIND`, `CAPAP` | `C_IBU`, `C_IND`, `C_APAP` |
| Unbound plasma concentration (internal diagnostic) | `UIBU`, `UIND`, `UAPAP` | `U_IBU`, `U_IND`, `U_APAP` |
| Unbound effect-site concentration (internal, bespoke) | `UIBUE`, `UINDE` | `U_IBU_EFF`, `U_IND_EFF` |
| Ductal peroxidase EC50, acetaminophen | `IC50_APAP` | `EC50_APAP` |
| Ductal peroxidase Emax, acetaminophen (new, was implicit) | *(none)* | `EMAX_APAP = 1.0` |
| Ductal peroxidase Hill coefficient, acetaminophen (new, was implicit) | *(none)* | `GAMMA_APAP = 1.0` |
| Compound's ductal effect (bespoke — see above) | *(none, buried in `occD`/`IPEROX`)* | `EFFECT_IBU`, `EFFECT_IND`, `EFFECT_APAP` |
| Volumes, clearances, Q, KA, F, KE0, MW, FU, Ki (all three drugs) | `V1_IBU`, `V2_IBU`, `Q_IBU`, `CL_IBU0`, `CLMAT_IBU`, `FU_IBU`, `MW_IBU`, `KE0_IBU`, `KI_IBU`, … | unchanged — already stem-scoped, matching the convention |
| Renal/gut/platelet site-specific Ki/IC50 (`KI_IBU_K`, `KI_IND_G`, `IC50_APAP_K`, `KI_IBU_PLT`, …) | unchanged | unchanged — deliberately out of `EFFECT_<STEM>` scope, see above |

Every parameter *value* is copied verbatim from the original — nothing
invented, nothing dropped, per the guide's "never invent or default a PK
parameter" rule (the two new parameters, `EMAX_APAP`/`GAMMA_APAP`, encode
values the original's formula already implied, not new pharmacology).

## A pre-existing upstream defect found while verifying, fixed in-scope

Logged as `translations/UPSTREAM_ISSUES.md` #55. The original does not
compile under mrgsolve 2.0.1 for two independent reasons, **both located
inside the three refactored compounds' own PK/PD declarations** (not
shared scaffolding), so per the guide's "when the original doesn't compile
at all" policy, item 4, both were fixed directly as part of this refactor
rather than being left as a bare workaround:

1. Two `$PARAM @annotated` entries (`KE0_IND`'s and `IC50_APAP_K`'s) wrap
   their description onto a second, non-`//`-prefixed line, which
   mrgsolve's annotated-block parser reads as a malformed entry rather
   than continuation text (`Error: improper annotation format`). Fixed by
   folding each onto one line (no value or wording change).
2. Three multi-declarator `double` lines — the ibuprofen/indomethacin/
   acetaminophen central/peripheral-concentration setup lines
   (`double c1i = IBU1/V1_IBU, c2i = IBU2/V2_IBU;` and its IND/APAP
   siblings) — drop the `double` type from the second declarator under
   mrgsolve 2.0.1's `$ODE` preprocessing (`'c2i' was not declared in this
   scope`), the same defect class as issue #37 (`diabetic-ketoacidosis`).
   Fixed by splitting each into two separate `double` statements, same
   values and order.

Both fixes are syntax-only and non-numeric; the exact-match verification
below is against an in-memory copy of the *original* patched with these
same two fixes (the tracked `pda_mrgsolve_model.R` itself is untouched and
still does not compile, per the never-edit-upstream rule), so the
verification result directly demonstrates neither fix changed any
behavior.

## qspserver compatibility

- `C_IBU`, `C_IND`, `C_APAP`, `EFFECT_IBU`, `EFFECT_IND`, `EFFECT_APAP` are
  `$ODE`/`$TABLE`-local `double`s (via the shared `PDA_ALG` macro),
  **not** `$PARAM` entries. Tried the `$PARAM`-with-`=0`-default route
  first per the guide's item 2; it does not compile under mrgsolve 2.0.1
  for the same structural reason already documented in five-plus prior
  calibration runs (`dka`, `ktx`, `sepsis`, `thyroid-eye-disease`,
  `polymyalgia-rheumatica`, `rheumatoid-arthritis`, `primary-sclerosing-
  cholangitis`, …): `$PARAM` values are passed into `$ODE` as
  `const double&`, so assigning to them is a compile error. Discoverability
  is instead satisfied the normal way for this repo's pattern: all six
  names are declared in `$CAPTURE @annotated` (confirmed present in
  `/model_manifest`'s `outputPaths`, see below) under their own bare names
  (this file's `$CAPTURE` is a separate `@annotated` block, not a fused
  `$ODE`/`$TABLE capture NAME = EXPR;` idiom, so no self-shadowing risk
  the way the `dka` file had).
- All 51 `$CMT` compartments (nine renamed, 42 unchanged) and all 181
  `$PARAM` entries (three new: `EC50_APAP`, `EMAX_APAP`, `GAMMA_APAP`;
  everything else unchanged) are confirmed present via `/model_manifest`.
- `$SET`/dosing: the model's own dosing is entirely NM-TRAN-style `ev()`
  bolus/infusion events built by the R wrapper (`rx_ibuprofen`,
  `rx_indomethacin`, `rx_acetaminophen`), not `$SET`-driven or a
  continuous-infusion `$PARAM` input — the API's `dosing` field drives it
  directly and correctly, confirmed by the verification runs below (which
  use exactly this mechanism).

## Diff scope confirmation

`diff -u pda_mrgsolve_model.R pda_mrgsolve_model_refactored.R` touches
exactly: the `$PARAM` COX-enzymology block (rename + 2 new lines + the 2
annotation-continuation fixes), the nine drug `$CMT` lines, the `PDA_ALG`
macro's drug-concentration/COX-inhibition section (renames plus the
`EFFECT_*` restructuring described above), two non-COX-vasoconstriction
lines (`vcind`/`vcibu`), the `NAPQI` line, the three drugs' `$ODE` PK
blocks (renames plus the 2 multi-declarator splits), the `bftgt`/
`indprot`/`HSIP` lines that read a drug concentration directly, the
`$CAPTURE` block (renames plus 3 new `EFFECT_*` entries), and the R
wrapper's `cmt=` string literals plus one `d$CIBU`→`d$C_IBU` reference in
`maturation_table()`. Every other line — all non-drug compartments and
their `$MAIN`/`$ODE` physiology (prostanoid, contraction, wall-O2,
haemodynamics, lung, kidney/gut/platelet/bilirubin, outcome hazards), the
16-scenario definitions, and every scenario-analysis helper function — is
byte-for-byte identical.

## Verification

**Method.** Both DSL blocks (`pda_code <- '...'`) were mechanically
extracted from the original and the refactored sibling (R single-quoted
string, so each `\\` in the source was unescaped to a single `\` before
use — the raw quoted text is not directly compilable C++ macro syntax
without this step) and posted to the qspserver `mrgsolve_api` container
(`http://localhost:8007`, confirmed healthy), using `POST /model_manifest`
to confirm both compile and expose the expected parameters/outputs, and
`POST /run_simulation` to compare every shared `$CAPTURE`d output across
seven scenarios built directly from dosing values already defined in the
original file's own `SCENARIOS` list (not invented ones), run one at a
time (never concurrently) given the API's `max_concurrent_jobs: 2` and its
recent history of crashing under load:

1. **S1** — expectant management, no dosing (baseline/machinery check).
2. **S2** — ibuprofen 10-5-5 mg/kg from postnatal hour 48 (`rx_ibuprofen(48)`).
3. **S6** — indomethacin 0.2-0.1-0.1 mg/kg from hour 48 (`rx_indomethacin(48)`).
4. **S8** — IV acetaminophen 15 mg/kg q6h × 3 d from hour 48 (`rx_acetaminophen(48)`).
5. **S9** — acetaminophen (as S8) with chorioamnionitis (`SEPSIS=1`) — the
   scenario that specifically exercises the peroxide-tone-scaled
   `ic50D`/`EC50_APAP` path.
6. **S11** — ibuprofen + acetaminophen combination, both regimens merged
   and time-sorted (`ev_bind(rx_ibuprofen(48), rx_acetaminophen(48))`) —
   exercises `ICOXD`'s `ICHAN × IPEROX` cross-term with both non-zero
   simultaneously.
7. **S13** — indomethacin + early hydrocortisone (`rx_indomethacin(48)`,
   `HCORT=1`) — exercises the SIP/IVH hazard terms that read `EFF_IND`
   (originally `INDE`) directly, outside the COX-inhibition macro.

Each scenario was run over a 240 h (10-day) window at delta=1 for the
main comparison; S2 was additionally re-run over the full 90-day (2160 h)
window used by the original's own `run_scenarios()` at delta=3, to confirm
the shorter window wasn't hiding a longer-horizon divergence. All shared
`$CAPTURE`d columns plus all 51 compartments were compared (renamed
1:1 per the naming table above; `EFFECT_IBU`/`EFFECT_IND`/`EFFECT_APAP`
have no original-side counterpart and are additions, not compared).

**Result: exact match. Maximum absolute and relative deviation observed =
0.0 across all seven scenarios, both time windows, and every shared
`$CAPTURE`d/compartment output.** This is the expected outcome for a pure
structural reorganization/rename with no Hill-refitting, per the guide's
tolerance table — including for the bespoke `EFFECT_IBU`+`EFFECT_IND`
competitive-occupancy split, confirming the algebraic identity holds
exactly, not just approximately.

**Solver-step budget: not a limiting factor for this model.** Unlike some
prior files in this fork (e.g. `distal-renal-tubular-acidosis`, which needs
`maxsteps=5e6` against the API's default 20000), this model's full 90-day
window at delta=3 (724 output points) ran to completion under the API's
default settings without a step-count error — the 240 h window used for
the bulk of the comparisons was chosen for turnaround speed, not because
the API couldn't handle the full window.

**What was not directly exercised: ibuprofen and indomethacin
co-administered simultaneously.** None of the original file's 16 named
scenarios doses both drugs at once (they are alternative NSAID choices in
this model, not a combination regimen), so the `EFFECT_IBU + EFFECT_IND`
summation inside `occD` is only exercised with one addend at zero across
all seven verified scenarios — with the other addend at zero, the sum
correctly reduces to the single-drug case, and the combination line itself
is a direct, one-line copy of the original's `occD` expression (see "The
Hill interface" above), but a scenario with *both* drugs simultaneously
nonzero was not run since the guide directs using the file's own defined
scenarios rather than inventing new dosing. Flagged here rather than
silently assumed equivalent for that untested combination.

## Anything else worth flagging

- This file's own header already states it was built by comparing against
  an independent Python/LSODA reference implementation
  (`pda_reference_model.py`) during original development, and lists six
  real defects found and fixed *in the original's own development history*
  (a 56×-too-large ductal resistance gain, a lung-water/PVR feedback loop,
  a wrong tone-combination functional form, a threshold fitted only to
  untreated data, a wrong-magnitude wall-O2 form, and a same-day-window
  delta-creatinine bug). None of that history is part of this refactor's
  scope (it predates this pass and is orthogonal to the PK/PD naming
  work) — noted here only so a reviewer doesn't mistake this refactor's
  silence on those six items as having missed them.
- `pda_reference_model.py` and `pda_shiny_app.R`/`pda_shiny_app_en.R` are
  untouched upstream-tracked files, not part of this refactor's
  deliverables, and were not read for correctness beyond confirming they
  don't need updating (they reference the original `pda_mrgsolve_model.R`,
  not the new `_refactored.R` sibling, which is expected — a `_refactored`
  sibling is never wired into the original's own app/reference chain).
