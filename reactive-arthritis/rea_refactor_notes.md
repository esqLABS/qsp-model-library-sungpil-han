# Refactor notes — `reactive-arthritis/rea_mrgsolve_model.R`

Scope: all **five** compounds this file gives their own PK compartments and
that have existing rows in `driver-patches/data/compound_perturbation_census.md`
— **NSAID** (naproxen), **SSZ** (Sulfasalazine), **MTX** (Methotrexate),
**TNFI** (a TNF inhibitor, modeled on etanercept — original stem `TNFi`),
and **IL17I** (an IL-17 inhibitor, modeled on secukinumab — original stem
`IL17i`) — all classified "Redirect concentration (clean single site)".
This file models no other pharmacological compound.

## Headline finding: the original does not compile at all — unrelated to any of the five compounds' own PK

Before any renaming could be attempted, `POST /model_manifest` on the
untouched original (via the qspserver `mrgsolve_api` container,
`http://localhost:8007`) failed to build. Both defects are logged in full,
with exact error text, as `translations/UPSTREAM_ISSUES.md` #82:

1. **`$CMT` (`@annotated`) and a separate `$INIT` block jointly redeclare
   all 26 compartments** — mrgsolve 2.0.1 treats `$INIT`'s `NAME = value`
   list as its own compartment declaration, so combining it with `$CMT`
   redeclares every compartment twice ("Duplicated model names"). Same
   defect class as `cervical-cancer`/`urolithiasis`/etc.
2. **`$ODE`'s `dxdt_CHRON_PATH` reads a bare identifier `chlamydia`** that
   is never declared anywhere — not in `$PARAM`, not as a `$CMT`
   compartment, and never set via `param()`/`mrgsim(param=...)` in any of
   the R wrapper's 5 scenarios. A missing declaration, not a name
   collision (`'chlamydia' was not declared in this scope`). The
   surrounding comment ("Chlamydia: small fraction seeds persistent
   synovial compartment") suggests `PATH` may have been the intended read,
   but nothing in the file recovers that intent unambiguously.

Per the fork guide's settled policy for this situation
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"),
both fixes are **syntax-only, non-numeric, non-behavioral**, and applied
directly in the delivered `rea_mrgsolve_model_refactored.R`/`.cpp` (not
just a throwaway scratch copy):

1. `$INIT` deleted; its 26 assignments moved into `$MAIN` using the modern
   `<CMT>_0 = value;` idiom, values unchanged, under the refactor's renamed
   PK compartments (e.g. `GUT_NSAID_0 = 0;`).
2. `chlamydia : 0` added to `$PARAM` — matching the value the original's
   own R wrapper already implicitly used everywhere (it is never set to
   anything else in any of the 5 scenarios), so nothing is invented beyond
   the original's own default-by-omission.

Neither change touches a single equation or parameter value belonging to
any of the five scoped compounds. The tracked `rea_mrgsolve_model.R` is
completely untouched and still carries both defects exactly as written;
the workaround lives only in the `_refactored.R`/`.cpp` siblings.

## Archetype per compound

**All five compounds: Archetype 3 minus peripheral (depot + central,
linear, no back-flux).** Every one of NSAID/SSZ/MTX/TNFI/IL17I in the
original is a plain first-order-absorption, first-order-elimination
one-compartment PK block:

```
dxdt_<STEM>_DEPOT = -ka_<STEM> * <STEM>_DEPOT;
dxdt_<STEM>_C     =  ka_<STEM> * F_<STEM> * <STEM>_DEPOT - ke_<STEM> * <STEM>_C;
```

renamed, per compound, to:

```
dxdt_GUT_<STEM>  = -KA_<STEM> * GUT_<STEM>;
dxdt_CENT_<STEM> =  KA_<STEM> * F_<STEM> * GUT_<STEM> - KE_<STEM> * CENT_<STEM>;
```

**Why `KE_<STEM>` and not `CL_<STEM>`/`V1_<STEM>`:** none of the five
compounds' central compartments in the original are divided by a volume
before use — `NSAID_C`/`SSZ_C`/`MTX_C`/`TNFi_C`/`IL17i_C` are themselves
the concentration state, with elimination governed by a single micro rate
constant (`ke_<stem>`), not a `CL/V` pair. Forcing a `CL_<STEM>`/`V1_<STEM>`
split here would require inventing a volume the original never had. This
is the same situation, and the same naming decision, as
`KE_LOP`/`KE_OND`/`KE_RIF` in
`bile-acid-diarrhea/bam_mrgsolve_model_refactored.R` ("no volume in
original; C_LOP is the raw amount") — `KE_<STEM>` is used as a bespoke
rate-constant name (the guide's table has no dedicated slot for a
combined absorption/elimination micro-constant), and `C_<STEM>` is a
direct alias of `CENT_<STEM>`, not `CENT_<STEM>/V1_<STEM>`.

**NSAID's `Vd_NSAID` (renamed `V1_NSAID`) is dead, preserved as such.** The
original declares `Vd_NSAID : 0.12 : NSAID apparent volume (L/kg, ...)` in
`$PARAM` but never reads it anywhere in `$ODE`/`$TABLE` — `NSAID_C` is
still governed purely by `ka_NSAID`/`F_NSAID`/`ke_NSAID`, exactly like the
other four compounds. `V1_NSAID` is carried over, renamed, still unused
— not invented, not dropped, per the guide's "don't invent, don't drop the
original's own quantities" treatment of dead parameters elsewhere in this
corpus (e.g. `RTOT_TCZ` in `rheumatoid-arthritis`, `CENT_TOBRA` in
`bronchiectasis`).

## TMDD check (IL17I and TNFI)

Both IL17I (modeled on secukinumab) and TNFI (modeled on etanercept) are
antibody biologics, so this refactor explicitly checked for target-mediated
drug disposition. **Neither is modeled as TMDD in this file.** There is no
free-receptor compartment, no drug-receptor complex compartment, and no
`kon`/`koff` binding kinetics anywhere in `$CMT` or `$ODE` for either
compound — each is a plain first-order depot→central PK feeding a plain
`Emax*C/(EC50+C)` ratio, identical in shape to the small-molecule
compounds (NSAID/SSZ/MTX) in the same file. Archetype 3-minus-peripheral
applies to all five compounds uniformly; Archetype 4 (TMDD) was not
needed for any of them.

## Naming

| Role | NSAID | SSZ | MTX | TNFI | IL17I |
|---|---|---|---|---|---|
| Depot | `GUT_NSAID` (was `NSAID_DEPOT`) | `GUT_SSZ` (was `SSZ_DEPOT`) | `GUT_MTX` (was `MTX_DEPOT`) | `GUT_TNFI` (was `TNFi_DEPOT`) | `GUT_IL17I` (was `IL17i_DEPOT`) |
| Central / exposed conc. | `CENT_NSAID` (was `NSAID_C`) | `CENT_SSZ` (was `SSZ_C`) | `CENT_MTX` (was `MTX_C`) | `CENT_TNFI` (was `TNFi_C`) | `CENT_IL17I` (was `IL17i_C`) |
| Absorption rate | `KA_NSAID` (was `ka_NSAID`) | `KA_SSZ` (was `ka_SSZ`) | `KA_MTX` (was `ka_MTX`) | `KA_TNFI` (was `ka_TNFi`) | `KA_IL17I` (was `ka_IL17i`) |
| Bioavailability | `F_NSAID` (unchanged) | `F_SSZ` (unchanged) | `F_MTX` (unchanged) | `F_TNFI` (was `F_TNFi`) | `F_IL17I` (was `F_IL17i`) |
| Elimination (bespoke rate const., no volume) | `KE_NSAID` (was `ke_NSAID`) | `KE_SSZ` (was `ke_SSZ`) | `KE_MTX` (was `ke_MTX`) | `KE_TNFI` (was `ke_TNFi`) | `KE_IL17I` (was `ke_IL17i`) |
| Dead volume param (unused in `$ODE`, same as original) | `V1_NSAID` (was `Vd_NSAID`) | — (none in original) | — | — | — |
| Hill EC50 | `EC50_NSAID` (was `IC50_NSAID`) | `EC50_SSZ` (was `IC50_SSZ`) | `EC50_MTX` (was `IC50_MTX`) | `EC50_TNFI` (was `IC50_TNFi`) | `EC50_IL17I` (was `IC50_IL17i`) |
| Hill Emax | `EMAX_NSAID` (was `Emax_NSAID`) | `EMAX_SSZ` (was `Emax_SSZ`) | `EMAX_MTX` (was `Emax_MTX`) | `EMAX_TNFI` (was `Emax_TNFi`) | `EMAX_IL17I` (was `Emax_IL17i`) |
| Hill gamma (new) | `GAMMA_NSAID = 1` | `GAMMA_SSZ = 1` | `GAMMA_MTX = 1` | `GAMMA_TNFI = 1` | `GAMMA_IL17I = 1` |
| Exposed concentration | `C_NSAID` (new; `= CENT_NSAID`) | `C_SSZ` (new) | `C_MTX` (new) | `C_TNFI` (new) | `C_IL17I` (new) |
| Effect on disease | `EFFECT_NSAID` (was `E_NSAID`) | `EFFECT_SSZ` (was `E_SSZ`) | `EFFECT_MTX` (was `E_MTX`) | `EFFECT_TNFI` (was `E_TNFi`) | `EFFECT_IL17I` (was `E_IL17i`) |

All parameter *values* are copied verbatim from the original — nothing
invented, nothing defaulted, nothing dropped, aside from the two
new-but-non-invented `GAMMA_<STEM> = 1` additions (below) and the
build-fix `chlamydia = 0` addition (unrelated to any compound).

## The Hill interface: five renames-only, no fitting anywhere in this file

Every one of the original's five effect terms was already a plain ratio
with no Hill exponent, e.g.:

```
double E_NSAID = Emax_NSAID * NSAID_C / (IC50_NSAID + NSAID_C + 1e-12);
```

Rename only — pulled out as `EMAX_<STEM>`/`EC50_<STEM>` with the original's
own values, plus a new explicit `GAMMA_<STEM> = 1` (the original had no
exponent at all). `pow(x, 1) == x` makes

```
double EFFECT_NSAID = EMAX_NSAID * pow(C_NSAID, GAMMA_NSAID)
                      / (pow(EC50_NSAID, GAMMA_NSAID) + pow(C_NSAID, GAMMA_NSAID) + 1e-12);
```

arithmetically identical to the original ratio, including the original's
own `+1e-12` guard against division-by-zero at `C=0`, kept verbatim in all
five effect terms. No `nls()` fitting was needed or performed for any
compound in this file.

### Each compound's effect stays a separate named variable, combined only where the disease equations already combined them

`EFFECT_MTX` appears in three separate disease terms exactly as `E_MTX`
did in the original — `Th1_in`, `Th17_in` (scaled ×0.8), `TNF_prod` (scaled
×0.6), and `IL6_prod` (scaled ×0.4) — a rename in place, not a
restructuring; `TNF_prod` likewise multiplies `(1-EFFECT_TNFI)`,
`(1-EFFECT_SSZ*0.5)`, and `(1-EFFECT_MTX*0.6)` together exactly as the
original combined `E_TNFi`/`E_SSZ`/`E_MTX`. Per the guide's "never collapse
several drugs into one shared Hill term" rule, no two compounds' effects
were merged; each keeps its own `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>`
triplet.

## qspserver compatibility

Followed the precedent set by the refactors before this one (e.g.
`bile-acid-diarrhea/bam_refactor_notes.md`): `C_<STEM>` and `EFFECT_<STEM>`
are plain `double` locals computed in `$ODE` (assigning to a
`$PARAM`-declared symbol from `$ODE` is a compile error in this mrgsolve
build), and discoverability is satisfied instead by listing every one of
them in a `$CAPTURE @annotated` block. Confirmed via `POST
/model_manifest` on the qspserver `mrgsolve_api` container that:

- `outputPaths` includes all 5 exposed concentrations (`C_NSAID`, `C_SSZ`,
  `C_MTX`, `C_TNFI`, `C_IL17I`), all 5 effect terms (`EFFECT_NSAID`,
  `EFFECT_SSZ`, `EFFECT_MTX`, `EFFECT_TNFI`, `EFFECT_IL17I`), and the
  original's own output names (`NSAID_obs`, `SSZ_obs`, `MTX_obs`,
  `TNFi_obs`, `IL17i_obs`, `VAS_pain`, `CRP_obs`, ...), kept for backward
  output comparability (74 parameters total, 59 outputs total).
- `parameters` includes every renamed PK parameter and every renamed Hill
  parameter (`EC50_<STEM>`/`EMAX_<STEM>`/`GAMMA_<STEM>` for all five
  compounds) with the original's numeric defaults.
- Compartment ordering is unchanged (position 19–28 in `$CMT`, same as
  the original's compartments 19–28), so 1-based compartment-number
  dosing (as `POST /run_simulation`'s `dosing` field uses) is unaffected
  by the rename.
- `rea_mrgsolve_model_refactored.cpp` (the bare-DSL extraction required
  for `model_content`) is confirmed byte-identical to the quoted `rea_code
  <- '...'` string in `rea_mrgsolve_model_refactored.R`.

**One incidental parse trap found and fixed while drafting this file:** an
early draft's `$PROB` free-text commentary began two lines with a literal
`[refactor]`/`[build-fix]` tag at column 1. mrgsolve's block parser reads
any line starting with `[` as a candidate block-header switch, so those
lines silently became an "invalid blocks found: REFACTOR" warning (not a
hard failure, but wrong) until reworded to "Refactor note:"/"Build-fix
note:" — a purely cosmetic authoring fix in this fork's own commentary
text, not a defect in the original file, so not logged to
`UPSTREAM_ISSUES.md`.

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
regimens through both an in-memory, syntax-fixed copy of the untouched
original (same two build-fixes as above, applied only for the comparison
run — the checked-in `rea_mrgsolve_model.R` was never touched) and
`rea_mrgsolve_model_refactored.cpp`, via the qspserver `mrgsolve_api`
service (`POST /run_simulation`), comparing every shared output
point-by-point. Requests were spaced ~2–3s apart per the guide's
concurrency-safety note.

**Scenario 5's own regimen (NSAID+SSZ+TNFI combination) — SSZ 1000mg
q12h continuous, NSAID 500mg q12h for the first 84 days, ETN 50mg weekly
starting day 84**, shortened to a 100-day (2400h) window from the
original's full 365-day horizon per the guide's solver-budget allowance
(no step-count issues were actually hit at this horizon; shortened mainly
to keep the ~370-event dosing set and request round-trip fast under the
API's `max_concurrent_jobs: 2`): `NSAID_obs`, `SSZ_obs`, `TNFi_obs`,
`VAS_pain`, `swollen_jt`, `CRP_obs`, `TNF_obs`, `IL17_obs`, `IL6_obs`,
`SYNOV_obs`, `CARTDMG_pct`, `TH1_obs`, `TH17_obs`, `TREG_obs` — **max abs
diff = 0 for every output, every time point (772 points)**. This exercises
NSAID, SSZ, and TNFI simultaneously, plus their combined effect on every
disease/immune endpoint in the file.

**None of the file's own 5 scenarios ever dose MTX or IL17i** — both
compounds have full PK/effect equations wired into the disease core, and
the R wrapper even defines `make_mtx_events()`/`make_il17i_events()`
helper functions with sensible clinical defaults (15mg/wk SC MTX; 300mg
loading + monthly SC secukinumab, matching the `$PROB`-commented doses),
but neither helper is ever actually called by any of Scenarios 1–5 (see
"Also found, not a build defect" below). Verification for these two
compounds therefore used the original's own (unused) helper-function
defaults directly, not an invented dosing scenario:

**MTX — `make_mtx_events()`'s own default (15mg SC weekly), shortened to
56 days / 9 doses**: `MTX_obs`, `TH1_obs`, `TH17_obs`, `TREG_obs`,
`TNF_obs`, `IL6_obs`, `VAS_pain`, `SYNOV_obs` — **max abs diff = 0 for
every output, every time point (234 points)**.

**IL17I — `make_il17i_events()`'s own default (300mg loading ×5 doses
then monthly), shortened to 90 days / 7 doses**: `IL17i_obs`, `IL17_obs`,
`TH17_obs`, `VAS_pain`, `SYNOV_obs` — **max abs diff = 0 for every output,
every time point (368 points)**. One early run of this exact comparison
came back with the refactored side pinned at 0 for `IL17i_obs` and `NaN`
for `VAS_pain`/`IL17_obs` while the request itself returned HTTP 200 —
re-running the identical request pair (~3s apart) reproduced a clean
exact match immediately, confirming this was the container's own known
intermittent instability under load (same class of flakiness documented
in `bronchiectasis/bex_refactor_notes.md` and
`translations/UPSTREAM_ISSUES.md` #64), not a property of either model.

**Result: exact match (max abs diff = 0) for every output, every
scenario, every compound.** This is the expected outcome for Archetype 3
minus peripheral applied uniformly to all five compounds — pure PK
reorganization plus `GAMMA=1` Hill renames, with no fitting anywhere in
this file.

## Also found, not a build defect

The original's own R wrapper defines `make_mtx_events()` and
`make_il17i_events()` helper functions (with sensible, correctly-commented
clinical defaults) but never actually invokes either one in any of its 5
named scenarios (`out_s1`...`out_s5`) — MTX and IL17i's PK/effect math is
fully wired into the disease equations and would behave correctly if
dosed, but the shipped script never exercises either compound. Not fixed
(adding a new scenario to the original is out of scope for a rename-only
refactor); logged here and in `UPSTREAM_ISSUES.md` #82 since it directly
affected how these two compounds had to be verified (via the original's
own unused helper defaults rather than one of its named scenarios, see
above).

## Anything else worth flagging

- Compartment ordering is unchanged (`GUT_NSAID` is still compartment 19,
  `GUT_IL17I` is still compartment 27, `CENT_IL17I` is still compartment
  28, etc.), so 1-based compartment-number dosing is unaffected between
  the original and the refactored file.
- Outside the DSL block, the surrounding R script needed five `cmt =`
  string updates in the dosing-event helper functions (`make_nsaid_events`,
  `make_ssz_events`, `make_mtx_events`, `make_tnfi_events`,
  `make_il17i_events`, each now targeting the renamed `GUT_<STEM>` depot)
  and one column-name-list update in Scenario 3's `init()` call (the 10
  raw PK compartment names it restores state into after the antibiotic
  phase). No other line of the R wrapper (scenario logic, plotting
  functions, summary table, calibration checks) differs from the
  original.
- `s3_state` (an intermediate variable built in Scenario 3) is dead code
  in the original — assigned but never read afterward — and needed no
  renaming since its own column-exclusion list only names `_obs`-suffixed
  aliases, which are unchanged; left exactly as dead as it was.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`reactive-arthritis | IL17I`, `reactive-arthritis | MTX`, `reactive-arthritis
| NSAID`, `reactive-arthritis | SSZ`, and `reactive-arthritis | TNFI` rows.
