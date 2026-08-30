# Refactor notes — `familial-mediterranean-fever/fmf_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
(`FORK_WORKFLOW_GUIDE.md`, Part 2), all **three** PK-bearing compounds this
file models were rewritten, per their existing rows in
`driver-patches/data/compound_perturbation_census.md` (all three classified
"Redirect concentration (clean single site)"): **Colchicine (COL)**,
**Anakinra (ANA)**, and **Canakinumab (CANA)**. Rilonacept has no PK
compartment of its own in the original (`RILO_dose` is a bare covariate
multiplied directly into an Emax term, the same shape as Perindopril in
`abdominal-aortic-aneurysm/aaa_mrgsolve_model.R`) and is therefore not a
refactor target; it is byte-identical to the original. Every disease-side
equation (PYRIN inflammasome, IL-1b/IL-18/SAA/CRP, neutrophils, attack
severity, amyloidosis/eGFR) is untouched.

## Archetype per compound

**Colchicine (COL): Archetype 3 (depot + central + peripheral, linear) plus
a bespoke leukocyte target-site accumulation compartment.** The original
already models colchicine intracellular accumulation as a real fourth
compartment (`LEU_COL`, unchanged name -- already conformed), fed from
plasma by first-order uptake/release (`k_leu_on`/`k_leu_off`, renamed
`K_LEU_ON_COL`/`K_LEU_OFF_COL`) rather than as an algebraic partition of
plasma concentration. Both of the file's own colchicine effect terms read
`Cl_col_conc` (the leukocyte compartment, renamed `C_COL`), never
`Cp_col` (plasma, renamed `C_COL_PLASMA`, kept only for reporting -- no PD
equation reads it). This is exactly the guide's stated exception ("two
[exposed concentrations] only when a genuinely different tissue site
matters"), the same shape as `TISSUE_DOXY`/`C_DOXY` in
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`.

**Anakinra (ANA): Archetype 3 minus peripheral (depot + central, linear).**
Two compartments in the original (`SC_ANA`, `CENT_ANA`), first-order SC
absorption, linear elimination. Renamed depot `SC_ANA` -> `GUT_ANA`
(the fork's depot role is always `GUT_<STEM>` regardless of physiological
route, matching the precedent in `rheumatoid-arthritis/ra_refactor_notes.md`
which renamed tocilizumab's own `SC_TCZ` -> `GUT_TCZ` the same way).

**Canakinumab (CANA): Archetype 3 (depot + central + peripheral, linear).**
Three compartments in the original (`SC_CANA`, `CENT_CANA`, `PERI_CANA`),
first-order SC absorption, linear elimination, no TMDD structure. Renamed
depot `SC_CANA` -> `GUT_CANA`, same reasoning as ANA above.

## Anakinra and canakinumab: comment claims TMDD, code does not implement it

The task instructions specifically asked to check canakinumab (an
anti-IL-1beta antibody) for TMDD/receptor-binding kinetics, since large
biologics like this very often are modeled with genuine receptor-binding
ODEs elsewhere in this corpus (e.g. tocilizumab's `REC_FREE_TCZ`/
`COMPLEX_TCZ` in `rheumatoid-arthritis/ra_mrgsolve_model.R`). Checking the
actual equations, not the comments:

- The file's own calibration-notes header comment block states "Anakinra:
  ... TMDD-based binding to IL-1R1" and "Canakinumab: ... TMDD to free
  IL-1b".
- The actual code for **both** compounds is plain linear PK: `SC_ANA`/
  `CENT_ANA` and `SC_CANA`/`CENT_CANA`/`PERI_CANA` are ordinary first-order
  compartments (`dxdt_CENT_ANA = ka_ana*F_ana*SC_ANA - (CL_ana/V_ana)*
  CENT_ANA`, etc.) -- there is no free-receptor compartment, no
  drug-receptor complex compartment, no `KON`/`KOFF`/`RTOT` anywhere in
  either compound's block, and neither compound's own effect term
  (`E_ana`, `E_cana`) is anything but a plain concentration-driven
  Emax/EC50 ratio (`Emax*C/(IC50+C)`).

Per the guide's explicit instruction to keep a mechanistically rich PK
model as rich as the original (never flatten it) -- the mirror image also
holds: never invent receptor-binding richness a file's code does not
actually contain, even where its own comment claims it. Both compounds
were refactored as the linear PK they are actually coded as (Archetype 3
variants), not upgraded to Archetype 4. This comment/code mismatch is
logged as its own item in `translations/UPSTREAM_ISSUES.md` (see below),
since it is exactly the class of discrepancy ("a PK model that is actually
TMDD when a comment implies otherwise" -- here inverted: the comment
implies TMDD, the code is not TMDD) the fork's shared rules call out for
logging.

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `ka_col` | `KA_COL` | 1.2 (1/h) | absorption rate |
| `F_col` | `F_COL` | 0.45 | bioavailability |
| `CL_col` | `CL_COL` | 18.0 (L/h) | clearance |
| `V1_col` | `V1_COL` | 120.0 (L) | central volume |
| `Q_col` | `Q_COL` | 60.0 (L/h) | inter-cpt CL |
| `V2_col` | `V2_COL` | 480.0 (L) | peripheral volume |
| `k_leu_on` | `K_LEU_ON_COL` | 2.5 (1/h) | leukocyte uptake (bespoke role) |
| `k_leu_off` | `K_LEU_OFF_COL` | 0.08 (1/h) | leukocyte release (bespoke role) |
| `Cp_col` (local) | `C_COL_PLASMA` | -- | plasma conc., reporting only |
| `Cl_col_conc` (local) | `C_COL` | -- | **the exposed concentration** (leukocyte) |
| `IC50_col_neu` | `EC50_COL_NEU` | 8.0 (ng/mL) | Hill EC50 |
| `Emax_col_neu` | `EMAX_COL_NEU` | 0.85 | Hill Emax |
| -- (none) | `GAMMA_COL_NEU` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| `IC50_col_PYRIN` | `EC50_COL_PYRIN` | 5.0 (ng/mL) | Hill EC50 |
| `Emax_col_PYRIN` | `EMAX_COL_PYRIN` | 0.60 | Hill Emax |
| -- (none) | `GAMMA_COL_PYRIN` (new) | 1 | Hill exponent |
| `E_col_neu` | `EFFECT_COL_NEU` | -- | colchicine effect on neutrophil migration |
| `E_col_PYRIN` | `EFFECT_COL_PYRIN` | -- | colchicine effect on PYRIN/ASC |
| `SC_ANA` (cmt) | `GUT_ANA` | -- | depot (was SC route, see above) |
| `CENT_ANA` | unchanged | -- | already conformed |
| `ka_ana` | `KA_ANA` | 0.40 (1/h) | absorption rate |
| `F_ana` | `F_ANA` | 0.95 | bioavailability |
| `CL_ana` | `CL_ANA` | 1.8 (L/h) | clearance |
| `V_ana` | `V1_ANA` | 8.5 (L) | volume |
| `Cp_ana` (local) | `C_ANA` | -- | **the exposed concentration** |
| `IC50_ana` | `EC50_ANA` | 50.0 (ng/mL) | Hill EC50 |
| `Emax_ana` | `EMAX_ANA` | 0.95 | Hill Emax |
| -- (none) | `GAMMA_ANA` (new) | 1 | Hill exponent |
| `E_ana` | `EFFECT_ANA` | -- | anakinra effect on IL-1 signaling |
| `SC_CANA` (cmt) | `GUT_CANA` | -- | depot (was SC route, see above) |
| `CENT_CANA`, `PERI_CANA` | unchanged | -- | already conformed |
| `ka_cana` | `KA_CANA` | 0.012 (1/h) | absorption rate |
| `F_cana` | `F_CANA` | 0.70 | bioavailability |
| `CL_cana` | `CL_CANA` | 0.18 (L/h) | clearance |
| `V1_cana` | `V1_CANA` | 4.5 (L) | central volume |
| `Q_cana` | `Q_CANA` | 0.4 (L/h) | inter-cpt CL |
| `V2_cana` | `V2_CANA` | 3.0 (L) | peripheral volume |
| `Cp_cana` (local) | `C_CANA` | -- | **the exposed concentration** |
| `IC50_cana` | `EC50_CANA` | 10.0 (ug/mL) | Hill EC50 |
| `Emax_cana` | `EMAX_CANA` | 0.98 | Hill Emax |
| -- (none) | `GAMMA_CANA` (new) | 1 | Hill exponent |
| `E_cana` | `EFFECT_CANA` | -- | canakinumab effect on IL-1b neutralization |

All parameter *values* are copied verbatim from the original. `USE_COL`/
`USE_ANA`/`USE_CANA` (the file's own per-compound dosing switches) are
unchanged -- they gate whether each compound's effect term is active at
all, orthogonal to the naming convention's own role table, same treatment
this fork has given other dosing switches (`RILO_dose` here, `TCZIN`/
`NIMPO` in prior runs).

## Hill interface: rename, not a fit -- for all three compounds

Every one of the five effect terms this file computes (`E_col_neu`,
`E_col_PYRIN`, `E_ana`, `E_cana`, plus the derived `IL1b_block`) was
already exactly `Emax*C/(IC50+C+1e-10)` in the original -- a plain Hill
ratio with an implicit exponent of 1 and a `1e-10` divide-by-zero guard.
Per the guide's rename rule, `EMAX_<STEM>`/`EC50_<STEM>` were pulled out
with the original's own values, `GAMMA_<STEM> = 1` was added for each (no
compound had an explicit Hill exponent), and the `1e-10` guard was kept in
the same position (`pow(EC50, GAMMA) + pow(C, GAMMA) + 1e-10`) so the
rewritten expression is bit-identical to the original for GAMMA=1, not an
approximation of it. No `nls()` fit was needed or performed for any of the
three compounds.

**Colchicine has two separate downstream effect terms, not one.** The
original computes distinct Emax/IC50 pairs for neutrophil-migration
inhibition (`IC50_col_neu`/`Emax_col_neu`) and PYRIN/ASC inhibition
(`IC50_col_PYRIN`/`Emax_col_PYRIN`) from the same leukocyte concentration.
These are genuinely two different physiological actions of the same drug
(not "several drugs on a combined multi-drug expression," the case the
guide's "combine only where disease equations actually use them" targets),
so they were kept as two independent named terms, `EFFECT_COL_NEU` and
`EFFECT_COL_PYRIN`, following the same pattern as Encaleret's two effect
terms in `hypoparathyroidism/hypopt_refactor_notes.md`.

## A pre-existing upstream build defect, unrelated to any of the three compounds

Neither `fmf_mrgsolve_model.R` nor a pure rename of it compiles under
mrgsolve 2.0.1. Five layered defects were found, none inside any of the
three refactored compounds' own blocks (all five sit in disease-PD scaffolding
shared by the whole file) -- confirmed reproducing identically from the
untouched original via the qspserver `mrgsolve_api` container
(`POST /model_manifest`):

1. **`$CMT` (bare, unannotated) and `$INIT` jointly redeclare all 22
   compartments** -- same defect class as issues #27/#34/#36/#42/#48:
   `Error in validObject(.Object): invalid class "mrgmod" object:
   Duplicated model names: GUT_COL CENT_COL ... eGFR`.
2. **Six disease-baseline `$PARAM` values collide with their own
   compartment's auto-generated `<CMT>_0` init symbol**, once defect 1's
   fix (moving `$INIT` into `$MAIN` via the modern `<CMT>_0 = value;`
   idiom) is applied: `IL18_0`, `SAA_0`, `CRP_0`, `Neu_circ_0`,
   `Neu_tis_0`, and `eGFR_0` are ordinary baseline-value parameters used
   only to seed the same-named compartment's initial condition (`IL18 =
   IL18_0` etc. in the original `$INIT`) -- but mrgsolve auto-generates an
   init-condition symbol also named `IL18_0` for compartment `IL18` (and
   symmetrically for the other five), so the two same-named symbols
   collide, the same class of incidental collision as `MMP9_0`/`MMP2_0` in
   `abdominal-aortic-aneurysm/aaa_refactor_notes.md`. (`IL1b_0`, used the
   same way to seed `IL1b_mat`, does *not* collide -- the compartment is
   named `IL1b_mat`, not `IL1b`.)
3. **`$CAPTURE` repeats eleven compartment names already in `$CMT`**
   (`IL1b_mat IL18 SAA CRP Neu_circ Neu_tis Att_sev AA_dep eGFR ASC
   Casp1`) -- same defect class as issues #44/#45/#46/#48:
   `compartment should not be in $CAPTURE: IL1b_mat,IL18,SAA,CRP,
   Neu_circ,Neu_tis,Att_sev,AA_dep,eGFR,ASC,Casp1`.
4. **`$MAIN` and `$ODE` both declare an identically-named local
   `double k_phos_eff = ...;`** -- mrgsolve 2.0.1 hoists every block-local
   `double` into the same anonymous-namespace scope, so the two
   declarations collide: `error: redefinition of 'double
   {anonymous}::k_phos_eff'`. The `$MAIN` copy is dead code (never read
   after being computed, in either block or capture); only the `$ODE`
   copy is actually used by `k_phos_drug`.
5. **`$MAIN`'s `if(NEWIND <= 1) { _nid++; }` references `_nid`**, an
   internal mrgsolve counter symbol not provided under this build:
   `error: '_nid' was not declared in this scope`. `_nid` appears nowhere
   else in the file -- this block is entirely dead code (its increment is
   never read, captured, or used to affect any equation).

**Verification workaround, applied to the delivered
`fmf_mrgsolve_model_refactored.R` per this fork's settled build-compat
policy (`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at
all")**, since a delivered file that cannot be run through the API defeats
the point of the refactor:
1. The `$INIT` block was deleted; its 22 assignments were moved into
   `$MAIN` using the `<CMT>_0 = value;` idiom, values copied verbatim.
2. The six colliding baseline parameters were renamed with a `_BASE`
   suffix (`IL18_0_BASE`, `SAA_0_BASE`, `CRP_0_BASE`, `Neu_circ_0_BASE`,
   `Neu_tis_0_BASE`, `eGFR_0_BASE`), and the new `$MAIN` init line for
   each reads the renamed parameter (`IL18_0 = IL18_0_BASE;` etc.) --
   same numeric value, same discoverability as a covariate, just a
   symbol rename plus a location move, not a new invented parameter.
   `IL1b_0` (no collision) is untouched.
3. The eleven duplicated names were dropped from `$CAPTURE` -- mrgsolve
   already reports every compartment's trajectory regardless of
   `$CAPTURE`, so nothing about what is reported changes.
4. The dead, unused `$MAIN` copy of `double k_phos_eff = ...;` was
   deleted; the `$ODE` copy (the one actually used) is untouched.
5. The dead `if(NEWIND <= 1) { _nid++; }` block was deleted.

None of these five fixes touches a single numeric value, parameter,
compartment initial condition, or equation -- confirmed by the exact-match
verification results below, which is the actual proof (not just the
description above) that these are syntax-only changes. The untouched
original `fmf_mrgsolve_model.R` still carries all five defects forward
unfixed, per the never-edit-upstream rule. Logged as **upstream issue #49**
in `translations/UPSTREAM_ISSUES.md` (also covering the anakinra/
canakinumab comment-vs-code TMDD mismatch as a second, non-build-blocking
item in the same entry).

## Verification

**Method.** Both files' embedded mrgsolve `fmf_model_code <- '...'` DSL
blocks were mechanically extracted (regex on the assignment, verbatim
quoted text), the original's copy was additionally patched with the
build-compat workaround above (in-memory only, for this verification run
-- never written back to the tracked `fmf_mrgsolve_model.R`), and both
were POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`/model_manifest` then `/run_simulation`), which
compiles and runs each DSL block directly with mrgsolve 2.0.1
server-side -- no local R/mrgsolve install used.

**A shared-infrastructure interruption occurred mid-verification and is
recorded here for transparency, since it is not a defect in this model.**
The qspserver `mrgsolve_api` container's on-disk build cache is shared
across concurrent sessions (several other disease-model refactors were
running against the same container at the same time, evidenced by
`UPSTREAM_ISSUES.md` entries #45-#48 logged during this same window). A
stale, incompatible-R-version cache-index file
(`pkpd_indirect/mrgmod_cache.RDS`, dated a week before this session) made
every model compile fail with `ReadItem: unknown type 47`, regardless of
model content -- confirmed by testing a trivial, freshly-unique one-line
model, which failed identically. This surfaced mid-run as a single
scenario (Anakinra) returning wildly divergent, non-physiological output
between the two models; that specific result was diagnosed as corrupted
API output (not a refactor bug) and discarded rather than reported. Once
the coordinator removed the stale cache-index file, all four scenarios
below were (re-)run cleanly against a healthy service and every result
below is from that clean run.

**Results: exact match on all four of the file's own dosing scenarios.**
Each was run for the scenario's own full stated duration/dosing
(`end = 8760 h`, `delta = 24 h`, 367 timepoints), no shortening needed --
this model's own step count stayed well within the API's default
`maxsteps` budget throughout.

1. **Scenario 2, "Colchicine 0.5 mg BID"** (500,000 ng into
   `GUT_COL`/cmt 1, q12h x731 doses, `USE_COL=1`): **exact match**, max
   abs/rel deviation **0.0** across all 33 shared `$CAPTURE`/compartment
   outputs.
2. **Scenario 3, "Colchicine 1.0 mg QD"** (1,000,000 ng into `GUT_COL`,
   q24h x366 doses, `USE_COL=1`): **exact match**, max abs/rel deviation
   **0.0** across all 33 shared outputs.
3. **Scenario 4, "Anakinra 100 mg SC QD"** (100,000 ng into
   `SC_ANA`/`GUT_ANA`, cmt 5, q24h x366 doses, `USE_ANA=1`): **exact
   match**, max abs/rel deviation **0.0** across all 33 shared outputs.
   (This is the scenario that briefly showed a corrupted, non-reproducible
   mismatch during the infrastructure interruption above; the clean rerun
   is unambiguous.)
4. **Scenario 5, "Canakinumab 150 mg SC Q8W"** (150 mg into
   `SC_CANA`/`GUT_CANA`, cmt 7, q1344h x7 doses, `USE_CANA=1`): **exact
   match**, max abs/rel deviation **0.0** across all 33 shared outputs.

(Scenario 1, "No Treatment," is trivially identical by construction --
all three `USE_*` switches at 0, no dosing -- and was not separately
POSTed.) Compartment indices (`GUT_COL`=1 ... `PERI_CANA`=9) are identical
between the original and the refactored model, since only compartment
*names* changed, never their `$CMT` declaration order. This is the
expected, and achieved, result for three compounds that are all pure
structural renames (Archetypes 3/3-minus-peripheral, no Hill-fitting) per
the guide's tolerance table for Archetypes 1-3: "expect a near-exact
match... anything beyond floating-point-scale deviation means a bug" --
here the match is exact, not merely near-exact.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`KA_COL`, `F_COL`, `CL_COL`, `V1_COL`, `Q_COL`, `V2_COL`,
`K_LEU_ON_COL`, `K_LEU_OFF_COL`, `EC50_COL_NEU`, `EMAX_COL_NEU`,
`GAMMA_COL_NEU`, `EC50_COL_PYRIN`, `EMAX_COL_PYRIN`, `GAMMA_COL_PYRIN`,
`KA_ANA`, `F_ANA`, `CL_ANA`, `V1_ANA`, `EC50_ANA`, `EMAX_ANA`,
`GAMMA_ANA`, `KA_CANA`, `F_CANA`, `CL_CANA`, `V1_CANA`, `Q_CANA`,
`V2_CANA`, `EC50_CANA`, `EMAX_CANA`, `GAMMA_CANA` all appear in the
manifest's `parameters` with their original numeric defaults. `C_COL`,
`C_COL_PLASMA`, `C_ANA`, `C_CANA`, `EFFECT_COL_NEU`, `EFFECT_COL_PYRIN`,
`EFFECT_ANA`, `EFFECT_CANA` are state-derived (computed in `$ODE` from
compartment concentrations), so -- per the same reasoning established for
`EFFECT_NIM`/`EFFECT_ATR`/`EFFECT_TCZ` in prior refactors in this fork --
they cannot also be `$PARAM` entries; all eight appear in the manifest's
`outputPaths` via the file's own (renamed) `$CAPTURE` list, confirmed
discoverable. Unlike some prior files in this fork, no `o`-prefixed
`$TABLE` alias was needed here -- this file's original `$CAPTURE` already
worked by listing already-declared `$ODE` local names directly (the same
mechanism, just with old names), so the renamed identifiers slot into the
same mechanism with no extra indirection.

No separate `.cpp` extraction file was left behind -- extraction was
in-memory only, used to build the verification requests above and then
discarded, per the workflow guide.

## Anything else flagged

- No compound other than Colchicine, Anakinra, and Canakinumab was
  touched. Rilonacept's parameters (`RILO_dose`, `IC50_rilo`,
  `Emax_rilo`) are byte-identical to the original.
- The R-side scenario list, `run_scenario()`, the severity-comparison
  loop, the summary-metrics table, every plot, and the sensitivity-analysis
  block were updated only where they name a renamed compartment or
  parameter (`"SC_ANA"` -> `"GUT_ANA"`, `"SC_CANA"` -> `"GUT_CANA"` in the
  `ev()` calls; `Cp_col`/`Cl_col_conc` -> `C_COL_PLASMA`/`C_COL` in
  `p_col_pk`; `IC50_col_neu` -> `EC50_COL_NEU` in the sensitivity loop) --
  same dosing amounts, same timing, same plot logic throughout.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`familial-mediterranean-fever | ANA`, `familial-mediterranean-fever |
CANA`, and `familial-mediterranean-fever | COL` rows.
