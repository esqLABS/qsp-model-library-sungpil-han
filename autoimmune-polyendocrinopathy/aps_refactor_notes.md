# Refactor notes — `autoimmune-polyendocrinopathy/aps_mrgsolve_model.R`

Compounds refactored: **Cyclosporine A (CSA)**, **Hydrocortisone (HC)**,
and the **JAK inhibitor tofacitinib (JAKI)** — the three rows in
`driver-patches/data/compound_perturbation_census.md` classified
`autoimmune-polyendocrinopathy | CSA`, `| HC`, and `| JAKI`, all "Redirect
concentration (clean single site)". Abatacept (`Drug_Aba`) and Rituximab
(`Drug_RTX`) are the file's other two immunomodulators and are completely
untouched — same compartments, parameters, and equations, byte-for-byte.

## A necessary build-compatibility fix, unrelated to any compound's PK

The untouched original does not compile under mrgsolve 2.0.1 at all:
`$CMT` names all 22 compartments (with trailing comment annotations) and a
separate `$INIT` block re-declares the same 22 names with starting values,
in identical order. mrgsolve treats `$INIT` as its own
compartment-declaring block, so having both is a duplicate declaration —
confirmed via `POST /model_manifest` against the untouched original alone
(`Duplicated model names: AIRE_func AutoT_pool ... Drug_HC`). This affects
every compartment in the file, not just CSA/HC/JAKI's three, and is the
same defect class already logged for `distal-renal-tubular-acidosis`
(#38) and `type1-diabetes` (#42).

Per the guide's settled policy for a non-compiling original
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"),
this is fixed directly in the delivered `aps_mrgsolve_model_refactored.R`
(not just a scratch workaround): the `$CMT` block is deleted outright,
`$INIT` is kept exactly as the original wrote it. Because `$CMT` and
`$INIT` already listed the same 22 names in the same order, this is a
pure syntax fix — every compartment keeps its name, its 1-based index, and
its original starting value; nothing numeric changes. Logged as
`translations/UPSTREAM_ISSUES.md` #51. The tracked `aps_mrgsolve_model.R`
itself is untouched and still carries the defect exactly as written.

## Archetype determination

All three compounds share the same shape: **a single PK compartment, no
depot, first-order elimination — Archetype 1, bespoke variant.** The
compartment's own state is already dimensionally a *concentration*
(ng/mL for CSA/JAKI, µg/dL for HC), not an amount: the original computes
each compound's zero-order input rate by dividing the dose term by a
volume (`ka_HC*F_HC*HC_dose/Vd_HC`, `ka_JAKi*F_JAKi*JAKi_dose/Vd_JAKi`) —
or, for CSA, by a hardcoded body weight of 70 — directly at the point
where the dose enters the ODE, rather than accumulating an amount that a
separate `C_<STEM> = CENT_<STEM>/V1_<STEM>` division would later convert.
Forcing an amount/volume split here (to match the archetype-1 template
literally) would require multiplying the input term back up by a volume
the original never applied, which is not a rename — it would change the
arithmetic. So the refactor keeps the single-compartment structure
(matches Archetype 1) but treats the state as the concentration directly:
`C_<STEM> = CENT_<STEM>` is a bare alias, not a division. This is
disclosed as bespoke per the guide's "None of these fit" fallback, rather
than forcing a division that isn't there in the original.

One genuine inconsistency in the original, preserved exactly rather than
fixed: `Vd_CsA` (renamed `V1_CSA`) is declared as a `$PARAM` but never
used anywhere in the original's own `CsA_input`/`CsA_elim` arithmetic
(which uses a hardcoded `70` instead of `Vd_CsA`) — `V1_CSA` stays
declared-but-inert in the refactor too, exactly matching the original's
own dead parameter.

There is no depot compartment, `ev()`/`data_set()` event, or continuous
zero-order infusion variable computed independently of a dose parameter
for any of the three compounds — each compound's `<STEM>_input` is a
plain algebraic function of its own `<STEM>_dose` `$PARAM`, gated
`(<STEM>_dose > 0) ? ... : 0`, evaluated fresh every `$ODE` step. This
does not match the guide's "continuous infusion input, not an event" edge
case either (that case is for a genuinely time-varying window function
like `win()`; here the input is a flat rate for as long as the dose
parameter is nonzero) — it is simply a linear one-compartment PK model
with its dosing entering as a constant per-day rate rather than a bolus.
Preserved exactly. (The R-side `simulate_scenario()` additionally layers
an `ev(..., rate=-2)` bolus/infusion event on top of this for each
compound; that mechanism is untouched except for renaming its `cmt=`
target string to match the renamed compartments — see "Diff scope" below.)

### Naming applied

| Role | CSA | HC | JAKI |
|---|---|---|---|
| Compartment (state = concentration directly) | `CENT_CSA` (was `Drug_CsA`) | `CENT_HC` (was `Drug_HC`) | `CENT_JAKI` (was `Drug_JAKi`) |
| Bioavailability | `F_CSA` (was `F_CsA`) | `F_HC` (unchanged) | `F_JAKI` (was `F_JAKi`) |
| Absorption rate | `KA_CSA` (was `ka_CsA`) | `KA_HC` (was `ka_HC`) | `KA_JAKI` (was `ka_JAKi`) |
| Elimination rate constant | `CL_CSA` (was `k_CsA_clear`) | `CL_HC` (was `k_HC_clear`) | `CL_JAKI` (was `k_JAKi_clear`) |
| Volume term in the input formula | `V1_CSA` (was `Vd_CsA` — unused, see above) | `V1_HC` (was `Vd_HC`) | `V1_JAKI` (was `Vd_JAKi`) |
| Dose driver | `CsA_dose` (unchanged) | `HC_dose` (unchanged) | `JAKi_dose` (unchanged) |
| Exposed concentration | `C_CSA` | `C_HC` | `C_JAKI` |
| Hill EC50 | `EC50_CSA` (was `IC50_CsA`) | — (no Hill; see below) | `EC50_JAKI` (was `IC50_JAKi`) |
| Hill Emax | `EMAX_CSA = 1.0` (new) | — | `EMAX_JAKI = 1.0` (new) |
| Hill gamma | `GAMMA_CSA = 1.0` (new) | — | `GAMMA_JAKI = 1.0` (new) |
| Effect | `EFFECT_CSA` (was `E_CsA`) | — (see below) | `EFFECT_JAKI` (was `E_JAKi`) |

Dose driver parameters (`CsA_dose`/`HC_dose`/`JAKi_dose`) are left
unchanged, same precedent as `kidney-transplant-rejection`'s
`TCZON`/`TCZSTART`/`TCZN`: they are not one of the naming convention's PK
structural roles, and the R-side `make_scenario()`/`simulate_scenario()`
functions already key off these exact names.

## Hill interface

### CSA and JAKI — rename, not a refit

Both original ratios are already the bare `C/(IC50+C)` shape with an
implicit `Emax=1`, `gamma=1` (no separate Emax constant, no exponent):

```
double E_CsA  = (Drug_CsA  > 0) ? Drug_CsA  / (IC50_CsA  + Drug_CsA)  : 0;
double E_JAKi = (Drug_JAKi > 0) ? Drug_JAKi / (IC50_JAKi + Drug_JAKi) : 0;
```

`EMAX_CSA`/`EMAX_JAKI`/`GAMMA_CSA`/`GAMMA_JAKI` are new named parameters
(all `= 1.0`), reproducing the original ratios exactly, not a fit:

```
EFFECT_CSA  = (C_CSA  > 0) ? EMAX_CSA*pow(C_CSA,  GAMMA_CSA)  / (pow(EC50_CSA,  GAMMA_CSA)  + pow(C_CSA,  GAMMA_CSA))  : 0;
EFFECT_JAKI = (C_JAKI > 0) ? EMAX_JAKI*pow(C_JAKI, GAMMA_JAKI) / (pow(EC50_JAKI, GAMMA_JAKI) + pow(C_JAKI, GAMMA_JAKI)) : 0;
```

Each is consumed in exactly the places the original's `E_CsA`/`E_JAKi`
were, with the downstream weighting literals **left untouched** (these
weights are not part of either compound's own Emax — they are separate,
per-target-organ physiological ceilings the disease equations apply where
they consume the shared ratio, same pattern already established for
tocilizumab in `kidney-transplant-rejection`):

- `Immuno_suppress = (1 - 0.85*EFFECT_CSA) * (1 - 0.80*E_Aba) * (1 - 0.80*E_RTX) * (1 - 0.70*EFFECT_JAKI)` ($MAIN)
- `beta_attack = ... * (1 - EFFECT_CSA*0.6)` ($ODE, pancreatic beta cells — CSA only)
- `adrenal_attack = ... * (1 - EFFECT_JAKI*0.4)` ($ODE, adrenal gland — JAKI only)
- `thy_attack = ... * (1 - EFFECT_JAKI*0.4)` ($ODE, thyroid — JAKI only)

Only the `E_CsA`/`E_JAKi` tokens were touched on these lines; `E_Aba`,
`E_RTX`, and every other term on them belong to the two untouched
compounds and are unchanged.

### HC — bespoke: no Hill term is fabricated

Hydrocortisone has no `Emax`/`EC50` anywhere in the original. Its action
is a direct, additive replacement of the hormone the disease itself fails
to produce (`HC_exogenous = Drug_HC`, then `Cortisol_eff = Cortisol_c +
HC_exogenous`, and `$TABLE`'s `cortisol_total = Cortisol_c + Drug_HC`) —
not a saturating receptor-occupancy effect. Per the guide's "a clean,
non-standard structure beats a standard structure that's wrong," no
`EFFECT_HC` is invented. `C_HC` is exposed (the guide's requirement 2/4 —
the single point an external covariate could later substitute in) and
consumed exactly where the original consumed `Drug_HC`: `HC_exogenous`,
`Cortisol_eff`, and the `$TABLE` `cortisol_total` capture. The census's
listed target (glucocorticoid receptor) describes HC's real pharmacology,
but the original model represents that pharmacology as mass-balance
replacement, not receptor kinetics — this refactor doesn't add kinetics
the original never had.

## A second necessary technical fix, introduced by this refactor's own new `$CAPTURE` entries

`C_CSA`, `EFFECT_CSA`, `C_JAKI`, `EFFECT_JAKI`, and `C_HC` must be
`$CAPTURE`d (via `capture NAME = expr;` in `$TABLE`) for qspserver
discoverability (guide requirement 4). Under mrgsolve 2.0.1, a `double
NAME = ...;` declaration for one of these names in **`$MAIN`** collides
with mrgsolve's own auto-declared member for that same `$CAPTURE`d name —
confirmed empirically via `POST /model_manifest`:

```
84:11: error: redefinition of 'capture {anonymous}::C_CSA'
   84 |   capture C_CSA;
28:10: note: 'double {anonymous}::C_CSA' previously declared here
```

This extends the same underlying mrgsolve-2.0.1 auto-promotion mechanism
already documented for `distal-renal-tubular-acidosis` (`$ODE`-scoped
collision) to `$MAIN` as well. Fix: in `$MAIN`, `C_CSA`/`EFFECT_CSA`/
`C_JAKI`/`EFFECT_JAKI`/`C_HC` are assigned **bare** (no `double`), relying
on mrgsolve's own auto-declaration from being `$CAPTURE`d in `$TABLE`;
`$TABLE` keeps its explicit `capture NAME = expr;` form, recomputed
independently from the compartment state (`CENT_CSA`/`CENT_JAKI`/
`CENT_HC`), matching this file's own existing convention of recomputing
every `$TABLE` quantity from state rather than reusing a `$MAIN`/`$ODE`
local (`cortisol_total`, `CaCorr`, etc. already do this). Confirmed
empirically both ways via `/model_manifest` (fails with `double` in
`$MAIN`, builds cleanly without it) — this is purely a declaration-site
fix, introduced by this refactor's own new capture entries (the original
never captured these names), so it is not logged to
`UPSTREAM_ISSUES.md`, per the same reasoning `distal-renal-tubular-acidosis`
used for its own `$ODE`-scoped version of the same mechanism.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted DSL
(`http://localhost:8007`):

- `F_HC`, `KA_HC`, `CL_HC`, `V1_HC`, `HC_dose`, `F_CSA`, `KA_CSA`,
  `CL_CSA`, `V1_CSA`, `CsA_dose`, `EC50_CSA`, `EMAX_CSA`, `GAMMA_CSA`,
  `F_JAKI`, `KA_JAKI`, `CL_JAKI`, `V1_JAKI`, `JAKi_dose`, `EC50_JAKI`,
  `EMAX_JAKI`, `GAMMA_JAKI` all appear in the manifest's `parameters`
  (fixed values, unchanged from the original).
- `CENT_CSA`, `CENT_HC`, `CENT_JAKI` (compartments), plus `C_CSA`,
  `C_HC`, `C_JAKI`, `EFFECT_CSA`, `EFFECT_JAKI` (via `$CAPTURE`) all
  appear in the manifest's `outputPaths` (34 total, including the disease
  compartments/captures and the untouched Abatacept/Rituximab states).
- The quoted `code <- '...'` block inside the delivered
  `aps_mrgsolve_model_refactored.R` was confirmed **byte-identical** to
  the DSL text actually exercised in every verification call below (a
  Python re-extraction + string-equality check, run after the file was
  written) — the DSL was never hand-edited after being tested. No
  separate `.cpp` file is kept as a deliverable; the extraction was a
  verification-only step.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`), comparing the (build-compat-fixed, see above)
original against `aps_mrgsolve_model_refactored.R`, reproducing the
original's own scenarios `2_HRT_Only` (HC only), `3_HRT_CsA` (HC+CSA),
and `6_HRT_JAKi` (HC+JAKI) — `AIRE_mut_sev=0.90`, `HC_dose=20`, plus each
scenario's own `CsA_dose=3.5`/`JAKi_dose=10` — over the scenarios' own
full 5-year duration (`end=1825`, `delta=7`, 261 points; no shortening
needed, this model did not hit the API's default solver step-count
budget at this duration). Every shared disease-side output was compared
(`AIRE_func`, `AutoT_pool`, `Treg_pool`, all four `AutoAb_*`,
`Adrenal_fn`, `Cortisol_c`, `PTG_fn`, `PTH_plasma`, `Ca_serum`,
`Beta_mass`, `Insulin_p`, `Glucose_p`, `Thyroid_fn`, `TSH_plasma`,
`FT4_plasma`, `cortisol_total`, `HbA1c_est`, `CaCorr`, `T3_est`,
`ACTH_est`, `APS_components`), plus each compound's own PK state under
its old/new name (`Drug_CsA`/`CENT_CSA`, `Drug_HC`/`CENT_HC`,
`Drug_JAKi`/`CENT_JAKI`).

**An unrelated qspserver infrastructure outage interrupted this
verification pass mid-way and is recorded here for transparency, since it
briefly affected every session sharing this container, not just this
one.** While probing build variants, one `/run_simulation` call crashed
the API's R subprocess worker with a glibc heap-corruption signal
(`free(): invalid next size (fast)`); this corrupted a *shared* mrgsolve
compiled-model cache index (`mrgmod_cache.RDS`, used by every inline
`model_content` submission on that container, not keyed to this model
specifically) and left the container itself briefly unresponsive to
`docker restart`/`docker kill`. Once the container was recreated and the
one corrupted cache-index file was removed (a stale index only — no
compiled `.so` artifacts or other sessions' data were touched), both
`/model_manifest` and `/run_simulation` returned to normal immediately,
with no changes needed to either model file.

**Result: exact match, not just near-exact.** Across all three scenarios,
**maximum absolute difference was exactly 0.0** for every one of the 20
shared disease-side outputs and all three compounds' own PK states
(`Drug_CsA` vs `CENT_CSA`, `Drug_HC` vs `CENT_HC`, `Drug_JAKi` vs
`CENT_JAKI`), across the full 261-point, 5-year time grid, in every
scenario. This is the expected outcome for Archetype 1 (pure structural
reorganization plus a bare rename of an already-Emax=1/gamma=1 Hill
ratio, no fit) per the guide's tolerance table. `EFFECT_CSA` ranged
0–0.5 in `3_HRT_CsA` and `EFFECT_JAKI` ranged 0–0.0223 in `6_HRT_JAKi`,
confirming the renamed Hill arithmetic is actively engaged (not trivially
zero) while still matching the original bit-for-bit.

## Anything else flagged

- No compound other than CSA, HC, and JAKI was touched. Abatacept and
  Rituximab — their compartment, parameters, and equations — are
  byte-identical to `aps_mrgsolve_model.R`.
- **Diff scope in the R harness (outside the DSL block):** the
  `simulate_scenario()` function's three `ev(cmt="Drug_HC"/"Drug_CsA"/
  "Drug_JAKi", ...)` calls were updated to `cmt="CENT_HC"/"CENT_CSA"/
  "CENT_JAKI"` so the refactored sibling's own R scenario-driver code
  stays internally consistent with the renamed compartments — a rename of
  the dosing target only, not a behavioral change (same compartment, same
  1-based index, same dose amounts/timing). No other R-side code
  (plotting, summary tables) references these compartments by name.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`autoimmune-polyendocrinopathy | CSA`, `| HC`, and `| JAKI` rows.

## Discoverability fix

A corpus-wide discoverability audit found `C_CSA` and `C_JAKI` were not
written as single contiguous `double C_<STEM> = <expr>;` statements anywhere
in the file. Both were bare-assigned once in `$MAIN` (`C_CSA = CENT_CSA;` /
`C_JAKI = CENT_JAKI;`) and exposed separately via `capture C_CSA =
CENT_CSA;` / `capture C_JAKI = CENT_JAKI;` in `$TABLE` — a legitimate,
working pattern (not a bug), but not literal-text-discoverable by tooling
that regexes for `double C_<STEM> = ...;`.

**First attempt (per the standing instruction for this pattern) failed to
compile.** Adding a *new*, separate `double C_CSA = CENT_CSA;` /
`double C_JAKI = CENT_JAKI;` line in `$TABLE` immediately before the
existing `capture C_CSA = ...;` / `capture C_JAKI = ...;` lines — leaving
those `capture` lines untouched, as instructed — was tried first (including
a variant wrapped in its own nested `{ ... }` C++ block, in case ordinary
brace scoping would avoid the clash). Both attempts failed identically when
posted to qspserver's `mrgsolve_api` `/model_manifest`:

```
81:11: error: redefinition of ‘capture {anonymous}::C_CSA’
   81 |   capture C_CSA;
79:10: note: ‘double {anonymous}::C_CSA’ previously declared here
   79 |   double C_CSA;
```

mrgsolve's own DSL parser auto-declares a same-named class member for
*every* `double NAME = ...;` and every `capture NAME = ...;` it finds
anywhere in the block's text — it is a text-scan, not a C++-scope-aware
pass, so wrapping the `double` line in braces does not help; two
declaration mechanisms for the same name is a genuine build failure, not a
tolerance issue. This matches this file's own pre-existing `$MAIN` comment
warning about the identical collision for a `double` re-declaration there.

**Actual fix applied:** converted `C_CSA`/`C_JAKI` from the inline
`capture NAME = expr;` mechanism to a genuine `double NAME = expr;`
declaration in `$TABLE` (identical formula, identical value — just declared
rather than captured-with-auto-declare), and added an explicit `$CAPTURE`
block at the end of the DSL (this file previously had none — it relied
entirely on inline `capture` statements) listing `C_CSA C_JAKI` so both
remain reported outputs:

```
$TABLE
...
double C_CSA = CENT_CSA;
double C_JAKI = CENT_JAKI;
capture C_HC  = CENT_HC;
capture EFFECT_CSA = ...;
capture EFFECT_JAKI = ...;
...
$CAPTURE
C_CSA C_JAKI
```

`C_HC` and every other `capture ... = ...;` line (`EFFECT_CSA`,
`EFFECT_JAKI`, `cortisol_total`, etc.) are untouched — out of scope for this
fix and still using the original mechanism; `$MAIN`'s bare `C_CSA =
CENT_CSA;` / `C_JAKI = CENT_JAKI;` assignments (lines ~202, ~206) are also
untouched and continue to work because `$CAPTURE`-listed names are
auto-declared and shared across `$MAIN`/`$ODE`/`$TABLE` the same way inline
`capture` names were.

**Verification:** `Rscript -e 'parse(...)'` succeeds (an early draft had a
straight apostrophe in "mrgsolve's" inside a new comment, which broke the
R string — caught by `parse()` and fixed with a curly apostrophe per this
guide's Part 2 rule 4, before any qspserver call). Extracted DSL posted to
qspserver's `mrgsolve_api` `/model_manifest` compiled cleanly with both
`C_CSA` and `C_JAKI` listed in `outputPaths`, and `EC50_CSA`/`EC50_JAKI` in
the parameter manifest. `/run_simulation` with a simple bolus dose into
`CENT_CSA` (cmt 19, amt 250) and `CENT_JAKI` (cmt 22, amt 5) at t=0, 0–30
(delta 1), run against both the pre-edit (`git show HEAD:...`) and
post-edit DSL, produced **bit-identical** `CENT_CSA`, `CENT_JAKI`, `C_CSA`,
`C_JAKI`, `EFFECT_CSA`, and `EFFECT_JAKI` columns (max abs diff = 0 across
the full time grid). No numeric behavior changed by either the discovered
build conflict or its fix.
