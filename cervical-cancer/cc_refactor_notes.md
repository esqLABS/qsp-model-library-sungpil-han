# Refactor notes — `cervical-cancer/cc_mrgsolve_model.R`

Five compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md`, all classified
"Redirect concentration (clean single site)": **Bevacizumab (BEV)**,
**Cisplatin (CIS)**, **Paclitaxel (PAC)**, **Pembrolizumab (PEM)**, and
**TV-ADC / tisotumab vedotin (TVADC)**. These are the only compounds
modeled in this file — nothing else was touched.

## Archetypes

- **CIS, PAC, BEV, PEM** — all four are **Archetype 2** (no depot,
  two-compartment, linear elimination), confirmed against the actual
  equations: `dxdt_CENT = -(CL+Q)/V1*CENT + Q/V2*PERI`, `dxdt_PERI =
  Q/V1*CENT - Q/V2*PERI`, IV bolus dosing directly into the central
  compartment. **No TMDD anywhere in this file** — Bevacizumab and
  Pembrolizumab (the two compounds explicitly checked for it, per real-
  world target-mediated disposition/checkpoint biology) both have plain
  linear PK with no receptor/complex compartment. VEGF neutralization is a
  first-order mass-action sink on a separate `VEGF` biomarker compartment
  (`VEGF_bind_BEV = kbind_BEV * C_BEV * VEGF`, was `BEV_effect`), not a
  dynamic receptor-complex system; anti-angiogenic tumor-growth inhibition
  is a separate, empirical Emax/EC50 ratio directly on `C_BEV`. Checkpoint
  inhibition is not concentration-driven at all in this model — see
  "No concentration-driven `EFFECT_PEM`" below.
- **TVADC (tisotumab vedotin)** — **bespoke**: an Archetype-2 conjugate PK
  base (`CENT_TVADC`/`PERI_TVADC`) plus a third, non-archetypal compartment
  (`MMAE_TVADC`, was `MMAE_free`) for the released free cytotoxic payload.
  This doesn't match any of the guide's four archetypes — it isn't TMDD
  (no receptor pool, no `KON`/`KOFF`/`RTOT`), and the payload compartment
  isn't a peripheral distribution volume; it's fed by a named fraction
  (`FR_MMAE_TVADC`, was a hardcoded `0.4`) of the conjugate's own clearance
  flux, with its own independent first-order decay (`k_dec_MMAE_TVADC`).
  This genuinely reflects ADC pharmacology (antibody-drug conjugate PK
  disposition is distinct from the released small-molecule payload's own
  disposition) and was kept as-is rather than forced into a mold — a
  clean bespoke structure beats a standard structure that's wrong.

## Renaming (values unchanged from the original)

| Original | Refactored |
|---|---|
| `CIS_C1` / `CIS_C2` | `CENT_CIS` / `PERI_CIS` |
| `PAC_C1` / `PAC_C2` | `CENT_PAC` / `PERI_PAC` |
| `BEV_C1` / `BEV_C2` | `CENT_BEV` / `PERI_BEV` |
| `PEMBRO_C1` / `PEMBRO_C2` | `CENT_PEM` / `PERI_PEM` |
| `TV_ADC_C1` / `TV_ADC_C2` | `CENT_TVADC` / `PERI_TVADC` |
| `MMAE_free` | `MMAE_TVADC` |
| `CL_TV`/`V1_TV`/`Q_TV`/`V2_TV` | `CL_TVADC`/`V1_TVADC`/`Q_TVADC`/`V2_TVADC` |
| `k_dec_MMAE` | `k_dec_MMAE_TVADC` |
| `Pt_DNA` | `ADDUCT_CIS` |
| `BEV_effect` | `VEGF_bind_BEV` (rename only, identical formula) |
| `eff_ICI` | `EFFECT_PEM` (rename only, identical formula) |
| `pac_eff` (`PAC_C1/(PAC_C1+100.0)`) | `EFFECT_PAC` (Hill form made explicit) |
| `MMAE_free/(MMAE_free+1.0)` inline ratio | `EFFECT_TVADC` (Hill form made explicit) |
| `BEV_C1/(BEV_C1+10.0)` inline ratio | `EFFECT_BEV` (Hill form made explicit) |

New (not in the original, for the named Hill interface; `EMAX=1`/`GAMMA=1`
make explicit a shape the original already had implicitly as a plain
ratio — **rename, not a refit**): `EMAX_PAC=1`, `EC50_PAC=100.0`,
`GAMMA_PAC=1`; `EMAX_BEV=1`, `EC50_BEV=10.0`, `GAMMA_BEV=1`;
`EMAX_TVADC=1`, `EC50_TVADC=1.0`, `GAMMA_TVADC=1`. Also new:
`FR_MMAE_TVADC=0.4`, naming a literal `0.4` fraction that was hardcoded
inline in the original's `dxdt_MMAE_free` line (value unchanged).

`CL_CIS`/`V1_CIS`/`Q_CIS`/`V2_CIS`, `CL_PAC`/`V1_PAC`/`Q_PAC`/`V2_PAC`,
`CL_BEV`/`V1_BEV`/`Q_BEV`/`V2_BEV`, and `CL_PEM`/`V1_PEM`/`Q_PEM`/`V2_PEM`
already matched the naming convention and are unchanged.

**Stem harmonization (`PEMBRO` -> `PEM`):** the original used two
different stems for the same drug — `$PARAM` already named
`CL_PEM`/`V1_PEM`/`Q_PEM`/`V2_PEM`, while `$CMT` named the compartments
`PEMBRO_C1`/`PEMBRO_C2`. This is exactly the naming chaos this refactor
exists to remove (per the guide's own tocilizumab example) — harmonized
to the shorter `PEM` stem, matching what the file's own `$PARAM` block
already used, so `CENT_PEM`/`PERI_PEM`.

**Stem disambiguation (`TV` -> `TVADC`):** the original's PK parameters for
tisotumab vedotin used the bare stem `TV` (`CL_TV`, `V1_TV`, `Q_TV`,
`V2_TV`), which reads confusingly next to the unrelated `TV` = tumor-volume
compartment in the same file (no actual symbol collision, since the params
are `CL_TV` etc., not literally `TV` — but genuinely confusing to a reader
and to the naming convention's `C_<STEM>`/`EFFECT_<STEM>` pattern).
Renamed to the `TVADC` stem throughout.

## Exposed concentration: `C_<STEM>` is a direct alias, not `CENT/V1`

The guide's canonical template computes the exposed concentration as
`C_TCZ = CENT_TCZ / V1_TCZ`. That does **not** apply here: this
original's own dosing functions already divide the bolus dose by `V1`
before adding it to the central compartment (e.g. `dose_cisplatin`:
`amt = dose_mg/V1`), and the original's own `$TABLE` captured the central
compartment directly as "the concentration" (`capture CIS_Conc = CIS_C1;`,
no further `/V1`). The `$CMT`-declared compartment is therefore already
concentration-scaled by construction, for all five compounds. Dividing by
`V1` again in the refactor would silently double-scale every PK output and
change the original's numeric behavior. So here, `C_<STEM> = CENT_<STEM>`
is a direct alias — confirmed correct by the verification below (exact
match against the original's own `*_Conc` outputs at every timepoint).

## No concentration-driven `EFFECT_PEM`

Pembrolizumab's disease effect in this model (`eff_ICI`, renamed
`EFFECT_PEM`) is `ICI_flag * CPS_high` — a binary treatment-presence ×
biomarker-eligibility switch. It is **not** a function of `C_PEM` at all;
the checkpoint-inhibition boost to CD8+ T cells (`ICI_effect`), the
ICI-enhanced tumor kill (`kill_ICI`), and the ICI-boosted HPV clearance
(`HPV_clear_eff`) all read `EFFECT_PEM`, never `C_PEM`. `C_PEM` is still
fully exposed (renamed, captured, discoverable via `/model_manifest`) for
future external driveability, but nothing in this file's own disease
equations currently consumes it. This looks like a deliberate modeling
simplification given the original's own comment that pembrolizumab
achieves "near-saturating receptor occupancy at 200mg q3w" (i.e., the
author treated exposure as effectively binary once dosed) rather than an
authoring error, so it is disclosed here, not "fixed" by inventing a new
concentration-dependent effect the original never had (out of scope for a
rename-only refactor, per the guide).

## Hill interface: renames, not fits

All disease-facing effect terms in this file were already plain
concentration/payload ratios — no ODE-solved receptor-occupancy kinetics
to approximate, so every one of these is a rename, confirmed exact by the
verification below:

- **PAC**: `pac_eff = PAC_C1/(PAC_C1+100.0)` is `1*C^1/(EC50^1+C^1)` with
  implicit `Emax=1`, `gamma=1`. `EMAX_PAC=1`/`GAMMA_PAC=1` were added to
  make this explicit; `EFFECT_PAC` computes the identical ratio on
  `C_PAC`.
- **BEV**: the inline `BEV_C1/(BEV_C1+10.0)` term inside `kill_bev` is
  likewise `1*C^1/(EC50^1+C^1)`. `EMAX_BEV=1`/`GAMMA_BEV=1` were added;
  `EFFECT_BEV` computes the identical ratio on `C_BEV`.
- **TVADC**: `adc_eff = MMAE_free/(MMAE_free+1.0)` is likewise
  `1*X^1/(EC50^1+X^1)`, but with `X = MMAE_TVADC` (the released payload),
  **not** `C_TVADC` (the conjugate concentration) — deliberately, since
  the whole point of modeling a separate payload compartment is that the
  cytotoxic effect is driven by free intratumoral MMAE, not circulating
  ADC. `EMAX_TVADC=1`/`EC50_TVADC=1.0`/`GAMMA_TVADC=1` were added;
  `EFFECT_TVADC` computes the identical ratio on `MMAE_TVADC`.
- **CIS**: `Pt_DNA` (renamed `ADDUCT_CIS`) is a linear, non-saturating
  turnover state (`dxdt_ADDUCT_CIS = k_adduct*C_CIS - k_repair*ADDUCT_CIS`),
  not a Hill/Emax ratio — there is no saturating nonlinearity here to
  approximate, so `EFFECT_CIS = ADDUCT_CIS` is a direct alias, exactly
  reproducing the original's `kill_Pt = k_kill_Pt * Pt_DNA * TV`. No
  `nls()` fitting was needed or attempted for any of the five compounds —
  none of them has an ODE-solved, saturating effect kinetics that would
  require it.

## `$GLOBAL` for `C_<STEM>`/`EFFECT_<STEM>` (not `$PARAM`)

Per the guide's qspserver compatibility requirement #2, these should
ideally live in `$PARAM` for direct `/model_manifest` discoverability.
Not done here, for the same reason already documented in the AMD,
membranous-nephropathy, and breast-cancer refactors: mrgsolve 2.0.1
compiles `$PARAM` members as **read-only references** inside `$ODE`, so a
value recomputed every timestep from state cannot also be declared in
`$PARAM`. All ten (`C_CIS`, `C_PAC`, `C_BEV`, `C_PEM`, `C_TVADC`,
`EFFECT_CIS`, `EFFECT_PAC`, `EFFECT_BEV`, `EFFECT_PEM`, `EFFECT_TVADC`)
are predeclared as `double`s in `$GLOBAL` and listed in `$CAPTURE`:
visible in every simulation's output columns and in `/model_manifest`'s
`outputPaths`, just not in its `parameters` list.

## A refactor-introduced staleness bug, found and fixed during verification

Early verification (before the fix below) found a real, non-floating-point
discrepancy: at a timestep where a dose event and a requested output row
coincide (`t=0`, cisplatin bolus), the qspserver `mrgsolve_api` run
reported that exact row using the state *after* the event but *before*
`$ODE` was re-invoked. `C_CIS` (a `$GLOBAL` double only ever written
inside `$ODE`) therefore read `0` instead of the correct `4.5333` at that
one row — a genuine bug introduced by computing the exposed concentration
in `$ODE` and only reading it back in `$TABLE`, something the original
never had a chance to exhibit because it computes `capture CIS_Conc =
CIS_C1;` by reading the bare compartment directly in `$TABLE`, every time,
regardless of `$ODE`'s call history.

**Fix:** every `C_<STEM>`/`EFFECT_<STEM>` is now *also* recomputed
directly from the live `$CMT` state inside `$TABLE` (redundantly with the
identical computation still needed inside `$ODE` for cross-compound
coupling during integration, e.g. `VEGF_bind_BEV`, `kill_Pt`). This closes
the gap entirely — confirmed by the verification below, which now shows
an exact match (including at every simultaneous dose/observation row) for
all five compounds across all six of the original's own scenarios. This
fix is disclosed here rather than in `UPSTREAM_ISSUES.md`, since it isn't
a defect in the original at all — it was introduced by an early draft of
this refactor's own `$GLOBAL`+`$ODE`-only design and fixed before delivery.

## Pre-existing upstream build defects (fixed syntax-only in the delivered file)

**The original does not compile under mrgsolve 2.0.1 at all**, for three
layered, build-blocking reasons plus one non-blocking numerical fragility,
all unrelated to any of the five refactored compounds' own PK math.
Logged in full, with exact error text, as
[`UPSTREAM_ISSUES.md` #76](../translations/UPSTREAM_ISSUES.md). Summary:

1. **`$CMT`+`$INIT` jointly redeclare all 19 compartments** (`Duplicated
   model names`). Fixed: the `$INIT` block is dropped; its 19 assignments
   become `<CMT>_0 = value;` lines in `$MAIN` (identical values, renamed
   compartment identifiers).
2. **Nine `$ODE` lines clamp a compartment state directly**
   (`if(MMAE_free<0) MMAE_free=0;` and eight more), which mrgsolve 2.0.1
   rejects as `assignment of read-only reference`. Confirmed (via an
   isolated minimal reproduction, and via #36/#73's prior findings) that
   this was **always a behavioral no-op**, in any mrgsolve version — only
   `dxdt_*` feeds the integrator, so reassigning the bare state inside
   `$ODE` can never affect the next solver step or the reported
   trajectory. All nine were deleted outright; confirmed dead, not a
   numeric change, by the exact-match verification below.
3. **Three `$PARAM` names collide with mrgsolve's auto-generated
   `<CMT>_0` init symbols**: `SCCAg_0`, `HPVload_0`, `CD8T_0` (all three
   unused dead parameters — the real `$INIT` values `8.0`/`5.0`/`1.0` are
   separate hardcoded literals, not sourced from these params at all) vs.
   the auto-init symbol for compartments `SCCAg`, `HPVload`, `CD8T`.
   Renamed to `SCCAg0`/`HPVload0`/`CD8T0` (no other usages existed).

All three fixes are syntax-only and non-numeric, applied **directly to
the delivered `cc_mrgsolve_model_refactored.R`**, per the guide's settled
policy for a non-compiling original — not just to a scratch copy, since a
"verified" deliverable nobody can run afterwards defeats the point of the
refactor. The checked-in `cc_mrgsolve_model.R` is left completely
untouched and still carries all three defects exactly as written. The
same three fixes were applied identically to an in-memory-only scratch
copy of the original purely so it could build for the comparison below.

4. **Non-blocking: `TV` underflows through zero and NaNs via
   `log(TV_max/TV)` once radiotherapy is active.** The original's own
   `S2` scenario (cisplatin CCRT, `RT_flag=1`) drives `TV` to `~3.39e-22`
   by `t=83` and `NaN` at `t=84`, propagating into every disease-side
   output thereafter. This is the scenario the (already-dead, per #2
   above) `if(TV<0.01) TV=0.01;` clamp was clearly meant to guard —
   meaning this scenario has apparently never actually run to completion,
   only appearing to work because the file never built at all before now.
   **Confirmed identical in both the original and the refactored file**:
   both NaN at the exact same timestep (`t=84`), with identical finite
   values through `t=83` at all 15 shared outputs. Not fixed (guarding it
   would change the original's own numeric behavior — out of scope for a
   rename-only refactor); the verification below reports the exact match
   over the finite prefix and the identical NaN onset for the remainder.

## Verification

**Method.** Both the (Defects-1-3-patched, in-memory-only) original and
`cc_mrgsolve_model_refactored.R`'s embedded DSL were extracted as bare
mrgsolve DSL text and run through the qspserver `mrgsolve_api` container
at `http://localhost:8007` (`POST /model_manifest`, `POST
/run_simulation`), spaced ~2s apart per request. `/model_manifest`
confirmed the refactored file compiles cleanly and lists every
`C_<STEM>`/`EFFECT_<STEM>` in `outputPaths`.

All six of the original file's own dosing scenarios were run, exactly as
defined in its own `dose_*()`/`ev_seq()` calls (dose amounts, `V1`-scaled
`amt`, `time` sequences all copied verbatim, translated to the API's
dosing-record convention — compartment indices are unchanged 1-based
positions, since the refactor renamed identifiers but did not add, drop,
or reorder any compartment), over the model's own 730-day horizon
(`end=730`, `delta=1`):

1. **S1** — Untreated (natural history)
2. **S2** — Cisplatin CCRT (RTOG-90-01 style, weekly cisplatin x6 + RT)
3. **S3** — CCRT + Pembrolizumab (KEYNOTE-A18)
4. **S4** — Chemo + Bevacizumab, recurrent/metastatic (GOG-240)
5. **S5** — Tisotumab vedotin monotherapy, recurrent/metastatic (innovaTV 301)
6. **S6** — Chemo + Bevacizumab + Pembrolizumab, recurrent/metastatic (KEYNOTE-826)

(Multi-drug scenarios' dosing records were time-sorted before submission —
a request-format requirement of the API's dosing path, not a change to
the original's own dosing logic, which builds each drug's event block
independently.)

All 15 shared `$CAPTURE`d outputs (`CIS_Conc`, `PAC_Conc`, `BEV_Conc`,
`PEMBRO_Conc`, `TVADC_Conc`, `MMAE_lvl`, `VEGF_free`, `PtDNA_rel`,
`RT_damage`, `TumorVol`, `SCCAg_lvl`, `HPV_rel`, `CD8T_rel`, `PDL1_rel`,
`TV_change`) were compared point-by-point across the full time grid for
all six scenarios.

**Result: exact match, max abs diff = 0.0, for every output, every
scenario** — S1/S4/S5/S6 over their full 731/763/751/798-point grids;
S2/S3 (the two RT-active scenarios) over the finite prefix up to the
shared NaN onset described in Defect 4 above (`t=0..83`, identical NaN
onset at `t=84` confirmed on both sides). No solver step-count issues or
timeouts were encountered at the model's own full 730-day horizon. This
is the expected outcome per the guide's tolerance rule for pure structural
reorganization with no Hill-fitting — every renamed/promoted effect term
computes the identical arithmetic on the identical (renamed) inputs, once
the refactor-introduced `$TABLE`-staleness bug described above was found
and fixed.

## Anything else worth flagging

- The `$MAIN` block's `EFFECT_PEM = ICI_flag * CPS_high;` assignment
  looks redundant with the identical line now also present in `$TABLE`
  (added as part of the staleness fix) — it is genuinely redundant (both
  compute the same flag-only expression, never state-dependent, so there
  is no staleness risk for this one specifically), but left in `$MAIN` as
  well since it was the original's own placement (`eff_ICI` was computed
  in `$MAIN`) and removing it would be an unnecessary additional change.
- The R-side scenario/plotting code was updated only where it referenced
  a renamed identifier: every `dose_*()` function's `cmt=` argument
  (`"CIS_C1"` -> `"CENT_CIS"`, etc.). All dosing amounts, intervals,
  `ev_seq()` calls, and all six scenario definitions are otherwise
  identical to the original.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`cervical-cancer | Bevacizumab`, `cervical-cancer | Cisplatin`,
`cervical-cancer | Paclitaxel total CL`, `cervical-cancer | Pembrolizumab`,
and `cervical-cancer | TV-ADC` rows.

## Discoverability fix

A corpus-wide discoverability audit found `C_PAC`, `C_BEV`, and `C_TVADC`
were never written as a single contiguous `double C_<STEM> = <expr>;`
statement anywhere in the file. Each was `$GLOBAL`-forward-declared
(`double C_CIS, C_PAC, C_BEV, C_PEM, C_TVADC;`) and bare-reassigned
*twice* — once in `$ODE`, and again in `$TABLE` (the reassignment already
present there specifically to fix the dose-instant reporting artifact
documented a few lines above it, per this file's own comment). Correct,
working code (not a bug), but not literal-text-discoverable by tooling
that regexes for `double C_<STEM> = ...;`, since neither reassignment
carried the `double` keyword.

**Naive fix (just prepend `double ` to the existing `$TABLE` reassignment)
was tried in a scratch copy first and failed to compile**, hitting the
identical mrgsolve auto-declare collision already diagnosed on
`breast-cancer/bc_mrgsolve_model_refactored.R` and
`celiac-disease/cd_mrgsolve_model_refactored.R` earlier in this batch:
prepending `double` to `C_PAC = CENT_PAC;` (etc.) in `$TABLE`, while the
`$GLOBAL` bare forward-declare still listed the same names, produced
"reference to 'C_PAC' is ambiguous" against qspserver's `mrgsolve_api`
`/model_manifest` — mrgsolve auto-declares a persistent, `$CAPTURE`-
visible class member for every `double NAME = ...;` initializing statement
found anywhere in the block, so a name already forward-declared in
`$GLOBAL` cannot also be (re)declared with `double` elsewhere.

**Fix applied:** removed `C_PAC`, `C_BEV`, `C_TVADC` from the `$GLOBAL`
bare forward-declare (`double C_CIS, C_PAC, C_BEV, C_PEM, C_TVADC;` →
`double C_CIS, C_PEM;`), then prepended `double ` to each of their
existing `$TABLE` reassignments, making that line each name's sole
declaration:

```
double C_PAC = CENT_PAC;
double C_BEV = CENT_BEV;
double C_TVADC = CENT_TVADC;
```

`C_CIS` and `C_PEM` are untouched (out of scope for this fix) — still
forward-declared in `$GLOBAL` and bare-reassigned (no `double`) in both
`$ODE` and `$TABLE`, exactly as before. The `$ODE` bare assignments for
`C_PAC`/`C_BEV`/`C_TVADC` and every downstream `$ODE`-side read
(`VEGF_bind_BEV`, `kill_Pt`, etc.) are untouched — they now simply target
the auto-declared member instead of the manually-declared one, with
identical storage semantics; the `$TABLE`-side recompute (the actual
fix for the dose-instant artifact) is unchanged in intent, only now
literal-text-discoverable as well.

**Verification:** `Rscript -e 'parse(...)'` succeeds with no error (no
straight apostrophes introduced). The DSL string was extracted and posted
to qspserver's `mrgsolve_api` at `localhost:8007`: `/model_manifest`
compiles cleanly and lists `C_PAC`, `C_BEV`, `C_TVADC` in `outputPaths`
alongside `EC50_PAC`/`EC50_BEV`/`EC50_TVADC` in `parameters`.
`/run_simulation` was run for two of the file's own scenarios against both
the pre-fix original DSL (`git show HEAD:...`) and the fixed DSL,
identical dosing:
- **S4-equivalent** (`dose_cisplatin(6, interval_d=21)` + `dose_paclitaxel(6)`
  + `dose_bevacizumab(0, n_doses=20)`, reproduced as three `DoseSpec`
  entries: `CENT_CIS` 4.5333 mg q21d ×6, `CENT_PAC` 46153.85 ng q21d ×6,
  `CENT_BEV` 360.82 mg q21d ×20; 730-day/daily-sampling horizon) — `TV`,
  `C_PAC`, `C_BEV`, `EFFECT_PAC`, `EFFECT_BEV` (734 time points) were
  numerically identical (max abs diff = 0) between the two runs.
- **S5** (`dose_tisotumab(0, n_doses=20)`, reproduced as `CENT_TVADC`
  41.79 mg q21d ×20; same horizon) — `TV`, `C_TVADC`, `EFFECT_TVADC`,
  `MMAE_TVADC` (732 time points) were numerically identical (max abs diff
  = 0) between the two runs.

Grep-confirmed `double C_PAC = `, `double C_BEV = `, `double C_TVADC = `
and `EC50_PAC`, `EC50_BEV`, `EC50_TVADC` all now appear in the file.
