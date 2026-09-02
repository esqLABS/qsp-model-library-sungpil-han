# FCS PK/PD refactor notes

`fcs_mrgsolve_model_refactored.R`, sibling of `fcs_mrgsolve_model.R`. Original
untouched. Every disease-side equation (two-limb TRL clearance, hepatic VLDL,
pancreatitis hazard, endpoints, safety) is copied verbatim; only the five
drugs' PK/PD blocks are reorganized into the fork's pluggable naming
convention.

## Compound identity confirmation (per the guide's warning about
classifier-artifact short codes)

The census (`driver-patches/data/compound_perturbation_census.md`) names
targets, not drugs, in its parenthetical short codes. Checked against the
actual code and `fcs_references.md`/`README.md`:

| Census label | Short code | Confirmed compound | Modality | Real/approved? |
|---|---|---|---|---|
| ANGPTL3 (EVI) | EVI | **Evinacumab** | Anti-ANGPTL3 monoclonal antibody | Yes (approved for HoFH; used off-label/investigationally in severe HTG/FCS) |
| ApoC-III (OLE) | OLE | **Olezarsen** | GalNAc3-conjugated ApoC-III antisense oligonucleotide | Yes (approved for FCS, Balance trial) |
| ApoC-III (VOL) | VOL | **Volanesorsen** | Naked 2'-MOE ApoC-III antisense oligonucleotide | Yes (approved in EU for FCS, APPROACH trial) — genuinely a *different* drug from olezarsen, not a duplicate. The model correctly keeps both as independent treatment arms (`FCS_scenario_volanesorsen` vs `FCS_scenario_olezarsen`), each with its own PK and its own `EC50_<STEM>`. |
| LPL protein (FIB) | FIB | **Fenofibrate** (fenofibric acid) | Oral PPAR-alpha agonist, fibrate class | Yes (approved, generic). The census's target label ("LPL protein") describes what the drug acts on (it induces LPL synthesis/transcription), not what the drug *is*; the compound itself is a real, approved oral small molecule, not a gene-therapy or protein construct. |

**A fifth real compound, plozasiran (PLO), has no census row at all** for
this disease directory. It is a genuine, fully-modeled ApoC-III-directed
GalNAc-siRNA (own SC depot → plasma → hepatic ASGPR uptake → RISC-loading
compartment; dosed in `FCS_scenario_plozasiran()`; used in
`FCS_limb_decomposition`, `FCS_trial_ledger`, `FCS_threshold_time`,
`FCS_galnac_dividend`, `FCS_diet_drug_interaction`, `FCS_scenario_binge`) —
this is a census omission, not a "this compound is out of scope" decision.
It has been refactored alongside the other four for consistency (all three
ApoC-III-directed drugs feed the same `dxdt_APOC3_M` term and belong to the
same naming family) and is reported below as a new census row.

## Archetype classification

| Compound | Archetype | Notes |
|---|---|---|
| Evinacumab (EVI) | **2** — two-compartment linear PK, no depot (dosed directly into `CENT_EVI`) | Not TMDD, despite being a monoclonal antibody: the original models its action as a plain algebraic neutralisation fraction on ANGPTL3 (`ANGFREE = ANGPTL3*(1 - IMAX_EVI*C/(IC50+C))`) — no free-receptor-pool compartment, no drug–receptor complex state, no target-mediated elimination term. A rename, not a refit. |
| Fenofibrate (FIB) | **3**, depot+central variant (no peripheral) | Oral absorption into `GUT_FIB`, single central compartment `CENT_FIB`. Original used a micro-constant (`KE_FIB`, elimination rate constant with no explicit `CL`); per the guide's Archetype 1–3 convention this was rewritten to `CL_FIB`/`V1_FIB` (`CL_FIB = KE_FIB × V_FIB` exactly: 0.0347 × 18 = 0.6246, so `CL_FIB/V1_FIB = 0.0347` reproduces `KE_FIB` bit-for-bit). Pure algebraic rename. |
| Volanesorsen (VOL) | **Bespoke** (depot + central + two one-way tissue sinks, no back-exchange) | `GUT_VOL → CENT_VOL → {LIV_VOL, SYS_VOL}`, each sink with its own first-order elimination. Doesn't fit Archetype 3's peripheral compartment (which exchanges bidirectionally); named `LIV_<STEM>`/`SYS_<STEM>` following the established corpus precedent for this same shape (`elevated-lipoprotein-a`'s `PEL_LIV`). |
| Olezarsen (OLE) | **Bespoke**, same shape as VOL | `GUT_OLE → CENT_OLE → {LIV_OLE, SYS_OLE}`. |
| Plozasiran (PLO) | **Bespoke**, same shape as VOL/OLE plus one extra downstream compartment | `GUT_PLO → CENT_PLO → LIV_PLO → RISC_PLO` (RISC-loading is a genuine extra linear compartment, not a receptor-binding/TMDD system — `KRISC_IN`/`KRISC_OUT` are plain first-order rate constants, no limited receptor pool, no saturable complex formation — so Archetype 4 does not apply here either). |

None of the five compounds required a Hill-curve fit (nls()); every
`EFFECT_<STEM>` below is an exact algebraic rename of the original's own
ratio, with `GAMMA_<STEM> = 1` introduced as an explicit parameter wherever
the original had no Hill exponent (all five).

## A structural nuance the naming convention doesn't fully anticipate:
"concentration" is a compartment amount, not amt/volume, for three compounds

Volanesorsen, olezarsen and plozasiran have **no plasma-volume parameter
anywhere in the original file**. Their disease-facing effect
(`EFF_VOL`/`EFF_OLE`/`EFF_SI` in the original) is computed directly from a
hepatic-tissue (or RISC-loaded) compartment **amount** (mg or AU), never
divided by a volume. `C_VOL`/`C_OLE`/`C_PLO` in the refactor are therefore
identity renames of that amount (`C_VOL = LIV_VOL;`, etc.) — not an invented
concentration, and not a unit conversion the original never had. This is
disclosed in-line in the DSL's own comments and in the `$CAPTURE` annotation
text (units given as mg/mg/AU respectively, not mg/L). Evinacumab and
fenofibrate both have a real plasma volume (`V1_EVI`, `V1_FIB`) so `C_EVI`
and `C_FIB` are true concentrations (mg/L).

## Renaming map (parameters)

| Original | Refactored | Change |
|---|---|---|
| `EMAX_ASO` (shared by VOL & OLE, value 6.0) | `EMAX_VOL = 6.0`, `EMAX_OLE = 6.0` | Split into per-stem copies of the identical original value, so each ApoC-III ASO is independently driveable per the guide's "never collapse several drugs into one shared Hill term" rule. Not a value change. |
| `EMAX_SI` | `EMAX_PLO` | Rename only. |
| `EC50_SI` | `EC50_PLO` | Rename only (needed so the `EC50_<STEM>` discoverability pairing actually matches stem `PLO`). |
| `EC50_VOL`, `EC50_OLE` | unchanged | Already matched the convention. |
| `IFIB_C3` | `EMAX_FIB_APOC3` | Rename only (matches the corpus's established multi-target-compound pattern, e.g. `lpa`'s `EMAX_STA_LDLR`/`EMAX_STA_PCSK9`). |
| `EFIB_LPL` | `EMAX_FIB_LPL` | Rename only. |
| `EC50_FIB` | unchanged | Already matched. |
| `IMAX_EVI` | `EMAX_EVI` | Rename only. |
| `IC50_EVI` | `EC50_EVI` | Rename — needed for the `EC50_<STEM>` discoverability pairing (the original's `IC50_` spelling would not have matched the driver-PK dashboard's regex). |
| `V_FIB` | `V1_FIB` | Rename only. |
| `KE_FIB` | `CL_FIB` | Algebraic rename: `CL_FIB = KE_FIB × V1_FIB` exactly (0.6246 = 0.0347 × 18); `CL_FIB/V1_FIB` reproduces `KE_FIB` bit-for-bit. |
| `V1_EVI`, `V2_EVI`, `CL_EVI`, `Q_EVI` | unchanged | Already matched the convention exactly. |
| new: `GAMMA_VOL`, `GAMMA_OLE`, `GAMMA_PLO`, `GAMMA_EVI`, `GAMMA_FIB` (all = 1) | — | Explicit Hill coefficients, all 1 (original had none, i.e. implicit 1). |
| `APOC3_REF`, `KOUT_M`, `KDEG_P`, `A_CKD_C3`, `KOUT_ANG`, `KOUT_LPL` | unchanged | Disease-side (not drug PK), out of refactor scope. |

## Renaming map (compartments)

| Original | Refactored |
|---|---|
| `VOL_SC`, `VOL_CP`, `VOL_LIV`, `VOL_SYS` | `GUT_VOL`, `CENT_VOL`, `LIV_VOL`, `SYS_VOL` |
| `OLE_SC`, `OLE_CP`, `OLE_LIV`, `OLE_SYS` | `GUT_OLE`, `CENT_OLE`, `LIV_OLE`, `SYS_OLE` |
| `PLO_SC`, `PLO_CP`, `PLO_LIV`, `PLO_RISC` | `GUT_PLO`, `CENT_PLO`, `LIV_PLO`, `RISC_PLO` |
| `EVI_C`, `EVI_P` | `CENT_EVI`, `PERI_EVI` |
| `FIB_GUT`, `FIB_C` | `GUT_FIB`, `CENT_FIB` |

Compartment **order/index is unchanged** (positions 26–41 in both files,
confirmed via `/model_manifest`), so the same numeric `cmt` index doses the
same physical compartment in both the original and the refactored model.

The five dosing helper functions (`volanesorsen()`, `olezarsen()`,
`plozasiran()`, `evinacumab()`, `fenofibrate()`) were updated to dose into
the renamed compartment names; every other R-side function (scenarios 1–10,
analyses 1–8, plotting, `FCS_run_all`) is untouched — none of them reference
a drug PK compartment directly, only the `$CAPTURE`d disease-side columns
(`TG`, `APOC3`, `FLUX1/2/3`, `HAZ_YR`, `PLT`, `ALT`, `XANTH`, …), which keep
their original names throughout.

## Multi-drug combination points (preserved exactly, per the guide's
"combine only at the point of use" rule)

- `dxdt_APOC3_M` sums `EFFECT_VOL + EFFECT_OLE + EFFECT_PLO` — each is
  computed independently from its own compound's own concentration/amount
  and its own `EC50_<STEM>`/`EMAX_<STEM>`; they are combined only at the one
  ODE line the original combined `EFF_VOL + EFF_OLE + EFF_SI`.
- `ASOSYS = SYS_VOL + SYS_OLE` and `ASOLIV = LIV_VOL + LIV_OLE` are a
  genuine shared downstream **safety** liability in the original (systemic/
  hepatic ASO burden driving platelet suppression and ALT rise respectively)
  — not the target-engagement Hill term, which stays per-compound. Renamed,
  not restructured; combined only where the original combined them (the
  `dxdt_PLT`/`dxdt_ALT` lines).

## A real, verification-failing bug found and fixed: the dose-instant
reporting artifact (evinacumab only)

Evinacumab is dosed as a **direct bolus into `CENT_EVI`** (no absorption
depot) — unlike the other four compounds, which are all dosed into an
upstream `GUT_`/depot compartment that builds their `C_<STEM>`-defining
state gradually. A first draft declared `C_EVI`/`EFFECT_EVI` as plain
`$ODE`-local doubles (per the guide's standard template) and captured them
directly. Live testing against the qspserver mrgsolve API found this exactly
reproduces the dose-instant reporting artifact this guide documents: at the
duplicate report row `/run_simulation` emits at the exact instant of a dose,
`C_EVI`/`EFFECT_EVI` read 0 (the pre-dose value) even though `CENT_EVI`
itself already shows the correct post-dose amount (1050 mg after a 15 mg/kg
× 70 kg dose) on that same row. Confirmed via a direct API call:

```
time  C_EVI    EFFECT_EVI  CENT_EVI
236   0        0           0
240   0        0           0      <- duplicate pre-dose row
240   0        0           1050   <- BUG (first draft): C_EVI/EFFECT_EVI still read 0
244   290.4    0.9371      1016.3
```

**Fix applied** (matching the original file's own pre-existing idiom, not
invented): the original already handles exactly this class of problem for
its own disease-side quantities by giving the `$TABLE`-recomputed version of
a quantity a different name from its `$ODE`-local counterpart
(`CL_LPL`→`CLLPL`, `GC3`→`GC3O`, `SATFRAC`→`SATO`, and, tellingly, its own
`CEVI` (ODE) → `CEVI2` (TABLE, captured)). The refactor mirrors this
precisely for evinacumab: the `$ODE`-internal variables that drive live
disease dynamics every solver substep (needed for `ANGFREE`, which feeds
`CL_LPL`/`KLPL` and therefore the TG ODEs directly, and so must stay in
`$ODE` per "keep a calculation in the block the original used it in") are
named `C_EVI_ODE`/`EFFECT_EVI_ODE`; the convention-named `C_EVI`/`EFFECT_EVI`
that get `$CAPTURE`d are recomputed fresh in `$TABLE` from post-dose state.
Confirmed live after the fix — the duplicate rows now read `0` then `300`
(`= 1050/3.5`, the correct post-dose concentration), matching the original's
own `CEVI2` behavior exactly (verification below shows `CEVI2 → C_EVI`
diff = 0 across the whole run).

Volanesorsen, olezarsen, plozasiran and fenofibrate were checked for the
same failure mode and confirmed **not** affected (their `C_<STEM>`-defining
compartment is always downstream of a depot/absorption step, so it never
receives a discontinuous jump at the exact dosing instant) — verified live,
no discontinuity at their respective dose rows. Their `C_<STEM>`/
`EFFECT_<STEM>` stay as single plain `$ODE`-local doubles, captured
directly, per the guide's standard template.

A related compile-time finding, worth recording for future refactors of
this corpus: **`$ODE` and `$TABLE` locals share one C++ scope** in this
mrgsolve build — declaring `double C_EVI = ...;` in both blocks under the
identical name fails with `error: redefinition of 'double {anonymous}::
C_EVI'`. This is presumably *why* the original author used `CEVI`/`CEVI2`
(different names) rather than reusing `CEVI` directly — not a stylistic
choice, a technical requirement. Anyone attempting the `$GLOBAL`-macro fix
this guide recommends elsewhere (`clostridioides-difficile-infection`)
should note that a macro named literally `C_<STEM>` would make the
mandatory discoverability statement `double C_<STEM> = <expr>;` impossible
to write (the macro would expand inside it), so the macro fix and the
discoverability contract are in tension for any file that needs it; the
`_ODE`/`$TABLE`-split approach used here satisfies both simultaneously.

## Verification (qspserver mrgsolve API, http://localhost:8007)

Both `fcs_mrgsolve_model.R`'s and `fcs_mrgsolve_model_refactored.R`'s quoted
DSL blocks were extracted verbatim and confirmed to compile and expose the
expected compartments/parameters via `/model_manifest`
(`Rscript -e 'parse(...)'` also passes for the `.R` file itself — 36
top-level statements in both the original and the refactored file).

`/model_manifest` confirms every discoverability requirement: for every
stem, `EC50_<STEM>` is in `$PARAM` and a literal `double C_<STEM> = <expr>;`
statement exists in the DSL body (`C_EVI`, `C_FIB`, `C_VOL`, `C_OLE`,
`C_PLO` and `EFFECT_EVI`, `EFFECT_VOL`, `EFFECT_OLE`, `EFFECT_PLO`,
`EFFECT_FIB_LPL`, `EFFECT_FIB_APOC3` all appear in `outputPaths`).

Six scenarios were run through both models via `/run_simulation`, identical
dosing (numeric `cmt` index — confirmed via `/model_manifest`'s
`compartments` list that both files order all 43 compartments identically,
so the same index addresses the same physical compartment in both), a
shortened 30-day window (10-day burn + 20-day treatment, `delta = 4 h`) per
this guide's solver-step-budget note, rather than the full 365/180-day
windows the named scenario functions default to:

1. **Natural history** (no drug) — sanity baseline.
2. **Volanesorsen** 300 mg SC weekly (`FCS_scenario_volanesorsen`'s regimen).
3. **Olezarsen** 80 mg SC monthly (`FCS_scenario_olezarsen`'s regimen).
4. **Plozasiran** 25 mg SC q90d (`FCS_scenario_plozasiran`'s regimen).
5. **Evinacumab** 15 mg/kg (70 kg → 1050 mg) q4w (`FCS_scenario_evinacumab`'s
   regimen, `lpl_null` genotype arm).
6. **Fenofibrate** 145 mg PO daily (`FCS_scenario_conventional`'s regimen).

Every `$CAPTURE`d disease-side output (`TG`, `TG_MMOL`, `CM_C`, `VLDL_C`,
`REM_C`, `APOC3`, `CLLPL`, `FLUX1`, `FLUX2`, `FLUX3`, `HAZ_YR`, `PROB_AP`,
`VISCI`, `LACT`, `ABOVE`) and every disease-side ODE state (all 25
non-drug compartments plus `CUMFAT`/`CUMIL6`) matched **exactly** across the
full time grid in all six scenarios (`max_abs_diff = 0`, `max_rel_diff = 0`
for every column, every scenario). Every renamed drug-PK compartment
matched its original counterpart exactly (`VOL_LIV→LIV_VOL`,
`VOL_SYS→SYS_VOL`, `OLE_LIV→LIV_OLE`, `OLE_SYS→SYS_OLE`,
`PLO_RISC→RISC_PLO`, `EVI_C→CENT_EVI`, `EVI_P→PERI_EVI`,
`FIB_C→CENT_FIB`, all `max_abs_diff = 0`), and — after the dose-instant fix
above — `CEVI2→C_EVI` also matched exactly (`max_abs_diff = 0`,
`max_rel_diff = 0`).

This is the expected result for a pure structural reorganization with no
Hill-curve fitting (all five compounds are Archetype 2/3/bespoke-linear, no
TMDD, no `nls()` fit anywhere): the guide's tolerance bar for this class of
refactor is "near-exact match... anything beyond floating-point-scale
deviation means a bug", and the achieved result is an exact 0, not merely
within tolerance.

## No upstream build defects found

`fcs_mrgsolve_model.R`'s DSL block compiled cleanly under mrgsolve 2.0.1 via
the qspserver API with no changes needed — no entry added to
`translations/UPSTREAM_ISSUES.md`.

## qspserver infrastructure note

One transient `all_loaded(x) is not TRUE` / `script_error` failure was hit
on `/model_manifest` and `/run_simulation` mid-session (consistent with the
concurrent-load cache contention `HANDOFF.md` describes other agents having
hit this session). It self-recovered on retry a few seconds later with no
intervention — unlike the persistent stale-`mrgmod_cache.RDS` failure mode
`HANDOFF.md` describes (which needed a `docker exec ... rm` fix from the
orchestrating session), this one cleared on its own. Reported here for
visibility; no action was needed and none was taken.

## Census update

`driver-patches/data/compound_perturbation_census.md` rows for
`familial-chylomicronemia-syndrome` updated with confirmed compound
identity, target/pathway and outcome; a new row added for plozasiran (PLO),
which the census had missed entirely.
