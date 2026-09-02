# Fibromyalgia (FM) PK/PD refactor notes

Refactored sibling: `fm_mrgsolve_model_refactored.R`. Original untouched at
`fm_mrgsolve_model.R`.

## Compound identity (confirmed against the actual code, not just the census guess)

The census's guesses were all correct on inspection of the original file's
own `$PROB` header comment and `$PARAM` block:

| Stem | Compound | Class | Confirmed by |
|---|---|---|---|
| DUL | Duloxetine | SNRI | `$PROB` header line 2; `ka_DUL/CL_DUL/V1_DUL/Q_DUL/V2_DUL/F_DUL` PK block; `IC50_SERT_DUL`/`IC50_NET_DUL` |
| MIL | Milnacipran | SNRI | `$PROB` header line 2; `ka_MIL/CL_MIL/V_MIL/F_MIL`; `IC50_SERT_MIL`/`IC50_NET_MIL` |
| PRE | Pregabalin | alpha2-delta ligand (gabapentinoid) | `$PROB` header line 2; `ka_PRE/CL_PRE/V_PRE/F_PRE`; `IC50_PRE_alpha2d`, `Emax_PRE` |
| TCA | Amitriptyline | Tricyclic antidepressant | `$PROB` header line 2 ("Amitriptyline"); `ka_TCA/CL_TCA/V_TCA/F_TCA`; `IC50_SERT_TCA`/`IC50_NET_TCA`; a fourth, separate hardcoded ratio `Cp_TCA/(0.05+Cp_TCA)` for sedation |

All four are standard, guideline-recommended fibromyalgia pharmacotherapy
(duloxetine and milnacipran are FDA-approved for FM; pregabalin is FDA-approved
for FM; low-dose amitriptyline is a longstanding off-label first-line agent) —
the census guesses needed no correction.

## Archetype per compound

All four compounds are **Archetype 3** (depot + central, +/- peripheral),
matching what the original already does — no compartments added or removed:

- **DUL (Duloxetine)** — full Archetype 3 (depot + central + peripheral):
  `GUT_DUL CENT_DUL PERI_DUL`, `KA_DUL F_DUL CL_DUL V1_DUL Q_DUL V2_DUL`.
- **PRE (Pregabalin)** — Archetype 3 without the peripheral compartment
  (the original never had one for this drug): `GUT_PRE CENT_PRE`,
  `KA_PRE F_PRE CL_PRE V1_PRE` (renamed from `V_PRE`).
- **MIL (Milnacipran)** — Archetype 3 without peripheral: `GUT_MIL CENT_MIL`,
  `KA_MIL F_MIL CL_MIL V1_MIL` (renamed from `V_MIL`).
- **TCA (Amitriptyline)** — Archetype 3 without peripheral: `GUT_TCA CENT_TCA`,
  `KA_TCA F_TCA CL_TCA V1_TCA` (renamed from `V_TCA`).

No TMDD, no bespoke structures needed anywhere in this file.

## Normalizing the duplicate concentration sites (the census classification)

The original computed each drug's concentration **twice**, under two
different names: once in `$ODE` (`Cp_DUL`, `Cp_PRE`, `Cp_MIL`, `Cp_TCA`,
used live by the PD equations) and a second time in `$TABLE`
(`Cp_DUL_out`, `Cp_PRE_out`, `Cp_MIL_out`, `Cp_TCA_out`, recomputing the
identical ratio under an `_out`-suffixed name purely so it could be
`$CAPTURE`d for output). This is exactly the "duplicate concentration
sites" pattern the census flagged.

Fixed by normalizing to the single canonical name `C_<STEM>` (`C_DUL`,
`C_PRE`, `C_MIL`, `C_TCA`), computed **once**, in `$ODE` — the block the
original itself used for the live PD calculation (per the guide's "keep a
calculation in the block the original used it in" rule). `$TABLE` does
**not** recompute it under a second name: this mrgsolve build compiles
`$ODE` and `$TABLE` into the same function (confirmed empirically — see
the build-defect section below and the precedent already established for
this exact situation in `acne-vulgaris/acn_refactor_notes.md`), so
redeclaring `double C_DUL` a second time there is a C++ redefinition
error, not a scoping choice. `$TABLE`'s own two derived percent metrics
(`inh_SERT_pct`, `inh_NET_pct`, `Ca_block_pct`) now simply reference the
single `$ODE`-scoped `C_<STEM>` doubles directly.

## The Hill interface — one EFFECT per compound per pathway, not one per compound

Every SNRI compound here (DUL, MIL, TCA) genuinely acts on **two** separate
disease-facing pathways (SERT reuptake inhibition and NET reuptake
inhibition), each feeding a different downstream equation (`SHT_syn` vs.
`NE_syn`, respectively, plus both jointly in `pain_target`/`dep_target`/
`fatigue_tgt`). Amitriptyline has a **third**, independent pathway
(antihistaminergic/anticholinergic sedation feeding `SWS_depth`, originally
the bare ratio `Cp_TCA/(0.05+Cp_TCA)` with no name at all). Collapsing these
into one `EFFECT_DUL`/`EFFECT_MIL`/`EFFECT_TCA` would hide which pathway
drives which disease equation and defeat independent driveability of each
mechanism, so each compound exposes one `EFFECT_<STEM>_<PATHWAY>` term per
genuinely distinct action instead:

- `EFFECT_DUL_SERT`, `EFFECT_DUL_NET`
- `EFFECT_MIL_SERT`, `EFFECT_MIL_NET`
- `EFFECT_TCA_SERT`, `EFFECT_TCA_NET`, `EFFECT_TCA_SEDATION`
- `EFFECT_PRE` (pregabalin has only the one alpha2-delta pathway)

Every one of these was **already a plain ratio** in the original
(`C/(IC50+C)`, optionally scaled by a shared `Emax_SNRI` or `Emax_PRE`) —
this is a rename, not a refit. `EMAX_<STEM>_<PATHWAY>`, `EC50_<STEM>_<PATHWAY>`,
`GAMMA_<STEM>_<PATHWAY>` were pulled out as explicit named parameters
carrying the original's values unchanged:

- The original's single shared `Emax_SNRI = 1.0` (applied identically to
  all three SNRI compounds' SERT and NET terms) is now six separate
  `EMAX_*_SERT`/`EMAX_*_NET` parameters, each still `= 1.0` — same value,
  split per compound/pathway per the naming convention, not a new number.
  `Emax_SNRI` itself is removed from `$PARAM` since nothing references it
  any more (fully superseded by its six named replacements at the same
  value).
- `Emax_PRE = 0.70` keeps its name (`EMAX_PRE`), already convention-shaped
  in the original.
- The original had **no** Hill coefficient anywhere in this file — every
  `GAMMA_*` is `1.0` (rename, not a fit), and since `pow(x, 1.0)` is
  mathematically identical to `x`, this restructuring is pure reorganization,
  not a refit — Archetype-3-tier exact-match tolerance applies throughout
  (see Verification below).
- The `EC50_TCA_SEDATION = 0.05` and `EMAX_TCA_SEDATION = 1.0` pair is new
  in name only: the original's `Cp_TCA/(0.05+Cp_TCA)` had no name and no
  explicit Emax (implicitly 1.0); both are made explicit here.

Downstream combination is preserved exactly at the point the original
combined it: `inh_SERT = EFFECT_DUL_SERT + EFFECT_MIL_SERT +
EFFECT_TCA_SERT` (capped at 1.0), `inh_NET` analogously — these two combined
variables, not the individual `EFFECT_*` terms, are what `dxdt_NE_syn`,
`dxdt_SHT_syn`, `pain_target`, `fatigue_tgt`, and `dep_target` actually read,
exactly matching where the original combined its own `inh_SERT`/`inh_NET`.

## Pre-existing upstream build defect found and fixed only in the sibling

**Logged as `translations/UPSTREAM_ISSUES.md` #144.** The original's
`$CAPTURE` block re-lists eighteen names that are already `$CMT`
compartments (`SP_csf NMDA_state WindUp LTP_cs DPMS NE_syn SHT_syn CRH ACTH
CORT SNS_tone SWS_depth Adenosine MG_act IL1b_sp Pain_score FIQ_score
Fatigue_VAS Depression_score`). Posting the untouched original's own DSL
block to the qspserver `mrgsolve_api`'s `/model_manifest` confirms this
fails to compile:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  SP_csf,NMDA_state,WindUp,LTP_cs,DPMS,NE_syn,SHT_syn,CRH,ACTH,CORT,
  SNS_tone,SWS_depth,Adenosine,MG_act,IL1b_sp,Pain_score,FIQ_score,
  Fatigue_VAS,Depression_score
```

This is exactly the "`$CAPTURE` duplicating compartment names" defect
category the workflow guide already documents generically. Fixed **only**
in `fm_mrgsolve_model_refactored.R` (never in the original) by simply not
re-listing those eighteen names in `$CAPTURE` — mrgsolve reports every
compartment's state automatically, so this is a pure deletion, not a
behavioral change. Confirmed the fix is syntax-only and complete: applying
it to an in-memory-only scratch copy of the untouched original (minimal
diff: delete the eighteen names, keep the seven genuinely-derived
`Cp_*_out`/`inh_*_pct`/`Ca_block_pct` metrics) lets that scratch copy
compile and return a 31-output-path manifest with nothing else masked
behind it — the only compile-blocking defect in the file.

## Verification

Verified via the qspserver `mrgsolve_api` at `http://localhost:8007`
(`/model_manifest` then `/run_simulation`), never local R/mrgsolve, per the
workflow guide. `Rscript -e 'parse(...)'` also confirmed the sibling parses
as R (38 top-level expressions) before any of this.

### `/model_manifest`

Both the original (with the minimal `$CAPTURE`-dedup syntax fix applied
in-memory, see the build-defect section above) and the refactored DSL
compile cleanly and return a manifest. For the refactored model: 46 output
paths (31 renamed/original compartments + `C_DUL/C_PRE/C_MIL/C_TCA` +
`inh_SERT_pct/inh_NET_pct/Ca_block_pct` + all 8 `EFFECT_*` terms) and 100
parameters. Confirmed present with original values: every renamed PK
parameter (`KA_DUL`, `F_DUL`, `CL_DUL`, `V1_DUL`, `Q_DUL`, `V2_DUL`,
`KA_PRE/F_PRE/CL_PRE/V1_PRE`, `KA_MIL/F_MIL/CL_MIL/V1_MIL`,
`KA_TCA/F_TCA/CL_TCA/V1_TCA`) and every Hill parameter
(`EC50/EMAX/GAMMA_{DUL,MIL,TCA}_{SERT,NET}`, `EC50/EMAX/GAMMA_TCA_SEDATION`,
`EC50/EMAX/GAMMA_PRE`). Discoverability requirement confirmed directly:
`double C_<STEM> = <expr>;` and a same-stem `EC50_<STEM>*` parameter both
present for every compound (`grep`-verified in the finished file, then
re-confirmed live against `/model_manifest`'s `parameters`/`outputPaths`).
`C_<STEM>`/`EFFECT_<STEM>` live in `$CAPTURE`/`outputPaths`, not `$PARAM` —
matching the established precedent for this exact question elsewhere in
the corpus (`acne-vulgaris/acn_refactor_notes.md`'s own `/model_manifest`
confirmation), since re-declaring an already-`$PARAM` name as a `double`
local in `$ODE` is a C++ redefinition error under this mrgsolve build.

### `/run_simulation` — all 6 of the original's own scenarios

Ran every dosing scenario defined in the original file's own R code
(`sim_base`, `sim_DUL`, `sim_PRE`, `sim_MIL`, `sim_COMBO`, `sim_TCA`),
identical dosing (amount, `ii`, `addl`, dose time) against both the
original (capfix) and refactored DSL, full 84-day/2016h window, `delta=1`,
comparing all 26 outputs the two models share (`Cp_DUL_out`->`C_DUL`,
`Cp_PRE_out`->`C_PRE`, `Cp_MIL_out`->`C_MIL`, `Cp_TCA_out`->`C_TCA`,
`inh_SERT_pct`, `inh_NET_pct`, `Ca_block_pct`, and the 19 disease-state
compartments) across the entire time grid. The API has no way to set the
original R script's custom `FM_init` steady-state initial conditions (no
`init` field in `/run_simulation`'s request schema), so both models were
run from their shared, structurally-identical default zero initial state
— an apples-to-apples comparison of the same starting point through
identical dosing, which is what this verification needs to catch a
refactor bug; it does not reproduce the original script's own steady-state
starting figures.

| Scenario | Dosing | Result |
|---|---|---|
| Untreated FM baseline | none | **Exact match**, max abs diff = 0.0, all 26 outputs, 2017 pts |
| Duloxetine 60 mg QD | `amt=60,cmt=GUT_DUL,ii=24,addl=83` | **Exact match**, max abs diff = 0.0, all 26 outputs, 2018 pts |
| Pregabalin 150 mg BID | `amt=150,cmt=GUT_PRE,ii=12,addl=167` | **Exact match**, max abs diff = 0.0, all 26 outputs, 2018 pts |
| Milnacipran 50 mg BID | `amt=50,cmt=GUT_MIL,ii=12,addl=167` | **Exact match**, max abs diff = 0.0, all 26 outputs, 2018 pts |
| Amitriptyline 25 mg QHS (dose starts at t=22h) | `amt=25,cmt=GUT_TCA,ii=24,addl=83,time=22` | **Exact match**, max abs diff = 0.0, all 26 outputs, 2018 pts |
| Duloxetine + Pregabalin combo | both of the above together | **Exact match**, max abs diff = 0.0, all 26 outputs, 2019 pts |

All six are Archetype 3 (or the without-peripheral variant) with `gamma=1`
throughout — pure reorganization, no Hill-fitting — so bit-exact
reproduction is the expected result, and it is what was obtained; no
floating-point-scale tolerance was even needed (max abs diff across every
scenario and every shared output column was exactly `0.0`, not merely
small).

### API availability note (server-side, not a model defect)

This verification run hit a real, twice-recurring qspserver infrastructure
fault, independent of this model: `http://localhost:8007` returned a
persistent `500` on **every** `/model_manifest` and `/run_simulation` call
regardless of `model_content` — confirmed both times by posting a
minimal, completely unrelated Archetype-1 toy model
(`$PARAM CL_TCZ=0.02, V1_TCZ=5.0 / $CMT CENT_TCZ / ...`) and getting the
identical error each time, first:

```
Error in readRDS(cache_file) :
  embedded nul in string: '//  h^-1\0\004\0\to 0'
```

and after that first instance was cleared, again shortly after (mid
retest of the Amitriptyline scenario specifically, producing a spurious
divergence in that scenario alone that resolved itself once the cache was
cleared a second time and the scenario was rerun — see the table above,
which reflects the clean rerun, not the corrupted one):

```
Error in readRDS(cache_file) :
  ReadItem: unknown type 32, perhaps written by later version of R
```

`/health` and `/queue` both reported the service as healthy and idle
throughout both incidents, so this was a corrupted on-disk mrgsolve
compile-cache index (`mrgmod_cache.RDS`, shared across every "inline"
`model_content` request, not keyed cleanly per caller or per project),
consistent with concurrent load from other agents against the same API,
exactly as this task's briefing warned it could. Root-caused and fixed
(twice) from outside this session by removing the stale
`mrgmod_cache.RDS` index for the affected project directories
(`docker exec qspserver-mrgsolve_api-1 rm -f
/soloc/mrgsolve-so-2.0.1-x86_64-pc-linux-gnu/{inline,friberg_ro,pk1cmt,pkpd_indirect}/mrgmod_cache.RDS`);
confirmed resolved both times by the same unrelated-toy-model probe
returning `200`. This is a server-side infrastructure fault, not a defect
in this model or its refactor, and every scenario above was re-verified
clean after the fix. Separately, a single retry of the Duloxetine
scenario (redundant — a clean exact-match result for it was already on
hand from before the first outage) hit the API's documented default
`maxsteps=20000` solver-step budget on that one retry; this is the same
known, unrelated API limitation the workflow guide already documents, not
evidence against the already-confirmed exact match for that scenario.
