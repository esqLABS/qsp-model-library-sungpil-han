# Refactor notes — `benign-prostatic-hyperplasia/bph_mrgsolve_model.R`

All four compounds in the file were refactored: **Dutasteride**,
**Finasteride**, **Tadalafil**, **Tamsulosin** — the four rows in
`driver-patches/data/compound_perturbation_census.md` classified
`benign-prostatic-hyperplasia | Dutasteride/Finasteride/Tadalafil/
Tamsulosin`, all "Redirect concentration (clean single site)". There are
no other compounds in this file, so nothing was left untouched at the
compound level (the disease-side PD — DHT/AR axis, prostate volume, IPSS/
Qmax/PVR/PSA — and the multi-drug combination arithmetic are unchanged,
only renamed and reorganized where a compound's own PK/effect terms
touch them).

## Archetype

**Archetype 3 (depot + central + peripheral, linear), for all four
compounds.** Each was already exactly this shape in the original: a
3-compartment oral PK model (`<STEM>_GUT`/`<STEM>_C`/`<STEM>_P`, CL/Q/V1/
V2/KA/F), no TMDD, no receptor-binding ODEs. Nothing was added or removed;
only renamed to the guide's convention.

| Role | Tamsulosin | Finasteride | Dutasteride | Tadalafil |
|---|---|---|---|---|
| Depot | `TAMS_GUT`→`GUT_TAMS` | `FINA_GUT`→`GUT_FINA` | `DUT_GUT`→`GUT_DUT` | `TAD_GUT`→`GUT_TAD` |
| Central | `TAMS_C`→`CENT_TAMS` | `FINA_C`→`CENT_FINA` | `DUT_C`→`CENT_DUT` | `TAD_C`→`CENT_TAD` |
| Peripheral | `TAMS_P`→`PERI_TAMS` | `FINA_P`→`PERI_FINA` | `DUT_P`→`PERI_DUT` | `TAD_P`→`PERI_TAD` |
| Absorption | `ka_tams`→`KA_TAMS` | `ka_fina`→`KA_FINA` | `ka_dut`→`KA_DUT` | `ka_tad`→`KA_TAD` |
| Bioavailability | `F_tams`→`F_TAMS` | `F_fina`→`F_FINA` | `F_dut`→`F_DUT` | `F_tad`→`F_TAD` |
| Clearance | `CL_tams`→`CL_TAMS` | `CL_fina`→`CL_FINA` | `CL_dut`→`CL_DUT` | `CL_tad`→`CL_TAD` |
| Central vol. | `V1_tams`→`V1_TAMS` | `V1_fina`→`V1_FINA` | `V1_dut`→`V1_DUT` | `V1_tad`→`V1_TAD` |
| Peripheral vol. | `V2_tams`→`V2_TAMS` | `V2_fina`→`V2_FINA` | `V2_dut`→`V2_DUT` | `V2_tad`→`V2_TAD` |
| Inter-cmt CL | `Q_tams`→`Q_TAMS` | `Q_fina`→`Q_FINA` | `Q_dut`→`Q_DUT` | `Q_tad`→`Q_TAD` |
| Exposed conc. | `CP_tams`→`C_TAMS` | `CP_fina`→`C_FINA` | `CP_dut`→`C_DUT` | `CP_tad`→`C_TAD` |
| EC50 | `IC50_tams`→`EC50_TAMS` | `IC50_fina`→`EC50_FINA` | `IC50_dut`→`EC50_DUT` | `IC50_tad`→`EC50_TAD` |
| Emax | `Emax_tams`→`EMAX_TAMS` (unchanged value, 0.95) | new `EMAX_FINA=1.0` | new `EMAX_DUT=1.0` | new `EMAX_TAD=1.0` |
| Gamma | new `GAMMA_TAMS=1.0` | new `GAMMA_FINA=1.0` | new `GAMMA_DUT=1.0` | new `GAMMA_TAD=1.0` |
| Compound effect | `occ_tams`→`EFFECT_TAMS` | `inh_fina`→`EFFECT_FINA` | `inh_dut`→`EFFECT_DUT` | `inh_pde5`→`EFFECT_TAD` |

All PK parameter *values* are copied verbatim from the original — nothing
invented or defaulted. `USE_TAMS`/`USE_FINA`/`USE_DUT`/`USE_TAD` (the
per-scenario on/off gating flags) are untouched, same names, same values.

## Hill interface: rename, not a fit — for all four

Every compound's effect term in the original was already a bare Emax/Imax
ratio, read directly by the downstream disease equations, with no further
transform on top:

- Tamsulosin: `occ_tams = Emax_tams * CP_tams / (IC50_tams + CP_tams)`
  (gated by `USE_TAMS > 0.5`) — already has an explicit `Emax_tams=0.95`.
- Finasteride: `inh_fina = CP_fina / (IC50_fina + CP_fina)` (gated by
  `USE_FINA`) — implicit `Emax=1`.
- Dutasteride: `inh_dut = CP_dut / (IC50_dut + CP_dut)` (gated by
  `USE_DUT`) — implicit `Emax=1`.
- Tadalafil: `inh_pde5 = CP_tad / (IC50_tad + CP_tad)` (gated by
  `USE_TAD`) — implicit `Emax=1`.

All four are exactly the guide's Hill-interface shape with `gamma=1`. Per
the guide ("If the original's effect term is already this shape... this
is a rename, not a refit"), each was rewritten to the explicit
`EMAX*pow(C,GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))` form with named
`EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>=1` parameters — Tamsulosin's
`EMAX_TAMS=0.95` is the original's own explicit value, renamed;
Finasteride/Dutasteride/Tadalafil's `EMAX_<STEM>=1.0` and all four
`GAMMA_<STEM>=1.0` are new named parameters making the previously-implicit
values explicit, per the guide's own instruction for this exact case
("`GAMMA_<STEM> = 1` if the original had no explicit Hill coefficient").
No curve fit was performed or needed; the resulting arithmetic is
identical to the original's bare ratio (`pow(x, 1) == x`).

## The multi-drug combination arithmetic (untouched, only its inputs renamed)

```
double inh_SRD5A2 = EFFECT_FINA * (1.0 - EFFECT_DUT) + EFFECT_DUT;  // was inh_fina/inh_dut
double inh_SRD5A1 = EFFECT_DUT;
double total_5AR_inh = f2_fraction * inh_SRD5A2 + f1_fraction * inh_SRD5A1;
...
double pde5_activity = 1.0 - EFFECT_TAD;  // was inh_pde5
...
dxdt_ALPHA1_OCC = k_alpha1_eq * (EFFECT_TAMS - ALPHA1_OCC);  // was occ_tams
```

`f2_fraction`/`f1_fraction` (SRD5A2/SRD5A1 tissue contribution split),
`total_5AR_inh`'s downstream use in the testosterone/DHT ODEs, and
`pde5_activity`'s use in the cGMP ODE are all byte-identical to the
original — only the four compound-effect inputs feeding these expressions
were renamed. Per the guide ("keep each compound's `EFFECT_<STEM>`
separate; combine them only at the point the disease equations actually
use them"), no drugs were collapsed into a shared term.

## `C_<STEM>`/`EFFECT_<STEM>` are `$ODE`-local (+ independently recomputed
`$TABLE`-local), not `$PARAM`

Same confirmed mrgsolve 2.0.1 constraint already documented in several
other `*_refactor_notes.md` in this repo (`autoimmune-polyendocrinopathy`,
`chronic-lymphocytic-leukemia`, the original `rheumatoid-arthritis`
tocilizumab calibration run): a name cannot simultaneously be a settable
`$PARAM` default and a value freshly assigned every `$ODE`/`$TABLE` step
in this build ("assignment of read-only reference"). All eight
(`C_TAMS`/`C_FINA`/`C_DUT`/`C_TAD`, `EFFECT_TAMS`/`EFFECT_FINA`/
`EFFECT_DUT`/`EFFECT_TAD`) are computed once in `$ODE` (bare assignment,
no `double` — see the build-defect section below for why) for use by the
disease equations, recomputed independently in `$TABLE` from raw
compartment state (this file's own pre-existing convention — every other
`$TABLE` quantity, e.g. `PV_change_pct`, is likewise re-derived from state
rather than re-using an `$ODE`/`$MAIN` local), and exposed via
`capture NAME = expr;`. Confirmed visible via `POST /run_simulation`'s
`outputs` selection and via `POST /model_manifest`'s `outputPaths` (48
entries, including all eight); confirmed **not** listed in
`/model_manifest`'s `parameters` (73 entries: every `EC50_<STEM>`/
`EMAX_<STEM>`/`GAMMA_<STEM>`/`CL_<STEM>`/`V1_<STEM>`/`V2_<STEM>`/
`Q_<STEM>`/`KA_<STEM>`/`F_<STEM>`/`USE_<STEM>` is present, none of the
eight `C_`/`EFFECT_` names). A future driver-patch redirect of, say,
`C_TAMS` will need to replace the `C_TAMS = ...;` line itself, not a
`param()`-style override.

## Two pre-existing build defects found (both fixed, syntax-only, in the delivered `_refactored.R`; logged as `translations/UPSTREAM_ISSUES.md` #59)

Neither is related to any of the four compounds' own PK/PD — both block
the file from compiling at all, original or refactored, until worked
around.

1. **The deprecated `_init_<CMT>` idiom.** `$MAIN`'s
   `if(NEWIND <= 1) { ... }` block set all 24 compartments' initial
   conditions with `_init_<CMT> = value;` (e.g. `_init_TAMS_GUT = 0;`,
   `_init_TEST_P = TEST0 * 50.0;`). mrgsolve 2.0.1 does not accept this
   symbol at all (`'_init_TAMS_GUT' was not declared in this scope`). Same
   defect class as issue #29 (`polymyalgia-rheumatica`). Fix: converted
   every line to the modern `<CMT>_0 = value;` idiom (24 lines) — no
   compartment, index, or starting value changed.
2. **`$ODE`/`$TABLE`-local `double`s colliding with same-named `capture`
   statements**, two separate instances:
   - The four original plasma-concentration locals (`CP_tams`, `CP_fina`,
     `CP_dut`, `CP_tad`) were declared `double NAME = ...;` in `$ODE` for
     that compound's own effect calculation, while `$TABLE` *also* used
     the identical bare name as an old-style
     `capture NAME = NAME_out;` target. mrgsolve 2.0.1 auto-promotes any
     `$CAPTURE`d name to a class member, colliding with the `$ODE`
     `double` of the same name (`redefinition of 'capture
     {anonymous}::CP_tams'`). This refactor's own new `C_TAMS`/`C_FINA`/
     `C_DUT`/`C_TAD` (and `EFFECT_TAMS`/`EFFECT_FINA`/`EFFECT_DUT`/
     `EFFECT_TAD`) hit the identical collision against their own new
     `$TABLE` captures. Fix (applied to both the four original names, in
     an in-memory-only verification copy, and to this refactor's own
     eight new names, in the delivered file): drop `double` from the
     `$ODE` declaration — a bare assignment to the member mrgsolve
     already promotes from the later `capture NAME = ...;` line (mrgsolve
     scans the whole model source for capture statements before compiling
     each block, so file order does not matter); confirmed empirically
     via `POST /model_manifest`.
   - Eight pre-existing `$TABLE`-only collisions (`DHT_inhibition_pct`,
     `DHT_PROST_inh_pct`, `PV_change_pct`, `IPSS_change`, `QMAX_change`,
     `PVR_change`, `PSA_change_pct`, `Alpha1_block_pct`), each computed as
     `double NAME = expr;` and then re-exposed as `capture NAME = NAME;`
     in the same block. Same defect class as issues #35 (`dengue`), #43
     (`disseminated-intravascular-coagulation`), #50 (`alcoholic-liver-
     disease`). Fix: renamed just the intermediate `double` local (append
     `_calc`); the exposed capture name and its arithmetic are unchanged.

Both fixes are syntax-only and change no numeric value; applied directly
in the delivered `bph_mrgsolve_model_refactored.R` (not just a throwaway
scratch copy), per the guide's settled policy. The checked-in original
(`bph_mrgsolve_model.R`) is untouched and still carries both defects
exactly as written. The same fixes (with the original's own
`CP_tams`/`CP_fina`/`CP_dut`/`CP_tad` names for the `$ODE`-vs-`$TABLE`
collision) were applied to an in-memory-only scratch copy of the original
purely to reach a buildable state for the `/run_simulation` comparison
below — never written back to the tracked original.

## Verification

Per the guide's mandatory protocol: ran the original file's own six
treatment scenarios through both the (build-defect-patched, for
buildability only) original and the refactored model via the qspserver
`mrgsolve_api` service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`), full 730-day duration with daily (`delta=24`)
output — matching the original R driver's own `sim_duration_days=730`,
`delta=24` exactly. No shortening was needed; the full 732-point (733 for
the two-drug combination scenario, which adds one extra distinct dose
time) grid ran to completion for every scenario, well within the API's
default solver step budget.

Scenarios run (dose amounts, route, and once-daily interval copied exactly
from the original file's own `scenario_list`, translated to the API's
`dosing`/`ii`/`addl` convention):

1. **Watchful Waiting** — no dosing (all `USE_*` flags 0).
2. **Tamsulosin 0.4 mg QD** — 400 µg into `GUT_TAMS` (cmt 1), q24h ×730.
3. **Finasteride 5 mg QD** — 5000 µg into `GUT_FINA` (cmt 4), q24h ×730.
4. **Dutasteride 0.5 mg QD** — 500 µg into `GUT_DUT` (cmt 7), q24h ×730.
5. **Combination (Dutasteride + Tamsulosin)** — both of the above doses
   together.
6. **Tadalafil 5 mg QD** — 5000 µg into `GUT_TAD` (cmt 10), q24h ×730.

Every `$CAPTURE`d output was compared point-by-point across the full time
grid for all six scenarios: all 24 compartment states (`GUT_TAMS`/
`CENT_TAMS`/`PERI_TAMS`/... and the 12 disease-side compartments,
unchanged names), the four exposed concentrations (`C_TAMS`/`C_FINA`/
`C_DUT`/`C_TAD`, renamed from `CP_tams`/`CP_fina`/`CP_dut`/`CP_tad`), the
four Hill effect terms (`EFFECT_TAMS`/`EFFECT_FINA`/`EFFECT_DUT`/
`EFFECT_TAD`, renamed from `occ_tams`/`inh_fina`/`inh_dut`/`inh_pde5`),
and every disease-side derived output (`DHT_inhibition_pct`,
`DHT_PROST_inh_pct`, `PV_change_pct`, `PV_cc`, `IPSS_score`,
`IPSS_change`, `QMAX_mL_s`, `QMAX_change`, `PVR_mL`, `PVR_change`,
`PSA_ngmL`, `PSA_change_pct`, `Alpha1_block_pct`, `AR_act_norm`,
`cGMP_norm`, `Inflam_idx`).

**Result: exact match, max abs diff = 0.0 for every output, every
scenario.** This is the expected outcome for a pure rename (Archetype 3,
no Hill-fitting) per the guide's own tolerance rule — since every
compound's PK is a straight compartment/parameter rename and every
effect term is an exact rename of an already-Hill-shaped ratio (`gamma=1`
throughout, `pow(x,1)==x`), there is no approximation to introduce error.
Sanity-checked directly as well: at steady state under Tamsulosin
monotherapy, `EFFECT_TAMS` (≈0.7996–0.7998, e.g. day ~10–15) tracks
`Alpha1_block_pct/100` (≈0.8003–0.8005) with the small, expected lag from
`ALPHA1_OCC`'s own first-order equilibration (`k_alpha1_eq=5/h` against
the instantaneous `EFFECT_TAMS` target) — consistent with the model's own
design, not an artifact of the refactor.

## Anything else worth flagging

- The `.cpp` DSL extraction used for every verification call above was
  confirmed byte-identical to the quoted `code <- '...'` block inside the
  delivered `bph_mrgsolve_model_refactored.R` (re-extracted and
  string-compared after the file was finalized) — the DSL was never
  hand-edited after being tested. No separate `.cpp` file is kept as a
  deliverable; the extraction was a verification-only step.
- Compartment ordering is unchanged (`GUT_TAMS` is still compartment 1,
  `GUT_FINA` compartment 4, `GUT_DUT` compartment 7, `GUT_TAD`
  compartment 10, etc.) since renaming was done in place without
  reordering — any external code addressing compartments by 1-based
  number rather than name is unaffected.
- The API's documented step-count limitation (`maxsteps` default well
  below what some models' own drivers request) was not encountered here;
  this model's dynamics are not stiff enough at the 730-day/daily-output
  scenario scale to hit it, so no scenario needed shortening.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`benign-prostatic-hyperplasia | Dutasteride`, `| Finasteride`,
`| Tadalafil`, and `| Tamsulosin` rows.
