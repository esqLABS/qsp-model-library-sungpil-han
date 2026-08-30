# Refactor notes — cytokine-release-syndrome / five compounds (ANA, DEX, RUX, SILT, TOCI)

Covers `cytokine-release-syndrome/crs_mrgsolve_model.R`. All five compounds
this file models were in scope (Anakinra, Dexamethasone, Ruxolitinib,
Siltuximab, Tocilizumab per their existing rows in
`driver-patches/data/compound_perturbation_census.md`) — there are no other
compounds in this file to leave untouched. Every disease-side compartment and
equation (CAR-T dynamics, T-cell/macrophage-derived cytokines, IL-6 downstream
signaling, endothelial/vascular response, CRS severity, organ damage) is
copied verbatim from the original into `crs_mrgsolve_model_refactored.R`;
only the five drugs' own PK blocks and their named effect terms were renamed
to the fork's convention.

## Archetype determined per compound

The task note flagged that Tocilizumab's archetype should be checked against
this file's own equations rather than assumed from a prior refactor of the
same drug elsewhere in the corpus (e.g. the TMDD archetype used in
`rheumatoid-arthritis/ra_refactor_notes.md`, or Archetype 1 in
`sepsis/sep_refactor_notes.md`). Checked directly: **none of the five
compounds in this file involve TMDD or any receptor-binding ODE system** —
all five are plain linear PK with an algebraic Emax effect term.

| Compound | Archetype | Structure |
|---|---|---|
| Tocilizumab (`TOCI`) | 2 | no depot, 2-cmt (central+peripheral), linear CL/Q/V |
| Siltuximab (`SILT`) | 1 | no depot, 1-cmt, linear elimination |
| Dexamethasone (`DEX`) | 2 | no depot, 2-cmt (central+peripheral), linear CL/Q/V |
| Ruxolitinib (`RUX`) | 3 | depot + central + peripheral, linear |
| Anakinra (`ANA`) | 3 minus peripheral | depot + central, linear |

Siltuximab was specifically checked for TMDD (anti-IL-6 antibody) and
Tocilizumab for TMDD (anti-IL-6 receptor antibody), per the task's explicit
instruction — **neither is TMDD in this file**: both are a single
concentration compartment feeding a plain `Emax*C/(EC50+C)` ratio, with no
free/complexed receptor pool anywhere in the ODE system. This is a genuinely
different structure from how tocilizumab is modeled in several *other* files
in this corpus (a reminder that the naming convention is keyed to what each
file actually does, not to the drug identity).

## A genuine structural quirk, preserved rather than "fixed" (all five drugs)

Every one of the five drug compartments in the original is dosed with the mg
dose pre-divided by its own volume of distribution *before* being injected
(see `make_doses()`'s `TOCI_conc`, `SILT_conc`, `DEX_conc`, `RUX_amt`,
`ANA_amt` conversions), and every effect term reads that same compartment
value directly, e.g.:

```r
double Inh_TOCI = (TOCI_C1 > 0) ? Emax_TOCI * TOCI_C1 / (EC50_TOCI + TOCI_C1) : 0.0;
```

There is no `amt / V1` division anywhere in `$ODE`; each drug's `V1_<STEM>`
only ever appears inside the elimination-rate ratio `CL_<STEM>/V1_<STEM>`
(and, for the two-compartment drugs, in `Q/V1`, `Q/V2`). So in this file
every drug compartment *is* the concentration by construction — same pattern
documented in `sepsis/sep_refactor_notes.md` (tocilizumab), and again in
`dengue/denv_refactor_notes.md` and `eosinophilic-esophagitis/eoe_refactor_notes.md`.
Applying the naming convention's literal template line
(`C_<STEM> = CENT_<STEM> / V1_<STEM>`) here would silently change every
downstream cytokine trajectory by a factor of that drug's own volume — the
ODE state trajectories are identical either way (only `CL/V1` and `Q/V1`,
`Q/V2` ratios enter the dynamics), but the *effect* term is not. Preserved
exactly: `C_TOCI = CENT_TOCI`, `C_SILT = CENT_SILT`, `C_DEX = CENT_DEX`,
`C_RUX = CENT_RUX`, `C_ANA = CENT_ANA` — all undivided.

**A second quirk, found and preserved rather than fixed:** `Emax_RUX_JAK`'s
and `Emax_ANA_IL1`'s corresponding bioavailability parameters, `F_RUX` and
`F_ANA`, are declared in `$PARAM` (0.95 each) but are **never applied to the
absorbed dose anywhere in `$ODE`** — `dxdt_RUX_C1`/`dxdt_ANA_C1` read
`ka_RUX * RUX_GUT`/`ka_ANA * ANA_SC` directly, with no `F_*` factor. This
means both drugs behave as if fully bioavailable (F=100%) regardless of the
declared value. Not a compile defect, just a modeling quirk of the original
— kept identically as `F_RUX`/`F_ANA` in the refactored `$PARAM` (renamed to
match convention, still unused in `$ODE`) rather than silently applying it
(which would change the model's numeric behavior) or silently dropping it
(which would remove a documented parameter the original declared).

## Naming map

| Original | Refactored |
|---|---|
| `TOCI_C1` / `TOCI_C2` | `CENT_TOCI` / `PERI_TOCI` |
| `Emax_TOCI` | `EMAX_TOCI` |
| `Inh_TOCI` | `EFFECT_TOCI` |
| `SILT_C1` | `CENT_SILT` |
| `Emax_SILT` | `EMAX_SILT` |
| `Inh_SILT` | `EFFECT_SILT` |
| `DEX_C1` / `DEX_C2` | `CENT_DEX` / `PERI_DEX` |
| `Emax_DEX_IL6` / `Emax_DEX_IFNg` | `EMAX_DEX_IL6` / `EMAX_DEX_IFNG` |
| `Inh_DEX_IL6` / `Inh_DEX_IFNg` | `EFFECT_DEX_IL6` / `EFFECT_DEX_IFNG` |
| `RUX_GUT` / `RUX_C1` / `RUX_C2` | `GUT_RUX` / `CENT_RUX` / `PERI_RUX` |
| `ka_RUX`, `Emax_RUX_JAK`, `EC50_RUX_JAK` | `KA_RUX`, `EMAX_RUX`, `EC50_RUX` |
| `Inh_RUX` | `EFFECT_RUX` |
| `ANA_SC` / `ANA_C1` | `GUT_ANA` / `CENT_ANA` |
| `ka_ANA`, `Emax_ANA_IL1` | `KA_ANA`, `EMAX_ANA` |
| `Inh_ANA` | `EFFECT_ANA` |

`CL_<STEM>`, `V1_<STEM>`, `V2_<STEM>`, `Q_<STEM>`, `EC50_<STEM>` already
matched the convention for all five compounds and are unchanged in value and
name. `GAMMA_TOCI`, `GAMMA_SILT`, `GAMMA_DEX`, `GAMMA_RUX`, `GAMMA_ANA` are
new, added explicitly at `1` (the original had no Hill exponent term for any
of the five). `use_<STEM>` switches and `<STEM>_dose_mg` R-side dose
parameters are unchanged (outside the naming convention's scope; also
unreferenced inside `$MAIN`/`$ODE`/`$TABLE` in both the original and the
refactor — the effect is gated entirely by whether/how much drug has been
dosed, same as the `useToci` observation in `sepsis/sep_refactor_notes.md`).

**Dexamethasone has two effect terms from one concentration**, not one:
`EFFECT_DEX_IL6` and `EFFECT_DEX_IFNG` both read `C_DEX`, sharing `EC50_DEX`
but each with its own `EMAX_DEX_*`, exactly as the original computed two
separate `Inh_DEX_*` values from the same `DEX_C1`. This isn't a second
"tissue site" in the guide's sense (same compartment, same concentration) —
it's one compound acting on two distinct downstream pathways (IL-6 and
IFN-γ suppression) with different potency ceilings. Handled as two named
`EFFECT_DEX_*` variables rather than forcing a single shared term, since
collapsing them would either lose the IFN-γ-specific `Emax` or invent a
combined one the original never had.

## Hill interface

All five compounds' effect terms were already exact Hill ratios in the
original (`Emax*C/(EC50+C)`, implicit `gamma=1`) — every `EFFECT_<STEM>`
below is a **rename, not a refit**:

```
double EFFECT_TOCI = (C_TOCI > 0) ?
    EMAX_TOCI * pow(C_TOCI, GAMMA_TOCI) / (pow(EC50_TOCI, GAMMA_TOCI) + pow(C_TOCI, GAMMA_TOCI)) : 0.0;
```

(and identically for `SILT`, `DEX_IL6`, `DEX_IFNG`, `RUX`, `ANA`). The
original's `(C > 0) ? ... : 0.0` zero-floor guard was preserved on every term
— with `GAMMA_<STEM>=1` this guard is not mathematically redundant once a
concentration could dip slightly negative from solver noise (`pow(x,1)=x`
even for negative `x`, unlike the original's explicit floor at zero), so
dropping it could have introduced a difference from the original at exactly
those noisy points. Kept for guaranteed exact equivalence rather than relying
on `pow()`'s behavior at negative arguments.

`Inh_IL6_sig` (the point where `EFFECT_TOCI`, `EFFECT_SILT`,
`EFFECT_DEX_IL6`, `EFFECT_RUX` are combined into one shared IL-6-signaling
inhibition term) is preserved as the original computed it, at its original
location right after the four individual effect terms — each compound's
`EFFECT_<STEM>` stays independently driveable up to that one combination
point, per the guide's "combine only at the point the disease equations
actually use them" rule.

## Build-compatibility fix (disclosed; syntax-only, non-numeric)

**The original `crs_mrgsolve_model.R` does not compile under mrgsolve
2.0.1** — logged as `translations/UPSTREAM_ISSUES.md` #75. `$CAPTURE` lists
ten names that are already `$CMT` compartments (`CRS_SEV`, `ORGAN_DMG`,
`CART_ACT`, `TUMOR`, `VASC_PERM`, `ENDO_ACT`, `MAC_ACT`, `STAT3`, `GMCSF`,
`IL2`), which mrgsolve 2.0.1's `validObject()` rejects outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: CRS_SEV,ORGAN_DMG,CART_ACT,TUMOR,VASC_PERM,ENDO_ACT,MAC_ACT,STAT3,GMCSF,IL2
```

Per the guide's settled policy, the fix — dropping those ten names from
`$CAPTURE` — was applied **directly to the delivered
`crs_mrgsolve_model_refactored.R`** (confirmed via `/model_manifest`'s
`outputPaths` that every one of the ten compartments is still reported
regardless, so this changes nothing about what is available, only what
compiles). The checked-in original, `crs_mrgsolve_model.R`, is untouched and
still carries the defect exactly as written; the same fix was applied only
to an in-memory-only copy of the original's DSL sent to the qspserver API
for comparison purposes, never saved to disk.

## Verification

**Method**: qspserver `mrgsolve_api` (`http://localhost:8007`),
`/model_manifest` (build-time check that `CL_/V1_/V2_/Q_/KA_/F_/EMAX_/EC50_/
GAMMA_` are declared in `$PARAM` for all five compounds — confirmed) and
`/run_simulation` (numeric comparison), extracting the bare DSL text from
both the original (with the build-compat fix above applied, in-memory only)
and the refactored file, run through identical dosing. Requests were spaced
~2s apart per the API's stated concurrency limits; no failures or crashes
encountered.

Four of the five compounds were verified against the **original file's own
scenarios** (from its `scenarios`/`make_doses()` R code):

| Compound(s) tested together | Scenario | Dosing | Result |
|---|---|---|---|
| Tocilizumab alone | 2 | 180.645 μg/mL-equiv bolus (1h infusion) Day 3, 90.32-equiv Day 5 | max abs/rel diff **0.0**, 19 shared outputs, 115 timepoints |
| Tocilizumab + Dexamethasone | 3 | TOCI Day 2 + Day 4 (half dose); DEX Day 2 pulse + daily ×5 (0.8×) | max abs/rel diff **0.0**, 19 shared outputs, 120 timepoints |
| Ruxolitinib + Dexamethasone (+ Tocilizumab, per the original's own scenario 4) | 4 | RUX BID Days 2–20; DEX Day 2 pulse (2×) + daily ×8; TOCI Day 2 | max abs/rel diff **0.0**, 19 shared outputs, 160 timepoints |
| Anakinra + Dexamethasone | 5 | ANA QD Days 2–20; DEX Day 2 pulse + daily ×6 (0.75×) | max abs/rel diff **0.0**, 19 shared outputs, 139 timepoints |

**Siltuximab has no scenario in the original file that doses it** —
`make_doses()` only handles scenarios 1–5 (Untreated / TOCI / TOCI+DEX /
RUX+DEX+TOCI / ANA+DEX); siltuximab's PK/PD block exists in the model but is
never exercised by any of the file's own five scenarios. Rather than
inventing an unrelated dosing regimen, siltuximab was verified with a
**bespoke bolus test using the file's own dose-conversion pattern** (the
already-defined but unused `SILT_dose_mg = 770` → `770*1000/5500 =
140.0` μg/mL-equivalent, dosed as a single 1h-infusion bolus at Day 2, the
same infusion-duration convention the other four IV drugs use): max abs/rel
diff **0.0**, 19 shared outputs, 114 timepoints.

**Result: exact match (max abs/rel deviation 0.0) for all five compounds**,
across every scenario tested — consistent with the expectation for
Archetypes 1–3 with a pure Hill-rename (no fitting): all five compounds'
`Inh_*`/`Emax*C/(EC50+C)` terms were already exact Hill ratios, so renaming
and reorganizing introduces no numeric difference.

No fitting was needed anywhere in this file (no TMDD, no ODE-solved
receptor kinetics) — the guide's "fit an approximation" path for Archetype 4
does not apply here for any of the five compounds.

## Manifest check

`/model_manifest` on the refactored DSL confirms every naming-convention
parameter is discoverable: `CL_TOCI, V1_TOCI, Q_TOCI, V2_TOCI, EMAX_TOCI,
EC50_TOCI, GAMMA_TOCI, CL_SILT, V1_SILT, EMAX_SILT, EC50_SILT, GAMMA_SILT,
CL_DEX, V1_DEX, Q_DEX, V2_DEX, EMAX_DEX_IL6, EC50_DEX, EMAX_DEX_IFNG,
GAMMA_DEX, KA_RUX, CL_RUX, V1_RUX, Q_RUX, V2_RUX, F_RUX, EMAX_RUX, EC50_RUX,
GAMMA_RUX, KA_ANA, CL_ANA, V1_ANA, F_ANA, EMAX_ANA, EC50_ANA, GAMMA_ANA` all
appear in `$PARAM`. `C_<STEM>`/`EFFECT_<STEM>` for all five compounds are
computed in `$ODE` (not `$PARAM` — following the established pattern in
every prior refactor checked, e.g. `sepsis/sep_mrgsolve_model_refactored.R`,
`dengue/denv_mrgsolve_model_refactored.R`,
`thyroid-eye-disease/ted_mrgsolve_model_refactored.R`) and exposed via
`$CAPTURE` (`C_TOCI EFFECT_TOCI C_SILT EFFECT_SILT C_DEX EFFECT_DEX_IL6
EFFECT_DEX_IFNG C_RUX EFFECT_RUX C_ANA EFFECT_ANA`), confirmed present in
`/model_manifest`'s `outputPaths`.
