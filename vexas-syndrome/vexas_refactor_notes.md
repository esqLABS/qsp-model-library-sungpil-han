# VEXAS syndrome — PK/PD refactor notes

Refactored file: `vexas_mrgsolve_model_refactored.R`
Original (untouched): `vexas_mrgsolve_model.R`
Compounds in scope (per `driver-patches/data/compound_perturbation_census.md`,
all classified "Redirect concentration (clean single site)"): **ANA**
(anakinra), **AZA** (azacitidine), **CAN** (canakinumab), **PRED**
(prednisone), **RUX** (ruxolitinib), **TOC** (tocilizumab).

## Archetype per compound

All six compounds use the **same PK structure** in the original: an SC or
oral absorption depot feeding a single central compartment, linear
elimination, no peripheral compartment. This is the guide's **Archetype 3
(depot + central, linear), minus the peripheral compartment** — i.e. the
`GUT_TCZ`/`CENT_TCZ` two-compartment variant the guide names explicitly
("drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a 2-compartment-no-depot variant" runs
the other direction here: this file keeps the depot but never had a
peripheral compartment to begin with).

| Compound | Original names | Refactored names | Archetype |
|---|---|---|---|
| PRED | `ka_PRED`, `V_PRED`, `CL_PRED`, `AGUT_PRED`, `CC_PRED` | `KA_PRED`, `V1_PRED`, `CL_PRED`, `GUT_PRED`, `CENT_PRED` | 3 minus peripheral |
| TOC  | `ka_TOC`, `V_TOC`, `CL_TOC`, `ADEP_TOC`, `CC_TOC` | `KA_TOC`, `V1_TOC`, `CL_TOC`, `GUT_TOC`, `CENT_TOC` | 3 minus peripheral |
| ANA  | `ka_ANA`, `V_ANA`, `CL_ANA`, `ADEP_ANA`, `CC_ANA` | `KA_ANA`, `V1_ANA`, `CL_ANA`, `GUT_ANA`, `CENT_ANA` | 3 minus peripheral |
| CAN  | `ka_CAN`, `V_CAN`, `CL_CAN`, `ADEP_CAN`, `CC_CAN` | `KA_CAN`, `V1_CAN`, `CL_CAN`, `GUT_CAN`, `CENT_CAN` | 3 minus peripheral |
| RUX  | `ka_RUX`, `V_RUX`, `CL_RUX`, `AGUT_RUX`, `CC_RUX` | `KA_RUX`, `V1_RUX`, `CL_RUX`, `GUT_RUX`, `CENT_RUX` | 3 minus peripheral |
| AZA  | `ka_AZA`, `V_AZA`, `CL_AZA`, `ADEP_AZA`, `CC_AZA` | `KA_AZA`, `V1_AZA`, `CL_AZA`, `GUT_AZA`, `CENT_AZA` | 3 minus peripheral |

### TOC and CAN: checked for TMDD, neither has it

The task brief for this batch flagged TOC and CAN specifically for a TMDD
check (TOC as "TMDD-lite" per the original file's own header comment; CAN
as an anti-IL-1beta antibody). Reading the actual `$ODE` block: neither
compound has a receptor compartment, a `KON`/`KOFF` binding pair, or any
ODE-solved occupancy dynamic. Both are:

```
double C_TOC = CC_TOC / V_TOC;
double occ_TOC = C_TOC / (Kd_TOC + C_TOC);
```
(and identically for CAN with its own `Kd_CAN`) — a plain algebraic
Michaelis-Menten-style occupancy ratio computed directly from plasma
concentration, not a dynamic receptor-binding system. This is **not**
TMDD despite the original's own "TMDD-lite" comment in its file header;
that header comment overstates what the equations actually do. Both are
handled as Archetype 3 minus peripheral (plain PK) with the Hill interface
built directly from `C_<STEM>`, not from any occupancy state variable
(there is none to build it from).

## The Hill interface — rename, not a refit

Every one of the five wired-up compounds (PRED, TOC, ANA, CAN, RUX) already
had its effect on the disease system as a plain algebraic ratio of the
shape `Emax * C / (EC50 + C)`, just with the `Emax` and `EC50` baked in as
literal numbers inside a larger expression instead of named parameters:

| Compound | Original inline expression | Refactored | Emax | EC50 |
|---|---|---|---|---|
| PRED | `(1.0 - 0.6*C_PRED/(0.05+C_PRED))` inside `NFkB_act` | `EFFECT_PRED`, `NFkB_act = 1 - EFFECT_PRED` | `EMAX_PRED=0.6` | `EC50_PRED=0.05` |
| TOC  | `occ_TOC*0.95` inside `IL6_eff` | `EFFECT_TOC` | `EMAX_TOC=0.95` | `EC50_TOC=Kd_TOC=0.5` |
| ANA  | `occ_ANA*0.90` inside `IL1_eff` | `EFFECT_ANA` | `EMAX_ANA=0.90` | `EC50_ANA=Kd_ANA=1.0` |
| CAN  | `occ_CAN*0.85` inside `IL1_eff` | `EFFECT_CAN` | `EMAX_CAN=0.85` | `EC50_CAN=Kd_CAN=0.3` |
| RUX  | `occ_RUX*0.70` inside `IFN_eff` and `dxdt_IFNa` | `EFFECT_RUX` | `EMAX_RUX=0.70` | `EC50_RUX=Kd_RUX=100` |

`GAMMA_<STEM> = 1` for all five (the original had no explicit Hill
exponent anywhere). **This is a pure rename, not a fit** — no `nls()` was
run and no R^2 is reported, because there was nothing to fit: pulling
`0.6`/`0.05`, `0.95`/`Kd_TOC`, etc. out of the inline expression into named
`EMAX_<STEM>`/`EC50_<STEM>` parameters and writing
`EMAX*pow(C,GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))` with `GAMMA=1` is
algebraically identical to the original expression term-for-term
(`pow(x,1) == x`). Confirmed by the exact-match verification below.

One simplification, disclosed: the original's `NFkB_act` used a defensive
ternary,
```
(1 + C_PRED/(0.05 + C_PRED) > 0 ? (1.0 - 0.6*C_PRED/(0.05 + C_PRED)) : 1.0)
```
Since `C_PRED >= 0` always, the ratio is always in `[0,1)`, so `1 + ratio`
is always `> 0` — the guard is always true and the `: 1.0` branch is dead
code. `NFkB_act = 1.0 - EFFECT_PRED` in the refactored file always takes
the same branch the original always took; this is an exact simplification,
not a behavioral change, and the verification below confirms it produces
identical output at every timepoint in every scenario, including the
prednisone-only and prednisone-combination scenarios.

### AZA is bespoke: PK renamed, but no `EFFECT_AZA` exists

Azacitidine does **not** fit the Hill-interface pattern above, and none was
invented. Grepping the original file for every use of `C_AZA` (`CC_AZA /
V_AZA`) shows it is computed and `$CAPTURE`d but never read by any other
equation. The disease-facing clone-suppression term in
`dxdt_VAF = (k_clone*VAF*(1-VAF)) - k_HSCT*VAF - k_Aza*VAF;` uses a
separate, free-standing `k_Aza` parameter (default `0.0`) that nothing in
the shown DSL ever sets from `C_AZA`. In other words: azacitidine's own PK
profile has **zero** effect on the modeled disease state as coded — dosing
`ADEP_AZA`/`GUT_AZA` moves `C_AZA` exactly as expected, but `VAF` and every
downstream biomarker are unaffected unless a user script overrides `k_Aza`
from outside this file. This is preserved exactly: `KA_AZA`/`F_AZA`/
`CL_AZA`/`V1_AZA`/`C_AZA` are renamed to convention as normal PK, `k_Aza`
is kept as-is, and **no** `EMAX_AZA`/`EC50_AZA`/`GAMMA_AZA`/`EFFECT_AZA`
was invented to paper over the gap — the guide is explicit that PK
parameter values (and, by the same logic, PD wiring) must come from the
original, never invented. Logged as a new defect,
`translations/UPSTREAM_ISSUES.md` entry **#86** (see below).

## Upstream defects found (logged, not fixed)

Two new entries appended to `translations/UPSTREAM_ISSUES.md` (previous
tail was #84; re-checked immediately before appending, no concurrent
writes found):

- **#85** — every one of the six compounds declares an SC/oral
  bioavailability parameter (`F_PRED`, `F_TOC`, `F_ANA`, `F_CAN`, `F_RUX`,
  `F_AZA`) that is never referenced in any `dxdt` equation anywhere in the
  file. Confirmed by grep (each name appears exactly once, in its own
  `$PARAM` declaration) and by `/model_manifest` on the untouched original
  (all six listed with their documented defaults). Practical effect: every
  dose in every one of the file's own scenarios reaches the central
  compartment as if `F=1`, regardless of the documented value (e.g.
  tocilizumab's stated 80% SC bioavailability is not actually applied).
  **Preserved bug-for-bug** in the refactored file (no `F_<STEM>`
  multiplier added to any `dxdt_CENT_<STEM>` equation) — confirmed
  numerically identical to the original in the verification below.
- **#86** — azacitidine's PK is fully disconnected from its own PD, as
  described in the AZA section above. Not fixed; `C_AZA` is still exposed
  as the redirect point per the naming convention, `k_Aza` is unchanged,
  and no effect term was invented.

A third, much lower-impact oddity was **not** escalated to
`UPSTREAM_ISSUES.md`, only noted here (same treatment the prior
tocilizumab refactor gave an analogous unused `RTOT_TCZ` parameter): the
original declares `IL6Rmax = 5.0` ("Apparent IL-6R density") in `$PARAM`
but never references it anywhere in `$ODE`/`$MAIN` — a single dead
parameter with no downstream numeric consequence, unlike the F_<STEM>
pattern above which is repeated identically across all six compounds and
directly affects simulated dose delivery. Preserved, renamed
`IL6RMAX_TOC`, in the refactored file for fidelity.

**No mrgsolve 2.0.1 build defect was found.** Both the original and the
refactored DSL compile cleanly through `POST /model_manifest` on the
qspserver `mrgsolve_api` container with no changes needed — this is not a
"doesn't compile" situation, so no syntax-only build-compat fix was
required or applied.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` on the refactored DSL through the
qspserver `mrgsolve_api` container (`http://localhost:8007`):

- `outputPaths` lists every renamed `$CMT` (`GUT_PRED`, `CENT_PRED`, ...,
  `GUT_AZA`, `CENT_AZA`) in the same order as the original's
  `AGUT_PRED`/`CC_PRED`/.../`ADEP_AZA`/`CC_AZA` — same compartment count
  (31 disease+PK states), same relative order, so 1-based `cmt` dosing
  indices transfer directly between original and refactored for
  verification purposes.
- `outputPaths` also lists all six `C_<STEM>` and all five `EFFECT_<STEM>`
  (no `EFFECT_AZA`, per above) plus `INFLAM`/`IL6_eff`/`IL1_eff`/`IFN_eff`.
- `parameters` lists every `KA_<STEM>`, `F_<STEM>`, `CL_<STEM>`,
  `V1_<STEM>` and every `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>` (for
  PRED/TOC/ANA/CAN/RUX) with their carried-over defaults — 96 parameters
  total, all discoverable.

Per the guide's `$PARAM`-vs-`$ODE` tension for `C_<STEM>`/`EFFECT_<STEM>`:
these remain plain `double` locals computed in `$ODE` (not declared in
`$PARAM`), matching the established pattern already used across this
repo's prior refactors (e.g. `abdominal-aortic-aneurysm/aaa_refactor_notes.md`,
`rheumatoid-arthritis/ra_mrgsolve_model_refactored.R`) — `$PARAM @annotated`
values compile to `const double&` in mrgsolve 2.0.1, so assigning to a
`$PARAM`-declared `C_<STEM>` from `$ODE` fails to build
(`assignment of read-only reference`). Discoverability is instead
satisfied via `$CAPTURE`, confirmed above.

## Verification (qspserver `mrgsolve_api`, `http://localhost:8007`)

The original file's own `scenarios()` function defines dosing for five of
the six compounds (PRED, TOC, ANA, RUX, AZA); it has **no** scenario that
doses CAN at all. Ran all five of the original's own scenarios, plus one
bespoke minimal CAN dosing test (150 mg SC on days 0/28/56, disclosed as
not from the original's own `scenarios()` since none exists), through
both the original and refactored DSL with identical dosing (`POST
/run_simulation`, requests spaced ~2.2s apart per the API's stated
concurrency limit):

| Scenario | Compounds | End / delta | Max relative deviation |
|---|---|---|---|
| 2_pred_taper | PRED | 2016h / 6h | 0.0 |
| 3_toci_lowGC | TOC + PRED | 2016h / 6h | 0.0 |
| 4_anakinra | ANA | 2016h / 6h | 0.0 |
| 5_ruxo_pred | RUX + PRED | 2016h / 6h | 0.0 |
| 6_azacitidine | AZA | 2016h / 6h | 0.0 |
| 7_canakinumab_bespoke | CAN | 2016h / 6h | 6.5e-13 (floating-point noise on `C_TOC`, which is analytically 0 in this CAN-only scenario and never rises above ~1e-22 in either model) |

Compared every `$CAPTURE`d output (all 19 disease-state `$CMT` variables
plus `C_PRED/C_TOC/C_ANA/C_CAN/C_RUX/C_AZA`, `INFLAM`, `IL6_eff`, `IL1_eff`,
`IFN_eff`) across the full time grid in all six runs. Every scenario
except the CAN-only bespoke run matched to `0.0` relative deviation
exactly; the one non-zero figure is floating-point noise on a
concentration that is numerically zero in both models (no tocilizumab
dosed in that scenario). This is the expected result for a pure structural
reorganization (Archetype 3 minus peripheral, no Hill-fitting) per the
guide's tolerance rule — anything beyond floating-point-scale deviation
would indicate a bug, and none was found.

Also confirmed (`3_toci_lowGC` scenario, the only one that doses two
compounds simultaneously through the shared HPA/inflammation pathway) that
the combined disease-state trajectories (VAF, MISF, ROS, ERST, IL1B, IL6,
TNFa, IFNa, CXCL8, CCL2, CRP, FER, HB, PLT, ANC, FEV, SKIN, VTE, HPA) all
match exactly, not just the drug concentrations — i.e. the rename did not
silently change how the two compounds' effects combine downstream.

## Notes for manual review

- The `"7_HSCT"` scenario in the original's own `scenarios()` is
  non-functional as written in both the original and the refactored file
  (an `evid=2` event on the `VAF` compartment with a comment admitting "use
  $TABLE override in user code" — no such override exists anywhere in the
  file). This is unrelated to any of the six compounds' PK/PD and was left
  exactly as-is (renamed nothing, since it references no PK compartment);
  not logged as a defect here since it is pre-existing, disclosed-as-a-stub
  behavior with a comment admitting it, not a silent one.
- Every parameter value in the refactored file is copied verbatim from the
  original — none invented, none defaulted, none dropped.
