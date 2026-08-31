# Refactor notes — `toxic-alcohol-poisoning/tap_mrgsolve_model.R` / fomepizole (FOM)

**Scope of this refactor**: only fomepizole's own PK compartments and its
one downstream effect (competitive inhibition of the shared alcohol
dehydrogenase, ADH, flux) were touched. Every other compound in this file
(methanol, ethylene glycol, ethanol, folinic acid, thiamine, pyridoxine,
haemodialysis/CRRT, exogenous bicarbonate, calcium) and the entire
disease network (formate/glycolate/glyoxylate/oxalate chemistry,
acid-base, CNS/optic/renal injury, hazards) is copied verbatim from the
original into `tap_mrgsolve_model_refactored.R`.

## The census classification for this row was wrong — corrected here

`driver-patches/data/compound_perturbation_census.md` classified FOM as
**"Delete PK compartment; concentration is itself the state"** — i.e. the
premise that fomepizole's own ODE compartment *is* the concentration
directly, with no amount/volume division to redirect. That is not what
the actual code does. In the original `$ODE`:

```
double C_FOM = A_FOM1  / V1FOM; double C_FOM2= A_FOM2  / V2FOM;
...
double abs_F  = KAFOM * A_GUTF;
double q_FOM  = 25.0 * (C_FOM - C_FOM2);
dxdt_A_GUTF = -abs_F;
dxdt_A_FOM1 = abs_F + fom_inf - q_FOM - el_fom - cl_hd_fom * C_FOM;
dxdt_A_FOM2 = q_FOM;
```

`A_GUTF`/`A_FOM1`/`A_FOM2` are ordinary amount states (mmol) with a real
absorption depot, a real central/peripheral exchange, and an explicit
`C_FOM = A_FOM1 / V1FOM` mass-to-concentration division — the textbook
shape the census's "Delete PK compartment" bucket is defined to exclude,
not an example of it. Per the fork guide's own header warning on this
table ("a second, independent reimplementation of this same classifier
produced materially different aggregate totals in the long tail...most
notably on `delete_compartment` classification") and the shared-rules
principle that a gate mismatch is a question to resolve, not a verdict to
force through: this row was **re-classified and refactored as Archetype 3
(depot + central + peripheral)**, with the original's own nonlinear
elimination term preserved as a disclosed bespoke variant (see below), not
linearised to fit the archetype template.

Recorded back into the census (see bottom of this file).

## Archetype: 3 (depot + central + peripheral), with a bespoke non-linear elimination term kept from the original

Topology matches Archetype 3 exactly (`GUT_FOM` -> `CENT_FOM` <->
`PERI_FOM`, first-order absorption, first-order inter-compartmental
exchange). The one respect in which it is **not** a clean Archetype 3 is
elimination: the guide's template assumes linear `CL_<STEM>/V1_<STEM>`
elimination, but the original eliminates fomepizole by **saturable,
autoinducible Michaelis-Menten kinetics** — a real, load-bearing feature
of fomepizole's actual clinical PK (it induces its own metabolism), not
an arbitrary embellishment:

```
double autoi  = 1.0 + AUTOIND / (1.0 + exp(-(SOLVERTIME - T_AUTO) / TAU_AUTO));
double el_fom = VMAX_FOM * autoi * C_FOM / (KM_FOM + C_FOM);
```

Per the guide's "match the archetype closest to what the original already
does; don't add or remove compartments (or, by the same logic, flatten
real non-linear elimination) to make it fit" — this term is kept exactly
as written, renamed only where it touches a renamed symbol (`C_FOM`
unchanged; `AUTOIND`/`T_AUTO`/`TAU_AUTO` unchanged, already effectively
FOM-scoped and not part of the guide's named-role table). This mirrors
`acne-vulgaris/acn_refactor_notes.md`'s handling of isotretinoin's
disclosed bespoke shared-Hill structure: reorganize and rename what the
convention actually covers, disclose what it doesn't rather than forcing
a linear fit that would silently change the model's autoinduction
behaviour.

## Renaming

| Original | Refactored |
|---|---|
| `A_GUTF` | `GUT_FOM` |
| `A_FOM1` | `CENT_FOM` |
| `A_FOM2` | `PERI_FOM` |
| `KAFOM` | `KA_FOM` |
| `V1FOM` / `V2FOM` (local `$ODE` doubles, `0.42*WT` / `0.28*WT`) | `V1_FOM` / `V2_FOM` (same computation, still local doubles — see below) |
| `25.0` (hardcoded inter-compartmental clearance in `q_FOM`) | `Q_FOM` (new named `$PARAM`, value `25.0`, unchanged) |
| `inh` (local `$ODE` double, `C_FOM/KI_FOM`) | `EFFECT_FOM` (the required named effect variable) |
| `C_FOM` (already correctly named) | unchanged, now the guide's required single redirect site, defined identically (`CENT_FOM / V1_FOM`) |

`VMAX_FOM`, `KM_FOM`, `KI_FOM`, `CL_HD_FOM`, `AUTOIND`, `T_AUTO`,
`TAU_AUTO` already matched or were already effectively FOM-scoped (no
other compound in this file uses them) and are unchanged. All parameter
*values* are copied verbatim from the original; `Q_FOM = 25.0` is a
promotion of an existing hardcoded literal into a named parameter with
the identical value (same treatment `abdominal-aortic-aneurysm/
aaa_refactor_notes.md` gave the Statin arm's hardcoded NF-kB/ROS ratios),
not an invented value.

**`V1_FOM`/`V2_FOM` were kept as local, WT-scaled `$ODE` doubles, not
promoted to `$PARAM`.** Every one of this file's seven modeled compounds
computes its own central/peripheral volumes the same way (`double V1M =
0.35*WT`, `double V1E = 0.38*WT`, ..., `double V1_FOM = 0.42*WT`) — this
is a whole-file structural convention, not a FOM-specific quirk, and
promoting only FOM's volumes to `$PARAM` would create an inconsistency
the rest of the file doesn't have. The guide's mandatory qspserver
discoverability requirement applies to `C_<STEM>`/`EFFECT_<STEM>`, not to
every intermediate PK constant, so this is compliant.

## The Hill interface: disclosed as linear, not forced into Emax/EC50/gamma shape

Fomepizole's "effect" on the disease is not a downstream Emax-style
response — it is a competitive-inhibition term inside the *shared* ADH
denominator that methanol, ethylene glycol, and ethanol all divide by:

```
double sM  = C_M  / KM_M;
double sE  = C_E  / KM_E;
double sET = C_ET / KM_ET;
capture EFFECT_FOM = C_FOM / KI_FOM;
double den = 1.0 + sM + sE + sET + EFFECT_FOM;
```

This is genuinely linear in `C_FOM` (no saturation of the inhibition
itself — what saturates is the *shared* flux `den` divides into), exactly
the situation the guide's Hill-interface section defers to when an
effect "emerges" in a non-Hill shape: `acne-vulgaris/acn_refactor_notes.md`
(SPI's AR-antagonism term) made the identical call — "kept as linear
C_SPI/KI_SPI_AR term ... not forced into a fabricated Hill." No `EMAX_FOM`/
`EC50_FOM`/`GAMMA_FOM` were added; naming it `EFFECT_FOM` and exposing it
as the single named quantity other code reads satisfies the guide's actual
requirement (one named, redirectable effect variable) without inventing a
shape the mechanism doesn't have. No curve-fitting was needed or
performed.

## Exposure via `capture`, not `$PARAM` or a bare `$CAPTURE` block

`C_FOM` and `EFFECT_FOM` are computed from live ODE state every step, so
per the same reasoning already established in `essential-thrombocythemia/
et_refactor_notes.md`, `x-linked-hypophosphatemia/xlh_refactor_notes.md`,
`mucopolysaccharidosis-type-1/mps1_refactor_notes.md`, and
`primary-sclerosing-cholangitis/psc_refactor_notes.md` for this engine,
they cannot also live in `$PARAM` (collides at the C++ stage). This file
has no pre-existing bare `$CAPTURE name-list` block at all (its own
`$TABLE` instead uses the modern mrgsolve `capture NAME = expr;` inline
syntax throughout), so the two lines that compute them were changed from
`double` to `capture` in place:

```
capture C_FOM = CENT_FOM / V1_FOM; double C_FOM2 = PERI_FOM / V2_FOM;
...
capture EFFECT_FOM = C_FOM / KI_FOM;
```

Confirmed via `POST /model_manifest`: both appear in `outputPaths`
(`C_FOM`, `EFFECT_FOM`), alongside `GUT_FOM`/`CENT_FOM`/`PERI_FOM` and
every renamed `$PARAM` (`KA_FOM`, `Q_FOM`, plus the unchanged
`VMAX_FOM`/`KM_FOM`/`KI_FOM`/`CL_HD_FOM`).

## Build-compat fix (pre-existing defect, unrelated to FOM) — `translations/UPSTREAM_ISSUES.md` #108

The original does not compile under mrgsolve 2.0.1 at all, for a reason
entirely outside fomepizole's own block: `$ODE`'s bicarbonate-titration
code declares a local named `_e`, and this build rejects any leading
underscore (`Error: leading underscore not allowed: _e`):

```
double _e  = (pH - PH_TARGET) / PH_TAU;
if (_e >  50.0) _e =  50.0;
if (_e < -50.0) _e = -50.0;
double bic_eff = bic_inf / (1.0 + exp(_e));
```

Per the guide's settled policy for this situation: **not fixed in
`tap_mrgsolve_model.R`** (left completely untouched); fixed directly in
the delivered `tap_mrgsolve_model_refactored.R` by renaming `_e` to
`bic_z` at all four occurrences (syntax-only, non-numeric — values and
control flow identical); logged as `translations/UPSTREAM_ISSUES.md` #108.

## R-side changes

`fom_events()`'s three dosing-event rows (`cmt = "A_FOM1"`) were updated
to `cmt = "CENT_FOM"`, matching the renamed compartment. No other R-side
code references fomepizole's compartments by name. `tap_run()`'s dynamic
`match(as.character(d$cmt), mrgsolve::cmt(mod))` resolves this
automatically; no other change was needed for the scenario library,
`tap_summary()`, `tap_verify()`, or the demonstration block.

## Verification

Per the guide's mandatory protocol, run via the qspserver `mrgsolve_api`
container (`http://localhost:8007`, `POST /model_manifest` and `POST
/run_simulation`) — no local R/mrgsolve was used. Both the original and
the refactored DSL were extracted from their respective `code <- '...'`
R-string blocks; the original was given the identical, minimal,
in-memory-only `_e` -> `bic_z` rename (never written to
`tap_mrgsolve_model.R`) purely so it would compile side by side with the
refactored model, per the guide's symmetric-workaround instruction.

**Scenario**: `M4_fom_hd` — the same scenario the original file's own
`tap_verify()` acceptance test uses (methanol 0.7 g/kg, fomepizole
loading dose at 8 h with q12h maintenance boosted to q4h on dialysis, and
a haemodialysis session at 9-15 h) — reproduced from `tap_scenarios()`
exactly: `A_GUTM` bolus of 1529.24 mmol at t=0; nine fomepizole doses at
t = 8, 13, 20, 32, 44, 56, 68, 80, 92 h (12.789/8.526/12.789 mmol per the
original's own `fom_events()` load/maintenance/HD-boost/autoinduction-
multiplier logic) into `CENT_FOM`; `HD1ON/OFF = 9/15`; bicarbonate
infusion 30 mmol/h from 8-26 h. Run 0-96 h (the scenario's own `tend`) at
`delta = 1 h` (97 points; well under the API's solver-step ceiling, so no
window shortening was needed).

**Result: exact match.** All 40 of the original's own shared `$TABLE`
outputs compared (`MEOH`, `FORMmM`, `FORMmgL`, `FCNS`, `FVIT`, `FOMug`,
`pHart`, `HCO3o`, `PACO2o`, `pHur`, `LACTo`, `NAo`, `AG`, `OG`, `DAG`,
`GAPSUM`, `OGoverAG`, `ADHDEN`, `INHFAC`, `FLUXM`, `FLUXE`, `FLUXLEFT`,
`FOMKIRAT`, `CAIONo`, `CATOTmg`, `SSPLAS`, `SSTUB`, `CAOXKID`, `GFRpct`,
`PTINJo`, `ATPo`, `PUTo`, `OPTICo`, `LOGMAR`, `COMA`, `PDEATH`, `PBLIND`,
`FATAL`, `THFo`, `MEOHOX`) had **max absolute/relative deviation = 0.0**
across every one of the 97 time points — the expected outcome for a pure
structural reorganization with a disclosed-linear (not fitted) effect
term, per the guide's tolerance table. This confirms `EFFECT_FOM`, as
actually consumed inside `den` and hence the whole ADH-competition
cascade, drives the ODE integration identically to the original's `inh`.

**One cosmetic, already-precedented artifact, not a correctness issue:**
the two newly-`capture`d diagnostic columns `C_FOM`/`EFFECT_FOM`
themselves (which have no counterpart to verify against in the original —
the original never exposed a raw `C_FOM`/`inh` output) report a
one-internal-step-stale value at the exact instant of each of the nine
dosing events (e.g. at t = 8 h, `C_FOM` reads 0 while the independently
`$TABLE`-recomputed `FOMug` = `cFM*MWFOM` already reflects the dose,
35.71 vs the expected 0; by the very next output row the two agree to
within JSON-rounding precision, and stay in agreement until the next dose
event). This is the same `$MAIN`/`$TABLE`(here, `$ODE`-`capture`)
evaluation-timing nuance already disclosed as a non-bug in
`sepsis/sep_refactor_notes.md` for `C_TCZ` vs `TOCI_C` — more visible here
only because this scenario doses fomepizole nine times instead of once.
It does not affect any of the 40 verified disease-relevant outputs above
(all exact matches at every point, including the dose timestamps
themselves), confirming the actual ODE state is correct throughout;
only the diagnostic capture column momentarily lags the compartment state
it derives from. Worth keeping in mind if a future patch reads `C_FOM`/
`EFFECT_FOM` directly from captured output at a timestep coincident with
a dose, rather than from live model state.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`toxic-alcohol-poisoning | FOM` row: perturbation method corrected from
"Delete PK compartment; concentration is itself the state" to "Redirect
concentration (clean single site)", target/pathway set to "Alcohol
dehydrogenase (competitive inhibitor)", and the last three columns filled
in with the redirect site, unit, and a summary of the archetype/
verification outcome above.
