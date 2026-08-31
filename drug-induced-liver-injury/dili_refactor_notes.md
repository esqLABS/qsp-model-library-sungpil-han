# Refactor notes — `drug-induced-liver-injury/dili_mrgsolve_model.R`

## Scope: three rows checked, one census entry corrected, one row missing from the census entirely

This file carries two rows in `driver-patches/data/compound_perturbation_census.md`:
**"Innate immunity (EL)"** and **NAC**. Per the task's specific flag, both
were checked against the actual code rather than trusted at face value, and
the file was also checked for a missing "culprit hepatotoxin" compound
(common in this model class) — one was found: **APAP**.

| Row | Real drug? | Evidence |
|---|---|---|
| "Innate immunity (EL)" | **No — not a compound at all.** No compartment, parameter, or identifier named `EL` exists anywhere in the file (confirmed by exact-name grep). `DAMP`, `KC`, `TNF`, `IL10` (this file's innate-immunity compartments) are all endogenous — produced by `KDAMP*death`, `KKC*DAMP*(1-KC)`, `KTNF*KC`, `KIL10P*KC` respectively — with no depot, no absorption parameter, and no `ev()`/dataset dosing route anywhere. "EL" is a substring fragment shared by three unrelated elimination-rate parameters that all sit inside that same disease section (`KDAMP_EL`, `KTNF_EL`, `KIL10_EL`); the classifier appears to have mis-extracted this fragment as if it were a compound stem. |
| NAC | **Yes.** Own compartment `NACC` (dosed via the file's own `nac_prescott()` 3-phase IV Prescott regimen builder, used in scenarios S5–S7, S9), own volume/clearance (`VNAC`/`CLNAC`), own Hill-shaped effect term (`enac = EMAX_NAC*CNAC/(EC50_NAC+CNAC)`), used in two disease equations (cysteine-supply boost and direct ROS scavenging). Census's own classification ("Redirect concentration inside #define macro") was already correct. |
| APAP | **Yes — missing from the census entirely.** The file's only drug-dosing route (`AGUT`→`AC`→`AP`/`ALIV`, dosed via the R wrapper's `oral()` helper, used in every one of the file's 13 scenarios) has no census row at all. This is the file's own "culprit hepatotoxin" — the model's entire thesis (reactive-metabolite flux vs. thiol resynthesis, the bistable JNK-Sab loop, BSEP-mediated cholestasis) is downstream of this compound's own hepatic bioactivation. Confirmed via the file's own constants (`MW_APAP <- 151.16`, `mw = MW_APAP` default in `oral()`, `CP_UGML` comment "ug/mL if MW = APAP") and calibration-anchor comments ("APAP clearance ~19-21 L/h", "VC ... [APAP ~0.9 L/kg]"). |

**Conclusion:** "Innate immunity (EL)" is corrected in the census to "Not a
compound — classifier artifact, corrected" (left completely untouched in
the code, since there is nothing to refactor). NAC is refactored as
originally classified. A new APAP row is added to the census and refactored.

## Archetype per compound

### APAP — Archetype 3 (depot + central + peripheral, linear) plus a bespoke 4th tissue compartment

The original models APAP with four PK-side compartments: `AGUT` (gut
depot), `AC` (central/plasma), `AP` (peripheral), `ALIV` (liver). The first
three are a clean Archetype-3 triad (`dxdt_AC` includes `-QD*CcP+QD*CpP`
peripheral exchange terms identical in shape to the guide's template). The
fourth, `ALIV`, is not a peripheral compartment in the ordinary sense — it
is a genuine second tissue site with its own hepatic-blood-flow exchange
(`QH*CcP - QH*CLIVP/KP`) and its own metabolism (phase-II conjugation +
CYP bioactivation consume `ALIV`, not `AC`). This is the same pattern as
`abdominal-aortic-aneurysm`'s `TISSUE_DOXY` and
`acute-intermittent-porphyria`'s `LIV_GIV`/`LIV_HEM` (see those files' own
refactor notes) — kept as its own compartment, renamed to the convention,
not collapsed into the peripheral compartment or discarded.

Renamed (all values unchanged from the original):

| Original | Refactored | Value | Role |
|---|---|---|---|
| `AGUT` (cmt) | `GUT_APAP` | -- | depot (oral) |
| `AC` (cmt) | `CENT_APAP` | -- | central (plasma) |
| `AP` (cmt) | `PERI_APAP` | -- | peripheral |
| `ALIV` (cmt) | `LIV_APAP` | -- | bespoke 4th tissue compartment (liver) |
| `KA` | `KA_APAP` | 1.2 (1/h) | absorption rate |
| `FA` | `F_APAP` | 0.88 | bioavailability x fraction absorbed |
| `CLR` | `CL_APAP` | 1.4 (L/h) | renal clearance of parent |
| `VC` | `V1_APAP` | 63.0 (L) | central volume |
| `VP` | `V2_APAP` | 25.0 (L) | peripheral volume |
| `QD` | `Q_APAP` | 25.0 (L/h) | inter-compartmental clearance |
| `KP` | `KP_APAP` | 1.0 | liver:plasma partition coefficient |
| `FULIV` | `FU_APAP` | 0.8 | unbound fraction in liver |
| `KI_BSEP` | `KI_APAP_BSEP` | 1e9 (uM) | BSEP-inhibition Ki (1e9 = inert by default) |
| `CcP` (macro) | `C_APAP` ($ODE-local) | -- | **exposed concentration 1/2: plasma** |
| `CpP` (macro) | `CPERI_APAP` ($ODE-local) | -- | peripheral concentration (internal only) |
| `CLIVP` (macro) | `CLIV_APAP` ($ODE-local) | -- | total hepatic concentration (internal only) |
| `CuP` (macro) | `CU_APAP` ($ODE-local) | -- | **exposed concentration 2/2: unbound hepatic** |
| `prod_RM` | `EFFECT_APAP` | -- | APAP's effect on disease (see below) |
| `1/(1+CuP/KI_BSEP)` (inline) | `EFFECT_APAP_BSEP` (named) | -- | APAP's 2nd effect (see below) |

**Two exposed concentrations, deliberately.** The guide allows a second
site "only when a genuinely different tissue site matters" — it does here:
plasma PK (`C_APAP`, the depot/central/peripheral triad) and hepatic
disposition (`CU_APAP`, unbound intrahepatic) are computed by different
equations and used by different consumers. `C_APAP` governs the
central<->peripheral<->liver exchange flux (pure PK). `CU_APAP` is what
every downstream toxicodynamic term actually reads (phase-II conjugation,
CYP bioactivation, BSEP inhibition) — it is the single point where an
external covariate could later substitute in for the whole reactive-
metabolite cascade.

**The metabolism-pathway parameters (`VMAX_SULT`, `KM_SULT`, `VMAX_UGT`,
`KM_UGT`, `APAPS`, `KSP`, `AUD`, `KSU`, `VMAX_CYP`, `KM_CYP`, `FCYP`) were
left named as in the original**, not given an `_APAP` suffix. They are not
part of the guide's naming-role table (no "conjugation cofactor" role
exists) and are already unambiguously scoped to APAP's own PK block (they
appear nowhere else in the file, and are used only via `CU_APAP` after this
refactor) — renaming them would add risk without resolving any naming
collision, since there is only one compound in this file that uses them.

**Not force-fit into a single Emax/EC50 Hill.** The guide's naming
convention calls for one named `EFFECT_<STEM>` function of the exposed
concentration. APAP's own header thesis is explicit that bioactivation is
deliberately *not* saturating over the dose range that matters
(`KM_CYP` comment: "near-linear, so it does NOT saturate while phase II
does; the shunt to reactive metabolite therefore rises with dose ON ITS
OWN") — this asymmetric saturation between phase-II conjugation and CYP
bioactivation is the entire mechanism behind the model's dose-threshold
behaviour. Collapsing `EFFECT_APAP` into a bounded `Emax*Xᵞ/(EC50ᵞ+Xᵞ)`
form would materially change what the model computes, not just rename it —
so per the guide's explicit allowance ("if a compound's structure genuinely
doesn't fit any of the guide's 4 archetypes, say so and handle it as
bespoke"), `EFFECT_APAP` is a straight rename of the original's own
`prod_RM = v_cyp / VLIV` (a production *flux*, uM/h, not a fractional
occupancy) — same Michaelis-Menten-on-`CU_APAP` shape the original already
had, values unchanged.

A second, genuinely separate effect exists: the original's `v_bsep` term
includes a competitive-inhibition factor, `1.0 / (1.0 + CuP/KI_BSEP)`. This
*is* extractable as its own named quantity without changing its shape
(same category as the linear competitive-antagonism term kept for SPI in
`acne-vulgaris/acn_refactor_notes.md`) — promoted to `EFFECT_APAP_BSEP`,
a straight rename, still not Hill-shaped (linear/competitive, not
Emax/EC50). By default (`KI_APAP_BSEP = 1e9`) this term is inert
(`EFFECT_APAP_BSEP ≈ 1`); the file's own `DRUG_B` archetype
(`KI_BSEP = 0.5` override) is what activates it, representing a real
BSEP-inhibiting hepatotoxin (troglitazone/bosentan-class) reusing the same
compartments with different parameter values — this parameter-swap pattern
(same DSL block, different `param()` overrides) is how the original file
itself represents multiple real drugs without multiple PK blocks; it is
unchanged by this refactor.

### NAC — Archetype 1 (no depot, single compartment, linear elimination)

```
dxdt_NACC = -CLNAC/VNAC * NACC;      // original
dxdt_CENT_NAC = -CL_NAC/V1_NAC * CENT_NAC;   // refactored (rename only)
```

Exact Archetype-1 match, dosed by IV infusion (the file's own 3-phase
Prescott regimen: 150 mg/kg/1h, 50 mg/kg/4h, 100 mg/kg/16h), not an
`ev()`/depot pattern — the guide's "continuous infusion input" edge case,
preserved unchanged (same `rate`/`amt` dosing-frame construction as the
original's `nac_prescott()`).

Renamed: `NACC`→`CENT_NAC`, `VNAC`→`V1_NAC`, `CLNAC`→`CL_NAC`. `EMAX_NAC`
and `EC50_NAC` already matched the convention (kept as-is, values
unchanged). Its concentration was inside a `#define` macro (`CNAC =
NACC/VNAC`), per the census's own classification — pulled out into a named
`$ODE`-local `double C_NAC = CENT_NAC / V1_NAC;`, as the task instructed.

`EFFECT_NAC` promoted from the original's inline `enac` variable:

```
double enac = EMAX_NAC * CNAC / (EC50_NAC + CNAC);                                    // original
double EFFECT_NAC = EMAX_NAC * pow(C_NAC, GAMMA_NAC)
                     / (pow(EC50_NAC, GAMMA_NAC) + pow(C_NAC, GAMMA_NAC));              // refactored
```

`GAMMA_NAC = 1.0` is new and explicit (the original had no Hill exponent —
this is algebraically identical since `pow(x,1)=x`, a rename not a fit).
`C_NAC` also drives a second, disclosed non-Hill term unchanged from the
original — direct ROS scavenging, `KROS_NAC * C_NAC * ROS` — the guide
permits one exposed concentration feeding more than one disease equation
as long as it isn't combining *multiple drugs'* effects into one term,
which this isn't.

## Upstream build defect found, fixed syntax-only, and disclosed

`POST /model_manifest` on the untouched original's own extracted DSL
(qspserver `mrgsolve_api`, `http://localhost:8007`) fails outright:

```
Error: improper annotation format
 input: : saturate while phase II does; the shunt to reactive
 context: parse annotated parameter block (PARAM)
Execution halted
```

Three `$PARAM @annotated` entries (`KM_CYP`, `STER_TAU`, `KBIL_ALT`) write
their description across three lines using a leading-`:` continuation
idiom, e.g.:

```
KM_CYP   : 6000.0 : Km CYP bioactivation (uM) — near-linear, so it does NOT
                  : saturate while phase II does; the shunt to reactive
                  : metabolite therefore rises with dose ON ITS OWN
```

mrgsolve 2.0.1's `@annotated` parser only accepts one `name : value :
description` triple per line and rejects the continuation. This is
unrelated to either compound's own PK block — all three parameters
(`KM_CYP`, `STER_TAU`, `KBIL_ALT`) are disease-side physiology (CYP
bioactivation kinetics, steroid-onset time constant, bilirubin's
non-hepatic elimination), not APAP's or NAC's PK.

**Not fixed in the original**, per the never-edit-upstream rule. Logged as
`translations/UPSTREAM_ISSUES.md` **#113**. Fixed *syntax-only* in
`dili_mrgsolve_model_refactored.R` (each 3-line description collapsed into
one line — no numeric or behavioral change), per the guide's settled
policy for a build-incompatible original, so the delivered refactor
actually compiles and runs through the qspserver API. Verification below
uses a **patched-original** baseline (the untouched original's own DSL,
with only this same 3-line-collapse applied and nothing else) as the
comparison target, confirming the fix is syntax-only.

## Verification (qspserver `mrgsolve_api`, `http://localhost:8007`)

Two of the original file's own scenarios (`scenarios()` in the R wrapper),
run through both the patched-original baseline and the refactored DSL via
`/run_simulation`, comparing every shared `$CAPTURE`/compartment output
across the full time grid:

**S5 "APAP 350 mg/kg + NAC at 8 h"** (single oral APAP dose at t=0 +
`nac_prescott(start=8)`, `end=336h`/`delta=0.25`, 1349 points): every
shared output matches to floating-point-solver noise —
`AGUT`/`AC`/`AP`/`ALIV`/`NACC` vs. their renamed counterparts max abs diff
1e-9 to 1e-16; the full disease cascade (`RM`, `GSH`, `CYS`, `ADD`, `ROS`,
`NRF2`, `MITO`, `ATP`, `JNK`, `BAH`, `BAP`, `HEP`, `HEPS`, `NECR`, `DAMP`,
`KC`, `TNF`, `IL10`, `TCELL`, `TREG`, `ALT`, `AST`, `TBIL`, `ALP`, `FV`,
`MIR`, `INR`, `RRATIO`, `HYLAW`, `REDOXG`) max abs diff 0 to ~1e-9 for
every variable except `FBIOACT` (see below). `EFFECT_NAC` at t=10h
independently verified against the closed-form Hill expression
(`1.6*C_NAC/(150+C_NAC)`) — exact match to displayed precision.

**S10 "Drug B (BSEP inhibitor) 8 mg/kg q12h x 28 d"** (`DRUG_B` archetype:
`KI_BSEP=0.5` genuinely activates the BSEP-inhibition pathway that S5 left
inert; 56 doses, `end=1400h`/`delta=0.25`, 5657 points, no solver-step
issues): every shared output matches to floating-point-solver noise (max
abs diff <=1e-4, most <=1e-10); peak ALT (140.9309), peak TBIL (0.9906),
and `HYLAW` (false throughout) are identical between original and
refactored.

**`FBIOACT` near-zero/near-zero artifact (disclosed, not a real
mismatch).** At isolated late timepoints in S5 (after ~260h, once the
single APAP dose is essentially fully cleared), `CU`/`CU_APAP` are
themselves at ODE-solver floating-point-noise level (~1e-13 to ~1e-14).
`FBIOACT`'s formula divides a numerator built from `CU` by a denominator
built from the same `CU` plus a `1e-12` floor constant — when the true
value of `CU` is comparable to that `1e-12` floor, which side of zero the
solver's own noise happens to land on dominates the ratio, producing
absolute swings up to ~4.6 between the two runs at those specific
timepoints even though the underlying state (`CU`/`CU_APAP`) differs by
~1e-13. Same category of artifact as documented for `HEMIN_PLASMA`/
`CENT_HEM` in `acute-intermittent-porphyria/aip_refactor_notes.md`. Not
adjusted or hidden — reported here per the guide's "a poor fit is a
finding, not a bug to hide" (though this isn't a fit issue, the same
transparency principle applies to a solver-noise artifact).

**qspserver `/model_manifest` discoverability.** Confirmed all renamed PK
parameters (`KA_APAP, F_APAP, CL_APAP, V1_APAP, V2_APAP, Q_APAP, KP_APAP,
FU_APAP, KI_APAP_BSEP, V1_NAC, CL_NAC, EMAX_NAC, EC50_NAC, GAMMA_NAC`)
appear in the manifest's `parameters` (134 total, fixed values unchanged
from the original). `C_APAP, CU_APAP, EFFECT_APAP, EFFECT_APAP_BSEP,
C_NAC, EFFECT_NAC` are state-dependent (derived from compartments in
`$ODE`), so per this fork's established precedent
(`acute-pancreatitis/ap_refactor_notes.md`,
`acute-intermittent-porphyria/aip_refactor_notes.md`) they cannot be
`$PARAM` defaults — confirmed discoverable instead in the manifest's
`outputPaths` (50 total: 5 renamed PK compartments + 28 disease
compartments + 11 original `$TABLE`/`$CAPTURE` outputs + the 6 new
`C_<STEM>`/`EFFECT_<STEM>` quantities).

No separate `.cpp` extraction was left in the repository — the deliverable
list is just the `_refactored.R` and `_refactor_notes.md` siblings plus the
census update. The quoted DSL block was extracted to a scratch file for
every API call above, confirmed byte-identical to the `dili_code <- '...'`
string in the delivered `dili_mrgsolve_model_refactored.R` (diffed
directly, `IDENTICAL`), then the scratch copy was deleted.

## R syntax check

`Rscript -e 'cat(length(parse("dili_mrgsolve_model_refactored.R")))'`
succeeds (24 top-level expressions parsed) — this is an R-syntax check
only (the DSL itself is a quoted string, opaque to R's parser), not a
local mrgsolve build, per the task's "do not use local mrgsolve/R for
verification" instruction; all mrgsolve-level verification above was done
exclusively through the qspserver API.
