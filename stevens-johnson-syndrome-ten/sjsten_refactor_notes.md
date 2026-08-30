# Refactor notes — `stevens-johnson-syndrome-ten/sjsten_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
(`FORK_WORKFLOW_GUIDE.md`, Part 2), all **seven** rows this file has in
`driver-patches/data/compound_perturbation_census.md` (all classified
"Redirect concentration (clean single site)") were rewritten: **Culprit
drug (DRUG)**, **Cyclosporine (CSA)**, **Etanercept (ETAN)**,
**Infliximab (INFL)**, **IVIG**, **JAK inhibitor (JAKI)**, and
**Prednisolone/methylprednisolone (PRED)**. Every disease-side equation
(antigen presentation kinetics, T-cell activation, cytokine/effector
pools, keratinocyte apoptosis, BSA detachment, SCORTEN, mortality hazard,
re-epithelialisation) is untouched, as is the therapeutic-plasma-exchange
(`SCEN_TPE`/`EFF_TPE_GNLY`) mechanism, which has no PK compartment and is
not one of this file's seven census rows.

## "Culprit drug" is not a mislabeled disease-state variable — it is a real, dosed PK compound

The task instructions specifically flagged this row for a sanity check:
is "Culprit drug" a generic placeholder/trigger with no real PK, or a
genuinely externally-dosed compound? Checking the code:

- `A_drug`/`A_drug_dep` (renamed `CENT_DRUG`/`GUT_DRUG`) are a genuine,
  independently-dosed one-compartment-with-depot PK pair: first-order
  absorption (`KA_drug`), a central volume (`V_drug`) and linear clearance
  (`CL_drug`), dosed via `ev(amt=400, ii=12, addl=2, cmt="A_drug_dep")` in
  the file's own trailing R driver example. This is a real dosing route,
  not an initial condition or a pure flag — the depot is genuinely
  drivable by an external covariate/dataset exactly like any other
  compound in this corpus.
- The file's own header comment and README (`stevens-johnson-syndrome-ten/README.md`,
  "Culprit drugs & PK" cluster) describe it as **illustrative**: "PK
  (illustrative: carbamazepine-like)", standing in generically for the
  class of real culprit drugs (carbamazepine, allopurinol,
  sulfamethoxazole, lamotrigine, etc.) rather than any one specific named
  compound. That is a real, if deliberately generic, modelling choice —
  unlike `neonatal-hyperbilirubinemia`'s "Bilirubin (SNMP)" or
  `prostate-cancer`'s "AR Signaling (DEG)" rows (both process-description
  names paired with the wrong compound identity), the census's "Culprit
  drug" label here is **accurate**: it names exactly what the compartment
  represents. No correction to the census's compound name was needed
  beyond adding the `(DRUG)` stem tag and this clarifying note.

**Conclusion: a real, drivable PK compound — refactored as such, stem `DRUG`.**

## Archetype per compound

**Culprit drug (DRUG): Archetype 3 minus peripheral (depot + central,
linear).** `GUT_DRUG`/`CENT_DRUG`, first-order absorption
(`KA_DRUG`), linear elimination (`CL_DRUG`/`V1_DRUG`), no peripheral
compartment. `C_DRUG = CENT_DRUG / V1_DRUG` — a genuine amount/volume
division, matching the original's own `Cdrug = A_drug / V_drug` exactly
(not the bespoke identity pattern below).

**Cyclosporine (CSA), Etanercept (ETAN), JAK inhibitor (JAKI): Archetype
3 minus peripheral (depot + central, linear), with one bespoke deviation
shared by all six treatment compounds — see next section.** Each has a
depot (`CSA_dep`/`ETAN_dep`/`JAKI_dep`, renamed `GUT_CSA`/`GUT_ETAN`/
`GUT_JAKI`) and a central compartment, first-order absorption, linear
elimination, no peripheral compartment.

**Infliximab (INFL), Prednisolone/methylprednisolone (PRED): Archetype 1
(single compartment, linear elimination, no depot), same bespoke
deviation.** Both are dosed directly into their own central compartment
(IV route) — no absorption compartment in the original.

**IVIG: Archetype 1 (single compartment, no depot) plus the guide's
"continuous infusion input, not an event" edge case, same bespoke
deviation.** `r_ivig = (SCEN_IVIG>0.5 && SOLVERTIME<DAYS_IVIG) ?
DOSE_IVIG*WT/V1_IVIG : 0.0` is a zero-order rate computed directly inside
`dxdt_CENT_IVIG`, gated by the model's own `SCEN_IVIG`/`DAYS_IVIG`
parameters rather than an `ev()`/dataset event — preserved exactly, per
the guide's explicit allowance for this pattern.

## The shared bespoke deviation: `C_<STEM>` is an identity of the compartment, not `CENT/V1`, for all six treatment compounds

The guide's Archetype 1/3 templates define `C_<STEM> = CENT_<STEM> /
V1_<STEM>`. This file's own effect terms for **all six treatment
compounds** (not the culprit drug) read the raw compartment value
directly, with no division by volume, even though `CL/V` is used
correctly as the elimination rate constant in each `dxdt_`:

```
double f_pred = 1.0 - (SCEN_PRED > 0.5 ? EFF_PRED_T * (PRED / (PRED + 0.05)) : 0.0);  // original
...
dxdt_PRED = - (CL_PRED/V_PRED)*PRED;                                                  // original
```

`PRED` (and `CSA`, `ETAN`, `INFL`, `IVIG`, `JAKI`) is read bare in every
one of the six effect ratios — never `PRED/V_PRED`. Since `dA/dt =
-(CL/V)*A` is a valid ODE regardless of whether the state is *labelled*
"amount" or "concentration" (the coefficient `CL/V` is just a rate
constant either way), and the file's own header states "concentrations =
mg/L", the most faithful reading is that the original author intended
these six central compartments to hold concentration directly (the same
convention used for `CSNMP`/`CPB`/`CUDCA` in
`neonatal-hyperbilirubinemia`), not amount. Converting to a literal
`CENT/V1` division would **change** the original's numeric behaviour
(rescaling each Hill term's effective magnitude by `V1_<STEM>`), not
merely rename it — which the guide explicitly rules out ("never flatten
... because it's inconvenient").

**Handled exactly as the `neonatal-hyperbilirubinemia` precedent
(`nhb_refactor_notes.md`):** all six central compartments
(`CENT_CSA`/`CENT_ETAN`/`CENT_INFL`/`CENT_IVIG`/`CENT_JAKI`/`CENT_PRED`)
are kept as concentration-integrating states, renamed only, and each
`C_<STEM>` is a plain identity of the compartment
(`double C_CSA = CENT_CSA;`, etc.), not a division. Disclosed here as a
bespoke deviation from the literal Archetype-1/3 template, per the
guide's explicit allowance ("if a compound's PK genuinely doesn't
resemble any archetype ... rename to the convention, isolate it from PD,
expose `C_<STEM>`, and note ... why"). The culprit drug (`DRUG`) is the
one exception — its own `Cdrug = A_drug/V_drug` in the original already
uses the standard division, so `C_DRUG = CENT_DRUG/V1_DRUG` is kept as a
genuine division, not an identity.

## The culprit drug's own effect term is bespoke: linear mass-action, not Hill

`Ag_HLA` (antigen-HLA complex formation, disease state, unchanged) is
driven by `KON_HLA*hla_amp*Cdrug` — a **linear** mass-action term in the
drug concentration, not a saturating `Emax*C/(EC50+C)` ratio. Per the
guide ("if a compound's structure genuinely doesn't fit any of the 4
archetypes, say so and handle it as bespoke"), no `EMAX_DRUG`/
`EC50_DRUG`/`GAMMA_DRUG` were invented — that would misrepresent a
linear driver as a saturable one, changing its shape at high
concentration. Instead, `EFFECT_DRUG = hla_amp * C_DRUG` is named as the
"one named function of concentration" the guide's naming section calls
for (satisfying the general "expose one named site" principle), while
`dxdt_Ag_HLA = KON_HLA*EFFECT_DRUG - KOFF_HLA*Ag_HLA` keeps `KON_HLA`
external to the named term, exactly reproducing the original's algebra.

## The other six compounds: rename, not a refit — all already Hill-shaped

Each of the six treatment compounds' original effect ratio was already
exactly `Emax*C/(EC50+C)` (implicit `gamma=1`, no explicit Hill exponent
in the original for any of them):

| Compound | Original inline term | Renamed |
|---|---|---|
| CSA  | `EFF_CSA_T * (CSA / (CSA + 0.05))`   | `EMAX_CSA=0.6`, `EC50_CSA=0.05`, `GAMMA_CSA=1` |
| PRED | `EFF_PRED_T * (PRED / (PRED + 0.05))`| `EMAX_PRED=0.5`, `EC50_PRED=0.05`, `GAMMA_PRED=1` |
| ETAN | `EFF_ETAN_TNF * (ETAN / (ETAN + 0.3))`| `EMAX_ETAN=0.8`, `EC50_ETAN=0.3`, `GAMMA_ETAN=1` |
| INFL | `EFF_INFL_TNF * (INFL / (INFL + 0.5))`| `EMAX_INFL=0.85`, `EC50_INFL=0.5`, `GAMMA_INFL=1` |
| IVIG | `EFF_IVIG_FASL * (IVIG / (IVIG + 0.5))`| `EMAX_IVIG=0.5`, `EC50_IVIG=0.5`, `GAMMA_IVIG=1` |
| JAKI | `EFF_JAKI_IFN * (JAKI / (JAKI + 0.1))`| `EMAX_JAKI=0.7`, `EC50_JAKI=0.1`, `GAMMA_JAKI=1` |

Each `EFFECT_<STEM> = EMAX_<STEM>*pow(C,GAMMA)/(pow(EC50,GAMMA)+pow(C,GAMMA))`
is a pure rename — no `nls()` fit performed, since none of the six needed
one. The original's `SCEN_<STEM> > 0.5 ? ... : 0.0` gating ternary is kept
exactly, wrapping `EFFECT_<STEM>` rather than being absorbed into it, so
the master on/off switch behaves identically to the original in every
case (including the edge case of a nonzero compartment concentration with
its scenario flag left off, which never arises under the file's own
scenarios but would behave identically to the original either way).

## Naming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `A_drug` (cmt) | `CENT_DRUG` | -- | central (amount) |
| `A_drug_dep` (cmt) | `GUT_DRUG` | -- | depot |
| `CL_drug` | `CL_DRUG` | 4.0 L/h | clearance |
| `V_drug` | `V1_DRUG` | 70 L | central volume |
| `KA_drug` | `KA_DRUG` | 0.5 1/h | absorption rate |
| `Cdrug` (local) | `C_DRUG` | -- | **the exposed concentration** (genuine amount/V) |
| -- (inline in `dxdt_Ag_HLA`) | `EFFECT_DRUG` (new name) | -- | bespoke linear driver, `hla_amp*C_DRUG` |
| `CSA` (cmt) | `CENT_CSA` | -- | central, concentration-state (bespoke) |
| `CSA_dep` (cmt) | `GUT_CSA` | -- | depot |
| `KA_CSA`, `CL_CSA`, `V_CSA` | `KA_CSA`, `CL_CSA`, `V1_CSA` | 1.4 /h, 5 L/h, 100 L | PK |
| `EFF_CSA_T` | `EMAX_CSA` | 0.6 | Hill Emax |
| -- (inline literal `0.05`) | `EC50_CSA` (new) | 0.05 | Hill EC50 |
| -- (none) | `GAMMA_CSA` (new) | 1.0 | Hill exponent |
| `CSA` (local, identity) | `C_CSA` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_csa`) | `EFFECT_CSA` (new name) | -- | T-activation inhibition |
| `ETAN` (cmt) | `CENT_ETAN` | -- | central, concentration-state (bespoke) |
| `ETAN_dep` (cmt) | `GUT_ETAN` | -- | depot |
| `KA_ETAN`, `CL_ETAN`, `V_ETAN` | `KA_ETAN`, `CL_ETAN`, `V1_ETAN` | 0.4 /d, 0.07 L/d, 7.6 L | PK |
| `EFF_ETAN_TNF` | `EMAX_ETAN` | 0.8 | Hill Emax |
| -- (inline literal `0.3`) | `EC50_ETAN` (new) | 0.3 | Hill EC50 |
| -- (none) | `GAMMA_ETAN` (new) | 1.0 | Hill exponent |
| `ETAN` (local, identity) | `C_ETAN` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_tnf_etan`) | `EFFECT_ETAN` (new name) | -- | TNF-alpha neutralisation |
| `INFL` (cmt) | `CENT_INFL` | -- | central, concentration-state (bespoke), no depot |
| `CL_INFL`, `V_INFL` | `CL_INFL`, `V1_INFL` | 0.27 L/d, 4 L | PK |
| `EFF_INFL_TNF` | `EMAX_INFL` | 0.85 | Hill Emax |
| -- (inline literal `0.5`) | `EC50_INFL` (new) | 0.5 | Hill EC50 |
| -- (none) | `GAMMA_INFL` (new) | 1.0 | Hill exponent |
| `INFL` (local, identity) | `C_INFL` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_tnf_infl`) | `EFFECT_INFL` (new name) | -- | TNF-alpha neutralisation |
| `IVIG` (cmt) | `CENT_IVIG` | -- | central, concentration-state (bespoke), no depot |
| `CL_IVIG`, `V_IVIG` | `CL_IVIG`, `V1_IVIG` | 0.05 L/d, 5 L | PK |
| `EFF_IVIG_FASL` | `EMAX_IVIG` | 0.5 | Hill Emax |
| -- (inline literal `0.5`) | `EC50_IVIG` (new) | 0.5 | Hill EC50 |
| -- (none) | `GAMMA_IVIG` (new) | 1.0 | Hill exponent |
| `IVIG` (local, identity) | `C_IVIG` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_fasl`) | `EFFECT_IVIG` (new name) | -- | sFasL neutralisation |
| `JAKI` (cmt) | `CENT_JAKI` | -- | central, concentration-state (bespoke) |
| `JAKI_dep` (cmt) | `GUT_JAKI` | -- | depot |
| `KA_JAKI`, `CL_JAKI`, `V_JAKI` | `KA_JAKI`, `CL_JAKI`, `V1_JAKI` | 6.0 /h, 22 L/h, 96 L | PK |
| `EFF_JAKI_IFN` | `EMAX_JAKI` | 0.7 | Hill Emax |
| -- (inline literal `0.1`) | `EC50_JAKI` (new) | 0.1 | Hill EC50 |
| -- (none) | `GAMMA_JAKI` (new) | 1.0 | Hill exponent |
| `JAKI` (local, identity) | `C_JAKI` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_ifn`) | `EFFECT_JAKI` (new name) | -- | IFN-gamma inhibition |
| `PRED` (cmt) | `CENT_PRED` | -- | central, concentration-state (bespoke), no depot |
| `CL_PRED`, `V_PRED` | `CL_PRED`, `V1_PRED` | 12 L/h, 50 L | PK |
| `EFF_PRED_T` | `EMAX_PRED` | 0.5 | Hill Emax |
| -- (inline literal `0.05`) | `EC50_PRED` (new) | 0.05 | Hill EC50 |
| -- (none) | `GAMMA_PRED` (new) | 1.0 | Hill exponent |
| `PRED` (local, identity) | `C_PRED` | -- | **the exposed concentration** (bespoke identity) |
| -- (inline in `f_pred`) | `EFFECT_PRED` (new name) | -- | T-activation inhibition |

All parameter *values* are copied verbatim from the original. `SCEN_TPE`/
`EFF_TPE_GNLY` (therapeutic plasma exchange — not a modelled compound) and
all patient covariates (`WT`, `AGE`, `HLA_RISK`, `PREDOSE`, `SCEN_WD`,
`DAY_WD`) are untouched.

## A design pitfall this refactor avoided (learned from the `neonatal-hyperbilirubinemia`/`familial-mediterranean-fever` precedent)

Every `C_<STEM>`/`EFFECT_<STEM>` is declared **exactly once**, as a plain
`double` local inside `[ODE]`, and captured directly in `[CAPTURE]` with
no re-declaration in `[TABLE]` — the same pattern this file's own
`Cdrug`/`SCORTEN`/`PredMort`/`Re_epi` already used successfully. This
sidesteps the `double`-redeclaration-across-blocks collision documented in
`nhb_refactor_notes.md` and `UPSTREAM_ISSUES.md` #48/#51, so no first
draft was needed here — confirmed correct on the first `/model_manifest`
call (see below).

## Verification

**Method.** Both files' bare mrgsolve DSL (this file has no `code <-
'...'` R wrapper — it is itself the DSL from `[PROB]` through
`[CAPTURE]`, with a trailing `/* ... */` C-style comment holding an
example R driver as documentation only) was mechanically extracted
(everything before the trailing `/*` comment, verbatim — 9,134 chars
original, 12,374 chars refactored) and POSTed to the local qspserver
`mrgsolve_api` service at `http://localhost:8007` (`/model_manifest` then
`/run_simulation`), which compiles and runs each DSL block directly with
mrgsolve 2.0.1 server-side — no local R/mrgsolve install used. Requests
were spaced ~2.5s apart and run sequentially (never more than one in
flight), respecting the service's `max_concurrent_jobs: 2` limit.
`POST /run_simulation`'s `dosing` field addresses compartments by
1-based index, not name; compartment order is identical between the two
files (only names changed), confirmed from `/model_manifest`'s
`outputPaths` (27 entries in both: 23 `$CMT` + 4 original / 17 refactored
`$CAPTURE` names, same 1–23 compartment order).

**The original compiles cleanly as-is** — no pre-existing mrgsolve 2.0.1
build defect, confirmed by the first `/model_manifest` call succeeding
before any refactor changes were made. (A separate, non-blocking logic
defect was found in the culprit-drug block and is disclosed below and in
`UPSTREAM_ISSUES.md` #98 — it does not prevent compilation.)

**Scenarios run — the file's own seven, not invented, plus one bespoke
addition for the one compound none of the seven exercises.** All at
`end=30, delta=0.1` (302 points), matching the file's own R driver
example exactly; no shortening needed (well inside the API's default
`maxsteps` budget):

1. **Supportive** (culprit drug only, 400mg q12h x3 doses into
   `GUT_DRUG`/`A_drug_dep`, cmt 2): **exact match, max abs diff 0.0**
   across all 27 shared outputs.
2. **IVIG** (`SCEN_IVIG=1`, internal zero-order infusion, no `ev()`
   needed): **max abs diff 1.0e-4**, appearing at only 5 of 302
   timepoints (t=8.6, 10.0, 10.3, 12.7, 20.7h) on `IVIG`/`CENT_IVIG`
   only — every other output matches exactly (0.0), and the last ~90
   timepoints of the run (t>21h onward, including the entire final
   segment to t=30) all match exactly (0.0) too, i.e. the deviation does
   not persist or grow. The magnitude (1e-4) is exactly the API's
   4-decimal display rounding granularity straddling a boundary (e.g.
   `48.6821` vs `48.682`) — consistent with genuine floating-point-scale
   noise from the `pow(x, 1.0)` calls the Hill rewrite introduces
   (mathematically identical to `x`, not bit-identical to the FPU)
   perturbing the adaptive LSODA stepper's step choice near a few points,
   the same phenomenon documented in `nhb_refactor_notes.md` for
   Stannsoporfin. Disclosed rather than smoothed over, per "a gate
   failure is a question, not a verdict."
3. **Cyclosporine** (`SCEN_CSA=1`, 100mg q0.5d x21 doses into
   `GUT_CSA`/`CSA_dep`, cmt 19, plus the culprit drug): **exact match,
   max abs diff 0.0**.
4. **Etanercept** (`SCEN_ETAN=1`, 25mg q3d x2 doses into
   `GUT_ETAN`/`ETAN_dep`, cmt 16, plus the culprit drug): **exact match,
   max abs diff 0.0**.
5. **Methylpred** (`SCEN_PRED=1`, 60mg q1d x6 doses into `CENT_PRED`/
   `PRED`, cmt 21, plus the culprit drug): **exact match, max abs diff
   0.0** — the refactored run's first attempt hit a transient
   container-level linker error (`bad reloc symbol index`, `ld` failure)
   unrelated to the DSL content; a clean retry a few seconds later
   (confirmed via `/health` that the service was idle) succeeded and
   reproduced exactly, consistent with the guide's warning that this
   service "has crashed under heavy concurrent load before."
6. **CSA_ETAN** (`SCEN_CSA=1, SCEN_ETAN=1`, both compounds' own dosing
   plus the culprit drug): **exact match, max abs diff 0.0**.
7. **JAK_inhib** (`SCEN_JAKI=1`, 5mg q0.5d x21 doses into `GUT_JAKI`/
   `JAKI_dep`, cmt 22, plus the culprit drug): **exact match, max abs
   diff ~1.0e-16** (floating-point noise only).
8. **Infliximab (bespoke, added)** — none of the file's own seven named
   scenarios doses infliximab at all (`SCEN_INFL`/`DOSE_INFL`/
   `EFF_INFL_TNF` are fully declared and correctly wired into `f_tnf_infl`,
   but the shipped R driver's `scenarios` list and `simulate()` helper
   never reference them). To verify INFL's renamed PK/effect block at a
   genuinely nonzero concentration (not just its default off-state),
   one additional scenario was constructed from the model's **own**
   declared `DOSE_INFL=5 mg/kg` (325mg for the file's own 65kg default
   patient), dosed as an IV bolus into `CENT_INFL`/`INFL` (cmt 18, no
   depot, matching how `PRED` is dosed directly in the file's own
   driver), q14d x2, with `SCEN_INFL=1`: **exact match, max abs diff
   0.0**, peak concentration 451.32 mg/L in both. Per the guide's
   verification rule ("not invented ones") this is disclosed as an
   addition beyond the file's own seven, following the precedent set for
   riluzole in `huntingtons-disease` (`UPSTREAM_ISSUES.md` #89), where
   the same situation (a fully-wired compound never exercised by any
   shipped scenario) was handled the same way.

Six of eight scenarios verify bit-for-bit exact (0.0), one verifies to
~1e-16 (pure floating-point noise), and one (IVIG) verifies to a
4-decimal-display-rounding-scale deviation (1e-4) isolated to 5/302
timepoints that self-corrects — consistent with the guide's tolerance
table for Archetypes 1–3 ("pure structural reorganization ... expect a
near-exact match").

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted, renamed DSL:
`CL_DRUG`, `V1_DRUG`, `KA_DRUG`, `KA_CSA`/`CL_CSA`/`V1_CSA`,
`EMAX_CSA`/`EC50_CSA`/`GAMMA_CSA`, `KA_ETAN`/`CL_ETAN`/`V1_ETAN`,
`EMAX_ETAN`/`EC50_ETAN`/`GAMMA_ETAN`, `CL_INFL`/`V1_INFL`,
`EMAX_INFL`/`EC50_INFL`/`GAMMA_INFL`, `CL_IVIG`/`V1_IVIG`,
`EMAX_IVIG`/`EC50_IVIG`/`GAMMA_IVIG`, `KA_JAKI`/`CL_JAKI`/`V1_JAKI`,
`EMAX_JAKI`/`EC50_JAKI`/`GAMMA_JAKI`, `CL_PRED`/`V1_PRED`,
`EMAX_PRED`/`EC50_PRED`/`GAMMA_PRED` all appear in the manifest's
`parameters` with their original numeric defaults. `C_DRUG`, `EFFECT_DRUG`,
`C_CSA`/`EFFECT_CSA`, `C_ETAN`/`EFFECT_ETAN`, `C_INFL`/`EFFECT_INFL`,
`C_IVIG`/`EFFECT_IVIG`, `C_JAKI`/`EFFECT_JAKI`, `C_PRED`/`EFFECT_PRED`
are state-derived (computed in `[ODE]` from compartment concentrations),
so — per the same reasoning established for `EFFECT_NIM`/`EFFECT_ATR`/
`EFFECT_TCZ`/`C_SNMP`/`C_ANA`/etc. in prior refactors in this fork —
they cannot also be `$PARAM` entries; all appear in the manifest's
`outputPaths` via the file's own (extended) `[CAPTURE]` list, confirmed
discoverable. `GUT_DRUG`/`CENT_DRUG`, `GUT_CSA`/`CENT_CSA`,
`GUT_ETAN`/`CENT_ETAN`, `CENT_INFL`, `CENT_IVIG`, `GUT_JAKI`/`CENT_JAKI`,
`CENT_PRED` also appear in `outputPaths` as ordinary compartments
(indices 1–2, 15, 16–23), unchanged order from the original.

No `.cpp` extraction file was left behind — extraction was in-memory
only, used to build the verification requests above and then discarded,
per the workflow guide.

## Anything else flagged

- **Dead-code finding (not fixed, disclosed and logged):** the culprit
  drug's own `[ODE]` block computes `drug_input =
  (PREDOSE>0.5 && (SCEN_WD<0.5 || SOLVERTIME<DAY_WD)) ? 1.0 : 0.0` but
  never uses it anywhere else in the file — `PREDOSE`, `SCEN_WD`, and
  `DAY_WD` have no effect on the simulated trajectory regardless of their
  values (confirmed by re-running the Supportive scenario with several
  different `PREDOSE`/`DAY_WD` overrides: bit-identical output every
  time). This is the same defect class as `UPSTREAM_ISSUES.md` #89
  (Huntington's riluzole's unused `Emax_riluzole`). Preserved verbatim
  (with an inline comment marking it as known-dead) rather than removed,
  since deleting it would be a silent behavioural-looking edit disguised
  as cleanup; logged as `UPSTREAM_ISSUES.md` #98.
- No compound other than the seven census rows above was touched.
  `SCEN_TPE`/`EFF_TPE_GNLY` (therapeutic plasma exchange, no PK
  compartment) and all disease-side compartments/equations are
  byte-identical to the original.
- The trailing `/* ... */` R driver example comment (documentation only,
  not executed as part of the DSL) was updated to reference the renamed
  compartments (`GUT_DRUG`, `CENT_IVIG`, `GUT_ETAN`, `GUT_CSA`,
  `CENT_PRED`, `GUT_JAKI`) so a reader copying it out still works against
  the refactored model; no other R-side code exists in this file to
  update (unlike files with a separate executable scenario list, this
  file's own example driver is inside a comment).
- All seven of the original's own named scenarios (`Supportive`, `IVIG`,
  `Cyclosporine`, `Etanercept`, `Methylpred`, `CSA_ETAN`, `JAK_inhib`)
  still reference the model the same way; none needed changes beyond the
  `cmt=` renames in the example driver comment above.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
seven `stevens-johnson-syndrome-ten` rows — Target/Pathway filled in for
all seven (previously `?` for six of them), and the "Culprit drug" row
given a `(DRUG)` stem tag and a clarifying note that, unlike some other
files' census rows, this name was already accurate (see above).
