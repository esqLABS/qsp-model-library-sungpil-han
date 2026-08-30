# Refactor notes — `prurigo-nodularis/pn_mrgsolve_model.R`

Compounds refactored: **Dupilumab (DUPI)**, **Gabapentin (GABA)**, the oral
**JAK1 inhibitor / abrocitinib (JAKI)**, **Nalbuphine ER (NAL)**, and
**Nemolizumab (NEMO)** — the five rows in
`driver-patches/data/compound_perturbation_census.md` classified
`prurigo-nodularis | DUPI/GABA/JAKI/NAL/NEMO`, all "Redirect concentration
inside #define macro". Every disease-side compartment and equation
(`BARRIER`, `TSLP`, `TH2`, `TH17`, `IL31`, `PSENS`, `IENFD`, `CSENS`,
`OPIOID`, `SCRATCH`, `NODULE`, `WINRS`, `IGA`, `SLEEP`, `DLQI`,
`CUM_SCRATCH`) is byte-for-byte identical to `pn_mrgsolve_model.R`; only the
five compounds' own PK compartments, exposed concentration, and Hill effect
terms were renamed and pulled out of `#define` macros.

## Why "inside #define macro" needed more than a rename

All five compounds' concentrations were computed as `$GLOBAL`
C-preprocessor macros (`#define CP_DUPI (C_DUPI/V_DUPI)`, `#define CP_NEMO
(C_NEMO/V_NEMO)`, `#define CP_GABA (C_GABA/V_GABA)`, `#define CP_NAL
(C_NAL/V2_NAL*1000.0)`, `#define CP_JAKI (C_JAKI/V_JAKI*1000.0)`), each
feeding its own effect macro (`#define EFF_DUPI (DUPI_EMAX*CP_DUPI/(CP_DUPI+
DUPI_EC50))`, etc.), which were then read directly inside `$ODE`'s disease
equations. A second complication specific to this file: the original's
**central-compartment amount** for every one of these five drugs is itself
named `C_<STEM>` (`C_DUPI`, `C_NEMO`, `C_GABA`, `C_NAL`, `C_JAKI` are `$CMT`
states, not concentrations) — the exact opposite of this guide's convention,
where `C_<STEM>` must be the exposed concentration and `CENT_<STEM>` the
amount. So this refactor needed two renames per compound, not one: the old
amount compartment `C_<STEM>` &rarr; `CENT_<STEM>`, and the old macro
`CP_<STEM>` &rarr; the new `C_<STEM>` (now a real, non-macro, state-dependent
quantity). `P_NAL` (peripheral amount) was renamed to `PERI_NAL` for the
same reason (`PERI_<STEM>` is the guide's convention).

Per the guide's qspserver-compatibility requirement 4, `C_<STEM>` and
`EFFECT_<STEM>` (replacing `CP_<STEM>`/`EFF_<STEM>`) must be `$CAPTURE`d so
they are discoverable outputs — which is incompatible with keeping them as
macros (a `$CAPTURE`d name that is also a `#define` risks the preprocessor
substituting mrgsolve's own auto-generated capture code, per the mechanism
documented in `hypoparathyroidism/hypopt_refactor_notes.md` and
`distal-renal-tubular-acidosis/drta_refactor_notes.md`). Resolution: `C_DUPI`,
`EFFECT_DUPI`, `C_NEMO`, `EFFECT_NEMO`, `C_GABA`, `EFFECT_GABA`, `C_NAL`,
`EFFECT_NAL`, `C_JAKI`, `EFFECT_JAKI` were pulled out of `$GLOBAL` entirely
(the `#define`s for all ten are deleted) and are instead computed as
ordinary `double` locals inside `$ODE`, immediately after each compound's own
PK block and before any disease equation reads them (each compound's PK
block already precedes the disease block that consumes its effect, so no
reordering of `$ODE` was needed, unlike the hypoparathyroidism case).

**A finding that differs from the `hypopt`/`drta` precedent's stated
mechanism**: those two files' notes say a bare `double NAME = ...;` inside
`$ODE` collides with mrgsolve's own auto-declaration when `NAME` is also
`$CAPTURE`d ("any bare `double NAME = expr;` inside `$ODE` is auto-promoted
by mrgsolve to a class member, and `$CAPTURE` naming that same identifier
requires mrgsolve to *also* declare a member of that exact name"), so both
files omit `double` on their `$ODE`-scoped, `$CAPTURE`d names. Tested
empirically here via `/model_manifest` (see Verification): **omitting
`double`** on `C_DUPI`/`EFFECT_DUPI`/etc. failed to compile (`'C_DUPI' was
not declared in this scope`) — the opposite problem, no auto-declaration
happened at all in this file's `$ODE` block. **Keeping the explicit
`double` prefix** on every one of these ten assignments compiled cleanly
with no duplicate-declaration error, and correctly appears in
`/model_manifest`'s `outputPaths`. So for this file, the ordinary
`$ODE`-scoped-local + `$CAPTURE` combination needed **no special-casing at
all** — a plain `double NAME = expr;` followed by listing `NAME` in
`$CAPTURE` is exactly what standard mrgsolve usage elsewhere in this same
file already does for `WINRS_PCT_IMPROVE`/`IGA_SUCCESS`/`RESPONDER` in
`$TABLE`. This is disclosed here because it directly contradicts the
generalization in the two cited precedents; the collision they describe is
real (confirmed independently in their own files) but is evidently not
universal across every `$ODE`/`$CAPTURE` combination in this build — worth a
flag for whoever standardizes this pattern further.

## Archetype

All five compounds fit **archetype 3, without a peripheral compartment**
(depot + central, linear elimination, first-order absorption, with
bioavailability multiplied directly into the absorption term inside the
ODE — exactly the guide's archetype-3 template, not the alternative
mrgsolve-`F_<CMT>`-dosing mechanism), **except nalbuphine ER, which is a
genuine archetype 3** (depot + central + peripheral) — the original already
modeled a real peripheral compartment for nalbuphine (`P_NAL`, now
`PERI_NAL`) via `Q`/`V2`/`V3`, so that structure is kept unchanged.

- **Dupilumab**: `GUT_DUPI`/`CENT_DUPI` (was `DEPOT_DUPI`/`C_DUPI`),
  `KA_DUPI`, `F_DUPI`, `CL_DUPI` unchanged; `V1_DUPI` (was `V_DUPI`, value
  4.6 L unchanged — renamed to match the guide's archetype-1/3 template,
  which uses `V1_<STEM>` even with no peripheral compartment).
- **Nemolizumab**: `GUT_NEMO`/`CENT_NEMO` (was `DEPOT_NEMO`/`C_NEMO`),
  `KA_NEMO`, `F_NEMO`, `CL_NEMO` unchanged; `V1_NEMO` (was `V_NEMO`, 4.6 L
  unchanged).
- **Gabapentin**: `GUT_GABA`/`CENT_GABA` (`GUT_GABA` unchanged, `C_GABA`
  &rarr; `CENT_GABA`), `KA_GABA`, `F_GABA`, `CL_GABA` unchanged; `V1_GABA`
  (was `V_GABA`, 58 L unchanged).
- **Oral JAK1 inhibitor**: `GUT_JAKI`/`CENT_JAKI` (`GUT_JAKI` unchanged,
  `C_JAKI` &rarr; `CENT_JAKI`), `KA_JAKI`, `F_JAKI`, `CL_JAKI` unchanged;
  `V1_JAKI` (was `V_JAKI`, 100 L unchanged). The scenario switch `JAKI_ON`
  (1 = drug pharmacologically active this run, concentration still
  simulated either way) is unchanged — it is a scenario flag, not a
  PK/Hill parameter covered by this guide's naming table, so it keeps its
  original name.
- **Nalbuphine ER**: `GUT_NAL` unchanged, `C_NAL` &rarr; `CENT_NAL`, `P_NAL`
  &rarr; `PERI_NAL`, `KA_NAL`, `CL_NAL`, `Q_NAL`, `F_NAL` unchanged. The
  original's volume naming was inverted relative to the guide's convention
  — `V2_NAL` was the **central** volume and `V3_NAL` the **peripheral**
  volume (guide: `V1_<STEM>`=central, `V2_<STEM>`=peripheral). Renamed to
  match: `V1_NAL` (was `V2_NAL`, 180 L unchanged) is now central,
  `V2_NAL` (was `V3_NAL`, 260 L unchanged) is now peripheral. Values and
  equations unchanged — this is a pure rename of which number means which
  volume, not a reparameterization.

No compound needed a bespoke structure; the archetype-3 (minus peripheral)
template already matched what the original did for four of the five, and
genuine archetype 3 matched nalbuphine exactly.

## TMDD check (dupilumab, nemolizumab)

Both dupilumab (anti-IL-4Rα) and nemolizumab (anti-IL-31RA) are monoclonal
antibodies against a cell-surface receptor, the class this guide flags to
check for target-mediated drug disposition. **Neither is modeled as TMDD in
this file.** There is no free-receptor or drug-receptor-complex compartment
for either drug anywhere in `pn_mrgsolve_model.R` — each is a single-
compartment linear PK model (`CL`/`V1`/`KA`/`F`), and each drug's receptor
blockade is captured purely as an algebraic Hill fraction of its own plasma
concentration (`EFF_DUPI`/`EFF_NEMO`, now `EFFECT_DUPI`/`EFFECT_NEMO`), never
as an ODE-solved binding system. So both are archetype 3 (no peripheral),
straight PK, with a rename-only Hill interface — not a curve-fit situation.
This is disclosed rather than silently treated as "obviously no TMDD",
since the guide explicitly asks to check.

## Hill interface: rename, not a fit, for all five

Every one of the five original effect macros was already exactly
`Emax·C/(C+EC50)` — an implicit-`gamma=1` Hill ratio — so per the guide's
rename rule, `EMAX_<STEM>`/`EC50_<STEM>`/`GAMMA_<STEM>=1` were pulled out as
explicit named parameters with the original's own values, and the macro
became a named `EFFECT_<STEM>` computed with `pow(...)`:

```
EFFECT_DUPI = EMAX_DUPI*pow(C_DUPI,GAMMA_DUPI)/(pow(EC50_DUPI,GAMMA_DUPI)+pow(C_DUPI,GAMMA_DUPI));
EFFECT_NEMO = EMAX_NEMO*pow(C_NEMO,GAMMA_NEMO)/(pow(EC50_NEMO,GAMMA_NEMO)+pow(C_NEMO,GAMMA_NEMO));
EFFECT_GABA = EMAX_GABA*pow(C_GABA,GAMMA_GABA)/(pow(EC50_GABA,GAMMA_GABA)+pow(C_GABA,GAMMA_GABA));
EFFECT_NAL  = EMAX_NAL*pow(C_NAL,GAMMA_NAL)/(pow(EC50_NAL,GAMMA_NAL)+pow(C_NAL,GAMMA_NAL));
EFFECT_JAKI = (JAKI_ON>0.5) ? (EMAX_JAKI*pow(C_JAKI,GAMMA_JAKI)/(pow(EC50_JAKI,GAMMA_JAKI)+pow(C_JAKI,GAMMA_JAKI))) : 0.0;
```

`EFFECT_JAKI`'s `JAKI_ON` gate is preserved verbatim from the original's
`EFF_JAKI` macro (concentration is still simulated regardless of the flag;
only the disease-facing effect is gated). Original parameter names
(`DUPI_EC50`/`DUPI_EMAX`, `NEMO_EC50`/`NEMO_EMAX`, `GABA_EC50`/`GABA_EMAX`,
`JAKI_EC50`/`JAKI_EMAX`, `NAL_EC50`/`NAL_EMAX`) were reordered to
`EC50_<STEM>`/`EMAX_<STEM>` per the guide's convention; all five numeric
values are unchanged. No fit was needed or performed for any compound.

Each `EFFECT_<STEM>` replaces the old `EFF_<STEM>` macro verbatim at every
one of its call sites in the disease `$ODE` block (`TH2`, `IL31`, `PSENS`,
`CSENS`, `OPIOID` tone correction, `NODULE` resolution) — no disease-side
coefficient (e.g. `2.0`/`1.5` multipliers on `EFFECT_DUPI`/`EFFECT_NEMO` in
`nodule_resolve`) was touched, since those are disease-model weighting
constants, not part of the compound's own PK/effect interface.

## qspserver `/model_manifest` discoverability

Confirmed via `POST /model_manifest` against the extracted DSL
(`http://localhost:8007`):

- `KA_DUPI, F_DUPI, CL_DUPI, V1_DUPI, EMAX_DUPI, EC50_DUPI, GAMMA_DUPI`;
  `KA_NEMO, F_NEMO, CL_NEMO, V1_NEMO, EMAX_NEMO, EC50_NEMO, GAMMA_NEMO`;
  `KA_GABA, F_GABA, CL_GABA, V1_GABA, EMAX_GABA, EC50_GABA, GAMMA_GABA`;
  `KA_NAL, CL_NAL, V1_NAL, V2_NAL, Q_NAL, F_NAL, EMAX_NAL, EC50_NAL,
  GAMMA_NAL`; `KA_JAKI, CL_JAKI, V1_JAKI, F_JAKI, EMAX_JAKI, EC50_JAKI,
  GAMMA_JAKI, JAKI_ON` — all 38 appear in the manifest's `parameters`
  (fixed values, unchanged from the original).
- `C_DUPI, EFFECT_DUPI, C_NEMO, EFFECT_NEMO, C_GABA, EFFECT_GABA, C_NAL,
  EFFECT_NAL, C_JAKI, EFFECT_JAKI` are state-dependent (derived from
  `CENT_<STEM>`), so they cannot be `$PARAM` defaults; all ten appear in the
  manifest's `outputPaths` (via `$CAPTURE`), confirmed discoverable.
- No separate `.cpp` extraction was left in the repository (the task's
  deliverable list is just the `_refactored.R` and `_refactor_notes.md`
  siblings plus the census update); the quoted DSL block was extracted
  verbatim to a scratch file for every API call below and confirmed
  byte-identical to the `pn_code <- '...'` string in the delivered
  `pn_mrgsolve_model_refactored.R` (diffed directly, `IDENTICAL`), then the
  scratch copy was deleted after verification.

## Verification

Per the guide's mandatory protocol, via the qspserver `mrgsolve_api`
service (`POST /model_manifest`, `POST /run_simulation`,
`http://localhost:8007`, requests spaced ~2s apart), comparing the untouched
original against `pn_mrgsolve_model_refactored.R`.

**No upstream defect found.** `pn_mrgsolve_model.R` compiled cleanly on the
first `/model_manifest` call, with the `CP_<STEM>` macros captured directly
as-is (its own pre-existing pattern happens not to hit the
macro/auto-declaration collision the guide warns about, since the original
never explicitly re-declares those names elsewhere) — nothing was logged to
`translations/UPSTREAM_ISSUES.md`.

**`/model_manifest`**: both the original and refactored DSL compiled
successfully; the refactored manifest additionally lists all 38 renamed/new
parameters above in `parameters` and all 10 renamed concentration/effect
quantities in `outputPaths`, as required.

**`/run_simulation`, all five compounds**: reproduced the dosing composition
of the original file's own R-side scenario builders (`ev_dupilumab()`,
`ev_nemolizumab()`, `ev_gabapentin()`, `ev_nalbuphine_er()`,
`ev_abrocitinib()`) exactly — same amounts, same intervals, same
compartment (1-based index, unchanged by the rename in every case: DUPI=1,
NEMO=3, GABA=5, NAL=7, JAKI=10) — run for the model's own full scenario
duration, the entire **24-week horizon (`end=4032h`)**, at `delta=24h`; no
shortening was needed for any of the five (no `maxsteps` issue encountered —
this is a smooth, non-stiff, linear-PK system with bounded logistic-ceiling
disease dynamics). The JAK1i (abrocitinib) scenario additionally set
`JAKI_ON=1` (matching the original's own "Abrocitinib 200mg QD oral
(off-label)" scenario parameters) so its gated effect term is actually
exercised (`EFFECT_JAKI` ranged 0–0.0004 over the run, confirming the gate
and Hill arithmetic are both active, not trivially zero — this small
magnitude reflects the original's own PK/potency calibration for
abrocitinib, `CL_JAKI=42 L/h` against `EC50_JAKI=150 ng/mL`, unchanged by
the refactor).

Every one of the model's 16 disease-side `$CAPTURE`/state outputs
(`BARRIER, TSLP, TH2, TH17, IL31, PSENS, IENFD, CSENS, OPIOID, SCRATCH,
NODULE, WINRS, IGA, SLEEP, DLQI, CUM_SCRATCH`) plus the 3 clinical-flag
outputs (`WINRS_PCT_IMPROVE, IGA_SUCCESS, RESPONDER`) plus each compound's
own PK states under their old/new names were compared per run:

- **Dupilumab**: `DEPOT_DUPI`/`GUT_DUPI`, `C_DUPI`/`CENT_DUPI`,
  `CP_DUPI`/`C_DUPI` — 181 time points.
- **Nemolizumab**: `DEPOT_NEMO`/`GUT_NEMO`, `C_NEMO`/`CENT_NEMO`,
  `CP_NEMO`/`C_NEMO` — 181 time points.
- **Gabapentin**: `GUT_GABA`/`GUT_GABA`, `C_GABA`/`CENT_GABA`,
  `CP_GABA`/`C_GABA` — 673 time points (504 dosing events add extra grid
  points).
- **Nalbuphine ER**: `GUT_NAL`/`GUT_NAL`, `C_NAL`/`CENT_NAL`,
  `P_NAL`/`PERI_NAL`, `CP_NAL`/`C_NAL` — 505 time points.
- **Oral JAK1 inhibitor**: `GUT_JAKI`/`GUT_JAKI`, `C_JAKI`/`CENT_JAKI`,
  `CP_JAKI`/`C_JAKI` — 337 time points.

**Result: exact match, not just near-exact, for all five compounds.** Every
compared output, in every run, matched with **maximum absolute difference
exactly 0.0** — no floating-point-scale residual, no API JSON-rounding
residual, bit-identical as returned. This is consistent with the guide's
expectation for a pure structural reorganization (archetypes 1–3, no
Hill-fitting): same math, reorganized, with no reparameterization needed
anywhere (unlike e.g. `distal-renal-tubular-acidosis`'s HCTZ, none of these
five compounds needed a micro-constant-to-CL/V conversion — the original
already used CL/V1/Q/V2-style parameters throughout, including nalbuphine
once its central/peripheral volume numbering was corrected).

## Anything else flagged

- No upstream defect found or logged in `translations/UPSTREAM_ISSUES.md`.
- The one genuine mechanism-level finding worth a reviewer's attention is
  the `$ODE`+`$CAPTURE`+explicit-`double` result described above (works
  fine here; the opposite of what two other refactors in this fork found
  necessary for their own files) — flagged in case it helps calibrate that
  precedent for future files.
- No compound other than DUPI, GABA, JAKI, NAL, NEMO was touched — this
  file has no other exogenous compound. Every disease-side compartment,
  parameter, and equation is byte-identical to `pn_mrgsolve_model.R`.
- The R-side dosing helpers `ev_dupilumab()` and `ev_nemolizumab()` were
  updated to dose into `"GUT_DUPI"`/`"GUT_NEMO"` (was `"DEPOT_DUPI"`/
  `"DEPOT_NEMO"`) so the refactored sibling's own R scenario-driver code
  stays internally consistent with the renamed compartments; this is a
  rename of the dosing target only (same compartment, same 1-based index,
  same dose amounts/timing), not a behavioral change. `ev_gabapentin()`,
  `ev_nalbuphine_er()`, and `ev_abrocitinib()` needed no change (`GUT_GABA`,
  `GUT_NAL`, `GUT_JAKI` were already named per convention in the original).
- All scratch/debug files used for API verification (`.cpp` DSL
  extractions, request/response JSON, the comparison script) were created
  under the session scratchpad directory only and deleted after use — none
  were left in `prurigo-nodularis/` or anywhere else in the repo.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`prurigo-nodularis | DUPI`, `prurigo-nodularis | GABA`, `prurigo-nodularis |
JAKI`, `prurigo-nodularis | NAL`, and `prurigo-nodularis | NEMO` rows.
