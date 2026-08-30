# Refactor notes — `mps1_mrgsolve_model.R`

**Scope of this pass.** Per the fork's PK/PD refactor spec
([`FORK_WORKFLOW_GUIDE.md`](../FORK_WORKFLOW_GUIDE.md), Part 2), only the two
compounds flagged in
[`driver-patches/data/compound_perturbation_census.md`](../driver-patches/data/compound_perturbation_census.md)
as "Redirect concentration (clean single site)" were rewritten: **Genistein**
(`GEN`) and **Laronidase, plasma** (`LARO`). Every other intervention in the
file — anti-laronidase antibody (ADA) formation, HSCT engraftment/chimerism,
CNS microglial replacement, investigational AAV9 CNS-directed gene therapy,
investigational lentiviral HSC gene therapy (LVGT) — is not a compound with
its own PK and was left completely untouched, confirmed by diff (see below).

## Archetype determination

### Laronidase — bespoke (two-compartment, one-way, no back-flux)

The original models laronidase PK as two compartments with **irreversible**
one-way transfer, not the guide's Archetype 2 (which has bidirectional
`Q/V1`/`Q/V2` exchange): `CENT_LARO` (plasma, mg) empties entirely into
`TISSUE_LARO` (tissue-retained active-enzyme pool, mg) via
`CL_LARO/V1_LARO`, and `TISSUE_LARO` alone decays via its own first-order
rate (`K_TISSUE_DECAY_LARO`). This reflects the real pharmacology described
in the file's own header comment: laronidase's fast plasma clearance
(t1/2 ~1.5–3.6 h) is receptor (M6PR)-mediated cellular uptake, not
elimination — the cleared drug is *routed into*, not lost from, the system,
and the tissue pool's much slower decay (functional t1/2 ~5 d) is what
sustains biochemical effect between weekly infusions. None of the guide's 4
archetypes has this "depot empties one-way into a second compartment with
its own independent decay, no exchange back" shape, so per the guide's "None
of these fit" clause this is handled as bespoke: renamed to the convention,
isolated from PD, single exposed concentration — but not forced into
Archetype 2's bidirectional-exchange shape.

**Two exposed concentrations, both genuinely necessary** (the guide's stated
exception: "two only when a genuinely different tissue site matters"):
- `C_LARO_PLASMA` (`= CENT_LARO/V1_LARO*1000`, ng/mL) — feeds the ADA
  immunogenicity sub-model's exposure signal only.
- `C_LARO` (`= TISSUE_LARO`, mg) — the actual driver of `EFFECT_LARO`
  (the therapeutic enzyme-replacement effect on systemic GAG access). This is
  the quantity the disease equations read, not `C_LARO_PLASMA`.

Both are new names for quantities the original already computed (as `LARO_CP`
and the bare `LARO_TISSUE` state respectively) — no new computation was
added, only naming and, for `C_LARO`, an explicit one-line alias.

### Genistein — Archetype 3 minus peripheral compartment

Oral depot (`GUT_GEN`) → central (`CENT_GEN`), no peripheral compartment,
first-order absorption and elimination. As in
`aneurysmal-subarachnoid-hemorrhage/sah_mrgsolve_model_refactored.R`'s
nimodipine ("central compartment stores concentration directly"), `CENT_GEN`
was already a **concentration** compartment in the original (`mg/L`, not an
amount divided by volume downstream) — preserved unchanged. The original's
micro-constant elimination (`KE_GEN`, a bare rate constant) was rewritten to
the corpus's preferred `CL/V` convention exactly as
`sah_refactor_notes.md` did for `CL_NIM`: **`CL_GEN := KE_GEN x V1_GEN =
0.12 x 25 = 3.0` (L/h)** is an exact arithmetic recombination of the
original's own two numbers, not an invented value — `CL_GEN/V1_GEN` reduces
to precisely `0.12` again, bit-identical to the original `KE_GEN` term
(confirmed by the verification below, not just asserted).

## Renaming applied

| Original | Refactored | Value | Role |
|---|---|---|---|
| `LARO_CENT` (cmt) | `CENT_LARO` | — | central (plasma) |
| `LARO_TISSUE` (cmt) | `TISSUE_LARO` | — | tissue-retained active-enzyme pool (bespoke role, not in the naming table) |
| `K_TISSUE_DECAY` | `K_TISSUE_DECAY_LARO` | 0.0058 (1/h) | tissue-pool decay rate |
| `LARO_CP` | `C_LARO_PLASMA` | — | plasma conc., ADA-exposure signal only |
| — (was the bare `LARO_TISSUE` state) | `C_LARO` (new alias) | — | **the exposed concentration** PD reads |
| `EMAX_ERT_SYS` | `EMAX_LARO` | 0.70 | Hill Emax |
| `EC50_ERT_TISSUE` | `EC50_LARO` | 1.5 (mg) | Hill EC50 |
| — (none; no explicit Hill term) | `GAMMA_LARO` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| `ERT_EFFECT_SYS` | `EFFECT_LARO` | — | compound's effect on disease (Hill x (1-ADA_INHIB), same formula) |
| `F_LARO_CENT` | `F_CENT_LARO` | — | bioavailability var, renamed to track `\$CMT` |
| `GEN_GUT` (cmt) | `GUT_GEN` | — | oral depot |
| `GEN_CENT` (cmt) | `CENT_GEN` | — | central (stores concentration directly) |
| `KE_GEN` + `V_GEN` | `CL_GEN` (new, `= KE_GEN x V1_GEN`) + `V1_GEN` | `CL_GEN = 3.0` (L/h), `V1_GEN = 25` (L) | clearance / volume |
| `F_GEN_GUT` | `F_GUT_GEN` | — | bioavailability var, renamed to track `\$CMT` |
| — (was the bare `CENT_GEN` state) | `C_GEN` (new alias) | — | **the exposed concentration** |
| `EMAX_GEN_SYNTH_RED` | `EMAX_GEN` | 0.25 | Hill Emax |
| `EC50_GEN` | `EC50_GEN` | 0.15 (mg/L) | Hill EC50 (name already conformed) |
| — (none) | `GAMMA_GEN` (new) | 1 | Hill exponent, "original had no explicit Hill term" |
| `GEN_SYNTH_RED` | `EFFECT_GEN` | — | compound's effect on disease |
| `KA_GEN`, `F_GEN` | unchanged | 0.55, 0.25 | already conformed to convention |

`ADA`, `K_ADA`, `ADA_TARGET`, `ADA_INHIB_MAX`, `EC50_ADA`, `ADA_INHIB`, and
all HSCT/AAV9/LVGT parameters, compartments, and equations are byte-identical
to the original — none is a compound's own PK, so none is in scope, and the
only edit touching the ADA block at all is updating its one reference to the
renamed `C_LARO_PLASMA` (was `LARO_CP`).

**Both Hill terms are rename-only, not refits.** Laronidase's
`EMAX_ERT_SYS * TISSUE_LARO/(EC50_ERT_TISSUE + TISSUE_LARO)` and genistein's
`EMAX_GEN_SYNTH_RED * GEN_CENT/(EC50_GEN + GEN_CENT)` were both already a
plain `Emax x X/(EC50+X)` ratio (implicit Hill coefficient of 1) — exactly
the guide's "already this shape: this is a rename, not a refit" case. Both
were rewritten to the explicit `pow(X, GAMMA)/(pow(EC50, GAMMA)+pow(X,
GAMMA))` form with `GAMMA_LARO = 1` / `GAMMA_GEN = 1`, which is
algebraically identical to the original ratio for any X — confirmed
bit-identical below, not just asserted from the algebra.

## Diff scope confirmation

`diff mps1_mrgsolve_model.R mps1_mrgsolve_model_refactored.R` touches
exactly: the Laronidase PK/PD `\$PARAM` groups, the Genistein PK/PD `\$PARAM`
groups, the four renamed `\$CMT` lines, the two `\$MAIN` bioavailability
lines, the Laronidase and Genistein `\$ODE` blocks (plus the one
`C_LARO_PLASMA` reference inside the ADA block and the three renamed
`EFFECT_LARO`/`EFFECT_GEN` usage sites downstream), the `\$CAPTURE` line, and
three R-side event/scenario-definition strings (`"LARO_CENT"` ->
`"CENT_LARO"` x2, `"GEN_GUT"` -> `"GUT_GEN"` x1) so the refactored file's own
scenario list still doses the right compartments by name. Nothing else in
the file differs — the header comment/citation block, the ADA/HSCT/AAV9/LVGT
`\$PARAM`/`\$CMT`/`\$ODE` content, all 10 scenario definitions' parameter
overrides, `run_scenario()`, and the calibration notes tail are
byte-identical.

## Verification

**Method.** Both files' embedded `code <- '...'` DSL blocks were mechanically
extracted (regex on the `mps1_code <- '...'` assignment, verbatim quoted
text) and POSTed to the local qspserver `mrgsolve_api` service at
`http://localhost:8007` (`POST /model_manifest` then `POST
/run_simulation`), which compiles and runs each block directly with mrgsolve
2.0.1 server-side. The extracted refactored `.cpp` text was confirmed
byte-identical to the `_refactored.R`'s own inline copy before use, then
discarded (no `.cpp` sibling left behind, per the guide). **No pre-existing
upstream compile defect was found** — both the original and the refactored
DSL compiled cleanly on the first attempt.

Four of the file's own 10 named scenarios were run through both models,
identical dosing/parameters, full 5-year horizon (`end=43800 h, delta=24 h`,
1828 timepoints) — **no shortened window was needed**; this model has no
step-count issue even under scenario 8's dense weekly-infusion +
daily-oral-dosing schedule (verified to complete against the API's default
solver budget):

1. **Scenario 1 — "Untreated natural history, severe Hurler"** (no dosing,
   `PHENOTYPE=1, HSCT_FLAG=0`): exercises the pure GAG-turnover/mortality
   equations with both drug compartments at zero.
2. **Scenario 2 — "Laronidase ERT, attenuated Hurler-Scheie"**
   (`laro_ev`: 7.54 mg IV weekly x260, 4 h infusion, `PHENOTYPE=0`):
   laronidase alone.
3. **Scenario 6 — "ERT, high ADA immunogenicity"** (as scenario 2 plus
   `ADA_TARGET=0.85`): exercises the `C_LARO_PLASMA` -> ADA -> `ADA_INHIB`
   -> `EFFECT_LARO` chain at its most active.
4. **Scenario 8 — "Genistein SRT adjunct to ERT"** (`laro_ev` + 130 mg
   genistein PO daily x1825, `PHENOTYPE=0`): both compounds active
   simultaneously.

**Result: exact match in all four scenarios. Maximum absolute and relative
deviation observed = 0.0 (bit-identical) across every shared `\$CAPTURE`d
output (`UGAG`, `LIVSPLEEN`, `VALVE`, `FVC`, `AHI`, `JOINTROM`, `CORNEA`,
`DQ`, `HEIGHTZ`, `HAZARD`, `SURVIVAL`, `ADA_INHIB`, `ENZ_ACCESS_SYS`, all
compartment states) over the full 1828-point time grid.** This is the
expected outcome for a pure structural reorganization with rename-only Hill
terms (no fitting performed, none needed) per the guide's tolerance table.

## qspserver `/model_manifest` discoverability

- `CL_LARO`, `V1_LARO`, `K_TISSUE_DECAY_LARO`, `EMAX_LARO`, `EC50_LARO`,
  `GAMMA_LARO`, `KA_GEN`, `CL_GEN`, `V1_GEN`, `F_GEN`, `EMAX_GEN`, `EC50_GEN`,
  `GAMMA_GEN` are all declared in `\$PARAM` and confirmed present with their
  defaults in the live `/model_manifest` response.
- `C_LARO`, `C_LARO_PLASMA`, `EFFECT_LARO`, `C_GEN`, `EFFECT_GEN` are all
  state-derived (recomputed every `\$ODE` step from compartment values, not
  fixed covariates) and are exposed via `\$CAPTURE` only, **not** `\$PARAM`.
  This was verified directly, not just assumed from precedent: a scratch
  variant adding e.g. `C_LARO : 0` to `\$PARAM` alongside the existing
  `double C_LARO = TISSUE_LARO;` in `\$ODE` failed to compile under
  mrgsolve 2.0.1 (`error: assignment of read-only reference 'C_LARO'`, and
  likewise for the other four) — confirming the same
  `\$PARAM`-vs-`\$ODE`-local mutual exclusion already documented in
  `amd_refactor_notes.md`, `sah_refactor_notes.md`, and `ted_refactor_notes.md`
  for this engine/plugin combination (here surfacing as a
  read-only-reference error under `\$PLUGIN autodec`, rather than the
  "redefinition" form those other files hit — same underlying constraint).
  This is a general mrgsolve/autodec behavior, not a defect in this file, so
  no `UPSTREAM_ISSUES.md` entry was added for it. All five remain fully
  discoverable via `/run_simulation`'s `outputs` selection through
  `\$CAPTURE`, just not via `/model_manifest`'s parameter list.
- `model_content` is pure DSL text extracted from the `code <- '...'`
  R-string wrapper (this file uses that pattern) — no `.cpp` sibling was
  left behind; the extraction was in-memory only, used to build the
  verification requests above and then discarded.

## Anything else flagged

- No pre-existing upstream compile defect was found in this file requiring
  an `UPSTREAM_ISSUES.md` entry.
- `CART_PENETRANCE` (fractional penetration of the *composite* systemic
  enzyme-access index into cartilage/cornea) was left unchanged and
  unrenamed — it modifies `ENZ_ACCESS_SYS` (which combines laronidase, HSCT,
  and LVGT contributions), not laronidase specifically, so it is out of
  scope exactly as `CLF`/`WT`-style shared covariates were left alone in
  prior calibration runs.
- `ADA_TARGET`, `ADA_INHIB_MAX`, `EC50_ADA`, `K_ADA` (the anti-laronidase
  immunogenicity sub-model) were left unrenamed. They are laronidase-specific
  in meaning but are not PK parameters in the naming convention's role table
  (no "immunogenicity" role listed) — only their one point of contact with
  the renamed PK (`C_LARO_PLASMA`, formerly `LARO_CP`) was updated.
