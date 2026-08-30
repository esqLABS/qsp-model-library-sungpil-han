# Refactor notes — `huntingtons-disease/hd_mrgsolve_model.R`

Scope: **TBZ (tetrabenazine), HTBZ ((alpha/beta)-dihydrotetrabenazine),
DTBZ (deutetrabenazine), VBZ (valbenazine), RILUZOLE, and TOMINERSEN
only** — their PK blocks and their downstream effect equations.
Branaplam (modeled only as an on/off flag applying a fixed 50%
multiplier to mHTT mRNA, `dose_branaplam`/`branaplam_eff`, no PK
compartment of its own) is untouched — it is not one of the six
compounds this task covers. All 10 disease-PD compartments and their
kinetics are untouched except for two disclosed build-compatibility
fixes (see below), which touch only the `$INIT`→`$MAIN` idiom, one dead
R-function call, and the addition of `capture` keywords — no formula in
any disease-PD equation was changed.

## HTBZ is TBZ's own active metabolite; DTBZ and VBZ are NOT chained to it

The task brief flagged this as something to check rather than assume,
per the sarcoidosis PRED/PREDL and urolithiasis ALLO/OXP precedent. The
code confirms a genuine parent-metabolite chain for **TBZ→HTBZ only**:

- `HTBZ_brain` has **no `GUT_`/depot compartment and no `ka_`/`F_`
  absorption parameters of its own** anywhere in the file.
- `HTBZ_brain`'s only inflow term is a conversion computed directly from
  TBZ's own saturable (Michaelis-Menten, CYP2D6) metabolic clearance:
  ```
  double R_met_TBZ = (Vmax_cyp * TBZ_plasma) / (km_cyp2d6 + TBZ_plasma); // CYP2D6
  dxdt_TBZ_plasma  = R_ka_TBZ * F_TBZ - R_CL_TBZ - R_met_TBZ;
  double R_HTBZ_in = R_met_TBZ * kp_HTBZ;
  dxdt_HTBZ_brain  = R_HTBZ_in - R_HTBZ_out;
  ```
  `R_met_TBZ` is subtracted from TBZ's own mass balance (real mass
  diversion, not just an informational read) and the same rate (scaled
  by `kp_HTBZ`) forms HTBZ — a genuine, if loosely parameterized,
  metabolic conversion.
- The VMAT2 disease-effect term reads `HTBZ_brain` (via the shared
  `Cp_VMAT2_inhib` sum, see below), **never** `TBZ_plasma` — TBZ itself
  has no direct pharmacodynamic effect in this model; only its
  metabolite HTBZ does.

**Checked the same way for DTBZ and VBZ, and found NOT to be chained to
HTBZ** — despite deutetrabenazine and valbenazine both being, in real
pharmacology, prodrugs of their own active metabolites (deuterium-substituted
dihydrotetrabenazine and NBI-98782, respectively):

- `DTBZ_brain` and `VBZ_brain` are fed from their **own** parent
  compartment (`DTBZ_plasma`, `VBZ_plasma`) via a `kp_*`-scaled rate, not
  from any shared metabolic pathway with TBZ/HTBZ. No state, parameter,
  or expression is shared between `HTBZ_brain` and `DTBZ_brain`/
  `VBZ_brain` anywhere in the file.
- Unlike TBZ's CYP2D6 term, `DTBZ_brain`'s and `VBZ_brain`'s inflow rates
  are **not** subtracted from their own parent's mass balance
  (`dxdt_DTBZ_plasma`/`dxdt_VBZ_plasma` only ever contain a `-R_CL_*`
  term) — i.e. these are one-way, non-mass-conserving "brain exposure"
  compartments of the *same* molecule, not a second chemical species.
- The three (HTBZ, DTBZ, VBZ) are nonetheless functionally
  interchangeable in the original's own VMAT2-inhibition term (see
  below) — all three feed one shared Hill curve — which is exactly why
  this looked, at first glance, like it might be a parent-metabolite
  triad. It is not: it is three independently-dosed drugs which happen
  to converge on the same pharmacological target (VMAT2) and share one
  (re-used, not re-derived) potency curve in the original.

This is flagged as a modeling simplification in the original (real
deutetrabenazine/valbenazine metabolite chemistry is not represented),
not a defect — the original never claims otherwise, so nothing is logged
in `UPSTREAM_ISSUES.md` for this specific point.

## Archetype per compound

**TBZ (tetrabenazine) — Archetype 3 minus the peripheral compartment**
(depot + central, linear elimination) **plus a parallel saturable
(Michaelis-Menten, CYP2D6) metabolic clearance pathway** that forms
HTBZ — not flattened to a plain archetype, since the saturable pathway
is a real part of the original's own pharmacology (same "don't simplify
away a mechanistically richer term" principle the guide applies to
TMDD):
```
dxdt_GUT_TBZ  = -R_ka_TBZ;
dxdt_CENT_TBZ = R_ka_TBZ*F_TBZ - R_CL_TBZ - R_met_TBZ;   // R_met_TBZ = Vmax*C/(Km+C)
```
No direct disease effect (`EFFECT_TBZ` does not exist) — matching the
original having no direct tetrabenazine PD effect (see above).

**HTBZ ((alpha/beta)-dihydrotetrabenazine) — Archetype 1** (no depot,
single compartment, linear elimination), fed by the conversion inflow
above:
```
dxdt_CENT_HTBZ = R_met_TBZ*FCONV_HTBZ - (CL_HTBZ/V1_HTBZ)*CENT_HTBZ;
```
`C_HTBZ = CENT_HTBZ` — the original never divides `HTBZ_brain` by
`Vc_HTBZ` anywhere it is actually read (`Cp_VMAT2_inhib += HTBZ_brain *
1000.0`, `double HTBZ_conc = HTBZ_brain;`); `Vc_HTBZ` (renamed
`V1_HTBZ`) is used only inside the elimination-rate-constant, never as a
concentration divisor. Preserved exactly — this is a "compartment IS
concentration" pattern, same class already documented for `POLY_MTX` in
`sarcoidosis/sarc_refactor_notes.md`.

**DTBZ (deutetrabenazine) — Archetype 1** (no depot, single compartment,
linear elimination) for its own central compartment, **plus a bespoke
one-way, non-mass-conserving brain-exposure compartment** (`PERI_DTBZ`)
that the original does not fit into Archetype 2 (that archetype requires
a bidirectional, mass-conserving `Q`/`V2` exchange; this one only ever
flows plasma→brain, with brain's own elimination never returning mass to
plasma):
```
dxdt_CENT_DTBZ = -R_CL_DTBZ;                              // dosed directly, no depot
dxdt_PERI_DTBZ = CENT_DTBZ*KP_DTBZ*(CL_DTBZ/V1_DTBZ)
                 - (CL_DTBZ/(V1_DTBZ*0.2))*PERI_DTBZ;
```
`KA_DTBZ`/`F_DTBZ` are declared in the original (`ka_DTBZ`/`F_DTBZ`) but
**never referenced by any `dxdt_` line** — the original's own comment
even says so (`// simplified: DTBZ_plasma = gut+plasma`) — preserved
unused, renamed only, same treatment as `ka_predl` in the sarcoidosis
precedent. `C_DTBZ = PERI_DTBZ` (compartment-is-concentration, same
pattern as HTBZ above).

**VBZ (valbenazine) — same structural pattern as DTBZ** (Archetype 1 +
bespoke one-way brain-exposure compartment):
```
dxdt_CENT_VBZ = -R_CL_VBZ;
dxdt_PERI_VBZ = CENT_VBZ*KP_VBZ*0.15 - (CL_VBZ/(V1_VBZ*0.15))*PERI_VBZ;
```
`KA_VBZ`/`F_VBZ` declared but unused, same treatment. `C_VBZ = PERI_VBZ`.

**TOMINERSEN — Archetype 1** (no depot, single compartment, linear
elimination):
```
dxdt_CENT_TOMINERSEN = -(CL_TOMINERSEN/V1_TOMINERSEN)*CENT_TOMINERSEN;
```
The task brief asked to check for "unusual PK" given the intrathecal
route — checked the actual equations and found none: this is
structurally the plainest archetype in the file, despite the CSF/IT
dosing context. `C_TOMINERSEN = CENT_TOMINERSEN` (compartment-is-
concentration, same pattern as above — the original's own `ASO_CSF =
tominersen_CSF` table line confirms this).

**RILUZOLE — Archetype 1** (no depot; `KA_RILUZOLE`/`F_RILUZOLE`
declared but unused, same treatment as DTBZ/VBZ) **plus the same bespoke
one-way brain-exposure compartment pattern**:
```
dxdt_CENT_RILUZOLE = -R_CL_RILUZOLE;
dxdt_PERI_RILUZOLE = CENT_RILUZOLE*KP_RILUZOLE*0.12
                     - (CL_RILUZOLE/(V1_RILUZOLE*0.12))*PERI_RILUZOLE;
```
`C_RILUZOLE = PERI_RILUZOLE`.

## Renaming

| Original | Refactored |
|---|---|
| `ka_TBZ`, `F_TBZ`, `CL_TBZ`, `Vc_TBZ` | `KA_TBZ`, `F_TBZ`, `CL_TBZ`, `V1_TBZ` |
| `km_cyp2d6`, `Vmax_cyp` | `KM_TBZ`, `VMAX_TBZ` |
| `TBZ_gut`, `TBZ_plasma` | `GUT_TBZ`, `CENT_TBZ` |
| `CL_HTBZ`, `Vc_HTBZ` | `CL_HTBZ`, `V1_HTBZ` |
| `kp_HTBZ` | `FCONV_HTBZ` (renamed for its actual role — a conversion scalar on the metabolic rate — not the "brain:plasma Kp" its own comment calls it; see above) |
| `HTBZ_brain` | `CENT_HTBZ` |
| `ka_DTBZ`, `F_DTBZ` | `KA_DTBZ`, `F_DTBZ` (unused in the original — see above) |
| `CL_DTBZ`, `Vc_DTBZ`, `kp_DTBZ` | `CL_DTBZ`, `V1_DTBZ`, `KP_DTBZ` |
| `DTBZ_plasma`, `DTBZ_brain` | `CENT_DTBZ`, `PERI_DTBZ` |
| `ka_VBZ`, `F_VBZ` | `KA_VBZ`, `F_VBZ` (unused) |
| `CL_VBZ`, `Vc_VBZ`, `kp_VBZ` | `CL_VBZ`, `V1_VBZ`, `KP_VBZ` |
| `VBZ_plasma`, `VBZ_brain` | `CENT_VBZ`, `PERI_VBZ` |
| `CL_tominersen`, `Vc_tominersen` | `CL_TOMINERSEN`, `V1_TOMINERSEN` |
| `tominersen_CSF` | `CENT_TOMINERSEN` |
| `ka_riluzole`, `F_riluzole` | `KA_RILUZOLE`, `F_RILUZOLE` (unused) |
| `CL_riluzole`, `Vc_riluzole`, `kp_riluzole` | `CL_RILUZOLE`, `V1_RILUZOLE`, `KP_RILUZOLE` |
| `riluzole_plasma`, `riluzole_brain` | `CENT_RILUZOLE`, `PERI_RILUZOLE` |
| `EC50_VMAT2`, `Emax_VMAT2`, `Hill_VMAT2` (one shared curve) | `EC50_HTBZ`/`EMAX_HTBZ`/`GAMMA_HTBZ`, `EC50_DTBZ`/`EMAX_DTBZ`/`GAMMA_DTBZ`, `EC50_VBZ`/`EMAX_VBZ`/`GAMMA_VBZ` (three duplicated copies, same values, see below) |
| `EC50_ASO`, `Emax_ASO` | `EC50_TOMINERSEN`, `EMAX_TOMINERSEN` |
| — | `GAMMA_TOMINERSEN = 1.0` (new; implicit in the original's `IMAX` plain-ratio shape) |
| `EC50_riluzole` | `EC50_RILUZOLE` |
| `Emax_riluzole = 0.45` (declared, **dead** — never read by the actual formula) | not carried forward under `EMAX_RILUZOLE`; see below |
| — | `EMAX_RILUZOLE = 0.03` (new; the value the original's own formula actually produces — see below), `GAMMA_RILUZOLE = 1.0` (new; implicit) |
| `aso_eff` | `EFFECT_TOMINERSEN` |
| `riluzole_motor` | `EFFECT_RILUZOLE` |
| `VMAT2_inh_effect` | kept (now computed from `EFFECT_HTBZ`/`EFFECT_DTBZ`/`EFFECT_VBZ` via Bliss independence, see below) |

All parameter *values* are copied verbatim from the original — nothing
invented except `GAMMA_TOMINERSEN`/`GAMMA_RILUZOLE` (both `1.0`, implicit
in the original's own plain-ratio Hill shapes) and `EMAX_RILUZOLE`
(discussed next).

## `EMAX_RILUZOLE` — the original's declared `Emax_riluzole` is dead; a different, hardcoded ceiling is what's actually used

Found while isolating riluzole's own block. The original declares
`Emax_riluzole = 0.45 // max glutamate/excitotoxicity reduction` in
`$PARAM`, but the only formula that computes riluzole's motor effect
never reads it:
```
double riluzole_motor = (riluzole_brain > 0.01) ?
    IMAX(riluzole_brain, EC50_riluzole, 0.20) * 0.15 : 0.0;
```
Grepping the whole file for `Emax_riluzole` confirms exactly one
occurrence (its own declaration) — the actual ceiling used is the
hardcoded `0.20` (inner `IMAX` ceiling) times the also-hardcoded `0.15`
outer scaling factor, `0.20*0.15=0.03`. Same defect class as
`UPSTREAM_ISSUES.md` #84 (vte's `EMAX_RIV`/`EMAX_WARF`) — logged as new
entry **#89**. Per that entry's precedent, `EMAX_RILUZOLE` is set to
`0.03` (the value the formula actually produces, correctly wired to
`EFFECT_RILUZOLE`), not to the original's dead `0.45`. `EFFECT_RILUZOLE`
preserves the concentration threshold gate (`C_RILUZOLE > 0.01`) exactly:
```
double EFFECT_RILUZOLE = (C_RILUZOLE > 0.01) ?
    EMAX_RILUZOLE * C_RILUZOLE / (EC50_RILUZOLE + C_RILUZOLE) : 0.0;
```
This is a rename/fold, not a fit — `0.20*0.15` is exact arithmetic, not
a curve fit — confirmed numerically identical in verification below.

## The VMAT2 Hill interface: split into three independent compounds, combined only at the point of use

The original combines **all three** VMAT2 inhibitors into **one shared
Hill curve**, applied to a summed concentration, gated by dose flags:
```
double Cp_VMAT2_inhib = 0.0;
if (dose_TBZ_flag > 0)  Cp_VMAT2_inhib += HTBZ_brain * 1000.0;
if (dose_DTBZ_flag > 0) Cp_VMAT2_inhib += DTBZ_brain * 1000.0;
if (dose_VBZ_flag > 0)  Cp_VMAT2_inhib += VBZ_brain * 1000.0;
double VMAT2_inh_effect = EMAX(Cp_VMAT2_inhib, EC50_VMAT2, Emax_VMAT2, Hill_VMAT2);
```
Per the guide ("keep each compound's `EFFECT_<STEM>` separate; combine
them only at the point the disease equations actually use them — never
collapse several drugs into one shared Hill term"), this is rewritten as
three independent terms (same shared parameter values, three copies,
each individually driveable), combined via Bliss independence only where
`VMAT2_inh_effect` is actually consumed (chorea reduction and dopamine
production):
```
double EFFECT_HTBZ = EMAX(C_HTBZ * 1000.0, EC50_HTBZ, EMAX_HTBZ, GAMMA_HTBZ);
double EFFECT_DTBZ = EMAX(C_DTBZ * 1000.0, EC50_DTBZ, EMAX_DTBZ, GAMMA_DTBZ);
double EFFECT_VBZ  = EMAX(C_VBZ  * 1000.0, EC50_VBZ,  EMAX_VBZ,  GAMMA_VBZ);
...
double VMAT2_inh_effect = 1.0 - (1.0-EFFECT_HTBZ)*(1.0-EFFECT_DTBZ)*(1.0-EFFECT_VBZ);
```
The three `dose_*_flag` gates are dropped from this calculation (no
longer read here) since concentration alone is now zero whenever the
corresponding flag would have been zero in every one of the original's
own scenarios (confirmed by inspection: no scenario doses more than one
of TBZ/DTBZ/VBZ, and none of the three compartments receives any
cross-contamination — unlike tominersen, see next section) — this is
required for pluggability (an external covariate driving `C_HTBZ`
wouldn't also set `dose_TBZ_flag`) and is numerically identical to the
original in every one of the seven named scenarios (confirmed in
verification below). `dose_TBZ_flag`/`dose_DTBZ_flag`/`dose_VBZ_flag`
are kept declared (for the R-side scenario setup, which still sets them)
but are now unread by any equation.

## `EFFECT_TOMINERSEN` deliberately keeps its dose-flag gate — a genuine original defect requires it

Unlike HTBZ/DTBZ/VBZ above, `EFFECT_TOMINERSEN` **preserves** the
original's `(dose_tominersen > 0) ? ... : 0.0` gate rather than reading
concentration alone:
```
double EFFECT_TOMINERSEN = (dose_tominersen > 0) ?
    IMAX(C_TOMINERSEN, EC50_TOMINERSEN, EMAX_TOMINERSEN) : 0.0;
```
Reason: the original's `make_dose_events()` doses **branaplam** (out of
this refactor's scope) into `cmt = 8` — the same compartment tominersen
itself uses (`tominersen_CSF`, renamed `CENT_TOMINERSEN`):
```
} else if (scenario == "Tominersen_Q8W") {
    e_aso <- ev(cmt = 8, amt = 120, ii = 1344, ...)
} else if (scenario == "Branaplam_Q1W") {
    e_bran <- ev(cmt = 8, amt = 50, ii = 168, ...)
```
This is a genuine, pre-existing cross-compound contamination bug —
logged as new `UPSTREAM_ISSUES.md` entry **#90**. During the
`Branaplam_Q1W` scenario, `dose_tominersen = 0` but `CENT_TOMINERSEN`
still receives spurious mass from branaplam's misrouted dose, which then
decays through tominersen's own elimination. If `EFFECT_TOMINERSEN` were
rewritten to read concentration alone (as done for HTBZ/DTBZ/VBZ), this
spurious mass would silently start suppressing mHTT mRNA production
during `Branaplam_Q1W` — a real numeric deviation from the original that
verification would (correctly) catch as a mismatch. Keeping the gate
preserves the original's exact behavior, including this bug, without
"fixing" it in the refactor (per the never-edit-upstream / don't-loosen-
verification-to-force-a-pass rule) — confirmed by the exact match on the
`Branaplam_Q1W` scenario in verification below.

## Build-compat fixes (mrgsolve 2.0.1) — logged as `UPSTREAM_ISSUES.md` #87, #88

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`)
that the **untouched original does not compile at all**, two independent
defects (full detail and error traces in `UPSTREAM_ISSUES.md` #87):

1. `$CMT`+`$INIT` jointly redeclare all 20 compartments (same class as
   issues #29/#34/#36/#59/#76/#81/#82/#83): `Duplicated model names: ...`.
2. `$ODE`'s UHDRS-TMS equation calls `mHTT_protein_reduction_f()` — an
   R-only function defined **outside** the quoted DSL string (a
   placeholder stub, `mHTT_protein_reduction_f <- function() 0.0`,
   several hundred lines below `hd_model <- '...'`) — undefined at
   compile time server-side since `model_content` is compiled with no
   surrounding R script.

**Fixes applied directly to the delivered `hd_mrgsolve_model_refactored.R`**
(not just a scratch copy), all syntax-only and non-numeric:

(a) The `$INIT` block deleted, its 20 assignments moved into `$MAIN` as
`<CMT>_0 = value;` under the refactor's renamed compartments (e.g.
`GUT_TBZ_0 = 0;`, `CENT_HTBZ_0 = 0;`), same values, no compartment added
or removed.

(b) `mHTT_protein_reduction_f()` deleted from the UHDRS-TMS `dxdt_` line
outright — confirmed a no-op, since the R function it referenced always
returned the constant `0.0` regardless of model state:
```
dxdt_UHDRS_TMS = chorea_drive - chorea_tx - EFFECT_RILUZOLE;
// (original's fourth term, - (kprog_TMS*0.3)*mHTT_protein_reduction_f(),
// removed — always evaluated to 0, see UPSTREAM_ISSUES.md #87)
```

A **third**, separate defect was found and fixed the same way, logged as
its own entry (**#88**) because it is a different class (the model
compiles, but produces an incomplete output rather than failing to
build): the original's entire `$TABLE` block is a set of bare `double`
declarations with **no** `capture`/`$CAPTURE` anywhere in the file, so
none of its 17 derived quantities are actually exposed as outputs
(confirmed empirically: `/model_manifest`'s `outputPaths` on the
compile-fixed original lists only the 20 raw compartments). This means
the original's own post-DSL R script — every `summarise()` column and
all six `ggplot()` calls, which reference `TMS`/`TFC`/`MSN_pct`/
`mHTT_total`/`BDNF_level` — would error at runtime on a nonexistent
column; as literally authored, the model's own downstream analysis has
never worked. Fixed by changing `double` to `capture` for all 17
pre-existing lines (no formula changed) plus new `_OUT`-suffixed
captures for this refactor's own `C_<STEM>`/`EFFECT_<STEM>` terms
(suffixed to avoid the `redefinition of capture` collision against the
identically-named `$ODE` locals — same pattern as
`sarcoidosis/sarc_refactor_notes.md` and
`urolithiasis/uri_refactor_notes.md`).

The checked-in original (`hd_mrgsolve_model.R`) is untouched and still
carries all three defects exactly as written; an identically
syntax-fixed scratch copy of the original (same three fixes, original
compartment/param/capture names, never committed) was built in-memory
purely to construct the verification comparison target.

## `$CAPTURE` naming for the six refactored compounds

`C_TBZ`, `C_HTBZ`, `EFFECT_HTBZ`, `C_DTBZ`, `EFFECT_DTBZ`, `C_VBZ`,
`EFFECT_VBZ`, `C_TOMINERSEN`, `EFFECT_TOMINERSEN`, `C_RILUZOLE`,
`EFFECT_RILUZOLE` are all computed as `$ODE`-local `double`s (needed
there for the model's own dynamics), so they are captured under an
`_OUT` suffix (`C_TBZ_OUT`, `C_HTBZ_OUT`, `EFFECT_HTBZ_OUT`, ...) rather
than a bare self-referential `capture C_TBZ = C_TBZ;`, which would be a
`redefinition of capture` under mrgsolve 2.0.1 — same pattern as the
sarcoidosis/urolithiasis precedents. Confirmed discoverable via
`/model_manifest`'s `outputPaths` (47 total, including all 17
pre-existing disease/PK-summary captures plus the 11 new `_OUT` ones).

## qspserver compatibility checklist

- `model_content` is pure mrgsolve DSL text (no R wrapper) — confirmed by
  extracting the quoted `hd_model <- '...'` block byte-for-byte (with
  `\'` unescaped back to `'`) and building it standalone via
  `POST /model_manifest`; the extracted text is byte-identical to the
  block actually embedded in `hd_mrgsolve_model_refactored.R`.
- Every `KA_`/`F_`/`CL_`/`V1_`/`KP_`/`KM_`/`VMAX_`/`FCONV_`/`EMAX_`/
  `EC50_`/`GAMMA_` parameter for all six compounds lives in `$PARAM`
  (confirmed present in the manifest's `parameters` list, 77 total).
- Every `C_<STEM>`/`EFFECT_<STEM>` for all six compounds is a
  discoverable `outputPath` (see above).
- `$SET end/delta` is not present in the original either — the R-side
  `run_scenario()` wrapper supplies `end`/`delta` explicitly; the
  verification below drives `time.end`/`time.delta` through the request.
- No R-only syntax inside the DSL block — the extracted text compiled
  standalone with no surrounding R script (this is exactly what defect
  #87's `mHTT_protein_reduction_f()` violated, and what the fix
  corrects).
- Solver step budget: the full 5-year (43,800h) horizon at `delta=24h`
  (1,826 output points), including TBZ's TID dosing (5,475 doses) and
  DTBZ's BID dosing (3,650 doses) over that window, ran to completion
  with no timeout or step-count error for every scenario tested — no
  shortening was needed.

## Verification

Per the guide's mandatory protocol: extracted the quoted DSL block from
both the build-compat-patched original (scratch copy only, never
committed) and the delivered `hd_mrgsolve_model_refactored.R`, confirmed
both build via `POST /model_manifest` on the qspserver `mrgsolve_api`
container (`http://localhost:8007`, confirmed healthy throughout,
requests spaced ~2s apart), then ran the original file's own seven named
dosing scenarios (from `make_dose_events()`) plus one bespoke scenario
for riluzole (which none of the original's own seven scenarios ever
doses) through `POST /run_simulation`, full 5-year (43,800h) horizon at
`delta=24h` (1,826 points), identical to the original's own
`run_scenario()` default:

| Scenario | Dosing (from the original's own `make_dose_events()`) | Result |
|---|---|---|
| 1. NaturalHistory | no dosing | exact match, max abs diff 0.0 |
| 2. TBZ_25mg | `GUT_TBZ` (cmt 1): 25/3 mg q8h x 5,475 doses | exact match, max abs diff 0.0 |
| 3. DTBZ_30mg | `CENT_DTBZ` (cmt 4): 15 mg q12h x 3,650 doses | exact match, max abs diff 0.0 |
| 4. VBZ_80mg | `CENT_VBZ` (cmt 6): 80 mg q24h x 1,825 doses | exact match, max abs diff 0.0 |
| 5. Tominersen_Q8W | `CENT_TOMINERSEN` (cmt 8): 120 mg q1344h x 33 doses | exact match, max abs diff 0.0 |
| 6. Branaplam_Q1W | `CENT_TOMINERSEN` (cmt 8, misrouted branaplam dose — see UPSTREAM_ISSUES #90): 50 mg q168h x 261 doses | exact match, max abs diff 0.0 (including the spurious `CENT_TOMINERSEN` trajectory) |
| 7. Combo_DTBZ_Tominersen | `CENT_DTBZ` + `CENT_TOMINERSEN` combined | exact match, max abs diff 0.0 |
| 8. Riluzole_bespoke (not an original scenario) | `CENT_RILUZOLE` (cmt 9): 50 mg q12h x 3,650 doses — arbitrary nonzero dose for an already-declared, previously-unexercised compartment | exact match, max abs diff 0.0 |

Every shared output was compared point-by-point across the full
1,826-point time grid: all 20 raw compartments (mapped through the
renames given above) and all 17 pre-existing `$TABLE` captures
(`HTBZ_conc`, `DTBZ_conc`, `VBZ_conc`, `ASO_CSF`, `mHTT_mRNA_rel`,
`mHTT_total`, `oligomer_pct`, `BDNF_level`, `DA_level`, `MSN_pct`,
`OxStress`, `Inflam`, `TMS`, `TFC`, `VMAT2_inh`, `chorea_red_pct`,
`cUHDRS`). **Result: exact match, max abs diff = 0.0, for every output in
all eight scenarios** — the expected result for pure structural
reorganization with no Hill-refitting (every one of the six compounds'
effect terms was already an exact Hill/plain-ratio shape in the
original, see Archetype section above; `EMAX_RILUZOLE`'s value is an
exact arithmetic fold, `0.20*0.15=0.03`, not a fit).

Sanity-checked the new `_OUT` captures are non-trivial: in the TBZ_25mg
scenario, `C_TBZ_OUT` stays near-zero at steady state (~1.8e-7 mg/L,
consistent with TBZ's fast CYP2D6 clearance to its metabolite) while
`C_HTBZ_OUT` rises to 0.444 and `EFFECT_HTBZ_OUT` saturates at 0.85
(=`EMAX_HTBZ`) by year 5, driving `chorea_red_pct` to ~85%.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, all
six `huntingtons-disease` rows (DTBZ, HTBZ, RILUZOLE, TBZ, TOMINERSEN,
VBZ) — HTBZ's `Target/Pathway` corrected from `?` to reflect the
parent-metabolite finding above; the other five rows' `Target/Pathway`
filled in from the actual code (VMAT2 for DTBZ/VBZ, glutamate-release
inhibition for RILUZOLE, none/PK-only-prodrug for TBZ, HTT mRNA/RNase-H1
for TOMINERSEN).

## Anything else flagged

- Branaplam (`dose_branaplam`, `branaplam_eff`, and the `Branaplam_Q1W`
  scenario's dosing mechanics) is out of this refactor's scope — untouched
  except that its dose event still lands in `CENT_TOMINERSEN` (renamed
  from `tominersen_CSF`), exactly as in the original; see UPSTREAM_ISSUES
  #90.
- The post-DSL R script required no changes to any dosing event's `cmt=`
  target (compartment declaration order — and therefore every numeric
  `cmt=` — is unchanged from the original) and no changes to any
  plotting/summary code (every `$TABLE` capture name referenced by
  `summarise()`/`ggplot()` — `TMS`, `TFC`, `MSN_pct`, `mHTT_total`,
  `BDNF_level`, `chorea_red_pct` — was kept unrenamed, only newly
  `capture`d). The only edits to the post-DSL script: `mcode()`'s model
  identifier changed from `"hd_qsp"` to `"hd_qsp_refactored"` (avoids an
  `mrgsolve` `soloc` cache collision if both files are compiled in the
  same R session; non-behavioral), and the now-orphaned R-level
  `mHTT_protein_reduction_f <- function() 0.0` placeholder (previously
  referenced only by the now-deleted dead DSL call) removed entirely
  rather than left dangling.
