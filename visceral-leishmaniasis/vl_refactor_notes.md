# Refactor notes — `visceral-leishmaniasis/vl_mrgsolve_model.R`

Scope: all **five** compound rows this file has in
`driver-patches/data/compound_perturbation_census.md` — **Amphotericin B
disposition (FRE)**, **Amphotericin B disposition (LIP)**, **Antimony
(SB)**, **Miltefosine**, and **Paromomycin (PM)** — all classified "Redirect
concentration (clean single site)". The file models exactly these four
drugs (five PK sub-systems once AmB's two formulations are counted
separately) and nothing else pharmacological, so there is no compound left
untouched to report.

## The duplicate "Miltefosine" census row: a true duplicate, not two entities

Checked the actual code (`$CMT`, `$PARAM`, `$ODE` section 2 "MILTEFOSINE").
There is exactly **one** miltefosine PK sub-system: a depot (`MIL_G`),
central (`MIL_C`), peripheral (`MIL_P`), and four macrophage-organ tissue
compartments (`MIL_SP`/`MIL_LI`/`MIL_BM`/`MIL_SK`) — one absorption route,
one clearance, one `EC50_M`/`EMAX_M`/`HILL_M` Hill term. No second dosing
regimen, no second physiological form, no parent+metabolite pair (contrast
with e.g. `acne-vulgaris`'s ISO/OXO parent+active-metabolite pair, which
*is* two real entities). This matches the census's own caveat calling out
this exact file as a case "worth a quick manual check" — it checked out as
the duplicate-classifier artifact the caveat warned about. The second
`Miltefosine` row has been removed from the census table; one row now
carries the full annotation.

## FRE vs LIP: genuinely two products, not a mislabeled single entity

Amphotericin B is the interesting case, and it is **not** analogous to the
Miltefosine duplicate. The original models two administered products that
diverge in disposition on purpose — this divergence is the file's own
stated "one structural commitment" (see the header comment: liposomal
encapsulation and free/deoxycholate amphotericin B send the same mg/kg dose
down two different physical routes, and that is the entire point of the
model, demonstrated in the file's own `vl_analysis_targeting()`).
Concretely, in the original:

- `A_LIP` (liposome-associated plasma pool) has its **own** volume
  (`V_LIP`) and clearance (`CL_LIP`), and is the dosing site for liposomal
  amphotericin B (`CMT_LAMB = 1`).
- `A_FRE`/`A_PER` (free/protein-bound central + peripheral) has its **own**
  volume/clearance/Q (`V_FRE`/`CL_FRE`/`Q_AMB`/`V_PER`), and is the dosing
  site for amphotericin B deoxycholate (`CMT_DAMB = 2`).
- `A_LIP`'s entire clearance (`lip_out`) splits into two fates: a fraction
  `FREL` **leaks** into `A_FRE` (so LIP feeds FRE), and the remainder is
  **phagocytosed directly** into the four macrophage-organ tissue pools
  (`AMB_SP`/`AMB_LI`/`AMB_BM`/`AMB_SK`) via the `mps` term — bypassing the
  free-plasma pool entirely. `A_FRE`'s concentration separately drives a
  **passive** Kp-scaled equilibration into those same four organs plus the
  renal cortex (`AMB_KID`).

So the two products are real, distinct PK sub-systems (distinct
compartments, distinct volumes/clearances, distinct dosing sites) — but
they are **not fully independent**: both deposit into the *same* downstream
tissue/kidney pools, because a macrophage (or a renal tubule cell) cannot
tell which formulation delivered the amphotericin B molecule sitting in
front of it. This is deliberate, real pharmacology (the entire "liposomal
encapsulation is a targeting device" argument the file's header makes), not
an artifact to collapse or an independence to force. The refactor keeps
this coupling exactly, renamed but mathematically untouched.

## Archetype per compound

**Amphotericin B — LIP: bespoke, no clean archetype match.** A single
compartment (`CENT_LIP`, was `A_LIP`) whose "elimination" splits
immediately into two simultaneous fates (a leak to `CENT_FRE` and a direct
phagocytic deposit into four tissue pools) rather than clearing to a single
next compartment or out of the system — none of the guide's four
archetypes describe this shape, so it is treated as its own thing, renamed
to convention and isolated as its own PK block per the guide's "none of
these fit" carve-out.

**Amphotericin B — FRE: archetype-3-like (central + peripheral, no depot;
IV so no `GUT`/`KA`/`F`), plus a genuinely distinct multi-tissue system.**
`CENT_FRE`+`PERI_FRE` (was `A_FRE`/`A_PER`) is a plain linear two-compartment
PK block. The five downstream tissue sites
(`TISSUE_SP_AMB`/`TISSUE_LI_AMB`/`TISSUE_BM_AMB`/`TISSUE_SK_AMB`/`TISSUE_KID_AMB`,
was `AMB_SP`/`AMB_LI`/`AMB_BM`/`AMB_SK`/`AMB_KID`) are the guide's "two [or
more] tissue sites matter" carve-out taken to its five-organ extreme —
each fed by both `C_FRE` (passive partition) and LIP's `mps` term.

**Antimony (SB): archetype-3-like** — `GUT_SB`+`CENT_SB`+`PERI_SB` (was
`SB5_D`/`SB5_C`/`SB_DEP`). `PERI_SB` is not a classic peripheral
compartment in the "extra volume of distribution" sense; it is a second,
genuinely distinct Sb(V) depot ("deep tissue depot") whose own
concentration (`C_SB_DEEP`, informational, was `C_SBD`) is what drives QTc
prolongation, while the exposed clean single site `C_SB` (plasma Sb(V), was
locally `C_SB5`) drives lipase/ALT. Plus four tissue sites (`TISSUE_SP_SB`
etc., was `SB3_SP` etc.) holding the *reduced* Sb(III) species that
actually kills the parasite (`KRED_SB` reduction from `C_SB`, `keff`
efflux). This structure (linear two-compartment plasma PK, with the actual
pharmacology one metabolic step downstream) most closely resembles
archetype 3 without a strict fit, so it is documented as archetype-3-like
rather than forced into an exact match.

**Miltefosine: archetype 3 (depot + central + peripheral), plus tissue —
clean fit, no coupling.** `GUT_MIL`+`CENT_MIL`+`PERI_MIL` (was
`MIL_G`/`MIL_C`/`MIL_P`) plus four tissue sites (`TISSUE_SP_MIL` etc., was
`MIL_SP` etc.), all fed purely by `C_MIL`. The one compound in this file
with zero cross-compound coupling.

**Paromomycin (PM): archetype 3, plus tissue.** `GUT_PM`+`CENT_PM`+`PERI_PM`
(was `PM_D`/`PM_C`/`PM_P`) plus six tissue sites: four macrophage organs
(`TISSUE_SP_PM` etc., was `PM_SP` etc. — using the original's slow
pinocytic uptake/egress kinetic form, `KIN_PM`\*`C_PM`\*`V` in, `KOUT_PM`\*
tissue out, a genuinely different shape from AmB/miltefosine's
equilibrium-restoring `KOUT*(Kp*C - tissue)` form, preserved unchanged) plus
renal cortex and cochlea (`TISSUE_KID_PM`, `TISSUE_COC_PM`, was `PM_KID`,
`PM_COC`).

## Naming

| Role | LIP | FRE | MIL | PM | SB |
|---|---|---|---|---|---|
| Depot | — (IV) | — (IV) | `GUT_MIL` (was `MIL_G`) | `GUT_PM` (was `PM_D`) | `GUT_SB` (was `SB5_D`) |
| Central | `CENT_LIP` (was `A_LIP`) | `CENT_FRE` (was `A_FRE`) | `CENT_MIL` (was `MIL_C`) | `CENT_PM` (was `PM_C`) | `CENT_SB` (was `SB5_C`) |
| Peripheral | — | `PERI_FRE` (was `A_PER`) | `PERI_MIL` (was `MIL_P`) | `PERI_PM` (was `PM_P`) | `PERI_SB` (deep depot, was `SB_DEP`) |
| Tissue (per organ, SP/LI/BM/SK) | `TISSUE_*_AMB` (shared with FRE) | `TISSUE_*_AMB` (shared with LIP) | `TISSUE_*_MIL` | `TISSUE_*_PM` | `TISSUE_*_SB` |
| Extra tissue | `TISSUE_KID_AMB` (shared) | `TISSUE_KID_AMB` (shared) | — | `TISSUE_KID_PM`, `TISSUE_COC_PM` | — |
| `CL_/V1_/V2_/Q_` | `CL_LIP`, `V1_LIP` (was `V_LIP`) | `CL_FRE`, `V1_FRE` (was `V_FRE`), `V2_FRE` (was `V_PER`), `Q_FRE` (was `Q_AMB`) | unchanged (`CL_MIL`, `V1_MIL`, `V2_MIL`, `Q_MIL`) | `V1_PM` (was `V_PM`); `CL_PM`/`V2_PM`/`Q_PM` unchanged | `V1_SB` (was `V_SB`); `CL_SB`/`V2_SB`/`Q_SB` unchanged |
| `KA_/F_` | — | — | unchanged | unchanged | unchanged |
| Exposed concentration | `C_LIP` (already named this in the original) | `C_FRE` (already named this in the original) | `C_MIL` (already named this in the original) | `C_PM` (already named this in the original) | `C_SB` (was locally `C_SB5`) |
| Hill: Emax/EC50/gamma | `EMAX_AMB`/`EC50_AMB`/`GAMMA_AMB` (was `EMAX_A`/`EC50_A`/`HILL_A`; shared with FRE) | same as LIP | `EMAX_MIL`/`EC50_MIL`/`GAMMA_MIL` (was `EMAX_M`/`EC50_M`/`HILL_M`) | `EMAX_PM`/`EC50_PM`/`GAMMA_PM` (was `EMAX_P`/`EC50_P`/`HILL_P`) | `EMAX_SB`/`EC50_SB`/`GAMMA_SB` (was `EMAX_S`/`EC50_S`/`HILL_S`) |
| Effect | `EFFECT_AMB_SP/LI/BM/SK` (shared with FRE) | same as LIP | `EFFECT_MIL_SP/LI/BM/SK` | `EFFECT_PM_SP/LI/BM/SK` | `EFFECT_SB_SP/LI/BM/SK` |

Other renames: `FREL`→`FREL_LIP`; `FSP_A`/`FLI_A`/`FBM_A`/`FSK_A` (LIP's
per-organ phagocytic split) →`FSP_LIP`/`FLI_LIP`/`FBM_LIP`/`FSK_LIP`;
`KOUT_A`→`KOUT_AMB`; `KPSP_A`/`KPLI_A`/`KPBM_A`/`KPSK_A`/`KPKID_A` (FRE's
passive partition coefficients) →`KPSP_FRE`/etc.; `KIN_KID`/`KOUT_KID`
→`KIN_KID_AMB`/`KOUT_KID_AMB`; `KPSP_M`/etc.→`KPSP_MIL`/etc., `KIN_M`
→`KIN_MIL`; `KIN_PMK`/`KOUT_PMK`→`KIN_KID_PM`/`KOUT_KID_PM`; `KIN_COC`/
`KOUT_COC`→`KIN_COC_PM`/`KOUT_COC_PM`; `KRED`→`KRED_SB`; `RFAC_A`/`RFAC_M`/
`RFAC_P`/`RFAC_S`(EC50 multiplier for the quiescent/persister subpopulation)
→`RFACQ_AMB`/`RFACQ_MIL`/`RFACQ_PM`/`RFACQ_SB`; `RES_M`/`RES_P`(acquired
resistance multiplier)→`RES_MIL`/`RES_PM` (`RES_SB`, the antimony efflux
resistance multiplier, was already stem-named). Every value is copied
verbatim from the original — nothing invented, nothing dropped, nothing
defaulted. Compartment **order is unchanged** (1-based `$CMT` position is
identical to the original throughout), so dosing by `cmt` number
(`CMT_LAMB=1`, `CMT_DAMB=2`, `CMT_MIL=9`, `CMT_PM=16`, `CMT_SB=25` in the
file's own R helpers) is unaffected.

The inert `SB5_P` compartment ("Spare compartment, kept for state-vector
parity", `dxdt_SB5_P = 0.0` always, read nowhere else in the file) is
preserved unchanged, renamed `SPARE_SB` — not removed, per "reorganize,
don't invent, don't drop."

## The Hill interface: four organs count as four tissue sites, twice over

Each of the four drugs kills amastigotes independently in each of four
organs (spleen/liver/marrow/skin), against two subpopulations (susceptible
and a slow-killed quiescent/persister form) — the guide's own "two [tissue
sites] only when a genuinely different tissue site matters" carve-out,
extended to this file's real four-organ, two-subpopulation structure. Per
drug per organ this yields:

- `EFFECT_<STEM>_<ORGAN>` — susceptible-population kill rate, e.g.
  `EFFECT_AMB_SP = EMAX_AMB * pow(cAsp, GAMMA_AMB) / (pow(EC50_AMB,
  GAMMA_AMB) + pow(cAsp, GAMMA_AMB))`, exactly the original's `kSa`
  computation for the spleen block, renamed. **Given `$CAPTURE`** for all
  16 (4 drugs × 4 organs) susceptible-form terms.
- `EFFECT_<STEM>_<ORGAN>Q` — the same drug's kill rate against the
  quiescent/persister subpopulation in that organ (EC50 penalised by
  `RFACQ_<STEM>`, exactly the original's `kRa`/etc.). Kept internal, **not**
  given its own `$CAPTURE` entry — it is a secondary population-specific
  variant of the same named effect, not an independently plumbable site
  (disclosed here rather than silently dropped).

All 32 terms (16 `S` + 16 `Q`) are a pure rename of the original's `kSa`/
`kSm`/`kSp`/`kSs`/`kRa`/`kRm`/`kRp`/`kRs` scratch variables — every one was
already exactly the `_EMX` Hill macro, so this is a rename, not a refit
(`GAMMA_<STEM>` carries the original's own `HILL_<X>` value, e.g.
`GAMMA_AMB = 2.0`, not a default `1`, since the original already had an
explicit Hill coefficient for every one of these four drugs).

**A drafting pitfall worth recording**: the first draft declared each
organ's four `EFFECT_*` variables as one combined `double a, b, c, d;`
predeclaration followed by separate assignment lines (mirroring the
original's own `kSa, kSm, kSp, kSs;` reuse-across-organs pattern). This
compiles fine as ordinary C++, but **mrgsolve's `$CAPTURE` could not see
those variables** — build failed against the qspserver `mrgsolve_api`
container with `'EFFECT_AMB_SP' was not declared in this scope` at the
`$TABLE`/capture step, for every `EFFECT_*` variable and *only* those (the
already-working `C_LIP`/`C_FRE`/`C_MIL`/`C_PM`/`C_SB`, each a single
`double NAME = expr;` statement, compiled and captured without issue).
Since each organ now has its own distinct name (no more cross-organ reuse),
every `EFFECT_*`/`EFFECT_*Q` declaration was rewritten as a single
`double NAME = expr;` statement, which resolved it. This is purely a fix to
this refactor's own first-draft syntax, discovered and corrected before
delivery — not an original-file defect, and not logged in
`UPSTREAM_ISSUES.md`.

## qspserver compatibility

Per the precedent already established in this batch (see e.g.
`abdominal-aortic-aneurysm/aaa_refactor_notes.md`): `$PARAM @annotated`
values compile to read-only `const double&`, so `C_<STEM>`/`EFFECT_<STEM>`
stay plain `double` locals computed in `$ODE` and are exposed for
discovery via `$CAPTURE` (which populates `/model_manifest`'s
`outputPaths`), not as overridable `$PARAM` entries. Confirmed via
`POST /model_manifest` against the qspserver `mrgsolve_api` container
(`http://localhost:8007`) that the refactored DSL compiles and that
`outputPaths` includes `C_LIP`, `C_FRE`, `C_MIL`, `C_PM`, `C_SB`, and all 16
`EFFECT_<STEM>_<ORGAN>` susceptible-form terms, and that `parameters`
includes every renamed `$PARAM` (all Hill parameters, all volumes/
clearances, all partition/rate constants listed in the naming table above).

The bare-DSL `.cpp` extraction (find `vl_code <- '...'`, pull the quoted
contents, unescape R's `\\`→`\`) needed for `model_content` was done as a
mechanical, throwaway verification step, not delivered as a checked-in
file — confirmed byte-identical to the quoted block inside
`vl_mrgsolve_model_refactored.R` before every verification run.

## The original compiles cleanly — no build-compatibility fix needed

Unlike several files in this batch, `vl_mrgsolve_model.R`'s DSL block
**compiles as-written** under mrgsolve 2.0.1 via the qspserver
`mrgsolve_api` container (confirmed via `POST /model_manifest` on the
untouched original, once the R string's `\\`→`\` unescaping was done
correctly — the DSL's only backslash usage, the `_EMX` macro's two
line-continuation backslashes, round-trips exactly). No syntax-only
build-compat fix was needed in `vl_mrgsolve_model_refactored.R`, and no new
entry was added to `translations/UPSTREAM_ISSUES.md` for a compile defect,
because there isn't one.

One pre-existing, non-blocking modelling gap was found and is disclosed but
**not** logged as a numbered `UPSTREAM_ISSUES.md` entry (it does not block
compilation or invalidate any scenario, unlike the build-defect class that
section is for): paromomycin's renal-cortex concentration (`cPkid`, was
`PM_KID`/`VkidA`) is computed in the original but never read by any
`dxdt_*` expression — no renal toxicity readout depends on it, unlike
amphotericin B's analogous `cAkid`, which does drive `TUBI`/`SCR`/`KSER`/
`MGSER`. Preserved as computed-but-unused, exactly as the original had it
(not wired up, not removed).

## Verification

Per the guide's mandatory protocol: ran the original file's own dosing
regimens (from `vl_regimens()`) through both the original and the
refactored DSL via the qspserver `mrgsolve_api` service (`POST
/model_manifest`, `POST /run_simulation`), comparing every shared
`$CAPTURE`d/compartment output point-by-point. Three scenarios were needed
to exercise all five compound rows (no single named scenario doses all
four drugs at once):

**Scenario `S13_lamb5_mil7`** (L-AmB 5 mg/kg × 1 IV, cmt 1, + miltefosine
2.5 mg/kg × 7 PO starting day 1, cmt 9; wt = 50 kg) — exercises **LIP** and
**MIL**. `end=240h, delta=6h` (10 days, covering the full dosing course).
**Result: exact match, max abs diff 0.0**, across all 65 shared outputs
(every PK compartment for all 4 drugs, all 8 amastigote states, all
toxicity/clinical states, all AUC integrals), 49 time points, no NaNs, no
step-count issues.

**Scenario `S16_ssg_pm17`** (SSG/antimony 20 mg/kg × 17 IM, cmt 25, +
paromomycin 15 mg/kg × 17 IM, cmt 16, both starting day 0; wt = 50 kg) —
exercises **SB** and **PM**. `end=432h, delta=6h` (18 days). **Result:
exact match, max abs diff 0.0**, 107 time points, no NaNs, no step-count
issues.

**Scenario `S04_damb_alt30`** (deoxycholate amphotericin B 1 mg/kg IV
q48h × 15, cmt 2; wt = 50 kg) — isolates **FRE** on its own (no other drug
dosed, so any LIP→FRE leak contribution is exactly zero and FRE's own
disposition is tested in isolation). `end=720h, delta=12h` (30 days).
**Result: exact match, max abs diff 0.0**, 76 time points, no NaNs, no
step-count issues.

No scenario needed a shortened window relative to what the original's own
`vl_regimens()` actually doses (the API's default solver budget handled all
three windows without a `maxsteps` error) — the guide's step-count caveat
did not apply here.

**Overall result: exact match everywhere, for every compound.** This is
the outcome the guide anticipates for pure structural reorganization plus
`gamma`-unchanged Hill renames (every Hill term here already carried the
original's own explicit coefficient, e.g. `HILL_A = 2.0`, `HILL_M = 1.5` —
no term needed a `gamma=1` default, and none needed a curve fit; nothing in
this refactor introduces any approximation).

## Anything else worth flagging

- The R-side wrapper code (`vl_simulate`, `vl_regimens`, `vl_scenario_par`,
  `vl_outcome`, `vl_pkdl_init`, all `vl_analysis_*` functions, the
  extinction-floor machinery) is **byte-identical** to the original in
  `vl_mrgsolve_model_refactored.R` — every name it references
  (`P_SP_S`, `PTOT_OUT`, `AUCFRE`, `CLIP_OUT`, `CMIL_OUT`, `TMEM`, `CD4`,
  `SPL`, `HGB`, etc.) was deliberately left unrenamed, because none of
  those are PK compartments/parameters for the five target compounds. Only
  the quoted `vl_code` DSL string changed.
- A subtlety confirmed identical in both models, not a refactor artifact:
  an `$ODE`-local concentration (e.g. `C_LIP`) reflects the **pre-dose**
  state at a dosing time point (computed before the event applies), while
  the corresponding `$TABLE` diagnostic (e.g. `CLIP_OUT`) reflects the
  **post-dose** state (computed after) — visible as `C_LIP=0`/`CLIP_OUT=50`
  at the same `t=0` row in scenario `S13`. This is inherent to mrgsolve's
  ODE→event→TABLE evaluation order and was already true of the original's
  own `C_LIP` local; the refactor changes nothing about it.
- Every organ compartment volume (`VSP`, `VLI`, `VBM`, `VSK`, `VKID`,
  `VCOC`) and every non-PK disease/immunity/clinical/toxicity-rate
  parameter (e.g. `KTUB`, `KWK`, `KQTC`, `KIMM`, `KTP`, `CD4SET`, ...) is
  untouched — these are shared physiology, not any one compound's PK, and
  are outside this refactor's scope.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the five
`visceral-leishmaniasis` rows (`Amphotericin B disposition (FRE)`,
`Amphotericin B disposition (LIP)`, `Antimony (SB)`, `Miltefosine`,
`Paromomycin (PM)`), with the duplicate `Miltefosine` row removed.
