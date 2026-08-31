# Refactor notes — `hypercalcemia-of-malignancy/mah_mrgsolve_model.R`

## Scope: the census's one listed row is not a drug; the file actually carries six

`driver-patches/data/compound_perturbation_census.md` carried exactly one row
for this file: **"Parathyroid (PTH)"**, classified "Delete PK compartment;
concentration is itself the state". Checked against the actual code, per the
task's specific flag (the row name itself already suggested an endogenous
hormone rather than a drug):

- `PTH` is a `$CMT` compartment (`plasma PTH (pmol/L)`) produced entirely by
  the model's own parathyroid gland dynamics — `dxdt_PTH = (S_PTH +
  PTH_ECT) / VD_PTH - KEL_PTH * PTH`, where `S_PTH` is the CaSR-driven
  secretion term and `PTH_ECT` is the *ectopic* (tumour-autonomous) secretion
  rate used by the "Ectopic PTH" aetiology arm. **No `ev()`/dataset dosing
  route targets `PTH` anywhere in the file's 23 scenarios or its R driver
  code** (confirmed by exact-name grep across the whole file). `PTHRP`
  (PTHrP, the actual humoral mediator of HHM) is likewise a pure disease-state
  compartment (`dxdt_PTHRP = PTHRP_BASE + PTHRP_S * TUM - PTHRP_KEL *
  PTHRP`), driven by tumour burden, never dosed.
- This is exactly the same pattern documented in `graves-disease`'s T3/T4/
  TSH/TRAb rows and `wilsons-disease`'s copper pools: a hormone or mediator
  that the disease itself produces and that a drug may act *through*, but
  that is never itself administered. **The census row is corrected below**
  (marked "not a drug — endogenous disease state") rather than force-fitting
  a PK compartment onto a variable with no dosing mechanism at all.

**The census had also undercounted this file, the same failure mode already
seen in `prostate-cancer` (1 listed row against 8 real drugs).** Reading the
model's own `$PARAM`/`$CMT` blocks and its 23 shipped scenarios turned up
**six** genuinely externally-dosed compounds, each with its own PK
compartment(s), elimination parameter(s), and a dosing event in at least one
of the file's own scenarios:

| Compound | Stem | Dosed in scenario(s) | Evidence |
|---|---|---|---|
| Zoledronic acid | `ZOL` | s05–s07, s09, s11–s12, s14–s22 | `ZOL_C`/`ZOL_P`/`ZOL_B`/`ZOL_D`/`ZOL_OC` — 5-compartment self-delivering PK, dosed IV bolus |
| Denosumab | `DMB` | s08, s10 | `DMB_SC`/`DMB_C`/`DMB_RK`, SC depot + TMDD-style RANKL sequestration |
| Calcitonin | `CTN` | s04, s06 | `CTN_SC`/`CTN_C`/`CTN_E` + receptor down-regulation `R_CT`, dosed SC |
| Prednisolone | `PRD` | s13 | `PRD_C`/`PRD_E`, dosed as a direct concentration increment |
| Cinacalcet | `CIN` | s14 | `CIN_C`, dosed as a direct concentration increment |
| Furosemide | `FUR` | s15, s16 | `FUR_C`, dosed as a direct concentration increment |

All six are refactored below. `PTH` (corrected) is the only census row that
needed a classification fix; no other row existed to correct or confirm.

## The original compiles cleanly — no build-compat fix needed

`POST /model_manifest` on the untouched original's own DSL (extracted
verbatim from `mah_code <- '...'`): HTTP 200, all 76 outputs and all 213
`$PARAM @annotated` parameters listed with no error. No entry was added to
`translations/UPSTREAM_ISSUES.md` — there is no defect to log.

## Archetypes — five different shapes for six compounds, none forced

### Zoledronate (ZOL) — bespoke, "the self-delivering drug"

None of the four reference archetypes fit. The original's own five-
compartment structure is central to the model's THESIS ("the uptake rate is
proportional to the very flux the drug abolishes") and is preserved exactly,
renamed only:

| Original | Refactored | Role |
|---|---|---|
| `ZOL_C` (cmt) | `CENT_ZOL` | central (dosed IV bolus) |
| `ZOL_P` (cmt) | `PERI_ZOL` | peripheral |
| `ZOL_B` (cmt) | `BONE_ZOL` | bound to hydroxyapatite bone surface |
| `ZOL_D` (cmt) | `BURIED_ZOL` | buried in matrix (terminal sink) |
| `ZOL_OC` (cmt) | `OC_ZOL` | intra-osteoclast drug — **the PD-driving quantity** |
| `ZOL_KEL` | `CL_ZOL` | renal elimination (GFR-scaled — a genuine, preserved renal covariate coupling) |
| `ZOL_K12`/`ZOL_K21` | `K12_ZOL`/`K21_ZOL` | central↔peripheral micro-rate constants |
| `ZOL_KB` | `KB_ZOL` | hydroxyapatite binding rate |
| `ZOL_BMAX` | `BMAX_ZOL` | bone-surface binding capacity |
| `ZOL_BURY` | `KBURY_ZOL` | burial into matrix |
| `ZOL_KUP` | `KUP_ZOL` | osteoclast uptake rate (flux-dependent — multiplies `J_res/J_RES0`) |
| `ZOL_KOUT_OC` | `KOUT_ZOL` | intracellular elimination |
| `ZOL_EMAX`/`ZOL_EC50` | `EMAX_ZOL_APO`/`EC50_ZOL_APO` | Hill — osteoclast apoptosis |
| `ZOL_IMAX`/`ZOL_IC50` | `EMAX_ZOL_RES`/`EC50_ZOL_RES` | Hill — per-cell resorption inhibition |
| — (none) | `GAMMA_ZOL_APO`, `GAMMA_ZOL_RES` (new) | Hill coefficients [math-implied =1, no explicit exponent in original] |

**Why `K12_ZOL`/`K21_ZOL` stay micro-rate-constants rather than becoming
`Q_ZOL`/`V1_ZOL`/`V2_ZOL`:** the guide's Archetype 2 says to always rewrite
`k12`/`k21` micro-constants to `CL`/`Q`/`V`, but that conversion presumes a
volume of distribution exists to convert with. This compound has none —
`ZOL_C`/`ZOL_P` are raw molar amounts (μmol) with no volume term anywhere in
the original, and inventing `V1_ZOL`/`V2_ZOL` purely to relabel `K12`/`K21`
as `Q_ZOL` would add parameters the original never had, for no structural
gain. Kept as disclosed micro-rate constants instead — a bespoke deviation,
not an oversight.

**Why `C_ZOL` is `OC_ZOL`, not `CENT_ZOL`:** both of zoledronate's effect
terms (`E_zol_apo`, `E_zol_act` in the original) are functions of
intra-osteoclast drug, never of the central/peripheral PK compartments —
`CENT_ZOL` only participates in the *delivery* mechanism (binding bone,
GFR-scaled renal clearance), not in the drug's actual pharmacodynamic
action. Exposing `C_ZOL = OC_ZOL` (the compartment the disease equations
actually read) matches the guide's own TMDD precedent, where the exposed
covariate is "`COMPLEX_TCZ/RTOT_TCZ`... not `C_TCZ` directly" whenever PD
acts through something other than raw central concentration.

**Two named effects, one shared concentration — same pattern as Graves'
PTU:** `EFFECT_ZOL_APO` (osteoclast apoptosis) and `EFFECT_ZOL_RES`
(per-cell resorption inhibition) are two structurally and numerically
distinct actions in the original (`ZOL_EMAX`/`ZOL_EC50` vs
`ZOL_IMAX`/`ZOL_IC50`, different values), both driven by `OC_ZOL`. Kept as
two separate named Hill interfaces rather than force-collapsed into one,
per the guide's explicit allowance for a structure that doesn't fit the
single-`EFFECT_<STEM>` mold.

### Denosumab (DMB) — Archetype 3 without peripheral, bespoke TMDD-adjacent binding

`DMB_C`/`CENT_DMB` is *itself* a concentration (pmol/L) in the original's own
`$CMT` annotation, not an amount — the depot's molar dose is converted by
dividing by `DMB_VD` once during absorption (`DMB_KA*DMB_SC*DMB_F/DMB_VD -
DMB_KEL*DMB_C`), exactly the same "compartment already IS the concentration"
bespoke pattern documented for Graves' MMI/PTU/PROP. `C_DMB = CENT_DMB` is
therefore a plain identity, not a division.

| Original | Refactored | Role |
|---|---|---|
| `DMB_SC` (cmt) | `GUT_DMB` | SC absorption depot |
| `DMB_C` (cmt) | `CENT_DMB` | central, itself a concentration (bespoke identity) |
| `DMB_RK` (cmt) | `COMPLEX_DMB` | denosumab–RANKL complex |
| `DMB_KA` | `KA_DMB` | absorption rate |
| `DMB_F` | `F_DMB` | bioavailability |
| `DMB_VD` | `V1_DMB` | volume (used only in the depot→concentration conversion) |
| `DMB_KEL` | `CL_DMB` | elimination rate constant |
| `DMB_KON`/`DMB_KOFF` | `KON_DMB`/`KOFF_DMB` | RANKL association/dissociation |
| `DMB_KINT` | `KINT_DMB` | complex internalisation |

**Why there is no `EFFECT_DMB`:** unlike the other five compounds, the
original has **no scalar effect term for denosumab at all** in its own
"DRUG EFFECT TERMS" section — denosumab's pharmacology is purely structural
sequestration, built directly into the free-RANKL balance equation
(`dxdt_RKL` includes `-bind_d + unbind_d` where `bind_d =
KON_DMB*RKL*C_DMB`). Manufacturing a scalar `EFFECT_DMB` Hill term would
misrepresent a genuine target-mediated-disposition mechanism as a simple
saturating effect, which the guide's "don't force a bad fit" principle
rules out. `C_DMB` alone satisfies requirement 2 (PD reads exactly one
named concentration — here, the binding reaction's own reactant).
**Requirement 1 (own delimited block) is necessarily incomplete for this
compound**, disclosed rather than hidden: `GUT_DMB`/`CENT_DMB`/`COMPLEX_DMB`
are DMB's own block, but the binding reaction inescapably also modifies
`RKL`, a disease-state compartment shared with the endogenous RANKL/OPG/
osteoclast axis — that cross-coupling **is** denosumab's mechanism of
action (RANKL neutralisation) and preserving it, not decoupling it, is what
"don't flatten a mechanistically rich model" requires. The guide's own
TMDD archetype has the identical property (drug binding modifies a receptor
pool read elsewhere); the only difference here is that the receptor pool
(`RKL`) is not a dedicated, drug-private compartment but the disease's own
free-RANKL state, so it is **not** renamed to `REC_FREE_DMB` — that would
misleadingly imply a private pool where none exists.

### Calcitonin (CTN) — Archetype 3 without peripheral, plus bespoke biophase + receptor tachyphylaxis

Same "compartment already is a concentration" pattern as DMB (`CTN_C`
annotated `IU/L`, `CTN_VD` used only in the depot conversion). No `F_CTN` —
the original has no bioavailability parameter for calcitonin at all
(implicitly F=1), so none is invented.

| Original | Refactored | Role |
|---|---|---|
| `CTN_SC` (cmt) | `GUT_CTN` | SC absorption depot |
| `CTN_C` (cmt) | `CENT_CTN` | central, itself a concentration |
| `CTN_E` (cmt) | `BIO_CTN` | osteoclast biophase (equilibration compartment) — **the PD-driving quantity** |
| `R_CT` (cmt) | `REC_CTN` | calcitonin receptor availability (0–1, tachyphylaxis state) |
| `CTN_KA` | `KA_CTN` | absorption rate |
| `CTN_KEL` | `CL_CTN` | elimination rate constant |
| `CTN_VD` | `V1_CTN` | volume (depot conversion only) |
| `CTN_KE0` | `KE0_CTN` | biophase equilibration rate |
| `CTN_EMAX`/`CTN_EC50` | `EMAX_CTN`/`EC50_CTN` | Hill |
| `CT_KIN`/`CT_KOUT` | `KIN_CTN`/`KOUT_CTN` | receptor resynthesis / down-regulation |
| — (none) | `GAMMA_CTN` (new) | Hill coefficient [math-implied =1] |

`C_CTN = BIO_CTN` (the biophase, not raw plasma) — the same "PD reads a
downstream site, not raw central" choice as ZOL's `OC_ZOL`, and explicitly
sanctioned by the guide's own TMDD precedent.

**Bespoke deviation, disclosed:** `EFFECT_CTN = EMAX_CTN * REC_CTN *
hilln(C_CTN, EC50_CTN, GAMMA_CTN)` is **not** a pure
`Emax·Xᵞ/(EC50ᵞ+Xᵞ)` — it carries an extra multiplicative factor,
`REC_CTN`, a *dynamic* receptor-availability state that itself decays under
sustained occupancy (tachyphylaxis) and only partially recovers. This is
disclosed rather than hidden: the original's own `E_ctn = CTN_EMAX * R_CT *
hillx(CTN_E, CTN_EC50)` already has this exact shape (a time-varying
effective ceiling, not a constant `EMAX`), and the model's own header
thesis is explicitly built on this ("calcitonin receptor availability R_CT
falls 1.00 → 0.26 → 0.11 over [48 h]"). Flattening `REC_CTN` out of
`EFFECT_CTN` to force a textbook Hill shape would delete the exact
mechanism (tachyphylaxis) the model exists to demonstrate — so `EFFECT_CTN`
is exposed as one named variable (satisfying the interface requirement)
whose formula is faithfully, not superficially, a rename of the original's.

### Prednisolone (PRD) — Archetype 1 (bespoke identity) + biophase, three separate Hill interfaces

`PRD_C` is dosed **directly as a concentration increment** ("the amount
given is the concentration increment D*F/Vd", per the original's own
in-code comment) — no depot, no volume conversion anywhere. Same bespoke
identity-mapping as Graves' MMI/PTU/PROP: `dxdt_CENT_PRD = -CL_PRD *
CENT_PRD` reproduces `dxdt_PRD_C = -PRD_KEL * PRD_C` exactly.

| Original | Refactored | Role |
|---|---|---|
| `PRD_C` (cmt) | `CENT_PRD` | central, dosed directly as a concentration bolus (bespoke identity) |
| `PRD_E` (cmt) | `BIO_PRD` | genomic biophase — **the PD-driving quantity** |
| `PRD_KEL` | `CL_PRD` | elimination rate constant |
| `PRD_KE0` | `KE0_PRD` | biophase equilibration rate |
| `PRD_MAC_EMAX`/`PRD_MAC_IC50` | `EMAX_PRD_MAC`/`EC50_PRD_MAC` | Hill — extrarenal 1α-hydroxylase suppression |
| `PRD_GUT`/`PRD_GUT_IC50` | `EMAX_PRD_GUT`/`EC50_PRD_GUT` | Hill — intestinal Ca absorption suppression |
| `PRD_OB`/`PRD_OB_IC50` | `EMAX_PRD_OB`/`EC50_PRD_OB` | Hill — osteoblast apoptosis |
| — (none) | `GAMMA_PRD_MAC`, `GAMMA_PRD_GUT`, `GAMMA_PRD_OB` (new) | Hill coefficients [math-implied =1 each] |

**Three named effects, one shared concentration:** prednisolone genuinely
acts at three pharmacologically distinct sites in the original, each with
its own `EMAX`/`IC50` (extrarenal 1α-hydroxylase, gut Ca absorption,
osteoblast apoptosis) — the same "several real, distinct affinities sharing
one drug concentration" situation as Graves' PTU (TPO inhibition vs
D1-deiodinase). Kept as three separate named interfaces
(`EFFECT_PRD_MAC`/`EFFECT_PRD_GUT`/`EFFECT_PRD_OB`), each an exact rename
(`EMAX`/`GAMMA=1` reproduce the original ratio exactly), rather than
force-collapsed into one term or silently dropping two of the three
mechanisms. `EFFECT_PRD_OB` — previously computed as a local `double
ob_prd_apo` inline inside the bone-cell section rather than in the "DRUG
EFFECT TERMS" section with the other five — was moved up alongside the
other `EFFECT_*` computations for consistency; this changes nothing
numerically (same expression, same inputs, same point of use), only where
in the file it is written.

### Cinacalcet (CIN) and Furosemide (FUR) — plain Archetype 1, bespoke identity, exact renames

Both are dosed directly as concentration increments exactly like
prednisolone (no depot, no volume, single compartment, single Hill effect
term, no biophase):

| Original | Refactored | Role |
|---|---|---|
| `CIN_C` (cmt) | `CENT_CIN` | central, dosed directly (bespoke identity) |
| `CIN_KEL` | `CL_CIN` | elimination rate constant |
| `CIN_EMAX`/`CIN_EC50` | `EMAX_CIN`/`EC50_CIN` | Hill — CaSR set-point left shift |
| — (none) | `GAMMA_CIN` (new) | Hill coefficient [math-implied =1] |
| `FUR_C` (cmt) | `CENT_FUR` | central, dosed directly (bespoke identity) |
| `FUR_KEL` | `CL_FUR` | elimination rate constant |
| `FUR_EMAX`/`FUR_IC50` | `EMAX_FUR`/`EC50_FUR` | Hill — TAL calcium-reabsorption block |
| — (none) | `GAMMA_FUR` (new) | Hill coefficient [math-implied =1] |

`FUR_DIU`/`FUR_NA` (furosemide's own maximal diuresis/natriuresis flux
magnitudes) are left untouched — they are downstream disease-flux
parameters scaled by `EFFECT_FUR/EMAX_FUR` (exactly reproducing the
original's `E_fur/FUR_EMAX` ratio), not part of the compound's own PK/Hill
role table.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Not attempted here as `$PARAM = 0` declarations — the read-only-reference
constraint this produces under mrgsolve 2.0.1 is already documented in
`pagets-disease/pbd_refactor_notes.md` and every refactor in this batch
since. Followed the same, actually-compiling precedent directly: all
fourteen new quantities (`C_ZOL`, `EFFECT_ZOL_APO`, `EFFECT_ZOL_RES`,
`C_DMB`, `C_CTN`, `EFFECT_CTN`, `C_PRD`, `EFFECT_PRD_MAC`,
`EFFECT_PRD_GUT`, `EFFECT_PRD_OB`, `C_CIN`, `EFFECT_CIN`, `C_FUR`,
`EFFECT_FUR`) are plain `double` locals computed in `$ODE`, listed in
`$CAPTURE`. Confirmed via `POST /model_manifest` on the refactored DSL: all
fourteen appear in `outputPaths`, and all forty-eight new/renamed
parameters (`CL_ZOL`, `K12_ZOL`, … `GAMMA_FUR`) appear in `parameters` with
their original numeric defaults — see the full list checked in Verification
below.

## Verification

**Method.** Extracted the bare DSL text from both `mah_mrgsolve_model.R`
(`mah_code <- '...'`, unmodified) and `mah_mrgsolve_model_refactored.R`
(same variable, refactored) and ran them through the local qspserver
`mrgsolve_api` (`http://localhost:8007`), `POST /model_manifest` then `POST
/run_simulation`, requests spaced ~2.5 s apart per the shared-service note.

**A request-contract limitation hit and worked around (not a defect):** the
API's `/run_simulation` contract has no `init`-override field, but every one
of the file's own aetiology definitions (`AETI` list in the R driver) sets a
non-zero tumour-burden initial condition via mrgsolve's `init(TUM = 0.30)`.
Worked around exactly as the `graves-disease` precedent describes: a `t=0`
dosing record adding `amt = 0.30` to the `TUM` compartment (default 0, per
`$MAIN`) reaches the scenario's own `init(TUM = 0.30)` exactly. Disclosed as
a verification-harness technique, not a change to either model.

**A second workaround, specific to this file's "segment" dosing style:**
the original's own R driver (`run_mah`/`merge_segs`) implements the
isotonic-saline schedule (`saline()`) as a piecewise-constant **parameter**
(`IV_RATE`) rather than a compartment dose, re-solving the ODE in segments
that hold `IV_RATE` fixed and restart from the prior segment's end state.
The API's mrgsolve engine supports exactly this through its `drivers` field
(`_runner.R`'s `.mrg_dataset()`, LOCF covariate injection): `drivers:
{"IV_RATE": {"time": [0, 12, 15, 19], "values": [0, 4.80, 2.40, 0], "method":
"constant"}}` reproduces the same step-function schedule as the R driver's
`saline(TX=12)` (hi=4.80 for 3 days, lo=2.40 for 4 days, then 0) in one
continuous run, without the R-side segment-restart machinery — confirmed
functionally equivalent by the exact-match result below, not merely assumed
so.

**Scenarios run — six of the file's own 23, chosen to exercise every one of
the six real compounds at least once (some scenarios dose two compounds
together), all with the file's own doses, aetiology parameters, and saline
schedule:**

1. **s04 — "HHM + saline + calcitonin 4 IU/kg SC q12h × 48 h"** (CTN):
   `ev_bolus(TX + c(0,0.5,1,1.5), GUT_CTN, 280)`.
2. **s05 — "HHM + saline + zoledronate 4 mg IV"** (ZOL):
   `ev_bolus(TX, CENT_ZOL, 14.70)`.
3. **s08 — "HHM + saline + denosumab 120 mg SC"** (DMB):
   `ev_bolus(TX, GUT_DMB, 816000)`.
4. **s13 — "Calcitriol (lymphoma) + saline + prednisolone 60 mg/day × 20 d"**
   (PRD): `ev_bolus(TX + 0:19, CENT_PRD, 1.0667)`.
5. **s14 — "Ectopic PTH + saline + cinacalcet 90 mg bd + zoledronate"**
   (CIN + ZOL together): `ev_bolus(TX + seq(0,29.5,by=0.5), CENT_CIN,
   56.25)` plus `ev_bolus(TX, CENT_ZOL, 14.70)`.
6. **s15 — "HHM + adequate saline + furosemide 40 mg q6h + ZA"**
   (FUR + ZOL together): `ev_bolus(TX + seq(0,2.75,by=0.25), CENT_FUR,
   3.3333)` plus `ev_bolus(TX, CENT_ZOL, 14.70)`.

Each was run for the scenario's own end time (60 days, the default for all
six) at `delta = 0.25` (finer than needed for the day-scale disease
dynamics but coarse enough to keep the shared service's per-request payload
light — no scenario needed the window shortened for a solver-step-count
reason; all six ran to their full 60-day horizon in both models).

**Result: exact match, max abs diff = 0.0, every one of 73 shared outputs,
every scenario** — compared point-by-point across each scenario's full time
grid (243–303 points depending on scenario). The 73 compared quantities per
scenario are (a) 57 disease-side outputs whose names are untouched by this
refactor (`Ca_tot, iCa, Ca_corr, Ca_UF, Pi_c, Mg_c, K_c, pHa, GFR, eGFR_ml,
FL, U_Ca, FE_Ca, Tm_dct, J_res, J_form, J_gut, Q_u, vol, vom, PS, QTc,
KDIGO, PTH, PTHRP, TUM, RKL, OPG, C_LO, OC, OB, OCP, OBP, CTX, P1NP,
CA_BONE, D125, D25, FGF23, AQP2, N_func, CAST, HCO3, A_Ca, V_ecf, SYM,
ADAPT, AUC_CA, MAC, CYT, A_P, A_Mg, A_Na, A_K, TGFB`, plus GFR/renal-flux
terms) and (b) 16 renamed pairs compared old-name-in-original against
new-name-in-refactored (`ZOL_C→CENT_ZOL, ZOL_P→PERI_ZOL, ZOL_B→BONE_ZOL,
ZOL_D→BURIED_ZOL, ZOL_OC→OC_ZOL, DMB_SC→GUT_DMB, DMB_C→CENT_DMB,
DMB_RK→COMPLEX_DMB, CTN_SC→GUT_CTN, CTN_C→CENT_CTN, CTN_E→BIO_CTN,
R_CT→REC_CTN, PRD_C→CENT_PRD, PRD_E→BIO_PRD, CIN_C→CENT_CIN,
FUR_C→CENT_FUR`) — genuinely `0.0` for every point in every scenario, not a
loosened tolerance. This is the expected outcome for six pure structural
reorganizations (all six compounds are renames, not fits — every
`GAMMA_*=1.0` reproduces `hillx()`'s implicit exponent exactly via
`hilln(x,k,1.0)`, which is bit-identical to `hillx(x,k)`).

`E_ctn`/`E_zol_act` were the only two of the original's seven drug-effect
terms it actually exposed via its own `$CAPTURE` (the other five —
`E_zol_apo, E_prd_mac, E_prd_gut, E_fur, E_cin` — were computed in `$ODE`
but never captured, so they are not retrievable from the *original*
compiled model at all, regardless of this refactor). Both of the two that
were capturable matched exactly against `EFFECT_CTN`/`EFFECT_ZOL_RES`
(included in the 16 renamed pairs above). For the five that were never
exposed by the original, correctness was instead confirmed two ways: (a)
indirectly — every downstream disease-state trajectory that depends on
them (`OC`/`MAC`/`U_Ca`/`FE_Ca` etc., all part of the 57 name-identical
outputs) matched the original to `0.0` across the full scenario, which is a
stronger test than an instantaneous local-variable comparison since it
integrates the effect term's numerical behaviour over the whole run; and
(b) directly, by spot-checking each newly-exposed `EFFECT_*` against its own
declared `EMAX`/`EC50`/`GAMMA` at a representative timepoint in the
refactored-only manifest output, e.g.:

- ZOL (s05, t=14): `C_ZOL = 0.4930`, `EFFECT_ZOL_APO = 2.1138 =
  3.4 × 0.4930/(0.30+0.4930)`, `EFFECT_ZOL_RES = 0.4974 = 0.80 ×
  0.4930/(0.30+0.4930)` — both exact.
- PRD (s13, t=14): `C_PRD = 0.0998`, `EFFECT_PRD_MAC = 0.5803 = 0.90 ×
  0.0998/(0.055+0.0998)`, `EFFECT_PRD_GUT = 0.2901 = 0.45 × …`,
  `EFFECT_PRD_OB = 0.2068 = 0.30 × 0.0998/(0.045+0.0998)` — all exact.
- FUR (s15, t=14): `C_FUR = 0.4786`, `EFFECT_FUR = 0.2424 = 0.85 ×
  0.4786/(1.2+0.4786)` — exact.
- CIN (s14, t=14): `C_CIN = 34.68`, `EFFECT_CIN = 0.1627 ≈ 0.28 ×
  34.68/(25+34.68)` — exact.
- CTN (s04, t=14, i.e. TX+2 ≈ 48 h post-therapy-start): `C_CTN = 1.4795`,
  `EFFECT_CTN = 0.0644 = 0.85 × REC_CTN × 1.4795/(0.60+1.4795)` implies
  `REC_CTN ≈ 0.107` — consistent with the model's own header thesis
  claiming `R_CT` falls to "0.11" by 48 h, an independent internal
  consistency check that the receptor-tachyphylaxis mechanism reproduces
  identically post-refactor.

`/model_manifest` on the refactored DSL additionally confirmed all
forty-eight renamed/new parameters (`CL_ZOL, K12_ZOL, K21_ZOL, KB_ZOL,
BMAX_ZOL, KBURY_ZOL, KUP_ZOL, KOUT_ZOL, EMAX_ZOL_APO, EC50_ZOL_APO,
GAMMA_ZOL_APO, EMAX_ZOL_RES, EC50_ZOL_RES, GAMMA_ZOL_RES, KA_DMB, F_DMB,
V1_DMB, CL_DMB, KON_DMB, KOFF_DMB, KINT_DMB, KA_CTN, CL_CTN, V1_CTN,
KE0_CTN, EMAX_CTN, EC50_CTN, GAMMA_CTN, KIN_CTN, KOUT_CTN, CL_PRD, KE0_PRD,
EMAX_PRD_MAC, EC50_PRD_MAC, GAMMA_PRD_MAC, EMAX_PRD_GUT, EC50_PRD_GUT,
GAMMA_PRD_GUT, EMAX_PRD_OB, EC50_PRD_OB, GAMMA_PRD_OB, CL_CIN, EMAX_CIN,
EC50_CIN, GAMMA_CIN, CL_FUR, EMAX_FUR, EC50_FUR, GAMMA_FUR`) appear in
`parameters` with their original numeric values, and all fourteen new/
renamed compartments and covariates (`CENT_ZOL, PERI_ZOL, BONE_ZOL,
BURIED_ZOL, OC_ZOL, GUT_DMB, CENT_DMB, COMPLEX_DMB, GUT_CTN, CENT_CTN,
REC_CTN, CENT_PRD, CENT_CIN, CENT_FUR, BIO_CTN, BIO_PRD, C_ZOL,
EFFECT_ZOL_APO, EFFECT_ZOL_RES, C_DMB, C_CTN, EFFECT_CTN, C_PRD,
EFFECT_PRD_MAC, EFFECT_PRD_GUT, EFFECT_PRD_OB, C_CIN, EFFECT_CIN, C_FUR,
EFFECT_FUR`) appear in `outputPaths`.

No scratch/debug artifacts were left in the repository — DSL extraction,
the six-scenario comparison harness, and all intermediate `.cpp`/JSON files
ran entirely from a session-local scratchpad directory outside the repo
tree, deleted after use.

## R driver consistency

Because six of this file's PK compartments were dosed by the R driver's
`ev_bolus()` calls using **compartment-name strings** (not numeric `cmt`
indices, unlike e.g. `graves-disease`), the same six renames were propagated
into the refactored sibling's own `SCEN` scenario list so the shipped R
driver still runs correctly end-to-end: `"ZOL_C"→"CENT_ZOL"` (16 call
sites), `"DMB_SC"→"GUT_DMB"` (2), `"CTN_SC"→"GUT_CTN"` (2),
`"PRD_C"→"CENT_PRD"` (1), `"CIN_C"→"CENT_CIN"` (1), `"FUR_C"→"CENT_FUR"`
(2) — confirmed by exact-count grep against the original before and after.
No other string in the R driver (aetiology parameter names, dose-constant
variable names, comments) needed to change, since none of those are
compartment references.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`: the
`Parathyroid (PTH)` row reclassified "not a drug — endogenous disease
state" (was previously left as a placeholder classification with no
resolution); six new rows added for Zoledronate, Denosumab, Calcitonin,
Prednisolone, Cinacalcet, and Furosemide, each marked refactored/verified
per this file.
