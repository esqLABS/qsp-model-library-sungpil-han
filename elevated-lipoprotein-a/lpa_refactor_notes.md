# Refactor notes — `elevated-lipoprotein-a/lpa_mrgsolve_model.R`

Six compounds refactored, per the existing rows in
`driver-patches/data/compound_perturbation_census.md` (all six classified
"Redirect concentration (clean single site)"): **EVO**, **MUV**, **NIA**,
**OBI**, **STA**, **ZIL**. The file also models pelacarsen (`PEL_*`, an
ASO, cmt 1-3) and a shared olpasiran/lepodisiran/zerlasiran siRNA block
(`SIR_*`, cmt 4-6) — neither has a census row in this file, both are
completely untouched.

## What each abbreviation actually turned out to be

The task brief's a-priori expectation (OBI = olpasiran/an siRNA, ZIL =
zilebesiran or pelacarsen) does **not** match this file's own code —
checked directly, not assumed, per the brief's own instruction to verify:

| Census abbr | Actual compound in this file | PK block |
|---|---|---|
| `EVO` | **Evolocumab** (anti-PCSK9 mAb) | `EVO_SC`/`EVO_CE`, cmt 9-10 |
| `MUV` | **Muvalaplin** (oral Lp(a)-assembly inhibitor) | `MUV_GU`/`MUV_CE`, cmt 7-8 |
| `NIA` | **Niacin ER** | `NIA_CE`, cmt 14 |
| `OBI` | **Obicetrapib** (oral CETP inhibitor — small molecule, NOT olpasiran) | `OBI_CE`, cmt 17 |
| `STA` | **Statin** (rosuvastatin-equivalent) | `STA_GU`/`STA_CE`, cmt 12-13 |
| `ZIL` | **Ziltivekimab** (anti-IL-6 ligand mAb — NOT zilebesiran/pelacarsen) | `ZIL_SC`/`ZIL_CE`, cmt 15-16 |

The actual siRNA (olpasiran/lepodisiran/zerlasiran, `SIR_*`) and ASO
(pelacarsen, `PEL_*`) compounds use different prefixes entirely and are
not among the six census rows for this file — out of scope, untouched.
Consequently **none of the six in-scope compounds is TMDD or an
siRNA/ASO**: every one of them is plain linear PK with no receptor-
binding ODEs, checked against the actual equations rather than assumed.

## Archetype per compound

All six are variants of the guide's Archetype 1/3 family — no peripheral
compartment anywhere, no TMDD:

- **EVO (evolocumab)** — Archetype 3 variant (depot + central, **no
  peripheral**): `GUT_EVO`→`CENT_EVO`, linear absorption/elimination.
- **MUV (muvalaplin)** — Archetype 3 variant (depot + central, no
  peripheral): `GUT_MUV`→`CENT_MUV`.
- **NIA (niacin)** — Archetype 1 (single compartment, no depot): the
  original doses `NIA_CE` directly (no gut compartment at all, i.e. an
  already-absorbed/IV-equivalent input), preserved as-is.
- **OBI (obicetrapib)** — Archetype 1 (single compartment, no depot):
  same as NIA, `OBI_CE` dosed directly, no gut compartment in the
  original.
- **STA (statin)** — Archetype 3 variant (depot + central, no
  peripheral): `STA_GU`→`STA_CE`.
- **ZIL (ziltivekimab)** — Archetype 3 variant (depot + central, no
  peripheral): `ZIL_SC`→`ZIL_CE`.

No compound needed a peripheral compartment, since none of the six had
one in the original — nothing was added or removed, only renamed.

### Why `KE_<STEM>`/`V1_<STEM>` rather than forcing `CL_<STEM>`

The original parametrizes every one of these six compounds' elimination
as a single micro-rate constant (`KEEVO`, `KEMUV`, `KESTA`, `KENIA`,
`KEZIL`, `KEOBI`) times a volume of distribution — never an explicit
`CL`/`V` pair. Converting this to `CL_<STEM> = KE_<STEM>_orig *
V_<STEM>_orig` as a new derived `$PARAM` literal was considered and
rejected: multiplying then dividing back by the same volume at run time
(`CL/V1`) is not guaranteed bit-identical to the original's direct
`KE*CENT` term (floating-point associativity), which would introduce
avoidable drift for a six-compound refactor whose whole point is a clean
verification. Instead, the elimination rate constant is kept and renamed
`KE_<STEM>` (exactly the same choice already made in this fork for
dengue's antiviral, `KEL_AV`, per `dengue/denv_refactor_notes.md`) —
still exactly one clearance-role parameter per compound, still fully
named to the guide's convention, but computing the *identical*
arithmetic as the original, not a re-derived approximation of it. `V1_<STEM>`
(exposed for `C_<STEM> = CENT_<STEM>/V1_<STEM>`) is the original's own
volume of distribution parameter, renamed only.

## Renaming map (values unchanged from the original)

| Original | Refactored |
|---|---|
| `MUV_GU`/`MUV_CE` | `GUT_MUV`/`CENT_MUV` |
| `EVO_SC`/`EVO_CE` | `GUT_EVO`/`CENT_EVO` |
| `STA_GU`/`STA_CE` | `GUT_STA`/`CENT_STA` |
| `NIA_CE` | `CENT_NIA` |
| `ZIL_SC`/`ZIL_CE` | `GUT_ZIL`/`CENT_ZIL` |
| `OBI_CE` | `CENT_OBI` |
| `KAMUV`/`KEMUV`/`VMUV`/`KIMUV` | `KA_MUV`/`KE_MUV`/`V1_MUV`/`EC50_MUV` |
| `KAEVO`/`KEEVO`/`VEVO`/`KONE` | `KA_EVO`/`KE_EVO`/`V1_EVO`/`K_EVO_PCSK9` |
| `KASTA`/`KESTA`/`VSTA`/`EC50STA` | `KA_STA`/`KE_STA`/`V1_STA`/`EC50_STA` |
| `ESRE`/`EPSK`/`ESTA` | `EMAX_STA_LDLR`/`EMAX_STA_PCSK9`/`EMAX_STA_LPA` |
| `KENIA`/`VNIA`/`EC50NIA` | `KE_NIA`/`V1_NIA`/`EC50_NIA` |
| `ENIA`/`ENIALDL` | `EMAX_NIA_LPA`/`EMAX_NIA_LDL` |
| `KAZIL`/`KEZIL`/`VZIL`/`IC50Z` | `KA_ZIL`/`KE_ZIL`/`V1_ZIL`/`EC50_ZIL` |
| `KEOBI`/`VOBI`/`EC50OBI` | `KE_OBI`/`V1_OBI`/`EC50_OBI` |
| `ECETPR`/`EOBI`/`EHDL` | `EMAX_OBI_LDLR`/`EMAX_OBI_LPA`/`EMAX_OBI_HDL` |
| `CMUV`/`CEVO`/`CSTA`/`CNIA`(uncaptured)/`CZIL`/`COBI` | `C_MUV`/`C_EVO`/`C_STA`/`C_NIA`/`C_ZIL`/`C_OBI` |
| `FSTAT`/`FNIAC`/`FOBIC` | `EFFECT_STA_FRAC`/`EFFECT_NIA_FRAC`/`EFFECT_OBI_FRAC` |

New (genuine new `$PARAM`, not derived from any existing value —
Hill-shape constants the original left implicit): `GAMMA_MUV = 1`,
`GAMMA_STA = 1`, `GAMMA_NIA = 1`, `GAMMA_ZIL = 1`, `GAMMA_OBI = 1`,
`EMAX_MUV = 1`. All parameter *values* otherwise are copied verbatim —
nothing invented or defaulted. **Note:** `C_NIA` was not captured at all
in the original (`$CAPTURE` lists `CPELL CMUV CEVO CSTA CZIL COBI ...` —
`CNIA` is conspicuously absent, an oversight in the original, confirmed
by grep). It is now captured in the refactored file (required for niacin
to be discoverable as a pluggable covariate) — a pure addition, not a
value change; niacin's own concentration was always computed identically
either way.

## Hill interface: renames plus exact algebraic derivations, no fitting

Every one of the six compounds' disease-facing effect terms was already
a plain Michaelis–Menten/Hill-shaped ratio in the original (implicit
`gamma = 1` throughout). No `nls()` fitting was needed or performed for
any of them.

- **Muvalaplin**: the original's competitive-blockade fraction
  `FMUV = 1/(1+C/Ki)` is algebraically `Ki/(Ki+C) = 1 - C/(Ki+C)` — an
  exact Hill ratio in disguise. `EFFECT_MUV = EMAX_MUV*C_MUV^GAMMA_MUV /
  (EC50_MUV^GAMMA_MUV + C_MUV^GAMMA_MUV)` with `EMAX_MUV=1`,
  `GAMMA_MUV=1`; `FMUV = 1 - EFFECT_MUV` used exactly where the original
  used `FMUV`.
- **Ziltivekimab**: same shape. Original `IL6EFF = IL6/(1+CZIL/IC50Z) =
  IL6*IC50Z/(IC50Z+CZIL)`. `EFFECT_ZIL = C_ZIL/(EC50_ZIL+C_ZIL)` (gamma=1);
  `IL6EFF = IL6*(1 - EFFECT_ZIL)`, algebraically identical.
- **Statin, niacin, obicetrapib**: each already computed a shared
  saturating fraction (`FSTAT`, `FNIAC`, `FOBIC` respectively) multiplied
  by a **different Emax weight at each of several downstream points of
  use** (statin: PCSK9 synthesis, LDLR synthesis, LPA transcription;
  niacin: LPA transcription, LDL production; obicetrapib: LDLR
  synthesis, LPA transcription, HDL-C). Renamed to one shared
  `EFFECT_<STEM>_FRAC` (the Hill ratio itself) feeding several
  independently-named `EFFECT_<STEM>_<TARGET>` terms — exactly the same
  pattern already used for statin in the `abdominal-aortic-aneurysm`
  refactor (`EFFECT_STAT_MMP9`/`NFKB`/`ROS`). Pure rename, `pow(x,1)==x`
  exactly per IEEE-754, so no drift from introducing `GAMMA_*=1`
  explicitly.
- **Evolocumab — bespoke, deliberately NOT forced into the Hill
  interface.** The original expresses evolocumab's entire disease-facing
  effect as a first-order mass-action clearance term added directly into
  the shared PCSK9 turnover ODE: `dxdt_PCSK9 = ... - KONE*CEVO*PCSK9`.
  This is not a saturating ratio anywhere in the original — it is a rate
  contribution proportional to concentration with no asymptote (PCSK9
  clearance keeps rising linearly with `C_EVO`, it never saturates).
  Per the guide's own rule ("if a compound's PK genuinely doesn't
  resemble any archetype, don't force it into one... a clean, non-standard
  structure beats a standard structure that's wrong"), the same principle
  is applied here to the *effect* term: `EFFECT_EVO = K_EVO_PCSK9*C_EVO`
  (renamed `KONE`→`K_EVO_PCSK9`), used at its point of use in the shared
  PCSK9 ODE (`- EFFECT_EVO*PCSK9`), with no `EMAX_EVO`/`EC50_EVO`
  introduced since none of the original's math needs one. Also note:
  `PCSK9` remains a **shared disease-network state** (also fed by
  statin's own `EFFECT_STA_PCSK9` production term) — it is not renamed to
  an evolocumab-only compartment, the same "combine only at the point of
  use, keep the shared state as disease state" principle already
  established in the `pompe-disease` (`ada_block`, `ag_drive`) and
  `chronic-lymphocytic-leukemia` (`MCL1_ADAPT`, `NK_ACT`) refactors.

## `$PARAM` vs `$CAPTURE` for `C_<STEM>`/`EFFECT_<STEM>`

Confirmed empirically against this exact build (a minimal test model with
a `$PARAM`-declared `C_TEST` reassigned in `$ODE` fails to compile:
`error: assignment of read-only reference 'C_TEST'`) — mrgsolve 2.0.1
compiles `$PARAM` members as read-only references inside `$ODE`, so a
quantity recomputed every timestep from state cannot also live in
`$PARAM`. Per the guide's qspserver compatibility requirement #2 these
would ideally live in `$PARAM`; this is the same, already-disclosed
deviation used throughout this fork (AMD, membranous-nephropathy,
breast-cancer, x-linked-hypophosphatemia, pompe-disease,
chronic-lymphocytic-leukemia): every `C_<STEM>`, `EFFECT_<STEM>`, and the
Hill-fraction intermediates are `double`s computed in `$ODE`, listed in
`$CAPTURE` instead — confirmed present in `/model_manifest`'s
`outputPaths` (99 total, up from the original's 84 — the increase is the
new `EFFECT_*` terms plus the previously-uncaptured `C_NIA`). Every PK
parameter and every new `GAMMA_*`/`EMAX_*_*` constant (none of which are
ever reassigned) *are* real `$PARAM` entries, confirmed present in
`/model_manifest`'s `parameters` list (174 total).

## No pre-existing build defect found

Unlike many other files in this batch, **the original compiles cleanly
under mrgsolve 2.0.1 as-is** — confirmed via `POST /model_manifest` on
the untouched original (200 OK, 84 output paths, 148 parameters). No
syntax-only build-compat fix was needed anywhere in the delivered
`_refactored.R`, and no new `UPSTREAM_ISSUES.md` entry was required for
this file.

## Verification

**Method.** Both the original's own quoted DSL string and the
`_refactored.R`'s embedded DSL were extracted as bare mrgsolve text (the
extracted refactored text confirmed byte-identical to the `_refactored.R`
file's own quoted string — 33,835 characters, exact match) and run
through the qspserver `mrgsolve_api` container (`POST /model_manifest`,
`POST /run_simulation`) at `http://localhost:8007`, requests spaced
~2.5s apart, respecting the API's 2-concurrent-job limit.

Rather than invoking one narrow named scenario from the original file's
own `sc` list per compound, dosing was built as the **superposition of
all six compounds' own real regimens** — the exact `amt`/`cmt`/`ii`
values the original file's own compartment-index comment and `ev_*`
definitions already specify (muvalaplin 240 mg q1d cmt 7, evolocumab
420 mg q28d cmt 9, statin 20 mg q1d cmt 12, niacin 2000 mg q1d cmt 14,
ziltivekimab 30 mg q28d cmt 15, obicetrapib 10 mg q1d cmt 17 — identical
to the original's `ev_muva`/`ev_evol`/`ev_stat`/`ev_niac`/`ev_zilt`/
`ev_obic`), run together over a shortened 90-day window (vs. the file's
own default 2-year `run_scenario` window — ample to exercise every
compound's dosing cadence: 4 doses of the two q28d drugs, ~90 daily doses
of the four q1d drugs). This is a genuine superposition of the file's own
dosing definitions, not an invented dose, and it is *more* thorough than
running scenarios one at a time since it also exercises every
cross-compound shared-state interaction (PCSK9 shared by statin+
evolocumab, LPA transcription shared by statin+niacin+obicetrapib+IL-6,
etc.) in a single run. Two independent runs were performed:

1. **Baseline dosing, default physiology** (`LPA0=250`, `NKIV2=12`,
   `FLDLRFN=1`, `IL6EXO=0` — the file's own defaults).
2. **HoFH + rheumatoid-arthritis background** (`FLDLRFN=0.05`,
   `IL6EXO=18` — the same parameter overrides the original file's own
   scenarios 18/19 use), same six-compound dosing — stress-tests the
   statin/evolocumab/PCSK9/LDLR interaction under near-zero LDL-receptor
   function and the ziltivekimab/IL-6 block under a high inflammatory
   background.

40 shared outputs were compared point-by-point across the full 91-point
time grid in both runs: all 6 renamed PK compartment pairs (`GUT_*`/
`CENT_*`), all 6 renamed concentrations (`C_*`), every shared disease
state (`LPA_P`, `LDL_P`, `VLDL_P`, `PCSK9`, `LDLR`, `MRNA`, `APOA_ER`,
`APOA_FR`, `OXPL`, `MONO`, `IL6`, `CRP`, `HDL_C`), and every downstream
read-out that any of the six compounds feeds (`LPA_NMOL`, `LPA_MASS_T`,
`ASSAY_MASS`, `ASSAY_APOA`, `LDLC_MEAS`, `APOB_TOT`, `HR_MACE`, `IL6EFF`,
`FIL6`, `TRANS`, `KNOCK`, `FMUV`, `KCATL`, `KCATD`, `ASSEM`).

**Result: exact match, max abs diff = 0.0 for every one of the 40 shared
outputs, at every time point, in both scenarios** (both the plain-dosing
run and the HoFH+RA parameter-override run). This includes the terms
that were algebraically re-derived rather than left as a literal
byte-for-byte copy (`FMUV`, `IL6EFF`, `TRANS`) — the guide's tolerance
rule anticipates floating-point-scale drift for these; none was measured
at the precision returned by the API. This is the expected outcome for
six compounds whose refactor is entirely rename-plus-exact-algebraic-
derivation, per the guide's tolerance rule for pure structural
reorganization (Archetypes 1/3, no Hill-fitting anywhere in this file).

`/model_manifest` on the refactored DSL confirms 99 output paths (up
from 84) and 174 parameters (up from 148), including every renamed PK
parameter and every new `GAMMA_*`/`EMAX_*_*` constant for all six
compounds.

## Anything else worth flagging

- **Compartment ordering is fully preserved** — no compartment was
  added, removed, or reordered (`GUT_MUV` is still compartment 7,
  `CENT_OBI` is still compartment 17, etc.), so the R driver script's
  existing `ev(... cmt = 7 ...)`-style dosing events needed **no
  changes at all** — they address compartments by number, not name, and
  the numbering is untouched. Only the compartment/dosing-index comment
  block (a comment, not code) was updated for the six in-scope compounds,
  to avoid the comment silently going stale; the pelacarsen/siRNA lines
  in that same comment block are untouched.
- The compiled model's cache name was changed to `"lpa_qsp_refactored"`
  (from `"lpa_qsp"`) so it does not share a `mcode_cache` compilation
  slot with the original — a standard, value-inert precaution used
  elsewhere in this fork (e.g. `cll_mrgsolve_model_refactored.R`'s
  `"CLL_QSP_refactored"`).
- The 20-scenario `sc` list, `run_scenario()`, the RUN 1-4 driver blocks,
  and the CALIBRATION NOTES block at the end of the file are all
  untouched — none of them reference the renamed compartments by name
  (dosing is by numeric `cmt=`, and the output columns they read —
  `ASSAY_MASS`, `CRP`, `LPA_NMOL`, etc. — are all unrenamed shared
  disease outputs), so they continue to work exactly as before.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`, the
`elevated-lipoprotein-a | EVO`, `| MUV`, `| NIA`, `| OBI`, `| STA`, and
`| ZIL` rows — replacing the compound names that were placeholders in the
existing rows with the real drug names confirmed against the code
(obicetrapib, not olpasiran, for `OBI`; ziltivekimab, not
zilebesiran/pelacarsen, for `ZIL`).
