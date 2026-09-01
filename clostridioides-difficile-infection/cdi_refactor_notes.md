# Refactor notes — `clostridioides-difficile-infection/cdi_mrgsolve_model.R` / vancomycin (VAN) and bezlotoxumab (BEZ)

**Scope of this refactor**: only vancomycin's own three GI-transit
compartments and bezlotoxumab's own three PK compartments (plus the
bezlotoxumab-TcdB complex state) were touched. Fidaxomicin/OP-1118,
metronidazole, rifaximin, ridinilazole, the index (precipitating)
antibiotic, and the live biotherapeutic/FMT arm are copied verbatim into
`cdi_mrgsolve_model_refactored.R` — same compartments, same parameter
names, same equations. The entire ecology/pathogen/host disease network
(microbiota guilds, bile-acid metabolism, the nutrient niche, the *C.
difficile* life cycle, PaLoc/toxin regulation minus the two renames noted
below, epithelium/barrier, immunity, and the clinical read-outs) is also
untouched.

## What the census's "Antibacterial" row actually is

`driver-patches/data/compound_perturbation_census.md` carried one row for
this file, `clostridioides-difficile-infection | Antibacterial |
Redirect concentration inside #define macro`. "Antibacterial" is not a
drug name anywhere in the model — it is the literal text of a `$PARAM`
section comment, `// ---------------- antibacterial PK
--------------------------------------`, which groups **six** distinct
real antibacterials under one heading: vancomycin, fidaxomicin (+ its
active metabolite OP-1118), metronidazole, rifaximin, ridinilazole, and
an abstract "index (precipitating) antibiotic". The census's classifier
evidently latched onto that section comment rather than a per-compound
name, so the single row cannot be mechanically resolved to one compound.

**Vancomycin (VAN) was selected as the drug behind that row**, on the
strength of several independent, consistent signals in the file itself
(per the guide's "read the actual code" rule, since the README/census is
not ground truth for this):

- It is the only antibacterial discussed by name in the model's own
  opening docstring ("Oral vancomycin sterilises the stool of vegetative
  *C. difficile* within days and cures ~80% of episodes...", the model's
  entire framing device).
- It is listed **first** in both places the six antibacterials are
  grouped (the `antibacterial PK` `$PARAM` section, and the `drug PD`
  section's kill-rate constants).
- It is used in more scenarios than any other antibacterial (8 of the
  file's 18: S06, S09, S11, S12, S13, S15, S16, S18), and scenario S06 is
  explicitly annotated `"standard of care first episode"`.

This is a reasoned selection, not a certainty the code can prove on its
own — fidaxomicin or metronidazole (the two alternatives the task brief
itself flagged as plausible) remain equally valid candidates for a
**separate**, future census row/refactor pass; they are left completely
untouched here rather than folded in speculatively.

## Bezlotoxumab (BEZ): the census classification is corrected here too

The census classified BEZ as **"Delete PK compartment; concentration is
itself the state."** That is half right and half misleading:

- `BEZ_GUT` (renamed `LUM_BEZ`) genuinely *is* tracked directly as a
  concentration (mg/L = ug/mL), with no volume division anywhere it is
  used (`ABEQ = pos(AB_IGG) + WBEZ*pos(BEZ_GUT)`, `BIND =
  KONB*pos(TCDB)*pos(BEZ_GUT)`) — so far, consistent with the census
  label.
- But that concentration is **not itself a dosed or independent state**.
  It is fed by a real, ordinary 2-compartment linear PK backbone
  (`BEZ_C`/`BEZ_P`, central/peripheral, `CLB`/`VCB`/`VPB`/`QB`) via
  first-order transudation, and it is **consumed by genuine
  target-mediated mass-action binding** to free TcdB (`BIND =
  KONB*pos(TCDB)*pos(BEZ_GUT)`, forming `TOX_CPLX` with an explicit
  off-rate `KOFFB` and complex clearance `KDEGCX`, and — disclosed in the
  original's own comment — a (deliberately negligible) stoichiometric
  drug loss term in `dxdt_BEZ_GUT` itself). That is TMDD-shaped kinetics,
  the opposite of "delete the PK compartment", and per the guide ("If the
  original genuinely models receptor-binding kinetics (TMDD)... keep
  that") it is kept, not simplified away.

So this is **not** a clean instance of any of the guide's four
archetypes: it is Archetype 2 (central/peripheral, no absorption depot —
the dose is a direct IV bolus into `CENT_BEZ`) **plus** a bespoke third,
genuinely distinct luminal/mucosal site (the guide's own "two
concentrations only when a genuinely different tissue site matters"
exception) **plus** real, disclosed, unmodified target-mediated drug
consumption at that third site. Recorded back into the census as
"Redirect concentration; bespoke 3-site structure (central/peripheral PK
+ luminal transudation site with real TMDD-style target consumption),"
not "Delete PK compartment."

One further nuance on the TMDD side: the guide's Archetype 4 template
assumes an isolated `REC_FREE_<STEM>`/`COMPLEX_<STEM>` pair with its own
conserved-total turnover. Here the "receptor" bezlotoxumab binds is
`TCDB` — the disease's own shared toxin state, independently produced
and degraded by C. difficile toxin regulation whether or not any
bezlotoxumab is present. Renaming `TCDB` itself to a bezlotoxumab-scoped
name would be wrong (it is read/written by many disease equations
unrelated to bezlotoxumab). Only `TOX_CPLX`, which exists *only* when
bezlotoxumab is present, is renamed to `COMPLEX_BEZ`; `TCDB` is left
exactly as it was. This is disclosed as a bespoke partial application of
Archetype 4, not a clean instance of it.

## Vancomycin: bespoke structure, not a clean 1-3 archetype either

Vancomycin's three original compartments (`VAN_D` -> `VAN_T` -> `VAN_COL`)
are **not** a depot-into-systemic-circulation chain — oral vancomycin for
CDI has essentially zero systemic absorption, and all three states are GI
lumen/content, ending at the colonic content whose concentration
(`CVAN = 1000*VAN_COL/VCOLG`, ug/g faeces) is what every PD term reads.
It is a three-stage first-order GI-transit chain sharing **one** rate
constant (`KTRV`, renamed `KTR_VAN`) across both transitions, terminating
in first-order faecal elimination (`KEXV`, renamed `KE_VAN`) — closest in
spirit to Archetype 3 (depot + "central") but with an extra transit stage
and no absorption into a true central/blood compartment at all. Renamed
per the guide's own naming spirit rather than forced into the fixed
four-role table:

| Original | Refactored | Role |
|---|---|---|
| `VAN_D` | `GUT_VAN` | dosing/absorption depot (unchanged role) |
| `VAN_T` | `GUT2_VAN` | **bespoke**: 2nd-stage GI transit compartment — no fixed-table role fits a drug with no systemic compartment at all |
| `VAN_COL` | `CENT_VAN` | the exposed, PD-facing site — "central" here means "site of action" (colonic content), not plasma |

`VCOLG` (colonic-content mass used for the faecal-concentration
conversion) and `KWASH` (diarrhoeal-washout multiplier) are **shared**
across every antibacterial in this file (fidaxomicin, metronidazole,
rifaximin, ridinilazole all use them too) and are deliberately left
unrenamed so the five out-of-scope compounds' equations stay untouched.

## Renaming tables

**Vancomycin PK** (`$PARAM`, `antibacterial PK` section):

| Original | Refactored |
|---|---|
| `KTRV` | `KTR_VAN` |
| `KEXV` | `KE_VAN` |
| `EC50V` | `EC50_VAN` |
| — (new) | `EMAX_VAN = 1.0` |
| — (new) | `GAMMA_VAN = 1.0` |

`EMAXV` (vancomycin's own maximal *C. difficile* kill-rate, 14.0/d) is
**intentionally left unrenamed**. It is a disease-side pathway weight —
exactly analogous to the untouched `EMAXF`/`EMAXM`/`EMAXD`/`EMAXR`/`EMAXT`
weights for the other, out-of-scope antibacterials, and to
`KV_SBA`/`KV_BUT`/`KV_BAC`/`KV_BIF`/`KV_ENC` (vancomycin's own per-guild
collateral-damage weights, also unrenamed for the same reason) — not part
of vancomycin's own PK/Hill-interface block per the guide's role table.

**Bezlotoxumab PK** (`$PARAM`, `bezlotoxumab PK` section and the toxin
section's TMDD constants):

| Original | Refactored |
|---|---|
| `CLB` | `CL_BEZ` |
| `VCB` | `V1_BEZ` |
| `VPB` | `V2_BEZ` |
| `QB` | `Q_BEZ` |
| `KTRAB` | `KTR_BEZ_LUM` |
| `KDGUT` | `KE_BEZ_LUM` |
| `KONB` | `KON_BEZ` |
| `KOFFB` | `KOFF_BEZ` |
| `KDEGCX` | `KDEG_BEZ` |
| — (new) | `EMAX_BEZ = 1.0` |
| — (new) | `GAMMA_BEZ = 1.0` |

`WBEZ` and `KIAB` are left unrenamed: both are disease-side constants
combining bezlotoxumab with the *endogenous* anti-toxin IgG response
(`ABEQ = pos(AB_IGG) + WBEZ*pos(LUM_BEZ)`), not bezlotoxumab's own PK.

**Compartments**: `VAN_D`->`GUT_VAN`, `VAN_T`->`GUT2_VAN`,
`VAN_COL`->`CENT_VAN`, `BEZ_C`->`CENT_BEZ`, `BEZ_P`->`PERI_BEZ`,
`BEZ_GUT`->`LUM_BEZ`, `TOX_CPLX`->`COMPLEX_BEZ`. `TCDA`, `TCDB`,
`TCDA_MUC`, `TCDB_MUC`, `CDT`, and every other compartment are unchanged.

## The Hill interface

**Vancomycin** — the original's `FV` macro (`CVAN/(EC50V+CVAN)`) is
already exactly the guide's canonical shape with an implicit `Emax=1,
gamma=1`. This is a rename/promotion, not a refit:

```
#define C_VAN      (1000.0 * CENT_VAN / VCOLG)
#define EFFECT_VAN (EMAX_VAN * pow(C_VAN, GAMMA_VAN) / (pow(EC50_VAN, GAMMA_VAN) + pow(C_VAN, GAMMA_VAN)))
```

`EFFECT_VAN` is substituted for every occurrence of `FV` in the disease
equations (`KILL_SBA`/`KILL_BUT`/`KILL_BAC`/`KILL_BIF`/`KILL_ENC`,
`ABXPRESS`, `KILLCD`) — same downstream weighting constants
(`KV_SBA`, ..., `EMAXV`), same math, only the name changed.

**Bezlotoxumab** — no single Emax/EC50 term exists in the original for
"bezlotoxumab's effect"; its real mechanism is the mass-action `BIND`
consumption of `TCDB`, kept completely unmodified. For discoverability, a
Hill-shaped `EFFECT_BEZ` was additionally derived by exact algebra from
the *existing* neutralization term, not fit to anything new. The original
computes a combined neutralizing-capacity fraction

```
double ABEQ = pos(AB_IGG) + WBEZ*pos(BEZ_GUT);
double fnAB = 1.0 / (1.0 + ABEQ/KIAB);
```

whose complement, `1 - fnAB = ABEQ/(KIAB+ABEQ)`, is already Hill-shaped
(`Emax=1, gamma=1`). Isolating just the bezlotoxumab-only contribution
(`WBEZ*LUM_BEZ`) from that combined term and re-deriving its own implied
EC50 by simple algebra (`X/((KIAB/WBEZ)+X) = (WBEZ*X)/(KIAB+WBEZ*X)`)
gives an exact, not fitted, rename:

```
#define C_BEZ      (pos(LUM_BEZ))
#define EC50_BEZ   (KIAB / WBEZ)
#define EFFECT_BEZ (EMAX_BEZ * pow(C_BEZ, GAMMA_BEZ) / (pow(EC50_BEZ, GAMMA_BEZ) + pow(C_BEZ, GAMMA_BEZ)))
```

**`EFFECT_BEZ` and `C_BEZ` are purely diagnostic/discoverability
additions** — the original's `ABEQ`/`fnAB`/`BIND`/`dxdt_TCDB`/
`dxdt_COMPLEX_BEZ` math is left **completely unmodified** (only the
renamed identifiers `LUM_BEZ`/`COMPLEX_BEZ`/`KON_BEZ`/`KOFF_BEZ`/
`KDEG_BEZ` are substituted in place). No R² is reported because no curve
was fit to anything: `EFFECT_BEZ` is an exact algebraic restatement of
values the original already computes, and every disease-facing quantity
that actually depends on bezlotoxumab (`ABEQ`, `fnAB`, `BIND`, `ABTOT`,
`RRI`) still reads `LUM_BEZ` directly, unchanged.

`C_BEZ` is deliberately the **luminal** concentration (the guide's "two
concentrations only when a genuinely different tissue site matters"
exception), not plasma — plasma is exposed separately as
`Cp_bez_plasma = CENT_BEZ/V1_BEZ`, informational only, read by no PD
term (mirrors `abdominal-aortic-aneurysm/aaa_refactor_notes.md`'s
`C_DOXY` (tissue) vs `Cp_doxy_plasma` split for the identical reason).

## A found-and-fixed reporting artifact: macros, not `$ODE` doubles, for `C_VAN`/`EFFECT_VAN`/`C_BEZ`/`EFFECT_BEZ`/`Cp_bez_plasma`

Most of this repo's other pluggable-PK refactors expose `C_<STEM>`/
`EFFECT_<STEM>` as plain `$ODE`-scope `double`s (or the `capture`
shorthand, e.g. `toxic-alcohol-poisoning/tap_refactor_notes.md`'s
`C_FOM`). That was the first approach tried here too, and it **compiles
and verifies correctly for every actual disease/PD output** — but a
direct comparison of the diagnostic columns themselves against a live
qspserver `/run_simulation` call surfaced a real, reproducible artifact,
confirmed with a minimal two-line repro model against the API: at the
duplicate report row mrgsolve emits for an instantaneous dose (pre-dose
and post-dose, same timestamp), a compartment's own state is correctly
already post-dose, but an `$ODE`-scope `double` derived from that state
is not recomputed for that extra row — it still holds the pre-dose value,
self-correcting only at the *next* real integration step. This is the
same phenomenon `sepsis/sep_refactor_notes.md` disclosed for `C_TCZ` and
`toxic-alcohol-poisoning/tap_refactor_notes.md` disclosed for `C_FOM`/
`EFFECT_FOM` (both left as a known, disclosed, non-correctness-affecting
cosmetic artifact). It shows up here specifically for bezlotoxumab
because its dose is a direct IV bolus straight into `CENT_BEZ` (no depot
buffers the jump) — vancomycin never triggers it, because its dose enters
`GUT_VAN`, two ODE transit stages upstream of `CENT_VAN`, so `CENT_VAN`
never has a within-timestep discontinuity to be stale about.

Rather than just disclose this as FOM's notes did, it was fixed here:
`C_VAN`, `EFFECT_VAN`, `C_BEZ`, `Cp_bez_plasma`, `EC50_BEZ`, and
`EFFECT_BEZ` are `$GLOBAL` `#define` macros (exactly like the original's
own `CVAN`/`FV`/`CBEZP` were), not `$ODE` doubles. A macro is re-expanded
textually wherever it is referenced, including inside `$TABLE`/`$CAPTURE`
at *their* execution time, so it always sees the current state — the
minimal repro model confirmed this directly (macro and `$TABLE`-scope
double both read the correct post-dose value on the duplicate row; the
`$ODE`-scope double did not). After the fix, **every output matches the
original exactly, including at the dosing-instant duplicate row** (see
Verification below) — there is no remaining artifact in this file, not
just a disclosed one.

## R-side changes

`ev_vanco()`, `ev_vanco_taper()` (`cmt = "VAN_D"` -> `"GUT_VAN"`,
5 occurrences), `ev_bezlo()` (`cmt = "BEZ_C"` -> `"CENT_BEZ"`), and
`RX_CMT` (used by `rx_end()` to find each scenario's last treatment day)
were updated to the renamed compartment names. `CVANo`/`CBEZo` (the
`$TABLE` columns the R post-processing code reads by name, e.g.
`pk$van_ugg <- mean(d$CVANo)`) are kept, now defined as `CVANo = C_VAN`
and `CBEZo = Cp_bez_plasma` — same values, same column names, so no other
R-side code needed to change.

## Build compatibility

`cdi_mrgsolve_model.R` compiles cleanly as-written under mrgsolve 2.0.1
(confirmed via a live `/model_manifest` call against the untouched
original: 240 parameters, 82 output paths, 61 compartments) — **no
pre-existing build defect for this file**, so no `UPSTREAM_ISSUES.md`
entry was needed, and no build-compat fix was required in the delivered
`_refactored.R`.

## qspserver `/model_manifest` discoverability

Confirmed via a live `/model_manifest` call against the refactored DSL:
244 parameters, 87 output paths, 61 compartments (same compartment count
as the original — renaming, not adding/removing states). Every renamed
`$PARAM` (`KTR_VAN`, `KE_VAN`, `EC50_VAN`, `EMAX_VAN`, `GAMMA_VAN`,
`CL_BEZ`, `V1_BEZ`, `V2_BEZ`, `Q_BEZ`, `KTR_BEZ_LUM`, `KE_BEZ_LUM`,
`KON_BEZ`, `KOFF_BEZ`, `KDEG_BEZ`, `EMAX_BEZ`, `GAMMA_BEZ`) appears in
`parameters`; every renamed/new output (`GUT_VAN`, `GUT2_VAN`,
`CENT_VAN`, `CENT_BEZ`, `PERI_BEZ`, `LUM_BEZ`, `COMPLEX_BEZ`, `C_VAN`,
`EFFECT_VAN`, `C_BEZ`, `EFFECT_BEZ`, `Cp_bez_plasma`) appears in
`outputPaths`.

## Verification

Per the guide's mandatory protocol, run via the qspserver `mrgsolve_api`
container (`http://localhost:8007`, `POST /model_manifest` and `POST
/run_simulation`; requests spaced ~2s apart, single-threaded, per the
task's concurrency guidance) — no local R/mrgsolve was used. The bare DSL
was extracted verbatim from each file's own `cdi_code <- '...'` R-string
(a pure text extraction; the extraction script and its output were
scratch files, deleted after use).

Two of the original file's own scenarios were run identically through
both models, end=90 (days), delta=0.25 — well within the API's default
solver-step budget, so no window shortening was needed:

- **S06 "vancomycin 125 qid x10d"** — index antibiotic days 0-6, spore
  inoculum day 2, vancomycin 125 mg qid days 9-19 (40 doses). Exercises
  vancomycin's full dosing chain (`GUT_VAN`->`GUT2_VAN`->`CENT_VAN`) and
  every guild/pathogen kill pathway `EFFECT_VAN` feeds.
- **S09 "vancomycin + bezlotoxumab"** — S06 plus a single 10 mg/kg (700
  mg) bezlotoxumab IV bolus on day 10. Exercises bezlotoxumab's full PK
  (`CENT_BEZ`/`PERI_BEZ`/`LUM_BEZ`), the TMDD-style `BIND`/`COMPLEX_BEZ`
  toxin-capture kinetics, and the recurrence-risk (`RRI`) pathway that
  reads bezlotoxumab's contribution.

**Result: exact match on both scenarios.** All 29 shared `$CAPTURE`
outputs compared at every time point (`CVANo`, `CBEZo`, `LGVEG`, `LGMUC`,
`LGSPL`, `LGSPB`, `LGSBA`, `fANAo`, `SHAN`, `SBA2o`, `TOXLo`, `SEVERE`,
`FULMIN`, `SYMPT`, `RRI`, `MB_SBA`, `CD_VEG`, `CD_MUC`, `TCDB`, `TCDA`,
`TCDB_MUC`, `EPI`, `EPI_PERM`, `IM_NEUT`, `STOOL`, `WBC`, `ALB`, `CRE`,
`AB_IGG`) — **max absolute difference = 0.0** across 364 (S06) and 365
(S09) time points each, including the exact duplicate row at the moment
of the bezlotoxumab bolus (`CBEZo` = 227.2727 in both models at that row,
after the macro fix described above; before the fix this row alone
differed by 227.2727 while every other row and every other output already
matched exactly — see the artifact section above for the diagnosis and
fix). This is the expected outcome for a pure structural reorganization
with an exact (non-fitted) Hill-interface rename, per the guide's
tolerance table.

## Census

Recorded in `driver-patches/data/compound_perturbation_census.md`:

- The `clostridioides-difficile-infection | Antibacterial` row is
  corrected to `clostridioides-difficile-infection | Vancomycin (VAN)`,
  target/pathway `Cell-wall synthesis (D-Ala-D-Ala)`, redirect site
  `C_VAN (faecal colonic concentration, NOT plasma) / EFFECT_VAN`, unit
  `ug/g (faecal)`, with a summary of the bespoke 3-stage GI-transit
  structure and the exact-match verification above.
- The `clostridioides-difficile-infection | Bezlotoxumab` row's
  perturbation method is corrected from "Delete PK compartment;
  concentration is itself the state" to "Redirect concentration; bespoke
  3-site structure (central/peripheral PK + luminal transudation site
  with real TMDD-style target consumption)", target/pathway `TcdB
  (toxin-neutralizing mAb)`, redirect site `C_BEZ (luminal/mucosal, NOT
  plasma) / EFFECT_BEZ; Cp_bez_plasma informational-only`, unit `mg/L
  (luminal)`.

## Discoverability audit (disclosed exception, no code change)

A corpus-wide discoverability audit (session after this refactor) flagged
`C_VAN`/`C_BEZ` as not matching the literal `double C_<STEM> = <expr>;`
pattern downstream tooling pattern-matches for. Investigated: both are
exposed via C preprocessor macros in `$GLOBAL` — `#define C_VAN (1000.0 *
CENT_VAN / VCOLG)` (~line 534) and `#define C_BEZ (pos(LUM_BEZ))` (~line
536) — not `double` statements.

A macro is structurally **incompatible** with also adding a literal
`double C_<STEM> = <expr>;` statement anywhere after the `#define`: the
preprocessor text-substitutes every subsequent occurrence of the token
(including inside a would-be declaration's own name), so
`double C_VAN = ...;` would expand to `double (1000.0 * CENT_VAN /
VCOLG) = ...;` — not valid C++. There is no block-scoped workaround (the
fix used elsewhere in this corpus this session, removing a `$GLOBAL` bare
`double` forward-declare and adding a `$TABLE`-local `double`, doesn't
apply here — the conflict is textual substitution, not a declaration
collision).

Left as-is: `C_VAN`/`C_BEZ` are genuinely, correctly exposed and
`$CAPTURE`d via this macro mechanism — always freshly evaluated at every
point of use, arguably more artifact-resistant than the `double`-based
pattern (no possibility of a stale cross-block read at all). They are
just not discoverable by a naive `double C_<STEM> = ` text grep. Flagged
for whoever owns the downstream discovery tool to decide whether to also
recognize `#define C_<STEM> (...)` as a valid discoverable form.
