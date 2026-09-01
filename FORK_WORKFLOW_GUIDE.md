# Fork workflow guide

This repository is a fork of an actively-developed upstream QSP disease-model
library (a new model committed almost daily). This is the single guide for
every kind of transformation this fork makes to that library — currently two:
**translating** a file to English, and **refactoring** a model's PK to be
pluggable with a named PK–PD interface. Read the shared rules below once, then
jump to whichever part applies.

## Shared rules (apply to every transformation in this fork)

**Never edit an upstream-tracked file.** An in-place edit becomes a merge
conflict on every upstream pull, forever. Every addition is a new sibling
file, named `<name>_<suffix>.<ext>` next to the original it derives from
(`_en` for a translation, `_refactored` for a PK rewrite) — never a separate
mirrored tree (see [`ESQlabs_curation_and_extension.md`](ESQlabs_curation_and_extension.md)
for why that was tried and abandoned) and never an edit to the original
itself, which stays the ground truth you verify against.

**Verify against the original, mechanically, every time.** Whatever a
transformation changes, something else must prove nothing unintended changed
alongside it — a token check, a re-run of the model, a diff of untouched
output. Never accept "it looks right" as the gate.

**A gate failure is a question, not a verdict, and it cuts both ways.**
Checking scripts have been wrong before (see the translation gates' bug
history) and an LLM-driven rewrite can be wrong in its own new ways. When
something doesn't match, find out which side is wrong and say so. **Never
reword a translation, or loosen a fit, or adjust a comparison, just to make a
gate pass** — that silently degrades the work and hides the bug from whoever
has to fix it next.

**Log what you find, don't fix it upstream.** A defect discovered in the
original (broken table, wrong count claim, a PK model that's actually
TMDD when a comment implies otherwise) goes in
[`translations/UPSTREAM_ISSUES.md`](translations/UPSTREAM_ISSUES.md), never
into the original file.

**The README is context, not ground truth for structure.** Useful for
compound identity, dosing rationale, and which trials calibrated which arm.
Checked directly and found unreliable for structural/mechanistic detail — one
model's README described its most mechanistically rich compartment only as
"PK compartments," with no mention that the code models full receptor-binding
kinetics. Read the actual code for how something is built.

---

# Part 1 — Translation

Use for any file still containing Korean. Two paths, chosen by what the file
mostly *is*.

## Path A — line-level substitution (structure-dominant files)

Use for `*_references.md`, `*.R`, `*.dot`, `*.py`, `*.txt`, `*.json`.

These files are mostly **not** Korean. A `*_references.md` is English
citations, PMIDs, and URLs with Korean section headers; a `*_shiny_app.R` is R
code with Korean UI labels. Retyping the whole file to translate a quarter of
it risks corrupting the parts that must not change, and there are 29,460
PMIDs in this repository. So only the Korean lines travel:

```bash
# 1. pull out the Korean lines as a numbered TSV
python3 translations/tools/lineio.py extract <path> -o /tmp/work.tsv

# 2. translate column 2 in place, keeping the line number and the tab

# 3. substitute them back into an otherwise byte-identical copy
python3 translations/tools/lineio.py apply <path> /tmp/work.tsv

# 4. confirm
python3 translations/tools/lineio.py verify <path>      # must print OK
```

`apply` refuses to write if a line number is missing or extra, or if a
translated line still contains Hangul. `verify` additionally confirms the
line count is unchanged and that **no line without Korean was altered** — so
a PMID, URL, or expression cannot drift through this path.

What must survive a Path A translation unchanged:

- PMIDs, DOIs, URLs, and PubMed links
- author names, journal names, years, article titles
- R and Python syntax: function names, arguments, operators, string quoting,
  `<-`, `%>%`, `$`, indentation
- Graphviz syntax: node and edge IDs, `->`, attribute names, `\n` inside a
  label, quoting, `shape=`, `fillcolor=`, `label=`
- variable, parameter, and compartment names (`KFI`, `MVG`, `CEM43`, `dTc/dt`, …)
- every number, unit, and comparison operator

Translate only the human-language text: section headings, prose, comments,
UI labels, axis titles, table cell text, and node label text.

**Korean is not always stored as Korean.** A few files write it as unicode
escapes (`\u` followed by four hex digits), which contain no Hangul codepoint
at all. `lineio` decodes those and treats them as Korean, so `extract` hands
you such a line like any other — but on screen it looks like ASCII gibberish.
Decode before translating:

```bash
python3 -c "import sys; print(sys.argv[1].encode().decode('unicode_escape'))" '<the line>'
```

Write plain English in place of the escapes; there is no need to re-escape.
Keep any escape that is *not* Hangul (em dash, middle dot) exactly as it was.

Two traps here, both hit for real:

- **Some of those escaped strings are garbled** — see `UPSTREAM_ISSUES.md`
  entry 17. Decode, work out the intended term, translate that, and say in
  your report which reading you used.
- **If the file already has a translation**, do not hand `apply` a TSV
  containing only the escape lines: it rebuilds the translation from the
  *original*, so every other Korean line would revert. Seed the TSV from the
  existing translation, which is line-aligned with the original, then
  override just the lines you are fixing.

### Situations that come up constantly — handle them the same way every time

**A heading that already carries an English gloss.** Many headings are
written `## 3. 유전적 소인 (Genetic susceptibility)`. Translate the Korean and
**drop the now-duplicated bracketed gloss** — never emit
`## 3. Genetic susceptibility (Genetic susceptibility)`. Where the author's
gloss is narrower than the Korean it glosses, keep the fuller translation
rather than the gloss.

**A table of contents whose anchors are built from Korean headings.**
Translating the heading changes its slug and silently breaks every link
pointing at it. **Regenerate each anchor from the translated heading** — the
slug is the heading lower-cased with markdown stripped, everything except
letters, digits, spaces, hyphens, and underscores removed, spaces turned into
hyphens. Checked by `check_links.py` against the file's own headings.

**A Graphviz label that pairs a Korean title with an English subtitle.**
Cluster and node labels are often written
`<b>6. 고환 축</b><br/><i>Testicular axis</i>` — all on one line. Collapse it:
translating both halves renders visible duplication in the image.

**A Korean passage the author already paraphrased in English underneath.**
The right handling differs by path:
- **Path A: leave both.** `apply`/`verify` depend on line alignment, so the
  translation carries the full version and the author's abridged one. Do not
  merge or trim; mention it in your report.
- **Path B: collapse them into one.** No alignment constraint here, and two
  English versions of the same paragraph reads as an error. Keep the fuller
  text (usually the translated Korean), drop the shorter paraphrase.

**The title case specifically, because it is the commonest.** Most
per-disease READMEs open with two `#` H1 lines — a Korean title and the
author's English title. **Collapse these into a single H1 carrying the union
of both** (abbreviation from the Korean line, expanded name from the English
one). A heading count of `n-1` is the correct outcome here, not drift. This
is **not** the same as an H1 followed by an `###` subtitle (a styling choice)
— there, translate the H1 and leave the subtitle alone.

### Is that Korean string a label or an identifier?

In `.R` and `.py` files a Korean string is sometimes a display label and
sometimes a key the code looks up. Do not decide by eye — **grep for it
before translating it**, and translate it only if every occurrence is a
display site:

- `choices = c("한국어 라벨" = "id")` — translate the label, never the `id`.
- `inputId=` / `outputId=` values: never translate.
- `data.frame` backtick column names are display headers **only** where the
  frame carries `check.names = FALSE` and the name is never subscripted
  anywhere. Without `check.names = FALSE`, `make.names()` replaces a space
  with a dot, so a Korean header that was one word becomes `Two.words` —
  use a single token or underscores instead. `tibble()` does not mangle.
- list names consumed as legend text or a table's first column are labels;
  list names matched against a literal elsewhere are identifiers.
- mrgsolve output column names and `$PARAM`/`$CMT` symbols are always
  identifiers.

### A label whose second line is an English gloss

Plot labels are often `경구 섭취\noral`. Collapse to the single English
label. Keep the `\n` where the second line carries something else (a
variable name, a unit, an abbreviation). The vector's element count must not
change.

### Korean groups numerals differently, and that is fine

Korean counts in 만 (10,000) and 억 (100,000,000), English in thousands and
millions. `130만 mL` is `1.3 million mL`: identical value, different digit
string. Translate naturally and let the token check flag it — record the
review in `translations/data/token_exceptions.tsv`. Do not contort the
English to keep the Korean digits.

### Trailing whitespace is handled for you

Two spaces at end of line is a markdown hard line break. `lineio.py apply`
reattaches each source line's trailing whitespace to its translated line; you
don't need to preserve it by hand, and shouldn't add it by hand either.

### Re-wrapping prose can invent markdown

A Path B translation rewraps every paragraph. Watch what each new line
*starts* with — `+ `, `- `, `15. `, `#` at the start of a wrapped line
silently becomes a list item or heading. Re-wrap so no continuation line
begins with a markdown block marker. Likewise a pipe inside a table cell
splits that cell even inside a code span; write `\|` and run
`check_tables.py --translations` if you touched one.

### Do not tidy notation

Reproduce the author's numeric and symbolic formatting exactly — `−34~−46%`
stays a tilde range, arrows (`→`, `↑`, `↓`) and multiplication signs stay as
written. Numbers and symbols stay byte-identical; only the words around them
change.

## Path B — whole-file translation (prose-dominant files)

Use for the per-disease `README.md`. Dense scientific prose where English
sentences don't line up with Korean sentences, so translate the file as a
document and write it to `<same path>` with `_en` inserted before the
extension — `heat-stroke/README.md` becomes `heat-stroke/README_en.md`.

Preserve the document as a document:

- keep the heading, table, and list structure
- **keep a code fence's *syntax* but translate the Korean *comments* inside
  it**, exactly as Path A does. "Keep code fences" means don't restructure or
  delete them, not "leave every character inside untouched." A run of 30+
  READMEs shared one boilerplate usage block whose comment lines were skipped
  for exactly this misreading. The whole-file Hangul grep must cover fence
  contents too.
- keep every file name, path, link, number, unit, gene name, drug name, and
  parameter name exactly as written, **except** a same-directory sibling's
  own file name once that sibling has its own `_en` translation —
  `source("x.R")` in a translated Shiny app must become `source("x_en.R")`.
  Applies to `source()` calls, `MODEL_FILE`-style path variables, and
  plain-text `## Files` directory listings alike.
- keep relative links working. A translation lives right next to its
  original, so a same-directory link needs no `../` at all — just the bare
  filename (`_en`-suffixed if that target is itself translated). Check
  mechanically: `python3 translations/tools/check_links.py <disease>/README_en.md`
- keep the emphasis the author used

## Register

Scientific writing for a pharmacometrics audience:

- British spelling (`haemoglobin`, `oedema`, `fibre`, `modelling`, `-isation`)
- standard clinical/pharmacometric terms, not literal glosses: 기계론적 지도 is
  "mechanistic map", 구획 is "compartment", 시나리오 is "scenario", 고정점 is
  "fixed point", 되먹임 is "feedback", 쌍안정 is "bistable", 용량-반응 is
  "dose-response"
- keep the author's voice — many files argue a specific structural claim;
  translate the argument, don't summarise or soften it
- do not add content, hedging, or corrections the author didn't write — this
  is a translation

## A hand-fix to a non-Korean line is not durable

`apply` reconstructs every non-Korean line from the pristine original, on
every run. A bare sibling link with zero Korean, hand-edited directly in the
`_en` output, does not survive a second `apply` on that file — it looks
fixed, then a later batch touching the same file's Korean content reverts it,
because `apply` sees the hand-edit as drift to remove. Fix it with
`python3 translations/tools/relink.py --write` instead, and run it **last**,
after every batch in a wave has finished, not mid-wave.

## Coverage and gates

```bash
python3 translations/tools/translation_status.py            # summary table
python3 translations/tools/translation_status.py --check     # residual Hangul
python3 translations/tools/render_maps.py                    # after any .dot
python3 translations/tools/check_tokens.py
python3 translations/tools/check_links.py
python3 translations/tools/check_tables.py --translations
```

A translated `.R` file must still parse (or, if it's pure mrgsolve DSL rather
than an R script, `parse()` failing is expected — test for the block marker
in its three spellings, `$PROB`/`[PROB]`/`[ PROB ]`, before concluding
anything is broken):

```bash
Rscript -e 'cat(length(parse("<path>_en.R")), length(parse("<path>.R")), "\n")'
python3 translations/tools/check_r_parses.py     # classifies every .R file
```

### A translated label can silently change a plot's colour

`ggplot2` orders factor levels alphabetically; an **unnamed**
`scale_*_manual(values = c(...))` assigns colours by that order, so
translating `무치료`/`치료` to `Untreated`/`Treated` can swap two colours
even though the legend still reads correctly. Check whether the vector is
named (translate names, order doesn't matter) or unnamed (choose wording that
preserves sort order where reasonable, and say so in your report when you
can't).

---

# Part 2 — PK/PD refactor (pluggable PK, named Hill interface)

Use when rewriting a model so its PK can later be driven by an external
covariate (see [`driver-patches/data/compound_perturbation_census.md`](driver-patches/data/compound_perturbation_census.md)
for the corpus-wide classification this work is based on). Produces
`<abbr>_mrgsolve_model_refactored.R`, original untouched, same as every other
transformation in this fork.

## What "plumbable PK" means, concretely

A model is plumbable when, for every compound it models:

1. That compound's PK lives in its own clearly-delimited block — its own
   compartments and `dxdt_` lines, not interleaved line-by-line with PD code.
2. It exposes **exactly one** concentration variable that PD equations read
   (two only when a genuinely different tissue site matters). That
   variable's *definition* is the single point where an external covariate
   could later be substituted in.
3. Naming follows the convention below, not whatever the original file used.
   (The corpus has six different names for "the central compartment" across
   17 files modeling the same drug, tocilizumab — `CENT_TCZ`, `TCZ_central`,
   `TCZC`, `A_TCZ1`, `C1_tcz`, `TOCZ_C`. This is the chaos being removed.)
4. The compound's effect on the disease is expressed as one named function
   of that concentration (or occupancy), not buried inside a combined
   multi-drug expression.

**Do not flatten a mechanistically rich PK model into a simpler one because
it's inconvenient.** If the original genuinely models receptor-binding
kinetics (TMDD), the rewrite keeps that — reorganized and renamed, never
less pharmacology than the original had.

## Naming convention

Use exactly this, for every compound, keyed to its existing abbreviation
stem (`TCZ`, `ADA`, `MTX`, …) — don't invent a new stem:

| Role | Pattern | Example |
|---|---|---|
| Absorption depot compartment | `GUT_<STEM>` | `GUT_TCZ` |
| Central compartment | `CENT_<STEM>` | `CENT_TCZ` |
| Peripheral compartment | `PERI_<STEM>` | `PERI_TCZ` |
| Free receptor (TMDD only) | `REC_FREE_<STEM>` | `REC_FREE_TCZ` |
| Drug–receptor complex (TMDD only) | `COMPLEX_<STEM>` | `COMPLEX_TCZ` |
| Clearance | `CL_<STEM>` | `CL_TCZ` |
| Central volume | `V1_<STEM>` | `V1_TCZ` |
| Peripheral volume | `V2_<STEM>` | `V2_TCZ` |
| Inter-compartmental CL | `Q_<STEM>` | `Q_TCZ` |
| Absorption rate | `KA_<STEM>` | `KA_TCZ` |
| Bioavailability | `F_<STEM>` | `F_TCZ` |
| **The exposed concentration** | `C_<STEM>` | `C_TCZ` |
| TMDD: on/off rate, receptor pool | `KON_<STEM>`, `KOFF_<STEM>`, `RTOT_<STEM>` | — |
| Hill: Emax, EC50, Hill coefficient | `EMAX_<STEM>`, `EC50_<STEM>`, `GAMMA_<STEM>` | — |
| The compound's effect on disease | `EFFECT_<STEM>` | `EFFECT_TCZ` |

Parameter *values* always come from the original file — never invent or
default a PK parameter, and don't add one the original didn't have.

## Reference archetypes

Worked examples to standardize toward, **not rigid molds to force every
compound into.** Derived by comparing how tocilizumab alone is modeled across
17 different disease files — the same real drug, six naming conventions,
five genuinely different structures. Match a compound to the archetype
closest to what the original already does; don't add or remove compartments
to make it fit one.

### Archetype 1 — no depot, single compartment, linear elimination

```
$PARAM CL_TCZ = 0.02, V1_TCZ = 5.0
$CMT CENT_TCZ
$ODE
double C_TCZ = CENT_TCZ / V1_TCZ;
dxdt_CENT_TCZ = -CL_TCZ * C_TCZ;
```

### Archetype 2 — no depot, two compartments, linear

Covers both CL/Q/V-style and micro-constant (k10/k12/k21) originals — always
rewrite to CL/Q/V (`k10=CL/V1`, `k12=Q/V1`, `k21=Q/V2`), the more
interpretable, more common convention across this corpus.

```
$PARAM CL_TCZ = 0.5, V1_TCZ = 4.0, V2_TCZ = 3.0, Q_TCZ = 0.6
$CMT CENT_TCZ PERI_TCZ
$ODE
double C_TCZ = CENT_TCZ / V1_TCZ;
dxdt_CENT_TCZ = -(CL_TCZ + Q_TCZ)/V1_TCZ * CENT_TCZ + Q_TCZ/V2_TCZ * PERI_TCZ;
dxdt_PERI_TCZ =  Q_TCZ/V1_TCZ * CENT_TCZ - Q_TCZ/V2_TCZ * PERI_TCZ;
```

### Archetype 3 — depot + central (+ peripheral), linear

```
$PARAM KA_TCZ = 0.1, F_TCZ = 0.8, CL_TCZ = 0.5, V1_TCZ = 4.0,
       V2_TCZ = 3.0, Q_TCZ = 0.6
$CMT GUT_TCZ CENT_TCZ PERI_TCZ
$ODE
double C_TCZ = CENT_TCZ / V1_TCZ;
dxdt_GUT_TCZ  = -KA_TCZ * GUT_TCZ;
dxdt_CENT_TCZ =  KA_TCZ*F_TCZ*GUT_TCZ - (CL_TCZ + Q_TCZ)/V1_TCZ*CENT_TCZ
                 + Q_TCZ/V2_TCZ*PERI_TCZ;
dxdt_PERI_TCZ =  Q_TCZ/V1_TCZ*CENT_TCZ - Q_TCZ/V2_TCZ*PERI_TCZ;
```

Drop `GUT_TCZ`/`KA_TCZ`/`F_TCZ` for a 2-compartment-no-depot variant.

### Archetype 4 — TMDD (target-mediated drug disposition)

Used when the original explicitly models receptor binding as its own dynamic
system, not just an algebraic ratio on concentration. Keep this whenever the
original has it — do not simplify to a plain PK archetype.

```
$PARAM KA_TCZ = 0.1, F_TCZ = 0.8, CL_TCZ = 0.5, V1_TCZ = 4.0,
       V2_TCZ = 3.0, Q_TCZ = 0.6,
       KON_TCZ = 0.01, KOFF_TCZ = 0.001, RTOT_TCZ = 2.0, KDEG_TCZ = 0.05
$CMT GUT_TCZ CENT_TCZ PERI_TCZ REC_FREE_TCZ COMPLEX_TCZ
$MAIN
REC_FREE_TCZ_0 = RTOT_TCZ;
$ODE
double C_TCZ = CENT_TCZ / V1_TCZ;
dxdt_GUT_TCZ    = -KA_TCZ * GUT_TCZ;
dxdt_CENT_TCZ   =  KA_TCZ*F_TCZ*GUT_TCZ - (CL_TCZ + Q_TCZ)/V1_TCZ*CENT_TCZ
                   + Q_TCZ/V2_TCZ*PERI_TCZ
                   - KON_TCZ*C_TCZ*REC_FREE_TCZ + KOFF_TCZ*COMPLEX_TCZ;
dxdt_PERI_TCZ   =  Q_TCZ/V1_TCZ*CENT_TCZ - Q_TCZ/V2_TCZ*PERI_TCZ;
dxdt_REC_FREE_TCZ = KDEG_TCZ*(RTOT_TCZ - REC_FREE_TCZ)
                    - KON_TCZ*C_TCZ*REC_FREE_TCZ + KOFF_TCZ*COMPLEX_TCZ;
dxdt_COMPLEX_TCZ  = KON_TCZ*C_TCZ*REC_FREE_TCZ - KOFF_TCZ*COMPLEX_TCZ
                    - KDEG_TCZ*COMPLEX_TCZ;
```

Here the Hill interface (below) is built from `COMPLEX_TCZ/RTOT_TCZ`
(fractional receptor occupancy), not from `C_TCZ` directly.

### Edge case — continuous infusion input, not an event

Some models dose via a zero-order input variable added directly into the
ODE, not an `ev()`/`data_set()` event. Preserve this if that's what the
original does — converting it to depot-style dosing would change how the
model's own scenarios must be run.

### None of these fit

If a compound's PK genuinely doesn't resemble any archetype, don't force it
into one. Rename to the convention, isolate it from PD, expose `C_<STEM>`,
and note in the refactor notes that this compound needed a bespoke structure
and why. A clean, non-standard structure beats a standard structure that's
wrong.

### Keep a calculation in the block the original used it in

Don't move a concentration/effect calculation from `$MAIN` into `$ODE` (or
the reverse) just because an archetype's worked example above shows it one
way — check where the *original* computed it and keep it there. mrgsolve
evaluates `$MAIN` once per reporting/dosing interval, using state from the
*start* of that interval; `$ODE` gets re-evaluated on every internal solver
substep. A value that depends on live state (most effect terms do) behaves
differently depending on which block it's declared in — moving it changes
the actual simulated trajectory, not just where the code lives. This has
independently bitten two refactors in this fork (`idiopathic-pulmonary-
fibrosis`, `hashimoto-thyroiditis`): both saw a real, verification-failing
divergence (up to ~35% on one output) from moving a `$MAIN`-placed
calculation into `$ODE`, and both were fixed by putting it back exactly
where the original had it. Treat "which block" as part of the original's
behavior to preserve, not an implementation detail free to normalize.

### The dose-instant reporting artifact — there is a real fix, not just a disclosure

Several refactors have hit the same qspserver/mrgsolve quirk: a `C_<STEM>`/
`EFFECT_<STEM>` declared as a `double` local inside `$ODE` reads a stale
(pre-dose) value on the duplicate report row `/run_simulation` emits at the
exact instant of a dose — because that row is generated from state that
hasn't yet had the dose applied when the `$ODE`-local is (re)computed. Most
refactors so far have just disclosed this as a cosmetic, self-healing
artifact (it disappears by the next timestep). `clostridioides-difficile-
infection` found the actual fix: declare `C_<STEM>`/`EFFECT_<STEM>` as
`$GLOBAL` macros (matching the original's own pre-existing convention for
similar derived quantities) instead of `$ODE` doubles — confirmed via a
minimal reproduction that this eliminates the artifact entirely rather than
just shrinking its window. Prefer this over disclosing-and-moving-on when
the original already has a `$GLOBAL`-macro precedent to match.

## The Hill interface

Every compound's effect on the disease system is one named variable:

```
double EFFECT_TCZ = EMAX_TCZ * pow(X_TCZ, GAMMA_TCZ)
                    / (pow(EC50_TCZ, GAMMA_TCZ) + pow(X_TCZ, GAMMA_TCZ));
```
where `X_TCZ` is `C_TCZ` for a plain-PK compound, or `COMPLEX_TCZ/RTOT_TCZ`
(fractional occupancy) for a TMDD compound.

**If the original's effect term is already this shape** (a plain ratio like
`C/(C+IC50)`): this is a rename, not a refit. Pull out `EMAX_<STEM>`,
`EC50_<STEM>`, `GAMMA_<STEM>` as explicit named parameters with the
original's values (`GAMMA_<STEM> = 1` if the original had no explicit Hill
coefficient).

**If the original's effect emerges from ODE-solved kinetics** (TMDD, or
anything you can't write down as one line): fit an approximation.
1. Hold the original model's PK/receptor system at steady state (or run long
   enough to approach it) across a range of concentrations spanning what the
   model's own scenarios actually produce — 15–20 log-spaced points is enough.
2. Record the resulting occupancy/effect at each point.
3. Fit `Emax·Xᵞ/(EC50ᵞ+Xᵞ)` (R's `nls()`; start `γ=1`, only let it vary if
   fixing it produces a visibly worse fit).
4. Report R² (or max relative residual) in the refactor notes. **A poor fit
   is a finding, not a bug to hide** — if you can't get a reasonable fit,
   say so and leave the compound's native calculation in place instead of
   forcing a bad approximation into the rewrite.

**Multiple drugs, one pathway**: keep each compound's `EFFECT_<STEM>`
separate; combine them only at the point the disease equations actually use
them. Never collapse several drugs into one shared Hill term — that defeats
independent driveability.

## When the original doesn't compile at all

Several files in this repo have pre-existing mrgsolve-2.0.1 build defects
totally unrelated to any compound's PK (bad `$PARAM` annotations, `$CMT`+
`$INIT` jointly redeclaring compartments, `$CAPTURE` duplicating compartment
names, `TIME`/`_init_<CMT>` idioms this build no longer accepts, name
collisions across `$PARAM`/`$ODE`/`$TABLE`, and more — see
`translations/UPSTREAM_ISSUES.md` #27 onward for the running list). This
comes up often enough at this repo's scale that it needs one settled answer,
not a fresh judgment call per file:

1. **Never edit the checked-in original** to work around it — that's still
   the shared never-edit-upstream rule. Log the defect as a new numbered
   entry in `UPSTREAM_ISSUES.md` (re-check the file's tail immediately
   before appending, since other agents may be writing concurrently — a
   numbering gap from a lost race is fine, see the note at the top of that
   file; don't let it block you).
2. **The delivered `_refactored.R` should still actually compile and run**
   through the qspserver API, because a "verified" deliverable nobody can
   run afterwards defeats the point of doing this refactor at all. So:
   apply the same **syntax-only, non-numeric, non-behavioral** fix that
   verification needed (the `<cmt>_0` idiom, deduplicating `$CAPTURE`,
   renaming a colliding declaration, adding a missing annotation field,
   etc.) directly into the `_refactored.R` sibling, not just into a
   throwaway scratch copy. This is a fork-owned new file, not the upstream
   original — giving it working syntax is not "fixing it upstream," and
   it's a normal, disclosed part of the refactor rather than a smuggled,
   silent one.
3. **Disclose it plainly in `_refactor_notes.md`**: name the exact defect,
   the exact syntax change made to work around it, and state explicitly
   that it changes nothing numeric — cite the verification result as proof.
   A reviewer should be able to tell, from the notes alone, that this was a
   build-compatibility fix and not a scope-creeping edit to the compound's
   own PK/PD.
4. If a defect is genuinely inside the scope compound's own block (not
   incidental), that's different — fix it as part of the refactor itself
   and say so, same as any other archetype decision.

(Earlier deliverables in this batch — the initial 5-model calibration set
plus part of the first scale-out batch — predate this settled answer and
still carry the original's build defect forward unfixed in their own
`_refactored.R`; each says so explicitly in its own notes, so the fix
needed is fully recoverable later without re-doing the verification work.
No need to redo them retroactively; apply this going forward.)

## Verification (mandatory, do not skip)

1. Run every dosing scenario already defined in the *original* file's own R
   code (not invented ones) through both the original and the
   `_refactored.R`, identical dosing.
2. Compare every `$CAPTURE`d output across the full time grid.
3. Tolerance:
   - Archetypes 1–3 and the edge case (pure structural reorganization, no
     Hill-fitting): expect a near-exact match. Same math, reorganized;
     anything beyond floating-point-scale deviation means a bug, not a
     tolerance to loosen.
   - Archetype 4 and any other ODE-derived interface: a stated, non-zero
     tolerance is expected — report the deviation and the fit quality
     together so a reviewer can judge the approximation.
4. A model that doesn't verify does not get committed. Report the mismatch;
   don't adjust the comparison to make it pass.

## Output artifacts

- `<abbr>_mrgsolve_model_refactored.R` — the rewritten sibling.
- `<abbr>_refactor_notes.md` — per model: archetype used per compound (or
  "bespoke" with justification), Hill-fit quality where applicable,
  verification result, anything flagged for manual review.
- Record the outcome against that compound's existing row in
  `driver-patches/data/compound_perturbation_census.md`, rather than
  starting a second parallel tracker.

## Structural template (mandatory) — a file that doesn't parse is invisible to everyone

Every `_refactored.R` **must** follow this exact shape. Not a stylistic
preference: a corpus-wide scan of all 106 refactored files found **69 that
failed to even parse as R** — not a qspserver issue, not a verification
issue, a file that `Rscript`/`parse()` chokes on before any tooling (this
guide's own checks, the qspserver API, the driver-PK dashboard, anything
written against this corpus in the future) can read it at all. Every one of
those 69 was recoverable by aligning to the template below; none needed any
change to the model's actual math.

```r
<stem>_code <- '
$PROB
...
$PARAM @annotated
...
$CMT @annotated
...
$MAIN
...
$ODE
...
$TABLE
...
$CAPTURE @annotated
...
'

mod <- mcode("<name>", <stem>_code)
```

**Non-negotiable structural rules:**

1. **Exactly one top-level string assignment**, `<stem>_code <- '`, opened
   with a single straight quote on its own statement, holding the *entire*
   mrgsolve spec ($PROB through $CAPTURE) as one contiguous R string.
   **Never** write `$BLOCK`/`[BLOCK]` markers as bare, unquoted text
   directly in the `.R` file — that is not valid R, full stop, regardless
   of whether mrgsolve's own `mread()` can load such a file directly from
   disk. Two files in this corpus (`stevens-johnson-syndrome-ten`,
   `thyroid-eye-disease`) used bracket-style `[BLOCK]` markers this way and
   were unreadable by anything that runs the file through R first. Four
   more (`distal-renal-tubular-acidosis`, `myopia-progression`,
   `hereditary-spherocytosis`, `mitral-regurgitation`) used `$BLOCK`
   markers correctly but never wrapped them in a string at all — same
   defect, same fix: wrap the whole spec in `<stem>_code <- '...'`.
2. **The string closes with a lone `'` on its own line**, immediately
   followed (after only blank lines/comments) by the compile call. A
   trailing stray `'` with nothing to close (a copy-paste artifact, seen
   live in `sarcoidosis`, sitting alone at end-of-file) produces an
   unterminated-string error that can look unrelated to its real cause —
   delete it, don't work around it.
3. **Exactly one compile call, after the string closes**: `mod <-
   mcode("<name>", <stem>_code)`. `mread(...)`, `mread_cache(...)`, and
   `mcode_cache(...)` are equally acceptable (some already-committed files
   use each of these), and namespace-qualified (`mrgsolve::mcode_cache(...)`)
   is fine — but the call must exist as real, live code (not only inside a
   `##`/`#`-commented usage example) and its code argument must be the
   string variable from rule 1, not a `...`-forwarded placeholder or a
   trailing named argument like `quiet = TRUE` sitting after it positionally.
4. **No straight apostrophe (`'`) anywhere inside the string body.** This
   is the single most common defect found (the root cause behind the
   69-file failure count above). It hides in more places than `//`
   comments:
   - `//` and `##` line comments (both conventions appear in this corpus)
   - multi-line `/* ... */` block comments
   - `$PROB`'s free-text prologue (no comment marker required there at
     all — the whole block is prose)
   - the free-text description column of a `$PARAM @annotated` /
     `$CMT @annotated` / `$CAPTURE @annotated` line (`NAME : value :
     description`, or `NAME : description` for $CMT/$CAPTURE) — and that
     description can itself contain a colon in prose (`"...same as in the
     original: there is no GUT_DAS..."`), so don't assume the *last* colon
     on the line is the delimiter
   Use a curly apostrophe (`'`, U+2019) instead — reads identically,
   never needs escaping, never breaks the string. If a straight apostrophe
   is already backslash-escaped (`\'`) somewhere, leave it exactly as is —
   it is already valid; do not "fix" it, that produces the invalid escape
   `\'` instead (caught live, twice, from an earlier over-eager fix pass).
5. **`/* */` is fine *inside* the string** (mrgsolve compiles that block as
   C++) but **must never appear in plain R code outside the string** — R
   has no block-comment syntax at all. Use `#` there (caught live in
   `takayasu-arteritis`: a helper R function used `/* 2h infusion */` after
   a live `ev(...)` call and failed to parse).
6. **Verify structurally, not just by eye, before considering a file
   done**: `Rscript -e 'parse("<abbr>_mrgsolve_model_refactored.R")'` must
   succeed with no error. This is necessary but not sufficient — it does
   not confirm the model *compiles* through mrgsolve, only that R can read
   the file at all. Treat a `parse()` failure as a hard blocker, the same
   tier as a failed verification run.

## What makes a compound's PK *discoverable* by downstream tooling

The naming convention above (`C_<STEM>`, `EC50_<STEM>`, `EFFECT_<STEM>`)
isn't only a readability convention — it's a **machine contract**. Tooling
built against this corpus (the driver-PK dashboard, and anything else
written after it) discovers which compounds a model exposes by pattern-
matching the source text directly, not by parsing C++ semantics. Concretely,
a compound is discovered only when **both** of these are literally present
in the string body:

- a statement matching `double C_<STEM> = <expr>;` — **one contiguous
  initializing statement**. A bare forward declaration that assigns the
  value later, `double C_MIT, OCC_MIT, EFFECT_MIT_ATP;` followed by
  `C_MIT = ...;` on its own line further down (seen live in
  `hereditary-spherocytosis`), is **not** matched — the compound becomes
  invisible to discovery even though the model itself is completely valid
  mrgsolve and compiles fine. Always write the exposed concentration as a
  single `double C_<STEM> = <expr>;` line, per the archetypes above.
- a same-stem `EC50_<STEM>` parameter appearing anywhere in the file
  (normally in `$PARAM`). Both must share the exact same `<STEM>` — this is
  what pairs a compound's exposed concentration with its potency parameter
  for automated discovery.

A model with no compound matching both conditions is not a bug — plenty of
disease models genuinely model no exogenous drug PK, or use it for
constructs that don't fit the Hill interface (a gene-therapy or antisense
construct has no meaningful "EC50", for instance) — but if a compound
*does* have both a concentration and a Hill-shaped potency parameter, write
them so this pattern actually finds them. When in doubt, grep the finished
`_refactored.R` yourself: `double C_<STEM> = ` and `EC50_<STEM>` should
both produce a hit for every compound meant to be driveable.

The compile call (rule 3 above) is discovered the same way: a call to
`mcode`/`mread`/`mread_cache`/`mcode_cache`, optionally namespace-qualified,
found anywhere in the file (it does not need to be assigned to a variable
named `mod`, or assigned at all) — but it does need to actually run, and
its code argument does need to resolve to the string from rule 1.

## qspserver compatibility requirements

A refactored model isn't finished just because it verifies locally against
the original — it also has to be usable through qspserver's mrgsolve REST
API (`app-mrgsolve-api`, service `mrgsolve_api` in
`qspserver/docker-compose.yaml`), since that API is the intended way this
fork's models actually get driven and simulated (no local R/mrgsolve
install required on the client side). That API's contract imposes a few
requirements a `_refactored.R` sibling must satisfy, beyond verification:

1. **`model_content` must be pure mrgsolve DSL text — no R wrapper.**
   `POST /run_simulation` takes `model_content` as an inline string that
   goes straight to `mrgsolve::mcode_cache()` server-side; it is not an R
   script. Every `<abbr>_mrgsolve_model.R` (and its `_refactored.R` sibling)
   is `library(mrgsolve); code <- '<DSL block>'; mod <- mcode(...)` — the
   DSL block is a quoted R string inside a larger script, not a standalone
   file. Producing valid `model_content` therefore needs one more
   mechanical step: extract the quoted block into its own
   `<abbr>_mrgsolve_model_refactored.cpp` (bare DSL, no R syntax around
   it). This is a straight text extraction (find the `code <- '...'`
   assignment, pull the quoted contents, save verbatim) — it should not
   change a single character of the DSL itself; if the extraction and the
   `_refactored.R`'s inline copy ever disagree, that is a bug in the
   extraction step, not a reason to hand-edit the `.cpp`.
2. **Exposed covariates must be discoverable via `/model_manifest`.**
   That endpoint introspects a model's `$PARAM`/`[PARAM]` block to list
   parameter names and defaults. Every `C_<STEM>` and `EFFECT_<STEM>` this
   guide's naming convention creates must live in `$PARAM` (with its `= 0`
   safe default), not be computed only in `$MAIN`/`$ODE` — otherwise a
   client has no way to learn the covariate's name except by reading the
   DSL text directly.
3. **`$SET end/delta` stays, but must not be load-bearing for API runs.**
   Keep it — it's what makes the model runnable standalone (matching the
   original's own behavior) — but `_runner.R`'s `.mrg_dataset()` builds the
   covariate dataset from the request's own `drivers` + `time` fields and
   does **not** consult the model's `$SET` when a `drivers` payload is
   present. A refactored model must produce correct output when time/data
   are supplied externally through the request, not only when its own
   `$SET` values are used. Don't rely on `$SET`'s `end`/`delta` to be
   respected by an API-driven run.
4. **Every output needed for verification or discovery must be in
   `$CAPTURE`.** The API's `outputs: [...]` selection can only return
   columns the model actually captures. If a compound's `C_<STEM>` or
   `EFFECT_<STEM>` (or the disease model's own primary readout) isn't
   captured, it's invisible to an API caller even though it exists inside
   the ODE system.
5. **No R-only syntax inside the DSL block.** Since `model_content` is
   compiled directly server-side with no R wrapper around it, anything
   that depends on the surrounding R script (e.g. R-side helper functions
   defined outside the quoted block and called from `$MAIN`/`$ODE`, or
   R constants interpolated into the string via `sprintf`/`glue` before
   quoting) will not exist at compile time through the API. Keep every
   value the DSL block needs — parameters, macros, helper functions —
   inside the block itself (`$GLOBAL`/`#define`/`$PARAM`), not stitched in
   from the R side.

None of this changes the verification discipline above — the `.cpp`
extraction is checked by confirming the extracted text is byte-identical
to the quoted block in `_refactored.R`, and calling `/model_manifest`
once per finished model as build-time insurance that every covariate this
guide's naming convention promises is actually declared in `$PARAM` — not
a separate parallel verification pass, just one more assertion alongside
the existing ones.

**Known API limitation: solver step count.** `/run_simulation` runs with a
default `maxsteps` well below what some of these models' own drivers
request (one file's own scenario needed `maxsteps=5e6`; the API's default
is 20000). If a scenario times out or errors on step count, shorten the
verification window (fewer days/hours, same dosing) rather than declaring
a mismatch — this is a solver-budget limit, not evidence the refactor is
wrong. Note the shortened window honestly in the `_refactor_notes.md`
rather than silently using a different scenario than the one named.
