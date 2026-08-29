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
