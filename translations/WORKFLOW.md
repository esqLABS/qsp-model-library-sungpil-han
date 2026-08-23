# Translation workflow

How a file gets translated in this fork. Read [`README.md`](README.md) first for the
rule about never editing upstream files.

There are two paths, chosen by what the file mostly *is*.

## Path A — line-level substitution (structure-dominant files)

Use for `*_references.md`, `*.R`, `*.dot`, `*.py`, `*.txt`, `*.json`.

These files are mostly **not** Korean. A `*_references.md` is English citations,
PMIDs, and URLs with Korean section headers; a `*_shiny_app.R` is R code with Korean
UI labels. Retyping the whole file to translate a quarter of it risks corrupting the
parts that must not change, and there are 29,460 PMIDs in this repository. So only
the Korean lines travel:

```bash
# 1. pull out the Korean lines as a numbered TSV
python3 translations/tools/lineio.py extract <path> -o /tmp/work.tsv

# 2. translate column 2 in place, keeping the line number and the tab

# 3. substitute them back into an otherwise byte-identical copy
python3 translations/tools/lineio.py apply <path> /tmp/work.tsv

# 4. confirm
python3 translations/tools/lineio.py verify <path>      # must print OK
```

`apply` refuses to write if a line number is missing or extra, or if a translated
line still contains Hangul. `verify` additionally confirms the line count is
unchanged and that **no line without Korean was altered** — so a PMID, URL, or
expression cannot drift through this path.

What must survive a Path A translation unchanged:

- PMIDs, DOIs, URLs, and PubMed links
- author names, journal names, years, article titles
- R and Python syntax: function names, arguments, operators, string quoting,
  `<-`, `%>%`, `$`, indentation
- Graphviz syntax: node and edge IDs, `->`, attribute names, `\n` inside a label,
  quoting, `shape=`, `fillcolor=`, `label=`
- variable, parameter, and compartment names (`KFI`, `MVG`, `CEM43`, `dTc/dt`, …)
- every number, unit, and comparison operator

Translate only the human-language text: section headings, prose, comments, UI
labels, axis titles, table cell text, and node label text.

**Korean is not always stored as Korean.** A few files write it as unicode escapes
(`\u` followed by four hex digits), which contain no Hangul codepoint at all.
`lineio` decodes those and treats them as Korean, so `extract` hands you such a line
like any other — but on screen it looks like ASCII gibberish. Decode before
translating:

```bash
python3 -c "import sys; print(sys.argv[1].encode().decode('unicode_escape'))" '<the line>'
```

Write plain English in place of the escapes; there is no need to re-escape. Keep any
escape that is *not* Hangul (em dash, middle dot) exactly as it was.

Two traps here, both hit for real:

- **Some of those escaped strings are garbled** — see `UPSTREAM_ISSUES.md` entry 17.
  Decode, work out the intended term, translate that, and say in your report which
  reading you used.
- **If the file already has a translation**, do not hand `apply` a TSV containing
  only the escape lines: it rebuilds the translation from the *original*, so every
  other Korean line would revert. Seed the TSV from the existing translation, which
  is line-aligned with the original, then override just the lines you are fixing.

### Two situations that come up constantly — handle them the same way every time

**A heading that already carries an English gloss.** Many headings are written
`## 3. 유전적 소인 (Genetic susceptibility)`. Translate the Korean and **drop the
now-duplicated bracketed gloss** — never emit
`## 3. Genetic susceptibility (Genetic susceptibility)`. Where the author's gloss is
narrower than the Korean it glosses, keep the fuller translation rather than the
gloss: `## 5. 담즙산 합성·되먹임·수송 (Bile acid feedback)` becomes
`## 5. Bile acid synthesis · feedback · transport`.

**A table of contents whose anchors are built from Korean headings.** Several files
open with a TOC of `#anchor` links slugged from their own headings, e.g.
`[3. 이차 지질 축적](#3-이차-지질-축적--리소좀-칼슘)`. Translating the heading changes
its slug and silently breaks every link pointing at it — nothing about the rendered
page looks wrong, the links just stop working. **Regenerate each anchor from the
translated heading.** The slug is the heading lower-cased with markdown stripped,
everything except letters, digits, spaces, hyphens, and underscores removed, and
spaces turned into hyphens — so a dropped `·` leaves a double hyphen. This is
checked: `check_links.py` validates same-document anchors against the file's own
headings.

**A Graphviz label that pairs a Korean title with an English subtitle.** Cluster and
node labels are often written
`<b>6. 고환 축</b><br/><i>Testicular axis</i>` — all on one line. Because it is one
line you may collapse it, and you should: translating both halves renders visible
duplication in the image. Keep the styling, drop the repetition. The same applies to
`ORGANISING THESIS (조직 원리)`, where the parenthetical merely glosses the English
beside it — drop the parenthetical.

**A Korean passage the author already paraphrased in English underneath.** Some
files carry a Korean note followed by the author's own shorter English version of
it. The right handling differs by path, because the constraint differs:

- **Path A: leave both.** `apply` cannot delete lines and `verify` depends on the
  line alignment, so the translation will carry the full version and the author's
  abridged one. Do not merge or trim. Mention it in your report and move on.
- **Path B: collapse them into one.** You are writing the whole document, so there
  is no alignment constraint, and emitting two English versions of the same
  paragraph reads as an error. Keep the fuller text (usually the translated Korean)
  and drop the author's shorter paraphrase. The duplication is an artefact of the
  file being bilingual, not something the author meant a reader to see twice.

**The title case specifically, because it is the commonest.** Most per-disease
READMEs open with two `#` H1 lines on consecutive lines — a Korean title and the
author's English title:

```
# 부신피질암 (ACC) QSP 모델
# Adrenocortical Carcinoma — Quantitative Systems Pharmacology Model
```

**Collapse these into a single H1 carrying the union of both**, here
`# Adrenocortical Carcinoma (ACC) — Quantitative Systems Pharmacology Model` — the
abbreviation from the Korean line, the expanded name from the English one. Do not
emit two H1s that say nearly the same thing, and do not preserve the heading count
for its own sake; a heading count of `n-1` is the correct outcome here and should be
reported as intended rather than as drift.

This is **not** the same as an H1 followed by an `###` subtitle, which is a styling
choice rather than a duplicate. There, translate the H1 and leave the subtitle
alone.

### Is that Korean string a label or an identifier?

In `.R` and `.py` files a Korean string is sometimes a display label and sometimes a
key the code looks up. Translating the second kind silently breaks the app. Do not
decide by eye — **grep for it before translating it**:

```bash
grep -rn '<the korean string>' <the file>
```

Translate it only if every occurrence is a display site. Specific cases in this
library:

- `choices = c("한국어 라벨" = "id")` — translate the label on the left, never the
  `id` on the right.
- `inputId=` / `outputId=` values: never translate.
- `data.frame` backtick column names are display headers **only** where the frame
  carries `check.names = FALSE` and the name is never subscripted (`df$name`,
  `df[["name"]]`) anywhere. Check, then translate.
- **Whether `check.names = FALSE` is present changes what you may write.** Without
  it, `data.frame()` pushes every name through `make.names()`, which replaces a
  space with a dot — so a Korean header that was one word becomes `Two.words` in
  the rendered table. In that case use a single token (`Item`, `Value`, `Metric`,
  `Endpoint`) or underscores, not a multi-word name. Only where
  `check.names = FALSE` is set may a translated header contain spaces. `tibble()`
  does not mangle, so it is unaffected.
- list names consumed as legend text or as a table's first column are labels; list
  names matched against a literal elsewhere are identifiers.
- mrgsolve output column names (`PNVOLWK`, `OCCGC`, …) and `$PARAM`/`$CMT` symbols
  are always identifiers.

### A label whose second line is an English gloss

Plot labels are often `경구 섭취\noral` — Korean, then a lower-case English gloss.
Collapse those to the single English label; do not emit "Oral intake / oral". Keep
the `\n` where the second line carries something else (a variable name, a unit, an
abbreviation). Either way the number of elements in the vector must not change.

### Korean groups numerals differently, and that is fine

Korean counts in 만 (10,000) and 억 (100,000,000), English in thousands and
millions. So `130만 mL` is `1.3 million mL` and `9700만` is `97 million`: the value
is identical and the digit string is not. Translate it into natural English and let
the token check flag it — that flag gets recorded in
`data/token_exceptions.tsv` after review. **Do not** contort the English to keep
the Korean digits, and do not write `1300000`.

Watch the reverse trap too: comma grouping means `12,800` tokenises as `800`. A
read-aloud script that spells numerals in Korean words (`백구십`) should become
English words (`one hundred and ninety`), not digits.

### Trailing whitespace is handled for you

Two spaces at end of line is a markdown hard line break, and several
`*_references.md` files use it to hold each citation's lines together. Editors strip
it silently. You do not need to preserve it: `lineio.py apply` reattaches each
source line's trailing whitespace run to its translated line, and `verify` reports
any that went missing. Do not add trailing whitespace by hand either — it will be
replaced by the source's.

### Re-wrapping prose can invent markdown

English runs longer than Korean, so a Path B translation rewraps every paragraph.
Watch what each new line now *starts* with: a wrapped line beginning `+ `, `- `,
`15. `, or `#` becomes a list item or a heading, silently changing the document
structure. Re-wrap so no continuation line begins with a markdown block marker.

Likewise, a pipe inside a table cell splits that cell even inside a code span --
GFM splits on pipes before parsing inline code. Write `\|` in a table cell, and run
`python3 translations/tools/check_tables.py --translations` if you touched one.

### Do not tidy notation

Reproduce the author's numeric and symbolic formatting exactly, even where a style
guide would change it: keep `−34~−46%` as a tilde range rather than normalising it
to an en dash, keep their arrows (`→`, `↑`, `↓`), their multiplication signs, and
their spacing inside units. The numbers and symbols must stay byte-identical; only
the words around them change.

## Path B — whole-file translation (prose-dominant files)

Use for the per-disease `README.md`. These are dense scientific prose in which
English sentences do not line up with Korean sentences, so translate the file as a
document and write it to `translations/en/<same path>`.

Preserve the document as a document:

- keep the heading structure, table structure, and list structure
- **keep a code fence's *syntax* — its commands, flags, file names, and
  punctuation — but translate the Korean *comments* inside it, exactly as Path A
  does.** "Keep code fences" means do not restructure or delete them, not "leave
  every character inside untouched." A run of 30+ READMEs shared one boilerplate
  usage block whose two comment lines — `# Shiny 대시보드:` and
  `# 기계론적 지도 렌더링` — were skipped for exactly this reason: the fence read as
  "code, don't touch," and the comment inside it went untranslated while
  everything else in the file was translated correctly. `grep -cP
  '[\x{AC00}-\x{D7A3}]' translations/en/<path>` must be 0 for the *whole* file,
  fences included — that is what would have caught this immediately, and is why
  it is listed as a required check and not a suggestion.
- keep every file name, path, link, number, unit, gene name, drug name, and
  parameter name exactly as written
- keep relative links working. **Count the levels — this is the easiest thing to get
  wrong.** A relative link resolves from the directory *containing* the file, so the
  number of `../` needed is the file's depth below the repository root:

  | translation | its directory | to reach the repo root |
  |---|---|---|
  | `translations/en/README.md` | `translations/en/` | `../../` |
  | `translations/en/<disease>/README.md` | `translations/en/<disease>/` | `../../../` |

  So from a per-disease README, a link to a sibling file in the original directory is
  `../../../<disease>/<file>` — three levels, not two. Link to the translation
  instead if that file has itself been translated. Then check it mechanically:

  ```bash
  python3 translations/tools/check_links.py translations/en/<disease>/README.md
  ```
- keep the emphasis the author used — these READMEs deliberately bold the claim
  each model is built around; do not flatten that

## Register

This is scientific writing for a pharmacometrics audience. Match it:

- British spelling, to match the terminology already used across the library
  (`haemoglobin`, `oedema`, `fibre`, `modelling`, `-isation`)
- standard clinical and pharmacometric terms, not literal glosses: 기계론적 지도 is
  a "mechanistic map", 구획 is a "compartment", 시나리오 is a "scenario",
  고정점 is a "fixed point", 되먹임 is "feedback", 쌍안정 is "bistable",
  용량-반응 is "dose-response"
- keep the author's voice. Many of these files argue a specific structural claim
  ("writing X as Y rather than Z derives W"). Translate the argument, do not
  summarise it or soften it into a description.
- do not add content, do not add hedging the author did not write, and do not
  correct the author's science even where it looks arguable — this is a
  translation

## A hand-fix to a non-Korean line is not durable

`apply` reconstructs every non-Korean line from the pristine original, on every
run. So if a line has **zero Korean** — a bare sibling link the original correctly
writes as `[file](file.R)`, which only needs the `../../../` prefix once it moves
three levels down — hand-editing it directly in `translations/en/` does not survive
a second `apply` on that same file. It looks fixed, `check_links` goes green, and
then a duplicate or later batch touching the same file's Korean content calls
`apply` again and silently reverts it, because from `apply`'s point of view the
hand-edit was drift from the original that its job is to remove.

This happened for real: an orchestrator hand-fixed three bare links, committed, and
a concurrently-running batch — assigned the same file by an earlier, now-stale
listing — re-ran `apply` from its own TSV and reverted all three, believing it was
restoring the correct state.

**So:** never hand-edit a bare link in a translated file while any other batch
might still touch that file's Korean content. Fix it with
`python3 translations/tools/relink.py --write` instead — same fragility in
principle, but it is the one tool everyone runs, so run it **last**, after every
translation batch in the wave has finished and committed, not mid-wave. If you
must hand-fix immediately, re-run `check_links.py` again after all other agents in
the wave have reported, since your fix may already be gone.

## When a gate fails, suspect the gate

The checks here are ordinary scripts and they have been wrong more than once — the
number check has already had three bugs, each of which reported honest translations
as damaged. So:

> **Never reword a translation to make a gate pass.** If a gate flags something you
> believe is correct, say so in your report and leave the translation as the source
> demands. Rewording to satisfy tooling silently degrades the text and hides the bug
> from whoever has to fix it.

A gate failure is a question, not a verdict. Read what it actually flagged, decide
whether the file or the checker is wrong, and report which.

## Coverage

```bash
python3 translations/tools/translation_status.py            # summary table
python3 translations/tools/translation_status.py --list todo
python3 translations/tools/translation_status.py --check     # residual Hangul
```

After translating any `*.dot`, render its images — the upstream `.svg`/`.png` still
carry Korean labels:

```bash
python3 translations/tools/render_maps.py
```

A translated source file must still parse, and must still have the same shape. Check
both, not just the first:

```bash
# R: parses, and yields the same number of top-level expressions as the original
Rscript -e 'cat(length(parse("translations/en/<path>")), length(parse("<path>")), "\n")'

# Python
python3 -m py_compile translations/en/<path>

# Graphviz
dot -Tsvg translations/en/<path> -o /dev/null
```

**Not every `.R` file here is an R script.** Two different things carry that
extension:

| Shape | First non-comment line | Gate |
|---|---|---|
| R script with the model in a quoted string | `library(mrgsolve)`, `mod <- mread(...)` | `parse()`, expression counts must match |
| pure mrgsolve DSL | `$PROB`, `$PARAM` at top level | `parse()` **cannot** work — `$` is not valid at R top level |

For a DSL file, `parse()` failing is correct behaviour and not a defect in either
the original or the translation. `lineio verify` is the gate there, and it is a
strong one: it proves every line without Korean is byte-identical, and all the
`$BLOCK` headers, parameter symbols, and values are ASCII, so they cannot have
moved.

So when `parse()` fails, find out which case you are in before reporting anything.
**The block marker has three spellings in this library** — `$PROB`, `[PROB]`, and
`[ PROB ]` with inner spaces — so test for all of them:

```bash
grep -nE '^[ \t]*(\$|\[[ \t]*)(PROB|PARAM|CMT|MAIN|ODE|TABLE)' <path> | head -1
```

A hit means DSL, and `parse()` failing is correct. Checking only `^\$` reports
dozens of healthy bracket-form files as broken.

Or just run the audit, which classifies every `.R` file in the repository:

```bash
python3 translations/tools/check_r_parses.py
```

If it is a script and the original fails too, that is a genuine upstream bug —
report it, and see `../UPSTREAM_ISSUES.md`, where two Shiny apps that cannot start
are already recorded.

### A translated label can silently change a plot's colours

`ggplot2` orders factor levels alphabetically, and `scale_*_manual(values = c(...))`
with an **unnamed** vector assigns colours in that order. So renaming
`무치료`/`치료` to `Untreated`/`Treated` can swap two colours, because the Korean and
English sort differently — the legend stays correct, but every line changes colour
and any claim in the prose about "the red curve" stops being true.

When you translate labels feeding a manual scale:

1. Check whether the `values =` vector is **named** (`c("Treated" = "#a00", …)`). If
   it is, translate the names to match and order does not matter.
2. If it is unnamed, choose English wording that preserves the sort positions where
   you reasonably can.
3. Where no natural term preserves the order, keep the correct clinical term and
   **say so in your report**. Contorting terminology to protect a colour is the
   wrong trade, but it must not be silent.

The root `README.md` is a third case: it is generated. See
[`tools/build_root_readme.py`](tools/build_root_readme.py).
