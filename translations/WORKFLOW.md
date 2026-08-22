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

### Two situations that come up constantly — handle them the same way every time

**A heading that already carries an English gloss.** Many headings are written
`## 3. 유전적 소인 (Genetic susceptibility)`. Translate the Korean and **drop the
now-duplicated bracketed gloss** — never emit
`## 3. Genetic susceptibility (Genetic susceptibility)`. Where the author's gloss is
narrower than the Korean it glosses, keep the fuller translation rather than the
gloss: `## 5. 담즙산 합성·되먹임·수송 (Bile acid feedback)` becomes
`## 5. Bile acid synthesis · feedback · transport`.

**A Graphviz label that pairs a Korean title with an English subtitle.** Cluster and
node labels are often written
`<b>6. 고환 축</b><br/><i>Testicular axis</i>` — all on one line. Because it is one
line you may collapse it, and you should: translating both halves renders visible
duplication in the image. Keep the styling, drop the repetition. The same applies to
`ORGANISING THESIS (조직 원리)`, where the parenthetical merely glosses the English
beside it — drop the parenthetical.

**A Korean passage the author already paraphrased in English underneath.** Some
files carry a Korean note followed by the author's own shorter English version of
it. Path A cannot delete lines, so the translation will hold both, the full one and
the author's abridged one. **Leave both.** Do not delete, merge, or trim either —
deduplicating is an editorial change, not a translation, and it would also break the
line alignment that `verify` depends on. Mention it in your report and move on.

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
- list names consumed as legend text or as a table's first column are labels; list
  names matched against a literal elsewhere are identifiers.
- mrgsolve output column names (`PNVOLWK`, `OCCGC`, …) and `$PARAM`/`$CMT` symbols
  are always identifiers.

### A label whose second line is an English gloss

Plot labels are often `경구 섭취\noral` — Korean, then a lower-case English gloss.
Collapse those to the single English label; do not emit "Oral intake / oral". Keep
the `\n` where the second line carries something else (a variable name, a unit, an
abbreviation). Either way the number of elements in the vector must not change.

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

- keep the heading structure, table structure, list structure, and code fences
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

The root `README.md` is a third case: it is generated. See
[`tools/build_root_readme.py`](tools/build_root_readme.py).
