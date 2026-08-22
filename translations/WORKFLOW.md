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

## Path B — whole-file translation (prose-dominant files)

Use for the per-disease `README.md`. These are dense scientific prose in which
English sentences do not line up with Korean sentences, so translate the file as a
document and write it to `translations/en/<same path>`.

Preserve the document as a document:

- keep the heading structure, table structure, list structure, and code fences
- keep every file name, path, link, number, unit, gene name, drug name, and
  parameter name exactly as written
- keep relative links working from `translations/en/<dir>/`: a link to a sibling
  file in the same disease directory becomes `../../<dir>/<file>`, unless that file
  has itself been translated, in which case link to the translation
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

The root `README.md` is a third case: it is generated. See
[`tools/build_root_readme.py`](tools/build_root_readme.py).
