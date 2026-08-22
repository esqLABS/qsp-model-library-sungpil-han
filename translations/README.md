# `translations/` — English translations for this fork

This repository is a **fork** of an upstream QSP model library whose text is
written mostly in Korean. Upstream is actively developed: a new disease model is
committed almost every day, and existing files (notably the root `README.md`
table) are rewritten in place.

## The rule

> **Never edit an upstream-tracked file. Every addition of ours lives in its own
> new file, inside this directory.**

An in-place edit of a Korean original would become a merge conflict on every
`git pull` from upstream, forever. A new file in a directory upstream does not
know about can never conflict.

Concretely:

* Translations go to `translations/en/<original path>`, mirroring the upstream
  layout exactly. The English version of `heat-stroke/README.md` is
  `translations/en/heat-stroke/README.md`.
* Do **not** create sibling files such as `heat-stroke/README.en.md`. Those sit
  inside upstream directories and make `git status` noisy after every pull.
* Translated Graphviz maps (`*_qsp_model.dot`) get their own re-rendered
  `.svg`/`.png` next to the translated `.dot` inside `translations/en/`.
* Tooling of ours goes in `translations/tools/`.

## Layout

```
translations/
  README.md                     this file
  tools/
    translation_status.py       coverage report: what still needs translating
  en/
    README.md                   translation of the root README.md
    CLAUDE.md                   translation of the root CLAUDE.md
    scripts/README.md
    <disease>/README.md
    <disease>/<abbr>_references.md
    <disease>/<abbr>_qsp_model.dot   (+ .svg / .png re-rendered from it)
    <disease>/<abbr>_shiny_app.R
    ...
```

Only files that actually contain Korean text are mirrored. A file that is
already fully English upstream is intentionally absent here — there is nothing
to translate, and mirroring it would create a second copy to keep in sync.

## Resuming the work

The translation is a long job — 893 files at the start — and it is designed to be
picked up mid-way by someone with no memory of the previous session. Everything
needed is on disk:

```bash
python3 translations/tools/translation_status.py     # what is done, per role
python3 translations/tools/batches.py                # what is left, as batches
python3 translations/tools/batches.py --role readme --index 1   # the next batch
```

Then follow [`WORKFLOW.md`](WORKFLOW.md) for the batch, and run the three gates
before committing:

```bash
python3 translations/tools/lineio.py verify <path>   # Path A files
python3 translations/tools/check_tokens.py           # add --loose for Path B files
python3 translations/tools/check_links.py
python3 translations/tools/render_maps.py            # after any .dot
```

Defects found in the originals go in [`UPSTREAM_ISSUES.md`](UPSTREAM_ISSUES.md), not
into the originals. Reviewed token differences go in
[`data/token_exceptions.tsv`](data/token_exceptions.tsv).

## Checking coverage

```bash
python3 translations/tools/translation_status.py
```

It reports, for every upstream file containing Hangul, whether a translation
exists under `translations/en/`, whether the translation is stale (the original
changed after the translation was last written), and whether any Hangul is left
in the translation itself.
