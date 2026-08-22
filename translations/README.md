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

## Checking coverage

```bash
python3 translations/tools/translation_status.py
```

It reports, for every upstream file containing Hangul, whether a translation
exists under `translations/en/`, whether the translation is stale (the original
changed after the translation was last written), and whether any Hangul is left
in the translation itself.
