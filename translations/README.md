# `translations/` — tooling and process docs for this fork's English translations

This repository is a **fork** of an upstream QSP model library whose text is
written mostly in Korean. Upstream is actively developed: a new disease model is
committed almost every day, and existing files (notably the root `README.md`
table) are rewritten in place.

## The rule

> **Never edit an upstream-tracked file. Every addition of ours is a new file.**

An in-place edit of a Korean original would become a merge conflict on every
`git pull` from upstream, forever. A new file upstream does not know about can
never conflict.

Concretely:

* A translation of `<original path>` lives at `<original path>` with `_en`
  inserted before the extension: the English version of
  `heat-stroke/README.md` is `heat-stroke/README_en.md`, next to the original,
  not in a separate mirrored tree. **This changed** partway through the
  project — see [`../ESQlabs_curation_and_extension.md`](../ESQlabs_curation_and_extension.md)
  for why, and for the risks that layout accepts in exchange for simpler links.
* Translated Graphviz maps (`*_qsp_model.dot`) get their own re-rendered
  `.svg`/`.png` — `*_qsp_model_en.svg`/`.png` — next to the translated `.dot`.
* Tooling of ours goes in `translations/tools/`; reviewed data goes in
  `translations/data/`. Those two directories, plus this file, `WORKFLOW.md`,
  and `UPSTREAM_ISSUES.md`, are the only things that still live under
  `translations/` — the translations themselves do not.

Only files that actually contain Korean text get a translation. A file that is
already fully English upstream has no `_en` sibling — there is nothing to
translate, and creating one would just be a second copy to keep in sync.

## Resuming the work

Coverage is currently 100% (894/894 tracked files), but upstream adds new
files continuously, so this will drift. Everything needed to pick it back up
is on disk:

```bash
python3 translations/tools/translation_status.py     # what is done, per role
python3 translations/tools/batches.py                # what is left, as batches
python3 translations/tools/batches.py --role readme --index 1   # the next batch
```

Then follow [`WORKFLOW.md`](WORKFLOW.md) for the batch, and run the gates
before committing:

```bash
python3 translations/tools/lineio.py verify <path>   # Path A files
python3 translations/tools/check_tokens.py           # mode chosen per file automatically
python3 translations/tools/check_links.py
python3 translations/tools/check_tables.py --translations
python3 translations/tools/render_maps.py            # after any .dot
python3 translations/tools/relink.py --write          # after any batch
```

Defects found in the originals go in [`UPSTREAM_ISSUES.md`](UPSTREAM_ISSUES.md), not
into the originals. Reviewed token differences go in
[`data/token_exceptions.tsv`](data/token_exceptions.tsv).

## A note on the translated scripts

`scripts/fix_readme_table_en.py` and the various `mkrefs_en.py` files are
translations, provided so an English reader can follow what the tooling does.
`fix_readme_table_en.py` is retargeted to operate on `README_en.md` (the
generated English gallery) rather than the upstream `README.md` it was copied
from — regenerate `README_en.md` with
[`tools/build_root_readme.py`](tools/build_root_readme.py) instead of running
a normaliser over it. The operative copy of `fix_readme_table.py` is
upstream's, and it is upstream's that `CLAUDE.md` tells a session to run.

## Checking coverage

```bash
python3 translations/tools/translation_status.py
```

It reports, for every upstream file containing Hangul, whether an `_en`
translation exists beside it, whether the translation is stale (the original
changed after the translation was last written), and whether any Hangul is left
in the translation itself.
