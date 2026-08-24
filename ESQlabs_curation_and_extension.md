# ESQlabs curation & extension of this fork

> This file is ESQlabs' own tracking document, not a translation of anything
> upstream. It records what we changed about the *structure* of this fork (as
> opposed to the *content* translation work, which is process-documented in
> [`translations/WORKFLOW.md`](translations/WORKFLOW.md) and
> [`translations/UPSTREAM_ISSUES.md`](translations/UPSTREAM_ISSUES.md)) and why,
> so a future session — human or Claude — can understand a structural decision
> without having to reconstruct the reasoning from git log archaeology.

## Status: migration complete

Translation content itself reached 100% coverage (894/894 tracked files) before
this migration started. What follows is a *layout* change: moving the finished
translations from a separate mirrored tree into same-directory siblings of the
upstream originals.

## Background: two layouts considered

**Layout A (original, 2026-08 through the 100%-completion milestone):** every
translation lives under `translations/en/<original path>`, mirroring the
upstream tree exactly. A translation of `heat-stroke/README.md` is
`translations/en/heat-stroke/README.md`.

**Layout B (this migration):** every translation lives directly next to its
original, named `<original name>_en.<ext>`. The translation of
`heat-stroke/README.md` is `heat-stroke/README_en.md`.

Layout A was chosen first because it gives the strongest possible guarantee
against ever conflicting with upstream: upstream does not know
`translations/` exists, full stop. Layout B trades some of that guarantee for
much simpler in-document links (same-directory bare filenames instead of
`../../../<dir>/<file>` back-references) and — the actual reason for the
switch — files that read naturally when browsing the repository directly on
GitHub, next to the file they translate, rather than in a parallel tree a
reader has to know to look for.

## Risks of Layout B, and how each is handled

1. **Filename collision with something upstream adds later.** If upstream
   ever created a file literally named e.g. `README_en.md` in a disease
   directory, our file of the same name would collide on the next `git pull`.
   Assessed as low-probability (upstream's own convention is Korean-first, and
   an `_en` suffix is not part of its naming scheme per `CLAUDE.md`'s File
   Naming table) but non-zero. **Mitigation:** none structural — this risk is
   accepted as the price of Layout B. If it ever happens, `git status` will
   surface the conflicting add immediately (a new untracked file with a name
   we already use), not a silent overwrite.
2. **`git status`/diff noise in upstream directories.** Every disease
   directory now visibly contains both the Korean original and our `_en`
   sibling. This was the explicit reason Layout A's rule was written
   (`translations/README.md`'s old text: "Do not create sibling files... make
   `git status` noisy after every pull"). **Accepted tradeoff**, not mitigated
   — this is the readability gain the migration is for.
3. **Upstream's own directory-scanning tooling.** `scripts/fix_readme_table.py`
   *does* glob every disease directory (`find_disk_dirs()` and `pick()`, both
   using patterns like `*_qsp*.dot` and `*_mrgsolve_model.R`), and after this
   migration such a glob matches both a disease's original file and its `_en`
   sibling. **Verified rather than assumed:** `pick()` takes
   `sorted(glob.glob(...))[0]` — the *first* match after an alphabetical
   sort — and `.` (0x2E) sorts before `_` (0x5F) in ASCII, so
   `abbr_qsp_model.png` always sorts ahead of `abbr_qsp_model_en.png` and
   `pick()` always resolves to the original regardless of which one glob finds
   first. Confirmed empirically: `python3 scripts/fix_readme_table.py --check`
   against the real, unmodified upstream `README.md` still reports `PASS`
   after the migration, exactly as before it. The **translated** copy,
   `scripts/fix_readme_table_en.py`, targets `README_en.md` at the repo root
   instead — not by design, but by a happy accident of the same-directory
   bare-substitution mechanism (see the migration script's docstring):
   `scripts/README.md` is a real, different file in the same directory as
   `fix_readme_table.py`, so the substitution correctly rewrote *that*
   sibling reference to `README_en.md` wherever the bare string `"README.md"`
   appeared in the translated script's text — including inside the line that
   actually meant "the root gallery two directories up", which happens to use
   the same literal string. Verified working with
   `python3 scripts/fix_readme_table_en.py --check`, which passed against the
   freshly generated `README_en.md`.
4. **Same-directory literal filename references inside translated scripts.**
   This is the one *content* risk the migration required fixing, not just
   tooling: `<abbr>_shiny_app.R` files `source()` their sibling
   `<abbr>_mrgsolve_model.R` by literal string (directly, or indirectly via a
   `MODEL_FILE <- "..."` variable). In Layout A this worked automatically
   because both translated files sat in the same mirrored subdirectory. In
   Layout B, the literal string must become `<abbr>_mrgsolve_model_en.R`,
   or the translated Shiny app would silently source the **Korean original**
   model file lying in the same upstream directory instead of its own
   translated sibling. 64 Shiny app files had this exact pattern at migration
   time (see Process log for the exact count and the verification method).
   Plain-text (non-hyperlink) filename mentions — e.g. a `## Files` directory
   listing in a README — have the identical risk and were swept the same way.
5. **In-document relative links change depth.** A file two levels under
   `translations/en/` needed `../../../` to reach the repo root; the same file
   one level under its own upstream directory needs only `../`. Every markdown
   link and HTML `href`/`src` attribute in every translated file needed its
   relative path recomputed, not just its filename adjusted. Handled
   mechanically (see Process log) by resolving each link against its *old*
   location, remapping the *target* through the same old→new table if the
   target was itself a translated file, then re-emitting the link relative to
   the file's *new* location. A file's own "this is a translation of `X`"
   attribution backlink is explicitly excluded from the target remap (it must
   keep pointing at the real Korean original, not at itself) — the same
   exclusion `translations/tools/relink.py` already applied under Layout A.
6. **Generated files.** `translations/en/README.md` (the translated root
   gallery) is generated by `translations/tools/build_root_readme.py`, not
   hand-translated. After migration it lives at `README_en.md` (repo root).
   The generator was retargeted to write there and to link disease-directory
   images at their new co-located paths.

## What did NOT need to change

* No translated **content** (prose, comments, code) required re-translation —
  the migration is a rename + relink + retarget operation, not a re-review of
  translation quality.
* Upstream files themselves are still never touched, still never will be —
  the fork rule from `translations/README.md` stands. Layout B still only
  ever *adds* new files; it changes *where* they are added, not the "never
  edit upstream" guarantee.

## Process log

Filled in as each step is executed and verified — this is the actual record
of what ran, not a plan written in advance and left stale.

- [x] Migration script written (`translations/tools/migrate_to_sibling_en.py`)
      and dry-run on two test directories, then the full tree, before the real
      run. Two real bugs were caught and fixed by the dry-run before anything
      was written to disk:
      1. The markdown-link regex, copied from `relink.py`'s established
         pattern, did not handle the "linked image" idiom
         `[![alt](imgtarget)](linktarget)`: it matched only the inner image
         link and left the outer link's target — which also needed its
         relative depth recomputed — untouched. Fixed by combining both into
         one pattern so a single pass consumes the whole nested construct.
      2. The bare-filename substitution (for `source("x.R")`-style same-
         directory references) did an unanchored substring replace, which
         corrupted an unrelated, already-correct link:
         `translations/README.md` was rewritten to `translations/README_en.md`
         because it contains `README.md` as a trailing substring, even though
         it points at a different file entirely (the un-migrated
         `translations/README.md` meta-doc, not the root gallery). Fixed with
         a path-boundary-aware regex (`(?<![\w/.-])name(?![\w-])`).
- [x] Full-tree migration executed: all 1098 files moved out of
      `translations/en/` to `<name>_en.<ext>` siblings; the now-empty
      `translations/en/` tree was removed. `git status` afterward showed
      exactly 1098 deletions and 1098 additions, matching expectation exactly.
- [x] Same-directory literal-filename references fixed as part of the same
      migration pass (the bare-substitution mechanism above) — this covered
      both `source()`/`MODEL_FILE`-indirection references in Shiny apps and
      plain-text `## Files` directory-listing mentions in READMEs, verified
      by spot-checking `mah_shiny_app.R`'s `.model <- file.path(.here, ...)`
      indirection and several READMEs' file tables before and after.
- [x] In-document link depth recomputed across all migrated files. One real
      gap surfaced only after the full run: a "linked image"
      `[![alt](https://img.shields.io/badge/...)](../local/target)` badge —
      external inner image, local outer link — didn't match either regex
      alternative (the image alternative excluded `https:` targets at the
      *regex* stage, so the whole nested construct fell through unmatched).
      4 files had this shape; their outer badge links were left pointing at
      the pre-migration depth and `check_links.py` correctly caught all of
      them. Fixed the script (scheme filtering now happens per-target inside
      the substitution function, after a broad match, not inside the regex)
      and hand-fixed the 4 already-migrated files directly, since
      `translations/en/` no longer existed to re-run the script against.
      Re-ran `check_links.py`: clean across all 524 files afterward.
- [x] `translations/tools/*.py` retargeted to the new layout: `lineio.py`
      (added a shared `en_path()` helper other tools now import),
      `translation_status.py`, `check_tokens.py`, `check_links.py`,
      `check_tables.py`, `check_r_parses.py`, `render_maps.py`, `relink.py`
      (simplified: no more depth recomputation needed, just a same-directory
      `_en` check), `batches.py` (no change needed — it only parses
      `translation_status.py`'s output, which is upstream-relative paths on
      both sides of this migration).
- [x] `build_root_readme.py` retargeted to write `README_en.md` at the repo
      root and link disease-directory images/files at their new co-located
      paths (no more `../../` prefix — the gallery now sits at the same
      level as every disease directory). `root_readme_head.md` and
      `root_readme_tail.md` had two more hardcoded `../../` links fixed the
      same way. Regenerated and spot-checked both a partially-translated row
      (`abdominal-aortic-aneurysm`, README only) and a fully-translated one
      (`achalasia`, every file) — each correctly links whichever of the
      original or the `_en` sibling actually exists.
- [x] Full gate suite re-run and green: `check_links.py` (524 files, OK after
      the badge-link fix above), `check_tokens.py` (894 files, OK),
      `check_tables.py --translations` (4 pre-existing ragged tables, same
      ones logged before the migration — no regression), `relink.py` (0 stale
      links), `translation_status.py` (894/894, 100%, unchanged from before
      the migration), `check_r_parses.py --include-translations` (840
      originals + 230 new `_en.R` files; the only 3 failures are the 2
      already-logged broken originals plus their 1 translated sibling that
      inherits the same defect verbatim, as intended — not a new fault).
      `render_maps.py` surfaced one new (pre-existing, not migration-caused)
      upstream defect: `retinitis-pigmentosa/rp_qsp_model.dot` crashes this
      machine's Graphviz exactly like the already-logged
      `visceral-leishmaniasis` case (`UPSTREAM_ISSUES.md` #20) — confirmed by
      reproducing the identical crash against the untouched original `.dot`
      directly. Logged as `UPSTREAM_ISSUES.md` #26.
- [x] `translations/README.md` and `WORKFLOW.md` updated to describe Layout B:
      the mirrored-tree layout, the three-level link-depth table, and the
      old `../../../` examples are gone, replaced with same-directory
      `_en`-suffix guidance and a note that a translated sibling's own file
      name (in `source()` calls, `MODEL_FILE`-style variables, and plain-text
      file listings) must also gain the `_en` suffix once translated.
- [x] Final commit
