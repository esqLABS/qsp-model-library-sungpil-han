# scripts/

> English translation of the upstream [`scripts/README.md`](README.md),
> maintained in this fork. The original is already in English apart from one
> cross-reference; see [`translations/README.md`](../translations/README.md) for why the
> translation lives here rather than as an edit to the original.

Maintenance tooling for the root `README.md` model gallery table.

## `fix_readme_table_en.py`

Normalizes the gallery table so it can never again break the way it did on
2026-07-02 (a stray unmatched `**` in one row's summary corrupted GFM's
rendering for every row after it; the trailing Disclaimer/References/License
section had also been lost).

```bash
# After adding a new model's row to README.md, before committing:
python3 scripts/fix_readme_table.py

# Lint-only: check README.md for the known failure modes, change nothing,
# exit 1 if anything is wrong:
python3 scripts/fix_readme_table.py --check
```

See the module docstring in the script for full details, and
`CLAUDE.md` → "README Update" for the row-writing rules this script
enforces (one-sentence summaries, no bold/citations in cells, capped
sub-titles, etc.). That section is translated in
[`translations/en/CLAUDE.md`](../CLAUDE_en.md).
