#!/usr/bin/env python3
"""Generate ``README_en.md`` (repo root) from the upstream root ``README.md``.

Most of a gallery row is mechanical: the directory paths, the ``<sub>`` subtitle
(already English upstream), and the four link labels. Only three things need
translating, and each has its own source of truth:

* the **category** column           -> ``translations/data/categories.tsv``
* the **English disease name**       -> taken from the row's own ``<sub>`` subtitle
* the **one-sentence summary**       -> ``translations/data/root_readme_summaries.tsv``

The surrounding prose is hand-translated once and lives in
``translations/data/root_readme_head.md`` and ``root_readme_tail.md``.

Doing it this way means upstream's daily new row costs exactly one new line in
the summaries TSV instead of a re-translation of the whole file.

Usage
-----
    python3 translations/tools/build_root_readme.py            # write the file
    python3 translations/tools/build_root_readme.py --check     # report gaps only
    python3 translations/tools/build_root_readme.py --todo      # dump untranslated
                                                                # rows as TSV stubs
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "translations" / "data"
OUT = REPO / "README_en.md"

ROW = re.compile(
    r"^\|\s*(\d+)\s*\|\s*([^|]+?)\s*\|\s*(.*?)\s*\|\s*(<a href.*?</a>)\s*\|\s*(.*?)\s*\|\s*$"
)
MODEL_CELL = re.compile(r"^\[\*\*(.+?)\*\*(?:<br><sub>(.*?)</sub>)?\]\((.+?)/\)$")
HANGUL = re.compile(r"[가-힣ᄀ-ᇿ㄰-㆏]")

# The four trailing links are identical in every row apart from the paths, so
# their Korean labels are mapped rather than translated per row.
LINK_LABELS = {"지도": "Map", "문헌": "Refs"}


def load_tsv(name: str) -> dict[str, str]:
    path = DATA / name
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        key, _, val = line.partition("\t")
        out[key.strip()] = val.strip()
    return out


def parse_rows(readme: str) -> list[dict]:
    rows = []
    for line in readme.splitlines():
        m = ROW.match(line)
        if not m:
            continue
        idx, cat, model, preview, summary = m.groups()
        mc = MODEL_CELL.match(model)
        if not mc:
            raise SystemExit(f"row {idx}: cannot parse model cell: {model[:120]}")
        ko_name, sub, directory = mc.groups()
        rows.append(
            {
                "idx": int(idx),
                "cat_ko": cat,
                "ko_name": ko_name,
                "sub": sub or "",
                "dir": directory,
                "preview": preview,
                "summary_ko": summary,
            }
        )
    return rows


def english_name(row: dict) -> str:
    """The English disease name, taken from the row's own subtitle.

    Upstream writes the subtitle as ``English Name . ABBR`` (with additional
    mechanism text on some rows), so the first bullet-separated segment is the
    name. Fall back to the directory name if the subtitle is unusable.
    """
    sub = row["sub"]
    if sub:
        first = sub.split("·")[0].strip()
        if first and not HANGUL.search(first):
            return first
    return row["dir"].replace("-", " ").title()


def trimmed_subtitle(row: dict, english: str) -> str:
    """The subtitle with the disease name removed when it merely repeats it.

    Upstream writes ``<sub>English Name . ABBR</sub>`` next to a Korean disease
    name, so once the name cell is English the first segment is redundant.
    """
    sub = row["sub"]
    if not sub:
        return sub
    parts = [p.strip() for p in sub.split("·")]
    if parts and parts[0].lower() == english.lower():
        parts = parts[1:]
    return " · ".join(parts) or english


def resolve(rel: str) -> str:
    """Resolve one repo-relative path as seen from the repo root.

    A path resolves to its ``<name>_en.<ext>`` sibling once that translation
    exists, so the English gallery wires itself to English content as the
    translation progresses, and falls back to the upstream original until then.
    Either way the path needs no ``../`` prefix: ``README_en.md`` lives at the
    repo root, the same level as every disease directory.
    """
    p = Path(rel)
    en_rel = (p.parent / f"{p.stem}_en{p.suffix}").as_posix()
    if (REPO / en_rel).exists():
        return en_rel
    return rel


def rewrite_paths(text: str) -> str:
    """Rewrite every relative markdown link and HTML attribute in ``text``."""
    text = re.sub(
        r'\]\((?!https?:|#|\.\./)([^)]+)\)',
        lambda m: "](" + resolve(m.group(1)) + ")",
        text,
    )
    text = re.sub(
        r'((?:href|src)=")(?!https?:|\.\./)([^"]+)"',
        lambda m: m.group(1) + resolve(m.group(2)) + '"',
        text,
    )
    return text


def build_summary_cell(row: dict, english: str) -> str:
    """English summary plus the row's four links, with paths rewritten."""
    ko = row["summary_ko"]
    _, _, links = ko.partition("<br>")
    if not links:
        links = ""
    for ko_label, en_label in LINK_LABELS.items():
        links = links.replace(f"] {ko_label}]", f"] {en_label}]")
        links = links.replace(f"{ko_label}](", f"{en_label}](")
    links = rewrite_paths(links)
    return f"{english}<br>{links}" if links else english


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--todo", action="store_true")
    args = ap.parse_args()

    readme = (REPO / "README.md").read_text(encoding="utf-8")
    rows = parse_rows(readme)
    cats = load_tsv("categories.tsv")
    summaries = load_tsv("root_readme_summaries.tsv")

    missing_cat = sorted({r["cat_ko"] for r in rows if r["cat_ko"] not in cats})
    missing_sum = [r for r in rows if r["dir"] not in summaries]

    if args.todo:
        for r in missing_sum:
            print(f"{r['dir']}\t{r['summary_ko'].partition('<br>')[0]}")
        return 0

    if missing_cat:
        print(f"categories.tsv is missing {len(missing_cat)}: {missing_cat}",
              file=sys.stderr)
    print(f"rows: {len(rows)}  summaries translated: "
          f"{len(rows) - len(missing_sum)}/{len(rows)}"
          f"  ({100 * (len(rows) - len(missing_sum)) / len(rows):.0f}%)")
    if missing_sum:
        print(f"untranslated summaries: {len(missing_sum)}"
              f" (first: {', '.join(r['dir'] for r in missing_sum[:5])})")

    if args.check:
        return 1 if (missing_cat or missing_sum) else 0

    counts = Counter(cats.get(r["cat_ko"], r["cat_ko"]) for r in rows)
    cat_line = " · ".join(f"{k} {v}" for k, v in counts.most_common())

    head = (DATA / "root_readme_head.md").read_text(encoding="utf-8")
    head = head.replace("{{N_MODELS}}", str(len(rows)))
    head = head.replace("{{CATEGORY_COUNTS}}", cat_line)
    tail = (DATA / "root_readme_tail.md").read_text(encoding="utf-8")

    out = [head, "", "| # | Category | Model | Preview | Summary & links |",
           "|---|------|------|----------|--------------|"]
    for r in rows:
        english = english_name(r)
        sub = trimmed_subtitle(r, english)
        cell_model = f"[**{english}**<br><sub>{sub}</sub>]({r['dir']}/)"
        summary_en = summaries.get(r["dir"], r["summary_ko"].partition("<br>")[0])
        cell_summary = build_summary_cell(r, summary_en)
        out.append(
            f"| {r['idx']} | {cats.get(r['cat_ko'], r['cat_ko'])} | {cell_model} "
            f"| {rewrite_paths(r['preview'])} | {cell_summary} |"
        )
    out += ["", tail]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(out), encoding="utf-8")
    print(f"wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
