#!/usr/bin/env python3
"""Coverage report for the English translations living next to their originals.

For every git-tracked file that contains Hangul, report whether a
``<name>_en.<ext>`` translation exists beside it, whether that translation
still contains Hangul, and whether the original has changed since the
translation was committed (stale).

Usage
-----
    python3 translations/tools/translation_status.py            # summary
    python3 translations/tools/translation_status.py --list todo # paths still to do
    python3 translations/tools/translation_status.py --list stale
    python3 translations/tools/translation_status.py --list residual
    python3 translations/tools/translation_status.py --check     # exit 1 if any
                                                                # translation has
                                                                # residual Hangul
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Hangul syllables, Jamo, and compatibility Jamo -- and the \uXXXX escape form,
# which some files use instead. Reusing lineio's matcher rather than a second
# regex here: when the two disagreed, this tool reported files as fully
# translated while nine lines of one Shiny app still rendered Korean.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lineio import HANGUL, en_path  # noqa: E402

# Files we never try to translate.
BINARY_SUFFIXES = {
    ".png", ".svg", ".jpg", ".jpeg", ".webp", ".woff2", ".npz",
    ".docx", ".pdf", ".ico", ".gif",
}


def tracked_files() -> list[str]:
    out = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8", check=True,
    ).stdout
    return [p for p in out.split("\0") if p]


def hangul_lines(path: Path) -> int:
    """Number of lines containing Hangul; 0 if unreadable as text."""
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return 0
    return sum(1 for line in text.splitlines() if HANGUL.search(line))


def last_commit_epoch(rel: str) -> int:
    """Author time of the last commit touching ``rel``; 0 if never committed."""
    out = subprocess.run(
        ["git", "log", "-1", "--format=%at", "--", rel],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8",
    ).stdout.strip()
    return int(out) if out.isdigit() else 0


def role(rel: str) -> str:
    if rel.endswith("_references.md"):
        return "references.md"
    if rel.endswith("README.md"):
        return "README.md"
    if rel.endswith(".md"):
        return "other .md"
    if "shiny" in rel.rsplit("/", 1)[-1] and rel.endswith(".R"):
        return "Shiny app .R"
    if rel.endswith(".R"):
        return "other .R"
    if rel.endswith(".dot"):
        return ".dot map"
    if rel.endswith(".py"):
        return ".py"
    return "other"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", choices=["todo", "done", "stale", "residual"])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any translation has residual Hangul")
    args = ap.parse_args()

    todo, done, stale, residual = [], [], [], []
    by_role: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0, 0])  # files,lines,doneF,doneL

    for rel in tracked_files():
        if rel.startswith("translations/"):
            continue
        src = REPO / rel
        if src.suffix.lower() in BINARY_SUFFIXES or not src.is_file():
            continue
        n = hangul_lines(src)
        if not n:
            continue

        r = role(rel)
        by_role[r][0] += 1
        by_role[r][1] += n

        dst = en_path(rel)
        if not dst.exists():
            todo.append((n, rel))
            continue

        by_role[r][2] += 1
        by_role[r][3] += n
        done.append((n, rel))

        left = hangul_lines(dst)
        if left:
            residual.append((left, rel))
        # A translation with no commit of its own is new, not stale: it was just
        # written and is about to be committed. Only compare commit times once
        # both sides actually have one, otherwise every fresh translation reports
        # as stale and the signal is useless.
        translated_at = last_commit_epoch(dst.relative_to(REPO).as_posix())
        if translated_at and last_commit_epoch(rel) > translated_at:
            stale.append((n, rel))

    if args.list:
        rows = {"todo": todo, "done": done, "stale": stale, "residual": residual}[args.list]
        for n, rel in sorted(rows, reverse=True):
            print(f"{n:6d}  {rel}")
        return 0

    print(f"{'role':<22} {'files':>7} {'done':>7} {'lines':>9} {'done':>9}  {'%':>5}")
    print("-" * 66)
    tf = td = tl = tdl = 0
    for r in sorted(by_role, key=lambda k: -by_role[k][1]):
        f, ln, df, dl = by_role[r]
        pct = 100.0 * dl / ln if ln else 100.0
        print(f"{r:<22} {f:>7} {df:>7} {ln:>9} {dl:>9}  {pct:>4.0f}%")
        tf += f; td += df; tl += ln; tdl += dl
    print("-" * 66)
    pct = 100.0 * tdl / tl if tl else 100.0
    print(f"{'TOTAL':<22} {tf:>7} {td:>7} {tl:>9} {tdl:>9}  {pct:>4.0f}%")

    if stale:
        print(f"\nstale (original newer than translation): {len(stale)}")
    if residual:
        print(f"residual Hangul inside translations: {len(residual)} file(s)")
        for n, rel in sorted(residual, reverse=True)[:10]:
            print(f"  {n:5d} lines  {en_path(rel).relative_to(REPO).as_posix()}")

    if args.check and residual:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
