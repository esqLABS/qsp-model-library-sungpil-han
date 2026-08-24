#!/usr/bin/env python3
"""Find markdown tables whose rows do not all have the same number of cells.

`CLAUDE.md` records an incident where a stray `**` corrupted the root README's
gallery table, and `scripts/fix_readme_table.py` exists to stop that recurring --
but only for the root README. The 343 per-disease READMEs have no such guard, and a
table that has lost a column renders as visible garbage.

The subtlety is GFM's actual cell-splitting rule, which is easy to get backwards:

* pipes are split **before** inline parsing, so a pipe inside a code span
  (`` `(1+|mu|)` ``) **does** split the cell and break the row;
* only a backslash-escaped pipe (`\\|`) does not.

Getting this wrong in either direction is worse than not checking: masking code
spans hides real breakage, and treating `\\|` as a separator invents breakage that
is not there. Both mistakes were made while writing this.

Usage
-----
    python3 translations/tools/check_tables.py                  # originals
    python3 translations/tools/check_tables.py --translations    # translations
    python3 translations/tools/check_tables.py <path> [<path>…]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

DELIM = re.compile(r"^\s*\|?[\s:-]*-[-\s:|]*\|")
SPLIT = re.compile(r"(?<!\\)\|")


def cells(line: str) -> int:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    return len(SPLIT.split(s))


def check(path: Path) -> list[tuple[int, str]]:
    try:
        lines = path.read_text(encoding="utf-8").split("\n")
    except (UnicodeDecodeError, OSError):
        return []
    out: list[tuple[int, str]] = []
    i = 0
    while i < len(lines):
        if (i + 1 < len(lines) and lines[i].lstrip().startswith("|")
                and DELIM.match(lines[i + 1])):
            hdr, dlm = cells(lines[i]), cells(lines[i + 1])
            if hdr != dlm:
                out.append((i + 2, f"delimiter has {dlm} cells, header has {hdr}"))
            j = i + 2
            while j < len(lines) and lines[j].strip() and "|" in lines[j]:
                if not lines[j].lstrip().startswith("|"):
                    out.append((j + 1, "body row does not start with a pipe"))
                elif cells(lines[j]) != hdr:
                    out.append((j + 1,
                                f"row has {cells(lines[j])} cells, header has {hdr}"))
                j += 1
            i = j
        else:
            i += 1
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--translations", action="store_true",
                    help="check the *_en.md translations instead of the originals")
    args = ap.parse_args()

    if args.paths:
        files = [Path(p).resolve() for p in args.paths]
    elif args.translations:
        files = sorted(p for p in REPO.rglob("*_en.md") if ".git" not in p.parts)
    else:
        files = sorted(p for p in REPO.glob("*/README.md")) \
            + sorted(REPO.glob("*/*_references.md")) + [REPO / "README.md"]

    bad = 0
    for f in files:
        problems = check(f)
        if not problems:
            continue
        bad += 1
        print(f"{f.relative_to(REPO).as_posix()}")
        for ln, why in problems[:4]:
            print(f"    line {ln}: {why}")
        if len(problems) > 4:
            print(f"    … and {len(problems) - 4} more")
    print(f"\n{bad} of {len(files)} file(s) have a ragged table")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
