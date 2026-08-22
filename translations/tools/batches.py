#!/usr/bin/env python3
"""Split the remaining translation work into batches for dispatch.

Reads the pending list from translation_status.py, groups by file role, and
prints batches of a chosen size so each can be handed to one translator.
Batch sizes differ by role because the cost per file does: a Path A file only
sends its Korean lines, while a per-disease README is translated whole.

Usage
-----
    python3 translations/tools/batches.py                 # per-role summary
    python3 translations/tools/batches.py --role readme   # list its batches
    python3 translations/tools/batches.py --role refs --batch 5
    python3 translations/tools/batches.py --role dot --index 3   # one batch
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# role -> (matcher, default files per batch, translation path)
ROLES: dict[str, tuple[object, int, str]] = {
    "refs":   (lambda p: p.endswith("_references.md"), 4, "A"),
    "readme": (lambda p: p.endswith("README.md"), 3, "B"),
    "shiny":  (lambda p: p.endswith(".R") and "shiny" in p.rsplit("/", 1)[-1], 3, "A"),
    "dot":    (lambda p: p.endswith(".dot"), 3, "A"),
    "rmodel": (lambda p: p.endswith(".R") and "shiny" not in p.rsplit("/", 1)[-1], 6, "A"),
    "misc":   (lambda p: p.endswith((".py", ".txt", ".json")), 6, "A"),
}


def pending() -> list[tuple[int, str]]:
    out = subprocess.run(
        [sys.executable, str(REPO / "translations/tools/translation_status.py"),
         "--list", "todo"],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8",
    ).stdout
    rows = []
    for line in out.splitlines():
        m = re.match(r"\s*(\d+)\s+(\S.*)$", line)
        if m:
            rows.append((int(m.group(1)), m.group(2).strip()))
    return rows


def classify(rows: list[tuple[int, str]]) -> dict[str, list[tuple[int, str]]]:
    """Assign each pending file to exactly one role, first match wins."""
    order = ["refs", "readme", "shiny", "dot", "rmodel", "misc"]
    buckets: dict[str, list[tuple[int, str]]] = {r: [] for r in order}
    for n, path in rows:
        for role in order:
            if ROLES[role][0](path):
                buckets[role].append((n, path))
                break
    return buckets


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--role", choices=list(ROLES))
    ap.add_argument("--batch", type=int, help="files per batch (default per role)")
    ap.add_argument("--index", type=int, help="print only this batch (1-based)")
    args = ap.parse_args()

    buckets = classify(pending())

    if not args.role:
        print(f"{'role':<8} {'path':<5} {'files':>6} {'lines':>8} {'batches':>8}")
        print("-" * 40)
        for role, rows in buckets.items():
            size = args.batch or ROLES[role][1]
            nb = (len(rows) + size - 1) // size
            print(f"{role:<8} {ROLES[role][2]:<5} {len(rows):>6} "
                  f"{sum(n for n, _ in rows):>8} {nb:>8}")
        print("-" * 40)
        print(f"{'TOTAL':<8} {'':<5} {sum(len(r) for r in buckets.values()):>6} "
              f"{sum(n for r in buckets.values() for n, _ in r):>8}")
        return 0

    rows = buckets[args.role]
    size = args.batch or ROLES[args.role][1]
    batches = [rows[i:i + size] for i in range(0, len(rows), size)]
    chosen = ([batches[args.index - 1]] if args.index
              else batches)
    start = args.index or 1
    for k, batch in enumerate(chosen, start):
        print(f"# batch {k}/{len(batches)}  role={args.role} "
              f"path={ROLES[args.role][2]}  "
              f"{sum(n for n, _ in batch)} Korean lines")
        for n, path in batch:
            print(f"{path}\t{n}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
