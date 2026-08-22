#!/usr/bin/env python3
"""Parse-check every R file in the repository and classify the failures.

Two different things carry a ``.R`` extension here, and only one of them is R:

* an **R script** — a Shiny app, or a model whose mrgsolve code sits inside a
  quoted string. ``parse()`` must succeed. If it does not, the file is broken.
* a **pure mrgsolve DSL** file — blocks at top level. ``parse()`` cannot succeed
  and its failure means nothing.

Telling them apart is the whole job, and it is easy to get wrong: the block marker
appears as ``$PROB``, ``[PROB]``, and ``[ PROB ]`` in different files, so a naive
``grep '^\\$'`` mis-reports the bracket forms as broken scripts.

Run this after an upstream pull to see whether a genuinely unparseable file has
appeared. It reads the originals, not the translations — a translation cannot
introduce this class of fault, because ``lineio verify`` proves it changed no line
that lacked Korean.

Usage
-----
    python3 translations/tools/check_r_parses.py
    python3 translations/tools/check_r_parses.py --include-translations
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

# Any mrgsolve block marker, in any of the three spellings the library uses.
BLOCK = re.compile(
    r"^[ \t]*(?:\$|\[[ \t]*)"
    r"(?:PROB|PARAM|CMT|MAIN|ODE|TABLE|SET|GLOBAL|PLUGIN|INIT|CAPTURE"
    r"|OMEGA|SIGMA|PRED|PKMODEL|NMXML|FIXED|THETA|ENV|BLOCK)\b",
    re.MULTILINE,
)

RSCRIPT_CANDIDATES = [
    r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe",
    "Rscript",
]

CHECKER = r"""
files <- readLines(commandArgs(TRUE)[1], warn = FALSE)
for (f in files) {
  r <- tryCatch({ parse(f); NA_character_ },
                error = function(e) conditionMessage(e))
  if (!is.na(r)) cat(sprintf("%s\t%s\n", f, gsub("[\r\n]+", " ", r)))
}
"""


def find_rscript() -> str | None:
    import shutil
    for cand in RSCRIPT_CANDIDATES:
        if Path(cand).exists() or shutil.which(cand):
            return cand
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-translations", action="store_true")
    args = ap.parse_args()

    rscript = find_rscript()
    if not rscript:
        print("Rscript not found; cannot parse-check", file=sys.stderr)
        return 2

    listed = subprocess.run(
        ["git", "ls-files", "*.R"],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8", check=True,
    ).stdout.split()
    files = [f for f in listed
             if args.include_translations or not f.startswith("translations/")]

    tmp = Path(tempfile.mkdtemp())
    (tmp / "files.txt").write_text("\n".join(files), encoding="utf-8")
    (tmp / "check.R").write_text(CHECKER, encoding="utf-8")

    out = subprocess.run(
        [rscript, str(tmp / "check.R"), str(tmp / "files.txt")],
        cwd=REPO, capture_output=True, text=True, encoding="utf-8",
    ).stdout

    dsl, broken = [], []
    for line in out.splitlines():
        if "\t" not in line:
            continue
        f, err = line.split("\t", 1)
        try:
            text = (REPO / f).read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        (dsl if BLOCK.search(text) else broken).append((f, err))

    print(f"checked {len(files)} R file(s)")
    print(f"  {len(dsl)} are mrgsolve DSL — parse() does not apply, not a fault")
    print(f"  {len(broken)} genuinely fail to parse\n")
    for f, err in broken:
        print(f"BROKEN {f}")
        print(f"       {err[:160]}")
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
