#!/usr/bin/env python3
"""Render the translated Graphviz maps under ``translations/en/``.

A translated ``*_qsp_model.dot`` needs its own ``.svg`` and ``.png``; the
upstream images still carry Korean node labels. Upstream renders PNGs at 150 dpi
(see ``CLAUDE.md``), so this matches that.

Usage
-----
    python3 translations/tools/render_maps.py              # render what is stale
    python3 translations/tools/render_maps.py --force      # re-render everything
    python3 translations/tools/render_maps.py <path.dot>    # just one
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"

# Graphviz is not always on PATH on Windows.
FALLBACK_DOT = [
    r"C:\Program Files\Graphviz\bin\dot.exe",
    r"C:\Program Files (x86)\Graphviz\bin\dot.exe",
    "/usr/bin/dot",
    "/usr/local/bin/dot",
]


def find_dot() -> str:
    found = shutil.which("dot")
    if found:
        return found
    for cand in FALLBACK_DOT:
        if os.path.exists(cand):
            return cand
    sys.exit("graphviz 'dot' not found; install it or put it on PATH")


def render(dot: str, src: Path, force: bool) -> tuple[int, int]:
    """Render ``src`` to .svg and .png beside it. Returns (made, skipped)."""
    made = skipped = 0
    for fmt, extra in (("svg", []), ("png", ["-Gdpi=150"])):
        out = src.with_suffix("." + fmt)
        if out.exists() and not force and out.stat().st_mtime >= src.stat().st_mtime:
            skipped += 1
            continue
        before = out.stat().st_mtime if out.exists() else 0
        proc = subprocess.run(
            [dot, f"-T{fmt}", *extra, str(src), "-o", str(out)],
            capture_output=True, text=True,
        )
        wrote = out.exists() and out.stat().st_mtime > before and out.stat().st_size
        if not wrote:
            print(f"FAIL {src.relative_to(REPO)} -> {fmt}: "
                  f"{proc.stderr.strip()[:200]}", file=sys.stderr)
            continue
        if proc.returncode != 0 or proc.stderr.strip():
            # dot writes the image and still complains: a missing font, or an
            # `init_rank` ranking problem in the graph. Several upstream maps do
            # this before translation too, so it is reported as a warning rather
            # than treated as a failed render -- see ../UPSTREAM_ISSUES.md.
            first = proc.stderr.strip().splitlines()[0] if proc.stderr.strip() else ""
            print(f"warn {src.relative_to(REPO)} -> {fmt}: {first[:160]}",
                  file=sys.stderr)
        made += 1
    return made, skipped


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*", help="specific .dot files (default: all)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    dot = find_dot()
    srcs = [Path(p).resolve() for p in args.paths] or sorted(EN.rglob("*.dot"))
    if not srcs:
        print("no translated .dot files yet")
        return 0

    total_made = total_skipped = 0
    for src in srcs:
        made, skipped = render(dot, src, args.force)
        total_made += made
        total_skipped += skipped
    print(f"rendered {total_made} image(s), {total_skipped} already current, "
          f"from {len(srcs)} .dot file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
