#!/usr/bin/env python3
"""Check that relative links inside the translations resolve to real files.

Every translated file sits deeper in the tree than the original it was made from,
so each relative link has to gain the right number of ``../`` segments. Getting
that count wrong is the single easiest mistake to make here and it produces links
that look plausible and go nowhere, so it is checked rather than trusted.

Usage
-----
    python3 translations/tools/check_links.py                 # everything
    python3 translations/tools/check_links.py <path> [<path>…] # specific files
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"

MD_LINK = re.compile(r"\]\((?!https?:|mailto:|#)([^)\s]+)")
HTML_ATTR = re.compile(r'(?:href|src)="(?!https?:|mailto:|#)([^"]+)"')

# Link targets we do not expect to exist on disk.
PLACEHOLDER = re.compile(r"[<>{}*]|^\.$")

FENCE = re.compile(r"^\s*(```|~~~).*?^\s*\1", re.MULTILINE | re.DOTALL)
INLINE_CODE = re.compile(r"`[^`\n]*`")


def strip_code(text: str) -> str:
    """Remove fenced blocks and inline code spans.

    A path inside a code sample is documentation, not a link — the shell snippets
    in these files are full of ``path/to/png`` and ``<disease>/<abbr>.dot``, and
    checking those as links only produces noise.
    """
    text = FENCE.sub("", text)
    return INLINE_CODE.sub("", text)


def targets(text: str) -> set[str]:
    text = strip_code(text)
    return set(MD_LINK.findall(text)) | set(HTML_ATTR.findall(text))


def check(path: Path) -> list[str]:
    """Return a list of unresolvable link targets in ``path``."""
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return []
    bad = []
    for t in sorted(targets(text)):
        if PLACEHOLDER.search(t):
            continue
        target = (path.parent / t.split("#")[0]).resolve()
        if not target.exists():
            bad.append(t)
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    args = ap.parse_args()

    if args.paths:
        files = [Path(p).resolve() for p in args.paths]
    else:
        files = sorted(
            p for p in EN.rglob("*")
            if p.is_file() and p.suffix.lower() in {".md", ".html"}
        )
    if not files:
        print("no translated files to check")
        return 0

    failed = 0
    for f in files:
        bad = check(f)
        rel = f.relative_to(REPO).as_posix()
        if bad:
            failed += 1
            print(f"FAIL {rel}: {len(bad)} unresolvable link(s)")
            for t in bad:
                print(f"       {t}")
    total = len(files)
    if failed:
        print(f"\n{failed} of {total} file(s) have broken links")
        return 1
    print(f"OK   all relative links resolve in {total} file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
