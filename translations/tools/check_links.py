#!/usr/bin/env python3
"""Check that relative links inside the translations resolve to real files.

Every translated file is a ``<name>_en.<ext>`` sibling next to the original it
was made from, and a link inside it can point at either one -- the untranslated
original, or another file's ``_en`` translation once that exists. Getting the
target or the relative path wrong produces a link that looks plausible and goes
nowhere, so it is checked rather than trusted.

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


HEADING = re.compile(r"^\s{0,3}(#{1,6})\s+(.*?)\s*#*\s*$", re.MULTILINE)
INLINE_MD = re.compile(r"\[([^\]]*)\]\([^)]*\)|[*_`~]")


def slug(heading: str) -> str:
    """GitHub's anchor slug for a heading.

    Lower-case, markdown stripped, anything that is not a letter, digit, space,
    hyphen, or underscore dropped, then spaces to hyphens. Close enough to
    GitHub's own rule to catch a heading that was translated while the anchor
    pointing at it was not.
    """
    text = INLINE_MD.sub(lambda m: m.group(1) or "", heading)
    text = "".join(
        c for c in text.lower()
        if c.isalnum() or c in " -_"
    )
    return text.strip().replace(" ", "-")


def anchors(text: str) -> set[str]:
    """Every anchor a reader could jump to in this document."""
    out: set[str] = set()
    counts: dict[str, int] = {}
    for _, heading in HEADING.findall(strip_code(text)):
        base = slug(heading)
        n = counts.get(base, 0)
        counts[base] = n + 1
        out.add(base if n == 0 else f"{base}-{n}")
    # Explicit HTML anchors are valid jump targets too.
    out |= set(re.findall(r'<a\s+[^>]*name="([^"]+)"', text))
    out |= set(re.findall(r'\bid="([^"]+)"', text))
    return out


def check(path: Path) -> list[str]:
    """Return a list of unresolvable link targets in ``path``."""
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return []
    bad = []
    own = anchors(text) if path.suffix.lower() == ".md" else set()
    for t in sorted(targets(text)):
        if PLACEHOLDER.search(t):
            continue
        file_part, _, frag = t.partition("#")
        if not file_part:
            # A same-document anchor. Translating a heading breaks these, and
            # nothing else checks them.
            if frag and frag not in own:
                bad.append(f"#{frag} (no such heading in this file)")
            continue
        target = (path.parent / file_part).resolve()
        if not target.exists():
            bad.append(t)
            continue
        # A fragment pointing into another markdown file is just as breakable as
        # one pointing into this file, and translating the target's headings is
        # what breaks it. Only check targets that are themselves translations: an
        # anchor into an untranslated original is still Korean-slugged and
        # correct.
        if frag and target.suffix.lower() == ".md":
            if not target.stem.endswith("_en"):
                continue
            if frag not in anchors(target.read_text(encoding="utf-8")):
                bad.append(f"{t} (no such heading in {file_part})")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    args = ap.parse_args()

    if args.paths:
        files = [Path(p).resolve() for p in args.paths]
    else:
        files = sorted(
            p for p in REPO.rglob("*_en.*")
            if p.is_file() and p.suffix.lower() in {".md", ".html"}
            and ".git" not in p.parts
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
