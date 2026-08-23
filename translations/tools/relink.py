#!/usr/bin/env python3
"""Re-point links in the translations at translated siblings once they exist.

The translation lands file by file and out of order, so a README translated on
Monday links three levels up to the Korean original of a map that gets translated
on Tuesday. Nothing is broken -- the link resolves -- but the reader of an English
README clicks through to a Korean-labelled image, and it silently stays that way
forever.

This walks every translated file, and wherever a link points at an upstream file
that now has a translation of its own, repoints it at the translation. It is
idempotent, so it can be run after every batch.

Usage
-----
    python3 translations/tools/relink.py            # report only
    python3 translations/tools/relink.py --write     # apply
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"

# Captures the link TEXT as well as the target. Several READMEs write a
# "parent library" backlink whose visible text is the literal relative path,
# e.g. [../../../README.md](../../../README.md) -- text and target identical.
# relink.py used to rewrite only the target half, leaving the visible text
# pointing at the stale path once the target moved. When text == old target,
# the text is now updated to the new target too.
MD_LINK = re.compile(r"\[([^\]]*)\]\((?!https?:|mailto:|#)([^)\s]+)\)")
HTML_ATTR = re.compile(r'((?:href|src)=")(?!https?:|mailto:|#)([^"]+)(")')

# Only these get repointed. A link to a .dot/.R/.py source is arguably better
# left on the original -- but no: the whole point is that the English reader
# wants the English file. Images especially, since a Korean map render is the
# most jarring.
INTERESTING = {".svg", ".png", ".dot", ".md", ".r", ".py", ".txt", ".json"}

# Files that are generated, not hand-written. build_root_readme.py already
# resolves its own links against what exists, and would overwrite anything done
# here on its next run.
GENERATED = {EN / "README.md"}


def rewrite_one(src_file: Path, target: str) -> str | None:
    """Return a replacement for ``target``, or None to leave it alone."""
    frag = ""
    if "#" in target:
        target, _, frag = target.partition("#")
        frag = "#" + frag
    if not target:
        return None

    resolved = (src_file.parent / target).resolve()

    # Never touch a file's own attribution backlink. Most translated files open
    # with "English translation of the upstream <link>", pointing at the very
    # original they were made from. Repointing that at the translation would make
    # the file cite itself, which is both wrong and silently circular.
    own_original = REPO / src_file.relative_to(EN)
    if resolved == own_original.resolve():
        return None
    # Already inside the translations tree? Nothing to do.
    try:
        resolved.relative_to(EN)
        return None
    except ValueError:
        pass
    try:
        rel = resolved.relative_to(REPO)
    except ValueError:
        return None

    translated = EN / rel
    if not translated.exists():
        return None

    new = os.path.relpath(translated, src_file.parent).replace(os.sep, "/")
    return new + frag


def process(path: Path, write: bool) -> list[tuple[str, str]]:
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return []
    changes: list[tuple[str, str]] = []

    def sub_md(m: re.Match) -> str:
        link_text, target = m.group(1), m.group(2)
        if Path(target.split("#")[0]).suffix.lower() not in INTERESTING:
            return m.group(0)
        new = rewrite_one(path, target)
        if new is None or new == target:
            return m.group(0)
        changes.append((target, new))
        # If the visible text was literally the old path, keep it literal.
        new_text = new if link_text == target else link_text
        return f"[{new_text}]({new})"

    def sub_attr(m: re.Match) -> str:
        pre, target, post = m.group(1), m.group(2), m.group(3)
        if Path(target.split("#")[0]).suffix.lower() not in INTERESTING:
            return m.group(0)
        new = rewrite_one(path, target)
        if new is None or new == target:
            return m.group(0)
        changes.append((target, new))
        return pre + new + post

    out = MD_LINK.sub(sub_md, text)
    out = HTML_ATTR.sub(sub_attr, out)
    if changes and write:
        path.write_text(out, encoding="utf-8", newline="\n")
    return changes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="apply the changes")
    args = ap.parse_args()

    files = sorted(
        p for p in EN.rglob("*")
        if p.is_file() and p.suffix.lower() in {".md", ".html"}
        and p not in GENERATED
    )
    total = 0
    touched = 0
    for f in files:
        changes = process(f, args.write)
        if not changes:
            continue
        touched += 1
        total += len(changes)
        print(f"{f.relative_to(REPO).as_posix()}: {len(changes)} link(s)")
        for old, new in changes[:6]:
            print(f"    {old}  ->  {new}")

    verb = "repointed" if args.write else "would repoint"
    print(f"\n{verb} {total} link(s) across {touched} file(s) "
          f"(scanned {len(files)})")
    if total and not args.write:
        print("re-run with --write to apply")
    return 0


if __name__ == "__main__":
    sys.exit(main())
