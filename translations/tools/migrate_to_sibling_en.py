#!/usr/bin/env python3
"""One-time migration: move translations out of the mirrored ``translations/en/``
tree and into ``<name>_en.<ext>`` siblings living directly next to the upstream
originals they translate.

Rationale, risks, and the checklist this script is part of are recorded in
``ESQlabs_curation_and_extension.md`` at the repo root. This script is kept
after running (not deleted) so the exact transform stays reproducible and
auditable, the same way the rest of ``translations/tools/`` is kept.

Two kinds of same-directory sibling reference need to move together with the
file, or a translated Shiny app would silently source the Korean original
lying in the same directory instead of its own translated model file:

1. Markdown/HTML links -- ``[text](target)``, ``href="target"``,
   ``src="target"`` -- resolved against the file's OLD location, remapped
   through the old-to-new table if the target is itself a migrated file, then
   re-emitted relative to the file's NEW location. A file's own "translation
   of X" attribution backlink is excluded from the remap -- it must keep
   pointing at the real original, not at itself -- mirroring the ``own_original``
   exclusion in ``relink.py``.
2. Bare literal filename mentions in the same directory -- ``source("x.R")``,
   ``MODEL_FILE <- "x.R"``, a plain-text ``## Files`` directory listing --
   caught with an exact-string substitution of each sibling's OLD basename to
   its NEW basename, scoped to files that were previously in the same
   directory. A file's own original basename is excluded from its own
   substitution table for the same reason as (1).

Usage
-----
    python3 translations/tools/migrate_to_sibling_en.py --dry-run [dir ...]
    python3 translations/tools/migrate_to_sibling_en.py [dir ...]

With no directory arguments, migrates every file under ``translations/en/``.
A directory argument restricts to translations whose *upstream* relative path
starts with that directory (e.g. ``achalasia``) -- useful for testing on one
directory before running the full tree.
"""

from __future__ import annotations

import argparse
import io
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"

TEXT_EXTS = {".md", ".r", ".py", ".txt", ".json", ".dot"}
BINARY_EXTS = {".svg", ".png"}

# A plain [text](target) link and a "linked image" [![alt](imgtarget)](target)
# are combined into ONE pattern so a single left-to-right non-overlapping
# re.sub() pass consumes the whole nested construct atomically. Matching them
# with two separate patterns (or two separate .sub() passes) leaves the outer
# link's target unprocessed: after a naive [text](target) pattern consumes
# "[![alt](imgtarget)" as `text`="![alt](imgtarget)"` up to the FIRST `]`, the
# trailing "](target)" has no preceding unmatched "[" left for a second pass
# to match against, so the outer target -- often a same-directory image whose
# relative depth also changed in this migration -- is silently left stale.
#
# Neither alternative excludes https:/mailto: targets at the REGEX stage
# (unlike relink.py's simpler patterns, which only ever see one target at a
# time). A shields.io badge is exactly [![alt](https://img.shields.io/...
# )](../local/page) -- an external inner image wrapped by a LOCAL outer link.
# Excluding https: targets from the image alternative's lookahead meant it
# never matched at all here (the inner target fails the lookahead), so the
# whole construct fell through to nothing, and the outer link -- the one that
# actually needed its relative depth recomputed -- was silently left stale.
# Four files had exactly this shape and went unnoticed until check_links.py's
# post-migration run caught them; scheme filtering now happens per-target
# inside remap_target(), after matching, not inside the regex.
LINK_PATTERN = re.compile(
    r"\[!\[(?P<alt>[^\]]*)\]\((?P<imgtarget>[^)\s]+)\)\]\((?P<linktarget>[^)\s]+)\)"
    r"|"
    r"\[(?P<text>[^\]]*)\]\((?P<target>[^)\s]+)\)"
)
HTML_ATTR = re.compile(r'((?:href|src)=")([^"]+)(")')
EXTERNAL = re.compile(r"^(?:https?:|mailto:)")


def new_rel(rel: Path) -> Path:
    """Upstream-relative old path -> upstream-relative new path."""
    new_name = f"{rel.stem}_en{rel.suffix}"
    return rel.parent / new_name


def build_rel_map() -> dict[Path, Path]:
    rels = [p.relative_to(EN) for p in EN.rglob("*") if p.is_file()]
    return {rel: new_rel(rel) for rel in rels}


def split_frag(target: str) -> tuple[str, str]:
    if "#" in target:
        t, _, frag = target.partition("#")
        return t, frag
    return target, ""


def resolve_and_remap(target: str, old_file: Path, rel_map: dict[Path, Path],
                       own_original_rel: Path) -> Path | None:
    t, _ = split_frag(target)
    if not t:
        return None  # pure same-document anchor, caller leaves it alone
    resolved = (old_file.parent / t).resolve()
    try:
        rel_from_en = resolved.relative_to(EN)
    except ValueError:
        rel_from_en = None
    if rel_from_en is not None:
        if rel_from_en == own_original_rel:
            # A translated file linking to what LOOKS like its own EN path
            # under the old tree shouldn't occur, but if it does, do not loop
            # back through the map -- treat as self and leave unmapped below.
            return REPO / own_original_rel
        return REPO / rel_map.get(rel_from_en, rel_from_en)
    try:
        rel_from_repo = resolved.relative_to(REPO)
    except ValueError:
        return None  # points outside the repo entirely -- leave alone
    # An upstream original (translated or not) keeps its location unchanged,
    # INCLUDING this file's own original -- explicit for clarity, matches the
    # fallthrough anyway, but named so the self-reference case reads plainly.
    return resolved


def rewrite_links(text: str, old_file: Path, new_file: Path,
                   rel_map: dict[Path, Path], own_original_rel: Path) -> str:
    def emit(new_abs: Path, frag: str) -> str:
        rel_path = os.path.relpath(new_abs, start=new_file.parent).replace(os.sep, "/")
        return rel_path + (f"#{frag}" if frag else "")

    def remap_target(target: str) -> str:
        t, frag = split_frag(target)
        if not t or EXTERNAL.match(t):
            return target
        new_abs = resolve_and_remap(target, old_file, rel_map, own_original_rel)
        if new_abs is None:
            return target
        return emit(new_abs, frag)

    def sub_md(m: re.Match) -> str:
        if m.group("imgtarget") is not None:
            alt = m.group("alt")
            new_img = remap_target(m.group("imgtarget"))
            new_link = remap_target(m.group("linktarget"))
            return f"[![{alt}]({new_img})]({new_link})"
        link_text, target = m.group("text"), m.group("target")
        new_target = remap_target(target)
        new_text = new_target if link_text == target else link_text
        return f"[{new_text}]({new_target})"

    def sub_attr(m: re.Match) -> str:
        pre, target, post = m.group(1), m.group(2), m.group(3)
        new_target = remap_target(target)
        if new_target == target:
            return m.group(0)
        return pre + new_target + post

    out = LINK_PATTERN.sub(sub_md, text)
    out = HTML_ATTR.sub(sub_attr, out)
    return out


def sibling_substitutions(rel: Path, rel_map: dict[Path, Path]) -> dict[str, str]:
    """old_basename -> new_basename for every OTHER migrated file that shares
    ``rel``'s upstream directory, excluding ``rel`` itself (own original)."""
    d = rel.parent
    subs = {}
    for other_rel, other_new in rel_map.items():
        if other_rel == rel or other_rel.parent != d:
            continue
        subs[other_rel.name] = other_new.name
    return subs


def apply_bare_substitutions(text: str, rel: Path, rel_map: dict[Path, Path]) -> str:
    subs = sibling_substitutions(rel, rel_map)
    for old_name, new_name in sorted(subs.items(), key=lambda kv: -len(kv[0])):
        # Bounded so "README.md" matches a bare mention but not the tail of
        # "translations/README.md" (a different, non-sibling file whose path
        # happens to end in the same basename) or of a longer identifier --
        # an unanchored substring replace already corrupted exactly this case
        # during testing, turning a correct link to translations/README.md
        # into one pointing at README_en.md.
        pattern = re.compile(r"(?<![\w/.-])" + re.escape(old_name) + r"(?![\w-])")
        text = pattern.sub(new_name, text)
    return text


def migrate_one(rel: Path, rel_map: dict[Path, Path], dry_run: bool) -> str:
    old_file = EN / rel
    new_file = REPO / rel_map[rel]
    ext = old_file.suffix.lower()

    if ext in BINARY_EXTS:
        data = old_file.read_bytes()
        if not dry_run:
            new_file.parent.mkdir(parents=True, exist_ok=True)
            new_file.write_bytes(data)
            old_file.unlink()
        return f"{rel}  ->  {rel_map[rel]}  (binary move)"

    text = old_file.read_text(encoding="utf-8")
    out = text
    if ext == ".md":
        out = rewrite_links(out, old_file, new_file, rel_map, rel)
    if ext in TEXT_EXTS:
        out = apply_bare_substitutions(out, rel, rel_map)

    changed = out != text
    if not dry_run:
        new_file.parent.mkdir(parents=True, exist_ok=True)
        with io.open(new_file, "w", encoding="utf-8", newline="\n") as f:
            f.write(out)
        old_file.unlink()
    return f"{rel}  ->  {rel_map[rel]}" + ("  (links/refs rewritten)" if changed else "")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="*", help="restrict to these upstream-relative directory prefixes")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    rel_map = build_rel_map()
    todo = sorted(rel_map)
    if args.dirs:
        prefixes = tuple(d.rstrip("/") + "/" for d in args.dirs)
        todo = [rel for rel in todo if str(rel).replace(os.sep, "/").startswith(prefixes)]

    if not todo:
        print("nothing matched")
        return 0

    for rel in todo:
        print(migrate_one(rel, rel_map, args.dry_run))

    print(f"\n{'would migrate' if args.dry_run else 'migrated'} {len(todo)} file(s)")

    if not args.dry_run and not args.dirs:
        # Full-tree run: the translations/en/ tree is now empty, remove it.
        remaining = [p for p in EN.rglob("*") if p.is_file()]
        if remaining:
            print(f"NOT removing translations/en/: {len(remaining)} file(s) remain "
                  f"(partial run?)")
        else:
            import shutil
            shutil.rmtree(EN)
            print("removed now-empty translations/en/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
