#!/usr/bin/env python3
"""Check that invariant tokens survive a translation unchanged.

A translation is allowed to change words. It is not allowed to change a PMID, a
DOI, a URL, or a measured value — and those are exactly the things that are
expensive to get wrong here and invisible on a read-through. This compares the
multiset of such tokens in the original against the translation and reports drift.

Only numbers that carry data are checked — decimals, exponents, and integers of
three or more digits. Small integers are skipped deliberately: they change form in
honest translation (제1형 becomes "type I", 2회 becomes "twice"), and checking them
buries real drift under false positives. Doses, thresholds, parameter values, and
cohort sizes all still fall inside the check.

Files translated by the whole-file path (Path B) reflow prose, so ``--loose``
relaxes the number check from a multiset comparison to a presence check while
keeping PMIDs, DOIs, and URLs strict.

Usage
-----
    python3 translations/tools/check_tokens.py                    # all translations
    python3 translations/tools/check_tokens.py <path> [<path>…]
    python3 translations/tools/check_tokens.py --loose <path>     # Path B files
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"
EXCEPTIONS = REPO / "translations" / "data" / "token_exceptions.tsv"


def load_exceptions() -> dict[str, set[str]]:
    """Reviewed, accepted token differences, keyed by repo-relative path.

    Some honest translations do change digits: Korean 1000분의 1 is "a
    thousandth" in idiomatic English. Recording those here keeps the gate green
    by default, so real drift in a later batch is not lost among known noise.
    """
    out: dict[str, set[str]] = {}
    if not EXCEPTIONS.exists():
        return out
    for line in EXCEPTIONS.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            out.setdefault(parts[0].strip(), set()).add(parts[1].strip())
    return out


PATTERNS = {
    "PMID": re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)"),
    "DOI": re.compile(r"\b10\.\d{4,9}/[-._;()/:A-Za-z0-9]+"),
    "URL": re.compile(r"https?://[^\s)\"'>\]]+"),
    # A number that carries DATA rather than grammar: a decimal, an exponent, or
    # an integer of three or more digits. Single- and double-digit integers are
    # excluded on purpose, because they legitimately change form in translation
    # -- Korean 제1형 becomes "type I", 3대 becomes "the three", 2회 becomes
    # "twice" -- and checking them buries the real drift in false positives.
    # Parameter values, doses, thresholds, and cohort sizes all survive the cut.
    #
    # The boundaries list ASCII rather than using \w: Python's \w matches Hangul,
    # so a \w boundary would skip the number in "37주" while counting it in
    # "37 weeks" and report every translated unit as invented. Excluding only
    # ASCII word characters still keeps identifiers like abc123 out.
    # The trailing boundary must allow a sentence-final full stop: with a plain
    # `(?![A-Za-z0-9_.])` the value in "the ratio is 0.90." is invisible, which
    # reports a number as missing when it is present. So reject a following digit
    # (to keep version strings like v1.2.3 out) but accept a following period.
    "number": re.compile(
        r"(?<![A-Za-z0-9_.])"
        r"(?:\d+\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+|\d{3,})"
        r"(?![A-Za-z0-9_])(?!\.\d)"
    ),
}


# Escape sequences that act as whitespace in the file's own rendering: Graphviz
# label breaks and the usual string escapes.
ESCAPE = re.compile(r"\\[nlrt]")


def tokens(text: str, kind: str) -> Counter:
    if kind == "number":
        # A value written straight after a label break -- "\n120-150 mg/kg" --
        # is preceded by the literal letter n, so the word boundary would reject
        # it and report the value as missing. These escapes are whitespace as far
        # as the reader is concerned, so treat them as whitespace here too.
        text = ESCAPE.sub(" ", text)
    return Counter(PATTERNS[kind].findall(text))


def check(rel: str, loose: bool, allowed: set[str] | None = None) -> list[str]:
    allowed = allowed or set()
    src, dst = REPO / rel, EN / rel
    if not dst.exists():
        return [f"no translation at translations/en/{rel}"]
    try:
        a = src.read_text(encoding="utf-8")
        b = dst.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError) as exc:
        return [f"unreadable: {exc}"]

    problems = []
    for kind in ("PMID", "DOI", "URL"):
        ca, cb = tokens(a, kind), tokens(b, kind)
        if ca != cb:
            lost = sorted((ca - cb).elements())[:6]
            gained = sorted((cb - ca).elements())[:6]
            detail = []
            if lost:
                detail.append(f"lost {lost}")
            if gained:
                detail.append(f"gained {gained}")
            problems.append(f"{kind}: " + "; ".join(detail))

    ca, cb = tokens(a, "number"), tokens(b, "number")
    if loose:
        missing = sorted((set(ca) - set(cb)) - allowed)
        if missing:
            problems.append(f"number: {len(missing)} value(s) absent "
                            f"e.g. {missing[:8]}")
    else:
        lost = [x for x in (ca - cb).elements() if x not in allowed]
        gained = [x for x in (cb - ca).elements() if x not in allowed]
        if lost or gained:
            detail = []
            if lost:
                detail.append(f"lost {sorted(lost)[:8]}")
            if gained:
                detail.append(f"gained {sorted(gained)[:8]}")
            problems.append("number: " + "; ".join(detail))
    return problems


def rels(paths: list[str]) -> list[str]:
    if paths:
        out = []
        for p in paths:
            q = Path(p)
            if q.is_absolute():
                q = q.resolve().relative_to(REPO)
            q = Path(q.as_posix().removeprefix("translations/en/"))
            out.append(q.as_posix())
        return out
    return sorted(
        p.relative_to(EN).as_posix()
        for p in EN.rglob("*")
        if p.is_file() and p.suffix.lower() not in
        {".png", ".svg", ".jpg", ".webp", ".pyc"}
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="*")
    ap.add_argument("--loose", action="store_true",
                    help="presence-only number check, for whole-file translations")
    args = ap.parse_args()

    exceptions = load_exceptions()
    failed = 0
    checked = 0
    for rel in rels(args.paths):
        if not (REPO / rel).exists():
            continue
        checked += 1
        problems = check(rel, args.loose, exceptions.get(rel))
        if problems:
            failed += 1
            print(f"FAIL {rel}")
            for p in problems:
                print(f"       {p}")
    if not checked:
        print("nothing to check")
        return 0
    if failed:
        print(f"\n{failed} of {checked} file(s) show token drift")
        return 1
    print(f"OK   invariant tokens preserved in {checked} file(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
