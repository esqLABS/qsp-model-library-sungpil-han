#!/usr/bin/env python3
"""Line-level translation for files that are structure-dominant.

Most non-README files in this library are mostly *not* Korean: a
``*_references.md`` is English citations, PMIDs, and URLs with Korean section
headers; a ``*_shiny_app.R`` is R code with Korean UI labels; a ``*_qsp_model.dot``
is Graphviz syntax with Korean node labels. Retyping such a file to translate it
risks corrupting exactly the parts that must not change — and there are 29,460
PMIDs in this repository.

So instead: ``extract`` pulls out only the lines containing Hangul, a translator
returns those lines, and ``apply`` substitutes them back into an otherwise
byte-identical copy under ``translations/en/``. Every line that was not extracted
is copied verbatim, so a PMID cannot be altered by this path.

Prose-dominant files (the per-disease ``README.md``) are translated whole
instead, because English prose does not line up with Korean prose line by line.

Usage
-----
    # 1. write the Korean lines of a file to a numbered TSV
    python3 translations/tools/lineio.py extract <path> [-o out.tsv]

    # 2. a translator edits column 2 of that TSV in place, then:
    python3 translations/tools/lineio.py apply <path> <tsv>

    # check an already-built translation is line-aligned and Hangul-free
    python3 translations/tools/lineio.py verify <path>

TSV format is ``<line number>\\t<text>``. Line numbers are 1-based and must match
the original exactly; a missing or extra number is an error, not a silent skip.
A tab inside the text itself is escaped as ``\\t`` on the way out and restored on
the way in, so the format stays two columns.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
EN = REPO / "translations" / "en"

HANGUL = re.compile(r"[가-힣ᄀ-ᇿ㄰-㆏]")


def rel_of(path: str) -> str:
    """Repo-relative POSIX path for ``path``, which may be absolute."""
    p = Path(path)
    if p.is_absolute():
        p = p.resolve().relative_to(REPO)
    return p.as_posix()


def read_lines(path: Path) -> list[str]:
    # newline="" keeps the file's own line endings out of the line text.
    return path.read_text(encoding="utf-8").split("\n")


def esc(s: str) -> str:
    return s.replace("\\", "\\\\").replace("\t", "\\t")


def unesc(s: str) -> str:
    out, i = [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s):
            nxt = s[i + 1]
            if nxt == "t":
                out.append("\t"); i += 2; continue
            if nxt == "\\":
                out.append("\\"); i += 2; continue
        out.append(s[i]); i += 1
    return "".join(out)


TRAILING_WS = re.compile(r"[ \t]*$")


def restore_trailing_ws(original: str, translated: str) -> str:
    """Give ``translated`` the trailing whitespace run that ``original`` had.

    Trailing whitespace is load-bearing in markdown: two spaces at end of line is
    a hard line break, and several ``*_references.md`` files use it to hold each
    citation's four lines together. Most editors -- and the Write tool a
    translator uses to save the TSV -- strip it silently, and ``verify`` cannot
    see the loss because it only compares lines that had no Korean. So rather
    than asking anyone to preserve invisible characters by hand, reattach them
    from the source. A translator who deliberately wants different trailing
    whitespace cannot express that here, which is the right trade: in this
    repository it is never intended.
    """
    if translated.strip() == "":
        return translated
    return TRAILING_WS.sub("", translated) + TRAILING_WS.search(original).group()


def cmd_extract(args: argparse.Namespace) -> int:
    rel = rel_of(args.path)
    lines = read_lines(REPO / rel)
    rows = [(i, ln) for i, ln in enumerate(lines, 1) if HANGUL.search(ln)]
    body = "".join(f"{i}\t{esc(ln)}\n" for i, ln in rows)
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(body, encoding="utf-8", newline="\n")
        print(f"{rel}: {len(rows)} Korean lines of {len(lines)} -> {out}")
    else:
        sys.stdout.write(body)
    return 0


def cmd_apply(args: argparse.Namespace) -> int:
    rel = rel_of(args.path)
    src = REPO / rel
    lines = read_lines(src)
    expected = {i for i, ln in enumerate(lines, 1) if HANGUL.search(ln)}

    repl: dict[int, str] = {}
    for lineno, raw in enumerate(
        Path(args.tsv).read_text(encoding="utf-8").splitlines(), 1
    ):
        if not raw.strip():
            continue
        num, tab, text = raw.partition("\t")
        if not tab or not num.strip().isdigit():
            print(f"{args.tsv}:{lineno}: not '<number><TAB><text>'", file=sys.stderr)
            return 2
        repl[int(num)] = unesc(text)

    missing = sorted(expected - repl.keys())
    extra = sorted(repl.keys() - expected)
    if missing:
        print(f"{rel}: {len(missing)} Korean line(s) not translated: "
              f"{missing[:10]}", file=sys.stderr)
    if extra:
        print(f"{rel}: {len(extra)} line number(s) that hold no Korean: "
              f"{extra[:10]}", file=sys.stderr)
    if (missing or extra) and not args.force:
        return 2

    left = [HANGUL.search(t) for t in repl.values() if HANGUL.search(t)]
    if left and not args.force:
        print(f"{rel}: {len(left)} translated line(s) still contain Hangul",
              file=sys.stderr)
        return 2

    for i, text in repl.items():
        lines[i - 1] = restore_trailing_ws(lines[i - 1], text)

    dst = EN / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(f"wrote translations/en/{rel} "
          f"({len(repl)} of {len(expected)} Korean lines replaced)")
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    rel = rel_of(args.path)
    src, dst = REPO / rel, EN / rel
    if not dst.exists():
        print(f"{rel}: no translation at translations/en/{rel}", file=sys.stderr)
        return 1
    a, b = read_lines(src), read_lines(dst)
    problems = []
    if len(a) != len(b):
        problems.append(f"line count {len(a)} -> {len(b)}")
    left = sum(1 for ln in b if HANGUL.search(ln))
    if left:
        problems.append(f"{left} line(s) still contain Hangul")
    drifted = sum(
        1 for x, y in zip(a, b) if not HANGUL.search(x) and x != y
    )
    if drifted:
        problems.append(f"{drifted} non-Korean line(s) changed")
    # Trailing whitespace is a markdown hard line break. It is invisible, most
    # editors strip it, and the non-Korean-line check above cannot see it going
    # missing from a translated line.
    lost_ws = sum(
        1 for x, y in zip(a, b)
        if TRAILING_WS.search(x).group() and not TRAILING_WS.search(y).group()
    )
    if lost_ws:
        problems.append(f"{lost_ws} line(s) lost trailing whitespace "
                        f"(markdown hard line break)")
    if problems:
        print(f"FAIL {rel}: " + "; ".join(problems))
        return 1
    print(f"OK   {rel}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    e = sub.add_parser("extract", help="write Korean lines as a numbered TSV")
    e.add_argument("path")
    e.add_argument("-o", "--out")
    e.set_defaults(fn=cmd_extract)

    a = sub.add_parser("apply", help="substitute translated lines back in")
    a.add_argument("path")
    a.add_argument("tsv")
    a.add_argument("--force", action="store_true",
                   help="write despite missing/extra lines or residual Hangul")
    a.set_defaults(fn=cmd_apply)

    v = sub.add_parser("verify", help="check an existing translation")
    v.add_argument("path")
    v.set_defaults(fn=cmd_verify)

    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
