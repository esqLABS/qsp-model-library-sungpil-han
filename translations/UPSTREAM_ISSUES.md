# Defects found in the upstream files while translating

Translating a file means reading every line of it, and running the gates means
parsing and rendering it. That turns up bugs in the originals. They are recorded
here rather than fixed, because of the rule in [`README.md`](README.md): we do not
edit upstream files.

Each entry states how it was confirmed to be upstream and not something the
translation introduced. Any of these could be sent upstream as a patch or an issue.

---

## 1. `acute-bacterial-meningitis/abm_mrgsolve_model.R` does not parse

**Line 560.** The whole mrgsolve model is held in a single-quoted R string, and this
comment line inside it contains an apostrophe:

```
// THE PRODUCT.  This single expression is the model's central claim: what
```

The apostrophe closes the string early, so `parse()` fails at `560:54` with
`unexpected symbol`. The file cannot be sourced as written.

**Confirmed upstream:** the original and the translation fail at the identical line
and column, and substituting `model's` → `models` in temp copies of *both* makes
both parse (27 expressions each). The line contains no Hangul, so `lineio verify`
proves the translation left it byte-identical.

**Fix upstream would be:** `model's` → `models`, or switch the enclosing string to
double quotes.

**Not to be confused with** the many `*_mrgsolve_model.R` files that are pure
mrgsolve DSL rather than R scripts — they open a `$PROB` or `$PARAM` block at top
level, so `parse()` fails on them by design and nothing is wrong. `aad_mrgsolve_model.R`
and `npc_mrgsolve_model.R` are of that kind; they were checked and are fine. The
distinguishing test is whether the first non-comment line is a `$BLOCK` or R code.
`abm_mrgsolve_model.R` above is a genuine script whose embedded string is broken by
an apostrophe, which is why it is a real defect and those are not.

---

## 2. Two Graphviz maps fail ranking

`niemann-pick-disease-type-c/npc_qsp_model.dot` and
`tardive-dyskinesia/td_qsp_model.dot` both make `dot` exit non-zero with:

```
Error: trouble in init_rank
```

`dot` still writes the `.svg` and `.png`, so the rendered images in the repository
look fine, but the layout engine is reporting a ranking problem — normally a cycle
involving cluster edges.

**Confirmed upstream:** running `dot -Tsvg` directly on the two untranslated Korean
originals produces the identical error.

---

## 3. Missing font on the build machine (environment, not a defect)

Many maps set `fontname="NanumGothic"`. Where that font is not installed, `dot`
emits `Pango-WARNING: couldn't load font "NanumGothic …", falling back to "Sans"`.
This is pre-existing and affects the originals equally. The attribute is a font
name, so the translation does not change it.

Worth considering upstream: now that the translated maps carry Latin text, a Latin
font stack would render them better. That is a change to the originals, so it is
not made here.

---

## 4. Possible authoring inconsistencies noted in passing

These are content questions for the author, not build failures. Translated as
written, per the rule that a translation does not correct its source.

| File | Observation |
|---|---|
| `niemann-pick-disease-type-c/npc_references.md` | Line 294 lists T5 among the validation targets not reproduced, but the table two lines above marks T5 as a 🎯 calibration target. |
| `toxic-alcohol-poisoning/README.md` | A cross-reference points at "§4" for the bicarbonate-space effect, but that subsection sits inside §①. |
| `toxic-alcohol-poisoning/README.md` | Prose says "37.5 at 1 hour" while the corresponding table row is t = 0.5 h. |
| `visceral-leishmaniasis/README.md` | The §3 separatrix table gives CD4/TMEM thresholds, but the committed `vl_reference_output.txt` prints "none" for every one of those rows. |
| `chronic-osteomyelitis/com_mrgsolve_model.R` | Line 118 states that the Korean prose is kept on the R-comment side of the file — a claim that is self-referentially false in the translated copy, and left as written. |

---

## 5. `allergic-bronchopulmonary-aspergillosis/README.md` — broken markdown table

**Line 280.** The A13 row's header cell contains an unescaped pipe inside a code
span:

```
**`E*|fixed`** (닫힌 형태)
```

GFM does not treat a `|` inside a code span as literal in a table row, so the header
splits into 5 cells while the delimiter row and every body row have 4. The rendered
table gains a spurious empty column.

**Fix upstream would be:** escape it as `\|`, or use `E*` and `fixed` without the
pipe.

## 6. Internal count inconsistencies

Numbers that contradict each other within one file. Translated as written.

| File | Observation |
|---|---|
| `allergic-bronchopulmonary-aspergillosis/README.md` | The header says a 43-compartment mrgsolve model, the Files table says 42-compartment, and a section heading says "mrgsolve 43 · Python 42". The prose alternates between "42-state model" and the 43-compartment claim for the same object. |
| `chemotherapy-induced-nausea-vomiting/README.md` | The table row says "fitted (8)" and lists 8 parameters, and the calibration section repeats 8 against 9 anchors — but the reported fitted values list 9, including `NAP`, which the failure discussion confirms was added to the fit in a second stage. |
| `chemotherapy-induced-nausea-vomiting/README.md` | Section-name drift: the section is titled "reported failures" but is referred to elsewhere as "the failures reported", and cross-referenced as "below" from a section that sits after it. |
| `adrenocortical-carcinoma/README.md` | "8 results" heads a section containing 10 subsections (`0.` plus `A.`–`I.`). |
| `adrenocortical-carcinoma/README.md` | Drug-name typo: 템로졸로마이드 for temozolomide (템 / 테모 transposition). Translated as "temozolomide". |

---

## 7. A corrupted PubMed link

**`immune-checkpoint-inhibitor-colitis/icic_references.md:493`** contains a PubMed
URL whose id has a stray CJK character in it:

```
<https://pubmed.ncbi.nlm.nih.gov/31872征>
```

The author noticed and added a correction on the following line
(`(정정 링크: …/31874109/)`), but the broken URL itself is still there. The line has
no Hangul, so the line-level path left it byte-identical and the same broken token
is now in the translation too.

**Checked whether this is systemic: it is not.** Scanning every tracked file for a
PubMed id that begins with a digit but contains a non-digit finds exactly this one
occurrence:

```bash
git ls-files | while read -r f; do
  grep -Hn -oP 'pubmed\.ncbi\.nlm\.nih\.gov/\d+[^\d/\s)>\]"'"'"']\S*' "$f" 2>/dev/null
done
```

The many `pubmed.ncbi.nlm.nih.gov/?term=…` URLs elsewhere are **not** corruption —
several `*_references.md` files state explicitly that their links are PubMed search
queries rather than id links.

## 8. Smaller factual and transliteration slips

Translated as written.

| File | Observation |
|---|---|
| `neonatal-hyperbilirubinemia/nhb_references.md:465` | Describes Guangdong as eastern China; it is in southern China. |
| `aneurysmal-subarachnoid-hemorrhage/sah_references.md` | Uses 헵토글로빈 for haptoglobin throughout; the standard Korean transliteration is 합토글로빈. Cosmetic. |
| `adrenocortical-carcinoma/acc_references.md:444` | `13 (RCT 무의차)` looks like a contraction of 유의차 없음; read as "no significant difference". |
| `toxic-alcohol-poisoning/tap_shiny_app.R:64` | A single string is half Korean, half English: `x = "시간 since ingestion (h)"`. |
| `toxic-alcohol-poisoning/tap_shiny_app.R:386` | A `data.frame()` column name contains a space without `check.names = FALSE`, so `make.names()` renders that header with an inserted dot in the running app. |
