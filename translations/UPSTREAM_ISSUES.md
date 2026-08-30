# Defects found in the upstream files while translating

Translating a file means reading every line of it, and running the gates means
parsing and rendering it. That turns up bugs in the originals. They are recorded
here rather than fixed, because of the rule in [`README.md`](README.md): we do not
edit upstream files.

Each entry states how it was confirmed to be upstream and not something the
translation introduced. Any of these could be sent upstream as a patch or an issue.

Entry numbers may skip (e.g. no #33) when two agents append concurrently during
a batch refactor run — each re-checks the file's tail before appending to avoid
overwriting the other, at the cost of an occasional gap. The numbering is not
meaningful beyond "later than," so a gap is cosmetic, not a sign of a missing
or deleted entry.

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

## 9. Numbers that disagree with their own surrounding text

Found by reading; translated as written. These are the kind a reader would trip
over, so they are worth the author's attention more than the transliteration slips
above.

| File | Observation |
|---|---|
| `retinopathy-of-prematurity/README.md` | The AXIS 7 heading says the measured suppression is 10-fold weaker; the paragraph directly beneath it computes 12-fold excess suppression. |
| `retinopathy-of-prematurity/README.md` | AXIS 3 says one parameter (`KNV`) was fitted, but the calibration table row lists `KNV, KPLV`, and fixed-defect 7 records `KPLV` being tuned from 700 to 1300. |
| `osteosarcoma/README.md` | Section 3 announces seven results with subsections A-G, but the third mechanism of G is promoted to its own `###` heading while (i) and (ii) stay inline, so the section renders eight `###` children. |
| `osteosarcoma/README.md` | Two calibration rows present a model value as matching a literature range it falls outside: MAP 5-year EFS 0.596 against 0.54-0.59, and metastatic-stage survival 0.322 against 0.20-0.30. |
| `allergic-bronchopulmonary-aspergillosis/README.md` | See entry 6 -- the same file also disagrees with itself on compartment count. |

---

## 10. Two Shiny apps cannot run at all

Parse-checking all 840 tracked `.R` files found exactly two that genuinely fail to
parse. Both are Shiny apps, so neither will start.

**`autoimmune-hepatitis/aih_shiny_app.R:498`** — a vector opened with `c(` is closed
with `]`:

```r
Parameter = c("Bioavailability", "6-TGN t1/2", "Therapeutic range", "Toxic range", "TPMT"],
```

**Fix:** `]` -> `)`.

**`vascular-dementia/vad_shiny_app.R:59`** — `Inf` is a reserved literal in R and
cannot be used as an argument name unquoted:

```r
Inf = numeric(n_days + 1),  MG  = numeric(n_days + 1),
```

**Fix:** backtick it as `` `Inf` `` throughout, or rename the column (it holds
inflammation, so `Infl` would read better and avoid the reserved word entirely).

Re-runnable with `python3 translations/tools/check_r_parses.py`. That tool exists
because the naive test gets this badly wrong: 97 of the 840 files are pure mrgsolve
DSL, where `parse()` failing is correct and meaningless, and the block marker is
spelled `$PROB`, `[PROB]`, **and** `[ PROB ]` in different files -- so a `grep '^\$'`
heuristic reports dozens of healthy files as broken. `abm_mrgsolve_model.R`
(entry 1) is a third case again: a real script whose embedded string is broken by an
apostrophe. It is listed separately because `parse()` catches it for a different
reason.

## 11. Colour mappings that already do not mean what the label says

`ggplot2` orders factor levels alphabetically, and a `scale_*_manual(values = c(...))`
given an **unnamed** vector maps colours by that order. Several files rely on the
Korean sort order, and in two places the Korean order already inverts the intent:

| File | Observation |
|---|---|
| `chronic-hepatitis-d/hdv_shiny_app.R:566-573` | 효능 (efficacy) and 비용 (cost) are mapped to green and olive by sort order, but 비용 sorts first, so *cost* gets the green and *efficacy* the olive -- the reverse of the evident intent, since olive is the bile-acid colour elsewhere in the same file. |
| `chronic-hepatitis-d/hdv_shiny_app.R:357-359` | `FTase 억제율` sorts ahead of both drug names, so lonafarnib takes the second colour although the title reads lonafarnib first. |
| `prolactinoma/prl_shiny_app.R:452-453` | A `data.frame()` without `check.names = FALSE` has Korean legend labels, so `make.names()` renders them as `보고된.값` and `참값...기준선.` in the running app. |

The translations reproduce the existing behaviour rather than silently correcting
it, except at `hdv_shiny_app.R:423`, where no natural English term for 복합반응
sorts after "HDV" -- there the clinically correct wording was kept and two colours
exchange places. That is recorded here rather than hidden.

## 12. Citation and reference-list problems

| File | Observation |
|---|---|
| `elevated-lipoprotein-a/lpa_references.md:29` | Attributes the MESA race-cutoff paper to "the model's `FANC` (race) covariate", but `FANC` exists nowhere else in that directory -- not in the `.R`, `.dot`, or `.py`. The citation justifies a parameter the model does not have. |
| `rosacea/ros_references.md` | About 13 entries (13, 16, 40, 43-45, 49, 53, 73, 77, 85, 90, 92) have no volume or page numbers and give the journal as a guess ("a Front Immunol / J Invest Dermatol-family review", "(related study)"). Only the PMID is authoritative for those. |
| `rosacea/ros_references.md` ref 52 | `J Clin Aesthet Dermatol. 2006;5(3):317-9` -- the year and the volume/issue cannot both be right for that journal. Worth checking against PMID 20725568. |
| `controlled-ovarian-stimulation/cos_reference_check.py:16` | The docstring says five structural defects were found, while `controlled-ovarian-stimulation/README.md:106` says eight and tabulates them. |
| `controlled-ovarian-stimulation/cos_reference_check.py:16` | Points at a README section named 검증에서 드러난 결함, but the README has no such heading -- only a table header 드러난 결함. |

## 13. A sentence that was never finished

**`necrotizing-enterocolitis/nec_reference_model.py:1271`**, and therefore
`nec_reference_output.txt:59`, ends:

```
...the signature of a logarithm.  이것이
```

`이것이` is a bare subject ("this is") with no predicate; the next statement opens a
new CAVEAT block. The translation renders it literally as `This is` rather than
inventing the missing clause.

## 14. Physically odd fitted values (flagged, not asserted as wrong)

`heat-stroke/hs_references.md` section 5 gives fitted lumped conductances of
UA 23.5 for ice water at 2 degC and UA 27.7 for cold water at 14 degC -- a larger
conductance for the *warmer* bath. It is arithmetically consistent with the cooling
rates quoted, and the author separately flags the 1.55 against 1.16-1.29 ratio as an
over-prediction, so this may be intentional. It is recorded because those two values
are the reason a single-UA calculation disagrees with the meta-analysis.

## 15. The node/edge counts quoted in READMEs are mostly stale

Many per-disease READMEs advertise their map as "N nodes / M edges / K clusters".
Counting with graphviz itself -- `dot -Tplain <file>.dot`, then counting the `node`
and `edge` lines, which is what the renderer actually resolves the graph to --
**40 of the 61 such claims disagree with the file they describe.**

Most gaps are small (one to twenty), which is the signature of a count written when
the map was authored and not updated as the map was edited. A few are too large for
that:

| Directory | Claim | graphviz |
|---|---|---|
| `intrahepatic-cholestasis-of-pregnancy` | 182 edges | 274 |
| `head-and-neck-squamous-cell-carcinoma` | 319 edges | 341 |
| `alcohol-use-disorder` | 403 edges | 362 |
| `takotsubo-syndrome` | 218 nodes / 374 edges | 230 / 391 |
| `fibrodysplasia-ossificans-progressiva` | 130 nodes / 203 edges | 143 / 183 |
| `pyruvate-kinase-deficiency` | 283 edges | 265 |

The ICP figure is the one worth a second look on its own: 182 edges against 191
nodes would make the graph disconnected, and the real count is 274.

**Caveat on the method, stated because it matters.** graphviz's edge count is not
always the number of `->` tokens in the source: `a -> {b c}` expands to two edges,
and a repeated declaration collapses. So "disagrees with graphviz" means the number
in the README does not describe the rendered graph -- not automatically that the
author miscounted at the time. A regex over the source agreed with the READMEs even
less often (34% against 39%), which is why the graphviz figure is the one quoted
here.

Reproduce with:

```bash
dot -Tplain <disease>/<abbr>_qsp_model.dot | grep -c '^node '
dot -Tplain <disease>/<abbr>_qsp_model.dot | grep -c '^edge '
```

## 16. Four READMEs contain a table that does not render

`CLAUDE.md` records the 2026-07-02 incident where stray markup corrupted the root
README's gallery table, and `scripts/fix_readme_table.py` guards against a repeat --
but only for the root README. The 343 per-disease READMEs have no such guard, and
four of them currently carry a table whose rows do not agree on cell count, so the
table renders as garbage:

| File | Line | Problem |
|---|---|---|
| `allergic-bronchopulmonary-aspergillosis/README.md` | 280-288 | Header has 5 cells, the delimiter row and all 7 body rows have 4. The header cell is ``**`E*|fixed`** (closed form)`` -- an unescaped pipe inside a code span. |
| `essential-tremor/README.md` | 334 | Row has 6 cells against a 4-cell header: ``  `(1+|mu|)`  `` and `(ADAPTF - KAF*P_RAW)` contain unescaped pipes. |
| `postpartum-depression/README.md` | 231-238 | Header has 6 cells, delimiter and body rows have 4. |
| `radiation-induced-lung-injury/README.md` | 350-361, 371-378 | Two tables whose 16 body rows have no leading pipe, so they do not parse as tables at all. |

**The fix in three of the four cases is `|` -> `\|`.** A pipe inside a code span
does not protect itself.

Re-runnable with `python3 translations/tools/check_tables.py`.

**Note on the detector, because getting it wrong is easy in both directions.** GFM
splits table cells *before* inline parsing, so a pipe inside backticks does split
the cell, while `\|` does not. My first version masked code spans and counted `\|`
as a separator, which simultaneously hid the essential-tremor breakage and invented
breakage in five healthy files that legitimately use `\|` in a header
(`carbon-monoxide-poisoning`, `chronic-hepatitis-d`,
`gastrointestinal-stromal-tumor`, `hypophosphatasia`,
`progressive-supranuclear-palsy`). Those five are fine. The count is 4 of 835, not
the 7 of 418 the wrong version reported.

## 17. Korean stored as unicode escapes, and garbled inside them

Five files write Korean as unicode escapes rather than as Korean characters:

| File | Lines |
|---|---|
| `chronic-osteomyelitis/com_shiny_app.R` | 9 (227, 230-236, 626) |
| `necrotizing-enterocolitis/nec_scenario_results.json` | 11 scenario labels |
| `chronic-insomnia-disorder/ins_shiny_app.R` | 5 |
| `complex-regional-pain-syndrome/crps_shiny_app.R` | 1 |
| `urea-cycle-disorders/ucd_shiny_app.R` | 1 |

That is not a defect in itself. What is worth reporting is that **several of the
escaped strings are garbled**, in a way that looks like the escaping went wrong
rather than the author mistyping. Decoded, in `com_shiny_app.R`:

| Decodes to | Almost certainly meant |
|---|---|
| 총 세기재 | 총 세균 (total bacteria) |
| 총 세군 | 총 세균 |
| 세햵내 | 세포내 (intracellular) |
| 거리곰 | 격리골 (sequestrum) |
| 피질곰 | 피질골 (cortical bone) |
| 관산종료 | 관찰종료 (end of observation) |

Every wrong character is a plausible single-codepoint slip and they all sit inside
escape sequences, which is what points at the encoding rather than the typing. Note
that 세균 comes out two *different* wrong ways in the same file. The translation uses
the intended term in each case; the readings are listed here so the author can
confirm them.

Two more things in that file:

* line 227's `metric = c(...)` is dead -- the column is overwritten wholesale on the
  next line;
* line 277 contains `리<b>ㅁ</b>에서의 농도`, a bold tag wrapped round a single bare
  jamo. Read as "the rim", from the contrast with 심부 농도 (deep concentration).

## 18. An inhibition edge labelled with the wrong error type

**`prosthetic-joint-infection/pji_qsp_model.dot:529`**

```
SUB_PER -> DX_CULT [arrowhead=tee, label="배양 음성 위양성"]
```

위양성 is *false positive*. But the edge is an inhibition (`arrowhead=tee`) from the
persister subpopulation into the culture node, and persisters suppressing culture
growth produce a culture-negative **false negative** (위음성) -- which is also what
"배양 음성" (culture-negative) in the same label says. The label appears to have the
error type inverted.

Translated literally as "culture-negative false positive" rather than corrected,
per the rule that a translation does not fix its source.

## 19. Where a translation genuinely lost something, and the gate caught it

Recorded as evidence the gates work, not as an upstream defect.

Collapsing the `GABAAR` node label in `postpartum-depression/ppd_qsp_model.dot`
dropped the author's `GABA&#8329;` HTML numeric entity (subscript A). `check_tokens`
reported `lost ['8329']` -- the entity's digits read as a number, which is a happy
accident of the number check rather than something it was designed for. The entity
was restored, the file re-applied, re-verified and re-rendered, and the check then
passed.

## 20. `visceral-leishmaniasis/vl_qsp_model.dot` crashes graphviz entirely

`dot -Tsvg` segfaults (signal 11, no output file produced) on this map, both the
untranslated original and the translation. Preceded by eight warnings of the same
shape:

```
Warning: SPLEEN was already in a rankset, deleted from cluster VL_QSP
Warning: LIVER was already in a rankset, deleted from cluster VL_QSP
... (MARROW, SKIN, ALIP, MILC, PMC, SB5C)
```

Each of those eight nodes is declared inside more than one `subgraph cluster_*`
block, which graphviz tolerates as a warning most of the time but appears not to
survive here -- the crash follows immediately after the eighth warning.

**Confirmed upstream:** identical warnings and the identical segfault on the
untranslated original, run directly with `dot -Tsvg <original> -o /dev/null`.
Unlike every other map in this library (which write an image despite an
`init_rank` warning), this one produces **no image at all**, before or after
translation.

**Fix upstream would be:** remove each of the eight nodes from every
`subgraph cluster_*` block but one.

## 21. `methemoglobinemia/README.md` claims five structures, table has six

Line 33 says "넣은 것은 다섯 개의 구조뿐이기 때문입니다" (only five structures were
put in), but the table immediately below it (lines 35-40) has six rows, each
pairing one structural assumption with one computed result. Translated literally
as "five" per the no-correction rule.

## 22. `prosthetic-joint-infection/README.md` heading says three, lists four

Line: "## 2. 세 개의 산수 (The three arithmetic pillars)" -- but the section
contains four numbered pillars (1-4). Translated literally as "The Three
Arithmetic Pillars" per the no-correction rule.

## 23. `endometriosis/README.md` says the map/model files are "planned," but they exist

The README text states the mechanistic-map (.dot/.svg/.png) and mrgsolve model are
"추후 추가 예정" (planned to be added later), but `endo_qsp_model.dot/.svg/.png`
and `endo_mrgsolve_model.R` already exist in the directory. The README is stale
relative to the actual file contents. Translated literally per the no-correction
rule.

## 24. `neurofibromatosis-type-1/README.md` glossary mistranslates 신경섬유종

A Korean/English glossary row pairs "신경섬유종 (뉴로피브로민)" with "Neurofibromin,"
but 신경섬유종 means "neurofibroma" (the tumor), not "neurofibromin" (the protein
whose loss causes it) -- an author slip in the original, not a translation
artifact. Recorded here since the glossary's Korean column was dropped from the
translation (see per-file note: a Korean/English quick-reference table cannot
survive as a table in an all-English document, so it was converted to a
single-column English term list, same terms, per the zero-residual-Hangul
requirement).

## 25. `chronic-pyelonephritis/README.md` Shiny tab says five, scenario table has seven

Tab (5)'s label says "5가지 항생제 전략 동시 비교" (comparing 5 antibiotic
strategies simultaneously), but the scenario table above it lists 7 scenarios.
Translated literally as "five" per the no-correction rule.

## 26. `retinitis-pigmentosa/rp_qsp_model.dot` also crashes graphviz, same class as #20

`dot -Tsvg retinitis-pigmentosa/rp_qsp_model.dot` produces the same "already in a
rankset, deleted from cluster" warning storm as issue #20
(`visceral-leishmaniasis`) -- about 20 nodes (`G1`-`G14`, `CL1`-`CL7`) declared
inside more than one `subgraph cluster_*` block -- and, on this machine's
Graphviz build, produces **no image at all** (exit code 2), for both the
original and the translated `.dot`.

Unlike `visceral-leishmaniasis`, the repository already has a committed
`rp_qsp_model.svg`/`.png` (dated before this session), so whoever built the
model originally had a Graphviz version that tolerates the warning and still
writes an image, as `render_maps.py`'s own docstring notes most maps do. This
machine's Graphviz does not. The translated `.dot` never got its own rendered
`rp_qsp_model_en.svg`/`.png` as a result -- confirmed via
`translations/tools/render_maps.py`, which reports it as a hard `FAIL`, not the
usual tolerated `warn`.

**Fix upstream would be:** the same as #20 -- remove each duplicated node from
every `subgraph cluster_*` block but one.

---

## 27. `kidney-transplant-rejection/ktx_mrgsolve_model.R` does not compile under mrgsolve 2.0.1

`$ODE` reuses the for-loop iterator name `k` in four separate, unrelated
top-level loops: the belatacept dosing sum (`for(int k = 0; k < 70; k++)
belasum += ...`), the anti-IL-6 (tocilizumab) dosing sum (`for(int k = 0; k <
(int)TCZN; k++) tczsum += ...`), the anti-CD38 dosing sum (two loops, `k < 5`
and `k < 4`), and the plasma-exchange sum (`for(int k = 0; k < (int)PLEXN;
k++) plexsum += ...`). Each loop's `k` is properly scoped to its own `for`
statement in the source, but mrgsolve's variable-hoisting preprocessor lifts
every `int k` declaration it finds to the top of the generated C++ function
without deduplicating, so the generated source declares `int k` five times in
one scope. GCC (tested: 14.3.0, the compiler mrgsolve invokes on this machine)
rejects that as `redefinition of 'int {anonymous}::k'` and the build fails
before any simulation can run.

**Confirmed upstream:** reproduced by `mread()`/`mcode()` on the untouched
original file alone, with no changes and no translation involved — the
generated `ktxdiag-mread-source.cpp` shows five consecutive `int k;` fields
in the hoisted-variable struct (lines 75, 79, 82, 83, 88 in the build attempted
here), one per loop, matching the compiler's error line numbers exactly.

**Why this matters here:** this file was the subject of a PK/PD refactor
(`ktx_mrgsolve_model_refactored.R`, tocilizumab's PK block only) whose mandatory
verification step requires actually building and running both the original and
the refactored model. Building either one, unmodified, fails with the error
above — it is unrelated to the refactor and reproduces identically from the
refactored file too, since that file copies these loops verbatim. Verification
was carried out against **in-memory-only** copies of both models' code with the
five reused `k` declarations mechanically renamed to be unique (`kbe`, `ktc`,
`kf1`, `kf2`, `kpl`) — a loop-variable rename only, no change to any loop bound,
increment, or body — applied identically to both models purely to make them
buildable for comparison. Neither the checked-in original nor the delivered
`_refactored.R` sibling was changed to work around this; both still contain the
original's `int k` loops as written, and so **neither currently builds against
mrgsolve 2.0.1** without the same workaround.

**Fix upstream would be:** rename each loop's `k` to a distinct name (or wrap
each loop body in its own `{ }` block, which does not appear to be enough by
itself since the hoisting pass looks for `int k` textually rather than
respecting brace scope — the safest fix is distinct names throughout `$ODE`).

---

## 28. `thyroid-eye-disease/ted_mrgsolve_model.R` has two independent mrgsolve-2.0.1 build breaks

Found while verifying a tocilizumab-only PK/PD refactor
(`ted_mrgsolve_model_refactored.R`); both reproduce identically from the
untouched original with no changes involved, and both are unrelated to
tocilizumab.

**(a) `[PARAM]` annotation missing its description field.** Two annotated
parameter lines have only `NAME : VALUE`, no third (description) field:

```
HILL_HEAR        : 1.8
HILL_CUSH        : 2.0
```

mrgsolve's annotated-block parser rejects this outright: `Error: improper
annotation format`. The model does not build at all until a description is
added to both lines.

**(b) `TIME` is not usable inside `$ODE` in this mrgsolve version.** The
original uses `TIME` twice inside `[ODE]` — `t_month = TIME/730.0 +
DURATION_MO;` and, inside `smoke_mult_now`, `exp(-TIME/4380.0)`. In mrgsolve
2.0.1, `TIME` expands to `self.time` (see `mrgsolve/base/modelheader.h`),
but `MRGSOLVE_ODE_SIGNATURE` (`mrgsolve/base/mrgsolv.h`) — the generated
`$ODE` function's parameter list — does not include a `self` parameter;
only the `$MAIN`/`$TABLE`/`$EVENT`/`$PREAMBLE` signatures do. Confirmed via
the generated `.cpp`: `error: 'self' was not declared in this scope`,
pointing exactly at the `TIME` macro expansion on the `t_month` line. The
correct in-`$ODE` accessor in this version is `SOLVERTIME` (`_ODETIME_[0]`).

**Confirmed upstream:** both reproduce from `mread()` on the untouched
original alone (via a temporary in-memory copy — see below), no translation
or refactor content involved.

**Why this matters here:** the refactor's mandatory verification step
requires actually building and running both the original and the
`_refactored.R` sibling. Neither builds, unmodified, past bug (a); fixing
(a) alone then hits bug (b). Verification was carried out against
**scratch, in-memory-only copies** of both models with (a) a description
string added to the two malformed `HILL_*` lines and (b) `TIME`→`SOLVERTIME`
substituted at the two `$ODE` usages — applied identically to both models,
changing no parameter value and no PK/PD equation, purely to make them
buildable for comparison. Neither the checked-in original nor the delivered
`_refactored.R` was changed to work around this; both still contain the
original's `TIME`-in-`$ODE` usage and malformed annotations exactly as
written, and so **neither currently builds against mrgsolve 2.0.1** (locally
via `mread()`/`mcode()`, or through the qspserver `mrgsolve_api` REST
service, which hit the identical two errors) without the same workaround.

**Fix upstream would be:** add a description field to the `HILL_HEAR`/
`HILL_CUSH` lines, and replace both `TIME` usages inside `[ODE]` with
`SOLVERTIME`.

---

## 29. `polymyalgia-rheumatica/pmr_mrgsolve_model.R` uses the deprecated `_init_<CMT>` idiom

Found while verifying a tocilizumab-only PK/PD refactor
(`pmr_mrgsolve_model_refactored.R`); reproduces identically from the
untouched original, unrelated to tocilizumab.

`$MAIN` sets initial conditions with the old `_init_<CMT> = value;` form
(e.g. `_init_CORT = ...`) for `CORT`, `IL6`, `SIL6R`, `CRP`, `ESR`, `BMD`,
`PMRAS`, `FLARE`. Local mrgsolve 2.0.1 (Windows, `R-4.6.0`) no longer
accepts this — `error: '_init_CORT' was not declared in this scope` —
confirmed with a minimal 1-compartment reproduction outside this file:
`_init_X = A0;` fails the same way, while the modern `X_0 = A0;` idiom
compiles and runs.

**Confirmed upstream:** reproduces from the untouched original alone, no
translation or refactor content involved.

**Note — this is a local-toolchain-specific break, not a universal one.**
Unlike issues #27 and #28, this file compiles and runs *unmodified* through
the qspserver `mrgsolve_api` container (Debian/GCC toolchain) — its
`/model_manifest` call on the original's byte-identical extracted DSL
succeeded with no patch applied. The `_init_<CMT>` idiom is accepted there.
So the "does this original still run" answer depends on which mrgsolve
build/toolchain is asked, not just the mrgsolve version.

**Why this matters here:** the refactor's mandatory verification step could
not use local mrgsolve for this reason, so verification was instead run
against the qspserver API, which needed no workaround for this issue.
Neither the checked-in original nor the delivered `_refactored.R` sibling
was changed — both still use `_init_<CMT>` exactly as written.

**Fix upstream would be:** replace each `_init_<CMT> = value;` line with the
modern `<CMT>_0 = value;` form.

---

## 30. `sepsis/sep_mrgsolve_model.R` has three independent mrgsolve-2.0.1 build breaks

Found while verifying a tocilizumab-only PK/PD refactor
(`sep_mrgsolve_model_refactored.R`); all three reproduce identically from
the untouched original, unrelated to tocilizumab.

1. **`$CAPTURE` lists compartment names directly** (e.g. `BACT`, `TNF`,
   `TOCI_C`, …) — this mrgsolve version rejects that outright: *"compartment
   should not be in $CAPTURE"*.
2. **Every compartment is declared via both `$CMT` (bare name) and `$INIT`
   (`name = value`)** — this mrgsolve version's `validObject()` flags that
   as *"Duplicated model names"*, even though declaring a compartment both
   ways used to be an unremarkable, common mrgsolve pattern.
3. **`$PARAM` names `IL6_0`, `IL10_0`, `PAI1_0` collide** at the generated-
   C++ stage with mrgsolve's auto-generated `<compartment>_0` initial-value
   symbols for the identically-named compartments `IL6`, `IL10`, `PAI1`,
   producing `conflicting declaration`/`redeclaration` compiler errors.

**Confirmed upstream:** reproduces from the untouched original alone, no
translation or refactor content involved; also reproduces identically
through the qspserver `mrgsolve_api` container, so this is a genuine
version-level incompatibility, not a local-toolchain quirk like issue #29.

**Why this matters here:** the refactor's mandatory verification step
requires building and running both the original and the `_refactored.R`
sibling; neither builds unmodified. Verification was carried out against
in-memory-only request payloads (never saved to either `.R` file) with:
compartment names dropped from `$CAPTURE`, the `$CMT` block removed in
favor of `$INIT` alone, and the three colliding `$PARAM` names renamed
`IL6_SS0`/`IL10_SS0`/`PAI1_SS0` (all usages updated to match) — applied
identically to both models, changing no numeric value. Neither the
checked-in original nor the delivered `_refactored.R` contains this
workaround; both still use the original's `$CAPTURE`/`$CMT`+`$INIT`/param-
naming exactly as written, and so **neither currently builds against
mrgsolve 2.0.1** without the same workaround.

**Fix upstream would be:** remove compartment names from `$CAPTURE`
(mrgsolve captures compartment amounts automatically), declare each
compartment via `$INIT` only (drop the redundant `$CMT` block), and rename
the three colliding `$PARAM` entries away from the `<compartment>_0`
pattern.

---

## 31. `rheumatoid-arthritis/ra_mrgsolve_model.R` does not compile under mrgsolve 2.0.1, and its `$MAIN`-scoped PD drivers report one interval late

Found while verifying a tocilizumab-only PK/PD refactor
(`ra_mrgsolve_model_refactored.R`); both issues reproduce identically from
the untouched original, unrelated to tocilizumab's own math.

**1. `$CAPTURE` duplicates compartment names.** `CRP`, `INFLAM`, `VDH`,
`CART`, `STAT3_P`, and `FLS_ACT` are `$CMT` compartments *and* are listed in
`$CAPTURE`. mrgsolve 2.0.1 refuses to build the model:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: CRP,INFLAM,VDH,CART,STAT3_P,FLS_ACT
```

**Confirmed upstream:** reproduces from the untouched original alone (no
translation/refactor content involved), independently on two separate
mrgsolve 2.0.1 installs — a local R install and the qspserver
`mrgsolve_api` container (Debian/GCC toolchain) — so unlike issue #29 this
is a genuine version-level break, not a local-toolchain quirk.

**2. Every disease-driving PK/PD variable is computed in `$MAIN`, not
`$ODE`/`$TABLE`.** `C_TCZ`, `C_ADA`, `C_MTX`, `C_BARI`, `OCC_TCZ`,
`JAK_INH`, `IL6_sig`, `TNF_INH_ADA`, `DRUG_ANTI_INFLAM`, `MTX_EFF_INFLAM`
are all `double`s declared in `$MAIN`. mrgsolve evaluates `$MAIN` once per
requested output/dosing interval, using the state as of the *start* of that
interval — so every one of these reported values is one interval stale
relative to the true continuous state. This is easy to miss at a fine
output grid but is large at the model's own `delta=24` (daily) scenarios:
after an IV tocilizumab bolus, `OCC_TCZ` reads exactly `0` in the very first
daily row even though the receptor is already ~96% occupied by then (a
finer, `delta=0.05` probe of the identical dose/model shows occupancy
already at 0.96 within the first hour); occupancy only "catches up" one row
later. Confirmed by comparing the same scenario at `delta=1` vs `delta=0.01`
against the underlying compartment states, and by an isolated version of
the same ODE with the double declarations moved to `$ODE`+`$TABLE`, which
tracks the true state with no lag.

**Why this matters here:** the refactor's mandatory verification step
requires building both the original and `_refactored.R`; the original does
not build unmodified (issue 1). Verification was carried out via the
qspserver `mrgsolve_api` service on an in-memory `$CAPTURE`-deduplicated
copy of both DSLs (compartment names removed from `$CAPTURE`; nothing else
changed — compartment values are always present in mrgsolve output
regardless of `$CAPTURE`, so this changes nothing about what is reported).
That same minimal, disclosed fix (and only that fix) was also applied to
the delivered `ra_mrgsolve_model_refactored.R`, since a delivered refactor
that cannot itself build would defeat the point of verifying it — see
`ra_refactor_notes.md`. `ra_mrgsolve_model.R` itself is untouched. Issue 2
(the `$MAIN` lag) was **not** fixed anywhere, including in the refactored
sibling: it is a pre-existing behavioral quirk of the whole model (all four
compounds, not just tocilizumab), and preserving it exactly is what let
the refactor's verification match the original to 0.0 absolute difference
rather than merely "close."

**Fix upstream would be:** remove `CRP`/`INFLAM`/`VDH`/`CART`/`STAT3_P`/
`FLS_ACT` from `$CAPTURE`; move the `$MAIN`-scoped double declarations
listed above into `$ODE` (so they update every solver evaluation) and/or
recompute them in `$TABLE` (so reporting reflects the state at the actual
output time) rather than `$MAIN`.

## 32. `myopia-progression/myp_mrgsolve_model.R` `$PARAM` block has a malformed annotation continuation line

Found while verifying a 7-methylxanthine ("MX")-only PK/PD refactor
(`myp_mrgsolve_model_refactored.R`); unrelated to MX's own math and
reproduces identically from the untouched original.

The `TRTPD` parameter's annotated description is split across two lines:

```
TRTPD    : 0.15  : imposed peripheral defocus from spectacles/CL (D; -ve = myopic)
                 : SV +0.15, PAL -0.30, MiSight -0.90, DIMS -1.20, HAL -2.20
```

mrgsolve 2.0.1's annotated-parameter-block parser requires every line to
have a `name : value : description` triplet; the second line has only a
leading `:` and a description with no name/value fields, which the parser
rejects outright:

```
Error: improper annotation format
 input: : SV +0.15, PAL -0.30, MiSight -0.90, DIMS -1.20, HAL -2.20
 context: parse annotated parameter block (PARAM)
Execution halted
```

**Confirmed upstream:** reproduces from the untouched original alone (no
refactor content involved), via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — the model does not build at all, for either the
original or the refactored file, until this is fixed.

**Verification workaround (in-memory only, not committed to either
file):** the two lines were merged into a single annotated line, moving the
per-device value list into the existing description's parentheses:

```
TRTPD    : 0.15  : imposed peripheral defocus from spectacles/CL (D; -ve = myopic; SV +0.15, PAL -0.30, MiSight -0.90, DIMS -1.20, HAL -2.20)
```

This is a pure text/formatting merge — no name, value, or number was
touched. It was applied identically to scratch copies of both
`myp_mrgsolve_model.R` and `myp_mrgsolve_model_refactored.R` purely so both
would build for the `/model_manifest` and `/run_simulation` comparison;
neither the tracked original nor the delivered `_refactored.R` was changed,
and both still contain the original's two-line form as written.

**Fix upstream would be:** either merge the continuation into `TRTPD`'s own
description field (as done for verification above), or give the second line
its own `name : value : description` triplet if it was intended as a
separate parameter.

## 34. `chronic-hypothyroidism/hypo_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT` and `$INIT` jointly redeclare every compartment, and `$CAPTURE` duplicates nine of them

Found while verifying a levothyroxine ("LT4")-only PK/effect refactor
(`hypo_mrgsolve_model_refactored.R`); both defects reproduce identically
from the untouched original and are unrelated to LT4's own math (LT3/
liothyronine, the file's other compound, is affected identically and was
not touched by the refactor).

**1. `$CMT @annotated` and `$INIT` both declare the same 15 compartments.**
The file uses `$CMT @annotated` to name and describe every compartment
(`A_TRH`, `A_TSH`, ... `Eff_BMD`), then immediately follows with a separate
`$INIT` block that assigns each of the same 15 names a starting value
(`A_TRH = 5.0`, etc.). mrgsolve 2.0.1 treats `$INIT` as its own
compartment-declaring block (an alternative to `$CMT`, not a companion to
it), so using both for the same names redeclares every compartment twice:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: 1: Duplicated model names: A_TRH A_TSH A_TT4 A_TT3 A_rT3 A_LT4_gut A_LT4_c A_LT4_p A_LT3_gut A_LT3_c Eff_HR Eff_LDL Eff_BMR Eff_Sym Eff_BMD
```

**2. `$CAPTURE` repeats nine compartment names.** `A_TSH`, `A_TT4`, `A_TT3`,
`A_rT3`, `Eff_HR`, `Eff_LDL`, `Eff_BMR`, `Eff_Sym`, and `Eff_BMD` are `$CMT`
compartments *and* are listed in `$CAPTURE` — the same defect shape as
issue #31 (`ra_mrgsolve_model.R`), independently present here:

```
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE: A_TSH,A_TT4,A_TT3,A_rT3,Eff_HR,Eff_LDL,Eff_BMR,Eff_Sym,Eff_BMD
```

**Confirmed upstream:** both errors reproduce from the untouched original
alone, via the qspserver `mrgsolve_api` container (`POST /model_manifest`)
— the model does not build at all, for either the original or the
refactored file, until both are fixed.

**Verification workaround (in-memory only, not committed to either
file):** (a) the `$INIT` block was deleted and its 15 assignments moved
into `$MAIN` using the modern `<CMT>_0 = value;` idiom (e.g.
`A_TRH_0 = 5.0;`), which declares no new compartment and changes no
numeric value; (b) the nine compartment names were removed from
`$CAPTURE` (compartment states are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported).
Both patches were applied identically to scratch copies of
`hypo_mrgsolve_model.R` and `hypo_mrgsolve_model_refactored.R` purely so
both would build for the `/model_manifest` and `/run_simulation`
comparison; neither the tracked original nor the delivered
`_refactored.R` was changed — both still contain the `$CMT`+`$INIT`
duplication and the nine-name `$CAPTURE` overlap exactly as written (under
the refactor's renamed `GUT_LT4`/`CENT_LT4`/`PERI_LT4` for the three LT4
compartments, since the refactor only renamed identifiers, not the defect
pattern itself). See `chronic-hypothyroidism/hypo_refactor_notes.md`.

**Fix upstream would be:** remove the `$INIT` block and set nonzero
initial conditions via `<CMT>_0 = value;` in `$MAIN` instead (or drop
`$CMT`'s annotations and rely on `$INIT` alone, whichever the author
prefers); remove `A_TSH`/`A_TT4`/`A_TT3`/`A_rT3`/`Eff_HR`/`Eff_LDL`/
`Eff_BMR`/`Eff_Sym`/`Eff_BMD` from `$CAPTURE`.

---

## 35. `dengue/denv_mrgsolve_model.R` — nine `$ODE`-local doubles collide with same-named `$TABLE` captures under mrgsolve 2.0.1

Found while verifying an antiviral (AV)-only PK/PD refactor
(`denv_mrgsolve_model_refactored.R`); unrelated to AV's own math and
reproduces identically from the untouched original.

`$ODE` declares nine local variables (`SV`, `CO`, `MAP`, `PP`, `Jv`, `Pser`,
`Jser`, `Dser`, `Hct`) used only to compute that section's own `dxdt_*`
expressions. `$TABLE` independently recomputes the same physiology from the
`$CMT` state under `_o`-suffixed names (`SV_o`, `CO_o`, `MAP_o`, `PP_o`,
`Jv_o`, `Pser_o`, `Jser_o`, `Dser_o`, `Hct_o`) and reports them via
`capture SV = SV_o;` etc. — i.e. the author already used distinct names
specifically to avoid a clash with the `$ODE` locals. Under mrgsolve 2.0.1,
that clash happens anyway: the compiler auto-promotes every bare
`double NAME = expr;` assignment found in `$ODE` to a reportable class
member (the same mechanism that lets `$ODE`-local doubles appear in output
without an explicit `$CAPTURE`), so `$TABLE`'s `capture SV = SV_o;` then
tries to redeclare a member that `$ODE`'s own `double SV = ...;` already
created:

```
130:11: error: redefinition of 'capture {anonymous}::Jv'
  130 |   capture Jv;
      |           ^~
66:10: note: 'double {anonymous}::Jv' previously declared here
   66 |   double Jv;
      |          ^~
```

(and identically for `Jser`, `Dser`, `Pser`, `Hct`, `SV`, `CO`, `MAP`, `PP`.)

**Confirmed upstream:** reproduces from the untouched original alone (no
refactor content involved), via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — the model does not build at all, for either the
original or the refactored file, until this is worked around.

**Verification workaround (in-memory only, not committed to either
file):** within the `$ODE` block only (never touching `$TABLE`, `$MAIN`, or
any numeric value), the nine colliding local names were suffixed
`_ode` (`SV`→`SV_ode`, `CO`→`CO_ode`, `MAP`→`MAP_ode`, `PP`→`PP_ode`,
`Jv`→`Jv_ode`, `Pser`→`Pser_ode`, `Jser`→`Jser_ode`, `Dser`→`Dser_ode`,
`Hct`→`Hct_ode`), with every downstream in-`$ODE` reference to each
(e.g. `PP` inside `Pc = (MAP + rratio*PVcvp)/(1+rratio)`, `CO`/`Hct` inside
the oxygen-delivery block) updated to match. `$TABLE`'s own `_o`-suffixed
recomputation and its `capture NAME = NAME_o;` lines were left completely
untouched, since they never referenced the `$ODE`-local bare names to begin
with. Applied identically to scratch copies of both
`denv_mrgsolve_model.R` and `denv_mrgsolve_model_refactored.R` purely so
both would build for the `/model_manifest` and `/run_simulation`
comparison; neither the tracked original nor the delivered `_refactored.R`
contains this workaround, and both still use the original's bare `$ODE`
local names exactly as written — so **neither currently builds against
mrgsolve 2.0.1** without the same rename.

**Fix upstream would be:** rename the nine `$ODE`-local doubles away from
the `$TABLE` capture names they collide with (e.g. adopt the same `_ode`
suffix used for verification above, or the reverse — rename `$TABLE`'s
`_o` locals to the bare names and drop `$ODE`'s bare locals), so mrgsolve's
auto-promoted `$ODE` members and `$TABLE`'s explicit `capture` declarations
never share a name.

---

## 36. `age-related-macular-degeneration/amd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` redeclare every compartment, and `$ODE` writes directly to two read-only compartment states

Found while verifying the anti-VEGF ("VIT", the file's only modeled
compound) PK/effect-interface refactor (`amd_mrgsolve_model_refactored.R`);
both defects reproduce identically from the untouched original.

**1. `$CMT` (bare, unannotated) and `$INIT` jointly redeclare all 20
compartments.** Same defect class as issue #34
(`chronic-hypothyroidism/hypo_mrgsolve_model.R`), independently present
here with plain (non-`@annotated`) `$CMT`: mrgsolve 2.0.1 treats `$INIT`'s
`NAME = value` list as its own compartment declaration rather than a
companion to `$CMT`, so every name is declared twice:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: DRUG_VIT DRUG_RET DRUG_SYS VEGF_FREE VEGF_BOUND VEGFR2_ACT ANG2_FREE ANG2_BOUND C3_LOCAL C5_LOCAL MAC_LOCAL RPE_NORM RPE_DAM LIPOFUSCIN DRUSEN CNV_AREA FLUID_EX GA_AREA BCVA_SCORE PR_FRAC
```

**2. `$ODE` assigns directly to two compartment states to clamp them at
zero** — `if(FLUID_EX < 0) FLUID_EX = 0;` and `if(PR_FRAC < 0) PR_FRAC =
0;`. mrgsolve 2.0.1 passes every `$CMT` state into the generated `$ODE`
function as a `const double&`, so writing to it is a hard C++ compile
error, not merely a warning:

```
439:27: error: assignment of read-only reference 'FLUID_EX'
  439 | if(FLUID_EX < 0) FLUID_EX = 0;
      |                  ~~~~~~~~~^~~
446:25: error: assignment of read-only reference 'PR_FRAC'
  446 | if(PR_FRAC < 0) PR_FRAC = 0;
      |                 ~~~~~~~~^~~
```

This is a defect class not previously logged in this file (distinct from
issues #31/#34's `$CAPTURE`-duplicates-a-compartment shape). It also looks
to have always been a no-op even where it once compiled: only `dxdt_*`
feeds mrgsolve's integrator, so reassigning the state variable itself
inside `$ODE` cannot affect the next solver step or the reported
trajectory either way — the clamp the author evidently intended
(preventing `FLUID_EX`/`PR_FRAC` from reporting as negative) was never
actually enforced by this line, compilable or not.

**Confirmed upstream:** both errors reproduce from the untouched original
alone, via the qspserver `mrgsolve_api` container (`POST
/model_manifest`) — the model does not build at all, for either the
original or the refactored file, until both are worked around.

**Verification workaround (in-memory only, not committed to either
file):** (a) the `$INIT` block was deleted and its 20 assignments moved
into a new `$MAIN` block using the modern `<CMT>_0 = value;` idiom (e.g.
`CENT_VIT_0 = 0;`), declaring no new compartment and changing no numeric
value; (b) the two illegal `if(... < 0) ... = 0;` lines were deleted
outright (see the no-op reasoning above — removing a statement that could
never have affected the integrated trajectory changes nothing numeric).
Both patches were applied identically to scratch copies of
`amd_mrgsolve_model.R` and `amd_mrgsolve_model_refactored.R` purely so
both would build for the `/model_manifest` and `/run_simulation`
comparison; neither the tracked original nor the delivered
`_refactored.R` was changed — both still contain the `$CMT`+`$INIT`
duplication and the two illegal compartment writes exactly as written
(under the refactor's renamed `CENT_VIT`/`PERI_VIT`/`SYS_VIT`, since the
refactor only renamed identifiers, not the defect pattern itself). See
`age-related-macular-degeneration/amd_refactor_notes.md`.

**Fix upstream would be:** remove the `$INIT` block and set initial
conditions via `<CMT>_0 = value;` in `$MAIN` instead; delete the two
`if(FLUID_EX < 0) FLUID_EX = 0;` / `if(PR_FRAC < 0) PR_FRAC = 0;` lines (or,
if the negative-value guard is actually wanted, implement it properly by
zeroing the *outflow term* in the affected `dxdt_` expression when the
state is at/near zero, e.g. `if(FLUID_EX <= 0 && Fluid_inflow <
Fluid_outflow) dxdt_FLUID_EX = 0;`).

## 37. `diabetic-ketoacidosis/dka_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: a multi-variable `double` declaration loses its type on all but the first name

Found while verifying the insulin PK/effect-interface refactor
(`dka_mrgsolve_model_refactored.R`); the defect reproduces identically from
the untouched original and has nothing to do with insulin.

**The pH bisection solve in `$ODE` (acid-base block) declares four working
variables on one line, only the first of which keeps its type after
mrgsolve's preprocessing:**

```c
double lo = 5.60, hi = 8.30, pHm, hres;
```

mrgsolve 2.0.1's `$ODE` preprocessor appears to hoist/rewrite `double`
declarations line by line rather than per comma-separated declarator, so the
generated C++ keeps `double lo = 5.60;` but drops the `double` from `hi`,
`pHm`, and `hres`, leaving them referenced with no declaration at all:

```
767:12: error: 'hi' was not declared in this scope
  767 | lo = 5.60, hi = 8.30, pHm, hres;
      |            ^~
767:23: error: 'pHm' was not declared in this scope
  767 | lo = 5.60, hi = 8.30, pHm, hres;
      |                       ^~~
767:28: error: 'hres' was not declared in this scope
  767 | lo = 5.60, hi = 8.30, pHm, hres;
      |                            ^~~~
```

Since `pH`/`HCO3`/`AG` and every downstream quantity in the file (the entire
acid-base, ketone, and renal blocks) are computed from this bisection loop's
result, the model cannot be built at all — for either the original or the
refactored file — until this is worked around. This is a new defect class
for this file (not the loop-counter-redeclaration class of issue #27, nor
the `$CMT`+`$INIT`/`$CAPTURE`-duplicate classes of issues #31/#34-36):
mrgsolve's handling of a single multi-declarator `double` statement itself,
independent of which names are involved.

**Confirmed upstream:** reproduces from the untouched original alone via the
qspserver `mrgsolve_api` container (`POST /model_manifest`).

**Verification workaround (in-memory only, not committed to either file):**
the one line above was split into four separate declarations —

```c
double lo = 5.60;
double hi = 8.30;
double pHm;
double hres;
```

— applied identically to scratch copies of `dka_mrgsolve_model.R` and
`dka_mrgsolve_model_refactored.R` purely so both would build for the
`/model_manifest` and `/run_simulation` comparison. Neither the tracked
original nor the delivered `_refactored.R` was changed; both still contain
the single-line multi-declaration exactly as written. See
`diabetic-ketoacidosis/dka_refactor_notes.md`.

**Fix upstream would be:** split the one `double lo = 5.60, hi = 8.30, pHm,
hres;` line into four separate `double` statements (or confirm, on a newer
mrgsolve, whether the multi-declarator form is actually supported and this
is version-specific).

---

## 38. `distal-renal-tubular-acidosis/drta_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: a valueless `#define` crashes the parameter-block parser

Found while verifying a thiazide (HCTZ) PK/effect-interface refactor
(`drta_mrgsolve_model_refactored.R`); the defect is in `$GLOBAL`, wholly
unrelated to HCTZ, and reproduces identically from the untouched original.

`$GLOBAL`'s first line is

```c
#define _MRG_DRTA_
```

a bare include-guard-style macro with no replacement value, and it is never
referenced anywhere else in the file (`grep -n "_MRG_DRTA_"` matches only
this one line) — apparently vestigial. Under mrgsolve 2.0.1, a `#define`
with no value token crashes the model's parameter-definition parser before
compilation is even attempted:

```
Error in FUN(X[[i]], ...) : subscript out of bounds
Calls: simcore_load_model ... pp_defs -> s_pick -> nonull -> unlist -> sapply -> lapply
Execution halted
```

Bisected by binary-searching the file's `$GLOBAL` block content down to this
one line: removing only this line (leaving the rest of `$GLOBAL`, and every
other block, untouched) makes the model build; restoring it alone reproduces
the crash. The error is unrelated to `$PARAM`'s own content — a trimmed
model keeping the real `$PARAM`/`$CMT` blocks but a one-line placeholder
`$GLOBAL` compiles fine, and the crash reappears the instant the bare
`#define` line is added back, with or without any other `$GLOBAL` content
present.

**Confirmed upstream:** reproduces from the untouched original alone (no
refactor content involved), via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — the model does not build at all, for either the
original or the refactored file, until this is worked around.

**Verification workaround (in-memory only, not committed to either file):**
the single line `#define _MRG_DRTA_` was deleted from scratch copies of both
`drta_mrgsolve_model.R` and `drta_mrgsolve_model_refactored.R` (plus, since
neither file's DSL block is wrapped in `code <- '...'`/`mcode()` — both are
raw mrgsolve block-marker text meant for direct `mread()`, like the
`thyroid-eye-disease` file in issue #28 — trimming each scratch copy to end
before its trailing `$ENV` block of R-only scenario-driver code, which is
not part of the DSL `/model_manifest`/`/run_simulation` compiles). No
numeric value or equation was touched by either change. Neither the tracked
original nor the delivered `_refactored.R` contains this workaround, and
both still contain the bare `#define` exactly as written — so **neither
currently builds against mrgsolve 2.0.1** without removing that one line.
See `distal-renal-tubular-acidosis/drta_refactor_notes.md`.

**Also noted during the same verification (not a defect in the checked-in
file, an engine-configuration mismatch):** every one of this model's own 28
scenarios is run, in its own R driver (`sim_scenario()`), with
`mrgsim(..., maxsteps = 5e6, hmax = 0.5)` — 250x mrgsolve's default
`maxsteps` of 20000, plus an explicit step-size cap. The qspserver
`mrgsolve_api` `/run_simulation` endpoint exposes no way to set either
option. Reproduced on the untouched, `#define`-patched original alone: an
undosed run exceeds mrgsolve's default 20000-step budget somewhere between
24h and 48h of simulated time, and a run carrying even a single dosing
event (any compartment) can exceed it within hours, both failing with
`[lsoda] 20000 steps taken before reaching tout` / `excess work done on
this call`. This is consistent with the file's own diurnal meal-forcing
term (`neap_rate()`, three raised-cosine humps per day) and its per-step
warm-started bisection (`solve_urine_pH()`) making the system numerically
demanding enough that the author's own driver anticipated needing a much
larger step budget than mrgsolve's default. Verification of the HCTZ
refactor below was therefore run over a 12h single-dose window (comfortably
under the failure threshold found by bisection), not the full 150-365 day
duration any of the model's own named scenarios actually use — see the
"Verification" section of `drta_refactor_notes.md` for the exact window and
why.

**Fix upstream would be:** delete the unused `#define _MRG_DRTA_` line (or
give it a value, e.g. `#define _MRG_DRTA_ 1`, if some future code is meant
to guard on it).

## 39. `acute-intermittent-porphyria/aip_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two independent, unrelated defects

Found while verifying a givosiran (GIV) + hemin (HEM) PK/effect-interface
refactor (`aip_mrgsolve_model_refactored.R`); both defects reproduce
identically from the untouched original via the qspserver `mrgsolve_api`
container and are unrelated to either compound's PK — they block the model
from building at all, refactored or not.

1. **`$PARAM @annotated` header, but the body uses plain (non-annotated)
   syntax.** Every line is `NAME = value // comment`, mrgsolve's ordinary
   `$PARAM` idiom — but the block is declared `$PARAM @annotated`, which
   requires the colon-delimited `NAME : description : (unit)` form instead.
   mrgsolve's own parser rejects the mismatch on the very first line:
   ```
   Error: improper annotation format
    input: KA_GIV    = 0.50
    context: parse annotated parameter block (PARAM)
   Execution halted
   ```
   (`$CMT @annotated` in the same file is written correctly, in the
   colon-delimited form, and builds fine — only `$PARAM` has the mismatch.)
   Verification workaround: dropped `@annotated` from the `$PARAM` header
   (`$PARAM @annotated` → `$PARAM`) in scratch copies of both models; no
   parameter name or value changed.

2. **`SOLVERTIME` (mrgsolve's `_ODETIME_[0]` macro) is used inside
   `$TABLE`, where this mrgsolve build does not expose it.** The `$TABLE`
   block recomputes the hormonal-trigger capture from scratch using
   `SOLVERTIME` (duplicating the `$ODE`-local `HORMTRIG` calculation
   instead of reading it), which fails at the C++ stage once defect #1
   above is worked around:
   ```
   modelheader.h:112:20: error: '_ODETIME_' was not declared in this scope
     112 | #define SOLVERTIME _ODETIME_[0]
   ...  in expansion of macro 'SOLVERTIME' (the HORM_TRIGGER capture line)
   ```
   Verification workaround: replaced the `$TABLE` capture's from-scratch
   recomputation with a direct reference to the already-identical `$ODE`
   local, i.e. `capture HORM_TRIGGER = HORMTRIG;` instead of re-deriving it
   from `SOLVERTIME`. This is a pure aliasing change — `HORMTRIG` in `$ODE`
   is computed with the exact same formula the `$TABLE` line duplicated —
   so it introduces no numeric difference, applied identically to scratch
   copies of both the original and the refactored model.

Applied identically to scratch copies of both `aip_mrgsolve_model.R` and
`aip_mrgsolve_model_refactored.R`, purely to get a working verification
harness; neither checked-in file was changed by either workaround — both
still contain the mismatched `$PARAM @annotated` header and the
`SOLVERTIME`-in-`$TABLE` line exactly as originally written. See
`acute-intermittent-porphyria/aip_refactor_notes.md`.

**Fix upstream would be:** either drop `@annotated` from `$PARAM` (since
the body is plain-style) or rewrite the body into the colon-delimited
annotated form; and replace the `$TABLE` block's `SOLVERTIME`-based
recomputation of `HORM_TRIGGER` with a reference to the `$ODE`-local
`HORMTRIG` (or otherwise move the calculation somewhere `_ODETIME_` is in
scope).

---

## 40. `eosinophilic-esophagitis/eoe_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: three layered defects in `$INIT`/`$CMT`/`$CAPTURE`

Found while verifying a mepolizumab (MEPO) + cendakimab (CENDA) PK/effect
refactor (`eoe_mrgsolve_model_refactored.R`); all three reproduce
identically from the untouched original via the qspserver `mrgsolve_api`
container and are unrelated to either compound's own PK/effect math — they
block the model from building at all, refactored or not.

1. **`$INIT @annotated` uses `NAME = value : description` syntax
   throughout** (all 18 lines, e.g. `BUD_ESO  = 0     : Initial budesonide
   esoph (mg/L)`). mrgsolve 2.0.1's annotated-block parser requires
   colon-delimited `NAME : value : description` — the same form the same
   file's own `$PARAM @annotated` and `$CMT @annotated` blocks already use
   correctly:
   ```
   Error: improper annotation format
    input: BUD_ESO  = 0     : Initial budesonide esoph (mg/L)
    context: parse annotated init block (INIT)
   Execution halted
   ```
   Confirmed with a minimal 1-compartment reproduction outside this file:
   `A  = 0     : Initial A (mg)` in an `$INIT @annotated` block fails
   identically; switching only the separator to `A : 0 : Initial A (mg)`
   parses. Not specific to any one line or to spacing — every `$INIT` line
   in the file uses the `=`-separated form and all are rejected the same
   way.
2. **Once (1) is worked around, `$CMT @annotated` and `$INIT` are revealed
   to jointly redeclare the same 18 compartment names** — the same defect
   class as issue #34 (`chronic-hypothyroidism`): mrgsolve 2.0.1 treats
   `$INIT` as its own compartment-declaring block, not a companion to
   `$CMT`:
   ```
   invalid class "mrgmod" object: Duplicated model names: BUD_ESO BUD_SYS
   DUP_SC DUP_C DUP_P MEPO_SC MEPO_C CENDA_GUT CENDA_C IL13 IL5 EOTAX3
   EOS_BL EOS_ESO MAST_ESO FIBRO IGE_TOT EPBAR
   ```
3. **Once (2) is worked around, `$CAPTURE` is revealed to list six
   compartment names directly** (`IL13`, `IL5`, `EOTAX3`, `EOS_ESO`,
   `MAST_ESO`, `EPBAR`) — the same defect class as issue #34's second half
   and issue #30 (`sepsis`):
   ```
   invalid class "mrgmod" object: compartment should not be in $CAPTURE:
   IL13,IL5,EOTAX3,EOS_ESO,MAST_ESO,EPBAR
   ```

**Confirmed upstream:** all three reproduce from the untouched original
alone, via the qspserver `mrgsolve_api` container — a genuine, three-layer
version-level incompatibility, not a refactor-introduced or local-toolchain
issue.

**Verification workaround (in-memory only, not committed to either
file):** applied identically to scratch copies of `eoe_mrgsolve_model.R`
and `eoe_mrgsolve_model_refactored.R`: (a) the `$INIT @annotated` block
deleted and its 18 assignments moved into `$MAIN` using the modern
`<CMT>_0 = value;` idiom (e.g. `BUD_ESO_0 = 0;`), which resolves both (1)
and (2) in one step since it removes the malformed, redeclaring block
entirely; (b) the six compartment names removed from `$CAPTURE`
(compartment states are always present in mrgsolve's output regardless of
`$CAPTURE`, so this changes nothing about what is reported). No numeric
value or equation was touched by either change. Neither the tracked
original nor the delivered `_refactored.R` contains this workaround — both
still contain the `$INIT` annotation-format bug, the `$CMT`+`$INIT`
duplication, and the six-name `$CAPTURE` overlap exactly as written (under
the refactor's renamed `GUT_MEPO`/`CENT_MEPO`/`GUT_CENDA`/`CENT_CENDA` for
the four MEPO/CENDA compartments, since the refactor only renamed
identifiers, not the defect pattern itself — the renamed `$INIT` lines keep
the same `NAME = value : description` syntax, including the original's own
missing space in `CENDA_GUT= 0` reproduced as `GUT_CENDA= 0`). See
`eosinophilic-esophagitis/eoe_refactor_notes.md`.

**Fix upstream would be:** rewrite `$INIT`'s 18 lines from `NAME = value :
description` to `NAME : value : description`, **or** remove the `$INIT`
block entirely and set initial conditions via `<CMT>_0 = value;` in `$MAIN`
instead (this also resolves defect 2); and remove `IL13`/`IL5`/`EOTAX3`/
`EOS_ESO`/`MAST_ESO`/`EPBAR` from `$CAPTURE`.

---

## 41. `neurofibromatosis-type-1/nf1_mrgsolve_model.R` does not compile under mrgsolve 2.0.1, and its `$MAIN` initial-condition idiom is load-bearing on that very defect

Found while verifying a mirdametinib (MIR) + selumetinib (SEL) PK/effect
refactor (`nf1_mrgsolve_model_refactored.R`); reproduces identically from
the untouched original via the qspserver `mrgsolve_api` container and is
unrelated to either compound's own PK/effect math — it blocks the model
from building at all, refactored or not.

1. **`$CAPTURE` lists ten names that duplicate `$CMT` compartment names**
   (`RESIST`, `OPG_VOL`, `CNF_BURDEN`, `PAIN`, `QOL`, `VISION`, `LVEF`,
   `DERM_AE`, `CPK_AE`, `GROWTHZ`) — the same defect class as issues #30,
   #34, #35, #36, and #40: mrgsolve 2.0.1 refuses to build any model whose
   `$CAPTURE` repeats a compartment name:
   ```
   invalid class "mrgmod" object: compartment should not be in $CAPTURE:
   RESIST,OPG_VOL,CNF_BURDEN,PAIN,QOL,VISION,LVEF,DERM_AE,CPK_AE,GROWTHZ
   ```
2. **Once (1) is worked around by itself, a second and distinct defect
   appears**, not present in the original as checked in (the original never
   reaches this code path because defect 1 aborts the build first): `$MAIN`
   sets every disease-compartment initial condition with bare
   compartment-name assignment (e.g. `PERK = PERK_BASE;`, `RESIST = 0;`,
   `PAIN = PAIN0;` — 13 lines total under `if (NEWIND <= 1) { ... }`)
   instead of the standard mrgsolve `<cmt>_0 = value;` idiom. With the
   duplicate names removed from `$CAPTURE`, this bare form fails to
   compile:
   ```
   184:14: error: assignment of read-only reference 'PERK'
     PERK       = PERK_BASE;
   ...(identically for RESIST, PN_PROLIF, PN_T1/T2/T3, OPG_VOL, CNF_BURDEN,
        PAIN, QOL, VISION, LVEF, DERM_AE, CPK_AE, GROWTHZ)
   ```
   **Confirmed with a minimal 1-compartment isolated reproduction** (`$CMT
   X`, `$MAIN if (NEWIND<=1) { X = XBASE; }`, `$ODE dxdt_X = -0.1*X;`):
   `X = XBASE;` compiles fine when `X` is *also* listed in `$CAPTURE`
   (fails only at the later `validObject` step, same as defect 1 above),
   but fails to compile with exactly the "read-only reference" error once
   `X` is removed from `$CAPTURE` — i.e. the file's own `$MAIN`
   initial-condition idiom, as written, only compiles *because of* defect 1
   being present. The two defects are entangled, not independent: fixing
   (1) alone (the minimum needed to pass `validObject`) breaks compilation
   a different way, so a working build requires both to be addressed
   together. Switching to `X_0 = XBASE;` in the same minimal model compiles
   regardless of `$CAPTURE` membership and was confirmed (via
   `/run_simulation` on the minimal model) to actually initialize the state
   to the intended value (`X` starts at 5, not 0, then decays correctly) —
   i.e. `<cmt>_0` is the semantically-correct idiom, not merely a
   workaround that happens to compile.

**Confirmed upstream:** both reproduce from the untouched original alone,
via the qspserver `mrgsolve_api` container, and the minimal reproduction
isolates defect 2 from any of `nf1`'s own PK/PD content — a genuine,
version-level incompatibility (and a fragile inter-dependency between two
separate parts of the DSL spec), not a refactor-introduced or
local-toolchain issue. Because defect 1 always fires first, it is not
possible to observe whether the *original, unmodified* file's `$MAIN` block
would have correctly set initial conditions had defect 1 not existed — the
original, as checked in, cannot be built under mrgsolve 2.0.1 at all, by
either defect's route.

**Verification workaround (applied to the delivered `_refactored.R`, not
to `nf1_mrgsolve_model.R`):** the ten duplicate names were removed from
`$CAPTURE` (compartment states are always present in mrgsolve's output
regardless of `$CAPTURE`, so nothing is lost) and all 13 `$MAIN`
initial-condition lines were switched to the `<cmt>_0` idiom (e.g.
`PERK_0 = PERK_BASE;`). Both changes are needed together per defect 2
above; neither changes a parameter value or a model dynamic — confirmed by
the refactor's own verification, which reproduced the (identically
patched) original's output on two of its own dosing scenarios to a max
absolute difference of 0.0 across the full time grid. The identical pair of
changes was also applied to an in-memory-only scratch copy of the
*original* file's DSL, used solely to build a working comparison target;
`nf1_mrgsolve_model.R` itself was never edited. See
`neurofibromatosis-type-1/nf1_refactor_notes.md`.

**Fix upstream would be:** remove `RESIST`/`OPG_VOL`/`CNF_BURDEN`/`PAIN`/
`QOL`/`VISION`/`LVEF`/`DERM_AE`/`CPK_AE`/`GROWTHZ` from `$CAPTURE`, **and**
rewrite the `$MAIN` initial-condition block's 13 lines from bare
`NAME = value;` to `NAME_0 = value;`.

---

## 43. `disseminated-intravascular-coagulation/dic_mrgsolve_model.R` — `$ODE`-local `FREETPA` collides with a same-named `$TABLE` capture under mrgsolve 2.0.1

Found while verifying an ATRA (ATR) + tranexamic acid (TXA) PK/PD refactor
(`dic_mrgsolve_model_refactored.R`); unrelated to either compound's own math
and reproduces identically from the untouched original.

`$ODE` declares `double FREETPA = TPAnM * (FLOC + (1.0 - FLOC) / (1.0 +
PAInM / KPAI)) * ANX;` (the free-tPA term feeding `vPLN`, computed from the
live ODE state each step). `$TABLE` independently recomputes the same
quantity for reporting, from the captured/report-time state, under the
*same* name: `capture FREETPA = TPAnMT * (FLOC + (1.0 - FLOC) / (1.0 +
PAInMT / KPAI)) * (1.0 + KANX * BLAST / (KBL50 + BLAST));` (using its own
`_T`-suffixed locals `TPAnMT`/`PAInMT`, but reporting under the bare
`FREETPA` name). This is the same collision class as issue #35
(`dengue`) and the `sah`/`ted`/`ktx` refactor notes' `EFFECT_<STEM>`
capture-naming workaround: mrgsolve 2.0.1 auto-promotes every bare
`double NAME = expr;` in `$ODE` to a reportable anonymous-namespace member,
so `$TABLE`'s `capture FREETPA = ...;` then tries to redeclare a member
`$ODE`'s own `double FREETPA = ...;` already created:

```
145:11: error: redefinition of 'capture {anonymous}::FREETPA'
  145 |   capture FREETPA;
      |           ^~~~~~~
83:10: note: 'double {anonymous}::FREETPA' previously declared here
   83 |   double FREETPA;
      |          ^~~~~~~
```

**Confirmed upstream:** reproduces from the untouched original alone (no
refactor content involved), via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — the model does not build at all, for either the
original or the refactored file, until this is worked around.

**Verification workaround (in-memory only, not committed to either
file):** within `$TABLE` only (never touching `$ODE`, `$MAIN`, or any
numeric value), the colliding capture was renamed `FREETPA_CAP` (right-hand
side untouched). Applied identically to scratch copies of both
`dic_mrgsolve_model.R` and `dic_mrgsolve_model_refactored.R` purely so both
would build for the `/model_manifest` and `/run_simulation` comparison;
neither the tracked original nor the delivered `_refactored.R` contains
this workaround, and both still use the original's bare `FREETPA` capture
name exactly as written — so **neither currently builds against mrgsolve
2.0.1** without the same rename. See
`disseminated-intravascular-coagulation/dic_refactor_notes.md`.

**Fix upstream would be:** rename `$TABLE`'s `capture FREETPA = ...;` to a
non-colliding name (e.g. `FREETPA_CAP`, or the `_o`/`_T` prefix convention
this same file's author already used for `TPAnMT`/`PAInMT`), or rename the
`$ODE`-local `double FREETPA` instead, so mrgsolve's auto-promoted `$ODE`
member and `$TABLE`'s explicit `capture` declaration never share a name.

---

## 42. `type1-diabetes/t1dm_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT` and `$INIT` jointly redeclare all 20 compartments

Found while verifying an insulin (INS) + teplizumab (TEP) PK/effect
refactor (`t1dm_mrgsolve_model_refactored.R`); reproduces identically from
the untouched original via the qspserver `mrgsolve_api` container and is
unrelated to either compound's own math — it blocks the model from
building at all, refactored or not. Same defect class as issues #34
(`chronic-hypothyroidism`) and #41 defect 1's sibling pattern
(`neurofibromatosis-type-1`), independently present here.

`$CMT` names all 20 compartments with trailing comment annotations
(`Bm // 1. Beta-cell mass (normalised)`, ... `MG_avg // 20. Running mean
glucose`), then a separate `$INIT` block assigns each of the same 20 names
a starting value (`Bm = 1.0`, ... `MG_avg = 90.0`). mrgsolve 2.0.1 treats
`$INIT` as its own compartment-declaring block (an alternative to `$CMT`,
not a companion to it), so using both for the same names redeclares every
compartment twice:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: Bm CTL Treg Cpep
  Gp Gt SC1 SC2 Ic Ip X_ins Gcg HbA1c Ctep CGM_comp APC_integral Qsto1
  Qsto2 Qgut MG_avg
```

**Confirmed upstream:** reproduces from the untouched original alone, via
`POST /model_manifest` — the model does not build at all, for either the
original or the refactored file, until this is worked around.

**Verification workaround (in-memory only, not committed to either
file):** the `$INIT` block was deleted and its 20 assignments moved into
`$MAIN` using the standard `<CMT>_0 = value;` idiom (e.g. `Bm_0 = 1.0;`),
which declares no new compartment and changes no numeric value. A first
attempt instead moved these into a new `$PARAM` block using mrgsolve's
`<CMT>_0 = value` initial-parameter idiom (tried because it would also let
the API's `parameters` field reproduce the original R scenarios' own
`init(...)` overrides, which `/run_simulation` has no other mechanism for
— see "What could not be verified" in
`type1-diabetes/t1dm_refactor_notes.md`); that attempt does **not**
compile, because mrgsolve 2.0.1 already auto-generates a mutable
`<CMT>_0` reference for every compartment inside `$MAIN`, and declaring
the same name again via `$PARAM` collides with it:

```
491:15: error: conflicting declaration 'const double& Gt_0'
  491 | const double& Gt_0 = _A_0_[5];
469:9: note: previous declaration as 'double& Gt_0'
  469 | double& Gt_0 = _THETA_[53];
```

so the plain `$MAIN`-assignment form (matching issue #34's precedent
exactly) was used instead. Both patches were applied identically to
scratch copies of `t1dm_mrgsolve_model.R` and
`t1dm_mrgsolve_model_refactored.R` purely so both would build for the
`/model_manifest` and `/run_simulation` comparison; neither the tracked
original nor the delivered `_refactored.R` was changed — both still
contain the `$CMT`+`$INIT` duplication exactly as written (under the
refactor's renamed `GUT1_INS`/`GUT2_INS`/`CENT_INS`/`PERI_INS`/`EFF_INS`/
`CENT_TEP` for the six INS/TEP compartments, since the refactor only
renamed identifiers, not the defect pattern itself). See
`type1-diabetes/t1dm_refactor_notes.md`.

**Fix upstream would be:** remove the `$INIT` block and set every initial
condition via `<CMT>_0 = value;` in `$MAIN` instead (or drop `$CMT`'s
comment annotations and rely on `$INIT` alone, whichever the author
prefers).

## 44. `abdominal-aortic-aneurysm/aaa_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two independent defects

Found while verifying a Doxycycline/Propranolol/Statin PK/effect refactor
(`aaa_mrgsolve_model_refactored.R`); both reproduce identically from the
untouched original via the qspserver `mrgsolve_api` container and are
unrelated to the three compounds' own math — they block the model from
building at all, refactored or not.

**Defect 1: `$CAPTURE` repeats ten compartment names already in `$CMT`.**
Same defect class as several other files in this corpus (e.g. the
`rheumatoid-arthritis` entry above). The original's
`$CAPTURE Cp_doxy Ct_doxy Cp_stat Cp_prop MMP9 MMP2 ELAST COLLAG VSMC ILT
DIAM TNF ROSO MAC SBP_mmHg Wall_Stress_idx Rupture_P MMP9_reduction` lists
`MMP9, MMP2, ELAST, COLLAG, VSMC, ILT, DIAM, TNF, ROSO, MAC` — all ten of
which are also `$CMT` compartments. mrgsolve 2.0.1 refuses to build any
model whose `$CAPTURE` repeats a compartment name:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  MMP9,MMP2,ELAST,COLLAG,VSMC,ILT,DIAM,TNF,ROSO,MAC
```

**Defect 2 (a new pattern, not previously logged): `$PARAM` names
`MMP9_0`/`MMP2_0` collide with mrgsolve's own auto-generated `<CMT>_0`
initial-condition symbol** for the compartments `MMP9`/`MMP2`. Unlike issue
#42's `Gt_0` case (an author's *deliberate but broken* attempt to set an
initial value via `$PARAM`), this is an *incidental* collision: `MMP9_0`
and `MMP2_0` are ordinary baseline-value parameters (`MMP9_0 : 1.00 :
Baseline active MMP-9`), used throughout the disease equations for
normalization (`MMP9/MMP9_0`, `MMP2/MMP2_0`). The author separately tried
to set each compartment's own initial condition in `$MAIN` with
`MMP9_0 = MMP9_0;` / `MMP2_0 = MMP2_0;` (self-assignments, presumably
intended as "compartment `MMP9` starts at the value of parameter
`MMP9_0`") — but mrgsolve auto-generates an init-condition symbol also
named `MMP9_0`/`MMP2_0` for compartments `MMP9`/`MMP2`, so the two
same-named symbols (one a `const double&` parameter, one a mutable
`double&` compartment-init reference) collide at the C++ level:

```
183:9: error: conflicting declaration 'double& MMP9_0'
  183 | double& MMP9_0 = _A_0_[13];
      |         ^~~~~~
129:15: note: previous declaration as 'const double& MMP9_0'
  129 | const double& MMP9_0 = _THETA_[29];
```

(and symmetrically for `MMP2_0`, plus further redeclaration errors later in
the generated C++ from the same root cause).

**Confirmed upstream:** both reproduce from the untouched original alone,
via `POST /model_manifest` — the model does not build at all, for either
the original or the refactored file, until both are worked around.

**Verification workaround (in-memory only, not committed to either
file):** (a) the ten duplicated compartment names were stripped from
`$CAPTURE` (compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported);
(b) the `$PARAM` symbols `MMP9_0`/`MMP2_0` were renamed to
`MMP9_0_BASE`/`MMP2_0_BASE` everywhere **except** the one `$MAIN` line that
assigns each compartment's own auto-generated init symbol
(`MMP9_0 = MMP9_0_BASE;`) — a pure symbol rename, zero numeric or
equation change. Both patches were applied identically to scratch copies of
`aaa_mrgsolve_model.R` and `aaa_mrgsolve_model_refactored.R` purely so both
would build for the `/model_manifest` and `/run_simulation` comparison;
neither the tracked original nor the delivered `_refactored.R` was
changed — both still contain the `$CAPTURE` duplication and the
`MMP9_0`/`MMP2_0` collision exactly as written (this part of the file is
disease-PD, not one of the three refactored compounds, so it is untouched
by the refactor either way). See `abdominal-aortic-aneurysm/aaa_refactor_notes.md`.

**Also noted during the same verification (not a build defect, a modeling
one — not fixed, not worked around, just observed):** the disease-PD rate
constants are declared in day⁻¹ (`k_MAC_in`, `k_MMP9_syn`, `k_diam_grow`,
...) while the PK rate constants are declared in h⁻¹ (`ka_d`, `CL_d`,
...), yet both integrate on the same (hour-valued) time axis — confirmed by
the original R script's own `sim_time <- seq(0, 8760, by=24)` /
`Day = time/24` conventions. None of the day⁻¹ rates are divided by 24
before use in `$ODE`, so the disease system (`MAC`→`MMP9`→`ELAST`/`DIAM`→
`WALL_STR`→`MAC`, a positive feedback loop) runs roughly 24× faster than
its labeled units imply and diverges to floating-point overflow (`NaN`)
well within the model's own stated 1-year (`end=8760`) scenario horizon —
reproduced identically in both the original and the refactored model
(same time index, every scenario tried, treated or untreated). Not logged
as its own numbered defect since it does not block building or running;
recorded here and in the refactor notes because it explains why one of the
three compounds' own monotherapy scenario (Doxycycline alone) has to be
verified over a much shorter window than the original's own 8760h.

**Fix upstream would be:** for defect 1, remove the ten duplicated names
from `$CAPTURE`; for defect 2, rename the `$PARAM` baseline parameters
(e.g. `MMP9_BASE`/`MMP2_BASE`) so they no longer collide with the
compartments' auto-generated `<CMT>_0` symbols; for the day/hour mismatch,
divide every day⁻¹ disease-PD rate constant by 24 before use (or declare
and integrate the whole model on a day-valued time axis instead).

## 45. `alzheimers-disease/ad_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` repeats ten compartment names already in `$CMT`

Found while verifying a Donepezil/Lecanemab/Memantine PK/effect refactor
(`ad_mrgsolve_model_refactored.R`); confirmed to reproduce identically from
the untouched original via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — it blocks the model from building at all,
refactored or not, and is unrelated to any of the three compounds' own PK.

Same defect class as issue #44 (`abdominal-aortic-aneurysm`) and several
other files in this corpus. The original's
`$CAPTURE MMSE ADAS_Cog CDR_SB CSF_Ab42 CSF_pTau181 AmyloidPET_CL
TauPET_SUVr NfL_plasma Cp_Donepezil Ccns_Donepezil Cp_Memantine
Cp_Lecanemab AChE_inhibition NMDAR_occupancy Proto_neutralized AB_MONO
AB_OLIGO AB_PROTO AB_PLAQUE TAU_SOL TAU_PHOS TAU_AGG NEURO_INFLAM ACH SYN`
lists `AB_MONO, AB_OLIGO, AB_PROTO, AB_PLAQUE, TAU_SOL, TAU_PHOS, TAU_AGG,
NEURO_INFLAM, ACH, SYN` — all ten of which are also `$CMT` compartments.
mrgsolve 2.0.1 refuses to build any model whose `$CAPTURE` repeats a
compartment name:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  AB_MONO,AB_OLIGO,AB_PROTO,AB_PLAQUE,TAU_SOL,TAU_PHOS,TAU_AGG,
  NEURO_INFLAM,ACH,SYN
```

**Verification workaround, applied to the delivered
`ad_mrgsolve_model_refactored.R` per this file's own settled build-compat
policy (see `FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile
at all")** — since these ten names are ordinary `$CMT` compartments,
mrgsolve already reports their trajectories in every simulation output
regardless of `$CAPTURE`, so the ten names were simply dropped from the
`$CAPTURE` list; nothing about what is reported changes, and no equation
or parameter value is touched. The untouched original
(`ad_mrgsolve_model.R`) still carries the defect forward unfixed, per the
never-edit-upstream rule.

**Also noted while verifying (not a build defect, and not fixed or worked
around — recorded here because it explains a gap in the Donepezil
verification run): under the Donepezil-monotherapy scenario copied
verbatim from the original file's own `simulate_scenario()` call (10 mg
QD, `duration_years = 3`, `tgrid <- seq(0, 3*8760, by = 8)`), `ACH`
abruptly turns `NaN` at `time = 23528` h (~2.7 years into the 3-year
run) — with no obviously divergent value in the run-up (`ACH` sits flat at
1.2349 for the preceding ~500 output points, `AChE_inhibition` is
essentially 0, `MMSE` is smoothly declining) — and every downstream
compartment/output (`SYN`, `MMSE`, `COGNITION`, ...) is `NaN` from that
point to the end of the 3-year horizon. This looks like an `lsoda`
internal-step failure (the same class of solver breakdown the API surfaces
elsewhere as "excess work done... increase maxsteps") rather than a real
analytic blow-up in the equations, but this was not root-caused further.
Reproduces at the identical time index with byte-identical pre-blowup
values in both the original and the refactored model (confirmed via
`POST /run_simulation`), so it is a pre-existing long-horizon solver
artifact, not something the refactor introduced. The Memantine and
Lecanemab scenarios, run to the same 3-year horizon, showed no such
blowup.

**Fix upstream would be:** remove the ten duplicated names from
`$CAPTURE`; separately, investigate the `lsoda` failure at `t≈23528h`
under the Donepezil-monotherapy scenario (e.g. tightening solver
tolerances or capping `maxsteps` explicitly in `$SET`).

## 46. `bronchiectasis/bex_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$TABLE` ends in bare, unheadered `capture` lines, one of which also duplicates `$CMT` names

Confirmed upstream via `POST /model_manifest` on the untouched original
alone (`http://localhost:8007`, qspserver `mrgsolve_api` container). Two
layered defects, both inside the `$TABLE` block:

**Defect 1 (a new pattern, not previously logged): the block has no
`$CAPTURE` header at all.** The original's last three lines are

```
capture FEV1 FEV1_pct log10_BACT IL8 NE MCC AD BIOFILM EXAC
capture Exac_active AZM_lung_conc Tobra_lung_conc Cipro_ELF_conc
capture NEUT_abs Sputum_purulence
```

— lowercase `capture`, with no preceding `$CAPTURE` (or any other) block
marker. Every other file in this repo that captures output uses the real
`$CAPTURE` directive; here the word is just three bare identifiers dangling
at the end of `$TABLE`, which mrgsolve tries to compile as a C++ function
call to an undeclared `capture(...)`, and fails:

```
493:14: error: expected initializer before 'FEV1_pct'
  493 | capture FEV1 FEV1_pct log10_BACT IL8 NE MCC AD BIOFILM EXAC
      |              ^~~~~~~~
```

**Defect 2 (the same pattern already logged for other files, e.g. #30,
#34, #41): once given a real `$CAPTURE` header, the first line is revealed
to also list six names — `IL8, NE, MCC, AD, BIOFILM, EXAC` — that
duplicate existing `$CMT` compartment names.** mrgsolve 2.0.1 refuses to
build any model whose `$CAPTURE` repeats a compartment name (compartments
are always present in the output regardless of `$CAPTURE`):

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  IL8,NE,MCC,AD,BIOFILM,EXAC
```

**Verification-only workaround, later promoted into the delivered
`_refactored.R` per the guide's settled policy for a non-compiling
original (`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at
all"):** (a) each bare `capture ...` line was given a real `$CAPTURE`
header (`capture` → `$CAPTURE`, three separate block headers, one per
line — no symbol added, removed, or reordered); (b) the six compartment
names were dropped from the first `$CAPTURE` line (`$CAPTURE FEV1
FEV1_pct log10_BACT`), since compartments are already in the output
regardless. Both changes only touch which line/symbols are captured —
zero equations or parameter values changed anywhere. Confirmed via
`POST /model_manifest` that the patched DSL compiles, and via
`POST /run_simulation` that the three refactored compounds' outputs
(`AZM_lung_conc`, `Tobra_lung_conc`, `Cipro_ELF_conc`, `FEV1_pct`,
`log10_BACT`, and every PK compartment) are byte-identical between the
identically-patched original and the refactored model, so the compat fix
introduces no numeric change. See
`bronchiectasis/bex_refactor_notes.md` for the full account; the tracked
`bex_mrgsolve_model.R` itself is untouched and still carries both defects
exactly as written.

**Also noted during the same verification (not a defect, an
observation):** the original's own comments label the three refactored
compounds' absorption/clearance/volume parameters in per-hour units
(`Ka_AZM = 0.5, // Absorption rate (h-1)`, similarly for tobramycin and
ciprofloxacin) while every one of the R script's own dosing scenarios and
`tgrid` calls (`times_long <- seq(0, 365*2, by = 1)`, `ii = 7/freq_AZM`,
etc.) run on a day-valued time axis — the same day/hour mismatch class
already flagged in issue #44 for `abdominal-aortic-aneurysm`. Unlike #44,
this does not visibly blow up any scenario within the horizons actually
run (natural history and each of the three compounds' own dosing regimens
were verified stable to at least 60 days), so it is recorded here as an
observation rather than a root cause of any solver failure — but it means
the absorption/elimination half-lives implied by the code are ~24x longer
than the comments' stated h⁻¹ values would suggest on an hour-valued
clock. Not fixed, not worked around, present identically in both the
original and the refactored model.

**Fix upstream would be:** replace the three bare `capture` lines with a
single `$CAPTURE` block (or several), and drop the six compartment names
duplicated from `$CMT`; separately, either relabel the PK rate constants'
units to d⁻¹ or divide them by 24 before use, so the stated half-lives
match the day-valued simulation clock.

## 47. `chronic-lymphocytic-leukemia/cll_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: four compounding defects

Found while verifying an Ibrutinib/Venetoclax/Obinutuzumab PK/effect
refactor (`cll_mrgsolve_model_refactored.R`); all four reproduce
identically from the untouched original via the qspserver `mrgsolve_api`
container and are unrelated to any of the three compounds' own math — they
block the model from building at all, refactored or not. Each defect had
to be fixed before the next one's error became visible.

**Defect 1: `$INIT @annotated` is not valid syntax for a plain
`NAME = value` init block.** `$INIT @annotated` expects each line in the
`NAME : value : description` annotated form; the original's `$INIT` lines
(`DEPOT_IB = 0`, `CENT_IB = 0`, ...) have no description field:

```
Error: improper annotation format
 input: DEPOT_IB  = 0
 context: parse annotated init block (INIT)
```

**Defect 2: `$CMT` (annotated) + `$INIT` jointly redeclare all 18
compartments.** Same defect class as issues #34, #36, #42, #44. Once
defect 1 is worked around (dropping `@annotated`), mrgsolve still refuses
to build because both blocks declare the same 18 names:

```
invalid class "mrgmod" object: 1: Duplicated model names: DEPOT_IB
CENT_IB BTK_FREE BTK_OCC DEPOT_VEN CENT_VEN PERI_VEN BCL2_FREE BCL2_OCC
CENT_OBI PERI_OBI CD20_FREE CD20_OCC ALC BM_CLL LN_CLL MCL1_ADAPT NK_ACT
```

**Defect 3 (a new incidental-collision instance, same pattern as #44's
`MMP9_0`/`MMP2_0`): the pre-existing `$PARAM` symbol `ALC_0` collides with
mrgsolve's own auto-generated `<CMT>_0` initial-condition symbol** for the
compartment `ALC`, once defect 2 is worked around by moving `$INIT` into
`$MAIN` (mrgsolve's `<CMT>_0 = value;` idiom). `ALC_0` is an ordinary
baseline-value parameter (`ALC_0 : 50.0 : Baseline ALC (x1e9/L); typical
symptomatic patient`), used in `$TABLE`'s `ALC_pch`/`PD_flag`
calculations (`(ALC - ALC_0)/ALC_0*100`, `ALC > ALC_0*1.5`) — an
incidental clash (its value, 50.0, happens to equal the `ALC` compartment's
own initial value, likely why the original author never noticed two
different things sharing one name):

```
150:9: error: conflicting declaration 'double& ALC_0'
  150 | double& ALC_0 = _A_0_[13];
      |         ^~~~~
117:15: note: previous declaration as 'const double& ALC_0'
  117 | const double& ALC_0 = _THETA_[27];
```

**Defect 4: `$CAPTURE` repeats five compartment names already in
`$CMT`.** Same defect class as #34 and #44:

```
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE:
ALC,BM_CLL,LN_CLL,MCL1_ADAPT,NK_ACT
```

**Confirmed upstream:** all four reproduce from the untouched original
alone, via `POST /model_manifest` — the model does not build at all, for
either the original or the refactored file, until all four are worked
around, in the order above (each fix's absence hides the next defect's
error).

**Fix applied in the delivered `_refactored.R` (per this guide's settled
policy for "when the original doesn't compile at all"):** (1) dropped
`@annotated` from `$INIT`'s header, moot after (2); (2) the 18 `$INIT`
lines were moved into `$MAIN` as `<CMT>_0 = value;`; (3) the colliding
`$PARAM` symbol was renamed `ALC_0` → `ALC_BASELINE` everywhere it is used
as a parameter (the `$PARAM` line itself, and both `$TABLE` reads),
leaving `$MAIN`'s `ALC_0 = 50.0;` to refer only to the compartment's own
auto-generated init symbol; (4) the five duplicated names were dropped
from `$CAPTURE` (compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so nothing about what is reported changes). All
four are syntax-only, zero numeric or equation change, and are disclosed
in `chronic-lymphocytic-leukemia/cll_refactor_notes.md` alongside the
verification results proving they change nothing behaviorally. The
tracked original (`cll_mrgsolve_model.R`) still contains all four defects,
unfixed, exactly as written — only the fork-owned `_refactored.R` sibling
carries the fix.

**Also observed (not a build defect, an R-driver-script inconsistency —
not fixed, not worked around, just noted):** `dose_events()`'s scenario-6
branch explicitly sets `ev_list[["IB6"]] <- NULL # no ibrutinib`, yet
`scenarios[["6_triplet"]]`'s label reads "Triplet IB+VEN+OBI" — tracing
the `%in%` conditions, ibrutinib is dosed only for scenario 1 or scenario
5, so scenario 6 never doses it and is dosing-identical to scenario 4
(`VEN` ramp + `OBI` cycles, no `IB`) despite the different label. Recorded
here and in the refactor notes because verifying scenario 4 also verifies
scenario 6's actual (as opposed to labeled) behavior.

**Fix upstream would be:** for defect 1, drop `@annotated` from `$INIT`
(or give every line a third `: description` field); for defect 2, remove
the separate `$INIT` block once `$CMT` carries the same information (or
vice versa); for defect 3, rename the baseline parameter (e.g.
`ALC_BASELINE`) so it no longer collides with the compartment's
auto-generated `<CMT>_0` symbol; for defect 4, remove the five duplicated
names from `$CAPTURE`; for the scenario-6 mislabeling, either dose
ibrutinib in that branch or relabel it to match what it actually runs.

---

## 48. `membranous-nephropathy/mn_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` lists eighteen compartment names

`$CAPTURE` at the end of the model lists eighteen names that are already
`$CMT` compartments (`RTX_cent RTX_peri RTX_CD20_bound TAC_blood CPx_cent
CPx_metab CD20_B Plasma_cells Anti_PLA2R1 IgG_deposit Complement_MAC
Podocyte_inj GBM_thick Proteinuria Serum_alb eGFR AngII Aldosterone`).
mrgsolve 2.0.1 validates the compiled model object and rejects this outright
before any simulation can run:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: RTX_cent,RTX_peri,RTX_CD20_bound,TAC_blood,CPx_cent,CPx_metab,CD20_B,Plasma_cells,Anti_PLA2R1,IgG_deposit,Complement_MAC,Podocyte_inj,GBM_thick,Proteinuria,Serum_alb,eGFR,AngII,Aldosterone
```

**Confirmed upstream:** reproduced via `POST /model_manifest` on the
untouched original file alone, no changes involved.

**Why this matters here:** this file was the subject of a PK/PD refactor
(`mn_mrgsolve_model_refactored.R`, rituximab/cyclophosphamide/tacrolimus PK
blocks only, per `FORK_WORKFLOW_GUIDE.md` Part 2) whose mandatory
verification step requires actually building and running both the original
and the refactored model. Per the guide's settled policy for this situation
("When the original doesn't compile at all"), the fix applied here is
syntax-only and non-numeric: **the eighteen compartment names were removed
from `$CAPTURE`** (mrgsolve always includes every compartment's state in its
output regardless of whether it also appears in `$CAPTURE` — confirmed by
diffing the `/model_manifest` `outputPaths` before and after the change,
which still lists all eighteen compartments — so this changes nothing about
what is reported, only what compiles). This fix was applied **directly to
the delivered `mn_mrgsolve_model_refactored.R`**, not just to a scratch
copy, per the guide's newer settled answer for this class of defect — the
checked-in original (`mn_mrgsolve_model.R`) was left untouched, still
carrying the defect exactly as written. See
`membranous-nephropathy/mn_refactor_notes.md` for full disclosure and the
verification result this fix enabled (bit-identical match between original
and refactored, once both are built with this same fix applied to
in-memory-only copies of the original for comparison purposes).

**Fix upstream would be:** remove the eighteen duplicated names from
`$CAPTURE`; nothing else about the model needs to change.

## 49. `familial-mediterranean-fever/fmf_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: five layered defects, plus a comment/code TMDD mismatch for two compounds

Found while verifying a Colchicine/Anakinra/Canakinumab PK/effect-interface
refactor (`fmf_mrgsolve_model_refactored.R`); every defect below reproduces
identically from the untouched original, confirmed via the qspserver
`mrgsolve_api` container (`POST /model_manifest`), and none sits inside any
of the three refactored compounds' own PK blocks — all five are in
disease-PD scaffolding shared by the whole file.

**Defect 1: `$CMT` (bare, unannotated) and `$INIT` jointly redeclare all 22
compartments** — same defect class as issues #27/#34/#36/#42/#48's sibling
class:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: GUT_COL CENT_COL
  PERI_COL LEU_COL SC_ANA CENT_ANA SC_CANA CENT_CANA PERI_CANA RhoA
  Pyrin_p ASC Casp1 IL1b_pro IL1b_mat IL18 SAA CRP Neu_circ Neu_tis
  Att_sev AA_dep eGFR
```

**Defect 2 (a new pattern, only surfaces once defect 1 is fixed with the
modern idiom): six disease-baseline `$PARAM` values collide with their own
compartment's auto-generated `<CMT>_0` init symbol.** `IL18_0`, `SAA_0`,
`CRP_0`, `Neu_circ_0`, `Neu_tis_0`, and `eGFR_0` are ordinary
baseline-value parameters used only to seed the same-named compartment's
initial condition in the original `$INIT` (`IL18 = IL18_0`, `SAA = SAA_0`,
etc.). Once `$INIT` is merged into `$MAIN` via `<CMT>_0 = value;` (the fix
for defect 1), mrgsolve's own auto-generated init symbol for compartment
`IL18` is *also* named `IL18_0`, so the parameter and the auto-generated
symbol collide — the same incidental-collision pattern as `MMP9_0`/
`MMP2_0` in issue #44 (`abdominal-aortic-aneurysm`). (`IL1b_0`, which
seeds `IL1b_mat` and is also read elsewhere in the disease equations,
does *not* collide — the compartment is named `IL1b_mat`, not `IL1b`.)

**Defect 3: `$CAPTURE` repeats eleven compartment names already in `$CMT`**
— same defect class as issues #44/#45/#46/#48:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  IL1b_mat,IL18,SAA,CRP,Neu_circ,Neu_tis,Att_sev,AA_dep,eGFR,ASC,Casp1
```

**Defect 4: `$MAIN` and `$ODE` both declare an identically-named local
`double k_phos_eff = ...;`.** mrgsolve 2.0.1 hoists every block-local
`double` into the same anonymous-namespace scope, so the two declarations
collide:

```
37:10: error: redefinition of 'double {anonymous}::k_phos_eff'
27:10: note: 'double {anonymous}::k_phos_eff' previously declared here
```

The `$MAIN` copy is dead code — its value is never read again in `$MAIN`,
never captured, and has no effect on any equation; only the `$ODE` copy
(independently computed, identical formula) is actually used, by
`k_phos_drug`.

**Defect 5: `$MAIN`'s `if(NEWIND <= 1) { _nid++; }` references `_nid`, an
internal counter symbol not provided under this build:**

```
196:3: error: '_nid' was not declared in this scope
```

`_nid` appears nowhere else in the file; this entire block is dead code
(the increment is never read, captured, or used to affect any equation or
initial condition).

**Confirmed upstream:** all five defects reproduce from the untouched
original alone, via `POST /model_manifest` — the model does not build at
all, for either the original or the refactored file, until all five are
worked around, and none is related to Colchicine, Anakinra, or
Canakinumab's own PK.

**Verification workaround, applied to the delivered
`fmf_mrgsolve_model_refactored.R` per this fork's settled build-compat
policy (`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at
all")** — a delivered file nobody can actually run defeats the point of
the refactor: (1) the `$INIT` block was deleted, its 22 assignments moved
into `$MAIN` via the modern `<CMT>_0 = value;` idiom, values copied
verbatim; (2) the six colliding baseline parameters were renamed with a
`_BASE` suffix (`IL18_0_BASE`, `SAA_0_BASE`, `CRP_0_BASE`,
`Neu_circ_0_BASE`, `Neu_tis_0_BASE`, `eGFR_0_BASE`), with the new `$MAIN`
init line for each reading the renamed parameter — same numeric value,
same discoverability as a covariate, just a symbol rename plus a location
move; (3) the eleven duplicated names were dropped from `$CAPTURE` —
mrgsolve already reports every compartment's trajectory regardless of
`$CAPTURE`, so nothing about what is reported changes; (4) the dead,
unused `$MAIN` copy of `double k_phos_eff = ...;` was deleted, the `$ODE`
copy left untouched; (5) the dead `if(NEWIND <= 1) { _nid++; }` block was
deleted. None of these five fixes touches a single numeric value,
parameter, initial condition, or equation — confirmed by the exact-match
verification result in `fmf_refactor_notes.md` (max abs/rel deviation 0.0
across all 33 shared outputs, on all four of the file's own dosing
scenarios). The untouched original (`fmf_mrgsolve_model.R`) still carries
all five defects forward unfixed, per the never-edit-upstream rule.

**Also noted (not a build defect, a modeling-documentation one — not
fixed, not worked around, just observed and disclosed):** the file's own
calibration-notes header comment states "Anakinra: ... TMDD-based binding
to IL-1R1" and "Canakinumab: ... TMDD to free IL-1b", but neither
compound's actual code implements anything of the kind — both are plain
one- or two-compartment linear PK (`SC_ANA`/`CENT_ANA` and `SC_CANA`/
`CENT_CANA`/`PERI_CANA`, ordinary first-order absorption/elimination, no
`REC_FREE`/`COMPLEX` compartments, no `KON`/`KOFF`/`RTOT` anywhere), and
each compound's own effect term is a plain concentration-driven Emax/EC50
ratio, not an occupancy computed from a receptor-binding ODE system. Per
the refactor guide's instruction to reproduce a model's actual
mechanistic richness (never flatten a genuine TMDD system, but equally
never invent one a file's code does not contain), both compounds were
refactored as the linear PK they are actually coded as. See
`familial-mediterranean-fever/fmf_refactor_notes.md` for the full
disclosure.

**Fix upstream would be:** remove the `$INIT` block and set initial
conditions via `<CMT>_0 = value;` in `$MAIN` instead; rename the six
colliding baseline parameters; remove the eleven duplicated names from
`$CAPTURE`; delete the dead `$MAIN` copy of `k_phos_eff` and the dead
`_nid` block; separately, either implement genuine TMDD for Anakinra/
Canakinumab if that mechanistic richness is actually wanted, or correct
the header comment to describe the linear PK that is actually there.

---

## 50. `alcoholic-liver-disease/ald_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two independent, unrelated defects

Found while verifying a G-CSF/NAC/Prednisolone PK/effect-interface refactor
(`ald_mrgsolve_model_refactored.R`); both reproduce identically from the
untouched original via the qspserver `mrgsolve_api` container and are
unrelated to the three refactored compounds' own math — they block the
model from building at all, refactored or not.

**Defect 1: a `$TABLE`-local double collides with a same-named `capture`
line in the same block.** `$TABLE` computes `double prob_d90 = 1.0 / (1.0
+ exp(-logit_d90));` and, four lines later, `capture prob_d90 = prob_d90;`
— the capture declaration tries to redefine the same symbol mrgsolve
already auto-promoted from the `double` assignment:

```
59:11: error: redefinition of 'capture {anonymous}::prob_d90'
   59 |   capture prob_d90;
      |           ^~~~~~~~
56:10: note: 'double {anonymous}::prob_d90' previously declared here
   56 |   double prob_d90;
      |          ^~~~~~~~
```

Same defect class as issue #35 (`dengue`) and #43 (`disseminated-
intravascular-coagulation`), here inside a single `$TABLE` block rather
than across `$ODE`/`$TABLE`. Every other line in this same `$TABLE`
already reports under a distinct `<name>_out` capture (`MELD_out`,
`DF_out`, `ALT_out`, ...) — `prob_d90` is the one output that was captured
under its own bare computed name instead.

**Defect 2 (a new pattern, not previously logged): fourteen `$ODE` lines
across twelve compartments assign directly to their own compartment state
to clamp it**, e.g. `if(AA < 0) AA = 0;`, `if(H < 0.01) H = 0.01;`,
`if(PRED_C < 0) PRED_C = 0;`. mrgsolve 2.0.1 passes every `$CMT` state into
the generated `$ODE` function as a `const double&`, so writing to it is a
hard compile error for every one of the fourteen:

```
376:15: error: assignment of read-only reference 'AA'
  376 | if(AA < 0) AA = 0;
      |            ~~~^~~
...(and identically for GSH, KC, IL1B, NEUT, H (x2), F (x2), PRED_C,
NAC_C, GCSF_C, PTX_C, ANK_C)
```

Same underlying mechanism as issue #36 (`age-related-macular-
degeneration`), at larger scale (14 lines vs. 2) and confirmed to be the
same kind of no-op even where it might once have compiled: nothing in
`$ODE` reads any of these twelve compartment names again after its own
clamp line within the same evaluation, so — as issue #36 reasoned — only
`dxdt_*` feeds mrgsolve's integrator, meaning a direct reassignment of the
state variable itself could never have altered the integrated trajectory
either way, compilable or not.

**Confirmed upstream:** both reproduce from the untouched original alone,
via the qspserver `mrgsolve_api` container (`POST /model_manifest`) — the
model does not build at all, for either the original or the refactored
file, until both are worked around.

**Fix applied directly to the delivered `ald_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation: (a) the fourteen illegal clamp lines were deleted outright (see
the no-op reasoning above — removing a statement that could never have
affected the integrated trajectory changes nothing numeric); (b) the one
colliding capture was renamed `prob_d90` → `prob_d90_out`, matching the
`_out` convention every other line in the same block already uses (a pure
rename — the underlying `double prob_d90` computation is untouched, and
`prob_d90_out`'s value is identical to what `prob_d90`'s capture would
have reported had it compiled). Both patches are unrelated to G-CSF/NAC/
Prednisolone's own PK/PD, which are otherwise the only substantive changes
in the refactored file. The checked-in original (`ald_mrgsolve_model.R`)
was left untouched and still carries both defects exactly as written; the
same two patches were also applied to an in-memory-only scratch copy of
the original so it could build for the `/run_simulation` comparison.
Verified: with both patches applied to each side identically, every
shared `$CAPTURE`d output (`MELD_out`, `ALT_out`, `BILI_out`, `INR_out`,
`H_out`, `GSH_out`, `ROS_out`, `KC_out`, `NEUT_out`, `F_out`, `TNF_out`,
`IL1B_out`, `prob_d90_out`, plus each compound's own PK) matched the
refactored file exactly (0.0 deviation) across three dosing scenarios. See
`alcoholic-liver-disease/ald_refactor_notes.md`.

**Fix upstream would be:** rename the `$TABLE` local `prob_d90` (or its
capture) so the two symbols do not collide; delete the fourteen no-op
compartment-state clamp lines (or, if the author's intent — preventing a
visibly negative/out-of-range reported value — is worth keeping, clamp a
`$TABLE`-side recomputation instead, since clamping the `$ODE` state
itself was never effective).

---

## 51. `autoimmune-polyendocrinopathy/aps_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT` and `$INIT` jointly redeclare all 22 compartments

Found while verifying a Cyclosporine A (CSA) + Hydrocortisone (HC) + JAK
inhibitor (JAKI) PK/effect refactor (`aps_mrgsolve_model_refactored.R`);
reproduces identically from the untouched original via the qspserver
`mrgsolve_api` container and is unrelated to any of the three compounds'
own math — it blocks the model from building at all, refactored or not.
Same defect class already logged for other files (#38
`distal-renal-tubular-acidosis`, #42 `type1-diabetes`), independently
present here.

`$CMT` names all 22 compartments with trailing comment annotations
(`AIRE_func // AIRE function index (0-1 scale, dimensionless)`, ...
`Drug_HC // Hydrocortisone plasma (µg/dL equivalent)`), then a separate
`$INIT` block assigns each of the same 22 names a starting value
(`AIRE_func = 0.0`, ... `Drug_HC = 0.0`), in identical order. mrgsolve
2.0.1 treats `$INIT` as its own compartment-declaring block (an
alternative to `$CMT`, not a companion to it), so using both for the same
names redeclares every compartment twice:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: AIRE_func
  AutoT_pool Treg_pool AutoAb_adren AutoAb_PTG AutoAb_beta AutoAb_thy
  Adrenal_fn Cortisol_c PTG_fn PTH_plasma Ca_serum Beta_mass Insulin_p
  Glucose_p Thyroid_fn TSH_plasma FT4_plasma Drug_CsA Drug_Aba Drug_RTX
  Drug_JAKi Drug_HC
```

**Confirmed upstream:** reproduces from the untouched original alone, via
`POST /model_manifest` — the model does not build at all, for either the
original or the refactored file, until this is worked around.

**Fix applied directly to the delivered `_refactored.R`, per the guide's
settled policy for a non-compiling original (`FORK_WORKFLOW_GUIDE.md`,
"When the original doesn't compile at all"):** the `$CMT` block was
deleted outright and `$INIT` kept as-is. Since `$CMT`'s and `$INIT`'s
compartment orders and names were already identical (confirmed by direct
comparison — `$INIT` lists the same 22 names in the same order as `$CMT`,
just with `= value` appended), `$INIT` alone declares every compartment
with its original initial value and the same 1-based compartment index as
before; this is a pure syntax fix (removes a redundant declaration,
changes no name, order, or numeric value) confirmed via
`POST /model_manifest` (compiles; `outputPaths` lists the same 22
compartments in the same order as the original's own `$CMT` block) and
via `POST /run_simulation` (exact match, max abs diff 0.0, across three
scenarios — see `autoimmune-polyendocrinopathy/aps_refactor_notes.md` for
the full numeric comparison). This is not a fix to CSA/HC/JAKI's own PK —
it is identically necessary for every one of the file's other 19
compartments too, and would have blocked verification of any compound in
this file, not just the three in scope for this pass.

**Fix upstream would be:** remove the `$CMT` block and rely on `$INIT`
alone (or drop `$INIT`'s initial values and set them via `<CMT>_0 =
value;` in `$MAIN` instead, keeping `$CMT`'s comment annotations —
whichever the author prefers).

## 52. `igg4-related-disease/igg4rd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `TGFB^2` uses R's `^` power operator, invalid in C++

Found while verifying a Rituximab/Prednisone/Dupilumab PK/effect refactor
(`igg4rd_mrgsolve_model_refactored.R`); confirmed to reproduce identically
from the untouched original via the qspserver `mrgsolve_api` container
(`POST /model_manifest`) — it blocks the model from building at all,
refactored or not, and is entirely unrelated to any of the three compounds'
own PK (it sits in the fibrosis/TGF-β disease block).

The `$ODE` line

```
double TGFB_eff = TGFB^2 / (TGFB_EC50_FIB^2 + TGFB^2);
```

uses `^`, R's exponentiation operator, three times. mrgsolve compiles the
DSL block directly to C++, where `^` is bitwise XOR and requires integral
operands, not `double`:

```
292:35: error: invalid operands of types 'const double' and 'double' to
binary 'operator^'
```

**Verification workaround, applied to the delivered
`igg4rd_mrgsolve_model_refactored.R` per this file's own settled
build-compat policy (see `FORK_WORKFLOW_GUIDE.md`, "When the original
doesn't compile at all")** — rewrote the line as
`pow(TGFB,2) / (pow(TGFB_EC50_FIB,2) + pow(TGFB,2))`, algebraically
identical for all real inputs (confirmed by the verification run's exact
0.0 max-abs-diff result on `Fibrosis_idx`/`ECM`/every other output across
all six of the original's own dosing scenarios). No parameter value or
disease equation's shape was touched. The untouched original
(`igg4rd_mrgsolve_model.R`) still carries the defect forward unfixed, per
the never-edit-upstream rule.

**Fix upstream would be:** replace all three `X^2` with `pow(X,2)` (or
`X*X`) on that one line.

## 53. `igg4-related-disease/igg4rd_mrgsolve_model.R` Scenario 2 (`GUT_PRED` taper) uses `rate = -2` without a modeled infusion-duration parameter, so the original's own R script would crash if actually run

Found while verifying the same Rituximab/Prednisone/Dupilumab refactor as
issue #52, while reproducing the original file's own six named dosing
scenarios (as required by `FORK_WORKFLOW_GUIDE.md`'s verification
protocol) via the qspserver `mrgsolve_api` container's `/run_simulation`.
Unrelated to PRED's own PK/effect equations (`GUT_PRED`/`CENT_PRED`,
`KA_PRED`/`CL_PRED`/`V_PRED`/`F_PRED`, `Eimmu_PRED`), which are untouched
by this defect and verify exactly (see the refactor notes).

Scenario 2 ("Prednisone 40mg/d × 4wk taper over 6mo") builds every one of
its 168 daily `ev()` records with `rate = -2`:

```r
lapply(0:27, function(d) ev(time = d, cmt = "GUT_PRED",
                             amt = 40 * 1e6 / 358.4 / 40.0,
                             rate = -2))
```

In mrgsolve/NM-TRAN convention, `RATE = -2` tells the solver "the infusion
duration for this dose is given by a modeled parameter `D_<compartment>`
(or `Dn` for numbered compartments)" — the model itself must define that
duration parameter. `igg4rd_mrgsolve_model.R` defines no such parameter for
`GUT_PRED` (or any compartment) anywhere in `$PARAM`/`$MAIN`. Reproducing
this exact dosing (same `time`/`amt`/`cmt`/`rate` values, via
`POST /run_simulation`'s `dosing` records, which mrgsolve builds into the
same `as.ev()` event set `ev()` would) against the **untouched, only
build-compat-patched** original throws, confirming the defect is intrinsic
to the scenario's own dosing spec, not to anything introduced by the
refactor or by the API translation:

```
Error in (function (cond)  :
  error in evaluating the argument 'x' in selecting a method for function
  'as.data.frame': [mrgsolve] modeled infusion duration D_CMT or Dn must be
  positive when dosing record RATE is set to -2.
```

This means Scenario 2, exactly as written in the original file, would
throw this same error if a user actually ran the original's own R script
top to bottom — it is not merely slow or step-limited (the guide's
"shorten the window" allowance does not apply; no window is short enough
to make an undefined modeled-duration parameter exist).

**Verification workaround (not a fix, disclosed, not applied to either
tracked file):** for the sole purpose of exercising `GUT_PRED`/`CENT_PRED`/
`EFFECT_PRED` under a plausible daily-taper dosing pattern, the same 168
`time`/`amt`/`cmt` triples were resubmitted with `rate = 0` (a plain bolus)
instead of `rate = -2`, identically against both the original and the
refactored model. Result: exact match, max abs diff 0.0, on every shared
output across the full 899-point time grid — see
`igg4-related-disease/igg4rd_refactor_notes.md`. This substitution changes
the *administration-route representation* of the verification run only
(bolus vs. an infusion the model cannot actually run); it does not touch,
and is not evidence about, PRED's own PK parameters or effect equation,
which are identical in both files by construction (pure rename, see notes).

**Fix upstream would be:** either drop `rate = -2` (use `rate = 0`,
ordinary bolus dosing, matching how every other compound in this file is
dosed) or add a modeled `D_GUT_PRED` (or `D2`, `GUT_PRED` being the second
declared compartment) parameter with a sensible infusion duration.

## 54. `primary-sclerosing-cholangitis/psc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two layered defects

Found while doing the PK/PD refactor of this file's three census-tracked
compounds -- Bezafibrate (BEZ), Obeticholic acid (OCA), Ursodeoxycholic
acid (UDCA) -- per `FORK_WORKFLOW_GUIDE.md` Part 2. Confirmed by POSTing
the untouched original's own embedded `psc_model_code` DSL block, verbatim,
to the qspserver `mrgsolve_api` container's `POST /model_manifest`
(`http://localhost:8007`) -- no local R/mrgsolve install used. Neither
defect is inside any of the three refactored compounds' own PK/PD logic;
both sit in shared scaffolding (the file's `$INIT` block and one disease
baseline parameter) that every compartment in the file, PK and PD alike,
goes through.

**Defect 1: `$INIT @annotated` is declared with plain, non-annotated
`NAME = value` entries.** mrgsolve's annotated-block parser requires
`NAME : value : description` triples wherever a block is tagged
`@annotated`; this file's `$INIT @annotated` block instead uses the same
bare `NAME = value // comment` syntax as a non-annotated block (e.g.
`UDCA_gut = 0`), giving:

```
Error: improper annotation format
 input: UDCA_gut     = 0
 context: parse annotated init block (INIT)
Execution halted
```

**Defect 2, once defect 1 is worked around by simply dropping the
`@annotated` tag: `$CMT @annotated` and `$INIT` jointly redeclare all 27
compartments.** Same defect class as issues #27/#34/#36/#42/#48/#49/#51 --
mrgsolve 2.0.1 treats an (annotated) `$CMT` block plus a same-named `$INIT`
list as two conflicting declarations of the same compartment:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: UDCA_gut
  UDCA_plasma UDCA_bile OCA_gut OCA_plasma OCA_bile NorUDCA_bile BEZ_plasma
  LPS GutBarrier FXR_act BilePool HydroIndex IL17A TNFa IL6 Treg_IL10
  Cholangio_health Senescence HSC_act Col1a1 LOXL2 ALP Bilirubin Fibroscan
  PortalPressure CCA_risk
Calls: simcore_load_model ... initialize -> callNextMethod -> .nextMethod -> validObject
Execution halted
```

Applying the modern `<CMT>_0 = value;` idiom in `$MAIN` (moving all 27
init assignments there and deleting `$INIT`) surfaces one further,
incidental collision: the file's own `$PARAM LOXL2_0 : 0.20 : Baseline
LOXL2 cross-linking activity` is a plain, otherwise-unused baseline
parameter, but its name collides with the init-condition symbol
`LOXL2_0` that this idiom auto-generates for compartment `LOXL2` (whose
own `$INIT` value is a different number, 0.35 -- the two were never the
same value, confirming `LOXL2_0` the parameter was never actually used to
seed `LOXL2` the compartment in the original, either). Same incidental-
collision class as `IL18_0`/`SAA_0`/etc. in issue #49.

**Verification workaround, applied to the delivered
`psc_mrgsolve_model_refactored.R` per this fork's settled build-compat
policy (`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at
all")**: the `$INIT` block was deleted and its 27 assignments moved into
`$MAIN` via `<CMT>_0 = value;` (values copied verbatim), and the colliding
`LOXL2_0` parameter was renamed `LOXL2_0_BASE` (same value, 0.20; still
unused). Neither fix touches a single number, equation, or PK/PD
behavior -- confirmed by the exact-match verification results in
`primary-sclerosing-cholangitis/psc_refactor_notes.md` (max abs diff 0.0
on one scenario, 1e-17 floating-point noise on the other, both against the
untouched original patched with only these same two syntax fixes
in-memory for the verification run). The untouched original
`psc_mrgsolve_model.R` still carries both defects forward unfixed, per the
never-edit-upstream rule.

## 55. `patent-ductus-arteriosus/pda_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two independent syntax defects, both inside the three refactored compounds' own PK blocks

Found while doing the PK/PD refactor of this file's three census-tracked
compounds -- Acetaminophen (APAP), Ibuprofen (IBU), Indomethacin (IND) --
per `FORK_WORKFLOW_GUIDE.md` Part 2. Confirmed by POSTing the untouched
original's own embedded `pda_code` DSL block, extracted verbatim (R
single-quoted string, backslash-escapes undone), to the qspserver
`mrgsolve_api` container's `POST /model_manifest` (`http://localhost:8007`)
-- no local R/mrgsolve install used. Unlike most prior entries in this
file, both defects sit *inside* the scope compounds' own declarations
(not shared scaffolding), so per `FORK_WORKFLOW_GUIDE.md`'s "when the
original doesn't compile" policy, item 4, both were fixed directly as part
of this refactor rather than treated as a separate generic workaround.

**Defect 1: two `$PARAM @annotated` entries wrap their description onto a
second line without a `//` comment marker.** mrgsolve's annotated-block
parser requires each parameter on its own `NAME : value : description`
line; an indented continuation line with no `//` is read as a malformed
entry, not as more description text of the previous line:

```
KE0_IND  :   0.12    : effect-site equilibration, indomethacin (1/h) - slow,
                       representing slowly-reversible tight COX binding
```
and, further down:
```
IC50_APAP_K: 260.0   : unbound IC50, acetaminophen at renal COX (uM) - renal
                       peroxide tone is high, so the peroxidase mechanism is weak
```
gives, on the first offending line reached:
```
Error: improper annotation format
 input: representing slowly-reversible tight COX binding
 context: parse annotated parameter block (PARAM)
Execution halted
```
Both are indomethacin's/acetaminophen's own parameters (`KE0_IND` is the
scope-compound IND's ductal effect-site rate constant; `IC50_APAP_K` is
scope-compound APAP's renal-COX potency), so both are inside the refactor's
own scope, not incidental. A whole-file scan for the same pattern (an
indented, non-`//`, non-blank continuation line following a `$PARAM`
entry) found exactly these two and no others.

**Defect 2, once defect 1 is worked around: three multi-declarator
`double` lines drop the type from every name after the first.** Same
defect class as issue #37 (diabetic-ketoacidosis) and the dka refactor's
own finding -- mrgsolve 2.0.1's preprocessing of `$ODE` strips `double`
from every comma-separated declarator but the first:

```
double c1i = IBU1 / V1_IBU, c2i = IBU2 / V2_IBU;
double c1n = IND1 / V1_IND, c2n = IND2 / V2_IND;
double c1a = APAP1 / V1_APAP, c2a = APAP2 / V2_APAP;
```
gives:
```
846:22: error: 'c2i' was not declared in this scope
850:22: error: 'c2n' was not declared in this scope
854:24: error: 'c2a' was not declared in this scope
```
All three lines are the central/peripheral concentration setup for the
three scope compounds' own two-compartment PK -- again inside scope, not
incidental. A whole-file regex scan for the same multi-declarator pattern
found exactly these three and no others.

**Fix applied, in both cases directly in the delivered
`pda_mrgsolve_model_refactored.R`** (per the "inside the scope compound's
own block" branch of the guide's policy, not just a generic build-compat
footnote): the two continuation lines were folded onto one line each (no
text or value changed), and each multi-declarator line was split into two
separate `double` statements (same values, same order). Neither fix
touches a single number, equation, or PK/PD behavior -- confirmed by the
exact-match verification results in `pda_refactor_notes.md` (max abs diff
0.0 across seven scenarios spanning all three compounds individually, a
combination dosing, a sepsis/peroxide-tone case, and both a 10-day and a
full 90-day window), run against the untouched original patched with only
these same two syntax fixes in-memory for the verification call. The
untouched original `pda_mrgsolve_model.R` still carries both defects
forward unfixed, per the never-edit-upstream rule.

## 56. `pagets-disease/pbd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` lists fifteen compartment names

`$CAPTURE` at the end of the model lists fifteen names that are already
`$CMT` compartments (`ZA_cen ZA_bon CTN_cen DMB_cen RANKL_free OPG_free
OCpre OC OBpre OB BMD bsALP NTX CTX_s Pain`). mrgsolve 2.0.1 validates the
compiled model object and rejects this outright before any simulation can
run:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: ZA_cen,ZA_bon,CTN_cen,DMB_cen,RANKL_free,OPG_free,OCpre,OC,OBpre,OB,BMD,bsALP,NTX,CTX_s,Pain
```

**Confirmed upstream:** reproduced via `POST /model_manifest` on the
untouched original file alone, no changes involved.

**Why this matters here:** this file was the subject of a PK/PD refactor
(`pbd_mrgsolve_model_refactored.R`, Calcitonin/Denosumab/Zoledronic Acid PK
blocks only, per `FORK_WORKFLOW_GUIDE.md` Part 2) whose mandatory
verification step requires actually building and running both the original
and the refactored model. Per the guide's settled policy for this situation
("When the original doesn't compile at all"), the fix applied here is
syntax-only and non-numeric: **the fifteen compartment names were removed
from `$CAPTURE`** (mrgsolve always includes every compartment's state in
its output regardless of whether it also appears in `$CAPTURE` — confirmed
by diffing the `/model_manifest` `outputPaths` before and after the change,
which still lists all fifteen compartments — so this changes nothing about
what is reported, only what compiles). This fix was applied **directly to
the delivered `pbd_mrgsolve_model_refactored.R`**, not just to a scratch
copy, per the guide's settled answer for this class of defect — the
checked-in original (`pbd_mrgsolve_model.R`) was left untouched, still
carrying the defect exactly as written. See
`pagets-disease/pbd_refactor_notes.md` for full disclosure and the
verification result this fix enabled (exact match between original and
refactored across three dosing scenarios, once both are built with this
same `$CAPTURE` fix applied to in-memory-only copies of the original for
comparison purposes).

**Fix upstream would be:** remove the fifteen duplicated names from
`$CAPTURE`; nothing else about the model needs to change.
forward unfixed, per the never-edit-upstream rule.

## 57. `breast-cancer/bc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` lists twelve compartment names

`$CAPTURE` at the end of the model lists twelve names that are already
`$CMT` compartments (`TUMOR CSC ER_SIGNAL CDK46_ACT HER2_SIGNAL PD_L1
CD8_EFF TREG E2_PLASMA AROMATASE Ki67 CA153`). mrgsolve 2.0.1 validates the
compiled model object and rejects this outright before any simulation can
run:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: TUMOR,CSC,ER_SIGNAL,CDK46_ACT,HER2_SIGNAL,PD_L1,CD8_EFF,TREG,E2_PLASMA,AROMATASE,Ki67,CA153
```

**Confirmed upstream:** reproduced via `POST /model_manifest` on the
untouched original file's own DSL, no changes involved.

**Why this matters here:** this file was the subject of a PK/PD refactor
(`bc_mrgsolve_model_refactored.R`, Letrozole/Olaparib/Palbociclib/
Trastuzumab PK+PD blocks only, per `FORK_WORKFLOW_GUIDE.md` Part 2) whose
mandatory verification step requires actually building and running both the
original and the refactored model. Per the guide's settled policy for this
situation ("When the original doesn't compile at all"), the fix applied
here is syntax-only and non-numeric: **the twelve compartment names were
removed from `$CAPTURE`** (mrgsolve always includes every compartment's
state in its output regardless of whether it also appears in `$CAPTURE` —
confirmed by diffing `/model_manifest`'s `outputPaths` before and after,
which still lists all twelve compartments — so this changes nothing about
what is reported, only what compiles). This fix was applied **directly to
the delivered `bc_mrgsolve_model_refactored.R`**, not just to a scratch
copy, per the guide's settled answer for this class of defect — the
checked-in original (`bc_mrgsolve_model.R`) was left untouched, still
carrying the defect exactly as written. See
`breast-cancer/bc_refactor_notes.md` for full disclosure and the
verification result this fix enabled (exact match, max abs diff 0.0,
between original and refactored across all six of the original's own
dosing scenarios, once both are built with this same `$CAPTURE` fix applied
to in-memory-only copies of the original for comparison purposes).

**Fix upstream would be:** remove the twelve duplicated names from
`$CAPTURE`; nothing else about the model needs to change.

## 58. `x-linked-hypophosphatemia/xlh_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two layered defects

Two independent build defects, both surfacing only when `POST
/model_manifest` actually tries to build the untouched original file's own
DSL:

1. **`$CAPTURE` lists fourteen compartment names.** The final line,
   `$CAPTURE BURO_CP FGF23_NEUT NPT2 TMPGFR PHOS CALCITRIOL PTH BSAP RSS
   HEIGHTZ_XLH AGV_CALC_XLH SIXMWT WOMAC UCACR NEPHROCALC PHOSORAL_SIG
   CALC_CENT`, repeats fourteen names that are already `$CMT` compartments
   (`NPT2 TMPGFR PHOS CALCITRIOL PTH BSAP RSS HEIGHTZ_XLH SIXMWT WOMAC UCACR
   NEPHROCALC PHOSORAL_SIG CALC_CENT`). mrgsolve 2.0.1 rejects this at
   object-validation time, before any C++ compilation is attempted:
   ```
   Error in validObject(.Object) :
     invalid class "mrgmod" object: compartment should not be in $CAPTURE: NPT2,TMPGFR,PHOS,CALCITRIOL,PTH,BSAP,RSS,HEIGHTZ_XLH,SIXMWT,WOMAC,UCACR,NEPHROCALC,PHOSORAL_SIG,CALC_CENT
   ```
2. **`$MAIN`'s `if (NEWIND <= 1) { ... }` block assigns directly to
   compartment names** (`NPT2 = NPT2_BASE;`, `TMPGFR = TMPGFR0;`, etc., 12
   assignments total) to set per-individual initial conditions. This defect
   only surfaces once defect 1 above is worked around — with `$PLUGIN
   autodec` in this mrgsolve build, a bare compartment name in `$MAIN` is a
   read-only reference, so the C++ compiler rejects every one of these
   twelve assignments:
   ```
   190:14: error: assignment of read-only reference 'NPT2'
     190 |   NPT2       = NPT2_BASE;
         |   ~~~~~~~~~~~^~~~~~~~~~~
   ... (identical error for TMPGFR, PHOS, CALCITRIOL, PTH, BSAP, RSS,
        HEIGHTZ_XLH, SIXMWT, WOMAC, UCACR, NEPHROCALC)
   make: *** [/usr/lib/R/etc/Makeconf:211: inline-mread-source.o] Error 1
   ```

**Confirmed upstream:** both reproduced via `POST /model_manifest` on the
untouched original file's own DSL, no changes involved — defect 1 first
(on the as-written file), then defect 2 (once defect 1 alone was patched in
an in-memory-only copy used purely to reach the next build stage).

**Why this matters here:** this file was the subject of a PK/PD refactor
(`xlh_mrgsolve_model_refactored.R`, Burosumab/Oral calcitriol/Oral phosphate
PK+PD blocks — the model's only three compounds — per
`FORK_WORKFLOW_GUIDE.md` Part 2) whose mandatory verification step requires
actually building and running both the original and the refactored model.
Per the guide's settled policy for this situation ("When the original
doesn't compile at all"), the fixes applied are syntax-only and non-numeric:
**the fourteen compartment names were removed from `$CAPTURE`** (mrgsolve
always includes every compartment's state in its output regardless of
`$CAPTURE` membership — confirmed by diffing `/model_manifest`'s
`outputPaths` before and after, which still lists all of them — so this
changes nothing about what is reported, only what compiles), and **the
twelve direct compartment assignments in `$MAIN` were switched to the
`<cmt>_0` initial-value idiom** (`NPT2_0 = NPT2_BASE;`, etc.), mrgsolve's
standard mechanism for setting a compartment's initial value from `$MAIN`,
which is numerically identical to a direct assignment at `NEWIND <= 1`.
Both fixes were applied **directly to the delivered
`xlh_mrgsolve_model_refactored.R`**, not just to a scratch copy, per the
guide's settled answer for this class of defect — the checked-in original
(`xlh_mrgsolve_model.R`) was left untouched, still carrying both defects
exactly as written. See `x-linked-hypophosphatemia/xlh_refactor_notes.md`
for full disclosure and the verification result these fixes enabled (exact
match, max abs diff 0.0, between original and refactored across all three
of the original's own dosing scenarios exercising all three compounds,
once both are built with these same two fixes applied to in-memory-only
copies of the original for comparison purposes).

**Fix upstream would be:** remove the fourteen duplicated names from
`$CAPTURE`, and rewrite the twelve `if (NEWIND <= 1) { CMT = value; }`
assignments in `$MAIN` using the `<cmt>_0` idiom; nothing else about the
model needs to change.

## 59. `benign-prostatic-hyperplasia/bph_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two independent, unrelated defects

Found while verifying a Dutasteride/Finasteride/Tadalafil/Tamsulosin PK/
effect-interface refactor (`bph_mrgsolve_model_refactored.R`); both
reproduce identically from the untouched original via the qspserver
`mrgsolve_api` container and are unrelated to any of the four refactored
compounds' own math — they block the model from building at all,
refactored or not.

**Defect 1: the deprecated `_init_<CMT>` idiom.** `$MAIN`'s
`if(NEWIND <= 1) { ... }` block sets all 24 compartments' initial
conditions with the old `_init_<CMT> = value;` form (`_init_TAMS_GUT = 0;`,
`_init_TEST_P = TEST0 * 50.0;`, etc.). mrgsolve 2.0.1 does not accept this
symbol at all:

```
251:3: error: '_init_TAMS_GUT' was not declared in this scope
  251 |   _init_TAMS_GUT = 0;
      |   ^~~~~~~~~~~~~~
... (identically for every one of the 24 compartments)
```

Same defect class already logged as issue #29 (`polymyalgia-rheumatica`),
here in a different file.

**Defect 2: `$ODE`/`$TABLE`-local `double`s collide with same-named
`capture` statements.** Two separate instances of this:

- Four `$ODE`-local plasma-concentration variables (`CP_tams`, `CP_fina`,
  `CP_dut`, `CP_tad` — one per compound) are declared `double NAME = ...;`
  inside `$ODE` for use in that compound's own effect calculation, and the
  identical bare name is *also* the target of an old-style
  `capture NAME = NAME_out;` statement in `$TABLE` (`capture CP_tams =
  CP_tams_out;`, etc.). mrgsolve 2.0.1 auto-promotes any `$CAPTURE`d name
  to a class member, which collides with the `double` declared for the
  same name in `$ODE`:
  ```
  86:11: error: redefinition of 'capture {anonymous}::CP_tams'
     86 |   capture CP_tams;
        |           ^~~~~~~
  27:10: note: 'double {anonymous}::CP_tams' previously declared here
     27 |   double CP_tams;
        |          ^~~~~~~
  ```
  (identically for `CP_fina`, `CP_dut`, `CP_tad`)
- Eight `$TABLE`-local doubles collide with a same-named `capture` line in
  the same block (`DHT_inhibition_pct`, `DHT_PROST_inh_pct`,
  `PV_change_pct`, `IPSS_change`, `QMAX_change`, `PVR_change`,
  `PSA_change_pct`, `Alpha1_block_pct` — each computed as
  `double NAME = expr;` and then re-exposed as `capture NAME = NAME;`):
  ```
  97:11: error: redefinition of 'capture {anonymous}::QMAX_change'
     97 |   capture QMAX_change;
        |           ^~~~~~~~~~~
  80:10: note: 'double {anonymous}::QMAX_change' previously declared here
     80 |   double QMAX_change;
        |          ^~~~~~~~~~~~
  ```
  (identically for the other seven names)

Same defect class already logged as issues #35 (`dengue`), #43
(`disseminated-intravascular-coagulation`), and #50 (`alcoholic-liver-
disease`) — here across `$ODE`/`$TABLE` for the four compound
concentrations, and within a single `$TABLE` block for the other eight.

**Confirmed upstream:** both reproduce from the untouched original alone,
via the qspserver `mrgsolve_api` container (`POST /model_manifest`) — the
model does not build at all, for either the original or the refactored
file, until both are worked around.

**Fix applied directly to the delivered `bph_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation: (a) every `_init_<CMT> = value;` line was converted to the
modern `<CMT>_0 = value;` idiom (24 lines; no compartment, index, or
starting value changed); (b) the eight `$TABLE`-only collisions were fixed
by renaming just the intermediate `double` local (append `_calc`; the
exposed capture name and its arithmetic are unchanged); (c) the four
`$ODE`-vs-`$TABLE` collisions were fixed by dropping `double` from the
`$ODE` declaration of this refactor's own `C_TAMS`/`C_FINA`/`C_DUT`/
`C_TAD` (and `EFFECT_TAMS`/`EFFECT_FINA`/`EFFECT_DUT`/`EFFECT_TAD`, which
hit the identical collision against their own new `$TABLE` captures) —
a bare assignment to the member mrgsolve already promotes from the later
`capture NAME = ...;` line, confirmed empirically via `POST
/model_manifest` (mrgsolve scans the whole model source for capture
statements before compiling each block, so file order does not matter).
All three fixes are syntax-only and change no numeric value. The
checked-in original (`bph_mrgsolve_model.R`) was left untouched and still
carries both defects exactly as written; the same fixes (using the
original's own `CP_tams`/`CP_fina`/`CP_dut`/`CP_tad` names for fix (c))
were applied to an in-memory-only scratch copy of the original so it
could build for the `/run_simulation` comparison. Verified: with these
fixes applied to each side, every shared `$CAPTURE`d output matched the
refactored file exactly (max abs diff 0.0) across all six of the
original's own treatment scenarios (Watchful Waiting, Tamsulosin,
Finasteride, Dutasteride, Combination DUT+TAMS, Tadalafil), full 730-day/
732-point daily grid, no shortening needed. See
`benign-prostatic-hyperplasia/bph_refactor_notes.md`.

**Fix upstream would be:** replace each `_init_<CMT> = value;` line with
the modern `<CMT>_0 = value;` idiom; rename the four `$ODE`-local
`CP_tams`/`CP_fina`/`CP_dut`/`CP_tad` variables (or their `$TABLE`
captures) so they no longer collide; rename the eight `$TABLE`-local
doubles (or their captures) so they no longer collide with each other.

## 60. `sarcoidosis/sarc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: the deprecated `_init_<CMT>` idiom, plus two `$PARAM` names that collide with mrgsolve's own auto-reserved `<CMT>_0` symbols, plus a pre-existing `TREG0`/`TREG_0` typo

Found while verifying a prednisone (PRED) / prednisolone (PREDL) /
methotrexate (MTX) PK/effect-interface refactor
(`sarc_mrgsolve_model_refactored.R`); reproduces identically from the
untouched original via the qspserver `mrgsolve_api` container and is
entirely inside the disease-PD initial-condition block (`MAC_ACT`, `TH1`,
`TREG`, `TNF`, `IFNG`, `IL12`, `GRAN`, `FIBR`, `ACE_BM`, `CALIT`,
`SERUM_CA`, `SIL2R`, `FVC_P`) — unrelated to any of the three refactored
compounds' own PK/effect math.

**Defect 1: the deprecated `_init_<CMT>` idiom.** `$MAIN` sets all 13
disease-compartment initial conditions with the old
`_init_<CMT> = value;` form (`_init_MAC_ACT = MAC_ACT0;`,
`_init_TH1 = TH1_0;`, etc., 13 lines total). mrgsolve 2.0.1 does not
accept this symbol at all:

```
246:1: error: '_init_MAC_ACT' was not declared in this scope
  246 | _init_MAC_ACT  = MAC_ACT0;
      | ^~~~~~~~~~~~~
... (identically for the other 12 compartments)
```

Same defect class already logged as issues #29, #41, #42, #59 (and
others), here in a different file.

**Defect 2: two `$PARAM` names collide with mrgsolve's own auto-reserved
per-compartment `<CMT>_0` initial-value symbol.** The original declares
baseline params `TH1_0` and `IL12_0` (matching compartments `TH1`/`IL12`
plus a `_0` suffix) and `TREG_0` (matching compartment `TREG` plus `_0`).
mrgsolve 2.0.1 auto-reserves exactly this `<CMT>_0` name per compartment
as the initial-value override slot, so a user `$PARAM` of the identical
name produces a hard conflicting-declaration error, independent of and in
addition to Defect 1:

```
207:9: error: conflicting declaration 'double& TREG_0'
  207 | double& TREG_0 = _A_0_[9];
      |         ^~~~~~
147:15: note: previous declaration as 'const double& TREG_0'
  147 | const double& TREG_0 = _THETA_[17];
```

(identically for `TH1_0`/compartment `TH1` and `IL12_0`/compartment
`IL12`; `TREG_0`'s collision recurs a second time later in the generated
source, at the point the disease `$ODE` block re-reads the parameter).

**Defect 3: a pre-existing `TREG0`/`TREG_0` typo, independent of Defects 1
and 2.** The original declares the baseline param as `TREG_0` (with an
underscore, unlike every sibling baseline param in the same block —
`MAC_ACT0`, `TNF0`, `IFNG0`, `GRAN0`, `FIBR0`, `ACE0`, `CALIT0`, `CA0`,
`SIL2R0`, `FVC_P0` — all of which omit the underscore) but then reads it
back as `TREG0` (no underscore) in three separate places: the `$MAIN`
initial-condition line (`_init_TREG = TREG0;`) and twice in the
macrophage-suppression logic in `$ODE`
(`double treg_ratio = TREG0 / (treg + 1e-9);` and
`double treg_suppress = 1.0 / (1.0 + treg / (TREG0 * 2.0));`). `TREG0` is
never declared anywhere in the file, so as literally written the model
could never have compiled under any mrgsolve version that requires every
referenced symbol to exist — this is a genuine authoring typo, not a
version-compatibility issue:

```
248:18: error: 'TREG0' was not declared in this scope; did you mean 'TREG'?
  248 | _init_TREG     = TREG0;
      |                  ^~~~~
      |                  TREG
```

**Confirmed upstream:** all three reproduce from the untouched original
alone, via the qspserver `mrgsolve_api` container (`POST /model_manifest`)
— the model does not build at all, for either the original or the
refactored file, until all three are addressed.

**Fix applied directly to the delivered `sarc_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation, all three syntax-only and non-numeric: (a) `TH1_0`/`TREG_0`/
`IL12_0` renamed to `TH10`/`TREG0`/`IL120` (matching the no-underscore
convention every other baseline param in the same block already uses,
which also happens to fix Defect 3's typo in the same stroke, since
`TREG0` — the name the code actually reads — becomes the declared name);
(b) all 13 `_init_<CMT> = value;` lines converted to the modern
`<CMT>_0 = value;` idiom (no compartment, index, or starting value
changed). The identical pair of fixes was applied to an in-memory-only
scratch copy of the original (never to `sarc_mrgsolve_model.R` itself,
which still carries all three defects exactly as written) so it could
build for the `/run_simulation` comparison. Verified: with these fixes
applied to each side, every shared output matched the refactored file to
floating-point noise (max abs diff ~1e-10 to ~1e-12, consistent with a
pure structural rename) across both of the original's own treatment
scenarios that exercise the refactored compounds (Scenario 2, prednisone
40 mg taper; Scenario 3, prednisone 20 mg + MTX 10 mg/wk), full 104-week/
2913-point 6-hourly grid, no shortening needed. See
`sarcoidosis/sarc_refactor_notes.md`.

**Fix upstream would be:** rename the `$PARAM` baseline params `TH1_0`,
`TREG_0`, `IL12_0` to any name that does not collide with mrgsolve's
`<CMT>_0` reservation (e.g. the no-underscore form every sibling baseline
already uses) and update the `TREG0`/`TREG_0` mismatch consistently;
replace each `_init_<CMT> = value;` line with the modern
`<CMT>_0 = value;` idiom.
