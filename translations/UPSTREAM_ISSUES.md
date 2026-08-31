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

---

## 61. `diabetic-nephropathy/dn_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$INIT` packs two assignments per line with no comma, and (once that is worked around) `$CMT`+`$INIT` jointly redeclare all 19 compartments while `$CAPTURE` duplicates 8 of them

Found while verifying an ACE-inhibitor (ACEI) / ARB / SGLT2-inhibitor
(SGLT) / finerenone (FINE) PK/PD refactor
(`dn_mrgsolve_model_refactored.R`). All three defects reproduce
identically from the untouched original and are unrelated to any of the
four compounds' own math — the disease-PD block and every drug's PK are
affected identically, dosed or not.

**1. `$INIT` puts two `name=value` assignments on one line with no comma
separator**, e.g.:

```
$INIT
GI_acei=0 CENT_acei=0
GI_arb=0  CENT_arb=0
GI_sglt2=0 CENT_sglt2=0
GI_fine=0 CENT_fine=0
BG=8.5
...
```

mrgsolve parses a plain (non-`@annotated`) `$INIT`/`$PARAM` block as an R
`list(...)` built from its lines; a line holding two assignments with no
comma between them is not valid R argument syntax, so the whole file
fails before any C++ compilation is attempted:

```
Error in parse(text = paste0("list(", x, ")")) :
  <text>:1:16: unexpected symbol
1: list(GI_acei=0 CENT_acei
                   ^
```

**2. Once (1) is worked around (commas inserted), `$CMT` (bare,
unannotated) and `$INIT` are revealed to jointly redeclare the same 19
compartment names** — the same defect class as issues #34/#36
(`chronic-hypothyroidism`, `age-related-macular-degeneration`), independently
present here with a plain, non-`@annotated` `$CMT` this time. mrgsolve
2.0.1 treats `$INIT` as its own compartment-declaring block, not a
companion to `$CMT`:

```
invalid class "mrgmod" object: 1: Duplicated model names: GI_acei CENT_acei
GI_arb CENT_arb GI_sglt2 CENT_sglt2 GI_fine CENT_fine BG AGE_cmpt AngII_cmpt
TGF_cmpt ROS_cmpt ECM_cmpt Pod_cmpt UACR_cmpt Tub_cmpt Fib_cmpt GFR_cmpt
```

**3. `$CAPTURE` repeats 8 compartment names already declared in `$CMT`**
(`BG`, `AGE_cmpt`, `AngII_cmpt`, `TGF_cmpt`, `ROS_cmpt`, `ECM_cmpt`,
`Pod_cmpt`, `Fib_cmpt`) — the same defect class as issues #31/#34/#36,
independently present here:

```
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE:
TGF_cmpt,ROS_cmpt,ECM_cmpt,Pod_cmpt,Fib_cmpt,AngII_cmpt,BG,AGE_cmpt
```

**Confirmed upstream:** all three reproduce from the untouched original
alone via the qspserver `mrgsolve_api` container (`POST /model_manifest`)
— the model does not build at all, for either the original or the
refactored file, until all three are addressed.

**Fix applied directly to the delivered `dn_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation, both syntax-only and non-numeric: (a) the `$INIT` block was
deleted and its 19 assignments moved into a new `$MAIN` block using the
modern `<CMT>_0 = value;` idiom (no compartment, order, or starting value
changed) — this fixes both defects (1) and (2) in one stroke; (b) the 8
duplicated names were removed from `$CAPTURE` (compartment states are
always present in mrgsolve's own output regardless of `$CAPTURE`, so
nothing is lost — all 8 are still read straight off the simulated output
data frame by the R script). The identical pair of fixes was applied to
an in-memory scratch copy of the original (never to
`dn_mrgsolve_model.R` itself, which still carries all three defects
exactly as written) so it could build for the `/run_simulation`
comparison. Verified: with these fixes applied to each side, every
shared output matched the refactored file exactly (max abs diff = 0.0)
across both the original's own triple-therapy scenario (S6: ACEi +
SGLT2i + Finerenone, full 730-day/109-point weekly grid) and its ARB
monotherapy scenario (S2, shortened to 60 days — see the step-count note
in `dn_refactor_notes.md`). See `diabetic-nephropathy/dn_refactor_notes.md`.

**Fix upstream would be:** remove the `$INIT` block and set the 19
initial conditions via `<CMT>_0 = value;` in `$MAIN` instead (or insert
commas between the packed `$INIT` assignments and drop `$CMT`'s bare
declarations, relying on `$INIT` alone — whichever the author prefers);
remove `BG`/`AGE_cmpt`/`AngII_cmpt`/`TGF_cmpt`/`ROS_cmpt`/`ECM_cmpt`/
`Pod_cmpt`/`Fib_cmpt` from `$CAPTURE`.

## 63. `copd/copd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 19 compartments, `$CAPTURE` duplicates 9 of them, and 6 `$PARAM` baseline names collide with mrgsolve's own auto-reserved `<CMT>_0` symbols

Found while verifying a LAMA (Tiotropium) + LABA (Salmeterol) + ICS
(Budesonide) + PDE4I (Roflumilast) PK/effect refactor
(`copd_mrgsolve_model_refactored.R`); reproduces identically from the
untouched original via the qspserver `mrgsolve_api` container and is
unrelated to any of the four compounds' own math — it blocks the model
from building at all, refactored or not.

**Defect 1: `$CMT` and `$INIT` jointly redeclare all 19 compartments.**
`$CMT @annotated` names all 19 compartments with descriptions
(`LAMA_LUNG : LAMA lung depot (ng)`, ... `AE_rate_ann : Annualized
exacerbation rate`), then a separate `$INIT` block assigns each of the
same 19 names a starting value (`LAMA_LUNG = 0, LAMA_C = 0, ...`), in
identical order. mrgsolve 2.0.1 treats `$INIT` as its own
compartment-declaring block (an alternative to `$CMT`, not a companion to
it), so using both for the same names redeclares every compartment twice.
Same defect class already logged for other files (#38, #42, #51, #61),
independently present here.

**Defect 2: `$CAPTURE` lists 9 names that are already `$CMT` compartments**
(`IL8 NE_sput CRP Eos FEV1 Emph PVR AE_cum AE_rate_ann`), which mrgsolve
2.0.1 rejects outright. Same defect class as #56/#57/#61.

**Defect 3: six `$PARAM` baseline names collide with mrgsolve's own
auto-reserved per-compartment `<CMT>_0` initial-value symbol.** The
original declares `IL8_0`, `CRP_0`, `Eos_0`, `FEV1_0`, `Emph_0`, and
`PVR_0` as ordinary baseline-value parameters (documented as "Baseline
sputum IL-8", "Baseline serum CRP", etc.), three of which
(`IL8_0`/`FEV1_0`/`PVR_0`) are also genuinely read inside `$ODE`
(`IL8/IL8_0` ratio terms, `FEV1_target = FEV1_0*(1+E_BD)`,
`dxdt_PVR` referencing `PVR_0` twice) — the other three
(`CRP_0`/`Eos_0`/`Emph_0`) are declared but never read anywhere in the
DSL body (dead; `Eos_0` is referenced only from the original's own R-side
`dose_response` scenario as a parameter override, which itself has no
effect on the simulated trajectory — see below). Because `$CMT` also
declares compartments named `IL8`, `CRP`, `Eos`, `FEV1`, `Emph`, and
`PVR`, mrgsolve 2.0.1 auto-reserves the identical `<CMT>_0` name per
compartment as the initial-value override slot, so all six user `$PARAM`
declarations produce a conflicting-declaration error once Defect 1 is
worked around — same defect class as #60 (sarcoidosis), here with six
collisions instead of three.

**Also noted (not a build defect, not fixed, disclosed only): the
original's own R-side `dose_response` block passes `Eos_0=350` as a
parameter override intending to simulate an eosinophilic-phenotype
scenario, but `$INIT` hardcodes `Eos = 200` as a literal, never reading
the `Eos_0` param — so this override has zero effect on the actual
simulated Eos trajectory in the checked-in original. This is a
pre-existing authoring bug unrelated to any of the four refactored
compounds, reproduced faithfully (not fixed) in both the build-fix scratch
copy and the delivered `_refactored.R`, both of which likewise hardcode
the same literal initial value. See `copd/copd_refactor_notes.md`.**

**Confirmed upstream:** all three defects reproduce from the untouched
original alone, via `POST /model_manifest` — the model does not build at
all, for either the original or the refactored file, until all three are
worked around.

**Fix applied directly to the delivered `copd_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for a
non-compiling original (`FORK_WORKFLOW_GUIDE.md`, "When the original
doesn't compile at all"), all three syntax-only and non-numeric: (a) the
`$INIT` block was deleted and its 8 non-default initial values moved into
a new `$MAIN` block using the modern `<CMT>_0 = value;` idiom (all-zero PK
compartments need no explicit statement, matching their implicit-0
default); (b) the 9 duplicated names were removed from `$CAPTURE`
(compartment states are always present in mrgsolve's own output
regardless of `$CAPTURE`, so nothing is lost); (c) the six colliding
baseline params were renamed `IL80`/`CRP0`/`Eos0`/`FEV10`/`Emph0`/`PVR0`
(dropping the underscore, same convention as #60's fix) and every read
site updated to match, freeing the `<CMT>_0` symbols for their
mrgsolve-reserved use. The identical set of fixes was applied to an
in-memory-only scratch copy of the original (never to
`copd_mrgsolve_model.R` itself, which still carries all three defects
exactly as written) so it could build for the `/run_simulation`
comparison. Verified: with these fixes applied to each side, every shared
output matched the refactored file exactly (max abs diff = 0.0, max
relative diff = 0.0) across all six of the original's own named dosing
scenarios (Placebo; LAMA monotherapy; LABA+LAMA; ICS+LABA; Triple
LAMA/LABA/ICS; Triple+Roflumilast), full 8760h/365-day horizon, no
solver step-count issues, plus one extra ICS-monotherapy check. See
`copd/copd_refactor_notes.md`.

**Fix upstream would be:** remove the `$CMT` block'\''s duplication with
`$INIT` (or drop `$INIT` and set values via `<CMT>_0 = value;` in
`$MAIN`, keeping `$CMT`'\''s annotations); remove
`IL8`/`NE_sput`/`CRP`/`Eos`/`FEV1`/`Emph`/`PVR`/`AE_cum`/`AE_rate_ann`
from `$CAPTURE`; rename `IL8_0`/`CRP_0`/`Eos_0`/`FEV1_0`/`Emph_0`/`PVR_0`
to any non-colliding name and update their read sites; either wire
`Eos_0` into `$INIT`'\''s `Eos` starting value or remove the dead
`dose_response` override, whichever the author intends.

---

## 62. `chronic-myeloid-leukemia/cml_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: four layered defects, plus an unrelated dead density-dependence term that lets leukemic/normal cell counts blow up to Inf/NaN around year 4.7 in every single-TKI scenario

Found while verifying an Asciminib (ASC) / Dasatinib (DAS) / Imatinib (IMT)
/ Nilotinib (NIL) PK/PD refactor (`cml_mrgsolve_model_refactored.R`); all
five findings reproduce identically from the untouched original via the
qspserver `mrgsolve_api` container and are unrelated to any of the four
compounds' own PK/PD math.

**Build defects (the model does not compile at all until all four are
addressed):**

1. **`$INIT @annotated` uses `NAME = value : description` instead of the
   required colon-delimited `NAME : value : description` form** — the same
   mismatch class as issue #39, here in `$INIT` rather than `$PARAM`:
   ```
   Error: improper annotation format
    input: y0  = 1e4    : HSC initial count
    context: parse annotated init block (INIT)
   Execution halted
   ```
2. **Once (1) is worked around, `$CMT @annotated` (bare names + descriptions)
   and `$INIT @annotated` (values + descriptions) jointly redeclare all 22
   compartments** — the same defect class as issues #34/#36/#61:
   ```
   invalid class "mrgmod" object: 1: Duplicated model names: y0 y1 y2 y3
   x0q x0p x1 x2 x3 xs xsp xr xrp Gut_imt Cp_imt Cic_imt Cp_das Cic_das
   Gut_nil Cp_nil Cic_nil Gut_asc Cp_asc
   ```
3. **`$CAPTURE` repeats seven compartment names already declared in
   `$CMT`** (`Cp_imt`, `Cic_imt`, `Cp_das`, `Cic_das`, `Cp_nil`, `Cic_nil`,
   `Cp_asc`) — the same defect class as issues #31/#34/#36/#45/#48/#56/#57:
   ```
   invalid class "mrgmod" object: 2: compartment should not be in
   $CAPTURE: Cp_imt,Cic_imt,Cp_das,Cic_das,Cp_nil,Cic_nil,Cp_asc
   ```
4. **`$TABLE` declares `double BCRPCT`, `double LOG_IS`, and
   `double Resist_frac` as ordinary locals, and *also* has an explicit
   `capture NAME = NAME;` line for each, while the same three bare names
   are listed in `$CAPTURE`** — mrgsolve 2.0.1 auto-promotes every bare
   `$TABLE`-local `double` to a reportable class member (the same mechanism
   behind issue #35's `$ODE`-local collisions, here triggered in `$TABLE`
   instead), so the `capture`-driven redeclaration collides with it at the
   C++ stage:
   ```
   68:11: error: redefinition of 'capture {anonymous}::BCRPCT'
      68 |   capture BCRPCT;
   59:10: note: 'double {anonymous}::BCRPCT' previously declared here
      59 |   double BCRPCT;
   ```
   (identically for `LOG_IS` and `Resist_frac`). Every *other* `$TABLE`
   capture in this file (`CHR`, `MMR`, `MR4`, `MR45`, `WBC`, `LSC_total`,
   `LSC_quies`, `E_imt_out`, `E_das_out`, `E_nil_out`, `E_asc_out`) already
   uses a distinct source-variable name (e.g. `capture CHR = CHR_reached;`),
   so only these three collide — the author evidently intended the same
   distinct-name pattern everywhere but missed these three.

**Confirmed upstream:** all four reproduce from the untouched original
alone via `POST /model_manifest` — the model does not build at all, for
either the original or the refactored file, until all four are addressed.

**Fix applied directly to the delivered `cml_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation, all syntax-only and non-numeric: (1) `$INIT`'s `=` delimiters
were changed to `:`; (2) the standalone `$CMT @annotated` block was
removed, folding its descriptions into `$INIT @annotated` (which already
carries every compartment's value) — no compartment, order, or starting
value changed; (3) the seven duplicated names were removed from
`$CAPTURE` (compartment states are always present in mrgsolve's own
output regardless of `$CAPTURE` membership — confirmed by diffing
`/model_manifest`'s `outputPaths` before/after: all renamed compartments
are still listed); (4) the three colliding `$TABLE` locals were renamed
`BCRPCT_calc`/`LOG_IS_calc`/`Resist_frac_calc` (and every internal
reference to them updated), with the `capture BCRPCT = BCRPCT_calc;`-style
lines and the external `$CAPTURE` names (`BCRPCT`, `LOG_IS`, `Resist_frac`)
left exactly as before — purely a local-variable rename, the reported
values are numerically identical to the original. The identical four fixes
were applied to an in-memory scratch copy of the original (never to
`cml_mrgsolve_model.R` itself, which still carries all four defects exactly
as written) so it could build for the `/run_simulation` comparison.
Verified: with these fixes applied to each side, every shared output
matched the refactored file exactly (max abs diff = 0.0) across all four
of the original's own single-TKI dosing scenarios — see
`chronic-myeloid-leukemia/cml_refactor_notes.md`.

**Fifth finding, behavioral rather than a build blocker, found during the
same verification run and left untouched (not a refactor concern, logged
here since it was found while confirming the refactor's exact-match
result):** the `$ODE` block computes a density-dependent suppression term,
`double suppress = fmax(0.0, 1.0 - Kcomp * Ntotal);`, intended to cap total
marrow cellularity, but it is only ever used as `suppress * 0.0` inside
`dxdt_y0` — multiplied by a literal zero, making it a complete no-op.
Nothing in the model actually enforces the stated carrying capacity
(`K_total = 1e12`). Consequently, in every one of the four single-TKI
monotherapy scenarios the original file's own R script defines (imatinib
400mg/day, dasatinib 100mg/day, nilotinib 300mg BID — asciminib 40mg BID
did not, within the 10-year horizon tested), leukemic/normal cell
compartments grow without bound and the simulation reaches `Inf`/`NaN`
around t=41064h (~4.68 years) — identically, at the same output row, in
both the original and the refactored model (confirmed: this is a
pre-existing dynamical defect reproduced by the refactor, not introduced
by it). See the verification section of `cml_refactor_notes.md` for the
exact shared blow-up index.

**Fix upstream would be:** either drop `@annotated` in favor of a single
merged `$INIT` (or plain `$CMT`+`$INIT` without descriptions) to avoid the
duplicate-declaration class; drop the seven compartment names from
`$CAPTURE`; rename the three colliding `$TABLE` locals; and, for the
suppression no-op, either remove the dead `suppress` variable or wire it
into every hematopoietic/leukemic `dxdt_*` term as the author's own
comment ("density-dependent competition, Lotka-Volterra type") implies was
intended.

---

## 64. `essential-thrombocythemia/et_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 16 compartments, `$CAPTURE` duplicates 2 of them, and `$ODE` reads five `DOSE_*` parameters never declared in `$PARAM`; plus a shared, independently-reproducing `pow()`-on-noisy-near-zero-concentration `NaN` fragility in all three fractional-Hill compounds

Found while verifying an Anagrelide (ANA) / Hydroxyurea (HU) / Pegylated
interferon (PIFN) / Ruxolitinib (RUX) PK/effect refactor
(`et_mrgsolve_model_refactored.R`); all findings reproduce identically from
the untouched original via the qspserver `mrgsolve_api` container.

**Build defects (the model does not compile at all until all three are
addressed):**

1. **`$CMT` (bare, unannotated) and `$INIT` jointly redeclare all 16
   compartments** (`HSC MKP MK PLT TPO JAK2 SPL HU_C HU_P ANA_C ANA_P
   RUX_C RUX_P pIFN_C pIFN_P RISK_T RISK_MF`) — the same defect class as
   issues #34/#36/#42/#44/#48/#49/#51/#61/#62/#63:
   ```
   invalid class "mrgmod" object: Duplicated model names: HSC MKP MK PLT
   TPO JAK2 SPL HU_C HU_P ANA_C ANA_P RUX_C RUX_P pIFN_C pIFN_P RISK_T
   RISK_MF
   ```
2. **`$CAPTURE` repeats two compartment names already declared in `$CMT`**
   (`RISK_T`, `RISK_MF`) — the same defect class as issues #31/#34/#36/#45/
   #48/#56/#57/#62:
   ```
   invalid class "mrgmod" object: compartment should not be in $CAPTURE:
   RISK_T,RISK_MF
   ```
3. **`$ODE` reads `DOSE_HU`, `DOSE_ANA`, `DOSE_RUX`, `DOSE_pIFN`, and
   `DOSE_ASA` directly in five `dxdt_*`/effect lines, but none of the five
   is ever declared in `$PARAM`** — a new defect class for this repo (not
   a name collision, a missing declaration): the original's own R wrapper
   only ever sets these five names via `param()`/`mrgsim(param = ...)` at
   the R side, which silently no-ops on an unknown name for some mrgsolve
   versions but is fatal at `mcode()`-time here, since the generated C++
   `$ODE` function references symbols that were never declared anywhere in
   the model:
   ```
   199:22: error: 'DOSE_HU' was not declared in this scope
   202:22: error: 'DOSE_ANA' was not declared in this scope
   205:22: error: 'DOSE_RUX' was not declared in this scope
   208:24: error: 'DOSE_pIFN' was not declared in this scope
   215:10: error: 'DOSE_ASA' was not declared in this scope
   ```
   (`SOLVERTIME`, also used in these same lines, is *not* flagged — it is
   this mrgsolve version's valid in-`$ODE` accessor for the solver's own
   advancing time, per issue #27's finding, `#define SOLVERTIME
   _ODETIME_[0]`.)

**Confirmed upstream:** all three reproduce from the untouched original
alone via `POST /model_manifest` — the model does not build at all, for
either the original or the refactored file, until all three are addressed.

**Fix applied directly to the delivered `et_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for this
situation, all syntax-only and non-numeric: (1) the standalone `$CMT`
block was removed, `$INIT` alone declaring every compartment (same fix as
issue #51); (2) `RISK_T`/`RISK_MF` were dropped from `$CAPTURE`
(compartment states are always present in mrgsolve's own output
regardless of `$CAPTURE` membership — confirmed via `/model_manifest`'s
`outputPaths`); (3) `DOSE_HU = 0`, `DOSE_ANA = 0`, `DOSE_RUX = 0`,
`DOSE_PIFN = 0`, `DOSE_ASA = 0` were added to `$PARAM`, matching the
defaults the original's own R wrapper already used
(`mod2 <- mod %>% param(DOSE_HU = 0, ...)`) — no value invented beyond
what the original itself already treated as the no-dose default. The
identical three fixes were applied to an in-memory scratch copy of the
original (never to `et_mrgsolve_model.R` itself, which still carries all
three defects exactly as written) so it could build for the
`/run_simulation` comparison. Verified: with these fixes applied to each
side, every shared output matched the refactored file to floating-point/
JSON-rounding precision across all seven of the original's own dosing
scenarios, up to the shared `NaN` fragility described next — see
`essential-thrombocythemia/et_refactor_notes.md`.

**Fourth finding, a numerical fragility rather than a build blocker, found
during the same verification run and left untouched (not a refactor
concern, logged here since it was found while confirming the refactor's
exact-match result):** three of the four compounds' Hill effect terms use
a non-integer exponent (`gam_HU = 1.20`, `gam_ANA = 1.50`, `gam_RUX =
1.80`; only `gam_pIFN = 1.00` is exempt). Each compound's own PK
compartment decays, over the original's 730-day simulation horizon, to
concentrations many orders of magnitude below any value with remaining
physiological meaning (sub-`1e-8`, against `EC50` values of `0.008`–`3.5`)
— at that scale, ordinary floating-point/adaptive-solver roundoff can
momentarily push the compartment's value slightly *negative* even though
the true solution never goes below zero. `pow(negative, non-integer
exponent)` is `NaN` in C++, and once any single `dxdt_*` evaluates to
`NaN`, the solver's own state vector is `NaN` from that point forward,
cascading into every other output (`PLT_count`, `MKP`, etc.) even though
those compartments' own trajectories were never physically unstable.
Confirmed present in the untouched original alone (no refactor changes
involved): in the four single-compound dosing scenarios the original's
own R script defines (HU 500 mg/d, HU 1500 mg/d, ANA 2 mg/d, RUX 20 mg/d
BID; PIFN 90 µg/week never triggers it, consistent with its integer
exponent), the original blows up to `NaN` at t = 158, 153, 33, and 74 days
respectively. The refactored model — mathematically equivalent PK,
reorganized into a real depot/central/peripheral ODE system instead of
the original's closed-form `exp()` term — hits the *same* fragility
independently, but at different times (t = 119, 129, 29, and 47 days
respectively for the same four scenarios), because the two formulations
take different floating-point paths to the same near-zero values. Every
output the two models share matches to floating-point/JSON-rounding
precision (`≤ 1e-4` absolute, at values of order `1`–`900`) right up to
each model's own independent blow-up point — see the verification section
of `et_refactor_notes.md` for the exact per-scenario numbers. This is a
pre-existing fragility of applying a fractional Hill exponent to a
compartment whose true value can be numerically indistinguishable from
zero, not a discrepancy introduced by the refactor.

**Fix upstream would be:** remove the `$CMT` block (or drop `$INIT` and
set values via `<CMT>_0` in `$MAIN` instead); drop `RISK_T`/`RISK_MF` from
`$CAPTURE`; declare `DOSE_HU`/`DOSE_ANA`/`DOSE_RUX`/`DOSE_pIFN`/`DOSE_ASA`
in `$PARAM`; and, for the `NaN` fragility, clamp each drug compartment's
concentration to a non-negative floor before it is raised to a fractional
power, e.g. `pow(fmax(HU_C, 0.0), gam_HU)`.

---

## 65. `myotonic-dystrophy/dm1_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 19 compartments, and (once that is worked around) a `$PARAM` baseline name collides with mrgsolve's own auto-reserved `<CMT>_0` symbol

Found while verifying a Mexiletine PK/effect refactor
(`dm1_mrgsolve_model_refactored.R`); reproduces identically from the
untouched original via the qspserver `mrgsolve_api` container and is
unrelated to mexiletine's own PK/PD math — it blocks the model from
building at all, refactored or not.

**Defect 1: `$CMT` and `$INIT` jointly redeclare all 19 compartments.**
`$CMT` names all 19 compartments with inline comments (`MEX_GUT // mg`,
`MEX_CENT // mg`, ... `FVC_PCT // % predicted FVC`), then a separate
`$INIT` block assigns each of the same 19 names a starting value
(`MEX_GUT = 0`, `MEX_CENT = 0`, ...), in identical order. mrgsolve 2.0.1
treats `$INIT` as its own compartment-declaring block (an alternative to
`$CMT`, not a companion to it), so using both for the same names
redeclares every compartment twice:
```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: MEX_GUT
  MEX_CENT MEX_PERI ASO_PLASMA ASO_MUSCLE ASO_NUCL CUG_FOCI MBNL1_FREE
  CUGBP1_ACT CLCN1_FETAL SERCA_FETAL INSR_FETAL MYOTONIA GRIP_STR
  MUSCLE_MASS PR_INT QTc_INT HOMA_IR FVC_PCT
```
Same defect class already logged for other files (#38, #42, #51, #61,
#63).

**Defect 2: a `$PARAM` baseline name collides with mrgsolve's own
auto-reserved per-compartment `<CMT>_0` initial-value symbol.** The
original declares `HOMA_IR_0 = 3.5` as an ordinary baseline-value
parameter ("DM1 baseline HOMA-IR"), genuinely read once in `$MAIN`
(`HOMA_IR_target = HOMA_IR_0 * (1.0 + 1.5*(INSR_FETAL - INSR_fetal_norm))`).
Because `$CMT`/`$INIT` also declare a compartment named `HOMA_IR`,
mrgsolve 2.0.1 auto-reserves the identical name `HOMA_IR_0` as that
compartment's initial-value override slot, so the user `$PARAM`
declaration produces a conflicting-declaration error once Defect 1 is
worked around:
```
155:9: error: conflicting declaration 'double& HOMA_IR_0'
  155 | double& HOMA_IR_0 = _A_0_[17];
131:15: note: previous declaration as 'const double& HOMA_IR_0'
  131 | const double& HOMA_IR_0 = _THETA_[52];
```
Same defect class as #59/#60/#63, here with a single collision (checked:
no other `$PARAM` name in this file matches any of the other 18
compartments' `<CMT>_0` form — `Muscle_mass_0` and `FVC_0` are
mixed/short-case names that do not match `MUSCLE_MASS_0`/`FVC_PCT_0`
exactly, and are therefore not affected).

**Confirmed upstream:** both defects reproduce from the untouched
original alone, via `POST /model_manifest` — the model does not build at
all, for either the original or the refactored file, until both are
worked around.

**Fix applied directly to the delivered `dm1_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for a
non-compiling original (`FORK_WORKFLOW_GUIDE.md`, "When the original
doesn't compile at all"), both syntax-only and non-numeric: (a) the
`$CMT` block was removed (its per-compartment comments folded into a
single documentation comment immediately above `$MAIN`), leaving `$INIT`
as the sole compartment declaration, same order, same 1-based compartment
numbers; (b) `HOMA_IR_0` was renamed `HOMA_IR_BASELINE` and its one read
site in `$MAIN` updated to match, freeing the `<CMT>_0` symbol for its
mrgsolve-reserved use. The identical two fixes were applied to an
in-memory-only scratch copy of the original (never to
`dm1_mrgsolve_model.R` itself, which still carries both defects exactly
as written) so it could build for the `/run_simulation` comparison.
Verified: with these fixes applied to each side, every shared output
matched the refactored file exactly (max abs diff = 0.0) across both a
no-dosing natural-history check and the original's own mexiletine 300 mg
TID (q8h) dosing regimen with `MEX_ON=1`, 30-day window (shortened from
the original's own 365-day scenario only for verification-request
convenience, not because of any solver step-count issue — no step-count
problem was encountered). See `myotonic-dystrophy/dm1_refactor_notes.md`.

**Fix upstream would be:** remove the `$CMT` block (or drop `$INIT` and
set values via `<CMT>_0 = value;` in `$MAIN` instead, keeping `$CMT`'s
comments); rename `HOMA_IR_0` to any non-colliding name (e.g.
`HOMA_IR_BASELINE`) and update its one read site.

## 66. `spinal-muscular-atrophy/sma_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$TABLE` re-declares `CMAP`/`HFMSE`/`RULM` a second time via self-referential `capture NAME = NAME;`; plus two unrelated, non-blocking reporting/dead-code quirks

Found while refactoring nusinersen/risdiplam/onasemnogene-abeparvovec
(Zolgensma) PK/PD per this file's `driver-patches/data/
compound_perturbation_census.md` row (see also the corrected-compound-
name note below). Reproduces identically from the untouched original via
the qspserver `mrgsolve_api` container's `/model_manifest`, unrelated to
any compound's own PK.

**Defect (build-blocking): `$TABLE` declares `double CMAP = ...;` and
then immediately writes `capture CMAP = CMAP;`.** mrgsolve's `capture
NAME = expr;` shorthand both declares `NAME` as a `double` *and*
registers it for output; when `NAME` is already a declared `double` in
the same block (here, because the report variable happens to share its
own capture name), the shorthand's own re-declaration collides with the
one already in scope:

```
65:11: error: redefinition of 'capture {anonymous}::CMAP'
   65 |   capture CMAP;
      |           ^~~~
56:10: note: 'double {anonymous}::CMAP' previously declared here
   56 |   double CMAP;
      |          ^~~~
```

and identically for `HFMSE` (declared as `double HFMSE = ...;`, captured
as `capture HFMSE = HFMSE;`) and `RULM` (`double RULM = ...;` /
`capture RULM = RULM;`). Three other derived variables in the same
`$TABLE` block — `FVC_pct`/`capture FVC`, `CHOP`/`capture CHOP_INTEND`,
and the file's `E7I_out`/`capture E7_inclusion` — do **not** hit this,
because their own local-variable name already differs from the capture
name they feed; only the three whose local name is spelled identically
to its own capture name collide. Same underlying mrgsolve mechanism as
the block-local-`double`-hoisting defects already logged for other files
(#48, #51 and others in the #45-#64 run), here triggered by a
self-referential `capture X = X;` rather than a cross-block name clash.

**Confirmed upstream:** reproduces from the untouched original alone via
`POST /model_manifest`, no refactor changes involved.

**Fix applied directly to the delivered
`sma_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"),
syntax-only and non-numeric: the three colliding locals were renamed
`CMAP_val`/`HFMSE_val`/`RULM_val` (the same pattern the original already
used successfully for `FVC_pct`/`CHOP`), and their `capture` lines
updated to read from the renamed local instead of re-declaring the
capture name itself (`capture CMAP = CMAP_val;`, etc.). The identical
rename was applied to an in-memory-only scratch copy of the original
(never to `sma_mrgsolve_model.R` itself, which still carries the defect
exactly as written) so it could build for the `/run_simulation`
comparison. Verified: with this fix applied to each side, all 13 of the
original's own `$CAPTURE` outputs matched the refactored file exactly
(max abs diff = 0.0) across all six of the original's own scenarios
(untreated natural history; nusinersen ENDEAR-style loading+maintenance;
risdiplam SUNFISH-style daily dosing, adult and pediatric weight-based;
Zolgensma single IV dose; nusinersen late-start). See
`spinal-muscular-atrophy/sma_refactor_notes.md`.

**Two additional findings, non-blocking (disclosed, not fixed):**

1. **The `E7_inclusion` capture output is not the same quantity that
   actually drives disease dynamics.** `$ODE`'s `E7I_current` (which
   integrates into `FL_SMN_mRNA`/`dSMN_mRNA`) combines the nusinersen and
   risdiplam Hill effects as `EFFECT_NUS + EFFECT_RIS -
   EFFECT_NUS*EFFECT_RIS`, then scales by `SMN2_copies/2`. `$TABLE`'s own
   `E7I_out` (captured as `E7_inclusion`) independently recomputes the
   *same two Hill terms from the same underlying concentrations* but
   combines them as a plain sum, `EFFECT_NUS + EFFECT_RIS`, with **no**
   product-complement term and **no** `SMN2_copies` scaling. The reported
   "exon-7 inclusion" time course therefore does not equal the exon-7
   inclusion rate the model itself uses to drive SMN mRNA transcription
   — a pre-existing internal inconsistency between a reporting variable
   and the actual mechanism, not introduced by this refactor. Preserved
   as-is in the refactored file (same formula, renamed variables only).
2. **A dead-code toggle around the Zolgensma transgene contribution.**
   `double SMN_from_ZOL = k_tg_prot > 0 ? A_tg_mRNA : 0.0;` is computed in
   `$ODE` but never read by any `dxdt_`/`$TABLE`/`$CAPTURE` expression —
   `SMN_synthesis` reads `A_tg_mRNA` directly, not `SMN_from_ZOL`. The
   parameter `k_tg_prot` ("Transgene protein synthesis rate") therefore
   has no effect on any output through this expression; whatever its
   value, disease dynamics are unaffected. Preserved verbatim (renamed
   `K_PROT_ZOL`/`SMN_from_ZOL`) in the refactored file, not fixed.

**Fix upstream would be:** rename the three colliding `$TABLE` locals
(e.g. `CMAP_val`/`HFMSE_val`/`RULM_val`) and update their `capture`
lines to match, as done here; separately, decide whether `E7I_out`
should be redefined to match `E7I_current`'s actual combination rule, and
whether `SMN_from_ZOL`/`k_tg_prot` should be wired into `SMN_synthesis`
or removed as dead code — both are modelling decisions for the original
author, out of scope for this syntax-only build fix.

## 67. `prostate-cancer/pc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: incomplete `$INIT` annotations, `$CMT`+`$INIT` jointly redeclaring all 33 compartments, `$CAPTURE` wholly duplicating `$CMT`, and an undefined `self.trt_leup` reference; plus a real (non-blocking) Degarelix/Denosumab volume mix-up found while refactoring all 8 of this file's real compounds

Found while refactoring this file's PK/PD per its
`driver-patches/data/compound_perturbation_census.md` row (see the
corrected-compound-identity note below — the census's one row for this
file, "AR Signaling (DEG)", turned out to name neither a process nor the
right compound cleanly, and the file actually contains **eight** distinct,
independently-dosed real drugs, only one of which the census had a row
for at all). Reproduces identically from the untouched original via the
qspserver `mrgsolve_api` container's `/model_manifest`, unrelated to any
compound's own PK except where noted.

**Defect 1 (build-blocking): incomplete `$INIT @annotated` entries.**
`$INIT @annotated` requires exactly three colon-separated fields per line
(`name : value : description`); 16 of this file's 33 `$INIT` lines carry
only `name : value`, e.g. `Leup_c    : 0.0` with no description field:

```
Error: improper annotation format
 input: Leup_c    : 0.0
 context: parse annotated init block (INIT)
```

**Defect 2 (build-blocking, same family as #45–#66): `$CMT`+`$INIT`
jointly redeclare all 33 compartments.** Once defect 1 is worked around
(by adding a placeholder description to each incomplete line), mrgsolve
rejects the model outright: `Duplicated model names: LH T DHT ... Den_c`
(all 33 names). Same underlying mechanism already logged for
`diabetic-nephropathy`, `copd`, `essential-thrombocythemia`, and
`myotonic-dystrophy` (#61, #63, #64, #65) — declaring a compartment in
both `$CMT` and `$INIT` is not accepted by this mrgsolve build.

**Defect 3 (build-blocking): `$CAPTURE` wholly duplicates `$CMT`.** All 18
names in the original's `$CAPTURE @annotated` block (`PSA`, `T`, `DHT`,
`TC_p`, `TC_q`, `CRPC_frac`, `ARv7_frac`, `OC`, `OB`, `BMD`, `BoneMets`,
`AKT_act`, `Enz_c`, `Abi_c`, `Doc_c`, `Leup_c`, `Deg_c`, `AR_nuc`) are
already `$CMT` compartments — none is a `$TABLE`-derived quantity. mrgsolve
2.0.1 rejects this: `compartment should not be in $CAPTURE: PSA,T,DHT,...`.

**Defect 4 (build-blocking, in-scope for Leuprolide, fixed as part of the
refactor rather than logged as a generic defect): undefined
`self.trt_leup`.** `$MAIN` reads `(NEWIND <= 1 || self.trt_leup == 0) ?
1.0 : (1.0 + GnRH_flare*Flare_eff)*(1.0 - 0.97*Leup_suppress)` to compute
Leuprolide's GnRH-agonist effect, but `trt_leup` is never declared
anywhere in the file (not in `$PARAM`, not set via any R-side
`param()`/`idata` override) — `self` (the per-individual `databox`) has
no such member, so this fails to compile: `'class databox' has no member
named 'trt_leup'`. Because the file never compiles at all, there is no
baseline "original runtime behaviour" for this branch to preserve; the
formula branch alone (dropping the ternary) is what the model actually
needs to run, and it is also what the surrounding comment ("leuprolide:
flare then desensitize") describes. It is also behaviourally inert
whenever Leuprolide is absent — with `Leup_c = 0` and `Flare_eff = 0` (the
model's own defaults, never dosed otherwise except via Leuprolide's own
depot), the formula already reduces to the neutral value the guarded
branch would have returned (1.0), confirmed by this file's own "untreated"
scenario reproducing the expected unsuppressed LH/T trajectory. Handled as
an in-scope Leuprolide archetype decision, not a generic upstream note —
see `prostate-cancer/pc_refactor_notes.md` for the full disclosure and the
verification that reproduces it.

**Fixes for defects 1–4 applied directly to the delivered
`pc_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original: the incomplete
`$INIT` lines, and then the whole `$INIT` block, are replaced by
`<CMT>_0 = value;` assignments in `$MAIN` (same values, same
compartments); `$CAPTURE` no longer repeats any `$CMT` name (mrgsolve
reports every compartment's state via `/model_manifest`'s `outputPaths`
regardless of `$CAPTURE`, confirmed); and the `self.trt_leup` ternary is
replaced by its unconditional formula branch. All syntax-only and
non-numeric except defect 4, which is disclosed above as a considered,
in-scope design decision rather than a value change. Verified: with the
identical defect-1–4 fixes applied to an in-memory-only scratch copy of
the original (never to `pc_mrgsolve_model.R` itself, which still carries
all four defects exactly as written), all 16 disease-side shared outputs
matched the refactored file to floating-point-scale precision (max abs
diff ranging 0.0–0.06 against output magnitudes of order 1–1000) across
all 7 of the original's own dosing scenarios, plus 3 additional
single-compound checks constructed for the three compounds (Degarelix,
Relugolix, Denosumab) that no shipped scenario doses at all. See
`prostate-cancer/pc_refactor_notes.md` for the full per-scenario table.

**One additional finding, non-blocking (disclosed, not fixed): Degarelix's
own PK divides by the wrong volume parameter.** `$PARAM` declares
`V_Deg = 1000.0` ("Volume of distribution (L)") for Degarelix, but
`dxdt_Deg_c = kDeg_abs*Deg_sc/V_Den - kDeg_elim*Deg_c` actually divides by
`V_Den` (Denosumab's own volume, 3.0 L) instead — an apparent copy-paste
error. `V_Deg` is declared but never referenced anywhere in `$ODE`; the
real, in-effect Degarelix central volume is 3.0 L, not 1000 L, making the
model's simulated Degarelix plasma concentration roughly 333x higher than
a correctly-parameterized 1000 L volume would produce, and coincidentally
tying Degarelix's kinetics to whatever value Denosumab's own `V_Den` is
overridden to at run time. Preserved as-is (numerically) in the refactored
file: `V1_DEG = 3.0` (the value actually used, not the dead `1000.0`), now
under Degarelix's own independent parameter name, decoupling it from
Denosumab's `V1_DEN` structurally while reproducing the original's actual
number exactly. See `prostate-cancer/pc_refactor_notes.md`.

**Fix upstream would be:** add a description field to the 16 incomplete
`$INIT` lines (or drop `$INIT` and set values via `<CMT>_0` in `$MAIN`
instead, as done here); remove the 18 compartment names from `$CAPTURE`;
either declare `trt_leup` (e.g. as a `$PARAM` treatment-arm flag,
consistent with how each scenario already toggles a different drug's
dosing) and decide what runtime behaviour was actually intended, or drop
the ternary as done here; and change `V_Den` to `V_Deg` in
`dxdt_Deg_c`'s denominator.

## 68. `long-covid/pasc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 26 compartments, `$CAPTURE` duplicates 4 of them, and `$TABLE` re-declares `FSS`/`MoCA`/`SF36_PCS` a second time via self-referential `capture NAME = NAME;`

Found while refactoring nirmatrelvir/metformin/sertraline/low-dose-
naltrexone (LDN) PK/PD per this file's `driver-patches/data/
compound_perturbation_census.md` rows (two of which were mislabeled
process-description phrases rather than compound names — see the
corrected-identity note in `long-covid/pasc_refactor_notes.md`).
Reproduces identically from the untouched original via the qspserver
`mrgsolve_api` container's `/model_manifest`, unrelated to any compound's
own PK.

**Defect 1 (build-blocking): `$CMT` and `$INIT` jointly redeclare all 26
compartments.** The original declares all 26 compartments in `$CMT
@annotated`, then repeats every one of the same 26 names in a separate,
plain (non-`@annotated`) `$INIT` block to set initial values
(`V_GUT = 0.0`, ..., `C_LDN = 0.0`). mrgsolve 2.0.1 rejects this outright:
`Duplicated model names: V_GUT V_PLASMA V_RES V_AG IFN CD8_exh Auto_Ab
IL6 TNF Fibrin Ddimer BBB Microglia Serotonin AutNom ROS MitoDmg Lactate
A_nirm C_nirm A_met C_met A_sert C_sert A_LDN C_LDN`. Same underlying
mechanism already logged for `diabetic-nephropathy`, `copd`,
`essential-thrombocythemia`, `myotonic-dystrophy`, and `prostate-cancer`
(#61, #63, #64, #65, #67).

**Defect 2 (build-blocking): `$CAPTURE` duplicates 4 of the 26`$CMT`
names.** A separate, bare `$CAPTURE C_nirm C_met C_sert C_LDN` line lists
four names that are already `$CMT` compartments — mrgsolve rejects this:
`compartment should not be in $CAPTURE: C_nirm,C_met,C_sert,C_LDN`. Same
family as #48/#56/#57/#58 and `prostate-cancer`'s Defect 3 (#67).

**Defect 3 (build-blocking): `$TABLE` declares `double FSS = ...;` /
`double MoCA = ...;` / `double SF36_PCS = ...;` and then immediately
writes `capture FSS = FSS;` / `capture MoCA = MoCA;` /
`capture SF36_PCS = SF36_PCS;`.** Exactly the `spinal-muscular-atrophy`
mechanism (#66): a self-referential `capture NAME = NAME;` collides with
the `double NAME` already declared earlier in the same shared member
namespace: `error: redefinition of 'capture {anonymous}::FSS'` (`note:
'double {anonymous}::FSS' previously declared here`), and identically for
`MoCA` and `SF36_PCS`. Four other `$TABLE`-derived variables in the same
block (`VO2max_pct`/`capture VO2max`, `POTS_HR_delta`/`capture POTS_HR`,
`NfL`/`capture NfL_pg`, `dyspnea`/`capture mMRC`) do **not** hit this,
because their own local name already differs from the capture name it
feeds; only the three whose local name is spelled identically to its own
capture name collide. `SF36_PCS`'s own formula additionally reads `FSS`
(the colliding local), so the rename has to propagate through that one
downstream reference too.

**Confirmed upstream:** all three defects reproduce from the untouched
original alone via `POST /model_manifest`, no refactor changes involved.

**Fixes applied directly to the delivered
`pasc_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original, syntax-only and
non-numeric:
1. The `$INIT` block is dropped; its 26 assignments become `<CMT>_0 = value;`
   lines in `$MAIN` (identical values, identical compartments, comments
   preserved).
2. The standalone `$CAPTURE C_nirm C_met C_sert C_LDN` line is dropped
   entirely — mrgsolve reports every compartment's own state via
   `/model_manifest`'s `outputPaths` regardless of `$CAPTURE`, confirmed;
   the refactor separately (in-scope, see the refactor notes) gives each
   of these four its own named `C_<STEM>`/`EFFECT_<STEM>` `$TABLE`
   capture anyway.
3. The three colliding `$TABLE` locals are renamed `FSS_calc`/`MoCA_calc`/
   `SF36_PCS_calc` (matching the same pattern the original already used
   successfully for `VO2max_pct`/`POTS_HR_delta`), `SF36_PCS_calc`'s own
   formula reads `FSS_calc` in place of `FSS`, and all three `capture`
   lines read from the renamed local instead of re-declaring the capture
   name itself (`capture FSS = FSS_calc;`, etc.).

All three fixes were applied identically to an in-memory-only scratch
copy of the original (never to `pasc_mrgsolve_model.R` itself, which
still carries all three defects exactly as written) so it could build for
the `/run_simulation` comparison. Verified: with these fixes applied to
each side, both of the original's own tested scenarios (S1 untreated
natural history; S7 full combination, all four compounds dosed at their
own full one-year regimens) matched the refactored file's 12 shared
disease-side outputs and all 4 renamed PK outputs **exactly** (max abs
diff = 0.0 at every one of 366/370 timepoints). See
`long-covid/pasc_refactor_notes.md`.

**One additional finding, non-blocking, not fixed (unrelated to any
compound):** `$CMT` declares `V_GUT` (a "viral depot" compartment) but no
`dxdt_V_GUT` line exists anywhere in `$ODE` — mrgsolve issues a
non-fatal audit warning (`Missing differential equation(s): --| missing:
dxdt_V_GUT`) and treats it as a permanent zero-derivative state. `V_GUT`
is never dosed or read by anything else in the model either, so this is
inert dead state, not a behavioural bug — the viral kinetics submodel
starts directly from `V_PLASMA`, never from an oral/gut viral depot.
Preserved as-is (unchanged) in the refactored file.

**Fix upstream would be:** drop the separate `$INIT` block (or add a
description field per line and keep it, per whichever idiom the author
prefers) and set initial values via `<CMT>_0` in `$MAIN` instead, as done
here; remove the four compartment names from `$CAPTURE`; rename the three
colliding `$TABLE` locals (e.g. `FSS_calc`/`MoCA_calc`/`SF36_PCS_calc`)
and update their `capture` lines and internal cross-reference to match, as
done here; and either give `V_GUT` a `dxdt_V_GUT` line or remove the
unused compartment.

---

## 69. `gaucher-disease/gcd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 26 compartments, `$CAPTURE` duplicates all 21 of them, and two `$PARAM` baseline names collide with mrgsolve's own auto-reserved `<CMT>_0` symbols

Found while refactoring ERT (imiglucerase/velaglucerase class)/eliglustat
(ELIS)/miglustat (MIGS)/venglustat (VENG) PK/PD per this file's
`driver-patches/data/compound_perturbation_census.md` rows. Reproduces
identically from the untouched original via the qspserver `mrgsolve_api`
container's `/model_manifest`, unrelated to any compound's own PK.

**Defect 1 (build-blocking): `$CMT` and `$INIT` jointly redeclare all 26
compartments.** The original declares all 26 compartments bare in `$CMT`,
then repeats every one of the same 26 names in a separate `$INIT` block
to set initial values (`ERT_C = 0, ERT_T = 0`, ..., `NFKB = 1.0`). mrgsolve
2.0.1 rejects this outright: `Duplicated model names: ERT_C ERT_T ELIS_GUT
ELIS_C MIGS_GUT MIGS_C VENG_GUT VENG_C GBA GC_MAC GC_SP GC_LV GC_BM GL1
LYSOGL1 CHITR FERRIT SV LV HGB PLT BMD OC OB IL6 NFKB`. Same underlying
mechanism already logged for `diabetic-nephropathy`, `copd`,
`essential-thrombocythemia`, `myotonic-dystrophy`, `prostate-cancer`, and
`long-covid` (#61, #63, #64, #65, #67, #68).

**Defect 2 (build-blocking): `$CAPTURE` duplicates all 21 of the
non-drug-PK `$CMT` names it lists.** The original's `$CAPTURE` block
lists `GBA GBA_PCT_NORMAL GC_MAC GC_SP GC_LV GC_BM GC_TOTAL GL1 LYSOGL1
CHITR FERRIT SV LV HGB PLT BMD OC OB IL6 NFKB ERT_C ERT_T ELIS_C MIGS_C
VENG_C ERT_Cplasma SRT_GCS_inh` — of these 25 names, 21 are `$CMT`
compartments (only `GBA_PCT_NORMAL`, `GC_TOTAL`, `ERT_Cplasma`, and
`SRT_GCS_inh` are genuine `$TABLE`-derived doubles). mrgsolve rejects
this: `compartment should not be in $CAPTURE: GBA,GC_MAC,GC_SP,GC_LV,
GC_BM,GL1,LYSOGL1,CHITR,FERRIT,SV,LV,HGB,PLT,BMD,OC,OB,IL6,NFKB,ERT_C,
ERT_T,ELIS_C,MIGS_C,VENG_C`. Same family as #48/#56/#57/#58/#61/#63/#67/#68.

**Defect 3 (build-blocking, surfaces only once Defect 1 is worked around
with the modern `<CMT>_0` idiom): two `$PARAM` baseline names collide
with mrgsolve's own auto-generated `<CMT>_0` initial-value symbols.**
`GL1_0` (baseline plasma GL-1) and `IL6_0` (baseline IL-6 composite
cytokine) are ordinary `$PARAM` values used nowhere in `$ODE`/`$TABLE`
(dead parameters, apparently left over from an earlier authoring pass),
but their names are spelled identically to the auto-generated init symbol
for the identically-named compartments `GL1` and `IL6`. Same incidental-
collision pattern already logged for `sepsis` (#30, `IL6_0`/`IL10_0`/
`PAI1_0` vs. compartments `IL6`/`IL10`/`PAI1` — the very same `IL6_0`
name, in fact), `chronic-lymphocytic-leukemia` (#44, `ALC_0`), `copd`
(#63), and `myotonic-dystrophy` (#65).

**Confirmed upstream:** all three defects reproduce from the untouched
original alone via `POST /model_manifest`, no refactor changes involved.

**Fixes applied directly to the delivered `gcd_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for a
non-compiling original, syntax-only and non-numeric:
1. The `$INIT` block is dropped; its 26 assignments become `<CMT>_0 = value;`
   lines in `$MAIN` (identical values, identical compartments).
2. `$CAPTURE` is trimmed to only the four genuine `$TABLE`-derived
   doubles (`GC_TOTAL`, `GBA_PCT_NORMAL`, `ERT_Cplasma`, `SRT_GCS_inh`)
   plus the refactor's own new named `C_<STEM>`/`EFFECT_<STEM>` doubles
   (in scope anyway, see `gcd_refactor_notes.md`) — every renamed
   compartment (`CENT_ERT`, `PERI_ERT`, `GUT_ELIS`, `CENT_ELIS`, etc.) is
   still reported automatically via `/model_manifest`'s `outputPaths`
   without being listed in `$CAPTURE`, confirmed.
3. `GL1_0` -> `GL1_SS0` and `IL6_0` -> `IL6_SS0` (all usages updated,
   including the R-side virtual-population tibble column that also
   carried the `GL1_0` name); both parameters remain unused dead
   parameters, so this is a pure rename with no numeric effect.

All three fixes were applied identically to an in-memory-only scratch
copy of the original (never to `gcd_mrgsolve_model.R` itself, which still
carries all three defects exactly as written) so it could build for the
`/run_simulation` comparison. Verified: with these fixes applied to each
side, all six of the original's own tested scenarios/dosing regimens
(S1 natural history; S2 imiglucerase; S4/S5 eliglustat EM/PM; S6
ERT+eliglustat combination; plus a bespoke miglustat/venglustat dosing
run built from the original's own unused `make_oral_doses()` helper and
its own `DOSE_MIGS`/`DOSE_VENG` defaults, since no named scenario in the
original doses either compound — see `gcd_refactor_notes.md`) matched the
refactored file's shared disease-side outputs and all 8 renamed/mapped PK
outputs **exactly** (max abs diff = 0.0 at every timepoint, 366-368
points per scenario). See `gaucher-disease/gcd_refactor_notes.md`.

**Fix upstream would be:** drop the separate `$INIT` block and set
initial values via `<CMT>_0` in `$MAIN` instead, as done here; remove the
21 compartment names from `$CAPTURE`, keeping only the genuine `$TABLE`
doubles; and rename `GL1_0`/`IL6_0` away from the `<compartment>_0`
pattern (or simply delete them, since neither is read anywhere).

## 70. `pseudogout/cppd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` duplicates all 12 disease-side `$CMT` compartment names

Found while verifying an anakinra (ANA) / colchicine (COLCH) / indomethacin
(INDO) / prednisolone (PRED) PK/effect-interface refactor
(`cppd_mrgsolve_model_refactored.R`); reproduces identically from the
untouched original via the qspserver `mrgsolve_api` container, before any
rename was applied.

**The defect.** `$CAPTURE` lists, alongside four genuine `$MAIN` locals
(`GR_occ E_COLCH_NLRP3 E_INDO_COX E_PRED_NFKB E_ANA`) and four plasma
concentrations (`Cp_COLCH Cp_INDO Cp_PRED Cp_ANA`, also `$MAIN` locals),
twelve names that are already `$CMT` compartments: `PPi_ext Cryst_cart
Cryst_SF NLRP3_act IL1b Neutrophil IL6 PGE2 CRP PainVAS CartInteg
LipoxA4`. mrgsolve 2.0.1 rejects this outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: PPi_ext,Cryst_cart,Cryst_SF,NLRP3_act,IL1b,Neutrophil,IL6,PGE2,CRP,PainVAS,CartInteg,LipoxA4
```

Confirmed on the untouched original alone via `POST /model_manifest` (no
renames involved). Same defect class already logged as issues #56, #61,
#63, #64, #65, #67 (and others), here in a different file.

**Fix applied, syntax-only, to the delivered `cppd_mrgsolve_model_refactored.R`
only** (never to `cppd_mrgsolve_model.R`, which still carries the defect
exactly as written): the twelve compartment names were removed from
`$CAPTURE`, leaving only the genuine `$MAIN`-local doubles (renamed to
`C_COLCH C_INDO C_PRED C_ANA EFFECT_COLCH_NLRP3 EFFECT_COLCH_NEUT
EFFECT_INDO EFFECT_PRED EFFECT_PRED_NFKB EFFECT_ANA` per the refactor).
Every compartment remains fully available in simulation output regardless
of `$CAPTURE` membership -- confirmed by diffing `/model_manifest`'s
`outputPaths` before/after, still lists all twelve under their own names
(renamed ones under their new PK names, the twelve disease-side ones
unchanged). This is purely a build-compatibility fix, not a value change:
verified by running all seven of the original file's own dosing scenarios
(untreated; colchicine loading; indomethacin; prednisolone taper;
anakinra; colchicine+indomethacin combo; 180-day colchicine prophylaxis)
through both the `$CAPTURE`-patched-only original and the fully renamed
`cppd_mrgsolve_model_refactored.R` via the qspserver `mrgsolve_api` --
every shared output matched exactly (max abs diff 0.0) in all seven. See
`pseudogout/cppd_refactor_notes.md` for the full disclosure.

**Fix upstream would be:** remove the twelve compartment names from
`$CAPTURE`, keeping only the eight genuine `$MAIN`-local doubles.

## 71. `pompe-disease/pompe_mrgsolve_model.R` does not compile under mrgsolve 2.0.1, and its `$TABLE` capture list is silently non-functional even when patched to compile

Confirmed upstream via `POST /model_manifest` and `POST /run_simulation` on
the untouched original alone (`http://localhost:8007`, qspserver
`mrgsolve_api` container), while refactoring this file's four Redirect
concentration compounds (Alglucosidase alfa, Avalglucosidase alfa,
Cipaglucosidase alfa, Miglustat). Two independent, layered defects:

**Defect 1: 11 `$PARAM @annotated` lines have no description field.**
mrgsolve 2.0.1's annotated-parameter parser requires the three-part
`name : value : description` form; these lines only have two parts
(`Q_AVAL    :   2.0`, `V2_AVAL   :   4.5`, `CL_CIPA   :  22`,
`V1_CIPA   :   3.6`, `Q_CIPA    :   1.8`, `V2_CIPA   :   4.5`,
`M6P_CIPA  :  25`, `DIAPH_LOSS:   0.0010`, `DIAPH_GAIN:   0.0006`,
`SMWT_MAX  : 600`, `SMWT_MIN  : 100`), which fails to build at all:

```
Error: improper annotation format
 input: Q_AVAL    :   2.0
 context: parse annotated parameter block (PARAM)
Execution halted
```

Seven of these eleven belong to Avalglucosidase alfa's and Cipaglucosidase
alfa's own PK blocks (in scope for the refactor); the remaining four
(`DIAPH_LOSS`, `DIAPH_GAIN`, `SMWT_MAX`, `SMWT_MIN`) are disease-side and
out of scope, but the whole file cannot compile at all until every one of
the eleven has a description added, so all eleven were fixed together as
one build-compat pass.

**Defect 2 (a new pattern, not previously logged): `$TABLE` ends in 13
bare `capture NAME;` lines — one bare identifier per line, each terminated
with a semicolon, no `$CAPTURE` header anywhere in the file.** Unlike the
`$CAPTURE`-duplicates-a-compartment-name defect logged repeatedly elsewhere
in this file (e.g. #30, #34, #41, #46, #57, #58, #64), and unlike
`bronchiectasis`'s bare `capture A B C` grouped-line variant (#46, which
fails to *compile*), this variant *compiles without any error* under this
mrgsolve build — but is silently non-functional at runtime: none of the 13
named quantities (`Cp_alglu`, `Cp_aval`, `Cp_cipa`, `Cp_mig`, `Cp_rtx`,
`SMWT`, `VENT_RISK`, `SF36_PCS`, `NTproBNP`, `EF_LV`, `CK`,
`tissue_supply`, `ada_block`) is actually retrievable as simulation output:

```
Cp_alglu is not a compartment or captured item
Cp_aval is not a compartment or captured item
...
Error in h(simpleError(msg, call)) :
  error in evaluating the argument 'x' in selecting a method for function
  'as.data.frame': invalid item in requested output
```

Confirmed the fix by converting the same 13 names to a real `$CAPTURE`
header (`$CAPTURE Cp_alglu Cp_aval Cp_cipa Cp_mig Cp_rtx SMWT VENT_RISK
SF36_PCS NTproBNP EF_LV CK tissue_supply ada_block`, no symbol added,
removed, or reordered): all 13 then appear correctly in `/run_simulation`
output, with values matching the model's own calibration-note expectations
(e.g. `SMWT` ~300-560 m, `Cp_alglu` rising and decaying with each q2w
dose). This means every diagnostic quantity this file's own author
intended to expose via `$TABLE` — plasma concentrations, 6MWT, ventilator
risk, SF-36, NT-proBNP, EF, CK, tissue enzyme supply, ADA neutralization —
has been silently unavailable to any caller of this model, with no error
ever surfacing to say so.

Both fixes are syntax-only and non-numeric, applied directly to the
delivered `pompe_mrgsolve_model_refactored.R` per the guide's settled
policy for a non-compiling original (`FORK_WORKFLOW_GUIDE.md`, "When the
original doesn't compile at all"): confirmed via `POST /run_simulation`
that all four refactored compounds' scenarios (no dosing; Alglucosidase
alfa alone; Avalglucosidase alfa alone; Cipaglucosidase alfa + Miglustat
combined) produce identical output (max abs diff 0.0 for every shared
quantity) between the identically-patched original and the fully renamed
refactored model. See `pompe-disease/pompe_refactor_notes.md` for the full
account; the tracked `pompe_mrgsolve_model.R` itself is untouched and still
carries both defects exactly as written.

**Fix upstream would be:** add a description field to the 11 incomplete
`$PARAM @annotated` lines, and replace the 13 bare `capture NAME;` lines
with a single `$CAPTURE` block (or several) listing the same 13 names.

## 72. `von-willebrand-disease/vwd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` duplicates 12 compartment names, and `$MAIN` sets 14 bare compartment names directly instead of using the `<CMT>_0` initial-value idiom; plus a DDAVP-specific fractional-`pow()` NaN fragility in the file's own repeat-dosing scenario

Found while refactoring this file's 4 compounds (DDAVP, recombinant VWF,
plasma-derived VWF/FVIII concentrate, tranexamic acid) per
`FORK_WORKFLOW_GUIDE.md` Part 2. Two build-blocking defects, both
unrelated to any compound's own PK:

**Defect 1: `$CAPTURE` duplicates 12 of the file's own `$CMT` names.**
`$CAPTURE VWF_AG_TOTAL VWF_RCO_TOTAL HMWM_EFFECTIVE ADAMTS13_ACT
FVIII_C_TOTAL PLT_COUNT BLEED_SCORE MENS_LOSS GI_LOSS HB NA_SERUM
THROMB_RISK WPB_STORE VWFPP DDAVP_CP RVWF_CONC PDVWF_CONC TXA_CP` lists
17 names, of which 12 (`ADAMTS13_ACT`, `PLT_COUNT`, `BLEED_SCORE`,
`MENS_LOSS`, `GI_LOSS`, `HB`, `NA_SERUM`, `THROMB_RISK`, `WPB_STORE`,
`VWFPP`, `DDAVP_CP`, `TXA_CP`) are also `$CMT` compartments (this build
already reports every compartment automatically, so listing one again in
`$CAPTURE` is rejected):

```
Error in validObject(.Object) :
  invalid class 'mrgmod' object: compartment should not be in $CAPTURE:
  ADAMTS13_ACT,PLT_COUNT,BLEED_SCORE,MENS_LOSS,GI_LOSS,HB,NA_SERUM,
  THROMB_RISK,WPB_STORE,VWFPP,DDAVP_CP,TXA_CP
```

Same category as #30/#34/#41/#46/#57/#58/#64/#70. The remaining 5 names
(`VWF_AG_TOTAL`, `VWF_RCO_TOTAL`, `HMWM_EFFECTIVE`, `FVIII_C_TOTAL`,
`RVWF_CONC`, `PDVWF_CONC` — 6, not `$CMT` names) are genuine `$ODE`-local
doubles and capture correctly once the 12 duplicates are removed.

**Defect 2 (once Defect 1 is worked around): `$MAIN` assigns 14 bare
compartment names directly inside `if (NEWIND <= 1) { ... }` instead of
using this build's `<CMT>_0` initial-value idiom:**

```
213:15: error: assignment of read-only reference 'WPB_STORE'
  213 |   WPB_STORE   = WPB0;
      |   ~~~~~~~~~~~~^~~~~~
```

(and 13 more identical errors for `VWFPP`, `VWF_AG`, `VWF_RCO`, `HMWM`,
`ADAMTS13_ACT`, `FVIII_C`, `PLT_COUNT`, `BLEED_SCORE`, `MENS_LOSS`,
`GI_LOSS`, `HB`, `NA_SERUM`, `THROMB_RISK`). Same category as the
`TIME`/`_init_<CMT>` idiom defects noted in `FORK_WORKFLOW_GUIDE.md`; the
fix is renaming each target to `<CMT>_0` (e.g. `WPB_STORE_0 = WPB0;`),
right-hand sides unchanged.

Both fixes are syntax-only and non-numeric, applied directly to the
delivered `vwd_mrgsolve_model_refactored.R` per the guide's settled policy
for a non-compiling original (`FORK_WORKFLOW_GUIDE.md`, "When the original
doesn't compile at all"); the tracked `vwd_mrgsolve_model.R` itself is
untouched and still carries both defects exactly as written.

**A third finding, DDAVP-specific and inside this refactor's own scope
(fixed as part of the refactor itself, not logged as a separate defect
per the guide's point 4, but noted here for completeness):** the file's
own `EFFECT_DDAVP`-equivalent Hill term (`DDAVP_CP` in the original, via
`pow(DDAVP_CP, HILL_DDAVP)` with `HILL_DDAVP = 1.3`, a genuinely
fractional exponent) NaNs when the decaying DDAVP concentration
numerically undershoots zero at ~1e-10 scale during adaptive-step
integration. Confirmed present in the *original* file too (patched only
for the two defects above): its own `3_Type1_DDAVP_Intranasal_Repeat`
scenario (300 mcg q12h x6 doses, Type-1 genotype) NaNs at t=167h. The
refactored file's algebraically-equivalent rewrite (`CL_DDAVP/V1_DDAVP`
in place of `KE_DDAVP`, plus one extra division for `C_DDAVP =
CENT_DDAVP/V1_DDAVP`) shifts the same floating-point knife-edge earlier,
to t=160h — not a structural error (see
`von-willebrand-disease/vwd_refactor_notes.md` for the full verification
proving numeric equivalence up to that point), but a real pre-existing
fragility. Guarded in the delivered refactored file with
`fmax(C_DDAVP, 0.0)` before the `pow()`; the tracked original still NaNs
past t=167h in this scenario, unfixed.

**Fix upstream would be:** remove the 12 duplicate names from `$CAPTURE`;
rename the 14 `$MAIN` assignment targets to `<CMT>_0`; and guard
`DDAVP_CP` (or replace `pow(DDAVP_CP, HILL_DDAVP)` with a safe
non-negative variant) before the fractional `pow()`.

## 73. `waldenstrom-macroglobulinemia/wm_mrgsolve_model.R` does not compile under mrgsolve 2.0.1 (`$ODE` writes directly to a read-only compartment state), and separately its `_F(CMT) = value;` bioavailability idiom crashes the solver at runtime even once patched to compile

Found while refactoring all four of this file's compounds (Ibrutinib,
Rituximab, Venetoclax, Zanubrutinib — `wm_mrgsolve_model_refactored.R`).
Two independent defects, confirmed to reproduce from the untouched
original alone via the qspserver `mrgsolve_api` container; both fixed
syntax-only, disclosed, in the delivered `_refactored.R` per the guide's
settled policy.

**1. `$ODE` assigns directly to a `$CMT` state to hold it at an algebraic
value:**
```
double total_tumor = LPC + PC;
BMInf = fmin(total_tumor / KMAX_BM, 1.0);
dxdt_BMInf = 0; // algebraic (set above)
```
mrgsolve 2.0.1 passes every `$CMT` state into the generated `$ODE`
function as a `const double&`, so writing to `BMInf` itself (not just
`dxdt_BMInf`) is a hard compile error:
```
424:9: error: assignment of read-only reference 'BMInf'
  424 |   BMInf = fmin(total_tumor / KMAX_BM, 1.0);
      |   ~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```
Same defect class as issue #36's `FLUID_EX`/`PR_FRAC` in
`age-related-macular-degeneration`, but not a no-op here: unlike that
case's dead clamp, this line is the *only* place `BMInf` is ever given a
non-zero value — `dxdt_BMInf` is always `0`, so the intent was plainly an
"instant algebraic" compartment, overwritten every step, which is exactly
what this mrgsolve build disallows. **Fix (in the delivered
`_refactored.R` only):** `BMInf` was removed from `$CMT` and its `$MAIN`
init line, and demoted to a plain `double BMInf = fmin(...)` local inside
`$ODE` (no `dxdt_` line needed, since it is no longer a state); the one
downstream `$ODE` read (the haemoglobin-suppression term) is unaffected
since it runs later in the same block. `$TABLE`'s `BMInf_pct` is
recomputed directly from `LPC`/`PC` with the identical formula (which
`$TABLE` can read as compartments), giving the same numeric value.
Confirmed via `POST /run_simulation`: with both this fix and defect 2's
fix applied identically to scratch copies of the original and the
refactored file, every scenario matches at max abs diff 0.0 (see
`waldenstrom-macroglobulinemia/wm_refactor_notes.md`).

**2. `_F(CMT) = value;` (used for all three oral compounds' bioavailability
— ibrutinib, zanubrutinib, venetoclax) compiles cleanly (`POST
/model_manifest` succeeds even with defect 1 unpatched, since manifest
introspection never executes `$MAIN`/`$ODE`) but crashes `POST
/run_simulation` outright, for *every* scenario tested, dosed or not:**
```
"stderr": "munmap_chunk(): invalid pointer\n"
"returncode": -6
```
Confirmed by bisection: removing the three `_F(...) = ...;` lines alone
(nothing else) eliminates the crash; a further minimal repro (a bare
2-compartment depot+central model with only `_F(GUT) = F1;` added, no
dosing at all) reproduces the identical signal-6 crash on its own, so
this is unrelated to dosing, to any of the four compounds' own math, and
to defect 1. Replacing `_F(GUT_IBR) = F_IBR;` (etc.) with the modern
`F_GUT_IBR = F_IBR;` idiom resolves it and applies bioavailability
identically — confirmed by comparing steady-state central-compartment
concentration with and without the substitution on the same minimal
repro (scales by exactly the `F` value either way, when the crash is
avoided). This appears to be the first time `_F()` syntax has been
logged in this corpus; worth checking for in any other file using it.

**Fix upstream would be:** delete the direct `BMInf = ...;` assignment
line (keep only the algebraic derivation, recomputed wherever `BMInf` is
read) and remove it from `$CMT`/`$MAIN`; replace all three
`_F(CMT) = value;` lines with `F_CMT = value;`.

## 74. `takayasu-arteritis/ta_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$PARAM`/`$CMT` C-style block comments break their non-annotated parsers, `$CMT`+`$INIT` jointly redeclare all 24 compartments, and `$TABLE`'s three multi-name `capture` lines are not valid capture syntax; plus a shared, pre-existing prednisolone `pow()` NaN fragility in every one of the file's own combination scenarios

Found while refactoring this file's four compounds (prednisone/
prednisolone `PRED`, tocilizumab `TCZ`, methotrexate `MTX`, infliximab
`IFX`) per `FORK_WORKFLOW_GUIDE.md` Part 2. Three independent
build-blocking defects, all unrelated to any compound's own PK, plus one
runtime fragility inside PRED's own PD term (in-scope, fixed as part of
the refactor itself per the guide's point 4 rather than merely logged —
see below).

**Defect 1: `$PARAM` and `$CMT` (both non-annotated blocks) use C-style
`/* ... */` block comments, which this build's non-annotated-block
parsers cannot handle.** `$PARAM`'s parser converts the raw block text to
an R `list(...)` call; the very first comment breaks the R parse:
```
Error in parse(text = paste0("list(", x, ")")) :
  <text>:1:8: unexpected '/'
1: list(  /
           ^
```
`$CMT`'s parser instead strips comment *contents* but leaves the `/*`/`*/`
delimiters as stray bare tokens, which are then read as (invalid,
duplicated) compartment names:
```
invalid class 'mrgmod' object: 1: Duplicated model names: /* */ /* */ ...
invalid class 'mrgmod' object: 2: Names without leading alpha character: /* */ ...
```
Confirmed by removing every `/* ... */` span (comments only, no name or
value touched) from both blocks, which clears both errors.

**Defect 2 (once Defect 1 is worked around): `$CMT` (bare, unannotated)
and `$INIT` (bare `name = value` block) jointly redeclare all 24
compartment names** — same category as #30/#34/#41/#46/#57/#58/#64/#70/
#72:
```
invalid class 'mrgmod' object: Duplicated model names: PRED_GUT PRED_C
PRED_P PRED_EFF TCZ_SC TCZ_C TCZ_P MTX_GUT MTX_C MTX_PG IFX_C IFX_P IL6
sIL6R IL6_cmplx TNF TH1 TH17 TREG VWI ST CRP PET VWT
```
Fix: move `$INIT`'s 24 initial values into a `$MAIN` block using the
standard `<CMT>_0 = value;` idiom, right-hand sides unchanged.

**Defect 3 (once Defects 1-2 are worked around): `$TABLE`'s three bare
`capture NAME1 NAME2 NAME3 ...;` statements (space-separated, no commas)
are not valid multi-name capture syntax under this mrgsolve build** —
they compile as ordinary C++ statements, and the C++ compiler chokes on
the second bare name:
```
723:19: error: expected ';' before 'NIH_SCORE'
  723 |   capture ESR_now NIH_SCORE ITAS_SCORE RESPONSE_FLAG
      |                   ^~~~~~~~~
```
Fix: consolidate the three `capture ...` lines into one ordinary
`$CAPTURE ESR_now NIH_SCORE ITAS_SCORE RESPONSE_FLAG CP_PRED CP_TCZ
CP_MTX CP_IFX Drug_inh_VWI Occ_TCZ Inh_PRED Inh_MTX Inh_IFX` block (same
14 names, same values).

All three fixes are syntax-only and non-numeric, applied directly to the
delivered `ta_mrgsolve_model_refactored.R` per the guide's settled policy
for a non-compiling original (`FORK_WORKFLOW_GUIDE.md`, "When the
original doesn't compile at all"); the tracked `ta_mrgsolve_model.R`
itself is untouched and still carries all three defects exactly as
written. Confirmed via the qspserver `mrgsolve_api` container
(`http://localhost:8007`) that an identically-patched scratch copy of the
original and the fully renamed refactored model produce max abs diff 0.0
across all 28 shared outputs (every `$CMT` state plus every `$TABLE`
capture) for all three of the file's own two-drug combination scenarios
(Scenario 3, prednisone+methotrexate; Scenario 4, prednisone+tocilizumab,
the TAKT regimen; Scenario 5, prednisone+infliximab), each run over the
full 365-day/8760h horizon the scenario itself specifies — see
`takayasu-arteritis/ta_refactor_notes.md`.

**A fourth finding, PRED-specific and inside this refactor's own scope
(fixed as part of the refactor itself, not logged as a separate defect,
per the guide's point 4, but noted here for completeness — same pattern
as issue #72's DDAVP fragility):** the original's prednisolone Hill term
(`Inh_PRED`, via `pow(PRED_EFF, hill_PRED)` with `hill_PRED = 1.5`, a
genuinely fractional exponent) NaNs once the tapering prednisone dose
decays `PRED_C`/`PRED_EFF` down to numerical noise level between doses
and the adaptive-step solver momentarily undershoots to a tiny *negative*
value (confirmed: `PRED_C = -8.5e-11` at t=2288h, immediately followed by
`NaN` in every compartment at t=2292h). Because every one of the file's
own named scenarios that includes prednisone (2 through 5 — i.e. every
combination scenario) runs the same taper, this reproduces identically
in the *original* file too (patched only for Defects 1-3 above) at
t=2292h (Scenario 4), t=1568h (Scenario 3), and t=1576h (Scenario 5) —
all well inside the intended 8760h horizon. Guarded in the delivered
refactored file with `fmax(C_PRED, 0.0)` (and, for consistency/
robustness against the same fragility, `fmax(...)` on `C_TCZ`/`C_MTX`/
`C_IFX` too — `GAMMA_TCZ=1.2` and `GAMMA_IFX=1.3` are likewise
non-integer and therefore equally NaN-exposed in principle, even though
none of the tested scenarios actually drove `CENT_TCZ`/`CENT_IFX`
negative before PRED's own blowup pre-empted them); the tracked original
still NaNs at the same three time points, unfixed. Verified: with the
guard applied, the refactored model matches the original at max abs diff
0.0 for every point before each scenario's blowup time, and then — unlike
the original — continues to produce finite, physically sensible output
(`NIH_SCORE`, `CRP`, etc. carrying on smoothly) through the full 365-day
horizon in all three scenarios.

**Fix upstream would be:** strip the `/* ... */` comments from `$PARAM`
and `$CMT`; rename the 24 `$INIT` targets to `<CMT>_0` inside `$MAIN`;
replace the three `capture ...` lines with one `$CAPTURE` block; and
guard `PRED_EFF` (or replace `pow(PRED_EFF, hill_PRED)` with a safe
non-negative variant) before the fractional `pow()`.

## 75. `cytokine-release-syndrome/crs_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` lists ten compartment names

Found while refactoring this file's five compounds (anakinra `ANA`,
dexamethasone `DEX`, ruxolitinib `RUX`, siltuximab `SILT`, tocilizumab
`TOCI`) per `FORK_WORKFLOW_GUIDE.md` Part 2. Single, isolated defect —
same category as #56/#57/#70.

`$CAPTURE` lists ten names that are already `$CMT` compartments
(`CRS_SEV`, `ORGAN_DMG`, `CART_ACT`, `TUMOR`, `VASC_PERM`, `ENDO_ACT`,
`MAC_ACT`, `STAT3`, `GMCSF`, `IL2`) alongside `$TABLE`-computed aliases
that are fine. mrgsolve 2.0.1 validates the compiled model object and
rejects this before any simulation can run:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE: CRS_SEV,ORGAN_DMG,CART_ACT,TUMOR,VASC_PERM,ENDO_ACT,MAC_ACT,STAT3,GMCSF,IL2
```

**Confirmed upstream:** reproduced via `POST /model_manifest` on the
untouched original file's own DSL, no changes involved.

**Why this matters here:** this file's mandatory refactor-verification
step (`crs_mrgsolve_model_refactored.R`, all five compounds' PK/PD blocks)
requires actually building and running both the original and the
refactored model through the qspserver `mrgsolve_api`. Per the guide's
settled policy for this situation, the fix applied is syntax-only and
non-numeric: **the ten compartment names were removed from `$CAPTURE`**
(mrgsolve still reports every compartment's state regardless of whether
it also appears in `$CAPTURE` — confirmed by checking `/model_manifest`'s
`outputPaths`, which lists all ten compartments both before and after the
fix, so nothing about what is reported changes, only what compiles). This
fix was applied directly to the delivered `crs_mrgsolve_model_refactored.R`
(and, for comparison purposes only, to an in-memory-only copy of the
original sent to the API — never to the checked-in
`crs_mrgsolve_model.R`, which is left untouched, still carrying the
defect exactly as written). See `crs_refactor_notes.md` for full
disclosure and the verification result this fix enabled (exact match,
max abs diff 0.0, across all five compounds' own/bespoke dosing
scenarios).

**Fix upstream would be:** remove the ten duplicated names from
`$CAPTURE`; nothing else about the model needs to change.

## 76. `cervical-cancer/cc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 19 compartments, nine `$ODE` clamps write directly to read-only compartment states, and three `$PARAM` baseline names collide with mrgsolve's own auto-reserved `<CMT>_0` symbols; plus a shared, pre-existing tumor-volume `log()` NaN fragility once RT is active

Found while refactoring this file's five compounds (Bevacizumab `BEV`,
Cisplatin `CIS`, Paclitaxel `PAC`, Pembrolizumab `PEM`, TV-ADC/tisotumab
vedotin `TVADC`) per `FORK_WORKFLOW_GUIDE.md` Part 2 and the existing rows
in `driver-patches/data/compound_perturbation_census.md`. All defects
reproduce identically from the untouched original via the qspserver
`mrgsolve_api` container, unrelated to any compound's own PK math.

**Defect 1 (build-blocking): `$CMT` (`@annotated`) and `$INIT` jointly
redeclare all 19 compartments.** The original names and describes all 19
compartments in `$CMT @annotated`, then repeats every one of the same 19
names in a separate `$INIT` block to set starting values (`CIS_C1=0,
CIS_C2=0, ..., PDL1_exp=1.0`). mrgsolve 2.0.1 treats `$INIT` as its own
compartment-declaring block rather than a companion to `$CMT`, so this is
rejected outright, same underlying mechanism as #34/#36/#58/#59/#67/#68/
#69 and others:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: CIS_C1 CIS_C2 PAC_C1 PAC_C2 BEV_C1 BEV_C2 PEMBRO_C1 PEMBRO_C2 TV_ADC_C1 TV_ADC_C2 MMAE_free VEGF Pt_DNA RT_SF TV SCCAg HPVload CD8T PDL1_exp
```

**Defect 2 (build-blocking, surfaces only once Defect 1 is worked around):
nine `$ODE` lines clamp a compartment state directly instead of only
writing `dxdt_*`.** `if(MMAE_free < 0) MMAE_free = 0;`, `if(VEGF < 0) VEGF
= 0;`, `if(Pt_DNA < 0) Pt_DNA = 0;`, `if(RT_SF < 0) RT_SF = 0;`, `if(CD8T <
0) CD8T = 0;`, `if(PDL1_exp < 0.1) PDL1_exp = 0.1;`, `if(TV < 0.01) TV =
0.01;`, `if(SCCAg < 0.1) SCCAg = 0.1;`, and `if(HPVload < 0) HPVload = 0;`.
mrgsolve 2.0.1 passes every `$CMT` state into the generated `$ODE`
function as a `const double&`, so a direct write is a hard compile error,
not a warning:

```
343:29: error: assignment of read-only reference 'MMAE_free'
346:19: error: assignment of read-only reference 'VEGF'
348:23: error: assignment of read-only reference 'Pt_DNA'
353:21: error: assignment of read-only reference 'RT_SF'
357:19: error: assignment of read-only reference 'CD8T'
359:29: error: assignment of read-only reference 'PDL1_exp'
371:18: error: assignment of read-only reference 'TV'
374:23: error: assignment of read-only reference 'SCCAg'
377:25: error: assignment of read-only reference 'HPVload'
```

Confirmed via an isolated two-line minimal reproduction (`$CMT A; $ODE
dxdt_A=-CL/V*A; if(A<0) A=0;`) that this is a general mrgsolve 2.0.1
property, not specific to this file: it fails identically with or without
a separate `$MAIN`-set initial condition for `A`. Same defect class as
#36 (`FLUID_EX`/`PR_FRAC`, age-related-macular-degeneration) and #73
(`BMInf`, waldenstrom-macroglobulinemia): only `dxdt_*` feeds mrgsolve's
integrator, so a bare compartment reassignment inside `$ODE` was already a
no-op in any mrgsolve version — confirmed here too (see Defect 4 below;
removing all nine clamps and re-verifying against the original with the
same nine clamps patched out on both sides gives an exact, 0.0-diff match
everywhere the trajectory is finite).

**Defect 3 (build-blocking, surfaces only once Defects 1-2 are worked
around): three `$PARAM` baseline names collide with mrgsolve's own
auto-generated `<CMT>_0` initial-value symbols.** `SCCAg_0` (baseline
SCC-Ag), `HPVload_0` (baseline HPV load), and `CD8T_0` (baseline CD8+ T)
are ordinary `$PARAM` values, each read nowhere in `$ODE`/`$MAIN`/`$TABLE`
(dead parameters — the actual `$INIT` starting values for `SCCAg`,
`HPVload`, and `CD8T` are separate hardcoded literals, `8.0`/`5.0`/`1.0`,
not sourced from these parameters at all), but their names are spelled
identically to the auto-generated init symbol for the identically-named
compartments `SCCAg`, `HPVload`, and `CD8T`. Same incidental-collision
pattern already logged for `sepsis` (#30), `chronic-lymphocytic-leukemia`
(#44), `copd` (#63), `myotonic-dystrophy` (#65), and `gaucher-disease`
(#69).

**Confirmed upstream:** all three defects reproduce from the untouched
original alone via `POST /model_manifest`, no refactor changes involved,
each surfacing only once the prior one is worked around.

**Fixes applied directly to the delivered `cc_mrgsolve_model_refactored.R`**
(not just a scratch copy), per the guide's settled policy for a
non-compiling original, syntax-only and non-numeric:
1. The `$INIT` block is dropped; its 19 assignments become `<CMT>_0 =
   value;` lines in `$MAIN` (identical values, identical compartments,
   under the refactor's renamed `CENT_CIS`/`PERI_CIS`/.../`ADDUCT_CIS`
   identifiers).
2. All nine `if(... ) X = value;` compartment clamps are deleted outright
   (confirmed dead — see Defect 2 reasoning and the verification below).
3. `SCCAg_0` -> `SCCAg0`, `HPVload_0` -> `HPVload0`, `CD8T_0` -> `CD8T0`
   (all usages updated — there were none besides the declaration itself);
   all three remain unused dead parameters, so this is a pure rename with
   no numeric effect.

All three fixes were applied identically to an in-memory-only scratch
copy of the original (never to `cc_mrgsolve_model.R` itself, which still
carries all three defects exactly as written) so it could build for the
`/run_simulation` comparison.

**Defect 4 (non-blocking, math fragility, not a build defect): tumor
volume underflows through zero and NaNs via `log(TV_max/TV)` once
radiotherapy is active.** The original's own `S2` scenario (cisplatin CCRT,
`RT_flag=1`, weekly cisplatin x6, no floor on `TV`) drives `TV` down
through double-precision denormal range (`~3.39e-22` at `t=83`) and then to
`NaN` at `t=84`, propagating into every disease-side output for the rest
of the 730-day run. This is not new: the `if(TV < 0.01) TV = 0.01;` clamp
in Defect 2 was clearly meant to guard exactly this, but (per Defect 2) it
was already a no-op even when it compiled, since only `dxdt_TV` feeds the
integrator — so this scenario has apparently never been runnable to
completion end-to-end, only appearing to work because the model never
built at all before now. **Confirmed identical in both the original
(patched only for Defects 1-3, to reach this point) and the refactored
file:** both NaN at the identical timestep (`t=84`), with identical
finite values up to and including `t=83`, at every compared output
(`CIS_Conc`, `PAC_Conc`, `BEV_Conc`, `PEMBRO_Conc`, `TVADC_Conc`,
`MMAE_lvl`, `VEGF_free`, `PtDNA_rel`, `RT_damage`, `TumorVol`, `SCCAg_lvl`,
`HPV_rel`, `CD8T_rel`, `PDL1_rel`, `TV_change`) — proof this is a
pre-existing fragility in the shared math, not something the refactor
introduced. Not fixed anywhere (guarding it would change the model's own
numeric behavior, out of scope for a rename-only refactor); disclosed in
full in `cc_refactor_notes.md`, which reports the comparison over the
finite prefix (`t=0..83`) for this scenario, exact match, and the
identical NaN onset for the remainder.

Verified: all six of the original's own scenarios (S1 natural history; S2
cisplatin CCRT; S3 CCRT+pembrolizumab; S4 chemo+bevacizumab GOG-240; S5
tisotumab vedotin monotherapy; S6 chemo+bevacizumab+pembrolizumab
KEYNOTE-826), run with identical dosing through both the Defects-1-3-patched
original and `cc_mrgsolve_model_refactored.R`, matched exactly (max abs
diff 0.0) at every one of the 15 shared captured outputs, over the full
731-751-point time grids (S2/S3's grids truncated only by Defect 4's NaN
onset, identical on both sides). See
`cervical-cancer/cc_refactor_notes.md`.

**Fix upstream would be:** drop the separate `$INIT` block and set initial
values via `<CMT>_0` in `$MAIN` instead; delete the nine dead compartment
clamps; rename `SCCAg_0`/`HPVload_0`/`CD8T_0` away from the
`<compartment>_0` pattern (or delete them, since none is read anywhere);
and add a real floor to `TV` inside `dxdt_TV`'s own computation (e.g.
clamping the `log(TV_max/TV)` argument, not reassigning `TV` itself) if
the RT-active scenario is meant to run to completion without NaN-ing.

## 77. `bile-acid-diarrhea/bam_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: two `$ODE` multi-declarator `double` lines drop the type from every name after the first

Found while refactoring this file's seven compounds (colesevelam `SEQ`,
elobixibat `ELO`, loperamide `LOP`, obeticholic acid `OCA`, ondansetron
`OND`, rifaximin `RIF`, tropifexor `TRO`) per `FORK_WORKFLOW_GUIDE.md`
Part 2. Same defect class as #37 (diabetic-ketoacidosis) and the
`pediatric-drug-interactions` finding logged there: mrgsolve 2.0.1's
`$ODE` preprocessing under `$PLUGIN autodec` strips `double` from every
comma-separated declarator but the first:

```
double J12 = kt*ILE1, J23 = kt*ILE2, J3c = kile*TRANS*ILE3;
double fb1 = ILE1/ILEtot, fb2 = ILE2/ILEtot, fb3 = ILE3/ILEtot;
```
and, further down:
```
double fc1 = CCA/tot_col, fc2 = CCDCA/tot_col;
double fc3 = CDCA/tot_col, fc4 = CLCA/tot_col;
```

gives:

```
711:16: error: 'J23' was not declared in this scope
712:20: error: 'fb2' was not declared in this scope
712:39: error: 'fb3' was not declared in this scope
743:20: error: 'fc2' was not declared in this scope
744:21: error: 'fc4' was not declared in this scope
```

Both pairs of lines are core disease-side scaffolding (ileal transit
fluxes across the three terminal-ileum tanks; colonic bile-acid pool
fractions used by the microbial-transformation and secretory-drive
blocks) — neither is inside any of the seven refactored compounds' own
PK/effect blocks, so this is incidental to the refactor's own scope, not
part of it. Confirmed upstream: reproduces from the untouched original
alone via `POST /model_manifest`, no refactor changes involved.

**Fix applied directly to the delivered `bam_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original, syntax-only
and non-numeric): each multi-declarator line was split into separate
`double` statements, same values, same order —

```
double J12 = kt*ILE1;
double J23 = kt*ILE2;
double J3c = kile*TRANS*ILE3;
double fb1 = ILE1/ILEtot;
double fb2 = ILE2/ILEtot;
double fb3 = ILE3/ILEtot;
```
and likewise for `fc1`–`fc4`. The same fix was applied in-memory only to
an untouched scratch copy of the original (never to the tracked
`bam_mrgsolve_model.R`, which still carries the defect exactly as written)
so it could build for the `/run_simulation` comparison. Neither fix
touches a single number, equation, or PK/PD behavior — confirmed by the
exact-match (max abs diff 0.0) verification results across all 17
scenarios in `bile-acid-diarrhea/bam_refactor_notes.md`.

**Fix upstream would be:** split each multi-declarator `double` line into
one statement per declarator (or avoid `$PLUGIN autodec` for these lines
by declaring each variable individually to begin with).

## 78. `narcolepsy/narc_mrgsolve_model.R` — `$CAPTURE` duplicates ten `$CMT` compartment names, fails to build

Found while refactoring this file's five compounds (sodium oxybate `OXY`,
modafinil `MOD`, pitolisant `PIT`, solriamfetol `SOL`, venlafaxine `VEN`)
per `FORK_WORKFLOW_GUIDE.md` Part 2. Same defect class as #30/#34/#41/#45/
#46: the model's `$CAPTURE` block lists ten names that are also `$CMT`
compartments —

```
$CAPTURE
CP_OXY_out CP_MOD_out CP_PIT_out CP_SOL_out CP_VEN_out
CSF_OREXIN OREXIN_NEURONS
EDS_NOW CAT_NOW SLEEP_LATENCY REM_PCT SLEEP_EFF CIRC_OUT
WAKE_STATE REM_STATE NREM_STATE VLPO_ACT SLEEP_P ADENOSINE
WAKE_LC WAKE_TMN WAKE_VTA WAKE_DRN
EFFECT_OXY_EDS EFFECT_MOD_EDS EFFECT_PIT_EDS EFFECT_SOL_EDS EFFECT_VEN_CAT
```

`WAKE_STATE, REM_STATE, NREM_STATE, VLPO_ACT, SLEEP_P, ADENOSINE, WAKE_LC,
WAKE_TMN, WAKE_VTA, WAKE_DRN` are all declared in `$CMT` earlier in the
same file. `POST /model_manifest` on the untouched original (via the
qspserver `mrgsolve_api` container) fails outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  WAKE_STATE,REM_STATE,NREM_STATE,VLPO_ACT,SLEEP_P,ADENOSINE,WAKE_LC,
  WAKE_TMN,WAKE_VTA,WAKE_DRN
```

Confirmed upstream: reproduces from the untouched original alone, unrelated
to any of the five refactored compounds' own PK/effect blocks (all ten
duplicated names belong to the disease/neuroscience side of the model).

**Fix applied directly to the delivered `narc_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original, syntax-only
and non-numeric): the ten compartment names are dropped from the
`$CAPTURE` line — compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported.
The tracked `narc_mrgsolve_model.R` is completely untouched and still
carries the defect exactly as written. Confirmed by the exact-match (max
abs diff 0.0) verification results across all five compounds in
`narcolepsy/narc_refactor_notes.md`.

**Fix upstream would be:** remove the ten compartment names from the
`$CAPTURE` line (they need no explicit listing — mrgsolve reports every
`$CMT` state automatically).


## 79. `ovarian-cancer/oc_mrgsolve_model.R` does not compile under mrgsolve 2.0.1 (three layered build defects), plus three non-blocking modeling defects found while refactoring its five compounds

Found while refactoring this file's five compounds (bevacizumab `BEV`,
carboplatin `CAR`, niraparib `NIRA`, olaparib `OLA`, paclitaxel `PAC`) per
`FORK_WORKFLOW_GUIDE.md` Part 2.

**Build defects (confirmed via `POST /model_manifest` on the untouched
original alone, no refactor changes involved):**

1. **`$CMT`+`$INIT` jointly redeclare all 18 compartments** (`Duplicated
   model names`), the same defect class as #76 (cervical-cancer):
   ```
   Error in validObject(.Object) :
     invalid class "mrgmod" object: Duplicated model names: CAR_C1 CAR_C2
     PAC_C1 PAC_C2 PAC_C3 OLA_gut OLA_C1 OLA_C2 NIRA_C1 NIRA_C2 BEV_C1
     BEV_C2 VEGF TV CA125 Pt_DNA CD8T HRD
   ```
2. **Two `$PARAM` baseline names collide with mrgsolve's own auto-reserved
   `<CMT>_0` symbols**: `CA125_0` (baseline CA-125, 300.0) collides with
   the auto-init symbol for compartment `CA125`; `CD8T_0` (baseline CD8T,
   1.0) collides with the auto-init symbol for compartment `CD8T`. Both
   parameters are otherwise dead — grepped and confirmed unused anywhere
   except a pair of already-broken `$MAIN` lines (`CA125_0_ = CA125_0;`,
   `CD8T_0_ = CD8T_0;`, note the stray trailing underscore on the target,
   which never actually reaches mrgsolve's real `CA125_0`/`CD8T_0`
   init-symbols even before the collision is considered — these two lines
   were always dead code).
3. **A `$ODE`-local `double VEGF_free` collides with `$TABLE`'s `capture
   VEGF_free`** (`redefinition of {anonymous}::VEGF_free`) — the same
   name is used both as a local alias inside `$ODE` (`double VEGF_free =
   VEGF;`, feeding `BEV_effect`) and as the output column name in
   `$TABLE`.
4. **Six `$ODE` lines clamp a compartment state directly** (`if(Pt_DNA<0)
   Pt_DNA=0;`, `if(HRD<0) HRD=0;`, `if(HRD>1) HRD=1;`, `if(CD8T<0)
   CD8T=0;`, `if(TV<0.01) TV=0.01;`, `if(CA125<1) CA125=1;`), which
   mrgsolve 2.0.1 rejects as `assignment of read-only reference` — the
   same defect class as #76's nine clamps, and by the same reasoning
   confirmed a behavioral no-op in any mrgsolve version (only `dxdt_*`
   feeds the integrator).

**Fix applied directly to the delivered `oc_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original, syntax-only
and non-numeric): the `$INIT` block is dropped and its 18 assignments
become `<CMT>_0 = value;` lines in `$MAIN`; `CA125_0`/`CD8T_0` are renamed
`CA125_BASE`/`CD8T_BASE` (dead elsewhere, so purely cosmetic); the
`$ODE`-local `VEGF_free` alias is removed (its one use, `BEV_effect`,
reads `VEGF` directly instead — an identity substitution, not a numeric
change); and all six dead clamps are deleted. The checked-in
`oc_mrgsolve_model.R` is left completely untouched and still carries all
four defects exactly as written. Confirmed by the exact-match
(floating-point-scale max abs diff, largest 0.0024 on a CA125 state of
order ~300) verification results across all six of the original's own
scenarios in `ovarian-cancer/oc_refactor_notes.md`.

**Non-blocking modeling defects (present in the original's own numerical
behavior, not build-related; preserved as-is in the refactor, not fixed,
per the guide's log-don't-fix rule):**

5. **Niraparib has no real absorption depot compartment.** `$PARAM`
   declares `ka_NIRA`/`F_NIRA` (an oral absorption rate and
   bioavailability, mirroring olaparib's real depot-based pattern:
   `dxdt_OLA_C1 = (F_OLA*ka_OLA*OLA_gut)/V1_OLA - ...`), but no
   `NIRA_gut` compartment is ever declared in `$CMT`, and
   `dose_niraparib()` doses directly into `NIRA_C1`. The corresponding
   term in `dxdt_NIRA_C1` reads `(F_NIRA*ka_NIRA*NIRA_C1)/V1_NIRA` —
   i.e. it uses the central compartment's own amount as if it were a
   separate depot's, creating a small stray positive-feedback term on
   `NIRA_C1` itself (`+0.73*0.36/537` is about `+0.00049` per time unit,
   versus a combined `-(CL_NIRA+Q_NIRA)/V1_NIRA` of about `-0.0392` —
   roughly a 1.2% relative offset, not catastrophic, but a genuine
   authoring defect, most likely a copy-paste from olaparib's depot
   pattern that dropped the depot state). Preserved exactly (renamed
   identifiers only) in the refactored file; confirmed via the exact-
   match verification that this reproduces identically on both sides.
6. **Paclitaxel's own header comment claims "Michaelis-Menten nonlinear"
   PK (citing Gianni 1995 JCO), but the actual `dxdt_PAC_C1` equation is
   pure linear elimination** — a constant `CL_PAC` divided by `V1_PAC`,
   no `Km`/`Vmax` saturation term anywhere in the block. The three-
   compartment linear structure is real (and was refactored faithfully,
   Archetype 2 extended with a second peripheral compartment), but the
   claimed nonlinear kinetics are not implemented.
7. **The olaparib/niraparib maintenance dosing schedule uses an
   hour-scaled numeric axis while the rest of the model (chemo cycles,
   Gompertz growth, CA-125 turnover, the R script's own `mrgsim(...,
   end=730, delta=1)` call, labeled "2-year simulation (days)") uses a
   day-scaled one, on the same raw `time` axis, with no unit
   reconciliation.** `dose_olaparib()`/`dose_niraparib()` compute event
   times as `start_d*24` / `(start_d+dur_d)*24`, e.g. scenario S4's
   `dose_olaparib(start_d=126, dur_d=604)` schedules its first dose at
   `time=3024` and its last at `time=17508` — on an axis where
   carboplatin is dosed at `time=21,42,...` (meant as days) and the
   model's own kill/growth kinetics (`kg=0.008`, `kdeg_CA125=0.03`, all
   "1/day") are calibrated for a ~730-unit run. Confirmed via `POST
   /run_simulation` against the qspserver `mrgsolve_api`: with S4's own
   dosing and `time: {end: 730, delta: 1}` requested, the returned time
   grid is unioned with the event times and actually runs to `t=17508`
   (about 24x the requested/intended horizon) — meaning no olaparib or
   niraparib dose ever lands within the model's own intended 0-730-day
   analysis window at all (first dose in S4 lands at t=3024, in S5 at
   t=3024), so the file's own `summarize_scenario()` (PFS/CA-125-nadir
   over what the header calls a 2-year run) would never actually see any
   PARPi maintenance effect within its intended window in either S4, S5,
   or S6, despite the header's own claim that these three scenarios are
   "calibrated to SOLO-1, PRIMA, PAOLA-1." Reproduces identically (exact
   match) on both the original and the refactored file, since neither
   the DSL nor the R-side `dose_*()` argument values were changed — only
   `cmt=` identifiers were renamed to match the refactor's renamed
   compartments.
8. **Bevacizumab has no mechanistic effect on tumor volume, CA-125, or
   CD8+ T cells anywhere in this file** — `BEV_effect` (renamed
   `VEGF_bind_BEV`) is a pure mass-action VEGF sink
   (`kbind_BEV*C_BEV*VEGF`) feeding only `dxdt_VEGF`; no downstream
   equation (`dxdt_TV`, `dxdt_CA125`, `dxdt_CD8T`) ever reads `VEGF` or
   any bevacizumab-derived quantity. Unlike `cervical-cancer/cc_mrgsolve_
   model.R` (#76), which has a real, if separate, anti-angiogenic
   `BEV_C1/(BEV_C1+10.0)` Hill-shaped kill term, this file's bevacizumab
   arms (S3, S6) carry chemotherapy's own tumor-kill effect only —
   bevacizumab contributes nothing to tumor-volume, CA-125, or T-cell
   trajectories despite being a named treatment arm and despite the
   header's calibration claim against PAOLA-1 (an anti-angiogenic-add-on
   trial). Confirmed by reading every `dxdt_TV`/`dxdt_CA125`/`dxdt_CD8T`
   line in the file; not fixed (inventing a new effect the original
   never had is out of scope for a rename-only refactor) — see "No
   `EFFECT_BEV`" in `ovarian-cancer/oc_refactor_notes.md`.

**Fix upstream would be:** drop the separate `$INIT` block and set
initial values via `<CMT>_0` in `$MAIN` instead; rename `CA125_0`/
`CD8T_0` away from the `<compartment>_0` pattern; rename the `$ODE`-local
`VEGF_free` alias; delete the six dead compartment clamps; add a real
`NIRA_gut` depot compartment for niraparib (or drop `ka_NIRA`/`F_NIRA` if
direct central dosing is intentional); either remove the "Michaelis-
Menten" claim from paclitaxel's header comment or actually implement
saturable clearance; rescale `dose_olaparib()`/`dose_niraparib()`'s event
times to the same day-scale axis as the rest of the model (drop the
`*24` factors, or switch the whole model to an hour-scale axis
consistently); and add a real anti-angiogenic tumor-kill term driven by
`VEGF`/`BEV` if bevacizumab is meant to affect tumor growth in this
model, as its own scenario labels and PAOLA-1 calibration claim it does.

## 80. `stable-angina/sa_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$TABLE` re-declares `Cp_BB`/`Cp_CCB`/`Cp_RAN`/`Cp_IVA`/`Cp_NIT` a second time via self-referential `capture NAME = NAME;`

Found while refactoring BB (bisoprolol)/CCB (amlodipine)/RAN
(ranolazine)/IVA (ivabradine)/NIT (ISMN nitrate) PK/PD per this file's
`driver-patches/data/compound_perturbation_census.md` rows. Reproduces
identically from the untouched original via the qspserver `mrgsolve_api`
container's `/model_manifest`, no refactor changes involved, unrelated to
any compound's own PK/PD shape.

**Defect (build-blocking): `$ODE` declares `double Cp_BB = CENTRAL_BB /
V1_BB * 1000.0;` (and identically for `Cp_CCB`, `Cp_RAN`, `Cp_IVA`,
`Cp_NIT`), then `$TABLE` immediately writes `capture Cp_BB = Cp_BB;`.**
Same underlying mrgsolve mechanism already logged for other files in this
run (#66, #68, #78): the `capture NAME = expr;` shorthand both declares
`NAME` as a `double` and registers it for output, so when `NAME` is
already a declared `double` in scope — here, because the report variable
happens to share its own capture name — the shorthand's redeclaration
collides with the one already there:

```
68:11: error: redefinition of 'capture {anonymous}::Cp_BB'
   68 |   capture Cp_BB;
      |           ^~~~~
30:10: note: 'double {anonymous}::Cp_BB' previously declared here
   30 |   double Cp_BB;
      |          ^~~~~
```

and identically for `Cp_CCB`, `Cp_RAN`, `Cp_IVA`, `Cp_NIT` (five
occurrences total, one per compound's own plasma concentration).

**Confirmed upstream:** reproduces from the untouched original alone via
`POST /model_manifest`, no refactor changes involved.

**Fix applied directly to the delivered
`sa_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original
(`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at all"),
syntax-only and non-numeric: the renamed concentration terms (`C_BB`,
`C_CCB`, `C_RAN`, `C_IVA`, `C_NIT`) and every `EFFECT_<STEM>` term are
predeclared as `double`s in `$GLOBAL`, assigned by plain statement
(without a re-declaring `double`) in `$ODE`, and listed bare (no
`capture NAME = expr;` shorthand) in a `$CAPTURE` block — the same
pattern already established in the AMD/membranous-nephropathy/sepsis/
SAH/AIP/breast-cancer refactors, which avoids the self-referential
collision entirely rather than reproducing it under the new names. A
syntax-fixed scratch copy of the original (renaming only the five
colliding `$ODE`-local doubles to `Cp_<STEM>_val`, leaving every
capture's own output name unchanged) was built in-memory for the
`/run_simulation` comparison; the checked-in `sa_mrgsolve_model.R` itself
was left completely untouched and still carries the defect exactly as
written. Verified: with this fix applied to each side, all 37 shared
outputs matched exactly (max abs diff = 0.0) across all six of the
original's own named dosing scenarios, its own bisoprolol dose-response
sweep (5 doses), and its own 2-year statin/no-statin plaque simulation —
no solver step-count issues at any of these horizons (up to 731 points).
See `stable-angina/sa_refactor_notes.md`.
## 81. `urolithiasis/uri_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$OMEGA`/`$SIGMA` `@block 0` misparses as a flat diagonal, `$CMT`+`$INIT` jointly redeclare all 17 compartments, and two `$ODE` lines write directly to read-only compartment state

Found while refactoring HCTZ (thiazide), ALLO/OXP (allopurinol/its own
active metabolite oxypurinol), KCIT (potassium citrate), and TAM
(tamsulosin) PK/PD per this file's
`driver-patches/data/compound_perturbation_census.md` rows. All three
defects reproduce identically from the untouched original via the
qspserver `mrgsolve_api` container's `/model_manifest`, no refactor
changes involved, unrelated to any of the four compounds' own PK/PD shape.

**1. `$OMEGA @block 0` / `$SIGMA @block 0` (build-blocking, found first).**
The original writes:
```
$OMEGA @block 0
@labels ETA_CA ETA_OX ETA_UA ETA_Stone
0.04
0.02 0.09
0.02 0.03 0.04
0.00 0.00 0.00 0.00

$SIGMA @block 0
@labels EPS_uCa EPS_uOx EPS_uUA
0.04 0.02 0.04
```
mrgsolve 2.0.1 does not accept a trailing numeric argument on `@block`;
rather than rejecting it outright, it silently falls back to reading
every number across the block as one flat sequence of independent
diagonal entries (10 numbers for `$OMEGA`, ignoring the intended
lower-triangular 4x4 grouping implied by the 1/2/3/4-value row shape),
which then fails validation against the 4 (`$OMEGA`) or 3 (`$SIGMA`)
declared `@labels`:
```
Error in validObject(.Object) :
  invalid class "omegalist" object: Length of labels (4) does not match the matrix rows (10).
```
**2. `$CMT` (`@annotated`) and `$INIT` jointly redeclare all 17
compartments** (`HCTZ1`/`HCTZ2`/`HCTZ3`/`ALLO1`/`ALLO2`/`OXP`/`KCIT1`/
`KCIT2`/`TAM1`/`TAM2`/`CA_PL`/`OX_PL`/`UA_PL`/`CIT_PL`/`STONE`/`INFL`/
`GFR_CUR`) — same defect class as issues #34/#36/etc: mrgsolve 2.0.1
treats `$INIT`'s `NAME = value` list as its own compartment declaration
rather than a companion to `$CMT`:
```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: HCTZ1 HCTZ2 HCTZ3 ALLO1 ALLO2 OXP KCIT1 KCIT2 TAM1 TAM2 CA_PL OX_PL UA_PL CIT_PL STONE INFL GFR_CUR
```
**3. `$ODE` writes directly to two compartment states to clamp them** —
`if(STONE < 0) STONE = 0;` and `if(GFR_CUR < 15) GFR_CUR = 15;`. Same
defect class as issue #36 (age-related-macular-degeneration): mrgsolve
2.0.1 passes every `$CMT` state into the generated `$ODE` function as a
`const double&`, so writing to it is a hard compile error:
```
492:21: error: assignment of read-only reference 'STONE'
495:26: error: assignment of read-only reference 'GFR_CUR'
```
As with #36, this is a no-op even where compilable — only `dxdt_*` feeds
mrgsolve's integrator, so reassigning the state variable itself inside
`$ODE` can never affect the next solver step or the reported trajectory,
compilable or not.

**Confirmed upstream:** all three reproduce from the untouched original
alone via `POST /model_manifest`; the model does not build at all, for
either the original or the refactored file, until worked around.

**Fix applied directly to the delivered
`uri_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original, all syntax-only and
non-numeric: (a) `$OMEGA @block 0` → `$OMEGA @block` (the 10 values now
parse as the intended 4x4 lower-triangular block, matching the 4
`@labels`, no value changed) and `$SIGMA @block 0` → `$SIGMA` (the 3
values now parse as the intended 3x3 diagonal, matching the 3 `@labels`,
no value changed); (b) the `$INIT` block deleted, its 17 assignments
moved into `$MAIN` using the modern `<CMT>_0 = value;` idiom under the
refactor's renamed compartments (e.g. `GUT_HCTZ_0 = 0;`), same values,
no compartment added or removed; (c) the two `if(...) ... = ...;` clamp
lines deleted outright (see the no-op reasoning above). An identically
syntax-fixed scratch copy of the untouched original (same three fixes,
original compartment/param names) was built in-memory, never committed,
purely as the verification comparison target; the checked-in
`uri_mrgsolve_model.R` itself was left completely untouched and still
carries all three defects exactly as written.

Verified: with this fix applied to each side, all 17 compartments and 19
deterministic `$CAPTURE` outputs matched exactly (max abs diff = 0.0)
across all six of the original's own named dosing scenarios (untreated,
HCTZ 25mg, K-citrate 60mEq, allopurinol 300mg, lifestyle+combo,
metabolic-syndrome combo), each run for the original's own full 5-year
(500-point) horizon — no solver step-count issues at this horizon. Three
`$TABLE` captures (`uCa_mgday`/`uOx_mgday`/`uUA_mgday`) carry the
original's own `$SIGMA`-drawn residual error (`* (1.0 + EPS(i))`) and are
therefore stochastic from one API call to the next by design (confirmed:
re-running the identical original twice gives different `uCa_mgday`
values but identical `StoneSize_mm`) — excluded from the exact-match
check for that reason alone, not because of any refactor-introduced
difference; every deterministic output, including the pre-noise
`uCa_24`-derived risk flags, matched exactly. See
`urolithiasis/uri_refactor_notes.md`.

**Fix upstream would be:** drop the stray `0` after `@block` in both
`$OMEGA` and `$SIGMA`; remove the `$INIT` block and set initial
conditions via `<CMT>_0 = value;` in `$MAIN` instead; delete the
`if(STONE < 0) STONE = 0;` / `if(GFR_CUR < 15) GFR_CUR = 15;` lines (or,
if the negative/floor guard is actually wanted, implement it by zeroing
the offending `dxdt_*` term instead of assigning the state directly).

## 82. `reactive-arthritis/rea_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 26 compartments, and `$ODE` reads a bare `chlamydia` identifier never declared anywhere in the file or its R wrapper

Found while refactoring NSAID, SSZ (Sulfasalazine), MTX (Methotrexate),
TNFi (renamed TNFI, a TNF inhibitor), and IL17i (renamed IL17I, an IL-17
inhibitor) PK/PD per this file's
`driver-patches/data/compound_perturbation_census.md` rows. Both defects
reproduce identically from the untouched original via the qspserver
`mrgsolve_api` container's `/model_manifest`, no refactor changes
involved, unrelated to any of the five compounds' own PK/PD shape.

**1. `$CMT` (`@annotated`) and `$INIT` jointly redeclare all 26
compartments** (`PATH CHRON_PATH INNATE NEUT MACRO TH1 TH17 TREG TNF_c
IL17_c IL6_c IL10_c IFNg_c SYNOV CARTDMG PAIN CRP_c JCOUNT NSAID_DEPOT
NSAID_C SSZ_DEPOT SSZ_C MTX_DEPOT MTX_C TNFi_DEPOT TNFi_C IL17i_DEPOT
IL17i_C`) — same defect class as issues #34/#36/#42/#44/#48/#49/#51/#61/
#62/#63/#64/#76/#81: mrgsolve 2.0.1 treats `$INIT`'s `NAME = value` list
as its own compartment declaration rather than a companion to `$CMT`:
```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: PATH CHRON_PATH
  INNATE NEUT MACRO TH1 TH17 TREG TNF_c IL17_c IL6_c IL10_c IFNg_c SYNOV
  CARTDMG PAIN CRP_c JCOUNT NSAID_DEPOT NSAID_C SSZ_DEPOT SSZ_C MTX_DEPOT
  MTX_C TNFi_DEPOT TNFi_C IL17i_DEPOT IL17i_C
```
**2. `$ODE`'s `dxdt_CHRON_PATH` reads a bare identifier `chlamydia` that
is never declared anywhere** — not in `$PARAM`, not as a `$CMT`
compartment, and not even set via `param()`/`mrgsim(param = ...)`
anywhere in the surrounding R wrapper's five scenarios. A missing
declaration, not a name collision — a new defect class for this file
(same class as issue #64's `DOSE_*` case): once fix #1 above is applied
in isolation, the model still fails, now on a hard C++ compile error:
```
415:27: error: 'chlamydia' was not declared in this scope
dxdt_CHRON_PATH = 0.002 * chlamydia - k_chron_decay * CHRON_PATH - k_abx_chlamyd * CHRON_PATH;
                          ^~~~~~~~~
```
The surrounding `$PROB`/comment text ("Chlamydia: small fraction seeds
persistent synovial compartment") suggests `PATH` (the pathogen-load
compartment) may have been the intended read, but nothing in the file or
its R wrapper ever assigns or references a `chlamydia` symbol of any
kind, so there is no way to recover original intent from the file alone.

**Confirmed upstream:** both reproduce from the untouched original alone
via `POST /model_manifest`; the model does not build at all, for either
the original or the refactored file, until both are worked around.

**Fix applied directly to the delivered
`rea_mrgsolve_model_refactored.R`** (not just a scratch copy), per the
guide's settled policy for a non-compiling original, both syntax-only and
non-numeric: (a) the `$INIT` block deleted, its 26 assignments moved into
`$MAIN` using the modern `<CMT>_0 = value;` idiom under the refactor's
renamed PK compartments (e.g. `GUT_NSAID_0 = 0;`), same values, no
compartment added or removed; (b) `chlamydia : 0` added to `$PARAM` —
matching the value the original's own R wrapper already implicitly used
everywhere (it is never set to anything else in any of the file's 5
scenarios), so this is not an invented value, just the original's own
default-by-omission made explicit and declarable. An identically
syntax-fixed scratch copy of the untouched original (same two fixes,
original compartment/param names) was built in-memory, never committed,
purely as the verification comparison target; the checked-in
`rea_mrgsolve_model.R` itself was left completely untouched and still
carries both defects exactly as written.

Verified: with this fix applied to each side, via `POST /run_simulation`
on the qspserver `mrgsolve_api` container, every shared output matched
exactly (max abs diff = 0.0): (1) the original's own Scenario 5 regimen
(SSZ 1000mg q12h continuous, NSAID 500mg q12h for the first 84 days, ETN
50mg weekly starting day 84), shortened to a 100-day window per the
guide's solver-budget allowance, 772 points across 14 outputs spanning
NSAID/SSZ/TNFi PK and every downstream disease endpoint; (2) MTX dosed
via the original's own (defined but never actually invoked in any of its
5 shipped scenarios) `make_mtx_events()` helper — 15mg/wk SC, shortened to
56 days/9 doses, 234 points across 8 outputs; (3) IL17i dosed via the
original's own, likewise-unused, `make_il17i_events()` helper — 300mg
loading + monthly, shortened to 90 days/7 doses, 368 points across 5
outputs. One `IL17i_obs` comparison run initially came back with the
refactored side's PK/effect outputs pinned at 0 and several downstream
disease outputs as `NaN`, alongside a `VAS_pain` value that looked like a
different simulation entirely — re-running the identical request pair a
few seconds apart reproduced a clean exact match, confirming this was the
container's own known intermittent instability under load (see e.g. the
bronchiectasis and essential-thrombocythemia entries' `mrgsim_api`
flakiness notes), not a property of either model. See
`reactive-arthritis/rea_refactor_notes.md`.

**Also found, not a build defect:** the original's own R wrapper defines
`make_mtx_events()` and `make_il17i_events()` helper functions (with
sensible clinical defaults matching the `$PROB`-commented doses, 15mg/wk
SC MTX and 300mg-load-then-monthly SC secukinumab) but never actually
calls either one in any of its 5 named scenarios — MTX and IL17i's PK/PD
math exists and is fully wired into the disease equations, but the
shipped script never exercises it. Not fixed (adding a new scenario to
the original is out of scope for a rename-only refactor) — flagged here
since it meant verification for these two compounds had to dose via the
original's own (unused) helper functions rather than one of its named
scenarios; see the "None of the file's own scenarios dose MTX or IL17i"
note in `reactive-arthritis/rea_refactor_notes.md`.

**Fix upstream would be:** drop the separate `$INIT` block and set
initial values via `<CMT>_0` in `$MAIN` instead; either declare
`chlamydia` in `$PARAM` (defaulting to, e.g., 0 or 1 depending on
intended semantics) or replace the bare `chlamydia` read in
`dxdt_CHRON_PATH` with the `PATH` compartment the surrounding comment
appears to intend; and either wire `make_mtx_events()`/
`make_il17i_events()` into a sixth/seventh scenario or remove them if
they are meant to stay illustrative-only.

## 83. `venous-thromboembolism/vte_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$PARAM`'s `VK0_red` computes its default from a sibling `$PARAM` name, `$CMT`+`$INIT` jointly redeclare all 18 compartments, and `$TABLE` re-declares six `$MAIN` doubles a second time via self-referential `capture NAME = NAME;`

Found while refactoring APIX (apixaban), DABI (dabigatran), ENOX
(enoxaparin), RIV (rivaroxaban), and WARF (warfarin) PK/PD per this
file's `driver-patches/data/compound_perturbation_census.md` rows. Three
independent defects, all unrelated to any compound's own PK/PD shape,
reproduce identically from the untouched original via the qspserver
`mrgsolve_api` container's `/model_manifest`:

1. `$PARAM VK0_red = VK0_ox * KIN_VK / KOUT_VK` — mrgsolve evaluates a
   `$PARAM` block's defaults as one R `list(...)` call with no access to
   sibling elements during that evaluation, so this fails with
   `object 'VK0_ox' not found`.
2. `$CMT` declares all 18 compartments (7 PK + 11 PD/factor-pool) and a
   separate `$INIT` block assigns every one of them an initial value —
   this build rejects the combination outright: `invalid class "mrgmod"
   object: Duplicated model names: RIV_GUT RIV_CENT ... FII_POOL`.
3. `$MAIN` computes `Cp_RIV`/`Cp_APIX`/`Cp_DABI`/`Cp_WARF`/`ANTI_XA`/`INR`
   as doubles, then `$TABLE` re-declares each one under the identical
   name (`capture Cp_RIV = Cp_RIV;`, etc.) — the same defect *class* as
   entry #43 (DIC's `FREETPA`), but at `$MAIN` scope rather than `$ODE`:
   `error: redefinition of 'capture {anonymous}::Cp_RIV'`.

See `venous-thromboembolism/vte_refactor_notes.md` for the syntax-only,
non-numeric fixes applied to the refactored sibling and the verification
proving all three are behaviorally inert (max relative deviation 0.0
across 16 shared outputs, in all 7 of the file's own dosing scenarios).

**Fix upstream would be:** replace `VK0_red`'s computed default with a
literal value; drop the separate `$INIT` block and set every initial
value via `<CMT>_0` in `$MAIN` instead; and delete the six redundant
self-referential `$TABLE capture` lines (the bare `$MAIN` doubles are
already reportable via `$CAPTURE` directly, with no `$TABLE` line needed).

## 84. `venous-thromboembolism/vte_mrgsolve_model.R`'s warfarin Hill term reuses rivaroxaban's own `EMAX_RIV` as warfarin's effect ceiling — no `EMAX_WARF` parameter exists anywhere in the file

Not a build defect (the file fails to compile for the unrelated reasons
in entry #83) but a genuine cross-compound parameter leak, found while
isolating warfarin's PK/PD block for the same refactor:

```c
double WARF_INH = (EMAX_RIV * pow(Cp_WARF, HILL_W)) /
                  (pow(IC50_WARF, HILL_W) + pow(Cp_WARF, HILL_W));
```

`EMAX_RIV` (rivaroxaban's own maximum-FXa-inhibition parameter, default
0.97) is read directly as warfarin's Hill ceiling; grepping the whole
file confirms no `EMAX_WARF` parameter is ever declared. This means
overriding `EMAX_RIV` for a rivaroxaban-only sensitivity analysis would
silently also change warfarin's factor-synthesis-inhibition curve, and it
buries warfarin's effect inside an expression that is not, in fact,
warfarin's own. Confirmed upstream by reading every `$PARAM` declaration
in the file (`grep -n "EMAX_WARF"` — no match) against every other
compound's `EMAX_<STEM>`, which are all present. Not fixed in the
original, per the never-edit-upstream rule; the refactored sibling adds a
new `EMAX_WARF = 0.97` parameter carrying the identical value so
`EFFECT_WARF` is numerically unchanged (confirmed in
`vte_refactor_notes.md`'s verification — no scenario in the file
overrides `EMAX_RIV`'s default) while being properly self-contained.

**Fix upstream would be:** declare a dedicated `EMAX_WARF` parameter
(clinically, warfarin's own ceiling on VKORC1 inhibition need not equal
rivaroxaban's ceiling on FXa inhibition — these are pharmacologically
unrelated targets that happen to share a numeric default in this file)
and reference it in `WARF_INH` in place of `EMAX_RIV`.

## 85. `pyoderma-gangrenosum/pg_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: twelve `$PARAM @annotated` lines are missing the required description field

Confirmed reproducing identically from the untouched original via the
qspserver `mrgsolve_api` container's `/model_manifest`:

```
Error: improper annotation format
 input: kout_Th17    : 0.05
 context: parse annotated parameter block (PARAM)
```

An `@annotated` `$PARAM` line requires three colon-delimited fields
(`name : value : description`); twelve lines in this file supply only
the first two, so mrgsolve's annotation parser rejects the whole block
on the first one it reaches, then (after that one is patched, for
diagnosis only) on the next, and so on. All twelve are ordinary disease-
side turnover-rate parameters, none inside any of the six refactored
compounds' own PK/PD blocks:

```
kout_Th17    : 0.05
kin_Treg     : 0.06
kout_Treg    : 0.06
kin_M1       : 0.20
kout_M1      : 0.20
kout_ROS     : 0.5
kout_CRP     : 0.3
k_Calp_syn   : 0.3
kout_Calp    : 0.3
IC50_IFX_TNF  : 0.10
IC50_UST_IL23 : 0.50
Emax_drug     : 0.92
```

(The last three sit inside the Infliximab/Ustekinumab/shared-Emax
parameters that the refactor renames anyway, so they were fixed as part
of the rename itself, each gaining a real description in the same edit
that renamed it to `EC50_IFX`/`EC50_UST`/`EMAX_<STEM>`; the first nine
are untouched disease-side parameters and got the settled build-compat
treatment below.)

**Verification workaround, applied to the delivered
`pg_mrgsolve_model_refactored.R` per this fork's settled build-compat
policy (`FORK_WORKFLOW_GUIDE.md`, "When the original doesn't compile at
all")**: a short, purely descriptive third field was appended to each of
the nine untouched lines (e.g. `kout_Th17    : 0.05    : Th17 clearance
(1/d)`), changing nothing about the parameter's name, value, or use in
any equation. Confirmed non-numeric by the exact-match verification in
`pg_refactor_notes.md` (max abs/rel deviation 0.0 across all 8 of the
file's own dosing scenarios). The untouched original
`pg_mrgsolve_model.R` still carries the defect forward unfixed, per the
never-edit-upstream rule.

**Fix upstream would be:** add a real description field to all twelve
`$PARAM` lines above.

## 86. `pyoderma-gangrenosum/pg_mrgsolve_model.R` reuses one shared `Emax_drug` parameter as the effect ceiling for all six of its compounds (Adalimumab, Anakinra, Cyclosporine, Infliximab, Prednisone, Ustekinumab)

Same defect class as entry #84 (VTE's `EMAX_RIV`/`EMAX_WARF` leak), found
while isolating each compound's own Hill effect term for the PK/PD
refactor:

```c
double EFF_ADA_TNF  = Emax_drug * CONC_ADA  / (IC50_ADA_TNF + CONC_ADA);
double EFF_IFX_TNF  = Emax_drug * CONC_IFX  / (IC50_IFX_TNF + CONC_IFX);
double EFF_ANA_IL1  = Emax_drug * CONC_ANA  / (IC50_ANA_IL1 + CONC_ANA);
double EFF_UST_IL23 = Emax_drug * CONC_UST  / (IC50_UST_IL23 + CONC_UST);
double EFF_CSA_Th17 = Emax_drug * CONC_CSA  / (IC50_CSA_Th17 + CONC_CSA);
double EFF_PRED     = Emax_drug * CONC_PRED / (IC50_PRED + CONC_PRED);
```

A single `Emax_drug = 0.92` parameter is read as the maximum-effect
ceiling for six pharmacologically unrelated targets (TNF-alpha blockade
by two different biologics, IL-1 axis blockade, IL-23 blockade,
calcineurin-mediated Th17 suppression, and glucocorticoid-receptor-
mediated NF-kB suppression). Grepping the whole file confirms no
compound-specific `Emax_<drug>` parameter is ever declared. This means a
sensitivity analysis or virtual-population override intended to change
only one compound's ceiling (e.g. a lower observed Adalimumab response
ceiling in a real cohort) would silently also change the other five
compounds' ceilings. Confirmed upstream by reading the full `$PARAM`
block (`grep -n "Emax"` — only the one shared declaration). Not fixed in
the original, per the never-edit-upstream rule; the refactored sibling
declares six independent `EMAX_ADA`/`EMAX_IFX`/`EMAX_ANA`/`EMAX_UST`/
`EMAX_CSA`/`EMAX_PRED` parameters, each carrying the identical `0.92`
default, so every `EFFECT_<STEM>` is numerically unchanged (confirmed in
`pg_refactor_notes.md`'s verification — no scenario in the file
overrides `Emax_drug`'s default) while each compound becomes
independently driveable, per the fork convention.

**Fix upstream would be:** declare six dedicated `Emax_<drug>`
parameters (clinically, each compound's own maximum achievable
suppression of its own target need not be identical to the others') and
reference each in its own compound's effect term in place of the shared
`Emax_drug`.

## 91. `vexas-syndrome/vexas_mrgsolve_model.R` declares an SC/oral bioavailability parameter for every one of its six PK compounds, but none of the six is ever referenced in any `dxdt` equation

Found while refactoring ANA (anakinra), AZA (azacitidine), CAN
(canakinumab), PRED (prednisone), RUX (ruxolitinib), and TOC
(tocilizumab) PK/PD per this file's
`driver-patches/data/compound_perturbation_census.md` rows. `$PARAM`
declares `F_PRED = 0.85`, `F_TOC = 0.80`, `F_ANA = 0.95`, `F_CAN = 0.66`,
`F_RUX = 0.95`, and `F_AZA = 0.89`, each commented as that drug's SC or
oral bioavailability. Grepping every occurrence of each name in the file
shows it appears exactly once — its own `$PARAM` declaration — and never
again in `$MAIN`/`$ODE`:

```
dxdt_AGUT_PRED = -ka_PRED*AGUT_PRED;
dxdt_CC_PRED   =  ka_PRED*AGUT_PRED - CL_PRED/V_PRED*CC_PRED;
```

is the actual depot-to-central transfer for every one of the six
compounds — no `F_<drug>` multiplier anywhere. Confirmed upstream via
`POST /model_manifest` on the untouched original through the qspserver
`mrgsolve_api` container (all six `F_*` parameters listed, each with its
documented default) plus a whole-file grep confirming a single
occurrence per name. The practical effect: every dose in every one of
the file's own `scenarios()` (prednisone taper, tocilizumab, anakinra,
ruxolitinib, azacitidine) reaches the central compartment as if `F=1`,
regardless of the value in the `F_<drug>` comment — e.g. tocilizumab's
documented 80% SC bioavailability is not actually applied, so simulated
exposure runs higher than the model's own stated PK. Not fixed in the
original, per the never-edit-upstream rule; preserved bug-for-bug
(unused parameters carried forward, no `F_<STEM>` multiplier added) in
`vexas_mrgsolve_model_refactored.R`, disclosed in
`vexas_refactor_notes.md`, and confirmed numerically identical to the
original in verification (max relative deviation ~6.5e-13, floating-point
noise near zero) via the qspserver `mrgsolve_api` container.

**Fix upstream would be:** multiply each depot-to-central transfer term
by its corresponding `F_<drug>` (e.g.
`dxdt_CC_PRED = ka_PRED*F_PRED*AGUT_PRED - CL_PRED/V_PRED*CC_PRED;`), which
would reduce simulated exposure for every compound to the documented
bioavailability.

## 92. `vexas-syndrome/vexas_mrgsolve_model.R`'s azacitidine PK is fully disconnected from its own PD: dosing `ADEP_AZA` changes `C_AZA` but nothing in the file reads `C_AZA` to affect the disease state

Found while refactoring the same six compounds as entry #91. The file
declares `k_Aza : 0.0 : Azacitidine clone-modulating rate (1/h)` in
`$PARAM` and uses it in `dxdt_VAF = (k_clone*VAF*(1-VAF)) - k_HSCT*VAF -
k_Aza*VAF;` — but `k_Aza` is a free-standing fixed parameter (default
0), never assigned from `C_AZA` (`= CC_AZA / V_AZA`) anywhere in
`$MAIN`/`$ODE`. Grepping the whole file for `C_AZA` shows it is computed
and captured but read nowhere else; grepping for `k_Aza` shows it is
read once (in `dxdt_VAF`) and never written except by its own `$PARAM`
default. Confirmed upstream via `POST /model_manifest` on the untouched
original (both `k_Aza` and `C_AZA`/`CC_AZA`/`ADEP_AZA` present as
separate, unconnected entries) through the qspserver `mrgsolve_api`
container, plus a run of the file's own `"6_azacitidine"` dosing
scenario showing `C_AZA` rising exactly as its PK dictates while `VAF`'s
trajectory is completely unaffected (`k_Aza` stays at its default 0
throughout). The practical effect: azacitidine dosing, as coded, has
**zero** effect on the modeled VEXAS clone burden or any downstream
biomarker unless a user manually overrides `k_Aza` from outside the
shown DSL — the drug's own PK profile is cosmetic. Not fixed in the
original, per the never-edit-upstream rule; no `EFFECT_AZA`/`EMAX_AZA`/
`EC50_AZA` was invented in `vexas_mrgsolve_model_refactored.R` to paper
over this gap — `C_AZA` is still exposed as the redirect point per the
naming convention, `k_Aza` is preserved exactly as-is, and the
disconnection is disclosed plainly in `vexas_refactor_notes.md`.

**Fix upstream would be:** replace the free-standing `k_Aza` read in
`dxdt_VAF` with a Hill function of `C_AZA` (e.g. `k_Aza_eff =
EMAX_AZA*C_AZA/(EC50_AZA+C_AZA)`), so that azacitidine's own dosing
schedule actually drives clonal suppression instead of a parameter that
only a user script can set.

## 87. `huntingtons-disease/hd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 20 compartments, and `$ODE` calls an R-only function that is undefined at compile time

Confirmed via the qspserver `mrgsolve_api` container (`http://localhost:8007`)
that the **untouched original does not compile at all**, two independent
defects, neither specific to any one of the six refactored compounds
(TBZ/HTBZ/DTBZ/VBZ/RILUZOLE/TOMINERSEN):

1. `$CMT`+`$INIT` jointly redeclare all 20 compartments (same class as
   issues #29/#34/#36/#59/#76/#81/#82/#83 — `$INIT`'s `NAME = value` list
   is its own compartment declaration under mrgsolve 2.0.1, not a
   companion to `$CMT`): `Duplicated model names: TBZ_gut TBZ_plasma
   HTBZ_brain DTBZ_plasma DTBZ_brain VBZ_plasma VBZ_brain tominersen_CSF
   riluzole_plasma riluzole_brain mHTT_mRNA mHTT_prot mHTT_oligo BDNF_cmt
   dopamine_cmt MSN_surv oxidative_idx neuroinflam_idx UHDRS_TMS TFC_cmt`.
2. A distinct, unrelated defect: the `$ODE` block's UHDRS-TMS equation
   calls `mHTT_protein_reduction_f()` —
   ```
   dxdt_UHDRS_TMS = chorea_drive - chorea_tx - riluzole_motor
                  - (kprog_TMS * 0.3) * (mHTT_protein_reduction_f());
   ```
   — but `mHTT_protein_reduction_f` is only ever defined as a plain R
   function **outside** the quoted DSL string (`mHTT_protein_reduction_f
   <- function() 0.0`, itself commented "R-level placeholder — mHTT
   reduction for DMT effect"), several hundred lines below `hd_model <-
   '...'` in the same `.R` file. Since `model_content` is compiled
   server-side as bare mrgsolve DSL text with no surrounding R script
   (see the fork guide's qspserver-compatibility requirement #5), this
   identifier is simply undefined in the generated C++ — a hard compile
   error independent of defect (1). New defect class, not seen in any
   prior entry: an R-only helper called from inside the DSL block rather
   than a param/constant interpolated via `sprintf`/`glue` before
   quoting.

**Fix applied directly to the delivered `hd_mrgsolve_model_refactored.R`**
(not just a scratch copy), both syntax-only and non-numeric: (a) the
`$INIT` block deleted, its 20 assignments moved into `$MAIN` as
`<CMT>_0 = value;` under the refactor's renamed compartments (e.g.
`GUT_TBZ_0 = 0;`, `CENT_HTBZ_0 = 0;`), same values, no compartment added
or removed; (b) the `mHTT_protein_reduction_f()` call deleted outright —
confirmed a no-op since the R function it referenced always returned the
constant `0.0` regardless of state, so removing `- (kprog_TMS*0.3)*0.0`
changes nothing numerically. The checked-in original
(`hd_mrgsolve_model.R`) is untouched and still carries both defects
exactly as written; an identically syntax-fixed scratch copy of the
original (same two fixes, original compartment/param names, never
committed) was built in-memory purely to construct the verification
comparison target — see `hd_refactor_notes.md`.

**Fix upstream would be:** convert the `$INIT` block to the `$MAIN`
`<CMT>_0` idiom, and either implement `mHTT_protein_reduction_f()`'s
intended DMT-driven mHTT-protein-reduction logic directly as DSL-local
C++ inside `$ODE`, or delete the dead term if it was never meant to do
anything beyond its own placeholder stub.

## 88. `huntingtons-disease/hd_mrgsolve_model.R`'s entire `$TABLE` block is never captured — none of its 17 derived quantities (including `TMS`/`TFC`/`MSN_pct`/`mHTT_total`/`BDNF_level`, which the file's own post-DSL R script plots and summarizes) are exposed as outputs

Found while adding `$CAPTURE`-equivalent output declarations for the six
refactored compounds' own `C_<STEM>`/`EFFECT_<STEM>` terms, per the fork
guide's qspserver-compatibility requirement #4. Every line in `$TABLE` is
written as a bare `double NAME = ...;` — no `capture` keyword and no
`$CAPTURE` block appear anywhere in the file. Confirmed via `POST
/model_manifest` on the untouched original (after fixing only the two
compile-blocking defects in issue #87, so it could build at all):
`outputPaths` lists exactly the 20 raw `$CMT` compartments and nothing
else — none of `HTBZ_conc`, `DTBZ_conc`, `VBZ_conc`, `ASO_CSF`,
`mHTT_mRNA_rel`, `mHTT_total`, `oligomer_pct`, `BDNF_level`, `DA_level`,
`MSN_pct`, `OxStress`, `Inflam`, `TMS`, `TFC`, `VMAT2_inh`,
`chorea_red_pct`, or `cUHDRS` is discoverable or returnable from
`/run_simulation`, even though every one of them is computed every
timestep. This is a broader, more severe variant of the `$CAPTURE`-
collision defect class already logged as issues #35/#43/#50/#59/#78/#80
(those are about a *duplicate* capture of an already-declared name
colliding; this file has *no* capture declaration at all for any derived
quantity). The practical effect: the original file's own post-DSL R
script — every column in its `summarise()` call (`TMS_baseline`,
`TMS_year5`, `TFC_change`, `MSN_yr5_pct`, `BDNF_yr5`, `mHTT_total_yr5`,
`chorea_red_pct`) and all six of its `ggplot()` calls (`y = TMS`,
`y = TFC`, `y = MSN_pct`, `y = mHTT_total`, `y = BDNF_level`) — would
error with "object not found" if actually run, because `mrgsim() %>%
as_tibble()` never contained those columns to begin with. As literally
authored, the model's own downstream analysis has never worked.

Not fixed in the original, per the never-edit-upstream rule. Fixed in
`hd_mrgsolve_model_refactored.R` by changing `double` to `capture` for
every one of the 17 pre-existing `$TABLE` lines (identical formulas, only
the declaration keyword changed) plus new `_OUT`-suffixed captures for
the refactor's own `C_TBZ`/`C_HTBZ`/`EFFECT_HTBZ`/`C_DTBZ`/`EFFECT_DTBZ`/
`C_VBZ`/`EFFECT_VBZ`/`C_TOMINERSEN`/`EFFECT_TOMINERSEN`/`C_RILUZOLE`/
`EFFECT_RILUZOLE` terms — disclosed in `hd_refactor_notes.md`, confirmed
discoverable via `/model_manifest`'s `outputPaths` (47 total) and
verified byte-for-byte reproducible against an identically-patched
scratch copy of the original.

**Fix upstream would be:** add `capture NAME;` (or an equivalent
`$CAPTURE` block) for every `$TABLE`-computed quantity the post-DSL R
script actually consumes, so the model's own analysis and plotting code
runs as written.

## 89. `huntingtons-disease/hd_mrgsolve_model.R`'s riluzole motor-benefit term never reads its own declared `Emax_riluzole` parameter — a hardcoded `0.20` is used instead

Found while isolating riluzole's PK/PD block for the same refactor. The
file declares `Emax_riluzole = 0.45 // max glutamate/excitotoxicity
reduction` in `$PARAM`, but the only place riluzole's motor effect is
computed reads a literal instead:
```
double riluzole_motor = (riluzole_brain > 0.01) ?
    IMAX(riluzole_brain, EC50_riluzole, 0.20) * 0.15 : 0.0;
```
Grepping the whole file for `Emax_riluzole` shows exactly one occurrence
— its own `$PARAM` declaration — never read anywhere else. This is the
same defect class already logged as issue #84 (vte's `EMAX_RIV`/
`EMAX_WARF`): a declared, physiologically-named ceiling parameter that a
reader would reasonably assume governs the effect it is named for, while
the formula that actually computes that effect uses an unrelated
hardcoded number instead. Confirmed upstream via `grep -n
"Emax_riluzole"` against the original (one match) and independently via
`/model_manifest` on the untouched (compile-fixed, see issue #87) original,
which lists `Emax_riluzole` as a live, overridable `$PARAM` entry that a
caller could reasonably expect to control riluzole's motor benefit but
which in fact does nothing.

Not fixed in the original, per the never-edit-upstream rule. In
`hd_mrgsolve_model_refactored.R`, `EMAX_RILUZOLE` is set to `0.03` — the
value the original's formula actually produces as its combined ceiling
(the hardcoded inner `0.20` times the also-hardcoded outer `0.15`
multiplier, `0.20*0.15=0.03`), matching issue #84's precedent of adding a
correctly-wired parameter under the canonical name carrying the value
that is actually used, rather than carrying forward the original's dead
`Emax_riluzole=0.45` under the `EMAX_RILUZOLE` slot. Disclosed in
`hd_refactor_notes.md`; verification (bespoke riluzole dosing scenario,
since none of the original's own seven named scenarios exercises
riluzole at all) confirms this reproduces the original's actual output
exactly (max abs diff 0.0).

**Fix upstream would be:** replace the hardcoded `0.20` in the
`IMAX(riluzole_brain, EC50_riluzole, 0.20)` call with `Emax_riluzole`
(and decide whether the outer `*0.15` scaling is still wanted, or should
be folded into a single, correctly-documented ceiling value).

## 90. `huntingtons-disease/hd_mrgsolve_model.R` doses branaplam into tominersen's own CSF compartment, contaminating an unrelated compound's PK state

Found while isolating tominersen's PK/PD block for the same refactor.
The file's `make_dose_events()` targets the CSF compartment
(`tominersen_CSF`, `$CMT` position 8) for **both** the tominersen
scenario and the entirely unrelated branaplam scenario:
```
} else if (scenario == "Tominersen_Q8W") {
    e_aso <- ev(cmt = 8, amt = 120, ii = 1344, ...)
    ...
} else if (scenario == "Branaplam_Q1W") {
    e_bran <- ev(cmt = 8, amt = 50, ii = 168, ...)
    ...
```
Branaplam has no PK compartment of its own anywhere in the file (its
effect is a bare on/off flag, `dose_branaplam`, multiplying a fixed
`branaplam_eff = 0.50`) — dosing it into `cmt = 8` deposits 50 mg every
168 h directly into the tominersen CSF compartment, which then decays via
tominersen's own `CL_tominersen/Vc_tominersen` elimination throughout the
`Branaplam_Q1W` scenario, even though `dose_tominersen = 0` for that
scenario. Confirmed upstream by grepping every `ev(cmt = ...)` call in
`make_dose_events()` (branaplam is the only scenario whose event targets
a compartment already owned by a different, correctly-modeled compound)
and by running the original's own `Branaplam_Q1W` scenario through the
qspserver `mrgsolve_api` container: `tominersen_CSF` (renamed `C_TOMINERSEN`
in the refactored sibling) rises from this spurious dosing exactly as
tominersen's own PK equation dictates, while `dose_tominersen` stays 0
throughout.

Not fixed in the original, per the never-edit-upstream rule. This is
also why `EFFECT_TOMINERSEN` in `hd_mrgsolve_model_refactored.R`
deliberately **keeps** the original's `(dose_tominersen > 0) ? ... : 0.0`
gate rather than driving purely off `C_TOMINERSEN` the way
`EFFECT_HTBZ`/`EFFECT_DTBZ`/`EFFECT_VBZ` now do (see
`hd_refactor_notes.md`) — dropping the gate would let this pre-existing
contamination bug silently start suppressing mHTT mRNA production during
the `Branaplam_Q1W` scenario, a numeric change the refactor must not
introduce. Verified: with the gate preserved, the `Branaplam_Q1W`
scenario matches the original exactly (max abs diff 0.0) despite the
spurious `CENT_TOMINERSEN` mass.

**Fix upstream would be:** give branaplam its own dosing compartment (or
model it as a flag-only intervention with no `ev()` at all, since its
effect is already flag-driven) instead of routing its dose through
tominersen's CSF compartment.

---

## 93. `graves-disease/gd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$PARAM F_MMI` collides with mrgsolve's own reserved `F_<compartment>` bioavailability symbol, and `$CAPTURE` re-lists eight `$CMT` compartment names directly

**Confirmed upstream**, via `POST /model_manifest` on the untouched original's
own DSL (extracted verbatim from `graves_model <- '...'`) through the local
qspserver `mrgsolve_api` container: HTTP 500,

```
invalid class "mrgmod" object: 1: reserved symbols in model names: F_MMI
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE:
MMI,PTU_C,RAI_THYR,PROP_C,TPO_inh,ThyMass,GO_act,BoneResor
```

Two independent defects:

1. **`F_MMI` reserved-symbol collision.** mrgsolve reserves `F_<cmt>` as the
   per-compartment bioavailability override for a compartment literally named
   `<cmt>`. The file declares a compartment named `MMI` (`$CMT ... MMI ...`)
   *and* a `$PARAM F_MMI = 0.93` (declared, and — separately — never actually
   referenced anywhere in `$ODE`, confirmed by exact-name grep). The two
   collide under mrgsolve 2.0.1's stricter reserved-symbol check.
2. **`$CAPTURE` duplicates eight compartment names.** `$CAPTURE ...
   MMI PTU_C RAI_THYR PROP_C TPO_inh ThyMass GO_act BoneResor` re-lists
   compartments that mrgsolve already exposes as output columns automatically
   — the same class of defect as `#78` (`narcolepsy`) and `#75`
   (`cytokine-release-syndrome`), just a different file.

**Fix upstream would be:** rename `F_MMI` to any non-colliding identifier (it
is unused regardless, so this changes no numeric behaviour), and drop the
second `$CAPTURE` line entirely (compartments do not need re-listing).

**Not fixed in the original**, per the never-edit-upstream rule. Fixed
directly in `gd_mrgsolve_model_refactored.R` as a disclosed, syntax-only,
non-numeric build-compatibility change: defect 1 disappears as a side effect
of this refactor's own naming-convention rename (the compartment becomes
`CENT_MMI`, no longer colliding with `F_MMI`); defect 2 is fixed by dropping
the redundant `$CAPTURE` line (replaced with the refactor's own new
`C_<STEM>`/`EFFECT_<STEM>` capture entries) — see `gd_refactor_notes.md` for
the full disclosure and the verification that proves this introduces no
numeric change.

---

## 94. `graves-disease/gd_mrgsolve_model.R` doses two flag-only interventions (exogenous levothyroxine, rituximab) into two other compounds' own PK compartments, contaminating their concentration states

Same class of defect as `#90` (`huntingtons-disease` branaplam-into-tominersen
contamination), found twice in this one file while refactoring its three
genuine drug compounds (MMI, PTU, propranolol):

1. **Scenario 5 ("Block-and-Replace")**: `br_lt4 <- ev(amt=0.10, ii=1.0,
   addl=730-1, cmt=5, time=90)` is commented `# LT4 start at day 90`, but
   `cmt=5` is `PROP_C` (propranolol's own PK compartment), not a dedicated
   levothyroxine pool — this model has **no LT4 compartment at all** despite
   the file's own header comment claiming `C5 = Levothyroxine (exogenous T4)
   pool` and a `USE_LT4` flag that is declared in `$PARAM` but never
   referenced anywhere in `$ODE` (confirmed by exact-name grep — `USE_LT4`
   only appears in its own `$PARAM` declaration and in the R-side
   `param(..., USE_LT4=1)` call for this one scenario). So this scenario's
   "LT4" dose is silently deposited into, and then cleared through, the
   propranolol PK equations instead.
2. **Scenario 6 ("Rituximab + MMI")**: `rtx_events <- ev(amt=1, cmt=1,
   time=0) + ev(amt=1, cmt=1, time=14)` doses two boluses into `cmt=1`
   (`MMI`), even though rituximab's entire effect in this model is a bare
   on/off flag (`RTX_eff = (USE_RTX > 0.5) ? RTX_Bcell : 0.0`) with no PK
   compartment of its own anywhere in the file. Confirmed upstream by
   grepping every `ev(cmt = ...)` call in the R driver code (rituximab and
   LT4 are the only two scenarios whose dosing targets a compartment already
   owned by a different, correctly-modeled compound) and by running the
   original's own `Rituximab + MMI` scenario through the qspserver
   `mrgsolve_api`: `MMI` (renamed `CENT_MMI` in the refactored sibling) rises
   transiently and decays via MMI's own `kel_MMI`/`CL_MMI` elimination from
   this spurious dosing, exactly as the compound's PK dictates, even though
   this scenario's own two `param(USE_MMI=1, ...)` intent is presumably to
   test rituximab's *own* effect, not to add extra antithyroid drug exposure.

Neither is one of this file's three census-row drugs' own defects — MMI,
PTU, and propranolol's PK/effect equations are otherwise unaffected — so
**not fixed as part of the refactor's own scope**: both scenarios are
preserved byte-for-behavior-identical in `gd_mrgsolve_model_refactored.R`
(same `cmt=` targets, same dosing), and the refactor's `Rituximab + MMI`
verification run confirms the spurious `CENT_MMI` contamination reproduces
identically between the original and the refactored DSL (max abs diff
~1e-27, floating-point-scale) — see `gd_refactor_notes.md`.

**Fix upstream would be:** give levothyroxine its own dosing compartment
(`GUT_LT4`/`CENT_LT4`) wired into a genuine exogenous-T4 contribution to
`T4_C`/`fT4_C`, and give rituximab either its own (flag-driven, no-`ev()`)
scenario or a dedicated dummy compartment that nothing else reads.

---

## 95. `cluster-headache/ch_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$PARAM @annotated` line `EMAX_TOPI` is missing the required description field

Confirmed via `POST /model_manifest` on the untouched original's own DSL
(extracted verbatim from `ch_code <- '...'`) through the local qspserver
`mrgsolve_api` container:

```
Error: improper annotation format
 input: EMAX_TOPI :  0.35
 context: parse annotated parameter block (PARAM)
```

An `@annotated` `$PARAM` line requires three colon-delimited fields
(`name : value : description`); this one line supplies only the first two
(`EMAX_TOPI :  0.35`, with no trailing `: description`), so mrgsolve's
annotation parser rejects the whole block at this line before reaching any
other parameter -- a single-line defect, confirmed by scanning every other
`$PARAM` line in the file for the same pattern (none found). This is inside
topiramate's own PK/PD block, one of the seven compounds this refactor
covers (Galcanezumab, Lithium, Prednisolone, Sumatriptan, Topiramate,
Verapamil PR, Zolmitriptan IN).

**Not fixed in the original** (`ch_mrgsolve_model.R`), per the never-edit-
upstream rule; it still carries the defect forward unfixed. Per the guide's
settled policy (a defect found inside the scope compound's own block is
folded into the refactor itself rather than treated as a separate build-
compat patch), the fix is applied directly in `ch_mrgsolve_model_refactored.R`
as part of topiramate's own naming-convention rename: a purely descriptive
third field was appended (`EMAX_TOPI :  0.35  : Max fractional attack-rate
reduction (preventive); ...`), changing nothing about the parameter's name
or value. Confirmed non-numeric: all seven of the original's own dosing
scenarios (S0-S6, 12-week horizon, matching `simulate_all()`'s own default
`weeks=12`/`delta=1`) were run through the qspserver `mrgsolve_api` against
both the untouched original (with only this one build-compat fix applied to
a throwaway verification copy, never to the checked-in file) and the
refactored DSL, and every shared `$CAPTURE`d output matched exactly (max abs
diff 0.0) -- see `ch_refactor_notes.md`.

**Fix upstream would be:** add a real description field to the `EMAX_TOPI`
line, e.g. `EMAX_TOPI :  0.35  : max fractional attack-rate reduction
(preventive)`.

---

## 96. `peripheral-arterial-disease/pad_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 20 compartments

Confirmed via `POST /model_manifest` on the untouched original's own DSL
(extracted verbatim from `pad_code <- '...'`) through the local qspserver
`mrgsolve_api` container:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: Duplicated model names: A_clopi C_clopi
  C_am A_asp C_asp A_tica C_tica A_riva C_riva A_cilo C_cilo A_atst C_atst
  Plt_agg Thrombin LDL_C ABI WalkDist Collat EF_idx PlaqueVol hsCRP
```

Same defect class as issues #29/#34/#36/#59/#76/#81/#82/#83: mrgsolve 2.0.1
treats `$INIT`'s plain `NAME = value` list as its own compartment
declaration, not a companion to `$CMT`, so having both blocks declare the
same 20 names is read as 20 duplicate compartments. Not specific to any of
this file's six modeled compounds (Clopidogrel, Aspirin, Ticagrelor,
Rivaroxaban, Cilostazol, Atorvastatin) — it sits entirely in the
boilerplate `$INIT` block shared by all of them plus the nine disease-PD
states.

**Not fixed in the original** (`pad_mrgsolve_model.R`), per the never-edit-
upstream rule; it still carries the defect forward unfixed. Per the guide's
settled policy, the fix is applied directly in
`pad_mrgsolve_model_refactored.R`: the `$INIT` block deleted, its 20
assignments moved into `$MAIN` as `<CMT>_0 = value;` under the refactor's
renamed compartments (e.g. `GUT_CLOP_0 = 0;`, `CENT_CLOPAM_0 = 0;`), same
values, no compartment added or removed — syntax-only, non-numeric.
Confirmed non-numeric: all seven of the original's own dosing scenarios
were run through the qspserver `mrgsolve_api` against both the untouched
original (with only this one build-compat fix applied to a throwaway
verification copy, never to the checked-in file) and the refactored DSL,
and every shared `$CAPTURE`d output matched exactly (max abs diff 0.0) over
a 60-day verification window — see `pad_refactor_notes.md` for why the
window is shortened.

**Fix upstream would be:** delete the `$INIT` block and set initial values
via `$MAIN`'s `<CMT>_0 = value;` idiom instead.

---

## 97. `peripheral-arterial-disease/pad_mrgsolve_model.R`'s LDL-C turnover parameters are inconsistent with their own stated half-life, causing LDL-C to grow without bound (thousands of mg/dL within months) rather than sit near its stated baseline

The `$PARAM` block declares:

```
k_ldl_prod: 0.0024    : LDL-C hepatic production rate (mg/dL/h)
LDL_base  : 130       : Baseline LDL-C (mg/dL)
k_ldl_cl  : 0.0000185 : LDL-C clearance rate (1/h; t1/2 ~55h)
```

and the turnover ODE is a plain zero-order-production /
first-order-elimination model:

```
double LDL_cl_total = k_ldl_cl * (1.0 + 3.0 * Eff_HMG) * LDL_C;
dxdt_LDL_C = k_ldl_prod * LDL_base - LDL_cl_total;
```

At the model's own baseline (`LDL_C = 130`, untreated so `Eff_HMG = 0`),
production is `k_ldl_prod * LDL_base = 0.0024*130 = 0.312 mg/dL/h` while
clearance is only `k_ldl_cl * LDL_C = 0.0000185*130 = 0.0024 mg/dL/h` —
production outweighs clearance by roughly 130-fold, so `LDL_C` is nowhere
near steady state at its own stated initial value and instead climbs
almost linearly (~7.4 mg/dL/day). Separately, `k_ldl_cl`'s own comment
claims a `t1/2 ~55h`, which would require `k_ldl_cl = ln(2)/55 ≈ 0.0126`,
not `0.0000185` (a value that alone implies a `t1/2` of roughly 4.3 years)
— the declared value and its own comment disagree by a factor of ~680,
consistent with a units/typo error rather than a deliberately slow
clearance.

Confirmed by direct simulation of the untouched original (build-compat-
fixed copy per issue #96, `mrgsolve_api`, no dosing): `LDL_C` reaches
~1005 mg/dL by day 121, ~3627 mg/dL by day 528, and the ODE solver fails
with a convergence error (`lsoda`, `istate -5`) or produces `NaN` shortly
after, well inside the model's own stated 2-year (730-day) default
simulation horizon and its own 180/365-day summary-table checkpoints —
i.e. the original's own published Figure 5 / summary table, if actually
run to completion, would report physiologically implausible LDL-C values
in the thousands of mg/dL rather than the intended statin-responsive
trajectory near 130 mg/dL. This is a disease-PD-side parameter defect, not
inside any of the six compounds' own PK/PD blocks (statin's `Eff_HMG` can
only ever quadruple the already-too-small clearance term, so it cannot
rescue the imbalance).

**Not fixed in the original** (`pad_mrgsolve_model.R`) or in
`pad_mrgsolve_model_refactored.R` — `k_ldl_cl`'s value is copied verbatim
in the refactor, per the never-invent-or-correct-a-parameter-value rule.
Worked around for verification purposes only by shortening the comparison
window to 60 days (safely before the divergence becomes numerically
severe) — see `pad_refactor_notes.md`.

**Fix upstream would be:** either raise `k_ldl_cl` to match its own stated
`t1/2 ~55h` (i.e. ~0.0126 /h) or lower `k_ldl_prod`/raise `LDL_base` so that
production and clearance balance near the intended 130 mg/dL baseline.

---

## 98. `stevens-johnson-syndrome-ten/sjsten_mrgsolve_model.R` computes a "continued culprit drug" toggle that is never actually used

Found while refactoring the culprit-drug PK block (stem `DRUG`) for this
file's PK/PD refactor. `[ODE]` opens with:

```
double drug_input = (PREDOSE > 0.5 && (SCEN_WD < 0.5 || SOLVERTIME < DAY_WD)) ? 1.0 : 0.0;
dxdt_A_drug_dep = -KA_drug*A_drug_dep;
dxdt_A_drug     =  KA_drug*A_drug_dep - (CL_drug/V_drug)*A_drug;
```

`drug_input` reads three patient/scenario covariates — `PREDOSE`
("continued culprit drug? 0/1"), `SCEN_WD` ("drug withdrawal toggle"), and
`DAY_WD` — that together look like they should gate a continued,
zero-order culprit-drug exposure after the withdrawal day (the file's own
comment above `[PARAM]` even lists "drug withdrawal (Supportive care
only)" as treatment scenario 1). But `drug_input` is never referenced
again anywhere in the file — grepping the whole DSL for `drug_input`
returns exactly the one line that declares it. The depot/central ODEs
right below it (`dxdt_A_drug_dep`/`dxdt_A_drug`) are driven purely by
first-order absorption from whatever bolus an external `ev()` puts into
`A_drug_dep`; `PREDOSE`, `SCEN_WD`, and `DAY_WD` have no effect on the
simulated trajectory at all, regardless of their values. Confirmed via
`/model_manifest` (all three list as live, overridable `$PARAM` entries)
and by re-running this file's own "Supportive" scenario (`SCEN_WD=1`,
the default) with `PREDOSE` and `DAY_WD` overridden to several different
values: the drug/antigen/T-cell/BSA/SCORTEN trajectory is bit-identical
in every case, confirming the three covariates are fully inert.

This is the same defect class already logged for other files in this
fork (e.g. #89's `Emax_riluzole`) — a declared, plausibly-named covariate
that a reader or API caller would reasonably expect to control the
model's behaviour, but which in fact does nothing.

**Not fixed in the original** (`sjsten_mrgsolve_model.R`), per the
never-edit-upstream rule; `drug_input`, `PREDOSE`, `SCEN_WD`, and
`DAY_WD` are all carried forward verbatim and unused in
`sjsten_mrgsolve_model_refactored.R` too, with an inline comment marking
the line as known-dead code rather than silently dropping it (dropping it
would be a silent behavioural edit disguised as a rename, since
`drug_input` computing to an unused value is not observably different
from not computing it at all, but the guide's rename-only discipline for
non-scope lines argues for leaving it untouched and disclosed instead).
See `stevens-johnson-syndrome-ten/sjsten_refactor_notes.md` for the full
disclosure.

**Fix upstream would be:** either wire `drug_input` into
`dxdt_A_drug`/`dxdt_A_drug_dep` as an actual continued zero-order exposure
term (consistent with the "continued culprit drug" comment), or remove
the three unused covariates and the dead computation if withdrawal is
meant to be modelled solely by stopping external dosing.

---

## 99. `type2-diabetes/t2dm_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CMT`+`$INIT` jointly redeclare all 25 compartments, `$CAPTURE` re-lists nine `$CMT` names directly, and `SOLVERTIME` is used inside `$MAIN` where it is out of scope

Confirmed via `POST /model_manifest` on the untouched original's own DSL
(extracted verbatim from `t2dm_model_code <- '...'`) through the local
qspserver `mrgsolve_api` container. Three independent, stacked defects,
each blocking the build in turn as the previous one is worked around:

**1. `$CMT`+`$INIT` duplicate declaration** (same defect class as issues
#29/#34/#36/#59/#76/#81/#82/#83/#96):

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: 1: Duplicated model names: MET_GUT MET_C
  MET_P EMPA_C SEMA_SC SEMA_C DPP4I_C SU_C INS_SC INS_C PIOG_C Gp Gt Ip
  X_action Gc GLP1 beta_mass IR_H IR_P FFA BW_t HbA1c_cmpt eGFR_cmpt UACR_cmpt
```

**2. `$CAPTURE` re-lists nine `$CMT` names directly** (same defect class as
issue #93):

```
invalid class "mrgmod" object: 2: compartment should not be in $CAPTURE:
  Gc,GLP1,IR_H,IR_P,FFA,BW_t,HbA1c_cmpt,eGFR_cmpt,UACR_cmpt
```

**3. `SOLVERTIME` used inside `$MAIN`.** After working around 1 and 2 in a
throwaway verification copy, the build still fails:

```
/usr/local/lib/R/site-library/mrgsolve/base/modelheader.h:112:20: error:
  '_ODETIME_' was not declared in this scope
  112 | #define SOLVERTIME _ODETIME_[0]
      |                    ^~~~~~~~~~~
303:65: note: in expansion of macro 'SOLVERTIME'
  303 | BW_sema_effect = USE_SEMA * Emax_sema_wt * (1 - exp(-ksema_wt * SOLVERTIME));
```

`_ODETIME_` is only declared inside the generated `$ODE` function; `$MAIN`
has no such variable in scope. Two lines in `$MAIN` reference it this way
(`BW_sema_effect`'s and `BW_piog_gain`'s exponential-decay weight-change
terms, in semaglutide's and pioglitazone's own effect blocks respectively)
— neither is inside any drug's PK block proper, both are downstream PD
terms reading elapsed time. `TIME` (mrgsolve's per-record simulation-time
covariate, valid in `$MAIN`) reads the identical value at every point this
model's own harness ever observes it (`mrgsim(end=, delta=)` generates one
observation record per requested output time, and `$MAIN` re-executes at
each one) — substituting `TIME` for `SOLVERTIME` is confirmed non-numeric
by the verification below.

**Not fixed in the original** (`t2dm_mrgsolve_model.R`), per the never-
edit-upstream rule; it still carries all three defects forward unfixed.
Per the guide's settled policy, all three fixes are applied directly in
`t2dm_mrgsolve_model_refactored.R`: the `$INIT` block deleted and its 25
assignments moved into `$MAIN` as `<CMT>_0 = value;` (under the refactor's
renamed compartments, e.g. `GUT_MET_0 = 0;`); the nine compartment names
removed from `$CAPTURE` (they remain in the output as ordinary compartment
states, just no longer double-declared); and both `SOLVERTIME` occurrences
replaced with `TIME`. Same values, no compartment or parameter added or
removed — syntax-only, non-numeric. Confirmed non-numeric: all six of the
original's own dosing scenarios (`metformin`, `combo_empa`, `combo_sema`,
`triple`, `insulin`, `dpp4i`, each with `build_doses()`'s own dose amounts/
intervals) plus the untreated baseline were run through the qspserver
`mrgsolve_api` against both the untouched original (with only these three
build-compat fixes applied to a throwaway verification copy, never to the
checked-in file) and the refactored DSL, full 365-day horizon at the
original's own `delta=6`, and every shared `$CAPTURE`d output and
compartment state matched exactly (max abs diff 0.0) in all seven runs —
see `t2dm_refactor_notes.md`.

**Fix upstream would be:** delete the `$INIT` block and set initial values
via `$MAIN`'s `<CMT>_0 = value;` idiom instead; drop the nine compartment
names from `$CAPTURE`; and change `SOLVERTIME` to `TIME` in `$MAIN`.

---

## 100. `type2-diabetes/t2dm_mrgsolve_model.R` models Insulin Degludec's own PK in full but never couples its concentration to any disease-PD term, and two of its seven advertised compounds (Glimepiride, Pioglitazone) are never dosed by any of the file's own scenarios

Found while refactoring all seven of this file's compounds' PK/PD blocks
into the fork's pluggable-PK naming convention.

**1. Insulin Degludec's concentration is computed and captured, but read by
nothing.** The file declares a full depot+central PK block for insulin
degludec (`INS_SC`→`INS_C`, `ka_ins`/`CL_ins`/`Vc_ins`), a `USE_INS` on/off
switch, a dedicated `Cins` concentration variable, and a "Metformin +
Insulin Degludec" treatment scenario with its own dosing — every
appearance of a seventh, fully-modeled antidiabetic drug. But grepping the
whole file for `Cins`/`INS_C`/`USE_INS` outside insulin's own PK block and
its own dosing/scenario wiring turns up nothing: no disease equation
(`EGP_val`, `Rd_val`, `X_action`, insulin secretion, or any other PD term)
ever reads `Cins`. The only place `CL_ins`/`Vc_ins` appear again is inside
`dxdt_Ip` (`(dIp_sec - CL_ins * Ip / Vc_ins) / Vi / BW_t * 10`), which
reuses those two constants as a generic clearance-rate denominator for
*endogenous* plasma insulin (`Ip`), applied unconditionally regardless of
`USE_INS` or any dose — i.e. exogenous insulin degludec's own simulated
plasma concentration has zero effect on glucose, insulin action, or any
other endpoint in this model; running the "Metformin + Insulin Degludec"
scenario changes nothing about the disease trajectory beyond what
metformin alone already does. Confirmed by running that scenario through
`mrgsolve_api` with `USE_INS` toggled 0 vs 1 (same dosing either way): the
`Gp`/`HbA1c_cmpt`/every other disease output is bit-identical between the
two runs, while `Cins`/`C_INS` itself rises and falls exactly as the PK
block dictates — the concentration exists and is captured, but is
provably inert.

**2. Glimepiride and Pioglitazone are fully parameterized (PK + a named
Hill effect term each: `E_su`/`Emax_su`/`EC50_su`; `E_piog_IR`/
`Emax_piog_IR`/`EC50_piog`) but neither compound is ever dosed.**
`build_doses()` only builds dosing for metformin, empagliflozin,
semaglutide, sitagliptin, and insulin degludec; there is no `su`/`piog`
branch, and none of the file's seven named scenario objects sets
`su=1`/`piog=1`. `USE_SU`/`USE_PIOG` therefore stay at their `$PARAM`
default of 0 in every scenario the file itself ever runs, and `SU_C`/
`PIOG_C` never receive a dose regardless. Confirmed by grepping
`build_doses()` and the `scenarios` list for `su`/`piog`/`SU_C`/`PIOG_C`:
the only matches are the dead compartment/parameter declarations
themselves. Two of the file's own header comment's seven advertised
compounds are therefore unreachable dead code paths in every scenario the
model ships with.

Neither defect is inside another compound's own block (unlike issue #94's
cross-compound dosing contamination) — each sits entirely within its own
named compound's block or the harness's own scenario list, so **not fixed
as part of the refactor's own scope**: `t2dm_mrgsolve_model_refactored.R`
preserves both exactly (Insulin Degludec's concentration remains exposed
as `C_INS` per the naming convention with no `EFFECT_INS` invented to fill
the gap; Glimepiride/Pioglitazone's Hill effect terms are renamed
`EFFECT_SU`/`EFFECT_PIOG_IR` and kept fully wired to their own PK, still
undosed by the shipped R harness). Verification for these two compounds
therefore could not reuse an existing scenario (per the guide's "run the
original's own scenarios" instruction) — a bespoke dosing pattern
(matching the style already used for the file's other direct-to-central
compounds) was built for verification purposes only, confirming the
renamed PK/PD blocks compute identically to the original's when actually
exercised; see `t2dm_refactor_notes.md`.

**Fix upstream would be:** for (1), add an exogenous-insulin contribution
to `X_action`/`Rd_val`/`EGP_val` (or wherever the model intends insulin
degludec's pharmacology to act) driven by `Cins`; for (2), add `su`/`piog`
dosing branches to `build_doses()` and at least one scenario that sets
`su=1`/`piog=1`, so both compounds' already-built PK/PD blocks are
actually exercised.

---

## 101. `acute-pancreatitis/ap_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` re-lists eighteen `$CMT` compartment names directly

Found while refactoring this file's eight compounds (indomethacin `IND`,
octreotide `OCT`, gabexate `GAB`, nafamostat `NAF`, ulinastatin `ULI`,
meropenem `MER`, anakinra `AKR`, fentanyl `FEN`) per
`FORK_WORKFLOW_GUIDE.md` Part 2. Same defect class as #30/#34/#41/#45/#46/
#78/#93: the model's `$CAPTURE` block lists names that are also `$CMT`
compartments —

```
$CAPTURE
TRYP TNF IL1 IL6 IL8 CRP NEC PERM GUT BT PF Cr BIL MAP GCS SOFA VAS
IND_C OCT_C GAB_C NAF_C ULI_C MER_C AKR_C FEN_C
SURVPROB BISAP MORT_HAZ
```

`TRYP, TNF, IL1, IL6, IL8, CRP, NEC, PERM, GUT, BT, PF, Cr, BIL, MAP, GCS,
SOFA, VAS, MORT_HAZ` (18 names) are all declared in `$CMT` earlier in the
same file (`IND_C`/`OCT_C`/etc. are `$GLOBAL` `#define` macros, not
compartments, so those are fine). `POST /model_manifest` on the untouched
original (via the qspserver `mrgsolve_api` container) fails outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  TRYP,TNF,IL1,IL6,IL8,CRP,NEC,PERM,GUT,BT,PF,Cr,BIL,MAP,GCS,SOFA,VAS,MORT_HAZ
```

Confirmed upstream: reproduces from the untouched original alone, unrelated
to any of the eight refactored compounds' own PK/effect blocks (all
eighteen duplicated names belong to the disease side of the model).

**Fix applied directly to the delivered `ap_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original, syntax-only
and non-numeric): the eighteen compartment names are dropped from the
`$CAPTURE` line — compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported.
The tracked `ap_mrgsolve_model.R` is completely untouched and still carries
the defect exactly as written. Confirmed by the exact-match (max abs diff
0.0, apart from one disclosed and unrelated single-row reporting artifact —
see `ap_refactor_notes.md`) verification results across all ten of the
file's own scenarios.

**Fix upstream would be:** remove the eighteen compartment names from the
`$CAPTURE` line (they need no explicit listing — mrgsolve reports every
`$CMT` state automatically).

## 102. `idiopathic-multicentric-castleman-disease/imcd_mrgsolve_model.R`'s ruxolitinib absorption ODE reads sirolimus's own `Ka_SIRO` parameter instead of its own declared `Ka_RUX`

Found while refactoring this file's nine compounds (siltuximab, tocilizumab,
sirolimus, rituximab, anakinra, ruxolitinib, CHOP-doxorubicin,
CHOP-cyclophosphamide, prednisone) per `FORK_WORKFLOW_GUIDE.md` Part 2. Same
class of defect as #90/#94 (one compound's dose/parameter silently
contaminating a different compound's own PK state), here a parameter
cross-reference rather than a dosing-target mixup. `$PARAM` declares both
`Ka_SIRO : 2.0 : 1/h` (sirolimus) and `Ka_RUX : 3.0 : 1/h` (ruxolitinib), but
`$ODE` reads:

```
dxdt_RUX_GUT = -Ka_SIRO * 24 * RUX_GUT;
dxdt_RUX_C   =  Ka_SIRO * 24 * RUX_GUT * F_RUX - k10_RUX * RUX_C;
```

`Ka_RUX` is declared in `$PARAM` but never referenced anywhere in `$ODE`
(confirmed by exact-name grep) — ruxolitinib's gut-to-central absorption is
actually governed by sirolimus's own 2.0/h rate, not its own declared 3.0/h.
Confirmed this changes real numeric behaviour (not a merely-cosmetic naming
slip): with `Ka_RUX` unused, ruxolitinib's own absorption half-life is set
by whatever value `Ka_SIRO` happens to hold, so changing sirolimus's
absorption rate in a future edit would silently also change ruxolitinib's,
an unintended coupling between two otherwise-independent drugs in the same
model.

**Not fixed here** (never edit an upstream-tracked file) — preserved
byte-for-behaviour-identical in `imcd_mrgsolve_model_refactored.R`: the
refactored file's `dxdt_GUT_RUX`/`dxdt_CENT_RUX` still read `KA_SIRO`, and
`KA_RUX` is carried forward as a declared-but-unused parameter (value `3.0`,
unchanged), exactly mirroring the original. Verified via qspserver
`mrgsolve_api`: a bespoke ruxolitinib-alone dosing scenario (20 mg PO BID —
no scenario in the original's own R code doses ruxolitinib at all) run
through both the original and the refactored DSL reproduces `Crux`/`C_RUX`
and every downstream disease output bit-for-bit (max abs/rel deviation
`0.0`), confirming the cross-wired absorption rate carries through
identically — see `imcd_refactor_notes.md`.

**Fix upstream would be:** replace `Ka_SIRO` with `Ka_RUX` in both of
ruxolitinib's own `$ODE` lines.

## 103. `idiopathic-multicentric-castleman-disease/imcd_mrgsolve_model.R` models Anakinra's own PK in full but never couples its concentration to any disease-PD term

Same class of finding as #100 (Insulin Degludec in `type2-diabetes`). The
file's header comment and `$CMT`/`$PARAM` blocks advertise anakinra
(IL-1 receptor antagonist) as one of nine treated compounds, and its full
depot+central PK is simulated (`ANA_SC`, `ANA_C`, `CL_ANA`, `V_ANA`, `F_ANA`,
`Ka_ANA`) and even scaled to a concentration for reporting (`Cana = ANA_C /
V_ANA * 1000`, listed in `$CAPTURE`) — but `Cana` is never read by any
disease/PD equation anywhere in `$ODE` (confirmed by exact-name grep: it
appears only in its own declaration and in `$CAPTURE`). Anakinra is dosed in
the original's own Scenario 7 ("Siltuximab + Sirolimus + Anakinra") and its
PK trajectory is fully simulated and reported, but has zero effect on IL-6,
CRP, lymph node size, or any other disease-side output in this model — a
"treatment" scenario that is pharmacokinetically but not pharmacodynamically
different from siltuximab+sirolimus alone.

**Not fixed here** — preserved as-is in `imcd_mrgsolve_model_refactored.R`:
`C_ANA` is renamed and exposed per the fork's naming convention (still
useful as a redirectable PK covariate for future work), but no `EFFECT_ANA`
term was invented, since the original has none to rename or refit. Verified
via qspserver `mrgsolve_api` (original Scenario 7 dosing): `Cana`/`C_ANA`
and every disease-side output match bit-for-bit between the original and the
refactored DSL, confirming anakinra's absence of PD coupling is reproduced
exactly, not silently patched over — see `imcd_refactor_notes.md`.

**Fix upstream would be:** either couple `Cana` into a real IL-1-blockade
mechanism (e.g. an `EFFECT_ANA` term on IL-6/CRP or plasmablast turnover) or
remove anakinra's PK block if it is not intended to have a modeled effect.

## 104. `mitral-regurgitation/mr_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: local `$ODE` scratch variables shadow `$TABLE`/`$CAPTURE` names, one local is declared twice in sibling loop scopes, a `$PARAM` name collides with mrgsolve's own auto-generated compartment-init symbol, and several multi-variable `double a=.., b=..;` declarations lose their `double` keyword

Four independent, unrelated build defects, all confirmed live against the
untouched original via the qspserver `mrgsolve_api` container
(`POST /model_manifest`):

1. **`$ODE`-local names shadow `$TABLE`-declared `$CAPTURE` names.** mrgsolve
   compiles `$ODE` and `$TABLE` into the same C++ scope, so `$TABLE`'s own
   `double NAME = g_NAME;` (the mechanism that actually populates a
   `$CAPTURE` output) collides with an *unrelated*, purely-internal `$ODE`
   scratch variable that happens to share the same name: `double Areq = Ann
   * (...)` (line 554) vs. `$TABLE`'s `double Areq = g_Areq;` (line 856);
   likewise `dPmv`, `Pes`, `RVol`, `SVtot`, `Ped`, `MAP`, `vwave`, `Ppcw`,
   `CI` — ten names total, each declared once inside the fast-block
   circulation solve in `$ODE` for a mid-computation quantity, and again in
   `$TABLE` to expose the corresponding `$CAPTURE` output. `gcc` reports
   `redefinition of 'double {anonymous}::Areq'` etc.
2. **One name is declared twice within `$ODE` itself, in sibling (not
   nested) block scopes**, which mrgsolve's build apparently does not
   scope the way ordinary C++ would: `double dP = c_s * m - Pmid;` inside
   the inner bisection `for` loop (line 587) and `double dP = c_s * Pe_ -
   Pmid;` immediately after that loop closes, still inside the outer `for`
   loop (line 593) — same "redefinition" class of error, confirming the
   build does not give the first `for` loop's body its own scope for this
   purpose.
3. **A `$PARAM` name collides with mrgsolve's own auto-generated
   compartment-initial-value symbol.** `$CMT` declares a compartment named
   `EROApri`, for which mrgsolve auto-generates a writable reference
   `EROApri_0` (the `$MAIN` idiom for setting that compartment's initial
   value). The original *also* declares its own `$PARAM` named literally
   `EROApri_0` ("Convenience copy of the initial degenerative orifice"),
   and `$MAIN` writes `EROApri_0 = EROApri_0;` — a self-assignment that
   was clearly meant to seed the compartment's initial value from the
   parameter, but the two meanings of the identical name collide: `gcc`
   reports `conflicting declaration 'double& EROApri_0'` (previously
   declared `const double&`) and `assignment of read-only reference
   'EROApri_0'`.
4. **Several multi-variable `double a = x, b = y;` declarations lose their
   `double` keyword during mrgsolve's build**, converting the whole
   statement into bare (undeclared) assignments: `double Pla = 0.2, Pla_hi
   = 140.0;` (line 569), `double Ved=0, Pes=0, SVf=0, RVol=0, SVtot=0,
   Ped=0, Q=0, MAP=0, Ppa=0, Psv=0, Vrved=0;` (line 570), `double SVt_ =
   60.0, Vd_ = 100.0, Pe_ = 100.0, SVf_ = 60.0, RV_ = 0.0, Pd_ = 8.0;`
   (line 575), and `double lo = 0.05, hi = 500.0;` (line 584) — confirmed
   because `gcc`'s error for each cascades to "was not declared in this
   scope" for every name after the first comma, and the compiler's own
   quoted source line for the error is missing the leading `double` that
   is plainly present in the checked-in file.

None of these four defects touch any of this file's eight drug PK/PD
blocks (furosemide, sacubitril, valsartan, beta-blocker, spironolactone/
canrenone, SGLT2 inhibitor, nitroprusside, dobutamine) — all four are
inside the shared haemodynamic/circulation and remodelling machinery that
every compound's effect ultimately reads from or writes into.

**Fix applied directly to the delivered `mr_mrgsolve_model_refactored.R`**
(syntax-only, non-numeric, per the guide's settled policy): (1) the ten
`$ODE`-local scratch variables that shadowed `$TABLE`/`$CAPTURE` names were
renamed with an `_l` suffix (`Areq_l`, `dPmv_l`, `Pes_l`, `RVol_l`,
`SVtot_l`, `Ped_l`, `MAP_l`, `vwave_l`, `Ppcw_l`, `CI_l`) everywhere they
are used as internal mid-computation quantities within `$ODE`; `$TABLE`'s
own declarations (the true `$CAPTURE`-populating ones) are untouched. (2)
The second `dP` declaration (line 593 and its one downstream use) was
renamed `dP2`. (3) The colliding `$PARAM` was renamed `EROApri_init`, and
`$MAIN` now reads `EROApri_0 = EROApri_init;` (the one usage-example
comment referencing the old parameter name was also updated for accuracy).
(4) The four multi-variable declarations were split into one `double X =
value;` statement per variable, each keeping its original initial value.
All four fixes are purely mechanical renames/declaration-splits with
identical numerics; see `mr_refactor_notes.md` for the full accounting and
the verification results (exact match, max abs diff 0.0, across two
scenarios exercising all eight compounds).

**Fix upstream would be:** apply the same four fixes to the checked-in
`mr_mrgsolve_model.R` directly.

---

## 105. `essential-tremor/et_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: a malformed `@annotated` line, a reserved parameter name, and an undeclared variable used in `$MAIN`

Found while refactoring this file's eleven compounds (propranolol `PRP`,
atenolol `ATN`, nadolol `NAD`, primidone `PRM`, phenobarbital `PB`, PEMA
`PEM`, topiramate `TOP`, gabapentin `GBP`, ethanol `ETH`, 1-octanol `OCT`,
T-type/Cav3.1 blocker `TTB`) per `FORK_WORKFLOW_GUIDE.md` Part 2. `POST
/model_manifest` on the untouched original (via the qspserver
`mrgsolve_api` container) fails with three independent defects, found one
at a time as each was fixed and the next surfaced:

1. **Malformed `$PARAM @annotated` line.** `EMAX_DBS: 0.90` has only two
   colon-separated fields (name, value); every other line in that block
   has three (name, value, description). mrgsolve 2.0.1's annotated-block
   parser requires exactly three fields per line and rejects the whole
   block:
   ```
   Error: improper annotation format
    input: EMAX_DBS: 0.90
    context: parse annotated parameter block (PARAM)
   ```
   Unrelated to any of the eleven compounds — `EMAX_DBS` belongs to the
   DBS/lesion (surgical) section of the model.

2. **`EPS` is a reserved word.** After fixing (1), the next failure is
   ```
   invalid class "mrgmod" object: Reserved words in model names: EPS
   ```
   `EPS` (`$PARAM ... EPS : 1e-4 : Noise seeding of the limit cycle`,
   used three times in `$ODE` as the oscillator's noise-seeding term)
   collides with a name mrgsolve reserves internally. Also unrelated to
   any of the eleven compounds — it is the tremor oscillator's own noise
   floor, not a drug parameter.

3. **`F_PRP` is used in `$MAIN` but never declared anywhere.** After fixing
   (1) and (2), the next failure is a C++ compile error:
   ```
   438:1: error: 'F_PRP' was not declared in this scope; did you mean 'F_PRM'?
     438 | F_PRP = F0_PRP + FMX_PRP*DPRP/(FD50_PRP + DPRP);
   ```
   `$MAIN` computes propranolol's dose-dependent first-pass-escape
   bioavailability into `F_PRP` (then copies it to `F_A_PRPG`, the
   mrgsolve bioavailability idiom for the gut depot), but `F_PRP` was
   never added to the file's own `$GLOBAL` block of shared `double`
   declarations (`PC0, G00, MU00, GH0, GV0, C_PRP, ...`) the way every
   other cross-block variable in this file is. Unlike (1) and (2), this
   one **is** inside a compound's own block (propranolol) rather than
   incidental to it.

Confirmed upstream: all three reproduce from the untouched original alone,
in this order, each only surfacing once the previous one is fixed.

**Fix applied directly to the delivered `et_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original — and, for
(3) specifically, per the guide's "if a defect is genuinely inside the
scope compound's own block... fix it as part of the refactor itself"
allowance, since it sits inside propranolol's own PK block): (1) the
missing description field was added, `EMAX_DBS: 0.90 : DBS max block
fraction`; (2) `EPS` was renamed `EPS_NOISE` everywhere (one declaration,
three uses); (3) `F_PRP` was added to the existing `$GLOBAL` double-
declaration list (`double PC0, G00, MU00, GH0, GV0, C_PRP, F_PRP;`). All
three are syntax-only, non-numeric, name/declaration-only changes — no
value, formula, or behaviour changed. Confirmed by exact-match (max abs
diff 0.0) verification across eleven scenarios (nine of the file's own
native dosing scenarios spanning all eleven compounds, plus two ad hoc
sanity doses for topiramate and gabapentin, which the original never
doses in any of its own scenarios) — see `et_refactor_notes.md`.

The tracked `et_mrgsolve_model.R` is completely untouched and still
carries all three defects exactly as written.

**Fix upstream would be:** apply the same three fixes (add the missing
annotation field, rename `EPS`, declare `F_PRP`) to the checked-in
`et_mrgsolve_model.R` directly.

## 106. `celiac-disease/cd_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: `$CAPTURE` re-lists four `$INIT` compartment names directly

Found while refactoring this file's one PK/PD-refactor census row (generic
"DRUG", identified as larazotide/AT-1001 per the file's own `$PARAM`
comments — see `cd_refactor_notes.md`). Same defect class as #30/#34/#41/
#45/#46/#78/#93/#101: the model's `$CAPTURE` block lists names that are
also compartments declared in `$INIT @annotated` —

```
$CAPTURE VH_CD_ratio Marsh_score Serology_pos Hgb_g_dL Ferritin_ug
         BMD_Tscore IEL_elevated GFD_flag Drug_Cp AntiTTG VH CrD AbsArea
```

`AntiTTG`, `VH`, `CrD`, `AbsArea` are all declared as compartments in the
file's own `$INIT @annotated` block. `POST /model_manifest` on the
untouched original (via the qspserver `mrgsolve_api` container) fails
outright:

```
Error in validObject(.Object) :
  invalid class "mrgmod" object: compartment should not be in $CAPTURE:
  AntiTTG,VH,CrD,AbsArea
```

Confirmed upstream: reproduces from the untouched original alone, unrelated
to the drug PK/PD block in scope for this refactor (all four duplicated
names belong to the disease side of the model — serology and
histopathology state variables).

**Fix applied directly to the delivered `cd_mrgsolve_model_refactored.R`**
(per the guide's settled policy for a non-compiling original, syntax-only
and non-numeric): the four compartment names are dropped from the
`$CAPTURE` line — compartments are always present in mrgsolve's output
regardless of `$CAPTURE`, so this changes nothing about what is reported.
The tracked `cd_mrgsolve_model.R` is completely untouched and still carries
the defect exactly as written. Confirmed by exact-match (max abs diff 0.0)
verification across all six of the file's own scenarios (untreated,
strict GFD, partial-GFD leak, and all three drug arms — larazotide,
ZED1227, AMG714 — full 1-year/366-point horizon) — see
`cd_refactor_notes.md`.

**Fix upstream would be:** remove the four compartment names (`AntiTTG`,
`VH`, `CrD`, `AbsArea`) from the `$CAPTURE` line (they need no explicit
listing — mrgsolve reports every `$INIT`/`$CMT` state automatically).

---

## 107. `hereditary-spherocytosis/hsph_mrgsolve_model.R` does not compile under mrgsolve 2.0.1: a stray duplicate `@annotated` marker, several multi-variable `double a=.., b=..;` declarations lose their `double` keyword, nine bare `N_1`..`N_9` locals are used with no declaration at all, and mrgsolve's cross-block variable hoisting collides ten local names between `$MAIN`/`$ODE`/`$TABLE`

Found while refactoring Mitapivat's PK/PD block per this file's
`driver-patches/data/compound_perturbation_census.md` row (classified
"Normalize duplicate concentration sites, then redirect"). All four defect
classes reproduce identically from the untouched original via the
qspserver `mrgsolve_api` container (`POST /model_manifest`) and are
unrelated to mitapivat's own PK/PD, which lives entirely inside its own
`$PARAM`/`$CMT`/`$ODE` lines.

1. **A stray, orphaned second `@annotated` marker sits inside the `$PARAM`
   block**, right after its last parameter (`k_don`) and before `$CMT`
   begins — a duplicate of the one already on the `$PARAM @annotated`
   header line, with no parameter name or value attached to it. mrgsolve
   2.0.1 rejects this outright, before any C++ compilation is even
   attempted: `Error: Found duplicated block option names.`
2. **Nine `N_1`..`N_9` locals are read/written in `$MAIN` with no
   declaration anywhere in the file** (`N_1 = 5.0*w[0]/wsum;` etc., then
   read back into `NA1_0 = N_1*Aent;` and friends): `error: 'N_1' was not
   declared in this scope; did you mean 'NB1'?`, and identically for
   `N_2`..`N_9`. Not part of either defect below — a plain missing
   declaration.
3. **Several multi-variable `double a = x, b = y, ...;` declarations lose
   the `double` keyword from every name after the first** during
   mrgsolve's build (the same mechanism already logged for
   `diabetic-ketoacidosis` (#37) and `mitral-regurgitation` (#104)):
   `double wsum = 0.0, w[NCOH];` (`$MAIN`); `double A[NCOH], V[NCOH],
   H[NCOH], Bq[NCOH], Dc[NCOH], sph[NCOH];`, `double pslow[NCOH],
   tauc[NCOH], Rcord[NCOH], kv[NCOH];`, `double vA[NCOH], vV[NCOH],
   vH[NCOH], vB[NCOH];`, `double hgeom[NCOH], hops[NCOH], hliv[NCOH],
   hlys[NCOH], hsen[NCOH], hz[NCOH];`, `double Ntot = 0.0, HBt = 0.0, VOLt
   = 0.0;`, `double destN = 0.0, destHb = 0.0, lysHb = 0.0, Wspl = 0.0,
   Rpool = 0.0;`, `double inN[NCOH], inA[NCOH], inV[NCOH], inH[NCOH],
   inB[NCOH];`, `double dN[NCOH], dNA[NCOH], dNV[NCOH], dNH[NCOH],
   dNB[NCOH];` (all `$ODE`); `double Nt = 0.0, HBs = 0.0, VLs = 0.0, dN_ =
   0.0, dHb_ = 0.0;`, `double wA = 0.0, wV = 0.0, wD_ = 0.0, wS = 0.0, wI =
   0.0, wR = 0.0;`, `double gG = 0.0, gO = 0.0, gS = 0.0, gL = 0.0, gY =
   0.0;`, `double Ai = NAq[i]/n, Vi = NVq[i]/n, Hi = NHq[i]/n, Bi =
   NBq[i]/n;`, `double Di = dcrit(Ai, Vi), si = Vi/vsph(Ai);` (all
   `$TABLE`) — fourteen lines total. `gcc` reports "was not declared in
   this scope" for every name after the first comma on each line, e.g.
   `error: 'V' was not declared in this scope`.
4. **mrgsolve hoists every `$MAIN`/`$ODE`/`$TABLE` block-local `double`/
   `int` declaration into one shared anonymous-namespace scope**, so the
   same scratch name used more than once anywhere across those three
   blocks — even inside separate, non-overlapping `if`/`for` braces —
   collides (same mechanism as #104's `Areq`/`dP` case, and the general
   pattern in #45's `k_phos_eff`): `i` (used as a `for`-loop counter once
   each in `$MAIN`, twice more in `$ODE`, once in `$TABLE` — five sites,
   one shared name); `ph` (`$ODE`'s transfusion-phase local and its
   mitapivat dosing-phase local, declared independently a few hundred
   lines apart); and `n`, `ex`, `vk`, `rho`, `ce`, `ig`, `fo`, `sl` (each
   declared once in `$ODE`'s per-cohort loop and again, independently, in
   `$TABLE`'s own duplicate per-cohort loop that recomputes the same
   quantities from live state). `gcc` reports, e.g., `error: redefinition
   of 'double {anonymous}::ph'` / `note: 'double {anonymous}::ph'
   previously declared here`.

**Fix applied directly to the delivered
`hsph_mrgsolve_model_refactored.R`** (per the guide's settled policy for a
non-compiling original, syntax-only and non-numeric, confirmed by exact
verification below): (1) the stray `@annotated` line deleted; (2) `double
N_1;` .. `double N_9;` declared (one per line) before first use; (3) every
multi-variable declaration listed above split into one `double NAME =
value;` statement per variable, values unchanged; (4) the five colliding
`i` for-loop counters renamed to five distinct names (`$MAIN`'s -> `im`;
`$ODE`'s first loop keeps `i`, its second and third renamed to `i2`/`i3`;
`$TABLE`'s renamed to `it`), `$ODE`'s transfusion-phase `ph` renamed
`ph_tx` (mitapivat's own kept as `ph`, being the compound in scope), and
`$TABLE`'s per-cohort-loop `n`/`ex`/`vk`/`rho`/`ce`/`ig`/`fo`/`sl` each
renamed with a `_t` suffix (`$ODE`'s copies of the same names left
unchanged). `$MAIN`'s and `$ODE`'s two independent, same-formula
`Aent` computations were also given distinct names (`$MAIN` keeps `Aent`,
`$ODE`'s renamed `Aent_ode`) since they hit the same hoisting collision.
All are pure symbol renames or declaration splits with identical
per-iteration math; see `hsph_refactor_notes.md` for the full accounting.
Verified via qspserver `mrgsolve_api`: exact match (max abs/rel deviation
0.0) across every one of the original's own `$CAPTURE`d outputs, over the
full 1400-day/1401-point time grid, in both a mitapivat-dosed (100 mg BID,
the file's own Shiny-app scenario value) and a no-drug run.

**Fix upstream would be:** delete the stray `@annotated` line; declare
`N_1`..`N_9` before use; split each multi-declarator line above into one
`double` statement per variable; and rename the eleven colliding scratch
names (`i`, used as five separate for-loop counters, plus `ph`, `n`,
`ex`, `vk`, `rho`, `ce`, `ig`, `fo`, `sl`, `Aent`) so no bare scratch name
is reused across `$MAIN`/`$ODE`/`$TABLE`.
