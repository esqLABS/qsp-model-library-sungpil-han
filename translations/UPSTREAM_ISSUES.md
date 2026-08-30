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
