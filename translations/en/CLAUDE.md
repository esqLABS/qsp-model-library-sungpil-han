# CLAUDE.md — QSP Disease Model Library

> English translation of the upstream [`CLAUDE.md`](../../CLAUDE.md), maintained in
> this fork. The original mixes Korean and English; this file is the whole document
> in English. See [`translations/README.md`](../README.md) for why it lives here
> rather than as an edit to the original.
>
> **This is a translation, not the operative instruction file.** Claude Code loads
> the upstream `CLAUDE.md` at the repository root; editing this copy changes
> nothing about how a session behaves.

## Project Purpose

This repository is a living library to which a disease-specific Quantitative
Systems Pharmacology (QSP) model is added automatically every day through
**Claude Code Routine (CCR)**. In each session Claude picks one disease, writes the
mechanistic map, the mrgsolve model, the Shiny app, and the references in full, then
commits and pushes directly to the `main` branch.

---

## Session Stop Requirements

The stop hook (`~/.claude/stop-hook-git-check.sh`) enforces that every session ends
with:
1. All changes committed
2. All commits pushed to the remote branch

The session will not terminate until these conditions are met.

---

## Guidelines for Adding a New Model

### Disease Selection
- Pick one from the disease list below (`category / disease name / pathogenesis`)
  that has **not been built yet**.
- Use today's date and what was built previously to pick from **a different category
  each day**.
- For diseases already added, consult the table in `README.md`.

### Directory Convention
- Lower case plus hyphens: `disease-name-in-english/` (e.g. `iga-nephropathy/`,
  `crohn-disease/`).
- No duplicate directories — check with `ls` before creating one.

### File Naming
| File | Convention |
|------|------|
| DOT  | `<abbr>_qsp_model.dot` or `<abbr>_qsp.dot` |
| SVG  | `<abbr>_qsp_model.svg` |
| PNG  | `<abbr>_qsp_model.png` (150 dpi, `dot -Tpng -Gdpi=150`) |
| R model | `<abbr>_mrgsolve_model.R` |
| Shiny | `<abbr>_shiny_app.R` (or `shiny_app/app.R`) |
| References | `<abbr>_references.md` |

### Model Quality Standards
1. **Mechanistic map (.dot)**: 100 or more nodes, at least 8 subgraph clusters, drug
   PK/PD included.
2. **mrgsolve model (.R)**:
   - At least 15 ODE compartments (drug PK + disease PD)
   - At least 5 treatment scenarios
   - Includes a note on parameter calibration against key clinical trial data
3. **Shiny app**: at least 6 tabs (patient profile · PK · key PD measures · clinical
   endpoints · scenario comparison · biomarkers).
4. **References**: at least 30 PubMed links, grouped by section.

### README Update

⚠️ **Preventing a repeat of the 2026-07-02 incident — read this.** Detailed
per-disease sections used to be appended below the table. They accumulated until the
file had swollen to thousands of lines, and eventually the citations and bold markup
left in some rows set with unmatched `**` pairs, **breaking the rendering of the
table itself**. The rules below exist to prevent that happening again and must be
followed without exception.

- **Never add a per-disease detail section below the table.** The README must end
  with "(1) intro sections → (2) model gallery table → (3) disclaimer / references /
  license". Detailed descriptions of a disease belong only in
  `<disease-dir>/README.md` (the per-directory README).
- When adding a new row to the table (after the last row; numbering is always
  consecutive):
  - **Write the summary cell as a single sentence**, roughly 190 characters or less.
  - **Do not use bold (`**`) or italics (`*`) inside the summary cell.** (The
    `[**Korean name**...]` in the model name cell is the exception — that is always
    exactly one pair.) If the count of `**` in a single cell becomes odd, the
    rendering of the entire table breaks from that row onward.
  - **Do not put citations, author names, years, journal names, or PMIDs in the
    summary cell.** (Anything like `(Author 2020 Journal[PMID 12345])` is
    forbidden — the evidence belongs in `<abbr>_references.md`.)
  - Keep the English subtitle (`<sub>...</sub>`) in the form `English Name · ABBR`
    only; do not append a string of mechanism and drug keywords.
  - Table image: use the format
    `<a href="path/to/svg"><img src="path/to/png" width="190" alt="ABBR"></a>`
    (follow the same pattern as the existing rows exactly).
- **After adding the row, and before committing, always run the cleanup script**:
  ```bash
  python3 scripts/fix_readme_table.py
  ```
  This script automatically normalises the formatting of the row just added (removing
  bold, removing citations, capping length, tidying the subtitle, normalising the
  category), renumbers the table, refreshes all the per-category counts and the
  "N models" wording, and restores the trailing sections (disclaimer / references /
  license) if they are missing. The script only normalises table formatting; it does
  not touch the content (disease descriptions) of other rows.
- After running the script, confirm that
  `python3 scripts/fix_readme_table.py --check` prints `PASS`, then commit. If it
  prints `FAIL` (one of: a break in the table numbering, a column count mismatch, an
  unmatched `**` pair, a broken link, or a missing trailing section), fix the cause
  and check again.
- Do not modify the content of other existing rows (formatting normalisation by the
  cleanup script excepted).

### Commit & Push
```bash
git add -A
git commit -m "Add <Disease> QSP model: mechanistic map, mrgsolve ODE, Shiny app, references"
git push -u origin HEAD
```
Do not create a PR; merge directly into `main`, or proceed as the user directs.

---

## Technology Stack

| Tool | Purpose |
|------|------|
| **Graphviz** (`dot`) | Rendering the mechanistic map (.dot to .svg/.png) |
| **mrgsolve** (R) | ODE-based PK/PD model |
| **Shiny** (R) | Interactive dashboard |
| **Claude Code Routine** | Daily automated model generation and commit |

---

## Disease List

Pick from the list below in each session. For entries already added, consult the
README table.

### Autoimmune Disease
Rheumatoid arthritis · Systemic lupus erythematosus · Sjögren syndrome · Systemic
sclerosis · Polymyositis · Dermatomyositis · Behçet disease · Ankylosing spondylitis
· Psoriatic arthritis · Reactive arthritis · Mixed connective tissue disease (MCTD)
· Antiphospholipid syndrome · Granulomatosis with polyangiitis (GPA) ·
Eosinophilic granulomatosis with polyangiitis (EGPA) · Microscopic polyangiitis
(MPA) · Polyarteritis nodosa (PAN) · Takayasu arteritis · Giant cell arteritis ·
IgA vasculitis · Relapsing polychondritis · Crohn disease · Ulcerative colitis ·
Autoimmune hepatitis · Primary biliary cholangitis (PBC) · Primary sclerosing
cholangitis (PSC) · Coeliac disease · Autoimmune pancreatitis · Pernicious anaemia ·
Type 1 diabetes · Hashimoto thyroiditis · Graves disease · Addison disease ·
Multiple sclerosis · Myasthenia gravis · Neuromyelitis optica (NMO) · Guillain-Barré
syndrome · CIDP · Autoimmune encephalitis · ITP · AIHA · Evans syndrome · Pemphigus
vulgaris · Bullous pemphigoid · Vitiligo · Alopecia areata · Goodpasture syndrome ·
**IgA nephropathy** ✓ · Adult-onset Still disease

### Chronic Disease
Type 2 diabetes · Dyslipidaemia · Metabolic syndrome · Obesity · Osteoporosis ·
Paget disease · Gout · Pseudogout · Polycystic ovary syndrome (PCOS) · Acromegaly ·
Primary hyperparathyroidism · Chronic hypothyroidism · **Essential hypertension** ✓
· Stable angina · Heart failure (HFrEF) · Heart failure (HFpEF) · Atrial
fibrillation · Peripheral arterial disease (PAD) · Chronic venous insufficiency ·
Abdominal aortic aneurysm · Hypertrophic cardiomyopathy (HCM) · Dilated
cardiomyopathy (DCM) · COPD · Bronchial asthma · Bronchiectasis · **Pulmonary
arterial hypertension (PAH)** ✓ · Idiopathic pulmonary fibrosis (IPF) ·
Pneumoconiosis · Obstructive sleep apnoea (OSA) · Sarcoidosis · GERD · Chronic
gastritis · Peptic ulcer · NAFLD · Liver cirrhosis · Chronic hepatitis B · Chronic
pancreatitis · Irritable bowel syndrome (IBS) · Diverticular disease ·
Cholelithiasis · Chronic kidney disease (CKD) · Minimal change nephrotic syndrome ·
FSGS · Membranous nephropathy · ADPKD · Benign prostatic hyperplasia (BPH) ·
Overactive bladder · Urolithiasis

(✓ = completed)
