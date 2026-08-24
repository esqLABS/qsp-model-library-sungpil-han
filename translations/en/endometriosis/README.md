# Endometriosis — QSP Model

## Overview

Endometriosis is a chronic inflammatory estrogen-dependent disease affecting approximately **10% of reproductive-age women** (approximately 190 million women worldwide). It is characterised by the presence of endometrial-like tissue outside the uterine cavity, most commonly on the ovaries, fallopian tubes, peritoneum, and rectovaginal septum. The disease causes significant morbidity through chronic pelvic pain, dysmenorrhea, dyspareunia, and infertility. Diagnosis is typically delayed 7–10 years from symptom onset.

The pathophysiology involves a complex interplay of retrograde menstruation, impaired immune surveillance, local estrogen biosynthesis via aromatase upregulation, progesterone resistance, neuroangiogenesis, and a self-sustaining PGE2–aromatase positive feedback loop. The QSP model in this directory captures these mechanisms quantitatively using an ODE-based framework spanning the HPO axis, lesion dynamics, inflammation, and pain signalling.

---

## Mechanistic Map

| Pathway | Key Nodes | Clinical Outcome |
|---|---|---|
| HPO Axis | GnRH → LH/FSH → E2 → negative feedback | Regulation of the ovulatory cycle, E2 supply |
| Retrograde Menstruation | Endometrial cells → peritoneal cavity → implantation | Lesion formation |
| Aromatase Upregulation | SF-1 ↑, PGE2 → CYP19A1 | Local E2 overproduction |
| Progesterone Resistance | PR-B ↓, HSD17B2 ↓ | Failure to suppress endometrial proliferation |
| PGE2 Loop | COX-2 → PGE2 → EP2/EP4 → aromatase | Sustained inflammation, promotes proliferation |
| NK Dysfunction | Reduced NK cell number/function | Failure to clear ectopic cells |
| Neural Invasion | NGF → TrkA → pain signalling | Allodynia, central sensitisation |
| Oxidative Stress | ROS → NF-κB → IL-6/IL-8 | Amplifies inflammation, promotes cell survival |
| Angiogenesis | VEGF → new vessel formation | Nutrient supply to lesions |
| Macrophage Infiltration | M2 macrophages → TGF-β | Immune suppression, fibrosis |

---

## Drug PK/PD Parameters

| Drug | Mechanism | Key PK Parameters | Clinical Effect |
|---|---|---|---|
| **Leuprolide depot 3.75 mg/month** | GnRH agonist → desensitisation → pituitary suppression | t½ = 3 h (active form); depot release ~4 weeks; Tmax ~4 h | E2 < 20 pg/mL (castrate level); flare-up over 1–2 weeks |
| **Elagolix 150 mg/day** | GnRH antagonist → immediate FSH/LH suppression (partial) | t½ = 4–6 h; Tmax = 1 h; linear PK | Partial E2 suppression (~50%); menstruation may be maintained |
| **Elagolix 200 mg BID** | GnRH antagonist → complete suppression | t½ = 4–6 h; BID dosing | Complete E2 suppression (< 12 pg/mL); induces amenorrhoea |
| **Dienogest 2 mg/day** | Selective progestin → endometrial atrophy, (partial) E2 suppression | t½ = 9–10 h; Tmax = 1.5 h; F ≈ 91% | ~60% reduction in E2; lesion atrophy; reduced pain |
| **Letrozole 2.5 mg/day + Add-back** | Non-steroidal aromatase inhibitor | t½ = 48 h; Tmax = 1 h; hepatic metabolism | >90% reduction in peripheral E2; add-back therapy required to prevent bone loss |
| **Combined OCP (EE/DRSP)** | Suppresses the hypothalamic-pituitary axis; endometrial atrophy | EE: t½ = 6–12 h; DRSP: t½ = 30 h | Partial reduction in E2; pain relief; endometrial thinning |
| **NSAIDs (Naproxen 500 mg BID)** | COX-1/2 inhibition → reduced PGE2 | t½ = 14 h; Tmax = 2–4 h | Pain relief; limited effect on lesion proliferation |

---

## Pathophysiology Table

| Pathway | Key Mechanism | Clinical Abnormality |
|---|---|---|
| Dysregulated HPO axis | Altered GnRH pulsatility → excess FSH/LH → excess ovarian E2 production | Excessive endometrial proliferation, irregular menstruation |
| Retrograde menstruation | Viable endometrial cells in menstrual blood → implantation and proliferation in the peritoneal cavity | Formation of peritoneal/ovarian/rectal lesions |
| Aromatase upregulation | Ectopic SF-1 expression → CYP19A1 induction → local E2 biosynthesis | Incomplete response to systemic therapy |
| Progesterone resistance | Reduced PR-B expression, MAPK hyperactivation → diminished progestin effect | Reduced response to progestin therapy |
| PGE2–Aromatase positive feedback | PGE2 → EP2/EP4 → cAMP → SF-1 → aromatase → E2 → COX-2 → increased PGE2 | Self-amplifying inflammation/proliferation loop |
| NK cell dysfunction | Reduced NK cell number, reduced cytotoxicity → survival of ectopic cells | Failure of immune clearance → lesion persistence |
| Neural invasion / NGF overexpression | Ectopic endometrium → NGF → invasion by TrkA-expressing sensory nerves → central sensitisation | Chronic pain, allodynia, hyperalgesia |
| Oxidative stress | Excess iron in peritoneal fluid, ROS → NF-κB activation | Increased IL-6, IL-8, TNF-α → lesion growth |

---

## Treatment Scenarios

| Scenario | Drug | Goal | Predicted Effect |
|---|---|---|---|
| 1 | No Treatment | Observe baseline disease course | Gradual lesion increase; persistent pain |
| 2 | Leuprolide depot (GnRH agonist) | Complete E2 suppression, lesion regression | E2 < 20 pg/mL; 40–60% lesion reduction; risk of BMD loss |
| 3 | Elagolix 150 mg/day | Partial E2 suppression, pain control | ~50% reduction in E2; improved dysmenorrhoea; minimal BMD impact |
| 4 | Elagolix 200 mg BID | Complete E2 suppression, potent pain control | E2 < 12 pg/mL; amenorrhoea; caution for BMD loss (use under 2 years) |
| 5 | Dienogest 2 mg/day | Long-term maintenance, endometrial atrophy | Partial E2 suppression; lesion reduction; relative preservation of BMD |
| 6 | Letrozole 2.5 mg/day + Add-back | Peripheral aromatase blockade + bone protection | >90% reduction in E2; treats recurrent/post-surgical residual lesions |
| 7 | Combined oral contraceptive pill (OCP) | Suppresses the hormonal cycle, pain relief | Moderate reduction in E2; suitable for long-term use; for patients not desiring pregnancy |

---

## File List

| File | Description |
|---|---|
| `endo_shiny_app.R` | Shiny interactive dashboard (6 tabs) |
| `endo_references.md` | Reference list (44 entries, grouped by section) |
| `README.md` | This file — model overview and parameter summary |

> **Note**: The DOT/SVG/PNG mechanistic map files and the mrgsolve ODE model file are planned to be added later.

---

## Key References

1. **Giudice LC, Kao LC.** "Endometriosis." *Lancet* 2004;364(9447):1789–1799. [PMID: 15488215](https://pubmed.ncbi.nlm.nih.gov/15488215/)

2. **Zondervan KT, Becker CM, Missmer SA.** "Endometriosis." *N Engl J Med* 2020;382(13):1244–1256. [PMID: 31777792](https://pubmed.ncbi.nlm.nih.gov/31777792/)

3. **Taylor HS et al.** "Treatment of Endometriosis-Associated Pain with Elagolix, an Oral GnRH Antagonist." *N Engl J Med* 2017;377(1):28–40. [PMID: 28985706](https://pubmed.ncbi.nlm.nih.gov/28985706/)

4. **Strowitzki T et al.** "Dienogest is as effective as leuprolide acetate in treating the painful symptoms of endometriosis: a 24-week, randomized, multicentre, open-label trial." *Hum Reprod* 2010;25(3):633–641. [PMID: 20188814](https://pubmed.ncbi.nlm.nih.gov/20188814/)

5. **Bulun SE et al.** "Molecular basis for treating endometriosis with aromatase inhibitors." *Hum Reprod Update* 2000;6(5):413–418. [PMID: 12065431](https://pubmed.ncbi.nlm.nih.gov/12065431/)
