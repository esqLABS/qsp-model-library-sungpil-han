# Post-Traumatic Stress Disorder (PTSD) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking trauma exposure
> and risk/resilience factors (FKBP5 x childhood-adversity interaction,
> peritraumatic dissociation, social support) to noradrenergic/glucocorticoid
> fear-memory consolidation, HPA axis dysregulation (enhanced glucocorticoid
> negative feedback → paradoxically low/normal cortisol), locus
> coeruleus-noradrenergic hyperarousal, and amygdala-hippocampus-vmPFC
> fear-extinction circuit failure (amygdala hyperreactivity, vmPFC
> hypoactivation, impaired extinction recall) driving the four DSM-5
> symptom clusters (intrusion, avoidance, negative cognition/mood,
> hyperarousal) and a composite CAPS-5-like endpoint — coupled to SSRI
> (sertraline/paroxetine), prazosin (alpha-1 antagonist, nightmares),
> ketamine/esketamine (rapid NMDA-antagonist extinction facilitation), and
> MDMA-assisted / trauma-focused psychotherapy (session-dose driven
> extinction-learning boost) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`ptsd_qsp_model_en.dot`](ptsd_qsp_model_en.dot) |
| 🖼️ Map (SVG)             | [`ptsd_qsp_model_en.svg`](ptsd_qsp_model_en.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`ptsd_qsp_model_en.png`](ptsd_qsp_model_en.png) |
| ⚙️ mrgsolve ODE model     | [`ptsd_mrgsolve_model.R`](ptsd_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`ptsd_shiny_app_en.R`](ptsd_shiny_app_en.R) |
| 📚 References             | [`ptsd_references_en.md`](ptsd_references_en.md) |

---

## 1. Disease in one paragraph

Post-traumatic stress disorder (PTSD) is a chronic neuropsychiatric disorder that
arises after exposure to an event involving a threat to life or serious physical
harm (combat, sexual assault, motor-vehicle accident, natural disaster, etc.). It
begins when noradrenaline- and glucocorticoid-mediated memory-consolidation
mechanisms centred on the amygdala fix an excessively strong fear memory in place.
Its core pathophysiology is the coexistence of **paradoxical dysregulation of the
HPA axis** — enhanced glucocorticoid receptor (GR) negative feedback that keeps
baseline cortisol low or normal despite acute-phase stress (the Yehuda hypothesis)
— with **locus coeruleus-noradrenergic hyperarousal**. Neuroimaging studies
consistently show a circuit dysfunction in which **amygdala hyperreactivity**
combines with **ventromedial prefrontal cortex (vmPFC) hypoactivation** to cause
top-down fear inhibition (extinction) to fail, and **hippocampal volume reduction**
from repeated stress impairs context-based fear discrimination, producing an
overgeneralisation in which threat is felt even in safe situations. These circuit
abnormalities manifest as the four DSM-5 symptom clusters — intrusive
re-experiencing, avoidance, negative alterations in cognition/mood, and
hyperarousal/reactivity — accompanied by sleep disturbance (nightmares, REM
fragmentation) and comorbidities such as depression, substance use disorder, and
suicide risk. First-line pharmacotherapy is the FDA-approved SSRIs
(**sertraline, paroxetine**), with **prazosin** used adjunctively for nightmares.
More recently, the rapid extinction-learning-facilitating effect of
**ketamine/esketamine** and **MDMA-assisted psychotherapy** — whose efficacy has
been confirmed in phase 3 trials — have emerged as a new axis complementing
standard exposure-based psychotherapies (prolonged exposure, cognitive processing
therapy, EMDR).

## 2. Mechanistic map clusters (15 clusters, 127 nodes)

1. Risk factors/gene-environment interaction (trauma type, FKBP5 x childhood
   trauma, 5-HTTLPR, peritraumatic dissociation, social support, resilience)
2. Acute stress response and fear acquisition (BLA threat detection,
   noradrenaline/cortisol-mediated memory enhancement, acute stress disorder,
   overgeneralisation)
3. HPA axis dysregulation (CRH-ACTH-cortisol, enhanced GR feedback, paradoxical
   hypocortisolism, FKBP5 demethylation)
4. Locus coeruleus-noradrenergic hyperarousal (LC firing, NE release,
   alpha1/beta receptors, heart rate/HRV, startle reflex, REM sympathetic surges)
5. Amygdala-hippocampus-vmPFC circuit dysfunction (BLA hyperreactivity, vmPFC
   hypoactivation, hippocampal atrophy/suppressed neurogenesis, insular
   hyperactivation)
6. Fear-memory acquisition/extinction neuroplasticity (consolidation,
   reconsolidation, extinction learning, BDNF-TrkB/mTOR, context-dependent
   relapse)
7. Neuroinflammation and oxidative stress (IL-6, TNF-α, CRP, microglia, the
   kynurenine pathway)
8. DSM-5 symptom clusters B-E (intrusion, avoidance, negative cognition/mood,
   hyperarousal, dissociative subtype)
9. Sleep disturbance (REM fragmentation, nightmares, insomnia, reduced slow-wave
   sleep)
10. Comorbidities and functional outcomes (depression, substance use disorder,
    suicide risk, chronic pain, cardiovascular risk, functional impairment)
11. SSRI (sertraline/paroxetine) PK/PD
12. Prazosin (alpha-1 blocker) PK/PD
13. Ketamine/esketamine PK/PD
14. MDMA-assisted psychotherapy and trauma-focused psychotherapy (PE/CPT/EMDR) PK/PD
15. Clinical assessment/biomarkers/endpoints (CAPS-5, PCL-5, remission criteria,
    HRV, startle habituation, fMRI)

## 3. mrgsolve model (24 ODE compartments)

* **Drug PK (10 compartments)** — sertraline gut depot/plasma (2), paroxetine
  gut depot/plasma (2), prazosin gut depot/plasma (2), ketamine plasma/effect
  compartment (2, biophase delay), MDMA gut depot/plasma (2).
* **Disease/PD (13 compartments)** — cortisol (CORTISOL), locus
  coeruleus/NE tone (NE_TONE), amygdala reactivity index (AMYG_REACT), vmPFC
  inhibitory tone (VMPFC_TONE), fear-memory strength (FEAR_MEM), extinction-memory
  strength (EXT_MEM), cumulative psychotherapy dose (THERAPY_CUM), severity of the
  four DSM-5 clusters (INTRUSION/AVOIDANCE/NEGCOG/HYPERAROUSE), a sleep-disturbance
  index (SLEEP_DIST), and a composite CAPS-5 score (CAPS5).
* **Time tracking (1 compartment)** — FX_WEEKS.
* Core design: **SSRIs** directly lower the amygdala-reactivity target
  (amyg_target), while **prazosin** lowers only the locus-coeruleus target
  (ne_target) and the sleep-disturbance target (sleep_target), so it is kept
  separate as a treatment that acts selectively on nocturnal symptoms.
  **Ketamine** acts through an effect compartment (representing the biophase and
  the delay of mTOR/synaptogenesis) to transiently accelerate the rate of
  extinction-memory formation (ext_drive), while **MDMA** is implemented with a
  dual action: it directly blunts the amygdala target during the acute processing
  session while simultaneously accelerating extinction-memory formation.
  **Cumulative psychotherapy dose (THERAPY_CUM)** increases through non-drug
  events (session-dose events) and linearly accelerates the rate of
  extinction-memory formation (ext_boost), mechanistically representing the
  complementary combined effect of medication and psychotherapy.

### 10 scenarios

| # | Scenario | Calibration basis |
|---|---|---|
| 1 | Natural course - moderate trauma (SEVERITY=1.0), untreated | Kessler 1995 Arch Gen Psychiatry |
| 2 | Natural course - severe/repeated trauma + FKBP5 + dissociative subtype (SEVERITY=1.5) | Binder 2008 JAMA; Klengel 2013 Nat Neurosci |
| 3 | Sertraline 100-200 mg/day | Brady 2000 JAMA |
| 4 | Paroxetine 20-50 mg/day | Marshall 2001 Am J Psychiatry |
| 5 | Prazosin 1-15 mg qhs adjunct | Raskind 2013 Am J Psychiatry; Raskind 2018 NEJM |
| 6 | Trauma-focused psychotherapy, weekly x12 (PE/CPT/EMDR) | Foa 2005 JAMA; Resick 2002 JCCP |
| 7 | Ketamine 0.5 mg/kg IV x6 (2 weeks) | Feder 2014 JAMA Psychiatry; Feder 2021 Am J Psychiatry |
| 8 | MDMA-assisted psychotherapy (3 sessions + 12 preparation/integration sessions) | Mitchell 2021 Nat Med; Mitchell 2023 Nat Med |
| 9 | High resilience + mild trauma (spontaneous remission) | Southwick 2014 Eur J Psychotraumatol |
| 10 | Combination: SSRI + weekly psychotherapy + prazosin | Krystal 2017 Biol Psychiatry (combination rationale) |

## 4. Shiny Dashboard (8 tabs)

1. **Patient profile** — set trauma severity, FKBP5 risk, dissociative subtype,
   resilience, simulation duration.
2. **PK** — plasma concentrations by drug (sertraline/paroxetine/prazosin/
   ketamine/MDMA).
3. **PD key measures (fear circuit/HPA)** — amygdala/vmPFC/fear-extinction
   memory, cortisol/NE tone.
4. **Clinical endpoints** — CAPS-5 trajectory (including the remission
   threshold line), trends by symptom cluster, time to remission.
5. **Scenario comparison** — overlaid comparison across multiple scenarios and a
   summary table.
6. **Biomarkers** — cortisol, NE tone, amygdala reactivity, vmPFC tone,
   extinction memory, sleep-disturbance index.
7. **Sleep/nocturnal symptoms** — sleep-disturbance index trajectory (to confirm
   the prazosin effect).
8. **References** — the full reference list.

## 5. How to Run

```bash
# 1) Render the mechanistic map
dot -Tsvg ptsd_qsp_model_en.dot -o ptsd_qsp_model_en.svg
dot -Tpng -Gdpi=150 ptsd_qsp_model_en.dot -o ptsd_qsp_model_en.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("ptsd_mrgsolve_model.R") %>% param(TRAUMA_SEVERITY = 1.0)
e_sert <- ev(amt = 100, cmt = "SERT_GUT", time = 0, ii = 24, addl = 180)
e_pe   <- ev(amt = 1, cmt = "THERAPY_CUM", time = 168, ii = 168, addl = 11)
out <- mod %>% ev(e_sert) %>% ev(e_pe) %>% mrgsim(end = 24*7*52, delta = 24)  # 52-week follow-up
plot(out, c("AMYG_REACT", "VMPFC_TONE", "EXT_MEM", "CAPS5"))

# 3) Run the Shiny dashboard
shiny::runApp("ptsd_shiny_app_en.R")
```

## 6. Key Clinical Calibration Basis

| Endpoint | Comparator | Basis |
|---|---|---|
| Lifetime prevalence and chronic course | National epidemiological survey (NCS) | Kessler 1995 Arch Gen Psychiatry |
| Sertraline RCT response rate | Randomised controlled trial | Brady 2000 JAMA |
| Paroxetine fixed-dose RCT | Randomised controlled trial | Marshall 2001 Am J Psychiatry |
| Prazosin nightmare/sleep improvement (initial) vs. no effect (large confirmatory trial, heterogeneous response) | 2 RCTs | Raskind 2013 Am J Psychiatry; Raskind 2018 NEJM |
| Rapid symptom reduction with single/repeated ketamine dosing | 2 RCTs | Feder 2014 JAMA Psychiatry; Feder 2021 Am J Psychiatry |
| MDMA-assisted psychotherapy phase 3 efficacy (2 confirmatory trials) | Randomised controlled phase 3 | Mitchell 2021 Nat Med; Mitchell 2023 Nat Med |
| Effect size of prolonged exposure/cognitive processing therapy | RCTs and meta-analyses | Foa 2005 JAMA; Watts 2013 J Clin Psychiatry |
| FKBP5 x childhood-trauma risk interaction | Gene-environment association studies | Binder 2008 JAMA; Klengel 2013 Nat Neurosci |
| Amygdala hyperreactivity/vmPFC hypoactivation circuit | Neuroimaging meta-analysis | Shin & Liberzon 2010 Neuropsychopharmacology |
| Hippocampal volume reduction | Neuroimaging meta-analysis | Woon 2010 Prog Neuropsychopharmacol Biol Psychiatry |

## 7. Model Validation Status

This container has no R/mrgsolve runtime installed (`Rscript` is absent), so the
mrgsolve model was completed through **literature-based parameter design and
self-review of the code (confirming all 24 compartments map 1:1 to `$CMT`/`$ODE`/
`$INIT`, and reviewing the logic that deliberately separates the sites of action
of SSRI/prazosin/ketamine/MDMA into the amygdala, locus coeruleus, effect
compartment, and acute-session effect respectively)**, but its numbers were not
verified by actual compilation and integration. The `.dot` file was actually
rendered with Graphviz `dot` installed via `apt-get install graphviz`, and SVG/PNG
generation was confirmed (127 nodes, 15 clusters). References were individually
verified by cross-checking author/year/journal information against PubMed via
WebSearch to confirm PMIDs. It is recommended to confirm the numerical integration
results by running the "How to Run" steps above in an environment with an
mrgsolve/R installation.
