# Atrial Fibrillation (AF) — QSP Model

## Overview

Atrial fibrillation (AF) is the most common sustained arrhythmia, occurring in approximately 1–2% of the global adult population and reaching 5–10% in those aged 65 and over. As an independent risk factor for stroke, heart failure, cognitive impairment, and death, its annual healthcare cost exceeds $26 billion in the United States alone.

| Item | Content |
|------|------|
| **ICD-10** | I48 |
| **Prevalence** | 33 million worldwide (as of 2010) |
| **Stroke risk** | 5-fold increase compared with those without AF |
| **5-year mortality** | ~50% when heart failure is present |
| **Core treatment strategy** | Rhythm control vs rate control + anticoagulation therapy |

---

## Key Pathophysiological Pathways

### 1. Initiation Mechanism

| Mechanism | Key molecule/structure | Clinical outcome |
|------|--------------|---------|
| Pulmonary vein ectopic beats | PV sleeve muscle → spontaneous depolarization | AF-triggering focus (85-95%) |
| Calcium overload | RyR2 dysfunction → delayed afterdepolarization (DAD) | Triggered activity |
| Early afterdepolarization (EAD) | ICaL reactivation → upon IKr blockade | Occurs especially under bradycardic conditions |
| Autonomic trigger | Simultaneous sympathetic/vagal activation | Early-morning AF onset pattern |
| Stretch-activated channels | SAC (stretch-activated channels) | Increased atrial pressure → AF |

### 2. Electrical Remodeling — "AF begets AF"

| Change | Mechanism | Result |
|------|------|------|
| AERP shortening | ICaL downregulation (Cav1.2↓) + IK1 upregulation | Wavelength shortening → re-entry maintenance |
| Reduced conduction velocity | Abnormal Cx40/Cx43 distribution | Increased electrical discontinuity |
| Loss of rate-dependent adaptation | IKs downregulation | Reduced AERP rate-dependence |
| Ion channel redistribution | Constitutive IKACh activation | Restoration of normal rhythm becomes difficult |

### 3. Structural Remodeling — Fibrosis

| Signalling pathway | Mediator | Final change |
|----------|--------|---------|
| Ang II overactivation | TGF-β1 → Smad2/3 → fibroblast activation | Atrial interstitial fibrosis |
| Oxidative stress | NADPH oxidase → ROS | Mitochondrial dysfunction |
| NLRP3 inflammasome | Release of IL-1β, IL-18 | Chronic atrial inflammation |
| EAT-secreted factors | FGF-2, IL-6, TGF-β | Epicardial fat → fibrosis |
| Atrial enlargement | Increased LAP → mechanical stress | Left atrial diameter >45mm |

### 4. Thromboembolism Pathway — Virchow's Triad

| Virchow's Triad | Manifestation in AF | Clinical outcome |
|----------------|------------|---------|
| Stasis | Slow blood flow within the left atrial appendage (LAA) | LAA thrombus formation |
| Endothelial dysfunction | Increased atrial endothelial vWF, TF expression | Platelet activation |
| Hypercoagulable state | Increased fibrinogen, D-dimer, thrombin generation | Promotion of clot formation |

---

## Drug PK/PD Parameters

### Rhythm Control

| Drug | Class | Key PK | Key mechanism | Clinical effect |
|------|------|---------|---------|---------|
| **Amiodarone** | Class III | t½=40–55 days, Vd=60L/kg | INa/ICaL/IKr/IKs/β blockade | 65% 1-year sinus rhythm maintenance |
| **Dronedarone** | Class III* | t½=13–19h, F=15% (high-fat meal) | Multichannel blockade (non-iodinated) | ATHENA: 24%↓ CV death/hospitalisation |
| **Flecainide** | Class Ic | t½=12–27h, F=85–90% | Potent INa blockade (use-dependent) | Only in the absence of structural heart disease |
| **Propafenone** | Class Ic | t½=2–10h (CYP2D6-dependent) | INa/ICaL/IKr + β blockade | ~57% sinus rhythm conversion rate |
| **Sotalol** | Class III + β | t½=12h, renal excretion | IKr blockade + β-blocker | TdP risk monitoring required |
| **Dofetilide** | Class III | t½=10h, CrCl-dependent | Pure IKr blockade | QTc prolongation monitoring |
| **Vernakalant** | Atrial-selective | IV, t½=3h | Atrial-selective INa/IKur blockade | 50–52% sinus rhythm conversion rate |

### Rate Control

| Drug | Class | Key PK | Mechanism | Target heart rate |
|------|------|---------|------|------------|
| **Metoprolol** | β1-blocker | t½=3–7h, F=40%, CL=65L/h | β1 blockade → AV conduction suppression | <110 bpm at rest |
| **Diltiazem** | Non-DHP CCB | t½=3.5–5h, F=40% | ICaL → AV conduction suppression | <130 bpm on exertion |
| **Digoxin** | Cardiac glycoside | t½=36–48h, renal excretion | Na/K ATPase → vagal tone↑ | Preferred with concomitant chronic HF |
| **Ivabradine** | If blocker | t½=2h, CYP3A4 | If (HCN4) blockade → HR↓ | Not indicated for arrhythmia (AF) |

### Anticoagulation

| Drug | Target | Key PK | Stroke reduction | Bleeding risk |
|------|------|---------|------------|---------|
| **Warfarin** | Vitamin K Cycle (VKOR) | t½=40h, CYP2C9/2C19 | 64% vs placebo | Requires INR monitoring |
| **Apixaban** | Factor Xa | t½=12h, F=50%, CL=3.3L/h | ARISTOTLE: 21% vs warfarin | 58%↓ intracranial haemorrhage |
| **Rivaroxaban** | Factor Xa | t½=5–13h, F=80% (high-fat meal) | ROCKET AF: non-inferior | Slightly↑ GI bleeding |
| **Edoxaban** | Factor Xa | t½=10–14h, F=62% | ENGAGE-AF: non-inferior | 18%↓ bleeding vs warfarin |
| **Dabigatran** | Thrombin (IIa) | t½=12–17h, P-gp/renal | RE-LY: 34%↓ at 150mg | 69%↓ intracranial haemorrhage |

---

## Stroke Risk Score (CHA₂DS₂-VASc Score)

| Item | Score |
|------|------|
| Congestive heart failure (C) | 1 |
| Hypertension (H) | 1 |
| Age ≥75 years (A₂) | 2 |
| Diabetes mellitus (D) | 1 |
| Prior stroke/TIA (S₂) | 2 |
| Vascular disease (V) | 1 |
| Age 65–74 years (A) | 1 |
| Female sex (Sc) | 1 |

| Score | Annual stroke risk | Recommendation |
|------|---------------|------|
| 0 | ~0% | Anticoagulation not required |
| 1 | ~1.3% | Consider (in men) |
| ≥2 | ≥2.2% | Anticoagulation therapy recommended |

---

## Model Files

| File | Description |
|------|------|
| [af_qsp_model.dot](af_qsp_model.dot) | Graphviz mechanistic map (128+ nodes, 12 clusters) |
| [af_qsp_model.svg](af_qsp_model.svg) | SVG vector image (scalable) |
| [af_qsp_model.png](af_qsp_model.png) | PNG raster image (150 dpi) |
| [af_mrgsolve_model.R](af_mrgsolve_model.R) | mrgsolve ODE model (20+ compartments, 6 scenarios) |
| [af_shiny_app_en.R](af_shiny_app_en.R) | Shiny dashboard (6 tabs) |
| [af_references_en.md](af_references_en.md) | 46 references (including PubMed links) |

---

## Key Treatment Scenarios (mrgsolve Simulation)

| Scenario | Drug | Key outcome variable | Supporting trial |
|---------|------|------------|------------|
| 1. Untreated new-onset AF | — | Increasing AF burden, fibrosis progression, increased stroke risk↑ | Baseline model |
| 2. Rate control | Metoprolol 50mg BID | Heart rate <110 bpm, AF burden unchanged | RACE II |
| 3. Rhythm control (amiodarone) | Amiodarone 200mg/day | 65% sinus rhythm maintenance, AERP normalisation | AFFIRM |
| 4. Anticoagulation alone | Apixaban 5mg BID | 21%↓ stroke risk vs warfarin | ARISTOTLE |
| 5. Rate control + anticoagulation | Metoprolol + apixaban | Rate control + stroke prevention | RACE II + ARISTOTLE |
| 6. Rhythm control + anticoagulation | Amiodarone + apixaban | Sinus rhythm maintenance + maximal stroke reduction | EAST-AFNET 4 |

---

## Shiny App Structure (6 Tabs)

| Tab | Content |
|----|------|
| **Patient profile** | CHA₂DS₂-VASc calculator, HAS-BLED score, patient characteristic settings |
| **Pharmacokinetics (PK)** | Amiodarone/apixaban/metoprolol plasma concentration trends |
| **AF dynamics (PD)** | AF burden (%), sinus rhythm maintenance rate, AERP dynamics |
| **Anticoagulation & stroke** | anti-FXa/INR trends, stroke risk reduction, bleeding risk |
| **Scenario comparison** | Side-by-side comparison of 6 treatment strategies |
| **Biomarkers** | NT-proBNP, CRP, D-dimer, troponin trends |

---

## Key References Summary

- **AFFIRM** (2002): rhythm vs rate control — no difference in mortality. Emphasises the importance of anticoagulation therapy.
- **EAST-AFNET 4** (2020): early rhythm control — 21%↓ CV death/stroke/hospitalisation.
- **ARISTOTLE** (2011): apixaban vs warfarin — 21%↓ stroke, 58%↓ intracranial haemorrhage.
- **CHA₂DS₂-VASc** (2010): standardises anticoagulation therapy at a score of ≥2.

See [af_references_en.md](af_references_en.md) for the full reference list (46 articles, including PubMed links).
