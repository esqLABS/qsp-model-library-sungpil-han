# Chronic Granulomatous Disease (CGD) — QSP Model

> **In one line.** CGD is a disease of a single number — the electron flux φ through
> NOX2. This model **calibrates only the phagosome of the normal neutrophil and
> predicts the disease.** Two things come out of that. First, there is no threshold
> inside the phagosome: K(φ) is convex and is still accelerating even at φ=1. The
> threshold observed in the clinic — "20% correction is enough" — is made not by the
> phagosome but by the **organisms sheltered inside the phagocytes that failed to
> kill them (the sheltered compartment)**. Second, the **single mean number the
> DHR-123 assay reports is not a sufficient statistic** — a carrier (mosaic) and a
> hypomorph (uniform) with the same mean of 0.20 have different infection risk, and
> which of the two is worse **changes sign** at a DHR mean of 0.238.

---

## What this directory contains

| File | Contents |
|---|---|
| [`cgd_qsp_model.dot`](../../../chronic-granulomatous-disease/cgd_qsp_model.dot) · [`.svg`](../../../chronic-granulomatous-disease/cgd_qsp_model.svg) · [`.png`](../../../chronic-granulomatous-disease/cgd_qsp_model.png) | Mechanistic map. **152 nodes, 17 clusters, 228 edges** |
| [`cgd_kernel.py`](../../../chronic-granulomatous-disease/cgd_kernel.py) | **Layer 1** — the 16-state kernel of phagosomal oxidative chemistry. Calibrated on the normal neutrophil only |
| [`cgd_python_reference.py`](../../../chronic-granulomatous-disease/cgd_python_reference.py) | **Layers 2 and 3** — infectious focus + whole patient, **53 ODEs**. The implementation that was actually run |
| [`cgd_analysis.py`](../../../chronic-granulomatous-disease/cgd_analysis.py) | The complete calibration, verification and scenario runner (sections A–L) |
| [`cgd_reference_output.txt`](../../../chronic-granulomatous-disease/cgd_reference_output.txt) | The full output of that run. The source of every number in this README |
| [`cgd_calibration.json`](../../../chronic-granulomatous-disease/cgd_calibration.json) | The calibrated constants (machine-injected into the R file) |
| [`cgd_mrgsolve_model.R`](cgd_mrgsolve_model.R) | mrgsolve 53-ODE model (an equation-by-equation mirror of the Python implementation) |
| [`cgd_shiny_app.R`](../../../chronic-granulomatous-disease/cgd_shiny_app.R) | Shiny dashboard, **12 tabs** |
| [`sync_r_params.py`](../../../chronic-granulomatous-disease/sync_r_params.py) | Machine injection of the parameters from JSON into R (`--check` verifies against drift) |
| [`cgd_references.md`](cgd_references.md) | **95 papers**, all looked up live through the PubMed E-utilities |
| [`mkrefs.py`](mkrefs.py) | The reference generator |

To reproduce:

```bash
python3 cgd_analysis.py            # kernel calibration → exposure calibration → sections A-L (about 4 min)
python3 sync_r_params.py --check   # check that the R file is in sync with the JSON
dot -Tsvg cgd_qsp_model.dot -o cgd_qsp_model.svg
python3 mkrefs.py --refresh        # re-query the references
```

---

## What is calibrated and what is prediction

This model has one design principle: **spend the calibration only on the normal cell
and on the exposure epidemiology, and leave the disease, the drug effects and the
thresholds entirely as predictions.**

### The calibrated numbers — nine in all

**The kernel (five, all measurements on the normal neutrophil)**

| Parameter | Target | Source |
|---|---|---|
| `AH0` | 70% of phagosomal H₂O₂ becomes HOCl | Winterbourn 2006 |
| `alk_gain` | peak phagosomal pH 7.80 | Segal 1981 |
| `kK_comp` | peak phagosomal K⁺ 0.50 M | Reeves 2002 |
| `kox` | the normal neutrophil kills 95% of S. aureus in 60 min | standard bactericidal assay |
| `kprot` | a φ=0 neutrophil kills 40% in 60 min | standard bactericidal assay |

Only the fifth comes from a measurement in CGD, and what it buys is exactly one
thing: **the size of the oxidase-independent killing arm.** Without that number the
model would have no way of knowing that the CGD neutrophil is not simply an empty
bag.

**Whole body (four, the exposure epidemiology)** — λ_b, μN_b, λ_f, μN_f. Each pair
was solved against two targets: the infection rate of untreated X-CGD and the
infection rate of a normal individual. The case-fatality rates (1% for bacteria, 20%
for invasive aspergillosis) were **fixed** from the literature, so mortality is a
prediction and not a target.

### And the predictions — nothing was fitted to any of these

| | Model | Actual |
|---|---|---|
| HOCl molecules per bacterium ingested | **7.3 × 10⁷** | ~10⁸ (Winterbourn 2006, by an entirely different route) |
| Co-trimoxazole prophylaxis | 1.90 → **0.752**/patient-year (**60% reduction**) | 1.90 → 0.83 (56%) — Margolis 1990 |
| Itraconazole prophylaxis | **52% reduction** | 86% — Gallin 2003 |
| Interferon gamma-1b | **72% reduction** | 67% — ICGDCSG 1991 |
| Triple prophylaxis | **0.240** severe infections/patient-year | 0.2–0.4 (registries) |
| Untreated mortality | **0.057**/patient-year | ~0.05 (the pre-prophylaxis era) |
| Mortality on triple prophylaxis | **0.0107**/patient-year | 0.01–0.02 (modern cohorts) |
| Correction threshold | **11.9%** halves mortality, **44.0%** removes 90% of it | remission above 20% chimerism, relapse below 10% |

Interferon gamma is **forbidden to touch φ** — because the 1991 trial could not show
a reproducible recovery of superoxide. And yet a 72% reduction comes out. That means
the drug has to work only in the oxidase-independent arm, and if so its benefit ought
to be **flat** with respect to residual ROS — which is exactly what section G
produces. The clinical trial had no power to test this prediction.

---

## Five things this model says that are not in the textbooks

### 1. The catalase dogma collapses under the arithmetic

The textbook: "CGD patients are infected by **catalase-positive** organisms.
Catalase-negative organisms hand their own H₂O₂ over to the phagosome and so kill
themselves."

That is **two claims about two different columns**.

- **Claim 1 (catalase removes the H₂O₂)**: computed as pseudo-first-order rate
  constants inside the phagosome, MPO gives 2.3 × 10⁴ /s and S. aureus catalase
  708 /s — catalase's share is **2.9%**. To take even half of it, the organism would
  have to carry **32 times** the catalase it has. And in any case this claim is about
  the *normal* phagosome. In the CGD phagosome there is no H₂O₂ to remove in the
  first place. Delete catalase from the model and the killing moves from 0.9500 to
  0.9525 at normal φ=1, and from 0.4000 to 0.4000 at φ=0. Messina 2002 actually made
  that mutant and it was still fully virulent in CGD mice.
- **Claim 2 (the organism arms the phagosome itself)**: the pneumococcus, scaled to a
  1.2 fL phagosome, puts out H₂O₂ at 4.6 mM/s — **more than NOX2 itself
  (1.0 mM/s).** Not one organism in the table with a self-arming index above 0.5 is a
  CGD pathogen.

**And yet the list of organisms cannot distinguish the two hypotheses.** Scored as
classifiers, the two rules are **exactly equally accurate (11/14)**. Because within
this table catalase status and peroxigenicity are almost perfectly confounded. The
decades of textbooks that have cited the organism list as evidence for the catalase
mechanism have been citing evidence with no discriminating power. What discriminates
is the arithmetic, and the arithmetic tells against catalase.

> **Catalase-negativity is a marker; peroxigenicity is the mechanism.**
> The discriminating prediction: an organism that is catalase-*negative* and does
> *not* make hydrogen peroxide should behave like a CGD pathogen, and an organism
> that is catalase-*positive* and does make hydrogen peroxide should be exempt.
> The dogma predicts the opposite in both cases.

### 2. There is no threshold inside the phagosome — and yet there is one in the clinic

Reading K(φ) in steps of 0.2:

```
    phi 0.0 -> 0.2:  K gains 0.138 log      ← the cheapest 20%
    phi 0.2 -> 0.4:  K gains 0.193 log
    phi 0.4 -> 0.6:  K gains 0.216 log
    phi 0.6 -> 0.8:  K gains 0.246 log
    phi 0.8 -> 1.0:  K gains 0.286 log      ← the most expensive 20%
```

**The first 20% of oxidase activity is worth the least.** That is the exact opposite
of the intuition behind gene therapy and mixed chimerism — that 10–20% correction is
enough. So the clinical threshold has to be made outside the phagosome. And it is:
correcting one cell *adds* one killer and at the same time **removes one hiding
place**, and the term that is convex is the sheltering term.

There is a place where the phagosome does have a threshold — **pH**. Alkalinisation
is a race between the oxidase's proton consumption and the V-ATPase, and below about
57% of normal flux the V-ATPase simply wins. In terms of pH and the alkaline
proteases, the phagosome of a hypomorph with 30% residual oxidase is
indistinguishable from that of a null patient — despite the DHR assay separating the
two cleanly.

### 3. The DHR mean is not a sufficient statistic, and which of the two is favoured changes sign

Let f be the fraction of neutrophils with normal oxidase and r the residual activity
of the remaining cells. DHR-123 reports only the single number `f·1 + (1−f)·r`. The
phagosome sees the two populations separately.

```
     DHR mean  mosaic surv  uniform surv    ratio  which is better
         0.05       0.5737        0.5978    0.960  MOSAIC
         0.20       0.4949        0.5161    0.959  MOSAIC
         0.30       0.4423        0.3996    1.107  uniform
         0.50       0.3372        0.2003    1.684  uniform
```

**The ordering turns over at a DHR mean of 0.238.** The crossing point is the
inflection point of K(φ), and neither more nor less than that. The reason this
matters is that the clinic operates below it — carriers, mixed chimerism after RIC
transplantation, and gene-therapy marking fractions are all in the 5–25% band. At the
level of the focus the difference is larger still: at a DHR mean of 0.20 the mosaic
patient tolerates **1.21 times** the inoculum of the uniform patient (log₁₀ N_crit
4.06 against 3.84).

**So "20% of normal" on a flow report does not tell you which patient this is.** The
shape of the histogram — bimodal or uniformly shifted — has to be read with it.

### 4. And yet for the inflammatory phenotype the DHR mean **is** sufficient

Running the same two patients for a year:

```
        X-CGD carrier 20%      IL-1b   4.93  granuloma  0.600  colitis  3.40
        X-CGD hypomorph 20%    IL-1b   4.93  granuloma  0.600  colitis  3.40
```

Completely identical. The reason is structural: **efferocytosis and inflammasome
suppression average over the population, whereas killing happens one phagosome at a
time.** The same assay number is a sufficient statistic for one phenotype and not for
the other.

### 5. Granulomas form without organisms, and antibiotics cannot touch them

The result of running a patient inoculated with nothing at all for a year:

```
    arm                              IL-1b  IL-10  apoptotic  granuloma  colitis
    healthy control                   2.35   4.46        9.2      0.079     0.45
    X-CGD null, no therapy            7.21   4.46       36.7      0.998     5.58
    X-CGD, co-trimox+itra             7.21   4.46       36.7      0.998     5.58   ← no change
    X-CGD, prednisolone 1 mg/kg       6.46   4.46       36.7      0.236     1.09
    X-CGD, anakinra 2 mg/kg           7.21   4.46       36.7      0.039     0.17
    X-CGD, HSCT 95% chimerism         2.44   4.46        9.5      0.097     0.55
    X-CGD, gene therapy 20%           4.94   4.46       23.0      0.603     3.41
```

The mechanism is not a level but a **loop**. The CGD neutrophil dies **without having
externalised the oxidised phosphatidylserine** that marks it for clearance (PS
externalisation is itself an oxidative event). The corpse that is not cleared is at
once a stimulus for IL-1β and a vanished source of the IL-10 that ought to switch the
lesion off. That is why anakinra abolishes the colitis and co-trimoxazole cannot —
there were no organisms in it to begin with.

Look at the IL-10 column: in CGD it **does not fall but rises** — because there are
far more corpses. What fails is IL-10 *per corpse*, and the model shows that that
shortfall on its own is not enough.

---

## Model structure — three layers, narrow interfaces

```
 LAYER 1  cgd_kernel.py            phagosomal chemistry, 16 states
          ├─ NADPH oxidase electron flux φ · R_ox 2.0 mM/s
          ├─ O2•− → (spontaneous + MPO-catalysed dismutation) → H2O2 → Cpd I → HOCl
          ├─ MPO 5-state cycle (MPO / Cpd I / Cpd II / Cpd III), Cl− depletion
          ├─ phagosomal pH (V-ATPase against the oxidase), K+ charge compensation, granule matrix
          └─ the two non-oxidative killing arms — opposite pH dependences
                  ↓  only this one thing crosses over: K(φ)
 LAYER 2  focus_rhs()              infectious focus, 4 states (B, Bi, N, NEC)
          └─ output: N_crit — the deterministic boundary between clearance and disease
                  ↓  how often the exposure distribution crosses that boundary
 LAYER 3  patient_rhs()            patient, 53 ODEs
          ├─ granulopoiesis · neutrophil traffic · monocytes/macrophages
          ├─ bacterial focus + sheltered compartment + necrosis, Aspergillus conidia/hyphae
          ├─ cytokine network · failed efferocytosis · granuloma/colitis/fibrosis
          ├─ PK of 7 drugs (co-trimoxazole · itraconazole · voriconazole · IFN-γ ·
          │   prednisolone · anakinra, including active metabolites)
          └─ chimerism · gene-therapy marking fraction · survival hazard function
```

**Why N_crit is the centre of the model.** A deterministic ODE cannot produce "0.83
events per patient-year". It either clears or goes septic, one or the other. The
common workaround is to bolt a fitted hazard function onto the side of the model, but
then that function does all the work and nothing is learnt. This model instead
computes a **critical inoculum** from the focus equations and leaves the infection
rate to how often a **fixed, patient-independent exposure distribution** crosses that
boundary:

```
    infection rate = λ · P( log10 N0 > log10 N_crit )
```

---

## Defects found and fixed in the course of running it

Every one of these would have survived had the model never been **run**.

1. **Bacteria were resurrected out of nothing.** An ODE has no notion of "one
   organism", so the 1e-30 left behind by integrator round-off, multiplied at
   μ_b = 4/d for 180 days, becomes exp(720) and hence 1e10 CFU. The first draft
   reported a 1e10 CFU abscess in a CGD patient who had been inoculated with nothing.
   Fixed by putting an `X/(X+1)` factor on every proliferating compartment, so that
   below one organism there is nothing to divide.
2. **Interferon gamma increased infections.** The boost had been written as
   `s**(1/b)`, and with s<1 and 1/b<1 survival *rises*. The model reported that IFN-γ
   multiplied the CGD infection rate by 2.5. Fixed by adding (b−1)·K₀
   **additively** to the log kill — which is also the exact formalisation of the claim
   that it "acts only in the oxidase-independent arm".
3. **Phagocytosis was a net amplifier.** The return coefficient
   s_enc·k_release/(k_release−μ_bi) came to 1.44 > 1 at φ=0, so the CGD neutrophil
   gave back more than it took in. Every CGD patient then goes septic within days of
   birth. Fixed by lowering the intraphagosomal replication rate from 0.35 to
   0.10/d, which brings it to 0.75 < 1.
4. **Healthy people developed granulomas.** Granuloma and colitis had been written to
   respond to the *absolute value* of IL-1β, so a normal baseline of 2 pg/mL on its
   own gave GRAN 3.04. Fixed to respond to the **increment** above baseline.
5. **The phagosomal chemistry was dominated by bacterial SOD.** Giving bacterial SOD
   access to the luminal superoxide made it account for 99% of the O2•− sink.
   Superoxide is a charged species and cannot cross the bacterial outer membrane —
   bacterial SOD defends against its own *endogenous* superoxide. Left at a default of
   0 and retained only as a sensitivity-analysis item.
6. **MPO's superoxide cycle was not returning the H₂O₂.** The MPO→Cpd III→MPO pair is
   effectively the phagosome's SOD (2 O2•− → H₂O₂ + O₂), but the first draft omitted
   the production of H₂O₂, thereby starving MPO of the very substrate it needs, and
   underestimated HOCl by about 40-fold.
7. **The AH pool ran out and the competition for chlorination disappeared.** Written
   as a consumable, the one-electron donor was exhausted immediately, fixing the HOCl
   yield at 0.90 whatever AH0 was set to. The tyrosine residues of 100 g/L of
   phagosomal protein are at the mM level, so it was changed to a fixed pool.
8. **An intraphagosomal doubling time of 30 minutes distorted the balance of the
   arms.** With the broth doubling time, replication ate up a third of the normal
   neutrophil's killing, and the absurdity emerged that switching the non-oxidative
   arm off collapsed normal killing from 95% to 7%. For a nutrient-limited acidic
   phagosome, three hours is a defensible figure.
9. **Without necrosis there was no critical inoculum.** The fate of the focus became
   independent of the inoculum and every exposure met the same end. A focus that has
   outrun neutrophil recruitment sequesters itself from both the neutrophils and the
   antibiotics — the clinical fact that a CGD abscess has to be drained.
10. **The PK units were out by a factor of 1000 in three places.** A dose in mg was
    being divided by a volume in L and then multiplied by 1000 again. All of it was
    unified on the convention "amounts in mg, volumes in L, concentrations in mg/L"
    (IFN-γ alone in µg/L = ng/mL).

---

## Honest failures, and what was not tested

- **The peak superoxide concentration is 16-fold lower than Winterbourn's estimate.**
  The two models genuinely disagree about whether the principal sink for O2•− is MPO
  or spontaneous dismutation. Nothing downstream in this model reads [O2•−] —
  everything reads the HOCl flux — so the error is contained, but it is a real
  disagreement and not a rounding difference.
- **The itraconazole prediction misses badly: 52% against an actual 86%.** It is the
  only bad one of the four predictive tests, and the direction is suggestive — the
  model gave the azoles inhibition of hyphal elongation only, and no blockade of
  germination. Adding a germination-blockade arm would almost certainly make it fit,
  but then it would be a fit and not a prediction.
- **The Kuhns gradient is under-predicted twofold.** The actual difference in survival
  between the top and bottom quartiles is about fivefold; the model gives twofold.
  Most probably because the residual-ROS distribution of the virtual population is
  bunched in the low range where K(φ) is flattest, so that three of the four quartiles
  are effectively the same patient.
- **The position of the alkalinisation threshold (φ≈0.57)** depends on two barely
  measured numbers (the resting phagosomal pH and the V-ATPase rate constant). The
  *existence* of the threshold is structural but its *position* is not, and no
  clinical claim in this file depends on that position.
- **Segal against Jankowski/Grinstein** — whether the CGD phagosome really fails to
  alkalinise is unresolved in the literature and unresolved here too. The model
  follows Segal, but prices the alternative in section I (normal killing
  0.950 → 0.877, killing at φ=0 unchanged).
- **ω (each organism's dependence on oxidant) is an assumption.** It is the softest
  input in the whole model and it governs the entire organism spectrum. The
  self-arming index of section B does not depend on ω, but the Winkelstein ordering
  does.
- **What was not tested at all: age.** Every patient in this file is a 30 kg child
  with fixed pharmacokinetics. CGD presents from infancy into the sixties, and drugs,
  marrow and exposure all change across that range.
- **The autosomal recessive genotypes were not tested as separate entities.** The fact
  that p47phox deficiency is milder in every registry is expressed in this model only
  as "a higher φ_res". If any part of the AR/XL difference is not residual ROS, this
  model cannot see it.
- **There was no R toolchain.** `cgd_mrgsolve_model.R` and `cgd_shiny_app.R` mirror
  the Python implementation that was run, equation by equation, and their parameters
  are injected mechanically by `sync_r_params.py` so they cannot drift — but the R
  syntax itself has not been verified by an interpreter.

---

## References

[`cgd_references.md`](cgd_references.md) — **95 papers**, in 12 sections. Not one of the
titles, journals, years, authors or PMIDs was written from memory; all of them were
looked up live through the NCBI E-utilities. Each entry also records **what the model
took from that paper**, so that a search result at odds with the intention is exposed
rather than hidden.

The central references: Winkelstein 2000 (the 368-patient registry) · Kuhns 2010
(residual ROS and survival) · Winterbourn 2006 (phagosomal chemistry) · Segal 1981
(phagosomal pH) · Reeves 2002 (K⁺/proteases) · Messina 2002 (catalase-negative
S. aureus) · Pericone 2000 (pneumococcal H₂O₂) · Margolis 1990 · Gallin 2003 ·
ICGDCSG 1991 (the three prophylaxis trials) · de Luca 2014 (anakinra) · Güngör 2014
(RIC transplantation) · Kohn 2020 (lentiviral gene therapy).

---

## ⚠️ Disclaimer

This is a semi-quantitative QSP model for educational and research purposes. It has
not been independently validated or certified and must not be used for real clinical
decision-making, prescribing, or regulatory submission.
