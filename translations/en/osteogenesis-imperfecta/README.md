# Osteogenesis Imperfecta (OI) — QSP Model

> Integrated Quantitative Systems Pharmacology model linking COL1A1/COL1A2
> (and recessive CRTAP/P3H1/PPIB/SERPINH1/FKBP10/WNT1) collagen-biosynthesis
> defects to osteoblast/osteocyte dysfunction, sclerostin-mediated Wnt
> suppression, RANKL/OPG-driven high bone turnover, and excess TGF-beta
> bioavailability — producing low bone mineral density and recurrent
> low-trauma fracture — coupled to bisphosphonate (pamidronate, zoledronic
> acid), denosumab (anti-RANKL), teriparatide (PTH1R anabolic, adult type I
> only), setrusumab (anti-sclerostin, investigational), and fresolimumab
> (anti-TGF-beta, investigational) PK/PD.

| Deliverable | File |
|---|---|
| 🗺️ Mechanistic map (DOT) | [`oi_qsp_model.dot`](oi_qsp_model.dot) |
| 🖼️ Map (SVG)             | [`oi_qsp_model.svg`](oi_qsp_model.svg) |
| 🖼️ Map (PNG, 150 dpi)    | [`oi_qsp_model.png`](oi_qsp_model.png) |
| ⚙️ mrgsolve ODE model     | [`oi_mrgsolve_model.R`](../../../osteogenesis-imperfecta/oi_mrgsolve_model.R) |
| 📊 Shiny dashboard        | [`oi_shiny_app.R`](oi_shiny_app.R) |
| 📚 References             | [`oi_references.md`](oi_references.md) |

---

## 1. Disease in One Paragraph

Osteogenesis imperfecta is the most common inherited bone-fragility disorder, caused by autosomal
dominant mutations (mostly glycine substitutions) in **COL1A1/COL1A2**, the genes encoding type I
collagen; it can also rarely arise from autosomal recessive mutations in the collagen
3-hydroxylation complex (CRTAP·P3H1·PPIB), chaperones (SERPINH1/HSP47, FKBP10), or WNT1. Quantitative
(null-allele) defects cause the relatively mild type I, while qualitative (structural) defects cause
a range of severities — from perinatally lethal type II, through type III with progressive
deformity, to moderate type IV (the Sillence classification). Misfolded collagen triggers ER stress
and the unfolded protein response (UPR), increasing osteoblast death, and the collapsed osteocyte
network oversecretes **sclerostin (SOST)**, suppressing Wnt/β-catenin signalling and further reducing
bone formation. At the same time, a **rising RANKL/OPG ratio** promotes osteoclast differentiation,
and **excess active TGF-beta** released from the disordered matrix impairs osteoblast maturation
while synergistically stimulating osteoclastogenesis. The result is a breakdown of the coupling
between bone formation and resorption, producing low bone mass with high-turnover bone metabolism,
recurrent low-trauma fractures, progressive long-bone bowing, vertebral compression fractures, and
extraskeletal manifestations including dentinogenesis imperfecta, blue-grey scleral discoloration,
hearing loss, and joint hypermobility. Standard treatment is intravenous bisphosphonates (cyclic
pamidronate, zoledronic acid), which induce osteoclast apoptosis to suppress bone resorption;
denosumab is tried in severe recessive forms, and teriparatide is used only in adult type I but is
contraindicated in children with open growth plates due to osteosarcoma risk. The anti-sclerostin
antibody setrusumab and the anti-TGF-beta antibody fresolimumab, targeting the disease mechanism
directly, are in clinical development.

## 2. Mechanistic Map Clusters (15 Clusters, 106 Nodes)

1. Genetic etiology (COL1A1/COL1A2 dominant · CRTAP/P3H1/PPIB 3-hydroxylation complex ·
   SERPINH1/FKBP10 chaperones · IFITM5/WNT1/SP7 recessive forms)
2. Collagen biosynthesis/secretion pathway (triple-helix formation, P4H/LH modification, HSP47, ER
   stress/UPR, procollagen cleavage, fibril assembly, LOX crosslinking)
3. Bone matrix/osteoblast dysfunction (matrix deposition defects, hypermineralisation, reduced bone
   strength, osteoblast apoptosis, osteocyte network collapse)
4. Bone remodelling coupling: WNT/sclerostin · RANKL-OPG (osteocyte SOST secretion, LRP5/6
   inhibition, rising RANKL/OPG ratio, osteoclast activation)
5. TGF-beta signalling hyperactivity (release from the latent matrix reservoir, SMAD2/3, impaired
   osteoblast maturation, increased osteoclastogenesis, muscle weakness)
6. Skeletal phenotype and the Sillence classification (types I-VIII, recurrent fractures, progressive
   bowing, vertebral compression fracture/scoliosis, short stature, basilar invagination)
7. Extraskeletal manifestations (dentinogenesis imperfecta, blue-grey sclerae, hearing loss, cardiac
   valve/aortic root dilation, joint hypermobility, easy bruising, restrictive lung disease)
8. Clinical endpoints (annual fracture rate, DXA BMD Z-score, bone turnover markers, growth
   Z-score, function/mobility, chronic pain/QoL, number of vertebral compression fractures)
9. Bisphosphonate PK/PD (pamidronate/zoledronic acid/oral agents, bone binding, FPPS inhibition,
   osteoclast apoptosis, acute-phase response, hypocalcaemia, ONJ/atypical fracture)
10. Denosumab PK/PD (anti-RANKL antibody, target-mediated drug disposition, rebound bone resorption
    on discontinuation)
11. Teriparatide/anabolic therapy (intermittent PTH1R stimulation, increased bone formation,
    contraindicated with open growth plates)
12. Anti-sclerostin antibody (setrusumab, the ASTEROID/ORBIT trials)
13. Surgical/rehabilitative interventions (intramedullary rodding, spinal fusion, physical therapy,
    vitamin D/calcium)
14. Integrated drug-disease feedback (bone remodelling unit activation frequency, coupling balance,
    bone mineral content trajectory)
15. (within cluster 5) The investigational anti-TGF-beta antibody fresolimumab

## 3. mrgsolve Model (22 ODE Compartments)

* **Drug PK (5 agents, 11 compartments)** — pamidronate (central/bone-bound, 2 compartments),
  zoledronic acid (central/bone-bound, 2 compartments), denosumab (SC depot/central, 2 compartments,
  including nonlinear Michaelis-Menten target-mediated elimination), teriparatide (SC depot/central,
  2 compartments), setrusumab (central/peripheral, 2-compartment linear mAb), fresolimumab (central,
  1-compartment linear mAb).
* **Disease/PD (11 compartments)** — sclerostin (SOST), RANKL, OPG, active TGF-beta (TGFB),
  osteoblast activity index (OB), osteoclast activity index (OC), bone mineral content index (BMC),
  bone formation marker P1NP, bone resorption marker CTX, cumulative fracture count (FX_CUM), height
  Z-score (HEIGHT_Z).
* Setrusumab lowers effective SOST signal via an Emax competitive-neutralisation model against
  sclerostin, and fresolimumab lowers effective TGF-beta the same way. Denosumab neutralises free
  RANKL via an Emax model, lowering the osteoclast-differentiation signal (the RANKL/OPG ratio). The
  amount of bisphosphonate bound to bone is converted into an osteoclast apoptosis probability
  (Emax, with a separate IC50 per drug — zoledronic acid is more potent than pamidronate at a lower
  dose). Fracture risk is integrated as an instantaneous hazard proportional to the BMC deficit
  relative to baseline, yielding the cumulative fracture count.

### 7 Scenarios

| # | Scenario | Calibration Basis |
|---|---|---|
| 1 | Natural history (untreated, type III/IV) | Marini 2017 Nat Rev Dis Primers, Rauch 2004 Lancet |
| 2 | Cyclic IV pamidronate (1 mg/kg/day x 3 days q3 months) | Glorieux 1998 NEJM |
| 3 | IV zoledronic acid (0.05 mg/kg every 6 months) | Vuorimies 2017, Barros 2012 |
| 4 | Denosumab SC (severe recessive form, off-label) | Hoyer-Kuhn 2014 JBMR |
| 5 | Teriparatide SC (adult type I only, 20 μg/day) | Orwoll 2014 JBMR |
| 6 | Setrusumab (anti-sclerostin, investigational) | Glorieux/ASTEROID · ORBIT trials |
| 7 | Fresolimumab (anti-TGF-beta, investigational) | Grafe 2014 Nat Med preclinical evidence |

## 4. Shiny Dashboard (8 Tabs)

1. **Patient profile** — Sillence classification severity, age, weight, growth-plate status.
2. **PK** — blood concentration by drug (pamidronate/zoledronic acid/denosumab/
   teriparatide/setrusumab/fresolimumab).
3. **Key PD metrics** — effective sclerostin/free RANKL, osteoblast/osteoclast activity index.
4. **Clinical endpoints** — BMC index (Z-score change), cumulative fracture count.
5. **Scenario comparison** — overlaid comparison of multiple scenarios plus a summary table.
6. **Biomarkers** — trends in the bone formation marker (P1NP) and bone resorption marker (CTX).
7. **Growth/skeleton** — height Z-score trajectory, teriparatide growth-plate contraindication
   warning.
8. **References** — the full reference list.

## 5. Usage

```bash
# 1) Render the mechanistic map
dot -Tsvg oi_qsp_model.dot -o oi_qsp_model.svg
dot -Tpng -Gdpi=150 oi_qsp_model.dot -o oi_qsp_model.png
```

```r
# 2) R/mrgsolve simulation
install.packages(c("mrgsolve","dplyr","tidyr","ggplot2","shiny","DT"))
library(mrgsolve)
mod <- mread("oi_mrgsolve_model.R")
e_pam <- ev(amt = 30, cmt = "PAM_CENT", ii = 24, addl = 2, rate = 30/2, time = 0) %>%
  ev_repeat(ii = 24*30*3, addl = 12)  # pamidronate, 3-month cycles x 4 years
out <- mod %>% ev(e_pam) %>% mrgsim(end = 24*365*4, delta = 24)
plot(out, c("BMC", "P1NP", "CTX", "FX_CUM"))

# 3) Launch the Shiny dashboard
shiny::runApp("oi_shiny_app.R")
```

## 6. Key Clinical Calibration Evidence

| Endpoint | Comparator | Evidence |
|---|---|---|
| Cyclic IV pamidronate, increased BMD | Improved lumbar spine BMD Z-score over 4 years | Glorieux 1998 NEJM |
| Overall bisphosphonate effect | Meta-analysis of fracture rate/BMD | Dwan 2016 Cochrane (CD005088) |
| Oral risedronate, no effect at 12 months | Fracture-prevention effect not demonstrated (PLUTO) | Bishop 2013 Lancet |
| Zoledronic acid vs. pamidronate non-inferiority | Twice-yearly vs. 3-monthly cycling | Vuorimies 2017 Horm Res Paediatr |
| Denosumab, severe recessive-form case series | Reduced fractures, rebound resorption on discontinuation | Hoyer-Kuhn 2014 JBMR |
| Teriparatide, adult type I | Increased BMD (only after growth-plate closure) | Orwoll 2014 JBMR |
| Anti-TGF-beta neutralisation, preclinical evidence | Normalised bone mass in oim/oim mice | Grafe 2014 Nat Med |

## 7. Model Validation Status

This container does not have an R/mrgsolve execution environment installed (no `Rscript`), so the
mrgsolve model was completed only through the stage of **literature-based parameter design and a
self-review of the code (dimension/limit checks)**; the numbers were not verified by actually
compiling and integrating the model. The `.dot` file was actually rendered into SVG/PNG and checked
using Graphviz `dot` (106 nodes, 15 clusters). It is recommended that anyone with an mrgsolve/R
environment run the model as described in "Usage" above to confirm the numerical integration
results.
