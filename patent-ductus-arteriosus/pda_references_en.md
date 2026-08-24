# Patent Ductus Arteriosus of Prematurity — references

> The source literature for the QSP model `pda_qsp_model.dot` · `pda_mrgsolve_model.R` · `pda_reference_model.py`.
> Every PMID was resolved individually through the NCBI E-utilities and its title, authors, journal and year confirmed
> (nothing here was transcribed from memory). The numbers quoted in the text were checked against the
> abstracts fetched directly for each paper, and that check did in fact turn up two errors,
> which were corrected — see the footnote to §1 below.
>
> All PMIDs were resolved individually through NCBI E-utilities and the numbers
> quoted below were checked against the fetched abstracts. That check caught two
> errors of recollection, both corrected — see the footnote to §1.

---

## 1. The key clinical trials used to calibrate or falsify the model

| # | Reference | PMID | Role in the model |
|---|------|------|-----------------|
| 1 | **Semberova J, et al.** Spontaneous Closure of Patent Ductus Arteriosus in Infants ≤1500 g. *Pediatrics* 2017;140(2):e20164258 | [28701390](https://pubmed.ncbi.nlm.nih.gov/28701390/) | **FITTED.** Of 280 infants under genuinely conservative management, 237 (85%) closed spontaneously before discharge. Median time to closure **71 days (<26+0 weeks) · 13 days (26+0–27+6) · 8 days (28+0–29+6) · 6 days (≥30 weeks)**. These four values determine `TAUSYN0` · `KTAUSYN` · `KINVGA` · `NET50` · `NETW`. |
| 2 | **Gupta S, et al.** (Baby-OSCAR) Trial of Selective Early Treatment of Patent Ductus Arteriosus with Ibuprofen. *N Engl J Med* 2024;390(4):314–325 | [38265644](https://pubmed.ncbi.nlm.nih.gov/38265644/) | **FITTED.** 653 infants, 23–28 weeks, large PDA within 72 hours. Primary outcome (death or moderate-to-severe BPD) **69.2% (ibuprofen) versus 63.5% (placebo)**, adjusted RR 1.09 (95% CI 0.98–1.20). Death 13.6% versus 10.3%. The three risk coefficients `B_GA` · `B_BUR` · `B_VENT` were fitted to both arms simultaneously. |
| 3 | **Gupta S, et al.** Two-year outcomes after selective early treatment of PDA with ibuprofen. *EClinicalMedicine* 2025 | [40896454](https://pubmed.ncbi.nlm.nih.gov/40896454/) | Two-year follow-up of Baby-OSCAR. Confirms whether the absence of short-term improvement persists in the long term. |
| 4 | **Hundscheid T, et al.** (BeNeDuctus) Expectant Management or Early Ibuprofen for Patent Ductus Arteriosus. *N Engl J Med* 2023;388(11):980–990 | [36477458](https://pubmed.ncbi.nlm.nih.gov/36477458/) | **PREDICTED.** 273 infants, median gestation 26 weeks · birth weight 845 g. Composite (NEC/BPD/death) **46.3% (expectant management) versus 63.5% (early ibuprofen)**, risk difference −17.2%p, non-inferiority met. BPD 33.3% versus 50.9%. No parameter was fitted to this. |
| 5 | **Schmidt B, et al.** (TIPP) Long-term effects of indomethacin prophylaxis in extremely-low-birth-weight infants. *N Engl J Med* 2001;344(26):1966–1972 | [11430325](https://pubmed.ncbi.nlm.nih.gov/11430325/) | **PREDICTED.** 1202 infants. Prophylactic indomethacin reduced PDA to **24% versus 50%** and severe IVH to **9% versus 13%**, but death/disability at 18 months was **47% versus 46%** (OR 1.1). Why the model's IVH term is built from two terms of opposite sign (platelet COX-1 inhibition = harm, germinal matrix vascular maturation = benefit). |
| 6 | **Clyman RI, et al.** (PDA-TOLERATE) An Exploratory Randomized Controlled Trial of Treatment of Moderate-to-Large Patent Ductus Arteriosus at 1 Week of Age. *J Pediatr* 2019;205:41–48 | [30340932](https://pubmed.ncbi.nlm.nih.gov/30340932/) | **PREDICTED.** <28 weeks, early routine treatment versus conservative management at 6–14 days. PDA at ligation/discharge 32% versus 39%, NEC 16% versus 19%, BPD 49% versus 53%, BPD/death 58% versus 57%. Rescue criteria met 31% versus 62%. |
| 7 | **Mitra S, et al.** Association of Placebo, Indomethacin, Ibuprofen, and Acetaminophen With Closure of Hemodynamically Significant PDA: Systematic Review and Network Meta-analysis. *JAMA* 2018;319(12):1221–1238 | [29584842](https://pubmed.ncbi.nlm.nih.gov/29584842/) | **FITTED.** The ranking of closure rates by drug determines `KI_IBU` · `KI_IND` · `IC50_APAP`. High-dose oral ibuprofen ranks highest. |
| 8 | **Buvaneswarran S, et al.** Active Treatment vs Expectant Management of Patent Ductus Arteriosus in Preterm Infants: A Meta-Analysis. *JAMA Pediatr* 2025 | [40423988](https://pubmed.ncbi.nlm.nih.gov/40423988/) | The most recent pooled evidence — external validation of the model's central claim that "the duct closes but the outcome does not improve". |
| 9 | **Mitra S, et al.** Early treatment versus expectant management of hemodynamically significant PDA. *Cochrane Database Syst Rev* 2025 | [40548426](https://pubmed.ncbi.nlm.nih.gov/40548426/) | Cochrane pooling of early treatment versus expectant management. |
| 10 | **Mitra S, et al.** Interventions for patent ductus arteriosus (PDA) in preterm infants: an overview of Cochrane Systematic Reviews. *Cochrane Database Syst Rev* 2023 | [37039501](https://pubmed.ncbi.nlm.nih.gov/37039501/) | A summary of the whole intervention landscape. |
| 11 | **Mitra S, et al.** Prophylactic cyclo-oxygenase inhibitor drugs for the prevention of morbidity and mortality in extremely preterm infants. *Arch Dis Child Fetal Neonatal Ed* 2024 | [37419686](https://pubmed.ncbi.nlm.nih.gov/37419686/) | Pooled prophylactic COX inhibitors — the basis for scenario S7. |
| 12 | **Hundscheid T, et al.** Conservative Management of Patent Ductus Arteriosus in Preterm Infants — A Systematic Review and Meta-Analyses. *Front Pediatr* 2021 | [33718300](https://pubmed.ncbi.nlm.nih.gov/33718300/) | The evidence base for conservative management. |

> **Footnote — the two errors the check found.** (i) The Baby-OSCAR primary outcome had been remembered as
> 69.4%, whereas the published value is **69.2%**. (ii) More importantly,
> Semberova's median time to spontaneous closure had been assumed to follow a gentle curve against
> gestational age, with the 26-week target set at 48 days, whereas the paper's actual value at 26+0–27+6 weeks was
> **13 days**. The **5.5-fold discontinuity** between <26 weeks (71 days) and 26–27 weeks (13 days) occurs over about
> a week of gestational age, and it coincides exactly with where this model places its structural closure
> threshold (the point at which the residual lumen left by TMAXGA crosses the closure criterion, about 25.5 weeks).
> Discovering that a wrong target value had been in use has, in the end,
> given independent support for the model's structure.

---

## 2. Normal physiology of the ductus arteriosus — prostaglandin-dependent patency

| # | Reference | PMID | Content |
|---|------|------|------|
| 13 | **Coceani F, Olley PM.** Lamb ductus arteriosus: effect of prostaglandin synthesis inhibitors on the muscle tone and the response to prostaglandin E2. *Prostaglandins* 1975 | [1135442](https://pubmed.ncbi.nlm.nih.gov/1135442/) | The original discovery. Established that ductal patency depends on local prostaglandin synthesis. |
| 14 | **Segi E, et al.** Patent ductus arteriosus and neonatal death in prostaglandin receptor EP4-deficient mice. *Biochem Biophys Res Commun* 1998 | [9600059](https://pubmed.ncbi.nlm.nih.gov/9600059/) | PDA and neonatal death in EP4 gene-deficient mice — the basis for the model's `EP4` state variable being the obligatory mediator of ductal patency. |
| 15 | **Rheinlaender C, et al.** Changing expression of cyclooxygenases and prostaglandin receptor EP4 during development of the human ductus arteriosus. *Pediatr Res* 2006 | [16857763](https://pubmed.ncbi.nlm.nih.gov/16857763/) | COX and EP4 expression change with development in the human ductus — the direct basis for the `FSYN` involution of local ductal synthesis and for `EP4GA`·`KEP4`. |
| 16 | **Trivedi DB, et al.** Attenuated cyclooxygenase-2 expression contributes to patent ductus arteriosus in preterm mice. *Pediatr Res* 2006 | [17065565](https://pubmed.ncbi.nlm.nih.gov/17065565/) | Reduced COX-2 expression contributes to PDA in preterm mice — the basis for the gestational-age-dependent prostanoid term. |
| 17 | **Momma K, Takeuchi H.** Constriction of fetal ductus arteriosus by nonsteroidal antiinflammatory drugs. *Adv Prostaglandin Thromboxane Leukot Res* 1983 | [6221638](https://pubmed.ncbi.nlm.nih.gov/6221638/) | NSAID-induced constriction of the fetal ductus — the classic pharmacodynamic evidence. |
| 18 | **Printz MP, et al.** Studies of pulmonary prostaglandin biosynthetic and catabolic enzymes as factors in ductus arteriosus patency. *Pediatr Res* 1984 | [6422431](https://pubmed.ncbi.nlm.nih.gov/6422431/) | Pulmonary prostaglandin synthetic and catabolic enzymes as factors in ductal patency — the basis for the `PGDH` state variable (pulmonary 15-PGDH). |
| 19 | **Takizawa T, et al.** Inhibitory effect of indomethacin on neonatal lung catabolism of prostaglandin E2. *J Toxicol Sci* 1996 | [8959648](https://pubmed.ncbi.nlm.nih.gov/8959648/) | Neonatal pulmonary catabolism of PGE2 and its inhibition — why the clearance of circulating PGE2 is written as the product of pulmonary blood flow and 15-PGDH. |
| 20 | **Ovalı F.** Molecular and Mechanical Mechanisms Regulating Ductus Arteriosus Closure in Preterm Infants. *Front Pediatr* 2020 | [32984222](https://pubmed.ncbi.nlm.nih.gov/32984222/) | A comprehensive review of the molecular and mechanical mechanisms of ductal closure in the preterm infant. |

---

## 3. Oxygen sensing and the smooth-muscle contractile apparatus

| # | Reference | PMID | Content |
|---|------|------|------|
| 21 | **Archer SL, et al.** O2 sensing in the human ductus arteriosus: redox-sensitive K+ channels are regulated by mitochondria-derived hydrogen peroxide. *Biol Chem* 2004 | [15134333](https://pubmed.ncbi.nlm.nih.gov/15134333/) | Oxygen sensing in the human ductus — mitochondria-derived H2O2 regulates redox-sensitive K+ channels. The `O2SENSE → MEMPOT → CAL` path on the map. |
| 22 | **Bentley RET, et al.** The molecular mechanisms of oxygen-sensing in human ductus arteriosus smooth muscle cells: A comprehensive transcriptome profile. *Genomics* 2021 | [34245829](https://pubmed.ncbi.nlm.nih.gov/34245829/) | The oxygen-sensing transcriptome of human ductal smooth muscle — the mechanistic background to `P50_O2` · `P50_GA`. |
| 23 | **Kajimoto H, et al.** Oxygen activates the Rho/Rho-kinase pathway and induces RhoB and ROCK-1 expression in human and rabbit ductus arteriosus. *Circulation* 2007 | [17353442](https://pubmed.ncbi.nlm.nih.gov/17353442/) | Oxygen activates the Rho/ROCK pathway — the `ROCK` (calcium sensitisation) node on the map. |
| 24 | **Hong Z, et al.** Role of store-operated calcium channels and calcium sensitization in normoxic contraction of the ductus arteriosus. *Circulation* 2006 | [16982938](https://pubmed.ncbi.nlm.nih.gov/16982938/) | Store-operated calcium channels and calcium sensitisation — the `CASMC` · `MLC20` path. |

---

## 4. Constriction ≠ closure — ductal wall hypoxia and anatomical remodelling (the model's bistable switch)

| # | Reference | PMID | Content |
|---|------|------|------|
| 25 | **Kajino H, et al.** Vasa vasorum hypoperfusion is responsible for medial hypoxia and anatomic remodeling in the newborn lamb ductus arteriosus. *Pediatr Res* 2002 | [11809919](https://pubmed.ncbi.nlm.nih.gov/11809919/) | **The key evidence for the model's structure.** After constriction, **vasa vasorum hypoperfusion** produces medial hypoxia and anatomical remodelling. Why `WALLO2` is written as two terms — "thin wall = diffusion suffices / thick wall = vasa-dependent, cut off by constriction". |
| 26 | **Clyman RI, et al.** VEGF regulates remodeling during permanent anatomic closure of the ductus arteriosus. *Am J Physiol Regul Integr Comp Physiol* 2002 | [11742839](https://pubmed.ncbi.nlm.nih.gov/11742839/) | VEGF regulates remodelling during permanent closure — the `VEGFD → TGFBD → REMOD` chain. |
| 27 | **Chorne N, et al.** Risk factors for persistent ductus arteriosus patency during indomethacin treatment. *J Pediatr* 2007 | [18035143](https://pubmed.ncbi.nlm.nih.gov/18035143/) | Risk factors for persistent patency during indomethacin treatment — the clinical correlate of reopening and of non-response. |
| 28 | **Keller RL, Clyman RI.** Persistent Doppler flow predicts lack of response to multiple courses of indomethacin in premature infants with recurrent patent ductus arteriosus. *Pediatrics* 2003 | [12949288](https://pubmed.ncbi.nlm.nih.gov/12949288/) | Prediction of non-response to repeated treatment — the basis for the second-course scenario (S12). |
| 29 | **Louis D, et al.** Factors associated with non-response to second course indomethacin for PDA treatment in preterm neonates. *J Matern Fetal Neonatal Med* 2018 | [28391737](https://pubmed.ncbi.nlm.nih.gov/28391737/) | Factors in non-response to a second course. |

---

## 5. Pharmacokinetics — the three drugs in the preterm infant

| # | Reference | PMID | Content |
|---|------|------|------|
| 30 | **Aranda JV, Thomas R.** Systematic review: intravenous ibuprofen in preterm newborns. *Semin Perinatol* 2006 | [16813969](https://pubmed.ncbi.nlm.nih.gov/16813969/) | A systematic review of intravenous ibuprofen. |
| 31 | **Van Overmeire B, et al.** Ibuprofen pharmacokinetics in preterm infants with patent ductus arteriosus. *Clin Pharmacol Ther* 2001 | [11673749](https://pubmed.ncbi.nlm.nih.gov/11673749/) | Ibuprofen PK in preterm infants. |
| 32 | **Hirt D, et al.** An optimized ibuprofen dosing scheme for preterm neonates with patent ductus arteriosus, based on a population pharmacokinetic and pharmacodynamic study. *Br J Clin Pharmacol* 2008 | [18307541](https://pubmed.ncbi.nlm.nih.gov/18307541/) | An optimal dosing scheme based on popPK/PD — the basis for **clearance maturation with postnatal age** (`CLMAT_IBU`) and for the higher doses needed when treatment is late. |
| 33 | **Yaffe SJ, et al.** The disposition of indomethacin in preterm babies. *J Pediatr* 1980 | [7441407](https://pubmed.ncbi.nlm.nih.gov/7441407/) | Disposition of indomethacin in preterm infants — `CL_IND0` · `V1_IND`. |
| 34 | **Allegaert K, et al.** Intravenous paracetamol (propacetamol) pharmacokinetics in term and preterm neonates. *Eur J Clin Pharmacol* 2004 | [15071761](https://pubmed.ncbi.nlm.nih.gov/15071761/) | Intravenous paracetamol PK in neonates — `CL_APAP0` · `V1_APAP`. |
| 35 | **Padavia F, et al.** Population Pharmacokinetics of Intravenous Paracetamol and Its Metabolites in Extreme Preterm Neonates. *Clin Pharmacokinet* 2024 | [39578300](https://pubmed.ncbi.nlm.nih.gov/39578300/) | popPK of paracetamol and its metabolites in extremely preterm neonates — the glucuronide/sulphate and CYP oxidation fractions (`FNAPQI`). |
| 36 | **Pacifici GM.** Clinical pharmacology of ibuprofen and indomethacin in preterm infants with patent ductus arteriosus. *Curr Pediatr Rev* 2014 | [25088343](https://pubmed.ncbi.nlm.nih.gov/25088343/) | A comparison of the clinical pharmacology of the two NSAIDs. |
| 37 | **Thibaut C, et al.** Effect of ibuprofen on bilirubin-albumin binding during the treatment of patent ductus arteriosus in preterm infants. *J Matern Fetal Neonatal Med* 2011 | [21815876](https://pubmed.ncbi.nlm.nih.gov/21815876/) | Ibuprofen's displacement of bilirubin from albumin binding — the `BILIDISP_IBU` · `BFREE` state variables. |
| 38 | **Diot C, et al.** Effect of ibuprofen on bilirubin-albumin binding in vitro at concentrations observed during treatment of patent ductus arteriosus. *Early Hum Dev* 2010 | [20472375](https://pubmed.ncbi.nlm.nih.gov/20472375/) | The extent of in vitro displacement at therapeutic concentrations. |

---

## 6. Acetaminophen — why "a different site" (a testable prediction of the model)

| # | Reference | PMID | Content |
|---|------|------|------|
| 39 | **Ouellet M, Percival MD.** Mechanism of acetaminophen inhibition of cyclooxygenase isoforms. *Arch Biochem Biophys* 2001 | [11370851](https://pubmed.ncbi.nlm.nih.gov/11370851/) | **The basis for the model's prediction.** Acetaminophen acts not at the arachidonic acid channel but at the **peroxidase site**, and is therefore **competitive with peroxide**. The direct basis for the `IC50_APAP × (1 + PEROX/KPEROX)` form. |
| 40 | **Graham GG, Scott KF.** Mechanism of action of paracetamol. *Am J Ther* 2005 | [15662292](https://pubmed.ncbi.nlm.nih.gov/15662292/) | A mechanism review — the concept that potency falls in tissues with a high peroxide tone. |
| 41 | **Ayoub SS.** Paracetamol (acetaminophen): A familiar drug with an unexplained mechanism of action. *Temperature (Austin)* 2021 | [34901318](https://pubmed.ncbi.nlm.nih.gov/34901318/) | An up-to-date account of the mechanistic controversy — cited in order to make explicit that the model's peroxide term is a hypothesis. |
| 42 | **Hammerman C, et al.** Ductal closure with paracetamol: a surprising new approach to patent ductus arteriosus treatment. *Pediatrics* 2011 | [22065264](https://pubmed.ncbi.nlm.nih.gov/22065264/) | The first report of closing the duct with paracetamol. |
| 43 | **Oncel MY, et al.** Oral paracetamol versus oral ibuprofen in the management of patent ductus arteriosus in preterm infants: a randomized controlled trial. *J Pediatr* 2014 | [24359938](https://pubmed.ncbi.nlm.nih.gov/24359938/) | A randomised comparison — similar closure rates. |
| 44 | **Ohlsson A, Shah PS.** Paracetamol (acetaminophen) for patent ductus arteriosus in preterm or low-birth-weight infants. *Cochrane Database Syst Rev* 2015 | [25758061](https://pubmed.ncbi.nlm.nih.gov/25758061/) | Cochrane pooling. |
| 45 | **Pranata R, et al.** The efficacy and safety of oral paracetamol versus oral ibuprofen for patent ductus arteriosus closure in preterm neonates. *Indian Heart J* 2020 | [32768013](https://pubmed.ncbi.nlm.nih.gov/32768013/) | A meta-analysis — similar efficacy, better safety. |

---

## 7. Organ-specific toxicity — the same enzyme, different organs

| # | Reference | PMID | Content |
|---|------|------|------|
| 46 | **Ohlsson A, et al.** Ibuprofen for the treatment of patent ductus arteriosus in preterm or low birth weight (or both) infants. *Cochrane Database Syst Rev* 2020 | [32045960](https://pubmed.ncbi.nlm.nih.gov/32045960/) | Ibuprofen versus indomethacin: similar closure rates, less NEC and less transient renal impairment with ibuprofen. The basis for `KI_*_K`·`KI_*_G` and for the `VC_IND` versus `VC_IBU` contrast. |
| 47 | **Pacifici GM.** Differential renal adverse effects of ibuprofen and indomethacin in preterm infants: a review. *Clin Pharmacol* 2014 | [25114597](https://pubmed.ncbi.nlm.nih.gov/25114597/) | The difference in renal adverse effects between the two NSAIDs — the comparison evidence for the `GFR_PGE2` · `SCR` · `UO` outputs. |
| 48 | **Patel J, et al.** Randomized double-blind controlled trial comparing the effects of ibuprofen with indomethacin on cerebral hemodynamics in preterm infants with patent ductus arteriosus. *Pediatr Res* 2000 | [10625080](https://pubmed.ncbi.nlm.nih.gov/10625080/) | **The direct basis for `VC_IND` versus `VC_IBU`.** Indomethacin acutely depresses cerebral haemodynamics; ibuprofen does not. |
| 49 | **Mosca F, et al.** Comparative evaluation of the effects of indomethacin and ibuprofen on cerebral perfusion and oxygenation in preterm infants with patent ductus arteriosus. *J Pediatr* 1997 | [9386657](https://pubmed.ncbi.nlm.nih.gov/9386657/) | A comparison of cerebral perfusion and oxygenation. |
| 50 | **Pezzati M, et al.** Effects of indomethacin and ibuprofen on mesenteric and renal blood flow in preterm infants with patent ductus arteriosus. *J Pediatr* 1999 | [10586177](https://pubmed.ncbi.nlm.nih.gov/10586177/) | Mesenteric and renal blood flow — the basis for `QMESREL` · `QRENREL`. |
| 51 | **Watterberg KL, et al.** Prophylaxis of early adrenal insufficiency to prevent bronchopulmonary dysplasia: a multicenter trial. *Pediatrics* 2004 | [15574629](https://pubmed.ncbi.nlm.nih.gov/15574629/) | Increased spontaneous intestinal perforation when early hydrocortisone is combined with indomethacin — `B_SIP_HC` and scenario S13. |
| 52 | **Ment LR, et al.** Randomized low-dose indomethacin trial for prevention of intraventricular hemorrhage in very low birth weight neonates. *J Pediatr* 1988 | [3373405](https://pubmed.ncbi.nlm.nih.gov/3373405/) | IVH prevention with low-dose indomethacin — the basis for the `IVH_IND_PROT` (germinal matrix vascular maturation) term. |
| 53 | **Sallmon H, et al.** Mature and immature platelets during the first week after birth and incidence of patent ductus arteriosus. *Cardiol Young* 2020 | [32340633](https://pubmed.ncbi.nlm.nih.gov/32340633/) | Platelets and PDA — the `TXA2` · `BTIME` path. |

---

## 8. Haemodynamics — the size of the shunt and its consequences

| # | Reference | PMID | Content |
|---|------|------|------|
| 54 | **Lindner W, et al.** Stroke volume and left ventricular output in preterm infants with patent ductus arteriosus. *Pediatr Res* 1990 | [2320395](https://pubmed.ncbi.nlm.nih.gov/2320395/) | Left ventricular output in PDA — the `QMAX_LV` · `QP` scale. |
| 55 | **Evans N, Kluckow M.** Assessment of ductus arteriosus shunt in preterm infants supported by mechanical ventilation: effect of interatrial shunting. *J Pediatr* 1994 | [7965434](https://pubmed.ncbi.nlm.nih.gov/7965434/) | Methodology for assessing the shunt. |
| 56 | **Fajardo MF, et al.** Effect of positive end-expiratory pressure on ductal shunting and systemic blood flow in preterm infants with patent ductus arteriosus. *Neonatology* 2014 | [24193163](https://pubmed.ncbi.nlm.nih.gov/24193163/) | The effect of PEEP on the shunt and on systemic blood flow — the `PEEPADJ` node on the map. |
| 57 | **Paradisis M, et al.** Randomized trial of milrinone versus placebo for prevention of low systemic blood flow in very preterm infants. *J Pediatr* 2009 | [18822428](https://pubmed.ncbi.nlm.nih.gov/18822428/) | Low systemic blood flow — the clinical reality of `QSDEF` (systemic perfusion deficit). |
| 58 | **Zonnenberg I, de Waal K.** The definition of a haemodynamic significant duct in randomized controlled trials: a systematic literature review. *Acta Paediatr* 2012 | [21913976](https://pubmed.ncbi.nlm.nih.gov/21913976/) | A literature review showing that RCTs have defined a "haemodynamically significant duct" differently from one another — why the `QSIG` threshold is carried as a state variable. |
| 59 | **El-Khuffash A, et al.** A Patent Ductus Arteriosus Severity Score Predicts Chronic Lung Disease or Death before Discharge. *J Pediatr* 2015 | [26474706](https://pubmed.ncbi.nlm.nih.gov/26474706/) | A PDA severity score predicts outcome — the basis for `PDABUR` (burden) and for the targeted treatment strategy (S16). |
| 60 | **Groves AM, et al.** Introduction to neonatologist-performed echocardiography. *Pediatr Res* 2018 | [30072808](https://pubmed.ncbi.nlm.nih.gov/30072808/) | Neonatal echocardiography — the `DALPA` · `LAAO` · `ABSDIAST` measurement nodes on the map. |
| 61 | **Levy PT, et al.** Application of Neonatologist Performed Echocardiography in the Assessment and Management of Neonatal Heart Disease. *Pediatr Res* 2018 | [30072802](https://pubmed.ncbi.nlm.nih.gov/30072802/) | Clinical application. |
| 62 | **MacLellan A, et al.** Fluid restriction for treatment of symptomatic patent ductus arteriosus in preterm infants. *Cochrane Database Syst Rev* 2024 | [39692231](https://pubmed.ncbi.nlm.nih.gov/39692231/) | Fluid restriction — the `FLUIDR` node on the map. |

---

## 9. PDA burden and outcome — BPD · NEC · IVH

| # | Reference | PMID | Content |
|---|------|------|------|
| 63 | **Clyman RI, et al.** Patent ductus arteriosus, tracheal ventilation, and the risk of bronchopulmonary dysplasia. *Pediatr Res* 2022 | [33790415](https://pubmed.ncbi.nlm.nih.gov/33790415/) | **The basis for the `PDABUR` × `VENTIX` product structure.** BPD risk is not whether the duct is present but the interaction between duration of exposure and mechanical ventilation. |
| 64 | **Nawaytou H, et al.** Patent ductus arteriosus and the risk of bronchopulmonary dysplasia-associated pulmonary hypertension. *Pediatr Res* 2023 | [36804505](https://pubmed.ncbi.nlm.nih.gov/36804505/) | PDA and BPD-associated pulmonary hypertension. |
| 65 | **Hamrick SEG, et al.** Patent Ductus Arteriosus of the Preterm Infant. *Pediatrics* 2020 | [33093140](https://pubmed.ncbi.nlm.nih.gov/33093140/) | A comprehensive review — from pathophysiology to the management controversy. |
| 66 | **Benitz WE; AAP Committee on Fetus and Newborn.** Patent Ductus Arteriosus in Preterm Infants. *Pediatrics* 2016 | [26672023](https://pubmed.ncbi.nlm.nih.gov/26672023/) | The American Academy of Pediatrics clinical report — the standard citation for the position that "there is no evidence that treatment improves outcome". |
| 67 | **Härkin P, et al.** Morbidities associated with patent ductus arteriosus in preterm infants. Nationwide cohort study. *J Matern Fetal Neonatal Med* 2018 | [28651469](https://pubmed.ncbi.nlm.nih.gov/28651469/) | Comorbidities in a nationwide cohort — the scale of the `H_NEC0` · `H_IVH0` baseline risks. |

---

## 10. Invasive treatment — ligation and transcatheter closure

| # | Reference | PMID | Content |
|---|------|------|------|
| 68 | **Weisz DE, et al.** Association of Patent Ductus Arteriosus Ligation With Death or Neurodevelopmental Impairment Among Extremely Preterm Infants. *JAMA Pediatr* 2017 | [28264088](https://pubmed.ncbi.nlm.nih.gov/28264088/) | The association of ligation with death or neurodevelopmental impairment. |
| 69 | **Sathanandam S, et al.** Consensus Guidelines for the Prevention and Management of Periprocedural Complications of Transcatheter Patent Ductus Arteriosus Closure in Extremely Low Birth Weight Infants. *Pediatr Cardiol* 2021 | [34195869](https://pubmed.ncbi.nlm.nih.gov/34195869/) | Consensus guidelines for transcatheter closure. |
| 70 | **Baruteau AE, et al.** The Transcatheter Closure of Patent Ductus Arteriosus in Extremely Low-Birth-Weight Infants: Technique and Outcomes. *J Cardiovasc Dev Dis* 2023 | [38132644](https://pubmed.ncbi.nlm.nih.gov/38132644/) | Technique and outcomes. |
| 71 | **Morray BH, et al.** 3-year follow-up of a prospective, multicenter study of the Amplatzer Piccolo Occluder for transcatheter patent ductus arteriosus closure. *J Perinatol* 2023 | [37587183](https://pubmed.ncbi.nlm.nih.gov/37587183/) | Three-year follow-up of the Piccolo. |
| 72 | **Duboue PM, et al.** Post-ligation cardiac syndrome after surgical versus transcatheter closure of patent ductus arteriosus. *Eur J Pediatr* 2024 | [38381375](https://pubmed.ncbi.nlm.nih.gov/38381375/) | Post-ligation cardiac syndrome — the `POSTLIG` node on the map. |
| 73 | **Serrano RM, et al.** Comparison of 'post-patent ductus arteriosus ligation syndrome' in premature infants after surgical versus percutaneous closure. *J Perinatol* 2020 | [31578421](https://pubmed.ncbi.nlm.nih.gov/31578421/) | A comparison between the two methods. |
| 74 | **Mosalli R, Alfaleh K.** Prophylactic surgical ligation of patent ductus arteriosus for prevention of mortality and morbidity in extremely low birth weight infants. *Cochrane Database Syst Rev* 2008 | [18254095](https://pubmed.ncbi.nlm.nih.gov/18254095/) | Cochrane on prophylactic ligation. |
| 75 | **Liebowitz M, Clyman RI.** Prophylactic Indomethacin Compared with Delayed Conservative Management of the Patent Ductus Arteriosus in Extremely Preterm Infants. *J Pediatr* 2017 | [28396025](https://pubmed.ncbi.nlm.nih.gov/28396025/) | Prophylactic indomethacin versus delayed conservative management — scenario S7 versus S1. |

---

## 11. Modelling tools

- **mrgsolve** (R, ODE-based PK/PD): <https://mrgsolve.org/>
- **QSP in R with mrgsolve**: <https://vantage-research.net/qsp-in-r/>
- **gPKPDviz** — an mrgsolve-based Shiny tool for PK/PD simulation:
  paper <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · code <https://github.com/Genentech/gPKPDviz/>
- **Graphviz** (rendering the mechanistic map): <https://graphviz.org/>

---

## ⚠️ Disclaimer

This model and this reference compilation are for **educational and research purposes**. The parameters were
taken from the published literature or are approximations fitted to the trial results of §1 above, and they have not
been independently validated or certified. **They must not be used for actual clinical decisions, prescribing, or regulatory submission.**
In patent ductus arteriosus of prematurity in particular, "should it be treated at all" is itself an
unresolved question, and this model is a tool not for resolving that uncertainty but for
**describing it quantitatively**.
