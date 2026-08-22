# Anthracycline-Induced Cardiotoxicity QSP model — references

The list of literature used to build the mechanistic map of anthracycline cardiotoxicity (`aic_qsp_model.dot`), the mrgsolve ODE model
(`aic_mrgsolve_model.R`), the Shiny dashboard (`aic_shiny_app.R`) and the independent numerical
verification script (`aic_reference_check.py`). 62 references in all.

> **A note on PMID verification**: of the 62 references below, **61 are direct PMID links
> looked up through the PubMed E-utilities (`esearch`/`esummary`), with title, authors, journal and year
> checked against the actual record**. No identifier was entered by guesswork. For the remaining
> one (no. 51, Slamon 2001) the automated look-up could not settle on a single record, so rather than guess a
> PMID it has been left as a **search query URL** — because a wrong PMID would cite an entirely
> different paper.
>
> **Important**: every figure marked "model:" in `README.md` and in the model files is a
> **computed result** of this model and not a value from the literature. The literature was used only as the
> basis for parameter calibration and as anchor values, and that correspondence is
> set out in the CALIBRATION ANCHORS section at the top of `aic_mrgsolve_model.R`.

---

## 1. Epidemiology · cumulative dose-response · natural history

The basis for calibrating the model's analysis A1 (the cumulative dose-incidence curve) and the late progression components (KDNH, FIB).

1. Von Hoff DD, et al. **Risk factors for doxorubicin-induced congestive heart
   failure.** *Ann Intern Med.* 1979. — the original source of the cumulative dose-heart failure
   relationship; the primary anchor for the model's A1. <https://pubmed.ncbi.nlm.nih.gov/496103/>
2. Swain SM, Whaley FS, Ewer MS. **Congestive heart failure in patients treated
   with doxorubicin: a retrospective analysis of three trials.** *Cancer.* 2003.
   — heart failure about 5% at 630 mg/m²; a modern re-analysis of the dose-response.
   <https://pubmed.ncbi.nlm.nih.gov/12767102/>
3. Cardinale D, et al. **Early detection of anthracycline cardiotoxicity and
   improvement with heart failure therapy.** *Circulation.* 2015. — n=2625,
   CTRCD 9%, 98% of it within one year, median 3.5 months. The model's time-scale anchor.
   <https://pubmed.ncbi.nlm.nih.gov/25948538/>
4. Cardinale D, et al. **Anthracycline-induced cardiomyopathy: clinical
   relevance and response to pharmacologic therapy.** *J Am Coll Cardiol.*
   2010. — LVEF recovery rate by the time treatment was started (64% within 2 months → ~0% beyond 6 months);
   the anchor for the model's A5 (the window of reversibility). <https://pubmed.ncbi.nlm.nih.gov/20117401/>
5. Lipshultz SE, et al. **Chronic progressive cardiac dysfunction years after
   doxorubicin therapy for childhood acute lymphoblastic leukemia.**
   *J Clin Oncol.* 2005. — continuing progression in childhood survivors; the basis of the irreversible-pool
   plus growth-mismatch concept. <https://pubmed.ncbi.nlm.nih.gov/15837978/>
6. Mulrooney DA, et al. **Major cardiac events for adult survivors of childhood
   cancer diagnosed between 1970 and 1999.** *BMJ.* 2020. — late cardiac events in a
   long-term cohort. <https://pubmed.ncbi.nlm.nih.gov/31941657/>
7. Cardinale D, Iacopo F, Cipolla CM. **Cardiotoxicity of anthracyclines.**
   *Front Cardiovasc Med.* 2020. — comprehensive clinical review.
   <https://pubmed.ncbi.nlm.nih.gov/32258060/>
8. Camilli M, et al. **Anthracyclines, diastolic dysfunction and the road to
   heart failure in cancer survivors.** *Prog Cardiovasc Dis.* 2024. — the diastolic
   dysfunction · fibrosis pathway. <https://pubmed.ncbi.nlm.nih.gov/39025347/>

## 2. The Top2b-dependent DNA damage arm — the evidence for the peak-driven arm

The basis of the model's `TOXN → DSB → P53` arm and of the Top2b degradation mechanism of dexrazoxane.

9. Zhang S, et al. **Identification of the molecular basis of doxorubicin-
   induced cardiotoxicity.** *Nat Med.* 2012. — cardiomyocyte-specific Top2b-null mice are
   resistant to cardiotoxicity; the direct basis for making Top2b an independent state variable in this model.
   <https://pubmed.ncbi.nlm.nih.gov/23104132/>
10. Lyu YL, et al. **Topoisomerase IIβ mediated DNA double-strand breaks:
    implications in doxorubicin cardiotoxicity and prevention by dexrazoxane.**
    *Cancer Res.* 2007. <https://pubmed.ncbi.nlm.nih.gov/17875725/>
11. Vejpongsa P, Yeh ETH. **Topoisomerase 2β: a promising molecular target for
    primary prevention of anthracycline-induced cardiotoxicity.**
    *Clin Pharmacol Ther.* 2014. <https://pubmed.ncbi.nlm.nih.gov/24091715/>
12. Deng S, et al. **Dexrazoxane may prevent doxorubicin-induced DNA damage via
    depleting both topoisomerase II isoforms.** *BMC Cancer.* 2014. — the basis for
    implementing the dexrazoxane effect in the model as proteolysis of Top2b rather than as iron chelation.
    <https://pubmed.ncbi.nlm.nih.gov/25406834/>
13. Martin E, et al. **Evaluation of the topoisomerase II-inactive
    bisdioxopiperazine ICRF-161 as a protectant against doxorubicin-induced
    cardiomyopathy.** *Toxicology.* 2009. — analogues that chelate iron but do not
    inhibit Top2 protect only weakly → the experimental basis for separating the mechanisms.
    <https://pubmed.ncbi.nlm.nih.gov/19010377/>
14. Sawyer DB, et al. **Mechanisms of anthracycline cardiac injury: can we
    identify strategies for cardioprotection?** *Prog Cardiovasc Dis.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20728697/>

## 3. Iron · ROS · mitochondria · ferroptosis (the iron-redox arm) — the evidence for the AUC-driven arm

The basis of the model's `TOXR → ROS → MITOD → death` arm, of the LIP feedback and of the subcritical loop gain.

15. Ichikawa Y, et al. **Cardiotoxicity of doxorubicin is mediated through
    mitochondrial iron accumulation.** *J Clin Invest.* 2014. — mitochondrial
    iron accumulation is the key; the model's LIP pool and FEAMP term.
    <https://pubmed.ncbi.nlm.nih.gov/24382354/>
16. Wallace KB, Sardão VA, Oliveira PJ. **Mitochondrial determinants of
    doxorubicin-induced cardiomyopathy.** *Circ Res.* 2020.
    <https://pubmed.ncbi.nlm.nih.gov/32213135/>
17. Fang X, et al. **Ferroptosis as a target for protection against
    cardiomyopathy.** *Proc Natl Acad Sci USA.* 2019.
    <https://pubmed.ncbi.nlm.nih.gov/30692261/>
18. Tadokoro T, et al. **Mitochondria-dependent ferroptosis plays a pivotal
    role in doxorubicin cardiotoxicity.** *JCI Insight.* 2023.
    <https://pubmed.ncbi.nlm.nih.gov/36946465/>
19. Minotti G, et al. **Anthracyclines: molecular advances and pharmacologic
    developments in antitumor activity and cardiotoxicity.**
    *Pharmacol Rev.* 2004. — the standard review of quinone redox cycling · metabolites · iron chemistry.
    <https://pubmed.ncbi.nlm.nih.gov/15169927/>

## 4. Metabolism · doxorubicinol · pharmacogenomics

The basis for putting the metabolite in a separate myocardial pool (CHM, t½ 28 d) in the model, and FM variation (analysis A7).

20. Olson RD, et al. **Doxorubicin cardiotoxicity may be caused by its
    metabolite, doxorubicinol.** *Proc Natl Acad Sci USA.* 1988. — the original source of the
    metabolite hypothesis; the model's WM weighting and CHM accumulation.
    <https://pubmed.ncbi.nlm.nih.gov/2897122/>
21. Forrest GL, et al. **Human carbonyl reductase overexpression in the heart
    advances the development of doxorubicin-induced cardiotoxicity in
    transgenic mice.** *Cancer Res.* 2000. — CBR overexpression brings the toxicity forward
    → the causal basis of the FM parameter. <https://pubmed.ncbi.nlm.nih.gov/11016643/>
22. Blanco JG, et al. **Anthracycline-related cardiomyopathy after childhood
    cancer: role of polymorphisms in carbonyl reductase genes.**
    *J Clin Oncol.* 2012. <https://pubmed.ncbi.nlm.nih.gov/22124095/>
23. Aminkeng F, et al. **A coding variant in RARG confers susceptibility to
    anthracycline-induced cardiotoxicity in childhood cancer.**
    *Nat Genet.* 2015. <https://pubmed.ncbi.nlm.nih.gov/26237429/>
24. Visscher H, et al. **Pharmacogenomic prediction of anthracycline-induced
    cardiotoxicity in children.** *J Clin Oncol.* 2012. — transporter variants such as
    SLC28A3. <https://pubmed.ncbi.nlm.nih.gov/21900104/>

## 5. Pharmacokinetics

The basis of the parameters of the model's PK block (three compartments + metabolite + liposomal carrier).

25. Speth PAJ, et al. **Clinical pharmacokinetics of doxorubicin.**
    *Clin Pharmacokinet.* 1988. — the CL · Vss · terminal t½ anchor.
    <https://pubmed.ncbi.nlm.nih.gov/3042244/>
26. Callies S, et al. **A population pharmacokinetic model for doxorubicin and
    doxorubicinol in patients.** *Cancer Chemother Pharmacol.* 2003. — a simultaneous
    parent-metabolite model; the basis of FM · CLM · VM in this model.
    <https://pubmed.ncbi.nlm.nih.gov/12647011/>
27. Gabizon A, Shmeeda H, Barenholz Y. **Pharmacokinetics of pegylated
    liposomal doxorubicin: review of animal and human studies.**
    *Clin Pharmacokinet.* 2003. — the long circulating half-life and low free-drug peak of PLD;
    the model's KEL_LIP · KREL_LIP. <https://pubmed.ncbi.nlm.nih.gov/12739982/>

## 6. Formulation · dosing schedule — the clinical test of the peak hypothesis

The key basis of the model's A3 (comparison of exposure metrics). The observation that the outcome differs while the
cumulative dose is identical is the clinical counterpart of the claim that the injury arm reads the peak and not the AUC.

28. Legha SS, et al. **Reduction of doxorubicin cardiotoxicity by prolonged
    continuous intravenous infusion.** *Ann Intern Med.* 1982. — the original source for
    continuous infusion. <https://pubmed.ncbi.nlm.nih.gov/7059060/>
29. van Dalen EC, et al. **Different dosage schedules for reducing
    cardiotoxicity in people with cancer receiving anthracycline
    chemotherapy.** *Cochrane Database Syst Rev.* 2016. — heart failure on continuous
    infusion RR≈0.27; the anchor for the model's A3/A8.
    <https://pubmed.ncbi.nlm.nih.gov/26938118/>
30. O'Brien MER, et al. **Reduced cardiotoxicity and comparable efficacy in a
    phase III trial of pegylated liposomal doxorubicin versus conventional
    doxorubicin for first-line treatment of metastatic breast cancer.**
    *Ann Oncol.* 2004. <https://pubmed.ncbi.nlm.nih.gov/14998846/>
31. Batist G, et al. **Reduced cardiotoxicity and preserved antitumor efficacy
    of liposome-encapsulated doxorubicin and cyclophosphamide compared with
    conventional doxorubicin in metastatic breast cancer.** *J Clin Oncol.*
    2001. — antitumour effect preserved + cardiotoxicity reduced = clinical evidence for a separation of the therapeutic index.
    <https://pubmed.ncbi.nlm.nih.gov/11230490/>

## 7. Dexrazoxane

The basis of the model's A4 (dose-reduction equivalence) and of results 5 and 6.

32. de Baat EC, et al. **Dexrazoxane for preventing or reducing cardiotoxicity
    in adults and children with cancer receiving anthracyclines.**
    *Cochrane Database Syst Rev.* 2022. — heart failure RR≈0.29; a model anchor.
    <https://pubmed.ncbi.nlm.nih.gov/36162822/>
33. Swain SM, et al. **Delayed administration of dexrazoxane provides
    cardioprotection for patients with advanced breast cancer treated with
    doxorubicin-containing therapy.** *J Clin Oncol.* 1997.
    <https://pubmed.ncbi.nlm.nih.gov/9193324/>
34. Lipshultz SE, et al. **Assessment of dexrazoxane as a cardioprotectant in
    doxorubicin-treated children with high-risk acute lymphoblastic
    leukaemia: long-term follow-up of a prospective, randomised, multicentre
    trial.** *Lancet Oncol.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20850381/>

## 8. Biomarkers · strain imaging

The basis of the model's A2 (lead time) and of the TNI/BNP/GLS observation equations.

35. Cardinale D, et al. **Prognostic value of troponin I in cardiac risk
    stratification of cancer patients undergoing high-dose chemotherapy.**
    *Circulation.* 2004. — the predictive power of troponin, above all its high negative predictive value.
    <https://pubmed.ncbi.nlm.nih.gov/15148277/>
36. Cardinale D, et al. **Prevention of high-dose chemotherapy-induced
    cardiotoxicity in high-risk patients by enalapril.** *Circulation.* 2006.
    — a troponin-guided early intervention strategy. <https://pubmed.ncbi.nlm.nih.gov/17101852/>
37. Lipshultz SE, et al. **Predictive value of cardiac troponin T in pediatric
    patients at risk for myocardial injury.** *Circulation.* 1997.
    <https://pubmed.ncbi.nlm.nih.gov/9355905/>
38. Ky B, et al. **Early increases in multiple biomarkers predict subsequent
    cardiotoxicity in patients with breast cancer treated with doxorubicin,
    taxanes, and trastuzumab.** *J Am Coll Cardiol.* 2014.
    <https://pubmed.ncbi.nlm.nih.gov/24291281/>
39. Oikonomou EK, et al. **Assessment of prognostic value of left ventricular
    global longitudinal strain for early prediction of chemotherapy-induced
    cardiotoxicity: a systematic review and meta-analysis.**
    *JAMA Cardiol.* 2019. — the predictive power of a relative fall in GLS; the basis of the model's GLS weighting.
    <https://pubmed.ncbi.nlm.nih.gov/31433450/>
40. Negishi K, et al. **Independent and incremental value of deformation
    indices for prediction of trastuzumab-induced cardiotoxicity.**
    *J Am Soc Echocardiogr.* 2013. — the basis of the 15% relative GLS threshold.
    <https://pubmed.ncbi.nlm.nih.gov/23562088/>
41. Cardinale D, et al. **Trastuzumab-induced cardiotoxicity: clinical and
    prognostic implications of troponin I evaluation.** *J Clin Oncol.* 2010.
    <https://pubmed.ncbi.nlm.nih.gov/20679614/>

## 9. Cardioprotective drug trials

The basis of the model's A9 (upstream vs downstream) and of results 7 and 8. The model places statins on the blockade of ROS
generation (effective only during exposure) and ACEi/BB/ARNI on the reversible-deficit · fibrosis arm (effective after exposure as well),
and that distinction is there to explain the conflicting results of the trials below.

42. Neilan TG, et al. **Atorvastatin for anthracycline-associated cardiac
    dysfunction: the STOP-CA randomized clinical trial.** *JAMA.* 2023.
    — a ≥10% fall in LVEF, 9% vs 22%; the model's A9 anchor.
    <https://pubmed.ncbi.nlm.nih.gov/37552303/>
43. Acar Z, et al. **Efficiency of atorvastatin in the protection of
    anthracycline-induced cardiomyopathy.** *J Am Coll Cardiol.* 2011.
    <https://pubmed.ncbi.nlm.nih.gov/21851890/>
44. Bosch X, et al. **Enalapril and carvedilol for preventing chemotherapy-
    induced left ventricular systolic dysfunction (OVERCOME).**
    *J Am Coll Cardiol.* 2013. <https://pubmed.ncbi.nlm.nih.gov/23583763/>
45. Heck SL, et al. **Prevention of cardiac dysfunction during adjuvant breast
    cancer therapy (PRADA): extended follow-up of a 2×2 factorial,
    randomized, placebo-controlled, double-blind clinical trial of
    candesartan and metoprolol.** *Circulation.* 2021. — that the effect is small and its
    persistence limited is the basis of the model's conservative downstream weighting.
    <https://pubmed.ncbi.nlm.nih.gov/33993702/>
46. Avila MS, et al. **Carvedilol for prevention of chemotherapy-related
    cardiotoxicity: the CECCY trial.** *J Am Coll Cardiol.* 2018. — LVEF
    the primary endpoint was negative but troponin fell → the basis for splitting BB in the model into partial ROS
    scavenging + a neurohormonal arm. <https://pubmed.ncbi.nlm.nih.gov/29540327/>
47. Cardinale D, et al. **Anthracycline-induced cardiotoxicity: a multicenter
    randomised trial comparing two strategies for guiding prevention with
    enalapril (ICOS-ONE).** *Eur J Cancer.* 2018. — prophylactic dosing and troponin-guided
    dosing are equivalent. <https://pubmed.ncbi.nlm.nih.gov/29567630/>
48. Livi L, et al. **Cardioprotective strategy for patients with nonmetastatic
    breast cancer who are receiving an anthracycline-based chemotherapy: the
    SAFE randomized clinical trial.** *JAMA Oncol.* 2021.
    <https://pubmed.ncbi.nlm.nih.gov/34436523/>
49. Gregorietti V, et al. **Use of sacubitril/valsartan in patients with
    cardiotoxicity and heart failure due to chemotherapy.**
    *Cardiooncology.* 2020. — ARNI as salvage therapy; the model's largest recovery lever.
    <https://pubmed.ncbi.nlm.nih.gov/33292750/>
50. Quagliariello V, et al. **The SGLT-2 inhibitor empagliflozin improves
    myocardial strain, reduces cardiac fibrosis and pro-inflammatory
    cytokines in non-diabetic mice treated with doxorubicin.**
    *Cardiovasc Diabetol.* 2021. <https://pubmed.ncbi.nlm.nih.gov/34301253/>

## 10. Trastuzumab · ErbB2 (the HER2 blockade interaction)

The basis of the model's A6 (interaction) and of result 9. ErbB2 blockade is implemented as a "repair-capacity deficit", and
the traditional Type I/II dichotomy is reinterpreted as the ratio of the two deficits.

51. Slamon DJ, et al. **Use of chemotherapy plus a monoclonal antibody against
    HER2 for metastatic breast cancer that overexpresses HER2.**
    *N Engl J Med.* 2001.
    <https://pubmed.ncbi.nlm.nih.gov/?term=Slamon+chemotherapy+plus+monoclonal+antibody+against+HER2+metastatic+breast+cancer+overexpresses+HER2+2001>
52. Romond EH, et al. **Trastuzumab plus adjuvant chemotherapy for operable
    HER2-positive breast cancer.** *N Engl J Med.* 2005. — NSABP B-31/N9831;
    the heart failure rate on concurrent administration. <https://pubmed.ncbi.nlm.nih.gov/16236738/>
53. Perez EA, et al. **Trastuzumab plus adjuvant chemotherapy for human
    epidermal growth factor receptor 2-positive breast cancer: planned joint
    analysis of overall survival from NSABP B-31 and NCCTG N9831.**
    *J Clin Oncol.* 2014. <https://pubmed.ncbi.nlm.nih.gov/25332249/>
54. Crone SA, et al. **ErbB2 is essential in the prevention of dilated
    cardiomyopathy.** *Nat Med.* 2002. — the ErbB2 survival signal is indispensable; the model's
    WTRD·WTR_REP·WTR_REG. <https://pubmed.ncbi.nlm.nih.gov/11984589/>
55. Ewer MS, Lippman SM. **Type II chemotherapy-related cardiac dysfunction:
    time to recognize a new entity.** *J Clin Oncol.* 2005.
    <https://pubmed.ncbi.nlm.nih.gov/15860848/>
56. Telli ML, et al. **Trastuzumab-related cardiotoxicity: calling into
    question the concept of reversibility.** *J Clin Oncol.* 2007. — a rebuttal of the
    "complete reversibility" of Type II; the basis for the model treating reversible/irreversible as a ratio rather than a dichotomy.
    <https://pubmed.ncbi.nlm.nih.gov/17687157/>

## 11. Cardiomyocyte renewal — the quantitative basis of irreversibility

The direct basis for setting KREG = 2×10⁻⁵/day (≈0.7 %/yr) in the model.

57. Bergmann O, et al. **Evidence for cardiomyocyte renewal in humans.**
    *Science.* 2009. — ¹⁴C dating; annual renewal rate of adult cardiomyocytes ~1%.
    <https://pubmed.ncbi.nlm.nih.gov/19342590/>
58. Bergmann O, et al. **Dynamics of cell generation and turnover in the human
    heart.** *Cell.* 2015. <https://pubmed.ncbi.nlm.nih.gov/26073943/>

## 12. Clinical guidelines · definitions · surveillance strategy

The basis of the CTRCD definition (a 10-point fall in LVEF & <50%), the 15% relative GLS threshold and the surveillance interval.

59. Lyon AR, et al. **2022 ESC Guidelines on cardio-oncology developed in
    collaboration with the European Hematology Association (EHA), the
    European Society for Therapeutic Radiology and Oncology (ESTRO) and the
    International Cardio-Oncology Society (IC-OS).** *Eur Heart J.* 2022.
    <https://pubmed.ncbi.nlm.nih.gov/36017575/>
60. Plana JC, et al. **Expert consensus for multimodality imaging evaluation of
    adult patients during and after cancer therapy: a report from the American
    Society of Echocardiography and the European Association of Cardiovascular
    Imaging.** 2014. — the standard for the CTRCD definition.
    <https://pubmed.ncbi.nlm.nih.gov/25239940/>
61. Curigliano G, et al. **Management of cardiac disease in cancer patients
    throughout oncological treatment: ESMO consensus recommendations.**
    *Ann Oncol.* 2020. <https://pubmed.ncbi.nlm.nih.gov/31959335/>

## 13. Preclinical · in vitro models

The basis of the structure of individual differences (inter-individual variation in KDROS · TOXN50).

62. Burridge PW, et al. **Human induced pluripotent stem cell-derived
    cardiomyocytes recapitulate the predilection of breast cancer patients to
    doxorubicin-induced cardiotoxicity.** *Nat Med.* 2016. — patient-specific susceptibility is
    reproduced at the cellular level → the basis of the assumption that inter-individual variation in the model lies in
    injury susceptibility and not in pharmacokinetics.
    <https://pubmed.ncbi.nlm.nih.gov/27089514/>

---

## How the literature maps onto the structure of the model

| Model structure | Corresponding references |
|---|---|
| Top2b made an independent state variable (peak-driven, DSB → p53 memory) | 9, 10, 11, 14 |
| Dexrazoxane = proteasomal degradation of Top2b (not iron chelation) | 12, 13, 32 |
| A slow residual pool + metabolite accumulation (the AUC-driven redox arm) | 15, 16, 19, 20, 21 |
| The metabolite in a separate myocardial pool (t½ 28 d) · FM variation | 20, 21, 22, 26 |
| The two arms read different exposure metrics (schedule and formulation effects) | 28, 29, 30, 31 |
| Irreversibility of myocyte loss (KREG ≈ 0.7 %/yr) | 57, 58 |
| A reversible functional deficit + the window of reversibility (the fibrosis clock) | 4, 36, 49 |
| LVEF masked by compensatory hypertrophy → GLS and troponin lead | 35, 37, 38, 39, 40 |
| The upstream (statin) vs downstream (ACEi/BB/ARNI) asymmetry | 42, 43, 44, 45, 46, 47, 48, 49 |
| ErbB2 blockade = a repair-capacity deficit (a supra-additive interaction) | 51, 52, 53, 54 |
| Reinterpretation of the Type I/II dichotomy (the ratio of the two deficits) | 55, 56 |
| The CTRCD definition and surveillance | 59, 60, 61 |
| Inter-individual variation placed in injury susceptibility | 23, 24, 62 |

## Where this model is more optimistic than the literature

Left on the record for the sake of honesty. Both items should be treated as an upper bound.

- **Pegylated liposomal doxorubicin**: the model's relative risk of CTRCD is ≈ 0.08, whereas the clinical heart
  failure RR in the Cochrane meta-analysis is 0.20 (refs 29 · 30 · 31). Because the free-drug
  peak is the dominant driver of the nuclear injury arm in the model, the benefit of encapsulation is overestimated.
- **Protectants in combination (dexrazoxane + statin + ACEi/BB)**: the model brings CTRCD down to
  almost zero at 480 mg/m², but no randomised trial testing this combination exists.
