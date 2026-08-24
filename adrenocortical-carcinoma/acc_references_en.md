# Adrenocortical carcinoma (ACC) QSP model — references
# Adrenocortical Carcinoma QSP Model — References

Every PMID was checked by direct lookup on PubMed at the time of writing (2026-08).
Links take the form `https://pubmed.ncbi.nlm.nih.gov/<PMID>/`.

> **Convention for citing the evidence behind model parameters.** The *"role in the model"* note on each
> entry below states which parameter or structure of `acc_mrgsolve_model_en.R` that reference supports.
> Couplings that were **ASSUMED** rather than measured directly in the literature are
> recorded as such. In particular, the **etoposide clearance induction coefficient** is not a value
> measured in a mitotane–etoposide interaction study but one derived from ① the magnitude of CYP3A4 induction (midazolam probe) and
> ② the fraction of etoposide metabolised by CYP3A4 (roughly 0.3–0.4 in the literature).

---

## 1. Pivotal clinical trials

1. **Fassnacht M, Terzolo M, Allolio B, Baudin E, Haak H, Berruti A, et al.; FIRM-ACT Study Group.**
   Combination chemotherapy in advanced adrenocortical carcinoma.
   *N Engl J Med.* 2012;366(23):2189-97. PMID: **22551107**
   <https://pubmed.ncbi.nlm.nih.gov/22551107/>
   *Role in the model:* The calibration standard for the two scenarios EDP-M vs Sz-M (07, 10).
   ORR 23.2% vs 9.2%, PFS 5.0 vs 2.1 months, OS 14.8 vs 12.0 months (no significant difference).
   The model's `SLPETO`/`SLPDOX`/`SLPCIS`/`SLPSZ` are set to reproduce this ordering.

2. **Terzolo M, Angeli A, Fassnacht M, Daffara F, Tauchmanova L, Conton PA, et al.**
   Adjuvant mitotane treatment for adrenocortical carcinoma.
   *N Engl J Med.* 2007;356(23):2372-80. PMID: **17554118**
   <https://pubmed.ncbi.nlm.nih.gov/17554118/>
   *Role in the model:* The basis for the adjuvant setting. The RFS endpoint.

3. **Terzolo M, Fassnacht M, Perotti P, Libe R, Kastelan D, Lacroix A, et al.**
   Adjuvant mitotane versus surveillance in low-grade, localised adrenocortical
   carcinoma (ADIUVO): an international, randomised, phase 3 trial.
   *Lancet Diabetes Endocrinol.* 2023;11(10):720-730. PMID: **37619579**
   <https://pubmed.ncbi.nlm.nih.gov/37619579/>
   *Role in the model:* The finding that adjuvant mitotane gives no benefit in the low-risk group —
   the basis for the structure in which adjuvant benefit varies with `KI67`/`ENSATst` in the model.

4. **Fassnacht M, Berruti A, Baudin E, Demeure MJ, Gilbert J, Haak H, et al.**
   Linsitinib (OSI-906) versus placebo for patients with locally advanced or
   metastatic adrenocortical carcinoma: a double-blind, randomised, phase 3 study.
   *Lancet Oncol.* 2015;16(4):426-35. PMID: **25795408**
   <https://pubmed.ncbi.nlm.nih.gov/25795408/>
   *Role in the model:* Scenario 17. The fact that IGF1R inhibition failed even though IGF2 is
   overexpressed in ~90% of cases — the model explains this as a combination of ① a small `APROLIF`
   (ACC growth is not limited by IGF) and ② the insulin–IR-A rescue pathway.

5. **Jones RL, Kim ES, Nava-Parada P, Alam S, Johnson FM, Stephens AW, et al.**
   Phase I study of intermittent oral dosing of the insulin-like growth factor-1
   and insulin receptors inhibitor OSI-906 in patients with advanced solid tumors.
   *Clin Cancer Res.* 2015;21(4):693-700. PMID: **25208878**
   <https://pubmed.ncbi.nlm.nih.gov/25208878/>
   *Role in the model:* The basis for linsitinib PK (`CLLIN`, `V1LIN`, `FLIN`) and for dose-limiting
   hyperglycaemia.

6. **Habra MA, Stephen B, Campbell M, Hess K, Tapia C, Xu M, et al.**
   Phase II clinical trial of pembrolizumab efficacy and safety in advanced
   adrenocortical carcinoma.
   *J Immunother Cancer.* 2019;7(1):253. PMID: **31533818**
   <https://pubmed.ncbi.nlm.nih.gov/31533818/>
   *Role in the model:* Scenarios 14–16. `SLPIMM`/`EMAXPD1` are set to reproduce an ORR of
   ~14–23%.

7. **Raj N, Zheng Y, Kelly V, Katz SS, Chou J, Do RKG, et al.**
   PD-1 Blockade in Advanced Adrenocortical Carcinoma.
   *J Clin Oncol.* 2020;38(1):71-80. PMID: **31644329**
   <https://pubmed.ncbi.nlm.nih.gov/31644329/>
   *Role in the model:* As above. The clinical observation that response is lower in patients with
   steroid excess is the starting point for result F.

8. **Berruti A, Terzolo M, Pia A, Angeli A, Dogliotti L.**
   Mitotane associated with etoposide, doxorubicin, and cisplatin in the
   treatment of advanced adrenocortical carcinoma.
   *Cancer.* 1998;83(10):2194-200. PMID: **9827725**
   <https://pubmed.ncbi.nlm.nih.gov/9827725/>
   *Role in the model:* The basis for the EDP-M regimen (etoposide 100 mg/m² d2-4, doxorubicin
   40 mg/m² d1, cisplatin 40 mg/m² d3-4, q28d).

9. **Khan TS, Imam H, Juhlin C, Skogseid B, Grondal S, Tibblin S, et al.**
   Streptozocin and o,p'DDD in the treatment of adrenocortical cancer patients:
   long-term survival in its adjuvant use.
   *Ann Oncol.* 2000;11(10):1281-7. PMID: **11106117**
   <https://pubmed.ncbi.nlm.nih.gov/11106117/>
   *Role in the model:* The basis for the Sz-M regimen (1 g d1-5 induction, then 2 g q3wk).

10. **Megerle F, Kroiss M, Hahner S, Fassnacht M, et al.**
    Mitotane Monotherapy in Patients With Advanced Adrenocortical Carcinoma.
    *J Clin Endocrinol Metab.* 2018;103(4):1686-1695. PMID: **29452402**
    <https://pubmed.ncbi.nlm.nih.gov/29452402/>
    *Role in the model:* Scenario 05. The benefit of monotherapy at low tumour burden —
    the basis for the model's `SLPMIT` (mitotane's own cytotoxicity is weak).

---

## 2. Mitotane pharmacokinetics · therapeutic drug monitoring — the heart of the model

11. **Arshad U, Taubert M, Kurlbaum M, Frechen S, Herterich S, Megerle F, et al.**
    Enzyme autoinduction by mitotane supported by population pharmacokinetic
    modelling in a large cohort of adrenocortical carcinoma patients.
    *Eur J Endocrinol.* 2018;179(5):287-297. PMID: **30087117**
    <https://pubmed.ncbi.nlm.nih.gov/30087117/>
    *Role in the model:* **The central evidence for this model.** ① Apparent volume of distribution ~6,086 L
    (interindividual variability **CV 81.5%**) → reproduced by `V1MIT + VLEANTIS + KPFAT*FATKG` = 5,966 L
    in the model. ② **BMI affects mitotane disposition** → the basis for making fat mass the
    determining variable of the depot volume. ③ **Linear enzyme autoinduction** of clearance →
    `FMIND`. ④ The high-dose start regimen and the recommendation of a **first TDM on day 16** → scenarios 01–04.

12. **Terzolo M, Baudin AE, Ardito A, Kroiss M, Leboulleux S, Daffara F, et al.**
    Mitotane levels predict the outcome of patients with adrenocortical
    carcinoma treated adjuvantly following radical resection.
    *Eur J Endocrinol.* 2013;169(3):263-70. PMID: **23704714**
    <https://pubmed.ncbi.nlm.nih.gov/23704714/>
    *Role in the model:* The therapeutic target concentration **≥14 mg/L** and the exposure–response relationship.
    The baseline for the model's `INWIN`/`TIW` outputs and for results A·B.

13. **Kerkhofs TM, Baudin E, Terzolo M, Allolio B, Chadarevian R, Mueller HH, et al.**
    Comparison of two mitotane starting dose regimens in patients with advanced
    adrenocortical carcinoma.
    *J Clin Endocrinol Metab.* 2013;98(12):4759-67. PMID: **24057287**
    <https://pubmed.ncbi.nlm.nih.gov/24057287/>
    *Role in the model:* **The core of result B.** The RCT finding that a high-dose start showed no
    significant difference in either mitotane concentration or adverse event rate. The model explains this as
    "dose was randomised but body composition was not stratified, plus the reported endpoint (time to
    target) censors exactly the patients the covariate harms most", and shows that this
    does not contradict the high-dose recommendation of the popPK study (no. 11).

14. **Kerkhofs TM, Derijks LJ, Ettaieb MH, Eekhoff EM, Neef C, Gelderblom H, et al.**
    Short-term variation in plasma mitotane levels confirms the importance of
    trough level monitoring.
    *Eur J Endocrinol.* 2014;171(6):677-83. PMID: **25201518**
    <https://pubmed.ncbi.nlm.nih.gov/25201518/>
    *Role in the model:* The observation that diurnal variation is limited → the basis for the structural
    decision to make the central compartment large (`V1MIT` 400 L), covering blood plus rapidly
    equilibrating lean tissue, which suppresses the diurnal amplitude.

15. **Puglisi S, Calabrese A, Basile V, Ceccato F, Scaroni C, Altieri B, et al.**
    Mitotane Concentrations Influence Outcome in Patients with Advanced
    Adrenocortical Carcinoma.
    *Cancers (Basel).* 2020;12(3):740. PMID: **32245135**
    <https://pubmed.ncbi.nlm.nih.gov/32245135/>
    *Role in the model:* The concentration–outcome relationship holds in advanced disease too.

16. **Puglisi S, Perotti P, Cosentini D, Fiorentini C, Basile V, et al.**
    Decreased Plasma Mitotane Levels Are Associated with Improved
    Recurrence-Free Survival in Adrenocortical Carcinoma Patients Treated
    Adjuvantly — target attainment analysis.
    *J Clin Med.* 2019;8(11):1850. PMID: **31684071**
    <https://pubmed.ncbi.nlm.nih.gov/31684071/>
    *Role in the model:* The relation in adjuvant therapy between reaching the target concentration and the risk of recurrence.

17. **Corso CR, Acco A, Bach C, Bonatto SJR, de Figueiredo BC, de Souza LM.**
    Pharmacological profile and effects of mitotane in adrenocortical carcinoma.
    *Br J Clin Pharmacol.* 2021;87(7):2698-2710. PMID: **33382119**
    <https://pubmed.ncbi.nlm.nih.gov/33382119/>
    *Role in the model:* The terminal half-life range (roughly 18–160 days), the extreme lipophilicity, and the
    o,p'-DDA / o,p'-DDE metabolites. The criterion for checking whether the model's pre-induction (86 days) /
    fully induced (42 days) half-lives fall inside this range.

18. **Hescot S, Paci A, Seck A, Slama A, Viengchareun S, Trabado S, et al.**
    Lipoprotein-Free Mitotane Exerts High Cytotoxic Activity in Adrenocortical
    Carcinoma.
    *J Clin Endocrinol Metab.* 2015;100(8):2890-2898. PMID: **26120791**
    <https://pubmed.ncbi.nlm.nih.gov/26120791/>
    *Role in the model:* The point that the fraction not bound to lipoprotein is what produces cytotoxicity —
    cell-level support for the depot view that "it is partitioning into the lipid phase, not the total
    plasma concentration, that determines the effect".

---

## 3. Mitotane's mechanism of action — adrenolysis

19. **Sbiera S, Leich E, Liebisch G, Sbiera I, Schirbel A, Wiemer L, et al.**
    Mitotane Inhibits Sterol-O-Acyl Transferase 1 Triggering Lipid-Mediated
    Endoplasmic Reticulum Stress and Apoptosis in Adrenocortical Carcinoma Cells.
    *Endocrinology.* 2015;156(11):3895-908. PMID: **26305886**
    <https://pubmed.ncbi.nlm.nih.gov/26305886/>
    *Role in the model:* The molecular mechanism of the adrenolytic pathway. The basis for the
    `SOAT1INH → OXYSTER → ERSTRESS → CASP → ADRLYSIS` chain in cluster 6 of the model and for
    `KLYS`/`EC50LYS`/`HLYS`.

20. **Weigand I, Altieri B, Lacombe AMF, Basile V, Kircher S, Landwehr LS, et al.**
    Expression of SOAT1 in Adrenocortical Carcinoma and Response to Mitotane
    Monotherapy.
    *J Clin Endocrinol Metab.* 2020;105(8):dgaa293. PMID: **32449514**
    <https://pubmed.ncbi.nlm.nih.gov/32449514/>
    *Role in the model:* Target (SOAT1) expression level is related to response → the basis for the model's
    `SOAT1exp` node and for the variability of predicted response.

21. **Smith DC, Kroiss M, Kebebew E, Habra MA, Chugh R, Schneider BJ, et al.**
    A phase 1 study of nevanimibe HCl, a novel SOAT1 inhibitor, in
    adrenocortical carcinoma.
    *Invest New Drugs.* 2020;38(5):1421-1429. PMID: **31984451**
    <https://pubmed.ncbi.nlm.nih.gov/31984451/>
    *Role in the model:* The contrasting branch in which a selective SOAT1 inhibitor stimulates only the
    adrenolytic axis (the `NEVANIM` node on the map) — the clinical basis for the counterfactual experiment
    "what happens if you get mitotane's benefit axis without the induction".

22. **Weigand I, Ronchi CL, Rizk-Rabin M, Dalmazi GD, Wild V, Bathon K, et al.**
    Active steroid hormone synthesis renders adrenocortical cells highly
    susceptible to type II ferroptosis induction.
    *Cell Death Dis.* 2020;11(3):192. PMID: **32184394**
    <https://pubmed.ncbi.nlm.nih.gov/32184394/>
    *Role in the model:* The point that cells actively synthesising steroid are more vulnerable —
    wired in the model as `CORT_SYN → FERROPT` (the mechanism by which a secreting tumour may be
    more sensitive to mitotane).

23. **Ronchi CL, Sbiera S, Volante M, Steinhauer S, Scott-Wild V, Altieri B, et al.**
    CYP2W1 is highly expressed in adrenal glands and is positively associated
    with the response to mitotane in adrenocortical carcinoma.
    *PLoS One.* 2014;9(8):e105855. PMID: **25144458**
    <https://pubmed.ncbi.nlm.nih.gov/25144458/>
    *Role in the model:* A pharmacogenomic marker predicting response (the `CYP2W1` node on the map).

24. **van Koetsveld PM, Creemers SG, Dogan F, Franssen GJH, de Herder WW, Hofland LJ, et al.**
    The Efficacy of Mitotane in Human Primary Adrenocortical Carcinoma Cultures.
    *J Clin Endocrinol Metab.* 2020;105(2):407-417. PMID: **31586196**
    <https://pubmed.ncbi.nlm.nih.gov/31586196/>
    *Role in the model:* Interindividual differences in in vitro sensitivity — the basis for putting
    response variability into the virtual population layer.

25. **Krüger AF, et al.**
    Adrenocortical Mitochondria-Associated Membranes and the Lipidoproteomic
    Response to Mitotane.
    *J Endocr Soc.* 2025;10(1):bvaf198. PMID: **41472891**
    <https://pubmed.ncbi.nlm.nih.gov/41472891/>
    *Role in the model:* The mitochondria–ER contact interface and lipid redistribution — the most recent
    evidence for the `MITOCHOND`/`ERSTRESS` coupling on the map.

---

## 4. CYP3A4 induction and drug interactions — harm axis #1

26. **Kroiss M, Quinkler M, Lutz WK, Allolio B, Fassnacht M.**
    Drug interactions with mitotane by induction of CYP3A4 metabolism in the
    clinical management of adrenocortical carcinoma.
    *Clin Endocrinol (Oxf).* 2011;75(5):585-91. PMID: **21883349**
    <https://pubmed.ncbi.nlm.nih.gov/21883349/>
    *Role in the model:* `EMAXIND`/`EC50IND` and the clinical breadth of the induction
    (the list of affected drugs in cluster 8 of the map).

27. **van Erp NP, Guchelaar HJ, Ploeger BA, Romijn JA, den Hartigh J, Gelderblom H.**
    Mitotane has a strong and a durable inducing effect on CYP3A4 activity.
    *Eur J Endocrinol.* 2011;164(4):621-6. PMID: **21220434**
    <https://pubmed.ncbi.nlm.nih.gov/21220434/>
    *Role in the model:* **The magnitude and persistence of the induction**, measured with a midazolam probe.
    The direct evidence for the behaviour in which ENZ rises up to ~4-fold and is sustained for months
    after dosing stops (result H).

28. **Theile D, Haefeli WE, Weiss J.**
    Effects of adrenolytic mitotane on drug elimination pathways assessed
    in vitro.
    *Endocrine.* 2015;49(3):842-53. PMID: **25542188**
    <https://pubmed.ncbi.nlm.nih.gov/25542188/>
    *Role in the model:* Co-induction of CYP2B6·UGT·SULT·P-gp alongside CYP3A4 →
    the `ENZ2B6`, `PGP`/`FPGPIND`, `DDI_LEVO` nodes.

29. **Kroiss M, Quinkler M, Johanssen S, van Erp NP, Lankheet N, Pöllinger A, et al.**
    Sunitinib in refractory adrenocortical carcinoma: a phase II, single-arm,
    open-label trial (SIRAC).
    *J Clin Endocrinol Metab.* 2012;97(10):3495-503. PMID: **22851488**
    <https://pubmed.ncbi.nlm.nih.gov/22851488/>
    *Role in the model:* Clinical evidence that mitotane induction markedly lowers the exposure of a
    co-administered TKI (`DDI_SUNI` on the map). A generalising case for result D.

---

## 5. Steroid replacement · binding proteins · assay pitfalls

30. **Chortis V, Taylor AE, Schneider P, Tomlinson JW, Hughes BA, O'Neil DM, et al.**
    Mitotane therapy in adrenocortical cancer induces CYP3A4 and inhibits
    5alpha-reductase, explaining the need for personalized glucocorticoid and
    androgen replacement.
    *J Clin Endocrinol Metab.* 2013;98(1):161-71. PMID: **23162091**
    <https://pubmed.ncbi.nlm.nih.gov/23162091/>
    *Role in the model:* **The direct evidence for result E.** Induction of cortisol clearance roughly
    doubles the hydrocortisone requirement (`FMHC`), and
    5α-reductase inhibition makes androgen replacement necessary in men (`INH5AR`).

31. **Nader N, Raverot G, Emptoz-Bonneton A, Déchaud H, Bonnay M, Baudin E, Pugeat M.**
    Mitotane has an estrogenic effect on sex hormone-binding globulin and
    corticosteroid-binding globulin in humans.
    *J Clin Endocrinol Metab.* 2006;91(6):2165-70. PMID: **16551731**
    <https://pubmed.ncbi.nlm.nih.gov/16551731/>
    *Role in the model:* **The other half of result E.** CBG·SHBG rise roughly
    two-fold (`EMAXCBG`, `EMAXSHBG`). The model shows quantitatively that this does not greatly change
    free cortisol but **raises the total/free ratio from 16.9 → 37.2**, misleading the total cortisol reading —
    that is, "the patient's error" (induction) and "the assay's error" (CBG) are of
    different kinds and point in the same direction.

32. **Daffara F, De Francia S, Reimondo G, Zaggia B, Aroasio E, Porpiglia F, et al.**
    Prospective evaluation of mitotane toxicity in adrenocortical cancer
    patients treated adjuvantly.
    *Endocr Relat Cancer.* 2008;15(4):1043-53. PMID: **18824557**
    <https://pubmed.ncbi.nlm.nih.gov/18824557/>
    *Role in the model:* The toxicity profile — raised GGT/ALT (`AALT`), lowered free T4
    (`AFT4`), gastrointestinal intolerance, and the effective universality of adrenal insufficiency (to be
    contrasted with the almost complete loss of adrenal cortex that `KLYS` produces).

---

## 6. Cortisol and immune escape — result F

33. **Landwehr LS, Altieri B, Schreiner J, Sbiera I, Weigand I, Kroiss M, et al.**
    Interplay between glucocorticoids and tumor-infiltrating lymphocytes on the
    prognosis of adrenocortical carcinoma.
    *J Immunother Cancer.* 2020;8(1):e000469. PMID: **32474412**
    <https://pubmed.ncbi.nlm.nih.gov/32474412/>
    *Role in the model:* **The direct evidence for result F.** Glucocorticoid excess depletes
    tumour-infiltrating lymphocytes and worsens prognosis. The model's
    `GR → GCIMSUP → TEFF(inhibition)` coupling and `GCIC50`/`HGC`.

34. **Baechle JJ, Hanna DN, Sekhar KR, Rathmell JC, Rathmell WK, Baregamian N.**
    Integrative computational immunogenomic profiling of cortisol-secreting
    adrenocortical carcinoma.
    *J Cell Mol Med.* 2021;25(21):10061-10072. PMID: **34664400**
    <https://pubmed.ncbi.nlm.nih.gov/34664400/>
    *Role in the model:* The immune-excluded phenotype of secreting ACC.

35. **Wu K, Liu X, Wang Y, Liu Y, Zhang Y, Wang Z, et al.**
    Discovery of a glucocorticoid receptor (GR) activity signature correlates
    with immune cell infiltration in adrenocortical carcinoma.
    *J Immunother Cancer.* 2023;11(10):e007528. PMID: **37793855**
    <https://pubmed.ncbi.nlm.nih.gov/37793855/>
    *Role in the model:* GR activity itself correlates with infiltration → the basis for the model using
    **GR occupancy (`GRO`)** rather than total cortisol as the driving variable of immune suppression.

36. **Landwehr LS, Schreiner J, Appenzeller S, Kircher S, Herterich S, Sbiera S, et al.**
    Expression and Prognostic Relevance of PD-1, PD-L1, and CTLA-4 Immune
    Checkpoints in Adrenocortical Carcinoma.
    *J Clin Endocrinol Metab.* 2024;109(9):2325-2334. PMID: **38415841**
    <https://pubmed.ncbi.nlm.nih.gov/38415841/>
    *Role in the model:* Checkpoint molecule expression and prognosis — the `PD1`/`CTLA4` nodes.

37. **Maier T, et al.**
    Wnt/beta-catenin pathway activation is associated with glucocorticoid
    secretion in adrenocortical carcinoma, but not directly with immune cell
    infiltration.
    *Front Endocrinol (Lausanne).* 2025;16:1502117. PMID: **40130164**
    <https://pubmed.ncbi.nlm.nih.gov/40130164/>
    *Role in the model:* The association between Wnt activity and the secreting phenotype
    (`WNT → STEROIDOG_PHEN` on the map), and the point that it is **not** connected directly to
    immune infiltration — the basis for the structural choice not to wire Wnt→immunity directly but only
    through cortisol as the mediator.

---

## 7. Genomics · prognosis

38. **Assié G, Letouzé E, Fassnacht M, Jouinot A, Luscap W, Barreau O, et al.**
    Integrated genomic characterization of adrenocortical carcinoma.
    *Nat Genet.* 2014;46(6):607-12. PMID: **24747642**
    <https://pubmed.ncbi.nlm.nih.gov/24747642/>
    *Role in the model:* Driver frequencies (cluster 1 of the map) and the molecular subtypes.

39. **Zheng S, Cherniack AD, Dewal N, Moffitt RA, Danilova L, Murray BA, et al.;
    Cancer Genome Atlas Research Network.**
    Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma.
    *Cancer Cell.* 2016;29(5):723-736. PMID: **27165744**
    <https://pubmed.ncbi.nlm.nih.gov/27165744/>
    *Role in the model:* TP53 ~21%, ZNRF3 ~19%, CTNNB1 ~16%, IGF2 overexpression
    ~90% — the figures on the genomic nodes of the map.

40. **Beuschlein F, Weigel J, Saeger W, Kroiss M, Wild V, Daffara F, et al.**
    Major prognostic role of Ki67 in localized adrenocortical carcinoma after
    complete resection.
    *J Clin Endocrinol Metab.* 2015;100(3):841-9. PMID: **25559399**
    <https://pubmed.ncbi.nlm.nih.gov/25559399/>
    *Role in the model:* The `KI67 → KGROW` coupling and the selection of candidates for adjuvant therapy.

41. **Arlt W, Biehl M, Taylor AE, Hahner S, Libé R, Hughes BA, et al.**
    Urine steroid metabolomics as a biomarker tool for detecting malignancy in
    adrenal tumors.
    *J Clin Endocrinol Metab.* 2011;96(12):3775-3784. PMID: **21917861**
    <https://pubmed.ncbi.nlm.nih.gov/21917861/>
    *Role in the model:* The `URSTER` node on the map. Why the model gives DHEAS·17-OHP·
    11-deoxycortisol their own outputs (the very pattern of substrate accumulation under enzyme
    blockade is diagnostic information).

---

## 8. Reviews · clinical guidelines

42. **Fassnacht M, Dekkers OM, Else T, Baudin E, Berruti A, de Krijger RR, et al.**
    European Society of Endocrinology Clinical Practice Guidelines on the
    management of adrenocortical carcinoma in adults, in collaboration with the
    European Network for the Study of Adrenal Tumors.
    *Eur J Endocrinol.* 2018;179(4):G1-G46. PMID: **30299884**
    <https://pubmed.ncbi.nlm.nih.gov/30299884/>
    *Role in the model:* The treatment algorithm, TDM targets, and steroid replacement recommendations —
    the clinical skeleton for the whole scenario design.

43. **Fassnacht M, Assié G, Baudin E, Eisenhofer G, de la Fouchardiere C, Haak HR, et al.**
    Adrenocortical carcinomas and malignant phaeochromocytomas: ESMO-EURACAN
    Clinical Practice Guidelines for diagnosis, treatment and follow-up.
    *Ann Oncol.* 2020;31(11):1476-1490. PMID: **32861807**
    <https://pubmed.ncbi.nlm.nih.gov/32861807/>
    *Role in the model:* The basis for the treatment sequence in advanced disease (EDP-M first line).

44. **Else T, Kim AC, Sabolch A, Raymond VM, Kandathil A, Caoili EM, et al.**
    Adrenocortical carcinoma.
    *Endocr Rev.* 2014;35(2):282-326. PMID: **24423978**
    <https://pubmed.ncbi.nlm.nih.gov/24423978/>
    *Role in the model:* A pathophysiology review — background for the map as a whole.

45. **Fassnacht M, Kroiss M, Allolio B.**
    Update in adrenocortical carcinoma.
    *J Clin Endocrinol Metab.* 2013;98(12):4551-64. PMID: **24081734**
    <https://pubmed.ncbi.nlm.nih.gov/24081734/>
    *Role in the model:* A review of clinical management.

46. **Allolio B, Fassnacht M.**
    Clinical review: Adrenocortical carcinoma: clinical update.
    *J Clin Endocrinol Metab.* 2006;91(6):2027-37. PMID: **16551738**
    <https://pubmed.ncbi.nlm.nih.gov/16551738/>
    *Role in the model:* The frequency of the hormone-hypersecreting phenotype (about 60% are functional) —
    the basis for the `SECRETOR` switch.

---

## 9. Couplings explicitly **ASSUMED** in the model

For transparency, every coupling not directly measured in the literature is collected here.

| Model element | Value | Nature of the evidence |
|---|---|---|
| `FM3A4ETO` (fraction of etoposide metabolised by CYP3A4) | 0.35 | **ASSUMED.** No mitotane–etoposide interaction study could be identified. Derived from the magnitude of CYP3A4 induction in nos. 26·27 and from the general pharmacological knowledge that etoposide is partly CYP3A4-metabolised and partly renally excreted. Consequence: at ENZ≈4, clearance about 2-fold → exposure about 1/2. Because this coefficient directly sets the size of result D, it is **the parameter that should be replaced by a measured value first**. |
| `KPFAT` (distribution volume per kg of fat, 203 L/kg) | 203 | **Partly assumed.** Back-calculated to reproduce the apparent Vss ~6,086 L of no. 11 and the BMI covariate effect. Not a measured tissue partition coefficient. |
| `AFAT` (cortisol excess → increased fat) | 0.35 | **ASSUMED.** Set approximately from the size of the body composition change in Cushing's syndrome. It sets the size of result C (the model predicts that this feedback barely operates during treatment and is concentrated in the pre-diagnosis period). |
| `IC50RB` < `IC50RA` (metabolic IR-B is more sensitive to linsitinib than tumour IR-A) | 0.35 vs 1.50 | **ASSUMED, mechanistically motivated.** The observation in no. 5 that dose-limiting hyperglycaemia appears before the antitumour effect is read as "the metabolic receptor is blocked first". It sets the size of the insulin rescue in result G. |
| `SLPIMM`, `EMAXPD1` | 0.010, 2.2 | **Calibrated values.** Set to reproduce the ORR ~14–23% of nos. 6·7. |
| `CNSTHR` = 24 mg/L (Hill EC50) | 24 | **Structural choice.** The EC50 was moved upwards so that the clinical alert level of 20 mg/L sits on the **rising limb** of the Hill curve rather than at its EC50. At 20 mg/L, about 17% of maximal injury. |
| Tumour growth `KGROW`, `TVMAX` | 0.0047, 20,000 mL | **Calibrated values.** Set so that the untreated volume doubling time is about 50–60 days. The actual growth rate of ACC is highly heterogeneous. |

---

## 10. Summary of how the model structure maps onto the literature

| Model result | Reference numbers | Character |
|---|---|---|
| Vss ~6,000 L, mostly fat | 11, 17, 18 | Directly supported |
| Therapeutic window 14–20 mg/L | 12, 15, 16 | Directly supported |
| A. Fat changes only τ, not Css | 11 (BMI covariate) | Literature + model deduction |
| B. Dose RCT vs unstratified body composition | 11 (CV 81.5%), 13 (RCT, no significant difference) | **A new explanation from the model** |
| C. Pre-diagnosis Cushing's enlarges the depot in advance | 46 (60% functional), 32 | **A new prediction from the model** |
| D. EDP-M metabolises its own partner | 26, 27, 28, 29 + ASSUMED coefficient | Mechanism established, magnitude assumed |
| E. Total cortisol lies | 30 (induction), 31 (CBG) | Directly supported, quantification from the model |
| F. The tumour premedicates itself with steroid | 33, 34, 35, 37 | Directly supported |
| G. The target is there but the drug failed | 4, 5 | Literature + a mechanistic explanation from the model |
| H. There is no washout | 27, 17 | Directly supported |

---

## 11. Tools

- **mrgsolve** — Baron KT. *mrgsolve: Simulate from ODE-Based Models.*
  <https://mrgsolve.org/>
- **Graphviz** — <https://graphviz.org/>
- **Shiny** — <https://shiny.posit.co/>
- Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO.
  Model of chemotherapy-induced myelosuppression with parameter consistency
  across drugs. *J Clin Oncol.* 2002;20(24):4713-21. PMID: **12488418**
  <https://pubmed.ncbi.nlm.nih.gov/12488418/>
  *Role in the model:* The 4-transit compartment structure for neutrophils (`PROLN`–`TR3N`–`CIRCN`),
  `MTTN`/`GAMN`.

---

> **Disclaimer.** This model and its reference list are for educational and research purposes. Do not use them for clinical
> decision-making, prescribing, or regulatory submission. The parameters are illustrative approximations, and
> separate fitting and validation against real patient data would be required.
