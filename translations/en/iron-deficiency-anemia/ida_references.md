# Iron Deficiency Anaemia — QSP model references

This document collects, section by section, the literature underlying the structure and
parameters of `ida_qsp_model.dot` / `ida_mrgsolve_model.R` / `ida_shiny_app.R`. The `→` mark on each
entry shows **which part** of the model that reference supports.

> ⚠️ This model is for educational and research purposes. Its parameters are approximations based on the literature and
> cannot be used directly in the care of an individual patient.

---

## 1. Epidemiology · diagnosis · clinical overview

1. Camaschella C. **Iron deficiency.** *Blood* 2019;133(1):30-39.
   <https://pubmed.ncbi.nlm.nih.gov/30401704/>
   → the overall pathophysiological skeleton, and the distinction between absolute and functional iron deficiency.
2. Camaschella C. **Iron-deficiency anemia.** *N Engl J Med* 2015;372(19):1832-1843.
   <https://pubmed.ncbi.nlm.nih.gov/25946282/>
   → the diagnostic algorithm, the ferritin <15 ng/mL criterion, and the treatment response to be expected.
3. Pasricha SR, Tye-Din J, Muckenthaler MU, Swinkels DW. **Iron deficiency.**
   *Lancet* 2021;397(10270):233-248.
   <https://pubmed.ncbi.nlm.nih.gov/33285139/>
   → classification by cause, and the criteria for choosing between oral and intravenous iron.
4. Kassebaum NJ, et al. **A systematic analysis of global anemia burden from 1990
   to 2010.** *Blood* 2014;123(5):615-624.
   <https://pubmed.ncbi.nlm.nih.gov/24297872/>
   → evidence for the scale: close to half the global anaemia burden is iron deficiency.
5. Peyrin-Biroulet L, Williet N, Cacoub P. **Guidelines on the diagnosis and
   treatment of iron deficiency across indications: a systematic review.**
   *Am J Clin Nutr* 2015;102(6):1585-1594.
   <https://pubmed.ncbi.nlm.nih.gov/26561626/>
   → the heterogeneity of the diagnostic thresholds (ferritin · TSAT) between indications.
6. Weiss G, Ganz T, Goodnough LT. **Anemia of inflammation.**
   *Blood* 2019;133(1):40-50.
   <https://pubmed.ncbi.nlm.nih.gov/30401705/>
   → the basis of the inflammation scenario (IL6_IN): the dissociation of ferritin↑ with TSAT↓.

## 2. The physiology of iron absorption — the enterocyte gate

7. Frazer DM, Anderson GJ. **The regulation of iron transport.**
   *Biofactors* 2014;40(2):206-214.
   <https://pubmed.ncbi.nlm.nih.gov/24132807/>
   → the DMT1–ferroportin axis and the role of enterocyte ferritin (the A_ENT / A_EFT structure).
8. Gunshin H, et al. **Cloning and characterization of a mammalian
   proton-coupled metal-ion transporter (DMT1).** *Nature* 1997;388(6641):482-488.
   <https://pubmed.ncbi.nlm.nih.gov/9242408/>
   → identification of the apical Fe²⁺ transporter — what VMAX_D / KM_D refer to.
9. McKie AT, et al. **An iron-regulated ferric reductase associated with the
   absorption of dietary iron (DCYTB).** *Science* 2001;291(5509):1755-1759.
   <https://pubmed.ncbi.nlm.nih.gov/11230685/>
   → the Fe³⁺→Fe²⁺ reduction step (dependence on gastric acid and ascorbate).
10. Donovan A, et al. **The iron exporter ferroportin/Slc40a1 is essential for
    iron homeostasis.** *Cell Metab* 2005;1(3):191-200.
    <https://pubmed.ncbi.nlm.nih.gov/16054062/>
    → the only route of basolateral export = the justification for the FPN_ENT factor.
11. Vulpe CD, et al. **Hephaestin, a ceruloplasmin homologue implicated in
    intestinal iron transport.** *Nat Genet* 1999;21(2):195-199.
    <https://pubmed.ncbi.nlm.nih.gov/9988272/>
    → the re-oxidation step needed to load transferrin.
12. Brasse-Lagnel C, et al. **Intestinal DMT1 cotransporter is down-regulated
    by hepcidin via proteasome internalization.**
    *Gastroenterology* 2011;140(4):1261-1271.
    <https://pubmed.ncbi.nlm.nih.gov/21199652/>
    → the direct apical action of hepcidin = the parameter `K_APIC`.
13. Hahn PF, Bale WF, Ross JF, Balfour WM, Whipple GH. **Radioactive iron
    absorption by gastro-intestinal tract: influence of anemia, anoxia, and
    antecedent feeding.** *J Exp Med* 1943;78(3):169-188.
    <https://pubmed.ncbi.nlm.nih.gov/19871320/>
    → the original source of the "mucosal block" — the historical basis of THESIS 1 of this model.
14. Hurrell R, Egli I. **Iron bioavailability and dietary reference values.**
    *Am J Clin Nutr* 2010;91(5):1461S-1467S.
    <https://pubmed.ncbi.nlm.nih.gov/20200263/>
    → the bioavailability of dietary iron (the range underlying `F_DIET` = 0.30).
15. Lynch S, et al. **Biomarkers of Nutrition for Development (BOND) — Iron
    review.** *J Nutr* 2018;148(suppl_1):1001S-1067S.
    <https://pubmed.ncbi.nlm.nih.gov/29878148/>
    → a synthesis of absorption rates, body iron distribution and biomarker reference values.

## 3. Hepcidin — the controller

16. Ganz T. **Systemic iron homeostasis.** *Physiol Rev* 2013;93(4):1721-1741.
    <https://pubmed.ncbi.nlm.nih.gov/24137020/>
    → a quantitative map of systemic iron flux: recycling 20-25 mg/day vs absorption 1-2 mg/day.
    The key basis of THESIS 3 of this model.
17. Nemeth E, et al. **Hepcidin regulates cellular iron efflux by binding to
    ferroportin and inducing its internalization.**
    *Science* 2004;306(5704):2090-2093.
    <https://pubmed.ncbi.nlm.nih.gov/15514116/>
    → hepcidin→ferroportin degradation = the mechanism behind `KDEG_FPE` / `KDEG_FPR`.
18. Nemeth E, et al. **IL-6 mediates hypoferremia of inflammation by inducing
    the synthesis of the iron regulatory hormone hepcidin.**
    *J Clin Invest* 2004;113(9):1271-1276.
    <https://pubmed.ncbi.nlm.nih.gov/15124018/>
    → `EMAX_IL6` / `KM_IL6` (the STAT3 pathway).
19. Andriopoulos B Jr, et al. **BMP6 is a key endogenous regulator of hepcidin
    expression and iron metabolism.** *Nat Genet* 2009;41(4):482-487.
    <https://pubmed.ncbi.nlm.nih.gov/19252486/>
    → sensing of stored iron (BMP6/SMAD) = `HEP_STORE_E` / `KM_HEP_LIV`.
20. Kautz L, et al. **Identification of erythroferrone as an erythroid regulator
    of iron metabolism.** *Nat Genet* 2014;46(7):678-684.
    <https://pubmed.ncbi.nlm.nih.gov/24880340/>
    → the basis of the ERFE state variable and of `KI_ERFE`.
21. Arezes J, et al. **Erythroferrone inhibits the induction of hepcidin by BMP6.**
    *Blood* 2018;132(14):1473-1477.
    <https://pubmed.ncbi.nlm.nih.gov/30097509/>
    → how ERFE sequesters BMP (the f_erfe multiplicative term in the model).
22. Finberg KE, et al. **Mutations in TMPRSS6 cause iron-refractory iron
    deficiency anemia (IRIDA).** *Nat Genet* 2008;40(5):569-571.
    <https://pubmed.ncbi.nlm.nih.gov/18408718/>
    → the `TMPRSS6` parameter and the IRIDA scenario.
23. Silvestri L, et al. **The serine protease matriptase-2 (TMPRSS6) inhibits
    hepcidin activation by cleaving membrane hemojuvelin.**
    *Cell Metab* 2008;8(6):502-511.
    <https://pubmed.ncbi.nlm.nih.gov/18976966/>
    → the matriptase-2 → HJV cleavage mechanism.
24. Ganz T, et al. **Immunoassay for human serum hepcidin.**
    *Blood* 2008;112(10):4292-4297.
    <https://pubmed.ncbi.nlm.nih.gov/18689548/>
    → the normal range of serum hepcidin (compare the model: healthy 6.4 ng/mL, IDA 0.32 ng/mL).
25. Girelli D, Nemeth E, Swinkels DW. **Hepcidin in the diagnosis of iron
    disorders.** *Blood* 2016;127(23):2809-2813.
    <https://pubmed.ncbi.nlm.nih.gov/27044621/>
    → the clinical use of hepcidin to predict oral iron responsiveness (NOTE E of the model).

## 4. Dosing interval · oral iron pharmacokinetics — the central thesis

26. Moretti D, et al. **Oral iron supplements increase hepcidin and decrease
    iron absorption from daily or twice-daily doses in iron-depleted young
    women.** *Blood* 2015;126(17):1981-1989.
    <https://pubmed.ncbi.nlm.nih.gov/26289639/>
    → direct evidence that the hepcidin rise after a 60 mg dose persists for 24 hours and reduces the absorption of a
    second dose given the same day. Model check: probe +4 h = 75.4 %, +24 h = 91.5 %.
27. Stoffel NU, et al. **Iron absorption from oral iron supplements given on
    consecutive versus alternate days and as single morning doses versus
    twice-daily split dosing in iron-depleted women: two open-label,
    randomised controlled trials.** *Lancet Haematol* 2017;4(11):e524-e533.
    <https://pubmed.ncbi.nlm.nih.gov/29032957/>
    → the fractional-absorption advantage of alternate-day dosing. Observed daily:alternate ratio ≈ 0.75,
    this model 0.87 (the model is conservative — stated as a limitation).
28. Stoffel NU, et al. **Iron absorption from supplements is greater with
    alternate day than with consecutive day dosing in iron-deficient anemic
    women.** *Haematologica* 2020;105(5):1232-1239.
    <https://pubmed.ncbi.nlm.nih.gov/31413088/>
    → the same conclusion in women who are anaemic (not simply depleted) — the basis for simulating from this
    model's IDA baseline state.
29. Kaundal R, Bhatia P, Jain A, et al. **Randomized controlled trial of
    twice-daily versus alternate-day oral iron therapy in the treatment of
    iron-deficiency anemia.** *Ann Hematol* 2020;99(1):57-63.
    <https://pubmed.ncbi.nlm.nih.gov/31811360/>
    → verification at the level of clinical outcome: alternate-day dosing produced an Hb response similar to
    twice-daily dosing with fewer adverse effects. The clinical counterpart of THESIS 2 — that total absorbed amount
    and fractional absorption rank in different orders.
30. Rimon E, et al. **Are we giving too much iron? Low-dose iron therapy is
    effective in octogenarians.** *Am J Med* 2005;118(10):1142-1147.
    <https://pubmed.ncbi.nlm.nih.gov/16194646/>
    → low doses (15-50 mg) give an Hb response similar to high doses — in keeping with the model's prediction
    that fractional absorption saturates (FIA 25.2 % at 30 mg).
31. Schrier SL. **So you know how to treat iron deficiency anemia.**
    *Blood* 2015;126(17):1971.
    <https://pubmed.ncbi.nlm.nih.gov/26494915/>
    → a commentary on the Moretti study: the argument for changing clinical prescribing practice.
32. Zimmermann MB, Hurrell RF. **Nutritional iron deficiency.**
    *Lancet* 2007;370(9586):511-520.
    <https://pubmed.ncbi.nlm.nih.gov/17693180/>
    → dietary interactions such as ascorbate and phytate.
33. Tolkien Z, Stecher L, Mander AP, Pereira DI, Powell JJ. **Ferrous sulfate
    supplementation causes significant gastrointestinal side-effects in
    adults: a systematic review and meta-analysis.**
    *PLoS One* 2015;10(2):e0117383.
    <https://pubmed.ncbi.nlm.nih.gov/25700159/>
    → the basis of the GI adverse-effect and adherence terms (`K_GI`, `EMAX_ADH`).
34. Gasche C, et al. **Ferric maltol is effective in correcting iron deficiency
    anemia in patients with inflammatory bowel disease (AEGIS trials).**
    *Inflamm Bowel Dis* 2015;21(3):579-588.
    <https://pubmed.ncbi.nlm.nih.gov/25545376/>
    → an alternative oral formulation (ferric maltol) — the oral pharmacology cluster of the map.

## 5. Intravenous iron — pharmacokinetics and safety

35. Auerbach M, Macdougall I. **The available intravenous iron formulations:
    History, efficacy, and toxicology.** *Hemodial Int* 2017;21 Suppl 1:S83-S92.
    <https://pubmed.ncbi.nlm.nih.gov/28371203/>
    → the carbohydrate shell, RES processing and free-iron release of each formulation (A_COL / F_COL_RES / K_COL).
36. Geisser P, Burckhardt S. **The pharmacokinetics and pharmacodynamics of
    iron preparations.** *Pharmaceutics* 2011;3(1):12-33.
    <https://pubmed.ncbi.nlm.nih.gov/24310424/>
    → the plasma half-life of FCM, 7-12 h (`K_COL` = 0.0693/h ≈ t½ 10 h).
37. Onken JE, et al. **A multicenter, randomized, active-controlled study to
    investigate the efficacy and safety of ferric carboxymaltose in patients
    with iron deficiency anemia (REPAIR-IDA).**
    *Transfusion* 2014;54(2):306-315.
    <https://pubmed.ncbi.nlm.nih.gov/23772856/>
    → the size of the Hb response to the FCM 750 mg×2 regimen (scenario 9 of the model).
38. Van Wyck DB, Mangione A, Morrison J, Hadley PE, Jehle JA, Goodnough LT.
    **Large-dose intravenous ferric carboxymaltose injection for iron
    deficiency anemia in heavy uterine bleeding.**
    *Transfusion* 2009;49(12):2719-2728.
    <https://pubmed.ncbi.nlm.nih.gov/19682342/>
    → high-dose FCM in IDA due to uterine bleeding — the same baseline patient as in this model.
39. Rognoni C, Venturini S, Meregaglia M, Marmifero M, Tarricone R.
    **Efficacy and safety of ferric carboxymaltose and other formulations in
    iron-deficient patients: a systematic review and network meta-analysis.**
    *Clin Drug Investig* 2016;36(3):177-194.
    <https://pubmed.ncbi.nlm.nih.gov/26692005/>
    → the similarity of the Hb response across formulations — the premise of THESIS 4 ("equivalent efficacy, different price").
40. Rampton D, et al. **Hypersensitivity reactions to intravenous iron:
    guidance for risk minimization and management.**
    *Haematologica* 2014;99(11):1671-1676.
    <https://pubmed.ncbi.nlm.nih.gov/25420283/>
    → CARPA and the Fishbane reaction (the safety nodes of the map).

## 6. FGF23 — the phosphate axis (the carboxymaltose phosphate penalty)

41. Wolf M, et al. **Effects of iron isomaltoside vs ferric carboxymaltose on
    hypophosphatemia in iron-deficiency anemia: two randomized clinical trials
    (PHOSPHARE-IDA).** *JAMA* 2020;323(5):432-443.
    <https://pubmed.ncbi.nlm.nih.gov/32016310/>
    → the main calibration basis of the FGF23-phosphate axis of this model: the incidence of phosphate <2.0 mg/dL is
    markedly higher with FCM and rare with derisomaltose. Model: FCM 1000 mg → phosphate nadir
    1.98 mg/dL (6.8 days), <2.0 mg/dL for 4.2 days; derisomaltose → 3.02 mg/dL, 0 days.
42. Wolf M, Koch TA, Bregman DB. **Effects of iron deficiency anemia and its
    treatment on fibroblast growth factor 23 and phosphate homeostasis in
    women.** *J Bone Miner Res* 2013;28(8):1793-1803.
    <https://pubmed.ncbi.nlm.nih.gov/23505057/>
    → the rise in iFGF23 is specific to FCM and peaks at 24-48 hours (`E_CLV`, `KOUT_CLV`).
43. Schouten BJ, Hunt PJ, Livesey JH, Frampton CM, Soule SG.
    **FGF23 elevation and hypophosphatemia after intravenous iron
    polymaltose: a prospective study.**
    *J Clin Endocrinol Metab* 2009;94(7):2332-2337.
    <https://pubmed.ncbi.nlm.nih.gov/19366850/>
    → an early observation of the difference between formulations.
44. Zoller H, Schaefer B, Glodny B. **Iron-induced hypophosphatemia: an
    emerging complication.** *Curr Opin Nephrol Hypertens* 2017;26(4):266-275.
    <https://pubmed.ncbi.nlm.nih.gov/28399017/>
    → the cleavage-inhibition hypothesis — the structure of the CLV state variable in the model.
45. Shimada T, et al. **FGF-23 is a potent regulator of vitamin D metabolism
    and phosphate homeostasis.** *J Bone Miner Res* 2004;19(3):429-435.
    <https://pubmed.ncbi.nlm.nih.gov/15040831/>
    → inhibition of NaPi-2a/2c + inhibition of CYP27B1 (`IMAX_FGF_P`, `IMAX_FGF_D`).
46. Klein K, et al. **Severe FGF23-based hypophosphataemic osteomalacia due to
    ferric carboxymaltose administration.** *BMJ Case Rep* 2018;2018:bcr-2017-222851.
    <https://pubmed.ncbi.nlm.nih.gov/29298794/>
    → the clinical consequence of prolonged hypophosphataemia (osteomalacia) — the OSTEOMAL node of the map.

## 7. Erythropoiesis · iron-restricted erythropoiesis

47. Koury MJ, Ponka P. **New insights into erythropoiesis: the roles of folate,
    vitamin B12, and iron.** *Annu Rev Nutr* 2004;24:105-131.
    <https://pubmed.ncbi.nlm.nih.gov/15189115/>
    → the inefficiency of erythropoiesis (apoptosis) under iron restriction — `KAPO_MAX`.
48. Rivella S. **Iron metabolism under conditions of ineffective
    erythropoiesis in beta-thalassemia.** *Blood* 2019;133(1):51-58.
    <https://pubmed.ncbi.nlm.nih.gov/30401707/>
    → iron recycling in ineffective erythropoiesis (in the model: apo_fe → A_RES).
49. Ponka P. **Tissue-specific regulation of iron metabolism and heme
    synthesis: distinct control mechanisms in erythroid cells.**
    *Blood* 1997;89(1):1-25.
    <https://pubmed.ncbi.nlm.nih.gov/8978272/>
    → the ALAS2 · mitoferrin · ferrochelatase pathway, and ZPP formation.
50. Shayeghi M, et al. **Identification of an intestinal heme transporter.**
    *Cell* 2005;122(5):789-801.
    <https://pubmed.ncbi.nlm.nih.gov/16143108/>
    → the haem iron absorption route (the HCP1 node of the map).
51. Brugnara C. **Iron deficiency and erythropoiesis: new diagnostic
    approaches.** *Clin Chem* 2003;49(10):1573-1578.
    <https://pubmed.ncbi.nlm.nih.gov/14500582/>
    → reticulocyte haemoglobin content (CHr) — the model's `CHR` output (IDA 20.5 pg).
52. Thomas C, Thomas L. **Biochemical markers and hematologic indices in the
    diagnosis of functional iron deficiency.**
    *Clin Chem* 2002;48(7):1066-1076.
    <https://pubmed.ncbi.nlm.nih.gov/12089176/>
    → the sTfR / log(ferritin) index (the `STFR` state variable).
53. Beguin Y. **Soluble transferrin receptor for the evaluation of
    erythropoiesis and iron status.** *Clin Chim Acta* 2003;329(1-2):9-22.
    <https://pubmed.ncbi.nlm.nih.gov/12589962/>
    → the normal range of sTfR and the size of its rise in iron deficiency (`STFR0`, `STFR_E`).

## 8. Non-erythroid iron function · symptoms

54. Haas JD, Brownlie T 4th. **Iron deficiency and reduced work capacity:
    a critical review of the research to determine a causal relationship.**
    *J Nutr* 2001;131(2S-2):676S-690S.
    <https://pubmed.ncbi.nlm.nih.gov/11160598/>
    → tissue iron (myoglobin, Fe-S enzymes) and exercise capacity — the basis of `W_FACIT_TISS`.
55. Allen RP, et al. **Restless legs syndrome/Willis-Ekbom disease
    pathophysiology.** *Sleep Med Clin* 2015;10(3):207-214.
    <https://pubmed.ncbi.nlm.nih.gov/26329429/>
    → the brain iron deficiency-dopamine axis → the IRLS endpoint.
56. Krayenbuehl PA, Battegay E, Breymann C, Furrer J, Schulthess G.
    **Intravenous iron for the treatment of fatigue in nonanemic,
    premenopausal women with low serum ferritin concentration.**
    *Blood* 2011;118(12):3222-3227.
    <https://pubmed.ncbi.nlm.nih.gov/21705493/>
    → repletion of tissue iron improves fatigue even without anaemia — the clinical counterpart of NOTE H of the
    model (tissue iron recovers later than, and differently from, Hb).
57. Falkingham M, et al. **The effects of oral iron supplementation on
    cognition in older children and adults: a systematic review and
    meta-analysis.** *Nutr J* 2010;9:4.
    <https://pubmed.ncbi.nlm.nih.gov/20100340/>
    → the cognitive endpoint.
58. Lozoff B, Beard J, Connor J, Barbara F, Georgieff M, Schallert T.
    **Long-lasting neural and behavioral effects of iron deficiency in
    infancy.** *Nutr Rev* 2006;64(5 Pt 2):S34-S43.
    <https://pubmed.ncbi.nlm.nih.gov/16770951/>
    → neurodevelopmental outcome in infants.

## 9. Indication-specific trials

59. Anker SD, et al. **Ferric carboxymaltose in patients with heart failure and
    iron deficiency (FAIR-HF).** *N Engl J Med* 2009;361(25):2436-2448.
    <https://pubmed.ncbi.nlm.nih.gov/19920054/>
    → functional improvement independent of whether anaemia is present — the clinical evidence for the tissue iron axis.
60. Ponikowski P, et al. **Ferric carboxymaltose for iron deficiency at
    discharge after acute heart failure (AFFIRM-AHF).**
    *Lancet* 2020;396(10266):1895-1904.
    <https://pubmed.ncbi.nlm.nih.gov/33197395/>
    → the readmission-reduction endpoint.
61. Muñoz M, et al. **International consensus statement on the peri-operative
    management of anaemia and iron deficiency.**
    *Anaesthesia* 2017;72(2):233-247.
    <https://pubmed.ncbi.nlm.nih.gov/27996086/>
    → the preoperative optimisation and transfusion-avoidance endpoints.
62. Pavord S, et al. **UK guidelines on the management of iron deficiency in
    pregnancy.** *Br J Haematol* 2020;188(6):819-830.
    <https://pubmed.ncbi.nlm.nih.gov/31578718/>
    → the requirement in pregnancy (+1000 mg) and the treatment threshold.
63. Dignass AU, Gasche C, et al. **European consensus on the diagnosis and
    management of iron deficiency and anaemia in inflammatory bowel diseases.**
    *J Crohns Colitis* 2015;9(3):211-222.
    <https://pubmed.ncbi.nlm.nih.gov/25518052/>
    → the recommendation to avoid oral iron in inflammatory bowel disease — the clinical counterpart of NOTE E of the model.
64. Munro MG, Mast AE, Powers JM, et al. **The relationship between heavy
    menstrual bleeding, iron deficiency, and iron deficiency anemia.**
    *Am J Obstet Gynecol* 2023;229(1):1-9.
    <https://pubmed.ncbi.nlm.nih.gov/36706856/>
    → the clinical basis for parameterising `VBLEED` (the >80 mL/cycle definition, and conversion to iron loss).

## 10. QSP · modelling methodology

65. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and
    physiologically-based pharmacokinetic modeling with mrgsolve: a
    hands-on tutorial.** *CPT Pharmacometrics Syst Pharmacol* 2019;8(12):883-893.
    <https://pubmed.ncbi.nlm.nih.gov/31652028/>
    → mrgsolve implementation conventions.
66. Parmar JH, Mendes P. **A computational model to understand mouse iron
    physiology and disease.** *PLoS Comput Biol* 2019;15(1):e1006680.
    <https://pubmed.ncbi.nlm.nih.gov/30608934/>
    → a precedent for an ODE model of systemic iron metabolism (the basis for comparing compartment structures).
67. Sarkar J, Potdar AA, Saidel GM. **Whole-body iron transport and metabolism:
    mechanistic, multi-scale model to improve treatment of anemia in chronic
    kidney disease.** *PLoS Comput Biol* 2018;14(4):e1006060.
    <https://pubmed.ncbi.nlm.nih.gov/29659573/>
    → a multiscale iron transport model — for comparison with this model's recycling/storage compartment design.
68. Enculescu M, et al. **Modelling systemic iron regulation during dietary
    iron overload and acute inflammation: role of hepcidin-independent
    mechanisms.** *PLoS Comput Biol* 2017;13(1):e1005322.
    <https://pubmed.ncbi.nlm.nih.gov/28068331/>
    → the separation of hepcidin-dependent from hepcidin-independent mechanisms — this model's design of two ferroportin clocks.
69. Ganz T, Nemeth E. **Hepcidin and iron homeostasis.**
    *Biochim Biophys Acta* 2012;1823(9):1434-1443.
    <https://pubmed.ncbi.nlm.nih.gov/22306005/>
    → a synthesis of quantitative parameters (hepcidin half-life, ferroportin recovery time).

---

## Calibration summary — the anchors and where they came from

| Model prediction | Observed / literature value | Supporting reference |
|---|---|---|
| Fractional absorption of a single 60 mg dose 22.0 % | about 20-22 % in iron-deficient women | 27, 28 |
| Hepcidin still 1.36-fold at 24 h | the rise persists to 24 h | 26 |
| +4 h probe absorption 75.4 % | absorption falls when two doses are given the same day | 26 |
| daily : alternate fractional-absorption ratio 0.87 | observed 0.75 (the model is conservative) | 27 |
| Serum iron peak 5.8 h | 2-5 h (the model is somewhat late) | 26 |
| Reticulocyte peak 12-14 days | 7-14 days | 2 |
| FCM 1000 mg → phosphate nadir 1.98 mg/dL (6.8 days) | nadir ~1.9-2.0 mg/dL, 7-14 days | 41, 42 |
| derisomaltose → no change in iFGF23 | a phenomenon specific to FCM | 41, 44 |
| Marrow iron consumption 13.3 mg/day (recycling ~85 %) | 20-25 mg/day (consistent once corrected for body size) | 16 |
| Hepcidin, normal 6.4 / IDA 0.32 ng/mL | normal 1-20, IDA below the limit of detection | 24, 25 |
| Mass conservation error ~1e-12 % | — | (internal audit, NOTE I) |

## Documented limitations

1. **The benefit of alternate-day dosing is smaller than observed.** The model's daily:alternate
      fractional-absorption ratio is 0.87, whereas the value Stoffel 2017 observed is about 0.75. That is, this model
      **underestimates** the refractory penalty. Holding the circulating half-life of hepcidin (2.5 h)
      fixed and generating the 24-hour memory from the recovery clock of enterocyte export capacity alone puts
      depth and persistence in conflict.
2. **The serum iron peak is late** (model 5.8 h vs observed 2-5 h): the absorption window (`KTR`) and
      basolateral export are connected in series, so the absorption delays accumulate.
3. **The IRIDA phenotype is mild.** The baseline Hb obtained with `TMPRSS6` = 0.30 is
      12.8 g/dL, milder than the moderate microcytic anaemia of real IRIDA. The qualitative conclusion of oral
      refractoriness (oral +0.32 vs +2.59 g/dL) still holds.
4. **The GI adverse-effect–adherence terms are semi-quantitative.** `K_GI` and `EMAX_ADH` were merely set to a
      size consistent with the adverse-effect incidence of the meta-analysis (33); they are not calibrated to predict
      the discontinuation rate of an individual patient.
5. It is calibrated to **a single body size** (a 60 kg woman). Children, pregnancy, CKD and so on require
      `BV_L`, `PV_L`, `VBLEED`, `IL6_IN` to be changed and then **re-equilibrated**
      (the "re-equilibrate baseline state" option in the Shiny app).
