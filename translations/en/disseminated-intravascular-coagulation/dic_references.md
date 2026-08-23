# Disseminated Intravascular Coagulation (DIC) — References
# Disseminated Intravascular Coagulation — Reference List

This list sets out, section by section, the evidence for the structure and parameters of
`dic_qsp_model.dot` / `dic_mrgsolve_model.R`. The PMID link on each entry goes to PubMed.
For which reference each equation of the model rests on, see the map in §12.

---

## 1. Definition · diagnostic criteria · epidemiology

1. Taylor FB Jr, Toh CH, Hoots WK, Wada H, Levi M. **Towards definition, clinical and laboratory criteria, and a scoring system for disseminated intravascular coagulation.** *Thromb Haemost.* 2001;86(5):1327-30. — the original ISTH overt-DIC score. Implemented directly as the model's `ISTH` output. [PMID 11816725](https://pubmed.ncbi.nlm.nih.gov/11816725/)
2. Levi M, Ten Cate H. **Disseminated intravascular coagulation.** *N Engl J Med.* 1999;341(8):586-92. [PMID 10451465](https://pubmed.ncbi.nlm.nih.gov/10451465/)
3. Levi M, van der Poll T. **Coagulation and sepsis.** *Thromb Res.* 2017;149:38-44. [PMID 27886531](https://pubmed.ncbi.nlm.nih.gov/27886531/)
4. Iba T, Nisio MD, Levy JH, Kitamura N, Thachil J. **New criteria for sepsis-induced coagulopathy (SIC) following the revised sepsis definition.** *BMJ Open.* 2017;7(9):e017046. — the SIC score. The reason fibrinogen was deliberately excluded is derived in the model. [PMID 28963294](https://pubmed.ncbi.nlm.nih.gov/28963294/)
5. Iba T, Levy JH, Warkentin TE, Thachil J, van der Poll T, Levi M. **Diagnosis and management of sepsis-induced coagulopathy and disseminated intravascular coagulation.** *J Thromb Haemost.* 2019;17(11):1989-94. [PMID 31410983](https://pubmed.ncbi.nlm.nih.gov/31410983/)
6. Gando S, Iba T, Eguchi Y, et al. **A multicenter, prospective validation of disseminated intravascular coagulation diagnostic criteria for critically ill patients: comparing current criteria.** *Crit Care Med.* 2006;34(3):625-31. — the JAAM criteria. [PMID 16521260](https://pubmed.ncbi.nlm.nih.gov/16521260/)
7. Gando S, Levi M, Toh CH. **Disseminated intravascular coagulation.** *Nat Rev Dis Primers.* 2016;2:16037. [PMID 27250996](https://pubmed.ncbi.nlm.nih.gov/27250996/)
8. Singer M, Deutschman CS, Seymour CW, et al. **The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3).** *JAMA.* 2016;315(8):801-10. — the SOFA definition. [PMID 26903338](https://pubmed.ncbi.nlm.nih.gov/26903338/)

## 2. Tissue factor · thrombin generation

9. Mackman N. **Role of tissue factor in hemostasis, thrombosis, and vascular development.** *Arterioscler Thromb Vasc Biol.* 2004;24(6):1015-22. [PMID 15117736](https://pubmed.ncbi.nlm.nih.gov/15117736/)
10. Osterud B, Bjorklid E. **Sources of tissue factor.** *Semin Thromb Hemost.* 2006;32(1):11-23. [PMID 16479458](https://pubmed.ncbi.nlm.nih.gov/16479458/)
11. Hemker HC, Giesen P, Al Dieri R, et al. **Calibrated automated thrombin generation measurement in clotting plasma.** *Pathophysiol Haemost Thromb.* 2003;33(1):4-15. — supplies the quantitative scale of the thrombin generation curve (peak 200-400 nM). [PMID 12853707](https://pubmed.ncbi.nlm.nih.gov/12853707/)
12. Mann KG, Butenas S, Brummel K. **The dynamics of thrombin formation.** *Arterioscler Thromb Vasc Biol.* 2003;23(1):17-25. [PMID 12524220](https://pubmed.ncbi.nlm.nih.gov/12524220/)
13. Hoffman M, Monroe DM 3rd. **A cell-based model of hemostasis.** *Thromb Haemost.* 2001;85(6):958-65. — the basis for the model holding `PSURF` (the activated-platelet phospholipid surface) as a separate state variable. [PMID 11434702](https://pubmed.ncbi.nlm.nih.gov/11434702/)

## 3. Endothelium · glycocalyx · thrombomodulin

14. Esmon CT. **The protein C pathway.** *Chest.* 2003;124(3 Suppl):26S-32S. [PMID 12970121](https://pubmed.ncbi.nlm.nih.gov/12970121/)
15. Faust SN, Levin M, Harrison OB, et al. **Dysfunction of endothelial protein C activation in severe meningococcal sepsis.** *N Engl J Med.* 2001;345(6):408-16. — direct histological evidence that TM/EPCR expression is lost in sepsis, so that protein C activation itself fails. The basis for the four-term product in the model (THR x TM x EPCR x PC). [PMID 11496851](https://pubmed.ncbi.nlm.nih.gov/11496851/)
16. Ait-Oufella H, Maury E, Lehoux S, Guidet B, Offenstadt G. **The endothelium: physiological functions and role in microcirculatory failure during severe sepsis.** *Intensive Care Med.* 2010;36(8):1286-98. [PMID 20443110](https://pubmed.ncbi.nlm.nih.gov/20443110/)
17. Chappell D, Jacob M, Becker BF, Hofmann-Kiefer K, Conzen P, Rehm M. **Expedition glycocalyx: a newly discovered "Great Barrier Reef".** *Anaesthesist.* 2008;57(10):959-69. [PMID 18773178](https://pubmed.ncbi.nlm.nih.gov/18773178/)
18. Ostrowski SR, Johansson PI. **Endothelial glycocalyx degradation induces endogenous heparinization in patients with severe injury and early traumatic coagulopathy.** *J Trauma Acute Care Surg.* 2012;73(1):60-6. — syndecan-1 shedding. The model's `GLX` and its capillary leak term. [PMID 22743373](https://pubmed.ncbi.nlm.nih.gov/22743373/)
19. Ono T, Mimuro J, Madoiwa S, et al. **Severe secondary deficiency of von Willebrand factor-cleaving protease (ADAMTS13) in patients with sepsis-induced disseminated intravascular coagulation.** *Blood.* 2006;107(2):528-34. [PMID 16189276](https://pubmed.ncbi.nlm.nih.gov/16189276/)

## 4. NETs · histones · thromboinflammation

20. Xu J, Zhang X, Pelayo R, et al. **Extracellular histones are major mediators of death in sepsis.** *Nat Med.* 2009;15(11):1318-21. — circulating histone concentrations (10-50 µg/mL) and cytotoxicity. The model's `HIST` and its downstream effects. [PMID 19855397](https://pubmed.ncbi.nlm.nih.gov/19855397/)
21. Fuchs TA, Brill A, Duerschmied D, et al. **Extracellular DNA traps promote thrombosis.** *Proc Natl Acad Sci USA.* 2010;107(36):15880-5. [PMID 20798043](https://pubmed.ncbi.nlm.nih.gov/20798043/)
22. Semeraro F, Ammollo CT, Morrissey JH, et al. **Extracellular histones promote thrombin generation through platelet-dependent mechanisms: involvement of platelet TLR2 and TLR4.** *Blood.* 2011;118(7):1952-61. — the basis for histones being the principal driving term for platelet activation in the model. [PMID 21673343](https://pubmed.ncbi.nlm.nih.gov/21673343/)
23. Ammollo CT, Semeraro F, Xu J, Esmon NL, Esmon CT. **Extracellular histones increase plasma thrombin generation by impairing thrombomodulin-dependent protein C activation.** *J Thromb Haemost.* 2011;9(9):1795-803. [PMID 21711444](https://pubmed.ncbi.nlm.nih.gov/21711444/)
24. Engelmann B, Massberg S. **Thrombosis as an intravascular effector of innate immunity.** *Nat Rev Immunol.* 2013;13(1):34-45. [PMID 23222502](https://pubmed.ncbi.nlm.nih.gov/23222502/)

## 5. The acute-phase response — CLOCK 1

25. Gabay C, Kushner I. **Acute-phase proteins and other systemic responses to inflammation.** *N Engl J Med.* 1999;340(6):448-54. — fibrinogen is a positive, antithrombin and protein C are negative acute-phase proteins. The direct basis for CLOCK 1 in the model. [PMID 9971870](https://pubmed.ncbi.nlm.nih.gov/9971870/)
26. Fuller GM, Zhang Z. **Transcriptional control mechanism of fibrinogen gene expression.** *Ann N Y Acad Sci.* 2001;936:469-79. — a 2-5 fold increase in fibrinogen synthesis by IL-6. The basis for `APOSF = 3`. [PMID 11460504](https://pubmed.ncbi.nlm.nih.gov/11460504/)
27. Dahlbäck B. **C4b-binding protein: a forgotten factor in thrombosis and hemostasis.** *Semin Thromb Hemost.* 2011;37(4):355-61. — a rise in C4BP → a fall in free protein S. The model's `ANEGPS`. [PMID 21805442](https://pubmed.ncbi.nlm.nih.gov/21805442/)
28. Collen D, Tytgat GN, Claeys H, Piessens R. **Metabolism and distribution of fibrinogen. I. Fibrinogen turnover in physiological conditions in humans.** *Br J Haematol.* 1972;22(6):681-700. — a fibrinogen half-life of about 100 h and a daily turnover of 2-3 g. The basis for `THFIB` and for the scale of the consumption rate. [PMID 5064500](https://pubmed.ncbi.nlm.nih.gov/5064500/)
29. Collen D, Schetz J, de Cock F, Holmer E, Verstraete M. **Metabolism of antithrombin III (heparin cofactor) in man: effects of venous thrombosis and of heparin administration.** *Eur J Clin Invest.* 1977;7(1):27-35. — an AT half-life of about 2.8 days, and the observation that giving heparin increases AT turnover. The model's `THAT` and its coupled heparin-AT consumption term. [PMID 65284](https://pubmed.ncbi.nlm.nih.gov/65284/)

## 6. Fibrinolysis — CLOCK 2

30. Bachmann F. **Plasminogen-plasmin enzyme system.** In: *Hemostasis and Thrombosis.* — reference concentrations such as plasminogen 2 µM and α2-antiplasmin 1 µM.
31. Kruithof EK, Tran-Thang C, Gudinchet A, et al. **Fibrinolysis in pregnancy: a study of plasminogen activator inhibitors.** *Blood.* 1987;69(2):460-6. [PMID 3099859](https://pubmed.ncbi.nlm.nih.gov/3099859/)
32. Mesters RM, Flörke N, Ostermann H, Kienast J. **Increase of plasminogen activator inhibitor levels predicts outcome of leukocytopenic patients with sepsis.** *Thromb Haemost.* 1996;75(6):902-7. — PAI-1 of 200-1000 ng/mL in sepsis. The quantitative basis for CLOCK 2. [PMID 8822583](https://pubmed.ncbi.nlm.nih.gov/8822583/)
33. Hermans PW, Hibberd ML, Booy R, et al. **4G/5G promoter polymorphism in the plasminogen-activator-inhibitor-1 gene and outcome of meningococcal disease.** *Lancet.* 1999;354(9178):556-60. — PAI-1 genotype changes prognosis = suggesting that the block on fibrinolysis is causal. [PMID 10470699](https://pubmed.ncbi.nlm.nih.gov/10470699/)
34. Bajzar L, Morser J, Nesheim M. **TAFI, or plasma procarboxypeptidase B, couples the coagulation and fibrinolytic cascades through the thrombin-thrombomodulin complex.** *J Biol Chem.* 1996;271(28):16603-8. — the same thrombin-TM complex activates protein C and TAFI at once. This is where the double-edged nature of rTM in the model comes from. [PMID 8663147](https://pubmed.ncbi.nlm.nih.gov/8663147/)
35. Boffa MB, Koschinsky ML. **Curiouser and curiouser: recent advances in measurement of thrombin-activatable fibrinolysis inhibitor (TAFI) and in understanding its molecular genetics, gene regulation, and biological roles.** *Clin Biochem.* 2007;40(7-8):431-42. — an active TAFIa half-life of about 8-15 minutes. The model's `KDTAFI`. [PMID 17331491](https://pubmed.ncbi.nlm.nih.gov/17331491/)
36. Mosnier LO, Bouma BN. **Regulation of fibrinolysis by thrombin activatable fibrinolysis inhibitor, an unstable carboxypeptidase B that unites the pathways of coagulation and fibrinolysis.** *Arterioscler Thromb Vasc Biol.* 2006;26(11):2445-53. [PMID 16960106](https://pubmed.ncbi.nlm.nih.gov/16960106/)

## 7. The coagulopathy of acute promyelocytic leukaemia (APL)

37. Menell JS, Cesarman GM, Jacovina AT, McLaughlin MA, Lev EA, Hajjar KA. **Annexin II and bleeding in acute promyelocytic leukemia.** *N Engl J Med.* 1999;340(13):994-1004. — annexin A2 on the surface of APL blasts amplifies tPA-dependent plasmin generation. The model's `KANX`. [PMID 10099140](https://pubmed.ncbi.nlm.nih.gov/10099140/)
38. Falanga A, Marchetti M. **Anticancer treatment and thrombosis.** *Thromb Res.* 2012;129(3):353-9. [PMID 22119389](https://pubmed.ncbi.nlm.nih.gov/22119389/)
39. Stein E, McMahon B, Kwaan H, Altman JK, Frankfurt O, Tallman MS. **The coagulopathy of acute promyelocytic leukaemia revisited.** *Best Pract Res Clin Haematol.* 2009;22(1):153-63. — a synthesis of the pathophysiology showing APL to be of the hyperfibrinolytic type. [PMID 19285282](https://pubmed.ncbi.nlm.nih.gov/19285282/)
40. Sanz MA, Fenaux P, Tallman MS, et al. **Management of acute promyelocytic leukemia: updated recommendations from an expert panel of the European LeukemiaNet.** *Blood.* 2019;133(15):1630-43. — antifibrinolytics are not routinely recommended, and a thrombotic risk is flagged when they are combined with ATRA. Consistent with the result of scenario Q in the model. [PMID 30803991](https://pubmed.ncbi.nlm.nih.gov/30803991/)
41. Avvisati G, ten Cate JW, Büller HR, Mandelli F. **Tranexamic acid for control of haemorrhage in acute promyelocytic leukaemia.** *Lancet.* 1989;2(8655):122-4. — a randomised trial of tranexamic acid in APL: no benefit. [PMID 2567895](https://pubmed.ncbi.nlm.nih.gov/2567895/)
42. Lo-Coco F, Avvisati G, Vignetti M, et al. **Retinoic acid and arsenic trioxide for acute promyelocytic leukemia.** *N Engl J Med.* 2013;369(2):111-21. — APL0406. Two-year EFS 97% vs 86%. [PMID 23841729](https://pubmed.ncbi.nlm.nih.gov/23841729/)
43. Tallman MS, Andersen JW, Schiffer CA, et al. **All-trans-retinoic acid in acute promyelocytic leukemia.** *N Engl J Med.* 1997;337(15):1021-8. [PMID 9321529](https://pubmed.ncbi.nlm.nih.gov/9321529/)
44. Breccia M, Latagliata R, Cannella L, et al. **Early hemorrhagic death before starting therapy in acute promyelocytic leukemia: association with high WBC count, late diagnosis and delayed treatment initiation.** *Haematologica.* 2010;95(5):853-4. [PMID 20015874](https://pubmed.ncbi.nlm.nih.gov/20015874/)

## 8. Anticoagulant therapy trials — the arms the model must reproduce

45. Bernard GR, Vincent JL, Laterre PF, et al. **Efficacy and safety of recombinant human activated protein C for severe sepsis (PROWESS).** *N Engl J Med.* 2001;344(10):699-709. — 28-day mortality 24.7% vs 30.8%. Model scenario I. [PMID 11236773](https://pubmed.ncbi.nlm.nih.gov/11236773/)
46. Ranieri VM, Thompson BT, Barie PS, et al. **Drotrecogin alfa (activated) in adults with septic shock (PROWESS-SHOCK).** *N Engl J Med.* 2012;366(22):2055-64. — 26.4% vs 24.2%. Led to withdrawal of the drug. Model scenario J. [PMID 22616830](https://pubmed.ncbi.nlm.nih.gov/22616830/)
47. Warren BL, Eid A, Singer P, et al. **Caring for the critically ill patient. High-dose antithrombin III in severe sepsis: a randomized controlled trial (KyberSept).** *JAMA.* 2001;286(15):1869-78. — overall 38.9% vs 38.7%; benefit only in the subgroup not given heparin as well; bleeding 22.0% vs 12.8%. Model scenarios F/G. [PMID 11597289](https://pubmed.ncbi.nlm.nih.gov/11597289/)
48. Vincent JL, Francois B, Zabolotskikh I, et al. **Effect of a recombinant human soluble thrombomodulin on mortality in patients with sepsis-associated coagulopathy: the SCARLET randomized clinical trial.** *JAMA.* 2019;321(20):1993-2002. — 26.8% vs 29.4%, p=0.32. Model scenario H. [PMID 31104069](https://pubmed.ncbi.nlm.nih.gov/31104069/)
49. Vincent JL, Ramesh MK, Ernest D, et al. **A randomized, double-blind, placebo-controlled, phase 2b study to evaluate the safety and efficacy of recombinant human soluble thrombomodulin, ART-123, in patients with sepsis and suspected disseminated intravascular coagulation.** *Crit Care Med.* 2013;41(9):2069-79. [PMID 23979365](https://pubmed.ncbi.nlm.nih.gov/23979365/)
50. Saito H, Maruyama I, Shimazaki S, et al. **Efficacy and safety of recombinant human soluble thrombomodulin (ART-123) in disseminated intravascular coagulation: results of a phase III, randomized, double-blind clinical trial.** *J Thromb Haemost.* 2007;5(1):31-41. [PMID 17059423](https://pubmed.ncbi.nlm.nih.gov/17059423/)
51. Zarychanski R, Abou-Setta AM, Kanji S, et al. **The efficacy and safety of heparin in patients with sepsis: a systematic review and metaanalysis.** *Crit Care Med.* 2015;43(3):511-8. — the real effect size of heparin. The yardstick by which the model is judged to overestimate this value. [PMID 25493972](https://pubmed.ncbi.nlm.nih.gov/25493972/)
52. Umemura Y, Yamakawa K, Ogura H, Yuhara H, Fujimi S. **Efficacy and safety of anticoagulant therapy in three specific populations with sepsis: a meta-analysis of randomized controlled trials.** *J Thromb Haemost.* 2016;14(3):518-30. — the stratified result showing benefit only in the DIC subgroup. [PMID 26670422](https://pubmed.ncbi.nlm.nih.gov/26670422/)
53. Yamakawa K, Umemura Y, Hayakawa M, et al. **Benefit profile of anticoagulant therapy in sepsis: a nationwide multicentre registry in Japan.** *Crit Care.* 2016;20(1):229. [PMID 27472991](https://pubmed.ncbi.nlm.nih.gov/27472991/)
54. Jaimes F, De La Rosa G, Morales C, et al. **Unfractioned heparin for treatment of sepsis: a randomized clinical trial (The HETRASE Study).** *Crit Care Med.* 2009;37(4):1185-96. [PMID 19242322](https://pubmed.ncbi.nlm.nih.gov/19242322/)

## 9. Heparin pharmacology · antithrombin · heparin resistance

55. Hirsh J, Anand SS, Halperin JL, Fuster V. **Guide to anticoagulant therapy: heparin.** *Circulation.* 2001;103(24):2994-3018. — the standard account that heparin accelerates the AT-thrombin reaction some 1000-4000 fold, and that it uses AT catalytically. The basis for the model writing this term as Michaelis-Menten in AT. [PMID 11413093](https://pubmed.ncbi.nlm.nih.gov/11413093/)
56. Olson ST, Björk I, Sheffer R, Craig PA, Shore JD, Choay J. **Role of the antithrombin-binding pentasaccharide in heparin acceleration of antithrombin-proteinase reactions.** *J Biol Chem.* 1992;267(18):12528-38. — the Kd of AT-heparin binding and the catalytic stoichiometry. The basis for `KMAT`. [PMID 1618758](https://pubmed.ncbi.nlm.nih.gov/1618758/)
57. Levy JH, Connors JM. **Heparin resistance — clinical perspectives and management strategies.** *N Engl J Med.* 2021;385(9):826-32. — a synthesis showing that heparin resistance is not solely a matter of AT deficiency, and that binding to PF4, histones, and acute-phase proteins is involved. The model's `KHEPN` term. [PMID 34437782](https://pubmed.ncbi.nlm.nih.gov/34437782/)
58. Lane DA, Denton J, Flynn AM, Thunberg L, Lindahl U. **Anticoagulant activities of heparin oligosaccharides and their neutralization by platelet factor 4.** *Biochem J.* 1984;218(3):725-32. [PMID 6721826](https://pubmed.ncbi.nlm.nih.gov/6721826/)
59. Hirsh J, Bauer KA, Donati MB, Gould M, Samama MM, Weitz JI. **Parenteral anticoagulants: ACCP evidence-based clinical practice guidelines.** *Chest.* 2008;133(6 Suppl):141S-159S. — source of the PK parameters for UFH, enoxaparin, fondaparinux, and argatroban. [PMID 18574264](https://pubmed.ncbi.nlm.nih.gov/18574264/)

## 10. Antifibrinolytics

60. CRASH-2 trial collaborators. **Effects of tranexamic acid on death, vascular occlusive events, and blood transfusion in trauma patients with significant haemorrhage (CRASH-2).** *Lancet.* 2010;376(9734):23-32. — benefit when given early, harm when given after 3 hours. The same structure as the sign reversal of TXA in the model. [PMID 20554319](https://pubmed.ncbi.nlm.nih.gov/20554319/)
61. WOMAN Trial Collaborators. **Effect of early tranexamic acid administration on mortality, hysterectomy, and other morbidities in women with post-partum haemorrhage (WOMAN).** *Lancet.* 2017;389(10084):2105-16. [PMID 28456509](https://pubmed.ncbi.nlm.nih.gov/28456509/)
62. HALT-IT Trial Collaborators. **Effects of a high-dose 24-h infusion of tranexamic acid on death and thromboembolic events in patients with acute gastrointestinal bleeding (HALT-IT).** *Lancet.* 2020;395(10241):1927-36. — no benefit, an increase in venous thromboembolism. [PMID 32563378](https://pubmed.ncbi.nlm.nih.gov/32563378/)
63. Pabinger I, Fries D, Schöchl H, Streif W, Toller W. **Tranexamic acid for treatment and prophylaxis of bleeding and hyperfibrinolysis.** *Wien Klin Wochenschr.* 2017;129(9-10):303-16. — the PK and IC50 of TXA. [PMID 28251277](https://pubmed.ncbi.nlm.nih.gov/28251277/)

## 11. Transfusion · replacement therapy · guidelines

64. Wada H, Thachil J, Di Nisio M, et al. **Guidance for diagnosis and treatment of disseminated intravascular coagulation from harmonization of the recommendations from three guidelines.** *J Thromb Haemost.* 2013;11(4):761-7. [PMID 23379279](https://pubmed.ncbi.nlm.nih.gov/23379279/)
65. Levi M, Toh CH, Thachil J, Watson HG. **Guidelines for the diagnosis and management of disseminated intravascular coagulation.** *Br J Haematol.* 2009;145(1):24-33. [PMID 19222477](https://pubmed.ncbi.nlm.nih.gov/19222477/)
66. Squizzato A, Hunt BJ, Kinasewitz GT, et al. **Supportive management strategies for disseminated intravascular coagulation. An international consensus.** *Thromb Haemost.* 2016;115(5):896-904. [PMID 26791023](https://pubmed.ncbi.nlm.nih.gov/26791023/)
67. Kozek-Langenecker SA, Ahmed AB, Afshari A, et al. **Management of severe perioperative bleeding: guidelines from the European Society of Anaesthesiology. First update 2016.** *Eur J Anaesthesiol.* 2017;34(6):332-95. — a rise of about 25-30 mg/dL per 1 g of fibrinogen concentrate. The dose conversion for `FIBINF` in the model. [PMID 28459785](https://pubmed.ncbi.nlm.nih.gov/28459785/)
68. Estcourt LJ, Birchall J, Allard S, et al. **Guidelines for the use of platelet transfusions.** *Br J Haematol.* 2017;176(3):365-94. — a rise of about 30-50 x10^9/L per adult apheresis platelet unit. [PMID 28009056](https://pubmed.ncbi.nlm.nih.gov/28009056/)
69. Levy JH, Sniecinski RM, Welsby IJ, Levi M. **Antithrombin: anti-inflammatory properties and clinical applications.** *Thromb Haemost.* 2016;115(4):712-28. — a rise of about 1.4% per 1 IU/kg of AT concentrate. The `ATC_LOAD` conversion in the model. [PMID 26676886](https://pubmed.ncbi.nlm.nih.gov/26676886/)
70. Iba T, Levy JH, Raj A, Warkentin TE. **Advance in the management of sepsis-induced coagulopathy and disseminated intravascular coagulation.** *J Clin Med.* 2019;8(5):728. [PMID 31121897](https://pubmed.ncbi.nlm.nih.gov/31121897/)

## 12. Obstetric, trauma and malignancy-associated DIC

71. Erez O, Mastrolia SA, Thachil J. **Disseminated intravascular coagulation in pregnancy: insights in pathophysiology, diagnosis and management.** *Am J Obstet Gynecol.* 2015;213(4):452-63. [PMID 25840271](https://pubmed.ncbi.nlm.nih.gov/25840271/)
72. Brohi K, Cohen MJ, Ganter MT, et al. **Acute coagulopathy of trauma: hypoperfusion induces systemic anticoagulation and hyperfibrinolysis.** *Ann Surg.* 2007;245(5):812-8. — the basis for traumatic coagulopathy being distinguished from DIC by way of the activated protein C pathway. [PMID 17457176](https://pubmed.ncbi.nlm.nih.gov/17457176/)
73. Moore HB, Moore EE, Gonzalez E, et al. **Hyperfibrinolysis, physiologic fibrinolysis, and fibrinolysis shutdown: the spectrum of postinjury fibrinolysis and relevance to antifibrinolytic therapy.** *J Trauma Acute Care Surg.* 2014;77(6):811-7. — quantification of the "fibrinolysis shutdown" concept. The clinical counterpart of CLOCK 2 in the model. [PMID 25051384](https://pubmed.ncbi.nlm.nih.gov/25051384/)
74. Levi M. **Disseminated intravascular coagulation in cancer patients.** *Best Pract Res Clin Haematol.* 2009;22(1):129-36. [PMID 19285280](https://pubmed.ncbi.nlm.nih.gov/19285280/)

## 13. QSP modelling methodology — systems modelling of coagulation

75. Hockin MF, Jones KC, Everse SJ, Mann KG. **A model for the stoichiometric regulation of blood coagulation.** *J Biol Chem.* 2002;277(21):18322-33. — the classic ODE model of the coagulation cascade. [PMID 11893748](https://pubmed.ncbi.nlm.nih.gov/11893748/)
76. Chelle P, Morin C, Montmartin A, Piot M, Cournil M, Tardy-Poncet B. **Evaluation and calibration of a mechanistic evolutionary model of blood coagulation.** *PLoS One.* 2018;13(3):e0193658. [PMID 29505602](https://pubmed.ncbi.nlm.nih.gov/29505602/)
77. Wajima T, Isbister GK, Duffull SB. **A comprehensive model for the humoral coagulation network in humans.** *Clin Pharmacol Ther.* 2009;86(3):290-8. — the standard reference for a coagulation-network model that can be coupled to clinical PK/PD. Many of the factor half-lives and baseline concentrations in this model agree with it. [PMID 19516255](https://pubmed.ncbi.nlm.nih.gov/19516255/)
78. Nayak S, Lee D, Patel-Hett S, et al. **Using a systems pharmacology model of the blood coagulation network to predict the effects of various therapies on biomarkers.** *CPT Pharmacometrics Syst Pharmacol.* 2015;4(7):396-405. [PMID 26312163](https://pubmed.ncbi.nlm.nih.gov/26312163/)
79. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** *J Pharmacokinet Pharmacodyn.* 2015;42:S84-5. — the simulation engine this repository uses.
80. Gulati A, Faed JM, Isbister GK, Duffull SB. **Application of a systems model of coagulation to the prediction of the effect of warfarin.** *CPT Pharmacometrics Syst Pharmacol.* 2014;3:e93. [PMID 24429593](https://pubmed.ncbi.nlm.nih.gov/24429593/)

---

## 14. Equation ↔ source map

A summary of which reference each structural choice in the model rests on.

| Model element | Supporting reference |
|---|---|
| ISTH / SIC / JAAM scores computed as **outputs** of the state vector | 1, 4, 6 |
| CLOCK 1 — fibrinogen is a positive, AT, PC and PS are negative acute-phase proteins | 25, 26, 27 |
| Fibrinogen half-life 100 h, daily turnover 2-3 g | 28 |
| AT half-life 65 h, heparin accelerates AT consumption | 29 |
| Protein C half-life 6 h (the anticoagulant factor that hits the floor first) | 14, 77 |
| CLOCK 2 — PAI-1 200-1000 ng/mL, t-PA 20-30 ng/mL | 32, 33 |
| The PAI-1 : t-PA molar ratio causally changes prognosis | 33 |
| annexin A2 → t-PA-dependent plasmin amplification in APL | 37 |
| α2-antiplasmin depletion is the gateway to systemic fibrinogenolysis | 30, 39 |
| APC generation = THR x TM x EPCR x PC (a product of four terms) | 14, 15, 23 |
| The same thrombin-TM complex also activates TAFI | 34, 36 |
| Active TAFIa half-life about 10 minutes | 35 |
| ART-123 activates protein C well but TAFI poorly | 34, 49, 50 |
| Heparin is a **catalyst** with respect to AT, so the AT term saturates | 55, 56 |
| Non-AT heparin resistance caused by PF4 and histones | 57, 58 |
| Histones are the principal driving term for platelet activation | 20, 22 |
| Loss of the glycocalyx → capillary leak of AT/PC | 17, 18 |
| Scenarios F/G (KyberSept) | 47 |
| Scenario H (SCARLET) | 48, 49, 50 |
| Scenarios I/J (PROWESS, PROWESS-SHOCK) | 45, 46 |
| Scenarios D/E, heparin (where the model **overestimates**) | 51, 52, 54 |
| Scenarios N/O/Q (APL, ATRA, antifibrinolytics) | 40, 41, 42, 43 |
| The sign reversal of TXA (early vs late, bleeding type vs shutdown type) | 60, 62, 73 |
| Blood-product dose → concentration conversion | 67, 68, 69 |
| Prior structures for ODE models of the coagulation network | 75, 76, 77, 78, 80 |

---

## 15. Where this model **disagrees** with the evidence

Recorded here honestly. Two of the predictions the model makes do not agree with the current
clinical evidence, and should be read as falsifiable hypotheses.

1. **The effect size of unfractionated heparin.** The model lowers 28-day mortality in septic DIC from 43.3% to
   28.3% (−15.0%p). In the actual meta-analyses (51, 52, 54) it is at best a few
   %p, and even that is confined to the DIC subgroup. The mechanisms most likely to be
   missing from the model are (i) heparin-induced bleeding into organs that are already damaged, and (ii) the fact
   that aPTT/anti-Xa monitoring is difficult in DIC, so that the target concentration is not in fact
   reached.
2. **AT concentrate combined with heparin.** The model predicts that the combination is better than AT alone, but
   in KyberSept (47) the benefit disappeared in the combination arm and bleeding increased. The model
   reproduces the direction (the additional benefit of AT shrinks and the bleeding burden grows when the two are combined)
   but not the magnitude.

Conversely, the point at which the model **falsified the author's prior hypothesis** is also recorded.
The explanation that "deposited fibrin is a stock and an anticoagulant touches only the rate, so giving it
late is useless" failed to account for PROWESS-SHOCK in this model.
Sweeping the start time of drotrecogin alfa across 0/6/12/24/30/48/72/96 hours gives
28-day mortality of 37.8 / 37.4 / 36.9 / 36.1 / 35.8 / 35.9 / 37.3 / 39.1 %,
so the benefit is **maximal at 24-48 hours** and only collapses beyond 72 hours.
Timing matters, but not in the direction of "too late".
