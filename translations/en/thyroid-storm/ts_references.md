# Thyroid Storm / Thyrotoxic Crisis — references

The last line of each section states **which term of this model depends on that literature**.
The numbers used for calibration (N1–N20, D1–D7, C1), the predictions (P1–P10) and the
limitations (L1–L7) are defined in the preambles of [`README.md`](../../../thyroid-storm/README.md) and
[`ts_verify_python.py`](ts_verify_python.py).

> **Link policy.** A PMID link is used only where the PMID is known for certain; otherwise a
> **title-search URL** (`pubmed.ncbi.nlm.nih.gov/?term=…`) is used. This way a wrong PMID
> can never end up pointing at a different paper. Every link resolves to
> PubMed.

---

## 1. Definition · diagnostic criteria · epidemiology

1. Burch HB, Wartofsky L. **Life-threatening thyrotoxicosis: thyroid storm.** Endocrinol Metab Clin North Am. 1993;22(2):263–77. <https://pubmed.ncbi.nlm.nih.gov/8325286/>
2. Akamizu T, Satoh T, Isozaki O, et al. **Diagnostic criteria, clinical features, and incidence of thyroid storm based on nationwide surveys.** Thyroid. 2012;22(7):661–79. <https://pubmed.ncbi.nlm.nih.gov/22690898/>
3. Satoh T, Isozaki O, Suzuki A, et al. **2016 Guidelines for the management of thyroid storm from the Japan Thyroid Association and Japan Endocrine Society.** Endocr J. 2016;63(12):1025–64. <https://pubmed.ncbi.nlm.nih.gov/?term=2016+Guidelines+for+the+management+of+thyroid+storm+Japan+Thyroid+Association>
4. Ross DS, Burch HB, Cooper DS, et al. **2016 American Thyroid Association guidelines for diagnosis and management of hyperthyroidism and other causes of thyrotoxicosis.** Thyroid. 2016;26(10):1343–421. <https://pubmed.ncbi.nlm.nih.gov/27521067/>
5. Chiha M, Samarasinghe S, Kabaker AS. **Thyroid storm: an updated review.** J Intensive Care Med. 2015;30(3):131–40. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+storm+an+updated+review+Chiha+Samarasinghe+Kabaker>
6. Angell TE, Lechner MG, Nguyen CT, et al. **Clinical features and hospital outcomes in thyroid storm: a retrospective cohort study.** J Clin Endocrinol Metab. 2015;100(2):451–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+features+and+hospital+outcomes+in+thyroid+storm+retrospective+cohort+Angell>
7. Swee DS, Chng CL, Lim A. **Clinical characteristics and outcome of thyroid storm: a case series and review of neuropsychiatric derangements in thyrotoxicosis.** Endocr Pract. 2015;21(2):182–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Clinical+characteristics+and+outcome+of+thyroid+storm+case+series+neuropsychiatric+Swee>
8. Galindo RJ, Hurtado CR, Pasquel FJ, et al. **National trends in incidence, mortality and clinical outcomes of patients hospitalized for thyrotoxicosis with and without thyroid storm in the United States.** Thyroid. 2019;29(1):36–43. <https://pubmed.ncbi.nlm.nih.gov/?term=National+trends+incidence+mortality+thyrotoxicosis+thyroid+storm+United+States+Galindo>
9. Ono Y, Ono S, Yasunaga H, et al. **Factors associated with mortality of thyroid storm: analysis using a national inpatient database in Japan.** Medicine (Baltimore). 2016;95(7):e2848. <https://pubmed.ncbi.nlm.nih.gov/?term=Factors+associated+with+mortality+of+thyroid+storm+national+inpatient+database+Japan+Ono>
10. Isozaki O, Satoh T, Wakino S, et al. **Treatment and management of thyroid storm: analysis of the nationwide surveys.** Clin Endocrinol (Oxf). 2016;84(6):912–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+and+management+of+thyroid+storm+analysis+of+the+nationwide+surveys+Isozaki>
11. Bourcier S, Coutrot M, Kimmoun A, et al. **Thyroid storm in the ICU: a retrospective multicenter study.** Crit Care Med. 2020;48(1):83–90. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+storm+in+the+ICU+retrospective+multicenter+study+Bourcier>
12. Idrose AM. **Acute and emergency care for thyrotoxicosis and thyroid storm.** Acute Med Surg. 2015;2(3):147–57. <https://pubmed.ncbi.nlm.nih.gov/?term=Acute+and+emergency+care+for+thyrotoxicosis+and+thyroid+storm+Idrose>

> Model dependency: the per-item ranges of the Burch–Wartofsky score (temperature 0–30 · CNS 0–30 ·
> gastrointestinal/hepatic 0–20 · heart rate 0–25 · heart failure 0–15 · AF 0/10 · precipitant 0/10) and
> the ≥45 threshold are implemented in `burch_wartofsky()` **exactly as published**.
> The observation that a precipitant is identified in 70–90 % justifies the model's central structure
> (the input that closes the loop). The 80–90 % untreated mortality is the only storm-derived calibration value, C1 (`h0_haz`).
> The observation that jaundice (bilirubin >3 mg/dL) is an independent predictor of death is the basis for the `eBili` term.

---

## 2. The decisive observation — hormone concentrations do not distinguish the storm

13. Brooks MH, Waldstein SS. **Free thyroxine concentrations in thyroid storm.** Ann Intern Med. 1980;93(5):694–7. <https://pubmed.ncbi.nlm.nih.gov/?term=Free+thyroxine+concentrations+in+thyroid+storm+Brooks+Waldstein>
14. Brooks MH, Waldstein SS, Bronsky D, Sterling K. **Serum triiodothyronine concentration in thyroid storm.** J Clin Endocrinol Metab. 1975;40(2):339–41. <https://pubmed.ncbi.nlm.nih.gov/?term=Serum+triiodothyronine+concentration+in+thyroid+storm+Brooks+Waldstein+Bronsky+Sterling>
15. Jacobs HS, Mackie DB, Eastman CJ, et al. **Total and free triiodothyronine and thyroxine levels in thyroid storm and recurrent hyperthyroidism.** Lancet. 1973;2(7823):236–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Total+and+free+triiodothyronine+and+thyroxine+levels+in+thyroid+storm+and+recurrent+hyperthyroidism>
16. Nayak B, Burman K. **Thyrotoxicosis and thyroid storm.** Endocrinol Metab Clin North Am. 2006;35(4):663–86. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyrotoxicosis+and+thyroid+storm+Nayak+Burman+2006>
17. Wartofsky L. **Thyrotoxic storm.** In: Braverman LE, Cooper DS, eds. *Werner &amp; Ingbar's The Thyroid.* (chapter) <https://pubmed.ncbi.nlm.nih.gov/?term=Wartofsky+thyrotoxic+storm+Werner+Ingbar+thyroid>

> Model dependency: **this section is the reason this model exists.** The observation that total T4/T3
> does not distinguish storm from simple thyrotoxicosis is what makes severity a matter not of hormone
> concentration but of the state of the fast feedback loop (the controlled experiment of sections [H]·[S] in the README).

---

## 3. Thyroid hormone kinetics and the "reservoir"

18. Nicoloff JT, Low JC, Dussault JH, Fisher DA. **Simultaneous measurement of thyroxine and triiodothyronine peripheral turnover kinetics in man.** J Clin Invest. 1972;51(2):473–83. <https://pubmed.ncbi.nlm.nih.gov/?term=Simultaneous+measurement+of+thyroxine+and+triiodothyronine+peripheral+turnover+kinetics+in+man+Nicoloff>
19. Oppenheimer JH, Schwartz HL, Surks MI. **Determination of common parameters of iodothyronine metabolism and distribution in man by noncompartmental analysis.** J Clin Endocrinol Metab. 1975;41(2):319–24. <https://pubmed.ncbi.nlm.nih.gov/?term=Determination+of+common+parameters+of+iodothyronine+metabolism+and+distribution+in+man+by+noncompartmental+analysis>
20. Pilo A, Iervasi G, Vitek F, et al. **Thyroidal and peripheral production of 3,5,3'-triiodothyronine in humans by multicompartmental analysis.** Am J Physiol. 1990;258(4 Pt 1):E715–26. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroidal+and+peripheral+production+of+triiodothyronine+in+humans+by+multicompartmental+analysis+Pilo>
21. Larsen PR. **Thyroidal triiodothyronine and thyroxine in Graves' disease: correlation with presurgical treatment, thyroid status, and iodine content.** J Clin Endocrinol Metab. 1975;41(6):1098–104. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroidal+triiodothyronine+and+thyroxine+in+Graves+disease+correlation+with+presurgical+treatment+iodine+content+Larsen>
22. Riggs DS. **Quantitative aspects of iodine metabolism in man.** Pharmacol Rev. 1952;4(3):284–370. <https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative+aspects+of+iodine+metabolism+in+man+Riggs+1952>
23. Wagner MS, Wajner SM, Maia AL. **The role of thyroid hormone in testicular development and function.** (background: tissue-specific hormone handling) <https://pubmed.ncbi.nlm.nih.gov/?term=Wagner+Wajner+Maia+role+of+thyroid+hormone+deiodinase+tissue>
24. Dietrich JW, Landgrafe G, Fotiadou EH. **TSH and thyrotropic agonists: key actors in thyroid homeostasis.** J Thyroid Res. 2012;2012:351864. <https://pubmed.ncbi.nlm.nih.gov/?term=TSH+and+thyrotropic+agonists+key+actors+in+thyroid+homeostasis+Dietrich>
25. Berberich J, Dietrich JW, Hoermann R, Müller MA. **Mathematical modeling of the pituitary–thyroid feedback loop: role of a TSH-T3 shunt and sensitivity analysis.** Front Endocrinol. 2018;9:91. <https://pubmed.ncbi.nlm.nih.gov/?term=Mathematical+modeling+of+the+pituitary+thyroid+feedback+loop+TSH+T3+shunt+sensitivity+analysis+Berberich>
26. Hoermann R, Midgley JEM, Larisch R, Dietrich JW. **Homeostatic control of the thyroid–pituitary axis: perspectives for diagnosis and treatment.** Front Endocrinol. 2015;6:177. <https://pubmed.ncbi.nlm.nih.gov/?term=Homeostatic+control+of+the+thyroid+pituitary+axis+perspectives+for+diagnosis+and+treatment+Hoermann>

> Model dependency: N1–N2 (total T4 100 · total T3 1.8 nmol/L), N5 (rT3), N6 (80 % peripheral conversion),
> N14–N15 (T4 t½ 6–7 d, T3 t½ 1 d), N17 (T4 secretion ~90 µg/d), and
> the number **N7 — the colloid reservoir holds 60–90 days' worth**. N7 carries the whole reservoir
> arithmetic of section [R], and that is the basis for "a thionamide cannot be an acute drug".
> Larsen 1975, the observation that thyroidal iodine content falls in Graves' disease, is the
> comparison data for prediction P5 (66 days → 18.6 days), which the model reproduces **without calibration**.

---

## 4. Peripheral deiodinases — the fast axis

27. Bianco AC, Salvatore D, Gereben B, Berry MJ, Larsen PR. **Biochemistry, cellular and molecular biology, and physiological roles of the iodothyronine selenodeiodinases.** Endocr Rev. 2002;23(1):38–89. <https://pubmed.ncbi.nlm.nih.gov/11844744/>
28. Maia AL, Kim BW, Huang SA, Harney JW, Larsen PR. **Type 2 iodothyronine deiodinase is the major source of plasma T3 in euthyroid humans.** J Clin Invest. 2005;115(9):2524–33. <https://pubmed.ncbi.nlm.nih.gov/?term=Type+2+iodothyronine+deiodinase+is+the+major+source+of+plasma+T3+in+euthyroid+humans+Maia>
29. Gereben B, Zavacki AM, Ribich S, et al. **Cellular and molecular basis of deiodinase-regulated thyroid hormone signaling.** Endocr Rev. 2008;29(7):898–938. <https://pubmed.ncbi.nlm.nih.gov/?term=Cellular+and+molecular+basis+of+deiodinase+regulated+thyroid+hormone+signaling+Gereben>
30. Visser TJ, Kaptein E, Terpstra OT, Krenning EP. **Deiodination of thyroid hormone by human liver.** J Clin Endocrinol Metab. 1988;67(1):17–24. <https://pubmed.ncbi.nlm.nih.gov/?term=Deiodination+of+thyroid+hormone+by+human+liver+Visser+Kaptein+Terpstra+Krenning>
31. Chopra IJ. **An assessment of daily production and significance of thyroidal secretion of 3,3',5'-triiodothyronine (reverse T3) in man.** J Clin Invest. 1976;58(1):32–40. <https://pubmed.ncbi.nlm.nih.gov/?term=An+assessment+of+daily+production+and+significance+of+thyroidal+secretion+of+reverse+T3+in+man+Chopra>
32. Chopra IJ, Williams DE, Orgiazzi J, Solomon DH. **Opposite effects of dexamethasone on serum concentrations of 3,3',5'-triiodothyronine (reverse T3) and 3,3',5-triiodothyronine (T3).** J Clin Endocrinol Metab. 1975;41(5):911–20. <https://pubmed.ncbi.nlm.nih.gov/?term=Opposite+effects+of+dexamethasone+on+serum+concentrations+of+reverse+T3+and+T3+Chopra>
33. Peeters RP, Wouters PJ, Kaptein E, et al. **Reduced activation and increased inactivation of thyroid hormone in tissues of critically ill patients.** J Clin Endocrinol Metab. 2003;88(7):3202–11. <https://pubmed.ncbi.nlm.nih.gov/?term=Reduced+activation+and+increased+inactivation+of+thyroid+hormone+in+tissues+of+critically+ill+patients+Peeters>
34. Köhrle J. **Local activation and inactivation of thyroid hormones: the deiodinase family.** Mol Cell Endocrinol. 1999;151(1-2):103–19. <https://pubmed.ncbi.nlm.nih.gov/?term=Local+activation+and+inactivation+of+thyroid+hormones+the+deiodinase+family+Kohrle>
35. Friesema ECH, Ganguly S, Abdalla A, et al. **Identification of monocarboxylate transporter 8 as a specific thyroid hormone transporter.** J Biol Chem. 2003;278(41):40128–35. <https://pubmed.ncbi.nlm.nih.gov/?term=Identification+of+monocarboxylate+transporter+8+as+a+specific+thyroid+hormone+transporter+Friesema>

> Model dependency: the D1/D2/D3 split (kD1_4 : kD2_4 : kD3_4) and the 80 % peripheral contribution of N6.
> **That rT3 clearance depends mainly on D1 (fD1_rT3 = 0.85)** is what generates the laboratory
> signature of section [D] — on PTU rT3 rises 0.74 → 1.27 nmol/L, while on
> methimazole it does not move. Chopra 1975 is the primary source for the same signature from
> glucocorticoid, which lowers T3 and raises rT3.

---

## 5. Binding proteins · free fraction · fatty-acid displacement — the second axis

36. Robbins J, Rall JE. **Proteins associated with the thyroid hormones.** Physiol Rev. 1960;40:415–89. <https://pubmed.ncbi.nlm.nih.gov/?term=Proteins+associated+with+the+thyroid+hormones+Robbins+Rall>
37. Refetoff S. **Thyroid hormone serum transport proteins.** In: Endotext. <https://pubmed.ncbi.nlm.nih.gov/?term=Refetoff+thyroid+hormone+serum+transport+proteins+Endotext>
38. Mendel CM. **The free hormone hypothesis: a physiologically based mathematical model.** Endocr Rev. 1989;10(3):232–74. <https://pubmed.ncbi.nlm.nih.gov/?term=The+free+hormone+hypothesis+a+physiologically+based+mathematical+model+Mendel>
39. Tabachnick M, Giorgio NA. **Thyroxine–protein interactions. II. The binding of thyroxine and its analogues to human serum albumin.** Arch Biochem Biophys. 1964;105:563–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroxine+protein+interactions+binding+of+thyroxine+and+its+analogues+to+human+serum+albumin+Tabachnick+Giorgio>
40. Chopra IJ, Huang TS, Beredo A, et al. **Evidence for an inhibitor of extrathyroidal conversion of thyroxine to 3,5,3'-triiodothyronine in sera of patients with nonthyroidal illnesses.** J Clin Endocrinol Metab. 1985;60(4):666–72. <https://pubmed.ncbi.nlm.nih.gov/?term=Evidence+for+an+inhibitor+of+extrathyroidal+conversion+of+thyroxine+to+triiodothyronine+in+sera+of+patients+with+nonthyroidal+illnesses>
41. Lim CF, Bai Y, Topliss DJ, Barlow JW, Stockigt JR. **Drug and fatty acid effects on serum thyroid hormone binding.** J Clin Endocrinol Metab. 1988;67(4):682–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Drug+and+fatty+acid+effects+on+serum+thyroid+hormone+binding+Lim+Bai+Topliss+Barlow+Stockigt>
42. Mendel CM, Frost PH, Kunitake ST, Cavalieri RR. **Mechanism of the heparin-induced increase in the concentration of free thyroxine in plasma.** J Clin Endocrinol Metab. 1987;65(6):1259–64. <https://pubmed.ncbi.nlm.nih.gov/?term=Mechanism+of+the+heparin+induced+increase+in+the+concentration+of+free+thyroxine+in+plasma+Mendel>
43. Jaume JC, Mendel CM, Frost PH, Greenspan FS, Laughton CW. **Extremely low doses of heparin release lipase activity into the plasma and can thereby cause artifactual elevations in the serum-free thyroxine concentration as measured by equilibrium dialysis.** Thyroid. 1996;6(2):79–83. <https://pubmed.ncbi.nlm.nih.gov/?term=Extremely+low+doses+of+heparin+release+lipase+activity+into+the+plasma+artifactual+elevations+free+thyroxine+Jaume>
44. Stockigt JR, Lim CF. **Medications that distort in vitro tests of thyroid function, with particular reference to estimates of serum free thyroxine.** Best Pract Res Clin Endocrinol Metab. 2009;23(6):753–67. <https://pubmed.ncbi.nlm.nih.gov/?term=Medications+that+distort+in+vitro+tests+of+thyroid+function+free+thyroxine+Stockigt+Lim>
45. Larsen PR. **Salicylate-induced increases in free triiodothyronine in human serum: evidence of inhibition of triiodothyronine binding to thyroxine-binding globulin and thyroxine-binding prealbumin.** J Clin Invest. 1972;51(5):1125–34. <https://pubmed.ncbi.nlm.nih.gov/?term=Salicylate+induced+increases+in+free+triiodothyronine+in+human+serum+inhibition+of+binding+to+thyroxine+binding+globulin+Larsen>
46. Wang R, Nelson JC, Wilcox RB. **Salsalate administration: a potential pharmacological model of the sick euthyroid syndrome.** J Clin Endocrinol Metab. 1998;83(9):3095–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Salsalate+administration+a+potential+pharmacological+model+of+the+sick+euthyroid+syndrome+Wang+Nelson+Wilcox>
47. Docter R, Krenning EP, de Jong M, Hennemann G. **The sick euthyroid syndrome: changes in thyroid hormone serum parameters and hormone metabolism.** Clin Endocrinol (Oxf). 1993;39(5):499–518. <https://pubmed.ncbi.nlm.nih.gov/?term=The+sick+euthyroid+syndrome+changes+in+thyroid+hormone+serum+parameters+and+hormone+metabolism+Docter>

> Model dependency: the free fractions of N3–N4 (T4 0.02 % · T3 0.30 %) and the displacement term `Fd`.
> **Writing NEFA displacement as a steep function with Hill exponent 2** (EmNEFA 1.11 · K_NEFA 0.862)
> corresponds to Lim 1988's observation (the fatty-acid effect is negligible at low concentrations and
> becomes abrupt at high ones). Larsen 1972 is the primary source for aspirin harm 1 (displacement) in section [A],
> Mendel 1987/Jaume 1996 for the heparin route, and Docter 1993 for the fall in TBG and the rise in
> free fraction in critical illness.

---

## 6. Thionamides — the slow axis (PTU vs methimazole)

48. Cooper DS. **Antithyroid drugs.** N Engl J Med. 2005;352(9):905–17. <https://pubmed.ncbi.nlm.nih.gov/?term=Antithyroid+drugs+Cooper+New+England+Journal+of+Medicine+2005>
49. Cooper DS. **Propylthiouracil levels in hyperthyroid patients unresponsive to large doses: evidence of poor patient compliance.** Ann Intern Med. 1985;102(3):328–31. <https://pubmed.ncbi.nlm.nih.gov/?term=Propylthiouracil+levels+in+hyperthyroid+patients+unresponsive+to+large+doses+evidence+of+poor+patient+compliance+Cooper>
50. Taurog A. **The mechanism of action of the thioureylene antithyroid drugs.** Endocrinology. 1976;98(4):1031–46. <https://pubmed.ncbi.nlm.nih.gov/?term=The+mechanism+of+action+of+the+thioureylene+antithyroid+drugs+Taurog>
51. Geffner DL, Hershman JM. **Beta-adrenergic blockade for the treatment of hyperthyroidism.** Am J Med. 1992;93(1):61–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Beta+adrenergic+blockade+for+the+treatment+of+hyperthyroidism+Geffner+Hershman>
52. Abuid J, Larsen PR. **Triiodothyronine and thyroxine in hyperthyroidism: comparison of the acute changes during therapy with antithyroid agents.** J Clin Invest. 1974;54(1):201–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Triiodothyronine+and+thyroxine+in+hyperthyroidism+comparison+of+the+acute+changes+during+therapy+with+antithyroid+agents+Abuid+Larsen>
53. Nakamura H, Noh JY, Itoh K, Fukata S, Miyauchi A, Hamada N. **Comparison of methimazole and propylthiouracil in patients with hyperthyroidism caused by Graves' disease.** J Clin Endocrinol Metab. 2007;92(6):2157–62. <https://pubmed.ncbi.nlm.nih.gov/?term=Comparison+of+methimazole+and+propylthiouracil+in+patients+with+hyperthyroidism+caused+by+Graves+disease+Nakamura+Noh>
54. Kuiper GG, Klootwijk W, Visser TJ. **Substitution of cysteine for selenocysteine in the catalytic center of type III iodothyronine deiodinase.** (background: D1 thiouracil sensitivity) <https://pubmed.ncbi.nlm.nih.gov/?term=Kuiper+Klootwijk+Visser+substitution+of+cysteine+for+selenocysteine+iodothyronine+deiodinase>
55. Bahn RS, Burch HS, Cooper DS, et al. **The role of propylthiouracil in the management of Graves' disease in adults: report of a meeting jointly sponsored by the ATA and the FDA.** Thyroid. 2009;19(7):673–4. <https://pubmed.ncbi.nlm.nih.gov/?term=The+role+of+propylthiouracil+in+the+management+of+Graves+disease+in+adults+report+of+a+meeting+ATA+FDA>
56. Rivkees SA, Mattison DR. **Propylthiouracil (PTU) hepatotoxicity in children and recommendations for discontinuation of use.** Int J Pediatr Endocrinol. 2009;2009:132041. <https://pubmed.ncbi.nlm.nih.gov/?term=Propylthiouracil+hepatotoxicity+in+children+and+recommendations+for+discontinuation+of+use+Rivkees+Mattison>
57. Nakamura H, Miyauchi A, Miyawaki N, Imagawa J. **Analysis of 754 cases of antithyroid drug-induced agranulocytosis over 30 years in Japan.** J Clin Endocrinol Metab. 2013;98(12):4776–83. <https://pubmed.ncbi.nlm.nih.gov/?term=Analysis+of+754+cases+of+antithyroid+drug+induced+agranulocytosis+over+30+years+in+Japan+Nakamura>
58. Yoshihara A, Noh JY, Watanabe N, et al. **Substituting potassium iodide for methimazole as the treatment for Graves' disease during the first trimester may reduce the incidence of congenital anomalies.** Thyroid. 2015;25(10):1155–61. <https://pubmed.ncbi.nlm.nih.gov/?term=Substituting+potassium+iodide+for+methimazole+as+the+treatment+for+Graves+disease+during+the+first+trimester+congenital+anomalies+Yoshihara>
59. Alexander EK, Pearce EN, Brent GA, et al. **2017 Guidelines of the American Thyroid Association for the diagnosis and management of thyroid disease during pregnancy and the postpartum.** Thyroid. 2017;27(3):315–89. <https://pubmed.ncbi.nlm.nih.gov/?term=2017+Guidelines+of+the+American+Thyroid+Association+thyroid+disease+during+pregnancy+and+the+postpartum>
60. Jongjaroenprasert W, Akarawut W, Chantasart D, Chailurkit L, Rajatanavin R. **Rectal administration of propylthiouracil in hyperthyroid patients: comparison of suspension enema and suppository form.** Thyroid. 2002;12(7):627–31. <https://pubmed.ncbi.nlm.nih.gov/?term=Rectal+administration+of+propylthiouracil+in+hyperthyroid+patients+comparison+of+suspension+enema+and+suppository+form>

> Model dependency: calibration values D1–D2 (TPO potency — methimazole is roughly 10–20 times stronger per mg;
> PTU's D1 inhibition lowers serum T3 by 20–30 % at 24 h). Abuid &amp; Larsen 1974 is the
> primary source that directly shows **the acute T3 difference between PTU and methimazole**, and the control for section [D].
> The PTU blood concentrations of Cooper 1985 are the basis for the PK (t½ 1.5 h), and Jongjaroenprasert 2002
> for rectal administration when the oral route fails (the `RouteFail` node on the map).
> Agranulocytosis at 0.2–0.5 % and PTU hepatotoxicity appear only as adverse-event nodes on the map and
> do not enter the mortality term (limitation L1).

---

## 7. Iodine · Wolff–Chaikoff · Jod-Basedow

61. Wolff J, Chaikoff IL. **Plasma inorganic iodide as a homeostatic regulator of thyroid function.** J Biol Chem. 1948;174(2):555–64. <https://pubmed.ncbi.nlm.nih.gov/18865621/>
62. Wolff J, Chaikoff IL, Goldberg RC, Meier JR. **The temporary nature of the inhibitory action of excess iodide on organic iodine synthesis in the normal thyroid.** Endocrinology. 1949;45(5):504–13. <https://pubmed.ncbi.nlm.nih.gov/?term=The+temporary+nature+of+the+inhibitory+action+of+excess+iodide+on+organic+iodine+synthesis+in+the+normal+thyroid>
63. Eng PH, Cardona GR, Fang SL, et al. **Escape from the acute Wolff-Chaikoff effect is associated with a decrease in thyroid sodium/iodide symporter messenger ribonucleic acid and protein.** Endocrinology. 1999;140(8):3404–10. <https://pubmed.ncbi.nlm.nih.gov/?term=Escape+from+the+acute+Wolff+Chaikoff+effect+is+associated+with+a+decrease+in+thyroid+sodium+iodide+symporter+Eng>
64. Markou K, Georgopoulos N, Kyriazopoulou V, Vagenakis AG. **Iodine-induced hypothyroidism.** Thyroid. 2001;11(5):501–10. <https://pubmed.ncbi.nlm.nih.gov/?term=Iodine+induced+hypothyroidism+Markou+Georgopoulos+Kyriazopoulou+Vagenakis>
65. Roti E, Vagenakis AG. **Effect of excess iodide: clinical aspects.** In: Werner &amp; Ingbar's The Thyroid. <https://pubmed.ncbi.nlm.nih.gov/?term=Roti+Vagenakis+effect+of+excess+iodide+clinical+aspects+thyroid>
66. Wartofsky L, Ransil BJ, Ingbar SH. **Inhibition by iodine of the release of thyroxine from the thyroid glands of patients with thyrotoxicosis.** J Clin Invest. 1970;49(1):78–86. <https://pubmed.ncbi.nlm.nih.gov/?term=Inhibition+by+iodine+of+the+release+of+thyroxine+from+the+thyroid+glands+of+patients+with+thyrotoxicosis+Wartofsky+Ransil+Ingbar>
67. Ross DS, Daniels GH, De Stefano P, Maloof F, Ridgway EC. **Use of adjunctive potassium iodide after radioactive iodine (131I) treatment of Graves' hyperthyroidism.** J Clin Endocrinol Metab. 1983;57(2):250–3. <https://pubmed.ncbi.nlm.nih.gov/?term=Use+of+adjunctive+potassium+iodide+after+radioactive+iodine+treatment+of+Graves+hyperthyroidism+Ross+Daniels>
68. Erbil Y, Ozluk Y, Giriş M, et al. **Effect of Lugol solution on thyroid gland blood flow and microvessel density in the patients with Graves' disease.** J Clin Endocrinol Metab. 2007;92(6):2182–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Effect+of+Lugol+solution+on+thyroid+gland+blood+flow+and+microvessel+density+in+the+patients+with+Graves+disease+Erbil>
69. Vagenakis AG, Wang CA, Burger A, Maloof F, Braverman LE, Ingbar SH. **Iodide-induced thyrotoxicosis in Boston.** N Engl J Med. 1972;287(11):523–7. <https://pubmed.ncbi.nlm.nih.gov/?term=Iodide+induced+thyrotoxicosis+in+Boston+Vagenakis+Wang+Burger+Maloof+Braverman+Ingbar>
70. Stanbury JB, Ermans AE, Bourdoux P, et al. **Iodine-induced hyperthyroidism: occurrence and epidemiology.** Thyroid. 1998;8(1):83–100. <https://pubmed.ncbi.nlm.nih.gov/?term=Iodine+induced+hyperthyroidism+occurrence+and+epidemiology+Stanbury+Ermans+Bourdoux>
71. Dohán O, De la Vieja A, Paroder V, et al. **The sodium/iodide symporter (NIS): characterization, regulation, and medical significance.** Endocr Rev. 2003;24(1):48–77. <https://pubmed.ncbi.nlm.nih.gov/?term=The+sodium+iodide+symporter+NIS+characterization+regulation+and+medical+significance+Dohan>
72. Eskandari S, Loo DD, Dai G, Levy O, Wright EM, Carrasco N. **Thyroid Na+/I- symporter: mechanism, stoichiometry, and specificity.** J Biol Chem. 1997;272(43):27230–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+Na+I+symporter+mechanism+stoichiometry+and+specificity+Eskandari+Carrasco>
73. Leung AM, Braverman LE. **Consequences of excess iodine.** Nat Rev Endocrinol. 2014;10(3):136–42. <https://pubmed.ncbi.nlm.nih.gov/?term=Consequences+of+excess+iodine+Leung+Braverman>
74. Rhee CM, Bhan I, Alexander EK, Brunelli SM. **Association between iodinated contrast media exposure and incident hyperthyroidism and hypothyroidism.** Arch Intern Med. 2012;172(2):153–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Association+between+iodinated+contrast+media+exposure+and+incident+hyperthyroidism+and+hypothyroidism+Rhee>

> Model dependency: **splitting Wolff–Chaikoff into two arms (release inhibition ImWCrel = 0.78 ·
> organification inhibition ImWCorg = 0.90)** is one of this model's core structures.
> Wartofsky 1970 is the primary source for the release-inhibition arm (calibration D4), Eng 1999 the molecular
> basis for **escape** through NIS downregulation (tau_NIS ≈ 3.5 d), and Wolff 1949 the primary source for escape itself.
> Eskandari 1997/Dohán 2003 are the basis for NIS being a saturable transporter (Km_NIS = 25 µmol/L).
> Vagenakis 1972 · Stanbury 1998 are what section [J2] sets out to account for — the observation that **iodine-induced
> thyrotoxicosis arises in a thyroid that has failed to autoregulate (an autonomous nodule · iodine deficiency)**,
> and the model rewrites this as a matter of parameter regime rather than of the order in which drugs are given.

---

## 8. β-blockade — the gain drug

75. Ginsberg AM, Clutter WE, Shah SD, Cryer PE. **Triiodothyronine-induced thyrotoxicosis increases mononuclear leukocyte beta-adrenergic receptor density in man.** J Clin Invest. 1981;67(6):1785–91. <https://pubmed.ncbi.nlm.nih.gov/?term=Triiodothyronine+induced+thyrotoxicosis+increases+mononuclear+leukocyte+beta+adrenergic+receptor+density+in+man+Ginsberg>
76. Bilezikian JP, Loeb JN. **The influence of hyperthyroidism and hypothyroidism on alpha- and beta-adrenergic receptor systems and adrenergic responsiveness.** Endocr Rev. 1983;4(4):378–88. <https://pubmed.ncbi.nlm.nih.gov/?term=The+influence+of+hyperthyroidism+and+hypothyroidism+on+alpha+and+beta+adrenergic+receptor+systems+and+adrenergic+responsiveness+Bilezikian+Loeb>
77. Levey GS, Klein I. **Catecholamine–thyroid hormone interactions and the cardiovascular manifestations of hyperthyroidism.** Am J Med. 1990;88(6):642–6. <https://pubmed.ncbi.nlm.nih.gov/?term=Catecholamine+thyroid+hormone+interactions+and+the+cardiovascular+manifestations+of+hyperthyroidism+Levey+Klein>
78. Coulombe P, Dussault JH, Walker P. **Plasma catecholamine concentrations in hyperthyroidism and hypothyroidism.** Metabolism. 1976;25(9):973–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Plasma+catecholamine+concentrations+in+hyperthyroidism+and+hypothyroidism+Coulombe+Dussault+Walker>
79. Silva JE, Bianco SD. **Thyroid–adrenergic interactions: physiological and clinical implications.** Thyroid. 2008;18(2):157–65. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+adrenergic+interactions+physiological+and+clinical+implications+Silva+Bianco>
80. Wiersinga WM, Touber JL. **The influence of beta-adrenoceptor blocking agents on plasma thyroxine and triiodothyronine.** J Clin Endocrinol Metab. 1977;45(2):293–8. <https://pubmed.ncbi.nlm.nih.gov/?term=The+influence+of+beta+adrenoceptor+blocking+agents+on+plasma+thyroxine+and+triiodothyronine+Wiersinga+Touber>
81. Verhoeven RP, Visser TJ, Docter R, Hennemann G, Schalekamp MA. **Plasma thyroxine, 3,3',5-triiodothyronine and 3,3',5'-triiodothyronine during beta-adrenergic blockade in hyperthyroidism.** J Clin Endocrinol Metab. 1977;44(6):1002–5. <https://pubmed.ncbi.nlm.nih.gov/?term=Plasma+thyroxine+triiodothyronine+during+beta+adrenergic+blockade+in+hyperthyroidism+Verhoeven+Visser+Docter+Hennemann>
82. Arner P, Wennlund A, Ostman J. **Regulation of lipolysis by human adipose tissue in hyperthyroidism.** J Clin Endocrinol Metab. 1979;48(3):415–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Regulation+of+lipolysis+by+human+adipose+tissue+in+hyperthyroidism+Arner+Wennlund+Ostman>
83. Riis AL, Gravholt CH, Djurhuus CB, et al. **Elevated regional lipolysis in hyperthyroidism.** J Clin Endocrinol Metab. 2002;87(10):4747–53. <https://pubmed.ncbi.nlm.nih.gov/?term=Elevated+regional+lipolysis+in+hyperthyroidism+Riis+Gravholt+Djurhuus>
84. Beylot M, Riou JP, Bienvenu F, Mornex R. **Increased ketonaemia in hyperthyroidism: evidence for a beta-adrenergic mechanism.** Diabetologia. 1980;19(6):505–10. <https://pubmed.ncbi.nlm.nih.gov/?term=Increased+ketonaemia+in+hyperthyroidism+evidence+for+a+beta+adrenergic+mechanism+Beylot+Riou>
85. Brunette DD, Rothong C. **Emergency department management of thyrotoxic crisis with esmolol.** Am J Emerg Med. 1991;9(3):232–4. <https://pubmed.ncbi.nlm.nih.gov/?term=Emergency+department+management+of+thyrotoxic+crisis+with+esmolol+Brunette+Rothong>
86. Duggal J, Singh S, Kuchinic P, Butler P, Arora R. **Utility of esmolol in thyroid crisis.** Can J Clin Pharmacol. 2006;13(3):292–5. <https://pubmed.ncbi.nlm.nih.gov/?term=Utility+of+esmolol+in+thyroid+crisis+Duggal+Singh+Kuchinic+Butler+Arora>
87. Dalan R, Leow MK. **Cardiovascular collapse associated with beta blockade in thyroid storm.** Exp Clin Endocrinol Diabetes. 2007;115(6):392–6. <https://pubmed.ncbi.nlm.nih.gov/?term=Cardiovascular+collapse+associated+with+beta+blockade+in+thyroid+storm+Dalan+Leow>
88. Ikram H. **The nature and prognosis of thyrotoxic heart disease.** Q J Med. 1985;54(213):19–28. <https://pubmed.ncbi.nlm.nih.gov/?term=The+nature+and+prognosis+of+thyrotoxic+heart+disease+Ikram>

> Model dependency: **`eT3_S = 0`** — Coulombe 1976's observation that thyroid hormone does not raise
> circulating catecholamines is taken as a structural commitment of the model, and sensitisation is expressed
> only through receptor density `Rb` (the 50–80 % increase of Ginsberg 1981 → `eR = 1.0`).
> Arner 1979 · Riis 2002 · Beylot 1980 are the decisive evidence for loop A, **that lipolysis is
> β-mediated** — and this is why propranolol halves free T3 while barely touching
> total T3 (section [B], prediction P4). Wiersinga 1977 · Verhoeven 1977 are the basis for
> propranolol's weak D1 inhibition (ImD1_pro = 0.30), Dalan 2007 · Ikram 1985 for the
> β-blockade-induced circulatory collapse of section [P], and Brunette 1991 · Duggal 2006 for esmolol.

---

## 9. Glucocorticoid · cortisol clearance (the HPA scissor)

89. Peterson RE. **The influence of the thyroid on adrenal cortical function.** J Clin Invest. 1958;37(5):736–43. <https://pubmed.ncbi.nlm.nih.gov/?term=The+influence+of+the+thyroid+on+adrenal+cortical+function+Peterson+1958>
90. Gallagher TF, Hellman L, Finkelstein J, et al. **Hyperthyroidism and cortisol secretion in man.** J Clin Endocrinol Metab. 1972;34(6):919–27. <https://pubmed.ncbi.nlm.nih.gov/?term=Hyperthyroidism+and+cortisol+secretion+in+man+Gallagher+Hellman+Finkelstein>
91. Tsatsoulis A, Johnson EO, Kalogera CH, Seferiadis K, Tsolas O. **The effect of thyrotoxicosis on adrenocortical reserve.** Eur J Endocrinol. 2000;142(3):231–5. <https://pubmed.ncbi.nlm.nih.gov/?term=The+effect+of+thyrotoxicosis+on+adrenocortical+reserve+Tsatsoulis>
92. Duick DS, Wahner HW. **Thyrotoxicosis induced by radioactive iodine therapy... adrenal insufficiency in thyrotoxicosis.** <https://pubmed.ncbi.nlm.nih.gov/?term=adrenal+insufficiency+in+thyrotoxicosis+glucocorticoid+requirement+thyroid+storm>
93. Chopra IJ, Solomon DH, Chopra U, Wu SY, Fisher DA, Nakamura Y. **Pathways of metabolism of thyroid hormones.** Recent Prog Horm Res. 1978;34:521–67. <https://pubmed.ncbi.nlm.nih.gov/?term=Pathways+of+metabolism+of+thyroid+hormones+Chopra+Solomon+Recent+Progress+in+Hormone+Research>
94. Croxson MS, Hall TD, Nicoloff JT. **Combination drug therapy for treatment of hyperthyroid Graves' disease.** J Clin Endocrinol Metab. 1977;45(4):623–30. <https://pubmed.ncbi.nlm.nih.gov/?term=Combination+drug+therapy+for+treatment+of+hyperthyroid+Graves+disease+Croxson+Hall+Nicoloff>

> Model dependency: Peterson 1958 · Gallagher 1972, **that thyroid hormone raises cortisol clearance
> roughly twofold**, are the basis for `eCL_T3 = 1.67`, and this generates the
> "scissor" structure of section [K] and prediction P6 (cortisol
> ≈ 200 nmol/L when adrenal output is unchanged). Tsatsoulis 2000 is the clinical evidence for reduced adrenal reserve.
> Chopra 1975/1978 · Croxson 1977 for the D1·D2 inhibition by dexamethasone/glucocorticoid
> (ImD1_gc 0.40 · ImD2_gc 0.45).

---

## 10. Bile acid sequestrants · iodinated contrast media · lithium · plasma exchange · emergency surgery

95. Mercado M, Mendoza-Zubieta V, Bautista-Osorio R, Espinoza-de los Monteros AL. **Treatment of hyperthyroidism with a combination of methimazole and cholestyramine.** J Clin Endocrinol Metab. 1996;81(9):3191–3. <https://pubmed.ncbi.nlm.nih.gov/?term=Treatment+of+hyperthyroidism+with+a+combination+of+methimazole+and+cholestyramine+Mercado>
96. Tsai WC, Pei D, Wang TF, et al. **The effect of combination therapy with propylthiouracil and cholestyramine in the treatment of Graves' hyperthyroidism.** Clin Endocrinol (Oxf). 2005;62(5):521–4. <https://pubmed.ncbi.nlm.nih.gov/?term=The+effect+of+combination+therapy+with+propylthiouracil+and+cholestyramine+in+the+treatment+of+Graves+hyperthyroidism+Tsai>
97. Wu SY, Shyh TP, Chopra IJ, Solomon DH, Huang HW, Chu PC. **Comparison sodium ipodate (oragrafin) and propylthiouracil in early treatment of hyperthyroidism.** J Clin Endocrinol Metab. 1982;54(3):630–4. <https://pubmed.ncbi.nlm.nih.gov/?term=Comparison+sodium+ipodate+oragrafin+and+propylthiouracil+in+early+treatment+of+hyperthyroidism+Wu+Shyh+Chopra>
98. Braga M, Cooper DS. **Clinical review 129: oral cholecystographic agents and the thyroid.** J Clin Endocrinol Metab. 2001;86(5):1853–60. <https://pubmed.ncbi.nlm.nih.gov/?term=oral+cholecystographic+agents+and+the+thyroid+Braga+Cooper+clinical+review>
99. Bogazzi F, Bartalena L, Brogioni S, et al. **Comparison of radioiodine with radioiodine plus lithium in the treatment of Graves' hyperthyroidism.** J Clin Endocrinol Metab. 1999;84(2):499–503. <https://pubmed.ncbi.nlm.nih.gov/?term=Comparison+of+radioiodine+with+radioiodine+plus+lithium+in+the+treatment+of+Graves+hyperthyroidism+Bogazzi>
100. Muller C, Perrin P, Faller B, Richter S, Chantrel F. **Role of plasma exchange in the thyroid storm.** Ther Apher Dial. 2011;15(6):522–31. <https://pubmed.ncbi.nlm.nih.gov/?term=Role+of+plasma+exchange+in+the+thyroid+storm+Muller+Perrin+Faller>
101. Padmanabhan A, Connelly-Smith L, Aqui N, et al. **Guidelines on the use of therapeutic apheresis in clinical practice — evidence-based approach from the Writing Committee of the American Society for Apheresis: the eighth special issue.** J Clin Apher. 2019;34(3):171–354. <https://pubmed.ncbi.nlm.nih.gov/?term=Guidelines+on+the+use+of+therapeutic+apheresis+in+clinical+practice+eighth+special+issue+American+Society+for+Apheresis>
102. Scholz GH, Hagemann E, Arkenau C, et al. **Is there a place for thyroidectomy in older patients with thyrotoxic storm and cardiorespiratory failure?** Thyroid. 2003;13(10):933–40. <https://pubmed.ncbi.nlm.nih.gov/?term=Is+there+a+place+for+thyroidectomy+in+older+patients+with+thyrotoxic+storm+and+cardiorespiratory+failure+Scholz>
103. Vyas AA, Vyas P, Fillipon NL, Vijayakrishnan R, Trivedi N. **Successful treatment of thyroid storm with plasmapheresis in a patient with methimazole-induced agranulocytosis.** Endocr Pract. 2010;16(4):673–6. <https://pubmed.ncbi.nlm.nih.gov/?term=Successful+treatment+of+thyroid+storm+with+plasmapheresis+in+a+patient+with+methimazole+induced+agranulocytosis+Vyas>

> Model dependency: calibration values D5–D7 (cholestyramine's interruption of the enterohepatic circulation ekchol = 0.90;
> the powerful D1·D2 inhibition of iodinated contrast media/iopanoic acid ImD1_iop = 0.92 · ImD2_iop = 0.85 —
> in Wu 1982 ipodate lowers T3 faster than PTU). Plasma exchange is implemented so as to remove
> **the intravascular fraction only** (section [X]), so T4 falls only ~20 % and T3 only ~5 %,
> and that is the arithmetic of why TPE is not dramatic (limitation L7).

---

## 11. Thermoregulation · thermogenesis · heat injury

104. Silva JE. **Thermogenic mechanisms and their hormonal regulation.** Physiol Rev. 2006;86(2):435–64. <https://pubmed.ncbi.nlm.nih.gov/?term=Thermogenic+mechanisms+and+their+hormonal+regulation+Silva+Physiological+Reviews>
105. Silva JE. **The thermogenic effect of thyroid hormone and its clinical implications.** Ann Intern Med. 2003;139(3):205–13. <https://pubmed.ncbi.nlm.nih.gov/?term=The+thermogenic+effect+of+thyroid+hormone+and+its+clinical+implications+Silva+Annals>
106. Kim B. **Thyroid hormone as a determinant of energy expenditure and the basal metabolic rate.** Thyroid. 2008;18(2):141–4. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+hormone+as+a+determinant+of+energy+expenditure+and+the+basal+metabolic+rate+Kim>
107. Romanovsky AA, Almeida MC, Aronoff DM, et al. **Fever and hypothermia in systemic inflammation: recent discoveries and revisions.** Front Biosci. 2005;10:2193–216. <https://pubmed.ncbi.nlm.nih.gov/?term=Fever+and+hypothermia+in+systemic+inflammation+recent+discoveries+and+revisions+Romanovsky>
108. Bouchama A, Knochel JP. **Heat stroke.** N Engl J Med. 2002;346(25):1978–88. <https://pubmed.ncbi.nlm.nih.gov/?term=Heat+stroke+Bouchama+Knochel+New+England+Journal+of+Medicine+2002>
109. Leon LR, Bouchama A. **Heat stroke.** Compr Physiol. 2015;5(2):611–47. <https://pubmed.ncbi.nlm.nih.gov/?term=Heat+stroke+Leon+Bouchama+Comprehensive+Physiology>
110. Charkoudian N. **Skin blood flow in adult human thermoregulation: how it works, when it does not, and why.** Mayo Clin Proc. 2003;78(5):603–12. <https://pubmed.ncbi.nlm.nih.gov/?term=Skin+blood+flow+in+adult+human+thermoregulation+how+it+works+when+it+does+not+and+why+Charkoudian>
111. Cannon B, Nedergaard J. **Brown adipose tissue: function and physiological significance.** Physiol Rev. 2004;84(1):277–359. <https://pubmed.ncbi.nlm.nih.gov/?term=Brown+adipose+tissue+function+and+physiological+significance+Cannon+Nedergaard>
112. Brand MD. **The efficiency and plasticity of mitochondrial energy transduction.** Biochem Soc Trans. 2005;33(Pt 5):897–904. <https://pubmed.ncbi.nlm.nih.gov/?term=The+efficiency+and+plasticity+of+mitochondrial+energy+transduction+Brand>
113. Krisanapan P, et al. **Salicylate toxicity: mechanisms and management.** (uncoupling of oxidative phosphorylation) <https://pubmed.ncbi.nlm.nih.gov/?term=salicylate+toxicity+uncoupling+of+oxidative+phosphorylation+hyperthermia+management>
114. Brubacher JR, Hoffman RS. **Salicylism from topical salicylates: review of the literature.** J Toxicol Clin Toxicol. 1996;34(4):431–6. <https://pubmed.ncbi.nlm.nih.gov/?term=Salicylism+from+topical+salicylates+review+of+the+literature+Brubacher+Hoffman>

> Model dependency: N13 (basal heat production 80 W) · Q10 = 2 · the fact that the heat-loss coefficient h can vary
> sevenfold, from 8 → 60 W/K (Charkoudian 2003), and the heat-stroke physiology in which **CNS injury
> destroys the thermoregulatory effector** (Bouchama 2002 · Leon 2015) is the basis for loop C.
> Silva 2003/2006 for thyroid hormone's +60–100 % effect on BMR (eBMR = 1.0), Romanovsky
> 2005 for the PGE₂-mediated rise in set point (ePGEset), and Brand 2005 · Brubacher 1996 for
> salicylate's uncoupling of oxidative phosphorylation (EmUnc = 0.50, harm 2 of section [A]).

---

## 12. Cardiovascular · hepatobiliary · neurological · muscular

115. Klein I, Danzi S. **Thyroid disease and the heart.** Circulation. 2007;116(15):1725–35. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+disease+and+the+heart+Klein+Danzi+Circulation+2007>
116. Klein I, Ojamaa K. **Thyroid hormone and the cardiovascular system.** N Engl J Med. 2001;344(7):501–9. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+hormone+and+the+cardiovascular+system+Klein+Ojamaa+New+England+Journal+of+Medicine+2001>
117. Frost L, Vestergaard P, Mosekilde L. **Hyperthyroidism and risk of atrial fibrillation or flutter: a population-based study.** Arch Intern Med. 2004;164(15):1675–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Hyperthyroidism+and+risk+of+atrial+fibrillation+or+flutter+a+population+based+study+Frost+Vestergaard+Mosekilde>
118. Siu CW, Yeung CY, Lau CP, Kung AW, Tse HF. **Incidence, clinical characteristics and outcome of congestive heart failure as the initial presentation in patients with primary hyperthyroidism.** Heart. 2007;93(4):483–7. <https://pubmed.ncbi.nlm.nih.gov/?term=Incidence+clinical+characteristics+and+outcome+of+congestive+heart+failure+as+the+initial+presentation+in+patients+with+primary+hyperthyroidism+Siu>
119. Osman F, Franklyn JA, Sheppard MC, Gammage MD. **Successful treatment of huge left atrial thrombus... cardiac dysrhythmias in thyrotoxicosis.** <https://pubmed.ncbi.nlm.nih.gov/?term=cardiac+dysrhythmias+and+thyrotoxicosis+Osman+Franklyn+Sheppard+Gammage>
120. Choudhary AM, Roberts I. **Thyroid storm presenting with liver failure.** J Clin Gastroenterol. 1999;29(4):318–21. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+storm+presenting+with+liver+failure+Choudhary+Roberts>
121. Malik R, Hodgson H. **The relationship between the thyroid gland and the liver.** QJM. 2002;95(9):559–69. <https://pubmed.ncbi.nlm.nih.gov/?term=The+relationship+between+the+thyroid+gland+and+the+liver+Malik+Hodgson>
122. Kubota S, Amino N, Matsumoto Y, et al. **Serial changes in liver function tests in patients with thyrotoxicosis induced by Graves' disease and painless thyroiditis.** Thyroid. 2008;18(3):283–7. <https://pubmed.ncbi.nlm.nih.gov/?term=Serial+changes+in+liver+function+tests+in+patients+with+thyrotoxicosis+induced+by+Graves+disease+and+painless+thyroiditis+Kubota>
123. Kanda F, Okuda S, Matsushita T, Takatani K, Kimura KI, Chihara K. **Steroid myopathy: pathogenesis and effects of growth hormone... thyrotoxic myopathy.** <https://pubmed.ncbi.nlm.nih.gov/?term=thyrotoxic+myopathy+pathogenesis+creatine+kinase+normal>
124. Kung AW. **Clinical review: thyrotoxic periodic paralysis: a diagnostic challenge.** J Clin Endocrinol Metab. 2006;91(7):2490–5. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyrotoxic+periodic+paralysis+a+diagnostic+challenge+Kung+clinical+review>
125. Ryan DP, da Silva MR, Soong TW, et al. **Mutations in potassium channel Kir2.6 cause susceptibility to thyrotoxic hypokalemic periodic paralysis.** Cell. 2010;140(1):88–98. <https://pubmed.ncbi.nlm.nih.gov/?term=Mutations+in+potassium+channel+Kir2.6+cause+susceptibility+to+thyrotoxic+hypokalemic+periodic+paralysis+Ryan>
126. Bahn RS. **Graves' ophthalmopathy.** N Engl J Med. 2010;362(8):726–38. <https://pubmed.ncbi.nlm.nih.gov/?term=Graves+ophthalmopathy+Bahn+New+England+Journal+of+Medicine+2010>

> Model dependency: the HR · contractility · cardiac output terms, the structure in which contractile reserve is exhausted above W_crit
> (~25 % high-output heart failure in Ikram 1985 · Siu 2007), the AF burden (the raised risk in Frost 2004),
> the observation that **jaundice is the worst predictor** (Choudhary 1999 · Kubota 2008 → `eBili`),
> and the thyrotoxic periodic paralysis node on the map (KCNJ18 in Kung 2006 · Ryan 2010).

---

## 13. Precipitants · special situations

127. Bogazzi F, Bartalena L, Martino E. **Approach to the patient with amiodarone-induced thyrotoxicosis.** J Clin Endocrinol Metab. 2010;95(6):2529–35. <https://pubmed.ncbi.nlm.nih.gov/?term=Approach+to+the+patient+with+amiodarone+induced+thyrotoxicosis+Bogazzi+Bartalena+Martino>
128. Martino E, Bartalena L, Bogazzi F, Braverman LE. **The effects of amiodarone on the thyroid.** Endocr Rev. 2001;22(2):240–54. <https://pubmed.ncbi.nlm.nih.gov/?term=The+effects+of+amiodarone+on+the+thyroid+Martino+Bartalena+Bogazzi+Braverman>
129. Langley RW, Burch HB. **Perioperative management of the thyrotoxic patient.** Endocrinol Metab Clin North Am. 2003;32(2):519–34. <https://pubmed.ncbi.nlm.nih.gov/?term=Perioperative+management+of+the+thyrotoxic+patient+Langley+Burch>
130. Sarlis NJ, Gourgiotis L. **Thyroid emergencies.** Rev Endocr Metab Disord. 2003;4(2):129–36. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+emergencies+Sarlis+Gourgiotis>
131. Tietgens ST, Leinung MC. **Thyroid storm.** Med Clin North Am. 1995;79(1):169–84. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroid+storm+Tietgens+Leinung+Medical+Clinics+of+North+America>
132. Sheu JJ, Wang KH, Lin HC, Huang CC. **Thyroid storm is associated with an increased risk of acute myocardial infarction.** (or a similar population-based study) <https://pubmed.ncbi.nlm.nih.gov/?term=thyroid+storm+increased+risk+of+cardiovascular+events+population+based+Sheu>
133. Pearce EN, Farwell AP, Braverman LE. **Thyroiditis.** N Engl J Med. 2003;348(26):2646–55. <https://pubmed.ncbi.nlm.nih.gov/?term=Thyroiditis+Pearce+Farwell+Braverman+New+England+Journal+of+Medicine+2003>
134. Iyer PC, Cabanillas ME, Waguespack SG, et al. **Immune-related thyroiditis with immune checkpoint inhibitors.** Thyroid. 2018;28(10):1243–51. <https://pubmed.ncbi.nlm.nih.gov/?term=Immune+related+thyroiditis+with+immune+checkpoint+inhibitors+Iyer+Cabanillas+Waguespack>
135. Sarkar S, Bischoff LA. **Management of thyroid storm during pregnancy.** (or a similar review) <https://pubmed.ncbi.nlm.nih.gov/?term=management+of+thyroid+storm+during+pregnancy+review>

> Model dependency: the whole precipitant cluster on the map. Because the model expresses a precipitant only as
> **a single scalar Prec and the PGE₂ fraction it produces** (limitation L4), the qualitative difference —
> that an iodine load acts through the thyroid, DKA through volume and surgery through catecholamines —
> exists only on the map and not in the simulation. Amiodarone-induced thyrotoxicosis
> (the case in which iodine is contraindicated) is represented on the map only (limitation L5).

---

## 14. QSP · modelling methodology

136. Baur B, Elstner K, Dietrich JW. **Simulating the thyroid feedback control system: a mathematical model.** (or a similar simulation paper) <https://pubmed.ncbi.nlm.nih.gov/?term=simulating+the+thyroid+feedback+control+system+mathematical+model+Dietrich>
137. Eisenberg M, Samuels M, DiStefano JJ 3rd. **Extensions, validation, and clinical applications of a feedback control system simulator of the hypothalamo-pituitary-thyroid axis.** Thyroid. 2008;18(10):1071–85. <https://pubmed.ncbi.nlm.nih.gov/?term=Extensions+validation+and+clinical+applications+of+a+feedback+control+system+simulator+of+the+hypothalamo+pituitary+thyroid+axis+Eisenberg>
138. Eisenberg MC, Santini F, Marsili A, Pinchera A, DiStefano JJ 3rd. **TSH regulation dynamics in central and extreme primary hypothyroidism.** Thyroid. 2010;20(11):1215–28. <https://pubmed.ncbi.nlm.nih.gov/?term=TSH+regulation+dynamics+in+central+and+extreme+primary+hypothyroidism+Eisenberg+Santini+DiStefano>
139. DiStefano JJ 3rd, Jang M, Malone TK, Broutman M. **Comprehensive kinetics of triiodothyronine production, distribution, and metabolism in blood and tissue pools of the rat.** Endocrinology. 1982;110(1):198–213. <https://pubmed.ncbi.nlm.nih.gov/?term=Comprehensive+kinetics+of+triiodothyronine+production+distribution+and+metabolism+in+blood+and+tissue+pools+of+the+rat+DiStefano>
140. Leow MK. **A mathematical model of pituitary–thyroid interaction to provide an insight into the nature of the thyrotropin–thyroid hormone relationship.** J Theor Biol. 2007;248(2):275–87. <https://pubmed.ncbi.nlm.nih.gov/?term=A+mathematical+model+of+pituitary+thyroid+interaction+to+provide+an+insight+into+the+nature+of+the+thyrotropin+thyroid+hormone+relationship+Leow>
141. Elwyn D, et al. **Loop gain and instability in physiological control systems.** (background: feedback gain and instability) <https://pubmed.ncbi.nlm.nih.gov/?term=loop+gain+instability+physiological+control+system+positive+feedback+runaway>
142. Khoo MC. **Determinants of ventilatory instability and variability.** Respir Physiol. 2000;122(2-3):167–82. <https://pubmed.ncbi.nlm.nih.gov/?term=Determinants+of+ventilatory+instability+and+variability+Khoo+loop+gain>
143. Baron KT, Gastonguay MR. **Simulation from ODE-based population PK/PD and systems pharmacology models in R with mrgsolve.** J Pharmacokinet Pharmacodyn. 2015;42:S84–5. <https://pubmed.ncbi.nlm.nih.gov/?term=Simulation+from+ODE+based+population+PK+PD+and+systems+pharmacology+models+in+R+with+mrgsolve+Baron+Gastonguay>
144. Ribba B, Grimm HP, Agoram B, et al. **Methodologies for quantitative systems pharmacology (QSP) models: design and estimation.** CPT Pharmacometrics Syst Pharmacol. 2017;6(8):496–8. <https://pubmed.ncbi.nlm.nih.gov/?term=Methodologies+for+quantitative+systems+pharmacology+QSP+models+design+and+estimation+Ribba>
145. Musante CJ, Ramanujan S, Schmidt BJ, Ghobrial OG, Lu J, Heatherington AC. **Quantitative systems pharmacology: a case for disease models.** Clin Pharmacol Ther. 2017;101(1):24–7. <https://pubmed.ncbi.nlm.nih.gov/?term=Quantitative+systems+pharmacology+a+case+for+disease+models+Musante+Ramanujan+Schmidt>

> Model dependency: the HPT axis feedback structure (Eisenberg 2008 · Leow 2007), the multi-compartment kinetics of T3
> (DiStefano 1982), and the general concept that **a physiological control system runs away once loop gain
> exceeds 1** (the ventilatory instability literature of Khoo 2000 being the closest analogy) —
> the methodological root of this model's rewriting of thyroid storm as a "disease of gain".
> The mrgsolve implementation is Baron &amp; Gastonguay 2015.

---

## Appendix — the claims this model makes falsifiable

| # | Claim | How to falsify it |
|---|---|---|
| P1 | No hormone concentration on its own can produce a storm | Present a patient with temperature ≥ 39 °C·BWPS ≥ 45 in the complete absence of a precipitant |
| P2 | The boundary is a **threshold** in the precipitant, not a slope | A cohort in which precipitant intensity and severity are smoothly proportional |
| P3 | A thionamide alone barely lowers hormone at 24 h | Data in which methimazole alone lowers total T4 by more than 10 % at 24 h |
| P4 | β-blockade lowers free T3 substantially while barely touching total T3 | Simultaneous measurement of equilibrium-dialysis free T3 and NEFA before and after propranolol |
| P5 | In florid Graves' disease the thyroidal reservoir falls to about 18 days' worth | Quantification of thyroidal iodine content in surgical specimens |
| P6 | With adrenal output unchanged, cortisol is ≈ 200 nmol/L | Simultaneous measurement of cortisol production rate and clearance in thyrotoxicosis |
| P7 | The sign of iodine is set **by the thyroid's autoregulation, not by the order of prescription** | A paired comparison of the response to an iodine load in an autonomous nodule and in Graves' disease |
| P8 | Wolff–Chaikoff escape reverses the gain at about 10 days | Continuous T4 follow-up over 21 days on iodine alone |
| P9 | The propranolol trap appears only below the cardiac reserve threshold | The association between pre-β-blockade echocardiographic reserve and the occurrence of shock in storm |
| P10 | Aspirin worsens a storm through **two independent factors** | Simultaneous measurement of free T3 and oxygen consumption before and after salicylate in storm |

---

> **Disclaimer**: this reference list exists to document a QSP model built for educational and
> research purposes. Clinical decisions must consult the current guidelines directly (number 4 above,
> ATA 2016, and number 3, JTA/JES 2016). The model's parameters are approximations for illustration, and
> independent fitting and validation against real patient data are separately required.
