# Male Hypogonadism QSP Model — References

**Male Hypogonadism · MHG**

This document sets out the sources of the parameters, structures, and calibration targets used in
`mhg_qsp_model.dot` (the mechanistic map), `mhg_mrgsolve_model.R` (the ODE model), and `mhg_shiny_app.R` (the dashboard).
For each item classified as "ANCHORED / FITTED / CALIBRATED" in the **CALIBRATION BLOCK** of the model file,
the reference below that it corresponds to is marked.

---

## 1. Diagnostic criteria and guidelines

1. Bhasin S, et al. **Testosterone Therapy in Men With Hypogonadism: An Endocrine
   Society Clinical Practice Guideline.** J Clin Endocrinol Metab. 2018;103(5):1715-1744.
   <https://pubmed.ncbi.nlm.nih.gov/29562364/>
   — the diagnostic algorithm, two morning fasting total T measurements, the 300 ng/dL threshold, the monitoring schedule.
   The rules that `MHG_diagnostic_frame()` takes as its object of validation.

2. Travison TG, et al. **Harmonized Reference Ranges for Circulating Testosterone
   Levels in Men of Four Cohort Studies in the United States and Europe.**
   J Clin Endocrinol Metab. 2017;102(4):1161-1173.
   <https://pubmed.ncbi.nlm.nih.gov/28324103/>
   — CDC-harmonised reference values (2.5th percentile 264 ng/dL for ages 19-39). The basis for the 300 ng/dL threshold.

3. Mulhall JP, et al. **Evaluation and Management of Testosterone Deficiency:
   AUA Guideline.** J Urol. 2018;200(2):423-432.
   <https://pubmed.ncbi.nlm.nih.gov/29601923/>

4. Salonia A, et al. **European Association of Urology Guidelines on Sexual and
   Reproductive Health — 2021 Update: Male Sexual Dysfunction.** Eur Urol.
   2021;80(3):333-357. <https://pubmed.ncbi.nlm.nih.gov/34183196/>

5. Rosner W, et al. **Position statement: Utility, limitations, and pitfalls in
   measuring testosterone: an Endocrine Society position statement.**
   J Clin Endocrinol Metab. 2007;92(2):405-413.
   <https://pubmed.ncbi.nlm.nih.gov/17090633/>
   — why immunoassays are inaccurate at low concentrations, and the LC-MS/MS reference method.

---

## 2. Binding equilibrium — the mathematical centre of this model

6. Vermeulen A, Verdonck L, Kaufman JM. **A critical evaluation of simple methods
   for the estimation of free testosterone in serum.** J Clin Endocrinol Metab.
   1999;84(10):3666-3672. <https://pubmed.ncbi.nlm.nih.gov/10523012/>
   — **ANCHORED**: `KS = 1.0e9 M^-1`, `KA = 3.6e4 M^-1`, and the quadratic equation that
   `free_T_pgmL()` in the model's `$GLOBAL` solves, itself.

7. Zakharov MN, et al. **A multi-step, dynamic allosteric model of testosterone's
   binding to sex hormone binding globulin.** Mol Cell Endocrinol. 2015;399:190-200.
   <https://pubmed.ncbi.nlm.nih.gov/25240469/>
   — an allosteric binding model for the SHBG dimer. It points out the limits of the simple Vermeulen model and
   shows that it can underestimate free T at low SHBG. This model adopted Vermeulen for computability and
   clinical convention, and that is a **deliberate simplification**.

8. Handelsman DJ. **Free Testosterone: Pumping up the Tires or Ending the Free
   Ride?** Endocr Rev. 2017;38(4):297-301.
   <https://pubmed.ncbi.nlm.nih.gov/28486701/>
   — a methodological critique of the calculation of free T. Why the model treats free T not as "the truth" but as
   "a frame different from total T".

9. Goldman AL, et al. **A Reappraisal of Testosterone's Binding in Circulation:
   Physiological and Clinical Implications.** Endocr Rev. 2017;38(4):302-324.
   <https://pubmed.ncbi.nlm.nih.gov/28673039/>

10. Dunn JF, Nisula BC, Rodbard D. **Transport of steroid hormones: binding of 21
    endogenous steroids to both testosterone-binding globulin and
    corticosteroid-binding globulin in human plasma.** J Clin Endocrinol Metab.
    1981;53(1):58-68. <https://pubmed.ncbi.nlm.nih.gov/7195404/>
    — the origin of the albumin binding constant.

11. Ly LP, Handelsman DJ. **Empirical estimation of free testosterone from
    testosterone and sex hormone-binding globulin immunoassays.** Eur J Endocrinol.
    2005;152(3):471-478. <https://pubmed.ncbi.nlm.nih.gov/15757866/>

---

## 3. SHBG regulation

12. Selby C. **Sex hormone binding globulin: origin, function and clinical
    significance.** Ann Clin Biochem. 1990;27(Pt 6):532-541.
    <https://pubmed.ncbi.nlm.nih.gov/2080856/>

13. Simó R, et al. **Novel insights in SHBG regulation and clinical implications.**
    Trends Endocrinol Metab. 2015;26(7):376-383.
    <https://pubmed.ncbi.nlm.nih.gov/26044465/>
    — **FITTED**: `EXP_INS` (insulin), `EXP_E2`, `THYR`. Regulation mediated by HNF-4α.

14. Selva DM, et al. **Monosaccharide-induced lipogenesis regulates the human
    hepatic sex hormone-binding globulin gene.** J Clin Invest. 2007;117(12):3979-3987.
    <https://pubmed.ncbi.nlm.nih.gov/17992261/>
    — monosaccharides/de novo lipogenesis → HNF-4α ↓ → SHBG ↓. The `DNL_SHBG` node of the map.

15. Wu FCW, et al. **Identification of late-onset hypogonadism in middle-aged and
    elderly men (EMAS).** N Engl J Med. 2010;363(2):123-135.
    <https://pubmed.ncbi.nlm.nih.gov/20554979/>
    — **FITTED**: the rise in SHBG with age, and the thresholds of the symptom-hormone relationship.

16. Feldman HA, et al. **Age trends in the level of serum testosterone and other
    hormones in middle-aged men: longitudinal results from the Massachusetts Male
    Aging Study.** J Clin Endocrinol Metab. 2002;87(2):589-598.
    <https://pubmed.ncbi.nlm.nih.gov/11836290/>
    — **ANCHORED**: total T falls about 1%/year, free T about 2-3%/year, SHBG rises about 1%/year
    (`AGE_S = 0.010`).

---

## 4. Testosterone production and clearance kinetics

17. Southren AL, et al. **Plasma production rates of testosterone in normal adult
    men and women.** J Clin Endocrinol Metab. 1965;25(11):1441-1450.
    <https://pubmed.ncbi.nlm.nih.gov/5842460/>
    — **ANCHORED**: a production rate of about 6 mg/day. The model's `KSPILL x ITT0 = 6000 ug/day` identity.

18. Vierhapper H, Nowotny P, Waldhäusl W. **Determination of testosterone
    production rates in men and women using stable isotope/dilution and mass
    spectrometry.** J Clin Endocrinol Metab. 1997;82(5):1492-1496.
    <https://pubmed.ncbi.nlm.nih.gov/9141539/>

19. Coviello AD, et al. **Effects of graded doses of testosterone on erythropoiesis
    in healthy young and older men.** J Clin Endocrinol Metab. 2008;93(3):914-919.
    <https://pubmed.ncbi.nlm.nih.gov/18160461/>
    — **FITTED**: the dose-dependent rise in Hb/Hct, larger in the elderly. The main calibration target of the erythrocyte module.

---

## 5. Intratesticular testosterone and spermatogenesis

20. Jarow JP, Chen H, Rosner W, Trentacoste S, Zirkin BR. **Assessment of the
    androgen environment within the human testis: minimally invasive method to
    obtain intratesticular fluid.** J Androl. 2001;22(4):640-645.
    <https://pubmed.ncbi.nlm.nih.gov/11451361/>
    — **ANCHORED**: the key observation that ITT is some tens of times the serum level.

21. Coviello AD, et al. **Low-dose human chorionic gonadotropin maintains
    intratesticular testosterone in normal men with testosterone-induced
    gonadotropin suppression.** J Clin Endocrinol Metab. 2005;90(5):2595-2602.
    <https://pubmed.ncbi.nlm.nih.gov/15713727/>
    — **FITTED**: ITT −94% on weekly T 200 mg; the dose response of hCG 125-500 IU EOD.
    The direct calibration targets of `MHG_ITT_collapse()`, `HCG_POT`, and `ITT50_S`.

22. Roth MY, et al. **Serum LH correlates highly with intratesticular steroid
    levels in normal men.** J Androl. 2010;31(2):138-145.
    <https://pubmed.ncbi.nlm.nih.gov/19726700/>
    — **ANCHORED**: quantification of the absolute ITT value and its correlation with serum LH. The basis for the model's `ITT0 = 700 nmol/L`
    and for the structure in which `ITT_IN` is a saturating function of LH_EQ.

23. Liu PY, et al. **Rate, extent, and modifiers of spermatogenic recovery after
    hormonal male contraception: an integrated analysis.** Lancet.
    2006;367(9520):1412-1420. <https://pubmed.ncbi.nlm.nih.gov/16650651/>
    — **FITTED**: recovery kinetics after discontinuation. Median 3.4 months to reach 20 M/mL, 67% at 6 months,
    90% at 12 months, 100% at 24 months. The comparison target of `MHG_recovery_curve()`.

24. Heller CG, Clermont Y. **Kinetics of the germinal epithelium in man.**
    Recent Prog Horm Res. 1964;20:545-575.
    <https://pubmed.ncbi.nlm.nih.gov/14285045/>
    — **ANCHORED**: a spermatogenic cycle of about 74 days (`TAU_SPG`).

25. World Health Organization. **WHO laboratory manual for the examination and
    processing of human semen, 6th edition.** 2021.
    <https://www.who.int/publications/i/item/9789240030787>
    — the lower reference limit for sperm concentration of 16 M/mL (revised from 15 M/mL in the 5th edition).

26. Ramasamy R, et al. **Testosterone supplementation versus clomiphene citrate
    for hypogonadism: an age matched comparison of satisfaction and efficacy.**
    J Urol. 2014;192(3):875-879. <https://pubmed.ncbi.nlm.nih.gov/24657837/>

27. Wenker EP, et al. **The Use of HCG-Based Combination Therapy for Recovery of
    Spermatogenesis after Testosterone Use.** J Sex Med. 2015;12(6):1334-1337.
    <https://pubmed.ncbi.nlm.nih.gov/25904023/>

---

## 6. Pharmacokinetics by formulation

28. Behre HM, Nieschlag E. **Testosterone buciclate (20 Aet-1) in hypogonadal men:
    pharmacokinetics and pharmacodynamics.** J Clin Endocrinol Metab.
    1992;75(5):1204-1210. <https://pubmed.ncbi.nlm.nih.gov/1430080/>

29. Snyder PJ, Lawrence DA. **Treatment of male hypogonadism with testosterone
    enanthate.** J Clin Endocrinol Metab. 1980;51(6):1335-1339.
    <https://pubmed.ncbi.nlm.nih.gov/7440698/>
    — **FITTED**: the peak-to-trough waveform of IM enanthate (`KA_IM`, `F_IM`).

30. Swerdloff RS, et al. **Long-term pharmacokinetics of transdermal testosterone
    gel in hypogonadal men.** J Clin Endocrinol Metab. 2000;85(12):4500-4510.
    <https://pubmed.ncbi.nlm.nih.gov/11134099/>
    — **FITTED**: the quasi-steady-state profile of the transdermal gel (`KA_GEL`, `F_GEL`).

31. Wang C, et al. **Long-term testosterone gel (AndroGel) treatment maintains
    beneficial effects on sexual function and mood, lean and fat mass, and bone
    mineral density in hypogonadal men.** J Clin Endocrinol Metab.
    2004;89(5):2085-2098. <https://pubmed.ncbi.nlm.nih.gov/15126525/>

32. Morgentaler A, et al. **Long acting testosterone undecanoate therapy in men
    with hypogonadism: results of a pharmacokinetic clinical study.** J Urol.
    2008;180(6):2307-2313. <https://pubmed.ncbi.nlm.nih.gov/18930277/>
    — **FITTED**: the very slow release of IM undecanoate (`KA_TU`).

33. Swerdloff RS, et al. **A New Oral Testosterone Undecanoate Formulation
    Restores Testosterone to Normal Concentrations in Hypogonadal Men.**
    J Clin Endocrinol Metab. 2020;105(8):2515-2531.
    <https://pubmed.ncbi.nlm.nih.gov/32427336/>
    — **FITTED**: the lymphatic absorption of oral TU, its food dependence, and its high peak (`F_ORAL`, `KA_ORAL`).

34. Kaminetsky JC, et al. **A 52-Week Study of Dose Adjusted Subcutaneous
    Testosterone Enanthate in Oil Self-Administered via Disposable Auto-Injector.**
    J Urol. 2019;201(3):587-594. <https://pubmed.ncbi.nlm.nih.gov/30273615/>
    — **FITTED**: the gentle waveform of the once-weekly SC autoinjector (Xyosted) (`KA_SC`, `F_SC`).

35. Kaminetsky J, et al. **Efficacy and Safety of Natesto, a Testosterone Nasal
    Gel, in Hypogonadal Men.** J Sex Med. 2017;14(9):1109-1116.
    <https://pubmed.ncbi.nlm.nih.gov/28781215/>
    — the observation that the nasal gel causes relatively less suppression of LH/FSH — a consequence of its short peak.

36. Kovac JR, et al. **Testosterone pellet implants: a review.** Asian J Androl.
    2015;17(6):938-943. <https://pubmed.ncbi.nlm.nih.gov/25994641/>

---

## 7. Erythrocytosis — the second convexity of the model

37. Bachman E, et al. **Testosterone induces erythrocytosis via increased
    erythropoietin and suppressed hepcidin: evidence for a new erythropoietin/
    hemoglobin set point.** J Gerontol A Biol Sci Med Sci. 2014;69(6):725-735.
    <https://pubmed.ncbi.nlm.nih.gov/24158761/>
    — **FITTED**: the hepcidin-EPO set-point structure of the model itself. `IMAX_H`, `SEPO`.

38. Guo W, et al. **Testosterone Administration Inhibits Hepcidin Transcription
    and Is Associated With Increased Iron Incorporation Into Red Blood Cells.**
    Aging Cell. 2013;12(2):280-291. <https://pubmed.ncbi.nlm.nih.gov/23399021/>
    — hepcidin suppression → ferroportin retained → iron availability increased. The `IRON_F` term.

39. Ohlander SJ, Varghese B, Pastuszak AW. **Erythrocytosis Following Testosterone
    Therapy.** Sex Med Rev. 2018;6(1):77-85.
    <https://pubmed.ncbi.nlm.nih.gov/28526632/>
    — **a major comparison target**: the observation that the intramuscular formulation causes erythrocytosis about three times
    more often than the transdermal one. Precisely the gap `MHG_convexity_decomposition()` sets out to explain.

40. Jones SD Jr, Dukovac T, Sangkum P, Yafi FA, Hellstrom WJG. **Erythrocytosis
    and Polycythemia Secondary to Testosterone Replacement Therapy in the Aging
    Male.** Sex Med Rev. 2015;3(2):101-112.
    <https://pubmed.ncbi.nlm.nih.gov/27784544/>

41. Roy CN, et al. **Association of Testosterone Levels With Anemia in Older Men:
    A Controlled Clinical Trial (T-Trials Anemia Trial).** JAMA Intern Med.
    2017;177(4):480-490. <https://pubmed.ncbi.nlm.nih.gov/28241237/>
    — **FITTED**: 54% correction of unexplained anaemia (placebo 15%). The basis for the structure in which hypogonadal men
    start from a lower Hct at baseline in the model.

---

## 8. Bone — oestrogen dependence

42. Finkelstein JS, et al. **Gonadal steroids and body composition, strength, and
    sexual function in men.** N Engl J Med. 2013;369(11):1011-1022.
    <https://pubmed.ncbi.nlm.nih.gov/24024838/>
    — **a key reference · FITTED**: goserelin background + graded T gel ± anastrozole.
    Lean mass and strength are T-dependent, **body fat is E2-dependent**, and sexual function depends on both.
    The design `MHG_finkelstein()` sets out to reproduce. The direct calibration targets of `SFAT_E` and `SLEAN_T`.

43. Finkelstein JS, et al. **Gonadal Steroid-Dependent Effects on Bone Turnover
    and Bone Mineral Density in Men.** J Clin Invest. 2016;126(3):1114-1125.
    <https://pubmed.ncbi.nlm.nih.gov/26901812/>
    — the E2-dominant effect on bone turnover and BMD. `E_SCL`, `E_OC`.

44. Snyder PJ, et al. **Effect of Testosterone Treatment on Volumetric Bone
    Density and Strength in Older Men With Low Testosterone (T-Trials Bone Trial).**
    JAMA Intern Med. 2017;177(4):471-479.
    <https://pubmed.ncbi.nlm.nih.gov/28241231/>
    — **FITTED**: lumbar trabecular vBMD about +7.5% in one year. The calibration target of `KFORM`/`KRES`.

45. Khosla S, et al. **Relationship of serum sex steroid levels and bone turnover
    markers with bone mineral density in men and women: a key role for bioavailable
    estrogen.** J Clin Endocrinol Metab. 1998;83(7):2266-2274.
    <https://pubmed.ncbi.nlm.nih.gov/9661593/>

46. Smith EP, et al. **Estrogen resistance caused by a mutation in the
    estrogen-receptor gene in a man.** N Engl J Med. 1994;331(16):1056-1061.
    <https://pubmed.ncbi.nlm.nih.gov/8090165/>
    — the natural experiment showing that oestrogen is essential to bone in men.

47. Snyder PJ, et al. **Testosterone Treatment and Fractures in Men with
    Hypogonadism (TRAVERSE fracture substudy).** N Engl J Med. 2024;390(3):203-211.
    <https://pubmed.ncbi.nlm.nih.gov/38231624/>
    — **a result the model does not reproduce**: clinical fractures 3.50% vs 2.46% (HR 1.43).
    BMD went up and yet fractures increased. This is why `MHG_trial_ledger()` prints this discrepancy rather than hiding it,
    and it exposes the limits of the very assumption that BMD is a surrogate.

---

## 9. Body composition and muscle

48. Bhasin S, et al. **Testosterone dose-response relationships in healthy young
    men.** Am J Physiol Endocrinol Metab. 2001;281(6):E1172-E1181.
    <https://pubmed.ncbi.nlm.nih.gov/11701431/>
    — **FITTED**: the dose-response of lean body mass. `SLEAN_T`, `KL50`.

49. Bhasin S, et al. **The effects of supraphysiologic doses of testosterone on
    muscle size and strength in normal men.** N Engl J Med. 1996;335(1):1-7.
    <https://pubmed.ncbi.nlm.nih.gov/8637535/>

50. Storer TW, et al. **Effects of Testosterone Supplementation for 3 Years on
    Muscle Performance and Physical Function in Older Men.** J Clin Endocrinol
    Metab. 2017;102(2):583-593. <https://pubmed.ncbi.nlm.nih.gov/27754805/>

---

## 10. Clinical outcome trials

51. Snyder PJ, et al. **Effects of Testosterone Treatment in Older Men
    (The Testosterone Trials).** N Engl J Med. 2016;374(7):611-624.
    <https://pubmed.ncbi.nlm.nih.gov/26886521/>
    — **FITTED**: sexual function improved (PDQ-Q4 +0.58 vs +0.10), but **the primary endpoints of vitality (FACIT-F)
    and walking distance were not met**. The basis for deliberately setting the PRO weights in the model small.

52. Resnick SM, et al. **Testosterone Treatment and Cognitive Function in Older
    Men With Low Testosterone and Age-Associated Memory Impairment.** JAMA.
    2017;317(7):717-727. <https://pubmed.ncbi.nlm.nih.gov/28241356/>
    — no cognitive improvement.

53. Budoff MJ, et al. **Testosterone Treatment and Coronary Artery Plaque Volume
    in Older Men With Low Testosterone.** JAMA. 2017;317(7):708-716.
    <https://pubmed.ncbi.nlm.nih.gov/28241355/>
    — an increase in non-calcified coronary plaque volume.

54. Lincoff AM, et al. **Cardiovascular Safety of Testosterone-Replacement Therapy
    (TRAVERSE).** N Engl J Med. 2023;389(2):107-117.
    <https://pubmed.ncbi.nlm.nih.gov/37326322/>
    — **a major comparison target**: n=5,246, MACE 7.0% vs 7.3% (HR 0.96, 95% CI 0.78-1.17,
    non-inferiority met). Atrial fibrillation was nonetheless increased at 3.5% vs 2.4%, pulmonary embolism at 0.9% vs 0.5%,
    and acute kidney injury at 2.3% vs 1.5%. A deterministic model has no event process and so does not predict this
    result; `MHG_trial_ledger()` states that limitation explicitly.

55. Vigen R, et al. **Association of testosterone therapy with mortality,
    myocardial infarction, and stroke in men with low testosterone levels.**
    JAMA. 2013;310(17):1829-1836. <https://pubmed.ncbi.nlm.nih.gov/24193080/>
    — the observational study that set off the controversy preceding TRAVERSE (much criticised on methodological grounds).

---

## 11. Pathophysiology by aetiology

56. Bojesen A, Juul S, Gravholt CH. **Prenatal and postnatal prevalence of
    Klinefelter syndrome: a national registry study.** J Clin Endocrinol Metab.
    2003;88(2):622-626. <https://pubmed.ncbi.nlm.nih.gov/12574191/>

57. Groth KA, et al. **Klinefelter Syndrome — A Clinical Update.** J Clin
    Endocrinol Metab. 2013;98(1):20-30.
    <https://pubmed.ncbi.nlm.nih.gov/23118429/>

58. Boehm U, et al. **Expert consensus document: European Consensus Statement on
    congenital hypogonadotropic hypogonadism — pathogenesis, diagnosis and
    treatment.** Nat Rev Endocrinol. 2015;11(9):547-564.
    <https://pubmed.ncbi.nlm.nih.gov/26194704/>
    — the genetics of Kallmann syndrome/nIHH (ANOS1, FGFR1, PROKR2, CHD7). Cluster ⑩ of the map.

59. Daniell HW. **Hypogonadism in men consuming sustained-action oral opioids.**
    J Pain. 2002;3(5):377-384. <https://pubmed.ncbi.nlm.nih.gov/14622741/>
    — **FITTED**: `OPI_SUP`. The high prevalence in chronic opioid users.

60. Grossmann M, Matsumoto AM. **A Perspective on Middle-Aged and Older Men With
    Functional Hypogonadism: Focus on Holistic Management.** J Clin Endocrinol
    Metab. 2017;102(3):1067-1075. <https://pubmed.ncbi.nlm.nih.gov/28359097/>
    — the concept of functional (reversible) hypogonadism. Cluster ⑪ of the map.

61. Corona G, et al. **Body weight loss reverts obesity-associated
    hypogonadotropic hypogonadism: a systematic review and meta-analysis.**
    Eur J Endocrinol. 2013;168(6):829-843.
    <https://pubmed.ncbi.nlm.nih.gov/23482592/>
    — **FITTED**: `WT_LOSS`. The extent of total T recovery as a function of the degree of weight loss.

62. Rastrelli G, et al. **Testosterone Replacement Therapy for Sexual Symptoms.**
    Sex Med Rev. 2019;7(3):464-475. <https://pubmed.ncbi.nlm.nih.gov/30926459/>

63. Rochira V, et al. **Hypogonadism in the HIV-infected man.** Endocrinol Metab
    Clin North Am. 2014;43(3):709-730.
    <https://pubmed.ncbi.nlm.nih.gov/25169563/>

---

## 12. Non-androgen strategies

64. Wheeler KM, et al. **A review of the role of human chorionic gonadotropin
    (hCG) in male reproduction.** Transl Androl Urol. 2019;8(Suppl 1):S6-S17.
    <https://pubmed.ncbi.nlm.nih.gov/31143669/>

65. Rambhatla A, et al. **The Role of Estrogen Modulators in Male Hypogonadism
    and Infertility.** Rev Urol. 2016;18(2):66-72.
    <https://pubmed.ncbi.nlm.nih.gov/27601965/>
    — a comparison of the mechanisms of clomiphene, enclomiphene, and anastrozole. Cluster ⑲ of the map.

66. Kim ED, et al. **Oral enclomiphene citrate raises testosterone and preserves
    sperm counts in obese hypogonadal men, unlike topical testosterone: restoration
    instead of replacement.** BJU Int. 2016;117(4):677-685.
    <https://pubmed.ncbi.nlm.nih.gov/26496621/>
    — the title itself says what the ⑥th analysis function of this model is getting at: replacement and
    restoration are different interventions.

67. Burnett-Bowie SA, et al. **Effects of aromatase inhibition on bone mineral
    density and bone turnover in older men with low testosterone levels.**
    J Clin Endocrinol Metab. 2009;94(12):4785-4792.
    <https://pubmed.ncbi.nlm.nih.gov/19820017/>
    — **FITTED**: the comparison target of `MHG_aromatase_cost()`. Aromatase inhibition raises total T
    while lowering lumbar BMD.

68. Dwyer AA, et al. **Trial of Recombinant Follicle-Stimulating Hormone Pretreatment
    for GnRH-Induced Fertility in Patients with Congenital Hypogonadotropic
    Hypogonadism.** J Clin Endocrinol Metab. 2013;98(11):E1790-E1795.
    <https://pubmed.ncbi.nlm.nih.gov/24037890/>

69. Jayasena CN, et al. **Twice-weekly administration of kisspeptin-54 for 8 weeks
    stimulates release of reproductive hormones in women with hypothalamic
    amenorrhea.** Clin Pharmacol Ther. 2010;88(6):840-847.
    <https://pubmed.ncbi.nlm.nih.gov/20980998/>
    — the concept of restarting the axis with a kisspeptin analogue (application in men is at the research stage).

---

## 13. Prostate safety — the saturation model

70. Morgentaler A, Traish AM. **Shifting the paradigm of testosterone and prostate
    cancer: the saturation model and the limits of androgen-dependent growth.**
    Eur Urol. 2009;55(2):310-320. <https://pubmed.ncbi.nlm.nih.gov/18838208/>
    — **FITTED**: the basis for the structural choice of setting `KPSA` below `DHT0` so that saturation
    occurs below the physiological range.

71. Boyle P, et al. **Endogenous and exogenous testosterone and the risk of
    prostate cancer and increased prostate-specific antigen (PSA): a
    meta-analysis.** BJU Int. 2016;118(5):731-741.
    <https://pubmed.ncbi.nlm.nih.gov/27313122/>

---

## 14. KNDy / the GnRH pulse generator

72. Lehman MN, Coolen LM, Goodman RL. **Minireview: kisspeptin/neurokinin B/
    dynorphin (KNDy) cells of the arcuate nucleus: a central node in the control
    of gonadotropin-releasing hormone secretion.** Endocrinology.
    2010;151(8):3479-3489. <https://pubmed.ncbi.nlm.nih.gov/20501670/>
    — the structure of cluster ① of the map.

73. Pitteloud N, et al. **Inhibition of luteinizing hormone secretion by
    testosterone in men requires aromatization for its pituitary but not its
    hypothalamic effects: evidence from the tandem study of normal and
    gonadotropin-releasing hormone-deficient men.** J Clin Endocrinol Metab.
    2008;93(3):784-791. <https://pubmed.ncbi.nlm.nih.gov/18073301/>
    — **FITTED**: `WE2 = 2.2`. Direct evidence that E2 is dominant in male feedback.

74. Veldhuis JD, et al. **Testosterone blunts feedback inhibition of gonadotropin
    secretion by exogenous estradiol in healthy older but not young men.**
    J Clin Endocrinol Metab. 2008;93(3):874-880.
    <https://pubmed.ncbi.nlm.nih.gov/18160462/>

75. Keenan DM, Veldhuis JD. **Pulsatility of hypothalamo-pituitary hormones: a
    challenge in quantification.** Physiology (Bethesda). 2016;31(1):34-50.
    <https://pubmed.ncbi.nlm.nih.gov/26674550/>
    — the background to this model's choice not to simulate LH pulses explicitly but to handle them as a mean
    concentration (and the price of that choice).

---

## 15. QSP methodology

76. Baron KT, Gastonguay MR. **mrgsolve: Simulate from ODE-Based Population PK/PD
    and Systems Pharmacology Models.** <https://mrgsolve.org/>

77. Musante CJ, et al. **Quantitative Systems Pharmacology: A Case for
    Disease Models.** Clin Pharmacol Ther. 2017;101(1):24-27.
    <https://pubmed.ncbi.nlm.nih.gov/27709613/>

78. Ribba B, et al. **Model-Informed Drug Development for Quantitative Systems
    Pharmacology.** CPT Pharmacometrics Syst Pharmacol. 2017;6(8):496-498.
    <https://pubmed.ncbi.nlm.nih.gov/28571121/>

---

## Appendix: What this model does not reproduce

An honest QSP model must state where it is wrong. The following are items where reproduction was deliberately
not attempted, or was attempted and failed.

| Item | Status | Reason |
|------|------|------|
| TRAVERSE MACE (HR 0.96) | not attempted | it is an event process across 5,246 people, and a deterministic single-patient model has no corresponding structure |
| TRAVERSE fractures (HR 1.43) | **reproduction failed (deliberately exposed)** | the model reproduces the rise in BMD, so it cannot predict an increase in fractures. This is a failure of the assumption that BMD is a surrogate, not a matter of adjusting a parameter (reference 47) |
| LH pulsatility | simplified | replaced by a mean concentration. The mechanism by which pulse frequency differentially regulates LHB/FSHB transcription is in the map only and not in the ODEs (reference 75) |
| SHBG allosteric binding | simplified | the simple Vermeulen model was adopted. It can underestimate free T at low SHBG (reference 7) |
| Between-subject variability (AR CAG repeats and the like) | not implemented | the model is deterministic, and a population distribution is used approximately in one place only, `HCT_SD` |
| `EC50_HEPC` | **a calibrated value, not a measured one** | the key assumption that explains the gap in erythrocytosis between formulations. `MHG_hepcidin_sensitivity()` prints the sensitivity of the conclusion to this value |

---

> **Disclaimer** — This model and this collation of the literature are for educational and research purposes.
> They must not be used directly for clinical decision-making, for prescribing, or for regulatory submission.
