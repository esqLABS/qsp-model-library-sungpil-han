# Controlled Ovarian Stimulation — References

A list of the literature underpinning this QSP model (mechanistic map · 65-ODE mrgsolve model ·
Shiny dashboard). **Every PMID was looked up individually through NCBI E-utilities (`esearch` → `esummary`)
and confirmed to match on first author, year, journal, and title.** Candidates that could not be resolved, or
whose title did not match, were dropped from the list.

The references used directly in the numerical calibration of the model are marked with ★, together with a note of
which parameter they underpin.

---

## 1. The FSH threshold/window concept and follicle selection

The whole of cluster 2 in this model (threshold distribution · selection · dominance) stands on this concept.
It is also why there is no parameter that specifies the number of follicles.

1. ★ Brown JB. Pituitary control of ovarian function—concepts derived from gonadotrophin therapy. *Aust N Z J Obstet Gynaecol* 1978. — [PMID 278588](https://pubmed.ncbi.nlm.nih.gov/278588/) — the origin of the threshold/window concept. The model's `T50` and `SIGT`
2. ★ Baird DT. A model for follicular selection and ovulation: lessons from superovulation. *J Steroid Biochem* 1987. — [PMID 3121918](https://pubmed.ncbi.nlm.nih.gov/3121918/) — the basis for the structure in which "exogenous FSH does not create follicles, it prevents atresia"
3. ★ Fauser BC, Van Heusden AM. Manipulation of human ovarian function: physiological concepts and clinical consequences. *Endocr Rev* 1997. — [PMID 9034787](https://pubmed.ncbi.nlm.nih.gov/9034787/) — the natural/stimulated cycle symmetry of the model
4. ★ Schipper I, et al. The follicle-stimulating hormone (FSH) threshold/window concept examined by different interventions with exogenous FSH during the follicular phase of the normal menstrual cycle. *J Clin Endocrinol Metab* 1998. — [PMID 9543158](https://pubmed.ncbi.nlm.nih.gov/9543158/) — the basis for `HT = 4` (the threshold Hill coefficient)
5. ★ Hillier SG. Current concepts of the roles of follicle stimulating hormone and luteinizing hormone in folliculogenesis. *Hum Reprod* 1994. — [PMID 8027271](https://pubmed.ncbi.nlm.nih.gov/8027271/) — the switch to LH dependence beyond 11 mm (`KLHA`, `DLH`)
6. Macklon NS, Fauser BC. Aspects of ovarian follicle development throughout life. *Horm Res* 1999. — [PMID 10725781](https://pubmed.ncbi.nlm.nih.gov/10725781/)
7. Mumusoglu S, et al. Initial and cyclic recruitment of ovarian follicles: a quarter-century update. *Reprod Biomed Online* 2025. — [PMID 41067165](https://pubmed.ncbi.nlm.nih.gov/41067165/)
8. Lyu Z, et al. Stochastic mechanism of dominant follicle selection: selection of one suppresses selection of others. *J R Soc Interface* 2026. — [PMID 42014065](https://pubmed.ncbi.nlm.nih.gov/42014065/) — a mechanism of dominance different from (stochastic rather than deterministic) this model's. Recorded as an alternative hypothesis

## 2. Steroidogenesis — the two-cell, two-gonadotrophin structure

Cluster 3. The multiplicative structure theca (LH) → androgen → granulosa cell (FSH) aromatase → E2 is what
produces "suppress LH deeply and E2 per follicle falls" in the model.

9. ★ Hillier SG, Whitelaw PF, Smyth CD. Follicular oestrogen synthesis: the 'two-cell, two-gonadotrophin' model revisited. *Mol Cell Endocrinol* 1994. — [PMID 8056158](https://pubmed.ncbi.nlm.nih.gov/8056158/) — `KTHECA = 0.40 IU/L`
10. Hillier SG. Gonadotropic control of ovarian follicular growth and development. *Mol Cell Endocrinol* 2001. — [PMID 11420129](https://pubmed.ncbi.nlm.nih.gov/11420129/)
11. Magoffin DA. The ovarian androgen-producing cells: a 2001 perspective. *Rev Endocr Metab Disord* 2002. — [PMID 11883104](https://pubmed.ncbi.nlm.nih.gov/11883104/)
12. Hattori N, et al. Epoxide hydrolase affects estrogen production in the human ovary. *Endocrinology* 2000. — [PMID 10965908](https://pubmed.ncbi.nlm.nih.gov/10965908/)

## 3. AMH — both a marker of reserve and a threshold factor

AMH plays two roles in the model: it is an observable surrogate for cohort size and, at the same time,
it is the term that raises the FSH threshold (`KAMHT`).

13. ★ Weenen C, et al. Anti-Müllerian hormone expression pattern in the human ovary: potential implications for initial and cyclic follicle recruitment. *Mol Hum Reprod* 2004. — [PMID 14742691](https://pubmed.ncbi.nlm.nih.gov/14742691/) — the fact that AMH is secreted only by small follicles = the model's `DAMH = 8 mm`
14. ★ Pellatt L, et al. Granulosa cell production of anti-Müllerian hormone is increased in polycystic ovaries. *J Clin Endocrinol Metab* 2007. — [PMID 17062765](https://pubmed.ncbi.nlm.nih.gov/17062765/) — the `KAMHT` pathway of the PCOS phenotype
15. Dewailly D, et al. The physiology and clinical utility of anti-Müllerian hormone in women. *Hum Reprod Update* 2014. — [PMID 24430863](https://pubmed.ncbi.nlm.nih.gov/24430863/)
16. La Marca A, et al. Anti-Müllerian hormone (AMH) as a predictive marker in assisted reproductive technology (ART). *Hum Reprod Update* 2010. — [PMID 19793843](https://pubmed.ncbi.nlm.nih.gov/19793843/)
17. ★ Broer SL, et al. AMH and AFC as predictors of excessive response in controlled ovarian hyperstimulation: a meta-analysis. *Hum Reprod Update* 2011. — [PMID 20667894](https://pubmed.ncbi.nlm.nih.gov/20667894/) — definition of the high-responder group as AFC 25–32
18. Broekmans FJ, et al. The antral follicle count: practical recommendations for better standardization. *Fertil Steril* 2010. — [PMID 19589513](https://pubmed.ncbi.nlm.nih.gov/19589513/) — definition of `AFC` (2–9 mm)

## 4. Gonadotrophin pharmacokinetics

The basis for the rFSH PK of the model (F 0.80 · V/F 26 L · multiple-dose t½ 42 h → 11.9 IU/L at 150 IU/day)
and for the corifollitropin compartment.

19. ★ le Cotonnec JY, et al. Clinical pharmacology of recombinant human follicle-stimulating hormone (FSH). I. Comparative pharmacokinetics with urinary human FSH. *Fertil Steril* 1994. — [PMID 8150109](https://pubmed.ncbi.nlm.nih.gov/8150109/) — `KAF`·`KELFX`·`VF`
20. ★ Duijkers IJ, et al. Single dose pharmacokinetics and effects on follicular growth and serum hormones of a long-acting recombinant FSH preparation (FSH-CTP) in healthy pituitary-suppressed females. *Hum Reprod* 2002. — [PMID 12151425](https://pubmed.ncbi.nlm.nih.gov/12151425/) — corifollitropin t½ 69 h
21. Fauser BC, et al. Advances in recombinant DNA technology: corifollitropin alfa, a hybrid molecule with sustained follicle-stimulating activity and reduced injection frequency. *Hum Reprod Update* 2009. — [PMID 19182099](https://pubmed.ncbi.nlm.nih.gov/19182099/)
22. Corifollitropin Alfa Dose-finding Study Group. A randomized dose-response trial of a single injection of corifollitropin alfa to sustain multifollicular growth during controlled ovarian stimulation. *Hum Reprod* 2008. — [PMID 18684735](https://pubmed.ncbi.nlm.nih.gov/18684735/) — calibration of `POTCO` (FSH-equivalent potency)
23. Olsson H, Sandström R, Grundemar L. Different pharmacokinetic and pharmacodynamic properties of recombinant follicle-stimulating hormone (rFSH) derived from a human cell line compared with rFSH from a non-human cell line. *J Clin Pharmacol* 2014. — [PMID 24800998](https://pubmed.ncbi.nlm.nih.gov/24800998/) — follitropin delta

## 5. GnRH antagonists — pharmacokinetics and surge blockade

24. ★ Oberyé JJ, et al. Pharmacokinetic and pharmacodynamic characteristics of ganirelix (Antagon/Orgalutran). Part I. Absolute bioavailability of 0.25 mg of ganirelix after a single subcutaneous injection in healthy female volunteers. *Fertil Steril* 1999. — [PMID 10593371](https://pubmed.ncbi.nlm.nih.gov/10593371/) — `FBIOANT 0.91`, Cmax 11.2 ng/mL, tmax 1.1 h
25. ★ Oberyé JJ, et al. Pharmacokinetic and pharmacodynamic characteristics of ganirelix (Antagon/Orgalutran). Part II. Dose-proportionality and gonadotropin suppression after multiple doses of ganirelix in healthy female volunteers. *Fertil Steril* 1999. — [PMID 10593372](https://pubmed.ncbi.nlm.nih.gov/10593372/) — LH suppression of 70–80% = the model's 77% (`IC50ANT`, `FANTL`)
26. ★ Ganirelix Dose-finding Study Group. A double-blind, randomized, dose-finding study to assess the efficacy of the gonadotrophin-releasing hormone antagonist ganirelix (Org 37462) to prevent premature luteinizing hormone surges in women undergoing ovarian stimulation with recombinant follicle stimulating hormone (Puregon). *Hum Reprod* 1998. — [PMID 9853849](https://pubmed.ncbi.nlm.nih.gov/9853849/) — the basis for choosing 0.25 mg
27. Lambalk CB, et al. GnRH antagonist versus long agonist protocols in IVF: a systematic review and meta-analysis accounting for patient type. *Hum Reprod Update* 2017. — [PMID 28903472](https://pubmed.ncbi.nlm.nih.gov/28903472/)
28. ★ Kuang Y, et al. Medroxyprogesterone acetate is an effective oral alternative for preventing premature luteinizing hormone surges in women undergoing controlled ovarian hyperstimulation for in vitro fertilization. *Fertil Steril* 2015. — [PMID 25956370](https://pubmed.ncbi.nlm.nih.gov/25956370/) — PPOS arm (`PPOS`, `KP4S`)

## 6. The trigger — where this model's central claim is tested

hCG and GnRH agonists use the same receptor, but in the model the two drugs diverge because maturation (the
leading edge) and VEGF (the area) are read through different kernels.

29. ★ Humaidan P, et al. GnRH agonist (buserelin) or hCG for ovulation induction in GnRH antagonist IVF/ICSI cycles: a prospective randomized study. *Hum Reprod* 2005. — [PMID 15760966](https://pubmed.ncbi.nlm.nih.gov/15760966/) — luteal phase deficiency with an agonist trigger. The model reproduces only the direction (−41%) and underpredicts the magnitude
30. ★ Kolibianakis EM, et al. A lower ongoing pregnancy rate can be expected when GnRH agonist is used for triggering final oocyte maturation instead of HCG in patients undergoing IVF with GnRH antagonists. *Hum Reprod* 2005. — [PMID 15979994](https://pubmed.ncbi.nlm.nih.gov/15979994/)
31. ★ Fauser BC, et al. Endocrine profiles after triggering of final oocyte maturation with GnRH agonist after cotreatment with the GnRH antagonist ganirelix during ovarian hyperstimulation for in vitro fertilization. *J Clin Endocrinol Metab* 2002. — [PMID 11836309](https://pubmed.ncbi.nlm.nih.gov/11836309/) — peak LH and duration of the agonist flare (`AMPA`, `SLHMAX`)
32. ★ Youssef MA, et al. Gonadotropin-releasing hormone agonist versus HCG for oocyte triggering in antagonist-assisted reproductive technology. *Cochrane Database Syst Rev* 2014. — [PMID 25358904](https://pubmed.ncbi.nlm.nih.gov/25358904/) — meta-analysis of reduced OHSS and lower pregnancy rates after fresh transfer
33. ★ Griffin D, et al. Dual trigger of oocyte maturation with gonadotropin-releasing hormone agonist and low-dose human chorionic gonadotropin to optimize live birth rates in high responders. *Fertil Steril* 2012. — [PMID 22480822](https://pubmed.ncbi.nlm.nih.gov/22480822/) — the dual trigger arm
34. ★ Shapiro BS, et al. Comparison of "triggers" using leuprolide acetate alone or in combination with low-dose human chorionic gonadotropin. *Fertil Steril* 2011. — [PMID 21550042](https://pubmed.ncbi.nlm.nih.gov/21550042/)
35. Humaidan P, et al. GnRH agonist for triggering of final oocyte maturation: time for a change of practice? *Hum Reprod Update* 2011. — [PMID 21450755](https://pubmed.ncbi.nlm.nih.gov/21450755/)
36. Humaidan P, et al. GnRHa trigger and individualized luteal phase hCG support according to ovarian response to stimulation: two prospective randomized controlled multi-centre studies in IVF patients. *Hum Reprod* 2013. — [PMID 23753114](https://pubmed.ncbi.nlm.nih.gov/23753114/)
37. Yding Andersen C, et al. Improving the luteal phase after ovarian stimulation: reviewing new options. *Reprod Biomed Online* 2014. — [PMID 24656557](https://pubmed.ncbi.nlm.nih.gov/24656557/)
38. ★ Beckers NG, et al. Nonsupplemented luteal phase characteristics after the administration of recombinant human chorionic gonadotropin, recombinant luteinizing hormone, or gonadotropin-releasing hormone (GnRH) agonist to induce final oocyte maturation in in vitro fertilization patients after ovarian stimulation with recombinant follicle-stimulating hormone and GnRH antagonist cotreatment. *J Clin Endocrinol Metab* 2003. — [PMID 12970285](https://pubmed.ncbi.nlm.nih.gov/12970285/) — luteal phase P4 profiles. The model's `KP4CL`, `KSATCL`, and `KLYSM`

## 7. Oocyte maturation and the time of retrieval — the interval between two clocks

39. ★ Park JY, et al. EGF-like growth factors as mediators of LH action in the ovulatory follicle. *Science* 2004. — [PMID 14726596](https://pubmed.ncbi.nlm.nih.gov/14726596/) — the molecular basis for the structure in which maturation is "a leading-edge signal plus autonomous progression"
40. ★ Mansour RT, Aboulghar MA, Serour GI. Study of the optimum time for human chorionic gonadotropin–ovum pickup interval in in vitro fertilization. *J Assist Reprod Genet* 1994. — [PMID 7633170](https://pubmed.ncbi.nlm.nih.gov/7633170/) — the 34–38 h window. The model's `TMAT 0.85 d`, `TRUP 1.55 d`, and `RUPX 0.90 d`
41. Dozortsev DI, Diamond MP. Luteinizing hormone-independent rise of progesterone as the physiological trigger of the ovulatory gonadotropins surge in the human. *Fertil Steril* 2020. — [PMID 32741458](https://pubmed.ncbi.nlm.nih.gov/32741458/) — an alternative trigger hypothesis that the model does not adopt
42. Jia Z, et al. Effects of C-type natriuretic peptide on meiotic arrest and developmental competence of oocytes. *Sci Rep* 2020. — [PMID 33106527](https://pubmed.ncbi.nlm.nih.gov/33106527/) — cGMP meiotic arrest (cluster 8 of the map)

## 8. OHSS — VEGF-mediated vascular permeability

In the model OHSS is not a severity scale but is computed as the result of (granulosa cell mass × area of LHCGR
occupancy).

43. ★ McClure N, et al. Vascular endothelial growth factor as capillary permeability agent in ovarian hyperstimulation syndrome. *Lancet* 1994. — [PMID 7913160](https://pubmed.ncbi.nlm.nih.gov/7913160/) — the origin of the VEGF → permeability axis
44. ★ Gómez R, et al. Vascular endothelial growth factor receptor-2 activation induces vascular permeability in hyperstimulated rats, and this effect is prevented by dopamine receptor-2 activation. *Endocrinology* 2002. — [PMID 12399430](https://pubmed.ncbi.nlm.nih.gov/12399430/) — the site of action of cabergoline (`EMAXCAB`)
45. ★ Soares SR, et al. Targeting the vascular endothelial growth factor system to prevent ovarian hyperstimulation syndrome. *Hum Reprod Update* 2008. — [PMID 18385260](https://pubmed.ncbi.nlm.nih.gov/18385260/)
46. ★ Golan A, et al. Ovarian hyperstimulation syndrome: an update review. *Obstet Gynecol Surv* 1989. — [PMID 2660037](https://pubmed.ncbi.nlm.nih.gov/2660037/) — the grade definitions that the model computes (Hct 45%/55%)
47. ★ Practice Committee of the American Society for Reproductive Medicine. Prevention of moderate and severe ovarian hyperstimulation syndrome: a guideline. *Fertil Steril* 2024. — [PMID 38099867](https://pubmed.ncbi.nlm.nih.gov/38099867/)
48. ★ Papanikolaou EG, et al. Incidence and prediction of ovarian hyperstimulation syndrome in women undergoing gonadotropin-releasing hormone antagonist in vitro fertilization cycles. *Fertil Steril* 2006. — [PMID 16412740](https://pubmed.ncbi.nlm.nih.gov/16412740/) — risk prediction based on follicle number
49. ★ Steward RG, et al. Oocyte number as a predictor for ovarian hyperstimulation syndrome and live birth: an analysis of 256,381 in vitro fertilization cycles. *Fertil Steril* 2014. — [PMID 24462057](https://pubmed.ncbi.nlm.nih.gov/24462057/) — the simultaneous relationship between oocyte number, OHSS, and birth rate. Corresponds to sweep D of the model
50. ★ Devroey P, Polyzos NP, Blockeel C. An OHSS-Free Clinic by segmentation of IVF treatment. *Hum Reprod* 2011. — [PMID 21828116](https://pubmed.ncbi.nlm.nih.gov/21828116/) — the agonist trigger + freeze-all arm
51. ★ Alvarez C, et al. Dopamine agonist cabergoline reduces hemoconcentration and ascites in hyperstimulated women undergoing assisted reproduction. *J Clin Endocrinol Metab* 2007. — [PMID 17456571](https://pubmed.ncbi.nlm.nih.gov/17456571/) — the calibration target for the cabergoline arm
52. Korhonen KV, et al. C-reactive protein response is higher in early than in late ovarian hyperstimulation syndrome. *Eur J Obstet Gynecol Reprod Biol* 2016. — [PMID 27865939](https://pubmed.ncbi.nlm.nih.gov/27865939/) — the early/late distinction
53. Youssef MA, et al. Volume expanders for the prevention of ovarian hyperstimulation syndrome. *Cochrane Database Syst Rev* 2016. — [PMID 27577848](https://pubmed.ncbi.nlm.nih.gov/27577848/)
54. Bassiouny YA, et al. Randomized trial of combined cabergoline and coasting in preventing ovarian hyperstimulation syndrome. *Int J Gynaecol Obstet* 2018. — [PMID 29055130](https://pubmed.ncbi.nlm.nih.gov/29055130/) — coasting arm
55. Wu D, et al. Comparison of the effectiveness of various medicines in the prevention of ovarian hyperstimulation syndrome. *Front Endocrinol (Lausanne)* 2022. — [PMID 35154015](https://pubmed.ncbi.nlm.nih.gov/35154015/)
56. Salari E, et al. Comparative study of cabergoline and hydroxychloroquine to prevent ovarian hyperstimulation syndrome. *J Ovarian Res* 2025. — [PMID 40442814](https://pubmed.ncbi.nlm.nih.gov/40442814/)
57. Timmons D, et al. Ovarian hyperstimulation syndrome: a review for emergency clinicians. *Am J Emerg Med* 2019. — [PMID 31097257](https://pubmed.ncbi.nlm.nih.gov/31097257/) — severe complications (thrombosis, ARDS)
58. Braat DD, et al. Maternal death related to IVF in the Netherlands 1984–2008. *Hum Reprod* 2010. — [PMID 20488805](https://pubmed.ncbi.nlm.nih.gov/20488805/)

## 9. Progesterone elevation on the day of trigger, and freeze-all

In the model, P4 on the day of trigger is not a separate event but a fourth reading of granulosa cell mass.

59. ★ Venetis CA, et al. Progesterone elevation and probability of pregnancy after IVF: a systematic review and meta-analysis of over 60,000 cycles. *Hum Reprod Update* 2013. — [PMID 23827986](https://pubmed.ncbi.nlm.nih.gov/23827986/) — the 1.5 ng/mL threshold
60. ★ Bosch E, et al. Circulating progesterone levels and ongoing pregnancy rates in controlled ovarian stimulation cycles for in vitro fertilization: analysis of over 4000 cycles. *Hum Reprod* 2010. — [PMID 20539042](https://pubmed.ncbi.nlm.nih.gov/20539042/)
61. Labarta E, et al. Endometrial receptivity is affected in women with high circulating progesterone levels at the end of the follicular phase: a functional genomics analysis. *Hum Reprod* 2011. — [PMID 21540246](https://pubmed.ncbi.nlm.nih.gov/21540246/)
62. Shapiro BS, et al. Evidence of impaired endometrial receptivity after ovarian stimulation for in vitro fertilization: a prospective randomized trial. *Fertil Steril* 2011. — [PMID 21737072](https://pubmed.ncbi.nlm.nih.gov/21737072/)
63. Roque M, et al. Fresh versus elective frozen embryo transfer in IVF/ICSI cycles: a systematic review and meta-analysis of reproductive outcomes. *Hum Reprod Update* 2019. — [PMID 30388233](https://pubmed.ncbi.nlm.nih.gov/30388233/)

## 10. Dose-response and individualised dosing

The clinical literature to which sweep A of the model (a 3-fold dose → oocytes +15%, ascites +125%) corresponds.

64. ★ Nyboe Andersen A, et al. Individualized versus conventional ovarian stimulation for in vitro fertilization: a multicenter, randomized, controlled, assessor-blinded, phase 3 noninferiority trial (ESTHER-1). *Fertil Steril* 2017. — [PMID 27912901](https://pubmed.ncbi.nlm.nih.gov/27912901/) — the calibration target for the standard arm (AFC 12 · 150 IU · ~10 oocytes)
65. Nelson SM, et al. Individualized versus conventional ovarian stimulation for in vitro fertilization (authors' reply). *Fertil Steril* 2024. — [PMID 38266801](https://pubmed.ncbi.nlm.nih.gov/38266801/)
66. Ngwenya O, et al. Individualised gonadotropin dose selection using markers of ovarian reserve for women undergoing IVF/ICSI. *Cochrane Database Syst Rev* 2024. — [PMID 38174816](https://pubmed.ncbi.nlm.nih.gov/38174816/)
67. van Tilborg TC, et al. Individualized versus standard FSH dosing in women starting IVF/ICSI: an RCT. Part 1: The predicted poor responder (OPTIMIST). *Hum Reprod* 2017. — [PMID 29121326](https://pubmed.ncbi.nlm.nih.gov/29121326/) — the result that raising the dose does not change the outcome
68. Arce JC, et al. Antimüllerian hormone in gonadotropin releasing-hormone antagonist cycles: prediction of ovarian response and cumulative treatment outcome in good-prognosis patients. *Fertil Steril* 2013. — [PMID 23394782](https://pubmed.ncbi.nlm.nih.gov/23394782/)
69. Blockeel C, et al. Follicular phase endocrine characteristics during ovarian stimulation and GnRH antagonist cotreatment for IVF. *J Clin Endocrinol Metab* 2011. — [PMID 21307142](https://pubmed.ncbi.nlm.nih.gov/21307142/) — trajectories of endogenous FSH, LH, and E2 during stimulation

## 11. Oocyte number · ploidy · cumulative live birth rate

The coefficients of the embryo chain in the model (2PN → blastocyst → euploid → CLBR) all came from here.

70. ★ Sunkara SK, et al. Association between the number of eggs and live birth in IVF treatment: an analysis of 400,135 treatment cycles. *Hum Reprod* 2011. — [PMID 21558332](https://pubmed.ncbi.nlm.nih.gov/21558332/) — the saturating curve of oocyte number against birth rate
71. ★ Drakopoulos P, et al. Conventional ovarian stimulation and single embryo transfer for IVF/ICSI. How many oocytes do we need to maximize cumulative live birth rates after utilization of all fresh and frozen embryos? *Hum Reprod* 2016. — [PMID 26724797](https://pubmed.ncbi.nlm.nih.gov/26724797/) — the reference corresponding to sweep D of the model
72. ★ Franasiak JM, et al. The nature of aneuploidy with increasing age of the female partner: a review of 15,169 consecutive trophectoderm biopsies evaluated with comprehensive chromosomal screening. *Fertil Steril* 2014. — [PMID 24355045](https://pubmed.ncbi.nlm.nih.gov/24355045/) — the euploid fraction curve `0.85/(1+exp((AGE−38.2)/4))`
73. Goldman RH, et al. Predicting the likelihood of live birth for elective oocyte cryopreservation: a counseling tool for physicians and patients. *Hum Reprod* 2017. — [PMID 28166330](https://pubmed.ncbi.nlm.nih.gov/28166330/) — an independent construction of the oocyte → birth chain
74. Toftager M, et al. Cumulative live birth rates after one ART cycle including all subsequent frozen-thaw cycles in 1050 women. *Hum Reprod* 2017. — [PMID 28130435](https://pubmed.ncbi.nlm.nih.gov/28130435/)
75. Totonchi M, et al. Preimplantation genetic screening and the success rate of in vitro fertilization: a three-years study on Iranian population. *Cell J* 2021. — [PMID 32347040](https://pubmed.ncbi.nlm.nih.gov/32347040/)

## 12. PCOS — the mechanism of the high-responder group

The model does not designate PCOS as a "high responder". Two terms alone — AMH raising the threshold, and LH tone
arresting follicles before they acquire LHCGR — produce both the 6–9 mm arrest (the polycystic morphology) and the
excessive response to stimulation.

76. ★ Jonard S, Dewailly D. The follicular excess in polycystic ovaries, due to intra-ovarian hyperandrogenism, may be the main culprit for the follicular arrest. *Hum Reprod Update* 2004. — [PMID 15073141](https://pubmed.ncbi.nlm.nih.gov/15073141/) — the mechanism of follicular arrest
77. ★ Willis DS, et al. Premature response to luteinizing hormone of granulosa cells from anovulatory women with polycystic ovary syndrome: relevance to mechanism of anovulation. *J Clin Endocrinol Metab* 1998. — [PMID 9814480](https://pubmed.ncbi.nlm.nih.gov/9814480/) — the model's `KLHARR` and `NLHARR`

## 13. Poor responders · LH supplementation

78. ★ Ferraretti AP, et al. ESHRE consensus on the definition of 'poor response' to ovarian stimulation for in vitro fertilization: the Bologna criteria. *Hum Reprod* 2011. — [PMID 21505041](https://pubmed.ncbi.nlm.nih.gov/21505041/)
79. ★ Humaidan P, et al. The novel POSEIDON stratification of 'Low prognosis patients in Assisted Reproductive Technology' and its proposed marker of successful outcome. *F1000Res* 2016. — [PMID 28232864](https://pubmed.ncbi.nlm.nih.gov/28232864/) — the poor-responder arm (AFC 5 · T50 12)
80. ★ Mochtar MH, et al. Recombinant luteinizing hormone (rLH) and recombinant follicle stimulating hormone (rFSH) for ovarian stimulation in IVF/ICSI cycles. *Cochrane Database Syst Rev* 2017. — [PMID 28537052](https://pubmed.ncbi.nlm.nih.gov/28537052/) — the reference corresponding to the model's prediction that LH supplementation fails to increase oocyte number

## 14. Luteal phase support

81. Fatemi HM, et al. An update of luteal phase support in stimulated IVF cycles. *Hum Reprod Update* 2007. — [PMID 17626114](https://pubmed.ncbi.nlm.nih.gov/17626114/)
82. Fatemi HM. Simplifying luteal phase support in stimulated assisted reproduction cycles. *Fertil Steril* 2018. — [PMID 30396544](https://pubmed.ncbi.nlm.nih.gov/30396544/)
83. ★ Tavaniotou A, Devroey P. Effect of human chorionic gonadotropin on luteal luteinizing hormone concentrations in natural cycles. *Fertil Steril* 2003. — [PMID 12969719](https://pubmed.ncbi.nlm.nih.gov/12969719/) — the direct basis for progesterone–LH negative feedback (`KP4LH`)

## 15. Mathematical models of the menstrual cycle and of follicular dynamics (prior models)

This model stands in the line of the family below; the differences are (i) that follicles are discretised into
threshold-quantile slots, (ii) that the LHCGR signal is split across three kernels, and (iii) that an OHSS fluid
compartment is included.

84. ★ Reinecke I, Deuflhard P. A complex mathematical model of the human menstrual cycle. *J Theor Biol* 2007. — [PMID 17448501](https://pubmed.ncbi.nlm.nih.gov/17448501/)
85. ★ Röblitz S, et al. A mathematical model of the human menstrual cycle for the administration of GnRH analogues. *J Theor Biol* 2013. — [PMID 23206386](https://pubmed.ncbi.nlm.nih.gov/23206386/) — a prior model that includes the administration of GnRH analogues
86. Röblitz S, et al. A computational systems biology view on the role of the menstrual cycle. *J Mol Endocrinol* 2026. — [PMID 42306999](https://pubmed.ncbi.nlm.nih.gov/42306999/)
87. ★ Clark LR, Schlosser PM, Selgrade JF. Multiple stable periodic solutions in a model for hormonal control of the menstrual cycle. *Bull Math Biol* 2003. — [PMID 12597121](https://pubmed.ncbi.nlm.nih.gov/12597121/) — bistability of the cycle model

---

## Citation summary

| Category | Papers |
|------|------|
| ★ Used directly in parameter calibration | 48 |
| Background · mechanism · verification | 39 |
| **Total (every PMID individually confirmed)** | **87** |

## Where the model disagrees with the literature (deliberately left in)

Recorded here honestly. What follows are predictions produced by the structure of the model, left in a
falsifiable form.

1. **Luteal letrozole**: in the model it lowers E2 by 79% but changes ascites by less than 4%.
   This follows directly from a structure that treats E2 as a marker and VEGF as the mediator. Should randomised
   data establish that letrozole reduces OHSS, the omission of an E2-dependent permeability term is falsified.
2. **Spontaneous ovulation after an agonist trigger**: because the rupture integral (`RUPX`) is never reached, the model
   predicts almost no ovulation even at 44 hours. This implies that a more forgiving retrieval time than with hCG
   is permissible, and it has not been confirmed clinically.
3. **coasting**: because it works only by killing granulosa cell mass in small follicles, OHSS cannot be reduced
   without a loss of oocytes (−35%).
4. **Luteal phase deficiency after an agonist trigger**: the direction is reproduced but the magnitude (−41%) is smaller than in the clinic
   (references 29 and 30). The missing mechanism is presumed to be a qualitative defect in luteinisation itself.
5. **Cabergoline**: more conservative than the RR of 0.38 in the meta-analysis (ascites −17%).
6. **Dominant follicle diameter in the natural cycle**: it stops at 16.9 mm (observed 20–22 mm). This is because the E2
   concentration that fires the surge is reached before the diameter is.
