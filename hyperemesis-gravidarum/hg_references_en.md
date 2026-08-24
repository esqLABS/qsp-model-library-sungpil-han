# Hyperemesis Gravidarum (HG) — QSP Model References

The reference base for the hyperemesis gravidarum QSP model.
Every PMID was cross-checked against author, journal, year, and title through NCBI E-utilities (2026-07-29).

Ordering principle: the literature that can support or refute the **structural claims** of the model is placed first.
The preamble to each section states which equation or parameter of the model that group of references went into.

---

## 1. The core mechanism — the GDF15·GFRAL axis and the "fold-change" hypothesis

> Structural assumption 1 of the model (`SP`·`ALPHA`·`TAU_SP`, cluster 3) comes from here.
> Fejzo 2024 is decisive: the **amount** of fetally derived GDF15 produced and the mother's **susceptibility**
> each contribute independently to risk, susceptibility is set by pre-pregnancy GDF15 exposure, and
> in the mouse, prior exposure to GDF15 induces **desensitisation**.

| # | Reference | PMID |
|---|------|------|
| 1 | Fejzo M, et al. **GDF15 linked to maternal risk of nausea and vomiting during pregnancy.** *Nature* 2024;625:760-767. | [38092039](https://pubmed.ncbi.nlm.nih.gov/38092039/) |
| 2 | Fejzo MS, et al. **Placenta and appetite genes GDF15 and IGFBP7 are associated with hyperemesis gravidarum.** *Nat Commun* 2018. | [29563502](https://pubmed.ncbi.nlm.nih.gov/29563502/) |
| 3 | Fejzo M, et al. **Multi-ancestry genome-wide association study of severe pregnancy nausea and vomiting.** *Nat Genet* 2026. | [41981316](https://pubmed.ncbi.nlm.nih.gov/41981316/) |
| 4 | Mullican SE, et al. **GFRAL is the receptor for GDF15 and the ligand promotes weight loss in mice and nonhuman primates.** *Nat Med* 2017. | [28846097](https://pubmed.ncbi.nlm.nih.gov/28846097/) |
| 5 | Breit SN, et al. **The GDF15-GFRAL Pathway in Health and Metabolic Disease: Friend or Foe?** *Annu Rev Physiol* 2021. | [33228454](https://pubmed.ncbi.nlm.nih.gov/33228454/) |
| 6 | Hes C, et al. **GDNF family receptor alpha-like (GFRAL) expression is restricted to the caudal brainstem.** *Mol Metab* 2025. | [39608751](https://pubmed.ncbi.nlm.nih.gov/39608751/) |
| 7 | Coll AP, et al. **GDF15 mediates the effects of metformin on body weight and energy balance.** *Nature* 2020. | [31875646](https://pubmed.ncbi.nlm.nih.gov/31875646/) |
| 8 | Alshaikh ABA, et al. **Hyperemesis gravidarum revisited: from GDF15 biology to precision multimodal therapy.** *Naunyn Schmiedebergs Arch Pharmacol* 2026. | [41942591](https://pubmed.ncbi.nlm.nih.gov/41942591/) |

**Where it entered the model.** References 4 and 6, on the basis that GFRAL is expressed only in the area postrema/NTS,
restrict the output of the fold-change detector so that it enters the hindbrain node alone (cluster 4).
Reference 7 is the sole basis for the metformin → raised baseline GDF15 pathway (`MET_EMAX_GDF`).

---

## 2. The time course of GDF15 in pregnancy — "the concentration keeps rising and the symptoms disappear"

> The most awkward observation the model has to explain. Maternal GDF15 does not fall after the first
> trimester (it rises, if anything, to term), yet symptoms peak at 9-11 weeks of gestation and
> disappear by 16-20 weeks. A concentration-based model cannot do this.

| # | Reference | PMID |
|---|------|------|
| 9 | Marjono AB, et al. **Macrophage inhibitory cytokine-1 in gestational tissues and maternal serum in normal and pre-eclamptic pregnancy.** *Placenta* 2003. | [12495665](https://pubmed.ncbi.nlm.nih.gov/12495665/) |
| 10 | Tong S, et al. **Serum concentrations of macrophage inhibitory cytokine 1 (MIC 1) as a predictor of miscarriage.** *Lancet* 2004. | [14726168](https://pubmed.ncbi.nlm.nih.gov/14726168/) |
| 11 | Kaitu'u-Lino TJ, et al. **Plasma MIC-1 and PAPP-A levels are decreased among women presenting to an early pregnancy assessment unit...** *PLoS One* 2013. | [24069146](https://pubmed.ncbi.nlm.nih.gov/24069146/) |
| 12 | Chen Q, et al. **Serum levels of GDF15 are reduced in preeclampsia...** *Cytokine* 2016. | [27173615](https://pubmed.ncbi.nlm.nih.gov/27173615/) |
| 13 | Lyu C, et al. **Insufficient GDF15 expression predisposes women to unexplained recurrent pregnancy loss by impairing extravillous trophoblast invasion.** *Cell Prolif* 2023. | [37272232](https://pubmed.ncbi.nlm.nih.gov/37272232/) |
| 14 | Yang SH, et al. **GDF15 promotes trophoblast invasion and pregnancy success via the BMPR1A/BMPR2/p-SMAD1 pathway.** *Life Sci* 2025. | [40157640](https://pubmed.ncbi.nlm.nih.gov/40157640/) |

**Caution.** References 13 and 14 are evidence that GDF15 is **required** for trophoblast invasion, and this is
the basis for the fetal safety warning (the `FETAL_SAF` node) attached to the anti-GDF15 antibody scenario
of the model. On fold-change logic alone, blocking the ligand is the best move; but that ligand is also
something the placenta needs in the first trimester.

---

## 3. Three independent lines of evidence that pre-pregnancy GDF15 exposure is protective

> Where the strongest prediction of the model is tested. β-thalassaemia (a chronic elevation),
> metformin (about a 2-fold rise), smoking (about a 1/3 rise) — all three exposures act **before**
> pregnancy, all three are protective, and the strength of protection lines up in the order of the size of the rise.
> The model reproduces that ordering from `ALPHA` alone (a prediction, not a fit).

| # | Reference | PMID |
|---|------|------|
| 15 | Sharma N, et al. **Prepregnancy metformin use associated with lower risk of severe nausea and vomiting of pregnancy and hyperemesis gravidarum.** *Am J Obstet Gynecol* 2025. — metformin aRR 0.29 (95% CI 0.12-0.71), smoking aRR 0.51 (0.30-0.86), n=5414 | [40588059](https://pubmed.ncbi.nlm.nih.gov/40588059/) |
| 16 | Ranjbaran R, et al. **GDF-15 negatively regulates excess erythropoiesis and its overexpression is involved in erythroid hyperplasia.** *Exp Cell Res* 2020. | [33164866](https://pubmed.ncbi.nlm.nih.gov/33164866/) |
| 17 | Piolatto A, et al. **GDF15 as a Marker of Ineffective Erythropoiesis and Erythroid Expansion in Thalassemia: a Clinical Perspective.** *Clin Lab* 2026. | [42159121](https://pubmed.ncbi.nlm.nih.gov/42159121/) |

---

## 4. Definition · epidemiology · natural history — PUQE and the Windsor criteria

> The band boundaries of the `puqe24()` function, the HG case definition (`hg_case()`), and the natural history fitting
> targets (peak at 9-11 weeks, resolution by 16 weeks) come from here.

| # | Reference | PMID |
|---|------|------|
| 18 | Fejzo MS, et al. **Nausea and vomiting of pregnancy and hyperemesis gravidarum.** *Nat Rev Dis Primers* 2019. | [31515515](https://pubmed.ncbi.nlm.nih.gov/31515515/) |
| 19 | Jansen LAW, et al. **The Windsor definition for hyperemesis gravidarum: a multistakeholder international consensus definition.** *Eur J Obstet Gynecol Reprod Biol* 2021. | [34555550](https://pubmed.ncbi.nlm.nih.gov/34555550/) |
| 20 | Koren G, et al. **Validation studies of the Pregnancy Unique-Quantification of Emesis (PUQE) scores.** *J Obstet Gynaecol* 2005. | [16147725](https://pubmed.ncbi.nlm.nih.gov/16147725/) |
| 21 | Koren G, et al. **Measuring the severity of nausea and vomiting of pregnancy; a 20-year perspective on the use of PUQE.** *J Obstet Gynaecol* 2021. | [32811235](https://pubmed.ncbi.nlm.nih.gov/32811235/) |
| 22 | Yilmaz T, et al. **Psychometric properties of the Pregnancy-Unique Quantification of Emesis (PUQE-24) Scale.** *J Obstet Gynaecol* 2022. | [35253594](https://pubmed.ncbi.nlm.nih.gov/35253594/) |
| 23 | Jansen LAW, et al. **Diagnosis and treatment of hyperemesis gravidarum.** *CMAJ* 2024. | [38621783](https://pubmed.ncbi.nlm.nih.gov/38621783/) |
| 24 | Ioannidou P, et al. **Predictive factors of Hyperemesis Gravidarum: a systematic review.** *Eur J Obstet Gynecol Reprod Biol* 2019. | [31126753](https://pubmed.ncbi.nlm.nih.gov/31126753/) |
| 25 | Beyene GA, et al. **Prevalence and determinants of hyperemesis gravidarum among pregnant women in Ethiopia: systematic review and meta-analysis.** *PLoS One* 2024. | [39625915](https://pubmed.ncbi.nlm.nih.gov/39625915/) |
| 26 | Niemeijer MN, et al. **Diagnostic markers for hyperemesis gravidarum: a systematic review and metaanalysis.** *Am J Obstet Gynecol* 2014. | [24530975](https://pubmed.ncbi.nlm.nih.gov/24530975/) |
| 27 | Vadakekut ES, et al. **Hyperemesis Gravidarum.** *StatPearls* (updated 2026). | [30422512](https://pubmed.ncbi.nlm.nih.gov/30422512/) |
| 28 | Nurmi M, et al. **Recurrence patterns of hyperemesis gravidarum.** *Am J Obstet Gynecol* 2018. | [30121224](https://pubmed.ncbi.nlm.nih.gov/30121224/) |
| 29 | Fassett MJ, et al. **Hyperemesis Gravidarum: Risk of Recurrence in Subsequent Pregnancies.** *Reprod Sci* 2023. | [36163577](https://pubmed.ncbi.nlm.nih.gov/36163577/) |
| 30 | Fejzo MS, et al. **Recurrence risk of hyperemesis gravidarum.** *J Midwifery Womens Health* 2011. | [21429077](https://pubmed.ncbi.nlm.nih.gov/21429077/) |

---

## 5. Clinical guidelines

| # | Reference | PMID |
|---|------|------|
| 31 | Nelson-Piercy C, et al. **The Management of Nausea and Vomiting in Pregnancy and Hyperemesis Gravidarum (RCOG Green-top Guideline No. 69).** *BJOG* 2024. | [38311315](https://pubmed.ncbi.nlm.nih.gov/38311315/) |
| 32 | Committee on Practice Bulletins—Obstetrics. **ACOG Practice Bulletin No. 189: Nausea and Vomiting of Pregnancy.** *Obstet Gynecol* 2018. | [29266076](https://pubmed.ncbi.nlm.nih.gov/29266076/) |
| 33 | Clark SM, et al. **Inpatient Management of Hyperemesis Gravidarum.** *Obstet Gynecol* 2024. | [38301258](https://pubmed.ncbi.nlm.nih.gov/38301258/) |
| 34 | Spinosa D, et al. **Management Considerations for Recalcitrant Hyperemesis.** *Obstet Gynecol Surv* 2020. | [31999353](https://pubmed.ncbi.nlm.nih.gov/31999353/) |

---

## 6. Randomised controlled trials — the quantitative benchmarks for validating the model

> The numbers in this table go straight into the VALIDATION table of `hg_reference_impl.py`.
> **Reference 35 (the VOMIT trial) is the starting point of the whole model**: ondansetron, the drug the
> guidelines recommend, managed only −0.51 (non-significant) while mirtazapine gave −1.86, and the difference
> opened up after day 4. Only two parameters (`W_VAG`, `E0`) were fitted to this trial;
> the remaining trials (36-39) are predictions.

| # | Reference | Key figures | PMID |
|---|------|-----------|------|
| 35 | Ostenfeld A, et al. **Mirtazapine or ondansetron for hyperemesis gravidarum: a randomized placebo-controlled trial (VOMIT).** *Am J Obstet Gynecol* 2026. | ΔPUQE-24 day 2: mirtazapine −1.86 (95% CI −3.61 to −0.12); ondansetron −0.51 (−2.32 to 1.30); n=59, 7 hospitals, Denmark | [41478546](https://pubmed.ncbi.nlm.nih.gov/41478546/) |
| 36 | Koren G, et al. **Effectiveness of delayed-release doxylamine and pyridoxine for nausea and vomiting of pregnancy: a randomized placebo controlled trial.** *Am J Obstet Gynecol* 2010. | ΔPUQE −4.8 vs −3.9 (P=.006), moderate NVP | [20843504](https://pubmed.ncbi.nlm.nih.gov/20843504/) |
| 37 | Guttuso T Jr, et al. **Effect of gabapentin on hyperemesis gravidarum: a double-blind, randomized controlled trial.** *Am J Obstet Gynecol MFM* 2021. | 52% greater reduction in PUQE on days 5-7 (95% CI 16-88), n=21 | [33451591](https://pubmed.ncbi.nlm.nih.gov/33451591/) |
| 38 | Maina A, et al. **Transdermal clonidine in the treatment of severe hyperemesis: a pilot randomised control trial (CLONEMESI).** *BJOG* 2014. | PUQE improvement 95% CI 0.43-3.24; systolic blood pressure −6 mmHg | [24684734](https://pubmed.ncbi.nlm.nih.gov/24684734/) |
| 39 | Yost NP, et al. **A randomized, placebo-controlled trial of corticosteroids for hyperemesis due to pregnancy.** *Obstet Gynecol* 2003. | readmission 34% vs 35% (P=.89) — **a negative result**, n=110 | [14662211](https://pubmed.ncbi.nlm.nih.gov/14662211/) |
| 40 | Smith C, et al. **A randomized controlled trial of ginger to treat nausea and vomiting in pregnancy.** *Obstet Gynecol* 2004. | | [15051552](https://pubmed.ncbi.nlm.nih.gov/15051552/) |
| 41 | McParlin C, et al. **Treatments for Hyperemesis Gravidarum and Nausea and Vomiting in Pregnancy: a Systematic Review.** *JAMA* 2016. | | [27701665](https://pubmed.ncbi.nlm.nih.gov/27701665/) |
| 42 | O'Donnell A, et al. **Treatments for hyperemesis gravidarum and nausea and vomiting in pregnancy: a systematic review and economic assessment.** *Health Technol Assess* 2016. | | [27731292](https://pubmed.ncbi.nlm.nih.gov/27731292/) |
| 43 | Abramowitz A, et al. **Treatment options for hyperemesis gravidarum.** *Arch Womens Ment Health* 2017. | | [28070660](https://pubmed.ncbi.nlm.nih.gov/28070660/) |

---

## 7. Pharmacokinetics · receptor pharmacology — the inputs to the node-location rule

> Structural assumption 2 of the model. The Ki, protein binding, brain penetration ratio, and PK of each drug were all taken
> from the literature below and were not fitted. The effect sizes for doxylamine, gabapentin,
> clonidine, and promethazine are therefore **predictions**.

| # | Reference | PMID |
|---|------|------|
| 44 | Elkomy MH, et al. **Ondansetron pharmacokinetics in pregnant women and neonates.** *Clin Pharmacol Ther* 2015. | [25670522](https://pubmed.ncbi.nlm.nih.gov/25670522/) |
| 45 | Siu SS, et al. **Placental transfer of ondansetron during early human pregnancy.** *Clin Pharmacokinet* 2006. | [16584287](https://pubmed.ncbi.nlm.nih.gov/16584287/) |
| 46 | Lemon LS, et al. **Ondansetron Exposure Changes in a Pregnant Woman.** *Pharmacotherapy* 2016. | [27374186](https://pubmed.ncbi.nlm.nih.gov/27374186/) |
| 47 | Davis R, Whittington R. **Mirtazapine: a review of its pharmacology and therapeutic potential.** *CNS Drugs* 1996. | [26071050](https://pubmed.ncbi.nlm.nih.gov/26071050/) |
| 48 | Andrews PL, et al. **The abdominal visceral innervation and the emetic reflex: pathways, pharmacology, and plasticity.** *Can J Physiol Pharmacol* 1990. | [2178756](https://pubmed.ncbi.nlm.nih.gov/2178756/) |
| 49 | Rudd JA, et al. **The involvement of TRPV1 in emesis and anti-emesis.** *Temperature* 2015. | [27227028](https://pubmed.ncbi.nlm.nih.gov/27227028/) |
| 50 | Wang XF, et al. **A meta-analysis of olanzapine for the prevention of chemotherapy-induced nausea and vomiting.** *Sci Rep* 2014. | [24770591](https://pubmed.ncbi.nlm.nih.gov/24770591/) |
| 51 | Langley-DeGroot M, et al. **Olanzapine in the treatment of refractory nausea and vomiting.** *J Pain Palliat Care Pharmacother* 2015. | [26095486](https://pubmed.ncbi.nlm.nih.gov/26095486/) |
| 52 | Grimison P, et al. **Oral Cannabis Extract for Secondary Prevention of Chemotherapy-Induced Nausea and Vomiting.** *J Clin Oncol* 2024. | [39151115](https://pubmed.ncbi.nlm.nih.gov/39151115/) |

---

## 8. Drug safety

| # | Reference | PMID |
|---|------|------|
| 53 | Huybrechts KF, et al. **Association of Maternal First-Trimester Ondansetron Use With Cardiac Malformations and Oral Clefts in Offspring.** *JAMA* 2018. | [30561479](https://pubmed.ncbi.nlm.nih.gov/30561479/) |
| 54 | Picot C, et al. **Risk of malformation after ondansetron in pregnancy: an updated systematic review and meta-analysis.** *Birth Defects Res* 2020. | [32420702](https://pubmed.ncbi.nlm.nih.gov/32420702/) |

---

## 9. Fluid and electrolytes — hypochloraemic metabolic alkalosis

> Half of structural assumption 3 of the model. Because chloride deficiency blocks renal bicarbonate excretion,
> the alkalosis is not corrected by anything else you give unless chloride is replaced
> (`cl_avail` gating). That bicarbonaturia forces sodium excretion is the mechanism of the
> hyponatraemia (`NA_HCO3_COUP`).

| # | Reference | PMID |
|---|------|------|
| 55 | Signorelli GC, et al. **Dietary Chloride Deficiency Syndrome: Pathophysiology, History, and Systematic Literature Review.** *Nutrients* 2020. | [33182508](https://pubmed.ncbi.nlm.nih.gov/33182508/) |
| 56 | Adrogué HJ, Madias NE. **Diagnosis and Management of Hyponatremia: a Review.** *JAMA* 2022. | [35852524](https://pubmed.ncbi.nlm.nih.gov/35852524/) |

---

## 10. Thiamine and Wernicke encephalopathy — the slowest clock

> The other half of structural assumption 3 of the model, and the most clinically actionable
> prediction: **glucose-containing fluid without thiamine carries a higher risk of Wernicke's than giving no
> fluid at all.** In the model P(WE) 15.7% vs 0.0%.

| # | Reference | PMID |
|---|------|------|
| 57 | Oudman E, et al. **Wernicke's encephalopathy in hyperemesis gravidarum: a systematic review.** *Eur J Obstet Gynecol Reprod Biol* 2019. | [30889425](https://pubmed.ncbi.nlm.nih.gov/30889425/) |
| 58 | Erick M. **Gestational malnutrition, hyperemesis gravidarum, and Wernicke's encephalopathy: what is missing?** *Nutr Clin Pract* 2022. | [36250744](https://pubmed.ncbi.nlm.nih.gov/36250744/) |
| 59 | Maslin K, Dean C. **Nutritional consequences and management of hyperemesis gravidarum: a narrative review.** *Nutr Res Rev* 2022. | [34526158](https://pubmed.ncbi.nlm.nih.gov/34526158/) |
| 60 | Stokke G, et al. **Hyperemesis gravidarum, nutritional treatment by nasogastric tube feeding: a 10-year retrospective cohort study.** *Acta Obstet Gynecol Scand* 2015. | [25581215](https://pubmed.ncbi.nlm.nih.gov/25581215/) |
| 61 | Hsu JJ, et al. **Nasogastric enteral feeding in the management of hyperemesis gravidarum.** *Obstet Gynecol* 1996. | [8752236](https://pubmed.ncbi.nlm.nih.gov/8752236/) |

---

## 11. Endocrine — cross-stimulation of the TSH receptor by hCG

> Because hCG and GDF15 come from the same cell (the syncytiotrophoblast), the model predicts that biochemical
> thyrotoxicosis follows the **fetal production axis** (`TROPH_GAIN`) and not the maternal susceptibility
> axis (`SENS`). That means two women with the same PUQE can have entirely different TSH, and it
> is testable with paired GDF15/TSH measurements.

| # | Reference | PMID |
|---|------|------|
| 62 | Rodien P, et al. **Abnormal stimulation of the thyrotrophin receptor during gestation.** *Hum Reprod Update* 2004. | [15073140](https://pubmed.ncbi.nlm.nih.gov/15073140/) |
| 63 | Albaar MT, Adam JMF. **Gestational transient thyrotoxicosis.** *Acta Med Indones* 2009. | [19390130](https://pubmed.ncbi.nlm.nih.gov/19390130/) |
| 64 | Iijima S. **Pitfalls in the assessment of gestational transient thyrotoxicosis.** *Gynecol Endocrinol* 2020. | [32301638](https://pubmed.ncbi.nlm.nih.gov/32301638/) |
| 65 | Zimmerman CF, et al. **Thyroid Storm Caused by Hyperemesis Gravidarum.** *AACE Clin Case Rep* 2022. | [35602873](https://pubmed.ncbi.nlm.nih.gov/35602873/) |

---

## 12. Gastrointestinal factors and the differential diagnosis

| # | Reference | PMID |
|---|------|------|
| 66 | Ng QX, et al. **A meta-analysis of the association between Helicobacter pylori infection and hyperemesis gravidarum.** *Helicobacter* 2018. | [29178407](https://pubmed.ncbi.nlm.nih.gov/29178407/) |
| 67 | Sorensen CJ, et al. **Cannabinoid Hyperemesis Syndrome: Diagnosis, Pathophysiology, and Treatment — a Systematic Review.** *J Med Toxicol* 2017. | [28000146](https://pubmed.ncbi.nlm.nih.gov/28000146/) |

---

## 13. Obstetric, fetal, and long-term paediatric outcomes

| # | Reference | PMID |
|---|------|------|
| 68 | Vandraas KF, et al. **Hyperemesis gravidarum and birth outcomes — a population-based cohort study of 2.2 million births in the Norwegian Birth Registry.** *BJOG* 2013. | [24021026](https://pubmed.ncbi.nlm.nih.gov/24021026/) |
| 69 | Nijsten K, et al. **Long-term health outcomes of children born to mothers with hyperemesis gravidarum: a systematic review and meta-analysis.** *Am J Obstet Gynecol* 2022. | [35367190](https://pubmed.ncbi.nlm.nih.gov/35367190/) |
| 70 | Pont S, et al. **Long-term health, neurodevelopmental, and educational outcomes of children born to mothers with hyperemesis gravidarum: a population-based sibling-design study.** *Am J Obstet Gynecol* 2025. | [40064411](https://pubmed.ncbi.nlm.nih.gov/40064411/) |
| 71 | Auger N, et al. **Hyperemesis gravidarum and the risk of offspring morbidity: a longitudinal cohort study.** *Eur J Pediatr* 2024. | [38884821](https://pubmed.ncbi.nlm.nih.gov/38884821/) |
| 72 | Koren G, et al. **The protective effects of nausea and vomiting of pregnancy against adverse fetal outcome — a systematic review.** *Reprod Toxicol* 2014. | [24893173](https://pubmed.ncbi.nlm.nih.gov/24893173/) |
| 73 | Boskovic R, et al. **Is lack of morning sickness teratogenic? A prospective controlled study.** *Birth Defects Res A* 2004. | [15329830](https://pubmed.ncbi.nlm.nih.gov/15329830/) |

**A tension with the model.** References 72 and 73 are the observation that pregnancies **with** NVP have
better outcomes. That means one cannot assume a treatment that abolishes the GDF15 signal is itself harmless,
and the model does not resolve this tension — it only states it explicitly in
cluster 16.

---

## 14. Mechanism-based therapeutics — targeting GDF15/GFRAL (all untested in pregnancy)

> The largest prediction of the model (the anti-GDF15 antibody) and its largest unknown are both here.
> The human data for ponsegromab in cancer cachexia show that ligand blockade
> can be safe and effective, but **there is no evidence whatever for use in the first trimester**, and
> references 13 and 14 raise a concern pointing the opposite way.

| # | Reference | PMID |
|---|------|------|
| 74 | Groarke JD, et al. **Ponsegromab for the Treatment of Cancer Cachexia.** *N Engl J Med* 2024. | [39282907](https://pubmed.ncbi.nlm.nih.gov/39282907/) |
| 75 | Groarke JD, et al. **Phase 2 study of the efficacy and safety of ponsegromab in patients with cancer cachexia: PROACC-1 study design.** *J Cachexia Sarcopenia Muscle* 2024. | [38500292](https://pubmed.ncbi.nlm.nih.gov/38500292/) |
| 76 | Crawford J, et al. **A Phase Ib First-In-Patient Study Assessing the Safety, Tolerability, Pharmacokinetics, and Pharmacodynamics of Ponsegromab in Participants with Cancer Cachexia.** *Clin Cancer Res* 2024. | [37982848](https://pubmed.ncbi.nlm.nih.gov/37982848/) |
| 77 | Lee BY, et al. **GDNF family receptor alpha-like antagonist antibody alleviates chemotherapy-induced cachexia in melanoma-bearing mice.** *J Cachexia Sarcopenia Muscle* 2023. | [37017344](https://pubmed.ncbi.nlm.nih.gov/37017344/) |

---

## Summary — which reference determined which parameter

| Parameter | Fitted/fixed | Basis |
|----------|-----------|------|
| `ALPHA`, `TAU_SP` | **fitted** (natural history) | #9, #18 — peak at 9-11 weeks, resolution by 16 weeks, GDF15 rises to term |
| `W_VAG` | **fitted** (ondansetron) | #35 — ΔPUQE −0.51 |
| `E0` | **fitted** (mirtazapine) | #35 — ΔPUQE −1.86 |
| `R_H1`,`R_A2`,`R_A2D`,… | fixed (pharmacology) | #47, #48 — Ki by receptor class and the NTS delivery gate |
| PK/Ki/Kp of each drug | fixed (literature) | #44-#47 |
| `MET_EMAX_GDF` | fixed | #7 — metformin → about a 2-fold GDF15 |
| `TOBACCO_F` | fixed | #15 — smoking gives a smaller rise than metformin |
| `THAL_FOLD` | fixed | #16, #17 |
| `THI0`, `KEL_THI`, `THI_PER_CHO` | fixed (physiology) | #57, #58 |
| `GJ_CL`, `GJ_H`, `KREN_HCO3`, `NA_HCO3_COUP` | fixed (physiology) | #55, #56 |
| `AH_HCG` | fixed | #62, #63 |
| The PUQE band boundaries | fixed (instrument definition) | #20, #21 |

**Predictions (results that were not fitted)**: the protection from β-thalassaemia, the pre-pregnancy protection and
in-pregnancy futility of metformin, the intermediate protection from smoking, the effect size of doxylamine, the superiority of gabapentin,
the effect size of clonidine, the negative corticosteroid result, thyrotoxicosis appearing on the fetal production axis
alone, and the harm of glucose without thiamine.

**A stated failure**: metoclopramide is predicted to be effectively useless, whereas in the trials it is
much like promethazine. This is the easiest place to refute the model.

---

## Tooling references

- mrgsolve: <https://mrgsolve.org/>
- gPKPDviz (a mrgsolve-based Shiny app): <https://pmc.ncbi.nlm.nih.gov/articles/PMC10941578/> · <https://github.com/Genentech/gPKPDviz/>
- Graphviz: <https://graphviz.org/>
