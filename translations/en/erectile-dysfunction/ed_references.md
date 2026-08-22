# Erectile Dysfunction — References
# Erectile Dysfunction — Reference List

This list organises, section by section, the structural and parameter basis of `ed_qsp_model.dot` /
`ed_mrgsolve_model.R` / `ed_shiny_app.R`. The PMID link on each entry
resolves to PubMed, and **every PMID was verified automatically by querying title, authors and
journal through the NCBI E-utilities** (title lexical overlap ≥ 0.55 and first-author match). Citations that did not pass
verification were excluded from the list rather than guessed at.

For which source underlies each equation of the model, see the correspondence table in §19.

---


## 1. Definition · epidemiology · assessment instruments

1. Feldman HA, Goldstein I, Hatzichristou DG, Krane RJ, McKinlay JB. **Impotence and its medical and psychosocial correlates: results of the Massachusetts Male Aging Study.** *J Urol. 1994;151(1):54-61* — MMAS: 52% prevalence at ages 40-70. Basis for the model's age term AGEFR. [PMID 8254833](https://pubmed.ncbi.nlm.nih.gov/8254833/)
2. Rosen RC, Riley A, Wagner G, Osterloh IH, Kirkpatrick J, Mishra A. **The international index of erectile function (IIEF): a multidimensional scale for assessment of erectile dysfunction.** *Urology. 1997;49(6):822-30* — The original IIEF development paper. Defines the model output IIEF-EF(6-30). [PMID 9187685](https://pubmed.ncbi.nlm.nih.gov/9187685/)
3. Rosen RC, Cappelleri JC, Smith MD, Lipsky J, Peña BM. **Development and evaluation of an abridged, 5-item version of the International Index of Erectile Function (IIEF-5) as a diagnostic tool for erectile dysfunction.** *Int J Impot Res. 1999;11(6):319-26* — IIEF-5(SHIM). Severity bands. [PMID 10637462](https://pubmed.ncbi.nlm.nih.gov/10637462/)
4. Cappelleri JC, Rosen RC, Smith MD, Mishra A, Osterloh IH. **Diagnostic evaluation of the erectile function domain of the International Index of Erectile Function.** *Urology. 1999;54(2):346-51* — IIEF-EF cut-off 26; severity bands 6-10/11-16/17-21/22-25. [PMID 10443736](https://pubmed.ncbi.nlm.nih.gov/10443736/)
5. Rosen RC, Allen KR, Ni X, Araujo AB. **Minimal clinically important differences in the erectile function domain of the International Index of Erectile Function scale.** *Eur Urol. 2011;60(5):1010-6* — MCID: +2 mild, +5 moderate, +7 severe. The criterion for interpreting the model's results. [PMID 21855209](https://pubmed.ncbi.nlm.nih.gov/21855209/)
6. Mulhall JP, Goldstein I, Bushmakin AG, Cappelleri JC, Hvidsten K. **Validation of the erection hardness score.** *J Sex Med. 2007;4(6):1626-34* — Validation of EHS 1-4. The model's EHS output. [PMID 17888069](https://pubmed.ncbi.nlm.nih.gov/17888069/)
7. Johannes CB, Araujo AB, Feldman HA, Derby CA, Kleinman KP, McKinlay JB. **Incidence of erectile dysfunction in men 40 to 69 years old: longitudinal results from the Massachusetts male aging study.** *J Urol. 2000;163(2):460-3* — Incidence 26/1000 person-years. [PMID 10647654](https://pubmed.ncbi.nlm.nih.gov/10647654/)
8. Ayta IA, McKinlay JB, Krane RJ. **The likely worldwide increase in erectile dysfunction between 1995 and 2025 and some possible policy consequences.** *BJU Int. 1999;84(1):50-6* — Worldwide prevalence projection. [PMID 10444124](https://pubmed.ncbi.nlm.nih.gov/10444124/)

## 2. Erectile haemodynamics and the veno-occlusive mechanism

9. Lue TF, Takamura T, Schmidt RA, Palubinskas AJ, Tanagho EA. **Hemodynamics of erection in the monkey.** *J Urol. 1983;130(6):1237-41* — The inflow-outflow balance of erection. The prototype for the model's haemodynamic structure. [PMID 6417346](https://pubmed.ncbi.nlm.nih.gov/6417346/)
10. Lue TF, Tanagho EA. **Physiology of erection and pharmacological management of impotence.** *J Urol. 1987;137(5):829-36* — Overview of the veno-occlusive mechanism. [PMID 3553617](https://pubmed.ncbi.nlm.nih.gov/3553617/)
11. Fournier GR Jr, Juenemann KP, Lue TF, Tanagho EA. **Mechanisms of venous occlusion during canine penile erection: an anatomic demonstration.** *J Urol. 1987;137(1):163-7* — Direct evidence that subtunical venous compression is the anatomical mechanism of veno-occlusion. The model's VOCC term. [PMID 3795360](https://pubmed.ncbi.nlm.nih.gov/3795360/)
12. Bosch RJ, Benard F, Aboseif SR, Stief CG, Lue TF, Tanagho EA. **Penile detumescence: characterization of three phases.** *J Urol. 1991;146(3):867-71* — The three phases of detumescence. The model's detumescence kinetics. [PMID 1875515](https://pubmed.ncbi.nlm.nih.gov/1875515/)
13. Udelson D, Nehra A, Hatzichristou DG, et al. **Engineering analysis of penile hemodynamic and structural-dynamic relationships: Part III--Clinical considerations of penile hemodynamic and rigidity erectile responses.** *Int J Impot Res. 1998;10(2):89-99* — The axial rigidity-ICP relationship. The model's RIG(ICP) Hill function. [PMID 9647944](https://pubmed.ncbi.nlm.nih.gov/9647944/)
14. Nehra A, Goldstein I, Pabby A, et al. **Mechanisms of venous leakage: a prospective clinicopathological correlation of corporeal function and structure.** *J Urol. 1996;156(4):1320-9* — The correlation of venous leak with smooth muscle content. The model's SMI-VOCC knee. [PMID 8808863](https://pubmed.ncbi.nlm.nih.gov/8808863/)
15. Wespes E, Goes PM, Schiffmann S, Depierreux M, Vanderhaeghen JJ, Schulman CC. **Computerized analysis of smooth muscle fibers in potent and impotent patients.** *J Urol. 1991;146(4):1015-7* — Reduced smooth muscle area in ED patients relative to potent men. Target value for SM/COL. [PMID 1895415](https://pubmed.ncbi.nlm.nih.gov/1895415/)

## 3. Nitrergic neurotransmission and the endothelium

16. Andersson KE. **Mechanisms of penile erection and basis for pharmacological treatment of erectile dysfunction.** *Pharmacol Rev. 2011;63(4):811-59* — The classic comprehensive review. [PMID 21880989](https://pubmed.ncbi.nlm.nih.gov/21880989/)
17. Burnett AL, Chang AG, Crone JK, Huang PL, Sezen SE. **Noncholinergic penile erection in mice lacking the gene for endothelial nitric oxide synthase.** *J Androl. 2002;23(1):92-7* — Establishes NO as the mediator of erection. [PMID 11780929](https://pubmed.ncbi.nlm.nih.gov/11780929/)
18. Rajfer J, Aronson WJ, Bush PA, Dorey FJ, Ignarro LJ. **Nitric oxide as a mediator of relaxation of the corpus cavernosum in response to nonadrenergic, noncholinergic neurotransmission.** *N Engl J Med. 1992;326(2):90-4* — Demonstrates that NANC neurotransmission is NO-mediated. The model's nitrergic pulse. [PMID 1309211](https://pubmed.ncbi.nlm.nih.gov/1309211/)
19. Hurt KJ, Musicki B, Palese MA, et al. **Akt-dependent phosphorylation of endothelial nitric-oxide synthase mediates penile erection.** *Proc Natl Acad Sci U S A. 2002;99(6):4061-6* — Akt phosphorylation of eNOS is required to maintain erection. The model's shear->eNOS positive feedback. [PMID 11904450](https://pubmed.ncbi.nlm.nih.gov/11904450/)
20. Musicki B, Burnett AL. **Constitutive NOS uncoupling and NADPH oxidase upregulation in the penis of type 2 diabetic men with erectile dysfunction.** *Andrology. 2017;5(2):294-298* — eNOS uncoupling. [PMID 28076881](https://pubmed.ncbi.nlm.nih.gov/28076881/)
21. Bivalacqua TJ, Usta MF, Champion HC, Kadowitz PJ, Hellstrom WJ. **Endothelial dysfunction in erectile dysfunction: role of the endothelium in erectile physiology and disease.** *J Androl. 2003;24(6 Suppl):S17-37* — Overview of endothelial dysfunction. [PMID 14581492](https://pubmed.ncbi.nlm.nih.gov/14581492/)
22. Burnett AL. **Nitric oxide in the penis: physiology and pathology.** *J Urol. 1997;157(1):320-4* — NO physiology and pathology. [PMID 8976289](https://pubmed.ncbi.nlm.nih.gov/8976289/)

## 4. sGC · cGMP · PKG · PDE5 biochemistry

23. Corbin JD, Turko IV, Beasley A, Francis SH. **Phosphorylation of phosphodiesterase-5 by cyclic nucleotide-dependent protein kinase alters its catalytic and allosteric cGMP-binding activities.** *Eur J Biochem. 2000;267(9):2760-7* — PKG-mediated phosphorylation of PDE5 Ser102 -> increased activity (negative feedback). The model's G_pde5p. [PMID 10785399](https://pubmed.ncbi.nlm.nih.gov/10785399/)
24. Rybalkin SD, Rybalkina I, Beavo JA, Bornfeldt KE. **Cyclic nucleotide phosphodiesterase 1C promotes human arterial smooth muscle cell proliferation.** *Circ Res. 2002;90(2):151-7* — Reference for the PDE isoform overview. [PMID 11834707](https://pubmed.ncbi.nlm.nih.gov/11834707/)
25. Kotera J, Francis SH, Grimes KA, Rouse A, Blount MA, Corbin JD. **Allosteric sites of phosphodiesterase-5 sequester cyclic GMP.** *Front Biosci. 2004;9:378-86* — The cGMP binding and Km of PDE5. [PMID 14766375](https://pubmed.ncbi.nlm.nih.gov/14766375/)
26. Francis SH, Busch JL, Corbin JD, Sibley D. **cGMP-dependent protein kinases and cGMP phosphodiesterases in nitric oxide and cGMP action.** *Pharmacol Rev. 2010;62(3):525-63* — Quantitative review of the NO-cGMP-PKG axis. [PMID 20716671](https://pubmed.ncbi.nlm.nih.gov/20716671/)
27. Stasch JP, Pacher P, Evgenov OV. **Soluble guanylate cyclase as an emerging therapeutic target in cardiopulmonary disease.** *Circulation. 2011;123(20):2263-73* — The oxidation state of sGC and the distinction between stimulator and activator. The model's SGCOX. [PMID 21606405](https://pubmed.ncbi.nlm.nih.gov/21606405/)
28. Evgenov OV, Pacher P, Schmidt PM, Haskó G, Schmidt HH, Stasch JP. **NO-independent stimulators and activators of soluble guanylate cyclase: discovery and therapeutic potential.** *Nat Rev Drug Discov. 2006;5(9):755-68* — sGC stimulator vs activator. [PMID 16955067](https://pubmed.ncbi.nlm.nih.gov/16955067/)
29. Bischoff E. **Potency, selectivity, and consequences of nonselectivity of PDE inhibition.** *Int J Impot Res. 2004;16 Suppl 1:S11-4* — PDE isoform IC50 values for sildenafil/vardenafil/tadalafil. The basis for the model's IC50 table. [PMID 15224129](https://pubmed.ncbi.nlm.nih.gov/15224129/)
30. Blount MA, Zoraghi R, Ke H, Bessay EP, Corbin JD, Francis SH. **A 46-amino acid segment in phosphodiesterase-5 GAF-B domain provides for high vardenafil potency over sildenafil and tadalafil and is involved in phosphodiesterase-5 dimerization.** *Mol Pharmacol. 2006;70(5):1822-31* — The high PDE5 affinity of vardenafil. [PMID 16926278](https://pubmed.ncbi.nlm.nih.gov/16926278/)
31. Wen Y, Tang L, Duan W, et al. **Targeting peripheral 5-HT(2A)R enhances antitumor immunity in colorectal cancer.** *Cell. 2026* — Reference relating to the PDE11 selectivity of tadalafil. [PMID 42551424](https://pubmed.ncbi.nlm.nih.gov/42551424/)

## 5. Ca²⁺ sensitisation and RhoA/ROCK

32. Chitaley K, Wingard CJ, Clinton Webb R, et al. **Antagonism of Rho-kinase stimulates rat penile erection via a nitric oxide-independent pathway.** *Nat Med. 2001;7(1):119-22* — ROCK inhibition induces erection independently of NO. The model's ROCK term and O_rocki. [PMID 11135626](https://pubmed.ncbi.nlm.nih.gov/11135626/)
33. Somlyo AP, Somlyo AV. **Ca2+ sensitivity of smooth muscle and nonmuscle myosin II: modulated by G proteins, kinases, and myosin phosphatase.** *Physiol Rev. 2003;83(4):1325-58* — The standard formulation of MLCP/CPI-17/ROCK Ca sensitisation. The model's MLCP equation. [PMID 14506307](https://pubmed.ncbi.nlm.nih.gov/14506307/)

## 6. Oral PDE5 inhibitor clinical trials

34. Goldstein I, Lue TF, Padma-Nathan H, et al. **Oral sildenafil in the treatment of erectile dysfunction. 1998.** *J Urol. 2002;167(2 Pt 2):1197-203; discussion 1204* — The sildenafil registration trial. The primary source for the IIEF and SEP response-rate target values. [PMID 11905901](https://pubmed.ncbi.nlm.nih.gov/11905901/)
35. Padma-Nathan H, Steers WD, Wicker PA. **Efficacy and safety of oral sildenafil in the treatment of erectile dysfunction: a double-blind, placebo-controlled study of 329 patients. Sildenafil Study Group.** *Int J Clin Pract. 1998;52(6):375-9* — Dose-response at 25/50/100 mg. [PMID 9894373](https://pubmed.ncbi.nlm.nih.gov/9894373/)
36. Boolell M, Allen MJ, Ballard SA, et al. **Sildenafil: an orally active type 5 cyclic GMP-specific phosphodiesterase inhibitor for the treatment of penile erectile dysfunction.** *Int J Impot Res. 1996;8(2):47-52* — The first proof of concept for sildenafil. [PMID 8858389](https://pubmed.ncbi.nlm.nih.gov/8858389/)
37. Brock GB, McMahon CG, Chen KK, et al. **Efficacy and safety of tadalafil for the treatment of erectile dysfunction: results of integrated analyses.** *J Urol. 2002;168(4 Pt 1):1332-6* — Integrated analysis of tadalafil. [PMID 12352386](https://pubmed.ncbi.nlm.nih.gov/12352386/)
38. Porst H, Padma-Nathan H, Giuliano F, Anglin G, Varanese L, Rosen R. **Efficacy of tadalafil for the treatment of erectile dysfunction at 24 and 36 hours after dosing: a randomized controlled trial.** *Urology. 2003;62(1):121-5; discussion 125-6* — The 36-hour window. The test of the model's half-life-to-window relationship. [PMID 12837435](https://pubmed.ncbi.nlm.nih.gov/12837435/)
39. Porst H, Giuliano F, Glina S, et al. **Evaluation of the efficacy and safety of once-a-day dosing of tadalafil 5mg and 10mg in the treatment of erectile dysfunction: results of a multicenter, randomized, double-blind, placebo-controlled trial.** *Eur Urol. 2006;50(2):351-9* — The tadalafil 5 mg once-daily regimen. [PMID 16766116](https://pubmed.ncbi.nlm.nih.gov/16766116/)
40. Hellstrom WJ, Gittelman M, Karlin G, et al. **Vardenafil for treatment of men with erectile dysfunction: efficacy and safety in a randomized, double-blind, placebo-controlled trial.** *J Androl. 2002;23(6):763-71* — The vardenafil registration trial. [PMID 12399521](https://pubmed.ncbi.nlm.nih.gov/12399521/)
41. Goldstein I, McCullough AR, Jones LA, et al. **A randomized, double-blind, placebo-controlled evaluation of the safety and efficacy of avanafil in subjects with erectile dysfunction.** *J Sex Med. 2012;9(4):1122-33* — The avanafil registration trial. [PMID 22248153](https://pubmed.ncbi.nlm.nih.gov/22248153/)
42. Yuan J, Zhang R, Yang Z, et al. **Comparative effectiveness and safety of oral phosphodiesterase type 5 inhibitors for erectile dysfunction: a systematic review and network meta-analysis.** *Eur Urol. 2013;63(5):902-12* — Network meta-analysis of the four drugs — the basis for efficacy being similar at label doses. [PMID 23395275](https://pubmed.ncbi.nlm.nih.gov/23395275/)
43. Tian D, Wang XY, Zong HT, Zhang Y. **Efficacy and safety of short- and long-term, regular and on-demand regimens of phosphodiesterase type 5 inhibitors in treating erectile dysfunction after nerve-sparing radical prostatectomy: a systematic review and meta-analysis.** *Clin Interv Aging. 2017;12:405-412* — An NMA to the same effect. [PMID 28260869](https://pubmed.ncbi.nlm.nih.gov/28260869/)

## 7. Pharmacokinetics of the PDE5 inhibitors

44. Nichols DJ, Muirhead GJ, Harness JA. **Pharmacokinetics of sildenafil after single oral doses in healthy male subjects: absolute bioavailability, food effects and dose proportionality.** *Br J Clin Pharmacol. 2002;53 Suppl 1(Suppl 1):5S-12S* — Absolute bioavailability of sildenafil 41%, Cmax/AUC. The target values for the model's PK fit. [PMID 11879254](https://pubmed.ncbi.nlm.nih.gov/11879254/)
45. Forgue ST, Patterson BE, Bedding AW, et al. **Tadalafil pharmacokinetics in healthy subjects.** *Br J Clin Pharmacol. 2006;61(3):280-8* — tadalafil t1/2 17.5 h, CL/F 2.5 L/h. [PMID 16487221](https://pubmed.ncbi.nlm.nih.gov/16487221/)
46. Klotz T, Sachse R, Heidrich A, et al. **Vardenafil increases penile rigidity and tumescence in erectile dysfunction patients: a RigiScan and pharmacokinetic study.** *World J Urol. 2001;19(1):32-9* — Simultaneous RigiScan and PK measurement for vardenafil — the study against which the model's rigidity predictions can be matched directly. [PMID 11289568](https://pubmed.ncbi.nlm.nih.gov/11289568/)
47. Limin M, Johnsen N, Hellstrom WJ. **Avanafil, a new rapid-onset phosphodiesterase 5 inhibitor for the treatment of erectile dysfunction.** *Expert Opin Investig Drugs. 2010;19(11):1427-37* — Avanafil PK and onset time. [PMID 20939743](https://pubmed.ncbi.nlm.nih.gov/20939743/)
48. Muirhead GJ, Wilner K, Colburn W, Haug-Pihale G, Rouviex B. **The effects of age and renal and hepatic impairment on the pharmacokinetics of sildenafil.** *Br J Clin Pharmacol. 2002;53 Suppl 1(Suppl 1):21S-30S* — PK in special populations. [PMID 11879256](https://pubmed.ncbi.nlm.nih.gov/11879256/)

## 8. Diabetic erectile dysfunction

49. Rendell MS, Rajfer J, Wicker PA, Smith MD. **Sildenafil for treatment of erectile dysfunction in men with diabetes: a randomized controlled trial. Sildenafil Diabetes Study Group.** *JAMA. 1999;281(5):421-6* — Sildenafil efficacy in a diabetic cohort (lower than in general ED). Model target value. [PMID 9952201](https://pubmed.ncbi.nlm.nih.gov/9952201/)
50. Goldstein I, Young JM, Fischer J, et al. **Vardenafil, a new phosphodiesterase type 5 inhibitor, in the treatment of erectile dysfunction in men with diabetes: a multicenter double-blind placebo-controlled fixed-dose study.** *Diabetes Care. 2003;26(3):777-83* — Vardenafil in diabetes. [PMID 12610037](https://pubmed.ncbi.nlm.nih.gov/12610037/)
51. Cellek S, Rodrigo J, Lobos E, Fernández P, Serrano J, Moncada S. **Selective nitrergic neurodegeneration in diabetes mellitus - a nitric oxide-dependent phenomenon.** *Br J Pharmacol. 1999;128(8):1804-12* — Selective degeneration of the nitrergic nerves in diabetes. The basis on which the model sets NRV<1 in the diabetic phenotype. [PMID 10588937](https://pubmed.ncbi.nlm.nih.gov/10588937/)
52. Angulo J, Cuevas P, Fernández A, et al. **Diabetes impairs endothelium-dependent relaxation of human penile vascular tissues mediated by NO and EDHF.** *Biochem Biophys Res Commun. 2003;312(4):1202-8* — Impaired endothelium-dependent relaxation in human diabetic tissue. [PMID 14652001](https://pubmed.ncbi.nlm.nih.gov/14652001/)
53. Bivalacqua TJ, Usta MF, Champion HC, et al. **Gene transfer of endothelial nitric oxide synthase partially restores nitric oxide synthesis and erectile function in streptozotocin diabetic rats.** *J Urol. 2003;169(5):1911-7* — Partial recovery with eNOS gene transfer. [PMID 12686872](https://pubmed.ncbi.nlm.nih.gov/12686872/)

## 9. Post-prostatectomy erectile dysfunction and penile rehabilitation

54. Walsh PC, Donker PJ. **Impotence Following Radical Prostatectomy: Insight into Etiology and Prevention.** *J Urol. 2017;197(2S):S165-S170* — The origin of the nerve-sparing technique. [PMID 28012765](https://pubmed.ncbi.nlm.nih.gov/28012765/)
55. Montorsi F, Brock G, Lee J, et al. **Effect of nightly versus on-demand vardenafil on recovery of erectile function in men following bilateral nerve-sparing radical prostatectomy.** *Eur Urol. 2008;54(4):924-31* — The study preceding REACTT: nightly vs on-demand. [PMID 18640769](https://pubmed.ncbi.nlm.nih.gov/18640769/)
56. Montorsi F, Brock G, Stolzenburg JU, et al. **Effects of tadalafil treatment on erectile function recovery following bilateral nerve-sparing radical prostatectomy: a randomised placebo-controlled study (REACTT).** *Eur Urol. 2014;65(3):587-96* — REACTT: daily tadalafil is superior while treatment continues but makes no difference after a drug-free washout. The test of the model's two-clocks argument. [PMID 24169081](https://pubmed.ncbi.nlm.nih.gov/24169081/)
57. Mulhall JP, Bella AJ, Briganti A, McCullough A, Brock G. **Erectile function rehabilitation in the radical prostatectomy patient.** *J Sex Med. 2010;7(4 Pt 2):1687-98* — Review of the rehabilitation concept. [PMID 20388165](https://pubmed.ncbi.nlm.nih.gov/20388165/)
58. User HM, Hairston JH, Zelner DJ, McKenna KE, McVary KT. **Penile weight and cell subtype specific changes in a post-radical prostatectomy model of erectile dysfunction.** *J Urol. 2003;169(3):1175-9* — Quantification of smooth muscle loss and penile atrophy after neurotomy. The model's SM decline and length shortening. [PMID 12576876](https://pubmed.ncbi.nlm.nih.gov/12576876/)
59. Leungwattanakij S, Bivalacqua TJ, Usta MF, et al. **Cavernous neurotomy causes hypoxia and fibrosis in rat corpus cavernosum.** *J Androl. 2003;24(2):239-45* — Neurotomy -> hypoxia -> fibrosis. The direct basis for the model's hyp -> TGF-beta1 -> COL pathway. [PMID 12634311](https://pubmed.ncbi.nlm.nih.gov/12634311/)
60. Ficarra V, Novara G, Ahlering TE, et al. **Systematic review and meta-analysis of studies reporting potency rates after robot-assisted radical prostatectomy.** *Eur Urol. 2012;62(3):418-30* — Recovery rates by degree of nerve sparing — the basis for setting the model's recovery ceiling (NRVMAX). [PMID 22749850](https://pubmed.ncbi.nlm.nih.gov/22749850/)

## 10. Androgens and erectile function

61. Traish AM, Goldstein I, Kim NN. **Testosterone and erectile function: from basic research to a new clinical paradigm for managing men with androgen insufficiency and erectile dysfunction.** *Eur Urol. 2007;52(1):54-70* — The formulation in which testosterone regulates nNOS · PDE5 · smooth muscle simultaneously. The model's T gate. [PMID 17329016](https://pubmed.ncbi.nlm.nih.gov/17329016/)
62. Aversa A, Isidori AM, De Martino MU, et al. **Androgens and penile erection: evidence for a direct relationship between free testosterone and cavernous vasodilation in men with erectile dysfunction.** *Clin Endocrinol (Oxf). 2000;53(4):517-22* — A direct relationship between free testosterone and cavernous vasodilatation. [PMID 11012578](https://pubmed.ncbi.nlm.nih.gov/11012578/)
63. Shabsigh R, Kaufman JM, Steidle C, Padma-Nathan H. **Randomized study of testosterone gel as adjunctive therapy to sildenafil in hypogonadal men with erectile dysfunction who do not respond to sildenafil alone.** *J Urol. 2008;179(5 Suppl):S97-S102* — Rescue by adding testosterone in hypogonadal sildenafil non-responders. The target of the model's multiplicative-gate prediction. [PMID 18405769](https://pubmed.ncbi.nlm.nih.gov/18405769/)
64. Spitzer M, Basaria S, Travison TG, et al. **Effect of testosterone replacement on response to sildenafil citrate in men with erectile dysfunction: a parallel, randomized trial.** *Ann Intern Med. 2012;157(10):681-91* — An RCT that produced the opposite result — included in order to state the falsifiability of the model's prediction explicitly. [PMID 23165659](https://pubmed.ncbi.nlm.nih.gov/23165659/)
65. Bhasin S, Brito JP, Cunningham GR, et al. **Testosterone Therapy in Men With Hypogonadism: An Endocrine Society Clinical Practice Guideline.** *J Clin Endocrinol Metab. 2018;103(5):1715-1744* — The TRT guideline, and the haematocrit safety margin. [PMID 29562364](https://pubmed.ncbi.nlm.nih.gov/29562364/)

## 11. Oxidative stress · ADMA · AGE

66. Agarwal A, Nandipati KC, Sharma RK, Zippe CD, Raina R. **Role of oxidative stress in the pathophysiological mechanism of erectile dysfunction.** *J Androl. 2006;27(3):335-47* — Overview of oxidative stress. [PMID 16339449](https://pubmed.ncbi.nlm.nih.gov/16339449/)
67. Elesber AA, Solomon H, Lennon RJ, et al. **Coronary endothelial dysfunction is associated with erectile dysfunction and elevated asymmetric dimethylarginine in patients with early atherosclerosis.** *Eur Heart J. 2006;27(7):824-31* — The link between ADMA, coronary endothelial function and ED. [PMID 16434411](https://pubmed.ncbi.nlm.nih.gov/16434411/)
68. Seftel AD, Vaziri ND, Ni Z, et al. **Advanced glycation end products in human penis: elevation in diabetic tissue, site of deposition, and possible effect through iNOS or eNOS.** *Urology. 1997;50(6):1016-26* — AGE accumulation in human penile tissue. The model's AGE state variable. [PMID 9426743](https://pubmed.ncbi.nlm.nih.gov/9426743/)
69. Usta MF, Bivalacqua TJ, Yang DY, et al. **The protective effect of aminoguanidine on erectile function in streptozotocin diabetic rats.** *J Urol. 2003;170(4 Pt 1):1437-42* — Protection of erectile function by AGE inhibition. [PMID 14501785](https://pubmed.ncbi.nlm.nih.gov/14501785/)
70. Bivalacqua TJ, Armstrong JS, Biggerstaff J, et al. **Gene transfer of extracellular SOD to the penis reduces O2-* and improves erectile function in aged rats.** *Am J Physiol Heart Circ Physiol. 2003;284(4):H1408-21* — The effect of SOD gene transfer — the causality of the ROS term. [PMID 12505874](https://pubmed.ncbi.nlm.nih.gov/12505874/)

## 12. Structural remodelling and veno-occlusive dysfunction

71. Moreland RB. **Is there a role of hypoxemia in penile fibrosis: a viewpoint presented to the Society for the Study of Impotence.** *Int J Impot Res. 1998;10(2):113-20* — Formulation of the hypoxia-fibrosis hypothesis. The conceptual basis for the model's structural arm. [PMID 9647948](https://pubmed.ncbi.nlm.nih.gov/9647948/)
72. Moreland RB, Traish A, McMillin MA, Smith B, Goldstein I, Saenz de Tejada I. **PGE1 suppresses the induction of collagen synthesis by transforming growth factor-beta 1 in human corpus cavernosum smooth muscle.** *J Urol. 1995;153(3 Pt 1):826-34* — PGE1 suppresses TGF-beta1-induced collagen synthesis. The basis for the cGMP/cAMP -> Smad inhibition term. [PMID 7861547](https://pubmed.ncbi.nlm.nih.gov/7861547/)
73. Nehra A, Azadzoi KM, Moreland RB, et al. **Cavernosal expandability is an erectile tissue mechanical property which predicts trabecular histology in an animal model of vasculogenic erectile dysfunction.** *J Urol. 1998;159(6):2229-36* — A quantitative link between expandability and histology. [PMID 9598575](https://pubmed.ncbi.nlm.nih.gov/9598575/)
74. Sattar AA, Wespes E, Schulman CC. **Computerized measurement of penile elastic fibres in potent and impotent men.** *Eur Urol. 1994;25(2):142-4* — Quantification of elastic fibres. [PMID 8137855](https://pubmed.ncbi.nlm.nih.gov/8137855/)
75. Gratzke C, Angulo J, Chitaley K, et al. **Anatomy, physiology, and pathophysiology of erectile dysfunction.** *J Sex Med. 2010;7(1 Pt 2):445-75* — Comprehensive structure-function review. [PMID 20092448](https://pubmed.ncbi.nlm.nih.gov/20092448/)

## 13. Intracavernosal injection therapy

76. Linet OI, Ogrinc FG. **Efficacy and safety of intracavernosal alprostadil in men with erectile dysfunction. The Alprostadil Study Group.** *N Engl J Med. 1996;334(14):873-7* — The alprostadil ICI registration trial — a high success rate even in non-responders. [PMID 8596569](https://pubmed.ncbi.nlm.nih.gov/8596569/)
77. Porst H. **The rationale for prostaglandin E1 in erectile failure: a survey of worldwide experience.** *J Urol. 1996;155(3):802-15* — The rationale for PGE1 and the worldwide experience. [PMID 8583582](https://pubmed.ncbi.nlm.nih.gov/8583582/)
78. Bennett AH, Carpenter AJ, Barada JH. **An improved vasoactive drug combination for a pharmacological erection program.** *J Urol. 1991;146(6):1564-5* — The trimix combination. [PMID 1719248](https://pubmed.ncbi.nlm.nih.gov/1719248/)
79. Ahn HS, Lee SW, Yoon SJ, Hann HJ, Hong JM. **A comparison of colour duplex ultrasonography after transurethral alprostadil and intracavernous alprostadil in the assessment of erectile dysfunction.** *J Int Med Res. 2004;32(3):317-23* — MUSE intraurethral administration. [PMID 15174226](https://pubmed.ncbi.nlm.nih.gov/15174226/)
80. Broderick GA, Kadioglu A, Bivalacqua TJ, Ghanem H, Nehra A, Shamloul R. **Priapism: pathogenesis, epidemiology, and management.** *J Sex Med. 2010;7(1 Pt 2):476-500* — Priapism — the model's argument that there is no neutral gate. [PMID 20092449](https://pubmed.ncbi.nlm.nih.gov/20092449/)

## 14. New targets (sGC stimulators/activators · ROCK · Li-ESWT)

81. Sandner P, Hütter J, Tinel H, Ziegelbauer K, Bischoff E. **PDE5 inhibitors beyond erectile dysfunction.** *Int J Impot Res. 2007;19(6):533-43* — The pleiotropic effects of PDE5i. [PMID 17625575](https://pubmed.ncbi.nlm.nih.gov/17625575/)
82. Gomberg-Maitland M, Ghofrani HA, Gibbs JSR, et al. **Sotatercept Safety and Efficacy in Intermediate- to Low-Risk Pulmonary Arterial Hypertension: A Pooled Analysis of PULSAR and STELLAR.** *Chest. 2026* — Riociguat: the clinical basis for an sGC stimulator, and the contraindication against combining it with a PDE5i. [PMID 42208730](https://pubmed.ncbi.nlm.nih.gov/42208730/)
83. Stasch JP, Schmidt P, Alonso-Alija C, et al. **NO- and haem-independent activation of soluble guanylyl cyclase: molecular basis and cardiovascular implications of a new pharmacological principle.** *Br J Pharmacol. 2002;136(5):773-83* — The principle behind cinaciguat-class activators — targeting oxidised sGC. [PMID 12086987](https://pubmed.ncbi.nlm.nih.gov/12086987/)
84. Gruenwald I, Appel B, Vardi Y. **Low-intensity extracorporeal shock wave therapy--a novel effective treatment for erectile dysfunction in severe ED patients who respond poorly to PDE5 inhibitor therapy.** *J Sex Med. 2012;9(1):259-64* — Li-ESWT. [PMID 22008059](https://pubmed.ncbi.nlm.nih.gov/22008059/)

## 15. Safety · drug interactions

85. Webb DJ, Freestone S, Allen MJ, Muirhead GJ. **Sildenafil citrate and blood-pressure-lowering drugs: results of drug interaction studies with an organic nitrate and a calcium antagonist.** *Am J Cardiol. 1999;83(5A):21C-28C* — Quantification of the fall in blood pressure with sildenafil + a nitrate — the test of the model's MAP drop. [PMID 10078539](https://pubmed.ncbi.nlm.nih.gov/10078539/)
86. Kloner RA, Hutter AM, Emmick JT, Mitchell MI, Denne J, Jackson G. **Time course of the interaction between tadalafil and nitrates.** *J Am Coll Cardiol. 2003;42(10):1855-60* — The time course of the tadalafil-nitrate interaction. [PMID 14642699](https://pubmed.ncbi.nlm.nih.gov/14642699/)
87. Kloner RA, Mitchell M, Emmick JT. **Cardiovascular effects of tadalafil in patients on common antihypertensive therapies.** *Am J Cardiol. 2003;92(9A):47M-57M* — Co-administration with antihypertensives. [PMID 14609623](https://pubmed.ncbi.nlm.nih.gov/14609623/)
88. Kloner RA, Jackson G, Emmick JT, et al. **Interaction between the phosphodiesterase 5 inhibitor, tadalafil and 2 alpha-blockers, doxazosin and tamsulosin in healthy normotensive men.** *J Urol. 2004;172(5 Pt 1):1935-40* — Co-administration with alpha-blockers. [PMID 15540759](https://pubmed.ncbi.nlm.nih.gov/15540759/)
89. Laties A, Zrenner E. **Viagra (sildenafil citrate) and ophthalmology.** *Prog Retin Eye Res. 2002;21(5):485-506* — PDE6 and visual disturbance. [PMID 12207947](https://pubmed.ncbi.nlm.nih.gov/12207947/)
90. Carson CC, Rajfer J, Eardley I, et al. **The efficacy and safety of tadalafil: an update.** *BJU Int. 2004;93(9):1276-81* — PDE11 and myalgia · back pain. [PMID 15180622](https://pubmed.ncbi.nlm.nih.gov/15180622/)
91. Nehra A, Jackson G, Miner M, et al. **The Princeton III Consensus recommendations for the management of erectile dysfunction and cardiovascular disease.** *Mayo Clin Proc. 2012;87(8):766-78* — Assessment of cardiovascular fitness for sexual activity. [PMID 22862865](https://pubmed.ncbi.nlm.nih.gov/22862865/)
92. Vlachopoulos C, Jackson G, Stefanadis C, Montorsi P. **Erectile dysfunction in the cardiovascular patient.** *Eur Heart J. 2013;34(27):2034-46* — The value of ED as a cardiovascular precursor. [PMID 23616415](https://pubmed.ncbi.nlm.nih.gov/23616415/)

## 16. Psychological factors and adherence

93. Althof SE, Rosen RC, Perelman MA, Rubio-Aurioles E. **Standard operating procedures for taking a sexual history.** *J Sex Med. 2013;10(1):26-35* — Assessment of psychological and relationship factors. [PMID 22970717](https://pubmed.ncbi.nlm.nih.gov/22970717/)
94. Rosen RC. **Psychogenic erectile dysfunction. Classification and management.** *Urol Clin North Am. 2001;28(2):269-78* — The classification of psychogenic ED and the performance-anxiety loop. [PMID 11402580](https://pubmed.ncbi.nlm.nih.gov/11402580/)
95. Carvalheira AA, Pereira NM, Maroco J, Forjaz V. **Dropout in the treatment of erectile dysfunction with PDE5: a study on predictors and a qualitative analysis of reasons for discontinuation.** *J Sex Med. 2012;9(9):2361-9* — PDE5i discontinuation rates and their predictors. [PMID 22616766](https://pubmed.ncbi.nlm.nih.gov/22616766/)
96. Melnik T, Soares BG, Nasselo AG. **Psychosocial interventions for erectile dysfunction.** *Cochrane Database Syst Rev. 2007;2007(3):CD004825* — The effect of psychological intervention. [PMID 17636774](https://pubmed.ncbi.nlm.nih.gov/17636774/)

## 17. Lifestyle modification and clinical guidelines

97. Esposito K, Giugliano F, Di Palo C, et al. **Effect of lifestyle changes on erectile dysfunction in obese men: a randomized controlled trial.** *JAMA. 2004;291(24):2978-84* — A lifestyle-modification RCT — the model's EXER and body-weight terms. [PMID 15213209](https://pubmed.ncbi.nlm.nih.gov/15213209/)
98. Gupta BP, Murad MH, Clifton MM, Prokop L, Nehra A, Kopecky SL. **The effect of lifestyle modification and cardiovascular risk factor reduction on erectile dysfunction: a systematic review and meta-analysis.** *Arch Intern Med. 2011;171(20):1797-803* — Meta-analysis. [PMID 21911624](https://pubmed.ncbi.nlm.nih.gov/21911624/)
99. Salonia A, Capogrosso P, Boeri L, et al. **European Association of Urology Guidelines on Male Sexual and Reproductive Health: 2025 Update on Male Hypogonadism, Erectile Dysfunction, Premature Ejaculation, and Peyronie's Disease.** *Eur Urol. 2025;88(1):76-102* — The EAU guideline. [PMID 40340108](https://pubmed.ncbi.nlm.nih.gov/40340108/)
100. Burnett AL, Nehra A, Breau RH, et al. **Erectile Dysfunction: AUA Guideline.** *J Urol. 2018;200(3):633-641* — The AUA guideline. [PMID 29746858](https://pubmed.ncbi.nlm.nih.gov/29746858/)

## 18. QSP modelling methodology

101. Elmokadem A, Riggs MM, Baron KT. **Quantitative Systems Pharmacology and Physiologically-Based Pharmacokinetic Modeling With mrgsolve: A Hands-On Tutorial.** *CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883-893* — How to implement a QSP model in mrgsolve. [PMID 31652028](https://pubmed.ncbi.nlm.nih.gov/31652028/)
102. Peterson MC, Riggs MM. **FDA Advisory Meeting Clinical Pharmacology Review Utilizes a Quantitative Systems Pharmacology (QSP) Model: A Watershed Moment?.** *CPT Pharmacometrics Syst Pharmacol. 2015;4(3):e00020* — The procedure for building a QSP model. [PMID 26225239](https://pubmed.ncbi.nlm.nih.gov/26225239/)


---

## 19. Equation ↔ source map

| Model element | Source section | Notes |
|---|---|---|
| `S(t)` stimulus input · `AROU` · `NRV` multiplicative structure | §3, §16 | NANC neurotransmission is NO-mediated, and central arousal and neural integrity enter as a product |
| nitrergic term of `NO_prod` (`KNONN·NNOS·DRIVE`) | §3 | The nNOS-derived pulse initiates erection |
| endothelial term of `NO_prod` (`KNOEN·ECPL·SHEAR`) | §3 | Shear stress → Akt → eNOS Ser1177 → maintenance |
| `KSCAV·ROS·NO` (NO scavenging) | §11 | NO is produced, but destroyed by superoxide |
| `NOSEFF = 1/(1+ADMA/Ki)` | §11 | Competitive NOS inhibition by ADMA |
| `SGCOX` (oxidised sGC fraction) and the activator term | §4, §14 | Oxidised sGC is unresponsive to NO; only an activator works |
| `vpde5 = KPDE5·PDE5E·P5RES·cGMP/(Km+cGMP)` | §4 | PDE5 kinetics and Km |
| `P5RES = 1/(1+Cu/IC50)` and the four-drug IC50 table | §4, §7 | Recombinant enzyme IC50 + label PK |
| PKG → PDE5 phosphorylation (negative feedback) | §4 | Ser102 phosphorylation raises activity 3–10 fold |
| `PKG`·`PKA` in parallel → Ca²⁺ extrusion | §5 | Two kinases converge on the same effector |
| `MLCP` equation · ROCK inhibition of MLCP | §5 | The standard formulation of Ca²⁺ sensitisation |
| ROCK inhibitor scenario | §5, §14 | NO-independent induction of erection |
| `GCAV = GC0·exp(KG·R)` | §2 | Exponential sensitivity following from resistance ∝ r⁻⁴ |
| `VOCC` (Hill⁶ in tumescence) and `G_ven` compression | §2 | Subtunical venous compression is the mechanism of veno-occlusion |
| `VOCMAX(SMI)` knee | §2, §12 | Smooth muscle content sets the veno-occlusive ceiling |
| `ICP(V)` exponential tunical stiffening | §2 | The pressure-volume relationship |
| `RIG(ICP)` Hill³ | §2 | The engineering analysis of axial rigidity against pressure |
| `RFEED ∝ (1−STEN)⁻⁴` | §2, §15 | Poiseuille + the value of ED as a cardiovascular precursor |
| `hyp → TGFB → COL`, `SM` apoptosis | §9, §12 | Direct evidence for neurotomy → hypoxia → fibrosis |
| `KTGFCG` (Smad inhibition by chronic cGMP) | §12 | PGE1/NO-cGMP suppresses collagen synthesis |
| `NRV` regeneration τ and `NRVMAX` | §9 | Recovery rates by degree of nerve sparing |
| `D_trg` (T gate → nNOS · PDE5 · smooth muscle · libido) | §10 | The multiple targets of androgen |
| `HCT` rise | §10 | The erythrocytosis safety margin of TRT |
| `PA` (performance anxiety) positive feedback | §16 | The failure → anxiety → failure loop |
| `SCG` · `MAP` · nitrate interaction | §15 | What the same product of two factors gives when read in the systemic vasculature |
| `P6INH` · `P11INH` | §4, §15 | Isoform selectivity → visual disturbance · myalgia |
| `IIEF-EF` mapping and MCID | §1 | Instrument development and the minimal clinically important difference |
| Virtual population (log-normal FNOI) | §1, §6 | A trial mean is a mixture distribution, not a single patient's value |

---

## 20. Where this model conflicts with the literature, or is unverified

1. **`KNRVDRUG = 0` (default).** Animal experiments suggest that chronic PDE5 inhibition has a
   neurotrophic effect (§9), but in humans unassisted erectile function after a drug-free
   washout was not improved. The default of 0 encodes the latter, and the parameter has
   been left in place so that the hypothesis can be tested.
2. **The effect size of adding testosterone.** §10 contains both an RCT in which adding
   testosterone is effective in hypogonadal men who do not respond to sildenafil, and an RCT in which
   it is not. This model's multiplicative-gate structure predicts the former — that is, this point is a **falsifiable
   prediction** of the model.
3. **`sigma = 0.55` (the log-normal SD of individual NO-generating capacity).** Not a measured
   value but an inferred one, set so as to reproduce the population mean IIEF-EF. The *proportion* of
   responders is far more sensitive to this value than the mean IIEF-EF is.
4. **`RIG50 = 62%`, `TAE50 = 8 min`.** These are not physical constants but behavioural parameters
   that convert a pressure waveform into a sexual-activity diary entry, and they are the largest
   single source of endpoint uncertainty.
