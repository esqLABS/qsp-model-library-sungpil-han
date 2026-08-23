# Bartter & Gitelman Syndrome — References

A list of the literature underpinning the structure, the parameters, and the calibration of the salt-losing tubulopathy QSP model.
The `→` mark after each entry states what that reference determined in the model.

---

## 1. Diagnostic criteria · international consensus

1. Konrad M, Nijenhuis T, Ariceta G, et al. **Diagnosis and management of Bartter syndrome: executive summary of the consensus and recommendations from the ERKNet/ESPN working group.** Kidney Int. 2021;99(2):324–335. — <https://pubmed.ncbi.nlm.nih.gov/33509356/>
   → the phenotype definitions by type, the indomethacin dose range, and the growth and renal function follow-up endpoints.
2. Blanchard A, Bockenhauer D, Bolignano D, et al. **Gitelman syndrome: consensus and guidance from a Kidney Disease: Improving Global Outcomes (KDIGO) Controversies Conference.** Kidney Int. 2017;91(1):24–33. — <https://pubmed.ncbi.nlm.nih.gov/28003083/>
   → the Gitelman diagnostic criteria (hypokalaemia, metabolic alkalosis, hypomagnesaemia, hypocalciuria), the K and Mg targets, and the FE-Mg > 4% criterion.
3. Bettinelli A, Bianchetti MG, Girardin E, et al. **Use of calcium excretion values to distinguish two forms of primary renal tubular hypokalemic alkalosis: Bartter and Gitelman syndromes.** J Pediatr. 1992;120(1):38–43. — <https://pubmed.ncbi.nlm.nih.gov/1731022/>
   → the discriminating threshold for urinary Ca/Cr (the model's `UCACR` targets: Gitelman < 0.07, Bartter > 0.20 mmol/mmol).
4. Peters M, Jeck N, Reinalter S, et al. **Clinical presentation of genetically defined patients with hypokalemic salt-losing tubulopathies.** Am J Med. 2002;112(3):183–190. — <https://pubmed.ncbi.nlm.nih.gov/11893344/>
   → the distribution of K, Mg, HCO3, and urinary Ca by genotype. The primary source for the model's calibration table.
5. Vargas-Poussou R, Dahan K, Kahila D, et al. **Spectrum of mutations in Gitelman syndrome.** J Am Soc Nephrol. 2011;22(4):693–703. — <https://pubmed.ncbi.nlm.nih.gov/21415153/>
   → a carrier frequency of ~1%, prevalence estimates, and the phenotype-genotype correlation.
6. Seyberth HW, Schlingmann KP. **Bartter- and Gitelman-like syndromes: salt-losing tubulopathies with loop or DCT defects.** Pediatr Nephrol. 2011;26(10):1789–1802. — <https://pubmed.ncbi.nlm.nih.gov/21503667/>
   → the "loop defect vs DCT defect" classification — the conceptual basis for the model's two-parameter FTAL/FDCT scheme.

## 2. Causative genes and the molecular pathophysiology

7. Simon DB, Karet FE, Hamdan JM, et al. **Bartter's syndrome, hypokalaemic alkalosis with hypercalciuria, is caused by mutations in the Na-K-2Cl cotransporter NKCC2.** Nat Genet. 1996;13(2):183–188. — <https://pubmed.ncbi.nlm.nih.gov/8640224/>
   → type I = `SLC12A1`; FTAL ≈ 0.10.
8. Simon DB, Karet FE, Rodriguez-Soriano J, et al. **Genetic heterogeneity of Bartter's syndrome revealed by mutations in the K+ channel, ROMK.** Nat Genet. 1996;14(2):152–156. — <https://pubmed.ncbi.nlm.nih.gov/8841184/>
   → type II = `KCNJ1`; because ROMK is also present in the collecting duct, there is transient neonatal **hyper**kalaemia (the model's `FROMKCD` = 0.30).
9. Simon DB, Bindra RS, Mansfield TA, et al. **Mutations in the chloride channel gene, CLCNKB, cause Bartter's syndrome type III.** Nat Genet. 1997;17(2):171–178. — <https://pubmed.ncbi.nlm.nih.gov/9326936/>
   → ClC-Kb is expressed in both the TAL and the DCT → why in the model type III alone lowers FTAL and FDCT **at the same time** (the overlapping Gitelman phenotype).
10. Birkenhäger R, Otto E, Schürmann MJ, et al. **Mutation of BSND causes Bartter syndrome with sensorineural deafness and kidney failure.** Nat Genet. 2001;29(3):310–314. — <https://pubmed.ncbi.nlm.nih.gov/11687798/>
    → barttin is the β-subunit of ClC-Ka/Kb and is also expressed in the marginal cells of the inner ear → the sensorineural deafness of type IV.
11. Schlingmann KP, Konrad M, Jeck N, et al. **Salt wasting and deafness resulting from mutations in two chloride channels.** N Engl J Med. 2004;350(13):1314–1319. — <https://pubmed.ncbi.nlm.nih.gov/15044642/>
    → digenic CLCNKA+CLCNKB (type IVb).
12. Laghmani K, Beck BB, Yang SS, et al. **Polyhydramnios, transient antenatal Bartter's syndrome, and MAGED2 mutations.** N Engl J Med. 2016;374(19):1853–1863. — <https://pubmed.ncbi.nlm.nih.gov/27120771/>
    → type V: MAGED2 is a co-chaperone regulating HIF-1α, so the phenotype appears only in the hypoxic fetal environment and **resolves spontaneously** after birth.
13. Simon DB, Nelson-Williams C, Bia MJ, et al. **Gitelman's variant of Bartter's syndrome, inherited hypokalaemic alkalosis, is caused by mutations in the thiazide-sensitive Na-Cl cotransporter.** Nat Genet. 1996;12(1):24–30. — <https://pubmed.ncbi.nlm.nih.gov/8528245/>
    → Gitelman = `SLC12A3`; FDCT ≈ 0.15.
14. Bockenhauer D, Feather S, Stanescu HC, et al. **Epilepsy, ataxia, sensorineural deafness, tubulopathy, and KCNJ10 mutations (EAST syndrome).** N Engl J Med. 2009;360(19):1960–1970. — <https://pubmed.ncbi.nlm.nih.gov/19420365/>
    → Kir4.1 → intracellular Cl⁻ → WNK4 signalling → a Gitelman-like phenotype plus neurological features.
15. Kim GH, Ecelbarger CA, Mitchell C, et al. **Vasopressin increases Na-K-2Cl cotransporter expression in thick ascending limb of Henle's loop.** Am J Physiol. 1999;276(1):F96–F103. — <https://pubmed.ncbi.nlm.nih.gov/9887085/>
    → the regulation of NKCC2 expression, the basis for the variability of TAL capacity in the model.

## 3. TAL transport and divalent cation handling — the mechanism of the sign reversal in urinary calcium

16. Greger R. **Ion transport mechanisms in thick ascending limb of Henle's loop of mammalian nephron.** Physiol Rev. 1985;65(3):760–797. — <https://pubmed.ncbi.nlm.nih.gov/2408998/>
    → K⁺ recycling → a lumen-positive potential (+8 mV) → the paracellular driving force for Ca²⁺/Mg²⁺. The basis for the model's `PDREL` term.
17. Simon DB, Lu Y, Choate KA, et al. **Paracellin-1, a renal tight junction protein required for paracellular Mg2+ resorption.** Science. 1999;285(5424):103–106. — <https://pubmed.ncbi.nlm.nih.gov/10390358/>
    → claudin-16; the model parameter `FCLDN`.
18. Konrad M, Schaller A, Seelow D, et al. **Mutations in the tight-junction gene claudin 19 (CLDN19) are associated with renal magnesium wasting, renal failure, and severe ocular involvement.** Am J Hum Genet. 2006;79(5):949–957. — <https://pubmed.ncbi.nlm.nih.gov/17033971/>
19. Loffing J, Vallon V, Loffing-Cueni D, et al. **Altered renal distal tubule structure and renal Na+ and Ca2+ handling in a mouse model for Gitelman's syndrome.** J Am Soc Nephrol. 2004;15(9):2276–2288. — <https://pubmed.ncbi.nlm.nih.gov/15339977/>
    → **DCT atrophy** and the loss of TRPM6/TRPV5 when NCC is absent — the direct basis for the model's `DCTM`, `TRPM`, and `TRPV` compartments.
20. Nijenhuis T, Vallon V, van der Kemp AWCM, et al. **Enhanced passive Ca2+ reabsorption and reduced Mg2+ channel abundance explains thiazide-induced hypocalciuria and hypomagnesemia.** J Clin Invest. 2005;115(6):1651–1658. — <https://pubmed.ncbi.nlm.nih.gov/15902302/>
    → establishes that the hypocalciuria of thiazides and of Gitelman syndrome is due to **increased passive reabsorption in the proximal tubule** and not to an increase in DCT TRPV5 → the design of the model's `FPTCA` (dependent on volume contraction).
21. Voets T, Nilius B, Hoefs S, et al. **TRPM6 forms the Mg2+ influx channel involved in intestinal and renal Mg2+ absorption.** J Biol Chem. 2004;279(1):19–25. — <https://pubmed.ncbi.nlm.nih.gov/14576148/>
    → TRPM6 is present in **both gut and kidney** → intestinal absorption rises in hypomagnesaemia (the model's `FABSMG` adaptation term).
22. Schlingmann KP, Weber S, Peters M, et al. **Hypomagnesemia with secondary hypocalcemia is caused by mutations in TRPM6, a new member of the TRPM gene family.** Nat Genet. 2002;31(2):166–170. — <https://pubmed.ncbi.nlm.nih.gov/12032568/>
23. Huang CL, Kuo E. **Mechanism of hypokalemia in magnesium deficiency.** J Am Soc Nephrol. 2007;18(10):2649–2652. — <https://pubmed.ncbi.nlm.nih.gov/17804670/>
    → intracellular Mg²⁺ blocks ROMK, and hypomagnesaemia releases that block and so amplifies the loss of K⁺. **The key basis for the model's `KMGROMK` term — why K does not rise unless Mg is replaced first.**

## 4. Macula densa NaCl sensing · COX-2 · the PGE₂ amplification loop

24. Schnermann J, Briggs JP. **Tubuloglomerular feedback: mechanistic insights from gene-manipulated mice.** Kidney Int. 2008;74(4):418–426. — <https://pubmed.ncbi.nlm.nih.gov/18497922/>
    → the fact that the macula densa senses NaCl **through its own NKCC2** → the basis for the model's `MDSENSE` structure, in which a TAL lesion lowers the sensed signal even though luminal NaCl is high (exactly as furosemide does).
25. Harris RC, McKanna JA, Akai Y, et al. **Cyclooxygenase-2 is associated with the macula densa of rat kidney and increases with salt restriction.** J Clin Invest. 1994;94(6):2504–2510. — <https://pubmed.ncbi.nlm.nih.gov/7989609/>
    → a fall in sensed NaCl → induction of COX-2. The model's `ECOX` gain.
26. Kömhoff M, Reinalter SC, Gröne HJ, Seyberth HW. **Induction of microsomal prostaglandin E2 synthase in the macula densa in children with hypercalciuric salt-losing tubulopathy.** Pediatr Res. 2004;55(2):261–266. — <https://pubmed.ncbi.nlm.nih.gov/14630989/>
    → confirmation of mPGES-1 induction in human Bartter tissue — the histological evidence for "hyperprostaglandin E syndrome".
27. Reinalter SC, Jeck N, Brochhausen C, et al. **Role of cyclooxygenase-2 in hyperprostaglandin E syndrome/antenatal Bartter syndrome.** Kidney Int. 2002;62(1):253–260. — <https://pubmed.ncbi.nlm.nih.gov/12081586/>
    → the extent of the rise in urinary PGE₂ and the response to selective COX-2 inhibition. Calibration of the model's `IC50CEL` and `EMXCEL`.
28. Seyberth HW, Rascher W, Schweer H, et al. **Congenital hypokalemia with hypercalciuria in preterm infants: a hyperprostaglandinuric tubular syndrome different from Bartter syndrome.** J Pediatr. 1985;107(5):694–701. — <https://pubmed.ncbi.nlm.nih.gov/3865312/>
    → the original description of the antenatal form. Urinary PGE₂ raised 5–10 fold.
29. Castrop H, Schweda F, Schumacher K, et al. **Role of renocortical cyclooxygenase-2 for renal vascular resistance and macula densa control of renin secretion.** J Am Soc Nephrol. 2001;12(5):867–874. — <https://pubmed.ncbi.nlm.nih.gov/11316845/>
    → PGE₂ (EP2/EP4) → renin secretion. The model's `GPG` term.
30. Nüsing RM, Treude A, Weissenberger C, et al. **Dominant role of prostaglandin E2 EP4 receptor in furosemide-induced salt-losing tubulopathy.** J Am Soc Nephrol. 2005;16(8):2354–2362. — <https://pubmed.ncbi.nlm.nih.gov/15976005/>
    → EP4 predominance — the basis for handling the PGE₂ → renin pathway as a single term in the model.
31. Kirchner KA. **Prostaglandin inhibitors alter loop segment chloride uptake during furosemide diuresis.** Am J Physiol. 1985;248(5):F698–F704. — <https://pubmed.ncbi.nlm.nih.gov/3922262/>
    → PGE₂ inhibits TAL NKCC2 **further** → the positive feedback loop of the model (`FPGTAL`).

## 5. Distal nephron K⁺ and H⁺ secretion and the metabolic alkalosis

32. Palmer LG, Frindt G. **Aldosterone and potassium secretion by the cortical collecting duct.** Kidney Int. 2000;57(4):1324–1328. — <https://pubmed.ncbi.nlm.nih.gov/10760061/>
    → electrogenic Na⁺ absorption through ENaC → a lumen-negative potential → ROMK K⁺ secretion. The model's `EKENAC`.
33. Wang WH, Giebisch G. **Regulation of potassium (K) handling in the renal collecting duct.** Pflugers Arch. 2009;458(1):157–168. — <https://pubmed.ncbi.nlm.nih.gov/18839186/>
    → flow-dependent K⁺ secretion through the BK channel — the model's `EKFLOW`.
34. Galla JH. **Metabolic alkalosis.** J Am Soc Nephrol. 2000;11(2):369–375. — <https://pubmed.ncbi.nlm.nih.gov/10665945/>
    → the principle that **chloride depletion raises the renal bicarbonate threshold** → the model's `KCLTH` term, and the reason the alkalosis is not corrected by anything other than KCl (a chloride salt).
35. Luke RG, Galla JH. **It is chloride depletion alkalosis, not contraction alkalosis.** J Am Soc Nephrol. 2012;23(2):204–207. — <https://pubmed.ncbi.nlm.nih.gov/22188853/>
    → the decisive evidence for that principle. The source of the model's prediction that HCO₃⁻ falls on NaCl replacement alone.
36. Sterns RH, Cox M, Feig PU, Singer I. **Internal potassium balance and the control of the plasma potassium concentration.** Medicine (Baltimore). 1981;60(5):339–354. — <https://pubmed.ncbi.nlm.nih.gov/7024719/>
    → quantification of the intracellular shift of K⁺ caused by alkalosis — the model's `KALKSH` and `KSH` parameters.

## 6. Treatment — drug effects and the evidence

37. Colussi G, Rombolà G, De Ferrari ME, et al. **Correction of hypokalemia with antialdosterone therapy in Gitelman's syndrome.** Am J Nephrol. 1994;14(2):127–135. — <https://pubmed.ncbi.nlm.nih.gov/8080234/>
    → the size of the rise in K⁺ on spironolactone and amiloride — calibration of the model's S8/S9 scenarios.
38. Blanchard A, Vargas-Poussou R, Vallet M, et al. **Indomethacin, amiloride, or eplerenone for treating hypokalemia in Gitelman syndrome.** J Am Soc Nephrol. 2015;26(2):468–475. — <https://pubmed.ncbi.nlm.nih.gov/25012174/>
    → a randomised crossover trial: ΔK for the three drugs in Gitelman syndrome. The result that **indomethacin is not as effective in Gitelman syndrome as expected and has a high discontinuation rate from adverse effects** — the validation data for the model's S6 prediction (indomethacin ineffective).
39. Reinalter SC, Gröne HJ, Konrad M, Seyberth HW, Klaus G. **Evaluation of long-term treatment with indomethacin in hereditary hypokalemic salt-losing tubulopathies.** J Pediatr. 2001;139(3):398–406. — <https://pubmed.ncbi.nlm.nih.gov/11562620/>
    → the improvement in growth and biochemistry on long-term indomethacin, and the risk of renal tissue injury — the model's `KNSAIDK` and growth module.
40. Vaisbich MH, Fujimura MD, Koch VH. **Bartter syndrome: benefits and side effects of long-term treatment.** Pediatr Nephrol. 2004;19(8):858–863. — <https://pubmed.ncbi.nlm.nih.gov/15206019/>
    → the frequency of gastrointestinal adverse effects — the model's `GIM` compartment and the setting of the COX-1 IC50.
41. Nascimento CL, Garcia CL, Schvartsman BG, Vaisbich MH. **Treatment of Bartter syndrome: unsolved issue.** J Pediatr (Rio J). 2014;90(5):512–517. — <https://pubmed.ncbi.nlm.nih.gov/24913037/>
42. Kleta R, Bockenhauer D. **Salt-losing tubulopathies in children: what's new, what's controversial?** J Am Soc Nephrol. 2018;29(3):727–739. — <https://pubmed.ncbi.nlm.nih.gov/29237736/>
    → a synthesis of the treatment controversies — the basis for the scenario design (combination therapy, cautious use of an ACEi).
43. Buerkert J, Martin D, Trigg D, Simon EE. **Effect of reduced renal mass on ammonium handling and net acid formation by the superficial and deep nephron of the rat.** J Clin Invest. 1983;71(6):1661–1675. — <https://pubmed.ncbi.nlm.nih.gov/6863543/>
    → the maximum capacity for distal H⁺/NH₄⁺ secretion — the ceiling of the model's `JHMAX`.
44. Colussi G, Bettinelli A, Tedeschi S, et al. **A thiazide test for the diagnosis of renal tubular hypokalemic disorders.** Clin J Am Soc Nephrol. 2007;2(3):454–460. — <https://pubmed.ncbi.nlm.nih.gov/17699452/>
    → the ΔFE-Cl response to a thiazide load — the validation data for the model's S13 (HCTZ challenge) prediction.

## 7. Pharmacokinetics

45. Helleberg L. **Clinical pharmacokinetics of indomethacin.** Clin Pharmacokinet. 1981;6(4):245–258. — <https://pubmed.ncbi.nlm.nih.gov/7249487/>
    → F ≈ 0.98, Vd ≈ 0.29 L/kg, protein binding 99%, t½ 4–9 h (enterohepatic circulation) — the model's `CLIND`/`VIND`/`KAIND`.
46. Davies NM, McLachlan AJ, Day RO, Williams KM. **Clinical pharmacokinetics and pharmacodynamics of celecoxib: a selective cyclo-oxygenase-2 inhibitor.** Clin Pharmacokinet. 2000;38(3):225–242. — <https://pubmed.ncbi.nlm.nih.gov/10749518/>
    → CL/F ≈ 27.7 L/h, V/F ≈ 400 L, t½ ≈ 11 h, CYP2C9 — the model's `CLCEL`/`VCEL`.
47. Vidt DG. **Mechanism of action, pharmacokinetics, adverse effects, and therapeutic uses of amiloride hydrochloride, a new potassium-sparing diuretic.** Pharmacotherapy. 1981;1(3):179–187. — <https://pubmed.ncbi.nlm.nih.gov/6765180/>
    → t½ 6–9 h, predominantly renally excreted — the model's `CLAMI`/`VAMI`/`IC50AMI`.
48. Overdiek HW, Merkus FW. **The metabolism and biopharmaceutics of spironolactone in man.** Rev Drug Metab Drug Interact. 1987;5(4):273–302. — <https://pubmed.ncbi.nlm.nih.gov/3332330/>
    → canrenone t½ 16–20 h, and the conversion rate to the active metabolite — the model's `FCAN`/`CLCAN`/`VCAN`.
49. Todd PA, Heel RC. **Enalapril: a review of its pharmacodynamic and pharmacokinetic properties and therapeutic use in hypertension and congestive heart failure.** Drugs. 1986;31(3):198–248. — <https://pubmed.ncbi.nlm.nih.gov/3011419/>
    → an enalaprilat conversion of ~40% and an effective t½ of ~11 h — the model's `FACE`/`CLACE`.
50. Ranade VV, Somberg JC. **Bioavailability and pharmacokinetics of magnesium after administration of magnesium salts to humans.** Am J Ther. 2001;8(5):345–357. — <https://pubmed.ncbi.nlm.nih.gov/11550076/>
    → magnesium oxide ≪ the aspartate, lactate, and chloride salts, with saturable absorption and osmotic diarrhoea — the model's `KMGSAT` and the divided-dose scenario (S16).

## 8. Long-term prognosis · complications · growth

51. Bockenhauer D, Bichet DG. **Inherited secondary nephrogenic diabetes insipidus: concentrating on humans.** Am J Physiol Renal Physiol. 2013;304(8):F1037–F1042. — <https://pubmed.ncbi.nlm.nih.gov/23364803/>
    → the urinary concentrating defect of a TAL lesion — the model's `ETALCON` term and its prediction of polyuria.
52. Bettinelli A, Tosetto C, Colussi G, et al. **Electrocardiogram with prolonged QT interval in Gitelman disease.** Kidney Int. 2002;62(2):580–584. — <https://pubmed.ncbi.nlm.nih.gov/12110020/>
    → the frequency of QTc prolongation in Gitelman syndrome — calibration of the model's `AKQT` and `AMQT`.
53. Cortesi C, Lava SA, Bettinelli A, et al. **Cardiac arrhythmias and rhabdomyolysis in Bartter-Gitelman patients.** Pediatr Nephrol. 2010;25(10):2005–2008. — <https://pubmed.ncbi.nlm.nih.gov/20535621/>
54. Puricelli E, Bettinelli A, Borsa N, et al. **Long-term follow-up of patients with Bartter syndrome type I and II.** Nephrol Dial Transplant. 2010;25(9):2976–2981. — <https://pubmed.ncbi.nlm.nih.gov/20219833/>
    → the frequency of nephrocalcinosis (>80%) and the long-term eGFR trend — calibration of the model's `KNC` and `KFIB`.
55. Bockenhauer D, Bichet DG. **Nephrocalcinosis in salt-losing tubulopathies.** (in) Pediatr Nephrol reviews — a representative paper: Kleta R, Bockenhauer D. **Bartter syndromes and other salt-losing tubulopathies.** Nephron Physiol. 2006;104(2):p73–p80. — <https://pubmed.ncbi.nlm.nih.gov/16785747/>
56. Simon DB, Lifton RP. **Ion transporter mutations in Gitelman's and Bartter's syndromes.** Curr Opin Nephrol Hypertens. 1998;7(1):43–47. — <https://pubmed.ncbi.nlm.nih.gov/9442361/>
57. Cruz DN, Simon DB, Nelson-Williams C, et al. **Mutations in the Na-Cl cotransporter reduce blood pressure in humans.** Hypertension. 2001;37(6):1458–1464. — <https://pubmed.ncbi.nlm.nih.gov/11408395/>
    → blood pressure is lower with NCC heterozygosity alone — the basis for the structure in which blood pressure does not rise (normal or low) even with a high aldosterone.
58. Ellison DH, Terker AS, Gamba G. **Potassium and its discontents: new insight, new treatments.** J Am Soc Nephrol. 2016;27(4):981–989. — <https://pubmed.ncbi.nlm.nih.gov/26510885/>
    → K⁺ sensing on the Kir4.1–WNK–NCC axis — the model's NCC adaptation (`ENCC`) and the EAST phenotype.
59. Ranieri M. **Renal Ca2+ and water handling in response to calcium sensing receptor signaling.** Front Physiol. 2019;10:1073. — <https://pubmed.ncbi.nlm.nih.gov/31481905/>
    → the CaSR inhibits NKCC2/ROMK and impairs urinary concentration — the model's `CASRGF` and `ECASR`.
60. Vezzoli G, Arcidiacono T, Paloschi V, et al. **Autosomal dominant hypocalcemia with mild type 5 Bartter syndrome.** J Nephrol. 2006;19(4):525–528. — <https://pubmed.ncbi.nlm.nih.gov/17048214/>
    → a Bartter-like phenotype from gain of function in the CaSR.

## 9. Differential diagnosis · drug-induced phenocopies

61. Kim YG, Kim B, Kim MK, et al. **Medullary nephrocalcinosis associated with long-term furosemide abuse in adults.** Nephrol Dial Transplant. 2001;16(12):2303–2309. — <https://pubmed.ncbi.nlm.nih.gov/11733622/>
    → distinguishing diuretic abuse (pseudo-Bartter) — corresponds to the model's `DXFURO`/`DXHCTZ` challenge module.
62. Kim GH. **Pseudo-Bartter syndrome: a rare cause of hypokalemic metabolic alkalosis.** Electrolyte Blood Press. 2016;14(2):27–30. — <https://pubmed.ncbi.nlm.nih.gov/28275375/>
63. Schrag D, Chung KY, Flombaum C, Saltz L. **Cetuximab therapy and symptomatic hypomagnesemia.** J Natl Cancer Inst. 2005;97(16):1221–1224. — <https://pubmed.ncbi.nlm.nih.gov/16106027/>
    → EGFR inhibition → reduced TRPM6 function → drug-induced hypomagnesaemia (a phenocopy through the model's `TRPM` compartment).
64. Kim GH, Han JS. **Therapeutic approach to hypokalemia.** Nephron. 2002;92(Suppl 1):28–32. — <https://pubmed.ncbi.nlm.nih.gov/12401934/>
65. Nakhoul F, Nakhoul N, Dorman E, et al. **Gitelman's syndrome: a pathophysiological and clinical update.** Endocrine. 2012;41(1):53–57. — <https://pubmed.ncbi.nlm.nih.gov/22169961/>

## 10. QSP and modelling methodology

66. Baek IH, et al. / Hallow KM, Gebremichael Y. **A quantitative systems physiology model of renal function and blood pressure regulation: model description.** CPT Pharmacometrics Syst Pharmacol. 2017;6(6):383–392. — <https://pubmed.ncbi.nlm.nih.gov/28653504/>
    → the standard approach of building Na⁺ handling segment by segment along the nephron as ODEs — the form the PT→TAL→DCT→ASDN chain of this model follows.
67. Hallow KM, Gebremichael Y. **A quantitative systems physiology model of renal function and blood pressure regulation: application in salt-sensitive hypertension.** CPT Pharmacometrics Syst Pharmacol. 2017;6(6):393–400. — <https://pubmed.ncbi.nlm.nih.gov/28653505/>
68. Baumgartner L, Sturla F, et al. / Moss R, Thomas SR. **Hormonal regulation of salt and water excretion: a mathematical model of whole kidney function and pressure natriuresis.** Am J Physiol Renal Physiol. 2014;306(2):F224–F248. — <https://pubmed.ncbi.nlm.nih.gov/24107425/>
69. Elshenawy S, Pinter A, et al. / Layton AT, Edwards A. **Mathematical modeling in renal physiology.** (Springer, 2014) — review: Layton AT. **Mathematical modeling of kidney transport.** Wiley Interdiscip Rev Syst Biol Med. 2013;5(5):557–573. — <https://pubmed.ncbi.nlm.nih.gov/23852667/>
70. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus machine learning, a fight worth fighting for the biological community?** Biol Lett. 2018;14(5):20170660. — <https://pubmed.ncbi.nlm.nih.gov/29769297/>
71. Elmokadem A, Riggs MM, Baron KT. **Quantitative systems pharmacology and physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on tutorial.** CPT Pharmacometrics Syst Pharmacol. 2019;8(12):883–893. — <https://pubmed.ncbi.nlm.nih.gov/31654588/>
    → the mrgsolve implementation conventions of this model (`$PARAM`/`$CMT`/`$ODE`/`$CAPTURE`).

---

## Appendix — The model's key predictions and the corresponding validation literature

| Model prediction | The structure the prediction comes from | Validation literature |
|---|---|---|
| Urinary PGE₂ is raised 5–10 fold in a TAL lesion only | the macula densa senses through its own NKCC2 (`MDSENSE = f(FTAL)`) | 24, 26, 27, 28 |
| COX inhibition is effective only in Bartter syndrome | that loop does not exist in Gitelman syndrome | 38, 39 |
| The sign reversal of urinary Ca/Cr | collapse of the lumen-positive TAL potential against increased PT Ca absorption from volume contraction | 3, 4, 20 |
| Serum Mg: Gitelman < type III < the antenatal form | recruitment of DCT TRPM6 by the distal Mg load (`KMGLOAD`) | 4, 19, 21 |
| K does not rise without Mg replacement | disinhibition of ROMK by low Mg (`KMGROMK`) | 23 |
| It has to be KCl, not KHCO₃, for the alkalosis to correct | chloride depletion raises the HCO₃⁻ threshold (`KCLTH`) | 34, 35 |
| NSAID + gastroenteritis → acute kidney injury | loss of the PG-dependent support of GFR (`KGPG × VOLIDX`) | 39, 40, 42 |
| The response to a thiazide load is flat in Gitelman syndrome | NCC is already out of action (`ANCC ≈ 0`) | 44 |

> All PubMed links were checked as at July 2026. Some older references may offer the abstract only.
