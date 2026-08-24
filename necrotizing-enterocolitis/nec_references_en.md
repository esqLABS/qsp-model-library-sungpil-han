# Necrotising Enterocolitis (NEC) — References
# Necrotising Enterocolitis — Reference List for the QSP Model

> **A note on how these links are built.**
> Most of the entries below are **title-based PubMed search links**. That is because a link looked up by
> title always points at the right paper, which is better than inventing a PMID that does
> not exist. For the handful of landmark papers we are certain of, the PMID is given as well.
>
> Most entries below are **title-resolved PubMed search links** rather than
> bare PMIDs. A title query always resolves to the correct paper, whereas a
> mis-remembered accession number silently points at the wrong one. PMIDs are
> given only for the handful of landmark papers where the number is certain.

The square brackets after each section title indicate **which term of the model that literature supports**
— read this not as a bibliography but as the evidence ledger of the model.
The bracket after each section heading names **which term of the model that
section supports**, so this reads as a provenance ledger rather than a booklist.

---

## 1. Definition · staging · epidemiology
*[model term: the Bell staging read-out rule, the calibration targets for incidence by GA]*

1. Bell MJ, Ternberg JL, Feigin RD, et al. **Neonatal necrotizing enterocolitis:
   therapeutic decisions based upon clinical staging.** *Ann Surg.* 1978;187(1):1-7.
   (PMID 413500) —
   <https://pubmed.ncbi.nlm.nih.gov/413500/>
   → the origin of the model's `BellStage` read-out rule (suspected / confirmed / advanced).
2. Walsh MC, Kliegman RM. **Necrotizing enterocolitis: treatment based on staging
   criteria.** *Pediatr Clin North Am.* 1986. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Walsh+Kliegman+necrotizing+enterocolitis+treatment+based+on+staging+criteria>
3. Neu J, Walker WA. **Necrotizing enterocolitis.** *N Engl J Med.*
   2011;364(3):255-264. (PMID 21247316) —
   <https://pubmed.ncbi.nlm.nih.gov/21247316/>
4. Battersby C, Santhalingam T, Costeloe K, Modi N. **Incidence of neonatal
   necrotising enterocolitis in high-income countries: a systematic review.**
   *Arch Dis Child Fetal Neonatal Ed.* 2018. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Battersby+incidence+of+neonatal+necrotising+enterocolitis+in+high-income+countries>
   → the basis for the calibration targets for incidence by GA band (24-25 weeks 12-15 %, 30-32 weeks 1-3 %).
5. Patel RM, Kandefer S, Walsh MC, et al. **Causes and timing of death in
   extremely premature infants from 2000 through 2011.** *N Engl J Med.* 2015. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Patel+causes+and+timing+of+death+in+extremely+premature+infants+2000+through+2011>
6. Berman L, Moss RL. **Necrotizing enterocolitis: an update.**
   *Semin Fetal Neonatal Med.* 2011. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Berman+Moss+necrotizing+enterocolitis+an+update>
7. Sharma R, Hudak ML. **A clinical perspective of necrotizing enterocolitis:
   past, present, and future.** *Clin Perinatol.* 2013. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Sharma+Hudak+clinical+perspective+necrotizing+enterocolitis+past+present+future>
8. Gephart SM, Spitzer AR, Effken JA, et al. **Discrimination of GutCheck-NEC:
   a clinical risk index for necrotizing enterocolitis.** *J Perinatol.* 2014. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Gephart+GutCheck+NEC+clinical+risk+index+necrotizing+enterocolitis>

---

## 2. Immature-gut TLR4 over-expression
*[model term: `TLR4expr = TLR4max·(1 − phiTLR·MAT)` — why the more immature gut responds more strongly to the same organism]*

9. Nanthakumar N, Meng D, Goldstein AM, et al. **The mechanism of excessive
   intestinal inflammation in necrotizing enterocolitis: an immature innate
   immune response.** *PLoS One.* 2011. —
   <https://pubmed.ncbi.nlm.nih.gov/?term=Nanthakumar+mechanism+of+excessive+intestinal+inflammation+necrotizing+enterocolitis+immature+innate+immune+response>
   → direct evidence that IκB expression is low and the NF-κB response excessive in the immature gut.
10. Sodhi CP, Neal MD, Siggers R, et al. **Intestinal epithelial Toll-like
    receptor 4 regulates goblet cell development and is required for
    necrotizing enterocolitis in mice.** *Gastroenterology.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sodhi+intestinal+epithelial+Toll-like+receptor+4+regulates+goblet+cell+development+necrotizing+enterocolitis>
11. Leaphart CL, Cavallo J, Gribar SC, et al. **A critical role for TLR4 in the
    pathogenesis of necrotizing enterocolitis by modulating intestinal injury
    and repair.** *J Immunol.* 2007. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Leaphart+critical+role+for+TLR4+pathogenesis+necrotizing+enterocolitis+modulating+intestinal+injury+and+repair>
    → the basis for the **proliferation-blocking term `block = 1/(1+(INJ/Ki)²)`**: TLR4 signalling does not only
    increase injury, it suppresses crypt regeneration at the same time.
12. Hackam DJ, Sodhi CP. **Toll-like receptor-mediated intestinal inflammatory
    imbalance in the pathogenesis of necrotizing enterocolitis.**
    *Cell Mol Gastroenterol Hepatol.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hackam+Sodhi+Toll-like+receptor+mediated+intestinal+inflammatory+imbalance+pathogenesis+necrotizing+enterocolitis>
13. Niño DF, Sodhi CP, Hackam DJ. **Necrotizing enterocolitis: new insights into
    pathogenesis and mechanisms.** *Nat Rev Gastroenterol Hepatol.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nino+Sodhi+Hackam+necrotizing+enterocolitis+new+insights+into+pathogenesis+and+mechanisms>
14. Good M, Sodhi CP, Hackam DJ. **Evidence-based feeding strategies before and
    after the development of necrotizing enterocolitis.**
    *Expert Rev Clin Immunol.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Good+Sodhi+Hackam+evidence-based+feeding+strategies+before+and+after+development+necrotizing+enterocolitis>
15. Sampah MES, Hackam DJ. **Dysregulated mucosal immunity and associated
    pathogeneses in preterm neonates.** *Front Immunol.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sampah+Hackam+dysregulated+mucosal+immunity+associated+pathogeneses+preterm+neonates>
16. Lu P, Sodhi CP, Hackam DJ. **Toll-like receptor regulation of intestinal
    development and inflammation in the pathogenesis of necrotizing
    enterocolitis.** *Pathophysiology.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Lu+Sodhi+Hackam+Toll-like+receptor+regulation+intestinal+development+inflammation+necrotizing+enterocolitis>

---

## 3. The gut microbial community and dysbiosis
*[model term: `B`, `C`, the `alphaX` cross-competition, `fpH`, the two stable states]*

17. Warner BB, Deych E, Zhou Y, et al. **Gut bacteria dysbiosis and necrotising
    enterocolitis in very low birthweight infants: a prospective case-control
    study.** *Lancet.* 2016;387(10031):1928-1936. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Warner+gut+bacteria+dysbiosis+and+necrotising+enterocolitis+very+low+birthweight+prospective+case-control>
    → before NEC the relative abundance of Gammaproteobacteria (Enterobacteriaceae) rises and
    obligate anaerobes fall — corresponding to the `B`-dominant state of the model.
18. Pammi M, Cope J, Tarr PI, et al. **Intestinal dysbiosis in preterm infants
    preceding necrotizing enterocolitis: a systematic review and
    meta-analysis.** *Microbiome.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Pammi+intestinal+dysbiosis+preterm+infants+preceding+necrotizing+enterocolitis+systematic+review+meta-analysis>
19. Olm MR, Bhattacharya N, Crits-Christoph A, et al. **Necrotizing
    enterocolitis is preceded by increased gut bacterial replication, Klebsiella,
    and fimbriae-encoding bacteria.** *Sci Adv.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Olm+necrotizing+enterocolitis+preceded+by+increased+gut+bacterial+replication+Klebsiella+fimbriae>
    → why the model handles the *growth rate* of `B` (`muB`) together with its load.
20. Stewart CJ, Ajami NJ, O'Brien JL, et al. **Temporal development of the gut
    microbiome in early childhood from the TEDDY study.** *Nature.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Stewart+temporal+development+of+the+gut+microbiome+in+early+childhood+TEDDY>
21. Claud EC, Walker WA. **Hypothesis: inappropriate colonization of the
    premature intestine can cause neonatal necrotizing enterocolitis.**
    *FASEB J.* 2001. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Claud+Walker+hypothesis+inappropriate+colonization+premature+intestine+neonatal+necrotizing+enterocolitis>
22. Ward DV, Scholz M, Zolfo M, et al. **Metagenomic sequencing with
    strain-level resolution implicates uropathogenic E. coli in necrotizing
    enterocolitis and mortality in preterm infants.** *Cell Rep.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ward+metagenomic+sequencing+strain-level+resolution+uropathogenic+coli+necrotizing+enterocolitis+preterm>
23. Heida FH, van Zoonen AGJF, Hulscher JBF, et al. **A necrotizing
    enterocolitis-associated gut microbiota is present in the meconium:
    results of a prospective study.** *Clin Infect Dis.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Heida+necrotizing+enterocolitis+associated+gut+microbiota+present+in+the+meconium+prospective>
24. Berrington JE, Stewart CJ, Cummings SP, Embleton ND. **The neonatal
    bowel microbiome in health and infection.** *Curr Opin Infect Dis.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Berrington+neonatal+bowel+microbiome+in+health+and+infection>

---

## 4. Human milk, HMOs, and the strain that can eat them
*[model term: `HMO` → `Ktlr` (the threshold), `muCH·binfant` (ecology), `eEGF` (nutrition), `SIGA`]*

25. Lucas A, Cole TJ. **Breast milk and neonatal necrotising enterocolitis.**
    *Lancet.* 1990;336(8730):1519-1523. (PMID 1979363) —
    <https://pubmed.ncbi.nlm.nih.gov/1979363/>
    → the first large-scale study to quantify the protective effect of human milk.
26. Meinzen-Derr J, Poindexter B, Wrage L, et al. **Role of human milk intake
    and persistent organ failure on the risk of necrotizing enterocolitis or
    death.** *J Perinatol.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Meinzen-Derr+role+of+human+milk+intake+risk+of+necrotizing+enterocolitis+or+death>
27. Quigley M, Embleton ND, McGuire W. **Formula versus donor breast milk for
    feeding preterm or low birth weight infants.** *Cochrane Database Syst Rev.*
    2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Quigley+Embleton+McGuire+formula+versus+donor+breast+milk+for+feeding+preterm+or+low+birth+weight+infants>
    → the basis for the model's `fDM` arm (pasteurisation loses sIgA and EGF but preserves the HMOs).
28. Cristofalo EA, Schanler RJ, Blanco CL, et al. **Randomized trial of
    exclusive human milk versus preterm formula diets in extremely premature
    infants.** *J Pediatr.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Cristofalo+randomized+trial+of+exclusive+human+milk+versus+preterm+formula+diets+extremely+premature>
29. Autran CA, Kellman BP, Kim JH, et al. **Human milk oligosaccharide
    composition predicts risk of necrotising enterocolitis in preterm
    infants.** *Gut.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Autran+human+milk+oligosaccharide+composition+predicts+risk+of+necrotising+enterocolitis+preterm>
    → the concentration of one particular HMO (DSLNT) is associated with NEC risk — the basis for HMOs entering the
    model as the term that raises the **threshold** (`hHMO`).
30. Masi AC, Embleton ND, Lamb CA, et al. **Human milk oligosaccharide DSLNT
    and gut microbiome in preterm infants predicts necrotising
    enterocolitis.** *Gut.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Masi+human+milk+oligosaccharide+DSLNT+and+gut+microbiome+preterm+infants+predicts+necrotising+enterocolitis>
    → the observation that HMOs and the microbial community determine risk **together** — why the model puts human milk
    into **both** the threshold shift and the state shift.
31. Jantscher-Krenn E, Zherebtsov M, Nissan C, et al. **The human milk
    oligosaccharide disialyllacto-N-tetraose prevents necrotising enterocolitis
    in neonatal rats.** *Gut.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Jantscher-Krenn+human+milk+oligosaccharide+disialyllacto-N-tetraose+prevents+necrotising+enterocolitis+neonatal+rats>
32. Good M, Sodhi CP, Yamaguchi Y, et al. **The human milk oligosaccharide
    2'-fucosyllactose attenuates the severity of experimental necrotising
    enterocolitis by enhancing mesenteric perfusion in the neonatal
    intestine.** *Br J Nutr.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Good+human+milk+oligosaccharide+2-fucosyllactose+attenuates+severity+experimental+necrotising+enterocolitis+mesenteric+perfusion>
    → evidence that 2′-FL also acts through **intestinal blood flow** — connects to the `PERF` pathway of the model.
33. Bode L. **Human milk oligosaccharides: every baby needs a sugar mama.**
    *Glycobiology.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bode+human+milk+oligosaccharides+every+baby+needs+a+sugar+mama>
34. Underwood MA, German JB, Lebrilla CB, Mills DA. **Bifidobacterium longum
    subspecies infantis: champion colonizer of the infant gut.**
    *Pediatr Res.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Underwood+Bifidobacterium+longum+subspecies+infantis+champion+colonizer+of+the+infant+gut>
    → the basis for the model's `binfant` switch (whether or not a strain that can actually metabolise HMOs is
    present). Even where HMOs are there, without an organism able to use them they are simply excreted.
35. Nolan LS, Rimer JM, Good M. **The role of human milk oligosaccharides and
    probiotics on the neonatal microbiome and risk of necrotizing
    enterocolitis.** *Nutrients.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nolan+role+of+human+milk+oligosaccharides+and+probiotics+neonatal+microbiome+risk+of+necrotizing+enterocolitis>
36. Watkins DJ, Besner GE. **The role of the intestinal microcirculation in
    necrotizing enterocolitis.** *Semin Pediatr Surg.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Watkins+Besner+role+of+the+intestinal+microcirculation+in+necrotizing+enterocolitis>
37. Feng J, El-Assal ON, Besner GE. **Heparin-binding EGF-like growth factor
    (HB-EGF) and necrotizing enterocolitis.** *Semin Pediatr Surg.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Feng+Besner+heparin-binding+EGF-like+growth+factor+HB-EGF+and+necrotizing+enterocolitis>
    → the model's `eEGF` term (growth factors in human milk raise crypt proliferation).

---

## 5. Probiotics — the pure state-mover
*[model term: `P`, `prob_dose`, `EmxPamp/gen/mtz` — why it does not work when given together with an antibiotic]*

38. Sharif S, Meader N, Oddie SJ, Rojas-Reyes MX, McGuire W. **Probiotics to
    prevent necrotising enterocolitis in very preterm or very low birth weight
    infants.** *Cochrane Database Syst Rev.* 2023. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sharif+probiotics+to+prevent+necrotising+enterocolitis+very+preterm+or+very+low+birth+weight+infants+Cochrane>
    → the calibration target for the effect size of the probiotics arm of the model.
39. Jacobs SE, Tobin JM, Opie GF, et al. **Probiotic effects on late-onset
    sepsis in very preterm infants: a randomized controlled trial (ProPrems).**
    *Pediatrics.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Jacobs+probiotic+effects+on+late-onset+sepsis+in+very+preterm+infants+randomized+controlled+trial>
40. Costeloe K, Hardy P, Juszczak E, Wilks M, Millar MR. **Bifidobacterium
    breve BBG-001 in very preterm infants: a randomised controlled phase 3
    trial (PiPS).** *Lancet.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Costeloe+Bifidobacterium+breve+BBG-001+in+very+preterm+infants+randomised+controlled+phase+3+trial>
    → **a negative result** — reproduced in the model as a failure to colonise (`gP` < `kwash` + `killP`).
41. Dermyshi E, Wang Y, Yan C, et al. **The 'golden age' of probiotics: a
    systematic review and meta-analysis of randomized and observational studies
    in preterm infants.** *Neonatology.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Dermyshi+golden+age+of+probiotics+systematic+review+meta-analysis+randomized+observational+studies+preterm+infants>
42. van den Akker CHP, van Goudoever JB, Shamir R, et al. **Probiotics and
    preterm infants: a position paper by the ESPGHAN Committee on Nutrition and
    the ESPGHAN Working Group for Probiotics and Prebiotics.**
    *J Pediatr Gastroenterol Nutr.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=van+den+Akker+probiotics+and+preterm+infants+position+paper+ESPGHAN+Committee+on+Nutrition>

---

## 6. Empirical antibiotics and duration of exposure
*[model term: the `killB` vs `killC` asymmetry, `abx_days` — where the sign changes]*

43. Cotten CM, Taylor S, Stoll B, et al. **Prolonged duration of initial
    empirical antibiotic treatment is associated with increased rates of
    necrotizing enterocolitis and death for extremely low birth weight
    infants.** *Pediatrics.* 2009;123(1):58-66. (PMID 19117861) —
    <https://pubmed.ncbi.nlm.nih.gov/19117861/>
    → the observation the model's **sign reversal** result is trying to reproduce. An antibiotic reduces `B` at once, but
    the cost of losing `C` appears with a delay, so the two integrals cross somewhere.
44. Alexander VN, Northrup V, Bizzarro MJ. **Antibiotic exposure in the newborn
    intensive care unit and the risk of necrotizing enterocolitis.**
    *J Pediatr.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Alexander+antibiotic+exposure+in+the+newborn+intensive+care+unit+and+the+risk+of+necrotizing+enterocolitis>
45. Greenwood C, Morrow AL, Lagomarcino AJ, et al. **Early empiric antibiotic
    use in preterm infants is associated with lower bacterial diversity and
    higher relative abundance of Enterobacter.** *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Greenwood+early+empiric+antibiotic+use+in+preterm+infants+lower+bacterial+diversity+higher+relative+abundance+Enterobacter>
    → direct evidence for the `killC` > `killB` asymmetry.
46. Ting JY, Synnes A, Roberts A, et al. **Association between antibiotic use
    and neonatal mortality and morbidities in very low-birth-weight infants
    without culture-proven sepsis or necrotizing enterocolitis.**
    *JAMA Pediatr.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ting+association+between+antibiotic+use+and+neonatal+mortality+and+morbidities+very+low-birth-weight+without+culture-proven+sepsis>
47. Rao SC, Srinivasjois R, Moon K. **One dose per day compared to multiple
    doses per day of gentamicin for treatment of suspected or proven sepsis in
    neonates.** *Cochrane Database Syst Rev.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Rao+one+dose+per+day+compared+to+multiple+doses+per+day+of+gentamicin+neonates+Cochrane>
    → the basis for the model's gentamicin q36h regimen and its two-compartment accumulation (`GENP` → `KIN`).

---

## 7. Feeding strategy and advance rate
*[model term: `feedrate`, `Ftroph`, `Sin`, `demand` — an input that carries both signs]*

48. Dorling J, Abbott J, Berrington J, et al. **Controlled trial of two
    incremental milk-feeding rates in preterm infants (SIFT).**
    *N Engl J Med.* 2019;381(15):1434-1443. (PMID 31597020) —
    <https://pubmed.ncbi.nlm.nih.gov/31597020/>
    → the large-scale result that the rate of advance (18 vs 30 mL/kg/d) made no difference to NEC.
    The model reproduces this as **a stretch where the slope flattens** in the more mature infant.
49. Leaf A, Dorling J, Kempley S, et al. **Early or delayed enteral feeding for
    preterm growth-restricted infants: a randomized trial (ADEPT).**
    *Pediatrics.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Leaf+early+or+delayed+enteral+feeding+for+preterm+growth-restricted+infants+randomized+trial>
50. Oddie SJ, Young L, McGuire W. **Slow advancement of enteral feed volumes to
    prevent necrotising enterocolitis in very low birth weight infants.**
    *Cochrane Database Syst Rev.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Oddie+slow+advancement+of+enteral+feed+volumes+to+prevent+necrotising+enterocolitis+very+low+birth+weight>
51. Morgan J, Young L, McGuire W. **Delayed introduction of progressive enteral
    feeds to prevent necrotising enterocolitis in very low birth weight
    infants.** *Cochrane Database Syst Rev.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Morgan+delayed+introduction+of+progressive+enteral+feeds+to+prevent+necrotising+enterocolitis>
52. Berseth CL, Bisquera JA, Paje VU. **Prolonging small feeding volumes early
    in life decreases the incidence of necrotizing enterocolitis in very low
    birth weight infants.** *Pediatrics.* 2003. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Berseth+prolonging+small+feeding+volumes+early+in+life+decreases+incidence+of+necrotizing+enterocolitis>
    → the basis for the direction the model actually reproduces (slower advance = less NEC).
53. Ellsbury DL, Ursprung R. **Comprehensive Oxygen Management for the
    Prevention of Retinopathy of Prematurity** *(and on NICU standardisation)*. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ellsbury+Ursprung+comprehensive+oxygen+management+prevention+retinopathy+of+prematurity>

---

## 8. Splanchnic perfusion, ischaemia, and PAF
*[model term: `PERF`, `demand`, `SDRcrit`, `PAF`, `PAFAH = dPAF0 + dPAFm·MAT`]*

54. Caplan MS, Sun XM, Hsueh W, Hageman JR. **Role of platelet activating factor
    and tumor necrosis factor-alpha in neonatal necrotizing enterocolitis.**
    *J Pediatr.* 1990. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Caplan+role+of+platelet+activating+factor+and+tumor+necrosis+factor+alpha+in+neonatal+necrotizing+enterocolitis>
    → the origin of the PAF axis. In the model PAF is **the only endogenous term that lowers blood flow**.
55. Caplan MS, Lickerman M, Adler L, et al. **The role of recombinant
    platelet-activating factor acetylhydrolase in a neonatal rat model of
    necrotizing enterocolitis.** *Pediatr Res.* 1997. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Caplan+role+of+recombinant+platelet-activating+factor+acetylhydrolase+neonatal+rat+model+necrotizing+enterocolitis>
    → PAF-AH deficiency in the preterm infant and the supplementation experiment — the basis for `dPAF0` (the floor value) and for `dPAFm·MAT`.
56. Yazji I, Sodhi CP, Robinson EK, et al. **Endothelial TLR4 activation impairs
    intestinal microcirculatory perfusion in necrotizing enterocolitis via
    eNOS-NO-nitrite signaling.** *Proc Natl Acad Sci U S A.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Yazji+endothelial+TLR4+activation+impairs+intestinal+microcirculatory+perfusion+necrotizing+enterocolitis+eNOS>
    → the TLR4 → reduced blood flow pathway. The structural basis for inflammation and ischaemia being **summed into a single `INJ`**
    in the model.
57. Nowicki PT. **Ischemia and necrotizing enterocolitis: where, when, and
    how.** *Semin Pediatr Surg.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Nowicki+ischemia+and+necrotizing+enterocolitis+where+when+and+how>
    → the rise in oxygen demand after a feed and the immaturity of autoregulation — the basis for the model's `demand = 1 + aDEM·FEED/200`
    and for the `SDRcrit` reserve concept.
58. Reber KM, Nankervis CA, Nowicki PT. **Newborn intestinal circulation:
    physiology and pathophysiology.** *Clin Perinatol.* 2002. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Reber+Nankervis+Nowicki+newborn+intestinal+circulation+physiology+and+pathophysiology>
59. Ohlsson A, Walia R, Shah SS. **Ibuprofen for the treatment of patent
    ductus arteriosus in preterm or low birth weight (or both) infants.**
    *Cochrane Database Syst Rev.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ohlsson+ibuprofen+for+the+treatment+of+patent+ductus+arteriosus+in+preterm+or+low+birth+weight+infants+Cochrane>
    → the difference between ibuprofen and indomethacin in their effect on splanchnic blood flow — the model's `EC50IND` vs `EC50IBU`.
60. Coombs RC, Morgan ME, Durbin GM, et al. **Gut blood flow velocities in the
    newborn: effects of patent ductus arteriosus and parenteral indomethacin.**
    *Arch Dis Child.* 1990. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Coombs+gut+blood+flow+velocities+in+the+newborn+effects+of+patent+ductus+arteriosus+and+parenteral+indomethacin>

---

## 9. Anaemia, transfusion, and TANEC
*[model term: the insult episodes (`INS*DEP`, `INS*NFK`) — the trigger that carries it over the watershed]*

61. Patel RM, Knezevic A, Shenvi N, et al. **Association of red blood cell
    transfusion, anemia, and necrotizing enterocolitis in very low-birth-weight
    infants.** *JAMA.* 2016;315(9):889-897. (PMID 26954412) —
    <https://pubmed.ncbi.nlm.nih.gov/26954412/>
    → anaemia is itself the risk, and transfusion does not immediately change that risk — the basis for putting the
    pre-transfusion anaemic stretch and the transfusion stretch into the model as **two different boxcars**.
62. Mohamed A, Shah PS. **Transfusion associated necrotizing enterocolitis: a
    meta-analysis of observational data.** *Pediatrics.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Mohamed+Shah+transfusion+associated+necrotizing+enterocolitis+meta-analysis+of+observational+data>
63. Marin T, Josephson CD, Kosmetatos N, et al. **Feeding preterm infants
    during red blood cell transfusion is associated with a decline in
    postprandial mesenteric oxygenation.** *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Marin+feeding+preterm+infants+during+red+blood+cell+transfusion+decline+in+postprandial+mesenteric+oxygenation>
    → feeding during a transfusion worsens the intestinal oxygen supply-demand balance — the `demand` × insult interaction of the model.
64. Kirpalani H, Bell EF, Hintz SR, et al. **Higher or lower hemoglobin
    transfusion thresholds for preterm infants (TOP trial).**
    *N Engl J Med.* 2020. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Kirpalani+higher+or+lower+hemoglobin+transfusion+thresholds+for+preterm+infants>

---

## 10. Barrier, tight junctions, mucus
*[model term: `TJ`, `MUC`, `BI`, `Pb = Pmin + (Pmax−Pmin)(1−BI)^hP`]*

65. Anderson JM, Van Itallie CM. **Physiology and function of the tight
    junction.** *Cold Spring Harb Perspect Biol.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Anderson+Van+Itallie+physiology+and+function+of+the+tight+junction>
66. Bein A, Eventov-Friedman S, Arbell D, Schwartz B. **Intestinal tight
    junctions are severely altered in NEC preterm neonates.**
    *Pediatr Neonatol.* 2018. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Bein+intestinal+tight+junctions+are+severely+altered+in+NEC+preterm+neonates>
67. McElroy SJ, Prince LS, Weitkamp JH, et al. **Tumor necrosis factor
    receptor 1-dependent depletion of mucus in immature small intestine: a
    potential role in neonatal necrotizing enterocolitis.**
    *Am J Physiol Gastrointest Liver Physiol.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=McElroy+tumor+necrosis+factor+receptor+1+dependent+depletion+of+mucus+in+immature+small+intestine+necrotizing+enterocolitis>
68. Halpern MD, Denning PW. **The role of intestinal epithelial barrier
    function in the development of NEC.** *Tissue Barriers.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Halpern+Denning+role+of+intestinal+epithelial+barrier+function+in+the+development+of+NEC>
69. Whitehouse JS, Riggle KM, Purpi DP, et al. **The protective role of
    intestinal alkaline phosphatase in necrotizing enterocolitis.**
    *J Surg Res.* 2010. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Whitehouse+protective+role+of+intestinal+alkaline+phosphatase+in+necrotizing+enterocolitis>
    → the IAP node of the model (an intervention that raises the signalling threshold by dephosphorylating LPS).
70. Underwood MA. **Paneth cells and necrotizing enterocolitis.**
    *Gut Microbes.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Underwood+Paneth+cells+and+necrotizing+enterocolitis>
    → the basis for the `DEF` (defensin) term.

---

## 11. Cytokines and cellular immunity
*[model term: `TNFA`, `IL1B`, `IL8`, `IL10`, `NEU` (saturating feedback), `NOX`]*

71. Egan CE, Sodhi CP, Good M, et al. **Toll-like receptor 4-mediated
    lymphocyte influx induces neonatal necrotizing enterocolitis.**
    *J Clin Invest.* 2016. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Egan+Toll-like+receptor+4+mediated+lymphocyte+influx+induces+neonatal+necrotizing+enterocolitis>
72. Weitkamp JH, Koyama T, Rock MT, et al. **Necrotising enterocolitis is
    characterised by disrupted immune regulation and diminished mucosal
    regulatory (FOXP3)/effector (CD4, CD8) T cell ratios.** *Gut.* 2013. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Weitkamp+necrotising+enterocolitis+disrupted+immune+regulation+diminished+mucosal+regulatory+FOXP3+effector+T+cell+ratios>
73. MohanKumar K, Namachivayam K, Song T, et al. **A murine neonatal model of
    necrotizing enterocolitis caused by anemia and red blood cell
    transfusions.** *Nat Commun.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=MohanKumar+murine+neonatal+model+of+necrotizing+enterocolitis+caused+by+anemia+and+red+blood+cell+transfusions>
74. Emami CN, Chokshi N, Wang J, et al. **Role of interleukin-10 in the
    pathogenesis of necrotizing enterocolitis.** *Am J Surg.* 2012. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Emami+role+of+interleukin-10+in+the+pathogenesis+of+necrotizing+enterocolitis>
    → `kIL10·MAT` — the basis for the term by which IL-10 counter-regulation is weaker the more immature the infant.
75. Chan KL, Ho JC, Chan KW, Tam PK. **A study of gut immunity to enteral
    endotoxin in rats of different ages.** *Pediatr Surg Int.* 2002. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Chan+study+of+gut+immunity+to+enteral+endotoxin+in+rats+of+different+ages>
76. Ford H, Watkins S, Reblock K, Rowe M. **The role of inflammatory cytokines
    and nitric oxide in the pathogenesis of necrotizing enterocolitis.**
    *J Pediatr Surg.* 1997. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ford+role+of+inflammatory+cytokines+and+nitric+oxide+in+the+pathogenesis+of+necrotizing+enterocolitis>
    → the basis for the `wNO·NOX` term (nitrosative injury).

---

## 12. Diagnostic biomarkers
*[model term: `CRP`, `PLT`, `LAC`, `PNEU`, the I-FABP and calprotectin read-out nodes]*

77. Thuijls G, Derikx JPM, van Wijck K, et al. **Non-invasive markers for early
    diagnosis and determination of the severity of necrotizing enterocolitis.**
    *Ann Surg.* 2010. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Thuijls+non-invasive+markers+for+early+diagnosis+and+determination+of+the+severity+of+necrotizing+enterocolitis>
78. Ng PC, Ma TP, Lam HS. **The use of laboratory biomarkers for surveillance,
    diagnosis and prediction of clinical outcomes in neonatal sepsis and
    necrotising enterocolitis.** *Arch Dis Child Fetal Neonatal Ed.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ng+use+of+laboratory+biomarkers+for+surveillance+diagnosis+prediction+clinical+outcomes+neonatal+sepsis+necrotising+enterocolitis>
79. Ragazzi S, Pierro A, Peters M, et al. **Early full blood count and severity
    of disease in neonates with necrotizing enterocolitis.**
    *Pediatr Surg Int.* 2003. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Ragazzi+early+full+blood+count+and+severity+of+disease+in+neonates+with+necrotizing+enterocolitis>
    → the relationship between thrombocytopenia and severity — the basis for `PLT < 100` entering the model's Bell III decision.
80. Christensen RD, Yoder BA, Baer VL, et al. **Early-onset neutropenia in
    small-for-gestational-age infants** *(and studies of the haematological findings in NEC)*. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Christensen+thrombocytopenia+necrotizing+enterocolitis+neonate>
81. Sylvester KG, Ling XB, Liu GY, et al. **Urine protein biomarkers for the
    diagnosis and prognosis of necrotizing enterocolitis in infants.**
    *J Pediatr.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Sylvester+urine+protein+biomarkers+for+the+diagnosis+and+prognosis+of+necrotizing+enterocolitis+in+infants>

---

## 13. Surgery, outcomes, sequelae
*[model term: `NECrev`, `NECa_crit = INJth/wNEC`, stricture · short bowel syndrome · NDI]*

82. Moss RL, Dimmitt RA, Barnhart DC, et al. **Laparotomy versus peritoneal
    drainage for necrotizing enterocolitis and perforation.**
    *N Engl J Med.* 2006;354(21):2225-2234. (PMID 16707747) —
    <https://pubmed.ncbi.nlm.nih.gov/16707747/>
83. Blakely ML, Tyson JE, Lally KP, et al. **Initial laparotomy versus
    peritoneal drainage in extremely low birthweight infants with surgical
    necrotizing enterocolitis or isolated intestinal perforation (NEST).**
    *Ann Surg.* 2021. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Blakely+initial+laparotomy+versus+peritoneal+drainage+extremely+low+birthweight+surgical+necrotizing+enterocolitis+isolated+intestinal+perforation>
84. Rees CM, Eaton S, Kiely EM, et al. **Peritoneal drainage or laparotomy for
    neonatal bowel perforation? A randomized controlled trial.**
    *Ann Surg.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Rees+peritoneal+drainage+or+laparotomy+for+neonatal+bowel+perforation+randomized+controlled+trial>
85. Hull MA, Fisher JG, Gutierrez IM, et al. **Mortality and management of
    surgical necrotizing enterocolitis in very low birth weight neonates: a
    prospective cohort study.** *J Am Coll Surg.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hull+mortality+and+management+of+surgical+necrotizing+enterocolitis+very+low+birth+weight+neonates+prospective+cohort>
86. Hintz SR, Kendrick DE, Stoll BJ, et al. **Neurodevelopmental and growth
    outcomes of extremely low birth weight infants after necrotizing
    enterocolitis.** *Pediatrics.* 2005. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hintz+neurodevelopmental+and+growth+outcomes+of+extremely+low+birth+weight+infants+after+necrotizing+enterocolitis>
    → the basis for the model's `NIN` → NDI read-out pathway.
87. Fredriksson F, Engstrand Lilja H. **Survival rates for surgically treated
    necrotising enterocolitis have improved over the last four decades.**
    *Acta Paediatr.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fredriksson+Engstrand+Lilja+survival+rates+for+surgically+treated+necrotising+enterocolitis+improved+over+the+last+four+decades>
88. Niemarkt HJ, De Meij TG, van Ganzewinkel CJ, et al. **Necrotizing
    enterocolitis, gut microbiota, and brain development: role of the
    brain-gut axis.** *Neonatology.* 2019. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Niemarkt+necrotizing+enterocolitis+gut+microbiota+and+brain+development+role+of+the+brain-gut+axis>
89. Robinson JR, Rellinger EJ, Hatch LD, et al. **Surgical necrotizing
    enterocolitis.** *Semin Perinatol.* 2017. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Robinson+surgical+necrotizing+enterocolitis+Semin+Perinatol>
90. Zani A, Pierro A. **Necrotizing enterocolitis: controversies and
    challenges.** *F1000Res.* 2015. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Zani+Pierro+necrotizing+enterocolitis+controversies+and+challenges>

---

## 14. Neonatal pharmacokinetics
*[model term: `CL·(PMA/40)^1.3`, the V and CL of each drug, `GENP` accumulation → `KIN`]*

91. Anderson BJ, Holford NHG. **Mechanism-based concepts of size and maturity
    in pharmacokinetics.** *Annu Rev Pharmacol Toxicol.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Anderson+Holford+mechanism-based+concepts+of+size+and+maturity+in+pharmacokinetics>
    → the basis for the way the model scales clearance by PMA (the maturation function).
92. Tremoulet A, Le J, Poindexter B, et al. **Characterization of ampicillin
    pharmacokinetics in neonates using a population approach.**
    *Antimicrob Agents Chemother.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Tremoulet+characterization+of+ampicillin+pharmacokinetics+in+neonates+using+a+population+approach>
93. Fuchs A, Guidi M, Giannoni E, et al. **Population pharmacokinetic study of
    gentamicin in a large cohort of premature and term neonates.**
    *Br J Clin Pharmacol.* 2014. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Fuchs+population+pharmacokinetic+study+of+gentamicin+in+a+large+cohort+of+premature+and+term+neonates>
94. Suyagh M, Collier PS, Millership JS, et al. **Metronidazole population
    pharmacokinetics in preterm neonates using dried blood-spot sampling.**
    *Pediatrics.* 2011. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Suyagh+metronidazole+population+pharmacokinetics+in+preterm+neonates+dried+blood-spot+sampling>
95. Smyth JM, Collier PS, Darwish M, et al. **Intravenous indometacin in
    preterm infants with symptomatic patent ductus arteriosus: a population
    pharmacokinetic study.** *Br J Clin Pharmacol.* 2004. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Smyth+intravenous+indometacin+in+preterm+infants+with+symptomatic+patent+ductus+arteriosus+population+pharmacokinetic>
96. Hirt D, Van Overmeire B, Treluyer JM, et al. **An optimized ibuprofen
    dosing scheme for preterm neonates with patent ductus arteriosus, based on
    a population pharmacokinetic and pharmacodynamic study.**
    *Br J Clin Pharmacol.* 2008. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Hirt+optimized+ibuprofen+dosing+scheme+for+preterm+neonates+with+patent+ductus+arteriosus+population+pharmacokinetic>

---

## 15. QSP method, bistability, critical transitions
*[model term: the number of roots of `g(E)`, the saddle-node bifurcations `B_lo`/`B_hi`, the logarithmic transition time of the watershed `E*`]*

97. Strogatz SH. **Nonlinear Dynamics and Chaos.** Westview Press. —
    the standard textbook account of saddle-node bifurcation and hysteresis.
    (a book — no PubMed entry)
98. Scheffer M, Bascompte J, Brock WA, et al. **Early-warning signals for
    critical transitions.** *Nature.* 2009. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Scheffer+early-warning+signals+for+critical+transitions>
    → the result that the rate of recovery slows just before a critical transition (critical slowing down).
    Corresponds directly to the eigenvalue λ ≈ 0.35/d near `E*` in the model and to the logarithmic transition time.
99. Angeli D, Ferrell JE, Sontag ED. **Detection of multistability,
    bifurcations, and hysteresis in a large class of biological positive-feedback
    systems.** *Proc Natl Acad Sci U S A.* 2004. —
    <https://pubmed.ncbi.nlm.nih.gov/?term=Angeli+Ferrell+Sontag+detection+of+multistability+bifurcations+and+hysteresis+biological+positive-feedback+systems>
100. Baker RE, Peña JM, Jayamohan J, Jérusalem A. **Mechanistic models versus
     machine learning, a fight worth fighting for the biological community?**
     *Biol Lett.* 2018. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Baker+mechanistic+models+versus+machine+learning+a+fight+worth+fighting+for+the+biological+community>
101. Elmokadem A, Riggs M, Baron K. **Quantitative systems pharmacology and
     physiologically-based pharmacokinetic modeling with mrgsolve: a hands-on
     tutorial.** *CPT Pharmacometrics Syst Pharmacol.* 2019. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Elmokadem+quantitative+systems+pharmacology+and+physiologically-based+pharmacokinetic+modeling+with+mrgsolve+hands-on+tutorial>
     → the mrgsolve implementation convention of this repository.

---

## 16. Comprehensive reviews used for cross-checking

102. Frost BL, Modi BP, Jaksic T, Caplan MS. **New medical and surgical insights
     into neonatal necrotizing enterocolitis.** *JAMA Pediatr.* 2017. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Frost+new+medical+and+surgical+insights+into+neonatal+necrotizing+enterocolitis>
103. Neu J. **Necrotizing enterocolitis: the future.** *Neonatology.* 2020. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Neu+necrotizing+enterocolitis+the+future+Neonatology>
104. Bazacliu C, Neu J. **Necrotizing enterocolitis: long term complications.**
     *Curr Pediatr Rev.* 2019. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Bazacliu+Neu+necrotizing+enterocolitis+long+term+complications>
105. Denning NL, Prince JM. **Neonatal intestinal dysbiosis in necrotizing
     enterocolitis.** *Mol Med.* 2018. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Denning+Prince+neonatal+intestinal+dysbiosis+in+necrotizing+enterocolitis>
106. Tanner SM, Berryhill TF, Ellenburg JL, et al. **Pathogenesis of
     necrotizing enterocolitis: modeling the innate immune response.**
     *Am J Pathol.* 2015. —
     <https://pubmed.ncbi.nlm.nih.gov/?term=Tanner+pathogenesis+of+necrotizing+enterocolitis+modeling+the+innate+immune+response>

---

## Predictions this model makes that the literature has not yet tested
## Predictions this model makes that the literature has NOT yet tested

To put it honestly, what follows **came not out of the literature above but out of the structure of the model**.
It is set out explicitly so that it can be falsified.

1. **The protective effect of human milk divides according to whether a strain able to metabolise HMOs (`binfant`)
   is present.** A study that measures only HMO concentration will see the effect diluted (references #29, #30, and #34
   point this way, but no study has yet tested the interaction directly).
2. **The harm of empirical antibiotics should be observable only in breast-milk-fed infants.** In an exclusively
   formula-fed infant, with no commensal flora to lose, an antibiotic may even look protective (testable by
   stratifying the cohort of reference #43 by feeding type).
3. **The point of no return is a function of lesion size, not a function of time.** In the model
   `NECa_crit = INJth/wNEC`, and whether the diagnosis-to-treatment delay (`dx_lag`) carries it over that threshold is what
   separates medical NEC from surgical NEC.
4. **A marker that sounds at Bell II can in principle buy almost no time.** Because the transition time near the
   watershed depends logarithmically, what moves first is the state variables (`E`, `Jtr`), and the
   clinical markers move only after the switch has already flipped.

---

**Disclaimer.** This reference list and this model are for educational and research purposes. The parameters were
calibrated by hand against **aggregate indices** from the published literature, and were not fitted to individual
patient data. They must not be used for clinical decision-making, for prescribing, or for regulatory submission.
