# Male- and Female-Pattern Hair Loss (Androgenetic Alopecia) — References

The literature underpinning this QSP model (mechanistic map · mrgsolve ODE · Shiny app).
**Every PMID was looked up individually through the NCBI E-utilities (`esearch` → `esummary`) to check that the author, year,
journal and title matched.** Papers used directly for the model's numerical calibration are marked ★.

---

## 1. Classification · epidemiology · natural history

1. Hamilton JB. Patterned loss of hair in man; types and incidence. *Ann N Y Acad Sci* 1951. — [PMID 14819896](https://pubmed.ncbi.nlm.nih.gov/14819896/)
2. Norwood OT. Male pattern baldness: classification and incidence. *South Med J* 1975. — [PMID 1188424](https://pubmed.ncbi.nlm.nih.gov/1188424/)
3. Ludwig E. Classification of the types of androgenetic alopecia (common baldness) occurring in the female sex. *Br J Dermatol* 1977. — [PMID 921894](https://pubmed.ncbi.nlm.nih.gov/921894/)
4. Sinclair R. Hair loss in women: medical and cosmetic approaches to increase scalp hair fullness. *Br J Dermatol* 2011. — [PMID 22171680](https://pubmed.ncbi.nlm.nih.gov/22171680/)
5. Cash TF. The psychosocial consequences of androgenetic alopecia: a review of the research literature. *Br J Dermatol* 1999. — [PMID 10583042](https://pubmed.ncbi.nlm.nih.gov/10583042/)
6. Ho CH, Sood T, Zito PM. Androgenetic Alopecia. *StatPearls*. — [PMID 28613674](https://pubmed.ncbi.nlm.nih.gov/28613674/)
7. Courtois M, Loussouarn G, Hourseau C, Grollier JF. Ageing and hair cycles. *Br J Dermatol* 1995. — [PMID 7756156](https://pubmed.ncbi.nlm.nih.gov/7756156/)
8. Ustuner ET. Baldness may be caused by the weight of the scalp: gravity as a proposed mechanism for hair loss. *Med Hypotheses* 2008. — [PMID 18667278](https://pubmed.ncbi.nlm.nih.gov/18667278/) *(an alternative hypothesis the model does not adopt; listed as something to be refuted)*

## 2. Genetics — the model's `GS` parameter

9. Richards JB, et al. Male-pattern baldness susceptibility locus at 20p11. *Nat Genet* 2008. — [PMID 18849991](https://pubmed.ncbi.nlm.nih.gov/18849991/)
10. Cobb JE, et al. Evidence for two independent functional variants for androgenetic alopecia around the androgen receptor gene. *Exp Dermatol* 2010. — [PMID 21073542](https://pubmed.ncbi.nlm.nih.gov/21073542/)
11. Heilmann-Heimbach S, et al. Meta-analysis identifies novel risk loci and yields systematic insights into the biology of male-pattern baldness. *Nat Commun* 2017. — [PMID 28272467](https://pubmed.ncbi.nlm.nih.gov/28272467/)
12. Sawaya ME, Shalita AR. Androgen receptor polymorphisms (CAG repeat lengths) in androgenetic alopecia, hirsutism, and acne. *J Cutan Med Surg* 1998. — [PMID 9677254](https://pubmed.ncbi.nlm.nih.gov/9677254/)

## 3. Androgen metabolism and the 5α-reductase isozymes

13. ★ Thigpen AE, et al. Tissue distribution and ontogeny of steroid 5α-reductase isozyme expression. *J Clin Invest* 1993. — [PMID 7688765](https://pubmed.ncbi.nlm.nih.gov/7688765/) — the basis of the model's `W_T1_LOC 0.30 / W_T2_LOC 0.70`
14. ★ Sawaya ME, Price VH. Different levels of 5α-reductase type I and II, aromatase, and androgen receptor in hair follicles of women and men with androgenetic alopecia. *J Invest Dermatol* 1997. — [PMID 9284093](https://pubmed.ncbi.nlm.nih.gov/9284093/) — `LOC5AR 0.30`, `AROM 1.6` in the female arm
15. Chen W, Zouboulis CC, et al. Cutaneous androgen metabolism: basic research and clinical perspectives. *J Invest Dermatol* 2002. — [PMID 12445184](https://pubmed.ncbi.nlm.nih.gov/12445184/)
16. Deplewski D, Rosenfield RL. Role of hormones in pilosebaceous unit development. *Endocr Rev* 2000. — [PMID 10950157](https://pubmed.ncbi.nlm.nih.gov/10950157/)
17. Randall VA. The hair follicle: a paradoxical androgen target organ. *Horm Res* 2000. — [PMID 11595812](https://pubmed.ncbi.nlm.nih.gov/11595812/)
18. Randall VA, et al. Androgen action in cultured dermal papilla cells from human hair follicles. *Skin Pharmacol* 1994. — [PMID 8003318](https://pubmed.ncbi.nlm.nih.gov/8003318/)
19. Hibberts NA, Howell AE, Randall VA. Balding hair follicle dermal papilla cells contain higher levels of androgen receptors than those from non-balding scalp. *J Endocrinol* 1998. — [PMID 9496234](https://pubmed.ncbi.nlm.nih.gov/9496234/)
20. Ryu HK, et al. Evaluation of androgens in the scalp hair and plasma of patients with male-pattern baldness before and after finasteride administration. *Br J Dermatol* 2006. — [PMID 16536818](https://pubmed.ncbi.nlm.nih.gov/16536818/)
21. ★ Dallob AL, et al. The effect of finasteride, a 5α-reductase inhibitor, on scalp skin testosterone and dihydrotestosterone concentrations. *J Clin Endocrinol Metab* 1994. — [PMID 8077349](https://pubmed.ncbi.nlm.nih.gov/8077349/)

## 4. Hair cycle biology — the model's 15 follicular compartments

22. Cotsarelis G, Sun TT, Lavker RM. Label-retaining cells reside in the bulge area of pilosebaceous unit. *Cell* 1990. — [PMID 2364430](https://pubmed.ncbi.nlm.nih.gov/2364430/)
23. Plikus MV, et al. Analyses of regenerative wave patterns in adult hair follicle populations. *Int J Dev Biol* 2009. — [PMID 19378257](https://pubmed.ncbi.nlm.nih.gov/19378257/)
24. ★ Whiting DA. Diagnostic and predictive value of horizontal sections of scalp biopsy specimens in male pattern androgenetic alopecia. *J Am Acad Dermatol* 1993. — [PMID 8496421](https://pubmed.ncbi.nlm.nih.gov/8496421/) — total follicle count and the distribution of miniaturisation
25. ★ Whiting DA. Possible mechanisms of miniaturization during androgenetic alopecia or pattern hair loss. *J Am Acad Dermatol* 2001. — [PMID 11511857](https://pubmed.ncbi.nlm.nih.gov/11511857/) — the direct basis of the model's "ratchet" structure
26. Guarrera M, Rebora A. Kenogen in female androgenetic alopecia. A longitudinal study. *Dermatology* 2005. — [PMID 15604539](https://pubmed.ncbi.nlm.nih.gov/15604539/)
27. ★ Van Neste D, et al. Finasteride increases anagen hair in men with androgenetic alopecia. *Br J Dermatol* 2000. — [PMID 11069460](https://pubmed.ncbi.nlm.nih.gov/11069460/) — change in the anagen fraction = the model's `T_A` term
28. Philpott MP, et al. Cultured human hair follicles and growth factors. *J Invest Dermatol* 1995. — [PMID 7738396](https://pubmed.ncbi.nlm.nih.gov/7738396/)

## 5. Dermal papilla paracrine signalling

29. ★ Kwack MH, et al. Dihydrotestosterone-inducible dickkopf 1 from balding dermal papilla cells causes apoptosis in follicular keratinocytes. *J Invest Dermatol* 2008. — [PMID 17657240](https://pubmed.ncbi.nlm.nih.gov/17657240/) — the model's `DKK1` state
30. Inui S, Itami S. Molecular basis of androgenetic alopecia: from androgen to paracrine mediators through dermal papilla. *J Dermatol Sci* 2011. — [PMID 21167691](https://pubmed.ncbi.nlm.nih.gov/21167691/)
31. Hibino T, Nishiyama T. Role of TGF-β2 in the human hair cycle. *J Dermatol Sci* 2004. — [PMID 15194142](https://pubmed.ncbi.nlm.nih.gov/15194142/)
32. Upton JH, et al. Oxidative stress-associated senescence in dermal papilla cells of men with androgenetic alopecia. *J Invest Dermatol* 2015. — [PMID 25647436](https://pubmed.ncbi.nlm.nih.gov/25647436/)
33. Trüeb RM. Oxidative stress and its impact on skin, scalp and hair. *Int J Cosmet Sci* 2021. — [PMID 34424547](https://pubmed.ncbi.nlm.nih.gov/34424547/)
34. Jaworsky C, Kligman AM, Murphy GF. Characterization of inflammatory infiltrates in male pattern alopecia. *Br J Dermatol* 1992. — [PMID 1390168](https://pubmed.ncbi.nlm.nih.gov/1390168/)

## 6. The prostaglandin arm — `PGD2`, `SETI`

35. ★ Garza LA, et al. Prostaglandin D2 inhibits hair growth and is elevated in bald scalp of men with androgenetic alopecia. *Sci Transl Med* 2012. — [PMID 22440736](https://pubmed.ncbi.nlm.nih.gov/22440736/) — `PGD_GAIN 0.45`
36. Nieves A, Garza LA. Does prostaglandin D2 hold the cure to male pattern baldness? *Exp Dermatol* 2014. — [PMID 24521203](https://pubmed.ncbi.nlm.nih.gov/24521203/)
37. Khidhir KG, et al. The prostamide-related glaucoma therapy, bimatoprost, offers a novel approach for treating scalp alopecias. *FASEB J* 2013. — [PMID 23104985](https://pubmed.ncbi.nlm.nih.gov/23104985/)
38. Santos Z, Avci P, Hamblin MR. Drug discovery for alopecia: gone today, hair tomorrow. *Expert Opin Drug Discov* 2015. — [PMID 25662177](https://pubmed.ncbi.nlm.nih.gov/25662177/)

## 7. Finasteride — mechanism of action and pharmacokinetics

39. ★ Faller B, Farley D, Nick H. Finasteride: a slow-binding 5α-reductase inhibitor. *Biochemistry* 1993. — [PMID 8389191](https://pubmed.ncbi.nlm.nih.gov/8389191/) — the direct basis of the model's quasi-irreversible trapping (`KDISS 30 d`)
40. ★ Steiner JF. Clinical pharmacokinetics and pharmacodynamics of finasteride. *Clin Pharmacokinet* 1996. — [PMID 8846625](https://pubmed.ncbi.nlm.nih.gov/8846625/) — `FIN_KE`, `FIN_V`, `FIN_F`
41. ★ **Drake L, et al. The effects of finasteride on scalp skin and serum androgen levels in men with androgenetic alopecia. *J Am Acad Dermatol* 1999.** — [PMID 10495374](https://pubmed.ncbi.nlm.nih.gov/10495374/) — **the primary target of the model's dose-response calibration**: scalp DHT −14.9 / −61.6 / −56.5 / −64.1 / −69.4%, serum DHT −49.5 / −68.6 / −71.4 / −72.2% (0.01–5 mg)
42. ★ Roberts JL, et al. Clinical dose ranging studies with finasteride, a type 2 5α-reductase inhibitor, in men with male pattern hair loss. *J Am Acad Dermatol* 1999. — [PMID 10495375](https://pubmed.ncbi.nlm.nih.gov/10495375/)

## 8. Finasteride — clinical trials

43. ★ **Kaufman KD, et al. Finasteride in the treatment of men with androgenetic alopecia. Finasteride Male Pattern Hair Loss Study Group. *J Am Acad Dermatol* 1998.** — [PMID 9777765](https://pubmed.ncbi.nlm.nih.gov/9777765/) — **the primary target for calibrating the model's clinical endpoint**: n=1553, baseline 876 hairs / 5.1 cm², +107 (1 year) · +138 (2 years) versus placebo
44. Leyden J, et al. Finasteride in the treatment of men with frontal male pattern hair loss. *J Am Acad Dermatol* 1999. — [PMID 10365924](https://pubmed.ncbi.nlm.nih.gov/10365924/)
45. ★ Price VH, et al. Changes in hair weight in men with androgenetic alopecia after treatment with finasteride (1 mg daily): three- and 4-year results. *J Am Acad Dermatol* 2006. — [PMID 16781295](https://pubmed.ncbi.nlm.nih.gov/16781295/) — the hair **weight** endpoint = the model's `HMI` (d²-weighted)
46. Kaufman KD, et al. Long-term treatment with finasteride 1 mg decreases the likelihood of developing further visible hair loss. *Eur J Dermatol* 2008. — [PMID 18573712](https://pubmed.ncbi.nlm.nih.gov/18573712/)
47. Olsen EA, et al. Global photographic assessment of men aged 18 to 60 years with male pattern hair loss receiving finasteride 1 mg or placebo. *J Am Acad Dermatol* 2012. — [PMID 22325459](https://pubmed.ncbi.nlm.nih.gov/22325459/)
48. Shapiro J, Kaufman KD. Use of finasteride in the treatment of men with androgenetic alopecia. *J Investig Dermatol Symp Proc* 2003. — [PMID 12894990](https://pubmed.ncbi.nlm.nih.gov/12894990/)
49. Sato A, Takeda A. Evaluation of efficacy and safety of finasteride 1 mg in 3177 Japanese men with androgenetic alopecia. *J Dermatol* 2012. — [PMID 21980923](https://pubmed.ncbi.nlm.nih.gov/21980923/)
50. ★ Price VH, et al. Lack of efficacy of finasteride in postmenopausal women with androgenetic alopecia. *J Am Acad Dermatol* 2000. — [PMID 11050579](https://pubmed.ncbi.nlm.nih.gov/11050579/) — the basis for introducing `ARIND` in the female arm

## 9. Dutasteride

51. ★ Clark RV, et al. Marked suppression of dihydrotestosterone in men with benign prostatic hyperplasia by dutasteride, a dual 5α-reductase inhibitor. *J Clin Endocrinol Metab* 2004. — [PMID 15126539](https://pubmed.ncbi.nlm.nih.gov/15126539/) — serum DHT ~−94%, `KON_DUT1/2`
52. ★ **Gubelin Harcha W, et al. A randomized, active- and placebo-controlled study of the efficacy and safety of different doses of dutasteride versus placebo and finasteride in the treatment of male subjects with androgenetic alopecia. *J Am Acad Dermatol* 2014.** — [PMID 24411083](https://pubmed.ncbi.nlm.nih.gov/24411083/) — at 24 weeks dutasteride 0.5 mg > finasteride 1 mg (p=0.003)

## 10. Minoxidil — mechanism and pharmacokinetics

53. ★ Buhl AE, Waldon DJ, et al. Minoxidil sulfate is the active metabolite that stimulates hair follicles. *J Invest Dermatol* 1990. — [PMID 2230218](https://pubmed.ncbi.nlm.nih.gov/2230218/) — the basis for the model treating `MXSF` as the active species
54. ★ Baker CA, et al. Minoxidil sulfation in the hair follicle. *Skin Pharmacol* 1994. — [PMID 7946376](https://pubmed.ncbi.nlm.nih.gov/7946376/) — follicular SULT1A1
55. ★ Goren A, et al. Clinical utility and validity of minoxidil response testing in androgenetic alopecia. *Dermatol Ther* 2015. — [PMID 25112173](https://pubmed.ncbi.nlm.nih.gov/25112173/) — separation of responders from non-responders (`SULT` 2.0 vs 0.25)
56. Messenger AG, Rundegren J. Minoxidil: mechanisms of action on hair growth. *Br J Dermatol* 2004. — [PMID 14996087](https://pubmed.ncbi.nlm.nih.gov/14996087/)
57. ★ Ferry JJ, et al. Relationship between contact time of applied dose and percutaneous absorption of minoxidil from a topical solution. *J Pharm Sci* 1990. — [PMID 2395092](https://pubmed.ncbi.nlm.nih.gov/2395092/) — systemic absorption of topical minoxidil ~1% (`MX_KSYS`)
58. Suchonwanit P, Thammarucha S, Leerunyakul K. Minoxidil and its use in hair disorders: a review. *Drug Des Devel Ther* 2019. — [PMID 31496654](https://pubmed.ncbi.nlm.nih.gov/31496654/)

## 11. Minoxidil — clinical trials

59. Olsen EA, Weiner MS, et al. Topical minoxidil in early male pattern baldness. *J Am Acad Dermatol* 1985. — [PMID 3900155](https://pubmed.ncbi.nlm.nih.gov/3900155/)
60. ★ **Olsen EA, et al. A randomized clinical trial of 5% topical minoxidil versus 2% topical minoxidil and placebo in the treatment of androgenetic alopecia in men. *J Am Acad Dermatol* 2002.** — [PMID 12196747](https://pubmed.ncbi.nlm.nih.gov/12196747/) — n=393, 48 weeks, 5% > 2% > placebo
61. Blume-Peytavi U, et al. Efficacy and safety of once-daily minoxidil foam 5% versus twice-daily minoxidil solution 2% in female pattern hair loss. *J Drugs Dermatol* 2016. — [PMID 27391640](https://pubmed.ncbi.nlm.nih.gov/27391640/)
62. Randolph M, Tosti A. Oral minoxidil treatment for hair loss: a review of efficacy and safety. *J Am Acad Dermatol* 2021. — [PMID 32622136](https://pubmed.ncbi.nlm.nih.gov/32622136/)
63. ★ Vañó-Galván S, et al. Safety of low-dose oral minoxidil for hair loss: a multicenter study of 1404 patients. *J Am Acad Dermatol* 2021. — [PMID 33639244](https://pubmed.ncbi.nlm.nih.gov/33639244/) — the frequency of hypertrichosis and oedema = the model's `HTR`, `MAP`
64. Jimenez-Cauhe J, et al. Safety of low-dose oral minoxidil treatment for hair loss: a systematic review and pooled-analysis of individual patient data. *Dermatol Ther* 2020. — [PMID 32757405](https://pubmed.ncbi.nlm.nih.gov/32757405/)

## 12. Topical finasteride · combination therapy

65. ★ Piraccini BM, et al. Efficacy and safety of topical finasteride spray solution for male androgenetic alopecia: a phase III, randomized, controlled clinical study. *J Eur Acad Dermatol Venereol* 2022. — [PMID 34634163](https://pubmed.ncbi.nlm.nih.gov/34634163/) — 0.25% topical; suppression of serum DHT is far shallower than with oral
66. Gupta AK, et al. Relative efficacy of minoxidil and the 5-α reductase inhibitors in androgenetic alopecia treatment of male patients: a network meta-analysis. *JAMA Dermatol* 2022. — [PMID 35107565](https://pubmed.ncbi.nlm.nih.gov/35107565/)
67. Adil A, Godwin M. The effectiveness of treatments for androgenetic alopecia: a systematic review and meta-analysis. *J Am Acad Dermatol* 2017. — [PMID 28396101](https://pubmed.ncbi.nlm.nih.gov/28396101/)
68. Rossi A, et al. Multi-therapies in androgenetic alopecia: review and clinical experiences. *Dermatol Ther* 2016. — [PMID 27424565](https://pubmed.ncbi.nlm.nih.gov/27424565/)
69. Piérard-Franchimont C, et al. Ketoconazole shampoo: effect of long-term use in androgenic alopecia. *Dermatology* 1998. — [PMID 9669136](https://pubmed.ncbi.nlm.nih.gov/9669136/)

## 13. Female pattern hair loss

70. Fabbrocini G, et al. Female pattern hair loss: a clinical, pathophysiologic, and therapeutic review. *Int J Womens Dermatol* 2018. — [PMID 30627618](https://pubmed.ncbi.nlm.nih.gov/30627618/)
71. Sinclair R, Wewerinke M, Jolley D. Treatment of female pattern hair loss with oral antiandrogens. *Br J Dermatol* 2005. — [PMID 15787815](https://pubmed.ncbi.nlm.nih.gov/15787815/)
72. Hoedemaker C, Sinclair R. Treatment of female pattern hair loss with a combination of spironolactone and minoxidil. *Australas J Dermatol* 2007. — [PMID 17222303](https://pubmed.ncbi.nlm.nih.gov/17222303/)

## 14. Safety

73. Irwig MS, Kolukula S. Persistent sexual side effects of finasteride for male pattern hair loss. *J Sex Med* 2011. — [PMID 21418145](https://pubmed.ncbi.nlm.nih.gov/21418145/)

---

## Summary of the quantitative targets used for model calibration

| Target | Source | Observed | Model |
|---|---|---|---|
| Baseline hair count in the target area (5.1 cm²) | Kaufman 1998 [9777765] | 876 | 876 (burn-in termination condition) |
| Finasteride 1 mg vs placebo, 1 year | Kaufman 1998 | +107 | **+101.3** |
| Finasteride 1 mg vs placebo, 2 years | Kaufman 1998 | +138 | **+139.2** |
| Scalp DHT, finasteride 0.2 mg, day 42 | Drake 1999 [10495374] | −56.5% | −57.5% |
| Scalp DHT, finasteride 1 mg, day 42 | Drake 1999 | −64.1% | **−63.6%** |
| Scalp DHT, finasteride 5 mg, day 42 | Drake 1999 | −69.4% | −67.4% |
| Serum DHT, finasteride 0.05 mg | Drake 1999 | −49.5% | −47.3% |
| Serum DHT, finasteride 1 mg | Drake 1999 | −71.4% | **−70.3%** |
| Serum DHT, dutasteride 0.5 mg | Clark 2004 [15126539] | ~−94% | −96.1% |
| Dutasteride > finasteride (24 weeks) | Gubelin Harcha 2014 [24411083] | significant (p=0.003) | +80.8 vs +63.4 |
| Rise in serum testosterone (finasteride) | Drake 1999 / label | +9~15% | +13.1% |
| Daily hair shedding on a normal scalp | textbook value | ~100 | 90 (= 100000×0.090/100 d) |
| Normal telogen fraction | textbook value | 9~14% | 9.0% |
| Serum DHT on topical finasteride 0.25% | Piraccini 2022 [34634163] | ~−25% | −27.4% |

### Where the model diverges from the literature (set down rather than hidden)

| Item | Literature | Model | Explanation |
|---|---|---|---|
| **Scalp** DHT on dutasteride | −51 ~ −79% | **−90.3%** | the 5AR-independent floor of the scalp compartment (`W_ALT` 6%) is small. Raising it fits dutasteride but breaks the better-measured finasteride 1 mg point. Reported as an over-prediction. |
| **Scalp** DHT on finasteride 0.05 mg | −61.6% | −43.1% | the raw data are themselves non-monotonic (0.05 mg suppresses more than 0.2 mg). Not reproducible with a mass-action model. The **serum** point at the same dose (−49.5%) fits well at −47.3%. |
| Setipiprant (CRTH2 blockade) | failed clinically | +46 hairs vs placebo | at `PGD_EFF 0.18` a detectable effect is predicted. Read the other way round, for compatibility with the negative trial the contribution of the PGD2 arm must be **5% or less** (`PGD_EFF` scan: 0.18→+46, 0.10→+28, 0.05→+15, 0.02→+6). |
| Rate of decline in the placebo arm | −21 (1 year) / −55 (2 years) | −14.7 / −29.3 | the model's untreated arm is linear, reality accelerates. |
