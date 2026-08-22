#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build tap_references.md from PubMed.

Every title, journal, year, author and PMID in the bibliography is resolved
LIVE from NCBI esearch + esummary.  Nothing is written from memory, because a
plausible-looking PMID written from memory is a fabricated citation.

Each entry is a (section, intent, query) triple.  The top RELEVANCE-ranked hit
is taken and printed together with its INTENT, so a retrieved paper that does
not match what it was asked for is visible at a glance instead of hidden.

E-utilities defaults to DATE order, not relevance, so sort=relevance is
mandatory: without it the "top hit" is simply the newest paper matching any
term.

    python3 mkrefs.py            # build (uses refs_raw.json cache if present)
    python3 mkrefs.py --refresh   # ignore the cache and re-query
"""
import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import OrderedDict

B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CACHE = "refs_raw.json"
OUT = "tap_references.md"


def _get(url):
    for k in range(6):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode(
                "utf-8", "replace")
        except Exception:
            time.sleep(2.5 + 3 * k)
    return ""


def esearch(term, n=1):
    u = B + ("esearch.fcgi?db=pubmed&retmode=json&sort=relevance&retmax=%d"
             "&term=%s" % (n, urllib.parse.quote(term)))
    try:
        return json.loads(_get(u))["esearchresult"]["idlist"]
    except Exception:
        return []


def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 120):
        chunk = pmids[i:i + 120]
        u = B + "esummary.fcgi?db=pubmed&retmode=json&id=" + ",".join(chunk)
        try:
            r = json.loads(_get(u))["result"]
        except Exception:
            continue
        for pid in chunk:
            if pid in r:
                out[pid] = r[pid]
        time.sleep(0.34)
    return out


# ---------------------------------------------------------------------------
#  (section, intent, query)
#  The INTENT column is the point: it records what the citation is being asked
#  to support, so that a mismatch is a visible defect rather than a silent one.
# ---------------------------------------------------------------------------
Q = [
 # ---------------------------------------------------------------- 1
 ("1 · Overview · epidemiology · guidelines",
  "modern review of methanol poisoning",
  "methanol poisoning review clinical"),
 ("1 · Overview · epidemiology · guidelines",
  "modern review of ethylene glycol poisoning",
  "ethylene glycol poisoning review diagnosis treatment"),
 ("1 · Overview · epidemiology · guidelines",
  "AACT/EAPCCT practice guideline on methanol",
  "American Academy Clinical Toxicology practice guidelines methanol poisoning"),
 ("1 · Overview · epidemiology · guidelines",
  "AACT practice guideline on ethylene glycol",
  "American Academy Clinical Toxicology practice guidelines ethylene glycol poisoning"),
 ("1 · Overview · epidemiology · guidelines",
  "toxic alcohols in the intensive care unit",
  "toxic alcohol poisoning intensive care management"),
 ("1 · Overview · epidemiology · guidelines",
  "Estonia 2001 mass methanol outbreak",
  "methanol mass poisoning outbreak Estonia 2001"),
 ("1 · Overview · epidemiology · guidelines",
  "Iran 2020 methanol outbreak during COVID-19",
  "methanol poisoning outbreak Iran COVID-19 hand sanitizer"),
 ("1 · Overview · epidemiology · guidelines",
  "Norway 2002-2004 methanol outbreak cohort",
  "methanol poisoning outbreak Norway prognostic factors cohort"),
 ("1 · Overview · epidemiology · guidelines",
  "Czech 2012 mass methanol poisoning outbreak",
  "Czech Republic mass methanol poisoning 2012 outbreak"),
 ("1 · Overview · epidemiology · guidelines",
  "global burden and prevention of methanol outbreaks",
  "methanol poisoning prevention strategies WHO"),

 # ---------------------------------------------------------------- 2
 ("2 · Alcohol dehydrogenase kinetics",
  "human ADH isoenzyme kinetics and substrate specificity",
  "human alcohol dehydrogenase isoenzyme kinetic properties substrate"),
 ("2 · Alcohol dehydrogenase kinetics",
  "kinetics of methanol oxidation by human ADH",
  "alcohol dehydrogenase methanol oxidation kinetics Km"),
 ("2 · Alcohol dehydrogenase kinetics",
  "4-methylpyrazole inhibition constant for alcohol dehydrogenase",
  "4-methylpyrazole alcohol dehydrogenase inhibition constant Ki"),
 ("2 · Alcohol dehydrogenase kinetics",
  "ethanol elimination kinetics Michaelis-Menten in man",
  "ethanol elimination kinetics Michaelis-Menten human"),
 ("2 · Alcohol dehydrogenase kinetics",
  "ADH1B/ADH1C polymorphism and alcohol metabolism rate",
  "ADH1B ADH1C polymorphism ethanol metabolism rate"),
 ("2 · Alcohol dehydrogenase kinetics",
  "catalase and CYP2E1 contribution to methanol oxidation",
  "methanol oxidation catalase"),
 ("2 · Alcohol dehydrogenase kinetics",
  "ethylene glycol oxidation by alcohol dehydrogenase",
  "ethylene glycol alcohol dehydrogenase oxidation glycolaldehyde"),

 # ---------------------------------------------------------------- 3
 ("3 · Methanol · formate toxicity mechanism",
  "formate is the metabolite responsible for methanol toxicity",
  "formate accumulation methanol toxicity mechanism"),
 ("3 · Methanol · formate toxicity mechanism",
  "formate inhibits cytochrome c oxidase",
  "formate inhibition cytochrome c oxidase mitochondrial"),
 ("3 · Methanol · formate toxicity mechanism",
  "folate dependence of formate oxidation and species differences",
  "folate dependent formate oxidation methanol species difference monkey"),
 ("3 · Methanol · formate toxicity mechanism",
  "blood formate concentrations in human methanol poisoning",
  "blood formate concentration human methanol poisoning outcome"),
 ("3 · Methanol · formate toxicity mechanism",
  "formate kinetics and elimination half-life in man",
  "formate elimination kinetics half-life human methanol"),
 ("3 · Methanol · formate toxicity mechanism",
  "methanol pharmacokinetics and elimination rate in poisoning",
  "methanol pharmacokinetics elimination rate poisoned patients"),
 ("3 · Methanol · formate toxicity mechanism",
  "methanol half-life during ADH blockade with fomepizole",
  "methanol elimination half-life fomepizole"),
 ("3 · Methanol · formate toxicity mechanism",
  "folinic acid / folate therapy for methanol poisoning",
  "folinic acid leucovorin methanol poisoning formate"),

 # ---------------------------------------------------------------- 4
 ("4 · Optic nerve · basal ganglia injury",
  "visual toxicity and retinal mechanism in methanol poisoning",
  "methanol optic neuropathy retina mechanism formate"),
 ("4 · Optic nerve · basal ganglia injury",
  "bilateral putaminal necrosis on imaging after methanol",
  "methanol poisoning putaminal necrosis MRI basal ganglia"),
 ("4 · Optic nerve · basal ganglia injury",
  "long-term visual and neurological sequelae after methanol",
  "methanol poisoning long term visual neurological sequelae follow-up"),
 ("4 · Optic nerve · basal ganglia injury",
  "delayed parkinsonism after methanol poisoning",
  "methanol poisoning parkinsonism"),
 ("4 · Optic nerve · basal ganglia injury",
  "vitreous and CSF formate compared with plasma",
  "vitreous humour formate methanol postmortem concentration"),
 ("4 · Optic nerve · basal ganglia injury",
  "optical coherence tomography of the optic nerve after methanol",
  "methanol optic neuropathy optical coherence tomography retinal nerve fiber"),

 # ---------------------------------------------------------------- 5
 ("5 · Ethylene glycol metabolism · oxalate",
  "glycolate is the principal acid in ethylene glycol poisoning",
  "glycolic acid concentration ethylene glycol intoxication"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "glycolate correlates with outcome better than ethylene glycol",
  "glycolate concentration correlation acidosis outcome ethylene glycol"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "calcium oxalate crystal nephropathy from ethylene glycol",
  "calcium oxalate crystal deposition renal ethylene glycol poisoning"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "glycolaldehyde and glyoxylate direct tubular cytotoxicity",
  "glycolaldehyde glyoxylate cytotoxicity renal tubular ethylene glycol"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "thiamine and pyridoxine as adjuncts in ethylene glycol poisoning",
  "pyridoxine thiamine ethylene glycol poisoning"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "ethylene glycol pharmacokinetics with and without ADH blockade",
  "ethylene glycol pharmacokinetics half-life fomepizole renal excretion"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "hypocalcaemia and QT prolongation in ethylene glycol poisoning",
  "hypocalcemia ethylene glycol poisoning"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "renal outcome and recovery after ethylene glycol AKI",
  "ethylene glycol poisoning acute kidney injury renal recovery outcome"),
 ("5 · Ethylene glycol metabolism · oxalate",
  "solubility product of calcium oxalate and urinary supersaturation",
  "calcium oxalate supersaturation solubility product urine crystallization"),

 # ---------------------------------------------------------------- 6
 ("6 · Fomepizole",
  "fomepizole efficacy and safety in methanol and ethylene glycol poisoning",
  "fomepizole efficacy safety methanol ethylene glycol poisoning trial"),
 ("6 · Fomepizole",
  "fomepizole pharmacokinetics and dosing in man",
  "fomepizole pharmacokinetics dosing healthy volunteers"),
 ("6 · Fomepizole",
  "fomepizole removal by haemodialysis and dose adjustment",
  "fomepizole hemodialysis clearance dosing during dialysis"),
 ("6 · Fomepizole",
  "fomepizole without dialysis for ethylene glycol poisoning",
  "fomepizole alone without hemodialysis ethylene glycol poisoning"),
 ("6 · Fomepizole",
  "fomepizole in paediatric toxic alcohol poisoning",
  "fomepizole pediatric poisoning case series children"),
 ("6 · Fomepizole",
  "fomepizole versus ethanol comparative effectiveness",
  "fomepizole ethanol comparison antidote"),
 ("6 · Fomepizole",
  "fomepizole shortage, cost and rationing in outbreaks",
  "fomepizole cost availability"),

 # ---------------------------------------------------------------- 7
 ("7 · Ethanol as antidote",
  "ethanol as an antidote: dosing and target concentration",
  "ethanol antidote methanol ethylene glycol dosing target concentration"),
 ("7 · Ethanol as antidote",
  "ethanol infusion during haemodialysis requires higher doses",
  "ethanol infusion hemodialysis dose requirement methanol antidote"),
 ("7 · Ethanol as antidote",
  "adverse effects of therapeutic ethanol infusion",
  "ethanol antidote adverse effects hypoglycemia"),
 ("7 · Ethanol as antidote",
  "co-ingested ethanol delays and mitigates methanol toxicity",
  "ethanol coingestion methanol poisoning"),

 # ---------------------------------------------------------------- 8
 ("8 · Extracorporeal removal",
  "EXTRIP recommendations for methanol",
  "EXTRIP extracorporeal treatment methanol poisoning recommendations"),
 ("8 · Extracorporeal removal",
  "EXTRIP recommendations for ethylene glycol",
  "EXTRIP extracorporeal treatment ethylene glycol poisoning workgroup"),
 ("8 · Extracorporeal removal",
  "haemodialysis clearance of methanol and formate",
  "formate clearance dialysis methanol kinetics"),
 ("8 · Extracorporeal removal",
  "continuous renal replacement therapy for toxic alcohols",
  "continuous venovenous hemodiafiltration ethylene glycol elimination"),
 ("8 · Extracorporeal removal",
  "rebound of toxic alcohol concentration after dialysis",
  "rebound hemodialysis ethylene glycol concentration"),
 ("8 · Extracorporeal removal",
  "when to stop dialysis in toxic alcohol poisoning",
  "duration hemodialysis methanol poisoning"),

 # ---------------------------------------------------------------- 9
 ("9 · Osmolal gap · anion gap",
  "osmolal gap in the diagnosis of toxic alcohol ingestion",
  "osmolal gap toxic alcohol ingestion diagnosis utility"),
 ("9 · Osmolal gap · anion gap",
  "the normal range of the osmolal gap is wide",
  "osmolal gap normal range"),
 ("9 · Osmolal gap · anion gap",
  "vapour-pressure versus freezing-point osmometry misses volatile alcohols",
  "vapor pressure osmometry freezing point depression volatile alcohol measurement"),
 ("9 · Osmolal gap · anion gap",
  "reciprocal relationship of osmolal and anion gaps over time",
  "osmolal gap decrease anion gap increase methanol"),
 ("9 · Osmolal gap · anion gap",
  "conversion factors from alcohol concentration to osmolality",
  "osmolal gap calculation conversion factor ethanol methanol ethylene glycol"),
 ("9 · Osmolal gap · anion gap",
  "high anion gap metabolic acidosis differential diagnosis",
  "high anion gap metabolic acidosis differential diagnosis approach"),
 ("9 · Osmolal gap · anion gap",
  "isopropanol poisoning: osmolal gap without acidosis",
  "isopropyl alcohol poisoning ketosis osmolal gap without acidosis"),
 ("9 · Osmolal gap · anion gap",
  "propylene glycol accumulation from sedative infusions",
  "propylene glycol toxicity lorazepam infusion lactic acidosis osmolal gap"),

 # ---------------------------------------------------------------- 10
 ("10 · Acid-base · ion trapping",
  "apparent bicarbonate space and its expansion in acidosis",
  "apparent bicarbonate space distribution volume metabolic acidosis"),
 ("10 · Acid-base · ion trapping",
  "sodium bicarbonate therapy in severe metabolic acidosis",
  "sodium bicarbonate therapy severe metabolic acidosis outcome"),
 ("10 · Acid-base · ion trapping",
  "urinary alkalinisation and non-ionic diffusion of weak acids",
  "urinary alkalinization poisoning elimination"),
 ("10 · Acid-base · ion trapping",
  "pH-dependent distribution of formate into tissue",
  "formate distribution pH methanol toxicity"),
 ("10 · Acid-base · ion trapping",
  "Winter's formula and expected respiratory compensation",
  "expected respiratory compensation metabolic acidosis Winter formula"),
 ("10 · Acid-base · ion trapping",
  "delta anion gap to delta bicarbonate ratio in organic acidosis",
  "delta anion gap delta bicarbonate ratio"),

 # ---------------------------------------------------------------- 11
 ("11 · Prognostic factors",
  "arterial pH predicts outcome in methanol poisoning",
  "arterial pH prognostic factor mortality methanol poisoning"),
 ("11 · Prognostic factors",
  "predictors of visual loss after methanol poisoning",
  "predictors visual outcome blindness methanol poisoning"),
 ("11 · Prognostic factors",
  "coma and seizure on presentation as predictors of death",
  "coma seizures presentation predictor mortality methanol poisoning"),
 ("11 · Prognostic factors",
  "methanol concentration is a poor predictor of severity",
  "serum methanol concentration outcome predictor"),
 ("11 · Prognostic factors",
  "prognostic scoring in ethylene glycol poisoning",
  "ethylene glycol poisoning mortality prognostic"),

 # ---------------------------------------------------------------- 12
 ("12 · Diagnosis · laboratory pitfalls",
  "glycolate interference with lactate assays: the lactate gap",
  "glycolate interference lactate assay lactate gap ethylene glycol"),
 ("12 · Diagnosis · laboratory pitfalls",
  "spurious creatinine elevation from glycolate",
  "glycolate interference creatinine assay falsely elevated ethylene glycol"),
 ("12 · Diagnosis · laboratory pitfalls",
  "gas chromatographic measurement of toxic alcohols",
  "gas chromatography measurement methanol ethylene glycol serum method"),
 ("12 · Diagnosis · laboratory pitfalls",
  "urinary fluorescence and calcium oxalate crystals are unreliable",
  "urine fluorescence ethylene glycol"),
 ("12 · Diagnosis · laboratory pitfalls",
  "point-of-care and enzymatic assays for ethylene glycol",
  "enzymatic assay ethylene glycol serum"),
 ("12 · Diagnosis · laboratory pitfalls",
  "diethylene glycol mass poisoning and its distinct metabolite",
  "diethylene glycol poisoning outbreak"),

 # ---------------------------------------------------------------- 13
 ("13 · QSP · PBPK · modelling",
  "physiologically based pharmacokinetic model of methanol",
  "physiologically based pharmacokinetic model methanol formate human"),
 ("13 · QSP · PBPK · modelling",
  "kinetic model of ethylene glycol and its metabolites",
  "pharmacokinetic model ethylene glycol glycolate simulation"),
 ("13 · QSP · PBPK · modelling",
  "simulation of antidote and dialysis strategies in toxic alcohol poisoning",
  "mathematical model methanol kinetics dialysis prediction"),
 ("13 · QSP · PBPK · modelling",
  "mrgsolve for ODE-based pharmacometric simulation",
  "mrgsolve R package simulation pharmacometrics"),
 ("13 · QSP · PBPK · modelling",
  "quantitative systems pharmacology model credibility and validation",
  "quantitative systems pharmacology model validation"),
 ("13 · QSP · PBPK · modelling",
  "toxicokinetic modelling to guide antidote dosing",
  "toxicokinetic modelling antidote dosing decision support poisoning"),
]


def author_str(a):
    if not a:
        return "—"
    names = [x.get("name", "") for x in a]
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return "%s, %s" % (names[0], names[1])
    return "%s et al." % names[0]


def year_of(rec):
    m = re.match(r"(\d{4})", rec.get("pubdate", "") or "")
    return m.group(1) if m else "n.d."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args()

    cache = {}
    if os.path.exists(CACHE) and not a.refresh:
        cache = json.load(open(CACHE))

    resolved = OrderedDict()
    for sec, intent, q in Q:
        key = q
        if key not in cache:
            ids = esearch(q, 1)
            if ids:
                cache[key] = ids
            time.sleep(0.45)
        ids = cache.get(key, [])
        if ids:
            resolved[(sec, intent, q)] = ids[0]
        else:
            resolved[(sec, intent, q)] = None
    json.dump(cache, open(CACHE, "w"), indent=1)

    pmids = sorted({v for v in resolved.values() if v})
    meta_file = "refs_meta.json"
    meta = json.load(open(meta_file)) if (os.path.exists(meta_file) and not a.refresh) else {}
    need = [p for p in pmids if p not in meta]
    if need:
        meta.update(esummary(need))
        json.dump(meta, open(meta_file, "w"), indent=1)

    # ---- assemble, de-duplicating within a section -----------------------
    secs = OrderedDict()
    seen = set()
    dupes = []
    for (sec, intent, q), pid in resolved.items():
        secs.setdefault(sec, [])
        if pid is None:
            secs[sec].append((intent, None, None))
            continue
        if pid in seen:
            dupes.append((sec, intent, pid))
        seen.add(pid)
        secs[sec].append((intent, pid, meta.get(pid, {})))

    n = sum(1 for s in secs.values() for e in s if e[1])
    lines = []
    lines.append("# Toxic Alcohol Poisoning (methanol · ethylene glycol) — References")
    lines.append("### Toxic Alcohol Poisoning · References")
    lines.append("")
    lines.append("**%d references in total** · 13 sections · every PMID·title·journal·year "
                 "was **queried live** from NCBI E-utilities by `mkrefs.py`, and none "
                 "of it was written from memory." % n)
    lines.append("")
    lines.append("The `intent` of each entry records **what that citation was asked "
                 "to support**. If a search result does not match the intent, this is so "
                 "that the mismatch is visible at a glance instead of being hidden.")
    lines.append("")
    lines.append("Regenerate:")
    lines.append("```bash")
    lines.append("python3 mkrefs.py --refresh    # re-query PubMed")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")
    for sec, entries in secs.items():
        lines.append("## %s" % sec)
        lines.append("")
        for intent, pid, rec in entries:
            if not pid:
                lines.append("- *(unresolved query — no PubMed hit)* — intent: %s" % intent)
                continue
            title = (rec.get("title") or "").strip().rstrip(".")
            src = rec.get("source", "")
            lines.append("- **%s** (%s, %s). %s  " %
                         (author_str(rec.get("authors")), src, year_of(rec), title))
            lines.append("  [PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/) · "
                         "intent: *%s*" % (pid, pid, intent))
        lines.append("")
    if dupes:
        lines.append("---")
        lines.append("")
        lines.append("## Duplicate query record (queries that resolved to an already-cited paper)")
        lines.append("")
        lines.append("These are cases where searches with different intents converged on "
                     "the same paper. They are left in rather than hidden.")
        lines.append("")
        for sec, intent, pid in dupes:
            lines.append("- [PMID %s](https://pubmed.ncbi.nlm.nih.gov/%s/) — %s / *%s*"
                         % (pid, pid, sec.split(" · ")[0], intent))
        lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Where these are used")
    lines.append("")
    lines.append("| Model component | Evidence section |")
    lines.append("|---|---|")
    lines.append("| ADH competitive inhibition (`KM_M`·`KM_E`·`KM_ET`·`KI_FOM`) | 2 |")
    lines.append("| Formate production · folate-dependent oxidation (`VMAX_THF`·`THF`) | 3 |")
    lines.append("| Ion trapping · CNS delivery (`PSF`·`fHA`) | 3, 4, 10 |")
    lines.append("| Optic nerve · basal ganglia damage threshold (`THR_OPTA`·`THR_PUT`) | 4 |")
    lines.append("| Glycolate · glyoxylate branching (`BR_OX`·`BR_TH`·`BR_PY`) | 5 |")
    lines.append("| CaOx supersaturation · Ksp · renal injury (`KSP`·`FP_MAX`) | 5 |")
    lines.append("| Fomepizole PK · dose during dialysis (`VMAX_FOM`·`CL_HD_FOM`) | 6 |")
    lines.append("| Ethanol antidote target concentration (`KM_ET`·`ETHHDMLT`) | 7 |")
    lines.append("| Dialysis clearance · rebound · stopping criteria (`CL_HD`·`CL_CRRT`) | 8 |")
    lines.append("| The two-gap clock · conversion factors (`OG`·`AG`·`GAPSUM`) | 9 |")
    lines.append("| Bicarbonate space · respiratory compensation (`bicspace`·Winter) | 10 |")
    lines.append("| Death · blindness hazard functions (`H_PH`·`H_BLIND`) | 11 |")
    lines.append("| Lactate gap · creatinine interference (`LACTGAP`) | 12 |")
    lines.append("| Model structure · verification methodology | 13 |")
    lines.append("")
    open(OUT, "w").write("\n".join(lines) + "\n")
    print("wrote %s: %d references, %d sections, %d duplicate resolutions"
          % (OUT, n, len(secs), len(dupes)))
    miss = [(s, i) for s, e in secs.items() for i, p, _ in e if not p]
    if miss:
        print("UNRESOLVED:")
        for s, i in miss:
            print("  ", s, "|", i)


if __name__ == "__main__":
    main()
