#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""abm_references.md generator.

Only the records actually returned by the NCBI E-utilities (esearch + esummary) are copied across.
Not one citation is written by hand — PMID·author·journal·year are all API response values.
"""
import json
import time
import urllib.request
import urllib.parse
import sys

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

SECTIONS = [
    ("Epidemiology · clinical presentation · prognosis", "Epidemiology, clinical presentation, prognosis", [
        "community-acquired bacterial meningitis adults prospective cohort",
        "pneumococcal meningitis adults outcome mortality",
        "bacterial meningitis clinical features diagnostic accuracy triad",
        "bacterial meningitis global burden epidemiology",
    ]),
    ("Pneumococcal virulence factors — the cargo that will be released", "Pneumococcal virulence factors", [
        "pneumolysin neurotoxicity meningitis",
        "pneumococcal cell wall inflammation cerebrospinal fluid",
        "LytA autolysin pneumococcus lysis",
        "pneumococcal capsule serotype invasive disease virulence",
        "pneumococcal hydrogen peroxide SpxB neuronal damage",
    ]),
    ("Invasion route and crossing of the blood-CSF barrier", "Invasion and crossing of the blood-CSF barrier", [
        "Streptococcus pneumoniae invasion blood brain barrier endothelial",
        "choroid plexus epithelium bacterial invasion meningitis",
        "pneumococcal nasopharyngeal colonization invasion bacteremia",
        "platelet activating factor receptor pneumococcal adherence",
    ]),
    ("Innate immune recognition — TLR · NOD · inflammasome", "Innate recognition", [
        "TLR2 pneumococcal meningitis innate immunity",
        "NLRP3 inflammasome pneumococcal meningitis interleukin-1",
        "NOD2 peptidoglycan central nervous system inflammation",
        "complement pneumococcal meningitis C5a",
        "microglia activation bacterial meningitis",
    ]),
    ("CSF cytokine network", "CSF cytokine network", [
        "cerebrospinal fluid tumor necrosis factor bacterial meningitis",
        "interleukin-1 beta cerebrospinal fluid meningitis outcome",
        "interleukin-6 interleukin-8 cerebrospinal fluid bacterial meningitis",
        "interleukin-10 cerebrospinal fluid bacterial meningitis",
    ]),
    ("Neutrophil influx and matrix metalloproteinases", "Neutrophil influx and matrix metalloproteinases", [
        "matrix metalloproteinase-9 cerebrospinal fluid bacterial meningitis",
        "neutrophil recruitment pleocytosis experimental meningitis",
        "MMP inhibitor experimental pneumococcal meningitis neuroprotection",
        "neutrophil extracellular traps cerebrospinal fluid meningitis",
    ]),
    ("Barrier permeability and the albumin quotient", "Barrier permeability and the albumin quotient", [
        "blood-cerebrospinal fluid barrier albumin quotient meningitis",
        "blood brain barrier permeability experimental bacterial meningitis",
        "tight junction claudin occludin bacterial meningitis",
        "VEGF cerebrospinal fluid bacterial meningitis barrier",
    ]),
    ("CSF hydrodynamics · intracranial pressure · compliance", "CSF hydrodynamics, intracranial pressure, compliance", [
        "intracranial pressure bacterial meningitis monitoring adults",
        "cerebrospinal fluid outflow resistance meningitis",
        "pressure volume index intracranial compliance",
        "brain edema bacterial meningitis mechanism",
        "hydrocephalus bacterial meningitis",
    ]),
    ("Cerebral perfusion · autoregulation · ischaemia", "Cerebral perfusion, autoregulation, ischaemia", [
        "cerebral blood flow autoregulation bacterial meningitis",
        "cerebral perfusion pressure meningitis outcome",
        "cerebral infarction pneumococcal meningitis",
        "transcranial Doppler bacterial meningitis vasculopathy",
    ]),
    ("CSF metabolism — glucose and lactate", "CSF metabolism: glucose and lactate", [
        "cerebrospinal fluid lactate bacterial meningitis diagnostic",
        "cerebrospinal fluid glucose ratio bacterial meningitis",
        "glucose transport blood brain barrier GLUT1 kinetics",
        "cerebrospinal fluid glucose consumption leukocytes",
    ]),
    ("Antibiotic-induced lysis and the inflammatory burst", "Antibiotic-induced lysis and the inflammatory burst", [
        "antibiotic induced release pneumococcal cell wall inflammation meningitis",
        "nonbacteriolytic antibiotic experimental pneumococcal meningitis",
        "cerebrospinal fluid cytokine increase after antibiotic meningitis",
        "rifampin versus ceftriaxone experimental meningitis inflammation",
        "daptomycin experimental pneumococcal meningitis non-lytic",
    ]),
    ("Antibiotic CSF penetration and PK", "Antibiotic CSF penetration and pharmacokinetics", [
        "ceftriaxone cerebrospinal fluid penetration meningitis pharmacokinetics",
        "vancomycin cerebrospinal fluid penetration meningitis",
        "rifampin cerebrospinal fluid penetration",
        "antibiotic pharmacokinetics cerebrospinal fluid penetration review",
        "ceftriaxone protein binding saturable pharmacokinetics",
        "dexamethasone effect vancomycin cerebrospinal fluid concentration",
    ]),
    ("Antibiotic PD — experimental meningitis models", "Antibiotic PD in experimental meningitis", [
        "rabbit model pneumococcal meningitis bactericidal activity cerebrospinal fluid",
        "pharmacodynamics antibiotic cerebrospinal fluid bactericidal rate meningitis",
        "time above MIC cerebrospinal fluid beta-lactam meningitis",
        "continuous infusion ceftriaxone meningitis",
    ]),
    ("Resistant pneumococci and combination therapy", "Resistant pneumococci and combination therapy", [
        "penicillin resistant pneumococcal meningitis treatment failure",
        "cephalosporin resistant Streptococcus pneumoniae meningitis vancomycin",
        "vancomycin plus ceftriaxone synergy pneumococcal meningitis",
        "PBP2x mosaic gene cephalosporin resistance pneumococcus",
    ]),
    ("Dexamethasone — trials and meta-analyses", "Dexamethasone: trials and meta-analyses", [
        "dexamethasone adults bacterial meningitis randomized trial",
        "corticosteroids acute bacterial meningitis Cochrane meta-analysis",
        "dexamethasone Haemophilus influenzae meningitis children hearing",
        "dexamethasone pneumococcal meningitis mortality subgroup",
        "adjunctive dexamethasone meningitis low income countries trial",
    ]),
    ("Dexamethasone — mechanism and timing of administration", "Dexamethasone: mechanism and timing", [
        "dexamethasone timing before antibiotic experimental meningitis",
        "glucocorticoid NF-kB transrepression mechanism",
        "dexamethasone cerebrospinal fluid inflammation experimental meningitis",
        "dexamethasone hippocampal apoptosis experimental meningitis",
    ]),
    ("Osmotherapy · intracranial pressure management", "Osmotherapy and ICP management", [
        "glycerol adjuvant therapy bacterial meningitis trial",
        "mannitol intracranial pressure treatment",
        "hypertonic saline intracranial pressure",
        "intracranial pressure targeted treatment bacterial meningitis outcome",
    ]),
    ("Hearing loss and cochlear injury", "Hearing loss and cochlear injury", [
        "hearing loss pneumococcal meningitis children",
        "labyrinthitis ossificans meningitis cochlear implantation",
        "cochlear injury experimental pneumococcal meningitis",
        "sensorineural hearing loss bacterial meningitis adults",
    ]),
    ("Cognitive sequelae and the hippocampal dentate gyrus", "Cognitive sequelae and the dentate gyrus", [
        "hippocampal apoptosis dentate gyrus experimental pneumococcal meningitis",
        "cognitive impairment survivors bacterial meningitis",
        "learning deficit experimental meningitis infant rat",
        "neuronal injury markers cerebrospinal fluid meningitis",
    ]),
    ("Delay to antibiotic and outcome", "Time to antibiotic and outcome", [
        "door to antibiotic time bacterial meningitis outcome",
        "delay antibiotic therapy adverse outcome bacterial meningitis",
        "prehospital antibiotic meningitis mortality",
    ]),
    ("Vaccination and serotype replacement", "Vaccination and serotype replacement", [
        "pneumococcal conjugate vaccine meningitis incidence",
        "serotype replacement pneumococcal conjugate vaccine",
        "pneumococcal vaccination adults invasive disease",
    ]),
    ("Diagnosis — CSF studies and molecular tests", "Diagnosis: CSF studies and molecular tests", [
        "multiplex PCR cerebrospinal fluid meningitis diagnosis",
        "cerebrospinal fluid analysis bacterial versus viral meningitis score",
        "cranial CT before lumbar puncture meningitis",
        "procalcitonin bacterial meningitis diagnosis",
    ]),
    ("Guidelines and reviews", "Guidelines and reviews", [
        "ESCMID guideline acute bacterial meningitis",
        "IDSA practice guidelines bacterial meningitis management",
        "acute bacterial meningitis review Lancet",
        "community-acquired bacterial meningitis review New England",
    ]),
    ("Quantitative modelling — QSP · PK/PD · systems biology", "Quantitative modelling: QSP, PK/PD, systems biology", [
        "mathematical model bacterial meningitis dynamics",
        "quantitative systems pharmacology infectious disease model",
        "pharmacokinetic pharmacodynamic model antibacterial killing bacteria",
        "mrgsolve simulation pharmacometrics",
        "physiologically based model central nervous system drug distribution",
        "intracranial pressure mathematical model cerebrospinal fluid dynamics",
    ]),
]

RETMAX = 8


def eget(path, params):
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode())
        except Exception as exc:                     # noqa: BLE001
            print(f"    retry {attempt+1}: {exc}", file=sys.stderr)
            time.sleep(2 + 2 * attempt)
    return None


def search(term):
    js = eget("esearch.fcgi", dict(db="pubmed", term=term, retmax=RETMAX,
                                   retmode="json", sort="relevance"))
    if not js:
        return []
    return js.get("esearchresult", {}).get("idlist", [])


def summary(pmids):
    if not pmids:
        return {}
    js = eget("esummary.fcgi", dict(db="pubmed", id=",".join(pmids), retmode="json"))
    return (js or {}).get("result", {})


def main():
    seen = set()
    out = []
    total = 0
    for ko, en, queries in SECTIONS:
        recs = []
        for q in queries:
            ids = search(q)
            time.sleep(0.4)
            res = summary(ids)
            time.sleep(0.4)
            for pid in ids:
                if pid in seen or pid not in res:
                    continue
                r = res[pid]
                title = (r.get("title") or "").strip().rstrip(".")
                if not title:
                    continue
                auth = r.get("authors") or []
                first = auth[0].get("name", "") if auth else ""
                al = f"{first} et al." if len(auth) > 1 else first
                journal = r.get("source", "")
                year = (r.get("pubdate", "") or "")[:4]
                recs.append((pid, al, title, journal, year))
                seen.add(pid)
            print(f"  {ko} / {q}: {len(ids)} hits", file=sys.stderr)
        recs.sort(key=lambda x: (x[4] or "0"), reverse=True)
        out.append((ko, en, recs))
        total += len(recs)

    L = []
    W = L.append
    W("# Acute Bacterial Meningitis (pneumococcus) — References")
    W("")
    W("**Acute Bacterial Meningitis (*Streptococcus pneumoniae*) — References**")
    W("")
    W(f"**{total} papers** in total.  Every entry was obtained by querying PubMed directly with "
      "the NCBI E-utilities (`esearch` + `esummary`) and copying **only the records returned** — "
      "PMID·author·journal·year are all API response values and no citation was written by hand.  The generating script is this directory's "
      "`fetch_refs.py`.")
    W("")
    W("| Section | Papers |")
    W("|------|------|")
    for i, (ko, en, recs) in enumerate(out, 1):
        anchor = f"{i}-{ko.replace(' ', '-').replace('·', '')}"
        W(f"| {i}. {ko} | {len(recs)} |")
    W("")
    W("---")
    W("")
    for i, (ko, en, recs) in enumerate(out, 1):
        W(f"## {i}. {ko}")
        W("")
        W(f"*{en}*")
        W("")
        for (pid, al, title, journal, year) in recs:
            W(f"- {al} **{title}.** *{journal}* {year}. "
              f"[PMID {pid}](https://pubmed.ncbi.nlm.nih.gov/{pid}/)")
        W("")
    W("---")
    W("")
    W("## Which literature each model parameter hangs on")
    W("")
    W("| Model element | Value | Section carrying the evidence |")
    W("|-----------|-----|------------------|")
    rows = [
        ("CSF bacterial density (at presentation)", "10⁵–10⁸ CFU/mL", "1, 13"),
        ("Pneumococcal CSF growth rate μ_max", "0.85 /h (doubling 49 min)", "13"),
        ("Ceftriaxone CSF penetration", "total conc 2–10 %", "12"),
        ("Ceftriaxone protein-binding saturation", "f_u 0.076 → 0.167", "12"),
        ("Vancomycin CSF penetration", "total conc 2–13 %", "12"),
        ("Rifampicin CSF penetration", "10–20 %, independent of inflammation", "12"),
        ("β-lactam bactericidal rate", "−0.5 ~ −0.7 log₁₀ CFU/mL/h", "13"),
        ("Vancomycin bactericidal rate", "−0.3 ~ −0.4 log₁₀ CFU/mL/h", "13"),
        ("Lytic versus non-lytic yield", "Y_rif ≈ 0.15 × Y_βlactam", "11"),
        ("CSF TNF-α", "100–1,000 pg/mL", "5"),
        ("CSF IL-6", "10⁴–10⁶ pg/mL", "5"),
        ("CSF MMP-9", "100–1,000 ng/mL", "6"),
        ("CSF white cells (pneumococcus)", "1,000–5,000 /µL", "1, 22"),
        ("CSF protein", "100–500 mg/dL", "1, 22"),
        ("CSF glucose ratio", "<0.4 (about 70 %)", "10, 22"),
        ("CSF lactate", ">3.5 mmol/L (commonly 6–12)", "10"),
        ("Albumin quotient Q_alb", "normal <8 → meningitis 30–100+", "7"),
        ("CSF formation rate Q_f", "0.35 mL/min (500 mL/day)", "8"),
        ("CSF outflow resistance R_out", "normal 6–10 mmHg/(mL/min)", "8"),
        ("Pressure-volume index PVI", "≈25 mL", "8"),
        ("Intracranial pressure", "normal 7–15 → 25–40 mmHg", "8"),
        ("Ischaemic threshold CBF", "<25 mL/100 g/min", "9"),
        ("Dexamethasone regimen", "0.15 mg/kg q6h × 4 days", "15"),
        ("Dexamethasone effect on mortality (pneumococcus)", "34 % → 14 %", "15"),
        ("Dexamethasone poor outcome", "52 % → 26 %", "15"),
        ("Severe deafness (paediatric Hib)", "~15 % → ~5 %", "15, 18"),
        ("Deafness (adult pneumococcus)", "any hearing loss 20–30 %", "18"),
        ("Cerebral infarction (pneumococcus)", "15–25 %", "9"),
        ("Antibiotic delay", "outcome worsens beyond 6 h", "20"),
        ("Cephalosporin-resistant MIC", "up to 4 mg/L", "14"),
    ]
    for a, b, c in rows:
        W(f"| {a} | {b} | {c} |")
    W("")
    with open("abm_references.md", "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"{total} papers saved", file=sys.stderr)


if __name__ == "__main__":
    main()
