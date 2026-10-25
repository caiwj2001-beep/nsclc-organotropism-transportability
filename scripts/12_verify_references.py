# -*- coding: utf-8 -*-
"""Verify every reference against PubMed (esummary) and dump machine-readable metadata."""
import json, os, time, urllib.parse, urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
PROXY = {"http": "http://127.0.0.1:7890", "https": "http://127.0.0.1:7890"}
OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "06_manuscript"))
os.makedirs(OUT, exist_ok=True)

PMIDS = ["37084736", "29858031", "37121400", "37771447", "35120664", "41589477",
         "41100253", "41747607", "42909660", "35100709", "28051122", "33188911",
         "41294679", "39506116", "22960745", "35598205", "35410334", "42085620",
         "40650339", "41743381"]


def get(url):
    last = None
    for _ in range(3):
        for up in (False, True):
            try:
                op = urllib.request.build_opener(urllib.request.ProxyHandler(PROXY if up else {}))
                op.addheaders = [("User-Agent", "Mozilla/5.0")]
                with op.open(url, timeout=60) as r:
                    return json.loads(r.read().decode("utf-8", "replace"))
            except Exception as e:
                last = e
                time.sleep(1.0)
    print("FAIL", last)
    return None


# Skoulidis STK11/LKB1 immunoresistance (find PMID first)
q = get(E + "esearch.fcgi?" + urllib.parse.urlencode(
    {"db": "pubmed", "term": "Skoulidis[Author] AND STK11 AND PD-1 inhibitor resistance KRAS",
     "retmax": 5, "retmode": "json"}))
extra = q["esearchresult"]["idlist"] if q else []
print("Skoulidis candidates:", extra)
allids = PMIDS + [i for i in extra[:1] if i not in PMIDS]

meta = {}
for i in range(0, len(allids), 10):
    chunk = allids[i:i + 10]
    raw = get(E + "esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "pubmed", "id": ",".join(chunk), "retmode": "json"}))
    if raw:
        meta.update(raw.get("result", {}))
    time.sleep(0.4)

rows = []
for pid in allids:
    a = meta.get(pid, {})
    if not a:
        print("MISSING", pid); continue
    au = [x["name"] for x in (a.get("authors") or [])]
    rows.append({
        "pmid": pid,
        "authors": au,
        "first3": ", ".join(au[:3]) + (", et al" if len(au) > 3 else ""),
        "title": a.get("title", "").rstrip("."),
        "journal": a.get("source", ""),
        "year": (a.get("pubdate") or "")[:4],
        "volume": a.get("volume", ""), "issue": a.get("issue", ""),
        "pages": a.get("pages", ""),
        "doi": a.get("elocationid", ""),
        "pubtype": a.get("pubtype", []),
    })

json.dump(rows, open(os.path.join(OUT, "references_verified.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
for r in rows:
    print(f"{r['pmid']:9s} {r['year']} {r['journal'][:26]:26s} {r['volume']}:{r['pages']:14s} {r['first3'][:38]:38s} {r['title'][:70]}")
print("\nwrote references_verified.json |", len(rows), "records")
