"""Rerunnable search of NCBI GEO and Synapse for public CTE brain expression datasets. Writes results/hits.json."""
import json, time, urllib.parse, urllib.request, datetime, pathlib
D = pathlib.Path(__file__).parent
UA = {"User-Agent": "CureCancer-research-search/1.0"}
def get(url, body=None):
    req = urllib.request.Request(url, data=json.dumps(body).encode() if body else None,
                                 headers={**UA, "Content-Type": "application/json", "Accept": "application/json"},
                                 method="POST" if body else "GET")
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)

GEO_QUERIES = ['"chronic traumatic encephalopathy"', 'CTE[All Fields] AND brain[All Fields] AND Homo sapiens[Organism]',
               '"repetitive head impact"', '"traumatic encephalopathy"', '"chronic traumatic encephalopathy" AND (expression profiling by high throughput sequencing[DataSet Type] OR expression profiling by array[DataSet Type])']
out = {"searched_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "geo": {}, "synapse": {}}
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for q in GEO_QUERIES:
    r = get(E + "esearch.fcgi?db=gds&retmode=json&retmax=200&term=" + urllib.parse.quote(q))["esearchresult"]
    ids = r.get("idlist", [])
    hits = []
    if ids:
        s = get(E + "esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids))["result"]
        for u in s.get("uids", []):
            x = s[u]
            hits.append({"accession": x.get("accession"), "type": x.get("entrytype"), "gdstype": x.get("gdstype"), "taxon": x.get("taxon"),
                         "n_samples": x.get("n_samples"), "title": x.get("title"), "pdat": x.get("pdat"),
                         "summary_snippet": (x.get("summary") or "")[:300]})
    out["geo"][q] = {"count": int(r.get("count", 0)), "hits": hits}
    time.sleep(0.5)
SYN = "https://repo-prod.prod.sagebase.org/repo/v1/search"
for q in ["chronic traumatic encephalopathy", "CTE RNA-seq", "traumatic encephalopathy brain"]:
    try:
        r = get(SYN, {"queryTerm": q.split(), "size": 50, "start": 0})
        out["synapse"][q] = {"found": r.get("found"), "hits": [{"id": h.get("id"), "name": h.get("name"), "node_type": h.get("node_type")} for h in r.get("hits", [])]}
    except Exception as e:
        out["synapse"][q] = {"error": f"{type(e).__name__}: {e}"}
(D / "results/hits.json").write_text(json.dumps(out, indent=2))
for q, v in out["geo"].items(): print("GEO", v["count"], q)
for q, v in out["synapse"].items(): print("SYN", v.get("found", v.get("error")), q)
