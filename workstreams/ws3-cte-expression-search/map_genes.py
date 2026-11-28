"""Map gene symbols to Ensembl gene IDs via the Ensembl REST API; writes gene_map.json (committed)."""
import json, urllib.request, pathlib
from genesets import ALL
req = urllib.request.Request("https://rest.ensembl.org/lookup/symbol/homo_sapiens", data=json.dumps({"symbols": ALL}).encode(),
                             headers={"Content-Type": "application/json", "Accept": "application/json"})
r = json.load(urllib.request.urlopen(req, timeout=60))
m = {s: r[s]["id"] for s in ALL if s in r}
pathlib.Path(__file__).with_name("gene_map.json").write_text(json.dumps(m, indent=1, sort_keys=True))
print(len(m), "of", len(ALL), "mapped; missing:", [s for s in ALL if s not in m])
