"""Build the pinned analysis table from the public cBioPortal API (open-tier TCGA PanCancer Atlas).

Writes data/glioma_mapt_idh_os.tsv (one row per primary tumour sample) and data/SHA256SUMS.
The table is derived, de-identified open-tier data; it is NOT committed to the repo (see README).
run.sh verifies the checksum so results are tied to one exact snapshot."""
import json, hashlib, pathlib, urllib.request, datetime

API = "https://www.cbioportal.org/api"
STUDIES = ["lgg_tcga_pan_can_atlas_2018", "gbm_tcga_pan_can_atlas_2018"]
MAPT, IDH1, IDH2 = 4137, 3417, 3418
OUT = pathlib.Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)

def call(path, body=None):
    req = urllib.request.Request(API + path, data=json.dumps(body).encode() if body is not None else None,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"},
                                 method="POST" if body is not None else "GET")
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)

def hotspot(entrez, change):
    return (entrez == IDH1 and change.startswith("R132")) or (entrez == IDH2 and (change.startswith("R172") or change.startswith("R140")))

rows, meta = [], {"fetched_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "api": API, "studies": {}}
for s in STUDIES:
    st = call(f"/studies/{s}")
    meta["studies"][s] = {"name": st.get("name"), "importDate": st.get("importDate")}
    expr = {d["sampleId"]: d["value"] for d in call(f"/molecular-profiles/{s}_rna_seq_v2_mrna/molecular-data/fetch",
                                                    {"sampleListId": f"{s}_all", "entrezGeneIds": [MAPT]})}
    seq = set(call(f"/sample-lists/{s}_sequenced/sample-ids"))
    muts = {}
    for m in call(f"/molecular-profiles/{s}_mutations/mutations/fetch?projection=SUMMARY",
                  {"sampleListId": f"{s}_all", "entrezGeneIds": [IDH1, IDH2]}):
        if hotspot(m["entrezGeneId"], m["proteinChange"]):
            muts.setdefault(m["sampleId"], set()).add(("IDH1" if m["entrezGeneId"] == IDH1 else "IDH2") + "_" + m["proteinChange"])
    clin = {}
    for c in call(f"/studies/{s}/clinical-data?clinicalDataType=PATIENT&projection=SUMMARY&pageSize=100000"):
        if c["clinicalAttributeId"] in ("OS_MONTHS", "OS_STATUS", "AGE"):
            clin.setdefault(c["patientId"], {})[c["clinicalAttributeId"]] = c["value"]
    for smp in sorted(expr):
        if not smp.endswith("-01"):          # primary tumour only
            continue
        pid = smp[:12]
        idh = ("mut" if smp in muts else "wt") if smp in seq else "NA"
        cl = clin.get(pid, {})
        rows.append([s.split("_")[0].upper(), smp, f"{expr[smp]}", idh, ";".join(sorted(muts.get(smp, []))) or "-",
                     cl.get("OS_MONTHS", "NA"), cl.get("OS_STATUS", "NA"), cl.get("AGE", "NA")])

tsv = OUT / "glioma_mapt_idh_os.tsv"
with open(tsv, "w") as fh:
    fh.write("cohort\tsample\tmapt_rsem\tidh_hotspot\tidh_variant\tos_months\tos_status\tage\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
(OUT / "snapshot_meta.json").write_text(json.dumps(meta, indent=2, sort_keys=True))
(OUT / "SHA256SUMS").write_text(f"{hashlib.sha256(tsv.read_bytes()).hexdigest()}  {tsv.name}\n")
print(f"wrote {len(rows)} rows; sha256 {hashlib.sha256(tsv.read_bytes()).hexdigest()}")
