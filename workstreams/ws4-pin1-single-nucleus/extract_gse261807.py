"""GSE261807: read one CellRanger filtered matrix (genes x nuclei, MatrixMarket).
Keeps per nucleus: library size (all genes), genes detected, counts of the pre-specified genes.
usage: extract_gse261807.py <dir> <GSM_prefix e.g. GSM8152382_S1_CTE> <out.tsv.gz>"""
import sys, re, numpy as np, pandas as pd
from genes import KEEP
d, pre, out = sys.argv[1:4]
gsm, sid, grp = re.match(r"(GSM\d+)_(S\d+)_(CTE|RHI|Control)$", pre).groups()
feat = pd.read_csv(f"{d}/{pre}_features.tsv.gz", sep="\t", header=None, usecols=[0, 1], names=["ens", "sym"])
bc = pd.read_csv(f"{d}/{pre}_barcodes.tsv.gz", header=None)[0].values
m = pd.read_csv(f"{d}/{pre}_matrix.mtx.gz", sep=" ", comment="%", header=None, names=["g", "c", "v"], skiprows=2, dtype=np.int32)
nc = len(bc)
lib = np.bincount(m.c - 1, weights=m.v, minlength=nc); ng = np.bincount(m.c - 1, minlength=nc)
df = pd.DataFrame({"group": grp, "donor": sid, "library": sid, "barcode": bc, "lib_size": lib, "n_genes": ng})
sym = feat.sym.values
for g in KEEP:
    gi = np.where(sym == g)[0]
    if len(gi) != 1: continue                      # absent or ambiguous symbol: skip, reported by the analysis
    s = m[m.g == gi[0] + 1]
    df[g] = np.bincount(s.c - 1, weights=s.v, minlength=nc)
df.to_csv(out, sep="\t", index=False, compression="gzip")
