"""GSE155114: stream one dense per-library counts TSV (rows = nuclei, columns = genes).
Keeps per nucleus: library size (sum over ALL genes), genes detected, counts of the pre-specified genes.
usage: extract_gse155114.py <counts.tsv> <out.tsv.gz>   (run_all.sh runs it over every file in parallel)"""
import sys, os, re, numpy as np, pandas as pd
from genes import KEEP
f, out = sys.argv[1], sys.argv[2]
grp, donor, lib = re.match(r"(CTE|Control)_(\d+)\.(\d+)\.counts\.tsv$", os.path.basename(f)).groups()
hdr = open(f).readline().rstrip("\n").split("\t")
keep = [g for g in KEEP if g in hdr]
parts = []
for ch in pd.read_csv(f, sep="\t", dtype={c: np.float32 for c in hdr[1:]}, index_col=0, chunksize=3000):
    d = ch[keep].copy(); d.insert(0, "n_genes", (ch.values > 0).sum(1)); d.insert(0, "lib_size", ch.values.sum(1)); parts.append(d)
df = pd.concat(parts); df.index.name = "barcode"; df = df.reset_index()
df.insert(0, "library", f"{grp}_{donor}.{lib}"); df.insert(0, "donor", f"{grp}_{donor}"); df.insert(0, "group", grp)
df.to_csv(out, sep="\t", index=False, compression="gzip")
