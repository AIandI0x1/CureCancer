"""WS4: does the CTE-associated PIN1 decrease (WS3, bulk) hold inside cell types? Pseudobulk per donor x cell type.
Datasets: GSE155114 (inDrops snRNA-seq, 8 CTE vs 8 control donors) and GSE261807 (10x snRNA-seq, CTE / RHI / control).
All thresholds below were fixed before group comparisons were run."""
import glob, json, pathlib, numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, norm
from statsmodels.stats.multitest import multipletests
from genes import PIN1, CELL_CYCLE, DNA_REPAIR, MARKERS, MARKER_GENES

D = pathlib.Path(__file__).parent
# nucleus QC. GSE261807 = CellRanger-filtered 10x matrices. GSE155114 = GEO dense matrices that still contain empty droplets
# (median 74 UMIs/barcode), so the cut must be low-ish: primary 500 UMI / 300 genes; 1000 / 400 is kept as a sensitivity analysis.
# NOTE (honesty): 1000/400 was my first pass for both; it left only ~10.6k GSE155114 nuclei, so the primary cut for that dataset was lowered
# to 500/300 AFTER seeing that pass. Both are reported.
QC = {"GSE261807": (1000, 400), "GSE155114": (500, 300), "GSE155114_strict": (1000, 400)}
Z_MIN, Z_GAP = 1.0, 0.5                 # marker assignment: top z >= 1 and gap to runner-up >= 0.5
MIN_NUCLEI = 50                          # per donor x cell type pseudobulk
MIN_DONORS = 3                           # per group
TYPES = list(MARKERS)
SETS = {"PIN1": PIN1, "CELL_CYCLE": CELL_CYCLE, "DNA_REPAIR": DNA_REPAIR}

def load(pattern, qc):
    min_lib, min_genes = QC[qc]
    df = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(str(D / pattern)))], ignore_index=True)
    return df[(df.lib_size >= min_lib) & (df.n_genes >= min_genes)].reset_index(drop=True)

def assign(df):
    norm_ = lambda g: np.log1p(df[g] / df.lib_size * 1e4)
    S = pd.DataFrame({t: np.mean([norm_(g) for g in gs if g in df], axis=0) for t, gs in MARKERS.items()})
    Z = (S - S.mean()) / S.std()
    srt = np.sort(Z.values, axis=1)
    ok = (srt[:, -1] >= Z_MIN) & ((srt[:, -1] - srt[:, -2]) >= Z_GAP)
    lab = np.where(ok, Z.idxmax(axis=1), "Unassigned")
    return lab, S

def qpois(c_case, c_ctl, t_case, t_ctl):
    """quasi-Poisson GLM, log2 rate ratio case vs control with library-size offset."""
    y = np.r_[c_case, c_ctl]; off = np.log(np.r_[t_case, t_ctl]); x = np.r_[np.ones(len(c_case)), np.zeros(len(c_ctl))]
    if y.sum() < 20: return None
    f = sm.GLM(y, sm.add_constant(x), family=sm.families.Poisson(), offset=off).fit(scale="X2")
    b, se = f.params[1], f.bse[1]; l2 = np.log2(np.e)
    return {"log2RR": b * l2, "se_log2": se * l2, "ci95": [(b - 1.96 * se) * l2, (b + 1.96 * se) * l2], "p": float(f.pvalues[1])}

def run(name, df, case, ctl="Control"):
    df = df.copy(); df["cell_type"], _ = assign(df)
    out = {"n_nuclei_qc": int(len(df)), "assigned_fraction": round(float((df.cell_type != "Unassigned").mean()), 3),
           "nuclei_by_type": df.cell_type.value_counts().to_dict(), "types": {}}
    g_all = PIN1 + CELL_CYCLE + DNA_REPAIR
    g_all = [g for g in g_all if g in df]
    out["genes_missing"] = [g for g in PIN1 + CELL_CYCLE + DNA_REPAIR if g not in df]
    rows = []
    for t in TYPES:
        sub = df[df.cell_type == t]
        pb = sub.groupby(["group", "donor"]).agg(n=("barcode", "size"), tot=("lib_size", "sum"), **{g: (g, "sum") for g in g_all}).reset_index()
        pb = pb[pb.n >= MIN_NUCLEI]
        a, b = pb[pb.group == case], pb[pb.group == ctl]
        res = {"n_donors": {case: int(len(a)), ctl: int(len(b))}, "median_nuclei_per_donor": {case: float(a.n.median()) if len(a) else None, ctl: float(b.n.median()) if len(b) else None}}
        if len(a) < MIN_DONORS or len(b) < MIN_DONORS:
            res["status"] = "skipped: fewer than %d donors with >= %d nuclei in a group" % (MIN_DONORS, MIN_NUCLEI); out["types"][t] = res; continue
        cpm = lambda p, g: np.log2(p[g] / p.tot * 1e6 + 1)
        res["genes"] = {}
        for g in g_all:
            q = qpois(a[g].values, b[g].values, a.tot.values, b.tot.values)
            if q is None: res["genes"][g] = {"status": "too few counts (<20 total)"}; continue
            mw = mannwhitneyu(cpm(a, g), cpm(b, g)).pvalue
            res["genes"][g] = {"log2RR": round(q["log2RR"], 3), "ci95": [round(x, 3) for x in q["ci95"]], "qpois_p": float(f"{q['p']:.3g}"),
                               "mwu_p": float(f"{mw:.3g}"), "total_counts": int(a[g].sum() + b[g].sum())}
        res["sets"] = {}
        for s, gs in list(SETS.items())[1:]:
            ok = [g for g in gs if g in res["genes"] and "log2RR" in res["genes"][g] and res["genes"][g]["total_counts"] >= 100]
            if len(ok) < 3: res["sets"][s] = {"status": "fewer than 3 genes with >=100 total counts"}; continue
            M = pd.concat([cpm(pb, g).rename(g) for g in ok], axis=1); Z = (M - M.mean()) / M.std().replace(0, np.nan)
            sc = Z.mean(axis=1)
            sa, sb = sc[pb.group == case], sc[pb.group == ctl]
            res["sets"][s] = {"n_genes_used": len(ok), "n_genes_lower_in_case": int(sum(res["genes"][g]["log2RR"] < 0 for g in ok)),
                              "mean_z_case": round(float(sa.mean()), 3), "mean_z_ctl": round(float(sb.mean()), 3), "mwu_p": float(f"{mannwhitneyu(sa, sb).pvalue:.3g}")}
        p1 = res["genes"].get("PIN1", {})
        if "log2RR" in p1:
            rows.append((t, p1["qpois_p"]))
            # donor-level table for transparency
            res["pin1_donor_cpm"] = {r.donor: round(float(r.PIN1 / r.tot * 1e6), 2) for r in pb.itertuples()}
        out["types"][t] = res
    if rows:
        adj = multipletests([p for _, p in rows], method="fdr_bh")[1]
        for (t, _), q in zip(rows, adj): out["types"][t]["genes"]["PIN1"]["bh_adj_p_across_types"] = float(f"{q:.3g}")
    return out

res = {"thresholds": {"qc_min_umi_min_genes": QC, "z_min": Z_MIN, "z_gap": Z_GAP, "min_nuclei_per_donor_type": MIN_NUCLEI, "min_donors_per_group": MIN_DONORS}}
d1 = load("data/n155/*.tsv.gz", "GSE155114"); res["GSE155114_CTE_vs_Control"] = run("GSE155114", d1, "CTE")
res["GSE155114_CTE_vs_Control_strictQC_sensitivity"] = run("GSE155114", load("data/n155/*.tsv.gz", "GSE155114_strict"), "CTE")
d2 = load("data/n261/*.tsv.gz", "GSE261807")
res["GSE261807_CTE_vs_Control"] = run("GSE261807", d2, "CTE")
res["GSE261807_RHI_vs_Control"] = run("GSE261807", d2, "RHI")
# fixed-effect meta of PIN1 log2RR (CTE vs control) across the two datasets, per cell type
meta = {}
for t in TYPES:
    xs = []
    for k in ["GSE155114_CTE_vs_Control", "GSE261807_CTE_vs_Control"]:
        g = res[k]["types"][t].get("genes", {}).get("PIN1", {})
        if "log2RR" in g:
            se = (g["ci95"][1] - g["ci95"][0]) / (2 * 1.96); xs.append((g["log2RR"], se))
    if len(xs) == 2:
        w = np.array([1 / s**2 for _, s in xs]); b = np.array([x for x, _ in xs]); m = (w * b).sum() / w.sum(); se = (1 / w.sum()) ** .5
        meta[t] = {"log2RR": round(float(m), 3), "ci95": [round(float(m - 1.96 * se), 3), round(float(m + 1.96 * se), 3)], "p": float(f"{2 * norm.sf(abs(m / se)):.3g}"),
                   "datasets_same_direction": bool(np.sign(xs[0][0]) == np.sign(xs[1][0]))}
res["PIN1_fixed_effect_meta_CTE_vs_control"] = meta
(D / "results").mkdir(exist_ok=True); (D / "results/results.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
