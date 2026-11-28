"""PIN1 CRISPR dependency by lineage, DepMap 26Q3 (portal CRISPRGeneEffect + Model).
Same metrics as the 22Q2 analysis in ../analyze.py. Needs ~2 GB RAM."""
import json, pathlib, numpy as np, pandas as pd
from scipy.stats import kruskal

D = pathlib.Path(__file__).parent
ROOT = D.parent
EFF = ROOT / "data/26Q3/CRISPRGeneEffect.csv"
MOD = ROOT / "data/26Q3/Model.csv"
REF = ["PIN1", "PIN4", "PPIA", "PCNA", "CDK1", "POLR2A", "RPL3"]

eff = pd.read_csv(EFF, index_col=0)
# columns are "SYMBOL (Entrez)"
sym = pd.Index([c.split(" (")[0] for c in eff.columns])
if sym.duplicated().any():
    raise SystemExit("duplicate gene symbols in CRISPRGeneEffect header: " + ",".join(sym[sym.duplicated()]))
eff.columns = sym
gene_median = eff.median(axis=0, skipna=True)

md = pd.read_csv(MOD, usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "CellLineName", "ModelType"])
md = md.rename(columns={"ModelID": "depmap_id", "OncotreeLineage": "lineage", "OncotreePrimaryDisease": "primary_disease"})

res = {
    "release": "DepMap 26Q3 portal files CRISPRGeneEffect.csv and Model.csv (downloaded 2026-11-28 after human verification)",
    "n_lines": int(eff.shape[0]),
    "n_genes": int(eff.shape[1]),
    "lineage_field": "Model.OncotreeLineage (26Q3 ontology; not the same labels as 22Q2 'lineage')",
    "reference_genes": {},
}
res["reference_genes_absent_from_release"] = [g for g in REF if g not in eff.columns]
for g in [g for g in REF if g in eff.columns]:
    v = eff[g].dropna()
    res["reference_genes"][g] = {
        "n_lines": int(v.size),
        "median": round(float(v.median()), 3),
        "q25": round(float(v.quantile(0.25)), 3),
        "q75": round(float(v.quantile(0.75)), 3),
        "min": round(float(v.min()), 3),
        "n_lt_-0.5": int((v < -0.5).sum()),
        "frac_lt_-0.5": round(float((v < -0.5).mean()), 3),
        "frac_lt_-1": round(float((v < -1).mean()), 3),
        "median_rank_pct_most_essential_first": round(float((gene_median < gene_median[g]).mean() * 100), 2),
    }

p = eff[["PIN1"]].rename(columns={"PIN1": "dependency"}).reset_index().rename(columns={eff.index.name or "index": "depmap_id"})
# the index name of a default csv first column is often None; pandas names it after the first header cell, which was empty -> "Unnamed: 0" avoided by index_col=0
if p.columns[0] != "depmap_id":
    p = p.rename(columns={p.columns[0]: "depmap_id"})
p = p.merge(md, on="depmap_id", how="left")
p = p.dropna(subset=["dependency"])

lin = p.groupby("lineage").dependency
tab = pd.DataFrame({
    "n": lin.size(),
    "median": lin.median().round(3),
    "frac_lt_-0.5": lin.apply(lambda s: (s < -0.5).mean()).round(3),
    "n_lt_-0.5": lin.apply(lambda s: int((s < -0.5).sum())),
    "min": lin.min().round(3),
})
tab = tab[tab.n >= 10].sort_values("median")
groups = [p.loc[p.lineage == l, "dependency"].values for l in tab.index]
kw = kruskal(*groups)
res["pin1_by_lineage_n>=10"] = tab.reset_index().to_dict(orient="records")
res["kruskal_wallis_pin1_across_lineages"] = {
    "n_lineages": int(len(groups)),
    "H": round(float(kw.statistic), 2),
    "p": float(f"{kw.pvalue:.3g}"),
}
top = p.nsmallest(10, "dependency")[["depmap_id", "CellLineName", "lineage", "primary_disease", "ModelType", "dependency"]].round(3)
res["pin1_most_dependent_lines_top10"] = top.to_dict(orient="records")
low = p[p.dependency < -0.5].sort_values("dependency")
cols = ["depmap_id", "CellLineName", "lineage", "primary_disease", "ModelType", "dependency"]
res["pin1_lines_lt_-0.5"] = low[cols].round(3).to_dict(orient="records")
res["n_lines_unmatched_to_Model"] = int(p["lineage"].isna().sum())
# Post-hoc: after seeing 14 lines < -0.5, we split by ModelType (not pre-specified).
mt = []
for t, g in p.groupby("ModelType"):
    v = g.dependency.dropna()
    mt.append({"ModelType": t, "n_lines": int(v.size), "median": round(float(v.median()), 3),
               "min": round(float(v.min()), 3), "n_lt_-0.5": int((v < -0.5).sum()),
               "frac_lt_-0.5": round(float((v < -0.5).mean()), 3)})
res["posthoc_pin1_by_ModelType"] = {"note": "Added after seeing the 14 lines < -0.5; not a pre-specified analysis.", "rows": mt}

# side-by-side with the committed 22Q2 results
old = json.loads((ROOT / "results/results.json").read_text())
cmp = {"releases": {"22Q2": old["release"], "26Q3": res["release"]},
       "n_genes": {"22Q2": old["n_genes"], "26Q3": res["n_genes"]},
       "n_lines_PIN1": {}, "genes": {}, "kruskal_wallis_pin1_across_lineages": {
           "22Q2": old["kruskal_wallis_pin1_across_lineages"],
           "26Q3": res["kruskal_wallis_pin1_across_lineages"],
           "note": "Lineage labels are not the same ontology (22Q2 DepMap lineage vs 26Q3 OncotreeLineage); H and p are not a paired test.",
       }}
genes = sorted(set(old["reference_genes"]) | set(res["reference_genes"]) | set(REF))
for g in genes:
    a, b = old["reference_genes"].get(g), res["reference_genes"].get(g)
    row = {"22Q2": a, "26Q3": b}
    if a and b:
        row["delta_median_26Q3_minus_22Q2"] = round(b["median"] - a["median"], 3)
        row["delta_rank_pct_26Q3_minus_22Q2"] = round(b["median_rank_pct_most_essential_first"] - a["median_rank_pct_most_essential_first"], 2)
    cmp["genes"][g] = row
cmp["n_lines_PIN1"] = {"22Q2": old["reference_genes"]["PIN1"]["n_lines"], "26Q3": res["reference_genes"]["PIN1"]["n_lines"]}
# 22Q2 results.json did not store an explicit overall min; recover from its top-10 list (most dependent line)
cmp["PIN1_min"] = {"22Q2": old["pin1_most_dependent_lines_top10"][0]["dependency"], "26Q3": res["reference_genes"]["PIN1"]["min"]}
cmp["PIN1_n_lt_-0.5"] = {
    "22Q2": int(round(old["reference_genes"]["PIN1"]["frac_lt_-0.5"] * old["reference_genes"]["PIN1"]["n_lines"])),
    "26Q3": res["reference_genes"]["PIN1"]["n_lt_-0.5"],
}
cmp["PIN1_median"] = {"22Q2": old["reference_genes"]["PIN1"]["median"], "26Q3": res["reference_genes"]["PIN1"]["median"]}
cmp["PIN1_median_rank_pct_most_essential_first"] = {
    "22Q2": old["reference_genes"]["PIN1"]["median_rank_pct_most_essential_first"],
    "26Q3": res["reference_genes"]["PIN1"]["median_rank_pct_most_essential_first"],
}

(D / "results").mkdir(exist_ok=True)
tab.to_csv(D / "results/pin1_by_lineage.csv")
(D / "results/results.json").write_text(json.dumps(res, indent=2))
(D / "results/comparison_22Q2_vs_26Q3.json").write_text(json.dumps(cmp, indent=2))
print(json.dumps({"26Q3": res, "comparison": cmp}, indent=2))
