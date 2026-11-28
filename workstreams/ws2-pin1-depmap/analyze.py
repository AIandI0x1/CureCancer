"""PIN1 CRISPR dependency by lineage, DepMap 22Q2 (Bioconductor depmap crispr_22Q2). Needs ~8 GB RAM."""
import json, pathlib, numpy as np, pandas as pd, pyreadr
from scipy.stats import kruskal

D = pathlib.Path(__file__).parent
cr = list(pyreadr.read_r(str(D / "data/crispr_22Q2.rda")).values())[0][["depmap_id", "gene_name", "dependency"]]
md = list(pyreadr.read_r(str(D / "data/metadata_22Q2.rda")).values())[0][["depmap_id", "lineage", "primary_disease"]]
gene_median = cr.groupby("gene_name")["dependency"].median()

REF = ["PIN1", "PIN4", "PPIA", "PCNA", "CDK1", "POLR2A", "RPL3"]
sub = cr[cr.gene_name.isin(REF)].merge(md, on="depmap_id", how="left")
del cr

res = {"release": "DepMap 22Q2 via Bioconductor ExperimentHub (EH7554 crispr_22Q2, EH7558 metadata_22Q2)",
       "n_genes": int(gene_median.size), "reference_genes": {}}
res["reference_genes_absent_from_release"] = [g for g in REF if g not in gene_median.index]
for g in [g for g in REF if g in gene_median.index]:
    v = sub[sub.gene_name == g].dependency.dropna()
    res["reference_genes"][g] = {"n_lines": int(v.size), "median": round(float(v.median()), 3),
        "q25": round(float(v.quantile(.25)), 3), "q75": round(float(v.quantile(.75)), 3),
        "frac_lt_-0.5": round(float((v < -0.5).mean()), 3), "frac_lt_-1": round(float((v < -1).mean()), 3),
        "median_rank_pct_most_essential_first": round(float((gene_median < gene_median[g]).mean() * 100), 2)}
p = sub[sub.gene_name == "PIN1"].dropna(subset=["dependency"])
lin = p.groupby("lineage").dependency
tab = pd.DataFrame({"n": lin.size(), "median": lin.median().round(3), "frac_lt_-0.5": lin.apply(lambda s: (s < -0.5).mean()).round(3),
                    "min": lin.min().round(3)})
tab = tab[tab.n >= 10].sort_values("median")
groups = [p[p.lineage == l].dependency.values for l in tab.index]
kw = kruskal(*groups)
res["pin1_by_lineage_n>=10"] = tab.reset_index().to_dict(orient="records")
res["kruskal_wallis_pin1_across_lineages"] = {"n_lineages": int(len(groups)), "H": round(float(kw.statistic), 2), "p": float(f"{kw.pvalue:.3g}")}
res["pin1_most_dependent_lines_top10"] = p.nsmallest(10, "dependency")[["depmap_id", "lineage", "primary_disease", "dependency"]].round(3).to_dict(orient="records")
(D / "results").mkdir(exist_ok=True)
tab.to_csv(D / "results/pin1_by_lineage.csv")
(D / "results/results.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
