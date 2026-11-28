"""H2 test on GSE319253 (bulk RNA-seq, superior frontal cortex, 6 CTE vs 6 control donors, 2 technical replicates each).
H2 prediction (inverse-comorbidity direction): PIN1, cell-cycle and DNA-repair expression LOWER in CTE than controls."""
import json, pathlib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
import statsmodels.api as sm
from genesets import PIN1, CELL_CYCLE, DNA_REPAIR, NEURON, ASTRO, MICROGLIA

D = pathlib.Path(__file__).parent
cnt = pd.read_csv(D / "data/GSE319253_featureCounts_count_matrix_human.txt.gz", sep="\t", index_col=0)
cnt.index = cnt.index.str.split(".").str[0]
gmap = json.loads((D / "gene_map.json").read_text())
logcpm = np.log2(cnt / cnt.sum() * 1e6 + 1)
donor = logcpm.T.groupby(lambda c: c.rsplit("_REP", 1)[0]).mean().T          # collapse technical replicates
grp = pd.Series(["CTE" if c.startswith("CTE") else "Control" for c in donor.columns], index=donor.columns)
cte, ctl = grp[grp == "CTE"].index, grp[grp == "Control"].index

def gene_stats(sym):
    v = donor.loc[gmap[sym]]
    return {"log2FC_CTE_vs_ctrl": round(float(v[cte].mean() - v[ctl].mean()), 3),
            "welch_p": float(f"{ttest_ind(v[cte], v[ctl], equal_var=False).pvalue:.3g}"),
            "mwu_p": float(f"{mannwhitneyu(v[cte], v[ctl]).pvalue:.3g}")}

def set_score(genes):
    m = donor.loc[[gmap[g] for g in genes]]
    z = m.sub(m.mean(1), axis=0).div(m.std(1).replace(0, np.nan), axis=0)
    return z.mean()

res = {"dataset": "GSE319253", "n_donors": {"CTE": int(len(cte)), "Control": int(len(ctl))},
       "library_size_millions": {c: round(float(cnt[c].sum() / 1e6), 2) for c in cnt.columns}, "genes": {}, "sets": {}}
for name, genes in [("PIN1", PIN1), ("CELL_CYCLE", CELL_CYCLE), ("DNA_REPAIR", DNA_REPAIR), ("NEURON", NEURON), ("ASTRO", ASTRO), ("MICROGLIA", MICROGLIA)]:
    for g in genes:
        res["genes"][g] = gene_stats(g)
    s = set_score(genes)
    res["sets"][name] = {"n_genes": len(genes), "mean_z_CTE": round(float(s[cte].mean()), 3), "mean_z_Control": round(float(s[ctl].mean()), 3),
                         "mwu_p_two_sided": float(f"{mannwhitneyu(s[cte], s[ctl]).pvalue:.3g}"),
                         "n_genes_lower_in_CTE": int(sum(res["genes"][g]["log2FC_CTE_vs_ctrl"] < 0 for g in genes))}
# Exploratory: PIN1 adjusted for neuronal marker score (composition)
neu = set_score(NEURON)
X = sm.add_constant(pd.DataFrame({"cte": (grp == "CTE").astype(float), "neuron_score": neu}))
fit = sm.OLS(donor.loc[gmap["PIN1"]], X).fit()
res["pin1_adjusted_for_neuron_score"] = {"coef_cte": round(float(fit.params["cte"]), 3), "ci95": [round(float(x), 3) for x in fit.conf_int().loc["cte"]],
                                         "p": float(f"{fit.pvalues['cte']:.3g}"), "coef_neuron_score": round(float(fit.params["neuron_score"]), 3)}
(D / "results").mkdir(exist_ok=True)
(D / "results/h2_results.json").write_text(json.dumps(res, indent=2))
print(json.dumps({k: v for k, v in res.items() if k != "library_size_millions"}, indent=2))
