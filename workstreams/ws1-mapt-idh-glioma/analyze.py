"""Pre-specified analysis: MAPT expression vs IDH hotspot status and overall survival (TCGA LGG+GBM)."""
import json, pathlib, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

D = pathlib.Path(__file__).parent
df = pd.read_csv(D / "data/glioma_mapt_idh_os.tsv", sep="\t", na_values=["NA"])
df["log2_mapt"] = np.log2(df["mapt_rsem"] + 1)
df["event"] = df["os_status"].str.startswith("1").astype(float)
df.loc[df["os_status"].isna(), "event"] = np.nan
res = {"n_samples": int(len(df)), "n_by_cohort_idh": {f"{k[0]}_{k[1]}": int(v) for k, v in df.groupby(["cohort", "idh_hotspot"]).size().items()}}

def expr_cmp(sub, label):
    m, w = sub[sub.idh_hotspot == "mut"].log2_mapt, sub[sub.idh_hotspot == "wt"].log2_mapt
    u = mannwhitneyu(m, w, alternative="two-sided")
    res[f"expr_{label}"] = {"n_mut": int(len(m)), "n_wt": int(len(w)),
        "median_log2_mut": round(float(m.median()), 3), "median_log2_wt": round(float(w.median()), 3),
        "median_diff_log2": round(float(m.median() - w.median()), 3),
        "auc_mut_gt_wt": round(float(u.statistic / (len(m) * len(w))), 3), "mannwhitney_p": float(f"{u.pvalue:.3g}")}

expr_cmp(df, "all_gliomas")
expr_cmp(df[df.cohort == "LGG"], "lgg_only")

sv = df.dropna(subset=["os_months", "event", "log2_mapt"]).copy()
sv = sv[sv.os_months > 0]
sv["mapt_z"] = (sv.log2_mapt - sv.log2_mapt.mean()) / sv.log2_mapt.std()

def cox(sub, covars, label):
    sub = sub.dropna(subset=covars)
    cph = CoxPHFitter().fit(sub[["os_months", "event"] + covars], "os_months", "event")
    r = cph.summary.loc["mapt_z"]
    res[f"cox_{label}"] = {"n": int(len(sub)), "events": int(sub.event.sum()), "covariates": covars,
        "hr_per_sd_mapt": round(float(r["exp(coef)"]), 3), "ci95": [round(float(r["exp(coef) lower 95%"]), 3), round(float(r["exp(coef) upper 95%"]), 3)],
        "p": float(f"{r['p']:.3g}")}

sv["idh_mut"] = (sv.idh_hotspot == "mut").astype(float)
sv.loc[sv.idh_hotspot.isna(), "idh_mut"] = np.nan
sv["gbm"] = (sv.cohort == "GBM").astype(float)
cox(sv, ["mapt_z"], "all_unadjusted")
cox(sv, ["mapt_z", "idh_mut", "age", "gbm"], "all_adj_idh_age_cohort")
cox(sv[sv.idh_hotspot == "mut"], ["mapt_z", "age"], "idh_mut_adj_age")
cox(sv[sv.idh_hotspot == "wt"], ["mapt_z", "age", "gbm"], "idh_wt_adj_age_cohort")

# Median split within IDH-mutant (descriptive)
mut = sv[sv.idh_hotspot == "mut"]
hi = mut.log2_mapt >= mut.log2_mapt.median()
lr = logrank_test(mut[hi].os_months, mut[~hi].os_months, mut[hi].event, mut[~hi].event)
res["logrank_idh_mut_median_split"] = {"n_high": int(hi.sum()), "n_low": int((~hi).sum()), "p": float(f"{lr.p_value:.3g}")}
_a = sv.dropna(subset=["age"]); rho = spearmanr(_a.log2_mapt, _a.age)
res["spearman_mapt_vs_age"] = {"n": int(len(_a)), "rho": round(float(rho.statistic), 3), "p": float(f"{rho.pvalue:.3g}")}

(D / "results").mkdir(exist_ok=True)
(D / "results/results.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
