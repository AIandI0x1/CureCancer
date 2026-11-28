import json, pathlib
r=json.load(open(pathlib.Path(__file__).parent/'results/results.json'))
for k in ["GSE155114_CTE_vs_Control","GSE155114_CTE_vs_Control_strictQC_sensitivity","GSE261807_CTE_vs_Control","GSE261807_RHI_vs_Control"]:
    x=r[k]; print("\n==",k, "nuclei after QC:",x['n_nuclei_qc'], "assigned:",x['assigned_fraction'], x['nuclei_by_type'], "missing genes:",x['genes_missing'])
    for t,v in x['types'].items():
        if 'genes' not in v: print(' ',t,v['n_donors'],v['status']); continue
        p=v['genes']['PIN1']
        print(f"  {t:16s} donors {v['n_donors']} med nuclei {v['median_nuclei_per_donor']}")
        print("     PIN1:",p)
        for s,d in v['sets'].items(): print("    ",s,d)
print("\nMETA",json.dumps(r['PIN1_fixed_effect_meta_CTE_vs_control'],indent=0))
