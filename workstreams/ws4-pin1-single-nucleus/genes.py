"""Pre-specified gene lists (fixed before any CTE-vs-control result was computed).
H2 sets are identical to WS3. Cell-type markers are used only to assign nuclei."""
PIN1 = ["PIN1"]
CELL_CYCLE = ["CCND1", "CDK4", "CDK6", "CCNE1", "CDK2", "CCNA2", "CCNB1", "CDK1", "E2F1", "PCNA", "MCM2", "MKI67"]
DNA_REPAIR = ["ATM", "BRCA1", "RAD51", "XRCC1", "OGG1", "ERCC1", "PARP1", "MLH1", "MSH2", "XRCC5"]
MARKERS = {
    "Excitatory": ["SLC17A7", "SLC17A6", "SATB2", "CAMK2A"],
    "Inhibitory": ["GAD1", "GAD2", "SLC6A1"],
    "Astrocyte":  ["AQP4", "SLC1A2", "ALDH1L1", "GFAP"],
    "Microglia":  ["CSF1R", "P2RY12", "DOCK8", "APBB1IP"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOBP", "MOG"],
}
MARKER_GENES = sorted({g for v in MARKERS.values() for g in v})
TARGETS = PIN1 + CELL_CYCLE + DNA_REPAIR
KEEP = sorted(set(TARGETS) | set(MARKER_GENES))
