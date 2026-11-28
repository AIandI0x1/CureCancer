# Pre-specified (written before group-level results were computed) gene sets for the H2 test.
PIN1 = ["PIN1"]
CELL_CYCLE = ["CCND1", "CDK4", "CDK6", "CCNE1", "CDK2", "CCNA2", "CCNB1", "CDK1", "E2F1", "PCNA", "MCM2", "MKI67"]
DNA_REPAIR = ["ATM", "BRCA1", "RAD51", "XRCC1", "OGG1", "ERCC1", "PARP1", "MLH1", "MSH2", "XRCC5"]
# Descriptive composition markers (not part of H2; used to judge confounding by cell loss / gliosis)
NEURON = ["RBFOX3", "SNAP25", "SYT1"]
ASTRO = ["GFAP", "AQP4"]
MICROGLIA = ["AIF1", "CD68"]
ALL = PIN1 + CELL_CYCLE + DNA_REPAIR + NEURON + ASTRO + MICROGLIA
