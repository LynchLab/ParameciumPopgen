import pandas as pd

# read files
etc = pd.read_csv("etc/ancestral_dS_icode5/ETC_avg_nuclear_genes.csv")
ribo = pd.read_csv("ribo/ancestral_dS_icode5/Ribo_avg_nuclear_genes.csv")

# fix species names in ribo file
species_map = {
    "PBIA": "biau",
    "PBUR": "burs",
    "PCAU": "caud",
    "PMMN": "mult",
    "PPRIM": "prim",
    "PSEX": "sexa",
    "PTET": "tetr"
}

ribo["species"] = ribo["species"].map(species_map)

# merge
merged = etc.merge(ribo, on="species", suffixes=("_etc", "_ribo"))

# weighted average using mean dS and gene counts
merged["weighted_dS"] = (
    merged["mean_etc"] * merged["count_etc"] +
    merged["mean_ribo"] * merged["count_ribo"]
) / (merged["count_etc"] + merged["count_ribo"])

result = merged[["species", "weighted_dS"]].sort_values("species")

# save
result.to_csv("weighted_avg_nuclear_ETC_RIBO_dS.csv", index=False)

print(result.to_string(index=False))
