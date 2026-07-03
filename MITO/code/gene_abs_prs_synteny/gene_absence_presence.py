#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import re

# ======== Input ========
INPUT_FILE = "gene_presence_absence_matrix_2.csv"
OUT_PNG    = "gene_presence_absence_heatmap.png"

# ======== Load presence/absence matrix ========
df = pd.read_csv(INPUT_FILE, index_col=0)

# --- Clean up genes and species ---
# Drop 'orfX' gene if present
if "orfX" in df.index:
    df = df.drop(index="orfX")

# Drop unwanted species: Tetrahymena, gigas, Pquin, and P. calkinski (case-insensitive)
drop_cols = [
    c for c in df.columns
    if ("tetrahymena" in c.lower())
    or ("gigas" in c.lower())
    or ("pquin" in c.lower())
    or ("calkinski" in c.lower())
]
if drop_cols:
    df = df.drop(columns=drop_cols)

# Force binary 0/1 values (robust to strings like "1", "0", "true", "false", etc.)
def to_bin(x):
    s = str(x).strip().lower()
    if s in {"1", "true", "t", "yes", "y"}:
        return 1
    elif s in {"0", "false", "f", "no", "n", ""}:
        return 0
    else:
        # Fallback: try numeric cast
        try:
            return 1 if float(s) != 0 else 0
        except Exception:
            return 0

df = df.applymap(to_bin)

# --- Robust species renaming (regex + longest-first to avoid substring collisions) ---
species_fullname_map = {
    "quindecaurelia": "P. quindecaurelia",
    "quadecaurelia": "P. quadecaurelia",
    "tredecaurelia": "P. tredecaurelia",
    "dodecaurelia": "P. dodecaurelia",
    "decaurelia": "P. decaurelia",
    "novaurelia": "P. novaurelia",
    "octaurelia": "P. octaurelia",
    "sexaurelia": "P. sexaurelia",
    "pentaurelia": "P. pentaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "biaurelia": "P. biaurelia",
    "primaurelia": "P. primaurelia",
    "sonneborni": "P. sonneborni",
    "jenningsi": "P. jenningsi",
    "Psche": "P. schewiakoffi",  # allows Psch, psch, etc.
    "caudatum": "P. caudatum",
    "multimicronucleatum": "P. multimicronucleatum",
    "multi": "P. multimicronucleatum",  # fallback short key
    "fokini": "P. fokini",
    "bursaria": "P. bursaria",
    "tribursaria": "P. tribursaria",
}

# Compile regexes with letter-boundaries (underscore/digits count as boundaries)
_compiled = [
    (k, re.compile(rf"(?i)(?:^|[^a-z]){re.escape(k)}(?:$|[^a-z])"))
    for k in sorted(species_fullname_map.keys(), key=len, reverse=True)
]

def rename_species(col: str) -> str:
    low = col.lower()
    for k, rx in _compiled:
        if rx.search(low):
            return species_fullname_map[k]
    return col

df.columns = [rename_species(c) for c in df.columns]

# Optional: detect duplicates after renaming (helps catch lingering collisions)
dup_counts = pd.Series(df.columns).value_counts()
dups = dup_counts[dup_counts > 1]
if not dups.empty:
    print("WARNING: Duplicate species columns after renaming:\n", dups)

# --- Desired species column ordering ---
aurelia_order = [
    "P. primaurelia", "P. biaurelia", "P. triaurelia", "P. tetraurelia",
    "P. pentaurelia", "P. sexaurelia", "P. octaurelia", "P. novaurelia",
    "P. decaurelia", "P. dodecaurelia", "P. tredecaurelia",
    "P. quadecaurelia", "P. quindecaurelia", "P. sonneborni",
    "P. jenningsi", "P. schewiakoffi"
]

group_order = [
    aurelia_order,
    ["P. caudatum"],
    ["P. multimicronucleatum"],
    ["P. fokini"],
    ["P. bursaria", "P. tribursaria"],
]

# Build final species order
ordered_cols = []
seen = set()
for group in group_order:
    for sp in group:
        if sp in df.columns and sp not in seen:
            ordered_cols.append(sp)
            seen.add(sp)

# Append any remaining unexpected species to the end
others = [c for c in df.columns if c not in seen]
final_cols = ordered_cols + others
df = df[final_cols]

# --- Desired gene order ---
gene_order = [
    "atp9",
    "cob",
    "cox1",
    "cox2",
    "nadh1-a",
    "nadh1-b",
    "nadh2",
    "nadh3",
    "nadh4",
    "nadh4L",
    "nadh5",
    "nadh6",
    "nadh7",
    "nadh9",
    "nadh10",
    "rpL2",
    "rpL6",
    "rpL14",
    "rpL16",
    "rpL32",
    "rpS3",
    "rpS12",
    "rpS13",
    "rpS14",
    "rpS19",
    "yejR",
    "ymf56",
    "ymf57",
    "ymf59",
    "ymf61",
    "ymf63",
    "ymf64",
    "ymf65",
    "ymf66",
    "ymf67",
    "ymf68",
    "ymf76",
    "ymf78",
    "ymf79",
    "ymf80",
    "ymf83",
    "ymf84",
    "orfL",
    "orfY",
    "orfZ",
]

# Keep only genes present in the dataframe, add extra genes at the end
ordered_genes = [g for g in gene_order if g in df.index]
extra_genes = [g for g in df.index if g not in ordered_genes]
final_rows = ordered_genes + extra_genes
df = df.loc[final_rows]

# ======== Plot ========
plt.figure(figsize=(16, 11))
# 0=absent (yellow), 1=present (steelblue)
cmap = sns.color_palette(["yellow", "steelblue"])

sns.heatmap(
    df,
    cmap=cmap,
    cbar=False,
    linewidths=1.5,
    linecolor="black",
    xticklabels=True,
    yticklabels=True
)

# Labels & Title
plt.xlabel("Species", fontsize=20,fontweight="bold")
plt.ylabel("Genes", fontsize=20,fontweight="bold")
plt.title("Presence and Absence of Genes in Paramecium Mitochondrial Genomes", fontsize=24, pad=20,fontweight="bold")

plt.xticks(rotation=45, ha="right", fontsize=14)
plt.yticks(fontsize=10)

# Custom Legend
legend_elements = [
    Patch(facecolor="yellow", edgecolor="black", label="Absent"),
    Patch(facecolor="steelblue", edgecolor="black", label="Present"),
]

plt.legend(
    handles=legend_elements,
    loc="upper right",
    bbox_to_anchor=(1.25, 1.0),
    title="Legend",
    frameon=True,
    fontsize=14,
    title_fontsize=14
)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.show()
