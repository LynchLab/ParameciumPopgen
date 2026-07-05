import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------
# Load data
# -----------------------------
file_path = "combined_gene_coords_table.csv"  # adjust path if needed
df = pd.read_csv(file_path)

# Clean names
df["sample"] = df["sample"].astype(str).str.strip()
df["gene_name"] = df["gene_name"].astype(str).str.strip()

# Remove tRNA / rRNA / misc RNA genes
df = df[~df["gene_name"].str.lower().str.contains(r"trna|rrna|rna", regex=True)].copy()

# Normalize each gene by its median length
df["normalized_length"] = df.groupby("gene_name")["length"].transform(
    lambda x: x / x.median()
)

# -----------------------------
# Define species groups
# -----------------------------
aurelia = [
    "tetr", "tria", "pent", "sexa", "octa", "dode", "tred",
    "biau", "deca", "jenn", "nova", "prima", "psch", "quad", "sonn"
]
others = ["caud", "mult", "foki", "burs", "tbur", "calk"]

# Map sample codes to proper species names
sample_name_map = {
    "caud": "P. caudatum",
    "burs": "P. bursaria",
    "tbur": "P. tribursaria",
    "mult": "P. multimicronucleatum",
    "foki": "P. fokini",
    "calk": "P. calkinski",  # spelling kept as-is from sample
}

df["species"] = df["sample"].map(sample_name_map).fillna(df["sample"])

aurelia_df = df[df["sample"].isin(aurelia)].copy()
others_df  = df[df["sample"].isin(others)].copy()

# -----------------------------
# Outgroup colors
# -----------------------------
outgroup_colors = {
    "P. tribursaria": "#88c96b",
    "P. bursaria": "#cc0000",
    "P. multimicronucleatum": "#9874f4",
    "P. fokini": "#ffe599",
    "P. caudatum": "#e69138",
    # add calkinsi if desired
    # "P. calkinski": "#000000",
}

# -----------------------------
# Plot function
# -----------------------------
def plot_gene_length_differences(gene_list, out_png, title, ylim):
    """
    Aurelia: gray boxplots + gray jitter points
    Outgroups: colored jitter points using outgroup_colors
    """
    genes = gene_list if gene_list is not None else sorted(df["gene_name"].unique())
    x = np.arange(len(genes))

    fig_w = max(10, 1.2 * len(genes))
    plt.figure(figsize=(fig_w, 6))
    ax = plt.gca()

    # --- Aurelia boxplots (gray) ---
    box_data = [
        aurelia_df.loc[aurelia_df["gene_name"] == g, "normalized_length"]
        .dropna()
        .values
        for g in genes
    ]
    ax.boxplot(
        box_data,
        positions=x,
        widths=0.6,
        showfliers=False,
        patch_artist=True,
        boxprops=dict(facecolor="lightgray", edgecolor="black", linewidth=0.8),
        medianprops=dict(color="black", linewidth=1.1),
        whiskerprops=dict(color="black", linewidth=0.8),
        capprops=dict(color="black", linewidth=0.8),
    )

    # RNG for reproducible jitter
    rng = np.random.default_rng(3)

    # --- Aurelia jitter points (gray) ---
    for i, g in enumerate(genes):
        vals = aurelia_df.loc[
            aurelia_df["gene_name"] == g, "normalized_length"
        ].dropna().values
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            color="gray",
            s=18,
            alpha=0.7,
            zorder=3,
        )

    # --- Outgroup colored points ---
    for sp, col in outgroup_colors.items():
        sub_sp = others_df[others_df["species"] == sp]
        for i, g in enumerate(genes):
            vals = sub_sp.loc[
                sub_sp["gene_name"] == g, "normalized_length"
            ].dropna().values
            if len(vals) == 0:
                continue
            jitter = rng.uniform(-0.15, 0.15, size=len(vals))
            ax.scatter(
                np.full(len(vals), i) + jitter,
                vals,
                color=col,
                s=40,
                edgecolor="black",
                linewidth=0.4,
                zorder=4,
            )

    # Reference line at median
    ax.axhline(1.0, linestyle="--", linewidth=1, color="black", alpha=0.6)

    # Labels & title
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=90)
    ax.set_ylabel("Relative gene length (length / median per gene)")
    ax.set_title(title, weight="bold")
    ax.set_ylim(ylim)

    # Legend
    handles = [Line2D([0], [0], marker="o", color="gray", linestyle="None", label="Aurelia")]
    for sp, col in outgroup_colors.items():
        handles.append(Line2D([0], [0], marker="o", color=col, linestyle="None", label=sp))

    ax.legend(
        handles=handles,
        title="Species",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=600)
    plt.close()
    print(f"✅ Saved: {out_png}")

# -----------------------------
# FIGURE 1: Selected genes (tighter scale, starts at 0.8)
# -----------------------------
selected_genes = ["cox2", "cob", "nadh1-b", "nadh10", "nadh5", "rpS19", "ymf66"]

plot_gene_length_differences(
    gene_list=selected_genes,
    out_png="gene_length_differences_selected_genes.png",
    title="Gene length differences",
    ylim=(0.9, 1.7),
)

# -----------------------------
# FIGURE 2: All genes (broader scale)
# -----------------------------
plot_gene_length_differences(
    gene_list=None,
    out_png="gene_length_differences_all_genes.png",
    title="Gene length differences (all mitochondrial genes)",
    ylim=(0.6, 2.0),
)
