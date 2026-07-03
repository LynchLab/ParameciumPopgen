import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# ===== INPUT =====
csv_path = "phi_full_summary.csv"   # adjust path if needed

# ===== OUTPUT (PNG) =====
output_png = "Recomb_tests_heatmaps_gene_vs_species_with_grid.png"

# ===== LOAD DATA =====
df = pd.read_csv(csv_path)

df["Species"] = df["Species"].astype(str).str.strip()
df["Gene"] = df["Gene"].astype(str).str.strip()

# ===== Species order you want (all species, exactly once) =====
species_order_full = [
    "sexaurelia",
    "fokini",
    "bursaria",
    "caudatum",
    "tribursaria",
    "multi",
    "biaurelia",
    "triaurelia",
    "tetraurelia",
    "quadecaurelia",
    "tredecaurelia",
    "primaurelia",
    "pquin",
]

# ===== TESTS =====
tests = [
    ("NSS", "NSS"),
    ("MaxChi2", "MaxChi2"),
    ("PHI_perm", "PHI_perm"),
    ("PHI_norm", "PHI_norm"),
]

# Convert test columns to numeric
for col, _ in tests:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# ===== Consistent gene order =====
genes_sorted = sorted(df["Gene"].unique())

# ===== Color map =====
# 0 = NA         → grey
# 1 = non-sig    → red
# 2 = significant→ green
cmap = ListedColormap(["grey", "red", "green"])

def classify_matrix(test_col):
    """
    Build Gene×Species matrix with:
      2 = significant (p < 0.05)
      1 = non-significant
      0 = NA
    """
    pivot = df.pivot_table(
        index="Gene",
        columns="Species",
        values=test_col,
        aggfunc="min",  # min p-value per Gene×Species
    )

    # reorder genes
    pivot = pivot.reindex(genes_sorted)

    # force all species in species_order_full (grey if missing)
    pivot = pivot.reindex(columns=species_order_full)

    species_order = species_order_full
    n_genes = len(genes_sorted)
    n_species = len(species_order)

    mat = np.zeros((n_genes, n_species), dtype=int)

    for i, gene in enumerate(genes_sorted):
        for j, sp in enumerate(species_order):
            val = pivot.loc[gene, sp]
            if pd.isna(val):
                mat[i, j] = 0
            elif val < 0.05:
                mat[i, j] = 2
            else:
                mat[i, j] = 1

    return mat, species_order

# ===== Create 4-panel heatmap figure =====
fig, axes = plt.subplots(1, 4, figsize=(28, 14), sharey=True)

for ax, (col, title) in zip(axes, tests):

    mat, species_order = classify_matrix(col)

    ax.imshow(mat, cmap=cmap, aspect="auto", vmin=0, vmax=2)

    n_genes = len(genes_sorted)
    n_species = len(species_order)

    ax.set_title(f"{title} (window = 100)", fontsize=16)

    ax.set_xticks(np.arange(n_species))
    ax.set_xticklabels([s[:4] for s in species_order], rotation=90, fontsize=10)

    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(genes_sorted, fontsize=10)

    ax.set_xlabel("Species (4 letters)")

    # ---- Black grid lines between all cells ----
    ax.set_xticks(np.arange(-0.5, n_species, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color="black", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

axes[0].set_ylabel("Gene")

# ===== Legend =====
legend_elements = [
    Patch(facecolor="green", edgecolor="black", label="Present (p < 0.05)"),
    Patch(facecolor="red",   edgecolor="black", label="Absent (non-sig)"),
    Patch(facecolor="grey",  edgecolor="black", label="Absent (NA / no data)"),
]
axes[0].legend(
    handles=legend_elements,
    title="Recombination signal",
    loc="upper left",
    bbox_to_anchor=(1.02, 1),
    fontsize="small",
    title_fontsize="small",
)

plt.tight_layout()
plt.savefig(output_png, dpi=600, bbox_inches="tight")
plt.close()

print(f"✅ Heatmap PNG with all species (grey if no data) saved as: {output_png}")
