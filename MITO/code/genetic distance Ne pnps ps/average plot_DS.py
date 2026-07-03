# FigS1_piS_barplot_Ne.py
# Synonymous diversity (πS) as a proxy for effective population size (Ne)
# Species ordered from low → high πS

import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv("poster images/paper figures/data/genetic distance Ne group pnps ps/average_dS_per_species.csv")

# Column names in YOUR file
species_col = "species"
y_col = "mean_dS"   # this is what your CSV actually contains

# -----------------------------
# Proper species labels (PNAS style)
# -----------------------------
species_labels = {
    "biaurelia": "P. biaurelia",
    "bursaria": "P. bursaria",
    "caudatum": "P. caudatum",
    "fokini": "P. fokini",
    "multi": "P. multimicronucleatum",
    "pquin": "P. quindecaurelia",
    "primaurelia": "P. primaurelia",
    "sexaurelia": "P. sexaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "tribursaria": "P. tribursaria"
}

# -----------------------------
# Your exact color map
# -----------------------------
color_map = {
    "tribursaria": "#88c96b",
    "bursaria": "#cc0000",
    "multi": "#9874f4",
    "fokini": "#ffe599",
    "caudatum": "#e69138",
    "sexaurelia": "#559738",
    "tetraurelia": "#0b5394",
    "pquin": "#d5ed8a",
    "biaurelia": "#ee79ac",
    "triaurelia": "#d5a6bd",
    "primaurelia": "#6fa8dc"
}

# -----------------------------
# Sort ascending (low → high), which corresponds to increasing Ne under pi_S ~ 2Neμ
# -----------------------------
df = df.sort_values(by=y_col, ascending=True).reset_index(drop=True)

# Build plot arrays
x_labels = [species_labels.get(s, s) for s in df[species_col]]
y_values = df[y_col].values
bar_colors = [color_map.get(s, "#999999") for s in df[species_col]]

# -----------------------------
# Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(7.2, 4.8))

ax.bar(
    x_labels,
    y_values,
    color=bar_colors,
    edgecolor="black",
    linewidth=0.6
)

# Axis labels + title (Ne framing)
ax.set_ylabel(r"Average synonymous diversity", fontsize=12)
ax.set_title(
    "Effective population size across Paramecium species",
    fontsize=13
)

# PNAS-clean look
ax.tick_params(axis="x", rotation=80, labelsize=10)
ax.tick_params(axis="y", labelsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

# -----------------------------
# Save
# -----------------------------

# -----------------------------
# Save
# -----------------------------
plt.savefig("poster images/paper figures/images/FigS1_Ne_inferred_from_piS.png", dpi=300)
plt.savefig("poster images/paper figures/images/FigS1_Ne_inferred_from_piS.pdf")
plt.show()
