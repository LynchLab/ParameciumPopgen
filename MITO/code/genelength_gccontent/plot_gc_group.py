import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================================================
# Median GC (%) values
# ============================================================
gc_data = [
    {"Species": "P. biaurelia",           "Protein-coding": 32.83, "rRNA": 35.97, "tRNA": 45.76, "Intergenic": 34.19},
    {"Species": "P. decaurelia",          "Protein-coding": 38.91, "rRNA": 36.37, "tRNA": 43.73, "Intergenic": 40.84},
    {"Species": "P. dodecaurelia",        "Protein-coding": 40.30, "rRNA": 36.52, "tRNA": 43.73, "Intergenic": 40.58},
    {"Species": "P. jenningsi",           "Protein-coding": 39.76, "rRNA": 36.80, "tRNA": 44.67, "Intergenic": 39.87},
    {"Species": "P. novaurelia",          "Protein-coding": 40.08, "rRNA": 35.97, "tRNA": 43.73, "Intergenic": 40.14},
    {"Species": "P. octaurelia",          "Protein-coding": 39.69, "rRNA": 35.77, "tRNA": 43.73, "Intergenic": 39.24},
    {"Species": "P. pentaurelia",         "Protein-coding": 39.63, "rRNA": 35.88, "tRNA": 43.73, "Intergenic": 39.85},
    {"Species": "P. quindecaurelia",      "Protein-coding": 33.34, "rRNA": 35.60, "tRNA": 44.44, "Intergenic": 29.07},
    {"Species": "P. primaurelia",         "Protein-coding": 35.85, "rRNA": 36.16, "tRNA": 44.28, "Intergenic": 39.31},
    {"Species": "P. schewiakoffi",        "Protein-coding": 37.61, "rRNA": 36.39, "tRNA": 43.73, "Intergenic": 37.49},
    {"Species": "P. quadecaurelia",       "Protein-coding": 38.71, "rRNA": 36.09, "tRNA": 43.73, "Intergenic": 39.67},
    {"Species": "P. sexaurelia",          "Protein-coding": 37.68, "rRNA": 35.99, "tRNA": 44.74, "Intergenic": 38.14},
    {"Species": "P. sonneborni",          "Protein-coding": 40.02, "rRNA": 36.43, "tRNA": 43.73, "Intergenic": 40.76},
    {"Species": "P. tetraurelia",         "Protein-coding": 41.54, "rRNA": 36.46, "tRNA": 44.67, "Intergenic": 41.23},
    {"Species": "P. triaurelia",          "Protein-coding": 35.65, "rRNA": 37.17, "tRNA": 44.46, "Intergenic": 34.42},
    {"Species": "P. caudatum",            "Protein-coding": 20.53, "rRNA": 32.41, "tRNA": 43.57, "Intergenic": 22.15},
    {"Species": "P. multimicronucleatum", "Protein-coding": 19.15, "rRNA": 30.78, "tRNA": 45.70, "Intergenic": 20.17},
    {"Species": "P. fokini",              "Protein-coding": 20.80, "rRNA": 31.47, "tRNA": 46.19, "Intergenic": 22.97},
    {"Species": "P. bursaria",            "Protein-coding": 14.70, "rRNA": 26.55, "tRNA": 39.65, "Intergenic": 14.89},
    {"Species": "P. tria-bursaria",       "Protein-coding": 14.55, "rRNA": 25.72, "tRNA": 40.56, "Intergenic":  9.65},
]

# ============================================================
# Groups and colors
# ============================================================
outgroups = [
    "P. caudatum",
    "P. multimicronucleatum",
    "P. fokini",
    "P. bursaria",
    "P. tria-bursaria",
]

outgroup_palette = {
    "P. caudatum": "#e69138",
    "P. multimicronucleatum": "#9874f4",
    "P. fokini": "#ffe599",
    "P. bursaria": "#cc0000",
    "P. tria-bursaria": "#88c96b",
}

# ============================================================
# Build dataframe
# ============================================================
df = pd.DataFrame(gc_data)
df["Group"] = df["Species"].apply(lambda s: "Outgroup" if s in outgroups else "Aurelia")

df_long = df.melt(
    id_vars=["Species", "Group"],
    value_vars=["Protein-coding", "Intergenic", "rRNA", "tRNA"],
    var_name="Region",
    value_name="GC"
).dropna()

region_order = ["Protein-coding", "Intergenic", "rRNA", "tRNA"]

# ============================================================
# Plot
# ============================================================
plt.figure(figsize=(9.5, 6))
ax = plt.gca()

# Aurelia boxplots
sns.boxplot(
    data=df_long[df_long["Group"] == "Aurelia"],
    x="Region", y="GC",
    order=region_order,
    color="lightgray",
    showfliers=False,
    width=0.65,
    zorder=1
)

# Aurelia points
sns.stripplot(
    data=df_long[df_long["Group"] == "Aurelia"],
    x="Region", y="GC",
    order=region_order,
    color="gray",
    size=4,
    jitter=0.12,
    alpha=0.7,
    zorder=2
)

# Outgroup points
sns.stripplot(
    data=df_long[df_long["Group"] == "Outgroup"],
    x="Region", y="GC",
    order=region_order,
    hue="Species",
    hue_order=outgroups,
    palette=outgroup_palette,
    size=7,
    jitter=0.12,
    dodge=True,
    edgecolor="black",
    linewidth=0.4,
    zorder=3
)

# ============================================================
# Styling
# ============================================================
ax.set_title(
    "Mitochondrial GC content variation",
    fontsize=13,
    weight="bold"
)
ax.set_xlabel("Genomic regions")
ax.set_ylabel("Median GC content (%)")
ax.set_ylim(0, 50)

# Legend
handles, labels = ax.get_legend_handles_labels()
aurelia_handle = Line2D([0], [0], marker="o", color="gray", linestyle="None", label="Aurelia")

ax.legend(
    [aurelia_handle] + handles,
    ["Aurelia"] + labels,
    title="Species",
    bbox_to_anchor=(1.02, 1),
    loc="upper left"
)

plt.tight_layout()
plt.savefig("median_gc_compartments_aurelia_vs_outgroups.png", dpi=600)
plt.close()

print("✅ Saved: median_gc_compartments_aurelia_vs_outgroups.png")
