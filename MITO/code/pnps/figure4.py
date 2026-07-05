import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
from matplotlib.lines import Line2D

# =============================
# PANEL A DATA (species-level)
# =============================
species_df = pd.DataFrame({
    "Species": [
        "biaurelia", "bursaria", "caudatum", "fokini",
        "primaurelia", "sexaurelia",
        "tetraurelia", "triaurelia", "tribursaria"
    ],
    "Mean_pN_pS": [
        0.26712627357559704,
        0.12799748995058846,
        0.10126902167113042,
        0.28768580992814785,
        0.06769503431414577,
        0.04808782536776819,
        0.12192411769988735,
        0.11489207795631895,
        0.14980295120543602
    ],
    "Mean_pS": [
        0.012181383869462506,
        0.02673186499325111,
        0.18390347246633298,
        0.07486924999279128,
        0.1518607012409594,
        0.2373432630233386,
        0.02037326153349239,
        0.024748386225031917,
        0.017854419415682266
    ]
})

# REMOVE multi
species_df = species_df[species_df["Species"] != "multi"].reset_index(drop=True)

# =============================
# COLOR MAP
# =============================
color_map = {
    "tribursaria": "#88c96b",
    "bursaria": "#cc0000",
    "multi": "#9874f4",
    "fokini": "#ffe599",
    "caudatum": "#e69138",
    "sexaurelia": "#559738",
    "tetraurelia": "#0b5394",
    "tredecaurelia": "#f1c232",
    "quadecaurelia": "#ccd9c7",
    "Psch": "#aaaaaa",
    "pquin": "#d5ed8a",
    "biaurelia": "#ee79ac",
    "triaurelia": "#d5a6bd",
    "primaurelia": "#6fa8dc"
}

species_df["Color"] = species_df["Species"].map(color_map)

rho, pval = spearmanr(species_df["Mean_pS"], species_df["Mean_pN_pS"])

# =============================
# INPUT FILES
# =============================
FIN_PNPS = "pnps_mean_summary_all_species_no_missing_ratio.csv"
FIN_DNDS = "M8_beta_p_q_omega.csv"

# =============================
# GENE SETS
# =============================
lgnames = {
    "ETC": ["nadh1-a","nadh1-b","nadh2","nadh3","nadh4","nadh4L","nadh5",
            "nadh6","nadh7","nadh9","nadh10","cob","cox1","cox2","atp9","yejR"],
    "RP":  ["rpL2","rpL6","rpL14","rpL16","rpL32",
            "rpS3","rpS12","rpS13","rpS14","rpS19"],
    "YMF": ["ymf56","ymf57","ymf59","ymf61","ymf63","ymf64","ymf65","ymf66",
            "ymf67","ymf68","ymf76","ymf78","ymf79","ymf80","ymf83","ymf84"]
}
gene2class = {g: cls for cls, genes in lgnames.items() for g in genes}

# =============================
# COLORS
# =============================
cols = {
    "ETC": "orange",
    "RP": "green",
    "YMF": "steelblue"
}

legend_labels = {
    "ETC": "Electron transport chain",
    "RP":  "Ribosomal proteins",
    "YMF": "Ciliate-specific"
}

# =============================
# pN/pS
# =============================
pnps = pd.read_csv(FIN_PNPS)
pnps["Gene"] = pnps["Gene"].astype(str).str.strip()
pnps["pN_over_pS"] = pd.to_numeric(pnps["pN_over_pS"], errors="coerce")
pnps = pnps.replace([np.inf, -np.inf], np.nan).dropna()
pnps = pnps[pnps["pN_over_pS"] > 0]

pnps_gene = (
    pnps[["Gene", "pN_over_pS"]]
    .groupby("Gene", as_index=False)
    .agg(pN_over_pS=("pN_over_pS", "median"))
)

# =============================
# dN/dS
# =============================
dnds = pd.read_csv(FIN_DNDS).rename(columns={"gene": "Gene"})
dnds["Gene"] = dnds["Gene"].astype(str).str.strip()
dnds["omega_mean"] = pd.to_numeric(dnds["omega_mean"], errors="coerce")
dnds = dnds.replace([np.inf, -np.inf], np.nan).dropna()
dnds = dnds[dnds["omega_mean"] > 0]

dnds_gene = (
    dnds[["Gene", "omega_mean"]]
    .groupby("Gene", as_index=False)
    .agg(omega_mean=("omega_mean", "mean"))
)

# =============================
# Merge + assign classes
# =============================
dat = pd.merge(pnps_gene, dnds_gene, on="Gene", how="inner")
dat["Class"] = dat["Gene"].map(gene2class)

# REMOVE unannotated genes
dat = dat.dropna(subset=["Class"])

group_order = ["ETC", "RP", "YMF"]

# =============================
# COMBINED FIGURE
# =============================
fig = plt.figure(figsize=(16, 6.5))
outer = GridSpec(1, 2, figure=fig, width_ratios=[1.2, 1.6], wspace=0.18)

# =============================
# LEFT PANEL
# =============================
ax1 = fig.add_subplot(outer[0, 0])

for _, row in species_df.iterrows():
    ax1.scatter(
        row["Mean_pS"], row["Mean_pN_pS"],
        s=100,
        color=row["Color"],
        edgecolor="black",
        linewidth=0.7,
        zorder=3
    )

coef = np.polyfit(species_df["Mean_pS"], species_df["Mean_pN_pS"], 1)
xfit = np.linspace(species_df["Mean_pS"].min(), species_df["Mean_pS"].max(), 200)
ax1.plot(xfit, np.polyval(coef, xfit), "--", color="black", linewidth=1.5, zorder=2)

ax1.text(
    0.98, 0.95,
    f"Spearman ρ = {rho:.2f}\nP = {pval:.3f}",
    transform=ax1.transAxes,
    ha="right", va="top",
    fontsize=11,
    color="black"
)

ax1.set_xlabel(r"πS", fontsize=20)
ax1.set_ylabel("πN/πS", fontsize=20,  color="black")
ax1.tick_params(axis="both", labelsize=12, colors="black")

for spine in ax1.spines.values():
    spine.set_linewidth(1.2)
    spine.set_color("black")

species_handles = [
    Line2D(
        [0], [0],
        marker='o',
        color='w',
        markerfacecolor=color_map[s],
        markeredgecolor='black',
        markersize=8,
        label=s
    )
    for s in species_df["Species"]
]

ax1.legend(
    handles=species_handles,
    loc="lower right",
    bbox_to_anchor=(1, 0.35),
    fontsize=12,
    frameon=True,
    edgecolor="black",
    facecolor="white"
)

# =============================
# RIGHT PANEL
# =============================
right = outer[0, 1].subgridspec(4, 4, hspace=0.05, wspace=0.05)

ax_top = fig.add_subplot(right[0, 0:3])
ax_scatter = fig.add_subplot(right[1:4, 0:3], sharex=ax_top)
ax_right = fig.add_subplot(right[1:4, 3], sharey=ax_scatter)
ax_scatter.tick_params(axis="both", labelsize=12)

# Scatter
for cls in group_order:
    sub = dat[dat["Class"] == cls]
    if len(sub) == 0:
        continue
    ax_scatter.scatter(
        sub["omega_mean"], sub["pN_over_pS"],
        color=cols[cls],
        label=legend_labels[cls],
        alpha=0.85,
        s=65
    )

ax_scatter.set_xlabel("dN/dS", fontsize=20, color="black")
ax_scatter.set_ylabel("πN/πS", fontsize=20,  color="black")
ax_scatter.legend(
    loc="upper right",
    fontsize=12,
    frameon=True,
    fancybox=False,
    edgecolor="black",
    framealpha=1.0
)

# =============================
# Top boxplots (dN/dS)
# =============================
data_top = [dat.loc[dat["Class"] == g, "omega_mean"].values for g in group_order]
bp_top = ax_top.boxplot(data_top, vert=False, patch_artist=True, showfliers=True)

for patch, g in zip(bp_top["boxes"], group_order):
    patch.set_facecolor(cols[g])
    patch.set_alpha(0.7)

for part in ["whiskers", "caps", "medians"]:
    for obj in bp_top[part]:
        obj.set_color("black")

ax_top.set_yticks([])

# =============================
# Right boxplots (pN/pS)
# =============================
bp_right = ax_right.boxplot(
    [dat[dat["Class"] == g]["pN_over_pS"] for g in group_order],
    vert=True,
    patch_artist=True,
    showfliers=True
)

for patch, g in zip(bp_right["boxes"], group_order):
    patch.set_facecolor(cols[g])
    patch.set_edgecolor("black")
    patch.set_linewidth(1.0)

for key in ["whiskers", "caps", "medians"]:
    for item in bp_right[key]:
        item.set_color("black")
        item.set_linewidth(1.0)

for spine in ax_right.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.2)
    spine.set_color("black")

ax_right.set_xticks([])
ax_right.tick_params(axis="y", labelleft=False, left=False)

plt.tight_layout()
plt.savefig("FINAL_combined_correct_median.png", dpi=300, bbox_inches="tight")
plt.show()