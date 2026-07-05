#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch

# =============================
# STYLE
# =============================
mpl.rcParams["hatch.linewidth"] = 1.5
plt.rcParams.update({
    "font.size": 18,
    "axes.linewidth": 1.2
})

# ============================================================
# PANEL A: S-hat mode values
# ============================================================
mode_data = [
    ("atp9", -6.77),
    ("cob", -5.91),
    ("cox1", -6.40),
    ("cox2", -4.99),
    ("nadh1-a", -6.82),
    ("nadh1-b", -5.85),
    ("nadh10", -6.95),
    ("nadh2", -4.65),
    ("nadh3", -4.88),
    ("nadh4", -5.01),
    ("nadh4L", -4.69),
    ("nadh5", -4.52),
    ("nadh6", -5.09),
    ("nadh7", -5.64),
    ("nadh9", -4.01),
    ("rpL14", -6.04),
    ("rpL16", -5.24),
    ("rpL2", -4.77),
    ("rpL32", -3.73),
    ("rpL6", -5.21),
    ("rpS12", -5.16),
    ("rpS13", -5.10),
    ("rpS14", -3.99),
    ("rpS19", -4.03),
    ("rpS3", -4.17),
    ("yejR", -3.96),
    ("ymf56", -5.81),
    ("ymf57", -5.24),
    ("ymf59", -3.54),
    ("ymf61", -3.91),
    ("ymf63", -4.76),
    ("ymf64", -3.25),
    ("ymf65", -4.61),
    ("ymf66", -5.03),
    ("ymf67", -4.15),
    ("ymf68", -5.80),
    ("ymf76", -3.06),
    ("ymf78", -4.64),
    ("ymf79", -4.86),
    ("ymf80", -4.59),
    ("ymf83", -3.76),
    ("ymf84", -3.46),
]

df = pd.DataFrame(mode_data, columns=["gene", "S_mode"])

# =============================
# Gene classes
# =============================
ETC = {
    "atp9", "cob", "cox1", "cox2",
    "nadh1-a", "nadh1-b", "nadh2", "nadh3", "nadh4",
    "nadh4L", "nadh5", "nadh6", "nadh7", "nadh9", "nadh10",
    "yejR"
}

RP = {
    "rpL2", "rpL6", "rpL14", "rpL16", "rpL32",
    "rpS3", "rpS12", "rpS13", "rpS14", "rpS19"
}

def gene_class(g):
    if g.startswith("ymf"):
        return "YMF"
    if g in RP:
        return "RP"
    if g in ETC:
        return "ETC"
    return "Other"

df["class"] = df["gene"].apply(gene_class)

cols = {
    "ETC": "orange",
    "RP": "green",
    "YMF": "steelblue",
    "Other": "lightcoral"
}

labels = {
    "ETC": "Electron transport chain",
    "RP": "Ribosomal proteins",
    "YMF": "Ciliate-specific",
    "Other": "Other"
}

order = ["ETC", "RP", "YMF", "Other"]
df["Class"] = pd.Categorical(df["class"], categories=order, ordered=True)

# Sort by functional group and S_mode
df = df.sort_values(["Class", "S_mode"]).reset_index(drop=True)

# ============================================================
# BARPLOT DATA
# Values are medians with asymmetric IQR error bars
# ============================================================
categories = ["ETC", "LSU", "SSU"]

# -----------------------------
# pN/pS medians
# UPDATED mitochondrial values:
# ETC = 0.0605 -0.0138 / +0.0321
# LSU = 0.0357 -0.0131 / +0.0044
# SSU = 0.0766 -0.0047 / +0.0209
# -----------------------------
mito_pnps = [0.0605, 0.0357, 0.0766]
nuc_pnps  = [0.0666, 0.0709, 0.0596]

# pN/pS error bars
mito_pnps_err_low  = [0.0138, 0.0131, 0.0047]
mito_pnps_err_high = [0.0321, 0.0044, 0.0209]

nuc_pnps_err_low   = [0.0259, 0.0216, 0.0143]
nuc_pnps_err_high  = [0.0286, 0.0210, 0.0263]

# -----------------------------
# dN/dS / omega medians
# -----------------------------
mito_dnds = [0.0335, 0.0294, 0.0654]
nuc_dnds  = [0.0355, 0.0300, 0.0339]

# dN/dS error bars
mito_dnds_err_low  = [0.0177, 0.0002, 0.0330]
mito_dnds_err_high = [0.0113, 0.0121, 0.0073]

nuc_dnds_err_low   = [0.0155, 0.0052, 0.0110]
nuc_dnds_err_high  = [0.0106, 0.0074, 0.0119]

cat_colors = {
    "ETC": "orange",
    "LSU": "#4caf50",
    "SSU": "#8bc34a"
}

x = np.arange(len(categories))
width = 0.34

error_kw = dict(
    ecolor="black",
    elinewidth=1.5,
    capsize=5,
    capthick=1.5
)

# ============================================================
# FIGURE
# ============================================================
fig, axes = plt.subplots(
    1, 3,
    figsize=(20, 7),
    gridspec_kw={"width_ratios": [2.3, 1.2, 1.2]}
)

# ============================================================
# PANEL A: S-hat
# ============================================================
ax = axes[0]

bar_colors = [cols[c] for c in df["class"]]

ax.bar(
    df["gene"],
    df["S_mode"],
    color=bar_colors,
    edgecolor="black",
    linewidth=0.8,
    width=0.8
)

ax.axhline(0, linestyle="--", color="black", linewidth=1)
ax.set_ylim(-7.5, 0)
ax.set_yticks([0, -1, -2, -3, -4, -5, -6, -7])

ax.set_ylabel("Modal scaled selection coefficient ($\\hat{S}$)", fontsize=20, color="black")
ax.set_xlabel("Genes", fontsize=24, color="black")

ax.tick_params(axis="x", rotation=90, labelsize=12, width=1.2)
ax.tick_params(axis="y", labelsize=16, width=1.2)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)

legend_patches = [
    Patch(facecolor=cols["ETC"], edgecolor="black", label=labels["ETC"]),
    Patch(facecolor=cols["RP"], edgecolor="black", label=labels["RP"]),
    Patch(facecolor=cols["YMF"], edgecolor="black", label=labels["YMF"]),
    Patch(facecolor=cols["Other"], edgecolor="black", label=labels["Other"]),
]

ax.legend(
    handles=legend_patches,
    loc="lower right",
    frameon=True,
    edgecolor="black",
    fontsize=11,
    title_fontsize=12,
    borderpad=0.35,
    labelspacing=0.25,
    handlelength=1.1,
    handletextpad=0.4
)

# ============================================================
# PANEL B: πN/πS with error bars
# ============================================================
ax = axes[1]

for i, cat in enumerate(categories):

    # Nuclear
    ax.bar(
        x[i] - width / 2,
        nuc_pnps[i],
        width,
        yerr=np.array([[nuc_pnps_err_low[i]], [nuc_pnps_err_high[i]]]),
        error_kw=error_kw,
        color=cat_colors[cat],
        edgecolor="black",
        linewidth=1.8,
        hatch="///"
    )

    # Mitochondrial
    ax.bar(
        x[i] + width / 2,
        mito_pnps[i],
        width,
        yerr=np.array([[mito_pnps_err_low[i]], [mito_pnps_err_high[i]]]),
        error_kw=error_kw,
        color=cat_colors[cat],
        edgecolor="black",
        linewidth=1.8
    )

ax.set_xticks(x)
ax.set_xticklabels(["ETC", "Ribo-LSU", "Ribo-SSU"], fontsize=18)
ax.set_ylabel("Median πN/πS", fontsize=24, color="black")
ax.set_ylim(0, 0.17)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

legend_elements = [
    Patch(facecolor="white", edgecolor="black", hatch="///", label="Nuclear"),
    Patch(facecolor="white", edgecolor="black", label="Mitochondrial")
]

ax.legend(
    handles=legend_elements,
    loc="upper left",
    frameon=True,
    edgecolor="black",
    fontsize=16
)

# ============================================================
# PANEL C: dN/dS with error bars
# ============================================================
ax = axes[2]

for i, cat in enumerate(categories):

    # Nuclear
    ax.bar(
        x[i] - width / 2,
        nuc_dnds[i],
        width,
        yerr=np.array([[nuc_dnds_err_low[i]], [nuc_dnds_err_high[i]]]),
        error_kw=error_kw,
        color=cat_colors[cat],
        edgecolor="black",
        linewidth=1.8,
        hatch="///"
    )

    # Mitochondrial
    ax.bar(
        x[i] + width / 2,
        mito_dnds[i],
        width,
        yerr=np.array([[mito_dnds_err_low[i]], [mito_dnds_err_high[i]]]),
        error_kw=error_kw,
        color=cat_colors[cat],
        edgecolor="black",
        linewidth=1.8
    )

ax.set_xticks(x)
ax.set_xticklabels(["ETC", "Ribo-LSU", "Ribo-SSU"], fontsize=18)
ax.set_ylabel("Median dN/dS", fontsize=24, color="black")
ax.set_ylim(0, 0.11)

for spine in ax.spines.values():
    spine.set_linewidth(1.2)

ax.legend(
    handles=legend_elements,
    loc="upper left",
    frameon=True,
    edgecolor="black",
    fontsize=16
)

# ============================================================
# FINAL
# ============================================================
plt.tight_layout()
plt.savefig("FINAL_clean_shared_legend_with_updated_mito_pnps_errorbars.png", dpi=300, bbox_inches="tight")
plt.show()