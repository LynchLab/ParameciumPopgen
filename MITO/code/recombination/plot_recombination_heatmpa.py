# plot_recombination_heatmap_all_tests_side_by_side.py
# ONE-PAGE SUPPLEMENT: 4 non-parametric recombination tests shown side-by-side
# - Panels: NSS, MaxChi2, PHI_perm, PHI_norm (auto-detected from column names)
# - Gene order: cox* → cob → nadh* → rpL* → rpS* → orf* → ymf* (then anything else)
# - Colors: green = present (p < 0.05), orange = absent (p >= 0.05), light grey = no data
# - Black borders between heatmap boxes + black-edged legend
#
# Fixes your error: converts p-value columns from strings to numeric before comparisons.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

# -----------------------------
# Input (EDIT THESE PATHS)
# -----------------------------
phi_path = r"poster images\paper figures\data\recombination\phi_full_summary.csv"
# Example Windows absolute path:
# phi_path = r"C:\Users\telka\OneDrive\Documents\omics\poster images\paper figures\data\recombination\phi_full_summary.csv"

# Output
out_png = r"poster images\paper figures\Images\FigS_recombination_tests_all_side_by_side.png"
out_pdf = r"poster images\paper figures\Images\FigS_recombination_tests_all_side_by_side.pdf"

# -----------------------------
# Load
# -----------------------------
phi = pd.read_csv(phi_path)

# -----------------------------
# Species labels (optional)
# -----------------------------
species_labels = {
    "biaurelia": "P. biaurelia",
    "bursaria": "P. bursaria",
    "caudatum": "P. caudatum",
    "multimicronucleatum": "P. multimicronucleatum",
    "multi": "P. multimicronucleatum",
    "fokini": "P. fokini",
    "quindecaurelia": "P. quindecaurelia",
    "pquin": "P. quindecaurelia",
    "primaurelia": "P. primaurelia",
    "sexaurelia": "P. sexaurelia",
    "tetraurelia": "P. tetraurelia",
    "triaurelia": "P. triaurelia",
    "tribursaria": "P. tribursaria"
}

# -----------------------------
# Helpers
# -----------------------------
def norm_species(s):
    s = str(s).strip()
    s = s.replace("Paramecium", "").strip()
    s = s.replace("P.", "").strip()
    s = s.replace(" ", "").replace(".", "")
    return s.lower()

def gene_group_rank(gene):
    g = str(gene).lower()
    if g.startswith("cox"):
        return (0, g)
    if g == "cob" or g.startswith("cytb"):
        return (1, g)
    if g.startswith("nadh"):
        return (2, g)
    if g.startswith("rpl"):
        return (3, g)
    if g.startswith("rps"):
        return (4, g)
    if g.startswith("orf"):
        return (5, g)
    if g.startswith("ymf"):
        return (6, g)
    return (7, g)

def find_col(df, must_contain, optional_contain=None):
    must = [x.lower() for x in (must_contain if isinstance(must_contain, (list, tuple)) else [must_contain])]
    opt  = [x.lower() for x in (optional_contain if optional_contain else [])]
    for c in df.columns:
        cl = c.lower()
        if all(m in cl for m in must) and all(o in cl for o in opt):
            return c
    return None

def to_numeric_clean(series):
    """
    Convert p-value-like strings to numeric safely.
    Handles: 'NA', 'na', '', 'None', 'nan', whitespace, and commas.
    Unparseable values -> NaN.
    """
    s = series.astype(str).str.strip()
    s = s.replace({
        "NA": np.nan, "na": np.nan, "NaN": np.nan, "nan": np.nan,
        "None": np.nan, "none": np.nan, "": np.nan
    })
    s = s.str.replace(",", "", regex=False)
    return pd.to_numeric(s, errors="coerce")

# -----------------------------
# Validate required columns
# -----------------------------
required = {"Species", "Gene"}
missing_req = [c for c in required if c not in phi.columns]
if missing_req:
    raise KeyError(f"Missing required columns in phi_full_summary.csv: {missing_req}\nColumns: {list(phi.columns)}")

# -----------------------------
# Standardize keys
# -----------------------------
phi2 = phi.copy()
phi2["species_key"] = phi2["Species"].apply(norm_species)
phi2["gene"] = phi2["Gene"].astype(str)

# -----------------------------
# Species order (preferred, then remaining)
# -----------------------------
preferred_species = [
    "biaurelia","bursaria","caudatum","fokini","primaurelia",
    "sexaurelia","tetraurelia","triaurelia","tribursaria","multimicronucleatum"
]
present_species = list(pd.unique(phi2["species_key"]))
species_order = [s for s in preferred_species if s in present_species] + \
                [s for s in present_species if s not in preferred_species]

# -----------------------------
# Gene order (as requested)
# -----------------------------
gene_order = sorted(pd.unique(phi2["gene"]), key=gene_group_rank)

# -----------------------------
# Detect the 4 p-value columns
# -----------------------------
col_nss     = find_col(phi2, "nss")
col_maxchi2 = find_col(phi2, ["maxchi"], ["2"]) or find_col(phi2, "maxchi2") or find_col(phi2, "maxchi")
col_phiperm = find_col(phi2, ["phi"], ["perm"]) or find_col(phi2, "phi_perm")
col_phinorm = find_col(phi2, ["phi"], ["norm"]) or find_col(phi2, "phi_norm")

test_cols = [
    ("NSS", col_nss),
    ("MaxChi2", col_maxchi2),
    ("PHI_perm", col_phiperm),
    ("PHI_norm", col_phinorm),
]

missing = [name for name, c in test_cols if c is None]
if missing:
    raise KeyError(
        f"Could not find columns for: {missing}\n"
        f"Available columns: {list(phi2.columns)}\n\n"
        "Tip: If your headers differ (e.g., 'MaxChi²_p' or 'PHI_perm_pvalue'), rename them or "
        "edit the find_col() logic."
    )

# -----------------------------
# Convert p-value columns to numeric (FIX for your error)
# -----------------------------
for _, pcol in test_cols:
    phi2[pcol] = to_numeric_clean(phi2[pcol])

# -----------------------------
# Build categorical matrices: 0=no data, 1=absent, 2=present
# -----------------------------
sig_thr = 0.05
mats = {}
for test_name, pcol in test_cols:
    piv = phi2.pivot_table(index="gene", columns="species_key", values=pcol, aggfunc="min")
    piv = piv.reindex(index=gene_order, columns=species_order)

    cat = pd.DataFrame(0, index=piv.index, columns=piv.columns, dtype=int)
    cat[(piv >= sig_thr) & (~piv.isna())] = 1
    cat[(piv <  sig_thr) & (~piv.isna())] = 2

    mats[test_name] = cat

# -----------------------------
# Plot: 1 page, 4 panels side-by-side
# -----------------------------
# Colors: grey/orange/green
cmap = ListedColormap(["#d9d9d9", "#f4a742", "#2ca25f"])
norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

fig, axes = plt.subplots(
    1, 4,
    figsize=(22, 10),
    gridspec_kw={"wspace": 0.12}
)

for ax, (test_name, _) in zip(axes, test_cols):
    cat = mats[test_name]
    ax.imshow(cat.values, aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")

    # x labels (species)
    ax.set_xticks(range(len(cat.columns)))
    ax.set_xticklabels([species_labels.get(s, s) for s in cat.columns],
                       rotation=45, ha="right", fontsize=10)

    # y labels (genes) only on first panel
    ax.set_yticks(range(len(cat.index)))
    if ax is axes[0]:
        ax.set_yticklabels(cat.index, fontsize=10)
        ax.set_ylabel("Gene", fontsize=12)
    else:
        ax.set_yticklabels([])

    ax.set_title(test_name, fontsize=13)

    # black grid between cells
    ax.set_xticks(np.arange(-0.5, cat.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, cat.shape[0], 1), minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=0.4)
    ax.tick_params(which="minor", bottom=False, left=False)

# shared x label
#fig.supxlabel("Species", fontsize=12, y=0.04)

# Legend (black border)
legend_handles = [
    Patch(facecolor="#2ca25f", edgecolor="black", label="Present (p < 0.05)"),
    Patch(facecolor="#f4a742", edgecolor="black", label="Absent (p ≥ 0.05)"),
    Patch(facecolor="#d9d9d9", edgecolor="black", label="No data"),
]
fig.legend(
    handles=legend_handles,
    loc="upper center",
    ncol=3,
    frameon=True,
    bbox_to_anchor=(0.5, 0.98),
    fontsize=11
)

fig.suptitle("Non-parametric recombination tests across Paramecium mitochondrial genes", fontsize=14, y=0.995)

# Save
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("Wrote:", out_png)
print("Wrote:", out_pdf)
print("Detected columns:",
      {name: col for name, col in test_cols})
