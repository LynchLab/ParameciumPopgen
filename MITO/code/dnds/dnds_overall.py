import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# =========================
# Load the data
# =========================
file_path = "M8_beta_p_q_omega.csv"
df = pd.read_csv(file_path)

# =========================
# Define functional group mapping
# =========================
functional_group_map = {
    # Electron Transport genes
    "atp9": "Electron Transport",
    "cob": "Electron Transport",
    "cox1": "Electron Transport",
    "cox2": "Electron Transport",
    "nadh1-a": "Electron Transport",
    "nadh1-b": "Electron Transport",
    "nadh10": "Electron Transport",
    "nadh2": "Electron Transport",
    "nadh3": "Electron Transport",
    "nadh4": "Electron Transport",
    "nadh4L": "Electron Transport",
    "nadh5": "Electron Transport",
    "nadh6": "Electron Transport",
    "nadh7": "Electron Transport",
    "nadh9": "Electron Transport",
    "yejR": "Electron Transport",

    # Ribosomal Proteins
    "rpL14": "Ribosomal Proteins",
    "rpL16": "Ribosomal Proteins",
    "rpL2":  "Ribosomal Proteins",
    "rpL32": "Ribosomal Proteins",
    "rpL6":  "Ribosomal Proteins",
    "rpS12": "Ribosomal Proteins",
    "rpS13": "Ribosomal Proteins",
    "rpS14": "Ribosomal Proteins",
    "rpS19": "Ribosomal Proteins",
    "rpS3":  "Ribosomal Proteins",

    # Ciliate-Specific Genes
    "ymf56": "Ciliate-Specific Genes",
    "ymf57": "Ciliate-Specific Genes",
    "ymf59": "Ciliate-Specific Genes",
    "ymf61": "Ciliate-Specific Genes",
    "ymf63": "Ciliate-Specific Genes",
    "ymf64": "Ciliate-Specific Genes",
    "ymf65": "Ciliate-Specific Genes",
    "ymf66": "Ciliate-Specific Genes",
    "ymf67": "Ciliate-Specific Genes",
    "ymf68": "Ciliate-Specific Genes",
    "ymf76": "Ciliate-Specific Genes",
    "ymf78": "Ciliate-Specific Genes",
    "ymf79": "Ciliate-Specific Genes",
    "ymf80": "Ciliate-Specific Genes",
    "ymf83": "Ciliate-Specific Genes",
    "ymf84": "Ciliate-Specific Genes"
}

# =========================
# Add functional group column
# =========================
df['functional_group'] = df['gene'].map(functional_group_map)

# Filter out rows with no mapping
df_filtered = df.dropna(subset=['functional_group'])

# =========================
# Define fixed functional group order
# =========================
group_order = ["Electron Transport", "Ribosomal Proteins", "Ciliate-Specific Genes"]
df_filtered['functional_group'] = pd.Categorical(df_filtered['functional_group'], categories=group_order, ordered=True)

# =========================
# Sort by functional group order and omega descending
# =========================
df_sorted = df_filtered.sort_values(by=['functional_group', 'omega'], ascending=[True, False])

# =========================
# Define custom colors
# =========================
palette_colors = {
    "Electron Transport": "#ff7f0e",     # Blue
    "Ribosomal Proteins": "#2ca02c",     # Orange
    "Ciliate-Specific Genes": "#1f77b4" # Green
}

# Assign colors
bar_colors = [palette_colors[group] for group in df_sorted['functional_group']]

# =========================
# Create the bar plot
# =========================
plt.figure(figsize=(14, 7))
plt.bar(df_sorted['gene'], df_sorted['omega'], color=bar_colors)

# Labels and title
plt.xlabel('Gene', fontsize=24)
plt.ylabel('dN/dS (Omega)', fontsize=24)
plt.title('dN/dS Values Across Genes by Functional Group', fontsize=30)

# Rotate x-axis labels for readability
plt.xticks(rotation=45, ha='right',fontsize = 16)
plt.yticks(fontsize=18)
# Add gridlines
plt.grid(axis='y', linestyle='--', alpha=0.7)

handles = [
    Patch(facecolor=palette_colors[g], edgecolor="black", label=g)
    for g in group_order
]

legend = plt.legend(
    handles=handles,
    title="Functional Group",
    fontsize=16,          # ↓ smaller legend text
    title_fontsize=18,    # slightly larger title
    loc="upper right",
    frameon=True,
    borderpad=0.5,
    labelspacing=0.4,
    handletextpad=0.5
)

legend.get_frame().set_edgecolor("black")
legend.get_frame().set_linewidth(1.2)
legend.get_frame().set_facecolor("white")
# Adjust layout
plt.tight_layout()

# Save the plot
plt.savefig("dnds_barplot_grouped_fixed_order.png", format="png")

# Show the plot
plt.show()

