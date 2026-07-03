import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg

# =========================
# INPUT FILES
# =========================
csv_file = "315_species_genetic_geographic_mito/315_geographic_cordinate_with_pop.csv"
tree_image_file = "tree_final.png"   # use a tightly cropped tree image
world_map_file = "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"

# =========================
# LOAD DATA
# =========================
df = pd.read_csv(csv_file)
df.columns = df.columns.str.replace("'", "", regex=False).str.strip()
df["species"] = df["species"].astype(str).str.replace("'", "", regex=False).str.strip()
df["latitude"] = df["latitude"].astype(float)
df["longitude"] = df["longitude"].astype(float)

# =========================
# MAP ABBREVIATED TO FULL SPECIES NAMES
# =========================
species_fullname_map = {
    "bursaria": "P. bursaria",
    "multi": "P. multimicronucleatum",
    "fokini": "P. fokini",
    "caudatum": "P. caudatum",
    "sexaurelia": "P. sexaurelia",
    "tetraurelia": "P. tetraurelia",
    "tredecaurelia": "P. tredecaurelia",
    "quadecaurelia": "P. quadecaurelia",
    "Psch": "P. schewiakoffi",
    "quindecaurelia": "P. quindecaurelia",
    "biaurelia": "P. biaurelia",
    "triaurelia": "P. triaurelia",
    "primaurelia": "P. primaurelia",
    "tribursaria": "P. tribursaria"
}

df["full_species"] = df["species"].map(species_fullname_map).fillna(df["species"])

# =========================
# FINAL COLOR MAP
# =========================
color_map = {
    "P. tribursaria": "#88c96b",
    "P. bursaria": "#cc0000",
    "P. multimicronucleatum": "#9874f4",
    "P. fokini": "#ffe599",
    "P. caudatum": "#e69138",
    "P. sexaurelia": "#559738",
    "P. tetraurelia": "#0b5394",
    "P. tredecaurelia": "#f1c232",
    "P. quadecaurelia": "#ccd9c7",
    "P. schewiakoffi": "#4c7f88",
    "P. quindecaurelia": "#d5ed8a",
    "P. biaurelia": "#ee79ac",
    "P. triaurelia": "#d5a6bd",
    "P. primaurelia": "#6fa8dc"
}

# =========================
# CREATE GEODATAFRAME
# =========================
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df["longitude"], df["latitude"]),
    crs="EPSG:4326"
)

# =========================
# ADD JITTER FOR OVERLAPPING POINTS
# =========================
np.random.seed(42)

counts = gdf.groupby(["latitude", "longitude"]).size().reset_index(name="count")
gdf = gdf.merge(counts, on=["latitude", "longitude"], how="left")

jitter_strength = 3.0

gdf["longitude_jitter"] = gdf.apply(
    lambda row: row["longitude"] + jitter_strength * (np.random.rand() - 0.5)
    if row["count"] > 1 else row["longitude"],
    axis=1
)

gdf["latitude_jitter"] = gdf.apply(
    lambda row: row["latitude"] + jitter_strength * (np.random.rand() - 0.5)
    if row["count"] > 1 else row["latitude"],
    axis=1
)

# =========================
# LOAD WORLD MAP
# =========================
world = gpd.read_file(world_map_file)

# =========================
# FIGURE SETUP
# Slightly taller overall
# =========================
fig, axes = plt.subplots(
    1, 2,
    figsize=(14, 9),
    gridspec_kw={"width_ratios": [1.85, 1.15]}
)

# =========================
# PANEL A: MAP
# =========================
ax_map = axes[0]
ax_map.set_box_aspect(0.48)   # slightly taller map panel

world.plot(ax=ax_map, color="#cfcfcf", edgecolor="white", linewidth=0.4)

for species, color in color_map.items():
    subset = gdf[gdf["full_species"] == species]
    if not subset.empty:
        ax_map.scatter(
            subset["longitude_jitter"],
            subset["latitude_jitter"],
            color=color,
            edgecolor="black",
            linewidth=0.7,
            s=65,
            alpha=0.95,
            zorder=3
        )

ax_map.set_xlim(-180, 180)
ax_map.set_ylim(-60, 85)

ax_map.set_xlabel("Longitude", fontsize=17)
ax_map.set_ylabel("Latitude", fontsize=17)
ax_map.tick_params(axis="both", labelsize=13, width=1.2, length=6)

for spine in ax_map.spines.values():
    spine.set_linewidth(1.2)

ax_map.text(
    0.01, 1.12, "A",
    transform=ax_map.transAxes,
    fontsize=20,
    fontweight="bold",
    va="top",
    ha="left"
)

# =========================
# PANEL B: TREE
# =========================
ax_tree = axes[1]

tree_img = mpimg.imread(tree_image_file)
ax_tree.imshow(tree_img)
ax_tree.axis("off")

ax_tree.text(
    0.01, 1.1, "B",
    transform=ax_tree.transAxes,
    fontsize=20,
    fontweight="bold",
    va="top",
    ha="left"
)

# =========================
# LAYOUT + SAVE
# =========================
plt.subplots_adjust(left=0.06, right=0.98, top=0.96, bottom=0.14, wspace=0.06)

plt.savefig("Figure1_map_tree_clean.pdf", format="pdf", bbox_inches="tight")
plt.savefig("Figure1_map_tree_clean.png", dpi=600, bbox_inches="tight")
plt.show()