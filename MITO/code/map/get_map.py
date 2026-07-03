import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# === Load your data ===
df = pd.read_csv("315_species_genetic_geographic_mito/315_geographic_cordinate_with_pop.csv")
df.columns = df.columns.str.replace("'", "").str.strip()
df['species'] = df['species'].str.replace("'", "").str.strip()
df['latitude'] = df['latitude'].astype(float)
df['longitude'] = df['longitude'].astype(float)

# === Mapping abbreviated species to full names ===
species_fullname_map = {
    'bursaria': 'P. bursaria',
    'multi': 'P. multimicronucleatum',
    'fokini': 'P. fokini',
    'caudatum': 'P. caudatum',
    'sexaurelia': 'P. sexaurelia',
    'tetraurelia': 'P. tetraurelia',
    'tredecaurelia': 'P. tredecaurelia',
    'quadecaurelia': 'P. quadecaurelia',
    'Psch': 'P. schewiakoffi',
    'quindecaurelia': 'P. quindecaurelia',
    'biaurelia': 'P. biaurelia',
    'triaurelia': 'P. triaurelia',
    'primaurelia': 'P. primaurelia',
    'tribursaria': 'P.tribursaria'
}

# === Final custom color map ===
color_map = {
    'P.tribursaria':'#88c96b',
    'P. bursaria': '#cc0000',
    'P. multimicronucleatum': '#9874f4',
    'P. fokini': '#ffe599',
    'P. caudatum': '#e69138',
    'P. sexaurelia': '#559738',
    'P. tetraurelia': '#0b5394',
    'P. tredecaurelia': '#f1c232',
    'P. quadecaurelia': '#ccd9c7',
    'P. schewiakoffi': '#aaaaaa',
    'P. quindecaurelia': '#d5ed8a',
    'P. biaurelia': '#ee79ac',
    'P. triaurelia': '#d5a6bd',
    'P. primaurelia': '#6fa8dc'

}

# === Custom legend order ===
legend_order = list(color_map.keys())

# === Map full names ===
df['full_species'] = df['species'].map(species_fullname_map).fillna(df['species'])

# === Create GeoDataFrame ===
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df['longitude'], df['latitude']),
    crs="EPSG:4326"
)

# === Add jitter ===
counts = gdf.groupby(['latitude', 'longitude']).size().reset_index(name='count')
gdf = gdf.merge(counts, on=['latitude', 'longitude'])
jitter_strength = 5
gdf['longitude_jitter'] = gdf.apply(
    lambda row: row['longitude'] + jitter_strength * (np.random.rand() - 0.5) if row['count'] > 1 else row['longitude'],
    axis=1
)
gdf['latitude_jitter'] = gdf.apply(
    lambda row: row['latitude'] + jitter_strength * (np.random.rand() - 0.5) if row['count'] > 1 else row['latitude'],
    axis=1
)

# === Load world map ===
world = gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip")

# === Plot ===
#fig, ax = plt.subplots(figsize=(16, 10))
#world.plot(ax=ax, color='lightgrey', edgecolor='white')

# === Plot setup with reduced size ===
fig, ax = plt.subplots(figsize=(12, 12))
world.plot(ax=ax, color='lightgrey', edgecolor='white')




# === Plot each species ===
for species in legend_order:
    subset = gdf[gdf['full_species'] == species]
    if not subset.empty:
        ax.scatter(
            subset['longitude_jitter'],
            subset['latitude_jitter'],
            color=color_map.get(species, '#999999'),
            label=species,
            edgecolor='black',
            s=40,
            alpha=0.9
        )

# === Ordered legend with italics ===
handles = [
    plt.Line2D([0], [0], marker='o', color='w', label=rf"$\it{{{sp}}}$",
               markerfacecolor=color_map[sp], markeredgecolor='black', markersize=10)
    for sp in legend_order if sp in gdf['full_species'].values
]
#ax.legend(handles=handles, title='Species', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title_fontsize=14)
# === Legend with adjusted anchor ===
ax.legend(
    handles=handles,
    title='Species',
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),  # Slightly closer to the plot
    fontsize=12,
    title_fontsize=14
)

# === Final touches ===
# ax.set_title("Species Sample Locations on World Map", fontsize=18)
# ax.set_xlabel("Longitude")
# ax.set_ylabel("Latitude")
plt.grid(True)

# Save
plt.tight_layout()
plt.savefig("species_map_final_colored_italic.pdf", format='pdf', bbox_inches='tight')
plt.savefig("species_map_final_colored_italic.png", dpi=300, bbox_inches='tight')  # <-- PNG export
plt.show()



