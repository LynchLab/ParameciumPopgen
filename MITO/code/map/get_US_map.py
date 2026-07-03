# Full script: US map (CONUS only), state abbrevs, jittered points,
# compact legend in lower-left

import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# 1) Load & clean your data
# -----------------------------
df = pd.read_csv("315_species_genetic_geographic_mito/USA_samples.csv")
df.columns = df.columns.str.replace("'", "").str.strip()
df['species']   = df['species'].str.replace("'", "").str.strip()
df['latitude']  = df['latitude'].astype(float)
df['longitude'] = df['longitude'].astype(float)

# -----------------------------
# 2) Species mapping & colors
# -----------------------------
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
legend_order = list(color_map.keys())

df['full_species'] = df['species'].map(species_fullname_map).fillna(df['species'])

# -----------------------------
# 3) GeoDataFrame + jitter
# -----------------------------
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df['longitude'], df['latitude']),
    crs="EPSG:4326"
)

# jitter only stacked coords (deterministic)
counts = gdf.groupby(['latitude', 'longitude']).size().reset_index(name='count')
gdf = gdf.merge(counts, on=['latitude', 'longitude'])
rng = np.random.default_rng(42)
jitter_strength = 0.35  # degrees (~15–20 km)
rand_lon = jitter_strength * (rng.random(len(gdf)) - 0.5)
rand_lat = jitter_strength * (rng.random(len(gdf)) - 0.5)
gdf['longitude_jitter'] = np.where(gdf['count'] > 1, gdf['longitude'] + rand_lon, gdf['longitude'])
gdf['latitude_jitter']  = np.where(gdf['count'] > 1, gdf['latitude']  + rand_lat, gdf['latitude'])

# -----------------------------
# 4) USA states, keep only CONUS
# -----------------------------
states = gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_1_states_provinces.zip")
usa_states = states[states['admin'] == 'United States of America'].copy()

exclude_names = {
    'Alaska', 'Hawaii', 'Puerto Rico', 'Guam',
    'United States Virgin Islands', 'American Samoa',
    'Commonwealth of the Northern Mariana Islands'
}
conus = usa_states[~usa_states['name'].isin(exclude_names)].copy()

# state abbreviations
state_abbrev = {
    "New York": "NY", "California": "CA", "Texas": "TX", "Florida": "FL",
    "Arizona": "AZ", "Massachusetts": "MA", "Illinois": "IL", "Ohio": "OH",
    "North Carolina": "NC", "South Carolina": "SC", "Pennsylvania": "PA",
    "New Jersey": "NJ", "Connecticut": "CT", "Maine": "ME", "New Hampshire": "NH",
    "Oregon": "OR", "Washington": "WA", "Colorado": "CO", "Kansas": "KS",
    "Missouri": "MO", "Louisiana": "LA", "Georgia": "GA", "Nevada": "NV",
    "Michigan": "MI", "Indiana": "IN", "Iowa": "IA", "Wisconsin": "WI",
    "Minnesota": "MN", "Nebraska": "NE", "South Dakota": "SD", "North Dakota": "ND",
    "Utah": "UT", "Wyoming": "WY", "Idaho": "ID", "Montana": "MT",
    "Oklahoma": "OK", "Arkansas": "AR", "Tennessee": "TN", "Kentucky": "KY",
    "Virginia": "VA", "West Virginia": "WV", "Maryland": "MD", "Delaware": "DE",
    "Rhode Island": "RI", "Vermont": "VT", "Alabama": "AL", "Mississippi": "MS",
    "District of Columbia": "DC"
}
conus['short_name'] = conus['name'].map(state_abbrev).fillna(conus['name'])
conus['label_pt']   = conus.representative_point()

# -----------------------------
# 5) Plot
# -----------------------------
fig, ax = plt.subplots(figsize=(16, 10), constrained_layout=True)

# Draw states
conus.plot(ax=ax, color='lightgrey', edgecolor='black', linewidth=0.8, zorder=1)

# State abbrev labels
for _, row in conus.iterrows():
    ax.text(row['label_pt'].x, row['label_pt'].y, row['short_name'],
            fontsize=9, ha='center', va='center',
            color='midnightblue', alpha=0.7, zorder=2)

# Sample points by species
for species in legend_order:
    subset = gdf[gdf['full_species'] == species]
    if len(subset):
        ax.scatter(
            subset['longitude_jitter'], subset['latitude_jitter'],
            color=color_map.get(species, '#999999'),
            edgecolor='black', s=50, alpha=0.9, zorder=3, label=species
        )

# -----------------------------
# Compact Legend at Lower Left
# -----------------------------
handles = [
    plt.Line2D([0], [0], marker='o', color='w',
               label=rf"$\it{{{sp}}}$",
               markerfacecolor=color_map[sp],
               markeredgecolor='black',
               markersize=7)  # smaller markers
    for sp in legend_order if (gdf['full_species'] == sp).any()
]
ax.legend(
    handles=handles,
    title='Species',
    loc='lower right',          # Lower left inside
    frameon=True,
    framealpha=0.9,
    fontsize=11,                # smaller text
    title_fontsize=13,
    borderpad=0.4,             # reduce border padding
    labelspacing=0.4,          # reduce space between labels
    handletextpad=0.5,          # reduce space between marker and text
    borderaxespad=0.4           # space between legend and axes edge
)

# -----------------------------
# Map extent for contiguous U.S.
# -----------------------------
ax.set_xlim(-125, -66)
ax.set_ylim(24, 50)
ax.set_aspect('equal', adjustable='box')

# Labels and grid
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.grid(True, alpha=0.4)

# Save & show
plt.savefig("species_map_CONUS_legend_compact.pdf", format='pdf', bbox_inches='tight')
plt.show()
