import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np

# === Load your data ===
df = pd.read_csv("data/table_S1.csv")
df['latitude'] = df['latitude'].astype(float)
df['longitude'] = df['longitude'].astype(float)
df = df[ ~df.isna().any(axis=1)]
df = df[['species','population','latitude','longitude']].drop_duplicates()
df['species'] = df.species.str.replace("_clade.+","",regex=True)
df['species'] = 'P. ' + df.species

# load color map
clrs = pd.read_csv("species_colors_pts.csv")
clrs = clrs[~clrs.species.str.contains("cladeB")]
clrs['species'] = clrs['species'].str.replace("_cladeA", "") 
clrs = clrs.iloc[ [9,6,8,7,5,4,0,10,3,1,2], :]
clrs['species'] = 'P. ' + clrs.species
# point types
clrs['ptype'] = 'o'
clrs.loc[clrs['pch'].eq(2),'ptype'] = '^'
clrs.set_index('species',inplace=True)


# === Custom legend order ===
legend_order = list(clrs.index)
legend_order = [ i.replace("_", " ") for i in legend_order]


# === Create GeoDataFrame ===
gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df['longitude'], df['latitude']),
    crs="EPSG:4326"
)


# === Load world map ===
world = gpd.read_file("https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip")

# drop "Antarctica" and "Seven seas" from the dataframe 
drop_idxs = world["CONTINENT"].isin([
    "Antarctica",
    "Seven seas (open ocean)"
])
world = world.drop(world[drop_idxs].index)

# === Plot setup ===
fig, ax = plt.subplots(figsize=(10, 20))
world.plot(ax=ax, color='lightgrey', edgecolor='white')

# turn off axis ticks
ax.set_xticks([])
ax.set_yticks([])


# === Plot each species ===
for species in legend_order:
    subset = gdf[ gdf.species == species]
    ax.scatter(
        subset['longitude'],
        subset['latitude'],
        color=clrs.loc[ species, 'hex'],
        marker=clrs.loc[ species, 'ptype'],
        label=species,
        edgecolor='black',
        linewidth=0.1,
        s=40,
        alpha=0.9
    ) 


# === Ordered legend with italics ===
handles = [
    plt.Line2D([0], [0], marker=clrs.loc[sp,'ptype'], color='w', label=rf"$\it{{{sp}}}$",
               markerfacecolor=clrs.loc[sp,'hex'], 
               markersize=10) for sp in legend_order]
ax.legend(handles=handles, title='', loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12, title_fontsize=14)
# === Legend with adjusted anchor ===
ax.legend(
     handles=handles,
     title='',
     loc='center left',
     bbox_to_anchor=(1.02, 0.5),  # Slightly closer to the plot
     fontsize=12,
     title_fontsize=14
)


plt.grid(True)

# Save
plt.tight_layout()
plt.savefig("images/figure_1A.png", dpi=300, bbox_inches='tight')
