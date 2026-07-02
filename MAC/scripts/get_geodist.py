import os, sys, pandas
from geopy.distance import geodesic

fin = sys.argv[1]
dout = sys.argv[2]

df = pandas.read_csv(fin)
# remove samples without GPS coordinates
df = df[ ~df.latitude.isna()]
# group data by species
spdfs = df.groupby(df.species)

def afunc(sdf):
	fout = os.path.join( dout, sdf.name+'.csv')
	nsam = sdf.shape[0]
	out = pandas.DataFrame( [ [ geodesic((sdf.iat[i,2],sdf.iat[i,3]),(sdf.iat[j,2],sdf.iat[j,3])).km for j in range(nsam)] for i in range(nsam)])
	out.index = sdf['sample']
	out.columns = sdf['sample']
	out.to_csv( fout)

spdfs.apply( afunc, include_groups=False)
