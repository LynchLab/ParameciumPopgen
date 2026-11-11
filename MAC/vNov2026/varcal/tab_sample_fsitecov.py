### generate a table of fraction of sites covered for each sample
import sys, os, pandas

indr = sys.argv[1] # dir of input files
fout = sys.argv[2] # output table 

samples = os.listdir(indr)
N = len(samples)

bcovs = pandas.concat( (pandas.read_csv( os.path.join( indr, i), 
	sep='\t')[['id','bcov']].groupby(['id']).mean().round(3)
for i in samples), axis=1)

# change column names
bcovs.columns = samples

# replace NaN with Zero
bcovs.fillna(0, inplace=True)

bcovs.to_csv( fout, sep='\t')

