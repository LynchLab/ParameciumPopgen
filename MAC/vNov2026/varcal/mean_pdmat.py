## take mean over pairwise distance matrices

import os, sys, pandas, re, multiprocessing, numpy

indr = sys.argv[1] # input dir of distance matrices
fout = sys.argv[2] # output table of mean pairwise distances
ncpu = int(sys.argv[3]) # number of CPUs

# list of input files
fls = os.listdir(indr)
nf = len(fls)

# function to load input matrices
def load_input(fname,dr=indr):
    df = pandas.read_csv( os.path.join( indr, fname), index_col=0)
    return df

# the function above is run parallely 
with multiprocessing.Pool(ncpu) as mp_pool:
	dfs = mp_pool.map( load_input, fls)

# pandas groupby is used to take mean over common pairs across files
out = pandas.concat(dfs).groupby(level=0).mean()
out = out.loc[ out.columns, :]
out.to_csv( fout)
