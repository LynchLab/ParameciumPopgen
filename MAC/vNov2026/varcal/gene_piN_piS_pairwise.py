### calculate genes' piN and piS as mean pairwise difference divided by no. of sites
## correct for breadth of coverage

import os, sys, pandas, multiprocessing

indr = sys.argv[1] # input directory of genes' pairwise comparison tables
fcns = sys.argv[2] # input site-count table
fbcv = sys.argv[3] # input coverage table
fout = sys.argv[4] # output table
ncpu = int(sys.argv[5])

# load input files of genes' pairwise differences among samples
fls = os.listdir(indr)
nf = len(fls)
genes = [ i.replace(".csv","") for i in fls]

# get average breadth of coverage for genes across samples
dbcv = pandas.read_csv( fbcv, sep='\t',index_col=0).mean(axis=1)
dcns = pandas.read_csv( fcns, index_col=0)
# only keep genes with coverage data
dcns = dcns.loc[ dbcv.index, :]
genes = [ i for i in genes if i in dbcv.index]
# add coverage info
dcns['cov'] = dbcv
# correct number of sites for coverage
dcns['syn_cor'] = dcns['syn'] * dcns['cov']
dcns['nsyn_cor'] = dcns['nsyn'] * dcns['cov']

def get_pNps(gene):
    fin = os.path.join( indr, gene+'.csv')
    df = pandas.read_csv(fin)
    # average pair-wise difference
    Pn = df['Pn'].mean()
    Ps = df['Ps'].mean()
    # normalize by no. of sites
    nN = dcns.loc[gene,'nsyn_cor']
    nS = dcns.loc[gene,'syn_cor']
    pN = Pn/nN
    pS = Ps/nS
    return [ gene, Pn, Ps, nN, nS, pN, pS, pN/pS]

with multiprocessing.Pool(ncpu) as mp_pool:
    out = mp_pool.map( get_pNps, genes)

out = pandas.DataFrame(out)
out.columns = ['gene', 'Dn', 'Ds', 'nonsyn_sites', 'syn_sites', 'piN', 'piS', 'pN/pS']
out.to_csv( fout, index=False)
