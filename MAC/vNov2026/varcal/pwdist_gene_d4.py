### Get distance matrices for different genes from SNPs_AA_GTs table generated from VCFs
### calculated from sites of a given degeneracy

import os, sys, pandas, multiprocessing

fin = sys.argv[1] # input snps_aa_gts table
dout = sys.argv[2] # output dir of distance matrices for genes
fgl = sys.argv[3] # genes' breadth of coverage across samples
fgp = sys.argv[4] # snp gene position info table
ddeg = sys.argv[5] # input dir of genes' degeneracy sequences
deg = int(sys.argv[6]) # degeneracy
ncpu = int(sys.argv[7]) # number of parallel processes

# create output dir if missing
not os.path.exists(dout) and os.mkdir(dout)

# load input table
df = pandas.read_csv(fin,sep='\t')
# genes' breadth of coverage data
dfgl = pandas.read_csv( fgl, sep='\t', index_col=0)
# gene position data
dfgp = pandas.read_csv( fgp, sep='\t')

# list genes' degeneracy files
dgls = os.listdir(ddeg)
dgls_genes = [ i.replace('.txt','') for i in dgls]

# no. of genes with deg
ngenes_deg = len(dgls)

# no. of samples
nsamples = dfgl.shape[1]

# only keep sites within genes in the nsites table
df = df[ df.feature_id.isin( dfgl.index)]
# only keep sites within genes with site-deg
df = df[ df.feature_id.isin( dgls_genes)]

# set chrom, pos, gene as index on gene-pos table
dfgp.set_index(['CHROM','POS','Feature_ID'], inplace=True)

# extract gene-pos info for sites in the AA GT table
def get_npos(X):
	return dfgp.loc[ (X.chrom,X.pos,X.feature_id), 'npos'].iloc[0]
# add a column to genotype table for mutated position wrt the gene
df['npos'] = df.apply( get_npos, axis=1)

# split genotype data by feature ids (genes)
gene_grps = df.groupby('feature_id')

### Functions ############
## get mismatch count (delta) between two diploid genotypes at a site
def delta_site_gt(gt1,gt2):
    # separate and sort characters excluding uncalled genotypes
    chars1 = sorted( [i for i in gt1.split('/') if i != "."])
    chars2 = sorted( [i for i in gt2.split('/') if i != "."])
    if len(chars1) < 2 or len(chars2) < 2:
        return
    
    # delta is zero if genotypes are the same after sorting characters
    if chars1 == chars2:
        delta = 0
    else:
        cfg = [ int( chars1[i] == chars2[j]) for i in range(2) for j in range(2)]
        sum_cfg = sum(cfg)
        if sum_cfg == 0:
            delta = 1
        else:
            delta = 0.5
    return delta

## generate pairwise distance matrix for the gene and save output
def write_gene_pdmat(gene):
    fout = os.path.join( dout, gene+".csv")
    
    if os.path.exists(fout):
        print("output exists! skipping %s\n"%gene)
        return
   
    # site-deg sequence
    fdeg = os.path.join( ddeg, gene+'.txt')
    with open(fdeg) as flob:
        degs = { (x+1):int(i) for x,i in enumerate(flob.read().rstrip('\n'))}
    
    # extract genotype data for the gene
    gdf = gene_grps.get_group(gene).copy()
    # add a column for site-degeneract to the gene-subset of genotype table
    gdf['deg'] = gdf.npos.apply( lambda x: degs[x])
    # no. of sites of the given deg in the reference gene
    nsites = sum( v==deg for v in degs.values())
    # skip the gene, if the reference has no sites of the given deg
    if nsites == 0:
        return
    
    # subset genotype data to the sites of the given deg
    gdf = gdf[ gdf['deg'].eq(deg)]
    # extract SNPs genotypes
    sdf = gdf.iloc[:,4:-2]
    sdf = sdf.apply( lambda x: x.apply( lambda y: y.split('|')[0]))
    samples = sdf.columns
    nsamples = len(samples)

    # initialize empty output dataframe
    outdf = pandas.DataFrame( index = samples, columns = samples)
    
    # do pairwise comparisons among samples  
    # if no site of the given degeneracy is present in the gene-subset of the genotype (mutated) table, then distance among all samples for that gene is Zero, 
    # except when there are no informative sites ( due to low coverage )
    if gdf.shape[0] == 0:
        for i in range(nsamples-1):
            for j in range(i+1,nsamples):
                # get number of info sites for this comparison
                ninfo = dfgl.loc[ gene, [ samples[i], samples[j]]].min()*nsites
                # skip the comparisons if no info sites
                if ninfo == 0: continue
                # otherwise the distance is zero, as no mutated sites seen
                outdf.iat[i,j] = 0
    else:
        for i in range(nsamples-1):
            for j in range(i+1,nsamples):
                # get number of info sites for this comparison
                ninfo = dfgl.loc[ gene, [ samples[i], samples[j]]].min()*nsites
                # skip the comparisons if no info sites
                if ninfo == 0: continue 
                # compute Delta over the gene length
                ndif = sdf.iloc[ :, [i,j]].apply( lambda x: delta_site_gt( x.iat[0], x.iat[1]), axis=1).sum()
                # divide total number of differences by the number of informative sites to get the distance
                outdf.iat[i,j] = ndif/ninfo

    outdf.to_csv(fout)
####################################################


# list of all genes with gene length info
genes = list( gene_grps.groups.keys())

# parallel run over genes
with multiprocessing.Pool(ncpu) as mp_pool:
    mp_pool.map( write_gene_pdmat, genes)

