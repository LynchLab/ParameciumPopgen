### Make tables of number of syn and non-syn differences between pairs of samples for genes
import os, sys, pandas, multiprocessing

fin = sys.argv[1] # input snps_aa_gts table
dout = sys.argv[2] # output dir of distance matrices for genes
ncpu = int(sys.argv[3]) # number of parallel processes

# create output dir if missing
not os.path.exists(dout) and os.mkdir(dout)

# load input table
df = pandas.read_csv( fin, sep='\t')
# split data by feature ids (genes)
gene_grps = df.groupby('feature_id')

### Functions ############
## get mismatch count (delta) between two diploid genotypes at a site
def delta_site_gt(gt1,gt2):
    # separate and sort characters excluding uncalled genotypes
    chars1 = sorted( [i for i in gt1.split('/') if i != "."])
    chars2 = sorted( [i for i in gt2.split('/') if i != "."])
    # if any of the two genotypes has less than 2 called alleles, then delta cannot be determined
    if len(chars1) < 2 or len(chars2) < 2:
        return
    
    # delta is zero if character sorted genotypes are the same
    if chars1 == chars2:
        delta = 0
    # otherwise, record the boolean vector of base-equality for comparisons between the two genotypes A & B
    # as [ A1B1, A1B2, A2B1, A2B2]
    else:
        cfg = [ int( chars1[i] == chars2[j]) for i in range(2) for j in range(2)]
        sum_cfg = sum(cfg)
        # if no two bases are equal between the genotypes, then genotypes are completely different
        if sum_cfg == 0:
            delta = 1
        # otherwise, they differ on 1 of the 2 chromosomes
        else:
            delta = 0.5
    return delta

## generate pairwise distance matrix for the gene and save output
def write_gene_pdmat(gene):

    fout = os.path.join( dout, gene+".csv")
    if os.path.exists(fout):
        print("output exists! skipping %s\n"%gene)
        return
    
    # extract gene data and remove unnecessary info
    gdf = gene_grps.get_group(gene).iloc[:,4:]
    nsites,nsamples = gdf.shape
    npairs = int(nsamples * (nsamples-1)/2)
    samples = gdf.columns
    # separate SNPs and AAs genotypes
    gts = gdf.apply( lambda x: x.apply( lambda y: y.split('|')))
    snps = gts.apply( lambda X: X.apply( lambda y: y[0]))
    aas = gts.apply( lambda X: X.apply( lambda y: y[1]))
    # initialize empty output dataframe
    outdf = pandas.DataFrame( index = range(npairs), columns = ["sample1","sample2","N","Pn","Ps"])
    # do pairwise comparisons among samples
    p = 0
    for i in range(nsamples-1):
        for j in range(i+1,nsamples):
            # compute Delta over the gene length for SNPs and AAs separately
            Dsnp = snps.iloc[:,[i,j]].apply( lambda x: delta_site_gt(x.iat[0],x.iat[1]),axis=1).sum()
            Daa = aas.iloc[:,[i,j]].apply( lambda x: delta_site_gt(x.iat[0],x.iat[1]),axis=1).sum()
            outdf.loc[p,:] = [samples[i],samples[j],Dsnp,Daa,Dsnp-Daa]
            p = p + 1

    outdf.to_csv( fout, index=False)
####################################################

# list of all genes with mutations
genes = list(gene_grps.groups.keys())

# parallel run over genes
with multiprocessing.Pool(ncpu) as mp_pool:
    mp_pool.map( write_gene_pdmat, genes)
