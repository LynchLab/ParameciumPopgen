### Create a table of variants in protein-coding genes with corresponding amino acid changes across samples

## require
import sys, pandas

fvcf = sys.argv[1] # input VCF file
fann = sys.argv[2] # TSV file of syn/nonsyn annotations
fout = sys.argv[3] # output table of amino-acid genotypes

## extract info from VCF
vcf = pandas.read_csv( fvcf, sep='\t', comment='#', header=None)

# no. of columns
nvar,ncol = vcf.shape

# find the index of GT:PL:DP column
ix = [ x for x,i in enumerate(vcf.iloc[0].str.contains('GT:PL:DP').values) if i is True][0]

# extract GT bit from sample columns
gts = vcf.iloc[:,ix+1:].apply( lambda X: X.apply( lambda y: y.split(':')[0]))

# number of samples
nsamples = gts.shape[1]

# extract columns : chrom, pos, ref, alt
vcf = vcf.iloc[ :, [0,1,3,4]]

vcf = pandas.concat( [ vcf, gts], axis = 1)

# add column names
vcf.columns = ['chrom','pos','ref','alt'] + [ 'sample_'+str(x+1) for x in range(nsamples)]

## extract info 
ann = pandas.read_csv( fann, sep='\t')

# extract columns : chrom, pos, ref, allele, feature_id, aa change
ann = ann[['CHROM','POS','REF','Allele','Feature_ID','HGVS.p','apos']]

# remove prefix
ann['HGVS.p'] = ann['HGVS.p'].str.replace('p.','')

# split columns by aa pos
aas = ann['HGVS.p'].str.split('[0-9]+').apply( lambda x: pandas.Series(x))

ann = pandas.concat([ ann, aas], axis = 1)
ann.drop( 'HGVS.p', axis=1, inplace=True)

# add column names
ann.columns = ['chrom','pos','ref','alt','feature_id','aapos','refa','alta']

# reset index
vcf.set_index(['chrom','pos'], inplace=True)
#ann.set_index(['chrom','pos'], inplace=True)

def num_to_nuc(X):
	nucs = [ X.ref ] + X.alt.split(',')
	nixs = { str(x):i for x,i in enumerate(nucs)}
	out = pandas.Series( X.ref, index = ['ref'])
	Z = X.iloc[2:].apply( lambda x: '/'.join([ nixs[k] if k in nixs.keys() else '.' for k in x.split('/')]))
	out = pandas.concat([ out, Z]) 
	return out

vcf = vcf.apply( num_to_nuc, axis=1)

# group annotations by chrom, pos and gene
grps = ann.groupby(['chrom','pos','feature_id'])

def dict_nuc_aa(X):
	# ref nuc and aa
	refn = X.ref.tolist()[0]
	refa = X.refa.tolist()[0]
	# alt nucs and aas
	n2a = X[['alt','alta']].set_index('alt').to_dict()['alta']
	# add ref to above dict
	n2a[refn] = refa
	return n2a

grp_na = grps.apply(dict_nuc_aa)
grp_na = grp_na.reset_index()

# only keep groups with valid amino-acid mapping
grp_na = grp_na[ grp_na[0].apply( lambda x: all( isinstance(v,str) for v in x.values() ) )]

def nuc_to_aa(X,df=vcf):
	Y = df.loc[ (X.chrom, X.pos), :]
	X['ref'] = Y.ref
	out = Y.iloc[1:].apply( lambda x: '|'.join( [ x, 
		'/'.join( X[0][n] if n in X[0].keys() 
			else '.' for n in x.split('/'))]))
	out = pandas.concat( [ X[['chrom','pos','feature_id','ref']], out])
	return out

out = grp_na.apply( nuc_to_aa, axis=1)
out.to_csv( fout, sep='\t', index=False)
