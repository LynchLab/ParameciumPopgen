### Get degeneracy of each site from NCBI reference assembly CDS file

## require
import os, sys, multiprocessing
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id as codon_table
from collections import Counter

## arguments
fin = sys.argv[1] # input cds from genomic
dout = sys.argv[2] # output directory of site-degeneracy tables
ctn = int(sys.argv[3]) # ncbi codon table number
ncpu = int(sys.argv[4]) # number of CPUs

# function to convert a list of nucleotides to a list of codons
n2c = lambda seq: [ ''.join(seq[x:x+3]) for x in range(0,len(seq),3)]

# create output dir if needed
not os.path.exists(dout) and os.mkdir(dout)

## Function ##
## Function to get degeneracy of each site in a codon
def sitedegen( ct=codon_table[ctn], nucs = 'ACGT'):
    """
    Get n-fold degeneracy of each site of every codon in a given codon table
    Arguments:
    - ct 	- codon table
    - nucs 	- string of allowed nucleotides ['ACGT']
    """
    ixs = range(3)
    table = ct.forward_table
    stops = ct.stop_codons
    codons = table.keys()
    out = { k:() for k in codons}

    # for each codon in the table
    for codon in codons:
        # corresponding amino acid
        aa = table[codon]
        # initialize list to hold fraction of synonymous changes at each site of a codon with zeros
        f = [ 0 for _ in ixs]
        # for every site in the codon
        for i in ixs:
            # produce all possible codons modifying this site
            altc = [ ''.join( [ n if x == i else codon[x] for x in ixs]) for n in nucs]
            # corresponding amino acids for above codons
            alta = [ table[c] for c in altc if c not in stops]
            # degeneracy = # codons coding for most frequent amino acid
            f[i] = Counter(alta).most_common(1)[0][1]
        out[codon] = f
    return out

# dict of codons and their degeneracy for each site
degen_tab = sitedegen()

# load records
recs = list( SeqIO.parse( fin, 'fasta'))

def get_dseq( rec, dtab=degen_tab, dr=dout):
    # gene id
    gid = rec.id
    # output path
    fout = os.path.join( dr, gid+'.txt')
    if os.path.exists(fout):
        return

    # split into codons
    codons = n2c(rec.seq)
    dlist = [ dtab[i] if i in dtab.keys() else list('000') for i in codons]
    dseq = ''.join( str(i) for x in dlist for i in x)
    with open(fout,'w') as flob:
    	flob.write(dseq+'\n')
    return

## Main
with multiprocessing.Pool(ncpu) as mp_pool:
	mp_pool.map( get_dseq, recs)
