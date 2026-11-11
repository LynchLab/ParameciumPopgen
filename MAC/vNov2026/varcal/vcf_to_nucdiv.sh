#!/bin/bash

## arguments
INDR=$1 # input dir path
REFDR=$2 # reference dir path

NCPU=4

## files that must already be present ######
FVCF=$INDR/snps_pass.vcf # input vcf
FANN=$INDR/snps_pass_ann.tsv # tsv file of all annotations 
FBCV=$INDR/gene_sample_bcov.tsv # file of genes' breadth of coverage across samples
SDDR=$REFDR/sitedeg # ref site deg dir path
FRNS=$REFDR/syn_nsyn_sitecount.csv # ref syn nsyn site count table

########################

## these files are created by the program ########
FSNS=$INDR/snps_syn_nsyn.tsv # tsv file of syn/non-syn annotations
FAGT=$INDR/snps_aa_gts.tsv # amino-acid genotypes
DDM4=$INDR/gene_pdmat_deg4 # output dir of distance matrices based on silent sites
FMD4=$INDR/mean_pdmat_deg4.csv # output file of mean pairwise distances over genes
DSNS=$INDR/pairwise_syn_nsyn_difs # dir of pairwise SYN NON-SYN dif counts for genes
FNSP=$INDR/gene_piN_piS_pairwise.csv # output file of genes pN/pS
#########################

echo "Start - `date`"

## main
# subset of annotated mutations that are synonymous/missense
if [ ! -s $FSNS ]; then
	head -n 1 $FANN > $FSNS
	awk -F '\t' '$5=="synonymous_variant" || $5=="missense_variant"' $FANN >> $FSNS
fi

# amino acid genotypes of above mutations
[ ! -s $FAGT ] && python vars_samples_aa.py $FVCF $FSNS $FAGT

# get matrices of pairwise distances based on 4-fold degenerate sites for each gene
[ ! -d $DDM4 ] && python pwdist_gene_d4.py $FAGT $DDM4 $FBCV $FSNS $SDDR 4 $NCPU
# take mean pairwise distance over genes
[ ! -s $FMD4 ] && python mean_pdmat.py $DDM4 $FMD4 $NCPU

# number of pairwise syn and non-syn differences by genes
[ ! -d $DSNS ] && python gene_pw_syn_nsyn_dif_vcf.py $FAGT $DSNS $NCPU
# pairwise piN/piS
[ ! -s $FNSP ] && python gene_piN_piS_pairwise.py $DSNS $FRNS $FBCV $FNSP $NCPU

conda deactivate

echo "Finish - `date`"
