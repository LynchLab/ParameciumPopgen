#!/bin/bash

WDR=$1 # working directory
FBED=$2 # input reference features bed file

## other arguments
MINQ=25 # only consider reads with minq mapping quality
MIND=3 # to report no. of bases with depth at least minc

## these files must already be present
FBAM=$WDR/bamlist_for_bed.txt # file listing BAM files

## these files are created in the analysis
DPAC=$WDR/gene_pa_by_cov # directory with sample directories for bedcov analysis
DGP=$WDR/genes_present # directory of lists of detected genes for each sample
FOUT=$WDR/gene_sample_bcov.tsv # output file with breadth of coverage for genes across samples

echo "Start: `date`"

## create directories
mkdir -p $DPAC
mkdir -p $DGP

# filter reads and calculate coverage for features in the reference BED file
while read fin; do 
	# sample name
	SAMPLE=$(basename $fin .bam); 
	# sample directory
	SDR=$DPAC/$SAMPLE;
	mkdir -p $SDR;
	# cov bed output file
	FCB=$SDR/cov.bed;
	if [ ! -s $FCB ]; then
		samtools bedcov -Q $MINQ -G UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP -j -d $MIND $FBED $fin > $FCB;
	fi;
	## gene presence inference
	FPT=$SDR/genes_present.txt;
	[ ! -s $FPT ] && Rscript infer_gene_presence.R $SDR;
	# save genes present file renamed by the sample name to a separate dir
	[ ! -L $DGP/$SAMPLE ] && ln -rs $FPT $DGP/$SAMPLE;
done < $FBAM

# pool results from across files into a single table
[ ! -s $FOUT ] && python tab_sample_fsitecov.py $DGP $FOUT
conda deactivate
echo "Finished: `date`"
