#!/bin/bash

## arguments
BAMLS=$1 # list of bam files to be used
FREF=$2 # reference fasta file 
WDR=$3 # working and output directory
MAXD=$4 # maximum depth ( > 2 x average sample depth of cov)

## aurelias treated as haploid
PLOIDY=1
# mapping quality
MAPQ=25
# base quality
BASEQ=30
# number of threads
NCPU=8 

echo -e "Input Arguments:-\n\tREF: $FREF\n\tWDR: $WDR\n\tBAMLS: $BAMLS\n\tMAXD: $MAXD\n\tBASEQ: $BASEQ\n\tMAPQ: $MAPQ\n\tPLOIDY: $PLOIDY\n"

## main
echo "`date`: Start.."

# index reference
if [ ! -s $FREF".fai" ]; then
    echo "Indexing reference"
    samtools faidx $FREF
fi

# calling variants
if [ ! -s $WDR/raw.bcf ]; then
    echo "running mpileup"
    bcftools mpileup --excl-flags UNMAP,SECONDARY,SUPPLEMENTARY,QCFAIL,DUP --redo-BAQ --min-BQ $BASEQ --min-MQ $MAPQ --per-sample-mF --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --threads $NCPU -d $MAXD -b $BAMLS -f $FREF -Ou | bcftools call --ploidy $PLOIDY -mv --threads $NCPU -Ob -o $WDR/raw.bcf
else
    echo "RAW BCF exists!"
fi


echo "`date`: Finished!"
