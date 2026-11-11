#!/bin/bash

## arguments
FRD=$1 # read file path prefix
FREF=$2 # reference genome file path
POUT=$3 # output files path prefix

NCPU=8

## main
echo "Start - `date`"

# create output directory if missing
mkdir -p $(dirname $POUT)

# index reference
if [ ! -s $FREF".bwt" ]; then
    echo "Indexing reference"
    bwa index $FREF
else
    echo "Index exists! Moving to next step.."
fi

# generate alignment, sort and index bam
if [ ! -s $POUT"_sorted.bam" ]; then
    echo "Aligning reads"
    bwa mem -t $NCPU $FREF $FRD"_1.fastq.gz" $FRD"_2.fastq.gz" | samtools sort -@ $NCPU -o $POUT"_sorted.bam" -
    samtools index -@ $NCPU $POUT"_sorted.bam"
else
    echo "Alignment exists!"
fi

# flagstat
if [ ! -s $POUT".flagstat" ]; then                                             
    echo "getting samtools FLAG statistics"                                                              
    samtools flagstat $POUT"_sorted.bam" > $POUT".flagstat"                              
else                                                                            
    echo "Flagstat file exists!"                                             
fi

# get coverage over reference
if [ ! -s $POUT"_cov.txt" ]; then
    echo "Sequencing Depth for every reference position"
    samtools depth $POUT"_sorted.bam" > $POUT"_cov.txt"
else
    echo "Coverage file exists!"
fi

echo "Finish - `date`"
