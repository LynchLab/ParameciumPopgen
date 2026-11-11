#!/bin/bash

## arguments
FCDS=$1 # input CDS FROM GENOMIC
ODR=$2  # output directory of degeneracy sequences
CTN=6 	# codon table

NCPU=16  # number of CPU cores

echo "Start - `date`"

## main
if [ ! -d $ODR ]; then
    echo "Sequence of degeneracy for each gene"
    python site_degeneracy.py $FCDS $ODR $CTN $NCPU
fi

echo "Finish - `date`"
conda deactivate
