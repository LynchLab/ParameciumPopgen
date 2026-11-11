#!/bin/bash

## arguments
FVCF=$1 # input vcf
REF=$2 # snpeff database

FANN=${FVCF/.vcf/_ann.vcf}
FTAB=${FANN/.vcf/.tsv}

echo "Start - `date`"
## main

# run VCF against the database
[ ! -s $FANN ] && java -Xmx8g -jar ~/opt/snpEff/snpEff.jar $REF $FVCF > $FANN
# extract vcf annotation into dataframe
[ ! -s $FTAB ] && python extract_vcf_annot.py $FANN $FTAB

echo "Finish - `date`"
conda deactivate
