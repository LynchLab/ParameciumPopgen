#!/bin/bash

WDR=$1

cd $WDR

# no. of samples
N=`bcftools query -l raw.bcf | wc -l`

# -g remove snps within 3 bp of indels, -G remove indels separated by less than 9 bp
bcftools filter -g3 -G9 raw.bcf -Ob -o variants_wo_clusters.bcf

# separate snps and indels
bcftools view -v snps variants_wo_clusters.bcf -Ob -o snps.bcf
bcftools view -v indels variants_wo_clusters.bcf -Ob -o indels.bcf

# header for quality stats table
lns=`for((i=0;i<$N;i++)); do echo "\tSP$i"; done`; lns=`echo $lns | sed 's/ //g'`; echo -e "QUAL\tDP\tMQB\tBQB\tRPB"$lns > snps_qual_stats.txt

# extract quality stats
bcftools query -f '%QUAL\t%DP\t%MQB\t%BQB\t%RPB[\t%SP]\n' snps.bcf >> snps_qual_stats.txt

# plot stats drb
Rscript ~/repos/pge_tol/plot_vcf_qual_stats_drbs.R snps_qual_stats.txt

cd $HOME
