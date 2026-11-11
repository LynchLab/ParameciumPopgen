#!/bin/bash
bcftools filter -i 'QUAL>=200 && INFO/DP<=1000 && RPB>1e-5 && MQB>1e-5 && BQB>1e-5' snps.bcf | bcftools filter -i 'FMT/SP<=2' > snps_pass.vcf
