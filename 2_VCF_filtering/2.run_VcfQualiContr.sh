#!/bin/sh 
set -euo pipefail

## Run vcftools filter process based on run_VcfQualiContr.sh
## Author Tianpeng Wang 2020-11


# INPUT AND OUTPUT

VCF_IN= napus_final_R289_SNP.fil30.rm.als.recode.vcf
VCF_OUT=napus_final_R289_SNP.fil30.rm.als.quali085.recode.vcf.gz


##MAF 0.05
MAC=3
MISS=0.85
QUAL=30
MIN_DEPTH=2
MAX_DEPTH=40


vcftools --gzvcf $VCF_IN \
    --mac $MAC \
    --max-missing $MISS \
    --minQ $QUAL \
    --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
    --minDP $MIN_DEPTH --maxDP $MAX_DEPTH \
    --recode --stdout  2>$VCF_OUT.log | bgzip -c > \$VCF_OUT
    
