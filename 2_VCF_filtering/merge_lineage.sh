#!/bin/bash
set -ueo pipefail

##############
#
# run bcftools isec and bcftools merge to merge 2 pops, require gzip & tbix pop files
# Author Tianpeng Wang 2020-11
#
#############

vcf_pop1=$1
vcf_pop2=$2

name_pop1=${vcf_pop1/.vcf.gz/}
name_pop2=${vcf_pop2/.vcf.gz/}

## run bcftools isec
bcftools isec -p dir ${vcf_pop1} ${vcf_pop2}

cd dir/

### basic stats
bcftools stats 0002.vcf > 0002.stats

## vcf_index.sh 0002
vcf_index.sh 0002.vcf
vcf_index.sh 0003.vcf

## bcftools merge 
bcftools merge -Oz -o ${name_pop1}_${name_pop2}.vcf.gz 0002.vcf.gz 0003.vcf.gz

## remove intermediate files
rm 0000.vcf; rm 0001.vcf; rm 0002.vcf; rm 0003.vcf

echo "isec & merge completed"
