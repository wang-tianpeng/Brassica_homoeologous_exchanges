#!/bin/bash 
set -euo pipefail

##############
#
# run_geno_fromExistVCF.sh $gvcf $exist_VCF $ref_genome
#
#############

gvcf=$1
exist_vcf=$2
ref_genome=$3

## new accessions 
nigra_name=${gvcf/.g.vcf.gz/}
out_name=${exist_vcf/.vcf.gz/}_${gvcf/.g.vcf.gz/}

## run gatk GenotypeGVCFs
gatk GenotypeGVCFs -V ${gvcf} -L ${exist_vcf} -all-sites true -R ${ref_genome} -O ${nigra_name}.select.vcf.gz

## run merge 
vcf-merge ${nigra_name}.select.vcf.gz ${exist_vcf} 2>${out_name}.log | bgzip -c > ${out_name}.lineageA.vcf.gz

## Only remain bi-allele SNPs 
vcftools --gzvcf ${out_name}.lineageA.vcf.gz --min-alleles 2 --max-alleles 2 --recode -c 2>>${out_name}.log | bgzip -c > ${out_name}.lineageA.ale.vcf.gz
