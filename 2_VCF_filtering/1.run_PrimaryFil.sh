#!/bin/sh 
set -euo pipefail

scripts=$HOME/scripts/7_vcf_filter/
file=napus_final_R289_SNP

gatk VariantFiltration \
	-V 6_vcf_raw/${file}.raw.vcf.gz \
	-O 7_snp_filter/${file}.fil30.vcf.gz \
	--filter-name FilterQual --filter-expression "QUAL < 30.0" \
	--filter-name FilterQD --filter-expression "QD < 2.0" \
	--filter-name FilterMQ --filter-expression "MQ < 30.0" \
	--filter-name FilterFS --filter-expression "FS > 60.0" 

vcftools --gzvcf 7_snp_filter/${file}.fil30.vcf.gz --remove-filtered-all --recode-INFO-all --recode \
	--min-alleles 2 --max-alleles 2 \
	 --stdout 2>7_snp_filter/${file}.fil30.rm.als.log | bgzip -c -@ 5 >7_snp_filter/${file}.fil30.rm.als.vcf.gz

cd 7_snp_filter/
tabix -p vcf -f ${file}.fil30.rm.als.vcf.gz
${scripts}cmd.VcfQualiContr.sh ${file}.fil30.rm.als.vcf.gz

cd ${file}.fil30.rm.als.recode_out/
Rscript ${scripts}plot_VcfQual.r -i ${file}.fil30.rm.als.recode
