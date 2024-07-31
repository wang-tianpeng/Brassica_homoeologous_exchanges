file=napus_final_R289_SNP.fil30.rm.als.quali08.recode_rapa_zs11_A_R199.SNP.raw.vcf.gz
file_name=${file/.vcf.gz/}

vcftools --gzvcf ${file} --bed zs11A_synSites.bed --remove LineageA_remove.id --recode -c 2>${file_name}.log |gzip -c >${file_name}.syn.vcf.gz

nohup vcftools --gzvcf napus_final_R289_SNP.fil30.rm.als.quali08.recode_ole_zs11C_R119_SNP.raw.vcf.gz --bed zs11C_synSites.bed --remove LineageC_remove.id --recode -c 2>napus_final_R289_SNP.fil30.rm.als.quali08.recode_ole_zs11C_R119_SNP.raw.syn.log |bgzip -c > napus_final_R289_SNP.fil30.rm.als.quali08.recode_ole_zs11C_R119_SNP.raw.syn.vcf.gz &

