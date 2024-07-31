chmod 755 run_geno_fromExistVCF.sh
nohup ./run_geno_fromExistVCF.sh nigra2.all.g.vcf.gz napus_final_R289_SNP.fil30.rm.als.quali08.recode_rapa_zs11_A_R199.SNP.raw.syn.vcf.gz ~/Database/brassica_napus/zs11.psudoA.genome.fa &

nohup ./run_geno_fromExistVCF.sh nigra2.all.g.vcf.gz napus_final_R289_SNP.fil30.rm.als.quali08.recode_ole_zs11C_R119_SNP.raw.syn.vcf.gz ~/Database/brassica_napus/zs11.psudoC.genome.fa &
