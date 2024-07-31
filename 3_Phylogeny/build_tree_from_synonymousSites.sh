#!/bin/bash
set -euo pipefail

#####
#
# filter Synonymous mutations & 4d sites & LD pruned sites for further phylogeny analysis
# snpeff config ref : napus_zs11
# Author : Tianpeng Wang
# 
#####

scripts=$HOME/scripts/8_phylogeny
pop=lineageA484_nigra2.ale.vcf.gz
name=${pop/.vcf.gz/}
sample=20N_wtp_017

### snp_annotation
perl ${scripts}/snpeff_annotation.pl ${pop} ${sample} napus_zs11

### filter synonymous mutants
perl -lne 'if(/^#/){print}else{print if $_=~/synonymous_variant/}' $name.snp.eff.vcf > $name.snp.syn.vcf

### Get the final Synonymous mutation sites.
zcat ${pop} | perl ${scripts}/filter_syn_from_snpeff.pl ${name}.snp.syn.vcf - > ${name}.snp.syn.final.vcf

### Convert the final vcf to the fasta file
${scripts}/vcf2phylip.py -i ${name}.snp.syn.final.vcf -f -m 10

### Construct ML phylogeny tree using iqtree with the BIC model set as TVM+R6
~/miniconda3/envs/wgs/bin/iqtree -s ${name}.snp.syn.final.min10.fasta -st DNA -m TVM+R6 -nt 30 -bb 1000 -pre ${name}.snp.syn.final &

## finally we get the phylogeny treefile. After carefully selection and polishment  of the tree file, we will use R package ggtree to visualise and plot.
