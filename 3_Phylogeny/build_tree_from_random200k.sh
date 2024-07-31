#!/bin/bash
set -euo pipefail

scripts=$HOME/scripts/8_phylogeny
file=lineageC407_nigra2.ale.vcf.gz
name=${file/.vcf.gz/}

plink --vcf ${file} --geno 0.1 --make-bed --out ${name} --aec --double-id

plink --bfile ${name} --thin-count 200000 --recode vcf-iid --out  ${name}.200k --aec --double-id --set-missing-var-ids @__#

${scripts}/vcf2phylip.py -i ${name}.200k.vcf -m 10 -f

~/miniconda3/envs/wgs/bin/iqtree -s ${name}.200k.min10.fasta -st DNA -m TVM+R6 -nt 30 -bb 1000 -pre ${name}.200k
