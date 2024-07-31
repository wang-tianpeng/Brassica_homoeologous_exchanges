#!/bin/sh 
set -euo pipefail

name=${1/.vcf.gz/}
mkdir ${name}_out
vcftools --gzvcf $1 --freq2 --out ${name}_out/$name &
vcftools --gzvcf $1 --depth --out ${name}_out/$name &
vcftools --gzvcf $1 --site-mean-depth --out ${name}_out/$name &
vcftools --gzvcf $1 --site-quality --out ${name}_out/$name &
vcftools --gzvcf $1 --missing-indv --out ${name}_out/$name &
vcftools --gzvcf $1 --missing-site --out ${name}_out/$name &
vcftools --gzvcf $1 --het --out ${name}_out/$name &
