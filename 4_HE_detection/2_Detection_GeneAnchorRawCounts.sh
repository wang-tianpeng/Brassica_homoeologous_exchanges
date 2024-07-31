#!/bin/sh 
set -euo pipefail

### 1. files preparation, including:

### 1.1 Brassica napus accessions Bam files after removing duplicates 
DirBam=../3_markdup

### 1.2 Bra and Bole Bed files
BedBra=Brapa_genes_only.bed
BedBol=Boleracea_genes_only.bed

### 1.3 reciprocal syntenic gene pair list
HeAnchor=Bnapus_ole2rapa_synortho.reciprocalmatch.txt


### 2. Use Samtools bedcov to calculate the raw counts of each homoeologous gene pair anchor
ls -1 ../3_markdup/*.id | cat |perl -pe 's/.*?dup\///; s/.markdup.bam//' > id.he.bnapus.acc

samtools bedcov $BedBra $DirBam/*.markdup.bam > Brapa_genes_only.cov
perl sort_by_reciprocalMatch.pl Brapa_genes_only.cov $HeAnchor >BnapusA_reciprocalMatch.cov

samtools bedcov $BedBol $DirBam/*.markdup.bam > Boleracea_genes_only.cov
perl sort_by_reciprocalMatch.pl Boleracea_genes_only.cov $HeAnchor >BnapusC_reciprocalMatch.cov


### Output the raw counts for both A and C homoeologous genes
### id.he.bnapus.acc; BnapusA_reciprocalMatch.cov; BnapusC_reciprocalMatch.cov