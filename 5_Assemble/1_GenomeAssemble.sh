#!/bin/sh 
set -euo pipefail

### 1. files preparation, including:

ONTreads=swede.ont.raw.merge.fastq.gz


### 2. Use canu to get the contig assembly
canu -d swede_canu_out_raw -p swede_canu_raw genomeSize=1g maxThreads=80 -nanopore-raw $ONTreads


### 3. Polish the raw contig genome with short reads
### by NextPolish
vim run.cfg sgs.fofn lgs.fofn
nohup nextPolish run.cfg 


# 4. build the chromosomes based on HiC data using ALLhic
## step0. prepare for the reference index
samtools faidx draft.asm.fasta 
bwa index -a bwtsw draft.asm.fasta

## step1. mapping the raw hic data
bwa aln -t 24 -f reads_R1.sai draft.asm.fasta reads_R1.fastq.gz
bwa aln -t 24 -f reads_R2.sai draft.asm.fasta reads_R2.fastq.gz
bwa sampe -f sample.bwa_aln.sam draft.asm.fasta reads_R1.sai reads_R2.sai reads_R1.fastq.gz reads_R2.fastq.gz
    ### Time consuming for reads data


## step2. remove the low-quality mapping reads, remain the uniq-paired reads
PreprocessSAMs.pl sample.bwa_aln.sam draft.asm.fasta MBOI &
        ## output: sample.bwa_aln.REduced.paired_only.bam
filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam
samtools view -bt draft.asm.fasta.fai sample.clean.sam > sample.clean.bam

## step3 build allele.ctg.table
ALLHiC_prune -i Allele.ctg.table -b sample.clean.bam -r draft.asm.fasta

# step4. group the contig based on K number
ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e AAGCTT -k 16  

# step5. rescue the unmapped contigs into  the groups
ALLHiC_rescue -b sample.clean.bam -r draft.asm.fasta -c prunning.clusters.txt -i prunning.counts_AAGCTT.txt

# step6. optimise and reorder the contigs in each group.
allhic extract sample.clean.bam draft.asm.fasta --RE AAGCTT
  
for i in group*.txt; do 
    allhic optimize $i sample.clean.clm &
done
      

# step7. finalise and convert tour files into the final HiC-adjusted fasta  
ALLHiC_build draft.asm.fasta
        # output the groups.asm.fasta; groups.agp

# step8. plot the heatmap
samtools faidx groups.asm.fasta
cut -f 1,2 groups.asm.fasta.fai > chrn.list
ALLHiC_plot -b sample.clean.bam -a groups.agp -l chrn.list -s 500k -o pdf

