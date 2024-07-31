#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
Usage:
	perl $0 -ref ref_genome.fasta -1 left_1.fas -2 right_2.fas -cpu <INT>(default:4)

This is Whole Genome Sequencing Upstream analysis pipeline with resuming function. This pipeline depends on the following softwares that can be run directly in teminal.

1. fastp (v0.20)
2. bwa (v0.7.17-r1188)
3. samtools (v1.9)
4. sambamba (0.7.0)
5. gatk (v4.1.2)
								version 0.3
								2019-11-21	
								Bug contact:'wangtianpeng19\@hotmail.com'

USAGE

die $usage if @ARGV==0;
my ($ref_genome,$left,$right,$cpu,$ref_genome_prefix);
GetOptions(
	"ref:s" => \$ref_genome,
	"1:s" => \$left,
	"2:s" => \$right,
	"cpu:i" => \$cpu,
);

# Preliminary Test : Detecting the softwares dependency
print STDERR "\nDetecting the software dependency\n\n";
##fastp
my $software= `fastp -h 2>&1`;
if($software=~/usage:/){
	print STDERR "fastp:\tOK\n";
}else{
	die "fastp\tfailed\n";
}

##bwa
my $software=`bwa 2>&1`;
if($software =~/Version/){
	print STDERR "bwa:\tOK\n";
}else{
	die "bwa\tfailed\n";
}

##samtools
my $software= `samtools 2>&1`;
if($software=~/Program: /){
	print STDERR "samtools\tOK\n";
}else{
	die "samtools\tfailed\n";
}

##sambamba
my $software= `sambamba 2>&1`;
if($software=~/Usage:/){
	print STDERR "sambamba\tOK\n";
}else{
	die "sambamba\tfailed\n\n";
}

##gatk4
my $software= `gatk -h`;
if ($software=~/Usage/){
	print STDERR "gatk\tOK\n";
}else{
	die "gatk\tfailed. GATK4 is recommended to be installed by conda.\n\n";
}

## step0 Detecting the reference genome bwa-mem Index
print "\n====================================\n";
print "Step0 Detecting Reference Genome BWA index files\n";
my $ref_dir=dirname($ref_genome);
unless(-e "$ref_genome.bwt"){
	my $command = "bwa index $ref_genome";
	print STDERR (localtime).": $command\n";
	system($command) == 0 or die "can not execute : $command\n";
#	print "$ref_genome.bwt BWA index is not ok\n";
}else{ print "Ready for $ref_genome BWA index\n";}
## step0_plus Detecting the GATK4 index for Reference genome
print "\n====================================\n";
print "Detecting gatk4 index for refence genome";
$ref_genome=~/(.*)\./;$ref_genome_prefix=$1;
unless (-e "$ref_genome_prefix.dict" && "$ref_genome.fai"){
	my $command = "gatk CreateSequenceDictionary -R $ref_genome";
	print STDERR (localtime).": $command\n";
	system ($command) == 0 or die "can not execute : $command\n";
	my $command = "samtools faidx $ref_genome";
	system ($command)==0 or die "failed to execute $command\n";
}else{ print "Ready for $ref_genome GATK dictionary index and Samtools index";}
## Detecting the samtools faidx for reference genome
print "\n====================================\n";
print "Detecting samtools faidx for reference genome";
unless(-e "$ref_genome.fai"){
	my $command = "samtools faidx $ref_genome";
	print STDERR (localtime).": $command\n";
	system($command) == 0 or die "can not execute : $command\n";
}else{
	print "Ready for $ref_genome samtools faidx\n";
}


$cpu = 4 unless (defined $cpu);
my $sample_dir=dirname($left);
basename($left)=~/(.*)_1\.(.*?)$/;
my $sample=$1; my $suffix=$2;

mkdir "00_logfile" unless -e "00_logfile";

## step1 From raw reads To clean reads; 
print "\n====================================\n"; 
print "Step1 Quality control and trimming using fastp\n"; 
mkdir "1_clean_data" unless -e "1_clean_data";
unless (-e "00_logfile/1.fastp.$sample.ok"){
	my $command = "fastp -i $sample_dir/$sample\_1.$suffix -I $sample_dir/$sample\_2.$suffix -o 1_clean_data/$sample\_1.clean.fq.gz -O 1_clean_data/$sample\_2.clean.fq.gz -z 4 -q 20 -u 30 -n 6 -c -w $cpu -h 1_clean_data/$sample.clean.html 2>00_logfile/logfile.$sample.log";
	print STDERR (localtime).":  $command\n\n";
	system($command) == 0 or die "can not execute : $command\n";
	open OUT, ">00_logfile/1.fastp.$sample.ok" or die "$!";
	close OUT;
}else{
	print STDERR "CMD(skipped): fastp  \n";
}

##step2 From clean reads to sorted bam;
print "\n====================================\n"; 
print "Step2 Map clean reads to RefGenome & Sort the bam file\n";
mkdir "2_map" unless -e "2_map";
unless (-e "00_logfile/2.map.$sample.ok"){
	my $command = "bwa mem -t $cpu -R '\@RG\\tID:$sample\\tLB:$sample\\tPL:ILLUMINA\\tSM:$sample' $ref_genome 1_clean_data/$sample\_1.clean.fq.gz 1_clean_data/$sample\_2.clean.fq.gz |samtools sort -@ $cpu -m 4G >2_map/$sample.sort.bam 2>>00_logfile/logfile.$sample.log";
	print STDERR (localtime).":  $command\n";
	system($command)==0 or die "cant execute :$command\n\n";
	open OUT, ">00_logfile/2.map.$sample.ok" or die "$!";
	close OUT
}else{
	print STDERR "CMD(skipped): bwa  \n";
}


##step3 From sorted bam to markdup bam
print "\n====================================\n"; 
print "Step3 Mark duplicates & Index bam\n";
mkdir "3_markdup" unless -e "3_markdup";
unless (-e "00_logfile/3.dup.$sample.ok"){
	my $command = "sambamba markdup --hash-table-size 786432 --overflow-list-size 600000 --tmpdir=. --sort-buffer-size=204800 -r -t $cpu 2_map/$sample.sort.bam 3_markdup/$sample.sort.markdup.bam 2>>00_logfile/logfile.$sample.log";
	print STDERR (localtime).": $command\n";
	system($command) == 0 or die "cannot execute: $command\n\n";
	open OUT, ">00_logfile/3.dup.$sample.ok" or die "$!";
	close OUT
}else{
	print STDERR "CMD(skipped): sambamba \n";
}

## STEP4 Get individual gVCF by using GATK haplotypeCaller
print "\n====================================\n"; 
print "STEP4 Get individual gVCF files by using GATK haplotypecaller\n";
mkdir "4_gvcf" unless -e "4_gvcf";
unless (-e "00_logfile/4_gvcf.$sample.ok"){
	my $command = "gatk HaplotypeCaller --emit-ref-confidence GVCF -R $ref_genome -I 3_markdup/$sample.sort.markdup.bam -O 4_gvcf/$sample.g.vcf.gz";
	print STDERR (localtime).": $command\n";
	system($command) == 0 or die "cannot execute : $command\n\n";
	open OUT, ">00_logfile/4.gvcf.$sample.ok" or die "$!";
	close OUT
}else{ print STDERR "CMD (skipped): gatk \n";}


#mkdir "$out_prefix.tmp" unless -e "$out_prefix.tmp";
#chdir "$out_prefix.tmp"



