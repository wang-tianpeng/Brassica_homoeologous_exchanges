#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;

my $usage = <<USAGE;
Usage:
	perl $0 -ref ref_genome.fasta -chr chr_name.id -dir 4_gvcf/ -pre name_prefix -cpu <INT>(based on the chrNum. default:2)

This is Whole genome sequencing pipeline STEP2 for variants calling of multi-samples.The pipeline can be paralleled based on multiple chromosomes chr_name.id. It starts from the result directory "4_gvcf/" of WGSpipeline_step1.pl and depends on the following softwares that can be run directly in teminal.

1. gatk (v4.1.2)
2. ParaFly


							Version 1.0
							2021-1-1
							Bug contact : 'wangtianpeng19\@hotmail.com'

USAGE

die $usage if @ARGV==0;
my($ref_genome, $gvcf_dir, $cpu, $prefix,$chr);
GetOptions(
	"pre:s" => \$prefix,
	"ref:s" => \$ref_genome,
	"chr:s" => \$chr,
	"dir:s" => \$gvcf_dir,
	"cpu:i" => \$cpu,
);

$cpu =4 unless (defined $cpu);
my @sample_file = glob "$gvcf_dir*.gz";
my $sample_count=@sample_file;
my ($sample_cmd, $prefix_name, @chr_array);

##chr_array
open CHR,"$chr" or die "$!";
while(<CHR>){
	chomp;
	push @chr_array,$_;
}
close CHR;

## gvcf_name;
foreach my $sample (@sample_file){
	$sample_cmd.="-V "."$sample"." ";
}
$prefix_name=$prefix."_R"."$sample_count";


print "\n==========================================\n";
print "Step5: Combine multi HaplotypeCaller GVCF files into a single GVCF file\n\n";
print "Step6: Joint genotyping on a multi-sample GVCF file from CombineGVCFs\n\n";

foreach my $chr_name (@chr_array){
	my $command = "gatk CombineGVCFs -L $chr_name -R $ref_genome $sample_cmd-O 5_gvcf_all/$prefix_name.all.$chr_name.g.vcf.gz && gatk GenotypeGVCFs -R $ref_genome -V 5_gvcf_all/$prefix_name.all.$chr_name.g.vcf.gz -O 6_vcf_raw/$prefix_name\_raw.$chr_name.vcf.gz";
	open OUT,">>5_6_gvcfCG.cmd" or die "$!";
	print OUT "$command\n";
	close OUT;
}

unless (-e "00_logfile/5_6.gvcf.ok"){
	mkdir "5_gvcf_all" unless -e "5_gvcf_all";
	mkdir "6_vcf_raw" unless -e "6_vcf_raw";
	my $command="ParaFly -c 5_6_gvcfCG.cmd -CPU $cpu";
	system($command) == 0 or die "cannot execute :$command\n";
	open OUT, ">00_logfile/5_6.gvcf.ok" or die "$!";
	close OUT;
}else{print STDERR "CMD skipped: gatk CombineGVCFs && GenotypeGVCFs";}


print "\n==========================================\n";
print "Step7: Merge different chr based vcf files";

my @vcf_chr_file = glob "6_vcf_raw/*.gz";
my $vcf_chr_cmd;
foreach my $vcf_chr (@vcf_chr_file){
	$vcf_chr_cmd.="-I "."$vcf_chr"." ";
}

unless (-e "00_logfile/7.merge.ok"){
	my $command = "gatk MergeVcfs $vcf_chr_cmd-O 6_vcf_raw/$prefix_name\_raw.all.vcf.gz";
	system($command) == 0 or die "cannot execute:$command\n";
	open OUT, ">00_logfile/7.merge.ok" or die "$!";
	close OUT;
}else{print STDERR "CMD skipped: gatk mergeVcfs";}



print "\n=======================================\n";
print "Step8: Select SNP & InDel variants";
	
my $command = "gatk SelectVariants -select-type SNP -V 6_vcf_raw/$prefix_name\_raw.all.vcf.gz -O 6_vcf_raw/$prefix_name\_SNP.raw.vcf.gz";
print STDERR "$command\n";
system ($command) == 0 or die "cannot execute : $command\n";
my $command = "gatk SelectVariants -select-type INDEL -V 6_vcf_raw/$prefix_name\_raw.all.vcf.gz -O 6_vcf_raw/$prefix_name\_INDEL.raw.vcf.gz";
print STDERR "$command\n"; system($command) ==0;
open OUT, ">00_logfile/8.select_vcf.ok" or die "$!";
close OUT;


