#!/usr/bin/env perl
use strict;

my $usage=<<USAGE;
	perl $0 napus_zs11.vcf random_sample ref 
USAGE

die $usage if @ARGV==0;

my ( $sample, $ref);
$sample=$ARGV[1];
$ref=$ARGV[2];

$ARGV[0]=~/(.*)\.vcf.gz/;#print"$ARGV[0]\n";
my $name=$1;
my $command = "vcftools --gzvcf $ARGV[0] --indv $sample --recode --out snpeff_$sample 2>snpeff_$sample.log";
#print "$command\n";
system($command)==0 or die "$!";
my $command = "java -Xmx100G -jar /data/mg1/wangtp/software/snpEff/snpEff.jar -c /data/mg1/wangtp/software/snpEff/snpEff.config $ref snpeff_$sample.recode.vcf > $name.snp.eff.vcf -csvStats $name.csv -stats $name.html";
#print "$command\n";
system($command)==0 or die "$!";
print "snpEFF annotation complete\n";


