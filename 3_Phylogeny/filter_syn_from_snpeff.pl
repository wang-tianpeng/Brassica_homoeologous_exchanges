#!/usr/bin/env perl

$usage=<<USAGE;
	perl $0 syn.vcf from.vcf
USAGE

die $usage if @ARGV==0;

open IN,$ARGV[0] or die;
while(<IN>){
	chomp;
	next if /^#/;
	@line=split;
	$hash{$line[0].$line[1]}=1;
}
close IN;


open IN,"$ARGV[1]"or die;
while(<IN>){
	chomp;
	if (/^#/){
		print"$_\n";
	}else{
		@F=split;
		print "$_\n"if (exists $hash{$F[0].$F[1]});
	}

}
e

