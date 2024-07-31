#!/usr/bin/env perl -w
## perl sort_by_reciprocalMatch.pl Brapa_genes_only.cov Bnapus_reciprocalMatch.txt

open IN1, "$ARGV[0]" or die;
open IN2, "$ARGV[1]" or die;

while (<IN1>){
	chomp;
	@F=split;
	$hash{$F[3]}=$_;
}

while (<IN2>){
	chomp;
	@F=split;
	if (exists $hash{$F[1]} ){
		print "$F[0]\t$F[1]\t$hash{$F[1]}\n";
	}elsif (exists $hash{$F[0]}){
		print "$F[0]\t$F[1]\t$hash{$F[0]}\n";
	}
}
