#!/usr/bin/env perl -w
## get_reciprocalMatch_synorth.pl

open IN1, "$ARGV[0]" or die;
open IN2, "$ARGV[1]" or die;

while(<IN1>){
	chomp;
	@F=split;
	$line="$F[4]_$F[0]";
	$hash{$line}=1;
}

while(<IN2>){
	chomp;
	@F=split;
	print"$_\n" if (exists $hash{"$F[0]_$F[4]"});
}
