#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2020-12-28)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tsv> <output_prefix>
Input: .tsv, 3+1+3 = 7 columns
	Col 1-3: junction's chr, pos, strand
	Col 4: read ID
	Col 5-7: bait-end's chr, pos, strand
Output:
	<output_prefix>.pos.bdg
	<output_prefix>.neg.bdg
".$version;

if(@ARGV < 2){
	die $usage;
}
my ($input_filename, $output_prefix) = @ARGV;

my %H = ();

my (@F, @B, @J);
open(IN, $input_filename) or die "Error: cannot open input file $input_filename\n";
while(<IN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	$H{join(":", @F[4..6])}++;
}
close(IN);


open(OUTpos, ">".$output_prefix.".pos.bdg") or die "Error: cannot open output file $output_prefix.pos.bdg\n";
open(OUTneg, ">".$output_prefix.".neg.bdg") or die "Error: cannot open output file $output_prefix.neg.bdg\n";
my (@A, $key);
foreach $key (sort {
			@A = split(/:/, $a);
			@B = split(/:/, $b);
			$A[0] cmp $B[0] || $A[1] <=> $B[1] || $A[2] <=> $B[2];
		} keys %H){
	@F=split(/:/, $key);
	if($F[2]>0){
		print OUTpos join("\t", $F[0], $F[1]-1, $F[1], $H{$key})."\n";
	}
	if($F[2]<0){
		print OUTneg join("\t", $F[0], $F[1]-1, $F[1], $H{$key})."\n";
	}
}
close(OUTpos);
close(OUTneg);

0;

