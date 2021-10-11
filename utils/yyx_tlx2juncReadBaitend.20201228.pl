#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2020-12-28)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.tlx>
Output: STDOUT (.tsv), 3+1+3 = 7 columns
	Col 1-3: junction's chr, pos, strand
	Col 4: read ID
	Col 5-7: bait-end's chr, pos, strand
".$version;

if(@ARGV < 1){
	die $usage;
}
my ($input_filename) = @ARGV;

my (@F, @B, @J);
open(IN, $input_filename) or die "Error: cannot open input file $input_filename\n";
while(<IN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	@B = @F[7..10];
	@J = @F[2..4];
	if($ B[3]>=0){
		@B = @B[0,2,3];
	}else{
		@B = @B[0,1,3];
	}
	print join("\t", @J, $F[0], @B)."\n";
}
close(IN);

0;

