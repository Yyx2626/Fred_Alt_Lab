#!/usr/bin/perl

use strict;
use warnings;

my $version = "
Author: Adam Yongxin Ye @ BCH
Version: 0.1.1 (2020-02-08)
";
my $usage = "Usage: $0 <c|u|v> <input.tlx> <ref.bed> <output.bed>
Input:	<input.tlx>
	<ref.bed>
Output: <output.bed>
".$version;

if(@ARGV < 4){
	die $usage;
}
my ($opt, $tlxfile, $bedfile, $output) = @ARGV;
if($opt ne "c" && $opt ne "u" && $opt ne "v"){
	die "<c|u|v> is not specified correctly
	c	count
	u	extract
	v	extract complement
";
}

my $cmd1 = "perl tlx2BED.pl $tlxfile /dev/stdout ";

my $cmd2 = "bedtools intersect -$opt  -a <( $cmd1 )  -b $bedfile  >$output";
if($opt eq "c"){
	$cmd2 = "bedtools intersect -$opt  -b <( $cmd1 )  -a $bedfile  >$output";
}

print STDERR "[CMD] $cmd2\n";
system("bash", "-c", $cmd2);
