#!/usr/bin/perl

use strict;
use warnings;

my $version = "
Author: Adam Yongxin Ye @ BCH
Version: 0.2.0 (2020-03-09)
";
my $usage = "Usage: $0 <c|u|v> <yyx_tlx2bed.py> <input.tlx> <should_strand> <ref.bed> <output_prefix>
Input:	<input.tlx>
	<ref.bed>
Output:
	<output.tlx.bed>
	[output.tlx.pos.bed]
	[output.tlx.neg.bed]
	<output.intersect.bed>
	[output.intersect.pos.bed]
	[output.intersect.neg.bed]
".$version;

if(@ARGV < 6){
	die $usage;
}
my ($opt, $yyx_tlx2bed_py, $tlxfile, $should_strand, $bedfile, $output_prefix) = @ARGV;
if($opt ne "c" && $opt ne "u" && $opt ne "v"){
	die "<c|u|v> is not specified correctly
	c	count
	u	extract
	v	extract complement
";
}
#if($strand == 0){
#	if($strand eq "+" || $strand eq "-" || $strand eq "*"){
#		# good
#	}elsif($strand eq "0"){
#		$strand = "*";
#	}else{
#		die "Error: cannot recognize strand=$strand, which should be one of + or >0, - or <0, * or 0\n";
#	}
#}elsif($strand > 0){
#	$strand = "+";
#}else{
#	$strand = "-";
#}


my @filenames = ();
$filenames[1] = $output_prefix . ".tlx.bed";
my $command = "cat $tlxfile | python $yyx_tlx2bed_py >".$filenames[1];
print STDERR "[CMD] $command\n";
system("bash", "-c", $command);

$filenames[6] = $output_prefix . ".intersect.bed";
$command = "bedtools intersect -$opt  -a ".$filenames[1]."  -b $bedfile  >".$filenames[6];
if($opt eq "c"){
	$command = "bedtools intersect -$opt  -b ".$filenames[1]."  -a $bedfile  >".$filenames[6];
}
print STDERR "[CMD] $command\n";
system("bash", "-c", $command);

if($should_strand){
	$filenames[2] = $output_prefix . ".tlx.pos.bed";
	$command = "cat ".$filenames[1].' | perl -ne \'@F=split/\t/; if($F[5]=~/[+]/){ print; }\' >'.$filenames[2];
	print STDERR "[CMD] $command\n";
	system("bash", "-c", $command);

	$filenames[3] = $output_prefix . ".tlx.neg.bed";
	$command = "cat ".$filenames[1].' | perl -ne \'@F=split/\t/; if($F[5]=~/-/){ print; }\' >'.$filenames[3];
	print STDERR "[CMD] $command\n";
	system("bash", "-c", $command);
	
	$filenames[7] = $output_prefix . ".intersect.pos.bed";
	$command = "bedtools intersect -$opt  -a ".$filenames[2]."  -b $bedfile  >".$filenames[7];
	if($opt eq "c"){
		$command = "bedtools intersect -$opt  -b ".$filenames[2]."  -a $bedfile  >".$filenames[7];
	}
	print STDERR "[CMD] $command\n";
	system("bash", "-c", $command);
	
	$filenames[8] = $output_prefix . ".intersect.neg.bed";
	$command = "bedtools intersect -$opt  -a ".$filenames[3]."  -b $bedfile  >".$filenames[8];
	if($opt eq "c"){
		$command = "bedtools intersect -$opt  -b ".$filenames[3]."  -a $bedfile  >".$filenames[8];
	}
	print STDERR "[CMD] $command\n";
	system("bash", "-c", $command);
}


0;

