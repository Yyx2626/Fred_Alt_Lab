#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2020-01-20)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.bdg|bw> <output.bdg> <multiply_factor> [chr] [start] [end]
".$version;

if(@ARGV < 3){
	die $usage;
}
my $input_filename = shift(@ARGV);
my $output_filename = shift(@ARGV);
my $multiply_factor = shift(@ARGV);

print STDERR "[DEBUG] input_filename=$input_filename\n";
print STDERR "[DEBUG] output_filename=$output_filename\n";
print STDERR "[DEBUG] multiply_factor=$multiply_factor\n";

my ($chr, $start, $end) = (undef, undef, undef);
if(@ARGV > 0){
	$chr = shift(@ARGV);
	$start = shift(@ARGV);
	$end = shift(@ARGV);
	print STDERR "[DEBUG] chr=$chr\n";
	print STDERR "[DEBUG] start=$start\n";
	print STDERR "[DEBUG] end=$end\n";
}

my (@F);

if($input_filename =~ /\.bw$/){
	if(defined($end)){
		open(IN, "bigWigToBedGraph -chrom=$chr -start=$start -end=$end  $input_filename  /dev/stdout | ") or die "Error: cannot input from bigWigToBedGraph $input_filename with $chr:$start-$end restriction\n";
	}else{
		open(IN, "bigWigToBedGraph  $input_filename  /dev/stdout | ") or die "Error: cannot input from bigWigToBedGraph $input_filename\n";
	}
}else{
	open(IN, $input_filename) or die "Error: cannot open $input_filename for input\n";
}
open(OUT, ">".$output_filename) or die "Error: cannot open $output_filename for output\n";

while(<IN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	if(defined($end)){
		if($F[0] eq $chr){
			if($F[2] >= $start && $F[1] <= $end){
				$F[3] *= $multiply_factor;
				print OUT join("\t", @F)."\n";
			}
		}
	}else{
		$F[3] *= $multiply_factor;
		print OUT join("\t", @F)."\n";
	}
}

close(OUT);
close(IN);


0;
