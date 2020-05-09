#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.2 (2020-03-19)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: cat <HTGTS_VDJ_annotated.tsv> | $0 <output_prefix> <mm9|mm9AJ>
	[allowed_possible_D_num (default: 0)] [ref_D_usage.tsv] [bait_IGHJ]
Input:
	[ref_D_usage.tsv]	two columns: D, usage(%)
Output:
	<output_prefix>.D_reannotated.tsv
	<output_prefix>.D_usage.tsv
".$version;

if(@ARGV < 2){
	die $usage;
}
my $output_prefix = shift(@ARGV);
my $mm9_mm9AJ = shift(@ARGV);

my $allowed_possible_D_num = 0;
if(@ARGV > 0){
	$allowed_possible_D_num = shift(@ARGV);
}

my $D_usage_pseudocount = 0.0001;
my $ref_D_usage_filename = undef;
my %D_usage = ();
my ($now_D, $now_usage);
my $total_usage = 0;
my (@F);
if(@ARGV > 0){
	$ref_D_usage_filename = shift(@ARGV);
	open(IN, $ref_D_usage_filename) or die "Error: cannot open ref_D_usage $ref_D_usage_filename for input\n";
	while(<IN>){
		if(/^\s*$/){  next;  }
		s/[\r\n]+$//;
		@F = split/\t/;
		$now_usage = $F[1] + $D_usage_pseudocount;
		$D_usage{$F[0]} = $now_usage;
		$total_usage += $now_usage;
	}
	close(IN);
	
	foreach $now_D (keys %D_usage){
		$D_usage{$now_D} /= $total_usage;
	}
}

my $bait_IGHJ = "IGHJ";
if(@ARGV > 0){
	$bait_IGHJ = shift(@ARGV);
}



my @output_D_order = qw/
IGHD1-3
IGHD3-1
IGHD1-1
IGHD2-3
IGHD2-2
IGHD3-3
IGHD1-2
IGHD2-1
IGHD2-4
IGHD2-14
IGHD2-9
IGHD2-10
IGHD2-11
IGHD3-2
IGHD4-1
/;
my @output_D_alias = qw/
DFL16.3
DST4.2
DFL16.1
DSP2.9
DSP2.3
DST4.3
DFL16.2
DSP2.5
DSP2.2a
DSP2.11
DSP2.2b
DSP2.8
DSP2.7
DST4
DQ52
/;
my %D2alias = ();
my ($i);
for($i=0; $i<@output_D_order; $i++){
	$D2alias{$output_D_order[$i]} = $output_D_alias[$i];
}

if($mm9_mm9AJ eq "mm9"){
	@output_D_order = qw/
Ighd5-1
Ighd3-1
Ighd1-1
Ighd6-1
Ighd2-3
Ighd6-2
Ighd2-4
Ighd5-2
Ighd2-5
Ighd5-3
Ighd5-7
Ighd2-6
Ighd5-4
Ighd5-8
Ighd2-7
Ighd5-5
Ighd2-8
Ighd5-6
Ighd3-2
Ighd4-1
/;
	%D2alias = (
"Ighd1-1" => "DFL16.1",
"Ighd3-2" => "DST4",
"Ighd4-1" => "DQ52"
);

}

push(@output_D_order, "-");
my %sum_D_reads = ();
for($i=0; $i<@output_D_order; $i++){
	$sum_D_reads{$output_D_order[$i]} = 0;
}



my $output_filename = $output_prefix . ".D_reannotated.tsv";
open(OUT, ">" . $output_filename) or die "Error: cannot open $output_filename for output\n";
my $headline = <STDIN>;
$headline =~ s/[\r\n]+$//;
my @fields = split(/\t/, $headline);
print OUT join("\t", @fields, "mid_D_reannotate")."\n";
my %field2idx = ();
for($i=0; $i<@fields; $i++){
	$field2idx{$fields[$i]} = $i;
}
my ($mid_D_annotate_str, @possible_mid_Ds, @possible_D_prob_vec, $now_sum, @str_vec);
my $total_reads = 0;
while(<STDIN>){
	s/[\r\n]+$//;
	@F = split/\t/;

	## For debug,
#	if($F[$field2idx{"Qname"}] eq "M01407:520:000000000-CCJHF:1:2105:2255:15067"){
#		print STDERR join("\t", @F)."\n";
#		print STDERR "'".$F[$field2idx{"bait_overlap_features"}]."'\n";
#		print STDERR "'".$F[$field2idx{"junction_overlap_features"}]."'\n";
#		print STDERR "'".$F[$field2idx{"mid_D_annotate"}]."'\n";
#	}
	if($F[$field2idx{"bait_overlap_features"}] =~ /^$bait_IGHJ/i){   # bait = IGHJ
		if($F[$field2idx{"junction_overlap_features"}] =~ /^IGHV/i){   # junction = IGHV
			$mid_D_annotate_str = $F[$field2idx{"mid_D_annotate"}];
			@possible_mid_Ds = unique(map{  s/\s*\(rC\)//g; $_;  } split(/, /, $mid_D_annotate_str));
#			print STDERR "[DEBUG] ".join("\t", $mid_D_annotate_str)."\n";
#			print STDERR "[DEBUG] ".join("\t", @possible_mid_Ds)."\n";
#			if($F[$field2idx{"Qname"}] eq "M01407:520:000000000-CCJHF:1:2105:2255:15067"){
#				print STDERR "'".$mid_D_annotate_str."'\n";
#				print STDERR "'".scalar(@possible_mid_Ds)."'\n";
#			}
			$total_reads++;
			if($mid_D_annotate_str eq "-" || @possible_mid_Ds > $allowed_possible_D_num){
				$mid_D_annotate_str = "-";
				$sum_D_reads{"-"} ++;
			}else{
				@possible_D_prob_vec = (1/scalar(@possible_mid_Ds)) x @possible_mid_Ds;   # evenly assigned, if no ref_D_usage is proviced
				if(defined($ref_D_usage_filename)){
					@possible_D_prob_vec = @D_usage{@possible_mid_Ds};
					$now_sum = sum(@possible_D_prob_vec);
					@possible_D_prob_vec = map {  $_ / $now_sum;  } @possible_D_prob_vec;
				}
				@str_vec = ();
				foreach $i (sort {
							$possible_D_prob_vec[$b] <=> $possible_D_prob_vec[$a];
						} 0..(@possible_mid_Ds-1)){
					$sum_D_reads{$possible_mid_Ds[$i]} += $possible_D_prob_vec[$i];
					if(exists($D2alias{$possible_mid_Ds[$i]})){
						push(@str_vec, $possible_mid_Ds[$i] . "(" . $D2alias{$possible_mid_Ds[$i]} . "):" . $possible_D_prob_vec[$i]);
					}else{
						push(@str_vec, $possible_mid_Ds[$i] . ":" . $possible_D_prob_vec[$i]);
					}
				}
				$mid_D_annotate_str = join(", ", @str_vec);
#				print STDERR $total_reads."\n";
			}
#			if($F[$field2idx{"Qname"}] eq "M01407:520:000000000-CCJHF:1:2105:2255:15067"){
#				print STDERR "'".$mid_D_annotate_str."'\n";
#			}
			print OUT join("\t", @F, $mid_D_annotate_str)."\n";
		}
	}
}
close(OUT);

#print STDERR $total_reads."\n";
if($total_reads < 0.0001){
	$total_reads = 0.0001
}
$output_filename = $output_prefix . ".D_usage.tsv";
open(OUT, ">" . $output_filename) or die "Error: cannot open $output_filename for output\n";
foreach $now_D (@output_D_order){
#	print STDERR join("'\t'", $now_D, $sum_D_reads{$now_D}, $total_reads)."\n";
	print OUT join("\t", $now_D, $sum_D_reads{$now_D}, $sum_D_reads{$now_D} / $total_reads)."\n";
}
close(OUT);

my @more_Ds = setdiff([keys %sum_D_reads], \@output_D_order);
if(@more_Ds > 0){
	print STDERR "Warning: some D existed in input file but are not output:\n";
	foreach (@more_Ds){
		print STDERR "\t".$_."\n";
	}
}


0;



sub unique{
	my %S = ();
	my @ans = ();
	foreach (@_){
		if(!exists($S{$_})){
			$S{$_} = 1;
			push(@ans, $_);
		}
	}
	return(@ans);
}

sub sum{
	my $ans = 0;
	foreach (@_){
		$ans += $_;
	}
	return($ans);
}

sub setdiff{
	my ($set1, $set2) = @_;
	my %S = ();
	my @ans = ();
	foreach (@$set1){
		$S{$_} = 1;
	}
	foreach (@$set2){
		$S{$_} = 0;
	}
	foreach (sort keys %S){
		if($S{$_} == 1){
			push(@ans, $_);
		}
	}
	return @ans;
}




