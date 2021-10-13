#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2021-06-16)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <input.sorted.bed> <overlapping_score_operator>
Input:  <input.sorted.bed>
Output:  STDOUT
Options:  <overlapping_score_operator> can be one of the following options:
	abs_max, abs_min, max, min, mean, bp_switch, comma_separated
".$version;


print STDERR "call bedops\t".`which bedops`;
print STDERR "call bedmap\t".`which bedmap`;
print STDERR "\n";


if(@ARGV < 2){
	die $usage;
}
my ($input_filename, $score_operator) = @ARGV;

if($score_operator eq "abs_max"){
}elsif($score_operator eq "abs_min"){
}elsif($score_operator eq "max"){
}elsif($score_operator eq "min"){
}elsif($score_operator eq "mean"){
}elsif($score_operator eq "bp_switch"){
}elsif($score_operator eq "comma_separated"){
}else{
	die "Error: unrecognized overlapping_score_operator = $score_operator
	should be one of the following options:
	abs_max, abs_min, max, min, mean, bp_switch, comma_separated\n";
}

my ($i, @F, @G);
open(IN, 'bedops --partition "'.$input_filename.'" | bedmap --delim "\t" --echo --echo-map-score - "'.$input_filename.'" | ') or die "Error: cannot split $input_filename by bedops and bedmap\n";
while(<IN>){
	s/[\r\n]+$//g;
	@F=split/\t/;
	@G=split(/;/, $F[3]);
	if(@G>0){
		if($score_operator eq "abs_max"){
			$F[3] = abs_max(@G);
		}elsif($score_operator eq "abs_min"){
			$F[3] = abs_min(@G);
		}elsif($score_operator eq "max"){
			$F[3] = max(@G);
		}elsif($score_operator eq "min"){
			$F[3] = min(@G);
		}elsif($score_operator eq "mean"){
			$F[3] = mean(@G);
		}elsif($score_operator eq "comma_separated"){
			$F[3] = join(",", @G);
		}elsif($score_operator eq "bp_switch"){
			for($i=$F[1]; $i<$F[2]; $i++){
				print join("\t", $F[0], $i, $i+1, $G[($i-$F[1]) % @G])."\n";
			}
		}
	}
	if($score_operator ne "bp_switch"){
		print join("\t", @F)."\n";
	}
}
close(IN);


0;



sub max{
	my $ans=undef;
	foreach (@_){
		if(!defined($ans) || $_ > $ans){
			$ans=$_;
		}
	}
	return $ans;
}

sub abs_max{
	my $ans=undef;
	foreach (@_){
		if(!defined($ans) || abs($_) > abs($ans)){
			$ans=$_;
		}
	}
	return $ans;
}

sub min{
	my $ans=undef;
	foreach (@_){
		if(!defined($ans) || $_ < $ans){
			$ans=$_;
		}
	}
	return $ans;
}

sub abs_min{
	my $ans=undef;
	foreach (@_){
		if(!defined(abs($ans)) || abs($_) < abs($ans)){
			$ans=$_;
		}
	}
	return $ans;
}

sub mean{
	my $S = 0;
	my $N = 0;
	foreach (@_){
		$S += $_;
		$N++;
	}
	return $S / $N;
}

