#!/usr/bin/perl

### ref: yyx_overlapping_bed_to_nonoverlapping_bdg.20210616.pl

use strict;
use warnings;

my $version = '
Version: 0.1.0 (2021-06-16)
Author: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: $0 <chrom.size> <output_prefix> <overlapping_score_operator>
Input:  STDIN
Output:  (skip if already exists)
	<output_prefix>.sorted.bdg
	<output_prefix>.nonoverlapping.sorted.bdg
	<output_prefix>.bw
Options:  <overlapping_score_operator> can be one of the following options:
	abs_max, abs_min, max, min, mean, bp_switch, comma_separated
".$version;


print STDERR "call sort\t".`which sort`;
print STDERR "call bedops\t".`which bedops`;
print STDERR "call bedmap\t".`which bedmap`;
print STDERR "call bedGraphToBigWig\t".`which bedGraphToBigWig`;
print STDERR "\n";


if(@ARGV < 3){
	die $usage;
}
my ($chrom_size_filename, $output_prefix, $score_operator) = @ARGV;

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

## sort input bdg
my $sorted_input_bed = $output_prefix.".sorted.bed";
my ($NR, @F);
if(!exist_file_or_dir($sorted_input_bed, "Skip sorting input bdg.")){
	print STDERR "Now sorting input bdg, and output to $sorted_input_bed ...\n";
	open(OUT, " | sort -k1,1 -k2,2n >$sorted_input_bed") or die "Error: cannot sort and output to $sorted_input_bed\n";
	$NR = 0;
	while(<STDIN>){
		$NR++;
		s/[\r\n]+$//;
		@F=split/\t/;
		print OUT join("\t", @F[0..2], $NR, $F[3])."\n";;
	}
	close(OUT);
}

## split overlapping if required
my $nonoverlapping_bdg = $output_prefix.".nonoverlapping.sorted.bdg";
my ($i, @G);
if(!exist_file_or_dir($nonoverlapping_bdg, "Skip splitting overlapping bed.")){
	print STDERR "Now splitting overlapping bed, and output to $nonoverlapping_bdg ...\n";
	open(OUT, ">".$nonoverlapping_bdg) or die "Error: cannot output to $nonoverlapping_bdg\n";
	open(IN, 'bedops --partition "'.$sorted_input_bed.'" | bedmap --delim "\t" --echo --echo-map-score - "'.$sorted_input_bed.'" | ') or die "Error: cannot split $sorted_input_bed by bedops and bedmap\n";
	while(<IN>){
#		print STDERR "[DEBUG] ".$_;
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
					print OUT join("\t", $F[0], $i, $i+1, $G[($i-$F[1]) % @G])."\n";
				}
			}
		}
		if($score_operator ne "bp_switch"){
			print OUT join("\t", @F)."\n";
		}
	}
	close(IN);
	close(OUT);
}

## bedGraphToBigWig
my $output_bw_filename = $output_prefix.".bw";
my $command = "bedGraphToBigWig $nonoverlapping_bdg $chrom_size_filename $output_bw_filename";
check_file_then_exec_command($output_bw_filename, $command, 1,1,0);


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

sub min{
	my $ans=undef;
	foreach (@_){
		if(!defined($ans) || $_ < $ans){
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



sub check_final_file_then_remove_intermediate_file{
	my ($final_filename, $intermediate_filename) = @_;
	my $command;
	if(ref($intermediate_filename) eq ""){
		$intermediate_filename = [$intermediate_filename];
	}
	if(exist_file_or_dir($final_filename, "So remove intermediate files...")){
		$command = "rm ".join(" ", @{$intermediate_filename});
		print STDERR "[PERL-SYSTEM] ".$command."\n";
		system($command);
	}
}


sub check_file_then_exec_command{
	my ($filename, $command, $should_time, $error_stop, $not_run) = @_;
	my $start_time = time();

	print STDERR "[PERL-SYSTEM] ".$command."\n";
	if(exist_file_or_dir($filename, "Skip this above command...")){
		return;
	}

	if(!(defined($not_run) && $not_run!=0)){
		if(system("/bin/bash", "-c", $command)!=0){
			if(defined($error_stop) && $error_stop!=0){
				die "Error: when exec last system command, return value = $?\n";
			}
		}
	}

	if(defined($should_time) && $should_time){
		check_elapsed_time($start_time);
	}
}

sub check_elapsed_time{
	my ($start_time, $end_time, $elapsed_time);
	my ($hour, $min, $sec, $day);
	$start_time = shift(@_);
	$end_time = time();
	$elapsed_time = $end_time - $start_time;
	$day = int($elapsed_time / (3600*24));
	$hour = int($elapsed_time % (3600*24) / 3600);
	$min = int($elapsed_time % 3600 / 60);
	$sec = $elapsed_time % 60;
	$elapsed_time = "";
	if($day>0){ $elapsed_time .= $day."day "; }
	if($hour>0){ $elapsed_time .= $hour."h"; }
	if($min>0){ $elapsed_time .= $min."min"; }
	if($sec>0 || $elapsed_time eq ""){ $elapsed_time .= $sec."s"; }
	print STDERR "[PERL-SYSTEM-TIME] ".$elapsed_time."\n";
}

sub exist_file_or_dir{
	my ($filenames, $str) = @_;
	my $returnValue = 0;
	if(ref($filenames) eq ""){
		$filenames = [$filenames];
	}
	if(!defined($str)){
		$str = "";
	}

	foreach my $filename (@{$filenames}){
		if(defined($filename) && -e $filename){
			if(-d $filename){
				if(! check_is_empty_dir($filename)){
					print STDERR "[CHECK-EXIST] Dir ".$filename." has already existed, and not empty. ".$str."\n";
					$returnValue = 1;
					last;
				}
			}elsif(-s $filename >= 100){
				print STDERR "[CHECK-EXIST] File ".$filename." has already existed. ".$str."\n";
				$returnValue = 1;
				last;
			}
		}
	}

	return $returnValue;
}


sub check_is_empty_dir{
	my ($dirname) = @_;
	my $dirHandle;
	my @contents;
	if(-d $dirname){
		if(opendir($dirHandle, $dirname)){
			@contents = readdir($dirHandle);
			closedir($dirHandle);
			if(scalar(@contents)>2){		# empty dir has . and ..
				return 0;		# not empty dir
			}else{
				return 1;		# is empty dir
			}
		}else{
			print STDERR ("Warning: cannot open dir $dirname\n");
			return -1;		# Cannot open dir
		}
	}else{
		print STDERR ("Warning: Not a dir is attemped to be checked\n");
		return -2;		# Not a dir
	}
}

