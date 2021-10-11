#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Author: Adam Yongxin Ye @ BCH
Version: 0.1.1 (2020-09-04)
';
my $usage = "Usage: $0 <input.tlx> <chromSize> <output_prefix>
	[should_strand (default:1)] [normalize_to (default:0)]
	[yyx_bdg_extract_multiply.pl (required for normalization)]
	[junction|prey (default:junction)]
".$version;

if(@ARGV < 3){
	die $usage;
}
my ($input_filename, $chromSize_filename, $output_prefix) = @ARGV;

my $should_strand = 1;
if(@ARGV > 3){
	$should_strand = $ARGV[3];
}
my $normalize_to = 0;
if(@ARGV > 4){
	$normalize_to = $ARGV[4];
}
my $yyx_bdg_extract_multiply_pl = "";
if(@ARGV > 5){
	$yyx_bdg_extract_multiply_pl = $ARGV[5];
}
if($normalize_to > 0 && $yyx_bdg_extract_multiply_pl eq ""){
	die "Error: need to provide the path of yyx_bdg_extract_multiply.pl for normalization\n";
}

my $retrieve_element = "junction";
if(@ARGV > 6){
	$retrieve_element = $ARGV[6];
}


my $start_time = time();
print STDERR "[PERL-START] ".scalar(localtime())."\n";

my @filenames = ();


### step 1. convert tlx to (junction) bed
$filenames[0] = $output_prefix . ".bed";
$filenames[1] = $output_prefix . ".pos.bed";
$filenames[2] = $output_prefix . ".neg.bed";

my $should_skip = 0;
if($should_strand){
	if(exist_file_or_dir($filenames[1], "") && exist_file_or_dir($filenames[2], "")){
		print STDERR "   ... Skip generating bed ...\n";
		$should_skip = 1;
	}else{
		open(OUT1, ">".$filenames[1]) or die "Error: cannot open ".$filenames[1]." for output\n";
		open(OUT2, ">".$filenames[2]) or die "Error: cannot open ".$filenames[1]." for output\n";
	}
}else{
	if(exist_file_or_dir($filenames[0], "Skip generating bed ...")){
		$should_skip = 1;
	}else{
		open(OUT, ">".$filenames[0]) or die "Error: cannot open ".$filenames[1]." for output\n";
	}
}
my $rowidx = undef;
if(!$should_skip){
	open(IN, $input_filename) or die "Error: cannot open tlx file $input_filename for input\n";
	my $headline = <IN>;
	$headline =~ s/[\r\n]+$//;
	my @F = split(/\t/, $headline);
	my %f2i = ();
	my $i;
	for($i=0; $i<@F; $i++){
		$f2i{$F[$i]} = $i;
	}
	my @required_fields = qw/Qname Rname Strand/;
	if($retrieve_element eq "junction"){
		push(@required_fields, "Junction");
	}elsif($retrieve_element eq "prey"){
		push(@required_fields, qw/Rstart Rend/);
	}else{
		die "Error: unknown retrieve $retrieve_element\n";
	}
	foreach (@required_fields){
		if(!exists($f2i{$_})){
			die "Error: cannot find column $_ in input tlx file\n";
		}
	}
	$rowidx = 1;
	my ($chr, $start,$end, $name, $strand);
	while(<IN>){
		s/[\r\n]+$//;
		if(/^\s*$/){ next; }
		@F = split/\t/;
		$chr = $F[$f2i{"Rname"}];
		if(!defined($chr)){ next; }
		$end = 0;
		$start = 0;
		if($retrieve_element eq "junction"){
			$end = $F[$f2i{"Junction"}];
			if(!defined($end)){ next; }
			$start = $end - 1;
			if($start < 0){ next; }
		}else{
			$end = $F[$f2i{"Rend"}];
			if(!defined($end)){ next; }
			$start = $F[$f2i{"Rstart"}] - 1;
			if($start < 0){ next; }
		}
		$name = $F[$f2i{"Qname"}];
		$strand = "+";
		if($F[$f2i{"Strand"}] < 0){
			$strand = "-";
		}
		if($should_strand){
			if($strand eq "+"){
				print OUT1 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
			}elsif($strand eq "-"){
				print OUT2 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
			}
		}else{
			print OUT join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
		}
		$rowidx++;
	}
	close(IN);
	if($should_strand){
		close(OUT1);
		close(OUT2);
	}else{
		close(OUT);
	}
}


### step 2. convert bed to bdg, by bedtools genomecov
my @file_last_idxes = (0);
if($should_strand){
	@file_last_idxes = (1, 2);
}
$filenames[10] = $output_prefix . ".sorted.bed";
$filenames[11] = $output_prefix . ".pos.sorted.bed";
$filenames[12] = $output_prefix . ".neg.sorted.bed";
$filenames[20] = $output_prefix . ".bdg";
$filenames[21] = $output_prefix . ".pos.bdg";
$filenames[22] = $output_prefix . ".neg.bdg";
my $command;
my $fIdx;
foreach $fIdx (@file_last_idxes){
	$command = "bedtools sort -i ".$filenames[$fIdx]." >".$filenames[10+$fIdx];
	check_file_then_exec_command([$filenames[10+$fIdx], $filenames[20+$fIdx]], $command, 1, 0, 0);
	$command = "bedtools genomecov -bg -i ".$filenames[10+$fIdx]." -g ".$chromSize_filename." >".$filenames[20+$fIdx];
	check_file_then_exec_command([$filenames[20+$fIdx]], $command, 1, 0, 0);
}


### step 3. normalize, if needed
my $multiply_factor;
if($normalize_to > 0){
	$filenames[30] = $output_prefix . ".norm_".$normalize_to.".bdg";
	$filenames[31] = $output_prefix . ".norm_".$normalize_to.".pos.bdg";
	$filenames[32] = $output_prefix . ".minus_norm_".$normalize_to.".neg.bdg";
	$should_skip = 0;
	if($should_strand){
		if(exist_file_or_dir($filenames[31], "") && exist_file_or_dir($filenames[32], "")){
			print STDERR "   ... Skip normalization ...\n";
			$should_skip = 1;
		}
	}else{
		if(exist_file_or_dir($filenames[30], "Skip normalization ...")){
			$should_skip = 1;
		}
	}
	if(!$should_skip){
		if(!defined($rowidx)){
			$rowidx = `grep -v "^\\s*\$" ".$input_filename." | wc -l` - 1;
		}
		$multiply_factor = $normalize_to / $rowidx;
		foreach $fIdx (@file_last_idxes){
			if($fIdx==2){
				$multiply_factor = -$multiply_factor;
			}
			$command = "perl ".$yyx_bdg_extract_multiply_pl." ".$filenames[20+$fIdx]." ".$filenames[30+$fIdx]." ".$multiply_factor;
			check_file_then_exec_command([$filenames[30+$fIdx]], $command, 1, 0, 0);
		}
	}
}else{
	foreach $fIdx (@file_last_idxes){
		$filenames[30+$fIdx] = $filenames[20+$fIdx];
		if($fIdx==2){
			$filenames[32] = $output_prefix . ".minus.neg.bdg";
			$command = "perl ".$yyx_bdg_extract_multiply_pl." ".$filenames[20+$fIdx]." ".$filenames[30+$fIdx]." -1";
			check_file_then_exec_command([$filenames[30+$fIdx]], $command, 1, 0, 0);
		}
	}
}


### step 4. convert bdg to bw, by bedGraphToBigWig
$filenames[40] = $output_prefix . ".bw";
$filenames[41] = $output_prefix . ".pos.bw";
$filenames[42] = $output_prefix . ".neg.bw";
foreach $fIdx (@file_last_idxes){
	$command = "bedGraphToBigWig ".$filenames[30+$fIdx]." ".$chromSize_filename." ".$filenames[40+$fIdx];
	check_file_then_exec_command([$filenames[40+$fIdx]], $command, 1, 0, 0);
}




check_elapsed_time($start_time);
print STDERR "[PERL-END] ".scalar(localtime())."\n";

0;




sub generateRandomString{
	my ($length) = @_;
	my $ans = "";
	my ($random, $random2);
	my ($i, $tmp);
	for($i=0; $i<$length; $i++){
		$random = rand();
		$random2 = $random * 2 - int($random * 2);
		if($random * 2 < 1){
			$tmp = int(ord('A') + $random2 * 26);
			if($tmp > ord('Z')){  $tmp = ord('Z'); }
			$tmp = chr($tmp);
			$ans .= $tmp;
		}else{
			$tmp = int(ord('a') + $random2 * 26);
			if($tmp > ord('z')){  $tmp = ord('z'); }
			$tmp = chr($tmp);
			$ans .= $tmp;
		}
	}
	return($ans);
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




