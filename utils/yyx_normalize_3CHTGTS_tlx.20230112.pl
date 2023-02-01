#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Author: Adam Yongxin Ye @ BCH
Version: 0.1.1 (2023-01-12)
';
my $usage = "Usage: $0 <input.tlx> <chromSize> <output_prefix>
	[signal_coordinate (default: all; example (Igk + nearby): chr6:64515000-73877000; or chr6)]
	[rm_bait_peak_coordinate (default: none; example: chr6:70675300-70675400)]
	[normalized_to (default: 1000000)] [remove_intermediate_files (default: 1)]

Required tools:
	bedtools
	bedGraphToBigWig

Intermediate output:
	<output_prefix>.  both/pos/neg  .bed , .bdg , .bw (if bedGraphToBigWig exists)
	<output_prefix>.rm_bait_peak.  both/pos/neg  .bed , .bdg , .bw (if rm_bait_peak_coordinate != none)

Final output:
	<output_prefix>.junction_count.txt
	<output_prefix>.norm_from_*_to_*.  both/pos/neg  .bdg or .bw (if bedGraphToBigWig exists)
	<output_prefix>.rm_bait_peak.norm_from_*_to_*.  both/pos/neg  .bdg or .bw (if rm_bait_peak_coordinate != none)
".$version;

if(@ARGV < 3){
	die $usage;
}
my $input_filename = shift(@ARGV);
my $chromSize_filename = shift(@ARGV);
my $output_prefix = shift(@ARGV);
print STDERR "[DEBUG] input_filename = '$input_filename'\n";
print STDERR "[DEBUG] chromSize_filename = '$chromSize_filename'\n";
print STDERR "[DEBUG] output_prefix = '$output_prefix'\n";

my $signal_coordinate = "all";
if(@ARGV > 0){
	$signal_coordinate = shift(@ARGV);
}
print STDERR "[DEBUG] signal_coordinate = '$signal_coordinate'\n";

my $rm_bait_peak_coordinate = "none";
if(@ARGV > 0){
	$rm_bait_peak_coordinate = shift(@ARGV);
}
print STDERR "[DEBUG] rm_bait_peak_coordinate = '$rm_bait_peak_coordinate'\n";

my $normalized_to = 1000000;
if(@ARGV > 0){
	$normalized_to = shift(@ARGV);
}
print STDERR "[DEBUG] normalized_to = '$normalized_to'\n";

my $should_remove_intermediate_files = 1;
if(@ARGV > 0){
	$should_remove_intermediate_files = shift(@ARGV);
}
print STDERR "[DEBUG] should_remove_intermediate_files = '$should_remove_intermediate_files'\n";


my @Signal_coord = parse_coordinate_string($signal_coordinate);
my @RmBait_coord = parse_coordinate_string($rm_bait_peak_coordinate);


my $start_time = time();
print STDERR "[PERL-START] ".scalar(localtime())."\n";


my $step_start_time = time();
print STDERR "\n[STEP1] convert tlx to (junction) bed, and count junctions\n";

### step 1. convert tlx to (junction) bed, and count junctions
my ($total_count, $signal_count, $rmBait_count) = parse_tlx_to_count_and_bed($input_filename, $output_prefix, $signal_coordinate, $rm_bait_peak_coordinate);
my $signalRmBait_count = $signal_count - $rmBait_count;

my $output_filename = $output_prefix . ".junction_count.txt";
if(open(OUT, ">".$output_filename)){
	print OUT join("\t", "total", $total_count)."\n";
	print OUT join("\t", "signal", $signal_count)."\n";
	print OUT join("\t", "bait", $rmBait_count)."\n";
	print OUT join("\t", "signal_rm_bait", $signalRmBait_count)."\n";
	close(OUT);
}else{
	print STDERR "Warning: cannot open $output_filename for output, so I skip output\n";
}

check_elapsed_time($step_start_time);
print STDERR "[STEP1-END] ".scalar(localtime())."\n";



$step_start_time = time();
print STDERR "\n[STEP2] convert bed to bdg and bw\n";

my $norm_str = "norm_from_" . $signalRmBait_count . "_to_" . $normalized_to;
my $multiply_factor = $normalized_to / $signalRmBait_count;

my @raw_file_prefixes = ();
my @norm_file_prefixes = ();
my ($strand, $ext);
foreach $strand (qw/both pos neg/){
	push(@raw_file_prefixes, $output_prefix . "." . $strand);
	push(@norm_file_prefixes, $output_prefix . "."  . $norm_str . "." . $strand);
}
if($RmBait_coord[0] ne "none"){
	foreach $strand (qw/both pos neg/){
		push(@raw_file_prefixes, $output_prefix . ".rm_bait_peak." . $strand);
		push(@norm_file_prefixes, $output_prefix . ".rm_bait_peak."  . $norm_str . "." . $strand);
	}
}

my $file_prefix;
my @intermediate_files = ();
my @final_files = ();
foreach $file_prefix (@raw_file_prefixes){
	foreach $ext (qw/bed bdg bw/){
		push(@intermediate_files, $file_prefix . "." . $ext);
	}
}
foreach $file_prefix (@norm_file_prefixes){
	foreach $ext (qw/bdg bw/){
		push(@final_files, $file_prefix . "." . $ext);
	}
}


### step 2. convert bed to bdg and bw
foreach $file_prefix (@raw_file_prefixes){
	convert_bed_to_bdg($file_prefix, $chromSize_filename);
	convert_bdg_to_bw($file_prefix, $chromSize_filename);
}

check_elapsed_time($step_start_time);
print STDERR "[STEP2-END] ".scalar(localtime())."\n";



### step 3. normalize, if needed
my $i;
if($normalized_to > 0){
	$step_start_time = time();
	print STDERR "\n[STEP3] normalize\n";

	for($i=0; $i<@raw_file_prefixes; $i++){
		bdg_extract_multiply($raw_file_prefixes[$i] . ".bdg", $norm_file_prefixes[$i] . ".bdg", "all", $multiply_factor);
		convert_bdg_to_bw($norm_file_prefixes[$i], $chromSize_filename);
	}
	
	check_elapsed_time($step_start_time);
	print STDERR "[STEP3-END] ".scalar(localtime())."\n";
}else{
	print STDERR "normalized_to = $normalized_to <= 0, so I skip normalization\n";
}


### step 4. remove intermediate files
if($should_remove_intermediate_files){
	check_final_file_then_remove_intermediate_file([@final_files], [@intermediate_files]);
	foreach $file_prefix (@norm_file_prefixes){
		check_final_file_then_remove_intermediate_file($file_prefix . ".bw", $file_prefix . ".bdg");
	}
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



sub parse_coordinate_string{
	my ($input_coordinate_str) = @_;
	if($input_coordinate_str =~ /all/i){
		return ("all", -1, -1);
	}
	if($input_coordinate_str =~ /none/i){
		return ("none", -1, -1);
	}
	my @F = split(/[:-]/, $input_coordinate_str);
	if(@F == 0){
		return ("none", -1, -1);
	}
	while(@F < 3){
		push(@F, -1);
	}
	return @F;
}

sub parse_tlx_to_count_and_bed{
	my ($input_filename, $output_prefix, $signal_coordinate, $rm_bait_peak_coordinate) = @_;

	my @Signal_coord = parse_coordinate_string($signal_coordinate);
	my @RmBait_coord = parse_coordinate_string($rm_bait_peak_coordinate);
	
	my @filenames = ();

	### step 1. convert tlx to (junction) bed
	$filenames[0] = $output_prefix . ".both.bed";
	$filenames[1] = $output_prefix . ".pos.bed";
	$filenames[2] = $output_prefix . ".neg.bed";

	$filenames[10] = $output_prefix . ".rm_bait_peak.both.bed";
	$filenames[11] = $output_prefix . ".rm_bait_peak.pos.bed";
	$filenames[12] = $output_prefix . ".rm_bait_peak.neg.bed";

	open(OUT, ">".$filenames[0]) or die "Error: cannot open ".$filenames[0]." for output\n";
	open(OUT1, ">".$filenames[1]) or die "Error: cannot open ".$filenames[1]." for output\n";
	open(OUT2, ">".$filenames[2]) or die "Error: cannot open ".$filenames[2]." for output\n";
	if($RmBait_coord[0] ne "none"){
		open(OUT10, ">".$filenames[10]) or die "Error: cannot open ".$filenames[10]." for output\n";
		open(OUT11, ">".$filenames[11]) or die "Error: cannot open ".$filenames[11]." for output\n";
		open(OUT12, ">".$filenames[12]) or die "Error: cannot open ".$filenames[12]." for output\n";
	}
	
	my $rowidx = undef;
	my $retrieve_element = "junction";
	
	my $total_count = 0;
	my $signal_count = 0;
	my $rmBait_count = 0;
	
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
		if(/^\s*$/){
			print STDERR "Warning: empty line $rowidx, so I skip it\n";
			$rowidx++;
			next;
		}
		@F = split/\t/;
		$chr = $F[$f2i{"Rname"}];
		if(!defined($chr) || $chr =~ /^\s*$/){
			print STDERR "Warning: no Rname info for line $rowidx, so I skip it\n";
			$rowidx++;
			next;
		}
		$end = 0;
		$start = 0;
		if($retrieve_element eq "junction"){
			$end = $F[$f2i{"Junction"}];
			if(!defined($end)){
				print STDERR "Warning: no Junction info for line $rowidx, so I skip it\n";
				$rowidx++;
				next;
			}
			$start = $end - 1;
			if($start < 0){
				print STDERR "Warning: Junction < 1 for line $rowidx, so I skip it\n";
				$rowidx++;
				next;
			}
		}else{
			$end = $F[$f2i{"Rend"}];
			$start = $F[$f2i{"Rstart"}];
			if(!defined($end)){
				print STDERR "Warning: no Rend info for line $rowidx, so I skip it\n";
				$rowidx++;
				next;
			}
			if(!defined($start)){
				print STDERR "Warning: no Rstart info for line $rowidx, so I skip it\n";
				$rowidx++;
				next;
			}
			$start = $start - 1;
			if($start < 0){
				print STDERR "Warning: Rstart < 1 for line $rowidx, so I skip it\n";
				$rowidx++;
				next;
			}
		}
		$name = $F[$f2i{"Qname"}];
		$strand = "+";
		if($F[$f2i{"Strand"}] < 0){
			$strand = "-";
		}
		if($strand eq "+"){
			print OUT1 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
		}elsif($strand eq "-"){
			print OUT2 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
		}
		print OUT join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
		
		$total_count++;
		if($Signal_coord[0] eq "all"){
			$signal_count++;
		}else{
			if($chr eq $Signal_coord[0] && ($Signal_coord[1] < 0 || $end >= $Signal_coord[1]) && ($Signal_coord[2] < 0 || $end <= $Signal_coord[2])){
				$signal_count++;
			}
		}
		if($RmBait_coord[0] ne "none"){
			if($chr eq $RmBait_coord[0] && ($RmBait_coord[1] < 0 || $end >= $RmBait_coord[1]) && ($RmBait_coord[2] < 0 || $end <= $RmBait_coord[2])){
				$rmBait_count++;
			}else{
				if($strand eq "+"){
					print OUT11 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
				}elsif($strand eq "-"){
					print OUT12 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
				}
				print OUT10 join("\t", $chr, $start, $end, $name, $rowidx, $strand)."\n";
			}
		}
		$rowidx++;
	}
	close(IN);
	close(OUT);
	close(OUT1);
	close(OUT2);
	if($RmBait_coord[0] ne "none"){
		close(OUT10);
		close(OUT11);
		close(OUT12);
	}
	
	return ($total_count, $signal_count, $rmBait_count);
}

sub convert_bed_to_bdg{
	my ($file_prefix, $chromSize_filename) = @_;
	
	### step 2. convert bed to bdg, by bedtools genomecov
	my $input_filename = $file_prefix . ".bed";
	my $sorted_filename = $file_prefix . ".sorted.bed";
	my $outBdg_filename = $file_prefix . ".bdg";
	
	my $command;
	$command = "bedtools sort -i ".$input_filename." >".$sorted_filename;
	check_file_then_exec_command([$sorted_filename, $outBdg_filename], $command, 1, 0, 0);
	$command = "bedtools genomecov -bg -i ".$sorted_filename." -g ".$chromSize_filename." >".$outBdg_filename;
	check_file_then_exec_command([$outBdg_filename], $command, 1, 0, 0);
	check_final_file_then_remove_intermediate_file($outBdg_filename, $sorted_filename);
}

### step 3. normalize, if needed
sub bdg_extract_multiply{
	my ($input_filename, $output_filename, $extract_coordinate, $multiply_factor) = @_;
	
	my ($chr, $start, $end) = parse_coordinate_string($extract_coordinate);
	if($chr eq "all"){
		($chr, $start, $end) = (undef, undef, undef);
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
}

sub convert_bdg_to_bw{
	my ($file_prefix, $chromSize_filename) = @_;
	### step 4. convert bdg to bw, by bedGraphToBigWig
	my $in_bdg_filename = $file_prefix . ".bdg";
	my $out_bw_filename = $file_prefix . ".bw";
	my $command = "bedGraphToBigWig ".$in_bdg_filename." ".$chromSize_filename." ".$out_bw_filename;
	check_file_then_exec_command([$out_bw_filename], $command, 1, 0, 0);
}


