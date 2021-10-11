#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

my $version = "
	Version: 0.1.0 (2012-12-02)
	Author: Adam_YeYongxin, CBI
";
my $usage = "
	Usage: $0 <files> ...
	Options:
		-d STR	set delimiter (default: \\t)
		-n INT	the number of head lines that will be checked
			(default: 4; -1 for all lines)
".$version;

my $delimiter = "\t";
my $num_check_lines = 4;
GetOptions(
	"d=s" => \$delimiter,
	"n=i" => \$num_check_lines,
);

#print STDERR "[DEBUG] delimiter='$delimiter', num_check_lines=$num_check_lines\n";
#die;

my @filenames;
if(scalar(@ARGV)<1){
	die $usage;
}else{
	@filenames = @ARGV;
}

my $filename;
my @F;
my $i;
foreach $filename (@filenames){
	if(open(INPUT, $filename)){
		print "# Now process file '$filename' :\n";
		$i = 0;
		while(<INPUT>){
			$i++;
			if($num_check_lines>=0 && $i>$num_check_lines){
				last;
			}
			@F = split/$delimiter/;
			print scalar(@F)."\n";
		}
		close(INPUT);
	}else{
		print STDERR "Warn: cannot open file '$filename' for input, so I just skip it\n";
	}
}

0;

