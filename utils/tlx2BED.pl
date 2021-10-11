#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Text::CSV;

my $tlxfile = shift;
my $bedfile = shift;


system("perl -pi -e 's/\r/\n/g' $tlxfile");

my $tlxfh = IO::File->new("<$tlxfile") or croak "Error: could not read file $tlxfile";
my $bedfh = IO::File->new(">$bedfile") or croak "Error: could not write file $bedfile";

my $csv = Text::CSV->new({sep_char=>"\t"});

my $header = $csv->getline($tlxfh);
$csv->column_names(map {lc} @$header);

while (my $read = $csv->getline_hr($tlxfh)) {
  next unless $read->{qname} =~ /\S/;
  my $strand = $read->{strand} == 1 ? "+" : "-";
  $bedfh->print(join("\t",$read->{rname},$read->{junction}-1,$read->{junction},$read->{qname},0,$strand)."\n");
}
