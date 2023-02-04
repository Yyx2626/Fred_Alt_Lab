#!/usr/bin/env perl

use strict;
use warnings;

use IO::File;
use Text::CSV;

my $tlxfile = shift;
my $bedfile = shift;
my $out = shift;

my %hash;

my $bedfh = IO::File->new("<$bedfile");
my $bedcsv = Text::CSV->new({sep_char => "\t"});
my @bedcolumns = qw(chr start end name);
$bedcsv->column_names(@bedcolumns);


while (my $bed = $bedcsv->getline_hr($bedfh)) {
  next unless defined $bed->{name} && $bed->{name} =~ /\S/;
  $hash{$bed->{name}} = 1;
}

my $tlxfh = IO::File->new("<$tlxfile");
my $tlxcsv = Text::CSV->new({sep_char => "\t"});
my $tlxcolumns = $tlxcsv->getline($tlxfh);
$tlxcsv->column_names(@$tlxcolumns);

my $outfh = IO::File->new(">$out");
$outfh->print(join("\t",@$tlxcolumns)."\n");

while (my $tlx = $tlxcsv->getline_hr($tlxfh)) {

  if (exists $hash{$tlx->{Qname}}) {
    $outfh->print(join("\t", @{$tlx}{@$tlxcolumns})."\n");
  }

}
