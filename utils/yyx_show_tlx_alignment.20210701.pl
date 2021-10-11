#!/usr/bin/perl

use strict;
use warnings;

my $version = '
Version: 0.1.1 (2021-07-01)
Authors: Adam Yongxin Ye @ BCH
';
my $usage = "Usage: cat <input.tlx> | perl $0 <ref.fa> <bait_slop>
	[junction_upstream] [junction_downstream]
	[anno.bed] [0-base(default:0)] [should_end_inclusive(default:0)]
".$version;

if(@ARGV < 2){
	die $usage;
}
my $ref_fasta_filename = $ARGV[0];
my $bait_slop = $ARGV[1];

my $junction_upstream = -1;
my $junction_downstream = -1;
if(@ARGV > 2){
	$junction_upstream = $ARGV[2];
}
if(@ARGV > 3){
	$junction_downstream = $ARGV[3];
}

my $anno_bed_filename = undef;
my %anno = ();
my @F;
if(@ARGV > 4){
	$anno_bed_filename = $ARGV[4];
}
my $bed_base = 0;
if(@ARGV > 5){
	$bed_base = $ARGV[5];
}
my $should_end_inclusive = 0;
if(@ARGV > 6){
	$should_end_inclusive = $ARGV[6];
}

if(defined($anno_bed_filename)){
	open(IN, $anno_bed_filename) or die "Error: cannot open anno_bed $anno_bed_filename\n";
	while(<IN>){
		s/[\r\n]+$//;
		@F = split/\t/;
		# convert to typical bed format
		$F[1] -= $bed_base;
		$F[2] -= $bed_base;
		$F[2] += $should_end_inclusive;
		$anno{$F[3]} = [@F];
	}
	close(IN);
}
my @anno_keys = sort keys %anno;



my %f2i = ();
$f2i{"Qname"} = 0;
$f2i{"Rname"} = 2;
$f2i{"Junction"} = 3;
$f2i{"Strand"} = 4;
$f2i{"Rstart"} = 5;
$f2i{"Rend"} = 6;
$f2i{"B_Rname"} = 7;
$f2i{"B_Rstart"} = 8;
$f2i{"B_Rend"} = 9;
$f2i{"B_Strand"} = 10;
$f2i{"B_Qstart"} = 11;
$f2i{"B_Qend"} = 12;
$f2i{"Qstart"} = 13;
$f2i{"Qend"} = 14;
$f2i{"Seq"} = 18;

my ($i);
my ($readSeq, $refBaitSeq, $refPreySeq, $refSeqMain, $refSeqSlop);
my ($chr, $start, $end, $strand, $bedStart, $bedEnd);
my ($readStart, $readEnd);
my (@baitAnno, @preyAnno);
while(<STDIN>){
	s/[\r\n]+$//;
	@F = split/\t/;
	if($F[0] eq "Qname"){
		%f2i = ();
		for($i=0; $i<@F; $i++){
			$f2i{$F[$i]} = $i;
		}
		next;
	}
	
	print join("\t", @F)."\n";
	
	$readSeq = $F[$f2i{"Seq"}];

	$chr = $F[$f2i{"B_Rname"}];
	$start = $F[$f2i{"B_Rstart"}] - 1;
	$end = $F[$f2i{"B_Rend"}];
	$strand = $F[$f2i{"B_Strand"}] > 0 ? "+" : "-";
	# slop towards bait end
	$bedStart = $start;
	$bedEnd = $end;
	if($strand eq "+"){
		$bedEnd += $bait_slop;
	}else{
		$bedStart -= $bait_slop;
	}
	## retrieve B_Rname:(B_Rstart-1)~(B_Rend+bait_slop) "+"
	##     or   B_Rname:(B_Rstart-1-bait_slop)~(B_Rend) "-"
	##  length = B_Rend - B_Rstart + 1 + bait_slop
	$refBaitSeq = `echo -e "$chr\t$bedStart\t$bedEnd\t\t\t$strand" | bedtools getfasta -tab -s -fi $ref_fasta_filename -bed -`;
	$refBaitSeq = (split(/\t/, $refBaitSeq))[1];
	$refBaitSeq =~ s/[\r\n]+$//;
#	if($strand eq "+"){
		$refSeqMain = substr($refBaitSeq, 0, $end-$start);
		$refSeqMain =~ tr/a-z/A-Z/;
		$refSeqSlop = substr($refBaitSeq, $end-$start);
		$refSeqSlop =~ tr/A-Z/a-z/;
		$refBaitSeq = $refSeqMain . $refSeqSlop;
#	}else{
#		$refSeqSlop = substr($refBaitSeq, 0, $bait_slop);
#		$refSeqSlop =~ tr/A-Z/a-z/;
#		$refSeqMain = substr($refBaitSeq, $bait_slop);
#		$refSeqMain =~ tr/a-z/A-Z/;
#		$refBaitSeq = $refSeqSlop . $refSeqMain;
#	}
	@baitAnno = retrieve_overlapping_anno_features($chr, $bedStart, $bedEnd);
	if($strand eq "-"){
		@baitAnno = map {
			$_->[1] = reverse_anno_str($_->[1]);
			$_;
		} @baitAnno;
	}
#	print STDERR "[DEBUG] '$refBaitSeq' refBaitSeq\n";
#	foreach (@baitAnno){
#		print STDERR "[DEBUG] '".$_->[1]."' ".$_->[0]."\n";
#	}
	
	$chr = $F[$f2i{"Rname"}];
	$start = $F[$f2i{"Rstart"}] - 1;
	$end = $F[$f2i{"Rend"}];
	$strand = $F[$f2i{"Strand"}] > 0 ? "+" : "-";
	# slop towards prey head
	$bedStart = $start;
	$bedEnd = $end;
	if($strand eq "-"){
		$bedEnd += $bait_slop;
	}else{
		$bedStart -= $bait_slop;
	}
	## retrieve Rname:(Rstart-1-bait_slop)~(Rend) "+"
	##     or   Rname:(Rstart-1)~(Rend+bait_slop) "-"
	##  length = B_Rend - B_Rstart + 1 + bait_slop
	$refPreySeq = `echo -e "$chr\t$bedStart\t$bedEnd\t\t\t$strand" | bedtools getfasta -tab -s -fi $ref_fasta_filename -bed -`;
	$refPreySeq = (split(/\t/, $refPreySeq))[1];
	$refPreySeq =~ s/[\r\n]+$//;
#	if($strand eq "-"){
#		$refSeqMain = substr($refPreySeq, 0, $end-$start);
#		$refSeqMain =~ tr/a-z/A-Z/;
#		$refSeqSlop = substr($refPreySeq, $end-$start);
#		$refSeqSlop =~ tr/A-Z/a-z/;
#		$refPreySeq = $refSeqMain . $refSeqSlop;
#	}else{
		$refSeqSlop = substr($refPreySeq, 0, $bait_slop);
		$refSeqSlop =~ tr/A-Z/a-z/;
		$refSeqMain = substr($refPreySeq, $bait_slop);
		$refSeqMain =~ tr/a-z/A-Z/;
		$refPreySeq = $refSeqSlop . $refSeqMain;
#	}
	@preyAnno = retrieve_overlapping_anno_features($chr, $bedStart, $bedEnd);
	if($strand eq "-"){
		@preyAnno = map {
			$_->[1] = reverse_anno_str($_->[1]);
			$_;
		} @preyAnno;
	}
#	print STDERR "[DEBUG] '$refPreySeq' refPreySeq\n";
#	foreach (@preyAnno){
#		print STDERR "[DEBUG] '".$_->[1]."' ".$_->[0]."\n";
#	}
	
	$readStart = 0;
	$readEnd = length($readSeq) - 1;
	if($junction_upstream >= 0){
		$readStart = $F[$f2i{"Qstart"}] - 1 - $junction_upstream;
		if($readStart < 0){
			$readStart = 0;
		}
	}
	if($junction_downstream >= 0){
		$readEnd = $F[$f2i{"Qstart"}] + $junction_downstream;
		if($readEnd > length($readSeq) - 1){
			$readEnd = length($readSeq) - 1;
		}
	}
	print_alignment($readSeq, $readStart, $readEnd, $bait_slop, $refBaitSeq, \@baitAnno, $refPreySeq, \@preyAnno, 0, @F);
	print_alignment($readSeq, $readStart, $readEnd, $bait_slop, $refBaitSeq, \@baitAnno, $refPreySeq, \@preyAnno, 1, @F);
	print "\n";
}


0;


sub get_match_mismatch_str{
	my ($query, $subjt, $caseSensitive) = @_;
	if(!defined($caseSensitive)){
		$caseSensitive = 0;
	}
	if(!$caseSensitive){
		$query =~ tr/a-z/A-Z/;
		$subjt =~ tr/a-z/A-Z/;
	}
	
	my @Q = split(//, $query);
	my @S = split(//, $subjt);
	my $i;
	my $ans = "";
	for($i=0; $i<@Q; $i++){
		if($i < @S){
			if($Q[$i] eq $S[$i]){
				$ans .= ".";
			}else{
				$ans .= $Q[$i];
			}
		}else{
			$ans .= $Q[$i];
		}
	}
	return $ans;
}

sub reverse_complement{
	my @ans = map { tr/ACGTacgt/TGCAtgca/; join("", reverse split//); } @_;
	if(@_ > 1){
		return @ans;
	}else{
		return $ans[0];
	}
}

sub print_alignment{
	my ($readSeq, $readStart, $readEnd, $baitSlop, $refBaitSeq, $baitAnnoList, $refPreySeq, $preyAnnoList, $should_rC, @F) = @_;
	
	my $outReadSeq = substr($readSeq, $readStart, $readEnd-$readStart);
	print "show read ".($readStart+1)."-$readEnd";
	if($should_rC){
		$outReadSeq = reverse_complement($outReadSeq);
		print " (reverse complement)";
	}
	print "\nread ".$outReadSeq."\n";
#	print STDERR "read'".$outReadSeq."'\n";
	
	my ($outRefSeq, $Qstart, $Qend);
	$Qstart = $F[$f2i{"B_Qstart"}] - 1;
	$Qend = $F[$f2i{"B_Qend"}] + $baitSlop;
#	$Qend = $Qstart + $F[$f2i{"B_Rend"}] - $F[$f2i{"B_Rstart"}] + $baitSlop;
#	print STDERR "$readStart $readEnd $Qstart $Qend\n";
	$outRefSeq = $refBaitSeq;
	if($should_rC){
#		print STDERR "befo'".$outRefSeq."'\n";
		if($readEnd > $Qend){
			$outRefSeq .= " " x ($readEnd - $Qend);
			$outRefSeq = reverse_complement($outRefSeq);
			$outRefSeq = substr($outRefSeq, 0, $readEnd-$readStart);
		}else{
			$outRefSeq = reverse_complement($outRefSeq);
			$outRefSeq = substr($outRefSeq, $Qend-$readEnd, $readEnd-$readStart);
		}
#		print STDERR "afte'".$outRefSeq."'\n";
	}else{
		if($readStart < $Qstart){
			$outRefSeq = (" " x ($Qstart-$readStart)) . $outRefSeq;
			$outRefSeq = substr($outRefSeq, 0, $readEnd-$readStart);
		}else{
			$outRefSeq = substr($outRefSeq, $readStart-$Qstart, $readEnd-$readStart);
		}
	}
	print "bait ".get_match_mismatch_str($outRefSeq, $outReadSeq, 1)."\n";
	
	my ($anno, $anno_title, $anno_str);
	foreach $anno (@$baitAnnoList){
		$anno_title = $anno->[0];
		$anno_str = $anno->[1];
		if($should_rC){
			if($readEnd > $Qend){
				$anno_str .= " " x ($readEnd - $Qend);
				$anno_str = reverse_anno_str($anno_str);
				$anno_str = substr($anno_str, 0, $readEnd-$readStart);
			}else{
				$anno_str = reverse_anno_str($anno_str);
				$anno_str = substr($anno_str, $Qend-$readEnd, $readEnd-$readStart);
			}
		}else{
			if($readStart < $Qstart){
				$anno_str = (" " x ($Qstart-$readStart)) . $anno_str;
				$anno_str = substr($anno_str, 0, $readEnd-$readStart);
			}else{
				$anno_str = substr($anno_str, $readStart-$Qstart, $readEnd-$readStart);
			}
		}
		print add_anno_title($anno_title, "     ".$anno_str."\n");
	}

	$Qstart = $F[$f2i{"Qstart"}] - 1 - $baitSlop;
	$Qend = $F[$f2i{"Qend"}];
	$outRefSeq = $refPreySeq;
	if($should_rC){
		if($readEnd > $Qend){
			$outRefSeq .= " " x ($readEnd - $Qend);
			$outRefSeq = reverse_complement($outRefSeq);
			$outRefSeq = substr($outRefSeq, 0, $readEnd-$readStart);
		}else{
			$outRefSeq = reverse_complement($outRefSeq);
			$outRefSeq = substr($outRefSeq, $Qend-$readEnd, $readEnd-$readStart);
		}
	}else{
		if($readStart < $Qstart){
			$outRefSeq = (" " x ($Qstart-$readStart)) . $outRefSeq;
			$outRefSeq = substr($outRefSeq, 0, $readEnd-$readStart);
		}else{
			$outRefSeq = substr($outRefSeq, $readStart-$Qstart, $readEnd-$readStart);
		}
	}
	print "prey ".get_match_mismatch_str($outRefSeq, $outReadSeq, 1)."\n";
	
	foreach $anno (@$preyAnnoList){
		$anno_title = $anno->[0];
		$anno_str = $anno->[1];
		if($should_rC){
			if($readEnd > $Qend){
				$anno_str .= " " x ($readEnd - $Qend);
				$anno_str = reverse_anno_str($anno_str);
				$anno_str = substr($anno_str, 0, $readEnd-$readStart);
			}else{
				$anno_str = reverse_anno_str($anno_str);
				$anno_str = substr($anno_str, $Qend-$readEnd, $readEnd-$readStart);
			}
		}else{
			if($readStart < $Qstart){
				$anno_str = (" " x ($Qstart-$readStart)) . $anno_str;
				$anno_str = substr($anno_str, 0, $readEnd-$readStart);
			}else{
				$anno_str = substr($anno_str, $readStart-$Qstart, $readEnd-$readStart);
			}
		}
		print add_anno_title($anno_title, "     ".$anno_str."\n");
	}
}


sub retrieve_overlapping_anno_features{
	my ($query_chr, $query_start, $query_end) = @_;
#	print STDERR "[DEBUG] query = " . join(" ", $query_chr, $query_start, $query_end)."\n";
	my @ans = ();
	my ($now_key, $F, $anno_str, $pos);
	my ($anno_chr, $anno_start, $anno_end);
	foreach $now_key (@anno_keys){
		($anno_chr, $anno_start, $anno_end) = @{$anno{$now_key}};
#		print STDERR "[DEBUG] anno  = " . join(" ", $anno_chr, $anno_start, $anno_end, $now_key)."\n";
		if($anno_chr eq $query_chr){
			if($anno_start >= $query_end || $anno_end <= $query_start){
				# no overlapping
#				print STDERR "[DEBUG] NO overlapping\n";
			}else{
				$anno_str = "";
				for($pos=$query_start; $pos<$query_end; $pos++){
					if($pos < $anno_start || $pos >= $anno_end){
						$anno_str .= " ";
					}elsif($pos==$anno_start){
						$anno_str .= "<";
					}elsif($pos==$anno_end-1){
						$anno_str .= ">";
					}else{
						$anno_str .= "-";
					}
				}
				push(@ans, [$now_key, $anno_str]);
#				print STDERR "[DEBUG] has overlapping\n";
#				print STDERR "[DEBUG] anno  = " . join(" ", $anno_chr, $anno_start, $anno_end, $now_key)."\n";
#				print STDERR "[DEBUG] $now_key '$anno_str'\n";
			}
		}
	}
	return @ans;
}

sub reverse_anno_str{
	my @ans = map { tr/<>/></; join("", reverse split//); } @_;
	if(@_ > 1){
		return @ans;
	}else{
		return $ans[0];
	}
}

sub add_anno_title{
	my ($anno_title, $anno_str) = @_;
	my @A = split(//, $anno_str);
	my @T = split(//, $anno_title);
	for($i=0; $i<@T && $i<@A; $i++){
		if($A[$i] eq " " || $A[$i] eq "-"){
			$A[$i] = $T[$i];
		}
	}
	return join("", @A);
}
	
	
	
	
