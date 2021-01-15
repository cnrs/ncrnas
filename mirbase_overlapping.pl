#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl $0 Day3_1_integrated_prediction_desc.xls mmu_mm9_circRNA.txt > mm9_circRNA.Day3_1.txt\n" unless (@ARGV == 2);


my %hash = ();

open (IN, "$ARGV[0]") || die "cannot not open $ARGV[0]\n";
while(<IN>){
	chomp;
	#Final_ID        Chr     Start   End
	#chr10:100038350|100046643       chr10   100038350       100046643
	if (/^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)/){
		my $rna = $1;
		my $c = $2;
		my $s = $3;
		my $e = $4;
		$hash{$rna}{c} = $c;
		$hash{$rna}{s} = $s;
		$hash{$rna}{e} = $e;
	}
	else{
		warn "not match: $_\n";
	}
}
close(IN);


open (IN, "$ARGV[1]") || die "cannot not open $ARGV[1]\n";
while(<IN>){
	chomp;
	#chrom start   end     strand
	#chr14   99472960        99473012        +
	if (/^(\S+)\s+(\d+)\s+(\d+)/){
		my $c = $1;
		my $s = $2;
		my $e = $3;
		
		foreach my $rna (keys %hash){
			my $rna_c = $hash{$rna}{c};
			my $rna_s = $hash{$rna}{s};
			my $rna_e = $hash{$rna}{e};
			
			next unless ($c eq $rna_c);
			
			if ($s >= $rna_s && $s <= $rna_e){
				print "$_\n";
				last;
			}
			elsif($e >= $rna_s && $e <= $rna_e){
				print "$_\n";
				last;
			}
			elsif($e >= $rna_e && $s <= $rna_s){
				print "$_\n";
				last;
			}
			
		}
	}
	else{
		warn "not match: $_\n";
	}
}
close(IN);

