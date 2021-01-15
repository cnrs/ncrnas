#!/usr/bin/perl -w
use warnings;
use strict;

die "Usage: perl $0 KX683219.1.fasta list.txt > s.fa\n" unless (@ARGV == 2);

my %hash = ();

{
	local $/ = '>';
	open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
	while(<IN>){
		chomp;
		my ($name, $sequence) = split (/\n/, $_, 2);
		next unless ($name && $sequence);
		
		$name =~ s/\s+.*$//g;
		$sequence =~ s/\s+//g;
			
		my $lens = length ($sequence);
		$hash{$name}{s} = $sequence;
		$hash{$name}{l} = $lens;
	}
	close(IN);
}

open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	#chr17   91100083        91103707        -       mmu_circ_0000835
	if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
		my $s_c = $1;
		my $s_s = $2;
		my $s_e = $3;
		my $strand = $4;
		my $s_n = $5;
		
		my $lens = $hash{$s_c}{l};
		my $sequence = $hash{$s_c}{s};
		warn "$s_e > $lens\n" if ($s_e > $lens);
		if ($s_e > $lens){
			$s_e = $lens;
		}
		
		my $new_s = substr ($sequence, $s_s - 1, $s_e - $s_s + 1);
		if ($strand eq "-"){
			$new_s =~ tr/atcgATCG/tagcTAGC/;
			$new_s = reverse ($new_s);
		}
		
		print ">$s_n\t$_\n$new_s\n";
		
	}
}
close(IN);




