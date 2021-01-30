#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 DAY3_vs_WT5.CIRCBASE.txt targets.tab > DAY3_vs_WT5.MIR_TARGET.txt\n" unless (@ARGV == 2);

my %hash = ();


open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	
	#chr10:104954369-104954678:+     4.80839329473504        0.0221094279605325      NA
	#chio log2FC pvalue circbaseid
	if (/^(\S+)\s+(\S+\s+\S+\s+\S+)/){
		$hash{$1} = $2;
	}
}
close(IN);

open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	#mmu-let-7a-2-3p chr10:41916271-41917707:-       161.00  -25.72
	if (/^\S+\s+(\S+)\s+/){
		my $circ = $1;
		if (exists $hash{$circ}){
			print "$_\t$hash{$circ}\n";
		}
	}
}
close(IN);
