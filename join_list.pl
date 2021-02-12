#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 DEG.DAY3_vs_WT5.txt targets.tab > DAY3_vs_WT5.MIR_TARGET.txt\n" unless (@ARGV == 2);

my %hash = ();


open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	#chr10:114826229-114827492:+|CIRC_000720 117.770808056798        45.0370542380077        95.9910566952797        74.828967385195 76.5842298302849        117.070529710647     36.3918143887116        38.0178607956739        41.3491000463199        71.4490467941019        1.19712541020518        0.363417181808364   3.29408038510531 0.000987442551384088    0.0744173316677463      intron (NM_001033261, intron 1 of 34)   intron (NM_001033261, intron 1 of 34)   4846    NM_001033261 216345  Mm.207621       NM_001033261    ENSMUSG00000034163      Zfc3h1  BC033596|Ccdc131|Psrc2  zinc finger, C3H1-type containing       protein-coding       NA
	#chr10:104954369-104954678:+     4.80839329473504        0.0221094279605325      NA
	#chio log2FC pvalue circbaseid
	my @array = split (/\t/, $_);
	my $circ = $array[0];
	my $lgfc = $array[11];
	my $gene = $array[24];
	my $type = $array[27];
	my $hits = $array[28];
	$circ =~ s/\|\S+$//g;
	
	$hash{$circ} = "$array[0]\t$lgfc\t$gene\t$type\t$hits";
}
close(IN);

print "miRNA\tcircRNA\tlog2FC\tsymbol\ttype\tcircbase\n";
open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	#mmu-let-7a-2-3p chr10:41916271-41917707:-       161.00  -25.72
	if (/^(\S+)\s+(\S+)\s+/){
		my $mirna = $1;
		my $circ  = $2;
		if (exists $hash{$circ}){
			print "$mirna\t$hash{$circ}\n";
		}
	}
}
close(IN);
