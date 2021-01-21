#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRI.ciri FIND_CIRC.candidates.bed /usr/local/db/ucsc/mouse/mm9.fa > circRNAs.fa\n" unless (@ARGV == 3);

my %hash = ();
open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	#circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio circRNA_type    gene_id strand  junction_reads_ID
	my @array = split (/\t/, $_);
	my $chromesome = $array[1];
	my $circrna_s  = $array[2];
	my $circrna_e  = $array[3];
	my $strand     = $array[10];
	next unless ($circrna_s =~ /\d+/);
	
	my $circ = "$circrna_s\-$circrna_e\:$strand";
	$hash{$chromesome}{$circ} = 1;
}
close(IN);

open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	# chrom start   end     name    n_reads strand  n_uniq  uniq_bridges    best_qual_left  best_qual_right tissues tiss_counts  edits   anchor_overlap  breakpoints     signal  strandmatch     category
	my @array = split (/\t/, $_);
	my $chromesome = $array[0];
	my $circrna_s  = $array[1];
	my $circrna_e  = $array[2];
	my $strand     = $array[5];
	next unless ($circrna_s =~ /\d+/);
	
	my $circ = "$circrna_s\-$circrna_e\:$strand";
	$hash{$chromesome}{$circ} = 1;
}
close(IN);

my %genome = ();
{
	local $/ = '>';
	open (IN, "$ARGV[2]") || die "cannot open $ARGV[2]\n";
	while(<IN>){
		chomp;
		my ($name, $sequence) = split (/\n/, $_, 2);
		next unless ($name && $sequence);
		$name =~ s/\s+.*$//g;
		$sequence =~ s/\s+//g;
		$genome{$name} = $sequence;
	}
	close(IN);
}

foreach my $chromesome (keys %hash) {
	my $sequence = $genome{$chromesome};
	foreach my $circ (keys %{$hash{$chromesome}}) {
		my ($circ_s, $circ_e, $strand) = $circ =~ /^(\d+)\-(\d+)\:(\S)$/;
		my $circrna_id = "$chromesome\:$circ_s\-$circ_e\:$strand";
		my $lens = ($circ_e - $circ_s + 1) / 1000;
		my $new_s = substr($sequence, $circ_s - 1, $circ_e - $circ_s + 1);
		print "\>$circrna_id\t$chromesome\t$circ_s\t$circ_e\t$circrna_id\t$lens\t$strand\n$new_s\n";
	}
}

# cd-hit-est -i circRNAs.fa -o circRNAs.C90.fa -c 0.9 -G 1 -M 64000 -T 36 