#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRC.GENECOUNT.txt /usr/local/db/ucsc/mouse/mm9.fa > circRNAs.fa\n" unless (@ARGV == 2);

my %hash = ();
open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	my $s = $_;
	#circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio circRNA_type    gene_id strand  junction_reads_ID
	my $circ       = "";
	my $chromesome = "";
	#next unless ($circrna_s =~ /\d+/);
	#chr10:10452447-10452958:+
	if ($s =~ /^(\S+)\:(\d+\-\d+\:\S)\s+/){
		$chromesome = $1;
		$circ       = $2;
		$hash{$chromesome}{$circ} = 1;
	}
	else{
		warn "not match: $s\n";
	}
	#warn "$chromesome\: $circ\n";
}
close(IN);

my %genome = ();
{
	local $/ = '>';
	open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
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