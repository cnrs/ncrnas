#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRI.ciri FIND_CIRC.candidates.bed > CIRC.bed\n" unless (@ARGV == 2);

my %hash = ();
my $lens_cut = 3000;
my $cnt = 0;

open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	#circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio circRNA_type    gene_id strand  junction_reads_ID
	my @array = split (/\t/, $_);
	#my $circrna_id = $array[0];
	my $chromesome = $array[1];
	my $circrna_s  = $array[2];
	my $circrna_e  = $array[3];
	my $strand     = $array[10];
	next unless ($circrna_s =~ /\d+/);
	next if ($circrna_e - $circrna_s + 1> $lens_cut);
	
	#my $circ_pos = "$chromesome\:$strand";
	push (@{$hash{$chromesome}}, [$circrna_s, $circrna_e, $strand]);
	#warn "$chromesome, $circrna_s, $circrna_e, $strand\n";
	$cnt ++;
}
close(IN);

warn "$cnt\n";
$cnt = 0;

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
	next if ($circrna_e - $circrna_s + 1> $lens_cut);
	
	my $circ_pos = "$chromesome\:$strand";
	push (@{$hash{$chromesome}}, [$circrna_s, $circrna_e, $strand]);
	#warn "$chromesome, $circrna_s, $circrna_e, $strand\n";
	$cnt ++;
}
close(IN);

warn "$cnt\n";

# merge circRNAs
# cover more that 90% sequences of shorter seuqences
my %circs = ();
my @tmp = ();

my @array = ("+", "-");
foreach my $strand (@array) {
	foreach my $chromesome (keys %hash) {
		foreach my $e (sort {$$a[0] <=> $$b[0]} @{$hash{$chromesome}}) {
			my ($circrna_s, $circrna_e, $str) = ($$e[0], $$e[1], $$e[2]);
			next unless ($strand eq $str);
			my $circ_id = "$circrna_s\-$circrna_e\:$strand";
			#warn "$circ_id\n";
			unless (exists $circs{$chromesome}){
				$circs{$chromesome}{$circ_id} = 1;
			}
			
			my $flag = 0;
			foreach my $cis (keys %{$circs{$chromesome}}){
				my ($c_s, $c_e, $c_r) = $cis =~ /^(\d+)\-(\d+)\:(\S)/;
				next unless ($strand eq $c_r);
				
				if (abs($circrna_s - $c_s) <= 10 && abs($circrna_e - $c_e) <= 10){
					$flag = 1;
					last;
				}
			}
			
			if ($flag == 0){
				$circs{$chromesome}{$circ_id} = 1;
			}
		}
	}
}

print "chromesome\tcirc_start\tcirc_end\tcircrna_id\tlength\tstrand\n";
foreach my $chromesome (keys %circs) {
	foreach my $e (keys %{$circs{$chromesome}}) {
		my ($circ_s, $circ_e, $strand) = $e =~ /^(\d+)\-(\d+)\:(\S)$/;
		my $circrna_id = "$chromesome\:$circ_s\-$circ_e\:$strand";
		my $lens = ($circ_e - $circ_s +1) / 1000;
		
		print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\t$lens\t$strand\n";
		#sleep (1);
	}
}

