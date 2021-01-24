#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRI.ciri FIND_CIRC.candidates.bed > CIRC.txt\n" unless (@ARGV == 2);

my %hash = ();
my $lens_cut = 5000;

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
	next if ($circrna_e - $circrna_s + 1> $lens_cut);
	
	my $circ_pos = "$chromesome\:$strand";
	push (@{$hash{$chromesome}}, [$circrna_s, $circrna_e, $strand]);
	#warn "$chromesome, $circrna_s, $circrna_e, $strand\n";
}
close(IN);

# merge circRNAs
# cover more that 90% sequences of shorter seuqences
my %circs = ();
my @tmp = ();

my @array = ("+", "-");
foreach my $strand (@array) {
	foreach my $chromesome (keys %hash) {
		my $circ_s_min = 0;
		my $circ_e_max = 0;
		
		foreach my $e (sort {$$a[0] <=> $$b[0]} @{$hash{$chromesome}}){
			my ($circrna_s, $circrna_e, $str_s) = ($$e[0], $$e[1], $$e[2]);
			# warn "$chromesome, $circrna_s, $circrna_e, $strand, $str_s\n";
			next unless ($strand eq $str_s);
			
			
			if ($circ_s_min == 0){
				$circ_s_min = $circrna_s;
			}
			if ($circ_e_max == 0){
				$circ_e_max = $circrna_e;
			}
			
			if ($circrna_s >= $circ_s_min && $circrna_s <= $circ_e_max){
				$circ_e_max = max($circ_e_max, $circrna_e);
			}
			
			if($circ_e_max < $circrna_s){
				my $circ_id = "$circ_s_min\-$circ_e_max\:$strand";
				#warn "$circ_id\n";
				$circs{$chromesome}{$circ_id} = 1;
				$circ_s_min = $circrna_s;
				$circ_e_max = $circrna_e;
			}
		}
		unless ($circ_s_min == 0){
			my $circ_id = "$circ_s_min\-$circ_e_max\:$strand";
			#warn "$circ_id\n";
			$circs{$chromesome}{$circ_id} = 1;
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

