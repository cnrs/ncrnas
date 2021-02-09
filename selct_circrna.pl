#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRI.ciri FIND_CIRC.candidates.bed > CIRC.bed\n" unless (@ARGV == 2);

my %hash = ();
my %circs = ();
my $lens_cut = 5000;
my $adjacent_distance = 20;
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
	next if ($circrna_e - $circrna_s + 1 > $lens_cut);
	
	my $circ = "$circrna_s\-$circrna_e";
	push (@{$hash{$chromesome}{$strand}}, [$circrna_s, $circrna_e]);
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
	next if ($circrna_e - $circrna_s + 1 > $lens_cut);
	
	& find_overlapping ($circrna_s, $circrna_e);
	$cnt ++;
}
close(IN);
warn "$cnt\n";
$cnt = 0;

& find_unis ();

sub find_overlapping {
	my ($c_s, $c_e) = @_;
	foreach my $chromesome (keys %hash) {
		foreach my $strand (keys %{$hash{$chromesome}}) {
			foreach my $e (sort {$$a[0] <=> $$b[0]} @{$hash{$chromesome}{$strand}}) {
				my ($circrna_s, $circrna_e) = ($$e[0], $$e[1]);
				
				my $circ = "$circrna_s\-$circrna_e";
				
				if (abs($circrna_s - $c_s) <= $adjacent_distance && abs($circrna_e - $c_e) <= $adjacent_distance){
					$circs{$chromesome}{$strand}{$circ} = 1;
				}
			}
		}
	}
}

sub find_unis {
	my %uni = ();
	$cnt = 0;
	print "chromesome\tcirc_start\tcirc_end\tcircrna_id\tlength\tstrand\n";
	foreach my $chromesome (sort keys %circs) {
		foreach my $strand (sort keys %{$circs{$chromesome}}) {
			foreach my $circ (sort keys %{$circs{$chromesome}{$strand}}) {
				my ($circ_s, $circ_e) = $circ =~ /^(\d+)\-(\d+)$/;
				my $circrna_id = "$chromesome\:$circ_s\-$circ_e\:$strand";
				my $lens = ($circ_e - $circ_s + 1) / 1000;
				
				unless (exists $uni{$chromesome}{$strand}){
					$uni{$chromesome}{$strand}{$circ} = 1;
					$cnt ++;
					my $CIRCID = sprintf ("CIRC_06d%", $cnt);
					print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\|$CIRCID\t$lens\t$strand\n";
				}
				else{
					my $flag = 0;
					foreach my $e (sort keys %{$uni{$chromesome}{$strand}}){
						next if ($e eq $circ);
						
						my ($c_s, $c_e) = $e =~ /^(\d+)\-(\d+)$/;
						if (abs($circ_s - $c_s) <= $adjacent_distance && abs($circ_e - $c_e) <= $adjacent_distance){
							$flag = 1;
							last;
						}
					}
					
					if ($flag == 0){
						$uni{$chromesome}{$strand}{$circ} = 1;
						$cnt ++;
						my $CIRCID = sprintf ("CIRC_06d%", $cnt);
						print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\|$CIRCID\t$lens\t$strand\n";
					}
				}
				
				#print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\t$lens\t$strand\n";
				#sleep (1);
			}
		}
	}
}


