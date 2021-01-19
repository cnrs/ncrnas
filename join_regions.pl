#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl $0 CIRI.ciri FIND_CIRC.candidates.bed > CIRC.txt\n" unless (@ARGV == 2);

my %hash = ();

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
	
	
	my $circ_pos = "$chromesome\:$strand";
	push (@{$hash{$circ_pos}}, [$circrna_s, $circrna_e])
	
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
	
	my $circ_pos = "$chromesome\:$strand";
	push (@{$hash{$circ_pos}}, [$circrna_s, $circrna_e, $chromesome, $strand])
}
close(IN);

my %circs = ();
print "chromesome\tcirc_start\tcirc_end\tcircrna_id\tlength\tstrand\n";
foreach my $circ_pos (keys %hash) {
	my $circ_s = 0;
	my $circ_e = 0;
	my ($chromesome, $strand) = $circ_pos =~ /^(\S+)\:(\S)$/;
	my $circrna_id = "";
	my $lens = 0;
	
	foreach my $ele (sort {$$a[0] <=> $$b[0]} @{$hash{$circ_pos}}){
		my $p_s = $$ele[0];
		my $p_e = $$ele[1];
		#my $chromesome = $$ele[2];
		#my $strand     = $$ele[3];
		
		if ($circ_s == 0){
			$circ_s = $p_s;
			$circ_e = $p_e;
		}
		
		if ($circ_e >= $p_s){
			if ($circ_e < $p_e){
				$circ_e = $p_e;
			}
		}
		elsif($circ_e < $p_s){
			$circrna_id = "$chromesome\:$circ_s\-$circ_e\:$strand"; 
			$circs{$circrna_id}{chromesome} = $chromesome;
			$circs{$circrna_id}{circ_s}     = $circ_s;
			$circs{$circrna_id}{circ_e}     = $circ_e;
			$circs{$circrna_id}{strand}     = $strand;
			$lens = int (($circ_e - $circ_s +1) / 1000);
			print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\t$lens\t$strand\n";
			#sleep (1);
			$circ_s = $p_s;
			$circ_e = $p_e;
		}
		#warn "$circ_pos\t$p_s\t$p_e\n";
		
	}
	
	$circrna_id = "$chromesome\:$circ_s\-$circ_e\:$strand"; 
	$circs{$circrna_id}{chromesome} = $chromesome;
	$circs{$circrna_id}{circ_s}     = $circ_s;
	$circs{$circrna_id}{circ_e}     = $circ_e;
	$circs{$circrna_id}{strand}     = $strand;
	$lens = int (($circ_e - $circ_s +1) / 1000);
	print "$chromesome\t$circ_s\t$circ_e\t$circrna_id\t$lens\t$strand\n";
	#sleep (1);
}

