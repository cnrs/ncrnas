#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;

die "Usage: perl $0 mmu_mm9_circRNA.txt CIRC.ANNO.txt > CIRC.ANNO.CIRCBASE.txt\n" unless (@ARGV == 2);

my %hash = ();

open (IN, "$ARGV[0]") || die "cannot not open $ARGV[0]\n";
while(<IN>){
	chomp;
	## chrom start   end     strand  circRNA ID      genomic length  spliced seq length      samples repeats annotation      best transcript gene symbol     circRNA study
	#chr14   99472960        99473012        +       mmu_circ_0000562        52      52      mm_ES   None    ALT_ACCEPTOR, ALT_DONOR, coding, INTERNAL, UTR3 NM_175265    6720463M24Rik   Memczak2013
	
	if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\S)\s+(\S+)\s+/){
		my $chromosome = $1;
		my $circ_s     = $2;
		my $circ_e     = $3;
		my $strand     = $4;
		my $circ_id    = $5;
		
		$hash{$circ_id}{c} = $chromosome;
		$hash{$circ_id}{s} = $circ_s;
		$hash{$circ_id}{e} = $circ_e;
		$hash{$circ_id}{r} = $strand;
	}
	else{
		warn "not match: $_\n";
	}
}
close(IN);


open (IN, "$ARGV[1]") || die "cannot not open $ARGV[1]\n";
while(<IN>){
	chomp;
	my $s = $_;
	$s =~ s/\s+$//g;
	print "$s";
	#chrom start   end     strand
	#chr10:10452447-10452958:+
	if ($s =~ /^GENEID\s+/){
		print "\tcircBase_ID";
	}
	elsif ($s =~ /^(\S+)\:(\d+)\-(\d+)\:(\S)\s+/){
		my $chromosome = $1;
		my $circ_s     = $2;
		my $circ_e     = $3;
		my $strand     = $4;
		my $circ_l = $circ_e - $circ_s + 1;
		
		my $circ = "";
		my %tmp = ();
		foreach my $circ_id (keys %hash){
			my $c = $hash{$circ_id}{c};
			my $s = $hash{$circ_id}{s};
			my $e = $hash{$circ_id}{e};
			my $r = $hash{$circ_id}{r};
			
			my $base_l = $e - $s + 1;
			my $min_l = min ($circ_l, $base_l);
			
			next unless ($chromosome eq $c && $strand eq $r);
			next if ($s > $circ_e || $circ_s > $e);
			my @array = sort {$a <=> $b} ($s, $e, $circ_s, $circ_e);
			my $overlap_l = $array[2] - $array[1] + 1;
			
			next unless ($overlap_l >= $min_l * 0.8);
			$tmp{$circ_id} = 1;
			#$circ .= "$circ_id\;";
		}
		
		$circ = join (";", sort keys %tmp);
		if ($circ eq ""){
			$circ = "NA";
		}
		
		print "\t$circ";
	}
	else{
		warn "not match: $_\n";
	}
	
	print "\n";
}
close(IN);

