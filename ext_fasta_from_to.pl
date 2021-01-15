#!/usr/bin/perl -w
use warnings;
use strict;

die "Usage: perl $0 KX683219.1.fasta 1 138767 + > s.fa\n" unless (@ARGV == 4);

my %hash = ();
my $s_s = $ARGV[1];
my $s_e = $ARGV[2];
my $strand = $ARGV[3];

local $/ = '>';
open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	my ($name, $sequence) = split (/\n/, $_, 2);
	next unless ($name && $sequence);
	
	$name =~ s/\s+.*$//g;
	$sequence =~ s/\s+//g;
	#warn "$name\n";
	#next unless (exists $hash{$name});
	
	my $lens = length ($sequence);
	warn "$s_e > $lens\n" if ($s_e > $lens);
	if ($s_e > $lens){
		$s_e = $lens;
	}
	
	my $new_s = substr ($sequence, $s_s - 1, $s_e - $s_s + 1);
	if ($strand eq "-"){
		$new_s =~ tr/atcgATCG/tagcTAGC/;
		$new_s = reverse ($new_s);
	}
	
	print ">$name\_$s_s\-$s_e\:$strand\n$new_s\n";
}
close(IN);
