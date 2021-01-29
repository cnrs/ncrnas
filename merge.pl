#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: perl $0 CIRC.bed.anno.xls GENECOUNTS.DAY3_vs_WT5.txt > GENECOUNTS.DAY3_vs_WT5.ANNO.txt\n" unless (@ARGV == 2);

my %hash = ();
my @array = ();
my $cnt = 0;

open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	my $s = $_;
	#circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio circRNA_type    gene_id strand  junction_reads_ID
	@array = ();
	@array = split (/\t/, $s);
	
	if ($cnt == 0){
		$cnt = scalar (@array);
		warn "cnt: $cnt\n";
		#sleep (3);
		for (my $i = 0; $i < $cnt; $i ++){
			warn "$i: $array[$i]\n";
		}
		sleep (1);
	}
	
	my $circ = $array[0];
	next if ($circ =~ s/^\s+$//g);
	
	my $anno = "";
	my $flag = 0;
	for (my $i = 7; $i < $cnt; $i ++){
		#chr8:71701654-71701920:+        chr8    71701655        71701920        +       0.267   NA      intron (NM_001370971, intron 1 of 3)    intron (NM_001370971, intron 1 of 3) 6310    NM_001370971
		#chr8:71701654-71701920:+        chr8    71701655        71701920        +       0.267   NA      intron (NM_001370971, intron 1 of 3)    intron (NM_001370971, intron 1 of 3) 6310    NM_199364       211134  Mm.409446       NM_199364       ENSMUSG00000036306      Lzts1   F37|FEZ1|PSD-Zip70      leucine zipper, putative tumor suppressor 1  protein-coding
		# ¶à×¢ÊÍ¼¸±é
		# perl /usr/local/.prog/anaconda/envs/chipseq/bin/annotatePeaks.pl CIRC.bed mm9 > CIRC.bed.anno.xls.1
		if ($array[$i] eq ""){
			$array[$i] = "NA";
			$flag = 1;
			warn "$i:    $array[$i]\n";
		}
		#warn "$i: $array[$i]\n";
		if ($i == $cnt - 1){
			$anno .= "$array[$i]";
		}
		else{
			$anno .= "$array[$i]\t";
		}
	}
	if ($flag == 1){
		#warn "######### $circ\t$anno\n";
		#next;
	}
	#$anno =~ s/\s+$//g;
	$hash{$circ} = $anno;
	#warn "$circ\t$anno\n";
	#sleep (1);
}
close(IN);

open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	my $s = $_;
	$s =~ s/\s*$//g;
	#GENEID  Day3_1  Day3_2  Day3_3  Day7_1  Day7_2  Day7_3  WT5_1   WT5_2   WT5_3   baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
	if ($s =~ /^(\S+)\s+(\S+.*\S+)$/){
		my $circ = $1;
		if (exists $hash{$circ}){
			print "$s\t$hash{$circ}\n";
		}
		else{
			warn "not exists: $circ\n";
		}
	}
	else{
		warn "not match: $_\n";
	}
}
close(IN);

# cnt: 19
# 0: PeakID (cmd=annotatePeaks.pl CIRC.bed mm9)
# 1: Chr
# 2: Start
# 3: End
# 4: Strand
# 5: Peak Score
# 6: Focus Ratio/Region Size
# 7: Annotation
# 8: Detailed Annotation
# 9: Distance to TSS
# 10: Nearest PromoterID
# 11: Entrez ID
# 12: Nearest Unigene
# 13: Nearest Refseq
# 14: Nearest Ensembl
# 15: Gene Name
# 16: Gene Alias
# 17: Gene Description
# 18: Gene Type

# grep 'chr10:24315446-24317824:-' CIRC.bed.anno.xls.1 
# grep 'chr10:24363384-24365084:+' CIRC.bed.anno.xls.1
# grep 'chr11:5093213-5093519:-' CIRC.bed.anno.xls.1
# grep 'chr14:8725386-8727201:+' CIRC.bed.anno.xls.1
# grep 'chr14:8747774-8748617:+' CIRC.bed.anno.xls.1
# grep 'chr14:8761514-8763660:+' CIRC.bed.anno.xls.1
# grep 'chr14:8767000-8768498:-' CIRC.bed.anno.xls.1
# grep 'chr14:8783418-8784060:+' CIRC.bed.anno.xls.1
# grep 'chr16:33461306-33461604:+' CIRC.bed.anno.xls.1
# grep 'chr16:33465332-33468269:+' CIRC.bed.anno.xls.1
# grep 'chr16:33502093-33502442:-' CIRC.bed.anno.xls.1
# grep 'chr17:35219772-35220449:-' CIRC.bed.anno.xls.1
# grep 'chr18:40684972-40686533:-' CIRC.bed.anno.xls.1
# grep 'chr18:40702037-40702494:+' CIRC.bed.anno.xls.1
# grep 'chr7:107072208-107073860:-' CIRC.bed.anno.xls.1
# grep 'chr7:107103624-107104050:-' CIRC.bed.anno.xls.1
# grep 'chr8:12866826-12867330:+' CIRC.bed.anno.xls.1
# grep 'chr8:87295682-87296872:-' CIRC.bed.anno.xls.1
# grep 'chr8:87311299-87312045:-' CIRC.bed.anno.xls.1
# grep 'chr9:122804391-122804718:+' CIRC.bed.anno.xls.1
