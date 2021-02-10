#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min sum maxstr minstr shuffle/;


die "Usage: perl $0 CIRC.sam CIRC.bed > CIRC.bed.good\n" unless (@ARGV == 2);

my %hash = ();
my @array = ();
my $cov_n = 2;
my $cov_p = 0.8;

open (IN, "$ARGV[0]") || die "cannot open $ARGV[0]\n";
while(<IN>){
	chomp;
	#E00603:107:HLKG3CCXY:2:1202:16193:48634_A__GATGGAACTGGTGGGCAACCAGAGTGGACCTTACAGAGGGTGACAATGGAAAACGTGAAGAGCAAGAAGACCCTACATTTTGTTGCAAATGTGTGTCTTT 16      chr1    4153955 32 20M      *       0       0       GGTTGCCCACCAGTTCCATC    iiiiiiiiiiiiiiiee```    AS:i:0  XS:i:-6 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:20 YT:Z:UU
	#bowtie2
	#QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL
	
	#https://blog.csdn.net/weixin_42670653/article/details/109185212
	#https://www.cnblogs.com/yahengwang/p/10695981.html
	#1�� ��paired-end����mate pair��һ��
	#2�� ˫ĩ�˱ȶԵ�һ��
	#4�� û�жԱ��ϲο�����
	#8����paired-end����mate pair��һ�������޷��ȶ��ϲο�������
	#16���ȶԵ��ο����еĸ�����
	#32�� ˫ĩ��reads����һ���ȶԵ��ο����еĸ�����
	#64������reads��mate1
	#128�� ����reads��mate2
	
	#	�ȶ���11���кͿ�ѡ�еĽ���
	#1  QNAME  �ȶԵ�������
	#2  FLAG   Bwise FLAG(�����ȶ����ͣ�pairing��strand��mate strand��)
	#3  RNAME  �ȶ��ϵĲο�������
	#4  POS    1-Based�ıȶ��ϵ�����ߵĶ�λ
	#5  MAPQ   �ȶ�����
	#6  CIGAR  Extended CIGAR string (��������MIDNSHP) �ȶԽ����Ϣ��ƥ���������ɱ���ӵȡ�
	#7  MRNM   ��ƥ�������һ�����У��ȶ��ϵĲο�������
	#8  MPOS   1-Based leftmost Mate POsition
	#9  ISIZE  ����Ƭ�γ���
	#10 SEQ    �Ͳο�������ͬһ�����ϵıȶ�����(���ȶԽ���ڸ������ϣ����������䷴���ظ�����)
	#11 QUAL   �ȶ����е�����(ASCII-33=Phred base quality)
	#12 ��ѡ���У���TAG��TYPE��VALUE����ʽ�ṩ�������Ϣ

	#SN1077:321:H2NH5BCXX:1:1103:1329:2233   99      S51     126036  60      133M    =       126036  -133    AAACCTTCGGGGCCTTGACTAAATTATCTTCTACGATCAAAATATCGTGGACTTCTGTCAGAATGTCCATCATATGGCCTGCGTTCATGAAAGAACTCTCTTGGAGGTGGATATCCGCCTGGCAAAGGAGAAT      DDDDDIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIIIIIIIIIHIIIIIIIIIIIIIHIIIIHIIIIIIIIIIGHIIIIIIIIIIIIIHIHFDHHIIIIIHHIIIIIIIIIHIIIHIIIIHII      AS:i:-12        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:13T10G108  YS:i:-12        YT:Z:CP NH:i:1
	#SN1077:321:H2NH5BCXX:1:1103:1329:2233   147     S51     126036  60      133M    =       126036  -133    AAACCTTCGGGGCCTTGACTAAATTATCTTCTACGATCAAAATATCGTGGACTTCTGTCAGAATGTCCATCATATGGCCTGCGTTCATGAAAGAACTCTCTTGGAGGTGGATATCCGCCTGGCAAAGGAGAAT      HHFHHIHIIIIIIIIIIIIIIIIIIIIIIIHFIIIIIIIIIIHHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHIHGIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDDDDD      AS:i:-12        XN:i:0  XM:i:2  XO:i:0  XG:i:0  NM:i:2  MD:Z:13T10G108  YS:i:-12        YT:Z:CP NH:i:1

	next if (/^\@HD\s+/);
	next if (/^\@SQ\s+/);
	next if (/^\@PG\s+/);
	
	@array = ();
	@array = split(/\t/, $_);
	my $flag       = $array[1];
	my $chromosome = $array[2];
	next if ($chromosome eq "*");
	
	my $c_s        = $array[3];
	my $insert_size = $array[8]; #bowtie2
	next if ($insert_size < 0);
	
	my $c_e = $c_s + $insert_size;
	
	#my $strand = "+";
	#if ($flag == 0){
	#	$strand = "-";
	#}
	#elsif ($flag == 16){
	#	$strand = "-";
	#}
	
	for (my $i = $c_s; $i <= $c_e; $i ++){
		unless (exists $hash{$chromosome}{$i}){
			$hash{$chromosome}{$i} = 1;
		}
		else{
			$hash{$chromosome}{$i} ++;
		}
	}
}
close(IN);

my $cnt = 0;
open (IN, "$ARGV[1]") || die "cannot open $ARGV[1]\n";
while(<IN>){
	chomp;
	#chromosome      circ_start      circ_end        circrna_id      length  strand
	#chr1    10035157        10037607        chr1:10035157-10037607:+|CIRC_000001    2.451   +
	
	my @array = split (/\t/, $_);
	my $chromosome = $array[0];
	my $circrna_s  = $array[1];
	my $circrna_e  = $array[2];
	my $strand     = $array[5];
	next unless ($circrna_s =~ /\d+/);
	
	my $c_num = 0;
	for (my $i = $circrna_s; $i <= $circrna_e; $i ++){
		unless (exists $hash{$chromosome}{$i}){
			$hash{$chromosome}{$i} = 0;
		}
		if ($hash{$chromosome}{$i} >= $cov_n){
			$c_num ++;
		}
	}
	
	if ($c_num / ($circrna_e - $circrna_s + 1) >= $cov_p){
		$cnt ++;
		#my $lens = ($circrna_e - $circrna_s + 1) / 1000;
		my $score = sprintf ("%.2f", $c_num / ($circrna_e - $circrna_s + 1));
		my $circrna_id = "$chromosome\:$circrna_s\-$circrna_e\:$strand";
		my $CIRCID = sprintf ("%06d", $cnt);
		print "$chromosome\t$circrna_s\t$circrna_e\t$circrna_id\|CIRC\_$CIRCID\t$score\t$strand\n";
	}
}
close(IN);
