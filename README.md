# circrna

bwa的使用流程
1.建立 Index
$ bwa index -a bwtsw mm9.fa

2.对reads进行mapping
cat Day3_1_1.fq  Day3_2_1.fq  Day3_3_1.fq  Day7_1_1.fq  Day7_2_1.fq  Day7_3_1.fq  WT5_1_1.fq  WT5_2_1.fq  WT5_3_1.fq > all_clean_1.fq & sleep 1s
cat Day3_1_2.fq  Day3_2_2.fq  Day3_3_2.fq  Day7_1_2.fq  Day7_2_2.fq  Day7_3_2.fq  WT5_1_2.fq  WT5_2_2.fq  WT5_3_2.fq > all_clean_2.fq & sleep 1s
bwa mem /usr/local/db/ucsc/mouse/mm9.fa all_clean_1.fq all_clean_2.fq > all_clean.sam

https://blog.csdn.net/weixin_43569478/article/details/108079100
1. BWA-backtrack 算法
对应的子命令为aln/samse/sample
单端数据用法如下：

bwa aln ref.fa reads.fq > aln_sa.sai
bwa samse ref.fa aln_sa.sai reads.fq > aln-se.sam

双端数据用法如下：

bwa aln ref.fa read1.fq > aln1_sa.sai
bwa aln ref.fa read2.fq > aln2_sa.sai
bwa sampe ref.fa aln_sa1.sai aln_sa2.sai read1.fq read2.fq > aln-pe.sam

2. BWA-SW 算法
对应的子命令为bwasw, 基本用法如下

bwa bwasw ref.fa reads.fq > aln-se.sam
bwa bwasw ref.fa read1.fq read2.fq > aln-pe.sam

3. BWA-MEM` 算法
对应的子命令为mem, 基本用法如下

bwa mem ref.fa reads.fq > aln-se.sam
bwa mem ref.fa read1.fq read2.fq > aln-pe.sam



如果是pair-end 数据（leftRead.fastq和rightRead.fastq）两个文件分别处理
$ bwa aln reference.fa leftRead.fastq leftRead.sai
$ bwa aln reference.fa rightRead.fastq rightRead.sai
如果是single-end 数据，则
$ bwa aln reference.fa singleRead.fastq singleRead.sai
3.将mapping输出文件*.sai处理成*.sam
如果是pair-end数据
$ bwa sampe -f pair-end.sam reference.fa leftRead.sai rightRead.sai leftRead.fastq rightread.fastq
如果是single-end数据
$ bwa samse -f single.sam reference.fa single.sai single.fastq


https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip


perl CIRI.pl -I in.sam -O output.ciri -F ref.fa

