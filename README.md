# circrna

conda create -n circrna python=2.7

conda activate circrna

conda install samtools hisat2 bwa bowtie2 pysam numpy





bwa的使用流程
1.建立 Index
$ bwa index -a bwtsw mm9.fa

2.对reads进行mapping

cat Day3_1_1.fq  Day3_2_1.fq  Day3_3_1.fq  Day7_1_1.fq  Day7_2_1.fq  Day7_3_1.fq  WT5_1_1.fq  WT5_2_1.fq  WT5_3_1.fq > all_clean_1.fq & sleep 1s

cat Day3_1_2.fq  Day3_2_2.fq  Day3_3_2.fq  Day7_1_2.fq  Day7_2_2.fq  Day7_3_2.fq  WT5_1_2.fq  WT5_2_2.fq  WT5_3_2.fq > all_clean_2.fq & sleep 1s

bwa mem -t 36 /usr/local/db/ucsc/mouse/mm9.fa all_clean_1.fq all_clean_2.fq > all_clean_s.sam

BWA用法：
https://blog.csdn.net/weixin_43569478/article/details/108079100

CIRI:
https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip

perl CIRI.pl -I in.sam -O output.ciri -F ref.fa

perl CIRI2.pl -T 36 -I all_clean.sam -O CIRI.ciri -F /usr/local/db/ucsc/mouse/mm9.fa

bowtie2-build --threads 36 mm9.fa mm9 

bowtie2 -p 36 --very-sensitive --score-min=C,-15,0 --mm -x /usr/local/db/ucsc/mouse/mm9 -q -1 all_clean_1.fq -2 all_clean_2.fq -S all_bowtie2.sam

samtools view -hbuS all_bowtie2.sam > all_bowtie2.sam.tmp

samtools sort -o all_bowtie2.bam all_bowtie2.sam.tmp

find_circ:
https://github.com/marvin-jens/find_circ

使用文档：
https://www.cnblogs.com/yanjiamin/p/11973687.html


 bowtie2 -p 16 --very-sensitive --score-min=C,-15,0 --mm -x /path/to/bowtie2_index -q -1 reads1.fq -2 reads2.fq | samtools view -hbuS - | samtools sort - -o output.bam
 

5.挑出没有比对上的序列，各取两头20bp短序列（anchor)
1 samtools view -hf 4 output.bam | samtools view -Sb - > unmapped.bam
2 /path/to/unmapped2anchors.py unmapped.bam | gzip > anchors.fq.gz
 

6.根据anchor比对基因组情况寻找潜在的circRNA

###find_circ.py参数介绍###

#--prefix参数指定的是spliced_sites.bed文件中第四列内容的前缀，建议指定为物种对应的三字母缩写，需要注意的是，在spliced_sites_bed中同时包含了环状RNA和线性RNA,环状RNA的名称用circ标识，线性RNA的名称用norm标识，这里设置为--prefix=hsa_
#--name参数会在生成的spliced_sites.bed文件中指定tissues列的名字
#--reads参数会生成包含spliced reads的fa文件
#--stats参数会生成包含数值统计信息的txt文件
1 bowtie2 -p 16 --reorder --mm  --score-min=C,-15,0 -q -x /path/to/bowtie2_index -U anchors.fq.gz | /path/to/find_circ.py --genome=/path/to/hg38.fa --prefix=hsa_ --name=my_test_sample --stats=<run folder>/stats.txt --reads=<run folder>/splice_reads.fa > <run folder>/spliced_sites.bed
 

###根据以下规则对结果进行筛选

1.根据关键词CIRCULAR筛选环状RNA

2.去除线粒体上的环状RNA

3.筛选unique junction reads数至少为2的环状RNA

4.去除断裂点不明确的环状RNA

5.过滤掉长度大于100kb的circRNA,这里的100kb为基因组长度，直接用环状RNA的头尾相减即可

1 grep CIRCULAR spliced_sites.bed | grep -v chrM | gawk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | /path/to/maxlength.py 100000 > find_circ.candidates.bed
 

7.分析多个样本

#如果有多个样本，需要分别用find_circ.py运行，然后将各自的结果合并
1 /path/to/merge_bed.py sample1.bed sample2.bed [...] >combined.bed


bowtie2 -p 36 --very-sensitive --score-min=C,-15,0 --mm -x /usr/local/db/ucsc/mouse/mm9 -q -1 all_clean_1.fq -2 all_clean_2.fq -S all_bowtie2.sam

samtools view -hbuS all_bowtie2.sam > all_bowtie2.sam.tmp

samtools sort -o all_bowtie2.bam all_bowtie2.sam.tmp


samtools view -hf 4 all_bowtie2.bam | samtools view -Sb - > unmapped.bam

unmapped2anchors.py unmapped.bam | gzip > anchors.fq.gz

bowtie2 -p 36 --score-min=C,-15,0 --reorder --mm -q -U anchors.fq.gz -x /usr/local/db/ucsc/mouse/mm9 | find_circ.py --genome=/usr/local/db/ucsc/mouse/mm9.fa --prefix=mm9_ --name=sample --stats=find_circ/stats.txt --reads=find_circ/spliced_reads.fa > find_circ/splice_sites.bed

