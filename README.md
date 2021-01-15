# circrna

conda create -n circrna python=2.7

conda activate circrna

conda install samtools hisat2 bwa bowtie2





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



bowtie2 -p 36 --very-sensitive --score-min=C,-15,0 --mm -x /usr/local/db/ucsc/mouse/mm9 -q -1 all_clean_1.fq -2 all_clean_2.fq -S all_bowtie2.sam

samtools view -hbuS all_bowtie2.sam > all_bowtie2.sam.tmp

samtools sort -o all_bowtie2.bam all_bowtie2.sam.tmp

find_circ.py --genome=/usr/local/db/ucsc/mouse/mm9.fa --prefix=mm9_ --name=my_test_sample --stats=find_circ/stats.txt --reads=find_circ/spliced_reads.fa > find_circ/splice_sites.bed




