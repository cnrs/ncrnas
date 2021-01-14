# circrna

bwa的使用流程
1.建立 Index
根据reference genome data(e.g. reference.fa) 建立 Index File
$ bwa index -a bwtsw human_hg18_ref.fa（human参考基因组19）
2.对reads进行mapping
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

