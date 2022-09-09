# Biochemical properties and progress in cancers of tRNA-derived fragments

```
conda create -n ncrnas
conda activate ncrnas
conda install curl=7.64.0 cutadapt=1.18 fastqc=0.11.9 ghostscript=9.22 libgit2=0.27.8 numpy=1.15.4 openssl=1.0.2u pandoc=2.10 pip=20.1.1 pigz=2.3.4 python=2.7.15 r=3.5.1 r-boot=1.3_20 r-checkpoint=0.4.4 r-class=7.3_14 r-cluster=2.0.7_1 r-codetools=0.2_15 r-curl=3.2 r-deployrrserve=9.0.0 r-doparallel=1.0.13 r-foreach=1.5.0 r-foreign=0.8_70 r-iterators=1.0.10 r-jsonlite=1.5 r-kernsmooth=2.23_15 r-lattice=0.20_35 r-mass=7.3_49 r-matrix=1.2_14 r-mgcv=1.8_23 r-microsoftr=3.5.0.108 r-nlme=3.1_137 r-nnet=7.3_12 r-png=0.1_7 r-r6=2.2.2 r-recommended=3.5.1 r-revoioq=10.0.0 r-revomods=11.0.0 r-revoutils=11.0.0 r-revoutilsmath=11.0.0 r-rpart=4.1_13 r-runit=0.4.26 r-spatial=7.3_11 r-survival=2.41_3 samtools=1.9 star=2.5.2b subread=2.0.1=hed695b0_0 trim-galore=0.4.1 
```

https://pubmed.ncbi.nlm.nih.gov/31674076/

tsRNAsearch
https://github.com/GiantSpaceRobot/tsRNAsearch

nextflow run tsRNAsearch --species mouse --input_dir tsRNAsearch/ExampleData --output_dir Results



tRF-seed-match:
https://github.com/ritianjiang/tRF-seed-match
tRNA-derived fragment is a kind of small-RNA, which may play the same role with miRNA. Here we proposed this repositry to perform the seed matching for tRFs to 3'UTR.

SVM_GA_tRF_targets:
https://github.com/cmuxiaoqiong/SVM_GA_tRF_targets
predicting the targets of tRNA-derived fragments (tRFs) by SVM and GA

ca-tRF
https://github.com/Load-Star/ca-tRF
the datasets and source code for investigating cancer-asssociated tRNA-derived fragments

tRF-Glu49-cervical-cancer
https://github.com/chenyouguosoochow/tRF-Glu49-cervical-cancer
tRNA-derived fragment tRF-Glu49 inhibits cell proliferation, migration and invasion in cervical cancer by targeting FGL1


https://www.ribobio.com/en/technical-resources-en/trfs/

tRFs（tRNA-Derived Fragments）是一类长度<40nt，来源于tRNA转录本的非编码小片段RNA。最近的研究表明，这些小RNA片段具有特定的生物学功能，如抑制基因表达、调节细胞凋亡和跨代表观遗传。通过对小RNA片段进行高通量测序和分析可鉴定tRFs，并且它们可对应到已知的tRNA基因。根据tRFs在初级或成熟tRNA转录本上的映射位置，可将其分为不同的类型：

      tiRNA (tiR)或tRNA halves：由应激(和饥饿)诱导的tRNA片段，通过在成熟tRNAs的反密码子环中特异性切割而产生，长约31-40个        碱基，在非应激条件下也可检测到。根据它们是否包含反密码子切割位点的5’或3’序列，tiR有两个亚类：5tiR从成熟tRNA的5’末端          开始并在反密码子环中终止，而3tiR从反密码子环开始并在成熟tRNA的3’末端终止。

      另一类tRNA衍生的小RNAs长约14-30nt，这些小RNAs与miRNAs类似，具有5’磷酸和3’羟基，并且由于其与miRNA的大小相似而          受到广泛关注。根据其映射位置，这些tRFs可分为三种类型：tRF-5、tRF-3和tRF-1。tRF-5和tRF-3分别从成熟tRNAs的5’和3’末端         产生，而tRF-1从初级tRNA转录本的3’末端产生。


https://pubmed.ncbi.nlm.nih.gov/27015120/

http://www.ebiotrade.com/custom/kangchen/160725/index.htm
tRF&tiRNA的定义及分类

       tRFs&tiRNAs是tRNA前体或成熟体在精确调控下由特定的核酸内切酶加工而成的。根据酶切位点的不同，可以分为两类：tiRNAs和tRFs（图2）。

       tiRNAs（亦称tRNA halves）是在多种应激条件下，由特定的核酸内切酶对成熟tRNA的反密码子环进行特异性切割产生，通常由29-50个核苷酸组成，分为5’-halves和3’-halves。
       5’-halves：源自成熟tRNA的5’部分。
       3’-halves：源自成熟tRNA的3’部分。

       tRFs (tRNA-derived fragments)，通常短于tiRNAs，约16-28个核苷酸；根据其在tRNA前体或成熟体上位置的不同，可分为4类：
       tRF-5：源自成熟tRNA的5’端，酶切发生在D-loop区域。
       tRF-3：对应于成熟tRNA的3’部分，包括CCA末端，酶切发生在T-loop区域。
       tRF-1：来源于前体tRNA的3’-UTR区域，3’末端含有多聚U序列。
       i-tRF：区别于以上3种类型，主要来自成熟tRNA的内部区域。
      
tRFdb
http://genome.bioch.virginia.edu/trfdb/
https://pubmed.ncbi.nlm.nih.gov/25270025/

A comprehensive repertoire of tRNA-derived fragments in prostate cancer
https://pubmed.ncbi.nlm.nih.gov/27015120/

GtRNAdb Genomic tRNA Database 
http://gtrnadb.ucsc.edu/

