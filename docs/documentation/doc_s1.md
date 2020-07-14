__主要的三部分模块都是可以独立运行的，只要提供相应的输入数据__

__1. RNA-seq数据处理__

####  使用
```
mkdir s1_pipeone_raw
cd s1_pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz" --genome hg38 --cleaned true
```
#### 参数
Require:
```
--reads  <string>     输入FASTQ文件, for example: "/home/reads/*_R{1,2}.fastq.gz"，
                      * 代表可以匹配的文件， {1,2} 代表双端测序

--genome <string>     使用的在 conf/igenomes.config 定义genome 版本
```
Optional:
```
--cleaned  <booloen>    true 或者 false FASTQ文件是否是质量控制过，如果是false 则用fastp进行质量控制。defalt [false]
--layout   <string>     paired 或者 single. 默认 [paired]
--threads  <int>       每一步使用的最大线程数。 默认 [8]
--maxForks <int>       并行多少步骤。 默认 [2]
--saveIntermediateFiles       是否保留中间文件到结果。 默认 [off]
--saveIntermediateHisat2Bam   是否保留hisat2 的bam 文件. 默认 [off]
-h --help     显示帮助
```


#### 输出

* __RNA-seq 经过不同的程序运行的结果__
    * 00_tables 经过不同的程序运行得到的表格
    * s1.*文件夹, 经过不同的程序运行目录
      * result 结果目录
      * work 工作目录
```
$ tree -L 2
.
├── 00_tables
│   ├── 00_rawdata
│   ├── s1.1_lncR_mRNA
│   ├── s1.3_APA-3TUR
│   ├── s1.4_retrotranscriptome
│   ├── s1.5_fusion
│   ├── s1.6_rnaEditing
│   ├── s1.7_alternative_splicing
│   └── s1.8_SNP
├── s1.1_lncRNA
│   ├── results
│   └── work
├── s1.3_APA-3TUR
│   ├── results
│   └── work
├── s1.4_retrotranscriptome
│   ├── results
│   └── work
├── s1.5_fusion
│   ├── results
│   └── work
├── s1.6_rnaEditing
│   ├── results
│   └── work
├── s1.7_alternative_splicing
│   ├── results
│   └── work
└── s1.8_SNP
    ├── results
    └── work

30 directories, 0 files
```

* __RNA-seq 经过不同的程序运行得到的表格: s1.*文件夹__
```
$ ls  00_tables
00_rawdata
s1.1_lncR_mRNA
s1.3_APA-3TUR
s1.4_retrotranscriptome
s1.5_fusion
s1.6_rnaEditing
s1.7_alternative_splicing
s1.8_SNP
```


* __RNA-seq 数据处理后得到的各个表格汇总，为下一步运行所需__
```
$ ls 00_tables/00_rawdata
APA_pau-distal-proximal.csv
lncR_gene.tpm.csv
merge_graphs_alt_3prime_C3.confirmed.psi.csv
merge_graphs_alt_5prime_C3.confirmed.psi.csv
merge_graphs_exon_skip_C3.confirmed.psi.csv
merge_graphs_intron_retention_C3.confirmed.psi.csv
merge_graphs_mult_exon_skip_C3.confirmed.psi.csv
merge_graphs_mutex_exons_C3.confirmed.psi.csv
prot_gene.tpm.csv
retro-FPKM-divide_totalMapReads.csv
RNA-editing-rate.csv
```

