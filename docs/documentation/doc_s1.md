__As long as the corresponding input data is provided, the three main modules can operate independently__

__Module 1: RNA-seq processing__

####  Usage
```
mkdir s1_pipeone_raw
cd s1_pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz" --genome hg38 --cleaned true
```
#### Options
Require:
```
--reads  <string>     FASTQ gzip files, for example: "/home/reads/*_R{1,2}.fastq.gz"

--genome <string>     Genome version defined in conf/igenomes.config
```
Optional:
```
--cleaned  <booloen>   true or false. Whether the FASTQ file has been quality controlled, if it is false, use fastp for quality control. defalt [false]
--layout   <string>     paired or single. Default [paired]
--threads  <int>       The maximum number of threads used in each step. Default [8]
--maxForks <int>       How many steps in parallel. Default [2]
--saveIntermediateFiles       Whether to keep the intermediate file to the result. Default [off]
--saveIntermediateHisat2Bam   Whether to keep the bam file of hisat2. Default [off]
-h --help     show  help message
```


#### Output

* __Files result from different programs applied to RNA-seq__
    * 00_tables Tables from different aspect of RNA-seq
    * s1.*  Program running directory
      * result result directory
      * work work directory
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

* __Tables from different aspect of RNA-seq: s1.*
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


* __Tables of different aspects of RNA-seq from one directory, which are required in the next module__
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