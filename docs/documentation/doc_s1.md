
__Note: As long as the corresponding input data is provided, the three main modules can be operate independently__

### __Module 1: RNA-seq processing__

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
--cleaned	<boolean>	true or false. defualt[true]
--layout	<string>	paired or single. defualt [paired]
--library	<string>	polyA or total. defualt [polyA]
--threads	<int>	    number of CPU process for each step. default [8]
--maxForks	<int>	    max forks number of parrallel. default [2]
--profile	<str>	    execution envirenment, e.g. docker, conda. defualt [docker]
--saveIntermediateFiles	save intermediate files defualt [off]
--update_GTF use customized GTF generated in step s1.1_lncRNA.nf instand of GENCODE GTF as input for step: s1.2_circRNA.quant.nf, s1.5_fusion.nf and s1.7_alternative_splicing.nf . defualt [off]
-h --help               print usage
```

#### Output

* __Files result from different programs applied to RNA-seq__
    * 00_tables Tables from different modal of RNA-seq
    * s1.*  Program running directory
      * result result directory
      * work work directory

```
$ tree -L 2
.
├── 00_tables
│   ├── 00_rawdata
│   ├── s1.1_lncR_mRNA
│   ├── s1.2_circRNA
│   ├── s1.3_APA-3TUR
│   ├── s1.4_retrotranscriptome
│   ├── s1.5_fusion
│   ├── s1.6_rnaEditing
│   ├── s1.7_alternative_splicing
│   └── s1.8_SNP
├── one_command.sh
├── s1.1_lncRNA
│   ├── results
│   └── work
├── s1.2_circRNA
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

```

* __Tables from different modal of RNA-seq: s1.*__
```
$ ls  00_tables
00_rawdata
s1.1_lncR_mRNA
s1.2_circRNA
s1.3_APA-3TUR
s1.4_retrotranscriptome
s1.5_fusion
s1.6_rnaEditing
s1.7_alternative_splicing
s1.8_SNP
```

* __Tables of different modals of RNA-seq in one directory__, which are required in module 2 and module3
```
$ tree 00_tables/00_rawdata 
00_tables/00_rawdata
├── APA_pau-distal-proximal.csv
├── circRNA_cpm.csv
├── lncR_gene.tpm.csv
├── merge_graphs_alt_3prime_C3.confirmed.psi.csv
├── merge_graphs_alt_5prime_C3.confirmed.psi.csv
├── merge_graphs_exon_skip_C3.confirmed.psi.csv
├── merge_graphs_intron_retention_C3.confirmed.psi.csv
├── merge_graphs_mult_exon_skip_C3.confirmed.psi.csv
├── merge_graphs_mutex_exons_C3.confirmed.psi.csv
├── prot_gene.tpm.csv
├── retro-FPKM-divide_totalMapReads.csv
├── RNA-editing-rate.csv
└── snp.geneBase.csv

0 directories, 13 files
```