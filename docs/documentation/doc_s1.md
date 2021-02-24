
__Note: As long as the corresponding input data is provided, these three main modules can be operate independently__

### __Module 1: RNA-seq processing__

####  Usage
```
mkdir s1
cd s1
baseDir=/your/path/to/PipeOne/
nextflow run ${baseDir}/s1_RNAseq.nf \
	-profile docker \
	--genome hg38 \
	--reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz"
```

#### Options
Require:
```
--reads  <string>   The FASTQ gzip files, for example: "/home/reads/*_R{1,2}.fastq.gz"

--genome <string>   Genome version defined in conf/igenomes.config
```

Optional:
```
--run_s1 choose programs to proceed use comma to seperate different programs. defualt ['1,2,3,4,5,6,7,8']
        1 represent task 'mRNA_lncRNA',
        2 represent task 'circRNA',
        3 represent task 'APA',
        4 represent task 'RetroTrans',
        5 represent task 'Fusion',
        6 represent task 'RNAediting',
        7 represent task 'AS',
        8 represent task 'SNP'.
        user could use the task number or task name, such as '1,2,RetroTrans,Fusion,RNAediting,AS,8'.
--cleaned   true or false. defualt[true]
--singleEnd paired or single. defualt [paired]
--library   <string>	polyA or total. defualt [polyA]
--max_cpus  <int>	    number of CPU process for each step. default [24]
--max_memory    <string>   max max_memory for processing. default [ 128.GB]
--max_time  <string>   number of hours for program to proceed. default [2400.h]
--maxForks  <int>	    max forks number of parrallel. default [2]
--saveIntermediateFiles save intermediate files defualt [off]
--update_GTF    use customized GTF generated in step s1.1_lncRNA.nf instand of GENCODE GTF as input for step: s1.2_circRNA.quant.nf, s1.5_fusion.nf and s1.7_alternative_splicing.nf . defualt [off]
-h --help               print usage
```

#### Output

* __Files result from different programs applied to RNA-seq__
    * __s1.*__, results of different submodules
    * __tables__, all tables from different modal of RNA-seq, which are required in module 2 and module3

```bash
-> % tree -L 2 results
.
├── pipeline_info
│   ├── pipeline_report.html
│   ├── pipeline_timeline.html
│   └── pipeline_trace.txt
├── s1.0_QC
│   └── fastp
├── s1.1_mRNA_lncRNA
│   ├── annotations_and_fasta
│   ├── coding_potential
│   ├── gffcompare
│   ├── hisat2
│   ├── novel_lncRNA
│   ├── reference_gtf_info
│   ├── salmon
│   ├── stringtie
│   └── taco
├── s1.2_circRNA
│   └── CIRIquant
├── s1.3_APA
│   ├── apa_3utr
│   ├── pau_results.filterPau-distal-proximal.txt
│   ├── pau_results.filterPau.txt
│   └── salmon
├── s1.4_RetroTrans
│   ├── bowtie2
│   └── telescope
├── s1.5_Fusion
│   └── arriba
├── s1.6_RNAediting
│   └── sprint
├── s1.7_AS
│   ├── spladder_out
│   ├── spladder_out_table
│   └── STAR
├── s1.8_SNP
│   ├── annovar_table
│   ├── variantAnnotateAnnovar
│   └── Variant_filtering
└── tables
    ├── APA_pau-distal-proximal.csv
    ├── circRNA_CPM.csv
    ├── lncR_gene.tpm.csv
    ├── merge_graphs_alt_3prime_C3.confirmed.psi.csv
    ├── merge_graphs_alt_5prime_C3.confirmed.psi.csv
    ├── merge_graphs_exon_skip_C3.confirmed.psi.csv
    ├── merge_graphs_intron_retention_C3.confirmed.psi.csv
    ├── merge_graphs_mult_exon_skip_C3.confirmed.psi.csv
    ├── merge_graphs_mutex_exons_C3.confirmed.psi.csv
    ├── prot_gene.tpm.csv
    ├── retro-FPKM-divide_totalMapReads.csv
    └── snp.geneBase.csv
```
