
#### Get test data

 `test_dat.7z`, Download one of the data sets below:

* [google drive](https://drive.google.com/drive/folders/1XX9NgpUTRj4llgJq6dGen__-qq4qJ-c0?usp=sharing)

* Baidu Cloud Disk:

	>Link: [https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA](https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA) 

	>Extraction code: __8fnl__




#### Module 1: RNA-seq processing
```
mkdir s1
cd s1
baseDir=/your/path/to/PipeOne/
nextflow run ${baseDir}/s1_RNAseq.nf  -resume -profile docker \
        --genome hg38 \
        --reads "../test_dat/s1_RNA-seq/*_{1,2}.fq.gz"
```

__Results__:
Tables for different aspects of RNA-seq:
```
$ ls -1 results/tables 
APA_pau-distal-proximal.csv
circRNA_CPM.csv
lncR_gene.tpm.csv
merge_graphs_alt_3prime_C3.confirmed.psi.csv
merge_graphs_alt_5prime_C3.confirmed.psi.csv
merge_graphs_exon_skip_C3.confirmed.psi.csv
merge_graphs_intron_retention_C3.confirmed.psi.csv
merge_graphs_mult_exon_skip_C3.confirmed.psi.csv
merge_graphs_mutex_exons_C3.confirmed.psi.csv
prot_gene.tpm.csv
retro-FPKM-divide_totalMapReads.csv
snp.geneBase.csv
```

#### Module 2: Feature Prioritization
We use these downloaded test tables as input data. 
__Note:__ _In real applications, the result tables of the previous step should be used as input_

```bash
mkdir s2
cd s2
nextflow run ${baseDir}/s2_RF.nf -profile docker \
    --rawdir ../test_dat/s2_tables/00_rawdata  \
    --sample_info ../test_dat/s2_tables/s1_sample_info-tumor-normal.csv \
    --gene_info ../test_dat/s2_tables/protein_coding_and_all_lncRNA.info.tsv

```

__Results__:

results/data/feature_importance.csv

results/data/feature_importance-addName.csv

results/data/discriminative_power_of_topk_feature.csv 

#### Module 3: Subtype Analysis
```bash
mkdir s3
cd s3
nextflow run ${baseDir}/s3_Subtype.nf -profile docker \
	--rawdir   ../test_dat/s3_subtype/00_rawdata/ \
	--clinical ../test_dat/s3_subtype/KIRP_cli.OS.csv \
	--test 
```
__Note:__  _The parameter `--test` is use for test running only, should not be use in real execution._

__Results__:
results/record_log_rank_test_pvalue.csv

results/FeatureSelection

__[Click for detailed documentation](../documentation/doc_s1)__
