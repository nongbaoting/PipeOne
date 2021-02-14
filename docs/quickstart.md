
#### Get test data

 `test_dat.7z`, Download one of the data sets below:

* [google drive](https://drive.google.com/drive/folders/1XX9NgpUTRj4llgJq6dGen__-qq4qJ-c0?usp=sharing)

* Baidu Cloud Disk:

	>Link: [https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA](https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA) 

	>Extraction code: __8fnl__


* Decompress the test data
`
7z x  test_dat.7z
`

#### Module 1: RNA-seq processing
```
mkdir s1_pipeone_raw
cd s1_pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz" --genome hg38 --cleaned true
```

Tables for different aspects of RNA-seq:

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

#### Module 2: Feature Prioritization
We use thess test tables as input data. in actual applications, one should use the result table of the previous step as input 

```bash
mkdir s2_randomForest
cd s2_randomForest
nextflow run /home/nbt2/pipe/dev/pipeOne-v2/s2_RF.nf  -profile docker \
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
mkdir s3_subtype
cd s3_subtype
nextflow run ~/pipe/dev/pipeOne-v2/s3_Subtype.nf -profile docker \
	--rawdir ../test_dat/s3_subtype/00_rawdata/ \
	--clinical ../test_dat/s3_subtype/KIRP_cli.OS.csv 

```
Note: We select the top 100 features with the largest variance in each table to facilitate the test program; in real operation, it should be a larger number, such as 1000

__[For detail](../documentation/doc_s1)__
