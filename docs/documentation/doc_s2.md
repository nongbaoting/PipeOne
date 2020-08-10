
__Module 2: Importance Features__

#### Run in One Command
```
mkdir s2_randomForest
cd s2_randomForest
nextflow run /your/path/to/PipeOne/s2_ml_randomForest.nf  -profile docker --sample_info ../test_dat/s2_tables/s1_sample_info-tumor-normal.csv --rawdir ../test_dat/s2_tables/00_rawdata --var_topK 1000  --threads 8
```
__Options__

```
--sample_info <str> The sample_info file must contains two columns, the column names are Sample and Group
--var_topK <int> top K most variance features. default 1000

```

1.sample_info
   
    The sample_info file must contains two columns, the column names are Sample and Group

    For Example

```
$ cat sample_info.csv
Sample,Group
sample_1,0
sample_2,0
sample_3,1
sample_4,1
sample_n,0
```
Note: There are only two types of group values, use 0 represents normal, 1 represents tumor

2.Directory containing data tables of various types, width is generated in the module 1

    The first row of each table is the sample ID, and the first column is the feature ID.

    For Example:

```
$ cat table.csv
feature_id,sample_1,sample_2,sample3...,sample_n
feature_1,value,value,value,value,...,value
feature_2,value,value,value,value,...,value
feature_3,value,value,value,value,...,value
feature_n,value,value,value,value,...,value
```

Note that each file is in csv comma separated format.




#### Run Module 2  step by step

__1. get top K variance features__
  
```
source activate pipeOne_ml
baseDir="/your/path/to/PipeOne/"
python3 ${baseDir}/bin/ML/proc_raw_data.py proc --rawdir 00_rawdata  --sample_info s1_sample_info-tumor-normal.csv --var_topk 1000 --tdir ./data/proc/
    
```
    --rawdir  RNA-seq result table, such as expression level (TPM value), RNA editing rate, fusion event, etc.
    --sample_info   sample information file
    --var_topk  top K most variance features
    --tdir  Output directory, default ./data/proc/

__2. Divide data into test set and training set__

```
python3 ${baseDir}/bin/ML/proc_raw_data.py train_test_split --indir ../data/proc  --sample_info s1_sample_info-tumor-normal.csv
```
    --indir The output of the previous step
    --sample_info   sample information file
    --test_size The percentage of the random test set. Default 0.3
    --random_state  Random seed, to get reproduce result. default 2

Output directory:
>train_dir Training data directory

>test_dir  Testing data directory

tables of the directory:
```
$ ls train_dir
dat_sample_clustering.csv
proc_APA_pau-distal-proximal.csv
proc_fusion_arriba_out.csv
proc_lncR_gene.tpm.csv
proc_merge_graphs_alt_3prime_C3.confirmed.psi.csv
proc_merge_graphs_alt_5prime_C3.confirmed.psi.csv
proc_merge_graphs_exon_skip_C3.confirmed.psi.csv
proc_merge_graphs_intron_retention_C3.confirmed.psi.csv
proc_merge_graphs_mult_exon_skip_C3.confirmed.psi.csv
proc_merge_graphs_mutex_exons_C3.confirmed.psi.csv
proc_prot_gene.tpm.csv
proc_retro-FPKM-divide_totalMapReads.csv
proc_RNA-editing-rate.csv
proc_snp.geneBase.csv
```

Or use a set of fixed sample

```
python3  ${baseDir}/bin/ML/python_code_2/subset_sample_nong.py subset data/proc ./data/train_dir train.csv
python3 /home/nbt2/pipe/PipeOne/bin/ML/python_code_2/subset_sample_nong.py subset --indir ./data/proc --sample_info test.csv --tdir ./data/test_dir 
```
    --indir The output of the previous step
    --sample_info   sample information file
    --tdir Output directory



__3. Run the main program__

```
python3 ${baseDir}/bin/ML/main_randomForest.py --threads 8 --train_dir ./data/train_dir --test_dir ./data/test_dir
```
    --threads   number of threads to use
    --train_dir Training data directory
    --test_dir  Testing data directory


__4. add gene info to  results (optional)__
```
python3 ${baseDir}/bin/ML/result_summary.py feature  --rf_res_fi data/feature_importance.csv  --ginfo_fi ../../s1.1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv
```
    --rf_res_fi Results of Random Foreast (data/feature_importance.csv)
    --ginfo_fi  The gene information file can be found in the results of the first step (s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv)

#### __Output files__:

* results/data/feature_importance.csv
    > random forest 的feature importance

* results/data/feature_importance-addName.csv

* results/data/discriminative_power_of_topk_feature.csv

    >module evaluation
    >
    >senitivity  = tp / (tp + fn)

    >specificity = tn / (tn + fp)

    >tp: true positive

    >fn: false negative

    >tn: true negative
    
    >fp: false positive
