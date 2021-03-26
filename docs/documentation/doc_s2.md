
### __Module 2: Feature Prioritization__

#### Run in One Command
```bash
mkdir s2
cd s2
baseDir=/your/path/to/PipeOne/
nextflow run ${baseDir}/s2_RF.nf -profile docker \
    --rawdir ../test_dat/s2_tables/00_rawdata  \
    --sample_info ../test_dat/s2_tables/s1_sample_info-tumor-normal.csv \
    --gene_info ../test_dat/s2_tables/protein_coding_and_all_lncRNA.info.tsv
```

#### Options
Require:
```bash
--rawdir <string> directory contain RNA-seq result tables produced in module 1
--sample_info <string> The sample_info file must contains two columns, the column names are Sample and Group
```

Optional:
```bash
--var_topK <int> top K most variance features for each table. default [1000]
--gene_info <string> The gene information file can be found in the results of the first step (s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv)
--test_size The percentage of the random test set. Default 0.25
--random_state  Random seed, to get reproduce result. default 2
```

* `--rawdir`, Directory containing data tables of various types, which is generated in the module 1. The first row of each table is the sample ID, and the first column is the feature ID.

For Example:
```bash
$ cat table.csv
feature_id,sample_1,sample_2,sample3...,sample_n
feature_1,value,value,value,value,...,value
feature_2,value,value,value,value,...,value
feature_3,value,value,value,value,...,value
...
feature_n,value,value,value,value,...,value
```
Note: _the format of each file in the directory is separated by  commas ( i.e. csv format )_

* `--sample_info`, The sample_info file must contains two columns, the column names are Sample and Group

For Example
```bash
$ cat sample_info.csv
Sample,Group
sample_1,0
sample_2,0
sample_3,1
sample_4,1
...
sample_n,0
```
Note: _there are only two types of group values, use 0 represents normal, 1 represents tumor_

* `--gene_info` information file at least contain two columns `gene_name`, `gene_id`. This parameter will convert gene_id to gene_name in the result.

#### __Output files__:

* `results/data/feature_importance*.csv`, random forest feature importance (feature weight)

* `results/data/feature_importance-addName.csv`, feature importance with gene name ( when `--gene_info` is provided ).

* `results/data/discriminative_power_of_topk_feature.csv`, module evaluation.
    * senitivity  = tp / (tp + fn)
    * specificity = tn / (tn + fp)
    * tp: true positive
    * fn: false negative
    * tn: true negative
    * fp: false positive
