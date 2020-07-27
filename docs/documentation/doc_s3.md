
__Module 3: Subtype Analysis__

#### Run in one command

```
mkdir s3_subtype
cd s3_subtype
nextflow run /your/path/to/PipeOne/s3_subtype.nf -resume -profile docker --rawdir ../test_dat/s3_subtype/00_rawdata --clinical ../test_dat/s3_subtype/KIRP_cli.OS.csv --var_topK 1000 
```
Required: 

    --rawdir <str>  directory contain RNA-seq result table, such as expression level (TPM value), RNA editing rate, fusion event, etc.
    --clinical <str> Clinical information file

Optional:
    
    --var_topK <int>  top K most variance features. default [1000]
    -- cluster_range <str> cluster of to test in Kmeans default ["3-8"]
    --threads <int> number of threads to use default [24]
    

--clinical Clinical information
The clinical  file must have three columns, the column names are Sample, Event, Time_to_event
```
$ head KIRP_cli.OS.csv
Sample,Event,Time_to_event
105248a5-eb2a-4336-b76c-cb759f45e45b_gdc_realn_rehead,0,214
e263e4a0-1489-40a9-8e89-9ee6aa19cdc7_gdc_realn_rehead,0,2298
be3d303c-ba1a-4cf8-a8a5-b5adcae05d14_gdc_realn_rehead,0,1795
dc5d11b5-742f-4740-acd0-2806962d9d1b_gdc_realn_rehead,1,1771
291bab1e-5d83-4b8f-9212-ecec1278ea1a_gdc_realn_rehead,0,3050
5779df14-6e54-4f72-97d3-5486e31ffc5d_gdc_realn_rehead,0,1731
a58cd3b2-3b69-4677-87f3-82b256bbcc48_gdc_realn_rehead,1,139
753f54d4-1994-4a3a-adeb-a37df28973d6_gdc_realn_rehead,0,2790
ff1a9e27-04bb-4732-8631-4ffad26c700e_gdc_realn_rehead,0,2294
```


#### Run step by step

1. select topK variance features and run NMF
```
source activate pipeOne_ml
baseDir=/path/to/PipeOne/

python3 ${baseDir}/bin/ML/python_code_2/proc_raw_data.py proc --rawdir 00_rawdata/ --sample_info sample.cli.csv --var_topk 1000

## defusion
python3 ${baseDir}/bin/ML/python_code_2/run_defusion.py --threads 24
python3 ${baseDir}/bin/ML/python_code_2/check_convergence.py 
```

__Options__

    --rawdir <str>  directory contain RNA-seq result table, such as expression level (TPM value), RNA editing rate, fusion event, etc.
    --clinical <str> Clinical information file
    --var_topK  top K most variance features. default [1000]
    --threads <int> number of threads to use default [24]

2. clustering and eval
```
python3 ${baseDir}/bin/ML/python_code_2/eval_cluster_num.py --cluster_range "3-8"
Rscript ${baseDir}/bin/ML/python_code_2/survival_eval.R ./data/sample.cli.csv ./clusters/surv_curve/ "3-8"
```

__Options__

    --cluster_range <str> cluster of to test in Kmeans default['3-8']

Script `survival_eval.R` need three mandotory inputs: __clinical information file__,  
__cluster_result__ which produce by script `eval_cluster_num.py`, and __cluster_range__ same as --cluster_range


3. select features
```
python3 ${baseDir}/bin/ML/python_code_2/select_topk_nong.py 
python3 ${baseDir}/bin/ML/python_code_2/find_best_RFparams.py 
```

__Options__

`select_topk_nong.py`

    --topK_importance <str>  top most importance features use to retraining. default ['50,100,200']
    --outdir <str> default [./data_randomForest]
    --cluster_survival_file <str>  default [record_log_rank_test_pvalue.csv].  when this option is provide program will select the clustering result with max silhoutte with value and log rand test p value is significance
    --cluster_file <str>  specify a cluster result file under ./clusters/surv_curve/ . default [None]. when this option is provide program will use the cluster

__Note:__ The `--cluster_survival_file` option is mutually exclusive with the `--cluster_file` option.

`find_best_RFparams.py`
    --ddir  input direcotry generate by last step. default ./data_randomForest
    --tdir output directory default ./FeatureSelection/



#### Results

__main result__
* `record_log_rank_test_pvalue.csv`, Survival difference of different classification results
* `FeatureSelection` random Forest results
    * `RF_best_params_settings_for_feature_selection.csv`   best Random Foerest parameter and accuracy
    * `feature*_importance.csv`  Random Foerest Feature importance
    * `RF_params_setting_record.txt` all Random Foerest parameter and accuracy

__Intermediate files__
* `NMF` results of Nonnegative matrix factorization (NMF) under different paramiters
    `weight_*` weight matrix W of NMF
    `X_*`      H matrix

* `clusters` cluster reuslt and survival curves
    * `clusters/eval_cluster_num/lowDim=*_alpha=*_gamma=*_clusters=*_clustering.csv` Kmeans cluster base of H matrix
    * `clusters/lowDim=*_alpha=*_gamma=*_silhouette_score.png` silhouette with of cluster results
    * `clusters/surv_curve/low_dim=2_alpha=0.01_gamma=0.00_clustering.pdf` survival plot of Kmeans cluster result which silhouette width is max

* `data_randomForest` selected data for running random Forest
    `top*` top K features in W matrix use for random forest training

* `rf_my_record.csv`