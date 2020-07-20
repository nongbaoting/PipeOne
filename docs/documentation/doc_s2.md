
__2. 寻找重要特征__

#### 需要
1.sample_info
   
    sample_info 必须有两列信息, 列名为 Sample 和Group

    例如：

```
$ cat sample_info.csv
Sample,Group
sample_1,0
sample_2,0
sample_3,1
sample_4,1
sample_n,0
```
注意：group的值只有两组

2.包含个类型数据的文件夹

    每个table的第一行是sample ID，第一列是feature。

    例如：

```
$ cat table.csv
feature_id,sample_1,sample_2,sample3...,sample_n
feature_1,value,value,value,value,...,value
feature_2,value,value,value,value,...,value
feature_3,value,value,value,value,...,value
feature_n,value,value,value,value,...,value
```

注意每个文件都是csv逗号分割格式。

#### 一步运行
```
mkdir s2_randomForest
cd s2_randomForest
nextflow run /your/path/to/PipeOne/s2_ml_randomForest.nf  -profile docker --sample_info ../test_dat/s2_tables/s1_sample_info-tumor-normal.csv --rawdir ../test_dat/s2_tables/00_rawdata --var_topK 1000  --threads 8
```



#### 分步骤运行
* 取 top N variance features
  
```
source activate pipeOne_ml
baseDir="/your/path/to/PipeOne/"
python3 ${baseDir}/bin/ML/proc_raw_data.py proc --rawdir 00_rawdata  --sample_info s1_sample_info-tumor-normal.csv --var_topk 1000 --tdir ./data/proc/
    
```
    --rawdir    RNA-seq各个结果的表格，例如表达水平（TPM值），RNA-editing rate, Fusion event 等
    --sample_info   样品信息
    --var_topk  各表格用多少feature 数，默认1000
    --tdir  输出目录， 默认 ./data/proc/

* 讲数据分为测试集和训练集
```
python3 ${baseDir}/bin/ML/proc_raw_data.py train_test_split --indir ../data/proc  --sample_info s1_sample_info-tumor-normal.csv
```
    --indir 上一步的输出结果
    --sample_info   样品信息
    --test_size 随机测试集百分比。 默认 0.3
    --random_state  随机种子，为得到重复的测试集

输出目录：
>train_dir 训练数据目录
>test_dir  测试数据目录

目录的内容:

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


或者用固定的sample
```
python3  ${baseDir}/bin/ML/python_code_2/subset_sample_nong.py subset data/proc ./data/train_dir train.csv
python3 /home/nbt2/pipe/PipeOne/bin/ML/python_code_2/subset_sample_nong.py subset --indir ./data/proc --sample_info test.csv --tdir ./data/test_dir 
```
    --indir 上一步的输出结果
    --sample_info   样品信息
    --tdir 输出目录



* 运行主程序

```
python3 ${baseDir}/bin/ML/main_randomForest.py --threads 8 --train_dir ./data/train_dir --test_dir ./data/test_dir
```
    --threads   用多少线程
    --train_dir 训练数据目录
    --test_dir 测试数据目录


* 整理结果（可选）
```
python3 ${baseDir}/bin/ML/result_summary.py feature  --rf_res_fi data/feature_importance.csv  --ginfo_fi ../../s1.1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv
```
    --rf_res_fi Random Foreast 的结果（data/feature_importance.csv）
    --ginfo_fi  基因的信息文件，可在第一步的结果里面找到（s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv）

#### __结果文件__:

* results/data/feature_importance.csv
    > random forest 的feature importance
* results/data/feature_importance-addName.csv

* results/data/discriminative_power_of_topk_feature.csv
    > 模型区分的两种的能力
    >senitivity  = tp / (tp + fn)
    >specificity = tn / (tn + fp)
    >tp: true positive
    >fn: false negative
    >tn: true negative
    >fp: false positive
