
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


* 运行主程序
```
python3 ${baseDir}/bin/ML/main_randomForest.py --threads 8
```
    --threads   用多少线程


* 整理结果（可选）
```
python3 ${baseDir}/bin/ML/result_summary.py feature  --rf_res_fi data/feature_importance.csv  --ginfo_fi ../../s1.1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv
```
    --rf_res_fi Random Foreast 的结果（data/feature_importance.csv）
    --ginfo_fi  基因的信息文件，可在第一步的结果里面找到（s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv）

#### __结果文件__:

* results/data/feature_importance.csv

* results/data/feature_importance-addName.csv

* results/data/discriminative_power_of_topk_feature.csv
