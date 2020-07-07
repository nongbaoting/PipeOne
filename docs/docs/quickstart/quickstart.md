

### 快速开始

#### 1. RNA-seq数据处理
```
mkdir pipeone_raw
cd pipeone_raw
bash /dsk2/who/nbt/pipe/lncRNA/pipeOne.sh --reads "../reads/testRaw/*_{1,2}.fq.gz" --genome hg38_124 --cleaned true --profile docker
```

#### 2. 寻找重要特征

```
mkdir pipeone_mlRF
cd pipeone_mlRF
nextflow run PipeOne/ml_randomForest.nf  --sample_info ../s1_sample_info-tumor-normal.csv --tables ../../00_rawdata --threads 24
```

结果文件:
results/data/feature_importance.csv
results/data/feature_importance-addName.csv
results/data/discriminative_power_of_topk_feature.csv

#### 3. 癌症亚型分析

```
mkdir pipeone_subtype
cd pipeone_subtype
nextflow run 

```

### [详细说明文档](../documentation)

