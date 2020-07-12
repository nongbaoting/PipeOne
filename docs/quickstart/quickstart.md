


#### 1. RNA-seq数据处理
```
mkdir pipeone_raw
cd pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../reads/testRaw/*_{1,2}.fq.gz" --genome hg38 --cleaned true --profile docker
```

#### 2. 寻找重要特征

```
mkdir pipeone_mlRF
cd pipeone_mlRF
nextflow run PipeOne/ml_randomForest.nf  --sample_info ../s1_sample_info-tumor-normal.csv --tables ../../00_rawdata --threads 24
```

__结果文件__:

results/data/feature_importance.csv

results/data/feature_importance-addName.csv

results/data/discriminative_power_of_topk_feature.csv 

#### 3. 癌症亚型分析


```
mkdir pipeone_subtype
cd pipeone_subtype
nextflow run PipeOne/s3_subtype.nf -resume --tables ../../00_rawdata --clinical KIRP_cli.OS.csv

```

__[详细使用说明](../../documentation/documentation)__

