

#### Get test data
* [Google Drive](https://drive.google.com/drive/folders/1XX9NgpUTRj4llgJq6dGen__-qq4qJ-c0?usp=sharing): test_dat.7z
* OR Baidu
  
    链接：https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA 

    提取码：8fnl


* 解压测试数据
`
7z x  test_dat.7z
`

#### 1. RNA-seq数据处理
```
mkdir s1_pipeone_raw
cd s1_pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz" --genome hg38 --cleaned true
```

RNA-seq不同方面的表格：
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

#### 2. 寻找重要特征
我们使用示例表格作为输入数据，实际应用中使用上一步骤的结果表格作为输入

输入文件
1. 表格目录
2. 样品信息，至少有两列：Sample,Group. Group 用数字标明组号，比如normal样品标0， 而tumor样品标1
```
mkdir s2_randomForest
cd s2_randomForest
nextflow run /your/path/to/PipeOne/s2_ml_randomForest.nf  -profile docker --sample_info ../test_dat/s2_tables/s1_sample_info-tumor-normal.csv --rawdir ../test_dat/s2_tables/00_rawdata --var_topK 1000  --threads 8
```

__结果文件__:

results/data/feature_importance.csv

results/data/feature_importance-addName.csv

results/data/discriminative_power_of_topk_feature.csv 

#### 3. 癌症亚型分析
```
mkdir s3_subtype
cd s3_subtype
nextflow run /your/path/to/PipeOne/s3_subtype.nf -resume -profile docker --rawdir ../test_dat/s3_subtype/00_rawdata --clinical ../test_dat/s3_subtype/KIRP_cli.OS.csv --var_topK 100 
```
注意： 我们选取各表格方差最大的前100个特征，方便测试程序；在真正运行是应该取更大的数，比如1000

__[详细使用说明](../../documentation/documentation)__

