__主要的三部分模块都是可以独立运行的，只要提供相应的输入数据__

### 1. RNA-seq数据处理
####  使用
```
mkdir s1_pipeone_raw
cd s1_pipeone_raw
bash /your/path/to/PipeOne/s1_pipeOne.sh --reads "../test_dat/s1_RNA-seq/*.R{1,2}.fastp.fq.gz" --genome hg38 --cleaned true
```
#### 参数
Require:
```
--reads  <string>     input Fastq files, for example: "/home/reads/*_{1,2}.fastq.gz"
--genome <string>     the genome profile your set in config files: conf/igenomes.config
```
Optional:
```
--cleaned  <booloen>   true or false
--layout   <string>     paired or single. defualt [paired]
--threads  <int>       number of CPU process for each steps
--maxForks <int>      max forks number of parrallel
--profile  <str>       execution envirenment. defualt [docker]
--saveIntermediateFiles       save intermediate files defualt [off]
--saveIntermediateHisat2Bam   save intermediate hisat2 BAM files. defualt [off]
-h --help     print usage
```


#### 输出
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


### 2. 寻找重要特征
需要
* sample_info, csv 格式
    >sample_info 必须有两列信息, 列名为 Sample 和Group

* 包含个类型数据的文件夹, 每种数据都是以csv格式

#### 分步骤运行



* 取 top N variance features
```
source activate pipeone
python3 /home/nbt2/proj/2020-pipeOne-test/python_code_2/proc_raw_data.py proc ../../00_rawdata/ ../s1_sample_info-tumor-normal.csv 3000
```
  
* 讲数据分为测试集和训练集
```
python3 /home/nbt2/proj/2020-pipeOne-test/python_code_2/proc_raw_data.py train_test_split data/proc ../s1_sample_info-tumor-normal.csv
```

* 运行主程序
```
python ./main.py
```

* 整理结果
```
python3 /home/nbt2/proj/2020-pipeOne-test/python_code_2/result_summary.py feature data/feature_importance.csv
```


## 3. 亚型分析