
2. 寻找重要特征
===================

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
