
__3. 亚型分析__

#### 运行

```
mkdir s3_subtype
cd s3_subtype
nextflow run /your/path/to/PipeOne/s3_subtype.nf -resume -profile docker --rawdir ../test_dat/s3_subtype/00_rawdata --clinical ../test_dat/s3_subtype/KIRP_cli.OS.csv --var_topK 100 
```

    --rawdir    RNA-seq各个结果的表格，例如表达水平（TPM值），RNA-editing rate, Fusion event 等
    --clinical  样品生存数据

--clinical 样品生存信息
clinical 文件必须有三列，列名分别为 Sample, Event, Time_to_event
```
$ head KIRP_cli.OS.csv
Run,Event,Time_to_event
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

#### 结果

* `record_log_rank_test_pvalue.csv`,  不同分类结果的生存差异