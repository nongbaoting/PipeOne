
### 准备软件
1. Java (version >= 1.7)
2. [Nextflow](https://www.nextflow.io/) (version >= 19.10.0)
3. 容器 Docker or [Singularity](https://singularity.lbl.gov/)
4. conda
   
### 安装

1.下载PipeOne Docker 镜像
   
```
docker pull nongbaoting/pipeone:latest
```


2.下载PipeOne

```
git clone https://github.com/nongbaoting/PipeOne.git
```

至此就安装完成了

__可选__

如果要单独运行随机森林和亚型分析，安装conda 环境即可
```
conda create -n pipeone python=3.6
source activate pipeone

cd PipeOne
pip install -r conf/requirement.txt
```



### 下载参考数据
[google drive](https://drive.google.com/file/d/1Z7E9DkNdG-zDyxhXx4jCUBFIzLNcHAag/view?usp=sharing)
[baidu 网盘]()
解压，构建索引
```
7z e hg38_ref.7z
cd hg38_ref
nextflow run /dsk2/who/nbt/pipe/PipeOne/prepare_ref.nf -resume 
```


### 配置
修改`PipeOne/conf/genomes.config`, 添加索引目录

`ref_directory = ""` 改为 `ref_directory = "/you/path/to/pipeOne_ref"`

