# PipeOne: an integrative RNA-seq analysis pipeline

### 安装前准备
1. Java (version >= 1.7)
2. [Nextflow](https://www.nextflow.io/) (version >= 19.10.0)
3. 容器 Docker or [Singularity](https://singularity.lbl.gov/)
   
### 安装
1. 下载PipeOne Docker 镜像
```
docker pull nongbaoting/pipeone:latest
```
* 或者 PipeOne Singularity 镜像
```
singularity pull --name nf-core-lncpipe.simg shub://nf-core/lncpipe
```

2. 下载软件包
```
git clone https://github.com/nongbaoting/PipeOne.git
```


### 软件参考基因组，注释文件

可直接下载打包好的全部参考数据


### 快速开始

#### 1. RNA-seq数据处理
```
bash /dsk2/who/nbt/pipe/lncRNA/pipeOne.sh --reads "../reads/testRaw/*_{1,2}.fq.gz" --genome hg38_124 --cleaned true --profile docker
```

#### 2. 寻找重要特征
__2.1__ 取各类型的数据的前N(默认1000)方差最大的特征
```
python 
```

### 详细说明文档()

### 引文
PipeOne 的论文请查阅文献: 