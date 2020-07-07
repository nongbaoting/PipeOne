
### 需要软件
1. Java (version >= 1.7)
2. [Nextflow](https://www.nextflow.io/) (version >= 19.10.0)
3. 容器 Docker or [Singularity](https://singularity.lbl.gov/)
4. conda
   
### 安装
1. 下载PipeOne Docker 镜像
```
docker pull nongbaoting/pipeone:latest
```


下载软件包, 并安装conda 环境（为找重要feature)
```
git clone https://github.com/nongbaoting/PipeOne.git
cd PipeOne
conda create -n pipeone python=3.6
source activate pipeone
pip install -r conf/requirement.txt
```



### 下载参考基因组，注释文件等

解压，构建索引
```
cd 
nextflow run 
```


### 配置


