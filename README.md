# 
## Install
```
pip install cpat

```

1. conda 安装

2. lncRNA 环境

`conda env create -f lncRNA_env.yml
conda env create  -n -f py3 base.yml
`

3. nextflow 安装


R 包
```
conda install r-ggplot2  r-dplyr r-reshape2 r-VennDiagram r-readr
conda install r-tidyverse r-pheatmap
conda install bioconductor-deseq2  bioconductor-tximport bioconductor-rhdf5
conda install bioconductor-vsn 
conda install bioconductor-clusterProfiler bioconductor-org.mm.eg.db bioconductor-org.hs.eg.db \
bioconductor-Rgraphviz bioconductor-topGO bioconductor-DOSE bioconductor-pathview
```

运行
```
# conda 
nextflow run /public/home/nong/pipeline/lncRNA/lnc_salmon.nf -profile conda -resume --executor local  --scratch false \
  --genome hg19  \
--reads "/public1/nong/projects/npc/cell_line/*_{1,2}.fastq.gz" --cleaned true


# docker
nextflow run /dsk2/who/nbt/pipe/lncRNA/lnc_salmon.nf -profile docker -resume --executor local  --scratch false \
  --genome hg19  \
--reads "/home/nbt2/proj/pipetest/reads/SRP142570/fastq/*.R{1,2}.fastp.fq.gz" --cleaned true


```


## docker 

```
conda install hisat2 stringtie bedops bedtools biopython \
 fastp hisat2 kallisto salmon samtools star taco  cpat  fire pysam
```



## pipeOne

### run all in one

conda
```
bash /dsk2/who/nbt/pipe/lncRNA/pipeOne.sh --reads "../reads/testRaw/*_{1,2}.fq.gz" --genome hg38_124 --cleaned true --profile standard
```

docker