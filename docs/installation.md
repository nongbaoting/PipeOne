
### 准备软件
1. Java (version >= 1.7)
2. [Nextflow](https://www.nextflow.io/) (version >= 19.10.0)
3. Docker
4. [conda](https://docs.conda.io/en/latest/miniconda.html)
5. [7z](https://www.7-zip.org/download.html)
   `sudo apt install p7zip-full p7zip-rar`
   
### 安装

1.下载PipeOne Docker 镜像
   
```
docker pull nongbaoting/pipeone:latest
```


2.下载PipeOne

```
git clone https://github.com/nongbaoting/PipeOne.git
```


如果只要单独运行随机森林和亚型分析，安装conda 环境即可
```
conda create -n pipeOne_ml python=3.6
source activate pipeOne_ml

cd PipeOne
pip install -r conf/requirement.txt
```



### 下载参考数据
[google drive](https://drive.google.com/drive/folders/1XX9NgpUTRj4llgJq6dGen__-qq4qJ-c0?usp=sharing)

Baidu 网盘:
>链接：https://pan.baidu.com/s/1gbZR1LJAmuT_fmFY1UJ7sA 
>提取码：8fnl

解压，构建索引
```
7z x hg38_ref.7z
cd hg38_ref
## 运行构建索引
main_code_path=/your/path/to/PipeOne/
nextflow run ${main_code_path}/s0_prepare_ref.nf -resume -profile docker --threads 4 
## 删掉中间文件
rm -rf work result
```


### 配置
修改程序配置文件`PipeOne/conf/genomes.config`, 添加索引目录

`ref_directory = ""` 改为 `ref_directory = "/your/path/to/hg38_ref"`

```
params {

  genomes {
   
	"hg38" {
		// 
		ref_directory           = "/your/path/to/hg38_ref"

        // 下面无需改
		species 				= "human"
		genome_build			= "hg38"

		// in the download files
		fasta   				= "${ref_directory}/hg38.fa"
		fasta_fai				= "${ref_directory}/hg38.fa.fai"
		gtf     				= "${ref_directory}/gencode.v32.gtf"
		ref						= "${ref_directory}/hg38_ref.txt"
		genecode_gtf  			= "${ref_directory}/gencode.v32.gtf"
		genecode_lncRNA_gtf		= "${ref_directory}/gencode.v32.long_noncoding_RNAs.gtf"
		
		lncpedia_gtf  			= "${ref_directory}/lncipedia_5_2_hg38.addquotes.gtf"
		repeat_gtf            	= "${ref_directory}/hg38_ucsc.repeat.gtf"
		retro_gtf             	= "${ref_directory}/telescope_HERV_rmsk.hg38.v2/transcripts.gtf"
		sprint_repeat 			= "${ref_directory}/sprint_hg38_repeat.txt"
		blacklisted   			= "${ref_directory}/arriba_db/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"
		cytobands 				= "${ref_directory}/arriba_db/cytobands_hg38_GRCh38_2018-02-23.tsv"
		proteinDomains 			= "${ref_directory}/arriba_db/protein_domains_hg38_GRCh38_2018-03-06.gff3"
		salmon_index_3UTR		= "${ref_directory}/apa_3UTR/gencode_v31/3utr"
		replace_SalmonIndex_ID 	= "${ref_directory}/apa_3UTR/gencode_v31/fa_id_relate.txt"
		apa3utr					= "${ref_directory}/apa_3UTR/gencode_v31/qapa_3utrs.gencode_V31.hg38.bed"
		utr_gtf					=  genecode_gtf
		annovar_data_dir		= "${ref_directory}/annovar_db"


		// generate by prepare_ref.nf
		star_index 	       		= "${ref_directory}/STAR_index"
		bwa_index		        = "${ref_directory}/bwa_index/bwa_index"
		bowtie2_index			= "${ref_directory}/bowtie2_index/bowtie2_base"
		sprint_index 			= "${ref_directory}/sprint_index/genome.fa"
		hisat2_index  			= "${ref_directory}/hisat2_index/genome"

		
		
	}
	  
   
  }

}

```

