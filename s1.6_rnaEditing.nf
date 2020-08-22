#!/usr/bin/env nextflow

params.reads = ""
params.sra = ""
params.cleaned = false

params.outdir = "results"
params.threads = 16
threads = params.threads
params.saveIntermediateFiles = false

params.sprint_index = params.genome ? params.genomes[ params.genome ].sprint_index ?:false :false
params.sprint_repeat = params.genome ? params.genomes[ params.genome ].sprint_repeat ?:false :false
params.snp = params.genome ? params.genomes[ params.genome ].snp ?:false :false
params.annovar_data_dir = params.genome ? params.genomes[ params.genome ].annovar_data_dir ?:false :false

sprint_repeat = file(params.sprint_repeat)


params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded = false


params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /single/ ){
	ifPaired = false
}

if(params.annovar_data_dir){
	Channel
		.fromPath("${params.annovar_data_dir}/*")
		.ifEmpty{exit 1, "annovar data not found: ${params.annovar_data_dir}"}
		.set{annovar_db }
}
annovar_BinDir = Channel.fromPath("${params.annovar_BinDir}/*.pl")

if(params.reads ){
	
	// defined single or paired reads !
	if(ifPaired ){
		ch = Channel
		.fromFilePairs(params.reads)
		.ifEmpty{error "Can not find any reads matching: ${params.reads}"}
		println "Paire-end mode"
	}else{
		Channel
			.fromPath(params.reads)
			.ifEmpty{error "Can not find any reads matching: ${params.reads}"}
			.map{file ->
				def key = file.name.toString().tokenize('.').get(0)
				return tuple(key,file)
				}
			.groupTuple()
			.set{ch}
		println "sing-end mode"
	}
	ch.into{ch;ch_print;}
	//ch_print.println()
	
	if (! params.cleaned ){
		process fastp{
			tag {id}
            maxForks 2
			
			publishDir "${params.outdir}/fastp/", mode: 'link',
				saveAs: {filename -> 
					if(filename =~ /fastp.fq.gz/) "clean/${filename}"
					else "report/${id}.${filename}" 
					}
			
			input:
			set id, file(reads) from ch
			
			output:
			set id, "${id}*.fastp.fq.gz" into read_ch
			set  "fastp.html", "fastp.json"
			
			script:
			if(! ifPaired ){
				println("fastp single-end!")
				"""
				if [ "${reads}" == "${id}.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 
				"""
			}else
			{
				
				"""
				if [ "${reads[0]}" == "${id}.R1.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean reads name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 
				"""
			}
				
		}
		
		//fastp_out=fastp_out.map{it ->a=[it[0],it[1 ..2 ] ] }
		

	}else{
		read_ch = ch
	}
	
		
}else if(params.sra && !params.reads){
	Channel
	.fromPath(params.sra)
	.ifEmpty{error "Can not find any SRA matching: ${params.sra}"}
	.map{file ->
		def key = file.name.toString().tokenize('.').get(0)
		return tuple(key,file)
		}
	.groupTuple()
	.set{ sra_ch}

	sra_ch.into{sra_ch; print_ch}
	//print_ch.println()
	process fastq_dump_and_fastp{
		tag {id}
		errorStrategy 'ignore'
		maxForks 1
		
		publishDir "${params.outdir}/fastp/", mode: 'link',
			saveAs: {filename -> 
				if(filename =~ /fastp.fq.gz/) "clean/${filename}"
				else "report/${id}.${filename}" 
				}
		
		input:
		set id, file(sra) from sra_ch
		
		output:
		set id, "${id}*fastp.fq.gz" into read_ch
		set "fastp.html", "fastp.json"
		
		script:
		if(!ifPaired ){
			"""
			fasterq-dump --split-3 -o ${id}  $sra
			fastp -i ${id}_1.fastq  -o ${id}.fastp.fq.gz -O -q 20 --thread  2
			"""
		}else{
			"""
			fasterq-dump --split-3 -o ${id}  $sra
			fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  2
			"""
		}
		
	}
	
	//read_ch=fastp_out.map{it ->a=[it[0],it[1 .. 2]] }
    }else{
	exit 1, "you should input --reads or --sra or --bam"
}


read_ch.into{ reads_print; reads_bwa; }
// reads_print.println()

if ( params.sprint_index  ){
   sprint_indices = Channel
        .fromPath("${params.sprint_index}*")
        .ifEmpty { exit 1, "sprint_index not found: ${params.sprint_index}" }
		sprint_base = params.sprint_index.split('/')[-1]
	
}

params.sprint_tmp = "false"

process sprint{
	publishDir "${params.outdir}/sprint/", mode: 'link',
		saveAs: {filename -> 
			if (filename =~ /SPRINT_identified_all.res/) "res/${filename}"
			else "tmp/${filename}"
			
		}
	
	input:
	set id, path(reads) from reads_bwa
	file "sprint_index/*" from sprint_indices.collect()
	path sprint_repeat
	
	output:
	path "${id}" optional true into sprint_out
	file "${id}.SPRINT_identified_all.res" into sprint_res
	
	"""
	set +u; source activate RnaEditing; set -u
	if [[ "${reads[0]}"  =~ gz\$|bzip2\$ ]]; then
		echo "decompressing fastq files ..."
		zcat  ${reads[0]} > read_1.fq
		zcat  ${reads[1]} > read_2.fq
	else
		echo "ln files ..."
		ln -s ${reads[0]}  read_1.fq
		ln -s ${reads[1]}  read_2.fq
	fi
	
	mkdir -m 775 -p tmp
	sprint main -rp ${sprint_repeat} -c 1 -p ${threads} -1 read_1.fq -2 read_2.fq sprint_index/${sprint_base} tmp bwa samtools
	cp tmp/SPRINT_identified_all.res ${id}.SPRINT_identified_all.res
	
	rm read_1.fq read_2.fq
	if [[ "${params.sprint_tmp}" == "true" ]];then
		ln -s tmp ${id}
	else
		rm -rf tmp
	fi
	
	"""
}

params.miniEditing_reads = 5

process merge_A2I {
	
	publishDir "${params.outdir}/sprint/Merge", mode: 'link'
		
	input:
	file "results/*"     from sprint_res.collect()
	file "annovar_db/*"  from annovar_db.collect()
	
	output:
	tuple "SPRINT_all.annovar.csv", "SPRINT_A2I.annovar.csv", "SPrint_A2I_table.annovar.csv" into merge_out
	
	"""
	set +u; source activate pipeOne_ml; set -u
	python3 ${baseDir}/bin/rnaEditing.py cat SPRINT_all.tsv SPRINT_A2I.tsv results
	
	# supported reads for A to I rna-editing is 5
	Rscript ${baseDir}/bin/rnaEditing_rate.R SPRINT_A2I.tsv SPrint_A2I_table.tsv 5

	# annovar annotate rnaEditing
	cut -f1-4 SPRINT_all.tsv |awk 'NR !=1 {split(\$4, a, ""); print \$1,\$3,\$3, a[1],a[2] }' OFS="\t" > rnaEditing_pos.txt
	annotate_variation.pl -geneanno -dbtype refGene -buildver hg38 rnaEditing_pos.txt annovar_db
	python3 ${baseDir}/bin/rnaEditing.py add_tab_anno rnaEditing_pos.txt.variant_function SPrint_A2I_table.tsv SPrint_A2I_table.annovar.csv
	python3 ${baseDir}/bin/rnaEditing.py add_res_anno rnaEditing_pos.txt.variant_function SPRINT_all.tsv SPRINT_all.annovar.csv
	python3 ${baseDir}/bin/rnaEditing.py add_res_anno rnaEditing_pos.txt.variant_function SPRINT_A2I.tsv SPRINT_A2I.annovar.csv
	"""
}


if (params.snp){
	snpFile = file(params.snp)
	
	process snpFilter{
		publishDir "${params.outdir}/sprint_filter/", mode: 'link'
		
		input:
		tuple  "SPRINT_all.annovar.csv", "SPRINT_A2I.annovar.csv", "SPrint_A2I_table.annovar.csv" from merge_out.collect()
		file "snpFile" from snpFile
		
		output:
		set  "SPRINT_all.annovar.filterSNP.csv", "SPRINT_A2I.annovar.filterSNP.csv","SPrint_A2I_table.annovar.filterSNP.tsv"
		
		"""
		python3 ${baseDir}/bin/rnaEditing.py filter_tab snpFile SPrint_A2I_table.annovar.csv SPrint_A2I_table.annovar.filterSNP.tsv
		python3 ${baseDir}/bin/rnaEditing.py filter_res snpFile SPRINT_all.annovar.csv SPRINT_all.annovar.filterSNP.csv
		python3 ${baseDir}/bin/rnaEditing.py filter_res snpFile SPRINT_A2I.annovar.csv SPRINT_A2I.annovar.filterSNP.csv
		"""
	}
}






