#!/usr/bin/env nextflow
nf_required_version = '0.31.1'
params.reads = ""
// "/public1/pub/guohh/other/mouse_fatty/00_rawdata/*.fastq.gz"
params.sra   = ""
params.hisat2_bam = ""
params.genome ="hg19"
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?:false :false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false
params.threads = 12

threads = params.threads
bwa_index = params.bwa_index
scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2}

params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /single/ ){
	ifPaired = false
}

if ( params.fasta ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}


if ( params.gtf ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}

if ( params.bwa_index  ){
    bwa_indices = Channel
        .fromPath("${params.bwa_index}.*")
        .ifEmpty { exit 1, "bwa index not found: ${params.bwa_index}" }
	bwa_base = params.bwa_index.split('/')[-1]
}

if ( params.hisat2_index ){
	Channel
		.fromPath("${params.hisat2_index}*")
		.ifEmpty { exit 1, "hisat2_index not found: ${params.hisat2_index}" }
		.set{hisat2_indices }
	hisat2_indices.into{hisat2_indices; hisat2_indices_2 }

	hisat2_base= params.hisat2_index.split('/')[-1]

}

// input reads or sra !
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
	ch_print.println()
	
	if (! params.cleaned ){
		process fastp{
			tag {id}
            maxForks 2
			
			publishDir "${params.outdir}/fastp/", mode: 'copy',
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
				println("fastp paired!")
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
		
		publishDir "${params.outdir}/fastp/", mode: 'copy',
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
	exit 1, "you should input --reads or --sra"
}


process cutadapt_fixLength {
	maxForks 10
	tag {id} 
	validExitStatus 0,141
	 
	input:
	set id, file(reads) from read_ch

	output:
	set id, "${id}.R1.fastp_trim.fastq","${id}.R2.fastp_trim.fastq" into trim_fq
	
	shell:
	'''
	#!/bin/bash
	read_len=$(zcat !{reads[0] } | head -n 2|tail -n 1 |wc -c )
	read_len=$(expr $read_len - 10)
	
	cutadapt -m $read_len -l $read_len -o !{id}.R1.fastp_trim.fastq -p !{id}.R2.fastp_trim.fastq !{reads[0]} !{reads[1]}
    '''
	
}

trim_fq.into{CIRIquant_ch; CIRI2_ch; print_ch;}


process CIRI_full{
	maxForks 10 
	
	tag {id} 
	errorStrategy 'ignore'
	publishDir "${params.outdir}/CIRI_full/", mode: 'copy'
	
	input:
	set id, file(reads_1), file(reads_2) from CIRI2_ch
	file "bwa_index/*" from bwa_indices.collect()
	file gtf
	file fasta
	
	output:
	set id, "${id}_output/CIRI-full_output/${id}_merge_circRNA_detail.anno", "${id}_output/CIRI-vis_out" into circ_full_result
	file "${id}_output/*" into res
	
	script:
	"""
	java -jar /dat1/apps/CIRI-full_v2.0/CIRI-full.jar Pipeline -1 ${reads_1 } -2 ${reads_2} -a ${gtf} -r /dat1/dat/ref/hg38/hg38+gencode.v32/bwa_index_ln/hg38.fa -d ${id}_output/ -o ${id} -t 16
	unset DISPLAY
	java -jar /dat1/apps/CIRI-full_v2.0/CIRI-vis.jar -i ${id}_output/CIRI-full_output/${id}_merge_circRNA_detail.anno -l ${id}_output/CIRI-AS_output/${id}_library_length.list -r /dat1/dat/ref/hg38/hg38+gencode.v32/bwa_index_ln/hg38.fa -d ${id}_output/CIRI-vis_out -min 1
	
	rm -rf ${id}_output/sam
    """
}



