#!/usr/bin/env nextflow

params.reads = ""
// "/public1/pub/guohh/other/mouse_fatty/00_rawdata/*.fastq.gz"
params.sra   = ""
params.hisat2_bam = ""
params.genome =""
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?:false :false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false
params.threads = 12
params.saveIntermediateFiles = false

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
	exit 1, "you should input --reads or --sra"
}



read_ch.into{CIRIquant_ch; hisat2_ch; print_ch;}

if( params.hisat2_bam) {
	 Channel
		.fromPath(params.hisat2_bam)
		.ifEmpty{error "Can not find any bam files  matching:${params.hisat2_bam}"}
		.map{ file ->
			key = file.name.toString().tokenize('.').get(0)
			return tuple(key, file)
			}
		.groupTuple()
		.set{ hisat2_bam }
}else{
	
	params.unstranded = true // defualt unstranded
	params.reverse_stranded = true // dUTP if stranded
	params.forward_stranded =false

	process hisat2 {
			
		tag {id}
		publishDir "${params.outdir}/hisat2/", mode: 'link',
			saveAs: {filename -> 
				if(filename =~ /bam/ && params.saveIntermediateFiles ) filename
				else if (filename =~/log/) "logs/${filename}"
				else null
				}
		
		input:
		set id, file(reads)   from hisat2_ch
		file "hisat2_index/*" from hisat2_indices.collect()
		
		output:
		set id, "${id}.hisat2.sortbycoordinate.bam" into hisat2_bam
		file "${id}.hisat2.log" into hisat2_log
		
		script:
		def rnastrandness = ''
        if (params.forward_stranded && !params.unstranded){
            rnastrandness = ifPaired ? '--rna-strandness FR' : '--rna-strandness F'
        } else if (params.reverse_stranded && !params.unstranded){
            rnastrandness = ifPaired ? '--rna-strandness RF' : '--rna-strandness R' 
        }
		
		if( ifPaired ){
			"""
			set +u; source activate lncRNA; set -u
			hisat2 -p $threads --dta  $rnastrandness  -x hisat2_index/$hisat2_base \
			-1 ${reads[0]} -2 ${reads[1]}   2> ${id}.hisat2.log | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam 
			samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
			"""
		}else{
			
			"""
			set +u; source activate lncRNA; set -u
			hisat2 -p $threads --dta  $rnastrandness -x  hisat2_index/$hisat2_base \
			-U  ${reads}   2> ${id}.hisat2.log  | samtools sort -@ 8 - -o ${id}.hisat2.sortbycoordinate.bam
			samtools index ${id}.hisat2.sortbycoordinate.bam
			conda deactivate
			"""
		}

	}

}

hisat2_bam
	.combine(CIRIquant_ch, by: 0)
	.set{CIRIquant_combine_ch }

process CIRIquant {
	tag {id }
	publishDir "${params.outdir}/CIRIquant/samples", mode: 'link',
		saveAs: {filename ->
				if(filename == id ) null
				else filename 
				}
	
	input:
	tuple id, file(bam), file(reads) from CIRIquant_combine_ch
	path fasta
	path gtf
	file "bwa_index/*"    from bwa_indices.collect()
	file "hisat2_index/*" from hisat2_indices_2.collect()
	
	output:
	file "${id}" into  CIRIquant_out
	set "${id}/circ/${id}.ciri", "${id}/gene/", "${id}/${id}.*" into CIRIquant_file
	
	"""
	set +u; source activate CIRIquant; set -u
	echo "name: ${params.genome}
	tools:
	  bwa: \$(which bwa)
	  hisat2:  \$(which hisat2)
	  stringtie: \$(which stringtie)
	  samtools: \$(which samtools)

	reference:
	  fasta: ${fasta}
	  gtf: ${gtf}
	  bwa_index: bwa_index/${bwa_base}
	  hisat_index: hisat2_index/${hisat2_base}
	" >config.yml
	
	samtools index ${bam}
	CIRIquant -t ${threads} \
          -1 ${reads[0] } \
          -2 ${reads[1] } \
          --config ./config.yml \
          -o ${id} \
          -p ${id} \
		  --bam ${bam} \
		  -t $threads
	"""
}




process merge_CIRIquant{

	publishDir "${params.outdir}/CIRIquant/", mode: 'link'
	
	input:
	path "CIRIquant/*" from CIRIquant_out.collect()

	output:
	tuple path("circRNA_info.csv"), path("circRNA_bsj.csv"), path("circRNA_ratio.csv"), path("gene_count_matrix.csv")
	tuple path("sample_gene.lst"), path("sample_psuedo.lst")
 	path "circRNA_cpm.csv"

	shell:
	''' 
	set +u; source activate CIRIquant; set -u
	find -L ./CIRIquant -name "*_out.gtf"|sort |awk '{print $3, $0}' FS="/" OFS="\t" > sample_gene.lst
	find -L ./CIRIquant -name "*.gtf"|awk '{print $3,"./CIRIquant/" $3 "/" $3 ".gtf", "C" }'  FS="/" OFS="\t" >sample_psuedo.lst

	prep_CIRIquant -i sample_psuedo.lst \
                 --lib library_info.csv \
                 --circ circRNA_info.csv \
                 --bsj circRNA_bsj.csv \
                 --ratio circRNA_ratio.csv

	prepDE.py -i sample_gene.lst
	Rscript  !{baseDir}/bin/circRNA_cpm.R
	'''
}


