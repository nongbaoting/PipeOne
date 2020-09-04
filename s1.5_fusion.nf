#!/usr/bin/env nextflow

params.arriba = "true"
params.star_fusion = "false"
params.update_GTF = false

params.reads = ""

params.genome =""
params.cleaned  = false
params.saveIntermediateFiles = false

params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false
params.blacklisted =  params.genome ? params.genomes[ params.genome ].blacklisted ?:false :false
params.proteinDomains =  params.genome ? params.genomes[ params.genome ].proteinDomains ?:false :false
params.cytobands =  params.genome ? params.genomes[ params.genome ].cytobands ?:false :false
params.ctat_dir =  params.genome ? params.genomes[ params.genome ].ctat_dir ?:false :false
params.threads = 8


blacklisted = file(params.blacklisted)
proteinDomains = file(params.proteinDomains)
cytobands = file(params.cytobands)

threads = params.threads

scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2; scripts_3; scripts_4}

if ( params.fasta ){
    fasta = file(params.fasta)
	fasta_fai = file("${params.fasta}.fai")
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}


if ( params.gtf && params.update_GTF == false ){
    gtf = file(params.gtf)
    if( !gtf.exists() ) exit 1, "GTF file not found: ${params.gtf}"
}else if(params.update_GTF == true){
	gtf = file("../s1.1_lncRNA/results/annotations_and_fasta/protein_coding_and_all_lncRNA.gtf")
	if( !gtf.exists() ) exit 1, "file: ../s1.1_lncRNA/results/annotations_and_fasta/protein_coding_and_all_lncRNA.gtf does not found\
	\nPlease check step s1.1_lncRNA has complete?"
}else{
	exit 1, "GTF file not found!"
}

if( params.star_index ){
    star_index = Channel
        .fromPath("${params.star_index}/*")
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    star_index.into{star_index; star_index_1; star_index_2}
}



params.layout = 'paired'
def ifPaired = true
if (params.layout =~ /single/ ){
	ifPaired = false
}


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
	exit 1, "you should input --reads or --sra or --bam"
}

read_ch.into{ reads_arriba; reads_star_fusion; reads_print}
//reads_print.println()

if (params.arriba == "true" ){
	
	process arriba{
		
		errorStrategy 'ignore'
		
		tag { id }
		publishDir "${params.outdir}/arriba/", mode: 'link',
		  saveAs: {filename ->
			if( filename =~ /.fusions.tsv/ ) "tsv/${filename}"
				else if(filename =~ /fusions.pdf/) "pdf/${filename}"
				else null
				}
		  
		input:
		file "STARIndex/*" from  star_index.collect()
		set  id, file(reads) from reads_arriba
		file gtf
		file fasta
		file blacklisted
		file cytobands
		file proteinDomains
		
		output:
		file "${id}.fusions.tsv" into arriba_out
		file "${id}.fusions.pdf" optional true
		
		shell:
		'''
		set +u; source activate fusion; set -u
		STAR \\
			--runThreadN !{threads} \\
			--genomeDir STARIndex --genomeLoad NoSharedMemory \\
			--readFilesIn !{reads} --readFilesCommand zcat \\
			--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \\
			--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \\
			--chimSegmentMin 10 --chimOutType WithinBAM SoftClip \\
			--chimJunctionOverhangMin 10 --chimScoreMin 1 \\
			--chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 \\
			--chimScoreSeparation 1 \\
			--alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 | \\
		arriba \\
			-x /dev/stdin \\
			-o !{id}.fusions.tsv -O fusions.discarded.tsv \\
			-a !{fasta} -g !{gtf} -b !{blacklisted} \\
			-T -P
			
		draw_fusions.R \\
			--fusions=!{id}.fusions.tsv \\
			--output=!{id}.fusions.pdf \\
			--annotation=!{gtf} \\
			--cytobands=!{cytobands} \\
			--proteinDomains=!{proteinDomains}
		'''
		
	}

	process merge_arriba{
		publishDir "${params.outdir}/arriba/", mode: 'link'
		
		input:
		path "arriba_out/*" from arriba_out.collect()
		
		output:
		path "fusion_arriba_out.tsv"
		
		"""
		python3 ${baseDir}/bin/fusion.py arriba_table fusion_arriba_out.tsv arriba_out/
		"""
	}

}


if(params.star_fusion == "true"){
	
	Channel
	.fromPath("${params.ctat_dir}/*", type: 'any')
	.ifEmpty { exit 1, "ctat_dir not found: ${params.ctat_dir}" }
	.set{CTAT_resource_lib}

		
	process star_fusion{
		errorStrategy 'ignore'
		
		tag { id }
		publishDir "${params.outdir}/star_fusion/", mode: 'link'
				
		input:
		set id , file(reads) from reads_star_fusion
		file "CTAT_resource_lib/*" from CTAT_resource_lib.collect()
		
		output:
		file "star_fusion_outdir" into star_fusion_out
		
		
		script:
		"""
		set +u; source activate fusion; set -u
		
		STAR-Fusion --genome_lib_dir CTAT_resource_lib \
				 --left_fq ${reads[0]} \
				 --right_fq ${reads[1]} \
				 --output_dir star_fusion_outdir
				 
		"""
	}

}



/*
process arriba {
	errorStrategy 'ignore'
	tag { id }
	publishDir "${params.outdir}/arriba/", mode: 'link',
      saveAs: {filename ->
		if( filename =~ /.fusions.tsv/ ) "tsv/${filename}"
		else if(filename =~ /fusions.pdf/) "pdf/${filename}"
        else null
			}
			
	input:
	set id, file(bam) from star_out_bam
	file gtf
    file fasta
    file blacklisted
	file cytobands
	file proteinDomains
	
	output:
	file "${id}.fusions.tsv" into arriba_out
	file "${id}.fusions.pdf" optional true
	
	script:
	"""
	/dat1/apps/arriba_v1.1.0/arriba \
		-x ${bam} \
		-o ${id}.fusions.tsv -O fusions.discarded.tsv \
		-a ${fasta} -g ${gtf} -b ${blacklisted} \
		-T -P
	/dat1/apps/arriba_v1.1.0/draw_fusions.R \
    --fusions=${id}.fusions.tsv \
    --output=${id}.fusions.pdf \
    --annotation=${gtf}\
    --cytobands=${cytobands} \
    --proteinDomains=$proteinDomains
	"""
	
}
*/