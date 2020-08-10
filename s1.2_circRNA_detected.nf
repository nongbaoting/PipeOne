#!/usr/bin/env nextflow

params.reads = ""
params.sra   = ""
params.single = ""
params.genome =""
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.ref = params.genome ? params.genomes[ params.genome ].ref ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?:false :false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false
params.repeat_gtf=  params.genome ? params.genomes[ params.genome ].repeat_gtf?:false :false
params.threads = 12
params.saveIntermediate =  false

threads = params.threads
bwa_index = params.bwa_index
scripts = Channel.fromPath("$baseDir/scripts/*")
scripts.into{scripts_1; scripts_2}

def check_file(myfile, par_name="your input file "){
    if ( myfile ){
        myfile_ = file(myfile)
        if( !myfile_.exists() ) exit 1, "File not found: ${myfile}"

    }else{exit 1, "value ${par_name} is null !"}
    return myfile_
}


repeat_gtf= check_file(params.repeat_gtf, par_name = "--repeat_gtf")

if ( params.fasta ){
    fasta = check_file(params.fasta)
	fasta_fai = check_file("${params.fasta}.fai")
}


if ( params.gtf ){
    gtf = check_file(params.gtf)

}
 


if ( params.ref ){
    ref = check_file(params.ref)
}


if ( params.bwa_index  ){
    bwa_indices = Channel
        .fromPath("${params.bwa_index}.*")
        .ifEmpty { exit 1, "bwa index not found: ${params.bwa_index}" }
	bwa_base = params.bwa_index.split('/')[-1]
}

if( params.star_index ){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    star_index.into{star_index; star_index_1; star_index_2}
}

// input reads or sra ! s
if(params.reads ){
	
	// defined single or paired reads !
	if(!params.single){
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
		println "your read file name must not have the pattern *.fastp.fq.gz which will cause an disaster"
		process fastp{
			tag {id}
			
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
			if(params.single){
				"""
				if [ "${reads}" = "${id}.fastp.fq.gz" ];
				then
					echo "input raw reads name and output clean name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 --thread  2
				"""
			}else
			{
				"""
				if [ "${reads[0]}" = "${id}.R1.fastp.fq.gz" ];
				then
					echo "input raw read name and output clean read name are identical, please check your input files!"
					exit 1
				fi
				
				fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  2
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
	print_ch.println()
	process fastq_dump_and_fastp{
		tag {id}
		
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
		if(params.single){
			"""
			fastq-dump --split-3  $sra
			fastp -i ${id}_1.fastq  -o ${id}.fastp.fq.gz -O -q 20 --thread  4
			"""
		}else{
			"""
			fastq-dump --split-3  $sra
			fastp -i ${id}_1.fastq -I ${id}_2.fastq -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  4
			"""
		}
		
	}
	
	//read_ch=fastp_out.map{it ->a=[it[0],it[1 .. 2]] }
}else{
	exit 1, "you should input --reads or --sra "
}

read_ch.into{CIRI2_ch; circexplorer2_ch; circRNA_finder_ch; find_circ_ch; dcc_ch; count_ch }

process dcc{
    maxForks 1
 	publishDir "${params.outdir}/DCC/", mode: 'copy',
    	saveAs: { filename ->
			if( (filename =~ /${id}.Aligned.sortedByCoord.out.bam/ || filename == "${id}.Log.final.out" || filename == "${id}.Chimeric.out.junction")  &&  params.saveIntermediate  )  "DCC_BAM/${filename}"
			else  if ( filename	=~ /Count|Circ/) "DCC/${filename}"
			else null
				}
   
   tag {id }
	
	label 'big_mem'	
  
	input:
	file "STARIndex" from  star_index_2.collect()
	set id, file(reads) from dcc_ch
	file fasta
	file fasta_fai
	file gtf
	file repeat_gtf
    
	output:
	tuple  id, "${id}.Chimeric.out.junction" into dcc_jun_out
	tuple "${id}.Aligned.sortedByCoord.out.bam*", "${id}.Log.final.out" optional true
	tuple "${id}.CircRNACount", "${id}.CircCoordinates", "${id}.LinearCount", "${id}.CircSkipJunctions" into dcc_result
	
	"""
	set +u; source activate circRNA; set -u
	mkdir -p ${id}/joint_mapping/ ${id}_mate1/joint_mapping/  ${id}_mate2/joint_mapping/ 
	STAR --runThreadN $threads \
       --genomeDir STARIndex \
       --outSAMtype BAM SortedByCoordinate \
	   --readFilesIn ${reads} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${id}/joint_mapping/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 
	if [ ${params.saveIntermediate} != false ]
		then
		ln -s ${id}/joint_mapping/Aligned.sortedByCoord.out.bam ${id}.Aligned.sortedByCoord.out.bam
		ln -s ${id}/joint_mapping/Aligned.sortedByCoord.out.bam.bai ${id}.Aligned.sortedByCoord.out.bam.bai
		ln -s ${id}/joint_mapping/Log.final.out                 ${id}.Log.final.out
	fi
	ln -s ${id}/joint_mapping/Chimeric.out.junction         ${id}.Chimeric.out.junction


	STAR --runThreadN $threads \
       --genomeDir STARIndex \
       --outSAMtype None \
       --readFilesIn ${reads[0]} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${id}_mate1/joint_mapping/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15 
	   
	STAR --runThreadN $threads \
       --genomeDir STARIndex \
       --outSAMtype None \
       --readFilesIn ${reads[1]} \
       --readFilesCommand zcat \
       --outFileNamePrefix ${id}_mate2/joint_mapping/ \
       --outReadsUnmapped Fastx \
       --outSJfilterOverhangMin 15 15 15 15 \
       --alignSJoverhangMin 15 \
       --alignSJDBoverhangMin 15 \
       --seedSearchStartLmax 30 \
       --outFilterMultimapNmax 20 \
       --outFilterScoreMin 1 \
       --outFilterMatchNmin 1 \
       --outFilterMismatchNmax 2 \
       --chimSegmentMin 15 \
       --chimScoreMin 15 \
       --chimScoreSeparation 10 \
       --chimJunctionOverhangMin 15
	   
	samtools index  ${id}/joint_mapping/Aligned.sortedByCoord.out.bam
	
	DCC ${id}/joint_mapping/Chimeric.out.junction -mt1 ${id}_mate1/joint_mapping/Chimeric.out.junction  \
	-mt2 ${id}_mate2/joint_mapping/Chimeric.out.junction  \
	-T $threads -D -N -R $repeat_gtf -an ${gtf} -Pi -F -M -fg -G -A ${fasta} -Nr 1 1
	
	ln -s CircRNACount      ${id}.CircRNACount 
	ln -s CircCoordinates   ${id}.CircCoordinates 
	ln -s LinearCount       ${id}.LinearCount 
	ln -s CircSkipJunctions ${id}.CircSkipJunctions 
	"""
}

dcc_jun_out.into{circRNA_finder_junc; circexplorer2_junc}

process circexplorer2{
	publishDir "${params.outdir}/circexplorer2/", mode: 'copy'
	
	input:
	set id, file(junc) from circexplorer2_junc
	file gtf
	file ref
	file fasta
	
	output:
	file "${id}.circularRNA_known.txt" into circexplorer2_result
	
	script:
	"""
	set +u; source activate circRNA; set -u
	CIRCexplorer2 parse -t STAR $junc > CIRCexplorer2_parse.log
	CIRCexplorer2 annotate -r $ref -g $fasta -b back_spliced_junction.bed -o ${id}.circularRNA_known.txt
	"""
}

process CIRI2{
	label 'big_mem'
	tag {id} 
	publishDir "${params.outdir}/CIRI2/", mode: 'copy'
	
	input:
	set id, file(reads) from CIRI2_ch
	file "bwa/*" from bwa_indices.collect()
	file gtf
	file fasta
	
	output:
	file "${id}.ciri2.txt" into circ2_result
	script:
	
	"""
	set +u; source activate circRNA; set -u
	bwa mem -t $threads -T 19 bwa/$bwa_base $reads > aln-bwa-mem.sam
	perl /opt/CIRI_v2.0.6/CIRI2.pl -T 8 -I aln-bwa-mem.sam -O ${id}.ciri2.txt -F ${fasta} -A ${gtf}
    rm aln-bwa-mem.sam
    """
}

process merge_circ{
	label 'py3'
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/merge_circ/", mode: 'copy'
	
	input:
	file "circ2/*" from circ2_result.collect()
	file "circexplorer2/*" from circexplorer2_result.collect()
	file "dcc/*" from dcc_result.collect()
	file "scripts/*" from scripts_1.collect()
	file  "ref.txt" from ref
	
    output:
    tuple path("all_circRNAs.hc.bed"), path('at_least_2_softs.bed') into circ_bed, circ_bed_plot
	// set "all_circ_anno.txt","exon_nums_size"
	
	"""
	merge_circRNA.py ciri2          ciri2.table.txt          circ2/*
	merge_circRNA.py circexplorer2  circexplorer2.table.txt  circexplorer2/*
	merge_circRNA.py dcc            dcc.table.txt            dcc/*CircRNACount 
    merge_circRNA.py dcc            dcc_linear.table.txt     dcc/*LinearCount
    merge_circRNA.py dcc_bed        dcc/*CircCoordinates
    
    merge_circRNA.py merge-tables all_circRNAs.table circexplorer2.table.txt ciri2.table.txt dcc.table.txt
    merge_circRNA.py merge-bed ciri2.bed dcc.bed circexplorer2.bed
    mkdir -p tmp
    cat ciri2.bed |cut -f1,2,3,6|awk '{print \$1,\$2,\$3,\$4}' OFS="__" FS='\t' >tmp/ciri2.bed.tmp
    cat dcc.bed |cut -f1,2,3,6|awk '{print \$1,\$2,\$3,\$4}' OFS="__" FS='\t' >tmp/dcc.bed.tmp
    cat circexplorer2.bed |cut -f1,2,3,6|awk '{print \$1,\$2,\$3,\$4}' OFS="__" FS='\t' > tmp/circexplorer2.bed.tmp
	
	Rscript ${baseDir}/scripts/venn.R
	#cat circ2/* >all_ciri2.txt
	#cat circexplorer2/* >all_ce2.txt
	#python3 scripts/anno_by_ciri2_ce2.py anno-ciri2-ce2  all_ciri2.txt all_ce2.txt all_circ_anno.txt
	python3 ${baseDir}/bin/bed.py select-common 1 ciri2_ce2.bed  ciri2.bed circexplorer2.bed
	python3 ${baseDir}/bin/bed.py venn ciri2_ce2.bed at_least_2.circRNA at_least_2_softs.bed
	
	#python3 scripts/stat.py exon-nums-size ref.txt at_least_2_softs.bed circexplorer2/*txt
	
	#mkdir -p files_cache__
	#ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
	#cp .command.sh files_cache__/command.sh
    #  */
	
	"""
}
