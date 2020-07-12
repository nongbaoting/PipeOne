#!/usr/bin/env nextflow

params.reads = ""
params.sra = ""

params.genome = ""
params.cleaned  = false
params.retro_gtf = params.genome ? params.genomes[ params.genome ].retro_gtf ?: false : false
params.bowtie2_index = params.genome ? params.genomes[ params.genome ].bowtie2_index ?:false :false
retro_gtf = file(params.retro_gtf)

params.layout = 'paired'
def ifPaired = true

params.outdir = "results"
params.threads = 16
threads = params.threads

params.unstranded = true // defualt unstranded
params.reverse_stranded = true // dUTP if stranded
params.forward_stranded = false

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
				
				fastp -i ${reads}  -o ${id}.fastp.fq.gz  -q 20 --thread  1
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
				
				fastp -i ${reads[0]} -I ${reads[1]} -o ${id}.R1.fastp.fq.gz -O ${id}.R2.fastp.fq.gz -q 20 --thread  1
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

read_ch.into{ hisat2_ch;reads_print; reads_star; reads_bowtie2}
//reads_print.println()


/*
process star{
	maxForks 2
	tag { id }
	
	publishDir "${params.outdir}/star/", mode: 'copy'
	
	input:
	set id,  file(reads) from reads_star
	
	output:
	set id, "${id}.Aligned.sortedByCoord.out.bam" into star_out
	
	script:
	"""
	STAR --genomeDir $starIdex --readFilesIn ${reads} --runThreadN $threads  --outFileNamePrefix ${id}. \\
	--chimSegmentMin 10  --outFilterMultimapNmax 100 --readFilesCommand zcat  \\
	--outSAMtype BAM SortedByCoordinate  >star.log
	
	#--genomeLoad LoadAndRemove --limitBAMsortRAM 40000000000
	
	"""
	
}

process telescope_star{
	maxForks 4
	tag {id}
	errorStrategy 'ignore'
	maxRetries 2
	
	publishDir "${params.outdir}/telescope_star/", mode: 'copy', 
		saveAs: {filename -> "${id}.${filename}"}
	
	input:
	set id, file(bams) from star_out
	
	output:
	file "telescope-telescope_report.tsv"
	
	script:
	
	"""
	samtools index ${bams}
	mkdir -p temp
	telescope assign ${bams} /dsk1/ref/hg38/telescope_HERV_rmsk.hg38.v2/transcripts.gtf \\
	--max_iter 200 --theta_prior 200000  --tempdir ./temp
	"""
	
}

*/


if ( params.bowtie2_index  ){
    bowtie2_indices = Channel
        .fromPath("${params.bowtie2_index}*")
        .ifEmpty { exit 1, "bowtie2_index not found: ${params.bowtie2_index}" }
		bowtie2_base = params.bowtie2_index.split('/')[-1]
	
}


process bowtie2{
	
	tag { id }
	maxForks 10
	
	publishDir "${params.outdir}/bowtie2/", mode: 'copy',
		saveAs: {filename -> 
			if(filename =~/bam/) "bams/${filename}"
			else if(filename =~/log/) "logs/${filename}"
			}
	
	input:
	set id,  file(reads) from reads_bowtie2
	file "bowtie2_index/*" from bowtie2_indices.collect()
	
	output:
	set id, "${id}.bowtie2.sortbycoordinate.bam" into bowtie2_out
	file "${id}.bowtie2.log"
	file "stat.${id}.log" into bowtie2_out_count
	
	script:
	if(ifPaired && reads.size() == 2 ){
		"""
		set +u; source activate telescope; set -u
		bowtie2  -p ${threads} -1  ${reads[0]} -2 ${reads[1]} -x bowtie2_index/${bowtie2_base} \\
			--no-mixed --very-sensitive-local -k 100 --score-min L,0,1.6 2> ${id}.bowtie2.log | \\
		samtools sort -@ 2 - -o ${id}.bowtie2.sortbycoordinate.bam
		samtools flagstat  ${id}.bowtie2.sortbycoordinate.bam >stat.${id}.log
		"""
	}else if (! ifPaired && reads.size() == 1 ){
		"""
		set +u; source activate telescope; set -u
		bowtie2  -p ${threads} -U ${reads} -x bowtie2_index/${bowtie2_base} \\
		 --very-sensitive-local -k 100 --score-min L,0,1.6 2> ${id}.bowtie2.log | \\
		samtools sort -@ 2 - -o ${id}.bowtie2.sortbycoordinate.bam 
		samtools flagstat  ${id}.bowtie2.sortbycoordinate.bam >stat.${id}.log
		"""
	}else{
		println("input reads or configure error!" )
		println("Your reads: " + reads)
		println("Your Configure ifPaired " + ifPaired)
		exit(1)
	}

}


process merge_MapReads{
	publishDir "${params.outdir}/bowtie2/", mode: 'copy'
	
	input:
	file "bowtie2/*" from  bowtie2_out_count.collect()
	
	output:
	file "bowtie2_properPaired_total.txt" into properPaired_total
	
	
	'''
	set +u; source activate telescope; set -u
	echo "Sample\tReads" >bowtie2_properPaired_total.txt

	find bowtie2 -name "stat*log" | \
		xargs -i grep -H "properly paired" {}| \
		cut -f1 -d ' '| \
		awk '{print $1,$2}' FS=":" OFS="\t"| \
		sed -e 's#bowtie2/stat.##' -e 's#.log##' >> bowtie2_properPaired_total.txt
	'''
}


process telescope_bowtie2 {
	maxForks 10
	
	tag {id}
	errorStrategy 'ignore'
	maxRetries 2
	
	publishDir "${params.outdir}/telescope/telescope/tsv/", mode: 'copy'
		
	input:
	set id, file(bams) from bowtie2_out
	file retro_gtf
	
	output:
	file "${id}.telescope_report.tsv" into telescope_out
	
	script:
	
	"""
	set +u; source activate telescope; set -u
	samtools index ${bams}
	mkdir -p temp
	telescope assign ${bams}  ${retro_gtf} --max_iter 200 --theta_prior 200000 --tempdir ./temp
	ln -s telescope-telescope_report.tsv ${id}.telescope_report.tsv
	"""
}

process merge_telescope{
	errorStrategy 'ignore'
	publishDir "${params.outdir}/telescope/", mode: 'copy'
	
	input:
	file retro_gtf
	file "out/*" from telescope_out.collect()
	file "bowtie2_properPaired_total.txt" from properPaired_total
	
	output:
	tuple "telescope.rawCount.tsv", "transcripts.info.tsv" ,"telescope.FPKM-divide_totalMapReads.csv", "telescope.FPKM-divide_totalRetroReads.csv"
	
	script:
	"""
	python3 ${baseDir}/bin/retro.py merge telescope.rawCount.tsv  out/
	python3 ${baseDir}/bin/gtf.py to_info ${retro_gtf} transcripts.info.tsv
	Rscript ${baseDir}/bin/retro_fpkm.R 
	"""
}
