

process bowtie2_tele {
	label 'bigCPU'
	tag { id }
	
	publishDir "${params.outdir_sub}/bowtie2/", mode: 'copy',
		saveAs: {filename -> 
			if(filename =~/bam/ &&  params.saveIntermediateFiles ) "bams/${filename}"
			else if(filename =~/log/) "logs/${filename}"
			else null	
			}
	
	input:
	tuple val(id),  path(reads)
    val bowtie2_base
	path "bowtie2_index/*" 
	
	output:
	tuple val(id), path("${id}.bowtie2.sortbycoordinate.bam") , emit: bam
	path "${id}.bowtie2.log", emit: log
	path "stat.${id}.log" , emit: stat_count
	
	script:
	if(! params.singleEnd && reads.size() == 2 ){
		"""
		set +u; source activate pipeOne_telescope; set -u
		bowtie2  -p ${task.cpus} -1  ${reads[0]} -2 ${reads[1]} -x bowtie2_index/${bowtie2_base} \\
			--no-mixed --very-sensitive-local -k 100 --score-min L,0,1.6 2> ${id}.bowtie2.log | \\
		samtools sort -@ 2 - -o ${id}.bowtie2.sortbycoordinate.bam
		samtools flagstat  ${id}.bowtie2.sortbycoordinate.bam >stat.${id}.log
		"""
	}else if ( params.singleEnd ){
		"""
		set +u; source activate pipeOne_telescope; set -u
		bowtie2  -p ${task.cpus} -U ${reads} -x bowtie2_index/${bowtie2_base} \\
		 --very-sensitive-local -k 100 --score-min L,0,1.6 2> ${id}.bowtie2.log | \\
		samtools sort -@ 2 - -o ${id}.bowtie2.sortbycoordinate.bam 
		samtools flagstat  ${id}.bowtie2.sortbycoordinate.bam >stat.${id}.log
		"""
	}else{
		println("input reads or configure error!" )
		println("Your reads: " + reads)
		println("Your Configure --singleEnd " + params.singleEnd )
		exit(1)
	}

}


process merge_bowtie2_MapReads{
	publishDir "${params.outdir_sub}/bowtie2/", mode: 'link'
	
	input:
	path "bowtie2/*" 
	
	output:
	path "bowtie2_properPaired_total.txt" 
	
	
	'''
	set +u; source activate pipeOne_telescope; set -u
	echo "Sample\tReads" >bowtie2_properPaired_total.txt

	find bowtie2 -name "stat*log" | \
		xargs -i grep -H "properly paired" {}| \
		cut -f1 -d ' '| \
		awk '{print $1,$2}' FS=":" OFS="\t"| \
		sed -e 's#bowtie2/stat.##' -e 's#.log##' >> bowtie2_properPaired_total.txt
	'''
}

