
params.sprint_tmp = false

process sprint{
	tag {id}
	publishDir "${params.outdir_sub}/sprint/", mode: 'link',
		saveAs: {filename -> 
			if (filename =~ /SPRINT_identified_all.res/) "res/${filename}"
			else "tmp/${filename}"
			
		}
	
	input:
	tuple val(id), path(reads) 
    val sprint_base
	path "sprint_index/*" 
	path sprint_repeat
	
	output:
	path "${id}" optional true 
	path "${id}.SPRINT_identified_all.res", emit: tsv
	
	"""
	set +u; source activate pipeOne_RnaEditing; set -u
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
	sprint main -rp ${sprint_repeat} -c 1 -p ${task.cpus} -1 read_1.fq -2 read_2.fq sprint_index/${sprint_base} tmp bwa samtools
	#sprint main -rp ${sprint_repeat} -c 1 -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1] } sprint_index/${sprint_base} tmp bwa samtools
	
	cp tmp/SPRINT_identified_all.res ${id}.SPRINT_identified_all.res
	
	rm read_1.fq read_2.fq
	if [[ "${params.sprint_tmp}" == "true" ]];then
		ln -s tmp ${id}
	else
		rm -rf tmp
	fi
	
	"""
}

params.miniEditing_reads = 10

process merge_sprint_A2I {
	
	publishDir "${params.outdir_sub}/sprint/Merge", mode: 'link'
		
	input:
	path "results/*"     
	
	output:
	tuple path("SPRINT_all.tsv"), path("SPRINT_A2I.tsv"), path("SPrint_A2I_table.tsv")
	
	"""
	set +u; source activate pipeOne_py3; set -u
	python3 ${baseDir}/bin/RNAseq/rnaEditing.py cat SPRINT_all.tsv SPRINT_A2I.tsv results
	
	# supported reads for A to I rna-editing is 5
	Rscript ${baseDir}/bin/RNAseq/rnaEditing_rate.R SPRINT_A2I.tsv SPrint_A2I_table.tsv ${params.miniEditing_reads}

	"""
}
