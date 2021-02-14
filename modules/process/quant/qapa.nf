
process salmon_APA {
	publishDir "${params.outdir_sub}/salmon/samples", mode: 'link',
		saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }

	tag {id } 
	input:
	tuple val(id), path(reads)
	path "transcripts_index/*" 
	
	
	output:
	path "${id}" 
	
	script:
	if(params.singleEnd){
					
		"""
		set +u; source activate pipeOne_apa; set -u
		salmon quant -i transcripts_index -l A -r $reads -p $task.cpus -o ${id}
		"""
	}else{
					
		"""
		set +u; source activate pipeOne_apa; set -u
		salmon quant -i transcripts_index -l A -1 ${reads[0] } -2 ${reads[1] } -p $task.cpus -o ${id}
		"""
	}
	
}

process qapa {
	
	publishDir "${params.outdir_sub}/apa_3utr/", mode: 'link'
	
	input:
	path utr_gtf
    path replace_SalmonIndex_ID
	path "salmon/*" 
	
	
	output:
	path "pau_results.txt" 
	
	
	"""
	set +u; source activate pipeOne_apa; set -u
	python3 ${baseDir}/bin/RNAseq/gtf.py to_gencode2biomark5 ${utr_gtf}  db_identifiers
	python3 ${baseDir}/bin/RNAseq/apa_3utr.py replace_SalmonIndex_ID ${replace_SalmonIndex_ID} salmon
	qapa quant --db db_identifiers salmon/*/quant_replace.sf > pau_results.txt
	"""
}



process qapa_filter {
	publishDir "${params.outdir_sub}/", mode: 'link'

	input:
	
	path apa_res // apa
	path TPM // TPM
	path gtf

	when: apa_res.size() >0

	output:
	path "pau_results.filterPau.txt" optional true
	path  "pau_results.filterPau-distal-proximal.txt", emit: pau // dis and pro

	script:
	"""
	set +u; source activate pipeOne_apa; set -u
	python3 ${baseDir}/bin/RNAseq/gtf.py to-info --gtf $gtf --out_tsv protein_coding_and_all_lncRNA.info.tsv
	Rscript --vanilla  ${baseDir}/bin/RNAseq/apa_3utr_filter.R $apa_res $TPM protein_coding_and_all_lncRNA.info.tsv
	"""

}