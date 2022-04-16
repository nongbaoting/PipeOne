process salmon_index {
		label 'local'

		label 'bigCPU'
		publishDir "${params.outdir_sub}/salmon/salmon_index/", mode: 'link', 
    		saveAs: { filename -> params.saveIntermediateFiles ? "$filename" : null }
		input:
		path fasta
		path gtf
		
		output:
		path "salmon_index/*" 
		
		"""
		set +u; source activate pipeOne_lncRNA; set -u
		gffread ${gtf} -g ${fasta} -w transcripts.fa -W
		salmon index -t transcripts.fa -i salmon_index -p $task.cpus
		"""
}

process salmon {
	label 'bigCPU'
	label 'local'

	publishDir "${params.outdir_sub}/salmon/samples", mode: 'copy',
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
		set +u; source activate pipeOne_lncRNA; set -u
		salmon quant -i transcripts_index -l A -r $reads -p $task.cpus -o ${id} --gcBias

		"""
	}else{
					
		"""
		set +u; source activate pipeOne_lncRNA; set -u
		salmon quant -i transcripts_index -l A -1 ${reads[0] } -2 ${reads[1] } -p $task.cpus -o ${id} 

		"""
	}
	
}

process salmon_merge {
	label 'local'
	
	publishDir "${params.outdir_sub}/salmon/", mode: 'link'
	
	input:
	path "ref_gtf" 
	path "samples/*" 
	
	output:
	path "salmon_gene_est_counts.tsv", emit: counts
	path "salmon_gene_tpm.keep.tsv", emit: keep_TPM
	path "salmon_gene_tpm.tsv", emit: TPM
	path "*tsv"
	path "gene_id.keep.txt" , emit: expr_gene_id
	
	"""
	set +u; source activate pipeOne_py3; set -u
	python3 ${baseDir}/bin/RNAseq/lncRNA.py get-txID-geneID ref_gtf protein_coding_and_all_lncRNA.txID_geneID.tsv
	Rscript ${baseDir}/bin/RNAseq/tximport_salmon.R
	# mkdir -p files_cache__
	# ls |grep -v files_cache__ |xargs -i cp -r {} files_cache__
	"""
}


workflow salmon_idx_cal {
	take:
	reads
	fasta
	gtf

	main:
	salmon_index(fasta, gtf)
	salmon(reads, salmon_index.out.collect() )
	salmon_merge(gtf, salmon.out.collect() )

	emit:
	expr_gene_id = salmon_merge.out.expr_gene_id
	counts = salmon_merge.out.counts
	keep_TPM = salmon_merge.out.keep_TPM
	TPM    = salmon_merge.out.TPM

}

