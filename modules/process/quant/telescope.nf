
process telescope_bowtie2 {
	
	tag {id}
	// errorStrategy 'ignore'
	maxRetries 2
	
	publishDir "${params.outdir_sub}/telescope/telescope/tsv/", mode: 'copy'
		
	input:
	tuple val(id), path(bams) 
	path retro_gtf
	
	output:
	path "${id}.telescope_report.tsv" 

	
	script:
	
	"""
	set +u; source activate pipeOne_telescope; set -u
	samtools index ${bams}
	mkdir -p temp
	telescope assign ${bams}  ${retro_gtf} --max_iter 200 --theta_prior 200000 --tempdir ./temp
	ln -s telescope-telescope_report.tsv ${id}.telescope_report.tsv
	"""
}

process merge_telescope{
	errorStrategy 'ignore'
	publishDir "${params.outdir_sub}/telescope/", mode: 'copy'
	
	input:
	path retro_gtf
    path "bowtie2_properPaired_total.txt" 
	path "tele_res/*" 
	
	
	output:
	tuple path("telescope.rawCount.tsv"), path("transcripts.info.tsv") 
    path "telescope.FPKM-divide_totalMapReads.tsv", emit: expr
	path "telescope.FPKM-divide_totalRetroReads.tsv" 

	script:
	"""
	set +u; source activate pipeOne_py3; set -u
	python3 ${baseDir}/bin/RNAseq/retro.py merge telescope.rawCount.tsv  tele_res/
	python3 ${baseDir}/bin/RNAseq/gtf.py to_info ${retro_gtf} transcripts.info.tsv
	Rscript ${baseDir}/bin/RNAseq/retro_fpkm.R 
	"""
}
