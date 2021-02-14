

process quant {
	tag {id }
	label 'bigCPU'
	label 'bigMEM'
	
	publishDir "${params.outdir_sub}/CIRIquant/samples", mode: 'link',
		saveAs: {filename ->
				if(filename == id ) null
				else filename 
				}
	
	input:
	tuple val(id), path(bam), path(reads) 
	path fasta
	path gtf
	val bwa_base
	path "bwa_index/*"
	val hisat2_base
	path "hisat2_index/*"
	
	output:
	path "${id}", emit: sample
	tuple path("${id}/circ/${id}.ciri"), path("${id}/gene/"), path("${id}/${id}.*") 
	
	"""
	set +u; source activate pipeOne_CIRIquant; set -u
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
	CIRIquant -t $task.cpus \
          -1 ${reads[0] } \
          -2 ${reads[1] } \
          --config ./config.yml \
          -o ${id} \
          -p ${id} \
		  --bam ${bam} 
	"""
}

process merge{

	publishDir "${params.outdir_sub}/CIRIquant/", mode: 'link'
	
	input:
	path "CIRIquant/*" 

	output:
	tuple path("circRNA_info.csv"), path("circRNA_bsj.csv"), path("circRNA_ratio.csv"),path("library_info.csv"), path("gene_count_matrix.csv"), emit: res_tab
	tuple path("sample_gene.lst"), path("sample_psuedo.lst")
 	

	shell:
	''' 
	set +u; source activate pipeOne_CIRIquant; set -u
	find -L ./CIRIquant -name "*_out.gtf"|sort |awk '{print $3, $0}' FS="/" OFS="\t" > sample_gene.lst
	find -L ./CIRIquant -name "*.gtf"|awk '{print $3,"./CIRIquant/" $3 "/" $3 ".gtf", "C" }'  FS="/" OFS="\t" >sample_psuedo.lst

	prep_CIRIquant -i sample_psuedo.lst \
                 --lib library_info.csv \
                 --circ circRNA_info.csv \
                 --bsj circRNA_bsj.csv \
                 --ratio circRNA_ratio.csv

	prepDE.py -i sample_gene.lst
	
	'''
}

process CPM {
	publishDir "${params.outdir_sub}/CIRIquant/", mode: 'link'

	input:
	tuple path("circRNA_info.csv"), path("circRNA_bsj.csv"), path("circRNA_ratio.csv"),path("library_info.csv"), path("gene_count_matrix.csv") 
	
	output:
	path "circRNA_cpm.tsv"

	"""
	set +u; source activate pipeOne_py3; set -u
	Rscript  ${baseDir}/bin/RNAseq/circRNA_cpm.R
	"""

}

workflow CIRIquant {
	take:
	reads
	bams
	fasta
	gtf
	bwa_base
	bwa_indices
	hisat2_base
	hisat2_indices

	main:
	bams
		.combine(reads, by: 0)
		.set{ id_bams_reads }
	quant(id_bams_reads, fasta, gtf, bwa_base, bwa_indices, hisat2_base, hisat2_indices )
	merge(quant.out.sample.collect() ) 
	CPM(merge.out.res_tab )

	emit:
	CPM.out
}