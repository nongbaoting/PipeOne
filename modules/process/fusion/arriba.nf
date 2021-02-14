process arriba{
	label 'bigMEM'
    label 'bigCPU'
    errorStrategy 'ignore'
    
    tag { id }
    publishDir "${params.outdir_sub}/arriba/", mode: 'link',
        saveAs: {filename ->
        if( filename =~ /.fusions.tsv/ ) "tsv/${filename}"
            else if(filename =~ /fusions.pdf/) "pdf/${filename}"
            else null
            }
        
    input:
    tuple  val(id), path(reads)
    path "STARIndex/*"
    path gtf
    path fasta
    path blacklisted
    path cytobands
    path proteinDomains
    
    output:
    path "${id}.fusions.tsv" , emit: tsv
    //path "${id}.fusions.pdf" optional true
    
    shell:
    '''
    set +u; source activate pipeOne_fusion; set -u
    STAR \
        --runThreadN !{task.cpus} \
        --genomeDir STARIndex --genomeLoad NoSharedMemory \
        --readFilesIn !{reads} --readFilesCommand zcat \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
        --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
        --chimSegmentMin 10 --chimOutType WithinBAM SoftClip \
        --chimJunctionOverhangMin 10 --chimScoreMin 1 \
        --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 | \
    arriba \
        -x /dev/stdin \
        -o !{id}.fusions.tsv -O fusions.discarded.tsv \
        -a !{fasta} -g !{gtf} -b !{blacklisted} \
        -T -P
        
    draw_fusions.R \
        --fusions=!{id}.fusions.tsv \
        --output=!{id}.fusions.pdf \
        --annotation=!{gtf} \
        --cytobands=!{cytobands} \
        --proteinDomains=!{proteinDomains}
    '''
		
	}

process merge_arriba {
    publishDir "${params.outdir_sub}/arriba/", mode: 'link'
    
    input:
    path "arriba_out/*" 
    
    output:
    path "fusion_arriba_out.tsv"
    
    """
    set +u; source activate pipeOne_fusion; set -u
    python3 ${baseDir}/bin/RNAseq/fusion.py arriba_table fusion_arriba_out.tsv arriba_out/
    """
}