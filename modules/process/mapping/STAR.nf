process star_1_pass_AS {
    label 'bigCPU'
    label 'bigMEM'
    publishDir "${params.outdir_sub}/STAR/star1pass", mode: 'link'

    tag {id}
    input:
    tuple val(id), path(reads) 
    path "STARIndex/*" 
    
    output:
    
    path "${id}.SJ.out.tab" 
        
    """
    set +u; source activate pipeOne_py3; set -u
    STAR 	--runThreadN ${task.cpus } \\
            --genomeDir STARIndex \\
            --outSAMtype None \\
            --readFilesIn ${reads} \\
            --outFileNamePrefix	${id}. \\
            --readFilesCommand zcat \\
            --outFilterMultimapNmax 50 \\
            --outFilterMultimapScoreRange 3 \\
            --outFilterScoreMinOverLread 0.7 \\
            --outFilterMatchNminOverLread 0.7 \\
            --outFilterMismatchNmax 10 \\
            --alignIntronMax 500000 \\
            --alignMatesGapMax 1000000 \\
            --sjdbScore 2
    
    """
    
}

	// Aligned.sortedByCoord.out.bam
	// Aligned.out.bam
	// Aligned.out.sam


process star_genomeGenerate_for_2pass{
    label 'bigCPU'
    label 'bigMEM'

    input:
    path fasta
    path gtf
    path "star_junction/*" 
    
    output:
    path "STARIndex_2pass/*" 
    
    """
    set +u; source activate pipeOne_py3; set -u
    find star_junction -type f |xargs -i cat {} > SJ.out.tab
    genomeDir=STARIndex_2pass
    mkdir \$genomeDir
    STAR --runMode genomeGenerate --genomeDir \$genomeDir \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
            --sjdbFileChrStartEnd SJ.out.tab  --runThreadN ${task.cpus }
    
    """
}


process star_2_pass_AS{
    label 'bigCPU'
    label 'bigMEM'
    tag {id}
    
    errorStrategy 'ignore'

    publishDir "${params.outdir_sub}/STAR/star2pass", mode: 'link'
    
    input:
    tuple  val(id), path(reads) 
    path "STARIndex/*" 
    
    output:
    tuple val(id), path("${id}.bam"), path("${id}.bam.bai" ), emit: id_bam_bai
    tuple val(id), path("${id}.bam"), emit: id_bam
    tuple path("${id}.bam"), path("${id}.bam.bai" ), emit: bam_bai

    
    
    """
    set +u; source activate pipeOne_py3; set -u
    STAR --runThreadN ${task.cpus } \\
        --genomeDir STARIndex \\
        --readFilesIn ${reads} \\
        --outFileNamePrefix ${id}. \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMstrandField intronMotif \\
        --outSAMattributes NH HI NM MD AS XS \\
        --outSAMheaderHD @HD VN:1.4 \\
        --outFilterMultimapNmax 50 \\
        --outFilterMultimapScoreRange 3 \\
        --outFilterScoreMinOverLread 0.7 \\
        --outFilterMatchNminOverLread 0.7 \\
        --outFilterMismatchNmax 10 \\
        --alignIntronMax 500000 \\
        --alignMatesGapMax 1000000 \\
        --sjdbScore 2 --limitBAMsortRAM 1443648494946
        
    mv ${id}.Aligned.sortedByCoord.out.bam ${id}.bam
    samtools index ${id}.bam
    
    """
}

workflow STAR_AS {
    take:
    reads
    fasta
    gtf
    star_indices

    main:
    star_1_pass_AS(reads, star_indices.collect() )
    star_genomeGenerate_for_2pass(fasta, gtf, star_1_pass_AS.out.collect() ) 
    star_2_pass_AS(reads, star_genomeGenerate_for_2pass.out.collect() )

    emit:
    id_bam = star_2_pass_AS.out.id_bam
    id_bam_bai = star_2_pass_AS.out.id_bam_bai
    bam_bai = star_2_pass_AS.out.bam_bai


}