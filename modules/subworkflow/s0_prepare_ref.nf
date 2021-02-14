#!/usr/bin/env nextflow
/*  @Date         : 2021/02/01 13:19:45
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/

nextflow.enable.dsl=2
// defined params

// include funtions must be placed after params !
include {Header; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
//include {  } from  './modules/modules.nf'

params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false


// bowtie2 index for s4_retrotranscriptome telescope 
workflow Prepare_References {
    fasta = check_file(params.fasta, '--fasta')
    gtf   = check_file(params.gtf, '--gtf')
    bowtie2_build(fasta)
    sprint_index(fasta, gtf)
    star_index(fasta, gtf)
    hisat2_build(fasta, gtf)
    bwa_index(fasta)
}

process bowtie2_build {
    label 'prepare'
    
    publishDir "./", mode: 'copy'
    input:
    path fasta

    output:
    path "bowtie2_index"

    """
    set +u; source activate pipeOne_telescope; set -u
    mkdir -p bowtie2_index
    bowtie2-build --threads ${task.cpus} ${fasta} bowtie2_index/bowtie2_base
    """
}

// bwa index for s2_circRNA, s6_RNAediting
process sprint_index {
    label 'prepare'

    publishDir "./", mode: 'copy',
        saveAs: {filename -> "$filename"}
    input:
    path "genome.fa" 
    path "genome.gtf"

    output:
    path "sprint_index/genome*"

    """
    set +u; source activate pipeOne_RnaEditing; set -u
    mkdir -p sprint_index; cd sprint_index
    cp ../genome.fa genome.fa
    sprint prepare -t ../genome.gtf genome.fa bwa
    """
    }

process star_index {
    label 'prepare'
    label 'bigMEM'
    label 'bigCPU'

    publishDir "./", mode: 'copy'
    input:
    path "genome.fa" 
    path "genome.gtf"  

    output:
    path "STAR_index"

    """
    set +u; source activate pipeOne_py3; set -u
    mkdir STAR_index; cd STAR_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir . --genomeFastaFiles ../genome.fa --sjdbGTFfile ../genome.gtf
    """
}

process hisat2_build {
    label 'prepare'
    label 'bigMEM'

    publishDir "./", mode: 'copy'
    input:
    path "genome.fa" 
    path "genome.gtf"

    output:
    path "hisat2_index"

    """
    set +u; source activate pipeOne_lncRNA; set -u
    # hisat2_extract_splice_sites.py genome.gtf >splice.txt
    # hisat2_extract_exons.py genome.gtf >exon.txt
    mkdir hisat2_index; cd hisat2_index
    # hisat2-build --exon ../exon.txt --ss ../splice.txt -p ${task.cpus}  ../genome.fa genome
    hisat2-build  -p ${task.cpus}  ../genome.fa genome
    """
}

process bwa_index {
    label 'prepare'
    publishDir "./", mode: 'copy'
    input:
    path fasta

    output:
    path "bwa_index"

    """
    set +u; source activate pipeOne_RnaEditing; set -u
    mkdir -p bwa_index
    bwa index -p bwa_index/bwa_index ${fasta}
    """
}

