#!/usr/bin/env nextflow

ref_directory = "/home/dat/ref/PipeOne/hg38+gencode.v32"
params.fasta  = "${ref_directory}/hg38.fa"
genecode_gtf  			= "${ref_directory}/gencode.v32.gtf"


params.threads = 48
threads = params.threads
println params.fasta

fasta = check_file(params.fasta)
gtf = check_file(params.gtf)


// bowtie2 index for s4_retrotranscriptome telescope 

process bowtie2_build {
    publishDir "./", mode: 'copy'
    input:
    file fasta

    output:
    "bowtie2_index"

    """
    set +u; source activate telescope; set -u
    mkdir -p bowtie2_index
    bowtie2-build --threads ${threads}  ${fasta} bowtie2_index/bowtie2_base
    """
}

// bwa index for s2_circRNA, s6_RNAediting

process sprint_index {
    publishDir "./", mode: 'copy'
    input:
    file "genome.fa" from fasta
    file "genome.gtf" from gtf 

    output:
    "sprint_index"

    """
    mkdir -p sprint_index
    sprint prepare -t genome.gtf genome.fa bwa
    """
    
}

/*

process bwa_index{
    publishDir "./", mode: 'copy'
    input:
    file fasta

    output:
    "bowtie2_index"

    """
    set +u; source activate RnaEditing; set -u
    mkdir -p bwa_index
    bwa index  --threads ${threads}  ${fasta} -p bwa_index/bwa_index
    """
    
}

process hisat2_build {

}

process sprint_build {
    
}
*/


def check_file(myfile){
    if ( myfile ){
        myfile_ = file(myfile)
        if( !myfile_.exists() ) exit 1, "File not found: ${myfile}"

    }else{exit 1, "value is null !"}

    return myfile_
}