#!/usr/bin/env nextflow

ref_directory = "./"
params.fasta  = "${ref_directory}/hg38.fa"
params.genecode_gtf = "${ref_directory}/gencode.v32.gtf"
params.threads = 12
threads = params.threads
println params.fasta

fasta = check_file(params.fasta)
gtf = check_file(params.genecode_gtf)

// bowtie2 index for s4_retrotranscriptome telescope 

process bowtie2_build {
    label 'prepare'
    label 'bigCPU'
    
    
    publishDir "./", mode: 'copy'
    input:
    file fasta

    output:
    file "bowtie2_index"

    """
    set +u; source activate pipeOne_telescope; set -u
    mkdir -p bowtie2_index
    bowtie2-build --threads ${threads}  ${fasta} bowtie2_index/bowtie2_base
    """
}
// bwa index for s2_circRNA, s6_RNAediting

process sprint_index {
    label 'prepare'
   

    publishDir "./", mode: 'copy',
        saveAs: {filename -> "$filename"}
    input:
    file "genome.fa" from fasta
    file "genome.gtf" from gtf 

    output:
    path "sprint_index/genome*"

    """
    set +u; source activate pipeOne_RnaEditing; set -u
    mkdir -p sprint_index; cd sprint_index
    cp ../genome.fa genome.fa
    sprint prepare -t ../genome.gtf genome.fa bwa
    """
    }

process star_index{
    label 'prepare'
    label 'bigMEM'

    publishDir "./", mode: 'copy'
    input:
    file "genome.fa" from fasta
    file "genome.gtf" from gtf 

    output:
    file "STAR_index"

    """
    set +u; source activate pipeOne_py3; set -u
    mkdir STAR_index; cd STAR_index
    STAR --runMode genomeGenerate --runThreadN ${threads} --genomeDir . --genomeFastaFiles ../genome.fa --sjdbGTFfile ../genome.gtf
    """
}

process hisat2_build {
    label 'prepare'
    label 'bigMEM'

    publishDir "./", mode: 'copy'
    input:
    file "genome.fa" from fasta
    file "genome.gtf" from gtf 

    output:
    file "hisat2_index"

    """
    set +u; source activate pipeOne_lncRNA; set -u
    # hisat2_extract_splice_sites.py genome.gtf >splice.txt
    # hisat2_extract_exons.py genome.gtf >exon.txt
    mkdir hisat2_index; cd hisat2_index
    # hisat2-build --exon ../exon.txt --ss ../splice.txt -p ${threads}  ../genome.fa genome
    hisat2-build  -p ${threads}  ../genome.fa genome
    """

}

process bwa_index{
    label 'prepare'
    publishDir "./", mode: 'copy'
    input:
    file fasta

    output:
    file "bwa_index"

    """
    set +u; source activate pipeOne_RnaEditing; set -u
    mkdir -p bwa_index
    bwa index -p bwa_index/bwa_index ${fasta}
    """
    
}


/*
process salmon_index{
    
    """
    ## apa
    qapa fasta -f $hg38_fa  qapa_3utrs.gencode_V31.hg38.bed qapa_3utrs.gencodeV31.hg38.fa
    python3 /dsk2/who/nbt/python/binPy3/fasta.py shorten_id qapa_3utrs.gencodeV31.hg38.fa qapa_3utrs.gencodeV31.hg38.short
    salmon index -t qapa_3utrs.gencodeV31.hg38.shortedID.fa -i 3utr -p 8

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