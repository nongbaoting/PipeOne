#!/usr/bin/env nextflow
/*  @Date         : 2021/02/01 00:01:44
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + '/s1.2_circRNA'
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.hisat2_idx = params.genome ? params.genomes[ params.genome ].hisat2_index ?:false :false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?: false : false

// include funtions must be placed after params !
include { hisat2; } from  "${baseDir}/modules/process/mapping/hisat2.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"
include { Header; check_file; input_reads; get_base_index;get_dir_files } from "${baseDir}/modules/functions.nf"
include {CIRIquant } from "${baseDir}/modules/process/circRNA/CIRIquant.nf"
//include {  } from  './modules/modules.nf'

def (hisat2_base, hisat2_indices) = get_base_index(params.hisat2_idx )
def (bwa_base, bwa_indices) = get_base_index( params.bwa_index )
//b2 = get_base_index( params.bwa_index )

def fasta    = check_file(params.fasta, '--fasta')

//println(b2)

workflow CircRNA_CIRIquant_fastq {
    take:
    reads
    gtf

    main:
    hisat2(reads,hisat2_base, hisat2_indices )
    CircRNA_CIRIquant_bam(reads, hisat2.out.bam, gtf )
    

}

workflow CircRNA_CIRIquant_bam {
    take:
    reads
    bam
    gtf

    main:
    CIRIquant(
            reads, 
            bam,
            fasta, 
            gtf,
            bwa_base, bwa_indices,
            hisat2_base, hisat2_indices
    )

    mark_feature(CIRIquant.out, "circRNA_CPM.csv", "circRNA")

}