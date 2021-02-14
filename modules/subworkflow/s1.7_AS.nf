#!/usr/bin/env nextflow
/*  @Date         : 2021/01/31 17:27:06
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + '/s1.7_AS'
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?:false :false

// include funtions must be placed after params !
include { Header; check_file; input_reads;
     get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include {STAR_AS } from "${baseDir}/modules/process/mapping/STAR.nf"
include {spladder;  } from "${baseDir}/modules/process/AS/spladder.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"
//include {  } from  './modules/modules.nf'


fasta= check_file(params.fasta, '--fasta')
star_indices = get_dir_files(params.star_index, '--star_index')

workflow AS {
    take: 
    reads
    gtf

    main: 
    STAR_AS(reads, fasta, gtf, star_indices.collect() )
    spladder(STAR_AS.out.id_bam_bai, gtf )


    emit:
    id_bam = STAR_AS.out.id_bam
}