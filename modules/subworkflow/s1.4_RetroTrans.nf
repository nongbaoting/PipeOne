#!/usr/bin/env nextflow
/*  @Date         : 2021/01/31 15:30:36
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + '/s1.4_RetroTrans'
params.retro_gtf = params.genome ? params.genomes[ params.genome ].retro_gtf ?: false : false
params.bowtie2_index = params.genome ? params.genomes[ params.genome ].bowtie2_index ?:false :false


// include funtions must be placed after params !
include {Header; check_file; input_reads; get_base_index } from "${baseDir}/modules/functions.nf"
include {bowtie2_tele; merge_bowtie2_MapReads } from "${baseDir}/modules/process/mapping/bowtie2.nf"
include {telescope_bowtie2;  merge_telescope }  from "${baseDir}/modules/process/quant/telescope.nf"
include { mark_feature } from "${baseDir}/modules/modules.nf"
//include {  } from  './modules/modules.nf'


retro_gtf = check_file(params.retro_gtf, '--retro_gtf')
(bowtie2_base, bowtie2_indices) = get_base_index(params.bowtie2_index)

workflow RetroTrans {
    take: reads

    main:
    bowtie2_tele(reads, bowtie2_base, bowtie2_indices.collect() )
    merge_bowtie2_MapReads(bowtie2_tele.out.stat_count.collect() )
    telescope_bowtie2(bowtie2_tele.out.bam, retro_gtf )
    merge_telescope(retro_gtf,merge_bowtie2_MapReads.out, telescope_bowtie2.out.collect()  )
    mark_feature(merge_telescope.out.expr, 'retro-FPKM-divide_totalMapReads.csv', "Retro")
    
}