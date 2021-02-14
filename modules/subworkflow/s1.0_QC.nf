#!/usr/bin/env nextflow
/*  @Date         : 2021/02/01 18:53:08
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir  + "/s1.0_QC"

// include funtions must be placed after params !
include {Header; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include {fastp; } from "${baseDir}/modules/process/QC/fastp.nf"


workflow QC {
    take:
    ch_reads

    main:
    if( ! params.cleaned ){
        fastp(ch_reads)
        reads = fastp.out.reads
    }else{
        reads = ch_reads
    }

    emit:
    reads = reads

}