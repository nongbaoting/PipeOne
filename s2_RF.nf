#!/usr/bin/env nextflow
/*  @Date         : 2021/02/02 16:26:07
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params

// include funtions must be placed after params !
include {Header; mkdir_tmp; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include { RF } from  "${baseDir}/modules/subworkflow/s2_RF.nf"

log.info Header()
mkdir_tmp('./tmp')

workflow {

    RF()

}