#!/usr/bin/env nextflow
/*  @Date         : 2021/02/02 19:12:02
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.test = false

// include funtions must be placed after params !
include {Header; mkdir_tmp; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include {Subtype} from  "${baseDir}/modules/subworkflow/s3_NMF.nf"

log.info Header()
mkdir_tmp('./tmp')

workflow {

    Subtype()
    
}