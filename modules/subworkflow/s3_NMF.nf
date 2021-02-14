#!/usr/bin/env nextflow
/*  @Date         : 2021/02/02 15:25:19
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/
nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + "/s3_Subtype"
params.rawdir = ""
params.clinical = ""
params.var_topK = 1000
params.cluster_range = "3-8"

// include funtions must be placed after params !
include {Header; check_file; input_reads; get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
//include {  } from  "${baseDir}/"
include {NMF } from "${baseDir}/modules/process/subtype/NMF.nf"
clinical    = check_file(params.clinical, '--clinical')
Channel
    .fromPath("${params.rawdir}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.rawdir}" }
    .set{ tables }

def myTmpDir = file('./tmp')
def result = myTmpDir.mkdir()

workflow Subtype {

    NMF(tables.collect(), clinical)


}