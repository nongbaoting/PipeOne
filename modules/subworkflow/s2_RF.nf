#!/usr/bin/env nextflow
/*  @Date         : 2021/02/02 12:59:21
    @Author       : Baoting Nong (nong55@foxmail.com)
    @Link         : https://github.com/nongbaoting
    @Version      : 
    @Description  : 
*/

nextflow.enable.dsl=2
// defined params
params.outdir_sub = params.outdir + "/s2_Features"
params.rawdir = ""
params.sample_info = ""
params.gene_info = ""
params.var_topK = 1000
params.test_size =  0.25
params.random_state = 2

// include funtions must be placed after params !
include { Header; check_file; input_reads; 
    get_base_index; get_dir_files } from "${baseDir}/modules/functions.nf"
include {proc; split_train_test_sample;
    RandomForest;add_gene_name } from "${baseDir}/modules/process/classification/RandomForest.nf"
//include {  } from  "${baseDir}/"
sample_info = check_file(params.sample_info, '--sample_info')

Channel
    .fromPath("${params.rawdir}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.rawdir}" }
    .set{tables}

workflow RF {
   
    main:
    proc(tables.collect(), sample_info)

    if(params.test_size ){

        split_train_test_sample( proc.out.collect(), sample_info )
        processed_data = split_train_test_sample.out.collect()
    }
    
    RandomForest( processed_data)
    if(params.gene_info){

        gene_info = check_file(params.gene_info, '--gene_info')
        add_gene_name(RandomForest.out.feature, gene_info)
    }
    

}