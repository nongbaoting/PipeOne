#!/usr/bin/env nextflow

params.rawdir = ""
params.var_topK = 1000
params.clinical = ""
params.threads = 24
params.cluster_range = "3-8"

clinical    = check_file(params.clinical)

Channel
    .fromPath("${params.rawdir}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.rawdir}" }
    .set{ tables }

process defusion {
    publishDir "./results",  mode: 'copy'
   
    input:
    file "00_rawdata/*" from tables.collect()
    file "sample.cli.csv" from clinical

    output:
    
    file "*" into next_step

    """
    set +u; source activate pipeOne_ml; set -u
    ## select topK variance features
    python3 ${baseDir}/bin/ML/python_code_2/proc_raw_data.py proc --rawdir 00_rawdata/ --sample_info sample.cli.csv --var_topk ${params.var_topK} 
    
    ## defusion
    python3 ${baseDir}/bin/ML/python_code_2/run_defusion.py --threads ${params.threads}  
    python3 ${baseDir}/bin/ML/python_code_2/check_convergence.py
    """
}

process clustering_and_eval{
    stageInMode 'copy'
    publishDir "./results",  mode: 'copy'

    input:
    file "*" from next_step.collect()
    file "data/sample.cli.csv" from clinical

    output:
    tuple "00_rawdata", "clusters",  "data" , "NMF", "record_log_rank_test_pvalue.csv", "record_log_rank_test_pvalue.log.txt" into cluster_res

    """
    set +u; source activate pipeOne_ml; set -u
    ## clustering and eval
    python3 ${baseDir}/bin/ML/python_code_2/eval_cluster_num.py --cluster_range ${params.cluster_range }
    Rscript ${baseDir}/bin/ML/python_code_2/survival_eval.R ./data/sample.cli.csv ./clusters/surv_curve/ ${params.cluster_range }
    """
}


process features_selection{
    stageInMode 'copy'
    publishDir "./results", mode: 'copy'

    input:
    file "*" from cluster_res.collect()
   
    output:
    file "*" into pre_ch

    """
    set +u; source activate pipeOne_ml; set -u
    ## select features
    python3 ${baseDir}/bin/ML/python_code_2/select_topk_nong.py
    python3 ${baseDir}/bin/ML/python_code_2/find_best_RFparams.py 
    """
}


def check_file(myfile){
    if ( myfile ){
        myfile_ = file(myfile)
        if( !myfile_.exists() ) exit 1, "File not found: ${myfile}"

    }else{exit 1, "value is null !"}

    return myfile_
}