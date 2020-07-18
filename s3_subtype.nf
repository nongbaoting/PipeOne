#!/usr/bin/env nextflow

params.rawdir = ""
params.var_topK = 1000
params.clinical = ""
params.threads = 4

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
    python3 ${baseDir}/bin/ML/python_code_2/run_defusion.py test_Params ${params.threads}
    python3 ${baseDir}/bin/ML/python_code_2/check_convergence.py

    
    """

}

process clustering_and_eval{
    stageInMode 'copy'
    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> if(filename =~ /data/ ) "./${filename}"}


    input:
    file "*" from next_step.collect()
    file "sample.cli.csv" from clinical

    output:
    file "data/*" into pre_ch

    """
    set +u; source activate pipeOne_ml; set -u
    ## clustering and eval
    python3 ${baseDir}/bin/ML/python_code_2/consistency_eval.py 3
    Rscript ${baseDir}/bin/ML/python_code_2/survival_eval.R
    
    ## select features
    python3 ${baseDir}/bin/ML/python_code_2/select_topk_nong.py select_topK
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