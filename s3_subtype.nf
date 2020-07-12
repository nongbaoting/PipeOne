#!/usr/bin/env nextflow


params.tables = ""

params.clinical = ""
params.threads = 4

clinical    = check_file(params.clinical)

Channel
    .fromPath("${params.tables}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.tables}" }
    .set{ tables }

process subtypes{

    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> "./${filename}"}

    input:
    file "00_rawdata/*" from tables.collect()
    file "sample.cli.csv" from clinical

    output:
    file "data/*" into pre_ch

    """
    set +u; source activate pipeOne_ml; set -u
    ## select topK variance features
    python3 ${baseDir}/bin/ML/python_code_2/proc_raw_data.py proc 00_rawdata/ sample.cli.csv
    
    ## defusion
    python3 ${baseDir}/bin/ML/python_code_2/run_defusion.py test_Params
    python3 ${baseDir}/bin/ML/python_code_2/check_convergence.py

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