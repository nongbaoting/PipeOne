#!/usr/bin/env nextflow

params.rawdir = ""
params.sample_info = ""
params.threads = 4
params.gene_info = ""
params.var_topK = 1000
params.test_size = 0.3
params.random_state = 2



sample_info = check_file(params.sample_info)

Channel
    .fromPath("${params.rawdir}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.rawdir}" }
    .set{tables}


process proc_and_split_train_test_sample {
    publishDir "./results",  mode: 'copy'

    input:
    file "s1_sample_info-tumor-normal.csv" from sample_info
    file "00_rawdata/*" from tables.collect()

    output:
    file "data/*" into pre_ch

    """
    set +u; source activate pipeOne_ml; set -u
    python3 ${baseDir}/bin/ML/proc_raw_data.py proc --rawdir 00_rawdata  --sample_info s1_sample_info-tumor-normal.csv --var_topk ${params.var_topK}
    python3 ${baseDir}/bin/ML/proc_raw_data.py train_test_split --indir data/proc \\
        --sample_info s1_sample_info-tumor-normal.csv \\
        --test_size ${params.test_size}   --random_state ${params.random_state}
    """
}

process RandomForest{
    
    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> "./${filename}"}

    input:
    file "data/*" from pre_ch.collect()

    output:
    file "data/feature_importance.csv" into rf_res
    file "data/*"

    """
    set +u; source activate pipeOne_ml; set -u
    python3 ${baseDir}/bin/ML/main_randomForest.py --threads ${params.threads }
    """

}

if( params.gene_info){
     
    gene_info = check_file(params.gene_info )

    process add_gene_name{
    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> "./${filename}" }

    input:
    file "data/feature_importance.csv" from  rf_res
    file "gene_info" from  gene_info
    
    output:
    file "data/feature_importance-addName.csv"
    
    """
    set +u; source activate pipeOne_ml; set -u
    python3 ${baseDir}/bin/ML/result_summary.py feature data/feature_importance.csv  gene_info
    """

    }

}



def check_file(myfile){
    if ( myfile ){
        myfile_ = file(myfile)
        if( !myfile_.exists() ) exit 1, "File not found: ${myfile}"

    }else{exit 1, "value is null !"}

    return myfile_
}