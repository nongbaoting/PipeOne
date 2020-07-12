#!/usr/bin/env nextflow

params.tables = ""
params.sample_info = ""
params.threads = 4
params.gene_info = ""

if ( params.sample_info ){
    sample_info = file(params.sample_info)
    if( !sample_info.exists() ) exit 1, "Sample file not found: ${params.sample_info}"
}else{exit 1, "No  sample file specified!"}

Channel
    .fromPath("${params.tables}/*.csv")
    .ifEmpty { exit 1, "files not found: ${params.tables}" }
    .set{tables}

process proc_and_split_train_test_sample{

    input:
    file "s1_sample_info-tumor-normal.csv" from sample_info
    file "00_rawdata/*" from tables.collect()

    output:
    file "data/*" into pre_ch

    """
    set +u; source activate pipeOne_ml; set -u
    python3 ${baseDir}/bin/ML/proc_raw_data.py proc --rawdir 00_rawdata  --sample_want s1_sample_info-tumor-normal.csv --var_topk 1000
    python3 ${baseDir}/bin/ML/proc_raw_data.py train_test_split data/proc s1_sample_info-tumor-normal.csv
    """
}

process RandomForest{
    
    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> "./${filename}"}

    input:
    file "data/*" from pre_ch.collect()

    output:
    tuple "data/feature_importance.csv", "data/discriminative_power_of_topk_feature.csv" into rf_res

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
    tuple "data/feature_importance.csv", "data/discriminative_power_of_topk_feature.csv" from  rf_res
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