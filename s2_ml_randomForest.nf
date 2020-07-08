#!/usr/bin/env nextflow

params.tables = ""
params.sample_info = ""
params.threads = 4

if ( params.sample_info ){
    sample_info = file(params.sample_info)
    if( !sample_info.exists() ) exit 1, "Sample file not found: ${params.sample_info}"
}else{exit 1, "No  sample file specified!"}

tables = Channel
        .fromPath("${params.tables}/*.csv")
        .ifEmpty { exit 1, "files not found: ${params.tables}" }

process proc_and_split_train_test_sample{
    input:
    file "s1_sample_info-tumor-normal.csv" from sample_info
    file "00_rawdata/*" from tables.collect()

    output:
    file "data/*" into pre_ch

    """
    set +u; source activate pipeone; set -u
    python3 ${baseDir}/ML/proc_raw_data.py proc --rawdir 00_rawdata  --sample_want s1_sample_info-tumor-normal.csv --var_topk 1000
    python3 ${baseDir}/ML/proc_raw_data.py train_test_split data/proc s1_sample_info-tumor-normal.csv
    """

}

process RandomForest{
    publishDir "./results",  mode: 'copy',
        saveAs: {filename -> "./${filename}"}

    input:
    file "data/*" from pre_ch.collect()

    output:
    tuple "data/feature_importance-addName.csv","data/feature_importance.csv", "data/discriminative_power_of_topk_feature.csv"

    """
    set +u; source activate pipeone; set -u
    python3 ${baseDir}/ML/main_randomForest.py --threads ${params.threads }
    python3 ${baseDir}/ML/result_summary.py feature data/feature_importance.csv
    """

}

