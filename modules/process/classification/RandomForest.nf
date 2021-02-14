
process proc_and_split_train_test_sample {
    publishDir "./results",  mode: 'copy'
    stageInMode 'copy'

    input:
    path "00_rawdata/*" 
    path "s1_sample_info-tumor-normal.csv" 

    output:
    path "data/*" 

    """
    set +u; source activate pipeOne_ml; set -u
    python3  ${baseDir}/bin/RF/proc_raw_data.py proc --rawdir 00_rawdata  --sample_info s1_sample_info-tumor-normal.csv --var_topk ${params.var_topK}
    python3  ${baseDir}/bin/RF/proc_raw_data.py train_test_split --indir data/proc \
        --sample_info ./data/common_sample_info.csv \
        --test_size ${params.test_size}   --random_state ${params.random_state}
    """
}

process proc {
    publishDir "./results",  mode: 'copy'
    stageInMode 'copy'

    input:
    path "00_rawdata/*" 
    path "s1_sample_info-tumor-normal.csv" 

    output:
    path "data/*"  

    """
    set +u; source activate pipeOne_ml; set -u
    python3  ${baseDir}/bin/RF/proc_raw_data.py proc --rawdir 00_rawdata  \\
        --sample_info s1_sample_info-tumor-normal.csv \\
        --var_topk ${params.var_topK}
    """
}

process split_train_test_sample {
    publishDir "${params.outdir_sub}",  mode: 'copy'
    stageInMode 'copy'

    input:
    path "data/*" 
    path "s1_sample_info-tumor-normal.csv" 

    output:
    path "data/*" 

    """
    set +u; source activate pipeOne_ml; set -u
    python3  ${baseDir}/bin/RF/proc_raw_data.py train_test_split --indir data/proc \
        --sample_info ./data/common_sample_info.csv \
        --test_size ${params.test_size}   --random_state ${params.random_state}
    """
}

process RandomForest {
    stageInMode 'copy'
    label "bigCPU"

    publishDir "${params.outdir_sub}",  mode: 'copy',
        saveAs: {filename -> "./${filename}"}

    input:
    path "data/*" 

    output:
    path "data/feature_importance.csv" , emit: feature
    path "data/*"

    """
    set +u; source activate pipeOne_ml; set -u
    python3  ${baseDir}/bin/RF/main_randomForest.py --threads ${task.cpus}
    """

}

process add_gene_name{
    stageInMode 'copy'

    publishDir "${params.outdir_sub}",  mode: 'copy',
        saveAs: {filename -> "./${filename}" }

    input:
    path "data/feature_importance.csv"
    path "gene_info" 

    output:
    path "data/feature_importance-addName.csv"

    """
    set +u; source activate pipeOne_ml; set -u
    python3  ${baseDir}/bin/RF/result_summary.py feature data/feature_importance.csv  gene_info
    """

}



