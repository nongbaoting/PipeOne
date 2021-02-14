

process defusion {
    label 'bigCPU'
    publishDir "${params.outdir_sub}",  mode: 'copy'
   
    input:
    path "00_rawdata/*" 
    path "sample.cli.csv" 

    output:
    
    path "*" 

    """
    set +u; source activate pipeOne_ml; set -u
    ## select topK variance features
    python3  ${baseDir}/bin/NMF/proc_raw_data.py proc --rawdir 00_rawdata/ --sample_info sample.cli.csv --var_topk ${params.var_topK} 
    
    ## defusion
    python3  ${baseDir}/bin/NMF/run_defusion.py --threads ${task.cpus}  --test ${params.test}
    python3  ${baseDir}/bin/NMF/check_convergence.py
    """
}

process clustering_and_eval{
    stageInMode 'copy'
    publishDir "${params.outdir_sub}",  mode: 'copy'

    input:
    path "*" 
    path "data/sample.cli.csv"

    output:
    tuple path("00_rawdata"),  path("clusters"),   path("data"), 
        path("NMF"),  path("record_log_rank_test_pvalue.csv"),  
        path("record_log_rank_test_pvalue.log.txt") 

    """
    set +u; source activate pipeOne_ml; set -u
    ## clustering and eval
    python3  ${baseDir}/bin/NMF/eval_cluster_num.py --cluster_range ${params.cluster_range }
    Rscript  ${baseDir}/bin/NMF/survival_eval.R ./data/sample.cli.csv ./clusters/surv_curve/ ${params.cluster_range }
    """
}


process features_selection {
    stageInMode 'copy'
    publishDir "${params.outdir_sub}", mode: 'copy'

    input:
    path "*" 
   
    output:
    path "*" 

    """
    set +u; source activate pipeOne_ml; set -u
    ## select features
    python3  ${baseDir}/bin/NMF/select_topk_nong.py
    python3  ${baseDir}/bin/NMF/find_best_RFparams.py 
    """
}


workflow NMF{

    take:
    tables
    clinical

    main:
    defusion(tables.collect(), clinical)
    clustering_and_eval( defusion.out.collect(), clinical )
    features_selection(clustering_and_eval.out.collect() )

    
}