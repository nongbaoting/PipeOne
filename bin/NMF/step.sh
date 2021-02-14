baseDir=/dat1/nbt2/pipe/PipeOne
set +u; source activate pipeOne_ml; set -u
## select topK variance features
#python3 ${baseDir}/bin/ML/python_code_2/proc_raw_data.py proc --rawdir 00_rawdata/ --sample_info sample.cli.csv --var_topk 1000

## defusion
#python3 ${baseDir}/bin/ML/python_code_2/run_defusion.py --threads 95
#python3 ${baseDir}/bin/ML/python_code_2/check_convergence.py

## clustering and eval

python3 ${baseDir}/bin/ML/python_code_2/eval_cluster_num.py --cluster_range 3-8
Rscript ${baseDir}/bin/ML/python_code_2/survival_eval.R ./data/sample.cli.csv ./clusters/surv_curve/ 3-8

   


## select features
python3 ${baseDir}/bin/ML/python_code_2/select_topk_nong.py
python3 ${baseDir}/bin/ML/python_code_2/find_best_RFparams.py 