
library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)
argv = commandArgs(trailingOnly=TRUE)
cli_fi = argv[1]
outdir = argv[2]
cluster_range = argv[3]

alpha <- c(0.01, 0.02, 0.1, 1, 5)
gamma <- c(0, 0.1, 1, 10, 100)
low_dim <- c(2, 3, 4, 5, 6, 7)
# cluster_num = c(2, 3, 4, 5, 6, 7, 8)
cr = as.numeric(str_split(cluster_range, '-', simplify=T) )
cluster_num = seq(cr[1], cr[2] )

res_dir = "clusters/eval_cluster_num"
# read data
#
record_file <- 'record_log_rank_test_pvalue.log.txt'
record_file_csv = "record_log_rank_test_pvalue.csv"
re_df = data.frame()
if (file.exists(record_file)){
    file.remove(record_file)
}

outdir = paste0(outdir, '/')
if (!file.exists(outdir)){
  dir.create(outdir)
}
cli = read_csv(cli_fi)
# if (!file.exists("./clusters/surv_curve/")){
#   dir.create("./clusters/surv_curve/")
# }

# cli_fi = "KIRP_cli.OS.csv"
# cli = read_csv(cli_fi)

fout <- file(record_file, open = "a")
for (low_dim_ in low_dim){
  
  for (alpha_ in alpha){
    for (gamma_ in gamma){
      for( cluster_num_ in cluster_num){
        #clinical_path <- paste0(res_dir, 
        #                      sprintf('lowDim=%d/lowDim=%d_alpha=%.2f_gamma=%.2f_clusters=%d_clustering.csv',
        #                              low_dim_, low_dim_, alpha_, gamma_, cluster_num_))
        clinical_path <- paste0(res_dir, 
                              sprintf('/lowDim=%d_alpha=%.2f_gamma=%.2f_clusters=%d_clustering.csv',
                                      low_dim_, alpha_, gamma_, cluster_num_ )  )
        if(! file.exists(clinical_path ) ) next                    
        clinical <- read.csv(clinical_path, header = T)
        clinical %>% left_join(cli) -> clinical

        surv_fit <- survfit(Surv(Time_to_event, Event)~Subtype, data=clinical)
        surv.diff <- survdiff(Surv(Time_to_event, Event)~Subtype, data=clinical)
        pvalue <- 1- pchisq (surv.diff$chisq[1], df= length(unique(clinical$Subtype)) -1 )
        pvalue = signif(pvalue, 2)

        outline <- sprintf("In lowDim=%d, alpha=%.2f_gamma=%.2f, clusters=%d, log rank test p=%f",
                          low_dim_, alpha_, gamma_, cluster_num_, pvalue)
        #a = data.frame(low_dim = c(low_dim_), alpha = c(alpha_), gamma = c( gamma_) , pvalue = c(pvalue) )
        a = data.frame(NMF_param = c(sprintf("lowDim=%d_alpha=%.2f_gamma=%.2f",
                          low_dim_, alpha_, gamma_ ) ) , n_clusters = c(cluster_num_) ,logRankTest_pvalue = c(pvalue) )
        re_df = rbind(re_df, a)
        print(outline)
        write(outline, fout, append = TRUE)
        
        surv_fit <- survfit(Surv(Time_to_event, Event)~Subtype, data=clinical)
        res <- ggsurvplot(surv_fit,
                          pval = T, conf.int = F, risk.table = T,
                          tables.theme = theme_cleantable())
        res$table <- res$table + theme(axis.line = element_blank())
        # res$plot <- res$plot + labs(title = "Survival Curves")

        #curve_path <- paste0('./data/surv_curve/',sprintf('low_dim=%d_alpha=%.2f_gamma=%.2f_clustering.png',low_dim_, alpha_, gamma_))
        #ggsave(curve_path, print(res))
        curve_path <- paste0(outdir, sprintf('low_dim=%d_alpha=%.2f_gamma=%.2f__clusters=%d_clustering.pdf',low_dim_, alpha_, gamma_, cluster_num_) )
        pdf(file = curve_path )
        print(res)
        dev.off()

        rm(a)
        rm(summary_fit)
        rm(pvalue)
        rm(clinical)
      }
      
    }
    }
    
  
}

close(fout)



silhoutte = read_csv("./clusters/eval_cluster_num/silhoutte_score_summary_2.csv")
as_tibble(re_df) %>% mutate(significant = logRankTest_pvalue <0.05 ) %>%
 left_join(silhoutte) %>% arrange(logRankTest_pvalue)  -> re_df
 
write_csv(re_df, record_file_csv)
