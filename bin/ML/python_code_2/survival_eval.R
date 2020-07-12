library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)
alpha <- c(0.01, 0.02, 0.1, 1, 5)
gamma <- c(0, 0.1, 1, 10, 100)
low_dim <- c(2, 3, 4, 5, 6, 7)

res_dir <- './results/'

# read data
#
record_file <- 'record_log_rank_test_pvalue.txt'
reocrd_file_csv = "record_log_rank_test_pvalue.csv"
re_df = data.frame()
if (file.exists(record_file)){
    file.remove(record_file)
}

if (!file.exists("./data/surv_curve/")){
  dir.create("./data/surv_curve")
}

fout <- file(record_file, open = "a")
for (low_dim_ in low_dim){
  
  for (alpha_ in alpha){
    for (gamma_ in gamma){
      clinical_path <- paste0(res_dir, 
                              sprintf('lowDim=%d/lowDim=%d_alpha=%.2f_gamma=%.2f_clustering.csv',
                                      low_dim_, low_dim_, alpha_, gamma_))
      
      clinical <- read.csv(clinical_path,header = T)
      fit <- coxph(Surv(time_to_event, event)~subtype, data=clinical)
      summary_fit <- summary(fit)
      pvalue <- summary_fit$logtest[3]
      outline <- sprintf("In low_dim=%d, alpha=%.2f_gamma=%.2f, log rank test p=%f",
                         low_dim_, alpha_, gamma_, pvalue)
      a = data.frame(low_dim = c(low_dim_), alpha = c(alpha_), gamma = c( gamma_) , pvalue = c(pvalue) )
      re_df = rbind(re_df, a)
      print(outline)
      write(outline, fout, append = TRUE)
      
      surv_fit <- survfit(Surv(time_to_event, event)~subtype, data=clinical)
      res <- ggsurvplot(surv_fit,
                        pval = T, conf.int = F, risk.table = T,
                        tables.theme = theme_cleantable())
      res$table <- res$table + theme(axis.line = element_blank())
      # res$plot <- res$plot + labs(title = "Survival Curves")

      curve_path <- paste0('./data/surv_curve/',
                           sprintf('low_dim=%d_alpha=%.2f_gamma=%.2f_clustering.png',
                                   low_dim_, alpha_, gamma_))
      ggsave(curve_path, print(res))
	  rm(a)
	  rm(summary_fit)
	  rm(pvalue)
    }
    }
    
  
}

close(fout)
re_df %>% arrange(pvalue) ->re_df

write_csv(re_df, reocrd_file_csv)
