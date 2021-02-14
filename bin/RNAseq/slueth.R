#source("http://bioconductor.org/biocLite.R")
#biocLite("devtools")    # only if devtools not yet installed
#biocLite("pachterlab/sleuth")
#biocLite("rhdf5")
#setwd('/public/home/nong/pro/npc/lncRNA_2/results/kallisto')
options(stringsAsFactors = F)
suppressMessages({
  library("sleuth")
})
library(readr)
library(dplyr)
wd = 'samples'
datainfo = read_tsv('datainfo.txt2')
s2c = data.frame(sample=datainfo$Sample,condition = datainfo$Tissue, path =file.path(wd,datainfo$Sample) )
s2c = tbl_df(s2c)
print(s2c)
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

#summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05)
head(sleuth_significant, 20)

#plot an example DE transcript result
p1 = plot_bootstrap(so, "ENST00000328933", units = "est_counts", color_by = "condition")
p2 = plot_pca(so, color_by = 'condition')

#Print out the plots created above and store in a single PDF file
pdf(file="SleuthResults.pdf")
print(p1)
print(p2)
dev.off()
library("tximport")
browseVignettes("tximport")
