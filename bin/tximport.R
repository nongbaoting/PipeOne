#!/usr/bin/env Rscript

library("tximport")
library(readr)
wd = 'samples'
tx2gene = read_tsv("protein_coding_and_all_lncRNA.txID_geneID.tsv")
sample_id <- dir(file.path(wd))
sample_id
kal_files = file.path(wd, sample_id,   "abundance.h5")
names(kal_files) = sample_id
txi <- tximport(kal_files  , type =  "kallisto", tx2gene = tx2gene)

write.table(txi$counts, 'kallisto_gene_est_counts.tsv',sep = '\t',quote = F)
write.table(txi$abundance, 'kallisto_gene_tpm.tsv', sep="\t",quote = F )

txi_tx <- tximport(kal_files, type =  "kallisto", txOut = TRUE)
write.table(txi_tx$counts, 'kallisto_tx_est_counts.tsv',sep = '\t',quote = F )
write.table(txi_tx$abundance, 'kallisto_tx_tpm.tsv', sep="\t",quote = F )


