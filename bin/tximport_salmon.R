#!/usr/bin/env Rscript
#http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon-sailfish
library("tximport")
library(readr)
wd = 'samples'
tx2gene = read_tsv("protein_coding_and_all_lncRNA.txID_geneID.tsv")
sample_id <- dir(file.path(wd))
sample_id
kal_files = file.path(wd, sample_id,   "quant.sf")
names(kal_files) = sample_id
txi <- tximport(kal_files  , type =  "salmon", tx2gene = tx2gene)

write.table(txi$counts, 'salmon_gene_est_counts.tsv',sep = '\t',quote = F)
write.table(txi$abundance, 'salmon_gene_tpm.tsv', sep="\t",quote = F )

txi_tx <- tximport(kal_files, type =  "salmon", txOut = TRUE)
write.table(txi_tx$counts, 'salmon_tx_est_counts.tsv',sep = '\t',quote = F )
write.table(txi_tx$abundance, 'salmon_tx_tpm.tsv', sep="\t",quote = F )


