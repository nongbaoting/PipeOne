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

write.table(txi$counts, 'salmon_gene_est_counts.tsv',sep = '\t', quote = F)
write.table(txi$abundance, 'salmon_gene_tpm.tsv', sep="\t", quote = F )

txi_tx <- tximport(kal_files, type =  "salmon", txOut = TRUE)
write.table(txi_tx$counts, 'salmon_tx_est_counts.tsv',sep = '\t', quote = F )
write.table(txi_tx$abundance, 'salmon_tx_tpm.tsv', sep="\t", quote = F )

# tx_abundance = txi_tx$abundance
# tx_keep =  rowSums(tx_abundance >= 0.01 ) >=2
# tx_abundance = tx_abundance[tx_keep, ]
# tx_id = rownames(tx_abundance)
# write(tx_id, "tx_id.keep.txt")

gene_abundance = txi$abundance
gene_keep =  rowSums(gene_abundance  >= 1 ) >= 2
gene_abundance = gene_abundance[gene_keep, ]
gene_id = rownames(gene_abundance)
write(gene_id, "gene_id.keep.txt")

tx2gene[tx2gene$geneID %in% gene_id, ] -> tx2gene_keep

tx_abundance = txi_tx$abundance
tx_tpm_keep   = tx_abundance[ rownames(tx_abundance) %in% tx2gene_keep$txID,  ]
tx_count_keep = txi_tx$counts[rownames(txi_tx$counts) %in% tx2gene_keep$txID, ]

tpm_keep   = txi$abundance[rownames(txi$abundance) %in% tx2gene_keep$geneID, ]
count_keep = txi$counts[rownames(txi$counts) %in% tx2gene_keep$geneID, ]

write.table(count_keep, 'salmon_gene_est_counts.keep.tsv',sep = '\t', quote = F)
write.table(tpm_keep, 'salmon_gene_tpm.keep.tsv', sep="\t", quote = F )

write.table(tx_count_keep, 'salmon_tx_est_counts.keep.tsv',sep = '\t', quote = F )
write.table(tx_tpm_keep, 'salmon_tx_tpm.keep.tsv', sep="\t", quote = F )