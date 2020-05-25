#!/usr/bin/env Rscript
options(stringsAsFactors = F)
arg = commandArgs(trailingOnly=TRUE)
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library(DESeq2)
library(dplyr)
library(pheatmap)
library(readr)
#setwd("/public/home/nong/pro/npc/lncRNA_2/results/dat_bak/deg_deseq2/feature_counts")
tab = arg[1]
#tab = "kallisto_gene_est_counts.tsv"
#tab = 'merged_gene_counts.txt'
tx_gene = read_tsv("protein_coding_and_all_lncRNA.txID_geneID.tsv")
lnc_tx = read.table('all_lncRNA.list',stringsAsFactors = F)
lnc_df = filter(tx_gene, txID %in% lnc_tx$V1 )

data =  read.table(tab)
dat_df =data.frame(id = row.names(data),data)
dat_df = tbl_df(dat_df)
datainfo = read.table('datainfo.txt',header = T)
rownames(datainfo) = datainfo$Sample
datainfo = with(datainfo, datainfo[order(Patient,Tissue),])

kp = row.names(data) %in% lnc_df$geneID
data = data[kp,]

keep =  rowSums(data >=5 ) >=2
data_f = data[keep,]

#dat_df %>% filter(id %in% ids$V1) %>% select(c("id", as.vector(datainfo$Sample) )) -> dat_remain

counts =data_f[, as.vector(datainfo$Sample)]
counts = round(counts)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = datainfo,
                              design = ~ Patient + Tissue)
dds <- DESeq(dds)
res <- results(dds)

res_da = as.data.frame(res)
write.csv(res_da, 'deseq2.csv')
res_da$id = rownames(res_da)
res_df = tbl_df(res_da)
deseq2_sig = res_df %>% filter( padj < 0.05, abs(log2FoldChange) >1, pvalue <0.01)


vsd <- vst(dds, blind = FALSE)
mat = assay(vsd)
mat = mat[row.names(mat) %in% deseq2_sig$id,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[,  c("Patient", "Tissue")])

pheatmap(mat, annotation_col = anno,show_rownames = F,filename= 'DEG_heatmap.tiff')
counts_sig = dat_df %>% filter(id %in% deseq2_sig$id)
deg_txid_gene = tx_gene %>% filter(geneID %in% deseq2_sig$id)
write_tsv(deg_txid_gene, 'deg_list.tsv')

count_sum =data.frame(counts, sum = rowSums(counts))
count_sum[order(count_sum$sum, decreasing = T),] ->count_sum
top_txid_gene = tx_gene %>% filter(geneID %in% rownames(count_sum[1:100,]) )
write_tsv(top_txid_gene, 'top_100_express_lncRNA_genes.txt')


