library(readr)
library(dplyr)
lncinfo = read_delim('all_lncRNA_info.tsv',delim ='\t')
tab = read.table('kallisto_tx_tpm.tsv')
tab = tab[row.names(tab) %in% lncinfo$tx_id,]
keep = rowSums(tab >10) >=12
table(keep)
tab_allE=tab[keep,]

tab_sum = data.frame(tab, sum = apply(tab, 1, sum ))

tab_sum =  tab_sum[order(tab_sum$sum,decreasing = T),]

keep = rowSums(tab_sum >10) >=13

tabb = tab_sum[keep,]
top_tpm = tabb[1:100,]
write.table(top_tpm, 'highest.lncRNA.tpm.xls',sep='\t',quote = F)



deg_txid_gene = read_delim('deg_list.tsv',delim = '\t')
names(deg_txid_gene) = c('tx_id','gene_id')
al = left_join(deg_txid_gene,lncinfo, by=c("tx_id",'gene_id'))
tx_tpm = data.frame(tx_id = row.names(tab),tab)
tx_tpm = tbl_df(tx_tpm)
al = left_join(al, tx_tpm)
write_csv(al,'deg_info_list_tpm.csv')


#top_100 = read_delim('top_expression_new_lncRNA.xls',delim = '\t')
#names(top_100) = c('tx_id','gene_id')
#top = left_join(top_100, lncinfo)
#top = left_join(top_100, tx_tpm)
#write.table(top, 'top_gene_info_new_lncRNA_tpm.xls', sep='\t')


fc = read.table('merged_gene_counts.txt')
fc = fc[row.names(fc) %in% lncinfo$gene_id,]
keep = rowSums(fc >10) >=12
table(keep)
fc_allE=fc[keep,]
fc_sum  = data.frame(fc, sum = apply(fc, 1, sum ))
fc_sum = fc_sum[order(fc_sum$sum,decreasing = T),]
keep = rowSums(fc_sum >10) >=13
table(keep)
fc_tabb = fc_sum[keep,]
top_fc = fc_tabb [1:100,]
write.table(top_fc, 'highest.lncRNA.feacturecount.xls',sep='\t',quote = F)

top_fc_info = lncinfo %>% filter( gene_id %in% row.names(top_fc) ) %>% filter(tx_id %in% row.names(top_tpm))

top_new_lncRNA = left_join(top_fc_info,tx_tpm)
write_excel_csv(top_new_lncRNA,"top_new_lncRNA.csv.tmp")
system("awk 'NR==1 || /^TU/' top_new_lncRNA.csv.tmp >top_new_lncRNA.csv")

