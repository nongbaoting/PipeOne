library(tidyverse)
library(data.table)
## input
options(stringsAsFactors = F)
totalMapReads = fread("bowtie2_properPaired_total.txt")
names(totalMapReads) = c("Sample", "Reads")
dat = fread("telescope.rawCount.tsv")
txInfo = read_tsv("transcripts.info.tsv")


## order column
gcols =c("tx_id", totalMapReads$Sample)
bwDt = dat[ tx_id != "__no_feature", ..gcols]

## order splice length
txInfo_len = txInfo %>% select(tx_id, splice_len)
tibble(tx_id = bwDt$tx_id) %>% left_join(txInfo_len) ->bwDt_txLen


## 
bwFPKM_tmp = t( t(bwDt[,-1]) / totalMapReads$Reads ) * 10^9 # divide total counts and plus 10^9
bwFPKM = bwFPKM_tmp/bwDt_txLen$splice_len             # divide len
as_tibble(bwFPKM ) %>% add_column( tx_id = bwDt$tx_id, .before = 1) -> bwFPKM
write_csv(bwFPKM, "telescope.FPKM-divide_totalMapReads.csv")


## fpkm use total retrotranscriptome reads as total reads
totalRetroReads = colSums(bwDt[,-1])
reFPKM_tmp = t( t(bwDt[,-1]) / totalRetroReads ) * 10^9 # divide total counts and plus 10^9
reFPKM = reFPKM_tmp / bwDt_txLen$splice_len                 # divide len
as_tibble(reFPKM ) %>% add_column( tx_id = bwDt$tx_id, .before = 1) -> reFPKM
write_csv(reFPKM, "telescope.FPKM-divide_totalRetroReads.csv")

## tpm
RPK = bwDt[,-1]/bwDt_txLen$splice_len * 1000 
perM = colSums(RPK) / (10^6)
tpm = t( t(RPK) /perM )
as_tibble(tpm ) %>% add_column( tx_id = bwDt$tx_id, .before = 1) -> tpm
write_csv(tpm , "telescope.TPM.csv")

