library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
inFile = args[1]
outfile = args[2]
Supporting_reads_cutoff = as.numeric(args[3])


res = fread(inFile)

res %>% filter( supporting_reads >=Supporting_reads_cutoff, type  == "AG" | type  == "TC" ) %>% 
  select(chrom, start, end, sample, AD, DP) -> resFilter
resFilter %>% group_by(chrom, start, end, sample ) %>% 
  mutate(AD_sum = sum(AD), DP_sum = sum(DP) ) %>% 
  mutate(editding_rate = AD_sum / DP_sum) %>% ungroup() ->resFilter

resFilter %>% mutate(id = str_c(chrom, start, end, sep=':') ) %>%
  select(id,sample,editding_rate ) %>% unique(.) ->resFilter

resFilter %>%  spread(key = sample , value = editding_rate) -> resFilter_tab

write_tsv(resFilter_tab, outfile)