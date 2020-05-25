

library(data.table)
library(tidyverse)


## input

apa = read_tsv("../s3_APA-3TUR/results/apa_3utr/pau_results.txt")
TPM = fread("../s1_lncRNA/results/salmon/salmon_gene_tpm.tsv")
allInfo = fread("../s1_lncRNA/results/novel_lncRNA/protein_coding_and_all_lncRNA.info.tsv")
names(TPM)[1] =  "gene_id"


## filter TPM >3 in 50% of samples
tpm4apa = TPM[ rowSums(TPM[,-1] >3) >= ( (dim(TPM)[2] -1 )* 0.9 ), ]
info4apa = allInfo %>% filter(gene_id %in% tpm4apa$gene_id )

apa %>% filter(Length <= 100 ) -> apa.exclude

## select tpm > 3 genes, apa length >100, no NA value
apa %>% filter(Gene_Name %in% info4apa$gene_name, ! Gene_Name %in% apa.exclude$Gene_Name ) %>% 
  filter(!grepl("^NA", APA_ID)) ->apa.pre

apa.pre %>% filter(!grepl("S$", APA_ID) ) ->apa.preDP


## write table

apa.pre %>% select(APA_ID, ends_with( '.PAU') ) ->apa.pre.pau
apa.preDP  %>% select(APA_ID, ends_with( '.PAU') ) ->apa.preDP.pau
write_tsv(apa.pre.pau, "../s3_APA-3TUR/results/apa_3utr/pau_results.filterPau.txt")
write_tsv(apa.preDP.pau, "../s3_APA-3TUR/results/apa_3utr/pau_results.filterPau-distal-proximal.txt")

system("mkdir -p s3_APA-3TUR")
write_tsv(apa.preDP.pau, "s3_APA-3TUR/pau_results.filterPau-distal-proximal.txt")

