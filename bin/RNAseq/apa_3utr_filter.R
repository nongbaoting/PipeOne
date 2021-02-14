
options(stringsAsFactors = F)
arg = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyverse, warn.conflicts = FALSE))

## input
apa = suppressMessages(read_tsv(arg[1]))
TPM = suppressWarnings(fread(arg[2]))
allInfo = fread( arg[3] )
names(TPM)[1] =  "gene_id"

## filter TPM >3 in 50% of samples
tpm4apa = TPM[ rowSums(TPM[,-1] >3) >= ( (dim(TPM)[2] -1 )* 0.5 ), ]
info4apa = allInfo %>% filter(gene_id %in% tpm4apa$gene_id )

apa %>% filter(Length <= 100 ) -> apa.exclude

## select tpm > 3 genes, apa length >100, no NA value
apa %>% filter(Gene_Name %in% info4apa$gene_name, ! Gene_Name %in% apa.exclude$Gene_Name ) %>% 
  filter(!grepl("^NA", APA_ID)) ->apa.pre

apa.pre %>% filter(!grepl("S$", APA_ID) ) ->apa.preDP

## write table
apa.pre %>% select(APA_ID, ends_with( '.PAU') ) ->apa.pre.pau
apa.preDP  %>% select(APA_ID, ends_with( '.PAU') ) ->apa.preDP.pau
names(apa.pre.pau) %>% str_replace(".PAU$", '') ->names(apa.pre.pau) 
names(apa.preDP.pau) %>% str_replace(".PAU$", '') ->names(apa.preDP.pau) 

write_tsv(apa.pre.pau, "./pau_results.filterPau.txt")
write_tsv(apa.preDP.pau, "./pau_results.filterPau-distal-proximal.txt")

