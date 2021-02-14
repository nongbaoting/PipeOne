library(data.table)
library(tidyverse)
reads_fi = "circRNA_bsj.csv"
libinfo_fi = "library_info.csv"
outfi = "circRNA_cpm.tsv"

reads = fread(reads_fi)
libinfo = read_csv(libinfo_fi)

libinfo %>% mutate(total_map = Mapped + Circular) -> libinfo
mysample = libinfo$Sample
reads_ord = reads[, ..mysample]

cpm = t( t(reads_ord)/ libinfo$total_map * 10^6 )
data.frame(circ_id = reads$circ_id, cpm, check.names = F) %>% as.data.table() ->cpm
write_tsv(cpm, path = outfi )

