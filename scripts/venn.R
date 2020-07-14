#!/usr/bin/env Rscript

#upload library
#install.packages("VennDiagram")
#setwd("/public/home/nong/pro/npc/process2/results/plot_dir")
arg = commandArgs(trailingOnly = TRUE)
read_cutoff = arg[1]
sample_cutfoff = arg[2]

library(data.table)
ciri2 = fread('ciri2.table.txt', sep="\t")
ce2 = fread('circexplorer2.table.txt', sep="\t")
dcc = fread('dcc.table.txt', sep="\t") 

ciri2_f = ciri2[rowSums(ciri2[,-1] >= read_cutoff) >=sample_cutfoff,]
ce2_f   =  ce2[rowSums(ce2[,-1] >= read_cutoff) >=sample_cutfoff,   ]
dcc_f   =  dcc[rowSums(dcc[,-1] >= read_cutoff) >=sample_cutfoff,   ]

ciri2Name = ciri2_f$V1
ce2Name = ce2_f$V1
dccName = dcc_f$V1

allName = c(ciri2Name,ce2Name,dccName)
allName = unique(allName)
ciri2_in = allName %in% ciri2Name
ce2_in = allName %in% ce2Name
dcc_in = allName %in% dccName
all_in = ciri2_in + ce2_in + dcc_in
table(all_in)
write(allName[all_in>=2], 'at_least_2.circRNA')
write(allName[all_in>=3], 'at_least_3.circRNA')


