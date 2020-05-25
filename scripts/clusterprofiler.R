#https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Dr.eg.db")
#biocLite("org.Mm.eg.db")
#biocLite("org.Hs.eg.db")
#biocLite("clusterProfiler")
#biocLite("Rgraphviz")
#biocLite("topGO")
#biocLite("DOSE")
#biocLite("pathview")
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DOSE)
library(Rgraphviz)
library(tools)

Argument = commandArgs(TRUE)
gene_list= Argument[1]
wd = Argument[2]
organism = Argument[3]
cmd = paste('mkdir -p ',wd);system(cmd)
wd = file_path_as_absolute(wd)


pvalueCutoff=0.05
qvalueCutoff=0.1
padjustCutoff = 0.1

x=read.table(gene_list)

setwd(wd)
org=data.frame( db=c('org.Hs.eg.db','org.Mm.eg.db'),sp=c('hsa','mmu'),stringsAsFactors = F)
row.names(org)= c('human','mouse')

orgDb=org[organism,'db']
sp=org[organism,'sp']

eg = bitr(row.names(x), fromType="SYMBOL", toType=c("ENTREZID","SYMBOL"),  OrgDb=orgDb)
foldchange = x[,1];  names(foldchange)=row.names(x)

enrichGO_plot=function(entrezID,ont,wd){

  ego= enrichGO(entrezID,  OrgDb=orgDb, ont = ont, pvalueCutoff = pvalueCutoff,
                  pAdjustMethod = "BH", qvalueCutoff = pvalueCutoff, minGSSize = 5,readable = T)
  setwd(wd)
  unlink(ont, recursive=TRUE)
  if(!dir.exists(ont)){ unlink(ont, recursive=TRUE);dir.create(ont)}
  setwd(ont)
  write.csv( as.data.frame(ego),'GOenrich.csv', row.names = F)
  pdf(file='barplot.pdf')
  myplot=barplot(ego,showCategory=10)
  print(myplot)
  dev.off()

  pdf(file='dotplot.pdf')
  myplot=dotplot(ego)
  print(myplot)
  dev.off()

  pdf(file='enrichMap.pdf')
  myplot=enrichMap(ego)
  print(myplot)
  dev.off()

  tryCatch({
    pdf(file='GOgraph.pdf')
    myplot=plotGOgraph(ego)
    print(myplot)
    dev.off()

  }, #cc??????
    error=function(e){cat(conditionMessage(e),"\n\n")},
    finally={print("error")})

  tryCatch({
    pdf(file="cnetplot.pdf")
    myplot=cnetplot(ego, categorySize="pvalue", foldChange=foldchange)
    print(myplot)
    dev.off()

  }, #cc??????
  error=function(e){cat(conditionMessage(e),"\n\n")},
  finally={print("error")} )

  setwd(wd)
}

turn2name = function(ids,eg_tbl){
  x=strsplit(ids,'/')[[1]]
  y=eg_tbl %>% filter(ENTREZID %in% x) %>% select(SYMBOL)
  yy=paste(y$SYMBOL,collapse = '/')
  return(yy)
}

wd_go=paste(wd,'GO',sep = '/')
if(!dir.exists(wd_go)){ dir.create(wd_go)}
enrichGO_plot(eg[,2],'MF',wd_go)
enrichGO_plot(eg[,2],'CC',wd_go)
enrichGO_plot(eg[,2],'BP',wd_go)

#kegg pathway
kk <- enrichKEGG(gene =  eg[,2], organism = sp,  pvalueCutoff = pvalueCutoff, 
                 pAdjustMethod = "BH", qvalueCutoff = qvalueCutoff , minGSSize = 3)
kks = as.data.frame(kk)
kks_copy =kks
setwd(wd)
if(!dir.exists('KEGG')){ dir.create('KEGG')}
setwd('KEGG')
library(dplyr)
eg_tbl= tbl_df(eg)
kks$geneID=sapply(kks$geneID, turn2name,eg)

write.csv(kks,'KEGG.csv',row.names = F)

detach("package:dplyr")
library("pathview")
b=50
if(b > length(kks$ID) ){
  b  = length(kks$ID)
}

for(id in kks$ID[1:b]){
  tryCatch({
    pathview(gene.data  = foldchange, pathway.id = id, gene.idtype = "SYMBOL", species = sp,
            limit      = list(gene=round( max(abs(foldchange)),2), cpd=1) ) }

      )
}
