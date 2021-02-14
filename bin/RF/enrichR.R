
bar_gokegg.5 =function(dat, plotname){
library(dplyr)
library(pheatmap)
library(ggplot2)
#library(ggthemes)
library(readr)

    dat %>% mutate(plog = - log10( P.value)) %>% filter( P.value<0.05 )  ->dat_st
    len = 10
    if(dim(dat_st)[1] <len ){
      len = dim(dat_st)[1]
    }
    dat_st = dat_st[1:len,]
    
    dat_st = arrange(dat_st, Combined.Score )
    dat_st$Term = factor(dat_st$Term, levels = dat_st$Term)
    p <- ggplot(data = dat_st, aes(x=Term, y=Combined.Score ) ) +  labs( y = "Enrichr Combined Score") +
      geom_bar(stat="identity", fill= '#FF6666' ) + 
      theme(axis.text.y = element_blank() ,axis.ticks.y  = element_blank()) +
      coord_flip() +
      geom_text(aes(x=Term,y=0.05, label=Term), color="black", size=8, hjust=0)   +   theme(plot.subtitle = element_text(vjust = 1), 
                                                                                            plot.caption = element_text(vjust = 1), 
                                                                                            axis.title = element_text(size = 20), 
                                                                                            axis.text = element_text(size = 18), 
                                                                                            
                                                                                            panel.background = element_rect(fill = NA)) +
      theme(plot.margin=unit(c(0.5,1.2,0.5,0),"cm") ,axis.title.y = element_blank() )
    #p <- ggplotly(p)
    #p <- ggplotly(p)
    p
    ggsave(filename = plotname,device = 'pdf', width=12, height=8)
}


my_enrichr = function(deGeneName, outdir){

library(dplyr)
library(pheatmap)
library(ggplot2)
#library(ggthemes)
library(readr)

library(enrichR)
 wd = getwd()
  cmd_mdkir = paste('mkdir -p ', outdir)
  system(cmd_mdkir)
  print(outdir)
  
  setwd(outdir)
  
  dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018','GO_Molecular_Function_2018', 'KEGG_2016', 'KEGG_2019_Human',
           'ChEA_2016','Transcription_Factor_PPIs', 'PPI_Hub_Proteins', 'TargetScan_microRNA_2017','miRTarBase_2017' )
  
  enriched <- enrichr(deGeneName, dbs)
  Sys.sleep(1)
  
  
  
  ## ----message=FALSE, warning=FALSE----------------------------------------
  
  
  for (db in names( enriched) ){
    print(db)
    
    out_result = paste(db, '.csv', sep='')
    write_csv(enriched[[db]], out_result )
    outfile = paste(db, '.pdf', sep='')
    bar_gokegg.5(enriched[[db]], outfile)
    
  }
  save(enriched,db, file="enriched.Rdata")
  setwd(wd)
}

argv = commandArgs(trailingOnly=TRUE)

print(argv)
library(tidyverse)
df = read_csv(argv[1])
print(df$gene_name)

my_enrichr(df$gene_name, "enrichR_out")