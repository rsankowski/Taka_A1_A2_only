library(tidyverse)
library(RaceID)
source('~/Documents/Single cell analysis/Advanced-plots/20190214-go_term-analysis.R')

date = Sys.Date()
load('data/sc.Robj')
source("/home/roman/Documents/Single cell analysis/Advanced-plots/20181025-sankowski-et-al-functions.R")

load("data/cluster-markers.RData") 
all.markers <- filter(all.markers, p_val_adj<0.05 )


enrich_up <- go_term_analysis_seurat(.df = all.markers)

dir.create('GO-terms')
dir.create('GO-terms/bp')
write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms.csv')

#mf terms
enrich_up_mf <- go_term_analysis_seurat(.df = all.markers,ontogeny = 'MF')

dir.create('GO-terms/mf')
write.csv(enrich_up_mf, 'GO-terms/mf/mf_GO_terms.csv')



#go term analysis of pseudo time
#mac
df <- read_csv('data/2019-05-29-nodes-stemid-vector-ctrl_mac.csv', col_names = c("X1", "Cluster", "GENEID"), skip = 1)

enrich_up <- go_term_analysis()

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-trajectory-mac.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',)

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-trajectory-mac.csv')

#GO terms across all nodes
df$Cluster <- rep(1, nrow(df))

rm(enrich_up)
#analysis
enrich_up <- go_term_analysis()

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-across-all-nodes-trajectory-mac.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',)

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-across-all-nodes-trajectory-mac.csv')


#modc
  df <- read_csv('data/2019-06-01-nodes-stemid-vector-ctrl_modc.csv', col_names = c("X1", "Cluster", "GENEID"), skip=1)
  
  enrich_up <- go_term_analysis()
  
  write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-trajectory-modc.csv')
  
  #mf
  enrich_up <- go_term_analysis(ontogeny = 'MF',)
  
  write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-trajectory-modc.csv')
  
  #GO terms across all nodes
  df$Cluster <- rep(1, nrow(df))
  
  #analysis
  enrich_up <- go_term_analysis()
  
  write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-across-all-nodes-trajectory-modc.csv')
  
  #mf
  enrich_up <- go_term_analysis(ontogeny = 'MF',)
  
  write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-across-all-nodes-trajectory-modc.csv')

  
  
  
                      #retain genes associated with certain go terms
                      #I downloaded te go term information from here: http://www.informatics.jax.org/vocab/gene_ontology/
                      #find genes in dataset
                      all_genes <- rownames(sc@ndata)[which(apply(as.matrix(sc@ndata) > 0,1,sum)>0)]
                      source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap.R')
                      
                      #before running loop please make sure to create the df object from plottin.R
                      for (i in list.files("data/go-terms-mirco")) {
                        genes <- read_tsv(paste0("data/go-terms-mirco/", i)) %>%
                          select(Symbol) %>%
                          unique %>%
                          unlist
                        
                        genes <- name2id(genes, rownames(sc@ndata))
                        present_genes <- genes[genes %in% all_genes]
                        
                        #plot_genes <-  as.character(unique(up_genes$GENEID))
                        svg(paste0('plots/tsne/', i, '.svg'), width = 8.57, height = 5.79)
                        pl <- plot_expmap(gene=present_genes, point_size = 2.2)
                        print(pl)
                        dev.off()
                        
                        svg(paste0('plots/tsne/', i, '-logsc.svg'), width = 8.57, height = 5.79)
                        pl <- plot_expmap(gene=present_genes, point_size = 2.2, logsc = T)
                        print(pl)
                        dev.off()
                        
                      }
                      