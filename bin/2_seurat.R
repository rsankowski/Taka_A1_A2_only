#RaceID4 
library(tidyverse)
library(viridis)
library(RaceID)
library(Seurat)
library(Matrix)
library(clustree)

date = Sys.Date()

load("data/prdata_with_mt.RData")
prdata <- prdata[, !grepl("EMP", colnames(prdata))]

prdata <- prdata[,!colnames(prdata) %in% read_csv("data/damaged-cells.csv")[[1]]] %>%
  CreateSeuratObject(project="a1_a2_emp", min.cells = 5, min.features = 500)

prdata$percent.mt <- PercentageFeatureSet(prdata,pattern="^mt-")
prdata$Celltype <- factor(ifelse(grepl("EMP", colnames(prdata)), "EMP", "A1_A2"), levels = c("EMP", "A1_A2"))

#normalize dataset 
all <- prdata %>%
  SCTransform(vars.to.regress = c("percent.mt"),
                          variable.features.n = 10000)
all <-all %>%
  RunPCA(features=grep("^(Kcnq1ot1|Gm|Rpl|Rps|Fos|Jun|mt-|Dusp|Zfp36|Hsp|Hbb|Hba)|(Rik)$", VariableFeatures(.), value=T, invert = T)) 

#
ElbowPlot(all)
all<- all %>% 
  RunUMAP(dims=1:15) %>%
  FindNeighbors(dims=1:15) %>%
  FindClusters(resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1, 1.5,2))

#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

pdf("plots/others/overview-cluster-resolutions.pdf")
clustree(all)
dev.off()

all<-FindClusters(all,resolution=2) #for overclustering resolution 3

FeaturePlot(all, "Mrc1")
FeaturePlot(all, "S100a9")
DimPlot(all, label = TRUE) + NoLegend()
DimPlot(all, label = TRUE, group.by = "Celltype")

save(all, file = "data/all.RData")

#find cluster markers
all.markers<-FindAllMarkers(all,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

save(all.markers, file = "data/diffgenes-seurat.RData")
write_csv(all.markers, "data/diffgenes-seurat.csv")

#write_csv(data.frame(ID=names(all@active.ident)[all@active.ident==2]), "data/damaged-cells.csv")
