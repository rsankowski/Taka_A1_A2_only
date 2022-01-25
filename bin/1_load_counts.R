# adjusted from Josip Herman
library(readr)
library(data.table)
library(tidyverse)
library(tools)
library(assertthat)

#create folders to structure project
dir.create("data")
dir.create("data/counts")
dir.create("GO-terms")
dir.create("GO-terms/bp")
dir.create("GO-terms/bp/tsne")
dir.create("GO-terms/mf")
dir.create("GO-terms/mf/tsne")
dir.create("plots")
dir.create("plots/heatmaps")
dir.create("plots/tsne")
dir.create("plots/others")

# Specify directories
file_paths <- list("data/counts/")
names(file_paths) <- file_paths

# Search for all files in the specified directories and extract files by a given extension
files_list            <- lapply(file_paths, list.files, recursive=T)
files_by_ext          <- lapply(files_list, function(x){x[endsWith(x, suffix=".coutt.csv")]} )

# Get complete paths to all files
all_file_paths        <- unlist(lapply(seq_along(files_by_ext), function(x) {  paste(names(files_by_ext[x]), files_by_ext[[x]], sep="") } ))
names(all_file_paths) <- lapply(strsplit(all_file_paths,split="/"), function(x) { sub(".coutt.csv","",x[length(x)]) } )

# Calculate md5sums to check for duplicated
md5sums      <- lapply(all_file_paths, function(x) {md5sum(x)} )

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(unname(unlist(lapply(strsplit(unlist(files_by_ext),split = "/"),tail,1))))) == 0)

####
#### LOADING
####
# Loading data using lapply
data_list   <- lapply(all_file_paths, function(x) {fread(x, header= T)} )

# Add dataset name prefix to all columns, Merge with remaining gene names
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))
  
}

# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- reduce(data_list, full_join, by = "GENEID")
data_list_cbind[is.na(data_list_cbind)]   <- 0

# Remove ERCC and mitochondrial genes
prdata <- data.frame(data_list_cbind)
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL

#save prdata with mt- cs <- colSums(prdata)
cs <- colSums(prdata)
prdata <- prdata[,names(cs)[cs > 499]]

#exclude not detected genes
prdata <- prdata[which(rowSums(prdata) >0),]

#save object
save(prdata, file =  "data/prdata_with_mt.RData")

prdata <- prdata[grep("^(mt-)",rownames(prdata),invert=TRUE),]

#remove Kcnq1ot1 and correlated genes
             cs <- cs[cs > 499]
          f <- t(prdata["Kcnq1ot1",])/cs < .02
          
          #define sigcor
          sigcor <- function(x,y,cthr=.4){
            if ( min(var(x),var(y)) == 0 ) return(NA)
            fit <- lm(x ~ y)
            pv <- as.data.frame(summary(fit)[4])[2,4]
            y <- as.data.frame(summary(fit)[4])[2,1]
            if ( is.na(pv) | is.na(y) ) return( NA )
            z <- sign(y)*sqrt(summary(fit)$r.square)
            if ( is.na(z) ) return(NA)
            if ( pv < .01 & abs(z) >= cthr ) return(z) else return(NA)
          }
          
          
          h <- rep(TRUE,nrow(prdata))
          for ( g in "Kcnq1ot1"){
            z <- apply(prdata,1,function(x,y) sigcor(x,y),y=t(prdata[g,]))
            h <- h & ( is.na(z) | z < .8 )
          }
          prdata <- prdata[h,]
          

save(prdata, file = "data/prdata.RData")
