##Download datasets from GEO and ArrayExpress for STACAS demo

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE))
   BiocManager::install("GEOquery")
if (!requireNamespace("DropletUtils", quietly = TRUE))
   BiocManager::install("DropletUtils")

if (!requireNamespace("Seurat", quietly = TRUE)) {
  BiocManager::install("multtest")
  install.packages("Seurat")
}

library(GEOquery)
library(DropletUtils)

########################################################################################################################
#Make sure we can write to this directory
datadir <- "data.STACAS.demo"
dir.create(datadir)

#########################################################################################################################
#1. Carmona et al. (GSE116390)
geo_acc <- "GSE116390"

destfile = paste0(datadir,"/GSE116390_aggregated_matrix.tar.gz")
download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE116390&format=file&file=GSE116390_aggregated_matrix.tar.gz",
              destfile = destfile)

untar(destfile, exdir=paste0(datadir,"/SJC"))

data <- Read10X(paste0(datadir,"/SJC"))

c <- gsub(pattern="(\\S+)", replacement="CA-\\1", colnames(data), perl=TRUE)
colnames(data) <- c
carmona.seurat <- CreateSeuratObject(counts=data, project = "Carmona_CD8", min.cells = 0, min.features = 50)

#########################################################################################################################
#2. Xiong et al. (E-MTAB-7919)
arrayexpress_acc <- "E-MTAB-7919"

destfile = paste0(datadir,"/E-MTAB-7919.1.zip")
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7919/E-MTAB-7919.processed.1.zip", destfile=destfile)
unzip(destfile, exdir=paste0(datadir,"/Xiong"))

destfile = paste0(datadir,"/E-MTAB-7919.2.zip")
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7919/E-MTAB-7919.processed.2.zip", destfile=destfile)
unzip(destfile, exdir=paste0(datadir,"/Xiong"))

destfile = paste0(datadir,"/E-MTAB-7919.3.zip")
download.file("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7919/E-MTAB-7919.processed.3.zip", destfile=destfile)
unzip(destfile, exdir=paste0(datadir,"/Xiong"))

data <- Read10X(paste0(datadir,"/Xiong"))

c <- gsub(pattern="(\\S+)", replacement="XI-\\1", colnames(data), perl=TRUE)
colnames(data) <- c
xiong.seurat <- CreateSeuratObject(counts=data, project = "Xiong_TILs", min.cells = 0, min.features = 50)

xiong.seurat$Sample <- factor(substring(colnames(xiong.seurat),21))
sampleSize <- table(xiong.seurat$Sample)

#######################################################################################################################
#3. Magen et al. (GSE124691)
geo_acc <- "GSE124691"

gse <- getGEO(geo_acc)

getGEOSuppFiles(geo_acc,baseDir=datadir)
system(paste0("gunzip ",datadir,"/",geo_acc,"/*Genes.tsv.gz")) 

accs <- c("GSM3543446","GSM3543449","GSM3543443","GSM3543447")
for (acc in accs) {
  system(paste0("mkdir ",datadir,"/",acc))
  getGEOSuppFiles(acc,baseDir=datadir) 
  
  fns <- list.files(paste0(datadir,"/",acc,"/"))
  for (fn in fns){
    system(paste0("mv ",paste0(datadir,"/",acc,"/",fn)," ",paste0(datadir,"/",acc,"/",sub(".*_","",fn,perl = T))))
  }
  system(paste0("gunzip ",datadir,"/",acc,"/*gz"))
  
  #add gene list
  system(paste0("cp ",datadir,"/",geo_acc,"/*Genes.tsv ", datadir,"/",acc,"/genes.tsv"))
}

geoAccToName <- gse[[1]]$title
names(geoAccToName) <-  gse[[1]]$geo_accession
geoAccToName

###Magen TILs
acc1 <- accs[1]
acc2 <- accs[2]

fname <- paste0(datadir, "/",acc1)
datasetID <- "Magen_CD4_TIL"
data <- Read10X(fname)
c <- gsub(pattern="(\\S+)", replacement=paste0("MA-\\1-",1), colnames(data), perl=TRUE)
colnames(data) <- c

data.seurat.1 <- CreateSeuratObject(counts = data, project = datasetID, min.cells = 0, min.features = 50)
data.seurat.1@meta.data$Sample <- geoAccToName[acc1]

fname <- paste0(datadir, "/",acc2)
datasetID <- "Magen_CD4_TIL"
data <- Read10X(fname)
c <- gsub(pattern="(\\S+)", replacement=paste0("MA-\\1-",2), colnames(data), perl=TRUE)
colnames(data) <- c

data.seurat.2 <- CreateSeuratObject(counts = data, project = datasetID, min.cells = 0, min.features = 50)
data.seurat.2@meta.data$Sample <- geoAccToName[acc2]

magen.til.seurat <- merge(data.seurat.1,data.seurat.2)

rm(data.seurat.1)
rm(data.seurat.2)

###Magen lymphnode
acc1 <- accs[3]
acc2 <- accs[4]

fname <- paste0(datadir, "/",acc1)
datasetID <- "Magen_CD4_dLN"
data <- Read10X(fname)
c <- gsub(pattern="(\\S+)", replacement=paste0("MA-\\1-",3), colnames(data), perl=TRUE)
colnames(data) <- c

data.seurat.1 <- CreateSeuratObject(counts = data, project = datasetID, min.cells = 0, min.features = 50)
data.seurat.1@meta.data$Sample <- geoAccToName[acc1]

fname <- paste0(datadir, "/",acc2)
datasetID <- "Magen_CD4_dLN"
data <- Read10X(fname)
c <- gsub(pattern="(\\S+)", replacement=paste0("MA-\\1-",4), colnames(data), perl=TRUE)
colnames(data) <- c

data.seurat.2 <- CreateSeuratObject(counts = data, project = datasetID, min.cells = 0, min.features = 50)
data.seurat.2@meta.data$Sample <- geoAccToName[acc2]

magen.dln.seurat <- merge(data.seurat.1,data.seurat.2)

rm(data.seurat.1)
rm(data.seurat.2)

raw.list=list("CD8_TILs"=carmona.seurat, "CD8CD4_TILs"=xiong.seurat, "CD4_TILs" = magen.til.seurat,  "CD4_dLN"=magen.dln.seurat)
saveRDS(raw.list, file="raw.list.demo.rds")


