library(STACAS)

# get test dataset
if (!requireNamespace("remotes")) install.packages("remotes")
if (!requireNamespace("SeuratData")) install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
options(timeout = max(300, getOption("timeout")))
InstallData("pbmcsca")
data("pbmcsca")

pbmcsca <- subset(pbmcsca, downsample=5000) # optionally take a sample for faster computation

# Integrate scRNA-seq datasets generated with different methods/technologies
pbmcsca.integrated <- NormalizeData(pbmcsca) |> SplitObject(split.by = "Method") |> Run.STACAS()
pbmcsca.integrated <- RunUMAP(pbmcsca.integrated, dims = 1:30) 

# Visualize
DimPlot(pbmcsca.integrated, group.by = c("Method","CellType")) 

# Semi-supervised integration taking into account cell labels/prior knowledge
pbmcsca.semisup <- NormalizeData(pbmcsca) |> SplitObject(split.by = "Method") |> Run.STACAS(cell.labels = "CellType")
pbmcsca.semisup <- RunUMAP(pbmcsca.semisup, dims = 1:30) 

# Visualize pre-integration
pbmcsca <- NormalizeData(pbmcsca) |> FindVariableFeatures() |> ScaleData() |> RunPCA() |> RunUMAP(dims = 1:30) 
DimPlot(pbmcsca, group.by = c("Method","CellType")) 

citation('pbmcsca.SeuratData')
