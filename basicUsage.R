library(STACAS)

# get test dataset
if (!requireNamespace("remotes")) install.packages("remotes")
if (!requireNamespace("SeuratData")) install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)
options(timeout = max(300, getOption("timeout")))
InstallData("pbmcsca")
data("pbmcsca")

# Integrate scRNA-seq data generated with different methods/technologies
#pbmcsca <- subset(pbmcsca, downsample=5000) # optionally take a sample for faster computation
pbmcsca <- pbmcsca |> NormalizeData() |> SplitObject(split.by = "Method") |> Run.STACAS() |> RunUMAP(dims = 1:30) 

# Visualize
DimPlot(pbmcsca, group.by = c("Method","CellType")) 

# Visualize pre-integration
#pbmcsca.pre <- pbmcsca |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA() |> RunUMAP(dims = 1:30) 
#DimPlot(pbmcsca.pre, group.by = c("Method","CellType")) 

citation('pbmcsca.SeuratData')
