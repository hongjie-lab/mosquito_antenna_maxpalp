suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)



#####################################

### 4. Neuron subcluster

normalizationMethod='LogNormalize'
# normalizationMethod='SCTransform' # not good for visualization
print(c('normalizationMethod: ', normalizationMethod))

receptorFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/3.checkReceptors')
neuronFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/4.neuronClusters')
dir.create(neuronFolderPath, showWarnings = FALSE)


#####################################

load( file = file.path(receptorFolderPath, 'srat.neuron.RData') )
picFolderPath=paste0(neuronFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)

print(neuronFolderPath)
print(picFolderPath)

head(srat.neuron@meta.data)


### Plot UMAP/tSNE

### (1)
sratDF <- NormalizeData(sratDF)
sratDF <- FindVariableFeatures(sratDF, selection.method = "vst", nfeatures = 2000)
sratDF <- ScaleData(sratDF, vars.to.regress = c('nCount_RNA'))
sratDF <- RunPCA(sratDF, npcs = 50, verbose = F)

sratDF <- FindNeighbors(sratDF, reduction = "pca", dims = 1:50)
sratDF <- FindClusters(sratDF, resolution = c(0.5, 0.8, 1, 2.5)[4])

### UMAP
sratDF <- RunUMAP(sratDF, reduction = "pca", dims = 1:50, metric = c("cosine", 'correlation', 'euclidean')[1]  )
sratDF <- RunTSNE(sratDF, reduction = "pca", dims = 1:50)

pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE', '.pdf')),
    width=10, height=5)
print( ( DimPlot(sratDF, reduction = "umap", label = TRUE) + theme(legend.position = "none")) + 
         ( DimPlot(sratDF, reduction = "tsne", label = TRUE) + theme(legend.position = "none")) )
dev.off()



srat.neuron <- NormalizeData(srat.neuron)
srat.neuron <- FindVariableFeatures(srat.neuron, selection.method = "vst", nfeatures = 2000)
srat.neuron <- ScaleData(srat.neuron, vars.to.regress = c('nCount_RNA'))

srat.neuron    <- RunPCA(srat.neuron, npcs = 50, verbose = F)

srat.neuron    <- RunUMAP(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- RunTSNE(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- FindNeighbors(srat.neuron, dims = 1:10, verbose = F)
srat.neuron    <- FindClusters(srat.neuron, verbose = T, resolution = c(0.5, 0.8, 2, 4)[2])


pdf(file = file.path(picFolderPath, 
                     paste0('UMAP_tSNE_', 'neuronOnly.pdf')),
    width=10, height=5)
print( ( DimPlot(srat.neuron, reduction = "umap", label = TRUE) + theme(legend.position = "none")) + 
        ( DimPlot(srat.neuron, reduction = "tsne", label = TRUE) + theme(legend.position = "none")) )
dev.off()

pdf(file = file.path(picFolderPath, 
                     paste0('nFeature_nCount_ptMt_ploidy.pdf')),
    width=12, height=3.5)
print(VlnPlot(sratDF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)) # ploidy difference in cluster #2?
dev.off()

print(VlnPlot(srat.neuron, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))


save.image(file = file.path(neuronFolderPath, 'neuralClusters.RData') )



