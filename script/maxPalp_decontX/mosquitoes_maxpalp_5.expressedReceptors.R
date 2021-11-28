suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
  library(ComplexHeatmap)
  library(dichromat)
  library(RColorBrewer)
  library(chorddiag)
})

theme_set(theme_cowplot())
projectFolder='/Users/tzuchiao/ProjectAnalysis/Leslie_Vosshall'
setwd(projectFolder)


#####################################

### 5. expressed receptors

normalizationMethod='LogNormalize'
print(c('normalizationMethod: ', normalizationMethod))

neuronFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/4.neuronClusters')

################################################################################################################################
load(file = file.path(neuronFolderPath, 'neuralClusters.RData'))

receptorFolderPath=paste0('Analysis/10282021_MaxPalp/', normalizationMethod, '/5.expressedReceptors')
dir.create(receptorFolderPath, showWarnings = FALSE)
picFolderPath=paste0(receptorFolderPath, '/pic')
dir.create(picFolderPath, showWarnings = FALSE)





### Check nFeature_RNA mtTx in the cluster level
p2 = VlnPlot(sratDF, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = 'seurat_clusters', pt.size = 0)

pdf(file = file.path(picFolderPath, 
                     paste0('nFeature_nCount_ptMt_', normalizationMethod, '_afterFiler.pdf')),
    width=10, height=5)
print(p2)
dev.off()



###
### coreceptor: not necessary for using neuron-specific cluster
pdf(file = file.path(picFolderPath, 
                     paste0('coreceptor_FeaturePlot_neuronOnly_batchCorrect.pdf')),
    width=15, height=10)
print(FeaturePlot(srat.neuron, reduction = 'tsne', features = c('Orco', 'Ir25a', 'Ir8a', 'Ir76b', 'Gr3'), ncol = 3, order = TRUE, pt.size = 1))
dev.off()


### Receptor list: for plotting all receptors
allGeneID.df<- as.character(read.table('SeuratFile/features.tsv', header = FALSE)$V1)
# head(allGeneID.df)

receptorL=c()
lapply(allGeneID.df, function(x) {
  # print(x[[1]])
  if ( isFALSE(startsWith(x, 'LOC'))  ) {
    if ( isFALSE(startsWith(x, 'MT'))  ) {
      receptorL <<- append(receptorL, x)
    }
  }
})


### Plot receptor using all MaxPalp cells
receptor.dot <- DotPlot(sratDF, features = receptorL) + scale_y_discrete(limits = rev) +
  xlab('') + ylab('') + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size = 9)) + scale_y_discrete(limits = rev)
# receptor.dot
receptor.dot.value <- receptor.dot$data %>% as.data.frame() %>%
  filter(id %in% neuronClusterIDs)

coReceptorL <- c('Orco', 'Ir25a', 'Ir76b', 'Ir8a', 'Gr3')

receptor.dot.value.noCoR <- receptor.dot.value %>% filter(! features.plot %in% coReceptorL)
receptor.dot.value.noCoR <- receptor.dot.value.noCoR %>% arrange(desc(avg.exp))

receptorExpSortL <- unique(receptor.dot.value.noCoR$features.plot)
receptorExpSortL <- c(as.character(receptorExpSortL), coReceptorL)
head(receptorExpSortL)
tail(receptorExpSortL)

clusterExpSortL <- unique(receptor.dot.value.noCoR$id) %>% as.character()
# clusterExpSortL


hex_codes1 <- hue_pal()(20)
scatterColors <- c('#66c2a5', '#fc8d62', '#8da0cb')

# relaxed criteria
pmain <- receptor.dot.value.noCoR %>%
  ggplot(aes(avg.exp, pct.exp, color=id)) +
  geom_point(size=0.5)+scale_x_log10()+
  geom_hline(yintercept=c(25,50)[1], color=scatterColors[3]) +
  xlab('Average expression') + ylab('Percentage of expressing cell')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = receptor.dot.value.noCoR, aes(x = avg.exp), fill='grey') + scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = receptor.dot.value.noCoR, aes(x = pct.exp), fill='grey') +
  geom_vline(xintercept = c(25,50)[1], color=scatterColors[3]) +
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath,
                     paste0('scatterPlot_avg.exp_pct.exp_allchemoreceptors_25Percent.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()


chemo25PerL <- receptor.dot.value.noCoR %>% filter(pct.exp >= 25) %>% 
  pull(features.plot) %>% as.character() %>%
  unique() %>% sort()
chemo25PerL


chemoreceptorInData <- 
  rownames(sratDF[["RNA"]]@data)[rownames(sratDF[["RNA"]]@data) %in% receptorL & (! rownames(sratDF[["RNA"]]@data) %in% coReceptorL)] %>% sort()
chemoreceptorCoRInData <- 
  rownames(sratDF[["RNA"]]@data)[rownames(sratDF[["RNA"]]@data) %in% receptorL ] %>% sort()

neuronMarkerL <- c('LOC5564848', 'LOC5570381', 'LOC5565901', 'LOC5570204')
# Epithelia, Muscle, Glia markers
OtherMarkerL <- c('LOC5564305', 'LOC5570152', 'LOC5565327',
                  'LOC5566990', 'LOC5580173', 'LOC5580231',
                  'LOC110678282', 'LOC5577632', 'LOC5566721')



### ordered heatmap

# neuronClusterIDs

sratDF@meta.data$annotation <- apply(sratDF@meta.data, 1, function(input){
  clusterID=input['seurat_clusters']
  
  if (clusterID %in% neuronClusterIDs) {
    return(paste('Neuron', clusterID, sep = '_'))
  } 
  else if (clusterID %in% c('1', '8')) {
    return('Epithelia')
  }
  else if (clusterID %in% c('14')) {
    return('Muscle')
  }
  else if (clusterID %in% c('15')) {
    return('Glia')
  }
  else {return('Others')}
})

sratDF@meta.data$annotation <- factor(sratDF@meta.data$annotation,
                                      levels = c('Others', 'Epithelia', 'Muscle', 'Glia', 'Neuron_2', 'Neuron_4', 'Neuron_5', 'Neuron_7'))
OtherMarkerL <- c('LOC5564305', 'LOC5570152', 'LOC5565327', 
                  'LOC5566990', 'LOC5580173', 'LOC5580231',
                  'LOC110678282', 'LOC5577632', 'LOC5566721')
coReceptorL <- c('Orco', 'Ir25a', 'Ir76b', 'Ir8a', 'Gr3')
chemo25PerL <- c('Gr1', 'Gr2', 'Or8', 'Or49')
# LOC5575210: nompC ortholog
# 
OtherNeuronMarkerL <- c('LOC5575210', 'LOC23687649')
heatmapGeneL <- c(OtherMarkerL, neuronMarkerL, coReceptorL, chemo25PerL, OtherNeuronMarkerL)


### Defined clusters

allReceptorMat <- sratDF[["RNA"]]@data[ heatmapGeneL, ] %>% as.matrix()
# allReceptorMat<- t(scale(t(allReceptorMat)))

quantileL <- quantile(allReceptorMat, c(0.1, 0.95))
quantileL <- append(quantileL, mean(quantileL), after = 1)
# col_fun = circlize::colorRamp2(quantileL, colorRampPalette(brewer.pal(9, "YlGnBu"))(3) ) 
# col_fun = circlize::colorRamp2(quantileL, c('blue', 'yellow') ) 
col_fun = circlize::colorRamp2(quantileL, c('navy', 'black', 'yellow') ) 

# col_fun = circlize::colorRamp2(quantileL, c("black", "yellow")) 
# colours=colorRampPalette( c('black', 'yellow'))(24)

cluster_anno<- sratDF@meta.data$annotation
left_anno <- c(replicate(3, 'Epithelia Markers'),
               replicate(3, 'Muscle Markers'),
               replicate(3, 'Glia Markers'),
               replicate(4, 'Neuron Markers'), 
               replicate(5, 'Coreceptors'),
               replicate(4, 'Dominant Chemoreceptors'),
               replicate(2, 'Mechanosensory Neuron Markers')
               )
left_anno <- factor(left_anno,
                    levels = c('Epithelia Markers', 'Muscle Markers', 'Glia Markers',
                               'Neuron Markers', 'Coreceptors', 'Dominant Chemoreceptors', 'Mechanosensory Neuron Markers'))



pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_allCells_mainMarkers_order1_navy-black-yellow.pdf')),
    width=20, height=6)

print(Heatmap(allReceptorMat, name = "Expression",
              column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
              column_title_rot = 0,
              height = unit(0.75, "npc"),
              
              col = col_fun, 
              show_column_names = FALSE, use_raster = FALSE,
              
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 10),
              
              row_names_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              cluster_row_slices = TRUE,
              
              
              row_split = left_anno,
              row_title_gp = gpar(fontsize = 10),
              row_title_rot = 0,
              
              
              
              top_annotation = 
                HeatmapAnnotation(
                  foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(sratDF@meta.data$annotation)) ) ))
                ),
              # left_annotation =
              #   rowAnnotation(
              #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )), 
              #     ),
              
))
dev.off()

### No others
### Defined clusters
heatmapGeneL <- c(OtherMarkerL, neuronMarkerL, coReceptorL, chemo25PerL, OtherNeuronMarkerL)
sratDF_annotated <- subset(sratDF, subset = annotation != 'Others')
sratDF_annotated <- subset(sratDF, subset = seurat_clusters %in% neuronClusterIDs)

allReceptorMat <- sratDF_annotated[["RNA"]]@data[ heatmapGeneL, ] %>% as.matrix()
# allReceptorMat<- t(scale(t(allReceptorMat)))

quantileL <- quantile(allReceptorMat, c(0.1, 0.95))
quantileL <- append(quantileL, mean(quantileL), after = 1)
# col_fun = circlize::colorRamp2(quantileL, colorRampPalette(brewer.pal(9, "YlGnBu"))(3) ) 
# col_fun = circlize::colorRamp2(quantileL, c('blue', 'yellow') ) 
col_fun = circlize::colorRamp2(quantileL, c('navy', 'black', 'yellow') ) 

cluster_anno<- sratDF_annotated@meta.data$annotation
left_anno <- c(replicate(3, 'Epithelia Markers'),
               replicate(3, 'Muscle Markers'),
               replicate(3, 'Glia Markers'),
               replicate(4, 'Neuron Markers'), 
               replicate(5, 'Coreceptors'),
               replicate(4, 'Dominant Chemoreceptors'),
               replicate(2, 'Mechanosensory Neuron Markers')
)
left_anno <- factor(left_anno,
                    levels = c('Epithelia Markers', 'Muscle Markers', 'Glia Markers',
                               'Neuron Markers', 'Coreceptors', 'Dominant Chemoreceptors', 'Mechanosensory Neuron Markers'))

pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_neuronCells_mainMarkers_order2_navy-black-yellow.pdf')),
    width=20, height=6)
# pdf(file = file.path(picFolderPath, 
    #                  paste0('heatmap_allCells_mainMarkers_order2_YlGnBu_test.pdf')),
    # width=20, height=6)
print(Heatmap(allReceptorMat, name = "Expression",
              column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
              column_title_rot = 0,
              height = unit(0.75, "npc"),
              
              col = col_fun, 
              show_column_names = FALSE,
              
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 10),
              
              row_names_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              cluster_row_slices = TRUE,
              
              
              row_split = left_anno,
              row_title_gp = gpar(fontsize = 10),
              row_title_rot = 0,              
              
              top_annotation = 
                HeatmapAnnotation(
                  foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(sratDF_annotated@meta.data$annotation)) ) ))
                ),
              # left_annotation =
              #   rowAnnotation(
              #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )), 
              #     ),
))
dev.off()


###
chemoGeneL <- c('Ir93a','Ir100a','Ir161','Ir41a')
length(chemoGeneL)

heatmapGeneL <- c(OtherMarkerL, neuronMarkerL, coReceptorL, chemo25PerL, OtherNeuronMarkerL, chemoGeneL)
sratDF_annotated <- subset(sratDF, subset = annotation != 'Others')
sratDF_annotated <- subset(sratDF, subset = seurat_clusters %in% neuronClusterIDs)

allReceptorMat <- sratDF_annotated[["RNA"]]@data[ heatmapGeneL, ] %>% as.matrix()
# allReceptorMat<- t(scale(t(allReceptorMat)))

quantileL <- quantile(allReceptorMat, c(0.1, 0.95))
# quantileL <- append(quantileL, mean(quantileL), after = 1)
col_fun = circlize::colorRamp2(quantileL, colorRampPalette(brewer.pal(9, "YlGnBu"))(3) ) 
col_fun = circlize::colorRamp2(quantileL, c('blue', 'yellow') ) 

cluster_anno<- sratDF_annotated@meta.data$annotation
left_anno <- c(replicate(3, 'Epithelia Markers'),
               replicate(3, 'Muscle Markers'),
               replicate(3, 'Glia Markers'),
               replicate(4, 'Neuron Markers'), 
               replicate(5, 'Coreceptors'),
               replicate(4, 'Dominant Chemoreceptors'),
               replicate(2, 'Mechanosensory Neuron Markers'),
               replicate(4, 'Other Chemoreceptors')
)
left_anno <- factor(left_anno,
                    levels = c('Epithelia Markers', 'Muscle Markers', 'Glia Markers',
                               'Neuron Markers', 'Coreceptors', 'Dominant Chemoreceptors', 
                               'Mechanosensory Neuron Markers', 'Other Chemoreceptors'))
# length(left_anno)

pdf(file = file.path(picFolderPath, 
                     paste0('heatmap_neuronCells_mainMarkers_order2_blue-yellow_moreChemo.pdf')),
    width=20, height=6)
# pdf(file = file.path(picFolderPath, 
#                  paste0('heatmap_allCells_mainMarkers_order2_YlGnBu_test.pdf')),
# width=20, height=6)
print(Heatmap(allReceptorMat, name = "Expression",
              column_split = cluster_anno, # column_split = factor(cluster_anno, levels(seurat.combined@meta.data$annotation2)),
              column_title_rot = 0,
              height = unit(0.75, "npc"),
              
              col = col_fun, 
              show_column_names = FALSE,
              
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              cluster_column_slices = TRUE,
              column_title_gp = gpar(fontsize = 10),
              
              row_names_gp = gpar(fontsize = 8),
              column_gap = unit(0.5, "mm"),
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              cluster_row_slices = TRUE,
              
              
              row_split = left_anno,
              row_title_gp = gpar(fontsize = 10),
              row_title_rot = 0,
              
              
              
              top_annotation = 
                HeatmapAnnotation(
                  foo = anno_block(gp = gpar(fill = scales::hue_pal()( length( levels(sratDF_annotated@meta.data$annotation)) ) ))
                ),
              # left_annotation =
              #   rowAnnotation(
              #     foo = anno_block(gp = gpar(fill = scales::hue_pal()( length(levels(left_anno)) ) )), 
              #     ),
))
dev.off()






##############################
##############################
### scatter plot & DotPlot


### matrix values
mat <- receptor.dot.value %>%
  select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = pct.exp) %>% as.data.frame()
row.names(mat) <- mat$id
row.names(mat)
mat <- mat[,-1] #drop gene column as now in rows

mat_avg.exp <- receptor.dot.value %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame()
row.names(mat_avg.exp) <- mat_avg.exp$id
row.names(mat_avg.exp)
mat_avg.exp <- mat_avg.exp[,-1] #drop gene column as now in rows


mat.T <- transpose(mat)
# t(mat)
colnames(mat.T) <- rownames(mat)
rownames(mat.T) <- colnames(mat)

mat.T$max.pct <- apply(mat.T, 1, function(input){
  # print(max(as.numeric(input)))
  return(max(as.numeric(input)))
})
mat.T$geneName <- row.names(mat.T)
head(mat.T)



## 25% of the cell
receptroMax25L <- mat.T %>% filter(max.pct >= 25) %>% pull(geneName)
receptorExpSort_max25L <- receptorExpSortL[receptorExpSortL %in% receptroMax25L]
receptor.dot.value.max25 <- receptor.dot.value %>% filter(features.plot %in% receptorExpSort_max25L)

dotplot <- receptor.dot.value.max25 %>%
  mutate(`% Expressing` = pct.exp,
         features.plot = factor(features.plot, levels =  receptorExpSortL),
         id = factor(id, levels = clusterExpSortL)
  ) %>% 
  ggplot(aes(x=features.plot, y = id, color = avg.exp, size = `% Expressing`)) + 
  geom_point() + 
  theme_half_open() +
  background_grid() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab('') +ylab('') +
  theme(axis.ticks = element_blank()) +
  scale_color_gradientn(colours = viridis::viridis(20),  name = 'Average\nexpression', oob = scales::squish, limits = c(1,20))+
  scale_y_discrete(position = "right", limits = rev) +
  scale_size(limits = c(25, 100)) 

pdf(file = file.path(picFolderPath, 
                     paste0('neuronOnly_DotPlot_receptor_25Percent.pdf')),
    width=5, height=4)
print(dotplot)
dev.off()



### Check expression of chemoreceptors
pdf(file = file.path(picFolderPath, 
                     paste0('FeaturePlot_Or8_Or49_Ir100a_Ir93a.pdf')),
    width=10, height=10)
print(FeaturePlot(sratDF, reduction = 'tsne', features = c('Or8', 'Or49', 'Ir100a', 'Ir93a'), ncol = 2, order = TRUE, pt.size = 1))
dev.off()







##########################################################
### Plotting scatter plot

### Orco vs. Ir25a
# head(srat.neuron@assays$RNA@counts[c('Orco','Ir25a'),])

srat.neuron$Orco_UMIs <- srat.neuron@assays$RNA@counts['Orco',]
srat.neuron$Ir25a_UMIs <- srat.neuron@assays$RNA@counts['Ir25a',]
srat.neuron$Gr3_UMIs <- srat.neuron@assays$RNA@counts['Gr3',]

srat.neuron$Orco_Exp <- srat.neuron@assays$RNA@data['Orco',]
srat.neuron$Ir25a_Exp <- srat.neuron@assays$RNA@data['Ir25a',]
srat.neuron$Gr3_Exp <- srat.neuron@assays$RNA@data['Gr3',]

### Exp: Orco vs. Ir25a
pmain <- srat.neuron@meta.data %>%
  ggplot( aes(Orco_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.5)
  # coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = srat.neuron@meta.data, aes(x = Orco_Exp), fill=scatterColors[1])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = srat.neuron@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Orco_Ir25a_normalizedExp.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()

### Exp: Gr3 vs. Ir25a
pmain <- srat.neuron@meta.data %>%
  ggplot( aes(Gr3_Exp, Ir25a_Exp) ) + 
  geom_point(size=0.5)
# coord_cartesian(xlim = c(0, 100))
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = srat.neuron@meta.data, aes(x = Gr3_Exp), fill=scatterColors[1])
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = srat.neuron@meta.data, aes(x = Ir25a_Exp), fill=scatterColors[2]) + 
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")

pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Gr3_Ir25a_normalizedExp.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()




### Gr3 vs. Ir25a in cluster level ==> not necessary
gr3Ir25a_avgExp.clusters.df <- 
  receptor.dot.value %>%
  filter(features.plot %in% c('Gr3', 'Ir25a')) %>%
  select(-pct.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = features.plot, values_from = avg.exp) %>% as.data.frame() 

pmain=  gr3Ir25a_avgExp.clusters.df %>%
  ggplot(aes(Gr3, Ir25a, color=id)) + 
  geom_point() +
  geom_text_repel(aes(label=id), hjust=0) + 
  theme(legend.position = "none") +
  scale_x_log10() + scale_y_log10()+
  xlab('Average expression of Gr3') + ylab('Average expression of Ir25a')
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = gr3Ir25a_avgExp.clusters.df, aes(x = Gr3), fill=scatterColors[1])+
  scale_x_log10()
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = gr3Ir25a_avgExp.clusters.df, aes(x = Ir25a), fill=scatterColors[2]) + 
  scale_x_log10()+
  coord_flip()
p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")


pdf(file = file.path(picFolderPath, 
                     paste0('scatterPlot_Gr3_Ir25a_normalizedExp_clusterLevel.pdf')),
    width=5, height=5)
print(ggdraw(p2))
dev.off()



##########################################################
### chord plot


# highReceptorRNA = srat.neuron@assays$RNA@counts[receptorExpSortL[1:20],]
# highReceptorRNA.data = srat.neuron@assays$RNA@data[receptorExpSortL[1:20],]

highReceptorRNA = as.matrix(srat.neuron@assays$RNA@counts[receptorExpSortL[1:20],])
highReceptorRNA.data = as.matrix(srat.neuron@assays$RNA@data[receptorExpSortL[1:20],])

highReceptorRNA.data = as.matrix(srat.neuron@assays$RNA@data[c(receptorExpSortL[1:20], 'Gr3'),])
highReceptorRNA.data = highReceptorRNA.data[-c(1, 4),] # Remove Gr1/2


highReceptorRNA.expCol=highReceptorRNA[,colSums(highReceptorRNA) != 0]
highReceptorRNA.data.expCol=highReceptorRNA.data[,colSums(highReceptorRNA.data) != 0]


# coexp.df <- data.frame(row.names = rownames(highReceptorRNA.data.expCol), )
coexp.df <- data.frame(matrix(ncol = 19, nrow = 19))
colnames(coexp.df) <- rownames(highReceptorRNA.data.expCol)
rownames(coexp.df) <- rownames(highReceptorRNA.data.expCol)
# coexp.df


### temp check srat.neuron
srat.neuron@meta.data$orig.ident
nrow(srat.neuron@meta.data)
levels(sratDF@meta.data$annotation)
nrow(sratDF@meta.data[sratDF@meta.data$annotation %in% c("Neuron_2", "Neuron_4", "Neuron_5", "Neuron_7"),])


# normalized expression
for (i in 1:19) {
  iExpressedCells.mx = highReceptorRNA.data.expCol[, highReceptorRNA.data.expCol[i,] >= 1 ]
  
  for (j in 1:19) {
    print(c('i = ', i))
    print(c('j = ', j))
    if (i==j) { coexp.df[i,j] <- 0 }
    else {
      lapply(1:length(ncol(highReceptorRNA.data.expCol)), function(x){
        print(c('i = ', i))
        print(c('j = ', j))
        iExpressedCells.mx = as.matrix(highReceptorRNA.data.expCol[, highReceptorRNA.data.expCol[i,] >= 1 ]) 
        jExpressedN=NCOL(iExpressedCells.mx[, iExpressedCells.mx[j,] >= 1]) 

        # print(jExpressedN)
        coexp.df[i,j]<<-jExpressedN

      })
    }
  }
}



# No cell# filter
chorddiag(as.matrix(coexp.df), groupColors = hex_codes1, 
          groupnamePadding = c( 15 ), 
          showTicks = TRUE, groupPadding=1, margin=120, )


coexp.df.filter10 <- coexp.df
coexp.df.filter10[coexp.df <= 10] <- 0
# check genes with 0 coexpression
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter10)) {
  if (sum(coexp.df.filter10[,i]) == 0) {
    # print(colnames(coexp.df.filter10)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter10 <- coexp.df.filter10[-zeroGenes_V, -zeroGenes_V]
chorddiag(as.matrix(coexp.df.filter10), groupColors = hue_pal()(ncol(coexp.df.filter10)), 
          groupnamePadding = 30, 
          showTicks = TRUE, groupPadding=1, margin=120)


coexp.df.filter20 <- coexp.df
coexp.df.filter20[coexp.df <= 20] <- 0
zeroGenes_V <- c()
for (i  in 1:ncol(coexp.df.filter20)) {
  if (sum(coexp.df.filter20[,i]) == 0) {
    print(colnames(coexp.df.filter20)[i])
    zeroGenes_V <- append(zeroGenes_V, i)
  }
}
coexp.df.filter20 <- coexp.df.filter20[-zeroGenes_V, -zeroGenes_V]
chorddiag(as.matrix(coexp.df.filter20), groupColors = hue_pal()(ncol(coexp.df.filter20)), 
          groupnamePadding = 30, 
          showTicks = TRUE, groupPadding=1, margin=120)





