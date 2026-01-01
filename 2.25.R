library(Seurat)
wt = Read10X("D:/bio/data/GSE242666_RAW/WT/")
tm60 = Read10X("D:/bio/data/GSE242666_RAW/60TM/")
wt = CreateSeuratObject(wt)
tm60 = CreateSeuratObject(tm60)                        

wt@meta.data$group = 'WT'
tm60@meta.data$group = "60TM"
data = merge(wt,tm60, add.cell.ids = c('WT','60TM'))
data = JoinLayers(data)

data = NormalizeData(data)
data = FindVariableFeatures(data)
data = ScaleData(data)
data <- RunPCA(data, verbose = FALSE)
data <- FindNeighbors(data, dims = 1:15)
data <- FindClusters(data, resolution = 0.5, verbose = FALSE)
data <- RunUMAP(data, dims = 1:15)
DimPlot(data, label = TRUE,reduction = "umap")
DimPlot(data,group.by = "group")
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
library(dplyr)
data.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
library(ggplot2)
DoHeatmap(data, features = top10$gene, slot = "count", disp.min = -1) + NoLegend()+theme(axis.text.y = element_text(size = 5.5))
cluster.id = c('0'='Microglia',
              '1' ='Astro',
              '2' ='Microglia',
              '3' ='Oligo',
              '4' ='Astro',
              '5' ='Astro',
              '6' ='Neuron',
              '7' ='Astro',
              '8' ='Oligo',
              '9' ='VECs',
              '10' ='NA',
              '11' ='Neuron',
              '12' ='Oligo',
              '13' ='Oligo',
             '14' ='BAMs',
             '15' = 'NA',
             '16' = 'NA',
             '17' = 'NA',
             '18' = 'NA',
             '19'='Pericytes',
             '20'='VECs',
             '21'='Oligo',
             '22'= 'Erythrocyte',
             '23' = 'NA')
cluster.id <- c('Microglia', 'Astro', 'Microglia', 'Oligo', 'Astro', 'Astro', 
                       'Neuron', 'Astro', 'Oligo', 'VECs', 'NA', 'Neuron', 'Oligo', 
                       'Oligo', 'BAMs', 'NA', 'NA', 'NA', 'NA', 'Pericytes', 'VECs', 
                       'Oligo', 'Erythrocyte', 'NA')
DimPlot(data, label = T)
names(cluster.id) <- levels(data)
data <- RenameIdents(data, cluster.id)
FeaturePlot(data, features = 'Sema5b',slot = 'data',split.by = 'group',order = T)
VlnPlot(data,features = 'Sema5b',group.by = 'celltype',split.by = 'group')
as = subset(data, subset = cell_type=='Astro')
FeaturePlot(as,features = 'Sema5b',order = T,split.by = 'group')
library(ggpubr)
VlnPlot(as,features = 'Sema5b', group.by = 'group')
data$cell_type <- Idents(data)
library(clustree)
as <- FindClusters(as, resolution = 0.5, verbose = FALSE)
DimPlot(as, group.by = "RNA_snn_res.0.5",split.by = "group",label = T)
DotPlot(as, features = "Sema5b")+facet_grid(cols = vars(RNA_snn_res.0.5))
AS6 = subset(as, subset = seurat_clusters=="6")
p5<-DotPlot(AS5,features="Sema5b",assay="RNA", group.by = "group",
            scale = T,cols=c("#f0dd93","#bd4335"))+
  labs(x="AC5")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position="none",
        panel.background=element_blank(),
        panel.spacing=unit(1,"line"),
        panel.border=element_rect(color="black",fill=NA,size=1),
        axis.title.y=element_blank())
library(patchwork)
p0+p1+p2+p3+p4+p5+plot_layout(ncol = 6)
VlnPlot(as, features = "Sema5b",slot = "data",group.by = "RNA_snn_res.0.5",split.by = "group")





library(plot1cell)
complex_dotplot_multiple(seu_obj = as, features = "Sema5b", group = "group")
