
library(Seurat)
library(ggplot2)
library(edgeR)
library(ggpubr)

metadata = read.csv("D:/bio/data/GSE228394_metadata",sep= '')
control <- Read10X("D:/bio/data/WT_control1/")
sus <- Read10X("D:/bio/data/WT_Sus/")
res <- Read10X("D:/bio/data/WT_Res/")
metadata$AllenV2[is.na(metadata$AllenV2)] <- 'AMY'
control <- CreateSeuratObject(counts = control)
sus <- CreateSeuratObject(counts = sus)
res <- CreateSeuratObject(counts = res)
data = AddMetaData(object = data, metadata = metadata)
data <- merge(control, y = c(sus, res), add.cell.ids = c("WT_control1", "WT_Sus", "WT_Res"))
data = JoinLayers(data)
seurat_cells <- colnames(data)
metadata_cells <- row.names(metadata)
cells_to_filter <- setdiff(seurat_cells, metadata_cells)
data <- subset(data, cells = seurat_cells[!(seurat_cells %in% cells_to_filter)])
data = AddMetaData(object = data, metadata = metadata)
DimPlot(data,group.by = 'BrainRegion')
# 归一化数据
data<- SCTransform(data,verbose=F)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data, features = VariableFeatures(object = data))
library(harmony)
data<-RunHarmony(data,"orig.ident")
data <- RunUMAP(data,dims = 1:40,reduction = 'pca') # 运行 UMAP

DimPlot(data,group.by = 'AllenV2')
FeaturePlot(data,features = 'Pkm',split.by = 'orig.ident', order = T)
VlnPlot(data,features = 'Pkm',group.by = "AllenV2",split.by = 'orig.ident',cols = c("red","orange","blue"))
# 设置默认图形宽度和高度（单位为英寸）
options(repr.plot.width = 6, repr.plot.height = 30)
DotPlot(data, features = 'Pkm',group.by = 'cluster1',split.by = 'sample',cols = c("red","orange","blue"))
RidgePlot(data, features = 'Pkm', group.by = 'sample')
con = subset(data, subset = sample=='CON')
count=GetAssayData(data, assay = 'RNA',layer = 'data')
count=log2(edgeR::cpm(count)+1)
# 使用 %in% 来筛选行
count = count[rownames(count) %in% c('Pkm',house.genes), ,drop = FALSE]
count = t(count)
count <- cbind(count,metadata)

'Pkm'%in% rownames(count)
'Actb' %in% rownames(count)

count <- as.data.frame(count)  # 将向量转换为数据框
count[] <- lapply(count, as.numeric) 
gene = "Pkm"
comparison = list(c("CON", "RES"), c("CON", "SUS"))
ggplot(data = count,aes(x=orig.ident,y=Pkm,fill=orig.ident))+ geom_boxplot() + theme(axis.text.x = element_text(size=0.1)) + 
  xlab(NULL) + ylab("Expression")+theme_bw()+facet_grid(.~AllenV2)+
  ggtitle('Pkm')+geom_jitter(size = 0.1,width = 0.1)+
  stat_compare_means(comparisons = comparison,aes(group=orig.ident),
                     method="wilcox.test",
                     label="p.signif",
                     label.y=by(count[,gene],count$AllenV2,max)+0.2*seq_along(comparison)-0.1,size=4.5)

house.genes <- c('Rrn18s','Actb','Gapdh','Pgk1','Ppia','Rpl13a','Rplp0','Arbp','B2m','Ywhaz','Sdha','Tfrc','Gusb','Hmbs','Hprt1','Tbp')
house.genes <- intersect(house.genes,rownames(data))
data <- PercentageFeatureSet(data,features = c('Pkm',house.genes),col.name = 'score')
table=table(metadata$AllenV2)
table = as.data.frame(table)
table$Freq = as.numeric(table$Freq)
ggplot(data = table, aes(x = Var1, y = Freq, fill = Freq)) +
  geom_bar(stat = "identity") +xlab(NULL) + ylab(NULL)+
  geom_text(aes(label = Freq), 
            position = position_stack(vjust = 0.5),  # 使文本居中
            size = 3, 
            color = "white")
