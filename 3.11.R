library(biomaRt)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(DESeq2)
count = read.table("D:/bio/data/GSE138614_countMatrix.txt/GSE138614_countMatrix.txt")

ensembl = useEnsembl(biomart = 'genes')
ensembl <- useMart("ensembl", host = "http://useast.ensembl.org") 
datasets <- listDatasets(ensembl)
dim(datasets)
head(datasets)
ensembl = useDataset(dataset ='hsapiens_gene_ensembl', mart = ensembl)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = 'hsapiens_gene_ensembl', 
                      mirror = "www")
count$gene_ids<-rownames(count)
count$gene_ids<- gsub("\\..*", "", count$gene_ids)
count=aggregate(.~gene_ids,mean,data=count)

gene_info<-getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                 filters="ensembl_gene_id",
                 values=count$gene_ids,
                 mart=ensembl)

count = count[count$gene_ids %in% gene_info$ensembl_gene_id,]
count = cbind(gene_info[,2],count)
colnames(count)[1] = "gene_name"

sample_info = read.delim("D:/bio/data/MS_sample_info.txt")
count = as.data.frame(t(count))
count_used = subset(matrix,gene_name%in%c("PLXNA2","PLXNA1","PLXNA4","PLXNA3","SEMA5B","SEMA5A"))
rownames(count_used)=count_used$gene_name
count_used = count_used[,-1]
count_used = as.data.frame(t(count_used))
count_used = cbind(count_used,sample_info)

sample_info$SampleType = gsub("MS.*","MS",sample_info$SampleType)
sample_info$SampleType = gsub("C.*","Control",sample_info$SampleType)
count_used <- as.data.frame(lapply(count_used, as.numeric))
dds <- DESeqDataSetFromMatrix(countData = count_used, colData = sample_info, design = ~Status )
dds = DESeq(dds)
res <- results(dds)
head(res)
vst<-assay(vst(dds,blind=FALSE))
vst = cbind(count$gene_name, vst)
vst = as.data.frame(vst)
colnames(vst)[1] <- "gene_name"
matrix <- aggregate(.~gene_name,data = vst,max)
rownames(matrix)=matrix$gene_name
count_used$Status = gsub("White.*","WM", count_used$Status)
count_used$Status = gsub("C.*","CA", count_used$Status)
count_used$Status = gsub("Normal.*","NAWM", count_used$Status)
count_used$Status = gsub("Active.*","AL", count_used$Status)
count_used$Status = gsub("R.*","RL", count_used$Status)
count_used$Status = gsub("Inactive.*","IL", count_used$Status)

comparisons = list(c("WM", "NAWM"),c("WM", "AL"),c("WM", "CA"),c("WM", "IL"),c("WM", "RL"))
ggplot(data = count_used,aes(x=Status,y=SEMA5B,fill=SampleType))+ geom_boxplot() + theme(axis.text.x = element_text(size=0.1)) + 
  xlab(NULL) + ylab("Expression")+theme_bw()+
  ggtitle('SEMA5B')+geom_jitter(size = 0.1,width = 0.1)+
  stat_compare_means(comparisons = comparisons,method="t.test",label="p.signif")
