library(ActivePathways)
library(qs)
library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
data = qread("D:/bio/data/NMO.integrated.48492cells.integrated.harmony.qs")
DimPlot(data,group.by = 'sctype_classification')
astro = subset(data,subset = sctype_classification=="Astrocytes")
DimPlot(astro, group.by = 'treatment')
astro = JoinLayers(astro)
diff <- FindMarkers(astro, group.by = "treatment",
                               ident.1 = "nmo",
                               ident.2 = "con",
                               only.pos = FALSE,
                               test.use = "wilcox",
                               assay = "RNA")

score = as.data.frame(diff[,1])
rownames(score)=rownames(diff)
names(score) = "p.val"
gmt = read.GMT("D:/bio/data/gene_symbol_to_go.gmt")
enriched_pathways <- ActivePathways(score, "D:/bio/data/gene_symbol_to_go.gmt",
                                    correction_method = "fdr",
                                    geneset_filter = c(10,750),
                                    cytoscape_file_tag = "enrichmentMap__")
score = as.matrix(score)
# 获取 Gene Symbol 到 GO term 的映射
gene_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")  # 获取所有 Gene Symbol

# 查询每个 Gene Symbol 对应的 GO term
go_annotations <- AnnotationDbi::select(org.Mm.eg.db,
                                        keys = gene_symbols,
                                        columns = c("GO","ONTOLOGY"),  # 返回 GO term 和 GO 类别（比如 Biological Process, Molecular Function）
                                        keytype = "SYMBOL")
go_bp <- go_annotations[go_annotations$ONTOLOGY == "BP", ]
go_annotations$TERM = Term(go_annotations$GO)
go_list <- split(go_annotations$SYMBOL, go_annotations$GO)

# 将 GO term 列表转换为 GMT 格式
go_ids = names(go_list)
go_term_names <- Term(go_ids)
gmt_list <- lapply(go_ids, function(go_id) {
  genes <- unique(go_list[[go_id]])
  term_name <- if (go_id %in% names(go_term_names)) go_term_names[[go_id]] else "NA"
  c(go_id, term_name, genes)
})


gmt_file <- "D:/bio/data/gene_symbol_to_go.gmt"
con <- file(gmt_file, "w")

for (line in gmt_list) {
  writeLines(paste(line, collapse = "\t"), con)
}

close(con)


setwd("D:/bio/data/")




library(dplyr)
# 按 GO term 分组，整理成 GMT 格式
gmt_data <- go_annotations %>%
  group_by(GO, TERM, ONTOLOGY) %>%  # 按 GO ID 和 Term 分组
  summarise(genes = paste(SYMBOL, collapse = "\t")) %>%  # 合并基因
  mutate(description = "")  # 添加空白描述列（可选）

# 查看整理后的数据
head(gmt_data)
gmt_lines <- paste(
  gmt_data$GO,              # GO ID
  gmt_data$TERM,      # 描述（可换成 gmt_data$TERM）
  gmt_data$genes,            # 基因列表（用 tab 分隔）
  sep = "\t"
)

# 写入文件
writeLines(gmt_lines, "gene_symbol_to_go.gmt")


sema_pathway = enriched_pathways[c(129,641,1546),]
sema_pathway = as.data.frame(sema_pathway)
write.csv(sema_pathway,file = 'sema_pathways.csv')
str(sema_pathway)
library(ggplot2)
ggplot(sema_pathway, aes(x = term_size,y= term_name))+
  geom_point(aes(size = term_size,color=adjusted_p_val))+
  scale_color_gradient(low='red',high = "grey")
  
FeaturePlot(object = astro,features = "Sema6a",split.by = "treatment")




Idents(astro) = astro$treatment
up_genes <- subset(diff, subset=diff$avg_log2FC > 0)
down_genes <- subset(diff, subset=diff$avg_log2FC < 0)
up_score = as.data.frame(up_genes[,1])
down_score = as.data.frame(down_genes[,1])
rownames(down_score)=rownames(down_genes)
names(down_score) = "p.val"
down_score = as.matrix(down_score)
enriched_pathways_up <- ActivePathways(up_score, "D:/bio/data/gene_symbol_to_go.gmt",

                                       cytoscape_file_tag = "up_pathways")
enriched_pathways_down <- ActivePathways(down_score, "D:/bio/data/gene_symbol_to_go.gmt",
                                       cytoscape_file_tag = "down_pathways")




enriched_pathways_up$group <- "up_regulated"
enriched_pathways_down$group <- "down_regulated"
enriched_pathways_all <- rbind(enriched_pathways_up,enriched_pathways_down)

library(dplyr)
pathways_filtered <- enriched_pathways_all %>% filter(term_name %in% enriched_pathways$term_name)

write.table(enriched_pathways_all, file = "output.txt", sep = "\t", quote = FALSE, row.names = FALSE)



enriched_pathways_all$term_name[duplicated(enriched_pathways_all$term_name)]

enriched_pathways_all$up_regulated <- ifelse(enriched_pathways_all$term_name %in% enriched_pathways_up$term_name, 1, 0)
enriched_pathways_all$down_regulated <- ifelse(enriched_pathways_all$term_name %in% enriched_pathways_down$term_name, 1, 0)

enriched_pathways_all <- enriched_pathways_all %>%
  distinct(term_id, .keep_all = TRUE)
subgroups = data.frame(term_id = enriched_pathways_all$term_id,
                       term_name = enriched_pathways_all$term_name,
                       up_regulated = enriched_pathways_all$up_regulated,
                       down_regulated = enriched_pathways_all$down_regulated)
utils::write.table(subgroups, 
                   file=paste0("All__", "subgroups.txt"), 
                   row.names=FALSE, 
                   sep="\t", 
                   quote=FALSE)


 all(enriched_pathways_up$term_id %in% subgroups$term_id)

enriched_pathways_all = enriched_pathways_all[,-"overlap"]
gmt_main <- gmt[enriched_pathways_all$term_id]
write.GMT(gmt_main,paste0("ALL__","pathways.gmt"))
