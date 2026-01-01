library(ggplot2)
library(ggpubr)
RL = read.table("D:/bio/data/GSE138614_WM_vs_RL.txt/GSE138614_WM_vs_RL.txt", header = T)
which(CA$Label=="SEMA5B")
RL$Label <- gsub("^[^;]*;", "", RL$Label)
rl = RL[6190,]
IL = read.table("C:/Users/LENOVO/Downloads/GSE138614_WM_vs_IL.txt/GSE138614_WM_vs_IL.txt", header = T)
IL$Label <- gsub("^[^;]*;", "", IL$Label)
il = IL[5526,]
NAWM = read.table("C:/Users/LENOVO/Downloads/GSE138614_WM_vs_NAWM.txt/GSE138614_WM_vs_NAWM.txt", header = T)
NAWM$Label = gsub("^[^;]*;", "", NAWM$Label)
nawm = NAWM[2561,]
AL = read.table("C:/Users/LENOVO/Downloads/GSE138614_WM_vs_AL.txt/GSE138614_WM_vs_AL.txt", header = T)
AL$Label = gsub("^[^;]*;", "", AL$Label)
al = AL[6340,]
CA = read.table("C:/Users/LENOVO/Downloads/GSE138614_WM_vs_CA.txt/GSE138614_WM_vs_CA.txt", header = T)
CA$Label = gsub("^[^;]*;", "", CA$Label)
ca = CA[3736,]
result = rbind(rl,il,nawm,al,ca)
rownames(result) = c("RL","IL","NAWM","AL","CA")
result$status = rownames(result)
result = result[c(5,3,4,1,2),]
ggbarplot(result, x = "status", y = "logFC",
              fill = "logFC",           # change fill color by mpg_level
              color = "white",            # Set bar border colors to white # jco journal color palett. see ?ggpar
              sort.val = "none",          # Sort the value in descending order
              sort.by.groups = FALSE,     # Don't sort inside each group
              x.text.angle = 90,       # Rotate vertically x axis texts
              ylab = "logFC",
              rotate = T,
              ggtheme = theme_minimal(),
              title = "SEMA5B"
)+theme(axis.text.y = element_text(size = 10))
