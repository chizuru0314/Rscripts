library(survival)
library(survminer)
fpkm = read.csv("D:/bio/data/TCGA-BLCA.htseq_fpkm.tsv.gz.csv")
survival = read.table("D:/bio/data/TCGA-BLCA.survival.tsv",header = T)
#处理数据
fpkm <- as.data.frame(fpkm)
fpkm$X ->rownames(fpkm)
fpkm <- fpkm[, -1]  # 如果 x 是第一列
fpkm<- t(fpkm)
rownames(fpkm) <- gsub("\\.", "-", rownames(fpkm))
fpkm = cbind(sample = rownames(fpkm), fpkm)
a = setdiff(survival$sample,rownames(fpkm))
survival <- survival[!survival$sample %in% a, ]
b = setdiff(rownames(fpkm),survival$sample)
fpkm <- fpkm[!rownames(fpkm) %in% b, ]
#合并两个数据
data <- merge(survival, fpkm, by = "sample")

#转换为数字型变量
gene = tail(colnames(data),295)
for (i in gene){
  data[[i]] <- as.numeric(as.character(data[[i]]))
}


gene_used = "A2M"#选择基因

##1.中值点
cutpoint = median(data[[gene_used]])
data$group <- ifelse(data[[gene_used]] >= cutpoint, "High", "Low")#分组
fit<-survfit(Surv(OS.time,OS)~group,data=data)
diff = survdiff(Surv(OS.time,OS)~group,data=data)#计算P值
pval_logrank = diff[['pvalue']]
#cox回归分析
data['group'] <- lapply(data['group'],factor)
cox_model<-coxph(formula = Surv(OS.time,OS)~group,data=data)
sum=summary(cox_model)
sum$coefficients
sum$conf.int
HR=sum$coefficients[,2]
pval_cox=sum$coefficients[,5]
ci.high = sum$conf.int[,4]
ci.low = sum$conf.int[,3]
plot = ggsurvplot(fit,
           data=data,
           title = 'Median',
           pval=F,
           conf.int=F,
           fun="pct",
           xlab="Time",
           palette=c("red","black"),
           legend.title=ggplot2::element_blank(),
           legend.labs=c("High","Low"),
           surv.median.line="hv",#添加中位生存时间
           tables.height=0.2,
           ggtheme=theme_bw())

plot$plot+annotate("text", x=4000, y = 85,
              label = paste0("cutpoint=",round(cutpoint,3),
                              '\n','pval(Logrank)=',round(pval_logrank,3),
                             '\n','pval(Cox)=',round(pval_cox,3),
                             '\n','HR=',round(HR,3),
                             '\n','95%CI=',round(ci.low,3),"~",round(ci.high,3)),
             size = 5, color = "black")
          


##2.最佳截断点
#获取最佳截断点
res.cut <- surv_cutpoint(
  data,                   # 数据集
  time = "OS.time",        # 生存时间的列名（字符串形式）
  event = "OS",            # 生存状态的列名（字符串形式）
  variables = gene,
  minprop = 0
)
summary = summary(res.cut)
summary[] <- (lapply(summary, function(x) round(x, 3)))
cutpoint = summary[gene_used,'cutpoint']
#分组
res.cat<-surv_categorize(res.cut)
res.cat[[gene_used]] <- as.factor(res.cat[[gene_used]])
fit<-survfit(Surv(OS.time,OS)~res.cat[[gene_used]],data=res.cat)#拟合生存分析
diff = survdiff(Surv(OS.time,OS)~res.cat[[gene_used]],data=res.cat)#计算P值
pval_logrank = diff[['pvalue']]
#cox回归分析
cox_model<-coxph(formula = Surv(OS.time,OS)~res.cat[[gene_used]],data=res.cat)
sum=summary(cox_model)
sum$coefficients
sum$conf.int
HR=sum$coefficients[,2]
pval_cox=sum$coefficients[,5]
ci.high = sum$conf.int[,4]
ci.low = sum$conf.int[,3]
plot = ggsurvplot(fit,
                  data=res.cat,
                  title = 'Optimal',
                  pval=F,
                  conf.int=F,
                  fun="pct",
                  xlab="Time",
                  palette=c("red","black"),
                  legend.title=ggplot2::element_blank(),
                  legend.labs=c("High","Low"),
                  surv.median.line="hv",#添加中位生存时间
                  tables.height=0.2,
                  ggtheme=theme_bw())

plot$plot+annotate("text", x=4000, y = 85,
                   label = paste0("cutpoint=",round(cutpoint,3),
                                  '\n','pval(Logrank)=',round(pval_logrank,3),
                                  '\n','pval(Cox)=',round(pval_cox,3),
                                  '\n','HR=',round(HR,3),
                                  '\n','95%CI=',round(ci.low,3),"~",round(ci.high,3)),
                   size = 5, color = "black")




##3.25%截断点
cutpoint = quantile(data[[gene_used]], 0.25, na.rm = TRUE)
data$group <- ifelse(data[[gene_used]]>= cutpoint, "High", "Low")
fit<-survfit(Surv(OS.time,OS)~group,data=data)
diff = survdiff(Surv(OS.time,OS)~group,data=data)#计算P值
pval_logrank = diff[['pvalue']]
#cox回归分析
data['group'] <- lapply(data['group'],factor)
cox_model<-coxph(formula = Surv(OS.time,OS)~group,data=data)
sum=summary(cox_model)
sum$coefficients
sum$conf.int
HR=sum$coefficients[,2]
pval_cox=sum$coefficients[,5]
ci.high = sum$conf.int[,4]
ci.low = sum$conf.int[,3]
plot = ggsurvplot(fit,
                  data=data,
                  title = 'Quartile(25%)',
                  pval=F,
                  conf.int=F,
                  fun="pct",
                  xlab="Time",
                  palette=c("red","black"),
                  legend.title=ggplot2::element_blank(),
                  legend.labs=c("High","Low"),
                  surv.median.line="hv",#添加中位生存时间
                  tables.height=0.2,
                  ggtheme=theme_bw())

plot$plot+annotate("text", x=4000, y = 85,
                   label = paste0("cutpoint=",round(cutpoint,3),
                                  '\n','pval(Logrank)=',round(pval_logrank,3),
                                  '\n','pval(Cox)=',round(pval_cox,3),
                                  '\n','HR=',round(HR,3),
                                  '\n','95%CI=',round(ci.low,3),"~",round(ci.high,3)),
                   size = 5, color = "black")



##4.75%截断点
cutpoint = quantile(data[[gene_used]], 0.75, na.rm = TRUE)
data$group <- ifelse(data[[gene_used]] >= cutpoint, "High", "Low")
fit<-survfit(Surv(OS.time,OS)~group,data=data)
diff = survdiff(Surv(OS.time,OS)~group,data=data)#计算P值
pval_logrank = diff[['pvalue']]
#cox回归分析
data['group'] <- lapply(data['group'],factor)
cox_model<-coxph(formula = Surv(OS.time,OS)~group,data=data)
sum=summary(cox_model)
sum$coefficients
sum$conf.int
HR=sum$coefficients[,2]
pval_cox=sum$coefficients[,5]
ci.high = sum$conf.int[,4]
ci.low = sum$conf.int[,3]
plot = ggsurvplot(fit,
                  data=data,
                  title = 'Quartile(75%)',
                  pval=F,
                  conf.int=F,
                  fun="pct",
                  xlab="Time",
                  palette=c("red","black"),
                  legend.title=ggplot2::element_blank(),
                  legend.labs=c("High","Low"),
                  surv.median.line="hv",#添加中位生存时间
                  tables.height=0.2,
                  ggtheme=theme_bw())

plot$plot+annotate("text", x=4000, y = 85,
                   label = paste0("cutpoint=",round(cutpoint,3),
                                  '\n','pval(Logrank)=',round(pval_logrank,3),
                                  '\n','pval(Cox)=',round(pval_cox,3),
                                  '\n','HR=',round(HR,3),
                                  '\n','95%CI=',round(ci.low,3),"~",round(ci.high,3)),
                   size = 5, color = "black")



