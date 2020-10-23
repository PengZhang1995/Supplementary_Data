#rm(list = ls())
#setwd("E:\\Desktop\\Supplementary_Data\\1.common DEGS\\GSE84005")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

library(GEOquery)
gset = getGEO('GSE84005',destdir = '.',getGPL = F,
              AnnotGPL = T)
gset = gset[[1]]         
expr = exprs(gset)        
pdata = pData(gset)       
gset@annotation           

probe = read.table(file = 'GPL5175-3188.txt',
                   sep = '\t',
                   quote = '',         
                   comment.char = '#',    
                   header = T,
                   fill = T,                     
                   stringsAsFactors = F) 
ids = probe[probe$gene_assignment != '',
            c(1,10)] 

library(dplyr)
colnames(ids)
expr = as.data.frame(expr)
#write.csv(expr,file = "expr.CSV", quote = F)
expr <- read.csv("expr.CSV")
exprSet = inner_join(ids,expr,by = 'ID')       
library(limma)
exprSet= avereps(exprSet[,-c(1,2)],             
                 ID = exprSet$gene_assignment)
exprSet = as.data.frame(exprSet)

pdf(file = 'Before_normalization.pdf')
par(pin=c(5.4, 3),family='serif',cex.axis=1.5)
p <- boxplot(exprSet,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F,ylim=c(1, 14))
title(main = list('GSE84005',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

library(limma)
normalized_expr = normalizeBetweenArrays(exprSet)
#normalized_expr=log2(normalized_expr)
pdf(file = 'Post_normalization.pdf')
par(pin=c(5.4, 3),family='serif',cex.axis=1.5)
p1 <- boxplot(normalized_expr,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F,ylim=c(1, 13))
title(main = list('GSE84005',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

group_list = pdata$source_name_ch1
#table(pdata$title) 
control = normalized_expr[,grep('^paired',group_list)]
tumor   = normalized_expr[,grep('^tumor',group_list)]
exprSet1 = cbind(tumor,control)
group_list = c(rep('tumor',ncol(tumor)),
               rep('normal',ncol(control)))
head(exprSet1)

data = exprSet1 

group_list = factor(group_list)
design <- model.matrix( ~0 + group_list)
colnames( design ) = levels(group_list)
rownames( design ) = colnames(data)
contrast.matrix <- makeContrasts( "tumor-normal", levels = design)

fit <- lmFit( data, design )
fit2 <- contrasts.fit( fit, contrast.matrix ) 
fit2 <- eBayes( fit2 )
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="alldiff.xls",sep="\t",quote=F)

logFoldChange=1
adjustP=0.05
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)

GSE84005_diff=diffSig
GSE84005_diff=GSE84005_diff[order(GSE84005_diff$logFC),]
GSE84005_diff=rbind(Gene=colnames(GSE84005_diff),GSE84005_diff)
write.table(GSE84005_diff,file="GSE84005_diff.txt",sep="\t",quote=F,col.names=F)

diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)

diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

hmExp=data[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)


xMax=max(-log10(allDiff$adj.P.Val))   
yMax=max(abs(allDiff$logFC))
library(ggplot2)
allDiff$change <- ifelse(allDiff$adj.P.Val < 0.05 & abs(allDiff$logFC) > 1,
                         ifelse(allDiff$logFC > 1,'UP','DOWN'),
                         'NOT')
table(allDiff$change)
pdf(file = 'volcano.pdf',width=6, height=4)
ggplot(data= allDiff, aes(x = -log10(adj.P.Val), y = logFC, color = change)) +
  geom_point(alpha=0.8, size = 1) +
  theme_bw() +
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"), #  标题居中
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+ # 网格线设置为空白
  geom_hline(yintercept= 0 ,linetype= 2 ) +
  scale_color_manual(name = "", 
                     values = c("red", "green", "black"),
                     limits = c("UP", "DOWN", "NOT")
  ) +
  xlim(0,xMax) + 
  ylim(-yMax,yMax) +
  labs(title = 'GSE84005', x = '-log10(adj.P.Val)', y = 'logFC',  colour="black")
dev.off()


top100exp = exprSet1[rownames(diffSig)[1:100],]
write.table(top100exp,file="top100exp.xls",sep="\t",quote=F)
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)

library(pheatmap)
pdf(file = 'heatmap.pdf',family='serif')
p=pheatmap(top100exp,
           main='GSE84005',#标题
           border_color =NA,#网格线
           cluster_rows = T,cluster_cols = T,#聚类
           #treeheight_col =0,#聚类树高度
           show_rownames = T,show_colnames = T,#行列名
           annotation_col = annotation_col,
           color = colorRampPalette(c("green", "black", "red"))(50),
           fontsize  = 5)
dev.off()
