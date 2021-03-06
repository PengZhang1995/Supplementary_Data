#rm(list = ls())
#setwd("E:\\Desktop\\Supplementary_Data\\1.common DEGS\\GSE84402")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

library(GEOquery)
gset = getGEO('GSE84402',destdir = '.',getGPL = F,
              AnnotGPL = T)
gset = gset[[1]]          
expr = exprs(gset)        
pdata = pData(gset)      
gset@annotation           

probe = read.table(file = 'GPL570-55999.txt',
                   sep = '\t',
                   quote = '',         
                   comment.char = '#',   
                   header = T,
                   fill = T,               
                   stringsAsFactors = F)  
ids = probe[probe$Gene.Symbol != '',
            c(1,11)] 

library(dplyr)
colnames(ids)
expr = as.data.frame(expr)
#write.csv(expr,file = "expr.CSV", quote = F)
expr <- read.csv("expr.CSV")
exprSet = inner_join(ids,expr,by = 'ID')       
library(limma)
exprSet= avereps(exprSet[,-c(1,2)],             
                 ID = exprSet$Gene.Symbol)
exprSet = as.data.frame(exprSet)
#exprSet=log2(exprSet+1)

pdf(file = 'Before_normalization.pdf')
par(pin=c(5, 3),family='serif',cex.axis=1.5)
p <- boxplot(exprSet,outline=FALSE,las=2,col = 'blue',xaxt = 'n',ann = F,ylim=c(0, 300))
title(main = list('GSE84402',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3.5,font = 2,cex = 1.5)
dev.off()

library(limma)
normalized_expr = normalizeBetweenArrays(exprSet)
normalized_expr=log2(normalized_expr)
pdf(file = 'Post_normalization.pdf')
par(pin=c(5, 3),family='serif',cex.axis=1.5)
p1 <- boxplot(normalized_expr,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F,ylim=c(-2.5, 11))
title(main = list('GSE84402',cex = 2 ,font = 2),
      xlab = list('Sample list',cex = 1.5,font = 2),
      ylab = '',line = 0.7)
mtext('Expression value',side = 2,padj = -3,font = 2,cex = 1.5)
dev.off()

group_list = pdata$source_name_ch1
control = normalized_expr[,grep('non-cancerous',group_list)]
tumor = normalized_expr[,grep('hepatocellular',group_list)]
exprSet1 = cbind(tumor,control)
group_list = c(rep('tumor',ncol(tumor)),
               rep('normal',ncol(control)))
head(exprSet1)
# exprSet1 = log2(exprSet1)  

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

GSE84402_diff=allDiff
GSE84402_diff=GSE84402_diff[order(GSE84402_diff$logFC),]
GSE84402_diff=rbind(Gene=colnames(GSE84402_diff),GSE84402_diff)
write.table(GSE84402_diff,file="GSE84402_diff.txt",sep="\t",quote=F,col.names=F)

logFoldChange=1
adjustP=0.05
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
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
  theme_bw(base_size = 18) +
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+ 
  geom_hline(yintercept= 0 ,linetype= 2 ) +
  scale_color_manual(name = "", 
                     values = c("red", "green", "black"),
                     limits = c("UP", "DOWN", "NOT")
  ) +
  xlim(0,xMax) + 
  ylim(-yMax,yMax) +
  labs(title = 'GSE84402', x = '-log10(adj.P.Val)', y = 'logFC',  colour="black")
dev.off()

top100exp = exprSet1[rownames(diffSig)[1:100],]
#write.table(top100exp,file="top100exp.xls",sep="\t",quote=F)
annotation_col = data.frame(group = group_list)
rownames(annotation_col) = colnames(exprSet1)
library(pheatmap)
pdf(file = 'heatmap.pdf',family='serif')
pheatmap(top100exp,
         main='GSE84402',
         border_color =NA,
         cluster_rows = T,cluster_cols = T,
         #treeheight_col =0,
         show_rownames = T,show_colnames = T,
         annotation_col = annotation_col,
         color = colorRampPalette(c("green", "black", "red"))(50),
         fontsize  = 5)
dev.off()
