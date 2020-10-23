#rm(list = ls())
#setwd("E:\\Desktop\\Supplementary_Data\\2.GO_KEGG")

all = read.table(file = 'DEGs_GO.TXT',sep = '\t',header = T,quote = '')
all_rt = all[all$PValue < 0.05,]
library(tidyr)
all_rt = separate(all_rt, Term, sep = "~",
                 into = c("ID", "Term"))

bp_df = all_rt[all_rt$Category == 'GOTERM_BP_DIRECT',]
bp_df = bp_df[order(bp_df$Count,decreasing = T),]
bp = bp_df[1:5,]

cc_df = all_rt[all_rt$Category == 'GOTERM_CC_DIRECT',]
cc_df = cc_df[order(cc_df$Count,decreasing = T),]
cc = cc_df[1:5,]

mf_df = all_rt[all_rt$Category == 'GOTERM_MF_DIRECT',]
mf_df = mf_df[order(mf_df$Count,decreasing = T),]
mf = mf_df[1:5,]

allGo = rbind(bp,cc,mf)
library(stringr)
table(allGo$Category)
allGo$Category = substr(allGo$Category,8,9)

library(ggpubr)
colnames(allGo)
p = ggbarplot(data = allGo,x = "Term",y = 'Count',
              
              x.text.angle = 45,
              fill = "Category",
              palette = c("cadetblue3","mediumslateblue","mediumorchid3"),
              sort.by.groups = T,
              xlab = '',
              ylab = "Target genes") +
  theme( 
    text=element_text(family="serif"),
    axis.text.x=element_text(color = "black",size = rel(1)),
    axis.text.y=element_text(color = "black",size = rel(1)),
    axis.title.x=element_text(color = "black",size = rel(1)),
    axis.title.y=element_text(color = "black",size = rel(1)),
    legend.text = element_text(color = "black",size = rel(1)),
    legend.title  = element_text(color = "black",size = rel(1)),
    plot.margin = unit(c(0.02,1.7,0.02,1.7),"cm")
  )
#ggpar(p,x.text.angle = 90)
ggsave(plot = p,'DEGs_GO.tiff',width = 8.8,height =6.4)
ggsave(plot = p,'DEGs_GO.pdf',width =8.8,height =6.4)

KEGG= read.table(file = 'DEGs_KEGG.txt',sep = '\t',header = T,quote = '')
keggSig = KEGG[KEGG$PValue < 0.01,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))

library(ggplot2)
p0 <- ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    y="Pathway name",
    title="Pathway enrichment")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(size=rel(1.0)),
        axis.title.y = element_blank()
  )
ggsave('DEGs_KEGG.tiff',p0,width = 8,height = 4)

KEGG= read.table(file = 'subnetwork1_KEGG.txt',sep = '\t',header = T,quote = '')
keggSig = KEGG[KEGG$PValue < 0.05,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))

library(ggplot2)
p1 <- ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    y="Pathway name",
    title="Pathway enrichment")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(size=rel(1.0)),
        axis.title.y = element_blank()
  )
ggsave('subnetwork1_KEGG.tiff',p1,width = 7,height = 4)


KEGG= read.table(file = 'subnetwork2_KEGG.txt',sep = '\t',header = T,quote = '')
keggSig = KEGG[KEGG$PValue < 0.05,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))

library(ggplot2)
p2 <- ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    y="Pathway name",
    title="Pathway enrichment")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(size=rel(1.0)),
        axis.title.y = element_blank()
  )
ggsave('subnetwork2_KEGG.tiff',p2,width = 5.5,height = 4)


KEGG= read.table(file = 'subnetwork3_KEGG.txt',sep = '\t',header = T,quote = '')
keggSig = KEGG[KEGG$PValue < 0.05,]
library(tidyr)
keggSig = separate(keggSig, Term, sep = ":",
                   into = c("ID", "Term"))

library(ggplot2)
p3 <- ggplot(keggSig,aes(x=Fold.Enrichment,y=Term)) + 
  geom_point(aes(size=Count,color=-1*log10(PValue)))+
  scale_colour_gradient(low="green",high="red")+
  labs(
    color=expression(-log[10](P.value)),
    size="Gene number",
    x="Fold enrichment",
    y="Pathway name",
    title="Pathway enrichment")+
  theme_bw()+
  theme(plot.title=element_text(hjust=0.5),
        axis.text = element_text(color = "black"),
        text=element_text(family="serif"),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.x = element_text(size=rel(1.0)),
        axis.title.y = element_blank()
  )
ggsave('subnetwork3_KEGG.tiff',p3,width = 7,height = 4)
