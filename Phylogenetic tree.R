###Phylogenetic tree


library(GUniFrac)
library(phyloseq) 
library(ggtreeExtra) 
library(ggstar) 
library(ggtree) 
library(treeio)
library(ggnewscale) 
library(tidytree)

tree<-read.newick("/Volumes/MAc KIOXIA Storage/gut/16s/tree/genus.phylogeny.tree",
                  node.label = "support") 

d<-read.csv("Genu-p1.csv",header=T) 

trs<-full_join(tree,d,by='label') 

tree@data$support<-ifelse(tree@data$support<50,NA,tree@data$support)


library(RColorBrewer)
library(scales)
getPalette = colorRampPalette(c(brewer.pal(8, "Accent"),brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(12, "Paired") ))(44)
g
show_col(getPalette,labels=T)


ggtree(trs, layout="circular", 
       #open.angle=10, #开口大小，circular时该条无效
       size=0.5, 
       branch.length="none" #末端对齐
)+
  #geom_tiplab(aes(color = Phylum),size=3, hjust= -0.06) + #显示物种信息，并设置颜色大小
  geom_tippoint(aes(color = Phylum), size=3) + #显示物种标识，并设置颜色大小
  #geom_text2(aes(subset=!isTip, label=node), hjust=-0.3, size=2, color="deepskyblue4")+ #显示节点支持率，并设置其位置、大小以及颜色
  geom_nodepoint(color="#FF6633", size=1.5)+ #显示节点标识及其颜色大小，alpha值为透明度
  #theme_tree2()  #显示坐标轴（绝对遗传距离）
  scale_color_manual(values = getPalette)