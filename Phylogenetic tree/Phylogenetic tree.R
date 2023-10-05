###Phylogenetic tree


library(GUniFrac)
library(phyloseq) 
library(ggtreeExtra) 
library(ggstar) 
library(ggtree) 
library(treeio)
library(ggnewscale) 
library(tidytree)



tree<-read.newick("Phylogeny.genus.tree",
                  node.label = "support") 
d<-read.csv("Phylogeny.lable.csv",header=T) 

trs<-full_join(tree,d,by='label') 

tree@data$support<-ifelse(tree@data$support<50,NA,tree@data$support)


library(RColorBrewer)
library(scales)
getPalette = colorRampPalette(c(brewer.pal(8, "Accent"),brewer.pal(12, "Paired"),brewer.pal(12, "Set3"),brewer.pal(12, "Paired") ))(44)

show_col(getPalette,labels=T)


ggtree(trs, layout="circular", 
       #open.angle=10, 
       size=0.5, 
       branch.length="none" 
)+
  # geom_tiplab(aes(color = Phylum),size=3, hjust= -0.06) +
  geom_tippoint(aes(color = Phylum), size=3) + 
  # geom_text2(aes(subset=!isTip, label=node), hjust=-0.3, size=2, color="deepskyblue4")+
  geom_nodepoint(color="#FF6633", size=1.5)+ 
  #theme_tree2() 
  scale_color_manual(values = getPalette)
