###PcOA


library(tidyverse)
library(data.table)
library(GUniFrac)
library(ape)
library(ade4)
library(ggplot2)



sample_info <- read.csv('sample_info.csv', row.names=1)
data <- read.csv('otu1.csv', header = T, row.names=1)
group <- sample_info$condition



###PCoA bray curtis
library(vegan)
df <- data
bray_dist<-vegdist(df,method = "bray")

df.pcoa<-cmdscale(bray_dist,k=(nrow(df)-1),eig=TRUE)
df.plot<-data.frame(df.pcoa$points)

library(plotly)
library(dplyr)
sum_eig <- sum(df.pcoa$eig)
eig_percent <- round(df.pcoa$eig/sum_eig*100,1)
x_label<-eig_percent[1]
y_label<-eig_percent[2]

Group <- group
df_pcs <- data.frame(df.pcoa$points, Group=Group)
Sample_Name <- rownames(df.pcoa) 

hull_group <- df_pcs %>%
  dplyr::mutate(Sample_Name = Sample_Name) %>%
  dplyr::group_by(Group) %>%
  dplyr::slice(chull(X1, X2))

ggplot(df.plot, aes(x = X1, y = X2, color = Group)) +
  geom_point(size = 1)+
  ggplot2::geom_polygon(data = hull_group, aes(fill = Group), alpha = 0.1, size = 0.5)+
  labs(x=paste0("PCoA1 (",x_label,"%)"),
       y=paste0("PCoA2 (",y_label,"%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  scale_color_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  scale_fill_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  geom_vline(xintercept = 0,lty="dashed", colour = '#444444', size = 0.5)+
  geom_hline(yintercept = 0,lty="dashed", colour = '#444444', size = 0.5)

  




set.seed(30)
bray.div <- adonis2(df ~ condition, data = sample_info, permutations = 999, method="bray")



library(pairwiseAdonis)
bray.pairwise.adonis <- pairwise.adonis(x=df, factors=sample_info$condition, sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)







###weighted UniFrac
library(GUniFrac)
library(phyloseq)

tree  <- read_tree("/Volumes/MAc KIOXIA Storage/gut/16s/a/OTU_final_phylogeny_tree")
otu.tab <- read.csv('otu1.csv', header = T, row.names=1)
sample_info <- read.csv('sample_info.csv', header = T, row.names=1)
groups <- factor(sample_info$condition)


comm <- otu.tab
ape::is.binary(tree)
tree <- ape::multi2di(tree)

otu.tab.rff <- Rarefy(otu.tab)$otu.tab.rff
unifracs <- GUniFrac(otu.tab.rff, tree, alpha=c(0, 0.5, 1))$unifracs 

dw <- unifracs[, , "d_1"] # Weighted UniFrac
du <- unifracs[, , "d_UW"] # Unweighted UniFrac
dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
d0 <- unifracs[, , "d_0"] # GUniFrac with alpha 0
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5



df.pcoa1<-cmdscale(dw,k=(nrow(otu.tab)-1),eig=TRUE) 
df.plot1<-data.frame(df.pcoa1$points)

sum_eig <- sum(df.pcoa1$eig)
eig_percent <- round(df.pcoa1$eig/sum_eig*100,1)
x_label<-eig_percent[1]
y_label<-eig_percent[2]

Group <- group
df_pcs1 <- data.frame(df.pcoa1$points, Group=Group)
Sample_Name <- rownames(df.pcoa1) 

hull_group <- df_pcs1 %>%
  dplyr::mutate(Sample_Name = Sample_Name) %>%
  dplyr::group_by(Group) %>%
  dplyr::slice(chull(X1, X2))


ggplot(df.plot1, aes(x = X1, y = X2, color = Group)) +
  geom_point(size = 1.3)+
  ggplot2::geom_polygon(data = hull_group, aes(fill = Group), alpha = 0.1, size = 0.5)+
  labs(x=paste0("PCoA1 (",x_label,"%)"),
       y=paste0("PCoA2 (",y_label,"%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  scale_color_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  scale_fill_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  geom_vline(xintercept = 0,lty="dashed", colour = '#444444', size = 0.5)+
  geom_hline(yintercept = 0,lty="dashed", colour = '#444444', size = 0.5)





set.seed(30)
wUniFrac.div <- adonis2(as.dist(dw) ~ groups)


wUniFrac.pairwise.adonis <- pairwise.adonis(x=as.dist(dw), factors=sample_info$condition, 
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)







###unweighted UniFrac
df.pcoa2<-cmdscale(du,k=(nrow(otu.tab)-1),eig=TRUE) 
df.plot2<-data.frame(df.pcoa2$points)

sum_eig <- sum(df.pcoa2$eig)
eig_percent <- round(df.pcoa2$eig/sum_eig*100,1)
x_label<-eig_percent[1]
y_label<-eig_percent[2]

Group <- group
df_pcs2 <- data.frame(df.pcoa2$points, Group=Group)
Sample_Name <- rownames(df.pcoa2) 

hull_group <- df_pcs2 %>%
  dplyr::mutate(Sample_Name = Sample_Name) %>%
  dplyr::group_by(Group) %>%
  dplyr::slice(chull(X1, X2))



ggplot(df.plot2, aes(x = X1, y = X2, color = Group)) +
  geom_point(size = 1.3)+
  ggplot2::geom_polygon(data = hull_group, aes(fill = Group), alpha = 0.1, size = 0.5)+
  labs(x=paste0("PCoA1 (",x_label,"%)"),
       y=paste0("PCoA2 (",y_label,"%)"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.5),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  scale_color_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  scale_fill_manual(values = c("#CC0033", "#999999", "#FF9933", "#0066CC"))+
  geom_vline(xintercept = 0,lty="dashed", colour = '#444444', size = 0.5)+
  geom_hline(yintercept = 0,lty="dashed", colour = '#444444', size = 0.5)




set.seed(30)

uUniFrac.div <- adonis2(as.dist(du) ~ groups)

library(pairwiseAdonis)
uUniFrac.pairwise.adonis <- pairwise.adonis(x=as.dist(du), factors=sample_info$condition, 
                                            p.adjust.m = "BH",
                                            reduce = NULL,
                                            perm = 999)