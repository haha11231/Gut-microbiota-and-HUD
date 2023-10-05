###Co-occurrence networks


library(ggClusterNet)
library(phyloseq)
library(sna)
library(tidyverse)
library(igraph)
library(ggrepel)
library(RColorBrewer)
library(scales)


getPalette <- colorRampPalette(c(brewer.pal(12, "Set3")))(24)
show_col(getPalette,labels=T)


metadata = read.csv("Co.occurrence.group.csv", row.names = 1)
datag = read.csv("Co.occurrence.test.data.control.csv", row.names=1)
taxonomy = read.csv("Co.occurrence.taxonomy.csv", row.names=1, header = T)




library(corrplot)
library(ggcorrplot)
d3 <- t(datag)
col <- cor(d3, method = "spearman")
col.p <- cor_pmat(d3, method = "spearman")
col.ci <- cor.mtest(d3, conf.level = 0.95)


ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(datag), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))#,
              # phy_tree(tree),
              # refseq(rep)
)


result = corMicro(ps = ps,
                  N = 0,
                  # method.scale = "RLC",
                  r.threshold=0.4, 
                  p.threshold=0.05, 
                  method = "spearman"
)




cor = result[[1]]
dim(cor)

ps_net = result[[3]]
otu_table = as.data.frame(t(vegan_otu(ps_net)))
tax_table = as.data.frame(vegan_tax(ps_net))



result2 <- model_Gephi.2(cor = cor,
                         method = "cluster_fast_greedy",
                         seed = 12
)


node = result2[[1]]



nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
head(nodes)



edge = edgeBuild(cor = cor,node = node)
head(edge)



colnames(edge)[8] = "cor"
head(nodes)




library(ggalt)
library(ggnewscale)
netClu = result2[[3]]
head(netClu)

row.names(netClu) = netClu$ID
nodeG = merge(nodes,netClu,by = "row.names",all  =FALSE)
dim(nodeG)
head(nodeG)
nodeG$group <- factor(nodeG$group) 




ggplot() + geom_curve(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = cor,curvature = 0.2),
                            data = edge, size = 0.5,alpha = 1) +
  geom_point(aes(X1, X2, fill = elements, size = mean), pch = 21, data = nodes) +
  geom_text(aes(X1, X2,label = elements),pch = 21, data = nodes) +
  scale_colour_manual(values = c("#C3D9EA","#F7BBBB")) +
  scale_fill_manual(values = getPalette, na.value = NA)+
  scale_size(range = c(15, 30)) +
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.background = element_rect(colour = NA)) +
  theme(panel.background = element_rect(fill = "white",  colour = NA)) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
