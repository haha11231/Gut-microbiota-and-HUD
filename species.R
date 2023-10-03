rm(list = ls()) 
getwd()
###设置工作文件夹
#setwd("/Users/a/Documents/adolscent RNA-seq/D/pvalue/NM_SMND")
#setwd("/Volumes/MAc KIOXIA Storage/R/adolscent RNA-seq/D/pvalue3_ok/NM_SMND")
setwd("/Volumes/MAc KIOXIA Storage/gut/16s/abundance stat")
dir()
library(tidyverse)
library(data.table)

#macbook12地址
#/Users/A/Desktop/16s/OPLSDA/addiction_withdrawal


# BiocManager::install('ropls')

library(ropls)

data <- read.csv("species_group.csv", header = T, row.names = 1)

# sample_info <- read.csv("sample_info.csv", header = T, row.names = 1)


###两组对比，OPLSDA，addiction vs withdrawal
# dataMatrix <- data[,2:157]

dataMatrix <- data[,2:ncol(data)]


sampleMetadata <- data[,1]

# genderFc = sampleMetadata[, "group"]

genderFc = sampleMetadata

# oplsda = opls(dataMatrix, genderFc, predI = 1, orthoI = NA)
plsda = opls(dataMatrix, genderFc)

vip <- getVipVn(plsda)
vip_select <- vip[vip > 1]    #通常以VIP值>1作为筛选标准
head(vip_select)


a <- t(dataMatrix)

vip_select <- cbind(a[names(vip_select), ], vip_select)
# colnames(vip_select)[71] <- 'VIP'
# names(vip_select)[71] <- 'VIP'
vip_select <- data.frame(vip_select)
vip_select <- vip_select[order(vip_select$vip_select, decreasing = TRUE), ]
head(vip_select)    #带注释的代谢物，VIP>1 筛选后，并按 VIP 降序排序


library(ggplot2)
# vip_select$cat = paste('A',1:nrow(vip_select), sep = '')

vip_select$cat = row.names(vip_select)
vip_select$cat <- factor(vip_select$cat, levels = vip_select$cat)

ggplot(vip_select, aes(factor(cat), vip_select)) +
  geom_segment(aes(x = cat, xend = cat,
                   y = 0, yend = vip_select)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = '', y = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())




ggplot(vip_select, aes(vip_select, factor(cat))) +
  geom_segment(aes(y = cat, yend = cat,
                   x = 0, xend = vip_select)) +
  geom_point(shape = 21, size = 5, color = '#008000' ,fill = '#008000') +
  geom_point(aes(1,2.5), color = 'white') +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = '', x = 'VIP value') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.text = element_text(color = 'black',size = 12, family = 'Arial', face = 'plain'),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'),
        axis.ticks.x = element_blank())





# df_l <- data %>% 
#   pivot_longer(cols = 2:157, names_to = "变量", values_to = "积分") %>% 
#   dplyr::mutate_if(is.character, as.factor)


df_l <- data %>% 
  pivot_longer(cols = 2:ncol(data), names_to = "变量", values_to = "积分") %>% 
  dplyr::mutate_if(is.character, as.factor)


str(df_l)

library(rstatix)

# K.test.dunn <- df_l %>% group_by(变量) %>% dunn_test(积分 ~ group)
# K.test.wilcox <- df_l %>% group_by(变量) %>% wilcox_test(积分 ~ group)
# K.test.wilcox <- K.test.wilcox[order(K.test.wilcox$p, decreasing = F), ]
# K.test.wilcox$p.adjust <- p.adjust(K.test.wilcox$p,method="fdr",n=length(K.test.wilcox$p))
# K.test.wilcox11 <- df_l %>% group_by(变量) %>% wilcox_test(积分 ~ group) %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj")


#秩和检验
# K.test.wilcox <- df_l %>% group_by(变量) %>% wilcox_test(积分 ~ group) %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj")
K.test.kruskal <- df_l %>% group_by(变量) %>% kruskal_test(积分 ~ group)  %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj")


df5 <- K.test.kruskal %>% filter(p.adj<0.05)

dir.create("species")

# write.table (vip_select, file ="./results/all_vip_select_Genus.csv", sep =",", row.names = F)
write.table (K.test.kruskal, file ="./species/all_kruskal_FDR_species.csv", sep =",", row.names = F)
# write.table (df4, file ="./results/all_OPLSDA_Genus.csv", sep =",", row.names = F)
write.table (df5, file ="./species/all_OPLSDA_species_select.csv", sep =",", row.names = F)

select <- as.character(df5$变量)

dunn.data <- data[,c("group", select)]


#dunn检验
# dunn.data <- read.csv("all_OPLSDA_Genus_dunn_data.csv", header = T, row.names = 1)

# df_l_dunn <- dunn.data %>% 
#   pivot_longer(cols = 2:94, names_to = "变量", values_to = "积分") %>% 
#   dplyr::mutate_if(is.character, as.factor)

df_l_dunn <- dunn.data %>% 
  pivot_longer(cols = 2:ncol(dunn.data), names_to = "变量", values_to = "积分") %>% 
  dplyr::mutate_if(is.character, as.factor)


str(df_l_dunn)

K.test.dunn <- df_l_dunn %>% group_by(变量) %>% dunn_test(积分 ~ group) %>% add_significance("p.adj")


library(tidyr)
K.test.dunn <- tidyr::unite(K.test.dunn, "comparison", group1, group2, sep = " vs. ", remove = FALSE) #对gtf_data中的group1与group2合并，以“ vs. ”连接，生成新的名为“comparison”的一列

write.table (K.test.dunn, file ="./species/all_OPLSDA_species_select_dunn.csv", sep =",", row.names = F)

##画图
# library(ggpubr)
# 
# pwc <- K.test.dunn %>% add_xy_position(x = "变量")
# 
# bxp <- ggboxplot(
#   df_l_dunn, x = "变量", y = "积分",
#   color = "group", palette = "jco",
#   # add = "jitter",
#   # facet.by = "变量" #按照变量分不同的面板
# )
# bxp
# 
# 
# bxp +
#   stat_pvalue_manual(pwc, hide.ns = TRUE) +
#   labs(
#     subtitle = get_test_label(K.test.dunn, detailed = TRUE),
#     caption = get_pwc_label(pwc),
#   )



###统计增加均数
summary_stats <- df_l %>% group_by(group, 变量) %>% get_summary_stats(积分, type = "mean_sd")
summary_stats <- summary_stats[,-3]

# b <- summary_stats[1:156,]
# c <- summary_stats[157:312,]
# bb <- summary_stats[313:468,]
# cc <- summary_stats[469:624,]


b <- summary_stats[1:(nrow(summary_stats)/4),]
c <- summary_stats[(nrow(summary_stats)/4+1):(nrow(summary_stats)/4*2),]
bb <- summary_stats[(nrow(summary_stats)/4*2+1):(nrow(summary_stats)/4*3),]
cc <- summary_stats[(nrow(summary_stats)/4*3+1):nrow(summary_stats),]

d1 <- merge(b, c, by = "变量", all = F)
d2 <- merge(bb, cc, by = "变量", all = F)
d <- merge(d1, d2, by = "变量", all = F)

K.test.kruskal.summary <- merge(d, K.test.kruskal, by.x = "变量", by.y = "变量", all = F)

write.table (K.test.kruskal.summary, file ="./species/all_OPLSDA_species_All_summary.csv", sep =",", row.names = F)


select1 <- data.frame(df5$变量)
df8 <- merge(select1, K.test.kruskal.summary, by.x = "df5.变量", by.y = "变量", all = F)

df8$logKrus <- -log10(df8$p.adj)

colnames(df8)[c(4,8,12,16)] <- c("Addiction", "Control", "Methadone", "Withdrawal")

df8$Addictionvs.Control <- log2(df8$Addiction/df8$Control)

df8$Withdrawalvs.Control <- log2(df8$Withdrawal/df8$Control)

df8$Methadonevs.Control <- log2(df8$Methadone/df8$Control)

df8$Addictionvs.Withdrawal <- log2(df8$Addiction/df8$Withdrawal)

df8$Addictionvs.Methadone <- log2(df8$Addiction/df8$Methadone)

df8$Withdrawalvs.Methadone <- log2(df8$Withdrawal/df8$Methadone)


write.table (df8, file ="./species/all_OPLSDA_species_select_summary.csv", sep =",", row.names = F)

df9 <- df8[,c(1,4,8,12,16, 27:32)]

colnames(df9)[6:11] <- c("Addiction vs. Control", "Control vs. Withdrawal", "Control vs. Methadone", "Addiction vs. Withdrawal", "Addiction vs. Methadone", "Methadone vs. Withdrawal")

###平均值
mean.data <- df9[,c(1,3,2,5,4)]
row.names(mean.data) <- mean.data$df5.变量
mean.data <- mean.data[-1]


df11 <- gather(df9, comparison, log2FC, 6:11)
df12 <- tidyr::unite(df11, "Genus_comparison", df5.变量, comparison, sep = "_", remove = FALSE)
K.test.dunn1 <- tidyr::unite(K.test.dunn, "Genus_comparison", 变量, comparison, sep = "_", remove = FALSE) #对gtf_data中的group1与group2合并，以“ vs. ”连接，生成新的名为“comparison”的一列
df13 <- merge(df12, K.test.dunn1, by.x = "Genus_comparison", by.y = "Genus_comparison", all = F)
df13$logDunn <- -log10(df13$p.adj)
write.table (df13, file ="./species/all_OPLSDA_species_select_dunn_summary.csv", sep =",", row.names = F)

library(pheatmap)
pheatmap(mean.data, 
         scale = "row", 
         # clustering_method = "average",  
         cluster_rows = F, cluster_cols=FALSE, 
         # gaps_row = 20, 
         # cellheight = 10, cellwidth = 10,
         # annotation_row = annotation_row, 
         # annotation_col = annotation_col, 
         # annotation_colors = ann_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         # breaks = c(seq(-2.4,1.6,by=0.4)),
         show_rownames=T)




ggplot(df13,
       aes(y = 变量, x = comparison.x))+
  geom_tile(aes(fill = log2FC), colour = 'white', size = 0.2)+
  # scale_fill_gradient2('deepskyblue1')+
  scale_fill_gradient2(low = 'deepskyblue1', high = 'deeppink2', mid = 'white')+
  # scale_fill_gradientn(colours = c('deepskyblue1','white','deeppink2'))+
  # scale_fill_gradientn(colours = colorRampPalette(colors = c("#2E9FDF", "white", "#FC4E07"))(10))+
  geom_text(aes(label = p.adj.signif), col = 'black', size = 6, vjust = 0.8)+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.background=element_rect(fill="white",color="black"), 
        axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(df13,
       aes(y = 变量, x = comparison.x))+
  geom_tile(aes(fill = log2FC), colour = 'white', size = 0.2)+
  # scale_fill_gradient2('deepskyblue1')+
  scale_fill_gradient2(low = 'darkgreen', high = 'orange', mid = 'white')+
  # scale_fill_gradientn(colours = c('deepskyblue1','white','deeppink2'))+
  # scale_fill_gradientn(colours = colorRampPalette(colors = c("#2E9FDF", "white", "#FC4E07"))(10))+
  geom_text(aes(label = p.adj.signif), col = 'black', size = 6, vjust = 0.8)+
  theme(panel.grid.major.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.background=element_rect(fill="white",color="black"), 
        axis.text.x = element_text(angle = 45, hjust = 1))



ggplot(df13,
       aes(y = factor(变量), x = comparison.x))+
  geom_point(aes(color = log2FC, size = logDunn))+
  theme_bw()+
  #theme(panel.grid =element_blank())+ #去掉网格
  #theme(panel.border = element_blank())+ #去掉外框
  #theme(axis.line = element_line(size=0.5, color = "black"))+ #添加坐标轴
  theme(axis.text.x=element_text(size=12, color = "black", 
                                 vjust = 0.5, 
                                 hjust = 1, 
                                 angle = 90),
        axis.text.y=element_text(size=12, color = "black"))+ #调整坐标轴字体大小、颜色
  theme(axis.title.x=element_blank())+
  theme(legend.key.height = unit(0.4,'cm'), legend.key.width = unit(0.4,'cm'),
        #legend.spacing = unit(2,'cm'),
        legend.text = element_text(colour = 'black',size = 12))+
  scale_color_gradientn(colours = c('blue','white','red'))



