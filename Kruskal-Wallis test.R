###Kruskal-Wallis test



library(tidyverse)
library(data.table)


data <- read.csv("Kruskal.Wallis.test.data.phylum.csv", header = T, row.names = 1)


df_l <- data %>% 
  pivot_longer(cols = 2:ncol(data), names_to = "变量", values_to = "积分") %>% 
  dplyr::mutate_if(is.character, as.factor)


str(df_l)

library(rstatix)

K.test.kruskal <- df_l %>% group_by(变量) %>% kruskal_test(积分 ~ group)  %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj")


df5 <- K.test.kruskal %>% filter(p.adj<0.05)

select <- as.character(df5$变量)

dunn.data <- data[,c("group", select)]



df_l_dunn <- dunn.data %>% 
  pivot_longer(cols = 2:ncol(dunn.data), names_to = "变量", values_to = "积分") %>% 
  dplyr::mutate_if(is.character, as.factor)


str(df_l_dunn)

K.test.dunn <- df_l_dunn %>% group_by(变量) %>% dunn_test(积分 ~ group) %>% add_significance("p.adj")


library(tidyr)
K.test.dunn <- tidyr::unite(K.test.dunn, "comparison", group1, group2, sep = " vs. ", remove = FALSE) 

summary_stats <- df_l %>% group_by(group, 变量) %>% get_summary_stats(积分, type = "mean_sd")
summary_stats <- summary_stats[,-3]




b <- summary_stats[1:(nrow(summary_stats)/4),]
c <- summary_stats[(nrow(summary_stats)/4+1):(nrow(summary_stats)/4*2),]
bb <- summary_stats[(nrow(summary_stats)/4*2+1):(nrow(summary_stats)/4*3),]
cc <- summary_stats[(nrow(summary_stats)/4*3+1):nrow(summary_stats),]

d1 <- merge(b, c, by = "变量", all = F)
d2 <- merge(bb, cc, by = "变量", all = F)
d <- merge(d1, d2, by = "变量", all = F)

K.test.kruskal.summary <- merge(d, K.test.kruskal, by.x = "变量", by.y = "变量", all = F)

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




df9 <- df8[,c(1,4,8,12,16, 27:32)]

colnames(df9)[6:11] <- c("Addiction vs. Control", "Control vs. Withdrawal", "Control vs. Methadone", "Addiction vs. Withdrawal", "Addiction vs. Methadone", "Methadone vs. Withdrawal")


mean.data <- df9[,c(1,3,2,5,4)]
row.names(mean.data) <- mean.data$df5.变量
mean.data <- mean.data[-1]


df11 <- gather(df9, comparison, log2FC, 6:11)
df12 <- tidyr::unite(df11, "Genus_comparison", df5.变量, comparison, sep = "_", remove = FALSE)
K.test.dunn1 <- tidyr::unite(K.test.dunn, "Genus_comparison", 变量, comparison, sep = "_", remove = FALSE) 
df13 <- merge(df12, K.test.dunn1, by.x = "Genus_comparison", by.y = "Genus_comparison", all = F)
df13$logDunn <- -log10(df13$p.adj)


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
