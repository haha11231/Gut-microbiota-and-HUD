###PLS-DA


library(tidyverse)
library(data.table)
library(ropls)

data <- read.csv("PLSDA.test.data.csv", header = T, row.names = 1)





dataMatrix <- data[,2:93]
sampleMetadata <- data[,1]



genderFc = sampleMetadata


plsda = opls(dataMatrix, genderFc)

vip <- getVipVn(plsda)
vip_select <- vip[vip > 1]    
head(vip_select)


a <- t(dataMatrix)

vip_select <- cbind(a[names(vip_select), ], vip_select)

vip_select <- data.frame(vip_select)
vip_select <- vip_select[order(vip_select$vip_select, decreasing = TRUE), ]
head(vip_select)    


library(ggplot2)

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




sample.score = plsda@scoreMN %>%  
  as.data.frame() %>%
  mutate(group = sampleMetadata,
         o1 = plsda@orthoScoreMN[,1]) 

Sample_Name <- rownames(sample.score) 

hull_group <- sample.score %>%
  dplyr::mutate(Sample_Name = Sample_Name) %>%
  dplyr::group_by(sample.score$group) %>%
  dplyr::slice(chull(p1, p2))


ggplot(sample.score, aes(p1, p2, color = group)) +
  geom_point(size = 1.3)+
  ggplot2::geom_polygon(data = hull_group, aes(fill = group), alpha = 0.1, linewidth = 0.5)+
  labs(x=paste0("Orthogonal T socre 1"),
       y=paste0("Orthogonal T socre 2"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))+
  scale_color_manual(values = c("#cc0033", "#00e990", "#ff9933", "#0066cc"))+
  scale_fill_manual(values = c("#cc0033", "#00e990", "#ff9933", "#0066cc"))+
  geom_vline(xintercept = 0,lty="dashed", colour = '#444444', linewidth = 0.5)+
  geom_hline(yintercept = 0,lty="dashed", colour = '#444444', linewidth = 0.5)



