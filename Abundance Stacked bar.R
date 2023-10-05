#abundance Stacked bar



library(ggplot2)
library(ggprism)
library(reshape)
library(ggalluvial)
library(ggsci)
library(ggthemes)
library(RColorBrewer)
library(scales)


df <- read.csv("Abundance.Stacked.test.phylum.csv", header = T)
names(df)[1:5] <- c("X", "Control", "Addiction","Withdrawal", "Methadone")
summary(df)


df1<-melt(df,id.vars='X')
names(df1)[1:2] <- c("X","group")  



getPalette = colorRampPalette(c(brewer.pal(9, "Set1"),brewer.pal(12, "Set3"),brewer.pal(8, "Accent"),brewer.pal(12, "Paired"), brewer.pal(8, "Dark2") ))(49)
show_col(getPalette,labels=T)

levels1 <- factor(df$X) 

ggplot(df1,
       aes(x = group, y = value,
           stratum = factor(X, levels = levels1),
           alluvium = factor(X, levels = levels1),
           fill = factor(X, levels = levels1),
           label = factor(X, levels = levels1))) +
  scale_fill_manual(values = getPalette) +
  geom_flow(stat = "alluvium", 
            lode.guidance = "frontback",
            #color = "darkgray",
            #alpha = .2,
            curve_type = "quintic",
            width = 0.55) +
  geom_stratum(width = 0.55, color = "black", size = 0.3) +
  #ggtitle("student curricula across several semesters")+
  labs(y="Relative Abundance (%)", size = 8, color = "black")+
  scale_y_continuous(expand = c(0,0))+
  #theme(legend.position = "bottom")+
  theme_bw()+
  theme(panel.grid =element_blank())+ 
  theme(panel.border = element_blank())+ 
  theme(axis.line = element_line(linewidth = 0.5, color = "black"))+ 
  theme(axis.text.x=element_text(size=8, color = "black", vjust = 1, hjust = 1, angle = 45),
        axis.text.y=element_text(size=8, color = "black"))+ 
  theme(axis.title.x=element_blank())+
  theme(legend.key.height = unit(0.4,'cm'), legend.key.width = unit(0.4,'cm'),
        #legend.spacing = unit(2,'cm'),
        legend.text = element_text(colour = 'black',size = 12))+
  guides(fill = guide_legend(title = "Phylum", ncol = 1, byrow = T, reverse = F))
