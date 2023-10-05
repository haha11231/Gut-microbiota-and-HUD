###Correlation analysis




library(ppcor)
data1 <- read.csv("Correlation.test.data.abundance.csv", header = T, row.names = 1)
data2 <- read.csv("Correlation.test.data.phenotype.csv", row.names = 1, header = T)


data3 <- data2


correlation1 <- function(class){
  y <- data3[,class]
  colnames <- colnames(data1)
  do.call(rbind, lapply(colnames, function(x){
    dd <- pcor.test(data1[,x], y, data3[ , !colnames(data3) %in% c(class)], method = 'spearman')
    data.frame(class = class, cellmarker = x, col = dd$estimate, p.value = dd$p.value)
  }))
}
class <- colnames(data3)
datax <- do.call(rbind, lapply(class, correlation1))

data <- datax

data$significant = ''
data$significant[which( (data$p.value < 0.05) & (data$p.value > 0.01) )] = '*'
data$significant[which( (data$p.value < 0.01) & (data$p.value > 0.001) )] = '**'
data$significant[which( (data$p.value < 0.001) & (data$p.value > 0.0001) )] = '***'
data$significant[which (data$p.value < 0.0001)] = '****'



library(ggplot2)

data$cellmarker <- factor(data$cellmarker, levels = c("Coprococcus", "Intestinibacter", "Flavonifractor", "Escherichia", "Clostridium_XlVa", "Turicibacter", "Terrisporobacter", "Megasphaera", "Dorea", "Actinomyces", "Butyricicoccus", "Paraprevotella", "Anaerostipes", "Streptococcus", "Megamonas", "Clostridium_XlVb", "Faecalibacterium", "Ruminococcus2", "Lactobacillus", "Romboutsia", "Roseburia", "Blautia", "Fusicatenibacter", "Weissella"))
data$class <- factor(data$class, levels = c("Age", "Education.level", "Income", "Spicy.foods", "Drink.tea", "Martital.status", "Employment", "Route.of.heroin.use", "Average.dose", "Age.of.onset", "Heroin.use.duration", "HAMD.score", "MTT.duration"))

ggplot(data, aes(class, cellmarker))+
  geom_tile(aes(fill = col), colour = 'white', linewidth = 1)+
  # scale_fill_gradient2(low = 'darkblue', mid = 'white', high = '#CC0000')+ #heroin
  # scale_fill_gradient2(low = 'deepskyblue2', mid = 'white', high = 'hotpink2')+ #withdrawal
  scale_fill_gradient2(low = 'seagreen', mid = 'white', high = 'orange')+ #methadone
  geom_text(aes(label = significant), col = 'black', size = 6, vjust = 0.8)+
  theme(panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(), 
        panel.background=element_rect(fill="white",color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))







library(corrplot)
library(ggcorrplot)

d3 <- data2


col <- cor(d3, method = "spearman")
col.p <- cor_pmat(d3, method = "spearman")
col.ci <- cor.mtest(d3, conf.level = 0.95)

library(RColorBrewer)
corrplot(corr = col, p.mat = col.p, method = "square", type = "full",
         #order = "hclust", hclust.method = 'average',  addrect = 4,
         col = colorRampPalette(c("seagreen", "white", "orange"))(100), #methadone
         #tl.pos="td",
         tl.col = 'black',
         #insig="blank",
         insig = 'label_sig',
         sig.level = c(.0001, .001, .01, .05), pch.cex = 1.5,
         #tl.srt = 45,
         diag = FALSE,
         cl.ratio = 0.07, cl.align = "l"
)


