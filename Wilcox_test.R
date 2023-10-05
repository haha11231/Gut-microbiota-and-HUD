#Wilcox test

library(tidyverse)
library(data.table)
library(rstatix)


data <- read.csv("Wilcox.test.data.phylum.csv", header = T, row.names = 1)
  
df_l <- data %>% 
    pivot_longer(cols = 2:ncol(data), names_to = "变量", values_to = "积分") %>% 
    dplyr::mutate_if(is.character, as.factor)
  
 
K.test.wilcox <- df_l %>% group_by(变量) %>% wilcox_test(积分 ~ group) %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj")
