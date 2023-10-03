###UpSet plot

library(tidyverse)
library(data.table)


data <- read.csv("pathway2.csv", header = T)
data <- data.frame(data)
colnames(data) <- c("ID", "Loc", "y")
data1 <- dcast(data=data, ID ~ Loc)


for (i in 1:ncol(data1)){
  data1[,i][is.na(data1[,i])] <- 0
}

data2 <- select(data1, 1,2,5,7,4,3,6)

data3 <- read.csv("data2.csv", header = T, row.names = 1)

library(UpSetR)

colnames(data3) <- c("Addiciton.vs.Control", "Methadone.vs.Control", "Withdrawal.vs.Control", "Addiciton.vs.Withdrawal", "Addiciton.vs.Methadone", "Withdrawal.vs.Methadone")


upset(data3,
      #nsets = 6,
      sets = c("Withdrawal.vs.Methadone", "Addiciton.vs.Methadone", "Addiciton.vs.Withdrawal", "Methadone.vs.Control", "Withdrawal.vs.Control", "Addiciton.vs.Control"),
      mb.ratio = c(0.65, 0.35), ###上下比例
      number.angles = 0, ###设置柱状图柱子上方数字的倾斜角度
      point.size = 6, #点大小
      line.size = 0.5, #线粗细
      mainbar.y.label = 'Pathway inersections (#)',
      sets.x.label = "Pathway per group (#)",
      text.scale = 2,
      order.by = c("degree", "freq"),
      decreasing = c(T, T),
      keep.order = T, #设置sets才可用此行参数，sets设置还需先设置data列名
      #empty.intersections = "on", #显示交集为空的组合
      sets.bar.color = 'darkgreen',
      main.bar.color = 'darkblue',
      matrix.color = 'darkred',
      mainbar.y.max = 33,
      matrix.dot.alpha = 0.5,
      set_size.show = T
      )