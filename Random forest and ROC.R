###Random forest and ROC


library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(tidyverse)
library(data.table)


df <- read.csv("Randomforest.test.data.csv", header = T, row.names = 1)

df2 <- data.frame(t(df))

df2$Sample <- row.names(df2)

df3 <- read.csv("Randomforest.test.data.group.csv", header = T)
df4 <- data.frame(df3[,c(1,3)])
colnames(df4) <- c("Sample", "group")
data1 <- merge(df4, df2, by = "Sample", all = F)
row.names(data1) <- data1$Sample
data <- data.frame(data1[,-1])

data=as.data.frame(lapply(data,as.numeric))
data$group <- df4$group
summary(data)



set.seed(6666)

trainlist <- caret::createDataPartition(data$group, 
                                        times = 1, 
                                        p = 0.6,  
                                        list = TRUE)

trainset <- data[trainlist$Resample1,]
testset <- data[-trainlist$Resample1,]




rf.train <- randomForest(as.factor(group) ~ ., 
                         data = trainset,
                         importance = TRUE, 
                         na.action = na.pass) 
rf.train

plot(rf.train, main = "randomforest origin")



rf.test <- predict(rf.train, 
                   newdata = testset,
                   type = "class") 
rf.test

t = table(rf.test, testset$group)
acc = sum(diag(t))/nrow(testset)*100
print(paste("模型准确率：", round(acc,4),sep = ""))


rf.cf <- caret::confusionMatrix(as.factor(rf.test), as.factor(testset$group)) 
rf.cf



rf.test2 <- predict(rf.train, 
                    newdata = testset,  
                    type = "prob") 

roc.rf <- roc(testset$group, rf.test2[,1])

auc(roc.rf)




pred <- ROCR::prediction(rf.test2[,2], testset$group)
perf <- performance(pred, "tpr", "fpr")
perf


x <- unlist(perf@x.values)
y <- unlist(perf@y.values)
plotdata <- data.frame(x,y)
names(plotdata) <- c("x", "y")

ggplot(plotdata) +
  geom_path(aes(x = x,
                y = y),
            colour = "#00416a",
            linewidth = 1) +
  labs(x = "False positive rate", y = "True positive rate")+
  # scale_color_gradient(name = "False positive rate",
  #                      low = "blue",
  #                      high = "red")+
  geom_segment(aes(y = 0, yend = 1, x = 0, xend = 1), color = "grey50", linetype = "dashed")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
