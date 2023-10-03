###Random forest


library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(tidyverse)
library(data.table)



df <- read.csv("all_OPLSDA_Genus.csv", header = T, row.names = 1)
df1 <- df %>% filter(p.adj<0.05)
df2 <- data.frame(t(df1))


df2$Sample <- row.names(df2)

df3 <- read.csv("Genus.all.group.csv", header = T)
df4 <- data.frame(df3[,1:2])
data1 <- merge(df4, df2, by = "Sample", all = F)
row.names(data1) <- data1$Sample
data <- data.frame(data1[,-1])

data=as.data.frame(lapply(data,as.numeric))
data$group <- as.factor(df4$group)
summary(data)
str(data)



data <- data


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

importance(rf.train)
varImpPlot(rf.train) 



rf.test <- predict(rf.train, 
                   newdata = testset,
                   type = "class") 
rf.test

t = table(rf.test, testset$group)
acc = sum(diag(t))/nrow(testset)*100


rf.cf <- caret::confusionMatrix(as.factor(rf.test), as.factor(testset$group)) 
rf.cf


