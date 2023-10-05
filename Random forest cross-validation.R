###Random forest cross-validation


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
df4 <- data.frame(df3[,c(1,2)])
colnames(df4) <- c("Sample", "group")
data1 <- merge(df4, df2, by = "Sample", all = F)
row.names(data1) <- data1$Sample
data <- data.frame(data1[,-1])



data=as.data.frame(lapply(data,as.numeric))
data$group <- as.factor(df4$group)
summary(data)
str(data)


trainlist <- caret::createDataPartition(data$group, 
                                        times = 1,
                                        p = 1,
                                        list = TRUE)

trainset <- data[trainlist$Resample1,]

rf.train <- randomForest(as.factor(group) ~ ., 
                         data = trainset,
                         importance = TRUE, 
                         ntree = 2000,
                         na.action = na.pass, 
                         keep.inbag = TRUE) 

rf.train

importance(rf.train) 
varImpPlot(rf.train) 


model_tuned <- tuneRF(
  x=trainset[-1], #define predictor variables
  y=trainset$group, #define response variable
  ntreeTry=2000,
  mtryStart=4, 
  stepFactor=1.5,
  improve=0.01,
  trace=FALSE #don't show real-time progress
)


otu_train.cv <- replicate(5, rfcv(trainset[-1], trainset$group, cv.fold = 10, 
                                  # scale="log", 
                                  step = 1.5), simplify = FALSE)

otu_train.cv

error.cv <- sapply(otu_train.cv, "[[","error.cv")

matplot(otu_train.cv[[1]]$n.var,
        cbind(rowMeans(error.cv),error.cv),
        type = "l",
        lend = par("lend"),
        #lwd = c(2, rep(1, ncol(error.cv))),
        #col = 1:6, lty = 1:6, 
        #log = "y",
        xlab = "Number of variables", ylab = "Cross-validation error"
)


otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)


library(ggplot2)

ggplot(otu_train.cv.mean, aes(Group.1, x)) +
  geom_line() +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of variables', y = 'Cross-validation error')
