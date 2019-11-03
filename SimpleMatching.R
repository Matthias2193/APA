library(ggplot2)
library(caret)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('Evaluation Methods.R')
source('X Model Approach.R')

library("glmnet")


set.seed(1234)

#Data import
email <- read.csv('Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$conversion <- email$spend <- email$segment <- NULL

estimation_set <- email[1:20000,]
train <- email[20001:50000,]
test <- email[50001:64000,]

men_treatment_model <- glm(as.factor(estimation_set$men_treatment)~., family = binomial(link='logit'), 
                           data = estimation_set[,1:8]) 
women_treatment_model <- glm(as.factor(estimation_set$women_treatment)~., family = binomial(link='logit'), 
                             data = estimation_set[,1:8]) 

train$men_probability <- predict(men_treatment_model,train[,1:8])
train$women_probability <- predict(women_treatment_model,train[,1:8])

test$men_probability <- predict(men_treatment_model,test[,1:8])
test$women_probability <- predict(women_treatment_model,test[,1:8])

test$best_treatment <- "NULL"
test$men_prediction <- 0 
test$women_prediction <- 0
test$control_prediction <- 0

target <- "visit"

start_time <- Sys.time()
for(x in 1:nrow(test)){
  best_treatment <- "NULL" 
  men_predicttion <- 0
  women_prediction <- 0
  control_prediction <- 0
  distances <- c()
  distances <- 
    apply(train,MARGIN = 1, function(y) ((test[x,"men_probability"]-as.numeric(y["men_probability"]))^2+
                                         (test[x,"women_probability"]-as.numeric(y["women_probability"]))^2+
                                         (test[x,"control_probability"]-as.numeric(y["control_probability"]))^2))
  # for(y in 1:nrow(train)){
  #   temp_dist <- (test[x,"men_probability"] - train[y,"men_prbobility"])^2
  #   temp_dist <- temp_dist + (test[x,"women_probability"] - train[y,"women_probability"])^2
  #   temp_dist <- temp_dist + (test[x,"control_probability"] - train[y,"control_probability"])^2
  #   distances <- c(distances,temp_dist)
  # }
  temp_data <- train[order(distances),][1:100,]
  test[x,"men_prediction"] <- mean(temp_data[temp_data[,"men_treatment"] == 1,target])
  test[x,"women_prediction"] <- mean(temp_data[temp_data[,"women_treatment"] == 1,target])
  test[x,"control_prediction"] <- mean(temp_data[temp_data[,"control"] == 1,target])
}
print(start_time-Sys.time())