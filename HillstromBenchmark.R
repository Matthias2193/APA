# This script evaluates several different models on the Hillstrom data set which is contained in this repository.
# The models evaluated are Contextual Treatment Selection (CTS), Causal Forest, Separate Model Approach with 
# Random Forest and a custom Tree/Random Forest with two different gain functions ("Simple" and "Frac). More 
# information about the custom Tree and Random Forest can be found under ModelImplementations/DecisionTreeImplementation.R
# The user can specify the parameter n_predictions. If it is set to one, each model is trained once on the 
# original data set and then evaluated. If n_predictions is greater than 1, there will be n_predictions iterations.
# For each iteration a bootstrap sample is taken as the new data and then the models are build. After the models
# have been built and the predictions have been made, the predictions are evaluated and the results ploted.

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)

source('ModelImplementations/DecisionTreeImplementation.R')
source('ModelImplementations/RzepakowskiTree.R')
source('Evaluation Methods.R')
source('ModelImplementations/CausalTree.R')
source('ModelImplementations/Separate Model Approach.R')
source('ModelImplementations/ContextualTreatmentSelection.R')
source('ModelImplementations/VisualizationHelper.R')
source("ModelImplementations/PredictionFunctions.R")


set.seed(1234)
remain_cores <- 1
#Data import and preprocessing
email <- read.csv('Data/Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)
email$temp_feature1 <- rnorm(nrow(email),0,1)
email$temp_feature2 <- rbinom(n=nrow(email), size=1, prob=0.5)
email$temp_feature2 <- as.factor(email$temp_feature2)
email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'
control <- "control"
treatments <- c('men_treatment','women_treatment')
temp_features <- c("","temp_feature1","temp_feature2")
size_vector <- c(50000,100000,200000)
treatment_list <- c()
original_email <- email
keep_cols <- c("recency","history_segment","history","mens","womens","zip_code",
               "newbie","channel","spend","control")
feature_cols <- c("recency","history_segment","history","mens","womens","zip_code",
                  "newbie","channel")
dom_time <- c()
cts_time <- c()
sma_time <- c()
causal_time <- c()

for(treatment in treatments){
  keep_cols <- c(keep_cols,treatment)
  for(f in temp_features){
    if(f != ""){
      remain_cols <- c(keep_cols,f)
      new_feature_cols <- c(feature_cols,f)
    } else{
      new_feature_cols <- feature_cols
    }
    for(size in size_vector){
      email <- original_email[,remain_cols]
      new_email <- sample_n(email,size = size,replace = T)
      train <- new_email
      test <- new_email[1,]
      est_list <- set_up_tests(train[,new_feature_cols],TRUE, max_cases = 10)
      #DOM
      start_time <- Sys.time()
      forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 500,n_features = 3,
                                             criterion = "frac", remain_cores = remain_cores)
      dom_time <- c(dom_time,difftime(Sys.time(), start_time, units='mins'))
      #CTS
      start_time <- Sys.time()
      cts_forest <- build_cts(response, control, treatment_list, train, 500, nrow(train), 5, 2, 100, parallel = TRUE,
                              remain_cores = remain_cores)
      cts_time <- c(cts_time,difftime(Sys.time(), start_time, units='mins'))
      #SMA
      start_time <- Sys.time()
      pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf", mtry = 3, ntree = 300)
      sma_time <- c(sma_time,difftime(Sys.time(), start_time, units='mins'))
      #Causal Forest
      start_time <- Sys.time()
      causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,ntree = 1000,
                                                    s_rule = "TOT", s_true = T)
      causal_time <- c(causal_time,difftime(Sys.time(), start_time, units='mins'))
    }
  }
}



