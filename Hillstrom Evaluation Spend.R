#Test Script

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('ModelImplementations/DecisionTreeImplementation.R')
source('ModelImplementations/RzepakowskiTree.R')
source('Evaluation Methods.R')
source('ModelImplementations/CausalTree.R')
source('ModelImplementations/Causal Forest.R')
source('ModelImplementations/Separate Model Approach.R')
source('ModelImplementations/ContextualTreatmentSelection.R')
source('ModelImplementations/VisualizationHelper.R')
source("ModelImplementations/PredictionFunctions.R")


set.seed(1234)

#Data import
email <- read.csv('Data/Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'
control <- "control"
treatment_list <- c('men_treatment','women_treatment')


#email$spend <- email$spend/max(email$spend)

folds <- createFolds(email$spend, k = 5, list = TRUE, returnTrain = FALSE)

original_email <- email

for(f in 1:1){
  
  email <- original_email[sample(nrow(original_email),nrow(original_email),replace = TRUE),]
  idx <- createDataPartition(y = email[ , response], p=0.2, list = FALSE)
  train <- email[-idx, ]
  
  test <- email[idx, ]
  
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  
  val <- train[p_idx,]
  train_val <- train[-p_idx,]
  
  
  start_time <- Sys.time()
  for(c in c("simple","max","frac")){
    print(c)
    # Single Tree
    raw_tree <- build_tree(train_val,0,100,treatment_list,response,control,test_list,criterion = c)
    pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control,criterion = c)
    pred <- predict.dt.as.df(pruned_tree, test)
    write.csv(pred, paste("Predictions/Hillstrom/tree_",c,as.character(f),".csv",sep = ""), row.names = FALSE)


    #Forest
    forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3,
                                    pruning = F, criterion = c)
    pred <- predict_forest_df(forest,test)
    write.csv(pred, paste("Predictions/Hillstrom/forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)

    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 100,n_features = 3,
                                           criterion = c)
    pred <- predict_forest_df(forest,test)
    write.csv(pred, paste("Predictions/Hillstrom/random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control)
  write.csv(causal_forest_pred, paste("Predictions/Hillstrom/causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # Causal Tree
  # causal_pred <- causalTreePredicitons(train, test, treatment_list, response, control)
  # write.csv(causal_pred, paste("Predictions/Hillstrom/causal_tree",as.character(f),".csv",sep = ""),
  #           row.names = FALSE)

  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf")
  write.csv(pred_sma_rf, paste("Predictions/Hillstrom/sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)

  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 100, nrow(train), 3, 0.15, 100, parallel = TRUE,
                          remain_cores = 1)

  pred <- predict_forest_df(cts_forest, test)
  write.csv(pred, paste("Predictions/Hillstrom/cts",as.character(f),".csv",sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}


start_time <- Sys.time()
folder <- "Predictions/Hillstrom/"
outcomes <- c()
p_treated <- c()
n_treated_df <- c()
n_predictions <- 1
for(model in c("tree","forest","random_forest","cts","sma rf","causal_forest")){
  if(sum(model == c("tree","forest","random_forest")) > 0){
    for(c in c("simple","max","frac")){
      for(f in 1:n_predictions){
        pred <- read.csv(paste(folder,model,"_",c,as.character(f),".csv",sep = ""))
        if(length(outcomes) == 0){
          outcomes <- c(new_expected_quantile_response(response,control,treatment_list,pred),
                        paste(model,"_",c,sep = ""))
          n_treated_df <- c(n_treated_decile(pred,control),paste(model,"_",c,sep = ""))
          p_treated <- c(perc_treated(pred,control),paste(model,"_",c,sep = ""))
        } else{
          outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),
                                       paste(model,"_",c,sep = "")))
          n_treated_df <- rbind(n_treated_df,c(n_treated_decile(pred,control),paste(model,"_",c,sep = "")))
          p_treated <- rbind(p_treated,c(perc_treated(pred,control),paste(model,"_",c,sep = "")))
        }
      }
    }
  } else{
    for(f in 1:n_predictions){
      pred <- read.csv(paste(folder,model,as.character(f),".csv",sep = ""))
      outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),model))
      n_treated_df <- rbind(n_treated_df,c(n_treated_decile(pred,control),model))
      p_treated <- rbind(p_treated,c(perc_treated(pred,control),model))
    }
  }
}
outcome_df <- data.frame(outcomes)
n_treated_df <- data.frame(n_treated_df)
perc_treated_df <- data.frame(p_treated)
colnames(outcome_df) <- c(0,10,20,30,40,50,60,70,80,90,100,"Model")
colnames(n_treated_df) <- c(0,10,20,30,40,50,60,70,80,90,100,"Model")
colnames(perc_treated_df) <- c("PercTreated","Model")
rownames(outcome_df) <- 1:nrow(outcome_df)
rownames(n_treated_df) <- 1:nrow(n_treated_df)
rownames(perc_treated_df) <- 1:nrow(perc_treated_df)
for(c in 1:11){
  outcome_df[,c] <- as.numeric(as.character(outcome_df[,c]))
  n_treated_df[,c] <- as.numeric(as.character(n_treated_df[,c]))
}
outcome_df[,12] <- as.character(outcome_df[,12])
n_treated_df[,12] <- as.character(n_treated_df[,12])
perc_treated_df[,1] <- as.numeric(as.character(perc_treated_df[,1]))
perc_treated_df[,2] <- as.character(perc_treated_df[,2])
print(difftime(Sys.time(),start_time,units = "mins"))



for(model in unique(outcome_df$Model)){
  temp_data <- outcome_df[outcome_df$Model == model,]
  visualize(temp_data)
}




