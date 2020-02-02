#Test Script


library(ggplot2)
library(caret)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('RzepakowskiTree.R')
source('Evaluation Methods.R')


set.seed(2193)

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

#email$spend <- email$spend/max(email$spend)

idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

test <- email[idx, ]

# # Partition training data for pruning
# p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
# 
# val <- train[p_idx,]
# train <- train[-p_idx,]
# 
treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE, max_cases = 10)


val <- train





forest <- parallel_build_random_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3, pruning = F)
pred <- predict_forest_df(forest,test)
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
exp_forest_random <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_forest_random <- new_expected_quantile_response(response,control,treatment_list,pred)


forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3, pruning = F)
pred <- predict_forest_df(forest,test)
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
exp_forest_all <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_forest_all <- new_expected_quantile_response(response,control,treatment_list,pred)

pred2 <- read.csv('Predictions/forest_Old.csv')
exp_forest_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_forest_old <- new_expected_quantile_response(response,control,treatment_list,pred)

#Forest
build_random_forest <- function(train_data, val_data,treatment_list,response,control,n_trees,n_features,
                         pruning,max_depth = 10){
  trees <- list()
  for(x in 1:n_trees){
    set.seed(rnorm(1))
    temp_train_data <- train_data[sample(nrow(train_data), nrow(train_data),replace = TRUE),]
    temp_tree <- build_tree(data = temp_train_data,0,treatment_list = treatment_list, 
                            test_list = test_list,target = response,control = control,
                            max_depth = max_depth,random = TRUE,n_features)
    if(pruning){
      temp_prune_tree <- simple_prune_tree(temp_tree,val_data[,chosen_cols], treatment_list, test_list, response, control)
      trees[[x]] <- temp_prune_tree
    } else{
      trees[[x]] <- temp_tree
    }
    
  }
  return(trees)
}



parallel_build_random_forest <- function(train_data, val_data,treatment_list,response,control,n_trees,n_features,
                                  pruning,max_depth = 10,remain_cores = 1){
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores-remain_cores)
  registerDoParallel(cl)
  trees <- foreach(x=1:n_trees) %dopar% {
    source('DecisionTreeImplementation.R')
    set.seed(rnorm(1))
    temp_train_data <- train_data[sample(nrow(train_data), nrow(train_data),replace = TRUE),]
    temp_tree <- build_tree(data = temp_train_data,0,treatment_list = treatment_list, 
                            test_list = test_list,target = response,control = control,
                            max_depth = max_depth,random = TRUE,n_features)
    return(temp_tree)
    if(pruning){
      temp_prune_tree <- simple_prune_tree(temp_tree,val_data[,chosen_cols], treatment_list,
                                           test_list, response,control = control)
      return(temp_prune_tree)
    } else{
      return(temp_tree)
    }
  }
  stopCluster(cl)
  return(trees)
}
