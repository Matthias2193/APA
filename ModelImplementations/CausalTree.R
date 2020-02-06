#In this script we adapt the causal tree and causal forest (https://github.com/susanathey/causalTree) 
#for multiple treatments.For each treatment we build one causal tree/forest to estimate the uplift for that treatment
#compared to control. Then we compare the eastimated uplift for all the treatments.

# install.packages("devtools")
# library(devtools) 
# install_github("susanathey/causalTree")
library(causalTree)


causalTreePredicitons <- function(train, test,treatment_list, response,control){
  pred <- data.frame(rep(0,nrow(test)))
  for(t in treatment_list){
    train_data <- train[train[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,setdiff(colnames(train_data),c(control,treatment_list))]
    
    test_data <- test[,setdiff(colnames(train_data),response)]
    
    tree <- causalTree(as.formula(paste(response, "~.")), data = train_data_new, treatment = train_data[,t], 
                       split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                       val = 5, cp = 0, minsize = 20, propensity = 0.5)

    opcp <- tree$cptable[,1][which.min(tree$cptable[,3])]
    pruned_tree <- rpart::prune(tree, opcp)
    
    assign(paste('predictions',t,sep = '_'),predict(tree,newdata = test_data))
    pred <- cbind(pred,eval(as.name(paste('predictions',t,sep = '_'))))
  }
  colnames(pred) <- c(control,treatment_list)
  for (t in treatment_list) {
    pred[,paste("uplift",t,sep = "_")] <- pred[t] - pred[control]
  }
  pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)
  pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
  pred[, "Outcome"] <- test[,response]
  return(pred)
}



causalForestPredicitons <- function(train,test,treatment_list,response,control){
  pred <- data.frame(rep(0,nrow(test)))
  for(t in treatment_list){
    train_data <- train[train[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,setdiff(colnames(train_data),c(control,treatment_list))]
    
    test_data <- test[,setdiff(colnames(train_data),response)]
    
    forest <- causalForest(as.formula(paste(response,paste(setdiff(colnames(train_data_new),response),collapse = "+"),
                                            sep = "~")), data = train_data_new, treatment = train_data[,t], 
                           split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                           minsize = 20, propensity = 0.5, mtry = 2, num.trees = 200, ncov_sample = 3, 
                           ncolx = ncol(train_data_new)-1)
    
    assign(paste('predictions',t,sep = '_'), predict(forest, newdata = test_data))
    pred <- cbind(pred,eval(as.name(paste('predictions',t,sep = '_'))))
  }
  colnames(pred) <- c(control,treatment_list)
  for (t in treatment_list) {
    pred[,paste("uplift",t,sep = "_")] <- pred[t] - pred[control]
  }
  pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)
  pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
  pred[, "Outcome"] <- test[,response]
  return(pred)
}






#Old Versions
oldcausalTreePredicitons <- function(train, test,treatment_list, response){
  for(t in treatment_list){
    train_data <- train
    train_data <- train_data[train_data[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,1:8]
    train_data_new[, response] <- train_data[, response]
    train_data_new[,t] <- train_data[,t]
    train_data <- train_data_new
    
    test_data <- test[,1:8]
    test_data[,t] <- test[,t]
    
    tree <- causalTree(as.formula(paste(response, "~.")), data = train_data, treatment = train_data[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                       cv.Honest = T, split.Bucket = F, xval = 5, cp = 0, minsize = 20, propensity = 0.5)
    
    opcp <- tree$cptable[,1][which.min(tree$cptable[,3])]
    pruned_tree <- rpart::prune(tree, opcp)
    
    assign(paste('predictions',t,sep = '_'),predict(pruned_tree,newdata = test_data))
  }
  
  pred <- data.frame(cbind(predictions_men_treatment,predictions_women_treatment))
  colnames(pred) <- treatment_list
  pred$control <- 0
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  return(pred)
}



oldcausalForestPredicitons <- function(train,test,treatment_list, response){
  for(t in treatment_list){
    train_data <- train
    train_data <- train_data[train_data[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,1:8]
    train_data_new[, response] <- train_data[, response]
    
    test_data <- test[,1:8]
    test_data[,t] <- test[,t]
    
    forest <- causalForest(as.formula(paste(response, "~ recency+history_segment+history+mens+womens+zip_code+newbie+channel")), data = train_data_new, treatment = train_data[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                           cv.Honest = T, split.Bucket = F, minsize = 20, propensity = 0.5, mtry = 2, num.trees = 200, ncov_sample = 3, ncolx = ncol(train_data_new)-1)
    
    assign(paste('predictions',t,sep = '_'), predict(forest, newdata = test_data))
  }
  
  pred <- data.frame(cbind(predictions_men_treatment,predictions_women_treatment))
  colnames(pred) <- treatment_list
  pred$control <- 0
  return(pred)
}