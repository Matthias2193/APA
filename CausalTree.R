# install.packages("devtools")
# library(devtools) 
# install_github("susanathey/causalTree")
 library(causalTree)

causalTreePredicitons <- function(train, test,treatment_list, response){
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
    
    ## TODO - include this ??
    
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

newcausalTreePredicitons <- function(train, test,treatment_list, response,control){
  pred <- data.frame(rep(0,nrow(test)))
  for(t in treatment_list){
    train_data <- train[train[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,setdiff(colnames(train_data),c(control,treatment_list))]
    
    test_data <- test[,setdiff(colnames(train_data),response)]
    
    tree <- causalTree(as.formula(paste(response, "~.")), data = train_data_new, treatment = train_data[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                       cv.Honest = T, split.Bucket = F, xval = 5, cp = 0, minsize = 20, propensity = 0.5)
    
    ## TODO - include this ??
    opcp <- tree$cptable[,1][which.min(tree$cptable[,3])]
    pruned_tree <- rpart::prune(tree, opcp)
    
    assign(paste('predictions',t,sep = '_'),predict(tree,newdata = test_data))
    pred <- cbind(pred,eval(as.name(paste('predictions',t,sep = '_'))))
  }
  colnames(pred) <- c(control,treatment_list)
  for (t in treatment_list) {
    pred[,paste("uplift",t,sep = "_")] <- pred[t] - pred[control]
  }
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, treatment_list], 1, which.max)]
  return(pred)
}
