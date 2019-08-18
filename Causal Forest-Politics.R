library(causalTree)

causalForestPredicitons_Politics <- function(train, test,treatment_list, response){
  for(t in treatment_list){
    train_data <- train
    train_data <- train_data[train_data[,setdiff(treatment_list,t)] == 0,]
    
    train_data_new <- train_data[,1:4]
    train_data_new[, response] <- train_data[, response]
    
    test_data <- test[,1:4]
    test_data[,t] <- test[,t]
    
    forest <- causalForest(as.formula(paste(response, "~ female+age+hh_size")), data = train_data_new, treatment = train_data[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                           cv.Honest = T, split.Bucket = F, minsize = 20, propensity = 0.5, mtry = 2, num.trees = 20, ncov_sample = 3, ncolx = ncol(train_data_new)-1)
    
    assign(paste('predictions',t,sep = '_'), predict(forest, newdata = test_data))
  }
  
  pred <- data.frame(cbind(predictions_treatment_CivicDuty, predictions_treatment_Self, predictions_treatment_Hawthorne,predictions_treatment_Neighbors))
  colnames(pred) <- treatment_list
  pred$control <- 0
  
  return(pred)
}

