library(causalTree)

causalForestPredicitons <- function(train, test,treatment_list, response){
  for(t in treatment_list){
    train_data <- train
    train_data <- train_data[train_data[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,1:8]
    train_data_new[, response] <- train_data[, response]
    #train_data <- train_data_new
    
    test_data <- test[,1:8]
    test_data[,t] <- test[,t]
    
    if(file.exists(paste(paste('models/forest',response,t,sep = '_'),'rda',sep='.'))){
      load(paste(paste('models/forest',response,t,sep = '_'),'rda',sep='.'))
    }
    else{
      forest <- causalForest(as.formula(paste(response, "~ recency+history_segment+history+mens+womens+zip_code+newbie+channel")), data = train_data_new, treatment = train_data[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                         cv.Honest = T, split.Bucket = F, minsize = 20, propensity = 0.5, mtry = 2, num.trees = 200, ncov_sample = 3, ncolx = ncol(train_data_new)-1)
      save(forest, file = paste(paste('models/forest',response,t,sep = '_'),'rda',sep='.'))
    }
    
    assign(paste('predictions',t,sep = '_'), predict(forest, newdata = test_data))
  }
  
  pred <- data.frame(cbind(predictions_men_treatment,predictions_women_treatment))
  colnames(pred) <- treatment_list
  pred$control <- 0
  pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
  pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
  
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  return(pred)
}

