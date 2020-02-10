# Implementation of the Seperate Model Approach for Multiple Treatments and different Base Learners
list.of.packages <- c("caret", "glmnet", "rpart")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(caret)
library(glmnet)
library(rpart)
library(randomForest)

########################################
# Decision Tree and Random Forest Models
########################################

dt_models <- function(train_data, response, prediction_method,treatments,control, test_data, model_type){  
  assignment <- predictions_to_treatment(test_data, treatment_list, control)
  outcome <- test_data[, response]
  for(t in c(treatments,control)){
    test_data[,t] <- NULL
  }
  test_data[, response] <- NULL
  
  for(t in c(treatments,control)){
    temp_train <- train_data[train_data[,t]==1,]
    
    for(s in c(treatments,control)){
      temp_train[,s] <- NULL
    }
    if(model_type == "dt"){
      dt <- rpart(as.formula(paste(response, "~.")), temp_train, method = prediction_method, 
                  control=rpart.control(minsplit=5,  cp=0.0001))
    } else if(model_type == "rf"){
      if(prediction_method == "class"){
        temp_train[,response] <- as.factor(temp_train[,response])
      }
      dt <- randomForest(as.formula(paste(response, "~.")), data = temp_train, mtry=3, ntree = 300)
    } else{
      print("No valid model type selected. Please choose 'dt' for a single tree or 'rf' for a random forest.")
      return()
    }
    
    if(prediction_method == "class"){
      if(!exists("predictions")){
        predictions <-data.frame(predict(dt, test_data, type = "prob")[ , 2])
        colnames(predictions) <- c(t)
      } else{
        predictions[ ,t] <- as.numeric(predict(dt, test_data, type = "prob")[ , 2] )
      }
      # Regression prediction
    } else if(prediction_method == "anova"){
      if(t == treatments[1]){
        predictions <- data.frame(predict(dt, test_data))
        colnames(predictions) <- c(t)
      } else{
        predictions[ , t] <- as.numeric(predict(dt, test_data))
      }
    }
    
  }
  
  colnames(predictions) <- c(treatment_list,control)
  predictions[ , "Treatment"] <- predictions_to_treatment(predictions, treatment_list, control)
  predictions[ , "Assignment"] <- assignment
  predictions[, "Outcome"] <- outcome
  for (t in treatment_list) {
    predictions[,paste("uplift",t,sep = "_")] <- predictions[t] - predictions[control]
  }
  
  
  return(predictions)
}
