#The predictions functions

#Prediction---- 
#Tree  
predict.dt.as.df <- function(tree,new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(new_data)
  temp_function <- function(x,node){
    type = node[['type']]
    while(type != 'leaf'){
      split = node[['split']]
      if(type_list[[names(split)]] == 'factor'){
        if(x[names(split)] == split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      } else{
        if(as.numeric(x[names(split)]) < split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      }
    }
    return(node[["results"]])
  }
  results <- data.frame(t(apply(new_data,1,temp_function,node=tree)))
  return(results)
}


#Takes predictions as input and returns just the name of the best treatment for each prediction

predictions_to_treatment <- function(pred,treatment_list,control){
  colnames(pred[,c(treatment_list,control)])[apply(pred[,c(treatment_list,control)],1,which.max)]
}

new_predictions_to_treatment <- function(pred,treatment_list,control){
  results <- colnames(pred[,c(treatment_list,control)])[apply(pred[,c(treatment_list,control)],1,which.max)]
  for(c in treatment_list){
    if(exists(temp_bool)){
      temp_bool <- temp_bool & pred[,c] == 0
    } else{
      temp_bool <- pred[,c] == 0
    }
  }
  if(sum(temp_bool) > 0){
    results[temp_bool] <- control
  }
}

FALSE + TRUE

#Forest

#Majority Vote
predict_forest_majority <- function(forest,test_data){
  predictions <- predictions_to_treatment(predict.dt(forest[[1]],test_data))
  for(x in 2:length(forest)){
    predictions <- cbind(predictions,predictions_to_treatment(predict.dt(forest[[x]],test_data)))
  }
  return(apply(data.frame(predictions), MARGIN = 1, FUN = forest_predictions_helper))
}

forest_predictions_helper <- function(preds){
  temp_names <- unique(unlist(preds))
  if(length(temp_names) == 1){
    return(temp_names)
  } else{
    return(names(which.max(table(unlist(preds)))))
  }
}


#Average
#Sequentially
predict_forest_average <- function(forest,test_data){
  predictions <- list()
  for(x in 1:length(forest)){
    predictions[[x]] <- predict.dt.as.df(forest[[x]],test_data)
  }
  final_predictions <- predictions[[1]]
  for(x in 2:length(predictions)){
    final_predictions <- final_predictions+predictions[[x]]
  }
  final_predictions <- final_predictions/length(predictions)
  return(final_predictions)
}

#Parallel
parallel_predict_forest_average <- function(forest,test_data,remain_cores = 1){
  predictions <- list()
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores-remain_cores)
  registerDoParallel(cl)
  predictions <- foreach(x = 1:length(forest)) %dopar%{
    source('DecisionTreeImplementation.R')
    tree <- forest[[x]]
    new_data <- test_data
    type_list <- sapply(new_data, class)
    names(type_list) = colnames(new_data)
    temp_function <- function(x,node){
      type = node[["type"]]
      while(type != 'leaf'){
        split = node[['split']]
        if(type_list[[names(split)]] == 'factor'){
          if(x[names(split)] == split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        } else{
          if(as.numeric(x[names(split)]) < split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        }
      }
      return(node[["results"]])
    }
    results <- data.frame(t(apply(new_data,1,temp_function,node=tree)))
    return(results)
  }
  stopCluster(cl)
  final_predictions <- predictions[[1]]
  for(x in 2:length(predictions)){
    final_predictions <- final_predictions+predictions[[x]]
  }
  final_predictions <- final_predictions/length(predictions)
  return(final_predictions)
}





predict_forest_df <- function(forest,test_data, parallel_pred = TRUE, remain_cores = 1){
  if(parallel_pred){
    return(parallel_predict_forest_average(forest,test_data,remain_cores))
  } else{
    return(predict_forest_average(forest,test_data))
  }
}