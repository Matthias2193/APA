#The predictions functions

#Prediction---- 
#Tree  
predict.dt.as.df <- function(tree, new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(new_data)
  results = data.frame()
  for(x in 1:nrow(new_data)){
    d = new_data[x,]
    type = 'root'
    node = tree
    while(type != 'leaf'){
      split = node[['split']]
      if(type_list[[names(split)]] == 'factor'){
        if(d[names(split)] == split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      } else{
        if(d[names(split)] < split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      }
    }
    results <- rbind(results , node[['results']])
    
    # in first round colnames need to be set T names
    if(x == 1){
      colnames(results) <- names(node[['results']])
    }
    
  }
  
  return(results)
}

predict.dt <- function(tree,new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(new_data)
  results = list()
  for(x in 1:nrow(new_data)){
    d = new_data[x,]
    type = 'root'
    node = tree
    while(type != 'leaf'){
      split = node[['split']]
      if(type_list[[names(split)]] == 'factor'){
        if(d[names(split)] == split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      } else{
        if(d[names(split)] < split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      }
    }
    results[[x]] <- node[['results']]
  }
  return(results)
}

predictions_to_df <- function(predictions){
  results <- data.frame()
  for(x in 1:length(predictions)){
    results <- rbind(results, predictions[[x]])
  }
  colnames(results) <- names(predictions[[1]])
  return(results)
}

#Takes predictions as input and returns just the name of the best treatment for each prediction

predictions_to_treatment <- function(pred){
  lapply(pred, predictions_to_treatment_helper)
}

predictions_to_treatment_helper <- function(x){
  names(x[match(max(x),x)])
}


#Forest
predict_forest_majority <- function(forest,test_data){
  predictions <- predictions_to_treatment(predict.dt(forest[[1]],test_data))
  for(x in 2:length(forest)){
    predictions <- cbind(predictions,predictions_to_treatment(predict.dt(forest[[x]],test_data)))
  }
  return(apply(data.frame(predictions), MARGIN = 1, FUN = forest_predictions_helper))
}

predict_forest_average <- function(forest,test_data){
  predictions <- list()
  for(x in 1:length(forest)){
    predictions[[x]] <- predict.dt(forest[[x]],test_data)
  }
  final_predictions <- list()
  for(x in 1:nrow(test_data)){
    temp <- unlist(predictions[[1]][x])
    for(y in 2:length(predictions)){
      temp <- temp + unlist(predictions[[y]][x])
    }
    temp <- temp/length(predictions)
    final_predictions[[x]] <- temp
  }
  return(final_predictions)
}
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
    results = list()
    for(y in 1:nrow(new_data)){
      d = new_data[y,]
      type = tree[['type']]
      node = tree
      while(type != 'leaf'){
        split = node[['split']]
        if(type_list[[names(split)]] == 'factor'){
          if(d[names(split)] == split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        } else{
          if(d[names(split)] < split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        }
      }
      results[[y]] <- node[['results']]
    }
    return(results)
  }
  stopCluster(cl)
  final_predictions <- list()
  for(x in 1:nrow(test_data)){
    temp <- unlist(predictions[[1]][x])
    for(y in 2:length(predictions)){
      temp <- temp + unlist(predictions[[y]][x])
    }
    temp <- temp/length(predictions)
    final_predictions[[x]] <- temp
  }
  return(final_predictions)
}


forest_predictions_helper <- function(preds){
  temp_names <- unique(unlist(preds))
  if(length(temp_names) == 1){
    return(temp_names)
  } else{
    return(names(which.max(table(unlist(preds)))))
  }
}

predict_forest_df <- function(forest,test_data){
  temp_pred <- predict_forest_average(forest,test_data)
  return(predictions_to_df(temp_pred))
}
parallel_predict_forest_df <- function(forest,test_data){
  temp_pred <- parallel_predict_forest_average(forest,test_data)
  return(predictions_to_df(temp_pred))
}
