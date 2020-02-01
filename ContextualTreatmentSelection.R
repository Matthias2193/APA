library(dplyr)
library(foreach)
library(doParallel)
source("DecisionTreeImplementation.R")
build_cts <- function(response, control, treatments, data, ntree, B, m_try, n_reg, min_split, parallel = TRUE,
                      remain_cores = 1){
  if(parallel){
    numCores <- detectCores()
    cl <- makePSOCKcluster(numCores-remain_cores)
    registerDoParallel(cl)
    trees <- foreach(x=1:ntree) %dopar% {
      source('ContextualTreatmentSelection.R')
      set.seed(x)
      for(t in c(treatments,control)){
        if(t == treatments[1]){
          temp_train_data <- sample_n(data[data[,t]==1,], B * (sum(data[,t]==1)/nrow(data)))
        } else{
          temp_train_data <- rbind(temp_train_data,sample_n(data[data[,t]==1,], B * (sum(data[,t]==1)/nrow(data))))
        }
        
      }
      return(build_cts_tree(response, control, treatments, temp_train_data, m_try, n_reg, 
                                  min_split, parent_predictions = NA))
      print(paste(as.character(x),"completed!",sep = " "))
    }
    stopCluster(cl)
    return(trees)
  }
  else{
    trees <- list()
    for(x in 1:ntree){
      for (t in c(treatments,control)){
        if(t == treatments[1]){
          temp_train_data <- sample_n(data[data[,t]==1,], B * (sum(data[,t]==1)/nrow(data)))
        } else{
          temp_train_data <- rbind(temp_train_data,sample_n(data[data[,t]==1,], B * (sum(data[,t]==1)/nrow(data))))
        }
        
      }
      trees[[x]] <- build_cts_tree(response, control, treatments, temp_train_data, m_try, n_reg, 
                                   min_split, parent_predictions = NA)
      
    }
    return(trees)
  }
}

build_cts_tree <- function(response, control, treatments, data, m_try, n_reg, min_split, 
                           min_gain = 0, parent_predictions = NA,depth = 0){
  node <- list()
  if(nrow(data) == 0){
    node[["type"]] <- "leaf"
    node[["results"]] <- parent_predictions
    return(node)
  }
  retain_cols <- c(treatments,control,response)
  sample_cols <- setdiff(colnames(data),retain_cols)
  temp_cols <- sample(sample_cols,m_try,replace = F)
  chosen_cols <- c(temp_cols,retain_cols)
  test_list<- set_up_tests(data[,temp_cols],TRUE)
  
  node <- list()
  results <- c()
  #Select split with maximum gain
  n_treatment_samples <- c()
  if(is.na(parent_predictions)){
    node[["type"]] <- "root"
    for(t in c(treatments,control)){
      results <- c(results, mean(data[data[,t] == 1,response]))
      n_treatment_samples <- c(n_treatment_samples,nrow(data[data[,t] == 1,]))
    }
    names(results) <- c(treatments,control)
    node[["results"]] <- results
  } else{
    node[["type"]] <- "node"
    for(t in c(treatments,control)){
      n_samples <- nrow(data[data[,t] == 1,])
      n_treatment_samples <- c(n_treatment_samples,n_samples)
      if(n_samples < min_split){
        results <- c(results, parent_predictions[[t]])
      } else{
        results <- c(results, 
                     (sum(data[data[,t] == 1,response])+parent_predictions[[t]]*n_reg)/(sum(data[,t] == 1)+n_reg))
      }
    }
  }
  names(results) <- c(treatments,control)
  names(n_treatment_samples) <- c(treatments,control)
  node[["results"]] <- results
  node[['n_samples']] <- nrow(data)
  node[['n_samplse_treatments']] <- n_treatment_samples
  terminate <- TRUE
  for(t in c(treatments,control)){
    if(sum(data[,t]==1) >= min_split){
      terminate = FALSE
    }
  }
  
  if(sum(data[,response] == data[1,response]) == nrow(data)){
    terminate = TRUE
  }
  
  if(terminate){
    node[["type"]] <- "leaf"
    return(node)
  }
  
  temp_split <- select_cts_split(node[["results"]], data, treatments, response, control, n_reg, min_split, 
                                   test_list, min_gain)
  
  if(temp_split == -1){
    node[["type"]] <- "leaf"
    return(node)
  }
  #Return a leaf, if there is no split with gain > 0
  
  node[['split']] <- temp_split
  if(names(temp_split) %in% names(test_list$categorical)){
    node[['left']] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)]==temp_split[[1]],], 
                                     m_try, n_reg, min_split, min_gain, node[["results"]],depth = depth + 1)
    node[['right']] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)]!=temp_split[[1]],], 
                                      m_try, n_reg, min_split, min_gain, node[["results"]],depth = depth + 1)
  } else{
    node[['left']] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)]<temp_split[[1]],], 
                                     m_try, n_reg, min_split, min_gain, node[["results"]],depth = depth + 1)
    node[['right']] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)]>=temp_split[[1]],], 
                                      m_try, n_reg, min_split, min_gain, node[["results"]],depth = depth + 1)
  }
  return(node)
}

select_cts_split <- function(parent_predictions, data, treatments, response, control, n_reg, min_split, 
                             test_list, min_gain){
  gain_list <- c()
  name_list <- c()
  if(length(test_list$categorical)>0){
    for(x in 1:length(test_list$categorical)){
      col_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        test_case <- test_list$categorical[[x]][y]
        new_name <- paste(col_name, as.character(test_case), sep = '@@')
        gain_list <- c(gain_list,cts_gain(test_case,treatments,control,response,data,'categorical',col_name,
                                          parent_predictions, n_reg, min_split))
        name_list <- c(name_list,new_name)
      }
    }
  }
  if(length(test_list$numerical)>0){
    for(x in 1:length(test_list$numerical)){
      col_name <- names(test_list$numerical[x])
      for(y in 1:length(test_list$numerical[[x]])){
        test_case <- test_list$numerical[[x]][y]
        new_name <- paste(col_name, as.character(test_case), sep = '@@')
        gain_list <- c(gain_list,cts_gain(test_case,treatments,control,response,data,'numerical',col_name,
                                          parent_predictions, n_reg, min_split))
        name_list <- c(name_list,new_name)
      }
    }
  }
  if(is.na(max(gain_list))){
    return(-1)
  }
  if(max(gain_list) > min_gain){
    temp_string <- name_list[match(max(gain_list),gain_list)]
    temp_result <- strsplit(temp_string,split='@@', fixed=TRUE)
    options(warn=-1)
    if(is.na(as.numeric(temp_result[[1]][2]))){
      result <- c(temp_result[[1]][2])
      names(result) <- temp_result[[1]][1]
      options(warn=0)
      return(result)
    } else{
      options(warn=0)
      result <- c(as.numeric(temp_result[[1]][2]))
      names(result) <- temp_result[[1]][1]
      return(result)
    }
  } else{
    return(-1)
  }
}

cts_gain <- function(test_case, treatments, control, response, data, test_type, col_name, parent_predictions, 
                     n_reg, min_split){
  max_left <- 0
  max_right <- 0
  max_root <- max(parent_predictions)
  if(test_type == "categorical"){
    data_left <- data[data[,col_name]!= test_case,]
    data_right <- data[data[,col_name] == test_case, ]
  } else{
    data_left <- data[data[,col_name] < test_case,]
    data_right <- data[data[,col_name] >= test_case, ]
  }
  p_left <- nrow(data_left)/nrow(data)
  p_right <- nrow(data_right)/nrow(data)
  for(t in c(treatments,control)){
    n_left <- nrow(data_left[data_left[,t] == 1,])
    n_right <- nrow(data_right[data_right[,t] == 1,])
    if(n_left < min_split){
      exp_left <- parent_predictions[[t]]
    } else{
      exp_left <- (sum(data_left[data_left[,t] == 1,response])+parent_predictions[[t]]*n_reg)/(sum(data_left[,t] ==1)+n_reg)
    }
    if(n_right < min_split){
      exp_right <- parent_predictions[[t]]
    } else{
      exp_right <- (sum(data_right[data_right[,t] == 1,response])+parent_predictions[[t]]*n_reg)/(sum(data_right[,t]==1)+n_reg)
    }
    max_left <- max(max_left, exp_left)
    max_right <- max(max_right, exp_right)
  }
  return(p_left * max_left + p_right * max_right - max_root)
}

check_split <- function(tree){
  if(tree[["type"]] == "leaf"){
    #print("Reached End")
  } else if((tree[["left"]][["n_samples"]] + tree[["right"]][["n_samples"]]) == tree[["n_samples"]]){
    if(tree[["left"]][["n_samples"]] == 0){
      print("Left no samples")
    } else if(tree[["right"]][["n_samples"]] == 0){
      print("Right no samples")
    }
    check_split(tree[["left"]])
    check_split(tree[["right"]])
  } else{
    print("Subsamples don't sum up!")
    break
  }
}
