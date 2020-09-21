#This script contains the implementation of a new decision tree for uplift modeling.

#Importing libraries
library(parallel)
library(foreach)
library(doParallel)


#Set up tests ----
#Creates a list of possible splits.
#Parameters:
#reduce_cases(Boolean): Allows to limit the number of possible splits for continuous covariates
#max_cases(Integer): The maximum number of possible splits to be examined for each covariate. These splits will
#be evenly distributed over the range of the covariate
set_up_tests <- function(x,reduce_cases,max_cases = 10){
  type_list <- sapply(x, class)
  categorical_splits = list()
  numerical_splits = list()
  for(n in colnames(x)){
    if(type_list[[n]] == 'factor'){
      categorical_splits[[n]]<-levels(x[,n])
    }
    else{
      temp_list <- sort(unique(x[,n]))
      r <- 1
      s <- 2
      final_list <- c()
      while(r < length(temp_list)){
        final_list <- c(final_list,round((temp_list[r]+temp_list[s])/2,1))
        r <- r+1
        s <- s+1
      }
      final_list <- unique(final_list)
      if((length(final_list)>max_cases)&& reduce_cases){
        final_list <- round(final_list)
        final_list <- unique(final_list)
      }
      if((length(final_list)>max_cases)&& reduce_cases){
        final_list <- quantile(x[,n],seq(0,1,1/max_cases))
        final_list <- round(final_list,2)
      }
      numerical_splits[[n]] <- final_list
    }
  }
  output = list()
  output[['categorical']] <- categorical_splits
  output[['numerical']] <- numerical_splits
  return(output)
}





#Split selection ----
#For each possible split the gain is calculated. The split with the highest gain is returned. If no split has a
#gain > 0, then -1 is returned indicating that no split is beneficial.
#Curretnly there are three criterions supported to calculate the gain. More explanation further down.
select_split <- function(test_list,treatment,control,target,temp_data,criterion,
                         min_split, parent_predictions){
  gain_list <- c()
  name_list <- c()
  if(length(test_list$categorical)>0){
    for(x in 1:length(test_list$categorical)){
      temp_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        t <- test_list$categorical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        if(criterion == "simple"){
          gain_list <- c(gain_list,simple_gain(t,treatment,control,target,temp_data,'categorical',temp_name,
                                               min_split, parent_predictions))
        } else if( criterion == "max"){
          gain_list <- c(gain_list,max_gain(t,treatment,control,target,temp_data,'categorical',temp_name,
                                            min_split, parent_predictions))
        } else{
          gain_list <- c(gain_list,frac_gain(t,treatment,control,target,temp_data,'categorical',temp_name,
                                             min_split, parent_predictions))
        }
        
        name_list <- c(name_list,new_name)
      }
    }
  }
  if(length(test_list$numerical)>0){
    for(x in 1:length(test_list$numerical)){
      temp_name <- names(test_list$numerical[x])
      for(y in 1:length(test_list$numerical[[x]])){
        t <- test_list$numerical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        if(criterion == "simple"){
          gain_list <- c(gain_list,simple_gain(t,treatment,control,target,temp_data,'numerical',temp_name,
                                               min_split, parent_predictions))
        } else if( criterion == "max"){
          gain_list <- c(gain_list,max_gain(t,treatment,control,target,temp_data,'numerical',temp_name,
                                            min_split, parent_predictions))
        } else{
          gain_list <- c(gain_list,frac_gain(t,treatment,control,target,temp_data,'numerical',temp_name,
                                             min_split, parent_predictions))
        }
        name_list <- c(name_list,new_name)
      }
    }
  }
  if(is.na(max(gain_list))){
    return(-1)
  }
  if(max(gain_list) > 0){
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
  }
  else{
    return(-1)
  }
}



#Simple
#The first of 3 possible ways to calculate the "gain". This option tries to maximize the difference in expected
#outcome between the treatments and control and between treatments.
#Doesn't take into accout current results or number of samples in each node after the split.
simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col,min_split,parent_predictions){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
  } else{
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
  }
  if((nrow(data) == 0) || nrow(data1) == 0 || nrow(data2) == 0 ){
    return(-1)
  }
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  #Here the gain is calculated
  for(t in treatments){
    for(s in treatments){
      if(sum(data1[,t]==1)<min_split){
        mean_t1 <- parent_predictions[[t]]
      }else{
        mean_t1 <- mean(data1[data1[,t] == 1,target])
      }
      if(sum(data[,s]==1)<min_split){
        mean_s1 <- parent_predictions[[s]]
      }else{
        mean_s1 <- mean(data1[data1[,s] == 1,target])
      }
      if(sum(data2[,t]==1)<min_split){
        mean_t2 <- parent_predictions[[t]]
      }else{
        mean_t2 <- mean(data2[data2[,t] == 1,target])
      }
      if(sum(data2[,s]==1)<min_split){
        mean_s2 <- parent_predictions[[s]]
      }else{
        mean_s2 <- mean(data2[data2[,s] == 1,target])
      }
      temp_gain <- (mean_t1 - mean_s1)^2 + (mean_t2 - mean_s2)^2
      gain <- gain + temp_gain
    }
  }
  if(is.na(gain)){
    gain = -1
  }
  return(gain)
}


#Frac
#Similar to "Simple" but makes sure the difference in outcome after the split is bigger than before. Also takes into
#account the number of samples in each node after the split.
frac_gain <- function(test_case, treatment, control, target, data, test_type, test_col,min_split, parent_predictions){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
  } else{
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
  }
  if((nrow(data) == 0) || nrow(data1) == 0 || nrow(data2) == 0 ){
    return(-1)
  }
  frac1 <- nrow(data1)/nrow(data)
  frac2 <- nrow(data2)/nrow(data)
  
  current_gain <- 0
  for(x in 1:(length(treatments)-1)){
    t <- treatments[x]
    for(y in (x+1):length(treatments)){
      s <- treatments[y]
      temp_gain <- (parent_predictions[[t]]-parent_predictions[[s]])^2
      current_gain <- current_gain + temp_gain
    }
  }
  
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  left_results <- c()
  right_results <- c()
  for(t in treatments){
    if(sum(data1[,t]==1)<min_split){
      left_results <- c(left_results,parent_predictions[[t]])
    }else{
      left_results <- c(left_results,mean(data1[data1[,t] == 1,target]))
    }
    if(sum(data2[,t]==1)<min_split){
      right_results <- c(parent_predictions[[t]],right_results)
    }else{
      right_results <- c(right_results,mean(data2[data2[,t] == 1,target]))
    }
  }
  for(x in 1:(length(treatments)-1)){
    for(y in (x+1):length(treatments)){
      temp_gain <- frac1*(left_results[[x]] - left_results[[y]])^2 + frac2*(right_results[[x]] - right_results[[y]])^2
      gain <- gain + temp_gain
    }
  }
  if(is.na(gain)){
    gain = -1
  }
  if(gain <= current_gain){
    gain = -1
  }
  return(gain)
}

#Max
#Simply tries to maximize the maximum expected outcome for all treatments and control.
max_gain <- function(test_case, treatment, control, target, data, test_type, test_col,min_split, parent_predictions){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
  } else{
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
  }
  if((nrow(data) == 0) || nrow(data1) == 0 || nrow(data2) == 0 ){
    return(-1)
  }
  frac1 <- nrow(data1)/nrow(data)
  frac2 <- nrow(data2)/nrow(data)
  
  current_results <- parent_predictions
  left_results <- c()
  right_results <- c()
  for(t in treatments){
    if(sum(data1[,t]==1)<min_split){
      left_results <- c(left_results, parent_predictions[[t]])
    }else{
      left_results <- c(left_results, mean(data1[data1[,t] == 1,target]))
    }
    if(sum(data2[,t]==1)<min_split){
      right_results <- c(right_results, parent_predictions[[t]])
    }else{
      right_results <- c(right_results, mean(data2[data2[,t] == 1,target]))
    }
    
  }
  current_gain <- (max(current_results) - sort(current_results,decreasing = T)[2])^2
  right_gain <- (max(right_results) - sort(right_results,decreasing = T)[2])^2
  left_gain <- (max(left_results) - sort(left_results,decreasing = T)[2])^2
  gain <- frac1*left_gain + frac2*right_gain - current_gain
  if(is.na(gain)){
    return(-1)
  } else if(gain < 0){
    return(-1)
  } else{
    return(gain)
  }
}




#Functions to build the tree ----

#Tree
#data: the training data
#depth: the current depth, can be ignored by the user, for internal pruposes
#treatment_list: a list with the names of all treatments
#target: the name of the response variable
#control: the name of the control 'treatment'
#test_list: a list of possible splits created by the 'set_up_tests' function
#random and n_features are used for building the random forest.
#criterion sepcifies which criterion to use to calculate the gain ("simple", "max", "frac")
build_tree <- function(data,depth,max_depth,treatment_list,target,control,test_list,random = FALSE,
                       n_features = 0,criterion = "simple",min_split = 0,parent_predictions = NULL){
  #Return leaf if current depth is max depth
  if(depth == max_depth){
    return(final_node(data,treatment_list,target,control))
  }
  #Create current node
  if(random){
    retain_cols <- c(treatment_list,control,target)
    sample_cols <- setdiff(colnames(data),retain_cols)
    temp_cols <- sample(sample_cols,n_features,replace = F)
    chosen_cols <- c(temp_cols,retain_cols)
    test_list<- set_up_tests(data[,temp_cols],TRUE)
  }
  treatment_names <- c(treatment_list,control)
  effects <- c()
  if(is.null(parent_predictions)){
    for(t in treatment_names){
      effects <- c(effects,mean(data[data[t]==1,target]))
    }
  } else{
    for(t in treatment_names){
      if(sum(data[,t])<min_split){
        effects <- c(effects,parent_predictions[[t]])
      }else{
        effects <- c(effects,mean(data[data[t]==1,target]))
      }
    }
  }
  names(effects) <- treatment_names
  node <- list()
  #Select split with maximum gain
  temp_split <- select_split(test_list = test_list, treatment = treatment_list, control, target,data,criterion,
                             min_split = min_split, parent_predictions = effects)
  #Return a leaf, if there is no split with gain > 0
  if(temp_split == -1){
    return(final_node(data,treatment_list,target,control))
  }
  #Construct current node
  node[['type']] = 'node'
  if(depth == 0){
    node[['type']] <- 'root'
  }
  #Number of training samples in current node
  node[['n_samples']] <- nrow(data)
  #The estimated effects for an observation in the current node, used for pruning
  
  node[['results']] <- effects
  #The current split
  node[['split']] <- temp_split
  #Split the data and create two new subtrees, each using one subset of the data
  if(names(temp_split) %in% names(test_list$categorical)){
    node[['left']] <- build_tree(data[data[names(temp_split)]==temp_split[[1]],],depth = depth+1,
                                 max_depth = max_depth,treatment_list =  treatment_list,target = target,
                                 control = control,test_list = test_list,random = random, n_features = n_features,
                                 criterion = criterion,min_split = min_split,parent_predictions = node[['results']])
    node[['right']] <- build_tree(data[data[names(temp_split)]!=temp_split[[1]],],depth = depth+1,
                                  max_depth = max_depth,treatment_list =  treatment_list,target = target,
                                  control = control,test_list = test_list,random = random, n_features = n_features,
                                  criterion = criterion,min_split = min_split,parent_predictions = node[['results']])
  } else{
    node[['left']] <- build_tree(data = data[data[names(temp_split)]<temp_split[[1]],],depth = depth+1,
                                 max_depth = max_depth,treatment_list =  treatment_list,target = target,
                                 control = control,test_list = test_list,random = random, n_features = n_features,
                                 criterion = criterion,min_split = min_split,parent_predictions = node[['results']])
    node[['right']] <- build_tree(data[data[names(temp_split)]>=temp_split[[1]],],depth = depth+1,
                                  max_depth = max_depth,treatment_list =  treatment_list,target = target,
                                  control = control,test_list = test_list,random = random, n_features = n_features,
                                  criterion = criterion,min_split = min_split,parent_predictions = node[['results']])
  }
  return(node)
}

#Used to create a leaf
final_node <- function(data,treatment_list,target,control){
  treatment_names <- c(treatment_list,control)
  if(nrow(data) == 0){
    effects <- rep(0,length(treatment_names))
  } else{
    effects <- c()
    for(t in treatment_names){
      temp_effect <- mean(data[data[t]==1,target])
      if(is.na(temp_effect)){
        effects <- c(effects,0)
      } else{
        effects <- c(effects,temp_effect)
      }
    }
  }
  names(effects) <- treatment_names
  node <- list()
  node[['type']] <- 'leaf'
  node[['results']] <- effects
  node[['n_samples']] <- nrow(data)
  return(node)
}


#Fucntion to build a random forest
#Two random forest specific parameters:
#n_trees: the number of trees in the forest
#n_features: the number of randomly selected covariates to use for each split 
build_random_forest <- function(train_data,treatment_list,response,control,n_trees,n_features,
                                max_depth = 100, criterion = "simple",min_split=0){
  trees <- list()
  for(x in 1:n_trees){
    print(x)
    set.seed(x)
    temp_train_data <- train_data[sample(nrow(train_data), nrow(train_data),replace = TRUE),]
    temp_tree <- build_tree(data = temp_train_data, depth = 0, max_depth = max_depth, treatment_list = treatment_list, 
                            target = response, control = control, test_list = test_list,
                            random = TRUE,n_features = n_features,criterion = criterion,min_split = min_split)
    trees[[x]] <- temp_tree
  }
  return(trees)
}

#Function to build a random forest with parallelization.
#remain_cores specifies how many cores should NOT be used for the process. 
parallel_build_random_forest <- function(train_data,treatment_list,response,control,n_trees,n_features,
                                         max_depth = 100,remain_cores = 1,criterion = "simple",min_split=0){
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores-remain_cores)
  registerDoParallel(cl)
  trees <- foreach(x=1:n_trees) %dopar% {
    source('ModelImplementations/DecisionTreeImplementation.R')
    set.seed(x)
    temp_train_data <- train_data[sample(nrow(train_data), nrow(train_data),replace = TRUE),]
    temp_tree <- build_tree(data = temp_train_data,0,treatment_list = treatment_list, 
                            test_list = test_list,target = response,control = control,
                            max_depth = max_depth,random = TRUE,n_features =  n_features, 
                            criterion =  criterion,min_split = min_split)
    return(temp_tree)
  }
  stopCluster(cl)
  return(trees)
}

#Pruning ----
#Takes a tree and prunes it with the help of a validation set.
#Returns a pruned tree.
simple_prune_tree <- function(tree, val_data, treatment_list, test_list, target,control,criterion = "simple"){
  new_tree <- assign_val_predictions(tree,val_data,treatment_list,test_list,target,control)
  pruned_tree <- simple_check_pruning(new_tree,val_data,target,control,treatment_list,criterion)
  return(pruned_tree)
}

simple_check_pruning <- function(node,val_data,target,control,treatments,criterion){
  if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
    if(criterion == "max"){
      return(max_pruning_helper(node,treatments,control))
    } else{
      return(simple_pruning_helper(node,treatments,control))
    }
  } else{
    if(node[['left']][['type']] != 'leaf'){
      node[['left']] <- simple_check_pruning(node[['left']],val_data,target,control,treatments,criterion)
    }
    if(node[['right']][['type']] != 'leaf'){
      node[['right']] <- simple_check_pruning(node[['right']],val_data,target,control,treatments,criterion)
    }
    if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
      if(criterion == "max"){
        return(max_pruning_helper(node,treatments,control))
      } else{
        return(simple_pruning_helper(node,treatments,control))
      }
    } else{
      return(node)
    }
  } 
}


#Given a tree and a validation data set this function assigns validation data to the nodes of the tree
#A helper function for the tuning process.
assign_val_predictions <- function(tree,val_data,treatment_list,test_list,target,control){
  treatment_names <- c(treatment_list,control)
  if(nrow(val_data) == 0){
    effects <- rep(0,length(treatment_names))
    names(effects) <- treatment_names
    tree[['val_samples']] <- nrow(val_data)
    for(n in treatment_names){
      tree[[n]] <- 0
    }
    tree[['val_predictions']] <- effects
  } else{
    effects <- c()
    for(t in c(treatment_list,control)){
      temp_effect <- mean(val_data[val_data[t]==1,target])
      if(is.na(temp_effect)){
        effects <- c(effects,0)
      } else{
        effects <- c(effects,temp_effect)
      }
    }
    names(effects) <- treatment_names
    tree[['val_samples']] <- nrow(val_data)
    for(n in treatment_names){
      tree[[n]] <- nrow(val_data[val_data[n] == 1,])
    }
    tree[['val_predictions']] <- effects
  }
  if(tree[['type']] != 'leaf'){
    if(nrow(val_data) == 0){
      tree[['left']] <- assign_val_predictions(tree[['left']],val_data,treatment_list,test_list,target,control)
      tree[['right']] <- assign_val_predictions(tree[['right']],val_data,treatment_list,test_list,target,
                                                control)
    }else{
      split_col <- names(tree[['split']])
      split_value <- tree[['split']][[1]]
      if(split_col %in% names(test_list$categorical)){
        data_left <- val_data[val_data[split_col] == split_value,]
        data_right <- val_data[val_data[split_col] != split_value,]
      } else{
        data_left <- val_data[val_data[split_col] < split_value,]
        data_right <- val_data[val_data[split_col] >= split_value,]
      }
      tree[['left']] <- assign_val_predictions(tree[['left']],data_left,treatment_list,test_list,target,control)
      tree[['right']] <- assign_val_predictions(tree[['right']],data_right,treatment_list,test_list,target,
                                                control)
    }
  }
  return(tree)
}


simple_pruning_helper <- function(node,treatments,control){
  #Check if we are already at the root
  if(node[['type']] == 'root'){
    return(node)
  }
  
  #The rows of the validation set, that ended up in the left and right leaf
  n_left <- node[['left']][['val_samples']]
  n_right <- node[['right']][['val_samples']]
  
  if(n_left == 0 || n_right == 0){
    node[['type']] <- 'leaf'
    return(node)
  }
  
  left_distance <- 0
  left_frac <- n_left/node[['val_samples']]
  left_results <- node[['left']][['val_predictions']]
  for(r in 1:(length(left_results)-1)){
    for(s in (r+1):length(left_results)){
      x <- left_results[r]
      y <- left_results[s]
      left_distance <- left_distance + (x-y)^2
    }
  }
  right_distance <- 0
  right_frac <- n_right/node[['val_samples']]
  right_results <- node[['right']][['val_predictions']]
  for(r in 1:(length(right_results)-1)){
    for(s in (r+1):length(right_results)){
      x <- right_results[r]
      y <- right_results[s]
      right_distance <- right_distance + (x-y)^2
    }
  }
  root_distance <- 0
  root_results <- node[['val_predictions']]
  for(r in 1:(length(root_results)-1)){
    for(s in (r+1):length(root_results)){
      x <- root_results[r]
      y <- root_results[s]
      root_distance <- root_distance + (x-y)^2
    }
  }
  sub_distance <- left_frac*left_distance + right_frac*right_distance
  
  if(is.nan(sub_distance) || is.nan(root_distance) || (sub_distance <= root_distance)){
    node[['type']] <- 'leaf'
    node[['left']] <- NULL
    node[['right']] <- NULL
    node[['split']] <- NULL
    return(node)
  } else{
    return(node)
  }
}

max_pruning_helper <- function(node,treatments,control){
  #Check if we are already at the root
  if(node[['type']] == 'root'){
    return(node)
  }
  
  #The rows of the validation set, that ended up in the left and right leaf
  temp_left_bool <- node[['left']][['val_samples']]
  temp_right_bool <- node[['right']][['val_samples']]
  
  if(temp_left_bool == 0 || temp_right_bool == 0){
    node[['type']] <- 'leaf'
    node[['left']] <- NULL
    node[['right']] <- NULL
    node[['split']] <- NULL
    return(node)
  }
  
  left_distance <- max(node[['left']][['val_predictions']])
  right_distance <- max(node[['right']][['val_predictions']])
  
  root_distance <- max(node[['val_predictions']])
  
  
  if(is.nan(left_distance) || is.nan(right_distance)|| is.nan(root_distance) || 
     (max(left_distance,right_distance) <= root_distance)){
    node[['type']] <- 'leaf'
    node[['left']] <- NULL
    node[['right']] <- NULL
    node[['split']] <- NULL
    return(node)
  } else{
    return(node)
  }
}