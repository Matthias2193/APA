#This script contains the implementation of the decision tree following the paper 
#'Decision trees for uplift modeling with single and multiple treatments'

#Importing libraries
set.seed(213)
library(parallel)
library(foreach)
library(doParallel)
source("PredictionFunctions.R")

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
      if((length(final_list)>1000)&& reduce_cases){
        final_list <- round(final_list)
      }
      if((length(final_list)>max_cases)&& reduce_cases){
        final_list <- quantile(x[,n],seq(0,1,1/max_cases))
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
#gain > 0, then -1 is returned indicating that no split is beneficial
select_split <- function(test_list,treatment,control,target,temp_data){
  gain_list <- c()
  name_list <- c()
  if(length(test_list$categorical)>0){
    for(x in 1:length(test_list$categorical)){
      temp_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        t <- test_list$categorical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,simple_gain(t,treatment,control,target,temp_data,'categorical',temp_name))
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
        gain_list <- c(gain_list,simple_gain(t,treatment,control,target,temp_data,'numerical',temp_name))
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
    if(is.na(as.numeric(temp_result[[1]][2]))){
      result <- c(temp_result[[1]][2])
      names(result) <- temp_result[[1]][1]
      return(result)
    }
    else{
      result <- c(as.numeric(temp_result[[1]][2]))
      names(result) <- temp_result[[1]][1]
      return(result)
    }
  }
  else{
    return(-1)
  }
}

#This method calculates the gain for a given split
simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    if((nrow(data) == 0) || nrow(data[data[test_col]==test_case,]) == 0 ||
       nrow(data[data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(data) == 0) || nrow(data[data[test_col]<test_case,]) == 0 ||
       nrow(data[data[test_col]>=test_case,]) == 0 ){
      return(-1)
    }
  }
  current_gain <- 0
  for(t in treatment){
    if(mean(data[data[,t]==1,target])>current_gain){
      current_gain <- mean(data[data[,t]==1,target])
    }
  }
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  if(test_type == 'categorical'){
    #First the data is split according to the given split
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
    frac1 <- nrow(data1)/nrow(data)
    frac2 <- nrow(data2)/nrow(data)
    #Here the gain is calculated
    left_gain <- 0
    right_gain <- 0
    for(t in treatment){
      left_gain <- max(left_gain,mean(data1[data1[,t]==1,target]))
      right_gain <- max(right_gain,mean(data2[data2[,t]==1,target]))
    }
    gain <- max(left_gain,right_gain)
    for(t in treatments){
      if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
        gain <- 0
      }
    }
  } else{
    #The same as above, but for numerical covariates
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
    frac1 <- nrow(data1)/nrow(data)
    frac2 <- nrow(data2)/nrow(data)
    left_gain <- 0
    right_gain <- 0
    for(t in treatment){
      left_gain <- max(left_gain,mean(data1[data1[,t]==1,target]))
      right_gain <- max(right_gain,mean(data2[data2[,t]==1,target]))
    }
    gain <- (frac1*left_gain+frac2*right_gain)
    # for(t in treatments){
    #   if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
    #     gain <- 0
    #   }
    # }
  }
  if(is.na(gain)){
    gain = -1
  }
  if(gain <= current_gain){
    gain = -1
  }
  return(gain)
}



new_simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    if((nrow(data) == 0) || nrow(data[data[test_col]==test_case,]) == 0 ||
       nrow(data[data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(data) == 0) || nrow(data[data[test_col]<test_case,]) == 0 ||
       nrow(data[data[test_col]>=test_case,]) == 0 ){
      return(-1)
    }
  }
  current_gain <- 0
  while(x < length(treatments)){
    t <- treatments[x]
    s <- treatments[x+1]
    temp_gain <- (mean(data[data[,t] == 1,target])-mean(data[data[,s] == 1,target]))^2
    current_gain <- current_gain + temp_gain
    x <- x+1
  }
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  if(test_type == 'categorical'){
    #First the data is split according to the given split
    data1 <- data[data[,test_col] == test_case,]
    frac1 <- nrow(data1)/nrow(data)
    data2 <- data[data[,test_col] != test_case,]
    frac2 <- nrow(data2)/nrow(data)
    #Here the gain is calculated
    x <- 1
    while(x < length(treatments)){
      t <- treatments[x]
      s <- treatments[x+1]
      temp_gain <- frac1*(mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2 +
        frac2*(mean(data2[data2[,t] == 1,target])-mean(data2[data2[,s] == 1,target]))^2
      gain <- gain + temp_gain
      x <- x+1
    }
    #Make sure that there are data points of each treatment in each subset of the data
    # for(t in treatments){
    #   if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
    #     gain <- 0
    #   }
    # }
  } else{
    #The same as above, but for numerical covariates
    data1 <- data[data[,test_col] < test_case,]
    frac1 <- nrow(data1)/nrow(data)
    data2 <- data[data[,test_col] >= test_case,]
    frac2 <- nrow(data2)/nrow(data)
    x <- 1
    while(x < length(treatments)){
      t <- treatments[x]
      s <- treatments[x+1]
      temp_gain <- frac1*(mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2 +
        frac2*(mean(data2[data2[,t] == 1,target])-mean(data2[data2[,s] == 1,target]))^2
      gain <- gain + temp_gain
      x <- x+1
    }
    # for(t in treatments){
    #   if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
    #     gain <- 0
    #   }
    # }
  }
  if(is.na(gain)){
    gain = -1
  }
  if(gain <= current_gain){
    gain = -1
  }
  return(gain)
}

old_simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
  treatments <- c(treatment, control)
  gain <- 0
  #First check if there is data in each subset after the data is split. If not return -1.
  if(test_type == 'categorical'){
    if((nrow(data) == 0) || nrow(data[data[test_col]==test_case,]) == 0 ||
       nrow(data[data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(data) == 0) || nrow(data[data[test_col]<test_case,]) == 0 ||
       nrow(data[data[test_col]>=test_case,]) == 0 ){
      return(-1)
    }
  }
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  if(test_type == 'categorical'){
    #First the data is split according to the given split
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
    #Here the gain is calculated
    for(t in treatments){
      for(s in treatments){
        temp_gain <- (mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2 +
          (mean(data2[data2[,t] == 1,target])-mean(data2[data2[,s] == 1,target]))^2
        gain <- gain + temp_gain
      }
    }
    #Make sure that there are data points of each treatment in each subset of the data
    for(t in treatments){
      if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
        gain <- 0
      }
    }
  } else{
    #The same as above, but for numerical covariates
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
    for(t in treatments){
      for(s in treatments){
        temp_gain <- (mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2 +
          (mean(data2[data2[,t] == 1,target])-mean(data2[data2[,s] == 1,target]))^2
        gain <- gain + temp_gain
      }
    }
    for(t in treatments){
      gain <- gain * (nrow(data1[data1[,t]==1,]))/nrow(data1)
      gain <- gain * (nrow(data2[data2[,t]==1,]))/nrow(data2)
    }
  }
  if(is.na(gain)){
    gain = -1
  }
  return(gain)
}







#Functions to build the tree ----

#Tree
#For a continuous target variable use divergence = 'EucDistance'. 
#For a binary categorical use 'binary_KL_divergence'
#depth: the current depth, can be ignored by the user, for internal pruposes
#treatment_list: a list with the names of all treatments
#target: the name of the response variable
#control: the name of the control 'treatment'
#test_list: a list of possible splits created by the 'set_up_tests' function
#criterion: 1 for Rzp-tree, 2 for simple tree
#alpha, l and g are parameters according to Rzepakowski paper (only necessare if criterion = 1)
build_tree <- function(data,depth,max_depth,treatment_list,target,control,test_list){
  #Return leaf if current depth is max depth
  if(depth == max_depth){
    return(final_node(data,treatment_list,target,control))
  }
  #Create current node
  node <- list()
  #Select split with maximum gain
  temp_split <- select_split(test_list = test_list, treatment = treatment_list, control, target,data)
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
  treatment_names <- c()
  effects <- c()
  for(t in treatment_list){
    treatment_names <- c(treatment_names,t)
    effects <- c(effects,mean(data[data[t]==1,target]))
  }
  treatment_names <- c(treatment_names,'control')
  effects <- c(effects,mean(data[data[control]==1,target]))
  names(effects) <- treatment_names
  node[['results']] <- effects
  #The current split
  node[['split']] <- temp_split
  #Split the data and create two new subtrees, each using one subset of the data
  if(names(temp_split) %in% names(test_list$categorical)){
    node[['left']] <- build_tree(data[data[names(temp_split)]==temp_split[[1]],],depth = depth+1,max_depth,
                                  treatment_list,target,control,test_list)
    node[['right']] <- build_tree(data[data[names(temp_split)]!=temp_split[[1]],],depth = depth+1,max_depth,
                                   treatment_list,target,control,test_list)
  }
  else{
    node[['left']] <- build_tree(data[data[names(temp_split)]<temp_split[[1]],],depth = depth+1,max_depth,
                                  treatment_list,target,control,test_list)
    node[['right']] <- build_tree(data[data[names(temp_split)]>=temp_split[[1]],],depth = depth+1,max_depth,
                                   treatment_list,target,control,test_list)
  }
  return(node)
}

#Used to creat a leaf
final_node <- function(data,treatment_list,target,control){
  treatment_names <- c()
  effects <- c()
  for(t in treatment_list){
    treatment_names <- c(treatment_names,t)
    temp_effect <- mean(data[data[t]==1,target])
    if(is.na(temp_effect)){
      effects <- c(effects,0)
    } else{
      effects <- c(effects,temp_effect)
    }
    
  }
  treatment_names <- c(treatment_names,control)
  temp_effect <- mean(data[data[control]==1,target])
  if(is.na(temp_effect)){
    effects <- c(effects,0)
  } else{
    effects <- c(effects,temp_effect)
  }
  names(effects) <- treatment_names
  node <- list()
  node[['type']] <- 'leaf'
  node[['results']] <- effects
  node[['n_samples']] <- nrow(data)
  return(node)
}

#Forest
build_forest <- function(train_data, val_data,treatment_list,response,control,n_trees,n_features,
                         criterion,pruning,max_depth = 10){
  trees <- list()
  retain_cols <- c(treatment_list,control,response)
  sample_cols <- setdiff(colnames(train_data),retain_cols)
  for(x in 1:n_trees){
    temp_cols <- sample(sample_cols,n_features,replace = F)
    chosen_cols <- c(temp_cols,retain_cols)
    test_list <- set_up_tests(train_data[,chosen_cols],TRUE)
    temp_tree <- build_tree(data = train_data[,chosen_cols],0,treatment_list = treatment_list, 
                             test_list = test_list, criterion = criterion,target = response,control = control,
                             max_depth = max_depth)
    if(pruning){
      temp_prune_tree <- simple_prune_tree(temp_tree,val_data[,chosen_cols], treatment_list, test_list, response, control)
      trees[[x]] <- temp_prune_tree
    } else{
      trees[[x]] <- temp_tree
    }
    
  }
  return(trees)
}

parallel_build_forest <- function(train_data, val_data,treatment_list,response,control,n_trees,n_features,pruning,max_depth = 10){
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores-1)
  registerDoParallel(cl)
  retain_cols <- c(treatment_list,control,response)
  sample_cols <- setdiff(colnames(train_data),retain_cols)
  trees <- foreach(x=1:n_trees) %dopar% {
    source('DecisionTreeImplementation.R')
    set.seed(x)
    temp_cols <- sample(sample_cols,n_features,replace = F)
    chosen_cols <- c(temp_cols,retain_cols)
    test_list <- set_up_tests(train_data[,chosen_cols],TRUE)
    temp_tree <- build_tree(data = train_data[,chosen_cols],0,treatment_list = treatment_list, 
                            test_list = test_list,target = response,control = control,
                            max_depth = max_depth)
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

#Pruning ----
#Takes a tree and prunes it with the help of a validation set.
#Returns a pruned tree.

simple_prune_tree <- function(tree, val_data, treatment_list, test_list, target,control){
  new_tree <- assign_val_predictions(tree,val_data,treatment_list,test_list,target,control)
  pruned_tree <- simple_check_pruning(new_tree,val_data,target,control,treatment_list)
  return(pruned_tree)
}

simple_check_pruning <- function(node,val_data,target,control,treatments){
  if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
    return(simple_pruning_helper(node,treatments,control))
  } else{
    if(node[['left']][['type']] != 'leaf'){
      node[['left']] <- simple_check_pruning(node[['left']],val_data,target,control,treatments)
    }
    if(node[['right']][['type']] != 'leaf'){
      node[['right']] <- simple_check_pruning(node[['right']],val_data,target,control,treatments)
    }
    if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
      return(simple_pruning_helper(node,treatments,control))
    } else{
      return(node)
    }
  } 
}



assign_val_predictions <- function(tree,val_data,treatment_list,test_list,target,control){
  if(nrow(val_data) == 0){
    treatment_names <- c(treatment_list,control)
    effects <- rep(0,length(treatment_names))
    names(effects) <- treatment_names
    tree[['val_samples']] <- nrow(val_data)
    for(n in treatment_names){
      tree[[n]] <- 0
    }
    tree[['val_predictions']] <- effects
  } else{
    treatment_names <- c()
    effects <- c()
    for(t in treatment_list){
      treatment_names <- c(treatment_names,t)
      temp_effect <- mean(val_data[val_data[t]==1,target])
      if(is.na(temp_effect)){
        effects <- c(effects,0)
      } else{
        effects <- c(effects,temp_effect)
        }
    }
    treatment_names <- c(treatment_names,control)
    temp_effect <- mean(val_data[val_data[control]==1,target])
    if(is.na(temp_effect)){
      effects <- c(effects,0)
    } else{
      effects <- c(effects,temp_effect)
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
  temp_left_bool <- node[['left']][['val_samples']]
  temp_right_bool <- node[['right']][['val_samples']]
  
  if(temp_left_bool == 0 || temp_right_bool == 0){
    node[['type']] <- 'leaf'
    node[['n_samples']] <- node[['left']][['n_samples']][1] + node[['right']][['n_samples']][1]
    return(node)
  }
  
  left_distance <- 0
  for(r in node[['left']][['val_predictions']]){
    for(s in node[['left']][['val_predictions']]){
      left_distance <- left_distance + (r-s)^2
    }
  }
  right_distance <- 0
  for(r in node[['right']][['val_predictions']]){
    for(s in node[['right']][['val_predictions']]){
      right_distance <- right_distance + (r-s)^2
    }
  }
  root_distance <- 0
  for(r in node[['val_predictions']]){
    for(s in node[['val_predictions']]){
      root_distance <- root_distance + (r-s)^2
    }
  }
  sub_distance <- (right_distance+left_distance)/2
  for (r in c(treatments,control)) {
    sub_distance <- sub_distance*(node[['left']][[r]]+node[['right']][[r]])/
      (node[['left']][['val_samples']]+node[['right']][['val_samples']])
    root_distance <- root_distance*node[[r]]/node[['val_samples']]
  }
  
  
  if(is.nan(sub_distance) || is.nan(root_distance)){
    return(node)
  }
  if(sub_distance <= root_distance){
    node[['type']] <- 'leaf'
    node[['n_samples']] <- node[['left']][['n_samples']][1] + node[['right']][['n_samples']][1]
    node[['left']] <- NULL
    node[['right']] <- NULL
    node[['split']] <- NULL
    return(node)
  } else{
    return(node)
  }
}


new_simple_pruning_helper <- function(node,treatments,control){
  #Check if we are already at the root
  if(node[['type']] == 'root'){
    return(node)
  }
  
  #The rows of the validation set, that ended up in the left and right leaf
  temp_left_bool <- node[['left']][['val_samples']]
  temp_right_bool <- node[['right']][['val_samples']]
  
  if(temp_left_bool == 0 || temp_right_bool == 0){
    node[['type']] <- 'leaf'
    node[['n_samples']] <- node[['left']][['n_samples']][1] + node[['right']][['n_samples']][1]
    return(node)
  }
  
  left_distance <- max(node[['left']][['val_predictions']])
  right_distance <- max(node[['right']][['val_predictions']])

  root_distance <- max(node[['val_predictions']])
  
  
  if(is.nan(left_distance) || is.nan(right_distance)|| is.nan(root_distance)){
    return(node)
  }
  if(max(left_distance,right_distance) <= root_distance){
    node[['type']] <- 'leaf'
    node[['n_samples']] <- node[['left']][['n_samples']][1] + node[['right']][['n_samples']][1]
    node[['left']] <- NULL
    node[['right']] <- NULL
    node[['split']] <- NULL
    return(node)
  } else{
    return(node)
  }
}


