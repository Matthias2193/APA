#This script contains the implementation of the decision tree following the paper 
#'Decision trees for uplift modeling with single and multiple treatments'

#Importing libraries
set.seed(213)

#Set up tests ----
set_up_tests <- function(x,reduce_cases,max_cases = 100){
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
        final_list <- unique(final_list)
        if((length(final_list)>1000)&& reduce_cases){
          final_list <- round(final_list)
        }
        if((length(final_list)>max_cases)&& reduce_cases){
          final_list <- seq(min(final_list), max(final_list),by = (max(final_list)-min(final_list))/max_cases)
        }
        r <- r+1
        s <- s+1
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

select_split <- function(a,l,g,divergence,test_list,treatment,control,target,temp_data,normalize,criterion = 1){
  gain_list <- c()
  name_list <- c()
  if(criterion == 1){
    for(x in 1:length(test_list$categorical)){
      temp_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        t <- test_list$categorical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,gain(a,l,g,divergence,t,
                                      treatment,control,target,temp_data,'categorical',temp_name,normalize))
        name_list <- c(name_list,new_name)
      }
    }
    for(x in 1:length(test_list$numerical)){
      temp_name <- names(test_list$numerical[x])
      for(y in 1:length(test_list$numerical[[x]])){
        t <- test_list$numerical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,gain(a,l,g,divergence,t,
                                      treatment,control,target,temp_data,'numerical',temp_name,normalize))
        name_list <- c(name_list,new_name)
      }
    }
  }
  if(criterion == 2){
    for(x in 1:length(test_list$categorical)){
      temp_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        t <- test_list$categorical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,simple_gain(t,treatment,control,target,temp_data,'categorical',temp_name))
        name_list <- c(name_list,new_name)
      }
    }
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
  if(criterion == 3){
    for(x in 1:length(test_list$categorical)){
      temp_name <- names(test_list$categorical[x])
      for(y in 1:length(test_list$categorical[[x]])){
        t <- test_list$categorical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,simple_gain_median(t,treatment,control,target,temp_data,'categorical',temp_name))
        name_list <- c(name_list,new_name)
      }
    }
    for(x in 1:length(test_list$numerical)){
      temp_name <- names(test_list$numerical[x])
      for(y in 1:length(test_list$numerical[[x]])){
        t <- test_list$numerical[[x]][y]
        new_name <- paste(temp_name, as.character(t), sep = '@@')
        gain_list <- c(gain_list,simple_gain_median(t,treatment,control,target,temp_data,'numerical',temp_name))
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

simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
  treatments <- c(treatment, control)
  gain <- 0
  if(test_type == 'categorical'){
    if((nrow(data) == 0) || nrow(data[data[test_col]==test_case,]) == 0 ||
       nrow(data[data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(data) == 0) || nrow(data[data[test_col]<test_case,]) == 0 ||
       nrow(data[data[test_col]!=test_case,]) >= 0 ){
      return(-1)
    }
  }
  if(test_type == 'categorical'){
    data1 <- data[data[,test_col] == test_case,]
    data2 <- data[data[,test_col] != test_case,]
    for(t in treatments){
      for(s in treatments){
        temp_gain <- (mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2
        gain <- gain + temp_gain
      }
    }
    for(t in treatments){
      gain <- gain * (nrow(data1[data1[,t]==1,]))/nrow(data1)
      gain <- gain * (nrow(data2[data2[,t]==1,]))/nrow(data2)
    }
  }
  if(test_type == 'numerical'){
    data1 <- data[data[,test_col] < test_case,]
    data2 <- data[data[,test_col] >= test_case,]
    temp_gain <- 0
    for(t in treatments){
      for(s in treatments){
        temp_gain <- (mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2
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


gain <- function(a,l,g,divergence, test_case,treatment,control,target,temp_data,test_type,test_col,normalize){
  if(test_type == 'categorical'){
    if((nrow(temp_data) == 0) || nrow(temp_data[temp_data[test_col]==test_case,]) == 0 ||
       nrow(temp_data[temp_data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(temp_data) == 0) || nrow(temp_data[temp_data[test_col]<test_case,]) == 0 ||
       nrow(temp_data[temp_data[test_col]>=test_case,]) == 0 ){
      return(-1)
    }
  }
  
  conditional <- conditional_divergence(a,l,g,divergence,test_case,treatment,control,target,temp_data,
                                        test_type,test_col)
  multiple <- multiple_divergence(a,l,g,divergence,treatment,control,target,temp_data)
  if(normalize){
    normalizer <- Normalization(a,temp_data,control,treatment,target,test_col,test_case,test_type,divergence)
    return((conditional-multiple)/normalizer)
  }
  else{
    return((conditional-multiple))
  }
}

multiple_divergence <- function(a,l,g,divergence,treatments,control,target,temp_data){
  divergence_function <- match.fun(divergence)
  multiple <- 0
  for(t in length(treatments)){
    multiple <- multiple + a*l[t]*divergence_function(temp_data[temp_data[treatments[t]]==1,target],
                                                      temp_data[temp_data[control]==1,target])
    between_treatments <- 0
    for(s in length(treatments)){
      between_treatments <- between_treatments + g[t,s]*divergence_function(
        temp_data[temp_data[treatments[t]]==1,target],temp_data[temp_data[treatments[s]]==1,target])
    }
    multiple <- multiple + (1-a)*between_treatments
  }
  if(is.na(multiple)){
    return(-1)
  }else{
    return(multiple)
  }
}

conditional_divergence <- function(a,l,g,divergence,test_case,treatments,control,target,temp_data,test_type,
                                   test_col){
  div <- 0
  if(test_type == 'categorical'){
    t <- test_case
    div <- div + nrow(temp_data[temp_data[test_col]==t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,temp_data[temp_data[test_col]==t,])
    div <- div + nrow(temp_data[temp_data[test_col]!=t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,temp_data[temp_data[test_col]!=t,])
  }
  else{
    t <- test_case
    div <- div + nrow(temp_data[temp_data[test_col]<t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,
                          temp_data[temp_data[test_col]<t,])
    div <- div + nrow(temp_data[temp_data[test_col]>=t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,
                          temp_data[temp_data[test_col]>=t,])
  }
  if(is.na(div)){
    return(-1)
  }else{
    return(div)
  }
}



#Divergence Measures

binary_KL_divergence <- function(x,y){
  p <- mean(x)
  q <- mean(y)
  temp_result <- (p*log(p/q))+((1-p)*log((1-p)/(1-q)))
  if(is.nan(temp_result)||is.infinite(temp_result)){
    return(0)
  }
  else{
    return(temp_result)
  }
}

EucDistance <- function(x,y){
  return(sqrt((mean(x) - mean(y)) ^ 2))
}

binary_Entropy <- function(prob_vec){
  temp_result <- 0
  for(x in prob_vec){
    temp_result <- temp_result + (x*log(x))
  }
  if(is.nan(temp_result)||is.infinite(temp_result)){
    return(0)
  }
  else{
    return(-temp_result)
  }
}

qini_coef <- function(prob_vec){
  temp_result <- 0
  for(x in prob_vec){
    temp_result <- temp_result + x^2
  }
  return(1-temp_result)
}

#Normalization ----
Normalization <- function(a,temp_data,control,treatments,target,test_col,test_case,test_type,divergence){
  n <- nrow(temp_data)
  nt <- nrow(temp_data[temp_data[,control] != 1,])/n
  nc <- nrow(temp_data[temp_data[,control] == 1,])/n
  divergence_function <- match.fun(divergence)
  if(divergence == "binary_KL_divergence"){
    temp_function <- match.fun("binary_Entropy")
  }else{
    temp_function <- match.fun("qini_coef")
  }
  if(test_type == 'categorical'){
    norm_factor <- a*temp_function(c(nt,nc))*
      divergence_function(nrow(temp_data[(temp_data[,control] != 1) & (temp_data[,test_col] == test_case),])/
                             nrow(temp_data[(temp_data[,control] != 1),]),
                           nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] == test_case),])/
                             nrow(temp_data[(temp_data[,control] == 1),]))
  } else{
    norm_factor <- a*temp_function(c(nt,nc))*
      divergence_function(nrow(temp_data[(temp_data[,control] != 1) & 
                                            (temp_data[,test_col] < test_case),])/
                             nrow(temp_data[(temp_data[,control] != 1),]),
                           nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] < test_case),])/
                             nrow(temp_data[(temp_data[,control] == 1),]))
  }
  
  for(t in treatments){
    nti <-nrow(temp_data[temp_data[,t] == 1,])
    nc <- nrow(temp_data[temp_data[,control] == 1,])
    pi <- nti/(nti+nc)
    pc <- nc/(nti+nc)
    if(test_type == 'categorical'){
      norm_factor <- norm_factor + (1-a) * temp_function(c(pi,pc)) * 
        divergence_function(nrow(temp_data[(temp_data[,t] == 1) & 
                                              (temp_data[,test_col] == test_case),])/
                               nrow(temp_data[(temp_data[,t] == 1),]),
                             nrow(temp_data[(temp_data[,control] == 1) &
                                              (temp_data[,test_col] == test_case),])/
                               nrow(temp_data[(temp_data[,control] == 1),]))
    } else{
      norm_factor <- norm_factor + (1-a) * temp_function(c(pi,pc)) * 
        divergence_function(nrow(temp_data[(temp_data[,t] == 1) & 
                                              (temp_data[,test_col] < test_case),])/
                               nrow(temp_data[(temp_data[,t] == 1),]),
                             nrow(temp_data[(temp_data[,control] == 1) &
                                              (temp_data[,test_col] < test_case),])/
                               nrow(temp_data[(temp_data[,control] == 1),]))
    }
    if(test_type == 'categorical'){
      pti <- nrow(temp_data[(temp_data[,t] == 1) & (temp_data[,test_col] == test_case),])/
        nrow(temp_data[(temp_data[,t] == 1),])
      # pti2 <- nrow(temp_data[(temp_data[,t] == 1) & (temp_data[,test_col] != test_case),])/
      #   nrow(temp_data[(temp_data[,t] == 1),])
    } else{
      pti <- nrow(temp_data[(temp_data[,t] == 1) & (temp_data[,test_col] < test_case),])/
        nrow(temp_data[(temp_data[,t] == 1),])
      # pti2 <- nrow(temp_data[(temp_data[,t] == 1) & (temp_data[,test_col] >= test_case),])/
      #   nrow(temp_data[(temp_data[,t] == 1),])
    }
    if(pti != 0){
      norm_factor <- norm_factor + nti/nrow(temp_data) * temp_function(c(pti,1-pti))
    }
  }
  if(test_type == 'categorical'){
    pc <- nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] == test_case),])/
      nrow(temp_data[(temp_data[,control] == 1),])
    # pc2 <- nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] != test_case),])/
    #   nrow(temp_data[(temp_data[,control] == 1),])
  } else{
    pc <- nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] < test_case),])/
      nrow(temp_data[(temp_data[,control] == 1),])
    # pc2 <- nrow(temp_data[(temp_data[,control] == 1) & (temp_data[,test_col] >= test_case),])/
    #   nrow(temp_data[(temp_data[,control] == 1),])
  }
  if(pc != 0){
    norm_factor <- norm_factor +nc/nrow(temp_data) * temp_function(c(pc,1-pc))
    # norm_factor <- norm_factor +nc/nrow(temp_data) * (-1) * pc2 * log(pc2)
  }
  norm_factor <- norm_factor  + 0.5
  if(norm_factor == 0 || is.na(norm_factor)){
    return(1)
  }
  return(norm_factor)
}

#Functions to build the tree ----

#Tree
#For a continuous target variable use divergence = 'EucDistance'. 
#For a binary categorical use 'binary_KL_divergence'
create_node <- function(data,depth,max_depth,treatment_list,target,control,test_list, alpha = 0.5,
                        l = c(0.5,0.5), g = matrix(0.25,nrow = 2, ncol = 2),
                        divergence = 'binary_KL_divergence',normalize = FALSE, criterion = 1){
  if(depth == max_depth){
    return(final_node(data,treatment_list,target,control))
  }
  for(t in treatment_list){
    if(nrow(data[data[t]==1,]) == 0){
      return(final_node(data,treatment_list,target,control))
    }
  }
  if(nrow(data[data[control]==1,]) == 0){
    return(final_node(data,treatment_list,target,control))
  }
  node <- list()
  temp_split <- select_split(alpha,l,g , divergence,test_list = test_list,
                             treatment = treatment_list,control,target,data,normalize,criterion)
  if(temp_split == -1){
    return(final_node(data,treatment_list,target,control))
  }
  node[['type']] = 'node'
  if(depth == 0){
    node[['type']] <- 'root'
  }
  node[['n_samples']] <- nrow(data)
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
  node[['split']] <- temp_split
  if(names(temp_split) %in% names(test_list$categorical)){
    node[['left']] <- create_node(data[data[names(temp_split)]==temp_split[[1]],],depth = depth+1,max_depth,
                                  treatment_list,target,control,test_list,alpha,l,g,divergence,normalize,
                                  criterion = criterion)
    node[['right']] <- create_node(data[data[names(temp_split)]!=temp_split[[1]],],depth = depth+1,max_depth,
                                   treatment_list,target,control,test_list,alpha,l,g,divergence,normalize,
                                   criterion = criterion)
  }
  else{
    node[['left']] <- create_node(data[data[names(temp_split)] < temp_split[[1]],],depth = depth+1,max_depth,
                                  treatment_list,target,control,test_list)
    node[['right']] <- create_node(data[data[names(temp_split)] >= temp_split[[1]],],depth = depth+1,max_depth,
                                   treatment_list,target, control,test_list)
  }
  return(node)
}

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
  treatment_names <- c(treatment_names,'control')
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
                         criterion,pruning,divergence = "binary_KL_divergence"){
  trees <- list()
  retain_cols <- c(treatment_list,control,response)
  sample_cols <- setdiff(colnames(train_data),retain_cols)
  for(x in 1:n_trees){
    temp_cols <- sample(sample_cols,n_features,replace = F)
    chosen_cols <- c(temp_cols,retain_cols)
    test_list <- set_up_tests(train_data[,chosen_cols],TRUE)
    temp_tree <- create_node(data = train_data[,chosen_cols],0,100,treatment_list = treatment_list, 
                             test_list = test_list, criterion = criterion,target = response,control = control,
                             divergence = divergence)
    if(pruning){
      temp_prune_tree <- prune_tree(temp_tree,val_data[,chosen_cols], treatment_list, test_list, response, control)
      trees[[x]] <- temp_prune_tree
    } else{
      trees[[x]] <- temp_tree
    }
    
  }
  return(trees)
}

#Pruning ----
#Takes a tree and prunes it with the help of a validation set.
#Returns a pruned tree.

prune_tree <- function(tree, val_data, treatment_list, test_list, target,control){
  new_tree <- assign_val_predictions(tree,val_data,treatment_list,test_list,target,control)
  pruned_tree <- check_pruning(new_tree,val_data,target,control,treatment_list)
  return(pruned_tree)
}

check_pruning <- function(node,val_data,target,control,treatments){
  if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
    return(pruning_helper(node,treatments))
  } else{
    if(node[['left']][['type']] != 'leaf'){
      node[['left']] <- check_pruning(node[['left']],val_data,target,control)
    }
    if(node[['right']][['type']] != 'leaf'){
      node[['right']] <- check_pruning(node[['right']],val_data,target,control)
    }
    if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
      return(pruning_helper(node,treatments))
    } else{
      return(node)
    }
  } 
}

pruning_helper <- function(node,treatments){
  #Check if we are already at the root
  if(node[['type']] == 'root'){
    return(node)
  }
  
  #Left
  temp_pred_left <- node[['left']][['results']]
  temp_left <- node[['left']][['n_samples']]
  best_treatment_left <- names(temp_pred_left[match(max(temp_pred_left[treatments]),temp_pred_left)])
  left_sign <- sign(temp_pred_left[[best_treatment_left]]-temp_pred_left[['control']])
  
  #Right
  temp_pred_right <- node[['right']][['results']]
  temp_right <- node[['right']][['n_samples']]
  best_treatment_right <- names(temp_pred_right[match(max(temp_pred_right[treatments]),temp_pred_right)])
  right_sign <- sign(temp_pred_right[[best_treatment_right]]-temp_pred_right[['control']])
  
  #Root
  temp_pred_root <- node[['results']]
  best_treatment_root <- names(temp_pred_root[match(max(temp_pred_root[treatments]),temp_pred_root)])
  root_sign <- sign(temp_pred_root[[best_treatment_root]]-temp_pred_root[['control']])
  
  
  #The rows of the validation set, that ended up in the left and right leaf
  temp_left_bool <- node[['left']][['val_samples']]
  temp_right_bool <- node[['right']][['val_samples']]
  
  if(temp_left_bool == 0 || temp_right_bool == 0){
    node[['type']] <- 'leaf'
    node[['n_samples']] <- node[['left']][['n_samples']][1] + node[['right']][['n_samples']][1]
    return(node)
  }
  
  n_treat_left <- node[['left']][[best_treatment_left]][1]
  n_control_left <- node[['left']][[control]][1]
  prob_treatment_left <- node[['left']][["val_predictions"]][[best_treatment_left]][1]
  prob_control_left <- node[['left']][["val_predictions"]][[control]][1]
  
  n_treat_right <- node[['right']][[best_treatment_right]][1]
  n_control_right <- node[['right']][[control]][1]
  prob_treatment_right <- node[['right']][["val_predictions"]][[best_treatment_right]][1]
  prob_control_right <- node[['right']][["val_predictions"]][[control]][1]
  
  n_treat_root <- node[[best_treatment_root]][1]
  n_control_root <-  node[[control]][1]
  prob_treatment_root <- node[["val_predictions"]][[best_treatment_root]][1]
  prob_control_root <- node[["val_predictions"]][[control]][1]
  
  
  d1 <-  (n_treat_left+n_control_left)/(n_treat_root+n_control_root)*
    left_sign*(prob_treatment_left-prob_control_left)
  d1 <-  d1 + (n_treat_right+n_control_right)/(n_treat_root+n_control_root)*
    right_sign*(prob_treatment_right-prob_control_right)
  
  d2 <- root_sign*(prob_treatment_root-prob_control_root)
  
  if(is.nan(d1) || is.nan(d2)){
    return(node)
  }
  if(d1 <= d2){
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


#Test Area----

# email <- read.csv('Email.csv')
# 
# email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
# email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
# email$control <- ifelse(email$segment=='No E-Mail',1,0)
# email$mens <- as.factor(email$mens)
# email$womens <- as.factor(email$womens)
# email$newbie <- as.factor(email$newbie)
# 
# email$visit <- email$spend <- email$segment <- NULL
# 
# response <- 'conversion'
# control <- 'control'
# 
# treatment_list <- c('men_treatment','women_treatment')
# test_list <- set_up_tests(email[,c("recency","history_segment","history","mens","womens","zip_code",
#                                    "newbie","channel")],TRUE)
# idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)
# 
# train <- email[-idx, ]
# 
# test <- email[idx, ]
# 
# test_tree <- create_node(email[1:50000,],0,100,treatment_list,'conversion','control',test_list,
#                          normalize  = TRUE)
# 
# start_time <- Sys.time()
# pruned_tree <- prune_tree_old(test_tree,email[50001:64000,],email[1:50000,],target = response)
# prune1_time <- Sys.time()
# pruned_tree2 <- prune_tree(test_tree,email[50001:64000,], treatment_list, test_list, response, control)
# prune2_time <- Sys.time()
# 
# time_old <- prune1_time - start_time
# time_new <- prune2_time - prune1_time
# print(time_old)
# print(time_old_updated)

#Old´Functions----
prune_tree_old <- function(tree, val_data, train_data, target){
  val_pred <- predict.dt(tree, val_data)
  train_pred <- predict.dt(tree, train_data)
  pruned_tree <- check_pruning_old(tree,train_pred,val_pred,val_data,target)
  while(check_tree_changes(pruned_tree)){
    pruned_tree <- type_subtrees(pruned_tree)
    val_pred <- predict.dt(pruned_tree, val_data)
    train_pred <- predict.dt(pruned_tree, train_data)
    pruned_tree <- check_pruning_old(pruned_tree,train_pred,val_pred,val_data,target)
  }
  return(pruned_tree)
}


check_pruning_old <- function(node,predictions,predictions_val,val_data,target){
  if(node[['left']][['type']] == 'leaf' && node[['right']][['type']] == 'leaf'){
    #Left
    if(node[['type']] == 'root'){
      return(node)
    }
    temp_pred_left <- node[['left']][['results']]
    temp_left <- plyr::compact(lapply(predictions, function(x) if(sum(x == temp_pred_left) == 3){x}))
    best_treatment_left <- names(temp_pred_left[match(max(temp_pred_left[1:2]),temp_pred_left)])
    left_sign <- sign(temp_pred_left[[best_treatment_left]]-temp_pred_left[['control']])
    
    #Right
    temp_pred_right <- node[['right']][['results']]
    temp_right <- plyr::compact(lapply(predictions, function(x) if(sum(x == temp_pred_right) == 3){x}))
    best_treatment_right <- names(temp_pred_right[match(max(temp_pred_right[1:2]),temp_pred_right)])
    right_sign <- sign(temp_pred_right[[best_treatment_right]]-temp_pred_right[['control']])
    
    #Root
    temp_pred_root <- (length(temp_left)*temp_left[[1]]+length(temp_right)*temp_right[[1]])/(length(temp_left)+
                                                                                               length(temp_right))
    best_treatment_root <- names(temp_pred_root[match(max(temp_pred_root[1:2]),temp_pred_root)])
    root_sign <- sign(temp_pred_root[[best_treatment_root]]-temp_pred_root[['control']])
    
    
    #The rows of the validation set, that ended up in the left and right leaf
    temp_left_bool <- unlist(lapply(predictions_val, function(x) (sum(x == temp_pred_left) == 3)))
    temp_right_bool <- unlist(lapply(predictions_val, function(x) (sum(x == temp_pred_right) == 3)))
    
    if(sum(temp_left_bool) == 0 || sum(temp_right_bool) == 0){
      result_node <- list()
      result_node[['type']] <- 'leaf'
      result_node[['results']] <- temp_pred_root
      predictions <- lapply(predictions, function(x)if(sum(x == temp_pred_right) == 3 || 
                                                       sum(x == temp_pred_left) == 3){temp_pred_root}
                            else{x})
      return(result_node)
    }
    
    temp_val <- val_data[temp_left_bool,]
    n_treat_left <- nrow(temp_val[temp_val[best_treatment_left] == 1,])
    n_control_left <- nrow(temp_val[temp_val['control'] == 1,])
    prob_treatment_left <- mean(temp_val[temp_val[best_treatment_left]== 1,target])
    prob_control_left <- mean(temp_val[temp_val['control'] == 1,target])
    
    temp_val <- val_data[temp_right_bool,]
    n_treat_right <- nrow(temp_val[temp_val[best_treatment_right] == 1,])
    n_control_right <- nrow(temp_val[temp_val[best_treatment_right] == 1,])
    prob_treatment_right <- mean(temp_val[temp_val[best_treatment_right] == 1,target])
    prob_control_right <- mean(temp_val[temp_val['control']== 1,target])
    
    
    temp_val <- val_data[(temp_left_bool+temp_right_bool) > 0,]
    n_treat_root <- nrow(temp_val[temp_val[best_treatment_root] == 1,])
    n_control_root <- nrow(temp_val[temp_val[best_treatment_root] == 1,])
    prob_treatment_root <- mean(temp_val[temp_val[best_treatment_root] == 1,target])
    prob_control_root <- mean(temp_val[temp_val['control'] == 1,target])
    
    
    d1 <-  (n_treat_left+n_control_left)/(n_treat_root+n_control_root)*
      left_sign*(prob_treatment_left-prob_control_left)
    d1 <-  d1 + (n_treat_right+n_control_right)/(n_treat_root+n_control_root)*
      right_sign*(prob_treatment_right-prob_control_right)
    
    d2 <- root_sign*(prob_treatment_root-prob_control_root)
    
    if(is.nan(d1) || is.nan(d2)){
      return(node)
    }
    if(d1 <= d2){
      result_node <- list()
      result_node[['type']] <- 'leaf'
      result_node[['results']] <- temp_pred_root
      predictions <- lapply(predictions, function(x)if(sum(x == temp_pred_right) == 3 || 
                                                       sum(x == temp_pred_left) == 3){temp_pred_root}
                            else{x})
      predictions_val <- lapply(predictions_val, function(x)if(sum(x == temp_pred_right) == 3 || 
                                                               sum(x == temp_pred_left) == 3){temp_pred_root}
                                else{x})
      return(result_node)
    } else{
      return(node)
    }
  } else{
    if(node[['left']][['type']] != 'leaf'){
      node[['left']] <- check_pruning_old(node[['left']],predictions,predictions_val,val_data,target)
    }
    if(node[['right']][['type']] != 'leaf'){
      node[['right']] <- check_pruning_old(node[['right']],predictions,predictions_val,val_data,target)
    }
  }
  return(node)
}

#Function that goes to the tree to identify subtrees. Used for pruning

type_subtrees <- function(tree){
  if(tree[['left']][['type']] == 'leaf' && tree[['right']][['type']] == 'leaf'){
    if(tree[['type']] != 'root'){
      tree[['type']] <- 'sub'
    }
  } else{
    if(tree[['left']][['type']] != 'leaf'){
      tree[['left']] <- type_subtrees(tree[['left']])
    }
    if(tree[['right']][['type']] != 'leaf'){
      tree[['right']] <- type_subtrees(tree[['right']])
    }
  }
  return(tree)
}


#Helper function for pruning

check_tree_changes <- function(tree){
  result <- FALSE
  if(tree[['left']][['type']] == 'leaf' && tree[['right']][['type']] == 'leaf' && tree[['type']] != 'sub'&& 
     tree[['type']] != 'root'){
    return(TRUE)
  } else{
    if(tree[['left']][['type']] != 'leaf'){
      if(check_tree_changes(tree[['left']]) == TRUE){
        return(TRUE)
      }
    }
    if(tree[['right']][['type']] != 'leaf'){
      if(check_tree_changes(tree[['right']]) == TRUE){
        return(TRUE)
      }
    }
  }
  return(result)
}
