#This script contains the implementation of the decision tree following the paper 'Decision trees for uplift modeling
#with single and multiple treatments'

#Importing libraries


#Set up tests
set_up_tests <- function(x,reduce_cases){
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




select_split <- function(a,l,g,divergence,test_list,treatment,control,target,temp_data){
  gain_list <- c()
  name_list <- c()
  for(x in 1:length(test_list$categorical)){
    temp_name <- names(test_list$categorical[x])
    for(y in 1:length(test_list$categorical[[x]])){
      t <- test_list$categorical[[x]][y]
      new_name <- paste(temp_name, as.character(t), sep = '@@')
      gain_list <- c(gain_list,gain(a,l,g,divergence,t,
                                    treatment,control,target,temp_data,'categorical',temp_name))
      name_list <- c(name_list,new_name)
    }
  }
  for(x in 1:length(test_list$numerical)){
    temp_name <- names(test_list$numerical[x])
    for(y in 1:length(test_list$numerical[[x]])){
      t <- test_list$numerical[[x]][y]
      new_name <- paste(temp_name, as.character(t), sep = '@@')
      gain_list <- c(gain_list,gain(a,l,g,divergence,t,
                                    treatment,control,target,temp_data,'numerical',temp_name))
      name_list <- c(name_list,new_name)
    }
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



gain <- function(a,l,g,divergence, test_case,treatment,control,target,temp_data,test_type,test_col){
  if(test_type == 'categorical'){
    if((nrow(temp_data) == 0) || nrow(temp_data[temp_data[test_col]==test_case,]) == 0 ||
       nrow(temp_data[temp_data[test_col]!=test_case,]) == 0 ){
      return(-1)
    }
  } else{
    if((nrow(temp_data) == 0) || nrow(temp_data[temp_data[test_col]<test_case,]) == 0 ||
       nrow(temp_data[temp_data[test_col]!=test_case,]) >= 0 ){
      return(-1)
    }
  }
  
  conditional <- conditional_divergence(a,l,g,divergence,test_case,treatment,control,target,temp_data,
                                        test_type,test_col)
  multiple <- multiple_divergence(a,l,g,divergence,treatment,control,target,temp_data)
  return(conditional-multiple)
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
  return(multiple)
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
  return(div)
}

binary_KL_divergence <- function(x,y){
  p <- sum(x)/length(x)
  q <- sum(y)/length(y)
  temp_result <- (p*log(p/q))+((1-p)*log((1-p)/1-q))
  if(is.nan(temp_result)||is.infinite(temp_result)){
    return(0)
  }
  else{
    return(temp_result)
  }
}

remove_split <- function(temp_test_list,temp_split){
  if(names(temp_split) %in% names(temp_test_list$categorical)){
    temp_test_list$categorical[[names(temp_split)]] <- 
      temp_test_list$categorical[[names(temp_split)]][-match(temp_split[[1]],
                                                           temp_test_list$categorical[[names(temp_split)]])]
  }
  else{
    temp_test_list$numerical[[names(temp_split)]] <- 
      temp_test_list$numerical[[names(temp_split)]][-match(temp_split[[1]],
                                                           temp_test_list$numerical[[names(temp_split)]])]
  }
  return(temp_test_list)
}




create_node <- function(data,depth,max_depth,treatment_list,target,control,test_list, alpha = 1,
                        l = c(0.5,0.5), g = matrix(0.25,nrow = 2, ncol = 2),
                        divergence = 'binary_KL_divergence'){
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
  if(depth == 0){
    node[['type']] <- 'root'
  }
  temp_split <- select_split(alpha,l,g , divergence,test_list = test_list,
                             treatment = treatment_list,control,target,data)
  if(temp_split == -1){
    return(final_node(data,treatment_list,target,control))
  }
  node[['type']] = 'node'
  node[['split']] <- temp_split
  if(names(temp_split) %in% names(test_list$categorical)){
    node[['left']] <- create_node(data[data[names(temp_split)]==temp_split[[1]],],depth = depth+1,max_depth,
                                  treatment_list,target,control,test_list)
    node[['right']] <- create_node(data[data[names(temp_split)]!=temp_split[[1]],],depth = depth+1,max_depth,
                                   treatment_list,target,control,test_list)
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
    effects <- c(effects,mean(data[data[t]==1,target]))
  }
  treatment_names <- c(treatment_names,'control')
  effects <- c(effects,mean(data[data[control]==1,target]))
  names(effects) <- treatment_names
  node <- list()
  node[['type']] <- 'leaf'
  node[['results']] <- effects
  return(node)
}

predict.dt <- function(tree,new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(email)
  type = 'root'
  node = tree
  results = list()
  for(x in 1:nrow(new_data)){
    d = new_data[x,]
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
        if(d[names(split)] <= split[[1]]){
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



####Test Area


treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(email[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE)
test_list$numerical$history <- test_list$numerical$history[1:100]


test_tree <- create_node(email[1:50000,],0,100,treatment_list,'visit','control',test_list)
 

temp_predictions = predict.dt(test_tree,email[50001:64000,])
