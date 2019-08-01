####################################
# Evaluation metrics
####################################

# As described in Zhao et al. 2017
expected_outcome <- function(eval_data){
  # Number of observations
  N <- nrow(eval_data)
  
  # only include points where the assigned treatment equals the predicted
  matching <- eval_data[eval_data$Treatment == eval_data$Assignment , ]
  
  p_c <- nrow(matching[matching$Treatment == 'control', ]) / nrow(eval_data[eval_data$Treatment == 'control', ]) 
  p_t <-  nrow(matching[matching$Treatment != 'control', ]) / nrow(eval_data[eval_data$Treatment != 'control', ])
  
  if(is.na(p_c) | p_c == 0) {
    p_c <- 1
  }
    
  if(is.na(p_t)| p_t == 0) {
    p_t <- 1
  }
  
  result <- (sum(matching$Outcome[matching$Treatment != "control"]) / p_t) + (sum(matching$Outcome[matching$Treatment == "control"]) / p_c )
  # Take AVG Outcome per customer
  result <- result / N
  
  return(result)
}


## Modified Uplift Curve by Zhao
# *Assumption* outcome for each treatment equals prediction of model
expected_percentile_response <- function(predictions){
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # Sum percentiles
  ret <- data.frame(matrix(ncol = 2, nrow = 0))
  
  control_level <- if (nrow(predictions[predictions$Assignment == 'control' , ]) == 0) "Control" else "control"
  
  for(x in seq(0,1, 0.05)){
    # for top x set T to max T
    predictions$T_index <- apply(predictions[, 1:2], 1, which.max)
    # Assign optimal treatment for all
    predictions$Treatment <- colnames(predictions)[predictions$T_index]
    
    
    ###quick fix
    predictions[is.na(predictions)] <- 0
    
    
    # For all who are not in top x Percentile assign Control Treatment
    # For implementation of Matthias lower case control...
    predictions$Treatment[predictions$max_uplift < quantile(predictions$max_uplift, prob=(1-x))] <- control_level
    
    
    # Calculate the Expected Response Value top Percentile
    ret <- rbind(ret, c(x, expected_outcome(predictions)) )
  }
  
  colnames(ret) <- c("Percentile", "Expected Outcome")
  
  return(ret)
}


# *Assumption* always the higher predicted treatment assigned
# Rzepakowski 2012
old_matching_evaluation <- function(predictions, control_level){
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # All treatments and control levels
  treatments <- levels(as.factor(predictions$Treatment))
  
  ret <- data.frame(Percentile=  seq(0,1, 0.05))
  
  # Iterate through treatements
  for(t in treatments){
    group <- predictions[predictions$Treatment == t , ]
    
    N <- nrow(group)
    
    
    ## Do not sort by max uplift, but uplift for this T
    # Sort by max uplift Treatment
    group <- group[order(-group$max_uplift) , ]
    
    outcomes <- c()
    # For each percentile
    for(x in seq(0,1, 0.05)){
      # AVG outcome of top percentile
      outcomes <- c(outcomes, mean(head(group$Outcome, x * N)))
    }     
    ret[ , t] <- outcomes
  }
  
  ret[is.na(ret)] <- 0
  # For each treatment uplift and dynamic curve as max
  curves <- data.frame(Percentile = ret$Percentile)
  
  treatments <- treatments[ !treatments %in% control_level]
  
  # in case never control assigned
  if(!control_level %in% colnames(ret)){
    ret[, control_level] <- 0
  }
  
  for(x in treatments){
    curves[ , x] <- ret[ , x] - ret[ , control_level]
  }
  
  if(ncol(curves) == 2){
    curves[ , "women_treatment"] <- 0
    curves[ , "max T"] <- curves[ , 2 ]
  } else {
    curves[ , "max T"] <- apply(curves[ , c(2 : ncol(curves))], 1, max)
  }
  
  return(curves)
}



## Input Uplift columns
matching_evaluation <- function(predictions, control_level){
  # score for each T individually
  c_group <- predictions[predictions$Assignment == control_level , ]
  N_c <- nrow(c_group)
  
  ret <- data.frame(Percentile=  seq(0,1, 0.05))
  
  
  treatments <- levels(as.factor(predictions$Assignment[predictions$Assignment != control_level]))
  
  for(t in treatments) {
    
    tmp <- predictions[predictions$Assignment == t, ]
    
    # score by uplift column of T
    tmp <- tmp[order(tmp[ , paste0("uplift_",t)] , decreasing = T) , ]
    c_tmp <- c_group[order(c_group[ , paste0("uplift_",t)] , decreasing = T) , ]
    
    N_t  <- nrow(tmp)
    
    outcomes <- c()
    # For each percentile
    for(x in seq(0,1, 0.05)){
      outcomes <- c(outcomes, mean(head(tmp$Outcome, x * N_t)) - mean(head(c_tmp$Outcome, x * N_c)))
    }     
    ret[ , t] <- outcomes
  }
  
  ret[is.na(ret)] <- 0 
  
  return(ret)
}

#new_matching_evaluation(predictions = causal_forest_pred, 'control')


# better ???

## Own naive evaluation approach
# *Assumption* outcome for each treatment equals prediction of model
naive_percentile_response <- function(predictions){
  N <- nrow(predictions)
  
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^Uplift",colnames(predictions))], 1 , max)
  
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # Sort by max uplift Treatment
  predictions <- predictions[order(-predictions$max_uplift) , ]
  
  # Sum percentiles
  ret <- data.frame(matrix(ncol = 2, nrow = 0))
  
  for(x in seq(0,1, 0.1)){
    # Top x % as treated outcome rest have control outcome
    sum_treated <- sum(head(predictions$max_treatment_outcome, N * x) )
    sum_control <- sum(tail(predictions$Control, N * (1 - x)) )
    
    ret <- rbind(ret, c(x, (sum_treated + sum_control) / N))
  }
  
  colnames(ret) <- c("% of Treated", "AVG Outcome")
  
  return(ret)
}


####old

# # As described in Zhao et al. 2017
# expected_outcome <- function(eval_data){
#   # Number of observations
#   N <- nrow(eval_data)
#   
#   ###TODO
#   # change back to Assignment ???
#   # t <- table(eval_data$Assignment)
#   # p_i <- data.frame(t/ nrow(eval_data))
#   # colnames(p_i) <- c("Assignment", "Freq")
#   
#   # only include points where the assigned treatment equals the predicted
#   matching <- eval_data[eval_data$Treatment == eval_data$Assignment , ]
#   #matching <- merge(matching, p_i, all.x = T )
#   
#   p_c <- nrow(matching[matching$Treatment == 'control', ]) / nrow(eval_data[eval_data$Treatment == 'control', ]) 
#   p_t <-  nrow(matching[matching$Treatment != 'control', ]) / nrow(eval_data[eval_data$Treatment != 'control', ])
#   
#   result <- (sum(matching$Outcome[matching$Treatment != "control"]) / p_t) + (sum(matching$Outcome[matching$Treatment == "control"]) / p_c )
#   # Take AVG Outcome per customer
#   result <- result / N
#   
#   #matching[matching$Assignment == 'control', ]
#   
#   # Expected value of response as AVG sum of outcomes / probability of Assignment
#   #res <- (sum(matching$Outcome / matching$Freq) / N)
#   
#   return(result)
# }
