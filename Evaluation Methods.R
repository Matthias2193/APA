####################################
# Evaluation metrics
####################################

##TODO
# is a second upscaling necessary?

# As described in Zhao et al. 2017
expected_outcome <- function(eval_data){
  # Number of observations
  N <- nrow(eval_data)
  
  # Frequencies of treatments z_i
  t <- table(eval_data$Assignment)
  p_i <- data.frame(t/ nrow(eval_data))
  colnames(p_i) <- c("Assignment", "Freq")
  
  # only include points where the assigned treatment equals the predicted
  matching <- eval_data[eval_data$Treatment == eval_data$Assignment , ]
  matching <- merge(matching, p_i, all.x = T )
  
  ##TODO
  # print how many have been matched for evaluation
  # print(paste("Total ",N))
  #print(aggregate(Outcome~Assignment,matching, length ))
  
  
  matching[matching$Assignment == 'control', ]
  
  # Expected value of response as AVG sum of outcomes / probability of Assignment
  res <- (sum(matching$Outcome / matching$Freq) / N)
  
  return(res)
}


# scaled_expected_outcome <- function(eval_data){
#   # Number of observations
#   N <- nrow(eval_data)
#   
#   # Frequencies of treatments z_i
#   t <- table(eval_data$Assignment)
#   p_i <- data.frame(t/ nrow(eval_data))
#   colnames(p_i) <- c("Assignment", "Freq")
#   
#   # only include points where the assigned treatment equals the predicted
#   matching <- eval_data[eval_data$Treatment == eval_data$Assignment , ]
#   
#   # also count how many have been matched...
#   # number matched for each T / this T in eval_data
#   aggregate(Outcome~Treatment, eval_data, length)
#   
#   
#   matching <- merge(matching, p_i, all.x = T )
#   
#   # Expected value of response as AVG sum of outcomes / probability of Assignment
#   res <- (sum(matching$Outcome / matching$Freq) / N)
#   
#   return(res)
# }

## Modified Uplift Curve by Zhao
# *Assumption* outcome for each treatment equals prediction of model
expected_percentile_response <- function(predictions){
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^Uplift",colnames(predictions))], 1 , max)
  
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # Sum percentiles
  ret <- data.frame(matrix(ncol = 2, nrow = 0))
  
  control_level <- if (nrow(predictions[predictions$Assignment == 'control' , ]) == 0) "Control" else "control"
  
  for(x in seq(0,1, 0.05)){
    # for top x set T to max T
    predictions$T_index <- apply(predictions[, 1:2], 1, which.max)
    # Assign optimal treatment for all
    predictions$Treatment <- colnames(predictions)[predictions$T_index]
    
    # For all who are not in top x Percentile assign Control Treatment
    # For implementation of Matthias lower case control...
    predictions$Treatment[predictions$max_uplift < quantile(predictions$max_uplift,prob=(1-x))] <- control_level
    
    # Calculate the Expected Response Value top Percentile
    ret <- rbind(ret, c(x, expected_outcome(predictions)) )
  }
  
  colnames(ret) <- c("Percentile", "Expected Outcome")
  
  return(ret)
}


# *Assumption* always the higher predicted treatment assigned
# Rzepakowski 2012
matching_evaluation <- function(predictions, control_level){
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^Uplift",colnames(predictions))], 1 , max)
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # All treatments and control levels
  treatments <- levels(as.factor(predictions$Treatment))
  
  ret <- data.frame(Percentile=  seq(0,1, 0.05))
  
  # Iterate through treatements
  for(t in treatments){
    group <- predictions[predictions$Treatment == t , ]
    
    N <- nrow(group)
    
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
