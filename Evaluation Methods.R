####################################
# Evaluation metrics
####################################

# As described in Zhao et al. 2017

new_expected_outcome <- function(pred,response,control,treatment_list){
  N <- nrow(pred)
  for (t in c(treatment_list, control)) {
    assign(paste(t,"prob",sep="_"),sum(pred$Assignment == t)/N)
  }
  results <- 0
  for(counter in 1:N){
    if(pred[counter,]$Assignment==pred[counter,]$Treatment){
      results <- results + as.numeric(pred[counter,]$Outcome/
                                        eval(as.name(paste(pred[counter,]$Treatment,"prob",sep = "_"))))
    }
  }
  return(results/N)
}



## Modified Uplift Curve by Zhao
# *Assumption* outcome for each treatment equals prediction of model

new_expected_quantile_response <- function(response,control,treatment_list,predictions){
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  sorted_predictions <- predictions[order(-predictions$max_uplift),]
  sorted_predictions$Treatment <- as.character(sorted_predictions$Treatment)
  sorted_predictions$Assignment <- as.character(sorted_predictions$Assignment)
  sorted_predictions$temp_pred <- sorted_predictions$Treatment
  n_tenth <- round(nrow(predictions)/10)
  sorted_predictions$Treatment <- control
  deciles <- c(new_expected_outcome(sorted_predictions,response,control,treatment_list))
  for (x in 1:9) {
    sorted_predictions$Treatment <- sorted_predictions$temp_pred
    sorted_predictions$Treatment[(x*n_tenth):nrow(sorted_predictions)] <- control
    deciles <- c(deciles,new_expected_outcome(sorted_predictions,response,control,treatment_list))
  }
  sorted_predictions$Treatment <- sorted_predictions$temp_pred
  deciles <- c(deciles,new_expected_outcome(sorted_predictions,response,control,treatment_list))
}




#Percent Matched
perc_matched <- function(predictions){
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  sorted_predictions <- predictions[order(-predictions$max_uplift),]
  sorted_predictions$Treatment <- as.character(sorted_predictions$Treatment)
  n_tenth <- round(nrow(predictions)/10)
  deciles <- c()
  for (x in 1:9) {
    if(x == 1){
      new_data <- sorted_predictions[1:n_tenth,]
      deciles <- c(deciles,sum(new_data$Treatment == new_data$Assignment)/nrow(new_data))
    }else{
      new_data <- sorted_predictions[((x-1)*n_tenth):(x*n_tenth),]
      deciles <- c(deciles,sum(new_data$Treatment == new_data$Assignment)/nrow(new_data))
    }
  }
  new_data <- sorted_predictions[(9*n_tenth):nrow(sorted_predictions),]
  deciles <- c(deciles,sum(new_data$Treatment == new_data$Assignment)/nrow(new_data))
  return(deciles)
}








# Old Evaluation Methods ----
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


new_expected_outcomeV2 <- function(pred,response,control,treatment_list){
  N <- nrow(pred)
  for (t in c(treatment_list, control)) {
    assign(paste(t,"prob",sep="_"),sum(pred$Assignment == t)/N)
  }
  results <- 0
  n_matches <- 0
  for(counter in 1:N){
    if(pred[counter,]$Assignment==pred[counter,]$Treatment){
      results <- results + as.numeric(pred[counter,]$Outcome/
                                        eval(as.name(paste(pred[counter,]$Treatment,"prob",sep = "_"))))
      n_matches <- n_matches + 1
    }
  }
  return(results/n_matches)
}

expected_percentile_response <- function(predictions){
  # Choose only the uplift columns
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  
  predictions$max_treatment_outcome <- apply(predictions[ , c(1: (length(levels(as.factor(predictions$Assignment))) - 1)  )], 1 , max)
  
  # Sum percentiles
  ret <- data.frame(matrix(ncol = 2, nrow = 0))
  
  control_level <- if (nrow(predictions[predictions$Assignment == 'control' , ]) == 0) "Control" else "control"
  
  for(x in seq(0,1, 0.1)){
    # for top x set T to max T
    predictions$T_index <- apply(predictions[, 1:2], 1, which.max)
    # Assign optimal treatment for all
    predictions$Treatment <- colnames(predictions)[predictions$T_index]
    
    # quick fix
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


new_expected_quantile_responseV2 <- function(response,control,treatment_list,predictions){
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  sorted_predictions <- predictions[order(-predictions$max_uplift),]
  n_tenth <- round(nrow(predictions)/10)
  deciles <- c(0)
  for (x in 1:9) {
    deciles <- c(deciles,new_expected_outcome(sorted_predictions[1:(x*n_tenth),],response,control,treatment_list))
  }
  deciles <- c(deciles,new_expected_outcome(sorted_predictions,response,control,treatment_list))
}

# Incremental Uplift Curve
uplift_curve <- function(predictions, control_level){
  # score for each T individually
  N_total <- nrow(predictions)
  
  c_group <- predictions[predictions$Assignment == control_level , ]
  N_c <- nrow(c_group)
  
  ret <- data.frame(Percentile=  seq(0,1, 0.1))
  
  treatments <- levels(as.factor(predictions$Assignment[predictions$Assignment != control_level]))
  
  
  # Uplift Curve for each treatment
  for(t in treatments) {
    
    tmp <- predictions[predictions$Assignment == t, ]
    
    # score by uplift column of T
    tmp <- tmp[order(tmp[ , paste0("uplift_",t)] , decreasing = T) , ]
    c_tmp <- c_group[order(c_group[ , paste0("uplift_",t)] , decreasing = T) , ]
    
    N_t  <- nrow(tmp)
    
    outcomes <- c()
    # For each decile
    for(x in seq(0,1, 0.1)){
      
      # Formula from Gutierrez 2017
      weighted_uplift <- (mean(head(tmp$Outcome, x * N_t)) - mean(head(c_tmp$Outcome, x * N_c)) ) * (N_total * x) 
      
      outcomes <- c(outcomes, weighted_uplift)
    }     
    ret[ , t] <- outcomes
    
  }
  
  # Combined Uplift Curve
  tmp <- predictions[predictions$Assignment != control_level, ]
  
  # score by max T
  tmp$max_uplift <- apply(tmp[ , grep("^uplift",colnames(tmp))], 1 , max)
  c_group$max_uplift <- apply(c_group[ , grep("^uplift",colnames(c_group))], 1 , max)
  
  tmp <- tmp[order(tmp[ , "max_uplift"] , decreasing = T) , ]
  c_tmp <- c_group[order(c_group[ , "max_uplift"] , decreasing = T) , ]
  
  N_t  <- nrow(tmp)
  
  outcomes <- c()
  # For each decile
  for(x in seq(0,1, 0.1)){
    
    # Formula from Gutierrez 2017
    weighted_uplift <- (mean(head(tmp$Outcome, x * N_t)) - mean(head(c_tmp$Outcome, x * N_c)) ) * (N_total * x) 
    
    outcomes <- c(outcomes, weighted_uplift)
  }     
  ret$combined <- outcomes
  
  
  ret[is.na(ret)] <- 0 
  
  return(ret)
}





# incremental Qini - Curve
qini_curve <- function(predictions, control_level){
  # score for each T individually
  c_group <- predictions[predictions$Assignment == control_level , ]
  N_c <- nrow(c_group)
  
  ret <- data.frame(Percentile=  seq(0,1, 0.1))
  
  treatments <- levels(as.factor(predictions$Assignment[predictions$Assignment != control_level]))
  
  treatments <- treatments[treatments != control_level]
  
  for(t in treatments) {
    print(t)
    tmp <- predictions[predictions$Assignment == t, ]
    
    # score by uplift column of T
    tmp <- tmp[order(tmp[ , paste0("uplift_",t)] , decreasing = T) , ]
    c_tmp <- c_group[order(c_group[ , paste0("uplift_",t)] , decreasing = T) , ]
    
    N_t  <- nrow(tmp)
    
    outcomes <- c()
    # For each decile
    for(x in seq(0,1, 0.1)){
      # Radcliffe 2007 
      # u = R_t - ((R_c * N_t) / N_c)
      qini_gain <- sum(head(tmp$Outcome, x * N_t)) - ((sum(head(c_tmp$Outcome, x * N_c)) * N_t) / N_c)
      
      outcomes <- c(outcomes, qini_gain)
    }     
    ret[ , t] <- outcomes
  }
  ret[is.na(ret)] <- 0 
  
  return(ret)
}
