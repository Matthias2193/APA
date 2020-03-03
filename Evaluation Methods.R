####################################
# Evaluation metrics
####################################

# This function calculates the expected outcome as described in Zhao et al. 2017
expected_outcome <- function(pred,response,control,treatment_list){
  N <- nrow(pred)
  for (t in c(treatment_list, control)) {
    assign(paste(t,"prob",sep="_"),sum(pred$Assignment == t)/N)
  }
  temp_function <- function(x){
    if(x["Assignment"]==x["Treatment"]){
      return(as.numeric(x["Outcome"])/eval(as.name(paste(x["Treatment"],"prob",sep = "_"))))
    } else{
      return(0)
    }
  }
  results <- sum(apply(pred,1,temp_function))
  return(results/N)
}


# Modified Uplift Curve by Zhao
new_expected_quantile_response <- function(response,control,treatment_list,predictions){
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  sorted_predictions <- predictions[order(-predictions$max_uplift),]
  sorted_predictions$Treatment <- as.character(sorted_predictions$Treatment)
  sorted_predictions$Assignment <- as.character(sorted_predictions$Assignment)
  sorted_predictions$temp_pred <- sorted_predictions$Treatment
  n_tenth <- round(nrow(predictions)/10)
  sorted_predictions$Treatment <- control
  deciles <- c(expected_outcome(sorted_predictions,response,control,treatment_list))
  for (x in 1:9) {
    sorted_predictions$Treatment <- sorted_predictions$temp_pred
    sorted_predictions$Treatment[(x*n_tenth):nrow(sorted_predictions)] <- control
    deciles <- c(deciles,expected_outcome(sorted_predictions,response,control,treatment_list))
  }
  sorted_predictions$Treatment <- sorted_predictions$temp_pred
  deciles <- c(deciles,expected_outcome(sorted_predictions,response,control,treatment_list))
  return(deciles)
}


# Given a prediction, this function checks how many what percentage of people where assigned a treatment 
# by the model
perc_treated <- function(pred, treatment_list, n_pred = NULL){
  perc <- c()
  for(t in treatment_list){
    if(is.null(n_pred)){
      perc <- c(perc, sum(pred$Treatment == t)/nrow(pred))
    } else{
      perc <- c(perc, sum(pred$Treatment == t)/n_pred)
    }
    
  }
  return(perc)
}

decile_perc_treated <- function(pred, treatment_list){
  pred$max_uplift <- apply(pred[ , grep("^uplift",colnames(pred))], 1 , max)
  sorted_predictions <- pred[order(-pred$max_uplift),]
  perc <- c()
  n_tenth <- round(nrow(pred)/10)
  for(x in 1:9){
    perc <- rbind(perc, cbind(perc_treated(sorted_predictions[1:(x*n_tenth),], treatment_list, nrow(pred)),
                              treatment_list, rep(x,length(treatment_list))))
  }
  perc <- rbind(perc,cbind(perc_treated(sorted_predictions, treatment_list),treatment_list,
                           rep(10,length(treatment_list))))
  return(perc)
}

# Given a prediction, this function calculates the percentage of people assigned a treatment by the model, for
# each decile.
n_treated_decile <- function(pred,control){
  pred$max_uplift <- apply(pred[ , grep("^uplift",colnames(pred))], 1 , max)
  sorted_predictions <- pred[order(-pred$max_uplift),]
  n_tenth <- round(nrow(pred)/10)
  treated <- c(0)
  for(x in 1:9){
    temp_data <- sorted_predictions[1:(n_tenth*x),]
    treated <- c(treated,sum(temp_data$Treatment != control))
  }
  temp_data <- sorted_predictions
  treated <- c(treated,sum(temp_data$Treatment != control))
  return(treated)
}


# Old Evaluation Methods ----
# Incremental Uplift Curve
uplift_curve <- function(predictions, control_level,treatments){
  # score for each T individually
  N_total <- nrow(predictions)
  
  c_group <- predictions[predictions$Assignment == control_level , ]
  N_c <- nrow(c_group)
  
  ret <- data.frame(Percentile=  seq(0,1, 0.1))

  
  # Uplift Curve for each treatment
  # for(t in treatments) {
  #   
  #   tmp <- predictions[predictions$Assignment == t, ]
  #   
  #   # score by uplift column of T
  #   tmp <- tmp[order(tmp[ , paste0("uplift_",t)] , decreasing = T) , ]
  #   c_tmp <- c_group[order(c_group[ , paste0("uplift_",t)] , decreasing = T) , ]
  #   
  #   N_t  <- nrow(tmp)
  #   
  #   outcomes <- c()
  #   # For each decile
  #   for(x in seq(0,1, 0.1)){
  #     
  #     # Formula from Gutierrez 2017
  #     weighted_uplift <- (mean(head(tmp$Outcome, x * N_t)) - mean(head(c_tmp$Outcome, x * N_c)) ) * (N_total * x) 
  #     
  #     outcomes <- c(outcomes, weighted_uplift)
  #   }     
  #   ret[ , t] <- outcomes
  #   
  # }
  # 
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
qini_curve <- function(predictions, control_level,treatments){
  # score for each T individually
  predictions$max_uplift <- apply(predictions[ , grep("^uplift",colnames(predictions))], 1 , max)
  c_group <- predictions[predictions$Assignment == control_level , ]
  N_c <- nrow(c_group)
  
  
  
  ret <- data.frame(Percentile=  seq(0,1, 0.1))
  tmp <- predictions[predictions$Assignment != control_level, ]
  
  # score by uplift column of T
  tmp <- tmp[order(tmp[ , "max_uplift"] , decreasing = T) , ]
  c_tmp <- c_group[order(c_group[ ,"max_uplift"] , decreasing = T) , ]
  
  N_t  <- nrow(tmp)
  
  outcomes <- c()
  # For each decile
  for(x in seq(0,1, 0.1)){
    # Radcliffe 2007 
    # u = R_t - ((R_c * N_t) / N_c)
    qini_gain <- sum(head(tmp$Outcome, x * N_t)) - ((sum(head(c_tmp$Outcome, x * N_c)) * N_t) / N_c)
    
    outcomes <- c(outcomes, qini_gain)
  }     
  ret$Values <- outcomes
  result <- ret
  result[is.na(result)] <- 0 
  return(result)
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
