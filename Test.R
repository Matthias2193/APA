#Test Script


library(ggplot2)
library(caret)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('RzepakowskiTree.R')
source('Evaluation Methods.R')


set.seed(2193)

#Data import
email <- read.csv('Data/Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'
control <- "control"

idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

test <- email[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE, max_cases = 10)



# Single Tree
start_time <- Sys.time()
raw_tree <- build_tree(train,0,100,treatment_list,response,control,test_list)
new_tree_time <- difftime(Sys.time(),start_time,units = "secs")
pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df_apply(pruned_tree, test)


### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)


write.csv(pred, 'Predictions/tree_Old.csv', row.names = FALSE)

pred <- read.csv('Predictions/tree_Old.csv')
exp_tree_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_tree_old <- new_expected_quantile_response(response,control,treatment_list,pred)



#Forest
forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3, pruning = F)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict_forest_df(forest,test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "New Assignment"] <- predictions_to_treatment(test, treatment_list, control)

write.csv(pred, 'Predictions/fores_Old.csv', row.names = FALSE)

pred <- read.csv('Predictions/fores_Old.csv')
exp_forest_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_forest_old <- new_expected_quantile_response(response,control,treatment_list,pred)



#Load original predictions
pred <- read.csv('Predictions/tree_Old.csv')
exp_tree_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_tree_old <- new_expected_quantile_response(response,control,treatment_list,pred)
pred <- read.csv('Predictions/fores_Old.csv')
exp_forest_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_forest_old <- new_expected_quantile_response(response,control,treatment_list,pred)




temp_df = data.frame(models=c("C_Forest","C_Tree","Simple_Tree","Simple_Forest","Naive Men", "Naive Women","SMA-DT","SMA-RF"),values = c(exp_outcome_c_forest,exp_outcome_c_tree,exp_outcome_simple,exp_outcome_simple_forest,exp_outcome_naive_men,exp_outcome_naive_women,exp_outcome_sma_dt,exp_outcome_sma_rf))

p<-ggplot(data=temp_df, aes(x=reorder(models, values), y=values)) +
  geom_bar(stat="identity") +
  labs(
    x = "Models",
    y = "Expected Outcome"
  )
p

temp_vec <- rep(seq(0,1,0.1),6)
name_vec <- c(rep("C_Forest",11),rep("C_Tree",11),rep("Simple_Tree",11),rep("Simple_Forest",11),rep("Naive Men",11),rep("Naive Women",11),rep("SMA-DT",11),rep("SMA-RF",11))
temp_results <- cbind(c(exp_inc_outcome_c_forest,exp_inc_outcome_c_tree,exp_inc_outcome_simple,exp_inc_outcome_simple_forest,exp_inc_outcome_naive_men,exp_inc_outcome_naive_women,exp_inc_outcome_sma_dt,exp_inc_outcome_sma_rf))
temp_df <- data.frame(cbind(temp_vec,temp_results,name_vec))
colnames(temp_df) <- c("Percentile","Expected_Outcome","Model")
temp_df$Expected_Outcome <- as.numeric(as.character(temp_df$Expected_Outcome))

p <-  melt(temp_df, id.vars = c("Percentile","Model")) %>% ggplot(aes(x = Percentile)) +
  geom_line(aes(y = value, group= Model, color= Model), size=0.5 ) +
  labs(
    color="Base Learner",
    title = "Model Comparison",
    y = "Avg. Spending per Person",
    x ="Amount of Treated"
  ) +
  scale_colour_brewer(palette = "Dark2") +
  theme_light()
p





### Test Gain Methods

#This method calculates the gain for a given split
new_simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
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
  current_gain <- 0
  for(t in treatment){
    if(mean(data[data[,t]==1,target])>current_gain){
      current_gain <- mean(data[data[,t]==1,target])
    }
  }
  
  #Here the gain is calculated
  left_gain <- 0
  right_gain <- 0
  for(t in treatment){
    left_gain <- max(left_gain,mean(data1[data1[,t]==1,target]))
    right_gain <- max(right_gain,mean(data2[data2[,t]==1,target]))
  }
  gain <- max(left_gain,right_gain)
  # gain <- (frac1*left_gain+frac2*right_gain)
  for(t in treatments){
    if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
      gain <- -1
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



new_simple_gain <- function(test_case, treatment, control, target, data, test_type, test_col){
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
    s <- treatments[x+1]
    temp_gain <- (mean(data[data[,t] == 1,target])-mean(data[data[,s] == 1,target]))^2
    current_gain <- current_gain + temp_gain
  }
  #The actual calculation of the gain
  #Here for a test of a categorical cavariate
  for(x in 1:(length(treatments)-1)){
    t <- treatments[x]
    s <- treatments[x+1]
    temp_gain <- frac1*(mean(data1[data1[,t] == 1,target])-mean(data1[data1[,s] == 1,target]))^2 +
      frac2*(mean(data2[data2[,t] == 1,target])-mean(data2[data2[,s] == 1,target]))^2
    gain <- gain + temp_gain
  }
  #Make sure that there are data points of each treatment in each subset of the data
  # for(t in treatments){
  #   if(nrow(data1[data1[,t]==1,])==0 || nrow(data2[data2[,t]==1,]) == 0){
  #     gain <- 0
  #   }
  # }
  if(is.na(gain)){
    gain = -1
  }
  if(gain <= current_gain){
    gain = -1
  }
  return(gain)
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








pred["predicted_treatment"] <- predictions_to_treatment(pred, treatment_list, control)


















#Benchmark old prediction vs new prediction

start_time <- Sys.time()
forest_pred <- parallel_predict_forest_df(forest, test)
forest_time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
forest_pred1 <- parallel_predict_forest_average_apply(forest,test)
forest_time2 <- difftime(Sys.time(),start_time,units = "secs")

if(sum(forest_pred==forest_pred1)==nrow(test)*3){
  print("New and old predictions are identical.")
  print(paste(c("The old time was:",forest_time1,"and the new time was:",forest_time2),collapse = " "))
  print(paste(c("The difference was:", as.double(forest_time1)-as.double(forest_time2),"seconds"),collapse = " "))
  paste(c("The new method takes ",round(as.double(forest_time2)/as.double(forest_time1),4)*100,
          "% of the time."),collapse = "")
}


start_time <- Sys.time() 
tree_pred1 <- predict.dt.as.df(pruned_tree,test)
tree_time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
tree_pred2 <- predict.dt.as.df_apply(pruned_tree,test)
tree_time2 <- difftime(Sys.time(),start_time,units = "secs")

if(sum(tree_pred1==tree_pred2)==nrow(test)*3){
  print("New and old predictions are identical.")
  print(paste(c("The old time was:",tree_time1,"and the new time was:",tree_time2),collapse = " "))
  print(paste(c("The difference was:", as.double(tree_time1)-as.double(tree_time2),"seconds"),collapse = " "))
  paste(c("The new method takes ",round(as.double(tree_time2)/as.double(tree_time1),4)*100,
          "% of the time."),collapse = "")
}