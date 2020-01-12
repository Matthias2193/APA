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
raw_tree <- build_tree(train,0,100,treatment_list,response,control,test_list)

pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)


### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]


write.csv(pred, 'Predictions/tree_Old.csv', row.names = FALSE)

pred <- read.csv('Predictions/tree_Old.csv')
exp_tree_old <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_tree_old <- new_expected_quantile_response(response,control,treatment_list,pred)



#Forest
forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 50,n_features = 3, pruning = F)

# add to the result df the outcome, assignment and calculate uplift for each T
start_time <- Sys.time()
pred <- parallel_predict_forest_df(forest, test)
time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
pred1 <- parallel_predict_forest_average_apply(forest,test)
time2 <- difftime(Sys.time(),start_time,units = "secs")
### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

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



start_time <- Sys.time() 
pred1 <- predict.dt.as.df(pruned_tree,test)
time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
pred2 <- predict.dt_apply(pruned_tree,test)
time2 <- difftime(Sys.time(),start_time,units = "secs")

sum(pred1==pred2)

predict.dt_apply <- function(tree,new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(new_data)
  temp_function <- function(x,node){
    type = 'root'
    while(type != 'leaf'){
      split = node[['split']]
      if(type_list[[names(split)]] == 'factor'){
        if(x[names(split)] == split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      } else{
        if(as.numeric(x[names(split)]) < split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      }
    }
    return(node[["results"]])
  }
  results <- data.frame(t(apply(new_data,1,temp_function,node=tree)))
  results <- split(results, seq(nrow(results)))
  return(results)
}

predict.dt.as.df_apply <- function(tree,new_data){
  type_list <- sapply(new_data, class)
  names(type_list) = colnames(new_data)
  temp_function <- function(x,node){
    type = 'root'
    while(type != 'leaf'){
      split = node[['split']]
      if(type_list[[names(split)]] == 'factor'){
        if(x[names(split)] == split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      } else{
        if(as.numeric(x[names(split)]) < split[[1]]){
          node = node[['left']]
          type = node[['type']]
        } else{
          node = node[['right']]
          type = node[['type']]
        }
      }
    }
    return(node[["results"]])
  }
  results <- data.frame(t(apply(new_data,1,temp_function,node=tree)))
  return(results)
}

parallel_predict_forest_average_apply <- function(forest,test_data){
  predictions <- list()
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores-1)
  registerDoParallel(cl)
  predictions <- foreach(x = 1:length(forest)) %dopar%{
    source('DecisionTreeImplementation.R')
    tree <- forest[[x]]
    new_data <- test_data
    type_list <- sapply(new_data, class)
    names(type_list) = colnames(new_data)
    temp_function <- function(x,node){
      type = 'root'
      while(type != 'leaf'){
        split = node[['split']]
        if(type_list[[names(split)]] == 'factor'){
          if(x[names(split)] == split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        } else{
          if(as.numeric(x[names(split)]) < split[[1]]){
            node = node[['left']]
            type = node[['type']]
          } else{
            node = node[['right']]
            type = node[['type']]
          }
        }
      }
      return(node[["results"]])
    }
    results <- data.frame(t(apply(new_data,1,temp_function,node=tree)))
    return(results)
  }
  stopCluster(cl)
  final_predictions <- predictions[[1]]
  for(x in 2:length(predictions)){
    final_predictions <- final_predictions+predictions[[x]]
  }
  final_predictions <- final_predictions/length(predictions)
  return(final_predictions)
}

for(x in 1:10){
  temp_df <- predict.dt.as.df_apply(forest[[x]],test)
  temp_df2 <- predict.dt.as.df(forest[[x]],test)
  if((sum(temp_df==predictions[[x]])==57600) && (sum(temp_df2==predictions[[x]])==57600)){
    print("Same")
  }
}