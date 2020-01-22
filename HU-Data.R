source("DecisionTreeImplementation.R")
source("RzepakowskiTree.R")
source("Evaluation Methods.R")
source("Causal Forest.R")
source("CausalTree.R")
library(caret)
set.seed(1234)
#Preprocessing---- 
hu_data <- read.csv("Data/explore_mt.csv",sep = ";")
# for(x in levels(hu_data$DeviceCategory)){
#   hu_data[x] <- ifelse(hu_data$DeviceCategory == x ,1,0)
# }
# hu_data$DeviceCategory <- NULL
for(x in levels(hu_data$multi_treat)){
  hu_data[x] <- ifelse(hu_data$multi_treat == x ,1,0)
}
tempfunction <- function(x){
  templist <- strsplit(x,",")
  new_string <- ""
  r <- 1
  for(s in templist[[1]]){
    if (r == 2) {
      new_string <- paste(new_string,s,sep = ".")
    } else{
      new_string <- paste(new_string,s,sep="")
    }
    r <- r+1
  }
  return(new_string)
}

hu_data$checkoutAmount <- as.numeric(lapply(as.character(hu_data$checkoutAmount), tempfunction))
hu_data$Number.of.seconds.between.last.and.previous.views <- 
  as.numeric(lapply(as.character(hu_data$Number.of.seconds.between.last.and.previous.views), tempfunction))
for (x in grep("^log.of",colnames(hu_data))) {
  hu_data[,x] <- as.numeric(lapply(as.character(hu_data[,x]), tempfunction))
}
for (x in colnames(hu_data[,16:155])) {
  if(is.na(mean(hu_data[[x]]))){
    print(x)
  }
}
# for (x in colnames(hu_data[,16:160])) {
#   print(x)
#   print(max(hu_data[[x]]) == 1 && min(hu_data[[x]]) == 0)
#   if(max(hu_data[[x]]) == 1 && min(hu_data[[x]]) == 0){
#     hu_data[,x] <- as.factor(hu_data[,x])
#   }
# }
# hu_data$ZipCode <- NULL

# Model Building----
response <- 'checkoutAmount'
control <- '0'


# Split into test and train data
idx <- createDataPartition(y = hu_data[ , response], p=0.2, list = FALSE)

train <- hu_data[-idx, ]

test <- hu_data[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.2, list = FALSE)

val <- train[p_idx,]
val_train <- train[-p_idx,]

treatment_list <- levels(hu_data$multi_treat)[2:7]
n_treatments <- length(treatment_list)
test_list <- set_up_tests(train[,colnames(train[,16:155])],TRUE,max_cases = 5)

#Single Tree
raw_tree <- build_tree(val_train,0,5,treatment_list,response,control,test_list)
pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)
tree_pred <-  predict.dt.as.df(pruned_tree, test)
colnames(tree_pred) <- c(treatment_list,control)
tree_pred[ , "Treatment"] <- colnames(tree_pred)[apply(tree_pred[, c(treatment_list,control)], 1, which.max)]
tree_pred[ , "Assignment"] <- colnames(test[, c(treatment_list,control)])[apply(test[, c(treatment_list,control)], 1, which.max)]
tree_pred[, "Outcome"] <- test[,response]
for (t in treatment_list) {
  tree_pred[,paste("uplift",t,sep = "_")] <- tree_pred[t] - tree_pred[control]
}

exp_outcome_simple <- new_expected_outcome(tree_pred,response,control,treatment_list) 
exp_inc_outcome_simple <- new_expected_quantile_response(response,control,treatment_list,tree_pred)
  
#Forest
forest <- parallel_build_forest(train,val,treatment_list,response,'0',n_trees = 4,n_features = 15, 
                                pruning = F,max_depth = 5)
forest_pred <- predict_forest_df(forest, test)
colnames(forest_pred) <- c(treatment_list,control)
forest_pred[ , "Treatment"] <- colnames(forest_pred)[apply(forest_pred[, c(treatment_list,control)], 1, which.max)]
forest_pred[ , "Assignment"] <- colnames(test[, c(treatment_list,control)])[apply(test[, c(treatment_list,control)], 1, which.max)]
forest_pred[, "Outcome"] <- test[,response]
for (t in treatment_list) {
  forest_pred[,paste("uplift",t,sep = "_")] <- forest_pred[t] - forest_pred[control]
}

exp_outcome_simple_forest <- new_expected_outcome(forest_pred,response,control,treatment_list) 
exp_inc_outcome_simple_forest <- new_expected_quantile_response(response,control,treatment_list,forest_pred)

#Random Forest
random_forest <- parallel_build_random_forest(train,val,treatment_list,response,'0',n_trees = 4,n_features = 15, 
                                pruning = F,max_depth = 5)
random_forest_pred <- predict_forest_df(random_forest, test)
colnames(random_forest_pred) <- c(treatment_list,control)
random_forest_pred[ , "Treatment"] <- colnames(random_forest_pred)[apply(random_forest_pred[, c(treatment_list,control)], 1, which.max)]
random_forest_pred[ , "Assignment"] <- colnames(test[, c(treatment_list,control)])[apply(test[, c(treatment_list,control)], 1, which.max)]
random_forest_pred[, "Outcome"] <- test[,response]
for (t in treatment_list) {
  random_forest_pred[,paste("uplift",t,sep = "_")] <- random_forest_pred[t] - random_forest_pred[control]
}

exp_outcome_random_forest <- new_expected_outcome(random_forest_pred,response,control,treatment_list) 
exp_inc_outcome_random_forest <- new_expected_quantile_response(response,control,treatment_list,random_forest_pred)


#Causal Tree
causal_tree_pred <- newcausalTreePredicitons(train,test,treatment_list,response,control)

#Causal Forest
causal_forest_pred <- newCausalForestPredicitons(train, test, treatment_list, response,control)
causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, treatment_list], 1, which.max)]
for (t in treatment_list) {
  causal_forest_pred[,paste("uplift",t,sep = "_")] <- causal_forest_pred[t] - causal_forest_pred[control]
}
exp_outcome_causal_forest <- new_expected_outcome(test,response,control,treatment_list,causal_forest_pred$Treatment)
exp_inc_outcome_causal_forest <- new_expected_quantile_response(test,response,control,treatment_list,causal_forest_pred)
