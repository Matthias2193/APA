source("DecisionTreeImplementation.R")
source("Evaluation Methods.R")
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
train <- train[-p_idx,]

treatment_list <- levels(hu_data$multi_treat)[2:7]
n_treatments <- length(treatment_list)
test_list <- set_up_tests(train[,colnames(train[,16:155])],TRUE,max_cases = 5)

#Single Tree
start_time <- Sys.time()
raw_tree <- build_tree(train,0,5,treatment_list,response,control,test_list,criterion = 2)
end_time <- Sys.time()
tree_time <- difftime(end_time,start_time)
print(tree_time)
pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)
tree_pred <-  predict.dt.as.df(pruned_tree, test)
tree_pred[ , "Treatment"] <- colnames(tree_pred)[apply(tree_pred[, treatment_list], 1, which.max)]
exp_outcome_simple_tree <- new_expected_outcome(test,response,control,treatment_list,tree_pred$Treatment) 
  
#Forest
start_time <- Sys.time()
forest <- parallel_build_forest(train,val,treatment_list,response,'0',n_trees = 50,n_features = 15,criterion = 2, 
                                pruning = F,l = rep(1/n_treatments,n_treatments),
                                g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments),max_depth = 5)
end_time <- Sys.time()
forest_time <- difftime(end_time,start_time)
print(forest_time)
forest_pred <- parallel_predict_forest_df(forest, test)
forest_pred[ , "Treatment"] <- colnames(forest_pred)[apply(forest_pred[, treatment_list], 1, which.max)]
exp_outcome_simple_forest <- new_expected_outcome(test,response,control,treatment_list,forest_pred$Treatment)
