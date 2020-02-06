
#####################################################################################
### Hillstorm Data
### https://blog.minethatdata.com/2008/03/minethatdata-e-mail-analytics-and-data.html
#####################################################################################

source('DecisionTreeImplementation.R')
source('Separate Model Approach.R')
source('Causal Forest.R')
source('Evaluation Methods.R')

set.seed(123)

#####################################################################################
### Conversion Prediction
#####################################################################################

#Data import
raw_email <- read.csv('Data/Email.csv')
response <- 'conversion'
control <- 'control'

#raw_email <-raw_email[1:10000 , ]

k <- 5
folds <- createFolds(raw_email[, response], k = k)
raw_email$fold <- 0

for(f in names(folds)){
  raw_email[folds[f][[1]] , "fold"] <- f
}

email <- raw_email

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit<- email$spend <- email$segment <- NULL

#################################

########
# Fit Rzp Tree
rzp_tree_exp_conv <- rzp_tree_up_conv <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list)
  
  pruned_tree <- prune_tree(raw_tree, val, treatment_list, test_list, response, 'control')
  
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
  
  # bind  matching and expected outcome evaluation
  rzp_tree_exp_conv <- rbind(rzp_tree_exp_conv, expected_percentile_response(pred))
  rzp_tree_up_conv <- rbind(rzp_tree_up_conv, uplift_curve(pred, "control"))
  
}
#aggregate by percentile r0esults of all models with mean 
rzp_tree_exp_conv <- aggregate(.~Percentile, rzp_tree_exp_conv, mean)
rzp_tree_up_conv <- aggregate(.~Percentile, rzp_tree_up_conv, mean)

write.csv(rzp_tree_exp_conv, 'CSV/rzp_tree_exp_conv.csv', row.names = FALSE)
write.csv(rzp_tree_up_conv, 'CSV/rzp_tree_up_conv.csv', row.names = FALSE)


######
# Tree simple Criterion
rzp_tree_exp_conv_simple <- rzp_tree_up_conv_simple <- NULL
for(f in names(folds)){
    train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train,0,100,treatment_list,response,'control',test_list,criterion = 2)
  
  pruned_tree <- raw_tree
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
  
  rzp_tree_exp_conv_simple <- rbind(rzp_tree_exp_conv_simple, expected_percentile_response(pred))
  rzp_tree_up_conv_simple <- rbind(rzp_tree_up_conv_simple, uplift_curve(pred, "control"))
  
}
rzp_tree_exp_conv_simple <- aggregate(.~Percentile, rzp_tree_exp_conv_simple, mean)
rzp_tree_up_conv_simple <- aggregate(.~Percentile, rzp_tree_up_conv_simple, mean)

write.csv(rzp_tree_exp_conv_simple, 'CSV/simple_tree_exp_conv.csv', row.names = FALSE)
write.csv(rzp_tree_up_conv_simple, 'CSV/simple_tree_up_conv.csv', row.names = FALSE)



## ## ## ## ## ## 
#Rzp Forest

rzp_forest_exp_conv <- rzp_forest_up_conv <- NULL
for(f in names(folds)){
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  rzp_forest <- build_forest(train_tree,val,treatment_list,response,control,n_trees = 20,n_features = 2,criterion = 1,
                             pruning = F,l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
  
  # add to the result df the outcome, assignment and calculate uplift for each T
  pred <- predict_forest_df(rzp_forest, test)
  
  ### Results Preparation to bring into equal format
  # Calculate Uplift for each T
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  rzp_forest_exp_conv <- rbind(rzp_forest_exp_conv, expected_percentile_response(pred))
  rzp_forest_up_conv <- rbind(rzp_forest_up_conv, uplift_curve(pred, "control"))
  
}
rzp_forest_exp_conv <- aggregate(.~Percentile, rzp_forest_exp_conv, mean)
rzp_forest_up_conv <- aggregate(.~Percentile, rzp_forest_up_conv, mean)

write.csv(rzp_forest_exp_conv, 'CSV/rzp_forest_exp_conv.csv', row.names = FALSE)
write.csv(rzp_forest_up_conv, 'CSV/rzp_forest_up_conv.csv', row.names = FALSE)


## ## ## ## ## ## 
#Simple Forest
simple_forest_exp_conv <- simple_forest_up_conv <- NULL
for(f in names(folds)){
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  simple_forest <- build_forest(train_tree,val,treatment_list,response,control,n_trees = 20,n_features = 2,criterion = 2,
                                pruning = F,l = rep(1/n_treatments,n_treatments),
                                g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
  
  # add to the result df the outcome, assignment and calculate uplift for each T
  pred <- predict_forest_df(simple_forest, test)
  
  ### Results Preparation to bring into equal format
  # Calculate Uplift for each T
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  simple_forest_exp_conv <- rbind(simple_forest_exp_conv, expected_percentile_response(pred))
  simple_forest_up_conv <- rbind(simple_forest_up_conv, uplift_curve(pred, "control"))
  
}
simple_forest_exp_conv <- aggregate(.~Percentile, simple_forest_exp_conv, mean)
simple_forest_up_conv <- aggregate(.~Percentile, simple_forest_up_conv, mean)

write.csv(simple_forest_exp_conv, 'CSV/simple_forest_exp_conv.csv', row.names = FALSE)
write.csv(simple_forest_up_conv, 'CSV/simple_forest_up_conv.csv', row.names = FALSE)



#######
# Causal Forest
c_forest_exp_conv <- c_forest_up_conv <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  treatment_list <- c('men_treatment','women_treatment')
  
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
  # Results Prep
  causal_forest_pred[ , "uplift_men_treatment"] <- causal_forest_pred[ , 1] - causal_forest_pred[ , 3]
  #
  causal_forest_pred[ , "uplift_women_treatment"] <- causal_forest_pred[ , 2] - causal_forest_pred[ , 3]
  causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, 1:3], 1, which.max)]
  
  causal_forest_pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  causal_forest_pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  #rbind  matching and expected outcome evaluation
  c_forest_exp_conv <- rbind(c_forest_exp_conv, expected_percentile_response(causal_forest_pred))
  c_forest_up_conv <- rbind(c_forest_up_conv, uplift_curve(causal_forest_pred, "control"))
}
c_forest_exp_conv <- aggregate(.~Percentile, c_forest_exp_conv, mean)
c_forest_up_conv <- aggregate(.~Percentile, c_forest_up_conv, mean)

write.csv(c_forest_exp_conv, 'CSV/c_forest_exp_conv.csv', row.names = FALSE)
write.csv(c_forest_up_conv, 'CSV/c_forest_up_conv.csv', row.names = FALSE)



####################################
####################################
# Import data and creat Train / Test splits
data <-  raw_email

data$visit <- NULL      
data$spend <- NULL
#data$conversion <- NULL      

treatment <- "segment"
response <-  "conversion"
control_level <- "No E-Mail"

treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]


#################
## Random Forest
sma_rf_conv_exp <- sma_rf_conv_up <- NULL
for(f in names(folds)){
  train_split <- data[data$fold != f , ]

  c_f <- c("fold")
  
  train_data <- list()
  for(x in treatments){
    t_data <- train_split[train_split[, treatment] == x, ]
    
    # Remove the treatment column from training data
    t_data <- t_data[ , names(t_data) != treatment]
    
    # Each Training Set as one element in list
    train_data <- append(train_data, list(x = t_data[t_data$fold != f , !names(t_data) %in% c_f]))
  }
  # Assign Treatment Names to Train Sets
  names(train_data) <- treatments
  # Add control group
  t_data <- train_split[train_split[, treatment] == control_level, ]
  t_data <- t_data[ , names(t_data) != treatment]
  train_data <- append(train_data, list(Control = t_data[t_data$fold != f, !names(t_data) %in% c_f]))
  
  # Test Data
  test_data <- data[data$fold == f , ]
  test_data$fold <- NULL
  
  rf_models_list <- rf_models(train_data, response, "class")
  
  rf_pred_class <- tree_sma_uplift(rf_models_list, test_data, response, treatment, control_level, "class")
  
  colnames(rf_pred_class)[4] <- "uplift_men_treatment"
  colnames(rf_pred_class)[5] <- "uplift_women_treatment"
  
  colnames(rf_pred_class)[1] <- "men_treatment"
  colnames(rf_pred_class)[2] <- "women_treatment"
  
  
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Mens E-Mail" , "Assignment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Womens E-Mail" , "Assignment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Control" , "Assignment"] <- "control"
  
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Mens E-Mail" , "Treatment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Womens E-Mail" , "Treatment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Control" , "Treatment"] <- "control"
  
  sma_rf_conv_exp <- rbind(sma_rf_conv_exp, expected_percentile_response(rf_pred_class))
  sma_rf_conv_up <- rbind(sma_rf_conv_up, uplift_curve(rf_pred_class, "control"))
}

write.csv(rf_pred_class,"rf_conv_pred.csv", row.names = F)

sma_rf_conv_exp <- aggregate(.~Percentile, sma_rf_conv_exp, mean)
sma_rf_conv_up <- aggregate(.~Percentile, sma_rf_conv_up, mean)

write.csv(sma_rf_conv_exp, 'CSV/sma_rf_exp_conv.csv', row.names = FALSE)
write.csv(sma_rf_conv_up, 'CSV/sma_rf_up_conv.csv', row.names = FALSE)

#####################################################################################
### Visit Prediction
#####################################################################################

#Data import
raw_email <- read.csv('Email.csv')
response <- 'visit'

k <- 5
folds <- createFolds(raw_email[, response], k = k)
raw_email$fold <- 0

for(f in names(folds)){
  raw_email[folds[f][[1]] , "fold"] <- f
}

email <- raw_email

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$conversion<- email$spend <- email$segment <- NULL

########
# Fit Rzp Tree
rzp_tree_exp_vis <- rzp_tree_up_vis <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list)
  
  pruned_tree <- pruned_tree <- prune_tree(raw_tree, val, treatment_list, test_list, response, 'control')
  
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
  
  # bind  matching and expected outcome evaluation
  rzp_tree_exp_vis <- rbind(rzp_tree_exp_vis, expected_percentile_response(pred))
  rzp_tree_up_vis <- rbind(rzp_tree_up_vis, uplift_curve(pred, "control"))
  
}
#aggregate by percentile r0esults of all models with mean 
rzp_tree_exp_vis <- aggregate(.~Percentile, rzp_tree_exp_vis, mean)
rzp_tree_up_vis <- aggregate(.~Percentile, rzp_tree_up_vis, mean)

write.csv(rzp_tree_exp_vis, 'CSV/rzp_tree_exp_vis.csv', row.names = FALSE)
write.csv(rzp_tree_up_vis, 'CSV/rzp_tree_up_vis.csv', row.names = FALSE)


######
# Tree simple Criterion
rzp_tree_exp_vis_simple <- rzp_tree_up_vis_simple <- NULL
for(f in names(folds)){
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train,0,100,treatment_list,response,'control',test_list,criterion = 2)
  
  pruned_tree <- raw_tree
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
  
  rzp_tree_exp_vis_simple <- rbind(rzp_tree_exp_vis_simple, expected_percentile_response(pred))
  rzp_tree_up_vis_simple <- rbind(rzp_tree_up_vis_simple, uplift_curve(pred, "control"))
  
}
rzp_tree_exp_vis_simple <- aggregate(.~Percentile, rzp_tree_exp_vis_simple, mean)
rzp_tree_up_vis_simple <- aggregate(.~Percentile, rzp_tree_up_vis_simple, mean)

write.csv(rzp_tree_exp_vis_simple, 'CSV/simple_tree_exp_vis.csv', row.names = FALSE)
write.csv(rzp_tree_up_vis_simple, 'CSV/simple_tree_up_vis.csv', row.names = FALSE)


#######
# Causal Forest
c_forest_exp_vis <- c_forest_up_vis <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  treatment_list <- c('men_treatment','women_treatment')
  
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
  # Results Prep
  causal_forest_pred[ , "uplift_men_treatment"] <- causal_forest_pred[ , 1] - causal_forest_pred[ , 3]
  #
  causal_forest_pred[ , "uplift_women_treatment"] <- causal_forest_pred[ , 2] - causal_forest_pred[ , 3]
  causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, 1:3], 1, which.max)]
  
  #levels(as.factor(causal_forest_pred[ , "Treatment"]))
  
  causal_forest_pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  causal_forest_pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  #rbind  matching and expected outcome evaluation
  c_forest_exp_vis <- rbind(c_forest_exp_vis, expected_percentile_response(causal_forest_pred))
  c_forest_up_vis <- rbind(c_forest_up_vis, uplift_curve(causal_forest_pred, "control"))
}
c_forest_exp_vis <- aggregate(.~Percentile, c_forest_exp_vis, mean)
c_forest_up_vis <- aggregate(.~Percentile, c_forest_up_vis, mean)

write.csv(c_forest_exp_vis, 'CSV/c_forest_exp_vis.csv', row.names = FALSE)
write.csv(c_forest_up_vis, 'CSV/c_forest_up_vis.csv', row.names = FALSE)


####################################
# SMA

# Import data and creat Train / Test splits
data <-  raw_email

data$spend <- NULL
data$conversion <- NULL      

treatment <- "segment"
response <-  "visit"
control_level <- "No E-Mail"

treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]


#################
## Random Forest
sma_rf_vis_exp <- sma_rf_vis_up <- NULL
for(f in names(folds)){
  train_split <- data[data$fold != f , ]
  
  c_f <- c("fold")
  
  train_data <- list()
  for(x in treatments){
    t_data <- train_split[train_split[, treatment] == x, ]
    
    # Remove the treatment column from training data
    t_data <- t_data[ , names(t_data) != treatment]
    
    # Each Training Set as one element in list
    train_data <- append(train_data, list(x = t_data[t_data$fold != f , !names(t_data) %in% c_f]))
  }
  # Assign Treatment Names to Train Sets
  names(train_data) <- treatments
  # Add control group
  t_data <- train_split[train_split[, treatment] == control_level, ]
  t_data <- t_data[ , names(t_data) != treatment]
  train_data <- append(train_data, list(Control = t_data[t_data$fold != f, !names(t_data) %in% c_f]))
  
  # Test Data
  test_data <- data[data$fold == f , ]
  test_data$fold <- NULL
  
  rf_models_list <- rf_models(train_data, response, "class")
  
  rf_pred_class <- tree_sma_uplift(rf_models_list, test_data, response, treatment, control_level, "class")
  
  colnames(rf_pred_class)[4] <- "uplift_men_treatment"
  colnames(rf_pred_class)[5] <- "uplift_women_treatment"
  
  colnames(rf_pred_class)[1] <- "men_treatment"
  colnames(rf_pred_class)[2] <- "women_treatment"
  
  
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Mens E-Mail" , "Assignment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Womens E-Mail" , "Assignment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Control" , "Assignment"] <- "control"
  
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Mens E-Mail" , "Treatment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Womens E-Mail" , "Treatment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Control" , "Treatment"] <- "control"
  
  sma_rf_vis_exp <- rbind(sma_rf_vis_exp, expected_percentile_response(rf_pred_class))
  sma_rf_vis_up <- rbind(sma_rf_vis_up, uplift_curve(rf_pred_class, "control"))
}

write.csv(rf_pred_class,"rf_vis_pred.csv", row.names = F)

sma_rf_vis_exp <- aggregate(.~Percentile, sma_rf_vis_exp, mean)
sma_rf_vis_up <- aggregate(.~Percentile, sma_rf_vis_up, mean)

write.csv(sma_rf_vis_exp, 'CSV/sma_rf_exp_vis.csv', row.names = FALSE)
write.csv(sma_rf_vis_up, 'CSV/sma_rf_up_vis.csv', row.names = FALSE)



#####################################################################################
### Spend Prediction
#####################################################################################
email <- raw_email

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'
control <- 'control'

########
# Fit Rzp Tree
rzp_tree_exp_spend <- rzp_tree_up_spend <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list, divergence = "EucDistance", criterion = 2)
  pruned_tree <- pruned_tree <- prune_tree(raw_tree, val, treatment_list, test_list, response, 'control')
  
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
  
  # bind  matching and expected outcome evaluation
  rzp_tree_exp_spend <- rbind(rzp_tree_exp_spend, expected_percentile_response(pred))
  rzp_tree_up_spend <- rbind(rzp_tree_up_spend, uplift_curve(pred, "control"))
  
}
#aggregate by percentile results of all models with mean 
rzp_tree_exp_spend <- aggregate(.~Percentile, rzp_tree_exp_spend, mean)
rzp_tree_up_spend <- aggregate(.~Percentile, rzp_tree_up_spend, mean)

write.csv(rzp_tree_exp_spend, 'CSV/rzp_tree_exp_spend.csv', row.names = FALSE)
write.csv(rzp_tree_up_spend, 'CSV/rzp_tree_up_spend.csv', row.names = FALSE)




######
# Tree simple Criterion
rzp_tree_exp_spend_simple <- rzp_tree_up_spend_simple <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list, divergence = "EucDistance", criterion = 2)
  
  pruned_tree <- raw_tree
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
  
  rzp_tree_exp_spend_simple <- rbind(rzp_tree_exp_spend_simple, expected_percentile_response(pred))
  rzp_tree_up_spend_simple <- rbind(rzp_tree_up_spend_simple, uplift_curve(pred, "control"))
  
}
rzp_tree_exp_spend_simple <- aggregate(.~Percentile, rzp_tree_exp_spend_simple, mean)
rzp_tree_up_spend_simple <- aggregate(.~Percentile, rzp_tree_up_spend_simple, mean)

write.csv(rzp_tree_exp_spend_simple, 'CSV/simple_tree_exp_spend.csv', row.names = FALSE)
write.csv(rzp_tree_up_spend_simple, 'CSV/simple_tree_up_spend.csv', row.names = FALSE)


######
# Rzp Forest
rzp_forest_exp_spend <- rzp_forest_up_spend <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  rzp_forest <- build_forest(train_tree,val,treatment_list,response,control,n_trees = 20,n_features = 2,criterion = 1,
                             pruning = F,l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
  
  # add to the result df the outcome, assignment and calculate uplift for each T
  pred <- predict_forest_df(rzp_forest, test)
  
  ### Results Preparation to bring into equal format
  # Calculate Uplift for each T
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  rzp_forest_exp_spend <- rbind(rzp_forest_exp_spend, expected_percentile_response(pred))
  rzp_forest_up_spend <- rbind(rzp_forest_up_spend, uplift_curve(pred, "control"))
  
}
rzp_forest_exp_spend <- aggregate(.~Percentile, rzp_forest_exp_spend, mean)
rzp_forest_up_spend <- aggregate(.~Percentile, rzp_forest_up_spend, mean)

write.csv(rzp_forest_exp_spend, 'CSV/rzp_forest_exp_spend.csv', row.names = FALSE)
write.csv(rzp_forest_up_spend, 'CSV/rzp_forest_up_spend.csv', row.names = FALSE)



######
# Simple Forest
simple_forest_exp_spend <- simple_forest_up_spend <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  val <- train[p_idx,]
  train_tree <- train[-p_idx,]
  
  treatment_list <- c('men_treatment','women_treatment')
  test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                          "newbie","channel")],TRUE)
  
  rzp_forest <- build_forest(train_tree,val,treatment_list,response,control,n_trees = 20,n_features = 2,criterion = 1,
                             pruning = F,l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
  
  # add to the result df the outcome, assignment and calculate uplift for each T
  pred <- predict_forest_df(rzp_forest, test)
  
  ### Results Preparation to bring into equal format
  # Calculate Uplift for each T
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  simple_forest_exp_spend <- rbind(simple_forest_exp_spend, expected_percentile_response(pred))
  simple_forest_up_spend <- rbind(simple_forest_up_spend, uplift_curve(pred, "control"))
  
}
simple_forest_exp_spend <- aggregate(.~Percentile, simple_forest_exp_spend, mean)
simple_forest_up_spend <- aggregate(.~Percentile, simple_forest_up_spend, mean)

write.csv(simple_forest_exp_spend, 'CSV/simple_forest_exp_spend.csv', row.names = FALSE)
write.csv(simple_forest_up_spend, 'CSV/simple_forest_up_spend.csv', row.names = FALSE)

#######
# Causal Forest
c_forest_exp_spend <- c_forest_up_spend <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
  # Results Prep
  causal_forest_pred[ , "uplift_men_treatment"] <- causal_forest_pred[ , 1] - causal_forest_pred[ , 3]
  causal_forest_pred[ , "uplift_women_treatment"] <- causal_forest_pred[ , 2] - causal_forest_pred[ , 3]
  causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, 1:3], 1, which.max)]
  
  causal_forest_pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  causal_forest_pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  #rbind  matching and expected outcome evaluation
  c_forest_exp_spend <- rbind(c_forest_exp_spend, expected_percentile_response(causal_forest_pred))
  c_forest_up_spend <- rbind(c_forest_up_spend, uplift_curve(causal_forest_pred, "control"))
}
c_forest_exp_spend <- aggregate(.~Percentile, c_forest_exp_spend, mean)
c_forest_up_spend <- aggregate(.~Percentile, c_forest_up_spend, mean)

write.csv(c_forest_exp_spend, 'CSV/c_forest_exp_spend.csv', row.names = FALSE)
write.csv(c_forest_up_spend, 'CSV/c_forest_up_spend.csv', row.names = FALSE)


####################################
####################################
# Import data and creat Train / Test splits
data <-  raw_email

data$visit <- NULL      
data$conversion <- NULL

treatment <- "segment"
response <-  "spend"
control_level <- "No E-Mail"

treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]


#################
## Random Forest
sma_rf_exp_spend <- sma_rf_up_spend <- NULL
for(f in names(folds)){
  train_split <- data[data$fold != f , ]
  
  train_data <- list()
  # make train_data list for each T
  for(x in treatments){
    t_data <- train_split[train_split[, treatment] == x, ]
    
    c_f <- c("fold")
    
    # Remove the treatment column from training data
    t_data <- t_data[ , names(t_data) != treatment]
    
    # Each Training Set as one element in list
    train_data <- append(train_data, list(x = t_data[t_data$fold != f, !names(t_data) %in% c_f]))
  }
  # Assign Treatment Names to Train Sets
  names(train_data) <- treatments
  # Add control group
  t_data <- train_split[train_split[, treatment] == control_level, ]
  t_data <- t_data[ , names(t_data) != treatment]
  train_data <- append(train_data, list(Control = t_data[t_data$fold != f, !names(t_data) %in% c_f]))
  
  # Test Data
  test_data <- data[data$fold == f , ]
  test_data$fold <- NULL
  
  rf_models_list <- rf_models(train_data, response, "anova")
  
  rf_pred_class <- tree_sma_uplift(rf_models_list, test_data, response, treatment, control_level, "anova")
  
  
  colnames(rf_pred_class)[4] <- "uplift_men_treatment"
  colnames(rf_pred_class)[5] <- "uplift_women_treatment"
  
  colnames(rf_pred_class)[1] <- "men_treatment"
  colnames(rf_pred_class)[2] <- "women_treatment"
  
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Mens E-Mail" , "Assignment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Womens E-Mail" , "Assignment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Assignment"] == "Control" , "Assignment"] <- "control"
  
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Mens E-Mail" , "Treatment"] <- "men_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Womens E-Mail" , "Treatment"] <- "women_treatment"
  rf_pred_class[rf_pred_class[ , "Treatment"] == "Control" , "Treatment"] <- "control"
  
  sma_rf_exp_spend <- rbind(sma_rf_exp_spend, expected_percentile_response(rf_pred_class))
  sma_rf_up_spend <- rbind(sma_rf_up_spend, uplift_curve(rf_pred_class, "control"))
}

sma_rf_exp_spend <- aggregate(.~Percentile, sma_rf_exp_spend, mean)
sma_rf_up_spend <- aggregate(.~Percentile, sma_rf_up_spend, mean)

write.csv(sma_rf_exp_spend, 'CSV/sma_rf_exp_spend.csv', row.names = FALSE)
write.csv(sma_rf_up_spend, 'CSV/sma_rf_up_spend.csv', row.names = FALSE)
