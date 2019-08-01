
source('DecisionTreeImplementation.R')
source('X Model Approach.R')
source('Causal Forest.R')
source('Evaluation Methods.R')

set.seed(123)

#####################################################################################
### Conversion Prediction
#####################################################################################

#Data import
raw_email <- read.csv('Email.csv')
response <- 'conversion'

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

#email$conversion 
email$visit<- email$spend <- email$segment <- NULL

########
# Fit Rzp Tree
rzp_tree_exp_conv <- rzp_tree_mat_conv <- NULL
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
  
  pruned_tree <- prune_tree(raw_tree,val,train_tree,target = response)
  #pruned_tree <- raw_tree
  
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
  rzp_tree_mat_conv <- rbind(rzp_tree_mat_conv, matching_evaluation(pred, "control"))
  
}
#aggregate by percentile r0esults of all models with mean 
rzp_tree_exp_conv <- aggregate(.~Percentile, rzp_tree_exp_conv, mean)
rzp_tree_mat_conv <- aggregate(.~Percentile, rzp_tree_mat_conv, mean)

write.csv(rzp_tree_exp_conv, 'CSV-Conv/rzp_tree_exp_conv.csv', row.names = FALSE)
write.csv(rzp_tree_mat_conv, 'CSV-Conv/rzp_tree_mat_conv.csv', row.names = FALSE)


######
# Tree simple Criterion
rzp_tree_exp_conv_simple <- rzp_tree_mat_conv_simple <- NULL
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
  rzp_tree_mat_conv_simple <- rbind(rzp_tree_mat_conv_simple, matching_evaluation(pred, "control"))
  
}
rzp_tree_exp_conv_simple <- aggregate(.~Percentile, rzp_tree_exp_conv_simple, mean)
rzp_tree_mat_conv_simple <- aggregate(.~Percentile, rzp_tree_mat_conv_simple, mean)

write.csv(rzp_tree_exp_conv_simple, 'CSV-Conv/rzp_tree_exp_conv_simple.csv', row.names = FALSE)
write.csv(rzp_tree_mat_conv_simple, 'CSV-Conv/rzp_tree_mat_conv_simple.csv', row.names = FALSE)


#######
# Causal Forest
c_forest_exp_conv <- c_forest_mat_conv <- NULL
for(f in names(folds)){
  # split into train and test
  train <- email[email$fold != f , ]
  test <- email[email$fold == f , ]
  train$fold <- test$fold <- NULL
  
  
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
  c_forest_exp_conv <- rbind(c_forest_exp_conv, expected_percentile_response(causal_forest_pred))
  c_forest_mat_conv <- rbind(c_forest_mat_conv, matching_evaluation(causal_forest_pred, "control"))
}
c_forest_exp_conv <- aggregate(.~Percentile, c_forest_exp_conv, mean)
c_forest_mat_conv <- aggregate(.~Percentile, c_forest_mat_conv, mean)

write.csv(c_forest_exp_conv, 'CSV-Conv/c_forest_exp_conv.csv', row.names = FALSE)
write.csv(c_forest_mat_conv, 'CSV-Conv/c_forest_mat_conv.csv', row.names = FALSE)

###
mean(causal_forest_pred$Outcome[causal_forest_pred$Assignment == 'control' ])

levels(as.factor(causal_forest_pred$Treatment))

min(causal_forest_pred$uplift_men_treatment)
min(causal_forest_pred$uplift_women_treatment)
## -> never negative uplift predicted....
###

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
sma_rf_conv_exp <- sma_rf_conv_mat <- NULL
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
  
  rf_pred_class <- dt_x_model_predictions(rf_models_list, test_data, response, treatment, control_level, "class")
  
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
  
  #predictions <- rf_pred_class
  
  ##
  sma_rf_conv_exp <- rbind(sma_rf_conv_exp, expected_percentile_response(rf_pred_class))
  sma_rf_conv_mat <- rbind(sma_rf_conv_mat, matching_evaluation(rf_pred_class, "control"))
}

write.csv(rf_pred_class,"rf_conv_pred.csv", row.names = F)

sma_rf_conv_exp <- aggregate(.~Percentile, sma_rf_conv_exp, mean)
sma_rf_conv_mat <- aggregate(.~Percentile, sma_rf_conv_mat, mean)

write.csv(sma_rf_conv_exp, 'CSV-Conv/sma_rf_conv_exp.csv', row.names = FALSE)
write.csv(sma_rf_conv_mat, 'CSV-Conv/sma_rf_conv_mat.csv', row.names = FALSE)


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


########
# Fit Rzp Tree
rzp_tree_exp_spend <- rzp_tree_mat_spend <- NULL
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
  pruned_tree <- prune_tree(raw_tree,val,train_tree,target = response)
  
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
  rzp_tree_mat_spend <- rbind(rzp_tree_mat_spend, matching_evaluation(pred, "control"))
  
}
#aggregate by percentile results of all models with mean 
rzp_tree_exp_spend <- aggregate(.~Percentile, rzp_tree_exp_spend, mean)
rzp_tree_mat_spend <- aggregate(.~Percentile, rzp_tree_mat_spend, mean)

write.csv(rzp_tree_exp_spend, 'CSV-Spend/rzp_tree_exp_spend.csv', row.names = FALSE)
write.csv(rzp_tree_mat_spend, 'CSV-Spend/rzp_tree_mat_spend.csv', row.names = FALSE)


######
# Tree simple Criterion
rzp_tree_exp_spend_simple <- rzp_tree_mat_spend_simple <- NULL
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
  rzp_tree_mat_spend_simple <- rbind(rzp_tree_mat_spend_simple, matching_evaluation(pred, "control"))
  
}
rzp_tree_exp_spend_simple <- aggregate(.~Percentile, rzp_tree_exp_spend_simple, mean)
rzp_tree_mat_spend_simple <- aggregate(.~Percentile, rzp_tree_mat_spend_simple, mean)

write.csv(rzp_tree_exp_spend_simple, 'CSV-Spend/rzp_tree_exp_spend_simple.csv', row.names = FALSE)
write.csv(rzp_tree_mat_spend_simple, 'CSV-Spend/rzp_tree_mat_spend_simple.csv', row.names = FALSE)


#######
# Causal Forest
c_forest_exp_spend <- c_forest_mat_spend <- NULL
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
  c_forest_mat_spend <- rbind(c_forest_mat_spend, matching_evaluation(causal_forest_pred, "control"))
}
c_forest_exp_spend <- aggregate(.~Percentile, c_forest_exp_spend, mean)
c_forest_mat_spend <- aggregate(.~Percentile, c_forest_mat_spend, mean)

write.csv(c_forest_exp_spend, 'CSV-Spend/c_forest_exp_spend.csv', row.names = FALSE)
write.csv(c_forest_mat_spend, 'CSV-Spend/c_forest_mat_spend.csv', row.names = FALSE)


####################################
####################################
# Import data and creat Train / Test splits
data <-  raw_email

data$visit <- NULL      
data$spendersion <- NULL

treatment <- "segment"
response <-  "spend"
control_level <- "No E-Mail"

treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]

#################
## Random Forest
sma_rf_spend_exp <- sma_rf_spend_mat <- NULL
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
  
  rf_pred_class <- dt_x_model_predictions(rf_models_list, test_data, response, treatment, control_level, "anova")
  
  
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
  
  
  
  sma_rf_spend_exp <- rbind(sma_rf_spend_exp, expected_percentile_response(rf_pred_class))
  sma_rf_spend_mat <- rbind(sma_rf_spend_mat, matching_evaluation(rf_pred_class, "control"))
}

sma_rf_spend_exp <- aggregate(.~Percentile, sma_rf_spend_exp, mean)
sma_rf_spend_mat <- aggregate(.~Percentile, sma_rf_spend_mat, mean)

write.csv(sma_rf_spend_exp, 'CSV-Spend/sma_rf_spend_exp.csv', row.names = FALSE)
write.csv(sma_rf_spend_mat, 'CSV-Spend/sma_rf_spend_mat.csv', row.names = FALSE)


#aggregate(conversion~fold+segment, raw_email, mean)
