set.seed(21)

source('DecisionTreeImplementation.R')
source('X Model Approach.R')
source('Causal Forest.R')

#####################################################################################
### Conversion Prediction
#####################################################################################

#Data import
email <- read.csv('Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$spend <- email$segment <- NULL

response <- 'conversion'

# Split into test and train data
idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

test <- email[idx, ]

write.csv(test, 'test_conv.csv', row.names = FALSE)

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train_tree <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE)

####
# Fit Rzp Tree

raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list)

pruned_tree <- prune_tree(raw_tree,val,train_tree,target = 'conversion') #conversion
# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'rzp tree pred.csv', row.names = FALSE)


######
# Tree simple Criterion
raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list,criterion = 2)

pruned_tree <- raw_tree

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)


### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'rzp tree simple pred.csv', row.names = FALSE)


#######
# Causal Forest
causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)

write.csv(causal_forest_pred, "causal forest conv pred.csv", row.names = FALSE)

####################################
####################################
# Import data and creat Train / Test splits
data <-  read.csv('Email.csv')

## Prepare data to only include independent variables, treatment assignment and dependent variable
# Possible for conversion or visit

data$visit <- NULL      
#data$conversion <- NULL      
data$spend <- NULL

## Set these values depending on the dataset 
# Treatment column in data
treatment <- "segment"
# Response column
response <-  "conversion"
# String Level of the control Level in Treatment column
control_level <- "No E-Mail"

treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]

train_split <- data[-idx , ]
train_data <- list()
# make train_data list for each T
for(x in treatments){
  t_data <- train_split[train_split[, treatment] == x, ]
  
  # Remove the treatment column from training data
  t_data <- t_data[ , names(t_data) != treatment]
  
  # Each Training Set as one element in list
  train_data <- append(train_data, list(x = t_data[-idx, ]))
}
# Assign Treatment Names to Train Sets
names(train_data) <- treatments
# Add control group
t_data <- train_split[train_split[, treatment] == control_level, ]
t_data <- t_data[ , names(t_data) != treatment]
train_data <- append(train_data, list(Control = t_data[-idx, ]))

# Test Data
test_data <- data[idx , ]

#################
## Random Forest
rf_models_list <- rf_models(train_data, response, "class")

rf_pred_class <- dt_x_model_predictions(rf_models_list, test_data, response, treatment, control_level, "class")

write.csv(rf_pred_class, "rf conv pred.csv", row.names = FALSE)


#####################################################################################
### Spend Prediction
#####################################################################################
email <- read.csv('Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'

train <- email[-idx, ]

test <- email[idx, ]

write.csv(test, 'test_spend.csv', row.names = FALSE)

val <- train[p_idx,]
train_tree <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE)

####
# Rzp Tree
raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list, divergence = "EucDistance")
pruned_tree <- prune_tree(raw_tree,val,train_tree,target = response) 

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'rzp tree spend pred.csv', row.names = FALSE)


####################
## Rzp Tree Simple - no pruning
pruned_tree <- raw_tree

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)

# Calculate Uplift for each T
pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'rzp tree simple spend pred.csv', row.names = FALSE)


#######
# Causal Forest
causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
write.csv(causal_forest_pred, 'causal forest spend pred.csv', row.names = FALSE)

#########################################
#########################################
data <-  read.csv('Email.csv')
data$visit <- NULL      
data$conversion <- NULL

## Set these values depending on the dataset 
treatment <- "segment"
response <-  "spend" 
control_level <- "No E-Mail"

# creates k + 1 train test sets (for each T and control model)
treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]

train_split <- data[-idx , ]
train_data <- list()
# make train_data list for each T
for(x in treatments){
  t_data <- train_split[train_split[, treatment] == x, ]
  
  # Remove the treatment column from training data
  t_data <- t_data[ , names(t_data) != treatment]
  
  # Each Training Set as one element in list
  train_data <- append(train_data, list(x = t_data[-idx, ]))
}
# Assign Treatment Names to Train Sets
names(train_data) <- treatments
# Add control group
t_data <- train_split[train_split[, treatment] == control_level, ]
t_data <- t_data[ , names(t_data) != treatment]
train_data <- append(train_data, list(Control = t_data[-idx, ]))

# Test Data
test_data <- data[idx , ]

rf_models_list <- rf_models(train_data, response, "anova")

rf_pred_spend <- dt_x_model_predictions(rf_models_list, test_data, response, treatment, control_level, "anova")

write.csv(rf_pred_spend, "rf spend pred.csv", row.names = FALSE)



