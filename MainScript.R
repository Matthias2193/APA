#The main script for data analysis and model evaluation


#Data import
email <- read.csv('Email.csv')
email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$segment <- NULL
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)
#insurance <- read.csv('Insurance.csv')


source('Evaluation Methods.R')

library(ggplot2)

####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')

library(caret)

response <- 'conversion'

# Split into test and train data
idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]
test <- email[idx, ]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE)

raw_tree <- create_node(train,0,100,treatment_list,'conversion','control',test_list)

# Partition training data for pruning
idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
# takes a lot of time..
pruned_tree <- prune_tree(raw_tree,train[idx,],email[-idx,],target = 'conversion')

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)


### Results Preparation to bring into equal format

# Calculate Uplift for each T
pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
pred[ , "Uplift - Womens E-Mail"] <- pred[ , 1] - pred[ , 3]

pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 12:14], 1, which.max) + 11]

# How often is each Treatment assigned by model
table(pred$Treatment)

# Expected outcome
# Formula from Zhao 2017
expected_outcome(pred)

# Expected Response per targeted customers
perc_eval_zhao <- expected_percentile_response(pred)
plot(perc_eval_zhao)


####################################################
# X Model Approach for continuous response
####################################################
set.seed(123)
source('X Model Approach.R')

data <-  read.csv('Email.csv')

## Prepare data to only include independent variables, treatment assignment and dependent variable
# Possible for conversion or visit

data$visit <- NULL      
data$conversion <- NULL

## Set these values depending on the dataset 
# Treatment column in data
treatment <- "segment"
# Response column
response <-  "spend" 
# String Level of the control Level in Treatment column
control_level <- "No E-Mail"

# Get the training and testing data as list
# creates k + 1 train test sets (for each T and control model)

train_test_list <- multiple_train_test_split(data, response, treatment, control_level, 0.7)
train_data <- train_test_list[1][[1]]
test_data <- train_test_list[[2]]



####################################################
# Decision Tree - X Model
####################################################

models_dt <- dt_models(train_data, response, "anova")

## Evaluate Models on Test data

# Calculated Upflift for each individual and Treatment and chosen Treatment
predictions_dt <- dt_x_model_predictions(models_dt, test_data, response, treatment, control_level, "anova")

# How often is each Treatment assigned by model
table(predictions_dt$Treatment)

# Expected outcome
# Formula from Zhao 2017
expected_outcome(predictions_dt)
# 1.32 AVG Spend

# With random assignment
mean(predictions_dt$Outcome)
# 0.99 AVG Spend

# Expected Response per targeted customers
perc_eval_zhao <- expected_percentile_response(predictions_dt)
plot(perc_eval_zhao)


# Rzepakowski 2012
# for each T uplift curve and one dynamic curve for combined treatments
curves <- matching_evaluation(predictions_dt)

library(reshape2)

melt(curves, id="Percentile") %>% ggplot(aes(x = Percentile, y = value, color = variable)) +
  geom_line() + 
  labs(
    y = "Uplift"
  ) +
  theme_light()
  
curves %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = curves[, colnames(curves)[4]]) ) + 
  labs(
    title = "Dynamic Treatment Curve - SMA DT",
    y = "Uplift"
  ) +
  theme_light()


####################################################
# X Model Approach for binary response
####################################################

# Import data and creat Train / Test splits
data <-  read.csv('Email.csv')

## Prepare data to only include independent variables, treatment assignment and dependent variable
# Possible for conversion or visit

data$visit <- NULL      
data$spend <- NULL

## Set these values depending on the dataset 
# Treatment column in data
treatment <- "segment"
# Response column
response <-  "conversion"  #"visit"
# String Level of the control Level in Treatment column
control_level <- "No E-Mail"

# Get the training and testing data as list
# creates k + 1 train test sets (for each T and control model)

train_test_list <- multiple_train_test_split(data, response, treatment, control_level, 0.7)
train_data <- train_test_list[1][[1]]
test_data <- train_test_list[[2]]


####################################################
# Logit - X Model
####################################################

# Creates k + 1 Logit models for binary response
models_logit <- logit_models(train_data, response)

## Evaluate Models on Test data
# Calculated Upflift for each individual and Treatment and chosen Treatment
predictions_logit <- logit_x_model_predictions(models_logit, test_data, response, treatment, control_level)

# How often is each Treatment assigned by model
table(predictions_logit$Treatment)

# Expected outcome on Predicted Outcomes
# Formula from Zhao 2017
expected_outcome(predictions_logit)
# 1.19 % response

# AVG outcome from random assignment in test set
mean(predictions_logit$Outcome)
# 0.94 % response


# Rzepakowski
perc_matching <- matching_evaluation(predictions_logit)

log_perc_rzpk <- matching_evaluation(predictions_logit)

# For each treatment
melt(log_perc_rzpk, id="Percentile") %>% ggplot(aes(x = Percentile, y = value, color = variable)) +
  geom_line() + 
  labs(
    y = "Uplift"
  ) +
  theme_light()

# Only for
log_perc_rzpk %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = log_perc_rzpk[, colnames(log_perc_rzpk)[4]]) ) + 
  labs(
    title = "Dynamic Treatment Curve - SMA Logit",
    y = "Uplift"
  ) +
  theme_light()


# Zhao 
logit_perc_eval_zhao <- expected_percentile_response(predictions_logit)

logit_perc_eval_zhao %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = logit_perc_eval_zhao[, colnames(logit_perc_eval_zhao)[2]]) ) + 
  labs(
    title = "Expected Outcome - SMA Logit",
    y = "Expected Response",
    x ="Amount of Treated"
  ) +
  theme_light()


####################################################
# Decision Tree - X Model
####################################################

models_clas_dt <- dt_models(train_data, response, "class")

## Evaluate Models on Test data

# Calculated Upflift for each individual and Treatment and chosen Treatment
predictions_clas_dt <- dt_x_model_predictions(models_clas_dt, test_data, response, treatment, control_level, "class")

# How often is each Treatment assigned by model
table(predictions_clas_dt$Treatment)

# Rzepakowski 2012
dt_perc_rzpk <- matching_evaluation(predictions_clas_dt)

dt_perc_rzpk %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = dt_perc_rzpk[, colnames(dt_perc_rzpk)[4]]) ) + 
  labs(
    title = "Dynamic Treatment Curve - SMA Logit",
    y = "Uplift"
  ) +
  theme_light()



# Zhao 2017
expected_outcome(predictions_clas_dt)
# 0.98 % response

# Expected Response per targeted customers
perc_eval_zhao <- expected_percentile_response(predictions_clas_dt)

perc_eval_zhao %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = perc_eval_zhao[, colnames(perc_eval_zhao)[2]]) ) + 
  labs(
    title = "Expected Outcome - SMA DT",
    y = "Expected Response",
    x ="Amount of Treated"
  ) +
  theme_light()


####################################################
# Decision Tree - X Model
####################################################

rf_models <- rf_models(train_data, response, "class")

rf_pred_class <- dt_x_model_predictions(rf_models, test_data, response, treatment, control_level, "class")

# Expected Response per targeted customers
rf_perc_eval_zhao <- expected_percentile_response(rf_pred_class)

rf_perc_eval_zhao %>% ggplot(aes(x = Percentile)) +
  geom_line(aes_string(y = rf_perc_eval_zhao[, colnames(rf_perc_eval_zhao)[2]]) ) + 
  labs(
    title = "Expected Outcome - SMA RF",
    y = "Expected Response",
    x ="Amount of Treated"
  ) +
  theme_light()


