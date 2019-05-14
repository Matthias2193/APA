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
insurance <- read.csv('Insurance.csv')


####################################################
# X Model Approach for binary response
####################################################
# Optional for same test train split
set.seed(21)

source('X Model Approach.R')

# Import data and creat Train / Test splits
data <-  read.csv('Email.csv')

## Prepare data to only include independent variables, treatment assignment and dependent variable
# First only look at conversion
data$visit <- NULL
data$spend <- NULL

## Set these values depending on the dataset 
# Treatment column in data
treatment <- "segment"
# Response column
response <- "conversion"
# String Level of the control Level in Treatment column
control_level <- "No E-Mail"

# Get the training and testing data as list
# creates k + 1 train test sets (for each T and control model)

train_test_list <- multiple_train_test_split(data, response, treatment, control_level, 0.7)
train_data <- train_test_list[1][[1]]
test_data <- train_test_list[[2]]

# Creates k + 1 Logit models for binary response
models_logit <- logit_models(train_data, response)


## Evaluate Models on Test data

# Calculated Upflift for each individual and Treatment and chosen Treatment
predictions <- logit_x_model_predictions(models_logit, test_data, response, treatment)

# Create merged test and prediction data for evaluation
eval_data <- merge_prediction_test_data(predictions, test_data, response, treatment, control_level)

# Expected outcome on Predicted Outcomes
# Formula from Zhao 2017
expected_outcome(eval_data)
# 1.23 % response

# AVG outcome from random assignment in test set
mean(eval_data$Outcome)
# 0.91 % response
