# Implementation of the Seperate Model Approach for Multiple Treatments and different Base Learners
list.of.packages <- c("dplyr", "caret", "uplift")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(caret)
library(glmnet)

# Get the data
data <- read.csv('Email.csv')

# First only look at conversion
data$visit <- NULL
data$spend <- NULL


###For Debugging
treatment <- "segment"
response <- "conversion"
control_level <- "No E-Mail"
############


# Splits the input data into train and test for each treatment and control
# Input: Randomized data; response column as string; treatment column as string; level of control group in treatment column
multiple_train_test_split <- function(data, response, treatment, control_level){  
  treatments <- levels(as.factor(data[, treatment]))
  # remove the control level from treatments
  treatments <- treatments[! treatments %in% c(control_level)]

  # Create train data for each model
  # Partition data by level and split each partition into test train..
  train_T <- list()
  test_data <- data.frame()
  
  for(x in treatments){
    t_data <- data[data[, treatment] == x, ]
    
    idx <- createDataPartition(y = t_data[ , response], p=0.3, list = FALSE)
  
    # Combine all test splits into one Test Dataset
    test_data <- rbind(test_data, t_data[idx, ] )
    
    # Remove the treatment column from training data
    t_data <- t_data[ , names(t_data) != treatment]
    
    # Each Training Set as one element in list
    train_T <- append(train_T, list(x = t_data[-idx, ]))
  }
  # Assign Treatment Names to Train Sets
  names(train_T) <- treatments
  
  # Control Data
  c_data <- data[data[, treatment] == control_level, ]
  
  idx <- createDataPartition(y = c_data[ , response], p=0.3, list = FALSE)
  
  # Add the control group to the test
  test_data <- rbind(test_data, c_data[idx, ])
  # Remove the treatment from train set
  c_data <- c_data[ , names(c_data) != treatment]
  # Add to train data
  train_T <- append(train_T, list(Control = c_data[-idx, ]))
  
  
  # List of training data for each input and list of test data
  return(list(Train = train_T, Test =test_data))
}

# Get the training and testing data as list
res <- multiple_train_test_split(data, "conversion", "segment", "No E-Mail")
##TODO
# why is the array of dataframes stored in here ??
training_list <- res[1]
train_data <- training_list[[1]]

test_data <- res[[2]]


# Returns list of Logit models for training data
# only applicable for binary response e.g. conversion
logit_models <- function(train_data, response){  
  # for each element in train data fir model and return list of models
  treatments <- names(train_data)
  
  models <- list()
  for(i in c(1: length(treatments))){
    train <- train_data[[i]]
    X <- model.matrix(as.formula(paste(response, "~.-1")) , train_data[[i]])
    y <- as.matrix(train[, response])
    # Scale the data
    X <- scale(X)
    
    # Fit the Logit model
    logit <- glmnet(X, as.factor(y), alpha=0, family = 'binomial', lambda = 0)
    
    models <- append(models, list(x= logit))
  }
  
  names(models) <- treatments

  # Named list with fitted models for each treatment and control
  return(models)
}

# Create the models for each treatment
models_logit <- logit_models(train_data, "conversion")

#coef(models_logit[[1]])


# Prepare Test Data for predcition
response <- "conversion"
treatment <- "segment"

# do not include the treatment assignment in the Test data
X_test <- model.matrix(as.formula(paste(response, "~.-1")) , test_data[ , names(test_data) != treatment])
y_test <- as.matrix(test_data[, response])
# Scale the data
X_test <- scale(X_test)

for (i in c(1 : length(models))) {
  if(i == 1){
    predictions <- data.frame(predict(models_logit[[i]], X_test, type="response"))
  } else{
    predictions[ , paste0(i)] <- predict(models_logit[[i]], X_test, type="response")  
  }
}
# Rename columns to treatment names
colnames(predictions) <- names(models)


k <- ncol(predictions)
# For each treatment calculate the Uplift as T_i - Control
for (i in c(1: k - 1)) {
  predictions[ , paste0("Uplift", names(models)[i])] <- predictions[ , i] - predictions[ , k]
}


##TODO
# choose predicted treatment by the model
# 


##TODO
# either include the treatment assignment or the ID to merge later on

# prediction for each treatment and control
# 
calculate_uplift <- function(){
  
  
  # 
  return()
  
}



### Evaluation

# Input: -individual, outcome, assigned treatment
#        -individial, predicted treatment

# As described in Zhao et al. 2017
expected_outcome <- function(){
  
  
  
}
