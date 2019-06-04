# Implementation of the Seperate Model Approach for Multiple Treatments and different Base Learners
list.of.packages <- c("caret", "glmnet", "rpart")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(caret)
library(glmnet)
library(rpart)
library(randomForest)

####################################
# X Model Data Training split
####################################

# Splits the input data into train and test for each treatment and control
# Input: Randomized data; response column as string; treatment column as string; level of control group in treatment column
multiple_train_test_split <- function(data, response, treatment, control_level, train_ratio){  
treatments <- levels(as.factor(data[, treatment]))
# remove the control level from treatments
treatments <- treatments[! treatments %in% c(control_level)]

# Create train data for each model
# Partition data by level and split each partition into test train..
train_T <- list()
test_data <- data.frame()

for(x in treatments){
  t_data <- data[data[, treatment] == x, ]
  
  
  idx <- createDataPartition(y = t_data[ , response], p=(1 - train_ratio), list = FALSE)

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

##TODO
# also include transformed single T case, --S treated or not treated


####################################
# Binary Logit models Prediction
####################################

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


logit_x_model_predictions <- function(models_logit, test_data, response, treatment, control_level){
# do not include the treatment assignment in the Test data
assignment <- test_data[ , treatment]
outcome <- test_data[, response]

# Scale the data
X_test <- model.matrix(as.formula(paste(response, "~.-1")) , test_data[ , names(test_data) != treatment])
X_test <- scale(X_test)

for (i in c(1 : length(models_logit))) {
  if(i == 1){
    predictions <- data.frame(predict(models_logit[[i]], X_test, type="response"))
  } else{
    predictions[ , paste0(i)] <- as.numeric(predict(models_logit[[i]], X_test, type="response") )
  }
}
# Rename columns to treatment names
colnames(predictions) <- names(models_logit)

k <- ncol(predictions)
# For each treatment calculate the Uplift as T_i - Control
for (i in c(1: k - 1)) {
  predictions[ , paste0("Uplift - ", names(models_logit)[i])] <- predictions[ , i] - predictions[ , k]
}

# choose predicted treatment by the model
predictions$T_index <- apply(predictions[, 1:k], 1, which.max)
predictions$Treatment <- colnames(predictions)[predictions$T_index]

# Add actual assignment 
predictions$Outcome <- outcome

predictions$Assignment <- assignment
predictions$Assignment <- ifelse(predictions$Assignment == control_level, "Control", as.character(predictions$Assignment))


return(predictions)
}


####################################
# Decision Tree model
####################################

dt_models <- function(train_data, response, prediction_method){  
  # for each element in train data fir model and return list of models
  treatments <- names(train_data)
  
  models <- list()
  for(i in c(1: length(treatments))){
    train <- train_data[[i]]
    
    # Grow the tree
    ## TODO how to set / adjust parameters ??
    dt <- rpart(as.formula(paste(response, "~.")), train, method = prediction_method, control=rpart.control(minsplit=10,  cp=0.001))
    
    models <- append(models, list(x= dt))
  }
  
  names(models) <- treatments
  
  # Named list with fitted models for each treatment and control
  return(models)
}


dt_x_model_predictions <- function(models_dt, test_data, response, treatment, control_level, prediction_class){
  assignment <- test_data[ , treatment]
  outcome <- test_data[, response]
  
  test_data[ , treatment] <- NULL
  
  # Check whether regression or class prediction due to different dim of predict() result
  if(prediction_class == "class"){
    for (i in c(1 : length(models_dt))) {
      if(i == 1){
        predictions <- data.frame(predict(models_dt[[i]], test_data, type = "prob")[ , 2])
      } else{
        predictions[ , paste0(i)] <- as.numeric(predict(models_dt[[i]], test_data, type = "prob")[ , 2] )
      }
    }  
  } else if(prediction_class == "anova"){
    for (i in c(1 : length(models_dt))) {
      if(i == 1){
        predictions <- data.frame(predict(models_dt[[i]], test_data))
      } else{
        predictions[ , paste0(i)] <- as.numeric(predict(models_dt[[i]], test_data) )
      }
    }
    
  }
  
  
  # Rename columns to treatment names
  colnames(predictions) <- names(models_dt)
  
  k <- ncol(predictions)
  # For each treatment calculate the Uplift as T_i - Control
  for (i in c(1: k - 1)) {
    predictions[ , paste0("Uplift - ", names(models_dt)[i])] <- predictions[ , i] - predictions[ , k]
  }
  
  # choose predicted treatment by the model
  predictions$T_index <- apply(predictions[, 1:3], 1, which.max)
  predictions$Treatment <- colnames(predictions)[predictions$T_index]
  
  predictions$Outcome <- outcome
  
  predictions$Assignment <- assignment
  predictions$Assignment <- ifelse(predictions$Assignment == control_level, "Control", as.character(predictions$Assignment))
  
  
  return(predictions)
}


#######################
# RF

rf_models <- function(train_data, response, prediction_method){  
  # for each element in train data fir model and return list of models
  
  if(prediction_method == "class"){
    treatments <- names(train_data)
  }
  
  models <- list()
  for(i in c(1: length(treatments))){
    train <- train_data[[i]]
    
    train[, response] <- as.factor(train[, response])
    
    ## TODO how to set / adjust parameters ??
    rf <- randomForest(as.formula(paste(response, "~.")), data = train, mtry=3, ntree = 350)
    
    models <- append(models, list(x= rf))
  }
  
  names(models) <- treatments
  
  # Named list with fitted models for each treatment and control
  return(models)
}


