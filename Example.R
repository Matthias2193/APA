#This script gives an example of how to use the implemented models, evaluate them and visualize the results


#Import the required libraries and files
library(caret)
source("src/Algorithm Implementations/DOM.R")
source("src/Algorithm Implementations/RzepakowskiTree.R")
source("src/Algorithm Implementations/ContextualTreatmentSelection.R")

source("src/Helper Functions/VisualizationHelper.R")
source("src/Helper Functions/Evaluation Methods.R")



set.seed(1234)

#The models are build using parallelizations. With this parameter users can decide how many cores should
#not be used for model building.
remain_cores <- 2

# Data import and preprocessing
email <- read.csv("Data/Email.csv")

#For categorical variables we need dummy variables with 0 and 1
email$men_treatment <- ifelse(email$segment == "Mens E-Mail", 1, 0)
email$women_treatment <- ifelse(email$segment == "Womens E-Mail", 1, 0)
email$control <- ifelse(email$segment == "No E-Mail", 1, 0)


email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

#Remove unused columns
email$visit <- email$conversion <- email$segment <- NULL

#The column names of the response, control and treatment variables
response <- "spend"
control <- "control"
treatment_list <- c("men_treatment", "women_treatment")

#Split the data into train and test set
idx <- createDataPartition(y = email[, response], p = 0.2, list = FALSE)
train <- email[-idx, ]
test <- email[idx, ]


#Difference in Outcome Model (DOM)

#Train the model
dom_forest <- dom_train(train, treatment_list, response, control,
  n_trees = 10, n_features = 3,
  criterion = "frac", remain_cores = remain_cores
)

#Make Predictions
pred <- predict_forest_df(dom_forest, test, treatment_list, control, remain_cores = remain_cores)

#Here two metrics are used for the evaluation.
#Expected outcome
dom_outcomes <- expected_outcome_curve(response, control, treatment_list, pred, model = "DOM")
#Qini
dom_qini <- qini_curve(pred, control, treatment_list, model_name = "DOM")


#The decision tree after Rzepakowski
rzp_forest <- rzp_train(
  train_data = train, treatment_list = treatment_list,
  response = response, control = control, n_trees = 10, n_features = 5,
  normalize = F, max_depth = 100, remain_cores = remain_cores,
  divergence = "binary_KL_divergence"
)

pred <- predict_forest_df(rzp_forest, test, treatment_list, control, remain_cores = remain_cores)

rzp_outcomes <- expected_outcome_curve(response, control, treatment_list, pred, model = "RZP")
rzp_qini <- qini_curve(pred, control, treatment_list, model_name = "RZP")


#Contextual Treatment Selection (CTS)
cts_forest <- cts_train(response, control, treatment_list, train,
                        ntree = 10, nrow(train), m_try = 4,
                        n_reg = 4, min_split = 10, parallel = TRUE, remain_cores = remain_cores
)
pred <- predict_forest_df(cts_forest, test, treatment_list, control, remain_cores = remain_cores)

cts_outcomes <- expected_outcome_curve(response, control, treatment_list, pred, model = "CTS")
cts_qini <- qini_curve(pred, control, treatment_list, model_name = "CTS")



#Visualize the evaluations 
visualize_outcome(rbind(dom_outcomes,rzp_outcomes,cts_outcomes),
                  ylabel = "Expected Conversion Probability per Person")

visualize_qini(rbind(dom_qini,rzp_qini,cts_qini), ylabel = "Cumulative Gained Conversion")

#By default all the results are displayed in one plot. 
#Using the "multiplot" parameter the results can displayed side by side in individual subplots

visualize_qini(rbind(dom_qini,rzp_qini,cts_qini), ylabel = "Cumulative Gained Conversion",
               multiplot = T)
